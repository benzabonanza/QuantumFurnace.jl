"""
Run-time cache for trajectory simulations
"""
struct TrajectoryWorkspace{T}
    jump_oft::Matrix{T}   # buffer for A_ω (or its dagger-partner is handled by swapping mul order)
    psi_tmp::Vector{T}    # generic tmp for matvec results (K0ψ, Aωψ, Uresψ, etc.)
    Rpsi::Vector{T}       # cache Rψ so K0ψ can be formed as ψ - α(Rψ) without a second matvec
end

function TrajectoryWorkspace(::Type{T}, dim::Int) where {T}
    return TrajectoryWorkspace{T}(
        zeros(T, dim, dim),  # jump_oft
        zeros(T, dim),       # psi_tmp
        zeros(T, dim),       # Rpsi
    )
end

"""
Per-operator Kraus data for Lie-Trotter splitting.
Each operator `a` has its own R_a, K0_a, U_residual_a, U_B_a.
"""
struct PerOperatorKraus{T}
    R::Matrix{T}                        # R^a for this operator (rescaled by 1/p_jump)
    K0::Matrix{T}                       # I - alpha * R^a
    U_residual::Matrix{T}               # sqrt(S^a) via eigendecomposition
    U_B::Union{Nothing, Matrix{T}}      # exp(-i * delta_eff * B^a), or nothing
end

struct TrajectoryFramework{T,D<:AbstractDomain}
    domain::D
    jumps::Vector{JumpOp}
    ham_or_trott::Union{HamHam, TrottTrott}
    config::AbstractThermalizeConfig{D}
    precomputed_data::Any  # NamedTuple from precompute_data, varies by domain

    # Per-operator Kraus data (Lie-Trotter splitting)
    per_operator::Vector{PerOperatorKraus{T}}
    n_jumps::Int

    # Step parameters
    delta::Float64           # original delta (for time stepping / num_steps)
    delta_eff::Float64       # delta / p_jump = delta * n_jumps (for per-operator Kraus probabilities)
    alpha::Float64           # α = 1 - sqrt(1-δ_eff)

    # Runtime workspace (mutated during stepping)
    ws::TrajectoryWorkspace{T}
end

function build_trajectoryframework(
    jumps::AbstractVector{<:JumpOp},
    ham_or_trott::Union{HamHam, TrottTrott},
    config::AbstractThermalizeConfig,
    precomputed_data,
    scratch::KrausScratch{ComplexF64},
    delta::Float64)

    dim = size(jumps[1].data, 1)
    n_jumps = length(jumps)
    p_jump = 1.0 / n_jumps
    delta_eff = delta / p_jump   # = delta * n_jumps

    @assert delta_eff < 1.0 "delta_eff = $(delta_eff) >= 1.0: delta=$(delta) * n_jumps=$(n_jumps) is too large for per-operator splitting"

    alpha = 1 - sqrt(1 - delta_eff)

    # For TrotterDomain, transform jump operators from Hamiltonian eigenbasis
    # to Trotter eigenbasis. The NUFFT prefactors use Trotter quasi-Bohr frequencies,
    # so the element-wise product A .* P requires A in the same basis.
    jumps_for_diss = if config.domain isa TrotterDomain && ham_or_trott isa TrottTrott
        transform_jumps_to_basis(jumps, ham_or_trott.eigvecs)
    else
        collect(JumpOp, jumps)
    end

    # Precompute per-operator coherent B terms (one per jump)
    # NOTE: precompute_coherent_total_B handles its own Trotter basis transform internally
    # (via B_trotter in coherent.jl), so we pass the ORIGINAL jumps here.
    per_op_U_B = Vector{Union{Nothing, Matrix{ComplexF64}}}(undef, n_jumps)
    if config.with_coherent
        @inbounds for a in 1:n_jumps
            single_jump = JumpOp[jumps[a]]  # Force Vector{JumpOp} for dispatch compatibility
            B_a = precompute_coherent_total_B(single_jump, ham_or_trott, config, precomputed_data)
            hermitianize!(B_a)
            per_op_U_B[a] = exp(-1im * delta_eff * Hermitian(B_a))
        end
    else
        fill!(per_op_U_B, nothing)
    end

    # Build per-operator Kraus data
    per_operator = Vector{PerOperatorKraus{ComplexF64}}(undef, n_jumps)

    @inbounds for a in 1:n_jumps
        # Compute R^a for single operator, rescaled by 1/p_jump
        precompute_R([jumps_for_diss[a]], ham_or_trott, config, precomputed_data, scratch)
        R_a = copy(scratch.R)
        R_a .*= (1.0 / p_jump)   # rescale: R_a = (1/p_jump) * sum_w rate2(w) * A_w' * A_w

        # K0_a = I - alpha * R_a
        K0_a = copy(R_a)
        K0_a .*= (-alpha)
        @inbounds for i in 1:dim
            K0_a[i,i] += 1
        end

        # S_a = (2*alpha - delta_eff)*R_a - alpha^2 * R_a^2
        mul!(scratch.tmp1, R_a, R_a)   # tmp1 := R_a^2
        s1 = 2 * alpha - delta_eff
        s2 = alpha * alpha
        @. scratch.tmp2 = s1 * R_a - s2 * scratch.tmp1

        # TFIX-04: PSD guard -- clamp negative eigenvalues to zero (silent fallback)
        hermitianize!(scratch.tmp2)
        S_herm = Hermitian(scratch.tmp2)
        eig = eigen(S_herm)
        eig.values .= max.(eig.values, 0.0)
        U_residual_a = Matrix{ComplexF64}(Diagonal(sqrt.(eig.values)) * eig.vectors')

        per_operator[a] = PerOperatorKraus(R_a, K0_a, U_residual_a, per_op_U_B[a])
    end

    ws = TrajectoryWorkspace(ComplexF64, dim)

    return TrajectoryFramework(
        config.domain,
        jumps_for_diss,
        ham_or_trott,
        config,
        precomputed_data,
        per_operator,
        n_jumps,
        delta,
        delta_eff,
        alpha,
        ws,
    )
end

"""
    precompute_R(domain, jumps, ham_or_trott, config, precomputed_data, scratch) -> Matrix{ComplexF64}

    Compute
        R = ∑_{k>0} L_k† L_k
    in the same basis as `jump.in_eigenbasis` (Hamiltonian eigenbasis for `HamHam`, Trotter basis for `TrottTrott`).

    Conventions are matched to `jump_contribution!(domain, ::AbstractThermalizeConfig, ...)`:
    - the weights are `rate2(ω) = base_prefactor * transition(ω)` (no extra `δ` factor),
    - for Hermitian jumps we iterate half-grid and add the mirrored negative-frequency partner explicitly.

    This returns `scratch.R` (Hermitianized).
"""
function precompute_R(
    jumps::AbstractVector{<:JumpOp},
    hamiltonian::HamHam,
    config::AbstractThermalizeConfig{EnergyDomain},
    precomputed_data,
    scratch::KrausScratch{ComplexF64},
    )
    dim = size(hamiltonian.data, 1)
    (; transition, gamma_norm_factor, energy_labels) = precomputed_data


    base_prefactor = config.w0 / (config.sigma * sqrt(2 * pi)) * gamma_norm_factor

    fill!(scratch.R, 0)

    @inbounds for jump in jumps
        if jump.hermitian
            # Half-grid (w_raw <= 0) and mirror partner at -w using Aω† as the Lindblad operator.
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)

                # Aω := A ∘ exp(-(w-ν)^2/(4σ^2))   (elementwise in eigenbasis)
                oft!(scratch.jump_oft, jump, w, hamiltonian, config.sigma)

                # Positive-frequency contribution: rate2(w) * (Aω† Aω)
                rate2_pos = base_prefactor * transition(w)
                mul!(scratch.LdagL, scratch.jump_oft', scratch.jump_oft)
                @. scratch.R += rate2_pos * scratch.LdagL

                if w > 1e-12
                    # Negative-frequency partner uses Lindblad op = Aω†:
                    # contribution is rate2(-w) * (Aω Aω†)
                    rate2_neg = base_prefactor * transition(-w)
                    mul!(scratch.LdagL, scratch.jump_oft, scratch.jump_oft')
                    @. scratch.R += rate2_neg * scratch.LdagL
                end
            end
        else
            # Non-Hermitian jump: full grid, no mirroring shortcut.
            for w in energy_labels
                oft!(scratch.jump_oft, jump, w, hamiltonian, config.sigma)

                rate2 = base_prefactor * transition(w)
                mul!(scratch.LdagL, scratch.jump_oft', scratch.jump_oft)
                @. scratch.R += rate2 * scratch.LdagL
            end
        end
    end

    # Numerical Hermitianization (R should be Hermitian PSD by construction).
    hermitianize!(scratch.R)
    return scratch.R
end


function precompute_R(
    jumps::AbstractVector{<:JumpOp},
    ham_or_trott::Union{HamHam, TrottTrott},
    config::AbstractThermalizeConfig{D},
    precomputed_data,
    scratch::KrausScratch{ComplexF64},
    ) where {D<:Union{TimeDomain, TrotterDomain}}
    dim = size(jumps[1].in_eigenbasis, 1)
    (; transition, gamma_norm_factor, energy_labels, oft_nufft_prefactors) = precomputed_data

    # Same weight as in jump_contribution!(::Union{TimeDomain,TrotterDomain}, ...), but without δ.
    base_prefactor = config.w0 * config.t0^2 * (config.sigma * sqrt(2 / pi)) / (2 * pi) * gamma_norm_factor

    fill!(scratch.R, 0)

    @inbounds for jump in jumps
        if jump.hermitian
            # Half-grid and explicit mirror, allocation-free (no filter/abs vector creation).
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)

                nufft_prefactor_matrix = prefactor_view(oft_nufft_prefactors, w)

                # Aω := A ∘ prefactor_matrix(ω)  (elementwise)
                @. scratch.jump_oft = jump.in_eigenbasis * nufft_prefactor_matrix

                rate2_pos = base_prefactor * transition(w)
                mul!(scratch.LdagL, scratch.jump_oft', scratch.jump_oft)
                @. scratch.R += rate2_pos * scratch.LdagL

                if w > 1e-12
                    rate2_neg = base_prefactor * transition(-w)
                    mul!(scratch.LdagL, scratch.jump_oft, scratch.jump_oft')
                    @. scratch.R += rate2_neg * scratch.LdagL
                end
            end
        else
            for w in energy_labels
                nufft_prefactor_matrix = prefactor_view(oft_nufft_prefactors, w)
                @. scratch.jump_oft = jump.in_eigenbasis * nufft_prefactor_matrix

                rate2 = base_prefactor * transition(w)
                mul!(scratch.LdagL, scratch.jump_oft', scratch.jump_oft)
                @. scratch.R += rate2 * scratch.LdagL
            end
        end
    end

    hermitianize!(scratch.R)
    return scratch.R
end

# Fast squared norm without sqrt
@inline _norm2(v::AbstractVector{ComplexF64}) = real(dot(v, v))

@inline function _accumulate_density_matrix!(
    rho_acc::Matrix{ComplexF64},
    psi::Vector{ComplexF64},
)
    # rho_acc += psi * psi'   (conjugate transpose, rank-1 update)
    mul!(rho_acc, psi, psi', one(ComplexF64), one(ComplexF64))
    return nothing
end

"""
    evolve_along_trajectory(psi0, fw::TrajectoryFramework, total_time) -> Vector{ComplexF64}

    Runs `num_steps = ceil(total_time/δ)` full δ-steps using `step_along_trajectory!`.
    Does not allocate inside the loop (beyond RNG internals); all workspace is in `fw.ws`.

    Normalization: each step normalizes in the branch logic (and after U_B if present).
"""
function _evolve_along_trajectory!(
    psi::Vector{ComplexF64},
    fw::TrajectoryFramework{ComplexF64},
    total_time::Float64,
    )::Vector{ComplexF64}

    delta = fw.delta
    @assert delta > 0
    @assert total_time ≥ 0

    num_steps = ceil(Int, total_time / delta)
    # ensure normalized input (optional, but makes probabilities consistent)
    psi_norm2 = _norm2(psi)
    rmul!(psi, 1.0 / sqrt(max(psi_norm2, eps(Float64))))

    @inbounds for _ in 1:num_steps
        step_along_trajectory!(psi, fw)  # dispatches to EnergyDomain or Time/Trotter variant
    end

    return psi
end

# Allocation-free measure-add (writes into `acc[:, save_idx]` by +=)
function _accumulate_measurements!(
    acc::AbstractMatrix{Float64},
    save_idx::Int,
    psi::Vector{ComplexF64},
    observables::Vector{Matrix{ComplexF64}},
    tmp::Vector{ComplexF64},
)
    @inbounds for i in eachindex(observables)
        mul!(tmp, observables[i], psi)                 # tmp := O_i ψ
        acc[i, save_idx] += real(dot(psi, tmp))        # += <ψ|O_i|ψ>
    end
    return nothing
end

"""
    run_trajectories(jumps, config, psi0, hamiltonian; trotter=nothing,
                     total_time=config.mixing_time, delta=config.delta,
                     ntraj=1, observables=nothing, save_every=1, store_states=false)

    Builds the TrajectoryFramework once, then runs `ntraj` trajectories.

    - If `observables === nothing`: returns final state(s) only.
    - If `observables` provided: returns time grid and trajectory-averaged ⟨O_i⟩(t).
"""
function run_trajectories(
    jumps::Vector{JumpOp},
    config::AbstractThermalizeConfig,
    psi0::Vector{ComplexF64},
    hamiltonian::HamHam;
    trotter::Union{TrottTrott,Nothing}=nothing,
    total_time::Float64 = config.mixing_time,
    delta::Float64 = config.delta,
    ntraj::Int = 1,
    observables::Union{Nothing, Vector{Matrix{ComplexF64}}} = nothing,
    save_every::Int = 1,
    store_states::Bool = false,
)

    validate_config!(config)
    print_press(config)

    @assert ntraj ≥ 1
    @assert save_every ≥ 1

    # Choose evolution object consistent with domain
    ham_or_trott = if config.domain isa TrotterDomain
        trotter === nothing && error("TrotterDomain requires `trotter`.")
        trotter
    else
        hamiltonian
    end

    precomputed_data = precompute_data(config, ham_or_trott)

    dim = size(hamiltonian.data, 1)
    rho_mean = zeros(ComplexF64, dim, dim)
    builder_scratch = KrausScratch(ComplexF64, dim)

    fw = build_trajectoryframework(jumps, ham_or_trott, config, precomputed_data, builder_scratch, delta)

    # ------------------------------------------------------------------
    # No measurements: just run and (optionally) store final ψ per traj
    # ------------------------------------------------------------------
    if observables === nothing
        if ntraj == 1 && !store_states
            psi = copy(psi0)
            _evolve_along_trajectory!(psi, fw, total_time)
            _accumulate_density_matrix!(rho_mean, psi)
            hermitianize!(rho_mean)
            return (framework = fw, psi = psi, rho_mean = rho_mean)
        end

        states = store_states ? Matrix{ComplexF64}(undef, dim, ntraj) : nothing
        psi = copy(psi0)

        @inbounds for trajectory in 1:ntraj
            copyto!(psi, psi0)
            _evolve_along_trajectory!(psi, fw, total_time)
            _accumulate_density_matrix!(rho_mean, psi)
            if store_states
                @views states[:, trajectory] .= psi
            end
        end

        rho_mean ./= ntraj
        hermitianize!(rho_mean)

        return (framework = fw, states = states, rho_mean = rho_mean)
    end

    # ------------------------------------------------------------------
    # With measurements: average ⟨O⟩ over trajectories, saved every `save_every`
    # ------------------------------------------------------------------
    delta = fw.delta
    num_steps = ceil(Int, total_time / delta)
    num_saves = div(num_steps, save_every) + 1
    num_obs   = length(observables)

    times = Vector{Float64}(undef, num_saves)
    @inbounds for s in 1:num_saves
        times[s] = (s - 1) * save_every * delta
    end

    mean_data = zeros(Float64, num_obs, num_saves)

    # reuse fw workspace vector as gemv buffer for measurement (safe between steps)
    tmp_meas = fw.ws.psi_tmp
    psi = copy(psi0)

    @inbounds for _ in 1:ntraj
        copyto!(psi, psi0)

        # normalize once per trajectory
        n2 = real(dot(psi, psi))
        rmul!(psi, 1.0 / sqrt(max(n2, eps(Float64))))

        # save at t=0
        _accumulate_measurements!(mean_data, 1, psi, observables, tmp_meas)

        save_idx = 1
        for step in 1:num_steps
            step_along_trajectory!(psi, fw)

            if step % save_every == 0
                save_idx += 1
                _accumulate_measurements!(mean_data, save_idx, psi, observables, tmp_meas)
            end
        end
        _accumulate_density_matrix!(rho_mean, psi)
    end

    mean_data ./= ntraj
    rho_mean ./= ntraj
    hermitianize!(rho_mean)

    return (framework = fw, times = times, measurements_mean = mean_data, rho_mean = rho_mean)
end

"""
    step_along_trajectory!(psi, fw::TrajectoryFramework)  for D ∈ {TimeDomain, TrotterDomain}

    Per-operator Lie-Trotter splitting: randomly select ONE operator a, apply that operator's
    CPTP channel (K0_a, U_res_a, jump outcomes for operator a only).

    # Per-operator channel structure (Chen 2023, adapted for Lie-Trotter splitting):
    #   Pick a ∈ {1,...,N_jumps} uniformly at random
    #   K0_a = I - alpha*R_a, where alpha = 1 - sqrt(1-delta_eff), delta_eff = delta*N_jumps
    #   K_{a,w} = sqrt(delta_eff * scaled_rate(w)) * L_{a,w}  (jump operators for operator a)
    #   U_res_a: U_res_a'*U_res_a = S_a  (residual for operator a)
    # Rates rescaled by 1/p_jump so net effect per unit time matches DM run_thermalization.
"""
function step_along_trajectory!(
    psi::Vector{ComplexF64},
    fw::TrajectoryFramework{ComplexF64,D},
    ) where {D<:Union{TimeDomain,TrotterDomain}}

    ws = fw.ws
    cfg = fw.config
    pd  = fw.precomputed_data

    # Pull hot fields into locals
    transition         = pd.transition
    energy_labels      = pd.energy_labels
    oft_prefactors     = pd.oft_nufft_prefactors
    gamma_norm_factor  = pd.gamma_norm_factor

    delta_eff = fw.delta_eff

    # Rate prefactor with 1/p_jump rescaling baked in
    # base_prefactor already includes gamma_norm_factor; multiply by n_jumps for 1/p_jump
    scaled_prefactor = cfg.w0 * cfg.t0^2 * (cfg.sigma * sqrt(2 / pi)) / (2 * pi) * gamma_norm_factor / (1.0 / fw.n_jumps)

    # ------------------------------------------------------------------
    # Random operator selection (Lie-Trotter splitting)
    # ------------------------------------------------------------------
    a = rand(1:fw.n_jumps)
    per_op = fw.per_operator[a]
    jump = fw.jumps[a]

    # ------------------------------------------------------------------
    # Apply per-operator coherent unitary FIRST (matches DM code ordering)
    # ------------------------------------------------------------------
    if per_op.U_B !== nothing
        mul!(ws.psi_tmp, per_op.U_B, psi)
        copyto!(psi, ws.psi_tmp)
        rmul!(psi, 1.0 / sqrt(max(_norm2(psi), eps(Float64))))
    end

    # ------------------------------------------------------------------
    # Compute probabilities from per-operator R_a, K0_a, U_residual_a
    # ------------------------------------------------------------------
    mul!(ws.Rpsi, per_op.R, psi)                    # R_a * psi
    expR = max(real(dot(psi, ws.Rpsi)), 0.0)        # <psi|R_a|psi>

    mul!(ws.psi_tmp, per_op.K0, psi)                # K0_a * psi
    p_nojump = _norm2(ws.psi_tmp)

    p_jump_total = delta_eff * expR                  # delta_eff, NOT delta

    mul!(ws.Rpsi, per_op.U_residual, psi)           # U_res_a * psi (reuse buffer)
    p_res = _norm2(ws.Rpsi)

    total_weight = p_nojump + p_res + p_jump_total
    total_weight = max(total_weight, 0.0)

    if abs(total_weight - 1.0) > 1e-6
        @warn "Normalization violation: sum = $(round(total_weight; digits=6))"
    end

    r = rand() * total_weight

    # ------------------------------------------------------------------
    # Branch: no-jump / residual / dissipative jump (for selected operator a only)
    # ------------------------------------------------------------------
    if r < p_nojump
        copyto!(psi, ws.psi_tmp)
        rmul!(psi, 1.0 / sqrt(max(p_nojump, eps(Float64))))

    elseif r < (p_nojump + p_res)
        copyto!(psi, ws.Rpsi)
        rmul!(psi, 1.0 / sqrt(max(p_res, eps(Float64))))

    else
        # Dissipative jump: sample one (ω,±) outcome for the SINGLE selected operator a
        target = r - p_nojump - p_res
        csum   = 0.0
        chosen = false
        last_norm2 = 0.0

        @inbounds if jump.hermitian
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)

                pref = prefactor_view(oft_prefactors, w)
                @. ws.jump_oft = jump.in_eigenbasis * pref

                # Positive-frequency
                mul!(ws.Rpsi, ws.jump_oft, psi)
                n2 = _norm2(ws.Rpsi)
                last_norm2 = n2
                p = delta_eff * (scaled_prefactor * transition(w)) * n2
                csum += p
                if csum >= target
                    copyto!(psi, ws.Rpsi)
                    rmul!(psi, 1.0 / sqrt(max(n2, eps(Float64))))
                    chosen = true
                    break
                end

                # Negative-frequency partner
                if w > 1e-12
                    mul!(ws.Rpsi, ws.jump_oft', psi)
                    n2 = _norm2(ws.Rpsi)
                    last_norm2 = n2
                    p = delta_eff * (scaled_prefactor * transition(-w)) * n2
                    csum += p
                    if csum >= target
                        copyto!(psi, ws.Rpsi)
                        rmul!(psi, 1.0 / sqrt(max(n2, eps(Float64))))
                        chosen = true
                        break
                    end
                end
            end
        else
            for w in energy_labels
                pref = prefactor_view(oft_prefactors, w)
                @. ws.jump_oft = jump.in_eigenbasis * pref

                mul!(ws.Rpsi, ws.jump_oft, psi)
                n2 = _norm2(ws.Rpsi)
                last_norm2 = n2
                p = delta_eff * (scaled_prefactor * transition(w)) * n2
                csum += p
                if csum >= target
                    copyto!(psi, ws.Rpsi)
                    rmul!(psi, 1.0 / sqrt(max(n2, eps(Float64))))
                    chosen = true
                    break
                end
            end
        end

        if !chosen
            copyto!(psi, ws.Rpsi)
            rmul!(psi, 1.0 / sqrt(max(last_norm2, eps(Float64))))
        end
    end
    return nothing
end

"""
    step_along_trajectory!(psi, fw::TrajectoryFramework)  for D == EnergyDomain

    Per-operator Lie-Trotter splitting for EnergyDomain.
    Same logic as Time/Trotter variant but A_{a,ω} is generated by `oft!(...)`.
"""
function step_along_trajectory!(
    psi::Vector{ComplexF64},
    fw::TrajectoryFramework{ComplexF64,EnergyDomain},
    )

    ws = fw.ws
    cfg = fw.config
    pd  = fw.precomputed_data

    transition         = pd.transition
    energy_labels      = pd.energy_labels
    gamma_norm_factor  = pd.gamma_norm_factor

    delta_eff = fw.delta_eff

    # Rate prefactor with 1/p_jump rescaling: base * n_jumps
    scaled_prefactor = cfg.w0 / (cfg.sigma * sqrt(2 * pi)) * gamma_norm_factor / (1.0 / fw.n_jumps)

    # ------------------------------------------------------------------
    # Random operator selection (Lie-Trotter splitting)
    # ------------------------------------------------------------------
    a = rand(1:fw.n_jumps)
    per_op = fw.per_operator[a]
    jump = fw.jumps[a]

    # ------------------------------------------------------------------
    # Apply per-operator coherent unitary FIRST
    # ------------------------------------------------------------------
    if per_op.U_B !== nothing
        mul!(ws.psi_tmp, per_op.U_B, psi)
        copyto!(psi, ws.psi_tmp)
        rmul!(psi, 1.0 / sqrt(max(_norm2(psi), eps(Float64))))
    end

    # ------------------------------------------------------------------
    # Compute probabilities from per-operator Kraus data
    # ------------------------------------------------------------------
    mul!(ws.Rpsi, per_op.R, psi)
    expR = max(real(dot(psi, ws.Rpsi)), 0.0)

    mul!(ws.psi_tmp, per_op.K0, psi)
    p_nojump = _norm2(ws.psi_tmp)

    p_jump_total = delta_eff * expR

    mul!(ws.Rpsi, per_op.U_residual, psi)
    p_res = _norm2(ws.Rpsi)

    total_weight = p_nojump + p_res + p_jump_total
    total_weight = max(total_weight, 0.0)

    if abs(total_weight - 1.0) > 1e-6
        @warn "Normalization violation: sum = $(round(total_weight; digits=6))"
    end

    r = rand() * total_weight

    # ------------------------------------------------------------------
    # Branch: no-jump / residual / dissipative jump (for selected operator a only)
    # ------------------------------------------------------------------
    if r < p_nojump
        copyto!(psi, ws.psi_tmp)
        rmul!(psi, 1.0 / sqrt(max(p_nojump, eps(Float64))))

    elseif r < (p_nojump + p_res)
        copyto!(psi, ws.Rpsi)
        rmul!(psi, 1.0 / sqrt(max(p_res, eps(Float64))))

    else
        target = r - p_nojump - p_res
        csum   = 0.0
        chosen = false
        last_norm2 = 0.0

        ham = fw.ham_or_trott
        @assert ham isa HamHam

        @inbounds if jump.hermitian
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)

                oft!(ws.jump_oft, jump, w, ham, cfg.sigma)

                mul!(ws.Rpsi, ws.jump_oft, psi)
                n2 = _norm2(ws.Rpsi)
                last_norm2 = n2
                p = delta_eff * (scaled_prefactor * transition(w)) * n2
                csum += p
                if csum >= target
                    copyto!(psi, ws.Rpsi)
                    rmul!(psi, 1.0 / sqrt(max(n2, eps(Float64))))
                    chosen = true
                    break
                end

                if w > 1e-12
                    mul!(ws.Rpsi, ws.jump_oft', psi)
                    n2 = _norm2(ws.Rpsi)
                    last_norm2 = n2
                    p = delta_eff * (scaled_prefactor * transition(-w)) * n2
                    csum += p
                    if csum >= target
                        copyto!(psi, ws.Rpsi)
                        rmul!(psi, 1.0 / sqrt(max(n2, eps(Float64))))
                        chosen = true
                        break
                    end
                end
            end
        else
            for w in energy_labels
                oft!(ws.jump_oft, jump, w, ham, cfg.sigma)

                mul!(ws.Rpsi, ws.jump_oft, psi)
                n2 = _norm2(ws.Rpsi)
                last_norm2 = n2
                p = delta_eff * (scaled_prefactor * transition(w)) * n2
                csum += p
                if csum >= target
                    copyto!(psi, ws.Rpsi)
                    rmul!(psi, 1.0 / sqrt(max(n2, eps(Float64))))
                    chosen = true
                    break
                end
            end
        end

        if !chosen
            copyto!(psi, ws.Rpsi)
            rmul!(psi, 1.0 / sqrt(max(last_norm2, eps(Float64))))
        end
    end
    return nothing
end