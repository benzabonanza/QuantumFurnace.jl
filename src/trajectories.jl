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

struct TrajectoryFramework{T,C,H,PD,D<:AbstractDomain}
    domain::D
    jumps::Vector{JumpOp}
    ham_or_trott::H
    config::C
    precomputed_data::PD

    # Coherent term (optional)
    B::Union{Nothing, Matrix{T}}
    U_B::Union{Nothing, Matrix{T}}

    # Dissipative Kraus skeleton
    R::Matrix{T}
    K0::Matrix{T}
    U_residual::Matrix{T}   # Cholesky factor U with S ≈ U' U (used as residual Kraus)

    # Step parameters
    delta::Float64
    alpha::Float64           # α = 1 - sqrt(1-δ)

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

    if config.with_coherent
        B_total = precompute_coherent_total_B(jumps, ham_or_trott, config, precomputed_data)
        B_total .= 0.5 .* (B_total .+ B_total')
        U_B = exp(-1im * delta * Hermitian(B_total))
    else
        B_total = nothing
        U_B = nothing
    end
    
    R = precompute_R(config.domain, jumps, ham_or_trott, config, precomputed_data, scratch)
    R = copy(scratch.R)

    alpha = 1 - sqrt(1 - delta)
    K0 = copy(R)
    K0 .*= (-alpha)
    @inbounds for i in 1:dim
        K0[i,i] += 1
    end

    mul!(scratch.tmp1, R, R)  # tmp1 := R²
    s1 = 2 * alpha - delta
    s2 = alpha * alpha
    @. scratch.tmp2 = s1 * R - s2 * scratch.tmp1   # tmp2 := S (up to roundoff)

    # Numerical symmetrization and tiny diagonal shift (matches thermalization code logic)
    scratch.tmp2 .= 0.5 .* (scratch.tmp2 .+ scratch.tmp2')
    eps_shift = 10 * eps(Float64)
    @inbounds for i in 1:dim
        scratch.tmp2[i,i] += eps_shift
    end

    S = copy(scratch.tmp2)  # store S explicitly (post-shift)

    # Factor for applying the residual Kraus to state vectors (K_res ≈ U where S ≈ U' U)
    cholesky_S = cholesky!(Hermitian(scratch.tmp2), check=false)
    U_residual = Matrix{ComplexF64}(cholesky_S.U)  # detach from scratch.tmp2

    ws = TrajectoryWorkspace(ComplexF64, dim)

    return TrajectoryFramework(
        config.domain,
        collect(jumps),
        ham_or_trott,
        config,
        precomputed_data,
        B_total,
        U_B,
        R,
        K0,
        U_residual,
        delta,
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
    ::EnergyDomain,
    jumps::AbstractVector{<:JumpOp},
    hamiltonian::HamHam,
    config::AbstractThermalizeConfig,
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
    scratch.R .= 0.5 .* (scratch.R .+ scratch.R')
    return scratch.R
end


function precompute_R(
    ::Union{TimeDomain, TrotterDomain},
    jumps::AbstractVector{<:JumpOp},
    ham_or_trott::Union{HamHam, TrottTrott},
    config::AbstractThermalizeConfig,
    precomputed_data,
    scratch::KrausScratch{ComplexF64},
    )
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

    scratch.R .= 0.5 .* (scratch.R .+ scratch.R')
    return scratch.R
end

# Fast squared norm without sqrt
@inline _norm2(v::AbstractVector{ComplexF64}) = real(dot(v, v))

@inline function _accumulate_density_matrix!(
    rho_acc::Matrix{ComplexF64},
    psi::Vector{ComplexF64},
)
    # rho_acc += psi * psi'   (conjugate transpose)
    BLAS.gerc!(one(ComplexF64), psi, psi, rho_acc)
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

    return nothing
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

    precomputed_data = precompute_data(config.domain, config, ham_or_trott)

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
            rho_mean .= 0.5 .* (rho_mean .+ rho_mean')
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
        rho_mean .= 0.5 .* (rho_mean .+ rho_mean')

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
    rho_mean .= 0.5 .* (rho_mean .+ rho_mean')

    return (framework = fw, times = times, measurements_mean = mean_data, rho_mean = rho_mean)
end

"""
    step_along_trajectory!(psi, fw::TrajectoryFramework)  for D ∈ {TimeDomain, TrotterDomain}

    Implements one Chen-faithful δ-step unraveling of the *full* CPTP Kraus map:
    - no-jump Kraus:      K0 = I - α R,   α = 1 - sqrt(1-δ)
    - jump Kraus family:  K_{a,ω} = √(δ rate2(ω)) A_{a,ω}   (generated on the fly via NUFFT prefactors)
    - residual Kraus:     K_res = U_residual  with  S = U_res† U_res = I - K0†K0 - δR
    Then applies the deterministic coherent unitary U_B (if present) after the Kraus update.
"""
function step_along_trajectory!(
    psi::Vector{ComplexF64},
    fw::TrajectoryFramework{ComplexF64,<:Any,<:Any,<:Any,D},
    ) where {D<:Union{TimeDomain,TrotterDomain}}

    ws = fw.ws
    cfg = fw.config
    pd  = fw.precomputed_data

    # Pull hot fields into locals (fewer dynamic dispatches / repeated getproperty)
    transition         = pd.transition
    energy_labels      = pd.energy_labels
    oft_prefactors     = pd.oft_nufft_prefactors
    gamma_norm_factor  = pd.gamma_norm_factor

    delta = fw.delta
    alpha = fw.alpha
    R = fw.R

    # Same rate prefactor used in precompute_R(::TimeDomain/::TrotterDomain) (no extra δ)
    base_prefactor = cfg.w0 * cfg.t0^2 * (cfg.sigma * sqrt(2 / pi)) / (2 * pi) * gamma_norm_factor

    # ------------------------------------------------------------------
    # Compute p_nojump and p_jump_total cheaply from one Rψ
    # ------------------------------------------------------------------
    mul!(ws.Rpsi, R, psi)                          # ws.Rpsi := Rψ
    expR = real(dot(psi, ws.Rpsi))                 # <ψ|R|ψ>, should be ≥ 0
    expR = max(expR, 0.0)

    # K0ψ = (I - αR)ψ = ψ - α(Rψ)
    copyto!(ws.psi_tmp, psi)
    @inbounds @. ws.psi_tmp = ws.psi_tmp - alpha * ws.Rpsi
    p_nojump = _norm2(ws.psi_tmp)

    # total dissipative jump probability (excluding residual completion)
    p_jump_total = delta * expR

    # ------------------------------------------------------------------
    # Residual probability p_res = ||U_residual ψ||^2
    #    Reuse ws.Rpsi buffer (Rψ no longer needed beyond expR + K0ψ)
    # ------------------------------------------------------------------
    mul!(ws.Rpsi, fw.U_residual, psi)              # ws.Rpsi := U_res ψ
    p_res = _norm2(ws.Rpsi)

    total_weight = p_nojump + p_res + p_jump_total
    total_weight = max(total_weight, 0.0)

    r = rand() * total_weight

    # ------------------------------------------------------------------
    # Deterministic coherent update (if enabled): ψ ← U_B ψ
    # ------------------------------------------------------------------
    if fw.U_B !== nothing
        mul!(ws.psi_tmp, fw.U_B, psi)
        copyto!(psi, ws.psi_tmp)
        rmul!(psi, 1.0 / sqrt(max(_norm2(psi), eps(Float64))))  # guard drift from expm / roundoff
    end

    # ------------------------------------------------------------------
    # Branch: no-jump / residual / dissipative jump
    # ------------------------------------------------------------------

    if r < p_nojump
        # ψ ← K0ψ / ||K0ψ||
        copyto!(psi, ws.psi_tmp)
        rmul!(psi, 1.0 / sqrt(max(p_nojump, eps(Float64))))

    elseif r < (p_nojump + p_res)
        # ψ ← U_res ψ / ||U_res ψ||
        copyto!(psi, ws.Rpsi)
        rmul!(psi, 1.0 / sqrt(max(p_res, eps(Float64))))

    else
        # Dissipative jump: sample one (a,ω,±) outcome by cumulative probability scan.
        target = r - p_nojump - p_res
        csum   = 0.0

        chosen = false
        last_norm2 = 0.0

        @inbounds for jump in fw.jumps
            if jump.hermitian
                # half-grid (w_raw <= 0) and explicit mirrored negative partner
                for w_raw in energy_labels
                    w_raw > 1e-12 && continue
                    w = abs(w_raw)

                    # Build A_{a,ω} = A ∘ P(ω) into the single d×d buffer.
                    pref = prefactor_view(oft_prefactors, w)
                    @. ws.jump_oft = jump.in_eigenbasis * pref

                    # Positive-frequency jump operator: A_{a,ω}
                    mul!(ws.Rpsi, ws.jump_oft, psi)              # reuse ws.Rpsi as candidate buffer
                    n2 = _norm2(ws.Rpsi)
                    last_norm2 = n2
                    p = delta * (base_prefactor * transition(w)) * n2
                    csum += p
                    if csum >= target
                        copyto!(psi, ws.Rpsi)
                        rmul!(psi, 1.0 / sqrt(max(n2, eps(Float64))))
                        chosen = true
                        break
                    end

                    # Negative-frequency partner uses operator A_{a,ω}† and weight transition(-w)
                    if w > 1e-12
                        mul!(ws.Rpsi, ws.jump_oft', psi)
                        n2 = _norm2(ws.Rpsi)
                        last_norm2 = n2
                        p = delta * (base_prefactor * transition(-w)) * n2
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
                # full grid for non-Hermitian jumps
                for w in energy_labels
                    pref = prefactor_view(oft_prefactors, w)
                    @. ws.jump_oft = jump.in_eigenbasis * pref

                    mul!(ws.Rpsi, ws.jump_oft, psi)
                    n2 = _norm2(ws.Rpsi)
                    last_norm2 = n2
                    p = delta * (base_prefactor * transition(w)) * n2
                    csum += p
                    if csum >= target
                        copyto!(psi, ws.Rpsi)
                        rmul!(psi, 1.0 / sqrt(max(n2, eps(Float64))))
                        chosen = true
                        break
                    end
                end
            end

            chosen && break
        end

        # Fallback if rounding prevents hitting target: normalize last candidate buffer.
        if !chosen
            # ws.Rpsi holds the last computed candidate state
            copyto!(psi, ws.Rpsi)
            rmul!(psi, 1.0 / sqrt(max(last_norm2, eps(Float64))))
        end
    end
    return nothing
end

"""
    step_along_trajectory!(psi, fw::TrajectoryFramework)  for D == EnergyDomain

    Same logic as the Time/Trotter variant, but A_{a,ω} is generated by `oft!(...)` (Gaussian in Bohr-frequency space),
    not via NUFFT prefactors.

    Probabilities:
    - p_nojump = ||(I-αR)ψ||^2
    - p_jump_total = δ <ψ|R|ψ>
    - p_res = ||U_res ψ||^2
    - per-outcome p(a,ω,+) = δ * rate2(ω) * ||A_{a,ω} ψ||^2
    - per-outcome p(a,ω,-) = δ * rate2(-ω) * ||A_{a,ω}† ψ||^2   (Hermitian jump case via explicit partner)
"""
function step_along_trajectory!(
    psi::Vector{ComplexF64},
    fw::TrajectoryFramework{ComplexF64,<:Any,<:Any,<:Any,EnergyDomain},
    )

    ws = fw.ws
    cfg = fw.config
    pd  = fw.precomputed_data

    transition         = pd.transition
    energy_labels      = pd.energy_labels
    gamma_norm_factor  = pd.gamma_norm_factor

    δ = fw.delta
    α = fw.alpha
    R = fw.R

    # Same base prefactor as precompute_R(::EnergyDomain) (no extra δ)
    base_prefactor = cfg.w0 / (cfg.sigma * sqrt(2 * pi)) * gamma_norm_factor

    # ------------------------------------------------------------------
    # Compute p_nojump and p_jump_total cheaply from one Rψ
    # ------------------------------------------------------------------
    mul!(ws.Rpsi, R, psi)                          # Rψ
    expR = real(dot(psi, ws.Rpsi))
    expR = max(expR, 0.0)

    # K0ψ = ψ - α Rψ
    copyto!(ws.psi_tmp, psi)
    @inbounds @. ws.psi_tmp = ws.psi_tmp - α * ws.Rpsi
    p_nojump = _norm2(ws.psi_tmp)

    p_jump_total = δ * expR

    # ------------------------------------------------------------------
    # Residual probability
    # ------------------------------------------------------------------
    mul!(ws.Rpsi, fw.U_residual, psi)              # U_res ψ
    p_res = _norm2(ws.Rpsi)

    total_weight = p_nojump + p_res + p_jump_total
    total_weight = max(total_weight, 0.0)

    r = rand() * total_weight

    # ------------------------------------------------------------------
    # Coherent unitary (if enabled)
    # ------------------------------------------------------------------
    if fw.U_B !== nothing
        mul!(ws.psi_tmp, fw.U_B, psi)
        copyto!(psi, ws.psi_tmp)
        rmul!(psi, 1.0 / sqrt(max(_norm2(psi), eps(Float64))))
    end

    # ------------------------------------------------------------------
    # Branch
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

        # Hamiltonian (or Trotter) object is needed to build A_{a,ω} in eigenbasis
        ham = fw.ham_or_trott
        @assert ham isa HamHam  # EnergyDomain should be paired with HamHam

        @inbounds for jump in fw.jumps
            if jump.hermitian
                # half-grid + explicit mirror partner
                for w_raw in energy_labels
                    w_raw > 1e-12 && continue
                    w = abs(w_raw)

                    # Build A_{a,ω} into ws.jump_oft buffer
                    oft!(ws.jump_oft, jump, w, ham, cfg.sigma)

                    # + branch: A_{a,ω}
                    mul!(ws.Rpsi, ws.jump_oft, psi)
                    n2 = _norm2(ws.Rpsi)
                    last_norm2 = n2
                    p = δ * (base_prefactor * transition(w)) * n2
                    csum += p
                    if csum >= target
                        copyto!(psi, ws.Rpsi)
                        rmul!(psi, 1.0 / sqrt(max(n2, eps(Float64))))
                        chosen = true
                        break
                    end

                    # - branch: A_{a,ω}† with transition(-w)
                    if w > 1e-12
                        mul!(ws.Rpsi, ws.jump_oft', psi)
                        n2 = _norm2(ws.Rpsi)
                        last_norm2 = n2
                        p = δ * (base_prefactor * transition(-w)) * n2
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
                # full grid
                for w in energy_labels
                    oft!(ws.jump_oft, jump, w, ham, cfg.sigma)

                    mul!(ws.Rpsi, ws.jump_oft, psi)
                    n2 = _norm2(ws.Rpsi)
                    last_norm2 = n2
                    p = δ * (base_prefactor * transition(w)) * n2
                    csum += p
                    if csum >= target
                        copyto!(psi, ws.Rpsi)
                        rmul!(psi, 1.0 / sqrt(max(n2, eps(Float64))))
                        chosen = true
                        break
                    end
                end
            end

            chosen && break
        end

        if !chosen
            copyto!(psi, ws.Rpsi)
            rmul!(psi, 1.0 / sqrt(max(last_norm2, eps(Float64))))
        end
    end
    return nothing
end

# """
#     evolve_and_measure_along_trajectory(psi0, fw::TrajectoryFramework, T, observables; save_every=1)

#     Runs a single trajectory for total time `T`, saving ⟨O_i⟩ at step 0 and then every `save_every` steps.
#     Uses the new `step_along_trajectory!` dispatch (EnergyDomain vs Time/Trotter).
#     Avoids allocations inside the loop except the `view` creation (cheap); can be removed by manual indexing.

#     Returns a NamedTuple: (psi, times, measurements).
# """
# function evolve_and_measure_along_trajectory(
#     psi0::Vector{ComplexF64},
#     fw::TrajectoryFramework{ComplexF64},
#     T::Float64,
#     observables::Vector{Matrix{ComplexF64}};
#     save_every::Int = 1,
#     )

#     δ = fw.delta
#     @assert δ > 0
#     @assert T ≥ 0
#     @assert save_every ≥ 1

#     num_steps = ceil(Int, T / δ)
#     num_saves = div(num_steps, save_every) + 1
#     num_obs   = length(observables)

#     data  = zeros(Float64, num_obs, num_saves)
#     times = zeros(Float64, num_saves)

#     psi = copy(psi0)
#     n2 = real(dot(psi, psi))
#     rmul!(psi, 1.0 / sqrt(max(n2, eps(Float64))))

#     # initial measurement
#     save_index = 1
#     times[save_index] = 0.0
#     measure!(view(data, :, save_index), psi, observables)

#     @inbounds for step in 1:num_steps
#         step_along_trajectory!(psi, fw)

#         if step % save_every == 0
#             save_index += 1
#             times[save_index] = step * δ
#             measure!(view(data, :, save_index), psi, observables)
#         end
#     end

#     return (psi = psi, times = times, measurements = data)
# end

# function measure!(measured_values, state, observable_list)
#     for i in eachindex(observable_list)
#         measured_values[i] = real(dot(state, observable_list[i], state))
#     end
# end


# struct KrausFramework{T}
#     kraus_H_eff::Matrix{T}                   # non-Hermitian no jump evolution
#     kraus_jumps::Vector{Matrix{T}}      # jump Kraus
#     R::Matrix{T}                    # sum(M^\dagger M) (without delta)
#     psi_temp::Vector{T}             # buffer for less allocations
#     delta::Float64                  # delta steps for trajectory
# end

# function krausframework(
#     H::AbstractMatrix{T}, 
#     kraus_jumps::Vector{Tuple{Float64, <:AbstractMatrix{T}}}, 
#     R::AbstractMatrix{T},
#     delta::Float64) where T

#     dim = size(H, 1)
    
#     kraus_jumps = [Matrix{T}(sqrt(delta) * rate * op) for (rate, op) in kraus_jumps]

#     # Construct H_eff and the 0th Kraus operator
#     H_eff = 1im * H + 0.5 * R
#     kraus_H_eff = Matrix{T}(I, dim, dim) - delta * H_eff

#     # Buffer for statevector
#     psi_temp = zeros(T, dim)

#     return KrausFramework(kraus_H_eff, kraus_jumps, R, psi_temp, delta)
# end


# function construct_gksl_lindbladian(
#     H::AbstractMatrix{ComplexF64}, 
#     kraus_jumps::Vector{Matrix{ComplexF64}}
#     )
#     dim = size(H, 1)

#     lindblad = zeros(ComplexF64, dim^2, dim^2)

#     # Add coherent part
#     vectorize_liouvillian_coherent!(lindblad, H)

#     # Add dissipative part
#     for kraus_jump in kraus_jumps
#         vectorize_liouv_diss_and_add!(lindblad, kraus_jump, 1.0)
#     end

#     return lindblad
# end

# function evolve_along_trajectory(psi0::Vector{ComplexF64}, fw::KrausFramework, total_time::Float64)::Vector{ComplexF64}

#     num_steps = round(Int, total_time / fw.delta)  # Number of trajectory steps

#     psi = copy(psi0)
#     for _ in 1:num_steps
#         step_along_the_trajectory!(psi, fw)
#     end

#     return psi
# end


# function step_along_the_trajectory!(psi::Vector{ComplexF64}, fw::KrausFramework)

#     # No jump
#     mul!(fw.psi_temp, fw.M0, psi)
#     prob_no_jump = norm(fw.psi_temp)^2

#     # Total jump probability
#     p_jump_total = fw.delta * real(dot(psi, fw.R, psi))

#     total_weight = prob_no_jump + p_jump_total

#     r = rand() * total_weight

#     if r < prob_no_jump
#         # No jump
#         copyto!(psi, fw.psi_temp)

#         # Force normalize
#         rmul!(psi, 1.0 / sqrt(prob_no_jump))
#     else  # Jump.
#         # Iterate through jumps and their probabilites till we find the winner
#         target_cummulative = r - prob_no_jump
#         current_cummulative = 0.0
        
#         for k in 1:length(fw.M_jumps)
#             mul!(fw.psi_temp, fw.M_jumps[k], psi)
#             prob_jump_k = norm(fw.psi_temp)^2
#             # println("Prob of jump k")
#             # println(prob_jump_k)

#             current_cummulative += prob_jump_k

#             if current_cummulative >= target_cummulative
#                 copyto!(psi, fw.psi_temp)

#                 # Force normalize
#                 rmul!(psi, 1.0 / sqrt(prob_jump_k))
#                 return
#             end
#         end

#         # If somehow we haven't picked any jumps, then use the last one.
#         copyto!(psi, fw.psi_temp)
#         normalize!(psi)
#     end
# end