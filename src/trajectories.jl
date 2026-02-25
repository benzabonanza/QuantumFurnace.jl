"""
Run-time cache for trajectory simulations
"""
struct TrajectoryWorkspace{T}
    jump_oft::Matrix{T}   # buffer for A_ω (or its dagger-partner is handled by swapping mul order)
    psi_tmp::Vector{T}    # generic tmp for matvec results (K0ψ, Aωψ, Uresψ, etc.)
    Rpsi::Vector{T}       # cache Rψ so K0ψ can be formed as ψ - α(Rψ) without a second matvec
    rho_acc::Matrix{T}    # density matrix accumulator for trajectory averaging
end

function TrajectoryWorkspace(::Type{T}, dim::Int) where {T}
    return TrajectoryWorkspace{T}(
        zeros(T, dim, dim),  # jump_oft
        zeros(T, dim),       # psi_tmp
        zeros(T, dim),       # Rpsi
        zeros(T, dim, dim),  # rho_acc
    )
end

"""
Result of a trajectory simulation run.
"""
struct TrajectoryResult{T}
    rho_mean::Matrix{T}
    n_trajectories::Int
    seed::Int
    times::Union{Nothing, Vector{Float64}}
    measurements_mean::Union{Nothing, Matrix{Float64}}
    convergence::Union{Nothing, ConvergenceData}
end

"""
Result of an observable-only trajectory run.

Unlike `TrajectoryResult`, density matrix reconstruction is optional (`rho_mean` may be `nothing`).
"""
struct ObservableTrajectoryResult{T}
    times::Vector{Float64}
    measurements_mean::Matrix{Float64}  # n_obs x n_saves
    n_trajectories::Int
    seed::Int
    rho_mean::Union{Nothing, Matrix{T}}  # nothing when reconstruct_dm=false

    # Inner constructor: explicit T (used when rho_mean is a Matrix{T})
    function ObservableTrajectoryResult{T}(
        times, measurements_mean, n_trajectories, seed, rho_mean
    ) where {T}
        new{T}(times, measurements_mean, n_trajectories, seed, rho_mean)
    end
end

# Outer constructor: infer T from rho_mean when it is a matrix
function ObservableTrajectoryResult(
    times::Vector{Float64}, measurements_mean::Matrix{Float64},
    n_trajectories::Int, seed::Int, rho_mean::Matrix{T},
) where {T}
    ObservableTrajectoryResult{T}(times, measurements_mean, n_trajectories, seed, rho_mean)
end

# Outer constructor: default T=ComplexF64 when rho_mean=nothing
function ObservableTrajectoryResult(
    times::Vector{Float64}, measurements_mean::Matrix{Float64},
    n_trajectories::Int, seed::Int, rho_mean::Nothing,
)
    ObservableTrajectoryResult{ComplexF64}(times, measurements_mean, n_trajectories, seed, nothing)
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

struct TrajectoryFramework{T,D<:AbstractDomain,F,P}
    domain::D
    jumps::Vector{JumpOp{Matrix{T}}}
    ham_or_trott::Union{HamHam, TrottTrott}
    config::Config{Thermalize, D}
    precomputed_data::Any  # NamedTuple from precompute_data, varies by domain

    # Per-operator Kraus data (Lie-Trotter splitting)
    per_operator::Vector{PerOperatorKraus{T}}
    n_jumps::Int

    # Step parameters
    delta::Float64           # original delta (for time stepping / num_steps)
    delta_eff::Float64       # CPTP channel time parameter (equals delta; R scaling handles 1/p_jump compensation)
    alpha::Float64           # α = 1 - sqrt(1-δ)

    # Hot-path fields with concrete types (avoid accessing abstract-typed config/precomputed_data in step loop)
    scaled_prefactor::Float64   # rate prefactor with 1/p_jump rescaling, domain-specific
    sigma::Float64              # cfg.sigma, needed by EnergyDomain oft!() call
    transition::F               # transition function (concrete closure type)
    energy_labels::Vector{Float64}
    oft_nufft_prefactors::P     # NUFFTPrefactors or Nothing (EnergyDomain uses oft! instead)
end

function TrajectoryWorkspace(fw::TrajectoryFramework{T}) where {T}
    dim = size(fw.per_operator[1].R, 1)
    TrajectoryWorkspace(T, dim)
end

function build_trajectoryframework(
    jumps::AbstractVector{<:JumpOp},
    ham_or_trott::Union{HamHam, TrottTrott},
    config::Config{Thermalize},
    precomputed_data,
    scratch::KrausScratch{<:Complex},
    delta::Real)

    CT = eltype(scratch.R)
    dim = size(jumps[1].data, 1)
    n_jumps = length(jumps)
    p_jump = 1.0 / n_jumps

    # The per-operator CPTP channel uses bare delta (NOT delta*n_jumps).
    # R_a is already scaled by 1/p_jump = n_jumps to compensate for random
    # operator selection. Using delta*n_jumps would double-count.
    # (The DM _finalize_kraus_step! uses the same approach: scaled R, bare delta.)
    @assert delta < 1.0 "delta = $(delta) >= 1.0: too large for CPTP channel"

    alpha = 1 - sqrt(1 - delta)

    # Convert to concrete element type for zero-allocation access in hot loop
    # (jumps arrive in the correct basis: trotter.eigvecs for TrotterDomain,
    #  hamiltonian.eigvecs for other domains -- basis selection is at the source)
    jumps_for_diss = convert(Vector{JumpOp{Matrix{CT}}}, collect(JumpOp, jumps))

    # Precompute per-operator coherent B terms (one per jump)
    # Jumps are already in the correct basis (trotter.eigvecs for TrotterDomain,
    # hamiltonian.eigvecs for other domains), so pass them directly.
    per_op_U_B = Vector{Union{Nothing, Matrix{CT}}}(undef, n_jumps)
    if with_coherent(config.construction)
        @inbounds for a in 1:n_jumps
            single_jump = JumpOp[jumps[a]]  # Force Vector{JumpOp} for dispatch compatibility
            B_a = _precompute_coherent_total_B(single_jump, ham_or_trott, config, precomputed_data)
            hermitianize!(B_a)
            # Coherent uses delta/p_jump = delta*n_jumps (matching DM coherent_unitaries scaling)
            per_op_U_B[a] = exp(-1im * (delta / p_jump) * Hermitian(B_a))
        end
    else
        fill!(per_op_U_B, nothing)
    end

    # Build per-operator Kraus data
    per_operator = Vector{PerOperatorKraus{CT}}(undef, n_jumps)

    @inbounds for a in 1:n_jumps
        # Compute R^a for single operator, rescaled by 1/p_jump
        _precompute_R([jumps_for_diss[a]], ham_or_trott, config, precomputed_data, scratch)
        R_a = copy(scratch.R)
        R_a .*= (1.0 / p_jump)   # rescale: R_a = (1/p_jump) * sum_w rate2(w) * A_w' * A_w

        # K0_a = I - alpha * R_a
        K0_a = copy(R_a)
        K0_a .*= (-alpha)
        @inbounds for i in 1:dim
            K0_a[i,i] += 1
        end

        # S_a = (2*alpha - delta)*R_a - alpha^2 * R_a^2
        mul!(scratch.tmp1, R_a, R_a)   # tmp1 := R_a^2
        s1 = 2 * alpha - delta
        s2 = alpha * alpha
        @. scratch.tmp2 = s1 * R_a - s2 * scratch.tmp1

        # TFIX-04: PSD guard -- clamp negative eigenvalues to zero (silent fallback)
        hermitianize!(scratch.tmp2)
        S_herm = Hermitian(scratch.tmp2)
        eig = eigen(S_herm)
        eig.values .= max.(eig.values, 0.0)
        U_residual_a = Matrix{CT}(Diagonal(sqrt.(eig.values)) * eig.vectors')

        per_operator[a] = PerOperatorKraus(R_a, K0_a, U_residual_a, per_op_U_B[a])
    end

    # Precompute scaled_prefactor for the hot path (avoids accessing abstract config/precomputed_data in step loop)
    gamma_norm_factor = precomputed_data.gamma_norm_factor
    scaled_prefactor = if config.domain isa EnergyDomain
        # EnergyDomain formula
        config.w0 / (config.sigma * sqrt(2 * pi)) * gamma_norm_factor / (1.0 / n_jumps)
    else
        # TimeDomain / TrotterDomain formula
        config.w0 * config.t0^2 * (config.sigma * sqrt(2 / pi)) / (2 * pi) * gamma_norm_factor / (1.0 / n_jumps)
    end

    # Extract hot-path fields from precomputed_data with concrete types
    transition_fn = precomputed_data.transition
    energy_labels_vec = Vector{Float64}(precomputed_data.energy_labels)
    # OFT NUFFT prefactors: present for Time/Trotter domains, absent for EnergyDomain
    oft_nufft_pref = if hasproperty(precomputed_data, :oft_nufft_prefactors)
        precomputed_data.oft_nufft_prefactors
    else
        nothing
    end

    return TrajectoryFramework(
        config.domain,
        jumps_for_diss,
        ham_or_trott,
        config,
        precomputed_data,
        per_operator,
        n_jumps,
        Float64(delta),
        Float64(delta),       # delta_eff field = delta (R scaling handles 1/p_jump)
        Float64(alpha),
        Float64(scaled_prefactor),
        Float64(config.sigma),
        transition_fn,
        energy_labels_vec,
        oft_nufft_pref,
    )
end

"""
    precompute_R(jumps, ham_or_trott, config, precomputed_data, scratch) -> Matrix{<:Complex}

    Compute
        R = ∑_{k>0} L_k† L_k
    in the same basis as `jump.in_eigenbasis` (Hamiltonian eigenbasis for `HamHam`, Trotter basis for `TrottTrott`).

    Conventions are matched to `jump_contribution!(domain, ::Config{Thermalize}, ...)`:
    - the weights are `rate2(ω) = base_prefactor * transition(ω)` (no extra `δ` factor),
    - for Hermitian jumps we iterate half-grid and add the mirrored negative-frequency partner explicitly.

    This returns `scratch.R` (Hermitianized).
"""
function _precompute_R(
    jumps::AbstractVector{<:JumpOp},
    hamiltonian::HamHam,
    config::Config{Thermalize, EnergyDomain},
    precomputed_data,
    scratch::KrausScratch{<:Complex},
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


function _precompute_R(
    jumps::AbstractVector{<:JumpOp},
    ham_or_trott::Union{HamHam, TrottTrott},
    config::Config{Thermalize, D},
    precomputed_data,
    scratch::KrausScratch{<:Complex},
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

                nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)

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
                nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)
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
@inline _norm2(v::AbstractVector{<:Complex}) = real(dot(v, v))

@inline function _accumulate_density_matrix!(
    rho_acc::Matrix{<:Complex},
    psi::Vector{<:Complex},
)
    # rho_acc += psi * psi'   (conjugate transpose, rank-1 update)
    CT = eltype(rho_acc)
    mul!(rho_acc, psi, psi', one(CT), one(CT))
    return nothing
end

# Allocation-free measure-add (writes into `acc[:, save_idx]` by +=)
function _accumulate_measurements!(
    acc::AbstractMatrix{<:Real},
    save_idx::Int,
    psi::Vector{<:Complex},
    observables::Vector{<:Matrix{<:Complex}},
    tmp::Vector{<:Complex},
)
    @inbounds for i in eachindex(observables)
        mul!(tmp, observables[i], psi)                 # tmp := O_i ψ
        acc[i, save_idx] += real(dot(psi, tmp))        # += <ψ|O_i|ψ>
    end
    return nothing
end

"""
    _partition_trajectories(range, n_chunks) -> Vector{UnitRange{Int}}

Partition a range into n_chunks approximately equal chunks. Earlier chunks get the remainder.
"""
function _partition_trajectories(range::UnitRange{Int}, n_chunks::Int)
    len = length(range)
    n_chunks = min(n_chunks, len)  # no more chunks than items
    base = div(len, n_chunks)
    remainder = rem(len, n_chunks)
    chunks = Vector{UnitRange{Int}}(undef, n_chunks)
    start = first(range)
    for i in 1:n_chunks
        chunk_size = base + (i <= remainder ? 1 : 0)
        chunks[i] = start:(start + chunk_size - 1)
        start += chunk_size
    end
    return chunks
end

"""
    _run_chunk_no_obs!(ws, fw, psi0, chunk, master_seed, total_time)

Run a chunk of trajectories without observables, accumulating density matrices in ws.rho_acc.
Each trajectory gets Xoshiro(master_seed + traj_id) for reproducibility.
Step loop is inlined (no intermediate _evolve_along_trajectory! wrapper).
"""
function _run_chunk_no_obs!(
    ws::TrajectoryWorkspace{<:Complex},
    fw::TrajectoryFramework{<:Complex},
    psi0::Vector{<:Complex},
    chunk::UnitRange{Int},
    master_seed::Int,
    total_time::Real,
)
    delta = fw.delta
    num_steps = ceil(Int, total_time / delta)
    psi = copy(psi0)
    for traj_id in chunk
        rng = Random.Xoshiro(master_seed + traj_id)
        copyto!(psi, psi0)
        # Normalize once per trajectory
        psi_norm2 = _norm2(psi)
        rmul!(psi, 1.0 / sqrt(max(psi_norm2, eps(Float64))))
        # Step loop (was _evolve_along_trajectory!)
        @inbounds for _ in 1:num_steps
            step_along_trajectory!(psi, fw, ws, rng)
        end
        _accumulate_density_matrix!(ws.rho_acc, psi)
    end
    return nothing
end

"""
    _run_chunk_with_obs!(ws, fw, psi0, chunk, master_seed, total_time,
                          observables, save_every, num_steps, num_saves, mean_data_local)

Run a chunk of trajectories with observable measurements, accumulating density matrices
in ws.rho_acc and measurements in mean_data_local.
"""
function _run_chunk_with_obs!(
    ws::TrajectoryWorkspace{<:Complex},
    fw::TrajectoryFramework{<:Complex},
    psi0::Vector{<:Complex},
    chunk::UnitRange{Int},
    master_seed::Int,
    total_time::Real,
    observables::Vector{<:Matrix{<:Complex}},
    save_every::Int,
    num_steps::Int,
    num_saves::Int,
    mean_data_local::Matrix{Float64},
)
    psi = copy(psi0)
    tmp_meas = ws.psi_tmp  # reuse workspace vector as gemv buffer

    for traj_id in chunk
        rng = Random.Xoshiro(master_seed + traj_id)
        copyto!(psi, psi0)

        # normalize once per trajectory
        n2 = real(dot(psi, psi))
        rmul!(psi, 1.0 / sqrt(max(n2, eps(Float64))))

        # save at t=0
        _accumulate_measurements!(mean_data_local, 1, psi, observables, tmp_meas)

        save_idx = 1
        for step in 1:num_steps
            step_along_trajectory!(psi, fw, ws, rng)
            if step % save_every == 0
                save_idx += 1
                _accumulate_measurements!(mean_data_local, save_idx, psi, observables, tmp_meas)
            end
        end
        _accumulate_density_matrix!(ws.rho_acc, psi)
    end
    return nothing
end

"""
    _run_chunk_obs_only!(ws, fw, psi0, chunk, master_seed, total_time,
                          observables, save_every, num_steps, num_saves, mean_data_local)

Run a chunk of trajectories with observable measurements but WITHOUT density matrix accumulation.
Each trajectory gets Xoshiro(master_seed + traj_id) for reproducibility.
"""
function _run_chunk_obs_only!(
    ws::TrajectoryWorkspace{<:Complex},
    fw::TrajectoryFramework{<:Complex},
    psi0::Vector{<:Complex},
    chunk::UnitRange{Int},
    master_seed::Int,
    total_time::Real,
    observables::Vector{<:Matrix{<:Complex}},
    save_every::Int,
    num_steps::Int,
    num_saves::Int,
    mean_data_local::Matrix{Float64},
)
    psi = copy(psi0)
    tmp_meas = ws.psi_tmp  # reuse workspace vector as gemv buffer

    for traj_id in chunk
        rng = Random.Xoshiro(master_seed + traj_id)
        copyto!(psi, psi0)

        # normalize once per trajectory
        n2 = real(dot(psi, psi))
        rmul!(psi, 1.0 / sqrt(max(n2, eps(Float64))))

        # save at t=0
        _accumulate_measurements!(mean_data_local, 1, psi, observables, tmp_meas)

        save_idx = 1
        for step in 1:num_steps
            step_along_trajectory!(psi, fw, ws, rng)
            if step % save_every == 0
                save_idx += 1
                _accumulate_measurements!(mean_data_local, save_idx, psi, observables, tmp_meas)
            end
        end
        # NO _accumulate_density_matrix! call here
    end
    return nothing
end

"""
    _run_batch_no_obs!(fw, psi0, ntraj, master_seed, total_time) -> Matrix{CT}

Run `ntraj` trajectories using pre-built `fw`, returning the averaged density matrix.
Handles serial vs multi-threaded dispatch internally. The master_seed is the base
seed; each trajectory gets Xoshiro(master_seed + traj_id) where traj_id is 1:ntraj.

This is the shared batch execution function used by `run_trajectories` (no-observables
path), `run_trajectories_convergence`, and `run_trajectories_adaptive`.
"""
function _run_batch_no_obs!(
    fw::TrajectoryFramework{<:Complex},
    psi0::Vector{<:Complex},
    ntraj::Int,
    master_seed::Int,
    total_time::Real,
)
    CT = eltype(psi0)
    dim = length(psi0)

    if ntraj > 1 && Threads.nthreads() > 1
        # Multi-threaded path
        nt = min(Threads.nthreads(), ntraj)
        chunks = _partition_trajectories(1:ntraj, nt)
        ws_per_task = [TrajectoryWorkspace(CT, dim) for _ in 1:length(chunks)]

        old_blas = BLAS.get_num_threads()
        BLAS.set_num_threads(1)
        try
            @sync for (idx, chunk) in enumerate(chunks)
                Threads.@spawn _run_chunk_no_obs!(
                    ws_per_task[idx], fw, psi0, chunk, master_seed, total_time)
            end
        finally
            BLAS.set_num_threads(old_blas)
        end

        rho_total = sum(ws.rho_acc for ws in ws_per_task)
        rho_result = rho_total ./ ntraj
        hermitianize!(rho_result)
    else
        # Serial path
        ws = TrajectoryWorkspace(CT, dim)
        _run_chunk_no_obs!(ws, fw, psi0, 1:ntraj, master_seed, total_time)
        rho_result = ws.rho_acc ./ ntraj
        hermitianize!(rho_result)
    end

    return rho_result
end

"""
    _build_framework_and_seed(jumps, config, psi0, hamiltonian; trotter, delta, seed)

One-time setup: validates config, chooses ham_or_trott, precomputes data, builds
TrajectoryFramework, and generates actual seed. Returns `(fw, actual_seed)`.

This is extracted so convergence/adaptive runners can build the framework ONCE
and reuse it across batches.
"""
function _build_framework_and_seed(
    jumps::Vector{JumpOp},
    config::Config{Thermalize},
    psi0::Vector{<:Complex},
    hamiltonian::HamHam;
    trotter::Union{TrottTrott,Nothing}=nothing,
    delta::Real = config.delta,
    seed::Union{Int,Nothing} = nothing,
)
    validate_config!(config)
    _print_press(config)

    CT = eltype(psi0)

    # Choose evolution object consistent with domain
    ham_or_trott = if config.domain isa TrotterDomain
        trotter === nothing && error("TrotterDomain requires `trotter`.")
        trotter
    else
        hamiltonian
    end

    precomputed_data = _precompute_data(config, ham_or_trott)

    dim = size(hamiltonian.data, 1)
    builder_scratch = KrausScratch(CT, dim)

    fw = build_trajectoryframework(jumps, ham_or_trott, config, precomputed_data, builder_scratch, delta)

    actual_seed = seed === nothing ? Int(rand(Random.RandomDevice(), UInt64) >> 1) : seed

    return fw, actual_seed
end

"""
    run_trajectories(jumps, config, psi0, hamiltonian; trotter=nothing,
                     total_time=config.mixing_time, delta=config.delta,
                     ntraj=1, observables=nothing, save_every=1, seed=nothing)

    Builds the TrajectoryFramework once, then runs `ntraj` trajectories.
    Returns a `TrajectoryResult` containing the averaged density matrix, trajectory count,
    and the RNG seed used (for reproducibility).

    - If `seed` is `nothing`, a random seed is generated from system entropy.
    - If `observables === nothing`: returns `TrajectoryResult` with `times=nothing`, `measurements_mean=nothing`.
    - If `observables` provided: returns `TrajectoryResult` with time grid and trajectory-averaged ⟨O_i⟩(t).
"""
function run_trajectories(
    jumps::Vector{JumpOp},
    config::Config{Thermalize},
    psi0::Vector{<:Complex},
    hamiltonian::HamHam;
    trotter::Union{TrottTrott,Nothing}=nothing,
    total_time::Real = config.mixing_time,
    delta::Real = config.delta,
    ntraj::Int = 1,
    observables::Union{Nothing, Vector{<:Matrix{<:Complex}}} = nothing,
    save_every::Int = 1,
    seed::Union{Int,Nothing} = nothing,
)

    @assert ntraj >= 1
    @assert save_every >= 1

    fw, actual_seed = _build_framework_and_seed(
        jumps, config, psi0, hamiltonian;
        trotter=trotter, delta=delta, seed=seed,
    )

    CT = eltype(psi0)
    dim = size(hamiltonian.data, 1)

    # ------------------------------------------------------------------
    # No measurements: use shared _run_batch_no_obs!
    # ------------------------------------------------------------------
    if observables === nothing
        rho_result = _run_batch_no_obs!(fw, psi0, ntraj, actual_seed, total_time)
        return TrajectoryResult(rho_result, ntraj, actual_seed, nothing, nothing, nothing)
    end

    # ------------------------------------------------------------------
    # With measurements: average <O> over trajectories, saved every `save_every`
    # ------------------------------------------------------------------
    ws = TrajectoryWorkspace(CT, dim)
    delta_step = fw.delta
    num_steps = ceil(Int, total_time / delta_step)
    num_saves = div(num_steps, save_every) + 1
    num_obs   = length(observables)

    times = Vector{Float64}(undef, num_saves)
    @inbounds for s in 1:num_saves
        times[s] = (s - 1) * save_every * delta_step
    end

    if ntraj > 1 && Threads.nthreads() > 1
        # Multi-threaded observable path
        nt = min(Threads.nthreads(), ntraj)
        chunks = _partition_trajectories(1:ntraj, nt)
        ws_per_task = [TrajectoryWorkspace(CT, dim) for _ in 1:length(chunks)]
        mean_data_per_task = [zeros(Float64, num_obs, num_saves) for _ in 1:length(chunks)]

        old_blas = BLAS.get_num_threads()
        BLAS.set_num_threads(1)
        try
            @sync for (idx, chunk) in enumerate(chunks)
                Threads.@spawn _run_chunk_with_obs!(
                    ws_per_task[idx], fw, psi0, chunk, actual_seed, total_time,
                    observables, save_every, num_steps, num_saves, mean_data_per_task[idx])
            end
        finally
            BLAS.set_num_threads(old_blas)
        end

        mean_data = sum(mean_data_per_task)
        mean_data ./= ntraj
        rho_total = sum(ws.rho_acc for ws in ws_per_task)
        rho_result = rho_total ./ ntraj
        hermitianize!(rho_result)
    else
        # Serial observable path
        mean_data = zeros(Float64, num_obs, num_saves)
        psi = copy(psi0)

        for traj_id in 1:ntraj
            rng_serial = Random.Xoshiro(actual_seed + traj_id)
            copyto!(psi, psi0)

            n2 = real(dot(psi, psi))
            rmul!(psi, 1.0 / sqrt(max(n2, eps(Float64))))

            _accumulate_measurements!(mean_data, 1, psi, observables, ws.psi_tmp)

            save_idx = 1
            for step in 1:num_steps
                step_along_trajectory!(psi, fw, ws, rng_serial)
                if step % save_every == 0
                    save_idx += 1
                    _accumulate_measurements!(mean_data, save_idx, psi, observables, ws.psi_tmp)
                end
            end
            _accumulate_density_matrix!(ws.rho_acc, psi)
        end

        mean_data ./= ntraj
        rho_result = ws.rho_acc ./ ntraj
        hermitianize!(rho_result)
    end

    return TrajectoryResult(rho_result, ntraj, actual_seed, times, mean_data, nothing)
end

"""
    run_observable_trajectories(jumps, config, psi0, hamiltonian;
        trotter=nothing, total_time=config.mixing_time, delta=config.delta,
        ntraj=1, observables, save_every=1, seed=nothing, reconstruct_dm=false)

Run trajectory simulations measuring time-resolved observables without per-trajectory
density matrix reconstruction. When `reconstruct_dm=true`, also accumulates the averaged
density matrix (using the existing `_run_chunk_with_obs!` path).

Returns an `ObservableTrajectoryResult` with `rho_mean=nothing` when `reconstruct_dm=false`.
"""
function run_observable_trajectories(
    jumps::Vector{JumpOp},
    config::Config{Thermalize},
    psi0::Vector{<:Complex},
    hamiltonian::HamHam;
    trotter::Union{TrottTrott,Nothing}=nothing,
    total_time::Real = config.mixing_time,
    delta::Real = config.delta,
    ntraj::Int = 1,
    observables::Vector{<:Matrix{<:Complex}},
    save_every::Int = 1,
    seed::Union{Int,Nothing} = nothing,
    reconstruct_dm::Bool = false,
)
    @assert ntraj >= 1
    @assert save_every >= 1

    fw, actual_seed = _build_framework_and_seed(
        jumps, config, psi0, hamiltonian;
        trotter=trotter, delta=delta, seed=seed,
    )

    CT = eltype(psi0)
    dim = size(hamiltonian.data, 1)

    delta_step = fw.delta
    num_steps = ceil(Int, total_time / delta_step)
    num_saves = div(num_steps, save_every) + 1
    num_obs   = length(observables)

    times = Vector{Float64}(undef, num_saves)
    @inbounds for s in 1:num_saves
        times[s] = (s - 1) * save_every * delta_step
    end

    if ntraj > 1 && Threads.nthreads() > 1
        # Multi-threaded path
        nt = min(Threads.nthreads(), ntraj)
        chunks = _partition_trajectories(1:ntraj, nt)
        ws_per_task = [TrajectoryWorkspace(CT, dim) for _ in 1:length(chunks)]
        mean_data_per_task = [zeros(Float64, num_obs, num_saves) for _ in 1:length(chunks)]

        old_blas = BLAS.get_num_threads()
        BLAS.set_num_threads(1)
        try
            if reconstruct_dm
                @sync for (idx, chunk) in enumerate(chunks)
                    Threads.@spawn _run_chunk_with_obs!(
                        ws_per_task[idx], fw, psi0, chunk, actual_seed, total_time,
                        observables, save_every, num_steps, num_saves, mean_data_per_task[idx])
                end
            else
                @sync for (idx, chunk) in enumerate(chunks)
                    Threads.@spawn _run_chunk_obs_only!(
                        ws_per_task[idx], fw, psi0, chunk, actual_seed, total_time,
                        observables, save_every, num_steps, num_saves, mean_data_per_task[idx])
                end
            end
        finally
            BLAS.set_num_threads(old_blas)
        end

        mean_data = sum(mean_data_per_task)
        mean_data ./= ntraj

        if reconstruct_dm
            rho_total = sum(ws.rho_acc for ws in ws_per_task)
            rho_result = rho_total ./ ntraj
            hermitianize!(rho_result)
            rho_mean = rho_result
        else
            rho_mean = nothing
        end
    else
        # Serial path
        ws = TrajectoryWorkspace(CT, dim)
        mean_data = zeros(Float64, num_obs, num_saves)
        psi = copy(psi0)

        for traj_id in 1:ntraj
            rng_serial = Random.Xoshiro(actual_seed + traj_id)
            copyto!(psi, psi0)

            n2 = real(dot(psi, psi))
            rmul!(psi, 1.0 / sqrt(max(n2, eps(Float64))))

            _accumulate_measurements!(mean_data, 1, psi, observables, ws.psi_tmp)

            save_idx = 1
            for step in 1:num_steps
                step_along_trajectory!(psi, fw, ws, rng_serial)
                if step % save_every == 0
                    save_idx += 1
                    _accumulate_measurements!(mean_data, save_idx, psi, observables, ws.psi_tmp)
                end
            end
            if reconstruct_dm
                _accumulate_density_matrix!(ws.rho_acc, psi)
            end
        end

        mean_data ./= ntraj

        if reconstruct_dm
            rho_result = ws.rho_acc ./ ntraj
            hermitianize!(rho_result)
            rho_mean = rho_result
        else
            rho_mean = nothing
        end
    end

    return ObservableTrajectoryResult{CT}(times, mean_data, ntraj, actual_seed, rho_mean)
end

"""
    step_along_trajectory!(psi, fw, ws, rng)  for D ∈ {TimeDomain, TrotterDomain}

    Per-operator Lie-Trotter splitting: randomly select ONE operator a, apply that operator's
    CPTP channel (K0_a, U_res_a, jump outcomes for operator a only).

    Arguments:
    - `psi`: state vector (modified in-place)
    - `fw`: read-only trajectory framework
    - `ws`: mutable workspace (scratch buffers)
    - `rng`: random number generator (explicit for thread safety and reproducibility)

    # Per-operator channel structure (Chen 2023, adapted for Lie-Trotter splitting):
    #   Pick a ∈ {1,...,N_jumps} uniformly at random
    #   K0_a = I - alpha*R_a, where alpha = 1 - sqrt(1-delta), R_a scaled by n_jumps
    #   K_{a,w} = sqrt(delta * scaled_rate(w)) * L_{a,w}  (jump operators for operator a)
    #   U_res_a: U_res_a'*U_res_a = S_a  (residual for operator a)
    # R_a rates rescaled by 1/p_jump; CPTP channel uses bare delta (matching DM).
"""
function step_along_trajectory!(
    psi::Vector{<:Complex},
    fw::TrajectoryFramework{<:Complex,D},
    ws::TrajectoryWorkspace{<:Complex},
    rng::AbstractRNG,
    ) where {D<:Union{TimeDomain,TrotterDomain}}

    # All hot-path data lives in concrete-typed fields of fw (no abstract access)
    delta = fw.delta
    scaled_prefactor = fw.scaled_prefactor
    transition = fw.transition
    energy_labels = fw.energy_labels
    oft_prefactors = fw.oft_nufft_prefactors

    # Select random operator (Vector elements are concrete-typed)
    a = rand(rng, 1:fw.n_jumps)
    @inbounds per_op = fw.per_operator[a]
    @inbounds jump = fw.jumps[a]

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

    p_jump_total = delta * expR                     # bare delta (R_a already scaled by n_jumps)

    mul!(ws.Rpsi, per_op.U_residual, psi)           # U_res_a * psi (reuse buffer)
    p_res = _norm2(ws.Rpsi)

    total_weight = p_nojump + p_res + p_jump_total
    total_weight = max(total_weight, 0.0)

    if abs(total_weight - 1.0) > 1e-6
        @warn "Normalization violation: sum = $(round(total_weight; digits=6))"
    end

    r = rand(rng) * total_weight

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

                pref = _prefactor_view(oft_prefactors, w)
                @. ws.jump_oft = jump.in_eigenbasis * pref

                # Positive-frequency
                mul!(ws.Rpsi, ws.jump_oft, psi)
                n2 = _norm2(ws.Rpsi)
                last_norm2 = n2
                p = delta * (scaled_prefactor * transition(w)) * n2
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
                    p = delta * (scaled_prefactor * transition(-w)) * n2
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
                pref = _prefactor_view(oft_prefactors, w)
                @. ws.jump_oft = jump.in_eigenbasis * pref

                mul!(ws.Rpsi, ws.jump_oft, psi)
                n2 = _norm2(ws.Rpsi)
                last_norm2 = n2
                p = delta * (scaled_prefactor * transition(w)) * n2
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
    step_along_trajectory!(psi, fw, ws, rng)  for D == EnergyDomain

    Per-operator Lie-Trotter splitting for EnergyDomain.
    Same logic as Time/Trotter variant but A_{a,ω} is generated by `oft!(...)`.
"""
function step_along_trajectory!(
    psi::Vector{<:Complex},
    fw::TrajectoryFramework{<:Complex,EnergyDomain},
    ws::TrajectoryWorkspace{<:Complex},
    rng::AbstractRNG,
    )

    # All hot-path data lives in concrete-typed fields of fw (no abstract access)
    delta = fw.delta
    scaled_prefactor = fw.scaled_prefactor
    sigma = fw.sigma
    transition = fw.transition
    energy_labels = fw.energy_labels

    # Select random operator (Vector elements are concrete-typed)
    a = rand(rng, 1:fw.n_jumps)
    @inbounds per_op = fw.per_operator[a]
    @inbounds jump = fw.jumps[a]

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

    p_jump_total = delta * expR                     # bare delta (R_a already scaled by n_jumps)

    mul!(ws.Rpsi, per_op.U_residual, psi)
    p_res = _norm2(ws.Rpsi)

    total_weight = p_nojump + p_res + p_jump_total
    total_weight = max(total_weight, 0.0)

    if abs(total_weight - 1.0) > 1e-6
        @warn "Normalization violation: sum = $(round(total_weight; digits=6))"
    end

    r = rand(rng) * total_weight

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

                oft!(ws.jump_oft, jump, w, ham, sigma)

                mul!(ws.Rpsi, ws.jump_oft, psi)
                n2 = _norm2(ws.Rpsi)
                last_norm2 = n2
                p = delta * (scaled_prefactor * transition(w)) * n2
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
                    p = delta * (scaled_prefactor * transition(-w)) * n2
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
                oft!(ws.jump_oft, jump, w, ham, sigma)

                mul!(ws.Rpsi, ws.jump_oft, psi)
                n2 = _norm2(ws.Rpsi)
                last_norm2 = n2
                p = delta * (scaled_prefactor * transition(w)) * n2
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