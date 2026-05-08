# ============================================================================
# Krylov Eigensolver: matrix-free spectral gap computation via KrylovKit
# ============================================================================
#
# Wraps KrylovKit.eigsolve around the Krylov matvec infrastructure (Phase 27-28)
# to provide a single-call API for computing Lindbladian spectral gaps without
# constructing the full dense Liouvillian.
#
# Two dispatch paths via config type:
#   - Config{Lindbladian}  -> Lindbladian eigsolve with :LR targeting
#   - Config{Thermalize}   -> CPTP channel eigsolve with :LM targeting
#
# Covers: MATVEC-06, KRYLOV-01 through KRYLOV-05


# ---------------------------------------------------------------------------
# Pre-flight memory guard
# ---------------------------------------------------------------------------

"""
    _check_krylov_memory(n_qubits, krylovdim)

Advisory pre-flight memory check for Krylov eigsolve.

Estimates memory usage as `krylovdim * 4^n_qubits * 16 * 1.5` bytes (KrylovKit
stores `krylovdim` vectors of length `dim^2 = 4^n`, each ComplexF64 = 16 bytes,
with a 1.5x overhead factor for internal Hessenberg matrix and temporaries).

Issues `@warn` if the estimate exceeds 80% of `Sys.free_memory()`. Does not error.
"""
function _check_krylov_memory(n_qubits::Int, krylovdim::Int)
    estimated_bytes = krylovdim * (4^n_qubits) * 16 * 1.5
    available = Sys.free_memory()
    if estimated_bytes > 0.8 * available
        est_gb = round(estimated_bytes / 1e9; digits=2)
        avail_gb = round(available / 1e9; digits=2)
        @warn "Krylov memory estimate $(est_gb) GB exceeds 80% of free memory $(avail_gb) GB. " *
              "Consider reducing krylovdim or num_qubits."
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Eigsolve with retry logic
# ---------------------------------------------------------------------------

"""
    _eigsolve_with_retry(f, x0, howmany, which; krylovdim=30, tol=1e-10, maxiter=100, max_retries=3, krylov_kwargs...)

Wrapper around `KrylovKit.eigsolve` with automatic retry on partial convergence.

On each retry, increases `krylovdim` by 50% (multiplied by 1.5, rounded up).
Issues `@warn` on each retry attempt.  Throws an error if all attempts fail.

Retry sequence (default): krylovdim = 30 -> 45 -> 68 -> 102.

# Arguments
- `f`: Linear map (callable `Vector -> Vector`)
- `x0`: Initial vector
- `howmany`: Number of eigenvalues requested
- `which`: Eigenvalue targeting (`:LR`, `:LM`, etc.)
- `krylovdim`: Initial Krylov subspace dimension (default 30)
- `tol`: Convergence tolerance (default 1e-10)
- `maxiter`: Maximum restart iterations (default 100)
- `max_retries`: Maximum number of retries (default 3)

# Returns
`(vals, vecs, info)` from KrylovKit.eigsolve
"""
function _eigsolve_with_retry(f, x0, howmany::Int, which::Symbol;
    krylovdim::Int=30, tol::Real=1e-10, maxiter::Int=100, max_retries::Int=3,
    krylov_kwargs...)

    current_krylovdim = krylovdim
    local vals, vecs, info

    for attempt in 1:(max_retries + 1)
        vals, vecs, info = eigsolve(f, x0, howmany, which,
            Arnoldi(; krylovdim=current_krylovdim, tol=tol, maxiter=maxiter, verbosity=0);
            krylov_kwargs...)

        if info.converged >= howmany
            return vals, vecs, info
        end

        if attempt <= max_retries
            new_krylovdim = ceil(Int, current_krylovdim * 1.5)
            @warn "KrylovKit: $(info.converged)/$(howmany) converged. " *
                  "Retrying with krylovdim=$new_krylovdim (attempt $(attempt+1)/$(max_retries+1))"
            current_krylovdim = new_krylovdim
        end
    end

    error("KrylovKit failed to converge: $(info.converged)/$(howmany) eigenvalues " *
          "after $(max_retries + 1) attempts (final krylovdim=$(current_krylovdim))")
end

# ---------------------------------------------------------------------------
# Note: _thermalize_to_liouv_config has been deleted. With unified Config{S,D,C,T},
# the Thermalize path passes Config directly (no conversion needed).
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# apply_delta_channel! -- Faithful Chen CPTP channel (Eq. 3.2)
# ---------------------------------------------------------------------------

"""
    apply_delta_channel!(ws, rho, config, hamiltonian) -> ws.scratch.rho_out

Apply the faithful Chen CPTP channel (Eq. 3.2) to a density matrix.

Uses precomputed channel matrices (K0, U_residual, U_coherent) from the workspace
(populated by the Config{Thermalize} constructor). The per-matvec computation is:

    1. Coherent rotation: rho_eff = U_coherent * rho * U_coherent' (if coherent enabled)
    2. Jump sandwich: rho_jump = delta * sum rate^2 * L * rho_eff * L'
    3. Assembly: E(rho) = K0 * rho_eff * K0' + rho_jump + U_res * rho_eff * U_res'

# Arguments
- `ws::Workspace{KrylovSpectrum}`: Pre-allocated workspace with channel fields populated
- `rho::Matrix{T}`: Input density matrix (dim x dim)
- `config::Config`: Configuration (for sandwich dispatch)
- `hamiltonian::HamHam`: Hamiltonian

# Returns
`ws.scratch.rho_out` containing E(rho).
"""
function apply_delta_channel!(
    ws::Workspace{KrylovSpectrum},
    rho::Matrix{T},
    config::Config,
    hamiltonian::HamHam,
) where {T<:Complex}
    sc = ws.scratch::KrylovScratch{T}
    K0 = ws.K0
    U_res = ws.U_residual
    U_coh = ws.U_coherent
    delta = ws.delta

    # 1. Coherent rotation: rho_eff = U_coh * rho * U_coh'
    #    Use sc.sandwich_out as scratch for rho_eff (safe: sandwich_out not used until sandwich loop)
    if U_coh !== nothing
        mul!(sc.sandwich_tmp, U_coh, rho)
        mul!(sc.sandwich_out, sc.sandwich_tmp, U_coh')
        rho_eff = sc.sandwich_out
    else
        rho_eff = rho
    end

    # 2. Accumulate jump sandwich: rho_jump = delta * sum rate^2 * L * rho_eff * L'
    fill!(sc.channel_rho_jump, 0)
    _accumulate_jump_sandwich!(sc.channel_rho_jump, ws, rho_eff, delta, config, hamiltonian)

    # 3. Need a safe copy of rho_eff before overwriting rho_out
    #    If U_coh !== nothing, rho_eff = sc.sandwich_out (not aliased with rho_out) -- safe
    #    If U_coh === nothing, rho_eff = rho (input arg) -- safe

    # Assembly: rho_out = K0 * rho_eff * K0' + rho_jump + U_res * rho_eff * U_res'
    mul!(sc.sandwich_tmp, K0, rho_eff)
    mul!(sc.rho_out, sc.sandwich_tmp, K0')

    sc.rho_out .+= sc.channel_rho_jump

    mul!(sc.sandwich_tmp, U_res, rho_eff)
    mul!(sc.rho_out, sc.sandwich_tmp, U_res', 1.0, 1.0)

    return sc.rho_out
end

# ---------------------------------------------------------------------------
# _accumulate_jump_sandwich!: domain-dispatched physics-convention sandwiches
# ---------------------------------------------------------------------------

"""
    _accumulate_jump_sandwich!(out, ws, rho, delta, config, hamiltonian) -> nothing

Accumulate `delta * sum rate^2 * L * rho * L'` (physics convention sandwich)
into `out`. Domain-dispatched for EnergyDomain.

Matches the rho_jump accumulation in `_jump_contribution!` for EnergyDomain
(jump_workers.jl) but operating on pre-rotated rho_eff.
"""
function _accumulate_jump_sandwich!(
    out::Matrix{T},
    ws::Workspace{KrylovSpectrum},
    rho::Matrix{T},
    delta::Real,
    config::Config{<:Any, EnergyDomain},
    hamiltonian::HamHam,
) where {T<:Complex}
    sc = ws.scratch::KrylovScratch{T}
    (; transition, gamma_norm_factor, energy_labels) = ws
    bohr_freqs = hamiltonian.bohr_freqs
    inv_4sigma2 = 1.0 / (4 * config.sigma^2)
    prefactor = ws.oft_domain_prefactor * gamma_norm_factor

    if Threads.nthreads() > 1 && length(energy_labels) >= OMEGA_THREAD_THRESHOLD
        return _accumulate_jump_sandwich_threaded_energy!(
            out, sc, rho, delta, ws.jump_eigenbases, ws.jump_hermitian,
            bohr_freqs, energy_labels, transition, prefactor, inv_4sigma2)
    end

    for (k, eigenbasis) in enumerate(ws.jump_eigenbases)
        is_herm = ws.jump_hermitian[k]
        if is_herm
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)
                oft!(sc.jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)
                rate2 = prefactor * transition(w)
                # Physics convention sandwich: delta * rate2 * L * rho * L'
                mul!(sc.sandwich_tmp, rho, sc.jump_oft')          # tmp = rho * L'
                mul!(out, sc.jump_oft, sc.sandwich_tmp, delta * rate2, 1.0)  # out += d*r2 * L * rho * L'
                if w > 1e-12
                    rate2_neg = prefactor * transition(-w)
                    # Neg freq: L_neg = L', sandwich = L' * rho * L
                    mul!(sc.sandwich_tmp, rho, sc.jump_oft)        # tmp = rho * L
                    mul!(out, sc.jump_oft', sc.sandwich_tmp, delta * rate2_neg, 1.0)
                end
            end
        else
            for w in energy_labels
                oft!(sc.jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)
                rate2 = prefactor * transition(w)
                mul!(sc.sandwich_tmp, rho, sc.jump_oft')
                mul!(out, sc.jump_oft, sc.sandwich_tmp, delta * rate2, 1.0)
            end
        end
    end
    return nothing
end

"""
    _accumulate_jump_sandwich!(out, ws, rho, delta, config, hamiltonian) -> nothing

TimeDomain / TrotterDomain version: same structure but uses NUFFT prefactor OFT.
"""
function _accumulate_jump_sandwich!(
    out::Matrix{T},
    ws::Workspace{KrylovSpectrum},
    rho::Matrix{T},
    delta::Real,
    config::Config{<:Any, D},
    hamiltonian::HamHam,
) where {T<:Complex, D<:Union{TimeDomain, TrotterDomain}}
    sc = ws.scratch::KrylovScratch{T}
    (; transition, gamma_norm_factor, energy_labels, oft_nufft_prefactors) = ws
    prefactor = ws.oft_domain_prefactor * gamma_norm_factor

    if Threads.nthreads() > 1 && length(energy_labels) >= OMEGA_THREAD_THRESHOLD
        return _accumulate_jump_sandwich_threaded_timetrot!(
            out, sc, rho, delta, ws.jump_eigenbases, ws.jump_hermitian,
            oft_nufft_prefactors.data, oft_nufft_prefactors.energy_to_index,
            energy_labels, transition, prefactor)
    end

    for (k, eigenbasis) in enumerate(ws.jump_eigenbases)
        is_herm = ws.jump_hermitian[k]
        if is_herm
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)
                nufft_pf = _prefactor_view(oft_nufft_prefactors, w)
                @. sc.jump_oft = eigenbasis * nufft_pf
                rate2 = prefactor * transition(w)
                mul!(sc.sandwich_tmp, rho, sc.jump_oft')
                mul!(out, sc.jump_oft, sc.sandwich_tmp, delta * rate2, 1.0)
                if w > 1e-12
                    rate2_neg = prefactor * transition(-w)
                    mul!(sc.sandwich_tmp, rho, sc.jump_oft)
                    mul!(out, sc.jump_oft', sc.sandwich_tmp, delta * rate2_neg, 1.0)
                end
            end
        else
            for w in energy_labels
                nufft_pf = _prefactor_view(oft_nufft_prefactors, w)
                @. sc.jump_oft = eigenbasis * nufft_pf
                rate2 = prefactor * transition(w)
                mul!(sc.sandwich_tmp, rho, sc.jump_oft')
                mul!(out, sc.jump_oft, sc.sandwich_tmp, delta * rate2, 1.0)
            end
        end
    end
    return nothing
end

"""
    _accumulate_jump_sandwich!(out, ws, rho, delta, config, hamiltonian) -> nothing

BohrDomain version: iterates over Bohr frequency buckets.

The physics-convention sandwich for BohrDomain is:
    rho_jump += delta * gamma_norm_factor * alpha_A * rho * A_nu2_dag
Matching jump_workers.jl:276-277.
"""
function _accumulate_jump_sandwich!(
    out::Matrix{T},
    ws::Workspace{KrylovSpectrum},
    rho::Matrix{T},
    delta::Real,
    config::Config{<:Any, BohrDomain},
    hamiltonian::HamHam,
) where {T<:Complex}
    sc = ws.scratch::KrylovScratch{T}
    bohr_alpha_fn = ws.bohr_alpha
    gamma_norm_factor = ws.gamma_norm_factor
    dim = size(rho, 1)

    bohr_keys = ws.bohr_keys
    bohr_keys_local = bohr_keys === nothing ? collect(keys(hamiltonian.bohr_dict)) : bohr_keys
    n_jumps = length(ws.jump_eigenbases)
    n_keys  = length(bohr_keys_local)
    if Threads.nthreads() > 1 && n_jumps * n_keys >= OMEGA_THREAD_THRESHOLD
        return _accumulate_jump_sandwich_threaded_bohr!(
            out, sc, rho, delta, ws.jump_eigenbases,
            bohr_alpha_fn, gamma_norm_factor,
            hamiltonian.bohr_freqs, hamiltonian.bohr_dict, bohr_keys_local)
    end

    for (k, eigenbasis) in enumerate(ws.jump_eigenbases)
        for nu_2 in bohr_keys_local
            # alpha_A = B_nu2
            @. sc.jump_oft = bohr_alpha_fn(hamiltonian.bohr_freqs, nu_2) * eigenbasis

            # Build A_nu2_dag: entrywise rho*A_nu2_dag via scatter (matches thermalization code)
            fill!(sc.sandwich_tmp, 0)
            indices = hamiltonian.bohr_dict[nu_2]
            @inbounds for idx in indices
                i = idx[1]; j = idx[2]
                v = conj(eigenbasis[i, j])
                @inbounds for p in 1:dim
                    sc.sandwich_tmp[p, i] += rho[p, j] * v  # sandwich_tmp = rho * A_nu2_dag
                end
            end

            # out += delta * gamma_norm_factor * alpha_A * (rho * A_nu2_dag)
            mul!(out, sc.jump_oft, sc.sandwich_tmp, delta * gamma_norm_factor, 1.0)
        end
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Threaded ω-loop variants for the channel jump sandwich (qf-in3 follow-up)
# ---------------------------------------------------------------------------
#
# Mirror of `_apply_lindbladian_threaded_*!` in `src/krylov_matvec.jl`, but for
# the additive accumulator pattern of `_accumulate_jump_sandwich!`. Uses
# `sc.task_scratches[idx].channel_rho_jump` as the per-thread accumulator
# (Thermalize Workspaces are constructed with `with_channel_rho_jump=true`),
# and `sc.work_list` as the pre-built `(jump_idx, label_idx)` schedule.
# Reduces the per-thread chunks into the caller's `out` buffer at the end.

function _accumulate_jump_sandwich_threaded_energy!(
    out::Matrix{T},
    sc::KrylovScratch{T},
    rho::Matrix{T},
    delta::Real,
    jump_eigenbases::Vector{Matrix{T}},
    jump_hermitian::Vector{Bool},
    bohr_freqs::AbstractMatrix{<:Real},
    energy_labels::Vector{Float64},
    transition,
    prefactor::Float64,
    inv_4sigma2::Float64,
) where {T<:Complex}
    work = sc.work_list
    _populate_lindblad_work_list!(work, jump_hermitian, energy_labels)
    n_work = length(work)
    n_work == 0 && return nothing

    pool = sc.task_scratches
    nt = min(Threads.nthreads(), n_work, length(pool))
    chunks = _partition_range(1:n_work, nt)

    old_blas = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        @sync for (idx, chunk) in enumerate(chunks)
            Threads.@spawn _accumulate_jump_sandwich_chunk_energy!(
                pool[idx], rho, delta, jump_eigenbases, jump_hermitian,
                bohr_freqs, energy_labels, transition, work, chunk,
                prefactor, inv_4sigma2)
        end
    finally
        BLAS.set_num_threads(old_blas)
    end

    @inbounds for idx in 1:length(chunks)
        out .+= pool[idx].channel_rho_jump
    end
    return nothing
end

function _accumulate_jump_sandwich_chunk_energy!(
    task_sc::KrylovScratch{T},
    rho::Matrix{T},
    delta::Real,
    jump_eigenbases::Vector{Matrix{T}},
    jump_hermitian::Vector{Bool},
    bohr_freqs::AbstractMatrix{<:Real},
    energy_labels::Vector{Float64},
    transition,
    work::Vector{Tuple{Int, Int}},
    chunk::UnitRange{Int},
    prefactor::Float64,
    inv_4sigma2::Float64,
) where {T<:Complex}
    fill!(task_sc.channel_rho_jump, 0)

    @inbounds for w_idx in chunk
        (k, li) = work[w_idx]
        eigenbasis = jump_eigenbases[k]
        is_herm = jump_hermitian[k]

        w_raw = energy_labels[li]
        # Hermitian fold: only `w_raw <= 1e-12` queued; OFT and rate use
        # `|w_raw|`. Non-Hermitian uses `w_raw` directly (signed).
        w = is_herm ? abs(w_raw) : w_raw

        oft!(task_sc.jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)
        rate2 = prefactor * transition(w)
        mul!(task_sc.sandwich_tmp, rho, task_sc.jump_oft')
        mul!(task_sc.channel_rho_jump, task_sc.jump_oft, task_sc.sandwich_tmp,
             delta * rate2, 1.0)

        if is_herm && w > 1e-12
            rate2_neg = prefactor * transition(-w)
            mul!(task_sc.sandwich_tmp, rho, task_sc.jump_oft)
            mul!(task_sc.channel_rho_jump, task_sc.jump_oft', task_sc.sandwich_tmp,
                 delta * rate2_neg, 1.0)
        end
    end
    return nothing
end

function _accumulate_jump_sandwich_threaded_timetrot!(
    out::Matrix{T},
    sc::KrylovScratch{T},
    rho::Matrix{T},
    delta::Real,
    jump_eigenbases::Vector{Matrix{T}},
    jump_hermitian::Vector{Bool},
    nufft_data::AbstractArray{T, 3},
    nufft_idx::AbstractDict,
    energy_labels::Vector{Float64},
    transition,
    prefactor::Float64,
) where {T<:Complex}
    work = sc.work_list
    _populate_lindblad_work_list!(work, jump_hermitian, energy_labels)
    n_work = length(work)
    n_work == 0 && return nothing

    pool = sc.task_scratches
    nt = min(Threads.nthreads(), n_work, length(pool))
    chunks = _partition_range(1:n_work, nt)

    old_blas = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        @sync for (idx, chunk) in enumerate(chunks)
            Threads.@spawn _accumulate_jump_sandwich_chunk_timetrot!(
                pool[idx], rho, delta, jump_eigenbases, jump_hermitian,
                nufft_data, nufft_idx, energy_labels, transition, work, chunk,
                prefactor)
        end
    finally
        BLAS.set_num_threads(old_blas)
    end

    @inbounds for idx in 1:length(chunks)
        out .+= pool[idx].channel_rho_jump
    end
    return nothing
end

function _accumulate_jump_sandwich_chunk_timetrot!(
    task_sc::KrylovScratch{T},
    rho::Matrix{T},
    delta::Real,
    jump_eigenbases::Vector{Matrix{T}},
    jump_hermitian::Vector{Bool},
    nufft_data::AbstractArray{T, 3},
    nufft_idx::AbstractDict,
    energy_labels::Vector{Float64},
    transition,
    work::Vector{Tuple{Int, Int}},
    chunk::UnitRange{Int},
    prefactor::Float64,
) where {T<:Complex}
    fill!(task_sc.channel_rho_jump, 0)

    @inbounds for w_idx in chunk
        (k, li) = work[w_idx]
        eigenbasis = jump_eigenbases[k]
        is_herm = jump_hermitian[k]

        w_raw = energy_labels[li]
        w = is_herm ? abs(w_raw) : w_raw
        prefactor_idx = is_herm ? nufft_idx[w] : li
        nufft_pf = @view nufft_data[:, :, prefactor_idx]
        @. task_sc.jump_oft = eigenbasis * nufft_pf
        rate2 = prefactor * transition(w)

        mul!(task_sc.sandwich_tmp, rho, task_sc.jump_oft')
        mul!(task_sc.channel_rho_jump, task_sc.jump_oft, task_sc.sandwich_tmp,
             delta * rate2, 1.0)

        if is_herm && w > 1e-12
            rate2_neg = prefactor * transition(-w)
            mul!(task_sc.sandwich_tmp, rho, task_sc.jump_oft)
            mul!(task_sc.channel_rho_jump, task_sc.jump_oft', task_sc.sandwich_tmp,
                 delta * rate2_neg, 1.0)
        end
    end
    return nothing
end

# --- BohrDomain channel sandwich threaded variant (qf-6af.7) ---
#
# Mirror of `_accumulate_jump_sandwich_threaded_timetrot!`, but with work
# decomposed over (jump_idx, bohr_key_idx) pairs. BohrDomain dissipator is
# `ρ_jump += δ·γ_nf · α(ν₂) · ρ · A_{ν₂}^†`, so the per-key inner work is
# (1) form `α·A_{eb}` (entrywise scaling), (2) scatter `rho · A_{ν₂}^†` via
# the precomputed `bohr_dict[ν₂]` index list, (3) `mul!(out, …, …, w, 1)`.

function _accumulate_jump_sandwich_threaded_bohr!(
    out::Matrix{T},
    sc::KrylovScratch{T},
    rho::Matrix{T},
    delta::Real,
    jump_eigenbases::Vector{Matrix{T}},
    bohr_alpha_fn,
    gamma_norm_factor::Real,
    bohr_freqs::AbstractMatrix{<:Real},
    bohr_dict,
    bohr_keys::Vector,
) where {T<:Complex}
    n_jumps = length(jump_eigenbases)
    n_keys  = length(bohr_keys)
    n_work  = n_jumps * n_keys
    n_work == 0 && return nothing

    pool = sc.task_scratches
    nt = min(Threads.nthreads(), n_work, length(pool))
    chunks = _partition_range(1:n_work, nt)
    n_chunks = length(chunks)

    old_blas = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        @sync for (idx, chunk) in enumerate(chunks)
            Threads.@spawn _accumulate_jump_sandwich_chunk_bohr!(
                pool[idx], rho, delta, jump_eigenbases,
                bohr_alpha_fn, gamma_norm_factor,
                bohr_freqs, bohr_dict, bohr_keys, n_keys, chunk)
        end
    finally
        BLAS.set_num_threads(old_blas)
    end

    @inbounds for idx in 1:n_chunks
        out .+= pool[idx].channel_rho_jump
    end
    return nothing
end

function _accumulate_jump_sandwich_chunk_bohr!(
    task_sc::KrylovScratch{T},
    rho::Matrix{T},
    delta::Real,
    jump_eigenbases::Vector{Matrix{T}},
    bohr_alpha_fn,
    gamma_norm_factor::Real,
    bohr_freqs::AbstractMatrix{<:Real},
    bohr_dict,
    bohr_keys::Vector,
    n_keys::Int,
    chunk::UnitRange{Int},
) where {T<:Complex}
    fill!(task_sc.channel_rho_jump, 0)
    dim = size(rho, 1)

    @inbounds for w_idx in chunk
        k       = ((w_idx - 1) ÷ n_keys) + 1
        key_idx = ((w_idx - 1) % n_keys) + 1
        eigenbasis = jump_eigenbases[k]
        nu_2 = bohr_keys[key_idx]

        @. task_sc.jump_oft = bohr_alpha_fn(bohr_freqs, nu_2) * eigenbasis

        fill!(task_sc.sandwich_tmp, 0)
        indices = bohr_dict[nu_2]
        @inbounds for idx in indices
            i = idx[1]; j = idx[2]
            v = conj(eigenbasis[i, j])
            @inbounds for p in 1:dim
                task_sc.sandwich_tmp[p, i] += rho[p, j] * v
            end
        end

        mul!(task_sc.channel_rho_jump, task_sc.jump_oft, task_sc.sandwich_tmp,
             delta * gamma_norm_factor, 1.0)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# krylov_spectral_gap: Lindbladian path (Config{Lindbladian})
# ---------------------------------------------------------------------------

"""
    krylov_spectral_gap(config::Config{Lindbladian}, hamiltonian, jumps; kwargs...) -> NamedTuple

Compute the Lindbladian spectral gap matrix-free using KrylovKit Arnoldi with `:LR` targeting.

Wraps the existing `apply_lindbladian!` matvec in a KrylovKit closure to find the leading
eigenvalues of the Lindbladian superoperator without constructing the full dense matrix.

The steady-state eigenvalue is near Re(lambda) ~ 0 (largest real part). The spectral gap
is `abs(real(lambda_2))` where lambda_2 is the second eigenvalue sorted by |Re(lambda)|.

# Arguments
- `config::Config{Lindbladian}`: Lindbladian configuration (EnergyDomain, TimeDomain, etc.)
- `hamiltonian::HamHam`: Hamiltonian with eigenbasis data
- `jumps::Vector{JumpOp}`: Jump operators

# Keyword Arguments
- `trotter::Union{TrottTrott, Nothing}=nothing`: Trotter object (required for TrotterDomain)
- `krylovdim::Int=30`: Krylov subspace dimension
- `howmany::Int=4`: Number of eigenvalues to compute
- `tol::Real=1e-10`: Convergence tolerance
- `max_retries::Int=3`: Maximum retry attempts on partial convergence
- `krylov_kwargs...`: Additional keyword arguments passed to KrylovKit

# Returns
A NamedTuple with Lindbladian eigenvalues, spectral gap, fixed point, and gap mode.
"""
function krylov_spectral_gap(
    config::Config{Lindbladian},
    hamiltonian::HamHam,
    jumps::Vector{JumpOp};
    trotter::Union{TrottTrott, Nothing}=nothing,
    krylovdim::Int=30,
    howmany::Int=4,
    tol::Real=1e-10,
    max_retries::Int=3,
    allow_unpaired_nonhermitian::Bool=false,
    krylov_kwargs...
)
    # Guards
    krylovdim > howmany || error("krylovdim ($krylovdim) must be > howmany ($howmany)")
    _check_krylov_memory(config.num_qubits, krylovdim)
    validate_jump_pairing(jumps; allow_unpaired_nonhermitian=allow_unpaired_nonhermitian)

    # Allocate workspace
    ws = Workspace(config, hamiltonian, jumps; trotter=trotter)

    # Dimensions
    dim = size(hamiltonian.data, 1)

    # Build matvec closure: Vector{ComplexF64} -> Vector{ComplexF64}
    function lindbladian_matvec(v::AbstractVector)
        rho = reshape(v, dim, dim)
        apply_lindbladian!(ws, rho, config, hamiltonian)
        return copy(vec(ws.scratch.rho_out))  # CRITICAL: copy to avoid aliasing (Pitfall 1)
    end

    # Initial vector: maximally mixed state
    x0 = vec(Matrix{ComplexF64}(I(dim) / dim))

    # Eigsolve with retry
    vals, vecs, info = _eigsolve_with_retry(
        lindbladian_matvec, x0, howmany, :LR;
        krylovdim=krylovdim, tol=tol, max_retries=max_retries,
        krylov_kwargs...)

    # Sort eigenvalues by |Re(lambda)| ascending (steady state first, then gap mode)
    perm = sortperm(vals; by=v -> abs(real(v)))

    eigenvalues_sorted = vals[perm]
    vecs_sorted = vecs[perm]

    # Extract fixed_point (eigenvector 1): reshape, hermitianize, trace-normalize
    fixed_point = reshape(vecs_sorted[1], dim, dim)
    hermitianize!(fixed_point)
    fixed_point ./= tr(fixed_point)

    # Extract gap_mode (eigenvector 2): reshape only
    gap_mode = reshape(vecs_sorted[2], dim, dim)

    # Spectral gap = abs(real(lambda_2))
    spectral_gap = abs(real(eigenvalues_sorted[2]))

    # Residual norms (reorder to match sorted eigenvalues)
    normres = Float64.(info.normres[perm])

    return (;
        eigenvalues = Complex{Float64}.(eigenvalues_sorted),
        spectral_gap,
        fixed_point = Complex{Float64}.(fixed_point),
        gap_mode = Complex{Float64}.(gap_mode),
        converged = info.converged,
        matvec_count = info.numops,
        num_restarts = info.numiter,
        normres,
        channel_eigenvalues = nothing,
        delta_used = nothing,
    )
end

# ---------------------------------------------------------------------------
# krylov_spectral_gap: Channel path (Config{Thermalize})
# ---------------------------------------------------------------------------

"""
    krylov_spectral_gap(config::Config{Thermalize}, hamiltonian, jumps; kwargs...) -> NamedTuple

Compute the Lindbladian spectral gap via the faithful Chen CPTP channel (Eq. 3.2),
using KrylovKit Arnoldi with `:LM` targeting.

The channel eigenvalues mu are related to Lindbladian eigenvalues by the first-order
approximation: `lambda_L = (mu - 1) / delta`. Since mu = exp(delta * lambda_L) + O(delta^2),
the conversion introduces O(delta) error. The steady state has mu ~ 1 (largest magnitude),
and the gap is recovered from the second eigenvalue after conversion.

The channel is CPTP and O(delta^2) accurate, matching what `run_thermalize`
actually implements via `_finalize_kraus_step!`.

# Arguments
- `config::Config{Thermalize}`: Thermalization configuration (provides delta)
- `hamiltonian::HamHam`: Hamiltonian with eigenbasis data
- `jumps::Vector{JumpOp}`: Jump operators

# Keyword Arguments
Same as the `Config{Lindbladian}` method.

# Returns
A NamedTuple with converted Lindbladian eigenvalues, spectral gap, fixed point,
gap mode, and the raw channel eigenvalues stored in `channel_eigenvalues`.
"""
function krylov_spectral_gap(
    config::Config{Thermalize},
    hamiltonian::HamHam,
    jumps::Vector{JumpOp};
    trotter::Union{TrottTrott, Nothing}=nothing,
    krylovdim::Int=30,
    howmany::Int=4,
    tol::Real=1e-10,
    max_retries::Int=3,
    allow_unpaired_nonhermitian::Bool=false,
    krylov_kwargs...
)
    # Guards
    krylovdim > howmany || error("krylovdim ($krylovdim) must be > howmany ($howmany)")
    _check_krylov_memory(config.num_qubits, krylovdim)
    validate_jump_pairing(jumps; allow_unpaired_nonhermitian=allow_unpaired_nonhermitian)

    # Get delta from config
    delta = config.delta

    # Allocate workspace using Config{Thermalize} constructor (precomputes channel matrices)
    ws = Workspace(config, hamiltonian, jumps; trotter=trotter)

    # Dimensions
    dim = size(hamiltonian.data, 1)

    # Build channel matvec closure using faithful Chen channel
    function channel_matvec(v::AbstractVector)
        rho = reshape(v, dim, dim)
        apply_delta_channel!(ws, rho, config, hamiltonian)
        return copy(vec(ws.scratch.rho_out))  # CRITICAL: copy to avoid aliasing
    end

    # Initial vector: maximally mixed state
    x0 = vec(Matrix{ComplexF64}(I(dim) / dim))

    # Eigsolve with :LM targeting (channel eigenvalues cluster near 1)
    vals, vecs, info = _eigsolve_with_retry(
        channel_matvec, x0, howmany, :LM;
        krylovdim=krylovdim, tol=tol, max_retries=max_retries,
        krylov_kwargs...)

    # Store raw channel eigenvalues before conversion
    channel_eigenvalues_raw = Complex{Float64}.(vals)

    # Convert channel eigenvalues to Lindbladian eigenvalues: lambda_L = (mu - 1) / delta
    lindblad_eigenvalues = (vals .- 1) ./ delta

    # Sort Lindbladian eigenvalues by |Re(lambda)| ascending (steady state first)
    perm = sortperm(lindblad_eigenvalues; by=v -> abs(real(v)))

    eigenvalues_sorted = lindblad_eigenvalues[perm]
    vecs_sorted = vecs[perm]
    channel_eigenvalues_sorted = channel_eigenvalues_raw[perm]

    # Extract fixed_point (eigenvector 1): reshape, hermitianize, trace-normalize
    fixed_point = reshape(vecs_sorted[1], dim, dim)
    hermitianize!(fixed_point)
    fixed_point ./= tr(fixed_point)

    # Extract gap_mode (eigenvector 2): reshape only
    gap_mode = reshape(vecs_sorted[2], dim, dim)

    # Spectral gap = abs(real(lambda_L_2))
    spectral_gap = abs(real(eigenvalues_sorted[2]))

    # Residual norms (reorder to match sorted eigenvalues)
    normres = Float64.(info.normres[perm])

    return (;
        eigenvalues = Complex{Float64}.(eigenvalues_sorted),
        spectral_gap,
        fixed_point = Complex{Float64}.(fixed_point),
        gap_mode = Complex{Float64}.(gap_mode),
        converged = info.converged,
        matvec_count = info.numops,
        num_restarts = info.numiter,
        normres,
        channel_eigenvalues = channel_eigenvalues_sorted,
        delta_used = Float64(delta),
    )
end

# ---------------------------------------------------------------------------
# run_krylov_spectrum: new public entry point (Phase 36)
# ---------------------------------------------------------------------------

"""
    run_krylov_spectrum(jumps, config, hamiltonian, trotter=nothing; krylov_kwargs...) -> KrylovSpectrumResults

Matrix-free spectral gap computation via KrylovKit, returning a `KrylovSpectrumResults` struct.

Wraps `krylov_spectral_gap` with the uniform positional signature
`(jumps, config, hamiltonian, trotter)` and adds config + timing metadata to the result.

Works for both `Config{Lindbladian}` (direct Lindbladian matvec with `:LR` targeting)
and `Config{Thermalize}` (CPTP channel matvec with `:LM` targeting and delta conversion).

# Arguments
- `jumps::Vector{JumpOp}`: Jump operators
- `config::Config`: Lindbladian or Thermalize configuration
- `hamiltonian::HamHam`: Hamiltonian with eigenbasis data
- `trotter::Union{TrottTrott, Nothing}=nothing`: Trotter object (required for TrotterDomain)

# Keyword Arguments
- `krylovdim::Int=30`: Krylov subspace dimension
- `howmany::Int=4`: Number of eigenvalues to compute
- `tol::Real=1e-10`: Convergence tolerance
- `max_retries::Int=3`: Maximum retry attempts on partial convergence
- `krylov_kwargs...`: Additional keyword arguments passed to KrylovKit

# Returns
`KrylovSpectrumResults` with eigenvalues, spectral gap, fixed point, gap mode,
convergence info, optional channel data, config, and timing metadata.
"""
function run_krylov_spectrum(
    jumps::Vector{JumpOp},
    config::Config{S,D,C,T},
    hamiltonian::HamHam,
    trotter::Union{TrottTrott, Nothing}=nothing;
    krylovdim::Int=30,
    howmany::Int=4,
    tol::Real=1e-10,
    max_retries::Int=3,
    allow_unpaired_nonhermitian::Bool=false,
    krylov_kwargs...
) where {S<:Union{Lindbladian, Thermalize}, D, C, T}

    t_start = time()

    # Delegate to existing krylov_spectral_gap (handles both Config{Lindbladian} and Config{Thermalize})
    # Note: krylov_spectral_gap uses (config, hamiltonian, jumps; ...) argument order
    krylov_result = krylov_spectral_gap(
        config, hamiltonian, jumps;
        trotter=trotter, krylovdim=krylovdim, howmany=howmany,
        tol=tol, max_retries=max_retries,
        allow_unpaired_nonhermitian=allow_unpaired_nonhermitian,
        krylov_kwargs...
    )

    wall_time = time() - t_start
    metadata = _capture_metadata(wall_time_seconds=wall_time)

    return KrylovSpectrumResults{Float64}(
        config,
        krylov_result.eigenvalues,
        krylov_result.spectral_gap,
        krylov_result.fixed_point,
        krylov_result.gap_mode,
        krylov_result.converged,
        krylov_result.matvec_count,
        krylov_result.num_restarts,
        krylov_result.normres,
        krylov_result.channel_eigenvalues,
        krylov_result.delta_used,
        metadata,
    )
end
