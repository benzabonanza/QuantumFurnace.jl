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
# Result struct
# ---------------------------------------------------------------------------

"""
    KrylovGapResult{T<:AbstractFloat}

Result of a matrix-free Krylov spectral gap computation.

Stores the leading Lindbladian eigenvalues, the spectral gap, the fixed-point
density matrix (steady state), and the gap mode operator.  For the channel path
(`Config{Thermalize}`), the raw channel eigenvalues and delta are also stored.

# Fields
- `eigenvalues`: Leading Lindbladian eigenvalues sorted by |Re(lambda)| ascending.
- `spectral_gap`: abs(real(eigenvalues[2])) -- the Lindbladian spectral gap.
- `fixed_point`: Steady-state density matrix (eigvec 1, Hermitianized, trace-normalized).
- `gap_mode`: Gap mode operator (eigvec 2, reshaped to dim x dim).
- `converged`: Number of converged eigenvalues reported by KrylovKit.
- `matvec_count`: Number of linear map applications (info.numops).
- `num_restarts`: Number of Krylov subspace rebuilds (info.numiter).
- `normres`: Residual norms for each eigenvalue.
- `channel_eigenvalues`: Raw channel eigenvalues before conversion (nothing for Lindbladian path).
- `delta_used`: Delta from Config{Thermalize} (nothing for Lindbladian path).
"""
struct KrylovGapResult{T<:AbstractFloat}
    eigenvalues::Vector{Complex{T}}
    spectral_gap::T
    fixed_point::Matrix{Complex{T}}
    gap_mode::Matrix{Complex{T}}
    converged::Int
    matvec_count::Int
    num_restarts::Int
    normres::Vector{T}
    channel_eigenvalues::Union{Nothing, Vector{Complex{T}}}
    delta_used::Union{Nothing, T}
end

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
    apply_delta_channel!(ws, rho, config, hamiltonian) -> ws.rho_out

Apply the faithful Chen CPTP channel (Eq. 3.2) to a density matrix.

Uses precomputed channel matrices (K0, U_residual, U_coherent) from the workspace
(populated by the Config{Thermalize} constructor). The per-matvec computation is:

    1. Coherent rotation: rho_eff = U_coherent * rho * U_coherent' (if coherent enabled)
    2. Jump sandwich: rho_jump = delta * sum rate^2 * L * rho_eff * L'
    3. Assembly: E(rho) = K0 * rho_eff * K0' + rho_jump + U_res * rho_eff * U_res'

# Arguments
- `ws::KrylovWorkspace{T}`: Pre-allocated workspace with channel fields populated
- `rho::Matrix{T}`: Input density matrix (dim x dim)
- `config::Config`: Configuration (for sandwich dispatch)
- `hamiltonian::HamHam`: Hamiltonian

# Returns
`ws.rho_out` containing E(rho).
"""
function apply_delta_channel!(
    ws::KrylovWorkspace{T},
    rho::Matrix{T},
    config::Config,
    hamiltonian::HamHam,
) where {T<:Complex}
    K0 = ws.channel_K0
    U_res = ws.channel_U_residual
    U_coh = ws.channel_U_coherent
    delta = ws.channel_delta

    # 1. Coherent rotation: rho_eff = U_coh * rho * U_coh'
    #    Use ws.LdagL as scratch for rho_eff (safe: LdagL not used until sandwich loop)
    if U_coh !== nothing
        mul!(ws.tmp1, U_coh, rho)
        mul!(ws.LdagL, ws.tmp1, U_coh')
        rho_eff = ws.LdagL
    else
        rho_eff = rho
    end

    # 2. Accumulate jump sandwich: rho_jump = delta * sum rate^2 * L * rho_eff * L'
    fill!(ws.channel_rho_jump, 0)
    _accumulate_jump_sandwich!(ws.channel_rho_jump, ws, rho_eff, delta, config, hamiltonian)

    # 3. Need a safe copy of rho_eff before overwriting rho_out
    #    If U_coh !== nothing, rho_eff = ws.LdagL (not aliased with rho_out) -- safe
    #    If U_coh === nothing, rho_eff = rho (input arg) -- safe
    #    So rho_eff is always safe to read after we write rho_out.

    # Assembly: rho_out = K0 * rho_eff * K0' + rho_jump + U_res * rho_eff * U_res'
    mul!(ws.tmp1, K0, rho_eff)
    mul!(ws.rho_out, ws.tmp1, K0')

    ws.rho_out .+= ws.channel_rho_jump

    mul!(ws.tmp1, U_res, rho_eff)
    mul!(ws.rho_out, ws.tmp1, U_res', 1.0, 1.0)

    return ws.rho_out
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
    ws::KrylovWorkspace{T},
    rho::Matrix{T},
    delta::Real,
    config::Config{<:Any, EnergyDomain},
    hamiltonian::HamHam,
) where {T<:Complex}
    (; transition, gamma_norm_factor, energy_labels) = ws.precomputed_data
    bohr_freqs = hamiltonian.bohr_freqs
    inv_4sigma2 = 1.0 / (4 * config.sigma^2)
    prefactor = ws.precomputed_data.domain_prefactor * gamma_norm_factor

    for (k, eigenbasis) in enumerate(ws.jump_eigenbases)
        is_herm = ws.jump_hermitian[k]
        if is_herm
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)
                oft!(ws.jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)
                rate2 = prefactor * transition(w)
                # Physics convention sandwich: delta * rate2 * L * rho * L'
                mul!(ws.tmp1, rho, ws.jump_oft')          # tmp1 = rho * L'
                mul!(out, ws.jump_oft, ws.tmp1, delta * rate2, 1.0)  # out += d*r2 * L * rho * L'
                if w > 1e-12
                    rate2_neg = prefactor * transition(-w)
                    # Neg freq: L_neg = L', sandwich = L' * rho * L
                    mul!(ws.tmp1, rho, ws.jump_oft)        # tmp1 = rho * L
                    mul!(out, ws.jump_oft', ws.tmp1, delta * rate2_neg, 1.0)
                end
            end
        else
            for w in energy_labels
                oft!(ws.jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)
                rate2 = prefactor * transition(w)
                mul!(ws.tmp1, rho, ws.jump_oft')
                mul!(out, ws.jump_oft, ws.tmp1, delta * rate2, 1.0)
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
    ws::KrylovWorkspace{T},
    rho::Matrix{T},
    delta::Real,
    config::Config{<:Any, D},
    hamiltonian::HamHam,
) where {T<:Complex, D<:Union{TimeDomain, TrotterDomain}}
    (; transition, gamma_norm_factor, energy_labels, oft_nufft_prefactors) = ws.precomputed_data
    prefactor = ws.precomputed_data.domain_prefactor * gamma_norm_factor

    for (k, eigenbasis) in enumerate(ws.jump_eigenbases)
        is_herm = ws.jump_hermitian[k]
        if is_herm
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)
                nufft_pf = _prefactor_view(oft_nufft_prefactors, w)
                @. ws.jump_oft = eigenbasis * nufft_pf
                rate2 = prefactor * transition(w)
                mul!(ws.tmp1, rho, ws.jump_oft')
                mul!(out, ws.jump_oft, ws.tmp1, delta * rate2, 1.0)
                if w > 1e-12
                    rate2_neg = prefactor * transition(-w)
                    mul!(ws.tmp1, rho, ws.jump_oft)
                    mul!(out, ws.jump_oft', ws.tmp1, delta * rate2_neg, 1.0)
                end
            end
        else
            for w in energy_labels
                nufft_pf = _prefactor_view(oft_nufft_prefactors, w)
                @. ws.jump_oft = eigenbasis * nufft_pf
                rate2 = prefactor * transition(w)
                mul!(ws.tmp1, rho, ws.jump_oft')
                mul!(out, ws.jump_oft, ws.tmp1, delta * rate2, 1.0)
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
    ws::KrylovWorkspace{T},
    rho::Matrix{T},
    delta::Real,
    config::Config{<:Any, BohrDomain},
    hamiltonian::HamHam,
) where {T<:Complex}
    (; alpha, gamma_norm_factor) = ws.precomputed_data
    dim = size(rho, 1)
    # Allocate A_nu2_dag buffer (one per call, acceptable for Bohr)
    A_nu2_dag = zeros(T, dim, dim)

    for (k, eigenbasis) in enumerate(ws.jump_eigenbases)
        for nu_2 in keys(hamiltonian.bohr_dict)
            # alpha_A = B_nu2
            @. ws.jump_oft = alpha(hamiltonian.bohr_freqs, nu_2) * eigenbasis

            # Build A_nu2_dag: entrywise rho*A_nu2_dag via scatter (matches thermalization code)
            fill!(ws.tmp1, 0)
            indices = hamiltonian.bohr_dict[nu_2]
            @inbounds for idx in indices
                i = idx[1]; j = idx[2]
                v = conj(eigenbasis[i, j])
                @inbounds for p in 1:dim
                    ws.tmp1[p, i] += rho[p, j] * v  # tmp1 = rho * A_nu2_dag
                end
            end

            # out += delta * gamma_norm_factor * alpha_A * (rho * A_nu2_dag)
            mul!(out, ws.jump_oft, ws.tmp1, delta * gamma_norm_factor, 1.0)
        end
    end
    return nothing
end

# ---------------------------------------------------------------------------
# krylov_spectral_gap: Lindbladian path (Config{Lindbladian})
# ---------------------------------------------------------------------------

"""
    krylov_spectral_gap(config::Config{Lindbladian}, hamiltonian, jumps; kwargs...) -> KrylovGapResult

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
`KrylovGapResult` with Lindbladian eigenvalues, spectral gap, fixed point, and gap mode.
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
    krylov_kwargs...
)
    # Guards
    krylovdim > howmany || error("krylovdim ($krylovdim) must be > howmany ($howmany)")
    _check_krylov_memory(config.num_qubits, krylovdim)

    # Allocate workspace
    ws = KrylovWorkspace(config, hamiltonian, jumps; trotter=trotter)

    # Dimensions
    dim = size(hamiltonian.data, 1)

    # Build matvec closure: Vector{ComplexF64} -> Vector{ComplexF64}
    function lindbladian_matvec(v::AbstractVector)
        rho = reshape(v, dim, dim)
        apply_lindbladian!(ws, rho, config, hamiltonian)
        return copy(vec(ws.rho_out))  # CRITICAL: copy to avoid aliasing (Pitfall 1)
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

    return KrylovGapResult{Float64}(
        Complex{Float64}.(eigenvalues_sorted),
        spectral_gap,
        Complex{Float64}.(fixed_point),
        Complex{Float64}.(gap_mode),
        info.converged,
        info.numops,
        info.numiter,
        normres,
        nothing,   # channel_eigenvalues
        nothing,   # delta_used
    )
end

# ---------------------------------------------------------------------------
# krylov_spectral_gap: Channel path (Config{Thermalize})
# ---------------------------------------------------------------------------

"""
    krylov_spectral_gap(config::Config{Thermalize}, hamiltonian, jumps; kwargs...) -> KrylovGapResult

Compute the Lindbladian spectral gap via the faithful Chen CPTP channel (Eq. 3.2),
using KrylovKit Arnoldi with `:LM` targeting.

The channel eigenvalues mu are related to Lindbladian eigenvalues by the first-order
approximation: `lambda_L = (mu - 1) / delta`. Since mu = exp(delta * lambda_L) + O(delta^2),
the conversion introduces O(delta) error. The steady state has mu ~ 1 (largest magnitude),
and the gap is recovered from the second eigenvalue after conversion.

The channel is CPTP and O(delta^2) accurate, matching what `run_thermalization`
actually implements via `_finalize_kraus_step!`.

# Arguments
- `config::Config{Thermalize}`: Thermalization configuration (provides delta)
- `hamiltonian::HamHam`: Hamiltonian with eigenbasis data
- `jumps::Vector{JumpOp}`: Jump operators

# Keyword Arguments
Same as the `Config{Lindbladian}` method.

# Returns
`KrylovGapResult` with converted Lindbladian eigenvalues, spectral gap, fixed point,
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
    krylov_kwargs...
)
    # Guards
    krylovdim > howmany || error("krylovdim ($krylovdim) must be > howmany ($howmany)")
    _check_krylov_memory(config.num_qubits, krylovdim)

    # Get delta from config
    delta = config.delta

    # Allocate workspace using Config{Thermalize} constructor (precomputes channel matrices)
    ws = KrylovWorkspace(config, hamiltonian, jumps; trotter=trotter)

    # Dimensions
    dim = size(hamiltonian.data, 1)

    # Build channel matvec closure using faithful Chen channel
    function channel_matvec(v::AbstractVector)
        rho = reshape(v, dim, dim)
        apply_delta_channel!(ws, rho, config, hamiltonian)
        return copy(vec(ws.rho_out))  # CRITICAL: copy to avoid aliasing
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

    return KrylovGapResult{Float64}(
        Complex{Float64}.(eigenvalues_sorted),
        spectral_gap,
        Complex{Float64}.(fixed_point),
        Complex{Float64}.(gap_mode),
        info.converged,
        info.numops,
        info.numiter,
        normres,
        channel_eigenvalues_sorted,
        Float64(delta),
    )
end
