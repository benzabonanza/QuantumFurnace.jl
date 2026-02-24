# ============================================================================
# Krylov Eigensolver: matrix-free spectral gap computation via KrylovKit
# ============================================================================
#
# Wraps KrylovKit.eigsolve around the Krylov matvec infrastructure (Phase 27-28)
# to provide a single-call API for computing Lindbladian spectral gaps without
# constructing the full dense Liouvillian.
#
# Two dispatch paths via config type:
#   - AbstractLiouvConfig     -> Lindbladian eigsolve with :LR targeting
#   - AbstractThermalizeConfig -> CPTP channel eigsolve with :LM targeting
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
(`AbstractThermalizeConfig`), the raw channel eigenvalues and delta are also stored.

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
- `delta_used`: Delta from ThermalizeConfig (nothing for Lindbladian path).
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
# ThermalizeConfig -> LiouvConfig conversion
# ---------------------------------------------------------------------------

"""
    _thermalize_to_liouv_config(tc::ThermalizeConfig) -> LiouvConfig

Build a `LiouvConfig` from a `ThermalizeConfig` by copying all shared Lindbladian
parameters and stripping `mixing_time` and `delta`.

Required because `apply_lindbladian!` dispatches on `AbstractLiouvConfig`, not
`AbstractThermalizeConfig`.
"""
function _thermalize_to_liouv_config(tc::ThermalizeConfig)
    LiouvConfig(
        num_qubits = tc.num_qubits,
        with_coherent = tc.with_coherent,
        with_linear_combination = tc.with_linear_combination,
        domain = tc.domain,
        beta = tc.beta,
        sigma = tc.sigma,
        gaussian_parameters = tc.gaussian_parameters,
        a = tc.a,
        b = tc.b,
        num_energy_bits = tc.num_energy_bits,
        t0 = tc.t0,
        w0 = tc.w0,
        eta = tc.eta,
        num_trotter_steps_per_t0 = tc.num_trotter_steps_per_t0,
    )
end

"""
    _thermalize_to_liouv_config(tc::ThermalizeConfigGNS) -> LiouvConfigGNS

Build a `LiouvConfigGNS` from a `ThermalizeConfigGNS`. GNS configs always have
`with_coherent = false`.
"""
function _thermalize_to_liouv_config(tc::ThermalizeConfigGNS)
    LiouvConfigGNS(
        num_qubits = tc.num_qubits,
        with_coherent = false,
        with_linear_combination = tc.with_linear_combination,
        domain = tc.domain,
        beta = tc.beta,
        sigma = tc.sigma,
        gaussian_parameters = tc.gaussian_parameters,
        a = tc.a,
        b = tc.b,
        num_energy_bits = tc.num_energy_bits,
        t0 = tc.t0,
        w0 = tc.w0,
        eta = tc.eta,
        num_trotter_steps_per_t0 = tc.num_trotter_steps_per_t0,
    )
end

# ---------------------------------------------------------------------------
# apply_delta_channel!
# ---------------------------------------------------------------------------

"""
    apply_delta_channel!(ws, rho, delta, config_liouv, hamiltonian) -> ws.rho_out

Compute the CPTP channel action E(rho) = rho + delta * L(rho).

First computes L(rho) via `apply_lindbladian!` (writing to `ws.rho_out`),
then overwrites `ws.rho_out` with `rho + delta * ws.rho_out`.

# Arguments
- `ws::KrylovWorkspace{T}`: Pre-allocated workspace
- `rho::Matrix{T}`: Input density matrix (dim x dim)
- `delta::Real`: Time step from ThermalizeConfig
- `config_liouv::AbstractLiouvConfig`: Lindbladian configuration (NOT ThermalizeConfig)
- `hamiltonian::HamHam`: Hamiltonian

# Returns
`ws.rho_out` containing E(rho).
"""
function apply_delta_channel!(
    ws::KrylovWorkspace{T},
    rho::Matrix{T},
    delta::Real,
    config_liouv::AbstractLiouvConfig,
    hamiltonian::HamHam,
) where {T<:Complex}
    # L(rho) -> ws.rho_out
    apply_lindbladian!(ws, rho, config_liouv, hamiltonian)
    # E(rho) = rho + delta * L(rho) -> ws.rho_out
    @. ws.rho_out = rho + delta * ws.rho_out
    return ws.rho_out
end

# ---------------------------------------------------------------------------
# krylov_spectral_gap: Lindbladian path (AbstractLiouvConfig)
# ---------------------------------------------------------------------------

"""
    krylov_spectral_gap(config::AbstractLiouvConfig, hamiltonian, jumps; kwargs...) -> KrylovGapResult

Compute the Lindbladian spectral gap matrix-free using KrylovKit Arnoldi with `:LR` targeting.

Wraps the existing `apply_lindbladian!` matvec in a KrylovKit closure to find the leading
eigenvalues of the Lindbladian superoperator without constructing the full dense matrix.

The steady-state eigenvalue is near Re(lambda) ~ 0 (largest real part). The spectral gap
is `abs(real(lambda_2))` where lambda_2 is the second eigenvalue sorted by |Re(lambda)|.

# Arguments
- `config::AbstractLiouvConfig`: Lindbladian configuration (EnergyDomain, TimeDomain, etc.)
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
    config::AbstractLiouvConfig,
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
# krylov_spectral_gap: Channel path (AbstractThermalizeConfig)
# ---------------------------------------------------------------------------

"""
    krylov_spectral_gap(config::AbstractThermalizeConfig, hamiltonian, jumps; kwargs...) -> KrylovGapResult

Compute the Lindbladian spectral gap via the CPTP channel E(rho) = rho + delta * L(rho),
using KrylovKit Arnoldi with `:LM` targeting.

The channel eigenvalues mu are related to Lindbladian eigenvalues by the exact linear
formula: `lambda_L = (mu - 1) / delta`. The steady state has mu ~ 1 (largest magnitude),
and the gap is recovered from the second eigenvalue after conversion.

# Arguments
- `config::AbstractThermalizeConfig`: Thermalization configuration (provides delta)
- `hamiltonian::HamHam`: Hamiltonian with eigenbasis data
- `jumps::Vector{JumpOp}`: Jump operators

# Keyword Arguments
Same as the `AbstractLiouvConfig` method.

# Returns
`KrylovGapResult` with converted Lindbladian eigenvalues, spectral gap, fixed point,
gap mode, and the raw channel eigenvalues stored in `channel_eigenvalues`.
"""
function krylov_spectral_gap(
    config::AbstractThermalizeConfig,
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

    # Convert ThermalizeConfig -> LiouvConfig for apply_lindbladian! dispatch
    config_liouv = _thermalize_to_liouv_config(config)

    _check_krylov_memory(config.num_qubits, krylovdim)

    # Get delta from config
    delta = config.delta

    # Allocate workspace using the LiouvConfig
    ws = KrylovWorkspace(config_liouv, hamiltonian, jumps; trotter=trotter)

    # Dimensions
    dim = size(hamiltonian.data, 1)

    # Build channel matvec closure
    function channel_matvec(v::AbstractVector)
        rho = reshape(v, dim, dim)
        apply_delta_channel!(ws, rho, delta, config_liouv, hamiltonian)
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
