# ============================================================================
# Gap Estimation API: single-call spectral gap estimation from trajectories
# ============================================================================
#
# Composes:
#   - build_preset_trajectory_observables (convergence.jl)
#   - run_observable_trajectories (trajectories.jl, Phase 22)
#   - fit_exponential_decay (fitting.jl, Phase 21)
# into a single estimate_spectral_gap function.

# ---------------------------------------------------------------------------
# SpectralGapResult struct
# ---------------------------------------------------------------------------

"""
    SpectralGapResult

Result of spectral gap estimation via trajectory-based exponential decay fitting.

Returned by [`estimate_spectral_gap`](@ref). Contains the best gap estimate,
per-observable fit details, and metadata for reproducibility.

# Fields
- `gap::Float64`: Gap estimate from the best-fit observable.
- `gap_ci::Tuple{Float64, Float64}`: 95% confidence interval on the gap.
- `gap_se::Float64`: Standard error on the gap.
- `best_observable::String`: Name of the observable selected as best.
- `best_r_squared::Float64`: R-squared of the best fit.
- `per_observable::Vector{FitResult}`: All fits, one per observable.
- `observable_names::Vector{String}`: Names matching `per_observable` order.
- `ntraj::Int`: Number of trajectories used.
- `total_time::Float64`: Total simulation time.
- `save_every::Int`: Save interval in steps.
- `seed::Int`: RNG seed used (for reproducibility).
- `skip_initial::Float64`: Fraction of initial data skipped for fitting.
"""
struct SpectralGapResult
    gap::Float64
    gap_ci::Tuple{Float64, Float64}
    gap_se::Float64
    best_observable::String
    best_r_squared::Float64
    per_observable::Vector{FitResult}
    observable_names::Vector{String}
    ntraj::Int
    total_time::Float64
    save_every::Int
    seed::Int
    skip_initial::Float64
end

# ---------------------------------------------------------------------------
# Internal helper: best-observable selection (GAP-03)
# ---------------------------------------------------------------------------

"""
    _select_best_observable(fits, names) -> (best_idx, best_name, best_r_squared)

Select the best observable for spectral gap estimation. The true spectral gap
is the slowest-decaying mode, so among fits with acceptable quality we pick the
one with the **smallest positive gap** (not the highest R-squared).

Selection criteria (in priority order):
1. **Primary:** Among fits where `converged && gap > 0 && R-squared > 0.8`,
   select the one with the smallest `gap`.
2. **Fallback 1:** If no fit passes the R-squared threshold, among fits where
   `converged && gap > 0`, select the one with the smallest `gap`.
3. **Fallback 2:** If still none, select the fit with the highest `R-squared`
   overall (for diagnostic purposes).
"""
function _select_best_observable(fits::Vector{FitResult}, names::Vector{String})
    # Primary: smallest gap among high-quality fits (converged, gap > 0, R² > 0.8)
    best_idx = 0
    best_gap = Inf

    for (i, fit) in enumerate(fits)
        if fit.converged && fit.gap > 0.0 && fit.r_squared > 0.8 && fit.gap < best_gap
            best_idx = i
            best_gap = fit.gap
        end
    end

    # Fallback 1: smallest gap among converged fits with positive gap (any R²)
    if best_idx == 0
        for (i, fit) in enumerate(fits)
            if fit.converged && fit.gap > 0.0 && fit.gap < best_gap
                best_idx = i
                best_gap = fit.gap
            end
        end
    end

    # Fallback 2: highest R-squared regardless (diagnostic)
    if best_idx == 0
        best_idx = argmax(fit.r_squared for fit in fits)
    end

    return best_idx, names[best_idx], fits[best_idx].r_squared
end

# ---------------------------------------------------------------------------
# Main API function
# ---------------------------------------------------------------------------

"""
    estimate_spectral_gap(jumps, config, psi0, hamiltonian; kwargs...) -> SpectralGapResult

Estimate the Lindbladian spectral gap from trajectory simulations. This is
the single-call API that composes observable construction, trajectory
simulation, exponential fitting, and best-observable selection.

# Arguments
- `jumps::Vector{JumpOp}`: Jump operators.
- `config::AbstractThermalizeConfig`: Thermalization configuration.
- `psi0::Vector{<:Complex}`: Initial state vector.
- `hamiltonian::HamHam`: Hamiltonian data.

# Keyword Arguments
- `observables`: Custom observable matrices (default: `nothing`, uses H + Mz bundle).
- `observable_names`: Names for custom observables (required when `observables` is provided).
- `ntraj::Int=1000`: Number of trajectories.
- `save_every::Int=10`: Save interval in steps.
- `total_time::Real=config.mixing_time`: Total simulation time.
- `delta::Real=config.delta`: Time step.
- `seed::Union{Int, Nothing}=nothing`: RNG seed (random if `nothing`).
- `trotter::Union{TrottTrott, Nothing}=nothing`: Trotter object (required for TrotterDomain).
- `skip_initial::Float64=0.0`: Fraction of initial data to skip for fitting.

# Returns
A [`SpectralGapResult`](@ref) containing the gap estimate, confidence interval,
per-observable fits, and simulation metadata.

# Example
```julia
psi0 = zeros(ComplexF64, dim); psi0[1] = 1.0
result = estimate_spectral_gap(jumps, config, psi0, hamiltonian;
    ntraj=1000, save_every=10, seed=42)
println("gap = \$(result.gap), best observable = \$(result.best_observable)")
```
"""
function estimate_spectral_gap(
    jumps::Vector{JumpOp},
    config::AbstractThermalizeConfig,
    psi0::Vector{<:Complex},
    hamiltonian::HamHam;
    observables::Union{Nothing, Vector{<:Matrix{<:Complex}}} = nothing,
    observable_names::Union{Nothing, Vector{String}} = nothing,
    ntraj::Int = 1000,
    save_every::Int = 10,
    total_time::Real = config.mixing_time,
    delta::Real = config.delta,
    seed::Union{Int, Nothing} = nothing,
    trotter::Union{TrottTrott, Nothing} = nothing,
    skip_initial::Float64 = 0.0,
)
    # 1. Build default observables if not provided
    if observables === nothing
        observables, observable_names = build_preset_trajectory_observables(
            hamiltonian, config.num_qubits; trotter=trotter)
    end

    # Validate: names must be provided with custom observables
    observable_names === nothing && throw(ArgumentError(
        "observable_names must be provided when custom observables are given"))
    length(observables) == length(observable_names) || throw(ArgumentError(
        "observables and observable_names must have same length " *
        "(got $(length(observables)) observables and $(length(observable_names)) names)"))

    # 2. Run observable-only trajectories
    traj_result = run_observable_trajectories(
        jumps, config, psi0, hamiltonian;
        observables=observables, save_every=save_every,
        ntraj=ntraj, total_time=total_time, delta=delta,
        seed=seed, trotter=trotter,
    )

    # 3. Fit each observable's time series
    fits = FitResult[]
    for i in eachindex(observables)
        obs_series = Float64.(traj_result.measurements_mean[i, :])
        fit = fit_exponential_decay(traj_result.times, obs_series;
                                     skip_initial=skip_initial)
        push!(fits, fit)
    end

    # 4. Select best observable (GAP-03)
    best_idx, best_name, best_r2 = _select_best_observable(fits, observable_names)
    best_fit = fits[best_idx]

    # 5. Build and return result
    return SpectralGapResult(
        best_fit.gap,
        best_fit.gap_ci,
        best_fit.gap_se,
        best_name,
        best_r2,
        fits,
        observable_names,
        traj_result.n_trajectories,
        Float64(total_time),
        save_every,
        traj_result.seed,
        skip_initial,
    )
end

# ---------------------------------------------------------------------------
# OverlapAnalysisResult struct
# ---------------------------------------------------------------------------

"""
    OverlapAnalysisResult

Result of eigenbasis overlap analysis for observable-Lindbladian coupling diagnostics.

Returned by [`eigenbasis_overlap_analysis`](@ref). Decomposes each observable into
the Lindbladian eigenbasis to quantify how strongly each observable couples to the
spectral gap mode (first excited eigenmode).

# Fields
- `eigenvalues::Vector{ComplexF64}`: Liouvillian eigenvalues sorted by |Re(lambda)|.
- `exact_gap::Float64`: `abs(real(eigenvalues[2]))` -- the exact spectral gap.
- `observable_names::Vector{String}`: Names of the observables analyzed.
- `overlap_coefficients::Matrix{ComplexF64}`: `n_obs x n_modes` matrix of expansion
  coefficients `c_k = <O|v_k> * alpha_k` where `v_k` is the k-th eigenmode and
  `alpha_k` is the initial state projection onto mode k.
- `gap_mode_overlap::Vector{Float64}`: `|c_2|` for each observable (mode 2 = gap mode).
- `relative_gap_overlap::Vector{Float64}`: `|c_2| / sum(|c_k| for k>=2)` for each
  observable, excluding the steady-state mode.
"""
struct OverlapAnalysisResult
    eigenvalues::Vector{ComplexF64}
    exact_gap::Float64
    observable_names::Vector{String}
    overlap_coefficients::Matrix{ComplexF64}
    gap_mode_overlap::Vector{Float64}
    relative_gap_overlap::Vector{Float64}
end

# ---------------------------------------------------------------------------
# Eigenbasis overlap analysis
# ---------------------------------------------------------------------------

"""
    eigenbasis_overlap_analysis(L, observables, observable_names, rho0; rho_beta=nothing) -> OverlapAnalysisResult

Decompose observables into the Lindbladian eigenbasis and measure coupling to
the spectral gap mode.

This is the key diagnostic for gap estimation: an observable that does not overlap
with the first excited eigenmode of L cannot estimate the spectral gap, regardless
of trajectory count. Larger `gap_mode_overlap` means the observable is better suited
for gap estimation.

# Arguments
- `L::Matrix{<:Complex}`: Full dense Liouvillian superoperator matrix.
- `observables::Vector{<:Matrix{<:Complex}}`: Observable matrices (in the same basis as L).
- `observable_names::Vector{String}`: Names matching each observable.
- `rho0::Matrix{<:Complex}`: Initial density matrix (in the same basis as L).

# Keyword Arguments
- `rho_beta::Union{Hermitian, Nothing}=nothing`: Gibbs (steady) state for proper
  steady-state subtraction. When provided, uses the exact biorthogonal formula
  `c_k = Tr[O R_k] * Tr[L_k^dagger (rho_0 - rho_beta)]` with explicit left and
  right eigenvectors. When `nothing`, uses the original simplified formula
  `c_k = vec(O)^H * v_k * (V \\ vec(rho0))_k` (v1.3 behavior).

# Returns
An [`OverlapAnalysisResult`](@ref) containing eigenvalues, overlap coefficients,
and per-observable gap mode coupling metrics.

# Mathematical Details
When `rho_beta` is provided, the time-dependent expectation value is decomposed as
`<O>(t) - <O>_beta = sum_{k>=2} c_k * exp(lambda_k * t)` where
`c_k = Tr[O R_k] * Tr[L_k^dagger (rho_0 - rho_beta)]`, with `R_k` being right
eigenvectors and `L_k` being left eigenvectors of the Liouvillian.

When `rho_beta` is `nothing`, uses `<O>(t) = sum_k c_k * exp(lambda_k * t)` where
`c_k = vec(O)^H * v_k * alpha_k` and `alpha_k = (V \\ vec(rho0))_k`.

# Example
```julia
L = Matrix(liouv_result.liouvillian)
obs, names = build_preset_trajectory_observables(ham, n)
rho0 = psi0 * psi0'
result = eigenbasis_overlap_analysis(L, obs, names, rho0; rho_beta=ham.gibbs)
println("Gap mode overlap: ", result.gap_mode_overlap)
```
"""
function eigenbasis_overlap_analysis(
    L::Matrix{<:Complex},
    observables::Vector{<:Matrix{<:Complex}},
    observable_names::Vector{String},
    rho0::Matrix{<:Complex};
    rho_beta::Union{Hermitian, Nothing}=nothing,
)
    # Full dense eigendecomposition
    F = eigen(L)
    perm = sortperm(abs.(real.(F.values)))
    lambda = F.values[perm]
    V_right = F.vectors[:, perm]

    n_obs = length(observables)
    n_modes = length(lambda)
    coeffs = zeros(ComplexF64, n_obs, n_modes)

    if rho_beta !== nothing
        # Exact biorthogonal formula: c_k = Tr[O R_k] * Tr[L_k^dagger (rho_0 - rho_beta)]
        # Left eigenvectors: rows of V_right^{-1}, as columns: (V_right^{-1})^T (transpose, not adjoint)
        V_inv = inv(V_right)
        V_left = transpose(V_inv)  # Left eigenvectors as columns (transpose, not conjugate transpose)

        rho_diff = reshape(rho0, :) - reshape(Matrix(rho_beta), :)
        dim = isqrt(size(V_right, 1))

        for k in 1:n_modes
            R_k_mat = reshape(V_right[:, k], dim, dim)  # Right eigenvector as dim x dim operator
            L_k_vec = V_left[:, k]                       # Left eigenvector (vectorized)
            # Factor 2: Tr[L_k^dagger (rho_0 - rho_beta)] = dot(L_k_vec, rho_diff)
            # Julia's dot conjugates first arg: dot(a,b) = conj(a)' * b
            factor2 = dot(L_k_vec, rho_diff)
            for (i, O) in enumerate(observables)
                # Factor 1: Tr[O_i R_k] = tr(O_i * R_k_mat)
                factor1 = tr(O * R_k_mat)
                coeffs[i, k] = factor1 * factor2
            end
        end
    else
        # Original simplified formula (v1.3 backward compatibility)
        # c[i, k] = vec(O_i)^H * v_k * alpha_k
        alpha = V_right \ reshape(rho0, :)
        for (i, O) in enumerate(observables)
            O_vec = reshape(O, :)
            for k in 1:n_modes
                # dot(a, b) in Julia computes conj(a)' * b = sum(conj.(a) .* b)
                o_k = dot(O_vec, V_right[:, k])
                coeffs[i, k] = o_k * alpha[k]
            end
        end
    end

    exact_gap = abs(real(lambda[2]))

    # Gap mode overlap: |c_2| for each observable (index 2 = gap mode in sorted eigenvalues)
    gap_overlap = Float64[abs(coeffs[i, 2]) for i in 1:n_obs]

    # Relative gap overlap: |c_2| / sum(|c_k| for k >= 2), excluding steady state (k=1)
    rel_overlap = Float64[]
    for i in 1:n_obs
        denom = sum(abs.(coeffs[i, 2:end]))
        if denom > 0.0
            push!(rel_overlap, abs(coeffs[i, 2]) / denom)
        else
            push!(rel_overlap, 0.0)
        end
    end

    return OverlapAnalysisResult(
        lambda, exact_gap, observable_names, coeffs, gap_overlap, rel_overlap,
    )
end
