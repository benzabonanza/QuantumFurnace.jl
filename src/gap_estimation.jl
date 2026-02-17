# ============================================================================
# Gap Estimation API: single-call spectral gap estimation from trajectories
# ============================================================================
#
# Composes:
#   - build_gap_estimation_observables (convergence.jl, Phase 20)
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
        observables, observable_names = build_gap_estimation_observables(
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
# CrossValidationResult struct
# ---------------------------------------------------------------------------

"""
    CrossValidationResult

Result of cross-validating a trajectory-fitted spectral gap against an exact
Liouvillian eigenvalue.

Returned by [`cross_validate_gap`](@ref). Compares the fitted gap from
[`SpectralGapResult`](@ref) against `abs(real(spectral_gap))` from the exact
eigendecomposition, and warns when the imaginary part is significant.

# Fields
- `fitted_gap::Float64`: Gap from trajectory-based exponential fitting.
- `exact_gap::Float64`: Gap from exact eigenvalue, computed as `abs(real(eigenvalue))`.
- `relative_error::Float64`: `|fitted - exact| / exact`. `Inf` if `exact == 0`.
- `absolute_error::Float64`: `|fitted - exact|`.
- `within_ci::Bool`: Whether `exact_gap` falls within the fitted confidence interval.
- `imaginary_ratio::Float64`: `|Im(eigenvalue)| / |Re(eigenvalue)|`. `Inf` if `Re == 0`.
- `imaginary_warning::Bool`: `true` if `imaginary_ratio > 0.1`.
"""
struct CrossValidationResult
    fitted_gap::Float64
    exact_gap::Float64
    relative_error::Float64
    absolute_error::Float64
    within_ci::Bool
    imaginary_ratio::Float64
    imaginary_warning::Bool
end

# ---------------------------------------------------------------------------
# cross_validate_gap: compare trajectory-fitted gap against exact eigenvalue
# ---------------------------------------------------------------------------

"""
    cross_validate_gap(estimated::SpectralGapResult, exact_result::LindbladianResult) -> CrossValidationResult

Cross-validate a trajectory-fitted spectral gap against the exact Liouvillian
eigenvalue from [`LindbladianResult`](@ref).

Extracts `exact_result.spectral_gap` and delegates to the `Complex` method.

# Arguments
- `estimated::SpectralGapResult`: Result from [`estimate_spectral_gap`](@ref).
- `exact_result::LindbladianResult`: Result from exact eigendecomposition.

# Returns
A [`CrossValidationResult`](@ref) with error metrics and imaginary-part warning.
"""
function cross_validate_gap(estimated::SpectralGapResult, exact_result::LindbladianResult)
    return cross_validate_gap(estimated, exact_result.spectral_gap)
end

"""
    cross_validate_gap(estimated::SpectralGapResult, exact_eigenvalue::Complex) -> CrossValidationResult

Cross-validate a trajectory-fitted spectral gap against an exact eigenvalue
provided directly as a complex number.

The exact gap is computed as `abs(real(exact_eigenvalue))` (locked decision:
the trajectory exponential fit captures only the real decay rate).

Emits `@warn` when `|Im(eigenvalue)/Re(eigenvalue)| > 0.1`, indicating that
a pure exponential fit may not capture oscillatory decay.

# Arguments
- `estimated::SpectralGapResult`: Result from [`estimate_spectral_gap`](@ref).
- `exact_eigenvalue::Complex`: Exact Liouvillian eigenvalue (spectral gap).

# Returns
A [`CrossValidationResult`](@ref) with error metrics and imaginary-part warning.

# Example
```julia
cv = cross_validate_gap(estimated, -0.5 + 0.1im)
println("relative error = \$(cv.relative_error), within CI = \$(cv.within_ci)")
```
"""
function cross_validate_gap(estimated::SpectralGapResult, exact_eigenvalue::Complex)
    exact_gap = abs(real(exact_eigenvalue))
    fitted_gap = estimated.gap

    # Handle edge case: Re(eigenvalue) == 0
    if real(exact_eigenvalue) == 0.0
        im_ratio = Inf
        relative_error = Inf
    else
        im_ratio = abs(imag(exact_eigenvalue)) / abs(real(exact_eigenvalue))
        relative_error = abs(fitted_gap - exact_gap) / exact_gap
    end

    absolute_error = abs(fitted_gap - exact_gap)
    within_ci = estimated.gap_ci[1] <= exact_gap <= estimated.gap_ci[2]
    im_warning = im_ratio > 0.1

    if im_warning
        @warn "Exact eigenvalue has significant imaginary part (|Im/Re| = $(round(im_ratio; digits=3))). " *
              "Pure exponential fit may not capture oscillatory decay." exact_eigenvalue=exact_eigenvalue
    end

    return CrossValidationResult(
        fitted_gap, exact_gap, relative_error, absolute_error,
        within_ci, im_ratio, im_warning,
    )
end
