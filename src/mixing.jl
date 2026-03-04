# ============================================================================
# Mixing Time Estimation: post-processing of trace distance convergence curves
# ============================================================================
#
# Given a ThermalizeResults (from run_thermalize), fits a decay model to the
# trace distance data and extracts the effective spectral gap, mixing time,
# and quality metrics.
#
# Supports two models:
#   :single  — d(t) = A * exp(-gap * t) + C  (default, backward compatible)
#   :biexp   — d(t) = A1 * exp(-g1 * t) + A2 * exp(-g2 * t) + C
#
# This is strictly post-processing -- it does NOT control or modify
# run_thermalize execution.

# ---------------------------------------------------------------------------
# MixingTimeEstimate struct
# ---------------------------------------------------------------------------

"""
    MixingTimeEstimate

Result of mixing time estimation from a `ThermalizeResults` trace distance curve.

Returned by [`estimate_mixing_time`](@ref). Contains the fitted spectral gap,
mixing time (actual and/or extrapolated), quality metrics, and the full
[`FitResult`](@ref) for advanced inspection.

# Fields
## Fit parameters (from FitResult slow mode)
- `fitted_gap::Float64`: Fitted decay rate (spectral gap estimate).
- `amplitude::Float64`: Fitted amplitude A (slow mode for biexp).
- `offset::Float64`: Fitted offset C (asymptotic value).
- `gap_ci::Tuple{Float64, Float64}`: Confidence interval on gap (lower, upper).
- `gap_se::Float64`: Standard error on gap.
- `r_squared::Float64`: Goodness of fit (1 - RSS/TSS).
- `converged::Bool`: Whether Levenberg-Marquardt optimization converged.

## Mixing time results
- `mixing_time::Float64`: Primary mixing time answer. If `extrapolate=true`, this is the
  extrapolated time (or `NaN` if extrapolation failed). If `target_epsilon` is provided
  without extrapolation, this is the actual time the data first crossed the target
  (or `NaN` if not reached). If neither, this is the total simulation time.
- `mixing_time_extrapolated::Union{Nothing, Float64}`: Extrapolated mixing time from
  the fitted model, or `nothing` if `extrapolate=false`.
- `mixing_time_actual::Union{Nothing, Float64}`: First time the trace distance
  crossed `target_epsilon` in the data, or `nothing` if not reached or not requested.
- `target_epsilon::Union{Nothing, Float64}`: The target trace distance used.

## Full fit result
- `fit_result::FitResult`: Complete single-exp fit result (or synthetic from biexp slow mode).

## Model info
- `model_used::Symbol`: `:single` or `:biexp`.
- `biexp_fit_result::Union{Nothing, BiexpFitResult}`: Full biexp fit, or `nothing` for single.
"""
struct MixingTimeEstimate
    # Fit parameters (from FitResult)
    fitted_gap::Float64
    amplitude::Float64
    offset::Float64
    gap_ci::Tuple{Float64, Float64}
    gap_se::Float64
    r_squared::Float64
    converged::Bool
    # Mixing time results
    mixing_time::Float64
    mixing_time_extrapolated::Union{Nothing, Float64}
    mixing_time_actual::Union{Nothing, Float64}
    target_epsilon::Union{Nothing, Float64}
    # Full fit for advanced users
    fit_result::FitResult
    # Model info (new in Phase 43)
    model_used::Symbol
    biexp_fit_result::Union{Nothing, BiexpFitResult}
end

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

"""
    _check_fit_quality(fit::FitResult, target_epsilon)

Issue `@warn` messages for quality gate violations. Does not throw.
"""
function _check_fit_quality(fit::FitResult, target_epsilon::Union{Nothing, Float64})
    if fit.r_squared < 0.95
        @warn "Fit R-squared = $(fit.r_squared) < 0.95. Single-exponential model may not describe the data well."
    end
    if target_epsilon !== nothing && fit.offset > 0.1 * target_epsilon
        @warn "Fit offset C = $(fit.offset) is large relative to target epsilon = $(target_epsilon). Extrapolation may be unreliable."
    end
    if !fit.converged
        @warn "Levenberg-Marquardt optimization did not converge. Fit results are unreliable."
    end
    if fit.gap_se > 0.5 * fit.gap && isfinite(fit.gap_se)
        @warn "Gap standard error ($(fit.gap_se)) exceeds 50% of fitted gap ($(fit.gap)). Confidence interval is very wide."
    end
    return nothing
end

"""
    _find_actual_mixing_time(times, dists, target_epsilon)

Find the first time in `times` where the trace distance `dists` drops to or
below `target_epsilon`. Returns `nothing` if `target_epsilon` is `nothing` or
if the target was never reached.
"""
function _find_actual_mixing_time(
    times::AbstractVector{<:Real},
    dists::AbstractVector{<:Real},
    target_epsilon::Union{Nothing, Float64},
)
    target_epsilon === nothing && return nothing
    for i in eachindex(dists)
        if dists[i] <= target_epsilon
            return Float64(times[i])
        end
    end
    return nothing  # target not reached in data
end

"""
    _extrapolate_mixing_time(fit::FitResult, target_epsilon)

Compute the extrapolated mixing time from the fitted model
`d(t) = A * exp(-gap * t) + C` by solving for `t` when `d(t) = epsilon`.

Returns `nothing` if extrapolation is not possible (missing target, non-positive
gap or amplitude, or offset exceeding target).
"""
function _extrapolate_mixing_time(fit::FitResult, target_epsilon::Union{Nothing, Float64})
    target_epsilon === nothing && return nothing
    fit.gap <= 0.0 && return nothing
    fit.amplitude <= 0.0 && return nothing

    # d(t) = A * exp(-gap * t) + C = epsilon
    # => A * exp(-gap * t) = epsilon - C
    # => t = -ln((epsilon - C) / A) / gap
    effective_target = target_epsilon - fit.offset
    effective_target <= 0.0 && return nothing   # offset exceeds target
    effective_target >= fit.amplitude && return nothing  # already below target at t=0

    return -log(effective_target / fit.amplitude) / fit.gap
end

"""
    _extrapolate_mixing_time_biexp(bifit::BiexpFitResult, target_epsilon)

Compute the extrapolated mixing time from the bi-exponential fitted model
`d(t) = A1*exp(-g1*t) + A2*exp(-g2*t) + C` by numerically solving
`d(t) = epsilon` using bisection via Roots.jl.

Returns `nothing` if extrapolation is not possible.
"""
function _extrapolate_mixing_time_biexp(bifit::BiexpFitResult, target_epsilon::Union{Nothing, Float64})
    target_epsilon === nothing && return nothing
    bifit.gap <= 0.0 && return nothing

    # Check: f(0) = A1 + A2 + C; need f(0) > target for a crossing to exist
    f0 = bifit.amplitude_fast + bifit.amplitude + bifit.offset
    f0 <= target_epsilon && return nothing  # already below target at t=0

    # Check: asymptotic value C must be below target
    bifit.offset >= target_epsilon && return nothing

    # Define the function to find root of: f(t) - epsilon = 0
    function biexp_residual(t)
        return bifit.amplitude_fast * exp(-bifit.gap_fast * t) +
               bifit.amplitude * exp(-bifit.gap * t) +
               bifit.offset - target_epsilon
    end

    # Upper bracket: use slow-mode estimate with 3x safety margin
    # From slow mode alone: t_slow = -ln((eps - C) / A_slow) / gap_slow
    eff = target_epsilon - bifit.offset
    if eff <= 0.0
        return nothing
    end
    t_slow_est = if bifit.amplitude > 0.0 && eff < bifit.amplitude
        -log(eff / bifit.amplitude) / bifit.gap
    else
        # Fallback: use a large upper bracket
        100.0 / bifit.gap
    end
    t_upper = max(t_slow_est * 3.0, 10.0 / bifit.gap)

    # Ensure f(t_upper) < target (bracket is valid)
    if biexp_residual(t_upper) > 0.0
        # Expand bracket
        t_upper *= 3.0
        if biexp_residual(t_upper) > 0.0
            return nothing  # cannot bracket
        end
    end

    try
        t_mix = Roots.find_zero(biexp_residual, (0.0, t_upper), Roots.Bisection())
        return t_mix
    catch
        return nothing
    end
end

"""
    _biexp_to_single_fit_result(bifit::BiexpFitResult) -> FitResult

Construct a synthetic `FitResult` from the slow-mode parameters of a
`BiexpFitResult` for backward compatibility with the `fit_result` field.
"""
function _biexp_to_single_fit_result(bifit::BiexpFitResult)
    return FitResult(
        bifit.gap,          # gap (slow mode)
        bifit.amplitude,    # amplitude (slow mode)
        bifit.offset,       # offset
        bifit.gap_ci,       # gap_ci (slow mode)
        bifit.gap_se,       # gap_se (slow mode)
        bifit.r_squared,    # r_squared (from full bi-exp fit)
        bifit.converged,    # converged
        bifit.residuals,    # residuals (from full bi-exp fit)
        bifit.times_used,   # times_used
        bifit.values_used,  # values_used
    )
end

# ---------------------------------------------------------------------------
# Main API
# ---------------------------------------------------------------------------

"""
    estimate_mixing_time(result::ThermalizeResults; kwargs...) -> MixingTimeEstimate

Estimate the mixing time from a `ThermalizeResults` trace distance curve by
fitting a decay model to the data.

This is a post-processing function: it operates on a completed simulation
result and does NOT control `run_thermalize` execution.

# Arguments
- `result::ThermalizeResults`: Completed thermalization simulation with
  `time_steps` and `trace_distances` fields.

# Keyword Arguments
- `skip_initial::Float64=0.2`: Fraction of initial data to skip (in [0, 1)).
  Early-time transients often deviate from single-exponential behavior.
- `target_epsilon::Union{Nothing, Float64}=nothing`: Target trace distance
  for mixing time. Required when `extrapolate=true`.
- `extrapolate::Bool=false`: If `true`, compute the extrapolated mixing time
  from the fitted model. Requires `target_epsilon`.
- `level::Float64=0.95`: Confidence level for the gap confidence interval.
- `model::Symbol=:single`: Fitting model to use.
  - `:single` — Single-exponential `A * exp(-gap * t) + C` (default).
  - `:biexp` — Bi-exponential `A1 * exp(-g1 * t) + A2 * exp(-g2 * t) + C`.
    Captures multi-timescale dynamics for more accurate offset estimation
    and extrapolation when the target epsilon is close to the floor.

# Returns
A [`MixingTimeEstimate`](@ref) containing the fitted gap, mixing time(s),
quality metrics, and the full [`FitResult`](@ref).

# Quality Gates
Issues `@warn` when (for `:single` model):
- R-squared < 0.95 (poor fit)
- Offset C > 0.1 * target_epsilon (unreliable extrapolation)
- Levenberg-Marquardt did not converge
- Gap standard error > 50% of fitted gap

# Example
```julia
result = run_thermalize(config, hamiltonian, jumps, rho0)
est = estimate_mixing_time(result; skip_initial=0.2, target_epsilon=0.01, extrapolate=true)
println("Spectral gap: \$(est.fitted_gap)")
println("Mixing time (extrapolated): \$(est.mixing_time)")
println("R-squared: \$(est.r_squared)")

# Bi-exponential model for improved extrapolation near floor:
est_bi = estimate_mixing_time(result; model=:biexp, target_epsilon=1e-4, extrapolate=true)
println("Biexp offset: \$(est_bi.offset)")
```
"""
function estimate_mixing_time(
    result::ThermalizeResults;
    skip_initial::Float64 = 0.2,
    target_epsilon::Union{Nothing, Float64} = nothing,
    extrapolate::Bool = false,
    level::Float64 = 0.95,
    model::Symbol = :single,
)
    times = result.time_steps
    dists = result.trace_distances

    # --- Input validation ---
    length(times) >= 10 || throw(ArgumentError(
        "Need at least 10 data points for mixing time estimation (got $(length(times)))"))

    if extrapolate && target_epsilon === nothing
        throw(ArgumentError("target_epsilon required when extrapolate=true"))
    end

    model in (:single, :biexp) || throw(ArgumentError(
        "model must be :single or :biexp (got :$model)"))

    # --- Compute actual mixing time (model-independent) ---
    t_mix_actual = _find_actual_mixing_time(times, dists, target_epsilon)

    if model == :single
        # --- Single-exponential path (original behavior) ---
        fit = fit_exponential_decay(Float64.(times), Float64.(dists);
            skip_initial=skip_initial, level=level)

        _check_fit_quality(fit, target_epsilon)

        t_mix_extrap = extrapolate ? _extrapolate_mixing_time(fit, target_epsilon) : nothing

        mixing_time = if extrapolate
            t_mix_extrap !== nothing ? t_mix_extrap : NaN
        elseif target_epsilon !== nothing
            t_mix_actual !== nothing ? t_mix_actual : NaN
        else
            Float64(last(times))
        end

        return MixingTimeEstimate(
            fit.gap, fit.amplitude, fit.offset,
            fit.gap_ci, fit.gap_se, fit.r_squared, fit.converged,
            mixing_time, t_mix_extrap, t_mix_actual, target_epsilon,
            fit, :single, nothing,
        )

    else  # model == :biexp
        # --- Bi-exponential path ---
        bifit = fit_biexponential_decay(Float64.(times), Float64.(dists);
            skip_initial=skip_initial, level=level)

        # Use bi-exponential extrapolation
        t_mix_extrap = extrapolate ? _extrapolate_mixing_time_biexp(bifit, target_epsilon) : nothing

        mixing_time = if extrapolate
            t_mix_extrap !== nothing ? t_mix_extrap : NaN
        elseif target_epsilon !== nothing
            t_mix_actual !== nothing ? t_mix_actual : NaN
        else
            Float64(last(times))
        end

        # Construct synthetic FitResult for backward compat
        synthetic_fit = _biexp_to_single_fit_result(bifit)

        return MixingTimeEstimate(
            bifit.gap, bifit.amplitude, bifit.offset,
            bifit.gap_ci, bifit.gap_se, bifit.r_squared, bifit.converged,
            mixing_time, t_mix_extrap, t_mix_actual, target_epsilon,
            synthetic_fit, :biexp, bifit,
        )
    end
end
