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
    estimate_mixing_time(times, distances; kwargs...) -> MixingTimeEstimate

Estimate the mixing time from raw `(times, distances)` vectors by fitting a
decay model and (optionally) extrapolating to a target epsilon.

This is the primary entry point for integrator-based pipelines (qf-lkb.2):
the output of [`lindblad_action_integrate`](@ref) and
[`discriminant_action_integrate`](@ref) provides `(t, distances)` directly,
and a NamedTuple forwarder method dispatches onto this one. The
[`ThermalizeResults`](@ref) overload also delegates here.

The extrapolator is metric-agnostic — `distances` may be trace distance
(L-mode), Frobenius distance / `χ` (K-mode), or any non-negative monotonically
decaying observable. The returned `mixing_time` is in the native metric of
the supplied data.

# Arguments
- `times::AbstractVector{<:Real}`: Time points (must match length of `distances`).
- `distances::AbstractVector{<:Real}`: Distance-to-equilibrium values at each time.

# Keyword Arguments
- `skip_initial::Real=0.2`: Fraction of initial data to skip (in [0, 1)).
  Early-time transients often deviate from single-exponential behavior.
- `target_epsilon::Union{Nothing, Real}=nothing`: Target distance for mixing
  time. Required when `extrapolate=true`.
- `extrapolate::Bool=false`: If `true`, compute the extrapolated mixing time
  from the fitted model. Requires `target_epsilon`.
- `level::Real=0.95`: Confidence level for the gap confidence interval.
- `model::Symbol=:biexp`: Fitting model to use. Defaults to `:biexp` here
  (vs `:single` on the `ThermalizeResults` overload, which preserves
  pre-Phase-43 behaviour). Bi-exp matches multi-timescale Liouvillian
  dynamics and gives 26%->0.001% error vs single-exp on real data.
  - `:single` — Single-exponential `A * exp(-gap * t) + C`.
  - `:biexp` — Bi-exponential `A1 * exp(-g1 * t) + A2 * exp(-g2 * t) + C`.

# Returns
A [`MixingTimeEstimate`](@ref) containing the fitted gap, mixing time(s),
quality metrics, and the full [`FitResult`](@ref).

# Example
```julia
res = lindblad_action_integrate(L_apply!, rho_0, sigma_beta, t_grid)
est = estimate_mixing_time(res.t, res.distances;
                            target_epsilon = 1e-3, extrapolate = true)
println("τ_mix = \$(est.mixing_time) (model=\$(est.model_used))")
```
"""
function estimate_mixing_time(
    times::AbstractVector{<:Real},
    distances::AbstractVector{<:Real};
    skip_initial::Real = 0.2,
    target_epsilon::Union{Nothing, Real} = nothing,
    extrapolate::Bool = false,
    level::Real = 0.95,
    model::Symbol = :biexp,
)::MixingTimeEstimate
    # --- Input validation ---
    length(times) == length(distances) || throw(ArgumentError(
        "times and distances must have the same length (got $(length(times)) and $(length(distances)))"))
    length(times) >= 10 || throw(ArgumentError(
        "Need at least 10 data points for mixing time estimation (got $(length(times)))"))

    if extrapolate && target_epsilon === nothing
        throw(ArgumentError("target_epsilon required when extrapolate=true"))
    end

    model in (:single, :biexp) || throw(ArgumentError(
        "model must be :single or :biexp (got :$model)"))

    # --- Coerce kwargs to the internal helpers' Float64 contract ---
    skip_initial_f = Float64(skip_initial)
    level_f        = Float64(level)
    target_eps_f   = target_epsilon === nothing ? nothing : Float64(target_epsilon)

    # --- Compute actual mixing time (model-independent) ---
    t_mix_actual = _find_actual_mixing_time(times, distances, target_eps_f)

    if model == :single
        # --- Single-exponential path (original behavior) ---
        fit = fit_exponential_decay(Float64.(times), Float64.(distances);
            skip_initial=skip_initial_f, level=level_f)

        _check_fit_quality(fit, target_eps_f)

        t_mix_extrap = extrapolate ? _extrapolate_mixing_time(fit, target_eps_f) : nothing

        mixing_time = if extrapolate
            t_mix_extrap !== nothing ? t_mix_extrap : NaN
        elseif target_eps_f !== nothing
            t_mix_actual !== nothing ? t_mix_actual : NaN
        else
            Float64(last(times))
        end

        return MixingTimeEstimate(
            fit.gap, fit.amplitude, fit.offset,
            fit.gap_ci, fit.gap_se, fit.r_squared, fit.converged,
            mixing_time, t_mix_extrap, t_mix_actual, target_eps_f,
            fit, :single, nothing,
        )

    else  # model == :biexp
        # --- Bi-exponential path ---
        bifit = fit_biexponential_decay(Float64.(times), Float64.(distances);
            skip_initial=skip_initial_f, level=level_f)

        # Use bi-exponential extrapolation
        t_mix_extrap = extrapolate ? _extrapolate_mixing_time_biexp(bifit, target_eps_f) : nothing

        mixing_time = if extrapolate
            t_mix_extrap !== nothing ? t_mix_extrap : NaN
        elseif target_eps_f !== nothing
            t_mix_actual !== nothing ? t_mix_actual : NaN
        else
            Float64(last(times))
        end

        # Construct synthetic FitResult for backward compat
        synthetic_fit = _biexp_to_single_fit_result(bifit)

        return MixingTimeEstimate(
            bifit.gap, bifit.amplitude, bifit.offset,
            bifit.gap_ci, bifit.gap_se, bifit.r_squared, bifit.converged,
            mixing_time, t_mix_extrap, t_mix_actual, target_eps_f,
            synthetic_fit, :biexp, bifit,
        )
    end
end

"""
    estimate_mixing_time(result::ThermalizeResults; kwargs...) -> MixingTimeEstimate

Estimate the mixing time from a `ThermalizeResults` trace distance curve by
fitting a decay model to the data.

This is a post-processing function: it operates on a completed simulation
result and does NOT control `run_thermalize` execution. Internally delegates
to the `(times, distances)` vector method on
`(result.time_steps, result.trace_distances)`.

# Arguments
- `result::ThermalizeResults`: Completed thermalization simulation with
  `time_steps` and `trace_distances` fields.

# Keyword Arguments
- `skip_initial::Real=0.2`: Fraction of initial data to skip (in [0, 1)).
  Early-time transients often deviate from single-exponential behavior.
- `target_epsilon::Union{Nothing, Real}=nothing`: Target trace distance
  for mixing time. Required when `extrapolate=true`.
- `extrapolate::Bool=false`: If `true`, compute the extrapolated mixing time
  from the fitted model. Requires `target_epsilon`.
- `level::Real=0.95`: Confidence level for the gap confidence interval.
- `model::Symbol=:single`: Fitting model to use. Default `:single` here
  preserves pre-Phase-43 behaviour for existing callers; pass `:biexp`
  explicitly for the multi-timescale model recommended for thesis-quality
  τ_mix near the floor (matches the new vector-method default).
  - `:single` — Single-exponential `A * exp(-gap * t) + C`.
  - `:biexp` — Bi-exponential `A1 * exp(-g1 * t) + A2 * exp(-g2 * t) + C`.

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
    skip_initial::Real = 0.2,
    target_epsilon::Union{Nothing, Real} = nothing,
    extrapolate::Bool = false,
    level::Real = 0.95,
    model::Symbol = :single,
)::MixingTimeEstimate
    return estimate_mixing_time(
        result.time_steps, result.trace_distances;
        skip_initial    = skip_initial,
        target_epsilon  = target_epsilon,
        extrapolate     = extrapolate,
        level           = level,
        model           = model,
    )
end

"""
    estimate_mixing_time(integrator_result::NamedTuple; kwargs...) -> MixingTimeEstimate

Convenience forwarder for the output of the matrix-free Lindbladian-action
integrators [`lindblad_action_integrate`](@ref) and
[`discriminant_action_integrate`](@ref) (qf-lkb.1). Dispatches onto the
`(times, distances)` vector method using `integrator_result.t` and
`integrator_result.distances`.

The integrator NamedTuple contract guarantees the field names `.t` and
`.distances`; any NamedTuple lacking those fields raises a `KeyError` at
field access. Default `model=:biexp` flows through from the vector method.

!!! warning "Not for the Krylov trajectory predictors"
    This (bi)exponential curve fit is for SAMPLED trajectory curves —
    integrator outputs and [`ThermalizeResults`](@ref). It must **not** be used
    on [`predict_lindbladian_trajectory`](@ref) /
    [`predict_channel_trajectory`](@ref) output: those return an exact
    bi-orthogonal eigendecomposition, and their τ_mix(ε) comes from the
    closed-form bisection [`eigenmode_mixing_time`](@ref), which is exact on the
    captured subspace and robust on flat tails where the LM fit degenerates
    (qf-3uj). A predictor NamedTuple (carrying `:R_modes` / `:eigenvalues`) is
    rejected here with an `ArgumentError` pointing to `eigenmode_mixing_time`.

# Example
```julia
res = lindblad_action_integrate(L_apply!, rho_0, sigma_beta, t_grid)
est = estimate_mixing_time(res; target_epsilon = 1e-3, extrapolate = true)
println("τ_mix = \$(est.mixing_time) (matvecs = \$(res.total_matvecs))")
```
"""
function estimate_mixing_time(integrator_result::NamedTuple; kwargs...)::MixingTimeEstimate
    # qf-3uj: refuse a Krylov trajectory-predictor result. predict_*_trajectory
    # return the spectral data (`:R_modes` / `:eigenvalues` / `:rho_inf` /
    # `:sigma_beta`); the matrix-free integrators return only
    # (t, distances, *_final, total_matvecs, all_converged, states). Keying on
    # `:R_modes` / `:eigenvalues` separates them cleanly, so the curve fit can
    # never silently stand in for the exact bisection on the predictor path.
    if haskey(integrator_result, :R_modes) || haskey(integrator_result, :eigenvalues)
        throw(ArgumentError(
            "estimate_mixing_time received a Krylov trajectory-predictor result " *
            "(it carries the spectral fields :R_modes / :eigenvalues). The " *
            "(bi)exponential curve fit must not be used on " *
            "predict_lindbladian_trajectory / predict_channel_trajectory output — " *
            "call `eigenmode_mixing_time(traj, target_epsilon)` (exact closed-form " *
            "bisection on the Krylov spectrum) instead. The curve-fit estimator " *
            "remains the right tool for the matrix-free integrator outputs " *
            "(lindblad_action_integrate / discriminant_action_integrate) and " *
            "ThermalizeResults."))
    end
    return estimate_mixing_time(
        integrator_result.t, integrator_result.distances;
        kwargs...,
    )
end

# ---------------------------------------------------------------------------
# Eigenmode τ_mix (qf-e4y.2): closed-form bisection on the Krylov spectral
# decomposition. Drops the bi-exp curve fit on the :krylov route — the
# eigendecomposition built by `predict_*_trajectory` already encodes the
# slow-mode amplitudes exactly, so the trace distance d(t) = ε equation
# reduces to a 1D bracketed bisection. Robust to LM degenerate-basin
# failures that haunt bi-exp on flat-tail cells.
# ---------------------------------------------------------------------------

const _EIGENMODE_ZERO_TOL = 1e-10  # |λ| below this counts as the steady mode

"""
    eigenmode_mixing_time(eigenvalues, c, R_modes, rho_inf, sigma_beta,
                          target_epsilon; t_upper, atol, max_iters,
                          eigenvalue_zero_tol)
        -> NamedTuple

Closed-form τ_mix(ε) from a Krylov bi-orthogonal eigendecomposition of `L`.

Bisects ``d(t) = \\| (\\rho_\\infty - \\sigma_\\beta) + \\sum_i c_i e^{\\lambda_i t} R_i \\|_1 / 2 = \\varepsilon``
on the continuous-time axis. The smallest ``|\\Re(\\lambda_i)|`` over the
non-steady eigenvalues (the steady mode `λ_1 ≈ 0` has `c_1 ≈ 0` by trace
preservation) is the spectral gap, which sets the bisection bracket.

Inputs are exactly the spectral data exposed by `predict_lindbladian_trajectory`
(and the channel analogue after `λ_eff = log(μ) / δ` conversion at the call
site).

# Arguments
- `eigenvalues::AbstractVector{<:Complex}`: Lindbladian eigenvalues. Steady
  mode (`|λ_1| < eigenvalue_zero_tol`) is automatically excluded from the
  gap and kept in the residual sum (with `c_1 ≈ 0` from the engine).
- `c::AbstractVector{<:Complex}`: biorthogonal coefficients
  ``c_i = \\langle L_i, \\rho_0 - \\rho_\\infty \\rangle_{HS}``.
- `R_modes::AbstractVector{<:AbstractMatrix}`: right-eigenvector matrices
  (one `d × d` complex matrix per captured Krylov mode).
- `rho_inf::AbstractMatrix`: captured steady state of `L` (the leading
  `R_modes[1]`, normalised to trace 1 by the engine).
- `sigma_beta::AbstractMatrix`: trace-distance reference (Gibbs state for
  the Lindbladian path; basis-aligned Gibbs for the channel path).
- `target_epsilon::Real`: the trace distance threshold ``\\varepsilon``.

# Keywords
- `t_upper::Real = 0.0`: bisection upper bracket. `0` triggers the heuristic
  `t_upper = max(10/gap, 50)` followed by one 3× expansion if needed.
- `atol::Real = 1e-3`: bisection tolerance (in `t` units). The bracket
  shrinks until `|t_high - t_low| < atol`.
- `max_iters::Int = 64`: hard cap on bisection steps (passed through to
  `Roots.find_zero`).
- `eigenvalue_zero_tol::Real = 1e-10`: threshold below which an eigenvalue
  is treated as the steady mode (excluded from the gap, included in the
  residual sum since the engine sets `c_1 = 0`).

# Returns
NamedTuple with:
- `mixing_time::Float64`: ``\\tau_{\\text{mix}}(\\varepsilon)``, or `Inf` if
  ``\\varepsilon`` is below the asymptotic floor
  ``\\| \\rho_\\infty - \\sigma_\\beta \\|_1 / 2`` (no crossing exists).
- `gap::Float64`: smallest ``|\\Re(\\lambda_i)|`` over non-steady modes.
- `floor_distance::Float64`: ``d(\\infty) = \\| \\rho_\\infty - \\sigma_\\beta \\|_1 / 2``.
- `source::Symbol`: `:extrapolated` (bisection found a crossing),
  `:floor` (target below floor; `mixing_time = Inf`), `:nan` (degenerate
  input — fewer than 2 eigenvalues, or the bracket failed even after one
  3× expansion).
- `n_evals::Int`: number of `d(t)` evaluations performed.

# Notes
- Defensive Hermitisation is applied to the residual matrix before
  `svdvals`, mirroring the existing predictor loop's convention.
- Complex eigenvalues come in conjugate pairs `(λ, λ̄)` with paired
  `(c, c̄)` and `(R, R̄)`, so the residual matrix is Hermitian to
  machine precision; the Hermitisation only suppresses round-off.
- The closed-form formula is exact on the Krylov-captured subspace.
  Truncation error mirrors that of the predictor itself — `all_converged`
  on the predictor flags it.
"""
function eigenmode_mixing_time(
    eigenvalues::AbstractVector{<:Complex},
    c::AbstractVector{<:Complex},
    R_modes::AbstractVector{<:AbstractMatrix},
    rho_inf::AbstractMatrix,
    sigma_beta::AbstractMatrix,
    target_epsilon::Real;
    t_upper::Real = 0.0,
    atol::Real = 1e-3,
    max_iters::Int = 64,
    eigenvalue_zero_tol::Real = _EIGENMODE_ZERO_TOL,
)::NamedTuple
    h = length(eigenvalues)
    h == length(c) == length(R_modes) || throw(ArgumentError(
        "eigenvalues, c, R_modes must have the same length (got $h, $(length(c)), $(length(R_modes)))"))
    target_epsilon > 0.0 || throw(ArgumentError(
        "target_epsilon must be positive (got $target_epsilon)"))

    # Floor: d(∞) = ‖rho_inf - sigma_beta‖_1 / 2 (in the Hermitian residual,
    # half the sum of singular values is the trace distance).
    floor_distance = sum(svdvals(rho_inf .- sigma_beta)) / 2

    # Spectral gap: smallest |Re(λ_i)| over non-steady modes.
    gap = Inf
    for i in 1:h
        abs(eigenvalues[i]) < eigenvalue_zero_tol && continue
        gi = abs(real(eigenvalues[i]))
        gi < gap && (gap = gi)
    end

    # Degenerate input: no non-steady modes captured.
    if !isfinite(gap) || gap <= 0.0
        return (
            mixing_time     = NaN,
            gap             = isfinite(gap) ? gap : NaN,
            floor_distance  = floor_distance,
            source          = :nan,
            n_evals         = 0,
        )
    end

    # Floor branch: target below asymptotic floor → no crossing.
    if floor_distance >= target_epsilon
        return (
            mixing_time     = Inf,
            gap             = gap,
            floor_distance  = floor_distance,
            source          = :floor,
            n_evals         = 0,
        )
    end

    d = size(rho_inf, 1)
    T = promote_type(eltype(rho_inf), eltype(sigma_beta), ComplexF64)
    rho_t = Matrix{T}(undef, d, d)
    floor_residual = Matrix{T}(rho_inf .- sigma_beta)
    eval_count = Ref(0)

    function d_at(t::Real)::Float64
        copyto!(rho_t, floor_residual)
        @inbounds for i in 1:h
            abs(eigenvalues[i]) < eigenvalue_zero_tol && continue
            phase = exp(eigenvalues[i] * t)
            rho_t .+= (c[i] * phase) .* R_modes[i]
        end
        # Defensive Hermitisation (mirrors predict_*_trajectory loop).
        @inbounds for j in 1:d, k in 1:d
            rho_t[k, j] = (rho_t[k, j] + conj(rho_t[j, k])) / 2
        end
        eval_count[] += 1
        return sum(svdvals(rho_t)) / 2
    end

    # Bisection bracket. The eigenmode formula is monotonically decreasing in
    # t on the slow tail, so [0, t_upper] with d(t_upper) < ε brackets a root.
    t_hi = float(t_upper)
    if t_hi <= 0.0
        t_hi = max(10.0 / gap, 50.0)
    end

    d_hi = d_at(t_hi)
    if d_hi > target_epsilon
        # One 3× expansion. The slow mode decays like e^{-gap * t}, so
        # log(d(t)/d(t_hi)) ≈ -gap (t - t_hi); 3× is enough margin in
        # practice for any well-behaved spectrum.
        t_hi *= 3.0
        d_hi = d_at(t_hi)
        if d_hi > target_epsilon
            return (
                mixing_time     = NaN,
                gap             = gap,
                floor_distance  = floor_distance,
                source          = :nan,
                n_evals         = eval_count[],
            )
        end
    end

    # d(0) ≥ target_epsilon must hold for a crossing on [0, t_hi]; if not,
    # the trajectory is below target at t=0 — degenerate (or trivially mixed).
    d_zero = d_at(0.0)
    if d_zero <= target_epsilon
        return (
            mixing_time     = 0.0,
            gap             = gap,
            floor_distance  = floor_distance,
            source          = :extrapolated,
            n_evals         = eval_count[],
        )
    end

    residual(t::Float64) = d_at(t) - float(target_epsilon)

    t_mix = try
        # `atol` is the t-axis tolerance ⇒ pass as `xatol` to Roots (their
        # `atol` kw governs f(x)=d(t)-ε, which would couple to the slow
        # mode's amplitude and hide stride O(ε / |d'|) ≫ atol).
        Roots.find_zero(residual, (0.0, t_hi), Roots.Bisection();
                         xatol=float(atol), maxiters=max_iters)
    catch
        return (
            mixing_time     = NaN,
            gap             = gap,
            floor_distance  = floor_distance,
            source          = :nan,
            n_evals         = eval_count[],
        )
    end

    return (
        mixing_time     = float(t_mix),
        gap             = gap,
        floor_distance  = floor_distance,
        source          = :extrapolated,
        n_evals         = eval_count[],
    )
end

"""
    eigenmode_mixing_time(traj::NamedTuple, target_epsilon; kwargs...) -> NamedTuple

Closed-form τ_mix(ε) directly from a Krylov trajectory-predictor result.

This is the **canonical** τ_mix extraction for
[`predict_lindbladian_trajectory`](@ref) and
[`predict_channel_trajectory`](@ref): it reads the bi-orthogonal
eigendecomposition the predictor already returned (`eigenvalues`, `c`,
`R_modes`, `rho_inf`, `sigma_beta`) and bisects `d(t) = ε` on it. Always prefer
this over a (bi)exponential curve fit of the sampled `(t, distances)` curve
(`estimate_mixing_time`) — the bisection is exact on the captured Krylov
subspace and robust on flat tails where the Levenberg–Marquardt fit degenerates
(qf-3uj). `estimate_mixing_time` actively refuses a predictor NamedTuple for
this reason.

For a **channel** trajectory (detected by the `delta_used` field) the channel
eigenvalues `μ` are converted to Lindbladian rates `λ_eff = log(μ) / δ`, so the
`exp(λ t)` machinery runs in physical time `t = k·δ` and matches the trajectory
curve (steady `μ ≈ 1 → λ ≈ 0`). For a **Lindbladian** trajectory the
eigenvalues are already rates and are used as-is.

`kwargs` (`t_upper`, `atol`, `max_iters`, `eigenvalue_zero_tol`) forward
unchanged to the spectral-data method.

# Example
```julia
traj = predict_lindbladian_trajectory(cfg, ham, jumps, rho_0, t_grid; krylovdim=40)
res  = eigenmode_mixing_time(traj, 1e-3)            # bisection, not a curve fit
println("τ_mix = \$(res.mixing_time)  (source=\$(res.source))")

ch   = predict_channel_trajectory(cfg_C, ham, jumps_C, rho_0, k_grid; trotter=trot)
resC = eigenmode_mixing_time(ch, 1e-3)              # μ → λ_eff conversion is automatic
```
"""
function eigenmode_mixing_time(traj::NamedTuple, target_epsilon::Real; kwargs...)::NamedTuple
    for f in (:eigenvalues, :c, :R_modes, :rho_inf, :sigma_beta)
        haskey(traj, f) || throw(ArgumentError(
            "eigenmode_mixing_time(traj::NamedTuple, …) expects a Krylov " *
            "trajectory-predictor result with fields :eigenvalues, :c, :R_modes, " *
            ":rho_inf, :sigma_beta (missing :$f). Pass the output of " *
            "predict_lindbladian_trajectory / predict_channel_trajectory."))
    end
    eigenvalues = if haskey(traj, :delta_used)
        # Channel: μ → λ_eff = log(μ)/δ (mirrors the predict_channel_trajectory
        # convention; at t = k·δ, exp(λ_eff·t) = μ^k exactly).
        δ = traj.delta_used
        ComplexF64[log(complex(μ)) / δ for μ in traj.eigenvalues]
    else
        ComplexF64.(traj.eigenvalues)
    end
    return eigenmode_mixing_time(
        eigenvalues, traj.c, traj.R_modes, traj.rho_inf, traj.sigma_beta,
        target_epsilon; kwargs...)
end
