# ============================================================================
# Exponential Decay Fitting: extract spectral gap from time-series data
# ============================================================================
#
# Model: y(t) = A * exp(-gap * t) + C
#
# Parameter ordering convention:
#   p[1] = A       (amplitude)
#   p[2] = gap     (decay rate / spectral gap)
#   p[3] = C       (offset / asymptotic value)

const _IDX_A   = 1
const _IDX_GAP = 2
const _IDX_C   = 3

# Exponential decay model for LsqFit (module-level, not inside any function)
_exp_decay_model(t, p) = @. p[_IDX_A] * exp(-p[_IDX_GAP] * t) + p[_IDX_C]

# ---------------------------------------------------------------------------
# FitResult struct
# ---------------------------------------------------------------------------

"""
    FitResult

Result of single-exponential decay fitting: `A * exp(-gap * t) + C`.

Returned by [`fit_exponential_decay`](@ref). All fitted parameters and
quality metrics are stored; the downstream consumer does not need access to
the raw `LsqFit.LsqFitResult`.

# Fields
- `gap::Float64`: Fitted decay rate (spectral gap estimate), constrained > 0.
- `amplitude::Float64`: Fitted amplitude A.
- `offset::Float64`: Fitted offset C (asymptotic value).
- `gap_ci::Tuple{Float64, Float64}`: Confidence interval on gap (lower, upper).
- `gap_se::Float64`: Standard error on gap.
- `r_squared::Float64`: Goodness of fit (1 - RSS/TSS). NOT clamped; can be negative for bad fits.
- `converged::Bool`: Whether Levenberg-Marquardt optimization converged.
- `residuals::Vector{Float64}`: Fit residuals (observed - model).
- `times_used::Vector{Float64}`: Time points used in fit (after `skip_initial`).
- `values_used::Vector{Float64}`: Data values used in fit (after `skip_initial`).
"""
struct FitResult
    gap::Float64
    amplitude::Float64
    offset::Float64
    gap_ci::Tuple{Float64, Float64}
    gap_se::Float64
    r_squared::Float64
    converged::Bool
    residuals::Vector{Float64}
    times_used::Vector{Float64}
    values_used::Vector{Float64}
end

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

"""
    _log_linear_initial_guess(times, values) -> Vector{Float64}

Estimate initial parameters `[A, gap, C]` from data using log-linear regression.

Algorithm:
1. Estimate plateau C from the last 20% of data.
2. Subtract C, filter for positive shifted values.
3. Linear regression on log(shifted) vs time to extract A and gap.

Falls back to heuristic guess if fewer than 3 valid points remain after filtering.
"""
function _log_linear_initial_guess(times::AbstractVector{<:Real}, values::AbstractVector{<:Real})
    n = length(times)

    # Step 1: Estimate plateau C from last 20% of data
    tail_start = max(1, n - div(n, 5))
    C_guess = sum(values[tail_start:n]) / (n - tail_start + 1)

    # Step 2: Shift data and filter for positive values
    y_shifted = values .- C_guess
    mask = y_shifted .> 1e-10

    if sum(mask) < 3
        # Fallback: heuristic initial guess
        A_guess = values[1] - values[end]
        gap_guess = 1.0 / max(times[end] - times[1], eps(Float64))
        return [A_guess, gap_guess, C_guess]
    end

    t_valid = times[mask]
    log_y = log.(y_shifted[mask])

    # Step 3: Linear regression: log(y - C) = log(A) - gap * t
    X = hcat(ones(length(t_valid)), t_valid)
    coeffs = X \ log_y

    A_guess = exp(coeffs[1])
    gap_guess = max(-coeffs[2], 1e-6)  # ensure positive

    return [A_guess, gap_guess, C_guess]
end

"""
    _compute_r_squared(values, residuals) -> Float64

Compute R-squared = 1 - RSS/TSS.  NOT clamped to [0,1]; negative values
indicate the model fits worse than a horizontal line at the mean.
"""
function _compute_r_squared(values::AbstractVector{<:Real}, residuals::AbstractVector{<:Real})
    ss_res = sum(residuals .^ 2)
    y_mean = sum(values) / length(values)
    ss_tot = sum((values .- y_mean) .^ 2)
    return 1.0 - ss_res / ss_tot
end

# ---------------------------------------------------------------------------
# Main fitting function
# ---------------------------------------------------------------------------

"""
    fit_exponential_decay(times, values; skip_initial=0.0, p0=nothing, level=0.95) -> FitResult

Fit a single-exponential decay model `A * exp(-gap * t) + C` to time-series
data using Levenberg-Marquardt optimization (via LsqFit.jl).

The spectral gap `gap` is constrained to be non-negative via parameter bounds.
An automatic log-linear initial guess is generated unless `p0` is provided.

# Arguments
- `times::AbstractVector{<:Real}`: Time points (must match length of `values`).
- `values::AbstractVector{<:Real}`: Observed values at each time point.

# Keyword Arguments
- `skip_initial::Float64=0.0`: Fraction of initial data to skip (in [0, 1)).
  Use this to exclude early-time transients from the fit.
- `p0::Union{Nothing, Vector{Float64}}=nothing`: Initial parameter guess `[A, gap, C]`.
  If `nothing`, an automatic log-linear estimate is computed.
- `level::Float64=0.95`: Confidence level for the gap confidence interval.

# Returns
A [`FitResult`](@ref) containing the fitted gap, amplitude, offset, confidence
interval, standard error, R-squared, convergence status, and residuals.

# Example
```julia
times = collect(0.0:0.1:10.0)
values = 2.0 .* exp.(-0.5 .* times) .+ 0.3
result = fit_exponential_decay(times, values)
println("gap = \$(result.gap), R2 = \$(result.r_squared)")
```
"""
function fit_exponential_decay(
    times::AbstractVector{<:Real},
    values::AbstractVector{<:Real};
    skip_initial::Float64 = 0.0,
    p0::Union{Nothing, Vector{Float64}} = nothing,
    level::Float64 = 0.95,
)
    # --- Input validation ---
    length(times) == length(values) ||
        throw(ArgumentError("times and values must have the same length (got $(length(times)) and $(length(values)))"))
    length(times) >= 4 ||
        throw(ArgumentError("need at least 4 data points for 3-parameter fit (got $(length(times)))"))
    0.0 <= skip_initial < 1.0 ||
        throw(ArgumentError("skip_initial must be in [0, 1) (got $skip_initial)"))

    # --- Apply skip_initial BEFORE fitting (critical -- see RESEARCH pitfall 4) ---
    start_idx = max(1, floor(Int, skip_initial * length(times)) + 1)
    times_fit = Float64.(times[start_idx:end])
    values_fit = Float64.(values[start_idx:end])

    length(times_fit) >= 4 ||
        throw(ArgumentError("fewer than 4 data points remain after skip_initial=$skip_initial (got $(length(times_fit)))"))

    # --- Generate initial guess ---
    p0_used = if p0 === nothing
        _log_linear_initial_guess(times_fit, values_fit)
    else
        Float64.(p0)
    end

    # --- Set bounds: gap > 0 (FIT-05) ---
    lower = [-Inf, 0.0, -Inf]
    upper = [Inf, Inf, Inf]

    # --- Fit using LsqFit.jl ---
    # CRITICAL: do NOT pass a weight vector (see RESEARCH pitfall 6 about corrupted covariance)
    fit = curve_fit(_exp_decay_model, times_fit, values_fit, p0_used;
                    lower=lower, upper=upper)

    # --- Extract results using StatsAPI names ---
    params = coef(fit)                            # [A, gap, C]
    gap     = params[_IDX_GAP]
    A       = params[_IDX_A]
    C       = params[_IDX_C]

    # Standard errors and confidence intervals may fail with SingularException
    # when the Jacobian is rank-deficient (e.g., flat time series after skip_initial).
    # In that case, report Inf/NaN to signal unreliable uncertainty estimates.
    gap_se, gap_ci = try
        se = stderror(fit)
        ci = confint(fit; level=level)
        se[_IDX_GAP], (ci[_IDX_GAP][1], ci[_IDX_GAP][2])
    catch e
        e isa LinearAlgebra.SingularException || rethrow(e)
        Inf, (-Inf, Inf)
    end

    resid   = residuals(fit)                       # Vector{Float64}
    conv    = fit.converged                        # Bool

    # --- Compute R-squared (NOT provided by LsqFit) ---
    r2 = _compute_r_squared(values_fit, resid)

    return FitResult(gap, A, C, gap_ci, gap_se, r2, conv, resid, times_fit, values_fit)
end

# ============================================================================
# Bi-Exponential Decay Fitting: captures multi-timescale Liouvillian dynamics
# ============================================================================
#
# Model: y(t) = A1 * exp(-g1 * t) + A2 * exp(-g2 * t) + C
#
# Parameter ordering convention:
#   p[1] = A1       (fast amplitude)
#   p[2] = g1       (fast decay rate)
#   p[3] = A2       (slow amplitude)
#   p[4] = g2       (slow decay rate / spectral gap estimate)
#   p[5] = C        (offset / asymptotic value)

const _BIEXP_IDX_A1   = 1
const _BIEXP_IDX_G1   = 2
const _BIEXP_IDX_A2   = 3
const _BIEXP_IDX_G2   = 4
const _BIEXP_IDX_C    = 5

# Bi-exponential decay model for LsqFit (module-level)
_biexp_decay_model(t, p) = @. p[_BIEXP_IDX_A1] * exp(-p[_BIEXP_IDX_G1] * t) +
                               p[_BIEXP_IDX_A2] * exp(-p[_BIEXP_IDX_G2] * t) +
                               p[_BIEXP_IDX_C]

# ---------------------------------------------------------------------------
# BiexpFitResult struct
# ---------------------------------------------------------------------------

"""
    BiexpFitResult

Result of bi-exponential decay fitting: `A1 * exp(-g1 * t) + A2 * exp(-g2 * t) + C`.

After fitting, modes are sorted so that `gap_fast >= gap` (g1 >= g2).
The slow mode (`gap`) is the spectral gap estimate used for extrapolation.

# Fields
- `gap::Float64`: Slow decay rate (spectral gap estimate), used for extrapolation.
- `gap_fast::Float64`: Fast decay rate.
- `amplitude::Float64`: Slow-mode amplitude A2.
- `amplitude_fast::Float64`: Fast-mode amplitude A1.
- `offset::Float64`: Fitted offset C (asymptotic value).
- `gap_ci::Tuple{Float64, Float64}`: Confidence interval on slow gap (lower, upper).
- `gap_se::Float64`: Standard error on slow gap.
- `r_squared::Float64`: Goodness of fit (1 - RSS/TSS).
- `converged::Bool`: Whether Levenberg-Marquardt optimization converged.
- `residuals::Vector{Float64}`: Fit residuals (observed - model).
- `times_used::Vector{Float64}`: Time points used in fit (after `skip_initial`).
- `values_used::Vector{Float64}`: Data values used in fit (after `skip_initial`).
"""
struct BiexpFitResult
    gap::Float64
    gap_fast::Float64
    amplitude::Float64
    amplitude_fast::Float64
    offset::Float64
    gap_ci::Tuple{Float64, Float64}
    gap_se::Float64
    r_squared::Float64
    converged::Bool
    residuals::Vector{Float64}
    times_used::Vector{Float64}
    values_used::Vector{Float64}
end

# ---------------------------------------------------------------------------
# Internal helpers for bi-exponential fitting
# ---------------------------------------------------------------------------

"""
    _biexp_initial_guess(times, values, single_fit::FitResult) -> Vector{Float64}

Generate initial parameters `[A1, g1, A2, g2, C]` for bi-exponential fit
by analyzing residuals of the single-exponential fit.

Strategy:
1. Use the single-exp fit as the slow mode (A2, g2, C).
2. Fit the residuals to estimate the fast mode (A1, g1).
"""
function _biexp_initial_guess(
    times::AbstractVector{<:Real},
    values::AbstractVector{<:Real},
    single_fit::FitResult,
)
    # Slow mode from single-exp fit
    A2_guess = single_fit.amplitude
    g2_guess = single_fit.gap
    C_guess  = single_fit.offset

    # Residuals of single-exp fit
    resids = single_fit.residuals
    t_used = single_fit.times_used

    # Estimate fast mode from early residuals
    # If residuals show a fast-decaying pattern, use it
    n_early = max(3, div(length(resids), 4))
    early_resids = resids[1:n_early]
    A1_guess = maximum(abs.(early_resids))
    if A1_guess < 1e-12
        A1_guess = 0.1 * abs(A2_guess)
    end

    # Estimate fast gap from residual decay rate
    # Expect fast mode ~ 3-10x the slow gap
    g1_guess = max(3.0 * g2_guess, 1.0)

    return [A1_guess, g1_guess, A2_guess, g2_guess, C_guess]
end

"""
    fit_biexponential_decay(times, values; skip_initial=0.0, p0=nothing, level=0.95) -> BiexpFitResult

Fit a bi-exponential decay model `A1 * exp(-g1 * t) + A2 * exp(-g2 * t) + C`
to time-series data using Levenberg-Marquardt optimization (via LsqFit.jl).

The model captures multi-timescale dynamics in Liouvillian evolution.
After fitting, modes are sorted so `g1 >= g2` (fast >= slow).
The slow mode `g2` is the spectral gap estimate.

# Arguments
- `times::AbstractVector{<:Real}`: Time points.
- `values::AbstractVector{<:Real}`: Observed values.

# Keyword Arguments
- `skip_initial::Float64=0.0`: Fraction of initial data to skip (in [0, 1)).
- `p0::Union{Nothing, Vector{Float64}}=nothing`: Initial guess `[A1, g1, A2, g2, C]`.
  If `nothing`, automatic initial guess via single-exp fit + residual analysis.
- `level::Float64=0.95`: Confidence level for the gap confidence interval.

# Returns
A [`BiexpFitResult`](@ref) with sorted modes (fast/slow).
"""
function fit_biexponential_decay(
    times::AbstractVector{<:Real},
    values::AbstractVector{<:Real};
    skip_initial::Float64 = 0.0,
    p0::Union{Nothing, Vector{Float64}} = nothing,
    level::Float64 = 0.95,
)
    # --- Input validation ---
    length(times) == length(values) ||
        throw(ArgumentError("times and values must have the same length (got $(length(times)) and $(length(values)))"))
    length(times) >= 8 ||
        throw(ArgumentError("need at least 8 data points for 5-parameter bi-exponential fit (got $(length(times)))"))
    0.0 <= skip_initial < 1.0 ||
        throw(ArgumentError("skip_initial must be in [0, 1) (got $skip_initial)"))

    # --- Apply skip_initial ---
    start_idx = max(1, floor(Int, skip_initial * length(times)) + 1)
    times_fit = Float64.(times[start_idx:end])
    values_fit = Float64.(values[start_idx:end])

    length(times_fit) >= 8 ||
        throw(ArgumentError("fewer than 8 data points remain after skip_initial=$skip_initial (got $(length(times_fit)))"))

    # --- Generate initial guess ---
    p0_used = if p0 !== nothing
        Float64.(p0)
    else
        # Fit single-exp first, then use residuals to seed bi-exp
        single_fit = fit_exponential_decay(times_fit, values_fit; skip_initial=0.0, level=level)
        _biexp_initial_guess(times_fit, values_fit, single_fit)
    end

    # --- Set bounds: both gaps >= 0, amplitudes unconstrained, offset unconstrained ---
    lower = [-Inf, 0.0, -Inf, 0.0, -Inf]
    upper = [Inf, Inf, Inf, Inf, Inf]

    # --- Fit using LsqFit.jl ---
    fit = curve_fit(_biexp_decay_model, times_fit, values_fit, p0_used;
                    lower=lower, upper=upper)

    # --- Extract raw parameters ---
    params = coef(fit)
    A1_raw = params[_BIEXP_IDX_A1]
    g1_raw = params[_BIEXP_IDX_G1]
    A2_raw = params[_BIEXP_IDX_A2]
    g2_raw = params[_BIEXP_IDX_G2]
    C_raw  = params[_BIEXP_IDX_C]

    # --- Sort so g1 >= g2 (fast >= slow) ---
    # Track which raw index is the slow gap for SE/CI extraction
    if g1_raw >= g2_raw
        # Already sorted: g1=fast, g2=slow
        gap_fast = g1_raw
        gap_slow = g2_raw
        amp_fast = A1_raw
        amp_slow = A2_raw
        slow_gap_raw_idx = _BIEXP_IDX_G2
    else
        # Swap: g2 was actually faster
        gap_fast = g2_raw
        gap_slow = g1_raw
        amp_fast = A2_raw
        amp_slow = A1_raw
        slow_gap_raw_idx = _BIEXP_IDX_G1
    end

    # --- SE and CI for the slow gap ---
    gap_se, gap_ci = try
        se = stderror(fit)
        ci = confint(fit; level=level)
        se[slow_gap_raw_idx], (ci[slow_gap_raw_idx][1], ci[slow_gap_raw_idx][2])
    catch e
        e isa LinearAlgebra.SingularException || rethrow(e)
        Inf, (-Inf, Inf)
    end

    resid = residuals(fit)
    conv  = fit.converged

    # --- R-squared ---
    r2 = _compute_r_squared(values_fit, resid)

    return BiexpFitResult(
        gap_slow, gap_fast, amp_slow, amp_fast, C_raw,
        gap_ci, gap_se, r2, conv, resid, times_fit, values_fit,
    )
end
