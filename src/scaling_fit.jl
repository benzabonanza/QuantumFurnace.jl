# ============================================================================
# Empirical (n, β) Scaling-Law Extraction for τ_mix  (qf-now)
# ============================================================================
#
# Two candidate forms, fitted in log τ-space via LsqFit.jl (LM-MLE):
#
#   M0 — separable power law:    log τ = c + x · log n + y · log β
#   M1 — power × Arrhenius:       log τ = c + x · log n + α · β
#
# Both reduce to linear regression in their respective parameters, so LM
# converges in a single iteration; using LsqFit gives uniform plumbing
# (stderror / confint / estimate_covar) with the rest of the package.
#
# Discrimination via AICc (Hurvich–Tsai small-sample correction); model
# weights via Burnham–Anderson Δ-AICc transform.
#
# Reference: drafts/scaling-analysis-research.md.

# --- Index conventions (parameter vectors are [c, x, slope] for both models) ---
const _SCALING_IDX_C     = 1
const _SCALING_IDX_X     = 2
const _SCALING_IDX_SLOPE = 3   # y for M0, α for M1

# --- Module-level model functions (LsqFit convention) ---
# `xdata` is an N×2 matrix; column 1 is log(n), column 2 is log(β) for M0 or β for M1.
_scaling_M0_model(xdata, p) = @. p[_SCALING_IDX_C] + p[_SCALING_IDX_X] * xdata[:, 1] + p[_SCALING_IDX_SLOPE] * xdata[:, 2]
_scaling_M1_model(xdata, p) = @. p[_SCALING_IDX_C] + p[_SCALING_IDX_X] * xdata[:, 1] + p[_SCALING_IDX_SLOPE] * xdata[:, 2]

# ---------------------------------------------------------------------------
# ScalingFit struct
# ---------------------------------------------------------------------------

"""
    ScalingFit

Result of an empirical (n, β) scaling-law fit produced by
[`fit_scaling`](@ref). One instance per candidate model. The fit is performed
in log τ-space.

# Model identification
- `model::Symbol`: `:M0` (separable power law) or `:M1` (power × Arrhenius).
- `param_names::NTuple{3, Symbol}`: `(:c, :x, :y)` for M0, `(:c, :x, :α)` for M1.

# Parameters (all indexed by `_SCALING_IDX_*`)
- `params::Vector{Float64}`: `[c, x, slope]` where slope = y (M0) or α (M1).
- `std_errors::Vector{Float64}`: Asymptotic standard errors from LsqFit
  covariance. `Inf` if the Jacobian was singular.
- `cis::Vector{Tuple{Float64, Float64}}`: Confidence intervals at `level`
  (default 95%). `(-Inf, Inf)` if SE estimation failed.
- `cov_matrix::Matrix{Float64}`: 3×3 parameter covariance.
- `corr_matrix::Matrix{Float64}`: 3×3 parameter correlation. `corr_matrix[2,3]`
  is the aliasing diagnostic — high absolute value means n-exponent and the
  β-side parameter are not independently identifiable (memo Pitfall 1).

# Model-comparison metrics
- `aicc::Float64`: AIC with Hurvich–Tsai small-sample correction.
- `log_likelihood::Float64`: Gaussian log-likelihood at MLE.
- `rss::Float64`: Residual sum of squares (in log τ).
- `sigma_residual::Float64`: MLE estimate of residual stdev.

# Provenance
- `n_data::Int`: Number of data points used in the fit.
- `converged::Bool`: Whether LM optimization converged.
- `n_values::Vector{Int}`: n values used (length `n_data`).
- `beta_values::Vector{Float64}`: β values used.
- `log_tau_observed::Vector{Float64}`: Observed log τ values.
- `log_tau_predicted::Vector{Float64}`: Model predictions at the data points.
- `residuals::Vector{Float64}`: `log_tau_observed - log_tau_predicted`.
"""
struct ScalingFit
    model::Symbol
    param_names::NTuple{3, Symbol}
    params::Vector{Float64}
    std_errors::Vector{Float64}
    cis::Vector{Tuple{Float64, Float64}}
    cov_matrix::Matrix{Float64}
    corr_matrix::Matrix{Float64}
    aicc::Float64
    log_likelihood::Float64
    rss::Float64
    sigma_residual::Float64
    n_data::Int
    converged::Bool
    n_values::Vector{Int}
    beta_values::Vector{Float64}
    log_tau_observed::Vector{Float64}
    log_tau_predicted::Vector{Float64}
    residuals::Vector{Float64}
end

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Compute AIC, AICc, and Gaussian log-likelihood for a fit with `n_data`
# points and `n_model_params` regression parameters. Treats σ² as an
# additional free parameter (Burnham–Anderson convention for NLS), so
# k = n_model_params + 1 in the AIC penalty.
function _scaling_aic_metrics(rss::Real, n_data::Integer, n_model_params::Integer)
    N = Int(n_data)
    k = n_model_params + 1                          # +1 for σ²
    if N <= k + 1
        # AICc denominator (N - k - 1) must be positive.
        # 3-param model + σ² ⇒ need N ≥ 6 for finite AICc.
        return (aicc = Inf, aic = Inf, log_likelihood = -Inf)
    end
    σ²_mle = rss / N
    log_L = -N / 2 * (log(2π) + log(σ²_mle) + 1.0)
    aic = 2k - 2 * log_L
    aicc = aic + 2k * (k + 1) / (N - k - 1)
    return (aicc = aicc, aic = aic, log_likelihood = log_L)
end

# Build a ScalingFit from a finished `curve_fit` LsqFitResult plus bookkeeping.
function _build_scaling_fit(
    model::Symbol,
    param_names::NTuple{3, Symbol},
    fit,
    n_vals::AbstractVector{<:Integer},
    beta_vals::AbstractVector{<:Real},
    log_τ::AbstractVector{<:Real},
    xdata::AbstractMatrix{<:Real},
    model_fn,
    level::Real,
)
    p = coef(fit)
    n_param = length(p)

    # SE / CI / covariance can fail for rank-deficient Jacobians (e.g., a
    # collinear (n, β) grid). LsqFit raises `SingularException` for an exactly
    # singular Jacobian and a plain `ErrorException` ("Covariance matrix is
    # negative …") when finite-difference numerics produce a non-PSD
    # covariance. Catch both and report Inf/NaN so the downstream consumer
    # can decide what to do.
    se, ci_vec, covmat = try
        se_v   = stderror(fit)
        ci_raw = confint(fit; level = level)            # already Vector{Tuple{Float64,Float64}}
        covm   = vcov(fit)                              # LsqFit ≥ 0.14; estimate_covar was deprecated
        se_v, Vector{Tuple{Float64, Float64}}(ci_raw), Matrix{Float64}(covm)
    catch e
        (e isa LinearAlgebra.SingularException || e isa ErrorException) || rethrow(e)
        fill(Inf, n_param),
        [(-Inf, Inf) for _ in 1:n_param],
        fill(NaN, n_param, n_param)
    end

    corrmat = if all(isfinite, covmat) && all(d -> d > 0, diag(covmat))
        s = sqrt.(diag(covmat))
        covmat ./ (s * s')
    else
        fill(NaN, n_param, n_param)
    end

    resid = residuals(fit)
    rss = sum(abs2, resid)
    metrics = _scaling_aic_metrics(rss, length(log_τ), n_param)
    σ_resid = length(log_τ) > 0 ? sqrt(rss / length(log_τ)) : NaN

    log_τ_pred = model_fn(xdata, p)

    return ScalingFit(
        model, param_names,
        Vector{Float64}(p), Vector{Float64}(se), ci_vec, covmat, corrmat,
        metrics.aicc, metrics.log_likelihood, rss, σ_resid,
        length(log_τ), fit.converged,
        Vector{Int}(n_vals), Vector{Float64}(beta_vals),
        Vector{Float64}(log_τ), Vector{Float64}(log_τ_pred),
        Vector{Float64}(resid),
    )
end

# ---------------------------------------------------------------------------
# Main fitting function
# ---------------------------------------------------------------------------

"""
    fit_scaling(n_vals, beta_vals, tau_vals; models=(:M0, :M1), level=0.95)
        -> Dict{Symbol, ScalingFit}

Fit one or more empirical scaling laws to a `(n, β, τ_mix)` sweep table.

The fit is performed in log τ-space. Both candidate models are linear in their
regression parameters, so LM converges in a single iteration; LsqFit is used
for its uniform plumbing with the rest of the package.

# Arguments
- `n_vals::AbstractVector{<:Integer}`: System-size values (qubit counts).
- `beta_vals::AbstractVector{<:Real}`: Inverse temperatures.
- `tau_vals::AbstractVector{<:Real}`: Mixing times τ_mix (must be positive).

# Keyword Arguments
- `models::NTuple{N, Symbol}=(:M0, :M1)`: Which models to fit. Currently
  `:M0` (separable power law `τ = C · n^x · β^y`) and `:M1`
  (power × Arrhenius `τ = C · n^x · exp(α·β)`).
- `level::Real=0.95`: Confidence level for parameter CIs.

# Returns
A `Dict{Symbol, ScalingFit}` keyed by model symbol. Use [`compare_models`](@ref)
to rank by AICc.

# Throws
- `ArgumentError` if the inputs have mismatched lengths, contain non-positive
  values where logarithm is required, or have fewer than 6 data points
  (3 regression params + σ² + 1 degree of freedom for AICc).

# Example
```julia
fits = fit_scaling([3,4,5,3,4,5,3,4,5], [5.,5.,5.,10.,10.,10.,20.,20.,20.],
                   [1.2, 2.4, 4.0, 3.5, 7.1, 13.0, 12.0, 25.0, 50.0])
println(formula_string(fits[:M0]))
println(compare_models(fits))
```
"""
function fit_scaling(
    n_vals::AbstractVector{<:Integer},
    beta_vals::AbstractVector{<:Real},
    tau_vals::AbstractVector{<:Real};
    models = (:M0, :M1),
    level::Real = 0.95,
)::Dict{Symbol, ScalingFit}

    # --- Input validation ---
    N = length(n_vals)
    length(beta_vals) == N || throw(ArgumentError(
        "n_vals and beta_vals must have the same length (got $N and $(length(beta_vals)))"))
    length(tau_vals) == N || throw(ArgumentError(
        "n_vals and tau_vals must have the same length (got $N and $(length(tau_vals)))"))
    N ≥ 6 || throw(ArgumentError(
        "need at least 6 data points for 3-parameter fit + σ² + AICc denominator (got $N)"))
    all(n -> n > 0, n_vals)    || throw(ArgumentError("all n must be positive"))
    all(b -> b > 0, beta_vals) || throw(ArgumentError("all β must be positive"))
    all(t -> t > 0 && isfinite(t), tau_vals) || throw(ArgumentError(
        "all τ_mix must be positive and finite (cannot take log of zero, negative, or non-finite)"))
    isempty(models) && throw(ArgumentError("models must be non-empty"))
    for m in models
        m in (:M0, :M1) || throw(ArgumentError("unknown model :$m (supported: :M0, :M1)"))
    end
    0 < level < 1 || throw(ArgumentError("level must be in (0, 1) (got $level)"))

    log_n = log.(Float64.(n_vals))
    log_β = log.(Float64.(beta_vals))
    log_τ = log.(Float64.(tau_vals))

    out = Dict{Symbol, ScalingFit}()

    if :M0 in models
        xdata = hcat(log_n, log_β)
        # Initial guess: c=0, x=1 (mild superlinear in n), y=1 (linear in β).
        p0 = [0.0, 1.0, 1.0]
        fit = curve_fit(_scaling_M0_model, xdata, log_τ, p0)
        out[:M0] = _build_scaling_fit(:M0, (:c, :x, :y), fit,
                                       n_vals, beta_vals, log_τ, xdata,
                                       _scaling_M0_model, level)
    end

    if :M1 in models
        xdata = hcat(log_n, Float64.(beta_vals))
        # Initial guess: c=0, x=1, α=0.1 (mild Arrhenius slope).
        p0 = [0.0, 1.0, 0.1]
        fit = curve_fit(_scaling_M1_model, xdata, log_τ, p0)
        out[:M1] = _build_scaling_fit(:M1, (:c, :x, :α), fit,
                                       n_vals, beta_vals, log_τ, xdata,
                                       _scaling_M1_model, level)
    end

    return out
end

"""
    fit_scaling(table::NamedTuple; kwargs...) -> Dict{Symbol, ScalingFit}

NamedTuple convenience wrapper. Expects fields `:n`, `:beta`, and `:tau_mix`.
"""
function fit_scaling(table::NamedTuple; kwargs...)
    haskey(table, :n) && haskey(table, :beta) && haskey(table, :tau_mix) || throw(ArgumentError(
        "table NamedTuple must have fields :n, :beta, :tau_mix (got $(keys(table)))"))
    return fit_scaling(table.n, table.beta, table.tau_mix; kwargs...)
end

"""
    fit_scaling(results::Vector{<:NamedTuple}; source_filter=(:extrapolated,), kwargs...)
        -> Dict{Symbol, ScalingFit}

Vector-of-NamedTuple convenience wrapper for the output of
[`sweep_mixing_times`](@ref) and [`sweep_channel_mixing`](@ref). Entries with
`mixing_time_source ∉ source_filter` or non-finite `mixing_time` are filtered
out before fitting.

`source_filter=(:extrapolated,)` (the default) keeps only the rigorously
extrapolated cells from the eigenmode τ_mix schema (qf-e4y); pass
`(:extrapolated, :floor)` to also include the gap-based bound, but be aware
that `:floor` cells have lower-bound semantics (the true τ_mix may be
infinite) and will bias the fit.
"""
function fit_scaling(
    results::Vector{<:NamedTuple};
    source_filter = (:extrapolated,),
    kwargs...,
)
    # Channel sweeps use :tau_mix / :tau_mix_source; Lindbladian sweeps use
    # :mixing_time / :mixing_time_source. Support both.
    function _get_tau(r)
        if haskey(r, :mixing_time)
            return r.mixing_time
        elseif haskey(r, :tau_mix)
            return r.tau_mix
        end
        return nothing
    end
    function _get_source(r)
        if haskey(r, :mixing_time_source)
            return r.mixing_time_source
        elseif haskey(r, :tau_mix_source)
            return r.tau_mix_source
        end
        return nothing
    end

    valid = NamedTuple[]
    for r in results
        haskey(r, :n) && haskey(r, :beta) || continue
        τ = _get_tau(r)
        τ === nothing && continue
        (τ isa Real && isfinite(τ) && τ > 0) || continue
        src = _get_source(r)
        # Strict policy: if the caller specified a source_filter, require the
        # source field to be present and match. Pass `source_filter = ()` to
        # disable source filtering entirely.
        if !isempty(source_filter)
            (src !== nothing && src in source_filter) || continue
        end
        push!(valid, r)
    end
    isempty(valid) && throw(ArgumentError(
        "no valid entries found (need :n, :beta, finite positive :mixing_time or :tau_mix" *
        " with :mixing_time_source/:tau_mix_source ∈ $source_filter)"))

    n_vals = [Int(r.n) for r in valid]
    β_vals = [Float64(r.beta) for r in valid]
    τ_vals = Float64[_get_tau(r) for r in valid]
    return fit_scaling(n_vals, β_vals, τ_vals; kwargs...)
end

# ---------------------------------------------------------------------------
# Predict on a single point
# ---------------------------------------------------------------------------

"""
    predict_scaling(fit::ScalingFit, n::Real, β::Real) -> Float64

Predicted τ_mix at `(n, β)` from the fitted model. Returns the linear-scale
prediction `exp(log τ̂)`.
"""
function predict_scaling(fit::ScalingFit, n::Real, β::Real)
    c     = fit.params[_SCALING_IDX_C]
    x     = fit.params[_SCALING_IDX_X]
    slope = fit.params[_SCALING_IDX_SLOPE]
    log_n = log(Float64(n))
    log_τ_pred = if fit.model === :M0
        c + x * log_n + slope * log(Float64(β))
    elseif fit.model === :M1
        c + x * log_n + slope * Float64(β)
    else
        throw(ArgumentError("unknown model :$(fit.model) in predict_scaling"))
    end
    return exp(log_τ_pred)
end

# ---------------------------------------------------------------------------
# Model comparison
# ---------------------------------------------------------------------------

"""
    aicc_weights(fits::Dict{Symbol, ScalingFit}) -> Dict{Symbol, Float64}

Compute AICc model weights (Burnham–Anderson). For a set of candidate models
with AICc values `{AICc_i}`,
``w_i = exp(-Δ_i / 2) / Σ_j exp(-Δ_j / 2)``
where ``Δ_i = AICc_i - AICc_{\\min}``.

Higher weight = stronger support. Rules of thumb (Burnham–Anderson 2002):
- `Δ ≤ 2`: substantial support for both / not distinguishable
- `2 < Δ ≤ 7`: considerably less support
- `Δ > 10`: essentially no support
"""
function aicc_weights(fits::Dict{Symbol, ScalingFit})
    isempty(fits) && return Dict{Symbol, Float64}()
    aicc_min = minimum(f.aicc for f in values(fits))
    raw = Dict(k => exp(-(f.aicc - aicc_min) / 2) for (k, f) in fits)
    Z = sum(values(raw))
    return Dict(k => v / Z for (k, v) in raw)
end

"""
    compare_models(fits::Dict{Symbol, ScalingFit}) -> NamedTuple

Rank fitted models by AICc and return a summary table.

# Returns
NamedTuple with fields:
- `ranked::Vector{Symbol}`: model symbols, best AICc first.
- `aicc::Vector{Float64}`: AICc values aligned with `ranked`.
- `delta_aicc::Vector{Float64}`: AICc differences from the best model.
- `weights::Vector{Float64}`: AICc model weights aligned with `ranked`.
"""
function compare_models(fits::Dict{Symbol, ScalingFit})
    isempty(fits) && throw(ArgumentError("fits dictionary is empty"))
    ordered = sort(collect(fits); by = kv -> kv[2].aicc)
    ranked = [kv[1] for kv in ordered]
    aiccs  = [kv[2].aicc for kv in ordered]
    aicc_min = aiccs[1]
    delta  = aiccs .- aicc_min
    w_dict = aicc_weights(fits)
    weights = [w_dict[k] for k in ranked]
    return (
        ranked     = ranked,
        aicc       = aiccs,
        delta_aicc = delta,
        weights    = weights,
    )
end

# ---------------------------------------------------------------------------
# Human-readable formula
# ---------------------------------------------------------------------------

"""
    formula_string(fit::ScalingFit) -> String

Produce a human-readable formula with ±1σ uncertainties propagated through
the `exp(c) → C` transform.

# Examples
- M0: `"τ_mix = (0.13 ± 0.02) · n^(2.30 ± 0.30) · β^(1.70 ± 0.20)"`
- M1: `"τ_mix = (0.13 ± 0.02) · n^(2.30 ± 0.30) · exp((0.1500 ± 0.0200)·β)"`
"""
function formula_string(fit::ScalingFit)
    c     = fit.params[_SCALING_IDX_C]
    x     = fit.params[_SCALING_IDX_X]
    slope = fit.params[_SCALING_IDX_SLOPE]
    σc    = fit.std_errors[_SCALING_IDX_C]
    σx    = fit.std_errors[_SCALING_IDX_X]
    σs    = fit.std_errors[_SCALING_IDX_SLOPE]

    # Propagate σ(c) → σ(C) through C = exp(c): σ_C ≈ C · σ_c (first-order).
    C = exp(c)
    σC = isfinite(σc) ? C * σc : Inf

    if fit.model === :M0
        return @sprintf("τ_mix = (%.3g ± %.2g) · n^(%.2f ± %.2f) · β^(%.2f ± %.2f)",
                        C, σC, x, σx, slope, σs)
    elseif fit.model === :M1
        return @sprintf("τ_mix = (%.3g ± %.2g) · n^(%.2f ± %.2f) · exp((%.4f ± %.4f)·β)",
                        C, σC, x, σx, slope, σs)
    else
        return "(unknown model :$(fit.model))"
    end
end

# ---------------------------------------------------------------------------
# Diagnostic-data helper (plot scripts consume this; src/ stays Plots-free)
# ---------------------------------------------------------------------------

"""
    scaling_fit_grid(fit::ScalingFit; n_grid=nothing, beta_grid=nothing)
        -> NamedTuple

Build a regular (n, β) grid covering the data range and evaluate the fitted
model on it. Returns a NamedTuple with the grids, the τ_mix prediction matrix,
plus a sparse `(observed, predicted, residual)` triple aligned with the
original data — everything a plotting script needs for a 2-panel data vs fit
diagnostic.

# Keyword Arguments
- `n_grid::Union{Nothing, AbstractVector}=nothing`: explicit n grid; defaults
  to `minimum(n_values):maximum(n_values)`.
- `beta_grid::Union{Nothing, AbstractVector}=nothing`: explicit β grid;
  defaults to a 21-point log-spaced grid spanning the β range.

# Returns
NamedTuple with fields:
- `n_grid::Vector{Int}`
- `beta_grid::Vector{Float64}`
- `tau_predicted::Matrix{Float64}`  — size `(length(n_grid), length(beta_grid))`
- `n_obs::Vector{Int}` / `beta_obs::Vector{Float64}` — original data coords
- `tau_obs::Vector{Float64}` — observed τ_mix at the data points (linear scale)
- `tau_pred_at_obs::Vector{Float64}` — model prediction at each data point
- `residuals_log::Vector{Float64}` — log τ residuals (positive ⇒ data > model)
"""
function scaling_fit_grid(
    fit::ScalingFit;
    n_grid::Union{Nothing, AbstractVector} = nothing,
    beta_grid::Union{Nothing, AbstractVector} = nothing,
)
    n_lo, n_hi = extrema(fit.n_values)
    β_lo, β_hi = extrema(fit.beta_values)

    n_g = n_grid === nothing ? collect(n_lo:n_hi) : collect(Int, n_grid)
    β_g = if beta_grid === nothing
        # Log-spaced 21-point grid; if the data range is a single β, fall
        # back to a single point so the matrix below is still well-defined.
        if β_hi ≈ β_lo
            [Float64(β_lo)]
        else
            exp.(range(log(β_lo), log(β_hi); length = 21))
        end
    else
        collect(Float64, beta_grid)
    end

    τ_pred_mat = Matrix{Float64}(undef, length(n_g), length(β_g))
    @inbounds for (i, nv) in enumerate(n_g), (j, βv) in enumerate(β_g)
        τ_pred_mat[i, j] = predict_scaling(fit, nv, βv)
    end

    # Per-datapoint predictions and residuals (already on the fit).
    tau_obs = exp.(fit.log_tau_observed)
    tau_pred_at_obs = exp.(fit.log_tau_predicted)

    return (
        n_grid          = n_g,
        beta_grid       = β_g,
        tau_predicted   = τ_pred_mat,
        n_obs           = fit.n_values,
        beta_obs        = fit.beta_values,
        tau_obs         = tau_obs,
        tau_pred_at_obs = tau_pred_at_obs,
        residuals_log   = fit.residuals,
    )
end
