# Phase 21: Exponential Fitting - Research

**Researched:** 2026-02-17
**Domain:** Nonlinear least-squares curve fitting / exponential decay analysis / LsqFit.jl
**Confidence:** HIGH

## Summary

Phase 21 implements single-exponential decay fitting (`A * exp(-gap * t) + C`) to extract spectral gap estimates from observable time-series data. The core dependency is LsqFit.jl v0.15, which provides Levenberg-Marquardt optimization with parameter bounds and confidence intervals via the StatsAPI interface. This phase operates entirely on synthetic data -- no trajectory dependency -- making it fully parallel with Phase 20.

The main technical challenges are: (1) automatic initial guess generation for the three-parameter model with offset, which cannot be linearized directly; (2) correct use of LsqFit.jl's `lower`/`upper` bounds to enforce gap > 0; and (3) computing R-squared manually since LsqFit.jl does not provide it natively. The log-linear initial guess approach works by estimating the plateau `C` from late-time data, subtracting it, then using linear regression on log-transformed residuals to get `A` and `gap`.

**Primary recommendation:** Create a single `fit_exponential_decay(times, values; skip_initial=0.0, ...)` function in a new `src/fitting.jl` file, returning a `FitResult` struct with gap, CI, R-squared, and standard error. Use LsqFit.jl's `curve_fit` with `lower=[0.0, 0.0, -Inf]` bounds to enforce gap > 0. Compute R-squared as `1 - RSS/TSS` manually.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| LsqFit.jl | 0.15 | Levenberg-Marquardt nonlinear least-squares fitting | Only pure-Julia NLLS library with built-in parameter bounds, confidence intervals, and StatsAPI integration. Standard choice in Julia ecosystem for curve fitting. |

### Supporting (already in project)
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| LinearAlgebra | stdlib | Matrix operations for log-linear initial guess | Always (already a dependency) |
| Statistics | stdlib | `mean()` for R-squared calculation and plateau estimation | Potentially useful but can be done with `sum()/length()` to avoid adding dependency |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| LsqFit.jl | LeastSquaresOptim.jl | Faster for large-scale problems, but no built-in confint/stderror; overkill for this use case |
| LsqFit.jl | Optim.jl (already in project) | Manual RSS minimization possible, but no Jacobian-based CIs or StatsAPI interface |
| Manual R-squared | StatsBase.jl | StatsBase already in test deps; but R-squared is trivial to compute inline |

**Installation:**
```julia
# In Julia REPL at project root:
using Pkg
Pkg.add("LsqFit")
```

Then add compat bound in Project.toml:
```toml
[deps]
LsqFit = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[compat]
LsqFit = "0.15"
```

**Note:** LsqFit.jl v0.15.1 requires Julia >= 1.6. QuantumFurnace.jl requires Julia >= 1.11. No compatibility issue.

## Architecture Patterns

### Recommended File Structure
```
src/
  fitting.jl           # NEW: FitResult struct + fit_exponential_decay + helpers
  convergence.jl       # Existing: observable builders (Phase 20)
  QuantumFurnace.jl    # Add: include("fitting.jl"), using LsqFit, exports

test/
  test_fitting.jl      # NEW: synthetic data tests for all FIT-* requirements
```

### Pattern 1: Result Struct with Quality Metrics
**What:** Return a dedicated struct from `fit_exponential_decay`, not a raw LsqFitResult.
**When to use:** Always. The downstream Phase 23 (`estimate_spectral_gap`) consumes this struct.
**Why:** LsqFitResult fields are implementation details. A clean `FitResult` struct provides a stable API boundary, allows adding R-squared (which LsqFit.jl does not provide), and follows the project's existing pattern (e.g., `ConvergenceData`, `TrajectoryResult`).

```julia
"""
    FitResult

Result of exponential decay fitting: A * exp(-gap * t) + C.

# Fields
- `gap::Float64`: Fitted decay rate (spectral gap estimate), constrained > 0.
- `amplitude::Float64`: Fitted amplitude A.
- `offset::Float64`: Fitted offset C (asymptotic value).
- `gap_ci::Tuple{Float64, Float64}`: 95% confidence interval on gap.
- `gap_se::Float64`: Standard error on gap.
- `r_squared::Float64`: Goodness of fit (1 - RSS/TSS).
- `converged::Bool`: Whether LM optimization converged.
- `residuals::Vector{Float64}`: Fit residuals.
- `times_used::Vector{Float64}`: Time points used in fit (after skip_initial).
- `values_used::Vector{Float64}`: Data values used in fit (after skip_initial).
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
```

### Pattern 2: Log-Linear Initial Guess (Three-Step Estimation)
**What:** Automatically estimate `[A, gap, C]` from data without user input.
**When to use:** Always as the default. The model `y = A * exp(-gap * t) + C` cannot be linearized directly because of the offset C.
**Algorithm:**

```
Step 1: Estimate C (plateau) from late-time data
  C_guess = mean(values[end-N:end])  where N = max(1, floor(0.2 * length))

Step 2: Subtract offset and take log
  y_shifted = values - C_guess
  # Filter: only use points where y_shifted > 0 (required for log)
  log_y = log.(y_shifted[mask])

Step 3: Linear regression on log_y vs times
  log_y = log(A) - gap * t
  Use least-squares: [log(A), -gap] = [ones(n) times] \ log_y
```

**Fallback:** If log-linear fails (e.g., non-monotone data, negative shifts), use heuristic defaults: `A = values[1] - values[end]`, `gap = 1.0 / (times[end] - times[1])`, `C = values[end]`.

### Pattern 3: Existing Module Conventions
**What:** Follow the project's established patterns for new modules.
**Details:**
- New file `src/fitting.jl` included in `src/QuantumFurnace.jl` after `convergence.jl`
- `using LsqFit` at the top of the main module (alongside other `using` statements)
- Export `fit_exponential_decay` and `FitResult` from the module
- Internal helper functions prefixed with `_` (e.g., `_log_linear_initial_guess`, `_compute_r_squared`)
- New test file `test/test_fitting.jl` included in `test/runtests.jl`

### Anti-Patterns to Avoid
- **Exposing LsqFitResult directly:** The downstream consumer should not need `using LsqFit` to interpret results. Wrap in `FitResult`.
- **Requiring manual initial guess:** The whole point of FIT-02 is auto-initialization. The `p0` keyword should exist for power users but never be required.
- **Using `estimate_covar` (deprecated):** Use `vcov(fit)` instead. The old name was deprecated in favor of StatsAPI conventions.
- **Using `confidence_interval` or `standard_error` (old API):** Use `confint(fit)` and `stderror(fit)` -- the StatsAPI names exported by LsqFit.jl v0.15.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Levenberg-Marquardt optimizer | Custom LM implementation | `LsqFit.curve_fit` | LM has many subtle edge cases (damping, step control, convergence criteria); LsqFit handles them |
| Confidence intervals | Manual Jacobian + covariance computation | `confint(fit; level=0.95)` | Requires correct handling of degrees of freedom, t-distribution quantiles, QR-based covariance |
| Standard errors | Manual sqrt(diag(vcov)) | `stderror(fit)` | Already scaled by degrees of freedom in LsqFit |
| Parameter bounds enforcement | Manual clamping in model | `lower=` / `upper=` kwargs to `curve_fit` | Bounds are integrated into LM step calculation, not just post-hoc clamping |

**Key insight:** R-squared is the ONE thing LsqFit.jl does NOT provide. It must be computed manually as `1 - sum(residuals.^2) / sum((values .- mean(values)).^2)`. This is trivial (3 lines) and not worth adding a dependency for.

## Common Pitfalls

### Pitfall 1: Log-Linear Guess Fails for Noisy Data Near Plateau
**What goes wrong:** When `values - C_guess` has negative or zero entries (noise pushes data below plateau), `log()` produces `NaN` or `-Inf`, corrupting the linear regression.
**Why it happens:** Real observable data has Monte Carlo noise. Near the plateau, `values[i] - C_guess` can be negative.
**How to avoid:** Filter to only use points where `y_shifted > threshold` (e.g., `threshold = 1e-10`). If fewer than 3 valid points remain, fall back to heuristic initial guess.
**Warning signs:** `NaN` in initial parameters, `curve_fit` fails to converge.

### Pitfall 2: Gap Bounds Must Match Parameter Vector Ordering
**What goes wrong:** If the model parameters are `[A, gap, C]` but bounds are `[gap_lower, A_lower, C_lower]`, the wrong parameter gets constrained.
**Why it happens:** LsqFit bounds are positional, mapping directly to `p0` indices. Easy to mix up.
**How to avoid:** Define a clear, documented parameter ordering convention. Use named constants: `const _IDX_A = 1; const _IDX_GAP = 2; const _IDX_C = 3`. Then `lower = [-Inf, 0.0, -Inf]` clearly maps gap (index 2) to lower bound 0.
**Warning signs:** Fit converges but returns negative gap, or amplitude is inexplicably bounded.

### Pitfall 3: confint Returns Tuples, Not Scalars
**What goes wrong:** `confint(fit)` returns `Vector{Tuple{Float64,Float64}}` -- one tuple per parameter. Code that expects a scalar CI width will fail.
**Why it happens:** The API returns the full interval `(lower, upper)` for each parameter at the requested confidence level.
**How to avoid:** Extract gap CI as `ci = confint(fit; level=0.95); gap_ci = ci[_IDX_GAP]` which gives `(lower_bound, upper_bound)`.
**Warning signs:** Type errors when trying to do arithmetic on confint output.

### Pitfall 4: skip_initial Must Be Applied Before Fitting, Not After
**What goes wrong:** If skip_initial truncates data after fitting, the fit uses the transient data and the window parameter has no effect.
**Why it happens:** Conceptual error in implementation order.
**How to avoid:** First compute `start_idx = max(1, floor(Int, skip_initial * length(times)) + 1)`, then pass `times[start_idx:end]` and `values[start_idx:end]` to `curve_fit`. The FitResult should store the truncated `times_used` and `values_used` for transparency.
**Warning signs:** Changing `skip_initial` has no effect on the fitted gap (Success Criterion 3 fails).

### Pitfall 5: R-Squared Can Be Negative for Bad Fits
**What goes wrong:** `R^2 = 1 - RSS/TSS` can be negative when the model fits worse than a horizontal line at the mean.
**Why it happens:** Bad initial guess or inappropriate model (data is not exponential decay).
**How to avoid:** Do NOT clamp R-squared to [0,1]. A negative R-squared is a useful diagnostic signal that the fit is bad. Report it honestly.
**Warning signs:** R-squared < 0 in test output. This is not a bug -- it means the fit failed.

### Pitfall 6: Passing Weight Vector of Ones Corrupts Covariance
**What goes wrong:** LsqFit.jl docs explicitly warn: "Passing vector of ones as the weight vector will cause mistakes in covariance estimation."
**Why it happens:** Internal normalization assumes non-trivial weights when weights are provided.
**How to avoid:** For unweighted fitting, do NOT pass a weights argument at all. Only use weights when you have actual variance estimates. For Phase 21, use the unweighted `curve_fit(model, times, values, p0; lower=lb, upper=ub)` signature.
**Warning signs:** Standard errors are wildly wrong (orders of magnitude off) when using uniform weights.

## Code Examples

Verified patterns from LsqFit.jl official documentation and source code:

### Basic Exponential Decay Fit with Bounds
```julia
# Source: LsqFit.jl README + API reference
using LsqFit

# Model: y = A * exp(-gap * t) + C
# Parameters: p = [A, gap, C]
model(t, p) = @. p[1] * exp(-p[2] * t) + p[3]

# Synthetic data
times = collect(0.0:0.1:10.0)
true_params = [2.0, 0.5, 0.3]  # A=2.0, gap=0.5, C=0.3
values = model(times, true_params) .+ 0.02 .* randn(length(times))

# Initial guess (from log-linear estimation)
p0 = [1.8, 0.4, 0.25]

# Enforce gap > 0 with parameter bounds
lb = [-Inf, 0.0, -Inf]   # lower bounds: only gap is bounded
ub = [Inf, Inf, Inf]      # upper bounds: none

fit = curve_fit(model, times, values, p0; lower=lb, upper=ub)
```

### Extracting All Quality Metrics from Fit Result
```julia
# Source: LsqFit.jl StatsAPI exports (v0.15)
using LsqFit

# Best-fit parameters
params = coef(fit)           # Vector{Float64} of [A, gap, C]
A_fit, gap_fit, C_fit = params

# Convergence check
converged = fit.converged     # Bool

# Standard errors (scaled by degrees of freedom)
se = stderror(fit)            # Vector{Float64}, one per parameter
gap_se = se[2]                # SE for gap parameter

# 95% confidence intervals
ci = confint(fit; level=0.95) # Vector{Tuple{Float64,Float64}}
gap_ci = ci[2]                # (lower, upper) for gap

# Residuals
resid = residuals(fit)        # same as fit.resid

# Covariance matrix (for advanced use)
cov_matrix = vcov(fit)        # Matrix{Float64}

# R-squared (NOT provided by LsqFit -- compute manually)
ss_res = sum(resid .^ 2)
y_mean = sum(values) / length(values)
ss_tot = sum((values .- y_mean) .^ 2)
r_squared = 1.0 - ss_res / ss_tot
```

### Log-Linear Initial Guess Algorithm
```julia
# Source: Standard numerical methods (not from LsqFit.jl)

function _log_linear_initial_guess(times::Vector{Float64}, values::Vector{Float64})
    n = length(times)

    # Step 1: Estimate plateau C from last 20% of data
    tail_start = max(1, n - div(n, 5))
    C_guess = sum(values[tail_start:n]) / (n - tail_start + 1)

    # Step 2: Shift data and filter for positive values
    y_shifted = values .- C_guess
    mask = y_shifted .> 1e-10  # avoid log(0) and log(negative)

    if sum(mask) < 3
        # Fallback: heuristic initial guess
        A_guess = values[1] - values[end]
        gap_guess = 1.0 / (times[end] - times[1])
        return [A_guess, gap_guess, C_guess]
    end

    t_valid = times[mask]
    log_y = log.(y_shifted[mask])

    # Step 3: Linear regression: log(y - C) = log(A) - gap * t
    # Build design matrix [ones, times]
    X = hcat(ones(length(t_valid)), t_valid)
    coeffs = X \ log_y  # least-squares solve

    A_guess = exp(coeffs[1])
    gap_guess = max(-coeffs[2], 1e-6)  # ensure positive

    return [A_guess, gap_guess, C_guess]
end
```

### Applying skip_initial Window Selection
```julia
# Source: Phase 21 requirement FIT-03

function _apply_skip_initial(times::Vector{Float64}, values::Vector{Float64},
                              skip_initial::Float64)
    @assert 0.0 <= skip_initial < 1.0 "skip_initial must be in [0, 1)"
    n = length(times)
    start_idx = max(1, floor(Int, skip_initial * n) + 1)
    return times[start_idx:end], values[start_idx:end]
end
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `estimate_covar(fit)` | `vcov(fit)` | LsqFit v0.15 | Old name deprecated; use StatsAPI convention |
| `standard_error(fit)` | `stderror(fit)` | LsqFit v0.15 | Renamed to match StatsAPI |
| `confidence_interval(fit, alpha)` | `confint(fit; level=0.95)` | LsqFit v0.15 | Note: old API used significance level alpha; new uses confidence level |
| Direct field access `fit.param` | `coef(fit)` | LsqFit v0.15 | StatsAPI accessor preferred (both work) |

**Deprecated/outdated:**
- `estimate_covar`: Deprecated in favor of `vcov`. Calling it emits a deprecation warning.
- `standard_error`: Renamed to `stderror` for StatsAPI compatibility.
- `confidence_interval`: Renamed to `confint`. Also note the parameter semantics changed from `alpha` (significance) to `level` (confidence).

## Open Questions

1. **Should FitResult include the raw LsqFitResult?**
   - What we know: Downstream Phase 23 only needs gap, CI, R-squared, and convergence status.
   - What's unclear: Whether advanced users might want access to the Jacobian or covariance matrix.
   - Recommendation: Do NOT include raw LsqFitResult in FitResult. Keep FitResult clean. If advanced access is needed later, it can be added in a future phase. Store only `residuals` as the most useful diagnostic.

2. **Should the model use `@.` broadcasting macro or explicit broadcasting?**
   - What we know: LsqFit.jl examples use `@.` macro: `model(t, p) = @. p[1] * exp(-p[2] * t) + p[3]`. The `t` argument can be a vector.
   - What's unclear: Whether `@.` has any performance implications vs explicit `.` operators.
   - Recommendation: Use `@.` macro for clarity, matching LsqFit.jl examples. Performance difference is negligible.

3. **Where in the module load order should fitting.jl be included?**
   - What we know: `convergence.jl` is included near the end of `QuantumFurnace.jl`. Fitting depends on LsqFit but not on any other QuantumFurnace modules.
   - Recommendation: Include `fitting.jl` right after `convergence.jl` and before `results.jl`, since Phase 23 will later connect fitting to convergence data. Add `using LsqFit` at the top of the module with other `using` statements.

## Sources

### Primary (HIGH confidence)
- LsqFit.jl Project.toml (v0.15.1): Julia >= 1.6, dependencies verified -- https://github.com/JuliaNLSolvers/LsqFit.jl/blob/master/Project.toml
- LsqFit.jl source code (curve_fit.jl): 6 method signatures, LsqFitResult struct verified -- https://github.com/JuliaNLSolvers/LsqFit.jl/blob/master/src/curve_fit.jl
- LsqFit.jl source code (LsqFit.jl): Exports verified (coef, confint, stderror, vcov, margin_error, dof, mse, nobs, rss, residuals, weights) -- https://github.com/JuliaNLSolvers/LsqFit.jl/blob/master/src/LsqFit.jl
- LsqFit.jl source code (levenberg_marquardt.jl): lower/upper bounds confirmed in function signature -- https://github.com/JuliaNLSolvers/LsqFit.jl/blob/master/src/levenberg_marquardt.jl
- LsqFit.jl README: curve_fit signature, bounds example, confint/stderror/vcov examples -- https://github.com/JuliaNLSolvers/LsqFit.jl/blob/master/README.md

### Secondary (MEDIUM confidence)
- LsqFit.jl Getting Started tutorial: Model definition, fit workflow -- https://julianlsolvers.github.io/LsqFit.jl/latest/getting_started/
- LsqFit.jl Tutorial: Weighted fitting, confidence intervals, standard errors, margin of error -- https://julianlsolvers.github.io/LsqFit.jl/latest/tutorial/
- LsqFit.jl maintenance discussion: Package is maintained but not highly active -- https://discourse.julialang.org/t/is-lsqfit-jl-abandoned/125192

### Tertiary (LOW confidence)
- Log-linear initial guess technique: Standard numerical methods, not verified from a single authoritative source but well-established in curve fitting literature.

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - LsqFit.jl source code and Project.toml directly verified. API confirmed from exports and function signatures.
- Architecture: HIGH - Follows established project patterns (ConvergenceData, TrajectoryResult structs). FitResult design driven by downstream Phase 23 requirements documented in ROADMAP.md.
- Pitfalls: HIGH - Weight-of-ones warning from official LsqFit.jl docs. Log-linear failure modes from numerical analysis fundamentals. Parameter ordering from positional kwargs design.
- Code examples: HIGH - Curve_fit, coef, confint, stderror, vcov signatures verified from LsqFit.jl source. R-squared formula is textbook statistics. Log-linear algorithm is standard.

**Research date:** 2026-02-17
**Valid until:** 2026-03-17 (LsqFit.jl is stable/slow-moving; 30-day validity appropriate)
