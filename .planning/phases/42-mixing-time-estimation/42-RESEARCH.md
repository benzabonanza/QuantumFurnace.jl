# Phase 42: Mixing Time Estimation - Research

**Researched:** 2026-03-01
**Domain:** Exponential curve fitting, spectral gap extraction, LsqFit.jl integration, Julia package management
**Confidence:** HIGH

## Summary

Phase 42 implements mixing time estimation as a post-processing function on `ThermalizeResults`. The core operation is fitting `d(t) = A*exp(-gap*t) + C` to the trace distance convergence curve produced by `run_thermalize`, extracting the effective spectral gap, and optionally extrapolating to a target epsilon. The existing `src/staging/fitting.jl` already contains a complete, tested `fit_exponential_decay` function using LsqFit.jl with Levenberg-Marquardt optimization, parameter bounds, confidence intervals, and R-squared computation. The `FitResult` struct and `_log_linear_initial_guess` helper are ready for promotion.

The phase requires three main work streams: (1) re-adding LsqFit.jl as an active dependency and promoting fitting.jl from staging to active source, (2) implementing `estimate_mixing_time` as a post-processing wrapper around `fit_exponential_decay` with a new `MixingTimeEstimate` struct, and (3) adding quality gates (R-squared thresholds, offset warnings) and the extrapolation formula. The DM trace distance curve is deterministic (no trajectory noise), making single-exponential fitting significantly more reliable than for trajectory-based observables.

**Primary recommendation:** Promote the existing `staging/fitting.jl` wholesale, add LsqFit.jl v0.15.x to Project.toml, then build `estimate_mixing_time(result::ThermalizeResults; ...)` as a thin wrapper that pre-processes the trace distance data, delegates to `fit_exponential_decay`, and post-processes the fit into a `MixingTimeEstimate`. The `extrapolate=true` mode computes the extrapolated mixing time from the fit result but does NOT control simulation execution (per the prior decision that mixing time estimation is post-processing, not embedded in `run_thermalize`).

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| LsqFit.jl | 0.15.x (latest: 0.15.1) | Levenberg-Marquardt nonlinear least-squares fitting | Only Julia package for bounded curve_fit with StatsAPI-compatible coef/stderror/confint/residuals. Already used in staging/fitting.jl. |
| LinearAlgebra (stdlib) | Julia 1.11+ | Matrix operations for trace distance | Already a dependency |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| Test (stdlib) | Julia 1.11+ | Unit tests | Test mixing time estimation |
| StableRNGs.jl | 1.x | Deterministic test RNG | Already in test extras |
| Printf (stdlib) | Julia 1.11+ | Warning messages | Already a dependency |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| LsqFit.jl | Optim.jl + manual Jacobian | LsqFit wraps the same algorithm with less boilerplate; provides confint/stderror out of the box |
| LsqFit.jl | Hand-rolled Gauss-Newton | Missing parameter bounds, covariance estimation, confidence intervals -- LsqFit handles all edge cases |

**Installation:**
```julia
# In Julia REPL from project root:
using Pkg
Pkg.add("LsqFit")
# Or manually add to [deps] in Project.toml:
# LsqFit = "2fda8390-95c7-5789-9bda-21331edee243"
# And to [compat]:
# LsqFit = "0.15"
```

## Architecture Patterns

### Recommended Project Structure
```
src/
  fitting.jl              # PROMOTED from src/staging/fitting.jl (FitResult, fit_exponential_decay, helpers)
  mixing.jl               # NEW: MixingTimeEstimate struct, estimate_mixing_time function, quality gates
  staging/
    fitting.jl            # REMOVED (promoted to src/fitting.jl)
    gap_estimation.jl     # STAYS in staging (trajectory-based, uses FitResult from promoted fitting.jl)
    log_sobolev.jl        # STAYS in staging (commented out, unrelated)
test/
  test_mixing.jl          # NEW: tests for estimate_mixing_time
  staging/
    test_fitting.jl       # STAYS (tests for fit_exponential_decay, now tests promoted code)
    test_gap_estimation.jl # STAYS
```

### Pattern 1: Post-Processing Wrapper
**What:** `estimate_mixing_time` takes a completed `ThermalizeResults` and returns a `MixingTimeEstimate`. It does NOT modify or control `run_thermalize`.
**When to use:** Always -- this is the core pattern for Phase 42.
**Example:**
```julia
# Source: Codebase pattern from staging/fitting.jl + FEATURES.md specification
function estimate_mixing_time(
    result::ThermalizeResults;
    skip_initial::Float64 = 0.2,
    target_epsilon::Union{Nothing, Float64} = nothing,
    extrapolate::Bool = false,
    level::Float64 = 0.95,
)
    times = result.time_steps
    dists = result.trace_distances

    # Validate sufficient data
    length(times) >= 10 || throw(ArgumentError(
        "Need at least 10 data points for mixing time estimation (got $(length(times)))"))

    # Delegate to existing fit_exponential_decay
    fit = fit_exponential_decay(Float64.(times), Float64.(dists);
        skip_initial=skip_initial, level=level)

    # Quality gate warnings
    _check_fit_quality(fit, target_epsilon)

    # Compute mixing times
    t_mix_actual = _find_actual_mixing_time(times, dists, target_epsilon)
    t_mix_extrap = extrapolate ? _extrapolate_mixing_time(fit, target_epsilon) : nothing

    return MixingTimeEstimate(
        fitted_gap = fit.gap,
        amplitude = fit.amplitude,
        offset = fit.offset,
        gap_ci = fit.gap_ci,
        gap_se = fit.gap_se,
        r_squared = fit.r_squared,
        converged = fit.converged,
        mixing_time = extrapolate ? t_mix_extrap : t_mix_actual,
        mixing_time_extrapolated = t_mix_extrap,
        mixing_time_actual = t_mix_actual,
        target_epsilon = target_epsilon,
        fit_result = fit,
    )
end
```

### Pattern 2: MixingTimeEstimate Struct (Result Wrapper)
**What:** A struct that wraps FitResult with mixing-time-specific fields.
**When to use:** Return type of `estimate_mixing_time`.
**Example:**
```julia
# Source: FEATURES.md specification + codebase pattern
struct MixingTimeEstimate
    # Fit parameters (delegated from FitResult)
    fitted_gap::Float64
    amplitude::Float64
    offset::Float64
    gap_ci::Tuple{Float64, Float64}
    gap_se::Float64
    r_squared::Float64
    converged::Bool
    # Mixing time results
    mixing_time::Float64                          # primary answer (actual or extrapolated)
    mixing_time_extrapolated::Union{Nothing, Float64}
    mixing_time_actual::Union{Nothing, Float64}   # from data (NaN/nothing if target not reached)
    target_epsilon::Union{Nothing, Float64}
    # Full fit result for advanced users
    fit_result::FitResult
end
```

### Pattern 3: Quality Gate Warnings via @warn
**What:** Use Julia's `@warn` macro for quality gate warnings (not errors).
**When to use:** When R-squared < 0.95, offset C is large, CI is wide.
**Example:**
```julia
function _check_fit_quality(fit::FitResult, target_epsilon::Union{Nothing, Float64})
    if fit.r_squared < 0.95
        @warn "Fit R-squared = $(fit.r_squared) < 0.95. Single-exponential model may not describe the data well. Consider increasing skip_initial or running longer."
    end
    if target_epsilon !== nothing && fit.offset > 0.1 * target_epsilon
        @warn "Fit offset C = $(fit.offset) is large relative to target epsilon = $(target_epsilon). Extrapolation may be unreliable."
    end
    if !fit.converged
        @warn "Levenberg-Marquardt optimization did not converge. Fit results are unreliable."
    end
    if fit.gap_se > 0.5 * fit.gap
        @warn "Gap standard error ($(fit.gap_se)) exceeds 50% of fitted gap ($(fit.gap)). Confidence interval is very wide."
    end
end
```

### Pattern 4: Extrapolation Formula
**What:** Compute mixing time to target epsilon from fitted parameters.
**When to use:** When `extrapolate=true` and `target_epsilon` is provided.
**Example:**
```julia
function _extrapolate_mixing_time(fit::FitResult, target_epsilon::Union{Nothing, Float64})
    target_epsilon === nothing && return nothing
    fit.gap <= 0.0 && return nothing
    fit.amplitude <= 0.0 && return nothing
    # d(t) = A * exp(-gap * t) + C
    # Solve A * exp(-gap * t_mix) + C = epsilon
    # => t_mix = -ln((epsilon - C) / A) / gap
    effective_target = target_epsilon - fit.offset
    effective_target <= 0.0 && return nothing  # offset already exceeds target
    return -log(effective_target / fit.amplitude) / fit.gap
end
```

### Anti-Patterns to Avoid
- **Embedding fitting in run_thermalize:** The prior decision explicitly states mixing time estimation is post-processing. Do NOT add fitting logic inside the DM evolution loop.
- **Using trajectory-based gap estimation for DM data:** The DM trace distance is a direct, noiseless signal. Do NOT route it through `estimate_spectral_gap` (which is designed for noisy trajectory-averaged observables).
- **Two-exponential fitting:** FEATURES.md explicitly lists this as an anti-feature for DM trace distance. Single-exp + skip_initial is sufficient for the noiseless DM curve.
- **Bootstrap uncertainty on DM data:** DM evolution is deterministic. Bootstrap is for trajectory sampling noise. Use LsqFit Jacobian-based confidence intervals.
- **Mutable MixingTimeEstimate struct:** Keep it immutable (plain struct, not mutable struct) following the FitResult pattern.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Nonlinear least-squares fitting | Manual Gauss-Newton or gradient descent | LsqFit.jl `curve_fit` | Handles parameter bounds, Jacobian estimation, convergence detection, covariance matrix |
| Confidence intervals on fit params | Manual Jacobian + t-distribution | LsqFit.jl `confint(fit; level=0.95)` | Correct degrees-of-freedom handling, handles singular Jacobian gracefully |
| Standard errors on fit params | Manual covariance extraction | LsqFit.jl `stderror(fit)` | Returns sqrt of diagonal of covariance matrix, properly scaled by residual variance |
| R-squared computation | Build from scratch | Existing `_compute_r_squared` in fitting.jl | Already handles edge cases (returns un-clamped R-squared, negative for bad fits) |
| Initial parameter guess | Random guessing | Existing `_log_linear_initial_guess` | Log-linear regression gives O(1) initial guess, avoids LM divergence |

**Key insight:** The `staging/fitting.jl` code was carefully designed to handle all edge cases (SingularException from stderror, fallback initial guess, un-clamped R-squared). Reuse it wholesale rather than rebuilding.

## Common Pitfalls

### Pitfall 1: Forgetting to add LsqFit.jl to ALL three Project.toml sections
**What goes wrong:** Adding to `[deps]` but forgetting `[compat]`, or adding to `[deps]` but the staging tests need it in `[extras]` too.
**Why it happens:** Julia's Pkg system has three sections: `[deps]` for runtime, `[compat]` for version bounds, `[extras]` for test-only deps.
**How to avoid:** LsqFit.jl is a runtime dependency (used in active source), so add to `[deps]` and `[compat]`. Since the tests call `fit_exponential_decay` through the module, it does not need to be in `[extras]`.
**Warning signs:** `Pkg.precompile()` fails or tests fail with `UndefVarError: curve_fit not defined`.

### Pitfall 2: Including fitting.jl in QuantumFurnace.jl without `using LsqFit`
**What goes wrong:** The fitting.jl code calls `curve_fit`, `coef`, `stderror`, `confint`, `residuals` which are provided by LsqFit.jl. Without `using LsqFit` in the module file, these are undefined.
**Why it happens:** The staging code was never included in the active module, so the `using LsqFit` was never added.
**How to avoid:** Add `using LsqFit` to `src/QuantumFurnace.jl` before the `include("fitting.jl")` line.
**Warning signs:** `UndefVarError: curve_fit not defined` at module load time.

### Pitfall 3: skip_initial applied TWICE (in estimate_mixing_time AND fit_exponential_decay)
**What goes wrong:** `estimate_mixing_time` pre-processes the data with skip_initial, then passes it to `fit_exponential_decay` which also applies skip_initial.
**Why it happens:** Both functions accept skip_initial. If the wrapper applies it and passes the keyword through, the data gets truncated twice.
**How to avoid:** Either (a) have `estimate_mixing_time` pass `skip_initial=0.0` to `fit_exponential_decay` after doing its own truncation, or (b) pass skip_initial through to `fit_exponential_decay` and do NO pre-processing in the wrapper.
**Warning signs:** Fit uses far fewer data points than expected; fitted gap is wrong for short runs.
**Recommendation:** Option (b) -- let `fit_exponential_decay` handle skip_initial (it already does it correctly). Do NOT duplicate the logic.

### Pitfall 4: Extrapolation with negative effective_target
**What goes wrong:** If `fit.offset >= target_epsilon`, then `effective_target = target_epsilon - fit.offset <= 0`, and `log(effective_target / A)` is undefined.
**Why it happens:** The offset C represents the asymptotic floor of the trace distance. If C >= epsilon, the system physically cannot reach epsilon (the fixed point accuracy limits it).
**How to avoid:** Check `effective_target > 0` before computing log. Return nothing and issue a warning.
**Warning signs:** `DomainError` on log of negative number.

### Pitfall 5: ThermalizeResults struct is immutable -- cannot extend fields
**What goes wrong:** MIX-05 says "ThermalizeResults extended with mixing_time, fitted_gap, and fit quality metrics." But ThermalizeResults is an immutable struct and cannot be modified without breaking serialization.
**Why it happens:** The struct was designed in Phase 36 with fixed fields. Adding fields breaks BSON serialization compatibility.
**How to avoid:** Do NOT modify ThermalizeResults. Instead, `MixingTimeEstimate` is a separate struct that wraps/references the ThermalizeResults data. MIX-05 should be interpreted as "the mixing time estimation output includes these fields" not "modify the existing struct." The metadata Dict can optionally store mixing_time info if the user wants to persist it.
**Warning signs:** BSON load_result fails on old files.

### Pitfall 6: Not exporting new public types and functions
**What goes wrong:** User cannot access `estimate_mixing_time`, `MixingTimeEstimate`, `FitResult`, `fit_exponential_decay` from `using QuantumFurnace`.
**Why it happens:** Promoting from staging means these were previously only accessible via `QuantumFurnace.fit_exponential_decay` (internal path). They need explicit export.
**How to avoid:** Add all new public names to the export list in `src/QuantumFurnace.jl`.
**Warning signs:** Tests work (they `using QuantumFurnace` internally) but user scripts fail.

### Pitfall 7: Staging tests use internal fixtures not available from active test suite
**What goes wrong:** `test/staging/test_fitting.jl` and `test/staging/test_gap_estimation.jl` reference `N3_*` fixtures and `make_config` from `test_helpers.jl`. When promoted, these tests must be run within the main test suite context.
**Why it happens:** Staging tests were designed to run in isolation or from the main runtests.jl which includes test_helpers.jl first.
**How to avoid:** Include the promoted test file in `runtests.jl` after `test_helpers.jl` is loaded, same as all other test files.
**Warning signs:** `UndefVarError: N3_HAM not defined` when running staging tests.

### Pitfall 8: Confusing "mixing time to convergence" with "mixing time to epsilon"
**What goes wrong:** MIX-03 says "actual number of steps/time to reach target trace distance." If no target_epsilon is provided, what does "actual mixing time" mean?
**Why it happens:** The convergence_cutoff in run_thermalize is hardcoded to 1e-5. The user may want a different target.
**How to avoid:** When `target_epsilon` is provided, find the first time_step where trace_distance <= target_epsilon. When not provided, use the last data point or the convergence_cutoff from run_thermalize metadata. Document clearly what "mixing time" means in each case.
**Warning signs:** User confusion about what the returned mixing_time means.

## Code Examples

Verified patterns from the existing codebase and official sources:

### Promoting a file from staging to active source
```julia
# In src/QuantumFurnace.jl, add:
using LsqFit   # NEW: before include("fitting.jl")

# In the include block, add:
include("fitting.jl")    # PROMOTED from staging
include("mixing.jl")     # NEW

# In the export block, add:
export fit_exponential_decay, FitResult
export estimate_mixing_time, MixingTimeEstimate
```

### LsqFit.jl API usage (verified from staging/fitting.jl + official docs)
```julia
# Source: src/staging/fitting.jl lines 189-210 (existing working code)
# Model function (module-level)
_exp_decay_model(t, p) = @. p[1] * exp(-p[2] * t) + p[3]

# Fitting with bounds
fit = curve_fit(_exp_decay_model, times, values, p0;
    lower=[-Inf, 0.0, -Inf],   # gap >= 0
    upper=[Inf, Inf, Inf])

# Extract results
params = coef(fit)           # Vector{Float64} of fitted parameters
resid = residuals(fit)       # Vector{Float64} of residuals
se = stderror(fit)           # Vector{Float64} of standard errors
ci = confint(fit; level=0.95) # Vector{Tuple{Float64,Float64}} confidence intervals
conv = fit.converged          # Bool
```

### Actual mixing time from data
```julia
function _find_actual_mixing_time(
    times::Vector{<:Real},
    dists::Vector{<:Real},
    target_epsilon::Union{Nothing, Float64},
)
    target_epsilon === nothing && return last(times)
    for i in eachindex(dists)
        if dists[i] <= target_epsilon
            return times[i]
        end
    end
    return nothing  # target not reached in data
end
```

### Adding LsqFit.jl to Project.toml
```toml
# [deps] section -- add:
LsqFit = "2fda8390-95c7-5789-9bda-21331edee243"

# [compat] section -- add:
LsqFit = "0.15"
```

### Test pattern for estimate_mixing_time
```julia
@testset "Mixing Time Estimation" begin
    # Generate synthetic ThermalizeResults with known exponential decay
    A_true = 1.5
    gap_true = 0.3
    C_true = 0.001
    times = collect(0.0:0.1:50.0)
    dists = A_true .* exp.(-gap_true .* times) .+ C_true

    # Create a mock ThermalizeResults
    config = make_config(Thermalize(), EnergyDomain(); mixing_time=50.0)
    mock_result = ThermalizeResults{Float64}(
        config,
        zeros(ComplexF64, DIM, DIM),  # dummy final_dm
        dists,
        times,
        Dict{Symbol, Any}(),
    )

    est = estimate_mixing_time(mock_result; skip_initial=0.1)
    @test isapprox(est.fitted_gap, gap_true; atol=0.01)
    @test est.r_squared > 0.99
    @test est.converged == true

    # Test extrapolation
    est_extrap = estimate_mixing_time(mock_result;
        skip_initial=0.1, target_epsilon=0.01, extrapolate=true)
    @test est_extrap.mixing_time_extrapolated !== nothing
    expected_tmix = -log((0.01 - C_true) / A_true) / gap_true
    @test isapprox(est_extrap.mixing_time_extrapolated, expected_tmix; rtol=0.05)
end
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| LsqFit in staging (dormant) | LsqFit promoted to active | Phase 42 | Fitting functions become part of public API |
| No mixing time output | MixingTimeEstimate struct | Phase 42 | Users get fitted gap + mixing time from ThermalizeResults |
| Manual trace distance inspection | Quality-gated exponential fit | Phase 42 | Automated, validated mixing time with confidence intervals |

**Deprecated/outdated:**
- `staging/fitting.jl` location: will be replaced by `src/fitting.jl`
- Commented STAGING exports in QuantumFurnace.jl: will become real exports

## Open Questions

1. **MIX-02 vs "post-processing" prior decision: What does "simulation stops early" mean?**
   - What we know: The prior decision says "Mixing time estimation is a post-processing function, not embedded in run_thermalize." But MIX-02 says "extrapolate=true keyword causes simulation to stop when fit is reliable."
   - What's unclear: Does MIX-02 require modifying run_thermalize, or does it mean "estimate_mixing_time returns the extrapolated time, and the user can choose to run shorter simulations"?
   - Recommendation: Interpret MIX-02 as post-processing. `extrapolate=true` computes the extrapolated mixing time from the fit on existing data. It does NOT control run_thermalize execution. The "simulation stops early" language likely refers to the user's workflow: they run a shorter simulation, call estimate_mixing_time with extrapolate=true, and get the predicted mixing time without running to full convergence. If literal early stopping is required, implement it as a separate wrapper function (e.g., `run_thermalize_with_extrapolation`) that calls run_thermalize in segments and periodically fits.

2. **MIX-05: Should ThermalizeResults be modified or should MixingTimeEstimate be separate?**
   - What we know: ThermalizeResults is an immutable struct with BSON serialization. Adding fields breaks backward compatibility.
   - What's unclear: Does MIX-05 require literal struct modification?
   - Recommendation: Keep ThermalizeResults unchanged. MixingTimeEstimate is a separate return type from estimate_mixing_time. If the user wants to persist mixing time info, they can store the MixingTimeEstimate separately or add to the ThermalizeResults metadata Dict before saving.

3. **Should gap_estimation.jl also be promoted from staging?**
   - What we know: gap_estimation.jl uses FitResult from fitting.jl. Once fitting.jl is promoted, gap_estimation.jl could also be promoted.
   - What's unclear: Is this in scope for Phase 42? The requirements only mention fitting.jl promotion.
   - Recommendation: Out of scope for Phase 42. gap_estimation.jl stays in staging. It can be promoted in a future phase that focuses on trajectory-based spectral gap estimation.

4. **What is the LsqFit.jl UUID for Project.toml?**
   - What we know: The UUID is `2fda8390-95c7-5789-9bda-21331edee243` (from Julia General Registry).
   - Confidence: HIGH (verified from the official JuliaNLSolvers/LsqFit.jl GitHub repository).

## Sources

### Primary (HIGH confidence)
- `src/staging/fitting.jl` -- Complete, tested implementation of `fit_exponential_decay`, `FitResult`, `_log_linear_initial_guess`, `_compute_r_squared`
- `src/staging/gap_estimation.jl` -- `estimate_spectral_gap` using FitResult (reference for composition pattern)
- `src/structs.jl` -- ThermalizeResults struct definition (lines 200-207)
- `src/furnace.jl` -- `run_thermalize` implementation (lines 143-256)
- `src/QuantumFurnace.jl` -- Module structure, exports, include order
- `src/results.jl` -- BSON serialization for ThermalizeResults
- `test/staging/test_fitting.jl` -- Complete test suite for fit_exponential_decay
- `test/test_helpers.jl` -- Test fixtures (N3_HAM, N3_JUMPS, make_config, etc.)
- `test/runtests.jl` -- Test inclusion pattern
- `Project.toml` -- Current dependencies (LsqFit NOT present)
- `.planning/research/FEATURES.md` -- Detailed specification for mixing time estimation (TS-04, DIFF-02, DIFF-03)
- `.planning/REQUIREMENTS.md` -- MIX-01 through MIX-08 requirements
- `.planning/ROADMAP.md` -- Phase 42 success criteria and dependencies

### Secondary (MEDIUM confidence)
- [LsqFit.jl GitHub](https://github.com/JuliaNLSolvers/LsqFit.jl) -- v0.15.1 (April 2025), latest version confirmed
- [LsqFit.jl API docs](https://julianlsolvers.github.io/LsqFit.jl/latest/api/) -- curve_fit signature, LsqFitResult type
- [LsqFit.jl Getting Started](https://julianlsolvers.github.io/LsqFit.jl/latest/getting_started/) -- confidence_interval, standard_error, estimate_covar functions
- [LsqFit.jl Tutorial](https://julianlsolvers.github.io/LsqFit.jl/latest/tutorial/) -- Levenberg-Marquardt algorithm, weight matrix support, bounds

### Tertiary (LOW confidence)
- None -- all findings verified against codebase or official docs

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - LsqFit.jl is the established Julia nonlinear fitting library, already used in staging code
- Architecture: HIGH - Pattern is clear from existing codebase (post-processing wrapper around fit_exponential_decay)
- Pitfalls: HIGH - All identified from direct code reading and prior phase experience
- LsqFit API: HIGH - Verified against existing working code in staging/fitting.jl and official docs
- Extrapolation formula: HIGH - Standard exponential decay algebra, documented in FEATURES.md
- Quality gates: HIGH - Thresholds specified in requirements (R^2 < 0.95, offset vs epsilon)
- MIX-02 interpretation: MEDIUM - Tension between requirement wording and prior decision, documented as open question

**Research date:** 2026-03-01
**Valid until:** 2026-04-01 (stable domain -- LsqFit API unlikely to change, codebase patterns well-established)
