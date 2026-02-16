# Stack Research: Spectral Gap Estimation from Trajectory Observable Decay

**Domain:** Nonlinear curve fitting, statistical inference on decay rates, spectral gap cross-validation for quantum Lindbladian simulation
**Researched:** 2026-02-16
**Confidence:** HIGH (LsqFit.jl verified via official docs + GitHub; Distributions.jl verified; integration points confirmed from codebase analysis)

## Scope

This stack research covers ONLY the additions needed for spectral gap estimation via trajectory-based observable decay fitting:

1. Exponential decay curve fitting: `f(t) = A * exp(-gap * t) + c`
2. Parameter uncertainty / confidence intervals on the fitted gap
3. Robust initial guess strategies for the Levenberg-Marquardt solver
4. Statistical tools for handling noisy trajectory-averaged data
5. Cross-validation against exact Liouvillian eigenvalues (Arpack)

It does NOT re-research the existing stack. See the v1.2 codebase STACK.md for Arpack, FINUFFT, LinearAlgebra, BSON, StableRNGs, HypothesisTests, StatsBase, threading, convergence monitoring, etc.

## Current Stack (Relevant Subset for This Feature)

These existing capabilities are directly used and need no changes:

| Existing Capability | Role in Spectral Gap Estimation | Status |
|---------------------|-------------------------------|--------|
| `run_trajectories(..., observables=..., save_every=N)` | Generates time series `<O_i>(t)` averaged over trajectories | Keep -- produces the raw data for fitting |
| `TrajectoryResult.measurements_mean` | Matrix{Float64} of shape (n_obs, n_saves) -- observable means at each save point | Keep -- this is the input to the fitter |
| `TrajectoryResult.times` | Vector{Float64} of time points corresponding to saves | Keep -- independent variable for fitting |
| `run_lindbladian(...)` -> `LindbladianResult.spectral_gap` | Exact spectral gap via Arpack shift-invert for n=4,6 | Keep -- ground truth for cross-validation |
| `build_convergence_observables(ham, n)` | Builds ZZ_ij and H observables in eigenbasis | Keep -- defines which observables to track |
| `_compute_gibbs_observable_values(gibbs, observables)` | Computes equilibrium values `<O_i>_gibbs` | Keep -- defines the asymptotic offset `c` in the fit |
| HypothesisTests (test extra) | OneSampleTTest for statistical validation | Keep -- useful for testing gap estimate consistency |
| StatsBase (test extra) | mean, std for trajectory statistics | Keep -- useful in test assertions |

## New Production Dependencies (src/)

### Core: LsqFit.jl -- Nonlinear Curve Fitting

| Technology | Version | Purpose | Why Recommended |
|------------|---------|---------|-----------------|
| LsqFit.jl | >= 0.15 | Levenberg-Marquardt nonlinear least squares fitting for `f(t) = A * exp(-gap * t) + c` | **The** standard Julia package for nonlinear curve fitting. Pure Julia. Provides `curve_fit` with parameter bounds (essential: gap > 0), weighted fitting (essential: weight by 1/variance at each time point), and built-in `confidence_interval` / `standard_error` for parameter uncertainty without needing a separate statistics package. Uses ForwardDiff for automatic Jacobian computation. Latest release v0.15.1 (April 2025), actively maintained by JuliaNLSolvers. |

**Confidence: HIGH** -- Verified via [official docs](https://julianlsolvers.github.io/LsqFit.jl/latest/), [GitHub v0.15.1 release](https://github.com/JuliaNLSolvers/LsqFit.jl), and [tutorial](https://julianlsolvers.github.io/LsqFit.jl/latest/tutorial/).

**Why LsqFit over alternatives:**

1. **Built-in parameter uncertainty** -- `confidence_interval(fit, alpha)` computes t-distribution-based CIs directly from the Jacobian at the solution. No need to write bootstrap code or pull in a separate statistics package for the primary use case.

2. **Parameter bounds** -- `curve_fit(model, t, y, p0; lower=[...], upper=[...])` constrains the gap to be positive and the amplitude to be physically reasonable. Essential because unconstrained Levenberg-Marquardt can converge to negative decay rates on noisy data.

3. **Weighted fitting** -- Pass `wt = 1 ./ variance_per_timepoint` to properly handle trajectory noise that varies with time (early time points have high variance from diverse initial conditions; late time points have low signal).

4. **Minimal dependency footprint** -- LsqFit depends on Distributions, ForwardDiff, NLSolversBase, StatsAPI, LinearAlgebra, Printf. Of these, ForwardDiff, NLSolversBase, StatsAPI, LinearAlgebra, and Printf are already in the Manifest.toml as transitive dependencies of Optim (which is a direct dep). Adding LsqFit only truly adds **LsqFit itself + Distributions.jl** as new resolved packages.

**Key API for this milestone:**

```julia
using LsqFit

# Model: observable deviation decays exponentially to equilibrium
# p = [A, gap, c] where A = amplitude, gap = decay rate, c = offset
model(t, p) = p[1] .* exp.(-p[2] .* t) .+ p[3]

# Fit with bounds: gap > 0, A can be any sign, c unconstrained
fit = curve_fit(model, times, obs_deviation, p0;
                lower = [-Inf, 0.0, -Inf],
                upper = [Inf, Inf, Inf])

# Extract gap estimate and uncertainty
gap_estimate = fit.param[2]
gap_stderr = standard_error(fit)[2]
gap_ci = confidence_interval(fit, 0.05)[2]  # 95% CI on gap

# Covariance matrix for full parameter correlation analysis
cov_matrix = estimate_covar(fit)
```

### Transitive: Distributions.jl -- Brought in by LsqFit

| Technology | Version | Purpose | Why Acceptable |
|------------|---------|---------|----------------|
| Distributions.jl | >= 0.25 | t-distribution quantiles for `confidence_interval()` in LsqFit; not used directly in QuantumFurnace code | LsqFit uses `quantile(TDist(dof), 1 - alpha/2)` internally to compute confidence intervals. This is a well-maintained, standard Julia statistics package (v0.25.123, Jan 2026). Compatible with Julia 1.11+. It brings some transitive deps (PDMats, StatsFuns, FillArrays, etc.) but these are lightweight and standard in the Julia ecosystem. |

**Confidence: HIGH** -- Verified via [Distributions.jl v0.25 docs](https://juliastats.org/Distributions.jl/v0.25/) generated with Julia 1.11.7. [v0.25.123 on Zenodo](https://zenodo.org/records/18145493).

## No New Test Dependencies Needed

The existing test extras are sufficient:

| Existing Test Dep | Role in Spectral Gap Tests |
|-------------------|---------------------------|
| HypothesisTests | `OneSampleTTest` to verify fitted gap matches exact Arpack gap within statistical tolerance |
| StatsBase | `mean`, `std` for computing reference statistics on fit residuals |
| StableRNGs | Deterministic seeding for reproducible trajectory data in gap estimation tests |

## How LsqFit Integrates with Existing Code

### Data Flow: Trajectory -> Fit -> Gap

```julia
# STEP 1: Run trajectories with observable measurement (EXISTING)
observables, obs_names = build_convergence_observables(hamiltonian, n_qubits)
result = run_trajectories(jumps, config, psi0, hamiltonian;
    observables = observables,
    save_every = save_every,
    ntraj = ntraj,
    seed = seed)

# result.times        :: Vector{Float64}      -- time grid
# result.measurements_mean :: Matrix{Float64}  -- (n_obs, n_saves)

# STEP 2: Compute deviation from equilibrium (NEW, uses existing helpers)
obs_gibbs = _compute_gibbs_observable_values(gibbs, observables)
# For observable i: deviation(t) = <O_i>(t) - <O_i>_gibbs

# STEP 3: Fit exponential decay (NEW, uses LsqFit)
using LsqFit

model(t, p) = p[1] .* exp.(-p[2] .* t) .+ p[3]

for i in 1:n_obs
    y = result.measurements_mean[i, :] .- obs_gibbs[i]
    # Or fit the raw data with c as free parameter:
    # model_with_offset(t, p) = p[1] * exp(-p[2] * t) + p[3]
    # y = result.measurements_mean[i, :]
    # p0 = [y[1] - y[end], initial_gap_guess, y[end]]

    p0 = _compute_initial_guess(result.times, y)
    fit = curve_fit(model, result.times, y, p0;
                    lower = [-Inf, 0.0, -Inf])
    gap_i = fit.param[2]
    gap_ci_i = confidence_interval(fit, 0.05)[2]
end

# STEP 4: Cross-validate against exact gap (EXISTING Arpack result)
liouv = run_lindbladian(jumps, liouv_config, hamiltonian; trotter=trotter)
exact_gap = abs(real(liouv.spectral_gap))
```

### Integration Points in the Codebase

| Integration Point | File | What Changes |
|-------------------|------|-------------|
| New function: `estimate_spectral_gap(result, gibbs, observables)` | New file: `src/spectral_gap.jl` | Takes TrajectoryResult + Gibbs state + observables, returns gap estimate + CI + per-observable fits |
| New struct: `SpectralGapResult` | `src/structs.jl` | Holds gap estimate, confidence interval, per-observable fit quality (R-squared, residual norm), fit parameters |
| Existing: `run_trajectories` with `observables` + `save_every` | `src/trajectories.jl` | No changes needed -- already produces the exact data format required |
| Existing: `run_lindbladian` -> `LindbladianResult.spectral_gap` | `src/furnace.jl` | No changes needed -- provides ground truth for validation |
| New export: `estimate_spectral_gap` | `src/QuantumFurnace.jl` | Add to exports |
| New `using LsqFit` | `src/QuantumFurnace.jl` | Add import |

## Robust Initial Guess Strategy

The Levenberg-Marquardt algorithm in LsqFit converges reliably IF the initial guess is reasonable. For noisy trajectory data, a poor initial guess causes convergence to local minima or failure. **This is the most likely pitfall**, not the fitting library itself.

### Recommended: Log-Linear Preprocessing for Initial Guess

```julia
function _compute_initial_guess(times, deviation)
    # deviation = <O>(t) - <O>_gibbs, should decay from ~A to ~0

    # Estimate offset c from late-time average (last 20% of data)
    n = length(deviation)
    c_guess = mean(deviation[max(1, round(Int, 0.8*n)):n])

    # Subtract offset, take abs to handle sign, take log
    shifted = deviation .- c_guess
    # Use only points where |shifted| > threshold to avoid log(~0)
    mask = abs.(shifted) .> 0.01 * maximum(abs.(shifted))
    if sum(mask) < 3
        # Fallback: crude guess
        A_guess = deviation[1] - c_guess
        gap_guess = 1.0 / (times[end] - times[1])
        return [A_guess, gap_guess, c_guess]
    end

    # Log-linear fit: log|shifted| = log|A| - gap * t
    log_y = log.(abs.(shifted[mask]))
    t_masked = times[mask]
    # Simple linear regression
    t_mean = mean(t_masked)
    y_mean = mean(log_y)
    slope = sum((t_masked .- t_mean) .* (log_y .- y_mean)) / sum((t_masked .- t_mean).^2)
    intercept = y_mean - slope * t_mean

    A_guess = sign(shifted[1]) * exp(intercept)
    gap_guess = max(-slope, 1e-6)  # Ensure positive

    return [A_guess, gap_guess, c_guess]
end
```

**Why log-linear for initial guess only (not final fit):**
- Log transform distorts error distribution (additive Gaussian noise on linear scale becomes non-Gaussian on log scale)
- Log transform cannot handle sign changes in the deviation (which happen when noise pushes the observable past equilibrium)
- But for initial guess, approximate is fine -- the nonlinear LM fit refines from there

**Confidence: HIGH** -- Standard numerical methods practice. The log-linear initial guess + nonlinear refinement pattern is used in scipy.optimize.curve_fit documentation, MATLAB's fit() documentation, and countless scientific computing references.

## Weighted Fitting for Trajectory Noise

Trajectory-averaged observable values have non-uniform variance along the time axis:

- **Early times (t << 1/gap):** High variance because trajectories haven't mixed -- each trajectory is near the initial state but random jumps create spread.
- **Late times (t >> 1/gap):** Low variance because trajectories have converged near the thermal equilibrium. The signal (deviation from Gibbs) is also small, so SNR is low.

### Weight Computation from Trajectory Batches

```julia
# Option A: If running multiple independent batches (e.g., from run_trajectories_convergence)
# Compute per-timepoint variance across batches
# wt = 1 ./ var_per_timepoint

# Option B: If running a single large batch with save_every
# Split trajectories into sub-batches, compute variance of means
# This requires modifying run_trajectories to return per-sub-batch data

# Option C (RECOMMENDED for first implementation): Unweighted fit
# Start with unweighted fit. Add weighting later if fit quality is poor.
# Reason: For n=4,6 with ntraj >= 10000, the trajectory average is smooth
# enough that unweighted fitting gives good results. Weighting is a
# refinement for when the gap estimate precision matters.
```

**Recommendation:** Start with unweighted fitting. The unweighted fit is simpler, and for the cross-validation against exact eigenvalues (the primary goal), the gap estimate only needs to be within ~10-20% of the exact value to demonstrate the method works. Add variance-weighted fitting as a refinement if unweighted residuals show systematic structure.

## Multi-Exponential Considerations

The observable deviation from equilibrium is not a pure single exponential. The exact expansion is:

```
<O>(t) - <O>_gibbs = sum_k c_k * exp(lambda_k * t)
```

where `lambda_k` are Liouvillian eigenvalues (all with Re(lambda_k) <= 0) and `c_k` depend on the overlap of the initial state and observable with the k-th eigenmode.

### Why Single Exponential Works for Gap Estimation

At long times (t >> 1/|lambda_3|), all faster modes have decayed, leaving:

```
<O>(t) - <O>_gibbs ≈ c_2 * exp(lambda_2 * t)   for large t
```

where `lambda_2` is the spectral gap (smallest non-zero eigenvalue in magnitude). **The gap dominates the late-time behavior**, which is precisely what we want to extract.

### Fitting Window Selection

Do NOT fit the entire time series. Instead:

1. **Skip early transient:** The first ~10% of the time series contains multi-exponential contributions from fast modes. Skip it.
2. **Skip late noise floor:** Once the deviation drops below the trajectory noise level, the data is pure noise. Stop there.
3. **Fit the middle region** where the single exponential dominates.

```julia
function _select_fitting_window(times, deviation; skip_fraction=0.1, noise_floor_ratio=0.05)
    n = length(times)
    t_start = round(Int, skip_fraction * n) + 1
    max_dev = maximum(abs.(deviation))
    t_end = findlast(abs.(deviation) .> noise_floor_ratio * max_dev)
    t_end = max(t_end, t_start + 5)  # Need at least a few points
    return t_start:t_end
end
```

**Confidence: HIGH** -- This is standard practice in spectroscopy, NMR relaxometry, and fluorescence lifetime imaging where exponential decay rates are extracted from noisy data with multi-exponential contributions.

### When Single Exponential Fails

If the observable has zero overlap with the gap eigenmode (`c_2 = 0`), the fitted rate will correspond to `lambda_3` or a later eigenvalue. This is why fitting multiple observables (ZZ correlations + energy) is essential -- at least one will have non-zero gap overlap.

**Detection:** Compare gap estimates across observables. If they agree, the estimate is the true gap. If they disagree, the smallest positive fitted rate (closest to zero) is the best gap estimate, because faster rates indicate the observable couples to higher modes.

## What NOT to Add

| Avoid | Why | What to Do Instead |
|-------|-----|-------------------|
| **Bootstrap.jl** for confidence intervals | LsqFit already provides analytically-derived `confidence_interval` from the Jacobian covariance matrix. Bootstrap resampling would require either (a) storing per-trajectory data (massive memory for ntraj=10000+), or (b) re-running trajectories (computationally prohibitive). The Jacobian-based CI is standard and appropriate for this application. | Use `confidence_interval(fit, alpha)` from LsqFit. If Jacobian-based CIs are insufficient (e.g., highly non-Gaussian residuals), implement a simple residual bootstrap on the fit residuals, which requires no extra package. |
| **Optim.jl for custom fitting** | Already a dependency, but LsqFit wraps the Levenberg-Marquardt algorithm with purpose-built curve fitting infrastructure (Jacobian management, covariance estimation, confidence intervals). Reimplementing this with raw Optim would be error-prone and redundant. | Use LsqFit which internally uses NLSolversBase (same optimization base as Optim). |
| **SignalDecomposition.jl / Prony method packages** | Prony and matrix pencil methods are theoretically superior for multi-exponential decomposition, but (a) no mature Julia package exists, (b) the single-exponential fit with window selection is sufficient for gap estimation, (c) Prony methods require equispaced samples (which we have) but are notoriously sensitive to noise without careful SVD-based regularization. | Fit single exponential with window selection. If multi-exponential decomposition is needed later, implement a simple SVD-based Prony method in ~50 lines rather than adding an external dependency. |
| **GLM.jl / StatsModels.jl** | These are for linear statistical models and generalized linear models. The exponential decay fit is inherently nonlinear. | Use LsqFit for nonlinear fitting. |
| **Measurements.jl** | Propagates uncertainties through calculations via dual numbers. Elegant but heavyweight for this use case where we only need CI on one parameter (the gap). | Use `standard_error(fit)[2]` from LsqFit directly. |
| **Multi-exponential fit (fitting 2+ decay rates)** | Tempting for extracting lambda_2 AND lambda_3. But 6+ parameter fits (A1, gap1, A2, gap2, A3, c) on noisy data are notoriously ill-conditioned. The gap/lambda_3 ratio determines identifiability; for Lindbladians this ratio is often < 2, making the two rates nearly indistinguishable from noisy data. | Fit single exponential with late-time window. Extract gap from the dominant late-time decay. Cross-validate against exact Arpack eigenvalues. |
| **Distributions.jl as direct dependency** | It is a transitive dependency of LsqFit. Do not import it directly in QuantumFurnace -- the confidence interval computation happens inside LsqFit. | Let LsqFit handle the t-distribution internally. If you ever need Distributions directly (e.g., for custom statistical tests), add it then. |

## Recommended Stack Summary

### Production Dependencies to ADD to `[deps]`

| Package | UUID | Purpose |
|---------|------|---------|
| LsqFit | `2fda8390-95c7-5789-9bda-21331edee243` | Nonlinear curve fitting for exponential decay gap estimation |

That is **one** new direct dependency. The transitive deps (Distributions, ForwardDiff, NLSolversBase, StatsAPI) are either already resolved or lightweight and standard.

### No New Test Dependencies

HypothesisTests and StatsBase (already in test extras) provide everything needed for validating gap estimates against exact eigenvalues.

## Installation

```julia
using Pkg
Pkg.activate(".")
Pkg.add("LsqFit")
```

Concrete `Project.toml` additions:

```toml
[deps]
# ... existing deps ...
LsqFit = "2fda8390-95c7-5789-9bda-21331edee243"

[compat]
# ... existing compat ...
LsqFit = "0.15"
```

**Note on compat range:** Pin to `"0.15"` (not `"0.14, 0.15"`) because v0.15 introduced important improvements to the Levenberg-Marquardt implementation and parameter bounds handling. The package requires Julia >= 1.6, well within the project's Julia >= 1.11 requirement.

## Alternatives Considered

| Recommended | Alternative | When to Use Alternative |
|-------------|-------------|-------------------------|
| LsqFit.jl `curve_fit` | Optim.jl + manual Jacobian | Only if you need a custom loss function (e.g., robust L1 loss instead of L2). For standard least squares, LsqFit is strictly better because it handles Jacobian computation, covariance estimation, and confidence intervals automatically. |
| LsqFit.jl `confidence_interval` | Manual bootstrap resampling | Only if residuals are severely non-Gaussian (would indicate a model misspecification rather than a statistics problem). For Monte Carlo trajectory data with sufficient ntraj, the central limit theorem ensures approximately Gaussian residuals. |
| Single exponential fit with window selection | Multi-exponential fit via Prony/ESPRIT | Only if you need to extract multiple eigenvalues simultaneously. For gap estimation alone, the single exponential approach is more robust and well-conditioned. Consider this for a future milestone if eigenvalue spectrum characterization is needed. |
| Unweighted initial fit | Variance-weighted fit from sub-batches | Add weighting as a refinement AFTER the unweighted approach is validated. Weighting requires either storing per-trajectory-batch data or modifying run_trajectories to return variance estimates. |
| LsqFit.jl | SciPy via PythonCall.jl | Never. Adding a Python dependency for curve fitting when a mature pure-Julia solution exists is unjustifiable overhead. |

## Stack Patterns by Use Case

**If estimating the gap for n=4 (dim=16, validation target):**
- Run ntraj >= 5000 with save_every = 10, total_time = 5/expected_gap
- Unweighted single exponential fit on late-time window
- Cross-validate against exact Arpack gap from `run_lindbladian`
- LsqFit confidence_interval gives error bars on the estimate

**If estimating the gap for n=6 (dim=64, validation target):**
- Run ntraj >= 10000 with save_every = 20, total_time = 5/expected_gap
- Same fitting approach as n=4
- Exact Arpack gap still available (64^2 = 4096 dim Liouvillian)
- May need weighted fit if signal-to-noise is poor at late times

**If estimating the gap for n=8+ (dim=256+, where Arpack is infeasible):**
- This is the production use case: trajectory-based gap estimation IS the method
- Run ntraj >= 50000 with save_every = 50
- Fit multiple observables (all ZZ correlations + energy)
- Take gap estimate as minimum fitted rate across observables with good fit quality
- No exact cross-validation available; rely on agreement across observables and CI overlap

## Version Compatibility

| Package | Compatible With | Notes |
|---------|-----------------|-------|
| LsqFit 0.15.1 | Julia >= 1.6 | Latest stable (April 2025). Pure Julia. No binary dependencies. |
| LsqFit 0.15.1 | Distributions >= 0.18 | LsqFit compat allows 0.18-0.25. Project will resolve 0.25.x. |
| LsqFit 0.15.1 | ForwardDiff >= 0.10 or 1.0 | Already in Manifest via Optim. No conflict. |
| LsqFit 0.15.1 | NLSolversBase 7.5 | Already in Manifest via Optim. No conflict. |
| LsqFit 0.15.1 | StatsAPI 1.x | Already in Manifest. No conflict. |
| Distributions 0.25.x | Julia 1.11, 1.12 | Verified: docs generated with Julia 1.11.7. |

**No version conflicts expected.** The shared transitive dependencies (ForwardDiff, NLSolversBase, StatsAPI) are already resolved in the Manifest.toml via Optim. Adding LsqFit will not trigger version conflicts.

## Sources

- [LsqFit.jl GitHub repository](https://github.com/JuliaNLSolvers/LsqFit.jl) -- v0.15.1 (April 2025), release history, dependency list. HIGH confidence.
- [LsqFit.jl official documentation](https://julianlsolvers.github.io/LsqFit.jl/latest/) -- API reference, tutorials. HIGH confidence.
- [LsqFit.jl Getting Started](https://julianlsolvers.github.io/LsqFit.jl/latest/getting_started/) -- confidence_interval, standard_error, estimate_covar API. HIGH confidence.
- [LsqFit.jl Tutorial](https://julianlsolvers.github.io/LsqFit.jl/latest/tutorial/) -- Exponential model example, parameter bounds, weighted fitting. HIGH confidence.
- [LsqFit.jl Project.toml](https://raw.githubusercontent.com/JuliaNLSolvers/LsqFit.jl/master/Project.toml) -- Dependencies: Distributions, ForwardDiff, NLSolversBase, StatsAPI, LinearAlgebra, Printf. Julia >= 1.6. HIGH confidence.
- [Distributions.jl v0.25 docs](https://juliastats.org/Distributions.jl/v0.25/) -- Generated with Julia 1.11.7, confirming compatibility. HIGH confidence.
- [Distributions.jl v0.25.123 on Zenodo](https://zenodo.org/records/18145493) -- Latest release Jan 2026. HIGH confidence.
- [HypothesisTests.jl parametric tests](https://juliastats.org/HypothesisTests.jl/stable/parametric/) -- OneSampleTTest for gap validation. HIGH confidence.
- [Bootstrap.jl](https://github.com/juliangehring/Bootstrap.jl) -- Evaluated but not recommended (LsqFit CI is sufficient). HIGH confidence evaluation.
- [Julia Discourse: weighted LsqFit](https://discourse.julialang.org/t/how-to-do-weighted-least-squares-fit-with-lsqfit-jl/61136) -- Weight parameter usage patterns. MEDIUM confidence.
- Standard exponential decay fitting methodology from NMR/spectroscopy literature -- log-linear initial guess, window selection, multi-exponential considerations. HIGH confidence (textbook methods).

---
*Stack research for: QuantumFurnace.jl Spectral Gap Estimation from Trajectory Observable Decay*
*Researched: 2026-02-16*
