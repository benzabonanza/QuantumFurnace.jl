# Feature Landscape: v1.3 Spectral Gap Estimation from Trajectory Observables

**Domain:** Spectral gap estimation for Lindbladian dynamics via exponential decay fitting of trajectory-averaged observables
**Researched:** 2026-02-16
**Confidence:** HIGH (physics well-established, codebase analysis thorough, fitting methodology standard)

## Scope

This research covers features needed for **v1.3 Mixing Time Estimation**: estimating the Lindbladian spectral gap (the real part of the second eigenvalue of the Liouvillian superoperator) from trajectory-based observable time series. The idea is physically grounded: under Lindbladian evolution, any observable expectation value `<O>(t) = tr(O rho(t))` decays exponentially toward its thermal equilibrium value at a rate determined by the spectral gap, provided the observable has nonzero overlap with the gap mode. By fitting `A * exp(-lambda * t) + C` to the trajectory-averaged observable time series, we extract `lambda` as an estimate of the spectral gap. Cross-validation against the exact Liouvillian eigenvalue (available for n=4,6 via `run_lindbladian`) establishes trust; the same fitting procedure then extends to n>8 where full Liouvillian construction is infeasible.

---

## Physics Foundation: Observable Decay and the Spectral Gap

### The Spectral Decomposition Argument

The Lindbladian generator L has eigenvalues {0, lambda_1, lambda_2, ...} where Re(lambda_i) < 0 for all i > 0 (for a primitive Lindbladian with unique steady state). The spectral gap is defined as:

```
Delta = min_i |Re(lambda_i)|  for lambda_i != 0
```

Under the semigroup evolution rho(t) = e^{Lt} rho(0), the density matrix decomposes as:

```
rho(t) = rho_ss + sum_i c_i * R_i * exp(lambda_i * t)
```

where `rho_ss` is the steady state (Gibbs state for our Lindbladian) and `R_i` are right eigenmatrices. For any observable O:

```
<O>(t) = tr(O rho(t)) = <O>_ss + sum_i c_i * tr(O R_i) * exp(lambda_i * t)
```

At long times, the slowest-decaying mode dominates:

```
<O>(t) ~ <O>_ss + A * exp(-Delta * t)     (for t >> 1/|Re(lambda_2)|)
```

where `A = c_1 * tr(O R_1)` depends on:
1. The overlap of the initial state with the gap mode (`c_1`)
2. The overlap of the observable with the gap mode (`tr(O R_1)`)

Both must be nonzero for the spectral gap to be visible in observable O's decay.

### Why This Works for QuantumFurnace

- The Lindbladian constructed by QuantumFurnace satisfies KMS detailed balance, which guarantees a unique steady state (the Gibbs state) and a real spectral gap for the symmetrized Lindbladian.
- The initial state is the maximally mixed state (or a random pure state), which has nonzero overlap with all eigenspaces, so `c_1 != 0`.
- The gap mode for spin chain Lindbladians typically has significant overlap with energy and magnetization observables (confirmed by `LindbladianResult.gap_mode` from `run_lindbladian`).

### Key Subtlety: Real vs. Complex Eigenvalues

The Lindbladian eigenvalues are generally complex: `lambda_i = -gamma_i + i * omega_i`. The real part `gamma_i` determines the decay rate; the imaginary part `omega_i` produces oscillations. For KMS-detailed-balanced Lindbladians, the eigenvalues come in conjugate pairs, but the gap eigenvalue may be real or complex depending on the Hamiltonian.

For the fitting:
- If the gap eigenvalue is real (`omega = 0`): pure exponential decay, `A * exp(-gamma * t) + C`
- If the gap eigenvalue is complex (`omega != 0`): damped oscillation, `A * exp(-gamma * t) * cos(omega * t + phi) + C`
- For KMS-detailed-balanced Lindbladians with real coupling to the gap mode (which is typical for energy and magnetization observables), the dominant contribution is real exponential decay.

The existing `LindbladianResult.spectral_gap` is complex, but for the Heisenberg chain with disordering field, `|Im(spectral_gap)| / |Re(spectral_gap)|` is typically small (< 0.1), so the pure exponential model is a good first approximation. A damped oscillation model should be available as a fallback.

---

## Table Stakes

Features users expect. Missing any of these means the milestone is incomplete.

| Feature | Why Expected | Complexity | Dependencies |
|---------|--------------|------------|--------------|
| **Observable-only trajectory runner** | Current `run_trajectories` with `observables` and `save_every` already returns time-resolved `<O>(t)`. But `run_trajectories_convergence` does batch-level DM reconstruction without time-resolved measurements. Need a runner that does time-resolved observable measurements without per-batch DM overhead -- just accumulate `<O>(t)` across many trajectories. | LOW | Existing `run_trajectories` with `observables` parameter; minor API design |
| **Total magnetization observable** | `sum_i Z_i` is the canonical order parameter for spin chains. Must be available alongside ZZ correlations and energy. Not currently in `build_convergence_observables`. Easy to add. | LOW | `build_convergence_observables`, `pad_term`, Pauli `Z` |
| **Single-exponential fit: `A * exp(-gap * t) + C`** | The core fitting model. Extracts the spectral gap from observable time series. Must handle: initial guess from data, weighted least squares with 1/variance weights, parameter bounds (gap > 0, C near Gibbs value). | MEDIUM | `LsqFit.jl` (new dependency), observable time series data |
| **Automatic initial guess for fitting** | Levenberg-Marquardt is sensitive to initial guess. Must derive good initial parameters from the data: `C_guess = mean(last 20% of data)` (plateau), `A_guess = data[1] - C_guess`, `gap_guess` from log-linear fit of `data - C_guess`. | LOW | Observable time series data |
| **Fit quality metrics** | R-squared, residual norm, and confidence intervals on the gap estimate. Without these, cannot assess whether the fit is trustworthy. LsqFit.jl provides `confidence_interval`, `standard_error`, `estimate_covar`. | LOW | LsqFit.jl fit result |
| **Cross-validation against exact Liouvillian gap** | For n=4 and n=6, compute the exact spectral gap via `run_lindbladian`, then compare against the trajectory-fitted gap. This is the validation that makes the method trustworthy for n>8. | MEDIUM | `run_lindbladian` (existing), fitted gap, matching configs |
| **`estimate_spectral_gap` function** | A single callable function that takes trajectory observable data and returns a gap estimate with confidence interval. This is the reusable API for larger systems. Signature: `estimate_spectral_gap(times, observable_data; kwargs...) -> (gap, confidence_interval, fit_result)` | MEDIUM | All fitting infrastructure |
| **Fitting window selection** | The exponential fit should exclude early transient (first ~10% of data where multi-exponential effects dominate) and use data from the "single-exponential regime" where only the gap mode contributes. Must detect or allow specifying the fitting window. | MEDIUM | Observable time series data, fit quality metrics |

## Differentiators

Features that add substantial value but are not strictly required for milestone sign-off.

| Feature | Value Proposition | Complexity | Dependencies |
|---------|-------------------|------------|--------------|
| **Multi-observable gap consistency check** | Fit the gap independently from energy, magnetization, and each ZZ correlation. If the same gap is recovered from all observables (within error bars), this strongly validates the method. Disagreement flags observables with poor gap mode overlap. | LOW | Multiple fitted gaps + comparison logic |
| **Damped oscillation model: `A * exp(-gamma * t) * cos(omega * t + phi) + C`** | Handles complex eigenvalues where the gap mode produces oscillatory decay. Falls back to pure exponential when omega is small. Paper completeness: shows the method handles the general case. | MEDIUM | LsqFit.jl, model selection (AIC/BIC or F-test) |
| **Automatic model selection (exponential vs. damped oscillation)** | Use Bayesian Information Criterion (BIC) or F-test to decide whether the oscillatory model is justified by the data. Prevents overfitting with the more complex model when the gap is purely real. | MEDIUM | Both fit models, statistical model comparison |
| **Observable overlap diagnostic** | Compute `|tr(O * gap_mode)|` for each observable (available for n=4,6 where the gap mode is known from `run_lindbladian`). Reports which observables have the strongest coupling to the gap mode. Guides observable selection for larger systems. | LOW | `LindbladianResult.gap_mode`, observable matrices |
| **Multi-exponential fit: `A1 * exp(-g1 * t) + A2 * exp(-g2 * t) + C`** | Resolves the first AND second decay modes from the data. More robust gap extraction at early times before the single-exponential regime. However, multi-exponential fitting is notoriously ill-conditioned -- keep this as optional/advanced. | HIGH | LsqFit.jl with 5+ parameters, careful regularization |
| **Bootstrap confidence intervals on gap** | Resample trajectory batches (block bootstrap), refit the gap from each resample, report the bootstrap confidence interval. More robust than the linearized LsqFit confidence interval for nonlinear models. | MEDIUM | Trajectory batch data, refitting loop |
| **Gap vs. beta scaling plot** | For each (n, beta) pair, plot the estimated gap vs. beta. Theory predicts the gap decreases with beta (slower mixing at lower temperatures). This is a key paper figure. | LOW | Multiple (n, beta) experiments with fitted gaps |
| **Gap vs. n scaling plot** | For fixed beta, plot the estimated gap vs. n (system size). Theory for 1D Heisenberg chains predicts system-size-independent gap at any finite temperature (proven in arXiv:2510.08533). This is a testable prediction. | LOW | Multiple (n, beta) experiments with fitted gaps |
| **Variance-weighted fitting** | Weight the least-squares fit by `1/var(O(t_i))` where the variance is estimated from trajectory-to-trajectory fluctuations. Downweights noisy early-time data and noisy late-time data near the plateau. | LOW | Per-time-point variance from trajectory batches |

## Anti-Features

Features to explicitly NOT build for this milestone.

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| **Prony's method for multi-exponential extraction** | Prony's method solves for exponential frequencies via polynomial root-finding, which is severely ill-conditioned for noisy data. It requires uniformly sampled data and the number of exponentials must be known a priori. The trajectory data is noisy (Monte Carlo variance ~1/sqrt(N_traj)); Prony would produce unstable results. | Use nonlinear least squares (LsqFit.jl) with single or double exponential models. NLLS is robust to noise when given good initial guesses and parameter bounds. |
| **Full Liouvillian construction for n>6** | The Liouvillian is dim^2 x dim^2 = 65536x65536 for n=8, requiring ~32 GB. For n=10, it is 1M x 1M -- completely infeasible. The entire point of trajectory-based gap estimation is to avoid this. | Use trajectory-based observable fitting for n>6. Cross-validate against exact Liouvillian only for n=4,6. |
| **Imaginary-time evolution or adiabatic state preparation for gap estimation** | Recent work (arXiv:2512.19288) proposes adiabatic gap estimation on quantum devices. This is irrelevant for classical simulation -- we have direct access to observable time series. | Fit the decay rate directly from real-time observable data. |
| **Autocorrelation-based gap estimation** | Autocorrelation time estimation requires the system to be at stationarity (already in the Gibbs state), then measures decorrelation. Our trajectories start far from equilibrium and evolve toward it. The decay-to-equilibrium approach is more natural and uses the same data we already collect. | Fit the observable's approach to equilibrium, not its fluctuations around equilibrium. |
| **Stochastic trace estimation for the Liouvillian gap** | Methods like the stochastic Lanczos quadrature can estimate eigenvalue densities without full matrix construction. These are sophisticated but require matrix-vector products with the Liouvillian, which we do not currently have as a separate primitive (only the full explicit matrix from `construct_lindbladian`). | Use trajectory-based fitting. If stochastic Liouvillian methods are needed later, they are a separate milestone. |
| **Real-time plotting during fitting** | Interactive plotting during the fit loop is a development convenience, not a feature. | Post-process and plot after the fit completes. Paper plots are made from saved data. |
| **GPU-accelerated fitting** | The fitting is a tiny computation (minutes of Levenberg-Marquardt on a few hundred data points). The bottleneck is trajectory generation, not fitting. | Run LsqFit.jl on CPU. All computation time is in trajectory sampling. |

---

## Feature Details

### Observable-Only Trajectory Runner

**Current state:** `run_trajectories` already supports `observables` and `save_every`. It returns `TrajectoryResult` with `measurements_mean::Matrix{Float64}` (n_obs x n_saves) and `times::Vector{Float64}`. This is already exactly what we need for observable time series.

**What is missing:** The existing `run_trajectories_convergence` and `run_trajectories_adaptive` do NOT support time-resolved observable measurements -- they track batch-level convergence metrics only. For spectral gap estimation, we need many trajectories with time-resolved `<O>(t)` data.

**Recommended approach:** Use `run_trajectories` directly with the `observables` and `save_every` parameters. No new runner function is needed. The existing code already:
1. Saves `<O_i>(t)` at every `save_every` steps across all trajectories
2. Averages across trajectories (both serial and multi-threaded paths)
3. Returns the time grid and measurement matrix

What IS needed:
- A convenience function `build_spectral_gap_observables(hamiltonian, num_qubits)` that returns [H, M_z, ZZ_1, ZZ_2, ...] with appropriate names
- Choose `save_every` and `total_time` appropriately: `save_every` should give ~200-500 time points for smooth decay curves; `total_time` should be long enough to see the observable reach its plateau (roughly `5/gap`)

**Complexity:** LOW. Mostly wiring existing functionality.

### Total Magnetization Observable

**What:** `M_z = sum_{i=1}^{n} Z_i` in the Hamiltonian eigenbasis (or Trotter eigenbasis for TrotterDomain).

**Implementation:** Add to `build_convergence_observables`:
```julia
# Total magnetization in eigenbasis
Mz_comp = sum(Matrix{ComplexF64}(pad_term([Z], num_qubits, i)) for i in 1:num_qubits)
Mz_eigen = V' * Mz_comp * V
push!(observables, Mz_eigen)
push!(names, "Mz")
```

**Why total magnetization specifically for the gap:**
- For the 1D Heisenberg chain with disorder, the total magnetization sector structure matters. The gap mode typically lives in the `M_z = 0` sector (same as the ground state for the antiferromagnet). However, the OBSERVABLE `M_z` couples to modes that CHANGE the total magnetization -- these are typically NOT the gap mode.
- Energy `<H>` is the better observable for the gap because the gap mode is an energy-like excitation (it describes the slowest approach to the thermal energy). The gap mode lives in the same symmetry sector as the Gibbs state, and `H` has diagonal matrix elements in the eigenbasis that directly couple to it.
- ZZ correlations `<Z_iZ_{i+1}>` also couple well to the gap mode because they probe the two-body correlations that change as the system thermalizes.
- Total magnetization `M_z` is included for completeness and as a consistency check. If `<M_z>_gibbs = 0` (as for the clean Heisenberg chain), then `M_z` couples only to odd-magnetization modes and may not see the gap at all. With disorder, `<M_z>_gibbs != 0`, so there is some coupling.

**Observable selection priority for gap estimation:**
1. **Energy `<H>`** -- strongest coupling to gap mode (diagonal in eigenbasis)
2. **ZZ correlations `<Z_iZ_{i+1}>`** -- strong coupling, site-resolved
3. **Total magnetization `M_z`** -- weaker coupling, included as consistency check

**Complexity:** LOW. One function addition.

### Single-Exponential Fit

**Model:** `f(t, p) = p[1] * exp(-p[2] * t) + p[3]`

where:
- `p[1] = A` -- amplitude (signed, can be positive or negative)
- `p[2] = gap` -- decay rate (must be positive)
- `p[3] = C` -- asymptotic value (the Gibbs expectation value)

**Implementation with LsqFit.jl:**
```julia
using LsqFit

function fit_spectral_gap(times, obs_data;
    t_start_frac=0.1, t_end_frac=1.0,
    obs_gibbs=nothing, weights=nothing)

    # Select fitting window
    n = length(times)
    i_start = max(1, round(Int, t_start_frac * n))
    i_end = round(Int, t_end_frac * n)
    t_fit = times[i_start:i_end]
    y_fit = obs_data[i_start:i_end]

    # Model
    model(t, p) = p[1] .* exp.(-p[2] .* t) .+ p[3]

    # Initial guess
    C_guess = mean(y_fit[end-div(length(y_fit),5):end])  # last 20%
    A_guess = y_fit[1] - C_guess
    # Log-linear estimate for gap
    y_shifted = abs.(y_fit .- C_guess) .+ 1e-15
    log_y = log.(y_shifted)
    gap_guess = -(log_y[end] - log_y[1]) / (t_fit[end] - t_fit[1])
    gap_guess = max(gap_guess, 1e-6)  # ensure positive

    p0 = [A_guess, gap_guess, C_guess]

    # Bounds: gap must be positive
    lower = [-Inf, 0.0, -Inf]
    upper = [Inf, Inf, Inf]

    fit = curve_fit(model, t_fit, y_fit, p0; lower=lower, upper=upper)

    return fit
end
```

**Key considerations:**
- **Parameter bounds:** The gap must be positive. LsqFit.jl supports box constraints via `lower` and `upper` keyword arguments.
- **Weights:** If per-time-point variance is available (from trajectory-to-trajectory fluctuations), use `wt = 1.0 ./ variance` in `curve_fit`.
- **Gibbs value constraint:** If the Gibbs expectation value is known (it is, from `tr(gibbs * O)`), we can fix `C` to the known value and fit only `A` and `gap`. This reduces the fit to 2 parameters and improves conditioning. Provide this as an option.

**Complexity:** MEDIUM. The fitting itself is straightforward with LsqFit.jl. The complexity is in the initial guess logic, fitting window selection, and edge case handling.

### Fitting Window Selection

**The problem:** The observable time series has three regimes:
1. **Early transient** (t < ~2/|lambda_2|): Multiple exponential modes contribute. The observable decays with a rate faster than the gap (because faster-decaying modes are still active).
2. **Single-exponential regime** (2/|lambda_2| < t < ~5/gap): Only the gap mode contributes. This is where the fit should be performed.
3. **Plateau** (t > 5/gap): The observable has reached its equilibrium value. The signal is noise. Including too much plateau data biases the fit toward C and loses sensitivity to the gap.

**Approach:** Default to excluding the first 10% and fitting the remaining data. Provide an optional `t_start` parameter. For automatic window selection:
1. Compute the running derivative `d<O>/dt` numerically.
2. Find where the decay rate stabilizes (i.e., `d/dt[log|<O> - C|]` becomes approximately constant). This marks the start of the single-exponential regime.
3. Find where `|<O>(t) - C| < noise_level`. This marks the end of useful data.

For the cross-validation at n=4,6, the exact `|lambda_2|` is known, so the fitting window can be set precisely: start at `t = 2/|lambda_2|`, end at `t = 5/gap`.

**Complexity:** MEDIUM. The automatic detection is the harder part; manual specification of the window is LOW complexity.

### Cross-Validation Against Exact Liouvillian

**Protocol:**
1. For n=4 and n=6 with identical configs (same beta, sigma, domain, delta, etc.):
   a. Run `run_lindbladian(jumps, liouv_config, hamiltonian)` to get exact `spectral_gap`
   b. Run `run_trajectories(jumps, therm_config, psi0, hamiltonian; observables=obs, save_every=k, ntraj=N, total_time=T)` to get time-resolved `<O>(t)`
   c. Fit exponential to each observable
   d. Compare `fitted_gap` vs `|Re(exact_gap)|`

2. Report:
   - Relative error: `|fitted_gap - exact_gap| / exact_gap`
   - The exact gap for reference
   - Confidence interval on fitted gap
   - Which observables give the best (closest to exact) gap estimate

**Expected precision:** With N_traj=10000 and T=5/gap, the fitted gap should agree with the exact gap to within ~5-10% for n=4 and ~10-20% for n=6 (more noise at larger dimensions due to higher-dimensional Hilbert space). This precision is sufficient for the paper's purpose (scaling analysis, not precision spectroscopy).

**Matching configs:** The `LiouvConfig` (for Liouvillian) and `ThermalizeConfig` (for trajectories) share all parameters except `mixing_time` and `delta`. Construct both from the same physical parameters to ensure consistency.

**Complexity:** MEDIUM. Mostly scripting -- the hard parts (trajectory simulation, Liouvillian construction) are already implemented.

### The `estimate_spectral_gap` Function

**Signature:**
```julia
function estimate_spectral_gap(
    times::Vector{Float64},
    measurements::Matrix{Float64};   # n_obs x n_times
    observable_names::Vector{String} = String[],
    gibbs_values::Union{Nothing, Vector{Float64}} = nothing,
    weights::Union{Nothing, Matrix{Float64}} = nothing,
    t_start_frac::Float64 = 0.1,
    model::Symbol = :exponential,  # :exponential or :damped_oscillation
) -> SpectralGapEstimate
```

**Return type:**
```julia
struct SpectralGapEstimate
    gap::Float64                          # best-estimate spectral gap
    gap_confidence::Tuple{Float64, Float64}  # 95% CI
    gap_per_observable::Vector{Float64}   # gap fitted per each observable
    best_observable_idx::Int              # which observable gave the best fit
    fit_results::Vector{Any}              # raw LsqFitResult per observable
    r_squared::Vector{Float64}           # R^2 per observable
    model_used::Symbol                    # :exponential or :damped_oscillation
end
```

**Logic:**
1. Fit the gap independently from each observable.
2. Rank by R-squared (best fit quality).
3. Report the gap from the observable with the best R-squared as the primary estimate.
4. Report the gap from all observables for consistency checking.
5. If gaps from multiple observables agree within their confidence intervals, report the weighted average as the final estimate.

**Complexity:** MEDIUM.

---

## Feature Dependencies

```
[build_spectral_gap_observables]  <-- adds M_z, extends build_convergence_observables
    |
    v
[Run trajectories with time-resolved observables]  <-- uses existing run_trajectories
    |
    v
[Single-exponential fit with LsqFit.jl]  <-- new dependency
    |           |
    v           v
[Fitting window selection]    [Fit quality metrics]
    |           |
    v           v
[estimate_spectral_gap function]  <-- combines all fitting
    |
    v
[Cross-validation vs exact Liouvillian]  <-- n=4,6 validation script
    |
    v
[Multi-observable consistency check]  (differentiator)
    |
    v
[Gap scaling plots: gap vs beta, gap vs n]  (differentiator)
```

### Dependency Notes

- **LsqFit.jl is the gate:** All fitting depends on this new dependency. Add it early.
- **Observable builder is independent of fitting:** Can build `build_spectral_gap_observables` immediately.
- **Cross-validation requires both paths:** Need both `run_lindbladian` (exact) and trajectory fitting to compare. These are independent and can be computed in parallel.
- **Model selection (exponential vs damped oscillation) is optional:** Start with pure exponential. Add damped oscillation only if the pure exponential gives poor fits.

---

## MVP Recommendation

### Must Complete

1. **Total magnetization observable** -- extend `build_convergence_observables` to include `M_z` (and rename/create `build_spectral_gap_observables`)
2. **Add LsqFit.jl dependency** -- needed for all fitting
3. **Single-exponential fit function** -- `A * exp(-gap * t) + C` with auto initial guess and parameter bounds
4. **Fitting window selection** -- at minimum, manual start fraction; ideally, simple auto-detection
5. **Fit quality metrics** -- R-squared, confidence interval on gap
6. **`estimate_spectral_gap` function** -- the reusable API that fits all observables and picks the best
7. **Cross-validation script** -- n=4 and n=6: exact gap vs fitted gap comparison

### Defer

- **Damped oscillation model** -- add only if pure exponential gives poor R-squared for some configs
- **Multi-exponential fit** -- ill-conditioned, only needed if early-time fitting is important
- **Bootstrap confidence intervals** -- linearized LsqFit CIs are sufficient for the paper
- **Automatic model selection (AIC/BIC)** -- overkill for single vs damped oscillation choice; visual inspection suffices
- **Observable overlap diagnostic** -- informative but not needed for the gap estimate itself

---

## Feature Prioritization Matrix

| Feature | Paper Value | Implementation Cost | Priority |
|---------|------------|---------------------|----------|
| Total magnetization observable | HIGH | LOW | P0 |
| LsqFit.jl dependency | CRITICAL | LOW | P0 |
| Single-exponential fit | CRITICAL | MEDIUM | P1 |
| Auto initial guess | CRITICAL | LOW | P1 |
| Fit quality metrics | CRITICAL | LOW | P1 |
| Fitting window selection (manual) | HIGH | LOW | P1 |
| estimate_spectral_gap function | CRITICAL | MEDIUM | P1 |
| Cross-validation n=4,6 | CRITICAL | MEDIUM | P2 |
| Multi-observable consistency | HIGH | LOW | P2 |
| Variance-weighted fitting | MEDIUM | LOW | P2 |
| Gap vs beta/n scaling plots | MEDIUM | LOW | P2 |
| Damped oscillation model | LOW-MEDIUM | MEDIUM | P3 |
| Fitting window auto-detection | LOW | MEDIUM | P3 |
| Observable overlap diagnostic | LOW | LOW | P3 |
| Bootstrap confidence intervals | LOW | MEDIUM | P3 |
| Multi-exponential fit | LOW | HIGH | P4 (defer) |

**Priority key:**
- P0: Enables all other work; do first
- P1: Core milestone features
- P2: Validates and strengthens the method
- P3: Nice to have; defer if time-constrained
- P4: Defer to future milestone

---

## Physics Reference: Observable Coupling to Gap Mode

### Which Observables See the Spectral Gap?

An observable O sees the spectral gap only if `tr(O * R_gap) != 0`, where `R_gap` is the gap eigenmatrix of the Lindbladian. For the 1D Heisenberg chain with disorder:

| Observable | Expected coupling to gap mode | Reasoning |
|------------|-------------------------------|-----------|
| Energy `<H>` | **STRONG** | `H` is diagonal in the eigenbasis. The gap mode is a density-matrix perturbation in the eigenbasis that describes the slowest redistribution of population among energy levels. The gap mode is "energy-like" -- it has large diagonal elements in the eigenbasis. `tr(H * R_gap)` is generically nonzero and large. |
| ZZ correlation `<Z_iZ_{i+1}>` | **STRONG** | The ZZ correlations probe two-body spin-spin interactions, which are the building blocks of the Hamiltonian. The gap mode, being related to the slowest energy redistribution, has significant structure in these two-body sectors. |
| Total magnetization `<M_z>` | **MODERATE to WEAK** | `M_z` commutes with the Heisenberg Hamiltonian's total spin (when disorder is absent). With disorder, this symmetry is broken, giving `M_z` some coupling. But the gap mode is primarily an energy excitation, not a magnetization excitation, so the coupling is weaker. |
| Single-site `<Z_i>` | **MODERATE** | Site-resolved magnetization has nontrivial coupling via the disorder field, but is a coarser probe than ZZ correlations. |

### Expected Spectral Gap Values

For the 1D Heisenberg chain simulated in QuantumFurnace (rescaled spectrum in [0, 0.45]):

| n | beta | Expected gap (rescaled units) | Mixing time ~1/gap | Notes |
|---|------|------------------------------|-------------------|-------|
| 4 | 5 | ~0.01-0.1 | ~10-100 | Warm, fast mixing |
| 4 | 10 | ~0.005-0.05 | ~20-200 | Moderate |
| 4 | 20 | ~0.001-0.01 | ~100-1000 | Cold, slow mixing |
| 6 | 10 | ~0.003-0.03 | ~30-300 | Larger system |
| 8 | 10 | ~0.001-0.01 | ~100-1000 | Trajectory-only regime |

Note: These are rough estimates. The actual gap depends on the specific disordered Hamiltonian instance, the domain (Energy/Time/Trotter), and the Lindbladian parameters (sigma, transition function). Cross-validation at n=4,6 will calibrate expectations.

### Required Simulation Parameters for Gap Estimation

To resolve the spectral gap from trajectory data:

- **`total_time`:** Must be at least `5/gap` for the observable to reach its plateau. For gap ~ 0.01, total_time ~ 500 in rescaled units.
- **`save_every`:** Should give ~200-500 time points. With `delta = 0.001` and `total_time = 500`, num_steps = 500000. Set `save_every = 1000-2500` for 200-500 data points.
- **`ntraj`:** More trajectories reduce noise, improving fit quality. For gap estimation: ntraj = 5000-20000 is a reasonable range. The fitted gap precision scales as ~1/sqrt(ntraj).
- **`delta`:** Must be small enough for simulation accuracy (existing tests confirm delta = 0.001 works well). No additional constraint from gap estimation.

---

## Sources

### Verified (HIGH confidence)
- QuantumFurnace.jl codebase -- Direct analysis of `trajectories.jl` (run_trajectories, _run_chunk_with_obs!, _accumulate_measurements!), `convergence.jl` (build_convergence_observables, run_trajectories_convergence), `furnace.jl` (run_lindbladian, LindbladianResult.spectral_gap), `structs.jl` (TrajectoryResult, LindbladianResult, ConvergenceData)
- [LsqFit.jl Documentation](https://julianlsolvers.github.io/LsqFit.jl/latest/tutorial/) -- curve_fit API, confidence_interval, estimate_covar, Levenberg-Marquardt algorithm
- [LsqFit.jl GitHub](https://github.com/JuliaNLSolvers/LsqFit.jl) -- MIT license, pure Julia, box constraints support

### Verified (MEDIUM confidence)
- [Spectral Gap and Exponential Decay of Correlations](https://arxiv.org/abs/math-ph/0507008) (Nachtergaele & Sims 2006) -- Mathematical foundation: spectral gap implies exponential decay of correlations
- [Fast Mixing of Quantum Spin Chains at All Temperatures](https://arxiv.org/html/2510.08533) -- System-size independent spectral gap for 1D chains at any finite temperature
- [Mixing Time of Open Quantum Systems via Hypocoercivity](https://arxiv.org/abs/2404.11503) -- Relationship between spectral gap, mixing time, and observable autocorrelation
- [Computationally Efficient Estimation of the Spectral Gap of a Markov Chain](https://arxiv.org/abs/1806.06047) -- UCPI algorithm for gap estimation from sample paths (classical analog)
- [Spectral Gap Estimation via Adiabatic Preparation](https://arxiv.org/html/2512.19288) -- Quantum device approach (not directly applicable but contextually relevant)
- [Universal Predictors for Mixing Time more than Liouvillian Gap](https://arxiv.org/html/2601.06256) -- Recent work on mixing time predictors beyond the gap

### Domain knowledge (HIGH confidence, established physics/mathematics)
- Lindbladian semigroup spectral decomposition: `rho(t) = rho_ss + sum c_i R_i exp(lambda_i t)` is standard open quantum systems theory
- Observable decay rate = real part of dominant non-zero eigenvalue (from spectral decomposition of `tr(O rho(t))`)
- Nonlinear least squares with Levenberg-Marquardt is the standard approach for exponential decay fitting in physics
- Heisenberg chain gap mode coupling to energy and spin-spin observables -- established condensed matter physics
- KMS-detailed-balanced Lindbladians have unique steady state and real spectral gap for the symmetrized generator

---
*Feature research for: v1.3 Spectral Gap Estimation from Trajectory Observables*
*Researched: 2026-02-16*
