# Feature Landscape: v1.4 Spectral Gap Refinement Diagnostics

**Domain:** Diagnostic infrastructure and improved estimators for Lindbladian spectral gap estimation
**Researched:** 2026-02-19
**Confidence:** HIGH (features specified by reference documents, existing codebase well-understood, numerical methods standard)

## Scope

This research covers features for **v1.4 Spectral Gap Refinement**: a comprehensive diagnostic and improved estimation pipeline that addresses the 7 catalogued error sources in trajectory-based spectral gap estimation. The v1.3 milestone delivered single-exponential fitting achieving 0.72% accuracy at n=4 but diagnosed zero-overlap physics limitations at n=6 and confirmed (Quick-32) that non-monotonic gap estimation is a fitting artifact, not a simulation bug. The v1.4 pipeline attacks these limitations head-on.

**Key constraint from reference documents:** All trajectory runs use `with_coherent=true`, 4 threads, Trotter steps=10 for OFT/B, delta is the Trotter parameter.

---

## Existing Features (Already Built in v1.3)

These are NOT part of v1.4 scope but are the foundation:

| Feature | Location | Status |
|---------|----------|--------|
| `estimate_spectral_gap` API | `gap_estimation.jl` | Shipped, single-exponential |
| `eigenbasis_overlap_analysis` | `gap_estimation.jl` | Shipped, 8 observables |
| `run_observable_trajectories` | `trajectories.jl` | Shipped, obs-only path |
| `fit_exponential_decay` | `fitting.jl` | Shipped, `A*exp(-gap*t)+C` |
| `SpectralGapResult`, `FitResult` | `gap_estimation.jl`, `fitting.jl` | Shipped |
| `build_preset_trajectory_observables` | `convergence.jl` | Shipped, 8 observables |
| Full Lindbladian construction + Arpack | `furnace.jl` | Shipped, shift-invert |
| Multi-threaded trajectory engine | `trajectories.jl` | Shipped, 4 threads |

---

## Table Stakes

Features that ARE the milestone. Without these, v1.4 delivers no value over v1.3.

| Feature | Why Expected | Complexity | Dependencies |
|---------|--------------|------------|--------------|
| **Effective rate plot `lambda_eff(t)`** | The single most useful model-free diagnostic. Reveals the three time regimes (transient, plateau, noise floor) without any fitting assumptions. Directly identifies the golden fitting window. This is the "Rosetta Stone" of the error diagnosis. | MEDIUM | Trajectory-averaged observable time series, steady-state subtraction |
| **Two-exponential fit with Prony initialization** | The root cause finding from Quick-32: single-exponential fitting is the bottleneck. Two-exponential `c1*exp(-g1*t) + c2*exp(-g2*t)` absorbs excited-state contamination into the second term, analogous to lattice QCD multi-state analysis. Prony two-point initialization provides robust starting values. | HIGH | `LsqFit.jl` (existing), effective rate plateau for init, bounded constraints |
| **Bootstrap resampling for trajectory error bars** | Linearized LsqFit confidence intervals assume Gaussian errors and local linearity -- poor assumptions for noisy exponential fits near the noise floor. Bootstrap over trajectories provides nonparametric error bars that correctly propagate the trajectory-level variance. | MEDIUM | Per-trajectory observable storage (new), resampling loop |
| **Automatic fitting window selection** | Manual `skip_initial=0.1` is fragile. t_max via SNR threshold (SNR > 3), t_min via gamma_1 stability test. The effective rate plot provides visual validation. | MEDIUM | Effective rate plot, trajectory variance estimates |
| **Richardson extrapolation for Trotter bias** | Quick-30 found Richardson was ineffective with single-exponential fitting. With improved two-exponential fitting, Richardson at two delta values eliminates the O(delta) Trotter bias, giving O(delta^2) residual. | LOW | Gap estimates at 2+ delta values, simple linear algebra |
| **Anti-Hermitian defect computation** | KMS similarity transform `D(rho,L) = rho^{-1/4} L[rho^{1/4}(.)rho^{1/4}] rho^{-1/4}`, decomposed as D = H + A. The ratio `||A|| / lambda_gap(H)` determines whether real-exponential fitting is valid or damped-oscillation model needed. | MEDIUM | Gibbs state eigendecomposition (existing), vectorized superoperator construction |
| **Delta-Sz symmetry sector labeling** | Label Lindbladian eigenvectors by the quantum number `Delta_Sz = Sz(E_i) - Sz(E_j)` they carry. Explains why certain observables see zero gap-mode overlap (Quick-25 through Quick-27 diagnosed this). | LOW | Exact eigenvectors from `run_lindbladian` (existing), Hamiltonian eigenbasis Sz values |

---

## Feature Details

### 1. Effective Rate Plot: `lambda_eff(t)`

**What it computes:**
```
Delta(t) = <O>(t) - <O>_ss
lambda_eff(t) = -(1/tau) * ln|Delta(t+tau) / Delta(t)|
```
where `tau = lag * dt` (default `lag=3` for noise smoothing).

**Why it works:** For a pure single-exponential `Delta(t) = c * exp(-lambda * t)`, the effective rate is exactly `lambda` for all t. For multi-exponential signals, `lambda_eff(t)` starts high (dominated by fast modes), decreases to a plateau at the true gap, then becomes erratic when the signal hits the noise floor. The plateau region IS the golden fitting window.

**What the plot shows:**
- Horizontal dashed line at exact `lambda_gap` (n=4,6 reference)
- Error bands from propagated trajectory variance: `sigma_lambda(t) ~ (1/tau) * sqrt(Var[Delta(t)]/|Delta(t)|^2 + Var[Delta(t+tau)]/|Delta(t+tau)|^2)`
- Vertical dashed lines at t_min (plateau entry) and t_max (noise floor)

**Implementation notes:**
- Sign-change guard: when `Delta(t)` and `Delta(t+tau)` have different signs, output NaN (signal has crossed zero due to noise or oscillation)
- Lag parameter: start with `lag=3`, adjustable. Larger lag smooths but loses time resolution.
- The reference document provides a complete Julia implementation snippet.

**Complexity:** MEDIUM. The computation itself is simple. The complexity lies in the error propagation, sensible defaults for lag, and the overlay plot generation.

**Dependencies:** Trajectory-averaged observable data (existing), Gibbs steady-state expectation values (existing via `_compute_gibbs_observable_values`).

---

### 2. Two-Exponential Fit with Prony Initialization

**Why single-exponential fails (from Quick-32 root cause):**
The signal is `Delta(t) = c1*exp(-g1*t) + c2*exp(-g2*t) + ...`. A single-exponential fit `A*exp(-gap*t)+C` finds a compromise rate between the true gap g1 and the faster modes g2, g3, ..., systematically overestimating the gap. Different delta values excite different spectral modes with different weights, explaining the non-monotonic delta-scaling found in Quick-30/31.

**Model:**
```
Delta(t) = c1 * exp(-g1 * t) + c2 * exp(-g2 * t)
```
with constraint `0 < g1 < g2`. The gap estimate is g1; the second term is a nuisance parameter that absorbs fast-mode contamination.

**Prony two-point initialization:**
Pick two time points t_a, t_b = t_a + tau in the signal (in the effective-rate plateau transition region). Form:
```
r1 = Delta(t_a + tau/2) / Delta(t_a)
r2 = Delta(t_b) / Delta(t_a + tau/2)
```
Solve the quadratic `z^2 - (r1+r2)*z + r1*r2 = 0` to get `z1, z2`. Then `g1 = -2*ln(z1)/tau`, `g2 = -2*ln(z2)/tau`. This gives initial guesses without iterative optimization.

**Alternative initialization from effective rate plot:**
```
g1_init = plateau value of lambda_eff(t)
g2_init = lambda_eff at early time (t ~ 0)
c1_init = Delta(t_mid) / exp(-g1_init * t_mid)
c2_init = Delta(0) - c1_init
```

**Bounded constraints (critical for convergence):**
```
g1 in [0.01, 3 * g1_init]
g2 in [g1, 10 * g2_init]
|c1|, |c2| < 10 * |Delta(0)|
```

**Validation:** At n=4,6, extracted g1 must agree with exact `lambda_gap` within 2-sigma bootstrap error bars.

**Complexity:** HIGH. Two-exponential fitting is notoriously ill-conditioned (nearly degenerate eigenvalues produce flat cost landscape). Requires both initialization strategies (Prony + effective-rate fallback), careful bounds, and bootstrap validation of convergence.

**Dependencies:** Effective rate plot (for initialization and window), LsqFit.jl (existing), bootstrap (for error bars).

---

### 3. Bootstrap Resampling for Error Bars

**Why needed:** The current `FitResult.gap_se` and `gap_ci` come from LsqFit's linearized covariance estimate `J^T J` inversion. This assumes: (a) Gaussian errors, (b) the Jacobian is approximately constant over the confidence region, (c) residuals are independent. None hold well for noisy exponential decay: the noise is not Gaussian (pure-state measurement outcomes), the model is nonlinear near boundaries, and consecutive time points from the same trajectory are correlated.

**Procedure:**
1. Store per-trajectory observable time series `O_i(t)` for i=1..N_traj (NEW: currently only the mean is stored)
2. For b = 1..N_boot (default 200):
   a. Draw N_traj trajectory indices with replacement
   b. Compute bootstrap-averaged `<O>^(b)(t)`
   c. Subtract steady-state value
   d. Compute `lambda_eff^(b)(t)` and/or run two-exponential fit
   e. Record g1^(b)
3. Report: mean, median, standard deviation, and [2.5%, 97.5%] percentiles of {g1^(b)}

**Key implementation decision:** Storing per-trajectory data. Currently `_accumulate_measurements!` sums in-place. For bootstrap, we need either:
- (a) Store all N_traj x n_saves x n_obs data in memory (for 20k trajectories, 500 saves, 8 observables: ~64 MB Float64 -- feasible)
- (b) Use block bootstrap: group trajectories into ~50 blocks, resample blocks (reduces storage to block-level means)

**Recommendation:** Option (a) for n=4,6 validation (affordable memory). Option (b) for larger n when N_traj is very large. Provide both.

**Complexity:** MEDIUM. The resampling loop is straightforward. The main complexity is in the new storage path for per-trajectory data.

**Dependencies:** Per-trajectory observable storage (new functionality in trajectories.jl), fitting functions.

---

### 4. Automatic Fitting Window Selection

**t_max selection (SNR threshold):**
```
SNR(t) = |Delta(t)| / sigma_Delta(t)
t_max = last time point where SNR(t) > 3
```
where `sigma_Delta(t) = sqrt(Var_traj[<psi|O|psi>] / N_traj)` is the standard error of the mean.

**t_min selection (stability test):**
Run two-exponential fit over [t_min, t_max] for a sequence of t_min values:
```
t_min = 0, dt, 2*dt, 3*dt, ...
```
Plot g1 vs t_min. Choose t_min as the smallest value where g1 has stabilized (does not change beyond its bootstrap error bar as t_min increases further).

**Why this is important:**
- Manual `skip_initial=0.1` is system-dependent. For fast-mixing systems, 10% skip is too much (wastes data). For slow-mixing, it is too little (transient still present).
- The t_min stability plot is itself a valuable diagnostic for the thesis: it demonstrates that the result is robust to fitting window choice.

**Complexity:** MEDIUM. t_max via SNR is simple. t_min stability requires running the fitter in a loop with different windows, which is computationally manageable but needs careful edge-case handling (when no stable plateau exists, when the window becomes too small).

**Dependencies:** Trajectory variance estimates, two-exponential fit, effective rate plot for validation.

---

### 5. Richardson Extrapolation for Trotter Bias

**The formula:**
```
gap_rich = 2 * gap(delta/2) - gap(delta)
```
with error `sigma_rich = sqrt(4 * sigma^2(delta/2) + sigma^2(delta))`.

**Why Quick-30 failed:** Quick-30 applied Richardson to single-exponential fits. Since single-exponential fitting introduces a fitting-model bias that varies non-monotonically with delta (different delta excites different spectral modes), Richardson amplified this fitting artifact rather than correcting Trotter bias.

**Why it should work with two-exponential fitting:** The two-exponential fit correctly isolates g1 (the true gap) from the fast-mode contamination, regardless of how delta affects the mode weights. The remaining bias IS the O(delta) Trotter error, which Richardson cancels.

**Required data:** Trajectory runs at 2+ delta values (e.g., delta=0.01 and delta=0.005) with identical everything else (same N_traj, same T, same observables, same psi0, same seed).

**Complexity:** LOW. The formula is trivial. The cost is in running trajectories at multiple delta values (but delta-convergence runs are already planned for diagnostics).

**Dependencies:** Two-exponential fit (must produce reliable per-delta gap estimates), bootstrap error bars.

---

### 6. Anti-Hermitian Defect Computation

**The KMS similarity transform:**
```
D(rho_beta, L) = rho_beta^{-1/4} * L[rho_beta^{1/4} (.) rho_beta^{1/4}] * rho_beta^{-1/4}
```
This is a superoperator. Vectorized, it becomes a d^2 x d^2 matrix D. Decompose:
```
H = (D + D^dagger) / 2     (Hermitian part)
A = (D - D^dagger) / 2     (anti-Hermitian part)
```

**What to report:**
1. `||A||` (operator norm of anti-Hermitian part)
2. `lambda_gap(H)` (spectral gap of the Hermitian part)
3. Ratio `||A|| / lambda_gap(H)` -- if < 0.01, real-exponential model is safe; if > 0.1, need damped-oscillation model
4. `|Im(lambda_2)| / |Re(lambda_2)|` for the leading eigenvalues -- direct check for oscillatory modes

**Implementation notes:**
- Computing `rho_beta^{1/4}`: Since `rho_beta` is diagonal in the energy eigenbasis, raise diagonal entries to the 1/4 power. For TrotterDomain, `rho_beta` is in the Trotter eigenbasis -- diagonalize it first, raise eigenvalues to 1/4, transform back.
- The vectorized D matrix is constructed by applying the map to each basis element of the d^2-dimensional vectorized operator space (same technique as `construct_lindbladian`).
- Size constraint: d^2 x d^2 is 4096x4096 for n=6 -- feasible.

**Complexity:** MEDIUM. The matrix fourth root is straightforward via eigendecomposition. The superoperator construction reuses the existing vectorization machinery from `qi_tools.jl` (_kron!, _vectorize_liouv_diss_and_add!). The main work is assembling the similarity-transformed superoperator.

**Dependencies:** Exact Lindbladian (existing), Gibbs state eigendecomposition (existing), vectorization utilities (existing).

---

### 7. Delta-Sz Symmetry Sector Labeling

**What to compute:**
For each energy eigenstate |E_i> of H, compute `Sz_tot(E_i) = <E_i| sum_j Z_j |E_i>`. For each Lindbladian eigenvector R_k (reshaped as a d x d matrix), determine the dominant `Delta_Sz = Sz(E_i) - Sz(E_j)` quantum number by computing:
```
weight(Delta_Sz) = sum_{i,j: Sz(i)-Sz(j)=Delta_Sz} |R_k[i,j]|^2
```
Label R_k by the Delta_Sz with the largest weight.

**Why this matters:**
Quick-25 through Quick-27 diagnosed that the n=6 gap mode lives in a symmetry sector that zero-overlap observables cannot access. Delta-Sz labeling makes this diagnosis systematic:
- If lambda_2 carries Delta_Sz = 0: it describes population decay. Diagonal observables (H, Z_i) can see it.
- If lambda_2 carries Delta_Sz != 0: it describes coherence decay. Only off-diagonal observables or symmetry-breaking observables (XZ_stagg) can see it.
- Explains why Quick-28 found that disorder breaks the symmetry protection.

**Complexity:** LOW. This is a post-processing analysis on existing eigendecomposition data. No new simulation infrastructure needed.

**Dependencies:** Exact Lindbladian eigenvectors (existing from `run_lindbladian`), Hamiltonian eigenbasis (existing).

---

## Differentiators

Features that strengthen the milestone but are not strictly required for the core diagnostic pipeline.

| Feature | Value Proposition | Complexity | Dependencies |
|---------|-------------------|------------|--------------|
| **Summary dashboard (7-panel figure)** | Single thesis-quality figure showing all diagnostics at once: spectrum, defect metrics, overlap coefficients, effective rate, delta-convergence, two-exp fit, t_min stability. Communicates the entire error analysis in one visual. | MEDIUM | All 7 table-stakes features |
| **External field comparison (h = 0.1J)** | Run full diagnostic with and without symmetry-breaking field. If the gap estimate changes dramatically while the exact gap changes mildly, symmetry sector restriction is confirmed as the dominant error source. | LOW | All diagnostics running, field parameter in Hamiltonian |
| **Multi-observable minimum-gap selector** | Extend the existing `_select_best_observable` to use the two-exponential fit and report the gap as `min_j(g1_j)` across all observables. Observables can only overestimate the gap. | LOW | Two-exponential fit per observable |
| **Damped-oscillation fit model** | `c * exp(-gamma*t) * cos(omega*t + phi)` for cases where anti-Hermitian defect is significant. Auto-selected when `|Im(lambda_2)/Re(lambda_2)| > 0.1`. | MEDIUM | Anti-Hermitian defect computation |
| **Observable overlap with left/right eigenvectors** | Extend `eigenbasis_overlap_analysis` to compute `c_k = Tr[O R_k] * Tr[L_k^dagger (rho0 - rho_beta)]` using both left and right eigenvectors (currently only uses right). More accurate coefficient decomposition. | LOW | Existing overlap analysis, left eigenvector computation |
| **Store leading 20-30 eigenvalues (not just 2)** | Currently `run_lindbladian` extracts only nev=2 eigenvalues. Expand to 20-30 for spectrum visualization and overlap coefficient analysis. | LOW | Modify Arpack `nev` parameter |

---

## Anti-Features

Features to explicitly NOT build in this milestone.

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| **Extend to n=8 via sparse Lindbladian** | The reference document explicitly defers this: "DEFERRED TO A LATER MILESTONE." n=8 requires 65536x65536 sparse Lindbladian, KrylovKit.jl methods, and a fast `apply_lindbladian!` function. Too much scope. | Validate at n=4,6 only. n=8 is a separate milestone. |
| **Matrix pencil / ESPRIT methods for spectral gap** | These are signal-processing alternatives to nonlinear least squares. While potentially more robust for multi-exponential extraction, they require uniform sampling (which we have) but add a new algorithmic paradigm with its own tuning parameters. The two-exponential fit with Prony init is sufficient. | Two-exponential fit with Prony initialization covers the same ground with less implementation complexity. |
| **Full GEVP (Generalized Eigenvalue Problem) approach** | Used in lattice QCD for multi-hadron spectroscopy. Requires building a matrix of correlation functions from multiple operators. Powerful but massively over-engineered for our 1-2 exponential extraction problem. | Two-exponential fit is the right level of sophistication. |
| **GPU-accelerated trajectory sampling** | The bottleneck for this milestone is understanding, not throughput. We need 20k trajectories at n=4,6 -- minutes on 4 CPU threads. GPU would add huge complexity for zero benefit at this system size. | Keep CPU-only. GPU is a future performance milestone. |
| **Adaptive delta selection** | Automatically choosing the optimal delta to balance Trotter bias vs. computational cost. This is an optimization problem that presumes the diagnostic pipeline already works. | Run at 2-3 fixed delta values. Use Richardson extrapolation. |
| **Real-time interactive plotting** | The reference document specifies static diagnostic figures for the thesis. Interactive plotting adds complexity (Makie.jl event handling, live updates) without scientific value. | Generate static plots post-hoc. Use CairoMakie for publication quality. |
| **Per-trajectory density matrix reconstruction during bootstrap** | Storing the full density matrix per trajectory would require N_traj x d^2 complex matrices. For 20k trajectories at n=6, that is 20000 x 4096 x 16 bytes = 1.3 GB. Unnecessary -- observable time series suffice. | Store per-trajectory observable values only (64 MB for 20k x 500 x 8). |
| **Weighted least squares with variance weights** | The reference document's effective rate approach and bootstrap make variance-weighted fitting less critical. The error propagation through lambda_eff already accounts for varying noise levels. Adding weights to the two-exponential fit complicates the already ill-conditioned optimization. | Use unweighted two-exponential fit with SNR-based window truncation instead. |

---

## Feature Dependencies

```
[Anti-Hermitian Defect Computation]          [Delta-Sz Symmetry Labels]
         |                                              |
         v                                              v
    Validates real-exp model              Explains zero-overlap at n=6
         |                                              |
         +-------> [Effective Rate Plot] <--------------+
                           |
                           v
              Identifies golden window
                    |              |
                    v              v
     [Auto Fitting Window]   [Two-Exponential Fit]
          |                        |
          +--------+   +---------+
                   v   v
          [Bootstrap Error Bars]
                    |
                    v
        [Richardson Extrapolation]  (requires gap at 2+ delta)
                    |
                    v
          [Validated Gap Estimate]
                    |
                    v
          [Summary Dashboard]
```

### Dependency Notes

- **Effective rate plot is the keystone:** It feeds initialization for the two-exponential fit, defines the fitting window, and provides the model-free reference against which all fits are validated. Build it first.
- **Anti-Hermitian defect and symmetry labels are independent diagnostics:** They can be computed in parallel with no dependency on each other or on the fitting pipeline. They answer "is the real-exponential model valid?" and "why does this observable miss the gap?" respectively.
- **Bootstrap requires per-trajectory storage:** This is a new data path in `trajectories.jl`. The current `_accumulate_measurements!` only stores the running sum. A new accumulation mode or return type is needed.
- **Richardson extrapolation depends on the two-exponential fit being reliable:** Quick-30 showed that Richardson amplifies fitting artifacts. Only apply Richardson after the two-exponential fit is validated.
- **Auto fitting window depends on effective rate AND two-exponential fit:** t_max is from SNR (effective rate infrastructure), t_min is from fit stability (requires running the fitter).
- **Summary dashboard depends on everything:** It is the last feature, composing all results into one figure.

---

## MVP Definition

### Launch With (v1.4 core)

- [x] **Effective rate plot** -- the model-free diagnostic that reveals everything
- [x] **Two-exponential fit with Prony init** -- the root cause fix from Quick-32
- [x] **Bootstrap error bars** -- trustworthy uncertainty quantification
- [x] **Auto fitting window (SNR + stability)** -- removes the manual skip_initial fragility
- [x] **Anti-Hermitian defect** -- validates the real-exponential model assumption
- [x] **Delta-Sz symmetry labels** -- explains the n=6 zero-overlap mystery
- [x] **Richardson extrapolation** -- eliminates O(delta) Trotter bias

### Add After Core Validation

- [ ] **Summary dashboard** -- after all diagnostics work, compose the thesis figure
- [ ] **External field comparison** -- test symmetry-breaking effect on gap estimation
- [ ] **Multi-observable minimum-gap selector** -- improve `_select_best_observable`
- [ ] **Store 20-30 eigenvalues** -- richer spectrum visualization

### Future Consideration (v1.5+)

- [ ] **n=8 sparse Lindbladian** -- separate milestone per reference document
- [ ] **Damped-oscillation fit model** -- only if anti-Hermitian defect is significant
- [ ] **GEVP / matrix pencil methods** -- only if two-exponential fit proves insufficient

---

## Feature Prioritization Matrix

| Feature | Scientific Value | Implementation Cost | Priority |
|---------|-----------------|---------------------|----------|
| Effective rate plot lambda_eff(t) | CRITICAL | MEDIUM | P0 |
| Two-exponential fit + Prony init | CRITICAL | HIGH | P0 |
| Per-trajectory observable storage | CRITICAL (enables bootstrap) | MEDIUM | P0 |
| Bootstrap error bars | HIGH | MEDIUM | P1 |
| Auto t_max (SNR threshold) | HIGH | LOW | P1 |
| Auto t_min (stability test) | HIGH | MEDIUM | P1 |
| Anti-Hermitian defect computation | HIGH | MEDIUM | P1 |
| Delta-Sz symmetry sector labels | HIGH | LOW | P1 |
| Richardson extrapolation | HIGH | LOW | P1 |
| Store 20-30 Lindbladian eigenvalues | MEDIUM | LOW | P2 |
| Summary dashboard (7-panel figure) | HIGH (thesis) | MEDIUM | P2 |
| External field comparison | MEDIUM | LOW | P2 |
| Multi-observable min-gap selector | MEDIUM | LOW | P2 |
| Left/right eigenvector overlap | LOW | LOW | P3 |
| Damped-oscillation fit model | LOW | MEDIUM | P3 |

**Priority key:**
- P0: Foundation that enables other features; build first
- P1: Core diagnostic pipeline
- P2: Validation, visualization, thesis figures
- P3: Only if needed based on diagnostic results

---

## Phase Structure Implications

Based on feature dependencies, the natural phase ordering is:

**Phase A (Foundation):** Per-trajectory storage, effective rate plot, anti-Hermitian defect, symmetry labels, expand eigenvalue extraction to 20-30. These are independent and can be built/validated separately.

**Phase B (Improved Estimator):** Two-exponential fit with Prony init, auto fitting window (SNR + stability), bootstrap error bars. These depend on Phase A outputs for initialization and validation.

**Phase C (Bias Elimination):** Richardson extrapolation at 2+ delta values, external field comparison. Depends on Phase B producing reliable per-delta gap estimates.

**Phase D (Synthesis):** Summary dashboard, final validated gap estimate at n=4,6, comparison table (exact vs estimated, with sigma discrepancy).

---

## Technical Notes for Implementation

### Per-Trajectory Storage Architecture

The current code in `_run_chunk_obs_only!` sums measurements in-place via `_accumulate_measurements!`. For bootstrap, we need one of:
1. **Full storage:** New `Matrix{Float64}` of size `(n_traj, n_obs, n_saves)`. Return alongside the mean.
2. **Streaming variance:** Welford's online algorithm to accumulate mean AND variance in one pass (no per-trajectory storage, but variance is needed for SNR calculation).
3. **Block storage:** Divide trajectories into B blocks of size N_traj/B. Store B block means. Bootstrap resamples blocks.

**Recommendation:** Use approach (1) for n=4,6 (small enough) and approach (3) as a fallback for larger systems.

### LsqFit.jl Two-Exponential Model

```julia
_two_exp_model(t, p) = @. p[1] * exp(-p[2] * t) + p[3] * exp(-p[4] * t)
# p = [c1, g1, c2, g2] with constraint g1 < g2
# Enforce via bounds: lower = [-Inf, 0.01, -Inf, 0.01], upper = [Inf, g1_ub, Inf, g2_ub]
```

The constraint `g1 < g2` is enforced indirectly: set `upper[2] = initial_g2` and `lower[4] = initial_g1` based on Prony initialization. If the fit swaps g1 and g2, sort post-hoc.

### Prony Two-Point Implementation

```julia
function prony_init(Delta, times; t_a_idx, tau_steps)
    tau = (times[t_a_idx + tau_steps] - times[t_a_idx])
    r1 = Delta[t_a_idx + div(tau_steps,2)] / Delta[t_a_idx]
    r2 = Delta[t_a_idx + tau_steps] / Delta[t_a_idx + div(tau_steps,2)]
    # Solve z^2 - (r1+r2)*z + r1*r2 = 0
    disc = (r1 + r2)^2 - 4 * r1 * r2
    disc < 0 && return nothing  # complex roots: fall back to effective rate init
    z1 = ((r1 + r2) + sqrt(disc)) / 2
    z2 = ((r1 + r2) - sqrt(disc)) / 2
    g1 = -2 * log(max(z1, 1e-10)) / tau
    g2 = -2 * log(max(z2, 1e-10)) / tau
    g1, g2 = minmax(g1, g2)  # ensure g1 < g2
    return (g1, g2)
end
```

### Observable Update

The reference document specifies updating the observable set to:
- Z1 (single-site Z on first site)
- X1 (single-site X on first site)
- Z1 * Z_{n/2} (two-point correlator)
- H (energy)
- Random traceless Hermitian (control)

This differs from the current 8-observable set. The `build_preset_trajectory_observables` function should be updated accordingly, or a separate builder provided for the diagnostics.

---

## Sources

### HIGH Confidence (codebase + reference documents)
- QuantumFurnace.jl codebase: `gap_estimation.jl`, `fitting.jl`, `trajectories.jl`, `convergence.jl`, `furnace.jl`, `qi_tools.jl`
- `supplementary-informations/spectral-gap-refinements-instructions.md` -- 5-part diagnostic pipeline specification
- `supplementary-informations/error_catalogue_spectral_gap_estimation.md` -- 7 error sources catalogued
- Quick-30/31/32 summaries: trajectory correctness confirmed, fitting is the bottleneck

### MEDIUM Confidence (established methods, verified by multiple sources)
- [Prony's Method (Wikipedia)](https://en.wikipedia.org/wiki/Prony's_method) -- exponential extraction via polynomial root-finding
- [LsqFit.jl Tutorial](https://julianlsolvers.github.io/LsqFit.jl/latest/tutorial/) -- Julia nonlinear least squares with bounds
- [Excited State Systematics in Lattice QCD](https://arxiv.org/abs/2104.05226) -- multi-state analysis for excited state contamination (analog problem)
- [Monte Carlo Wave-Function Method](https://www.sciencedirect.com/science/article/abs/pii/S0010465518304314) -- robust adaptive MCWF with convergence study
- [Reducing Circuit Depth in Lindblad Simulation via Step-Size Extrapolation](https://arxiv.org/html/2507.22341) -- Richardson extrapolation for Lindbladian Trotter error
- [Exponentially Reduced Circuit Depths Using Trotter Error Mitigation](https://doi.org/10.1103/kw39-yxq5) -- rigorous Richardson extrapolation performance guarantees
- [Symmetry Classification of Many-Body Lindbladians](https://dx.doi.org/10.1103/PhysRevX.13.031019) -- symmetry sector labeling framework
- [Chen et al. 2023 (arXiv:2303.18224)](https://arxiv.org/pdf/2303.18224) -- KMS detailed balance, Proposition II.3, similarity transform
- [Efficient Quantum Gibbs Samplers with KMS Detailed Balance](https://link.springer.com/article/10.1007/s00220-025-05235-3) -- theoretical framework for KMS Lindbladians
- [QuantumToolbox.jl Monte Carlo Solver](https://qutip.org/QuantumToolbox.jl/stable/users_guide/time_evolution/mcsolve) -- bootstrap/jackknife for trajectory error bars in QuTiP
- [Classical Shadow Nonparametric Bootstrap](https://arxiv.org/html/2511.09793) -- bootstrap methods for quantum state estimation

---
*Feature research for: v1.4 Spectral Gap Refinement Diagnostics*
*Researched: 2026-02-19*
