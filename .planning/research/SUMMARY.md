# Project Research Summary

**Project:** v1.3 Spectral Gap Estimation from Trajectory Observable Decay
**Domain:** Quantum trajectory simulation + nonlinear exponential decay fitting + statistical cross-validation
**Researched:** 2026-02-16
**Confidence:** HIGH

## Executive Summary

This milestone adds spectral gap estimation via exponential fitting to time-resolved observable decay from quantum trajectory simulations. The physics is well-established: for a Lindbladian with spectral gap λ, observable expectation values decay exponentially toward thermal equilibrium at rate λ. By running trajectory simulations with observable measurements at regular intervals, fitting `A * exp(-λ * t) + C` to the decay curve, and extracting λ, we estimate the spectral gap without constructing the full Liouvillian (which becomes infeasible at n≥8). Cross-validation against exact Liouvillian eigenvalues at n=4,6 establishes trust.

The recommended approach is minimal and well-constrained: add LsqFit.jl (one new dependency) for Levenberg-Marquardt curve fitting with parameter bounds and built-in confidence intervals, create a new `spectral_gap.jl` module that wraps existing trajectory primitives with exponential fitting logic, extend observable builders to include total magnetization, and implement cross-validation helpers. The existing trajectory infrastructure (`run_trajectories` with observables, `step_along_trajectory!`, measurement accumulation) is reused as-is. The architecture follows established patterns from `run_trajectories_convergence`.

The key risk is multi-exponential contamination: the observable decay is a sum of exponentials (one per Liouvillian eigenmode), not a pure single exponential. Fitting the full time series biases the gap estimate high. The mitigation is late-time fitting (skip the first 10-20% of data where fast modes dominate) combined with multi-observable consistency checks (fit energy, magnetization, and ZZ correlations independently; the smallest fitted rate with good quality is the true gap). Secondary risks include basis mismatch for new observables (critical pitfall from v1.2 quick-task-20 experience), noise floor at late times (requires sufficient trajectories: 5000+ for n=4,6, 10000+ for n=8), and confusing complex Liouvillian eigenvalues with real decay rates (must compare against `-real(spectral_gap)`, not `abs(spectral_gap)`).

## Key Findings

### Recommended Stack

LsqFit.jl (v0.15+) is the single new production dependency. It provides Levenberg-Marquardt nonlinear least squares with parameter bounds (essential: gap > 0), confidence intervals via t-distribution from Jacobian covariance, standard errors, and weighted fitting (for non-uniform trajectory noise). It is pure Julia, actively maintained by JuliaNLSolvers, and adds minimal dependency footprint: its transitive deps (ForwardDiff, NLSolversBase, StatsAPI) are already resolved via the existing Optim.jl dependency. The only truly new package is Distributions.jl (used internally by LsqFit for t-quantiles in confidence intervals), which is a standard, lightweight, well-maintained statistics package.

**Core technologies:**
- **LsqFit.jl v0.15+**: Nonlinear curve fitting for `A * exp(-gap * t) + C` — provides `curve_fit` with bounds, `confidence_interval`, `standard_error`, and automatic Jacobian computation. The standard Julia package for this task; purpose-built for curve fitting unlike raw Optim.jl.
- **Distributions.jl v0.25+ (transitive)**: t-distribution for confidence intervals — brought in by LsqFit, used internally, not directly imported by QuantumFurnace.
- **Existing trajectory infrastructure (reused)**: `run_trajectories` with observables, `TrajectoryFramework`, `step_along_trajectory!`, `_accumulate_measurements!` — all primitives work as-is; no changes needed to core simulation engine.

**Alternatives considered and rejected:**
- Optim.jl alone: already a dependency, but lacks curve-fitting infrastructure (confidence intervals, covariance estimation). Would require manual implementation of what LsqFit provides.
- Bootstrap.jl: LsqFit's Jacobian-based CIs are sufficient. Bootstrap would require storing per-trajectory data (massive memory) or re-running trajectories (prohibitively expensive).
- Prony's method / matrix pencil: theoretically superior for multi-exponential decomposition but (a) no mature Julia package, (b) notoriously noise-sensitive without SVD regularization, (c) requires equally-spaced samples. Single-exponential fit with window selection is more robust.

### Expected Features

**Must have (table stakes):**
- **Observable-only trajectory runner**: Run trajectories with time-resolved observable measurements (`<O>(t)` at `save_every` intervals) without per-trajectory DM reconstruction. Existing `run_trajectories` with observables already does this; just need to clarify/simplify the API.
- **Total magnetization observable**: `M_z = sum_i Z_i` in eigenbasis/Trotter basis, following the pattern of existing `build_convergence_observables`. Easy addition; one function.
- **Single-exponential fit with bounds**: Model `f(t) = A * exp(-gap * t) + C`, constrain `gap > 0`, auto-initialize from log-linear estimate. Core of the milestone.
- **Fit quality metrics**: R-squared, residual norm, confidence interval on gap. LsqFit provides these directly.
- **Cross-validation against exact Liouvillian**: For n=4,6, compare trajectory-fitted gap vs `run_lindbladian().spectral_gap`. This is the validation that makes the method credible for n≥8.
- **`estimate_spectral_gap` function**: Public API that orchestrates observable-trajectories + multi-observable fitting + best-fit selection + optional cross-validation. Returns `SpectralGapResult` with gap estimate, CI, per-observable results, and metadata.

**Should have (competitive):**
- **Multi-observable consistency check**: Fit gap from energy, M_z, and all ZZ correlations independently; report agreement. Minimal code, high scientific value.
- **Fitting window selection**: Skip early transient (first 10-20%) and late noise floor. Improves fit quality significantly with minimal complexity.
- **Variance-weighted fitting**: Weight by `1/var(O(t_i))` if per-time variance is available. LsqFit supports this via `wt` parameter. Optional refinement.
- **Gap vs beta scaling plot**: For paper figures. Low complexity, deferred to simulation scripts (not library code).

**Defer (v2+):**
- **Damped oscillation model**: `A * exp(-gamma * t) * cos(omega * t + phi) + C` for complex eigenvalues with significant imaginary part. Only needed if pure exponential fits fail (unlikely for KMS-balanced Lindbladians which have real eigenvalues).
- **Multi-exponential fit**: Extracts multiple decay rates simultaneously. Ill-conditioned; only needed if eigenvalue spectrum characterization (not just gap) is the goal.
- **Bootstrap confidence intervals**: Jacobian-based CIs from LsqFit are sufficient for validation. Bootstrap adds complexity without clear benefit for this use case.

### Architecture Approach

The design follows the "compose existing primitives, add minimal new code" principle. All trajectory simulation machinery is reused. The new milestone adds: (1) a new file `src/spectral_gap.jl` containing all gap estimation code (keeps the already-large `trajectories.jl` focused), (2) a variant of the trajectory inner loop `_run_chunk_obs_only!` that measures observables without per-trajectory DM accumulation (optional, for clarity and slight memory efficiency), (3) exponential fitting function `fit_exponential_decay` wrapping LsqFit.jl, (4) result struct `SpectralGapResult` co-located in `spectral_gap.jl`, (5) observable builder extensions in `convergence.jl` for total magnetization, and (6) cross-validation helper that compares fitted gap against `LindbladianResult.spectral_gap`.

**Major components:**
1. **Observable builders (convergence.jl)** — `build_total_magnetization(ham, n)` and `build_gap_estimation_observables(ham, n)` extend existing patterns. Transform observables to eigenbasis/Trotter basis via `V' * O * V` (critical to avoid basis mismatch pitfall).
2. **Trajectory runner variant (spectral_gap.jl)** — `run_observable_trajectories` and `_run_chunk_obs_only!` follow the `run_trajectories` pattern exactly but make DM reconstruction optional. Reuses all existing primitives: `_build_framework_and_seed`, `step_along_trajectory!`, `_accumulate_measurements!`, multi-threading via `_partition_trajectories`.
3. **Exponential fitting (spectral_gap.jl)** — `fit_exponential_decay(times, values)` uses LsqFit.jl `curve_fit` with: log-linear initial guess, parameter bounds `[gap >= 0]`, optional Gibbs value for offset, returns named tuple with gap, CI, SE, residual norm, converged flag.
4. **Top-level API (spectral_gap.jl)** — `estimate_spectral_gap(jumps, config, psi0, ham; observables, ntraj, save_every, exact_result)` runs trajectories, fits all observables, selects best fit (lowest residual, converged, gap > 0), optionally cross-validates, returns `SpectralGapResult`.
5. **Cross-validation helper (spectral_gap.jl)** — `cross_validate_gap(estimated, exact_result::LindbladianResult)` compares `estimated.gap` vs `abs(real(exact_result.spectral_gap))`, warns if `|Im/Re| > 0.1`, returns relative error.

**Integration points:**
- `src/spectral_gap.jl` (NEW): ~250 lines, all gap estimation logic
- `src/convergence.jl` (MODIFY): +2 functions (~40 lines) for observable builders
- `src/QuantumFurnace.jl` (MODIFY): `include("spectral_gap.jl")`, `using LsqFit`, export new API
- `Project.toml` (MODIFY): add LsqFit dependency and compat entry

**No changes needed:**
- `src/trajectories.jl`: all primitives reused as-is
- `src/furnace.jl`: `run_lindbladian` already provides `.spectral_gap` for cross-validation
- `src/structs.jl`: `SpectralGapResult` lives in `spectral_gap.jl` (not a cross-module type)

### Critical Pitfalls

1. **Confusing complex eigenvalue with real decay rate (CRITICAL)** — The Liouvillian `spectral_gap` is `Complex{T}`. Observable decay rate is `-real(spectral_gap)`, NOT `abs(spectral_gap)`. Cross-validation MUST use `abs(real(...))` for comparison. For KMS-balanced Lindbladians (exact BohrDomain), eigenvalues are real; for approximate domains (Energy, Time, Trotter), small imaginary parts can appear. Flag when `|Im/Re| > 0.1` (indicates oscillatory decay, not pure exponential). **Mitigation:** Define comparison quantity explicitly: `exact_gap = abs(real(liouv_result.spectral_gap))`, assert positive, warn on large imaginary part.

2. **Observable basis mismatch (CRITICAL)** — Trajectories evolve in eigenbasis/Trotter basis; observables MUST be transformed to the same basis via `V' * O * V`. This is the exact bug class from v1.2 quick-task-20 (0.83 gap instead of 0.0807 from mixing bases). Total magnetization `M_z = sum_i Z_i` is naturally defined in computational basis; must transform. **Mitigation:** Follow `build_convergence_observables` pattern exactly; add regression test `tr(gibbs * M_z_eigen) = sum_i <Z_i>_gibbs` to catch basis errors.

3. **Multi-exponential contamination (CRITICAL)** — Observable decay is `sum_k c_k * exp(-gamma_k * t)`, not a single exponential. Fitting the full time series gives a biased estimate (weighted average of all rates, pulled high by fast modes). **Mitigation:** (a) Late-time fitting only — skip first 10-20% where fast modes dominate; (b) Fit multiple observables independently, select the smallest fitted rate with good R-squared as the true gap; (c) Use exact steady-state value `<O>_ss = tr(gibbs * O)` as fixed parameter rather than fitting it.

4. **Noise floor at late times (CRITICAL)** — At `t >> 1/gap`, signal `|<O>(t) - <O>_ss| ~ exp(-gap * t)` drops below noise `sigma / sqrt(N_traj)`. Including noisy tail corrupts fit. For gap ~ 0.08, signal reaches noise floor at t ~ 50 with 1000 trajectories. **Mitigation:** (a) Pre-compute signal-to-noise from rough gap estimate, choose `N_traj` such that signal > 3*noise for at least 3 decay times; (b) Use weighted fitting by `1/var(t)` to downweight noisy regions; (c) Empirically determine noise floor and exclude time points where `|signal| < 2*noise`.

5. **Fitting sensitivity to initial guess (MODERATE)** — Levenberg-Marquardt converges to local minima. Poor initial `gap_0` causes convergence to wrong rate or failure. **Mitigation:** (a) Log-linear pre-estimate: fit `log|<O>(t) - C|` vs `t` to get initial gap; (b) Fix offset to known Gibbs value to reduce from 3-parameter to 2-parameter fit; (c) Use LsqFit bounds `lower=[..., 0.0, ...]` to constrain `gap >= 0`.

## Implications for Roadmap

Based on research, suggested phase structure:

### Phase 1: Observable Infrastructure
**Rationale:** Zero-dependency foundation. Observable builders have no external deps and enable all subsequent work. Following the same pattern as existing `build_convergence_observables` minimizes basis-mismatch risk. Can be tested in isolation against Gibbs state traces.

**Delivers:**
- `build_total_magnetization(ham, n)`
- `build_gap_estimation_observables(ham, n)` returning [H, M_z, ZZ_12, ZZ_23, ...]
- Basis transform regression tests

**Addresses:** Total magnetization observable (table stakes), basis mismatch prevention (critical pitfall 2)

**Avoids:** Pitfall 2 (basis mismatch) by following `build_convergence_observables` pattern exactly

**Complexity:** LOW (~40 lines in `convergence.jl`)

---

### Phase 2: Add LsqFit Dependency + Exponential Fitting
**Rationale:** The fitting logic can be developed and tested independently of trajectory simulation using synthetic exponential data. This validates the fitting methodology (initial guess, bounds, convergence) before integrating with noisy trajectory data. Can run in parallel with Phase 1.

**Delivers:**
- LsqFit.jl added to `Project.toml`
- `fit_exponential_decay(times, values; skip_initial, gibbs_value)` in new `spectral_gap.jl`
- Synthetic data tests: fit `y = 2.0 * exp(-0.5 * t) + 1.0 + noise`, verify recovery

**Uses:** LsqFit.jl for `curve_fit`, `confidence_interval`, `standard_error`

**Addresses:** Single-exponential fit with bounds (table stakes), fit quality metrics (table stakes)

**Avoids:** Pitfall 5 (initial guess sensitivity) via log-linear pre-estimate, Pitfall 3 (multi-exponential) via `skip_initial` parameter

**Complexity:** MEDIUM (~70 lines fitting function + tests)

---

### Phase 3: Observable-Only Trajectory Runner
**Rationale:** Depends on Phase 1 for observables but independent of Phase 2 (fitting). Follows `run_trajectories_convergence` template. The trajectory runner produces time-series data that Phase 4 will consume.

**Delivers:**
- `_run_chunk_obs_only!` (observable variant of trajectory loop)
- `run_observable_trajectories(jumps, config, psi0, ham; observables, save_every, ntraj, reconstruct_dm=false)`
- Returns `TrajectoryResult` with `measurements_mean` and `times`

**Uses:** Existing trajectory primitives (`_build_framework_and_seed`, `step_along_trajectory!`, `_accumulate_measurements!`, multi-threading)

**Implements:** Observable-only trajectory runner architecture component

**Addresses:** Observable-only trajectory runner (table stakes), correct time grid selection (pitfall 7 mitigation)

**Avoids:** Pitfall 4 (noise floor) by exposing `ntraj` and `total_time` parameters for signal-to-noise planning

**Complexity:** MEDIUM (~120 lines following existing patterns)

---

### Phase 4: Gap Estimation API + Result Struct
**Rationale:** Depends on Phases 2 and 3 (requires both fitting and trajectory runner). Integrates observable-trajectories with multi-observable fitting, implements best-fit selection logic, packages results.

**Delivers:**
- `SpectralGapResult` struct in `spectral_gap.jl`
- `estimate_spectral_gap(jumps, config, psi0, ham; observables, ntraj, save_every, exact_result)`
- Fits all observables, selects best (lowest residual, converged, gap > 0)
- Returns gap estimate + CI + per-observable results

**Implements:** Top-level API architecture component

**Addresses:** `estimate_spectral_gap` function (table stakes), multi-observable consistency check (differentiator)

**Avoids:** Pitfall 3 (multi-exponential) by comparing gap across observables (smallest rate with good fit is true gap)

**Complexity:** MEDIUM (~80 lines orchestration logic)

---

### Phase 5: Cross-Validation Helpers + n=4,6 Validation
**Rationale:** Depends on Phase 4 (needs `estimate_spectral_gap` API). Cross-validation establishes trust in the method by comparing against exact Liouvillian eigenvalues at n=4,6. This is the scientific validation step.

**Delivers:**
- `cross_validate_gap(estimated, exact_result::LindbladianResult)` helper
- Validation simulation script for n=4 and n=6
- Report: relative error, confidence interval overlap, per-observable gap comparison

**Addresses:** Cross-validation against exact Liouvillian gap (table stakes)

**Avoids:**
- Pitfall 1 (complex vs real) by using `abs(real(spectral_gap))` and warning on large imaginary part
- Pitfall 6 (wrong eigenvalue from Arpack) by requesting more eigenvalues and validating against full spectrum for n=4
- Pitfall 9 (wrong steady-state value) by using `liouv_result.fixed_point` for baseline subtraction in cross-validation tests

**Complexity:** MEDIUM (helper is ~30 lines; validation script is separate, not library code)

---

### Phase 6: Gap Scaling Studies (Optional, Deferred)
**Rationale:** Uses complete implementation from Phases 1-5. This is a simulation/paper phase, not library development. Can be done outside the main milestone delivery.

**Delivers:**
- Gap vs beta scaling plots
- Gap vs n scaling plots (tests system-size independence prediction from arXiv:2510.08533)
- Variance-weighted fitting refinement (if needed for precision)

**Addresses:** Gap scaling plots (differentiator), variance-weighted fitting (differentiator)

**Complexity:** LOW (simulation scripts using stable API)

---

### Phase Ordering Rationale

- **Phases 1 and 2 can be parallel:** Observable builders and exponential fitting are independent. Both are foundational.
- **Phase 3 depends on Phase 1:** Needs observables to measure during trajectories.
- **Phase 4 depends on Phases 2 and 3:** Integrates fitting with trajectory runner.
- **Phase 5 depends on Phase 4:** Cross-validation needs the complete API.
- **Phase 6 is post-delivery:** Paper figures and refinements happen after validation.

**Key dependency chain:**
```
Phase 1 (observables)  ─┐
                        ├─> Phase 3 (trajectories) ─┐
Phase 2 (fitting)     ──┼───────────────────────────┼─> Phase 4 (API) ──> Phase 5 (validation) ──> Phase 6 (scaling)
```

**Risk mitigation via ordering:**
- Address critical basis-mismatch pitfall (2) first via Phase 1 (observable patterns)
- Validate fitting methodology via synthetic data in Phase 2 before applying to noisy trajectory data
- Cross-validation (Phase 5) validates the entire pipeline before extending to n≥8

### Research Flags

**Phases with standard patterns (minimal additional research needed):**
- **Phase 1 (Observables):** Follows `build_convergence_observables` exactly. Basis transform pattern is established. Skip `/gsd:research-phase`.
- **Phase 2 (Fitting):** LsqFit.jl is well-documented, exponential fitting is standard numerics. Skip `/gsd:research-phase`.
- **Phase 3 (Trajectories):** Follows `run_trajectories_convergence` template. All primitives exist. Skip `/gsd:research-phase`.

**Phases needing light validation during planning (but not deep research):**
- **Phase 4 (API):** Multi-observable selection logic (pick best fit by R-squared, converged, gap > 0) is straightforward but may benefit from a quick-task to define the exact selection heuristic. Consider a quick-task for "best-fit selection criteria".
- **Phase 5 (Cross-validation):** The comparison metric (exact vs fitted) needs precision in handling complex eigenvalues and imaginary part warnings. Consider a quick-task to validate the comparison logic for n=4 BohrDomain vs TrotterDomain.

**No phase requires `/gsd:research-phase`:** All technical approaches are well-understood from domain research.

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | LsqFit.jl verified via official docs, GitHub, and tutorial. Dependency analysis confirms minimal footprint. All existing trajectory primitives are production-ready. |
| Features | HIGH | Physics foundation is textbook Lindbladian theory. Feature scope is well-defined by cross-validation at n=4,6. Observable builders follow established patterns. |
| Architecture | HIGH | Direct codebase analysis of all relevant files (trajectories.jl, convergence.jl, furnace.jl, structs.jl). New components follow existing patterns exactly. Integration points are clear. |
| Pitfalls | HIGH | Critical pitfalls identified from (a) codebase analysis (basis mismatch from v1.2 quick-task-20), (b) numerical analysis literature (multi-exponential fitting, noise floor), (c) quantum open systems theory (complex eigenvalues, decay rate vs oscillation frequency). Mitigations are concrete and testable. |

**Overall confidence:** HIGH

The research is grounded in direct codebase analysis (all source files read), verified external dependencies (LsqFit.jl documentation and GitHub), established physics (Lindbladian spectral decomposition, observable decay), and battle-tested numerical methods (Levenberg-Marquardt fitting, log-linear initial guess). The architecture reuses 95% of existing code. The new code (~250 lines in `spectral_gap.jl`, ~40 lines in `convergence.jl`) follows established patterns from `run_trajectories_convergence` and `build_convergence_observables`.

### Gaps to Address

**Multi-observable selection heuristic:** The "best fit" selection logic (from energy, M_z, and ZZ correlations) uses "lowest residual + converged + gap > 0" as the criterion. This is reasonable but may need refinement based on validation results. **Mitigation:** Phase 5 cross-validation will reveal if this heuristic is sufficient. If not, consider weighted average across observables with inverse-variance weighting, or select the observable with best R-squared rather than best residual norm.

**Exact eigenvalue for non-BohrDomain:** The research assumes `run_lindbladian` returns the correct spectral gap for Energy, Time, and TrotterDomain. Pitfall 6 flags that Arpack with shift-invert may misidentify the gap eigenvalue when eigenvalues have similar real parts. **Mitigation:** Phase 5 validation should compute full spectrum (via `eigen()`) for n=4 to verify Arpack's result. If discrepancy found, extend `run_lindbladian` to request more eigenvalues and sort properly.

**Signal-to-noise planning for n=8:** The recommended trajectory counts (5000+ for n=4,6, 10000+ for n=8) are rough estimates. The actual required `N_traj` depends on the specific Hamiltonian's spectral gap and observable variance. **Mitigation:** Phase 5 validation at n=4,6 will calibrate the signal-to-noise relationship. Use this to plan n=8 runs. Consider a pilot run (100 trajectories) before committing to full N_traj.

**Gibbs vs fixed-point for offset:** The exponential fit uses `<O>_ss` as the offset parameter. For exact KMS (BohrDomain), this is the Gibbs value. For approximate domains, it should be the Liouvillian fixed point. The research flags this as Pitfall 9. **Mitigation:** For cross-validation (Phase 5), use `liouv_result.fixed_point` to compute `<O>_ss`. For production use at n≥8 (where Liouvillian is unavailable), use Gibbs value and acknowledge domain approximation error. Document this in the function docstring.

## Sources

### Primary (HIGH confidence)

**Codebase (direct analysis):**
- `src/trajectories.jl` (927 lines) — `run_trajectories`, `_run_chunk_with_obs!`, `step_along_trajectory!`, `_accumulate_measurements!`, trajectory primitives
- `src/convergence.jl` (387 lines) — `build_convergence_observables`, `run_trajectories_convergence`, observable builders and basis transforms
- `src/furnace.jl` (163 lines) — `run_lindbladian`, `LindbladianResult.spectral_gap`, Arpack eigenvalue extraction
- `src/structs.jl` (358 lines) — `LindbladianResult`, `TrajectoryResult`, `ConvergenceData`, struct definitions
- `src/QuantumFurnace.jl` — module structure, exports, dependency imports
- `Project.toml` — current dependencies (Optim.jl, Arpack, etc.)
- `.planning/quick/20-debug-gns-trotterdomain-0-83-gap-suspect/20-SUMMARY.md` — documented basis mismatch bug (0.83 gap instead of 0.0807)

**External dependencies:**
- [LsqFit.jl Documentation](https://julianlsolvers.github.io/LsqFit.jl/latest/) — API reference, tutorials, curve_fit usage
- [LsqFit.jl GitHub v0.15.1](https://github.com/JuliaNLSolvers/LsqFit.jl) — Dependencies, Project.toml, release history, Julia compat
- [LsqFit.jl Tutorial](https://julianlsolvers.github.io/LsqFit.jl/latest/tutorial/) — Exponential model example, parameter bounds, confidence intervals
- [Distributions.jl v0.25 docs](https://juliastats.org/Distributions.jl/v0.25/) — Generated with Julia 1.11.7, confirming compatibility
- [Distributions.jl v0.25.123 on Zenodo](https://zenodo.org/records/18145493) — Latest release Jan 2026

### Secondary (HIGH confidence, domain knowledge)

**Quantum open systems theory:**
- Lindbladian spectral decomposition: `rho(t) = rho_ss + sum c_k R_k exp(lambda_k t)` — textbook result
- Observable decay rate = `-real(lambda)` where `lambda` is Liouvillian eigenvalue — standard open quantum systems
- KMS detailed balance implies real eigenvalues w.r.t. GNS inner product — Chen, Kastoryano, Gilyen (2025) arXiv:2311.09207

**Numerical methods:**
- Exponential fitting of sums of exponentials is ill-conditioned — classic numerical analysis result
- Log-linear initial guess for nonlinear exponential fitting — standard practice in spectroscopy, NMR
- Levenberg-Marquardt for nonlinear least squares — Nocedal & Wright, Numerical Optimization
- Late-time fitting extracts slowest decay mode — standard in quantum Monte Carlo gap estimation

**Physics literature:**
- [Sandvik (2011) "Excitation Gap from Optimized Correlation Functions in QMC"](https://ar5iv.labs.arxiv.org/html/1112.2269) — Signal-to-noise in extracting gaps from Monte Carlo, multi-exponential fitting window selection
- [Nachtergaele, Sims (2006) "Spectral Gap and Exponential Decay of Correlations"](https://link.springer.com/article/10.1007/s00220-006-0030-4) — Mathematical foundation: spectral gap implies exponential correlation decay
- [Fast Mixing of Quantum Spin Chains at All Temperatures](https://arxiv.org/html/2510.08533) — System-size independent gap for 1D chains at finite temperature (testable prediction for gap vs n scaling)

### Tertiary (MEDIUM confidence, methodological references)

- [HypothesisTests.jl parametric tests](https://juliastats.org/HypothesisTests.jl/stable/parametric/) — OneSampleTTest for gap validation (existing test dependency)
- [Mixing Time of Open Quantum Systems via Hypocoercivity](https://arxiv.org/abs/2404.11503) — Relationship between spectral gap, mixing time, and observable autocorrelation
- [Mori (2022) "Liouvillian analysis of relaxation time"](https://www2.yukawa.kyoto-u.ac.jp/~nqs2022/slide/4th/Mori.pdf) — Complex eigenvalue structure, decay rate vs oscillation frequency
- Julia Discourse on weighted LsqFit — Weight parameter usage patterns

---
*Research completed: 2026-02-16*
*Ready for roadmap: YES*
