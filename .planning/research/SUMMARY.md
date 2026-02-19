# Project Research Summary

**Project:** QuantumFurnace.jl v1.4 Spectral Gap Refinement Diagnostics
**Domain:** Quantum Lindbladian simulation — spectral gap estimation diagnostics and improved estimators
**Researched:** 2026-02-19
**Confidence:** HIGH

## Executive Summary

QuantumFurnace.jl v1.4 is a diagnostic and improved-estimation milestone built on top of a working v1.3 single-exponential fitting pipeline. The core problem, confirmed by Quick-30 through Quick-32, is that trajectory-based spectral gap estimation achieves ~0.72% accuracy at n=4 but suffers 37-49% error at n=6. Root cause analysis (Quick-32) conclusively separated two error sources: trajectory simulation errors are O(1e-3) and monotonically decrease with Trotter step size delta, while fitting model errors are O(10-50%) and non-monotonic in delta. The v1.4 milestone attacks the fitting bottleneck directly through seven interlocking diagnostic capabilities: effective rate plots, two-exponential fitting with Prony initialization, bootstrap error bars, automatic fitting window selection, anti-Hermitian defect computation, Delta-Sz symmetry sector labeling, and Richardson extrapolation.

The recommended approach is a pure additive layer on top of the existing simulation infrastructure. All new capabilities are post-hoc analyses of trajectory data and exact Lindbladian results — none require modifying the performance-critical trajectory step loop. The sole architectural exception is bootstrap error bars, which require a new batched trajectory runner variant that stores per-batch mean observable time series. The complete existing Julia dependency set (LsqFit.jl, LinearAlgebra stdlib, Arpack.jl, Plots.jl in extras, StatsBase in extras) is sufficient for all v1.4 features. Zero new production dependencies are needed.

The three critical risks are: (1) two-exponential fitting parameter non-identifiability when the two decay rates are within ~2x of each other, which requires two-stage initialization and a separation quality check before any downstream use of the fitted rates; (2) numerical blow-up of the rho^{-1/4} similarity transform at low temperature (beta > 5), which requires Gibbs eigenvalue spectrum truncation and must be validated against the BohrDomain ground truth; and (3) Richardson extrapolation amplifying fitting artifacts rather than correcting Trotter bias, which requires a monotonicity precondition test across 3+ delta values before any extrapolation is applied. These three pitfalls produce silently wrong results that look plausible, making proactive prevention essential.

## Key Findings

### Recommended Stack

The existing dependency set covers all v1.4 needs without exception. Two-exponential fitting uses LsqFit.jl's `curve_fit` with a 5-parameter model — the same API already used for single-exponential fitting in v1.3. Matrix fourth roots for the KMS similarity transform are computed via `LinearAlgebra.eigen(Hermitian(...))` exploiting the diagonal structure of the Gibbs state in the energy eigenbasis, with a numerically stable path for TrotterDomain via general eigendecomposition. Expanding the eigenvalue extraction from nev=2 to nev=30 requires only a parameter change to the existing Arpack.jl shift-invert call. Bootstrap resampling uses `StatsBase.sample(replace=true)` already in the test extras. Diagnostic figures use Plots.jl already in extras.

**Core technologies:**
- **LsqFit.jl (0.15):** Two-exponential fitting — same `curve_fit` API with 5-parameter model `c1*exp(-g1*t) + c2*exp(-g2*t)`; Prony two-point method provides robust starting values for the ill-conditioned optimization
- **LinearAlgebra (stdlib):** Matrix fourth root via `eigen(Hermitian(...))` exploiting diagonal Gibbs state structure; `opnorm` for anti-Hermitian defect norm; `Diagonal` for efficient similarity transforms in eigenbasis
- **Arpack.jl (0.5.4):** Leading 20-30 eigenvalue extraction via shift-invert; only change is increasing `nev` from 2 to 30 in the existing `eigs` call; KrylovKit.jl explicitly rejected (no shift-invert mode for interior eigenvalues)
- **StatsBase.jl (0.34, extras):** `sample(1:K, K; replace=true)` for batch-level bootstrap resampling; 10-line custom loop preferred over Bootstrap.jl for non-standard batch-level use case
- **Plots.jl (1, extras):** 7-panel diagnostic dashboard via `plot(layout=...)` — sufficient for thesis-quality diagnostic figures; CairoMakie.jl deferred (adds ~30 transitive deps for marginal quality gain)
- **Statistics (stdlib):** `mean`, `std` for bootstrap statistics and SNR computation — always available, no new dep entry needed

### Expected Features

The feature landscape is precisely specified by the reference documents. All seven table-stakes features are required for the milestone to deliver value over v1.3. The effective rate plot is the keystone feature: it feeds initialization data into the two-exponential fit, defines the fitting window, and provides model-free validation of all fitted results. The dependency chain flows from the effective rate plot outward to bootstrap, then to Richardson extrapolation.

**Must have (table stakes):**
- **Effective rate plot lambda_eff(t)** — model-free diagnostic `lambda_eff(t) = -(1/tau)*ln|Delta(t+tau)/Delta(t)|` that reveals multi-exponential time regimes and identifies the golden fitting window without fitting assumptions; keystone for all other features
- **Two-exponential fit with Prony initialization** — root cause fix from Quick-32; model `c1*exp(-g1*t) + c2*exp(-g2*t)` absorbs fast-mode contamination into the second term, analogous to lattice QCD multi-state analysis; Prony two-point method gives reliable initial rates without iteration
- **Bootstrap error bars (batch-level)** — nonparametric uncertainty quantification using K~100 batches of trajectories, resampling at batch level to avoid 640 MB+ per-trajectory storage; provides correct confidence intervals for the ill-conditioned two-exponential fit
- **Automatic fitting window selection (SNR + stability)** — removes fragile manual `skip_initial=0.1`; t_max via SNR > 3 threshold, t_min via gamma_1 stability test sweeping window starts
- **Anti-Hermitian defect computation** — KMS similarity transform `D = rho^{-1/4} L[rho^{1/4}(.)rho^{1/4}] rho^{-1/4}`, decomposed as H+A; ratio ||A||/lambda_gap(H) determines whether real-exponential fitting is appropriate
- **Delta-Sz symmetry sector labeling** — explains the n=6 zero-overlap mystery from Quick-25 through Quick-27; labels each Lindbladian eigenvector by the dominant Delta_Sz quantum number
- **Richardson extrapolation for Trotter bias** — formula `gap_rich = 2*gap(delta/2) - gap(delta)` eliminates O(delta) Trotter error; only viable after two-exponential fitting validates monotonic error structure

**Should have (add after core validation):**
- **Summary dashboard (7-panel figure)** — composes all diagnostics into one thesis-quality communicable result
- **External field comparison (h=0.1J)** — confirms symmetry sector restriction as dominant n=6 error source
- **Multi-observable minimum-gap selector** — extend `_select_best_observable` to use two-exponential g1 estimates

**Defer (v1.5+):**
- **n=8 sparse Lindbladian** — explicitly deferred in reference documents; requires KrylovKit.jl and new architecture
- **Damped-oscillation fit model** — only if anti-Hermitian defect proves significant in practice
- **GEVP / matrix pencil methods** — only if two-exponential fit proves insufficient

### Architecture Approach

The diagnostic layer is purely additive — a post-hoc analysis layer consuming existing result structs and trajectory data. Two new source files are added (`src/diagnostics.jl` for all structural diagnostics, `src/bootstrap.jl` for the batched trajectory runner and bootstrap analysis). Three existing files receive additions with minimal modification risk: `src/fitting.jl` gets `fit_two_exponential_decay()`, `src/trajectories.jl` gets `run_observable_trajectories_batched()`, and `src/QuantumFurnace.jl` gets include and export additions. All new result structs are flat immutable Julia structs following the existing convention (no methods on structs, plain field types).

**Major components:**
1. **`src/fitting.jl` (modified — new function)** — adds `fit_two_exponential_decay()` returning `TwoExpFitResult`; 5-parameter model with Prony + effective-rate fallback initialization; separation test (g2/g1 > 1.5); kept separate from existing `fit_exponential_decay` to preserve unchanged API
2. **`src/diagnostics.jl` (new, 500-700 LOC)** — all structural diagnostic functions: `compare_gap_to_exact`, `compute_effective_rates`, `anti_hermitian_defect` (with Gibbs spectrum truncation), `sector_gap_analysis` (with degeneracy detection), `trotter_convergence_sweep`, `run_gap_diagnostics` wrapper; all corresponding result structs
3. **`src/bootstrap.jl` (new, 300-400 LOC)** — `run_observable_trajectories_batched` (new trajectory runner storing per-batch means in `n_batches x n_obs x n_saves` 3D array), `bootstrap_spectral_gap`, `BootstrapGapResult` with percentile CI; `BatchedTrajectoryResult`
4. **Include order in QuantumFurnace.jl** — `fitting.jl` then `gap_estimation.jl` then `bootstrap.jl` then `diagnostics.jl`, following the dependency DAG

The critical architectural decision is batch-level bootstrap over per-trajectory storage: `n_batches x n_obs x n_saves` (~3 MB for 100 batches) versus `n_traj x n_obs x n_saves` (640 MB+ for 20k trajectories at n=6). This enables genuine bootstrap distributions with percentile confidence intervals at ~1% of the naive memory cost, using the block bootstrap methodology validated for independent batches.

### Critical Pitfalls

1. **Two-exponential fit parameter non-identifiability** — when g2/g1 < 1.5, the Levenberg-Marquardt Jacobian is ill-conditioned and `stderror(fit)` from LsqFit dramatically underestimates uncertainty; running from 5 different initial guesses gives 5 different (g1, g2) pairs with similar residuals. Prevention: two-stage initialization (single-exp tail fit for g1, residual fit for g2), explicit `g2 > g1` bounds in `curve_fit`, mandatory separation test with fallback to single-exponential tail fit when g2/g1 < 1.5. Never use LsqFit `stderror` for two-exponential uncertainty — use bootstrap exclusively.

2. **rho^{-1/4} blow-up at low temperature** — the KMS similarity transform amplifies Gibbs eigenvalues by the -1/4 power; at beta=10 for the n=6 Heisenberg chain, rho^{-1/4} diagonal entries reach ~10^5, drowning signal in Float64 numerical noise. Prevention: floor on Gibbs eigenvalues at ~1e-12 before taking the -1/4 power; report condition number max(d_k)/min(d_k); validate that BohrDomain defect is < 1e-10 after truncation (the physical ground truth). This truncation must be designed before any anti-Hermitian defect computations.

3. **Richardson extrapolation amplifying fitting artifacts** — Quick-30 proved the gap estimation error does NOT follow O(delta^p) structure due to fitting model bias; Richardson applied to poorly-initialized two-exponential fits can produce results worse than the finest-delta un-extrapolated estimate. Prevention: mandatory delta-convergence diagnostic (plot fitted gap vs 3+ delta values) before any extrapolation; gate extrapolation on monotonicity check; always report un-extrapolated estimates alongside extrapolated ones; flag when extrapolated gap falls outside [min, max] of un-extrapolated estimates.

4. **Bootstrap memory blow-up from per-trajectory storage** — naive implementation storing N_traj individual time series at n=6 uses 640 MB+ and makes bootstrap infeasible for N_traj > 10k. Prevention: batch-level bootstrap architecture from the start; `run_observable_trajectories_batched` must be the only supported bootstrap input path; the existing per-trajectory grand-mean runner is not extended.

5. **Effective rate lambda_eff(t) divergence at sign changes** — when trajectory noise causes Delta(t) to cross zero, the log-ratio diverges to NaN/Inf. Prevention: noise floor cutoff (exclude t where |Delta(t)| < 3*sigma_traj/sqrt(N)), sign-change guard returning NaN with companion boolean mask, and the steady-state value must be the Lindbladian fixed point (not Gibbs state) for TrotterDomain to prevent a systematic offset that corrupts all lambda_eff values.

## Implications for Roadmap

Based on the feature dependency graph in FEATURES.md and the build order in ARCHITECTURE.md, five phases are natural. The effective rate plot is the keystone and must exist before bootstrap (which needs to compute lambda_eff per bootstrap sample). The two-exponential fitter must exist before both the effective rate plot (feeds its initialization) and bootstrap (fits each bootstrap resample). Richardson extrapolation requires reliable per-delta gap estimates and is placed last among functional components.

### Phase 1: Two-Exponential Fitting Infrastructure
**Rationale:** Zero dependencies on new code; pure numerical function testable immediately with synthetic data. The `TwoExpFitResult` struct established here is consumed by every downstream phase. Validating convergence and separation quality checks here prevents silently wrong results from propagating through the entire pipeline. This is the lowest-risk high-value starting point.
**Delivers:** `fit_two_exponential_decay()` in `fitting.jl` with Prony two-point initialization, two-stage effective-rate fallback, separation test (g2/g1 > 1.5 with fallback to single-exponential), LsqFit `confint` for parameter uncertainty, `TwoExpFitResult` struct.
**Addresses:** Two-exponential fit with Prony init (P0 feature)
**Avoids:** Pitfall 1 (non-identifiability — separation test and two-stage init), Pitfall 8 (negative amplitude cancellation — cancellation ratio check)
**Research flag:** Skip — well-documented numerical method; Prony algorithm is textbook; LsqFit API verified from official docs.

### Phase 2: Exact Reference and Structural Diagnostics
**Rationale:** These diagnostics (exact gap comparison, anti-Hermitian defect, symmetry sector labels, expanded eigenvalue extraction to nev=30) are all independent of each other and of the fitting pipeline. They consume existing `construct_lindbladian` and `run_lindbladian` infrastructure without modification. Building them before the fitting pipeline provides the ground truth references needed to validate Phase 3 and 4 results — in particular, the Lindbladian fixed point computed here is the correct steady-state for lambda_eff in Phase 3.
**Delivers:** `compare_gap_to_exact()` with `GapComparisonResult` including Lindbladian fixed point; `anti_hermitian_defect()` with Gibbs spectrum truncation and condition number reporting; `sector_gap_analysis()` with `SectorAnalysisResult` and near-degeneracy detection; `run_lindbladian` nev increased to 30; `_thermalize_to_liouv_config()` helper. New file `src/diagnostics.jl` starts here.
**Addresses:** Anti-Hermitian defect (P1), Delta-Sz symmetry labels (P1), store 20-30 eigenvalues (P2)
**Avoids:** Pitfall 2 (rho^{-1/4} blow-up — spectrum truncation built in from the start), Pitfall 7 (symmetry label ambiguity — degeneracy detection built in), Pitfall 12 (steady-state mismatch — fixed point available for Phase 3)
**Research flag:** Targeted audit needed for `_thermalize_to_liouv_config` field mapping — a field-by-field verification of ThermalizeConfig vs LiouvConfig struct definitions is needed before implementation. Skip research-phase for the rest (standard linear algebra).

### Phase 3: Effective Rate Plot and Automatic Window Selection
**Rationale:** The effective rate plot is the keystone diagnostic that provides model-free ground truth for all fitting validation. It must exist before the bootstrap pipeline because bootstrap computes lambda_eff curves per bootstrap sample. Automatic window selection is bundled here because it depends directly on the lambda_eff infrastructure and is used by Phase 1's two-exponential fit for t_min determination.
**Delivers:** `compute_effective_rates()` in `diagnostics.jl` with NaN masking, noise floor cutoff at 3*sigma, midpoint time axis convention, sign-change guard; `EffectiveRateResult` struct; SNR-based t_max selector (SNR > 3); t_min stability test (g1 plateau detection over window-start sweep); steady-state value defaulting to Phase 2's Lindbladian fixed point.
**Uses:** Phase 1's `fit_two_exponential_decay` for t_min stability loop; Phase 2's Lindbladian fixed point for correct steady-state subtraction
**Avoids:** Pitfall 3 (lambda_eff divergence — noise floor cutoff and sign-change guard), Pitfall 6 (silent window selection failure — multiple-window comparison), Pitfall 12 (steady-state mismatch — fixed point from Phase 2), Pitfall 14 (wrong time axis — midpoint convention)
**Research flag:** Skip — the lambda_eff computation is 15 lines of arithmetic; window selection logic is fully specified in reference documents.

### Phase 4: Batched Bootstrap and Richardson Extrapolation
**Rationale:** Bootstrap requires the new batched trajectory runner, which is the highest-risk new code (touches trajectories.jl) and has the most significant memory implications. Placing it after Phases 1-3 means all analytical diagnostics are available as validation tools when debugging bootstrap outputs. Richardson extrapolation is bundled here because it is a one-line formula once reliable per-delta gap estimates exist from the two-exponential fitter.
**Delivers:** `run_observable_trajectories_batched()` with `n_batches x n_obs x n_saves` 3D batch storage in new `src/bootstrap.jl`; `bootstrap_spectral_gap()` with percentile CI (not normal approximation); `BootstrapGapResult` including per-bootstrap-sample lambda_eff confidence bands; Richardson extrapolation with mandatory monotonicity precondition gate; `BatchedTrajectoryResult` struct.
**Uses:** Phase 1 `fit_two_exponential_decay` per bootstrap resample; Phase 3 `compute_effective_rates` for per-sample lambda_eff curves
**Implements:** Bootstrap architecture component from ARCHITECTURE.md
**Avoids:** Pitfall 4 (Richardson failure — monotonicity gate required before extrapolation), Pitfall 5 (memory blow-up — batch-level storage from the start), Pitfall 9 (bootstrap lambda_eff bands — per-batch time series enables pointwise CI), Pitfall 13 (degenerate bootstrap samples — minimum B=100 batches)
**Research flag:** Targeted audit needed for per-batch seeding arithmetic — the `master_seed + batch_idx * ntraj_per_batch` scheme must be verified against the existing `Xoshiro(seed + traj_id)` seeding in `trajectories.jl` before implementation to ensure batch independence. Skip research-phase for bootstrap resampling logic itself (textbook block bootstrap).

### Phase 5: Integration, Dashboard, and Validation
**Rationale:** After all diagnostic components are independently validated, the final phase composes them into the `run_gap_diagnostics()` convenience wrapper and the 7-panel summary dashboard. This is the last phase because the summary figure must display results from all prior phases simultaneously and the external field comparison requires all diagnostics to be operational.
**Delivers:** `run_gap_diagnostics()` wrapper returning `DiagnosticReport` bundling all prior result structs; 7-panel Plots.jl dashboard (spectrum, defect metrics, overlap coefficients, effective rate, delta-convergence, two-exp fit overlay, t_min stability); external field comparison (h=0.1J symmetry-breaking validation); multi-observable minimum-gap selector using two-exponential g1 estimates; final validated gap table for n=4,6 comparing exact vs estimated with sigma discrepancy.
**Addresses:** Summary dashboard (P2), external field comparison (P2), multi-observable min-gap selector (P2)
**Research flag:** Skip — Plots.jl multi-panel layout uses documented API already validated in existing docs/ usage; composition of existing structs is mechanical.

### Phase Ordering Rationale

- Phase 1 first because `TwoExpFitResult` is a dependency of Phases 3, 4, and 5, and it is the easiest to test in isolation with synthetic data without any trajectory infrastructure.
- Phase 2 before Phase 3 because the Lindbladian fixed point from Phase 2's `compare_gap_to_exact` is needed as the correct steady-state value for lambda_eff computation in Phase 3. Using the wrong steady-state (Gibbs instead of fixed point for TrotterDomain) corrupts all effective rate plots.
- Phase 3 before Phase 4 because the batched trajectory runner (Phase 4) needs to compute per-bootstrap-sample lambda_eff curves using Phase 3's `compute_effective_rates`. Having the lambda_eff infrastructure ready first avoids implementing it twice.
- Phase 4 before Phase 5 because the dashboard plots bootstrap confidence bands and Richardson-extrapolated estimates that are produced in Phase 4.
- This ordering also matches the ARCHITECTURE.md suggested build order: arch phases 1 (two-exp fitting) maps to roadmap Phase 1; arch phases 2-5 (effective rates, exact reference, defect, symmetry) map to roadmap Phases 2-3; arch phases 6-7 (batched runner, bootstrap) map to roadmap Phase 4; arch phase 9 (integration) maps to roadmap Phase 5.

### Research Flags

Phases needing targeted verification during planning:
- **Phase 2 (config conversion):** The `_thermalize_to_liouv_config` helper must map all fields between `ThermalizeConfig` and `LiouvConfig` correctly. A field-by-field audit of both struct definitions in `src/` is needed before implementation. Silent field mismatch means running diagnostics against a wrong Lindbladian. This is a 30-minute verification task, not a full research-phase.
- **Phase 4 (batched trajectory seeding):** The per-batch seed arithmetic must produce statistically independent batches, compatible with the existing `Xoshiro(seed + traj_id)` scheme. Incorrect seeding makes bootstrap confidence intervals invalid. This is also a targeted audit of the seeding code in `trajectories.jl`, not a full research-phase.

Phases with standard patterns (skip research-phase):
- **Phase 1:** Prony algorithm is textbook; LsqFit bounded fitting is documented from official sources with verified API.
- **Phase 3:** lambda_eff computation is ~15 lines of arithmetic; window selection logic is fully specified in reference documents.
- **Phase 5:** Plots.jl multi-panel dashboard uses documented API already validated in existing `docs/` usage; composition is mechanical.

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | All packages verified against official Julia docs; zero new dependencies confirmed by cross-checking each v1.4 feature against existing stdlib and extras; version compatibility verified for all API calls |
| Features | HIGH | Feature set precisely specified in `spectral-gap-refinements-instructions.md` and `error_catalogue_spectral_gap_estimation.md`; priorities validated by Quick-30/31/32 empirical data with exact gap comparison |
| Architecture | HIGH | Based on direct analysis of all 26 source files and 19 test files; complete data flow traced through existing codebase; all integration points identified; build order validated against dependency DAG |
| Pitfalls | HIGH | Critical pitfalls grounded in Quick-30/31/32 empirical results plus established numerical analysis literature; mitigations are concrete, testable, and cross-referenced to specific code locations |

**Overall confidence:** HIGH

### Gaps to Address

- **Config field mapping (ThermalizeConfig to LiouvConfig):** The `_thermalize_to_liouv_config` helper needs a field-by-field audit against the actual struct definitions before Phase 2 implementation. This is mechanical verification but the consequence of a silent mismatch is running all structural diagnostics against a wrong Lindbladian with no error signal.
- **Per-batch seeding arithmetic:** The exact seeding scheme for `run_observable_trajectories_batched` needs explicit design verification that batches are statistically independent. Correlated batches silently invalidate bootstrap confidence intervals. Audit the `Xoshiro(seed + traj_id)` logic in `trajectories.jl` before writing the batched runner.
- **Two-exponential identifiability at n=4:** The n=4 Heisenberg chain has g1=0.173 and g2~0.35 (ratio ~2x), right at the boundary of identifiability per the numerical literature. Whether two-stage Prony initialization reliably recovers g1 at this ratio should be validated early in Phase 1 with actual n=4 trajectory data (not just synthetic). If not, profile likelihood over a g1 grid may be needed.
- **TrotterDomain fixed point for n=6:** The effective rate computation requires the Lindbladian fixed point for correct steady-state subtraction. For n=6, the Liouvillian is 4096x4096 (dense eigen takes ~30s but is feasible). Confirming `liouv_result.fixed_point` is available from `run_lindbladian` for TrotterDomain with nev=30 should be verified in Phase 2 before Phase 3 depends on it.

## Sources

### Primary (HIGH confidence)
- QuantumFurnace.jl codebase — direct analysis of all 26 source files and 19 test files; all integration points and data flows traced
- `supplementary-informations/spectral-gap-refinements-instructions.md` — 5-part diagnostic pipeline specification; feature and architecture ground truth
- `supplementary-informations/error_catalogue_spectral_gap_estimation.md` — 7 catalogued error sources; pitfall grounding
- Quick-30/31/32 summaries — empirical validation of root cause (fitting model error >> simulation error); Richardson extrapolation failure mechanism confirmed
- Julia LinearAlgebra documentation — eigen, Hermitian, opnorm, Diagonal; matrix fourth root via eigenbasis approach verified
- LsqFit.jl official documentation — curve_fit API, bounds, confidence intervals; two-exponential model approach verified
- Arpack.jl documentation — eigs with sigma (shift-invert) and nev parameters; KrylovKit.jl rejection verified by explicit "no shift-invert" documentation quote
- KrylovKit.jl eigenvalue documentation — confirms no shift-invert mode; interior eigenvalue finding not supported
- StatsBase.jl documentation — `sample` function with `replace` keyword; batch bootstrap pattern verified

### Secondary (MEDIUM confidence)
- Prony's method (Wikipedia + SIAM J. Sci. Comput. paper) — two-point exponential initialization; mathematical foundation sound and standard
- Richardson extrapolation for Lindbladian Trotter error (arXiv:2507.22341, DOI:10.1103/kw39-yxq5) — theoretical support for extrapolation when error structure is correct
- Symmetry classification of many-body Lindbladians (PhysRevX.13.031019) — symmetry sector labeling framework; Delta-Sz approach validated
- KMS detailed balance and similarity transform (Chen et al. arXiv:2303.18224) — rho^{-1/4} transform mathematical foundation
- Bootstrap resampling theory (Efron 1979; Efron and Tibshirani 1993) — block bootstrap validity for independent blocks; minimum ~50 blocks for reliable percentile CIs
- Parameter identifiability in two-exponential models (BMC Systems Biology 2015) — condition number growth as function of rate separation ratio; g2/g1 > 1.5 threshold grounded here

### Tertiary (MEDIUM-LOW confidence)
- Julia Discourse sparse eigenvalue comparison (Arpack vs KrylovKit vs ArnoldiMethod) — community discussion supporting Arpack choice; not official source
- Plots.jl vs Makie.jl ecosystem comparison (Julia Discourse) — community discussion supporting Plots.jl for diagnostic purposes; not official benchmark

---
*Research completed: 2026-02-19*
*Ready for roadmap: yes*
