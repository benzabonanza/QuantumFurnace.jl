# Roadmap: QuantumFurnace.jl

## Milestones

- ✅ **v1.0 Trajectories** -- Phases 1-5 (shipped 2026-02-14)
- ✅ **v1.1 Reduce** -- Phases 6-11 (shipped 2026-02-15)
- ✅ **v1.2 Multi-threading** -- Phases 12-19 (shipped 2026-02-16)
- ✅ **v1.3 Mixing Time Estimation** -- Phases 20-25 (shipped 2026-02-19)
- 🚧 **v1.4 Spectral Gap Refinement** -- Phases 26-30 (in progress)

## Phases

<details>
<summary>✅ v1.0 Trajectories (Phases 1-5) -- SHIPPED 2026-02-14</summary>

- [x] Phase 1: Foundation and Compilation (1/1 plans) -- completed 2026-02-13
- [x] Phase 2: Trajectory Bug Fixes (2/2 plans) -- completed 2026-02-14
- [x] Phase 3: DM Reference Test Suite (3/3 plans) -- completed 2026-02-14
- [x] Phase 4: Trajectory Cross-Validation (2/2 plans) -- completed 2026-02-14
- [x] Phase 5: Statistical Validation and Regression (2/2 plans) -- completed 2026-02-14

Full details: [milestones/v1.0-ROADMAP.md](milestones/v1.0-ROADMAP.md)

</details>

<details>
<summary>✅ v1.1 Reduce (Phases 6-11) -- SHIPPED 2026-02-15</summary>

- [x] Phase 6: Dead Code Pruning (2/2 plans) -- completed 2026-02-15
- [x] Phase 7: DRY Refactoring (2/2 plans) -- completed 2026-02-15
- [x] Phase 8: Struct Simplification (3/3 plans) -- completed 2026-02-15
- [x] Phase 9: Type Parameterization (3/3 plans) -- completed 2026-02-15
- [x] Phase 10: API Surface Cleanup (3/3 plans) -- completed 2026-02-15
- [x] Phase 11: Allocation Optimization (3/3 plans) -- completed 2026-02-15

Full details: [milestones/v1.1-ROADMAP.md](milestones/v1.1-ROADMAP.md)

</details>

<details>
<summary>✅ v1.2 Multi-threading (Phases 12-19) -- SHIPPED 2026-02-16</summary>

- [x] Phase 12: Workspace Refactor (2/2 plans) -- completed 2026-02-15
- [x] Phase 13: Multi-Threaded Trajectory Engine (2/2 plans) -- completed 2026-02-16
- [x] Phase 14: GNS Trajectory Path (1/1 plans) -- completed 2026-02-16
- [x] Phase 15: Data Architecture (2/2 plans) -- completed 2026-02-16
- [x] Phase 16: Convergence Tracking (2/2 plans) -- completed 2026-02-16
- [x] Phase 17: Adaptive Sampling (2/2 plans) -- completed 2026-02-16
- [x] Phase 18: KMS-vs-GNS Experiments (1/1 plans) -- completed 2026-02-16
- [x] Phase 19: Logic Simplification (3/3 plans) -- completed 2026-02-16

Full details: [milestones/v1.2-ROADMAP.md](milestones/v1.2-ROADMAP.md)

</details>

<details>
<summary>✅ v1.3 Mixing Time Estimation (Phases 20-25) -- SHIPPED 2026-02-19</summary>

- [x] Phase 20: Observable Infrastructure (1/1 plans) -- completed 2026-02-17
- [x] Phase 21: Exponential Fitting (1/1 plans) -- completed 2026-02-17
- [x] Phase 22: Observable-Only Trajectory Runner (1/1 plans) -- completed 2026-02-17
- [x] Phase 23: Gap Estimation API (1/1 plans) -- completed 2026-02-17
- [x] Phase 24: Cross-Validation (3/3 plans) -- completed 2026-02-17
- [x] Phase 25: Spectral Gap Validation Overhaul (3/3 plans) -- completed 2026-02-18

Full details: [milestones/v1.3-ROADMAP.md](milestones/v1.3-ROADMAP.md)

</details>

### 🚧 v1.4 Spectral Gap Refinement (In Progress)

**Milestone Goal:** Build comprehensive diagnostics and improved estimation methods for spectral gap estimation, validated against exact results at n=4,6, to produce reliable gap estimates for larger systems.

**Constraints:**
- All trajectory runs use `with_coherent=true`
- 4 threads for trajectory simulations
- Trotter steps = 10 for OFT/B
- Delta is the relevant Trotter parameter
- Storage: ~40 GB free -- warn before large data operations (batch-level bootstrap ~3 MB, fine)

- [ ] **Phase 26: Exact Reference and Structural Diagnostics** -- Exact Lindbladian eigenvalues, fixed point, anti-Hermitian defect, symmetry sector labels, and observable overlap analysis
- [ ] **Phase 27: Two-Exponential Fitting Infrastructure** -- Robust two-exponential decay fitting with Prony initialization and separation quality checks
- [ ] **Phase 28: Effective Rate Plot and Automatic Window Selection** -- Model-free lambda_eff(t) diagnostic with SNR-based t_max and stability-based t_min selection
- [ ] **Phase 29: Batched Bootstrap and Richardson Extrapolation** -- Batch-level bootstrap uncertainty quantification and delta-convergence extrapolation with monotonicity gate
- [ ] **Phase 30: Integration, Dashboard, and Validation** -- 7-panel diagnostic dashboard, external field comparison, and final validated gap table at n=4,6

## Phase Details

### Phase 26: Exact Reference and Structural Diagnostics
**Goal**: Researchers can extract exact Lindbladian spectral data and structural diagnostics (anti-Hermitian defect, symmetry sectors, observable overlaps) as ground truth for validating trajectory-based estimates
**Depends on**: Phase 25 (existing run_lindbladian infrastructure)
**Requirements**: DIAG-01, DIAG-02, DIAG-03, DIAG-04, DIAG-05, DIAG-06
**Success Criteria** (what must be TRUE):
  1. Leading 20-30 Lindbladian eigenvalues (real + imaginary parts) are extracted via Arpack shift-invert at n=4 and n=6, sorted by proximity to zero
  2. Lindbladian fixed point (lambda_1=0 eigenvector) is computed and its trace distance to the Gibbs state is reported (should be ~0 for BohrDomain, ~1e-8 for TrotterDomain)
  3. Anti-Hermitian defect is computed via KMS similarity transform with Gibbs spectrum truncation, and the defect ratio ||A||/lambda_gap(H) determines whether real-exponential fitting is appropriate
  4. Delta-Sz symmetry sector labels are assigned to each eigenvector with near-degeneracy detection, explaining the n=6 zero-overlap mystery
  5. Observable overlap coefficients c_k are computed for leading modes using exact left/right eigenvectors, identifying which observables couple to the gap mode
**Plans**: 2 plans

Plans:
- [ ] 26-01-PLAN.md -- Diagnostics core: result structs, DIAG-01 through DIAG-06 functions, run_exact_diagnostics bundle, tests
- [ ] 26-02-PLAN.md -- Observable set replacement (6 canonical) + eigenbasis_overlap_analysis corrected formula

### Phase 27: Two-Exponential Fitting Infrastructure
**Goal**: Researchers can fit trajectory observable decay to a two-exponential model and trust the separation quality assessment before using the fitted rates
**Depends on**: Phase 25 (v1.3 single-exponential fitting infrastructure), Phase 26 (exact diagnostics inform expected rate separation)
**Requirements**: FIT-01, FIT-02, FIT-03
**Success Criteria** (what must be TRUE):
  1. `fit_two_exponential_decay()` returns a `TwoExpFitResult` with two rates (g1, g2), two amplitudes (c1, c2), confidence intervals, and R-squared
  2. Prony two-point initialization produces starting values that converge reliably (no manual tuning needed for n=4 trajectory data)
  3. When g2/g1 < 1.5, the function automatically falls back to single-exponential tail fit and reports the fallback in the result
  4. Synthetic two-exponential data with known rates (g1=0.17, g2=0.5) is recovered within confidence intervals
**Plans**: TBD

Plans:
- [ ] 27-01: TBD
- [ ] 27-02: TBD

### Phase 28: Effective Rate Plot and Automatic Window Selection
**Goal**: Researchers can compute model-free effective decay rates lambda_eff(t) and automatically select the fitting window (t_min, t_max) without manual parameter tuning
**Depends on**: Phase 27 (two-exponential fit for t_min stability loop), Phase 26 (Lindbladian fixed point for correct steady-state subtraction)
**Requirements**: RATE-01, RATE-02, RATE-03, RATE-04
**Success Criteria** (what must be TRUE):
  1. `compute_effective_rates()` returns lambda_eff(t) with NaN masking at sign changes and noise floor cutoff, using the Lindbladian fixed point as steady-state (not Gibbs state for TrotterDomain)
  2. Bootstrap error bars on lambda_eff(t) are computed via trajectory batch resampling, producing pointwise confidence bands
  3. SNR-based t_max selection (SNR > 3 threshold) automatically truncates the noisy tail of the decay curve
  4. Stability-based t_min selection (plateau detection over window-start sweep using two-exponential g1) removes the early-time fast-mode contamination window
**Plans**: TBD

Plans:
- [ ] 28-01: TBD
- [ ] 28-02: TBD

### Phase 29: Batched Bootstrap and Richardson Extrapolation
**Goal**: Researchers can quantify uncertainty on spectral gap estimates via bootstrap resampling and correct Trotter bias via Richardson extrapolation with mandatory quality gates
**Depends on**: Phase 27 (two-exponential fit per bootstrap resample), Phase 28 (effective rate computation per bootstrap sample)
**Requirements**: BOOT-01, BOOT-02, BOOT-03, RICH-01, RICH-02, RICH-03
**Success Criteria** (what must be TRUE):
  1. `run_observable_trajectories_batched()` stores per-batch mean observable time series in an n_batches x n_obs x n_saves 3D array (~3 MB for 100 batches at n=6)
  2. Batch-level bootstrap with B>=100 resamples produces percentile confidence intervals (not normal approximation) for the spectral gap, with the full distribution (median, mean, SE) reported per observable
  3. Delta-convergence sweep at 3+ delta values with identical parameters (same N_traj, t_final, observable, seed) produces gap estimates that can be compared for monotonicity
  4. Richardson extrapolation is gated by a mandatory monotonicity precondition; un-extrapolated estimates are always reported alongside extrapolated values
**Plans**: TBD

Plans:
- [ ] 29-01: TBD
- [ ] 29-02: TBD
- [ ] 29-03: TBD

### Phase 30: Integration, Dashboard, and Validation
**Goal**: Researchers can run a single diagnostic command and get a validated gap table with a 7-panel summary figure, including external field comparison confirming symmetry sector analysis
**Depends on**: Phase 26, 27, 28, 29 (all prior phases)
**Requirements**: VAL-01, VAL-02, VAL-03, VAL-04
**Success Criteria** (what must be TRUE):
  1. 7-panel diagnostic dashboard (spectrum, defect metrics, overlap coefficients, effective rate, delta-convergence, two-exp fit overlay, t_min stability) is generated for a given system configuration
  2. External field comparison (h=0 vs h=0.1J) at n=4,6 confirms that symmetry breaking improves gap-mode coupling for the n=6 case
  3. Final validated gap table at n=4,6 reports exact vs estimated gap, bootstrap CI, and sigma discrepancy -- with n=4 accuracy maintained at <1%
  4. Multi-observable minimum-gap selector uses two-exponential g1 estimates to pick the best observable, improving on v1.3's single-exponential selection
**Plans**: TBD

Plans:
- [ ] 30-01: TBD
- [ ] 30-02: TBD
- [ ] 30-03: TBD

## Progress

**Execution Order:**
Phases execute in numeric order: 26 -> 27 -> 28 -> 29 -> 30

| Phase | Milestone | Plans Complete | Status | Completed |
|-------|-----------|----------------|--------|-----------|
| 1. Foundation and Compilation | v1.0 | 1/1 | Complete | 2026-02-13 |
| 2. Trajectory Bug Fixes | v1.0 | 2/2 | Complete | 2026-02-14 |
| 3. DM Reference Test Suite | v1.0 | 3/3 | Complete | 2026-02-14 |
| 4. Trajectory Cross-Validation | v1.0 | 2/2 | Complete | 2026-02-14 |
| 5. Statistical Validation and Regression | v1.0 | 2/2 | Complete | 2026-02-14 |
| 6. Dead Code Pruning | v1.1 | 2/2 | Complete | 2026-02-15 |
| 7. DRY Refactoring | v1.1 | 2/2 | Complete | 2026-02-15 |
| 8. Struct Simplification | v1.1 | 3/3 | Complete | 2026-02-15 |
| 9. Type Parameterization | v1.1 | 3/3 | Complete | 2026-02-15 |
| 10. API Surface Cleanup | v1.1 | 3/3 | Complete | 2026-02-15 |
| 11. Allocation Optimization | v1.1 | 3/3 | Complete | 2026-02-15 |
| 12. Workspace Refactor | v1.2 | 2/2 | Complete | 2026-02-15 |
| 13. Multi-Threaded Trajectory Engine | v1.2 | 2/2 | Complete | 2026-02-16 |
| 14. GNS Trajectory Path | v1.2 | 1/1 | Complete | 2026-02-16 |
| 15. Data Architecture | v1.2 | 2/2 | Complete | 2026-02-16 |
| 16. Convergence Tracking | v1.2 | 2/2 | Complete | 2026-02-16 |
| 17. Adaptive Sampling | v1.2 | 2/2 | Complete | 2026-02-16 |
| 18. KMS-vs-GNS Experiments | v1.2 | 1/1 | Complete | 2026-02-16 |
| 19. Logic Simplification | v1.2 | 3/3 | Complete | 2026-02-16 |
| 20. Observable Infrastructure | v1.3 | 1/1 | Complete | 2026-02-17 |
| 21. Exponential Fitting | v1.3 | 1/1 | Complete | 2026-02-17 |
| 22. Observable-Only Trajectory Runner | v1.3 | 1/1 | Complete | 2026-02-17 |
| 23. Gap Estimation API | v1.3 | 1/1 | Complete | 2026-02-17 |
| 24. Cross-Validation | v1.3 | 3/3 | Complete | 2026-02-17 |
| 25. Spectral Gap Validation Overhaul | v1.3 | 3/3 | Complete | 2026-02-18 |
| 26. Exact Reference and Structural Diagnostics | v1.4 | 0/TBD | Not started | - |
| 27. Two-Exponential Fitting Infrastructure | v1.4 | 0/TBD | Not started | - |
| 28. Effective Rate Plot and Window Selection | v1.4 | 0/TBD | Not started | - |
| 29. Batched Bootstrap and Richardson Extrapolation | v1.4 | 0/TBD | Not started | - |
| 30. Integration, Dashboard, and Validation | v1.4 | 0/TBD | Not started | - |
