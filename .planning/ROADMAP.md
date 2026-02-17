# Roadmap: QuantumFurnace.jl

## Milestones

- ✅ **v1.0 Trajectories** -- Phases 1-5 (shipped 2026-02-14)
- ✅ **v1.1 Reduce** -- Phases 6-11 (shipped 2026-02-15)
- ✅ **v1.2 Multi-threading** -- Phases 12-19 (shipped 2026-02-16)
- 🚧 **v1.3 Mixing Time Estimation** -- Phases 20-24 (in progress)

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

### 🚧 v1.3 Mixing Time Estimation (In Progress)

**Milestone Goal:** Estimate the Lindbladian spectral gap from trajectory-based observable decay, cross-validated against exact Liouvillian eigenvalues for small systems.

- [x] **Phase 20: Observable Infrastructure** - Total magnetization and combined gap estimation observables in both eigenbasis and Trotter basis -- completed 2026-02-17
- [x] **Phase 21: Exponential Fitting** - LsqFit.jl-based single-exponential decay fitting with auto-initialization, window selection, and quality metrics -- completed 2026-02-17
- [x] **Phase 22: Observable-Only Trajectory Runner** - Time-resolved observable measurement from trajectories without mid-simulation DM reconstruction -- completed 2026-02-17
- [x] **Phase 23: Gap Estimation API** - Unified `estimate_spectral_gap` function orchestrating trajectories, multi-observable fitting, and best-estimate selection -- completed 2026-02-17
- [ ] **Phase 24: Cross-Validation** - Verify trajectory-fitted gap against exact Liouvillian spectral gap for n=4 and n=6

## Phase Details

### Phase 20: Observable Infrastructure
**Goal**: Users can construct all observables needed for spectral gap estimation, correctly transformed to the simulation basis
**Depends on**: Nothing (first phase of v1.3; builds on existing `build_convergence_observables` pattern from v1.2)
**Requirements**: OBS-01, OBS-02, OBS-03
**Success Criteria** (what must be TRUE):
  1. User can call `build_total_magnetization(ham, n)` and get M_z in Hamiltonian eigenbasis, verified by `tr(gibbs * M_z) = sum_i <Z_i>_gibbs`
  2. User can call `build_total_magnetization(ham, n; trotter)` and get M_z in Trotter eigenbasis for TrotterDomain simulations
  3. User can call `build_gap_estimation_observables(ham, n)` and receive all gap-relevant observables (H, M_z, ZZ correlations) in one call, ready to pass to trajectory runner
  4. All observables pass basis-transform regression tests (trace against Gibbs state matches expected physical values)
**Plans**: 1 plan -- completed 2026-02-17

Plans:
- [x] 20-01-PLAN.md -- Implement build_total_magnetization and build_gap_estimation_observables with exports and regression tests

### Phase 21: Exponential Fitting
**Goal**: Users can fit exponential decay curves to time-series data and extract decay rates with confidence intervals, validated on synthetic data before touching real trajectories
**Depends on**: Nothing (parallel with Phase 20; uses LsqFit.jl on synthetic data, no trajectory dependency)
**Requirements**: FIT-01, FIT-02, FIT-03, FIT-04, FIT-05
**Success Criteria** (what must be TRUE):
  1. User can call `fit_exponential_decay(times, values)` and recover a known decay rate from synthetic data (`A * exp(-gap * t) + C + noise`) to within the returned confidence interval
  2. Fitting auto-generates initial guess via log-linear estimate -- user does not need to provide initial parameters for typical exponential decay data
  3. User can specify `skip_initial` fraction to exclude early transient, and the fitted gap changes measurably (demonstrating window selection works)
  4. Fit result includes R-squared, confidence interval on gap, and standard error; gap is always constrained positive (gap > 0 enforced via parameter bounds)
  5. LsqFit.jl is added to Project.toml with proper compat bounds
**Plans**: 1 plan -- completed 2026-02-17

Plans:
- [x] 21-01-PLAN.md -- Implement FitResult struct, fit_exponential_decay with LsqFit.jl, and synthetic data test suite

### Phase 22: Observable-Only Trajectory Runner
**Goal**: Users can run trajectory simulations that measure time-resolved observables efficiently, without the overhead of per-trajectory density matrix reconstruction
**Depends on**: Phase 20 (needs observable builders to provide observables for measurement)
**Requirements**: TRAJ-01, TRAJ-02, TRAJ-03
**Success Criteria** (what must be TRUE):
  1. User can call `run_observable_trajectories(jumps, config, psi0, ham; observables, save_every, ntraj)` and receive time-resolved observable means `<O>(t)` at regular intervals
  2. Observable-only runner produces identical observable time series as the existing `run_trajectories` with observables (cross-check for correctness)
  3. User can optionally reconstruct the averaged density matrix at the end of the run (but not during), via `reconstruct_dm=true`
  4. Multi-threaded execution works with the same per-thread workspace/RNG pattern as existing `run_trajectories`
**Plans**: 1 plan -- completed 2026-02-17

Plans:
- [x] 22-01-PLAN.md -- Implement run_observable_trajectories with ObservableTrajectoryResult struct and cross-validation tests

### Phase 23: Gap Estimation API
**Goal**: Users can estimate the Lindbladian spectral gap from a single function call that orchestrates trajectory simulation, multi-observable fitting, and best-estimate selection
**Depends on**: Phase 21 (fitting), Phase 22 (trajectory runner)
**Requirements**: GAP-01, GAP-02, GAP-03
**Success Criteria** (what must be TRUE):
  1. User can call `estimate_spectral_gap(jumps, config, psi0, ham; ...)` and receive a `SpectralGapResult` containing the gap estimate, confidence interval, and per-observable fit details
  2. `SpectralGapResult` struct contains: gap estimate, CI, per-observable gaps, R-squared values, best observable name, and fit metadata (ntraj, total_time, save_every)
  3. Best observable is automatically selected by fit quality (R-squared, convergence status, gap > 0), and the selection is explainable from the per-observable results
**Plans**: 1 plan -- completed 2026-02-17

Plans:
- [x] 23-01-PLAN.md -- Implement SpectralGapResult, estimate_spectral_gap orchestration, and integration tests with SMALL system

### Phase 24: Cross-Validation
**Goal**: Users can verify that trajectory-fitted spectral gap agrees with exact Liouvillian eigenvalues, establishing trust in the method for systems where exact diagonalization is infeasible
**Depends on**: Phase 23 (needs complete `estimate_spectral_gap` API)
**Requirements**: VAL-01, VAL-02, VAL-03
**Success Criteria** (what must be TRUE):
  1. User can call `cross_validate_gap(estimated, exact_result)` comparing fitted gap against `abs(real(exact_result.spectral_gap))`, with relative error reported
  2. Cross-validation warns when the imaginary part of the exact eigenvalue is significant (`|Im/Re| > 0.1`), indicating oscillatory decay that pure exponential fitting cannot capture
  3. Validation script demonstrates gap estimation agreement for n=4 and n=6 Heisenberg chains, with fitted gap within confidence interval of exact gap (or within documented tolerance)
**Plans**: TBD

Plans:
- [ ] 24-01: TBD

## Progress

**Execution Order:**
Phases 20 and 21 can execute in parallel (independent). Then 22 (depends on 20), then 23 (depends on 21+22), then 24 (depends on 23).

```
Phase 20 (observables)  ---\
                            +--> Phase 22 (trajectories) --\
Phase 21 (fitting)      ---+-----------------------------+--> Phase 23 (API) --> Phase 24 (validation)
```

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
| 24. Cross-Validation | v1.3 | 0/? | Not started | - |
