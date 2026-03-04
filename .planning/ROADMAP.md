# Roadmap: QuantumFurnace.jl

## Milestones

- ✅ **v1.0 Trajectories** -- Phases 1-5 (shipped 2026-02-14)
- ✅ **v1.1 Reduce** -- Phases 6-11 (shipped 2026-02-15)
- ✅ **v1.2 Multi-threading** -- Phases 12-19 (shipped 2026-02-16)
- ✅ **v1.3 Mixing Time Estimation** -- Phases 20-25 (shipped 2026-02-19)
- ✅ **v1.4 Spectral Gap Refinement** -- Phase 26 (partial, shipped 2026-02-20)
- ✅ **v1.5 Krylov Gap Estimation** -- Phases 27-32 (shipped 2026-02-25)
- ✅ **v2.0 Restructure** -- Phases 33-38 (shipped 2026-02-28)
- ✅ **v2.1 Speedup & Mixing Time** -- Phases 39-43 (shipped 2026-03-04)

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

<details>
<summary>✅ v1.4 Spectral Gap Refinement (Phase 26) -- PARTIAL, SHIPPED 2026-02-20</summary>

- [x] Phase 26: Exact Reference and Structural Diagnostics (2/2 plans) -- completed 2026-02-19

Phases 27-30 deferred (two-exponential fitting, effective rates, bootstrap, dashboard).

Full details: [milestones/v1.4-ROADMAP.md](milestones/v1.4-ROADMAP.md)

</details>

<details>
<summary>✅ v1.5 Krylov Gap Estimation (Phases 27-32) -- SHIPPED 2026-02-25</summary>

- [x] Phase 27: Core Matvec Infrastructure (2/2 plans) -- completed 2026-02-20
- [x] Phase 28: Domain Extension (2/2 plans) -- completed 2026-02-20
- [x] Phase 29: Eigensolver Integration (2/2 plans) -- completed 2026-02-24
- [x] Phase 30: Cross-Validation (2/2 plans) -- completed 2026-02-24
- [x] Phase 31: Scaling Benchmarks (2/2 plans) -- completed 2026-02-24
- [x] Phase 32: Krylov Simulator Speedup (2/2 plans) -- completed 2026-02-25

Full details: [milestones/v1.5-ROADMAP.md](milestones/v1.5-ROADMAP.md)

</details>

<details>
<summary>✅ v2.0 Restructure (Phases 33-38) -- SHIPPED 2026-02-28</summary>

- [x] Phase 33: Type Foundation (4/4 plans) -- completed 2026-02-25
- [x] Phase 34: Code Deduplication (2/2 plans) -- completed 2026-02-26
- [x] Phase 35: Workspace and Channel Consolidation (2/2 plans) -- completed 2026-02-27
- [x] Phase 36: API and Results (4/4 plans) -- completed 2026-02-27
- [x] Phase 37: File Organization and Dead Code (2/2 plans) -- completed 2026-02-27
- [x] Phase 38: Test Cleanup (5/5 plans) -- completed 2026-02-28

Full details: [milestones/v2.0-ROADMAP.md](milestones/v2.0-ROADMAP.md)

</details>

### v2.1 Speedup & Mixing Time (Complete, shipped 2026-03-04)

**Milestone Goal:** Optimize run_thermalize performance via per-jump precomputation and multi-threaded BLAS/omega-loops; add mixing time estimation via exponential fit on trace distance convergence curve; bi-exponential fitting for improved extrapolation accuracy.

- [x] **Phase 39: Per-Jump Precomputation** - Eliminate per-step eigendecomposition by precomputing K0, U_residual, and U_coherent per jump at simulation start -- completed 2026-03-01
- [x] **Phase 40: Save Every** - Control trace distance computation frequency to decouple observation cost from step cost -- completed 2026-03-01
- [x] **Phase 41: Threading** - Multi-threaded BLAS and optional omega-loop parallelism for DM thermalization -- completed 2026-03-01
- [x] **Phase 42: Mixing Time Estimation** - Exponential fit on trace distance convergence curve with extrapolation and quality gates -- completed 2026-03-01
- [x] **Phase 43: Bi-Exponential Fitting** - Bi-exponential decay model for improved mixing time extrapolation (<5% error vs 26% single-exp) -- completed 2026-03-04

## Phase Details

### Phase 39: Per-Jump Precomputation
**Goal**: Users run_thermalize without redundant per-step eigendecomposition -- all CPTP channel data is precomputed once per jump at simulation start, producing identical results faster
**Depends on**: Phase 38 (v2.0 complete)
**Requirements**: PRECOMP-01, PRECOMP-02, PRECOMP-03, PRECOMP-04
**Success Criteria** (what must be TRUE):
  1. run_thermalize for Energy, Time, and Trotter domains produces trace distance results matching the pre-precomputation baseline within floating-point tolerance (atol < 1e-12)
  2. No eigen() or _build_cptp_channel call occurs inside the per-step hot loop -- all channel data comes from precomputed vectors of K0s and U_residuals
  3. BohrDomain run_thermalize continues to work correctly with general speedups applied but without per-Bohr-frequency precomputation
  4. Existing test suite passes with zero or near-zero regression test baseline updates (regenerated if O(1e-14) level shifts occur)
**Plans:** 2 plans

Plans:
- [x] 39-01-PLAN.md -- Precomputation infrastructure (BohrDomain _precompute_R, _precompute_per_jump_channels, _accumulate_rho_jump!, _apply_precomputed_channel!)
- [x] 39-02-PLAN.md -- Integration and validation (refactor run_thermalize, CPTP test extension, full test suite)

### Phase 40: Save Every
**Goal**: Users control how often trace distance to the Gibbs state is computed during run_thermalize, reducing observation overhead for long simulations while preserving backward compatibility
**Depends on**: Phase 39
**Requirements**: SAVE-01, SAVE-02
**Success Criteria** (what must be TRUE):
  1. run_thermalize accepts a save_every keyword that controls trace distance computation frequency (e.g., save_every=10 computes trace distance every 10 steps)
  2. Default save_every=1 produces identical results to pre-change behavior (full backward compatibility)
  3. ThermalizeResults contains aligned time_steps and trace_distances arrays whose lengths match ceil(n_steps / save_every) (plus or minus 1 for initial/final step handling)
**Plans:** 1 plan

Plans:
- [x] 40-01-PLAN.md -- Add save_every keyword to run_thermalize with gated trace distance computation and behavioral tests

### Phase 41: Threading
**Goal**: Users get multi-threaded BLAS speedup on DM thermalization automatically, with optional omega-loop parallelism for large frequency grids, without data races or BLAS thread state leaks
**Depends on**: Phase 39 (allocation-free hot path required for safe threading)
**Requirements**: THREAD-01, THREAD-02, THREAD-03, THREAD-04, THREAD-05
**Success Criteria** (what must be TRUE):
  1. run_thermalize with multi-threaded BLAS produces results matching single-threaded execution within expected floating-point tolerance (atol < 1e-10 for accumulated multi-step)
  2. BLAS thread count is always restored after run_thermalize returns, even if an error occurs (try/finally pattern verified)
  3. Threaded omega-loop rho_jump accumulation (when enabled) produces results matching serial accumulation within floating-point tolerance
  4. BohrDomain benefits from threading where applicable (Bohr bucket iteration, sandwich accumulation)
  5. No BLAS thread oversubscription occurs when run_thermalize is called independently or in sequence with trajectory functions
**Plans:** 2 plans

Plans:
- [x] 41-01-PLAN.md -- BLAS thread management for run_thermalize (try/finally save/restore, DM threading tests)
- [x] 41-02-PLAN.md -- Omega-loop parallelism with per-task accumulators (_accumulate_rho_jump!, _precompute_R for all 4 domains)

### Phase 42: Mixing Time Estimation
**Goal**: Users can estimate mixing time from a run_thermalize trace distance curve via exponential fit, with optional early stopping via extrapolation, quality gates on fit reliability, and LsqFit.jl re-integrated as an active dependency
**Depends on**: Phase 40 (save_every produces time-series format needed for fitting)
**Requirements**: MIX-01, MIX-02, MIX-03, MIX-04, MIX-05, MIX-06, MIX-07, MIX-08
**Success Criteria** (what must be TRUE):
  1. estimate_mixing_time(result::ThermalizeResults) returns a MixingTimeEstimate with fitted_gap, mixing_time, R_squared, and confidence intervals
  2. With extrapolate=false (default), mixing_time reflects the actual number of steps to reach target trace distance from the full simulation
  3. With extrapolate=true, simulation stops early when the exponential fit is reliable, and the returned mixing_time is the extrapolated number of steps to target epsilon
  4. Quality gate warnings fire when R^2 < 0.95 or when offset C is large relative to target epsilon
  5. LsqFit.jl is restored as an active dependency in Project.toml and fitting.jl is promoted from src/staging/ to active source
**Plans:** 2 plans

Plans:
- [x] 42-01-PLAN.md -- LsqFit.jl promotion + MixingTimeEstimate implementation (fitting.jl promoted, mixing.jl created)
- [x] 42-02-PLAN.md -- Comprehensive tests and full validation (fitting tests promoted, mixing tests created)

### Phase 43: Bi-Exponential Fitting
**Goal**: Reduce mixing time extrapolation error from ~26% to <5% by adding bi-exponential decay model with explicit :biexp keyword, preserving backward compatibility
**Depends on**: Phase 42 (mixing time estimation)
**Success Criteria** (what must be TRUE):
  1. All existing tests pass unchanged
  2. Bi-exp extrapolation achieves <5% error on synthetic bi-exp data
  3. `model=:single` (default) preserves exact current behavior
**Plans:** 1 plan

Plans:
- [x] PLAN.md -- BiexpFitResult, fit_biexponential_decay, model=:biexp in estimate_mixing_time, Roots.Bisection extrapolation

## Progress

**Execution Order:**
Phases execute in numeric order: 39 -> 40 -> 41 -> 42 -> 43

| Phase | Milestone | Plans Complete | Status | Completed |
|-------|-----------|----------------|--------|-----------|
| 1-5 | v1.0 | 10/10 | Complete | 2026-02-14 |
| 6-11 | v1.1 | 16/16 | Complete | 2026-02-15 |
| 12-19 | v1.2 | 15/15 | Complete | 2026-02-16 |
| 20-25 | v1.3 | 10/10 | Complete | 2026-02-18 |
| 26 | v1.4 | 2/2 | Complete | 2026-02-19 |
| 27-32 | v1.5 | 12/12 | Complete | 2026-02-25 |
| 33-38 | v2.0 | 19/19 | Complete | 2026-02-28 |
| 39. Precomputation | v2.1 | 2/2 | Complete | 2026-03-01 |
| 40. Save Every | v2.1 | 1/1 | Complete | 2026-03-01 |
| 41. Threading | v2.1 | 2/2 | Complete | 2026-03-01 |
| 42. Mixing Time | v2.1 | 2/2 | Complete | 2026-03-01 |
| 43. Bi-Exp Fitting | v2.1 | 1/1 | Complete | 2026-03-04 |
