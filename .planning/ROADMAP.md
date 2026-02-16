# Roadmap: QuantumFurnace.jl

## Milestones

- ✅ **v1.0 Trajectories** -- Phases 1-5 (shipped 2026-02-14)
- ✅ **v1.1 Reduce** -- Phases 6-11 (shipped 2026-02-15)
- 🚧 **v1.2 Multi-threading** -- Phases 12-18 (in progress)

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

### 🚧 v1.2 Multi-threading (In Progress)

**Milestone Goal:** Multi-threaded trajectory sampling engine with GNS comparison path, enabling paper-ready KMS-vs-GNS convergence experiments across temperature regimes and system sizes.

- [x] **Phase 12: Workspace Refactor** - Separate mutable workspace from shared framework to enable thread-safe trajectory stepping -- completed 2026-02-15
- [x] **Phase 13: Multi-Threaded Trajectory Engine** - Parallel trajectory sampling with per-thread workspaces, BLAS control, and deterministic seeding -- completed 2026-02-16
- [x] **Phase 14: GNS Trajectory Path** - Verify and test GNS (approximate detailed balance) trajectory simulation end-to-end -- completed 2026-02-16
- [ ] **Phase 15: Data Architecture** - Experiment result serialization and round-trip via BSON
- [ ] **Phase 16: Convergence Tracking** - Batch-level trace distance and per-observable convergence metrics during trajectory sampling
- [ ] **Phase 17: Adaptive Sampling** - Convergence-driven trajectory batching with automatic stopping and hard cap
- [ ] **Phase 18: KMS-vs-GNS Experiments** - Parameter sweep experiments comparing KMS and GNS across system sizes and temperatures

## Phase Details

### Phase 12: Workspace Refactor
**Goal**: Trajectory stepping accepts explicit workspace and RNG arguments, enabling safe concurrent execution from a single shared framework
**Depends on**: Phase 11 (v1.1 codebase)
**Requirements**: THRD-01
**Success Criteria** (what must be TRUE):
  1. `step_along_trajectory!` accepts explicit `TrajectoryWorkspace` and `AbstractRNG` arguments instead of pulling them from the framework
  2. Two independent workspaces can step trajectories from the same `TrajectoryFramework` without interfering with each other
  3. All existing trajectory tests pass unchanged (backward compatibility preserved)
**Plans**: 2 plans
  - [x] 12-01-PLAN.md -- Core struct/function refactor (remove ws from framework, add ws+rng params, TrajectoryResult)
  - [x] 12-02-PLAN.md -- Test migration and workspace independence verification

### Phase 13: Multi-Threaded Trajectory Engine
**Goal**: Users can run thousands of trajectories in parallel across CPU threads with reproducible results
**Depends on**: Phase 12
**Requirements**: THRD-02, THRD-03, THRD-04
**Success Criteria** (what must be TRUE):
  1. Multi-threaded trajectory sampling distributes N trajectories across available Julia threads, each with its own workspace and RNG
  2. BLAS thread count is set to 1 before parallel execution and restored afterward, preventing thread oversubscription
  3. Given the same master seed and thread count, multi-threaded trajectory results are bitwise identical across runs
  4. Multi-threaded execution at n=4 with 4+ threads is faster than serial execution (no performance regression from threading overhead)
**Plans**: 2 plans
  - [x] 13-01-PLAN.md -- Core threaded engine (allocation guard, partition helper, @sync/@spawn dispatch, BLAS control)
  - [x] 13-02-PLAN.md -- Threading correctness and performance tests (determinism, BLAS restore, serial-threaded agreement, speedup)

### Phase 14: GNS Trajectory Path
**Goal**: GNS (approximate detailed balance, no coherent B term) trajectory simulation works end-to-end and produces physically correct results
**Depends on**: Phase 12
**Requirements**: GNS-01, GNS-02
**Success Criteria** (what must be TRUE):
  1. `ThermalizeConfigGNS` dispatches through the full trajectory pipeline (workspace allocation, stepping, density matrix accumulation) without error
  2. GNS trajectories produce valid density matrices (Hermitian, unit trace, positive semidefinite) at every checkpoint
  3. Averaged GNS trajectory density matrix converges toward the GNS approximate fixed point (not exact Gibbs, but within documented approximation bound for the given sigma)
**Plans**: 1 plan
  - [x] 14-01-PLAN.md -- GNS test infrastructure + full trajectory validation suite (fixed point, CPTP, convergence, DM validity)

### Phase 15: Data Architecture
**Goal**: Experiment results (configs, convergence curves, observables, density matrices) are persistable and reproducible from saved files
**Depends on**: Phase 13
**Requirements**: DATA-01, DATA-02
**Success Criteria** (what must be TRUE):
  1. An experiment result containing config, convergence history, observable time series, and final density matrix can be saved to a BSON file
  2. A saved experiment result can be loaded back and all fields match the original (configs, matrices, scalar values) to machine precision
  3. Saved files include sufficient metadata (seed, thread count, Julia version, timestamp) to reproduce or contextualize the result
**Plans**: TBD

### Phase 16: Convergence Tracking
**Goal**: Trajectory sampling reports trace distance to Gibbs and per-observable values at batch checkpoints, giving users visibility into convergence progress
**Depends on**: Phase 13
**Requirements**: CONV-01, CONV-02, CONV-03
**Success Criteria** (what must be TRUE):
  1. Trace distance between the running average density matrix and the Gibbs state is computed and recorded at each batch checkpoint
  2. Nearest-neighbor correlation `<Z_iZ_{i+1}>` is tracked per batch and its value converges as trajectory count increases
  3. Energy expectation `<H>` is tracked per batch and converges toward the thermal equilibrium value
  4. Convergence data (trace distance curve, observable curves) is accessible programmatically after a trajectory run completes
**Plans**: TBD

### Phase 17: Adaptive Sampling
**Goal**: Trajectory sampling automatically runs enough batches to reach convergence without the user specifying a fixed trajectory count
**Depends on**: Phase 15, Phase 16
**Requirements**: CONV-04, CONV-05
**Success Criteria** (what must be TRUE):
  1. Adaptive mode runs trajectory batches and stops when relative change in tracked metrics is below 1% for 3 consecutive batches
  2. A hard maximum trajectory cap prevents infinite loops even when convergence is slow
  3. Adaptive mode returns the same result structure as fixed-count mode (convergence history, final density matrix, observables)
**Plans**: TBD

### Phase 18: KMS-vs-GNS Experiments
**Goal**: Paper-ready comparison data showing KMS and GNS convergence behavior across system sizes and inverse temperatures
**Depends on**: Phase 14, Phase 17
**Requirements**: EXPT-01, EXPT-02, EXPT-03, EXPT-04
**Success Criteria** (what must be TRUE):
  1. Matched KMS and GNS(sigma=1/beta) experiments run on the same Hamiltonian, beta, and delta, producing comparable convergence data
  2. GNS experiments also run at sigma=0.5/beta, giving a second comparison point for cost-accuracy tradeoff
  3. The full parameter grid (n=4,6,8 x beta=5,10,20 x {KMS, GNS@1/beta, GNS@0.5/beta}) executes using TrotterDomain with pre-built Hamiltonians
  4. Results demonstrate that KMS achieves lower final trace distance to Gibbs than GNS at sigma=1/beta (validating exact vs approximate detailed balance)
  5. All 18 experiment results are saved to BSON and loadable for analysis
**Plans**: TBD

## Progress

**Execution Order:**
Phases execute in numeric order: 12 -> 13 -> 14 -> 15 -> 16 -> 17 -> 18
(Phase 14 depends only on 12, so may overlap with 13 in practice)

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
| 15. Data Architecture | v1.2 | 0/TBD | Not started | - |
| 16. Convergence Tracking | v1.2 | 0/TBD | Not started | - |
| 17. Adaptive Sampling | v1.2 | 0/TBD | Not started | - |
| 18. KMS-vs-GNS Experiments | v1.2 | 0/TBD | Not started | - |
