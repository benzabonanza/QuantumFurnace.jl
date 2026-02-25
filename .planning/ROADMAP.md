# Roadmap: QuantumFurnace.jl

## Milestones

- ✅ **v1.0 Trajectories** -- Phases 1-5 (shipped 2026-02-14)
- ✅ **v1.1 Reduce** -- Phases 6-11 (shipped 2026-02-15)
- ✅ **v1.2 Multi-threading** -- Phases 12-19 (shipped 2026-02-16)
- ✅ **v1.3 Mixing Time Estimation** -- Phases 20-25 (shipped 2026-02-19)
- ✅ **v1.4 Spectral Gap Refinement** -- Phase 26 (partial, shipped 2026-02-20)
- ✅ **v1.5 Krylov Gap Estimation** -- Phases 27-32 (shipped 2026-02-25)
- 🚧 **v2.0 Restructure** -- Phases 33-38 (in progress)

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

### 🚧 v2.0 Restructure (In Progress)

**Milestone Goal:** Major codebase restructure -- redesign Config type hierarchy for extensibility (KMS/GNS/DLL), eliminate code duplication across simulation paths, consolidate workspaces, reorganize files, and slim down tests. Prepare architecture for DLL construction, error estimation, and gate complexity features.

- [x] **Phase 33: Type Foundation** - Define Config{S,D,C,T} hierarchy with simulation/construction singleton types, hard swap all call sites -- completed 2026-02-25
- [ ] **Phase 34: Code Deduplication** - Extract domain_prefactor(), foreach_frequency(), and unified oft!() replacing 16+ copy-pasted patterns
- [ ] **Phase 35: Workspace and Channel Consolidation** - Merge KrylovWorkspace + KrausScratch + LindbladianWorkspace; unify R/K0/U_residual computation paths
- [ ] **Phase 36: API and Results** - Define 4 clean run_* entry points with matching Result structs and save capability
- [ ] **Phase 37: File Organization and Dead Code** - Rename src/ files to PRE/MID/POST grouping, move staging code, remove @distributed, update exports
- [ ] **Phase 38: Test Cleanup** - Consolidate test helpers, add @info printouts, review dubious thresholds

## Phase Details

### Phase 33: Type Foundation
**Goal**: All simulation code dispatches on a single `Config{S,D,C,T}` struct, eliminating 4 duplicate config types and enabling future DLL construction via a new singleton type
**Depends on**: Nothing (first phase of v2.0)
**Requirements**: TYPE-01, TYPE-02, TYPE-03, TYPE-04, TYPE-05, TYPE-06
**Success Criteria** (what must be TRUE):
  1. `Config{Lindbladian,EnergyDomain,KMS,Float64}` can be constructed and used in place of `LiouvConfig{EnergyDomain,Float64}` for all 4 domain types
  2. `Config{Thermalize,TimeDomain,GNS,Float64}` can be constructed and used in place of `ThermalizeConfigGNS{TimeDomain,Float64}` for all 4 domain types
  3. `with_coherent` is derived from construction type (KMS -> true, GNS -> false) at compile time, not stored as a field
  4. All existing tests pass using either the new `Config{S,D,C,T}` type or backward-compatible aliases (`LiouvConfig`, `ThermalizeConfig`, `LiouvConfigGNS`, `ThermalizeConfigGNS`)
  5. Adding a `DLL` construction type requires only `struct DLL <: AbstractConstruction end` plus dispatch methods, with zero changes to existing code
**Plans**: 4 plans

Plans:
- [x] 33-01-PLAN.md -- Define type hierarchies, Config struct, with_coherent trait, update exports
- [x] 33-02-PLAN.md -- Migrate core dispatch (energy_domain, bohr_domain, furnace_utensils, coherent, misc_tools, results)
- [x] 33-03-PLAN.md -- Migrate simulation pipeline (furnace, jump_workers, krylov, trajectories, convergence, gap_estimation)
- [x] 33-04-PLAN.md -- Migrate tests, simulations, experiments, playground; run full test suite

### Phase 34: Code Deduplication
**Goal**: The 16+ copy-pasted prefactor formulas, hermitian half-grid branching patterns, and OFT variants are each single-source functions dispatched on domain type
**Depends on**: Phase 33 (new Config type needed for dispatch signatures)
**Requirements**: DEDUP-01, DEDUP-02, DEDUP-03
**Success Criteria** (what must be TRUE):
  1. `domain_prefactor(config, gamma_norm_factor)` is called from all 5 files (jump_workers, krylov_workspace, krylov_matvec, krylov_eigsolve, trajectories) instead of inline formula computation
  2. `foreach_frequency()` iterator replaces the 16 hermitian half-grid branching patterns with zero allocation overhead verified by existing allocation tests
  3. A single `oft!()` function (with domain dispatch) replaces both `oft!` and `_krylov_oft!`; `time_oft!`/`trotter_oft!` remain as clearly-marked test/debug utilities
  4. All regression, allocation, and Krylov cross-validation tests pass with identical numerical results
**Plans**: TBD

Plans:
- [ ] 34-01: TBD
- [ ] 34-02: TBD

### Phase 35: Workspace and Channel Consolidation
**Goal**: Workspace types are consolidated with unified naming, and R/K0/U_residual CPTP channel computation uses shared helper functions with correct per-jump vs summed semantics
**Depends on**: Phase 34 (deduplication helpers used by workspace constructors)
**Requirements**: WORK-01, WORK-02
**Success Criteria** (what must be TRUE):
  1. `KrylovWorkspace`, `KrausScratch`, and `LindbladianWorkspace` are consolidated into a unified workspace struct (TrajectoryWorkspace stays separate)
  2. Shared CPTP channel helper functions compute R, K0, U_residual correctly: per-jump (R^a, K0^a, U_residual^a) for DM/Trajectory paths, summed (R_total, K0_total, U_residual_total) for Krylov path
  3. Zero-allocation hot-path invariants are preserved: existing allocation tests pass with 0 bytes in all domains
  4. Workspace independence test passes (no shared mutable state between workspaces)
**Plans**: TBD

Plans:
- [ ] 35-01: TBD
- [ ] 35-02: TBD

### Phase 36: API and Results
**Goal**: Four clean public entry points (`run_lindblad`, `run_thermalize`, `run_krylov_spectrum`, `run_trajectory`) each return a typed Result struct with optional BSON save capability
**Depends on**: Phase 35 (consolidated workspaces used by run_* functions)
**Requirements**: WORK-03, WORK-04, WORK-05
**Success Criteria** (what must be TRUE):
  1. `run_lindblad()`, `run_thermalize()`, `run_krylov_spectrum()`, `run_trajectory()` are the 4 public entry points, each dispatching on `Config{S,D,C,T}`
  2. Each entry point returns a typed Result struct (`LindbladResults`, `ThermalizeResults`, `KrylovSpectrumResults`, `TrajectoryResults`) containing the config and metadata (git hash, timestamp)
  3. Passing `save_path="path.bson"` to any `run_*` function serializes the result to BSON with companion .txt, and the result round-trips correctly via `load_result`
  4. Simulation scripts in `simulations/` demonstrate all 4 entry points with working examples
**Plans**: TBD

Plans:
- [ ] 36-01: TBD
- [ ] 36-02: TBD

### Phase 37: File Organization and Dead Code
**Goal**: Source files are renamed for clarity with PRE/MID/POST logical grouping, dead code is removed, staging code is separated, and the module export list matches the new structure
**Depends on**: Phase 36 (API must be stable before renaming files)
**Requirements**: ORG-01, ORG-02, ORG-03, ORG-07, ORG-08, ORG-09
**Success Criteria** (what must be TRUE):
  1. All src/ files are renamed with clear names reflecting their role in the PRE (types, physics)/MID (simulation runners)/POST (results, analysis) architecture -- flat directory, no subdirectories
  2. Gap estimation/fitting/convergence code is moved to a staging area separated from active source
  3. `@distributed` dead code and `using Distributed` import are removed from furnace.jl; SharedArrays import stays
  4. Diagnostics remains as a separate analysis module, not folded into the Lindblad simulation path
  5. Module export list is reorganized by simulation type (Lindblad, Thermalize, Krylov, Trajectory) with clean groupings
**Plans**: TBD

Plans:
- [ ] 37-01: TBD
- [ ] 37-02: TBD

### Phase 38: Test Cleanup
**Goal**: Test infrastructure is consolidated with parametrized helpers, informative output, and validated thresholds so that tests are maintainable and trustworthy
**Depends on**: Phase 37 (tests cleaned after code is finalized)
**Requirements**: ORG-04, ORG-05, ORG-06
**Success Criteria** (what must be TRUE):
  1. Config factory functions in test_helpers.jl are parametrized by system size and construction type, eliminating duplicate setup patterns across test files
  2. Every `@testset` block prints `@info` showing what is being tested and key numerical results (trace distances, gaps, allocation counts)
  3. Previously dubious test thresholds are reviewed and either tightened with documented rationale or relaxed with explanation of why the original threshold was wrong
  4. All tests pass and the full test suite runs without regressions from the restructure
**Plans**: TBD

Plans:
- [ ] 38-01: TBD
- [ ] 38-02: TBD

## Progress

**Execution Order:**
Phases execute in numeric order: 33 -> 34 -> 35 -> 36 -> 37 -> 38

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
| 26. Exact Reference and Structural Diagnostics | v1.4 | 2/2 | Complete | 2026-02-19 |
| 27. Core Matvec Infrastructure | v1.5 | 2/2 | Complete | 2026-02-20 |
| 28. Domain Extension | v1.5 | 2/2 | Complete | 2026-02-20 |
| 29. Eigensolver Integration | v1.5 | 2/2 | Complete | 2026-02-24 |
| 30. Cross-Validation | v1.5 | 2/2 | Complete | 2026-02-24 |
| 31. Scaling Benchmarks | v1.5 | 2/2 | Complete | 2026-02-24 |
| 32. Krylov Simulator Speedup | v1.5 | 2/2 | Complete | 2026-02-25 |
| 33. Type Foundation | v2.0 | 4/4 | Complete | 2026-02-25 |
| 34. Code Deduplication | v2.0 | 0/TBD | Not started | - |
| 35. Workspace and Channel Consolidation | v2.0 | 0/TBD | Not started | - |
| 36. API and Results | v2.0 | 0/TBD | Not started | - |
| 37. File Organization and Dead Code | v2.0 | 0/TBD | Not started | - |
| 38. Test Cleanup | v2.0 | 0/TBD | Not started | - |
