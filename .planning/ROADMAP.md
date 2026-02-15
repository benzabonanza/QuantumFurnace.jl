# Roadmap: QuantumFurnace.jl

## Milestones

- ✅ **v1.0 Trajectories** -- Phases 1-5 (shipped 2026-02-14)
- ✅ **v1.1 Reduce** -- Phases 6-11 (shipped 2026-02-15)

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

### ✅ v1.1 Reduce (Complete)

**Milestone Goal:** Refactor and simplify the codebase -- prune dead code, simplify structs, introduce type parameterization, clean up public API, remove duplication, and minimize allocations in core simulation paths. All 224 tests must pass after every phase.

- [x] **Phase 6: Dead Code Pruning** - Remove ~930 lines of commented code, ~35 unused functions, and dead structs -- completed 2026-02-15
- [x] **Phase 7: DRY Refactoring** - Extract shared patterns into reusable helpers -- completed 2026-02-15
- [x] **Phase 8: Struct Simplification** - Reduce config field duplication, fix HamHam initialization, make TrottTrott immutable -- completed 2026-02-15
- [x] **Phase 9: Type Parameterization** - Parameterize core structs on element type for F64/F32 flexibility -- completed 2026-02-15
- [x] **Phase 10: API Surface Cleanup** - Remove dead exports, internalize implementation details, export missing public functions -- completed 2026-02-15
- [x] **Phase 11: Allocation Optimization** - Eliminate unnecessary allocations in hot paths -- completed 2026-02-15

## Phase Details

### Phase 6: Dead Code Pruning
**Goal**: Codebase contains only live, reachable code -- no commented-out blocks, no unused functions, no dead structs
**Depends on**: Phase 5 (v1.0 complete)
**Requirements**: PRUNE-01, PRUNE-02, PRUNE-03
**Success Criteria** (what must be TRUE):
  1. No commented-out code blocks remain in source files (the ~930 lines across 9 files are removed)
  2. All exported and internal functions are reachable from either the public API or the test suite
  3. LindbladianJumpCaches, LiouvLiouv structs and non-mutating `oft` wrapper no longer exist in the codebase
  4. All 224 existing tests pass with no regressions
**Plans**: 2 plans
- [x] 06-01-PLAN.md -- Remove ~930 lines of commented-out code blocks across 10 source files
- [x] 06-02-PLAN.md -- Remove ~35 unused functions, 3 dead structs, and update exports

### Phase 7: DRY Refactoring
**Goal**: Repeated code patterns are extracted into single-source helpers -- Hermitianization, CPTP channel application, coherent unitary application, and Trotter basis transforms each have one canonical implementation
**Depends on**: Phase 6
**Requirements**: DRY-01, DRY-02, DRY-03, DRY-04
**Success Criteria** (what must be TRUE):
  1. A single `hermitianize!` helper replaces all 8+ inline Hermitianization patterns
  2. CPTP channel application (K0/residual/Cholesky sequence) exists as one shared function, called from 3 sites in jump_workers.jl
  3. Coherent unitary application exists as one shared function, called from 3 sites in jump_workers.jl
  4. Trotter basis transform of jumps has one canonical implementation shared between furnace.jl and trajectories.jl
  5. All 224 existing tests pass with no regressions
**Plans**: 2 plans
- [x] 07-01-PLAN.md -- Extract hermitianize! helper (DRY-01) and transform_jumps_to_basis helper (DRY-04)
- [x] 07-02-PLAN.md -- Extract apply_coherent_unitary! (DRY-03) and apply_cptp_channel! (DRY-02) helpers in jump_workers.jl

### Phase 8: Struct Simplification
**Goal**: Core data structures are minimal and correct -- config struct field duplication reduced, HamHam fully initialized in constructor, TrottTrott immutable with correct field types
**Depends on**: Phase 7
**Requirements**: STRUCT-01, STRUCT-02, STRUCT-03
**Success Criteria** (what must be TRUE):
  1. Config structs share common fields through composition or shared constructors -- no 12-field duplication across 4 types (separate types preserved for KMS/GNS/future Ding extensibility)
  2. HamHam constructor computes bohr_freqs and gibbs state directly -- no Nothing-typed fields, no two-step initialization pattern
  3. TrottTrott is an immutable struct and num_trotter_steps_per_t0 is typed as Int (not Float64)
  4. All 224 existing tests pass with no regressions
**Plans**: 3 plans
Plans:
- [x] 08-01-PLAN.md -- Config struct deduplication + sentinel cleanup + TrottTrott immutability (STRUCT-01, STRUCT-03)
- [x] 08-02-PLAN.md -- HamHam initialization redesign (STRUCT-02)
- [x] 08-03-PLAN.md -- TrajectoryFramework type param cleanup + domain dispatch refactoring

### Phase 9: Type Parameterization
**Goal**: Core structs are parameterized on element type, enabling future Float32 paths without changing calling code
**Depends on**: Phase 8
**Requirements**: TYPE-01, TYPE-02, TYPE-03
**Success Criteria** (what must be TRUE):
  1. HamHam is parameterized as `HamHam{T<:AbstractFloat}` and all its numeric fields use type T
  2. LindbladianWorkspace is parameterized on element type, consistent with HamHam's type parameter
  3. Config structs are parameterized on float type for tolerance and step-size fields
  4. Existing Float64 call sites work without modification (T defaults to or infers Float64)
  5. All 224 existing tests pass with no regressions
**Plans**: 3 plans
Plans:
- [x] 09-01-PLAN.md -- Parameterize HamHam{T} and TrottTrott{T} on element type
- [x] 09-02-PLAN.md -- Parameterize Config structs, LindbladianWorkspace, JumpOp, KrausScratch, NUFFTPrefactors
- [x] 09-03-PLAN.md -- Propagate T through simulation function signatures and verify full pipeline

### Phase 10: API Surface Cleanup
**Goal**: Public API exposes exactly what users and researchers need -- building blocks for pedagogy are exported, implementation details are internal
**Depends on**: Phase 9
**Requirements**: API-01, API-02, API-03
**Success Criteria** (what must be TRUE):
  1. Non-mutating `oft` wrapper and dead struct names are no longer exported
  2. ~18 implementation-detail exports (workspaces, precompute helpers, internal dispatch functions) are removed from the public API
  3. `trace_distance_h` is exported and accessible for convergence analysis workflows
  4. All 224 existing tests pass with no regressions
**Plans**: 3 plans
Plans:
- [x] 10-01-PLAN.md -- Reorganize export block: add physics exports, remove internal exports, confirm dead exports absent
- [x] 10-02-PLAN.md -- Apply _ prefix to internal function definitions and intra-file call sites
- [x] 10-03-PLAN.md -- Update cross-file call sites and test qualified access for _-prefixed names

### Phase 11: Allocation Optimization
**Goal**: Core simulation hot paths avoid unnecessary heap allocations -- sparse matrices, Diagonal wrappers, filter intermediates, and redundant basis transforms are eliminated or precomputed
**Depends on**: Phase 10
**Requirements**: ALLOC-01, ALLOC-02, ALLOC-03, ALLOC-04
**Success Criteria** (what must be TRUE):
  1. coherent_bohr inner loop does not allocate sparse matrices per iteration (A_nu matrices pre-allocated or precomputed)
  2. B_time and B_trotter closures compute phase rotations in-place without allocating Diagonal wrappers
  3. Time/Trotter thermalize hot path in jump_workers.jl avoids the abs.(filter(...)) intermediate allocation
  4. B_trotter multi-jump variant precomputes Trotter basis transforms instead of recomputing in the inner loop
  5. All 224 existing tests pass with no regressions
**Plans**: 3 plans
Plans:
- [x] 11-01-PLAN.md -- Eliminate sparse allocations in B_bohr and filter+abs in jump_workers
- [x] 11-02-PLAN.md -- Eliminate Diagonal wrappers and redundant basis transforms in B_time/B_trotter
- [x] 11-03-PLAN.md -- Add @allocated regression tests for all optimized hot paths

## Progress

**Execution Order:**
Phases execute in numeric order: 6 -> 7 -> 8 -> 9 -> 10 -> 11

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
