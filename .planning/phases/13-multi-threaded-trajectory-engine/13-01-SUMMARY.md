---
phase: 13-multi-threaded-trajectory-engine
plan: 01
subsystem: trajectories
tags: [threading, @spawn, BLAS, Xoshiro, zero-allocation, function-barrier]

# Dependency graph
requires:
  - phase: 12-workspace-refactor
    provides: "TrajectoryFramework read-only, explicit ws/rng workspace independence"
provides:
  - "_partition_trajectories helper for dividing work across threads"
  - "_run_chunk_no_obs! and _run_chunk_with_obs! chunk workers"
  - "Threaded dispatch in run_trajectories with per-task workspace and Xoshiro RNG"
  - "BLAS.set_num_threads(1) with try/finally restore pattern"
  - "Zero-allocation step_along_trajectory! via concrete-typed TrajectoryFramework fields"
affects: [14-convergence-checking, 15-adaptive-batching]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Function barrier pattern for type-stable access to abstract-typed struct fields"
    - "Concrete-typed hot-path fields in TrajectoryFramework (scaled_prefactor, transition, energy_labels, oft_nufft_prefactors)"
    - "Per-trajectory Xoshiro(seed + traj_id) seeding for serial/threaded reproducibility"
    - "@sync/@spawn with per-task workspace for thread-safe trajectory execution"

key-files:
  created: []
  modified:
    - "src/trajectories.jl"
    - "test/test_allocation.jl"

key-decisions:
  - "Added concrete-typed fields (F,P type params) to TrajectoryFramework to eliminate allocation in step_along_trajectory! hot path"
  - "Changed jumps field from Vector{JumpOp} (abstract) to Vector{JumpOp{Matrix{T}}} (concrete) for zero-allocation vector access"
  - "Per-trajectory Xoshiro(actual_seed + traj_id) seeding in both serial and threaded paths for reproducibility"
  - "Precomputed scaled_prefactor stored in framework to avoid accessing abstract config/precomputed_data in hot loop"

patterns-established:
  - "Function barrier for Any/abstract-typed struct fields: pass field as argument to force specialization"
  - "Allocation test wrapper: @allocated in a function (not top-level scope) for accurate measurement"
  - "Chunk-based parallel execution: _partition_trajectories + per-task workspace + @sync/@spawn"

# Metrics
duration: 27min
completed: 2026-02-16
---

# Phase 13 Plan 01: Multi-Threaded Trajectory Engine Summary

**Zero-allocation step_along_trajectory! with @sync/@spawn threaded dispatch, per-task Xoshiro RNG, and BLAS thread control for parallel trajectory sampling**

## Performance

- **Duration:** 27 min
- **Started:** 2026-02-16T04:46:59Z
- **Completed:** 2026-02-16T05:14:03Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments
- step_along_trajectory! verified as zero-allocation on the hot path (was 1376 bytes before fix)
- Multi-threaded trajectory engine with automatic dispatch when ntraj>1 and nthreads()>1
- BLAS.set_num_threads(1) with try/finally restore prevents oversubscription
- Both no-observable and observable code paths support threaded execution
- Serial path uses consistent per-trajectory Xoshiro seeding matching threaded path

## Task Commits

Each task was committed atomically:

1. **Task 1: Add allocation regression test for step_along_trajectory!** - `0205a17` (test + fix)
2. **Task 2: Implement multi-threaded trajectory engine** - `0bb37a2` (feat)

## Files Created/Modified
- `src/trajectories.jl` - Added concrete-typed hot-path fields to TrajectoryFramework (F,P type params for transition function and NUFFT prefactors), _partition_trajectories helper, _run_chunk_no_obs! and _run_chunk_with_obs! chunk workers, threaded dispatch in run_trajectories, per-trajectory Xoshiro seeding
- `test/test_allocation.jl` - Added step_along_trajectory! zero-allocation regression test with function-wrapped measurement

## Decisions Made
- Added type parameters F and P to TrajectoryFramework{T,D,F,P} to store transition function and NUFFT prefactors with concrete types, eliminating dynamic dispatch in the hot loop
- Changed jumps vector from Vector{JumpOp} (abstract element type) to Vector{JumpOp{Matrix{T}}} (concrete) to prevent allocation on element access
- Per-trajectory seeding with Xoshiro(actual_seed + traj_id) replaces the old shared-RNG serial path, ensuring serial and threaded paths produce equivalent results (up to FP reduction order)
- Allocation test wrapped in a helper function to avoid Julia's top-level scope @allocated measurement artifacts

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed 1376 bytes allocation in step_along_trajectory! hot path**
- **Found during:** Task 1 (allocation test)
- **Issue:** step_along_trajectory! allocated 1376 bytes per call due to accessing abstract-typed fields fw.config (AbstractThermalizeConfig{D}), fw.precomputed_data (Any), and fw.jumps (Vector{JumpOp} with abstract element type)
- **Fix:** Added concrete-typed fields to TrajectoryFramework: scaled_prefactor::Float64, sigma::Float64, transition::F, energy_labels::Vector{Float64}, oft_nufft_prefactors::P. Changed jumps to Vector{JumpOp{Matrix{T}}}. Refactored step_along_trajectory! to read from concrete fields instead of abstract ones.
- **Files modified:** src/trajectories.jl
- **Verification:** Allocation test confirms 0 bytes; all 247 tests pass
- **Committed in:** 0205a17 (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (1 bug)
**Impact on plan:** Essential fix -- without zero-allocation step function, parallel scaling would be destroyed by GC contention. No scope creep.

## Issues Encountered
None -- the allocation fix was anticipated by the plan ("If the test fails because step_along_trajectory! allocates, investigate and fix the allocation source before proceeding").

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Multi-threaded trajectory engine complete, ready for Plan 02 (thread safety and performance testing)
- TrajectoryFramework now has concrete-typed hot-path fields, future plans should use fw.transition, fw.energy_labels, fw.oft_nufft_prefactors instead of fw.precomputed_data
- Serial path seeding changed from shared RNG to per-trajectory Xoshiro(seed+traj_id) -- existing seed-based tests updated and passing

## Self-Check: PASSED

All files and commits verified:
- src/trajectories.jl: FOUND
- test/test_allocation.jl: FOUND
- Commit 0205a17: FOUND
- Commit 0bb37a2: FOUND
- 13-01-SUMMARY.md: FOUND

---
*Phase: 13-multi-threaded-trajectory-engine*
*Completed: 2026-02-16*
