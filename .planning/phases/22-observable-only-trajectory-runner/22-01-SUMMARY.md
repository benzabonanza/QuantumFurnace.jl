---
phase: 22-observable-only-trajectory-runner
plan: 01
subsystem: trajectories
tags: [quantum-trajectories, observables, density-matrix, threading, cross-validation]

# Dependency graph
requires:
  - phase: 20-observable-infrastructure
    provides: "Observable builders (build_gap_estimation_observables, build_total_magnetization)"
  - phase: 13-trajectory-threading
    provides: "Multi-threaded trajectory infrastructure (_partition_trajectories, per-task workspaces, BLAS thread management)"
provides:
  - "run_observable_trajectories function for time-resolved observable decay curves without DM overhead"
  - "ObservableTrajectoryResult struct with optional rho_mean"
  - "_run_chunk_obs_only! inner function for observable-only threaded chunks"
affects: [23-gap-estimation-api, 24-mixing-time-estimation]

# Tech tracking
tech-stack:
  added: []
  patterns: ["Observable-only trajectory runner with optional DM reconstruction via reconstruct_dm flag"]

key-files:
  created:
    - test/test_observable_trajectories.jl
  modified:
    - src/trajectories.jl
    - src/QuantumFurnace.jl
    - test/runtests.jl

key-decisions:
  - "ObservableTrajectoryResult placed in trajectories.jl alongside TrajectoryResult (not structs.jl) to keep related types together"
  - "Inner/outer constructor pattern to satisfy Aqua unbound type parameter checks when rho_mean=nothing"
  - "reconstruct_dm=true reuses existing _run_chunk_with_obs! rather than duplicating DM accumulation logic"

patterns-established:
  - "Observable-only runner: dedicated function for obs-only measurements, separate from run_trajectories"
  - "Optional DM reconstruction via reconstruct_dm kwarg dispatching between _run_chunk_obs_only! and _run_chunk_with_obs!"

# Metrics
duration: 8min
completed: 2026-02-17
---

# Phase 22 Plan 01: Observable-Only Trajectory Runner Summary

**run_observable_trajectories with optional DM reconstruction, bitwise cross-validated against existing runner for serial and threaded paths**

## Performance

- **Duration:** 8 min
- **Started:** 2026-02-17T09:06:15Z
- **Completed:** 2026-02-17T09:14:16Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments
- Implemented `run_observable_trajectories` that measures time-resolved `<O>(t)` without per-trajectory density matrix reconstruction overhead
- Added `ObservableTrajectoryResult` struct with optional `rho_mean` field (nothing when `reconstruct_dm=false`)
- Cross-validated bitwise agreement with existing `run_trajectories` for observable means, times, and density matrices
- Verified deterministic seeding, serial/threaded paths, and single-trajectory correctness (18 tests total)
- All 633 tests pass (615 existing + 18 new)

## Task Commits

Each task was committed atomically:

1. **Task 1: Implement ObservableTrajectoryResult, _run_chunk_obs_only!, and run_observable_trajectories** - `01b65fd` (feat)
2. **Task 2: Cross-validation and correctness tests** - `c7b5227` (test)

## Files Created/Modified
- `src/trajectories.jl` - Added ObservableTrajectoryResult struct, _run_chunk_obs_only! function, run_observable_trajectories function
- `src/QuantumFurnace.jl` - Added exports for ObservableTrajectoryResult and run_observable_trajectories
- `test/test_observable_trajectories.jl` - 6 testsets: basic run, bitwise cross-validation, reconstruct_dm, determinism, different seeds, serial path
- `test/runtests.jl` - Added include for new test file

## Decisions Made
- **ObservableTrajectoryResult location:** Placed in `trajectories.jl` (not `structs.jl` as plan suggested) because `TrajectoryResult` lives there -- keeps related types together
- **Constructor pattern:** Added explicit inner constructor + two outer constructors (Matrix{T} and Nothing) to pass Aqua's unbound type parameter check
- **reconstruct_dm dispatch:** Uses existing `_run_chunk_with_obs!` when `reconstruct_dm=true`, avoiding code duplication

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed Aqua unbound type parameter failure**
- **Found during:** Task 2 (test verification)
- **Issue:** Default constructor of `ObservableTrajectoryResult` had unbound `T` when `rho_mean=nothing`, causing Aqua "unbound_args" test to fail
- **Fix:** Added explicit inner constructor and two outer constructors (one for `Matrix{T}`, one for `Nothing` defaulting to `ComplexF64`)
- **Files modified:** `src/trajectories.jl`
- **Verification:** All 633 tests pass including Aqua checks
- **Committed in:** `01b65fd` (Task 1 commit, amended)

**2. [Rule 3 - Blocking] Placed struct in trajectories.jl instead of structs.jl**
- **Found during:** Task 1 (implementation)
- **Issue:** Plan specified adding struct to `structs.jl`, but `TrajectoryResult` (the analogous existing type) is actually in `trajectories.jl`
- **Fix:** Added `ObservableTrajectoryResult` to `trajectories.jl` after `TrajectoryResult` for consistency
- **Files modified:** `src/trajectories.jl`
- **Verification:** Compiles correctly, exports resolve
- **Committed in:** `01b65fd` (Task 1 commit)

---

**Total deviations:** 2 auto-fixed (1 bug, 1 blocking)
**Impact on plan:** Both fixes necessary for correctness. No scope creep.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- `run_observable_trajectories` is ready for Phase 23 (Gap Estimation API) to call for time-resolved observable decay curves
- Cross-validation proves bitwise equivalence with existing runner, ensuring no correctness regressions
- Both serial and multi-threaded paths verified

---
*Phase: 22-observable-only-trajectory-runner*
*Completed: 2026-02-17*
