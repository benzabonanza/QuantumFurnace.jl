---
phase: 12-workspace-refactor
plan: 01
subsystem: trajectories
tags: [julia, workspace, rng, thread-safety, struct-refactor]

# Dependency graph
requires:
  - phase: 02-trajectory-bug-fixes
    provides: "TrajectoryFramework with per-operator Kraus stepping"
  - phase: 11-allocation-perf
    provides: "Optimized allocation-free step paths"
provides:
  - "Read-only TrajectoryFramework (no ws field)"
  - "TrajectoryWorkspace with rho_acc accumulator"
  - "TrajectoryResult struct for structured returns"
  - "Explicit ws and rng parameters on step_along_trajectory!"
  - "Explicit ws and rng parameters on _evolve_along_trajectory!"
  - "Seed capture in run_trajectories for reproducibility"
affects: [13-thread-pool, 14-batch-convergence, 15-seeding]

# Tech tracking
tech-stack:
  added: []
  patterns: [explicit-rng-threading, workspace-separation, structured-result-type]

key-files:
  created: []
  modified:
    - src/trajectories.jl
    - src/QuantumFurnace.jl
    - test/test_trajectory_fixes.jl
    - test/test_regression.jl
    - test/trajectory_validation/run_trajectory_validation.jl
    - test/trajectory_validation/run_convergence_tests.jl

key-decisions:
  - "TrajectoryResult includes optional times and measurements_mean fields for observable path"
  - "store_states kwarg dropped from run_trajectories (no tests use it, TrajectoryResult does not carry states)"
  - "Seed auto-generated via Int(rand(RandomDevice(), UInt64) >> 1) for positive Int values"
  - "TrajectoryWorkspace(fw) convenience constructor placed after TrajectoryFramework definition for Julia ordering"

patterns-established:
  - "Explicit RNG: all rand() calls in step functions use rand(rng, ...) for thread safety"
  - "Workspace separation: framework is read-only, workspace is mutable, passed explicitly"
  - "Structured returns: run_trajectories returns TrajectoryResult instead of NamedTuple"

# Metrics
duration: 7min
completed: 2026-02-15
---

# Phase 12 Plan 01: Workspace Refactor Summary

**Separated mutable TrajectoryWorkspace from read-only TrajectoryFramework with explicit RNG threading and TrajectoryResult return type**

## Performance

- **Duration:** 7 min
- **Started:** 2026-02-15T19:16:10Z
- **Completed:** 2026-02-15T19:23:45Z
- **Tasks:** 2
- **Files modified:** 6

## Accomplishments
- TrajectoryFramework is now fully read-only during stepping (ws field removed)
- All 4 rand() call sites in step functions now use explicit rng::AbstractRNG
- TrajectoryResult struct provides structured, reproducible returns with seed capture
- All 231 existing tests pass with updated call signatures

## Task Commits

Each task was committed atomically:

1. **Task 1: Refactor structs and step function signatures** - `10dabd9` (feat)
2. **Task 2: Update exports and verify compilation** - `025c97b` (feat)

## Files Created/Modified
- `src/trajectories.jl` - Removed ws from TrajectoryFramework, added rho_acc to TrajectoryWorkspace, added TrajectoryResult struct, refactored all step/evolve/run functions to take explicit ws and rng
- `src/QuantumFurnace.jl` - Added TrajectoryResult to exports
- `test/test_trajectory_fixes.jl` - Updated 3 step_along_trajectory! call sites to 4-arg form
- `test/test_regression.jl` - Updated 2 trajectory regression loops to use explicit ws/rng
- `test/trajectory_validation/run_trajectory_validation.jl` - Updated single-step crossval and run_trajectories calls
- `test/trajectory_validation/run_convergence_tests.jl` - Updated reference computation and batch loops

## Decisions Made
- TrajectoryResult includes optional `times` and `measurements_mean` fields (Union{Nothing, ...}) to keep a single return type for both observable and non-observable paths
- Dropped `store_states` kwarg entirely from run_trajectories since no tests use it and TrajectoryResult does not carry per-trajectory state vectors
- Used `Int(rand(Random.RandomDevice(), UInt64) >> 1)` for seed generation to ensure positive Int values
- Placed TrajectoryWorkspace(fw) convenience constructor after TrajectoryFramework struct definition to satisfy Julia's top-to-bottom evaluation order

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Updated test call sites to new 4-argument signatures**
- **Found during:** Task 1 (Refactor structs and step function signatures)
- **Issue:** Plan listed only src/trajectories.jl for Task 1, but 4 test files had 8 call sites using the old 2-argument step_along_trajectory!(psi, fw) signature that would fail at runtime
- **Fix:** Updated all test call sites to create local TrajectoryWorkspace and Xoshiro RNG, then call step_along_trajectory!(psi, fw, ws, rng)
- **Files modified:** test/test_trajectory_fixes.jl, test/test_regression.jl, test/trajectory_validation/run_trajectory_validation.jl, test/trajectory_validation/run_convergence_tests.jl
- **Verification:** All 231 tests pass
- **Committed in:** 10dabd9 (Task 1 commit)

**2. [Rule 1 - Bug] Fixed TrajectoryWorkspace(fw) constructor ordering**
- **Found during:** Task 1 (Refactor structs and step function signatures)
- **Issue:** Convenience constructor TrajectoryWorkspace(fw::TrajectoryFramework{T}) was placed before TrajectoryFramework struct definition, causing UndefVarError at module load time
- **Fix:** Moved constructor to after TrajectoryFramework struct definition
- **Verification:** Module compiles successfully
- **Committed in:** 10dabd9 (Task 1 commit)

---

**Total deviations:** 2 auto-fixed (1 blocking, 1 bug)
**Impact on plan:** Both auto-fixes were necessary for correct compilation and test passing. No scope creep.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- TrajectoryFramework is now read-only and safe for concurrent access from multiple threads
- TrajectoryWorkspace can be independently allocated per-thread
- Explicit RNG parameters enable reproducible per-thread seeding
- Ready for Phase 12 Plan 02 (test updates) and Phase 13 (thread pool implementation)

---
*Phase: 12-workspace-refactor*
*Completed: 2026-02-15*
