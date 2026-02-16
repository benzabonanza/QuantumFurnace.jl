---
phase: 19-logic-simplification
plan: 02
subsystem: quantum-simulation
tags: [trajectory, call-chain, framework, convergence, threading]

# Dependency graph
requires:
  - "19-01: Domain-aware JumpOp construction at the source"
provides:
  - "Flattened 3-level trajectory call chain (public API -> batch execution -> step_along_trajectory!)"
  - "_build_framework_and_seed for one-time framework setup"
  - "_run_batch_no_obs! for shared serial/threaded batch execution"
  - "TrajectoryFramework built once in convergence/adaptive runners (not per batch)"
affects: [19-03]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "One-time framework building via _build_framework_and_seed, reused across batches"
    - "Shared _run_batch_no_obs! handles serial vs threaded dispatch for all callers"
    - "Step loop inlined into chunk runners (no wrapper function)"

key-files:
  created: []
  modified:
    - "src/trajectories.jl"
    - "src/convergence.jl"

key-decisions:
  - "Inline _evolve_along_trajectory! into chunk runners rather than keeping as thin wrapper"
  - "Extract _build_framework_and_seed as shared setup function for framework + seed"
  - "Create _run_batch_no_obs! that encapsulates serial/threaded dispatch for batch execution"
  - "Convergence runners call _run_batch_no_obs! directly instead of run_trajectories"

patterns-established:
  - "_build_framework_and_seed is the single entry point for TrajectoryFramework construction"
  - "_run_batch_no_obs! is the single entry point for batch trajectory execution without observables"

# Metrics
duration: 7min
completed: 2026-02-16
---

# Phase 19 Plan 02: Flatten Trajectory Call Chain Summary

**3-level call chain with _run_batch_no_obs! shared across run_trajectories, convergence, and adaptive runners -- framework built once, not per batch**

## Performance

- **Duration:** 7 min
- **Started:** 2026-02-16T14:26:20Z
- **Completed:** 2026-02-16T14:33:46Z
- **Tasks:** 2/2
- **Files modified:** 2

## Accomplishments
- Deleted _evolve_along_trajectory! (trivial 13-line wrapper), step loop inlined into chunk runners
- Extracted _build_framework_and_seed for one-time setup reuse across convergence/adaptive batches
- Created _run_batch_no_obs! encapsulating serial/threaded dispatch, used by all three public runners
- Convergence runners (run_trajectories_convergence, run_trajectories_adaptive) no longer rebuild TrajectoryFramework per batch
- Call chain reduced from 5 levels to 3: public API -> _run_batch_no_obs! -> step_along_trajectory!

## Task Commits

Each task was committed atomically:

1. **Task 1: Inline _evolve_along_trajectory! into chunk runners and extract framework building** - `aac622d` (refactor)
2. **Task 2: Lift framework building out of convergence runners** - `8906084` (refactor)

## Files Created/Modified
- `src/trajectories.jl` - Deleted _evolve_along_trajectory!, inlined step loop into _run_chunk_no_obs! and serial path, added _build_framework_and_seed and _run_batch_no_obs!, refactored run_trajectories to use both
- `src/convergence.jl` - Refactored run_trajectories_convergence and run_trajectories_adaptive to build framework once via _build_framework_and_seed and call _run_batch_no_obs! per batch

## Decisions Made
- Inlined _evolve_along_trajectory! rather than keeping it as a thinner wrapper -- the function added no abstraction value and removing it eliminates one call level
- Created _run_batch_no_obs! as a separate function (not merging into _run_chunk_no_obs!) because it handles the serial vs threaded dispatch decision
- Convergence runners call _run_batch_no_obs! directly instead of going through run_trajectories, avoiding redundant config validation and framework building per batch

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- Flattened call chain complete, all 539 tests pass
- Public API signatures unchanged (run_trajectories, run_trajectories_convergence, run_trajectories_adaptive)
- Ready for plan 03 (result struct simplification)

## Self-Check: PASSED

All 2 modified files verified present. All 2 task commits verified in git log. 539/539 tests pass.

---
*Phase: 19-logic-simplification*
*Completed: 2026-02-16*
