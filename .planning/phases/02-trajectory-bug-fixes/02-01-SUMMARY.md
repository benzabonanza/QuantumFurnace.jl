---
phase: 02-trajectory-bug-fixes
plan: 01
subsystem: simulation
tags: [trajectories, kraus, cptp, eigendecomposition, psd-guard, normalization]

# Dependency graph
requires:
  - phase: 01-foundation-and-compilation
    provides: "Compilation fixes, test infrastructure (TEST_HAM, TEST_JUMPS, KrausScratch, DIM, TEST_DELTA)"
provides:
  - "Fixed step_along_trajectory! with correct U_B ordering (TFIX-02)"
  - "Normalization warning in both trajectory variants (TFIX-03)"
  - "PSD-guarded build_trajectoryframework via eigendecomposition (TFIX-04)"
  - "Cross-checked channel structure documentation (TFIX-05)"
  - "make_test_trotter() and TEST_TROTTER for TrotterDomain tests"
  - "20 single-step trajectory fix tests"
affects: [02-02 (CPTP verification), phase-04 (trajectory validation)]

# Tech tracking
tech-stack:
  added: []
  patterns: [eigendecomposition-for-psd-guard, normalization-check-via-warn]

key-files:
  created:
    - test/test_trajectory_fixes.jl
  modified:
    - src/trajectories.jl
    - test/test_helpers.jl
    - test/runtests.jl

key-decisions:
  - "EnergyDomain had U_B ordering bug; Time/Trotter variant was already correct"
  - "Used TrottTrott() constructor directly instead of nonexistent create_trotter() function"
  - "Used Random.seed!() for test reproducibility since step_along_trajectory! uses global rand()"

patterns-established:
  - "PSD guard via eigen + clamp: replace cholesky!(check=false) with eigen + max(eigenvalues, 0)"
  - "TFIX comment convention: # TFIX-NN: description"

# Metrics
duration: 5min
completed: 2026-02-14
---

# Phase 2 Plan 1: Trajectory Bug Fixes Summary

**Fixed U_B ordering bug in EnergyDomain step_along_trajectory!, added normalization warning and PSD eigendecomposition guard, with 20 single-step verification tests**

## Performance

- **Duration:** 5 min
- **Started:** 2026-02-14T08:55:29Z
- **Completed:** 2026-02-14T09:00:42Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments
- Fixed critical U_B ordering bug in EnergyDomain step_along_trajectory! (TFIX-02) -- U_B was applied after probability computation, causing inconsistent state
- Added normalization warning (@warn) to both trajectory variants when total_weight deviates from 1.0 by more than 1e-6 (TFIX-03)
- Replaced Cholesky decomposition with eigendecomposition PSD guard in build_trajectoryframework (TFIX-04) -- negative eigenvalues clamped to zero silently
- Added channel structure documentation cross-referencing Chen 2023 Theorem III.1 (TFIX-05)
- Created 20 new tests covering all four TFIX requirements across domains (42 total tests, up from 22)
- Added make_test_trotter() helper and TEST_TROTTER constant for TrotterDomain tests

## Task Commits

Each task was committed atomically:

1. **Task 1: Fix trajectory bugs in src/trajectories.jl** - `26c37ed` (fix)
2. **Task 2: Create single-step bug fix tests** - `9689c35` (test)

## Files Created/Modified
- `src/trajectories.jl` - Fixed U_B ordering (EnergyDomain), added normalization warning (both variants), replaced Cholesky with eigen PSD guard, added TFIX-05 channel documentation
- `test/test_trajectory_fixes.jl` - 4 testsets: TFIX-02 (U_B ordering), TFIX-03 (normalization), TFIX-04 (PSD guard), TFIX-05 (probability decomposition)
- `test/test_helpers.jl` - Added make_test_trotter() and TEST_TROTTER constant
- `test/runtests.jl` - Added include for test_trajectory_fixes.jl

## Decisions Made
- **EnergyDomain only had the bug:** The Time/Trotter variant already had U_B at the top of step_along_trajectory!. Only the EnergyDomain variant needed the U_B block moved. Both variants received the TFIX-02 comment and TFIX-03 normalization warning for consistency.
- **Used TrottTrott() constructor directly:** The plan referenced `create_trotter()` which is exported but never defined. Used `TrottTrott(TEST_HAM, T0, NUM_TROTTER_STEPS_PER_T0)` instead.
- **Used Random.seed!() instead of StableRNG:** step_along_trajectory! uses `rand()` (global RNG) with no rng keyword argument. Used `Random.seed!()` for test reproducibility instead of the plan's `StableRNG(42)` approach.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] create_trotter() function does not exist**
- **Found during:** Task 1 (make_test_trotter helper)
- **Issue:** Plan specified `create_trotter(TEST_HAM, T0, NUM_TROTTER_STEPS_PER_T0)` but `create_trotter` is exported but never defined in the codebase
- **Fix:** Used `TrottTrott(TEST_HAM, T0, NUM_TROTTER_STEPS_PER_T0)` constructor directly
- **Files modified:** test/test_helpers.jl
- **Verification:** Tests compile and pass
- **Committed in:** 26c37ed (Task 1 commit)

**2. [Rule 3 - Blocking] step_along_trajectory! has no rng keyword argument**
- **Found during:** Task 2 (test creation)
- **Issue:** Plan's test template used `step_along_trajectory!(psi, fw; rng=rng)` with `StableRNG`, but the function uses global `rand()` with no rng parameter
- **Fix:** Used `Random.seed!(42)` before each step call for reproducibility
- **Files modified:** test/test_trajectory_fixes.jl
- **Verification:** Tests pass deterministically
- **Committed in:** 9689c35 (Task 2 commit)

**3. [Rule 1 - Bug] Time/Trotter variant U_B ordering was already correct**
- **Found during:** Task 1 (TFIX-02 fix)
- **Issue:** Plan stated BOTH variants had the U_B ordering bug, but the Time/Trotter variant (lines 447-451) already had U_B at the top before probability computation. Only the EnergyDomain variant (lines 635-639) had U_B in the wrong position.
- **Fix:** Only moved U_B in the EnergyDomain variant. Added TFIX-02 comment to both variants for documentation consistency.
- **Files modified:** src/trajectories.jl
- **Verification:** Both variants now have U_B at top; all tests pass
- **Committed in:** 26c37ed (Task 1 commit)

---

**Total deviations:** 3 auto-fixed (1 bug clarification, 2 blocking)
**Impact on plan:** All deviations were necessary for correctness. No scope creep. The key insight is that the Time/Trotter variant was already correct -- only EnergyDomain needed the U_B fix.

## Issues Encountered
None -- all tasks completed without unexpected problems beyond the documented deviations.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Trajectory bugs fixed; ready for CPTP verification (Plan 02)
- TEST_TROTTER available for TrotterDomain CPTP tests
- All 42 tests pass (22 Phase 1 + 20 Phase 2 Plan 1)
- build_trajectoryframework now produces reliable U_residual via PSD-guarded eigendecomposition

## Self-Check: PASSED

All files verified present, all commit hashes found in git log.

---
*Phase: 02-trajectory-bug-fixes*
*Completed: 2026-02-14*
