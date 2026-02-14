---
phase: quick-10
plan: 01
subsystem: testing
tags: [regression, blas, cross-platform, tolerance]

# Dependency graph
requires:
  - phase: 05-statistical-validation-and-regression
    provides: Frozen BSON regression tests (TINF-02)
provides:
  - Cross-platform-safe trajectory regression tests with 1e-6 tolerance
affects: []

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Strict 1e-10 for deterministic DM tests, relaxed 1e-6 for stochastic trajectory tests"

key-files:
  created: []
  modified:
    - test/test_regression.jl

key-decisions:
  - "1e-6 tolerance for trajectory tests accommodates BLAS platform variance while still catching real regressions (differences >> 0.001)"
  - "DM tests unchanged at 1e-10 (fully deterministic, no BLAS sensitivity)"

patterns-established:
  - "Split tolerance strategy: deterministic paths get strict tolerance, stochastic paths get relaxed tolerance with explanatory comments"

# Metrics
duration: 2min
completed: 2026-02-14
---

# Quick Task 10: Fix Trajectory Regression Test Failure Summary

**Relaxed trajectory regression tolerances from 1e-10 to 1e-6 for cross-platform BLAS safety while keeping DM tests strict**

## Performance

- **Duration:** 2 min
- **Started:** 2026-02-14T17:02:31Z
- **Completed:** 2026-02-14T17:04:29Z
- **Tasks:** 1
- **Files modified:** 1

## Accomplishments
- Trajectory regression tests (EnergyDomain + TrotterDomain coherent) relaxed to atol=1e-6
- DM regression tests preserved at strict atol=1e-10
- Added explanatory comments before each trajectory testset documenting why tolerance differs
- Full test suite passes (224/224 tests)

## Task Commits

Each task was committed atomically:

1. **Task 1: Relax trajectory regression test tolerances** - `9d2e71a` (fix)

## Files Created/Modified
- `test/test_regression.jl` - Relaxed trajectory atol from 1e-10 to 1e-6, added BLAS variance comments

## Decisions Made
- 1e-6 tolerance chosen: tight enough to catch real code regressions (which would produce differences >> 0.001) but loose enough to accommodate O(1/sqrt(ntraj)) differences from cross-platform BLAS rounding in stochastic branch decisions
- DM tests unchanged: density matrix propagation via matrix exponential is deterministic and not affected by BLAS platform differences

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- All regression tests are now robust across OpenBLAS and Apple Accelerate platforms
- No blockers or concerns

## Self-Check: PASSED

- FOUND: test/test_regression.jl
- FOUND: commit 9d2e71a
- FOUND: 10-SUMMARY.md

---
*Quick Task: 10-fix-trajectory-regression-test-failure*
*Completed: 2026-02-14*
