---
phase: 21-fix-test-errors-after-removing-transform
plan: 01
subsystem: testing
tags: [exports, cleanup, allocation-tests]

# Dependency graph
requires:
  - phase: 19-logic-simplification
    provides: "transform_jumps_to_basis removed from qi_tools.jl, TEST_TROTTER_JUMPS constant in test_helpers.jl"
provides:
  - "Clean module exports with no dangling references"
  - "Allocation tests using pre-built TEST_TROTTER_JUMPS"
affects: []

# Tech tracking
tech-stack:
  added: []
  patterns: []

key-files:
  created: []
  modified:
    - src/QuantumFurnace.jl
    - test/test_allocation.jl

key-decisions:
  - "Used TEST_TROTTER_JUMPS from test_helpers.jl (identical transform, already available)"

patterns-established: []

# Metrics
duration: 2min
completed: 2026-02-16
---

# Quick Task 21: Fix Test Errors After Removing transform_jumps_to_basis

**Removed dangling export and test references to deleted transform_jumps_to_basis function, restoring all 539 tests to passing**

## Performance

- **Duration:** 2 min
- **Started:** 2026-02-16T15:17:36Z
- **Completed:** 2026-02-16T15:20:01Z
- **Tasks:** 1
- **Files modified:** 2

## Accomplishments
- Removed `transform_jumps_to_basis` from module export list in `src/QuantumFurnace.jl`
- Removed `transform_jumps_to_basis` from test import in `test/test_allocation.jl`
- Replaced `transform_jumps_to_basis()` call with pre-built `TEST_TROTTER_JUMPS` constant
- All 539 tests pass (0 failures, 0 errors)

## Task Commits

Each task was committed atomically:

1. **Task 1: Remove transform_jumps_to_basis from exports and fix allocation test** - `b2e4123` (fix)

## Files Created/Modified
- `src/QuantumFurnace.jl` - Removed `transform_jumps_to_basis` from QI Tools export block
- `test/test_allocation.jl` - Removed import and replaced function call with TEST_TROTTER_JUMPS constant

## Decisions Made
- Used `TEST_TROTTER_JUMPS` from `test/test_helpers.jl` which is constructed via `make_test_system(; trotter=TEST_TROTTER).jumps` -- identical transform to the removed function call

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Codebase clean: no remaining references to `transform_jumps_to_basis` in src/ or test/
- All 539 tests passing

## Self-Check: PASSED

- All modified files exist on disk
- Commit b2e4123 verified in git log
- No remaining references to transform_jumps_to_basis in src/ or test/
- 539/539 tests pass

---
*Quick Task: 21-fix-test-errors-after-removing-transform*
*Completed: 2026-02-16*
