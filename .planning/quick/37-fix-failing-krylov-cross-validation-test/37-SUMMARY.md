---
phase: 37-fix-failing-krylov-cross-validation-test
plan: 01
subsystem: testing
tags: [krylov, spectral-gap, cross-validation, convergence-order, chen-channel]

# Dependency graph
requires:
  - phase: quick-36
    provides: "Faithful Chen CPTP channel replacing Euler channel"
provides:
  - "All Krylov cross-validation tests passing with correct O(delta) thresholds"
  - "Accurate docstrings for channel-to-Lindbladian eigenvalue conversion"
affects: [krylov-eigsolve, cross-validation]

# Tech tracking
tech-stack:
  added: []
  patterns: []

key-files:
  created: []
  modified:
    - src/krylov_eigsolve.jl
    - test/test_krylov_crossvalidation.jl

key-decisions:
  - "Convergence threshold lowered from >= 1.5 to >= 0.9 to match O(delta) error of faithful Chen channel"
  - "Docstring corrected from 'exact linear formula' to 'first-order approximation' with O(delta) error explanation"

patterns-established: []

# Metrics
duration: 2min
completed: 2026-02-25
---

# Quick Task 37: Fix Failing Krylov Cross-Validation Test Summary

**Corrected XVAL-03 convergence order thresholds from >= 1.5 to >= 0.9 and fixed misleading docstrings after Euler-to-Chen channel migration**

## Performance

- **Duration:** 2 min
- **Started:** 2026-02-25T09:29:21Z
- **Completed:** 2026-02-25T09:31:38Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments
- Fixed 8 failing XVAL-03 tests by lowering convergence order threshold from >= 1.5 to >= 0.9
- Corrected docstring in `krylov_spectral_gap` to describe first-order approximation (not "exact linear formula")
- Updated comments and docstrings in test file to explain O(delta) convergence from faithful Chen channel
- Verified all 16 cross-validation tests pass (XVAL-01, XVAL-02, XVAL-03) with no regressions

## Task Commits

Each task was committed atomically:

1. **Task 1: Fix docstring in krylov_eigsolve.jl** - `01be229` (fix)
2. **Task 2: Fix XVAL-03 test thresholds and comments** - `1036b35` (fix)

## Files Created/Modified
- `src/krylov_eigsolve.jl` - Corrected docstring for `krylov_spectral_gap(config::AbstractThermalizeConfig, ...)` to describe O(delta) approximation error
- `test/test_krylov_crossvalidation.jl` - Updated `run_le_convergence` docstring, XVAL-03 section comments, and all 4 domain thresholds from 1.5 to 0.9

## Decisions Made
- Convergence threshold >= 0.9 chosen because measured orders are ~0.98-1.02 (O(delta) convergence), providing comfortable margin
- Docstring wording "first-order approximation" with explicit `mu = exp(delta * lambda_L) + O(delta^2)` formula chosen for mathematical precision

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- All Krylov cross-validation tests green
- No blockers for future work

---
*Quick Task: 37-fix-failing-krylov-cross-validation-test*
*Completed: 2026-02-25*

## Self-Check: PASSED

All files exist, all commits verified.
