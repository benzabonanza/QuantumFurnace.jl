---
phase: 03-dm-reference-test-suite
plan: 01
subsystem: testing
tags: [lindbladian, detailed-balance, gibbs-state, domain-hierarchy, eigendecomposition]

# Dependency graph
requires:
  - phase: 01-foundation-and-compilation
    provides: "test infrastructure (test_helpers.jl, HamHam, JumpOp fixtures)"
  - phase: 02-trajectory-bug-fixes
    provides: "corrected construct_lindbladian with U_B ordering fix and PSD guard"
provides:
  - "3-qubit test fixture (make_small_test_system, SMALL_HAM, SMALL_JUMPS, SMALL_GIBBS)"
  - "DMTST-01: BohrDomain detailed balance verification (Gibbs as fixed point)"
  - "DMTST-02: Domain error hierarchy verification (bohr <= energy <= time <= trotter)"
  - "Lindbladian fixed-point extraction pattern via eigen() for small dense matrices"
affects: [03-02-PLAN, 03-03-PLAN, 04-trajectory-validation]

# Tech tracking
tech-stack:
  added: []
  patterns: ["eigen() for small Liouvillian eigendecomposition instead of Arpack eigs", "QuantumFurnace.trace_distance_h for unexported function access"]

key-files:
  created: [test/test_dm_detailed_balance.jl]
  modified: [test/test_helpers.jl, test/runtests.jl]

key-decisions:
  - "Use eigen() instead of Arpack eigs for 64x64 and 256x256 Liouvillians (small dense, deterministic)"
  - "Use QuantumFurnace.trace_distance_h since function is not exported from module"
  - "Hierarchy tolerance 1e-12 for numerical noise in domain distance comparisons"

patterns-established:
  - "Lindbladian fixed-point extraction: eigen -> argmin |Re(lambda)| -> reshape -> Hermitianize -> normalize"
  - "3-qubit test system available via SMALL_HAM/SMALL_JUMPS/SMALL_GIBBS constants"

# Metrics
duration: 5min
completed: 2026-02-14
---

# Phase 03 Plan 01: Detailed Balance and Domain Hierarchy Summary

**BohrDomain Lindbladian fixed point matches Gibbs state to 1.6e-15 trace distance; domain error hierarchy bohr <= energy <= time <= trotter verified on 4-qubit system**

## Performance

- **Duration:** 5 min
- **Started:** 2026-02-14T09:32:13Z
- **Completed:** 2026-02-14T09:37:44Z
- **Tasks:** 1
- **Files modified:** 3

## Accomplishments
- DMTST-01: BohrDomain with coherent term B produces exact Gibbs state as fixed point (trace distance 1.6e-15, threshold 1e-10) on 3-qubit Heisenberg system
- DMTST-02: Domain error hierarchy verified: bohr (3.4e-16) <= energy (2.3e-14) <= time (6.3e-14) <= trotter (0.76) on 4-qubit system
- Added make_small_test_system() 3-qubit fixture with 9 Pauli jump operators for future test reuse
- All 49 tests pass (45 existing + 4 new), zero regressions

## Task Commits

Each task was committed atomically:

1. **Task 1: Add 3-qubit fixture and implement detailed balance + hierarchy tests** - `79ccb1a` (test)

**Plan metadata:** pending (docs: complete plan)

## Files Created/Modified
- `test/test_helpers.jl` - Added make_small_test_system() factory and SMALL_* constants for 3-qubit system
- `test/test_dm_detailed_balance.jl` - DMTST-01 (Bohr detailed balance) and DMTST-02 (domain error hierarchy) testsets
- `test/runtests.jl` - Added include for test_dm_detailed_balance.jl

## Decisions Made
- Used `eigen()` from LinearAlgebra instead of Arpack `eigs` for eigendecomposition of small dense Liouvillians (64x64 for 3-qubit, 256x256 for 4-qubit). This is deterministic and avoids Arpack convergence issues.
- Qualified `trace_distance_h` as `QuantumFurnace.trace_distance_h` since the function is not exported from the module (discovered during testing).
- Used `<= ... + 1e-12` tolerance for hierarchy comparisons to allow for tiny numerical noise in domain distance ordering.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Qualified unexported function trace_distance_h**
- **Found during:** Task 1 (test execution)
- **Issue:** `trace_distance_h` is defined in `qi_tools.jl` but not exported from QuantumFurnace module. Test file could not resolve the function.
- **Fix:** Used `QuantumFurnace.trace_distance_h` qualified name in test_dm_detailed_balance.jl
- **Files modified:** test/test_dm_detailed_balance.jl
- **Verification:** Tests pass with qualified call
- **Committed in:** 79ccb1a (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (1 blocking)
**Impact on plan:** Minimal -- function exists but is not in the public API export list. Qualified access is the standard Julia pattern for internal functions.

## Issues Encountered
None beyond the deviation documented above.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- DM reference ground truth established for BohrDomain (exact) and all four domains (hierarchy)
- 3-qubit fixture available for Plan 03-02 (DM step error scaling) and Plan 03-03 (Aqua)
- Fixed-point extraction pattern established for reuse in future DM tests

## Self-Check: PASSED

All files verified present, commit 79ccb1a confirmed in history, all content markers (make_small_test_system, DMTST-01, DMTST-02, runtests include) verified.

---
*Phase: 03-dm-reference-test-suite*
*Completed: 2026-02-14*
