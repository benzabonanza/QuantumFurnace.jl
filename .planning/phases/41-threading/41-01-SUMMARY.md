---
phase: 41-threading
plan: 01
subsystem: threading
tags: [blas, threading, density-matrix, try-finally, linearalgebra]

# Dependency graph
requires:
  - phase: 39-precomputation
    provides: "per-jump precomputed channels used in run_thermalize hot loop"
  - phase: 40-save-every
    provides: "save_every gating in run_thermalize (loop structure wrapping target)"
provides:
  - "BLAS try/finally save/restore in run_thermalize with explicit multi-threaded enablement"
  - "DM thermalization threading tests (restoration, agreement, isolation)"
affects: [41-02-omega-threading, mixing-time-estimation]

# Tech tracking
tech-stack:
  added: []
  patterns: ["BLAS save/restore try/finally in DM path (opposite of trajectory: sets BLAS threads high instead of 1)"]

key-files:
  created: []
  modified:
    - src/furnace.jl
    - test/test_threading.jl

key-decisions:
  - "Explicit BLAS.set_num_threads(Threads.nthreads()) inside try block -- enables multi-threaded BLAS even if caller reduced it"
  - "Only hot loop wrapped in try/finally -- precomputation and result construction remain outside"

patterns-established:
  - "DM BLAS pattern: save -> set high -> try loop finally restore (opposite of trajectory pattern which sets to 1)"

# Metrics
duration: 7min
completed: 2026-03-01
---

# Phase 41 Plan 01: BLAS Threading Summary

**BLAS try/finally wrapping run_thermalize hot loop with explicit multi-threaded BLAS enablement and DM-specific threading tests**

## Performance

- **Duration:** 7 min
- **Started:** 2026-03-01T11:44:33Z
- **Completed:** 2026-03-01T11:51:21Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments
- run_thermalize hot loop wrapped in BLAS save/restore with explicit multi-threaded BLAS enablement (THREAD-03)
- BLAS thread count guaranteed restored after run_thermalize returns, even on error (THREAD-05)
- DM BLAS thread restoration verified across Energy/Time/Trotter domains
- Serial-BLAS and multi-threaded-BLAS agreement verified within atol=1e-10
- DM-trajectory BLAS isolation confirmed (no thread count leak between run_thermalize and run_trajectories)
- Full test suite: 1170 tests pass (7 new, 0 regressions)

## Task Commits

Each task was committed atomically:

1. **Task 1: Add BLAS try/finally to run_thermalize** - `a4c898c` (feat)
2. **Task 2: Add DM thermalization threading tests** - `45dbcfe` (test)

## Files Created/Modified
- `src/furnace.jl` - Added BLAS save/restore try/finally around run_thermalize hot loop with Threads.nthreads() enablement
- `test/test_threading.jl` - Added 3 DM threading test blocks: restoration, serial-threaded agreement, DM-trajectory isolation

## Decisions Made
- Explicit BLAS.set_num_threads(Threads.nthreads()) inside try block ensures multi-threaded BLAS even if caller has reduced BLAS threads (e.g., trajectory engine sets BLAS to 1)
- Only the hot loop is wrapped in try/finally -- precomputation and result construction stay outside (matches plan spec)
- DM BLAS pattern is the inverse of trajectory pattern: trajectory sets BLAS=1 for Julia-level parallelism, DM sets BLAS=nthreads for BLAS-level parallelism

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed TrotterDomain test to use correct trotter-basis jumps**
- **Found during:** Task 2 (DM threading tests)
- **Issue:** Plan used N3_JUMPS (hamiltonian eigenbasis) for TrotterDomain test, but TrotterDomain requires N3_TROTTER_JUMPS (trotter eigenbasis) and N3_TROTTER object
- **Fix:** Used N3_TROTTER_JUMPS and passed N3_TROTTER as the trotter argument in the TrotterDomain iteration of the restoration test
- **Files modified:** test/test_threading.jl
- **Verification:** All 1170 tests pass including TrotterDomain restoration test
- **Committed in:** 45dbcfe (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (1 bug)
**Impact on plan:** Auto-fix necessary for correctness. Without it, TrotterDomain test would assertion-error. No scope creep.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- BLAS threading for DM path complete and tested
- Ready for Plan 02 (omega-loop threading or additional threading optimizations)
- run_thermalize now automatically benefits from multi-threaded BLAS for dim >= 64 systems

## Self-Check: PASSED

All files exist, all commits verified.

---
*Phase: 41-threading*
*Completed: 2026-03-01*
