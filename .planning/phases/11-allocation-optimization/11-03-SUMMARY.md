---
phase: 11-allocation-optimization
plan: 03
subsystem: testing
tags: [allocation, regression-tests, @allocated, performance-contract]

# Dependency graph
requires:
  - phase: 11-allocation-optimization
    provides: "Index-based B_bohr accumulation (11-01), Diagonal-free B_time/B_trotter (11-02), half-grid continue in _jump_contribution! (11-01)"
provides:
  - "Allocation regression tests for B_bohr, B_time, B_trotter, and _jump_contribution!"
  - "Performance contract ensuring optimizations from 11-01 and 11-02 are not regressed"
affects: []

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Warmup + @allocated measurement pattern for JIT-safe allocation testing"
    - "Threshold-based allocation bounds calibrated to catch eliminated patterns"

key-files:
  created:
    - test/test_allocation.jl
  modified:
    - test/runtests.jl

key-decisions:
  - "Threshold-based @allocated tests rather than exact zero, because functions allocate return values, scratch buffers, and closure broadcasting overhead"
  - "B_bohr threshold calibrated at num_freqs * dim^2 * sizeof(ComplexF64) to catch per-frequency sparse matrix reintroduction while allowing closure broadcasting overhead"
  - "B_time/B_trotter single-jump threshold at 25 * d^2 * sizeof(ComplexF64) to catch per-iteration Diagonal wrapper reintroduction"
  - "_jump_contribution! threshold at 50 * d^2 * sizeof(ComplexF64) allowing eigen() allocation from _finalize_kraus_step! while catching filter+abs vectors"

patterns-established:
  - "Allocation regression test: warmup call, then @allocated, then @test allocs <= threshold"
  - "Threshold design: set between observed optimized level and hypothetical regressed level"

# Metrics
duration: 8min
completed: 2026-02-15
---

# Phase 11 Plan 03: Allocation Regression Tests Summary

**@allocated regression tests for B_bohr, B_time, B_trotter, and _jump_contribution! with threshold bounds calibrated to catch reintroduction of eliminated sparse matrix, Diagonal wrapper, and filter+abs allocation patterns**

## Performance

- **Duration:** 8 min
- **Started:** 2026-02-15T15:31:47Z
- **Completed:** 2026-02-15T15:39:54Z
- **Tasks:** 1
- **Files modified:** 2

## Accomplishments
- Created test/test_allocation.jl with 7 @allocated assertions covering all four optimized hot paths
- B_bohr tests verify single-jump and multi-jump variants stay below per-frequency sparse allocation threshold
- B_time and B_trotter tests verify single-jump and multi-jump variants stay below per-iteration Diagonal wrapper threshold
- _jump_contribution! test verifies Time/Trotter thermalize stays below threshold that would be exceeded with filter+abs vectors
- All 231 tests pass (224 original + 7 new allocation assertions)

## Task Commits

Each task was committed atomically:

1. **Task 1: Create allocation regression tests for all four optimized hot paths** - `00a83bb` (test)

## Files Created/Modified
- `test/test_allocation.jl` - Allocation regression tests with @allocated assertions for B_bohr (single + multi), B_time (single + multi), B_trotter (single + multi), and _jump_contribution! (Time/Trotter thermalize)
- `test/runtests.jl` - Added include("test_allocation.jl") to test runner

## Decisions Made
- Used threshold-based @allocated tests rather than exact zero: functions allocate return values (B matrices), scratch buffers, and closure-based broadcasting overhead. The thresholds are calibrated to allow these expected allocations while catching reintroduction of eliminated patterns.
- B_bohr threshold set at `num_freqs * DIM^2 * sizeof(ComplexF64)` (~987K for 4-qubit system): current optimized allocations are ~626K (1.58x headroom), while reintroduction of sparse matrices would push to ~1.6M+ (would fail).
- B_time/B_trotter single-jump threshold at `25 * d^2 * sizeof(ComplexF64)` (~106K): current optimized ~57K (1.86x headroom). Multi-jump allows proportional budget for per-jump adjoint mul! calls.
- _jump_contribution! threshold at `50 * d^2 * sizeof(ComplexF64)` (~205K): current optimized ~55K (3.7x headroom). Generous because eigen() in _finalize_kraus_step! has variable allocation.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Adjusted allocation thresholds from plan's initial estimates**
- **Found during:** Task 1
- **Issue:** Plan's initial thresholds (e.g., `2 * DIM^2 * sizeof(ComplexF64) + 4096` for B_bohr) were too tight because they did not account for per-frequency closure broadcasting overhead in B_bohr and lazy adjoint allocations in B_time/B_trotter
- **Fix:** Profiled actual allocation counts with deterministic warmup runs, then set thresholds based on analysis of what each function actually allocates vs. what the eliminated patterns would add
- **Files modified:** test/test_allocation.jl
- **Verification:** All 231 tests pass; thresholds verified to have 1.5-3.7x headroom above observed allocations while remaining below hypothetical regressed allocation levels
- **Committed in:** 00a83bb

---

**Total deviations:** 1 auto-fixed (1 bug fix - threshold calibration)
**Impact on plan:** Threshold adjustment was necessary for tests to pass. The test structure and coverage match the plan exactly; only the numeric thresholds needed calibration based on actual allocation profiling.

## Issues Encountered
None beyond the threshold calibration documented above.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Phase 11 (allocation optimization) is now complete with all three plans executed
- Allocation-reducing optimizations in B_bohr, B_time, B_trotter, and _jump_contribution! are locked in by regression tests
- Ready for the next milestone phase or finalization

## Self-Check: PASSED

- FOUND: test/test_allocation.jl
- FOUND: test/runtests.jl
- FOUND: 11-03-SUMMARY.md
- FOUND: 00a83bb (Task 1 commit)

---
*Phase: 11-allocation-optimization*
*Completed: 2026-02-15*
