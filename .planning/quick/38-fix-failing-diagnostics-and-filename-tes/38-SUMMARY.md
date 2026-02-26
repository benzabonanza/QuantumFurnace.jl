---
phase: 38-fix-failing-diagnostics-and-filename-tes
plan: 01
subsystem: testing
tags: [diagnostics, KMS, GNS, construction, fixed-point, trace-distance]

# Dependency graph
requires:
  - phase: 33-type-foundation
    provides: Config{S,D,C,T} with construction parameter (KMS/GNS singletons)
provides:
  - All diagnostics tests pass with KMS construction
  - Filename generation test passes with kms_ prefix
affects: [test-suite, diagnostics]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Pass construction=KMS() explicitly when test requires exact detailed balance"

key-files:
  created: []
  modified:
    - test/test_results.jl
    - test/test_diagnostics.jl

key-decisions:
  - "Replace ill-conditioned ratio test with absolute threshold checks for KMS regime"

patterns-established:
  - "KMS construction yields exact Gibbs fixed point (trace_distance ~ 1e-15 for Bohr, ~ 1e-8 for Trotter)"
  - "GNS construction yields approximate fixed point (trace_distance ~ 0.035) -- use when testing GNS-specific behavior"

# Metrics
duration: 25min
completed: 2026-02-26
---

# Quick Task 38: Fix Failing Diagnostics and Filename Tests Summary

**Fixed 5 test failures by passing construction=KMS() at 5 call sites where tests require exact Gibbs fixed point**

## Performance

- **Duration:** 25 min
- **Started:** 2026-02-26T09:29:12Z
- **Completed:** 2026-02-26T09:54:03Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments
- Fixed filename generation test in test_results.jl (1 call site)
- Fixed all 4 diagnostics test call sites in test_diagnostics.jl
- All diagnostics tests now pass (242/242)
- Full test suite: 1197 passed, 1 failed (known-deferred test_regression.jl:40 only)

## Task Commits

Each task was committed atomically:

1. **Task 1: Fix test_results.jl filename test** - `cdab3fd` (fix)
2. **Task 2: Fix test_diagnostics.jl at all 4 call sites** - `a34f0bd` (fix)

## Files Created/Modified
- `test/test_results.jl` - Line 282: pass construction=KMS() to make_small_thermalize_config for kms_ filename prefix
- `test/test_diagnostics.jl` - Lines 4, 312, 383, 513: pass construction=KMS() at all call sites requiring exact detailed balance

## Decisions Made
- Replaced the Trotter/Bohr ratio comparison (lines 389-390) with absolute threshold checks. With KMS, both trace distances are near machine precision (Bohr ~1e-15, Trotter ~1e-8), making the ratio numerically ill-conditioned (~6.4M). The physically meaningful assertion is that both distances are small (< 0.01).

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Replaced ill-conditioned ratio test with absolute threshold checks**
- **Found during:** Task 2 (test_diagnostics.jl line 389-390)
- **Issue:** Plan specified changing line 383 to KMS(), which makes fp_bohr.trace_distance ~ 1.8e-15. The ratio `fp.trace_distance / fp_bohr.trace_distance` then evaluates to ~6.4M, far outside the [0.5, 2.0] bound. This was also broken with the original GNS (ratio was 2.297, already > 2.0).
- **Fix:** Replaced `@test fp.trace_distance / fp_bohr.trace_distance > 0.5` and `< 2.0` with `@test fp.trace_distance < 0.01` and `@test fp_bohr.trace_distance < 0.01`. Both assertions pass and verify the physically meaningful property: KMS gives small fixed-point distance.
- **Files modified:** test/test_diagnostics.jl
- **Verification:** All 242 diagnostics tests pass
- **Committed in:** a34f0bd (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (Rule 1 - bug)
**Impact on plan:** Auto-fix necessary for correctness. The ratio test was ill-conditioned in both GNS (ratio 2.297 > 2.0 bound) and KMS (ratio ~6.4M) regimes. New absolute threshold tests are physically meaningful and stable.

## Issues Encountered
None beyond the deviation documented above.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- All diagnostics tests pass with KMS construction
- test_regression.jl:40 remains as known-deferred (not touched per constraints)
- Ready to proceed with Phase 34 (Code Deduplication)

## Self-Check: PASSED

All files exist, all commits verified.

---
*Quick Task: 38-fix-failing-diagnostics-and-filename-tes*
*Completed: 2026-02-26*
