---
phase: quick-5
plan: 01
subsystem: testing
tags: [trotter, error-thresholds, dm-scaling, regression-safety]

# Dependency graph
requires:
  - phase: quick-3
    provides: "Fixed OFT consistency test basis transformation (DMTST-06)"
  - phase: quick-4
    provides: "NUFFT OFT consistency test (DMTST-06b)"
provides:
  - "Tight TrotterDomain error thresholds catching regressions"
affects: [phase-04, trajectory-validation]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Threshold = 2x measured error for physics-limited tests (DMTST-05)"
    - "Threshold = 1e-5 for numerically exact Trotter OFT tests (DMTST-06, 06b)"

key-files:
  created: []
  modified:
    - test/test_dm_scaling.jl

key-decisions:
  - "DMTST-05 threshold 0.02 (not 1e-5) because B-term Trotter error ~0.011 is a real physical approximation"
  - "DMTST-06 and DMTST-06b thresholds 1e-5 since measured OFT Trotter errors are ~1.5e-8"

patterns-established:
  - "Trotter threshold tiers: physics-limited (2x margin) vs numerically-exact (orders of magnitude margin)"

# Metrics
duration: 1min
completed: 2026-02-14
---

# Quick Task 5: Tighten TrotterDomain Error Thresholds Summary

**Replaced three loose 0.1 Trotter error thresholds with tight values (0.02 for B-term, 1e-5 for OFT) based on measured errors**

## Performance

- **Duration:** 1 min
- **Started:** 2026-02-14T11:14:04Z
- **Completed:** 2026-02-14T11:15:11Z
- **Tasks:** 1
- **Files modified:** 1

## Accomplishments
- Tightened DMTST-05 B-term threshold from 0.1 to 0.02 (measured error ~0.011, 2x safety margin)
- Tightened DMTST-06 OFT threshold from 0.1 to 1e-5 (measured error ~1.5e-8, 680x margin)
- Tightened DMTST-06b NUFFT OFT threshold from 0.1 to 1e-5 (measured error ~1.5e-8, 680x margin)
- All 18 tests pass with the tightened thresholds

## Task Commits

Each task was committed atomically:

1. **Task 1: Tighten TrotterDomain error thresholds** - `eb2569b` (feat)

## Files Created/Modified
- `test/test_dm_scaling.jl` - Updated three Trotter error thresholds and their comments

## Decisions Made
- DMTST-05 threshold set to 0.02 (not 1e-5) because the B-term Trotter approximation error (~0.011) is a real physical effect, not a numerical artifact
- DMTST-06 and DMTST-06b thresholds set to 1e-5 since the measured OFT Trotter errors (~1.5e-8) are orders of magnitude smaller, indicating numerical precision rather than physics-limited accuracy

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- All DM scaling tests have tight thresholds, providing strong regression detection
- Ready for Phase 4 (trajectory validation) which depends on these DM reference tests

---
*Quick task: 5*
*Completed: 2026-02-14*
