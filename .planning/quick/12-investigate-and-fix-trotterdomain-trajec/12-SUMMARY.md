---
phase: quick-12
plan: 01
subsystem: testing
tags: [regression, trajectory, DM-evolution, platform-portability, Lie-Trotter]

# Dependency graph
requires:
  - phase: 05-02
    provides: "BSON-based regression test infrastructure"
provides:
  - "Platform-portable trajectory regression tests comparing against DM evolution"
  - "Cleaned up reference generator (DM-only, no trajectory BSON)"
affects: [regression-tests, trajectory-validation]

# Tech tracking
tech-stack:
  added: []
  patterns: ["DM-based trajectory validation via exp(delta*L) computed at test time"]

key-files:
  created: []
  modified:
    - test/test_regression.jl
    - test/reference/generate_references.jl

key-decisions:
  - "delta=0.01 (not 0.1) to keep Lie-Trotter splitting bias within atol=0.05"
  - "atol=0.05 accommodates O(delta) splitting bias + O(1/sqrt(N)) statistical noise"
  - "Trajectory BSON files deleted since no test references them"

patterns-established:
  - "DM-based trajectory comparison: compute exp(delta*L) at test time as platform-portable reference"

# Metrics
duration: 8min
completed: 2026-02-14
---

# Quick Task 12: Replace Frozen Trajectory BSON with DM-Based Comparison Summary

**Platform-portable trajectory regression via DM evolution computed at test time with delta=0.01 and atol=0.05**

## Performance

- **Duration:** 8 min
- **Started:** 2026-02-14T17:26:15Z
- **Completed:** 2026-02-14T17:34:13Z
- **Tasks:** 2
- **Files modified:** 3 (1 modified, 2 deleted)

## Accomplishments
- Trajectory regression tests now compare against DM evolution computed fresh at test time, eliminating all platform-dependence
- Removed 2 stale trajectory BSON reference files that caused failures across 3 consecutive quick tasks (10, 11, 12)
- Cleaned up generate_references.jl to only generate the 2 DM reference files still in use
- All 224 tests pass including both trajectory regression tests at atol=0.05

## Task Commits

Each task was committed atomically:

1. **Task 1: Replace frozen trajectory BSON comparison with DM-based comparison** - `895e1b2` (feat)
2. **Task 2: Clean up generate_references.jl and remove stale trajectory BSON files** - `394daf9` (chore)

## Files Created/Modified
- `test/test_regression.jl` - Trajectory regression tests now compute exp(delta*L) as reference instead of loading BSON
- `test/reference/generate_references.jl` - Removed trajectory generator function, constants, and Random import
- `test/reference/energy_traj_reference.bson` - Deleted (no longer used)
- `test/reference/trotter_coherent_traj_reference.bson` - Deleted (no longer used)

## Decisions Made
- **delta=0.01 instead of plan's delta=0.1:** At delta=0.1, the systematic Lie-Trotter splitting bias between trajectory (per-operator Kraus product formula) and DM (full exp(delta*L)) is ~0.10, exceeding the planned atol=0.05. At delta=0.01 the bias is ~0.01, well within tolerance. Verified across 5 different seeds that max element-wise error stays around 0.01.
- **atol=0.05 tolerance:** Accommodates ~0.01 splitting bias + ~0.03 statistical noise (O(1/sqrt(1000))). Gives ~2.5x headroom. Still catches real regressions which produce O(1) errors.
- **Kept `using BSON` in test_regression.jl:** DM regression tests still load frozen BSON reference files.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Changed delta from 0.1 to 0.01 to account for Lie-Trotter splitting bias**
- **Found during:** Task 1 (initial test run failed)
- **Issue:** Plan specified delta=0.1 but the systematic Lie-Trotter splitting error between trajectories (per-operator product formula) and DM (full Liouvillian exp(delta*L)) is O(delta) ~ 0.10, exceeding the atol=0.05 tolerance. Verified error is ~0.108 at 50k trajectories (not statistical).
- **Fix:** Used delta=0.01 where the splitting bias is ~0.01, comfortably within atol=0.05. The plan's core approach (compare trajectories against DM) is correct; only the delta needed adjustment.
- **Files modified:** test/test_regression.jl
- **Verification:** Tested with 5 different seeds at delta=0.01/ntraj=1000, max error consistently ~0.01. Full test suite passes (224/224).
- **Committed in:** 895e1b2 (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (1 bug)
**Impact on plan:** Essential correction. The plan's approach was correct but delta=0.1 was incompatible with atol=0.05 due to Lie-Trotter splitting. No scope creep.

## Issues Encountered
None beyond the delta adjustment documented in deviations.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Trajectory regression tests are now fully platform-portable
- The recurring failure pattern (quick-10, quick-11, quick-12) should be permanently resolved
- DM regression tests remain unchanged at atol=1e-10

## Self-Check: PASSED

- All created/modified files verified to exist on disk
- Both trajectory BSON files confirmed deleted
- Both task commits (895e1b2, 394daf9) verified in git log
- Full test suite passes (224/224)

---
*Phase: quick-12*
*Completed: 2026-02-14*
