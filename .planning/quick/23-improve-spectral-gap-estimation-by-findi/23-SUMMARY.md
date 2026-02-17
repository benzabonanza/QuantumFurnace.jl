---
phase: quick-23
plan: 01
subsystem: gap-estimation
tags: [spectral-gap, observable-selection, exponential-fitting, cross-validation]

# Dependency graph
requires:
  - phase: v1.3
    provides: "estimate_spectral_gap, _select_best_observable, cross_validate_gap"
  - phase: quick-22
    provides: "Corrected trajectory CPTP channel (bare delta, R_a scaled by n_jumps)"
provides:
  - "Smallest-gap-among-good-fits selection criterion for _select_best_observable"
  - "Tightened validation pass criterion [0.8, 1.5] replacing [1.0, 3.0]"
affects: [gap-estimation, cross-validation, experiments]

# Tech tracking
tech-stack:
  added: []
  patterns: ["smallest-gap selection: filter by quality then pick minimum gap"]

key-files:
  created: []
  modified:
    - src/gap_estimation.jl
    - test/test_gap_estimation.jl
    - experiments/validate_gap_estimation.jl

key-decisions:
  - "Smallest gap among good fits (R-squared > 0.8) replaces highest R-squared selection"
  - "Tightened validation pass criterion from [1.0, 3.0] to [0.8, 1.5]"
  - "n=4 residual factor improved from ~1.6x to ~1.17x; n=6 at ~1.46x (both pass)"

patterns-established:
  - "Observable selection: filter by convergence+positivity+quality, then minimize gap"

# Metrics
duration: 6min
completed: 2026-02-17
---

# Quick Task 23: Improve Spectral Gap Estimation Summary

**Smallest-gap selection criterion for _select_best_observable, reducing n=4 residual factor from ~1.6x to ~1.17x**

## Performance

- **Duration:** 6 min
- **Started:** 2026-02-17T15:33:21Z
- **Completed:** 2026-02-17T15:40:10Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments
- Replaced highest-R-squared observable selection with smallest-gap-among-good-fits criterion
- n=4 residual factor improved from ~1.6x to ~1.17x (closest to true spectral gap)
- n=6 residual factor at ~1.46x, both systems pass tightened [0.8, 1.5] criterion
- Added test proving smallest-gap selection beats highest-R-squared when two valid fits exist

## Task Commits

Each task was committed atomically:

1. **Task 1: Fix _select_best_observable and add test** - `c5f6c68` (fix)
2. **Task 2: Update validation script and delete diagnostic** - `fa95a0c` (chore)

## Files Created/Modified
- `src/gap_estimation.jl` - Replaced _select_best_observable: smallest gap among fits with R-squared > 0.8
- `test/test_gap_estimation.jl` - Added smallest-gap selection test; updated existing test comment
- `experiments/validate_gap_estimation.jl` - Tightened pass criterion [0.8, 1.5], updated comments
- `experiments/diagnose_observable_overlap.jl` - Deleted (untracked temporary diagnostic)

## Decisions Made
- Selection threshold R-squared > 0.8 (not 0.9): allows more observables into the candidate pool while still filtering noise fits
- Three-tier fallback: (1) smallest gap with R-squared > 0.8, (2) smallest gap with any R-squared, (3) highest R-squared overall
- Pass criterion [0.8, 1.5] chosen based on observed residual factors: n=4 at 1.17, n=6 at 1.46

## Deviations from Plan

None - plan executed exactly as written.

Note: The diagnostic file `experiments/diagnose_observable_overlap.jl` was never tracked in git, so its deletion did not produce a git change. This is expected behavior, not a deviation.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Gap estimation now uses physically correct observable selection
- Cross-validation passes with tightened criteria for both n=4 and n=6
- No blockers for future work

## Self-Check: PASSED

- All 3 modified files exist on disk
- Diagnostic file confirmed deleted
- Commit c5f6c68 (Task 1) found in git log
- Commit fa95a0c (Task 2) found in git log
- All 63 unit tests pass
- Validation script passes OVERALL: PASS for both n=4 and n=6

---
*Phase: quick-23*
*Completed: 2026-02-17*
