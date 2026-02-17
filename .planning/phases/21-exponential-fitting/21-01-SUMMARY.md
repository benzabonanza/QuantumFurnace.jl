---
phase: 21-exponential-fitting
plan: 01
subsystem: fitting
tags: [lsqfit, exponential-decay, curve-fitting, spectral-gap, levenberg-marquardt]

# Dependency graph
requires: []
provides:
  - "FitResult struct with gap, CI, R-squared, standard error, convergence status"
  - "fit_exponential_decay function with auto log-linear initial guess"
  - "skip_initial window selection for transient exclusion"
  - "LsqFit.jl dependency with compat bounds"
affects: [23-gap-estimation-api]

# Tech tracking
tech-stack:
  added: [LsqFit.jl v0.15]
  patterns: [result-struct-with-quality-metrics, log-linear-initial-guess, parameter-bounds-enforcement]

key-files:
  created:
    - src/fitting.jl
    - test/test_fitting.jl
  modified:
    - Project.toml
    - src/QuantumFurnace.jl
    - test/runtests.jl

key-decisions:
  - "Use _IDX_A=1, _IDX_GAP=2, _IDX_C=3 parameter ordering with named constants"
  - "R-squared NOT clamped to [0,1] -- negative values are valid diagnostics for bad fits"
  - "No weight vector passed to curve_fit (avoids corrupted covariance per LsqFit docs)"

patterns-established:
  - "FitResult struct pattern: wrap LsqFit output in clean API boundary (downstream does not need using LsqFit)"
  - "Log-linear initial guess with fallback: estimate plateau, subtract, regress on log-shifted data"

# Metrics
duration: 5min
completed: 2026-02-17
---

# Phase 21 Plan 01: Exponential Fitting Summary

**Single-exponential decay fitting (A*exp(-gap*t)+C) via LsqFit.jl with auto log-linear initial guess, bounded gap >= 0, skip_initial window selection, and comprehensive quality metrics**

## Performance

- **Duration:** 5 min
- **Started:** 2026-02-17T08:39:04Z
- **Completed:** 2026-02-17T08:44:28Z
- **Tasks:** 2
- **Files modified:** 5

## Accomplishments
- FitResult struct and fit_exponential_decay function with full quality metrics (gap, CI, R-squared, SE, convergence)
- Automatic log-linear initial guess generation with fallback for degenerate data
- skip_initial parameter measurably improves gap estimate when early transient is present
- 9 testsets covering all FIT-01 through FIT-05 requirements, all 610 tests pass (no regressions)

## Task Commits

Each task was committed atomically:

1. **Task 1: Add LsqFit.jl dependency and implement FitResult + fit_exponential_decay** - `0eafd22` (feat)
2. **Task 2: Add comprehensive synthetic data tests for all FIT-* requirements** - `e50f29c` (test)

## Files Created/Modified
- `src/fitting.jl` - FitResult struct, fit_exponential_decay function, _log_linear_initial_guess, _compute_r_squared helpers
- `test/test_fitting.jl` - 9 testsets for all FIT-* requirements with synthetic data
- `Project.toml` - LsqFit.jl added to [deps] and [compat] with "0.15" bound
- `src/QuantumFurnace.jl` - Added using LsqFit, include("fitting.jl"), export fit_exponential_decay/FitResult
- `test/runtests.jl` - Added include("test_fitting.jl")

## Decisions Made
- Used named index constants (_IDX_A, _IDX_GAP, _IDX_C) for parameter ordering clarity and safety
- R-squared computed as 1 - RSS/TSS without clamping; negative values serve as valid diagnostic for bad fits
- No weight vector passed to curve_fit to avoid corrupted covariance (per LsqFit.jl documentation warning)
- LsqFit compat bound set to "0.15" (permissive) rather than "0.15.1" (pinned)

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- fit_exponential_decay and FitResult are exported and ready for Phase 23 (Gap Estimation API)
- All 610 tests pass including new fitting tests; no regressions in existing tests
- LsqFit.jl v0.15 installed and precompiled

## Self-Check: PASSED

All files verified present: src/fitting.jl, test/test_fitting.jl, 21-01-SUMMARY.md
All commits verified: 0eafd22, e50f29c

---
*Phase: 21-exponential-fitting*
*Completed: 2026-02-17*
