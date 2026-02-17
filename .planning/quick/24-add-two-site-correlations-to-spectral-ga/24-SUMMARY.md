---
phase: quick-24
plan: 01
subsystem: convergence
tags: [spectral-gap, observables, two-site-correlations, gap-estimation]

# Dependency graph
requires:
  - phase: 20-gap-observables
    provides: "build_gap_estimation_observables with H + Mz bundle"
  - phase: quick-23
    provides: "smallest-gap observable selection criterion"
provides:
  - "build_gap_estimation_observables returning 5 observables: H, Mz, XX_avg, YY_avg, ZZ_avg"
  - "Per-bond averaged nearest-neighbor two-site correlations for gap estimation"
  - "n=6 validation script for measuring correlation observable improvement"
affects: [gap-estimation, spectral-gap, validation]

# Tech tracking
tech-stack:
  added: []
  patterns: [per-bond-averaged-correlations]

key-files:
  created:
    - "experiments/run_gap_validation.jl"
  modified:
    - "src/convergence.jl"
    - "src/fitting.jl"
    - "test/test_gap_estimation.jl"
    - "test/test_convergence.jl"

key-decisions:
  - "XX_avg, YY_avg, ZZ_avg use per-bond averaging (divide by n) matching Mz per-site normalization"
  - "Correlations use pad_term with periodic=true for all nearest-neighbor bonds"
  - "SingularException in LsqFit handled gracefully with Inf/(-Inf,Inf) fallback"

patterns-established:
  - "Per-bond averaged correlation observable pattern: sum over sites, divide by n, transform to eigenbasis"

# Metrics
duration: 11min
completed: 2026-02-17
---

# Quick Task 24: Add Two-Site Correlations to Spectral Gap Estimation Summary

**Added XX_avg, YY_avg, ZZ_avg per-bond averaged nearest-neighbor correlations to gap estimation observable bundle (now 5 observables total), with n=6 validation script**

## Performance

- **Duration:** 11 min
- **Started:** 2026-02-17T15:51:45Z
- **Completed:** 2026-02-17T16:02:38Z
- **Tasks:** 2
- **Files modified:** 5

## Accomplishments
- `build_gap_estimation_observables` now returns 5 observables: H, Mz, XX_avg, YY_avg, ZZ_avg
- Each two-site correlation is per-bond averaged over all nearest-neighbor bonds (periodic), correctly transformed to the selected eigenbasis (Hamiltonian or Trotter)
- All 725 tests pass, including updated convergence and gap estimation test suites
- Standalone validation script created for n=6 gap estimation with 10k trajectories

## Task Commits

Each task was committed atomically:

1. **Task 1: Add XX_avg, YY_avg, ZZ_avg to build_gap_estimation_observables** - `509d540` (feat)
2. **Task 2: Update tests for 5 observables and create validation script** - `595aba3` (feat)

## Files Created/Modified
- `src/convergence.jl` - Added XX_avg, YY_avg, ZZ_avg construction loop to build_gap_estimation_observables
- `src/fitting.jl` - Added SingularException handling in stderror/confint extraction
- `test/test_gap_estimation.jl` - Updated assertions from 2 to 5 observables, added length check
- `test/test_convergence.jl` - Updated eigenbasis and Trotter basis tests for 5-observable bundle with correlation verification
- `experiments/run_gap_validation.jl` - New standalone n=6 validation script with per-observable fit summary

## Decisions Made
- XX_avg, YY_avg, ZZ_avg use per-bond averaging (divide by num_qubits) matching existing Mz per-site normalization pattern. This keeps amplitudes system-size-independent while preserving decay rates.
- Used mutable push!/append! pattern instead of vcat for observable construction to support the correlation loop cleanly.
- SingularException in LsqFit.stderror is handled gracefully with Inf standard error and (-Inf, Inf) confidence interval, rather than crashing. This is appropriate since the fit itself may still converge even if the covariance matrix is singular.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Updated test_convergence.jl for 5 observables**
- **Found during:** Task 2 (test execution)
- **Issue:** test_convergence.jl had hardcoded expectations of 2 observables and names == ["H", "Mz"] in both eigenbasis and Trotter basis testsets
- **Fix:** Updated both testsets to expect 5 observables with ["H", "Mz", "XX_avg", "YY_avg", "ZZ_avg"], added per-bond correlation verification tests
- **Files modified:** test/test_convergence.jl
- **Verification:** All 725 tests pass
- **Committed in:** 595aba3 (Task 2 commit)

**2. [Rule 1 - Bug] Fixed SingularException in fit_exponential_decay**
- **Found during:** Task 2 (test execution)
- **Issue:** With 5 observables and skip_initial=0.2, one correlation observable produced near-flat data causing LsqFit.stderror to throw SingularException(1) from singular Jacobian
- **Fix:** Added try-catch around stderror/confint calls, returning Inf/(-Inf,Inf) on SingularException
- **Files modified:** src/fitting.jl
- **Verification:** skip_initial test now passes; all 725 tests pass
- **Committed in:** 595aba3 (Task 2 commit)

---

**Total deviations:** 2 auto-fixed (1 blocking, 1 bug)
**Impact on plan:** Both auto-fixes necessary for test correctness. No scope creep.

## Issues Encountered
None beyond the auto-fixed deviations above.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- The validation script `experiments/run_gap_validation.jl` is ready for the user to run and evaluate whether the additional correlation observables reduce the ~1.46x residual factor at n=6
- All existing functionality preserved; the 5-observable bundle is backward-compatible since estimate_spectral_gap auto-selects the best observable

---
*Quick task: 24-add-two-site-correlations*
*Completed: 2026-02-17*
