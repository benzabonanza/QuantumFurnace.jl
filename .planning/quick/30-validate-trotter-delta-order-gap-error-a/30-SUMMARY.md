---
phase: quick-30
plan: 01
subsystem: validation
tags: [spectral-gap, trotter, richardson-extrapolation, delta-scaling, disordered-heisenberg]

# Dependency graph
requires:
  - phase: quick-28
    provides: "Disordered Heisenberg gap estimation baseline (delta=0.01)"
  - phase: quick-29
    provides: "Observable gap-mode coupling analysis"
provides:
  - "Empirical evidence that gap estimation error is NOT O(delta) in Trotter step size"
  - "Richardson extrapolation tested and shown ineffective (1.0x improvement)"
  - "Per-observable delta-scaling data for 8 observables at 3 delta values"
affects: [spectral-gap-estimation, trotter-discretization, gap-accuracy]

# Tech tracking
tech-stack:
  added: []
  patterns: [delta-override-per-run, richardson-extrapolation-formula]

key-files:
  created:
    - experiments/validate_gap_delta_scaling.jl
  modified: []

key-decisions:
  - "Gap estimation error is NOT O(delta) -- error/delta varies 96x across delta=0.1,0.01,0.001"
  - "Richardson extrapolation provides only 1.0x improvement -- not useful"
  - "Dominant error source is NOT Trotter discretization but systematic observable bias"
  - "YY_avg shows smallest and most stable error across delta values (3.6-6.1%)"

patterns-established:
  - "delta keyword override in estimate_spectral_gap for multi-delta comparison"

# Metrics
duration: 12min
completed: 2026-02-18
---

# Quick Task 30: Delta-Scaling Validation & Richardson Extrapolation Summary

**Gap estimation error does NOT scale as O(delta) with Trotter step size; Richardson extrapolation provides no improvement (1.0x), indicating dominant error is systematic observable bias not Trotter discretization**

## Performance

- **Duration:** ~12 min
- **Started:** 2026-02-18T19:53:39Z
- **Completed:** 2026-02-18T20:05:39Z
- **Tasks:** 2
- **Files modified:** 1

## Accomplishments

- Ran n=4 disordered Heisenberg gap estimation at delta=0.1, 0.01, 0.001 (20k trajectories each)
- Demonstrated error/delta ratio varies by 96x across delta values (NOT O(delta))
- Richardson extrapolation tested for two pairs -- only 1.0x improvement (ineffective)
- Per-observable analysis reveals YY_avg has most stable and smallest error (2-6%)

## Key Numerical Results

### Best-Observable (Mz_stagg) Gap Estimates

| delta  | est_gap    | exact_gap  | error      | rel_err% | error/delta |
|--------|------------|------------|------------|----------|-------------|
| 0.1000 | 0.10597865 | 0.17298091 | -0.06700226| 38.73%   | -0.670      |
| 0.0100 | 0.08870246 | 0.17298091 | -0.08427845| 48.72%   | -8.428      |
| 0.0010 | 0.10869142 | 0.17298091 | -0.06428949| 37.17%   | -64.289     |

Error/delta varies from 0.67 to 64.3 (96x spread). O(delta) scaling requires constant ratio.

### Richardson Extrapolation

| Pair           | gap_rich   | error_rich | rel_err% | improvement |
|----------------|------------|------------|----------|-------------|
| (0.01, 0.1)   | 0.08678288 | 0.08619803 | 49.83%   | 1.0x        |
| (0.001, 0.01) | 0.11091242 | 0.06206849 | 35.88%   | 1.0x        |

Richardson extrapolation is ineffective because the error is not dominated by O(delta) term.

### Per-Observable Analysis (Highlights)

| Observable | error/delta (0.1) | error/delta (0.01) | error/delta (0.001) | Stable? |
|------------|-------------------|--------------------|---------------------|---------|
| H          | +0.449            | +3.697             | +40.125             | No      |
| Mz         | +0.951            | +9.213             | +90.259             | No      |
| XX_avg     | +0.001            | -1.288             | -5.305              | No      |
| YY_avg     | -0.056            | -0.615             | -3.853              | No      |
| ZZ_avg     | +1.013            | +8.866             | +92.910             | No      |
| Mz_stagg   | -0.670            | -8.428             | -64.289             | No      |
| Z1         | +1.296            | +11.623            | +118.755            | No      |

No observable shows O(delta) scaling. Error/delta grows roughly as 1/delta, suggesting the dominant error is delta-independent (constant bias).

### Per-Observable Richardson (Best Results)

| Observable | Rich (0.01,0.1) rel% | Rich (0.001,0.01) rel% |
|------------|----------------------|------------------------|
| YY_avg     | 3.59%                | 2.08%                  |
| XX_avg     | 8.28%                | 2.58%                  |
| H          | 20.86%               | 23.40%                 |

YY_avg and XX_avg have the smallest absolute errors and benefit marginally from Richardson.

## Task Commits

Each task was committed atomically:

1. **Task 1: Create delta-scaling validation script with Richardson extrapolation** - `502afe4` (feat)
2. **Task 2: Write summary documenting delta-scaling and Richardson results** - (this commit, docs)

## Files Created/Modified

- `experiments/validate_gap_delta_scaling.jl` - Delta-scaling validation script with Richardson extrapolation (266 lines)

## Decisions Made

- **Error is NOT O(delta):** error/delta varies by 96x (0.67 to 64.3) for Mz_stagg across delta=0.1, 0.01, 0.001. This rules out Trotter discretization as the dominant error source.
- **Richardson extrapolation is ineffective:** 1.0x improvement for both pairs. Since error is not O(delta), canceling the leading delta term has no effect.
- **Dominant error is systematic observable bias:** The error/delta ratio growing as ~1/delta means the error is approximately constant (delta-independent). This is consistent with the exponential fitting picking up a mode other than the true gap mode.
- **YY_avg is most accurate:** YY_avg achieves 2-6% relative error across all delta values, far better than Mz_stagg (37-49%). The "best observable" selection criterion (smallest gap) systematically picks Mz_stagg, which underestimates.
- **Non-monotonic error in delta:** Error at delta=0.01 (49%) is LARGER than at delta=0.1 (39%) and delta=0.001 (37%), ruling out simple O(delta) or O(delta^2) scaling.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None - script ran successfully for all three delta values.

## Key Insights

1. **Trotter step size is NOT the bottleneck:** Reducing delta from 0.1 to 0.001 (100x smaller) does not meaningfully reduce gap estimation error. The ~37-49% error at all delta values indicates a fundamental limitation in the observable-based estimation approach.

2. **Best-observable selection criterion may need revisiting:** The "smallest gap" criterion consistently selects Mz_stagg, which has large systematic bias. YY_avg at 2-6% error would be far more accurate but is not selected because its gap estimate (0.167-0.169) is not the smallest.

3. **Richardson extrapolation premise falsified:** Richardson extrapolation requires error = C*delta^p for some known p. The data shows error is essentially independent of delta, so no extrapolation order can help.

## Next Steps

- Consider revising best-observable selection to use accuracy-aware criteria (e.g., cross-validation between observables)
- Investigate why Mz_stagg systematically underestimates -- may be fitting a faster-decaying mode
- YY_avg appears to be the most reliable observable for disordered n=4 Heisenberg

## User Setup Required

None - no external service configuration required.

---
*Quick Task: 30-validate-trotter-delta-order-gap-error-a*
*Completed: 2026-02-18*
