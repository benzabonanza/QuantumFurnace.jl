---
phase: 43-biexp-fitting
plan: 01
subsystem: fitting
tags: [lsqfit, bi-exponential, mixing-time, roots, levenberg-marquardt, spectral-gap]

# Dependency graph
requires:
  - phase: 42-mixing-time-estimation
    provides: FitResult, estimate_mixing_time, MixingTimeEstimate
provides:
  - BiexpFitResult struct for multi-timescale decay fitting
  - fit_biexponential_decay function with mode sorting
  - estimate_mixing_time model=:biexp keyword for bi-exponential extrapolation
  - _extrapolate_mixing_time_biexp using Roots.Bisection numerical solving
affects: [mixing-time, fitting, diagnostics, scripts]

# Tech tracking
tech-stack:
  added: []
  patterns: [bi-exponential decay with mode sorting (fast/slow), residual-seeded initial guess, bisection root-finding for multi-exponential extrapolation]

key-files:
  created: []
  modified:
    - src/fitting.jl
    - src/mixing.jl
    - src/QuantumFurnace.jl
    - test/test_fitting.jl
    - test/test_mixing.jl
    - scripts/mixing_time_extrapolate_verify.jl

key-decisions:
  - "Explicit :biexp model keyword only, no :auto mode with AICc"
  - "Mode sorting after fit (g1 >= g2) with pre-swap SE/CI index tracking"
  - "Bi-exp extrapolation via Roots.Bisection instead of analytical formula"
  - "Synthetic FitResult from slow-mode params for backward-compatible fit_result field"

patterns-established:
  - "Multi-exponential fitting with post-fit mode sorting for consistent slow/fast separation"
  - "Residual analysis of single-exp fit to seed bi-exponential initial guess"

# Metrics
duration: 7min
completed: 2026-03-04
---

# Phase 43: Bi-Exponential Fitting Summary

**Bi-exponential decay fitting with Roots.Bisection extrapolation, reducing mixing time error from ~0.13% to <0.001% on synthetic data**

## Performance

- **Duration:** 7 min
- **Started:** 2026-03-04T09:34:15Z
- **Completed:** 2026-03-04T09:41:18Z
- **Tasks:** 6
- **Files modified:** 6

## Accomplishments
- BiexpFitResult struct and fit_biexponential_decay function with automatic initial guess via single-exp residual analysis
- estimate_mixing_time extended with model=:biexp keyword; uses Roots.Bisection for numerical extrapolation of A1*exp(-g1*t) + A2*exp(-g2*t) + C = epsilon
- Offset accuracy: bi-exp error 1.2e-18 vs single-exp error 1.8e-9 on synthetic data (BIEXP-02)
- Extrapolation accuracy: bi-exp error 6.8e-12 (<<5% threshold) vs single-exp error 0.13% (BIEXP-MIX-01)
- Full backward compatibility: default model=:single produces identical results, all 1273 tests pass (+27 new)

## Task Commits

Each task was committed atomically:

1. **Task 1: src/fitting.jl -- bi-exponential fitting** - `094184f` (feat)
2. **Task 2: src/mixing.jl -- model=:biexp support** - `13324e2` (feat)
3. **Task 3: src/QuantumFurnace.jl -- exports** - `94e2c50` (feat)
4. **Task 4: test/test_fitting.jl -- BIEXP-01/02/03** - `c82c63b` (test)
5. **Task 5: test/test_mixing.jl -- BIEXP-MIX-01/02** - `3b2d690` (test)
6. **Task 6: scripts/mixing_time_extrapolate_verify.jl -- demo update** - `6ccccf1` (feat)

## Files Created/Modified
- `src/fitting.jl` - Added BiexpFitResult struct, _biexp_decay_model, _biexp_initial_guess, fit_biexponential_decay (+220 lines)
- `src/mixing.jl` - Added model_used/biexp_fit_result fields to MixingTimeEstimate, _extrapolate_mixing_time_biexp, _biexp_to_single_fit_result, model keyword
- `src/QuantumFurnace.jl` - Added export for fit_biexponential_decay, BiexpFitResult
- `test/test_fitting.jl` - BIEXP-01 (data recovery), BIEXP-02 (offset accuracy), BIEXP-03 (skip_initial), edge case
- `test/test_mixing.jl` - BIEXP-MIX-01 (<5% extrapolation accuracy), BIEXP-MIX-02 (backward compat), edge case, MIX-05 field update
- `scripts/mixing_time_extrapolate_verify.jl` - Side-by-side single vs biexp comparison with error percentages

## Decisions Made
- Explicit `:biexp` model keyword only (no `:auto` mode with AICc) per user decision
- Post-fit mode sorting (g1 >= g2) with pre-swap index tracking for correct SE/CI extraction
- Bi-exponential extrapolation via Roots.Bisection rather than analytical formula (multi-exponential has no closed-form solution for t)
- Synthetic FitResult constructed from slow-mode params for backward-compatible fit_result field

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Bi-exponential fitting fully operational and tested
- All 1273 tests pass (1246 existing + 27 new)
- Acceptance criterion met: <5% extrapolation error on synthetic bi-exp data (actual: <0.001%)
- Script ready for real-data validation: `julia --project scripts/mixing_time_extrapolate_verify.jl`

---
*Phase: 43-biexp-fitting*
*Completed: 2026-03-04*
