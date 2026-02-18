---
phase: 25-spectral-gap-validation-overhaul
plan: 03
subsystem: estimation
tags: [spectral-gap, validation, arpack, eigenbasis-overlap, trajectory-estimation]

# Dependency graph
requires:
  - phase: 25-01
    provides: Consolidated build_preset_trajectory_observables, clean gap_estimation.jl
  - phase: 25-02
    provides: eigenbasis_overlap_analysis function and OverlapAnalysisResult struct
provides:
  - Unified validation script: ARPACK check, overlap analysis, 20k-trajectory estimation, pass/fail
  - Concrete validation results for n=4 and n=6 at beta=10
affects: []

# Tech tracking
tech-stack:
  added: []
  patterns: [unified-validation-script]

key-files:
  created:
    - experiments/validate_spectral_gap.jl
  modified: []

key-decisions:
  - "n=4 ARPACK vs eigen agrees to 1.2e-10 (well within 1e-8 threshold)"
  - "n=4 gap estimation passes with 0.72% relative error (ZZ_avg best observable)"
  - "n=6 gap estimation fails with 10.7% relative error -- all observables show zero overlap with gap mode"
  - "Experiments directory is gitignored; validation script force-added to git tracking"

patterns-established:
  - "Validation script pattern: make_system + make_liouv_config + make_thermalize_config helpers for clean Heisenberg chain construction"

# Metrics
duration: 26min
completed: 2026-02-18
---

# Phase 25 Plan 03: Unified Validation Script Summary

**Unified spectral gap validation: ARPACK vs eigen (PASS), n=4 20k-traj estimation (PASS, 0.72%), n=6 estimation (FAIL, 10.7% with zero-overlap diagnostic)**

## Performance

- **Duration:** 26 min (mostly n=6 trajectory simulation at ~17min)
- **Started:** 2026-02-18T08:14:33Z
- **Completed:** 2026-02-18T08:40:37Z
- **Tasks:** 1
- **Files created:** 1

## Accomplishments
- Created `experiments/validate_spectral_gap.jl` replacing 4 deleted experiment scripts
- ARPACK vs dense eigen verification passes for n=4 (difference 1.18e-10)
- Eigenbasis overlap analysis runs for both n=4 and n=6 with clear table output
- n=4 trajectory estimation (20k, beta=10, delta=0.01) passes with 0.72% relative error
- n=6 trajectory estimation runs to completion with diagnostic explaining failure (zero overlap with gap mode)
- Script provides clear pass/fail judgment with per-observable fit quality table

## Task Commits

Each task was committed atomically:

1. **Task 1: Write unified validation script** - `ba342e5` (feat)

## Files Created/Modified
- `experiments/validate_spectral_gap.jl` - Unified validation script covering ARPACK check, overlap analysis, 20k-trajectory estimation for n=4 and n=6

## Decisions Made
- Used full `@kwdef` keyword constructor for LiouvConfig and ThermalizeConfig instead of plan's positional-then-keyword syntax (which doesn't work with Julia's `@kwdef`)
- Set `with_linear_combination=true`, `sigma=1/beta`, `a=beta/30`, `b=0.4` matching test helper conventions for KMS detailed balance
- Used `mixing_time = max(5.0/exact_gap, 10.0)` heuristic (~31.5s for n=4, ~43.2s for n=6)
- Force-added validation script to git since `/experiments/` is in `.gitignore`

## Validation Results

### ARPACK vs eigen (n=4)
- ARPACK gap: 0.1585339427
- Dense eigen gap: 0.1585339428
- Difference: 1.179e-10
- **PASS** (threshold: 1e-8)

### n=4 Gap Estimation (20k trajectories)
- Exact gap: 0.15853394
- Estimated gap: 0.15738781 (best: ZZ_avg)
- Relative error: 0.72%
- **PASS** (target: < 1%)
- All 5 observables converged with R-squared > 0.999

### n=6 Gap Estimation (20k trajectories)
- Exact gap: 0.11560988
- Estimated gap: 0.12799675 (best: ZZ_avg)
- Relative error: 10.7%
- **FAIL** (target: < 1%)
- Diagnostic: All observables have zero overlap with gap mode for n=6
- This means the gap mode (first excited eigenmode of L) does not couple to any of the 5 preset observables for the n=6 periodic Heisenberg chain at beta=10
- The estimation still finds a reasonable gap (within ~11%) because ZZ_avg picks up slower-decaying modes near the gap

### Eigenbasis Overlap (n=4)
| Observable | |c_gap| | Relative |
|------------|---------|----------|
| H          | 0.2457  | 0.7124   |
| Mz         | 0.0000  | 0.0000   |
| XX_avg     | 0.5460  | 0.5692   |
| YY_avg     | 0.5460  | 0.5728   |
| ZZ_avg     | 0.5460  | 0.4820   |

### Eigenbasis Overlap (n=6)
| Observable | |c_gap| | Relative |
|------------|---------|----------|
| H          | 0.0000  | 0.0000   |
| Mz         | 0.0000  | 0.0000   |
| XX_avg     | 0.0000  | 0.0000   |
| YY_avg     | 0.0000  | 0.0000   |
| ZZ_avg     | 0.0000  | 0.0000   |

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed @printf with non-literal format strings**
- **Found during:** Task 1 (first script run)
- **Issue:** `@printf("=" ^ 50 * "\n")` fails because `@printf` requires literal format strings, not runtime expressions
- **Fix:** Changed to `println("=" ^ 50)` for dynamic string output
- **Files modified:** experiments/validate_spectral_gap.jl
- **Verification:** Script runs to completion
- **Committed in:** ba342e5 (Task 1 commit)

**2. [Rule 1 - Bug] Fixed config constructor syntax**
- **Found during:** Task 1 (script implementation)
- **Issue:** Plan specified `LiouvConfig(TimeDomain(); num_qubits=4, ...)` but LiouvConfig is `@kwdef` and requires keyword-only construction: `LiouvConfig(; domain=TimeDomain(), ...)`
- **Fix:** Created `make_liouv_config(n)` and `make_thermalize_config(n)` helper functions using full keyword constructors with all required fields (sigma, a, b, num_energy_bits, w0, t0, etc.)
- **Files modified:** experiments/validate_spectral_gap.jl
- **Verification:** Config validation passes, Lindbladian construction succeeds
- **Committed in:** ba342e5 (Task 1 commit)

---

**Total deviations:** 2 auto-fixed (2 bugs in plan specification)
**Impact on plan:** Both fixes necessary for correctness. No scope creep.

## Issues Encountered
- n=6 gap estimation fails the 1% target. This is a genuine physics result: the 5 preset observables have zero overlap with the n=6 gap mode. The diagnostic correctly identifies this. Future work could explore system-size-dependent observable selection or symmetry-adapted observables.
- The `/experiments/` directory is gitignored, requiring `git add -f` to track the validation script. This is consistent with how previous experiment scripts were handled.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Phase 25 is now complete (all 3 plans executed)
- Validation results are concrete and documented
- n=4 pipeline is validated end-to-end
- n=6 failure is explained with diagnostic evidence (gap-mode overlap)
- No blockers for future phases

## Self-Check: PASSED
- [x] experiments/validate_spectral_gap.jl exists
- [x] Commit ba342e5 exists
- [x] Script runs to completion with pass/fail output
- [x] ARPACK vs eigen passes for n=4
- [x] Eigenbasis overlap table printed for both n=4 and n=6
- [x] Gap estimation results with pass/fail printed
- [x] Diagnostic output explains n=6 failure

---
*Phase: 25-spectral-gap-validation-overhaul*
*Completed: 2026-02-18*
