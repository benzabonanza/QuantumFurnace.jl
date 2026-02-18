---
phase: 27-add-symmetry-breaking-observable
plan: 01
subsystem: physics
tags: [observable, spectral-gap, SU2, symmetry-breaking, heisenberg-chain]

# Dependency graph
requires:
  - phase: quick-26
    provides: "7-observable bundle with Mz_stagg/Z1; SU(2) gap-mode protection diagnosis"
provides:
  - "8-observable bundle with XZ_stagg (staggered nearest-neighbor XZ correlation)"
  - "Empirical confirmation that n=6 gap mode is protected beyond simple SU(2) breaking"
affects: [gap-estimation, observable-design, symmetry-analysis]

# Tech tracking
tech-stack:
  added: []
  patterns: [staggered-asymmetric-pauli-correlation]

key-files:
  created: []
  modified:
    - src/convergence.jl
    - test/test_gap_estimation.jl
    - test/test_convergence.jl
    - experiments/validate_spectral_gap.jl

key-decisions:
  - "XZ_stagg has |c_gap|=0 for n=6 despite breaking SU(2) -- gap mode protection is stronger than spin-rotation symmetry alone"
  - "Keep TARGET_REL_ERROR_N6 at 12% -- XZ_stagg did not improve n=6 estimation"
  - "n=4 unaffected (0.72% error, ZZ_avg best) -- no regression from adding 8th observable"

patterns-established:
  - "Asymmetric Pauli pair observable: pad_term([X, Z], ...) with staggered sign for mixed-operator correlations"

# Metrics
duration: 24min
completed: 2026-02-18
---

# Quick Task 27: Add Symmetry-Breaking Observable Summary

**Added XZ_stagg staggered XZ correlation as 8th observable; confirmed n=6 gap mode protected beyond SU(2) (|c_gap|=0 for all 8 observables)**

## Performance

- **Duration:** 24 min
- **Started:** 2026-02-18T11:01:01Z
- **Completed:** 2026-02-18T11:25:30Z
- **Tasks:** 3
- **Files modified:** 4

## Accomplishments
- Added XZ_stagg = sum((-1)^i X_i Z_{i+1}) / n as 8th observable in build_preset_trajectory_observables
- Updated all test assertions (7 -> 8 observables) across test_gap_estimation.jl and test_convergence.jl; all 696 tests pass
- Validated with 20k trajectories: n=4 at 0.72% (PASS), n=6 at 10.71% (PASS under 12% threshold)
- Confirmed XZ_stagg has zero gap-mode overlap for n=6 -- SU(2) breaking via asymmetric Pauli pair is insufficient

## Task Commits

Each task was committed atomically:

1. **Task 1: Add XZ_stagg observable** - `e6d0694` (feat)
2. **Task 2: Update test assertions 7->8** - `66bab95` (test)
3. **Task 3: Run validation, confirm results** - `a0f4c62` (feat)

## Files Created/Modified
- `src/convergence.jl` - Added XZ_stagg construction block and updated docstring (7->8 observables)
- `test/test_gap_estimation.jl` - Updated count assertions, name lists, overlap dimensions (7->8)
- `test/test_convergence.jl` - Updated count assertions, name lists, added XZ_stagg construction verification test
- `experiments/validate_spectral_gap.jl` - Updated TARGET_REL_ERROR_N6 comment to document XZ_stagg result

## Decisions Made
- **XZ_stagg has zero n=6 gap-mode overlap:** Despite breaking SU(2) spin-rotation symmetry and having k=pi momentum component, XZ_stagg cannot couple to the n=6 gap mode. The protection mechanism is stronger than simple SU(2) -- likely related to the 3-fold degeneracy structure of the gap eigenspace and additional discrete symmetries of the Heisenberg chain.
- **Keep 12% threshold for n=6:** No improvement from XZ_stagg, so TARGET_REL_ERROR_N6 remains at 0.12. The ~10.7% error comes from ZZ_avg (same as before adding XZ_stagg).
- **n=4 XZ_stagg has zero overlap too:** For n=4, XZ_stagg has |c_gap|=0.000000, consistent with the k=0 gap mode being in a different symmetry sector.

## Deviations from Plan
None - plan executed exactly as written. The plan anticipated that XZ_stagg might not improve n=6 (the action specified conditional threshold updates based on results).

## Issues Encountered
- Git index corruption after Task 1 commit (resolved by `rm .git/index && git reset`)
- experiments/ directory is gitignored; used `git add -f` (documented precedent from Phase 25-03)

## Physics Analysis

The n=6 Heisenberg chain gap mode remains inaccessible to all 8 observables:

| Observable | n=4 |c_gap| | n=6 |c_gap| | Notes |
|------------|-------------|-------------|-------|
| H          | 0.2457      | 0.0000      | Strong for n=4 |
| Mz         | 0.0000      | 0.0000      | Zero (symmetry) |
| XX_avg     | 0.5460      | 0.0000      | Strong for n=4 |
| YY_avg     | 0.5460      | 0.0000      | Strong for n=4 |
| ZZ_avg     | 0.5460      | 0.0000      | Strong for n=4, best estimator |
| Mz_stagg   | 0.0000      | 0.0000      | k=pi but SU(2)-symmetric |
| Z1         | 0.0000      | 0.0000      | All k sectors but SU(2)-symmetric |
| XZ_stagg   | 0.0000      | 0.0000      | k=pi + breaks SU(2), still zero |

The n=6 gap eigenspace (3-fold degenerate, k=pi sector) appears to be protected by a symmetry more specific than translational or SU(2). Possible candidates: parity (reflection symmetry) combined with spin-flip, or a higher symmetry of the isotropic Heisenberg model. Further investigation would require analyzing the gap eigenspace structure under all discrete symmetries of the periodic Heisenberg chain.

## Next Steps
- n=6 gap estimation at ~10.7% is the practical limit with current observable design
- Breaking the symmetry protection likely requires constructing observables that target the specific irrep of the gap eigenspace
- Alternative approach: use non-Hermitian or time-dependent observables

---
*Quick Task: 27-add-symmetry-breaking-observable*
*Completed: 2026-02-18*
