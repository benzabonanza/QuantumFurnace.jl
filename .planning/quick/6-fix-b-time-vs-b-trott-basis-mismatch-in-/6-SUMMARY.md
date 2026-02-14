---
phase: quick-6
plan: 01
subsystem: quantum-simulation
tags: [trotter, coherent-term, basis-transformation, eigenbasis]

requires:
  - phase: quick-3
    provides: "Pattern for trafo_from_eigen_to_trotter basis fix in OFT functions"
provides:
  - "Corrected B_trotter() single-jump and multi-jump with proper Trotter eigenbasis jump operators"
  - "Tightened DMTST-05 threshold from 0.02 to 1e-5"
affects: [phase-04, trajectory-validation]

tech-stack:
  added: []
  patterns: ["Transform jump.in_eigenbasis to Trotter eigenbasis via trafo_from_eigen_to_trotter before mixing with Trotter time evolution operators"]

key-files:
  created: []
  modified:
    - src/coherent.jl
    - test/test_dm_scaling.jl

key-decisions:
  - "Only fixed active B_trotter() functions; left coherent_term_trotter() unfixed as it is dead code (all call sites commented out)"

patterns-established:
  - "Trotter eigenbasis consistency: any function using Trotter time evolution Diagonal(eigvals_t0.^n) must have jump operators in Trotter eigenbasis, not Hamiltonian eigenbasis"

duration: 2min
completed: 2026-02-14
---

# Quick Task 6: Fix B_trotter() Basis Mismatch Summary

**Fixed B_trotter() basis mismatch by transforming jump operators from H-eigenbasis to Trotter eigenbasis, reducing dist_bohr_trott from ~0.011 to ~1.2e-10**

## Performance

- **Duration:** 2 min
- **Started:** 2026-02-14T11:28:32Z
- **Completed:** 2026-02-14T11:30:47Z
- **Tasks:** 1
- **Files modified:** 2

## Accomplishments
- Fixed basis mismatch bug in both single-jump and multi-jump B_trotter() functions in src/coherent.jl
- Reduced dist_bohr_trott from ~0.011 to ~1.2e-10 (5 orders of magnitude improvement)
- Tightened DMTST-05 test threshold from 0.02 to 1e-5
- All 76 tests pass with no regressions

## Task Commits

Each task was committed atomically:

1. **Task 1: Fix basis mismatch in both B_trotter() functions and tighten DMTST-05 threshold** - `2334a5e` (fix)

## Files Created/Modified
- `src/coherent.jl` - Added trafo_from_eigen_to_trotter basis transformation for jump operators in both B_trotter() variants
- `test/test_dm_scaling.jl` - Tightened DMTST-05 threshold from 0.02 to 1e-5

## Decisions Made
- Only fixed active B_trotter() functions; coherent_term_trotter() has the same bug but is dead code (all call sites commented out) -- fixing it risks regressions with no benefit

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- All coherent term functions now produce consistent results across Bohr, Time, and Trotter domains
- B_trotter() fix matches the OFT fix pattern from quick task 3 (same trafo_from_eigen_to_trotter transformation)
- Ready for Phase 4 trajectory validation work

---
*Quick Task: 6*
*Completed: 2026-02-14*
