---
phase: 13-unify-residual-cholesky
plan: 01
subsystem: simulation
tags: [eigendecomposition, psd, cholesky, numerical-robustness, lindbladian]

# Dependency graph
requires:
  - phase: 07-dry-refactoring
    provides: apply_cptp_channel! helper extracted into jump_workers.jl
provides:
  - Unified eigendecomposition-based U_residual computation in apply_cptp_channel!
  - Eliminated fragile Cholesky + eps_shift hack from DM simulator path
affects: [jump_workers, trajectories, dm-simulation]

# Tech tracking
tech-stack:
  added: []
  patterns: [eigendecomposition with clamped eigenvalues for PSD square root]

key-files:
  created: []
  modified: [src/jump_workers.jl]

key-decisions:
  - "hermitianize!(scratch.tmp2) added before eigen to handle floating-point asymmetry in S matrix"

patterns-established:
  - "PSD square root: hermitianize! -> Hermitian -> eigen -> clamp negatives -> Diagonal(sqrt) * vectors'"

# Metrics
duration: 2min
completed: 2026-02-15
---

# Quick Task 13: Unify Residual Cholesky Computation Summary

**Replaced fragile Cholesky + eps_shift PSD guard in apply_cptp_channel! with eigendecomposition-based clamped square root matching trajectories.jl**

## Performance

- **Duration:** 2 min
- **Started:** 2026-02-15T07:47:19Z
- **Completed:** 2026-02-15T07:49:30Z
- **Tasks:** 1
- **Files modified:** 1

## Accomplishments
- Replaced `cholesky!(Hermitian(scratch.tmp2), check=false)` with `eigen(Hermitian(scratch.tmp2))` + clamped eigenvalues
- Removed the `eps_shift = 10 * eps(Float64)` diagonal hack that could still produce NaN for non-PSD matrices
- DM simulator (jump_workers.jl) and trajectory simulator (trajectories.jl) now use identical PSD square root approach
- All 224 tests pass with identical results

## Task Commits

Each task was committed atomically:

1. **Task 1: Replace Cholesky with eigendecomposition in apply_cptp_channel!** - `2d5ca52` (feat)

## Files Created/Modified
- `src/jump_workers.jl` - Replaced Cholesky-based U_residual with eigendecomposition approach in apply_cptp_channel!

## Decisions Made
- Added `hermitianize!(scratch.tmp2)` before `eigen` call to ensure floating-point symmetry in S matrix, matching the trajectories.jl pattern where this was already done

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Both DM and trajectory simulators now share the same numerically robust PSD square root pattern
- No further unification needed for this computation path
- Ready for continued v1.1 phase execution

## Self-Check: PASSED

- FOUND: src/jump_workers.jl
- FOUND: commit 2d5ca52
- FOUND: 13-SUMMARY.md

---
*Quick Task: 13-unify-residual-cholesky*
*Completed: 2026-02-15*
