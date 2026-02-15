---
phase: 07-dry-refactoring
plan: 01
subsystem: core-simulation
tags: [julia, DRY, hermitianize, basis-transform, refactoring]

# Dependency graph
requires:
  - phase: 06-dead-code-pruning
    provides: clean codebase with no dead code
provides:
  - hermitianize! in-place helper in qi_tools.jl
  - transform_jumps_to_basis helper in qi_tools.jl
  - all 15 inline Hermitianization patterns replaced
  - all 3 inline Trotter basis transform comprehensions replaced
affects: [07-02, 08-naming-conventions]

# Tech tracking
tech-stack:
  added: []
  patterns: [hermitianize! for numerical Hermiticity enforcement, transform_jumps_to_basis for eigenbasis transforms]

key-files:
  created: []
  modified:
    - src/qi_tools.jl
    - src/jump_workers.jl
    - src/trajectories.jl
    - src/furnace.jl
    - src/log_sobolev.jl

key-decisions:
  - "hermitianize! modifies rho_next in-place then copyto! for evolving_dm sites (semantically equivalent to original broadcast pattern)"
  - "coherent.jl B_trotter single-jump transforms left untouched (different pattern from DRY-04)"

patterns-established:
  - "hermitianize!(A): canonical way to enforce Hermiticity after numerical accumulation"
  - "transform_jumps_to_basis(jumps, eigvecs): canonical way to transform JumpOp vectors to a new eigenbasis"

# Metrics
duration: 4min
completed: 2026-02-15
---

# Phase 7 Plan 1: Hermitianize and Basis Transform Helpers Summary

**Extracted hermitianize! and transform_jumps_to_basis helpers to qi_tools.jl, replacing 15 inline Hermitianization patterns and 3 Trotter basis transform comprehensions across 5 source files**

## Performance

- **Duration:** 4 min
- **Started:** 2026-02-15T07:18:33Z
- **Completed:** 2026-02-15T07:22:14Z
- **Tasks:** 1
- **Files modified:** 5

## Accomplishments
- Created `hermitianize!(A)` in-place helper replacing 15 inline `0.5 .* (A .+ A')` and `(A + A') / 2` patterns
- Created `transform_jumps_to_basis(jumps, eigvecs)` helper replacing 3 identical JumpOp comprehensions
- Updated call sites across jump_workers.jl (6), trajectories.jl (7), furnace.jl (3), log_sobolev.jl (1)
- All 224 tests pass with zero regressions

## Task Commits

Each task was committed atomically:

1. **Task 1: Create hermitianize! and transform_jumps_to_basis helpers, replace all call sites** - `c8c4478` (feat)

## Files Created/Modified
- `src/qi_tools.jl` - Added hermitianize! and transform_jumps_to_basis function definitions
- `src/jump_workers.jl` - Replaced 6 inline Hermitianization patterns (3x scratch.R, 3x evolving_dm via rho_next)
- `src/trajectories.jl` - Replaced 7 inline patterns (2x scratch.R, 1x scratch.tmp2, 3x rho_mean, 1x B_a) and 1 Trotter transform
- `src/furnace.jl` - Replaced 1 steady_state_dm Hermitianization and 2 Trotter transform comprehensions
- `src/log_sobolev.jl` - Replaced 1 v2_mat Hermitianization pattern

## Decisions Made
- For `evolving_dm .= 0.5 .* (scratch.rho_next .+ scratch.rho_next')` sites: used `hermitianize!(scratch.rho_next); copyto!(evolving_dm, scratch.rho_next)` which is semantically identical
- coherent.jl B_trotter functions left untouched: they transform a single jump's data matrix (returning a plain Matrix), not a vector of JumpOps -- different pattern from the DRY-04 duplication

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- DRY-01 and DRY-04 satisfied; ready for 07-02 (remaining DRY refactoring)
- hermitianize! and transform_jumps_to_basis available as canonical helpers for any future code

## Self-Check: PASSED

All files verified present. Commit c8c4478 verified in git log.

---
*Phase: 07-dry-refactoring*
*Completed: 2026-02-15*
