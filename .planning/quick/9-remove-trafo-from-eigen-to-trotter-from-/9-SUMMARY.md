---
phase: quick-9
plan: 01
subsystem: quantum-simulation
tags: [julia, TrottTrott, basis-transform, refactoring, eigvecs]

# Dependency graph
requires:
  - phase: quick-8
    provides: "TrotterDomain Gibbs fixed point basis fix with eigvecs transforms"
provides:
  - "TrottTrott struct without redundant trafo_from_eigen_to_trotter field"
  - "All call sites use trotter.eigvecs directly for basis transforms"
affects: [trotter-domain, lindbladian, trajectories, coherent-terms]

# Tech tracking
tech-stack:
  added: []
  patterns: ["Direct eigvecs transform: trotter.eigvecs' * jump.data * trotter.eigvecs"]

key-files:
  created: []
  modified:
    - src/trotter_domain.jl
    - src/furnace.jl
    - src/coherent.jl
    - src/trajectories.jl
    - src/ofts.jl
    - test/test_dm_scaling.jl

key-decisions:
  - "Use trotter.eigvecs' * j.data * trotter.eigvecs (computational -> Trotter) instead of U * j.in_eigenbasis * U' (H-eigen -> Trotter)"
  - "Test basis transforms use U_t2e = trotter.eigvecs' * ham.eigvecs for Trotter-to-H-eigen conversion"

patterns-established:
  - "Direct eigvecs pattern: always transform from computational basis via trotter.eigvecs, never via intermediate H-eigenbasis"

# Metrics
duration: 3min
completed: 2026-02-14
---

# Quick Task 9: Remove trafo_from_eigen_to_trotter Summary

**Removed redundant TrottTrott field, replacing all H-eigen-to-Trotter transforms with direct eigvecs products**

## Performance

- **Duration:** 3 min
- **Started:** 2026-02-14T15:28:36Z
- **Completed:** 2026-02-14T15:31:46Z
- **Tasks:** 2
- **Files modified:** 6

## Accomplishments
- Removed `trafo_from_eigen_to_trotter` field from TrottTrott struct (5 fields remaining)
- Updated all 5 source call sites to use `trotter.eigvecs' * j.data * trotter.eigvecs`
- Updated 3 test sections (DMTST-05, DMTST-06, DMTST-06b) to use explicit eigvec products
- All 220 tests pass with identical numerical results

## Task Commits

Each task was committed atomically:

1. **Task 1: Remove field from struct/constructor and update all source call sites** - `b9217bf` (feat)
2. **Task 2: Update test code and verify all tests pass** - `5646187` (test)

## Files Created/Modified
- `src/trotter_domain.jl` - Removed trafo field from struct and constructor
- `src/furnace.jl` - Updated construct_lindbladian and run_thermalization jump transforms
- `src/coherent.jl` - Updated B_trotter single-jump and multi-jump variants
- `src/trajectories.jl` - Updated build_trajectoryframework jump transforms
- `src/ofts.jl` - Updated comment referencing the removed field
- `test/test_dm_scaling.jl` - Updated DMTST-05, DMTST-06, DMTST-06b basis transforms

## Decisions Made
- Used `trotter.eigvecs' * j.data * trotter.eigvecs` (computational basis to Trotter eigenbasis directly) instead of the previous two-step path through H-eigenbasis. Mathematically equivalent but simpler.
- For test back-transforms (Trotter-eigen to H-eigen), used `U_t2e = trotter.eigvecs' * ham.eigvecs` explicitly, which is exactly the matrix that was stored in the removed field.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None

## Next Phase Readiness
- TrottTrott struct is now cleaner with 5 fields
- All basis transforms are explicit and use available eigenvectors directly
- No downstream code depends on the removed field

---
*Quick Task: 9*
*Completed: 2026-02-14*
