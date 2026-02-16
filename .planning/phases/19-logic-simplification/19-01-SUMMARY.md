---
phase: 19-logic-simplification
plan: 01
subsystem: quantum-simulation
tags: [jumpop, trotter, basis-transform, lindbladian, coherent]

# Dependency graph
requires: []
provides:
  - "Domain-aware JumpOp construction at the source (trotter.eigvecs for TrotterDomain)"
  - "TEST_TROTTER_JUMPS and SMALL_TROTTER_JUMPS test constants"
  - "Zero internal transform_jumps_to_basis callers in src/"
affects: [19-02, 19-03]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Basis selection at JumpOp construction source instead of downstream transform"
    - "Domain-conditional jumps selection in test loops: `domain isa TrotterDomain ? *_TROTTER_JUMPS : *_JUMPS`"

key-files:
  created: []
  modified:
    - "test/test_helpers.jl"
    - "test/test_gns_trajectory.jl"
    - "test/test_regression.jl"
    - "test/test_cptp.jl"
    - "test/test_dm_detailed_balance.jl"
    - "test/test_dm_scaling.jl"
    - "test/reference/generate_references.jl"
    - "test/trajectory_validation/run_trajectory_validation.jl"
    - "test/trajectory_validation/run_convergence_tests.jl"
    - "experiments/run_sweep.jl"
    - "src/trajectories.jl"
    - "src/furnace.jl"
    - "src/coherent.jl"

key-decisions:
  - "Use trotter.eigvecs for TrotterDomain JumpOp basis, hamiltonian.eigvecs for all others"
  - "Keep transform_jumps_to_basis as public API export, remove only internal callers"
  - "Pre-compute SMALL_TROTTER_JUMPS and TEST_TROTTER_JUMPS constants alongside existing ones"

patterns-established:
  - "JumpOp in_eigenbasis field contains domain-appropriate basis from construction"
  - "TrotterDomain test call sites use *_TROTTER_JUMPS constants"

# Metrics
duration: 8min
completed: 2026-02-16
---

# Phase 19 Plan 01: Eliminate Redundant Jump Basis Transforms Summary

**Domain-aware JumpOp construction at the source eliminates 6 internal transform_jumps_to_basis calls across 3 src files, verified by all 539 tests passing**

## Performance

- **Duration:** 8 min
- **Started:** 2026-02-16T14:15:15Z
- **Completed:** 2026-02-16T14:23:47Z
- **Tasks:** 3/3
- **Files modified:** 13

## Accomplishments
- Test helpers produce both hamiltonian-basis and trotter-basis JumpOps via optional `trotter` kwarg
- All TrotterDomain test call sites updated to use pre-built trotter-basis JumpOps (8 test files)
- All 6 internal `transform_jumps_to_basis` calls removed from `trajectories.jl`, `furnace.jl`, `coherent.jl`
- `build_trotter_system()` in sweep script uses `trotter.eigvecs` for JumpOp construction
- `transform_jumps_to_basis` remains exported as public API utility

## Task Commits

Each task was committed atomically:

1. **Task 1: Add trotter-basis JumpOp construction to test helpers** - `4c20ade` (feat)
2. **Task 2: Update TrotterDomain test call sites to use trotter-basis JumpOps** - `8bb3db6` (feat)
3. **Task 3: Remove internal transform_jumps_to_basis calls and update sweep script** - `e54aa56` (refactor)

## Files Created/Modified
- `test/test_helpers.jl` - Added optional `trotter` kwarg to make_test_system/make_small_test_system; added TEST_TROTTER_JUMPS/SMALL_TROTTER_JUMPS constants
- `test/test_gns_trajectory.jl` - 4 TrotterDomain call sites use SMALL_TROTTER_JUMPS
- `test/test_regression.jl` - 3 TrotterDomain call sites use SMALL_TROTTER_JUMPS
- `test/test_cptp.jl` - 1 TrotterDomain call site uses TEST_TROTTER_JUMPS
- `test/test_dm_detailed_balance.jl` - Domain loop selects jumps by domain type
- `test/test_dm_scaling.jl` - Uses TEST_TROTTER_JUMPS[1] instead of manual transform_jumps_to_basis
- `test/reference/generate_references.jl` - Selects jumps by domain type
- `test/trajectory_validation/run_trajectory_validation.jl` - Helper + 2 explicit sites use SMALL_TROTTER_JUMPS
- `test/trajectory_validation/run_convergence_tests.jl` - Helper selects jumps by domain type
- `experiments/run_sweep.jl` - build_trotter_system uses trotter.eigvecs for JumpOp construction
- `src/trajectories.jl` - Removed conditional transform in build_trajectoryframework
- `src/furnace.jl` - Removed conditional transforms in construct_lindbladian and run_thermalization
- `src/coherent.jl` - Removed transform_jumps_to_basis from all 3 precompute functions

## Decisions Made
- Used `trotter.eigvecs` for TrotterDomain and `hamiltonian.eigvecs` for all other domains (per user LOCKED decision)
- Kept `transform_jumps_to_basis` as public API export -- only removed internal automatic callers
- Pre-computed trotter-basis JumpOp constants at test include time (same pattern as existing constants)

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- All internal transform_jumps_to_basis usage eliminated from src/
- Test infrastructure provides both basis variants via constants
- Ready for plans 02 and 03 (further logic simplification)

## Self-Check: PASSED

All 13 modified files verified present. All 3 task commits verified in git log. 539/539 tests pass.

---
*Phase: 19-logic-simplification*
*Completed: 2026-02-16*
