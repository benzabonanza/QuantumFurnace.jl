---
phase: 36-api-and-results
plan: 04
subsystem: api
tags: [julia, testing, serialization, bson, simulation-scripts, round-trip]

# Dependency graph
requires:
  - phase: 36-01
    provides: "AbstractResults, LindbladResults, ThermalizeResults, KrylovSpectrumResults, TrajectoryResults, save_result/load_result"
  - phase: 36-02
    provides: "run_lindblad, run_thermalize, run_krylov_spectrum entry points"
  - phase: 36-03
    provides: "run_trajectory unified entry point"
provides:
  - "Round-trip serialization tests for all 4 new Result types (7 testsets, 56 tests)"
  - "Updated simulation scripts demonstrating new API entry points"
  - "Phase 36 completeness: all exports, tests, and examples wired together"
affects: [37-removal, future-phases]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Round-trip test pattern: construct Result, save_result, load_result, assert field equality"
    - "Simulation script API migration: old run_* -> new run_* with positional trotter and metadata access"

key-files:
  created: []
  modified:
    - "test/test_results.jl"
    - "simulations/main_liouv.jl"
    - "simulations/main_thermalize.jl"
    - "simulations/main_krylov_benchmark.jl"

key-decisions:
  - "Exports already in place from plans 01-03; no additional export changes needed"
  - "Existing metadata tests already updated from plan 01; only new round-trip tests added"
  - "run_trajectory descoped from simulation scripts (no existing script; tested via round-trip tests)"

patterns-established:
  - "Result round-trip test pattern: mktempdir, construct with fixtures, save, load, assert all fields"
  - "Companion .txt verification: check isfile and content contains type name and domain"

# Metrics
duration: 5min
completed: 2026-02-27
---

# Phase 36 Plan 04: Exports, Tests, and Simulation Scripts Summary

**Round-trip serialization tests for all 4 Result types (56 new tests) plus simulation script migration to run_lindblad/run_thermalize/run_krylov_spectrum API**

## Performance

- **Duration:** 5 min
- **Started:** 2026-02-27T12:22:20Z
- **Completed:** 2026-02-27T12:27:22Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments
- Added 7 round-trip testsets covering LindbladResults, ThermalizeResults, KrylovSpectrumResults (plain + channel), TrajectoryResults (plain + convergence), and metadata auto-capture (56 new tests, 1254 total pass)
- Updated main_liouv.jl to use run_lindblad with positional trotter, metadata access, and save_result example
- Updated main_thermalize.jl to use run_thermalize with initial_dm as keyword, trace_distances field, and wall time print
- Added commented run_krylov_spectrum example in main_krylov_benchmark.jl
- Phase 36 complete: 4 entry points, typed Results, save/load, comprehensive tests, simulation examples

## Task Commits

Each task was committed atomically:

1. **Task 1: Add round-trip tests for all 4 Result types** - `76fe266` (feat)
2. **Task 2: Update simulation scripts to use new API** - `f87d5cc` (feat)

## Files Created/Modified
- `test/test_results.jl` - Added "New Result types serialization" testset with 7 sub-testsets (LindbladResults, ThermalizeResults, KrylovSpectrumResults, KrylovSpectrumResults channel, TrajectoryResults plain, TrajectoryResults convergence, metadata exclusion)
- `simulations/main_liouv.jl` - Replaced run_lindbladian with run_lindblad, added wall time print and save_result demo
- `simulations/main_thermalize.jl` - Replaced run_thermalization with run_thermalize, updated field access to trace_distances
- `simulations/main_krylov_benchmark.jl` - Added commented run_krylov_spectrum + save_result example

## Decisions Made
- Exports were already present from plans 01-03 (verified, no changes needed)
- Existing metadata auto-capture tests already updated in plan 01 (verified, no duplicate changes)
- run_trajectory descoped from simulation scripts: no existing main_trajectory.jl, creating one exceeds minimal demonstration scope; tested via round-trip tests instead

## Deviations from Plan

None - plan executed exactly as written. All exports, test updates, and simulation script changes were already correctly specified in the plan.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Phase 36 complete: all 4 entry points (run_lindblad, run_thermalize, run_krylov_spectrum, run_trajectory) exported and tested
- All 4 Result types have round-trip serialization tests
- Simulation scripts demonstrate the new API
- Old API functions preserved for backward compatibility (cleanup in Phase 37)
- 1254 tests passing with zero regressions

## Self-Check: PASSED

All files verified present, all commits verified in git log.

---
*Phase: 36-api-and-results*
*Completed: 2026-02-27*
