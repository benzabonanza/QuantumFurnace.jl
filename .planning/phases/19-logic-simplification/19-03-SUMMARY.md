---
phase: 19-logic-simplification
plan: 03
subsystem: simulation
tags: [structs, result-types, convergence, serialization, bson]

# Dependency graph
requires:
  - phase: 19-02
    provides: "Flattened call chain with _build_framework_and_seed, _run_batch_no_obs!"
provides:
  - "LindbladianResult replacing HotSpectralResults (4 fields, no hamiltonian/config baggage)"
  - "DMSimulationResult replacing HotAlgorithmResults (3 fields, no hamiltonian/config baggage)"
  - "TrajectoryResult.convergence field embedding ConvergenceData"
  - "Single-value returns from run_trajectories_convergence and run_trajectories_adaptive"
  - "Backward-compatible BSON serialization for new struct layouts"
affects: [experiments, results, convergence]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Embed convergence data inside result struct instead of returning tuple"
    - "Slim result structs: carry only output data, not input config/hamiltonian"

key-files:
  created: []
  modified:
    - "src/structs.jl"
    - "src/trajectories.jl"
    - "src/convergence.jl"
    - "src/furnace.jl"
    - "src/results.jl"
    - "src/QuantumFurnace.jl"
    - "src/log_sobolev.jl"
    - "test/test_results.jl"
    - "test/test_convergence.jl"
    - "test/test_threading.jl"
    - "test/test_allocation.jl"
    - "experiments/run_sweep.jl"

key-decisions:
  - "ConvergenceData struct moved from convergence.jl to structs.jl for include-order safety (TrajectoryResult needs it before convergence.jl is loaded)"
  - "data -> liouvillian, evolved_dm -> final_dm, distances_to_gibbs -> trace_distances field renames for clarity"
  - "Convergence runners return single TrajectoryResult with .convergence field (not tuple)"
  - "test_threading.jl _evolve_along_trajectory! references replaced with inline step loop (function was deleted in 19-02)"

patterns-established:
  - "Result structs carry only output data, not input references"
  - "Optional sub-results embedded as Union{Nothing, T} fields"

# Metrics
duration: 7min
completed: 2026-02-16
---

# Phase 19 Plan 03: Simplify Result Struct Hierarchy Summary

**LindbladianResult/DMSimulationResult replacing heavy result structs, ConvergenceData embedded in TrajectoryResult with single-value returns**

## Performance

- **Duration:** 7 min
- **Started:** 2026-02-16T14:36:38Z
- **Completed:** 2026-02-16T14:44:03Z
- **Tasks:** 3
- **Files modified:** 12

## Accomplishments
- HotSpectralResults renamed to LindbladianResult (4 fields: liouvillian, fixed_point, gap_mode, spectral_gap -- no hamiltonian/config/trotter baggage)
- HotAlgorithmResults renamed to DMSimulationResult (3 fields: final_dm, trace_distances, time_steps -- no hamiltonian/config/trotter baggage)
- TrajectoryResult gains convergence::Union{Nothing, ConvergenceData} field
- Convergence/adaptive runners return single TrajectoryResult (not tuple)
- BSON backward compatibility preserved (missing convergence key loads as nothing)
- All 539 tests pass

## Task Commits

Each task was committed atomically:

1. **Task 1a: Rename and slim result structs** - `332c54d` (refactor)
2. **Task 1b: Embed ConvergenceData in TrajectoryResult and update serialization** - `8746fd4` (feat)
3. **Task 2: Update test files and sweep script for new struct shapes** - `a2befc6` (test)

## Files Created/Modified
- `src/structs.jl` - LindbladianResult and DMSimulationResult replacing old heavy structs; ConvergenceData moved here for include-order safety
- `src/trajectories.jl` - TrajectoryResult with 6th convergence field, constructors updated
- `src/convergence.jl` - Single-value returns from run_trajectories_convergence and run_trajectories_adaptive; ConvergenceData struct removed (now in structs.jl)
- `src/furnace.jl` - Updated result construction for LindbladianResult and DMSimulationResult
- `src/results.jl` - Updated _trajectory_to_dict and _dict_to_experiment for convergence field with backward compatibility
- `src/QuantumFurnace.jl` - Updated exports (LindbladianResult, DMSimulationResult)
- `src/log_sobolev.jl` - Updated compute_LSI_alpha2 type signature and field access (result.data -> result.liouvillian)
- `test/test_results.jl` - TrajectoryResult constructors gain 6th arg
- `test/test_convergence.jl` - All tuple destructuring changed to .convergence field access
- `test/test_threading.jl` - Replaced deleted _evolve_along_trajectory! with inline step loop
- `test/test_allocation.jl` - Removed _evolve_along_trajectory! import
- `experiments/run_sweep.jl` - Single-value return from run_trajectories_adaptive

## Decisions Made
- ConvergenceData struct moved from convergence.jl to structs.jl: TrajectoryResult in trajectories.jl (included at line 97) references ConvergenceData, but convergence.jl was included at line 103. Moving the struct definition to structs.jl (line 85) resolves the ordering dependency.
- Field renames chosen for clarity: `data` -> `liouvillian` (describes the matrix), `evolved_dm` -> `final_dm` (emphasizes it's the final state), `distances_to_gibbs` -> `trace_distances` (shorter, standard term)
- test_threading.jl references to `_evolve_along_trajectory!` (deleted in 19-02) replaced with inline step loop matching the pattern in `_run_chunk_no_obs!`

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] ConvergenceData include-order dependency**
- **Found during:** Task 1b (Embed ConvergenceData in TrajectoryResult)
- **Issue:** TrajectoryResult in trajectories.jl (included before convergence.jl) references ConvergenceData type, causing UndefVarError
- **Fix:** Moved ConvergenceData struct and backward-compatible constructor from convergence.jl to structs.jl (included before trajectories.jl)
- **Files modified:** src/structs.jl, src/convergence.jl
- **Verification:** Module compiles successfully, all 539 tests pass
- **Committed in:** 8746fd4 (Task 1b commit)

**2. [Rule 3 - Blocking] Stale _evolve_along_trajectory! references in test_threading.jl and test_allocation.jl**
- **Found during:** Task 2 (Update test files)
- **Issue:** test_threading.jl called _evolve_along_trajectory! which was deleted in Phase 19 Plan 02; test_allocation.jl imported it
- **Fix:** Replaced calls with inline step loop (normalize + loop step_along_trajectory!); removed import
- **Files modified:** test/test_threading.jl, test/test_allocation.jl
- **Verification:** All 539 tests pass including threading tests
- **Committed in:** a2befc6 (Task 2 commit)

---

**Total deviations:** 2 auto-fixed (2 blocking issues)
**Impact on plan:** Both auto-fixes necessary for compilation/correctness. No scope creep.

## Issues Encountered
None beyond the auto-fixed deviations above.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Phase 19 (logic-simplification) is now complete (all 3 plans executed)
- Clean result struct hierarchy with per-method types
- All 539 tests pass

---
*Phase: 19-logic-simplification*
*Completed: 2026-02-16*
