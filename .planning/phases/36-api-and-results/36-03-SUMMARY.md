---
phase: 36-api-and-results
plan: 03
subsystem: api
tags: [julia, trajectories, entry-point, keyword-dispatch, convergence, adaptive]

# Dependency graph
requires:
  - phase: 36-01
    provides: "TrajectoryResults struct with metadata field"
  - phase: 35-02
    provides: "Workspace-based trajectory infrastructure (_build_framework_and_seed, _run_batch_no_obs!)"
provides:
  - "run_trajectory unified entry point consolidating 4 trajectory functions"
  - "_run_trajectory_with_obs internal helper for observable path"
  - "_run_trajectory_convergence internal helper for fixed-batch convergence"
  - "_run_trajectory_adaptive internal helper for adaptive early stopping"
affects: [36-04, testing, playground-scripts]

# Tech tracking
tech-stack:
  added: []
  patterns: ["keyword-driven mode dispatch", "internal helper delegation to existing functions", "NamedTuple bridge between old and new return types"]

key-files:
  created: []
  modified:
    - src/trajectories.jl
    - src/convergence.jl
    - src/QuantumFurnace.jl

key-decisions:
  - "Convergence helpers placed in convergence.jl (not trajectories.jl) since they reference functions defined there"
  - "Julia late binding resolves cross-file function calls at runtime, allowing run_trajectory in trajectories.jl to call helpers in convergence.jl"
  - "_run_trajectory_with_obs returns NamedTuple bridge to decouple from ObservableTrajectoryResult internal structure"

patterns-established:
  - "Unified entry point pattern: single function with keyword-driven dispatch to internal helpers"
  - "NamedTuple bridge: internal helpers return NamedTuples that the entry point wraps into typed Results"

# Metrics
duration: 4min
completed: 2026-02-27
---

# Phase 36 Plan 03: Unified run_trajectory Entry Point Summary

**Unified run_trajectory with keyword-driven dispatch across default/observable/convergence/adaptive modes, returning TrajectoryResults with wall_time metadata**

## Performance

- **Duration:** 4 min
- **Started:** 2026-02-27T12:15:37Z
- **Completed:** 2026-02-27T12:20:11Z
- **Tasks:** 1
- **Files modified:** 3

## Accomplishments
- Consolidated 4 existing trajectory functions (run_trajectories, run_observable_trajectories, run_trajectories_convergence, run_trajectories_adaptive) into a single `run_trajectory` entry point
- Keyword-driven mode switching: default (plain DM), observables, convergence, and adaptive modes
- All modes return `TrajectoryResults` with wall_time metadata via `_capture_metadata`
- psi0 as required keyword argument; uniform positional signature (jumps, config, hamiltonian, trotter)
- Convergence modes auto-build observables via `build_preset_trajectory_observables` when none provided
- Gibbs reference auto-computed from `hamiltonian.gibbs` (Trotter-transformed when applicable)

## Task Commits

Each task was committed atomically:

1. **Task 1: Create the unified run_trajectory function** - `32e2b08` (feat)

## Files Created/Modified
- `src/trajectories.jl` - Added `_run_trajectory_with_obs` helper and `run_trajectory` entry point (182 lines)
- `src/convergence.jl` - Added `_run_trajectory_convergence` and `_run_trajectory_adaptive` helpers (84 lines)
- `src/QuantumFurnace.jl` - Added `run_trajectory` to module exports

## Decisions Made
- Convergence helpers placed in `convergence.jl` rather than `trajectories.jl` because they reference `build_preset_trajectory_observables`, `_compute_gibbs_observable_values`, and other convergence functions defined there
- Julia's late binding (runtime method resolution) means `run_trajectory` in `trajectories.jl` can safely call `_run_trajectory_convergence` from `convergence.jl` even though it's included later
- `_run_trajectory_with_obs` returns a NamedTuple rather than the raw `ObservableTrajectoryResult`, providing a clean interface between the old return type and the new `TrajectoryResults` wrapper

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- `run_trajectory` is ready for use alongside `run_lindblad` and `run_thermalize` from 36-02
- Plan 36-04 (run_krylov_spectrum) completes the public API surface
- All 4 existing trajectory functions preserved for backward compatibility

## Self-Check: PASSED

All files verified present. Commit 32e2b08 confirmed in git log. SUMMARY.md created.

---
*Phase: 36-api-and-results*
*Completed: 2026-02-27*
