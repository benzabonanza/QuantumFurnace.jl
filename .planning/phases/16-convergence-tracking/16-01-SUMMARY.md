---
phase: 16-convergence-tracking
plan: 01
subsystem: simulation
tags: [convergence, trajectory, trace-distance, observables, batch-runner, density-matrix]

# Dependency graph
requires:
  - phase: 13-threading
    provides: "run_trajectories with multi-threaded batch execution and deterministic seeding"
  - phase: 15-data-architecture
    provides: "ExperimentResult, TrajectoryResult, Dict-based BSON serialization pattern"
provides:
  - "ConvergenceData struct for batch-level convergence metrics"
  - "run_trajectories_convergence batch runner with non-overlapping seed management"
  - "build_convergence_observables: ZZ and H observables in eigenbasis"
  - "build_convergence_observables_trotter: ZZ and H observables in Trotter eigenbasis"
  - "_gibbs_in_trotter_basis: Gibbs state transform for TrotterDomain"
  - "_convergence_to_dict / _dict_to_convergence: Dict serialization for BSON"
affects: [17-adaptive-stopping, 18-kms-gns-comparison]

# Tech tracking
tech-stack:
  added: []
  patterns: [batch-accumulate-measure loop, eigenbasis observable construction, seed-offset batching]

key-files:
  created: [src/convergence.jl]
  modified: [src/results.jl, src/QuantumFurnace.jl]

key-decisions:
  - "Observables built in eigenbasis (not computational) to match run_trajectories rho output"
  - "Separate function (run_trajectories_convergence) instead of modifying existing run_trajectories"
  - "Scalar metrics only (no density matrix snapshots) for O(n_batches) memory"
  - "convergence.jl included before results.jl so ConvergenceData type available to serialization"
  - "Pre-allocated storage (not push!) for convergence data arrays for performance"

patterns-established:
  - "Batch-accumulate-measure: run batches via run_trajectories, accumulate rho_acc, measure metrics after each batch"
  - "Seed offset batching: batch_seed = actual_seed + n_total ensures non-overlapping trajectory seeds"
  - "Eigenbasis observable construction: V' * O_comp * V transform for all observables"

# Metrics
duration: 3min
completed: 2026-02-16
---

# Phase 16 Plan 01: Convergence Tracking Infrastructure Summary

**Batch-level convergence monitoring with ConvergenceData struct, eigenbasis ZZ/H observable builders, and run_trajectories_convergence batch runner wrapping existing trajectory infrastructure**

## Performance

- **Duration:** 3 min
- **Started:** 2026-02-16T09:48:38Z
- **Completed:** 2026-02-16T09:52:23Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments
- ConvergenceData struct storing scalar convergence metrics (trace distance, per-observable values) at batch checkpoints
- Observable builders for nearest-neighbor ZZ correlations and energy H, correctly transformed to eigenbasis (EnergyDomain) or Trotter eigenbasis (TrotterDomain)
- Batch convergence runner that wraps existing run_trajectories with non-overlapping seed management, producing (TrajectoryResult, ConvergenceData) tuple
- Dict-based serialization for ConvergenceData following the same pattern as ExperimentResult (forward-compatible with get/default)
- All 364 existing tests pass -- zero regressions

## Task Commits

Each task was committed atomically:

1. **Task 1: ConvergenceData struct, observable builders, and run_trajectories_convergence** - `96b32f3` (feat)
2. **Task 2: ConvergenceData Dict serialization and module integration** - `3ad335c` (feat)

## Files Created/Modified
- `src/convergence.jl` - ConvergenceData struct, build_convergence_observables, build_convergence_observables_trotter, _gibbs_in_trotter_basis, _compute_gibbs_observable_values, run_trajectories_convergence
- `src/results.jl` - _convergence_to_dict, _dict_to_convergence for BSON serialization
- `src/QuantumFurnace.jl` - include(convergence.jl) before results.jl, export public API symbols

## Decisions Made
- Observables built in eigenbasis (not computational basis) because run_trajectories returns rho_mean in eigenbasis for EnergyDomain; using computational basis would give wrong tr(rho * O) values
- Separate run_trajectories_convergence function rather than modifying existing run_trajectories to avoid API bloat (the existing observables parameter serves a different purpose: time-resolved within-trajectory measurements)
- Scalar-only ConvergenceData (no density matrix snapshots) to keep memory O(n_batches) per research Pitfall 11
- convergence.jl included before results.jl in module file so ConvergenceData type is available when results.jl defines serialization functions
- Pre-allocated arrays (Vector/Matrix with undef) instead of push!/hcat for convergence data storage

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Convergence tracking infrastructure complete, ready for Phase 16 Plan 02 (tests)
- Phase 17 (adaptive stopping) can consume ConvergenceData to implement convergence detection
- Phase 18 (KMS-vs-GNS comparison) can use convergence curves for paper figures

## Self-Check: PASSED

All files and commits verified:
- src/convergence.jl: FOUND
- src/results.jl: FOUND
- src/QuantumFurnace.jl: FOUND
- 16-01-SUMMARY.md: FOUND
- 96b32f3 (Task 1): FOUND
- 3ad335c (Task 2): FOUND

---
*Phase: 16-convergence-tracking*
*Completed: 2026-02-16*
