---
phase: 35-workspace-and-channel-consolidation
plan: 02
subsystem: api
tags: [julia, parametric-structs, workspace, trajectory, type-stability, per-operator-kraus, threading]

# Dependency graph
requires:
  - phase: 35-workspace-and-channel-consolidation
    plan: 01
    provides: "Unified Workspace{S,D,C,T,SC} struct with Krylov, Lindbladian, Thermalize variants"
provides:
  - "Workspace{Trajectory,D,C,T} replacing TrajectoryWorkspace + TrajectoryFramework + PerOperatorKraus"
  - "TrajectoryScratch nested sub-struct for mutable trajectory buffers"
  - "_build_trajectory_workspace factory function for trajectory workspace construction"
  - "_copy_workspace_for_thread helper for per-thread workspace copies"
  - "Flattened per-operator Kraus vectors (Rs, K0s, U_residuals, U_Bs)"
affects: [36-channel-and-simulation-pipeline, 37-file-rename-and-cleanup]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "TrajectoryScratch{CT} nested sub-struct for thread-local mutable buffers"
    - "_build_trajectory_workspace factory avoiding Workspace constructor dispatch conflict"
    - "_copy_workspace_for_thread sharing immutable data with independent scratch"
    - "Flat per-operator vectors (Rs, K0s, U_residuals, U_Bs) replacing PerOperatorKraus struct"
    - "3-arg step_along_trajectory!(psi, ws, rng) replacing 4-arg (psi, fw, ws, rng)"

key-files:
  created: []
  modified:
    - "src/structs.jl"
    - "src/trajectories.jl"
    - "src/QuantumFurnace.jl"
    - "src/krylov_workspace.jl"
    - "src/furnace.jl"
    - "src/convergence.jl"
    - "test/test_allocation.jl"
    - "test/test_workspace_independence.jl"
    - "test/test_compilation.jl"
    - "test/test_cptp.jl"
    - "test/test_threading.jl"
    - "test/test_trajectory_fixes.jl"
    - "test/test_gns_trajectory.jl"
    - "test/test_regression.jl"
    - "test/trajectory_validation/run_convergence_tests.jl"
    - "test/trajectory_validation/run_trajectory_validation.jl"

key-decisions:
  - "_build_trajectory_workspace factory function avoids dispatch conflict with Workspace(Config{Thermalize}) in krylov_workspace.jl"
  - "Per-operator Kraus data stored as flat Vector{Matrix{CT}} (Rs, K0s, U_residuals, U_Bs) eliminating PerOperatorKraus struct"
  - "TrajectoryScratch holds only mutable buffers (jump_oft, psi_tmp, Rpsi, rho_acc)"
  - "Dead fields eliminated: config, precomputed_data, delta_eff from old TrajectoryFramework"
  - "_copy_workspace_for_thread shares immutable physics data but creates fresh TrajectoryScratch per thread"
  - "step_along_trajectory! simplified from 4-arg (psi, fw, ws, rng) to 3-arg (psi, ws, rng)"

# Metrics
duration: 15min
completed: 2026-02-27
---

# Phase 35 Plan 02: Trajectory Workspace Consolidation Summary

**Consolidated TrajectoryWorkspace, TrajectoryFramework, and PerOperatorKraus into Workspace{Trajectory,D,C,T} with flattened per-operator Kraus vectors, nested TrajectoryScratch, and _build_trajectory_workspace factory function**

## Performance

- **Duration:** 15 min
- **Started:** 2026-02-27T09:10:06Z
- **Completed:** 2026-02-27T09:25:06Z
- **Tasks:** 2
- **Files modified:** 16

## Accomplishments
- Consolidated three trajectory-specific types (TrajectoryWorkspace, TrajectoryFramework, PerOperatorKraus) into unified Workspace{Trajectory,D,C,T,TrajectoryScratch{CT}}
- Eliminated 4 dead fields: config, precomputed_data, delta_eff, domain from old TrajectoryFramework
- Flattened per-operator Kraus data from nested PerOperatorKraus struct to flat vectors (Rs, K0s, U_residuals, U_Bs) on workspace
- Simplified step_along_trajectory! from 4-arg to 3-arg form (psi, ws, rng)
- Created _copy_workspace_for_thread helper ensuring thread-safe workspace copies share immutable data
- Added trajectory-specific fields to Workspace struct (ham_or_trott, n_jumps, scaled_prefactor, sigma, Rs, K0s, U_residuals, U_Bs)
- Updated all 16 source and test files, all 1199 tests pass with identical numerical results

## Task Commits

Each task was committed atomically:

1. **Task 1: Define TrajectoryScratch and build Workspace{Trajectory} constructor** - `6cf17a6` (feat)
2. **Task 2: Update all trajectory callers and tests for Workspace{Trajectory}** - `b5e3668` (feat)

## Files Created/Modified
- `src/structs.jl` - Added TrajectoryScratch sub-struct and trajectory-specific fields to Workspace struct
- `src/trajectories.jl` - Complete rewrite: deleted TrajectoryWorkspace/TrajectoryFramework/PerOperatorKraus, added _build_trajectory_workspace factory, _copy_workspace_for_thread, updated step_along_trajectory! to 3-arg form, updated all run_* functions
- `src/QuantumFurnace.jl` - Removed TrajectoryFramework/build_trajectoryframework exports, added TrajectoryScratch export
- `src/krylov_workspace.jl` - Updated both Workspace constructors to pass nothing for new trajectory fields
- `src/furnace.jl` - Updated Workspace{Lindbladian} construction for new field count
- `src/convergence.jl` - Updated convergence/adaptive runners from fw to ws naming
- `test/test_allocation.jl` - Replaced TrajectoryWorkspace/build_trajectoryframework with new patterns
- `test/test_workspace_independence.jl` - Workspace independence using _copy_workspace_for_thread, ws.scratch.psi_tmp access
- `test/test_compilation.jl` - _build_trajectory_workspace replaces build_trajectoryframework, Workspace{Trajectory} type checks
- `test/test_cptp.jl` - CPTP completeness via ws.K0s[a]/ws.Rs[a]/ws.U_residuals[a] flat access
- `test/test_threading.jl` - Per-thread workspaces via _copy_workspace_for_thread
- `test/test_trajectory_fixes.jl` - 3-arg step_along_trajectory!, ws.U_Bs/ws.K0s/ws.Rs access
- `test/test_gns_trajectory.jl` - _build_trajectory_workspace for GNS trajectory tests
- `test/test_regression.jl` - Trajectory regression using new workspace pattern
- `test/trajectory_validation/run_convergence_tests.jl` - Convergence tests using new workspace
- `test/trajectory_validation/run_trajectory_validation.jl` - Cross-validation using new workspace

## Decisions Made
- **Factory function over constructor**: `_build_trajectory_workspace` avoids dispatch ambiguity with `Workspace(Config{Thermalize}, ...)` in krylov_workspace.jl. Both match Config{Thermalize} and Julia cannot disambiguate.
- **Flat per-operator vectors**: Rs, K0s, U_residuals, U_Bs stored as Vector{Matrix{CT}} directly on workspace. No nested PerOperatorKraus indirection. Hot-path access is `ws.Rs[a]` instead of `ws.per_operator[a].R`.
- **TrajectoryScratch for thread safety**: Each thread gets its own TrajectoryScratch via `_copy_workspace_for_thread`, which aliases all immutable fields (Rs, K0s, jumps, etc.) but creates fresh mutable buffers.
- **Dead field elimination**: config (callers have it), precomputed_data (absorbed into flat fields), delta_eff (always equals delta), domain (available via type parameter D) -- none needed in hot path.
- **3-arg step simplification**: The old 4-arg pattern (psi, fw, ws, rng) split read-only framework from mutable scratch. With Workspace{Trajectory} containing both, the split is no longer needed.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None - clean execution. The plan's dispatch constraint warning about _build_trajectory_workspace vs Workspace constructor was well-founded and followed exactly.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- Unified Workspace{S,D,C,T,SC} now covers ALL simulation paths: Krylov, Lindbladian, Thermalize (DM), and Trajectory
- Phase 35 workspace consolidation complete (Plans 01 + 02)
- Ready for Phase 36 (channel and simulation pipeline) or Phase 37 (file rename and cleanup)

## Self-Check: PASSED

- All 16 modified files: FOUND
- Commit 6cf17a6 (Task 1): FOUND
- Commit b5e3668 (Task 2): FOUND
- TrajectoryScratch in structs.jl: FOUND
- _build_trajectory_workspace in trajectories.jl: FOUND
- _copy_workspace_for_thread in trajectories.jl: FOUND
- TrajectoryWorkspace struct deleted: CONFIRMED (0 matches in src/)
- TrajectoryFramework struct deleted: CONFIRMED (0 matches in src/, 1 comment only)
- PerOperatorKraus struct deleted: CONFIRMED (0 matches in src/)
- build_trajectoryframework deleted: CONFIRMED (0 matches in src/, 1 comment only)
- All 1199 tests pass: CONFIRMED

---
*Phase: 35-workspace-and-channel-consolidation*
*Completed: 2026-02-27*
