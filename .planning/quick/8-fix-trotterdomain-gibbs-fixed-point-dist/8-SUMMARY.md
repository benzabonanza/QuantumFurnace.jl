---
phase: quick-8
plan: 01
subsystem: physics
tags: [lindbladian, trotter, basis-transform, gibbs, fixed-point]

# Dependency graph
requires:
  - phase: quick-7
    provides: per-operator Lie-Trotter trajectory splitting
provides:
  - Trotter-basis-transformed jumps for dissipative Liouvillian and trajectory code
  - TrotterDomain Gibbs fixed point distance ~1e-8 (was ~0.004)
affects: [trajectory-validation, dm-reference-tests]

# Tech tracking
tech-stack:
  added: []
  patterns: [trotter-basis-transform-before-nufft-dissipative-loop]

key-files:
  created: []
  modified:
    - src/furnace.jl
    - src/trajectories.jl
    - test/trajectory_validation/run_trajectory_validation.jl

key-decisions:
  - "JumpOp[] typed comprehension to preserve Vector{JumpOp} for TrajectoryFramework compatibility"
  - "Original jumps preserved for coherent B computation (handles own transform internally)"
  - "Trotter basis transform applied at three entry points: construct_lindbladian, run_thermalization, build_trajectoryframework"

patterns-established:
  - "Trotter basis transform: U * A * U' where U = trotter.trafo_from_eigen_to_trotter before any NUFFT element-wise product"

# Metrics
duration: 7min
completed: 2026-02-14
---

# Quick Task 8: Fix TrotterDomain Gibbs Fixed Point Distance Summary

**Trotter eigenbasis transform for dissipative jump operators reduces TrotterDomain Gibbs distance from ~0.004 to ~9e-9**

## Performance

- **Duration:** 7 min
- **Started:** 2026-02-14T15:09:31Z
- **Completed:** 2026-02-14T15:17:02Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments
- Fixed basis mismatch: jump.in_eigenbasis (Hamiltonian eigenbasis) is now transformed to Trotter eigenbasis before the NUFFT element-wise product A .* P
- TrotterDomain Liouvillian fixed point trace distance to Gibbs dropped from ~0.004 to ~9e-9
- Tightened test thresholds: fixed point check 0.01 -> 1e-4, trajectory check 0.02 -> 0.015
- All 220 Pkg.test() tests pass, TVAL-06 trajectory validation passes with tightened thresholds

## Task Commits

Each task was committed atomically:

1. **Task 1: Transform jump operators to Trotter eigenbasis** - `e0ef0fc` (feat)
2. **Task 2: Tighten TrotterDomain test thresholds** - `e4b550c` (feat)

## Files Created/Modified
- `src/furnace.jl` - Added jumps_for_diss with Trotter basis transform in construct_lindbladian and run_thermalization
- `src/trajectories.jl` - Added jumps_for_diss with Trotter basis transform in build_trajectoryframework, stored in framework for step_along_trajectory! use
- `test/trajectory_validation/run_trajectory_validation.jl` - Tightened thresholds and updated comments to reflect corrected domain approximation error

## Decisions Made
- Used `JumpOp[...]` typed comprehension (not bare `[...]`) to ensure `Vector{JumpOp}` type matches `TrajectoryFramework` struct field expectation
- Kept original `jumps` for `precompute_coherent_total_B` and `precompute_coherent_unitary_terms` calls -- these functions handle their own Trotter basis transform internally via `B_trotter()` in coherent.jl
- Transform applied at all three entry points: `construct_lindbladian`, `run_thermalization`, and `build_trajectoryframework`

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Fixed Vector type mismatch in JumpOp comprehension**
- **Found during:** Task 1 (Transform jump operators)
- **Issue:** `[JumpOp(...) for j in jumps]` produced `Vector{JumpOp{Matrix{ComplexF64}}}` but `TrajectoryFramework` requires `Vector{JumpOp}` (unparameterized)
- **Fix:** Used `JumpOp[JumpOp(...) for j in jumps]` to force `Vector{JumpOp}` element type
- **Files modified:** src/furnace.jl, src/trajectories.jl
- **Verification:** CPTP TrotterDomain test passes, all 220 tests pass
- **Committed in:** e0ef0fc (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (1 blocking)
**Impact on plan:** Necessary type fix for Julia's parametric type system. No scope creep.

## Issues Encountered
None beyond the type mismatch documented above.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- TrotterDomain now matches Energy/TimeDomain accuracy level for Gibbs fixed point
- Trajectory validation TVAL-06 confirms convergence with corrected error budget
- Ready for Phase 5 (regression suite) or further trajectory analysis

---
*Quick Task: 8-fix-trotterdomain-gibbs-fixed-point-dist*
*Completed: 2026-02-14*
