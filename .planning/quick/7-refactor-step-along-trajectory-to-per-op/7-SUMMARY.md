---
phase: quick-7
plan: 01
subsystem: trajectories
tags: [quantum-trajectories, lie-trotter, per-operator, kraus, cptp]

# Dependency graph
requires:
  - phase: 02-trajectory-bug-fixes
    provides: Correct trajectory stepping (TFIX-02/03/04/05), CPTP completeness test
provides:
  - PerOperatorKraus struct for per-operator Lie-Trotter splitting
  - Refactored TrajectoryFramework with per_operator vector and delta_eff
  - Per-operator step_along_trajectory! matching DM run_thermalization structure
affects: [05-regression, trajectory-validation]

# Tech tracking
tech-stack:
  added: []
  patterns: [per-operator-lie-trotter-splitting, rate-rescaling-by-inv-p-jump]

key-files:
  created: []
  modified:
    - src/trajectories.jl
    - test/test_cptp.jl
    - test/test_compilation.jl
    - test/test_trajectory_fixes.jl

key-decisions:
  - "delta_eff = delta * N_jumps for per-operator rate rescaling; alpha derived from delta_eff"
  - "Per-operator coherent B computed via precompute_coherent_total_B with single-jump vector"
  - "JumpOp[jumps[a]] used to force Vector{JumpOp} type for dispatch compatibility with coherent_bohr"

patterns-established:
  - "Per-operator pattern: fw.per_operator[a].R/K0/U_residual/U_B instead of single fw.R/K0/U_residual/U_B"
  - "Rate rescaling: 1/p_jump baked into per-operator R_a and jump probability prefactor"

# Metrics
duration: 10min
completed: 2026-02-14
---

# Quick Task 7: Per-Operator Lie-Trotter Trajectory Splitting Summary

**Refactored trajectory simulator from single combined CPTP channel to per-operator random selection matching DM run_thermalization structure**

## Performance

- **Duration:** 10 min
- **Started:** 2026-02-14T14:33:51Z
- **Completed:** 2026-02-14T14:43:54Z
- **Tasks:** 3
- **Files modified:** 4

## Accomplishments
- Per-operator Kraus data (PerOperatorKraus struct) with R_a, K0_a, U_residual_a, U_B_a for each jump operator
- step_along_trajectory! randomly selects one operator per step (both EnergyDomain and Time/TrotterDomain variants)
- CPTP completeness verified for each per-operator channel individually across all three domains (39 test assertions)
- Full test suite passes: 220 tests, 0 failures

## Task Commits

Each task was committed atomically:

1. **Task 1: Add PerOperatorKraus struct and refactor TrajectoryFramework** - `7be1ade` (feat)
2. **Task 2: Refactor step_along_trajectory! for per-operator branching** - `7080046` (feat)
3. **Task 3: Update CPTP and trajectory tests for per-operator verification** - `0c879e7` (feat)

## Files Created/Modified
- `src/trajectories.jl` - PerOperatorKraus struct, refactored TrajectoryFramework, per-operator build and step functions
- `test/test_cptp.jl` - Per-operator CPTP completeness test (K0_a'K0_a + delta_eff*R_a + U_res_a'U_res_a = I)
- `test/test_compilation.jl` - Updated to test per_operator U_B presence instead of single fw.B/fw.U_B
- `test/test_trajectory_fixes.jl` - Updated all field references from fw.R/K0/U_residual to per_operator access

## Decisions Made
- delta_eff = delta * N_jumps for per-operator rate rescaling (matches DM run_thermalization's 1/p_jump convention)
- Per-operator coherent B computed by calling precompute_coherent_total_B with single-jump vector, not precompute_coherent_unitary_terms (avoids needing HamHam separately for TrotterDomain)
- JumpOp[jumps[a]] syntax used to force Vector{JumpOp} type for dispatch compatibility with coherent_bohr/B_time/B_trotter which have concrete Vector{JumpOp} signatures

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Fixed JumpOp vector type dispatch for coherent term computation**
- **Found during:** Task 3 (test verification)
- **Issue:** `[jumps[a]]` created `Vector{JumpOp{Matrix{ComplexF64}}}` which didn't match `coherent_bohr`'s `Vector{JumpOp}` signature
- **Fix:** Used `JumpOp[jumps[a]]` to force correct container type
- **Files modified:** src/trajectories.jl
- **Verification:** All compilation tests pass, coherent term computed correctly
- **Committed in:** 0c879e7 (Task 3 commit)

**2. [Rule 3 - Blocking] Updated test_compilation.jl and test_trajectory_fixes.jl for new struct fields**
- **Found during:** Task 3 (test verification)
- **Issue:** Tests referenced removed fields (fw.R, fw.K0, fw.U_residual, fw.B, fw.U_B) causing 10 test errors
- **Fix:** Updated all field references to use fw.per_operator[a].R/K0/U_residual/U_B
- **Files modified:** test/test_compilation.jl, test/test_trajectory_fixes.jl
- **Verification:** Full test suite passes (220 tests)
- **Committed in:** 0c879e7 (Task 3 commit)

---

**Total deviations:** 2 auto-fixed (2 blocking)
**Impact on plan:** Both auto-fixes necessary for correctness. No scope creep. Plan anticipated test_cptp.jl changes but not the other two test files.

## Issues Encountered
None beyond the blocking issues documented above.

## Next Phase Readiness
- Trajectory simulator now uses per-operator Lie-Trotter splitting matching DM code structure
- Phase 5 (regression) can proceed with this as foundation
- Trajectory cross-validation tests (run_trajectory_validation.jl) were NOT updated -- they only call step_along_trajectory! and build_trajectoryframework without accessing internal fields, so they work as-is. However, they will now validate the per-operator splitting behavior instead of the combined channel behavior.

## Self-Check: PASSED

All files, commits, and key artifacts verified:
- 4/4 files found (src/trajectories.jl, test/test_cptp.jl, test/test_compilation.jl, test/test_trajectory_fixes.jl)
- 3/3 commits found (7be1ade, 7080046, 0c879e7)
- PerOperatorKraus struct defined (1 occurrence)
- per_operator field used (6 occurrences in trajectories.jl, 3 in test_cptp.jl)
- rand(1:fw.n_jumps) present in both step variants (2 occurrences)

---
*Quick Task: 7-refactor-step-along-trajectory-to-per-op*
*Completed: 2026-02-14*
