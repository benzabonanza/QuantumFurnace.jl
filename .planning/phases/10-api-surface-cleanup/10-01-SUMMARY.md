---
phase: 10-api-surface-cleanup
plan: 01
subsystem: api
tags: [exports, api-surface, julia-module]

# Dependency graph
requires:
  - phase: 06-dead-code-pruning
    provides: "Dead exports (LindbladianJumpCaches, LiouvLiouv, LindbladWorkspace, non-mutating oft) already removed"
  - phase: 09-type-parameterization
    provides: "Generalized type signatures for exported functions"
provides:
  - "Organized export block with labeled groups"
  - "Physics building block functions exported (trace_distance_h, create_f, trotterize, etc.)"
  - "Internal workspaces and precompute helpers de-exported"
affects: [10-02, 10-03]

# Tech tracking
tech-stack:
  added: []
  patterns: ["Labeled export groups in module file", "Internal functions accessed via QuantumFurnace.name in tests"]

key-files:
  created: []
  modified:
    - src/QuantumFurnace.jl
    - test/test_compilation.jl
    - test/test_cptp.jl
    - test/test_dm_scaling.jl
    - test/test_trajectory_fixes.jl
    - test/test_regression.jl
    - test/test_dm_detailed_balance.jl
    - test/trajectory_validation/run_trajectory_validation.jl
    - test/trajectory_validation/run_convergence_tests.jl

key-decisions:
  - "De-exported internal types and functions require qualified QuantumFurnace.name access in tests (Rule 3 deviation)"

patterns-established:
  - "Export block organized into labeled groups: Types, Simulation, QI Tools, Gibbs, Transition functions, Coherent terms, Pauli/Trotter, Config validation, Log-Sobolev, OFT"
  - "Internal workspace types (OFTCaches, KrausScratch, etc.) accessed via qualified QuantumFurnace.name"

# Metrics
duration: 6min
completed: 2026-02-15
---

# Phase 10 Plan 01: Export Curation Summary

**Reorganized export block with labeled groups, exported 17 physics building block functions, de-exported 12 internal workspace/precompute items**

## Performance

- **Duration:** 6 min
- **Started:** 2026-02-15T13:11:16Z
- **Completed:** 2026-02-15T13:17:25Z
- **Tasks:** 2
- **Files modified:** 9

## Accomplishments
- Export block reorganized into 10 labeled groups (Types: Simulation, Types: Domains, QI Tools, etc.)
- 17 physics functions newly exported: trace_distance_h, trace_distance_nh, trace_norm_h, trace_norm_nh, fidelity, frobenius_norm, is_density_matrix, random_density_matrix, hermitianize!, transform_jumps_to_basis, create_f, create_f_gauss, create_alpha_gauss, check_alpha_skew_symmetry, trotterize, group_hamiltonian_terms, pauli_string_to_matrix
- 12 internal items de-exported: OFTCaches, NUFFTPrefactors, KrausScratch, TrajectoryWorkspace, prepare_oft_nufft_prefactors, prefactor_view, precompute_coherent_terms, precompute_coherent_total_B, precompute_R, precompute_data, jump_contribution!, generate_filename
- API-01 confirmed: dead/deprecated exports (LindbladianJumpCaches, LiouvLiouv, LindbladWorkspace, non-mutating oft) were already removed in Phase 6
- All 224 tests pass including Aqua.jl quality checks

## Task Commits

Each task was committed atomically:

1. **Task 1: Reorganize export block** - `9eecd80` (feat)
2. **Task 2: Update tests for unqualified trace_distance_h** - `6c2243e` (refactor)

## Files Created/Modified
- `src/QuantumFurnace.jl` - Reorganized export block with labeled groups, new physics exports, removed internal exports
- `test/test_compilation.jl` - Qualified precompute_data and KrausScratch references
- `test/test_cptp.jl` - Qualified precompute_data and KrausScratch references
- `test/test_dm_scaling.jl` - Qualified precompute_data, OFTCaches references
- `test/test_trajectory_fixes.jl` - Qualified precompute_data and KrausScratch references
- `test/test_regression.jl` - Qualified precompute_data and KrausScratch references
- `test/test_dm_detailed_balance.jl` - Unqualified trace_distance_h calls
- `test/trajectory_validation/run_trajectory_validation.jl` - Unqualified trace_distance_h, qualified precompute_data/KrausScratch
- `test/trajectory_validation/run_convergence_tests.jl` - Unqualified trace_distance_h, qualified precompute_data/KrausScratch

## Decisions Made
- De-exported internal types/functions require qualified `QuantumFurnace.name` access in test files (consistent with the convention that tests accessing internals must use qualified names)

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Qualified de-exported function calls in test files**
- **Found during:** Task 1 (Export block reorganization)
- **Issue:** Tests used `precompute_data`, `KrausScratch`, and `OFTCaches` without qualification; de-exporting them caused 16 test errors
- **Fix:** Added `QuantumFurnace.` prefix to all unqualified references to de-exported functions in 6 test files (test_compilation.jl, test_cptp.jl, test_dm_scaling.jl, test_trajectory_fixes.jl, test_regression.jl, trajectory validation files)
- **Files modified:** test/test_compilation.jl, test/test_cptp.jl, test/test_dm_scaling.jl, test/test_trajectory_fixes.jl, test/test_regression.jl, test/trajectory_validation/run_trajectory_validation.jl, test/trajectory_validation/run_convergence_tests.jl
- **Verification:** All 224 tests pass
- **Committed in:** 9eecd80 (Task 1 commit) and 6c2243e (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (1 blocking)
**Impact on plan:** Necessary for test correctness after de-exporting internal functions. No scope creep.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Export block curated and organized, ready for Plan 02 (internal function _ prefix renaming) and Plan 03 (documentation)
- De-exported items still accessible via `QuantumFurnace.name` for power users and tests

## Self-Check: PASSED

- SUMMARY.md exists
- Commit 9eecd80 exists (Task 1)
- Commit 6c2243e exists (Task 2)
- src/QuantumFurnace.jl exists with Public API marker
- All 224 tests pass

---
*Phase: 10-api-surface-cleanup*
*Completed: 2026-02-15*
