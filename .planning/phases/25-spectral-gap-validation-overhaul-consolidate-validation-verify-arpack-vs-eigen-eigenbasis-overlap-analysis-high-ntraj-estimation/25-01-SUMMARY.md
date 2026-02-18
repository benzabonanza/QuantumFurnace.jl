---
phase: 25-spectral-gap-validation-overhaul
plan: 01
subsystem: estimation
tags: [observable-builders, spectral-gap, convergence, cleanup]

# Dependency graph
requires:
  - phase: 24-cross-validation
    provides: CrossValidationResult and cross_validate_gap (now removed)
  - phase: 20-observable-builders
    provides: build_gap_estimation_observables (now renamed)
provides:
  - Single observable builder function: build_preset_trajectory_observables
  - Clean codebase without cross-validation wrappers or old experiment scripts
affects: [25-02, 25-03]

# Tech tracking
tech-stack:
  added: []
  patterns: [single-builder-pattern]

key-files:
  created: []
  modified:
    - src/convergence.jl
    - src/gap_estimation.jl
    - src/QuantumFurnace.jl
    - test/test_convergence.jl
    - test/test_gap_estimation.jl

key-decisions:
  - "Single observable builder (build_preset_trajectory_observables) replaces 4 old builders"
  - "Mz construction inlined into the single builder (was delegated to deleted build_total_magnetization)"
  - "CrossValidationResult and cross_validate_gap removed from source and exports"

patterns-established:
  - "Single builder pattern: all trajectory observables come from build_preset_trajectory_observables"

# Metrics
duration: 7min
completed: 2026-02-18
---

# Phase 25 Plan 01: Consolidation Summary

**Consolidated 4 observable builders into single build_preset_trajectory_observables, removed cross-validation code and 4 experiment scripts**

## Performance

- **Duration:** 7 min
- **Started:** 2026-02-18T07:55:17Z
- **Completed:** 2026-02-18T08:02:34Z
- **Tasks:** 2
- **Files modified:** 5 source/test files, 4 deleted

## Accomplishments
- Deleted 4 old experiment scripts (validate_gap_estimation.jl, run_gap_validation.jl, run_sweep.jl, eigenmode_decomposition.jl)
- Removed CrossValidationResult struct and cross_validate_gap functions (2 methods)
- Consolidated build_convergence_observables, build_convergence_observables_trotter, build_total_magnetization, and build_gap_estimation_observables into single build_preset_trajectory_observables
- Updated all test assertions for new observable ordering (H at index 1, not last)
- All 653 tests pass after consolidation

## Task Commits

Each task was committed atomically:

1. **Task 1: Delete old experiment scripts and remove cross-validation code** - `04b56eb` (chore)
2. **Task 2: Consolidate observable builders and update all references** - `36b22cd` (feat)

## Files Created/Modified
- `src/convergence.jl` - Single builder build_preset_trajectory_observables with inlined Mz
- `src/gap_estimation.jl` - Removed CrossValidationResult, cross_validate_gap; updated internal call
- `src/QuantumFurnace.jl` - Updated exports (removed 6 old, added 1 new)
- `test/test_convergence.jl` - Deleted 4 testsets, renamed calls, fixed index assertions
- `test/test_gap_estimation.jl` - Deleted cross-validation testset, renamed builder calls
- `experiments/validate_gap_estimation.jl` - DELETED
- `experiments/run_gap_validation.jl` - DELETED
- `experiments/run_sweep.jl` - DELETED
- `experiments/eigenmode_decomposition.jl` - DELETED

## Decisions Made
- Inlined Mz construction directly into build_preset_trajectory_observables rather than keeping build_total_magnetization as private helper. Simpler code with no loss of functionality.
- Updated test Mz cross-checks to use inline reference construction instead of calling deleted function. Preserves correctness assertion.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
- `eigenmode_decomposition.jl` was not tracked by git (untracked file), required `rm` instead of `git rm`. Resolved by using regular file deletion.
- `run_sweep.jl` had local modifications, required `git rm -f` to force removal. Expected since it was listed in git status at session start.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Codebase is clean with single observable builder pattern
- Ready for Plan 02 (ARPACK vs eigen verification) and Plan 03 (eigenbasis overlap analysis)
- No blockers

## Self-Check: PASSED
- [x] src/convergence.jl exists and contains build_preset_trajectory_observables
- [x] src/gap_estimation.jl exists and calls build_preset_trajectory_observables
- [x] src/QuantumFurnace.jl exports build_preset_trajectory_observables
- [x] Commit 04b56eb exists
- [x] Commit 36b22cd exists
- [x] Zero references to deleted functions in src/ and test/
- [x] All 653 tests pass

---
*Phase: 25-spectral-gap-validation-overhaul*
*Completed: 2026-02-18*
