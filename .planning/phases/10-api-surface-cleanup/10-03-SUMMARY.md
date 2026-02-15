---
phase: 10-api-surface-cleanup
plan: 03
subsystem: api
tags: [internal-functions, underscore-prefix, cross-file-calls, julia-module]

# Dependency graph
requires:
  - phase: 10-api-surface-cleanup
    plan: 02
    provides: "All ~45 internal function definitions prefixed with _ across 14 source files"
provides:
  - "All cross-file call sites updated to _-prefixed function names"
  - "All test qualified access updated to _-prefixed function names"
  - "Phase 10 API surface cleanup fully complete"
affects: []

# Tech tracking
tech-stack:
  added: []
  patterns: ["Cross-file calls consistently use _ prefix convention"]

key-files:
  created: []
  modified:
    - src/furnace.jl
    - src/jump_workers.jl
    - src/trajectories.jl
    - src/furnace_utensils.jl
    - src/coherent.jl
    - test/test_dm_scaling.jl
    - test/test_compilation.jl
    - test/test_cptp.jl
    - test/test_regression.jl
    - test/test_trajectory_fixes.jl
    - test/trajectory_validation/run_trajectory_validation.jl
    - test/trajectory_validation/run_convergence_tests.jl

key-decisions:
  - "Updated docstrings in coherent.jl to match _-prefixed function names"
  - "Fixed trajectory_validation scripts from 3-arg to 2-arg _precompute_data call"

patterns-established:
  - "All non-exported functions use _ prefix in both definitions and call sites"
  - "Test qualified access uses QuantumFurnace._func_name for internal functions"

# Metrics
duration: 8min
completed: 2026-02-15
---

# Phase 10 Plan 03: Cross-file Call Site Update Summary

**Updated all cross-file call sites and test qualified access to _-prefixed internal function names, completing Phase 10 API surface cleanup**

## Performance

- **Duration:** 8 min
- **Started:** 2026-02-15T13:29:34Z
- **Completed:** 2026-02-15T13:38:12Z
- **Tasks:** 2
- **Files modified:** 12

## Accomplishments
- All cross-file call sites in 5 source files updated to _-prefixed names (~30 call sites)
- All test qualified access patterns updated across 7 test files
- All 224 tests pass including Aqua.jl quality checks
- Phase 10 API surface cleanup is now complete: export list curated (10-01), definitions renamed (10-02), cross-file calls updated (10-03)

## Task Commits

Each task was committed atomically:

1. **Task 1: Update cross-file call sites in source files** - `6b8ee6e` (feat)
2. **Task 2: Update test files for _-prefixed internal function names** - `f5b29a4` (feat)

## Files Created/Modified
- `src/furnace.jl` - Updated _precompute_data, _print_press, _precompute_coherent_total_B, _vectorize_liouvillian_coherent!, _jump_contribution!, _precompute_coherent_unitary_terms
- `src/jump_workers.jl` - Updated _vectorize_liouv_diss_and_add! (7 sites), _vectorize_liouvillian_coherent! (3 sites), _prefactor_view (3 sites)
- `src/trajectories.jl` - Updated _precompute_data, _print_press, _precompute_coherent_total_B, _prefactor_view (4 sites)
- `src/furnace_utensils.jl` - Updated _create_energy_labels, _truncate_energy_labels, _pick_alpha, _truncate_time_labels_for_oft, _prepare_oft_nufft_prefactors, _compute_truncated_func, _compute_b_minus, _compute_b_plus, _compute_b_plus_metro, _compute_b_plus_smooth
- `src/coherent.jl` - Updated docstrings to match _-prefixed function names
- `test/test_dm_scaling.jl` - Updated _precompute_data, _create_energy_labels, _truncate_time_labels_for_oft, _time_oft!, _trotter_oft!, _prefactor_view
- `test/test_compilation.jl` - Updated _precompute_data
- `test/test_cptp.jl` - Updated _precompute_data
- `test/test_regression.jl` - Updated _precompute_data
- `test/test_trajectory_fixes.jl` - Updated _precompute_data
- `test/trajectory_validation/run_trajectory_validation.jl` - Updated _precompute_data (fixed 3-arg -> 2-arg)
- `test/trajectory_validation/run_convergence_tests.jl` - Updated _precompute_data (fixed 3-arg -> 2-arg)

## Decisions Made
- Updated docstrings in coherent.jl to match the _-prefixed function names (consistency)
- Fixed trajectory_validation standalone scripts that used an old 3-arg `precompute_data(domain, config, ham)` call to the current 2-arg `_precompute_data(config, ham)` signature

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed trajectory_validation scripts using obsolete 3-arg precompute_data**
- **Found during:** Task 2 (test file updates)
- **Issue:** `run_trajectory_validation.jl` and `run_convergence_tests.jl` used `precompute_data(domain, config, ham_or_trott)` -- an obsolete 3-arg form that no longer exists
- **Fix:** Updated to 2-arg `_precompute_data(config, ham_or_trott)` since config already contains the domain
- **Files modified:** test/trajectory_validation/run_trajectory_validation.jl, test/trajectory_validation/run_convergence_tests.jl
- **Verification:** These are standalone scripts not in the test suite; function signature matches current API
- **Committed in:** f5b29a4 (Task 2 commit)

**2. [Rule 2 - Missing Critical] Updated additional test files not mentioned in plan**
- **Found during:** Task 2 (test file updates)
- **Issue:** Plan only listed test_dm_scaling.jl, but test_compilation.jl, test_cptp.jl, test_regression.jl, and test_trajectory_fixes.jl also used `QuantumFurnace.precompute_data` qualified access
- **Fix:** Updated all test files to use `QuantumFurnace._precompute_data`
- **Files modified:** test/test_compilation.jl, test_cptp.jl, test_regression.jl, test_trajectory_fixes.jl
- **Verification:** All 224 tests pass
- **Committed in:** f5b29a4 (Task 2 commit)

**3. [Rule 2 - Missing Critical] Updated coherent.jl docstrings**
- **Found during:** Task 1 (source file updates)
- **Issue:** Docstrings for _precompute_coherent_total_B and _precompute_coherent_unitary_terms still showed old un-prefixed names
- **Fix:** Updated docstring function name headers to match _-prefixed definitions
- **Files modified:** src/coherent.jl
- **Verification:** Module loads without errors
- **Committed in:** 6b8ee6e (Task 1 commit)

---

**Total deviations:** 3 auto-fixed (1 bug fix, 2 missing critical)
**Impact on plan:** All auto-fixes necessary for correctness. No scope creep.

## Issues Encountered
- Plan referenced `construct_disordering_terms` in `trotter_domain.jl` at line ~243, but this function only exists in `hamiltonian.jl` (intra-file, already handled in Plan 10-02). No action needed.
- `pick_transition` calls were correctly left unchanged as it is an exported function (confirmed from Plan 10-02 decision).

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Phase 10 (API Surface Cleanup) is fully complete
- All exported functions accessible via `using QuantumFurnace` without qualification
- All internal functions use _ prefix convention and require `QuantumFurnace._func_name` for access
- All 224 tests pass with Aqua.jl quality checks
- Ready for Phase 11 or any subsequent work

## Self-Check: PASSED

- SUMMARY.md exists
- Commit 6b8ee6e exists (Task 1)
- Commit f5b29a4 exists (Task 2)
- All 12 modified files exist on disk
- Module loads without errors
- All 224 tests pass including Aqua.jl quality checks

---
*Phase: 10-api-surface-cleanup*
*Completed: 2026-02-15*
