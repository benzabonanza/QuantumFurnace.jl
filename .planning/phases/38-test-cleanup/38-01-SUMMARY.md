---
phase: 38-test-cleanup
plan: 01
subsystem: testing
tags: [julia, test-infrastructure, factory-pattern, config-consolidation]

# Dependency graph
requires:
  - phase: 37-nufft-allocation
    provides: "Stabilized test suite with all 1057 tests passing"
provides:
  - "Unified make_config(sim, domain; kwargs...) factory replacing 8 old factories"
  - "Unified make_test_system(; num_qubits=N) factory replacing 2 old factories"
  - "N3_* globals replacing SMALL_* across all test files"
  - "ALL_DOMAINS constant for parametric test loops"
  - "trajectory_validation integrated into runtests.jl behind QUANTUMFURNACE_FULL_TESTS gate"
affects: [38-02, 38-03, 38-04, 38-05]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "make_config(sim, domain; num_qubits, construction, delta, mixing_time) for all Config construction in tests"
    - "N3_ prefix for 3-qubit test globals (N3_HAM, N3_JUMPS, N3_GIBBS, N3_DIM, N3_TROTTER, N3_TROTTER_JUMPS)"
    - "ALL_DOMAINS = [EnergyDomain(), TimeDomain(), TrotterDomain(), BohrDomain()] for parametric sweeps"

key-files:
  created: []
  modified:
    - test/test_helpers.jl
    - test/runtests.jl
    - test/test_allocation.jl
    - test/test_compilation.jl
    - test/test_convergence.jl
    - test/test_cptp.jl
    - test/test_diagnostics.jl
    - test/test_dm_detailed_balance.jl
    - test/test_dm_scaling.jl
    - test/test_gns_trajectory.jl
    - test/test_krylov_crossvalidation.jl
    - test/test_krylov_eigsolve.jl
    - test/test_krylov_matvec.jl
    - test/test_observable_trajectories.jl
    - test/test_regression.jl
    - test/test_results.jl
    - test/test_threading.jl
    - test/test_trajectory_fixes.jl
    - test/test_workspace_independence.jl
    - test/staging/test_gap_estimation.jl
    - test/trajectory_validation/run_trajectory_validation.jl
    - test/trajectory_validation/run_convergence_tests.jl
    - test/reference/generate_references.jl

key-decisions:
  - "Unified make_config uses Config(; ...) with leading semicolon for keyword-only NamedTuple splatting"
  - "Default construction=KMS() in make_config matches most common usage; GNS callers pass explicitly"
  - "Deleted local n=6 factories in test_krylov_crossvalidation.jl; unified make_config handles all num_qubits values"
  - "Trajectory validation gated at runtests.jl level to save compilation time during normal testing"

patterns-established:
  - "make_config(sim, domain; kwargs...): single parametrized factory for all test Config construction"
  - "make_test_system(; num_qubits, trotter): single parametrized factory for test system creation"
  - "N3_ prefix: consistent naming for 3-qubit precomputed test fixtures"
  - "QUANTUMFURNACE_FULL_TESTS env gate: slow test gating pattern in runtests.jl"

# Metrics
duration: 18min
completed: 2026-02-28
---

# Phase 38 Plan 01: Test Infrastructure Consolidation Summary

**Unified 8 config factories into single make_config(sim, domain; kwargs...), renamed SMALL_* to N3_*, deleted 7 dead test files, integrated trajectory_validation into runtests.jl**

## Performance

- **Duration:** ~18 min
- **Started:** 2026-02-28T08:50:00Z
- **Completed:** 2026-02-28T09:08:00Z
- **Tasks:** 3
- **Files modified:** 30

## Accomplishments
- Consolidated 8 config factory functions (make_liouv_config, make_liouv_config_gns, make_thermalize_config, make_small_liouv_config, make_small_liouv_config_gns, make_small_thermalize_config, make_small_thermalize_config_gns, plus duplicate) into single `make_config(sim, domain; num_qubits, construction, delta, mixing_time)`
- Consolidated 2 system factory functions (make_test_system, make_small_test_system) into single `make_test_system(; num_qubits, trotter)`
- Renamed all SMALL_* globals to N3_* across 22 test files (net -625 lines of code)
- Deleted 7 dead test files in test/old_tests/ (-516 lines)
- Integrated trajectory_validation into runtests.jl behind QUANTUMFURNACE_FULL_TESTS env gate
- All 1057 tests pass with zero regressions

## Task Commits

Each task was committed atomically:

1. **Task 1: Consolidate test_helpers.jl factories and globals** - `b066e89` (feat)
2. **Task 2: Migrate all test files to new factories and N3_* globals** - `a58434c` (refactor)
3. **Task 3: Delete old_tests, integrate trajectory_validation, update runtests.jl** - `f298283` (chore)

## Files Created/Modified
- `test/test_helpers.jl` - Unified make_config, make_test_system, N3_* globals, ALL_DOMAINS constant
- `test/runtests.jl` - Added trajectory_validation integration behind env gate
- `test/test_allocation.jl` - Migrated to make_config API and N3_* globals
- `test/test_compilation.jl` - Migrated to make_config API
- `test/test_convergence.jl` - Migrated to make_config API
- `test/test_cptp.jl` - Migrated to make_config API
- `test/test_diagnostics.jl` - Migrated to make_config API and N3_* globals
- `test/test_dm_detailed_balance.jl` - Replaced inline Config() with make_config, N3_* globals
- `test/test_dm_scaling.jl` - Migrated to make_config API
- `test/test_gns_trajectory.jl` - Migrated to make_config API and N3_* globals
- `test/test_krylov_crossvalidation.jl` - Deleted 3 local n=6 factories, migrated to unified API
- `test/test_krylov_eigsolve.jl` - Migrated to make_config API
- `test/test_krylov_matvec.jl` - Migrated to make_config API with explicit construction=GNS()
- `test/test_observable_trajectories.jl` - Migrated to make_config API and N3_* globals
- `test/test_regression.jl` - Migrated to make_config API and N3_* globals
- `test/test_results.jl` - Migrated to make_config API and N3_* globals
- `test/test_threading.jl` - Migrated to N3_* globals
- `test/test_trajectory_fixes.jl` - Migrated to make_config API
- `test/test_workspace_independence.jl` - Migrated to N3_* globals
- `test/staging/test_gap_estimation.jl` - Migrated to make_config API and N3_* globals
- `test/trajectory_validation/run_trajectory_validation.jl` - Migrated to make_config API and N3_* globals
- `test/trajectory_validation/run_convergence_tests.jl` - Migrated to make_config API and N3_* globals
- `test/reference/generate_references.jl` - Migrated to make_config API and N3_* globals
- `test/old_tests/` (7 files deleted) - B_test.jl, ham_test.jl, kossakowski_test.jl, log_sobolev_test.jl, time_tests.jl, trajectory_test.jl, trott_test.jl

## Decisions Made
- Used `Config(; ...)` with leading semicolon for keyword-only NamedTuple splatting (Julia @kwdef struct requires this syntax when splatting NamedTuples)
- Default `construction=KMS()` in make_config matches the majority of test callers; GNS-specific tests pass `construction=GNS()` explicitly
- Deleted local n=6 factory functions in test_krylov_crossvalidation.jl rather than keeping them; the unified make_config(; num_qubits=6) handles this cleanly
- Gated trajectory_validation at the runtests.jl level (not just within the individual files) to avoid compilation overhead during normal test runs

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Config keyword splatting syntax**
- **Found during:** Task 1 (test_helpers.jl consolidation)
- **Issue:** Using `Config(... therm_kw...)` caused NamedTuple values to be interpreted as positional Float64 arguments, producing MethodError
- **Fix:** Changed to `Config(; ... therm_kw...)` with leading semicolon to force keyword-only splatting
- **Files modified:** test/test_helpers.jl
- **Verification:** `make_config(Thermalize(), EnergyDomain())` produces valid Config
- **Committed in:** b066e89 (Task 1 commit)

**2. [Rule 1 - Bug] replace_all corrupted function definition in test_krylov_crossvalidation.jl**
- **Found during:** Task 2 (test file migration)
- **Issue:** Bulk replace of `make_n6_liouv_config(` also modified the function definition line `function make_n6_liouv_config(domain; ...)` into invalid Julia syntax
- **Fix:** Deleted all 3 local n=6 factory function definitions (now redundant) and fixed the replacement at call sites
- **Files modified:** test/test_krylov_crossvalidation.jl
- **Verification:** All 1057 tests pass
- **Committed in:** a58434c (Task 2 commit)

**3. [Rule 1 - Bug] Double semicolons in test_krylov_crossvalidation.jl**
- **Found during:** Task 2 (test file migration)
- **Issue:** replace_all of `make_n6_test_system(` to `make_test_system(; num_qubits=6,` created `make_test_system(; num_qubits=6,; trotter=...)` with invalid double semicolons
- **Fix:** Manually fixed to single semicolon syntax
- **Files modified:** test/test_krylov_crossvalidation.jl
- **Verification:** All 1057 tests pass
- **Committed in:** a58434c (Task 2 commit)

---

**Total deviations:** 3 auto-fixed (3 bugs)
**Impact on plan:** All auto-fixes were necessary for correct Julia syntax. No scope creep.

## Issues Encountered
None beyond the auto-fixed deviations above.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Unified make_config factory ready for parametric threshold work in plans 02-05
- ALL_DOMAINS constant enables clean domain sweeps
- N3_* naming consistent and searchable
- trajectory_validation properly integrated for full test runs

## Self-Check: PASSED

All 23 modified/created files verified present on disk. All 3 task commits (b066e89, a58434c, f298283) verified in git log. test/old_tests/ directory confirmed deleted.

---
*Phase: 38-test-cleanup*
*Completed: 2026-02-28*
