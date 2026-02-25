---
phase: 33-type-foundation
plan: 04
subsystem: types
tags: [julia, parametric-types, config, migration, test-suite]

# Dependency graph
requires:
  - phase: 33-03
    provides: "Config{S,D,C,T} dispatch migration in all src/ pipeline files"
provides:
  - "Complete codebase migration -- zero old type references in any .jl or .ipynb file"
  - "Full test suite validation confirming Config{S,D,C,T} works end-to-end"
  - "Phase 33 success criteria met: unified type system operational"
affects: [34-file-restructure, 35-api-cleanup]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Config factory functions use construction=KMS()/GNS() kwarg instead of with_coherent Bool"
    - "isa checks use Config{Lindbladian}, Config{Thermalize}, Config{S,<:Any,GNS} patterns"

key-files:
  created: []
  modified:
    - "test/test_helpers.jl (factory functions: make_liouv_config, make_thermalize_config, etc.)"
    - "test/test_compilation.jl (isa checks)"
    - "test/test_results.jl (isa checks with construction type parameter)"
    - "test/test_krylov_crossvalidation.jl (n6 factory functions)"
    - "simulations/main_thermalize.jl, main_liouv.jl, main_krylov_benchmark.jl"
    - "experiments/*.jl (7 experiment scripts)"
    - "docs/src/literate/tutorial_thermalize.jl"
    - "src/krylov_eigsolve.jl, src/krylov_workspace.jl (comment-only)"

key-decisions:
  - "Playground files migrated on disk but not committed (gitignored)"
  - "Pre-existing test_diagnostics.jl failures accepted (numerical threshold, not migration-related)"

patterns-established:
  - "Config(sim=Lindbladian(), domain=..., construction=KMS(), ...) replaces LiouvConfig(...)"
  - "Config(sim=Thermalize(), domain=..., construction=KMS(), ...) replaces ThermalizeConfig(...)"
  - "construction=GNS() replaces with_coherent=false in factory kwargs"
  - "isa Config{Lindbladian, <:Any, GNS} replaces isa LiouvConfigGNS"

# Metrics
duration: 23min
completed: 2026-02-25
---

# Phase 33 Plan 04: Consumer Migration Summary

**Complete codebase migration of 37+ files from 4-type config system to unified Config{S,D,C,T}, validated by full test suite (1187 tests pass)**

## Performance

- **Duration:** 23 min
- **Started:** 2026-02-25T22:16:06Z
- **Completed:** 2026-02-25T22:39:09Z
- **Tasks:** 2
- **Files modified:** 37 (23 test + 14 non-test tracked files, plus 5 gitignored playground files on disk)

## Accomplishments

- Migrated all 23 test files: factory functions return Config with sim/construction singletons, isa checks use Config{Lindbladian}/Config{Thermalize}
- Migrated 3 simulation scripts, 7 experiment scripts, 2 doc files, and 2 src comment references
- Migrated 5 gitignored playground files on disk (including Jupyter notebook)
- Verified zero references to old type names (LiouvConfig, ThermalizeConfig, etc.) remain in any .jl or .ipynb file
- Full test suite passes: 1187 passed (7 pre-existing failures in test_diagnostics.jl unrelated to migration)

## Task Commits

Each task was committed atomically:

1. **Task 1: Migrate test files** - `5cb331c` (feat)
2. **Task 2: Migrate non-src/non-test files and validate** - `4114811` (feat)

## Files Created/Modified

### Task 1: Test files (23 files)
- `test/test_helpers.jl` - Factory functions: make_liouv_config, make_thermalize_config, make_small_* variants
- `test/test_compilation.jl` - isa LiouvConfig -> isa Config{Lindbladian}, isa ThermalizeConfig -> isa Config{Thermalize}
- `test/test_results.jl` - Round-trip isa checks with construction type: Config{Lindbladian,<:Any,KMS}, Config{Thermalize,<:Any,GNS}
- `test/test_krylov_crossvalidation.jl` - make_n6_liouv_config and make_n6_thermalize_config factories
- `test/test_krylov_eigsolve.jl` - Testset name references
- `test/test_dm_detailed_balance.jl` - Inline LiouvConfig -> Config(sim=Lindbladian(), ...)
- `test/test_gns_trajectory.jl` - References to GNS config types
- `test/old_tests/trajectory_test.jl` - Old config constructors
- `test/old_tests/B_test.jl` - Old config constructors
- `test/old_tests/time_tests.jl` - Old config constructors
- `test/old_tests/kossakowski_test.jl` - Old config constructors
- `test/test_krylov_matvec.jl` - Config type references
- `test/test_trajectory_fixes.jl` - Config type references
- `test/test_observable_trajectories.jl` - Config type references
- `test/test_gap_estimation.jl` - Config type references
- `test/test_workspace_independence.jl` - Config type references
- `test/test_regression.jl` - Config type references
- `test/test_threading.jl` - Config type references
- `test/test_allocation.jl` - Config type references
- `test/test_convergence.jl` - Config type references
- `test/trajectory_validation/run_convergence_tests.jl` - Config type references
- `test/trajectory_validation/run_trajectory_validation.jl` - Config type references
- `test/reference/generate_references.jl` - Config type references

### Task 2: Non-test files (14 tracked + 5 gitignored)
- `simulations/main_thermalize.jl` - ThermalizeConfig -> Config(sim=Thermalize(), construction=KMS())
- `simulations/main_liouv.jl` - LiouvConfig -> Config(sim=Lindbladian(), construction=KMS())
- `simulations/main_krylov_benchmark.jl` - Factory function updated with construction kwarg
- `experiments/validate_spectral_gap.jl` - Both factory functions migrated
- `experiments/diagnose_gap_momentum.jl` - Factory function migrated
- `experiments/validate_gap_delta_scaling.jl` - Both factory functions migrated
- `experiments/validate_gap_delta_scaling_v2.jl` - Both factory functions migrated
- `experiments/validate_gap_xx_stagg.jl` - Both factory functions migrated
- `experiments/validate_gap_disordered.jl` - Both factory functions migrated
- `experiments/investigate_delta_scaling_bug.jl` - Both factory functions migrated
- `docs/src/literate/tutorial_thermalize.jl` - ThermalizeConfig -> Config(sim=Thermalize())
- `docs/src/generated/tutorial_thermalize.ipynb` - Generated notebook updated
- `src/krylov_eigsolve.jl` - Comment: ThermalizeConfig -> Config{Thermalize}
- `src/krylov_workspace.jl` - Two comments: ThermalizeConfig -> Config{Thermalize}
- `playground/vegyessali.jl` (gitignored) - Config constructors updated
- `playground/unraveling.jl` (gitignored) - Vector{LiouvConfig} -> Vector{Config}
- `playground/src_playground.jl` (gitignored) - config.with_coherent -> with_coherent(config.construction)
- `playground/trotter_match.jl` (gitignored) - Config constructor updated
- `playground/jl_playground.ipynb` (gitignored) - Notebook cell updated

## Decisions Made

- **Playground files not committed:** All 5 playground files are in .gitignore. Migrations applied on disk for developer convenience but cannot be version-controlled.
- **Pre-existing test failures accepted:** test_diagnostics.jl has 7 failures (numerical threshold: 0.035 < 0.01 assertion). Verified pre-existing by running test suite on stashed (pre-migration) code -- same 7 failures. Not caused by type migration.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 2 - Missing Critical] Migrated docs tutorial files not listed in plan**
- **Found during:** Task 2 (comprehensive grep sweep)
- **Issue:** `docs/src/literate/tutorial_thermalize.jl` and `docs/src/generated/tutorial_thermalize.ipynb` contained ThermalizeConfig references but were not listed in the plan's files_modified
- **Fix:** Migrated both files using same patterns as other files
- **Files modified:** docs/src/literate/tutorial_thermalize.jl, docs/src/generated/tutorial_thermalize.ipynb
- **Verification:** grep confirms zero old type references in docs/
- **Committed in:** 4114811 (Task 2 commit)

**2. [Rule 2 - Missing Critical] Updated src/ comment references not listed in plan**
- **Found during:** Task 2 (comprehensive grep sweep)
- **Issue:** src/krylov_eigsolve.jl and src/krylov_workspace.jl had ThermalizeConfig in code comments/docstrings
- **Fix:** Updated comment references to Config{Thermalize}
- **Files modified:** src/krylov_eigsolve.jl, src/krylov_workspace.jl
- **Verification:** grep confirms zero old type references in src/
- **Committed in:** 4114811 (Task 2 commit)

---

**Total deviations:** 2 auto-fixed (2 missing critical -- files not in plan's explicit list)
**Impact on plan:** Both fixes necessary for completeness. Without them, old type references would remain in docs and src comments, violating the "zero references" success criterion.

## Issues Encountered

- **Gitignored directories:** experiments/ and playground/ are in .gitignore but files are tracked (force-added previously). Required `git add -f` to stage experiment files. Playground files could not be committed at all since they were never tracked.
- **Jupyter notebook editing:** .ipynb files required python3 JSON manipulation scripts for reliable editing, since they are JSON files with embedded Julia code.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Phase 33 (Type Foundation) is COMPLETE. All 4 plans executed successfully.
- The unified `Config{S,D,C,T}` parametric type system is fully operational across the entire codebase.
- Ready for Phase 34 (File Restructure) which will reorganize src/ files for clarity.
- Adding new construction types (e.g., `DLL`) requires only `struct DLL <: AbstractConstruction end` plus dispatch methods.

## Self-Check: PASSED

- All key files exist on disk
- Both task commits found (5cb331c, 4114811)
- Zero old type name references in tracked .jl files

---
*Phase: 33-type-foundation*
*Completed: 2026-02-25*
