---
phase: 08-struct-simplification
plan: 03
subsystem: structs
tags: [julia, type-parameters, dispatch, domain-dispatch, trajectory-framework]

# Dependency graph
requires:
  - phase: 08-01
    provides: Config structs with Union{..., Nothing} fields and domain type parameter
provides:
  - Simplified TrajectoryFramework{T,D} with 2 type parameters
  - Domain dispatch via config type parameter instead of separate domain argument
  - Clean function signatures without redundant domain positional args
affects: [09-type-parameterization]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Domain dispatch via config type parameter: f(config::AbstractConfig{D}) where {D<:...}"
    - "TrajectoryFramework{T,D} with only essential type parameters for dispatch"

key-files:
  created: []
  modified:
    - src/trajectories.jl
    - src/furnace_utensils.jl
    - src/jump_workers.jl
    - src/misc_tools.jl
    - src/furnace.jl

key-decisions:
  - "Reduce TrajectoryFramework from {T,C,H,PD,D} to {T,D} -- C,H,PD not used for dispatch"
  - "Domain dispatch uses config type parameter: f(config::AbstractConfig{D}) instead of f(::DomainType, config)"
  - "Config domain::D field retained for runtime isa checks and display purposes"
  - "Collapsed duplicate Time/Trotter NUFFT prefactor code block in precompute_data"

patterns-established:
  - "Domain dispatch pattern: f(config::AbstractConfig{D}) where {D<:Union{TimeDomain, TrotterDomain}}"
  - "Config type parameter carries domain info -- no need to pass domain as separate arg"

# Metrics
duration: 12min
completed: 2026-02-15
---

# Phase 8 Plan 3: TrajectoryFramework Simplification and Domain Dispatch Refactor Summary

**TrajectoryFramework reduced from 5 to 2 type params {T,D}, domain dispatch refactored to use config type parameter across 13 dispatch functions**

## Performance

- **Duration:** 12 min
- **Started:** 2026-02-15T09:17:31Z
- **Completed:** 2026-02-15T09:29:52Z
- **Tasks:** 3
- **Files modified:** 13

## Accomplishments
- Reduced TrajectoryFramework type parameters from {T,C,H,PD,D} to {T,D}, eliminating 3 unused dispatch parameters
- Refactored 13 domain-dispatched function definitions to use config type parameter instead of separate domain argument
- Updated all call sites (source and test files) to remove explicit domain arguments
- Collapsed duplicate TimeDomain/TrotterDomain NUFFT prefactor code block in precompute_data
- All 224 tests pass with zero regressions

## Task Commits

Each task was committed atomically:

1. **Task 1: Simplify TrajectoryFramework type parameters** - `83f89b5` (refactor)
2. **Task 2a: Refactor core dispatch function signatures** - `a3678f8` (refactor)
3. **Task 2b: Update all call sites to remove explicit domain arguments** - `93d5c46` (refactor)

## Files Created/Modified
- `src/trajectories.jl` - TrajectoryFramework{T,D} struct, step_along_trajectory! signatures, precompute_R definitions, build/run call sites
- `src/furnace_utensils.jl` - precompute_labels and precompute_data dispatch on config type param
- `src/jump_workers.jl` - jump_contribution! (6 variants) dispatch on config type param
- `src/misc_tools.jl` - _collect_config_errors! (4 variants) dispatch on config type param, validate_config! call site
- `src/furnace.jl` - construct_lindbladian and run_thermalization call sites
- `test/test_compilation.jl` - Updated precompute_data call sites
- `test/test_cptp.jl` - Updated precompute_data call sites
- `test/test_dm_scaling.jl` - Updated precompute_data call sites
- `test/test_regression.jl` - Updated precompute_data call sites
- `test/test_trajectory_fixes.jl` - Updated precompute_data call sites
- `test/old_tests/{B_test,time_tests,trajectory_test}.jl` - Updated legacy test call sites

## Decisions Made
- Reduced TrajectoryFramework from {T,C,H,PD,D} to {T,D}: only T (element type) and D (domain) are used for dispatch. C, H, PD are accessed at runtime only.
- Config struct's domain::D field retained for runtime isa checks (e.g., `config.domain isa TrotterDomain` in furnace.jl) and display purposes (string(typeof(config.domain)) in misc_tools.jl).
- Collapsed the duplicate if/elseif TimeDomain/TrotterDomain code block in precompute_data's last variant since both branches had identical NUFFT prefactor computation.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Updated test files to match new function signatures**
- **Found during:** Task 2b (update call sites)
- **Issue:** Test files (test_dm_scaling.jl, test_regression.jl, test_compilation.jl, test_cptp.jl, test_trajectory_fixes.jl) still used old `precompute_data(config.domain, config, ...)` pattern
- **Fix:** Updated all test call sites to use new `precompute_data(config, ...)` pattern
- **Files modified:** 8 test files
- **Verification:** All 224 tests pass
- **Committed in:** 93d5c46 (Task 2b commit)

---

**Total deviations:** 1 auto-fixed (1 blocking)
**Impact on plan:** Test file updates were necessary for correctness -- plan focused on source files but test files also had call sites. No scope creep.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Phase 08 (struct simplification) is now complete with all 3 plans executed
- TrajectoryFramework has clean {T,D} parameterization ready for Phase 9 (type parameterization)
- Domain dispatch pattern via config type param established for all domain-specific functions
- All 224 tests pass, codebase ready for next phase

## Self-Check: PASSED

All source files exist. All 3 commit hashes verified. All 224 tests pass.

---
*Phase: 08-struct-simplification*
*Completed: 2026-02-15*
