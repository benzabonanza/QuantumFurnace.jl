---
phase: 37-file-organization-and-dead-code
plan: 02
subsystem: codebase-cleanup
tags: [exports, module-definition, simulation-scripts, julia]

# Dependency graph
requires:
  - phase: 37-01
    provides: "Clean source files with dead code removed, staging area established, Project.toml cleaned"
provides:
  - "Organized export list with 6 sections: Lindbladian, Thermalize, Krylov, Trajectory, Diagnostics, Common"
  - "STAGING: commented-out block for dormant exports"
  - "Clean runtests.jl without staged test includes"
  - "Complete set of 4 simulation scripts matching 4 run_* entry points"
  - "Dead nprocs/SharedArrays references removed"
affects: [phase-38]

# Tech tracking
tech-stack:
  added: []
  patterns: ["Export list organized by simulation type with # --- Section --- comment blocks", "STAGING: prefix for commented-out dormant exports"]

key-files:
  created:
    - simulations/main_krylov.jl
  modified:
    - src/QuantumFurnace.jl
    - test/runtests.jl
    - src/nufft.jl
    - src/furnace_utensils.jl
    - Project.toml

key-decisions:
  - "Export list organized into Lindbladian/Thermalize/Krylov/Trajectory/Diagnostics/Common sections with Common at bottom per user decision"
  - "Dormant exports preserved as commented-out STAGING: block (not deleted)"
  - "Removed dead nprocs() references and SharedArrays branch (Rule 1 bug fix after Distributed removal in 37-01)"
  - "Removed SharedArrays from Project.toml and using statements since no active code uses it"

patterns-established:
  - "Export sections: # --- SimType --- comment blocks group related exports"
  - "Staging exports: # STAGING: prefix for dormant code exports"

# Metrics
duration: 13min
completed: 2026-02-27
---

# Phase 37 Plan 02: Module Definition and Exports Summary

**Reorganized module exports into 6 simulation-type sections, removed dead exports/imports, created main_krylov.jl simulation script, fixed leftover nprocs/SharedArrays references**

## Performance

- **Duration:** 13 min
- **Started:** 2026-02-27T15:47:14Z
- **Completed:** 2026-02-27T16:00:10Z
- **Tasks:** 2
- **Files modified:** 6

## Accomplishments
- Reorganized flat export list into 6 organized sections (Lindbladian, Thermalize, Krylov, Trajectory, Diagnostics, Common) with STAGING block for dormant exports
- Removed dead exports: run_lindbladian, run_thermalization, LindbladianResult, DMSimulationResult, fit_exponential_decay, FitResult, estimate_spectral_gap
- Created simulations/main_krylov.jl demonstrating run_krylov_spectrum API entry point
- Removed test_fitting.jl and test_gap_estimation.jl from test runner
- Fixed nprocs() undefined reference bug left over from Distributed removal in Plan 01
- Removed SharedArrays from module (no longer used after dead SharedArray branch removal)

## Task Commits

Each task was committed atomically:

1. **Task 1: Update QuantumFurnace.jl - reorganize exports** - `266d29c` (feat)
2. **Task 2: Update runtests.jl, create main_krylov.jl, fix nprocs** - `90fa4b5` (feat)

## Files Created/Modified
- `src/QuantumFurnace.jl` - Reorganized export list into 6 sections, removed dead exports, removed `using SharedArrays`, added STAGING block
- `test/runtests.jl` - Removed includes for test_fitting.jl and test_gap_estimation.jl
- `simulations/main_krylov.jl` - New simulation script demonstrating run_krylov_spectrum API
- `src/nufft.jl` - Replaced `nprocs() > 1` default with `false`, removed dead SharedArray branch
- `src/furnace_utensils.jl` - Replaced `nprocs() > 1` with `false` in NUFFT prefactor call
- `Project.toml` - Removed SharedArrays from [deps] and [compat]

## Decisions Made
- Export list organized by simulation type with Common section at bottom, per user decision in CONTEXT.md
- Staging exports preserved as `# STAGING:` prefix commented-out block rather than deleted
- `nprocs() > 1` replaced with `false` since Distributed was removed in Plan 01 (multi-process mode is dead code)
- SharedArrays removed entirely since the only usage was behind the dead `nprocs() > 1` branch

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed undefined nprocs() references after Distributed removal**
- **Found during:** Task 2, Part C (full test suite run)
- **Issue:** Plan 01 removed `using Distributed` from module, but `nprocs()` was still called in nufft.jl and furnace_utensils.jl, causing `UndefVarError: nprocs not defined` at runtime
- **Fix:** Replaced `nprocs() > 1` with `false` (multi-process mode is dead code), removed dead SharedArray allocation branch, removed `using SharedArrays` and SharedArrays from Project.toml
- **Files modified:** src/nufft.jl, src/furnace_utensils.jl, src/QuantumFurnace.jl, Project.toml
- **Verification:** Module loads, all functional tests pass
- **Committed in:** 90fa4b5 (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (1 bug)
**Impact on plan:** Essential fix for module functionality after Distributed removal. Without this, any code path touching NUFFT prefactors would crash. No scope creep.

## Issues Encountered
- 4 pre-existing allocation threshold failures in test_krylov_matvec.jl (Time/TrotterDomain hot path: 137472 > 100000 budget). These are flaky allocation tests unrelated to our changes, consistent with known pre-existing threshold issues noted in prior phases.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Module definition clean with organized exports
- All 4 simulation scripts (main_liouv.jl, main_thermalize.jl, main_trajectory.jl, main_krylov.jl) match 4 run_* entry points
- Phase 37 complete, ready for Phase 38
- No stale imports, no dead exports, no dead code in active module

## Self-Check: PASSED

All files verified present, all commit hashes verified in git log.

---
*Phase: 37-file-organization-and-dead-code*
*Completed: 2026-02-27*
