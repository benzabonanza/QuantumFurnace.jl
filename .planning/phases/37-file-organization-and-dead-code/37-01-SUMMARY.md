---
phase: 37-file-organization-and-dead-code
plan: 01
subsystem: codebase-cleanup
tags: [dead-code, staging, project-toml, julia]

# Dependency graph
requires:
  - phase: 36-api-and-results
    provides: "New entry points (run_lindblad, run_thermalize, run_krylov_spectrum, run_trajectory) that replace old ones"
provides:
  - "Clean furnace.jl with only active entry points (construct_lindbladian, run_lindblad, run_thermalize)"
  - "Clean structs.jl with only active structs (no DMSimulationResult, LindbladianResult, LSIFramework)"
  - "Staging area src/staging/ and test/staging/ for dormant gap estimation code"
  - "Clean Project.toml without Distributed, LsqFit, Optim deps"
affects: [37-02-module-definition, phase-38]

# Tech tracking
tech-stack:
  added: []
  patterns: ["src/staging/ directory for dormant code excluded from module"]

key-files:
  created:
    - src/staging/gap_estimation.jl
    - src/staging/fitting.jl
    - src/staging/log_sobolev.jl
    - test/staging/test_gap_estimation.jl
    - test/staging/test_fitting.jl
  modified:
    - src/furnace.jl
    - src/structs.jl
    - src/convergence.jl
    - src/krylov_eigsolve.jl
    - src/QuantumFurnace.jl
    - Project.toml

key-decisions:
  - "Removed using Distributed/LsqFit/Optim and include() calls for moved/deleted files from QuantumFurnace.jl as blocking fix (Rule 3)"
  - "Left stale export lines (LindbladianResult, DMSimulationResult, run_lindbladian, etc.) for Plan 02 export reorganization"

patterns-established:
  - "staging/ directory: dormant code moved here is excluded from module includes and test suite"

# Metrics
duration: 4min
completed: 2026-02-27
---

# Phase 37 Plan 01: Dead Code and Staging Summary

**Removed old entry points (run_lindbladian, run_thermalization) and dead structs, moved gap estimation code to src/staging/, deleted placeholder files, cleaned Project.toml of stale deps**

## Performance

- **Duration:** 4 min
- **Started:** 2026-02-27T15:40:03Z
- **Completed:** 2026-02-27T15:44:26Z
- **Tasks:** 2
- **Files modified:** 10

## Accomplishments
- Deleted run_lindbladian() and run_thermalization() old entry points from furnace.jl (125 lines removed)
- Deleted DMSimulationResult, LindbladianResult, LSIFramework dead structs from structs.jl (61 lines removed)
- Moved gap_estimation.jl, fitting.jl, log_sobolev.jl to src/staging/ with git history preserved
- Moved test_gap_estimation.jl, test_fitting.jl to test/staging/
- Deleted empty placeholder files errors.jl and kraus.jl
- Removed Distributed, LsqFit, Optim from Project.toml [deps] and [compat]
- Cleaned 4 stale comments referencing old function names across convergence.jl, krylov_eigsolve.jl, furnace.jl, structs.jl

## Task Commits

Each task was committed atomically:

1. **Task 1: Delete dead code, clean stale comments** - `a30a29c` (feat)
2. **Task 2: Move staging files, delete placeholders, clean Project.toml** - `e89ce9a` (feat)

## Files Created/Modified
- `src/furnace.jl` - Removed run_lindbladian (36 lines) and run_thermalization (71 lines), cleaned stale comment
- `src/structs.jl` - Removed DMSimulationResult, LindbladianResult, LSIFramework structs, updated ThermalizeScratch docstring
- `src/convergence.jl` - Updated Gibbs helper docstring reference from run_thermalization to run_thermalize
- `src/krylov_eigsolve.jl` - Updated channel docstring reference from run_thermalization to run_thermalize
- `src/QuantumFurnace.jl` - Removed using Distributed/Optim/LsqFit and include() for deleted/moved files
- `Project.toml` - Removed Distributed, LsqFit, Optim from [deps] and [compat]
- `src/staging/gap_estimation.jl` - Moved from src/ (git mv)
- `src/staging/fitting.jl` - Moved from src/ (git mv)
- `src/staging/log_sobolev.jl` - Moved from src/ (git mv)
- `test/staging/test_gap_estimation.jl` - Moved from test/ (git mv)
- `test/staging/test_fitting.jl` - Moved from test/ (git mv)
- `src/errors.jl` - Deleted (empty placeholder)
- `src/kraus.jl` - Deleted (empty placeholder)

## Decisions Made
- Removed `using Distributed`, `using Optim`, `using LsqFit` and `include()` calls for moved/deleted files from QuantumFurnace.jl as part of Task 2. Without this, the module would fail to load after removing the deps from Project.toml and deleting/moving the files. This was a Rule 3 blocking fix.
- Left stale export lines (LindbladianResult, DMSimulationResult, run_lindbladian, run_thermalization, fit_exponential_decay, estimate_spectral_gap) in the export list for Plan 02 to handle as part of the comprehensive export reorganization.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Removed using/include statements for deleted/moved files from QuantumFurnace.jl**
- **Found during:** Task 2 (Move staging files, delete placeholders, clean Project.toml)
- **Issue:** Removing Distributed/LsqFit/Optim from Project.toml [deps] while leaving `using Distributed`, `using Optim`, `using LsqFit` in the module file, and leaving `include("errors.jl")`, `include("kraus.jl")`, `include("log_sobolev.jl")`, `include("fitting.jl")`, `include("gap_estimation.jl")` would cause the module to fail to load.
- **Fix:** Removed the 3 `using` statements and 5 `include()` calls from QuantumFurnace.jl.
- **Files modified:** src/QuantumFurnace.jl
- **Verification:** Module file now references only existing deps and files
- **Committed in:** e89ce9a (Task 2 commit)

**2. [Rule 3 - Blocking] Rebuilt corrupt git index**
- **Found during:** Task 1 commit
- **Issue:** Git index was corrupt (bad index file sha1 signature), preventing commits
- **Fix:** Removed corrupt index file and ran `git reset` to rebuild
- **Files modified:** .git/index (internal)
- **Verification:** Subsequent git operations succeed normally
- **Committed in:** N/A (git internal fix)

---

**Total deviations:** 2 auto-fixed (2 blocking)
**Impact on plan:** Both auto-fixes essential for module loadability and git operation. No scope creep.

## Issues Encountered
- Git index corruption required index rebuild before first commit could succeed. Fixed by removing corrupt index and running `git reset`.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Source files cleaned, staging area established
- Ready for Plan 02: Module definition and export list reorganization
- Export list still contains references to deleted types/functions (LindbladianResult, DMSimulationResult, run_lindbladian, etc.) -- Plan 02 will clean these

---
*Phase: 37-file-organization-and-dead-code*
*Completed: 2026-02-27*
