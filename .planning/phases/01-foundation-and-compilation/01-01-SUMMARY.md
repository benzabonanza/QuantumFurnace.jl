---
phase: 01-foundation-and-compilation
plan: 01
subsystem: testing, compilation
tags: [julia, trajectory, kraus, project-toml, test-infrastructure, fixtures]

# Dependency graph
requires: []
provides:
  - "Compilable trajectory code (build_trajectoryframework works with/without coherent term)"
  - "Clean Project.toml with dev deps in [extras] only"
  - "Test infrastructure: fixtures (TEST_HAM, TEST_JUMPS, TEST_GIBBS), tolerances, config factories"
  - "Working Pkg.test() with 22 passing tests"
affects: [02-bug-fixes-and-cptp, 03-dm-tests-and-aqua, 04-trajectory-validation, 05-regression-suite]

# Tech tracking
tech-stack:
  added: [StableRNGs, HypothesisTests, StatsBase, Aqua]
  patterns: [test-helpers-pattern, config-factory-pattern, tolerance-tiers]

key-files:
  created:
    - test/runtests.jl
    - test/test_helpers.jl
    - test/test_compilation.jl
  modified:
    - Project.toml
    - Manifest.toml
    - src/trajectories.jl
    - src/QuantumFurnace.jl
    - src/furnace.jl

key-decisions:
  - "Relaxed stdlib compat bounds to support Julia 1.11+ (original were Julia 1.12-only)"
  - "Load Hamiltonian via @__DIR__ path in tests instead of load_hamiltonian() to work under Pkg.test() temp dir"
  - "TOL_DELTA constant C=5.0 per Claude's discretion for unraveling error tolerance"
  - "TEST_DELTA = 0.01 for compilation smoke tests"

patterns-established:
  - "test-helpers: All shared test state in test/test_helpers.jl, included before testsets"
  - "config-factories: make_liouv_config() and make_thermalize_config() with locked physical params"
  - "tolerance-tiers: TOL_EXACT (1e-12), TOL_QUADRATURE (1e-6), TOL_DELTA(delta) = 5*delta"

# Metrics
duration: 8min
completed: 2026-02-13
---

# Phase 1 Plan 1: Foundation and Compilation Summary

**Fixed 5 trajectory compilation bugs, separated dev/test deps in Project.toml, and built test infrastructure with 22 passing smoke tests via Pkg.test()**

## Performance

- **Duration:** 8 min
- **Started:** 2026-02-13T17:32:26Z
- **Completed:** 2026-02-13T17:41:02Z
- **Tasks:** 2
- **Files modified:** 8

## Accomplishments
- Fixed all 5 compilation bugs blocking trajectory code (dangling `where T`, invalid `trotter` kwarg, uninitialized `B_total`, broken export block, furnace.jl kwarg)
- Cleaned Project.toml: moved 6 dev-only deps to [extras], added 4 test deps (StableRNGs, HypothesisTests, StatsBase, Aqua)
- Created complete test infrastructure with shared fixtures, tolerance tiers, and config factory functions
- All 22 smoke tests pass: module loading, build_trajectoryframework with/without coherent, fixture validity, tolerance definitions, config factories

## Task Commits

Each task was committed atomically:

1. **Task 1: Fix compilation bugs and clean Project.toml** - `e138250` (fix)
2. **Task 2: Create test infrastructure and compilation smoke tests** - `2e2abb1` (feat)

## Files Created/Modified
- `src/trajectories.jl` - Fixed 3 bugs: dangling where T, invalid trotter kwarg, missing B_total=nothing
- `src/QuantumFurnace.jl` - Fixed export block (missing comma, removed dead export), removed ClusterManagers using
- `src/furnace.jl` - Fixed invalid trotter kwarg in precompute_coherent_total_B call
- `Project.toml` - Dev deps to [extras], test deps added, compat bounds relaxed for Julia 1.11
- `Manifest.toml` - Regenerated for clean dependency resolution
- `test/runtests.jl` - Pkg.test() entry point
- `test/test_helpers.jl` - Shared fixtures (TEST_HAM/JUMPS/GIBBS), tolerances, config factories
- `test/test_compilation.jl` - 22 smoke tests for compilation and fixture availability

## Decisions Made
- Relaxed stdlib compat bounds (LinearAlgebra, Pkg, SparseArrays, etc.) from Julia-1.12-only to "1.11, 1.12" to support Julia 1.11 runtime
- Used `@__DIR__` path resolution for Hamiltonian loading in tests instead of `load_hamiltonian()` which relies on `Pkg.project().path` (breaks under Pkg.test() temp directory)
- Set TOL_DELTA constant C=5.0 and TEST_DELTA=0.01 for smoke test parameters

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Regenerated Manifest.toml for clean dependency resolution**
- **Found during:** Task 1 (Project.toml cleanup)
- **Issue:** Existing Manifest.toml referenced unregistered packages (IrrationalConstants), stale from previous Julia version
- **Fix:** Deleted Manifest.toml and ran Pkg.instantiate() to regenerate
- **Files modified:** Manifest.toml
- **Verification:** Pkg.resolve() and Pkg.instantiate() succeed
- **Committed in:** e138250 (Task 1 commit)

**2. [Rule 3 - Blocking] Relaxed stdlib compat bounds for Julia 1.11**
- **Found during:** Task 1 (Project.toml cleanup)
- **Issue:** Compat entries for SparseArrays, LinearAlgebra, Pkg were pinned to 1.12.0 which is incompatible with Julia 1.11 runtime
- **Fix:** Changed to "1.11, 1.12" ranges; relaxed third-party deps to major version ranges
- **Files modified:** Project.toml
- **Verification:** Pkg.instantiate() succeeds on Julia 1.11.3
- **Committed in:** e138250 (Task 1 commit)

**3. [Rule 3 - Blocking] Fixed Hamiltonian load path for Pkg.test() environment**
- **Found during:** Task 2 (test infrastructure creation)
- **Issue:** load_hamiltonian() uses Pkg.project().path which points to a temp directory during Pkg.test(), causing FileNotFoundError for BSON file
- **Fix:** Used `dirname(@__DIR__)` to resolve source tree path and load BSON directly
- **Files modified:** test/test_helpers.jl
- **Verification:** Pkg.test() passes with Hamiltonian loaded correctly
- **Committed in:** 2e2abb1 (Task 2 commit)

**4. [Rule 1 - Bug] Removed AbstractDomain type annotation from test factory functions**
- **Found during:** Task 2 (test infrastructure creation)
- **Issue:** AbstractDomain is not exported from QuantumFurnace, causing UndefVarError in test
- **Fix:** Removed type annotation from domain parameter (duck typing works fine)
- **Files modified:** test/test_helpers.jl
- **Verification:** Pkg.test() passes
- **Committed in:** 2e2abb1 (Task 2 commit)

---

**Total deviations:** 4 auto-fixed (1 bug, 3 blocking)
**Impact on plan:** All auto-fixes necessary for correct execution in the runtime environment. No scope creep.

## Issues Encountered
None beyond the auto-fixed deviations above.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Trajectory compilation works end-to-end (with and without coherent term)
- Test infrastructure ready for all subsequent phases (fixtures, tolerances, config factories)
- Pkg.test() pipeline established -- future phases just add include() calls to runtests.jl
- Phase 2 (Bug Fixes and CPTP) can proceed immediately

---
*Phase: 01-foundation-and-compilation*
*Completed: 2026-02-13*
