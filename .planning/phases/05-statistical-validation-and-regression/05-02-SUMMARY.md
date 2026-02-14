---
phase: 05-statistical-validation-and-regression
plan: 02
subsystem: testing
tags: [bson, regression, frozen-reference, trajectory, density-matrix]

# Dependency graph
requires:
  - phase: 04-trajectory-cross-validation
    provides: "Verified trajectory and DM evolution correctness for all domains"
  - phase: 03-dm-reference-test-suite
    provides: "DM test infrastructure, test_helpers.jl with SMALL system fixtures"
provides:
  - "Frozen BSON reference data for EnergyDomain and TrotterDomain (DM + trajectory)"
  - "Always-on regression tests (TINF-02) comparing fresh computation against frozen references"
  - "Generator script for regenerating references if intentional algorithm changes are made"
affects: []

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "BSON with symbol keys for Julia-idiomatic access (d[:rho] not d[\"rho\"])"
    - "Frozen reference pattern: generate once, commit binary, compare on every test run"

key-files:
  created:
    - test/reference/generate_references.jl
    - test/reference/energy_dm_reference.bson
    - test/reference/energy_traj_reference.bson
    - test/reference/trotter_coherent_dm_reference.bson
    - test/reference/trotter_coherent_traj_reference.bson
    - test/test_regression.jl
  modified:
    - test/runtests.jl

key-decisions:
  - "Symbol keys in BSON (not string keys) for idiomatic Julia d[:rho] access pattern"
  - "Tolerance 1e-10 for regression comparison (allows floating-point accumulation across platforms)"
  - "Trajectory seed=12345, ntraj=1000 for deterministic reference (distinct from Phase 4 seeds)"

patterns-established:
  - "Frozen BSON reference: save Matrix(rho) not Hermitian(rho) to avoid BSON round-trip issues"
  - "Regression test reads seed/ntraj from BSON metadata, not hardcoded in test"

# Metrics
duration: 5min
completed: 2026-02-14
---

# Phase 5 Plan 2: Regression Tests Summary

**Frozen BSON reference data with always-on regression tests verifying DM and trajectory numerical stability at 1e-10 tolerance (TINF-02)**

## Performance

- **Duration:** 5 min
- **Started:** 2026-02-14T16:26:12Z
- **Completed:** 2026-02-14T16:31:00Z
- **Tasks:** 2
- **Files modified:** 7

## Accomplishments
- 4 frozen BSON reference files generated for EnergyDomain DM/trajectory and TrotterDomain+coherent DM/trajectory
- Generator script for one-time reference data creation (reproducible via `julia --project test/reference/generate_references.jl`)
- 4 always-on regression tests in `test/test_regression.jl` verifying fresh computation matches frozen references
- All regression tests run in Pkg.test() in under 6 seconds
- TINF-02 requirement fully satisfied

## Task Commits

Each task was committed atomically:

1. **Task 1: Create reference data generator and generate frozen BSON files** - `2fee822` (feat)
2. **Task 2: Create always-on regression test and add to runtests.jl** - `7fac61e` (feat)

## Files Created/Modified
- `test/reference/generate_references.jl` - One-time generator script for frozen BSON reference data
- `test/reference/energy_dm_reference.bson` - Frozen DM reference for EnergyDomain 3-qubit system
- `test/reference/energy_traj_reference.bson` - Frozen trajectory reference for EnergyDomain 3-qubit system
- `test/reference/trotter_coherent_dm_reference.bson` - Frozen DM reference for TrotterDomain+coherent 3-qubit system
- `test/reference/trotter_coherent_traj_reference.bson` - Frozen trajectory reference for TrotterDomain+coherent 3-qubit system
- `test/test_regression.jl` - Always-on regression tests (TINF-02) with 4 sub-tests
- `test/runtests.jl` - Updated to include test_regression.jl

## Decisions Made
- Used symbol keys (`:rho`, `:delta`, etc.) in BSON dictionaries for idiomatic Julia access
- Tolerance 1e-10 for regression comparison (per user decision, allows floating-point accumulation)
- Seed 12345 with 1000 trajectories for trajectory references (distinct from Phase 4 seeds 42 and 123)
- Save `Matrix(rho)` not `Hermitian(rho)` to ensure BSON round-trip fidelity

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed BSON key type: string keys to symbol keys**
- **Found during:** Task 1 (generate BSON references)
- **Issue:** Using `Dict("rho" => ...)` with string keys caused `d[:rho]` access to fail with KeyError
- **Fix:** Changed to `Dict(:rho => ...)` with symbol keys for Julia-standard BSON access
- **Files modified:** test/reference/generate_references.jl
- **Verification:** BSON round-trip validated with `d[:rho]` access
- **Committed in:** 2fee822 (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (1 bug)
**Impact on plan:** Essential fix for BSON interoperability between generator and regression tests. No scope creep.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Phase 5 complete: both statistical validation (plan 1) and regression tests (plan 2) are done
- All 224 tests pass in Pkg.test() under 48 seconds
- Regression tests protect against future numerical drift

## Self-Check: PASSED

All 8 files verified present. Both task commits (2fee822, 7fac61e) confirmed in git history.

---
*Phase: 05-statistical-validation-and-regression*
*Completed: 2026-02-14*
