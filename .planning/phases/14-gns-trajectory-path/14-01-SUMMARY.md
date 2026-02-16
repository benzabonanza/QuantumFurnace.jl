---
phase: 14-gns-trajectory-path
plan: 01
subsystem: testing
tags: [gns, lindbladian, trajectory, density-matrix, detailed-balance, cptp]

# Dependency graph
requires:
  - phase: 13-multi-threaded-engine
    provides: "TrajectoryFramework, run_trajectories, step_along_trajectory! with workspace/RNG separation"
provides:
  - "GNS config factory functions (make_small_liouv_config_gns, make_small_thermalize_config_gns)"
  - "GNS trajectory validation test suite (test_gns_trajectory.jl)"
  - "GNS-to-Gibbs approximation gap baseline (EnergyDomain: 0.081, BohrDomain: 0.035)"
  - "Fix for LiouvConfigGNS/ThermalizeConfigGNS @kwdef construction bug"
affects: [18-kms-vs-gns-comparison]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "GNS config factory pattern matching KMS factories with with_coherent=false"
    - "Outer constructor pattern for @kwdef structs with custom inner constructors"

key-files:
  created:
    - test/test_gns_trajectory.jl
  modified:
    - test/test_helpers.jl
    - test/runtests.jl
    - src/structs.jl

key-decisions:
  - "GNS-to-Gibbs approximation gap at sigma=0.1: EnergyDomain=0.081, BohrDomain=0.035"
  - "Trajectory convergence at ntraj=1000, delta=0.01, mixing_time=5.0: trace distance 0.029 to GNS fixed point"
  - "Fixed @kwdef outer constructor bug for GNS config structs (pre-existing, never triggered before)"

patterns-established:
  - "GNS test factories: make_small_liouv_config_gns(domain) and make_small_thermalize_config_gns(domain; delta, mixing_time)"
  - "GNS trajectory convergence pattern: compute Lindbladian fixed point as reference, compare trajectory DM to that (not Gibbs)"

# Metrics
duration: 5min
completed: 2026-02-16
---

# Phase 14 Plan 01: GNS Trajectory Path Summary

**GNS trajectory validation with 4 testsets covering Lindbladian fixed point, CPTP completeness, trajectory convergence (0.029 trace distance), and BohrDomain detailed balance, plus GNS-to-Gibbs gap baseline for Phase 18**

## Performance

- **Duration:** 5 min
- **Started:** 2026-02-16T06:57:51Z
- **Completed:** 2026-02-16T07:02:47Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments
- GNS Lindbladian fixed point extracted and validated as proper density matrix (Hermitian, trace 1, PSD) for both EnergyDomain and BohrDomain
- Trajectory-averaged density matrix converges to GNS fixed point with trace distance 0.029 (under 0.05 threshold) at ntraj=1000, delta=0.01, mixing_time=5.0
- GNS-to-Gibbs approximation gap documented: EnergyDomain=0.081, BohrDomain=0.035 (Phase 18 baseline)
- Fixed pre-existing bug where LiouvConfigGNS and ThermalizeConfigGNS were unconstructable via @kwdef

## Task Commits

Each task was committed atomically:

1. **Task 1: Add GNS config factory functions** - `dfe36ef` (feat) - also fixed GNS struct outer constructors
2. **Task 2: Create GNS trajectory validation test suite** - `e02c456` (test) - 4 testsets, all 284 tests pass

## Files Created/Modified
- `test/test_gns_trajectory.jl` - GNS trajectory validation: Lindbladian fixed point, CPTP, convergence, BohrDomain
- `test/test_helpers.jl` - Added make_small_liouv_config_gns and make_small_thermalize_config_gns factories
- `test/runtests.jl` - Include test_gns_trajectory.jl in test suite
- `src/structs.jl` - Added outer constructors for LiouvConfigGNS and ThermalizeConfigGNS

## Decisions Made
- GNS-to-Gibbs gap at sigma=0.1, beta=10.0, n=3: EnergyDomain=0.081, BohrDomain=0.035. These serve as Phase 18 baselines for the two-sigma comparison.
- Trajectory convergence parameters: ntraj=1000, delta=0.01 (delta_eff=0.09), mixing_time=5.0 (500 steps). Achieved 0.029 trace distance to GNS fixed point, well under the 0.05 threshold.
- Used `result.rho_mean` (actual TrajectoryResult field name) instead of `result.density_matrix` (plan's incorrect reference).

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed LiouvConfigGNS and ThermalizeConfigGNS unconstructable via @kwdef**
- **Found during:** Task 1 (GNS config factory functions)
- **Issue:** The @kwdef macro generates an outer constructor that calls the struct positionally, but GNS structs only had inner constructors (LiouvConfigGNS{D,T}(...)), not outer constructors (LiouvConfigGNS(...)). The @kwdef-generated code could not infer type parameters, making GNS configs impossible to construct via keyword arguments.
- **Fix:** Added outer constructor functions for both LiouvConfigGNS and ThermalizeConfigGNS that infer D from domain and T from beta, then forward to the inner constructor with validation.
- **Files modified:** src/structs.jl
- **Verification:** Both factory functions construct valid configs; all 284 tests pass
- **Committed in:** dfe36ef (Task 1 commit)

**2. [Rule 1 - Bug] Corrected result.density_matrix to result.rho_mean**
- **Found during:** Task 2 (test suite creation)
- **Issue:** Plan referenced `result.density_matrix` but the actual TrajectoryResult struct field is `rho_mean`
- **Fix:** Used correct field name `result.rho_mean` in test assertions
- **Files modified:** test/test_gns_trajectory.jl
- **Verification:** Tests access the correct field and pass
- **Committed in:** e02c456 (Task 2 commit)

---

**Total deviations:** 2 auto-fixed (2 bugs)
**Impact on plan:** Both fixes necessary for correctness. The @kwdef bug was pre-existing and would have blocked any GNS config usage. No scope creep.

## Issues Encountered
None - all tests passed on first run with the planned parameters (ntraj=1000, delta=0.01, mixing_time=5.0).

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- GNS trajectory path fully validated with documented approximation gaps
- Phase 18 (KMS vs GNS comparison) has baseline metrics: EnergyDomain gap=0.081, BohrDomain gap=0.035
- All 284 tests pass including new GNS tests, zero regressions
- No blockers or concerns

---
*Phase: 14-gns-trajectory-path*
*Completed: 2026-02-16*
