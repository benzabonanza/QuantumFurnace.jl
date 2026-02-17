---
phase: 20-observable-infrastructure
plan: 01
subsystem: simulation
tags: [quantum-observables, magnetization, spectral-gap, basis-transform, convergence]

# Dependency graph
requires:
  - phase: 16-convergence-tracking
    provides: "build_convergence_observables pattern, ConvergenceData struct, Gibbs helpers"
provides:
  - "build_total_magnetization(ham, n; trotter) for per-site M_z in eigenbasis or Trotter basis"
  - "build_gap_estimation_observables(ham, n; trotter) for H + M_z bundle"
  - "Regression tests (testsets 19-22) with Gibbs trace verification in both bases"
affects: [22-trajectory-runner, 23-gap-estimation-api, 24-cross-validation]

# Tech tracking
tech-stack:
  added: []
  patterns: ["unified trotter keyword argument instead of separate _trotter suffix functions"]

key-files:
  created: []
  modified:
    - src/convergence.jl
    - src/QuantumFurnace.jl
    - test/test_convergence.jl

key-decisions:
  - "M_z = sum(Z_i)/n per-site normalization (locked decision)"
  - "Single function with trotter keyword argument, not separate _trotter suffix"
  - "build_gap_estimation_observables internally calls build_total_magnetization (DRY)"
  - "No ZZ correlations in gap estimation bundle (deferred decision)"
  - "Trotter H test uses explicit transform verification instead of off-diagonal magnitude check"

patterns-established:
  - "Unified trotter keyword: build_func(ham, n; trotter=nothing) pattern for dual-basis support"
  - "Gap estimation observable bundle: H first, then M_z in return tuple"

# Metrics
duration: 8min
completed: 2026-02-17
---

# Phase 20 Plan 01: Observable Infrastructure Summary

**M_z = sum(Z_i)/n and H + M_z gap estimation bundle in both Hamiltonian eigenbasis and Trotter eigenbasis, with Gibbs trace regression tests**

## Performance

- **Duration:** 8 min
- **Started:** 2026-02-17T08:10:27Z
- **Completed:** 2026-02-17T08:19:06Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments
- `build_total_magnetization` constructs per-site M_z in computational basis, transforms to selected eigenbasis
- `build_gap_estimation_observables` bundles H + M_z for spectral gap estimation, correctly handling H diagonality in eigenbasis vs full transform in Trotter basis
- 4 new testsets (19-22) cover both functions in both bases with analytical Gibbs trace verification and cross-basis consistency checks
- All 578 tests pass (including 22 new assertions from the 4 testsets)

## Task Commits

Each task was committed atomically:

1. **Task 1: Implement build_total_magnetization and build_gap_estimation_observables** - `d3b57bc` (feat)
2. **Task 2: Add regression tests for magnetization and gap estimation observables** - `424ccb8` (test)

## Files Created/Modified
- `src/convergence.jl` - Added build_total_magnetization and build_gap_estimation_observables functions with docstrings
- `src/QuantumFurnace.jl` - Added exports for both new functions
- `test/test_convergence.jl` - Added testsets 19-22 covering both bases with Gibbs trace verification

## Decisions Made
- Used `vcat([H], mz_obs)` for proper `Vector{Matrix}` concatenation (Julia `[H; mz_obs]` fails with dimension mismatch for matrix-in-vector types)
- Trotter basis H test verifies `H == V_T' * H_data * V_T` rather than asserting off-diagonal magnitude, because the test system's Trotter eigenbasis is very close to the Hamiltonian eigenbasis (Trotter error ~1e-9)

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed Vector concatenation in build_gap_estimation_observables**
- **Found during:** Task 2 (test execution revealed error in Task 1 code)
- **Issue:** `[H; mz_obs]` attempted vertical matrix concatenation instead of Vector{Matrix} concatenation, causing DimensionMismatch
- **Fix:** Changed to `vcat([H], mz_obs)` and `vcat(["H"], mz_names)` for proper 1D vector concatenation
- **Files modified:** src/convergence.jl
- **Verification:** All tests pass after fix
- **Committed in:** d3b57bc (amended Task 1 commit)

**2. [Rule 1 - Bug] Fixed Trotter H diagonality test assertion**
- **Found during:** Task 2 (testset 22 execution)
- **Issue:** Plan assumed H would have large off-diagonal elements in Trotter basis, but the test system's Trotter approximation is very accurate (off-diagonal ~3e-9 << 1e-6 threshold)
- **Fix:** Changed assertion from `maximum(abs.(H - diagm(diag(H)))) > 1e-6` to verifying H matches the explicit basis transform `V_T' * H_data * V_T`
- **Files modified:** test/test_convergence.jl
- **Verification:** Test now passes and correctly validates the transform
- **Committed in:** 424ccb8 (Task 2 commit)

---

**Total deviations:** 2 auto-fixed (2 bug fixes)
**Impact on plan:** Both fixes necessary for correctness. No scope creep.

## Issues Encountered
None beyond the auto-fixed deviations above.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Observable infrastructure complete for Phase 22 (trajectory runner) and Phase 23 (gap estimation API)
- `build_total_magnetization` and `build_gap_estimation_observables` are exported and tested in both eigenbasis and Trotter basis
- Return format matches existing `(Vector{Matrix{ComplexF64}}, Vector{String})` convention expected by `run_trajectories_convergence` and `run_trajectories_adaptive`

---
*Phase: 20-observable-infrastructure*
*Completed: 2026-02-17*
