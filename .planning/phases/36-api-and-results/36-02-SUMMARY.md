---
phase: 36-api-and-results
plan: 02
subsystem: api
tags: [julia, entry-points, lindblad, thermalize, krylov, spectral-gap, results]

# Dependency graph
requires:
  - phase: 36-01
    provides: "LindbladResults, ThermalizeResults, KrylovSpectrumResults structs and _capture_metadata"
  - phase: 35-workspace-consolidation
    provides: "Unified Workspace{S,D,C,T}, Config{S,D,C,T}, existing run_lindbladian/run_thermalization/krylov_spectral_gap"
provides:
  - "run_lindblad(jumps, config, hamiltonian, trotter=nothing) -> LindbladResults"
  - "run_thermalize(jumps, config, hamiltonian, trotter=nothing; initial_dm, rng, rescale_by_inv_prob) -> ThermalizeResults"
  - "run_krylov_spectrum(jumps, config, hamiltonian, trotter=nothing; krylov_kwargs...) -> KrylovSpectrumResults"
  - "Uniform positional signature (jumps, config, hamiltonian, trotter) across all 3 entry points"
  - "Wall time capture via _capture_metadata in all 3 entry points"
affects: [36-03-PLAN, 36-04-PLAN, 37-removal]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Uniform positional signature (jumps, config, hamiltonian, trotter) for all run_* entry points"
    - "run_* wraps existing internal logic and returns typed Result struct with config + metadata"
    - "Old functions untouched -- removal deferred to Phase 37"

key-files:
  created: []
  modified:
    - "src/furnace.jl"
    - "src/krylov_eigsolve.jl"
    - "src/QuantumFurnace.jl"

key-decisions:
  - "run_thermalize defaults initial_dm to maximally mixed I/d (keyword arg) rather than requiring it positionally"
  - "run_krylov_spectrum uses Union{Lindbladian, Thermalize} type constraint to support both config types"
  - "run_krylov_spectrum reorders args internally to match krylov_spectral_gap(config, hamiltonian, jumps) signature"

patterns-established:
  - "Thin wrapper pattern: new run_* functions inline existing logic or delegate, adding timing + Result struct wrapping"
  - "Metadata capture at entry point: t_start = time() at top, _capture_metadata(wall_time_seconds=...) before return"

# Metrics
duration: 4min
completed: 2026-02-27
---

# Phase 36 Plan 02: Entry Points Summary

**Three new public entry points (run_lindblad, run_thermalize, run_krylov_spectrum) with uniform positional signature returning typed Result structs with timing metadata**

## Performance

- **Duration:** 4 min
- **Started:** 2026-02-27T12:15:43Z
- **Completed:** 2026-02-27T12:19:44Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments
- Added run_lindblad wrapping existing Arpack eigs logic, returning LindbladResults with sorted eigenvalues, fixed point, gap mode, spectral gap, and timing metadata
- Added run_thermalize with keyword initial_dm defaulting to maximally mixed I/d, returning ThermalizeResults with final DM, trace distances, time steps, and timing metadata
- Added run_krylov_spectrum wrapping krylov_spectral_gap with uniform argument order, returning KrylovSpectrumResults for both Config{Lindbladian} and Config{Thermalize}
- All three use uniform positional signature: run_*(jumps, config, hamiltonian, trotter=nothing)
- Existing run_lindbladian, run_thermalization, and krylov_spectral_gap remain untouched
- All 1198 existing tests pass with zero regressions

## Task Commits

Each task was committed atomically:

1. **Task 1: Create run_lindblad and run_thermalize entry points in furnace.jl** - `b7262e9` (feat)
2. **Task 2: Create run_krylov_spectrum entry point in krylov_eigsolve.jl** - `8d0d4eb` (feat)

## Files Created/Modified
- `src/furnace.jl` - Added run_lindblad (returns LindbladResults) and run_thermalize (returns ThermalizeResults) after existing functions
- `src/krylov_eigsolve.jl` - Added run_krylov_spectrum (returns KrylovSpectrumResults) after existing krylov_spectral_gap functions
- `src/QuantumFurnace.jl` - Added exports for run_lindblad, run_thermalize, run_krylov_spectrum

## Decisions Made
- run_thermalize accepts initial_dm as keyword argument defaulting to maximally mixed state I/d (research recommendation: keeps positional signature uniform across all 4 entry points)
- run_krylov_spectrum constrains S<:Union{Lindbladian, Thermalize} to support both config types via a single function
- run_krylov_spectrum internally reorders arguments to match krylov_spectral_gap's (config, hamiltonian, jumps; ...) signature while presenting the uniform (jumps, config, hamiltonian, trotter) external API

## Deviations from Plan

None - plan executed exactly as written. Both tasks were already partially implemented from a prior session; this execution verified correctness, ran the full test suite, and committed atomically.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- 3 of 4 entry points complete (run_lindblad, run_thermalize, run_krylov_spectrum)
- Plan 03 (run_trajectory) can proceed -- it wraps existing trajectory infrastructure
- Plan 04 (comprehensive tests) can test all entry points and save/load round-trips
- All old API functions remain functional for backward compatibility

## Self-Check: PASSED

All files verified present, all commits verified in git log.

---
*Phase: 36-api-and-results*
*Completed: 2026-02-27*
