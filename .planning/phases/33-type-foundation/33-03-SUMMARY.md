---
phase: 33-type-foundation
plan: 03
subsystem: types
tags: [julia, parametric-types, type-parameter-dispatch, config-migration, pipeline-migration]

# Dependency graph
requires:
  - "33-01: Config{S,D,C,T} struct, simulation/construction singletons, with_coherent trait"
  - "33-02: Dispatch migration of 6 source files (energy_domain, bohr_domain, furnace_utensils, coherent, misc_tools, results)"
provides:
  - "All 8 simulation pipeline files migrated from old 4-type system to Config{S,D,C,T}"
  - "_thermalize_to_liouv_config deleted -- unified Config eliminates need for conversion"
  - "KrylovWorkspace accepts Config{Lindbladian} and Config{Thermalize} directly"
  - "apply_lindbladian!/apply_adjoint_lindbladian! dispatch on Config{Lindbladian, D}"
  - "krylov_spectral_gap dispatches on Config{Lindbladian} and Config{Thermalize}"
affects: [33-04, 34, 35, 36]

# Tech tracking
tech-stack:
  added: []
  patterns: [unified-config-passthrough, simulation-type-dispatch, domain-only-dispatch]

key-files:
  created: []
  modified:
    - src/furnace.jl
    - src/jump_workers.jl
    - src/trajectories.jl
    - src/convergence.jl
    - src/gap_estimation.jl
    - src/krylov_workspace.jl
    - src/krylov_matvec.jl
    - src/krylov_eigsolve.jl

key-decisions:
  - "construct_lindbladian constrained to Config{Lindbladian} (only called from run_lindbladian)"
  - "_accumulate_R_total! uses Config{<:Any, D} -- domain-only dispatch since both sim types need it"
  - "_accumulate_jump_sandwich! uses Config{<:Any, D} -- called from channel path with Config{Thermalize}"
  - "apply_lindbladian!/apply_adjoint_lindbladian! constrained to Config{Lindbladian, D} -- only used in Lindbladian path"
  - "apply_delta_channel! accepts Config (unconstrained) -- receives Config{Thermalize} from channel path"

patterns-established:
  - "Config{Lindbladian, D} for Lindbladian-specific functions (construct, matvec, eigsolve)"
  - "Config{Thermalize, D} for thermalization-specific functions (run_thermalization, trajectories, convergence)"
  - "Config{<:Any, D} for shared helpers that serve both paths (_accumulate_R_total!, _accumulate_jump_sandwich!)"
  - "with_coherent(config.construction) for all coherent-term boolean checks (no config.with_coherent field access)"

# Metrics
duration: 7min
completed: 2026-02-25
---

# Phase 33 Plan 03: Pipeline Migration Summary

**Migrated 8 simulation pipeline files to Config{S,D,C,T} dispatch and deleted _thermalize_to_liouv_config (~60 lines of dead conversion code)**

## Performance

- **Duration:** 7 min
- **Started:** 2026-02-25T22:05:59Z
- **Completed:** 2026-02-25T22:13:34Z
- **Tasks:** 2
- **Files modified:** 8

## Accomplishments
- Replaced all AbstractLiouvConfig, AbstractThermalizeConfig, and AbstractConfig references across 8 pipeline files
- Deleted both _thermalize_to_liouv_config methods (~60 lines) -- unified Config eliminates the need for Thermalize-to-Liouvillian conversion
- KrylovWorkspace Thermalize constructor now passes Config directly to _precompute_data and _precompute_coherent_total_B
- config.with_coherent replaced with with_coherent(config.construction) in all remaining call sites

## Task Commits

Each task was committed atomically:

1. **Task 1: Migrate furnace.jl, jump_workers.jl, trajectories.jl, convergence.jl, gap_estimation.jl** - `da976ac` (feat)
2. **Task 2: Migrate krylov_workspace.jl, krylov_matvec.jl, krylov_eigsolve.jl and delete _thermalize_to_liouv_config** - `301d7a5` (feat)

## Files Created/Modified
- `src/furnace.jl` - run_lindbladian Config{Lindbladian,D,C,Tc}, construct_lindbladian Config{Lindbladian}, run_thermalization Config{Thermalize,D,C,Tc}
- `src/jump_workers.jl` - All _jump_contribution! methods: Config{Lindbladian,D} for Liouvillian, Config{Thermalize,D} for thermalization
- `src/trajectories.jl` - TrajectoryFramework config field, build_trajectoryframework, _build_framework_and_seed, run_trajectories, run_observable_trajectories, _precompute_R all use Config{Thermalize}; config.with_coherent -> with_coherent(config.construction)
- `src/convergence.jl` - run_trajectories_convergence and run_trajectories_adaptive accept Config{Thermalize}
- `src/gap_estimation.jl` - estimate_spectral_gap accepts Config{Thermalize}
- `src/krylov_workspace.jl` - KrylovWorkspace(Config{Lindbladian}), KrylovWorkspace(Config{Thermalize}), _accumulate_R_total! uses Config{<:Any,D}, _thermalize_to_liouv_config deleted
- `src/krylov_matvec.jl` - apply_lindbladian! and apply_adjoint_lindbladian! dispatch on Config{Lindbladian,D}
- `src/krylov_eigsolve.jl` - krylov_spectral_gap Config{Lindbladian} and Config{Thermalize}, apply_delta_channel! accepts Config, _accumulate_jump_sandwich! uses Config{<:Any,D}

## Decisions Made
- construct_lindbladian constrained to Config{Lindbladian} since it is only called from run_lindbladian (diagnostics can relax this constraint if needed later)
- _accumulate_R_total! and _accumulate_jump_sandwich! use Config{<:Any, D} (domain-only constraint) because both Lindbladian and Thermalize paths call these helpers
- apply_delta_channel! accepts bare Config (no type parameter constraint) since it receives Config{Thermalize} from channel path and Config{Lindbladian} could work too

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- All 8 pipeline files now use Config{S,D,C,T} dispatch exclusively
- Combined with Plans 01-02, a total of 16 source files have been migrated
- Plan 04 will complete the migration of test files
- Module still will NOT compile -- test files and any remaining callers still reference old types

## Self-Check: PASSED

- FOUND: src/furnace.jl
- FOUND: src/jump_workers.jl
- FOUND: src/trajectories.jl
- FOUND: src/convergence.jl
- FOUND: src/gap_estimation.jl
- FOUND: src/krylov_workspace.jl
- FOUND: src/krylov_matvec.jl
- FOUND: src/krylov_eigsolve.jl
- FOUND: commit da976ac (Task 1)
- FOUND: commit 301d7a5 (Task 2)

---
*Phase: 33-type-foundation*
*Completed: 2026-02-25*
