---
phase: 33-type-foundation
plan: 02
subsystem: types
tags: [julia, parametric-types, type-parameter-dispatch, config-migration, trait-function]

# Dependency graph
requires:
  - "33-01: Config{S,D,C,T} struct, simulation/construction singletons, with_coherent trait"
provides:
  - "All 6 dispatch-heavy source files migrated from old 4-type system to Config{S,D,C,T}"
  - "2-way KMS/GNS dispatch via Config type parameter instead of 4-way concrete type dispatch"
  - "with_coherent accessed via trait function everywhere (no field access)"
  - "Serialization/deserialization using Config with construction/sim singletons"
affects: [33-03, 33-04, 34, 35, 36]

# Tech tracking
tech-stack:
  added: []
  patterns: [type-parameter-dispatch, trait-based-property-access, singleton-based-serialization]

key-files:
  created: []
  modified:
    - src/energy_domain.jl
    - src/bohr_domain.jl
    - src/furnace_utensils.jl
    - src/coherent.jl
    - src/misc_tools.jl
    - src/results.jl

key-decisions:
  - "GNS coherent validation check removed from validate_config! -- type system now enforces via trait"
  - "with_coherent kept in serialized dict for backward compat but computed via trait on write"
  - "_reconstruct_config builds Config with explicit sim/construction singletons instead of 4 old constructors"
  - "DLL support added to serialization type tags (config_type field) for forward compatibility"

patterns-established:
  - "Config{<:Any, <:Any, KMS/GNS} for construction-specific dispatch (2 methods, not 4)"
  - "Config{Lindbladian}/Config{Thermalize} for simulation-specific dispatch"
  - "Config{<:Any, D} for domain-specific dispatch"
  - "config.construction isa GNS for runtime type checks (replaces isa Union{...} patterns)"
  - "with_coherent(config.construction) for all coherent-term boolean checks"

# Metrics
duration: 5min
completed: 2026-02-25
---

# Phase 33 Plan 02: Dispatch Migration Summary

**Migrated 6 source files from 4-way concrete type dispatch to 2-way Config type-parameter dispatch with trait-based with_coherent access**

## Performance

- **Duration:** 5 min
- **Started:** 2026-02-25T21:58:35Z
- **Completed:** 2026-02-25T22:04:00Z
- **Tasks:** 2
- **Files modified:** 6

## Accomplishments
- Replaced all 4-method dispatch sites (pick_transition, _pick_alpha) with 2-method Config{<:Any,<:Any,KMS/GNS} dispatch
- Migrated all function signatures across 6 files from old type names (LiouvConfig, ThermalizeConfig, AbstractConfig, etc.) to Config
- Replaced all config.with_coherent field accesses with with_coherent(config.construction) trait calls
- Updated serialization to use config.construction/config.sim for type detection instead of isa checks
- Removed GNS coherent validation check from validate_config! (type system enforces it)
- Updated _reconstruct_config to return Config with appropriate singletons instead of 4 old constructors

## Task Commits

Each task was committed atomically:

1. **Task 1: Migrate energy_domain.jl and bohr_domain.jl dispatch** - `bb72953` (feat)
2. **Task 2: Migrate furnace_utensils.jl, coherent.jl, misc_tools.jl, results.jl** - `204c6f9` (feat)

## Files Created/Modified
- `src/energy_domain.jl` - pick_transition 2-method dispatch, _pick_transition_kms/gns Config signatures, _truncate_energy_labels Config
- `src/bohr_domain.jl` - B_bohr (both), _pick_f, _pick_alpha 2-method dispatch, _pick_alpha_kms/gns Config signatures
- `src/furnace_utensils.jl` - _precompute_labels/data Config signatures, with_coherent trait call, _select_b_plus_calculator KMS dispatch
- `src/coherent.jl` - _precompute_coherent_total_B/unitary_terms/terms Config signatures, with_coherent trait calls
- `src/misc_tools.jl` - validate_config!, _collect_config_errors!, _generate_filename, _print_press all updated to Config
- `src/results.jl` - ExperimentResult{Config}, _config_to_dict, _reconstruct_config, _write_companion_txt, filename/dir functions

## Decisions Made
- GNS coherent check removed entirely from validate_config! rather than converted -- the type system now makes it impossible to construct a GNS config with with_coherent=true
- with_coherent kept in serialized Dict output for backward compatibility with existing BSON files, but computed via trait on write
- DLL added to type tag serialization (config_type field) for forward compatibility even though DLL construction is not yet implemented
- domain removed from _dict_to_config_kwargs since _reconstruct_config passes it explicitly to Config constructor

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- All 6 dispatch-heavy files now use Config exclusively
- Module still will NOT compile -- remaining files (liouvillian.jl, thermalize.jl, trajectory.jl, krylov.jl, etc.) still reference old types
- Plans 03-04 will complete the migration of remaining callers and test files

## Self-Check: PASSED

- FOUND: src/energy_domain.jl
- FOUND: src/bohr_domain.jl
- FOUND: src/furnace_utensils.jl
- FOUND: src/coherent.jl
- FOUND: src/misc_tools.jl
- FOUND: src/results.jl
- FOUND: commit bb72953 (Task 1)
- FOUND: commit 204c6f9 (Task 2)

---
*Phase: 33-type-foundation*
*Completed: 2026-02-25*
