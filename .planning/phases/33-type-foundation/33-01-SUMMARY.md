---
phase: 33-type-foundation
plan: 01
subsystem: types
tags: [julia, parametric-types, singleton-dispatch, type-hierarchy, trait-function]

# Dependency graph
requires: []
provides:
  - "Config{S,D,C,T} unified parametric struct replacing 4 old config types"
  - "AbstractSimulation hierarchy with Lindbladian, Thermalize, KrylovSpectrum, Trajectory singletons"
  - "AbstractConstruction hierarchy with KMS, GNS, DLL singletons"
  - "with_coherent trait function dispatching on construction type"
affects: [33-02, 33-03, 33-04, 34, 35, 36, 37, 38]

# Tech tracking
tech-stack:
  added: []
  patterns: [singleton-type-dispatch, trait-function-for-boolean-property, unified-parametric-config]

key-files:
  created: []
  modified:
    - src/structs.jl
    - src/QuantumFurnace.jl

key-decisions:
  - "Field ordering: singletons first, then system params, physics params, grid params, thermalize-specific last"
  - "with_coherent removed as field, replaced by trait function dispatching on construction singleton"
  - "Outer constructor infers S,D,C from singleton args and T from beta kwarg"
  - "Optional fields (grid params, mixing_time, delta) default to nothing for cross-simulation compatibility"

patterns-established:
  - "Singleton dispatch: define abstract type + concrete singletons, use as both field and type parameter"
  - "Trait function: with_coherent(::KMS) = true pattern for compile-time boolean properties"
  - "Config dispatch: Config{Lindbladian,...} for simulation-specific, Config{<:Any,<:Any,KMS} for construction-specific"

# Metrics
duration: 2min
completed: 2026-02-25
---

# Phase 33 Plan 01: Type Foundation Summary

**Unified Config{S,D,C,T} struct with simulation/construction singleton hierarchies and with_coherent trait, replacing 4 duplicate config types**

## Performance

- **Duration:** 2 min
- **Started:** 2026-02-25T21:54:26Z
- **Completed:** 2026-02-25T21:56:20Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments
- Replaced 4 duplicate config structs (LiouvConfig, LiouvConfigGNS, ThermalizeConfig, ThermalizeConfigGNS) and 3 abstract supertypes with single Config{S,D,C,T}
- Defined AbstractSimulation hierarchy with 4 singletons: Lindbladian, Thermalize, KrylovSpectrum, Trajectory
- Defined AbstractConstruction hierarchy with 3 singletons: KMS, GNS, DLL
- Added with_coherent trait function: KMS->true, GNS->false, DLL->true (placeholder)
- Updated module exports to reference all new types and remove all old type names

## Task Commits

Each task was committed atomically:

1. **Task 1: Define type hierarchies, Config struct, and with_coherent trait** - `27b6669` (feat)
2. **Task 2: Update module exports** - `4e6a8c5` (feat)

## Files Created/Modified
- `src/structs.jl` - New type hierarchies, Config{S,D,C,T} struct, with_coherent trait; removed 4 old config structs and 3 abstract supertypes (net -62 lines)
- `src/QuantumFurnace.jl` - Updated export block: Config, AbstractSimulation, Lindbladian, Thermalize, KrylovSpectrum, Trajectory, AbstractConstruction, KMS, GNS, DLL, with_coherent

## Decisions Made
- Field ordering follows semantic grouping: singletons -> system params -> physics params -> grid params -> thermalize-specific
- Outer constructor captures S,D,C from singleton field values and T from beta kwarg to handle @kwdef 4-parameter inference
- Optional fields (grid params, mixing_time, delta) keep nothing defaults as functional necessities for cross-simulation use
- Comprehensive docstring adapted from LiouvConfig, updated to document all type parameters and the trait-based with_coherent derivation

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- Config{S,D,C,T} struct ready for dispatch migration in plans 02-04
- Module will NOT compile yet -- callers still reference old types (expected, migration continues in plan 02)
- All singleton types exported and available for dispatch patterns

## Self-Check: PASSED

- FOUND: src/structs.jl
- FOUND: src/QuantumFurnace.jl
- FOUND: .planning/phases/33-type-foundation/33-01-SUMMARY.md
- FOUND: commit 27b6669 (Task 1)
- FOUND: commit 4e6a8c5 (Task 2)

---
*Phase: 33-type-foundation*
*Completed: 2026-02-25*
