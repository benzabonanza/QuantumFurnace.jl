---
phase: quick-17
plan: 01
subsystem: structs
tags: [julia, kwdef, constructors, config-structs]

# Dependency graph
requires:
  - phase: 08-struct-simplification
    provides: "Config struct hierarchy with AbstractConfig{D,T}"
  - phase: 09-type-parameterization
    provides: "Type-parameterized config structs with T<:AbstractFloat"
provides:
  - "All 4 config structs using consistent @kwdef pattern"
  - "Eliminated standalone keyword constructor boilerplate for GNS configs"
affects: []

# Tech tracking
tech-stack:
  added: []
  patterns: ["@kwdef + inner constructor + bridge outer for validated structs"]

key-files:
  created: []
  modified: [src/structs.jl]

key-decisions:
  - "Bridge outer constructor needed for @kwdef with inner-constructor-only structs"

patterns-established:
  - "@kwdef + inner constructor + bridge: pattern for structs needing validation with keyword construction"

# Metrics
duration: 1min
completed: 2026-02-15
---

# Quick Task 17: Simplify Config Constructors Summary

**Converted GNS config structs from 3-part constructor boilerplate to @kwdef with inner constructor validation, reducing 32 lines across LiouvConfigGNS and ThermalizeConfigGNS**

## Performance

- **Duration:** 1 min
- **Started:** 2026-02-15T14:04:51Z
- **Completed:** 2026-02-15T14:06:16Z
- **Tasks:** 2
- **Files modified:** 1

## Accomplishments
- LiouvConfigGNS converted from struct + inner + standalone keyword constructor to @kwdef + inner + bridge (16 lines removed)
- ThermalizeConfigGNS converted using same pattern (16 lines removed)
- All 4 config structs (LiouvConfig, LiouvConfigGNS, ThermalizeConfig, ThermalizeConfigGNS) now use consistent @kwdef pattern
- with_coherent=true validation preserved at construction time for both GNS configs
- Keyword API identical: same field names, same defaults, same required fields

## Task Commits

Each task was committed atomically:

1. **Task 1: Convert LiouvConfigGNS to @kwdef with inner constructor** - `ada4592` (feat)
2. **Task 2: Convert ThermalizeConfigGNS to @kwdef with inner constructor** - `be8597b` (feat)

## Files Created/Modified
- `src/structs.jl` - All 4 config structs now use @kwdef; GNS configs simplified from 3-part to @kwdef+inner+bridge pattern

## Decisions Made
- Bridge outer constructor pattern: @kwdef generates `StructName(args...)` (unparameterized) but inner constructor only defines `StructName{D,T}(args...)`. The bridge function routes unparameterized positional calls to the parameterized inner constructor, which is required for @kwdef keyword construction to work with inner-constructor-only structs.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Config struct cleanup complete
- All 4 config structs now follow consistent @kwdef pattern
- No blockers

## Self-Check: PASSED

- [x] src/structs.jl exists and contains 4 @kwdef config structs
- [x] Commit ada4592 found (Task 1: LiouvConfigGNS)
- [x] Commit be8597b found (Task 2: ThermalizeConfigGNS)
- [x] Package loads without errors
- [x] Both GNS configs reject with_coherent=true at construction time

---
*Quick Task: 17-simplify-config-constructors*
*Completed: 2026-02-15*
