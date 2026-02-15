---
phase: 08-struct-simplification
plan: 01
subsystem: structs
tags: [julia, kwdef, union-types, sentinel-replacement, immutable-struct]

# Dependency graph
requires:
  - phase: 08-02
    provides: "Redesigned HamHam struct with fully-initialized fields (no Union{..., Nothing})"
provides:
  - "Config structs with Union{..., Nothing} defaults replacing -1 sentinels"
  - "_build_common_fields() shared constructor helper"
  - "GNS inner constructor enforcement for with_coherent=false"
  - "Immutable TrottTrott with Int-typed num_trotter_steps_per_t0"
  - "load_hamiltonian_bson() for legacy BSON deserialization"
affects: [08-03, config-validation, trotter-domain, api-cleanup]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Union{T, Nothing} = nothing for optional config fields instead of sentinel defaults"
    - "Inner constructor enforcement for invariant validation on GNS configs"
    - "Manual keyword constructor for structs with inner constructors (Julia @kwdef limitation)"
    - "BSON.parse + manual field reconstruction for legacy struct deserialization"

key-files:
  created: []
  modified:
    - "src/structs.jl"
    - "src/misc_tools.jl"
    - "src/trotter_domain.jl"
    - "src/hamiltonian.jl"
    - "src/QuantumFurnace.jl"
    - "test/test_helpers.jl"
    - "simulations/main_liouv.jl"
    - "simulations/main_thermalize.jl"

key-decisions:
  - "GNS structs use manual keyword constructor + inner constructor (not @kwdef) to enforce with_coherent=false invariant"
  - "TrottTrott.bohr_freqs name kept (not renamed to quasi_bohr_freqs) to preserve polymorphic access with HamHam"
  - "Added load_hamiltonian_bson for legacy BSON compat after 08-02 changed HamHam struct"

patterns-established:
  - "isnothing(field) || field <= 0 pattern for validating optional numeric config fields"
  - "filter(p -> p[2] !== nothing, params) for excluding unset optional fields from display"

# Metrics
duration: 11min
completed: 2026-02-15
---

# Phase 08 Plan 01: Config Struct Deduplication and TrottTrott Immutability Summary

**Config sentinel defaults replaced with Union{..., Nothing}, GNS constructors enforce invariants, TrottTrott made immutable with Int field type**

## Performance

- **Duration:** 11 min
- **Started:** 2026-02-15T08:58:06Z
- **Completed:** 2026-02-15T09:09:34Z
- **Tasks:** 2
- **Files modified:** 8

## Accomplishments
- Replaced all 5 sentinel default fields (-1 / -1.) with Union{..., Nothing} = nothing across 4 config structs
- Added _build_common_fields() shared constructor helper centralizing default logic
- Added inner constructor enforcement on LiouvConfigGNS and ThermalizeConfigGNS to reject with_coherent=true
- Made TrottTrott immutable and changed num_trotter_steps_per_t0 from Float64 to Int
- Normalized ThermalizeConfig.gaussian_parameters type to match other config structs
- Added load_hamiltonian_bson() for legacy BSON file compatibility
- All 224 tests pass with zero regressions

## Task Commits

Each task was committed atomically:

1. **Task 1: Deduplicate config structs and replace sentinel defaults** - `011e706` (feat)
2. **Task 2: Make TrottTrott immutable with correct field types** - `1224dfb` (feat)

## Files Created/Modified
- `src/structs.jl` - Config structs with Union{..., Nothing} defaults, _build_common_fields helper, GNS inner constructors
- `src/misc_tools.jl` - Updated validation checks (isnothing guards) and print_press filters (nothing instead of -1.0)
- `src/trotter_domain.jl` - TrottTrott changed to immutable struct, num_trotter_steps_per_t0 typed as Int
- `src/hamiltonian.jl` - Added load_hamiltonian_bson() and finalize_hamham() backward compat
- `src/QuantumFurnace.jl` - Export list updated with finalize_hamham, load_hamiltonian_bson
- `test/test_helpers.jl` - Uses _load_test_hamiltonian for legacy BSON (from prior 08-02 work)
- `simulations/main_liouv.jl` - Updated to use load_hamiltonian with beta kwarg (from prior 08-02 work)
- `simulations/main_thermalize.jl` - Updated to use load_hamiltonian with beta kwarg (from prior 08-02 work)

## Decisions Made
- **GNS constructors use manual keyword constructors:** Julia @kwdef with inner constructors generates outer constructors that call positional args without the type parameter, causing MethodError. Solution: define struct without @kwdef, write explicit keyword constructor that delegates to inner constructor.
- **TrottTrott.bohr_freqs name kept:** Renaming to quasi_bohr_freqs would break polymorphic `ham_or_trott.bohr_freqs` access in furnace_utensils.jl. Per CONTEXT.md decision: keep as bohr_freqs.
- **Added load_hamiltonian_bson:** The 08-02 commit changed HamHam struct (removed Union{..., Nothing} fields), making BSON.load of legacy files fail. Added a BSON.parse-based loader that reconstructs HamHam from raw field data.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Legacy BSON deserialization failure after HamHam struct change**
- **Found during:** Task 1 (test suite execution)
- **Issue:** The prior 08-02 commit changed HamHam to remove Union{..., Nothing} fields for bohr_freqs/bohr_dict/gibbs. BSON.load of legacy .bson files fails because stored `nothing` values cannot convert to the new non-Nothing field types.
- **Fix:** Added load_hamiltonian_bson() in hamiltonian.jl using BSON.parse + manual field reconstruction via HamHam(NamedTuple, beta) constructor. Updated test_helpers.jl to use this approach. Added finalize_hamham backward-compat entry point.
- **Files modified:** src/hamiltonian.jl, src/QuantumFurnace.jl, test/test_helpers.jl
- **Verification:** All 224 tests pass
- **Committed in:** 011e706 (Task 1 commit)

**2. [Rule 1 - Bug] @kwdef + inner constructor incompatibility for GNS structs**
- **Found during:** Task 1 (GNS constructor enforcement)
- **Issue:** @kwdef generates an outer constructor that calls `StructName(args...)` without the type parameter, but inner constructors are defined as `StructName{D}(args...)`. The unparameterized call has no matching method.
- **Fix:** Replaced @kwdef on GNS structs with plain struct + explicit keyword outer constructor that infers D from the domain argument and delegates to the inner constructor.
- **Files modified:** src/structs.jl
- **Verification:** GNS construction with default with_coherent=false works; with_coherent=true correctly throws error
- **Committed in:** 011e706 (Task 1 commit)

**3. [Rule 1 - Bug] ThermalizeConfig gaussian_parameters type inconsistency**
- **Found during:** Task 1 (config struct update)
- **Issue:** ThermalizeConfig used `Tuple{Union{Float64, Nothing}, Union{Float64, Nothing}}` while other 3 configs used `Union{Tuple{Float64, Float64}, Tuple{Nothing, Nothing}}`. Inconsistent types could cause dispatch issues.
- **Fix:** Normalized ThermalizeConfig gaussian_parameters to `Union{Tuple{Float64, Float64}, Tuple{Nothing, Nothing}}` matching the other config types.
- **Files modified:** src/structs.jl
- **Committed in:** 011e706 (Task 1 commit)

---

**Total deviations:** 3 auto-fixed (2 bugs, 1 blocking)
**Impact on plan:** All auto-fixes necessary for correctness. The BSON fix addressed a pre-existing incompatibility from 08-02; the @kwdef fix addressed a Julia language limitation; the type normalization fixed an inconsistency. No scope creep.

## Issues Encountered
- Plan 08-02 was executed before 08-01 (out of order), leaving uncommitted changes in test_helpers.jl and simulation files. These were included in Task 1 commit to maintain a working codebase.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- Config structs ready for type parametrization cleanup (08-03)
- Domain dispatch refactoring can proceed with the clean Union{..., Nothing} field types
- TrajectoryFramework type parameter simplification (08-03) is independent of these changes

## Self-Check: PASSED

- All 8 modified files exist on disk
- Commits 011e706 and 1224dfb verified in git log
- Zero sentinel defaults (= -1) in src/structs.jl
- Zero sentinel checks (!= -1 / == -1) in src/misc_tools.jl
- 224/224 tests pass

---
*Phase: 08-struct-simplification*
*Completed: 2026-02-15*
