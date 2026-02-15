---
phase: 09-type-parameterization
plan: 02
subsystem: core-structs
tags: [type-parameterization, AbstractFloat, generics, Config, LindbladianWorkspace, JumpOp, KrausScratch, NUFFTPrefactors]

# Dependency graph
requires:
  - phase: 09-type-parameterization
    plan: 01
    provides: HamHam{T} and TrottTrott{T} parameterized structs
provides:
  - Config structs (LiouvConfig, LiouvConfigGNS, ThermalizeConfig, ThermalizeConfigGNS) parameterized as {D,T}
  - LindbladianWorkspace{T} with Complex{T} buffers
  - OFTCaches{T} with Complex{T} fields
  - JumpOp widened to accept any AbstractMatrix{<:Complex}
  - HotAlgorithmResults{D,T} and HotSpectralResults{D,T} referencing HamHam{T}/TrottTrott{T}
  - KrausScratch constructor accepting generic Complex{T} types
  - NUFFTPrefactors{T,A} parameterized on AbstractFloat T
affects: [09-03, downstream-simulation-functions]

# Tech tracking
tech-stack:
  added: []
  patterns: [config-struct-dual-parameterization-{D,T}, workspace-parameterization-on-T, NUFFT-Float64-promotion]

key-files:
  created: []
  modified:
    - src/structs.jl
    - src/kraus.jl
    - src/nufft.jl

key-decisions:
  - "AbstractConfig hierarchy carries {D,T} -- backward compatible since AbstractConfig{D} matches AbstractConfig{D,T} where T"
  - "JumpOp widened to AbstractMatrix{<:Complex} rather than adding explicit T parameter for simplicity"
  - "GNS keyword constructors default beta/sigma to Float64 literals for ergonomic backward compatibility"
  - "NUFFTPrefactors promotes all inputs to Float64 for FINUFFT, converts results back to Complex{T}"

patterns-established:
  - "Dual parameterization {D,T}: Config structs carry both domain D and float precision T"
  - "Convenience constructor pattern: StructName(args...) = StructName{Float64}(args...) for all parameterized workspace/cache structs"
  - "NUFFT Float64 promotion: FINUFFT requires Float64, so generic inputs are promoted and results downcast to T"

# Metrics
duration: 6min
completed: 2026-02-15
---

# Phase 9 Plan 02: Config/Workspace Type Parameterization Summary

**Config structs parameterized as {D,T}, LindbladianWorkspace/OFTCaches/KrausScratch/NUFFTPrefactors generalized for generic float precision**

## Performance

- **Duration:** 6 min
- **Started:** 2026-02-15T10:37:52Z
- **Completed:** 2026-02-15T10:44:16Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments
- All config structs (LiouvConfig, LiouvConfigGNS, ThermalizeConfig, ThermalizeConfigGNS) and their abstract type hierarchy parameterized on {D, T<:AbstractFloat}
- LindbladianWorkspace, OFTCaches parameterized on T with Complex{T} buffer matrices and Float64 convenience constructors
- JumpOp widened from AbstractMatrix{ComplexF64} to AbstractMatrix{<:Complex} for generic complex type support
- Result structs (HotAlgorithmResults, HotSpectralResults) reference HamHam{T} and TrottTrott{T}
- KrausScratch constructor generalized to accept any Complex{T} type instead of ComplexF64-only
- NUFFTPrefactors{T,A} parameterized on T with Float64 promotion for FINUFFT computation
- All 224 existing tests pass unchanged with Float64 inference from literal values

## Task Commits

Each task was committed atomically:

1. **Task 1: Parameterize Config structs and LindbladianWorkspace** - `270fcff` (feat)
2. **Task 2: Update KrausScratch and NUFFTPrefactors for generic types** - `e28236b` (feat)

**Plan metadata:** [pending] (docs: complete plan)

## Files Created/Modified
- `src/structs.jl` - All config structs, abstract types, LindbladianWorkspace, OFTCaches, JumpOp, result structs parameterized on T
- `src/kraus.jl` - KrausScratch constructor generalized from ComplexF64-only to any Complex{T}
- `src/nufft.jl` - NUFFTPrefactors{T,A} struct, _unique_with_invmap generalized, prepare_oft_nufft_prefactors accepts generic floats with Float64 promotion, prefactor_view untyped omega

## Decisions Made
- **AbstractConfig{D,T} backward compatibility**: Adding T as second type parameter is backward compatible in Julia because existing dispatch on `AbstractConfig{D}` implicitly means `AbstractConfig{D,T} where T`. No downstream dispatch signatures needed changing.
- **JumpOp widening approach**: Chose `AbstractMatrix{<:Complex}` rather than adding explicit float type parameter `{T, M<:AbstractMatrix{Complex{T}}}`. Simpler, still type-safe, and existing `Matrix{ComplexF64}` construction works unchanged.
- **GNS keyword constructor defaults**: `beta::T = 1.0` and `sigma::T = 0.1` in GNS keyword constructors to allow type inference from Float64 literals when no explicit type is given.
- **NUFFT Float64 promotion**: FINUFFT library only supports Float64. For generic T inputs, promote to Float64 for computation, convert results back to Complex{T}. Pragmatic approach since NUFFT is a precomputation step, not a hot path.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Complete struct-level type chain: HamHam{T} -> Config{D,T} -> LindbladianWorkspace{T} -> Results{D,T}
- KrausScratch and NUFFTPrefactors ready for generic T construction
- Plan 03 (simulation function generics) can now propagate T through function signatures and computation paths
- All existing Float64 paths work without modification

## Self-Check: PASSED

- src/structs.jl exists and contains parameterized types
- src/kraus.jl exists and contains generic KrausScratch constructor
- src/nufft.jl exists and contains parameterized NUFFTPrefactors
- Commit 270fcff (Task 1) found in git log
- Commit e28236b (Task 2) found in git log
- All 224 tests pass

---
*Phase: 09-type-parameterization*
*Completed: 2026-02-15*
