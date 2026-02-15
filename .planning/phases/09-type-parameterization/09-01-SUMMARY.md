---
phase: 09-type-parameterization
plan: 01
subsystem: core-structs
tags: [type-parameterization, AbstractFloat, generics, HamHam, TrottTrott, precision]

# Dependency graph
requires:
  - phase: 08-struct-simplification
    provides: Simplified HamHam and TrottTrott struct definitions, _gibbs_in_eigen helper
provides:
  - HamHam{T<:AbstractFloat} parameterized struct with precision kwarg support
  - TrottTrott{T<:AbstractFloat} parameterized struct inferring T from HamHam
  - Generic _gibbs_in_eigen, create_bohr_dict, rescaling_and_shift_factors
  - Mixed-precision policy enforcement (downward error, upward promotion)
affects: [09-02, 09-03, downstream-simulation-functions]

# Tech tracking
tech-stack:
  added: []
  patterns: [type-parameterization-on-AbstractFloat, precision-kwarg-pattern, Float64-compute-then-convert]

key-files:
  created: []
  modified:
    - src/hamiltonian.jl
    - src/trotter_domain.jl
    - src/bohr_domain.jl
    - src/misc_tools.jl

key-decisions:
  - "HamHam data field kept as Matrix{Complex{T}} (not Hermitian wrapper) for downstream compatibility"
  - "TrottTrott Trotter computation stays in Float64 internally, converts to T at constructor level"
  - "group_hamiltonian_terms and compute_U_group convert terms to Float64 for Pauli matrix compatibility"
  - "create_bohr_dict generalized with zero(T) instead of 0.0 literal"
  - "NamedTuple constructor infers T from eigvals eltype; beta widened to Real for flexibility"

patterns-established:
  - "precision=Float64 kwarg: top-level constructors accept precision=Type{T} defaulting to Float64"
  - "Compute-in-Float64: Trotter and sparse Hamiltonian construction always use Float64 intermediates, convert at boundaries"
  - "Type inference chain: HamHam{T} -> TrottTrott{T} (internals infer from HamHam)"

# Metrics
duration: 7min
completed: 2026-02-15
---

# Phase 9 Plan 01: Core Struct Type Parameterization Summary

**HamHam{T} and TrottTrott{T} parameterized on AbstractFloat with precision=Float32 kwarg, mixed-precision policy, and Float64 backward compatibility**

## Performance

- **Duration:** 7 min
- **Started:** 2026-02-15T10:27:56Z
- **Completed:** 2026-02-15T10:34:49Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments
- HamHam parameterized as HamHam{T<:AbstractFloat} with all numeric fields using T/Complex{T}
- TrottTrott parameterized as TrottTrott{T<:AbstractFloat}, inferring T from HamHam
- Mixed-precision policy enforced: precision=Float32 with Float64 data throws ArgumentError; Float32 data to Float64 silently promotes
- All 224 existing tests pass unchanged with default HamHam{Float64}/TrottTrott{Float64} construction
- BSON loader naturally produces HamHam{Float64} without code changes

## Task Commits

Each task was committed atomically:

1. **Task 1: Parameterize HamHam{T} struct and constructors** - `93aa0d0` (feat)
2. **Task 2: Parameterize TrottTrott{T} and update BSON loader** - `b9cd7df` (feat)

**Plan metadata:** [pending] (docs: complete plan)

## Files Created/Modified
- `src/hamiltonian.jl` - HamHam{T} struct, _gibbs_in_eigen{T}, constructors with precision kwarg, rescaling_and_shift_factors generalized
- `src/trotter_domain.jl` - TrottTrott{T} struct, constructor with Float64 compute + T conversion, group_hamiltonian_terms generalized
- `src/bohr_domain.jl` - create_bohr_dict generalized to Matrix{T} with zero(T) key
- `src/misc_tools.jl` - BSON loader docstrings updated to reflect HamHam{Float64} return type

## Decisions Made
- **HamHam.data as Matrix{Complex{T}}**: Kept as Matrix rather than Hermitian wrapper because all downstream code calls `size(hamiltonian.data)` which works with Matrix. The Hermitian wrapping is stripped by Julia's constructor conversion.
- **Trotter Float64 compute path**: trotterize2 and helper functions always compute in Float64 since Pauli matrices (X, Y, Z) are ComplexF64 constants. Results are converted to T at TrottTrott constructor boundary. This avoids deep changes to Pauli constant types.
- **group_hamiltonian_terms generic**: Local variable type annotations changed from hardcoded ComplexF64/Float64 to Complex{T}/T to match HamHam{T} field types.
- **compute_U_group/expm_pauli_padded bridge**: Terms are explicitly converted to Vector{Matrix{ComplexF64}} before passing to expm_pauli_padded, which expects ComplexF64 (Pauli matrix arithmetic).
- **NamedTuple constructor beta::Real**: Widened from Float64 to Real to allow beta inference alongside T from eigvals.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Converted terms to Float64 in trotterize2 and compute_U_group**
- **Found during:** Task 2 (TrottTrott parameterization)
- **Issue:** expm_pauli_padded requires Vector{Matrix{ComplexF64}} but HamHam{T} stores Vector{Matrix{Complex{T}}}. For T=Float64 this is fine, but for generic T the types wouldn't match.
- **Fix:** Added explicit `Vector{Matrix{ComplexF64}}(term)` conversions in trotterize2, compute_U_group, and trotterize before calling expm_pauli_padded.
- **Files modified:** src/trotter_domain.jl
- **Verification:** All 224 tests pass
- **Committed in:** b9cd7df (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (1 blocking)
**Impact on plan:** Auto-fix necessary for type-safe Float64 Trotter computation with generic HamHam{T}. No scope creep.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- HamHam{T} and TrottTrott{T} provide the type parameter infrastructure for Plan 02 (Config/Workspace parameterization) and Plan 03 (simulation function generics)
- Downstream code using `::HamHam` dispatch annotations works unchanged since HamHam{Float64} <: HamHam
- The precision=Float32 kwarg is accepted but Float32 end-to-end paths require Plan 02/03 to parameterize configs and workspaces

## Self-Check: PASSED

- All 4 modified files exist on disk
- SUMMARY.md created at expected path
- Commit 93aa0d0 (Task 1) found in git log
- Commit b9cd7df (Task 2) found in git log
- All 224 tests pass

---
*Phase: 09-type-parameterization*
*Completed: 2026-02-15*
