# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-14)

**Core value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers
**Current focus:** v1.1 Reduce -- Phase 10: API Surface Cleanup (in progress)

## Current Position

Phase: 10 of 11 (API Surface Cleanup) -- IN PROGRESS
Plan: 1 of 3 in current phase -- COMPLETE
Status: Executing 10-01
Last activity: 2026-02-15 - Completed 10-01 (Export Curation)

Progress: [##################..] 90% (v1.0 complete, v1.1 7/8 phases)

## Performance Metrics

**Velocity:**
- Total plans completed: 25 (v1.0: 10, v1.1: 12, quick: 3)
- Average duration: --
- Total execution time: --

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| v1.0 phases 1-5 | 10 | -- | -- |
| 06-dead-code-pruning | 2 | 16min | 8min |
| 07-dry-refactoring | 2 | 8min | 4min |
| quick-13 | 1 | 2min | 2min |
| quick-15 | 1 | 2min | 2min |
| 08-struct-simplification | 3 | 36min | 12min |
| 09-type-parameterization | 3/3 | 22min | 7.3min |
| quick-16 | 1 | 2min | 2min |
| 10-api-surface-cleanup | 1/3 | 6min | 6min |

*Updated after each plan completion*

## Accumulated Context

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.

- 06-01: Preserved inline single-line debug assertions within active functions
- 06-01: Preserved 4-line non-Hermitian alternative in coherent.jl as design documentation
- 06-01: linearmaps_liouv.jl excluded per CONTEXT.md explicit keep list
- 06-02: Preserved all qi_tools.jl functions except are_we_tp (pedagogical/user-facing value)
- 06-02: Left errors.jl as placeholder for Phase 10 API cleanup
- 06-02: Preserved coherent.jl support functions (compute_b_minus, etc.) used by furnace_utensils.jl
- 07-01: hermitianize! modifies rho_next in-place then copyto! for evolving_dm sites
- 07-01: coherent.jl B_trotter single-jump transforms left untouched (different pattern)
- 07-02: apply_cptp_channel! expects scratch.R pre-Hermitianized; hermitianize! remains at call site
- 07-02: apply_coherent_unitary! marked @inline for zero-overhead nothing dispatch
- quick-13: hermitianize!(scratch.tmp2) added before eigen to handle floating-point asymmetry in S matrix
- 08-01: GNS structs use manual keyword constructor + inner constructor (not @kwdef) for with_coherent enforcement
- 08-01: TrottTrott.bohr_freqs name kept for polymorphic access with HamHam
- 08-01: Added load_hamiltonian_bson for legacy BSON compat after 08-02 changed HamHam struct
- 08-02: Inline Gibbs computation via _gibbs_in_eigen helper (avoids circular dependency)
- 08-02: BSON.parse + raise_recursive for legacy deserialization (no BSON file re-saving)
- 08-02: load_hamiltonian requires beta kwarg (breaking change, clean break per CONTEXT.md)
- 08-02: find_ideal_heisenberg returns NamedTuple for composable HamHam construction
- 08-03: TrajectoryFramework reduced from {T,C,H,PD,D} to {T,D} -- only dispatch-relevant params kept
- 08-03: Domain dispatch via config type param: f(config::AbstractConfig{D}) instead of f(::DomainType, config)
- 08-03: Config domain::D field retained for runtime isa checks and display
- 09-01: HamHam.data kept as Matrix{Complex{T}} (not Hermitian wrapper) for downstream compatibility
- 09-01: TrottTrott Trotter computation stays Float64 internally, converts to T at constructor level
- 09-01: group_hamiltonian_terms generalized with Complex{T}/T local variables to match HamHam{T}
- 09-01: create_bohr_dict generalized with zero(T) key for generic AbstractFloat support
- 09-01: NamedTuple constructor beta widened to Real for type inference flexibility
- 09-02: AbstractConfig{D,T} backward compat -- existing AbstractConfig{D} dispatch matches {D,T} where T
- 09-02: JumpOp widened to AbstractMatrix{<:Complex} rather than adding explicit T parameter
- 09-02: GNS keyword constructors default beta/sigma to Float64 literals for ergonomic backward compat
- 09-02: NUFFTPrefactors promotes inputs to Float64 for FINUFFT, converts results back to Complex{T}
- 09-03: Function signatures use <:Complex for dispatch, eltype() for internal allocations
- 09-03: Cross-struct T mismatch check at run_lindbladian and run_thermalization entry points
- 09-03: TrajectoryFramework step parameters (delta, delta_eff, alpha) stay Float64 for numerical precision
- 09-03: Domain helper functions (create_alpha, create_f, etc.) widened from Float64 to Real
- quick-16: Removed convenience constructor since only construction site already uses explicit {T}
- 10-01: De-exported internal types/functions require qualified QuantumFurnace.name access in tests

### Pending Todos

None

### Blockers/Concerns

None

### Quick Tasks Completed

| # | Description | Date | Commit | Directory |
|---|-------------|------|--------|-----------|
| 13 | Unify residual Cholesky computation: compare cholesky() vs eigendecomposition approaches and use the more robust one in both DM and trajectory simulators | 2026-02-15 | 5bd9dbe | [13-unify-residual-cholesky-computation-comp](./quick/13-unify-residual-cholesky-computation-comp/) |
| 15 | Remove unused _build_common_fields() helper from src/structs.jl | 2026-02-15 | adf5398 | [15-remove-unused-build-common-fields-helper](./quick/15-remove-unused-build-common-fields-helper/) |
| 16 | Remove LindbladianWorkspace default Float64 convenience constructor | 2026-02-15 | 1dba871 | [16-defer-lindbladianworkspace-construction-](./quick/16-defer-lindbladianworkspace-construction-/) |

## Session Continuity

Last session: 2026-02-15
Stopped at: Completed 10-01-PLAN.md (Export Curation)
Resume file: None
