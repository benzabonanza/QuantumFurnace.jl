# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-14)

**Core value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers
**Current focus:** v1.1 Reduce -- Phase 8: Struct Simplification (in progress)

## Current Position

Phase: 8 of 11 (Struct Simplification) -- IN PROGRESS
Plan: 2 of 3 in current phase
Status: Executing phase 08
Last activity: 2026-02-15 - Completed 08-02: HamHam initialization redesign

Progress: [##############......] 70% (v1.0 complete, v1.1 4/6 phases in progress)

## Performance Metrics

**Velocity:**
- Total plans completed: 18 (v1.0: 10, v1.1: 7, quick: 1)
- Average duration: --
- Total execution time: --

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| v1.0 phases 1-5 | 10 | -- | -- |
| 06-dead-code-pruning | 2 | 16min | 8min |
| 07-dry-refactoring | 2 | 8min | 4min |
| quick-13 | 1 | 2min | 2min |
| 08-struct-simplification | 3 | 24min | 8min |

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

### Pending Todos

None

### Blockers/Concerns

None

### Quick Tasks Completed

| # | Description | Date | Commit | Directory |
|---|-------------|------|--------|-----------|
| 13 | Unify residual Cholesky computation: compare cholesky() vs eigendecomposition approaches and use the more robust one in both DM and trajectory simulators | 2026-02-15 | 5bd9dbe | [13-unify-residual-cholesky-computation-comp](./quick/13-unify-residual-cholesky-computation-comp/) |

## Session Continuity

Last session: 2026-02-15
Stopped at: Completed 08-02-PLAN.md (HamHam initialization redesign)
Resume file: None
