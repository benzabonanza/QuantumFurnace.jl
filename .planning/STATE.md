# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-14)

**Core value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers
**Current focus:** v1.1 Reduce -- Phase 7: DRY Refactoring (Plan 1 complete)

## Current Position

Phase: 7 of 11 (DRY Refactoring)
Plan: 1 of 2 in current phase
Status: Plan 07-01 complete
Last activity: 2026-02-15 -- Completed 07-01 (hermitianize! and basis transform helpers)

Progress: [############........] 59% (v1.0 complete, v1.1 2/6 phases in progress)

## Performance Metrics

**Velocity:**
- Total plans completed: 13 (v1.0: 10, v1.1: 3)
- Average duration: --
- Total execution time: --

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| v1.0 phases 1-5 | 10 | -- | -- |
| 06-dead-code-pruning | 2 | 16min | 8min |
| 07-dry-refactoring | 1 | 4min | 4min |

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

### Pending Todos

None

### Blockers/Concerns

None

## Session Continuity

Last session: 2026-02-15
Stopped at: Completed 07-01-PLAN.md
Resume file: None
