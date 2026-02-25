# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-25)

**Core value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers
**Current focus:** v2.0 Restructure -- Phase 33: Type Foundation

## Current Position

Phase: 33 of 38 (Type Foundation)
Plan: 1 of 4 in current phase
Status: Executing
Last activity: 2026-02-25 -- Completed 33-01 (type hierarchies, Config struct, exports)

Progress: [##░░░░░░░░] ~4% (v2.0, 1/~24 plans)

## Performance Metrics

**Velocity:**
- Total plans completed: 88 (v1.0: 10, v1.1: 16, quick: 24, v1.2: 12, cleanup: 3, v1.3: 10, v1.4: 2, v1.5: 12)

**By Milestone:**

| Milestone | Phases | Plans | Timeline |
|-----------|--------|-------|----------|
| v1.0 Trajectories | 1-5 | 10 | 2026-02-13 to 2026-02-14 |
| v1.1 Reduce | 6-11 | 16 (+5 quick) | 2026-02-15 |
| v1.2 Multi-threading | 12-19 | 15 (+3 quick) | 2026-02-15 to 2026-02-16 |
| v1.3 Mixing Time | 20-25 | 10 (+11 quick) | 2026-02-17 to 2026-02-18 |
| v1.4 Spectral Gap Refinement | 26 | 2 (+1 quick) | 2026-02-19 to 2026-02-20 |
| v1.5 Krylov Gap Estimation | 27-32 | 12 (+3 quick) | 2026-02-20 to 2026-02-25 |
| v2.0 Restructure | 33-38 | 1/~24 | 2026-02-25 to ... |

## Accumulated Context

### Decisions

All decisions logged in PROJECT.md Key Decisions table.

- [33-01] Field ordering: singletons first, then system/physics/grid/thermalize-specific params
- [33-01] with_coherent derived from construction type via trait, not stored as field
- [33-01] Outer constructor infers S,D,C from singleton args and T from beta

### Key Constraints for v2.0

- Keep time_oft!/trotter_oft! as test utilities (don't delete)
- Keep OFTCaches for testing reference
- R/K0/U_residual: per-jump (R^a) for DM/Trajectory, summed (R_total) for Krylov -- don't mix
- Config is `Config{S,D,C,T}` where C = Construction (KMS, GNS, DLL future)
- Stay flat in src/ (no subdirectories) -- rename files for clarity
- Diagnostics stays as separate module
- SharedArrays stays, only @distributed code is dead
- Bohr domain may resist full unification (different loop structure)

### Pending Todos

None

### Blockers/Concerns

None

## Session Continuity

Last session: 2026-02-25
Stopped at: Completed 33-01-PLAN.md (type hierarchies, Config struct, exports)
Resume file: None
