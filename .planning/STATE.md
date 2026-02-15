# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-15)

**Core value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers
**Current focus:** v1.2 Multi-threading -- Phase 12: Workspace Refactor

## Current Position

Phase: 12 of 18 (Workspace Refactor) -- first phase of v1.2
Plan: 0 of TBD in current phase
Status: Ready to plan
Last activity: 2026-02-15 -- Roadmap created for v1.2 Multi-threading (phases 12-18)

Progress: [██████████████████████░░░░░░░░] 73% (26/TBD plans -- v1.0: 10, v1.1: 16+5q, v1.2: 0/TBD)

## Performance Metrics

**Velocity:**
- Total plans completed: 31 (v1.0: 10, v1.1: 16, quick: 5)

**By Milestone:**

| Milestone | Phases | Plans | Timeline |
|-----------|--------|-------|----------|
| v1.0 Trajectories | 1-5 | 10 | 2026-02-13 to 2026-02-14 |
| v1.1 Reduce | 6-11 | 16 (+5 quick) | 2026-02-15 |
| v1.2 Multi-threading | 12-18 | TBD | 2026-02-15 to ... |

## Accumulated Context

### Decisions

All decisions logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- Workspace refactor (THRD-01) is the gate for all threading work
- BSON for serialization (existing dependency, no JLD2)
- BLAS.set_num_threads(1) mandatory before threaded execution
- TaskLocalRNG with seed = master_seed + trajectory_id for reproducibility
- Primary observable: nearest-neighbor correlations <Z_iZ_{i+1}>
- Adaptive = batch convergence (relative change <1% for 3 consecutive batches, hard cap N_max)

### Pending Todos

None

### Blockers/Concerns

None

## Session Continuity

Last session: 2026-02-15
Stopped at: v1.2 roadmap created, ready to plan Phase 12
Resume file: None
