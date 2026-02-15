# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-15)

**Core value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers
**Current focus:** v1.2 Multi-threading -- Phase 12: Workspace Refactor

## Current Position

Phase: 12 of 18 (Workspace Refactor) -- first phase of v1.2
Plan: 2 of 2 in current phase (COMPLETE)
Status: Phase 12 Complete
Last activity: 2026-02-15 -- Completed 12-02-PLAN.md (test migration + workspace independence)

Progress: [██████████████████████░░░░░░░░] 76% (28/TBD plans -- v1.0: 10, v1.1: 16+5q, v1.2: 2/TBD)

## Performance Metrics

**Velocity:**
- Total plans completed: 33 (v1.0: 10, v1.1: 16, quick: 5, v1.2: 2)

**By Milestone:**

| Milestone | Phases | Plans | Timeline |
|-----------|--------|-------|----------|
| v1.0 Trajectories | 1-5 | 10 | 2026-02-13 to 2026-02-14 |
| v1.1 Reduce | 6-11 | 16 (+5 quick) | 2026-02-15 |
| v1.2 Multi-threading | 12-18 | 2/TBD | 2026-02-15 to ... |

## Accumulated Context

### Decisions

All decisions logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- Phase 12 complete: workspace independence verified, all 246 tests pass with explicit ws/rng
- Workspace refactor (THRD-01) complete: TrajectoryFramework is read-only, workspace passed explicitly
- TrajectoryResult struct with optional times/measurements_mean fields (single return type)
- store_states kwarg dropped from run_trajectories (unused, not in TrajectoryResult)
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
Stopped at: Completed 12-02-PLAN.md (Phase 12 complete)
Resume file: None
