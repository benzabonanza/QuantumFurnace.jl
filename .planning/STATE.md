# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-15)

**Core value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers
**Current focus:** v1.2 Multi-threading -- Phase 13 Plan 01 complete, Plan 02 remaining

## Current Position

Phase: 13 of 18 (Multi-Threaded Trajectory Engine)
Plan: 1/2 complete
Status: Plan 01 (multi-threaded engine) complete. Plan 02 (thread safety testing) ready.
Last activity: 2026-02-16 -- Phase 13 Plan 01 executed (2 tasks, 2 commits)

Progress: [██████████████████████░░░░░░░░] 78% (29/TBD plans -- v1.0: 10, v1.1: 16+5q, v1.2: 3/TBD)

## Performance Metrics

**Velocity:**
- Total plans completed: 34 (v1.0: 10, v1.1: 16, quick: 5, v1.2: 3)

**By Milestone:**

| Milestone | Phases | Plans | Timeline |
|-----------|--------|-------|----------|
| v1.0 Trajectories | 1-5 | 10 | 2026-02-13 to 2026-02-14 |
| v1.1 Reduce | 6-11 | 16 (+5 quick) | 2026-02-15 |
| v1.2 Multi-threading | 12-18 | 3/TBD | 2026-02-15 to ... |

## Accumulated Context

### Decisions

All decisions logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- Phase 13 Plan 01: Added F,P type params to TrajectoryFramework{T,D,F,P} for zero-allocation hot path
- Phase 13 Plan 01: Changed jumps from Vector{JumpOp} to Vector{JumpOp{Matrix{T}}} for concrete element access
- Phase 13 Plan 01: Per-trajectory Xoshiro(seed + traj_id) seeding replaces shared RNG in serial path
- Phase 13 Plan 01: Precomputed scaled_prefactor stored in framework to avoid abstract config access in hot loop
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

Last session: 2026-02-16
Stopped at: Completed 13-01-PLAN.md. Plan 02 (thread safety testing) ready for execution.
Resume file: None
