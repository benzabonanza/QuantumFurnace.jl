# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-15)

**Core value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers
**Current focus:** v1.2 Multi-threading -- Phase 13 complete, Phase 14 ready

## Current Position

Phase: 13 of 18 (Multi-Threaded Trajectory Engine) -- COMPLETE
Plan: 2/2 complete
Status: Phase 13 complete. All threading tests pass (257 tests with 2 threads). Ready for Phase 14.
Last activity: 2026-02-16 -- Phase 13 Plan 02 executed (2 tasks, 2 commits)

Progress: [███████████████████████░░░░░░░] 79% (30/TBD plans -- v1.0: 10, v1.1: 16+5q, v1.2: 4/TBD)

## Performance Metrics

**Velocity:**
- Total plans completed: 35 (v1.0: 10, v1.1: 16, quick: 5, v1.2: 4)

**By Milestone:**

| Milestone | Phases | Plans | Timeline |
|-----------|--------|-------|----------|
| v1.0 Trajectories | 1-5 | 10 | 2026-02-13 to 2026-02-14 |
| v1.1 Reduce | 6-11 | 16 (+5 quick) | 2026-02-15 |
| v1.2 Multi-threading | 12-18 | 4/TBD | 2026-02-15 to ... |

## Accumulated Context

### Decisions

All decisions logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- Phase 13 Plan 02: Performance test uses ntraj=200 x mixing_time=5.0 (not 50 x 1.0) to amortize threading overhead on dim=8
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
Stopped at: Completed 13-02-PLAN.md. Phase 13 complete. Phase 14 (convergence checking) ready.
Resume file: None
