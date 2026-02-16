# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-15)

**Core value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers
**Current focus:** v1.2 Multi-threading -- Phase 16 complete, ready for Phase 17

## Current Position

Phase: 16 of 18 (Convergence Tracking) -- COMPLETE
Plan: 2/2 complete
Status: Phase 16 complete. Convergence tracking infrastructure + comprehensive tests (470 total). Ready for Phase 17 (Adaptive Stopping).
Last activity: 2026-02-16 - Completed 16-02: Convergence tracking tests

Progress: [████████████████████████████░░] 89% (40/TBD plans -- v1.0: 10, v1.1: 16+5q, v1.2: 9/TBD)

## Performance Metrics

**Velocity:**
- Total plans completed: 40 (v1.0: 10, v1.1: 16, quick: 5, v1.2: 9)

**By Milestone:**

| Milestone | Phases | Plans | Timeline |
|-----------|--------|-------|----------|
| v1.0 Trajectories | 1-5 | 10 | 2026-02-13 to 2026-02-14 |
| v1.1 Reduce | 6-11 | 16 (+5 quick) | 2026-02-15 |
| v1.2 Multi-threading | 12-18 | 9/TBD | 2026-02-15 to ... |

## Accumulated Context

### Decisions

All decisions logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- Phase 16 Plan 02: 470 total tests (364 existing + 106 new convergence tests), 10 testsets
- Phase 16 Plan 02: Integration test uses 1000 trajectories (200x5 batches) with mixing_time=5.0 for convergence proof
- Phase 16 Plan 02: Convergence assertion uses generous "last < first" tolerance (not strict monotonic)
- Phase 16 Plan 01: Observables built in eigenbasis (not computational) to match run_trajectories rho output basis
- Phase 16 Plan 01: Separate run_trajectories_convergence function (not modifying existing run_trajectories) to avoid API bloat
- Phase 16 Plan 01: Scalar-only ConvergenceData (no density matrix snapshots) for O(n_batches) memory
- Phase 16 Plan 01: convergence.jl included before results.jl so ConvergenceData available to serialization functions
- Phase 16 Plan 01: Pre-allocated arrays for convergence data storage instead of push!/hcat
- Phase 15 Plan 02: 364 total tests (284 existing + 80 new) -- all round-trip, forward compat, integration tests pass
- Phase 15 Plan 02: Integration test uses 3-qubit SMALL system with 10 trajectories (fast, proves full pipeline)
- Phase 15 Plan 01: Dict-based BSON serialization for ExperimentResult (avoids parametric struct pitfalls)
- Phase 15 Plan 01: TrajectoryResult embedded as field, domain singletons stored as strings
- Phase 15 Plan 01: LibGit2 for git hash capture, Dates for timestamps (both stdlibs)
- Phase 15 Plan 01: Hamiltonian params store only reproduction-relevant subset (no eigendecomposition)
- Quick-20: GNS-to-Gibbs gap at sigma=0.1: EnergyDomain=0.081, TrotterDomain=0.081, BohrDomain=0.035 (Phase 18 baselines)
- Quick-20: TrotterDomain Lindbladian is in Trotter eigenbasis; must transform fixed point to energy eigenbasis before comparing to Gibbs
- Phase 14 Plan 01: Trajectory convergence params: ntraj=1000, delta=0.01, mixing_time=5.0 -> 0.029 trace distance
- Phase 14 Plan 01: Fixed @kwdef outer constructor bug for GNS config structs (pre-existing)
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

### Quick Tasks Completed

| # | Description | Date | Commit | Directory |
|---|-------------|------|--------|-----------|
| 19 | Fix failing test after EnergyDomain to TrotterDomain rename in GNS trajectory | 2026-02-16 | 0683acc | [19-fix-failing-test-after-energydomain-to-t](./quick/19-fix-failing-test-after-energydomain-to-t/) |
| 20 | Fix basis mismatch in GNS TrotterDomain tests (0.83 -> 0.08 gap) | 2026-02-16 | 0161e01 | [20-debug-gns-trotterdomain-0-83-gap-suspect](./quick/20-debug-gns-trotterdomain-0-83-gap-suspect/) |

## Session Continuity

Last session: 2026-02-16
Stopped at: Completed 16-02-PLAN.md. Phase 16 complete (both plans). Ready for Phase 17 (Adaptive Stopping).
Resume file: None
