# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-03-01)

**Core value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers
**Current focus:** v2.1 Speedup & Mixing Time -- Phase 40 complete, ready for Phase 41 (BLAS/Omega Threading)

## Current Position

Phase: 40 of 42 (Save Every) -- COMPLETE
Plan: 1 of 1 in current phase (COMPLETE)
Status: Phase 40 complete -- save_every keyword implemented and tested
Last activity: 2026-03-01 -- Phase 40 plan 01 executed (2 tasks, 1158 tests pass)

Progress: [██████████] 100%

## Performance Metrics

**Velocity:**
- Total plans completed: 111 (v1.0: 10, v1.1: 16, quick: 26, v1.2: 12, cleanup: 3, v1.3: 10, v1.4: 2, v1.5: 12, v2.0: 19, v2.1: 3)

**By Milestone:**

| Milestone | Phases | Plans | Timeline |
|-----------|--------|-------|----------|
| v1.0 Trajectories | 1-5 | 10 | 2026-02-13 to 2026-02-14 |
| v1.1 Reduce | 6-11 | 16 (+5 quick) | 2026-02-15 |
| v1.2 Multi-threading | 12-19 | 15 (+3 quick) | 2026-02-15 to 2026-02-16 |
| v1.3 Mixing Time | 20-25 | 10 (+11 quick) | 2026-02-17 to 2026-02-18 |
| v1.4 Spectral Gap Refinement | 26 | 2 (+1 quick) | 2026-02-19 to 2026-02-20 |
| v1.5 Krylov Gap Estimation | 27-32 | 12 (+3 quick) | 2026-02-20 to 2026-02-25 |
| v2.0 Restructure | 33-38 | 19 (+2 quick) | 2026-02-25 to 2026-02-28 |
| v2.1 Speedup & Mixing Time | 39-42 | 3/TBD | 2026-03-01 to present |

## Accumulated Context

### Decisions

All decisions logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- v2.1 scope: Per-jump precomputation, save_every, BLAS/omega threading, mixing time estimation
- BohrDomain: NO per-Bohr-frequency precomputation (frequency count grows too fast); general speedups and threading only
- Omega-loop threading: Optional/deferred within Phase 41 -- BLAS threading alone provides meaningful speedup
- Mixing time estimation: Post-processing function, not embedded in run_thermalize
- _precompute_per_jump_channels stores K0s/U_residuals only (no Rs) -- DM path does not need raw R matrices
- BohrDomain _precompute_R uses precomputed bohr_is/bohr_js with fallback to hamiltonian.bohr_dict
- jump_weight_scaling precomputed before hot loop as gamma_norm_factor / p_jump
- B_bohr/B_time/B_trotter signatures use AbstractVector{<:JumpOp} for Julia type invariance correctness
- save_every gates observation only (trace_distance_h, push, printf, convergence); physics runs every step unconditionally
- recorded_steps array replaces num_steps-based time_steps construction for save_every support
- Convergence cutoff checked only at save points (coarser detection acceptable for reduced observation cost)

### Pending Todos

None

### Blockers/Concerns

None -- Phase 40 complete, no outstanding blockers.

## Session Continuity

Last session: 2026-03-01
Stopped at: Completed 40-01-PLAN.md (save_every feature). Ready for Phase 41 (BLAS/Omega Threading).
Resume file: None
