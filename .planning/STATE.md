# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-03-01)

**Core value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers
**Current focus:** v2.1 Speedup & Mixing Time -- Phase 42 complete (all plans executed, 1246 tests pass)

## Current Position

Phase: 42 of 42 (Mixing Time Estimation)
Plan: 2 of 2 in current phase (COMPLETE)
Status: Phase 42 complete -- all mixing time estimation tests pass (1246 total, +65 new)
Last activity: 2026-03-01 -- Phase 42 plan 02 executed (2 tasks, full test suite green)

Progress: [██████████] 100%

## Performance Metrics

**Velocity:**
- Total plans completed: 115 (v1.0: 10, v1.1: 16, quick: 26, v1.2: 12, cleanup: 3, v1.3: 10, v1.4: 2, v1.5: 12, v2.0: 19, v2.1: 7)

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
| v2.1 Speedup & Mixing Time | 39-42 | 7 | 2026-03-01 |

## Accumulated Context

### Decisions

All decisions logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- v2.1 scope: Per-jump precomputation, save_every, BLAS/omega threading, mixing time estimation
- BohrDomain: NO per-Bohr-frequency precomputation (frequency count grows too fast); general speedups and threading only
- Omega-loop threading: Implemented in Phase 41-02 -- OMEGA_THREAD_THRESHOLD=50, per-task ThermalizeScratch, BLAS=1 during Julia-level parallelism
- Mixing time estimation: Post-processing function, not embedded in run_thermalize
- MixingTimeEstimate is separate struct (not ThermalizeResults modification) -- preserves BSON compatibility
- skip_initial passed through to fit_exponential_decay (not applied twice) -- Research pitfall 3 option b
- Primary mixing_time: extrapolated when extrapolate=true, actual crossing when target_epsilon, total sim time otherwise
- Mixing time tests: synthetic ThermalizeResults via helper function with zero DM for clean unit testing
- Integration test uses N3 (3-qubit) with 5s mixing_time -- proves end-to-end API works
- _precompute_per_jump_channels stores K0s/U_residuals only (no Rs) -- DM path does not need raw R matrices
- BohrDomain _precompute_R uses precomputed bohr_is/bohr_js with fallback to hamiltonian.bohr_dict
- jump_weight_scaling precomputed before hot loop as gamma_norm_factor / p_jump
- B_bohr/B_time/B_trotter signatures use AbstractVector{<:JumpOp} for Julia type invariance correctness
- save_every gates observation only (trace_distance_h, push, printf, convergence); physics runs every step unconditionally
- recorded_steps array replaces num_steps-based time_steps construction for save_every support
- Convergence cutoff checked only at save points (coarser detection acceptable for reduced observation cost)
- DM BLAS pattern: save -> set_num_threads(nthreads) -> try loop finally restore (inverse of trajectory pattern which sets BLAS=1)
- Only hot loop wrapped in try/finally -- precomputation and result construction remain outside

### Pending Todos

None

### Blockers/Concerns

None -- Phase 42 complete. v2.1 milestone fully delivered.

## Session Continuity

Last session: 2026-03-01
Stopped at: Completed 42-02-PLAN.md (mixing time estimation tests). Phase 42 and v2.1 milestone complete.
Resume file: None
