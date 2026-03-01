# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-03-01)

**Core value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers
**Current focus:** v2.1 Speedup & Mixing Time -- Phase 39 (Per-Jump Precomputation)

## Current Position

Phase: 39 of 42 (Per-Jump Precomputation)
Plan: 1 of 2 in current phase
Status: In progress
Last activity: 2026-03-01 -- Completed Plan 39-01 (precomputation infrastructure)

Progress: [█████░░░░░] 50%

## Performance Metrics

**Velocity:**
- Total plans completed: 109 (v1.0: 10, v1.1: 16, quick: 26, v1.2: 12, cleanup: 3, v1.3: 10, v1.4: 2, v1.5: 12, v2.0: 19)

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
| v2.1 Speedup & Mixing Time | 39-42 | 1/TBD | 2026-03-01 to present |

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

### Pending Todos

None

### Blockers/Concerns

- BohrDomain _precompute_R has no existing codebase precedent -- requires careful extraction of rho-independent R part from _jump_contribution!
- Regression test baselines may need regeneration after precomputation (O(1e-14) level eigendecomposition FP variation)

## Session Continuity

Last session: 2026-03-01
Stopped at: Completed 39-01-PLAN.md (precomputation infrastructure). Ready for 39-02.
Resume file: None
