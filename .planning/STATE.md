# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-25)

**Core value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers
**Current focus:** v2.0 Restructure -- Phase 34: Code Deduplication (Plan 01 complete, Plan 02 next)

## Current Position

Phase: 34 of 38 (Code Deduplication) -- IN PROGRESS
Plan: 1/2 complete
Status: Plan 01 complete (domain_prefactor + OFT unification), Plan 02 next
Last activity: 2026-02-26 - Completed 34-01: Extract domain_prefactor and unify oft!

Progress: [#####░░░░░] ~21% (v2.0, 5/~24 plans)

## Performance Metrics

**Velocity:**
- Total plans completed: 95 (v1.0: 10, v1.1: 16, quick: 26, v1.2: 12, cleanup: 3, v1.3: 10, v1.4: 2, v1.5: 12, v2.0: 5)

**By Milestone:**

| Milestone | Phases | Plans | Timeline |
|-----------|--------|-------|----------|
| v1.0 Trajectories | 1-5 | 10 | 2026-02-13 to 2026-02-14 |
| v1.1 Reduce | 6-11 | 16 (+5 quick) | 2026-02-15 |
| v1.2 Multi-threading | 12-19 | 15 (+3 quick) | 2026-02-15 to 2026-02-16 |
| v1.3 Mixing Time | 20-25 | 10 (+11 quick) | 2026-02-17 to 2026-02-18 |
| v1.4 Spectral Gap Refinement | 26 | 2 (+1 quick) | 2026-02-19 to 2026-02-20 |
| v1.5 Krylov Gap Estimation | 27-32 | 12 (+3 quick) | 2026-02-20 to 2026-02-25 |
| v2.0 Restructure | 33-38 | 5/~24 | 2026-02-25 to ... |

## Accumulated Context

### Decisions

All decisions logged in PROJECT.md Key Decisions table.

- [33-01] Field ordering: singletons first, then system/physics/grid/thermalize-specific params
- [33-01] with_coherent derived from construction type via trait, not stored as field
- [33-01] Outer constructor infers S,D,C from singleton args and T from beta
- [33-02] GNS coherent validation removed from validate_config! -- type system enforces via trait
- [33-02] with_coherent kept in serialized Dict for backward compat, computed via trait on write
- [33-02] _reconstruct_config builds Config with explicit sim/construction singletons
- [33-03] construct_lindbladian constrained to Config{Lindbladian} (only caller is run_lindbladian)
- [33-03] _accumulate_R_total! and _accumulate_jump_sandwich! use Config{<:Any,D} for both sim types
- [33-03] _thermalize_to_liouv_config deleted -- unified Config eliminates conversion need
- [33-04] Gitignored playground files migrated on disk but not committed
- [33-04] Pre-existing test_diagnostics.jl failures (7) accepted -- numerical threshold, not migration-related
- [Phase 38]: Replace ill-conditioned Trotter/Bohr ratio test with absolute threshold checks for KMS regime
- [Phase quick-39]: Use kron(conj(J),J) for L*rho*L' Lindblad convention; simplify adjoint G to swap(G_left,G_right) for all domains
- [34-01] domain_prefactor returns domain-only scalar; callers compose with gamma_norm_factor/jump_weight_scaling
- [34-01] Old JumpOp-based oft! deleted; unified oft!(out, eigenbasis, bohr_freqs, energy, inv_4sigma2) is single source
- [34-01] inv_4sigma2 computed at call site, not stored in precomputed_data

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

### Quick Tasks Completed

| # | Description | Date | Commit | Directory |
|---|-------------|------|--------|-----------|
| 38 | Fix failing diagnostics and filename tests by passing construction=KMS() where tests assume exact Gibbs fixed point | 2026-02-26 | d1a2577 | [38-fix-failing-diagnostics-and-filename-tes](./quick/38-fix-failing-diagnostics-and-filename-tes/) |
| 39 | Fix Krylov Lindblad operator convention to L*rho*L' matching Thermalize path | 2026-02-26 | b4b1cec | [39-fix-krylov-lindblad-operator-convention-](./quick/39-fix-krylov-lindblad-operator-convention-/) |

## Session Continuity

Last session: 2026-02-26
Stopped at: Completed 34-01-PLAN.md (domain_prefactor + OFT unification)
Resume file: None
