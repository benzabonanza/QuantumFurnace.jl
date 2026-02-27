# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-25)

**Core value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers
**Current focus:** v2.0 Restructure -- Phase 36 complete, ready for Phase 37: File Organization and Dead Code

## Current Position

Phase: 36 of 38 (API and Results)
Plan: 4/4 complete (Phase 36 COMPLETE)
Status: Completed 36-04: Exports, tests, and simulation scripts
Last activity: 2026-02-27 - Completed 36-04: Round-trip tests for all 4 Result types + simulation script API migration

Progress: [######░░░░] ~42% (v2.0, 13/~24 plans)

## Performance Metrics

**Velocity:**
- Total plans completed: 102 (v1.0: 10, v1.1: 16, quick: 26, v1.2: 12, cleanup: 3, v1.3: 10, v1.4: 2, v1.5: 12, v2.0: 12)

**By Milestone:**

| Milestone | Phases | Plans | Timeline |
|-----------|--------|-------|----------|
| v1.0 Trajectories | 1-5 | 10 | 2026-02-13 to 2026-02-14 |
| v1.1 Reduce | 6-11 | 16 (+5 quick) | 2026-02-15 |
| v1.2 Multi-threading | 12-19 | 15 (+3 quick) | 2026-02-15 to 2026-02-16 |
| v1.3 Mixing Time | 20-25 | 10 (+11 quick) | 2026-02-17 to 2026-02-18 |
| v1.4 Spectral Gap Refinement | 26 | 2 (+1 quick) | 2026-02-19 to 2026-02-20 |
| v1.5 Krylov Gap Estimation | 27-32 | 12 (+3 quick) | 2026-02-20 to 2026-02-25 |
| v2.0 Restructure | 33-38 | 13/~24 | 2026-02-25 to ... |

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
- [34-02] _accumulate_sandwich_adj! as canonical name for L'*rho*L sandwich (symmetric with _accumulate_sandwich!)
- [34-02] _build_cptp_channel returns NamedTuple (; K0, U_residual, alpha); callers destructure what they need
- [34-02] scratch.K0 field now dead in KrausScratch (struct removal deferred to Phase 35)
- [35-01] Workspace{S,D,C,T,SC} uses 5th type parameter SC for concrete scratch typing on hot paths
- [35-01] Union{Nothing, Function} boxing accepted at ~300 bytes/call (MATVEC_ALLOC_BUDGET=512)
- [35-01] _TransitionWrap{F} pattern for function barrier dispatch on transition closures
- [35-01] KrylovWorkspace kept as const alias for backward compatibility
- [35-01] JumpOp[jump] typed vector literal for Julia invariant parameterization
- [35-02] _build_trajectory_workspace factory avoids dispatch conflict with Workspace(Config{Thermalize}) constructor
- [35-02] Per-operator Kraus data as flat Vector{Matrix{CT}} (Rs, K0s, U_residuals, U_Bs) eliminating PerOperatorKraus struct
- [35-02] TrajectoryScratch holds only mutable buffers; _copy_workspace_for_thread shares immutable data
- [35-02] step_along_trajectory! simplified from 4-arg (psi, fw, ws, rng) to 3-arg (psi, ws, rng)
- [36-01] config_kind tag changed from "liouv" to "lindbladian" for new saves; backward compat with "liouv" preserved
- [36-01] _result_to_dict uses multiple dispatch (one method per concrete Result type)
- [36-01] _trajectory_to_dict_new suffix avoids name clash with existing _trajectory_to_dict

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

Last session: 2026-02-27
Stopped at: Completed 36-04-PLAN.md (Exports, tests, simulation scripts). Phase 36 COMPLETE. Ready for Phase 37.
Resume file: None
