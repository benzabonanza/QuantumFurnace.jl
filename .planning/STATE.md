# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-25)

**Core value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers
**Current focus:** v2.0 Restructure -- Phase 38 in progress

## Current Position

Phase: 38 of 38 (Test Cleanup)
Plan: 4/5 complete
Status: Completed 38-02: @info output and threshold rationale for 8 shorter test files (filling gap)
Last activity: 2026-02-28 - Completed 38-02: Added 38 @info lines and threshold rationale to 8 shorter test files (compilation, cptp, dm_detailed_balance, trajectory_fixes, regression, allocation, workspace_independence, threading)

Progress: [######░░░░] ~50% (v2.0, 19/~24 plans)

## Performance Metrics

**Velocity:**
- Total plans completed: 108 (v1.0: 10, v1.1: 16, quick: 26, v1.2: 12, cleanup: 3, v1.3: 10, v1.4: 2, v1.5: 12, v2.0: 18)

**By Milestone:**

| Milestone | Phases | Plans | Timeline |
|-----------|--------|-------|----------|
| v1.0 Trajectories | 1-5 | 10 | 2026-02-13 to 2026-02-14 |
| v1.1 Reduce | 6-11 | 16 (+5 quick) | 2026-02-15 |
| v1.2 Multi-threading | 12-19 | 15 (+3 quick) | 2026-02-15 to 2026-02-16 |
| v1.3 Mixing Time | 20-25 | 10 (+11 quick) | 2026-02-17 to 2026-02-18 |
| v1.4 Spectral Gap Refinement | 26 | 2 (+1 quick) | 2026-02-19 to 2026-02-20 |
| v1.5 Krylov Gap Estimation | 27-32 | 12 (+3 quick) | 2026-02-20 to 2026-02-25 |
| v2.0 Restructure | 33-38 | 18/~24 | 2026-02-25 to ... |

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
- [37-01] src/staging/ directory for dormant gap estimation code excluded from module includes and test suite
- [37-01] using/include cleanup for deleted/moved files done in Task 2 as Rule 3 blocking fix (module would not load otherwise)
- [37-02] Export list organized into Lindbladian/Thermalize/Krylov/Trajectory/Diagnostics/Common sections with Common at bottom
- [37-02] Dormant exports kept as # STAGING: commented-out block (not deleted)
- [37-02] Dead nprocs()/SharedArrays removed -- multi-process mode fully dead after Distributed removal
- [38-01] Unified make_config(sim, domain; kwargs...) uses Config(; ...) with leading semicolon for keyword-only NamedTuple splatting
- [38-01] Default construction=KMS() in make_config; GNS callers pass construction=GNS() explicitly
- [38-01] Trajectory validation gated at runtests.jl level (QUANTUMFURNACE_FULL_TESTS) to save compilation time
- [38-03] All println() in test_dm_scaling.jl replaced with structured @info keyword args
- [38-03] test_observable_trajectories.jl and test_results.jl get zero @info (all exact/structural checks)
- [38-03] Per-iteration @info for ratio tests; summary @info for CPTP loop tests
- [38-04] Classify tests as structural vs numerical: isapprox/abs/norm -> @info; isa/==/haskey -> skip
- [38-04] Loop-summary @info pattern: track max_err across iterations, emit single @info after loop
- [38-04] All existing thresholds confirmed appropriate with O(DIM^n * eps) error analysis -- no changes needed
- [38-02] CPTP threshold 1e-10 kept: algebraic identity error scales as DIM^2 * eps ~ 3e-13, giving ~300x margin
- [38-02] Allocation @info includes both allocs_bytes and threshold for direct comparison in test output
- [38-02] PSD guard eigenvalue threshold -1e-14 documented for FP rounding in Hermitian eigvals

### Key Constraints for v2.0

- Keep time_oft!/trotter_oft! as test utilities (don't delete)
- Keep OFTCaches for testing reference
- R/K0/U_residual: per-jump (R^a) for DM/Trajectory, summed (R_total) for Krylov -- don't mix
- Config is `Config{S,D,C,T}` where C = Construction (KMS, GNS, DLL future)
- Stay flat in src/ for active code (src/staging/ is for dormant code only) -- rename files for clarity
- Diagnostics stays as separate module
- SharedArrays removed (was only used behind dead Distributed branches)
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

Last session: 2026-02-28
Stopped at: Completed 38-02-PLAN.md (@info and threshold rationale for 8 shorter test files). 38-05 remaining.
Resume file: None
