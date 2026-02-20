# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-20)

**Core value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers
**Current focus:** v1.5 Krylov Gap Estimation -- Phase 28 (Other-Domain Dispatch)

## Current Position

Phase: 28 of 31 (Other-Domain Dispatch)
Plan: 1 of 2 complete
Status: In Progress
Last activity: 2026-02-20 -- Completed 28-01 (Time/Trotter Domain Matvec)

Progress: [████░░░░░░] 40% (v1.5 phases 27-31, phase 28 plan 1/2 complete)

## Performance Metrics

**Velocity:**
- Total plans completed: 74 (v1.0: 10, v1.1: 16, quick: 18, v1.2: 12, cleanup: 3, v1.3: 10, v1.4: 2, v1.5: 3)

**By Milestone:**

| Milestone | Phases | Plans | Timeline |
|-----------|--------|-------|----------|
| v1.0 Trajectories | 1-5 | 10 | 2026-02-13 to 2026-02-14 |
| v1.1 Reduce | 6-11 | 16 (+5 quick) | 2026-02-15 |
| v1.2 Multi-threading | 12-19 | 15 (+3 quick) | 2026-02-15 to 2026-02-16 |
| v1.3 Mixing Time | 20-25 | 10 (+11 quick) | 2026-02-17 to 2026-02-18 |
| v1.4 Spectral Gap Refinement | 26 | 2 (+1 quick) | 2026-02-19 to 2026-02-20 |
| v1.5 Krylov Gap Estimation | 27-31 | 3 | 2026-02-20 to -- |

## Accumulated Context

### Decisions

All decisions logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- Dense eigen() for exact diagnostics (v1.4): provides the reference for Krylov cross-validation
- Canonical 6-observable set (v1.4): available for trajectory-based comparison if needed
- Biorthogonal overlap formula (v1.4): left+right eigenvectors for diagnostics
- KrylovWorkspace{T,PD} dual type params for zero-overhead NamedTuple access (v1.5, 27-01)
- BLAS.gemm!/axpy! instead of mul!/broadcast for zero-allocation hot path (v1.5, 27-02)
- Concrete-typed jump_eigenbases in KrylovWorkspace to avoid JumpOp abstract field boxing (v1.5, 27-02)
- Separate _accumulate_adjoint_dissipator! preserving {L'L, rho} anticommutator (v1.5, 27-02)
- _prefactor_view + broadcast multiply for zero-alloc NUFFT OFT in Krylov hot path (v1.5, 28-01)
- _measure_matvec_allocs helper to avoid @testset soft-scope variable boxing (v1.5, 28-01)

### Deferred from v1.4

- FIT-01/02/03: Two-exponential fitting with Prony initialization
- RATE-01/02/03/04: Effective rate lambda_eff(t) and automatic window selection
- BOOT-01/02/03: Batched bootstrap uncertainty quantification
- RICH-01/02/03: Richardson extrapolation with monotonicity gate
- VAL-01/02/03/04: Diagnostic dashboard and final validation

### Deferred from v1.5

- BIEIG-01: bieigsolve for left+right eigenvectors simultaneously
- SECTOR-01: Sector-resolved gap computation
- SCALE-01: n=10, n=12 production runs on cluster
- ADAPT-01: Adaptive krylovdim auto-increase on partial convergence

### Pending Todos

None

### Blockers/Concerns

None

## Session Continuity

Last session: 2026-02-20
Stopped at: Completed 28-01-PLAN.md (Time/Trotter Domain Matvec) -- ready for 28-02 (BohrDomain)
Resume file: None
