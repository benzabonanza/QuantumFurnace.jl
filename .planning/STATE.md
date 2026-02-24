# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-20)

**Core value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers
**Current focus:** v1.5 Krylov Gap Estimation -- Phase 30 (Cross-validation)

## Current Position

Phase: 30 of 31 (Cross-validation)
Plan: 0 of ? complete
Status: Ready
Last activity: 2026-02-24 -- Completed 29-02 (Eigensolver Integration: krylov eigsolve tests)

Progress: [██████░░░░] 60% (v1.5 phases 27-31, phases 27-29 complete)

## Performance Metrics

**Velocity:**
- Total plans completed: 78 (v1.0: 10, v1.1: 16, quick: 19, v1.2: 12, cleanup: 3, v1.3: 10, v1.4: 2, v1.5: 6)

**By Milestone:**

| Milestone | Phases | Plans | Timeline |
|-----------|--------|-------|----------|
| v1.0 Trajectories | 1-5 | 10 | 2026-02-13 to 2026-02-14 |
| v1.1 Reduce | 6-11 | 16 (+5 quick) | 2026-02-15 |
| v1.2 Multi-threading | 12-19 | 15 (+3 quick) | 2026-02-15 to 2026-02-16 |
| v1.3 Mixing Time | 20-25 | 10 (+11 quick) | 2026-02-17 to 2026-02-18 |
| v1.4 Spectral Gap Refinement | 26 | 2 (+1 quick) | 2026-02-19 to 2026-02-20 |
| v1.5 Krylov Gap Estimation | 27-31 | 6 | 2026-02-20 to -- |

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
- Kron-derived 2op dissipator formula (B_dag'*rho*A' sandwich) instead of physics convention (v1.5, 28-02)
- Dedicated _accumulate_adjoint_dissipator_2op! with all-N BLAS flags for real-valued BohrDomain (v1.5, 28-02)
- kron(A,B)vec(X)=vec(B*X*A^T) convention consistently applied to all matvec terms (quick-35)
- Anticommutator uses (L'L)^T, coherent uses i[B^T, rho] -- matches dense vectorization exactly (quick-35)
- KrylovKit Arnoldi (not Lanczos) for non-Hermitian Lindbladian eigsolve (v1.5, 29-01)
- copy(vec(ws.rho_out)) in KrylovKit closure to prevent aliasing of Krylov basis vectors (v1.5, 29-01)
- Exact linear lambda_L = (mu-1)/delta for channel-to-Lindbladian eigenvalue conversion (v1.5, 29-01)
- 50% krylovdim increase per retry (30->45->68->102) for convergence recovery (v1.5, 29-01)
- Dense extract_leading_eigendata as ground truth for Krylov eigsolve accuracy tests (v1.5, 29-02)
- rtol=1e-6 for Lindbladian path, rtol=1e-3 for channel path O(delta^2) error (v1.5, 29-02)
- Back-computation mu = 1 + delta*lambda_L for channel eigenvalue conversion verification (v1.5, 29-02)

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

### Quick Tasks Completed

| # | Description | Date | Commit | Directory |
|---|-------------|------|--------|-----------|
| 35 | Fix Krylov matvec dissipator convention to match dense kron-based convention for complex jump operators | 2026-02-20 | ca813c9 | [35-fix-krylov-matvec-dissipator-convention-](./quick/35-fix-krylov-matvec-dissipator-convention-/) |

## Session Continuity

Last session: 2026-02-24
Stopped at: Completed 29-02-PLAN.md (Eigensolver Integration: krylov eigsolve tests) -- Phase 29 complete, ready for Phase 30
Resume file: None
