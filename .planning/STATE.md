# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-20)

**Core value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers
**Current focus:** Phase 32 (Krylov Simulator Speedup) COMPLETE

## Current Position

Phase: 32 of 32 (Some speedup for the Krylov simulator)
Plan: 2 of 2 complete
Status: Phase Complete
Last activity: 2026-02-25 -- Plan 32-02 complete (dead code removal of legacy Euler + dissipator functions)

Progress: [██████████] 100% (Phase 32: plan 2 of 2 done)

## Performance Metrics

**Velocity:**
- Total plans completed: 85 (v1.0: 10, v1.1: 16, quick: 21, v1.2: 12, cleanup: 3, v1.3: 10, v1.4: 2, v1.5: 11)

**By Milestone:**

| Milestone | Phases | Plans | Timeline |
|-----------|--------|-------|----------|
| v1.0 Trajectories | 1-5 | 10 | 2026-02-13 to 2026-02-14 |
| v1.1 Reduce | 6-11 | 16 (+5 quick) | 2026-02-15 |
| v1.2 Multi-threading | 12-19 | 15 (+3 quick) | 2026-02-15 to 2026-02-16 |
| v1.3 Mixing Time | 20-25 | 10 (+11 quick) | 2026-02-17 to 2026-02-18 |
| v1.4 Spectral Gap Refinement | 26 | 2 (+1 quick) | 2026-02-19 to 2026-02-20 |
| v1.5 Krylov Gap Estimation | 27-31 | 9 | 2026-02-20 to 2026-02-24 |

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
- First-order approximation lambda_L = (mu-1)/delta for channel-to-Lindbladian eigenvalue conversion with O(delta) error (v1.5, 29-01; corrected quick-37)
- 50% krylovdim increase per retry (30->45->68->102) for convergence recovery (v1.5, 29-01)
- Dense extract_leading_eigendata as ground truth for Krylov eigsolve accuracy tests (v1.5, 29-02)
- rtol=1e-6 for Lindbladian path, rtol=1e-3 for channel path O(delta^2) error (v1.5, 29-02)
- Back-computation mu = 1 + delta*lambda_L for channel eigenvalue conversion verification (v1.5, 29-02)
- Physics convention (L*rho*L') for channel sandwiches, kron convention for Lindbladian matvec (quick-36)
- Precompute R_total, K0, U_residual, U_coherent at workspace construction for faithful Chen channel (quick-36)
- alpha_chen = 1-sqrt(1-delta) naming to avoid shadowing BohrDomain alpha function (quick-36)
- Relaxed channel eigsolve rtol to 2e-3 due to faithful channel O(delta^2) eigenvalue mapping error (quick-36)
- atol=1e-8 for n=4 Krylov vs dense cross-validation (KrylovKit tol=1e-10 provides margin) (v1.5, 30-01)
- L-vs-E convergence order >= 1.5 hard assertion with deltas [0.1, 0.01, 0.001] (v1.5, 30-01)
- compare_krylov_dense shared helper pattern for all domain/balance cross-validation tests (v1.5, 30-01)
- atol=1e-6 for n=6 cross-validation (looser than n=4 due to larger Krylov subspace approximation error at dim^2=4096) (v1.5, 30-02)
- No GNS at n=6 per locked decision (n=4 sufficient for balance-type correctness) (v1.5, 30-02)
- Dedicated n6_trotter and n6_trotter_sys for TrotterDomain (separate eigenbasis from Hamiltonian) (v1.5, 30-02)
- krylovdim=50 sufficient at n=6 EnergyDomain (1e-15 error vs dense), used for n=6 and n=7 (v1.5, 31-01)
- @elapsed for timing, separate @allocated call for allocation measurement (not BenchmarkTools) (v1.5, 31-01)
- BenchmarkRow struct for clean data pipeline from measurement to report generation (v1.5, 31-01)
- Scaling assertion [3.5,12.0] not [3.5,4.5]: O(8^n) BLAS gemm per matvec gives b~10 for EnergyDomain (v1.5, 31-02)
- Log-space linear regression for power-law fit across 4+ orders of magnitude timing data (v1.5, 31-02)
- EnergyDomain n=10 feasible on cluster (~111h, ~1.3 GB); n=12 infeasible (~12000h) (v1.5, 31-02)
- TrotterDomain n=10+ infeasible due to ~34 GB NUFFT prefactors; needs on-the-fly NUFFT (v1.5, 31-02)
- Store 4 separate G matrices (G_left, G_right, G_left_adj, G_right_adj) for correct BohrDomain adjoint (32-01)
- BohrDomain adjoint uses conj(R_total) not R_total^T because R_total is non-Hermitian for Bohr (32-01)
- Pre-transpose G matrices at construction so hot path uses gemm!('N','N',...) for zero-allocation (32-01)

### Roadmap Evolution

- Phase 32 added: Some speedup for the Krylov simulator

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
| 36 | Replace Euler apply_delta_channel! with faithful Chen CPTP channel (Eq. 3.2) | 2026-02-24 | f761e79 | [36-fix-apply-delta-channel-to-use-faithful-](./quick/36-fix-apply-delta-channel-to-use-faithful-/) |
| 37 | Fix failing XVAL-03 Krylov cross-validation tests (convergence order threshold 1.5->0.9) | 2026-02-25 | 1036b35 | [37-fix-failing-krylov-cross-validation-test](./quick/37-fix-failing-krylov-cross-validation-test/) |

## Session Continuity

Last session: 2026-02-25
Stopped at: Completed quick task 37 -- fix failing XVAL-03 convergence tests
Resume file: None
