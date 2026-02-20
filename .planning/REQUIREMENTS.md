# Requirements: QuantumFurnace.jl

**Defined:** 2026-02-20
**Core Value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers

## v1.5 Requirements

Requirements for Krylov Gap Estimation milestone. Each maps to roadmap phases.

### Matvec — Matrix-Free Lindbladian/Channel Action

- [ ] **MATVEC-01**: Matrix-free `apply_lindbladian!` computes L(rho) directly via dissipator formula for EnergyDomain (KMS + GNS)
- [ ] **MATVEC-02**: Matrix-free `apply_lindbladian!` computes L(rho) for TimeDomain using NUFFT prefactors (KMS + GNS)
- [ ] **MATVEC-03**: Matrix-free `apply_lindbladian!` computes L(rho) for TrotterDomain in Trotter basis (KMS + GNS)
- [ ] **MATVEC-04**: Matrix-free `apply_lindbladian!` computes L(rho) for BohrDomain using Bohr frequency buckets (KMS + GNS)
- [ ] **MATVEC-05**: Coherent term -i[B, rho] included in L(rho) for KMS configs (with_coherent=true)
- [ ] **MATVEC-06**: Matrix-free `apply_delta_channel!` computes E(rho) = (I + delta*L)(rho) for all domains (KMS + GNS)
- [ ] **MATVEC-07**: Adjoint Lindbladian L'(rho) for all domains (swap A <-> A' in dissipator)
- [ ] **MATVEC-08**: KrylovWorkspace pre-allocates all scratch matrices for zero-alloc matvec hot path
- [ ] **MATVEC-09**: Round-trip test validates ||L_dense * vec(rho) - vec(L_matvec(rho))|| < 1e-12 at n=4

### Krylov — Eigensolver Integration

- [ ] **KRYLOV-01**: `krylov_spectral_gap()` API wrapping KrylovKit eigsolve with :LR targeting for Lindbladian
- [ ] **KRYLOV-02**: `krylov_spectral_gap()` also supports :LM targeting for CPTP channel eigenvalues
- [ ] **KRYLOV-03**: KrylovGapResult struct with eigenvalues, spectral_gap, convergence info, matvec count, fixed point
- [ ] **KRYLOV-04**: Pre-flight memory estimation (krylovdim * 4^n * 16 * 1.5 bytes) with warning if exceeds threshold
- [ ] **KRYLOV-05**: Convergence assertion: info.converged >= howmany, with fallback to increased krylovdim

### Validation — Cross-Validation Framework

- [ ] **XVAL-01**: Krylov gap matches dense eigen() gap to < 1e-8 at n=4 for all 4 domains (KMS)
- [ ] **XVAL-02**: Krylov gap matches dense eigen() gap to < 1e-6 at n=6 for all 4 domains (KMS)
- [ ] **XVAL-03**: Krylov L gap vs Krylov E gap consistency: gap_L ~ -log(|lambda_2(E)|)/delta
- [ ] **XVAL-04**: Krylov KMS vs GNS gap comparison at n=4 matches existing results

### Benchmark — Scaling and Resource Estimation

- [ ] **BENCH-01**: Timing benchmarks at n=3,4,5,6 (n=7 if feasible) with 4 BLAS threads
- [ ] **BENCH-02**: Memory usage measurement at each n
- [ ] **BENCH-03**: Scaling model fit (time proportional to 4^n) with extrapolation to n=10,12
- [ ] **BENCH-04**: Per-matvec timing breakdown: BLAS vs precompute vs Krylov overhead

## Future Requirements

### Deferred from v1.5

- **BIEIG-01**: bieigsolve for left+right eigenvectors simultaneously (needed for biorthogonal overlap diagnostics at n>6)
- **SECTOR-01**: Sector-resolved gap computation (symmetry projection infrastructure)
- **SCALE-01**: n=10, n=12 production runs on cluster (requires n=8 empirical parameter tuning)
- **ADAPT-01**: Adaptive krylovdim auto-increase on partial convergence

### Deferred from v1.4

- FIT-01/02/03: Two-exponential fitting with Prony initialization
- RATE-01/02/03/04: Effective rate lambda_eff(t) and automatic window selection
- BOOT-01/02/03: Batched bootstrap uncertainty quantification
- RICH-01/02/03: Richardson extrapolation with monotonicity gate
- VAL-01/02/03/04: Diagnostic dashboard and final validation

## Out of Scope

| Feature | Reason |
|---------|--------|
| GPU acceleration of Krylov matvec | Data transfer overhead dominates many small dim x dim matrix multiplications |
| Arnoldi-Lindblad time-evolution approach | Higher complexity; only needed if :LR targeting fails (unlikely) |
| n>12 qubit support | Memory and time scale exponentially; 12 qubits is the practical limit for full Lindbladian simulation |
| Sparse Lindbladian storage | Full Lindbladian is dense for these constructions; sparsity not exploitable |
| MKL BLAS substitution | OpenBLAS sufficient; MKL optimization deferred to performance tuning |

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| (populated during roadmap creation) | | |

**Coverage:**
- v1.5 requirements: 22 total
- Mapped to phases: 0
- Unmapped: 22

---
*Requirements defined: 2026-02-20*
*Last updated: 2026-02-20 after initial definition*
