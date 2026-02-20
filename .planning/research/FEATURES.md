# Feature Landscape: Krylov-Based Lindbladian Spectral Gap Estimation

**Domain:** Matrix-free Krylov eigensolving for Lindbladian spectral analysis
**Researched:** 2026-02-20
**Confidence:** HIGH (core algorithms well-understood, codebase deeply analyzed, KrylovKit.jl API verified via official docs)

## Scope

This research covers features for a new milestone adding **matrix-free Krylov-based spectral gap estimation** to QuantumFurnace.jl. The goal is to compute leading Lindbladian eigenvalues (steady state at lambda=0, gap mode at lambda_1) for systems up to 12 qubits without constructing the dense d^2 x d^2 superoperator matrix.

**Key constraint:** The Lindbladian superoperator at n=12 qubits would be a 16,777,216 x 16,777,216 matrix (~2 petabytes dense). Matrix-free Krylov methods reduce this to storing ~k vectors of size d^2 = 16,777,216, where k is the Krylov dimension (typically 30-100).

---

## Existing Features (Already Built)

These are NOT part of this milestone but form the foundation:

| Feature | Location | Status | Relevance |
|---------|----------|--------|-----------|
| `construct_lindbladian()` dense | `furnace.jl` | Shipped | Returns dense d^2 x d^2 matrix; n<=6 only |
| `run_lindbladian()` with Arpack | `furnace.jl` | Shipped | Shift-invert via `eigs(L, sigma=eps)` for 2 eigenvalues |
| `extract_leading_eigendata()` | `diagnostics.jl` | Shipped | Dense `eigen(L)` for n_modes=20 with left+right eigenvectors |
| `run_exact_diagnostics()` bundle | `diagnostics.jl` | Shipped | Full diagnostic suite: eigendata, fixed point, defect, overlap, Sz |
| `estimate_spectral_gap()` trajectory | `gap_estimation.jl` | Shipped | Exponential fitting from observable time series |
| `_jump_contribution!()` Kraus step | `jump_workers.jl` | Shipped | Per-jump CPTP channel application to density matrix |
| `_vectorize_liouv_diss_and_add!()` | `qi_tools.jl` | Shipped | Vectorized Lindbladian accumulation (dense target) |
| Commented-out `LinearMap` approach | `linearmaps_liouv.jl` | Abandoned | Early attempt; marked "slow" due to calling time_oft repeatedly |
| `TrajectoryFramework` with `PerOperatorKraus` | `trajectories.jl` | Shipped | Precomputed per-jump Kraus data for fast channel application |
| `LinearMaps` package | `Project.toml` | Available | Already a dependency; provides `LinearMap{ComplexF64}` type |

---

## Table Stakes

Features that ARE the milestone. Without these, no Krylov-based spectral gap computation is possible.

| Feature | Why Expected | Complexity | Dependencies |
|---------|--------------|------------|--------------|
| **Matrix-free `apply_lindbladian!(out, rho)` function** | The core linear map that KrylovKit needs. Computes L(rho) without forming the d^2 x d^2 matrix. Must reuse existing `_jump_contribution!` and `_vectorize_liouv_diss_and_add!` infrastructure but operating on a single vectorized rho input. | HIGH | Existing jump infrastructure, precomputed data, coherent terms |
| **KrylovKit.jl `eigsolve` integration** | The actual eigenvalue computation. Must target eigenvalues near zero (the spectral gap) for a non-Hermitian operator. | MEDIUM | Matrix-free linear map, eigenvalue targeting strategy |
| **Eigenvalue targeting strategy for gap** | Krylov methods find extremal eigenvalues, not interior ones. For the Lindbladian, eigenvalues near zero (gap) are NOT extremal -- they are interior. Requires either (a) spectral transformation or (b) the Arnoldi-Lindblad time-evolution approach. | HIGH | Deep understanding of Lindbladian spectrum structure |
| **Cross-validation against dense `eigen()` for n<=6** | Must agree with `extract_leading_eigendata()` results to machine precision (or known tolerance). This is the correctness gate. | MEDIUM | Both dense and Krylov paths working |
| **Memory-bounded operation for n<=12** | At n=12, d^2 = 16M complex entries per vector = 256 MB per Krylov vector. With k=50 Krylov vectors: ~12.8 GB. Must stay within ~32 GB RAM. | MEDIUM | Careful Krylov dimension tuning, pre-allocated buffers |

---

## Feature Details

### 1. Matrix-Free `apply_lindbladian!` (LINEAR MAP)

**What it computes:**
Given a vectorized density matrix `v = vec(rho)`, compute `w = vec(L(rho))` where:
```
L(rho) = sum_k (A_k rho A_k^dagger - 1/2 {A_k^dagger A_k, rho}) + coherent_term
```

**Two formulation variants needed:**

**(a) Continuous-time Lindbladian `L(rho)`:**
Computes the Lindbladian action. Eigenvalues of L have Re(lambda) <= 0. The spectral gap is `|Re(lambda_2)|` where lambda_2 is the eigenvalue with second-smallest `|Re(lambda)|`.

**(b) Discrete-time CPTP channel `E(rho) = (I + delta*L)(rho)`:**
The channel that the trajectory simulation actually uses. Eigenvalues mu of E satisfy |mu| <= 1, with mu=1 for the fixed point. The gap is `1 - |mu_2|`. This is the formulation faithful to Chen's algorithm.

**Why both matter:** The continuous-time Lindbladian is the theoretical object; the discrete-time channel is what the trajectory simulation implements. Cross-checking both against dense diagonalization validates the full pipeline.

**Implementation approach:**
The existing `_jump_contribution!` for Liouvillian construction (lines ~1-135 of `jump_workers.jl`) accumulates into a dense `L_target` matrix using `_vectorize_liouv_diss_and_add!()` which does kronecker products. For matrix-free operation, we need to instead compute `L(rho)` directly as operator-level multiplications:

```julia
# For each jump operator A_k with rate gamma_k:
# L_k(rho) = gamma_k * (A_k * rho * A_k' - 0.5 * (A_k' * A_k * rho + rho * A_k' * A_k))
# Plus coherent: L_coh(rho) = -i * [B_total, rho]
```

This uses dim x dim matrix multiplications (cost O(d^3)) per jump per energy label, compared to the d^2 x d^2 dense matrix (cost O(d^4) to build, O(d^4) per matvec). The matrix-free approach is O(d^3 * n_jumps * n_energies) per application.

**Critical subtlety:** The existing `_jump_contribution!` for Liouvillian construction operates differently from the Kraus step `_jump_contribution!` for thermalization. The Liouvillian version vectorizes into a superoperator matrix. The matrix-free version must implement the same physics as the Liouvillian version but without vectorization -- directly computing `L(rho)` as a dim x dim matrix.

**For the CPTP channel variant:** Reuse the existing `_jump_contribution!` for `AbstractThermalizeConfig` which already applies E(rho) to a density matrix in-place (via `_finalize_kraus_step!`). The Krylov interface just wraps this: `f(v) = vec(E(reshape(v, dim, dim)))`.

**Complexity:** HIGH. Not because any single computation is hard, but because:
1. Must handle all four domains (Bohr, Energy, Time, Trotter) and both KMS/GNS balance types
2. Must precompute per-jump data (OFT, NUFFT prefactors, coherent unitaries) once, not per Krylov iteration
3. Must manage workspace allocation carefully -- KrylovKit will call the linear map O(k * maxiter) times
4. The existing `linearmaps_liouv.jl` was abandoned as "slow" -- the new implementation must be faster by precomputing everything upfront

**Dependencies:** All existing jump worker infrastructure, precomputed data framework, coherent term computation.

---

### 2. Eigenvalue Targeting Strategy (SPECTRAL TRANSFORM)

This is the hardest algorithmic decision in the milestone.

**The problem:** KrylovKit.jl's `eigsolve` with Arnoldi targets extremal eigenvalues -- eigenvalues on the periphery of the spectrum. For a Lindbladian:
- The spectrum spans from Re(lambda) = 0 (steady state) down to Re(lambda) ~ -||L|| (fastest-decaying modes)
- The spectral gap eigenvalue lambda_2 is near zero: Re(lambda_2) ~ -0.01 to -0.1 for typical systems
- The most negative eigenvalue could be Re(lambda_min) ~ -100 or worse
- So the gap eigenvalue is INTERIOR, not extremal

**Three viable strategies, in order of preference:**

**(a) Use `:LR` (largest real part) targeting -- RECOMMENDED**
All Lindbladian eigenvalues have Re(lambda) <= 0. The steady state has Re(lambda_1) = 0. The gap mode has the next-largest real part. So `:LR` (largest real part) targets exactly the right part of the spectrum: the eigenvalues closest to zero from below.

Request `howmany=5` to get the steady state + a few gap candidates. The steady state is identified by Re(lambda) closest to zero; the gap mode is the next one.

**Why this works:** The steady state is guaranteed to be the rightmost eigenvalue (Re=0), and the gap mode is next. This makes it an EXTREMAL eigenvalue problem in the real-part sense. KrylovKit's Arnoldi handles this natively.

**Confidence:** HIGH. This is the standard approach for Lindbladian spectral gap computation. The Lindbladian spectrum structure (all Re(lambda) <= 0 with lambda=0 at the boundary) makes the gap eigenvalue extremal when sorted by real part.

**(b) Arnoldi-Lindblad time-evolution approach (alternative)**
Instead of applying L directly, apply the dynamical map E_dt = exp(L*dt) for a chosen timestep dt. The eigenvalues of E_dt are mu_k = exp(lambda_k * dt). The largest-magnitude eigenvalue is mu=1 (steady state), and the gap mode has |mu_2| close to 1.

Use `eigsolve(E_dt_map, x0, howmany, :LM)` to find largest-magnitude eigenvalues, then recover lambda_k = log(mu_k)/dt.

**Advantages:** Naturally targets long-lived modes. The time-evolution map is the physical operation.
**Disadvantages:** Requires choosing dt (too small: all mu_k near 1, hard to distinguish; too large: gap mode has decayed to noise). More expensive per iteration (time evolution rather than single L application).

**Confidence:** MEDIUM. The Arnoldi-Lindblad paper (Huybrechts & Roscilde, Quantum 2022) demonstrates this works, but the implementation complexity is higher and the dt selection is system-dependent.

**(c) Manual shift-invert via GMRES (fallback only)**
Define the shifted operator (L - sigma*I)^{-1} and find its largest eigenvalues (which correspond to eigenvalues of L closest to sigma). Apply the inverse via iterative GMRES solve inside each Krylov iteration.

**Disadvantages:** Nested iteration (GMRES inside Arnoldi). Expensive. Convergence of the outer iteration depends on GMRES accuracy. Much more complex to implement.
**When needed:** Only if strategy (a) fails to converge, which is unlikely given the Lindbladian spectrum structure.

**Recommendation:** Start with strategy (a) `:LR` targeting. Fall back to (b) Arnoldi-Lindblad if convergence is poor for larger systems. Do not implement (c) unless (a) and (b) both fail.

**Complexity:** LOW for strategy (a), HIGH for strategy (b), VERY HIGH for strategy (c).

**Dependencies:** Matrix-free linear map.

---

### 3. KrylovKit.jl eigsolve Integration (KRYLOV INTERFACE)

**Core API call:**
```julia
using KrylovKit

vals, vecs, info = eigsolve(
    apply_L,            # function: Vector{ComplexF64} -> Vector{ComplexF64}
    x0,                 # starting vector: vec(rho_0), size d^2
    howmany,            # number of eigenvalues (5 for gap analysis)
    :LR;                # largest real part
    issymmetric=false,  # Lindbladian is NOT Hermitian
    krylovdim=50,       # Krylov subspace dimension
    maxiter=200,        # maximum restarts
    tol=1e-10,          # convergence tolerance
    verbosity=1,        # print convergence info
)
```

**Starting vector choice:**
Use `x0 = vec(rho_gibbs) + eps * vec(random_perturbation)` where rho_gibbs is the Gibbs state. This seeds the Krylov subspace with a vector rich in both the steady-state component (lambda=0) and the gap mode.

Alternative: `x0 = vec(I/d)` (maximally mixed state). Has equal weight in all symmetry sectors, so can find the gap mode even if it lives in a symmetry sector that rho_gibbs does not overlap with.

**Return value processing:**
1. Sort returned eigenvalues by `|Re(lambda)|`
2. Identify steady state: the eigenvalue with smallest `|Re(lambda)|`
3. Identify gap mode: next eigenvalue
4. Extract spectral gap: `|Re(lambda_gap)|`
5. Check convergence: `info.converged >= howmany`

**Eigenvector extraction (n<=6 only):**
Reshape each eigenvector `vecs[k]` to a dim x dim matrix. For n<=6, also compute left eigenvectors via `bieigsolve` for biorthogonal overlap analysis. For n>6, eigenvectors may be too expensive to store/analyze.

**Krylov dimension tuning:**
- n=4 (d^2=256): krylovdim=30, maxiter=100 (fast, for testing)
- n=6 (d^2=4096): krylovdim=50, maxiter=200 (moderate)
- n=8 (d^2=65536): krylovdim=60, maxiter=300 (main target)
- n=10 (d^2=1M): krylovdim=80, maxiter=500 (aggressive)
- n=12 (d^2=16M): krylovdim=100, maxiter=1000 (memory-limited)

**Tolerance management:**
- For cross-validation against dense eigen(): tol=1e-10 (tight)
- For production gap estimation at n>6: tol=1e-6 (relaxed, faster convergence)
- User-configurable with sensible defaults

**Complexity:** MEDIUM. The integration itself is straightforward once the linear map exists. The complexity is in parameter tuning and convergence verification.

**Dependencies:** Matrix-free linear map, KrylovKit.jl (new dependency).

---

### 4. Cross-Validation Against Dense eigen() (VALIDATION)

**What to validate:**
At n=4 and n=6, compare Krylov results against existing `extract_leading_eigendata()`:
1. Spectral gap: `|gap_krylov - gap_dense| / gap_dense < tol`
2. Steady-state eigenvalue: `|Re(lambda_0_krylov)| < tol`
3. Gap mode eigenvector overlap: `|<v_krylov | v_dense>| > 1 - tol`
4. Fixed point trace distance: `trace_distance(rho_krylov, rho_dense) < tol`

**Test matrix:**

| n | Domain | Balance | Expected gap | Validation target |
|---|--------|---------|--------------|-------------------|
| 4 | Bohr | KMS | ~0.011 | 1e-8 agreement |
| 4 | Trotter | KMS | ~0.010 | 1e-6 agreement (Trotter error) |
| 6 | Bohr | KMS | ~0.003 | 1e-8 agreement |
| 6 | Trotter | KMS | ~0.003 | 1e-6 agreement |
| 4 | Bohr | GNS | varies | 1e-8 agreement |

**Cross-validation against trajectory estimate:**
Also compare Krylov gap vs trajectory-based `estimate_spectral_gap()` at n=4. The trajectory estimate has ~1% accuracy (from v1.3 validation), so agreement within ~5% confirms consistency.

**Complexity:** MEDIUM. The tests themselves are straightforward. The setup (constructing Lindbladian both ways, running both solvers) takes development time.

**Dependencies:** Both dense diagnostics and Krylov solver working.

---

### 5. Memory-Bounded Operation (RESOURCE MANAGEMENT)

**Memory scaling formula:**

| n | dim | d^2 | Per-vector (ComplexF64) | k=50 vectors | k=100 vectors |
|---|-----|-----|------------------------|--------------|---------------|
| 4 | 16 | 256 | 4 KB | 200 KB | 400 KB |
| 6 | 64 | 4,096 | 64 KB | 3.2 MB | 6.4 MB |
| 8 | 256 | 65,536 | 1 MB | 50 MB | 100 MB |
| 10 | 1,024 | 1,048,576 | 16 MB | 800 MB | 1.6 GB |
| 12 | 4,096 | 16,777,216 | 256 MB | 12.8 GB | 25.6 GB |

**Additional memory per linear map application:**
- Per-jump workspace: ~3 dim x dim matrices = 3 * dim^2 * 16 bytes
- At n=12: 3 * 4096^2 * 16 = ~768 MB per jump (but reused across iterations)
- Coherent term: 1 dim x dim matrix
- Total workspace overhead at n=12: ~1 GB

**Total memory budget at n=12:** ~14 GB (Krylov vectors) + ~1 GB (workspace) = ~15 GB

**Pre-allocation strategy:**
All workspace matrices must be allocated once before the Krylov iteration begins. KrylovKit calls the linear map function hundreds of times -- any allocation inside the function multiplied by iteration count would be catastrophic.

**Memory profiling feature:**
Before running eigsolve, compute and report estimated memory:
```julia
function estimate_krylov_memory(n_qubits, krylov_dim, n_jumps)
    d = 2^n_qubits
    d2 = d^2
    krylov_bytes = krylov_dim * d2 * sizeof(ComplexF64)
    workspace_bytes = (3 * n_jumps + 2) * d^2 * sizeof(ComplexF64)
    total = krylov_bytes + workspace_bytes
    return (krylov_gb=total/1e9, krylov_vectors_gb=krylov_bytes/1e9, workspace_gb=workspace_bytes/1e9)
end
```

**Complexity:** LOW for estimation, MEDIUM for pre-allocation. The main effort is ensuring zero allocations in the hot path.

**Dependencies:** Matrix-free linear map with pre-allocated workspace.

---

## Differentiators

Features that strengthen the milestone beyond the minimum viable Krylov solver.

| Feature | Value Proposition | Complexity | Dependencies |
|---------|-------------------|------------|--------------|
| **Both L and E formulations** | Cross-check continuous vs discrete: gap from L should equal `-log(lambda_2(E))/delta` up to O(delta^2). Validates the Trotter discretization. | MEDIUM | Matrix-free map for both L and E |
| **Adaptive Krylov dimension** | Auto-increase krylovdim if convergence stalls. Start at 30, double if `info.converged < howmany` after maxiter. Prevents user from having to guess. | LOW | KrylovKit integration |
| **Timing benchmarks with scaling extrapolation** | Time per matvec at n=4,6,8, fit scaling law t(n) = a * 4^n + b, extrapolate to n=10,12. Report estimated wall-clock time before starting large runs. | LOW | Working solver at multiple sizes |
| **Per-domain apply_lindbladian! specialization** | BohrDomain can exploit the Bohr dictionary for sparse accumulation. EnergyDomain uses OFT. TrotterDomain uses NUFFT prefactors. Each has different cost per application. Specializing maximizes performance. | HIGH | All four domain worker implementations |
| **bieigsolve for left+right eigenvectors** | KrylovKit's `bieigsolve` with `BiArnoldi` returns both left and right eigenvectors simultaneously. Enables biorthogonal overlap analysis without dense eigen(). Would extend existing `compute_overlap_coefficients()` to n>6. | MEDIUM | KrylovKit bieigsolve API |
| **Sector-resolved gap computation** | Apply `apply_lindbladian!` restricted to a symmetry sector (e.g., Delta_Sz=0 subspace). The sector-restricted Lindbladian is smaller, so Krylov converges faster and the gap within that sector is directly obtained. | HIGH | Symmetry sector projection infrastructure |
| **Resource estimation report** | Before running: memory estimate, time estimate, expected number of matvecs. User decides go/no-go. | LOW | Memory formula + timing benchmarks |

---

## Anti-Features

Features to explicitly NOT build in this milestone.

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| **Sparse Lindbladian matrix construction** | Building a sparse d^2 x d^2 matrix defeats the purpose of matrix-free. At n=12, even sparse storage with ~100 nonzeros per row would be ~26 GB. The matrix-free approach is specifically chosen to avoid this. | Matrix-free `apply_lindbladian!` function. |
| **Full spectrum computation (all d^2 eigenvalues)** | At n=8, d^2=65536 eigenvalues. Dense diag is O(d^6). Even Krylov with howmany=1000 is wasteful. We need only the ~5-20 leading eigenvalues. | Krylov with `howmany=5` for gap, `howmany=20` for extended diagnostics. |
| **GPU-accelerated Krylov** | CUDA adds massive complexity (GPU memory management, kernel launches, PCIe transfer overhead). The bottleneck is dim x dim matrix multiplications which BLAS already handles well on CPU. At n=12, dim=4096, so BLAS gemm at 4096x4096 is fast. | CPU BLAS with appropriate thread count. |
| **Custom Krylov implementation** | KrylovKit.jl is mature, well-tested, and handles the Arnoldi restart logic correctly. Reimplementing Arnoldi would be error-prone and slower. | Use KrylovKit.jl eigsolve directly. |
| **Shift-invert with explicit LU factorization** | LU of a d^2 x d^2 matrix at n=12 is impossible (petabytes). Shift-invert requires solving (L-sigma*I)x=b, which for matrix-free L requires an iterative inner solve. Too complex for this milestone. | Use `:LR` targeting which avoids the need for shift-invert entirely. |
| **Floquet (periodically driven) Lindbladian support** | The Arnoldi-Lindblad paper extends to Floquet systems but QuantumFurnace.jl has no periodic driving. Pure overhead. | Stick with time-independent Lindbladian. |
| **Distributed (multi-node) Krylov** | At n=12, ~15 GB fits in a single workstation. Distributed computing adds communication overhead and code complexity without benefit until n>14. | Single-node, multi-threaded BLAS. |
| **Left eigenvectors for n>6** | At n>6, storing left eigenvectors doubles memory. The gap estimate only needs eigenvalues. Overlap analysis (which needs left eigenvectors) is only scientifically needed at n<=6 for cross-validation. | `eigsolve` (right only) for n>6, `bieigsolve` for n<=6. |
| **Trajectory-Krylov hybrid (Arnoldi-Lindblad with trajectory simulator)** | Using the trajectory simulator to approximate exp(L*dt) inside a Krylov iteration would couple stochastic noise into the eigenvalue estimate. The matrix-free L application is deterministic and exact. | Deterministic `apply_lindbladian!` for Krylov; trajectory-based estimation remains a separate, complementary approach. |

---

## Feature Dependencies

```
[Precomputed Jump Data]  ------>  [Matrix-Free apply_lindbladian!]
  (existing infrastructure)              |
                                         |
                              +----------+----------+
                              |                     |
                     [Continuous L(rho)]    [Discrete E(rho)]
                              |                     |
                              +----------+----------+
                                         |
                              [KrylovKit eigsolve Integration]
                                         |
                              +----------+----------+
                              |                     |
                   [Eigenvalue Targeting]    [Memory Estimation]
                   (use :LR for gap)        (pre-flight check)
                              |
                              v
                   [Cross-Validation n<=6]
                   (dense eigen vs Krylov)
                              |
                     +--------+--------+
                     |                 |
              [Timing Benchmarks]  [bieigsolve for
               (n=4,6,8 scaling)    left+right eigvecs]
                     |                 |
                     v                 v
              [Resource Estimate   [Overlap Analysis
               for n=10,12]         at n>6 (optional)]
                     |
                     v
              [Production Runs n=8,10,12]
```

### Dependency Notes

- **Matrix-free `apply_lindbladian!` is the keystone:** Everything depends on it. Build first, validate thoroughly.
- **Continuous L and discrete E are independent formulations:** Can be built in parallel. Both wrap the same jump infrastructure differently.
- **`:LR` targeting removes the need for shift-invert:** This is the critical insight. Lindbladian spectrum structure (all Re(lambda) <= 0) means the gap eigenvalue IS the second-largest real part. No spectral transformation needed.
- **Cross-validation must precede production runs:** Any gap estimate at n>6 is only trustworthy if the Krylov solver is validated against dense results at n<=6.
- **bieigsolve is optional:** Only needed for overlap analysis. The gap estimate itself only needs right eigenvectors.
- **Timing benchmarks enable resource estimation:** Without timing data at n=4,6,8, the extrapolation to n=10,12 is pure guesswork.

---

## MVP Definition

### Launch With (Core Krylov Gap Estimation)

- [ ] **Matrix-free `apply_lindbladian!(out, rho)` for BohrDomain** -- the simplest domain, cleanest implementation
- [ ] **KrylovKit `eigsolve` with `:LR` targeting** -- find steady state + gap mode
- [ ] **KrylovResult struct** -- eigenvalues, gap, convergence info, timing
- [ ] **Cross-validation at n=4 BohrDomain** -- Krylov gap matches dense gap to 1e-8
- [ ] **Memory estimation function** -- pre-flight check before large runs

### Add After Core Validation

- [ ] **TrotterDomain support** -- extend to the domain used in actual quantum algorithm simulations
- [ ] **EnergyDomain and TimeDomain support** -- complete domain coverage
- [ ] **Discrete CPTP channel formulation** -- E(rho) = (I + delta*L)(rho) as alternative linear map
- [ ] **Cross-validation at n=6** -- confirm at the boundary of dense feasibility
- [ ] **n=8 production run** -- first result beyond dense regime
- [ ] **Timing benchmarks and scaling extrapolation** -- predictive resource estimation

### Future Consideration (Subsequent Milestones)

- [ ] **n=10, n=12 production runs** -- may need memory optimization or reduced krylovdim
- [ ] **bieigsolve for biorthogonal overlap analysis at n>6** -- extends diagnostic power
- [ ] **Sector-resolved gap computation** -- leverages symmetry for faster convergence
- [ ] **Adaptive Krylov dimension** -- auto-tuning for unknown spectrum structure
- [ ] **Integration with deferred v1.4 features** -- combine Krylov gap with effective rate plots, bootstrap, Richardson

---

## Feature Prioritization Matrix

| Feature | Scientific Value | Implementation Cost | Priority |
|---------|-----------------|---------------------|----------|
| Matrix-free `apply_lindbladian!` (BohrDomain) | CRITICAL | HIGH | P0 |
| KrylovKit `eigsolve` with `:LR` | CRITICAL | MEDIUM | P0 |
| Cross-validation n=4 BohrDomain | CRITICAL | MEDIUM | P0 |
| Memory estimation function | HIGH | LOW | P0 |
| KrylovResult struct with convergence info | HIGH | LOW | P0 |
| Matrix-free `apply_lindbladian!` (TrotterDomain) | HIGH | MEDIUM | P1 |
| Cross-validation n=6 (both domains) | HIGH | MEDIUM | P1 |
| Discrete CPTP channel formulation | HIGH | MEDIUM | P1 |
| n=8 production run | HIGH | LOW (if P0/P1 done) | P1 |
| EnergyDomain / TimeDomain support | MEDIUM | MEDIUM | P1 |
| Timing benchmarks + scaling extrapolation | MEDIUM | LOW | P2 |
| Resource estimation report | MEDIUM | LOW | P2 |
| bieigsolve for left+right eigenvectors | MEDIUM | MEDIUM | P2 |
| Sector-resolved gap computation | LOW (this milestone) | HIGH | P3 |
| Adaptive Krylov dimension | LOW | LOW | P3 |

**Priority key:**
- P0: Foundation without which nothing works; build first
- P1: Core milestone deliverables; production capability
- P2: Polish and diagnostics; improves confidence and usability
- P3: Only if needed based on results

---

## Phase Structure Implications

Based on feature dependencies, the natural phase ordering is:

**Phase A (Matrix-Free Foundation):**
- `apply_lindbladian!` for BohrDomain (simplest)
- Precomputed workspace allocation
- Unit test: L(rho_gibbs) = 0 (steady state property)
- Unit test: Tr[L(rho)] = 0 (trace preservation)

**Phase B (Krylov Integration):**
- KrylovKit `eigsolve` wrapper with `:LR` targeting
- KrylovResult struct
- Memory estimation
- Cross-validation at n=4 BohrDomain

**Phase C (Domain Extension):**
- TrotterDomain apply_lindbladian!
- EnergyDomain and TimeDomain apply_lindbladian!
- Discrete CPTP channel formulation
- Cross-validation at n=6 all domains

**Phase D (Production):**
- n=8 Krylov runs with timing
- Scaling benchmarks
- Resource estimation for n=10,12
- GNS balance type support

**Phase ordering rationale:**
- A before B: the linear map must exist before the eigensolver can use it
- B before C: validate the eigensolver on the simplest domain before adding complexity
- C before D: ensure all domains work before running expensive production computations
- Each phase has clear deliverables and can be cross-validated independently

---

## Technical Notes for Implementation

### KrylovKit.jl API Details (Verified from Official Docs)

```julia
# Install
using Pkg; Pkg.add("KrylovKit")

# Basic eigsolve
vals, vecs, info = eigsolve(f, x0, howmany, :LR;
    issymmetric=false,
    krylovdim=50,
    maxiter=200,
    tol=1e-10,
    verbosity=1
)

# Check convergence (CRITICAL: no automatic warning on failure)
@assert info.converged >= howmany "Only $(info.converged) of $howmany eigenvalues converged"

# bieigsolve for left+right eigenvectors
vals, (vecs_right, vecs_left), (info_right, info_left) = bieigsolve(
    f, x0, howmany, :LR;
    krylovdim=50,
    maxiter=200,
    tol=1e-10
)
```

**Key KrylovKit facts:**
- Accepts ANY callable `f(x) -> y` as the linear map -- perfect for matrix-free
- Returns eigenvalues as `Vector{ComplexF64}`, eigenvectors as `Vector{Vector{ComplexF64}}`
- `info.converged` gives number of converged eigenvalues -- MUST check this
- `info.normres` gives residual norms for each eigenvalue
- Default `krylovdim=30` and `tol=1e-12` -- may need adjustment
- `:LR` sorts by largest real part -- exactly what we need for Lindbladian gap

### Matrix-Free Linear Map Signature

```julia
struct KrylovLindbladianMap{T, D<:AbstractDomain}
    jumps::Vector{JumpOp{Matrix{Complex{T}}}}
    ham_or_trott::Union{HamHam{T}, TrottTrott}
    config::AbstractLiouvConfig{D,T}
    precomputed_data::Any
    workspace::KrylovWorkspace{T}  # NEW: pre-allocated buffers
end

# Make it callable (KrylovKit interface)
function (map::KrylovLindbladianMap)(v::AbstractVector{ComplexF64})
    rho = reshape(v, map.workspace.dim, map.workspace.dim)
    out = map.workspace.output  # pre-allocated
    fill!(out, 0)
    # ... apply L(rho) into out ...
    return vec(out)
end
```

### Lindbladian Application Without Vectorization

The core computation for each jump operator A_w at frequency w:
```julia
# L_w(rho) = gamma_w * (A_w * rho * A_w' - 0.5 * (A_w' * A_w * rho + rho * A_w' * A_w))
# Cost: 3 matrix multiplications of size dim x dim = O(dim^3) per frequency per jump
mul!(tmp1, A_w, rho)         # A_w * rho
mul!(out, tmp1, A_w', gamma_w, 1.0)  # += gamma_w * A_w * rho * A_w'

mul!(tmp1, A_w', A_w)        # A_w' * A_w (can precompute!)
mul!(out, tmp1, rho, -0.5*gamma_w, 1.0)  # -= 0.5 * gamma_w * (A_w'A_w) * rho
mul!(out, rho, tmp1, -0.5*gamma_w, 1.0)  # -= 0.5 * gamma_w * rho * (A_w'A_w)
```

**Optimization: precompute `A_w' * A_w` for all w.** This is a dim x dim matrix per frequency per jump. At n=8 with 3 jump paulis * 8 sites * ~2048 energy labels: potentially thousands of precomputed matrices. Must evaluate memory vs recomputation tradeoff per system size.

### Starting Vector Selection

```julia
# Option 1: Gibbs state + perturbation (biases toward gap mode in population sector)
x0 = vec(Matrix(gibbs)) + 1e-6 * randn(ComplexF64, d^2)

# Option 2: Maximally mixed + perturbation (all symmetry sectors equally)
x0 = vec(Matrix{ComplexF64}(I(d)/d)) + 1e-6 * randn(ComplexF64, d^2)

# Option 3: Random positive-definite (no bias)
A = randn(ComplexF64, d, d)
rho0 = A * A'; rho0 /= tr(rho0)
x0 = vec(rho0)
```

Recommendation: Option 2 (maximally mixed + perturbation) for robustness. It has nonzero projection onto all symmetry sectors, avoiding the n=6 zero-overlap problem diagnosed in v1.3.

---

## Comparison with Existing Approaches

| Approach | n range | Accuracy | Wall time (n=8) | Memory (n=8) | Status |
|----------|---------|----------|-----------------|--------------|--------|
| Dense `eigen(L)` | n<=6 | Machine precision | N/A (impossible) | N/A | Existing |
| Arpack shift-invert `eigs(L, sigma)` | n<=6 (dense L) | 1e-12 | ~1 min | ~34 GB (dense L) | Existing |
| Trajectory exponential fitting | any n | ~1% | ~5 min (20k traj) | ~100 MB | Existing |
| **Krylov `eigsolve` (this milestone)** | n<=12 | ~1e-8 | ~10 min (est.) | ~50 MB | Planned |

The Krylov approach fills the gap: trajectory fitting is inaccurate but scales well; dense methods are precise but limited to n<=6. Krylov gives near-machine-precision eigenvalues at n=8-12 where dense methods fail.

---

## Sources

### HIGH Confidence (verified via official docs and codebase)
- [KrylovKit.jl eigenvalue problems documentation](https://jutho.github.io/KrylovKit.jl/stable/man/eig/) -- eigsolve API, `:LR` targeting, convergence info
- [KrylovKit.jl available algorithms](https://jutho.github.io/KrylovKit.jl/stable/man/algorithms/) -- Arnoldi, BiArnoldi, Lanczos parameters
- [KrylovKit.jl GitHub repository](https://github.com/Jutho/KrylovKit.jl) -- source code, release history
- QuantumFurnace.jl codebase: `jump_workers.jl`, `furnace.jl`, `diagnostics.jl`, `qi_tools.jl`, `linearmaps_liouv.jl`, `trajectories.jl`

### MEDIUM Confidence (established methods, multiple sources agree)
- [Arnoldi-Lindblad time evolution paper (Huybrechts & Roscilde, Quantum 2022)](https://quantum-journal.org/papers/q-2022-02-10-649/) -- Arnoldi on dynamical map for Lindbladian spectrum
- [Arnoldi-Lindblad GitHub implementation](https://github.com/DHuybrechts/Arnoldi-Lindblad-time-evolution) -- Python/QuTiP reference implementation
- [Tensor Network Framework for Lindbladian Spectra (arXiv:2509.07709)](https://arxiv.org/abs/2509.07709) -- complex-time Krylov methods for Lindbladian eigenvalues
- [Bi-Lanczos for Lindbladian Krylov complexity (JHEP 2023)](https://link.springer.com/article/10.1007/JHEP12(2023)066) -- bi-Lanczos algorithm for non-Hermitian Lindbladian
- [Lindbladian Wikipedia](https://en.wikipedia.org/wiki/Lindbladian) -- spectrum structure: Re(lambda) <= 0, complex conjugate pairs
- [BifurcationKit.jl matrix-free eigensolver discussion](https://discourse.julialang.org/t/matrix-free-eigensolver-in-bifurcationkit-jl/92780) -- shift-invert via GMRES for KrylovKit

### LOW Confidence (need validation during implementation)
- Time-per-matvec estimates for n=8,10,12 -- based on O(d^3) scaling extrapolation, not measured
- Krylov dimension requirements at n>8 -- system-dependent, needs empirical tuning
- Memory estimates at n=12 -- formulaic, not accounting for Julia GC overhead or BLAS workspace

---
*Feature research for: Krylov-Based Lindbladian Spectral Gap Estimation*
*Researched: 2026-02-20*
