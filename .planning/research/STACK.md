# Technology Stack: Matrix-Free Krylov Spectral Gap Estimation

**Project:** QuantumFurnace.jl -- Krylov-based Lindbladian spectral gap at up to 12 qubits
**Researched:** 2026-02-20
**Confidence:** HIGH (primary recommendations verified via official docs and GitHub releases)

---

## Scope

This STACK.md covers ONLY the additions needed for matrix-free Krylov eigensolving of the Lindbladian superoperator at scales where the full dim^2 x dim^2 Liouvillian matrix cannot be stored (n > 6 qubits). It does NOT re-research the existing stack (LinearAlgebra, Arpack, LsqFit, FINUFFT, etc.) which remains unchanged.

The milestone target: compute the spectral gap `|Re(lambda_2)|` of the Lindbladian or the leading eigenvalue of the CPTP channel `E = I + delta*L` for systems up to 12 qubits (dim = 4096, vectorized dim^2 = 16,777,216).

---

## Critical Design Insight: Lindbladian Spectrum Structure

All Lindbladian eigenvalues have `Re(lambda) <= 0`. The steady state has `lambda_1 = 0`, and the spectral gap is `|Re(lambda_2)|` where `lambda_2` is the eigenvalue with the second-largest (least negative) real part. This means:

- **lambda_1 = 0 is the EXTREMAL eigenvalue in the `:LR` (largest real part) direction.**
- **lambda_2 is the SECOND extremal eigenvalue in that same direction.**

Therefore, Krylov methods with `:LR` targeting can find the spectral gap WITHOUT shift-invert. This is the key enabler: no need for matrix factorization, no need for solving linear systems. Standard Arnoldi with `:LR` will converge to the steady state first, then the gap mode second. We need `howmany=2` (or a few more for robustness).

---

## Recommended Stack

### New Production Dependency: KrylovKit.jl

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| **KrylovKit.jl** | 0.10.2 | Matrix-free Krylov eigensolving via Arnoldi for non-Hermitian Lindbladian | Pure Julia, accepts ANY callable as linear map (no wrapping needed), supports `:LR` for Lindbladian spectrum, `bieigsolve` provides left eigenvectors for biorthogonal decomposition, `exponentiate` for matrix exponential action |

**Confidence: HIGH** -- Version verified via [GitHub releases](https://github.com/Jutho/KrylovKit.jl/releases) (v0.10.2, October 11, 2025). Documentation verified via [official docs](https://jutho.github.io/KrylovKit.jl/stable/) (generated October 10, 2025).

**Why KrylovKit.jl over alternatives:**

| Criterion | KrylovKit.jl | ArnoldiMethod.jl | Arpack.jl (existing) |
|-----------|-------------|-------------------|---------------------|
| Matrix-free (callable) | YES -- any callable `f(x)` | YES -- needs `mul!(y,A,x)` | NO for shift-invert (needs explicit matrix for factorization) |
| `:LR` targeting | YES | YES | YES (but `which=:LR` in Arpack is different from shift-invert `sigma`) |
| Left eigenvectors | YES via `bieigsolve`/`BiArnoldi` | NO | NO (Arpack.jl wraps only partial Arpack -- no left eigvec support) |
| `exponentiate(t, A, x)` | YES -- computes `exp(tA)x` matrix-free | NO | NO |
| Pure Julia | YES | YES | NO (Fortran wrapper) |
| Memory control | `krylovdim` parameter controls subspace size | `maxdim` parameter | `ncv` parameter |
| Actively maintained | YES (Oct 2025 release) | YES (Feb 2025 release v0.4.0) | Maintenance-mode |
| AD/differentiability | YES (gradient support) | NO | NO |

**The decisive factor is `:LR` + callable interface + `bieigsolve`.** ArnoldiMethod.jl is also excellent (pure Julia, `mul!`-based, good memory control via `ArnoldiWorkspace`), but it lacks `bieigsolve` for left eigenvectors and `exponentiate` for time evolution. Since we need left eigenvectors for the biorthogonal overlap diagnostics (DIAG-05), KrylovKit is the clear winner.

### New Development Dependency: TimerOutputs.jl

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| **TimerOutputs.jl** | 0.5 | Hierarchical timing/memory profiling of Krylov iterations | `@timeit` sections for matvec cost, orthogonalization, convergence -- critical for tuning `krylovdim` and understanding where time goes at 10-12 qubits |

**Confidence: HIGH** -- Standard Julia profiling package. Already used by ITensors.jl and other large Julia numerical packages.

BenchmarkTools.jl is already in `[extras]` for microbenchmarks. TimerOutputs.jl complements it for coarse-grained profiling of multi-minute Krylov runs where `@benchmark` is impractical.

### Existing Dependencies: No Changes Required

| Existing Dep | Role in Krylov Milestone | Change Needed |
|---|---|---|
| **LinearAlgebra** (stdlib) | `mul!` for BLAS-3 matrix products in matvec, `BLAS.set_num_threads()` for thread control | NONE |
| **LinearMaps.jl** (3) | Already in `[deps]`. Could optionally wrap the matvec as `LinearMap` for Arpack compatibility, but KrylovKit does NOT need it (accepts bare callables) | NONE (keep for backward compat with existing dense Arpack path) |
| **Arpack.jl** (0.5.4) | Keep for dense Lindbladian at n<=6. The existing `run_lindbladian` with `eigs(L, sigma=shift)` remains the gold standard for small systems | NONE |
| **FINUFFT.jl** (3) | NUFFT prefactor computation for Time/Trotter domains -- used INSIDE the matvec | NONE |
| **SparseArrays** (stdlib) | Sparse Bohr-frequency indexing inside matvec | NONE |

---

## How KrylovKit.jl Integrates with QuantumFurnace

### The Matrix-Free Linear Map

KrylovKit's `eigsolve` accepts a plain Julia function `f(x) -> y` as the linear operator. The Lindbladian action `L(rho)` on a vectorized density matrix is:

```julia
function apply_lindbladian_matvec(v::Vector{ComplexF64},
                                   jumps, config, ham_or_trott, precomputed_data)
    dim = isqrt(length(v))
    rho = reshape(v, dim, dim)      # View, no allocation
    d_rho = zeros(ComplexF64, dim, dim)  # Preallocated workspace

    # Accumulate: d_rho = sum_k L_k(rho) = sum_k [A_k rho A_k' - 0.5{A_k'A_k, rho}]
    for jump in jumps
        for (w, A_w) in frequency_components(jump, config, precomputed_data)
            # A_w rho A_w'
            tmp = A_w * rho * A_w'       # Two BLAS-3 mul! calls
            d_rho .+= gamma_w .* tmp
            # -0.5 * {A_w'A_w, rho}
            LdL = A_w' * A_w
            d_rho .-= 0.5 .* gamma_w .* (LdL * rho + rho * LdL)
        end
    end

    # Add coherent part if applicable: -i[B, rho]
    if B !== nothing
        d_rho .-= 1im .* (B * rho - rho * B)
    end

    return vec(d_rho)
end
```

**Key point:** Each matvec costs O(n_jumps * n_freqs * dim^2 * dim) for the matrix multiplications. At 12 qubits (dim=4096), each `mul!` is a 4096x4096 complex matrix multiply = O(dim^3) = O(68 billion) flops. This is expensive but feasible with BLAS-3.

**The existing code already has the building blocks:**
- `_jump_contribution!` for `LiouvConfig` (Liouvillian construction) accumulates `A_w rho A_w' - 0.5{...}` into a vectorized matrix via `_vectorize_liouv_diss_and_add!`.
- `_jump_contribution!` for `ThermalizeConfig` (DM simulator) applies the CPTP channel to a density matrix directly.
- The Krylov matvec is a hybrid: it applies the Lindbladian to a density matrix (like the DM simulator) but returns the result as a vector (like the Liouvillian constructor).

The new code extracts the "accumulate L(rho)" logic from the existing `_jump_contribution!` functions, removing the Kraus/TP-correction overhead (which is only for the DM simulator's CPTP channel).

### The eigsolve Call

```julia
using KrylovKit

# Construct closure capturing all precomputed data
function make_lindbladian_map(jumps, config, ham_or_trott, precomputed_data)
    dim = size(ham_or_trott.data, 1)
    # Preallocate workspace (reused across matvec calls)
    ws = KrylovMatvecWorkspace(dim)  # new struct with scratch matrices

    function L_matvec(v::Vector{ComplexF64})
        return apply_lindbladian_action!(ws, v, jumps, config, ham_or_trott, precomputed_data)
    end
    return L_matvec
end

L = make_lindbladian_map(jumps, config, ham_or_trott, precomputed_data)
x0 = randn(ComplexF64, dim^2)  # Random initial vector

# Find 2-5 eigenvalues with largest real part (steady state + gap)
vals, vecs, info = eigsolve(L, x0, 5, :LR;
    krylovdim = 30,       # Krylov subspace size
    maxiter = 300,         # Max restarts
    tol = 1e-8,            # Convergence tolerance
    eager = true,          # Check convergence after each expansion
    verbosity = 1,         # Print convergence info
)

# vals[1] ~ 0 (steady state), vals[2] = spectral gap eigenvalue
spectral_gap = abs(real(vals[2]))
```

### Left Eigenvectors via bieigsolve

For the overlap diagnostics (DIAG-05), we need biorthogonal left eigenvectors:

```julia
# BiArnoldi: simultaneous left and right eigenvectors
vals, (vecs_right, vecs_left), (info_r, info_l) = bieigsolve(
    L, x0, x0, 5, :LR;
    krylovdim = 30,
    maxiter = 300,
    tol = 1e-8,
)
# vecs_left satisfy L' * w = conj(lambda) * w
# Biorthogonality: dot(vecs_left[i], vecs_right[j]) ~ delta_{ij}
```

**Important:** `bieigsolve` requires BOTH `L(v)` and `L'(w)` (adjoint action). The adjoint Lindbladian `L'(rho)` has a known form:
`L'(rho) = sum_k gamma_k * (A_k' rho A_k - 0.5{A_k'A_k, rho})` (swap A and A'). This is straightforward to implement alongside the forward map.

---

## Memory Budget Analysis

### Per-Vector Cost

At n qubits: dim = 2^n, vectorized = dim^2 = 4^n.

| Qubits | dim | dim^2 | Bytes per vector (ComplexF64) | Bytes |
|--------|-----|-------|------------------------------|-------|
| 6 | 64 | 4,096 | 32 KB | Trivial |
| 8 | 256 | 65,536 | 512 KB | Trivial |
| 10 | 1,024 | 1,048,576 | 8 MB | Manageable |
| 11 | 2,048 | 4,194,304 | 32 MB | Significant |
| 12 | 4,096 | 16,777,216 | 128 MB | Heavy |

### KrylovKit Memory Model

The Arnoldi method stores a Krylov basis of `krylovdim` vectors plus a few scratch vectors. Total Krylov memory:

```
memory_krylov = (krylovdim + ~5) * bytes_per_vector
```

For `bieigsolve`, double this (both left and right Krylov bases):

```
memory_bieigsolve = 2 * (krylovdim + ~5) * bytes_per_vector
```

| Qubits | krylovdim=30 eigsolve | krylovdim=30 bieigsolve | krylovdim=50 eigsolve |
|--------|----------------------|------------------------|----------------------|
| 8 | 18 MB | 36 MB | 28 MB |
| 10 | 288 MB | 576 MB | 448 MB |
| 11 | 1.1 GB | 2.2 GB | 1.8 GB |
| 12 | 4.5 GB | 9.0 GB | 7.0 GB |

### Matvec Workspace

The matvec itself needs scratch matrices for the dim x dim operations:

```
workspace = ~6 * dim * dim * 16 bytes  # 6 scratch matrices of size dim x dim
```

| Qubits | dim | Workspace |
|--------|-----|-----------|
| 8 | 256 | 6 MB |
| 10 | 1,024 | 96 MB |
| 12 | 4,096 | 1.5 GB |

### Total Memory Budget

| Qubits | Krylov (k=30) | Workspace | NUFFT Prefactors | Total | Fits in... |
|--------|---------------|-----------|------------------|-------|------------|
| 8 | 18 MB | 6 MB | ~50 MB | ~75 MB | Laptop (16 GB) |
| 10 | 288 MB | 96 MB | ~200 MB | ~600 MB | Laptop (16 GB) |
| 11 | 1.1 GB | 384 MB | ~400 MB | ~2 GB | Workstation (64 GB) |
| 12 | 4.5 GB | 1.5 GB | ~800 MB | ~7 GB | Cluster node (512 GB) |

**Verdict:** 12 qubits with `krylovdim=30` needs ~7 GB for `eigsolve`, ~12 GB for `bieigsolve`. Well within the 512 GB cluster budget. Even a 64 GB workstation handles 11 qubits comfortably.

**Note on KrylovKit memory overhead:** [Issue #9](https://github.com/Jutho/KrylovKit.jl/issues/9) reports KrylovKit allocates ~20x more than Arpack for equivalent problems due to its flexible vector type design. For our case, the overhead is in the orthogonalization, not the vector storage. The 20x factor applies to small problems; for our large vectors, the orthogonalization overhead (O(krylovdim^2 * dim^2)) is dwarfed by vector storage. Real overhead is likely 1.5-2x the theoretical minimum. Budget accordingly.

---

## BLAS Threading Strategy

### The Problem

Each Krylov matvec performs multiple dim x dim matrix multiplications via BLAS-3 (`mul!`). Julia's BLAS (OpenBLAS or MKL) uses its own thread pool. If Julia also uses multiple threads, BLAS threads and Julia threads compete for CPU cores, causing oversubscription and cache thrashing.

### The Solution

For Krylov eigensolving (inherently sequential outer loop), set BLAS threads high and Julia threads to 1:

```julia
using LinearAlgebra
BLAS.set_num_threads(N_physical_cores)  # Let BLAS parallelize the mul! calls
# Julia threading is NOT used in the Krylov loop itself
```

KrylovKit's Arnoldi loop is sequential: matvec -> orthogonalize -> check convergence -> repeat. The parallelism lives inside each matvec (BLAS-3 matrix multiply).

| System | BLAS threads | Julia threads | Rationale |
|--------|-------------|---------------|-----------|
| Laptop (4 cores) | 4 | 1 | BLAS owns all cores for mat-mul |
| Workstation (16 cores) | 16 | 1 | Same -- Krylov loop is sequential |
| Cluster node (64 cores) | 32-64 | 1 | OpenBLAS saturates at ~32 threads for large matrices; MKL scales better |

**At 12 qubits (4096x4096 matrices):** Each BLAS-3 `mul!` is O(68G) flops. With 64 BLAS threads on a cluster node, each `mul!` takes ~1-2 seconds. A matvec with ~100 frequency components and 3 jumps performs ~1000 `mul!` calls = ~1000-2000 seconds per matvec. With `krylovdim=30` and 100 restarts, this is ~1-2 days of compute. This motivates:

1. **Precomputation of jump-OFT products** to reduce the number of `mul!` calls per matvec.
2. **Exploiting Hermitian symmetry** of jumps (half the frequency grid).
3. **Starting with smaller systems** (n=8-10) to validate before scaling.

---

## Preconditioning and Convergence Acceleration

### Why NOT Shift-Invert

Shift-invert requires solving `(L - sigma*I) \ v` at each Krylov step, which requires either:
- Factorizing the dim^2 x dim^2 matrix (impossible for n > 6), OR
- Iteratively solving with GMRES (each GMRES iteration needs a matvec, creating a nested Krylov loop that is expensive and numerically delicate).

Since `:LR` targeting works directly (lambda_1=0 is extremal), shift-invert is unnecessary.

### Convergence Considerations

The gap between `lambda_1 = 0` and `lambda_2` determines convergence rate. If the spectral gap is small (slow mixing), the ratio `|lambda_2| / |lambda_3|` may be close to 1, requiring more Krylov iterations. Typical convergence:

- **Well-separated gap:** `krylovdim=30`, `maxiter=100` should suffice.
- **Clustered eigenvalues near 0:** Increase `krylovdim` to 50-80 and `maxiter` to 500. Memory cost grows linearly.
- **Very small gap (slow mixing):** The Krylov method may struggle if `|lambda_2|` is much smaller than `|lambda_end|` (the most negative eigenvalue). In this regime, consider working with the CPTP channel `E = I + delta*L` instead, where `E`'s leading eigenvalue is 1 (maps to steady state) and the second eigenvalue is `1 + delta*lambda_2` (maps to gap). For small delta, the gap eigenvalue of E is closer to the extremal eigenvalue 1, improving Krylov convergence.

### Alternative: Power Iteration on the CPTP Channel

For the CPTP channel `E(rho)`, the leading eigenvalue is 1 (largest magnitude). This makes `:LM` targeting natural. The spectral gap of L maps to `|1 - mu_2|/delta` where `mu_2` is the second eigenvalue of E. This approach:

- Uses the EXISTING `_jump_contribution!` for `ThermalizeConfig` directly (the DM simulator step IS the channel application).
- Converges faster when the gap is very small (because `mu_2 ~ 1` is near the extremal eigenvalue `mu_1 = 1`).
- BUT: introduces discretization error (delta-dependent). Richardson extrapolation can correct this.

**Recommendation:** Implement both L-based (`eigsolve` with `:LR`) and E-based (`eigsolve` with `:LM`) approaches. Use L-based as the primary method (exact, no delta dependence). Use E-based as a fallback/cross-validation (faster convergence for small gaps, but needs delta extrapolation).

---

## Alternatives Considered

| Category | Recommended | Alternative | Why Not |
|----------|-------------|-------------|---------|
| Krylov solver | **KrylovKit.jl** (0.10.2) | ArnoldiMethod.jl (0.4.0) | No `bieigsolve` (left eigvecs), no `exponentiate`. ArnoldiMethod has better memory control (`ArnoldiWorkspace` pre-allocation) but lacks critical features. |
| Krylov solver | KrylovKit.jl | Arpack.jl (0.5.4, existing) | Arpack needs explicit matrix for shift-invert. Cannot do matrix-free. Keep for n<=6 dense path. |
| Krylov solver | KrylovKit.jl | Krylov.jl (JSO ecosystem) | Krylov.jl focuses on linear systems (GMRES, CG, etc.), not eigenvalue problems. No `eigsolve`. |
| Linear map wrapper | Direct callable (KrylovKit native) | LinearMaps.jl (existing) | KrylovKit accepts bare functions. No need to wrap. LinearMaps.jl stays for Arpack backward compat. |
| Profiling | TimerOutputs.jl | BenchmarkTools.jl (existing) | BenchmarkTools for micro-benchmarks, TimerOutputs for profiling multi-minute Krylov runs with section-level breakdown. Complementary. |
| BLAS | OpenBLAS (Julia default) | MKL.jl | MKL may be faster for large matrices on Intel CPUs, but adds proprietary dependency. Benchmark first with OpenBLAS, switch to MKL only if profiling shows BLAS is bottleneck. |

---

## Dependencies to NOT Add

| Package | Why Evaluated | Why NOT Adding |
|---------|---------------|----------------|
| **ArnoldiMethod.jl** | Pure Julia alternative to KrylovKit | Missing `bieigsolve` and `exponentiate`. If we only needed forward eigsolve, ArnoldiMethod would be competitive. But left eigenvectors for DIAG-05 overlap analysis are essential. |
| **ArnoldiMethodTransformations.jl** | Shift-invert for ArnoldiMethod | Shift-invert is unnecessary (`:LR` works for Lindbladian). Would add complexity for no benefit. |
| **Krylov.jl** | JSO linear solver ecosystem | No eigenvalue solver. Wrong tool for this problem. |
| **ExponentialUtilities.jl** | Matrix exponential action `expmv` | KrylovKit.jl already provides `exponentiate` which does the same thing. No need for a second package. |
| **MKL.jl** | Potentially faster BLAS for Intel | Premature optimization. Benchmark with OpenBLAS first. MKL adds non-trivial binary dependency. |
| **CUDA.jl** / **GPUArrays.jl** | GPU acceleration of matvec | The matvec is many small (dim x dim) BLAS-3 calls, not one giant matmul. GPU data transfer overhead would dominate. Defer to a future GPU milestone if needed. |
| **Distributed.jl** (for Krylov) | Distributed matvec | The Krylov loop is sequential. Distributing the matvec across nodes adds MPI-like complexity. The parallelism should be within BLAS. |
| **IncompleteLU.jl** / **ILUZero.jl** | Preconditioners for shift-invert GMRES | Shift-invert is not needed. These are for sparse linear system preconditioning. |

---

## Project.toml Changes

### Production Dependencies [deps]: Add KrylovKit

```toml
[deps]
# ... existing deps unchanged ...
KrylovKit = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"

[compat]
# ... existing compat unchanged ...
KrylovKit = "0.10"
```

### Development/Profiling Dependencies [extras]: Add TimerOutputs

```toml
[extras]
# ... existing extras unchanged ...
TimerOutputs = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"

[compat]
# ... existing compat unchanged ...
TimerOutputs = "0.5"
```

### Test Targets: Add TimerOutputs

```toml
[targets]
test = ["Test", "BenchmarkTools", "Debugger", "StableRNGs", "HypothesisTests", "StatsBase", "Aqua", "TimerOutputs"]
```

---

## Installation

```bash
# In Julia REPL:
using Pkg
Pkg.add("KrylovKit")
Pkg.add("TimerOutputs")  # for benchmarking scripts
```

---

## Integration Points with Existing Code

### New Files (src/)

| New File | Dependencies Used | Purpose |
|----------|-------------------|---------|
| `src/krylov_matvec.jl` | LinearAlgebra (`mul!`), existing domain logic | Matrix-free `L(rho)` and `L'(rho)` actions, with preallocated workspace |
| `src/krylov_eigsolve.jl` | KrylovKit (`eigsolve`, `bieigsolve`) | Krylov-based spectral gap estimation: wrapper around KrylovKit with Lindbladian-specific defaults |
| `src/krylov_channel.jl` | KrylovKit (`eigsolve`) | CPTP channel eigensolve (`E(rho)` with `:LM`), delta extrapolation |

### Modified Files (src/)

| File | Change | Stack Impact |
|------|--------|--------------|
| `src/QuantumFurnace.jl` | `using KrylovKit`, new `include()` lines, new exports | One new `using` statement |
| `src/structs.jl` | New `KrylovMatvecWorkspace` struct, new `KrylovGapResult` result struct | No new deps |
| `src/diagnostics.jl` | Optional: adapt `extract_leading_eigendata` to accept Krylov results alongside dense eigen | No new deps |

### New Files (experiments/)

| New File | Dependencies Used | Purpose |
|----------|-------------------|---------|
| `experiments/krylov/benchmark_matvec.jl` | TimerOutputs, BenchmarkTools | Profile single matvec cost across qubit counts |
| `experiments/krylov/validate_against_dense.jl` | Existing dense eigen, new Krylov | Cross-validate Krylov gap vs dense gap at n<=6 |
| `experiments/krylov/scaling_study.jl` | TimerOutputs | Measure time and memory scaling n=4..12 |

---

## Verification Plan

### Step 1: Cross-Validate at n=4,6 (Dense Reference Exists)

```julia
# Dense reference (existing)
L_dense = construct_lindbladian(jumps, config, hamiltonian)
F = eigen(L_dense)
gap_dense = abs(real(sort(F.values, by=x->abs(real(x)))[2]))

# Krylov
L_map = make_lindbladian_map(jumps, config, ham_or_trott, precomputed_data)
vals, _, _ = eigsolve(L_map, x0, 5, :LR; krylovdim=30, tol=1e-10)
gap_krylov = abs(real(vals[2]))

@assert abs(gap_dense - gap_krylov) / gap_dense < 1e-6
```

### Step 2: Scale to n=8 (Krylov Only)

At n=8 (dim=256, dim^2=65536): dense Lindbladian is 65536 x 65536 = 32 GB. Too large for dense, but trivial for Krylov (each vector is 512 KB). This is the first "Krylov-only" regime.

### Step 3: Scale to n=10, 12

Systematic benchmarks tracking: time per matvec, number of matvecs to convergence, total wall time, peak memory.

---

## Sources

- [KrylovKit.jl official documentation](https://jutho.github.io/KrylovKit.jl/stable/) -- API reference, algorithm descriptions. HIGH confidence. Generated October 10, 2025.
- [KrylovKit.jl GitHub releases](https://github.com/Jutho/KrylovKit.jl/releases) -- v0.10.2, October 11, 2025. HIGH confidence.
- [KrylovKit.jl eigenvalue problems](https://jutho.github.io/KrylovKit.jl/stable/man/eig/) -- `eigsolve`, `bieigsolve` API with `:LR`/`:SR`/`:LM` selectors. Explicit note: "since no (shift-and)-invert is used, this will only be successful if you somehow know that eigenvalues close to zero are also close to the periphery of the spectrum." For Lindbladians, lambda=0 IS the periphery under `:LR`. HIGH confidence.
- [KrylovKit.jl available algorithms](https://jutho.github.io/KrylovKit.jl/stable/man/algorithms/) -- Arnoldi parameters: `krylovdim` (default 30), `maxiter` (default 100), `tol` (default 1e-12). BiArnoldi for `bieigsolve`. HIGH confidence.
- [KrylovKit.jl memory allocation issue #9](https://github.com/Jutho/KrylovKit.jl/issues/9) -- Reports ~20x memory overhead vs Arpack for small problems due to flexible vector type design. MEDIUM confidence (2018 issue, may have improved).
- [ArnoldiMethod.jl documentation](https://julialinearalgebra.github.io/ArnoldiMethod.jl/dev/) -- `partialschur` API, `ArnoldiWorkspace` for pre-allocation. v0.4.0 (Feb 2025). HIGH confidence.
- [ArnoldiMethod.jl releases](https://github.com/JuliaLinearAlgebra/ArnoldiMethod.jl/releases) -- v0.4.0 released February 22, 2025. HIGH confidence.
- [BifurcationKit eigensolver docs](https://github.com/bifurcationkit/BifurcationKitDocs.jl/blob/main/docs/src/eigensolver.md) -- Documents how to implement custom shift-invert with KrylovKit GMRES. MEDIUM confidence (third-party docs).
- [Krylov.jl performance tips](https://jso.dev/Krylov.jl/dev/tips/) -- BLAS threading recommendations: use physical cores, not logical. HIGH confidence.
- [Julia BLAS threading issue #27962](https://github.com/JuliaLang/julia/issues/27962) -- `BLAS.set_num_threads()` overhead and interaction with Julia threads. HIGH confidence.
- [TimerOutputs.jl GitHub](https://github.com/KristofferC/TimerOutputs.jl) -- `@timeit` macro for hierarchical profiling. HIGH confidence.
- [QuantumToolbox.jl](https://qutip.org/QuantumToolbox.jl/stable/users_guide/steadystate) -- Uses `eigsolve` for Lindbladian steady state. Validates the `:SR`/`:LR` approach. MEDIUM confidence (different codebase but same physics).
- [Arnoldi-Lindblad time evolution paper](https://arxiv.org/pdf/2109.01648) -- Discusses Arnoldi methods applied to Lindbladian operators, validates the Krylov approach for open quantum systems. HIGH confidence (peer-reviewed).

---

*Stack research for: QuantumFurnace.jl Krylov-based Lindbladian spectral gap estimation*
*Researched: 2026-02-20*
