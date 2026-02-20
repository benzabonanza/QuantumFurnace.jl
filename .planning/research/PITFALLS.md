# Domain Pitfalls: Krylov-Based Lindbladian Spectral Gap Estimation

**Domain:** Adding matrix-free Krylov eigensolving (KrylovKit.jl) to existing Lindbladian simulator (QuantumFurnace.jl) for spectral gap estimation at n=8..12 qubits
**Researched:** 2026-02-20
**Confidence:** HIGH (codebase analysis of furnace.jl/gap_estimation.jl/diagnostics.jl/qi_tools.jl + KrylovKit.jl official documentation + Arpack.jl existing usage patterns + numerical linear algebra literature)

**Relationship to prior research:** This document covers pitfalls specific to the Krylov eigensolving milestone. The v1.4 PITFALLS.md (2026-02-19) covered trajectory-based diagnostics pitfalls. This document focuses on the distinct challenges of integrating matrix-free Krylov methods into the existing dense-eigendecomposition and trajectory infrastructure.

---

## Critical Pitfalls

Mistakes that produce wrong spectral gap values, cause out-of-memory crashes, or require architectural rewrites.

---

### Pitfall 1: Targeting Interior Eigenvalues Without Shift-Invert (Wrong Eigenvalues Returned)

**What goes wrong:**
The Lindbladian L has eigenvalue 0 (steady state) and all other eigenvalues with Re(lambda) < 0. The spectral gap eigenvalue lambda_1 has the smallest nonzero |Re(lambda)|, making it an **interior** eigenvalue -- neither extremal in magnitude nor in real part. Calling `KrylovKit.eigsolve(L_matvec, v0, 2, :LR)` targets the eigenvalues with largest real part. This correctly finds lambda_0 = 0 (steady state), but the second eigenvalue returned is the one with the **second-largest** real part -- which IS the spectral gap eigenvalue lambda_1 **only if** Re(lambda_1) > Re(lambda) for all other eigenvalues. For the Lindbladian, this is exactly the definition of the spectral gap, so `:LR` is the correct `which` parameter. However, there is a subtle trap: if eigenvalues exist with Re(lambda) very close to Re(lambda_1) but different imaginary parts, Arnoldi may converge to those instead, or oscillate between candidates without converging.

The CPTP channel E = I + delta*L has a different spectrum. Its eigenvalue 1 (steady state) is the largest in magnitude, and the spectral gap eigenvalue has |mu_1| closest to 1 among non-steady-state eigenvalues. Using `:LM` correctly targets eigenvalue 1 first, then mu_1 second. But if delta is too large, some channel eigenvalues may have |mu| > 1 (violating complete positivity), and `:LM` returns those spurious eigenvalues instead.

**Why it happens:**
KrylovKit.jl does NOT implement shift-and-invert mode natively. The `EigSorter` and `ClosestTo` options sort the Krylov subspace eigenvalues but do NOT change the operator being applied. The documentation explicitly states: "no (shift-and)-invert is used." This means finding interior eigenvalues requires either (a) the target eigenvalues being naturally extremal under the chosen sorting, or (b) manually constructing a shifted-inverted operator `(L - sigma*I)^{-1}`, which requires solving a linear system at every Krylov iteration -- defeating the matrix-free advantage.

The existing code in `furnace.jl` (line 17) uses Arpack with `sigma=1e-9*(1+1im)` for shift-invert mode on a dense matrix. This approach requires the dense matrix to exist, which is impossible at n=12.

**How to avoid:**
1. **For the Lindbladian L:** Use `which=:LR` (largest real part). The spectral gap eigenvalue has the largest real part among all nonzero eigenvalues (it is the "least negative" real part). Request `howmany=4` or more eigenvalues to have a safety margin, then post-filter to find the one with smallest |Re(lambda)| among non-steady-state eigenvalues.
2. **For the CPTP channel E:** Use `which=:LM` (largest magnitude). The steady state eigenvalue 1 is the largest, and the spectral gap eigenvalue has the second-largest magnitude. Verify that no eigenvalue has |mu| > 1 + epsilon (CPTP violation).
3. **For L with eigenvalue clustering near the gap:** Increase `krylovdim` substantially (50-100 for n=8, 100-200 for n>=10). More Krylov vectors resolve clustered eigenvalues better.
4. **Validate against dense results at n=4,6.** The existing `extract_leading_eigendata` function (diagnostics.jl) returns exact eigenvalues via dense `eigen(L)`. Any Krylov gap estimate at n=4 or n=6 that disagrees with the dense gap by more than the requested tolerance is a bug.

**Warning signs:**
- `info.converged < howmany` -- fewer eigenvalues converged than requested.
- The Krylov gap at n=4 disagrees with `extract_leading_eigendata` gap by more than 1e-6.
- Multiple eigenvalues returned with nearly identical real parts but different imaginary parts (degenerate multiplet was not fully resolved).
- Re(lambda) > 0 for any returned eigenvalue of L (unphysical: Lindbladian eigenvalues must have Re <= 0).

**Phase to address:** Initial Krylov integration phase. This is the first thing to get right. Wrap in a validation test that compares Krylov vs dense at n=4.

**Memory/time context:** No additional memory cost -- this is about choosing the right algorithm parameters, not about memory.

---

### Pitfall 2: Krylov Subspace Memory Blowup at Large System Sizes

**What goes wrong:**
The Arnoldi iteration stores `krylovdim` vectors, each of length dim^2 = 4^n (ComplexF64 = 16 bytes per element). The memory cost is:

| n (qubits) | dim | dim^2 (vector length) | Bytes per vector | krylovdim=30 | krylovdim=100 |
|-------------|-----|-----------------------|------------------|--------------|---------------|
| 4 | 16 | 256 | 4 KB | 120 KB | 400 KB |
| 6 | 64 | 4,096 | 64 KB | 1.9 MB | 6.4 MB |
| 8 | 256 | 65,536 | 1 MB | 30 MB | 100 MB |
| 10 | 1,024 | 1,048,576 | 16 MB | 480 MB | 1.6 GB |
| 12 | 4,096 | 16,777,216 | 256 MB | 7.5 GB | 25 GB |

At n=12 with krylovdim=100, the Krylov subspace alone consumes 25 GB. Add the scratch vectors for the matvec (Kraus operators, density matrix buffers) and orthogonalization workspace, and the total easily exceeds 32 GB. On a 512 GB cluster node this is fine, but on a 4-thread laptop for benchmarking, n=10 with krylovdim=100 already consumes 1.6 GB just for the Krylov basis.

Additionally, KrylovKit.jl has been documented to allocate ~20x more memory than expected due to internal copies during the orthonormalization process (see GitHub issue #9). Each restart cycle may trigger significant temporary allocations.

**Why it happens:**
The Arnoldi method inherently requires storing the full Krylov subspace for orthogonalization. The Krylov-Schur variant used by KrylovKit.jl keeps a portion of the subspace after thick restarts, but the peak memory during each cycle is `krylovdim` full vectors. At n=12, each vector is 256 MB, so even modest krylovdim values create large memory footprints.

The second factor is that KrylovKit.jl's internal `basistransform!` function in `orthonormal.jl` creates temporary copies during the modified Gram-Schmidt process with iterative refinement (the default orthogonalizer `ModifiedGramSchmidtIR()`).

**How to avoid:**
1. **Compute memory budget before running.** Before launching eigsolve, compute `mem_estimate = krylovdim * 4^n * 16 * 1.5` (the 1.5x factor accounts for KrylovKit internals). If this exceeds available RAM, reduce krylovdim or abandon the run.
2. **Start with small krylovdim and increase.** Use krylovdim=30 as default (KrylovKit's default). If convergence fails, increase to 50, then 100. Log memory usage at each step.
3. **Use maxiter to compensate for small krylovdim.** More restarts with a smaller subspace can converge with less peak memory. Set maxiter=200-500 when krylovdim is kept small.
4. **For n=12, target krylovdim=30-50.** At krylovdim=50, the Krylov subspace is ~12.5 GB. With scratch space and Julia overhead, total memory should be under 20 GB. This is feasible on a 512 GB node.
5. **Profile with `Base.summarysize` or `@allocated`.** Wrap the eigsolve call and measure peak allocation to catch memory surprises early.
6. **Consider Float32 for exploratory runs at n=12.** Halves the Krylov subspace memory. But see Pitfall 7 on precision -- only use this for rough gap estimates, not for validation.

**Warning signs:**
- Julia GC activity spikes (visible via `GC.gc_live_bytes()` or `--track-allocation=user`).
- System swap usage increases (check `/proc/meminfo` or `free -h`).
- OOM kill on Linux (the process simply disappears).
- Eigsolve taking hours at n=12 when it took minutes at n=10 (may indicate thrashing rather than computation).

**Phase to address:** Memory budgeting phase. Must be established BEFORE attempting n>=10 runs. Build a `memory_estimate(n_qubits, krylovdim)` utility function.

---

### Pitfall 3: Vectorization Convention Mismatch (Column-Stacking vs Row-Stacking)

**What goes wrong:**
The existing dense Lindbladian construction in `qi_tools.jl` uses the column-stacking (Watrous) convention for vectorization:

```
vec(rho) = [rho[:,1]; rho[:,2]; ...]   (column-major, Julia default)
L * vec(rho) = vec(L_super(rho))
```

where the superoperator satisfies:
```
L_super(rho) = sum_k (A_k rho A_k^dagger - 0.5 {A_k^dagger A_k, rho})
```

The vectorized form uses `kron(A, conj(A))` for the sandwich term (line 77 of qi_tools.jl) and `kron(A^dagger A, I)` and `kron(I, (A^dagger A)^T)` for the anticommutator terms (lines 80-81). This is the standard column-stacking convention where `vec(AXB) = (B^T kron A) vec(X)`.

When implementing the matrix-free matvec, the operator must apply the Lindbladian to a vectorized density matrix and return the vectorized result. The crucial operation is:

```julia
rho = reshape(v, dim, dim)     # un-vectorize
result = L_action(rho)          # apply Lindbladian
return vec(result)              # re-vectorize
```

If `reshape` and `vec` use different conventions (e.g., if someone writes `reshape(v, dim, dim)'` to get row-major or uses `permutedims`), the result is silently wrong. Julia's `vec` and `reshape` are both column-major by default, so this is consistent **as long as nobody transposes**.

The specific trap: the existing Kraus-based step function (`_jump_contribution!` for thermalization in jump_workers.jl) operates on density matrices directly (not on vectorized form). When wrapping this for Krylov matvec, the temptation is to copy the Kraus application logic. But the Kraus operators are applied as `A_w * rho * A_w^dagger` (sandwich), while the vectorized Lindbladian applies `kron(A_w, conj(A_w)) * vec(rho)`. These are equivalent ONLY under column-major vectorization.

**Why it happens:**
The existing codebase has two representation pathways that have never needed to interoperate:
1. Dense Lindbladian (dim^2 x dim^2 matrix) built via `construct_lindbladian` using `_kron!` and `_vectorize_liouv_diss_and_add!`. This uses column-major vec.
2. Trajectory stepping (`step_along_trajectory!`) which works with state vectors psi, never with vectorized density matrices.

The new Krylov path introduces a third pathway: matrix-free matvec that must be consistent with pathway (1) for validation but uses the computational approach of pathway (2). Mixing conventions between these pathways produces wrong eigenvalues.

**How to avoid:**
1. **Use Julia's native `vec()` and `reshape(v, dim, dim)` everywhere.** Never use `permutedims`, `transpose`, or `PermutedDimsArray` in the vectorization/unvectorization step.
2. **Write a round-trip test:** Build the dense Lindbladian L_dense for n=4. For 10 random density matrices rho, verify that `L_dense * vec(rho) == vec(L_matvec(rho))` to machine precision. This is the single most important validation test.
3. **Put the vectorization convention in a docstring and NEVER deviate.** The convention is: `vec(rho)` stacks columns, `reshape(vec_rho, dim, dim)` recovers the density matrix.
4. **Be explicit about the Kraus-to-Lindbladian relationship.** The Lindbladian action is `L(rho) = sum_k (A_k rho A_k^dagger) - rho` (after normalization). The CPTP channel is `E(rho) = sum_k K_k rho K_k^dagger`. These are related by `L = (E - I) / delta`. Implement the matvec as the Lindbladian action, NOT as the channel action, to match the dense L matrix.

**Warning signs:**
- The round-trip test fails: `norm(L_dense * vec(rho) - vec(L_matvec(rho))) > 1e-10`.
- The Krylov steady-state eigenvector, when reshaped to a density matrix, is not Hermitian.
- The spectral gap from Krylov disagrees with dense by more than the tolerance but the steady state eigenvalue is correct (indicating the matvec is partially correct but has a sign or transpose error).

**Phase to address:** Core matvec implementation phase. The round-trip test must pass before ANY Krylov eigensolving is attempted.

---

### Pitfall 4: CPTP Channel vs Lindbladian Spectrum Confusion

**What goes wrong:**
The project has two operator representations:
1. **Lindbladian L**: eigenvalues are {0, lambda_1, lambda_2, ...} with Re(lambda_k) < 0. The spectral gap is |Re(lambda_1)|.
2. **CPTP channel E = I + delta*L**: eigenvalues are {1, mu_1, mu_2, ...} with |mu_k| < 1 (for small enough delta). The spectral gap is 1 - |mu_1|.

These are related by `mu_k = 1 + delta * lambda_k`, so `|Re(lambda_1)| = (1 - Re(mu_1)) / delta` approximately. But this relationship holds exactly only in the limit delta -> 0. For finite delta:
- The channel may not be CPTP if delta is too large (eigenvalues escape the unit disk).
- The channel eigenvalues have a delta-dependent shift: `mu_k = 1 + delta * lambda_k` is exact for the relationship between the generators, but the physical channel eigenvalues also acquire O(delta^2) corrections from the Trotterized evolution.

The dangerous mistake: computing the Krylov eigenvalues of the CPTP channel E (which is what the existing trajectory framework naturally provides as a matvec) and then forgetting to convert back to Lindbladian eigenvalues, or converting with the wrong delta, or using the Trotterized channel eigenvalues as if they were exact Lindbladian eigenvalues.

**Why it happens:**
The existing trajectory code (`step_along_trajectory!` in trajectories.jl) implements a single step of the CPTP channel. It is natural to wrap this as a matvec for the channel E, not for the Lindbladian L. But the diagnostics infrastructure (`extract_leading_eigendata`, `eigenbasis_overlap_analysis`) works with the Lindbladian L. Mixing these two representations leads to factor-of-delta errors or comparison bugs.

**How to avoid:**
1. **Decide up front: Krylov on L or Krylov on E.** The Lindbladian matvec computes `L(rho)` directly. The channel matvec computes `E(rho)`. They have different spectra.
2. **Recommendation: Implement Krylov on L (Lindbladian) for consistency with dense diagnostics.** The dense code uses `eigen(L)`. The Krylov code should target the same operator. Build the Lindbladian matvec from scratch using the Kraus operators: `L(rho) = (1/delta) * (E(rho) - rho)`.
3. **If using the channel E:** Target `:LM` (largest magnitude) to find eigenvalue 1 and the spectral gap eigenvalue. Convert: `lambda_k = (mu_k - 1) / delta`. Document this conversion prominently.
4. **Cross-validate:** At n=4, compute both the dense Lindbladian eigenvalues and the dense channel eigenvalues. Verify the relationship `mu_k = 1 + delta * lambda_k` holds to expected precision.
5. **Watch the delta dependence.** The spectral gap from the channel converges to the Lindbladian gap only as delta -> 0. If using the channel, run at multiple delta values and extrapolate (Richardson extrapolation from v1.4 diagnostics).

**Warning signs:**
- The gap from Krylov-on-channel is delta-dependent (changes significantly when delta is halved).
- The gap from Krylov-on-channel disagrees with dense L gap by more than O(delta) relative error.
- Eigenvalues of the channel have |mu| > 1 (CPTP violation, delta too large).

**Phase to address:** Matvec design phase. The L-vs-E decision must be made and documented before implementation. The recommended approach is to implement `L(rho) = (E(rho) - rho) / delta` to reuse existing Kraus infrastructure.

---

### Pitfall 5: Basis Mismatch Between Krylov Matvec and Dense Diagnostics

**What goes wrong:**
The existing codebase operates in different bases depending on the domain:
- **BohrDomain**: Hamiltonian eigenbasis. Jump operators stored as `jump.in_eigenbasis = eigvecs' * jump.data * eigvecs`.
- **TrotterDomain**: Trotter eigenbasis. Jump operators transformed via `trotter.eigvecs' * jump.data * trotter.eigvecs`.
- **Diagnostics** (`run_exact_diagnostics`): Takes the Lindbladian L and works in whichever basis L was constructed in.

The Krylov matvec must operate in the **same basis** as the dense Lindbladian was constructed in. If the matvec applies Kraus operators in the computational basis while the dense L was built in the eigenbasis, the eigenvalues will be correct (eigenvalues are basis-independent) but the eigenvectors will be in different bases, making cross-validation of eigenvectors impossible.

More dangerously: if the matvec uses jump operators from `jump.data` (computational basis) but the Hamiltonian evolution uses `hamiltonian.eigvecs` (eigenbasis), the intermediate steps may be inconsistent. The existing `construct_lindbladian` in `furnace.jl` carefully selects the basis via `ham_or_trott` (line 48-53) and uses `jump.in_eigenbasis`. The matrix-free path must replicate this basis selection exactly.

**Why it happens:**
The basis selection logic is scattered across multiple files: `furnace.jl` chooses `ham_or_trott`, `furnace_utensils.jl` precomputes domain-specific data, `jump_workers.jl` applies jumps in the chosen basis. When building a new matrix-free path, it is easy to accidentally use `jump.data` instead of `jump.in_eigenbasis`, or to forget that the TrotterDomain requires `trotter.eigvecs` instead of `hamiltonian.eigvecs`.

**How to avoid:**
1. **Reuse the existing precomputation pipeline.** Call `_precompute_data(config, ham_or_trott)` exactly as `construct_lindbladian` does. This ensures the same domain-specific setup.
2. **The matrix-free function should accept the same (jumps, config, hamiltonian; trotter) signature** as `construct_lindbladian`. Internally, it should replicate the basis selection: `ham_or_trott = config.domain isa TrotterDomain ? trotter : hamiltonian`.
3. **Validate at n=4 in BOTH BohrDomain and TrotterDomain.** Build the dense L for BohrDomain, compare with Krylov matvec for BohrDomain. Repeat for TrotterDomain. Both must match.
4. **Use `jump.in_eigenbasis` exclusively in the matvec, never `jump.data`.** The `.data` field is in computational basis and is only used for constructing `.in_eigenbasis`.

**Warning signs:**
- Round-trip test passes for BohrDomain but fails for TrotterDomain (or vice versa).
- Eigenvectors from Krylov, when reshaped and compared to dense eigenvectors, are rotated by a unitary transformation (basis mismatch).
- The fixed point from Krylov does not match the Gibbs state in the expected basis.

**Phase to address:** Matvec implementation phase. Include TrotterDomain tests from the beginning, not as an afterthought.

---

### Pitfall 6: Non-Convergence of Arnoldi for Clustered or Degenerate Eigenvalues

**What goes wrong:**
The Lindbladian of a system with symmetries (e.g., SZ conservation in the Heisenberg chain) has eigenvalue multiplets -- groups of eigenvalues with identical or nearly identical real parts but different imaginary parts. The existing diagnostics code (`detect_multiplets` in diagnostics.jl) already handles this for dense eigendecomposition. But Krylov methods struggle with degenerate or clustered eigenvalues:

1. **Degenerate eigenvalues:** If lambda_1 and lambda_2 have the same value, the Arnoldi iteration may return a random linear combination of their eigenvectors, or fail to converge to both.
2. **Near-degenerate eigenvalues:** If |lambda_1 - lambda_2| is very small relative to the Krylov subspace dimension, the Krylov-Schur restart may not separate them. The residual norm for their Schur vectors may stagnate above the tolerance.
3. **Complex conjugate pairs:** For a real Lindbladian (which this is NOT -- the Lindbladian is complex), eigenvalues come in conjugate pairs. For the complex Lindbladian here, this is not guaranteed, but near-conjugate pairs may still appear and cause convergence difficulties.

The n=4 Heisenberg chain has well-separated eigenvalues (the exact diagnostics show distinct multiplets). But at n=8+, the eigenvalue density near the spectral gap increases, and multiplets become denser. The probability of near-degeneracy at the gap boundary increases with system size.

**Why it happens:**
The Arnoldi iteration builds a Krylov subspace by repeatedly applying the operator. When two eigenvalues are close, the Krylov subspace cannot distinguish their eigenvectors until enough iterations have been performed. The number of iterations needed grows as ~O(1/gap_between_eigenvalues). With thick restarts, the Krylov-Schur algorithm can partially mitigate this, but it still struggles when the gap between eigenvalues is smaller than the achievable residual norm.

KrylovKit.jl has a documented issue (#23, #38) where degenerate eigenvalues cause convergence failure. The v0.10 release added `BlockLanczos` for Hermitian problems with degeneracies, and `BiArnoldi` for non-Hermitian problems, but the standard Arnoldi `eigsolve` still has this limitation.

**How to avoid:**
1. **Request more eigenvalues than needed.** Instead of `howmany=2`, use `howmany=6` or `howmany=10`. This gives the Krylov subspace room to resolve nearby eigenvalues. Post-filter to find the spectral gap.
2. **Increase krylovdim.** The rule of thumb is `krylovdim >= 3 * howmany` (KrylovKit documentation). For near-degenerate problems, use `krylovdim >= 5 * howmany`.
3. **Always check `info.converged`.** KrylovKit does NOT warn if fewer eigenvalues converge than requested (prior to v0.10 defaults). Always assert `info.converged >= howmany` and handle the failure case.
4. **Use multiple random starting vectors.** If convergence fails from one starting vector, try others. Different starting vectors project differently onto the eigenspaces and may help.
5. **Consider BiArnoldi** (KrylovKit v0.10+) if both left and right eigenvectors are needed, which the existing diagnostics infrastructure requires for overlap analysis.
6. **For the spectral gap specifically:** The gap eigenvalue must have converged. Even if some higher eigenvalues fail to converge, the gap estimate may still be reliable if the first few eigenvalues converged.

**Warning signs:**
- `info.converged == 0` or `info.converged == 1` when 2+ were requested.
- Running eigsolve twice with different random starting vectors gives different eigenvalues.
- The returned eigenvalues change significantly when krylovdim is increased by 50%.
- Residual norms plateau above the tolerance and do not decrease with more iterations.

**Phase to address:** Krylov convergence tuning phase. After basic matvec is validated, systematically test convergence at n=6,8,10 and build adaptive parameter selection.

---

### Pitfall 7: Numerical Precision Loss in Matrix-Free Matvec

**What goes wrong:**
The matrix-free Lindbladian action computes `L(rho) = (1/delta) * (sum_k K_k rho K_k^dagger - rho)` or equivalently the Lindblad dissipator directly. Each evaluation involves:
1. Applying Kraus/jump operators to a density matrix (matrix-matrix multiplications).
2. Subtracting the original density matrix (cancellation).
3. Summing over all jump operators and frequency components.

The subtraction in step 2 is catastrophic cancellation when `E(rho)` is close to `rho` (which it always is for small delta). If `delta = 0.01` and `E(rho) - rho` has relative magnitude delta, then dividing by delta recovers order-1 quantities, but the absolute precision is reduced by a factor of delta. In Float64, this means losing ~2 digits of precision for delta=0.01.

More seriously, if the matvec accumulates contributions from many jump operators and frequency components (as the existing `construct_lindbladian` does), rounding errors accumulate. For n=12 with 3*n=36 jump operators and ~2^11 frequency components each, the total number of floating-point operations per matvec is O(dim^2 * n_jumps * n_freqs) = O(16M * 36 * 2048) ~ 10^12 FLOPs. Rounding errors at the 10^{-16} level accumulate as sqrt(10^12) * 10^{-16} = 10^{-10}. This limits the achievable Krylov tolerance to ~10^{-10} at best for n=12.

Additionally, each Krylov iteration applies the matvec once and orthogonalizes against all previous Krylov vectors. Orthogonalization errors compound over iterations. With krylovdim=50 iterations, the accumulated orthogonality loss can be O(krylovdim * eps_matvec).

**Why it happens:**
The matrix-free approach trades memory for computation. Instead of storing the dim^2 x dim^2 Lindbladian matrix (impossible at n>=8), we recompute the action each time. But each recomputation is subject to rounding errors, and these errors are not identical across iterations (unlike a stored matrix, where the matrix-vector product has consistent rounding behavior). This means the linear operator seen by Arnoldi is not exactly the same at each iteration -- it has a tiny fluctuating perturbation from rounding errors.

The catastrophic cancellation in `(E(rho) - rho) / delta` is particularly insidious because it looks like the code is working (the magnitudes are correct) but the noise floor is elevated.

**How to avoid:**
1. **Implement the Lindbladian action directly, NOT as `(E(rho) - rho) / delta`.** Compute `L(rho) = sum_k (A_k rho A_k^dagger - 0.5 * {A_k^dagger A_k, rho})` directly, without the E-I-divide-by-delta pathway. This avoids the catastrophic cancellation. The existing `_vectorize_liouv_diss_and_add!` function computes exactly this in vectorized form; the matrix-free version should compute the same thing in operator form.
2. **Validate the matvec precision.** At n=4, compute `norm(L_dense * v - L_matvec(v)) / norm(L_dense * v)` for multiple random vectors. This relative error should be at most ~10^{-12} at n=4. Track how it degrades with n.
3. **Set the Krylov tolerance appropriately.** Do not set `tol=1e-12` if the matvec has 10^{-10} precision. A reasonable tolerance is `tol = max(1e-8, estimated_matvec_precision * 100)`.
4. **Use the `ModifiedGramSchmidtIR()` orthogonalizer** (KrylovKit default). The iterative refinement step re-orthogonalizes to combat loss of orthogonality. Do NOT switch to plain `ModifiedGramSchmidt()` to save time.
5. **Consider using the CPTP channel E directly with `:LM`** if the Lindbladian direct computation proves too noisy. The channel avoids the subtraction-and-division step.

**Warning signs:**
- The residual norm reported by KrylovKit plateaus at 10^{-8} or worse despite requesting 10^{-12}.
- At n=4, the Krylov eigenvalues match dense eigenvalues to only 6-8 digits instead of 12+.
- The Krylov eigenvectors are not orthogonal to each other (check `abs(dot(v1, v2)) > 1e-6`).
- The computed steady-state eigenvalue is not close to 0 (for L) or 1 (for E).

**Phase to address:** Matvec validation phase. Quantify the precision at each system size before trusting the eigenvalue results.

---

## Moderate Pitfalls

---

### Pitfall 8: BLAS Thread Contention During Krylov Iteration

**What goes wrong:**
The Krylov iteration involves dense matrix-vector products (BLAS Level 2) and matrix-matrix products (BLAS Level 3) for orthogonalization. If Julia is started with multiple threads (`-t 4`) and BLAS also uses multiple threads (default on most systems), the total thread count becomes `Julia_threads * BLAS_threads`. On a 4-core laptop, this means 16 threads competing for 4 cores, causing severe performance degradation from context switching and cache thrashing.

The existing trajectory code (`test_threading.jl`) already handles this correctly: it sets `BLAS.set_num_threads(1)` during multi-threaded trajectory sampling and restores afterward. But the Krylov iteration is a different use pattern -- it is inherently serial (each iteration depends on the previous), so Julia threading is not used for the eigensolve itself. However, the matvec (applying the Lindbladian) may benefit from BLAS threads for the matrix-matrix multiplications `A * rho * A'`.

**Why it happens:**
At n=12, `rho` is a 4096x4096 matrix. The Kraus application `A_w * rho * A_w'` involves two dim x dim matrix multiplications, each O(dim^3) = O(64 billion) FLOPs. BLAS multithreading is genuinely helpful here -- a single matvec at n=12 takes ~seconds with 1 BLAS thread and ~hundreds of ms with 4 BLAS threads. But if the Krylov iteration also launches Julia threads for something (e.g., parallel Kraus operator application), the two threading levels conflict.

**How to avoid:**
1. **Set BLAS threads = physical cores during Krylov iteration.** The matvec is compute-bound on BLAS operations. Let BLAS use all cores. Do not use Julia `@threads` for the matvec.
2. **Do NOT parallelize the inner Krylov loop.** KrylovKit's Arnoldi iteration is inherently serial. The only parallelism opportunity is inside the matvec.
3. **Profile the matvec at n=8 and n=10.** Measure wall-clock time per matvec with 1, 2, 4 BLAS threads. Choose the optimal setting.
4. **Use `BLAS.set_num_threads()` with a restore pattern** (as in the existing trajectory code). Wrap the eigsolve call in a try-finally that restores the original BLAS thread count.
5. **For the cluster (512 GB, many cores):** Use BLAS threads = 8-16 for the matvec. Do not launch multiple Krylov eigensolves in parallel (memory constraint dominates).

**Warning signs:**
- Krylov iteration is slower at n=10 than expected based on n=8 scaling (should scale as ~16x for dim doubling, not 100x).
- CPU utilization is >100% per core (visible in `top` or `htop`) indicating thread oversubscription.
- System load average >> number of physical cores.

**Phase to address:** Performance tuning phase, after correctness is established.

---

### Pitfall 9: Wrong Starting Vector for Krylov Iteration

**What goes wrong:**
KrylovKit's `eigsolve(f, v0, howmany, which)` requires an initial vector `v0` to seed the Krylov subspace. If `v0` has zero or negligible overlap with the target eigenspace, convergence is extremely slow or fails entirely. For the spectral gap, we need overlap with both the steady-state eigenvector (lambda=0) and the gap eigenvector (lambda_1).

The steady-state eigenvector is `vec(rho_beta)` (the vectorized Gibbs state). The gap eigenvector is unknown a priori. If `v0` is orthogonal to the gap eigenvector, the Krylov subspace will not contain any component of that eigenmode, and Arnoldi will converge to higher eigenvalues instead.

A common naive choice is `v0 = randn(ComplexF64, dim^2)`, which has overlap with all eigenmodes in expectation. This works but may converge slowly if the random vector has unusually small overlap with the gap mode.

**Why it happens:**
For n=12, dim^2 = 16M. A random vector has overlap ~O(1/sqrt(dim^2)) = O(1/4096) with any particular eigenvector. This is usually sufficient for the largest-real-part eigenvalues to emerge, but convergence may take many iterations. For the gap eigenvalue, which is the SECOND-largest real part, convergence requires the first eigenvalue to be resolved first (so the Krylov subspace can "deflate" it), and then the second to emerge from the remaining subspace.

**How to avoid:**
1. **Use `v0 = vec(rho_0 - rho_beta)` where `rho_0` is a physically motivated initial state** (e.g., maximally mixed state or the all-up state). This removes the steady-state component and enhances the gap-mode component.
2. **At minimum, use `v0 = randn(ComplexF64, dim^2)` normalized.** This has nonzero overlap with all eigenmodes with probability 1.
3. **For the channel E:** Use `v0 = vec(rho_0)` where `rho_0` is NOT the Gibbs state. The Gibbs state is the eigenvalue-1 eigenvector, and using it as v0 would miss the gap eigenvalue entirely.
4. **If convergence is slow, try multiple starting vectors.** Run eigsolve 3 times with different random seeds and take the result with the smallest residual.

**Warning signs:**
- `info.converged == 1` but `info.converged < howmany` (only the steady state converged).
- All returned eigenvalues are near 0 (for L) or near 1 (for E) -- the algorithm found the steady state but nothing else.
- The number of matvec applications (`info.numops`) is near `maxiter * krylovdim` (exhausted all iterations without convergence).

**Phase to address:** Initial Krylov implementation phase. Use the physically motivated starting vector from the beginning.

---

### Pitfall 10: Comparing Krylov vs Dense Results Incorrectly

**What goes wrong:**
Cross-validation requires comparing Krylov eigenvalues/eigenvectors with dense results at n=4,6. Several comparison mistakes are common:

1. **Phase ambiguity in eigenvectors.** Eigenvectors are defined only up to a complex phase. `v_krylov = exp(i*theta) * v_dense` is perfectly valid. Comparing `norm(v_krylov - v_dense)` will show a large difference even though the eigenvectors are the same.
2. **Ordering ambiguity.** The Krylov eigensolver returns eigenvalues in the order determined by `which`, while dense `eigen()` returns them in an implementation-specific order. The existing diagnostics code sorts by `abs.(real.(eigenvalues))` (line 173 of diagnostics.jl), but the Krylov results need the same sorting before comparison.
3. **Degenerate eigenvalue subspaces.** If two eigenvalues are equal (within tolerance), their eigenvectors span a 2D subspace. Any basis of this subspace is valid. Comparing individual eigenvectors is meaningless; instead, compare the subspace projection.
4. **Different left/right eigenvector normalization.** KrylovKit returns right eigenvectors normalized to unit norm. The dense code computes left eigenvectors via `inv(V_right)` (diagnostics.jl line 183). The biorthonormality `V_left' * V_right = I` may have a different normalization convention than KrylovKit's right eigenvectors.

**Why it happens:**
Dense eigendecomposition (LAPACK `eigen`) and Krylov methods use fundamentally different algorithms. They produce equivalent but differently parametrized results. Naive element-wise comparison fails.

**How to avoid:**
1. **Compare eigenvalues only for validation.** Sort both sets by real part (largest first), then compare `abs(lambda_krylov[k] - lambda_dense[k]) < tol` for the first few eigenvalues.
2. **For eigenvector comparison, use the residual norm.** For right eigenvector v: `residual = norm(L * v - lambda * v) / norm(v)`. This is phase-independent and should be near machine epsilon for dense, and near the Krylov tolerance for Krylov.
3. **For subspace comparison with degeneracies:** Compute the projection matrix `P = V * V'` for the degenerate subspace from both methods. Compare `norm(P_krylov - P_dense)`.
4. **Use the existing `extract_leading_eigendata` as ground truth.** Its eigenvalues are sorted consistently. Compare Krylov eigenvalues against `eigen_result.eigenvalues[1:howmany]`.
5. **Build a dedicated `validate_krylov_vs_dense(n_qubits, config, ...)` function** that performs all these comparisons and returns a pass/fail with diagnostic messages.

**Warning signs:**
- Eigenvalue comparison shows large differences for some eigenvalues but not others (ordering mismatch).
- Eigenvector comparison shows `norm(v_krylov - v_dense) ~ 2` (phase flip, not a real difference).
- Validation passes at n=4 but fails at n=6 (indicates a real bug, not a comparison artifact).

**Phase to address:** Validation framework phase, implemented alongside the core Krylov solver.

---

### Pitfall 11: Convergence Criteria -- When Is the Gap Estimate "Good Enough"?

**What goes wrong:**
The Krylov tolerance `tol` controls the residual norm of the Schur decomposition, NOT the absolute error in the eigenvalue. For non-Hermitian operators, the relationship between residual norm and eigenvalue error depends on the **condition number of the eigenvalue**, which can be much larger than 1 for non-normal operators.

For the Lindbladian (which is non-normal due to the anti-Hermitian part), the eigenvalue condition number `kappa(lambda_k) = 1 / |<l_k, r_k>|` where `l_k` and `r_k` are the left and right eigenvectors. The eigenvalue error satisfies `|lambda_exact - lambda_krylov| <= kappa(lambda_k) * residual_norm`. If `kappa(lambda_1) = 100` (moderate non-normality) and `residual_norm = 1e-8`, the eigenvalue error can be as large as `1e-6`.

For the spectral gap, which is the DIFFERENCE between two eigenvalues (lambda_0 = 0 and lambda_1), the gap error is bounded by the sum of the individual eigenvalue errors.

**Why it happens:**
The KMS similarity transform analysis (DIAG-03/04 in diagnostics.jl) already quantifies non-normality via the anti-Hermitian defect ratio. At low temperature (beta=10, n=6), the defect ratio can be O(1) or larger, meaning eigenvalue condition numbers are elevated. This problem worsens at larger system sizes and lower temperatures.

**How to avoid:**
1. **Request tighter tolerance than you need.** If you need the gap to 4 significant figures, set `tol=1e-10` (with the expectation that eigenvalue error is ~kappa * tol ~ 100 * 1e-10 = 1e-8, giving ~8 digits).
2. **Compute the eigenvalue condition number from the Krylov result.** KrylovKit returns right eigenvectors. If using BiArnoldi, you also get left eigenvectors. Compute `kappa = 1 / abs(dot(l_k, r_k))` and use it to bound the eigenvalue error.
3. **Use the residual norm as a diagnostic, not as the accuracy guarantee.** Report both the residual norm and the estimated eigenvalue accuracy (`kappa * residual`).
4. **Cross-validate with trajectory-based gap estimation** (the existing `estimate_spectral_gap`). The trajectory estimate has different error characteristics (statistical, not algebraic) and provides an independent check.
5. **For the gap specifically:** Compute at least 4-6 eigenvalues. If eigenvalues 2 and 3 have similar real parts, the gap may be a degenerate multiplet. Report the gap with uncertainty bounds.

**Warning signs:**
- The Krylov gap changes by more than 10% when tol is tightened from 1e-8 to 1e-10 (indicates the gap was not converged at the looser tolerance).
- The eigenvalue condition number `kappa` > 1000 (severe non-normality, gap accuracy may be poor).
- The Krylov gap disagrees with the trajectory-based gap by more than the trajectory confidence interval.

**Phase to address:** Gap estimation API phase, when wrapping the raw eigsolve into a user-facing function with error bounds.

---

## Minor Pitfalls

---

### Pitfall 12: Forgetting to Normalize the Krylov Eigenvector Before Reshaping to Density Matrix

**What goes wrong:**
The eigenvector from KrylovKit has unit norm in the vector sense (`norm(v) = 1`). When reshaped to a density matrix via `reshape(v, dim, dim)`, the result has `tr(rho) != 1` in general (because trace is not the same as vector norm). The existing code in `furnace.jl` (lines 26-27) correctly handles this:
```julia
steady_state_dm = reshape(steady_state_vec, size(hamiltonian.data))
hermitianize!(steady_state_dm)
steady_state_dm ./= tr(steady_state_dm)
```
Forgetting the `./= tr(...)` step means the density matrix is not normalized, which corrupts any downstream computation (trace distances, fidelity, expectation values).

**How to avoid:** Always normalize after reshaping. Copy the existing pattern from `furnace.jl` and `diagnostics.jl`.

**Phase to address:** Eigenvector post-processing, part of the initial implementation.

---

### Pitfall 13: KrylovKit Version Incompatibility

**What goes wrong:**
KrylovKit v0.10 introduced breaking API changes: `BiArnoldi`, `BlockLanczos`, and changes to the `ConvergenceInfo` struct. Code written against v0.8 may not work with v0.10 and vice versa. The `eigsolve` return signature also changed subtly (the info struct gained new fields).

**How to avoid:** Pin the KrylovKit version in `Project.toml`. Test against the pinned version. Add it to `[compat]` with a specific version range.

**Phase to address:** Dependency management, first thing when adding KrylovKit to the project.

---

### Pitfall 14: Not Exploiting Hermiticity of the Density Matrix in the Krylov Vector

**What goes wrong:**
The Krylov vectors are ComplexF64 vectors of length dim^2. But the density matrix they represent is Hermitian: `rho[i,j] = conj(rho[j,i])`. This means only dim*(dim+1)/2 real parameters are independent, not dim^2 complex parameters. Ignoring this means the Krylov subspace has twice the effective dimension it needs, wasting memory and iterations.

However, the Lindbladian is non-Hermitian, so the Krylov iteration (Arnoldi) cannot exploit the density matrix Hermiticity. The operator `L(rho)` does NOT preserve Hermiticity in general (it does physically, but the action on the vectorized form does not reduce to a real-symmetric problem). So this is NOT something that can be fixed by switching to Lanczos.

**How to avoid:** Accept that Arnoldi on the full dim^2-dimensional space is necessary. Do not attempt to project onto the "Hermitian subspace" during Krylov iteration -- this would break the Arnoldi recurrence.

**Phase to address:** No action needed. This is a "do not attempt this optimization" warning.

---

## Technical Debt Patterns

Shortcuts that seem reasonable but create long-term problems.

| Shortcut | Immediate Benefit | Long-term Cost | When Acceptable |
|----------|-------------------|----------------|-----------------|
| Using `(E(rho) - rho) / delta` for the Lindbladian matvec instead of direct Lindblad formula | Reuses existing Kraus infrastructure with minimal code | Catastrophic cancellation loses ~2 digits per factor of delta; limits achievable precision | Never for the Lindbladian. Acceptable only if eigensolving the channel E directly. |
| Hardcoding krylovdim=30 for all system sizes | Simple API, no parameter tuning needed | Convergence failure at n>=10 where eigenvalue density is higher; silent wrong results if info.converged not checked | Only for n<=8 where eigenvalues are well-separated |
| Skipping TrotterDomain validation | Faster development (only test BohrDomain) | Basis mismatch bugs surface late when TrotterDomain is needed for physical predictions | Never -- TrotterDomain is the physically relevant domain |
| Storing dense Kraus operators for the matvec | Simple code, each operator is a dim x dim matrix | 36 operators * dim^2 * 16 bytes = 36 * 256 MB = 9 GB at n=12 | Acceptable up to n=10 (~576 MB). At n=12, need lazy/on-the-fly Kraus computation. |

---

## Integration Gotchas

Common mistakes when connecting Krylov results to the existing QuantumFurnace infrastructure.

| Integration Point | Common Mistake | Correct Approach |
|-------------------|----------------|------------------|
| Krylov vs `extract_leading_eigendata` | Comparing eigenvalues without matching sort order | Sort both by `abs(real(lambda))`, then compare pairwise |
| Krylov eigenvectors vs `eigenbasis_overlap_analysis` | Passing Krylov vectors into overlap analysis without reshaping | Reshape to density matrix, ensure same basis as L was constructed in |
| Krylov gap vs `estimate_spectral_gap` (trajectory) | Expecting exact agreement | Trajectory estimate has O(delta) bias + statistical noise; agreement to ~10% is good |
| KrylovKit + existing Arpack dependency | Both active, different APIs, version conflicts | Consider replacing Arpack with KrylovKit for the dense path too (Arpack only needed in `furnace.jl`) |
| Dense diagnostics at n=6 for validation | Running `run_exact_diagnostics` at n=8+ (4096^2 dense matrix) | n=6 is the ceiling for dense; n=8+ is Krylov-only territory |

---

## Performance Traps

Patterns that work at small scale but fail as system size grows.

| Trap | Symptoms | Prevention | When It Breaks |
|------|----------|------------|----------------|
| Dense Kraus operators stored in memory | OOM at n=12 | Compute Kraus operators on-the-fly using NUFFT/OFT; only store jump.in_eigenbasis | n=12: 36 operators * 4096^2 * 16 bytes ~ 9 GB |
| krylovdim=100 at n=12 | 25 GB Krylov basis + scratch = OOM on most machines | Start with krylovdim=30, increase only with memory headroom | n=12 with krylovdim>50 exceeds 16 GB for Krylov basis alone |
| Full orthogonalization at n=12 | Each orthogonalization pass is O(krylovdim * dim^2) = O(50 * 16M) ~ 800M operations | Accept: this is intrinsic to Arnoldi. Reduce krylovdim if too slow. | n=12 with krylovdim=50: ~40 seconds per Krylov iteration |
| Matvec recomputing everything from scratch | Matvec takes minutes at n=12 (36 jumps * 2048 frequencies * dim^2 work) | Precompute and cache Kraus operators if memory allows; or precompute partial sums | n=10+: matvec time dominates total runtime |
| Trajectory validation at n=12 | 1000 trajectories * many steps * large dim = weeks of compute | Use Krylov for gap estimation, trajectory only for sanity checks at n<=8 | n=10+: trajectory approach becomes impractical for gap estimation |

---

## "Looks Done But Isn't" Checklist

Things that appear complete but are missing critical pieces.

- [ ] **Matvec implementation:** Often missing the coherent term (B operator). The existing `construct_lindbladian` adds it separately via `_precompute_coherent_total_B` and `_vectorize_liouvillian_coherent!`. The matrix-free path must include this term if `config.with_coherent == true`.
- [ ] **Convergence check:** `eigsolve` returns results even when `info.converged == 0`. Always check `info.converged >= howmany`.
- [ ] **Gap extraction:** Getting lambda_0 and lambda_1 right requires sorting by `abs(real(lambda))`, not by `abs(lambda)` or by `real(lambda)` (the latter would put the most negative eigenvalue first).
- [ ] **Eigenvector normalization:** Unit-norm vector != trace-1 density matrix. Must normalize by trace after reshaping.
- [ ] **TrotterDomain support:** BohrDomain works != TrotterDomain works. The basis selection, precomputation, and Kraus operator source differ.
- [ ] **Memory estimation:** "Works at n=8" does NOT mean "works at n=10". Compute memory budget explicitly: `krylovdim * 4^n * 16 * 1.5` bytes for the Krylov subspace alone.
- [ ] **Multiple eigenvalue request:** Requesting only `howmany=2` may miss the spectral gap if there are near-degenerate eigenvalues clustered near lambda_0. Request 4-10.

---

## Recovery Strategies

When pitfalls occur despite prevention, how to recover.

| Pitfall | Recovery Cost | Recovery Steps |
|---------|---------------|----------------|
| Wrong eigenvalues (Pitfall 1) | LOW | Change `which` parameter; increase `howmany`; validate against dense at n=4 |
| Memory blowup (Pitfall 2) | LOW | Kill process; reduce krylovdim; add memory estimation guard |
| Vectorization mismatch (Pitfall 3) | MEDIUM | Rewrite the matvec; requires re-running all validation tests |
| L vs E confusion (Pitfall 4) | MEDIUM | Rewrite the matvec to use the correct operator; rerun validation |
| Basis mismatch (Pitfall 5) | HIGH | Debug which basis is wrong; may require restructuring the matvec to match domain logic |
| Non-convergence (Pitfall 6) | LOW | Increase krylovdim and maxiter; try multiple starting vectors |
| Precision loss (Pitfall 7) | MEDIUM | Switch from (E-I)/delta to direct Lindblad formula; requires new matvec implementation |
| BLAS contention (Pitfall 8) | LOW | Set BLAS.set_num_threads() appropriately |
| Bad starting vector (Pitfall 9) | LOW | Change v0 to physically motivated choice |
| Wrong comparison (Pitfall 10) | LOW | Fix the comparison logic; eigenvalue sort order |

---

## Pitfall-to-Phase Mapping

How roadmap phases should address these pitfalls.

| Pitfall | Prevention Phase | Verification |
|---------|------------------|--------------|
| 1: Wrong eigenvalue targeting | Krylov API design | Test: Krylov gap matches dense gap at n=4 to 10 digits |
| 2: Memory blowup | Memory budget utilities | Test: `memory_estimate(12, 50)` returns < 512 GB; OOM guard in eigsolve wrapper |
| 3: Vectorization mismatch | Core matvec implementation | Test: `norm(L_dense * vec(rho) - vec(L_matvec(rho))) < 1e-12` for 10 random rho at n=4 |
| 4: L vs E confusion | Matvec design document | Test: explicit formula for L(rho) vs E(rho) in docstring; eigenvalue conversion documented |
| 5: Basis mismatch | Matvec implementation | Test: round-trip test passes for both BohrDomain and TrotterDomain at n=4 |
| 6: Non-convergence | Convergence tuning phase | Test: `info.converged >= howmany` asserted; adaptive krylovdim fallback |
| 7: Precision loss | Matvec validation | Test: relative precision of matvec quantified at n=4,6,8; tol set accordingly |
| 8: BLAS contention | Performance tuning | Test: eigsolve benchmarked with different BLAS thread counts; optimal choice documented |
| 9: Bad starting vector | API design | Test: default v0 = vec(rho_0 - rho_beta); convergence in <100 matvec applications at n=8 |
| 10: Wrong comparison | Validation framework | Test: `validate_krylov_vs_dense(n=4)` passes; comparison handles phase, ordering, degeneracy |
| 11: Convergence criteria | Gap estimation API | Test: gap reported with error bounds; kappa(lambda) computed |
| 12: Eigenvector normalization | Post-processing utilities | Test: `abs(tr(reshape(v, dim, dim))) - 1 > eps` triggers normalization |
| 13: Version incompatibility | Project.toml setup | KrylovKit version pinned in [compat]; CI tests pass against pinned version |
| 14: Hermiticity exploitation | Documentation | Docstring explicitly states: do NOT project onto Hermitian subspace during Arnoldi |

---

## Memory Budget Reference Table

Essential reference for phase planning. All values for ComplexF64 (16 bytes/element).

| n | dim | dim^2 | 1 vector | krylovdim=30 | krylovdim=50 | krylovdim=100 | Dense L matrix |
|---|-----|-------|----------|--------------|--------------|---------------|----------------|
| 4 | 16 | 256 | 4 KB | 120 KB | 200 KB | 400 KB | 1 MB |
| 6 | 64 | 4,096 | 64 KB | 1.9 MB | 3.2 MB | 6.4 MB | 256 MB |
| 8 | 256 | 65,536 | 1 MB | 30 MB | 50 MB | 100 MB | 64 GB |
| 10 | 1,024 | 1,048,576 | 16 MB | 480 MB | 800 MB | 1.6 GB | 16 TB |
| 12 | 4,096 | 16,777,216 | 256 MB | 7.5 GB | 12.5 GB | 25 GB | 4 PB |

**Implication:** Dense eigendecomposition is feasible only for n<=6 (256 MB matrix). Arpack shift-invert is feasible at n<=7 (requires dense matrix or sparse factorization). Krylov matrix-free is the ONLY option for n>=8.

**Kraus operator storage** (for the matvec):

| n | n_jumps (3 Paulis * n sites) | Operator size | Total Kraus storage |
|---|------------------------------|---------------|---------------------|
| 8 | 24 | 256x256 = 1 MB | 24 MB |
| 10 | 30 | 1024x1024 = 16 MB | 480 MB |
| 12 | 36 | 4096x4096 = 256 MB | 9 GB |

At n=12, storing all Kraus operators in memory alongside the Krylov subspace is feasible on the 512 GB node but tight. Total: 9 GB (Kraus) + 12.5 GB (Krylov, krylovdim=50) + 5 GB (scratch) ~ 27 GB.

---

## Sources

- [KrylovKit.jl Eigenvalue Problems Documentation](https://jutho.github.io/KrylovKit.jl/stable/man/eig/) -- eigsolve API, `which` parameters, convergence semantics
- [KrylovKit.jl Available Algorithms](https://jutho.github.io/KrylovKit.jl/stable/man/algorithms/) -- Arnoldi defaults: krylovdim=30, maxiter=100, tol=1e-12
- [KrylovKit.jl GitHub Issues #9, #23, #38](https://github.com/Jutho/KrylovKit.jl/issues/9) -- memory allocation concerns, degenerate eigenvalue issues
- [KrylovKit.jl Releases](https://github.com/Jutho/KrylovKit.jl/releases) -- v0.10 changes: BiArnoldi, BlockLanczos
- [Arpack.jl Standard Eigen Decomposition](https://julialinearalgebra.github.io/Arpack.jl/stable/eigs/) -- shift-invert sigma parameter
- [Arnoldi iteration failure modes (Embree, Virginia Tech)](https://personal.math.vt.edu/embree/66909.pdf) -- convergence proof gaps for non-Hermitian case
- [Implicitly Restarted Arnoldi Methods (STFC report RAL-TR-97-058)](https://epubs.stfc.ac.uk/manifestation/1345/RAL-TR-97-058.pdf) -- restart strategies
- [Tensor Network Framework for Lindbladian Spectra (arXiv:2509.07709)](https://arxiv.org/html/2509.07709) -- Krylov methods for Lindbladian spectral gap
- [Krylov complexity in open quantum systems (arXiv:2303.04175)](https://arxiv.org/html/2303.04175v2) -- bi-Lanczos for non-Hermitian Lindbladians
- [Julia BLAS threading discussion](https://discourse.julialang.org/t/julia-threads-vs-blas-threads/8914) -- thread contention patterns
- [ITensors.jl Multithreading guide](https://itensor.github.io/ITensors.jl/dev/Multithreading.html) -- BLAS thread management best practices
- [Three approaches for representing Lindblad dynamics (arXiv:1510.08634)](https://arxiv.org/pdf/1510.08634) -- vectorization conventions
- QuantumFurnace.jl codebase analysis: `furnace.jl` (shift-invert Arpack usage), `qi_tools.jl` (vectorization convention), `diagnostics.jl` (dense eigendecomposition), `jump_workers.jl` (per-domain Lindbladian construction), `trajectories.jl` (CPTP channel stepping), `test_threading.jl` (BLAS thread management)

---
*Pitfalls research for: Krylov-based Lindbladian spectral gap estimation in QuantumFurnace.jl*
*Researched: 2026-02-20*
