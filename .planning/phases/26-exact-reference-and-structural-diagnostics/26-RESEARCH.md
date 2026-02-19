# Phase 26: Exact Reference and Structural Diagnostics - Research

**Researched:** 2026-02-19
**Domain:** Lindbladian spectral analysis, KMS detailed balance diagnostics, symmetry sectors
**Confidence:** HIGH

## Summary

Phase 26 builds the exact diagnostic infrastructure that Phase 25 showed is critically needed: the n=6 Heisenberg chain has zero gap-mode overlap with all v1.3 preset observables (documented in Phase 25 VERIFICATION), confirming that symmetry sector analysis and a proper overlap formula with left/right eigenvectors are essential.

The codebase already has all core building blocks: `construct_lindbladian()` builds the full dense Lindbladian, `run_lindbladian()` extracts 2 eigenvalues via Arpack shift-invert, `eigenbasis_overlap_analysis()` does dense eigendecomposition with overlap computation (but uses a simplified `V\vec(rho0)` formula instead of left eigenvectors), and `log_sobolev.jl` demonstrates the sigma^{1/4} sandwich computation pattern needed for the KMS similarity transform. The phase needs to extend these into a proper multi-eigenvalue extraction, left+right eigenvector computation, correct overlap formula, anti-Hermitian defect computation, and symmetry sector labeling.

**Primary recommendation:** Build six standalone diagnostic functions (one per DIAG requirement) plus a bundling `run_exact_diagnostics()` that returns a single `ExactDiagnosticsResult` struct. Place all new code in a new `src/diagnostics.jl` file included after `gap_estimation.jl`. Replace the existing `build_preset_trajectory_observables` observable set with the new canonical set. Use dense `eigen()` for both left and right eigenvectors at n=4,6 (Arpack cannot compute left eigenvectors, and dense eigen is feasible at 4096x4096).

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

#### Defect ratio interpretation
- Advisory warning only -- report ||A||/lambda_gap(H) and flag if high, but do NOT gate or block downstream fitting
- Threshold value to be determined from supplementary information papers and numerical evidence at n=4,6 (Claude's discretion)
- KMS similarity transform uses a fixed internal spectrum truncation threshold (no user-facing parameter) -- chosen to be robust for our n=4,6 cases
- Report scalar ratio only (no full spectrum decomposition of anti-Hermitian part A) -- full eigen() infeasible above n=6
- Imaginary-to-real eigenvalue ratios |Im(lambda_k)/Re(lambda_k)| also reported for leading modes (per supplementary info Task 1.2)

#### Observable overlap scope
- Canonical observable set from supplementary info + Mz_stagg:
  - Z_1 (Pauli-Z on first site)
  - X_1 (Pauli-X on first site)
  - Z_1 * Z_{floor(n/2)} (two-point correlator)
  - H (the Hamiltonian itself)
  - Random traceless Hermitian matrix (control, fixed seed for reproducibility)
  - Mz_stagg (staggered magnetization, user addition)
- Replace existing v1.3 observable set (XX_avg, YY_avg, ZZ_avg, XZ_stagg, etc.) with this new canonical set in build_preset_observables
- Compute c_k for best-coupling observables across all three initial states:
  - |0>^n (all spins up)
  - |+>^n (all spins in X-plus state)
  - I/2^n (maximally mixed state)
- 20 leading modes for overlap coefficient computation
- c_k = Tr[O R_k] * Tr[L_k^dagger (rho_0 - rho_beta)] per supplementary info formula

#### Symmetry sector reporting
- General Sz-based analysis (works for any Hamiltonian that conserves total Sz: Heisenberg, Ising, XXZ, etc.)
- Delta_Sz labels assigned to each Lindbladian eigenvector based on density matrix support structure

#### Eigenvalue extraction defaults
- Default: 20 leading eigenvalues via Arpack shift-invert (nearest to zero)
- Store eigenvalues as Vector{ComplexF64} (natural representation, caller decomposes)
- Extract both left AND right eigenvectors upfront in DIAG-01 (both needed for overlap computation)
- API: standalone functions for each diagnostic PLUS a `run_exact_diagnostics()` bundle that calls them all and returns a single result struct

### Claude's Discretion
- Defect ratio warning threshold (determined from supplementary papers and numerical evidence)
- Near-degeneracy detection and multiplet grouping approach for symmetry sectors
- Handling of mixed-sector eigenvectors (dominant sector + purity fraction vs threshold-based labeling)
- Dashboard visualization choices (color-by-sector, marker shapes for the spectrum plot)
- Fixed internal truncation threshold for KMS similarity transform

### Deferred Ideas (OUT OF SCOPE)
- n=8 sparse Lindbladian diagonalization via KrylovKit.jl -- explicitly deferred to v2 (EVAL-01)
- Damped-oscillation fit model c*exp(-gamma*t)*cos(omega*t+phi) -- deferred to v2 (EVAL-02)
- GEVP / matrix pencil methods for multi-eigenvalue extraction -- deferred to v2 (EVAL-03)
</user_constraints>

## Standard Stack

### Core

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| LinearAlgebra (stdlib) | Julia 1.11+ | Dense `eigen()` for full eigendecomposition | Already used throughout codebase; needed for left/right eigenvectors |
| Arpack.jl | 0.5.4 | Shift-invert sparse eigenvalue extraction (20 leading) | Already in Project.toml; used by `run_lindbladian()` |
| LsqFit.jl | 0.15 | Not used in Phase 26 but downstream | Already in Project.toml |
| Random (stdlib) | Julia 1.11+ | Seed-controlled random traceless Hermitian observable | Already used throughout |

### Supporting

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| SparseArrays (stdlib) | Julia 1.11+ | Sparse Lindbladian construction (already used) | Already in stack |
| Printf (stdlib) | Julia 1.11+ | Diagnostic output formatting | Already in stack |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Dense `eigen()` for left eigenvectors | Arpack `eigs(L')` | Arpack cannot compute left eigenvectors directly; `eigs(L')` gives right eigenvectors of L^T which are conjugates of left eigenvectors of L, but normalizing biorthogonality is error-prone. Dense eigen at n=6 (4096x4096) takes ~seconds, well within budget. Use dense. |
| `opnorm(A)` for anti-Hermitian defect norm | Full eigendecomposition of A | User decision: scalar ratio only, no full spectrum of A. `opnorm()` computes the largest singular value (operator 2-norm) -- exactly what's needed. |

### No New Dependencies

No new packages needed. All required functionality is already in the dependency tree.

## Architecture Patterns

### Recommended File Structure

```
src/
  diagnostics.jl          # NEW: all Phase 26 diagnostic functions
  convergence.jl          # MODIFIED: replace observable set in build_preset_trajectory_observables
  gap_estimation.jl       # MODIFIED: update eigenbasis_overlap_analysis to use left+right eigenvectors
  QuantumFurnace.jl       # MODIFIED: include diagnostics.jl, update exports
```

### Pattern 1: Result Struct + Standalone Functions + Bundle

**What:** Each DIAG requirement maps to a standalone function returning its own result, plus a single `run_exact_diagnostics()` that calls all of them and returns a unified `ExactDiagnosticsResult`.

**When to use:** Always -- this is a locked decision from CONTEXT.md.

**Example:**
```julia
# Individual diagnostic functions
struct EigenDecompositionResult
    eigenvalues::Vector{ComplexF64}        # Leading 20, sorted by |Re(λ)|
    right_eigenvectors::Matrix{ComplexF64} # dim^2 x n_modes
    left_eigenvectors::Matrix{ComplexF64}  # dim^2 x n_modes (rows of L^{-1})
    spectral_gap::Float64                  # -Re(λ_2)
end

function extract_leading_eigendata(L::Matrix{ComplexF64}; n_modes::Int=20) -> EigenDecompositionResult
    # Dense eigen, sort by |Re(λ)|, extract leading n_modes
end

# Bundle function
function run_exact_diagnostics(L, hamiltonian, gibbs; ...) -> ExactDiagnosticsResult
    eigen_result = extract_leading_eigendata(L)
    fp_result = compute_fixed_point_distance(L, gibbs)
    defect_result = compute_anti_hermitian_defect(L, gibbs)
    # ...
    return ExactDiagnosticsResult(eigen_result, fp_result, defect_result, ...)
end
```

### Pattern 2: Eigendecomposition Strategy -- Dense for n<=6

**What:** Use dense `eigen(L)` for the full Lindbladian at n=4,6. This gives all 4^n eigenvalues and right eigenvectors. Left eigenvectors come from `eigen(L')` (transpose, not adjoint -- left eigenvectors of L are right eigenvectors of L^T). Then biorthonormalize.

**Why not Arpack for 20 modes:** Arpack returns right eigenvectors only. The overlap formula `c_k = Tr[O R_k] * Tr[L_k^dagger (rho_0 - rho_beta)]` requires both left (`L_k`) and right (`R_k`) eigenvectors. For n=4 (256x256) and n=6 (4096x4096), dense eigen is fast enough (~seconds). We take the leading 20 modes after sorting.

**Important subtlety:** For a non-normal matrix, the left and right eigenvectors from separate `eigen()` calls on L and L' need to be matched and biorthonormalized: `L_k^dagger R_j = delta_{kj}`. Julia's `eigen()` returns right eigenvectors. If `F = eigen(L)`, then `F.vectors` are right eigenvectors. For left eigenvectors: `G = eigen(transpose(L))` gives `G.vectors` where `G.vectors[:, k]` is the left eigenvector of L (since `L^T w = lambda w` implies `w^T L = lambda w^T`, so `w` is a left eigenvector of L). We biorthonormalize: `L_k_normalized = L_k / dot(L_k, R_k)`.

**Alternative (simpler):** Since we have the full dense eigen decomposition `F = eigen(L)` with `V = F.vectors`, we can compute left eigenvectors as rows of `V^{-1}` (i.e., `inv(V)`). This is the mathematically correct biorthogonal dual: if `L V = V D`, then `V^{-1} L = D V^{-1}`, so rows of `V^{-1}` are the left eigenvectors, already biorthonormalized with the right eigenvectors. This approach avoids a second eigendecomposition.

**Recommendation:** Use the `inv(V)` approach. It is simpler, avoids eigenvalue matching between two decompositions, and is automatically biorthonormalized. For n=6 (4096x4096), `inv(V)` is a single LAPACK call and takes ~seconds.

```julia
F = eigen(L)
perm = sortperm(abs.(real.(F.values)))
eigenvalues = F.values[perm]
V_right = F.vectors[:, perm]         # Right eigenvectors as columns
V_left = inv(V_right)'               # Left eigenvectors as columns (rows of V^{-1}, transposed)
# Verification: V_left' * V_right should be ≈ I
```

### Pattern 3: KMS Similarity Transform (from log_sobolev.jl)

**What:** The anti-Hermitian defect computation needs `D = rho^{-1/4} L[rho^{1/4} (.) rho^{1/4}] rho^{-1/4}` as a full matrix. The codebase already computes `sigma^{1/4}` in `log_sobolev.jl` (lines 29-31). Reuse this pattern.

**Implementation approach:** The similarity transform D is a superoperator. To build its matrix representation, apply it to each basis element of the vectorized operator space. For operator X (represented as a dim x dim matrix):

1. Compute `rho^{1/4}` and `rho^{-1/4}` from eigendecomposition of the Gibbs state (truncating eigenvalues below threshold for numerical stability of the inverse)
2. For each basis vector `e_j` of C^{dim^2}: reshape to dim x dim matrix `X_j`, compute `Y_j = rho^{-1/4} L[rho^{1/4} X_j rho^{1/4}] rho^{-1/4}`, vectorize back, store as column j of D_matrix
3. Actually, this is expensive (dim^2 matrix-vector products). More efficiently: if we have the full Lindbladian matrix L, then D_matrix = (rho^{-1/4} kron rho^{-1/4}) * L * (rho^{1/4} kron rho^{1/4}).

**Wait -- that's not quite right.** The similarity transform acts on the vectorized operator as:
```
D|X>> = |rho^{-1/4} L[rho^{1/4} X rho^{1/4}] rho^{-1/4}>>
```
In vectorized form, `|AXB>> = (B^T kron A)|X>>`. So:
- The inner sandwich `rho^{1/4} X rho^{1/4}` becomes `(rho^{1/4T} kron rho^{1/4})|X>>` = `(rho^{1/4} kron rho^{1/4})|X>>` since rho is Hermitian with real diagonal in eigenbasis
- Then L acts on this
- Then outer `rho^{-1/4} (.) rho^{-1/4}` becomes `(rho^{-1/4} kron rho^{-1/4})`

So: `D_matrix = kron(rho^{-1/4}, rho^{-1/4}) * L * kron(rho^{1/4}, rho^{1/4})`

**But** we must be careful about which basis we are in. The Lindbladian L is constructed in the Hamiltonian eigenbasis (or Trotter eigenbasis for TrotterDomain). The Gibbs state is diagonal in the Hamiltonian eigenbasis. So rho^{1/4} is diagonal with entries `(exp(-beta*E_i)/Z)^{1/4}`. The kron products are diagonal matrices when rho^{1/4} is diagonal, making this very efficient.

**Truncation for numerical stability:** For `rho^{-1/4}`, small eigenvalues of rho produce large entries. The user decided on a fixed internal threshold. Based on the supplementary info guidance and the physics (beta=10, n<=6), eigenvalues of the Gibbs state at the high-energy end are exp(-beta*E_max)/Z, which for the rescaled Heisenberg spectrum [0, 0.45] with beta=10 is exp(-4.5)/Z ~ very small but not catastrophically so. Recommendation: truncate eigenvalues below `eps_trunc = 1e-12` before taking the -1/4 power. This is conservative and will not affect the physics at n=4,6.

### Pattern 4: Symmetry Sector Labeling via Delta_Sz

**What:** For each Lindbladian eigenvector (reshaped as a dim x dim matrix), determine its dominant `Delta_Sz = Sz(E_i) - Sz(E_j)` quantum number based on which density matrix elements have the largest weight.

**Implementation:**

1. Compute total Sz eigenvalue for each Hamiltonian eigenstate: `Sz_vals[i] = <E_i|Sz_tot|E_i>` where `Sz_tot = sum_j Z_j / 2`. Since the Hamiltonian eigenstates are known, compute `Sz_tot` in the eigenbasis and take diagonal.

2. For each Lindbladian eigenvector R_k (reshaped as dim x dim matrix M_k), compute the weight in each (i,j) element: `|M_k[i,j]|^2`. The Delta_Sz of element (i,j) is `Sz_vals[i] - Sz_vals[j]`.

3. Group weights by Delta_Sz value. The dominant Delta_Sz is the one with the largest total weight.

4. Report: (dominant Delta_Sz, purity = fraction of weight in dominant sector)

**Near-degeneracy detection:** Recommendation: group Delta_Sz values that are within `1e-6` of each other (to handle floating point). For the Heisenberg chain, Sz eigenvalues are half-integers, so Delta_Sz values are integers -- no near-degeneracy issues from this. However, eigenvalues of L can be near-degenerate (multiplets). Group eigenvalues that are within `|lambda_i - lambda_j| / max(|lambda_i|, 1e-10) < 0.01` relative distance as a multiplet. Report multiplet structure.

### Pattern 5: Observable Construction -- New Canonical Set

**What:** Replace the 8-observable set in `build_preset_trajectory_observables` (H, Mz, XX_avg, YY_avg, ZZ_avg, Mz_stagg, Z1, XZ_stagg) with the new 6-observable canonical set (Z1, X1, Z1*Z_{n/2}, H, random traceless Hermitian, Mz_stagg).

**Impact analysis:** The function `build_preset_trajectory_observables` is called from:
- `estimate_spectral_gap()` in `gap_estimation.jl` (line 158)
- `run_trajectories_convergence()` callers pass their own observables
- `test/test_convergence.jl` (14 occurrences)
- `experiments/validate_spectral_gap.jl`

The tests reference specific observable names ("H", "Mz", "XX_avg", etc.) and check array sizes. These will need updating.

**Random traceless Hermitian construction:**
```julia
# Fixed seed for reproducibility across runs
rng = StableRNG(12345)  # or MersenneTwister(12345)
R = randn(rng, ComplexF64, dim, dim)
R = (R + R') / 2  # Hermitianize
R = R - tr(R) / dim * I(dim)  # Make traceless
R = R / opnorm(R)  # Normalize (operator norm = 1)
```

### Anti-Patterns to Avoid

- **Computing left eigenvectors via separate `eigen(L')`:** This requires matching eigenvalues between two decompositions (fragile for near-degenerate eigenvalues). Use `inv(V)` instead.
- **Using Arpack for left eigenvectors:** Arpack does not support left eigenvectors. Do not attempt `eigs(L'; ...)`.
- **Ignoring biorthonormalization:** The overlap formula requires biorthonormal left/right pairs. If using `inv(V)`, this is automatic. If using two separate eigen calls, explicit normalization is needed.
- **Applying rho^{-1/4} without truncation:** The Gibbs state has exponentially small eigenvalues at high energies. Without truncation, `rho^{-1/4}` overflows or produces numerically meaningless results.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Matrix fourth root | Custom Newton iteration | `eigvals^(1/4)` via eigendecomposition | Gibbs state is diagonal in eigenbasis; just raise diagonal entries to power |
| Operator norm of anti-Hermitian part | Manual singular value computation | `opnorm(A_matrix)` from LinearAlgebra | Standard LAPACK-backed function; returns largest singular value (= operator 2-norm) |
| Biorthogonal left/right eigenvectors | Two separate eigendecompositions + matching | `inv(V)` from single `eigen(L)` | Mathematically correct, avoids eigenvalue matching bugs |
| Dense Kronecker product for KMS transform | Manual loop over basis vectors | Since `rho^{1/4}` is diagonal in eigenbasis, use `Diagonal(d) kron Diagonal(d)` = `Diagonal(kron(d, d))` | Diagonal kron is O(dim^2) storage, not O(dim^4) |

**Key insight:** Because the Gibbs state is diagonal in the Hamiltonian eigenbasis, and the Lindbladian is constructed in this basis, the KMS similarity transform reduces to diagonal scaling: `D_matrix = Diagonal(d_inv) * L * Diagonal(d)` where `d = kron(rho_quarter_diag, rho_quarter_diag)` and `d_inv = kron(rho_inv_quarter_diag, rho_inv_quarter_diag)`. This is O(dim^4) multiplication (two diagonal matrix-matrix multiplies on a dim^2 x dim^2 matrix) but avoids forming a dense kron product.

## Common Pitfalls

### Pitfall 1: Left Eigenvector Normalization

**What goes wrong:** Computing overlap coefficients `c_k = Tr[O R_k] * Tr[L_k^dagger (rho_0 - rho_beta)]` with unnormalized left eigenvectors gives wrong magnitudes. The existing `eigenbasis_overlap_analysis` uses `alpha = V \ vec(rho0)` which is equivalent to `inv(V) * vec(rho0)`, so the left eigenvectors are implicitly `rows of inv(V)`. However, it then computes `c_k = dot(vec(O), V[:, k]) * alpha[k]` which is `Tr[O^dagger R_k] * (V^{-1} vec(rho0))_k`. For Hermitian O this gives `Tr[O R_k] * Tr[L_k^dagger rho_0]` -- close but not exactly the supplementary info formula which uses `(rho_0 - rho_beta)` not `rho_0`.

**Why it happens:** The existing code subtracts steady-state contribution by relying on `lambda_1 = 0` (so `exp(0*t) = 1` cancels between O(t) and O_ss). The supplementary info formula makes the subtraction explicit.

**How to avoid:** Use the exact formula: `c_k = Tr[O R_k] * Tr[L_k^dagger (rho_0 - rho_beta)]` for k >= 2. This properly handles the case where lambda_1 is not exactly 0 (numerical noise) and is consistent with the supplementary info.

**Warning signs:** `c_1` (steady-state mode) should be nearly zero when using `(rho_0 - rho_beta)`. If it's large, the subtraction is wrong.

### Pitfall 2: Basis Mismatch in KMS Transform

**What goes wrong:** The Lindbladian is in the Hamiltonian eigenbasis (for BohrDomain/TimeDomain) or Trotter eigenbasis (for TrotterDomain). The Gibbs state `hamiltonian.gibbs` is in the Hamiltonian eigenbasis. If working in TrotterDomain, a basis change is needed.

**Why it happens:** The codebase has two paths: Hamiltonian eigenbasis and Trotter eigenbasis. The Phase 26 context says "Fixed point trace distance should be ~0 for BohrDomain, ~1e-8 for TrotterDomain."

**How to avoid:** For Phase 26, focus on BohrDomain/TimeDomain where L and rho_beta are in the same basis. If TrotterDomain is needed, use `_gibbs_in_trotter_basis()` from convergence.jl. Document which basis each matrix is in.

**Warning signs:** Trace distance between fixed point and Gibbs state is unexpectedly large.

### Pitfall 3: Vectorization Convention (Row-Major vs Column-Major)

**What goes wrong:** Julia uses column-major storage. `vec(M)` stacks columns. The Lindbladian `L` is built using `kron(A, conj(B))` in the Watrous convention (seen in `qi_tools.jl`). When reshaping eigenvectors back to density matrix form, `reshape(v, dim, dim)` must use the same convention as the Lindbladian construction.

**Why it happens:** Different references use different vectorization conventions. The codebase uses the `|A kron B^* |X>>` convention consistently (Watrous convention).

**How to avoid:** Always use `reshape(v, dim, dim)` to convert eigenvectors back to operators. Test by verifying that the fixed-point eigenvector, reshaped and normalized, matches the Gibbs state.

**Warning signs:** The reshaped fixed-point density matrix is not Hermitian or has wrong trace.

### Pitfall 4: Degenerate Eigenvalues in Symmetry Sector Labeling

**What goes wrong:** For the isotropic Heisenberg chain, SU(2) symmetry creates multiplets of degenerate Lindbladian eigenvalues. The eigenvectors within a degenerate subspace are not uniquely determined -- any linear combination is valid. This means a single eigenvector might span multiple Delta_Sz sectors.

**Why it happens:** `eigen()` (LAPACK) returns an arbitrary basis within each degenerate subspace. Two runs with slightly different inputs may give different eigenvectors for the same eigenvalue.

**How to avoid:** For near-degenerate eigenvalues (within relative tolerance 1%), group them as a multiplet. Report the Delta_Sz distribution across the entire multiplet subspace rather than per-vector. Within a multiplet, compute the weight-averaged Delta_Sz profile.

**Warning signs:** Eigenvectors with low "purity" (less than 80% weight in a single sector). These are likely members of a multiplet with an arbitrary basis choice.

### Pitfall 5: Defect Ratio Interpretation at Finite Temperature

**What goes wrong:** The supplementary info says: if `||A||/lambda_gap(H) < 0.01`, the Lindbladian is effectively normal. But `lambda_gap(H)` here refers to the gap of the *Hermitian part* of the similarity-transformed Lindbladian, not the Hamiltonian spectral gap. These are different quantities.

**Why it happens:** Notation overload between `lambda_gap(H)` (Hamiltonian gap) and `lambda_gap(H_D)` (Hermitian part of similarity-transformed Lindbladian).

**How to avoid:** Compute `lambda_gap(Hermitian_part_of_D)` explicitly. The anti-Hermitian defect ratio is `||A|| / lambda_gap(H_D)` where `H_D = (D + D')/2`. This is what the supplementary info Task 1.2 specifies.

**Warning signs:** Defect ratio wildly different from what the eigenvalue Im/Re ratios suggest.

## Code Examples

### DIAG-01: Leading Eigenvalue Extraction

```julia
# Source: Codebase pattern from furnace.jl run_lindbladian() + extension
function extract_leading_eigendata(L::Matrix{ComplexF64}; n_modes::Int=20)
    # Dense eigendecomposition (feasible for n<=6: 4096x4096)
    F = eigen(L)

    # Sort by proximity to zero (|Re(lambda)|)
    perm = sortperm(abs.(real.(F.values)))
    eigenvalues = F.values[perm[1:n_modes]]
    V_right = F.vectors[:, perm[1:n_modes]]

    # Left eigenvectors: rows of V^{-1} for the selected modes
    # Full inverse then select rows -- or use the full V and invert
    V_full = F.vectors[:, perm]
    V_inv = inv(V_full)  # Left eigenvectors as rows
    V_left_conj = V_inv[1:n_modes, :]'  # Columns = left eigenvectors (conjugated for biorthogonality)

    spectral_gap = abs(real(eigenvalues[2]))

    return EigenDecompositionResult(eigenvalues, V_right, V_left_conj, spectral_gap)
end
```

### DIAG-02: Fixed Point Distance

```julia
# Source: Pattern from furnace.jl run_lindbladian() lines 24-28
function compute_fixed_point_distance(eigen_result, gibbs::Hermitian)
    # Fixed point = eigenvector for lambda_1 ≈ 0
    fp_vec = eigen_result.right_eigenvectors[:, 1]
    dim = isqrt(length(fp_vec))
    fp_dm = reshape(fp_vec, dim, dim)

    # Normalize: make Hermitian, trace = 1
    hermitianize!(fp_dm)
    fp_dm ./= tr(fp_dm)

    dist = trace_distance_h(Hermitian(fp_dm), gibbs)
    return (fixed_point=fp_dm, trace_distance=dist)
end
```

### DIAG-03/04: Anti-Hermitian Defect via KMS Transform

```julia
# Source: Supplementary info Task 1.2 + log_sobolev.jl sigma^{1/4} pattern
function compute_anti_hermitian_defect(L::Matrix{ComplexF64}, gibbs::Hermitian;
                                        eps_trunc::Float64=1e-12)
    dim = size(gibbs, 1)
    dim2 = dim^2

    # Compute rho^{1/4} and rho^{-1/4} diagonal entries (Gibbs is diagonal in eigenbasis)
    gibbs_diag = real.(diag(Matrix(gibbs)))
    gibbs_diag_safe = max.(gibbs_diag, eps_trunc)  # Truncate for stability

    rho_quarter = gibbs_diag_safe .^ 0.25
    rho_inv_quarter = gibbs_diag_safe .^ (-0.25)

    # KMS similarity transform: D = kron(rho^{-1/4}, rho^{-1/4}) * L * kron(rho^{1/4}, rho^{1/4})
    # Since rho^{1/4} is diagonal, kron of diagonals is diagonal
    d_right = kron(rho_quarter, rho_quarter)       # Vector of length dim^2
    d_left  = kron(rho_inv_quarter, rho_inv_quarter)

    # D_ij = d_left[i] * L_ij * d_right[j]
    D_matrix = d_left .* L .* d_right'  # Broadcasting: d_left is column, d_right' is row

    # Hermitian/anti-Hermitian decomposition
    H_part = (D_matrix + D_matrix') / 2
    A_part = (D_matrix - D_matrix') / 2

    # Metrics
    A_norm = opnorm(A_part)  # Operator 2-norm (largest singular value)

    # Gap of Hermitian part
    H_eigenvalues = eigvals(Hermitian(H_part))
    sorted_H = sort(H_eigenvalues, by=x->abs(x))
    H_gap = abs(sorted_H[2])  # Second smallest |eigenvalue|

    defect_ratio = A_norm / H_gap

    return (A_norm=A_norm, H_gap=H_gap, defect_ratio=defect_ratio, D_matrix=D_matrix)
end
```

### DIAG-05: Observable Overlap Coefficients (Correct Formula)

```julia
# Source: Supplementary info Task 1.3
function compute_overlap_coefficients(
    eigen_result,  # from DIAG-01
    observables::Vector{<:Matrix{<:Complex}},
    rho0::Matrix{<:Complex},
    rho_beta::Hermitian;
    n_modes::Int=20
)
    dim = size(rho0, 1)
    rho_diff = rho0 - Matrix(rho_beta)  # rho_0 - rho_beta

    n_obs = length(observables)
    coeffs = zeros(ComplexF64, n_obs, n_modes)

    for k in 1:n_modes
        R_k = reshape(eigen_result.right_eigenvectors[:, k], dim, dim)
        L_k = reshape(eigen_result.left_eigenvectors[:, k], dim, dim)

        # c_k = Tr[O R_k] * Tr[L_k^dagger (rho_0 - rho_beta)]
        lk_factor = tr(L_k' * rho_diff)

        for (i, O) in enumerate(observables)
            ok_factor = tr(O * R_k)
            coeffs[i, k] = ok_factor * lk_factor
        end
    end

    return coeffs  # n_obs x n_modes
end
```

### DIAG-06: Delta_Sz Symmetry Labels

```julia
# Source: Supplementary info Task 4.1
function compute_sz_labels(eigen_result, hamiltonian::HamHam; n_modes::Int=20)
    dim = size(hamiltonian.data, 1)
    n = Int(log2(dim))

    # Total Sz operator in eigenbasis
    Sz_comp = zeros(ComplexF64, dim, dim)
    for site in 1:n
        Sz_comp .+= Matrix{ComplexF64}(pad_term([Z], n, site)) / 2
    end
    V = hamiltonian.eigvecs
    Sz_eigen = V' * Sz_comp * V
    sz_vals = real.(diag(Sz_eigen))  # Sz eigenvalue for each energy eigenstate

    labels = Vector{NamedTuple}(undef, n_modes)
    for k in 1:n_modes
        M_k = reshape(eigen_result.right_eigenvectors[:, k], dim, dim)
        weights = abs2.(M_k)  # |M_k[i,j]|^2

        # Compute Delta_Sz for each (i,j) entry
        delta_sz_map = Dict{Float64, Float64}()  # Delta_Sz => total weight
        for j in 1:dim, i in 1:dim
            w = weights[i, j]
            w < 1e-14 && continue
            dsz = round(sz_vals[i] - sz_vals[j]; digits=6)
            delta_sz_map[dsz] = get(delta_sz_map, dsz, 0.0) + w
        end

        total_weight = sum(values(delta_sz_map))
        dominant_dsz, dominant_weight = findmax(delta_sz_map)
        purity = dominant_weight / max(total_weight, 1e-30)

        labels[k] = (delta_sz=dominant_dsz, purity=purity, sector_weights=delta_sz_map)
    end

    return labels
end
```

## State of the Art

| Old Approach (v1.3) | Current Approach (Phase 26) | When Changed | Impact |
|--------|---------|--------|--------|
| Dense `eigen()` + `V\rho0` for overlap (no left eigenvectors) | Dense `eigen()` + `inv(V)` for left eigenvectors + correct formula | Phase 26 | Fixes zero-overlap mystery at n=6; enables proper c_k computation |
| 8 observables (H, Mz, XX_avg, YY_avg, ZZ_avg, Mz_stagg, Z1, XZ_stagg) | 6 canonical observables (Z1, X1, Z1*Z_{n/2}, H, random traceless, Mz_stagg) | Phase 26 | Matches supplementary info prescription; adds X1 and two-point correlator |
| Only right eigenvectors, no symmetry info | Right + left eigenvectors, Delta_Sz sector labels | Phase 26 | Explains which modes are accessible from which initial states/observables |
| No anti-Hermitian defect analysis | Full KMS similarity transform + defect ratio | Phase 26 | Determines if real exponential fits are appropriate |
| 2 eigenvalues from Arpack | 20 leading eigenvalues (dense eigen then select) | Phase 26 | Reveals multiplet structure, enables multi-mode overlap analysis |

## Discretion Recommendations

### Defect Ratio Warning Threshold

**Recommendation: 0.1**

From the supplementary info, Chen et al. (2023) Proposition II.3 states that if `||A||/lambda_gap(H_D) <= 1/100`, the Hermitian gap controls the decay. The supplementary info also says if the ratio is `>= 0.1`, non-normality effects cause oscillatory transients. Use `0.1` as the advisory warning threshold. Below 0.1 = "real exponential fitting appropriate"; above 0.1 = "warning: consider oscillatory model". This aligns with the `|Im(lambda_k)/Re(lambda_k)| > 0.1` threshold also mentioned.

**Confidence: MEDIUM** -- Based on supplementary info guidance, not independent verification. Should be validated empirically at n=4,6.

### KMS Transform Truncation Threshold

**Recommendation: 1e-12**

For the rescaled Heisenberg spectrum [0, 0.45] at beta=10, the smallest Gibbs eigenvalue is approximately exp(-10 * 0.45) / Z ~ exp(-4.5) / Z. With dim=64 (n=6), Z ~ sum of exp(-10*E_i), this is of order 1e-3 to 1e-4. So rho^{-1/4} involves raising 1e-4 to the -0.25 power ~ 5.6, which is perfectly well-behaved. The 1e-12 truncation is extremely conservative and will never trigger at n=4,6. It serves as a safety net only.

**Confidence: HIGH** -- Direct computation from known spectrum bounds.

### Near-Degeneracy Detection

**Recommendation:** Two eigenvalues lambda_i, lambda_j are in the same multiplet if `|lambda_i - lambda_j| / max(|lambda_i|, |lambda_j|, 1e-10) < 0.01` (1% relative distance). Within a multiplet, report the collective Delta_Sz distribution.

**Confidence: MEDIUM** -- The 1% threshold is a heuristic. For the Heisenberg chain, SU(2) multiplets should have exactly degenerate eigenvalues (up to numerical noise), so even a much tighter threshold would work. The 1% is conservative to handle near-degeneracies from approximate symmetries.

### Mixed-Sector Eigenvector Handling

**Recommendation:** Report dominant sector + purity fraction. An eigenvector is "pure" if purity > 0.95 (95% weight in one Delta_Sz sector). Otherwise, label it as "mixed" and report the top-2 sectors with their weights. This is more informative than a hard threshold.

**Confidence: MEDIUM** -- The 95% threshold is a heuristic but physically motivated: exact symmetry eigenstates would have 100% purity.

### Dashboard Visualization

**Recommendation:** For the spectrum plot (Panel A of the supplementary info dashboard):
- X-axis: Re(lambda_k), Y-axis: Im(lambda_k)
- Color by Delta_Sz sector (using a categorical colormap)
- Shape: circles for pure-sector modes (purity > 0.95), diamonds for mixed-sector modes
- Size proportional to max overlap coefficient |c_k| across all observables

This directly addresses the n=6 mystery: the plot will visually show that the gap mode lives in a sector that no observable can reach.

**Confidence: MEDIUM** -- Visualization choices; could be refined.

## Open Questions

1. **Arpack vs Dense Eigen for DIAG-01**
   - What we know: The locked decision says "20 leading eigenvalues via Arpack shift-invert", but left eigenvectors require dense eigen. Arpack cannot compute left eigenvectors.
   - What's unclear: Should we use Arpack for just the eigenvalues and dense eigen for the full decomposition needed by DIAG-05? Or just use dense eigen for everything?
   - Recommendation: Use dense `eigen()` for everything at n=4,6. The CONTEXT says "Arpack shift-invert" as the method, but the actual constraint is "both left AND right eigenvectors upfront in DIAG-01". Since Arpack cannot provide left eigenvectors, dense eigen is the only feasible approach. The Arpack call in the existing `run_lindbladian()` extracts only 2 modes and can remain as-is for backward compatibility. The new diagnostic function does a separate dense decomposition. This is consistent with the existing `eigenbasis_overlap_analysis()` which already uses dense `eigen()`.

2. **Observable Set Replacement Impact on Tests**
   - What we know: 14 test occurrences reference the old observable set by name. Changing `build_preset_trajectory_observables` will break them.
   - What's unclear: Should the observable builder be renamed or versioned? Should old tests be updated in Phase 26 or a separate task?
   - Recommendation: Update `build_preset_trajectory_observables` in-place (same function, new observable set) and update all tests in the same plan. The function contract (returns observables + names) stays the same; only the specific observables change. Tests that check specific names ("XX_avg") must be updated. This is a locked decision from CONTEXT.md.

3. **Initial State for Random Traceless Hermitian**
   - What we know: Need a reproducible random traceless Hermitian matrix as a control observable.
   - What's unclear: Which RNG seed? Using `StableRNGs.StableRNG(12345)` would be ideal for cross-platform reproducibility, but StableRNGs is only in test dependencies. `MersenneTwister(12345)` from stdlib `Random` is available in the main module.
   - Recommendation: Use `Random.MersenneTwister(12345)` for the random observable construction. Julia's MersenneTwister is deterministic for a given seed within the same Julia version, which suffices for reproducibility. If cross-version reproducibility becomes important, move StableRNGs to main deps.

## Sources

### Primary (HIGH confidence)
- Codebase inspection: `src/furnace.jl` (run_lindbladian, construct_lindbladian), `src/gap_estimation.jl` (eigenbasis_overlap_analysis), `src/log_sobolev.jl` (sigma^{1/4} pattern), `src/qi_tools.jl` (vectorization helpers), `src/convergence.jl` (build_preset_trajectory_observables), `src/structs.jl` (data types)
- Phase 25 VERIFICATION: Documents n=6 zero-overlap problem as genuine physics, confirms ARPACK works at n=4
- Supplementary information: `spectral-gap-refinements-instructions.md` (Task 1.1-1.3, 4.1 formulas)
- Supplementary information: `error_catalogue_spectral_gap_estimation.md` (Error sources 1-7, overlap formula)

### Secondary (MEDIUM confidence)
- [Arpack.jl documentation](https://arpack.julialinearalgebra.org/stable/eigs/) -- confirmed `eigs` does not support left eigenvectors
- [Julia LinearAlgebra docs](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/) -- `eigen()`, `opnorm()`, `eigvals()` API

### Tertiary (LOW confidence)
- Defect ratio threshold 0.1: derived from supplementary info interpretation of Chen et al. (2023) Proposition II.3, needs empirical validation at n=4,6

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- no new dependencies, all tools already in Project.toml
- Architecture: HIGH -- clear patterns from existing codebase (log_sobolev sigma^{1/4}, eigenbasis_overlap_analysis dense eigen, furnace.jl Lindbladian construction)
- Pitfalls: HIGH -- based on direct codebase inspection and Phase 25 documented issues
- Discretion areas: MEDIUM -- warning thresholds need empirical validation, visualization choices are subjective

**Research date:** 2026-02-19
**Valid until:** 2026-03-19 (stable domain; no fast-moving dependencies)
