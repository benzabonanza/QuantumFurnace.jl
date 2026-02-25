# Phase 32: Krylov Simulator Speedup - Research

**Researched:** 2026-02-25
**Domain:** Julia BLAS performance optimization for Lindbladian matvec
**Confidence:** HIGH

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

#### Precomputed effective Hamiltonian for Lindbladian matvec
- Currently each dissipator term computes L'L (1 GEMM) and anticommutator (2 GEMMs) per term, totaling 5N GEMMs for N = jumps x energy_labels
- Since Krylov always applies the FULL Lindbladian (unlike trajectory which samples individual jumps), precompute R_total = sum_i scalar_i * L_i'L_i at workspace construction
- Define G_left = i*B^T - 0.5*R_total^T and G_right = -i*B^T - 0.5*R_total^T (the non-Hermitian effective Hamiltonian in kron convention)
- Per-matvec becomes: rho_out = G_left*rho + rho*G_right + sandwiches (2 GEMMs constant + 2 GEMMs per sandwich term)
- This reduces from 5N GEMMs to 2 + 2N GEMMs -- roughly 2.5x fewer BLAS calls
- R_total accumulation already exists (`_accumulate_R_total!` in krylov_workspace.jl) for the channel path -- reuse for Lindbladian path
- Apply to all domains: EnergyDomain, TimeDomain, TrotterDomain, BohrDomain
- Similarly precompute for adjoint: adjoint just swaps G_left and G_right

#### Delete legacy Euler apply_delta_channel!
- Remove the 5-argument `apply_delta_channel!(ws, rho, delta, config_liouv, hamiltonian)` Euler approximation
- It is faulty: doesn't retain correct O(delta^2) error properties of Chen's CPTP algorithm
- Also delete the corresponding test: "apply_delta_channel! legacy Euler" testset in test_krylov_eigsolve.jl

#### Energy label count is small (truncated)
- Energy labels are truncated from 2^num_energy_bits to just a few dozen relevant labels
- The optimization still matters: saving 3 GEMMs per iteration x ~30 iterations x ~100 matvecs

### Claude's Discretion
- Whether to store G_left and G_right as separate fields or compute G_right = -conj(G_left) on the fly
- Exact scratch matrix reuse strategy for the sandwich-only loop
- Whether BohrDomain needs a separate precomputed matrix or can share the same R_total pattern

### Deferred Ideas (OUT OF SCOPE)
None -- discussion stayed within phase scope
</user_constraints>

## Summary

Phase 32 is a performance optimization that restructures the Lindbladian matvec hot path to precompute aggregate matrices at workspace construction time, reducing the per-matvec GEMM count from 5N to 2+2N (where N = jumps x energy_labels). The optimization exploits the fact that the Krylov path always applies the FULL Lindbladian (summing over all jumps and frequencies), unlike the trajectory path which samples individual jumps. The L'L product and anticommutator terms can therefore be precomputed and fused into a single "effective Hamiltonian" matrix G_left/G_right applied once per matvec, leaving only the sandwich terms (which depend on rho) in the per-term inner loop.

The second deliverable is dead code removal: deleting the legacy 5-argument Euler `apply_delta_channel!` which was superseded by the faithful Chen CPTP channel in Quick-36 but retained for backward compatibility.

The codebase is well-structured for this change. The `_accumulate_R_total!` helper already exists for all 4 domains (EnergyDomain, TimeDomain, TrotterDomain, BohrDomain) and is currently used only by the ThermalizeConfig workspace constructor for the channel path. Reusing it for the Lindbladian path is straightforward.

**Primary recommendation:** Add `G_left` and `G_right` fields to `KrylovWorkspace`, compute them once at construction by calling the existing `_accumulate_R_total!`, then refactor all 8 `apply_lindbladian!`/`apply_adjoint_lindbladian!` functions to use a 2-GEMM effective Hamiltonian term + sandwich-only inner loop.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| LinearAlgebra.BLAS | Julia stdlib | `gemm!`, `axpy!` for zero-allocation GEMM | Already used throughout the hot path; the only way to avoid mul! wrapper allocations |
| KrylovKit | >=0.7 | Arnoldi eigsolve consuming the optimized matvec | Existing dependency; not changed in this phase |

### Supporting
No new dependencies. This phase modifies only existing files.

## Architecture Patterns

### Current File Structure (affected files)
```
src/
├── krylov_workspace.jl    # KrylovWorkspace struct + constructors + _accumulate_R_total!
├── krylov_matvec.jl       # apply_lindbladian!, apply_adjoint_lindbladian!, _accumulate_dissipator!
└── krylov_eigsolve.jl     # apply_delta_channel! (legacy Euler to delete), krylov_spectral_gap
test/
├── test_krylov_matvec.jl  # Round-trip correctness + allocation tests (must keep passing)
└── test_krylov_eigsolve.jl # Legacy Euler testset to delete
```

### Pattern 1: Current Dissipator Per-Term GEMM Cost (to be optimized)

**What:** Each call to `_accumulate_dissipator!(out, L_op, rho, scalar, ws)` performs 5 GEMMs:
1. `L'L -> ws.LdagL` (1 GEMM: `gemm!('C','N',...)`)
2. `(L'L)^T * rho -> ws.tmp1` (1 GEMM: `gemm!('T','N',...)`)
3. `rho * (L'L)^T -> ws.tmp1` (1 GEMM: `gemm!('N','T',...)`)
4. `conj(L) * rho -> ws.tmp1` (1 GEMM: `gemm!('N','N',...)`)
5. `(conj(L)*rho) * L^T -> ws.LdagL` (1 GEMM: `gemm!('N','T',...)`)

Plus 3 `axpy!` calls for accumulation. For N dissipator terms, this is 5N GEMMs total.

**After optimization:** The per-term inner loop only computes the sandwich:
1. `conj(L) * rho -> ws.tmp1` (1 GEMM)
2. `(conj(L)*rho) * L^T -> ws.LdagL` (1 GEMM)

Plus 1 `axpy!`. The anticommutator is absorbed into the precomputed G_left/G_right:
1. `G_left * rho -> ws.tmp1` (1 GEMM, computed once)
2. `rho * G_right -> ws.tmp1` (1 GEMM, computed once)

Total: 2 + 2N GEMMs (down from 5N).

### Pattern 2: G_left / G_right Precomputation

**What:** At workspace construction, compute:
```
R_total = sum_i scalar_i * L_i'L_i     (already exists via _accumulate_R_total!)
G_left  = i*B^T - 0.5*R_total^T        (non-Hermitian effective Hamiltonian, left action)
G_right = -i*B^T - 0.5*R_total^T       (right action, note: G_right = -conj(G_left) only if B and R_total are Hermitian)
```

The Lindbladian action becomes:
```
L(rho) = G_left * rho + rho * G_right + sum_i scalar_i * conj(L_i) * rho * L_i^T
```

For the adjoint:
```
L*(rho) = G_right * rho + rho * G_left + sum_i scalar_i * L_i^T * rho * conj(L_i)
```
i.e., adjoint just swaps G_left and G_right (and swaps the sandwich direction).

**Source:** Derived from the existing code in `krylov_matvec.jl`:
- Coherent term (lines 228-235): `i*B^T*rho - i*rho*B^T`
- Anticommutator (lines 50-56 in `_accumulate_dissipator!`): `-0.5*(L'L)^T*rho - 0.5*rho*(L'L)^T`
- These two combine linearly across all terms into G_left * rho + rho * G_right.

### Pattern 3: KrylovWorkspace Extension

**What:** Add two new fields to the struct for the precomputed effective Hamiltonian matrices.

Current struct has 16 fields. Adding G_left and G_right brings it to 18. The fields are `Union{Nothing, Matrix{T}}` to handle the case where both B_total and R_total might be zero (though this never happens in practice for non-trivial systems).

```julia
struct KrylovWorkspace{T<:Complex, PD<:NamedTuple}
    # ... existing fields ...

    # Precomputed effective Hamiltonian for Lindbladian matvec (Phase 32)
    G_left::Union{Nothing, Matrix{T}}     # i*B^T - 0.5*R_total^T
    G_right::Union{Nothing, Matrix{T}}    # -i*B^T - 0.5*R_total^T
end
```

### Pattern 4: Sandwich-Only Inner Loop

**What:** After extracting the anticommutator into G_left/G_right, the inner loop over energy labels only computes the sandwich term. This is a simplified version of `_accumulate_dissipator!`.

For EnergyDomain/TimeDomain/TrotterDomain (1-op dissipator):
```julia
function _accumulate_sandwich!(out, L_op, rho, scalar, ws)
    CT = one(eltype(out))
    ZT = zero(eltype(out))
    # conj(L) * rho -> tmp1
    @. ws.tmp2 = conj(L_op)
    BLAS.gemm!('N', 'N', CT, ws.tmp2, rho, ZT, ws.tmp1)
    # (conj(L)*rho) * L^T -> LdagL (or tmp2, since we can reuse)
    BLAS.gemm!('N', 'T', CT, ws.tmp1, L_op, ZT, ws.LdagL)
    BLAS.axpy!(eltype(out)(scalar), ws.LdagL, out)
    return nothing
end
```

For the Hermitian half-grid negative-frequency partner (L = A(w)', sandwich = L^T * rho * conj(L)):
```julia
function _accumulate_sandwich_adj_L!(out, L_op, rho, scalar, ws)
    CT = one(eltype(out))
    ZT = zero(eltype(out))
    @. ws.tmp2 = conj(L_op)
    BLAS.gemm!('T', 'N', CT, L_op, rho, ZT, ws.tmp1)    # L^T * rho
    BLAS.gemm!('N', 'N', CT, ws.tmp1, ws.tmp2, ZT, ws.LdagL)  # L^T * rho * conj(L)
    BLAS.axpy!(eltype(out)(scalar), ws.LdagL, out)
    return nothing
end
```

### Pattern 5: BohrDomain Two-Operator Dissipator

**What:** BohrDomain uses a different dissipator structure with two distinct operators A and B_dag. The anticommutator base product is `B_dag * A` (not `L'L`), and R_total = sum gamma_norm_factor * (B_dag * A).

The existing `_accumulate_R_total!` for BohrDomain (krylov_workspace.jl lines 282-314) already computes this correctly:
```julia
mul!(R, A_nu2_dag, jump_oft, gamma_norm_factor, 1.0)  # R += gamma_norm_factor * B_dag * A
```

So G_left/G_right precomputation works the same way: R_total accumulates the anticommutator base, and the per-matvec sandwich uses only `_accumulate_dissipator_2op!` with the anticommutator terms removed.

The BohrDomain sandwich-only function:
```julia
function _accumulate_sandwich_2op!(out, A, B_dag, rho, scalar, ws)
    CT = one(eltype(out))
    ZT = zero(eltype(out))
    # B_dag^T * rho -> tmp1
    BLAS.gemm!('T', 'N', CT, B_dag, rho, ZT, ws.tmp1)
    # B_dag^T * rho * A^T -> tmp2
    BLAS.gemm!('N', 'T', CT, ws.tmp1, A, ZT, ws.tmp2)
    BLAS.axpy!(eltype(out)(scalar), ws.tmp2, out)
    return nothing
end
```

### Anti-Patterns to Avoid

- **Storing only G_left and computing G_right = -conj(G_left) on the fly:** This identity only holds when both B_total and R_total are real or Hermitian. B_total is Hermitian (by construction), and R_total is Hermitianized after accumulation. In kron convention with the transpose, `-conj(i*B^T - 0.5*R_total^T) = -(-i*conj(B^T)) + 0.5*conj(R_total^T) = i*conj(B)^T + 0.5*conj(R_total)^T`. Since B is Hermitian, `conj(B)^T = B` (NOT `B^T`), and similarly for R_total. So `conj(G_left) = -i*B + 0.5*R_total` which equals `-G_right` only if we compute the transpose again. The relationship is: `G_right = -transpose(conj(G_left))` = `-G_left'`. **Recommendation: Store both G_left and G_right as separate fields.** The memory cost is negligible (two dim x dim matrices) and avoids subtle conjugation errors. Computing `G_right = -adjoint(G_left)` at construction time is cleaner than deriving it on the fly.

- **Modifying the sandwich scratch buffer strategy carelessly:** The current `_accumulate_dissipator!` uses `ws.LdagL` as scratch for the L'L product, then repurposes it for the sandwich result (since L'L is no longer needed after the anticommutator). After the optimization, `ws.LdagL` is fully available for scratch in the sandwich-only loop. **The existing scratch layout (tmp1, tmp2, LdagL) is sufficient.** No new scratch matrices needed.

- **Attempting to fuse BohrDomain R_total with the A_nu2_dag allocation:** The BohrDomain `apply_lindbladian!` currently allocates `A_nu2_dag = zeros(T, dim, dim)` per matvec call. This allocation exists in the current code and is "acceptable for Bohr" per the existing comment. Phase 32 should not try to fix this -- it's a separate concern.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| R_total accumulation | New accumulation code | Existing `_accumulate_R_total!` (3 domain methods) | Already handles Hermitian half-grid, NUFFT prefactors, Bohr buckets correctly |
| Hermitianization | Manual symmetrize | Existing `hermitianize!` | Already used in workspace constructor for channel R_total |
| Dense matrix exponential | Custom exp(-i*delta*B) | `exp(Hermitian(...))` | Already used for U_coherent in ThermalizeConfig constructor |

**Key insight:** The `_accumulate_R_total!` functions already exist in `krylov_workspace.jl` and handle all 4 domains. They were written for the ThermalizeConfig channel path (Quick-36). The Lindbladian path now simply reuses the same functions at construction time.

## Common Pitfalls

### Pitfall 1: Kron Convention Transpose in G_left/G_right
**What goes wrong:** The kron convention used throughout is `kron(A,B)vec(X) = vec(B*X*A^T)`. This means the anticommutator `{(L'L)^T, rho}` un-vectorizes to `(L'L)^T * rho + rho * (L'L)^T` -- note the transposes. So R_total^T (not R_total) appears in G_left/G_right.
**Why it happens:** Mixing up `R_total` (physics convention, Hermitian) with `R_total^T` (kron convention, what appears in the matvec).
**How to avoid:** Define G_left = i*B^T - 0.5*R_total^T and G_right = -i*B^T - 0.5*R_total^T. Since R_total is Hermitian, R_total^T = conj(R_total). Use `transpose(R_total)` explicitly when building G_left/G_right.
**Warning signs:** Round-trip tests failing with errors that scale with the coherent or dissipator strength.

### Pitfall 2: Adjoint Swap Direction
**What goes wrong:** Getting the adjoint G_left_adj/G_right_adj wrong. The forward has G_left on the left and G_right on the right. The adjoint SWAPS them: G_left_adj = G_right, G_right_adj = G_left.
**Why it happens:** This comes from two sign flips that cancel in the anticommutator but swap in the coherent term:
- Coherent: forward = `i*B^T*rho - i*rho*B^T`, adjoint = `-i*B^T*rho + i*rho*B^T`
- Anticommutator: forward = `-0.5*(R^T*rho + rho*R^T)`, adjoint = same (self-adjoint for Hermitian R)
- So forward G_left = `i*B^T - 0.5*R^T`, forward G_right = `-i*B^T - 0.5*R^T`
- Adjoint G_left_adj = `-i*B^T - 0.5*R^T` = G_right, adjoint G_right_adj = `i*B^T - 0.5*R^T` = G_left
**How to avoid:** Verify with the adjoint duality test: `tr(X' * L(Y)) == tr(L*(X)' * Y)`.
**Warning signs:** Adjoint round-trip tests fail but forward tests pass.

### Pitfall 3: Hermitian Half-Grid Negative Frequency Anticommutator
**What goes wrong:** For Hermitian operators on the half-grid, the negative-frequency partner uses `L_neg = L'`, so `L_neg'L_neg = LL'` (not `L'L`). The R_total accumulation must include both `rate_w * L'L` AND `rate_{-w} * LL'`.
**Why it happens:** Forgetting that the Hermitian half-grid optimization means negative frequencies use a different product.
**How to avoid:** The existing `_accumulate_R_total!` already handles this correctly (see krylov_workspace.jl lines 218-223 for EnergyDomain). Since we reuse the existing function, this pitfall is automatically avoided.
**Warning signs:** Errors that appear only for Hermitian jump operators.

### Pitfall 4: BohrDomain Adjoint Sandwich
**What goes wrong:** The BohrDomain adjoint dissipator uses `_accumulate_adjoint_dissipator_2op!` which has a fundamentally different structure from simply swapping A and B_dag. The adjoint sandwich is `conj(B_dag) * rho * conj(A)` while the forward is `B_dag^T * rho * A^T`.
**Why it happens:** The HS adjoint of `rho -> X*rho*Y` is `rho -> X^H*rho*Y^H`, not `rho -> Y*rho*X`.
**How to avoid:** Keep separate `_accumulate_sandwich_2op!` and `_accumulate_adjoint_sandwich_2op!` functions. Do NOT try to unify them.
**Warning signs:** BohrDomain adjoint duality test fails.

### Pitfall 5: KrylovWorkspace Struct Field Ordering
**What goes wrong:** Julia structs have positional constructors. Adding new fields (G_left, G_right) to `KrylovWorkspace` breaks all existing constructor call sites that use positional arguments.
**Why it happens:** The inner constructor in the struct definition and the outer constructors pass fields positionally.
**How to avoid:** Add G_left and G_right after the channel fields (at the end), update both outer constructors (LiouvConfig and ThermalizeConfig), and ensure the positional order matches.
**Warning signs:** `MethodError: no method matching KrylovWorkspace(...)` at construction time.

### Pitfall 6: B_total = nothing Case
**What goes wrong:** When `with_coherent = false` (GNS configs), `B_total` is `nothing`. G_left and G_right must still be computed from R_total alone: `G_left = G_right = -0.5*R_total^T`.
**Why it happens:** Branching on `B_total !== nothing` is already done in the existing apply_lindbladian! code. The precomputation must handle this case.
**How to avoid:** In the constructor, compute `G_left = -0.5 * transpose(R_total)` when `B_total === nothing`, and `G_left = 1im * transpose(B_total) - 0.5 * transpose(R_total)` otherwise.
**Warning signs:** GNS round-trip tests fail or G computation throws a TypeError on `nothing`.

### Pitfall 7: Zero-Allocation Guarantee
**What goes wrong:** The current matvec achieves zero allocations (verified by `@allocated` tests in test_krylov_matvec.jl). The optimized code must maintain this property.
**Why it happens:** Using `transpose()` or `adjoint()` wrapper types with `BLAS.gemm!` can cause boxing allocations. The current code avoids this by using `'T'` and `'C'` character flags directly.
**How to avoid:** Apply G_left and G_right directly via `BLAS.gemm!('N', 'N', ...)` since they are pre-transposed at construction time. Store G_left as the already-transposed matrix to avoid runtime transpose wrappers.
**Warning signs:** `@allocated` tests fail with 176-byte allocations (the signature of Adjoint/Transpose wrapper boxing).

## Code Examples

### Example 1: G_left/G_right Computation at Construction

```julia
# In KrylovWorkspace constructor (LiouvConfig path)
# After computing precomputed_data, B_total, jump_eigenbases, jump_hermitian...

# Compute R_total (reuse existing helper)
R_total = zeros(CT, dim, dim)
_accumulate_R_total!(R_total, jump_eigenbases, jump_hermitian,
                     precomputed_data, config, ham_or_trott)
hermitianize!(R_total)

# Build G_left and G_right (store pre-transposed for zero-alloc BLAS)
# Convention: G_left * rho means gemm!('N','N', ..., G_left, rho, ...)
# where G_left is stored as (i*B^T - 0.5*R_total^T) ready to left-multiply rho
R_total_T = transpose(R_total)  # since R is Hermitian, R^T = conj(R)
if B_total !== nothing
    B_T = transpose(B_total)
    G_left  = 1im .* B_T .- 0.5 .* R_total_T
    G_right = -1im .* B_T .- 0.5 .* R_total_T
else
    G_left  = -0.5 .* R_total_T
    G_right = -0.5 .* R_total_T  # same as G_left when no coherent
end
# Convert to Matrix{CT} for type stability
G_left  = Matrix{CT}(G_left)
G_right = Matrix{CT}(G_right)
```

### Example 2: Optimized apply_lindbladian! (EnergyDomain)

```julia
function apply_lindbladian!(
    ws::KrylovWorkspace{T},
    rho::Matrix{T},
    config::AbstractLiouvConfig{EnergyDomain},
    hamiltonian::HamHam,
) where {T<:Complex}
    (; transition, gamma_norm_factor, energy_labels) = ws.precomputed_data
    bohr_freqs = hamiltonian.bohr_freqs
    inv_4sigma2 = 1.0 / (4 * config.sigma^2)
    CT = one(T)
    ZT = zero(T)

    # Effective Hamiltonian: rho_out = G_left * rho + rho * G_right
    BLAS.gemm!('N', 'N', CT, ws.G_left, rho, ZT, ws.rho_out)     # G_left * rho
    BLAS.gemm!('N', 'N', CT, rho, ws.G_right, CT, ws.rho_out)     # + rho * G_right

    # Sandwich-only loop: sum_i scalar_i * conj(L_i) * rho * L_i^T
    prefactor = (config.w0 / (config.sigma * sqrt(2 * pi))) * gamma_norm_factor
    for (k, eigenbasis) in enumerate(ws.jump_eigenbases)
        is_herm = ws.jump_hermitian[k]
        if is_herm
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)
                _krylov_oft!(ws.jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)
                scalar_w = prefactor * transition(w)
                _accumulate_sandwich!(ws.rho_out, ws.jump_oft, rho, scalar_w, ws)
                if w > 1e-12
                    scalar_neg = prefactor * transition(-w)
                    _accumulate_sandwich_adj_L!(ws.rho_out, ws.jump_oft, rho, scalar_neg, ws)
                end
            end
        else
            for w in energy_labels
                _krylov_oft!(ws.jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)
                scalar_w = prefactor * transition(w)
                _accumulate_sandwich!(ws.rho_out, ws.jump_oft, rho, scalar_w, ws)
            end
        end
    end
    return ws.rho_out
end
```

### Example 3: Optimized apply_adjoint_lindbladian! (EnergyDomain)

```julia
function apply_adjoint_lindbladian!(
    ws::KrylovWorkspace{T},
    rho::Matrix{T},
    config::AbstractLiouvConfig{EnergyDomain},
    hamiltonian::HamHam,
) where {T<:Complex}
    # ... same setup ...
    CT = one(T)
    ZT = zero(T)

    # Adjoint effective Hamiltonian: swap G_left and G_right
    BLAS.gemm!('N', 'N', CT, ws.G_right, rho, ZT, ws.rho_out)    # G_right * rho (was G_left)
    BLAS.gemm!('N', 'N', CT, rho, ws.G_left, CT, ws.rho_out)     # + rho * G_left (was G_right)

    # Adjoint sandwich loop: sum_i scalar_i * L_i^T * rho * conj(L_i)
    # ... uses _accumulate_adjoint_sandwich! instead of _accumulate_sandwich! ...
end
```

## Detailed GEMM Accounting

### Current Cost (per matvec, N dissipator terms)

| Operation | GEMMs per term | Count | Total GEMMs |
|-----------|---------------|-------|-------------|
| Coherent: B^T*rho, rho*B^T | 1 each | 1 | 2 |
| _accumulate_dissipator! | 5 | N | 5N |
| **Total** | | | **2 + 5N** |

### After Optimization

| Operation | GEMMs per term | Count | Total GEMMs |
|-----------|---------------|-------|-------------|
| G_left*rho, rho*G_right | 1 each | 1 | 2 |
| _accumulate_sandwich! | 2 | N | 2N |
| **Total** | | | **2 + 2N** |

### Savings: 3N GEMMs removed per matvec

For a typical system with 12 jumps and ~30 energy labels (Hermitian half-grid ~ 15 effective labels per jump = ~180 terms), this removes ~540 GEMMs per matvec. Over ~100 Krylov iterations, that's ~54,000 fewer BLAS calls.

## Existing _accumulate_R_total! Analysis

The R_total accumulation is already implemented for all 4 domains in `krylov_workspace.jl`:

### EnergyDomain (lines 190-235)
- Physics convention: `R += rate^2 * L'L` using `mul!(LdagL, jump_oft', jump_oft)`
- Handles Hermitian half-grid: adds both `rate_w * L'L` and `rate_{-w} * L*L'`
- Uses same OFT computation as the matvec

### TimeDomain/TrotterDomain (lines 237-280)
- Same structure, but uses NUFFT prefactors instead of Gaussian OFT
- Prefactor formula matches the matvec exactly

### BohrDomain (lines 282-314)
- Different structure: `R += gamma_norm_factor * (A_nu2_dag * alpha_A)`
- This is the anticommutator base `B_dag * A` in the two-operator dissipator
- The same R_total serves for G_left/G_right because the anticommutator in `_accumulate_dissipator_2op!` is `{(B_dag*A)^T, rho}` which matches `{R_total^T, rho}`

**All 3 methods use `mul!` (not `BLAS.gemm!`) for the accumulation because they run once at construction, not per-matvec. This is acceptable.**

## Discretion Recommendations

### G_left / G_right Storage: Store Both as Separate Fields

**Recommendation: Store both G_left and G_right as separate `Matrix{T}` fields in `KrylovWorkspace`.**

Reasons:
1. **Correctness:** The relationship `G_right = -G_left'` (adjoint) or `G_right = -conj(G_left)` requires careful verification of the conjugation/transpose interaction with the kron convention. Storing both avoids this complexity entirely.
2. **Performance:** Two dim x dim matrices cost negligible memory compared to the workspace (which already has 5 scratch matrices). At dim=16 (4 qubits), each is 2 KB. At dim=64 (6 qubits), each is 32 KB.
3. **Zero-allocation:** Pre-stored matrices can be passed directly to `BLAS.gemm!` with `'N','N'` flags. Computing on-the-fly would require temporary allocation or wrapper types.
4. **When B_total is nothing:** G_left = G_right = -0.5*R_total^T, so we could share one pointer. But storing separate copies simplifies the code and costs nothing.

### Scratch Matrix Reuse: Use Existing tmp1, tmp2, LdagL

**Recommendation: The existing 3 scratch matrices (tmp1, tmp2, LdagL) are sufficient for the sandwich-only loop.**

The sandwich computation for the 1-op case needs:
- `tmp2` <- `conj(L_op)` (elementwise, in-place)
- `tmp1` <- `conj(L)*rho` via `gemm!('N','N', tmp2, rho -> tmp1)`
- `LdagL` <- `conj(L)*rho*L^T` via `gemm!('N','T', tmp1, L_op -> LdagL)`
- `axpy!(scalar, LdagL, out)`

This matches the existing tail of `_accumulate_dissipator!` (lines 59-62 of krylov_matvec.jl). Since we're removing the L'L and anticommutator computations that previously used these buffers, the same buffers are available.

For the 2-op BohrDomain sandwich:
- `tmp1` <- `B_dag^T * rho` via `gemm!('T','N', B_dag, rho -> tmp1)`
- `tmp2` <- `B_dag^T * rho * A^T` via `gemm!('N','T', tmp1, A -> tmp2)`
- `axpy!(scalar, tmp2, out)`

Again, 2 scratch matrices suffice.

### BohrDomain R_total: Share the Same Pattern

**Recommendation: BohrDomain uses the same `_accumulate_R_total!` and the same G_left/G_right precomputation pattern.**

The Bohr R_total already accumulates `B_dag * A` which is the correct anticommutator base for the two-operator dissipator. The G_left/G_right definitions are identical:
- `G_left = i*B^T - 0.5*R_total^T`
- `G_right = -i*B^T - 0.5*R_total^T`

The per-matvec sandwich loop changes from `_accumulate_dissipator_2op!` (5 GEMMs) to `_accumulate_sandwich_2op!` (2 GEMMs), same reduction ratio.

## Adjoint Dissipator Anticommutator Verification

A critical verification: the adjoint dissipator has the SAME anticommutator as the forward, so the same R_total works for both.

**Forward:** `D_L(rho) = conj(L)*rho*L^T - 0.5*(L'L)^T*rho - 0.5*rho*(L'L)^T`
**Adjoint:** `D_L*(rho) = L^T*rho*conj(L) - 0.5*(L'L)^T*rho - 0.5*rho*(L'L)^T`

The anticommutator `{(L'L)^T, rho}` is identical. Only the sandwich changes direction. This is confirmed in the code: `_accumulate_adjoint_dissipator!` (lines 119-147) uses the same `BLAS.gemm!('C','N',...)` for L'L and the same anticommutator GEMMs.

For the BohrDomain 2-op case:
**Forward:** uses `(B_dag*A)^T` in anticommutator
**Adjoint:** uses `conj(B_dag*A)` in anticommutator (line 428 of krylov_matvec.jl)

Wait -- this is a **critical difference**! The BohrDomain adjoint uses `conj(B_dag*A)`, NOT `(B_dag*A)^T`. So the adjoint anticommutator for BohrDomain is different from the forward.

Let me verify: `(B_dag*A)^T` vs `conj(B_dag*A)`. For Hermitian product, these would be different (transpose vs conjugate). However, note that `conj(M) = (M^T)^*` and `M^T = conj(M)^* `, so `conj(M) = M^T` only for real M.

Looking at `_accumulate_adjoint_dissipator_2op!` (lines 413-446):
- Line 428: `@. ws.tmp2 = conj(ws.LdagL)` where LdagL = B_dag*A
- Lines 431-432: `conj(B_dag*A) * rho` and `rho * conj(B_dag*A)`

vs `_accumulate_dissipator_2op!` (lines 365-393):
- Lines 380-385: `(B_dag*A)^T * rho` and `rho * (B_dag*A)^T`

These ARE different: `(M)^T` vs `conj(M)`.

**This means for BohrDomain, the adjoint needs a separate G_left_adj/G_right_adj, OR we need to accept that BohrDomain cannot fully precompute the anticommutator for the adjoint using the same R_total.**

Actually, let's reconsider. The relationship: `conj(M) = conj(B_dag * A)` while the forward uses `(B_dag * A)^T`. Since `conj(M) = transpose(M')` and `M^T = transpose(M)`, these are different by a conjugation.

For the forward G_right: contains `-0.5*R_total^T` where R_total = sum gamma_norm_factor * B_dag * A.
For the adjoint: the anticommutator uses `conj(R_total)` instead of `R_total^T`.

Since R_total is **not necessarily Hermitian** for BohrDomain (it's `sum B_dag * A` which is a general matrix), `R_total^T != conj(R_total)` in general.

**Resolution:** We need to define separate adjoint versions for BohrDomain:
- Forward: `G_left = i*B^T - 0.5*R_total^T`, `G_right = -i*B^T - 0.5*R_total^T`
- Adjoint for Energy/Time/Trotter: `G_left_adj = G_right`, `G_right_adj = G_left` (since R_total is Hermitian so R^T = conj(R))
- Adjoint for BohrDomain: `G_left_adj = -i*B^T - 0.5*conj(R_total)`, `G_right_adj = i*B^T - 0.5*conj(R_total)`

Wait, but we need to be more careful. Let me re-derive:

Forward coherent: `i*B^T*rho - i*rho*B^T`
Adjoint coherent: `-i*B^T*rho + i*rho*B^T`

Forward anticommutator (Energy/Time/Trotter): `-0.5*(L'L)^T * rho - 0.5*rho*(L'L)^T`, R_total Hermitian, so R^T = conj(R)
Adjoint anticommutator (Energy/Time/Trotter): same as forward (anticommutator is self-adjoint for Hermitian R)

Forward anticommutator (Bohr): `-0.5*(B_dag*A)^T * rho - 0.5*rho*(B_dag*A)^T`
Adjoint anticommutator (Bohr): `-0.5*conj(B_dag*A) * rho - 0.5*rho*conj(B_dag*A)`

So for Bohr adjoint, G_left_adj = `-i*B^T - 0.5*conj(R_total)` and G_right_adj = `i*B^T - 0.5*conj(R_total)`.

**Options:**
1. Store 4 matrices: G_left, G_right, G_left_adj, G_right_adj (universal, handles all cases)
2. Store 2 matrices + derive adjoint: works for Energy/Time/Trotter (R Hermitian), needs special handling for Bohr
3. For Bohr, store R_total separately and compute adjoint anticommutator from conj(R_total)

**Recommendation:** Store 4 matrices (G_left, G_right, G_left_adj, G_right_adj). For Energy/Time/Trotter domains where R is Hermitian, G_left_adj = G_right and G_right_adj = G_left automatically. For BohrDomain, they differ. The memory cost is 4 dim x dim matrices = negligible. This is the simplest correct approach.

Actually, looking more carefully: is R_total for BohrDomain really non-Hermitian? Let me check:

`R_total += gamma_norm_factor * A_nu2_dag * alpha_A`

where `A_nu2_dag[j,i] = conj(eigenbasis[i,j])` and `alpha_A = alpha(bohr_freqs, nu_2) * eigenbasis`. The product `A_nu2_dag * alpha_A` is the sum over Bohr buckets of `B_dag * A` terms. This is generally not Hermitian.

However, the existing `_accumulate_R_total!` for the ThermalizeConfig constructor follows this with `hermitianize!(R_total)` (line 380 of krylov_workspace.jl). But that hermitianization is for the Chen channel (which requires a Hermitian R_total for the PSD guard). For the Lindbladian path, the R_total as accumulated is the correct non-Hermitian matrix.

**Key insight:** The forward Lindbladian uses `(B_dag*A)^T` and the adjoint uses `conj(B_dag*A)`. For the forward, R_total^T is correct. For the adjoint, conj(R_total) is correct. These are the same if and only if R_total is Hermitian (since `conj(H) = H^T` for Hermitian H).

For BohrDomain, R_total is NOT Hermitian in general. So we need separate adjoint matrices.

**Final recommendation for BohrDomain:** Store `G_left`, `G_right`, `G_left_adj`, `G_right_adj` as workspace fields. Compute all four at construction. For non-Bohr domains where R is Hermitian, G_left_adj = G_right and G_right_adj = G_left (verified, not assumed). For Bohr, compute independently.

## Scope of Changes

### Files to Modify
1. **`src/krylov_workspace.jl`**: Add G_left, G_right, G_left_adj, G_right_adj fields to struct; compute them in both constructors (LiouvConfig and ThermalizeConfig)
2. **`src/krylov_matvec.jl`**: Replace `_accumulate_dissipator!` calls with G_left/G_right + `_accumulate_sandwich!`; same for adjoint variants and all domains
3. **`src/krylov_eigsolve.jl`**: Delete the legacy 5-arg `apply_delta_channel!` method (lines 208-227)
4. **`test/test_krylov_eigsolve.jl`**: Delete "apply_delta_channel! legacy Euler" testset (lines 53-68)
5. **`src/QuantumFurnace.jl`**: Possibly update export if apply_delta_channel! still has 4-arg form

### Files NOT Modified
- `test/test_krylov_matvec.jl`: All existing tests should pass unchanged (verifying optimization correctness)
- `test/test_krylov_crossvalidation.jl`: All cross-validation tests should pass unchanged
- `src/coherent.jl`: No changes needed
- `src/furnace_utensils.jl`: No changes needed

### New Functions
- `_accumulate_sandwich!`: Forward 1-op sandwich only (2 GEMMs)
- `_accumulate_sandwich_adj_L!`: Forward 1-op negative-frequency sandwich (2 GEMMs)
- `_accumulate_adjoint_sandwich!`: Adjoint 1-op sandwich (2 GEMMs)
- `_accumulate_adjoint_sandwich_adj_L!`: Adjoint 1-op negative-frequency sandwich (2 GEMMs)
- `_accumulate_sandwich_2op!`: Forward BohrDomain sandwich (2 GEMMs)
- `_accumulate_adjoint_sandwich_2op!`: Adjoint BohrDomain sandwich (2 GEMMs)

### Deleted Functions/Methods
- `apply_delta_channel!(ws, rho, delta, config_liouv, hamiltonian)` -- the 5-arg Euler form

### Deleted Tests
- "apply_delta_channel! legacy Euler" testset in test_krylov_eigsolve.jl

## Open Questions

1. **Should the old `_accumulate_dissipator!` and friends be deleted entirely or kept?**
   - What we know: After the optimization, `_accumulate_dissipator!`, `_accumulate_dissipator_adj_L!`, `_accumulate_adjoint_dissipator!`, `_accumulate_adjoint_dissipator_adj_L!`, `_accumulate_dissipator_2op!`, `_accumulate_adjoint_dissipator_2op!` are no longer called by the matvec.
   - What's unclear: Are they used anywhere else (e.g., by the channel `_accumulate_jump_sandwich!` path)? The channel path uses its own sandwich logic in `_accumulate_jump_sandwich!` which does NOT call `_accumulate_dissipator!`. So these functions become dead code.
   - Recommendation: Delete them in this phase to keep the codebase clean. They were only ever called by the matvec path.

2. **BohrDomain adjoint R_total non-Hermiticity verification**
   - What we know: The mathematical analysis shows `conj(R_total) != R_total^T` for non-Hermitian R_total.
   - What's unclear: In practice, is the BohrDomain R_total close enough to Hermitian that the difference is negligible?
   - Recommendation: Do NOT assume Hermiticity. Store 4 matrices and compute all correctly. Validate with existing adjoint duality tests.

## Sources

### Primary (HIGH confidence)
- `src/krylov_matvec.jl` -- Current dissipator implementations, GEMM counts, kron convention
- `src/krylov_workspace.jl` -- KrylovWorkspace struct, `_accumulate_R_total!` for all domains, ThermalizeConfig constructor
- `src/krylov_eigsolve.jl` -- Legacy Euler `apply_delta_channel!`, faithful Chen channel, eigsolve API
- `test/test_krylov_matvec.jl` -- Round-trip correctness tests and zero-allocation checks
- `test/test_krylov_eigsolve.jl` -- Legacy Euler test to delete, Chen channel test to keep

### Secondary (MEDIUM confidence)
- `.planning/phases/32-some-speedup-for-the-krylov-simulator/32-CONTEXT.md` -- User decisions and discretion areas
- `.planning/ROADMAP.md` -- Phase 32 requirements and success criteria

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- No new dependencies, pure refactoring of existing BLAS calls
- Architecture: HIGH -- Direct code analysis of all affected functions, GEMM counting verified line-by-line
- Pitfalls: HIGH -- BohrDomain adjoint non-Hermiticity discovered through code analysis (Pitfall 4 + adjoint verification section), kron convention transpose pitfall identified from existing comments
- BohrDomain adjoint: HIGH -- Code analysis clearly shows `conj(B_dag*A)` vs `(B_dag*A)^T` in the two adjoint variants

**Research date:** 2026-02-25
**Valid until:** 2026-03-25 (stable internal optimization, no external dependencies)
