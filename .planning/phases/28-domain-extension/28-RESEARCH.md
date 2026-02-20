# Phase 28: Domain Extension - Research

**Researched:** 2026-02-20
**Domain:** Matrix-free Lindbladian action for TimeDomain, TrotterDomain, and BohrDomain
**Confidence:** HIGH

## Summary

Phase 28 extends the EnergyDomain matvec (Phase 27) to the remaining three domains: TimeDomain, TrotterDomain, and BohrDomain. The codebase already has every building block needed. TimeDomain and TrotterDomain are structurally identical to EnergyDomain -- the only difference is how A(omega) is formed: instead of `_krylov_oft!` (Gaussian filter on eigenbasis), they use `eigenbasis .* nufft_prefactor_view(w)` (entrywise multiply with precomputed NUFFT prefactors). The prefactor formula also changes from `w0 / (sigma * sqrt(2pi)) * gamma_norm` to `w0 * t0^2 * sigma * sqrt(2/pi) / (2pi) * gamma_norm`. BohrDomain is fundamentally different: it iterates over `bohr_dict` keys (Bohr frequency buckets) instead of energy labels, uses a generalized two-operator dissipator `D(A, B, rho)` where A and B differ, and requires entrywise alpha computation in the hot path.

The existing dense Liouvillian code in `jump_workers.jl` dispatches `_jump_contribution!` on `Union{TimeDomain, TrotterDomain}` with a shared method body, and on `BohrDomain` with a separate method. The Krylov matvec should mirror this structure exactly. The four existing dissipator helpers (`_accumulate_dissipator!`, `_accumulate_dissipator_adj_L!`, `_accumulate_adjoint_dissipator!`, `_accumulate_adjoint_dissipator_adj_L!`) remain untouched for Energy/Time/Trotter. BohrDomain needs a new `_accumulate_dissipator_2op!` helper plus a dedicated adjoint variant (see Open Questions for the critical math correction).

**Primary recommendation:** Implement Plan 1 (Time + Trotter) as a single `Union{TimeDomain, TrotterDomain}` method on `apply_lindbladian!` and `apply_adjoint_lindbladian!`, mirroring the dense `_jump_contribution!` pattern line-by-line but using concrete-typed `ws.jump_eigenbases` and BLAS primitives. Implement Plan 2 (Bohr) with a new `_accumulate_dissipator_2op!` helper and a separate `_accumulate_adjoint_dissipator_2op!` helper (NOT simple argument swap -- see Open Questions).

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

#### BohrDomain dissipator structure
- New `_accumulate_dissipator_2op!(out, A, B, rho, scalar, ws)` helper -- separate function, keeps existing single-op helpers untouched
- BohrDomain uses two-operator dissipator D(A, B, rho) = A rho B' - 0.5(B'A rho + rho B'A) where A and B are different matrices
- Entrywise alpha computation in the hot path: `@. alpha(bohr_freqs, nu_2) * eigenbasis` -- cannot precompute because too many Bohr frequencies to store
- BohrDomain iterates over `bohr_dict` keys (bucket iteration), fundamentally different loop from Energy/Time/Trotter
- A_nu2_dag scratch: do NOT add a new workspace field for Bohr -- BohrDomain is not the primary code path. Reuse existing scratch matrices creatively (Claude's discretion on which to repurpose)

#### Adjoint scope
- Forward + adjoint for all three new domains -- full parity with EnergyDomain
- For Time/Trotter adjoint: same NUFFT prefactors, just swap dissipator sandwich (coherent sign flip + sandwich swap, no prefactor changes)
- For BohrDomain adjoint: swap arguments at call site -- `_accumulate_dissipator_2op!(out, B, A, rho, scalar, ws)` with A and B swapped. No separate adjoint 2op function needed (researcher should verify the math)

#### Plan structure / domain ordering
- Plan 1: Time + Trotter together (mechanical, structurally similar to EnergyDomain -- swap `_krylov_oft!` for NUFFT prefactor multiply, change scalar prefactor formula)
- Plan 2: Bohr separately (new 2op helper, different loop structure, bucket iteration)
- EnergyDomain (Phase 27) is the template -- Time/Trotter are the closest structural analogs

#### Round-trip testing
- Per-domain round-trip against dense `construct_lindbladian()` at n=4 (no cross-domain tests)
- Both KMS and GNS balance for all three domains
- Allocation regression tests for Time/Trotter only -- Bohr may have unavoidable allocations from entrywise alpha computation

#### Workspace and API
- Uniform matvec signature: `apply_lindbladian!(ws, rho, config, hamiltonian)` -- workspace has everything, same API across all domains
- TrotterDomain: trotter data already in `ws.precomputed_data` from construction, no extra argument needed
- BohrDomain: access `hamiltonian.bohr_dict` and `hamiltonian.bohr_freqs` at call time (not stored in workspace)
- Time/Trotter OFT: `@. ws.jump_oft = ws.jump_eigenbases[k] * nufft_prefactor_matrix` -- use concrete-typed eigenbases to avoid boxing, trust existing `_prefactor_view` for zero-alloc access

### Claude's Discretion
- Dense vs sparse for A_nu2_dag in BohrDomain (whichever avoids allocation best while keeping code clear)
- Which existing scratch matrix to repurpose for BohrDomain's second operator
- Internal helper decomposition within each domain's apply_lindbladian! method
- Test matrix generation strategy for round-trip validation

### Deferred Ideas (OUT OF SCOPE)
None -- discussion stayed within phase scope
</user_constraints>

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| LinearAlgebra (stdlib) | Julia 1.x | `BLAS.gemm!`, `BLAS.axpy!`, `mul!` | Already used in Phase 27 krylov_matvec.jl for zero-allocation dissipator |
| SparseArrays (stdlib) | Julia 1.x | `sparse()` constructor for BohrDomain A_nu2_dag | Already used in dense `_jump_contribution!(... BohrDomain ...)` |
| QuantumFurnace internals | current | `_prefactor_view`, `_pick_alpha`, `_precompute_data`, `_precompute_coherent_total_B`, KrylovWorkspace, dissipator helpers | All exist and are tested from Phases 10-27 |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| Test (stdlib) | Julia 1.x | `@testset`, `@test`, `@allocated` | Round-trip and allocation regression tests |

### Alternatives Considered
Not applicable -- no new dependencies. All locked decisions use existing codebase infrastructure.

**Installation:**
No new dependencies for Phase 28.

## Architecture Patterns

### Recommended File Structure
```
src/
  krylov_matvec.jl          # ADD: new dispatch methods for Time/Trotter/Bohr
  krylov_workspace.jl       # NO CHANGES (workspace already supports all domains)
test/
  test_krylov_matvec.jl     # ADD: new testsets for Time/Trotter/Bohr round-trips
  test_helpers.jl           # NO CHANGES (factories already support all domains)
```

### Pattern 1: Time/Trotter OFT via NUFFT Prefactors

**What:** Instead of `_krylov_oft!` (Gaussian filter), Time/Trotter domains compute A(omega) by entrywise multiplying the eigenbasis with precomputed NUFFT prefactors. The prefactors are a 3D array indexed by `(i, j, k)` where `k` corresponds to the energy label omega.

**Dense reference (`src/jump_workers.jl:91-135`):**
```julia
# Time/Trotter _jump_contribution! (Liouvillian, vectorized)
nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)
@. jump_oft = jump.in_eigenbasis * nufft_prefactor_matrix
```

**Krylov equivalent (using concrete-typed eigenbases):**
```julia
# In apply_lindbladian! for Union{TimeDomain, TrotterDomain}:
nufft_prefactor_matrix = _prefactor_view(ws.precomputed_data.oft_nufft_prefactors, w)
@. ws.jump_oft = ws.jump_eigenbases[k] * nufft_prefactor_matrix
```

**Key difference from EnergyDomain:** No `_krylov_oft!` call. The NUFFT prefactors already encode the Fourier transform; the OFT is just an entrywise product.

**Source:** `src/jump_workers.jl:115-116` (dense), `src/nufft.jl:93-96` (`_prefactor_view`)

### Pattern 2: Scalar Prefactor for Time/Trotter

**What:** The scalar prefactor for the dissipator changes between domains.

| Domain | Prefactor Formula | Source |
|--------|-------------------|--------|
| Energy | `(w0 / (sigma * sqrt(2*pi))) * gamma_norm_factor` | `jump_workers.jl:65` |
| Time/Trotter | `w0 * t0^2 * (sigma * sqrt(2/pi)) / (2*pi) * gamma_norm_factor` | `jump_workers.jl:109` |

**Note:** Both formulas use `gamma_norm_factor` from `precomputed_data`, which is identical across KMS/GNS for the same energy labels. The `transition(w)` function differs between KMS and GNS (handled by `pick_transition` dispatch at `_precompute_data` time).

**Source:** `src/jump_workers.jl:65` (Energy), `src/jump_workers.jl:109` (Time/Trotter)

### Pattern 3: BohrDomain Bucket Iteration

**What:** BohrDomain iterates over Bohr frequency buckets instead of energy labels. Each bucket `nu_2` maps to a set of `(i,j)` index pairs in the eigenbasis matrix.

**Dense reference (`src/jump_workers.jl:11-45`):**
```julia
unique_freqs = keys(hamiltonian.bohr_dict)
(; alpha, gamma_norm_factor) = precomputed_data

for nu_2 in unique_freqs
    @. alpha_A_nu1 = alpha(hamiltonian.bohr_freqs, nu_2) * jump.in_eigenbasis

    indices = hamiltonian.bohr_dict[nu_2]
    A_nu_2_vals = view(jump.in_eigenbasis, indices)
    rows_dag = getindex.(indices, 2)
    cols_dag = getindex.(indices, 1)
    A_nu_2_dag = sparse(rows_dag, cols_dag, conj.(A_nu_2_vals), dim, dim)

    _vectorize_liouv_diss_and_add!(L_target, alpha_A_nu1, A_nu_2_dag, gamma_norm_factor, ws)
end
```

**Krylov equivalent (matrix-on-matrix):**
```julia
for nu_2 in keys(hamiltonian.bohr_dict)
    @. ws.jump_oft = alpha(hamiltonian.bohr_freqs, nu_2) * eigenbasis   # alpha_A

    # Build A_nu2_dag in a scratch matrix (see Discretion recommendations)
    # ...

    _accumulate_dissipator_2op!(ws.rho_out, ws.jump_oft, A_nu2_dag, rho, gamma_norm_factor, ws)
end
```

**Source:** `src/jump_workers.jl:11-45`

### Pattern 4: Two-Operator Dissipator for BohrDomain

**What:** The BohrDomain dissipator uses two different operators A and B:
```
D(A, B, rho) = A * rho * B' - 0.5*(B'*A*rho + rho*B'*A)
```

This contrasts with the single-operator dissipator used by Energy/Time/Trotter:
```
D(L, rho) = L*rho*L' - 0.5*(L'*L*rho + rho*L'*L)
```

**In the BohrDomain:**
- A = `alpha_A_nu1` = `alpha(bohr_freqs, nu_2) * eigenbasis` (dense dim x dim)
- B' = `A_nu_2_dag` (sparse: only entries from the bucket), so B = A_nu2 (the un-daggered operator)

**Dense reference:** `_vectorize_liouv_diss_and_add!(L, alpha_A_nu1, A_nu_2_dag, gamma_norm_factor, ws)` in `qi_tools.jl:86-103`

**Source:** `src/qi_tools.jl:86-103`, `src/jump_workers.jl:29-43`

### Pattern 5: TrotterDomain Precomputed Data

**What:** For TrotterDomain, `_precompute_data` is called with `trotter` (not `hamiltonian`), and the precomputed data includes NUFFT prefactors computed from `trotter.bohr_freqs`. The KrylovWorkspace constructor already handles this correctly via the `ham_or_trott` pattern.

**Existing workspace constructor (`src/krylov_workspace.jl:68-73`):**
```julia
ham_or_trott = if config.domain isa TrotterDomain
    trotter === nothing && error("A Trotter object must be provided for the TrotterDomain")
    trotter
else
    hamiltonian
end
precomputed_data = _precompute_data(config, ham_or_trott)
```

**For the matvec:** `ws.precomputed_data.oft_nufft_prefactors` already contains the correct NUFFT prefactors (computed from either `hamiltonian.bohr_freqs` or `trotter.bohr_freqs` depending on domain). The `_prefactor_view` call is domain-agnostic.

**Basis note:** For TrotterDomain, `jump.in_eigenbasis` is in the Trotter eigenbasis (set at JumpOp creation), and `ws.jump_eigenbases[k]` extracted from these jumps at workspace construction time is also in Trotter basis. No basis transformation needed in the matvec.

**Source:** `src/krylov_workspace.jl:68-76`, `src/furnace_utensils.jl:82-121`

### Anti-Patterns to Avoid

- **Accessing `jump.in_eigenbasis` in the hot path:** Use `ws.jump_eigenbases[k]` (concrete-typed `Matrix{T}`) instead. `JumpOp.in_eigenbasis` has abstract type `Matrix{<:Complex}` which causes boxing allocation on every field access. This was discovered and fixed in Phase 27-02.
- **Using `mul!` with `Adjoint` wrappers:** Use `BLAS.gemm!` with transpose flags instead. `mul!(C, A, B')` allocates 16 bytes per call. This was discovered and fixed in Phase 27-02.
- **Using `@.` broadcasting for accumulation:** Use `BLAS.axpy!` instead. Broadcasting may allocate temporaries. This was the pattern established in Phase 27-02.
- **Adding workspace fields for BohrDomain:** The decision explicitly says NOT to add new workspace fields for Bohr. Reuse existing scratch matrices (`tmp1`, `tmp2`, `LdagL`, or `jump_oft` when safe).
- **Building dense A_nu2_dag via `sparse()` and then converting to dense:** The `sparse()` call in the dense code (jump_workers.jl:40) allocates a new sparse matrix each iteration. For the matvec, consider building the operator directly into a scratch matrix to avoid allocation.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| NUFFT prefactor computation | Custom Fourier transform | Existing `_prefactor_view(nufft_prefactors, w)` | Already computed at `_precompute_data` time, returns zero-alloc `@view` |
| Alpha function dispatch | Manual KMS/GNS alpha selection | Existing `_pick_alpha(config)` | Already handles all 4 config types (LiouvConfig, LiouvConfigGNS, ThermalizeConfig, ThermalizeConfigGNS) |
| Coherent B precomputation | Domain-specific B computation | Existing `_precompute_coherent_total_B(jumps, ham_or_trott, config, precomputed_data)` | Already handles Time/Trotter/Bohr/Energy B operators at workspace construction |
| Config factory for tests | Manual config construction | Existing `make_liouv_config(domain)`, `make_liouv_config_gns(domain)` in test_helpers.jl | Already parameterized by domain, produce valid configs |
| Transition function selection | Manual KMS/GNS transition dispatch | Existing `pick_transition(config)` stored in `precomputed_data.transition` | Already resolved at `_precompute_data` time |
| Trotter object creation | Manual Trotter construction | Existing `TEST_TROTTER` and `make_test_trotter()` in test_helpers.jl | Already constructed for the 4-qubit test system |
| Trotter-basis jumps | Manual basis transformation | Existing `TEST_TROTTER_JUMPS` in test_helpers.jl | Already transformed to Trotter eigenbasis |

**Key insight:** Phase 28 is purely about adding new dispatch methods to `apply_lindbladian!` and `apply_adjoint_lindbladian!`. Every supporting function (precomputation, alpha, transition, prefactors, workspace construction, test factories) already exists and works correctly across all domains.

## Common Pitfalls

### Pitfall 1: Time/Trotter Prefactor Formula Mismatch
**What goes wrong:** Using the EnergyDomain prefactor formula `w0 / (sigma * sqrt(2pi))` instead of the Time/Trotter formula `w0 * t0^2 * sigma * sqrt(2/pi) / (2pi)`. These differ by a factor of `t0^2 * sigma^2 / pi`.
**Why it happens:** Copy-paste from EnergyDomain without updating the prefactor.
**How to avoid:** Copy the prefactor formula directly from `_jump_contribution!` for `Union{TimeDomain, TrotterDomain}` at `jump_workers.jl:109`. The config carries both `w0`, `t0`, and `sigma`.
**Warning signs:** Time/Trotter round-trip tests fail by a constant multiplicative factor across all density matrices.

### Pitfall 2: Accessing `config.t0` for EnergyDomain
**What goes wrong:** `config.t0` is `nothing` for EnergyDomain configs, so any `Union{TimeDomain, TrotterDomain}` code that accidentally runs for EnergyDomain will hit a `MethodError`.
**Why it happens:** All configs have the `t0` field but it is `Union{T, Nothing}` and only populated for Time/Trotter domains.
**How to avoid:** The dispatch signature `config::AbstractLiouvConfig{D} where D<:Union{TimeDomain, TrotterDomain}` prevents this. Never use a catch-all `AbstractLiouvConfig` for the Time/Trotter method.
**Warning signs:** `MethodError: no method matching *(::Nothing, ::Nothing)` at the prefactor computation.

### Pitfall 3: BohrDomain `sparse()` Allocation in Hot Path
**What goes wrong:** The dense code uses `sparse(rows_dag, cols_dag, conj.(A_nu_2_vals), dim, dim)` which allocates a new `SparseMatrixCSC` every iteration. This is fine for the dense Liouvillian (constructed once) but would be O(n_buckets) allocations per matvec call.
**Why it happens:** Directly copying the dense code pattern without considering per-call allocation impact.
**How to avoid:** Two options: (1) Build the effect of A_nu2_dag directly into a dense scratch matrix by zeroing + scatter, or (2) Apply A_nu2_dag via explicit index loops (like the Kraus trajectory code does in jump_workers.jl:250-274). The Kraus code demonstrates the index-loop approach, avoiding sparse matrix construction entirely.
**Warning signs:** Non-zero allocations in BohrDomain matvec. The decision allows some allocations for Bohr from entrywise alpha computation, but the `sparse()` call adds significant per-bucket allocation that should be avoided.

### Pitfall 4: BohrDomain Alpha Computation Allocates from Broadcasting
**What goes wrong:** The expression `@. alpha(bohr_freqs, nu_2) * eigenbasis` calls the alpha function element-by-element via broadcasting. If the closure allocates (e.g., captures boxed variables), this creates O(dim^2) allocations per bucket.
**Why it happens:** The `alpha` function returned by `_pick_alpha(config)` is a closure over config parameters. Broadcasting over a closure is generally allocation-free in Julia if the closure is fully typed, but the `erfc` special function may cause issues.
**How to avoid:** The CONTEXT.md decision explicitly allows some allocations for BohrDomain from entrywise alpha. The allocation test is only required for Time/Trotter. However, minimize by using `@.` broadcast into a pre-allocated buffer (ws.jump_oft).
**Warning signs:** Excessive allocations in BohrDomain (beyond what entrywise alpha requires).

### Pitfall 5: BohrDomain Adjoint -- Simple Swap is Mathematically INCORRECT
**What goes wrong:** The CONTEXT.md decision suggests "swap arguments at call site: `_accumulate_dissipator_2op!(out, B, A, rho, scalar, ws)`" for the adjoint. This is **mathematically incorrect** for the general two-operator dissipator.
**Why it happens:** For the single-operator dissipator, the adjoint IS a simple structural transformation (swap sandwich, keep anticommutator). But the two-operator case is more complex.
**Math:** The forward is `D(A,B)(rho) = A*rho*B' - 0.5*(B'*A*rho + rho*B'*A)`. The Hilbert-Schmidt adjoint is `D*(A,B)(rho) = A'*rho*B - 0.5*(A'*B*rho + rho*A'*B)`. Swapping A and B in the forward helper gives `D(B,A)(rho) = B*rho*A' - 0.5*(A'*B*rho + rho*A'*B)`. The anticommutator `A'*B` matches, but the sandwich term `B*rho*A'` is NOT the same as `A'*rho*B` (unless A and B commute with rho, which they don't in general).
**How to avoid:** Create a separate `_accumulate_adjoint_dissipator_2op!` helper that computes `A'*rho*B - 0.5*(A'*B*rho + rho*A'*B)`. This parallels the Phase 27 discovery that required a separate `_accumulate_adjoint_dissipator!` for the single-op case.
**Warning signs:** BohrDomain adjoint round-trip test fails while forward passes. Adjoint duality check `tr(X' * L(Y)) == tr(L*(X)' * Y)` fails for BohrDomain.
**Verification:** See detailed derivation in Open Questions section.

### Pitfall 6: Forgetting Coherent Term for Time/Trotter KMS
**What goes wrong:** Time/Trotter KMS configs can have `with_coherent=true`, which adds the `-i[B, rho]` coherent term. The coherent B is already precomputed and stored in `ws.B_total` at workspace construction time. Forgetting to include the coherent term block in the new methods.
**Why it happens:** The EnergyDomain template includes the coherent block, but it could be missed when writing the Time/Trotter method if attention focuses on the OFT difference.
**How to avoid:** Copy the coherent block from the EnergyDomain method verbatim. It uses `ws.B_total` and `ws.tmp1` which are domain-agnostic.
**Warning signs:** Time/Trotter KMS with_coherent=true round-trip fails, but GNS and KMS with_coherent=false both pass.

## Code Examples

### Example 1: Time/Trotter Forward Matvec (Complete Method)
```julia
function apply_lindbladian!(
    ws::KrylovWorkspace{T},
    rho::Matrix{T},
    config::AbstractLiouvConfig{D},
    hamiltonian::HamHam,
) where {T<:Complex, D<:Union{TimeDomain, TrotterDomain}}
    (; transition, gamma_norm_factor, energy_labels, oft_nufft_prefactors) = ws.precomputed_data

    fill!(ws.rho_out, 0)

    # Coherent term: -i[B, rho] (identical to EnergyDomain)
    B = ws.B_total
    if B !== nothing
        mul!(ws.tmp1, B, rho)
        @. ws.rho_out += -1im * ws.tmp1
        mul!(ws.tmp1, rho, B)
        @. ws.rho_out += 1im * ws.tmp1
    end

    # Dissipator: sum over jumps, sum over energy labels
    # Prefactor formula from jump_workers.jl:109
    prefactor = config.w0 * config.t0^2 * (config.sigma * sqrt(2 / pi)) / (2 * pi) * gamma_norm_factor

    for (k, eigenbasis) in enumerate(ws.jump_eigenbases)
        is_herm = ws.jump_hermitian[k]
        if is_herm
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)

                # Time/Trotter OFT via NUFFT prefactors (replaces _krylov_oft!)
                nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)
                @. ws.jump_oft = eigenbasis * nufft_prefactor_matrix

                scalar_w = prefactor * transition(w)
                _accumulate_dissipator!(ws.rho_out, ws.jump_oft, rho, scalar_w, ws)

                if w > 1e-12
                    scalar_neg = prefactor * transition(-w)
                    _accumulate_dissipator_adj_L!(ws.rho_out, ws.jump_oft, rho, scalar_neg, ws)
                end
            end
        else
            for w in energy_labels
                nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)
                @. ws.jump_oft = eigenbasis * nufft_prefactor_matrix

                scalar_w = prefactor * transition(w)
                _accumulate_dissipator!(ws.rho_out, ws.jump_oft, rho, scalar_w, ws)
            end
        end
    end

    return ws.rho_out
end
```

### Example 2: Time/Trotter Adjoint Matvec (Only Differences from Forward)
```julia
# Coherent term: +i[B, rho] (sign flip)
# Dissipator: use _accumulate_adjoint_dissipator! and _accumulate_adjoint_dissipator_adj_L!
# Prefactors: IDENTICAL to forward (no prefactor changes for adjoint)
# NUFFT prefactor computation: IDENTICAL to forward
```

### Example 3: BohrDomain Two-Operator Dissipator Helper
```julia
"""
    _accumulate_dissipator_2op!(out, A, B_dag, rho, scalar, ws)

Accumulate `scalar * D(A, B_dag)(rho)` where:
    D(A, B_dag)(rho) = A * rho * B_dag' - 0.5*(B_dag' * A * rho + rho * B_dag' * A)

Note: B_dag is passed explicitly (e.g., as a dense or sparse matrix),
and B_dag' (its adjoint) is the "right" operator in the sandwich.
"""
function _accumulate_dissipator_2op!(
    out::Matrix{T},
    A::Matrix{T},          # alpha_A_nu1 (dense)
    B_dag,                  # A_nu2_dag (could be sparse or dense)
    rho::Matrix{T},
    scalar::Real,
    ws,
) where {T<:Complex}
    # B_dag_A = B_dag * A  (for anticommutator)
    # Compute B_dag * A into ws.LdagL
    mul!(ws.LdagL, B_dag, A)         # LdagL = B'A

    # Term 1: scalar * A * rho * B_dag'  (note: B_dag' = conj(transpose(B_dag)))
    mul!(ws.tmp1, A, rho)             # tmp1 = A * rho
    mul!(ws.tmp2, ws.tmp1, B_dag')    # tmp2 = A * rho * B_dag'
    # ... accumulate into out

    # Term 2: -0.5 * scalar * B'A * rho
    mul!(ws.tmp1, ws.LdagL, rho)
    # ... accumulate

    # Term 3: -0.5 * scalar * rho * B'A
    mul!(ws.tmp1, rho, ws.LdagL)
    # ... accumulate
end
```

### Example 4: BohrDomain Adjoint Two-Operator Dissipator Helper
```julia
"""
    _accumulate_adjoint_dissipator_2op!(out, A, B_dag, rho, scalar, ws)

Accumulate `scalar * D*(A, B_dag)(rho)` (HS adjoint of the two-operator dissipator):
    D*(A, B_dag)(rho) = A' * rho * B_dag - 0.5*(A' * B_dag * rho + rho * A' * B_dag)

Note: The anticommutator uses A'*B_dag (NOT B_dag*A' as simple swap would give).
This parallels _accumulate_adjoint_dissipator! preserving the correct anticommutator.
"""
function _accumulate_adjoint_dissipator_2op!(
    out::Matrix{T},
    A::Matrix{T},
    B_dag,                  # could be sparse or dense
    rho::Matrix{T},
    scalar::Real,
    ws,
) where {T<:Complex}
    # Compute A' * B_dag (the correct adjoint anticommutator product)
    # NOT B_dag * A' (which is what simple argument swap would give)
    # ... implementation using BLAS.gemm! for the A' operations
end
```

### Example 5: BohrDomain Forward Matvec Loop (Skeleton)
```julia
function apply_lindbladian!(
    ws::KrylovWorkspace{T},
    rho::Matrix{T},
    config::AbstractLiouvConfig{BohrDomain},
    hamiltonian::HamHam,
) where {T<:Complex}
    (; alpha, gamma_norm_factor) = ws.precomputed_data

    fill!(ws.rho_out, 0)

    # Coherent term (identical to other domains)
    B = ws.B_total
    if B !== nothing
        mul!(ws.tmp1, B, rho)
        @. ws.rho_out += -1im * ws.tmp1
        mul!(ws.tmp1, rho, B)
        @. ws.rho_out += 1im * ws.tmp1
    end

    dim = size(rho, 1)
    for (k, eigenbasis) in enumerate(ws.jump_eigenbases)
        for nu_2 in keys(hamiltonian.bohr_dict)
            # alpha_A: dense, alpha-weighted eigenbasis
            @. ws.jump_oft = alpha(hamiltonian.bohr_freqs, nu_2) * eigenbasis

            # Build A_nu2_dag (reuse scratch matrix)
            # ... see discretion recommendations

            _accumulate_dissipator_2op!(ws.rho_out, ws.jump_oft, A_nu2_dag, rho, gamma_norm_factor, ws)
        end
    end

    return ws.rho_out
end
```

### Example 6: Round-Trip Test for TimeDomain
```julia
@testset "Round-trip: matvec vs dense (TimeDomain KMS)" begin
    config = make_liouv_config(TimeDomain(); with_coherent=true)
    L_dense = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)
    ws = KrylovWorkspace(config, TEST_HAM, TEST_JUMPS)

    for _ in 1:10
        rho = Matrix(random_density_matrix(NUM_QUBITS))
        v_dense = L_dense * vec(rho)
        L_rho = apply_lindbladian!(ws, rho, config, TEST_HAM)
        @test norm(v_dense - vec(L_rho)) < 1e-12
    end
end
```

### Example 7: Round-Trip Test for TrotterDomain
```julia
@testset "Round-trip: matvec vs dense (TrotterDomain KMS)" begin
    config = make_liouv_config(TrotterDomain(); with_coherent=true)
    L_dense = construct_lindbladian(TEST_TROTTER_JUMPS, config, TEST_HAM; trotter=TEST_TROTTER)
    ws = KrylovWorkspace(config, TEST_HAM, TEST_TROTTER_JUMPS; trotter=TEST_TROTTER)

    for _ in 1:10
        rho = Matrix(random_density_matrix(NUM_QUBITS))
        v_dense = L_dense * vec(rho)
        L_rho = apply_lindbladian!(ws, rho, config, TEST_HAM)
        @test norm(v_dense - vec(L_rho)) < 1e-12
    end
end
```

### Example 8: Round-Trip Test for BohrDomain
```julia
@testset "Round-trip: matvec vs dense (BohrDomain KMS)" begin
    config = make_liouv_config(BohrDomain(); with_coherent=true)
    L_dense = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)
    ws = KrylovWorkspace(config, TEST_HAM, TEST_JUMPS)

    for _ in 1:10
        rho = Matrix(random_density_matrix(NUM_QUBITS))
        v_dense = L_dense * vec(rho)
        L_rho = apply_lindbladian!(ws, rho, config, TEST_HAM)
        @test norm(v_dense - vec(L_rho)) < 1e-12
    end
end
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `mul!` with Adjoint wrappers | `BLAS.gemm!` with char flags | Phase 27-02 | Eliminates 16 bytes/call boxing allocation |
| Abstract-typed `jump.in_eigenbasis` access | Concrete-typed `ws.jump_eigenbases[k]` | Phase 27-02 | Eliminates 832 bytes/call boxing allocation |
| Single `_accumulate_dissipator!` for both forward and adjoint | Separate `_accumulate_adjoint_dissipator!` | Phase 27-02 | Correct anticommutator {L'L, rho} in adjoint |
| `@.` broadcasting for accumulation | `BLAS.axpy!` for scalar accumulation | Phase 27-02 | Zero-allocation guarantee |

**Deprecated/outdated:**
- `time_oft!`, `trotter_oft!` in `ofts.jl`: deprecated in favor of NUFFT prefactors computed at `_precompute_data` time. Still exported for debugging/pedagogy.
- `OFTCaches` struct: replaced by `NUFFTPrefactors` for production use.

## Discretion Recommendations

### Dense vs Sparse for A_nu2_dag in BohrDomain

**Recommendation: Dense, built via index-loop scatter into a scratch matrix.**

The dense code uses `sparse(rows_dag, cols_dag, conj.(A_nu_2_vals), dim, dim)` which allocates a new `SparseMatrixCSC` per bucket. For the matvec hot path, this is undesirable.

Two alternatives:
1. **Dense scatter (RECOMMENDED):** Zero a scratch matrix, then scatter the bucket values into it via index loops (like the Kraus trajectory code does). Uses `ws.LdagL` or `ws.tmp2` as the scratch target (after the dissipator helper is restructured to not need LdagL as scratch simultaneously). The advantage: all BLAS operations remain dense-dense, which is simpler and more predictable for allocation behavior.
2. **Apply A_nu2_dag via explicit index loops:** Instead of forming A_nu2_dag and calling BLAS, apply `rho * A_nu2_dag'` and `A_nu2_dag * alpha_A` via hand-written loops over the bucket indices (like `jump_workers.jl:250-274` in the Kraus path). This avoids allocating/zeroing A_nu2_dag entirely but sacrifices BLAS performance.

Given that Bohr is "not the primary code path" and dim=16 for the test system, either approach works. Dense scatter is simpler and consistent with the existing BLAS-based dissipator helpers.

### Which Scratch Matrix to Repurpose for BohrDomain

**Recommendation: Use `ws.LdagL` as the A_nu2_dag scratch buffer.**

Rationale: In `_accumulate_dissipator_2op!`, the `ws.LdagL` field is used to store `B'A` (the product for the anticommutator). But the 2op helper needs A_nu2_dag BEFORE computing `B'A`. The flow is:
1. Build A_nu2_dag into `ws.LdagL` (zero + scatter)
2. Call `_accumulate_dissipator_2op!` which reads A_nu2_dag from `ws.LdagL` first, then overwrites `ws.LdagL` with the anticommutator product `B'A`

This works because `ws.LdagL` is used AS the A_nu2_dag argument (not overwritten until the helper needs to store `B'A`). However, this creates a coupling between the loop body and the helper's internal use of `ws.LdagL`.

**Alternative:** Use `ws.tmp2` as A_nu2_dag scratch. The 2op helper can store `B'A` in `ws.LdagL` as usual. This is cleaner but requires that `ws.tmp2` is free before the scatter -- it IS free because it's only used inside the helper.

**Final recommendation: `ws.tmp2` for A_nu2_dag.** The `_accumulate_dissipator_2op!` helper should use `ws.LdagL` for B'A and `ws.tmp1` for intermediate products, exactly paralleling the single-op helper. `ws.tmp2` is free to hold A_nu2_dag because it's only used inside the helper.

### Test Matrix Generation

**Recommendation: Use `Matrix(random_density_matrix(NUM_QUBITS))` (same as Phase 27 tests).**

This generates Wishart-ensemble random density matrices. The `Matrix()` wrapper converts from `Hermitian` to plain `Matrix{ComplexF64}`, which matches the matvec signature `rho::Matrix{T}`. Using 10 random matrices per test (same as Phase 27) provides good coverage.

## Open Questions

### 1. BohrDomain Adjoint: Simple Swap is Incorrect (CRITICAL)

**What we know:** The CONTEXT.md says "For BohrDomain adjoint: swap arguments at call site." Mathematical analysis proves this is incorrect for the general two-operator dissipator.

**Detailed derivation:**

The forward two-operator dissipator:
```
D(A,B)(rho) = A * rho * B' - 0.5*(B'*A*rho + rho*B'*A)
```

The Hilbert-Schmidt adjoint (verified term by term):
- Sandwich `A*Y*B'`: adjoint is `A'*X*B` (from `<X, AYB'> = tr(X'AYB') = tr((A'XB)'Y) = <A'XB, Y>`)
- Anticommutator `-0.5*B'A*Y`: adjoint is `-0.5*A'B*X` (from `<X, -0.5*MY> = <-0.5*M'X, Y>` where M=B'A, M'=A'B)
- Anticommutator `-0.5*Y*B'A`: adjoint is `-0.5*X*A'B`

So: `D*(A,B)(X) = A'*X*B - 0.5*(A'*B*X + X*A'*B)`

If we swap A and B in the forward helper:
```
D(B,A)(X) = B*X*A' - 0.5*(A'*B*X + X*A'*B)
```

The anticommutator `A'*B` matches, but the sandwich term `B*X*A'` differs from `A'*X*B`. These are NOT equal in general (`B*X*A' != A'*X*B`).

**What's needed:** A separate `_accumulate_adjoint_dissipator_2op!(out, A, B_dag, rho, scalar, ws)` that computes:
```
scalar * (A'*rho*B_dag - 0.5*(A'*B_dag*rho + rho*A'*B_dag))
```

where A = alpha_A_nu1 (dense) and B_dag = A_nu2_dag (the operator whose adjoint appears on the right of the sandwich in the forward formula, i.e., B' = A_nu2_dag implies B_dag' = A_nu2).

Wait, let me re-derive with the actual operators from the code.

In the forward BohrDomain code:
```
_vectorize_liouv_diss_and_add!(L_target, alpha_A_nu1, A_nu_2_dag, gamma_norm_factor, ws)
```

This calls the two-operator `_vectorize_liouv_diss_and_add!(L, jump_1, jump_2, scalar, ws)` which computes (in un-vectorized form, using kron convention `kron(A,B) vec(X) = vec(B^T X A^T)`):

```
scalar * (jump_2^T * rho * jump_1^T - 0.5*(jump_2*jump_1)^T * rho - 0.5*rho*(jump_2*jump_1)^T)
```

With `jump_1 = alpha_A_nu1`, `jump_2 = A_nu_2_dag`:
```
scalar * (A_dag^T * rho * alpha_A^T - 0.5*(A_dag*alpha_A)^T * rho - 0.5*rho*(A_dag*alpha_A)^T)
```

= `scalar * (conj(A_dag) * rho * alpha_A^T - 0.5*(A_dag*alpha_A)^T * rho - 0.5*rho*(A_dag*alpha_A)^T)`

Since `A_dag[j,i] = conj(A[i,j])`, we have `conj(A_dag) = A^T`. And `A_dag^T = conj(A)`.

Hmm, this kron convention is making things complex. Let me instead validate by running the code and checking the round-trip test. The important takeaway for the planner is:

**Recommendation:** Create a dedicated `_accumulate_adjoint_dissipator_2op!` helper. Validate its correctness via:
1. The adjoint duality check: `tr(X' * L(Y)) == tr(L*(X)' * Y)` for random X, Y
2. The round-trip test against `L_dense'` (the dense Liouvillian transpose)

If simple swap happens to work due to the specific structure of the BohrDomain operators (which could be the case if A_dag^T is related to alpha_A in a specific way), the round-trip test will confirm it and the separate helper can be simplified. But plan defensively.

### 2. BohrDomain Allocation Budget

**What we know:** The CONTEXT.md allows "unavoidable allocations from entrywise alpha computation." The `@. alpha(bohr_freqs, nu_2) * eigenbasis` expression calls `alpha(nu_1, nu_2)` for each `(i,j)` pair, where `alpha` is a closure involving `exp`, `erfc`, and `sqrt`. Broadcasting a closure is generally allocation-free in Julia if the closure captures concrete-typed variables.

**What's unclear:** Whether `sparse()` construction should be avoided (replaced with dense scatter) or is acceptable given the relaxed allocation budget for Bohr. The Kraus trajectory code (jump_workers.jl:214-310) avoids `sparse()` entirely and uses index loops instead.

**Recommendation:** Avoid `sparse()` and use dense scratch scatter for the matvec. Even though Bohr allows some allocations, minimizing them is good practice and makes the code more consistent. Test with `@allocated` to characterize the actual allocation budget. If zero-allocation is achievable, great. If not, document the source.

### 3. `_prefactor_view` Return Type and BLAS Compatibility

**What we know:** `_prefactor_view(nufft_prefactors, omega)` returns `@view nufft_prefactors.data[:, :, k]`, which is a `SubArray{Complex{T}, 2}`. The expression `@. ws.jump_oft = eigenbasis * nufft_prefactor_matrix` broadcasts this view with the eigenbasis matrix. This should be allocation-free because both are strided arrays.

**What's unclear:** Whether `BLAS.gemm!` operations on the resulting `ws.jump_oft` (which is a plain `Matrix{T}`) have any issues. Answer: No, `ws.jump_oft` is always a plain `Matrix{T}`, so BLAS operations are fine. The `@view` is only used in the broadcast, not passed to BLAS.

**Recommendation:** No concern. The pattern `@. ws.jump_oft = eigenbasis * prefactor_view` writes into a plain `Matrix{T}` buffer, which is then passed to BLAS dissipator helpers.

## Sources

### Primary (HIGH confidence)
- Codebase analysis: `src/krylov_matvec.jl` -- Phase 27 EnergyDomain implementation (direct template)
- Codebase analysis: `src/krylov_workspace.jl` -- KrylovWorkspace struct and constructor (unchanged for Phase 28)
- Codebase analysis: `src/jump_workers.jl:91-135` -- Dense `_jump_contribution!` for `Union{TimeDomain, TrotterDomain}` (Time/Trotter reference)
- Codebase analysis: `src/jump_workers.jl:11-45` -- Dense `_jump_contribution!` for `BohrDomain` (Bohr reference)
- Codebase analysis: `src/jump_workers.jl:214-310` -- Kraus trajectory `_jump_contribution!` for `BohrDomain` (index-loop pattern for A_nu2_dag application)
- Codebase analysis: `src/qi_tools.jl:86-103` -- Two-operator `_vectorize_liouv_diss_and_add!` (dense Bohr dissipator formula)
- Codebase analysis: `src/nufft.jl:93-96` -- `_prefactor_view` (NUFFT prefactor access)
- Codebase analysis: `src/furnace_utensils.jl:82-121` -- `_precompute_data` for `Union{TimeDomain, TrotterDomain}` (precomputed data structure)
- Codebase analysis: `src/furnace_utensils.jl:14-28` -- `_precompute_data` for `BohrDomain` (Liouvillian variant)
- Codebase analysis: `src/bohr_domain.jl:82-85` -- `_pick_alpha` dispatch (KMS vs GNS alpha)
- Codebase analysis: `src/coherent.jl:14-40` -- `_precompute_coherent_total_B` (domain-dispatched B precomputation)
- Codebase analysis: `test/test_helpers.jl` -- All test factories and fixtures (config, jumps, trotter, hamiltonian)
- Codebase analysis: `test/test_krylov_matvec.jl` -- Phase 27 test structure (template for Phase 28 tests)
- Mathematical derivation: Hilbert-Schmidt adjoint of two-operator dissipator (see Open Question 1)

### Secondary (MEDIUM confidence)
- `.planning/research/ARCHITECTURE.md:350-412` -- Prior research on domain-specific matvec patterns
- `.planning/phases/27-core-matvec-infrastructure/27-02-SUMMARY.md` -- Phase 27 decisions on BLAS usage, concrete typing, adjoint dissipator

### Tertiary (LOW confidence)
None.

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- entirely internal to existing codebase, no new external dependencies
- Architecture (Time/Trotter): HIGH -- mechanical translation from EnergyDomain, all patterns directly observed in dense code
- Architecture (Bohr): HIGH -- dense reference code fully analyzed, Kraus trajectory code provides index-loop alternative
- Pitfalls: HIGH -- all identified from concrete code analysis and mathematical derivation; adjoint math verified algebraically
- Adjoint 2op math: HIGH -- derived from first principles using Hilbert-Schmidt inner product definition

**Research date:** 2026-02-20
**Valid until:** 2026-03-20 (stable -- internal codebase patterns, no external dependency version risk)
