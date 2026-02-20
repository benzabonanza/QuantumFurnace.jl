# Phase 27: Core Matvec Infrastructure - Research

**Researched:** 2026-02-20
**Domain:** Matrix-free Lindbladian action (superoperator-vector product) for quantum open systems
**Confidence:** HIGH

## Summary

Phase 27 implements matrix-free application of the Lindbladian superoperator L to a density matrix rho, producing L(rho) without forming the full dim^2 x dim^2 superoperator. This is the foundation for KrylovKit eigsolve integration (Phase 29). The codebase already has all the mathematical building blocks: `construct_lindbladian()` in `furnace.jl` builds the dense superoperator by calling `_jump_contribution!` in `jump_workers.jl`, which uses `_vectorize_liouv_diss_and_add!` to accumulate kron-product terms. The matvec equivalent replaces the kron-product vectorization with direct matrix operations: instead of `L_target += scalar * kron(A, conj(A))`, we compute `out += scalar * A * rho * A'`.

The existing `TrajectoryWorkspace` and `KrausScratch` structs demonstrate the project's workspace pattern: pre-allocate all scratch matrices at construction time, pass workspace to in-place functions, achieve zero allocations on the hot path. The new `KrylovWorkspace` follows this pattern exactly. The dense `_jump_contribution!` for `AbstractLiouvConfig{EnergyDomain}` (jump_workers.jl:47-89) is the direct template: same loop structure (sites outer, frequencies inner), same `oft!()` call to build A(omega), same Hermitian half-grid optimization, same prefactor computation. The matvec just applies A(omega) to rho directly instead of vectorizing into a superoperator matrix.

**Primary recommendation:** Implement `apply_lindbladian!` and `apply_adjoint_lindbladian!` as thin wrappers around a shared dissipator accumulation loop that mirrors `_jump_contribution!(L_target, jump, hamiltonian, config::AbstractLiouvConfig{EnergyDomain}, ...)` line-by-line, replacing `_vectorize_liouv_diss_and_add!(L_target, jump_oft, scalar, ws)` with `out += scalar * (A*rho*A' - 0.5*(A'A*rho + rho*A'A))` using pre-allocated scratch matrices.

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

#### Adjoint interface
- Separate function `apply_adjoint_lindbladian!` (not keyword flag) -- explicit, no ambiguity
- Shared KrylovWorkspace between forward and adjoint -- same scratch matrices, used differently
- Leverage single-site Hermitian Pauli structure: A^a are Hermitian, so A^a' is lazy/free. Kraus jumps A(omega) are not Hermitian but have symmetry with negative omegas -- use existing patterns from DM simulator code
- Memory-efficient: use lazy `'` adjoint where possible, don't store explicit adjoints of Hermitian operators
- Design signatures now for KrylovKit compatibility -- `apply_lindbladian!` should be easily wrappable into a closure `rho -> L(rho)` for KrylovKit's `eigsolve`

#### Workspace scope
- KrylovWorkspace designed for all domains from the start (Energy, Time, Trotter, Bohr) -- Phase 28 adds dispatch methods, not workspace restructuring
- Tied to (config, ham) at construction -- precomputes everything once, one workspace per problem instance
- Scratch matrix pattern: follow existing DM simulator workspace pattern (TrajectoryWorkspace or equivalent)
- Jump operator storage: match DM simulator -- dense matrices for A(omega) in EnergyDomain; for Time/TrotterDomain store NUFFT matrices and create A(omega) on the fly via entrywise multiplication

#### Coherent term control
- `with_coherent` read from config object (not keyword arg) -- config already carries this information
- Coherent correction B precomputed at workspace construction, using existing functions that precompute full B and its action (or B^a and its action)
- GNS balance: silently skip coherent term (no error if GNS config encountered)
- Adjoint coherent: automatic sign flip -- `apply_adjoint_lindbladian!` uses +i[B, rho] instead of -i[B, rho], same B matrix. B is approximately Hermitian (up to Trotter and quadrature errors)

#### Dispatch pattern
- Dispatch on domain type via config parametrization: `apply_lindbladian!(... config::SamplerConfig{EnergyDomain} ...)` -- matches existing DM simulator dispatch
- KMS vs GNS dispatch: follow DM simulator pattern with different config names (one with GNS, one without)
- Inner loop order: match existing `construct_lindbladian` pattern (sites outer, frequencies inner) -- consistency, easier to validate
- Output convention: determined by KrylovKit eigsolve requirements -- research phase should check what signature KrylovKit expects and design accordingly

### Claude's Discretion
- Exact number of scratch matrices in the workspace pool
- Cache-friendliness optimizations within the loop (as long as loop order matches existing pattern)
- Internal helper function decomposition
- Test matrix generation strategy for round-trip validation

### Deferred Ideas (OUT OF SCOPE)
None -- discussion stayed within phase scope
</user_constraints>

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| LinearAlgebra (stdlib) | Julia 1.x | `mul!`, `adjoint`, `Hermitian`, `rmul!` | Already used throughout codebase; all matrix operations use BLAS via this |
| QuantumFurnace internals | current | `oft!`, `_precompute_data`, `_precompute_coherent_total_B`, `pick_transition`, workspace patterns | The reference implementation we mirror |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| KrylovKit.jl | latest (not yet added) | `eigsolve` consumer of the matvec closure | Phase 29 adds dependency; Phase 27 designs the interface for it |
| Test (stdlib) | Julia 1.x | `@testset`, `@test`, `@allocated` | Round-trip validation and allocation regression tests |

### Alternatives Considered
Not applicable -- all decisions locked. No external libraries needed for Phase 27 itself.

**Installation:**
No new dependencies for Phase 27. KrylovKit.jl is added in Phase 29.

## Architecture Patterns

### Recommended File Structure
```
src/
  krylov_workspace.jl       # KrylovWorkspace struct + constructor
  krylov_matvec.jl          # apply_lindbladian!, apply_adjoint_lindbladian!
  structs.jl                # (existing) domain types, config types
  ...
test/
  test_krylov_matvec.jl     # Round-trip tests, allocation tests
```

### Pattern 1: Workspace Pre-allocation (from TrajectoryWorkspace)

**What:** All scratch matrices allocated once at construction, stored in a typed struct, passed to every in-place function call.

**When to use:** Any hot-path function that will be called O(thousands) of times (KrylovKit calls the matvec O(krylovdim * maxiter) times).

**Reference implementation (`src/trajectories.jl:1-18`):**
```julia
struct TrajectoryWorkspace{T}
    jump_oft::Matrix{T}   # buffer for A_omega
    psi_tmp::Vector{T}    # generic tmp for matvec results
    Rpsi::Vector{T}       # cache R*psi
    rho_acc::Matrix{T}    # density matrix accumulator
end

function TrajectoryWorkspace(::Type{T}, dim::Int) where {T}
    return TrajectoryWorkspace{T}(
        zeros(T, dim, dim),
        zeros(T, dim),
        zeros(T, dim),
        zeros(T, dim, dim),
    )
end
```

**KrylovWorkspace design (recommended):**
```julia
struct KrylovWorkspace{T<:Complex}
    # Precomputed data (immutable after construction)
    precomputed_data::Any          # NamedTuple from _precompute_data()
    B_total::Union{Nothing, Matrix{T}}  # Precomputed coherent B (nothing for GNS)

    # Scratch matrices for dissipator accumulation (dim x dim)
    jump_oft::Matrix{T}       # A(omega) buffer
    tmp1::Matrix{T}           # scratch for A*rho, rho*A', etc.
    tmp2::Matrix{T}           # scratch for A'*A, etc.
    LdagL::Matrix{T}          # scratch for L'L accumulation
    rho_out::Matrix{T}        # output accumulator (zeroed at start of each matvec)
end
```

### Pattern 2: Dense Dissipator Loop (from `_jump_contribution!` for EnergyDomain Liouvillian)

**What:** The dense Liouvillian builder iterates: for each jump operator, for each energy label (half-grid for Hermitian jumps), compute A(omega) via `oft!()`, then accumulate the dissipator contribution. The matvec follows the identical loop but applies the dissipator formula directly to rho.

**Dense reference (`src/jump_workers.jl:47-89`):**
```julia
# Dense: accumulates into dim^2 x dim^2 superoperator matrix
prefactor = (config.w0 / (config.sigma * sqrt(2*pi))) * gamma_norm_factor
for w_raw in energy_labels
    w_raw > 1e-12 && continue
    w = abs(w_raw)
    oft!(jump_oft, jump, w, hamiltonian, config.sigma)
    scalar_w = prefactor * transition(w)
    # Vectorized dissipator: kron(A, conj(A)) - 0.5*kron(A'A, I) - 0.5*kron(I, (A'A)^T)
    _vectorize_liouv_diss_and_add!(L_target, jump_oft, scalar_w, ws)
    if w > 1e-12  # negative freq partner: L = A(w)', i.e. adjoint
        scalar_neg = prefactor * transition(-w)
        _vectorize_liouv_diss_and_add!(L_target, jump_oft', scalar_neg, ws)
    end
end
```

**Matvec equivalent:**
```julia
# Matvec: applies dissipator directly to rho, accumulates into rho_out
prefactor = (config.w0 / (config.sigma * sqrt(2*pi))) * gamma_norm_factor
for w_raw in energy_labels
    w_raw > 1e-12 && continue
    w = abs(w_raw)
    oft!(ws.jump_oft, jump, w, hamiltonian, config.sigma)
    scalar_w = prefactor * transition(w)
    _accumulate_dissipator!(ws.rho_out, ws.jump_oft, rho, scalar_w, ws)
    if w > 1e-12
        scalar_neg = prefactor * transition(-w)
        _accumulate_dissipator!(ws.rho_out, ws.jump_oft', rho, scalar_neg, ws)
    end
end
```

Where `_accumulate_dissipator!` computes: `out += scalar * (L * rho * L' - 0.5 * (L'L * rho + rho * L'L))`

### Pattern 3: Dispatch on Domain Type via Config Parametrization

**What:** The project dispatches domain-specific behavior via the `D` type parameter in `AbstractLiouvConfig{D}`. The same function name has multiple methods for different domains.

**Reference (`src/jump_workers.jl`):**
```julia
function _jump_contribution!(L, jump, ham, config::AbstractLiouvConfig{EnergyDomain}, ...)
function _jump_contribution!(L, jump, ham, config::AbstractLiouvConfig{BohrDomain}, ...)
function _jump_contribution!(L, jump, ham, config::AbstractLiouvConfig{D}, ...) where D<:Union{TimeDomain,TrotterDomain}
```

**For Phase 27 (EnergyDomain only):**
```julia
function apply_lindbladian!(ws::KrylovWorkspace, rho, config::AbstractLiouvConfig{EnergyDomain}, ham)
```

Phase 28 adds:
```julia
function apply_lindbladian!(ws::KrylovWorkspace, rho, config::AbstractLiouvConfig{D}, ham) where D<:Union{TimeDomain,TrotterDomain}
function apply_lindbladian!(ws::KrylovWorkspace, rho, config::AbstractLiouvConfig{BohrDomain}, ham)
```

### Pattern 4: KMS vs GNS Dispatch

**What:** GNS configs (`LiouvConfigGNS`) always have `with_coherent=false` and use a different transition function. The existing code uses separate config types (`LiouvConfig` vs `LiouvConfigGNS`) that both are subtypes of `AbstractLiouvConfig`. Dispatch on the config type happens via `pick_transition()` and the `with_coherent` check.

**Reference (`src/energy_domain.jl:3-6`):**
```julia
pick_transition(config::LiouvConfig) = _pick_transition_kms(config)
pick_transition(config::LiouvConfigGNS) = _pick_transition_gns(config)
```

**For Phase 27:** The matvec function accepts `config::AbstractLiouvConfig{EnergyDomain}`, which naturally covers both KMS and GNS. The `_precompute_data()` and `pick_transition()` functions already handle the dispatch. The coherent term is skipped when `config.with_coherent == false` (which includes all GNS configs).

### Pattern 5: KrylovKit Closure Interface

**What:** KrylovKit's `eigsolve` accepts `f(x) -> y` (not `f!(y, x)`). The function takes a vector and returns a vector. For density matrices, the "vector" is `vec(rho)` (a `Vector{ComplexF64}` of length dim^2).

**KrylovKit API (from official docs):**
```julia
vals, vecs, info = eigsolve(f, x0, howmany, :LR; issymmetric=false, krylovdim=50, tol=1e-10)
```

Where `f` is a function `Vector{ComplexF64} -> Vector{ComplexF64}`.

**Closure design (Phase 29 wraps Phase 27):**
```julia
# Phase 29 will create this closure:
function make_lindbladian_map(ws::KrylovWorkspace, config, ham)
    dim = size(ham.data, 1)
    return v -> begin
        rho = reshape(v, dim, dim)
        apply_lindbladian!(ws, rho, config, ham)
        return vec(ws.rho_out)
    end
end
```

**Critical note:** KrylovKit calls `f(x)` and expects a NEW vector as output (not mutation of x). The closure must return `vec(ws.rho_out)` (or a copy if KrylovKit might store it). Investigation of KrylovKit source shows it stores the output vectors in the Krylov basis, so the closure should return `copy(vec(ws.rho_out))` or, better, use a fresh output each time. However, allocating a new vector per call defeats zero-allocation. The standard approach is to let KrylovKit handle vector storage -- the returned `vec(ws.rho_out)` will be copied by KrylovKit into its internal Krylov basis. The workspace `rho_out` is overwritten on each call, which is fine because KrylovKit copies the return value.

**Alternative (avoiding vec/reshape overhead):** KrylovKit supports custom vector types. The density matrix itself could serve as the "vector" if we define `LinearAlgebra.mul!`, `LinearAlgebra.rmul!`, `LinearAlgebra.axpy!`, `LinearAlgebra.axpby!`, `LinearAlgebra.dot`, and `LinearAlgebra.norm` for it. However, using `vec()` with plain `Vector{ComplexF64}` is simpler and the reshape/vec are zero-cost (views, not copies). **Recommend: use vec/reshape approach.** This is what the abandoned `linearmaps_liouv.jl` code used.

### Anti-Patterns to Avoid
- **Allocating inside the loop:** Never call `zeros()`, `similar()`, or `copy()` inside the energy-label loop. All scratch must be pre-allocated in KrylovWorkspace.
- **Building the superoperator matrix:** The whole point is to avoid forming dim^2 x dim^2. Never call `_vectorize_liouv_diss_and_add!` or `_kron!` in the matvec path.
- **Redundant eigendecompositions:** Unlike `_finalize_kraus_step!` (which does `eigen()` for PSD clamping), the Lindbladian matvec has no eigen() calls. The matvec is purely matrix multiply + add.
- **Storing explicit adjoints of Hermitian operators:** For Hermitian `jump.in_eigenbasis`, `jump.in_eigenbasis'` is lazy (Julia returns an `Adjoint` wrapper, not a copy). Use this directly in `mul!`.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Dissipator formula application | Manual index loops for L*rho*L' | `mul!()` chaining with scratch matrices | BLAS3 is 10-100x faster than manual loops for dim=16+ |
| Energy-domain jump A(omega) | Custom Gaussian filter | Existing `oft!(out, jump, w, ham, sigma)` | Already optimized, tested, allocation-free |
| Coherent B precomputation | Custom B calculation | Existing `_precompute_coherent_total_B(jumps, ham, config, precomputed_data)` | Handles all domains, all balance types, already scaled |
| Precomputed data (transition, labels, gamma_norm) | Manual precomputation | Existing `_precompute_data(config, ham)` | Already handles all domain x balance-type combinations |
| Config validation | Manual checks | Existing `validate_config!(config)` | Comprehensive validation of all config fields |
| Half-grid optimization for Hermitian jumps | Custom frequency filtering | Existing `jump.hermitian` flag + `w_raw > 1e-12 && continue` pattern | Exact pattern from dense Liouvillian builder, halves iteration count |

**Key insight:** Phase 27 is NOT building new physics. It is translating the existing dense vectorization formula `L_target += scalar * kron(A, conj(A)) - ...` into the equivalent matrix-on-matrix formula `out += scalar * (A * rho * A' - ...)`. Every computational ingredient already exists.

## Common Pitfalls

### Pitfall 1: Sign Convention in Coherent Term
**What goes wrong:** The coherent term is `-i[B, rho] = -i*B*rho + i*rho*B` for the forward Lindbladian, and `+i[B, rho]` for the adjoint. Getting the sign wrong produces a Lindbladian that doesn't match the dense reference.
**Why it happens:** The dense vectorization uses `_kron!(L, B, Id, -1im)` and `_kron!(L, Id, B^T, +1im)`, which translates to `-i*B*rho + i*rho*B`. Easy to flip signs when translating.
**How to avoid:** Write the coherent term as a separate helper. Verify: `coherent_contribution = -1im * (B * rho - rho * B)` for forward, `+1im * (B * rho - rho * B)` for adjoint.
**Warning signs:** Round-trip test fails only when `with_coherent=true` but passes when `with_coherent=false`.

### Pitfall 2: Adjoint Dissipator -- Wrong Swap Direction
**What goes wrong:** The adjoint Lindbladian swaps L with L' in the dissipator: `D*(rho) = sum L' rho L - 0.5*(L L' rho + rho L L')`. Note that the anticommutator uses `L L'` (not `L' L`). Getting this wrong produces wrong adjoint.
**Why it happens:** The adjoint of `D(rho) = L rho L' - 0.5*(L'L rho + rho L'L)` is `D*(rho) = L' rho L - 0.5*(L L' rho + rho L L')`. The swap is symmetric: wherever you see L, put L', and vice versa.
**How to avoid:** Factor out the dissipator application into a helper `_accumulate_dissipator!(out, L_op, rho, scalar, ws)` that takes the Lindblad operator as argument. For forward: pass `L_op = A(omega)`. For adjoint: pass `L_op = A(omega)'`. The helper applies `out += scalar * (L_op * rho * L_op' - 0.5*(L_op' * L_op * rho + rho * L_op' * L_op))` regardless.
**Warning signs:** Adjoint round-trip test fails. Check by computing `tr(X' * L(Y))` vs `tr(L*(X)' * Y)` for random X, Y.

### Pitfall 3: Hermitian Half-Grid Negative Frequency Partner
**What goes wrong:** For Hermitian jumps, the negative-frequency contribution uses `L = A(omega)'` (the adjoint). When the jump is Hermitian, `A^a' = A^a` (lazy, free). But `A(omega)` is NOT Hermitian (it's `A^a` element-wise multiplied by a Gaussian filter that breaks Hermiticity). So `A(omega)' != A(omega)`.
**Why it happens:** Confusion between the Hermiticity of the single-site Pauli operator A^a and the frequency-projected jump A(omega).
**How to avoid:** For the negative-frequency partner at `-omega`: the Lindblad operator is `A(omega)'` (adjoint of the positive-frequency A(omega)). This is implemented via `jump_oft'` (lazy adjoint view) in the existing dense code. Use the same `ws.jump_oft'` in the matvec.
**Warning signs:** Round-trip test fails for KMS configs (which have non-trivial transition(-w)).

### Pitfall 4: GNS Transition vs KMS Transition Mismatch
**What goes wrong:** Using the wrong transition function for GNS configs produces wrong prefactors.
**Why it happens:** `pick_transition(config::LiouvConfig)` returns the KMS-shifted transition, while `pick_transition(config::LiouvConfigGNS)` returns the unshifted GNS transition. The transition function is part of `precomputed_data`, so as long as `_precompute_data` is called correctly, this is handled.
**How to avoid:** Call `_precompute_data(config, ham)` once at workspace construction. The transition function is stored in `precomputed_data.transition` and is correct for whatever config type was passed.
**Warning signs:** KMS round-trip passes but GNS round-trip fails (or vice versa).

### Pitfall 5: Forgetting to Zero the Output Accumulator
**What goes wrong:** If `ws.rho_out` is not zeroed at the start of each `apply_lindbladian!` call, the output accumulates across multiple KrylovKit iterations.
**Why it happens:** The workspace is reused across calls. KrylovKit calls the matvec hundreds of times.
**How to avoid:** `fill!(ws.rho_out, 0)` at the very start of `apply_lindbladian!`, before any accumulation.
**Warning signs:** First KrylovKit eigenvalue is correct but subsequent ones are garbage.

### Pitfall 6: mul! with Adjoint Views and Beta Accumulation
**What goes wrong:** `mul!(C, A, B', alpha, beta)` (5-arg mul!) may allocate if `B'` is a lazy `Adjoint` wrapper and BLAS cannot dispatch to `gemm` directly.
**Why it happens:** Julia's `mul!` dispatches to BLAS `gemm` for `Matrix * Adjoint{Matrix}` but the path may go through a fallback that allocates.
**How to avoid:** Use the 3-arg `mul!(tmp, A, B')` followed by `rho_out .+= scalar .* tmp`, or use `mul!(rho_out, A, B', scalar, one(eltype(rho_out)))` which chains into BLAS gemm correctly for `StridedMatrix` types. Test with `@allocated`.
**Warning signs:** `@allocated` shows non-zero allocations in the matvec hot path.

## Code Examples

### Example 1: Dissipator Helper (Core Building Block)
```julia
"""
    _accumulate_dissipator!(out, L_op, rho, scalar, ws) -> nothing

Accumulate scalar * D_L(rho) into `out`:
    D_L(rho) = L * rho * L' - 0.5 * (L'L * rho + rho * L'L)

Uses ws.tmp1, ws.tmp2, ws.LdagL as scratch.
"""
function _accumulate_dissipator!(
    out::Matrix{<:Complex},
    L_op::AbstractMatrix{<:Complex},
    rho::AbstractMatrix{<:Complex},
    scalar::Real,
    ws,  # KrylovWorkspace or similar with tmp1, tmp2, LdagL fields
)
    # L'L -> ws.LdagL
    mul!(ws.LdagL, L_op', L_op)

    # Term 1: L * rho * L'
    mul!(ws.tmp1, L_op, rho)         # tmp1 = L * rho
    mul!(ws.tmp2, ws.tmp1, L_op')    # tmp2 = L * rho * L'
    @. out += scalar * ws.tmp2

    # Term 2: -0.5 * L'L * rho
    mul!(ws.tmp1, ws.LdagL, rho)     # tmp1 = L'L * rho
    @. out -= 0.5 * scalar * ws.tmp1

    # Term 3: -0.5 * rho * L'L
    mul!(ws.tmp1, rho, ws.LdagL)     # tmp1 = rho * L'L
    @. out -= 0.5 * scalar * ws.tmp1

    return nothing
end
```

### Example 2: Forward EnergyDomain Matvec (Skeleton)
```julia
function apply_lindbladian!(
    ws::KrylovWorkspace{CT},
    rho::AbstractMatrix{<:Complex},
    config::AbstractLiouvConfig{EnergyDomain},
    hamiltonian::HamHam,
) where {CT}
    (; transition, gamma_norm_factor, energy_labels) = ws.precomputed_data

    fill!(ws.rho_out, 0)

    # Coherent term: -i[B, rho]
    B = ws.B_total
    if B !== nothing
        # out += -i*B*rho + i*rho*B
        mul!(ws.tmp1, B, rho)
        @. ws.rho_out += -1im * ws.tmp1
        mul!(ws.tmp1, rho, B)
        @. ws.rho_out += 1im * ws.tmp1
    end

    # Dissipator: sum over jumps, sum over energy labels
    prefactor = (config.w0 / (config.sigma * sqrt(2 * pi))) * gamma_norm_factor

    for jump in ws.jumps
        if jump.hermitian
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)
                oft!(ws.jump_oft, jump, w, hamiltonian, config.sigma)

                scalar_w = prefactor * transition(w)
                _accumulate_dissipator!(ws.rho_out, ws.jump_oft, rho, scalar_w, ws)

                if w > 1e-12
                    scalar_neg = prefactor * transition(-w)
                    _accumulate_dissipator!(ws.rho_out, ws.jump_oft', rho, scalar_neg, ws)
                end
            end
        else
            for w in energy_labels
                oft!(ws.jump_oft, jump, w, hamiltonian, config.sigma)
                scalar_w = prefactor * transition(w)
                _accumulate_dissipator!(ws.rho_out, ws.jump_oft, rho, scalar_w, ws)
            end
        end
    end

    return ws.rho_out
end
```

### Example 3: Adjoint Matvec (Differences from Forward)
```julia
function apply_adjoint_lindbladian!(
    ws::KrylovWorkspace{CT},
    rho::AbstractMatrix{<:Complex},
    config::AbstractLiouvConfig{EnergyDomain},
    hamiltonian::HamHam,
) where {CT}
    (; transition, gamma_norm_factor, energy_labels) = ws.precomputed_data

    fill!(ws.rho_out, 0)

    # Adjoint coherent term: +i[B, rho] (sign flip vs forward)
    B = ws.B_total
    if B !== nothing
        mul!(ws.tmp1, B, rho)
        @. ws.rho_out += 1im * ws.tmp1   # +i instead of -i
        mul!(ws.tmp1, rho, B)
        @. ws.rho_out -= 1im * ws.tmp1   # -i instead of +i
    end

    # Adjoint dissipator: swap L <-> L'
    # D*(rho) = L' rho L - 0.5*(L L' rho + rho L L')
    prefactor = (config.w0 / (config.sigma * sqrt(2 * pi))) * gamma_norm_factor

    for jump in ws.jumps
        if jump.hermitian
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)
                oft!(ws.jump_oft, jump, w, hamiltonian, config.sigma)

                scalar_w = prefactor * transition(w)
                # Forward uses L=A(w), adjoint uses L=A(w)' (i.e. pass A(w)' as L_op)
                _accumulate_dissipator!(ws.rho_out, ws.jump_oft', rho, scalar_w, ws)

                if w > 1e-12
                    scalar_neg = prefactor * transition(-w)
                    # Forward uses L=A(w)', adjoint uses L=(A(w)')' = A(w)
                    _accumulate_dissipator!(ws.rho_out, ws.jump_oft, rho, scalar_neg, ws)
                end
            end
        else
            for w in energy_labels
                oft!(ws.jump_oft, jump, w, hamiltonian, config.sigma)
                scalar_w = prefactor * transition(w)
                _accumulate_dissipator!(ws.rho_out, ws.jump_oft', rho, scalar_w, ws)
            end
        end
    end

    return ws.rho_out
end
```

### Example 4: Round-Trip Test Pattern
```julia
@testset "Round-trip: matvec vs dense" begin
    config = make_liouv_config(EnergyDomain())
    L_dense = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)
    ws = KrylovWorkspace(config, TEST_HAM, TEST_JUMPS)

    for _ in 1:10
        rho = random_density_matrix(DIM)

        # Dense reference
        v_dense = L_dense * vec(rho)

        # Matvec
        L_rho = apply_lindbladian!(ws, rho, config, TEST_HAM)
        v_matvec = vec(L_rho)

        @test norm(v_dense - v_matvec) < 1e-12
    end
end
```

### Example 5: Allocation Test Pattern (from existing test_allocation.jl)
```julia
@testset "apply_lindbladian! zero allocations" begin
    config = make_liouv_config(EnergyDomain())
    ws = KrylovWorkspace(config, TEST_HAM, TEST_JUMPS)
    rho = random_density_matrix(DIM)

    function _measure_matvec_allocs(ws, rho, config, ham)
        # Warmup
        apply_lindbladian!(ws, rho, config, ham)
        # Measure
        return @allocated apply_lindbladian!(ws, rho, config, ham)
    end

    allocs = _measure_matvec_allocs(ws, rho, config, TEST_HAM)
    @test allocs == 0
end
```

## Scratch Matrix Count Analysis (Claude's Discretion)

The dissipator formula `D_L(rho) = L*rho*L' - 0.5*(L'L*rho + rho*L'L)` requires:

**Minimum scratch matrices needed (dim x dim):**
1. `jump_oft` -- A(omega) buffer (written by `oft!`)
2. `LdagL` -- L'*L product (reused across anticommutator terms)
3. `tmp1` -- general scratch (L*rho, L'L*rho, rho*L'L, etc.)
4. `tmp2` -- second scratch (L*rho*L' before accumulation)
5. `rho_out` -- output accumulator

**Total: 5 scratch matrices** for the dissipator loop.

The coherent term `-i[B, rho]` reuses `tmp1` (for B*rho and rho*B), so no additional scratch needed.

**Comparison with existing workspaces:**
- `LindbladianWorkspace`: 5 matrices (Id, jump_tmp, jump_conj, jump_dag_jump, jump2_jump1) -- for dense vectorization
- `KrausScratch`: 8 matrices -- more complex due to K0, R accumulation, rho_jump, etc.
- `TrajectoryWorkspace`: 2 matrices + 2 vectors -- simpler (state-vector ops, not density-matrix)

**Recommendation: 5 scratch matrices.** This is the minimum needed for zero-allocation dissipator application. No Id matrix is needed (the matvec doesn't use kron products). The B_total matrix is stored separately as precomputed data, not counted as scratch.

## Adjoint Lindbladian: Mathematical Details

The Lindbladian L and its adjoint L* (Hilbert-Schmidt adjoint) satisfy:
```
tr(X' * L(Y)) = tr(L*(X)' * Y)    for all X, Y
```

If the forward Lindbladian is:
```
L(rho) = -i[B, rho] + sum_k gamma_k * (L_k * rho * L_k' - 0.5*{L_k'*L_k, rho})
```

Then the adjoint is:
```
L*(rho) = +i[B, rho] + sum_k gamma_k * (L_k' * rho * L_k - 0.5*{L_k*L_k', rho})
```

**Key differences:**
1. Coherent term sign flips: `-i[B, .]` becomes `+i[B, .]`
2. Dissipator swaps L_k with L_k': `L*rho*L'` becomes `L'*rho*L`, and `L'L` becomes `LL'`

**For Hermitian A^a:** `(A^a)' = A^a`, so `A(omega)' = (A^a o G(omega))' = (A^a)' o conj(G(omega)) = A^a o conj(G(omega))`. The `oft!` function computes `A(omega) = A^a o G(omega)` where G is real and positive (Gaussian filter), so `A(omega)' = A^a o G(omega) = A(omega)` only when the Bohr frequencies are zero. In general `A(omega)' != A(omega)`.

**Implementation:** The `_accumulate_dissipator!` helper takes the Lindblad operator as an argument. For the adjoint, we simply pass the adjoint of what we'd pass for the forward. This is clean because Julia's `'` on a matrix returns a lazy `Adjoint` wrapper.

## KrylovWorkspace Constructor Design

**Precomputation at construction (tied to config + ham):**
1. `_precompute_data(config, ham)` -> `precomputed_data` (energy_labels, transition, gamma_norm_factor, and domain-specific data like NUFFT prefactors)
2. `_precompute_coherent_total_B(jumps, ham, config, precomputed_data)` -> `B_total` (or nothing for GNS/with_coherent=false)
3. Store a reference to the jump operators (needed in the matvec loop)
4. Allocate all scratch matrices once

**Constructor signature:**
```julia
function KrylovWorkspace(
    config::AbstractLiouvConfig,
    hamiltonian::HamHam,
    jumps::Vector{JumpOp};
    trotter::Union{TrottTrott, Nothing}=nothing,
)
```

This mirrors `_build_framework_and_seed()` from trajectories.jl but for the Krylov context. The workspace is designed for all domains from the start -- Phase 28 just adds dispatch methods for `apply_lindbladian!`.

## Open Questions

1. **`mul!` with `Adjoint` views and zero-allocation guarantee**
   - What we know: Julia's `mul!(C, A, B')` dispatches to BLAS `gemm` with `transB='C'` for strided matrices, which should be allocation-free. The 5-arg `mul!(C, A, B', alpha, beta)` also works.
   - What's unclear: Whether passing `jump_oft'` (lazy Adjoint wrapper) to `_accumulate_dissipator!` as the `L_op` argument causes any allocation in the `mul!(ws.LdagL, L_op', L_op)` call (which becomes `mul!(ws.LdagL, (A')', A')` = `mul!(ws.LdagL, A, A')`). This double-adjoint should resolve to the original matrix, but need to verify with `@allocated`.
   - Recommendation: Test with `@allocated` during implementation. If double-adjoint causes allocation, unwrap manually.

2. **Random density matrix generation for tests**
   - What we know: `random_density_matrix(num_qubits)` exists in `qi_tools.jl` and returns a `Hermitian` matrix. It takes `num_qubits` (not dim), using Wishart ensemble (A*A'/tr(A*A')).
   - Resolved: Use `random_density_matrix(NUM_QUBITS)` in tests. Note the argument is num_qubits, not dim. The function returns `Hermitian` -- may need `Matrix()` wrapper if the matvec function signature requires `Matrix{Complex}`.

3. **Whether `vec(ws.rho_out)` allocates**
   - What we know: `vec()` on a `Matrix` returns a `Base.ReshapedArray` (a view, not a copy) in Julia. This should be allocation-free.
   - What's unclear: Whether KrylovKit's internal operations (axpy!, norm, dot) work on `ReshapedArray` without allocating.
   - Recommendation: This is a Phase 29 concern. For Phase 27, the matvec returns `ws.rho_out` (a Matrix). Phase 29 wraps it in the vec/reshape closure.

## Sources

### Primary (HIGH confidence)
- Codebase analysis: `src/jump_workers.jl` -- dense `_jump_contribution!` for all domains and config types (reference loop structure)
- Codebase analysis: `src/qi_tools.jl` -- `_vectorize_liouv_diss_and_add!` and `_kron!` (dense vectorization formula that matvec replaces)
- Codebase analysis: `src/trajectories.jl` -- `TrajectoryWorkspace`, `TrajectoryFramework`, `build_trajectoryframework()` (workspace pattern reference)
- Codebase analysis: `src/kraus.jl` -- `KrausScratch` (alternative workspace pattern reference)
- Codebase analysis: `src/coherent.jl` -- `_precompute_coherent_total_B()` (coherent B precomputation, reused directly)
- Codebase analysis: `src/furnace_utensils.jl` -- `_precompute_data()` for all domains (precomputed data, reused directly)
- Codebase analysis: `src/energy_domain.jl` -- `pick_transition()`, KMS vs GNS dispatch
- Codebase analysis: `src/ofts.jl` -- `oft!()` (jump operator frequency projection, reused directly)
- Codebase analysis: `src/structs.jl` -- domain types, config types, workspace struct patterns
- Codebase analysis: `src/furnace.jl` -- `construct_lindbladian()` (dense reference for round-trip test)
- Codebase analysis: `test/test_allocation.jl` -- zero-allocation test pattern
- Codebase analysis: `test/test_helpers.jl` -- test system setup, config factories

### Secondary (MEDIUM confidence)
- [KrylovKit.jl official documentation](https://jutho.github.io/KrylovKit.jl/stable/man/eig/) -- eigsolve API signature, function-based interface `f(x) -> y`
- [KrylovKit.jl GitHub](https://github.com/Jutho/KrylovKit.jl) -- confirms `eigsolve` takes callable, returns (vals, vecs, info)
- `.planning/research/FEATURES.md` -- prior research on KrylovKit integration strategy, eigenvalue targeting

### Tertiary (LOW confidence)
None.

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- entirely internal to existing codebase, no new external dependencies
- Architecture: HIGH -- direct translation of existing dense patterns into matvec form, reference implementations fully analyzed
- Pitfalls: HIGH -- all identified from concrete code analysis, not theoretical concerns
- KrylovKit interface: MEDIUM -- verified from official docs, but closure allocation behavior needs empirical validation in Phase 29

**Research date:** 2026-02-20
**Valid until:** 2026-03-20 (stable -- internal codebase patterns, no external dependency version risk)
