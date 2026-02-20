# Architecture Patterns: Krylov-based Lindbladian Spectral Gap Estimation

**Domain:** Matrix-free Krylov eigensolving for Lindbladian spectral gap in QuantumFurnace.jl
**Researched:** 2026-02-20
**Confidence:** HIGH (direct analysis of all 26 source files, existing architecture patterns fully traced, KrylovKit API verified via official docs)

## Recommended Architecture

### Design Principle: The Lindbladian Action, Not the Channel, Not the Matrix

There are three options for what to wrap as a KrylovKit linear map:

1. **Full Liouvillian matrix L** -- build the dim^2 x dim^2 matrix, pass to KrylovKit as an `AbstractMatrix`. Defeats the purpose of matrix-free methods.

2. **CPTP channel E_delta** -- wrap the existing DM thermalization step (`_jump_contribution!` for `AbstractThermalizeConfig`). Eigenvalues cluster around 1.0 (since E_delta(rho) ~ rho + delta*L(rho)), making Krylov convergence poor.

3. **Lindbladian generator L** -- compute L(rho) directly from the dissipator formula without forming the full matrix. Eigenvalues have Re(lambda) <= 0; spectral gap = |Re(lambda_2)|.

**Recommendation: Option 3 -- wrap the Lindbladian generator L as a matrix-free action.** This is the mathematically clean approach. KrylovKit's Arnoldi with `:SR` (smallest real part, i.e., most negative) directly targets the spectral gap. The existing `_precompute_data()` and coherent B infrastructure is reused verbatim; only the per-iteration action needs new code.

### High-Level Architecture

```
User API                          Internal                              KrylovKit
---------                         --------                              ---------

krylov_spectral_gap(              build_lindbladian_action(             eigsolve(
  jumps, config, hamiltonian;       jumps, config, hamiltonian;           lindbladian_action,
  trotter=nothing,                  trotter=nothing                       x0,
  n_eigenvalues=6,              ) -> function f(v) -> w                   howmany,
  krylov_dim=30,                     |                                    which=:SR,
  tol=1e-10                          |  (closure captures                 tol=tol
) -> KrylovGapResult                 |   precomputed_data,              )
                                     |   jumps, config,                 -> (vals, vecs, info)
                                     |   scratch buffers)
                                     |
                                     v
                                  For each call f(v):
                                    rho = reshape(v, dim, dim)
                                    fill!(d_rho, 0)
                                    for jump in jumps:
                                      d_rho += dissipator(rho, jump, ...)
                                    d_rho += -i[B_total, rho]
                                    return copy(vec(d_rho))
```

### Component Boundaries

| Component | Responsibility | Location | Status |
|-----------|---------------|----------|--------|
| `krylov_spectral_gap()` | User-facing API: validates, precomputes, calls KrylovKit, packages results | New: `src/krylov.jl` | NEW |
| `build_lindbladian_action()` | Constructs the closure that applies L to vec(rho); captures precomputed data and scratch | New: `src/krylov.jl` | NEW |
| `_lindbladian_jump_action!()` | Applies ONE jump's Lindbladian contribution to rho (4 domain-dispatched methods) | New: `src/krylov.jl` | NEW |
| `_accumulate_dissipator!()` | Shared helper: rate * (A rho A^dag - 0.5 {A^dag A, rho}) | New: `src/krylov.jl` | NEW |
| `KrylovWorkspace` | Scratch buffers for the Lindbladian action (reused across Krylov iterations) | New: `src/krylov.jl` | NEW |
| `KrylovGapResult` | Result type: eigenvalues, gap, fixed point, convergence info | Add to `src/structs.jl` | NEW |
| `_precompute_data()` | Domain-specific precomputation (transition, energy labels, NUFFT prefactors) | Existing: `src/furnace_utensils.jl` | REUSE as-is |
| `_precompute_coherent_total_B()` | Precomputes total B operator for coherent term | Existing: `src/coherent.jl` | REUSE as-is |
| `oft!()` | Oscillatory Fourier Transform for EnergyDomain A_w computation | Existing: `src/ofts.jl` | REUSE as-is |
| `_prefactor_view()` | NUFFT prefactor matrix lookup for Time/TrotterDomain | Existing: `src/nufft.jl` | REUSE as-is |

### Where New Code Lives

**Single new file: `src/krylov.jl`** -- approximately 300-450 lines containing:
- `KrylovWorkspace` struct (~15 lines)
- `_accumulate_dissipator!()` and `_accumulate_dissipator_adjoint!()` shared helpers (~40 lines)
- `_lindbladian_jump_action!` for BohrDomain (~60 lines)
- `_lindbladian_jump_action!` for EnergyDomain (~50 lines)
- `_lindbladian_jump_action!` for TimeDomain/TrotterDomain (~50 lines)
- `build_lindbladian_action()` (~40 lines)
- `krylov_spectral_gap()` public API (~80 lines)
- Utility helpers (~30 lines)

**Modified files:**
- `src/structs.jl` -- add `KrylovGapResult` struct (~25 lines)
- `src/QuantumFurnace.jl` -- add `include("krylov.jl")` after `linearmaps_liouv.jl`, add exports
- `Project.toml` -- add KrylovKit dependency

**NO modifications to:** `furnace.jl`, `furnace_utensils.jl`, `coherent.jl`, `jump_workers.jl`, `diagnostics.jl`, `trajectories.jl`, `ofts.jl`, `nufft.jl`, or any other existing source file.

### Data Flow

```
krylov_spectral_gap(jumps, config, hamiltonian; trotter, n_eigenvalues, krylov_dim, tol)
  |
  +-- validate_config!(config)
  +-- ham_or_trott = (TrotterDomain ? trotter : hamiltonian)
  +-- precomputed_data = _precompute_data(config, ham_or_trott)          # REUSE existing
  +-- B_total = _precompute_coherent_total_B(jumps, ham_or_trott, ...)   # REUSE existing
  +-- ws = KrylovWorkspace(ComplexF64, dim)                              # NEW workspace
  |
  +-- lindbladian_action = build_lindbladian_action(
  |     jumps, ham_or_trott, config, precomputed_data, B_total, ws)
  |     |
  |     +-- Returns closure: function(v::AbstractVector{ComplexF64}) -> Vector{ComplexF64}
  |           |
  |           +-- rho = reshape(v, dim, dim)      # reshape view, no copy
  |           +-- fill!(ws.d_rho, 0)
  |           +-- for jump in jumps:
  |           |     _lindbladian_jump_action!(ws.d_rho, rho, jump, ham_or_trott,
  |           |                               config, precomputed_data, ws)
  |           +-- if B_total !== nothing:
  |           |     _coherent_action!(ws.d_rho, rho, B_total, ws)  # -i[B, rho]
  |           +-- return copy(vec(ws.d_rho))      # MUST be a fresh vector for KrylovKit
  |
  +-- x0 = randn(ComplexF64, dim^2)   # random starting vector
  +-- vals, vecs, info = KrylovKit.eigsolve(
  |     lindbladian_action, x0, n_eigenvalues, :SR;
  |     krylovdim=krylov_dim, tol=tol,
  |     issymmetric=false, ishermitian=false)
  |
  +-- perm = sortperm(abs.(real.(vals)))           # sort by proximity to zero
  +-- spectral_gap = abs(real(vals[perm[2]]))
  +-- fixed_point = reshape(vecs[perm[1]], dim, dim); normalize
  +-- gap_mode = reshape(vecs[perm[2]], dim, dim)
  +-- fixed_point_distance = trace_distance_h(Hermitian(fixed_point), gibbs)
  +-- Return KrylovGapResult(...)
```

## Patterns to Follow

### Pattern 1: Domain-Dispatched Lindbladian Action (Mirrors Existing `_jump_contribution!`)

**What:** Write `_lindbladian_jump_action!` as domain-dispatched methods that mirror the mathematical STRUCTURE of the existing `_jump_contribution!` for `AbstractLiouvConfig` in `jump_workers.jl`, but compute L(rho) on a density matrix instead of assembling the dim^2 x dim^2 Liouvillian.

**When:** Every Krylov iteration calls this for each jump operator.

**Why:** The existing `_jump_contribution!` for `AbstractLiouvConfig` builds the dim^2 x dim^2 matrix by calling `_vectorize_liouv_diss_and_add!` (which uses `_kron!` to place `A kron conj(A)` etc. into L). We need the ACTION of that matrix on a vector. Mathematically: instead of building L and computing L*v, we compute L(rho) directly where rho = reshape(v, dim, dim).

The mapping from vectorized Liouvillian terms to density matrix operations is:

| Vectorized term (in `_jump_contribution!`) | Direct action on rho |
|--------------------------------------------|---------------------|
| `_kron!(L, A, conj(A), rate)` adds rate*(A kron conj(A)) | `d_rho += rate * A * rho * A'` |
| `_kron!(L, A'A, I, -0.5*rate)` | `d_rho -= 0.5*rate * A'A * rho` |
| `_kron!(L, I, (A'A)^T, -0.5*rate)` | `d_rho -= 0.5*rate * rho * A'A` |
| `_kron!(L, B, I, -1im)` (coherent) | `d_rho -= 1im * B * rho` |
| `_kron!(L, I, B^T, +1im)` (coherent) | `d_rho += 1im * rho * B` |

**Example (EnergyDomain):**
```julia
function _lindbladian_jump_action!(
    d_rho::Matrix{<:Complex},
    rho::Matrix{<:Complex},
    jump::JumpOp,
    hamiltonian::HamHam,
    config::AbstractLiouvConfig{EnergyDomain},
    precomputed_data,
    ws::KrylovWorkspace,
)
    (; transition, gamma_norm_factor, energy_labels) = precomputed_data
    prefactor = config.w0 / (config.sigma * sqrt(2 * pi)) * gamma_norm_factor

    if jump.hermitian
        for w_raw in energy_labels
            w_raw > 1e-12 && continue
            w = abs(w_raw)
            oft!(ws.A_w, jump, w, hamiltonian, config.sigma)
            _accumulate_dissipator!(d_rho, rho, ws.A_w, prefactor * transition(w), ws)
            if w > 1e-12
                _accumulate_dissipator_adjoint!(d_rho, rho, ws.A_w,
                    prefactor * transition(-w), ws)
            end
        end
    else
        for w in energy_labels
            oft!(ws.A_w, jump, w, hamiltonian, config.sigma)
            _accumulate_dissipator!(d_rho, rho, ws.A_w, prefactor * transition(w), ws)
        end
    end
end
```

### Pattern 2: Shared Dissipator Accumulation Helpers

**What:** Two small helper functions that accumulate the standard Lindblad dissipator term.

**When:** Called from every domain's `_lindbladian_jump_action!`.

**Why:** All four domains share the same dissipator formula `D(A, rho) = A rho A^dag - 0.5{A^dag A, rho}`. Only the computation of A (the filtered/projected jump operator) differs per domain. Factoring out the dissipator accumulation eliminates code duplication across domains.

```julia
"""
Accumulate: d_rho += rate * (A * rho * A' - 0.5 * A'A * rho - 0.5 * rho * A'A)
"""
function _accumulate_dissipator!(
    d_rho::Matrix{<:Complex},
    rho::Matrix{<:Complex},
    A::AbstractMatrix{<:Complex},
    rate::Real,
    ws::KrylovWorkspace,
)
    # A * rho * A'
    mul!(ws.tmp1, A, rho)
    mul!(ws.tmp2, ws.tmp1, A')
    d_rho .+= rate .* ws.tmp2

    # -0.5 * A'A * rho - 0.5 * rho * A'A
    mul!(ws.LdagL, A', A)
    mul!(ws.tmp1, ws.LdagL, rho)
    d_rho .-= (0.5 * rate) .* ws.tmp1
    mul!(ws.tmp1, rho, ws.LdagL)
    d_rho .-= (0.5 * rate) .* ws.tmp1
    return nothing
end

"""
Accumulate dissipator with A' as the Lindblad operator (negative-frequency partner):
d_rho += rate * (A' * rho * A - 0.5 * AA' * rho - 0.5 * rho * AA')
"""
function _accumulate_dissipator_adjoint!(
    d_rho::Matrix{<:Complex},
    rho::Matrix{<:Complex},
    A::AbstractMatrix{<:Complex},
    rate::Real,
    ws::KrylovWorkspace,
)
    # A' * rho * A
    mul!(ws.tmp1, A', rho)
    mul!(ws.tmp2, ws.tmp1, A)
    d_rho .+= rate .* ws.tmp2

    # -0.5 * AA' * rho - 0.5 * rho * AA'
    mul!(ws.LdagL, A, A')
    mul!(ws.tmp1, ws.LdagL, rho)
    d_rho .-= (0.5 * rate) .* ws.tmp1
    mul!(ws.tmp1, rho, ws.LdagL)
    d_rho .-= (0.5 * rate) .* ws.tmp1
    return nothing
end
```

### Pattern 3: Closure-Based Linear Map for KrylovKit

**What:** KrylovKit accepts `f(x) -> y` where x and y are vectors. We build a closure that captures all precomputed data and scratch space, making the closure itself zero-allocation on the hot path (except the final `copy(vec(d_rho))`).

**When:** Constructing the linear map for `KrylovKit.eigsolve`.

**Why:** KrylovKit calls the linear map O(krylov_dim * n_restarts) times. Each call must be fast. Closures in Julia are efficient when the captured variables are type-stable.

**Critical detail:** KrylovKit internally stores the returned vectors to build the Krylov subspace. The closure MUST return a new vector (or a copy) each time, not a view into the same buffer. Use `copy(vec(ws.d_rho))` to produce a fresh allocation. This is the one unavoidable allocation per Krylov iteration.

```julia
function build_lindbladian_action(
    jumps::Vector{JumpOp},
    ham_or_trott::Union{HamHam, TrottTrott},
    config::AbstractLiouvConfig,
    precomputed_data,
    B_total::Union{Nothing, Matrix{<:Complex}},
    ws::KrylovWorkspace,
)
    dim = ham_or_trott isa HamHam ? size(ham_or_trott.data, 1) :
                                     size(ham_or_trott.eigvecs, 1)

    function action(v::AbstractVector{<:Complex})
        rho = reshape(v, dim, dim)
        fill!(ws.d_rho, 0)

        for jump in jumps
            _lindbladian_jump_action!(ws.d_rho, rho, jump,
                ham_or_trott, config, precomputed_data, ws)
        end

        if B_total !== nothing
            # Coherent: -i[B, rho] = -i*B*rho + i*rho*B
            mul!(ws.tmp1, B_total, rho)
            ws.d_rho .-= 1im .* ws.tmp1
            mul!(ws.tmp1, rho, B_total)
            ws.d_rho .+= 1im .* ws.tmp1
        end

        return copy(vec(ws.d_rho))
    end

    return action
end
```

### Pattern 4: Reuse Existing Precomputation Verbatim

**What:** Call the exact same `_precompute_data()` and `_precompute_coherent_total_B()` functions that `construct_lindbladian()` and `run_thermalization()` use.

**When:** Setting up the Krylov solve.

**Why:** The precomputed data (transition functions, gamma_norm_factor, energy labels, NUFFT prefactors, coherent B) is domain-specific and already battle-tested across 26 phases of development. Duplicating or reimplementing this logic would introduce bugs and diverge from the validated Liouvillian construction.

```julia
# In krylov_spectral_gap():
precomputed_data = _precompute_data(config, ham_or_trott)   # exact same call as furnace.jl
B_total = _precompute_coherent_total_B(jumps, ham_or_trott, config, precomputed_data)
```

### Pattern 5: Config Acceptance -- LiouvConfig Only

**What:** `krylov_spectral_gap()` accepts `AbstractLiouvConfig` (not `AbstractThermalizeConfig`).

**When:** API design.

**Why:** The Krylov solver computes the Lindbladian's eigenvalues, which are independent of the thermalization time step delta and mixing time. These are parameters of the CPTP channel approximation, not the generator. Using `AbstractLiouvConfig` makes this explicit and avoids confusion.

The existing `run_lindbladian()` already follows this pattern. A user who has a `ThermalizeConfig` can construct the corresponding `LiouvConfig` by dropping `mixing_time` and `delta`. A helper function for this conversion is recommended (see Pattern in existing `_thermalize_to_liouv_config` pattern from diagnostics research).

## Anti-Patterns to Avoid

### Anti-Pattern 1: Forming the Full Liouvillian Then Using KrylovKit

**What:** Building the full dim^2 x dim^2 dense Liouvillian matrix and then passing it to KrylovKit.

**Why bad:** At 12 qubits, the Liouvillian is 16M x 16M = 2 petabytes. The existing `construct_lindbladian()` works for n<=6 (dim^2 = 4096) but is fundamentally unscalable. Using KrylovKit on a dense matrix provides no advantage over the existing `eigs()` approach.

**Instead:** Use the matrix-free action described above.

### Anti-Pattern 2: Wrapping the DM Channel (E_delta) Instead of the Lindbladian (L)

**What:** Using the existing `_jump_contribution!(evolving_dm, ...)` for `AbstractThermalizeConfig` as the linear map.

**Why bad:**
1. Channel eigenvalues are mu_i = 1 + delta * lambda_i + O(delta^2). For delta=0.01, all eigenvalues cluster near 1.0. KrylovKit's Arnoldi must resolve |mu_1 - mu_2| ~ delta * gap, which is tiny.
2. The Kraus channel includes nonlinear PSD clamping (K0, U_residual, `_finalize_kraus_step!`) that makes it NOT a true linear operator. The eigenvalues of the approximate channel differ from the Lindbladian's.
3. The delta-dependence means the "gap" depends on the time step, which is a nuisance parameter.

**Instead:** Compute L(rho) directly from the dissipator formula.

### Anti-Pattern 3: Allocating Inside the Lindbladian Action Closure

**What:** Creating new matrices/vectors inside the closure that gets called per Krylov iteration.

**Why bad:** At 12 qubits, a single dim x dim ComplexF64 allocation is 256 MB and triggers GC pressure. KrylovKit calls the action hundreds of times.

**Instead:** Pre-allocate all scratch in `KrylovWorkspace` and reuse. The only allocation per call is `copy(vec(ws.d_rho))` (128 MB at 12 qubits) which is unavoidable since KrylovKit stores the result.

### Anti-Pattern 4: Modifying the Existing `_jump_contribution!` Methods

**What:** Adding Lindbladian-action logic to the existing `_jump_contribution!` methods in `jump_workers.jl`.

**Why bad:** The existing methods serve two distinct purposes (Liouvillian matrix assembly for `AbstractLiouvConfig`, DM channel application for `AbstractThermalizeConfig`). Adding a third mode (density-matrix-level Lindbladian action) would require a third dispatch dimension or runtime branching, making the code harder to maintain.

**Instead:** Write new `_lindbladian_jump_action!` methods in `src/krylov.jl`. They share the same mathematical structure but have clean, purpose-built signatures without the vectorization/Kronecker machinery of the Liouvillian methods or the Kraus/K0/residual machinery of the channel methods.

## Domain-Specific Architecture Details

### How Each Domain Affects the Krylov Computation

All four domains share the same Lindbladian structure:
```
L(rho) = sum_jumps sum_w [ rate(w) * D(A_w, rho) ] - i[B_total, rho]
where D(A, rho) = A rho A' - 0.5 {A'A, rho}
```

The domains differ in HOW A_w is computed and what the prefactor/rate is.

| Domain | A_w computation | Basis for rho | B computation | Precomputed data |
|--------|----------------|---------------|---------------|-----------------|
| **Bohr** | Sparse Bohr-bucket projection: `alpha(bohr_freqs, nu) .* jump.in_eigenbasis` | H eigenbasis | `B_bohr()` | `alpha` fn, `bohr_dict`, `gamma_norm_factor` |
| **Energy** | Gaussian filter: `oft!(A_w, jump, w, hamiltonian, sigma)` | H eigenbasis | `B_bohr()` | `transition` fn, `energy_labels`, `gamma_norm_factor` |
| **Time** | NUFFT prefactors: `jump.in_eigenbasis .* nufft_prefactor_view(w)` | H eigenbasis | `B_time()` | `transition`, `energy_labels`, `oft_nufft_prefactors`, `b_minus/b_plus` |
| **Trotter** | NUFFT prefactors: `jump.in_eigenbasis .* nufft_prefactor_view(w)` | Trotter eigenbasis | `B_trotter()` | Same as Time but using `trotter.eigvecs/eigvals_t0` |

**Basis handling:** For Bohr/Energy/Time domains, all density matrices live in the Hamiltonian eigenbasis. For TrotterDomain, they live in the Trotter eigenbasis. This is handled automatically because `jump.in_eigenbasis` is already constructed in the correct basis at JumpOp creation time. The Krylov vectors (length dim^2) represent vectorized density matrices in whichever basis the domain uses. No explicit basis transformation is needed in the Krylov action.

### BohrDomain: Distinct Loop Structure

BohrDomain does NOT iterate over an energy grid. Instead it iterates over Bohr frequency buckets from `hamiltonian.bohr_dict`. For each bucket nu2, it constructs:
- `alpha_A_{nu1}[i,j] = alpha(bohr_freqs[i,j], nu2) * jump.in_eigenbasis[i,j]` (Kossakowski-weighted jump, dense matrix)
- `A_{nu2}` is a SPARSE operator: only the entries (i,j) where `bohr_freqs[i,j] == nu2` are nonzero

The dissipator for each (nu1, nu2) pair is: `D(alpha_A_{nu1}, A_{nu2}^dag, rho)`. This is a generalized dissipator where the two operators differ: `alpha_A_{nu1} * rho * A_{nu2} - 0.5 * A_{nu2} * alpha_A_{nu1} * rho - 0.5 * rho * A_{nu2} * alpha_A_{nu1}`.

This matches the existing two-operator form `_vectorize_liouv_diss_and_add!(L, jump_1, jump_2_dag, scalar, ws)` used in the Bohr Liouvillian construction.

The Krylov action for BohrDomain needs a generalized dissipator helper that takes two different operators:

```julia
function _accumulate_dissipator_two_ops!(
    d_rho::Matrix{<:Complex},
    rho::Matrix{<:Complex},
    A1::AbstractMatrix{<:Complex},   # alpha_A_{nu1}
    A2_dag::AbstractMatrix{<:Complex}, # A_{nu2}^dag (sparse)
    rate::Real,
    ws::KrylovWorkspace,
)
    # rate * (A1 * rho * A2_dag' - 0.5 * A2_dag' * A1 * rho - 0.5 * rho * A2_dag' * A1)
    # Note: A2_dag' = A2 (adjoint of the adjoint)
    mul!(ws.tmp1, A1, rho)
    mul!(ws.tmp2, ws.tmp1, A2_dag')     # A1 * rho * A2
    d_rho .+= rate .* ws.tmp2

    mul!(ws.tmp1, A2_dag', A1)           # A2 * A1
    mul!(ws.tmp2, ws.tmp1, rho)
    d_rho .-= (0.5 * rate) .* ws.tmp2
    mul!(ws.tmp2, rho, ws.tmp1)
    d_rho .-= (0.5 * rate) .* ws.tmp2
end
```

Wait -- examining the existing BohrDomain Liouvillian code more carefully. The `_jump_contribution!` for BohrDomain in `jump_workers.jl` lines 11-45 does:

```julia
alpha_A_nu1 = alpha(bohr_freqs, nu2) .* jump.in_eigenbasis  # full dim x dim
A_nu_2_dag = sparse(rows, cols, conj(vals), dim, dim)         # sparse projection
_vectorize_liouv_diss_and_add!(L, alpha_A_nu1, A_nu_2_dag, gamma_norm_factor, ws)
```

So the dissipator is: `gamma_norm * (alpha_A * rho * A_nu2 - 0.5 * A_nu2 * alpha_A * rho - 0.5 * rho * A_nu2 * alpha_A)`. The density-matrix action version needs to construct `A_nu2_dag` as a sparse matrix (or apply it element-by-element). Given that BohrDomain is the "highest approximation" domain and typically used only for small systems, constructing the sparse matrix is fine.

### Energy/Time/Trotter: Shared Loop Structure

These three domains share the same loop: iterate over energy labels, compute A_w, accumulate the standard single-operator dissipator. They differ only in how A_w is formed:
- **Energy:** `oft!(A_w, jump, w, hamiltonian, sigma)` -- in-place Gaussian filter
- **Time/Trotter:** `A_w[i,j] = jump.in_eigenbasis[i,j] * nufft_prefactors[i,j,k(w)]` -- elementwise multiply with precomputed NUFFT

Time and Trotter share nearly identical Krylov action code. The only difference is implicit: `ham_or_trott` is a `HamHam` for TimeDomain and a `TrottTrott` for TrotterDomain, and `jump.in_eigenbasis` is already in the correct basis.

## Result Type Design

### New: KrylovGapResult

```julia
"""
    KrylovGapResult

Result of Krylov-based Lindbladian spectral gap estimation.

# Fields
- `spectral_gap::Float64`: |Re(lambda_2)|, the Lindbladian spectral gap.
- `eigenvalues::Vector{ComplexF64}`: Leading eigenvalues sorted by |Re(lambda)|.
- `fixed_point::Matrix{ComplexF64}`: Normalized density matrix from lambda_1 eigenvector.
- `fixed_point_distance::Float64`: Trace distance from fixed point to Gibbs state.
- `gap_mode::Matrix{ComplexF64}`: Density-matrix-shaped eigenvector for lambda_2.
- `converged::Int`: Number of converged eigenvalues (from KrylovKit info).
- `residual_norms::Vector{Float64}`: Residual norms for each eigenvalue.
- `n_matvecs::Int`: Number of Lindbladian action evaluations (from info.numops).
- `krylov_dim::Int`: Krylov subspace dimension used.
- `tol::Float64`: Convergence tolerance used.
"""
struct KrylovGapResult
    spectral_gap::Float64
    eigenvalues::Vector{ComplexF64}
    fixed_point::Matrix{ComplexF64}
    fixed_point_distance::Float64
    gap_mode::Matrix{ComplexF64}
    converged::Int
    residual_norms::Vector{Float64}
    n_matvecs::Int
    krylov_dim::Int
    tol::Float64
end
```

**Relationship to existing types:**
- `LindbladianResult` stores the full Liouvillian matrix + 2 eigenvalues. `KrylovGapResult` stores NO matrix but more eigenvalues.
- `ExactDiagnosticsResult` provides 6 diagnostic outputs from dense eigen. `KrylovGapResult` provides spectral data from matrix-free Krylov.
- `SpectralGapResult` is trajectory-based (stochastic). `KrylovGapResult` is Krylov-based (deterministic, exact up to tolerance).

## Workspace Design

```julia
struct KrylovWorkspace{T<:Complex}
    d_rho::Matrix{T}     # output accumulator: L(rho)
    A_w::Matrix{T}       # filtered jump operator at frequency w
    LdagL::Matrix{T}     # L^dag L or similar product scratch
    tmp1::Matrix{T}      # general scratch 1
    tmp2::Matrix{T}      # general scratch 2
end

function KrylovWorkspace(::Type{T}, dim::Int) where {T<:Complex}
    Zm() = zeros(T, dim, dim)
    return KrylovWorkspace{T}(Zm(), Zm(), Zm(), Zm(), Zm())
end
```

**Memory budget at key qubit counts:**

| n_qubits | dim | dim^2 | KrylovWorkspace (5 matrices) | Krylov subspace (kd=20) | Total |
|----------|-----|-------|------------------------------|------------------------|-------|
| 4 | 16 | 256 | 20 KB | 40 KB | <1 MB |
| 6 | 64 | 4,096 | 320 KB | 640 KB | ~1 MB |
| 8 | 256 | 65,536 | 5 MB | 10 MB | ~15 MB |
| 10 | 1,024 | 1,048,576 | 80 MB | 160 MB | ~250 MB |
| 12 | 4,096 | 16,777,216 | 1.3 GB | 2.6 GB | ~4 GB |

At n=12 with krylov_dim=20, the total memory is approximately 4 GB, feasible on a 16+ GB workstation.

## Validation Strategy

### Cross-Validation Against Dense Eigen (n=4,6)

For small systems where the dense Liouvillian is feasible, validate that:
1. `krylov_result.spectral_gap` matches `exact_diagnostics.eigen.spectral_gap` to within `tol`.
2. `krylov_result.fixed_point_distance` is near zero.
3. All Krylov eigenvalues match the leading dense eigenvalues.

This leverages the existing `construct_lindbladian()` + `run_exact_diagnostics()` infrastructure.

```julia
# Validation test structure:
for domain in [BohrDomain(), EnergyDomain(), TimeDomain(), TrotterDomain()]
    L_dense = construct_lindbladian(jumps, liouv_config, hamiltonian; trotter=trotter)
    exact = run_exact_diagnostics(L_dense, hamiltonian, gibbs; ...)

    krylov = krylov_spectral_gap(jumps, liouv_config, hamiltonian; trotter=trotter)

    @test isapprox(krylov.spectral_gap, exact.eigen.spectral_gap; rtol=1e-6)
    @test krylov.fixed_point_distance < 1e-6
end
```

### Linearity Verification

Verify that the Lindbladian action is actually linear by checking `f(a*x + b*y) == a*f(x) + b*f(y)` for random x, y and scalars a, b. This catches bugs where scratch buffers are incorrectly reused or where the action has inadvertent state.

## Scalability Considerations

| Concern | n=4 (dim=16) | n=6 (dim=64) | n=8 (dim=256) | n=10 (dim=1024) | n=12 (dim=4096) |
|---------|-------------|-------------|--------------|----------------|----------------|
| One L(rho) eval (3 jumps, EnergyDomain) | <1ms | ~5ms | ~200ms | ~5s | ~2min (est.) |
| Full Krylov solve (kd=20, 2 restarts) | <1s | ~1s | ~10s | ~5min | ~2hr (est.) |
| Full dense Liouvillian | instant | ~2s | ~30s (2GB) | impossible | impossible |
| Dense eigen of L | instant | ~1s | ~30min | impossible | impossible |
| Krylov advantage factor | none | ~2x | ~100x | infinity | infinity |

**Key bottleneck at large n:** Each L(rho) involves O(n_jumps * n_energy_labels) matrix multiplications of size dim x dim. The cost per mul! is O(dim^3) via BLAS. With truncated energy labels (~200-400) and 36 jumps (3 Paulis x 12 qubits), each L(rho) requires ~10,000 mul! calls at dim=4096.

**Estimated times with BLAS threading (8 cores):**
- Single mul!(dim=4096): ~2ms with multi-threaded BLAS
- Per L(rho): ~20 seconds
- Full Krylov (40 iterations): ~13 minutes

This is the breakthrough: matrix-free Krylov at n=10-12 takes minutes to hours, while forming the full Liouvillian is impossible (terabytes of memory).

## Build Order (Incremental Testing)

### Phase 1: Core Infrastructure and EnergyDomain (Test at n=4)

1. Add KrylovKit to `Project.toml`
2. `KrylovWorkspace` struct
3. `_accumulate_dissipator!()` and `_accumulate_dissipator_adjoint!()` helpers
4. `_lindbladian_jump_action!` for EnergyDomain
5. `build_lindbladian_action()` (EnergyDomain path)
6. `krylov_spectral_gap()` API
7. `KrylovGapResult` in `src/structs.jl`
8. Test: Krylov gap matches dense eigen for n=4 EnergyDomain (with and without coherent)

**Why EnergyDomain first:** Simplest A_w computation (direct `oft!`), no NUFFT dependency, no Bohr-bucket iteration.

### Phase 2: Remaining Domains (Test at n=4,6)

9. `_lindbladian_jump_action!` for TimeDomain/TrotterDomain (NUFFT prefactors)
10. `_lindbladian_jump_action!` for BohrDomain (Bohr buckets, generalized dissipator)
11. Coherent term action integration
12. Cross-domain validation at n=4 and n=6 for all 4 domains
13. TrotterDomain-specific test (verify basis handling)

### Phase 3: Scaling and Polish (n=8+)

14. Test at n=8 (first size with genuine Krylov advantage)
15. Performance profiling and optimization
16. BLAS threading configuration
17. Krylov parameter guidance (krylov_dim, tol defaults)
18. Integration with existing diagnostics pipeline

## Dependency: KrylovKit.jl

KrylovKit.jl is NOT currently in `Project.toml`. It must be added as a dependency.

**Why KrylovKit over Arpack (already a dependency):**
1. KrylovKit accepts `f(x) -> y` directly; Arpack requires wrapping in `LinearMap`.
2. KrylovKit's Krylov-Schur algorithm with thick restarts is more modern than Arpack's IRAM.
3. KrylovKit returns richer convergence info.
4. KrylovKit works with arbitrary Julia vector types.

**Historical context:** The existing `linearmaps_liouv.jl` contains a commented-out prototype that attempted matrix-free Liouvillian maps via Arpack. The comment says "Sadly, slow. Misery." This was because the prototype recomputed OFT per Krylov iteration WITHOUT the NUFFT prefactor cache (which was added later in the project). With current precomputed NUFFT prefactors, the matrix-free approach should be dramatically faster.

**The LinearMaps dependency is already present** in `Project.toml` but is unused (the code in `linearmaps_liouv.jl` is fully commented out). KrylovKit does not need `LinearMaps` since it accepts plain functions.

## Sources

- [KrylovKit.jl documentation](https://jutho.github.io/KrylovKit.jl/stable/)
- [KrylovKit.jl eigsolve API](https://jutho.github.io/KrylovKit.jl/stable/man/eig/)
- [KrylovKit.jl GitHub repository](https://github.com/Jutho/KrylovKit.jl)
- Existing codebase: `src/jump_workers.jl` (Liouvillian assembly + channel application patterns)
- Existing codebase: `src/furnace.jl` (run_lindbladian entry point, precomputation flow)
- Existing codebase: `src/linearmaps_liouv.jl` (abandoned matrix-free prototype with historical context)
- Existing codebase: `src/diagnostics.jl` (dense exact diagnostics for validation reference)
- Existing codebase: `src/furnace_utensils.jl` (_precompute_data domain dispatch)
- Existing codebase: `src/qi_tools.jl` (_vectorize_liouv_diss_and_add! for understanding L(rho) semantics)
