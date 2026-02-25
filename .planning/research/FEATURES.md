# Feature Landscape: QuantumFurnace.jl Codebase Restructure

**Domain:** Refactoring a Julia quantum Gibbs sampling simulation package
**Researched:** 2026-02-25
**Confidence:** HIGH (based on exhaustive reading of all 28 source files, all structs, all duplication sites)

---

## Scope

This document defines the **concrete architectural features** for restructuring QuantumFurnace.jl. Every recommendation below is grounded in exact field names, function signatures, and line counts from the current codebase. The restructure touches:

1. **Config type hierarchy** -- 4 structs to 1 parametric struct
2. **Workspace consolidation** -- 5 workspace types to 2-3
3. **Prefactor/OFT/channel deduplication** -- 16+ copy-pasted code blocks to single-source functions
4. **Result struct cleanup** -- 4+ result types with uniform save capability
5. **File renaming** -- PRE/MID/POST logical grouping in flat `src/`
6. **Test deduplication** -- factory patterns replacing ad-hoc setup code

---

## Current State: What Exists

### Config Structs (4 structs, 90% identical fields)

| Struct | Lines | Unique Fields | Purpose |
|--------|-------|---------------|---------|
| `LiouvConfig{D,T}` | 73-88 | none (base) | Lindbladian spectral analysis |
| `LiouvConfigGNS{D,T}` | 99-114 | `with_coherent=false` enforced | GNS variant of above |
| `ThermalizeConfig{D,T}` | 143-162 | `mixing_time`, `delta` | DM thermalization simulation |
| `ThermalizeConfigGNS{D,T}` | 174-192 | `with_coherent=false` + `mixing_time`, `delta` | GNS variant of above |

All 4 share these 14 fields verbatim: `num_qubits`, `with_coherent`, `with_linear_combination`, `domain`, `beta`, `sigma`, `gaussian_parameters`, `a`, `b`, `num_energy_bits`, `t0`, `w0`, `eta`, `num_trotter_steps_per_t0`.

The abstract hierarchy is: `AbstractConfig{D,T}` > `AbstractLiouvConfig{D,T}` / `AbstractThermalizeConfig{D,T}`.

### Workspace Types (5 types)

| Workspace | Location | Fields | Used By |
|-----------|----------|--------|---------|
| `LindbladianWorkspace{T}` | structs.jl:21 | Id, jump_tmp, jump_conj, jump_dag_jump, jump2_jump1 (5 matrices) | `construct_lindbladian` vectorized Liouvillian |
| `KrausScratch{T}` | kraus.jl:1 | jump_oft, LdagL, R, rho_jump, K0, tmp1, tmp2, rho_next (8 matrices) | `_jump_contribution!` for DM thermalization |
| `KrylovWorkspace{T,PD}` | krylov_workspace.jl:43 | precomputed_data, B_total, jumps, jump_eigenbases, jump_hermitian, jump_oft, tmp1, tmp2, LdagL, rho_out + 9 channel/G fields (23 fields) | `apply_lindbladian!`, `apply_delta_channel!` |
| `TrajectoryWorkspace{T}` | trajectories.jl:4 | jump_oft, psi_tmp, Rpsi, rho_acc (3 matrices + 1 vector) | `step_along_trajectory!` |
| `OFTCaches{T}` | structs.jl:347 | prefactors, U, temp_op (deprecated) | Legacy OFT functions only |

### Result Types (5+ types)

| Result | Location | Fields | Serializable |
|--------|----------|--------|-------------|
| `LindbladianResult{T}` | structs.jl:259 | liouvillian, fixed_point, gap_mode, spectral_gap | No (manual) |
| `DMSimulationResult{T}` | structs.jl:239 | final_dm, trace_distances, time_steps | No (manual) |
| `KrylovGapResult{T}` | krylov_eigsolve.jl:40 | eigenvalues, spectral_gap, fixed_point, gap_mode, converged, matvec_count, num_restarts, normres, channel_eigenvalues, delta_used | No |
| `TrajectoryResult{T}` | trajectories.jl:23 | rho_mean, n_trajectories, seed, times, measurements_mean, convergence | Via ExperimentResult |
| `ObservableTrajectoryResult{T}` | trajectories.jl:37 | times, measurements_mean, n_trajectories, seed, rho_mean | No |
| `ExperimentResult{C,T}` | results.jl:18 | config, trajectory_result, hamiltonian_params, metadata | Yes (BSON) |
| `SpectralGapResult` | gap_estimation.jl:37 | gap, gap_ci, gap_se, best_observable, per_observable, ... | No |

### Duplication Sites (measured)

| Pattern | Occurrences | Files |
|---------|-------------|-------|
| EnergyDomain prefactor formula | 8 | jump_workers, krylov_workspace, krylov_matvec, krylov_eigsolve, trajectories |
| Time/TrotterDomain prefactor formula | 8 | same 5 files |
| Hermitian half-grid branch (`w_raw > 1e-12 && continue`) | 16 | 5 files |
| `hermitianize!(scratch.R)` after R accumulation | 5 | jump_workers, trajectories |
| R accumulation loop (Energy) | 3 | jump_workers (2x), krylov_workspace |
| R accumulation loop (Time/Trotter) | 3 | jump_workers (2x), krylov_workspace |
| `_thermalize_to_liouv_config` field-by-field copy | 2 | krylov_workspace.jl |
| Sandwich function variants | 6 | krylov_matvec.jl (4 are just transposes of 2) |
| `_precompute_coherent_*` near-identical domain branches | 3 functions | coherent.jl |
| `_select_b_plus_calculator` (KMS-only, no GNS guard) | 1 | furnace_utensils.jl |

---

## Table Stakes

Features that MUST be in this restructure. Without these, the refactor is incomplete and the codebase remains duplicated.

### TS-01: Unified Config Struct -- `SimConfig{S,D,DB,T}`

| Aspect | Detail |
|--------|--------|
| **Why Expected** | 4 structs with 14 identical fields is the #1 maintenance burden. Adding DLL variant would require 2 more structs (6 total). |
| **Complexity** | HIGH -- touches every function signature, all dispatch, serialization, tests |
| **Notes** | This is the keystone feature; everything else depends on it |

**Concrete design:**

```julia
# --- Simulation type tags ---
abstract type AbstractSimType end
struct Lindblad <: AbstractSimType end    # was "Liouv"
struct Thermalize <: AbstractSimType end
struct KrylovSpectrum <: AbstractSimType end
struct Trajectory <: AbstractSimType end

# --- Detailed balance type tags ---
abstract type AbstractDBType end
struct KMS <: AbstractDBType end
struct GNS <: AbstractDBType end
# future: struct DLL <: AbstractDBType end

# --- Single config struct ---
@kwdef struct SimConfig{S<:AbstractSimType, D<:AbstractDomain, DB<:AbstractDBType, T<:AbstractFloat}
    # === Physics (shared by ALL configs) ===
    num_qubits::Int
    beta::T
    sigma::T
    domain::D

    # === Linear combination parameters ===
    with_linear_combination::Bool
    gaussian_parameters::Tuple{Union{T,Nothing}, Union{T,Nothing}} = (nothing, nothing)
    a::Union{T, Nothing} = nothing
    b::Union{T, Nothing} = nothing

    # === Discretization (required for Energy/Time/Trotter, nothing for Bohr) ===
    num_energy_bits::Union{Int, Nothing} = nothing
    t0::Union{T, Nothing} = nothing
    w0::Union{T, Nothing} = nothing
    eta::Union{T, Nothing} = nothing
    num_trotter_steps_per_t0::Union{Int, Nothing} = nothing

    # === Thermalization-specific (only for S <: Union{Thermalize, Trajectory}) ===
    mixing_time::Union{T, Nothing} = nothing
    delta::Union{T, Nothing} = nothing
end
```

**Key decisions:**
- `with_coherent` is **removed as a field**. It is derived: `has_coherent(::SimConfig{S,D,KMS}) = true`, `has_coherent(::SimConfig{S,D,GNS}) = false`. This eliminates the "GNS configs must have with_coherent=false" runtime guard entirely.
- `mixing_time` and `delta` become optional fields (Union with Nothing). They are required for `Thermalize`/`Trajectory` simulation types and validated at construction via an inner constructor or `validate_config!`.
- Four type parameters because: `S` selects simulation runner, `D` selects domain dispatch, `DB` selects detailed balance physics (alpha, transition, coherent), `T` is numeric precision.

**Dispatch implications:**
- `run_lindblad(config::SimConfig{Lindblad})` replaces `run_lindbladian(config::AbstractLiouvConfig)`
- `run_thermalize(config::SimConfig{Thermalize})` replaces `run_thermalization(config::AbstractThermalizeConfig)`
- `pick_transition(config::SimConfig{S,D,KMS})` replaces `pick_transition(config::LiouvConfig)` and `pick_transition(config::ThermalizeConfig)`
- `pick_transition(config::SimConfig{S,D,GNS})` replaces `pick_transition(config::LiouvConfigGNS)` and `pick_transition(config::ThermalizeConfigGNS)`
- `_precompute_data(config::SimConfig{S,EnergyDomain})` unifies all energy domain precomputation regardless of sim type
- `has_coherent(::SimConfig{S,D,KMS}) where {S,D}` = true, eliminates the `with_coherent` flag

**DLL extensibility:** Adding DLL requires only `struct DLL <: AbstractDBType end` plus new `pick_transition(::SimConfig{S,D,DLL})` and `_pick_alpha(::SimConfig{S,D,DLL})` methods. Zero changes to existing code.

**Backward compatibility for serialization:** The `_config_to_dict` / `_reconstruct_config` in results.jl already uses string tags `"KMS"/"GNS"` and `"liouv"/"thermalize"`. Extend with `"lindblad"/"krylov_spectrum"/"trajectory"` for `S`, add `"DLL"` for `DB`.

---

### TS-02: Prefactor Deduplication -- `domain_prefactor(config, gamma_norm_factor)`

| Aspect | Detail |
|--------|--------|
| **Why Expected** | 16 occurrences of 2 formulas across 5 files. Bug in one = silent divergence. |
| **Complexity** | LOW -- pure extraction, no logic changes |
| **Notes** | Must be done BEFORE any other inner-loop refactoring |

**Concrete design:**

```julia
# In a new file: src/prefactors.jl (PRE layer)

"""
    domain_prefactor(config::SimConfig{S,EnergyDomain}, gamma_norm_factor) -> Float64

Energy domain: w0 / (sigma * sqrt(2pi)) * gamma_norm_factor
"""
@inline function domain_prefactor(
    config::SimConfig{S,EnergyDomain}, gamma_norm_factor::Real
) where {S}
    return config.w0 / (config.sigma * sqrt(2 * pi)) * gamma_norm_factor
end

"""
    domain_prefactor(config::SimConfig{S,D}, gamma_norm_factor) -> Float64

Time/Trotter domain: w0 * t0^2 * sigma * sqrt(2/pi) / (2pi) * gamma_norm_factor
"""
@inline function domain_prefactor(
    config::SimConfig{S,D}, gamma_norm_factor::Real
) where {S, D<:Union{TimeDomain, TrotterDomain}}
    return config.w0 * config.t0^2 * (config.sigma * sqrt(2 / pi)) / (2 * pi) * gamma_norm_factor
end
```

Every site that currently computes the prefactor inline becomes a call to `domain_prefactor(config, gamma_norm_factor)`. For trajectory framework hot paths where config access is avoided, the value is precomputed once at framework construction and stored as `scaled_prefactor::Float64` (already the pattern in `TrajectoryFramework`).

---

### TS-03: Hermitian Half-Grid Iterator -- `foreach_energy_pair`

| Aspect | Detail |
|--------|--------|
| **Why Expected** | 16 copies of the hermitian branching pattern with half-grid iteration |
| **Complexity** | MEDIUM -- needs careful closure design for zero allocations |
| **Notes** | The biggest single source of code duplication |

**Concrete design:**

```julia
"""
    foreach_energy_pair(f_pos, f_neg, energy_labels, is_hermitian)

Iterate energy labels with hermitian half-grid optimization.
- If is_hermitian: iterate w <= 0 only, call f_pos(w) and f_neg(w) for w > 0.
- If !is_hermitian: iterate all w, call f_pos(w) only.

Replaces 16 copies of the `if jump.hermitian ... for w_raw ... w_raw > 1e-12 && continue` pattern.
"""
@inline function foreach_energy_pair(
    f_pos, f_neg, energy_labels::AbstractVector{<:Real}, is_hermitian::Bool
)
    if is_hermitian
        @inbounds for w_raw in energy_labels
            w_raw > 1e-12 && continue
            w = abs(w_raw)
            f_pos(w)
            if w > 1e-12
                f_neg(w)
            end
        end
    else
        @inbounds for w in energy_labels
            f_pos(w)
        end
    end
    return nothing
end
```

This is a **do-block pattern** -- callers pass closures that capture their local variables. The closures must be typed (not `Any`) for zero-allocation. Julia's compiler specializes on closure types, so this is allocation-free when the closures are non-boxing.

**Example migration (from krylov_matvec.jl EnergyDomain `apply_lindbladian!`):**

```julia
# BEFORE (17 lines per jump):
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

# AFTER (8 lines per jump):
for (k, eigenbasis) in enumerate(ws.jump_eigenbases)
    foreach_energy_pair(energy_labels, ws.jump_hermitian[k]) do w
        compute_jump_oft!(ws.jump_oft, eigenbasis, oft_source, w)
        scalar_w = prefactor * transition(w)
        _accumulate_sandwich!(ws.rho_out, ws.jump_oft, rho, scalar_w, ws)
    end do w
        scalar_neg = prefactor * transition(-w)
        _accumulate_sandwich_adj_L!(ws.rho_out, ws.jump_oft, rho, scalar_neg, ws)
    end
end
```

Actually, Julia's do-block syntax only supports one closure. A better pattern is:

```julia
for (k, eigenbasis) in enumerate(ws.jump_eigenbases)
    foreach_energy_pair(energy_labels, ws.jump_hermitian[k],
        w -> begin  # positive frequency action
            compute_jump_oft!(ws.jump_oft, eigenbasis, oft_source, w)
            _accumulate_sandwich!(ws.rho_out, ws.jump_oft, rho, prefactor * transition(w), ws)
        end,
        w -> begin  # negative frequency action (only for hermitian mirrors)
            _accumulate_sandwich_adj_L!(ws.rho_out, ws.jump_oft, rho, prefactor * transition(-w), ws)
        end
    )
end
```

The key constraint is that `ws.jump_oft` must already be computed by `f_pos` before `f_neg` reads it. The current code relies on this sequential ordering, and `foreach_energy_pair` preserves it by calling `f_pos(w)` then `f_neg(w)` in that order.

---

### TS-04: Unified OFT Computation -- `compute_jump_oft!`

| Aspect | Detail |
|--------|--------|
| **Why Expected** | `oft!`, `_krylov_oft!`, and inline NUFFT prefactor multiply are 3 ways to do the same thing |
| **Complexity** | LOW -- merge to single dispatched function |
| **Notes** | `_krylov_oft!` is slightly more efficient than `oft!` (uses `inv_4sigma2` precomputed) |

**Concrete design:**

```julia
"""
    compute_jump_oft!(out, eigenbasis, bohr_freqs, w, inv_4sigma2)

EnergyDomain OFT: Gaussian filter in eigenbasis.
    out[i,j] = eigenbasis[i,j] * exp(-(w - bohr_freqs[i,j])^2 * inv_4sigma2)
"""
@inline function compute_jump_oft!(
    out::Matrix{T}, eigenbasis::Matrix{T},
    bohr_freqs::Matrix{<:Real}, w::Real, inv_4sigma2::Real
) where {T<:Complex}
    @. out = eigenbasis * exp(-(w - bohr_freqs)^2 * inv_4sigma2)
    return nothing
end

"""
    compute_jump_oft!(out, eigenbasis, nufft_prefactors, w)

Time/TrotterDomain OFT: NUFFT prefactor multiply.
    out[i,j] = eigenbasis[i,j] * prefactor_view(nufft_prefactors, w)[i,j]
"""
@inline function compute_jump_oft!(
    out::Matrix{T}, eigenbasis::Matrix{T},
    nufft_prefactors::NUFFTPrefactors, w::Real
) where {T<:Complex}
    pf = _prefactor_view(nufft_prefactors, w)
    @. out = eigenbasis * pf
    return nothing
end
```

The old `oft!(out, jump::JumpOp, ...)` is replaced because it accesses `jump.in_eigenbasis` (an abstractly-typed field causing boxing). The new version takes the concrete-typed `eigenbasis` matrix directly.

The old `time_oft!` and `trotter_oft!` using `OFTCaches` are deprecated (marked so in code already). Remove them entirely; they exist only in `ofts.jl` for legacy test support. The tests can use `compute_jump_oft!` with the NUFFT path.

---

### TS-05: Workspace Consolidation

| Aspect | Detail |
|--------|--------|
| **Why Expected** | 5 workspace types with overlapping buffer names is confusing and wastes memory |
| **Complexity** | HIGH -- touches every simulation path |
| **Notes** | Must preserve zero-allocation hot paths |

**Concrete target: 3 workspace types**

**1. `SimWorkspace{T}` (replaces `KrausScratch` + `LindbladianWorkspace`)**

The DM-level workspace for `run_lindblad` and `run_thermalize`. Both need dim x dim scratch matrices.

```julia
struct SimWorkspace{T<:Complex}
    # Scratch matrices (dim x dim)
    jump_oft::Matrix{T}      # A(omega) buffer
    tmp1::Matrix{T}           # general scratch
    tmp2::Matrix{T}           # general scratch
    LdagL::Matrix{T}          # L'L product
    R::Matrix{T}              # R accumulator (Thermalize only, zeroed per step)
    rho_jump::Matrix{T}       # jump sandwich accumulator (Thermalize only)
    K0::Matrix{T}             # I - alpha*R (Thermalize only)
    rho_next::Matrix{T}       # output buffer

    # Vectorized Liouvillian extras (Lindblad only)
    Id::Matrix{T}             # identity
    jump_conj::Matrix{T}      # conj(jump) for vectorization
end
```

`LindbladianWorkspace` currently has 5 fields (Id, jump_tmp, jump_conj, jump_dag_jump, jump2_jump1). `KrausScratch` has 8 fields. The overlap is: `jump_tmp` = `jump_oft`, `jump_dag_jump` = `LdagL`. Combined: ~10 unique buffers. The `Lindblad`-specific ones (`Id`, `jump_conj`) are only used in vectorized Liouvillian construction, not in the hot path, so the small overhead of carrying them is acceptable.

Constructor:
```julia
function SimWorkspace(::Type{T}, dim::Int; vectorized::Bool=false) where {T<:Complex}
    # ... allocate all, skip Id/jump_conj if !vectorized
end
```

**2. `KrylovWorkspace{T,PD}` (keep as-is, just remove channel duplication)**

The Krylov workspace is already well-designed for its purpose. The only change is that `_accumulate_R_total!` and `_accumulate_jump_sandwich!` should call the unified `compute_jump_oft!` and `foreach_energy_pair` from the shared PRE layer instead of inlining them.

The workspace should also drop the `jumps::Vector{JumpOp}` field (only used for external inspection) and keep only `jump_eigenbases`/`jump_hermitian` which are the concrete-typed hot-path fields.

**3. `TrajectoryWorkspace{T}` (keep as-is)**

Already minimal: `jump_oft`, `psi_tmp`, `Rpsi`, `rho_acc`. No changes needed. The trajectory hot path operates on state vectors, not density matrices, so it correctly has different buffer shapes.

**Remove: `OFTCaches{T}`** -- deprecated, used only by `time_oft!` and `trotter_oft!` which are being removed.

---

### TS-06: Unified R Accumulation -- `accumulate_R!`

| Aspect | Detail |
|--------|--------|
| **Why Expected** | R computation appears in 3 places: `_precompute_R` (trajectories.jl x2 domains), `_accumulate_R_total!` (krylov_workspace.jl x3 domains). Same physics. |
| **Complexity** | MEDIUM -- domain dispatch needed |
| **Notes** | Uses `compute_jump_oft!` and `foreach_energy_pair` from TS-03/TS-04 |

**Concrete design:**

```julia
"""
    accumulate_R!(R, jump_eigenbases, jump_hermitian, precomputed_data, config, ham_or_trott)

Accumulate R_total = sum_k sum_w rate^2(w) * L_kw' * L_kw into R.
Dispatches on domain via config type parameter.
"""
function accumulate_R!(
    R::Matrix{T},
    eigenbases::Vector{Matrix{T}},
    hermitian_flags::Vector{Bool},
    precomputed_data,
    config::SimConfig{S,EnergyDomain},
    hamiltonian::HamHam,
) where {T<:Complex, S}
    (; transition, gamma_norm_factor, energy_labels) = precomputed_data
    prefactor = domain_prefactor(config, gamma_norm_factor)
    bohr_freqs = hamiltonian.bohr_freqs
    inv_4sigma2 = 1.0 / (4 * config.sigma^2)

    jump_oft = zeros(T, size(R)...)
    LdagL = zeros(T, size(R)...)

    for (k, eigenbasis) in enumerate(eigenbases)
        foreach_energy_pair(energy_labels, hermitian_flags[k],
            w -> begin
                compute_jump_oft!(jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)
                mul!(LdagL, jump_oft', jump_oft)
                @. R += prefactor * transition(w) * LdagL
            end,
            w -> begin
                mul!(LdagL, jump_oft, jump_oft')
                @. R += prefactor * transition(-w) * LdagL
            end
        )
    end
    hermitianize!(R)
    return R
end
```

The `_precompute_R` in trajectories.jl (which takes a `KrausScratch`) wraps this by pointing `eigenbases = [j.in_eigenbasis for j in jumps]`. The Krylov version already has `ws.jump_eigenbases`. Both call the same inner `accumulate_R!`.

---

### TS-07: Sandwich Function Consolidation

| Aspect | Detail |
|--------|--------|
| **Why Expected** | 6 sandwich functions in krylov_matvec.jl, but only 2 unique computations |
| **Complexity** | LOW -- algebraic identity elimination |
| **Notes** | The adjoint functions are literally identical to their non-adjoint counterparts |

**Current situation:**
- `_accumulate_sandwich!(out, L, rho, s, ws)` = `s * conj(L) * rho * L^T`
- `_accumulate_sandwich_adj_L!(out, L, rho, s, ws)` = `s * L^T * rho * conj(L)`
- `_accumulate_adjoint_sandwich!(out, L, rho, s, ws)` = `s * L^T * rho * conj(L)` -- **identical to `_adj_L`**
- `_accumulate_adjoint_sandwich_adj_L!(out, L, rho, s, ws)` = `s * conj(L) * rho * L^T` -- **identical to base `_sandwich`**

**Consolidation:** Keep exactly 2 functions:

```julia
"""Accumulate s * conj(L) * rho * L^T (forward sandwich / adjoint neg-freq)"""
@inline function sandwich_forward!(out, L, rho, s, ws)

"""Accumulate s * L^T * rho * conj(L) (forward neg-freq / adjoint sandwich)"""
@inline function sandwich_adjoint!(out, L, rho, s, ws)
```

The 2-operator Bohr domain variant `_accumulate_sandwich_2op!` stays separate since it takes two different operators A and B_dag.

---

### TS-08: Result Struct Pattern with Save

| Aspect | Detail |
|--------|--------|
| **Why Expected** | User explicitly requested 4 result types with optional save capability |
| **Complexity** | MEDIUM -- new save infrastructure, backward-compatible BSON |
| **Notes** | Replaces the current `ExperimentResult` which only wraps trajectory results |

**Concrete design -- 4 result types:**

```julia
# Each result type has a `save` kwarg that triggers BSON serialization.
# Common metadata is captured uniformly.

struct LindbladResult{T<:AbstractFloat}
    liouvillian::Matrix{Complex{T}}
    fixed_point::Matrix{Complex{T}}
    gap_mode::Matrix{Complex{T}}
    spectral_gap::Complex{T}
    config::SimConfig
    metadata::Dict{Symbol, Any}
end

struct ThermalizeResult{T<:AbstractFloat}
    final_dm::Matrix{Complex{T}}
    trace_distances::Vector{T}
    time_steps::Vector{T}
    config::SimConfig
    metadata::Dict{Symbol, Any}
end

struct KrylovSpectrumResult{T<:AbstractFloat}
    eigenvalues::Vector{Complex{T}}
    spectral_gap::T
    fixed_point::Matrix{Complex{T}}
    gap_mode::Matrix{Complex{T}}
    converged::Int
    matvec_count::Int
    num_restarts::Int
    normres::Vector{T}
    channel_eigenvalues::Union{Nothing, Vector{Complex{T}}}
    delta_used::Union{Nothing, T}
    config::SimConfig
    metadata::Dict{Symbol, Any}
end

struct TrajectoryResult{T<:AbstractFloat}
    rho_mean::Matrix{Complex{T}}
    n_trajectories::Int
    seed::Int
    times::Union{Nothing, Vector{Float64}}
    measurements_mean::Union{Nothing, Matrix{Float64}}
    convergence::Union{Nothing, ConvergenceData}
    config::SimConfig
    metadata::Dict{Symbol, Any}
end
```

**Save protocol:** Each `run_*` function gets an optional `save_path::Union{Nothing,String}=nothing` kwarg. If non-nothing, the result is serialized via:

```julia
function save_result(result::Union{LindbladResult, ThermalizeResult, KrylovSpectrumResult, TrajectoryResult}, path::String)
    d = _result_to_dict(result)
    mkpath(dirname(path))
    BSON.bson(path, d)
    _write_companion_txt(result, replace(path, ".bson" => ".txt"))
    return path
end
```

The `ExperimentResult` wrapper becomes unnecessary -- each result type carries its own config and metadata.

---

### TS-09: File Renaming -- PRE/MID/POST Grouping

| Aspect | Detail |
|--------|--------|
| **Why Expected** | Flat `src/` with 28 files needs a logical naming scheme |
| **Complexity** | LOW -- pure renames, no logic changes |
| **Notes** | Keep flat directory, use numeric prefix for include order |

**File rename mapping:**

```
# PRE: Physics building blocks (included first)
constants.jl           -> 01_constants.jl          (unchanged content)
structs.jl             -> 02_types.jl              (SimConfig, domains, JumpOp, workspaces)
hamiltonian.jl         -> 03_hamiltonian.jl         (HamHam, TrottTrott, Pauli, Trotter)
qi_tools.jl            -> 04_qi_tools.jl            (trace_distance, fidelity, hermitianize!)
misc_tools.jl          -> 05_misc_tools.jl          (load_hamiltonian, config validation, labels)
nufft.jl               -> 06_nufft.jl               (NUFFTPrefactors)
NEW                    -> 07_prefactors.jl          (domain_prefactor, compute_jump_oft!, foreach_energy_pair)
energy_domain.jl       -> 08_transitions.jl         (pick_transition, create_alpha -- merge energy/bohr domain tools)
bohr_domain.jl         -> 08_transitions.jl         (merged into above)
time_domain.jl         -> 08_transitions.jl         (merged into above)
trotter_domain.jl      -> 09_trotter.jl             (trotterize, TrottTrott constructor)
coherent.jl            -> 10_coherent.jl            (B operators, coherent precomputation)

# MID: Simulation runners (the 4 run_* functions + their inner helpers)
ofts.jl                -> REMOVE (deprecated, logic absorbed into 07_prefactors.jl)
kraus.jl               -> ABSORB into 02_types.jl (just a struct definition)
jump_workers.jl        -> 20_jump_workers.jl        (vectorized Liouvillian _jump_contribution!)
furnace_utensils.jl    -> 21_precompute.jl          (_precompute_data, _precompute_labels)
furnace.jl             -> 22_run_lindblad.jl         (run_lindblad, construct_lindbladian)
                       -> 23_run_thermalize.jl       (run_thermalize, split from furnace.jl)
krylov_workspace.jl    -> 24_krylov_workspace.jl     (KrylovWorkspace constructor)
krylov_matvec.jl       -> 25_krylov_matvec.jl        (apply_lindbladian!, apply_adjoint!)
krylov_eigsolve.jl     -> 26_run_krylov_spectrum.jl  (krylov_spectral_gap -> run_krylov_spectrum)
trajectories.jl        -> 27_run_trajectory.jl       (run_trajectory, TrajectoryFramework, step!)

# POST: Results and analysis
results.jl             -> 40_results.jl              (result types, save/load, BSON serialization)
convergence.jl         -> 41_convergence.jl          (ConvergenceData, adaptive trajectory runs)
fitting.jl             -> 42_fitting.jl              (exponential decay fitting)
gap_estimation.jl      -> 43_gap_estimation.jl       (trajectory-based gap estimation)
diagnostics.jl         -> 44_diagnostics.jl          (exact diagnostics suite)
log_sobolev.jl         -> 45_log_sobolev.jl          (LSI framework)
errors.jl              -> 50_errors.jl               (custom error types)
```

**Include order in `QuantumFurnace.jl`:** The numeric prefix guarantees correct include order. Files in the same prefix group (e.g., 08_*) can be included in any order since they are in the same dependency tier.

**Unfinished trajectory-gap code:** The `SpectralGapResult` and `estimate_spectral_gap` in gap_estimation.jl, along with `ObservableTrajectoryResult`, move to a `staging/` directory outside `src/`. They are not part of the 4-runner architecture.

---

### TS-10: `_thermalize_to_liouv_config` Elimination

| Aspect | Detail |
|--------|--------|
| **Why Expected** | With unified config, there is no "conversion" needed |
| **Complexity** | FREE -- falls out of TS-01 |
| **Notes** | Currently 2 functions x 14 field-by-field copies |

With `SimConfig{S,D,DB,T}`, the Krylov workspace constructor for a thermalize config simply reads the Lindbladian-relevant fields directly from the same struct. No conversion needed. The dispatch on `AbstractLiouvConfig` vs `AbstractThermalizeConfig` becomes dispatch on `SimConfig{Lindblad}` vs `SimConfig{Thermalize}`, and the shared physics parameters are accessed identically.

---

## Differentiators

Features that make the restructure excellent rather than merely adequate.

### DIFF-01: Unified Dissipator Loop

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| Single `dissipator_loop!` function with callback | Replaces 6 domain-specific loops across jump_workers, krylov_matvec, krylov_eigsolve | HIGH | The "holy grail" deduplication |

**Concept:** A single inner loop that iterates over jumps and energies, computing the OFT for each, and calls back with the Lindblad operator L_w and scalar rate for each (w, jump) pair. The callback implements what to do with L_w:
- Vectorized Liouvillian: `_vectorize_liouv_diss_and_add!(L_target, L_w, scalar, ws)`
- Krylov sandwich: `_accumulate_sandwich!(rho_out, L_w, rho, scalar, ws)`
- DM thermalization: accumulate into `R` and `rho_jump`
- Trajectory step: compute jump probability and sample

**Risk:** The DM thermalization and trajectory paths have fundamentally different per-iteration logic (DM accumulates R+rho_jump, trajectory samples). Forcing them into a single callback may be awkward. Better to have `dissipator_loop!` handle the "iterate jumps x energies, compute L_w" part, and let callers handle per-L_w actions. This still eliminates the 16x hermitian branching duplication.

---

### DIFF-02: Config Factory Functions for Tests

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| `make_test_config(; sim=Lindblad(), domain=EnergyDomain(), db=KMS(), ...)` | Eliminates 20+ ad-hoc config constructions across test files | LOW | Pure convenience |

```julia
# In test/test_helpers.jl:
function make_test_config(;
    sim::AbstractSimType = Lindblad(),
    domain::AbstractDomain = EnergyDomain(),
    db::AbstractDBType = KMS(),
    num_qubits = NUM_QUBITS,
    beta = BETA,
    sigma = SIGMA,
    kwargs...
)
    SimConfig(
        sim_type = sim, domain = domain, db_type = db,
        num_qubits = num_qubits, beta = beta, sigma = sigma,
        with_linear_combination = false,
        num_energy_bits = 6, t0 = ..., w0 = ...,
        kwargs...
    )
end
```

---

### DIFF-03: Precomputed Data as Typed Struct

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| Replace NamedTuple precomputed_data with typed structs | Eliminates `hasproperty` checks, enables dispatch | MEDIUM | Optional but clean |

Currently `_precompute_data` returns ad-hoc NamedTuples with different fields per domain. The `TrajectoryFramework` constructor checks `hasproperty(precomputed_data, :oft_nufft_prefactors)`. With typed structs:

```julia
struct EnergyPrecomputed{T,F}
    transition::F
    gamma_norm_factor::T
    energy_labels::Vector{T}
end

struct TimePrecomputed{T,F,P}
    transition::F
    gamma_norm_factor::T
    energy_labels::Vector{T}
    oft_nufft_prefactors::P
    b_minus::Union{Nothing, Dict}
    b_plus::Union{Nothing, Dict}
end
```

This makes `_precompute_data` return type-stable and enables dispatch on precomputed data type in the dissipator loop.

---

## Anti-Features

Features to explicitly NOT build during this restructure.

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| Module/subdirectory splitting of `src/` | Julia's module system adds compilation overhead for internal modules; flat `src/` with numeric prefixes is idiomatic for packages of this size | Use numeric file prefixes for logical grouping |
| Generic abstract result type with parametric fields | Each simulation produces fundamentally different output (Liouvillian matrix vs trace distances vs eigenvalues vs trajectories); forcing them into one type is Procrustean | Keep 4 concrete result types, share metadata capture and save protocol via functions |
| Automatic config validation via inner constructor | Inner constructors with validation make testing painful (cannot construct partial/invalid configs for edge case tests) | Keep `validate_config!` as explicit pre-flight check in each `run_*` function |
| Trait-based dispatch instead of type parameters | Julia traits (Holy traits) add complexity without benefit here; `SimConfig{S,D,DB,T}` with 4 type params already gives precise dispatch | Use parametric dispatch directly |
| Observable/convergence tracking in `run_thermalize` | The DM thermalization is for small-scale validation; observable tracking belongs in trajectory simulations | Keep `run_thermalize` minimal; trajectory code handles observables |
| Moving `TrajectoryFramework` construction into workspace | The framework pre-computes per-operator Kraus data (eigendecompositions) which is expensive and should be explicit | Keep `build_trajectoryframework` as explicit construction step |

---

## Feature Dependencies

```
TS-01 (Config)
  |
  +---> TS-02 (Prefactors) -- needs new config type for dispatch
  |       |
  |       +---> TS-03 (Half-grid iterator) -- no config dependency, but logically grouped
  |       |
  |       +---> TS-04 (OFT unification) -- no config dependency, but logically grouped
  |       |
  |       +---> TS-06 (R accumulation) -- uses TS-02, TS-03, TS-04
  |
  +---> TS-05 (Workspace) -- needs new config type
  |
  +---> TS-07 (Sandwich) -- independent, can parallelize
  |
  +---> TS-08 (Results) -- needs new config type for embedding
  |
  +---> TS-09 (File rename) -- do LAST, after all logic changes
  |
  +---> TS-10 (Config conversion elimination) -- FREE with TS-01

TS-03 + TS-04 can be done BEFORE TS-01 as pure extractions.
TS-07 can be done BEFORE TS-01 (pure algebraic simplification).
```

**Critical path:** TS-01 -> TS-05 -> TS-06 -> TS-08 -> TS-09

**Parallelizable early wins:** TS-02, TS-03, TS-04, TS-07 (can be done before TS-01 using current types)

---

## MVP Recommendation

### Phase 1 (Foundation -- do first, unblocks everything):

1. **TS-02: Prefactor extraction** -- 30 min, zero risk, immediate deduplication
2. **TS-03: Half-grid iterator** -- 1 hr, eliminates 16 code blocks
3. **TS-04: OFT unification** -- 30 min, kills 3 OFT variants
4. **TS-07: Sandwich consolidation** -- 30 min, kills 4 functions

These 4 can be done on the CURRENT type system. They reduce code volume by ~400 lines before touching any struct definitions. This is the safest starting point.

### Phase 2 (Type system -- the big change):

5. **TS-01: Unified Config** -- the hardest single change. Do it in one focused session. Write the struct, update all dispatch sites, update serialization, update tests. The config change touches every file.
6. **TS-10: Config conversion elimination** -- free with TS-01.

### Phase 3 (Workspace and results):

7. **TS-05: Workspace consolidation** -- after config is stable
8. **TS-06: R accumulation unification** -- uses new workspaces
9. **TS-08: Result struct pattern** -- after config and workspace are stable

### Phase 4 (Cleanup):

10. **TS-09: File rename** -- absolutely last. Git blame preservation via `git mv`.

### Defer:

- **DIFF-01: Unified dissipator loop** -- too risky for this milestone. The 4 simulation paths have enough structural differences that a single callback-based loop would be fragile. The prefactor + half-grid + OFT deduplication (TS-02/03/04) captures 80% of the benefit at 20% of the risk.
- **DIFF-03: Typed precomputed data** -- nice-to-have, defer to a follow-up milestone.
- **Trajectory-gap code staging** -- move `estimate_spectral_gap`, `SpectralGapResult`, `ObservableTrajectoryResult`, and the `trajectory_validation/` test directory to `staging/` at the end.

---

## Runner Function Signatures (Target API)

```julia
# The 4 public runner functions after restructure:

function run_lindblad(
    jumps::Vector{JumpOp},
    config::SimConfig{Lindblad, D, DB, T},
    hamiltonian::HamHam;
    trotter::Union{TrottTrott, Nothing} = nothing,
    save_path::Union{Nothing, String} = nothing,
) -> LindbladResult{T}

function run_thermalize(
    jumps::Vector{JumpOp},
    config::SimConfig{Thermalize, D, DB, T},
    evolving_dm::Matrix{Complex{T}},
    hamiltonian::HamHam;
    trotter::Union{TrottTrott, Nothing} = nothing,
    rng::AbstractRNG = Random.default_rng(),
    save_path::Union{Nothing, String} = nothing,
) -> ThermalizeResult{T}

function run_krylov_spectrum(
    config::SimConfig{KrylovSpectrum, D, DB, T},
    hamiltonian::HamHam,
    jumps::Vector{JumpOp};
    trotter::Union{TrottTrott, Nothing} = nothing,
    krylovdim::Int = 30,
    howmany::Int = 4,
    tol::Real = 1e-10,
    save_path::Union{Nothing, String} = nothing,
) -> KrylovSpectrumResult{T}

function run_trajectory(
    jumps::Vector{JumpOp},
    config::SimConfig{Trajectory, D, DB, T},
    psi0::Vector{Complex{T}},
    hamiltonian::HamHam;
    trotter::Union{TrottTrott, Nothing} = nothing,
    ntraj::Int = 1,
    observables::Union{Nothing, Vector{<:Matrix}} = nothing,
    save_every::Int = 1,
    seed::Union{Int, Nothing} = nothing,
    save_path::Union{Nothing, String} = nothing,
) -> TrajectoryResult{T}
```

Note: `run_krylov_spectrum` accepts both `SimConfig{KrylovSpectrum,D,DB}` (Lindbladian path) and can be called with a thermalize-like config by constructing `SimConfig{KrylovSpectrum,...}` with `delta` populated. The dispatch between Lindbladian eigsolve and channel eigsolve is determined by whether `config.delta` is nothing.

---

## Sources

- Exhaustive reading of all 28 source files in `/Users/bence/code/QuantumFurnace.jl/src/`
- User's architecture notes in `supplementary-informations/quantumfurnace-structure.md`
- Current test infrastructure in `test/test_helpers.jl` and `test/runtests.jl`
- Julia language documentation on parametric types, multiple dispatch, and `@kwdef`
- Existing research files in `.planning/research/` from the Krylov milestone
