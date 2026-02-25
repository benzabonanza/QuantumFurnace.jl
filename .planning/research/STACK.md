# Technology Stack: QuantumFurnace.jl Codebase Restructure

**Project:** QuantumFurnace.jl -- Parametric type hierarchy, workspace consolidation, code deduplication
**Researched:** 2026-02-25
**Confidence:** HIGH (patterns verified against Julia official docs, SciML patterns, and codebase analysis)

---

## Scope

This STACK.md covers Julia language patterns, idioms, and tools for restructuring the QuantumFurnace.jl codebase (8,312 LOC src + 5,071 LOC test). It does NOT re-research the existing dependency stack (KrylovKit, FINUFFT, Arpack, etc.), which remains unchanged.

Focus areas:
1. **Parametric type hierarchy** for Config{S,D,DB,T} -- encoding Simulation, Domain, DetailedBalance, and Float type
2. **Workspace consolidation** -- reducing 5 workspace structs to 2-3
3. **Code deduplication** -- eliminating 8+ copies of prefactor formulas and sandwich loops
4. **File organization** -- grouping 28 source files into logical sections

No new package dependencies are needed. This is a pure refactoring milestone using Julia's existing type system, multiple dispatch, and module organization.

---

## 1. Parametric Config Type Hierarchy

### Current Problem

Four config structs (`LiouvConfig`, `LiouvConfigGNS`, `ThermalizeConfig`, `ThermalizeConfigGNS`) duplicate 13 identical fields. The only differences are:
- `ThermalizeConfig*` adds `mixing_time::T` and `delta::T`
- `*GNS` forces `with_coherent = false` and uses unshifted transition weights
- Dispatch currently uses `Union{LiouvConfig, ThermalizeConfig}` in 10+ places

Functions like `pick_transition`, `_pick_alpha`, `B_bohr`, `_select_b_plus_calculator` all manually dispatch on concrete config types with near-identical bodies.

### Recommended Pattern: Three-Level Parametric Hierarchy

Use Julia's parametric abstract types to encode orthogonal axes (Simulation type, Domain, Detailed Balance flavor) as type parameters, enabling dispatch on any axis independently.

**Confidence: HIGH** -- This is the standard Julia pattern used by DiffEqBase.jl (algorithms as parametric types), KrylovKit.jl (Arnoldi algorithm types), and recommended in Julia Performance Tips.

```julia
# ---- Abstract hierarchy ----
# Detailed Balance flavors (singleton types, like Domain types)
abstract type AbstractDBFlavor end
struct KMSDB <: AbstractDBFlavor end   # Exact KMS detailed balance
struct GNSDB <: AbstractDBFlavor end   # Approx GNS detailed balance
# Future: struct DLLDB <: AbstractDBFlavor end

# Simulation mode (what kind of output/computation)
abstract type AbstractSimMode end
struct Lindbladian <: AbstractSimMode end    # Dense Liouvillian construction + eigen
struct Thermalize <: AbstractSimMode end     # DM stepping (weak measurement emulation)
struct KrylovSpectrum <: AbstractSimMode end # Matrix-free Krylov eigsolve
struct Trajectory <: AbstractSimMode end     # Monte Carlo state-vector sampling

# Config hierarchy: parametric on all axes
abstract type AbstractConfig{S<:AbstractSimMode, D<:AbstractDomain, DB<:AbstractDBFlavor, T<:AbstractFloat} end
```

**Key design decision: Simulation and DB as type parameters, NOT struct fields.**

Why type parameters:
- The compiler specializes functions for each concrete combination at compile time
- Zero-cost dispatch: `f(::AbstractConfig{Lindbladian, EnergyDomain, KMSDB})` compiles to a single method body
- Type stability: the compiler knows the exact config type throughout the call chain
- Cannot accidentally mix incompatible combinations at runtime

Why NOT Val-based dispatch for DB:
- `Val{:KMS}` requires explicit `Val` wrapping at every call site
- Type parameters propagate automatically through the type system
- Val dispatch adds latency (each new value is a new type) with no benefit over parametric types here

### Concrete Config Struct Design

```julia
"""
Core physical parameters shared by ALL simulation modes.
Stored as a separate struct for composition (not inheritance).
"""
@kwdef struct PhysicsParams{D<:AbstractDomain, T<:AbstractFloat}
    num_qubits::Int
    domain::D
    beta::T
    sigma::T
    with_linear_combination::Bool
    # Linear combination parameters (concrete Union -- isbits optimization applies)
    gaussian_parameters::Tuple{T, T} = (zero(T), zero(T))
    a::T = zero(T)
    b::T = zero(T)
    eta::T = zero(T)
    # Discretization parameters (required for Energy/Time/Trotter, ignored for Bohr)
    num_energy_bits::Int = 0
    t0::T = zero(T)
    w0::T = zero(T)
    num_trotter_steps_per_t0::Int = 0
end

"""
    SimConfig{S, D, DB, T}

Unified simulation configuration. The type parameters encode:
- S: simulation mode (Lindbladian, Thermalize, KrylovSpectrum, Trajectory)
- D: domain (BohrDomain, EnergyDomain, TimeDomain, TrotterDomain)
- DB: detailed balance flavor (KMSDB, GNSDB)
- T: float precision (Float64, Float32)
"""
struct SimConfig{S<:AbstractSimMode, D<:AbstractDomain, DB<:AbstractDBFlavor, T<:AbstractFloat}
    physics::PhysicsParams{D, T}
    with_coherent::Bool
    # Thermalize/Trajectory-specific (nothing for Lindbladian/KrylovSpectrum)
    mixing_time::Union{T, Nothing}
    delta::Union{T, Nothing}
end
```

**Critical: The `Union{T, Nothing}` for mixing_time/delta is acceptable here** because:
1. These fields are accessed ONCE during setup, never in the hot path
2. Julia's isbits union optimization applies to `Union{Float64, Nothing}` (stored inline, no heap allocation)
3. The compiler can union-split on 2 concrete types efficiently

**Confidence: HIGH** -- Julia docs explicitly state isbits Union optimization handles `Union{T, Nothing}` efficiently for isbits `T`.

### Accessor Functions (The Julian Interface Pattern)

Instead of directly accessing fields, define accessor functions on the abstract type. This decouples callers from struct layout and enables future changes.

```julia
# Physics accessors -- work on any config
num_qubits(c::SimConfig) = c.physics.num_qubits
domain(c::SimConfig) = c.physics.domain
beta(c::SimConfig) = c.physics.beta
sigma(c::SimConfig) = c.physics.sigma
w0(c::SimConfig) = c.physics.w0
t0(c::SimConfig) = c.physics.t0
num_energy_bits(c::SimConfig) = c.physics.num_energy_bits

# Simulation-mode accessors
mixing_time(c::SimConfig{S}) where {S<:Union{Thermalize, Trajectory}} = c.mixing_time
delta_step(c::SimConfig{S}) where {S<:Union{Thermalize, Trajectory}} = c.delta
with_coherent(c::SimConfig) = c.with_coherent

# DB-flavor query (compile-time known)
is_gns(::SimConfig{S,D,GNSDB}) where {S,D} = true
is_gns(::SimConfig) = false
```

### Dispatch on Individual Axes

The real power: dispatch on exactly the axis that matters.

```julia
# Dispatch on DB flavor only (transition function differs)
function pick_transition(c::SimConfig{S, D, KMSDB}) where {S, D}
    # KMS-shifted transition weight
    _pick_transition_kms(c)
end
function pick_transition(c::SimConfig{S, D, GNSDB}) where {S, D}
    # Unshifted GNS transition weight
    _pick_transition_gns(c)
end

# Dispatch on Domain only (prefactor formula differs)
function compute_prefactor(c::SimConfig{S, EnergyDomain}, gnf::Real) where {S}
    (w0(c) / (sigma(c) * sqrt(2pi))) * gnf
end
function compute_prefactor(c::SimConfig{S, D}, gnf::Real) where {S, D<:Union{TimeDomain, TrotterDomain}}
    w0(c) * t0(c)^2 * (sigma(c) * sqrt(2/pi)) / (2pi) * gnf
end

# Dispatch on Simulation mode (what to return/compute)
function run_simulation(c::SimConfig{Lindbladian}, ...)
    # Dense Liouvillian + Arpack eigen
end
function run_simulation(c::SimConfig{Thermalize}, ...)
    # DM stepping with Kraus channel
end
function run_simulation(c::SimConfig{KrylovSpectrum}, ...)
    # Matrix-free Krylov eigsolve
end
function run_simulation(c::SimConfig{Trajectory}, ...)
    # Monte Carlo trajectory sampling
end

# Dispatch on DB + coherent interaction
function _precompute_coherent_total_B(c::SimConfig{S, D, GNSDB}, ...) where {S, D}
    return nothing  # GNS never has coherent term
end
function _precompute_coherent_total_B(c::SimConfig{S, D, KMSDB}, ...) where {S, D}
    c.with_coherent || return nothing
    # ... compute B ...
end
```

### Backward Compatibility: Type Aliases

Maintain existing API names as type aliases during transition:

```julia
const LiouvConfig{D, T} = SimConfig{Lindbladian, D, KMSDB, T}
const LiouvConfigGNS{D, T} = SimConfig{Lindbladian, D, GNSDB, T}
const ThermalizeConfig{D, T} = SimConfig{Thermalize, D, KMSDB, T}
const ThermalizeConfigGNS{D, T} = SimConfig{Thermalize, D, GNSDB, T}
```

### Convenience Constructors

```julia
function LiouvConfig(; num_qubits, domain::D, beta::T, sigma::T, with_coherent=true,
                       with_linear_combination=false, kwargs...) where {D, T}
    physics = PhysicsParams{D, T}(; num_qubits, domain, beta, sigma,
                                    with_linear_combination, kwargs...)
    SimConfig{Lindbladian, D, KMSDB, T}(physics, with_coherent, nothing, nothing)
end

function ThermalizeConfig(; num_qubits, domain::D, beta::T, sigma::T,
                            mixing_time::T, delta::T, with_coherent=true,
                            with_linear_combination=false, kwargs...) where {D, T}
    physics = PhysicsParams{D, T}(; num_qubits, domain, beta, sigma,
                                    with_linear_combination, kwargs...)
    SimConfig{Thermalize, D, KMSDB, T}(physics, with_coherent, mixing_time, delta)
end
```

### What This Eliminates

| Current Duplication | After Restructure |
|---|---|
| 4 config struct definitions (52 duplicate field lines) | 1 SimConfig + 1 PhysicsParams |
| `_thermalize_to_liouv_config()` (2 methods, 40 lines) | Unnecessary -- dispatch on S parameter directly |
| `pick_transition` (4 methods with Union types) | 2 methods dispatching on DB axis |
| `_pick_alpha` (4 methods with Union types) | 2 methods dispatching on DB axis |
| `_pick_f` (Union{LiouvConfig, ThermalizeConfig} argument) | Dispatches on DB, domain irrelevant |
| `_select_b_plus_calculator` (Union type argument) | 1 method, DB axis dispatches to correct formula |

---

## 2. Workspace Consolidation

### Current Problem

Five workspace structs with overlapping scratch buffers:

| Workspace | Fields | Used By | Unique Fields |
|---|---|---|---|
| `LindbladianWorkspace` | Id, jump_tmp, jump_conj, jump_dag_jump, jump2_jump1 | Dense Liouvillian construction | Id, jump_conj, jump2_jump1 |
| `KrausScratch` | jump_oft, LdagL, R, rho_jump, K0, tmp1, tmp2, rho_next | DM thermalization | R, rho_jump, K0, rho_next |
| `KrylovWorkspace` | jump_oft, tmp1, tmp2, LdagL, rho_out + precomputed data | Krylov matvec | rho_out, precomputed_data, G_left/G_right, channel fields |
| `TrajectoryWorkspace` | jump_oft, psi_tmp, Rpsi, rho_acc | Trajectory stepping | psi_tmp, Rpsi, rho_acc |
| `OFTCaches` | prefactors, U, temp_op | Debugging OFT | prefactors, U |

Common scratch across all: `jump_oft` (dim x dim), `tmp1` (dim x dim), `LdagL` (dim x dim).

### Recommended Pattern: Composition of Scratch + Precomputed Data

Use a common scratch buffer struct (pure preallocated arrays) composed into simulation-specific workspaces that also hold precomputed data.

**Confidence: HIGH** -- This is the standard Julia composition pattern recommended by the community. Avoids inheritance (which Julia does not support for concrete types), keeps scratch buffers concrete-typed, and allows the inner scratch to be shared.

```julia
"""
Common scratch buffers for dim x dim complex matrix operations.
Used by ALL simulation modes. Contains only preallocated arrays,
no precomputed physics data.
"""
struct MatrixScratch{T<:Complex}
    jump_oft::Matrix{T}   # A(omega) buffer
    tmp1::Matrix{T}       # generic scratch
    tmp2::Matrix{T}       # generic scratch
    LdagL::Matrix{T}      # L'*L product
end

function MatrixScratch(::Type{T}, dim::Int) where {T<:Complex}
    Zm() = zeros(T, dim, dim)
    MatrixScratch{T}(Zm(), Zm(), Zm(), Zm())
end

"""
Scratch buffers for density-matrix based simulations (Thermalize, dense Lindbladian).
Extends MatrixScratch with DM-specific buffers.
"""
struct DMScratch{T<:Complex}
    core::MatrixScratch{T}  # Composition, not inheritance
    R::Matrix{T}            # R_total accumulator
    rho_jump::Matrix{T}     # Jump sandwich accumulator
    K0::Matrix{T}           # Kraus K0 operator
    rho_next::Matrix{T}     # Next DM step
end

function DMScratch(::Type{T}, dim::Int) where {T<:Complex}
    Zm() = zeros(T, dim, dim)
    DMScratch{T}(MatrixScratch(T, dim), Zm(), Zm(), Zm(), Zm())
end

"""
Scratch buffers for trajectory simulations (state-vector based).
"""
struct TrajScratch{T<:Complex}
    core::MatrixScratch{T}  # Shared matrix scratch
    psi_tmp::Vector{T}      # State-vector scratch
    Rpsi::Vector{T}         # R*psi cache
    rho_acc::Matrix{T}      # DM accumulator for averaging
end

function TrajScratch(::Type{T}, dim::Int) where {T<:Complex}
    TrajScratch{T}(MatrixScratch(T, dim), zeros(T, dim), zeros(T, dim), zeros(T, dim, dim))
end
```

### Forwarding Access (Zero-Overhead)

Use `@inline` property access forwarding so callers can access `ws.jump_oft` regardless of workspace type:

```julia
# Forward core scratch access -- compiler inlines to direct field access
@inline jump_oft(ws::DMScratch) = ws.core.jump_oft
@inline tmp1(ws::DMScratch) = ws.core.tmp1
@inline tmp2(ws::DMScratch) = ws.core.tmp2
@inline LdagL(ws::DMScratch) = ws.core.LdagL

# Or use an abstract interface with getproperty overloading:
# (This pattern is used by ComponentArrays.jl and similar)
@inline Base.getproperty(ws::DMScratch, s::Symbol) = begin
    s === :core && return getfield(ws, :core)
    s === :R && return getfield(ws, :R)
    s === :rho_jump && return getfield(ws, :rho_jump)
    s === :K0 && return getfield(ws, :K0)
    s === :rho_next && return getfield(ws, :rho_next)
    # Forward to core scratch
    return getproperty(getfield(ws, :core), s)
end
```

**Recommended approach: explicit accessor functions over `getproperty` overloading.** The `getproperty` approach is clever but makes code harder to grep and can confuse IDE tools. Accessor functions are explicit and the compiler inlines them to zero overhead.

However, for the codebase as-is where callers use `ws.jump_oft` extensively, the simpler approach is to NOT use composition for the hot-path scratch and instead just define unified workspace structs with all needed fields:

### Pragmatic Approach (Recommended)

```julia
"""
    SimWorkspace{T, PD}

Unified workspace for Krylov/DM/Lindbladian simulations.
Contains precomputed physics data AND scratch buffers.
PD is the concrete NamedTuple type of precomputed_data (varies by domain).
"""
struct SimWorkspace{T<:Complex, PD<:NamedTuple}
    # Precomputed physics (immutable after construction)
    precomputed_data::PD
    B_total::Union{Nothing, Matrix{T}}
    jumps::Vector{JumpOp}
    jump_eigenbases::Vector{Matrix{T}}
    jump_hermitian::Vector{Bool}

    # Core scratch (shared by all simulation paths)
    jump_oft::Matrix{T}
    tmp1::Matrix{T}
    tmp2::Matrix{T}
    LdagL::Matrix{T}
    rho_out::Matrix{T}       # Krylov output / generic output

    # DM-specific (nothing for Lindbladian-only path)
    R::Union{Nothing, Matrix{T}}
    rho_jump::Union{Nothing, Matrix{T}}
    K0::Union{Nothing, Matrix{T}}
    rho_next::Union{Nothing, Matrix{T}}

    # Phase 32 effective Hamiltonian (precomputed for optimized matvec)
    G_left::Union{Nothing, Matrix{T}}
    G_right::Union{Nothing, Matrix{T}}
    G_left_adj::Union{Nothing, Matrix{T}}
    G_right_adj::Union{Nothing, Matrix{T}}

    # Channel matrices (precomputed for ThermalizeConfig path)
    channel_K0::Union{Nothing, Matrix{T}}
    channel_U_residual::Union{Nothing, Matrix{T}}
    channel_U_coherent::Union{Nothing, Matrix{T}}
    channel_rho_jump::Union{Nothing, Matrix{T}}
    channel_delta::Union{Nothing, Float64}
end
```

**Why Union{Nothing, T} is acceptable for these fields:**
1. None of these `Union` fields are accessed in hot loops -- they are checked ONCE during setup or dispatch boundary
2. Julia's small-union optimization handles `Union{Nothing, Matrix{T}}` efficiently (tag byte + pointer, no boxing)
3. Functions that use these fields are behind a dispatch barrier (the caller already knows the simulation mode)

**What this consolidates:**
- `KrylovWorkspace` + `KrausScratch` + `LindbladianWorkspace` -> `SimWorkspace` (one struct, constructed with appropriate nothing fields per simulation mode)
- `TrajectoryWorkspace` stays separate (state-vector scratch is fundamentally different from density-matrix scratch)
- `OFTCaches` is eliminated (debug-only, inline the logic)

### Dense Lindbladian Path

The `LindbladianWorkspace` fields (`Id`, `jump_conj`, `jump2_jump1`) are only used by the vectorized Liouvillian construction. These can be local allocations inside `construct_lindbladian` since that function runs once per setup, not in a hot loop. No workspace struct needed.

---

## 3. Code Deduplication via Dispatch

### Current Problem: Prefactor Formula Duplication

The same prefactor computation appears in 8+ places:

```
EnergyDomain:     prefactor = (w0 / (sigma * sqrt(2pi))) * gamma_norm_factor
Time/Trotter:     prefactor = w0 * t0^2 * (sigma * sqrt(2/pi)) / (2pi) * gamma_norm_factor
```

Locations: `_jump_contribution!` (3 Liouvillian + 3 Thermalize), `_accumulate_R_total!` (3 methods), `_precompute_R` (2 methods), `apply_lindbladian!` (3 methods), `apply_adjoint_lindbladian!` (3 methods), `build_trajectoryframework` (1 method with if/else), `step_along_trajectory!` (2 methods).

### Recommended Pattern: Domain-Dispatched Helper Functions

Extract domain-specific computations into small `@inline` helper functions. The compiler will inline them, producing identical machine code to the current hand-duplicated version, but the source code lives in one place.

**Confidence: HIGH** -- Julia `@inline` annotation combined with dispatch on singleton types is the standard zero-overhead abstraction pattern. The compiler specializes and inlines these at compile time.

```julia
# ---- Domain-specific prefactor computation ----
# Called once per frequency loop, inlined to zero overhead

@inline function compute_rate_prefactor(
    ::EnergyDomain, w0::Real, t0::Real, sigma::Real, gnf::Real
)
    (w0 / (sigma * sqrt(2pi))) * gnf
end

@inline function compute_rate_prefactor(
    ::Union{TimeDomain, TrotterDomain}, w0::Real, t0::Real, sigma::Real, gnf::Real
)
    w0 * t0^2 * (sigma * sqrt(2/pi)) / (2pi) * gnf
end

# ---- Domain-specific OFT computation ----
# Replaces both oft!() and the NUFFT prefactor multiply

@inline function compute_jump_oft!(
    out::Matrix{T}, eigenbasis::Matrix{T},
    ::EnergyDomain, bohr_freqs::Matrix{<:Real},
    energy::Real, inv_4sigma2::Real, ::Nothing
) where {T<:Complex}
    @. out = eigenbasis * exp(-(energy - bohr_freqs)^2 * inv_4sigma2)
    return nothing
end

@inline function compute_jump_oft!(
    out::Matrix{T}, eigenbasis::Matrix{T},
    ::Union{TimeDomain, TrotterDomain}, ::Matrix{<:Real},
    energy::Real, ::Real, nufft_prefactors
) where {T<:Complex}
    pf = _prefactor_view(nufft_prefactors, energy)
    @. out = eigenbasis * pf
    return nothing
end
```

### Hermitian Symmetry Loop Deduplication

The "iterate half-grid for Hermitian jumps, full grid otherwise" pattern is duplicated in every frequency loop (Liouvillian, Krylov, DM, Trajectory). Extract this into a unified iteration pattern:

```julia
"""
    foreach_frequency(f_positive!, f_negative!, jump_herm, energy_labels)

Unified frequency iteration. For Hermitian jumps, iterates half-grid
(w <= 0) and calls both f_positive!(w) and f_negative!(w) for w > 0.
For non-Hermitian jumps, iterates full grid calling f_positive!(w) only.

This eliminates the duplicated if/else Hermitian branching across all
simulation modes.
"""
@inline function foreach_frequency(
    f_positive!::F1, f_negative!::F2,
    is_hermitian::Bool, energy_labels::AbstractVector{<:Real}
) where {F1, F2}
    if is_hermitian
        @inbounds for w_raw in energy_labels
            w_raw > 1e-12 && continue
            w = abs(w_raw)
            f_positive!(w)
            if w > 1e-12
                f_negative!(w)
            end
        end
    else
        @inbounds for w in energy_labels
            f_positive!(w)
        end
    end
    return nothing
end
```

**Why closures for f_positive!/f_negative!:** Julia specializes on the closure type, so the compiler inlines the closure body. The function parameter types `F1, F2` are concrete closure types, not `Function` (which would be abstract). This produces the same machine code as the hand-written loop.

**Confidence: HIGH** -- Verified against Julia Performance Tips: "Annotate values taken from untyped locations... Write 'type-stable' functions."

### Sandwich Loop Deduplication

The forward/adjoint sandwich accumulation pattern is currently 4 nearly-identical functions (`_accumulate_sandwich!`, `_accumulate_sandwich_adj_L!`, `_accumulate_adjoint_sandwich!`, `_accumulate_adjoint_sandwich_adj_L!`). The only difference is the order of conj/transpose operations.

```julia
"""
    _accumulate_sandwich!(out, L_op, rho, scalar, ws, ::Val{:forward})
    _accumulate_sandwich!(out, L_op, rho, scalar, ws, ::Val{:adjoint})

Unified sandwich: conj(L)*rho*L^T (forward) or L^T*rho*conj(L) (adjoint).
Val dispatch is eliminated at compile time.
"""
@inline function _accumulate_sandwich!(
    out::Matrix{T}, L_op::Matrix{T}, rho::Matrix{T},
    scalar::Real, ws, ::Val{:forward}
) where {T<:Complex}
    CT = one(T); ZT = zero(T)
    @. ws.tmp2 = conj(L_op)
    BLAS.gemm!('N', 'N', CT, ws.tmp2, rho, ZT, ws.tmp1)
    BLAS.gemm!('N', 'T', CT, ws.tmp1, L_op, ZT, ws.LdagL)
    BLAS.axpy!(T(scalar), ws.LdagL, out)
end

@inline function _accumulate_sandwich!(
    out::Matrix{T}, L_op::Matrix{T}, rho::Matrix{T},
    scalar::Real, ws, ::Val{:adjoint}
) where {T<:Complex}
    CT = one(T); ZT = zero(T)
    @. ws.tmp2 = conj(L_op)
    BLAS.gemm!('T', 'N', CT, L_op, rho, ZT, ws.tmp1)
    BLAS.gemm!('N', 'N', CT, ws.tmp1, ws.tmp2, ZT, ws.LdagL)
    BLAS.axpy!(T(scalar), ws.LdagL, out)
end
```

**Why Val here and not for DB:** Val dispatch is appropriate for binary forward/adjoint switching because:
1. The caller always knows at compile time whether it is forward or adjoint
2. Only 2 values, and they are structural (operation direction), not domain-semantic
3. The Val is passed as a constant, so the compiler eliminates the dispatch entirely

### Unified apply_lindbladian! via Domain Dispatch

The 6 `apply_lindbladian!` / `apply_adjoint_lindbladian!` methods (2 per domain: Energy, Time/Trotter, Bohr) share 90% of their code. The only differences are:
1. How jump_oft is computed (Gaussian filter vs NUFFT prefactor vs alpha-weighted)
2. The prefactor formula
3. BohrDomain uses a 2-operator sandwich with A_nu2_dag

```julia
function _apply_lindbladian_impl!(
    ws::SimWorkspace{T}, rho::Matrix{T},
    config::SimConfig{S, D}, hamiltonian::HamHam,
    ::Val{direction}
) where {T<:Complex, S, D<:Union{EnergyDomain, TimeDomain, TrotterDomain}, direction}
    (; transition, gamma_norm_factor, energy_labels) = ws.precomputed_data
    CT = one(T); ZT = zero(T)

    # Select G matrices based on direction (forward vs adjoint)
    G_L = direction === :forward ? ws.G_left : ws.G_left_adj
    G_R = direction === :forward ? ws.G_right : ws.G_right_adj

    # Effective Hamiltonian
    BLAS.gemm!('N', 'N', CT, G_L, rho, ZT, ws.rho_out)
    BLAS.gemm!('N', 'N', CT, rho, G_R, CT, ws.rho_out)

    # Domain-dispatched prefactor (inlined)
    prefactor = compute_rate_prefactor(domain(config), w0(config), t0(config),
                                        sigma(config), gamma_norm_factor)

    # Unified frequency loop
    for (k, eigenbasis) in enumerate(ws.jump_eigenbases)
        foreach_frequency(
            ws.jump_hermitian[k], energy_labels,
            # Positive frequency
            w -> begin
                compute_jump_oft!(ws.jump_oft, eigenbasis, domain(config), ...)
                s = prefactor * transition(w)
                _accumulate_sandwich!(ws.rho_out, ws.jump_oft, rho, s, ws, Val(direction))
            end,
            # Negative frequency (Hermitian partner)
            w -> begin
                s_neg = prefactor * transition(-w)
                neg_dir = direction === :forward ? :forward_adj_L : :adjoint_adj_L
                _accumulate_sandwich!(ws.rho_out, ws.jump_oft, rho, s_neg, ws, Val(neg_dir))
            end
        )
    end
    return ws.rho_out
end

# Public API: thin wrappers
apply_lindbladian!(ws, rho, config, ham) = _apply_lindbladian_impl!(ws, rho, config, ham, Val(:forward))
apply_adjoint_lindbladian!(ws, rho, config, ham) = _apply_lindbladian_impl!(ws, rho, config, ham, Val(:adjoint))
```

### What This Eliminates

| Duplicated Code | Lines Saved | Mechanism |
|---|---|---|
| Prefactor formula (8 locations) | ~40 lines | `compute_rate_prefactor` dispatch |
| OFT computation (6 locations) | ~30 lines | `compute_jump_oft!` dispatch |
| Hermitian half-grid branching (8 locations) | ~80 lines | `foreach_frequency` |
| 4 sandwich variants | ~60 lines | `_accumulate_sandwich!` with Val |
| 6 apply_(adjoint_)lindbladian! | ~350 lines | Unified `_apply_lindbladian_impl!` |
| `_thermalize_to_liouv_config` converters | ~40 lines | Direct S-parameter dispatch |
| Total estimated | ~600 lines | |

---

## 4. File Organization

### Current State

28 source files in flat `src/` directory. The main module file (`QuantumFurnace.jl`) includes them in a dependency-ordered sequence. No subdirectories, no submodules.

### Recommended Pattern: Grouped Includes with Section Comments

**Do NOT use subdirectories or submodules.** Julia packages conventionally keep source files in flat `src/`. Submodules add namespace complexity (`using .SubModule`) and can cause precompilation issues. For an 8K-LOC package, flat `src/` with logical grouping in the include order is the right scale.

**Confidence: HIGH** -- This is standard Julia package convention. KrylovKit.jl (1,800 LOC), DiffEqBase.jl (15K+ LOC), and ITensors.jl all use flat `src/` with grouped includes.

```julia
module QuantumFurnace

# ===== Dependencies =====
using LinearAlgebra, SparseArrays, Random
using Printf, Dates
using BSON, Arpack, FINUFFT, KrylovKit
using LsqFit, QuadGK, Optim, Roots
using SpecialFunctions: erfc
using ProgressMeter, DataStructures
using Base.Threads

# ===== Exports =====
# (grouped by domain, listed here for clarity)

# ===== Type System =====
include("types/domains.jl")         # AbstractDomain, BohrDomain, etc.
include("types/db_flavors.jl")      # KMSDB, GNSDB
include("types/sim_modes.jl")       # Lindbladian, Thermalize, etc.
include("types/config.jl")          # PhysicsParams, SimConfig, constructors, accessors
include("types/results.jl")         # LindbladianResult, DMSimulationResult, etc.

# ===== Core Physics =====
include("physics/constants.jl")
include("physics/hamiltonian.jl")   # HamHam, TrottTrott, Gibbs state
include("physics/transition.jl")    # pick_transition (KMS/GNS dispatch)
include("physics/kossakowski.jl")   # create_alpha, create_alpha_gns (Bohr domain)
include("physics/coherent.jl")      # B operators (Bohr, Time, Trotter)
include("physics/nufft.jl")         # NUFFT prefactor computation
include("physics/precompute.jl")    # _precompute_data, _precompute_labels

# ===== Workspaces =====
include("workspace.jl")             # SimWorkspace, TrajWorkspace, constructors

# ===== Simulation Engines =====
include("engines/oft.jl")           # compute_jump_oft!, frequency iteration
include("engines/sandwich.jl")      # _accumulate_sandwich! (unified forward/adjoint)
include("engines/lindbladian.jl")   # construct_lindbladian, dense L matrix
include("engines/krylov.jl")        # apply_lindbladian!, apply_adjoint_lindbladian!
include("engines/thermalize.jl")    # run_thermalization, DM stepping
include("engines/trajectory.jl")    # TrajectoryFramework, step_along_trajectory!

# ===== Analysis =====
include("analysis/convergence.jl")
include("analysis/fitting.jl")
include("analysis/gap_estimation.jl")
include("analysis/diagnostics.jl")
include("analysis/log_sobolev.jl")

# ===== Utilities =====
include("util/qi_tools.jl")         # trace_distance, fidelity, etc.
include("util/misc_tools.jl")       # Pauli matrices, padding, etc.
include("util/validation.jl")       # validate_config!, _print_press
include("util/results_io.jl")       # save/load ExperimentResult

end # module
```

**Note on subdirectory `include()` paths:** Julia resolves `include()` paths relative to the including file. Since all includes are in `QuantumFurnace.jl` which lives in `src/`, paths like `"types/config.jl"` resolve to `src/types/config.jl`. This is standard and works correctly.

### File Consolidation

Several current files can be merged:

| Current Files | Merged Into | Rationale |
|---|---|---|
| `structs.jl` (358 LOC) | Split into `types/config.jl`, `types/results.jl`, `workspace.jl` | Currently mixes config types, result types, and workspace types |
| `energy_domain.jl` (165), `bohr_domain.jl` (184), `time_domain.jl` (19), `trotter_domain.jl` (206) | `physics/transition.jl` + `physics/kossakowski.jl` | Domain files contain transition functions and alpha computations, not domain definitions |
| `kraus.jl` (14 LOC) | `workspace.jl` | Just the KrausScratch struct definition |
| `errors.jl` (1 LOC) | `util/validation.jl` | Single error type |
| `furnace.jl` + `furnace_utensils.jl` | `engines/lindbladian.jl` + `physics/precompute.jl` | Furnace = Lindbladian construction + precomputation |
| `jump_workers.jl` (461) + `krylov_matvec.jl` (586) | `engines/sandwich.jl` + `engines/krylov.jl` | Shared sandwich logic extracted, domain-specific parts unified |
| `krylov_workspace.jl` (498) + `krylov_eigsolve.jl` (568) | `workspace.jl` + `engines/krylov.jl` | Workspace construction separated from eigensolve logic |
| `ofts.jl` (110) | `engines/oft.jl` | Kept as-is but with unified dispatch |

### File Count Reduction

Current: 28 files, no subdirectories
Proposed: ~22 files in 5 subdirectories (`types/`, `physics/`, `engines/`, `analysis/`, `util/`)

The slight reduction in file count is less important than the logical grouping. A developer looking for "how is the transition function computed?" goes to `physics/transition.jl`, not `energy_domain.jl` (which contains transition functions but also Kossakowski matrices and Bohr dictionary utilities).

---

## 5. Patterns to Avoid

### Anti-Pattern 1: Abstract Type in Hot-Path Fields

**Current problem (already known and partially addressed):** `JumpOp.in_eigenbasis::Matrix{<:Complex}` has an abstract element type parameter, causing boxing allocations when accessed in loops. The existing code already works around this by extracting `jump_eigenbases::Vector{Matrix{CT}}` with concrete `CT` in workspace constructors.

**Rule:** Every field accessed in an inner loop must have a fully concrete type known at compile time. Use parametric struct fields, not abstract type annotations.

```julia
# BAD: abstract element type causes boxing
struct JumpOp{T <: AbstractMatrix{<:Complex}}
    in_eigenbasis::Matrix{<:Complex}   # <:Complex is abstract!
end

# GOOD: concrete element type, parametric
struct JumpOp{T<:Complex}
    in_eigenbasis::Matrix{T}
end
```

### Anti-Pattern 2: Type Parameters that Prevent Specialization

**Trap:** Using `AbstractConfig` as a function argument type instead of the parametric `SimConfig{S,D,DB,T}`.

```julia
# BAD: loses all type information, prevents specialization
function foo(config::AbstractConfig, ...)
    if config.domain isa EnergyDomain  # Runtime check!
        ...
    end
end

# GOOD: dispatch on the type parameters
function foo(config::SimConfig{S, EnergyDomain, DB, T}, ...) where {S, DB, T}
    # Compiler knows domain at compile time
end
```

### Anti-Pattern 3: Generated Functions for This Use Case

`@generated` functions are sometimes suggested for eliminating branching. **Do NOT use them here.** The codebase's branching is already on type parameters (Domain, DB flavor), which the compiler handles through normal dispatch specialization. `@generated` functions add complexity, are harder to debug, and provide no performance benefit when dispatch already eliminates the branches.

**Confidence: HIGH** -- Julia docs warn: "do not use a generated function if you can write normally dispatched function that does the same thing."

### Anti-Pattern 4: Over-Parameterizing the Config Type

Adding a type parameter for every boolean flag (e.g., `Config{S, D, DB, WithCoherent, WithLinearComb, T}`) creates a combinatorial explosion of compiled specializations. Each new type parameter doubles the number of compiled method bodies.

**Rule of thumb:** Parameterize on axes that have DIFFERENT CODE PATHS (Domain, DB, SimMode). Keep boolean flags as runtime values for flags that only affect scalar computations (with_coherent, with_linear_combination). The `with_coherent` flag is checked once per setup to decide whether to compute B; it does not change the hot-path loop structure.

### Anti-Pattern 5: Macro-Heavy Struct Generation

Packages like Mixers.jl or custom `@def_config` macros can auto-generate struct fields from templates. **Avoid this.** The resulting code is hard to navigate (field definitions are hidden behind macros), hard to debug (error messages refer to generated code), and unnecessary when composition solves the duplication problem.

---

## 6. Tools for Verification

### Existing Tools (No New Dependencies)

| Tool | Purpose | Usage |
|---|---|---|
| `@code_warntype` | Verify type stability of refactored functions | `@code_warntype apply_lindbladian!(ws, rho, config, ham)` |
| `@allocated` | Verify zero-allocation hot paths after refactoring | `@allocated step_along_trajectory!(psi, fw, ws, rng)` |
| `Test.@testset` | Regression testing against existing behavior | Compare output matrices before/after refactoring |
| `Aqua.jl` (existing) | Package quality checks (ambiguities, unbound args) | `Aqua.test_all(QuantumFurnace)` -- catches dispatch ambiguities from new type hierarchy |

### Key Verification Steps

1. **Type stability check** after introducing `SimConfig{S,D,DB,T}`:
   ```julia
   config = SimConfig{Lindbladian, EnergyDomain, KMSDB, Float64}(...)
   @code_warntype pick_transition(config)  # Must show concrete return type
   @code_warntype compute_rate_prefactor(domain(config), ...)  # Must show Float64
   ```

2. **Allocation check** on refactored hot paths:
   ```julia
   # Must still report 0 allocations
   @allocated apply_lindbladian!(ws, rho, config, hamiltonian)
   @allocated step_along_trajectory!(psi, fw, ws, rng)
   ```

3. **Numerical regression** against saved reference values:
   ```julia
   # Before refactor: save reference outputs
   BSON.bson("reference.bson", liouv=construct_lindbladian(...), gap=krylov_spectral_gap(...))

   # After refactor: compare
   ref = BSON.load("reference.bson")
   new_liouv = construct_lindbladian(...)
   @test norm(new_liouv - ref[:liouv]) / norm(ref[:liouv]) < 1e-12
   ```

4. **Aqua.jl dispatch ambiguity check**:
   ```julia
   # New type hierarchy must not introduce method ambiguities
   Aqua.test_ambiguities(QuantumFurnace)
   ```

---

## Alternatives Considered

| Category | Recommended | Alternative | Why Not |
|---|---|---|---|
| Config hierarchy | Parametric `SimConfig{S,D,DB,T}` | Separate structs with shared base | Julia has no struct inheritance. Shared abstract base + duplicate fields is the current (broken) pattern. |
| Config hierarchy | Parametric types for S,D,DB | Holy Traits for DB/SimMode | Unnecessary complexity. Type parameters on the config struct propagate naturally through dispatch. Traits add an indirection layer (`DBTrait(config)`) with no benefit when the config already carries the information as a type parameter. |
| Config hierarchy | Type parameters | Val-based dispatch | Val requires explicit wrapping at every call site. Type parameters are automatic and checked at construction time. |
| Workspace consolidation | Union{Nothing, T} fields | Separate workspace per simulation mode | Current approach (5 separate structs) is what we are trying to fix. Union fields are acceptable for setup-time-only fields. |
| Workspace consolidation | Single SimWorkspace | Composition (MatrixScratch embedded) | Composition requires accessor forwarding or getproperty overloading. For this codebase size, a single flat struct with Union fields is simpler and equally performant. |
| Code deduplication | Dispatch on Domain singleton | if/elseif chains on domain | Runtime branching inside hot loops. Dispatch eliminates the branch at compile time. |
| Code deduplication | `foreach_frequency` with closures | Macro-generated loops | Macros are harder to debug and understand. Closures with concrete type parameters are zero-overhead in Julia. |
| File organization | Flat src/ with subdirectories | Julia submodules | Submodules add namespace complexity (using .SubModule), can cause precompilation issues, and are not the Julia convention for single-package codebases under 20K LOC. |
| File organization | Subdirectories in src/ | Fully flat src/ (current) | 28 files in flat directory is navigable but unstructured. 22 files in 5 subdirectories groups by concern without adding module boundaries. |

---

## Summary of Stack Decisions

| Decision | Choice | Rationale |
|---|---|---|
| New dependencies | **None** | Pure refactoring using Julia's type system |
| Config type design | `SimConfig{S, D, DB, T}` with 4 type parameters | Encodes all orthogonal axes; enables zero-cost dispatch on any axis |
| DB flavor encoding | Singleton types (`KMSDB`, `GNSDB`) as type parameter | Compile-time dispatch, extensible for future DLL variant |
| Workspace design | Single `SimWorkspace` with Union{Nothing,T} for optional fields | Consolidates 4 workspace structs; Union fields only accessed at setup time |
| Deduplication mechanism | `@inline` helper functions + domain dispatch | Zero-overhead; compiler inlines and specializes |
| Hermitian loop pattern | `foreach_frequency` with closure parameters | Eliminates 8 copies of half-grid iteration; closures are zero-cost in Julia |
| File organization | Subdirectories in flat src/ (`types/`, `physics/`, `engines/`, `analysis/`, `util/`) | Logical grouping without module boundaries |
| Backward compatibility | Type aliases (`const LiouvConfig = SimConfig{Lindbladian,...}`) | Existing tests and scripts continue to work during transition |

---

## Sources

- [Julia Performance Tips](https://docs.julialang.org/en/v1/manual/performance-tips/) -- Type stability, avoiding abstract containers, function barriers. HIGH confidence.
- [Julia Types Documentation](https://docs.julialang.org/en/v1/manual/types/) -- Parametric types, abstract type hierarchy, type parameter constraints. HIGH confidence.
- [Julia Methods Documentation](https://docs.julialang.org/en/v1/manual/methods/) -- Multiple dispatch, parametric methods, method specialization. HIGH confidence.
- [Julia isbits Union Optimization](https://docs.julialang.org/en/v1/devdocs/isbitsunionarrays/) -- Union{Nothing, T} stored inline for isbits T. HIGH confidence.
- [Performance of Fields with Union Type](https://discourse.julialang.org/t/performance-of-fields-with-union-type/88863) -- Explicit unions enable union-splitting; abstract types do not. HIGH confidence.
- [Types vs Traits for Dispatch](https://discourse.julialang.org/t/types-vs-traits-for-dispatch/46296) -- Type parameters sufficient for most cases; traits for orthogonal extension. MEDIUM confidence (community discussion).
- [Composition and Inheritance: The Julian Way](https://discourse.julialang.org/t/composition-and-inheritance-the-julian-way/11231) -- Composition over inheritance; accessor functions for interface. MEDIUM confidence (community discussion).
- [Holy Traits Pattern](https://ahsmart.com/pub/holy-traits-design-patterns-and-best-practice-book/) -- Tim Holy Trait Trick: zero-cost compile-time dispatch via singleton types. HIGH confidence.
- [Multiple Dispatch Designs](http://ucidatascienceinitiative.github.io/IntroToJulia/Html/DispatchDesigns) -- Duck typing, hierarchies, traits comparison. MEDIUM confidence.
- [ValSplit.jl](https://github.com/ztangent/ValSplit.jl) -- Val dispatch compiles to zero-overhead switch. HIGH confidence (benchmark data in repo).
- [DifferentialEquations.jl Architecture](https://openresearchsoftware.metajnl.com/articles/jors.151) -- Algorithm types as parametric dispatch targets; common solve interface. HIGH confidence (peer-reviewed paper).
- [Best Practices for Structuring Larger Projects](https://discourse.julialang.org/t/best-practices-for-structuring-larger-projects/2652) -- Flat src/ with grouped includes is standard for single packages. MEDIUM confidence (community discussion).
- [Julia Modules Documentation](https://docs.julialang.org/en/v1/manual/modules/) -- Submodule semantics, include path resolution. HIGH confidence.

---

*Stack research for: QuantumFurnace.jl codebase restructure*
*Researched: 2026-02-25*
