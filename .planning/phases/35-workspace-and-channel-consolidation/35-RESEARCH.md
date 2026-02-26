# Phase 35: Workspace and Channel Consolidation - Research

**Researched:** 2026-02-26
**Domain:** Julia struct consolidation, parametric dispatch, quantum simulation workspaces
**Confidence:** HIGH

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

#### Unified Workspace Struct Design
- Parameterized `Workspace{S,D,C,T}` matching Config's type parameters
- S expanded to 4 simulation singletons: `Lindbladian`, `Thermalize`, `Krylov`, `Trajectory`
- Each Workspace{S,...} has different fields appropriate to its simulation path
- Nested structure: precomputed/immutable physics data as flat fields, mutable scratch buffers bundled into a single nested `Scratch` struct field
- `Workspace{Krylov,...}` handles both Krylov-Lindbladian (G_left/G_right) and Krylov-Thermalize (CPTP channel) modes as a single type, using `Union{Nothing,...}` for mode-specific fields
- Constructor: `Workspace(config::Config{S,D,C,T}, hamiltonian, jumps)` -- type constructor, dispatch on Config
- **Critical constraint:** Only ONE constructor per workspace variant (inner OR outer, never both) -- multiple constructors caused bugs previously
- Domain-specific precomputed data (NUFFT for TimeDomain/TrotterDomain) absorbed into the workspace parameterization, replacing the generic `precomputed_data::NamedTuple`

#### CPTP Channel Computation API
- `_build_cptp_channel(R, delta)` stays as a single pure function -- callers pass R_total (summed) for Krylov or R^a (per-jump) for DM/Trajectory
- Function doesn't know or care whether input R is summed or per-jump
- Coherent B computation stays separate from `_build_cptp_channel` -- clean physics separation
- Eliminate single-jump variants of `B_time`, `B_trotter`, etc. -- keep only the `jumps` (vector) variant; callers wrap single jump as `[jump]`
- Same deduplication for `_precompute_coherent_B` and `_precompute_coherent_unitary` -- vector-of-jumps only
- K0, U_residual, alpha stored as direct flat fields on the workspace (not wrapped in NamedTuple or sub-struct)
- G_left, G_right, G_left_adj, G_right_adj as direct workspace fields with `Union{Nothing,...}` for non-Krylov-Lindbladian workspaces

#### TrajectoryWorkspace Integration
- `Workspace{Trajectory,D,C,T}` follows same pattern as other simulation types
- PerOperatorKraus data absorbed: per-jump R^a, K0^a, U_residual^a, U_B^a stored as vectors of matrices (e.g., `K0s::Vector{Matrix{T}}`) -- flattened, not nested in PerOperatorKraus
- Only mutable scratch buffers (jump_oft, psi_tmp, Rpsi, rho_acc) go into the nested Scratch struct

#### Field Naming
- Aim for physics-descriptive names for scratch buffers (not generic tmp1/tmp2)
- No `channel_` prefix on CPTP fields: `K0`, `U_residual` directly on workspace
- Claude to audit dead/redundant fields across all three workspaces during research (KrausScratch.K0 already known dead from Phase 34-02)

### Claude's Discretion
- Exact Scratch struct field names (descriptive but Claude determines the best names per simulation path)
- Whether domain-specific precomputed data gets its own type parameter or is handled via dispatch in constructors
- How to handle the Identity matrix (currently a field in LindbladianWorkspace -- may become computed or shared)
- Optimal number of scratch matrices per simulation type

### Deferred Ideas (OUT OF SCOPE)
- Adding `Krylov` and `Trajectory` as new Config simulation singleton types may overlap with Phase 36 (API and Results). The workspace phase defines the singletons; Phase 36 wires them to `run_*` entry points.
</user_constraints>

## Summary

This phase consolidates four separate workspace types (`KrylovWorkspace`, `KrausScratch`, `LindbladianWorkspace`, `TrajectoryWorkspace` + `TrajectoryFramework` + `PerOperatorKraus`) into a unified parametric `Workspace{S,D,C,T}` struct hierarchy. The current codebase has significant field duplication and naming inconsistency across workspaces (e.g., KrylovWorkspace has `channel_K0`/`channel_U_residual`/`channel_rho_jump` while KrausScratch has `K0`/`rho_jump`; TrajectoryFramework duplicates hot-path fields from precomputed_data). A unified design eliminates this redundancy while preserving the zero-allocation hot-path guarantees.

The research reveals that the consolidation is primarily a structural refactoring with well-defined boundaries. The key risk is breaking the zero-allocation invariant in the matvec hot path (currently verified by 8 separate `@allocated == 0` tests across 4 domains x 2 directions). The main complexity lies in the fact that `Workspace{Krylov,...}` must serve two distinct roles (Lindbladian matvec via G_left/G_right, and Thermalize channel via K0/U_residual/U_coherent), while `Workspace{Trajectory,...}` absorbs the per-operator data that currently lives in `TrajectoryFramework` + `PerOperatorKraus`.

**Primary recommendation:** Execute in 2 plans: Plan 01 consolidates KrylovWorkspace + KrausScratch + LindbladianWorkspace into `Workspace{Lindbladian|Thermalize|Krylov,...}`, and Plan 02 consolidates TrajectoryWorkspace + TrajectoryFramework + PerOperatorKraus into `Workspace{Trajectory,...}`. Both plans preserve zero-allocation invariants and run the full test suite.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Julia | 1.x | Language | Project language |
| LinearAlgebra | stdlib | BLAS, eigen, etc. | Hot-path matrix operations |
| KrylovKit | pkg | Arnoldi eigsolve | Matrix-free spectral gap |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| FINUFFT | pkg | NUFFT prefactors | TimeDomain/TrotterDomain OFT |
| Test | stdlib | @test, @allocated | Allocation regression tests |

No new libraries needed. This is purely a structural refactoring of existing code.

## Architecture Patterns

### Current Workspace Landscape (What Exists)

```
src/structs.jl:
  LindbladianWorkspace{T}    -- 5 fields (Id, 4 scratch matrices)
                             -- Used ONLY by construct_lindbladian (dense Liouvillian)

src/kraus.jl:
  KrausScratch{T}            -- 8 fields (jump_oft, LdagL, R, rho_jump, K0*, tmp1, tmp2, rho_next)
                             -- Used by run_thermalization (DM evolution)
                             -- K0 field is DEAD (confirmed Phase 34-02)
                             -- Also used as builder scratch for TrajectoryFramework

src/krylov_workspace.jl:
  KrylovWorkspace{T,PD}      -- 21 fields total:
                             -- precomputed_data::PD (domain-specific NamedTuple)
                             -- B_total, jumps, jump_eigenbases, jump_hermitian
                             -- jump_oft, tmp1, tmp2, LdagL, rho_out (5 scratch)
                             -- channel_K0, channel_U_residual, channel_U_coherent,
                                channel_rho_jump, channel_delta (5 channel, Nothing for Lindbladian)
                             -- G_left, G_right, G_left_adj, G_right_adj (4 effective Hamiltonian)
                             -- 2 constructors: Config{Lindbladian}, Config{Thermalize}

src/trajectories.jl:
  TrajectoryWorkspace{T}     -- 4 fields (jump_oft, psi_tmp, Rpsi, rho_acc)
  PerOperatorKraus{T}        -- 4 fields (R, K0, U_residual, U_B)
  TrajectoryFramework{T,D,F,P} -- 14 fields (domain, jumps, ham_or_trott, config,
                                precomputed_data, per_operator, n_jumps, delta, delta_eff,
                                alpha, scaled_prefactor, sigma, transition, energy_labels,
                                oft_nufft_prefactors)
```

### Target Architecture

```
Workspace{S<:AbstractSimulation, D<:AbstractDomain, C<:AbstractConstruction, T<:AbstractFloat}
  where S in {Lindbladian, Thermalize, Krylov, Trajectory}

Each variant has:
  - Flat immutable physics fields (precomputed at construction)
  - Nested mutable Scratch struct (simulation-path-specific buffers)

Workspace{Lindbladian, D, C, T}:
  -- Used by construct_lindbladian (dense Liouvillian)
  -- Fields: Id, scratch::LiouvillianScratch{T}

Workspace{Thermalize, D, C, T}:
  -- Used by run_thermalization (DM Kraus steps)
  -- Fields: coherent_unitaries, scratch::ThermalizeScratch{T}

Workspace{Krylov, D, C, T}:
  -- Used by apply_lindbladian!, apply_delta_channel!, krylov_spectral_gap
  -- Fields: jump_eigenbases, jump_hermitian, B_total,
             G_left, G_right, G_left_adj, G_right_adj (Union{Nothing,...}),
             K0, U_residual, U_coherent (Union{Nothing,...}),
             delta (Union{Nothing,...}),
             [domain-specific precomputed data: transition, gamma_norm_factor,
              energy_labels, oft_domain_prefactor, alpha, ...]
             scratch::KrylovScratch{T}

Workspace{Trajectory, D, C, T}:
  -- Used by step_along_trajectory!, run_trajectories
  -- Fields: jumps, n_jumps, delta, alpha, scaled_prefactor, sigma,
             transition, energy_labels, oft_nufft_prefactors,
             K0s, U_residuals, Rs, U_Bs (vectors of matrices, flattened from PerOperatorKraus),
             scratch::TrajectoryScratch{T}
```

### Pattern 1: Parametric Struct with Simulation-Type Dispatch

**What:** Julia parametric struct where the first type parameter `S` determines which fields exist, using `Union{Nothing,...}` for optional fields.

**When to use:** When a single struct serves multiple roles with overlapping but not identical field sets.

**Current example (KrylovWorkspace):**
```julia
# Already uses this pattern for channel fields:
struct KrylovWorkspace{T<:Complex, PD<:NamedTuple}
    # ... common fields ...
    channel_K0::Union{Nothing, Matrix{T}}  # populated only for Thermalize
    # ... etc ...
end
```

**Target pattern:**
```julia
struct Workspace{S<:AbstractSimulation, D<:AbstractDomain, C<:AbstractConstruction, T<:AbstractFloat}
    # Fields vary by S -- Union{Nothing,...} for optional
    # Only ONE outer constructor per S variant
end

# Single dispatch constructor:
function Workspace(config::Config{Lindbladian,D,C,T}, hamiltonian, jumps; trotter=nothing) where {D,C,T}
    # Build Workspace{Lindbladian,D,C,T}
end
```

### Pattern 2: Nested Scratch Structs

**What:** Separate mutable scratch buffers from immutable precomputed data by nesting mutable buffers in a dedicated struct.

**Why:** Clean separation makes it obvious which data is read-only (physics) vs mutable (per-call scratch). Enables independent workspace instances for thread safety.

**Example:**
```julia
struct KrylovScratch{T<:Complex}
    jump_oft::Matrix{T}
    sandwich_tmp::Matrix{T}    # was tmp1
    sandwich_out::Matrix{T}    # was LdagL
    rho_out::Matrix{T}
end

struct TrajectoryScratch{T<:Complex}
    jump_oft::Matrix{T}
    psi_tmp::Vector{T}
    Rpsi::Vector{T}
    rho_acc::Matrix{T}
end
```

### Pattern 3: Vector-of-Jumps Only API for Coherent B Functions

**What:** Eliminate single-jump overloads of B_time, B_trotter, B_bohr. Callers wrap single jump as `[jump]`.

**Current state:** 6 functions (2 each for B_time, B_trotter, B_bohr -- single-jump and vector variants).

**Target:** 3 functions (vector-of-jumps only). Callers do `B_time([jump], ...)` instead of `B_time(jump, ...)`.

**Call sites that use single-jump:**
- `_precompute_coherent_unitary` in coherent.jl (lines 76, 85, 94): iterates per-jump, calls single-jump B
- `build_trajectoryframework` in trajectories.jl (line 141): wraps as `JumpOp[jumps[a]]`

### Anti-Patterns to Avoid
- **Multiple constructors for same type variant:** The user explicitly flagged this as a past source of bugs. Each workspace variant (Lindbladian, Thermalize, Krylov, Trajectory) must have exactly ONE constructor.
- **Storing mutable and immutable data at same nesting level:** Makes it unclear what can be mutated. Use Scratch sub-struct.
- **Generic field names (tmp1, tmp2):** Prefer physics-descriptive names that indicate purpose.

## Field Audit: Dead and Redundant Fields

### Confirmed Dead Fields
| Field | Location | Evidence | Action |
|-------|----------|----------|--------|
| `KrausScratch.K0` | kraus.jl:5 | Never written after Phase 34-02 extracted `_build_cptp_channel`. `_finalize_kraus_step!` no longer writes to `scratch.K0`. | Remove |

### Potentially Redundant Fields
| Field | Location | Analysis | Recommendation |
|-------|----------|----------|----------------|
| `LindbladianWorkspace.Id` | structs.jl:22 | Only used by `_vectorize_liouvillian_coherent!` and `_vectorize_liouv_diss_and_add!` in the dense Liouvillian construction path. Could be computed inline as `Matrix{CT}(I, dim, dim)` -- it's a construction-time-only cost, not hot path. | **Compute inline** -- the dense Liouvillian path allocates `dim^2 x dim^2` anyway, so one `dim x dim` identity is negligible. Remove from workspace. |
| `KrylovWorkspace.jumps` | krylov_workspace.jl:48 | Stored "for external access" per docstring. Only used by callers via `ws.jumps`. Hot path uses `ws.jump_eigenbases` instead. With unified workspace, jumps could be a field on all variants. | **Keep on Krylov workspace** for external access pattern |
| `TrajectoryFramework.config` | trajectories.jl:83 | Stored but not accessed in hot path (all hot-path fields are extracted to concrete-typed flat fields). Only used at construction time. | **Remove from workspace** -- callers already have config |
| `TrajectoryFramework.ham_or_trott` | trajectories.jl:82 | Used by EnergyDomain step_along_trajectory! to access `ham.bohr_freqs`. | **Keep** -- needed for EnergyDomain oft! calls |
| `TrajectoryFramework.precomputed_data::Any` | trajectories.jl:84 | Abstract-typed, not accessed in hot path. All needed fields extracted to concrete-typed flat fields. | **Remove** -- all needed data already extracted |
| `TrajectoryFramework.delta_eff` | trajectories.jl:92 | Always equals `delta` (comment says "delta_eff field = delta"). | **Remove** -- use `delta` directly |
| `KrylovWorkspace.channel_delta` | krylov_workspace.jl:66 | Only used in `apply_delta_channel!` where it's read as `ws.channel_delta`. With unified naming, becomes just `delta`. | **Rename to `delta`** on Krylov workspace |

### Fields Per Simulation Path (After Consolidation)

#### Workspace{Lindbladian,...} -- Dense Liouvillian Construction
Immutable:
- (none needed -- Id computed inline)

Scratch:
- `jump_tmp::Matrix{CT}` -- OFT output / alpha-weighted jump
- `jump_conj::Matrix{CT}` -- elementwise conjugate scratch
- `jump_dag_jump::Matrix{CT}` -- L'*L product
- `jump2_jump1::Matrix{CT}` -- Bohr domain mixed product

#### Workspace{Thermalize,...} -- DM Kraus Steps
Immutable:
- `coherent_unitaries::Union{Nothing, Vector{Matrix{CT}}}` -- per-jump exp(-i*delta*B_k)

Scratch:
- `jump_oft::Matrix{CT}` -- OFT output buffer
- `LdagL::Matrix{CT}` -- L'*L accumulation
- `R::Matrix{CT}` -- R accumulator (zeroed each step)
- `rho_jump::Matrix{CT}` -- jump sandwich accumulator
- `sandwich_tmp::Matrix{CT}` -- was tmp1
- `rho_work::Matrix{CT}` -- was tmp2
- `rho_next::Matrix{CT}` -- next-step result

(Note: dead `K0` field removed; total reduced from 8 to 7 scratch fields)

#### Workspace{Krylov,...} -- Matrix-Free Lindbladian/Channel
Immutable (physics):
- `jump_eigenbases::Vector{Matrix{CT}}` -- concrete-typed eigenbasis matrices
- `jump_hermitian::Vector{Bool}` -- hermitian flags
- `B_total::Union{Nothing, Matrix{CT}}` -- coherent B (Nothing for GNS)
- `G_left::Union{Nothing, Matrix{CT}}` -- effective Hamiltonian left
- `G_right::Union{Nothing, Matrix{CT}}` -- effective Hamiltonian right
- `G_left_adj::Union{Nothing, Matrix{CT}}` -- adjoint left
- `G_right_adj::Union{Nothing, Matrix{CT}}` -- adjoint right
- `K0::Union{Nothing, Matrix{CT}}` -- CPTP no-event operator (Thermalize only)
- `U_residual::Union{Nothing, Matrix{CT}}` -- CPTP residual (Thermalize only)
- `U_coherent::Union{Nothing, Matrix{CT}}` -- coherent unitary (Thermalize+coherent only)
- `delta::Union{Nothing, Float64}` -- step size (Thermalize only)

Immutable (domain-specific precomputed data -- absorbed from NamedTuple):
- `transition::F` -- transition function (concrete closure type)
- `gamma_norm_factor::Float64`
- `energy_labels::Vector{Float64}`
- `oft_domain_prefactor::Float64`
- `alpha::Union{Nothing, Function}` -- BohrDomain alpha function
- `oft_nufft_prefactors::Union{Nothing, P}` -- Time/TrotterDomain NUFFT data
- `bohr_keys::Union{Nothing, Vector{...}}` -- BohrDomain Thermalize
- `bohr_is::Union{Nothing, Vector{Vector{Int}}}` -- BohrDomain Thermalize
- `bohr_js::Union{Nothing, Vector{Vector{Int}}}` -- BohrDomain Thermalize
- `b_minus::Union{Nothing, Dict{...}}` -- Time/Trotter coherent
- `b_plus::Union{Nothing, Dict{...}}` -- Time/Trotter coherent

Scratch:
- `jump_oft::Matrix{CT}` -- OFT output
- `sandwich_tmp::Matrix{CT}` -- was tmp1 (BLAS gemm scratch)
- `sandwich_out::Matrix{CT}` -- was LdagL (sandwich result)
- `rho_out::Matrix{CT}` -- output accumulator
- `channel_rho_jump::Union{Nothing, Matrix{CT}}` -- channel jump sandwich (Thermalize only)

Also stored (for external access, not hot path):
- `jumps::Vector{JumpOp}` -- reference to jump operators

#### Workspace{Trajectory,...} -- State-Vector Trajectory Simulation
Immutable (physics):
- `jumps::Vector{JumpOp{Matrix{CT}}}` -- concrete-typed jumps
- `ham_or_trott::Union{HamHam, TrottTrott}`
- `n_jumps::Int`
- `delta::Float64`
- `alpha::Float64`
- `scaled_prefactor::Float64`
- `sigma::Float64`
- `transition::F` -- concrete closure
- `energy_labels::Vector{Float64}`
- `oft_nufft_prefactors::P` -- NUFFT data (Nothing for EnergyDomain)
- `Rs::Vector{Matrix{CT}}` -- per-jump R^a (flattened from PerOperatorKraus)
- `K0s::Vector{Matrix{CT}}` -- per-jump K0^a
- `U_residuals::Vector{Matrix{CT}}` -- per-jump U_residual^a
- `U_Bs::Vector{Union{Nothing, Matrix{CT}}}` -- per-jump coherent unitary

Scratch:
- `jump_oft::Matrix{CT}`
- `psi_tmp::Vector{CT}`
- `Rpsi::Vector{CT}`
- `rho_acc::Matrix{CT}`

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| CPTP channel matrices | Inline K0/U_residual formulas | `_build_cptp_channel(R, delta)` | Already extracted in Phase 34-02; single source of truth |
| Coherent B computation | Per-simulation-type B calls | `_precompute_coherent_B(jumps, ...)` | Already the canonical API |
| PSD square root | Manual eigendecomposition | The one in `_build_cptp_channel` | Handles clamping of negative eigenvalues |

**Key insight:** The main consolidation work is structural (moving fields between structs, renaming, eliminating dead code). The physics computation functions (`_build_cptp_channel`, `_precompute_coherent_B`, `_accumulate_R_total!`, etc.) are already in good shape from Phase 34.

## Common Pitfalls

### Pitfall 1: Breaking Zero-Allocation Hot Path
**What goes wrong:** Changing field types or struct layout introduces type instability, causing runtime allocations in the matvec/step loops.
**Why it happens:** Julia's type inference is sensitive to `Union{Nothing, ...}` fields and abstract type parameters. If the compiler cannot narrow the type at a hot-path access site, it boxes the value.
**How to avoid:**
1. Keep all hot-path fields concrete-typed (no `Union` in hot path access)
2. For `Union{Nothing, Matrix{T}}` fields: the hot path must branch on `=== nothing` BEFORE accessing the matrix, so the compiler narrows the type within the branch.
3. Run `@allocated == 0` tests after every structural change.
**Warning signs:** `@allocated` tests fail with small (16-176 byte) allocations -- indicates boxing.

### Pitfall 2: Accidentally Aliased Scratch Buffers
**What goes wrong:** Two different operations share the same scratch buffer within a single call, leading to data corruption.
**Why it happens:** When consolidating scratch fields, it's tempting to reduce buffer count. But within a single matvec call, buffers like `sandwich_tmp` and `sandwich_out` must not alias because they're used simultaneously by BLAS.gemm!.
**How to avoid:** Trace every BLAS call in the hot path to verify no two concurrent reads/writes share a buffer.
**Warning signs:** Non-deterministic numerical errors, results that differ by small amounts between runs.

### Pitfall 3: Multiple Constructors for Same Variant
**What goes wrong:** Two constructors compute overlapping fields differently, leading to subtle initialization bugs.
**Why it happens:** Natural temptation to have both `Workspace{Krylov}(Config{Lindbladian})` and `Workspace{Krylov}(Config{Thermalize})` as inner constructors.
**How to avoid:** Hard constraint from user: only ONE constructor per workspace variant. Use a single outer constructor with dispatch on Config's simulation type parameter.
**Warning signs:** Fields being `nothing` when they shouldn't be, or vice versa.

### Pitfall 4: Domain-Specific Precomputed Data Type Instability
**What goes wrong:** Absorbing the `precomputed_data::NamedTuple` into flat fields creates type instability for domain-varying fields (e.g., `alpha` function exists only for BohrDomain, `oft_nufft_prefactors` exists only for Time/Trotter).
**Why it happens:** `Union{Nothing, ...}` on frequently-accessed fields can hurt performance.
**How to avoid:** Hot-path code already branches on domain type (`D` parameter) before accessing domain-specific fields. The compiler specializes on `D` and eliminates dead branches. Keep this pattern.
**Warning signs:** Type inference warnings, allocation in domain-specific hot paths.

### Pitfall 5: Forgetting to Update Exports and Test References
**What goes wrong:** Renaming `KrylovWorkspace` to `Workspace` breaks downstream imports and test files.
**Why it happens:** Many test files explicitly reference `KrylovWorkspace`, `TrajectoryWorkspace`, `KrausScratch`, `build_trajectoryframework` etc.
**How to avoid:** Search all `export` statements in `QuantumFurnace.jl`, all `using QuantumFurnace:` in tests, and all direct type references. Update systematically.
**Warning signs:** `UndefVarError` at test time.

## Code Examples

### Example 1: Unified Workspace Struct Definition

```julia
# Simulation singletons (add Krylov to existing set in structs.jl)
struct Krylov <: AbstractSimulation end

# Scratch sub-structs (one per simulation path)
struct KrylovScratch{T<:Complex}
    jump_oft::Matrix{T}       # OFT output buffer
    sandwich_tmp::Matrix{T}   # BLAS gemm scratch (was tmp1)
    sandwich_out::Matrix{T}   # sandwich result (was LdagL)
    rho_out::Matrix{T}        # output accumulator
    channel_rho_jump::Union{Nothing, Matrix{T}}  # Thermalize-only jump scratch
end

struct Workspace{S<:AbstractSimulation, D<:AbstractDomain, C<:AbstractConstruction, T<:AbstractFloat}
    # Fields determined by S (details vary per variant)
    # ... see field listings above ...
end
```

### Example 2: Single Outer Constructor with Config Dispatch

```julia
# Krylov workspace -- single constructor, dispatches on config.sim
function Workspace(
    config::Config{S,D,C,T},
    hamiltonian::HamHam,
    jumps::Vector{JumpOp};
    trotter::Union{TrottTrott, Nothing}=nothing,
) where {S<:Union{Lindbladian,Thermalize}, D, C, T}
    # ... common setup ...
    ham_or_trott = _resolve_ham_or_trott(config, hamiltonian, trotter)
    precomputed_data = _precompute_data(config, ham_or_trott)
    B_total = _precompute_coherent_B(jumps, ham_or_trott, config, precomputed_data)
    # ... build all fields ...
    # Key: Lindbladian sets channel fields to nothing; Thermalize populates them
end
```

### Example 3: Eliminating Single-Jump B Variants

```julia
# BEFORE (coherent.jl has both):
B_time(jump::JumpOp, ...)   # single-jump
B_time(jumps::Vector{JumpOp}, ...)  # multi-jump

# AFTER: only vector variant remains
B_time(jumps::Vector{JumpOp}, ...)

# Callers that had single jump wrap it:
# _precompute_coherent_unitary (coherent.jl):
B_time([jump], hamiltonian, b_minus, b_plus, t0, beta, sigma)
```

### Example 4: Flattening PerOperatorKraus into Workspace{Trajectory}

```julia
# BEFORE:
per_operator::Vector{PerOperatorKraus{T}}  # nested struct
per_operator[a].K0  # access pattern

# AFTER:
K0s::Vector{Matrix{T}}           # flat vector of matrices
U_residuals::Vector{Matrix{T}}
Rs::Vector{Matrix{T}}
U_Bs::Vector{Union{Nothing, Matrix{T}}}

K0s[a]  # access pattern (same indexing, no struct indirection)
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| 4 separate workspace types | Unified `Workspace{S,D,C,T}` | Phase 35 (this phase) | Single struct, clear dispatch |
| `precomputed_data::NamedTuple` generic | Flat typed fields on workspace | Phase 35 (this phase) | Better type stability, clearer API |
| Single + vector B variants | Vector-only B functions | Phase 35 (this phase) | Reduced code duplication |
| `channel_` prefixed fields | Direct physics names | Phase 35 (this phase) | Cleaner naming |

**Deprecated/outdated:**
- `KrylovWorkspace{T,PD}`: Replaced by `Workspace{Krylov,D,C,T}`
- `KrausScratch{T}`: Replaced by `ThermalizeScratch{T}` (nested in Workspace{Thermalize})
- `LindbladianWorkspace{T}`: Replaced by `Workspace{Lindbladian,D,C,T}` with `LiouvillianScratch`
- `TrajectoryWorkspace{T}`: Replaced by scratch sub-struct of `Workspace{Trajectory}`
- `TrajectoryFramework{T,D,F,P}`: Absorbed into `Workspace{Trajectory,D,C,T}`
- `PerOperatorKraus{T}`: Flattened into vectors on `Workspace{Trajectory}`

## Recommendations for Claude's Discretion Items

### Domain-Specific Precomputed Data: Type Parameter vs Dispatch

**Recommendation: Use dispatch in constructors, NOT an additional type parameter.**

Rationale: The current `PD<:NamedTuple` type parameter on KrylovWorkspace is the source of 176-byte boxing allocations that required workaround function barriers in tests. Absorbing the precomputed data fields directly into the workspace eliminates this. The `D` type parameter from Config already provides domain dispatch -- constructors can branch on `D` to populate domain-specific fields. Fields not relevant to a domain are set to `nothing`.

The alternative (adding a `PD` type parameter) would mean `Workspace{S,D,C,T,PD}` with 5 type parameters, which is unwieldy and provides no hot-path benefit since the compiler already specializes on `D`.

### Identity Matrix Handling

**Recommendation: Compute inline at construction time, do not store as field.**

The `Id` matrix in `LindbladianWorkspace` is only used by dense Liouvillian vectorization helpers (`_vectorize_liouvillian_coherent!`, `_vectorize_liouv_diss_and_add!`). These are called during `construct_lindbladian`, which already allocates a `dim^2 x dim^2` output matrix. Adding one `dim x dim` identity allocation is negligible. Compute it at the call site with `Matrix{CT}(I, dim, dim)`.

### Scratch Field Names (Per Simulation Path)

**Recommended names:**

| Simulation | Field | Description | Replaces |
|------------|-------|-------------|----------|
| Krylov | `jump_oft` | OFT-filtered jump operator | `jump_oft` (unchanged) |
| Krylov | `sandwich_tmp` | BLAS intermediate for sandwich | `tmp1` |
| Krylov | `sandwich_out` | Sandwich product result | `LdagL` |
| Krylov | `rho_out` | Output density matrix | `rho_out` (unchanged) |
| Krylov | `channel_rho_jump` | Channel jump accumulator | `channel_rho_jump` (keep for clarity) |
| Thermalize | `jump_oft` | OFT-filtered jump | `scratch.jump_oft` (unchanged) |
| Thermalize | `LdagL` | L'*L accumulation | Keep as-is (physics name) |
| Thermalize | `R` | R accumulator | Keep as-is |
| Thermalize | `rho_jump` | Jump sandwich accumulator | Keep as-is |
| Thermalize | `sandwich_tmp` | mul! intermediate | `tmp1` |
| Thermalize | `rho_work` | Working density matrix | `tmp2` |
| Thermalize | `rho_next` | Next-step result | Keep as-is |
| Trajectory | `jump_oft` | OFT buffer | unchanged |
| Trajectory | `psi_tmp` | State vector scratch | unchanged |
| Trajectory | `Rpsi` | R*psi cache | unchanged |
| Trajectory | `rho_acc` | Density matrix accumulator | unchanged |
| Lindbladian | `jump_tmp` | Jump operator scratch | unchanged |
| Lindbladian | `jump_conj` | Conjugate scratch | unchanged |
| Lindbladian | `jump_dag_jump` | L'*L product | unchanged |
| Lindbladian | `jump2_jump1` | Bohr mixed product | unchanged |

### Optimal Number of Scratch Matrices

| Simulation Path | Current | After Consolidation | Change |
|----------------|---------|---------------------|--------|
| Lindbladian | 5 (incl Id) | 4 (Id removed) | -1 |
| Thermalize | 8 (incl dead K0) | 7 | -1 (dead K0 removed) |
| Krylov | 5 + 1 channel | 5 (channel_rho_jump moves to Union on scratch) | 0 net |
| Trajectory | 4 | 4 | 0 |

## Impact Surface Analysis

### Source Files Requiring Modification

| File | Current Types Used | Changes Needed |
|------|-------------------|----------------|
| `src/structs.jl` | `LindbladianWorkspace` definition | Remove `LindbladianWorkspace`, add `Krylov` singleton |
| `src/kraus.jl` | `KrausScratch` definition | Replace with `ThermalizeScratch` or remove file |
| `src/krylov_workspace.jl` | `KrylovWorkspace` definition + 2 constructors | Replace with `Workspace{Krylov}` + merged constructor |
| `src/krylov_matvec.jl` | `KrylovWorkspace` in all matvec signatures | Update to `Workspace{Krylov}` |
| `src/krylov_eigsolve.jl` | `KrylovWorkspace` construction | Update constructor calls |
| `src/jump_workers.jl` | `LindbladianWorkspace`, `KrausScratch` in signatures | Update to new workspace types |
| `src/furnace.jl` | `LindbladianWorkspace`, `KrausScratch` construction | Update to new workspace constructors |
| `src/furnace_utensils.jl` | `_precompute_data` (unchanged), `_build_cptp_channel` (unchanged) | Possibly no changes |
| `src/trajectories.jl` | `TrajectoryWorkspace`, `TrajectoryFramework`, `PerOperatorKraus` | Major restructure into `Workspace{Trajectory}` |
| `src/coherent.jl` | `B_time(jump)`, `B_trotter(jump)` single-jump variants | Remove single-jump overloads |
| `src/bohr_domain.jl` | `B_bohr(ham, jump)` single-jump variant | Remove single-jump overload |
| `src/QuantumFurnace.jl` | All exports | Update exported names |

### Test Files Requiring Modification

| File | Types Referenced | Changes Needed |
|------|-----------------|----------------|
| `test/test_krylov_matvec.jl` | `KrylovWorkspace` (26 refs) | Update all refs |
| `test/test_allocation.jl` | `KrausScratch`, `TrajectoryWorkspace` | Update refs + verify allocations |
| `test/test_workspace_independence.jl` | `TrajectoryWorkspace`, `TrajectoryFramework` | Update refs |
| `test/test_compilation.jl` | `KrausScratch`, `TrajectoryFramework` | Update refs |
| `test/test_cptp.jl` | `KrausScratch`, `build_trajectoryframework` | Update refs |
| `test/test_krylov_eigsolve.jl` | `KrylovWorkspace` | Update refs |
| `test/test_threading.jl` | `TrajectoryWorkspace` | Update refs |
| `test/test_regression.jl` | `KrausScratch` | Update refs |
| `test/test_trajectory_fixes.jl` | `TrajectoryWorkspace`, `TrajectoryFramework` | Update refs |
| `test/test_gns_trajectory.jl` | `TrajectoryFramework` | Update refs |

## Critical Success Criteria Verification Plan

1. **Struct consolidation:** `KrylovWorkspace`, `KrausScratch`, `LindbladianWorkspace` no longer exist in source. `Workspace{S,D,C,T}` replaces all.
2. **CPTP semantics:** `_build_cptp_channel(R, delta)` unchanged. Per-jump R^a used by DM/Trajectory, summed R_total by Krylov. Verified by existing CPTP tests.
3. **Zero-allocation hot path:** `@allocated == 0` for:
   - `apply_lindbladian!` (Energy, Time, Trotter, Bohr domains)
   - `apply_adjoint_lindbladian!` (same 4 domains)
   - `step_along_trajectory!` (Time/Trotter domains)
4. **Workspace independence:** Independent workspace instances share no mutable state.
5. **All 1198+ tests pass** with identical numerical results.

## Open Questions

1. **Workspace{Krylov} for BohrDomain Thermalize**
   - What we know: BohrDomain Thermalize uses extra precomputed fields (`bohr_keys`, `bohr_is`, `bohr_js`) not present in Lindbladian BohrDomain precomputed_data.
   - What's unclear: Whether these fields belong on Workspace{Krylov} (since krylov_eigsolve can target Config{Thermalize, BohrDomain}) or only on Workspace{Thermalize}.
   - Recommendation: Include on Workspace{Krylov} as `Union{Nothing, ...}`, populated when config is Thermalize+BohrDomain.

2. **Whether to keep `KrylovWorkspace` as a type alias**
   - What we know: `KrylovWorkspace` is exported and used in tests, user code, and docstrings. Renaming to `Workspace{Krylov}` is a breaking change.
   - Recommendation: Add `const KrylovWorkspace{T} = Workspace{Krylov,...}` type alias for backward compatibility, or deprecation warning. This is a Phase 36 (API) concern but worth noting.

3. **Compilation time impact**
   - What we know: A single parametric struct with 4 type parameters will generate more specializations than 4 separate structs.
   - What's unclear: Whether this meaningfully increases compile time for the test suite.
   - Recommendation: Monitor. If compile time becomes an issue, the Scratch sub-structs can be moved to separate files for incremental compilation.

## Sources

### Primary (HIGH confidence)
- Direct codebase analysis of all workspace types, their fields, constructors, and usage sites
- Phase 34-02 SUMMARY confirming `KrausScratch.K0` is dead
- Existing allocation tests (`test_allocation.jl`, `test_krylov_matvec.jl`) providing zero-allocation baselines
- User decisions in CONTEXT.md providing locked architectural constraints

### Secondary (MEDIUM confidence)
- Julia documentation on parametric types and Union type narrowing behavior (from training data, consistent with observed codebase patterns)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - no new libraries, purely structural refactoring
- Architecture: HIGH - thorough audit of all 4 workspace types, all field usages, all constructors, all test references
- Pitfalls: HIGH - informed by existing allocation tests, known boxing issues (PD parameter), and user's explicit "one constructor" constraint
- Field audit: HIGH - every field traced to usage sites; dead K0 confirmed by Phase 34-02

**Research date:** 2026-02-26
**Valid until:** 2026-03-26 (stable codebase, no external dependencies changing)
