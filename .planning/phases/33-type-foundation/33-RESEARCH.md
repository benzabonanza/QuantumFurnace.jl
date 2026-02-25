# Phase 33: Type Foundation - Research

**Researched:** 2026-02-25
**Domain:** Julia type system refactoring -- singleton type hierarchies, parametric structs, trait-based dispatch
**Confidence:** HIGH

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- Struct: `Config{S,D,C,T}` -- single-letter type params (S=Simulation, D=Domain, C=Construction, T=Float)
- Simulation hierarchy: `AbstractSimulation` with subtypes `Lindbladian`, `Thermalize`, `KrylovSpectrum`, `Trajectory`
- Domain hierarchy: keep existing `AbstractDomain` with subtypes `BohrDomain`, `EnergyDomain`, `TimeDomain`, `TrotterDomain`
- Construction hierarchy: `AbstractConstruction` with subtypes `KMS`, `GNS`, `DLL`
- All three singletons stored as both fields AND type parameters: `sim::S`, `domain::D`, `construction::C`
- Hard swap -- remove old type names immediately, no aliases, no deprecation warnings
- Update all call sites in src/ and test/ in this phase
- All 4 simulation types defined in Phase 33 (not deferred to Phase 36)
- Dispatch approach: use full `Config{S,D,C,T}` where methods need config fields; dispatch on extracted singletons (e.g., `::KMS`, `::Lindbladian`) where only the type tag matters -- Claude's discretion per call site
- Trait function: `with_coherent(::KMS) = true`, `with_coherent(::GNS) = false`, `with_coherent(::DLL) = true` (placeholder)
- Remove `with_coherent` field from Config -- derived from construction type
- Single constructor with all keyword arguments: `Config(sim=Lindbladian(), domain=EnergyDomain(), construction=KMS(), N=4, beta=1.0, ...)`
- No convenience constructors, no defaults -- explicit is better
- No shortcut aliases

### Claude's Discretion
- Per-method decision on dispatching full Config vs extracted singleton
- Exact field ordering in Config struct
- How to handle the Bohr domain's different loop structure during migration
- Whether to introduce helper accessor functions (e.g., `simulation(config)`, `construction(config)`)

### Deferred Ideas (OUT OF SCOPE)
- DLL-specific behavior (no linear combinations, no energy labels) -- future phase when DLL construction is implemented
- KrylovSpectrum and Trajectory run_* entry points -- Phase 36 (types defined here, methods added there)
</user_constraints>

## Summary

Phase 33 replaces 4 duplicate config types (`LiouvConfig`, `ThermalizeConfig`, `LiouvConfigGNS`, `ThermalizeConfigGNS`) with a single parametric `Config{S,D,C,T}` struct, where S is a simulation singleton, D is a domain singleton, C is a construction singleton, and T is the float type. This eliminates approximately 130 lines of duplicated struct definitions and replaces ~80 `isa` checks and 4-way dispatch patterns with clean type-parameter dispatch.

The codebase currently has three distinct branching axes baked into the type system: simulation kind (Liouvillian vs Thermalize), domain (4 domains), and construction (KMS vs GNS). The first axis is currently encoded by using different struct types (`LiouvConfig` vs `ThermalizeConfig`), the third by using different struct types (`*Config` vs `*ConfigGNS`), and the second is already a type parameter. The refactoring unifies all three into type parameters on a single struct.

The key technical challenge is that the KMS/GNS distinction is currently detected via runtime `isa` checks in ~11 places across `results.jl`, `misc_tools.jl`, and `validate_config!`. Post-refactoring, these become type-parameter dispatch on the `C` parameter or trait function queries. The `with_coherent` field (a `Bool`) gets replaced by a compile-time trait `with_coherent(::KMS) = true`, `with_coherent(::GNS) = false`. All existing runtime checks of `config.with_coherent` become calls to `with_coherent(config.construction)`.

**Primary recommendation:** Define the three type hierarchies and the unified Config struct in `structs.jl`, then systematically migrate dispatch sites file-by-file, running tests after each file to catch regressions immediately.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Julia Base | 1.10+ | Parametric types, multiple dispatch, `@kwdef` | Native language features, zero dependencies |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| Julia `@kwdef` macro | Base | Keyword argument constructors | Config struct constructor |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Singleton types | Val{:KMS} etc. | Singletons are more readable, extensible, and allow method dispatch with `::KMS` syntax |
| Trait functions | Holy traits pattern (abstract type dispatch) | Simple trait functions are cleaner for single boolean queries; Holy traits better for grouping |

## Architecture Patterns

### Current Type Hierarchy (BEFORE)
```
AbstractConfig{D, T}
  AbstractLiouvConfig{D, T}
    LiouvConfig{D, T}         # KMS, Lindbladian
    LiouvConfigGNS{D, T}      # GNS, Lindbladian
  AbstractThermalizeConfig{D, T}
    ThermalizeConfig{D, T}    # KMS, Thermalize
    ThermalizeConfigGNS{D, T} # GNS, Thermalize
```

### New Type Hierarchy (AFTER)
```
AbstractSimulation
  Lindbladian
  Thermalize
  KrylovSpectrum
  Trajectory

AbstractDomain           # Already exists
  BohrDomain             # Already exists
  EnergyDomain           # Already exists
  TimeDomain             # Already exists
  TrotterDomain          # Already exists

AbstractConstruction
  KMS
  GNS
  DLL                    # Placeholder

Config{S <: AbstractSimulation, D <: AbstractDomain, C <: AbstractConstruction, T <: AbstractFloat}
```

### Pattern 1: Singleton Type Hierarchies
**What:** Define abstract type + concrete singleton subtypes for each axis.
**When to use:** When you need multiple dispatch on a conceptual category with a fixed set of variants.
**Example:**
```julia
abstract type AbstractSimulation end
struct Lindbladian    <: AbstractSimulation end
struct Thermalize     <: AbstractSimulation end
struct KrylovSpectrum <: AbstractSimulation end
struct Trajectory     <: AbstractSimulation end

abstract type AbstractConstruction end
struct KMS <: AbstractConstruction end
struct GNS <: AbstractConstruction end
struct DLL <: AbstractConstruction end
```

### Pattern 2: Unified Config Struct
**What:** Single parametric struct with all fields, type parameters encode the axes.
**When to use:** When 4+ struct types share 90%+ of their fields.
**Example:**
```julia
@kwdef struct Config{S <: AbstractSimulation, D <: AbstractDomain, C <: AbstractConstruction, T <: AbstractFloat}
    # Singleton fields (also encoded as type parameters)
    sim::S
    domain::D
    construction::C

    # Core physics parameters
    num_qubits::Int
    with_linear_combination::Bool
    beta::T
    sigma::T
    gaussian_parameters::Union{Tuple{T, T}, Tuple{Nothing, Nothing}} = (nothing, nothing)
    a::Union{T, Nothing} = nothing
    b::Union{T, Nothing} = nothing

    # Grid parameters
    num_energy_bits::Union{Int, Nothing} = nothing
    t0::Union{T, Nothing} = nothing
    w0::Union{T, Nothing} = nothing
    eta::Union{T, Nothing} = nothing
    num_trotter_steps_per_t0::Union{Int, Nothing} = nothing

    # Thermalize-specific (nothing for non-Thermalize simulations)
    mixing_time::Union{T, Nothing} = nothing
    delta::Union{T, Nothing} = nothing
end
```
**Important note on defaults:** The user decided "no defaults, explicit is better." However, the current `@kwdef` structs DO use defaults for optional fields like `a`, `b`, `num_energy_bits`, etc. (they default to `nothing`). The "no defaults" decision means no smart defaults for *required* parameters (sim, domain, construction, num_qubits, beta, sigma, etc.). The `nothing` defaults for truly optional fields (grid params that only some domains need) should be preserved, as they are functional necessities, not convenience shortcuts.

### Pattern 3: Trait Functions for Construction-Dependent Behavior
**What:** Free functions dispatching on the construction singleton to determine behavior.
**When to use:** When a boolean property is inherent to the construction type.
**Example:**
```julia
with_coherent(::KMS) = true
with_coherent(::GNS) = false
with_coherent(::DLL) = true   # placeholder

# Usage in existing code (BEFORE):
config.with_coherent || return nothing

# Usage (AFTER):
with_coherent(config.construction) || return nothing
```

### Pattern 4: Dispatch Migration Strategies
**What:** Converting `isa` checks and Union-typed dispatch to type-parameter dispatch.

#### Strategy A: Direct type-parameter dispatch (for physics functions)
```julia
# BEFORE:
pick_transition(config::LiouvConfig) = _pick_transition_kms(config)
pick_transition(config::ThermalizeConfig) = _pick_transition_kms(config)
pick_transition(config::LiouvConfigGNS) = _pick_transition_gns(config)
pick_transition(config::ThermalizeConfigGNS) = _pick_transition_gns(config)

# AFTER:
pick_transition(config::Config{<:Any, <:Any, KMS}) = _pick_transition_kms(config)
pick_transition(config::Config{<:Any, <:Any, GNS}) = _pick_transition_gns(config)
# Or equivalently, extract and dispatch on singleton:
pick_transition(config::Config) = _pick_transition(config.construction, config)
_pick_transition(::KMS, config) = _pick_transition_kms(config)
_pick_transition(::GNS, config) = _pick_transition_gns(config)
```

#### Strategy B: Abstract type dispatch (for simulation-kind dispatch)
```julia
# BEFORE:
function run_lindbladian(jumps, config::AbstractLiouvConfig, ...)
function run_thermalization(jumps, config::AbstractThermalizeConfig, ...)

# AFTER (Phase 33 -- keep existing function names, just change signatures):
function run_lindbladian(jumps, config::Config{Lindbladian}, ...)
function run_thermalization(jumps, config::Config{Thermalize}, ...)
```

#### Strategy C: Runtime `isa` replacement (for serialization/display)
```julia
# BEFORE:
d[:config_type] = (config isa Union{LiouvConfigGNS, ThermalizeConfigGNS}) ? "GNS" : "KMS"
d[:config_kind] = (config isa AbstractThermalizeConfig) ? "thermalize" : "liouv"

# AFTER:
d[:config_type] = config.construction isa GNS ? "GNS" : config.construction isa KMS ? "KMS" : "DLL"
d[:config_kind] = config.sim isa Thermalize ? "thermalize" : config.sim isa Lindbladian ? "liouv" : string(typeof(config.sim))
# Or better: define name traits
_construction_name(::KMS) = "KMS"
_construction_name(::GNS) = "GNS"
_construction_name(::DLL) = "DLL"
```

### Pattern 5: ThermalizeConfig -> LiouvConfig Elimination
**What:** The `_thermalize_to_liouv_config` conversion functions become unnecessary.
**Why:** With unified Config, you can simply dispatch on the domain type parameter instead.
**Example:**
```julia
# BEFORE: krylov_workspace.jl line 190-232
function _thermalize_to_liouv_config(tc::ThermalizeConfig)
    LiouvConfig(num_qubits=tc.num_qubits, ...) # Copy 14 fields
end

# AFTER: not needed. Functions that currently require AbstractLiouvConfig
# should instead accept Config{<:Any, D, C, T} and dispatch on D.
# The simulation type parameter is irrelevant for Lindbladian construction.
```

### Anti-Patterns to Avoid
- **Keeping `AbstractLiouvConfig` / `AbstractThermalizeConfig` as intermediate types:** These become redundant. Replace with `Config{Lindbladian,...}` and `Config{Thermalize,...}` dispatch.
- **Adding `with_coherent` back as a field "just in case":** The whole point is to derive it from construction type. If a function needs the boolean, call `with_coherent(config.construction)`.
- **Union types for mixed dispatch:** `Union{Config{Lindbladian,...}, Config{Thermalize,...}}` is a code smell. Use `Config{<:Any, D, C, T}` where D,C,T are the constraints that matter.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Keyword constructor | Manual positional constructor | `@kwdef` | Handles defaults, named args, type inference automatically |
| Type-checking on construction | Runtime `if/else` in constructor | Julia's type parameter dispatch | Compile-time guarantees, no runtime overhead |
| Config validation for GNS coherent | Runtime `isa` check + error | Parametric dispatch + trait | `with_coherent(::GNS) = false` is checked at the type level |

## Common Pitfalls

### Pitfall 1: Breaking `@kwdef` Type Inference
**What goes wrong:** `@kwdef` may not auto-generate the outer constructor correctly when there are 4 type parameters.
**Why it happens:** `@kwdef` generates an outer constructor that calls the inner one positionally. With 4 type params, the auto-generated outer constructor must infer S, D, C from the `sim`, `domain`, `construction` field values and T from `beta`.
**How to avoid:** Write an explicit outer constructor that infers type parameters:
```julia
function Config(; sim::S, domain::D, construction::C, beta::T, kwargs...) where {S,D,C,T}
    Config{S,D,C,T}(; sim=sim, domain=domain, construction=construction, beta=beta, kwargs...)
end
```
**Warning signs:** "MethodError: no method matching Config(...)" at test time.

### Pitfall 2: Mixing Time / Delta Fields for Non-Thermalize Configs
**What goes wrong:** `Config{Lindbladian,...}` will have `mixing_time` and `delta` fields that are meaningless.
**Why it happens:** Unified struct includes ALL fields from all 4 old structs.
**How to avoid:** Make them `Union{T, Nothing} = nothing` with defaults. Functions that need them (like `run_thermalization`) dispatch on `Config{Thermalize,...}` and access them directly. The `validate_config!` function should check that Thermalize configs have non-nothing values.
**Warning signs:** `nothing` errors at runtime when accessing `config.delta` on a Lindbladian config.

### Pitfall 3: `B_bohr` Only Defined for KMS
**What goes wrong:** `B_bohr` currently only accepts `Union{LiouvConfig, ThermalizeConfig}` (KMS types). If dispatch changes, GNS configs might accidentally reach it.
**Why it happens:** GNS has `with_coherent = false`, so the coherent code path short-circuits before calling `B_bohr`. This remains safe because `with_coherent(::GNS) = false`.
**How to avoid:** Ensure all coherent code paths check `with_coherent(config.construction)` BEFORE any construction-type-specific dispatch on B operators. The existing pattern (`config.with_coherent || return nothing`) already does this; just update to `with_coherent(config.construction) || return nothing`. Additionally, update `B_bohr` signature from `Union{LiouvConfig, ThermalizeConfig}` to `Config{<:Any, <:Any, KMS}`.
**Warning signs:** MethodError on `B_bohr` for GNS configs.

### Pitfall 4: `_select_b_plus_calculator` is KMS-Only
**What goes wrong:** This function on line 123 of `furnace_utensils.jl` dispatches on `Union{LiouvConfig, ThermalizeConfig}` (KMS). It is called inside the `with_coherent` branch, so it never runs for GNS.
**How to avoid:** Update signature to `Config{<:Any, <:Any, KMS}`. The guard `with_coherent(config.construction) || ...` in the caller ensures GNS/DLL never reach it.

### Pitfall 5: `_thermalize_to_liouv_config` Elimination Cascading Failures
**What goes wrong:** The KrylovWorkspace ThermalizeConfig constructor converts to LiouvConfig for dispatch. Removing this conversion requires updating all downstream functions.
**Why it happens:** `_precompute_data`, `_accumulate_R_total!`, `_precompute_coherent_total_B` currently dispatch on `AbstractLiouvConfig`. With unified Config, they need to accept any Config with the right domain parameter.
**How to avoid:** Change their signatures from `AbstractLiouvConfig{D}` to `Config{<:Any, D, <:Any}` (i.e., don't constrain on simulation type). Then the ThermalizeConfig workspace constructor can pass the original config directly.

### Pitfall 6: Results Serialization/Deserialization Breaks
**What goes wrong:** `_config_to_dict` and `_reconstruct_config` use `isa` checks against the old types. BSON files saved with old types won't load after the change.
**Why it happens:** The serialization format encodes `config_type` ("KMS"/"GNS") and `config_kind` ("liouv"/"thermalize") as strings in the Dict.
**How to avoid:** The Dict format already stores type info as strings, not Julia types. Update `_reconstruct_config` to map ("KMS", "thermalize") -> `Config(sim=Thermalize(), construction=KMS(), ...)`. Also add `sim` and `construction` string tags to `_config_to_dict`. Old BSON files will still load because the string-based reconstruction logic is independent of Julia's type system.

### Pitfall 7: Forgetting to Update Exports
**What goes wrong:** Module exports reference old type names that no longer exist.
**Why it happens:** Hard swap means no aliases.
**How to avoid:** Update `QuantumFurnace.jl` exports to replace `LiouvConfig, LiouvConfigGNS, ThermalizeConfig, ThermalizeConfigGNS` with `Config, AbstractSimulation, Lindbladian, Thermalize, KrylovSpectrum, Trajectory, AbstractConstruction, KMS, GNS, DLL, with_coherent`.

## Code Examples

### Example 1: Complete Config Struct Definition
```julia
# src/structs.jl

# --- Simulation hierarchy ---
abstract type AbstractSimulation end
struct Lindbladian    <: AbstractSimulation end
struct Thermalize     <: AbstractSimulation end
struct KrylovSpectrum <: AbstractSimulation end
struct Trajectory     <: AbstractSimulation end

# --- Construction hierarchy ---
abstract type AbstractConstruction end
struct KMS <: AbstractConstruction end
struct GNS <: AbstractConstruction end
struct DLL <: AbstractConstruction end

# --- Trait: with_coherent ---
with_coherent(::KMS) = true
with_coherent(::GNS) = false
with_coherent(::DLL) = true  # placeholder

# --- Unified Config ---
@kwdef struct Config{S <: AbstractSimulation, D <: AbstractDomain, C <: AbstractConstruction, T <: AbstractFloat}
    sim::S
    domain::D
    construction::C

    num_qubits::Int
    with_linear_combination::Bool
    beta::T
    sigma::T
    gaussian_parameters::Union{Tuple{T, T}, Tuple{Nothing, Nothing}} = (nothing, nothing)
    a::Union{T, Nothing} = nothing
    b::Union{T, Nothing} = nothing
    num_energy_bits::Union{Int, Nothing} = nothing
    t0::Union{T, Nothing} = nothing
    w0::Union{T, Nothing} = nothing
    eta::Union{T, Nothing} = nothing
    num_trotter_steps_per_t0::Union{Int, Nothing} = nothing

    # Thermalize-specific
    mixing_time::Union{T, Nothing} = nothing
    delta::Union{T, Nothing} = nothing
end
```

### Example 2: Outer Constructor with Type Inference
```julia
function Config(;
    sim::S, domain::D, construction::C,
    beta::T,
    kwargs...
) where {S <: AbstractSimulation, D <: AbstractDomain, C <: AbstractConstruction, T <: AbstractFloat}
    Config{S, D, C, T}(; sim=sim, domain=domain, construction=construction, beta=beta, kwargs...)
end
```

### Example 3: Migrating pick_transition Dispatch
```julia
# BEFORE (4 methods):
pick_transition(config::LiouvConfig) = _pick_transition_kms(config)
pick_transition(config::ThermalizeConfig) = _pick_transition_kms(config)
pick_transition(config::LiouvConfigGNS) = _pick_transition_gns(config)
pick_transition(config::ThermalizeConfigGNS) = _pick_transition_gns(config)

# AFTER (2 methods):
pick_transition(config::Config{<:Any, <:Any, KMS}) = _pick_transition_kms(config)
pick_transition(config::Config{<:Any, <:Any, GNS}) = _pick_transition_gns(config)
```

### Example 4: Migrating _pick_alpha Dispatch
```julia
# BEFORE (4 methods):
_pick_alpha(config::LiouvConfig) = _pick_alpha_kms(config)
_pick_alpha(config::ThermalizeConfig) = _pick_alpha_kms(config)
_pick_alpha(config::LiouvConfigGNS) = _pick_alpha_gns(config)
_pick_alpha(config::ThermalizeConfigGNS) = _pick_alpha_gns(config)

# AFTER (2 methods):
_pick_alpha(config::Config{<:Any, <:Any, KMS}) = _pick_alpha_kms(config)
_pick_alpha(config::Config{<:Any, <:Any, GNS}) = _pick_alpha_gns(config)
```

### Example 5: Migrating run_lindbladian / run_thermalization Signatures
```julia
# BEFORE:
function run_lindbladian(jumps, config::AbstractLiouvConfig{D,Tc}, hamiltonian::HamHam{Th}; ...) where {D, Tc, Th}
function run_thermalization(jumps, config::AbstractThermalizeConfig{D,Tc}, dm, hamiltonian::HamHam{Th}; ...) where {D, Tc, Th}

# AFTER:
function run_lindbladian(jumps, config::Config{Lindbladian,D,C,Tc}, hamiltonian::HamHam{Th}; ...) where {D, C, Tc, Th}
function run_thermalization(jumps, config::Config{Thermalize,D,C,Tc}, dm, hamiltonian::HamHam{Th}; ...) where {D, C, Tc, Th}
```

### Example 6: Migrating with_coherent Checks
```julia
# BEFORE (coherent.jl line 21):
config.with_coherent || return nothing

# AFTER:
with_coherent(config.construction) || return nothing
```

### Example 7: Migrating validate_config! GNS Check
```julia
# BEFORE (misc_tools.jl line 125):
if (config isa Union{LiouvConfigGNS, ThermalizeConfigGNS}) && config.with_coherent
    push!(errors, "GNS configs must have with_coherent=false")
end

# AFTER: This check becomes unnecessary -- with_coherent(::GNS) is always false.
# The validation that GNS doesn't use coherent is encoded in the type system.
# Remove this check entirely.
```

### Example 8: Test Helper Migration
```julia
# BEFORE:
function make_liouv_config(domain; with_coherent=true)
    LiouvConfig(num_qubits=NUM_QUBITS, with_coherent=with_coherent, ...)
end

# AFTER:
function make_liouv_config(domain; construction=KMS())
    Config(sim=Lindbladian(), domain=domain, construction=construction,
           num_qubits=NUM_QUBITS, with_linear_combination=true, ...)
end
```

## Detailed File-by-File Migration Audit

### Files That Define Config Types (must change)
| File | What Changes | Complexity |
|------|-------------|------------|
| `src/structs.jl` | Remove 4 old struct defs + 2 abstract types, add 3 hierarchies + Config struct + trait function | HIGH -- foundation of everything |
| `src/QuantumFurnace.jl` | Update exports | LOW |

### Files That Dispatch on Concrete Config Types (must change)
| File | Lines Affected | Pattern Used |
|------|---------------|--------------|
| `src/energy_domain.jl` | 6 lines: `pick_transition` 4-way dispatch, `_pick_transition_kms/gns` Union types | Strategy A |
| `src/bohr_domain.jl` | 10 lines: `B_bohr` Union types, `_pick_f` Union, `_pick_alpha` 4-way, `_pick_alpha_kms/gns` Union types | Strategy A |
| `src/furnace_utensils.jl` | 2 lines: `_select_b_plus_calculator` Union type | Strategy A |
| `src/misc_tools.jl` | 11 lines: `_generate_filename` 2x, `validate_config!` isa checks, `_print_press` 2x isa checks | Strategy C |
| `src/results.jl` | 8 lines: `_config_to_dict` isa, `_reconstruct_config`, `_write_companion_txt`, `_generate_experiment_filename`, `_default_results_dir` | Strategy C |
| `src/krylov_workspace.jl` | 4 lines: `_thermalize_to_liouv_config` 2x (DELETE), ThermalizeConfig workspace constructor | Strategy B + deletion |

### Files That Dispatch on Abstract Config Types (signature changes only)
| File | Lines Affected | Change Required |
|------|---------------|-----------------|
| `src/furnace.jl` | 3 signatures: `run_lindbladian`, `construct_lindbladian`, `run_thermalization` | `AbstractLiouvConfig` -> `Config{Lindbladian,...}`, `AbstractThermalizeConfig` -> `Config{Thermalize,...}` |
| `src/jump_workers.jl` | 6 signatures: `_jump_contribution!` 3x Liouv + 3x Thermalize | Same pattern |
| `src/krylov_matvec.jl` | 6 signatures: `apply_lindbladian!` etc. | `AbstractLiouvConfig{D}` -> `Config{<:Any, D, <:Any}` |
| `src/krylov_eigsolve.jl` | ~12 signatures + comments | `AbstractLiouvConfig` / `AbstractThermalizeConfig` -> `Config{...}` |
| `src/krylov_workspace.jl` | 5 signatures: workspace constructors + `_accumulate_R_total!` 3x | Same pattern |
| `src/coherent.jl` | 3 signatures + `config.with_coherent` checks | Signature + trait call |
| `src/trajectories.jl` | ~6 signatures | `AbstractThermalizeConfig` -> `Config{Thermalize,...}` or `Config{<:Any, D,...}` |
| `src/convergence.jl` | 2 signatures | `AbstractThermalizeConfig` -> `Config{Thermalize,...}` |
| `src/gap_estimation.jl` | 1 signature | Same |
| `src/furnace_utensils.jl` | 4 signatures: `_precompute_labels`, `_precompute_data` | `AbstractConfig{D}` -> `Config{<:Any, D, ...}` |

### Test Files (must change)
| File | Lines Affected | Change Required |
|------|---------------|-----------------|
| `test/test_helpers.jl` | 8 factory functions | Rewrite to use `Config(sim=..., construction=..., ...)` |
| `test/test_compilation.jl` | 2 `isa` checks | `isa LiouvConfig` -> `isa Config{Lindbladian}` etc. |
| `test/test_results.jl` | 6 `isa` checks + testset names | Same pattern |
| `test/test_gns_trajectory.jl` | Uses factories only | No direct changes if factories updated |
| `test/test_krylov_crossvalidation.jl` | 2 factory functions | Rewrite like test_helpers |
| `test/test_dm_detailed_balance.jl` | 1 direct `LiouvConfig()` call | Update to `Config()` |
| `test/old_tests/*.jl` | 4 files with direct `LiouvConfig()` calls | Update or leave (they're old tests) |

### Non-src Files (experiments, simulations, playground)
| Directory | Files Affected | Change Required |
|-----------|---------------|-----------------|
| `simulations/` | 3 files | Update `LiouvConfig`/`ThermalizeConfig` constructors |
| `experiments/` | 7 files | Update `LiouvConfig`/`ThermalizeConfig` constructors |
| `playground/` | 5 files | Update constructors |

## Key Design Decisions for Claude's Discretion

### 1. Field Ordering in Config Struct
**Recommendation:** Group by semantics, singletons first.
```
sim, domain, construction,         # Type-encoding singletons
num_qubits, with_linear_combination,  # System parameters
beta, sigma, gaussian_parameters,     # Physics parameters
a, b, eta,                            # Linear combination params
num_energy_bits, t0, w0,              # Grid parameters
num_trotter_steps_per_t0,             # Trotter-specific
mixing_time, delta                    # Thermalize-specific
```
**Rationale:** Follows the pattern of most-general-first (every config has sim/domain/construction) to most-specific-last (only Thermalize has mixing_time/delta). This matches how users mentally construct configs.

### 2. Dispatch Strategy Per Call Site
**Recommendation:** Use type-parameter dispatch (`Config{<:Any, <:Any, KMS}`) for physics functions that already have separate `_kms`/`_gns` implementations. Use extracted singleton dispatch (`with_coherent(config.construction)`) for boolean trait checks. Use runtime `isa` only for serialization code.

**Specific recommendations:**
| Function | Recommended Dispatch |
|----------|---------------------|
| `pick_transition` | `Config{<:Any, <:Any, KMS}` / `Config{<:Any, <:Any, GNS}` |
| `_pick_alpha` | Same |
| `B_bohr` | `Config{<:Any, <:Any, KMS}` (only KMS has B) |
| `_pick_f` | `Config{<:Any, <:Any, KMS}` (only KMS has f) |
| `_select_b_plus_calculator` | `Config{<:Any, <:Any, KMS}` |
| `_precompute_data` (Bohr, Liouv) | `Config{<:Any, BohrDomain}` with `Lindbladian`-compatible sim |
| `_precompute_data` (Bohr, Therm) | `Config{<:Any, BohrDomain}` with `Thermalize` sim (has extra bohr_keys) |
| `run_lindbladian` | `Config{Lindbladian, D, C, T}` |
| `run_thermalization` | `Config{Thermalize, D, C, T}` |
| `_jump_contribution!` (Liouv) | `Config{Lindbladian, D}` |
| `_jump_contribution!` (Therm) | `Config{Thermalize, D}` |
| `krylov_spectral_gap` (Liouv) | `Config{Lindbladian}` |
| `krylov_spectral_gap` (Therm) | `Config{Thermalize}` |
| `validate_config!` | `Config` (generic) |
| `_config_to_dict` | `Config` (generic, use field access) |
| `_print_press` | Separate `Config{Lindbladian}` and `Config{Thermalize}` (Thermalize prints extra fields) |

### 3. Helper Accessor Functions
**Recommendation:** Introduce a minimal set:
```julia
simulation(config::Config) = config.sim
construction(config::Config) = config.construction
```
**Rationale:** These read slightly better than `config.sim` and `config.construction` in dispatch chains. They also provide a stable API if the field name changes later. However, this is a very minor convenience and could be omitted without harm.

### 4. Bohr Domain's Different Loop Structure
**Recommendation:** Keep the existing split between `_precompute_data(config::..., BohrDomain)` for Lindbladian vs Thermalize paths. The Thermalize path precomputes `bohr_keys`, `bohr_is`, `bohr_js` for index caching. Post-migration:
```julia
# Lindbladian: simpler precompute (no index caching needed)
function _precompute_data(config::Config{Lindbladian, BohrDomain}, ham_or_trott)
    # ... alpha, gamma_norm_factor only
end

# Thermalize: richer precompute (with Bohr index caching)
function _precompute_data(config::Config{Thermalize, BohrDomain}, hamiltonian::HamHam)
    # ... alpha, gamma_norm_factor, bohr_keys, bohr_is, bohr_js
end
```
Actually the current split is `AbstractLiouvConfig{BohrDomain}` vs `AbstractThermalizeConfig{BohrDomain}`, which maps cleanly. The other domains (`EnergyDomain`, `TimeDomain/TrotterDomain`) use `AbstractConfig{D}` and don't split on simulation type, so they become `Config{<:Any, D}`.

### 5. Handling `_thermalize_to_liouv_config`
**Recommendation:** Delete both `_thermalize_to_liouv_config` methods entirely. The KrylovWorkspace constructor for Thermalize configs currently converts to LiouvConfig for:
1. `_precompute_data` dispatch
2. `_precompute_coherent_total_B` dispatch
3. `_accumulate_R_total!` dispatch

With unified Config, all these functions accept `Config{<:Any, D, C, T}` (not constrained on simulation type), so the original Thermalize config can be passed directly.

## Open Questions

1. **Should `_pick_transition_kms` and `_pick_transition_gns` accept `Config` or just extract the fields they need?**
   - What we know: They only use `config.beta`, `config.sigma`, `config.a`, `config.b`, `config.gaussian_parameters`, `config.with_linear_combination`. They don't use simulation type at all.
   - What's unclear: Whether it's cleaner to pass the full Config or destructure into fields.
   - Recommendation: Keep passing Config. The functions are not hot-path and the readability of `config.beta` is better than positional arguments.

2. **Should old BSON files be loadable after migration?**
   - What we know: `_reconstruct_config` uses string tags ("KMS"/"GNS", "liouv"/"thermalize") to select the constructor. The Dict format is independent of Julia types.
   - What's unclear: Whether there are any BSON files that encode the Julia struct type directly (not via Dict).
   - Recommendation: Update `_reconstruct_config` to return `Config(sim=..., construction=...)` instead of the old type names. The string-based format should be backward-compatible. Test with existing BSON files if available.

3. **What about `AbstractConfig` -- keep or remove?**
   - What we know: `AbstractConfig{D, T}` is used in 5 function signatures (`_precompute_labels`, `_precompute_data`, `_truncate_energy_labels`, `validate_config!`, `_collect_config_errors!`, `ExperimentResult`).
   - Recommendation: Remove `AbstractConfig` and replace with `Config` directly. Since Config is now the only concrete config type, there's no need for an abstract supertype. The parametric dispatch `Config{<:Any, D, <:Any, T}` replaces `AbstractConfig{D, T}`. Similarly remove `AbstractLiouvConfig` and `AbstractThermalizeConfig`.

## Sources

### Primary (HIGH confidence)
- Direct codebase analysis of all 28 source files in `src/`
- Direct analysis of all 30+ test files in `test/`
- `supplementary-informations/quantumfurnace-structure.md` -- user's architectural vision document

### Secondary (MEDIUM confidence)
- Julia documentation on parametric types, `@kwdef`, multiple dispatch (well-established language features)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- pure Julia type system, no external dependencies
- Architecture: HIGH -- direct codebase analysis of every affected file
- Pitfalls: HIGH -- identified through systematic audit of all dispatch sites and `isa` checks
- Migration scope: HIGH -- complete file-by-file audit with line counts

**Research date:** 2026-02-25
**Valid until:** 2026-03-25 (stable -- no external dependency changes expected)
