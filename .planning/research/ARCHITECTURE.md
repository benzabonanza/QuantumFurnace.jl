# Architecture Patterns: v2.2 Hamiltonian Simulation Time Counter

**Domain:** Analytical cost counting for quantum Gibbs sampler algorithm
**Researched:** 2026-03-04
**Confidence:** HIGH (based on exhaustive reading of all source files in existing codebase)

## Recommended Architecture

The Ham sim time counter is a **pure analytical computation module** that sits alongside the existing post-processing layer (fitting.jl, mixing.jl). It does NOT run simulations -- it computes the quantum algorithm's cost from parameters. This makes it architecturally clean: no workspace buffers, no scratch memory, no hot-loop concerns.

### High-Level Data Flow

```
HamHam (rescaling_factor, data size)
    |
    +-- Scalar params (r, delta, beta, sigma, transition weight type)
    |       |
    |       +-- MixingTimeEstimate (mixing_time field)  [optional convenience]
    |               |
    v               v
compute_simulation_time(ham, r, delta, mixing_time; beta, transition_weight=:smooth_metro, ...)
    |
    v
SimulationTimeBudget (result struct)
    |-- oft_time          (per step)
    |-- b_time            (per step)
    |-- per_step_time     (per step total)
    |-- total_steps
    |-- total_time
    |-- grid_points, energy_range, etc.
```

### Component Boundaries

| Component | Responsibility | Communicates With |
|-----------|---------------|-------------------|
| `SimulationTimeBudget` (struct) | Immutable result container for all time counts | Returned by compute functions, consumed by callers |
| `compute_simulation_time` (main API) | Orchestrate full time budget computation from parameters | Reads HamHam fields, calls internal helpers |
| `_qpe_grid_info` (internal) | Compute 2^r QPE grid point count and energy range | Used by OFT and B time helpers |
| `_oft_hamiltonian_time` (internal) | Compute Ham sim time for OFT (dissipative channels) per step | Uses rescaling_factor, r, transition weight |
| `_b_hamiltonian_time` (internal) | Compute Ham sim time for coherent B term per step | Uses rescaling_factor, beta, sigma |
| Transition weight reuse | Leverage existing `pick_transition` or compute weights analytically | Reads from energy_domain.jl transition weight functions |

### New File: `src/simulation_time.jl`

Single new file, approximately 200-350 lines. No modifications to existing source files except:
1. `src/QuantumFurnace.jl`: Add `include("simulation_time.jl")` line after `include("mixing.jl")` and add exports
2. No changes to HamHam, Config, Workspace, or any existing struct

This follows the established pattern: `fitting.jl` and `mixing.jl` were each added as standalone post-processing files with their own structs, included at the end of the module, exported from the main module file.

## Patterns to Follow

### Pattern 1: Immutable Result Struct (match FitResult, MixingTimeEstimate)

**What:** All computation outputs stored in an immutable struct with descriptive field names and full docstring.
**When:** Always -- this is the universal pattern in QuantumFurnace for computed results.
**Why:** FitResult (10 fields), BiexpFitResult (12 fields), MixingTimeEstimate (14 fields) all follow this. No mutable state, full data capture.

```julia
struct SimulationTimeBudget
    # Per-step costs (Ham sim time units)
    oft_time::Float64
    b_time::Float64
    per_step_time::Float64
    # Step count
    total_steps::Int
    mixing_time::Float64
    delta::Float64
    # Total cost
    total_time::Float64
    # QPE grid info
    r::Int
    grid_points::Int
    energy_range::Float64
    # Parameters used (provenance)
    n::Int
    beta::Float64
    rescaling_factor::Float64
    transition_weight::Symbol
    with_coherent::Bool
end
```

### Pattern 2: Pure Function + Keywords (match estimate_mixing_time)

**What:** Main API function takes a physics object + keyword arguments, dispatches internally on keyword values.
**When:** For the top-level `compute_simulation_time` API.
**Why:** `estimate_mixing_time(result; model=:single, ...)` established this pattern. Keywords for optional behavior, positional args for required data.

```julia
function compute_simulation_time(
    ham::HamHam,
    r::Int,
    delta::Float64,
    mixing_time::Float64;
    beta::Float64,
    sigma::Float64 = 1.0 / beta,
    transition_weight::Symbol = :smooth_metro,
    with_coherent::Bool = true,
    a::Float64 = 0.0,
    b::Float64 = 0.0,
) :: SimulationTimeBudget
```

### Pattern 3: Convenience Overload for MixingTimeEstimate

**What:** Allow passing a MixingTimeEstimate directly to extract mixing_time.
**When:** Common workflow: user runs `estimate_mixing_time` then wants simulation time budget.
**Why:** Avoids manual field extraction. Follows Julia multiple dispatch convention.

```julia
function compute_simulation_time(
    ham::HamHam, r::Int, delta::Float64,
    est::MixingTimeEstimate;
    kwargs...
)
    return compute_simulation_time(ham, r, delta, est.mixing_time; kwargs...)
end
```

### Pattern 4: Internal Helpers with Underscore Prefix

**What:** All internal computation helpers prefixed with `_`, not exported.
**When:** Always for non-public functions.
**Why:** Established codebase convention: `_check_fit_quality`, `_extrapolate_mixing_time`, `_log_linear_initial_guess`, `_create_energy_labels`.

## Anti-Patterns to Avoid

### Anti-Pattern 1: Coupling to Config Struct

**What:** Making `compute_simulation_time` accept a `Config` object.
**Why bad:** Config carries `sim::S`, `domain::D`, `construction::C` type parameters that are irrelevant to time counting. The counting function is independent of any specific simulation mode -- it counts the cost of the *ideal quantum algorithm*, not the classical simulation. Also, Config requires many fields (num_energy_bits, t0, w0, etc.) that are simulation parameters, not algorithm cost parameters.
**Instead:** Accept HamHam + scalar parameters directly.

### Anti-Pattern 2: Adding Fields to Existing Structs

**What:** Extending MixingTimeEstimate or ThermalizeResults with simulation time fields.
**Why bad:** Single responsibility violation. MixingTimeEstimate is about curve fitting; SimulationTimeBudget is about algorithm cost. Bundling them creates coupling where changes to one affect the other.
**Instead:** Separate `SimulationTimeBudget` struct from a standalone function.

### Anti-Pattern 3: Making This Domain-Dispatched

**What:** Creating separate `compute_simulation_time` methods for EnergyDomain, TrotterDomain, etc.
**Why bad:** The time counting is analytical. It counts the cost of the quantum algorithm (which runs on a quantum computer), not the cost of the classical simulation (which runs on the workstation). The quantum algorithm's OFT cost is the same regardless of which classical approximation domain was used.
**Instead:** A single function with transition weight type as a keyword.

### Anti-Pattern 4: Mutating HamHam or Creating Workspace

**What:** Storing results in mutable containers or creating scratch buffers.
**Why bad:** Time counting is pure arithmetic. It reads HamHam fields (rescaling_factor, data size) but never mutates anything. No matrices are constructed; no eigendecompositions happen.
**Instead:** Pure function returning immutable result struct.

## Integration Points with Existing Code

### Dependencies (what the new code reads from existing components)

| Existing Component | Fields/Functions Used | Why |
|-------------------|----------------------|-----|
| `HamHam{T}` | `rescaling_factor`, `size(ham.data, 1)` | Norm scaling determines physical time; dimension determines n |
| `MixingTimeEstimate` | `mixing_time` field | Convenience overload to chain with mixing time estimation |
| Transition weight logic (energy_domain.jl) | `pick_transition(config, w)` or analytical formulas | Compute transition weight values to sum over QPE grid |

### Critical Observation: HamHam Does Not Store beta or n Directly

- **n** (num_qubits) is inferred via `Int(log2(size(ham.data, 1)))` -- established pattern in `_extract_hamiltonian_params`
- **beta** is baked into `ham.gibbs` but NOT stored as a field on HamHam

**Design decision:** `compute_simulation_time` must accept `beta` as a required keyword argument. This is consistent with how the codebase works -- `HamHam(raw, beta)` requires beta at construction, `load_hamiltonian(type, n; beta=...)` requires it at load time. Beta is always an external parameter.

### Transition Weight Reuse Strategy

The existing `pick_transition(config, w)` dispatches on Config construction type (KMS/GNS) and uses Config fields (a, b, beta, sigma, gaussian_parameters). For the time counter, we have two options:

**Option A: Reuse `pick_transition` by constructing a minimal Config**
- Pro: No code duplication
- Con: Requires constructing a fake Config just for transition evaluation; Config validation would need to be bypassed

**Option B: Accept a transition weight function directly**
- Pro: Clean separation; no Config dependency
- Con: Caller must know how to construct the function

**Option C (recommended): Accept weight type as Symbol + parameters, implement lightweight dispatch internally**
- Accept `transition_weight::Symbol` (`:gaussian`, `:metropolis`, `:smooth_metro`)
- Accept relevant parameters (a, b, beta, sigma, gaussian_parameters) as keywords
- Internally call the correct formula
- Pro: Self-contained; no Config dependency; clear API
- Con: Small amount of formula duplication (but these are one-line formulas)

This is the cleanest approach. The formulas in `pick_transition` are compact analytical expressions that can be called directly without constructing Config objects.

### No Modifications to Existing Structs or Functions

The new module is purely additive:
- **New file:** `src/simulation_time.jl`
- **Modified file:** `src/QuantumFurnace.jl` (add `include` + `export` lines)
- **New test file:** `test/test_simulation_time.jl`
- **Modified test file:** `test/runtests.jl` (add `include` line)

Zero changes to HamHam, Config, Workspace, or any Result struct.

## Data Flow Detail: The Physics of Time Counting

### What Is Being Counted

The quantum Gibbs sampler algorithm runs on a quantum computer. Each step of the algorithm requires:
1. **OFT (Operator Fourier Transform):** For each energy grid point w, simulate the Hamiltonian for time proportional to the QPE circuit depth. This creates the frequency-filtered jump operators A(w).
2. **B (Coherent correction):** For KMS detailed balance, a coherent unitary correction that involves further Hamiltonian simulation.

The "Hamiltonian simulation time" is the total time the quantum computer spends simulating `e^{iHt}` across all steps of the algorithm. This is the primary cost metric for the quantum algorithm.

### QPE Grid

With `r` estimating qubits:
- Grid has `N = 2^r` points
- Energy spacing: `w0` (determined by desired resolution)
- Time unit: `t0 = 2*pi / (N * w0)` (Fourier relation)
- Each OFT evaluation at energy `w` sums over time grid points: the Hamiltonian is simulated for total time proportional to `N * t0` = `2*pi / w0`
- In rescaled units, the Hamiltonian has norm ~0.45, so physical simulation time scales with `rescaling_factor`

### Per-Step OFT Cost

For each thermalization step:
- Sum over truncated energy grid points `w` that have non-negligible transition weight
- At each `w`: Ham sim time = number of time grid points * t0
- Weight by transition function value (Gaussian/Metro/SmoothMetro)
- The per-step OFT cost depends on how many energy grid points are "active" (non-negligible weight)

### Per-Step B Cost

For KMS with coherent term:
- B involves convolutions `b_minus` and `b_plus` over the time grid
- Each time label `t` requires Ham sim for duration `t` (or `s * beta` for `b_plus`)
- The B cost is determined by the number of non-negligible time labels in b_minus and b_plus

### Total Steps

```
total_steps = ceil(Int, mixing_time / delta)
```

## Suggested Build Order

### Phase 1: Struct + Grid Utilities (zero physics dependencies)

1. Define `SimulationTimeBudget` struct with full docstring
2. Implement `_qpe_grid_info(r, w0)` returning `(grid_points, energy_range, t0)`
3. Implement `_num_qubits(ham)` helper (just `Int(log2(size(ham.data, 1)))`)
4. Unit tests for struct construction and grid arithmetic

### Phase 2: OFT Time Counting (core cost computation)

1. Implement `_compute_transition_weight(w, transition_weight, beta, sigma, a, b, gaussian_params)` -- lightweight dispatch on Symbol
2. Implement `_oft_hamiltonian_time(ham, r, w0, beta, sigma, transition_weight, ...)` -- sum over energy grid
3. Must handle: `:gaussian`, `:metropolis`, `:smooth_metro` transition weight types
4. Unit tests: verify OFT time scales correctly with r, n; verify against hand calculations

### Phase 3: B Term Time Counting (coherent correction cost)

1. Implement `_b_hamiltonian_time(ham, r, w0, beta, sigma, ...)` -- coherent term cost
2. Returns 0.0 when `with_coherent=false`
3. Unit tests: verify 0 when disabled; verify positive and plausible when enabled

### Phase 4: Integration API + Validation

1. Implement main `compute_simulation_time(ham, r, delta, mixing_time; ...)` orchestrator
2. Add `MixingTimeEstimate` convenience overload
3. Wire into module: add include + exports to `src/QuantumFurnace.jl`
4. Add `include("test_simulation_time.jl")` to `test/runtests.jl`
5. Integration tests: end-to-end from HamHam to SimulationTimeBudget
6. Validation: compare against known analytical results or paper formulas

**Rationale for ordering:**
- Phase 1 is pure infrastructure with no physics -- safe to build and test independently
- Phase 2 is the core: OFT dominates the time budget in practice
- Phase 3 adds the optional coherent correction
- Phase 4 ties everything together; can validate the full pipeline

## Sources

- Direct source code analysis of all files in `src/` directory of QuantumFurnace.jl
- Existing architecture documentation: `.planning/codebase/ARCHITECTURE.md`
- Existing patterns from fitting.jl (Phase 42) and mixing.jl (Phase 42-43)
- Milestone context: v2.2 Hamiltonian Simulation Time Counter
