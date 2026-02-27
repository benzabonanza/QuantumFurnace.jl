# Phase 36: API and Results - Research

**Researched:** 2026-02-27
**Domain:** Julia API design, struct serialization, BSON round-trip
**Confidence:** HIGH

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

#### Entry point signatures
- Keep separate positional args: `run_*(jumps, config, hamiltonian; kwargs...)` -- Config stays focused on simulation parameters, not problem definition
- Hamiltonian and trotter are both positional args at the same level (both are eigenbasis objects, used interchangeably as `ham_or_trott`): `run_*(jumps, config, hamiltonian, trotter; kwargs...)`
- Trotter defaults to `nothing` when not used
- Trajectory-specific kwargs (ntraj, total_time, delta, seed, save_every, observables) stay as flat keyword args, no options struct

#### Trajectory API consolidation
- All 4 current trajectory functions (`run_trajectories`, `run_observable_trajectories`, `run_trajectories_convergence`, `run_trajectories_adaptive`) collapse into a single `run_trajectory`
- `convergence=true` keyword triggers convergence tracking; Gibbs reference comes from `hamiltonian.gibbs` (already available in HamHam struct)
- `adaptive=true` keyword (requires `convergence=true`) enables adaptive early stopping with additional kwargs (convergence_threshold, patience, etc.)
- DM reconstruction logic: no observables = reconstruct DM at save points (mid-run); with observables = track observables mid-run, reconstruct DM once at the end as final average

#### Result struct contents
- 4 typed Result structs: `LindbladResults`, `ThermalizeResults`, `KrylovSpectrumResults`, `TrajectoryResults`
- `ExperimentResult` is removed entirely -- replaced by the 4 typed structs
- `AbstractResults` as common supertype if needed for dispatch
- Each Result stores the full `Config{S,D,C,T}` for reproducibility
- Metadata in all results: git hash, timestamp, wall_time_seconds, thread count (no Julia version)
- `LindbladResults`: spectral data only (eigenvalues, fixed_point, gap_mode, spectral_gap) -- no full Liouvillian matrix
- `ThermalizeResults`: includes trace distances over time for convergence-to-Gibbs plotting
- `TrajectoryResults`: convergence data and observable time series are `Union{Nothing, ...}` -- only populated when those features were used

#### Save/load workflow
- Standalone `save_result(result, path)` function -- not a keyword on run_* functions
- Every `save_result` call always creates a companion `.txt` file with: date, config summary (n_qubits, beta, domain, construction), key results
- `load_result(path)` auto-detects the result type from a tag stored in the BSON -- returns the correct typed Result
- Clean break with old `ExperimentResult` .bson files -- no backward compatibility needed

### Claude's Discretion
- Internal serialization format (Dict-based vs direct struct BSON)
- Companion .txt content formatting and level of detail per result type
- Whether `AbstractResults` is actually needed or if duck typing suffices
- Auto-filename generation patterns (if any)
- How metadata is captured internally (_capture_metadata helper pattern)

### Deferred Ideas (OUT OF SCOPE)
None -- discussion stayed within phase scope
</user_constraints>

## Summary

Phase 36 consolidates the QuantumFurnace public API from 7+ entry points and inconsistent result types into 4 clean entry points (`run_lindblad`, `run_thermalize`, `run_krylov_spectrum`, `run_trajectory`) each returning a typed Result struct. The existing codebase already has well-factored internals -- the key work is reorganizing the public surface and unifying the serialization layer.

The codebase currently has these entry points: `run_lindbladian` (furnace.jl:1), `run_thermalization` (furnace.jl:95), `run_trajectories` (trajectories.jl:566), `run_observable_trajectories` (trajectories.jl:680), `run_trajectories_convergence` (convergence.jl:162), `run_trajectories_adaptive` (convergence.jl:289), and `krylov_spectral_gap` (krylov_eigsolve.jl:389,490). The trajectory consolidation is the most complex task -- collapsing 4 functions into one `run_trajectory` with keyword-driven mode switching. The Lindbladian, Thermalize, and KrylovSpectrum paths each need new Result structs but their internal logic is already clean.

**Primary recommendation:** Define the 4 Result structs first, then build the 4 `run_*` entry points wrapping existing internals, then implement `save_result`/`load_result` with Dict-based BSON (continuing the proven pattern from `ExperimentResult`), and finally update `simulations/` scripts. The trajectory consolidation should be done as a careful composition of the 3 existing chunk runners (`_run_batch_no_obs!`, `_run_chunk_with_obs!`, `_run_chunk_obs_only!`).

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| BSON.jl | 0.3 | Result serialization to disk | Already used, proven round-trip in `results.jl` |
| LibGit2 | stdlib | Git hash capture for metadata | Already used in `_capture_git_hash()` |
| Dates | stdlib | Timestamp capture for metadata | Already used in `_capture_metadata()` |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| Pkg | stdlib | Project root path resolution | For default save directories |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| BSON Dict-based | Direct struct BSON | Direct struct BSON breaks when struct definitions change; Dict-based is proven safer for forward/backward compat in this codebase |
| `AbstractResults` supertype | Duck typing | Supertype is better: enables `save_result(::AbstractResults)` dispatch and `load_result` return type annotation |

**Recommendation (Claude's Discretion): Use `AbstractResults` supertype.**
Reason: The `save_result` function needs to dispatch on a common type to avoid 4 separate method definitions. `load_result` return type annotation `::AbstractResults` documents intent. Cost is one abstract type definition (trivial). Duck typing would work but provides no dispatch benefit for the serialization layer.

## Architecture Patterns

### Current Codebase Structure (relevant files)
```
src/
  QuantumFurnace.jl     # Module definition, exports
  structs.jl            # Config, JumpOp, DMSimulationResult, LindbladianResult,
                        #   ConvergenceData, TrajectoryResult, ObservableTrajectoryResult,
                        #   Workspace, scratch structs
  results.jl            # ExperimentResult, save/load, Dict-based BSON, companion .txt
  furnace.jl            # run_lindbladian, run_thermalization, construct_lindbladian
  trajectories.jl       # run_trajectories, run_observable_trajectories, step_along_trajectory!
  convergence.jl        # run_trajectories_convergence, run_trajectories_adaptive
  krylov_eigsolve.jl    # krylov_spectral_gap (2 dispatch paths)
  krylov_workspace.jl   # Workspace constructors for KrylovSpectrum
```

### Pattern 1: New Result Structs in structs.jl

**What:** Define `AbstractResults` and 4 concrete Result structs alongside the existing types in `structs.jl`.

**Why here:** `structs.jl` is included early in the module (line 110 of QuantumFurnace.jl), before `results.jl`. All type definitions live here. The new Result structs need to be available before `results.jl` defines serialization.

**Struct design:**
```julia
abstract type AbstractResults end

struct LindbladResults{T<:AbstractFloat} <: AbstractResults
    config::Config
    eigenvalues::Vector{Complex{T}}    # leading eigenvalues (from Arpack eigs)
    fixed_point::Matrix{Complex{T}}
    gap_mode::Matrix{Complex{T}}
    spectral_gap::Complex{T}
    metadata::Dict{Symbol, Any}
end

struct ThermalizeResults{T<:AbstractFloat} <: AbstractResults
    config::Config
    final_dm::Matrix{Complex{T}}
    trace_distances::Vector{T}
    time_steps::Vector{T}
    metadata::Dict{Symbol, Any}
end

struct KrylovSpectrumResults{T<:AbstractFloat} <: AbstractResults
    config::Config
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
    metadata::Dict{Symbol, Any}
end

struct TrajectoryResults{T<:AbstractFloat} <: AbstractResults
    config::Config
    rho_mean::Matrix{Complex{T}}
    n_trajectories::Int
    seed::Int
    # Observable data (Nothing when no observables used)
    times::Union{Nothing, Vector{Float64}}
    measurements_mean::Union{Nothing, Matrix{Float64}}
    # Convergence data (Nothing when convergence=false)
    convergence::Union{Nothing, ConvergenceData}
    metadata::Dict{Symbol, Any}
end
```

### Pattern 2: run_* Entry Points Wrap Existing Internals

**What:** Each `run_*` function wraps existing internal logic, times it, captures metadata, and returns the typed Result.

**Key insight:** The existing internal functions (`construct_lindbladian`, `_build_framework_and_seed`, `_run_batch_no_obs!`, etc.) are already well-factored. The new entry points are thin wrappers that add timing and metadata capture.

**Example pattern:**
```julia
function run_lindblad(jumps::Vector{JumpOp}, config::Config{Lindbladian,D,C,Tc},
                      hamiltonian::HamHam{Th}, trotter::Union{TrottTrott, Nothing}=nothing
                      ) where {D, C, Tc<:AbstractFloat, Th<:AbstractFloat}
    # ... validation, timing wrapper ...
    t_start = time()
    # ... existing logic from current run_lindbladian ...
    wall_time = time() - t_start
    metadata = _capture_metadata(wall_time_seconds=wall_time)
    return LindbladResults(config, eigenvalues, fixed_point, gap_mode, spectral_gap, metadata)
end
```

### Pattern 3: Trajectory Consolidation via Keyword Dispatch

**What:** Single `run_trajectory` function dispatches to different internal code paths based on keyword arguments.

**Design:**
```julia
function run_trajectory(
    jumps::Vector{JumpOp},
    config::Config{Thermalize},
    hamiltonian::HamHam,
    trotter::Union{TrottTrott, Nothing}=nothing;
    ntraj::Int = 1,
    total_time::Real = config.mixing_time,
    delta::Real = config.delta,
    seed::Union{Int, Nothing} = nothing,
    save_every::Int = 1,
    observables::Union{Nothing, Vector{<:Matrix{<:Complex}}} = nothing,
    observable_names::Union{Nothing, Vector{String}} = nothing,
    convergence::Bool = false,
    adaptive::Bool = false,
    # Convergence-specific kwargs
    batch_size::Int = 1000,
    n_batches::Int = 10,
    n_max::Int = 20_000,
    convergence_threshold::Float64 = 0.01,
    patience::Int = 3,
    min_batches::Int = 5,
    window_size::Int = 3,
)
```

**Dispatch logic:**
1. `adaptive=true` (requires `convergence=true`): delegates to adaptive convergence path
2. `convergence=true, adaptive=false`: delegates to batch convergence path
3. `observables !== nothing`: delegates to observable trajectory path (with DM reconstruction at end)
4. Default: plain trajectory path (DM reconstruction at save points)

**Gibbs reference for convergence:** `hamiltonian.gibbs` is used directly (no separate gibbs kwarg). For TrotterDomain, compute `_gibbs_in_trotter_basis(hamiltonian, trotter)` internally.

### Pattern 4: Dict-Based BSON Serialization (Proven Pattern)

**What:** Continue the existing `_experiment_to_dict` / `_dict_to_experiment` pattern from `results.jl` for all 4 Result types.

**Why:** The codebase already proves this works. Direct struct BSON (storing Julia type tags) breaks when struct definitions change across versions. Dict-based serialization with string type tags is robust.

**Type tag:** Store a `:result_type` key in the Dict: `"lindblad"`, `"thermalize"`, `"krylov_spectrum"`, `"trajectory"`. `load_result` reads this tag and dispatches to the correct reconstruction function.

### Pattern 5: Companion .txt Formatting Per Result Type

**What:** Each Result type gets a specialized companion `.txt` section showing its key results.

**Recommendation (Claude's Discretion):**
```
=== QuantumFurnace [LindbladResults|ThermalizeResults|...] ===

Date:       2026-02-27_12:00:00
Git:        abc123
Threads:    4
Wall time:  1.5 s

--- Config ---
Type:       KMS
Domain:     EnergyDomain
n_qubits:   4
beta:       10.0

--- Results ---
[Type-specific section]
```

LindbladResults: spectral_gap, fixed_point dimension
ThermalizeResults: final trace distance, number of steps
KrylovSpectrumResults: spectral_gap, matvec_count, converged eigenvalues
TrajectoryResults: n_trajectories, seed, rho_mean dimension, convergence status if applicable

### Anti-Patterns to Avoid
- **Storing the full Liouvillian in LindbladResults:** Decision is locked -- spectral data only. The Liouvillian is dim^2 x dim^2 which is prohibitively large at scale.
- **Making save_path a keyword on run_*:** Decision is locked -- standalone `save_result` function.
- **Backward compatibility with ExperimentResult .bson files:** Decision is locked -- clean break.
- **Storing Julia version in metadata:** Decision is locked -- only git hash, timestamp, wall_time_seconds, thread count.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Git hash capture | Custom shell exec | `_capture_git_hash()` (existing, uses LibGit2) | Already works, handles errors gracefully |
| BSON serialization | Custom binary format | BSON.jl with Dict-based conversion | Already proven in ExperimentResult round-trip |
| Hermitian wrap for BSON | Store raw Matrix | `Matrix(hermitian)` unwrap before save | BSON cannot serialize `Hermitian` wrapper; existing code already does this |
| Timestamp formatting | Custom string concat | `Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH:MM:SS")` | Already used in `_capture_metadata()` |

**Key insight:** Almost all serialization infrastructure already exists in `results.jl`. The work is extending the pattern to 4 types, not inventing new infrastructure.

## Common Pitfalls

### Pitfall 1: BSON Cannot Serialize Julia Type Wrappers
**What goes wrong:** `Hermitian`, `Diagonal`, `Adjoint` wrappers cause BSON deserialization failures when struct definitions change.
**Why it happens:** BSON stores Julia type tags; if the type changes between save/load, deserialization fails.
**How to avoid:** Always unwrap to `Matrix()` before saving: `Matrix(result.fixed_point)`, `Matrix(gibbs)`. The existing code already does this in `_trajectory_to_dict`.
**Warning signs:** `MethodError` or `UndefVarError` on `BSON.load`.

### Pitfall 2: Tuple-to-Array Conversion in BSON
**What goes wrong:** BSON stores Julia tuples as arrays. On load, `(1.0, 2.0)` becomes `[1.0, 2.0]`.
**Why it happens:** BSON format limitation.
**How to avoid:** Convert back in reconstruction: `if gp isa AbstractVector; kwargs[:gaussian_parameters] = (gp[1], gp[2]); end`. Already handled in `_dict_to_config_kwargs`.
**Warning signs:** Config reconstruction fails with type mismatch.

### Pitfall 3: Breaking the run_trajectories Seed Scheme
**What goes wrong:** Changing seed generation or per-trajectory seed assignment breaks reproducibility.
**Why it happens:** Trajectory results depend on exact `Xoshiro(master_seed + traj_id)` per trajectory.
**How to avoid:** The new `run_trajectory` must use identical seeding as existing functions: `_build_framework_and_seed` generates `actual_seed`, each trajectory gets `Xoshiro(actual_seed + traj_id)`.
**Warning signs:** Same seed produces different density matrices before/after refactor.

### Pitfall 4: Forgetting to Update Exports
**What goes wrong:** New types and functions are defined but not exported; user code cannot access them.
**Why it happens:** `QuantumFurnace.jl` has manual export lists.
**How to avoid:** Update the export block in `QuantumFurnace.jl` with: new Result types, new run_* functions, `save_result`, `load_result`, `AbstractResults`. Also remove deprecated exports.
**Warning signs:** `UndefVarError` when calling from user code.

### Pitfall 5: Config Type Mismatch Between run_trajectory and Internals
**What goes wrong:** The unified `run_trajectory` function receives `Config{Thermalize}` but internal convergence functions expect specific gibbs/observable arguments.
**Why it happens:** `run_trajectories_convergence` requires explicit `gibbs::Hermitian` kwarg. The new API gets Gibbs from `hamiltonian.gibbs`.
**How to avoid:** Compute gibbs internally: `gibbs = config.domain isa TrotterDomain ? _gibbs_in_trotter_basis(hamiltonian, trotter) : hamiltonian.gibbs`.
**Warning signs:** Missing keyword argument errors or wrong-basis Gibbs comparison.

### Pitfall 6: Metadata Capture Must Exclude Julia Version
**What goes wrong:** Including `julia_version` in metadata when the decision explicitly excludes it.
**Why it happens:** Existing `_capture_metadata` includes `:julia_version`.
**How to avoid:** Modify `_capture_metadata` or use a new helper that omits Julia version. The locked decision specifies: git hash, timestamp, wall_time_seconds, thread count only.
**Warning signs:** Extra field in saved BSON.

## Code Examples

### Example 1: Current run_lindbladian Signature (to be replaced)
```julia
# Current (src/furnace.jl:1-40)
function run_lindbladian(jumps::Vector{JumpOp}, config::Config{Lindbladian,D,C,Tc},
    hamiltonian::HamHam{Th};
    trotter::Union{TrottTrott, Nothing}=nothing) where {D, C, Tc, Th}
    # ... returns LindbladianResult (has liouvillian field -- removed in new API)
end
```

### Example 2: Current ExperimentResult Save Pattern (to be adapted)
```julia
# Current (src/results.jl:48-55)
function _experiment_to_dict(result::ExperimentResult)
    return Dict{Symbol, Any}(
        :config             => _config_to_dict(result.config),
        :trajectory         => _trajectory_to_dict(result.trajectory_result),
        :hamiltonian_params => result.hamiltonian_params,
        :metadata           => result.metadata,
    )
end
```

### Example 3: Auto-detect Result Type on Load
```julia
function load_result(path::String)
    d = BSON.load(path)
    tag = d[:result_type]  # "lindblad", "thermalize", "krylov_spectrum", "trajectory"
    if tag == "lindblad"
        return _dict_to_lindblad_results(d)
    elseif tag == "thermalize"
        return _dict_to_thermalize_results(d)
    elseif tag == "krylov_spectrum"
        return _dict_to_krylov_spectrum_results(d)
    elseif tag == "trajectory"
        return _dict_to_trajectory_results(d)
    else
        error("Unknown result type: $tag")
    end
end
```

### Example 4: run_trajectory Keyword Dispatch Logic
```julia
function run_trajectory(jumps, config, hamiltonian, trotter=nothing;
                        convergence=false, adaptive=false, observables=nothing, ...)
    # Validation
    adaptive && !convergence && error("adaptive=true requires convergence=true")

    # Build workspace once
    ws, actual_seed = _build_framework_and_seed(jumps, config, psi0, hamiltonian; ...)

    # Gibbs reference for convergence
    gibbs = if convergence
        config.domain isa TrotterDomain ?
            _gibbs_in_trotter_basis(hamiltonian, trotter) : hamiltonian.gibbs
    else
        nothing
    end

    t_start = time()
    if adaptive
        # ... adaptive convergence path (from run_trajectories_adaptive)
    elseif convergence
        # ... batch convergence path (from run_trajectories_convergence)
    elseif observables !== nothing
        # ... observable trajectory path
    else
        # ... plain trajectory path
    end
    wall_time = time() - t_start

    metadata = _capture_metadata(wall_time_seconds=wall_time)
    return TrajectoryResults(config, rho_mean, ntraj, actual_seed, times, measurements,
                             convergence_data, metadata)
end
```

### Example 5: run_trajectory Needs psi0 as Argument
**Critical observation:** The current trajectory functions take `psi0::Vector{<:Complex}` as a positional argument. The CONTEXT.md signature is `run_*(jumps, config, hamiltonian, trotter; kwargs...)` which does NOT include `psi0`. However, trajectories REQUIRE an initial state vector. There are two options:
1. Add `psi0` as a keyword argument with a default (e.g., computational basis |0>)
2. Add `psi0` as a positional argument between config and hamiltonian

**Recommendation:** Since the other 3 entry points don't need `psi0`, and trajectories always need it, make `psi0` a required keyword argument of `run_trajectory`. This keeps the positional signature uniform across all 4 entry points.

### Example 6: run_thermalize Needs initial_dm
**Critical observation:** Similarly, `run_thermalization` currently takes `evolving_dm::Matrix{<:Complex}` as a positional argument. The new `run_thermalize` needs an initial density matrix. Same recommendation: make it a required keyword argument or use a default (maximally mixed state I/d).

## Mapping: Current Functions to New API

| Current Function | New API | Notes |
|-----------------|---------|-------|
| `run_lindbladian` (furnace.jl) | `run_lindblad` | Returns `LindbladResults` (no Liouvillian matrix) |
| `run_thermalization` (furnace.jl) | `run_thermalize` | Returns `ThermalizeResults` with trace distances |
| `krylov_spectral_gap` (krylov_eigsolve.jl) | `run_krylov_spectrum` | Wraps existing, returns `KrylovSpectrumResults` |
| `run_trajectories` (trajectories.jl) | `run_trajectory` | Plain mode (no obs, no convergence) |
| `run_observable_trajectories` | `run_trajectory` with `observables=...` | Observable mode |
| `run_trajectories_convergence` | `run_trajectory` with `convergence=true` | Convergence mode |
| `run_trajectories_adaptive` | `run_trajectory` with `convergence=true, adaptive=true` | Adaptive mode |
| `ExperimentResult` + `save_experiment` | `save_result` + typed Results | Clean break |
| `load_experiment` | `load_result` | Auto-detects type via tag |

## Existing Types to Keep, Remove, or Rename

| Type | Action | Reason |
|------|--------|--------|
| `LindbladianResult` | **Replace** with `LindbladResults` | Different fields (no liouvillian), different name |
| `DMSimulationResult` | **Replace** with `ThermalizeResults` | Adds config, metadata, renames fields |
| `KrylovGapResult` | **Replace** with `KrylovSpectrumResults` | Adds config, metadata |
| `TrajectoryResult` | **Keep as internal** or merge into `TrajectoryResults` | Used internally for batch aggregation |
| `ObservableTrajectoryResult` | **Remove** | Absorbed into unified `TrajectoryResults` |
| `ExperimentResult` | **Remove** | Replaced by 4 typed Results |
| `ConvergenceData` | **Keep** | Embedded in `TrajectoryResults.convergence` |

## Metadata Capture

**Recommendation (Claude's Discretion):** Modify the existing `_capture_metadata` to match the locked decision (no Julia version):

```julia
function _capture_metadata(;
    n_threads::Int = Threads.nthreads(),
    wall_time_seconds::Union{Float64, Nothing} = nothing,
)
    return Dict{Symbol, Any}(
        :timestamp         => Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH:MM:SS"),
        :git_hash          => _capture_git_hash(),
        :n_threads         => n_threads,
        :wall_time_seconds => wall_time_seconds,
    )
end
```

The existing `_capture_git_hash()` function is kept as-is.

## Auto-Filename Generation

**Recommendation (Claude's Discretion):** Keep the existing `_generate_experiment_filename` pattern but adapt it:

```julia
function _generate_result_filename(result::AbstractResults)
    cfg = result.config
    type_str = _result_type_tag(result)  # "lindblad", "thermalize", etc.
    db_str = cfg.construction isa GNS ? "gns" : cfg.construction isa KMS ? "kms" : "dll"
    domain_str = lowercase(replace(string(typeof(cfg.domain)), "Domain" => ""))
    n_str = "n$(cfg.num_qubits)"
    beta_str = "beta$(round(Int, cfg.beta))"
    date_str = Dates.format(Dates.now(), dateformat"yyyymmdd")
    return "$(type_str)_$(db_str)_$(n_str)_$(beta_str)_$(domain_str)_$(date_str).bson"
end
```

## Open Questions

1. **psi0 / initial_dm in run_trajectory and run_thermalize**
   - What we know: `run_trajectory` needs `psi0::Vector{<:Complex}`, `run_thermalize` needs `evolving_dm::Matrix{<:Complex}`. The CONTEXT.md signature is `run_*(jumps, config, hamiltonian, trotter; kwargs...)`.
   - What's unclear: Should psi0/initial_dm be a required keyword? A positional argument? Have a default?
   - Recommendation: Make `psi0` a required keyword of `run_trajectory`. For `run_thermalize`, use a default of maximally mixed state I/d (matching the existing simulation scripts). This keeps positional signatures uniform: `run_*(jumps, config, hamiltonian, trotter; ...)`.

2. **What happens to construct_lindbladian?**
   - What we know: `construct_lindbladian` is a public export that constructs the full Liouvillian without spectral analysis. It's used in tests and benchmarks.
   - What's unclear: The phase scope says 4 entry points. Should `construct_lindbladian` remain as a lower-level utility?
   - Recommendation: Keep `construct_lindbladian` as an internal (or semi-public) utility. The 4 entry points are the public API; `construct_lindbladian` can remain for testing/benchmarking but is not one of the 4 primary entry points.

3. **LindbladResults eigenvalues field**
   - What we know: The current `run_lindbladian` uses Arpack `eigs(liouv, nev=2)` which returns only 2 eigenvalues. The locked decision says "spectral data only (eigenvalues, fixed_point, gap_mode, spectral_gap)".
   - What's unclear: Should `eigenvalues` store just the 2 from Arpack, or more?
   - Recommendation: Store whatever Arpack returns (currently 2). The field type `Vector{Complex{T}}` handles any count.

4. **Convergence observable_names in run_trajectory**
   - What we know: `run_trajectories_convergence` and `run_trajectories_adaptive` require `observable_names::Vector{String}` kwarg.
   - What's unclear: When `convergence=true` is passed to `run_trajectory`, should it auto-build preset observables or require the user to pass them?
   - Recommendation: Auto-build using `build_preset_trajectory_observables` when `convergence=true` and `observables` is `nothing`. If the user passes custom observables, require `observable_names` too.

## Sources

### Primary (HIGH confidence)
- `/Users/bence/code/QuantumFurnace.jl/src/results.jl` -- Current ExperimentResult serialization (Dict-based BSON pattern)
- `/Users/bence/code/QuantumFurnace.jl/src/structs.jl` -- All current type definitions, Config, Workspace, scratch structs
- `/Users/bence/code/QuantumFurnace.jl/src/furnace.jl` -- Current `run_lindbladian`, `run_thermalization`
- `/Users/bence/code/QuantumFurnace.jl/src/trajectories.jl` -- Current trajectory functions and workspace builder
- `/Users/bence/code/QuantumFurnace.jl/src/convergence.jl` -- Current convergence and adaptive trajectory functions
- `/Users/bence/code/QuantumFurnace.jl/src/krylov_eigsolve.jl` -- Current `krylov_spectral_gap`
- `/Users/bence/code/QuantumFurnace.jl/src/QuantumFurnace.jl` -- Module exports
- `/Users/bence/code/QuantumFurnace.jl/test/test_results.jl` -- ExperimentResult round-trip tests
- `/Users/bence/code/QuantumFurnace.jl/test/test_helpers.jl` -- Test fixtures and config factories
- `/Users/bence/code/QuantumFurnace.jl/simulations/` -- All 3 simulation scripts showing current call patterns

### Secondary (MEDIUM confidence)
- BSON.jl v0.3 behavior with Dict serialization (verified via existing codebase usage)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- All dependencies already in Project.toml, patterns already proven in codebase
- Architecture: HIGH -- Direct codebase analysis, all functions read and understood
- Pitfalls: HIGH -- All identified from actual codebase patterns and existing test failures

**Research date:** 2026-02-27
**Valid until:** 2026-03-27 (stable domain, internal refactor only)
