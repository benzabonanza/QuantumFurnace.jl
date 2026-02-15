# Architecture: Multi-Threaded Trajectory Engine, GNS Path, Adaptive Sampling, Experiments

**Domain:** Multi-threaded trajectory sampling + GNS trajectory path + convergence experiments for QuantumFurnace.jl
**Researched:** 2026-02-15
**Confidence:** HIGH (direct codebase analysis of all source files + Julia threading documentation)

## System Overview: New Components Integrated With Existing Architecture

```
EXISTING (unchanged)                        NEW / MODIFIED (this milestone)
=================                           ================================

+-- structs.jl ----------------------------+
|  HamHam{T}, TrottTrott{T}               |   READ-ONLY by all threads
|  Config{D,T}: ThermalizeConfig,          |
|    ThermalizeConfigGNS                   |
|  JumpOp, HotAlgorithmResults             |
+------------------------------------------+

+-- furnace_utensils.jl --------------------+
|  _precompute_data(config, ham_or_trott)  |   Called ONCE before spawning threads
|  NUFFTPrefactors, transition functions   |   Result is READ-ONLY shared data
+------------------------------------------+

+-- trajectories.jl ------------------------+   +-- NEW: parallel_trajectories.jl --------+
|  TrajectoryFramework{T,D}                |   |  run_trajectories_parallel(...)         |
|  build_trajectoryframework(...)          |   |    - Builds ONE framework               |
|  step_along_trajectory!(psi, fw)         |   |    - Clones per-thread workspace        |
|  run_trajectories(...)  [single-thread]  |   |    - Spawns @threads loop               |
|  _evolve_along_trajectory!(psi, fw, t)   |   |    - Reduces rho_mean + observables     |
+------------------------------------------+   |    - Returns TrajectoryResults          |
                                               +----------------------------------------+

+-- energy_domain.jl -----------------------+
|  pick_transition(ThermalizeConfigGNS)    |   Already dispatches correctly for GNS
|    -> _pick_transition_gns(config)       |   No changes needed
+------------------------------------------+

+-- coherent.jl ----------------------------+   +-- NEW: adaptive_sampling.jl -------------+
|  _precompute_coherent_total_B(...)       |   |  AdaptiveBatchManager                   |
|  with_coherent enforced false for GNS    |   |    - Runs batches of N_batch trajs      |
|  -> B term naturally skipped             |   |    - Evaluates convergence criteria      |
+------------------------------------------+   |    - Decides continue/stop               |
                                               |  Convergence criteria:                   |
                                               |    - Trace distance variance             |
                                               |    - Observable variance                 |
                                               +------------------------------------------+

                                               +-- NEW: experiment_runner.jl --------------+
                                               |  ExperimentSpec (parameter grid)         |
                                               |  ExperimentResult (config+data+metadata) |
                                               |  run_kms_vs_gns_experiment(...)          |
                                               |  save/load via JLD2                      |
                                               +------------------------------------------+
```

## Question 1: Thread-Safe Read-Only Access to Precomputed Data

### The Problem

`TrajectoryFramework` currently bundles read-only precomputed data (per-operator Kraus matrices, precomputed_data NamedTuple, config, jumps) with mutable per-trajectory workspace (`TrajectoryWorkspace` containing `jump_oft`, `psi_tmp`, `Rpsi`). Multiple threads cannot share a single `TrajectoryWorkspace` because `step_along_trajectory!` mutates it in-place.

### Recommended Architecture: Immutable Framework + Per-Thread Workspace

```julia
# EXISTING (no changes):
struct TrajectoryFramework{T,D<:AbstractDomain}
    domain::D
    jumps::Vector{JumpOp}
    ham_or_trott::Union{HamHam, TrottTrott}
    config::AbstractThermalizeConfig{D}
    precomputed_data::Any
    per_operator::Vector{PerOperatorKraus{T}}
    n_jumps::Int
    delta::Float64
    delta_eff::Float64
    alpha::Float64
    ws::TrajectoryWorkspace{T}  # <-- This is the problem for threading
end
```

The cleanest approach: **separate the mutable workspace from the immutable framework**.

```julia
# NEW: Thread-local workspace, allocated per thread
struct ThreadTrajectoryState{T}
    ws::TrajectoryWorkspace{T}   # jump_oft, psi_tmp, Rpsi buffers
    psi::Vector{T}               # Current state vector (reused across trajectories)
    rho_local::Matrix{T}         # Thread-local density matrix accumulator
    obs_local::Union{Nothing, Matrix{Float64}}  # Thread-local observable accumulator
    rng::Random.AbstractRNG      # Per-thread RNG for reproducibility
end

function ThreadTrajectoryState(::Type{T}, dim::Int, n_obs::Int, seed::Int) where {T}
    ThreadTrajectoryState{T}(
        TrajectoryWorkspace(T, dim),
        zeros(T, dim),
        zeros(T, dim, dim),
        n_obs > 0 ? zeros(Float64, n_obs, num_saves) : nothing,
        Random.MersenneTwister(seed),
    )
end
```

### Why This Works

1. **`TrajectoryFramework` fields are all immutable or read-only after construction:**
   - `per_operator::Vector{PerOperatorKraus{T}}` -- Kraus matrices R, K0, U_residual, U_B are precomputed and never mutated during stepping
   - `precomputed_data` -- NUFFTPrefactors, transition function, energy_labels are read-only
   - `jumps` -- JumpOp data and in_eigenbasis are read-only
   - `config` -- immutable struct
   - `ham_or_trott` -- immutable after construction

2. **`TrajectoryWorkspace` fields are the ONLY mutable state:**
   - `jump_oft` -- scratch buffer overwritten each step
   - `psi_tmp` -- scratch buffer overwritten each step
   - `Rpsi` -- scratch buffer overwritten each step

3. **Julia's memory model guarantees**: Read-only access to shared immutable data across threads requires no synchronization. The `NUFFTPrefactors.data` (3D array) and `per_operator` Kraus matrices are allocated once and only read during `step_along_trajectory!`.

### BLAS Thread Management

Critical: `step_along_trajectory!` calls `mul!` (BLAS gemv/gemm) inside each trajectory step. If Julia threads > 1 AND BLAS threads > 1, nested parallelism causes thread oversubscription and performance collapse.

```julia
# At the start of run_trajectories_parallel:
BLAS.set_num_threads(1)  # Disable BLAS multithreading
# ... run @threads loop ...
# Optionally restore after: BLAS.set_num_threads(n_blas_original)
```

This is standard practice for Julia multi-threaded scientific computing. The existing codebase already imports `Base.Threads` (line 18 of `QuantumFurnace.jl`).

### Per-Thread RNG Strategy

Julia's `TaskLocalRNG` (default since Julia 1.7) is per-task, not per-thread, which is suitable for `@spawn` but not ideal for `@threads` reproducibility. For reproducible parallel Monte Carlo:

```julia
# Create independent RNGs with jump-ahead (guaranteed non-overlapping streams)
rngs = [Random.MersenneTwister(seed + i) for i in 1:nthreads()]
```

Then pass `rngs[threadid()]` into the trajectory loop. The existing codebase uses `rand()` directly in `step_along_trajectory!`, so thread-local RNG requires either:

- **Option A (recommended):** Modify `step_along_trajectory!` to accept an `rng` argument: `rand(rng, 1:fw.n_jumps)` instead of `rand(1:fw.n_jumps)`. This is a small signature change with big payoff for reproducibility.
- **Option B:** Use `@threads :static` scheduling (guarantees thread IDs are stable) and seed task-local RNG per thread before the loop.

**Recommendation: Option A.** It is explicit, testable, and matches the pattern already used in `run_thermalization` (which already accepts `rng::AbstractRNG`).

## Question 2: GNS Trajectory Path Integration

### Current GNS Dispatch Architecture

The existing code already has a clean KMS-vs-GNS dispatch via config type:

```
Config type hierarchy:
  AbstractThermalizeConfig{D,T}
    |-- ThermalizeConfig{D,T}      (KMS line)
    |-- ThermalizeConfigGNS{D,T}   (GNS line)

Dispatch points (already working):
  1. pick_transition(config::ThermalizeConfigGNS)  -> _pick_transition_gns(config)
     Returns UNSHIFTED gamma_tilde(w) satisfying KMS condition directly
  2. _pick_alpha(config::ThermalizeConfigGNS)      -> _pick_alpha_gns(config)
     Returns GNS alpha function for BohrDomain
  3. validate_config!(config): enforces with_coherent=false for GNS
  4. _precompute_data(config, ham_or_trott): uses pick_transition(config),
     so GNS transition weight flows automatically via dispatch
  5. Coherent B term: _precompute_coherent_total_B returns nothing when
     config.with_coherent==false (enforced for GNS configs)
```

### What Needs to Change for GNS Trajectory Path

**Almost nothing.** The trajectory code in `trajectories.jl` dispatches on `AbstractThermalizeConfig{D}` for the domain type, and the GNS-specific behavior flows through `_precompute_data` and `pick_transition` which already handle `ThermalizeConfigGNS` correctly.

Here is the exhaustive list of functions that must accept `ThermalizeConfigGNS`:

| Function | File | Current Dispatch | GNS Status |
|----------|------|-----------------|------------|
| `run_trajectories()` | `trajectories.jl` | `config::AbstractThermalizeConfig` | **Already works** -- accepts any subtype |
| `build_trajectoryframework()` | `trajectories.jl` | `config::AbstractThermalizeConfig` | **Already works** |
| `_precompute_data()` | `furnace_utensils.jl` | Dispatches on domain `D` | **Already works** -- calls `pick_transition(config)` which dispatches on GNS |
| `pick_transition()` | `energy_domain.jl` | Has explicit `ThermalizeConfigGNS` method | **Already works** |
| `_precompute_R()` | `trajectories.jl` | `config::AbstractThermalizeConfig{D}` | **Already works** -- uses `precomputed_data.transition` |
| `step_along_trajectory!()` | `trajectories.jl` | `fw::TrajectoryFramework{T,D}` | **Already works** -- uses `fw.precomputed_data.transition` |
| `_precompute_coherent_total_B()` | `coherent.jl` | `config::AbstractConfig` | **Already works** -- returns `nothing` when `with_coherent==false` |
| `validate_config!()` | `misc_tools.jl` | `config::AbstractConfig` | **Already works** -- enforces `with_coherent=false` for GNS |
| `_select_b_plus_calculator()` | `furnace_utensils.jl` | `config::Union{LiouvConfig, ThermalizeConfig}` | **N/A** -- only called when `with_coherent==true`, which GNS configs cannot be |

### The One Required Change

The `_select_b_plus_calculator` function only dispatches on `Union{LiouvConfig, ThermalizeConfig}` (KMS types). If someone bypasses validation and passes a GNS config with `with_coherent=true` (which the inner constructor blocks), it would error. **This is NOT a bug -- the inner constructor prevents it.** No code change needed.

### Verification Strategy for GNS Trajectory Path

The right way to verify GNS trajectories work is:

1. Create a `ThermalizeConfigGNS` with valid parameters
2. Run `run_trajectories(jumps, config_gns, psi0, ham; ntraj=N)`
3. Compare trajectory-averaged rho against:
   - The GNS Liouvillian fixed point (from `run_lindbladian` with `LiouvConfigGNS`)
   - The Gibbs state (with expected domain approximation error)
4. Confirm that `per_op.U_B === nothing` for all operators (no coherent term)

```julia
# Smoke test for GNS trajectory path:
gns_config = ThermalizeConfigGNS(
    num_qubits=4, with_coherent=false,
    with_linear_combination=true, domain=EnergyDomain(),
    beta=10.0, sigma=0.1, a=1/30, b=0.4,
    num_energy_bits=12, w0=0.05, t0=T0,
    num_trotter_steps_per_t0=10,
    mixing_time=10.0, delta=0.01
)
result = run_trajectories(jumps, gns_config, psi0, ham; ntraj=1000)
# result.framework.per_operator[1].U_B === nothing  # Confirm no B term
```

## Question 3: Adaptive Sampling Architecture

### Design: Batch-Based Convergence Testing

Adaptive sampling means "keep running trajectories until convergence criteria are met." The architecture should separate three concerns:

1. **Batch execution** -- run N trajectories in parallel and accumulate results
2. **Convergence evaluation** -- assess whether current results are converged
3. **Decision logic** -- continue, stop, or adjust batch size

```julia
"""
Manages adaptive trajectory sampling with convergence detection.

Usage:
    manager = AdaptiveBatchManager(fw, psi0;
        batch_size=100, max_trajectories=100_000,
        convergence_criterion=TraceDistanceVariance(rtol=0.01))

    while !is_converged(manager)
        run_batch!(manager)
    end

    result = finalize(manager)
"""
mutable struct AdaptiveBatchManager{T, C<:ConvergenceCriterion}
    # Immutable setup
    fw::TrajectoryFramework{T}
    psi0::Vector{T}
    gibbs::Union{Nothing, Hermitian}  # Target state (if known, for convergence check)
    observables::Union{Nothing, Vector{Matrix{T}}}

    # Mutable state
    rho_running::Matrix{T}              # Running sum of |psi><psi| across all trajectories
    obs_running::Union{Nothing, Matrix{Float64}}  # Running sum of <O_i> per trajectory
    n_completed::Int                     # Total trajectories completed so far

    # Convergence tracking
    criterion::C
    convergence_history::Vector{Float64}  # Value of criterion after each batch
    batch_size::Int
    max_trajectories::Int

    # Thread-local state (created lazily on first batch)
    thread_states::Union{Nothing, Vector{ThreadTrajectoryState{T}}}
end
```

### Convergence Criteria (Trait Pattern)

```julia
abstract type ConvergenceCriterion end

"""Stop when trace distance between consecutive batch-averages is below rtol."""
struct TraceDistanceStability <: ConvergenceCriterion
    rtol::Float64       # Relative tolerance
    window::Int         # Number of consecutive batches that must be stable
end

"""Stop when standard error of per-observable means is below atol."""
struct ObservableVariance <: ConvergenceCriterion
    atol::Float64       # Absolute tolerance on standard error
end

"""Stop when batch-averaged rho is within threshold of target Gibbs state."""
struct GibbsProximity <: ConvergenceCriterion
    threshold::Float64  # Trace distance threshold
end

# Evaluate convergence after each batch
function evaluate_convergence(
    crit::TraceDistanceStability,
    manager::AdaptiveBatchManager,
)::Bool
    length(manager.convergence_history) < crit.window && return false
    recent = manager.convergence_history[end-crit.window+1:end]
    return all(v -> v < crit.rtol, recent)
end
```

### Thread Coordination for Adaptive Batches

Each batch runs a fixed number of trajectories in parallel, then synchronizes to evaluate convergence:

```
Batch 1:  Thread1[N/nthreads trajs] | Thread2[...] | ... | ThreadK[...]
          |                          |                     |
          +--- reduce into rho_running (atomic add or per-thread + merge) ---+
          |
          Evaluate convergence criterion
          |
Batch 2:  Thread1[...] | Thread2[...] | ... | ThreadK[...]
          |
          ...
          |
Done:     Return averaged rho_mean = rho_running / n_completed
```

The reduction after each batch can use a simple sequential merge of per-thread accumulators (cost: K * dim^2 additions, negligible compared to trajectory cost).

### Online Mean Computation

For numerical stability, use the Welford online algorithm for the running mean:

```julia
function accumulate_batch!(manager::AdaptiveBatchManager, batch_rho::Matrix)
    manager.n_completed += batch_size
    # Simple running sum (divide by n_completed at query time)
    manager.rho_running .+= batch_rho
end

function current_mean(manager::AdaptiveBatchManager)
    return manager.rho_running ./ manager.n_completed
end
```

For variance tracking (needed by `TraceDistanceStability`), track batch-level means:

```julia
# After each batch:
batch_mean = batch_rho ./ batch_size
overall_mean = current_mean(manager)
# Convergence metric: trace_distance(batch_mean, overall_mean)
push!(manager.convergence_history, trace_distance_h(
    Hermitian(batch_mean), Hermitian(overall_mean)))
```

## Question 4: Per-Observable Tracking

### During Trajectory vs Post-Average

**During trajectory (recommended).** Per-observable expectation values must be computed on the pure state `psi` within each trajectory step, then averaged. You cannot compute `<O>` from the averaged `rho_mean` after the fact if you want time-resolved observable curves, because the trajectory-average of `<psi|O|psi>` at intermediate times gives the time-resolved expectation, while `tr(O * rho_mean)` only gives the final-time expectation.

The existing `run_trajectories` already implements this correctly via `_accumulate_measurements!`:

```julia
# FROM trajectories.jl (existing, lines 405-437):
for _ in 1:ntraj
    copyto!(psi, psi0)
    # normalize once per trajectory
    save_idx = 1
    for step in 1:num_steps
        step_along_trajectory!(psi, fw)
        if step % save_every == 0
            save_idx += 1
            _accumulate_measurements!(mean_data, save_idx, psi, observables, tmp_meas)
        end
    end
end
mean_data ./= ntraj  # Average over trajectories
```

### Thread-Safe Observable Accumulation

Each thread accumulates into its own `obs_local::Matrix{Float64}` (dimensions: `n_observables x num_saves`). After the `@threads` loop, merge:

```julia
# Merge thread-local accumulators
for state in thread_states
    mean_data .+= state.obs_local
end
mean_data ./= ntraj
```

This avoids any atomic operations or locks during the hot loop. The merge cost is `O(nthreads * n_obs * num_saves)`, which is negligible.

### Observable Choice for KMS-vs-GNS Experiments

For the paper experiments, the relevant observables are:
1. **Trace distance to Gibbs**: `trace_distance_h(Hermitian(rho_traj), gibbs)` -- computed post-average from `rho_mean`
2. **Per-site magnetization**: `<Z_i>` for each qubit site -- measured per-trajectory at each save point
3. **Energy**: `<H>` -- measured per-trajectory, tracks thermalization
4. **Purity**: `tr(rho^2)` -- computed post-average from `rho_mean`

Items 2-3 require during-trajectory measurement. Items 1 and 4 can be computed from the final `rho_mean`.

## Question 5: Experiment Results Data Model

### Recommended Structure

```julia
"""
Complete record of a single experiment run: configuration + results + metadata.
Designed for serialization via JLD2.
"""
struct ExperimentResult{D, T<:AbstractFloat}
    # --- Configuration (what was run) ---
    config::AbstractThermalizeConfig{D,T}
    hamiltonian_id::String          # e.g. "heis_disordered_periodic_n4"
    db_type::Symbol                 # :KMS or :GNS
    num_qubits::Int
    beta::T
    domain::D

    # --- Trajectory results ---
    rho_mean::Matrix{Complex{T}}    # Trajectory-averaged density matrix
    n_trajectories::Int             # Number of trajectories used
    wall_time_seconds::Float64      # Total wall-clock time

    # --- Observable time series (if measured) ---
    times::Union{Nothing, Vector{Float64}}
    observable_means::Union{Nothing, Matrix{Float64}}  # (n_obs x n_times)
    observable_names::Union{Nothing, Vector{String}}

    # --- Convergence diagnostics ---
    trace_dist_to_gibbs::T          # Final trace distance to Gibbs
    trace_dist_to_fixed_point::Union{Nothing, T}  # If Liouvillian FP was computed
    convergence_history::Union{Nothing, Vector{Float64}}  # Batch convergence values

    # --- Metadata ---
    timestamp::String               # ISO 8601
    julia_version::String
    n_threads::Int
    git_commit::Union{Nothing, String}
    seed::Int
end
```

### Serialization Strategy

Use **JLD2** (already HDF5-compatible, pure Julia, handles arbitrary Julia types):

```julia
using JLD2

function save_experiment(result::ExperimentResult, path::String)
    jldsave(path; result=result)
end

function load_experiment(path::String)::ExperimentResult
    return load(path, "result")
end
```

JLD2 natively serializes complex Julia types including parametric structs, so `ExperimentResult{EnergyDomain, Float64}` round-trips correctly. The existing codebase uses BSON for Hamiltonian serialization (which has known fragility with struct evolution -- see `_load_hamiltonian_bson`), so JLD2 is a deliberate upgrade.

### File Organization for Experiments

```
experiments/
  kms_vs_gns/
    n4_beta10_energy/
      kms_ntraj10000_seed42.jld2
      gns_ntraj10000_seed42.jld2
    n6_beta10_energy/
      kms_ntraj10000_seed42.jld2
      gns_ntraj10000_seed42.jld2
    summary.jld2  # Aggregated comparison data
```

### Separation from Source

Experiment scripts should live in `experiments/` (not `src/`), with the new source components (`parallel_trajectories.jl`, `adaptive_sampling.jl`) in `src/`. The experiment runner uses the library API:

```
src/                        experiments/
  parallel_trajectories.jl    kms_vs_gns.jl  (script: defines grid, calls library)
  adaptive_sampling.jl        analyze_results.jl  (script: loads JLD2, makes plots)
  experiment_types.jl         plot_convergence.jl
```

## Question 6: Parameter Sweep Organization (n=4,6,8)

### Recommended: Flat Parameter Grid with Job Specification

```julia
"""
Specification for a single experiment point in a parameter sweep.
"""
struct ExperimentSpec{D<:AbstractDomain, T<:AbstractFloat}
    num_qubits::Int
    beta::T
    domain::D
    db_type::Symbol   # :KMS or :GNS
    with_coherent::Bool
    n_trajectories::Int
    seed::Int
    # ... other parameters from ThermalizeConfig
end

"""
Generate the full parameter grid for KMS-vs-GNS experiments.
"""
function kms_vs_gns_grid(;
    qubit_counts = [4, 6, 8],
    betas = [10.0],
    domains = [EnergyDomain()],
    db_types = [:KMS, :GNS],
    n_trajectories = 10_000,
    base_seed = 42,
)::Vector{ExperimentSpec}
    specs = ExperimentSpec[]
    seed_counter = base_seed

    for n in qubit_counts
        for beta in betas
            for domain in domains
                for db in db_types
                    with_coh = (db == :KMS)  # KMS can have B; GNS cannot
                    push!(specs, ExperimentSpec(
                        n, beta, domain, db, with_coh,
                        n_trajectories, seed_counter))
                    seed_counter += 1
                end
            end
        end
    end
    return specs
end
```

### Why Flat Grid, Not Nested Loops

A flat grid of `ExperimentSpec` objects has several advantages over nested `for n in [4,6,8]; for db in [:KMS,:GNS]; ...`:

1. **Resumability:** Each spec has a unique identity. If the experiment crashes at spec 7/12, you can resume from spec 8.
2. **Progress tracking:** `ProgressMeter` over a flat vector gives clear ETA.
3. **Selective re-running:** Failed or suspicious specs can be re-run individually.
4. **Parallelism flexibility:** The grid can be split across multiple Julia processes (Distributed) or run serially.

### Execution Pattern

```julia
function run_experiment_grid(specs::Vector{ExperimentSpec}, output_dir::String)
    results = ExperimentResult[]

    for (i, spec) in enumerate(specs)
        @info "Running experiment $i/$(length(specs))" spec.num_qubits spec.db_type

        # Build Hamiltonian (expensive for n=8, dim=256)
        ham = load_hamiltonian("heis", spec.num_qubits; beta=spec.beta)
        jumps = build_standard_jumps(ham, spec.num_qubits)
        psi0 = fill(ComplexF64(1.0), 2^spec.num_qubits) / sqrt(2^spec.num_qubits)

        # Build config (KMS or GNS)
        config = build_config(spec)

        # Run with parallel trajectories + adaptive sampling
        t0 = time()
        result_data = run_trajectories_parallel(
            jumps, config, psi0, ham;
            ntraj=spec.n_trajectories,
            seed=spec.seed,
        )
        wall_time = time() - t0

        # Package result
        result = ExperimentResult(
            config=config,
            hamiltonian_id="heis_disordered_periodic_n$(spec.num_qubits)",
            db_type=spec.db_type,
            # ... fill remaining fields ...
        )

        # Save immediately (crash resilience)
        save_experiment(result, joinpath(output_dir, filename(spec)))
        push!(results, result)
    end
    return results
end
```

### Scaling Considerations for n=4,6,8

| n | dim | Precompute cost | Per-trajectory step cost | Memory per thread |
|---|-----|----------------|------------------------|-------------------|
| 4 | 16  | ~1s (NUFFT)    | ~100us (12 jumps, ~200 energies) | ~50 KB |
| 6 | 64  | ~10s (NUFFT)   | ~1ms (18 jumps, ~200 energies) | ~800 KB |
| 8 | 256 | ~100s (NUFFT)  | ~20ms (24 jumps, ~200 energies) | ~12 MB |

For n=8, each trajectory step involves `mul!` on 256x256 matrices. With 10K trajectories at 1000 steps each, that is 10M steps * 20ms = ~55 hours single-threaded. With 8 threads: ~7 hours. With adaptive sampling stopping early: potentially much less.

**Implication:** n=8 experiments should:
1. Use adaptive sampling aggressively (stop when converged)
2. Use the largest reasonable batch size (amortize reduction overhead)
3. Consider fewer trajectories if convergence is fast
4. Save intermediate results after each batch

## Component Boundaries: New vs Modified

### New Files

| File | Purpose | Key Types/Functions |
|------|---------|-------------------|
| `src/parallel_trajectories.jl` | Multi-threaded trajectory runner | `ThreadTrajectoryState`, `run_trajectories_parallel`, `step_along_trajectory!(psi, fw, rng)` |
| `src/adaptive_sampling.jl` | Adaptive batch manager + convergence criteria | `AdaptiveBatchManager`, `ConvergenceCriterion`, `run_batch!`, `is_converged` |
| `src/experiment_types.jl` | Data model for experiment results | `ExperimentResult`, `ExperimentSpec`, `save_experiment`, `load_experiment` |
| `experiments/kms_vs_gns.jl` | Experiment script (not in module) | `kms_vs_gns_grid`, `run_experiment_grid` |

### Modified Files

| File | What Changes | Why |
|------|-------------|-----|
| `src/trajectories.jl` | Add `rng::AbstractRNG` parameter to `step_along_trajectory!` (default: `Random.default_rng()`) | Thread-safe RNG for parallel execution |
| `src/QuantumFurnace.jl` | Add `include` for new files, export new public API | Module registration |
| `Project.toml` | Add `JLD2` dependency | Experiment result serialization |

### Unchanged Files

Every other source file remains unchanged. The GNS trajectory path works through existing dispatch without modification. The key insight is that `ThermalizeConfigGNS <: AbstractThermalizeConfig`, so all trajectory functions that dispatch on `AbstractThermalizeConfig` already accept GNS configs. The GNS-specific behavior (unshifted transition, no B term) is delivered by:
1. `pick_transition(config::ThermalizeConfigGNS)` returning the unshifted gamma
2. `config.with_coherent == false` causing B term precomputation to be skipped

## Data Flow: End-to-End for a KMS-vs-GNS Experiment

```
[ExperimentSpec(n=4, db=:KMS)]
    |
    v
[load_hamiltonian("heis", 4; beta=10.0)]  -->  HamHam{Float64}
    |
    v
[build_standard_jumps(ham, 4)]  -->  Vector{JumpOp} (12 jumps)
    |
    v
[build_config(spec)]  -->  ThermalizeConfig{EnergyDomain, Float64}
    |                       (with_coherent=true for KMS)
    v
[_precompute_data(config, ham)]  -->  NamedTuple (transition, energy_labels, gamma_norm_factor)
    |
    v
[build_trajectoryframework(jumps, ham, config, precomputed, scratch, delta)]
    |                       -->  TrajectoryFramework (per_operator with U_B matrices)
    v
[run_trajectories_parallel(jumps, config, psi0, ham; ntraj=10000)]
    |
    |  BLAS.set_num_threads(1)
    |
    |  Create nthreads() ThreadTrajectoryState instances
    |
    |  @threads for batch in 1:ntraj
    |      state = thread_states[threadid()]
    |      copyto!(state.psi, psi0)
    |      _evolve_along_trajectory!(state.psi, fw, total_time, state.rng)
    |      _accumulate_density_matrix!(state.rho_local, state.psi)
    |      (optional: _accumulate_measurements!(...))
    |  end
    |
    |  Reduce: rho_mean = sum(state.rho_local for state in thread_states) / ntraj
    |
    v
[ExperimentResult(config, rho_mean, trace_dist_to_gibbs, ...)]
    |
    v
[save_experiment(result, "kms_n4_beta10.jld2")]


[ExperimentSpec(n=4, db=:GNS)]
    |
    v
... same pipeline but with ThermalizeConfigGNS ...
    |
    |  pick_transition dispatches to _pick_transition_gns (unshifted gamma)
    |  config.with_coherent=false -> per_op.U_B === nothing (no B term)
    |
    v
[ExperimentResult(config_gns, rho_mean_gns, ...)]
    |
    v
[save_experiment(result, "gns_n4_beta10.jld2")]


COMPARISON:
    load both results
    compare trace_dist_to_gibbs: KMS vs GNS
    compare observable time series
    compare convergence rates
```

## Suggested Build Order

Based on dependency analysis:

### Phase 1: RNG Signature Change (enables all subsequent work)

**Modify:** `step_along_trajectory!` to accept `rng::AbstractRNG = Random.default_rng()`.

This is a 2-line change per domain variant (replace `rand()` -> `rand(rng)`, `rand(1:fw.n_jumps)` -> `rand(rng, 1:fw.n_jumps)`). No behavioral change for existing callers.

**Files modified:** `trajectories.jl`
**Depends on:** Nothing
**Enables:** Phase 2 (parallel execution)

### Phase 2: Multi-Threaded Trajectory Engine

**New:** `parallel_trajectories.jl` with `ThreadTrajectoryState` and `run_trajectories_parallel`.

**Files new:** `src/parallel_trajectories.jl`
**Depends on:** Phase 1 (RNG parameter)
**Enables:** Phase 4 (experiments)

### Phase 3: GNS Trajectory Verification

**New:** Tests confirming GNS trajectory path works. No source changes needed.

**Files new:** `test/test_gns_trajectory.jl` or section in existing test file
**Depends on:** Nothing (GNS dispatch already works)
**Enables:** Phase 4 (experiments need both KMS and GNS)

### Phase 4: Experiment Data Model + Serialization

**New:** `experiment_types.jl` with `ExperimentResult`, `ExperimentSpec`, save/load.

**Files new:** `src/experiment_types.jl`
**Depends on:** Nothing (pure data types)
**Enables:** Phase 5 (adaptive sampling), Phase 6 (experiments)

### Phase 5: Adaptive Sampling

**New:** `adaptive_sampling.jl` with `AdaptiveBatchManager`, convergence criteria.

**Files new:** `src/adaptive_sampling.jl`
**Depends on:** Phase 2 (parallel trajectories), Phase 4 (result types)
**Enables:** Phase 6 (experiments use adaptive sampling)

### Phase 6: KMS-vs-GNS Experiment Runner

**New:** `experiments/kms_vs_gns.jl` script.

**Files new:** `experiments/kms_vs_gns.jl`, `experiments/analyze_results.jl`
**Depends on:** All previous phases
**Enables:** Paper results

## Anti-Patterns to Avoid

### Anti-Pattern: Shared Mutable State in @threads Loop

**What goes wrong:** Two threads write to the same `rho_mean` accumulator matrix without synchronization. Data race produces incorrect results silently (no crash, just wrong numbers).

**Prevention:** Each thread has its own `rho_local` accumulator. Merge after the loop completes. Never use `@threads` with a shared mutable accumulator unless using `Threads.Atomic` or locks.

### Anti-Pattern: Global RNG in Parallel Code

**What goes wrong:** All threads share `Random.default_rng()`, producing correlated random streams and non-reproducible results. In Julia >= 1.7, TaskLocalRNG is per-task, but `@threads` tasks may share RNG state depending on scheduler.

**Prevention:** Explicitly construct per-thread RNGs with independent seeds. Pass as argument to `step_along_trajectory!`.

### Anti-Pattern: Saving Results Only at the End

**What goes wrong:** A 10-hour n=8 experiment crashes at hour 9. All results lost.

**Prevention:** Save after each experiment point. Use `jldsave` with unique filenames per spec. Implement resume logic (check if output file exists, skip if so).

### Anti-Pattern: Monolithic Experiment Script

**What goes wrong:** Experiment logic, trajectory engine, data types, and analysis all in one 500-line script. Cannot test components independently. Cannot reuse trajectory engine without experiment boilerplate.

**Prevention:** Library code in `src/` (importable, testable), experiment scripts in `experiments/` (use the library).

## Sources

- Direct codebase analysis of all 23 source files in QuantumFurnace.jl `src/`
- Direct analysis of all test files in `test/`
- [Julia Multi-Threading Documentation](https://docs.julialang.org/en/v1/manual/multi-threading/)
- [Julia for HPC: Multithreading](https://enccs.github.io/julia-for-hpc/multithreading/) -- BLAS thread management patterns
- [JLD2.jl Documentation](https://juliaio.github.io/JLD2.jl/stable/) -- struct serialization, HDF5 compatibility
- [Julia Random Numbers Documentation](https://docs.julialang.org/en/v1/stdlib/Random/) -- TaskLocalRNG, per-thread RNG
- [Reproducible Multithreaded Monte Carlo Discussion](https://discourse.julialang.org/t/reproducible-multithreaded-monte-carlo-task-local-random/35269) -- per-task RNG patterns
- [BLAS Thread Overhead Issue](https://github.com/JuliaLang/julia/issues/43292) -- BLAS.set_num_threads(1) for multi-threaded code
- [ITensors.jl Multithreading Guide](https://itensor.github.io/ITensors.jl/dev/Multithreading.html) -- BLAS thread disabling pattern
- [Carlo.jl Framework (2025)](https://arxiv.org/abs/2408.03386) -- Monte Carlo simulation framework design patterns
- [DifferentialEquations.jl Parallel Monte Carlo](https://docs.sciml.ai/release-4.6/features/monte_carlo.html) -- multi-threaded trajectory patterns

---
*Architecture research for: Multi-threaded trajectory engine + GNS path + adaptive sampling + experiments (v1.1 milestone)*
*Researched: 2026-02-15*
