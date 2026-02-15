# Stack Research: v1.1 Multi-Threaded Trajectory Engine & KMS-vs-GNS Experiments

**Domain:** Multi-threaded Monte Carlo trajectory sampling, adaptive convergence, experiment data management for quantum Lindbladian simulation
**Researched:** 2026-02-15
**Confidence:** HIGH (core recommendations use Julia stdlib threading + verified pure-Julia packages; one MEDIUM item noted)

## Scope

This stack research covers ONLY the additions needed for the v1.1 milestone:
1. Multi-threaded trajectory sampling with shared read-only precomputed data
2. Thread-safe RNG patterns (per-task seeding for reproducibility)
3. Adaptive convergence-based early stopping per observable
4. Experiment data serialization (parameter sweeps, convergence curves)
5. BLAS thread management for trajectory parallelism

It does NOT re-research the existing stack (Arpack, FINUFFT, LinearAlgebra, BSON, StableRNGs, HypothesisTests, etc.). See the v1.0 STACK.md for those.

## Current Stack (Relevant Subset)

These existing dependencies are directly used and need no changes:

| Dependency | Role in This Milestone | Status |
|------------|----------------------|--------|
| LinearAlgebra (stdlib) | mul!, dot, BLAS for matrix ops in trajectory steps | Keep -- must manage BLAS.set_num_threads(1) |
| Random (stdlib) | rand() calls in step_along_trajectory! via TaskLocalRNG | Keep -- Julia's default RNG is already task-local since 1.7 |
| Base.Threads (stdlib) | Already imported in QuantumFurnace.jl | Keep -- extend usage |
| StableRNGs (test extra) | Reproducible per-trajectory seeding in tests | Keep in test target |
| Statistics (stdlib) | mean(), std() for convergence monitoring | Keep -- use more heavily |
| BSON | Hamiltonian loading (load_hamiltonian) | Keep as-is, not for new data |

## New Production Dependencies (src/)

### Core: Multi-Threading Infrastructure

| Technology | Version | Purpose | Why Recommended |
|------------|---------|---------|-----------------|
| *No new package needed* | -- | Julia threading via `Threads.@spawn` + manual reduction | The existing `using Base.Threads` in QuantumFurnace.jl is sufficient. The trajectory parallelism pattern is embarrassingly parallel: spawn N_traj independent tasks, each with its own workspace clone, accumulate results via per-thread partial sums merged after `@sync`. This does NOT need OhMyThreads.jl -- see Alternatives Considered below. |

**Confidence: HIGH** -- Julia's task-based threading with `@spawn`/`@sync` is the standard approach for embarrassingly parallel Monte Carlo. Confirmed in [Julia multi-threading docs](https://docs.julialang.org/en/v1/manual/multi-threading/) and used by QuantumToolbox.jl, SequentialMonteCarlo.jl, and MCIntegration.jl.

### Core: Data Serialization

| Technology | Version | Purpose | Why Recommended |
|------------|---------|---------|-----------------|
| JLD2.jl | >= 0.5 | Save/load experiment results (convergence curves, parameter sweeps, density matrices) | Pure Julia, no C library dependency. Preserves Julia types exactly -- critical for saving `NamedTuple` results, `Matrix{ComplexF64}` density matrices, and nested experiment metadata without manual serialization. HDF5-compatible format means results can also be read from Python/MATLAB if needed for paper figures. Already the de facto standard for Julia scientific computing data persistence. Replaces any ad-hoc BSON usage for experiment data. |

**Confidence: HIGH** -- JLD2 v0.6.3 (latest release Nov 2024) verified via [GitHub releases](https://github.com/JuliaIO/JLD2.jl/releases). API confirmed via [official docs](https://juliaio.github.io/JLD2.jl/stable/). Pure Julia, actively maintained by JuliaIO org.

### Supporting: Experiment Result Tables

| Technology | Version | Purpose | Why Recommended |
|------------|---------|---------|-----------------|
| DataFrames.jl | >= 1.7 | Tabular experiment results (parameter sweep tables, convergence summaries) | For the KMS-vs-GNS paper experiments, you need to sweep over (n_qubits, beta, sigma, domain, kms_vs_gns) and record metrics (trace_distance, spectral_gap, convergence_time, n_traj_to_converge). A DataFrame is the natural structure for this. Enables filtering, grouping, and CSV export for paper tables. |
| CSV.jl | >= 0.10 | Export experiment tables to CSV for paper figures/LaTeX | Companion to DataFrames. Write `CSV.write("results.csv", df)` for direct consumption by plotting tools or LaTeX table generators. Lightweight. |

**Confidence: HIGH** -- DataFrames.jl is the standard tabular data package in Julia, verified via [official docs](https://dataframes.juliadata.org/stable/). CSV.jl is its standard I/O companion.

## No New Production Dependencies Needed for Threading

This is the critical insight: **Julia's stdlib threading is sufficient for this milestone.**

### Why stdlib Threading Is Enough

The trajectory sampling loop in `run_trajectories` (trajectories.jl:383-395) is embarrassingly parallel:

```julia
# Current sequential loop:
@inbounds for trajectory in 1:ntraj
    copyto!(psi, psi0)
    _evolve_along_trajectory!(psi, fw, total_time)
    _accumulate_density_matrix!(rho_mean, psi)
end
```

The parallel version needs:
1. **Per-task workspace clones** (TrajectoryWorkspace + psi vector) -- cheap to allocate
2. **Per-task partial rho accumulation** -- merge after all tasks complete
3. **Shared read-only data** -- TrajectoryFramework fields (per_operator, precomputed_data, config) are immutable after construction

This maps directly to Julia's `@sync`/`@spawn` pattern with per-task buffers and a final reduction. No locks needed because each task writes only to its own buffers.

### Threading Pattern for run_trajectories

```julia
function run_trajectories_threaded(fw, psi0, ntraj; nchunks=Threads.nthreads())
    dim = length(psi0)
    CT = eltype(psi0)

    # Divide trajectories into chunks
    chunk_sizes = _divide_work(ntraj, nchunks)

    # Per-chunk partial results (allocated before spawning)
    partial_rhos = [zeros(CT, dim, dim) for _ in 1:nchunks]

    @sync for (chunk_id, chunk_size) in enumerate(chunk_sizes)
        Threads.@spawn begin
            # Per-task workspace (NOT shared)
            ws = TrajectoryWorkspace(CT, dim)
            psi = Vector{CT}(undef, dim)

            # Per-task RNG: Julia's TaskLocalRNG is automatically task-local
            # Seed deterministically from chunk_id for reproducibility
            Random.seed!(fw.base_seed + chunk_id)

            local_rho = partial_rhos[chunk_id]
            for _ in 1:chunk_size
                copyto!(psi, psi0)
                _evolve_along_trajectory!(psi, fw_with_ws(fw, ws), total_time)
                _accumulate_density_matrix!(local_rho, psi)
            end
        end
    end

    # Merge: sum partial rhos (sequential, cheap)
    rho_mean = sum(partial_rhos) ./ ntraj
    hermitianize!(rho_mean)
    return rho_mean
end
```

### Why TrajectoryFramework Needs Refactoring

The current `TrajectoryFramework` struct embeds a mutable `ws::TrajectoryWorkspace{T}` field. This is the **single blocker** for thread safety -- the framework is otherwise read-only. The fix is straightforward:

1. Remove `ws` from `TrajectoryFramework` (make it a separate argument)
2. Each spawned task creates its own `TrajectoryWorkspace`
3. `step_along_trajectory!(psi, fw, ws)` takes workspace as explicit argument

This is a small refactor (adding one argument to 2-3 functions) with large payoff.

## Thread-Safe RNG Strategy

### Julia's TaskLocalRNG (No Package Needed)

Since Julia 1.7, the default RNG (`Random.default_rng()`) returns a `TaskLocalRNG` -- each `@spawn`ed task gets its own RNG state, deterministically seeded from the parent task's state. This means:

- `rand()` inside a spawned task is automatically thread-safe
- No data races on RNG state between tasks
- Deterministic if the parent task's RNG is seeded before spawning

**Confidence: HIGH** -- Verified in [Julia Random docs](https://docs.julialang.org/en/v1/stdlib/Random/) and confirmed by [JuliaCon 2025 talk on TaskLocalRNG](https://pretalx.com/juliacon-2025/talk/ZNBEAN/).

### Reproducibility Pattern

For reproducible multi-threaded experiments:

```julia
# Seed parent task before spawning
Random.seed!(experiment_seed)

# Each @spawn inherits a deterministic child RNG state
# Results are reproducible IF trajectory count per chunk is fixed
@sync for chunk_id in 1:nchunks
    Threads.@spawn begin
        # TaskLocalRNG is already seeded deterministically
        for _ in 1:chunk_size
            r = rand()  # thread-safe, deterministic
        end
    end
end
```

**Caveat:** Reproducibility requires fixed `nchunks` (i.e., fixed thread count). Changing `julia -t 4` to `julia -t 8` will change results because task seeding depends on spawn order. This is acceptable for experiments -- document the thread count alongside results.

### StableRNGs for Cross-Version Test Reproducibility

For regression tests that must produce identical results across Julia versions, continue using `StableRNG(seed)` per task:

```julia
Threads.@spawn begin
    rng = StableRNG(base_seed + chunk_id)
    # Pass rng explicitly or copy! to task-local RNG
end
```

This is test-only. Production code uses `TaskLocalRNG` (zero overhead, no extra package).

## BLAS Thread Management

### The Problem

`step_along_trajectory!` calls `mul!(ws.psi_tmp, per_op.K0, psi)` -- matrix-vector products via BLAS. By default, BLAS uses multiple threads per `mul!` call. When Julia also uses multiple threads for trajectory parallelism, you get thread oversubscription: N_julia_threads x N_blas_threads total threads competing for N_cores physical cores.

### The Solution

Set `BLAS.set_num_threads(1)` before launching parallel trajectories. The system sizes in QuantumFurnace (dim = 2^n for n=4,6,8,12, i.e., dim = 16, 64, 256, 4096) are small enough that single-threaded BLAS is optimal:

- **dim <= 256 (n <= 8):** BLAS parallelism has zero benefit. Overhead of thread management exceeds gains for matrices this small.
- **dim = 4096 (n = 12):** Marginal benefit from BLAS threads, but trajectory-level parallelism is far more efficient (thousands of independent trajectories vs. parallelizing a single 4096x4096 multiply).

```julia
function run_trajectories_threaded(...)
    old_blas_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)  # Prevent oversubscription
    try
        # ... spawn trajectory tasks ...
    finally
        BLAS.set_num_threads(old_blas_threads)  # Restore
    end
end
```

**Confidence: HIGH** -- This pattern is used by [ITensors.jl](https://itensor.github.io/ITensors.jl/dev/Multithreading.html) and recommended in [Julia for HPC course](https://enccs.github.io/julia-for-hpc/multithreading/). The [Julia Discourse thread on BLAS vs Julia threads](https://discourse.julialang.org/t/julia-threads-vs-blas-threads/8914) confirms this is the standard approach.

## Adaptive Convergence Monitoring

### No External Package Needed

The convergence monitoring for adaptive trajectory count requires tracking running mean and variance of per-observable expectation values. This is a Welford online algorithm -- 10 lines of code, no package dependency:

```julia
mutable struct OnlineMeanVar{T}
    n::Int
    mean::T
    M2::T  # sum of squared deviations
end

function update!(s::OnlineMeanVar, x)
    s.n += 1
    delta = x - s.mean
    s.mean += delta / s.n
    delta2 = x - s.mean
    s.M2 += delta * delta2
end

variance(s::OnlineMeanVar) = s.n > 1 ? s.M2 / (s.n - 1) : zero(s.M2)
stderr(s::OnlineMeanVar) = sqrt(variance(s) / s.n)
```

For adaptive stopping: check `stderr(stat) / abs(stat.mean) < rtol` every `check_interval` trajectories.

### Why NOT OnlineStats.jl

OnlineStats.jl provides Mean, Variance, and many other streaming statistics. However:
- It is a heavyweight dependency (many transitive deps) for something that is literally 10 lines of code
- The package is designed for general-purpose streaming analytics, not scientific convergence monitoring
- You need per-observable tracking (a vector of `OnlineMeanVar` structs), which is trivial to implement but awkward with OnlineStats types

**Decision:** Implement Welford's algorithm inline. No package needed.

**Confidence: HIGH** -- Welford's algorithm is numerically stable and well-understood. Used by NumPy, Rust's statistical crates, and OnlineStats.jl itself internally.

## What NOT to Add for This Milestone

| Avoid | Why | What to Do Instead |
|-------|-----|-------------------|
| OhMyThreads.jl | Adds a dependency for syntactic sugar. The `tmapreduce` pattern is elegant but the trajectory parallelism here is simple enough that `@sync`/`@spawn` with manual chunk allocation is clearer and has zero dependency cost. OhMyThreads' `TaskLocalValue` pattern is useful but `TrajectoryWorkspace` allocation in the spawn body achieves the same thing. | Use `@sync`/`@spawn` with per-task workspace allocation. Consider OhMyThreads only if you later add more complex parallel patterns (e.g., dynamic load balancing for variable-length trajectories). |
| Distributed.jl (multi-process) | Already in `[deps]` but using it for trajectory parallelism adds complexity (serialization of frameworks, separate memory spaces) for zero benefit on a single node. Multi-threading is strictly better for shared-memory trajectory parallelism because the precomputed data (NUFFT prefactors, per-operator Kraus matrices) is read-only and shared without copying. | Use `Threads.@spawn` for single-node parallelism. Distributed is relevant only for multi-node cluster runs (separate milestone). |
| OnlineStats.jl | 10 lines of Welford's algorithm replaces the entire package for this use case. OnlineStats brings transitive deps (OrderedCollections, AbstractTrees, etc.) for features you will never use. | Implement `OnlineMeanVar` struct inline. |
| DrWatson.jl | Excellent for managing large experiment projects with file naming conventions, git tagging, etc. However, QuantumFurnace is a library package, not an experiment-management project. The experiment scripts that use QuantumFurnace should use DrWatson; the package itself should not depend on it. | Let experiment scripts (outside the package) use DrWatson if desired. The package provides `save_experiment_results(path, results)` using JLD2 directly. |
| HDF5.jl | JLD2 already produces HDF5-compatible files without requiring the HDF5 C library. Adding HDF5.jl brings a binary dependency (libhdf5) that complicates installation, especially on clusters. | Use JLD2, which is pure Julia and HDF5-compatible. If raw HDF5 is needed later (e.g., for interop with specific tools), add it then. |
| Arrow.jl | Faster than JLD2 for flat tabular data but cannot serialize Julia-specific types like `Matrix{ComplexF64}` or nested NamedTuples. The experiment results mix tabular (parameter sweep) and non-tabular (density matrices, convergence curves) data. | Use JLD2 for full results (types preserved), CSV for human-readable tables. |
| Transducers.jl / FLoops.jl | Predecessor of OhMyThreads with heavier API. FLoops is in maintenance mode. Transducers API is powerful but overkill for this use case. | Use `@sync`/`@spawn`. |
| SharedArrays (stdlib) | Already in `[deps]` for NUFFT prefactors with Distributed. For multi-threaded (not multi-process) code, regular `Array` is shared between threads automatically (same memory space). SharedArrays adds overhead for no benefit in the threading context. | Regular Julia `Array` for shared read-only data in threaded code. |
| Atomics / lock-based accumulation | Julia `Threads.Atomic` only works on primitive types (Float64, Int64), not matrices. Lock-based matrix accumulation adds contention. Per-task partial sums with final merge is both simpler and faster. | Per-task partial rho accumulation, merge after `@sync`. |

## Recommended Stack Summary

### Production Dependencies to ADD to `[deps]`

| Package | UUID | Purpose |
|---------|------|---------|
| JLD2 | `033835bb-8acc-5ee8-8aed-2f6b5d7c090b` | Experiment result serialization |

### Production Dependencies to ADD to `[deps]` (Experiment Scripts)

| Package | UUID | Purpose |
|---------|------|---------|
| DataFrames | `a93c6f00-e57d-5684-b7b6-d8193f3e46c0` | Parameter sweep result tables |
| CSV | `336ed68f-0bac-5ca0-87d4-7b16caf5d00b` | CSV export for paper tables/figures |

**Note:** DataFrames and CSV could alternatively live only in experiment scripts outside the package. If the package itself provides a `sweep_parameters(...)` API that returns a DataFrame, they belong in `[deps]`. If experiments are standalone scripts, they go in the script's own environment.

### No New Test Dependencies

The existing test extras (StableRNGs, HypothesisTests, StatsBase, Aqua) are sufficient for testing the threaded trajectory engine.

## Installation

```julia
using Pkg
Pkg.activate(".")

# Add production dependency
Pkg.add("JLD2")

# For experiment scripts (decide if in-package or script-level)
Pkg.add(["DataFrames", "CSV"])
```

Concrete `Project.toml` additions:

```toml
[deps]
# ... existing deps ...
JLD2 = "033835bb-8acc-5ee8-8aed-2f6b5d7c090b"

# Optional: only if sweep API is part of the package
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"

[compat]
# ... existing compat ...
JLD2 = "0.5, 0.6"
DataFrames = "1"
CSV = "0.10"
```

## How Each Addition Integrates with Existing Code

### JLD2 Integration

```julia
using JLD2

# Save full experiment result
function save_experiment(path::String, result::NamedTuple)
    jldsave(path;
        rho_mean = result.rho_mean,
        measurements = result.measurements_mean,
        times = result.times,
        config = result.config_summary,   # Dict, not the full config (avoid type issues across versions)
        metadata = Dict(
            "n_qubits" => config.num_qubits,
            "ntraj" => ntraj,
            "julia_version" => string(VERSION),
            "nthreads" => Threads.nthreads(),
            "timestamp" => Dates.now(),
        )
    )
end

# Load back
data = load(path)
rho = data["rho_mean"]  # Matrix{ComplexF64} -- types preserved exactly
```

### Threading Integration with Existing TrajectoryFramework

The key refactor: separate workspace from framework.

```julia
# BEFORE (current): workspace embedded in framework
struct TrajectoryFramework{T,D}
    # ... read-only fields ...
    ws::TrajectoryWorkspace{T}  # MUTABLE -- prevents sharing
end

# AFTER: workspace is a separate argument
struct TrajectoryFramework{T,D}
    # ... read-only fields only ...
    # ws removed
end

# step_along_trajectory! gains a ws argument:
function step_along_trajectory!(psi, fw, ws)
    # Uses fw (shared read-only) and ws (per-task mutable)
end

# Backward compatibility wrapper:
function step_along_trajectory!(psi, fw_with_ws::TrajectoryFrameworkLegacy)
    step_along_trajectory!(psi, fw_with_ws.framework, fw_with_ws.ws)
end
```

### BLAS Thread Control Integration

```julia
# In the threaded trajectory runner:
function run_trajectories(jumps, config, psi0, hamiltonian;
                          ntraj=1, nthreads=Threads.nthreads(), ...)
    # ... build framework ...

    if ntraj > 1 && nthreads > 1
        _run_trajectories_threaded(fw, psi0, ntraj, nthreads; ...)
    else
        _run_trajectories_sequential(fw, psi0, ntraj; ...)
    end
end

function _run_trajectories_threaded(fw, psi0, ntraj, nthreads; ...)
    old_blas = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        # ... @sync/@spawn pattern ...
    finally
        BLAS.set_num_threads(old_blas)
    end
end
```

### Adaptive Convergence Integration

```julia
function run_trajectories_adaptive(fw, psi0, observables;
                                    rtol=1e-3, max_traj=100_000,
                                    check_every=100)
    n_obs = length(observables)
    stats = [OnlineMeanVar(0, 0.0, 0.0) for _ in 1:n_obs]
    rho_acc = zeros(ComplexF64, dim, dim)

    ntraj = 0
    while ntraj < max_traj
        # Run a batch of trajectories (potentially threaded)
        for _ in 1:check_every
            psi = copy(psi0)
            _evolve_along_trajectory!(psi, fw, ws, total_time)
            _accumulate_density_matrix!(rho_acc, psi)
            for i in 1:n_obs
                val = real(dot(psi, observables[i] * psi))
                update!(stats[i], val)
            end
            ntraj += 1
        end

        # Check convergence: all observables within tolerance
        converged = all(i -> stderr(stats[i]) / max(abs(stats[i].mean), 1e-15) < rtol,
                       1:n_obs)
        converged && break
    end

    rho_mean = rho_acc ./ ntraj
    return (rho_mean=rho_mean, stats=stats, ntraj=ntraj, converged=converged)
end
```

## Alternatives Considered

| Recommended | Alternative | When to Use Alternative |
|-------------|-------------|-------------------------|
| `@sync`/`@spawn` manual chunking | OhMyThreads.jl `tmapreduce` | If you add dynamic load balancing (e.g., early-stopped trajectories where some finish faster). OhMyThreads' `:greedy` scheduler handles uneven workloads well. For fixed-length trajectories (current case), manual chunking is simpler. |
| JLD2 for full results | BSON (already in deps) | Never for new data. BSON is used only for legacy Hamiltonian loading. JLD2 is faster, preserves types better, and produces HDF5-compatible files. |
| JLD2 for full results | HDF5.jl | Only if you need HDF5 features that JLD2 doesn't expose (chunked datasets, compression filters, MPI-IO). For experiment results at this scale, JLD2 is sufficient and avoids the C library dependency. |
| JLD2 + CSV for experiment data | Arrow.jl | Only if you need columnar format for very large tabular datasets (millions of rows). For parameter sweep tables with hundreds of rows, CSV is simpler and human-readable. |
| Inline Welford algorithm | OnlineStats.jl | Only if you need dozens of different streaming statistics (quantiles, histograms, etc.). For mean + variance + standard error, inline code is better. |
| Manual chunk allocation | FLoops.jl / Transducers.jl | Never. FLoops is in maintenance mode. Transducers is powerful but the learning curve is not justified for this use case. |
| TaskLocalRNG (Julia stdlib) | Per-task StableRNG instances | Only in tests where cross-Julia-version reproducibility is required. For production runs, TaskLocalRNG is zero-overhead and deterministic within a Julia version. |
| DataFrames + CSV | Printf-based table output | If you want zero dependencies and your experiment scripts are standalone. Printf tables work but cannot be filtered, sorted, or joined programmatically. |

## Stack Patterns by Use Case

**If running experiments on a single workstation (n=4,6,8):**
- Use `julia -t auto` (all cores)
- Set `BLAS.set_num_threads(1)`
- JLD2 for results, CSV for summary tables
- Adaptive convergence with rtol=1e-3 to avoid over-sampling

**If running on a cluster (n=12, many parameter points):**
- Use `julia -t N` where N = cores per node
- Each SLURM job runs one parameter point
- JLD2 per job, merge results post-hoc with a collection script
- Consider DrWatson in the experiment management scripts (not in the package)
- Distributed.jl remains relevant for multi-node runs within a single job

**If preparing paper figures:**
- Use CSV export from DataFrames for plotting in Julia (Plots.jl/Makie.jl) or Python (matplotlib)
- JLD2 for archiving full density matrices and convergence data
- Save metadata (git commit, Julia version, thread count) alongside results for reproducibility

## Version Compatibility

| Package | Compatible With | Notes |
|---------|-----------------|-------|
| JLD2 >= 0.5 | Julia >= 1.6 | Pure Julia. v0.6.3 is latest (Nov 2024). `[compat]` should be `"0.5, 0.6"` to accept both major-minor ranges. No known issues with Julia 1.11/1.12. |
| DataFrames >= 1.7 | Julia >= 1.6 | Stable API since 1.0. Actively maintained. |
| CSV >= 0.10 | Julia >= 1.6 | Lightweight. Works with DataFrames seamlessly. |
| Base.Threads (stdlib) | Julia >= 1.3 (basic), 1.7+ (TaskLocalRNG) | Julia 1.11+ required per Project.toml `[compat]`. TaskLocalRNG fully stable. |
| Random (stdlib) | Julia >= 1.7 (task-local RNG) | Default RNG is task-local. `Random.seed!()` seeds the current task's RNG. |

## Cleanup Opportunities

While adding stack for this milestone:

| Item | Current | Should Be | Impact |
|------|---------|-----------|--------|
| SharedArrays in `[deps]` | Used by NUFFT for Distributed | Keep but do NOT use for threading | Avoid confusion: threading shares memory automatically |
| Distributed in `[deps]` | Used by `#! uncomment for multi-threads` comment in furnace.jl | Keep for future cluster use | The TODO comment in `construct_lindbladian` uses `@distributed`, which is multi-process, not multi-thread. Do not confuse with `Threads.@spawn`. |
| BSON in `[deps]` | Used for legacy Hamiltonian loading | Keep, but use JLD2 for all new data | Do not use BSON for new experiment results |

## Sources

- [Julia Multi-Threading documentation](https://docs.julialang.org/en/v1/manual/multi-threading/) -- `@spawn`, `@sync`, task-based parallelism. HIGH confidence.
- [Julia Random stdlib documentation](https://docs.julialang.org/en/v1/stdlib/Random/) -- TaskLocalRNG, per-task seeding, `Random.seed!` behavior. HIGH confidence.
- [JLD2.jl official documentation](https://juliaio.github.io/JLD2.jl/stable/) -- API, HDF5 compatibility, type preservation. HIGH confidence.
- [JLD2.jl GitHub releases](https://github.com/JuliaIO/JLD2.jl/releases) -- v0.6.3 (Nov 2024), latest stable. HIGH confidence.
- [OhMyThreads.jl documentation](https://juliafolds2.github.io/OhMyThreads.jl/stable/) -- TaskLocalValue, tmapreduce API, scheduler options. Evaluated but not recommended. HIGH confidence.
- [OhMyThreads.jl thread-safe storage](https://juliafolds2.github.io/OhMyThreads.jl/stable/literate/tls/tls/) -- TaskLocalValue pattern for per-task buffers. HIGH confidence.
- [ITensors.jl multithreading guide](https://itensor.github.io/ITensors.jl/dev/Multithreading.html) -- BLAS.set_num_threads(1) pattern for Julia threading + BLAS coexistence. HIGH confidence.
- [Julia Discourse: BLAS vs Julia threads](https://discourse.julialang.org/t/julia-threads-vs-blas-threads/8914) -- Thread oversubscription problem and solutions. HIGH confidence.
- [Julia Discourse: Reproducible multithreaded Monte Carlo](https://discourse.julialang.org/t/reproducible-multithreaded-monte-carlo-task-local-random/35269) -- Per-task RNG patterns. HIGH confidence.
- [JuliaCon 2025: Fixing Julia's task-local RNG](https://pretalx.com/juliacon-2025/talk/ZNBEAN/) -- TaskLocalRNG design, DotMix algorithm. HIGH confidence.
- [QuantumToolbox.jl Monte Carlo solver](https://qutip.org/QuantumToolbox.jl/stable/users_guide/time_evolution/mcsolve) -- EnsembleThreads() pattern for parallel trajectories. HIGH confidence.
- [Julia Discourse: JLD2 vs Arrow performance](https://discourse.julialang.org/t/why-jld2-jl-is-40x-slower-than-arrow-jl/122217) -- JLD2 slower for flat tables, but preserves Julia types. Informed JLD2+CSV dual strategy. MEDIUM confidence.
- [DataFrames.jl documentation](https://dataframes.juliadata.org/stable/) -- Tabular data API. HIGH confidence.
- [StableRNGs.jl GitHub](https://github.com/JuliaRandom/StableRNGs.jl) -- Cross-version reproducible RNG for tests. HIGH confidence.
- [Welford's online algorithm (Wikipedia)](https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm) -- Numerically stable streaming variance. HIGH confidence (textbook algorithm).

---
*Stack research for: QuantumFurnace.jl v1.1 Multi-Threaded Trajectory Engine & KMS-vs-GNS Experiments*
*Researched: 2026-02-15*
