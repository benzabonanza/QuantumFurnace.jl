# Phase 17: Adaptive Sampling - Research

**Researched:** 2026-02-16
**Domain:** Convergence-driven adaptive trajectory batching with automatic stopping and hard cap (Julia, QuantumFurnace.jl internals)
**Confidence:** HIGH

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

#### Convergence criteria
- **Stop trigger: trace distance only** -- observables (correlations, energy) are tracked but do NOT gate stopping
- **Relative change computed via windowed average** -- compare average of last K=3 batches vs previous K=3 batches
- **3 consecutive stable checks required** -- relative change must be below threshold for 3 consecutive batch checkpoints
- **Minimum 5 batches before stopping can trigger** -- prevents premature stopping from lucky early batches

#### Non-convergence handling
- **Return with `converged=false` flag** -- no warning logged, no error thrown; user checks the flag programmatically
- **Include diagnostics** -- final relative change value, number of consecutive stable batches reached, total batches run; helps user decide next steps
- **Extend existing ConvergenceData struct** -- add `converged::Bool` and stop diagnostics fields directly to ConvergenceData (no new wrapper type)
- **Return immediately when converged** -- no extra confirmation batches after convergence criteria are met

#### Batch sizing strategy
- **Fixed batch size throughout** -- every batch has the same number of trajectories
- **Default 200 trajectories per batch** -- matches Phase 16 convergence tracking structure
- **Batch size configurable** -- user can pass batch_size kwarg, default 200
- **Wraps `run_trajectories_convergence`** -- adaptive function calls the existing Phase 16 function in a loop with increasing trajectory counts, reusing existing batch infrastructure

#### Default parameters (all configurable)
- **N_max = 20,000 trajectories** (100 batches of 200) -- generous ceiling for large systems at high beta
- **Relative change threshold = 0.01** (1%) -- user can tighten (0.001) or loosen (0.05)
- **Patience = 3 consecutive stable batches** -- user can increase for more confidence
- **Minimum batches = 5** -- user can adjust burn-in period
- **Batch size = 200** -- user can adjust per-batch trajectory count

### Claude's Discretion
- How to accumulate windowed averages efficiently
- Internal loop structure for wrapping run_trajectories_convergence
- Exact field names and types for diagnostics in ConvergenceData
- How to handle the seed progression across adaptive batches

### Deferred Ideas (OUT OF SCOPE)
None -- discussion stayed within phase scope
</user_constraints>

## Summary

Phase 17 adds an adaptive stopping layer on top of the Phase 16 convergence tracking infrastructure. The core idea: instead of requiring the user to specify a fixed number of batches, the system monitors trace distance relative change across batches and stops automatically when convergence is detected (or a hard maximum is reached). This is a thin orchestration layer -- the heavy lifting (trajectory execution, density matrix accumulation, trace distance measurement) already lives in `run_trajectories_convergence`.

The implementation has two main parts: (1) extending the `ConvergenceData` struct with convergence status and diagnostic fields (`converged`, `final_relative_change`, `consecutive_stable_batches`, `total_batches`), and (2) creating a new `run_trajectories_adaptive` function that calls `run_trajectories_convergence` in a loop of single-batch increments, checking the windowed-average convergence criterion after each batch. The windowed average compares the mean trace distance of the last K=3 batches against the mean of the previous K=3 batches. Convergence is declared when this relative change stays below the threshold for 3 consecutive checks, with a minimum of 5 batches before stopping can trigger.

The key design decision is how the adaptive function wraps `run_trajectories_convergence`. Since the Phase 16 function runs a fixed number of batches in a loop and returns accumulated results, the adaptive function should NOT call `run_trajectories_convergence` repeatedly with different n_batches. Instead, the adaptive function should implement its own batch loop that reuses the internal logic of `run_trajectories_convergence` (call `run_trajectories` per batch, accumulate rho, measure metrics) while adding the convergence check. This avoids redundant computation and keeps seed management clean.

**Primary recommendation:** Create `run_trajectories_adaptive` as a new function in `convergence.jl` that implements its own batch loop (mirroring `run_trajectories_convergence` internals) with convergence checking. Extend `ConvergenceData` with 4 new fields for adaptive diagnostics. The function should be a thin layer: ~80 lines of new code, mostly the convergence check logic and the adaptive loop wrapper around the same `run_trajectories` + accumulate + measure pattern already in Phase 16.

## Standard Stack

### Core

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Existing `convergence.jl` | Phase 16 | `ConvergenceData` struct, `run_trajectories_convergence` pattern, observable builders | The adaptive layer extends this directly |
| Existing `trajectories.jl` | Phase 13 | `run_trajectories` for per-batch execution | Called in the batch loop, handles threading/seeding |
| Existing `qi_tools.jl` | N/A | `trace_distance_h`, `hermitianize!` | Checkpoint measurements (same as Phase 16) |
| `LinearAlgebra` (stdlib) | Julia 1.11+ | `Hermitian`, `tr`, matrix ops | Standard Julia scientific computing |
| `Random` (stdlib) | Julia 1.11+ | `RandomDevice`, `Xoshiro` | Seed management for adaptive batches |

### Supporting

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| Existing `results.jl` | Phase 15 | `_convergence_to_dict`, `_dict_to_convergence` | Must be updated for new ConvergenceData fields |
| Existing `constants.jl` | N/A | Pauli matrices | Observable construction (unchanged from Phase 16) |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Own batch loop in adaptive function | Calling `run_trajectories_convergence` repeatedly | Repeated calls waste computation: each call re-runs all previous batches since it starts from scratch. Own loop avoids this. |
| Extending ConvergenceData | New AdaptiveResult wrapper | CONTEXT.md locks decision: extend ConvergenceData directly (no new wrapper type) |
| Windowed average of trace distances | Raw batch-to-batch relative change | Windowed average is smoother, avoids false convergence from noisy single batches (locked decision) |

## Architecture Patterns

### Current State (After Phase 16)

```
run_trajectories_convergence(jumps, config, psi0, ham;
    gibbs, observables, observable_names,
    batch_size, n_batches, seed, ...)
  -> (TrajectoryResult, ConvergenceData)

ConvergenceData:
    batch_sizes::Vector{Int}
    cumulative_n_traj::Vector{Int}
    trace_distances::Vector{Float64}
    observable_names::Vector{String}
    observable_values::Matrix{Float64}      # n_obs x n_checkpoints
    observable_gibbs_values::Vector{Float64}
```

The Phase 16 function runs a fixed `n_batches` loop. All storage is pre-allocated to `n_batches` length. There is no convergence checking or early stopping.

### Target State (After Phase 17)

```
ConvergenceData (EXTENDED):
    batch_sizes::Vector{Int}
    cumulative_n_traj::Vector{Int}
    trace_distances::Vector{Float64}
    observable_names::Vector{String}
    observable_values::Matrix{Float64}          # n_obs x n_checkpoints
    observable_gibbs_values::Vector{Float64}
    converged::Bool                              # NEW: did adaptive stopping trigger?
    final_relative_change::Float64               # NEW: relative change at termination
    consecutive_stable_batches::Int              # NEW: how many consecutive stable checks achieved
    total_batches::Int                           # NEW: number of batches actually run

run_trajectories_adaptive(jumps, config, psi0, hamiltonian;
    gibbs, observables, observable_names,
    batch_size=200, n_max=20_000,
    convergence_threshold=0.01, patience=3,
    min_batches=5, window_size=3,
    seed=nothing, trotter=nothing,
    total_time=config.mixing_time, delta=config.delta)
  -> (TrajectoryResult, ConvergenceData)
```

The adaptive function returns the **same** `(TrajectoryResult, ConvergenceData)` tuple as the fixed-count function, but with the new diagnostic fields populated. When called with fixed-count mode (`run_trajectories_convergence`), the new fields default to `converged=false`, `final_relative_change=NaN`, `consecutive_stable_batches=0`, `total_batches=n_batches`.

### Pattern 1: Adaptive Batch Loop with Convergence Check

**What:** Run batches one at a time, checking convergence after each batch (after the minimum burn-in). Use push!-based dynamic arrays instead of pre-allocated fixed arrays, since the total batch count is unknown.

**When to use:** Always for the adaptive function. The fixed-count function retains pre-allocated arrays.

**Pseudocode:**
```julia
function run_trajectories_adaptive(jumps, config, psi0, hamiltonian; kwargs...)
    max_batches = div(n_max, batch_size)

    # Dynamic storage (unknown final length)
    trace_dists = Float64[]
    obs_values_list = Vector{Float64}[]   # will hcat at the end
    cum_n_traj = Int[]
    batch_sizes_vec = Int[]

    rho_acc = zeros(CT, dim, dim)
    n_total = 0
    consecutive_stable = 0
    last_relative_change = NaN

    for batch_idx in 1:max_batches
        # Run one batch (same pattern as run_trajectories_convergence)
        batch_seed = actual_seed + n_total
        result = run_trajectories(jumps, config, psi0, hamiltonian;
            ntraj=batch_size, seed=batch_seed, ...)

        # Accumulate and measure (identical to Phase 16)
        rho_acc .+= result.rho_mean .* batch_size
        n_total += batch_size
        rho_running = rho_acc ./ n_total
        hermitianize!(rho_running)
        td = trace_distance_h(Hermitian(rho_running), gibbs)
        push!(trace_dists, td)
        # ... push obs values ...

        # Check convergence (only after min_batches)
        if batch_idx >= min_batches
            rel_change = _compute_windowed_relative_change(trace_dists, window_size)
            last_relative_change = rel_change

            if rel_change < convergence_threshold
                consecutive_stable += 1
            else
                consecutive_stable = 0
            end

            if consecutive_stable >= patience
                converged = true
                break
            end
        end
    end

    # Build ConvergenceData with diagnostic fields
    # ...
end
```

### Pattern 2: Windowed Relative Change Computation

**What:** Compare the mean trace distance of the last `W` batches against the mean of the previous `W` batches. The relative change is `|mean_recent - mean_previous| / max(mean_previous, eps)`.

**When to use:** At each convergence check point (every batch after `min_batches`).

**Key design detail:** The window size `W` defaults to 3 (matching the `patience` parameter). The check requires at least `2*W` data points, so `min_batches` must be >= `2*W`. With defaults (`min_batches=5`, `window_size=3`), the first check happens at batch 6 (need 6 data points for two windows of 3). Actually -- `min_batches=5` with `window_size=3` means at batch 5 we have 5 data points. Two windows of 3 need indices [1,2,3] and [3,4,5] with overlap, or [1,2] and [3,4,5] without. The cleanest approach: require `min_batches >= 2 * window_size` and use non-overlapping windows: `previous = trace_dists[end-2W+1:end-W]`, `recent = trace_dists[end-W+1:end]`. If `min_batches < 2 * window_size`, silently clamp `min_batches` to `2 * window_size`.

**Implementation:**
```julia
function _compute_windowed_relative_change(
    trace_dists::Vector{Float64},
    window_size::Int
)
    n = length(trace_dists)
    if n < 2 * window_size
        return Inf  # not enough data, return large value (won't trigger convergence)
    end

    recent = @view trace_dists[n - window_size + 1 : n]
    previous = @view trace_dists[n - 2*window_size + 1 : n - window_size]

    mean_recent = sum(recent) / window_size
    mean_previous = sum(previous) / window_size

    # Relative change: |new - old| / max(|old|, eps)
    return abs(mean_recent - mean_previous) / max(abs(mean_previous), eps(Float64))
end
```

### Pattern 3: Seed Management Across Adaptive Batches

**What:** Same pattern as Phase 16: `batch_seed = actual_seed + n_total` before each batch. Since each batch adds exactly `batch_size` trajectories, the seeds form a contiguous non-overlapping sequence regardless of how many batches run.

**When to use:** Always. The seed progression is deterministic: for a given `seed`, the adaptive function produces identical results to running a fixed-count function with the same total number of batches (since the per-batch seeds are identical).

**Key insight:** This means `run_trajectories_adaptive` with seed=42 that converges after 15 batches produces *exactly* the same trace distances for batches 1-15 as `run_trajectories_convergence` with seed=42, n_batches=15. The seed management is the same; the only difference is when the loop stops.

### Pattern 4: ConvergenceData Struct Extension (Backward Compatible)

**What:** Add 4 new fields to `ConvergenceData`. Since ConvergenceData uses a positional constructor (no @kwdef), the new fields go at the end with default-value outer constructors for backward compatibility.

**Implementation:**
```julia
struct ConvergenceData
    # Existing fields (Phase 16)
    batch_sizes::Vector{Int}
    cumulative_n_traj::Vector{Int}
    trace_distances::Vector{Float64}
    observable_names::Vector{String}
    observable_values::Matrix{Float64}
    observable_gibbs_values::Vector{Float64}
    # New fields (Phase 17)
    converged::Bool
    final_relative_change::Float64
    consecutive_stable_batches::Int
    total_batches::Int
end

# Backward-compatible constructor (Phase 16 callers pass 6 args)
function ConvergenceData(
    batch_sizes, cumulative_n_traj, trace_distances,
    observable_names, observable_values, observable_gibbs_values
)
    return ConvergenceData(
        batch_sizes, cumulative_n_traj, trace_distances,
        observable_names, observable_values, observable_gibbs_values,
        false, NaN, 0, length(batch_sizes),
    )
end
```

This means the existing Phase 16 `run_trajectories_convergence` function does not need to change -- its 6-argument ConvergenceData constructor call still works, and the new fields get sensible defaults (`converged=false`, `final_relative_change=NaN`, `consecutive_stable_batches=0`, `total_batches=n_batches`).

### Pattern 5: Dynamic Storage with Final Conversion

**What:** Since the adaptive function does not know the final batch count, use `push!`-based vectors for trace distances, observable values, etc. At the end, convert to the fixed-size arrays expected by ConvergenceData.

**When to use:** In the adaptive function only. The fixed-count function retains pre-allocated arrays.

**For observable_values:** Use a `Vector{Vector{Float64}}` during accumulation, then `hcat(obs_values_list...)` or `reduce(hcat, obs_values_list)` to produce the `n_obs x n_batches` Matrix at the end.

### Anti-Patterns to Avoid

- **Calling `run_trajectories_convergence` in a loop with increasing n_batches:** Each call re-runs all previous batches from scratch. For 50 batches, this runs 1+2+3+...+50 = 1275 batch-equivalents instead of 50. The adaptive function must implement its own single-pass batch loop.
- **Using `run_trajectories_convergence` with n_batches=1 in a loop:** While this avoids redundant computation, it creates a new `rho_acc` accumulator on every call and discards it. The running average would need to be maintained externally, duplicating the logic in `run_trajectories_convergence`. Better to implement the loop directly.
- **Checking convergence on raw batch-to-batch change:** Single-batch trace distances can be noisy (especially early on). The windowed average smooths this out (locked decision).
- **Modifying `run_trajectories_convergence` to add adaptive logic:** Keep the fixed-count function simple and predictable. The adaptive function is a separate entry point with different parameters.
- **Pre-allocating storage to `max_batches` length:** For `n_max=20_000` with `batch_size=200`, that is 100 batches. Pre-allocating 100 entries is fine, but most runs converge much earlier. Use push!-based dynamic storage to avoid confusion about which entries are valid.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Per-batch trajectory execution | Custom trajectory loop | `run_trajectories(...)` with `ntraj=batch_size` | Handles threading, seeding, workspace management |
| Trace distance computation | Custom matrix norm | `trace_distance_h(Hermitian(rho), gibbs)` | Already handles Hermitianization, eigendecomposition |
| Observable expectation values | Manual matrix multiplication | `real(tr(rho_running * obs))` pattern from Phase 16 | Consistent with existing convergence tracking |
| Density matrix accumulation | Custom accumulator | Same `rho_acc .+= result.rho_mean .* batch_size` pattern | Proven in Phase 16 integration tests |
| Hermitianization | Inline `(A+A')/2` | `hermitianize!(A)` from qi_tools.jl | In-place, consistent with codebase |
| Seed management | Custom seed scheme | `batch_seed = actual_seed + n_total` (Phase 16 pattern) | Non-overlapping, deterministic, tested |

**Key insight:** Phase 17 introduces only one genuinely new computation: the windowed relative change of trace distances. Everything else (batch execution, accumulation, measurement, seed management) is identical to Phase 16. The adaptive function is ~30 lines of convergence-checking logic wrapped around ~50 lines of batch-loop code copied from `run_trajectories_convergence`.

## Common Pitfalls

### Pitfall 1: Windowed Average Requires Enough Batches

**What goes wrong:** If `min_batches` is set lower than `2 * window_size`, the windowed relative change computation does not have enough data points and either crashes (index out of bounds) or returns garbage.

**Why it happens:** Two non-overlapping windows of size `W` need at least `2W` data points. With defaults `window_size=3`, this means at least 6 batches.

**How to avoid:** Clamp `min_batches` to be at least `2 * window_size` at the start of the adaptive function. If the user passes `min_batches=3` with `window_size=3`, silently set `min_batches=6`. Alternatively, return `Inf` from the relative change function when there are not enough data points (safer, no silent clamping).

**Warning signs:** Convergence triggers at batch 5 with suspiciously good relative change values.

### Pitfall 2: ConvergenceData Struct Extension Breaks Existing Code

**What goes wrong:** Adding new fields to `ConvergenceData` changes the positional constructor. All existing code that constructs a `ConvergenceData` with 6 positional arguments breaks because the constructor now expects 10 arguments.

**Why it happens:** Julia structs have a single auto-generated inner constructor that takes ALL fields positionally.

**How to avoid:** Add a backward-compatible outer constructor that accepts 6 arguments and fills in defaults for the 4 new fields. The existing `run_trajectories_convergence` and all tests continue to work without changes. The inner constructor (10 args) is used only by the new adaptive function.

**Warning signs:** `MethodError: no method matching ConvergenceData(::Vector{Int}, ...)` errors when running existing Phase 16 code.

### Pitfall 3: Dict Serialization Must Handle New Fields

**What goes wrong:** The existing `_convergence_to_dict` and `_dict_to_convergence` in `results.jl` do not include the new fields. Saving a ConvergenceData from adaptive mode loses the diagnostic information. Loading an old ConvergenceData (without new fields) crashes on missing keys.

**Why it happens:** The serialization functions were written for Phase 16's 6-field struct.

**How to avoid:** Update `_convergence_to_dict` to include the 4 new fields. Update `_dict_to_convergence` to use `get(d, :converged, false)` etc. with defaults, ensuring backward compatibility when loading files saved before Phase 17.

**Warning signs:** Missing `:converged` key error when loading old BSON files; adaptive diagnostics silently lost when saving.

### Pitfall 4: Dynamic Array Concatenation for Observable Values

**What goes wrong:** Using `hcat` repeatedly inside the loop (`obs_values = hcat(obs_values, batch_obs)`) creates a new matrix allocation on every batch. For 100 batches, this is 100 matrix copies with growing sizes.

**Why it happens:** Julia matrices cannot be appended in-place.

**How to avoid:** Collect batch observable values in a `Vector{Vector{Float64}}` during the loop, then do a single `reduce(hcat, obs_values_list)` at the end. Or pre-allocate a large matrix and track the used column count.

**Warning signs:** Unexpected GC pressure or slowdown in the adaptive loop, especially with many batches.

### Pitfall 5: Relative Change Denominator Near Zero

**What goes wrong:** If the mean of the previous window is very close to zero (trace distance nearly converged), dividing by it amplifies noise and can produce a relative change >> 1, resetting the consecutive stable counter.

**Why it happens:** Trace distance approaches zero as rho_running -> gibbs. At very small values (< 1e-10), the relative change becomes numerically unreliable.

**How to avoid:** Use `max(abs(mean_previous), eps(Float64))` as the denominator. This prevents division by zero and keeps the relative change finite. At trace distance ~ eps, the system is effectively converged. Alternatively, add an absolute threshold floor: if `mean_recent < absolute_floor` (e.g., 1e-8), declare converged regardless of relative change.

**Warning signs:** Convergence counter resets near the end of a long run when trace distances are very small.

### Pitfall 6: Adaptive Function Must Export and Be Tested

**What goes wrong:** The new function exists in `convergence.jl` but is not exported in `QuantumFurnace.jl`, or is not covered by tests.

**Why it happens:** Forgetting to add exports and tests for new public API.

**How to avoid:** Add `run_trajectories_adaptive` to the export list in `QuantumFurnace.jl`. Write tests covering: convergence detection (converged=true path), hard cap (converged=false path), determinism, backward compatibility of ConvergenceData, serialization round-trip with new fields.

**Warning signs:** `run_trajectories_adaptive` not found as exported symbol.

## Code Examples

### Extended ConvergenceData Struct

```julia
# Source: convergence.jl (Phase 16 struct, extended for Phase 17)
struct ConvergenceData
    batch_sizes::Vector{Int}
    cumulative_n_traj::Vector{Int}
    trace_distances::Vector{Float64}
    observable_names::Vector{String}
    observable_values::Matrix{Float64}      # n_obs x n_checkpoints
    observable_gibbs_values::Vector{Float64} # <O_i>_gibbs reference values
    # Phase 17: Adaptive diagnostics
    converged::Bool
    final_relative_change::Float64
    consecutive_stable_batches::Int
    total_batches::Int
end

# Backward-compatible constructor (6 args -> 10 fields with defaults)
function ConvergenceData(
    batch_sizes::Vector{Int},
    cumulative_n_traj::Vector{Int},
    trace_distances::Vector{Float64},
    observable_names::Vector{String},
    observable_values::Matrix{Float64},
    observable_gibbs_values::Vector{Float64},
)
    ConvergenceData(
        batch_sizes, cumulative_n_traj, trace_distances,
        observable_names, observable_values, observable_gibbs_values,
        false, NaN, 0, length(batch_sizes),
    )
end
```

### Windowed Relative Change Helper

```julia
# Source: New for Phase 17
"""
    _windowed_relative_change(trace_dists, window_size) -> Float64

Compute the relative change between two consecutive non-overlapping windows
of trace distance values. Returns Inf if not enough data.
"""
function _windowed_relative_change(trace_dists::Vector{Float64}, window_size::Int)
    n = length(trace_dists)
    n < 2 * window_size && return Inf

    mean_recent = sum(@view trace_dists[n - window_size + 1 : n]) / window_size
    mean_previous = sum(@view trace_dists[n - 2*window_size + 1 : n - window_size]) / window_size

    return abs(mean_recent - mean_previous) / max(abs(mean_previous), eps(Float64))
end
```

### Adaptive Function Skeleton

```julia
# Source: New for Phase 17 (wrapping Phase 16 batch-loop pattern)
function run_trajectories_adaptive(
    jumps::Vector{JumpOp},
    config::AbstractThermalizeConfig,
    psi0::Vector{<:Complex},
    hamiltonian::HamHam;
    gibbs::Hermitian,
    observables::Vector{<:Matrix{<:Complex}},
    observable_names::Vector{String},
    batch_size::Int = 200,
    n_max::Int = 20_000,
    convergence_threshold::Float64 = 0.01,
    patience::Int = 3,
    min_batches::Int = 5,
    window_size::Int = 3,
    seed::Union{Int,Nothing} = nothing,
    trotter::Union{TrottTrott,Nothing} = nothing,
    total_time::Real = config.mixing_time,
    delta::Real = config.delta,
)
    actual_seed = seed === nothing ? Int(rand(Random.RandomDevice(), UInt64) >> 1) : seed
    max_batches = cld(n_max, batch_size)  # ceiling division

    # Ensure min_batches allows at least 2 full windows
    effective_min = max(min_batches, 2 * window_size)

    CT = eltype(psi0)
    dim = length(psi0)
    n_obs = length(observables)

    rho_acc = zeros(CT, dim, dim)
    n_total = 0
    obs_gibbs = _compute_gibbs_observable_values(gibbs, observables)

    # Dynamic storage
    trace_dists = Float64[]
    obs_values_list = Vector{Float64}[]
    cum_n_traj = Int[]
    batch_sizes_vec = Int[]

    converged = false
    consecutive_stable = 0
    last_relative_change = NaN

    for batch_idx in 1:max_batches
        batch_seed = actual_seed + n_total

        result = run_trajectories(
            jumps, config, psi0, hamiltonian;
            trotter=trotter, total_time=total_time,
            delta=delta, ntraj=batch_size, seed=batch_seed,
        )

        rho_acc .+= result.rho_mean .* batch_size
        n_total += batch_size
        rho_running = rho_acc ./ n_total
        hermitianize!(rho_running)

        td = trace_distance_h(Hermitian(rho_running), gibbs)
        push!(trace_dists, td)

        batch_obs = [real(tr(rho_running * observables[i])) for i in 1:n_obs]
        push!(obs_values_list, batch_obs)
        push!(cum_n_traj, n_total)
        push!(batch_sizes_vec, batch_size)

        # Convergence check (only after burn-in)
        if batch_idx >= effective_min
            rel_change = _windowed_relative_change(trace_dists, window_size)
            last_relative_change = rel_change

            if rel_change < convergence_threshold
                consecutive_stable += 1
            else
                consecutive_stable = 0
            end

            if consecutive_stable >= patience
                converged = true
                break
            end
        end
    end

    # Final density matrix
    rho_final = rho_acc ./ n_total
    hermitianize!(rho_final)

    # Build observable values matrix from collected vectors
    obs_values = reduce(hcat, obs_values_list)

    conv_data = ConvergenceData(
        batch_sizes_vec, cum_n_traj, trace_dists,
        observable_names, obs_values, obs_gibbs,
        converged, last_relative_change, consecutive_stable, length(trace_dists),
    )

    traj_result = TrajectoryResult(rho_final, n_total, actual_seed, nothing, nothing)

    return traj_result, conv_data
end
```

### Updated Dict Serialization (Backward Compatible)

```julia
# Source: results.jl (updated for Phase 17 fields)
function _convergence_to_dict(conv::ConvergenceData)
    return Dict{Symbol, Any}(
        :batch_sizes              => conv.batch_sizes,
        :cumulative_n_traj        => conv.cumulative_n_traj,
        :trace_distances          => conv.trace_distances,
        :observable_names         => conv.observable_names,
        :observable_values        => conv.observable_values,
        :observable_gibbs_values  => conv.observable_gibbs_values,
        :converged                => conv.converged,
        :final_relative_change    => conv.final_relative_change,
        :consecutive_stable_batches => conv.consecutive_stable_batches,
        :total_batches            => conv.total_batches,
    )
end

function _dict_to_convergence(d::Dict)
    return ConvergenceData(
        d[:batch_sizes],
        d[:cumulative_n_traj],
        d[:trace_distances],
        d[:observable_names],
        d[:observable_values],
        get(d, :observable_gibbs_values, Float64[]),
        get(d, :converged, false),
        get(d, :final_relative_change, NaN),
        get(d, :consecutive_stable_batches, 0),
        get(d, :total_batches, length(d[:batch_sizes])),
    )
end
```

## State of the Art

| Phase 16 (Fixed-count) | Phase 17 (Adaptive) | Difference |
|---|---|---|
| User specifies `n_batches` | User specifies `n_max` (hard cap) | Automatic stopping eliminates guesswork |
| Always runs all n_batches | Stops early when converged | Saves compute time for well-behaved systems |
| No convergence judgment | `converged::Bool` flag in result | User gets programmatic convergence status |
| No diagnostics | `final_relative_change`, `consecutive_stable_batches` | User can debug non-convergence |
| Pre-allocated arrays | Push-based dynamic arrays | Handles variable-length runs |
| ConvergenceData has 6 fields | ConvergenceData has 10 fields | Backward-compatible via outer constructor |

**Scope boundary:**
- Phase 16 = fixed-count convergence monitoring (CONV-01, CONV-02, CONV-03) -- COMPLETE
- Phase 17 = adaptive stopping based on convergence criteria (CONV-04, CONV-05)
- Phase 18 = KMS-vs-GNS experiments using adaptive sampling

## Open Questions

1. **Should `run_trajectories_adaptive` validate that `n_max` is a multiple of `batch_size`?**
   - What we know: `max_batches = cld(n_max, batch_size)` rounds up. With `n_max=20_000` and `batch_size=200`, that is exactly 100 batches. With `n_max=1_000` and `batch_size=300`, that is 4 batches (1200 trajectories, exceeding n_max by 200).
   - What's unclear: Whether exceeding `n_max` slightly is acceptable, or whether the last batch should be truncated.
   - Recommendation: Use ceiling division and accept the slight overshoot. Truncating the last batch would break the fixed-batch-size invariant (locked decision: "fixed batch size throughout"). Document this in the docstring.

2. **Should the Phase 16 `run_trajectories_convergence` function be refactored to share code with the adaptive function?**
   - What we know: Both functions have nearly identical batch-loop internals (run batch, accumulate, measure). The key difference is pre-allocated vs. dynamic storage and the convergence check.
   - What's unclear: Whether code duplication or abstraction is preferable.
   - Recommendation: Accept moderate code duplication (~40 lines of shared loop body). The two functions serve different purposes and the adaptive function has enough differences (dynamic storage, convergence checking, different parameters) that a shared helper would be over-engineered. The Phase 16 function remains untouched except for the ConvergenceData struct extension (backward compatible).

## Sources

### Primary (HIGH confidence)
- **QuantumFurnace.jl codebase** -- Direct analysis of:
  - `src/convergence.jl` (236 lines): ConvergenceData struct, run_trajectories_convergence batch loop, observable builders
  - `src/trajectories.jl` (914 lines): run_trajectories, seed management, threading
  - `src/results.jl` (434 lines): _convergence_to_dict, _dict_to_convergence serialization
  - `src/structs.jl` (317 lines): TrajectoryResult, config structs
  - `src/qi_tools.jl`: trace_distance_h, hermitianize!
  - `test/test_convergence.jl` (350 lines): 10 testsets, 106 assertions covering Phase 16
  - `test/test_helpers.jl` (380 lines): shared test fixtures, config factories
- **Phase 16 planning docs** -- 16-01-PLAN.md (architecture decisions), 16-RESEARCH.md (convergence patterns, pitfalls), 16-VERIFICATION.md (10/10 truths verified)
- **Phase 17 CONTEXT.md** -- Locked decisions on convergence criteria, non-convergence handling, batch sizing, default parameters

### Secondary (MEDIUM confidence)
- **ROADMAP.md** -- Phase 17 success criteria (3 items), dependency on Phase 15 + 16
- **REQUIREMENTS.md** -- CONV-04 (adaptive sampling), CONV-05 (hard cap)
- **STATE.md** -- Prior decisions affecting Phase 17 (adaptive = batch convergence, relative change <1% for 3 consecutive batches, hard cap N_max)

### Tertiary (LOW confidence)
- None -- all findings verified against codebase analysis

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- No new dependencies; adaptive layer is pure Julia orchestration over existing Phase 16 components
- Architecture: HIGH -- Single new function + 4 struct fields + serialization update. Pattern directly mirrors Phase 16 batch loop with convergence check added.
- Pitfalls: HIGH -- All identified from codebase analysis (backward-compat constructor, serialization update, windowed average edge cases, dynamic arrays)
- Code examples: HIGH -- All derived from existing codebase patterns (convergence.jl batch loop, results.jl serialization, qi_tools.jl trace distance)

**Research date:** 2026-02-16
**Valid until:** 60 days (stable domain; no external dependency drift; all code is internal)
