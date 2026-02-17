# Phase 22: Observable-Only Trajectory Runner - Research

**Researched:** 2026-02-17
**Domain:** Quantum trajectory simulation with observable-only measurement (no per-trajectory DM reconstruction)
**Confidence:** HIGH

## Summary

Phase 22 adds a new trajectory runner `run_observable_trajectories` that measures time-resolved observables `<O>(t)` from quantum trajectories without the per-trajectory overhead of density matrix reconstruction. The existing `run_trajectories` already has a "with observables" path (lines 571-639 of `trajectories.jl`) that accumulates BOTH observable measurements AND the `|psi><psi|` outer product at each trajectory's end. The new runner eliminates the `rho_acc` accumulation by default, reducing the per-trajectory memory cost from `O(dim^2)` to `O(num_obs * num_saves)` and removing the rank-1 `mul!(rho_acc, psi, psi', ...)` call per trajectory.

The implementation is straightforward because all building blocks exist. The step loop (`step_along_trajectory!`), threading pattern (`@sync`/`Threads.@spawn` with per-task workspaces), observable measurement (`_accumulate_measurements!`), and framework setup (`_build_framework_and_seed`) are all production-ready from Phases 12-13. The new function is essentially a refactored version of the existing observable path that: (1) drops `rho_acc` from the hot loop, (2) adds an optional `reconstruct_dm=true` flag to re-enable it when needed, and (3) returns a result type that may or may not contain a density matrix.

**Primary recommendation:** Create `run_observable_trajectories` as a thin wrapper that reuses `_build_framework_and_seed`, `step_along_trajectory!`, and `_accumulate_measurements!`, with a new `_run_chunk_obs_only!` inner function that omits density matrix accumulation. Use a new `ObservableTrajectoryResult` struct (or extend `TrajectoryResult` with a `nothing` rho_mean) to hold the observable-only output.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| LinearAlgebra (stdlib) | Julia 1.11+ | `dot`, `mul!`, `rmul!`, matrix ops | Already used in all trajectory code |
| Random (stdlib) | Julia 1.11+ | `Xoshiro` per-trajectory RNG | Existing threading/determinism pattern |
| Base.Threads (stdlib) | Julia 1.11+ | `@spawn`, `@sync`, `nthreads()` | Existing multi-threaded pattern |

### Supporting
No new dependencies needed. All building blocks exist in the codebase.

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| New `ObservableTrajectoryResult` struct | Reuse `TrajectoryResult` with `rho_mean=nothing` | `TrajectoryResult.rho_mean` is typed `Matrix{T}`, not `Union{Nothing, Matrix{T}}`, so it cannot hold `nothing` without breaking the type. A new struct is cleaner than changing an established type. |
| New `run_observable_trajectories` function | Add `observable_only=true` kwarg to `run_trajectories` | Separate function is clearer; the existing `run_trajectories` already has two code paths (no-obs vs with-obs), adding a third would increase complexity. A dedicated function keeps each function's contract simple. |

**Installation:** No new packages needed.

## Architecture Patterns

### Recommended Project Structure
```
src/
  trajectories.jl    # ADD new function here (same file as existing trajectory runners)
  structs.jl         # ADD new ObservableTrajectoryResult struct here
  QuantumFurnace.jl  # ADD exports for new function and struct
test/
  test_observable_trajectories.jl  # NEW test file for Phase 22
  runtests.jl        # ADD include for new test file
```

### Pattern 1: Observable-Only Workspace (No rho_acc)

**What:** Create a lightweight workspace struct (or simply reuse `TrajectoryWorkspace` with the understanding that `rho_acc` goes unused) for the observable-only path.
**When to use:** When `reconstruct_dm=false` (the default).
**Recommendation:** Reuse `TrajectoryWorkspace` as-is. The `rho_acc` field is `dim x dim` (e.g., 16x16 = 256 complex numbers for 4 qubits), which is negligible compared to framework memory. The hot-loop savings come from not calling `_accumulate_density_matrix!`, not from skipping the allocation.

**Rationale:** Creating a separate workspace type would require duplicating the workspace-creation logic and threading dispatch. The `rho_acc` allocation is one-time per thread (not per trajectory) and is small. The performance win comes from skipping the `mul!(rho_acc, psi, psi', one(CT), one(CT))` BLAS call at the end of every trajectory.

### Pattern 2: _run_chunk_obs_only! Inner Function

**What:** A new chunk runner that measures observables but does NOT call `_accumulate_density_matrix!`.
**When to use:** Default behavior of `run_observable_trajectories`.

**Example:**
```julia
# Source: modeled after _run_chunk_with_obs! (trajectories.jl lines 391-429)
function _run_chunk_obs_only!(
    ws::TrajectoryWorkspace{<:Complex},
    fw::TrajectoryFramework{<:Complex},
    psi0::Vector{<:Complex},
    chunk::UnitRange{Int},
    master_seed::Int,
    total_time::Real,
    observables::Vector{<:Matrix{<:Complex}},
    save_every::Int,
    num_steps::Int,
    num_saves::Int,
    mean_data_local::Matrix{Float64},
)
    psi = copy(psi0)
    tmp_meas = ws.psi_tmp  # reuse workspace vector as gemv buffer

    for traj_id in chunk
        rng = Random.Xoshiro(master_seed + traj_id)
        copyto!(psi, psi0)

        n2 = real(dot(psi, psi))
        rmul!(psi, 1.0 / sqrt(max(n2, eps(Float64))))

        # save at t=0
        _accumulate_measurements!(mean_data_local, 1, psi, observables, tmp_meas)

        save_idx = 1
        for step in 1:num_steps
            step_along_trajectory!(psi, fw, ws, rng)
            if step % save_every == 0
                save_idx += 1
                _accumulate_measurements!(mean_data_local, save_idx, psi, observables, tmp_meas)
            end
        end
        # NO _accumulate_density_matrix! call here
    end
    return nothing
end
```

### Pattern 3: Optional DM Reconstruction at End

**What:** When `reconstruct_dm=true`, accumulate `|psi><psi|` at the END of each trajectory (same as existing code), but this is opt-in, not default.
**When to use:** When the user needs the averaged density matrix for downstream analysis (e.g., trace distance to Gibbs).

**Implementation:** Use the existing `_run_chunk_with_obs!` when `reconstruct_dm=true`, and `_run_chunk_obs_only!` when `reconstruct_dm=false`. This avoids code duplication while giving the user control.

### Pattern 4: ObservableTrajectoryResult Struct

**What:** A dedicated result type for observable-only runs.
**Why needed:** The existing `TrajectoryResult` has `rho_mean::Matrix{T}` (not `Union{Nothing, Matrix{T}}`), so it cannot represent "no DM was computed." Changing the type would break all existing callers.

```julia
struct ObservableTrajectoryResult{T}
    times::Vector{Float64}
    measurements_mean::Matrix{Float64}  # n_obs x n_saves
    n_trajectories::Int
    seed::Int
    rho_mean::Union{Nothing, Matrix{T}}  # nothing when reconstruct_dm=false
end
```

### Pattern 5: Reuse Existing Framework and Threading

**What:** The function reuses `_build_framework_and_seed` for one-time setup, `_partition_trajectories` for thread load balancing, and the `BLAS.set_num_threads(1)` / `@sync`+`@spawn` pattern for multi-threaded execution.

**Example threading structure:**
```julia
if ntraj > 1 && Threads.nthreads() > 1
    nt = min(Threads.nthreads(), ntraj)
    chunks = _partition_trajectories(1:ntraj, nt)
    ws_per_task = [TrajectoryWorkspace(CT, dim) for _ in 1:length(chunks)]
    mean_data_per_task = [zeros(Float64, num_obs, num_saves) for _ in 1:length(chunks)]

    old_blas = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        @sync for (idx, chunk) in enumerate(chunks)
            Threads.@spawn _run_chunk_obs_only!(
                ws_per_task[idx], fw, psi0, chunk, actual_seed, total_time,
                observables, save_every, num_steps, num_saves, mean_data_per_task[idx])
        end
    finally
        BLAS.set_num_threads(old_blas)
    end

    mean_data = sum(mean_data_per_task)
    mean_data ./= ntraj
else
    # Serial path
    ...
end
```

### Anti-Patterns to Avoid
- **Modifying `TrajectoryResult` type signature:** Changing `rho_mean::Matrix{T}` to `Union{Nothing, Matrix{T}}` would break all existing callers, serialization, and tests. Use a new struct instead.
- **Duplicating the step loop:** Do NOT copy the `step_along_trajectory!` call or the framework setup. Reuse the existing functions.
- **Allocating per-trajectory measurement buffers:** The existing `_accumulate_measurements!` is allocation-free (writes into pre-allocated `acc[:, save_idx]`). Do not create per-trajectory measurement arrays.
- **Creating a new workspace type:** The existing `TrajectoryWorkspace` works perfectly. The `rho_acc` field costs negligible memory and simply goes unused.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Framework setup | Custom init logic | `_build_framework_and_seed(jumps, config, psi0, ham; ...)` | One-time setup already handles all domain types, Kraus precomputation, seed generation |
| Step evolution | Custom step loop | `step_along_trajectory!(psi, fw, ws, rng)` | Already optimized, domain-dispatched, tested across Phases 1-5 |
| Observable measurement | Custom expectation value | `_accumulate_measurements!(acc, idx, psi, obs, tmp)` | Allocation-free, already tested in existing observable path |
| Thread partitioning | Custom load balancing | `_partition_trajectories(range, n_chunks)` | Handles remainder distribution, already tested in Phase 13 |
| Per-thread workspace | Custom buffer management | `TrajectoryWorkspace(CT, dim)` per task | Thread-safe by construction (no shared mutable state) |

**Key insight:** This phase is 95% composition of existing building blocks. The only truly new code is: (1) `_run_chunk_obs_only!` (a simplified copy of `_run_chunk_with_obs!` minus the DM accumulation), (2) `ObservableTrajectoryResult` struct, and (3) `run_observable_trajectories` function that wires everything together.

## Common Pitfalls

### Pitfall 1: Seed Sequence Must Match Existing Runner for Cross-Validation
**What goes wrong:** Observable time series from `run_observable_trajectories(seed=42)` does not match `run_trajectories(seed=42, observables=obs)` despite identical parameters.
**Why it happens:** The seed-to-trajectory mapping `Xoshiro(master_seed + traj_id)` must be identical. If the traj_id range differs (e.g., starting from 0 instead of 1, or using batch offsets), seeds diverge.
**How to avoid:** Use the exact same `for traj_id in chunk` pattern with `Xoshiro(master_seed + traj_id)` where `traj_id` ranges over `1:ntraj`. This is identical to the existing `_run_chunk_with_obs!`.
**Warning signs:** Cross-validation test fails: observable means from new runner do not match existing runner for same seed, ntraj, observables.

### Pitfall 2: Division by ntraj Must Happen AFTER Thread Summation
**What goes wrong:** Each thread divides its local `mean_data_local` by its chunk size instead of the total `ntraj`.
**Why it happens:** Misunderstanding the accumulation pattern. Each thread accumulates raw sums, and the global average is computed after merging all threads.
**How to avoid:** Follow the existing pattern: each chunk accumulates raw sums in `mean_data_local`, then `sum(mean_data_per_task) ./ ntraj` computes the global average.
**Warning signs:** Observable means are off by a factor related to thread count.

### Pitfall 3: psi_tmp Reuse Conflict Between Measurement and Step
**What goes wrong:** `_accumulate_measurements!` uses `ws.psi_tmp` as a temporary buffer, but `step_along_trajectory!` also uses `ws.psi_tmp` for the K0*psi computation.
**Why it happens:** Both functions share the workspace.
**How to avoid:** This is NOT actually a problem because `_accumulate_measurements!` is called BETWEEN steps (not concurrently), and `step_along_trajectory!` leaves `psi_tmp` in an undefined state anyway. The existing `_run_chunk_with_obs!` already uses this pattern safely (line 405: `tmp_meas = ws.psi_tmp`). Follow the same pattern.
**Warning signs:** None -- this is a false alarm. Noting it here because it looks like a bug but is not.

### Pitfall 4: num_saves Calculation Must Match Existing Code
**What goes wrong:** Off-by-one in the number of save points, causing observable time series length mismatch in cross-validation.
**Why it happens:** The existing code uses `num_saves = div(num_steps, save_every) + 1` (the `+1` accounts for the t=0 measurement before any steps). If the new code uses a different formula, lengths differ.
**How to avoid:** Copy the exact calculation from `run_trajectories` (lines 576-577):
```julia
num_steps = ceil(Int, total_time / delta_step)
num_saves = div(num_steps, save_every) + 1
```
And measure at t=0 (before the step loop starts), then at every `save_every`-th step.
**Warning signs:** `length(result.times)` differs between old and new runner.

### Pitfall 5: BLAS Thread Restoration in Error Path
**What goes wrong:** If an exception occurs during threaded execution, BLAS threads remain set to 1 for the rest of the session.
**Why it happens:** Missing `try/finally` around the BLAS thread manipulation.
**How to avoid:** Always use the existing pattern:
```julia
old_blas = BLAS.get_num_threads()
BLAS.set_num_threads(1)
try
    @sync for (idx, chunk) in enumerate(chunks)
        Threads.@spawn ...
    end
finally
    BLAS.set_num_threads(old_blas)
end
```
This is already standard in the codebase (trajectories.jl lines 457-466, 593-602).

### Pitfall 6: Export and Module Registration
**What goes wrong:** New function/struct not accessible from user code.
**Why it happens:** Forgetting to add `export run_observable_trajectories, ObservableTrajectoryResult` to `QuantumFurnace.jl`.
**How to avoid:** Add exports to the trajectory block (line 42-43 of QuantumFurnace.jl).

### Pitfall 7: Forgetting to normalize psi0 at start of each trajectory
**What goes wrong:** Unnormalized initial state causes wrong expectation values.
**Why it happens:** Skipping the `rmul!(psi, 1.0 / sqrt(max(n2, eps(Float64))))` step.
**How to avoid:** Copy the normalization pattern from existing chunk runners.

## Code Examples

Verified patterns from the existing codebase:

### Existing Observable Measurement (Allocation-Free)
```julia
# Source: trajectories.jl lines 317-329
function _accumulate_measurements!(
    acc::AbstractMatrix{<:Real},
    save_idx::Int,
    psi::Vector{<:Complex},
    observables::Vector{<:Matrix{<:Complex}},
    tmp::Vector{<:Complex},
)
    @inbounds for i in eachindex(observables)
        mul!(tmp, observables[i], psi)                 # tmp := O_i * psi
        acc[i, save_idx] += real(dot(psi, tmp))        # += <psi|O_i|psi>
    end
    return nothing
end
```

### Existing Time Grid Construction
```julia
# Source: trajectories.jl lines 580-583
times = Vector{Float64}(undef, num_saves)
@inbounds for s in 1:num_saves
    times[s] = (s - 1) * save_every * delta_step
end
```

### Existing Multi-Threaded Observable Pattern
```julia
# Source: trajectories.jl lines 585-607
nt = min(Threads.nthreads(), ntraj)
chunks = _partition_trajectories(1:ntraj, nt)
ws_per_task = [TrajectoryWorkspace(CT, dim) for _ in 1:length(chunks)]
mean_data_per_task = [zeros(Float64, num_obs, num_saves) for _ in 1:length(chunks)]

old_blas = BLAS.get_num_threads()
BLAS.set_num_threads(1)
try
    @sync for (idx, chunk) in enumerate(chunks)
        Threads.@spawn _run_chunk_with_obs!(
            ws_per_task[idx], fw, psi0, chunk, actual_seed, total_time,
            observables, save_every, num_steps, num_saves, mean_data_per_task[idx])
    end
finally
    BLAS.set_num_threads(old_blas)
end
mean_data = sum(mean_data_per_task)
mean_data ./= ntraj
```

### Existing Chunk Runner with Observables
```julia
# Source: trajectories.jl lines 391-429
function _run_chunk_with_obs!(ws, fw, psi0, chunk, master_seed, total_time,
                               observables, save_every, num_steps, num_saves,
                               mean_data_local)
    psi = copy(psi0)
    tmp_meas = ws.psi_tmp

    for traj_id in chunk
        rng = Random.Xoshiro(master_seed + traj_id)
        copyto!(psi, psi0)
        n2 = real(dot(psi, psi))
        rmul!(psi, 1.0 / sqrt(max(n2, eps(Float64))))

        _accumulate_measurements!(mean_data_local, 1, psi, observables, tmp_meas)

        save_idx = 1
        for step in 1:num_steps
            step_along_trajectory!(psi, fw, ws, rng)
            if step % save_every == 0
                save_idx += 1
                _accumulate_measurements!(mean_data_local, save_idx, psi, observables, tmp_meas)
            end
        end
        _accumulate_density_matrix!(ws.rho_acc, psi)  # <-- THIS IS WHAT WE SKIP
    end
end
```

### API Signature for New Function
```julia
# Modeled after run_trajectories (trajectories.jl lines 538-550)
function run_observable_trajectories(
    jumps::Vector{JumpOp},
    config::AbstractThermalizeConfig,
    psi0::Vector{<:Complex},
    hamiltonian::HamHam;
    trotter::Union{TrottTrott,Nothing}=nothing,
    total_time::Real = config.mixing_time,
    delta::Real = config.delta,
    ntraj::Int = 1,
    observables::Vector{<:Matrix{<:Complex}},
    save_every::Int = 1,
    seed::Union{Int,Nothing} = nothing,
    reconstruct_dm::Bool = false,
)
```

## Cross-Validation Test Strategy

The critical correctness test is: for the same seed, ntraj, observables, and save_every, the observable time series from `run_observable_trajectories` MUST be bitwise identical to those from `run_trajectories` with observables.

```julia
# Cross-validation pattern
result_old = run_trajectories(jumps, config, psi0, ham;
    observables=obs, save_every=5, ntraj=20, seed=42)
result_new = run_observable_trajectories(jumps, config, psi0, ham;
    observables=obs, save_every=5, ntraj=20, seed=42)

@test result_new.measurements_mean == result_old.measurements_mean  # bitwise
@test result_new.times == result_old.times                          # bitwise
```

With `reconstruct_dm=true`:
```julia
result_dm = run_observable_trajectories(jumps, config, psi0, ham;
    observables=obs, save_every=5, ntraj=20, seed=42, reconstruct_dm=true)
@test result_dm.rho_mean == result_old.rho_mean  # bitwise
```

## Performance Characteristics

### What is saved by dropping DM accumulation

The `_accumulate_density_matrix!` call performs `mul!(rho_acc, psi, psi', one(CT), one(CT))` -- a BLAS rank-1 update (`ZHER` or `ZGERU`) on a `dim x dim` complex matrix. Cost per trajectory:
- 3-qubit (dim=8): 128 complex multiply-adds -- negligible
- 4-qubit (dim=16): 512 complex multiply-adds -- small
- 6-qubit (dim=64): 8192 complex multiply-adds -- starts to matter
- 8-qubit (dim=256): 131072 complex multiply-adds -- significant

For the target use case (gap estimation with many trajectories on larger systems), removing this call saves meaningful compute time. But for the 3-4 qubit test systems, the savings are negligible. Tests should focus on correctness, not performance.

### Memory savings

Per-thread `rho_acc` is `dim x dim` complex. With `reconstruct_dm=false`:
- The `TrajectoryWorkspace` still allocates it (trivial cost), but it is never written to.
- The main memory savings come from not needing to sum `rho_acc` across threads at the end.

## State of the Art

| Old Approach (run_trajectories with obs) | New Approach (run_observable_trajectories) | Impact |
|-------------------------------------------|-------------------------------------------|--------|
| Always computes `rho_mean` even when only observables needed | Skips DM accumulation by default | Faster for large dim, same correctness |
| Returns `TrajectoryResult` with mandatory `rho_mean` | Returns `ObservableTrajectoryResult` with optional `rho_mean` | Clearer intent in API |
| Single function handles 3 modes (no-obs, obs, obs-only) | Dedicated function for obs-only | Simpler per-function contract |

## Open Questions

1. **Should `run_observable_trajectories` reuse `TrajectoryResult` or define a new struct?**
   - What we know: `TrajectoryResult.rho_mean` is `Matrix{T}` (not Union{Nothing, ...}). Changing it would break existing code and BSON serialization.
   - Recommendation: Define a new `ObservableTrajectoryResult` struct. It is lightweight (5 fields) and avoids all backward-compatibility concerns.

2. **Should the result include a `convergence` field?**
   - What we know: The existing `TrajectoryResult` has `convergence::Union{Nothing, ConvergenceData}` for the batched convergence runner. Phase 22 does NOT implement batched convergence -- it is a single-shot runner.
   - Recommendation: No convergence field in `ObservableTrajectoryResult`. If batched observable convergence is needed later, it would be a separate function (`run_observable_trajectories_convergence` or similar).

3. **Should the time grid be returned as Float64 or match the config's precision type?**
   - What we know: The existing `run_trajectories` returns `times::Vector{Float64}` regardless of config type. The `measurements_mean` is always `Matrix{Float64}` because `real(dot(psi, tmp))` produces `Float64`.
   - Recommendation: Keep `Float64` for both `times` and `measurements_mean`, matching the existing convention.

## Sources

### Primary (HIGH confidence)
- `src/trajectories.jl` -- Complete existing implementation including `run_trajectories`, `_run_chunk_with_obs!`, `_run_chunk_no_obs!`, `_accumulate_measurements!`, `_accumulate_density_matrix!`, `_partition_trajectories`, `TrajectoryWorkspace`, `TrajectoryResult`, `TrajectoryFramework`, `_build_framework_and_seed`, `step_along_trajectory!`
- `src/structs.jl` -- `TrajectoryWorkspace`, `TrajectoryResult`, `TrajectoryFramework`, `PerOperatorKraus` struct definitions
- `src/convergence.jl` -- `build_gap_estimation_observables`, `build_total_magnetization` (Phase 20 output, used as observable source)
- `src/QuantumFurnace.jl` -- Module exports (line 42-48)
- `test/test_threading.jl` -- Deterministic multi-threaded trajectory tests, serial-threaded agreement tests
- `test/test_helpers.jl` -- Shared test fixtures (SMALL_HAM, SMALL_JUMPS, make_small_thermalize_config, etc.)

### Secondary (MEDIUM confidence)
- None needed -- all information from codebase inspection

### Tertiary (LOW confidence)
- None

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- no new dependencies, all building blocks exist in codebase
- Architecture: HIGH -- direct composition of existing, tested building blocks
- Pitfalls: HIGH -- all pitfalls derived from existing code patterns and known threading gotchas
- Code examples: HIGH -- all examples copied from verified, working source code

**Research date:** 2026-02-17
**Valid until:** No expiry -- codebase-internal research, no external dependency version concerns
