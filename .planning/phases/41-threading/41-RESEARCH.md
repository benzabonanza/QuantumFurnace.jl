# Phase 41: Threading - Research

**Researched:** 2026-03-01
**Domain:** Multi-threaded BLAS and omega-loop parallelism for `run_thermalize` DM evolution
**Confidence:** HIGH

## Summary

Phase 41 adds two complementary threading strategies to `run_thermalize`: (1) enabling multi-threaded BLAS for the dense matrix multiplications in the hot loop (the "free lunch" -- no code changes to the math, just thread count control), and (2) optional omega-loop parallelism for the `_accumulate_rho_jump!` frequency summation with per-task accumulators (the "structured parallelism" -- requires careful accumulator isolation).

The codebase already has a mature threading pattern in the trajectory engine (`trajectories.jl`): `BLAS.get_num_threads()` save, `BLAS.set_num_threads(1)` for trajectory-level parallelism, `@sync/@spawn` with per-task workspace copies, and `try/finally` restoration. The DM thermalization path (`run_thermalize`) currently runs single-threaded with whatever BLAS thread count the user has. Phase 41 inverts the trajectory pattern: instead of setting BLAS threads to 1 for Julia-level parallelism, it explicitly enables multi-threaded BLAS for the DM path since there is no Julia-level parallelism competing. For the optional omega-loop threading (THREAD-01/02), the pattern flips: BLAS is set to 1, frequency iterations are partitioned across tasks with per-task scratch accumulators, and results are summed after synchronization.

The BohrDomain (THREAD-04) has a different inner loop structure (sparse Bohr bucket iteration rather than dense frequency grid) that makes omega-loop threading more nuanced -- the bucket sizes are highly variable and the inner accumulation uses sparse index patterns rather than BLAS `mul!`. Threading BohrDomain means parallelizing the Bohr bucket iteration in `_accumulate_rho_jump!` with per-task `rho_jump` accumulators, and potentially the `_precompute_R` Bohr loop as well.

**Primary recommendation:** Implement BLAS thread enablement (THREAD-03/05) first as it provides immediate speedup with minimal code changes, then add omega-loop parallelism (THREAD-01/02/04) as a second plan with the per-task accumulator pattern from the trajectory engine.

## Standard Stack

### Core

No new libraries required. This phase uses only Julia stdlib threading facilities already imported by the project.

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| `LinearAlgebra` (stdlib) | Julia 1.11+ | `BLAS.get_num_threads()`, `BLAS.set_num_threads()`, `mul!` | Already used throughout; BLAS thread control is the standard Julia API |
| `Base.Threads` (stdlib) | Julia 1.11+ | `Threads.nthreads()`, `Threads.@spawn`, `@sync` | Already imported in `QuantumFurnace.jl` line 16 |

### Supporting

No new supporting libraries needed. The existing `ThermalizeScratch` struct provides all necessary scratch buffers.

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| `@sync/@spawn` per-task pattern | `@threads :static` with threadid-indexed buffers | PSA from Julia team explicitly discourages threadid()-based patterns; task migration makes it unsafe. The codebase already uses `@sync/@spawn` in trajectories.jl |
| Manual chunk partitioning | OhMyThreads.jl `tmap`/`treduce` | Adds external dependency; the existing `_partition_trajectories` utility is simple and battle-tested |
| Per-task scratch allocation | Task-local storage (`Base.task_local_storage`) | Per-task scratch via explicit chunk allocation is simpler, already proven in trajectory engine, and avoids TLS overhead |

## Architecture Patterns

### Current run_thermalize Hot Loop (Post Phase 39/40)

```
run_thermalize (furnace.jl)
  |-- _precompute_per_jump_channels(...)      # K0s, U_residuals precomputed
  |-- for step in 1:num_steps                 # HOT LOOP (single-threaded)
  |     |-- _apply_coherent_unitary!(...)     # 2x mul! (BLAS gemm)
  |     |-- _accumulate_rho_jump!(...)        # omega loop: N_w iterations, each with 2x mul!
  |     |-- _apply_precomputed_channel!(...)  # 4x mul! (BLAS gemm) + hermitianize
  |     |-- if step % save_every == 0: trace_distance_h(...)
```

### Pattern 1: BLAS Thread Enablement (THREAD-03/05)

**What:** Wrap `run_thermalize` hot loop in BLAS thread save/restore to allow multi-threaded BLAS during DM evolution.
**When to use:** Always (no downside when DM path runs alone).
**Key insight:** Unlike the trajectory engine which sets BLAS to 1 (because it does Julia-level parallelism across trajectories), `run_thermalize` has a single sequential hot loop where each step's `mul!` calls can benefit from multi-threaded BLAS.

```julia
# Pattern: BLAS thread enablement for DM path
function run_thermalize(...)
    # ... precomputation (already exists) ...

    # Enable multi-threaded BLAS for the DM hot loop
    old_blas = BLAS.get_num_threads()
    try
        # The hot loop runs serially over steps, but each step's mul! calls
        # benefit from multi-threaded BLAS (gemm on dim x dim matrices).
        # No need to change BLAS thread count -- just ensure it's restored.
        for step in 1:num_steps
            # ... existing hot loop body (unchanged) ...
        end
    finally
        BLAS.set_num_threads(old_blas)
    end

    # ... post-loop result construction (unchanged) ...
end
```

**Critical detail:** The `try/finally` ensures BLAS thread count is restored even if an error occurs inside the loop (e.g., user interrupt, numerical failure). This matches the trajectory engine pattern at `trajectories.jl` lines 565-573.

### Pattern 2: Omega-Loop Parallelism with Per-Task Accumulators (THREAD-01/02)

**What:** Parallelize the frequency summation inside `_accumulate_rho_jump!` by partitioning the `energy_labels` across tasks, each with its own scratch buffers, then summing the per-task `rho_jump` accumulators.
**When to use:** When `length(energy_labels)` is large enough to amortize task spawn overhead (empirically: > ~50 frequencies with dim >= 64).

```julia
# Pattern: Per-task accumulator for omega-loop
function _accumulate_rho_jump_threaded!(
    scratch::ThermalizeScratch,
    evolving_dm::Matrix{<:Complex},
    jump::JumpOp,
    ham_or_trott,
    config::Config{Thermalize, EnergyDomain},
    precomputed_data;
    jump_weight_scaling::Real,
)
    energy_labels = precomputed_data.energy_labels
    n_labels = length(energy_labels)
    nt = min(Threads.nthreads(), n_labels)

    if nt <= 1 || n_labels < OMEGA_THREAD_THRESHOLD
        # Fall back to serial
        return _accumulate_rho_jump!(scratch, evolving_dm, jump, ham_or_trott,
                                     config, precomputed_data;
                                     jump_weight_scaling=jump_weight_scaling)
    end

    CT = eltype(evolving_dm)
    dim = size(evolving_dm, 1)
    chunks = _partition_energy_labels(energy_labels, nt)

    # Per-task scratch: each task needs its own rho_jump, jump_oft, sandwich_tmp
    task_scratches = [ThermalizeScratch(CT, dim) for _ in 1:length(chunks)]

    old_blas = BLAS.get_num_threads()
    BLAS.set_num_threads(1)  # Prevent BLAS oversubscription inside threaded loop
    try
        @sync for (idx, chunk_labels) in enumerate(chunks)
            Threads.@spawn _accumulate_rho_jump_chunk!(
                task_scratches[idx], evolving_dm, jump, ham_or_trott,
                config, precomputed_data, chunk_labels;
                jump_weight_scaling=jump_weight_scaling)
        end
    finally
        BLAS.set_num_threads(old_blas)
    end

    # Reduce: sum per-task rho_jump into scratch.rho_jump
    fill!(scratch.rho_jump, 0)
    for ts in task_scratches
        scratch.rho_jump .+= ts.rho_jump
    end

    return nothing
end
```

**Critical details:**
- Each task gets its own `ThermalizeScratch` (same pattern as `_copy_workspace_for_thread` in trajectory engine)
- BLAS is set to 1 during the omega-loop parallelism to prevent oversubscription
- The `evolving_dm` is read-only during `_accumulate_rho_jump!` (no write race)
- Only `rho_jump` accumulation is parallelized; the `_apply_precomputed_channel!` call that mutates `evolving_dm` remains serial

### Pattern 3: BohrDomain Parallel Bucket Iteration (THREAD-04)

**What:** Parallelize the Bohr frequency bucket iteration in `_accumulate_rho_jump!` for BohrDomain.
**Key difference from Energy/Time/Trotter:** The BohrDomain inner loop iterates over Bohr frequency buckets (variable-size sparse index sets), not a uniform energy grid. Each bucket computes `B_{nu_2} = alpha(bohr_freqs, nu_2) * A`, then accumulates `rho_jump += scaled_delta * B_{nu_2} * (rho * A_{nu_2}^dag)`.

```julia
# Pattern: Per-task Bohr bucket accumulation
# Each task processes a subset of Bohr frequency keys
# Per-task scratch provides separate rho_jump, jump_oft, sandwich_tmp
# After @sync, sum per-task rho_jump into the main scratch.rho_jump
```

**Threading considerations for BohrDomain:**
- Bucket sizes are highly variable (diagonal bucket has `dim` entries; most off-diagonal buckets have 1-2 entries)
- The inner accumulation uses manual index loops, not BLAS `mul!` (except the final `mul!` for `rho_jump += B * sandwich_tmp`)
- Load balancing: simple chunk partitioning may leave some tasks idle if bucket sizes are very uneven. For now, use equal-count partitioning of `bohr_keys` -- the per-bucket cost is dominated by the `dim`-sized inner loops which are similar across buckets.

### Anti-Patterns to Avoid

- **Mutating `evolving_dm` from multiple tasks:** The omega-loop reads `evolving_dm` (for `rho * Aw^dag` sandwiches). It must NOT be modified during parallel accumulation. Only `_apply_precomputed_channel!` writes to `evolving_dm`, and it runs after the parallel accumulation completes.
- **Using `threadid()` to index scratch buffers:** Julia PSA (2023) explicitly discourages this. Use per-task allocation with `@spawn`, exactly as the trajectory engine does.
- **Forgetting BLAS thread restoration on error:** Always use `try/finally` -- never bare `BLAS.set_num_threads()` before/after.
- **Setting BLAS threads to 1 when not doing Julia-level parallelism:** The DM hot loop (Pattern 1) benefits from multi-threaded BLAS. Only set BLAS to 1 when doing omega-loop parallelism (Pattern 2/3).
- **Threading the `_precompute_per_jump_channels` loop:** This runs once at startup and iterates over `n_jumps` (typically 3*n_qubits = 12-36). Each iteration calls `_precompute_R` which has its own omega loop. Threading the jump loop would require per-task scratch AND is dominated by `_build_cptp_channel` (which calls `eigen()`). The startup cost is negligible compared to the hot loop; thread the hot loop, not the startup.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| BLAS thread save/restore | Custom thread count tracking | `BLAS.get_num_threads()` + `BLAS.set_num_threads()` + `try/finally` | Already validated in trajectory engine (3 call sites); any custom wrapper would duplicate stdlib |
| Chunk partitioning for omega-loop | Custom partitioning function | Adapt `_partition_trajectories` from `trajectories.jl` (lines 400-413) | Already battle-tested for trajectory parallelism; same pattern applies to energy_labels |
| Per-task scratch allocation | Thread-local storage or global scratch pool | Allocate `ThermalizeScratch(CT, dim)` per task in a vector before `@sync` | Same pattern as `ws_per_task = [_copy_workspace_for_thread(ws) for _ in 1:length(chunks)]` in trajectories.jl |
| Thread count decision logic | Hardcoded thread counts | `nt = min(Threads.nthreads(), n_items)` with a threshold constant | Same guard as trajectory engine: `if ntraj > 1 && Threads.nthreads() > 1` |

**Key insight:** Every threading pattern needed for Phase 41 already exists in the trajectory engine. The DM thermalization path needs the same infrastructure applied to a different parallelism axis (frequency labels instead of trajectories).

## Common Pitfalls

### Pitfall 1: BLAS Thread Oversubscription

**What goes wrong:** If `run_thermalize` enables multi-threaded BLAS (THREAD-03) and then the omega-loop threading (THREAD-01) also leaves BLAS enabled, each spawned task's `mul!` calls will spawn BLAS threads, leading to `Threads.nthreads() * BLAS.get_num_threads()` total threads competing for cores.

**Why it happens:** Multi-threaded BLAS and Julia-level `@spawn` threading are two independent parallelism layers. They must not both be active simultaneously.

**How to avoid:** When doing omega-loop parallelism, set `BLAS.set_num_threads(1)` inside the `try` block. When doing serial DM evolution with BLAS parallelism, leave BLAS at its default. The decision is: "Is the hot loop doing Julia-level parallelism on this iteration?" If yes, BLAS=1. If no, BLAS=default.

**Warning signs:** Wall time increases with thread count; CPU usage exceeds 100% on all cores; `perf stat` shows excessive context switches.

### Pitfall 2: Data Race on scratch.rho_jump

**What goes wrong:** If multiple tasks accumulate into the same `scratch.rho_jump` matrix concurrently, the `+=` operations race and produce incorrect results (silently -- no crash, just wrong numbers).

**Why it happens:** `mul!(scratch.rho_jump, A, B, alpha, 1.0)` reads and writes `scratch.rho_jump`. If two tasks call this on the same matrix simultaneously, the BLAS `gemm!` writes interleave unpredictably.

**How to avoid:** Each task MUST have its own `ThermalizeScratch` instance. The per-task scratches are allocated before `@sync` and summed after. This is the exact pattern used in the trajectory engine (`ws_per_task`).

**Warning signs:** Non-deterministic results when `Threads.nthreads() > 1`; results differ between serial and threaded runs by more than FP accumulation order differences.

### Pitfall 3: BLAS Thread Count Not Restored After Error

**What goes wrong:** If `run_thermalize` throws an error (e.g., NaN in trace distance, user interrupt), the BLAS thread count remains modified, affecting all subsequent computations in the Julia session.

**Why it happens:** Without `try/finally`, an exception bypasses the `BLAS.set_num_threads(old_blas)` call.

**How to avoid:** Always use `try/finally` pattern. The trajectory engine uses this pattern in 3 places (lines 567-573, 686-693, 789-804 of trajectories.jl). Copy it exactly.

**Warning signs:** After a failed `run_thermalize`, subsequent trajectory runs or BLAS operations behave differently (slower or faster than expected).

### Pitfall 4: FP Accumulation Order Differences

**What goes wrong:** When accumulating `rho_jump` in parallel (sum of per-task results) vs serial (single accumulation loop), the floating-point results differ due to different addition order. The tolerance in the success criteria (atol < 1e-10 for multi-step) must account for this.

**Why it happens:** Floating-point addition is not associative. `(a + b) + c != a + (b + c)` in general. With `dim^2` complex values being accumulated across `N_w` frequency iterations, the order difference per step is O(N_w * dim^2 * eps) ~ O(100 * 10000 * 1e-16) ~ O(1e-10). Across `num_steps`, this compounds.

**How to avoid:** Accept that threaded and serial results will differ at O(eps * N_w * dim^2 * num_steps) level. The success criterion already allows atol < 1e-10 for accumulated multi-step. For verification, compare single-step threaded vs serial (should agree within O(1e-13)) and multi-step (should agree within O(1e-10)).

**Warning signs:** Differences larger than 1e-8 suggest a bug, not FP ordering.

### Pitfall 5: Hermitian Jump Half-Grid Partitioning

**What goes wrong:** For Hermitian jumps, the Energy/Time/Trotter omega loops iterate only the half-grid (w_raw <= 0) and explicitly add the mirrored negative-frequency partner. If the energy_labels are partitioned across tasks, each task must handle its own half-grid filter correctly.

**Why it happens:** The half-grid optimization `w_raw > 1e-12 && continue` skips positive frequencies for Hermitian jumps. If labels are naively partitioned, some tasks may get all-positive chunks and do no work, while others get all the work.

**How to avoid:** Either (a) pre-filter energy_labels to the half-grid before partitioning (so all tasks get actual work), or (b) partition the full grid and let each task filter -- simpler but with load imbalance. Option (a) is better for balanced work distribution.

**Warning signs:** Some tasks finish instantly while others take the full time; overall speedup is much less than expected.

### Pitfall 6: Omega-Loop Threading Threshold Too Aggressive

**What goes wrong:** For small systems (dim < 32, few energy labels), the overhead of task spawning and scratch allocation exceeds the parallelism benefit, making threaded execution slower than serial.

**Why it happens:** `Threads.@spawn` has overhead (~1-10 microseconds per spawn). If each frequency iteration takes < 1 microsecond (small matrices), spawning tasks is a net loss.

**How to avoid:** Add a threshold constant (e.g., `const OMEGA_THREAD_THRESHOLD = 50`) below which the serial path is always used. The threshold should be tuned empirically but a conservative default of 50 frequency labels is reasonable. For BLAS threading (THREAD-03), there is no such threshold -- BLAS internally decides whether to thread based on matrix size.

**Warning signs:** Regression in small-system benchmarks when threading is enabled.

## Code Examples

Verified patterns from the existing codebase:

### Existing BLAS Save/Restore Pattern (trajectories.jl)

```julia
# Source: trajectories.jl lines 565-573
old_blas = BLAS.get_num_threads()
BLAS.set_num_threads(1)
try
    @sync for (idx, chunk) in enumerate(chunks)
        Threads.@spawn _run_chunk_no_obs!(
            ws_per_task[idx], psi0, chunk, master_seed, total_time)
    end
finally
    BLAS.set_num_threads(old_blas)
end
```

### Existing Per-Task Workspace Copy (trajectories.jl)

```julia
# Source: trajectories.jl lines 561-563
nt = min(Threads.nthreads(), ntraj)
chunks = _partition_trajectories(1:ntraj, nt)
ws_per_task = [_copy_workspace_for_thread(ws) for _ in 1:length(chunks)]
```

### Existing Partition Utility (trajectories.jl)

```julia
# Source: trajectories.jl lines 400-413
function _partition_trajectories(range::UnitRange{Int}, n_chunks::Int)
    len = length(range)
    n_chunks = min(n_chunks, len)
    base = div(len, n_chunks)
    remainder = rem(len, n_chunks)
    chunks = Vector{UnitRange{Int}}(undef, n_chunks)
    start = first(range)
    for i in 1:n_chunks
        chunk_size = base + (i <= remainder ? 1 : 0)
        chunks[i] = start:(start + chunk_size - 1)
        start += chunk_size
    end
    return chunks
end
```

### Existing Per-Task Result Reduction (trajectories.jl)

```julia
# Source: trajectories.jl lines 576-578
rho_total = sum(ws.scratch.rho_acc for ws in ws_per_task)
rho_result = rho_total ./ ntraj
hermitianize!(rho_result)
```

### Current _accumulate_rho_jump! EnergyDomain (jump_workers.jl)

```julia
# Source: jump_workers.jl lines 453-501
# This is the function whose inner omega-loop would be parallelized.
# Each iteration does:
#   1. oft!(scratch.jump_oft, ...)           -- elementwise, no BLAS
#   2. mul!(scratch.sandwich_tmp, rho, Aw')  -- BLAS gemm (reads rho, writes scratch)
#   3. mul!(scratch.rho_jump, Aw, tmp, ...)  -- BLAS gemm (accumulates into rho_jump)
# The rho read is safe (immutable during accumulation).
# The rho_jump write is the data race risk that per-task accumulators solve.
```

### Current _apply_precomputed_channel! (furnace_utensils.jl)

```julia
# Source: furnace_utensils.jl lines 212-233
# This function contains 4x mul! (BLAS gemm) on dim x dim matrices.
# These are the primary beneficiaries of multi-threaded BLAS (THREAD-03).
# At dim=256 (8 qubits), each gemm is 256x256 complex = well within
# BLAS multi-thread sweet spot.
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `threadid()` indexed buffers | Per-task allocation via `@spawn` return | Julia PSA July 2023 | Prevents data races from task migration |
| `@threads :static` for reductions | `@sync/@spawn` with explicit partitioning | Julia 1.7+ task migration | Composable, no deadlock risk |
| No BLAS thread control | Save/restore via `try/finally` | This codebase v1.2 (Phase 12-19) | Prevents thread count leaks; trajectory engine pattern |

**Julia version notes:**
- Julia 1.11 (current compat target): stable `@spawn`, stable `BLAS.set_num_threads`, dynamic task scheduling
- `Base.Threads` is fully stable and the recommended threading API
- No changes expected in Julia 1.12 that affect this threading model

## Two-Tier Threading Strategy

The phase requirements describe two distinct parallelism levels that must NOT be active simultaneously:

### Tier 1: Multi-Threaded BLAS (THREAD-03, THREAD-05)

- **Benefit:** Every `mul!` call in the hot loop (`_apply_coherent_unitary!`, `_accumulate_rho_jump!`, `_apply_precomputed_channel!`) benefits from BLAS-level parallelism
- **Implementation:** Minimal -- just save/restore BLAS thread count with `try/finally`
- **Matrix sizes that benefit:** At dim=64 (6 qubits), BLAS threading provides modest speedup. At dim=256+ (8+ qubits), significant speedup. At dim < 32, BLAS internally remains single-threaded regardless
- **No code changes to math:** The `mul!` calls are already there; BLAS threads are a runtime configuration

### Tier 2: Omega-Loop Parallelism (THREAD-01, THREAD-02, THREAD-04)

- **Benefit:** The omega-loop in `_accumulate_rho_jump!` iterates `N_w` times (typically 50-500 for Energy/Time/Trotter). Parallelizing across frequencies provides speedup when `N_w * dim^2` work exceeds spawn overhead
- **Implementation:** Significant -- requires per-task scratch allocation, chunk partitioning of energy_labels, and post-accumulation reduction
- **Trade-off:** BLAS must be set to 1 during omega-loop parallelism to avoid oversubscription. So we lose BLAS parallelism on the individual `mul!` calls within the loop, but gain Julia-level parallelism across iterations
- **When beneficial:** When `N_w` is large and `dim` is moderate (so individual `mul!` are cheap but there are many). For large `dim` with few `N_w`, Tier 1 (BLAS threading) is better
- **Decision:** User-configurable or heuristic threshold. Conservative default: serial omega-loop, multi-threaded BLAS. Enable omega-loop threading only when `N_w > OMEGA_THREAD_THRESHOLD`

### Mutual Exclusion

At any given moment during a step:
- Either BLAS is multi-threaded (Tier 1 active, Tier 2 inactive)
- Or Julia tasks are running in parallel with BLAS=1 (Tier 2 active, Tier 1 inactive)

This is the same fundamental constraint the trajectory engine already manages.

## BohrDomain Threading Specifics (THREAD-04)

BohrDomain's `_accumulate_rho_jump!` (jump_workers.jl lines 569-623) has a different structure:

```
for (k, nu_2) in pairs(bohr_keys)           # Outer: iterate Bohr frequency buckets
    @. scratch.jump_oft = alpha(...) * A     # Elementwise (no BLAS)
    fill!(scratch.sandwich_tmp, 0)
    for t in eachindex(is)                   # Inner: sparse bucket indices
        for p in 1:dim                       # Manual accumulation loop
            scratch.sandwich_tmp[p, i] += evolving_dm[p, j] * v
        end
    end
    mul!(scratch.rho_jump, jump_oft, sandwich_tmp, ...)  # One BLAS call per bucket
end
```

**Threading approach for BohrDomain:**
1. Partition `bohr_keys` across tasks (same as partitioning `energy_labels`)
2. Each task gets its own scratch with private `rho_jump`, `jump_oft`, `sandwich_tmp`
3. Set BLAS to 1 during parallel bucket iteration
4. After `@sync`, sum per-task `rho_jump` into main scratch

**BohrDomain-specific considerations:**
- The `sandwich_tmp` computation (`rho * A_{nu_2}^dag`) is the main compute cost per bucket -- it involves a `dim x dim` sparse-times-dense operation
- The final `mul!` is a single BLAS gemm per bucket -- with BLAS=1, this is fast for moderate dim
- Load imbalance: the diagonal bucket (nu=0) has `dim` entries while most off-diagonal buckets have 1-2 entries. However, the `dim`-sized inner loop (`for p in 1:dim`) runs for ALL entries, so per-bucket cost is proportional to `n_entries * dim`
- For the `_precompute_R` BohrDomain (THREAD-02), the same partitioning applies -- identical loop structure, just no `evolving_dm` dependency

## Open Questions

1. **Should omega-loop threading be opt-in via a keyword argument?**
   - What we know: The phase description says "optional omega-loop parallelism" and THREAD-01/02 say "runs multi-threaded"
   - What's unclear: Whether this should be always-on (with a heuristic threshold) or controlled by a user-facing keyword (e.g., `thread_omega::Bool=false`)
   - Recommendation: Use a heuristic threshold internally (no user keyword). If `Threads.nthreads() > 1` AND `length(energy_labels) > OMEGA_THREAD_THRESHOLD`, use threaded omega-loop. Otherwise serial. This avoids API surface expansion while providing automatic speedup for large problems.

2. **What is the optimal omega-loop threading threshold?**
   - What we know: Task spawn overhead is ~1-10 microseconds. A single `mul!` on a 64x64 complex matrix takes ~5-20 microseconds (without BLAS threading). At 100 frequency labels with dim=64, the serial loop takes ~2-4 ms.
   - What's unclear: Exact threshold that balances spawn overhead vs parallelism benefit
   - Recommendation: Start with `OMEGA_THREAD_THRESHOLD = 50` and validate empirically. This is a module-level constant that can be tuned later.

3. **Should per-task ThermalizeScratch allocation be cached across steps?**
   - What we know: Currently, each step would allocate new per-task scratches if omega-loop threading is enabled. Across `num_steps` (thousands), this creates GC pressure.
   - What's unclear: Whether allocation overhead matters compared to the loop work
   - Recommendation: Allocate per-task scratches ONCE before the hot loop (alongside K0s/U_residuals precomputation), then reuse across steps. This mirrors the trajectory engine's `ws_per_task` allocation before the `@sync` block.

4. **How to handle BLAS thread count when run_thermalize is called from user code that has its own BLAS settings?**
   - What we know: The `try/finally` pattern guarantees restoration to the caller's BLAS thread count
   - What's unclear: Whether the user expects BLAS threads to be active during `run_thermalize` or controls it externally
   - Recommendation: Document that `run_thermalize` manages BLAS threads internally and restores the caller's setting on return. Do not change the default BLAS thread count -- just ensure restoration.

## Sources

### Primary (HIGH confidence)
- `src/trajectories.jl` lines 559-574, 677-694, 780-805 -- existing BLAS save/restore + `@sync/@spawn` pattern (3 proven implementations)
- `src/trajectories.jl` lines 154-178 -- `_copy_workspace_for_thread` per-task workspace isolation pattern
- `src/trajectories.jl` lines 396-413 -- `_partition_trajectories` chunk partitioning utility
- `src/furnace.jl` lines 144-250 -- current `run_thermalize` implementation (Phase 39/40 state)
- `src/jump_workers.jl` lines 442-623 -- `_accumulate_rho_jump!` for all 4 domains (the code to be threaded)
- `src/furnace_utensils.jl` lines 212-233 -- `_apply_precomputed_channel!` (BLAS gemm beneficiary)
- `test/test_threading.jl` -- existing trajectory threading tests (BLAS restoration, serial-threaded agreement, determinism)
- `.planning/REQUIREMENTS.md` -- THREAD-01 through THREAD-05 specifications
- `.planning/ROADMAP.md` -- Phase 41 success criteria

### Secondary (MEDIUM confidence)
- [Julia Multi-Threading documentation](https://docs.julialang.org/en/v1/manual/multi-threading/) -- `@threads`, `@spawn`, data race avoidance
- [Julia PSA: Thread-local state is no longer recommended](https://julialang.org/blog/2023/07/PSA-dont-use-threadid/) -- per-task accumulation pattern, avoiding `threadid()`
- [Julia issue #49455](https://github.com/JuliaLang/julia/issues/49455) -- multi-threaded parallel calls to `mul!` slower than serial (confirms need to set BLAS=1 during Julia-level parallelism)
- [ENCCS Julia HPC Multithreading guide](https://enccs.github.io/julia-for-hpc/multithreading/) -- BLAS thread management best practices

### Tertiary (LOW confidence)
- None -- all findings verified with primary or secondary sources

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- no new libraries, uses existing Julia stdlib threading APIs already in the codebase
- Architecture: HIGH -- all threading patterns directly copy existing trajectory engine code (3+ validated implementations); only the parallelism axis changes (frequencies instead of trajectories)
- Pitfalls: HIGH -- all identified from direct code analysis and verified Julia threading documentation; FP tolerance pitfall validated by existing `test_threading.jl` serial-threaded agreement test
- BohrDomain threading: MEDIUM -- BohrDomain's sparse bucket structure is less naturally parallelizable than the dense frequency grid; load balancing may need empirical tuning

**Research date:** 2026-03-01
**Valid until:** 2026-04-01 (stable -- internal threading changes, Julia stdlib threading API is stable in 1.11+)
