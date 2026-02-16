# Phase 13: Multi-Threaded Trajectory Engine - Research

**Researched:** 2026-02-16
**Domain:** Julia multi-threading for embarrassingly parallel quantum trajectory sampling
**Confidence:** HIGH

## Summary

Phase 13 implements multi-threaded trajectory sampling in QuantumFurnace.jl. The core problem is embarrassingly parallel: each trajectory evolves an independent state vector `psi` through `step_along_trajectory!` calls using a private `TrajectoryWorkspace` and RNG, reading from a shared `TrajectoryFramework`. Phase 12 already separated the mutable workspace from the read-only framework and made all functions accept explicit `ws` and `rng` arguments, so the threading infrastructure is ready.

The implementation requires three coordinated changes: (1) a new `run_trajectories_threaded` function (or a `threaded=true` kwarg on the existing `run_trajectories`) that distributes N trajectories across Julia threads using `Threads.@spawn` with per-task workspace/RNG/psi, (2) BLAS thread management via `BLAS.set_num_threads(1)` save/restore around the parallel section, and (3) deterministic seeding using `Xoshiro(master_seed + trajectory_id)` so that each trajectory's random stream is independent of thread count and scheduling. The per-thread density matrix accumulators are reduced after all tasks complete.

**Primary recommendation:** Use `@sync`/`Threads.@spawn` with per-task closures (not `@threads` with `threadid()` indexing) to avoid task migration issues. Seed each trajectory with `Xoshiro(master_seed + traj_id)` for bitwise reproducibility given the same thread count. Set `BLAS.set_num_threads(1)` before and restore after. Accumulate density matrices per-task locally, then reduce.

## Standard Stack

### Core

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| `Base.Threads` | Julia 1.11+ | `@spawn`, `@sync`, `nthreads()`, `Atomic` | Julia's built-in threading primitives; composable task-based parallelism since Julia 1.3, mature since 1.7+ |
| `Random` (stdlib) | Julia 1.11+ | `Xoshiro(seed)`, `AbstractRNG` | Per-trajectory deterministic RNG; Xoshiro256++ is Julia's default PRNG, fast with small memory footprint |
| `LinearAlgebra.BLAS` | Julia 1.11+ | `BLAS.set_num_threads(1)`, `BLAS.get_num_threads()` | Thread count control for OpenBLAS to prevent oversubscription |

### Supporting

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| `StableRNGs` | 1.x (test dep) | Cross-Julia-version reproducible RNG | Only if bitwise reproducibility across Julia versions is needed; not required for same-version reproducibility |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| `@sync`/`@spawn` per-batch | `@threads :static` | `@threads :static` is simpler but pins tasks to threads, preventing load balancing; `@spawn` is more flexible and avoids `threadid()` fragility |
| `@threads :static` | `@threads :dynamic` (default) | `:dynamic` may reassign iterations to different threads, making `threadid()`-based workspace indexing unsafe; `:static` guarantees stable thread assignment but cannot be nested |
| `Xoshiro(master_seed + traj_id)` | `TaskLocalRNG` with parent seeding | TaskLocalRNG seeding depends on task spawn order, which may vary with thread count; explicit Xoshiro per trajectory is fully independent of scheduling |

## Architecture Patterns

### Current State (After Phase 12)

```
TrajectoryFramework{T,D}    (READ-ONLY during stepping)
  domain, jumps, ham_or_trott, config, precomputed_data
  per_operator::Vector{PerOperatorKraus{T}}
  n_jumps, delta, delta_eff, alpha

TrajectoryWorkspace{T}       (MUTABLE, per-trajectory)
  jump_oft::Matrix{T}   # dim x dim scratch
  psi_tmp::Vector{T}    # dim scratch
  Rpsi::Vector{T}       # dim scratch
  rho_acc::Matrix{T}    # dim x dim density matrix accumulator

step_along_trajectory!(psi, fw, ws, rng)  # 4-arg, explicit ws/rng
run_trajectories(...)                      # serial, returns TrajectoryResult
```

### Target State (After Phase 13)

```
run_trajectories(...)
  |
  +--> ntraj <= 1 or nthreads() == 1?
  |      YES: serial path (current code, unchanged)
  |
  +--> Multi-threaded path:
         1. old_blas = BLAS.get_num_threads()
         2. BLAS.set_num_threads(1)
         3. Partition trajectories into chunks (one per @spawn task)
         4. @sync for chunk in chunks
              @spawn begin
                ws_local = TrajectoryWorkspace(CT, dim)
                psi_local = copy(psi0)
                rng_local = Xoshiro(master_seed + first_traj_in_chunk)
                for traj_id in chunk
                    rng_traj = Xoshiro(master_seed + traj_id)
                    copyto!(psi_local, psi0)
                    _evolve_along_trajectory!(psi_local, fw, ws_local, rng_traj, total_time)
                    _accumulate_density_matrix!(ws_local.rho_acc, psi_local)
                end
              end
            end
         5. Reduce: rho_mean = sum(ws.rho_acc for ws in all_workspaces) / ntraj
         6. BLAS.set_num_threads(old_blas)
         7. return TrajectoryResult(rho_mean, ntraj, master_seed, ...)
```

### Pattern 1: Per-Task Workspace with @spawn

**What:** Each spawned task creates its own `TrajectoryWorkspace`, `psi` vector, and per-trajectory `Xoshiro` RNG inside a closure. The shared `TrajectoryFramework` is captured by reference (read-only).

**When to use:** Always for the multi-threaded trajectory engine.

**Example:**
```julia
# Source: Julia docs + ParallelMCWF.jl pattern + QuantumFurnace Phase 12 architecture
function _run_trajectories_threaded(fw, psi0, ntraj, master_seed, total_time)
    CT = eltype(psi0)
    dim = length(psi0)
    nt = Threads.nthreads()

    # Partition trajectories into chunks, one per thread
    chunks = _partition(1:ntraj, nt)

    # Per-task workspace + accumulator
    ws_per_task = [TrajectoryWorkspace(CT, dim) for _ in 1:length(chunks)]

    @sync for (task_idx, chunk) in enumerate(chunks)
        Threads.@spawn begin
            ws = ws_per_task[task_idx]
            psi = copy(psi0)
            for traj_id in chunk
                rng = Random.Xoshiro(master_seed + traj_id)
                copyto!(psi, psi0)
                _evolve_along_trajectory!(psi, fw, ws, rng, total_time)
                _accumulate_density_matrix!(ws.rho_acc, psi)
            end
        end
    end

    # Reduce
    rho_total = sum(ws.rho_acc for ws in ws_per_task)
    rho_mean = rho_total ./ ntraj
    hermitianize!(rho_mean)
    return rho_mean
end
```

### Pattern 2: Deterministic Seeding (master_seed + trajectory_id)

**What:** Each trajectory gets `Xoshiro(master_seed + traj_id)` regardless of which thread runs it. This makes the random stream for trajectory `i` depend only on `master_seed` and `i`, not on thread count or scheduling.

**When to use:** Always for reproducible multi-threaded runs.

**Key insight:** The seed derivation `master_seed + traj_id` is simple, deterministic, and produces independent Xoshiro streams (Xoshiro256++ has excellent independence across close seeds). No need for hash-based seed derivation or explicit stream splitting.

**Caveat:** Bitwise reproducibility requires the same `nthreads()`. Different thread counts produce different chunk assignments, and while each trajectory's RNG stream is identical, the accumulation order of density matrices differs. Since floating-point addition is not associative, `sum(rho_1 + rho_2 + ... + rho_N)` can differ in the last few ULP depending on grouping. For practical purposes (atol ~1e-14), this is negligible, but if strict bitwise matching is required across thread counts, the reduction must use a canonical ordering.

### Pattern 3: BLAS Thread Save/Restore

**What:** Save the current BLAS thread count, set it to 1, run the parallel section, restore the original count.

**When to use:** At the entry point of any multi-threaded trajectory execution.

**Example:**
```julia
# Source: ITensors.jl, Krylov.jl, Julia HPC best practices
function with_single_blas_thread(f)
    old = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        return f()
    finally
        BLAS.set_num_threads(old)
    end
end
```

### Pattern 4: Observable Accumulation in Threaded Mode

**What:** When `observables !== nothing`, each task also accumulates measurement data locally (`mean_data_local`), then all are summed after the parallel section.

**When to use:** When `run_trajectories` is called with `observables` parameter in multi-threaded mode.

**Key difference from serial:** In serial mode, a single `mean_data` matrix is accumulated across all trajectories. In threaded mode, each task has its own `mean_data_local` and the final result is `sum(mean_data_locals) / ntraj`.

### Anti-Patterns to Avoid

- **Using `threadid()` to index workspace arrays:** Tasks can migrate between threads (Julia >= 1.7), making `threadid()` unreliable. Use per-task closures with `@spawn` instead.
- **Locking during density matrix accumulation:** Using a lock to protect shared `rho_acc` serializes the accumulation step and destroys parallelism. Use per-task local accumulators and reduce after.
- **Using `Random.seed!()` for per-thread seeding:** `Random.seed!()` modifies the task-local RNG, not a named RNG object. It cannot produce independent streams across tasks spawned from the same parent.
- **Cloning `TrajectoryFramework` per thread:** The framework contains large read-only data (per-operator Kraus matrices: ~18 MB for n=8). Only the workspace (~600 KB for n=8) needs per-task copies.
- **Disabling GC during parallel section without confirming zero allocations:** If any code in the hot loop allocates, `GC.enable(false)` causes unbounded memory growth.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Task partitioning | Manual chunk assignment loop | Simple `_partition` helper or `Iterators.partition` | Off-by-one errors in chunk boundaries; `Iterators.partition` handles remainders correctly |
| Thread-safe accumulation | Lock-based shared accumulator | Per-task local accumulators + post-hoc reduction | Locks serialize the hot path; local accumulation is allocation-free and contention-free |
| Independent RNG streams | Custom stream-splitting / jumping | `Xoshiro(master_seed + traj_id)` | Xoshiro256++ has provably independent streams for nearby seeds; no need for explicit jump-ahead |
| BLAS thread management | Manual `set_num_threads` calls without restore | `try/finally` pattern or helper function | Forgetting to restore BLAS threads after an error leaves the process in degraded state |

**Key insight:** The threading pattern for embarrassingly parallel Monte Carlo is well-established. The novel part is integrating it with the existing `TrajectoryFramework`/`TrajectoryWorkspace` separation, not inventing new concurrency primitives.

## Common Pitfalls

### Pitfall 1: BLAS Thread Oversubscription

**What goes wrong:** OpenBLAS spawns its own thread pool. With `julia -t 8` and BLAS using 8 threads, each `mul!(ws.psi_tmp, per_op.K0, psi)` call in `step_along_trajectory!` launches 8 BLAS threads. 8 Julia threads x 8 BLAS threads = 64 OS threads competing for 8 cores. Result: slower than single-threaded.

**Why it happens:** The trajectory step uses `mul!` for matrix-vector products (gemv). For dim=16 (n=4) to dim=256 (n=8), these are too small for BLAS threading to help -- kernel launch overhead dominates.

**How to avoid:** `BLAS.set_num_threads(1)` before entering the parallel section. Restore afterward with `try/finally`.

**Warning signs:** Multi-threaded runs are slower than serial. `htop` shows many more threads than expected.

### Pitfall 2: False Sharing on Adjacent Accumulators

**What goes wrong:** Per-task density matrix accumulators allocated in a tight `Vector{Matrix}` may share cache lines at their boundaries, causing false sharing that destroys parallel scaling for small matrices (n=4, dim=16, matrix = 4 KB).

**Why it happens:** Julia's allocator does not guarantee cache-line alignment between array allocations.

**How to avoid:** Use the `@spawn`-per-batch pattern where each task allocates its workspace independently (in the spawned closure), naturally separating memory. For n>=6 (dim>=64, matrix = 64 KB), false sharing is negligible.

**Warning signs:** Scaling plateaus beyond 2-4 threads for n=4 but works fine for n=8.

### Pitfall 3: Floating-Point Reduction Order Breaks Bitwise Reproducibility

**What goes wrong:** Even with deterministic per-trajectory seeds, different thread counts produce different chunk-to-task assignments. The density matrices are accumulated in chunk order within each task, then reduced across tasks. Since FP addition is not associative, `(a + b) + c != a + (b + c)` at ULP level. Result: bitwise different `rho_mean` across thread counts.

**Why it happens:** The mathematical result is identical to ~1e-14 relative error, but strict `==` comparison fails.

**How to avoid:** For the success criterion "bitwise identical given same master seed and thread count," this is automatically satisfied because the same chunk assignment produces the same accumulation order. For cross-thread-count comparison, use `isapprox` with `atol=1e-13`. The success criteria specify "same thread count," so this pitfall does not block phase completion but should be documented.

**Warning signs:** Tests using exact equality (`==`) on `rho_mean` across different thread counts fail.

### Pitfall 4: Observable Path Memory Scaling

**What goes wrong:** With observables, each task needs a local `mean_data` matrix (n_obs x n_saves). For n=8 with 17 observables and 10000 save points, this is 17 x 10000 x 8 bytes = 1.3 MB per task. With 64 threads, ~85 MB total. Manageable, but must be accounted for.

**Why it happens:** The observable accumulation cannot share a single `mean_data` across threads (data race).

**How to avoid:** Pre-allocate per-task `mean_data` arrays alongside the workspace. Include them in the memory budget calculation.

**Warning signs:** Unexpected memory growth with high thread counts in the observable path.

### Pitfall 5: GC Stop-the-World Pauses

**What goes wrong:** Julia's GC is stop-the-world. If trajectory steps allocate (even small amounts), GC frequency increases with thread count, destroying scaling.

**Why it happens:** The step function is designed to be allocation-free, but subtle allocations can creep in (e.g., `Adjoint` wrappers, closure boxing).

**How to avoid:** Run `@allocated step_along_trajectory!(psi, fw, ws, rng)` to verify zero allocations before parallelizing. The existing `test_allocation.jl` guards against allocation regressions but does not cover the trajectory step function itself -- a new allocation test for `step_along_trajectory!` should be added.

**Warning signs:** `@time` reports significant GC time (>5% of total). Scaling curve shows diminishing returns beyond 4 threads.

## Code Examples

### Complete Multi-Threaded run_trajectories Pattern

```julia
# Source: Codebase analysis + Julia threading docs + community patterns
function run_trajectories(jumps, config, psi0, hamiltonian;
    trotter=nothing, total_time=config.mixing_time, delta=config.delta,
    ntraj=1, observables=nothing, save_every=1, seed=nothing)

    # ... existing setup code (build framework, etc.) ...

    actual_seed = seed === nothing ? Int(rand(Random.RandomDevice(), UInt64) >> 1) : seed

    if ntraj > 1 && Threads.nthreads() > 1
        return _run_trajectories_threaded(fw, psi0, ntraj, actual_seed,
            total_time, observables, save_every)
    else
        return _run_trajectories_serial(fw, psi0, ntraj, actual_seed,
            total_time, observables, save_every)
    end
end

function _run_trajectories_threaded(fw, psi0, ntraj, master_seed,
    total_time, observables, save_every)

    CT = eltype(psi0)
    dim = length(psi0)
    nt = min(Threads.nthreads(), ntraj)

    # BLAS thread control
    old_blas = BLAS.get_num_threads()
    BLAS.set_num_threads(1)

    try
        chunks = _partition_trajectories(1:ntraj, nt)
        ws_per_task = [TrajectoryWorkspace(CT, dim) for _ in 1:nt]

        if observables === nothing
            # No-observable path
            @sync for (idx, chunk) in enumerate(chunks)
                Threads.@spawn _run_chunk_no_obs!(
                    ws_per_task[idx], fw, psi0, chunk, master_seed, total_time)
            end
            rho_total = sum(ws.rho_acc for ws in ws_per_task)
        else
            # Observable path -- each task also accumulates measurements
            # ... similar pattern with per-task mean_data ...
        end

        rho_mean = rho_total ./ ntraj
        hermitianize!(rho_mean)
        return TrajectoryResult(rho_mean, ntraj, master_seed, nothing, nothing)
    finally
        BLAS.set_num_threads(old_blas)
    end
end

function _run_chunk_no_obs!(ws, fw, psi0, chunk, master_seed, total_time)
    psi = copy(psi0)
    for traj_id in chunk
        rng = Random.Xoshiro(master_seed + traj_id)
        copyto!(psi, psi0)
        _evolve_along_trajectory!(psi, fw, ws, rng, total_time)
        _accumulate_density_matrix!(ws.rho_acc, psi)
    end
end
```

### Chunk Partitioning Helper

```julia
function _partition_trajectories(range, n_chunks)
    len = length(range)
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

### Deterministic Seeding Verification Test

```julia
@testset "Deterministic multi-threaded results" begin
    # Run twice with same seed and thread count
    result1 = run_trajectories(..., ntraj=100, seed=42)
    result2 = run_trajectories(..., ntraj=100, seed=42)
    @test result1.rho_mean == result2.rho_mean  # bitwise identical
end
```

### BLAS Thread Management Test

```julia
@testset "BLAS thread restoration" begin
    old = BLAS.get_num_threads()
    result = run_trajectories(..., ntraj=10, seed=42)
    @test BLAS.get_num_threads() == old  # restored after run
end
```

### Performance Comparison Test

```julia
@testset "Threading speedup" begin
    if Threads.nthreads() >= 4
        t_serial = @elapsed run_trajectories(..., ntraj=100, seed=42)
        # Note: need to measure threaded time properly
        t_threaded = @elapsed run_trajectories(..., ntraj=100, seed=42)
        @test t_threaded < t_serial  # threaded is faster (no regression)
    end
end
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Workspace embedded in framework | Workspace passed explicitly | Phase 12 (2026-02-15) | Enables per-thread workspace without framework duplication |
| Global `rand()` in stepping | Explicit `rng::AbstractRNG` | Phase 12 (2026-02-15) | Enables per-trajectory deterministic RNG |
| `threadid()` workspace indexing | `@spawn` with per-task closures | Julia 1.7+ best practice | Avoids task migration issues |
| `@threads` for embarrassingly parallel | `@sync`/`@spawn` for chunked batches | Julia community consensus ~2022+ | Better composability, no nesting restrictions |
| Shared accumulator with lock | Per-task local accumulator + reduce | Standard Monte Carlo pattern | Eliminates contention, enables linear scaling |

## Implementation Considerations

### Memory Budget (n=8, dim=256)

| Component | Size | Shared? | Count |
|-----------|------|---------|-------|
| TrajectoryFramework (per_operator) | ~18 MB (12 ops x 4 matrices x 256^2 x 16 bytes) | Shared (1 copy) | 1 |
| NUFFT prefactors | ~500 MB (256^2 x ~500 labels x 16 bytes) | Shared (1 copy) | 1 |
| TrajectoryWorkspace per task | ~1.1 MB (2 x 256^2 + 2 x 256) x 16 bytes | Per-task | nthreads |
| psi vector per task | ~4 KB (256 x 16 bytes) | Per-task | nthreads |
| **Total per-task overhead** | **~1.1 MB** | | |
| **Total for 64 threads** | **~70 MB** | | |

Memory is dominated by the shared precomputed data, not per-thread state. This is the correct architecture.

### Performance Expectations

For n=4 (dim=16):
- Each trajectory step: ~microseconds (small gemv)
- 100 steps x 100 trajectories = ~1 ms serial
- Expected speedup with 4 threads: ~3-3.5x (some overhead from task creation)

For n=8 (dim=256):
- Each trajectory step: ~100 microseconds (dim=256 gemv is larger)
- 1000 steps x 1000 trajectories = ~100 seconds serial
- Expected speedup with 8 threads: ~6-7x (good scaling, BLAS at 1 thread)

### API Design Decision

**Option A:** Add `threaded::Bool=false` kwarg to existing `run_trajectories`:
- Pro: Single entry point, simple API
- Con: Internal branching complexity

**Option B:** New `run_trajectories_threaded(...)` function:
- Pro: Clean separation, easier to test independently
- Con: API duplication

**Recommendation:** Option A with internal dispatch. The high-level API should be:
```julia
run_trajectories(...; ntraj=1000, seed=42)
# Automatically uses threads if nthreads() > 1 and ntraj > 1
```

No new keyword needed. The function internally checks `Threads.nthreads()` and dispatches. This keeps the API minimal and matches the user's expectation: "I set Julia threads, things go faster."

### Testing Strategy

1. **Determinism test:** Same seed + same nthreads => bitwise identical results (run twice, compare `rho_mean`)
2. **BLAS restoration test:** `BLAS.get_num_threads()` same before and after `run_trajectories`
3. **Correctness test:** Multi-threaded result matches serial result (same seed, `isapprox` with tight tolerance)
4. **Allocation test:** `step_along_trajectory!` allocates zero bytes (guards against GC regression)
5. **Performance test:** Multi-threaded faster than serial for ntraj>=100 with nthreads>=4 (CI may not have 4 threads; guard with `if Threads.nthreads() >= 4`)

### Thread Count Sensitivity

Tests must handle single-threaded CI environments gracefully. All multi-threading tests should be wrapped in `if Threads.nthreads() > 1` guards. The serial fallback path must remain functional and tested independently.

## Open Questions

1. **Automatic vs explicit threading**
   - What we know: The success criteria state "distributes N trajectories across available Julia threads." This implies automatic detection.
   - What's unclear: Should there be a way to force serial execution even with multiple threads? (e.g., for debugging)
   - Recommendation: Auto-detect with `nthreads()`, but document that `julia -t 1` forces serial mode. No explicit kwarg needed for Phase 13.

2. **Observable path complexity**
   - What we know: The `observables !== nothing` code path in `run_trajectories` has a more complex accumulation loop (saves measurements at intervals).
   - What's unclear: Whether the observable path should be threaded in Phase 13 or deferred.
   - Recommendation: Thread both paths (no-observable and observable). The observable path is a straightforward extension of the no-observable pattern: each task gets a local `mean_data` matrix in addition to `ws.rho_acc`, reduced after the parallel section.

3. **Allocation test for step function**
   - What we know: `test_allocation.jl` tests B_bohr, B_time, B_trotter, and `_jump_contribution!` but does NOT test `step_along_trajectory!`.
   - What's unclear: Whether `step_along_trajectory!` is truly zero-allocation on all code paths.
   - Recommendation: Add an allocation regression test for `step_along_trajectory!` before implementing threading. If it allocates, fix first.

## Sources

### Primary (HIGH confidence)
- **QuantumFurnace.jl codebase** -- Direct analysis of `src/trajectories.jl` (761 lines), `src/structs.jl` (312 lines), `src/QuantumFurnace.jl` (95 lines), `src/furnace.jl` (172 lines), `src/kraus.jl` (15 lines), `src/furnace_utensils.jl` (136 lines), `test/test_workspace_independence.jl` (91 lines), `test/test_allocation.jl` (131 lines)
- **Phase 12 RESEARCH.md and VERIFICATION.md** -- Workspace separation complete, 4-arg step signature verified, all 246 tests pass
- **Julia Multi-Threading Documentation** -- https://docs.julialang.org/en/v1/base/multi-threading/ -- `@threads`, `@spawn`, `@sync`, `nthreads()`, scheduling options
- **Julia Random Documentation** -- https://docs.julialang.org/en/v1/stdlib/Random/ -- `Xoshiro`, `AbstractRNG`, TaskLocalRNG behavior, per-task seeding
- **Prior milestone research** -- `.planning/research/FEATURES.md`, `.planning/research/PITFALLS.md` -- threading patterns, BLAS interaction, memory budget, pitfalls catalog

### Secondary (MEDIUM confidence)
- **Julia Discourse: Reproducible multithreaded Monte Carlo** -- https://discourse.julialang.org/t/reproducible-multithreaded-monte-carlo-task-local-random/35269 -- Per-task RNG patterns, practical approaches
- **ParallelMCWF.jl** -- https://github.com/Z-Denis/ParallelMCWF.jl -- `:threads` parallelism for MCWF trajectories in QuantumOptics.jl
- **QuantumToolbox.jl mcsolve** -- https://qutip.org/QuantumToolbox.jl/stable/users_guide/time_evolution/mcsolve -- EnsembleThreads pattern for parallel trajectory Monte Carlo
- **Julia HPC Multithreading** -- https://enccs.github.io/julia-for-hpc/multithreading/ -- BLAS oversubscription, `BLAS.set_num_threads(1)` pattern
- **ITensors.jl Multithreading Guide** -- BLAS.set_num_threads pattern used by production HPC Julia packages

### Tertiary (LOW confidence)
- None -- all findings verified against official docs or codebase analysis

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- Julia stdlib only, no new dependencies needed
- Architecture: HIGH -- Phase 12 laid the groundwork; pattern is well-established in Julia quantum computing ecosystem
- Pitfalls: HIGH -- Comprehensive pitfall catalog from `.planning/research/PITFALLS.md` already identifies all major risks (BLAS oversubscription, workspace sharing, RNG reproducibility, false sharing, GC pressure)
- Code examples: HIGH -- Patterns derived from codebase analysis and verified Julia documentation

**Research date:** 2026-02-16
**Valid until:** Indefinite (Julia threading primitives are stable; no external dependency drift)
