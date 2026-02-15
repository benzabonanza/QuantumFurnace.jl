# Domain Pitfalls: Multi-Threaded Trajectory Engine, GNS Comparison, and Adaptive Sampling

**Domain:** Adding multi-threaded trajectory sampling, GNS-vs-KMS convergence experiments, adaptive sampling termination, and convergence tracking to an existing Julia quantum Gibbs sampling package
**Researched:** 2026-02-15
**Confidence:** HIGH (codebase-grounded analysis + verified Julia documentation + community-documented pitfalls)

---

## Critical Pitfalls

Mistakes that cause rewrites, silently wrong results, or major performance regressions.

---

### Pitfall 1: BLAS Thread Oversubscription -- Julia Threads x OpenBLAS Threads = CPU Starvation

**What goes wrong:**
OpenBLAS (Julia's default BLAS) spawns its own thread pool for matrix operations. When you run `julia -t 8` and OpenBLAS also uses 8 threads, every `mul!(ws.psi_tmp, per_op.K0, psi)` call in `step_along_trajectory!` launches 8 BLAS threads, and 8 Julia threads each do this concurrently, creating 64 OS threads competing for 8 cores. The result is *slower* than single-threaded execution. This is documented in [Julia issue #49455](https://github.com/JuliaLang/julia/issues/49455) where multi-threaded `mul!` was 73% slower than serial.

**Why it happens:**
The current codebase uses `mul!` extensively in the trajectory hot path (lines 491-507 of `trajectories.jl`: `mul!(ws.Rpsi, per_op.R, psi)`, `mul!(ws.psi_tmp, per_op.K0, psi)`, `mul!(ws.Rpsi, per_op.U_residual, psi)`). These are matrix-vector products (`gemv!`) which call into BLAS. For n=4 (dim=16), BLAS threading overhead exceeds the computation time. For n=8 (dim=256), each `gemv!` does 256x256 work -- still small enough that BLAS threading hurts more than helps when Julia threads handle the outer parallelism.

**Consequences:**
- Multi-threaded trajectory runs are slower than single-threaded.
- Performance puzzlement leads to removing threading ("it doesn't help"), when the fix is simply `BLAS.set_num_threads(1)`.
- At n=12 (dim=4096), BLAS threading on individual `gemv!` calls *might* help, but only if Julia threads are few. The crossover is system-dependent.

**Prevention:**
1. **Set `BLAS.set_num_threads(1)` at the entry point of the multi-threaded trajectory engine.** This is the standard pattern used by ITensors.jl, Krylov.jl, and other Julia HPC packages.
2. Save and restore the original BLAS thread count: `old_blas_threads = BLAS.get_num_threads(); BLAS.set_num_threads(1); ... ; BLAS.set_num_threads(old_blas_threads)`.
3. For n=12 on cluster, benchmark the crossover: try `BLAS.set_num_threads(1)` with `julia -t 64` vs `BLAS.set_num_threads(4)` with `julia -t 16`. The optimal split depends on the ratio of trajectory-level parallelism to per-trajectory BLAS work.
4. The codebase already has a TODO comment at `jump_workers.jl:464`: `#TODO: test it; set BLAS threads to 1, let julia threads be more.` -- this confirms the author already identified the issue.

**Detection:**
- `htop` shows many more threads than expected.
- Multi-threaded runs are *slower* than single-threaded (anti-scaling).
- `BLAS.get_num_threads()` returns > 1 when Julia is started with `-t N`.

**Phase to address:** Phase 1 (Multi-threaded trajectory engine) -- must be set before any performance benchmarking.

---

### Pitfall 2: Shared Mutable TrajectoryWorkspace Across Threads -- Data Race in the Hot Path

**What goes wrong:**
`TrajectoryFramework` contains a single mutable `ws::TrajectoryWorkspace{T}` (defined at `trajectories.jl:46`) with shared buffers `jump_oft`, `psi_tmp`, and `Rpsi`. The current `run_trajectories` function (line 383) uses a single `fw` and loops sequentially over trajectories. If this loop is naively parallelized with `@threads for trajectory in 1:ntraj`, all threads write to the same `ws.psi_tmp`, `ws.Rpsi`, and `ws.jump_oft` concurrently. This is a data race that produces silently wrong results (corrupted state vectors, wrong jump probabilities, incorrect density matrices).

**Why it happens:**
The `TrajectoryFramework` struct was designed for single-threaded execution. The workspace is embedded in the framework to avoid per-trajectory allocation. When parallelizing, the natural instinct is to parallelize the outer loop without realizing the workspace is shared mutable state.

**Consequences:**
- Corrupted trajectory results that *look* reasonable (density matrix is still approximately Hermitian, trace approximately 1) but converge to wrong steady state.
- Non-deterministic behavior that changes with thread count, timing, and system load.
- Extremely difficult to diagnose because the race condition is probabilistic and the corruption is small per-step.

**Prevention:**
1. **Create per-thread workspaces.** Allocate `nthreads()` copies of `TrajectoryWorkspace` and index by `threadid()` or use `@spawn`-based task parallelism with task-local workspaces.
2. Better pattern: allocate a `Vector{TrajectoryWorkspace}` of length `nthreads()` and pass `ws_pool[threadid()]` into the step function.
3. **Do NOT share `psi` across threads.** Each thread must have its own state vector. The current `psi = copy(psi0)` on line 384 must move inside the parallelized loop, creating per-thread copies.
4. The read-only data in `TrajectoryFramework` (per_operator, jumps, precomputed_data, config) is safe to share -- only `ws` is mutable.
5. Preferred Julia pattern: use `Threads.@spawn` with closures that capture per-task workspace, rather than `@threads` with `threadid()` indexing (which is fragile if tasks migrate between threads).

**Detection:**
- Results change when thread count changes (even with same RNG seed per thread).
- `@assert` on state vector norm fails sporadically.
- ThreadSanitizer (not available for Julia) would catch this, but manual code review is the primary defense.

**Phase to address:** Phase 1 (Multi-threaded trajectory engine) -- core design decision, must be right from the start.

---

### Pitfall 3: Global RNG in step_along_trajectory! -- Non-Reproducible and Non-Thread-Safe Random Sampling

**What goes wrong:**
The current `step_along_trajectory!` calls `rand(1:fw.n_jumps)` and `rand() * total_weight` (lines 483, 517, 627, 661 of `trajectories.jl`) using Julia's global default RNG. When multiple threads call `rand()` concurrently:
1. The global RNG's state is task-local (since Julia 1.7), so concurrent tasks get different streams. But the stream assignment depends on task scheduling order, which is non-deterministic. Result: the same seed produces different results on different runs.
2. Reproducibility across different thread counts is impossible because task creation order changes the parent RNG state (documented in [Julia issue #49522](https://github.com/JuliaLang/julia/issues/49522)).

**Why it happens:**
Julia's `TaskLocalRNG` is designed so that `rand()` is thread-safe (no race conditions), but it is NOT reproducible across different threading configurations. Each spawned task gets a deterministically-derived child RNG from the parent, but the order of task spawning depends on the parallel decomposition, which depends on `nthreads()`.

**Consequences:**
- Cannot reproduce results from a paper by changing thread count (e.g., laptop has 4 threads, cluster has 64).
- Cannot bisect a bug by varying thread count while holding the random stream fixed.
- Test suite becomes flaky when CI runners have varying thread counts.
- The existing use of `StableRNGs` in tests (see `Project.toml` test deps) provides cross-version stability but does NOT solve multi-threading reproducibility.

**Prevention:**
1. **Pass explicit per-trajectory RNG objects** derived from a master seed: `rng_per_traj = [StableRNG(master_seed + i) for i in 1:ntraj]` or `Xoshiro(master_seed + i)`.
2. Modify `step_along_trajectory!` to accept an `rng::AbstractRNG` parameter: `a = rand(rng, 1:fw.n_jumps)` instead of `a = rand(1:fw.n_jumps)`.
3. Seed derivation must be deterministic and independent of thread count. The pattern `seed_i = hash(master_seed, i)` ensures trajectory `i` always uses the same stream regardless of which thread runs it.
4. Do NOT use `Random.seed!(seed)` + `TaskLocalRNG` for multi-threaded code. This seeds only the current task's RNG, not child tasks.
5. For paper-quality results, store the per-trajectory seeds alongside the results so any trajectory can be replayed exactly.

**Detection:**
- Run the same simulation twice with same master seed but different `-t` values. If results differ, reproducibility is broken.
- Check that `run_trajectories(..., ntraj=1000, seed=42)` on `-t 1` gives the same `rho_mean` as on `-t 8`.

**Phase to address:** Phase 1 (Multi-threaded trajectory engine) -- API design must include RNG parameters from the start.

---

### Pitfall 4: False Sharing on Per-Thread Density Matrix Accumulation

**What goes wrong:**
The natural multi-threading pattern accumulates trajectory results into a shared `rho_mean` matrix using atomic operations or a lock: `lock(lk); rho_mean .+= psi * psi'; unlock(lk)`. Even without a lock, if per-thread partial sums are stored in adjacent memory (e.g., `rho_partials = [zeros(dim,dim) for _ in 1:nthreads()]`), the matrices for different threads may share cache lines (64 bytes = 4 ComplexF64 values). When thread 1 writes to `rho_partials[1][end, end]` and thread 2 writes to `rho_partials[2][1, 1]`, false sharing forces cache line invalidation between cores, destroying parallel speedup.

**Why it happens:**
Julia's default allocator does not guarantee cache-line alignment between array allocations. For small matrices (n=4, dim=16, 16x16 ComplexF64 = 4 KB), different per-thread accumulator matrices may land within the same cache line at their boundaries. The effect is most severe when `dim` is small (n=4: 4 KB per matrix, cache line covers ~1.5% of the matrix).

**Consequences:**
- Multi-threaded accumulation shows sub-linear or no speedup despite embarrassingly parallel work.
- Performance varies unpredictably with `nthreads()` and allocation patterns.
- The issue is invisible in profiling (no locks, no contention visible) -- only cache miss counters reveal it.

**Prevention:**
1. **Use per-thread accumulators with a final reduction** instead of writing to shared memory:
   ```julia
   rho_partials = [zeros(ComplexF64, dim, dim) for _ in 1:nthreads()]
   @threads for i in 1:ntraj
       # ... run trajectory ...
       rho_partials[threadid()] .+= psi * psi'
   end
   rho_mean = sum(rho_partials) ./ ntraj
   ```
2. For n>=6 (dim>=64, matrix = 64 KB), false sharing is negligible because each matrix spans many cache lines. Focus prevention effort on n=4 benchmarks.
3. Padding between per-thread allocations (allocate each in a separate array with 64-byte alignment) eliminates the issue but is rarely needed for matrices >= 4 KB.
4. Prefer the `@spawn`-per-batch pattern where each task accumulates its own local `rho_batch` and returns it, then the main thread sums the batch results. This naturally separates memory.

**Detection:**
- `perf stat` (Linux) shows high L1d cache miss rate.
- Speedup plateaus or degrades beyond 2-4 threads for small systems.
- Speedup is fine for n=8 but poor for n=4 with the same code.

**Phase to address:** Phase 1 (Multi-threaded trajectory engine) -- design the accumulation pattern correctly.

---

### Pitfall 5: KMS-vs-GNS "Comparison" Using Different Sigma Values Without Realizing It Changes the Physics

**What goes wrong:**
The KMS transition function (`_pick_transition_kms` in `energy_domain.jl:9-46`) uses a shifted argument `w + beta*sigma^2/2` in the exponential, while the GNS transition function (`_pick_transition_gns` in `energy_domain.jl:64-103`) uses the unshifted argument `w`. This means that for the same `(beta, sigma, a, b)` parameters, the two lines have different effective temperature-dependent behavior. A "comparison" that uses the same sigma for both lines is comparing apples to oranges: the KMS line at `sigma=0.1` and the GNS line at `sigma=0.1` are implementing physically different Lindbladians with different spectral gaps, different mixing times, and different fixed-point approximation errors.

**Why it happens:**
The mathematical difference is by design (KMS-DB uses the shifted gamma, GNS-DB uses the unshifted gamma satisfying KMS condition). But when setting up comparison experiments, it is natural to use the same config parameters for both lines: `LiouvConfig(sigma=0.1, ...)` vs `LiouvConfigGNS(sigma=0.1, ...)`. The sigma appears as the same parameter but its role in the transition function differs.

**Consequences:**
- Paper results that claim "KMS converges faster than GNS" (or vice versa) may be an artifact of comparing at different effective resolutions.
- Referee asks "did you control for the effective smoothing parameter?" and the answer is no.
- Results do not reproduce when using a different Hamiltonian because the sigma-shift effect depends on the spectral gap.

**Prevention:**
1. **Define a "fair comparison" protocol** before running experiments:
   - Option A: **Same sigma** -- accept that the Lindbladians are different and the comparison is "at fixed algorithmic parameter sigma, which line converges faster?"
   - Option B: **Matched spectral gap** -- tune sigma_GNS so that the GNS Lindbladian has the same spectral gap as the KMS Lindbladian. This is fair but expensive (requires computing the spectral gap for each).
   - Option C: **Same computational budget** -- run both for the same wall-clock time (or same number of CPTP map applications) and compare trace distance to Gibbs. This is the most operationally meaningful comparison.
2. **Document which comparison protocol is used** in every experiment.
3. **Plot both the fixed-point approximation error** (`||rho_fixedpoint - gibbs||`) and the **convergence rate** separately. KMS may have a better fixed point (exact with coherent term) but slower mixing, or vice versa.

**Detection:**
- Large discrepancy in mixing times between KMS and GNS that reverses when sigma is changed.
- GNS fixed point is further from Gibbs than KMS (expected, since GNS omits the coherent correction), but the trajectory comparison does not account for this.
- Results flip when switching Hamiltonian.

**Phase to address:** Phase for KMS-vs-GNS experiments -- must define protocol before running any comparison.

---

### Pitfall 6: Comparing KMS-vs-GNS Convergence Without Separating Fixed-Point Error from Mixing Rate

**What goes wrong:**
The KMS Lindbladian (with coherent correction) has the exact Gibbs state as its fixed point. The GNS Lindbladian (without coherent correction, enforced by `ThermalizeConfigGNS` constructor constraint `with_coherent && error(...)`) has an *approximate* fixed point. When measuring "convergence to Gibbs state" as trace distance over time, the KMS line eventually reaches `distance ~ 0` while the GNS line saturates at `distance ~ epsilon_fixedpoint > 0`. This makes KMS look superior, but the mixing *rate* (spectral gap) might actually favor GNS for certain systems.

**Why it happens:**
The GNS construction (CKBG23) satisfies the KMS detailed balance condition in the transition weights but omits the Lamb shift (coherent B term). This means its fixed point is the Gibbs state only in the `sigma -> 0` limit. For finite sigma, the fixed point deviates from Gibbs by `O(sigma)`. The KMS construction (CKG23) adds the coherent correction to achieve exact detailed balance. Comparing convergence to Gibbs mixes two effects: (1) how fast the Lindbladian mixes (spectral gap), and (2) how close the fixed point is to Gibbs (approximation error).

**Consequences:**
- The paper cannot make claims about relative mixing rates if the convergence metric conflates mixing with fixed-point accuracy.
- At small sigma (where both fixed points are close to Gibbs), the comparison is fair for mixing rate. At large sigma, the GNS line is penalized by its worse fixed point.
- Reviewers familiar with the Chen et al. papers will flag this immediately.

**Prevention:**
1. **Separate the two metrics:**
   - **Mixing rate:** Compute the spectral gap of each Lindbladian via `run_lindbladian`. Report `gap_KMS` vs `gap_GNS` directly.
   - **Fixed-point quality:** Compute `||rho_fixedpoint_KMS - gibbs||` and `||rho_fixedpoint_GNS - gibbs||` via eigenanalysis of the Lindbladian.
   - **Trajectory convergence:** Measure trace distance to *each line's own fixed point*, not to the Gibbs state. This isolates the mixing rate from the fixed-point error.
2. Plot convergence as `log(||rho(t) - rho_fixedpoint||)` vs `t`. The slope gives the mixing rate. The y-intercept is the initial distance. The x-intercept is meaningless if measured against Gibbs for GNS.
3. For the paper, include a panel showing `||rho_fixedpoint - gibbs||` vs sigma for both lines. This demonstrates the cost of omitting the coherent correction.

**Detection:**
- GNS convergence curve flattens at a nonzero trace distance while KMS reaches machine precision.
- The "plateau" level for GNS changes with sigma -- this is the fixed-point error, not a sampling artifact.
- Removing the sigma-dependence (setting sigma very small) makes both lines converge to the same point.

**Phase to address:** Phase for KMS-vs-GNS experiments -- analysis methodology, must precede result interpretation.

---

### Pitfall 7: Adaptive Sampling Terminates Prematurely Due to Autocorrelation in Sequential Trajectories

**What goes wrong:**
An adaptive stopping criterion checks whether the standard error of the trajectory-averaged observable has dropped below a threshold: `se = std(observations) / sqrt(n_traj)`. If `se < epsilon`, stop sampling. But consecutive trajectories started from the same initial state `psi0` and evolved for the same total time produce *correlated* final states when the mixing time is not much longer than the total evolution time. The naive standard error estimate assumes i.i.d. samples and underestimates the true uncertainty by a factor of `sqrt(tau_corr)`, where `tau_corr` is the integrated autocorrelation time.

**Why it happens:**
In the current `run_trajectories` (line 383-396), each trajectory starts from `psi0` and evolves independently. If the system has not fully mixed by `total_time`, the final states cluster around the transient-evolved state rather than sampling the full steady-state distribution. The samples ARE independent (different random streams), but they are all biased toward the same transient regime. The standard error of the mean correctly reflects sampling noise, but the *mean itself* is biased. The adaptive criterion detects low variance (all trajectories are close to each other) and terminates, not realizing that the trajectories are close to each other because they are all close to the wrong answer.

**Consequences:**
- Adaptive sampling terminates with a "converged" flag but the result is far from the steady state.
- The reported uncertainty is small but does not include the systematic bias from insufficient mixing.
- More trajectories do NOT fix this -- the bias is a property of `total_time`, not `ntraj`.

**Prevention:**
1. **Never use convergence of the trajectory average as evidence of convergence to the steady state.** The two are completely different:
   - "Trajectory average has converged" = sampling noise is small = you have enough trajectories for THIS total_time.
   - "Result is close to steady state" = total_time is long enough = you need to check against a known steady-state metric.
2. **Use an observable-based convergence criterion** instead of a statistical one: track `||rho(t) - rho(t-Delta)|| / ||rho(t-Delta)||` over macro time windows. Convergence to the steady state means this ratio drops below a threshold. This is what `run_thermalization` already does (checking distance to Gibbs at each step, line 143-150 in `furnace.jl`).
3. **For adaptive trajectory count at fixed total_time:** use batch means with gap. Run `B` batches of `M` trajectories each. Compute `rho_batch_b` for each batch. The adaptive criterion is `std([||rho_batch_b - rho_running_mean||^2 for b]) / sqrt(B) < epsilon`. This gives a proper estimate of the mean's uncertainty.
4. **For adaptive total_time at fixed trajectory count:** run trajectories to time `T`, compute `rho_avg(T)`. Then extend all trajectories to time `2T` and compute `rho_avg(2T)`. If `||rho_avg(2T) - rho_avg(T)|| < threshold`, the total time is sufficient. This is a time-doubling protocol.

**Detection:**
- Convergence criterion triggers very quickly (e.g., after 50 trajectories) for a system that is known to have a long mixing time.
- Increasing total_time changes the "converged" result significantly.
- The "converged" result depends on `psi0` (if truly at steady state, it should not).

**Phase to address:** Phase for adaptive sampling -- core algorithm design.

---

## Moderate Pitfalls

Mistakes that cause significant debugging time or subtle performance issues.

---

### Pitfall 8: Per-Thread Workspace Copies Balloon Memory at n=12

**What goes wrong:**
For n=12 qubits (dim=4096), a single `TrajectoryWorkspace` contains:
- `jump_oft`: 4096x4096 ComplexF64 = 256 MB
- `psi_tmp`: 4096 ComplexF64 = 64 KB
- `Rpsi`: 4096 ComplexF64 = 64 KB

Plus each `PerOperatorKraus` contains R, K0, U_residual (each 4096x4096 = 256 MB) and optionally U_B (256 MB). With 36 jump operators (3 Paulis x 12 sites), the per-operator data alone is `36 * 4 * 256 MB = 36 GB`. This is shared read-only across threads.

But if each of `T` threads needs its own workspace, that adds `T * 256 MB` for `jump_oft` alone. With 64 threads on a 512 GB cluster node, that is 16 GB of workspace. Manageable, but only if:
1. The precomputed data is shared (not copied per thread).
2. The NUFFTPrefactors 3D array (`dim x dim x n_energy_labels`) is shared (not copied).

The NUFFTPrefactors at n=12: `4096 * 4096 * n_labels * 16 bytes`. With `num_energy_bits=12` giving ~4096 energy labels (after truncation, maybe ~500), that is `4096 * 4096 * 500 * 16 = 128 GB`. This alone may exceed node memory.

**Prevention:**
1. **Audit memory before parallelizing.** Print total precomputed data size before entering the trajectory loop. For n=12, this is a critical path: if the prefactors do not fit in memory, multi-threading is moot.
2. Share all read-only data (`per_operator`, `precomputed_data`, `jumps`, `config`) via a single `TrajectoryFramework` instance. Only `ws` needs per-thread copies.
3. For the `jump_oft` buffer in the workspace: this is a temporary matrix used during the dissipative jump scan. At n=12, the scan iterates over energy labels and writes to `jump_oft` using `@. ws.jump_oft = jump.in_eigenbasis * pref`. This is the hot allocation. Consider using smaller buffers if the scan can be restructured (e.g., compute `||A_w psi||^2` without materializing the full A_w matrix).
4. The `SharedArray` support already in `nufft.jl` (line 51-53) is for `Distributed` parallelism, not `Threads`. For `Threads`, regular arrays are automatically shared. Do not accidentally create `nprocs() * nthreads()` copies.

**Detection:**
- `OutOfMemoryError` when starting multi-threaded run at n=12.
- Memory usage grows linearly with thread count when it should be constant (precomputed data) + linear (workspaces only).
- Swap usage on cluster node indicates memory pressure.

**Phase to address:** Phase 1 (Multi-threaded trajectory engine) -- memory budget must be established before scaling to n=12.

---

### Pitfall 9: GC Stop-the-World Pauses Destroy Multi-Thread Scaling

**What goes wrong:**
Julia's garbage collector is stop-the-world: when GC triggers, ALL threads pause until collection completes. If trajectory threads allocate frequently (e.g., creating temporary matrices, closures, or type-unstable intermediate values), GC pauses increase with thread count because more threads = more allocation = more frequent GC. This is documented in [Julia issue #33033](https://github.com/JuliaLang/julia/issues/33033) where multi-threaded allocation-heavy workloads showed significant slowdown as thread count increased.

**Why it happens:**
The trajectory step function (`step_along_trajectory!`) is *designed* to be allocation-free (all workspace is pre-allocated). But subtle allocations can creep in:
- `rand(1:fw.n_jumps)` -- allocates a `UnitRange` on each call.
- `@. ws.jump_oft = jump.in_eigenbasis * pref` -- if types are not inferred, creates temporary arrays.
- `mul!(ws.Rpsi, ws.jump_oft', psi)` -- the adjoint `ws.jump_oft'` creates a lazy `Adjoint` wrapper, which itself does not allocate the array but creates a small heap object.
- Closure captures in the transition function (e.g., `w -> exp(-beta * w / 2)` from `pick_transition`) -- each call to the closure is fine, but if the closure itself is not inferred, it boxes.

**Consequences:**
- Scaling from 1 to 8 threads gives 3x speedup instead of expected 7x.
- GC time (visible via `@time` or `GC.gc_stats()`) grows from 1% to 20% of wall time.
- Performance degrades non-uniformly: some trajectory batches take much longer than others (those that trigger GC).

**Prevention:**
1. **Profile allocations of the hot path before parallelizing.** Use `@allocated step_along_trajectory!(psi, fw)` (warmup + measure, as already done in `test_allocation.jl`). If allocations per step > 0, find and eliminate them.
2. Replace `rand(1:fw.n_jumps)` with `rand(rng, 1:fw.n_jumps)` using an explicit RNG to avoid potential task-local RNG lookup overhead.
3. Ensure all workspace buffers are `isbitstype` or concretely typed. Check with `@code_warntype step_along_trajectory!(psi, fw)`.
4. Consider `GC.enable(false)` during the trajectory loop and `GC.gc()` after completion, but ONLY if the loop is confirmed allocation-free. If it allocates, disabling GC causes unbounded memory growth.
5. The existing `test_allocation.jl` tests guard against allocation regressions -- ensure these tests cover the multi-threaded hot path as well.

**Detection:**
- `@time` reports significant GC time (> 5% of total).
- `@allocated` in the step function returns > 0 bytes.
- Scaling curve shows diminishing returns beyond 4 threads.

**Phase to address:** Phase 1 (Multi-threaded trajectory engine) -- allocation audit before parallelization.

---

### Pitfall 10: Thread-Unsafe FINUFFT Plan Execution

**What goes wrong:**
The NUFFT prefactors are precomputed in `_prepare_oft_nufft_prefactors` (`nufft.jl`) using a FINUFFT plan object. This plan is used only during precomputation (not in the trajectory loop). However, if precomputation is ever parallelized (e.g., precomputing for KMS and GNS configs concurrently), FINUFFT plan objects are NOT thread-safe -- concurrent `finufft_exec!` calls on different plans may share internal FFTW wisdom, causing race conditions.

**Why it happens:**
FINUFFT uses FFTW internally, and FFTW's plan creation/execution shares global state (the "wisdom" cache). Multiple concurrent FINUFFT executions can corrupt this shared state.

**Consequences:**
- Silent data corruption in precomputed NUFFT prefactors.
- Non-deterministic NaN or incorrect values in the prefactor matrix.
- Trajectory results are wrong in ways that look like physics bugs, not threading bugs.

**Prevention:**
1. **Keep precomputation single-threaded** (current behavior is correct -- `nthreads=1` on line 109 of `furnace_utensils.jl`).
2. If parallelizing precomputation for multiple configs (e.g., KMS and GNS in parallel), serialize the FINUFFT calls or use separate FFTW wisdom sessions.
3. The trajectory loop does NOT call FINUFFT (it uses precomputed prefactors via `_prefactor_view`), so this is only a concern during setup.

**Detection:**
- Incorrect trajectory results when precomputing KMS and GNS configs concurrently.
- FFTW-related segfaults or assertion failures.

**Phase to address:** Phase for KMS-vs-GNS experiments -- relevant if parallelizing experiment setup.

---

### Pitfall 11: Convergence Tracking Accumulates O(ntraj * dim^2) Memory

**What goes wrong:**
If convergence tracking stores the density matrix `rho_avg` at each checkpoint (e.g., every 100 trajectories), and there are 100,000 trajectories with 1000 checkpoints, this creates 1000 density matrices of size `dim x dim`. For n=8 (dim=256): `1000 * 256 * 256 * 16 bytes = 1 GB`. For n=12: `1000 * 4096 * 4096 * 16 = 256 GB`. This silently exhausts memory before the simulation completes.

**Why it happens:**
The natural convergence tracking pattern is: "store a snapshot every K trajectories so I can plot the convergence curve later." This works for small systems but does not scale.

**Consequences:**
- `OutOfMemoryError` during long runs, losing all partially computed results.
- Swap thrashing makes the simulation appear to hang.
- Results are lost because no intermediate output was saved to disk.

**Prevention:**
1. **Track scalar convergence metrics, not full density matrices.** Store `trace_distance(rho_avg, rho_target)` or `fidelity(rho_avg, rho_target)` at each checkpoint. This is O(1) per checkpoint.
2. **Keep only the current and previous `rho_avg`** in memory for computing convergence rate. Overwrite the previous at each checkpoint.
3. For post-hoc analysis, save checkpointed `rho_avg` to disk (BSON/JLD2) rather than keeping in memory.
4. For convergence curve plotting, store `(ntraj, trace_distance)` pairs: `Vector{Tuple{Int, Float64}}`, which is negligible memory.
5. The existing `run_thermalization` already follows this pattern: it stores `distances_to_gibbs::Vector{Float64}` (scalar per step), not density matrix snapshots.

**Detection:**
- Memory usage grows linearly during trajectory accumulation.
- GC becomes frequent and slow as heap grows.
- For n>=8, the simulation slows dramatically partway through.

**Phase to address:** Phase for convergence tracking -- data structure design.

---

### Pitfall 12: Batch Size Effect on Convergence Estimate Gives False Precision

**What goes wrong:**
Adaptive stopping uses batch means: group trajectories into batches of size `B`, compute `rho_batch` for each, then estimate uncertainty from the inter-batch variance. If `B` is too small (e.g., `B=10`), each batch `rho_batch` is very noisy, making the inter-batch variance large. This gives a *correct but conservative* stopping criterion (safe, but wastes computation). If `B` is too large (e.g., `B=10000`), there are few batches (e.g., 10 batches for 100K trajectories), and the inter-batch variance estimate has high uncertainty (few degrees of freedom in the chi-squared distribution), giving false precision -- the estimated standard error may be too small by a factor of 2-3.

**Why it happens:**
The number of batches `K` determines the degrees of freedom for the variance estimate: `var(batch_means)` has `K-1` degrees of freedom. With `K=5`, the 95% confidence interval for the variance spans a factor of ~6. The estimated standard error can be off by sqrt(6) ~ 2.4x in either direction. For `K=30`, the factor drops to ~1.5, which is acceptable.

**Consequences:**
- With too few batches: premature stopping with underestimated uncertainty (false confidence).
- With too many small batches: correct but very slow convergence detection (wasted computation).
- Published uncertainty bounds may be too tight by 2-3x.

**Prevention:**
1. **Use at least 20-30 batches** for the adaptive stopping criterion. This gives ~20-30 degrees of freedom and a variance estimate within ~50% of the true value.
2. **Batch size formula:** `B = max(ntraj_min / 30, 100)` where `ntraj_min` is the minimum expected trajectory count. Start with `B=100` and adjust.
3. **Report effective sample size (ESS)** alongside the batch estimate: `ESS = ntraj * var_naive / var_batch`. If `ESS << ntraj`, the batch means are correlated (possible if trajectories share some hidden state, though they should not).
4. Use the Student-t correction: the stopping criterion should use `t_alpha(K-1) * se_batch` rather than `z_alpha * se_batch` to account for the finite number of batches.

**Detection:**
- Batch variance estimate changes significantly when batch size is doubled (should be relatively stable).
- Re-running the same experiment gives significantly different uncertainty estimates.
- `ESS / ntraj < 0.5` indicates a problem with the batching.

**Phase to address:** Phase for adaptive sampling -- statistical methodology.

---

## Minor Pitfalls

Issues that cause confusion or small inefficiencies.

---

### Pitfall 13: mul! on Adjoint Views May Allocate Unexpectedly

**What goes wrong:**
In the dissipative branch of `step_along_trajectory!`, `mul!(ws.Rpsi, ws.jump_oft', psi)` creates a lazy `Adjoint` wrapper for `ws.jump_oft'`. While the wrapper itself is small, some BLAS dispatch paths for `Adjoint` inputs may allocate intermediate arrays. This is harmless for single-threaded code but creates GC pressure in multi-threaded code (see Pitfall 9).

**Prevention:**
- Use `mul!(ws.Rpsi, adjoint(ws.jump_oft), psi)` explicitly and verify with `@allocated` that it does not allocate.
- Alternatively, precompute `jump_oft_dag` as a separate buffer and fill it with `copy!(ws.jump_oft_dag, ws.jump_oft); conj!(ws.jump_oft_dag)` -- but this wastes a buffer for a minor optimization.
- Profile first: the allocation may already be zero (Julia's BLAS dispatch for `gemv!` with `Adjoint` is typically allocation-free for dense matrices).

**Phase to address:** Phase 1 (Allocation audit) -- check during allocation profiling.

---

### Pitfall 14: `threadid()` Indexing is Unsafe with Task Migration (Julia >= 1.7)

**What goes wrong:**
Using `Threads.threadid()` to index into per-thread workspace arrays assumes tasks stay on the same thread. Since Julia 1.7, tasks can migrate between threads during yield points. If a task starts on thread 3 (uses workspace 3), yields (e.g., GC pause, I/O), and resumes on thread 5, it now uses workspace 5, which may be in use by another task. This causes the same data-race as Pitfall 2.

**Prevention:**
1. Use `@spawn`-based parallelism where each task creates and owns its workspace in a closure. The workspace is stack-local to the task, not indexed by thread ID.
2. If using `@threads`, ensure the loop body never yields. For trajectory stepping, this is likely the case (pure computation, no I/O), but GC pauses can trigger migration.
3. Pin tasks to threads using `@threads :static` (Julia >= 1.7), which prevents migration. This is the simplest fix but limits scheduler flexibility.
4. Preferred pattern:
   ```julia
   @sync for batch in batches
       @spawn begin
           ws_local = TrajectoryWorkspace(CT, dim)
           psi_local = Vector{ComplexF64}(undef, dim)
           for i in batch
               # ... use ws_local, psi_local ...
           end
       end
   end
   ```

**Phase to address:** Phase 1 (Multi-threaded trajectory engine) -- architecture decision.

---

### Pitfall 15: Convergence Comparison Using Different Initial States

**What goes wrong:**
When comparing KMS-vs-GNS convergence rates, using a different initial state for each run (e.g., random pure states with different seeds) introduces variance in the convergence curve that masks the actual rate difference. The trajectory-averaged convergence is `E[||rho(t) - rho_ss||]`, which depends on the initial state.

**Prevention:**
1. **Use the same initial state for all comparison runs.** The standard choice is `|0...0>` (computational zero state) or the maximally mixed state `I/dim`.
2. **Use the maximally mixed state** for convergence rate comparison: it is rotationally symmetric (no basis-dependent artifacts) and its initial distance to Gibbs is `1 - 1/dim` (known analytically), making convergence curves comparable.
3. **Avoid using random pure states** for comparison experiments unless the comparison is averaged over initial states (which requires many runs).

**Phase to address:** Phase for KMS-vs-GNS experiments -- experimental protocol.

---

### Pitfall 16: FINUFFT nthreads Parameter Interacts with Julia Threading

**What goes wrong:**
The NUFFT precomputation (`_prepare_oft_nufft_prefactors` in `nufft.jl:59`) uses `nthreads=1` for the FINUFFT plan. If this is changed to use more FINUFFT threads while Julia threads are also active, FINUFFT's internal OpenMP parallelism conflicts with Julia's threading in the same way as BLAS (Pitfall 1). Additionally, FINUFFT linked against FFTW may use FFTW's own thread pool, creating a three-way thread contention: Julia threads x BLAS threads x FFTW threads.

**Prevention:**
1. Keep FINUFFT `nthreads=1` when using Julia multi-threading (current behavior is correct).
2. If precomputation is a bottleneck, parallelize across energy labels using Julia threads rather than FINUFFT internal parallelism.
3. Document the thread budget: `total_OS_threads <= num_cores` where `total = julia_threads * BLAS_threads * FINUFFT_threads`.

**Phase to address:** Phase 1 -- inherited from current design, no change needed unless precomputation is parallelized.

---

## Phase-Specific Warnings

| Phase Topic | Likely Pitfall | Mitigation |
|-------------|---------------|------------|
| Multi-threaded trajectory engine | BLAS oversubscription (P1), shared workspace race (P2), global RNG (P3) | Set `BLAS.set_num_threads(1)`, per-task workspace, explicit per-trajectory RNG |
| Multi-threaded trajectory engine | False sharing on accumulation (P4), GC pressure (P9) | Per-batch accumulation with final reduction, allocation audit before parallelizing |
| Multi-threaded trajectory engine | Memory at n=12 (P8), task migration (P14) | Memory budget audit, `@threads :static` or `@spawn`-per-batch pattern |
| GNS trajectory path | Sigma interpretation (P5), fixed-point vs mixing rate (P6) | Define comparison protocol, separate metrics |
| KMS-vs-GNS experiments | Fair comparison methodology (P5, P6), initial state (P15) | Same initial state, same computational budget, separate fixed-point error from mixing rate |
| Adaptive sampling | Premature convergence (P7), batch size effects (P12) | Observable-based convergence + statistical stopping, 20-30 batches minimum |
| Convergence tracking | Memory from snapshots (P11), autocorrelation (P7) | Scalar metrics only, time-doubling protocol for mixing verification |
| Precomputation parallelism | FINUFFT thread safety (P10, P16) | Keep FINUFFT single-threaded, serialize precomputation |

---

## Sources

### Julia Threading and BLAS Interaction
- [Julia Multi-Threading Documentation](https://docs.julialang.org/en/v1/manual/multi-threading/)
- [Multi-threaded `mul!` slower than serial (Julia issue #49455)](https://github.com/JuliaLang/julia/issues/49455)
- [ITensors.jl Multithreading Guide (BLAS.set_num_threads pattern)](https://itensor.github.io/ITensors.jl/dev/Multithreading.html)
- [Julia BLAS thread count defaults (issue #33409)](https://github.com/JuliaLang/julia/issues/33409)
- [Julia Threads vs BLAS threads (Discourse)](https://discourse.julialang.org/t/julia-threads-vs-blas-threads/8914)

### RNG and Reproducibility
- [TaskLocalRNG vs Xoshiro performance (Discourse)](https://discourse.julialang.org/t/why-is-tasklocalrng-faster-than-xoshiro-with-multiple-threads/74577)
- [Random.seed! reproducibility issue (#49522)](https://github.com/JuliaLang/julia/issues/49522)
- [Reproducible multithreaded Monte Carlo (Discourse)](https://discourse.julialang.org/t/reproducible-multithreaded-monte-carlo-task-local-random/35269)
- [StableRNGs.jl](https://github.com/JuliaRandom/StableRNGs.jl)

### GC and Allocation
- [Multi-threaded allocation benchmark shows GC slowdown (issue #33033)](https://github.com/JuliaLang/julia/issues/33033)
- [GC and threading performance (Discourse)](https://discourse.julialang.org/t/garbage-collection-and-threading/107250)
- [Julia Memory Management Documentation](https://docs.julialang.org/en/v1/manual/memory-management/)

### False Sharing
- [False sharing in multi-threading (Julia blog)](https://blog.jling.dev/blog/false_share/)
- [Julia for HPC: Multithreading](https://enccs.github.io/julia-for-hpc/multithreading/)

### Quantum Trajectory Methodology
- [Monte Carlo wave-function method: robust algorithm and convergence (Abdelhafez et al. 2019)](https://arxiv.org/abs/1803.08589)
- [QuTiP Monte Carlo Solver documentation](https://qutip.org/docs/4.5/guide/dynamics/dynamics-monte.html)
- [QuantumToolbox.jl Monte Carlo Solver](https://qutip.org/QuantumToolbox.jl/stable/users_guide/time_evolution/mcsolve)

### Adaptive MCMC and Convergence
- [Convergence diagnostics for MCMC (Roy 2019)](https://arxiv.org/pdf/1909.11827)
- [Adaptive MCMC survey (Liang et al. 2020)](https://asp-eurasipjournals.springeropen.com/articles/10.1186/s13634-020-00675-6)

### Codebase References
- `src/trajectories.jl` -- TrajectoryWorkspace, step_along_trajectory!, run_trajectories
- `src/structs.jl` -- ThermalizeConfig, ThermalizeConfigGNS, TrajectoryFramework
- `src/energy_domain.jl` -- _pick_transition_kms, _pick_transition_gns (sigma shift difference)
- `src/jump_workers.jl:464` -- existing TODO: "set BLAS threads to 1"
- `src/nufft.jl` -- NUFFTPrefactors, memory scaling
- `src/furnace_utensils.jl:109` -- FINUFFT nthreads=1
- `src/furnace.jl:101` -- run_thermalization RNG parameter
- `test/test_allocation.jl` -- existing allocation regression tests

---
*Pitfalls research for: QuantumFurnace.jl -- Multi-threaded trajectory engine, GNS comparison, adaptive sampling*
*Researched: 2026-02-15*
*Milestone: Multi-threaded trajectory engine + KMS-vs-GNS experiments*
