# Domain Pitfalls: v2.1 Speedup & Mixing Time Features

**Domain:** Adding per-jump precomputation, multi-threaded frequency loops, multi-threaded BLAS control, mixing time extrapolation, and save_every to an existing Lindbladian simulator (QuantumFurnace.jl)
**Researched:** 2026-03-01
**Confidence:** HIGH (grounded in direct codebase analysis + Julia threading documentation + numerical analysis literature)

**Relationship to prior research:** This document covers pitfalls specific to the v2.1 milestone. The v2.0 PITFALLS.md covered codebase restructure risks (abstract type boxing, BSON serialization, closure capture). Those pitfalls remain relevant but are not repeated here. This document covers NEW pitfalls introduced by adding performance and mixing time features to the existing v2.0 system.

---

## Critical Pitfalls

Mistakes that cause silent numerical errors, race conditions, or fundamentally wrong results.

---

### CRIT-01: BLAS.set_num_threads Is Global State -- Cannot Be Scoped Per Code Path

**What goes wrong:**
The existing trajectory engine calls `BLAS.set_num_threads(1)` before spawning Julia threads and restores it after (see `trajectories.jl` lines 500-508). This works because trajectory sampling is the ONLY code path running at that time. The v2.1 milestone introduces a second threading pattern: parallelizing the omega-loop within `run_thermalize` (DM evolution). The problem is that `BLAS.set_num_threads()` sets a **process-global** variable -- there is ONE OpenBLAS thread pool shared by all Julia threads. You cannot have "BLAS uses 8 threads for the DM engine" and "BLAS uses 1 thread for trajectory tasks" simultaneously in the same Julia process.

**Why it matters for v2.1:**
The design calls for three different BLAS threading strategies in the same codebase:
1. **Trajectory engine**: `BLAS.set_num_threads(1)` -- Julia threads handle parallelism, each doing small BLAS calls (existing, works)
2. **DM thermalize (no omega threading)**: `BLAS.set_num_threads(N)` -- single Julia thread, let BLAS handle parallelism (new feature)
3. **DM thermalize (with omega threading)**: `BLAS.set_num_threads(1)` -- Julia threads parallelize the omega loop, each doing small BLAS calls (new feature)

The danger: if a user calls `run_thermalize` and then `run_trajectories` in the same session, the BLAS thread count from one call persists into the next unless explicitly reset. The existing try/finally pattern in trajectories.jl handles this correctly for that path, but the NEW DM thermalize code must follow the same discipline.

**Consequences:**
- If DM thermalize sets `BLAS.set_num_threads(8)` and an error occurs before restore, subsequent trajectory calls with Julia threading will oversubscribe by 8x, causing ~10x slowdown.
- If omega-loop threading forgets to set `BLAS.set_num_threads(1)`, each Julia thread's `mul!` call will internally spawn BLAS threads, causing thread oversubscription (N_julia * N_blas threads competing).
- OpenBLAS with thread contention can produce **wrong numerical results** (not just slowness) in versions before 0.3.7 without USE_LOCKING=1. Julia ships OpenBLAS with locking enabled, but this is still a serialization bottleneck.

**Prevention:**
- Wrap ALL threading entry points in a `try/finally` that saves and restores `BLAS.get_num_threads()`. The existing pattern in trajectories.jl (lines 500-508) is the model:
  ```julia
  old_blas = BLAS.get_num_threads()
  BLAS.set_num_threads(1)
  try
      @sync for ...
          Threads.@spawn ...
      end
  finally
      BLAS.set_num_threads(old_blas)
  end
  ```
- For the "let BLAS handle parallelism" path (no Julia threading), do NOT change `BLAS.set_num_threads` at all -- just let the user's default apply. Only set it when Julia threads are active.
- Add a test: run trajectory engine, check `BLAS.get_num_threads()` is restored. Then run DM engine, check again. The existing `test_threading.jl` "BLAS thread restoration" testset is the model.
- Document clearly: "The omega-loop threaded DM path requires BLAS single-threaded, same as trajectories."

**Detection:**
- The existing `test_threading.jl` "BLAS thread restoration" test catches leaked BLAS state for trajectories. Write the equivalent for the DM engine.
- Performance benchmark: if `run_thermalize` is slower than expected despite having multiple Julia threads, suspect BLAS oversubscription.

**Which phase:** Must be addressed in the omega-loop threading phase. Design decision at the start.

---

### CRIT-02: Per-Jump Precomputed Eigendecomposition Becomes Stale Under Numerical Drift

**What goes wrong:**
The existing code computes `_build_cptp_channel(R_a, delta)` per-jump at trajectory workspace construction time (trajectories.jl lines 112-122). This calls `eigen(Hermitian(S))` to perform the PSD guard: negative eigenvalues of the residual matrix S are clamped to zero, then `U_residual = sqrt(diag(clamped_eigenvalues)) * eigvecs'` is computed. This precomputation is correct for the trajectory engine because delta is fixed and R is state-independent.

For the DM thermalize engine, the situation is different. The current `_finalize_kraus_step!` in jump_workers.jl (lines 172-192) calls `_build_cptp_channel(scratch.R, delta)` **every step**, recomputing the eigendecomposition of S from the current accumulated R. If v2.1 precomputes the eigendecomposition once (like the trajectory engine does), the key question is: does R depend on the evolving state?

**Analysis of the existing code:**
Looking at `_jump_contribution!` for `Config{Thermalize, ...}` (jump_workers.jl lines 194-290, 292-367, 369-440): the R accumulation uses `scratch.R` which is **reset to zero** at the start of each jump contribution (`fill!(scratch.R, 0)` at lines 219, 316, 395). The R matrix is built from `rate^2 * L'L` terms where L depends on jump operators and transition rates -- NOT on the evolving density matrix. This means R IS state-independent and CAN be precomputed.

**However**, there is a subtlety: in the BohrDomain path (lines 194-290), the inner loop builds `scratch.jump_oft = alpha(...) * jump.in_eigenbasis` and uses it for both R and rho_jump. The `rho_jump` part DOES depend on the evolving state, but R does not. These two computations are interleaved in the current code. Precomputing R requires untangling R from rho_jump.

**The real danger:**
The PSD guard (`max.(eig.values, 0.0)`) handles numerical noise where S should theoretically be PSD but has eigenvalues at -1e-16 due to floating-point arithmetic. If R is precomputed once and used for many steps:
1. The clamping is correct for the precomputed R.
2. But if delta changes between steps (e.g., adaptive delta), the precomputed U_residual is WRONG because `S = (2*alpha - delta)*R - alpha^2*R^2` depends on delta through `alpha = 1 - sqrt(1 - delta)`.
3. Even with fixed delta, if there are multiple jump operators and the code changes from "one R per step" to "per-jump R precomputed", the residual S per-jump might have different numerical properties than the summed R over all jumps in the current code.

**Prevention:**
- Verify that delta is truly fixed for the entire DM thermalize run (it is -- `config.delta` is a const field). This makes precomputation safe.
- When precomputing per-jump, verify that `S_a = (2*alpha - delta)*R_a - alpha^2*R_a^2` has eigenvalues >= -eps for each individual R_a. Add an assertion: `@assert minimum(eigen(Hermitian(S_a)).values) > -1e-10 "Per-jump S_a has unexpectedly negative eigenvalues"` during workspace construction.
- After precomputing, validate with a regression test: run one step with precomputed vs. recomputed channel and verify bitwise-identical results (or within machine epsilon).
- The memory cost of storing per-jump precomputed data is `N_jumps * (3 * dim^2)` complex matrices (K0, U_residual, R). For 3 qubits (dim=8) with 6 jumps: 6 * 3 * 64 * 16 bytes = 18 KB. For 6 qubits (dim=64) with 12 jumps: 12 * 3 * 4096 * 16 bytes = 2.4 MB. Manageable. For 10+ qubits this starts to matter.

**Detection:**
- Trace distance to Gibbs state diverging instead of converging after switching to precomputed channels.
- CPTP violation: `tr(rho) != 1.0` or `eigvals(rho)` having negative values beyond numerical noise after many steps.
- The existing `test_cptp.jl` should catch CPTP violations.

**Which phase:** Per-jump precomputation phase. Must include a regression test comparing precomputed vs. recomputed results.

---

### CRIT-03: Race Conditions in Omega-Loop Accumulator Matrices

**What goes wrong:**
The omega-loop in `_jump_contribution!` for the Thermalize path accumulates into `scratch.R` and `scratch.rho_jump` across frequency iterations. If this loop is parallelized with `Threads.@threads` or `@sync/@spawn`:

```julia
# WRONG: concurrent writes to shared scratch.R
@threads for w in energy_labels
    # ... compute L(w), rate(w) ...
    scratch.R .+= rate^2 * LdagL   # RACE CONDITION
    mul!(scratch.rho_jump, ...)     # RACE CONDITION
end
```

Matrix `.+=` is NOT atomic. Each element addition is a separate read-modify-write operation. Two threads writing to the same matrix element will lose one update.

**This is different from the trajectory threading:**
The trajectory engine avoids this by giving each thread its OWN workspace (`_copy_workspace_for_thread`). Each thread accumulates into its own `rho_acc`, then results are summed after all threads complete. The omega-loop threading is different because all frequency iterations contribute to the SAME R and rho_jump matrices within a single DM step.

**Prevention strategy -- per-thread partial accumulators:**
```julia
# CORRECT: each task accumulates into its own R, then merge
R_per_task = [zeros(CT, dim, dim) for _ in 1:ntasks]
rho_jump_per_task = [zeros(CT, dim, dim) for _ in 1:ntasks]

@sync for (task_idx, freq_chunk) in enumerate(chunks)
    Threads.@spawn begin
        for w in freq_chunk
            # accumulate into R_per_task[task_idx], rho_jump_per_task[task_idx]
        end
    end
end

# Merge (single-threaded, safe)
R_total = sum(R_per_task)
rho_jump_total = sum(rho_jump_per_task)
```

**Do NOT use Threads.threadid() for indexing:**
Julia tasks can migrate between OS threads. Using `threadid()` to index into pre-allocated buffers is wrong because two tasks could observe the same `threadid()` simultaneously. Julia's own PSA (July 2023) explicitly warns against this pattern. Use task-local storage or indexed by task identity instead.

**Memory cost of per-task accumulators:**
For `ntasks` tasks with dim-by-dim complex matrices: `ntasks * 2 * dim^2 * 16` bytes (R and rho_jump). With 8 tasks and dim=64: 8 * 2 * 4096 * 16 = 1 MB. With dim=256: 8 * 2 * 65536 * 16 = 16 MB. Acceptable for typical problem sizes.

**Additional subtlety -- scratch.jump_oft:**
The current code reuses `scratch.jump_oft` across frequency iterations: `@. scratch.jump_oft = jump.in_eigenbasis * nufft_prefactor_matrix`. Each thread needs its OWN `jump_oft` buffer too, plus any other scratch matrices used within the loop body (`LdagL`, `sandwich_tmp`). Essentially, each task needs a full `ThermalizeScratch` worth of buffers.

**Detection:**
- Non-deterministic trace distance results across runs with the same seed (race condition symptom).
- Trace distance plateauing at a non-zero value (lost accumulation updates).
- Compare threaded omega-loop result against serial loop result with `isapprox(..., atol=1e-14)`. Any disagreement beyond FP accumulation order differences indicates a race.

**Which phase:** Omega-loop threading phase. Must be the core design pattern.

---

### CRIT-04: False Sharing When Per-Task Accumulators Are Contiguous in Memory

**What goes wrong:**
If per-task accumulator matrices are allocated as views into a single large array (e.g., `big_buffer = zeros(dim, dim, ntasks)` with `R_per_task[i] = @view big_buffer[:, :, i]`), adjacent matrices share cache lines. When thread 1 writes to the last elements of its R matrix and thread 2 writes to the first elements of its R matrix, the CPU invalidates the shared cache line for both cores, causing cache thrashing.

On modern CPUs, cache lines are 64 bytes (Intel/AMD x86) or 128 bytes (Apple M-series). A ComplexF64 is 16 bytes, so 4 elements fit in an x86 cache line and 8 in an M-series cache line. The last row of one matrix and the first row of the next are almost certainly in the same cache line.

**This is specific to the omega-loop accumulator pattern:**
The trajectory engine does not have this problem because per-thread workspaces are independently allocated Julia Arrays (each has its own heap allocation, naturally cache-line-separated).

**Prevention:**
- Allocate each per-task matrix independently: `[zeros(CT, dim, dim) for _ in 1:ntasks]`. Julia's allocator naturally aligns to 16-byte boundaries, and individual heap allocations are typically cache-line-separated.
- Do NOT use a single 3D array with views into slices. The independent allocation pattern used by `_copy_workspace_for_thread` (trajectories.jl line 164) is correct.
- For very small matrices (dim <= 4), the entire matrix fits in 1-2 cache lines and false sharing is unavoidable. At these sizes, the omega-loop has too few iterations to benefit from threading anyway, so disable threading below a size threshold.

**Detection:**
- Threading speedup is less than expected (e.g., 1.5x with 8 threads instead of 4x+).
- Performance profiling shows high L1/L2 cache miss rates.
- Empirically: benchmark threaded vs. serial, compare against roofline.

**Which phase:** Omega-loop threading phase. An implementation detail, not a design decision.

---

### CRIT-05: Nested mul!/BLAS Calls Within Julia Threads Without BLAS Single-Threading

**What goes wrong:**
The omega-loop body contains `mul!` calls (BLAS GEMM):
```julia
mul!(scratch.LdagL, scratch.jump_oft', scratch.jump_oft)     # L'L
mul!(scratch.sandwich_tmp, evolving_dm, scratch.jump_oft')    # rho * L'
mul!(scratch.rho_jump, scratch.jump_oft, scratch.sandwich_tmp, ...)  # L * (rho * L')
```

If Julia threads are used for the omega-loop and `BLAS.set_num_threads` is NOT set to 1, each `mul!` call will internally launch BLAS threads. With N Julia threads each launching M BLAS threads, you get N*M total threads competing for CPU resources. This is called "oversubscription" and causes severe performance degradation (worse than serial due to context switching overhead).

**Specific danger for small matrices:**
For typical QuantumFurnace problem sizes (dim = 8 to 64), BLAS GEMM on matrices this small is actually SLOWER with multi-threading because the overhead of thread creation/synchronization dominates the computation. OpenBLAS has an internal threshold below which it runs single-threaded anyway, but this threshold varies by OpenBLAS version and is not guaranteed.

**Prevention:**
- Before the omega-loop's threaded section, call `BLAS.set_num_threads(1)`.
- After the threaded section, restore the original count.
- This is the SAME pattern as the trajectory engine (trajectories.jl lines 500-508). Extract it into a shared helper:
  ```julia
  function with_serial_blas(f)
      old = BLAS.get_num_threads()
      BLAS.set_num_threads(1)
      try
          return f()
      finally
          BLAS.set_num_threads(old)
      end
  end
  ```
- Add a threshold: only use Julia threading for the omega-loop when `length(energy_labels) >= min_freqs_for_threading` (suggest min_freqs_for_threading = 16-32). Below this, the serial loop is faster due to threading overhead.

**Detection:**
- The new DM thermalize threaded path is slower than the serial path despite having multiple threads.
- CPU utilization is 100% across all cores but throughput is low (oversubscription symptom).
- `perf stat` or similar shows high context switch counts.

**Which phase:** Omega-loop threading phase. Must be implemented together with the threading itself.

---

### CRIT-06: Hermitianize Happening Before vs After Merge Changes Numerical Result

**What goes wrong:**
The current code calls `hermitianize!(scratch.R)` after the omega-loop completes (jump_workers.jl line 285 for Bohr, line 362 for Energy, line 435 for Time/Trotter). This enforces Hermiticity of the accumulated R matrix. When threading the omega-loop:

- **Option A: Hermitianize each per-task R, then sum.** Result: `sum(hermitianize.(R_task_i))`
- **Option B: Sum per-task R, then hermitianize.** Result: `hermitianize(sum(R_task_i))`

These give the same mathematical result but different floating-point results due to rounding. Hermitianizing before summation rounds each partial accumulator independently, then sums rounded values. Hermitianizing after summation sums the raw values, then rounds once. The difference is O(ntasks * dim^2 * eps), which for ntasks=8, dim=64, eps=1e-16 gives ~3e-12 -- above the 1e-13 tolerance used in the serial-threaded agreement test.

**Prevention:**
- Use Option B (hermitianize after merge). This matches the serial code's behavior where hermitianize happens once after the full omega-loop.
- Update the serial-threaded agreement test to use a tolerance that accounts for floating-point accumulation order differences: `atol = ntasks * dim^2 * eps(Float64)`. The existing trajectory threading test uses `atol=1e-13` for this reason (test_threading.jl line 89).
- Do NOT aim for bitwise identical serial/threaded results -- FP addition is not associative. Aim for `isapprox` with a tight but realistic tolerance.

**Detection:**
- Serial-threaded agreement tests failing with small errors (1e-13 to 1e-11 range).
- Non-reproducible trace distance differences between serial and threaded DM evolution.

**Which phase:** Omega-loop threading phase. Must decide the hermitianize placement as part of the merge strategy.

---

## Moderate Pitfalls

Issues that cause performance problems or incorrect mixing time estimates but not silent data corruption.

---

### MOD-01: Exponential Fit to Trace Distance Curve Has Structural Model Mismatch

**What goes wrong:**
The fitting code (staging/fitting.jl) fits `y(t) = A * exp(-gap * t) + C` to trace distance time series. The underlying assumption is that trace distance decays as a single exponential governed by the spectral gap. This is only true asymptotically (after initial transients die out). The actual trace distance dynamics is:

`d(rho(t), rho_gibbs) ~ sum_k |c_k| * exp(Re(lambda_k) * t)`

where lambda_k are Liouvillian eigenvalues and c_k are overlap coefficients (the code already analyzes these in gap_estimation.jl's `eigenbasis_overlap_analysis`).

**When the single-exponential model breaks:**
1. **Burn-in period:** Early transients from multiple eigenmode contributions make the initial data points follow a DIFFERENT exponential (or sum of exponentials). The `skip_initial` parameter handles this, but choosing it wrong biases the gap estimate.
2. **Nearly degenerate eigenvalues:** If lambda_2 and lambda_3 have similar real parts, the trace distance follows a SUM of two exponentials. A single-exponential fit will return an average of the two rates, which is wrong.
3. **Plateau misidentification:** The `C` parameter captures the asymptotic plateau. If the simulation has not run long enough, the fit interprets the still-decaying curve as a plateau plus fast exponential, overestimating the gap.
4. **Non-monotone trace distance:** For some initial states, trace distance can temporarily INCREASE before decaying (the "initial slip" phenomenon in open quantum systems). The fitting code does not handle this.

**Prevention:**
- The existing `skip_initial` parameter is critical. Default should be at least 0.1 (skip first 10% of data). For mixing time estimation where the goal is extrapolation, skip_initial = 0.2-0.3 is safer.
- Add a diagnostic: compute R-squared for the fit. The existing code does this (`_compute_r_squared`). Set a hard threshold: if R^2 < 0.9, flag the fit as unreliable and do not use it for mixing time extrapolation.
- For mixing time estimation specifically, consider fitting to the LOG of trace distance instead: `log(d(t)) ~ -gap * t + const`. This is a LINEAR fit (no iterative solver needed) and is more numerically stable. The existing `_log_linear_initial_guess` already does this for the initial guess -- consider making it an alternative fitting mode.
- Validate the gap estimate against the Krylov spectral gap when available (`run_krylov_spectrum` gives the exact gap for small systems). The existing `eigenbasis_overlap_analysis` provides this cross-validation.

**Detection:**
- Gap estimate varies significantly with `skip_initial` (sensitivity analysis).
- R-squared < 0.9 indicates model mismatch.
- Gap estimate is negative or much larger than expected from Krylov results.
- The existing `FitResult.converged` flag catches Levenberg-Marquardt non-convergence.

**Which phase:** Mixing time extrapolation phase. Must include validation against known spectral gaps.

---

### MOD-02: Extrapolation Beyond the Fitted Range Produces Meaningless Mixing Time

**What goes wrong:**
Mixing time estimation requires extrapolating the fitted exponential to find when `d(rho(t), rho_gibbs) < epsilon`. If the simulation only ran for time T_sim and the fitted curve predicts convergence at T_mix >> T_sim, the extrapolation is unreliable because:
1. The fit was trained on data in [0, T_sim] but is being evaluated at T_mix which may be 10-100x larger.
2. Small errors in the gap estimate are amplified exponentially: if gap_true = 0.1 and gap_fit = 0.11 (10% error), the mixing time estimate is off by `log(1/eps) * (1/0.1 - 1/0.11)` which for eps=1e-5 is `11.5 * 0.91 = 10.5` time units -- a significant error.
3. The single-exponential model may be valid in [0, T_sim] but the actual dynamics at T_mix could be dominated by a different eigenmode.

**Prevention:**
- Compute the extrapolation ratio: `T_mix / T_sim`. Flag results where this exceeds a threshold (suggest: 5x). At >10x, the extrapolation is essentially meaningless.
- Report confidence intervals on the mixing time using the gap uncertainty: `T_mix = log(A/epsilon) / gap`, so `dT_mix/dgap = -log(A/epsilon) / gap^2`. With `gap_se` from the fit, the mixing time uncertainty is approximately `T_mix_se = gap_se * |dT_mix/dgap|`.
- Require that the fitted trace distance at the END of the simulation is at least 10x smaller than at the START. If the curve has barely decayed, the gap estimate is poorly constrained.
- The existing `FitResult.gap_ci` (confidence interval) provides the raw statistical uncertainty. Propagate this through to mixing time.

**Detection:**
- Mixing time estimates that vary wildly with small changes to `skip_initial` or `total_time`.
- Mixing time confidence interval spanning more than an order of magnitude.
- Extrapolation ratio > 10.

**Which phase:** Mixing time extrapolation phase.

---

### MOD-03: save_every Off-by-One Corrupts Downstream Time Grid

**What goes wrong:**
The trajectory engine computes the number of saves as `num_saves = div(num_steps, save_every) + 1` (trajectories.jl line 604). The +1 is for the initial state (step 0). The time grid is `times[s] = (s - 1) * save_every * delta_step` (line 609). The DM thermalize code does NOT currently have save_every -- it records every step:
```julia
trace_distances = [trace_distance_h(Hermitian(evolving_dm), gibbs)]  # initial
for step in 1:num_steps
    ...
    push!(trace_distances, dist)  # after each step
end
time_steps = collect(0.0:config.delta:(num_steps * config.delta))
```

Adding save_every to the DM thermalize code creates several off-by-one risks:

1. **Mismatched lengths:** If `num_steps` is not divisible by `save_every`, the last partial interval may or may not be saved, creating a mismatch between `trace_distances` and `time_steps`.
2. **Initial state inclusion:** The trajectory code saves the initial measurement at step 0. If the DM code skips the initial step, the arrays are shifted by one.
3. **Time grid construction:** `times = collect(0.0 : save_every*delta : num_steps*delta)` vs `times = [(s-1) * save_every * delta for s in 1:num_saves]` produce different results when `num_steps % save_every != 0`.

**Downstream impact:**
The fitting code (`fit_exponential_decay`) takes `times` and `values` vectors that must have matching lengths. The gap estimation code (`estimate_spectral_gap`) constructs these from trajectory results. If DM thermalize produces arrays with different length semantics, the fitting will either error (length mismatch) or silently use misaligned time-value pairs.

**Prevention:**
- Copy the EXACT formula from trajectories.jl: `num_saves = div(num_steps, save_every) + 1`. Always include the initial state. Always use `times[s] = (s - 1) * save_every * delta`.
- Add an assertion: `@assert length(trace_distances) == length(time_steps) "save_every off-by-one: got $(length(trace_distances)) trace distances but $(length(time_steps)) time steps"`.
- Write a dedicated test: set `save_every = 7, num_steps = 20`. Expected `num_saves = div(20, 7) + 1 = 3` (saves at steps 0, 7, 14). The step 20 result is NOT saved because 20 % 7 != 0.
- Consider also saving the FINAL step regardless of save_every (save at steps 0, 7, 14, 20). This avoids losing the last ~6 steps of data. But document this clearly and ensure the time grid matches.

**Detection:**
- `ArgumentError: times and values must have the same length` from `fit_exponential_decay`.
- Fitting produces a gap that does not match the visual slope of the trace distance plot (misaligned time grid).
- Test that constructs trace_distances with save_every and verifies array sizes.

**Which phase:** save_every implementation phase. Must be tested before mixing time extrapolation depends on it.

---

### MOD-04: GC Pauses Disrupting Threaded Omega-Loop Performance

**What goes wrong:**
Julia's garbage collector (GC) is stop-the-world: when a GC event triggers, ALL threads are paused. In the trajectory engine, this is mitigated by the zero-allocation hot path (`step_along_trajectory!` allocates exactly 0 bytes per step, verified in test_allocation.jl). The omega-loop threading is different because:

1. The omega-loop body calls `mul!` which is allocation-free, but also does `@. scratch.jump_oft = ...` broadcasting which MAY allocate depending on the broadcast fusion.
2. The `_finalize_kraus_step!` at the end of each DM step calls `_build_cptp_channel` which calls `eigen(Hermitian(S))` -- this allocates (eigenvalue/eigenvector arrays). With precomputation, this allocation moves to construction time.
3. Even if the omega-loop body is allocation-free, other Julia code running concurrently (logging, GC of previous allocations) can trigger a GC pause.

**Impact on omega-loop threading:**
A GC pause during the omega-loop stalls ALL threads, including the ones doing useful frequency computation. If threads are unbalanced (some finish their chunk before others), idle threads could trigger GC from their own unrelated work, pausing the still-busy threads.

**Prevention:**
- Precompute the eigendecomposition (CRIT-02) to move the `eigen()` allocation out of the hot path.
- Verify that the omega-loop body is allocation-free using `@allocated` in a function barrier (same pattern as test_allocation.jl).
- Use `GC.enable(false)` around the threaded section if GC pauses are problematic, but ONLY if the allocation within the section is bounded and small. Re-enable immediately after. This is a last resort.
- Prefer `Threads.@spawn` with `@sync` over `Threads.@threads` -- the `@spawn` pattern gives better control over task granularity and avoids the scheduler overhead of `@threads :dynamic`.

**Detection:**
- Inconsistent timing for the same workload across runs (GC timing is non-deterministic).
- `@time` shows non-zero GC percentage.
- Threading speedup is below expected (1.5x with 4 threads instead of 3x+).

**Which phase:** Omega-loop threading phase. Profile before optimizing.

---

### MOD-05: Per-Jump Precomputation Changes Test Expectations for DM Thermalize

**What goes wrong:**
The existing `run_thermalize` recomputes the CPTP channel (R, K0, U_residual) from scratch at every step inside `_jump_contribution!` -> `_finalize_kraus_step!` -> `_build_cptp_channel`. Moving to per-jump precomputed channels changes the computation path:

**Before (current):**
```
For each DM step:
  pick random jump a
  accumulate R from all frequencies using current jump a
  hermitianize R
  S = f(R, delta)
  eigen(S) -> clamp -> sqrt -> U_residual  # recomputed
  apply K0 * rho * K0' + rho_jump + U_residual * rho * U_residual'
```

**After (precomputed):**
```
At construction time:
  For each jump a:
    accumulate R_a from all frequencies
    hermitianize R_a
    S_a = f(R_a, delta)
    eigen(S_a) -> clamp -> sqrt -> U_residual_a  # precomputed

For each DM step:
  pick random jump a
  accumulate rho_jump from all frequencies using jump a (needs evolving_dm)
  apply K0_a * rho * K0_a' + rho_jump + U_residual_a * rho * U_residual_a'
```

The numerical results will differ at machine precision because:
1. Eigendecomposition is numerically sensitive: two eigendecompositions of the "same" matrix constructed through different floating-point computation paths produce eigenvectors that differ by O(kappa * eps) where kappa is the condition number.
2. The R matrix in the current code accumulates within a scratch buffer that may have different rounding than a freshly-allocated matrix.

**Prevention:**
- Do NOT try to maintain bitwise identical results between precomputed and recomputed paths. Instead:
  - Verify trace distance convergence is equivalent (final trace distance within 2x of each other).
  - Verify the CPTP property is maintained (trace preservation, complete positivity).
  - Compare at a meaningful tolerance: `isapprox(result_precomp.final_dm, result_recomp.final_dm, atol=1e-8)` for a long run.
- Existing regression tests (test_regression.jl) store reference BSON outputs from the CURRENT (recomputed) code. These must be regenerated after switching to precomputed channels, or tests must be updated with looser tolerances.
- Keep the recomputed path available (behind a flag or separate method) for validation.

**Detection:**
- Regression tests failing with small numerical differences (1e-14 to 1e-10 range).
- CPTP test failures after precomputation.

**Which phase:** Per-jump precomputation phase.

---

### MOD-06: LsqFit SingularException from Flat or Noisy Trace Distance Data

**What goes wrong:**
The existing fitting code (staging/fitting.jl, lines 199-208) already handles `SingularException` from `stderror(fit)` and `confint(fit)`:
```julia
gap_se, gap_ci = try
    se = stderror(fit)
    ci = confint(fit; level=level)
    se[_IDX_GAP], (ci[_IDX_GAP][1], ci[_IDX_GAP][2])
catch e
    e isa LinearAlgebra.SingularException || rethrow(e)
    Inf, (-Inf, Inf)
end
```

This is good. But the SingularException occurs when the Jacobian is rank-deficient at the solution. For mixing time estimation, several common scenarios produce rank-deficient Jacobians:

1. **Already converged:** If `skip_initial` removes the transient and the remaining data is essentially flat at the plateau value, all three parameters (A, gap, C) are poorly determined. The Jacobian for A and gap becomes nearly zero.
2. **Too few data points after skip_initial:** The code requires >= 4 points (fitting.jl line 174), but with save_every the actual number of data points could be small.
3. **Very noisy data:** DM evolution trace distances can be noisy for small systems where the CPTP approximation error (O(delta^2) per step) accumulates.

**Downstream impact for mixing time estimation:**
If `gap_se = Inf` and `gap_ci = (-Inf, Inf)`, the mixing time uncertainty is also infinite, making the estimate useless. The code needs a graceful fallback.

**Prevention:**
- Before fitting, check data quality: if `max(values) - min(values) < 1e-6`, skip fitting and report "already converged."
- After fitting, check `FitResult.converged && FitResult.gap > 0 && FitResult.r_squared > 0.8`. The existing `_select_best_observable` (gap_estimation.jl line 31-59) already does this for trajectory-based gap estimation. The DM mixing time path should use the same quality gates.
- If gap_se = Inf, report the gap estimate as "unreliable" in the result struct rather than silently using it for mixing time computation.
- Consider implementing a fallback: if single-exponential fit fails, use the log-linear estimate from `_log_linear_initial_guess` directly as a rough gap. It is less accurate but does not require a converged Levenberg-Marquardt solve.

**Detection:**
- `FitResult.gap_se == Inf` or `FitResult.gap_ci == (-Inf, Inf)`.
- `FitResult.converged == false`.
- `FitResult.r_squared < 0` (model worse than horizontal line).

**Which phase:** Mixing time extrapolation phase.

---

## Minor Pitfalls

Issues that cause inconvenience or suboptimal performance but not incorrect results.

---

### MIN-01: Threading Threshold Too Low -- Overhead Dominates for Small Systems

**What goes wrong:**
For 2-3 qubit systems (dim = 4-8), the omega-loop has `2^num_energy_bits` frequency points (typically 16-64). Spawning Julia tasks, distributing work, and merging accumulators has overhead (microseconds). The per-frequency computation is a few small matrix multiplies (dim x dim, so 4x4 or 8x8) which take nanoseconds. Threading overhead exceeds computation time, making the threaded path slower than serial.

**Prevention:**
- Add a minimum work threshold: `use_threading = length(energy_labels) * dim^2 >= THREADING_THRESHOLD`. Suggested threshold: `dim^2 >= 32^2 = 1024` (i.e., dim >= 32, or 5+ qubits). Below this, run the omega-loop serially.
- Make the threshold configurable via a keyword argument with a sensible default.
- The trajectory engine does not have this problem because each trajectory takes many steps (thousands of matrix multiplies), so the per-trajectory work always exceeds threading overhead.

**Which phase:** Omega-loop threading phase.

---

### MIN-02: save_every Breaks Existing Plotting and Analysis Scripts

**What goes wrong:**
The existing `ThermalizeResults` stores `trace_distances::Vector{T}` and `time_steps::Vector{T}` with one entry per DM step. Downstream analysis code (plotting, convergence diagnostics) may assume `length(trace_distances) == num_steps + 1`. Introducing save_every changes this to `length(trace_distances) == div(num_steps, save_every) + 1`, breaking any code that indexes by step number.

**Prevention:**
- Default save_every to 1 for backward compatibility. Users opt-in to reduced saving.
- Store `save_every` in the results struct so downstream code can reconstruct the step-to-save mapping.
- Update the results BSON serialization (results.jl) to include the new field.
- Add the `save_every` field as a Union{Nothing, Int} defaulting to nothing (meaning 1) for backward compatibility with old BSON files.

**Which phase:** save_every implementation phase.

---

### MIN-03: Memory Spike from Per-Jump Precomputed Matrices at Large Scale

**What goes wrong:**
Per-jump precomputation stores 3 dim-by-dim matrices (K0, U_residual, R) per jump operator. The number of jump operators scales linearly with the number of qubits (typically `n_qubits` Pauli operators). Memory usage:

| n_qubits | dim | n_jumps | Memory per jump | Total |
|----------|-----|---------|-----------------|-------|
| 3 | 8 | 6 | 3 KB | 18 KB |
| 4 | 16 | 8 | 12 KB | 96 KB |
| 6 | 64 | 12 | 192 KB | 2.3 MB |
| 8 | 256 | 16 | 3 MB | 48 MB |
| 10 | 1024 | 20 | 48 MB | 960 MB |

At 10 qubits, storing per-jump precomputed data approaches 1 GB. This is borderline acceptable but surprising for users who do not expect the precomputation to use more memory than the density matrix itself (16 MB at dim=1024).

**Prevention:**
- Add a warning when total precomputed memory exceeds a threshold (100 MB): `@warn "Per-jump precomputation allocating $(total_mb) MB. Consider using lazy recomputation for systems this large."`.
- Optionally provide a `precompute_jumps::Bool` flag in the API. When false, fall back to per-step recomputation (current behavior).
- For large systems, consider precomputing only K0 and R, and computing U_residual lazily if the memory cost is dominated by the three matrices.

**Which phase:** Per-jump precomputation phase.

---

## Phase-Specific Warning Summary

| Phase | Pitfall | Severity | Key Mitigation |
|-------|---------|----------|----------------|
| Per-jump precomputation | CRIT-02: Stale eigendecomposition | Critical | Verify R is state-independent; assert S eigenvalues > -eps |
| Per-jump precomputation | MOD-05: Changed numerical results | Moderate | Regression test with loose tolerance; keep recomputed path for validation |
| Per-jump precomputation | MIN-03: Memory at large scale | Minor | Warning + optional lazy fallback |
| Omega-loop threading | CRIT-01: BLAS global thread state | Critical | try/finally save/restore pattern; never mix BLAS threading strategies |
| Omega-loop threading | CRIT-03: Race conditions in accumulators | Critical | Per-task accumulator pattern; do NOT use threadid() |
| Omega-loop threading | CRIT-04: False sharing | Critical | Independently allocated matrices, not 3D array views |
| Omega-loop threading | CRIT-05: BLAS oversubscription | Critical | set_num_threads(1) before threading; min_freqs threshold |
| Omega-loop threading | CRIT-06: Hermitianize placement | Critical | Hermitianize after merge, matching serial semantics |
| Omega-loop threading | MOD-04: GC pauses | Moderate | Zero-allocation hot path; precompute eigen |
| Omega-loop threading | MIN-01: Threading threshold | Minor | Skip threading for dim < 32 |
| save_every | MOD-03: Off-by-one errors | Moderate | Copy trajectory formula; assertion on array lengths |
| save_every | MIN-02: Breaking downstream code | Minor | Default save_every=1; store in results |
| Mixing time extrapolation | MOD-01: Model mismatch | Moderate | R^2 threshold; skip_initial; cross-validate with Krylov |
| Mixing time extrapolation | MOD-02: Extrapolation reliability | Moderate | Extrapolation ratio limit; propagate gap_se to T_mix |
| Mixing time extrapolation | MOD-06: Singular Jacobian | Moderate | Quality gates; fallback to log-linear estimate |

---

## Recommended Phase Ordering Based on Pitfall Dependencies

1. **Per-jump precomputation** first: This is the foundation. It eliminates the per-step `eigen()` call, which is needed before threading the omega-loop (otherwise each thread's eigen call would allocate and trigger GC). Also provides the biggest single-threaded speedup.

2. **save_every** second: Simple feature, independent of threading. Produces the time series data needed for mixing time extrapolation. Off-by-one bugs are easier to debug in serial code.

3. **Omega-loop threading** third: Builds on precomputed channels (no eigen in hot path = allocation-free threading body). BLAS thread management pattern is well-established from trajectory engine.

4. **Mixing time extrapolation** last: Depends on save_every for efficient data collection. Depends on understanding the DM convergence behavior, which benefits from having the threaded DM engine for faster iteration.

---

## Sources

- [Julia Multi-Threading Documentation](https://docs.julialang.org/en/v1/manual/multi-threading/) (HIGH confidence -- official docs)
- [Julia PSA: Don't Use threadid()](https://julialang.org/blog/2023/07/PSA-dont-use-threadid/) (HIGH confidence -- official Julia blog)
- [Julia for HPC: Multithreading](https://enccs.github.io/julia-for-hpc/multithreading/) (MEDIUM confidence -- ENCCS/community resource)
- [BLAS Thread Count vs Julia Thread Count](https://discourse.julialang.org/t/blas-thread-count-vs-julia-thread-count/57197) (MEDIUM confidence -- Discourse discussion)
- [Document the interaction of Julia and BLAS threads, Issue #44201](https://github.com/JuliaLang/julia/issues/44201) (MEDIUM confidence -- open issue with developer discussion)
- [MKL.jl Issue #106: Warn about NUM_THREADS behavior](https://github.com/JuliaLinearAlgebra/MKL.jl/issues/106) (MEDIUM confidence -- MKL-specific but pattern applies)
- [Threads.threadid() and Task Migration -- Sharp Edge](https://discourse.julialang.org/t/sharp-edge-with-threads-threadid-and-task-migration/124550) (MEDIUM confidence -- community report)
- [OhMyThreads.jl: Thread-Safe Storage](https://juliafolds2.github.io/OhMyThreads.jl/stable/literate/tls/tls/) (MEDIUM confidence -- alternative threading library docs)
- [False Sharing in Multi-threading](https://blog.jling.dev/blog/false_share/) (MEDIUM confidence -- Julia-specific blog post)
- [Exponential Data Fitting, SDSU Research Report](https://www.csrc.sdsu.edu/research_reports/CSRSR2009-04.pdf) (MEDIUM confidence -- academic reference)
- [LsqFit.jl Tutorial](https://julianlsolvers.github.io/LsqFit.jl/latest/tutorial/) (MEDIUM confidence -- official LsqFit docs)
- Direct codebase analysis of QuantumFurnace.jl v2.0 src/ and test/ directories (HIGH confidence)
