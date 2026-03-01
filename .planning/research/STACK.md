# Stack Research: v2.1 Speedup & Mixing Time Features

**Domain:** Julia threading patterns, BLAS thread control, exponential fitting for quantum simulation
**Researched:** 2026-03-01
**Confidence:** HIGH (patterns verified against codebase analysis, Julia official docs, community best practices)

---

## Scope

This STACK.md covers stack additions and patterns needed for the v2.1 milestone features:

1. **Per-jump precomputation** of K0^a, U_residual^a, U_coherent^a (already exists in trajectory path -- extend to DM thermalization)
2. **Multi-threaded omega-loops** with thread-local accumulators for Lindbladian/Krylov/DM paths
3. **Multi-threaded BLAS** for DM thermalization steps (different strategy from trajectory parallelism)
4. **save_every parameter** for trace distance in run_thermalize
5. **Mixing time estimation** via exponential fit on trace distance convergence curve

This does NOT re-research the existing dependency stack (KrylovKit, FINUFFT, Arpack, etc.).

---

## 1. Threading Patterns for omega-Loop Parallelization

### Context: The Problem

The frequency loop in `_jump_contribution!` (Thermalize path) and `apply_lindbladian!`/`_apply_lindbladian_impl!` (Krylov matvec) iterates over 50-100 energy labels, performing a GEMM sandwich (`L*rho*L'`) at each frequency. This loop is currently serial. For 100 energy labels with dim=64 matrices, this is 100 sequential GEMM calls on small (64x64) matrices.

The trajectory path already parallelizes at the trajectory level (via `Threads.@spawn` with `_copy_workspace_for_thread`), not at the omega-loop level, because each trajectory is independent. For the DM/Krylov paths, there is only one density matrix, so trajectory-level parallelism does not apply. The omega-loop is the natural parallelization target.

### Recommendation: `@sync`/`Threads.@spawn` with Pre-Allocated Thread-Local Accumulators

**Do NOT use `@threads :static`** for the omega-loop. Use the existing pattern from trajectories: pre-allocate per-task scratch buffers, spawn tasks via `@sync`/`Threads.@spawn`, and sum results.

**Why `@spawn` over `@threads`:**
- `@threads :static` pins tasks to threads and prevents composable nesting. If a caller already has threading (e.g., future multi-system sweep), `@threads` inside `@threads` deadlocks.
- `@threads :dynamic` (Julia >= 1.11 default) is better, but still schedules one iteration per task, incurring spawn overhead for 50-100 short iterations.
- `@spawn` with chunking matches the existing trajectory pattern exactly, avoids introducing a second threading idiom, and allows explicit control over chunk granularity.

**Confidence: HIGH** -- The codebase already uses `@sync`/`Threads.@spawn` with per-task workspace copies for trajectories (see `_run_batch_no_obs!` at trajectories.jl:494-509). This is the Julia-recommended approach since the [PSA: Don't use threadid()](https://julialang.org/blog/2023/07/PSA-dont-use-threadid/) blog post.

### Pattern: Thread-Local Matrix Accumulators for omega-Loop

```julia
function _jump_contribution_threaded!(
    L_target::AbstractMatrix{<:Complex},    # or evolving_dm for Thermalize
    jump::JumpOp,
    ham_or_trott,
    config::Config{S, D},
    precomputed_data,
    ws;                                     # workspace or scratch
    kwargs...
) where {S, D}

    energy_labels = precomputed_data.energy_labels
    n_omega = length(energy_labels)

    # Only parallelize if enough work to amortize spawn cost
    # Rule of thumb: each omega iteration does 2-3 GEMM on dim x dim matrices
    # For dim=64, one GEMM ~ 5-10us, so each iteration ~ 15-30us
    # Spawn overhead ~ 1-5us. With 50 labels, serial is fine for dim<32.
    # For dim>=64 and n_omega>=20, threading wins.
    nt = min(Threads.nthreads(), n_omega)

    if nt <= 1 || n_omega < 20
        # Serial fallback -- current code unchanged
        _jump_contribution_serial!(L_target, jump, ...)
        return L_target
    end

    # Pre-allocate per-task accumulators (same shape as L_target or R/rho_jump)
    chunks = _partition_trajectories(1:n_omega, nt)  # reuse existing helper
    CT = eltype(L_target)
    dim = size(L_target, 1)  # or isqrt(size(L_target, 1)) for Liouvillian

    # Per-task: own scratch buffers + own accumulator
    task_accumulators = [zeros(CT, size(L_target)) for _ in 1:length(chunks)]
    task_scratches = [_make_omega_scratch(CT, dim) for _ in 1:length(chunks)]

    old_blas = BLAS.get_num_threads()
    BLAS.set_num_threads(1)  # Single-threaded BLAS within each task
    try
        @sync for (idx, chunk) in enumerate(chunks)
            Threads.@spawn _omega_chunk!(
                task_accumulators[idx],
                task_scratches[idx],
                jump, energy_labels, chunk, ...
            )
        end
    finally
        BLAS.set_num_threads(old_blas)
    end

    # Reduce: sum all task accumulators into L_target
    for acc in task_accumulators
        L_target .+= acc
    end

    return L_target
end
```

### Per-Task Scratch Struct for omega-Loop

The omega-loop needs these scratch buffers per iteration:
- `jump_oft` (dim x dim): the OFT-filtered jump operator A(omega)
- `LdagL` (dim x dim): the L'*L product (for R accumulation in Thermalize, or sandwich in Liouvillian)
- `sandwich_tmp` (dim x dim): intermediate GEMM result

These are the same buffers currently in `ThermalizeScratch` and `LiouvillianScratch`. For threading, each task needs its own copy:

```julia
struct OmegaLoopScratch{T<:Complex}
    jump_oft::Matrix{T}
    LdagL::Matrix{T}
    sandwich_tmp::Matrix{T}
end

function OmegaLoopScratch(::Type{CT}, dim::Int) where {CT<:Complex}
    Zm() = zeros(CT, dim, dim)
    OmegaLoopScratch{CT}(Zm(), Zm(), Zm())
end
```

**Memory cost:** 3 * dim^2 * 16 bytes per task. For dim=64 and 8 threads: 3 * 4096 * 16 * 8 = 1.5 MB. Negligible.

### When NOT to Parallelize the omega-Loop

| dim | n_omega | GEMM time per omega | Total serial | Spawn overhead | Verdict |
|-----|---------|---------------------|-------------|----------------|---------|
| 16  | 50      | ~1 us               | ~50 us      | ~5 us * 8      | Serial  |
| 32  | 50      | ~3 us               | ~150 us     | ~40 us         | Maybe   |
| 64  | 50      | ~15 us              | ~750 us     | ~40 us         | Thread  |
| 128 | 100     | ~100 us             | ~10 ms      | ~40 us         | Thread  |

**Recommendation:** Add a `parallel_omega::Bool` keyword (default: `dim >= 32 && n_omega >= 20`) or compute it automatically. The trajectory path already has this pattern (`ntraj > 1 && Threads.nthreads() > 1`).

**Confidence: MEDIUM** -- The dim/n_omega thresholds are rough estimates based on BLAS benchmark data. Actual thresholds should be tuned with the project's specific workloads.

---

## 2. BLAS Thread Control Strategy

### Context: Two Different Parallelism Regimes

The codebase has two fundamentally different parallelism scenarios:

**Regime A: Trajectory Parallelism (existing)**
- Many independent trajectories, each doing small sequential BLAS calls (gemv on dim-vectors, rank-1 updates)
- Current strategy: `BLAS.set_num_threads(1)` + Julia threads for trajectories
- This is CORRECT. Each trajectory does O(dim) work per BLAS call. Multi-threaded BLAS on such small operations adds overhead.

**Regime B: DM Thermalization / Krylov Matvec (v2.1 target)**
- Single density matrix, iterated over omega-loop
- Each omega iteration does `mul!(C, A, B)` on dim x dim matrices (O(dim^3) BLAS gemm)
- Two sub-strategies depending on dim:
  - **dim <= 64:** GEMM is fast regardless. Parallelize the omega-loop with single-threaded BLAS per task.
  - **dim >= 128:** Single GEMM is expensive. Can use multi-threaded BLAS for each GEMM within a serial omega-loop.

### Recommended Strategy: Adaptive BLAS Threading

```julia
function _select_blas_strategy(dim::Int, n_omega::Int)
    n_julia_threads = Threads.nthreads()
    if n_julia_threads == 1
        return :serial     # No threading at all
    elseif dim <= 64 && n_omega >= 20
        return :omega_parallel   # Multi-task omega-loop, BLAS.set_num_threads(1)
    elseif dim >= 128
        return :blas_parallel    # Serial omega-loop, BLAS.set_num_threads(n_julia_threads)
    else
        return :omega_parallel   # Default to omega-parallelism
    end
end
```

**Regime B1: omega-parallel (dim <= 64)**
- Set `BLAS.set_num_threads(1)` (same as trajectory path)
- Spawn N tasks, each processes a chunk of energy labels
- Each task accumulates into its own R/rho_jump buffer
- Sum accumulators at the end

**Regime B2: blas-parallel (dim >= 128)**
- Keep `BLAS.set_num_threads(Threads.nthreads())` (or physical cores)
- Run the omega-loop serially
- Each `mul!(C, A, B)` call benefits from multi-threaded GEMM internally
- No need for per-task scratch -- single-threaded Julia, multi-threaded BLAS

**Crossover point:** For OpenBLAS, multi-threaded GEMM starts beating serial at approximately dim >= 128 on modern CPUs. For dim 64-128, the benefit is marginal. The omega-parallel strategy is almost always better for dim <= 64 because the per-iteration GEMM is too small for BLAS threading overhead.

**Confidence: HIGH** -- This matches the ITensors.jl recommendation: "If you are using ITensors.jl on a single node with multi-threaded BLAS, you will likely see the best performance by setting the number of Julia threads to 1 and using BLAS threads for the dense linear algebra operations." The key insight is that Julia threading and BLAS threading compete for the same CPU cores.

### Critical: BLAS Thread Scoping

`BLAS.set_num_threads()` is a GLOBAL setting -- it affects all threads in the process. This is why the existing trajectory code uses try/finally:

```julia
old_blas = BLAS.get_num_threads()
BLAS.set_num_threads(1)
try
    # ... parallel work ...
finally
    BLAS.set_num_threads(old_blas)
end
```

**This pattern is mandatory.** There is no per-task BLAS thread control in Julia 1.11/1.12. The `BLAS.set_num_threads()` call is process-global and affects the entire BLAS library (OpenBLAS).

**Confidence: HIGH** -- Verified via [Julia docs](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/) and [GitHub issue #44201](https://github.com/JuliaLang/julia/issues/44201). Per-task BLAS thread control is a desired but unimplemented feature.

### DM Thermalization: Sequential Steps, Parallel Within

For `run_thermalize`, the outer step loop is inherently sequential (each step depends on the previous density matrix). The parallelism opportunities are:

1. **Within each step:** Parallelize the omega-loop inside `_jump_contribution!` (as described above)
2. **The CPTP channel construction** (`_build_cptp_channel`): Contains an eigendecomposition (`eigen(Hermitian(S))`). Julia's `eigen` for `Hermitian` matrices already uses multi-threaded LAPACK internally when BLAS threads > 1. For the omega-parallel regime, this means we need to temporarily restore BLAS threads for the eigendecomposition step.

**Per-jump precomputation avoids this problem entirely:** If K0^a, U_residual^a, U_coherent^a are precomputed once before the step loop (as the milestone requires), then `_build_cptp_channel` is not called during the hot loop. The eigendecomposition happens only during setup, where BLAS threading is at its default value.

**Confidence: HIGH** -- The trajectory workspace already precomputes these via `_build_trajectory_workspace` (see trajectories.jl:104-122). Extending this to the Thermalize DM path is a straightforward code reuse.

---

## 3. Per-Jump Precomputation for DM Thermalization

### What Already Exists (Trajectory Path)

The trajectory workspace (`_build_trajectory_workspace`, trajectories.jl:61-152) already precomputes per-operator:
- `Rs[a]`: The accumulated R^a matrix (R = sum of L_k'*L_k weighted by transition rates)
- `K0s[a]`: K0 = I - alpha * R^a
- `U_residuals[a]`: sqrt_psd of residual matrix S = (2*alpha - delta)*R - alpha^2*R^2
- `U_Bs[a]`: Per-operator coherent unitary exp(-i*delta/p_jump * B^a)

This is computed via `_precompute_R` (one method per domain) + `_build_cptp_channel`.

### What v2.1 Needs

The DM Thermalize path (`run_thermalize`, furnace.jl:143-223) currently:
1. Calls `_precompute_data` once for transition weights (shared across all jumps)
2. Precomputes coherent unitaries once
3. In the step loop: selects a random jump, calls `_jump_contribution!` which recomputes R and rho_jump from scratch every step
4. Then calls `_finalize_kraus_step!` which calls `_build_cptp_channel` every step

The per-step `_build_cptp_channel` call includes an eigendecomposition -- this is the main performance bottleneck for the DM path.

**v2.1 precomputation:** Reuse the trajectory path's precomputation logic. Build per-jump K0^a, U_residual^a, U_coherent^a once before the step loop. Then the step loop becomes:
1. Select random jump a
2. Apply coherent unitary (U_B^a, already precomputed)
3. Apply `K0^a * rho * K0^a' + rho_jump_a + U_res^a * rho * U_res^a'`

Step 3 still requires computing `rho_jump_a` (the dissipative sandwich sum over omega), but the eigendecomposition is eliminated from the hot path.

### Fitting into the Workspace Struct

The `Workspace{Thermalize}` struct in the current codebase does not have the per-jump vector fields (`Rs`, `K0s`, `U_residuals`, `U_Bs`). These fields exist on `Workspace{Trajectory}`. Two options:

**Option A (recommended): Add fields to Workspace or use a separate precomputed struct.**
The Workspace struct already has these fields declared (see structs.jl:408-411), they are just set to `nothing` for non-Trajectory simulation modes. Set them for Thermalize too.

**Option B: Build a lightweight NamedTuple during setup.**
```julia
precomputed_kraus = (
    K0s = Vector{Matrix{CT}}(...),
    U_residuals = Vector{Matrix{CT}}(...),
    U_Bs = Vector{Union{Nothing, Matrix{CT}}}(...),
)
```

**Recommendation: Option A.** The fields already exist in the Workspace struct. Setting them for Thermalize is the minimal-diff approach.

**Confidence: HIGH** -- Direct code analysis. The Workspace struct already has `K0s`, `U_residuals`, `U_Bs` fields; they just need to be populated for the Thermalize path.

---

## 4. Dependencies: What to Add/Change

### New Dependencies Required

| Dependency | Version | Purpose | Justification |
|------------|---------|---------|---------------|
| LsqFit.jl | 0.15+ (latest) | Exponential fitting for mixing time estimation | Was previously a dependency (removed in v2.0 Phase 37 cleanup when fitting.jl was staged). Now needs to be re-added as fitting.jl moves from staging to src. |

### No Other New Dependencies Needed

| Considered | Verdict | Why Not |
|------------|---------|---------|
| OhMyThreads.jl | Not needed | Provides `tmapreduce` and `TaskLocalValue`, but the codebase already has a working `@sync`/`@spawn` pattern with manual task-local scratch. Adding OhMyThreads would introduce a new dependency and a different threading idiom for no clear benefit. The existing pattern is well-understood and proven in the trajectory path. |
| FLoops.jl / Transducers.jl | Not needed | Over-engineered for this use case. The omega-loop reduction is a simple sum of matrix accumulators, not a complex data pipeline. |
| ThreadPinning.jl | Not needed | Useful for benchmarking BLAS/Julia thread interaction, but not a runtime dependency. Can use for development benchmarking without adding to Project.toml. |
| Optim.jl | Not needed | The mixing time estimation uses LsqFit's Levenberg-Marquardt, not general optimization. Optim.jl was previously used only by log_sobolev.jl (also staged). |

### LsqFit.jl Re-Addition

LsqFit was removed from `[deps]` in Project.toml during Phase 37 (v2.0 cleanup) because fitting.jl was moved to staging. For v2.1:

1. Re-add to `[deps]`: `LsqFit = "2a205ff7-8bc8-5f38-aade-5765a3247307"` (UUID from Julia registry)
2. Add to `[compat]`: `LsqFit = "0.15, 0.16"` (current latest is 0.15.0 as of 2025)
3. Re-add `using LsqFit` to `src/QuantumFurnace.jl`
4. Move `src/staging/fitting.jl` back to `src/fitting.jl` and re-include

The existing `fitting.jl` code is feature-complete for the mixing time estimation use case. It provides:
- `fit_exponential_decay(times, values; skip_initial, p0, level)` -> `FitResult`
- Automatic log-linear initial guess
- Levenberg-Marquardt with gap > 0 constraint
- R-squared, confidence intervals, standard error on gap
- Robust handling of SingularException from degenerate fits

**Confidence: HIGH** -- Direct code analysis. fitting.jl has 217 LOC with comprehensive error handling, tested in Phase 21.

---

## 5. Mixing Time Estimation: Fitting Strategy

### Model

The trace distance convergence curve from `run_thermalize` follows:

```
d(t) = A * exp(-gap * t) + C
```

where:
- `gap` is the spectral gap (exponential convergence rate)
- `A` is the initial amplitude (depends on initial state)
- `C` is the asymptotic offset (should be ~0 for exact convergence, small positive for numerical residual)

The mixing time is estimated as:

```
t_mix(epsilon) = -ln(epsilon / A) / gap
```

for a target trace distance threshold `epsilon`.

### Fitting Approach: Use Existing `fit_exponential_decay`

The existing `fit_exponential_decay` in `staging/fitting.jl` is exactly what is needed. The mixing time estimation API wraps it:

```julia
function estimate_mixing_time(
    trace_distances::Vector{<:Real},
    time_steps::Vector{<:Real};
    epsilon::Float64 = 1e-3,
    skip_initial::Float64 = 0.0,
    extrapolate::Bool = true,
)
    fit = fit_exponential_decay(time_steps, trace_distances; skip_initial)

    if !fit.converged || fit.gap <= 0 || fit.r_squared < 0.5
        @warn "Poor fit quality" gap=fit.gap r_squared=fit.r_squared converged=fit.converged
    end

    # Compute mixing time from fit
    if extrapolate && fit.gap > 0 && fit.amplitude > 0
        t_mix = -log(epsilon / fit.amplitude) / fit.gap
    else
        # Fallback: find last time trace distance < epsilon, or NaN if never
        idx = findlast(d -> d <= epsilon, trace_distances)
        t_mix = idx === nothing ? NaN : time_steps[idx]
    end

    return (; mixing_time=t_mix, fit=fit, epsilon, extrapolate)
end
```

### Why NOT a Multi-Exponential Fit

For some systems, the convergence curve has multiple exponential modes:
```
d(t) = A1 * exp(-gap1 * t) + A2 * exp(-gap2 * t) + C
```

**Recommendation: Do NOT implement multi-exponential fitting for v2.1.** Reasons:
1. Multi-exponential fitting is notoriously ill-conditioned (the Prony problem). Parameters are highly correlated.
2. The mixing time is dominated by the slowest mode (smallest gap), which the single-exponential fit captures when `skip_initial > 0` removes the fast transients.
3. The existing `skip_initial` parameter (available in `fit_exponential_decay`) handles early-time fast modes: skip the first 20-30% of data, fit the tail where the slowest mode dominates.
4. Adding model complexity without clear need violates the principle of minimal changes.

**Confidence: HIGH** -- This is standard practice in spectral gap estimation. The existing code in `gap_estimation.jl` (the `estimate_spectral_gap` function) already uses single-exponential fitting with skip_initial on trajectory observable time series.

### `extrapolate` Flag Design

When `extrapolate=true`:
- Uses the fitted gap to compute `t_mix = -ln(epsilon/A) / gap`
- Allows estimating mixing time even when the simulation has not reached the target epsilon
- Reports quality metrics (R-squared, confidence interval on gap) for the user to assess reliability

When `extrapolate=false`:
- Only reports the mixing time if the trace distance actually dropped below epsilon during the simulation
- Returns NaN if not converged (conservative)
- Useful as a ground-truth check when the simulation runs long enough

**Confidence: HIGH** -- Direct design decision, well-scoped.

---

## 6. `save_every` for Trace Distance in `run_thermalize`

### Current Behavior

`run_thermalize` (furnace.jl:143-223) computes trace distance at EVERY step:
```julia
for step in 1:num_steps
    # ... apply jump ...
    dist = trace_distance_h(Hermitian(evolving_dm), gibbs)
    push!(trace_distances, dist)
end
```

For large systems (dim=128, num_steps=100000), `trace_distance_h` computes eigenvalues of a dim x dim Hermitian matrix at each step. This is O(dim^3) per step -- potentially dominating the runtime.

### Implementation

Add a `save_every::Int = 1` parameter:
```julia
function run_thermalize(
    jumps, config, hamiltonian, trotter;
    initial_dm=nothing, rng=Random.default_rng(),
    rescale_by_inv_prob=true,
    save_every::Int = 1,  # NEW
)
    # ...
    trace_distances = [trace_distance_h(Hermitian(evolving_dm), gibbs)]

    for step in 1:num_steps
        # ... apply jump ...
        if step % save_every == 0 || step == num_steps
            dist = trace_distance_h(Hermitian(evolving_dm), gibbs)
            push!(trace_distances, dist)
        end
    end

    # Adjust time_steps to match saved points
    saved_steps = [0; [s for s in 1:num_steps if s % save_every == 0 || s == num_steps]]
    time_steps = saved_steps .* config.delta
    # ...
end
```

This is identical to the `save_every` pattern already used in `run_trajectories` (trajectories.jl:602-610) and `run_observable_trajectories` (trajectories.jl:706-713).

**No new dependencies or patterns needed.**

**Confidence: HIGH** -- Direct extension of existing `save_every` pattern from trajectory code.

---

## 7. Version Compatibility

### Julia Version Requirements

| Requirement | Min Version | Notes |
|-------------|-------------|-------|
| Julia | 1.11 | Already set in Project.toml `[compat]`. Required for `using Base.Threads`, current `@spawn` semantics, `@threads :dynamic` default. |
| Julia (recommended) | 1.12 | Improved BLAS thread affinity (BLAS threads respect CPU affinity settings), better task scheduler. Not required but beneficial for BLAS thread strategy. |

### Dependency Compatibility

| Package | Version | Compatible With | Notes |
|---------|---------|-----------------|-------|
| LsqFit | 0.15+ | Julia 1.11-1.12 | Last breaking change was 0.14->0.15 (API change to `coef`, `stderror`, `confint` functions). The staging code already uses the 0.15 API. |
| LinearAlgebra | 1.11, 1.12 | `BLAS.set_num_threads` and `BLAS.get_num_threads` | Stable API. No breaking changes between 1.11 and 1.12 for these functions. |

---

## 8. Alternatives Considered

| Category | Recommended | Alternative | Why Not |
|----------|-------------|-------------|---------|
| omega-loop threading | `@sync`/`@spawn` + per-task accumulators | `@threads :dynamic` | `@threads` spawns one task per iteration, 50-100 tasks for omega-loop is fine but mismatches existing codebase pattern. `@spawn` with chunking is already the codebase idiom and allows explicit control over task count. |
| omega-loop threading | `@sync`/`@spawn` + manual scratch | OhMyThreads.jl `tmapreduce` | Adds external dependency. `tmapreduce` is elegant but the existing manual pattern is well-tested and understood. No clear benefit for this use case. |
| BLAS thread control | Adaptive (omega-parallel vs blas-parallel) | Always `BLAS.set_num_threads(1)` | For large dim (>=128), serial omega-loop with multi-threaded BLAS inside each GEMM is faster than N-way omega-parallel with single-threaded BLAS, because GEMM parallelizes better at large N than task spawning. |
| Fitting library | LsqFit.jl (Levenberg-Marquardt) | Optim.jl with custom loss | LsqFit is purpose-built for curve fitting with confidence intervals and residual analysis. Optim.jl would require reimplementing covariance estimation, confidence intervals, etc. LsqFit was already used and tested in Phase 21. |
| Fitting model | Single exponential A*exp(-gap*t)+C | Multi-exponential sum | Multi-exponential is ill-conditioned (Prony problem). `skip_initial` handles fast transients. Single-exponential with tail fitting is the standard approach for spectral gap estimation. |
| Thread-local storage | Pre-allocated Vector of scratch structs | `TaskLocalValue` (from OhMyThreads) | `TaskLocalValue` is elegant but requires OhMyThreads dependency. Pre-allocated scratch vectors (indexed by task chunk) are the existing pattern and avoid any external dependency. |
| Precomputation target | Populate existing Workspace fields | New struct for precomputed Kraus data | Workspace already has `K0s`, `U_residuals`, `U_Bs` fields (currently `nothing` for Thermalize). Populating them is zero-struct-change. |

---

## 9. What NOT to Use

| Avoid | Why | Use Instead |
|-------|-----|-------------|
| `threadid()` for indexing scratch buffers | Task migration between threads makes threadid()-indexed arrays unsafe. This is a [known Julia pitfall](https://julialang.org/blog/2023/07/PSA-dont-use-threadid/). | Pre-allocate per-chunk scratch arrays indexed by chunk number (1:n_chunks), not thread ID. |
| `@threads :static` in library code | Prevents composable nesting. If any caller also uses `@threads`, the inner loop runs on a single thread. | `@sync`/`@spawn` for composable parallelism. |
| `Threads.Atomic` for matrix accumulation | Atomics work for scalars, not for matrix element-wise accumulation. Would need one atomic per matrix element (dim^2 atomics), massive overhead. | Per-task accumulators with post-loop reduction (sum). |
| Multi-threaded BLAS + Julia threading simultaneously | BLAS threads and Julia threads compete for the same CPU cores. Running 8 Julia tasks each calling 8-thread BLAS = 64 threads contending for 8 cores = massive slowdown. | Choose one: either omega-parallel with `BLAS.set_num_threads(1)`, or serial omega with multi-threaded BLAS. Never both. |
| `SharedArrays` or `Distributed` | Process-level parallelism. Overkill for within-process matrix operations. Adds IPC overhead. | `Threads.@spawn` for shared-memory threading. |
| `@generated` functions for branching on dim thresholds | The serial/parallel decision is a runtime choice based on problem size, not a compile-time type specialization. `@generated` is inappropriate here. | Runtime `if` branch on `dim` and `n_omega`. |

---

## 10. Stack Patterns by Variant

### If dim <= 32 (small systems):

- Use serial omega-loop everywhere (trajectory, DM, Krylov)
- Use `BLAS.set_num_threads(1)` during trajectory parallelism
- `save_every` can be 1 (trace distance is cheap for small dim)
- Threading overhead dominates for small matrices; do not parallelize omega-loop

### If 32 < dim <= 64 (medium systems):

- Parallelize omega-loop with `@spawn` + `BLAS.set_num_threads(1)` per task
- Use `save_every >= 10` for trace distance in `run_thermalize`
- Per-jump precomputation saves one eigendecomposition per step (significant at dim=64: eigen is O(dim^3) ~ 260K ops)

### If dim >= 128 (large systems):

- Use serial omega-loop with `BLAS.set_num_threads(Threads.nthreads())` for multi-threaded GEMM
- `save_every >= 100` for trace distance (eigen of 128x128 Hermitian ~ 0.1ms per call)
- Per-jump precomputation is critical: eigendecomposition per step would be the dominant cost
- Consider whether trajectory parallelism (existing) or DM thermalization is the primary use case

### If system has many jump operators (>10):

- Per-jump precomputation cost scales linearly with n_jumps
- Memory for K0s, U_residuals: 2 * n_jumps * dim^2 * 16 bytes
- For n_jumps=20, dim=64: 2 * 20 * 4096 * 16 = 2.5 MB -- still negligible
- The step loop randomly selects one jump per step, so precomputation is amortized well

---

## 11. Installation / Dependency Changes

```toml
# Project.toml additions:

[deps]
# Re-add (was removed in v2.0 Phase 37):
LsqFit = "2a205ff7-8bc8-5f38-aade-5765a3247307"

[compat]
# Add:
LsqFit = "0.15, 0.16"
```

```julia
# src/QuantumFurnace.jl additions:
using LsqFit

# Move from staging to active:
include("fitting.jl")

# New exports:
export fit_exponential_decay, FitResult
export estimate_mixing_time  # new v2.1 function
```

---

## Sources

- [Julia Multi-Threading Documentation](https://docs.julialang.org/en/v1/base/multi-threading/) -- `@spawn`, `@sync`, `@threads` semantics. HIGH confidence.
- [PSA: Don't use threadid()](https://julialang.org/blog/2023/07/PSA-dont-use-threadid/) -- Why pre-allocated task-indexed buffers over threadid(). HIGH confidence.
- [Julia 1.12 Highlights](https://julialang.org/blog/2025/10/julia-1.12-highlights/) -- BLAS thread affinity improvements. HIGH confidence.
- [Multithreading -- Julia for HPC (ENCCS)](https://enccs.github.io/julia-for-hpc/multithreading/) -- Nested parallelism, BLAS thread interaction. HIGH confidence.
- [ITensors.jl Multithreading](https://itensor.github.io/ITensors.jl/dev/Multithreading.html) -- BLAS vs Julia thread strategy for dense linear algebra. HIGH confidence.
- [OhMyThreads.jl Documentation](https://juliafolds2.github.io/OhMyThreads.jl/stable/) -- TaskLocalValue, tmapreduce patterns. HIGH confidence (evaluated and decided against for this project).
- [BLAS Thread Count Discussion (Julia Discourse)](https://discourse.julialang.org/t/ideal-number-of-blas-threads/79197) -- Physical vs logical cores for BLAS. MEDIUM confidence.
- [Julia GitHub Issue #44201: Document BLAS/Julia thread interaction](https://github.com/JuliaLang/julia/issues/44201) -- Confirms BLAS threads are process-global, no per-task control. HIGH confidence.
- [Julia GitHub Issue #43292: BLAS call overhead with multi-threading](https://github.com/JuliaLang/julia/issues/43292) -- BLAS overhead for small matrices in threaded code. MEDIUM confidence.
- [LsqFit.jl Documentation](https://julianlsolvers.github.io/LsqFit.jl/latest/) -- curve_fit API, Levenberg-Marquardt, confidence intervals. HIGH confidence.
- Codebase analysis: `trajectories.jl`, `jump_workers.jl`, `furnace.jl`, `furnace_utensils.jl`, `staging/fitting.jl`, `staging/gap_estimation.jl`, `structs.jl`, `Project.toml` -- Direct code inspection. HIGH confidence.

---

*Stack research for: QuantumFurnace.jl v2.1 Speedup & Mixing Time*
*Researched: 2026-03-01*
