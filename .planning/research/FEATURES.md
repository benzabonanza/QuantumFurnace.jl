# Feature Landscape: v2.1 Speedup & Mixing Time

**Domain:** Performance optimization and mixing time estimation for quantum Gibbs sampling simulator
**Researched:** 2026-03-01
**Confidence:** HIGH (based on exhaustive reading of all source files, staging code, supplementary papers, and error catalogue)

---

## Scope

This document defines the concrete features for the v2.1 milestone, which has two pillars:

1. **Performance:** Per-jump CPTP channel precomputation, multi-threaded frequency summation, BLAS thread management, and `save_every` for trace distance computation in `run_thermalize`.
2. **Mixing Time Estimation:** Fit the trace distance convergence curve `d(t) ~ A*exp(-gap*t) + C` to estimate mixing time, with optional extrapolation to a target epsilon without running full convergence.

All recommendations are grounded in the existing codebase (`src/furnace.jl`, `src/jump_workers.jl`, `src/trajectories.jl`, `src/staging/fitting.jl`), the theoretical papers (Chen et al. 2023 Theorem III.1, Proposition II.3; Ramkumar & Soleimanifar 2024), and the comprehensive error catalogue (`supplementary-informations/error_catalogue_spectral_gap_estimation.md`).

---

## Current State: What Exists and What's Missing

### The Hot Loop Problem in `run_thermalize`

The current `_jump_contribution!` for DM thermalization (jump_workers.jl, lines 194-440) recomputes everything from scratch **every single delta-step**:

1. Iterates over all Bohr frequencies (or energy labels) for the selected jump
2. Computes OFT: `scratch.jump_oft = eigenbasis .* gaussian_filter`
3. Accumulates `scratch.R += rate^2 * L'L` (sum over all frequencies)
4. Accumulates `scratch.rho_jump += delta * rate^2 * L * rho * L'` (sum over all frequencies)
5. Calls `_finalize_kraus_step!` which does eigendecomposition of S for PSD guard

Steps 1-3 are **independent of the current density matrix** -- they depend only on the jump operator and static parameters. Yet they are recomputed at every step. The trajectory simulator already solved this in `_build_trajectory_workspace` (trajectories.jl, lines 61-152) where `Rs`, `K0s`, and `U_residuals` are precomputed once per jump and reused.

### The `save_every` Gap

The trajectory simulator already has `save_every` (trajectories.jl, line 576+) controlling how often observable measurements are saved. But `run_thermalize` (furnace.jl, lines 143-223) computes trace distance at **every single step** via:

```julia
dist = trace_distance_h(Hermitian(evolving_dm), gibbs)
push!(trace_distances, dist)
@printf("Dist to Gibbs: %s\n", dist)
```

The `trace_distance_h` call involves eigendecomposition (O(d^3)) which is expensive relative to the CPTP channel application for large d.

### Existing Fitting Infrastructure

`src/staging/fitting.jl` already contains:
- `fit_exponential_decay(times, values; skip_initial, p0, level)` -- single-exponential A*exp(-gap*t)+C via LsqFit.jl
- `_log_linear_initial_guess` -- auto-initialization from data
- `FitResult` struct with gap, amplitude, offset, CIs, R-squared, convergence status

This is currently in `src/staging/` (dormant), not included in the active module. It was built for trajectory-based spectral gap estimation (fitting observable decay curves). The mixing time estimation feature needs to **reactivate and adapt** this for DM trace distance curves.

### What the Papers Say About Mixing Time

**Chen et al. 2023, Proposition II.3:** If the anti-Hermitian defect ratio is small (lambda_1(H)/lambda_gap(H) <= 1/100), then:

```
t_mix(L) <= 3 * ln(3 * ||rho^{-1/2}||) / lambda_gap(H)
```

This gives the theoretical mixing time bound in terms of the spectral gap and the Gibbs state condition number.

**Operational definition:** The mixing time to precision epsilon is `t_mix(epsilon) = inf{t : ||e^{tL}[rho_0] - rho_beta||_1 <= epsilon}`. For exponential convergence, `d(t) = d(0) * exp(-lambda_gap * t)`, so:

```
t_mix(epsilon) = (1/lambda_gap) * ln(d(0)/epsilon)
```

This is exactly what fitting gives us: extract `lambda_gap` from the trace distance curve, then extrapolate.

**Ramkumar & Soleimanifar 2024:** The mixing time bound is:

```
t_mix(L_beta) <= O(beta*||H|| + log(n)) / lambda_gap(L_beta)
```

The spectral gap `lambda_gap` depends on the jump operators, energy resolution, and temperature. For graph-local jumps on cyclic graphs, `lambda_gap = Theta(n^{-3})` (Theorem 2.1). For unitary design jumps on bounded-degree systems, `lambda_gap = Omega(1)` (Theorem 2.2).

**Key insight for feature design:** The trace distance `d(t) = ||rho(t) - rho_beta||_1` does NOT perfectly follow a single exponential in practice because:

1. **Multiple modes:** The exact decay is a sum of exponentials from all Lindbladian eigenvalues (error catalogue, Error 6)
2. **Burn-in transient:** Early steps may not follow exponential decay if the initial state has strong overlap with fast-decaying modes
3. **Delta-step Trotterization bias:** The effective decay rate has O(delta) systematic bias (error catalogue, Error 1)
4. **DM vs trajectory:** For DM evolution, there is no trajectory noise -- the trace distance is deterministic. This eliminates Errors 2, 3 from the error catalogue (those affect trajectory-based estimation only)

The DM trace distance curve is a **clean, deterministic signal** -- much better suited for exponential fitting than noisy trajectory-averaged observables. This is a key advantage of doing mixing time estimation from `run_thermalize` rather than from trajectories.

---

## Table Stakes

Features users expect. Missing = milestone feels incomplete.

### TS-01: Per-Jump CPTP Channel Precomputation

| Aspect | Detail |
|--------|--------|
| **Why Expected** | Trajectory simulator already does this. DM simulator recomputing R, K0, U_residual every step is an obvious inefficiency. |
| **Complexity** | MEDIUM |
| **Notes** | Most impactful single optimization for run_thermalize |

**What to build:**

Precompute `R^a`, `K0^a`, `U_residual^a`, and `U_B^a` (coherent unitary) once per jump operator before the main loop, following the exact pattern from `_build_trajectory_workspace` (trajectories.jl, lines 104-122):

```julia
# Precompute per-operator CPTP channel parts
Rs = Vector{Matrix{CT}}(undef, n_jumps)
K0s = Vector{Matrix{CT}}(undef, n_jumps)
U_residuals = Vector{Matrix{CT}}(undef, n_jumps)

for a in 1:n_jumps
    _precompute_R([jumps[a]], ham_or_trott, config, precomputed_data, scratch)
    R_a = copy(scratch.R)
    R_a .*= (1.0 / p_jump)  # rescale_by_inv_prob

    (; K0, U_residual) = _build_cptp_channel(R_a, config.delta)

    Rs[a] = R_a
    K0s[a] = K0
    U_residuals[a] = U_residual
end
```

**What changes in the hot loop:**

The current `_jump_contribution!` for `Config{Thermalize}` does three things:
1. Apply coherent unitary (already cached as `coherent_unitary_cache`)
2. Compute R and rho_jump via frequency summation (the expensive part)
3. Call `_finalize_kraus_step!` which builds K0/U_residual from R and applies the channel

With precomputation, step 2 reduces to looking up `Rs[idx]` and step 3 uses precomputed `K0s[idx]` and `U_residuals[idx]`. The **only per-step work** that depends on rho is:
- `rho_jump = delta * sum_w rate^2(w) * L_w * rho * L_w'` (the dissipative sandwich)
- `rho_next = K0 * rho * K0' + rho_jump + U_res * rho * U_res'` (the CPTP application)

The R and eigendecomposition (most expensive parts) become one-time costs.

**Critical subtlety -- rho_jump still needs frequency loop:**

The rho_jump term `sum_w L_w * rho * L_w'` depends on rho, so it cannot be fully precomputed. However, R can be precomputed since `R = sum_w L_w' * L_w` is rho-independent. So we need:

Option A: Keep the frequency loop for rho_jump, skip R accumulation, use precomputed K0/U_residual.
Option B: Precompute the individual `L_w` operators and store them, then the frequency loop becomes BLAS sandwiches only.

Option A is simpler and already captures the biggest win (eliminating the per-step eigendecomposition in `_finalize_kraus_step!`). Option B requires storing O(n_freq) matrices per jump, which may be memory-prohibitive.

**Recommendation: Option A.** Keep frequency loop for rho_jump only; precompute R, K0, U_residual. The eigendecomposition in `_build_cptp_channel` (the PSD square root of S) is the dominant cost being eliminated.

**Dependencies on existing code:**
- `_precompute_R` already exists for Energy and Time/Trotter domains (jump_workers.jl, lines 193-300+)
- `_build_cptp_channel` already exists (furnace_utensils.jl, lines 183-200)
- Coherent unitary precomputation already exists via `_precompute_coherent_unitary` (furnace.jl, lines 180-181)
- Pattern already proven in trajectory workspace construction

**Edge case:** BohrDomain has a different R accumulation pattern using `alpha(bohr_freqs, nu_2)` Kossakowski matrices. The precomputation needs separate handling for BohrDomain but follows the same principle.

---

### TS-02: Precomputed CPTP Application (Eliminate Per-Step Eigendecomposition)

| Aspect | Detail |
|--------|--------|
| **Why Expected** | `_finalize_kraus_step!` calls `_build_cptp_channel` which does `eigen(Hermitian(S))` every step. This O(d^3) eigendecomposition is the single most expensive per-step operation for d >= 32. |
| **Complexity** | LOW (follows directly from TS-01) |
| **Notes** | Once R is precomputed, K0 and U_residual are constants |

**What changes:**

Replace `_finalize_kraus_step!(evolving_dm, config.delta, scratch)` with a simpler function that uses precomputed K0 and U_residual:

```julia
function _apply_precomputed_channel!(
    evolving_dm::Matrix{<:Complex},
    K0::Matrix{<:Complex},
    U_residual::Matrix{<:Complex},
    rho_jump::Matrix{<:Complex},
    scratch::ThermalizeScratch{<:Complex},
)
    # rho_next = K0 * rho * K0' + rho_jump + U_res * rho * U_res'
    mul!(scratch.sandwich_tmp, K0, evolving_dm)
    mul!(scratch.rho_next, scratch.sandwich_tmp, K0')
    scratch.rho_next .+= rho_jump

    mul!(scratch.sandwich_tmp, U_residual, evolving_dm)
    mul!(scratch.rho_next, scratch.sandwich_tmp, U_residual', 1.0, 1.0)

    hermitianize!(scratch.rho_next)
    copyto!(evolving_dm, scratch.rho_next)
    return evolving_dm
end
```

This is 4 BLAS gemm calls (O(d^3) each) vs the current path which does those 4 calls PLUS an eigendecomposition. Net savings: one full eigendecomposition per step.

---

### TS-03: save_every for Trace Distance in run_thermalize

| Aspect | Detail |
|--------|--------|
| **Why Expected** | Trajectory simulator already has this pattern. DM simulator computing trace_distance_h every step is wasteful for long runs. |
| **Complexity** | LOW |
| **Notes** | Direct pattern reuse from trajectory simulator |

**What to build:**

Add a `save_every::Int=1` parameter to `run_thermalize`. The trace distance computation and printf only execute when `step % save_every == 0`:

```julia
function run_thermalize(
    jumps, config, hamiltonian, trotter=nothing;
    initial_dm=nothing, rng=Random.default_rng(),
    rescale_by_inv_prob=true,
    save_every::Int=1,              # NEW
)
    # ... setup ...
    num_steps = Int(ceil(config.mixing_time / config.delta))

    # Pre-allocate with correct size
    num_saves = div(num_steps, save_every) + 1
    trace_distances = Vector{T}(undef, num_saves)
    time_steps = Vector{T}(undef, num_saves)

    save_idx = 1
    trace_distances[1] = trace_distance_h(Hermitian(evolving_dm), gibbs)
    time_steps[1] = 0.0
    save_idx = 2

    for step in 1:num_steps
        # ... apply channel ...

        if step % save_every == 0
            dist = trace_distance_h(Hermitian(evolving_dm), gibbs)
            trace_distances[save_idx] = dist
            time_steps[save_idx] = step * config.delta
            save_idx += 1
            @printf("Step %d/%d, dist to Gibbs: %.6e\n", step, num_steps, dist)
            if dist < convergence_cutoff
                resize!(trace_distances, save_idx - 1)
                resize!(time_steps, save_idx - 1)
                break
            end
        end
    end

    # Trim to actual size
    resize!(trace_distances, save_idx - 1)
    resize!(time_steps, save_idx - 1)
    # ...
end
```

**Time steps alignment:** The returned `time_steps` vector should correspond to the saved trace distances. With `save_every=k`, time_steps are `[0, k*delta, 2k*delta, ...]`. This matches the trajectory simulator convention.

**Impact:** For n=6 (dim=64), `trace_distance_h` involves an eigendecomposition of a 64x64 matrix, roughly 0.1ms. With save_every=10, this reduces from (say) 10000 calls to 1000 calls, saving ~1 second. For n=8 (dim=256), the savings become significant.

---

### TS-04: Mixing Time Estimation via Exponential Fit

| Aspect | Detail |
|--------|--------|
| **Why Expected** | This is the primary scientific deliverable of the milestone. Without it, the speedup features are optimization-only with no new capability. |
| **Complexity** | MEDIUM |
| **Notes** | Reactivates and adapts staging/fitting.jl |

**What to build:**

A function that takes the output of `run_thermalize` (trace distances and time steps) and fits `d(t) = A*exp(-gap*t) + C` to extract the effective spectral gap, then optionally extrapolates to a target epsilon.

```julia
struct MixingTimeEstimate
    gap::Float64                          # fitted effective spectral gap
    amplitude::Float64                    # A in A*exp(-gap*t)+C
    offset::Float64                       # C (should be ~0 for trace distance to fixed point)
    gap_ci::Tuple{Float64, Float64}       # confidence interval on gap
    gap_se::Float64                       # standard error on gap
    r_squared::Float64                    # goodness of fit
    converged::Bool                       # LM convergence
    t_mix_fitted::Float64                 # mixing time to offset C from fitted data
    t_mix_extrapolated::Union{Nothing, Float64}  # extrapolated t_mix to target epsilon
    target_epsilon::Union{Nothing, Float64}       # the target epsilon
    fit_window::Tuple{Float64, Float64}           # (t_min, t_max) used for fitting
end
```

**The fitting function:**

```julia
function estimate_mixing_time(
    result::ThermalizeResults;
    skip_initial::Float64 = 0.2,    # skip first 20% as burn-in
    target_epsilon::Union{Nothing, Float64} = nothing,
    extrapolate::Bool = false,
) -> MixingTimeEstimate
```

**Burn-in detection (critical design decision):**

The trace distance curve `d(t) = ||rho(t) - rho_beta||_1` during DM thermalization has three regimes:

1. **Transient regime (burn-in):** Early steps where multiple Lindbladian modes contribute. The trace distance may even increase briefly if the initial state has poor overlap with the slow mode. The curve is not well-described by a single exponential.

2. **Exponential regime:** After the fast modes have decayed, the trace distance follows `d(t) ~ A*exp(-lambda_gap*t) + C` where lambda_gap is the spectral gap.

3. **Plateau regime:** The trace distance saturates at a floor set by the fixed-point accuracy (the deviation of the CPTP channel fixed point from the true Gibbs state). For KMS-BohrDomain this is ~1e-15; for TrotterDomain this is ~1e-8.

**Approaches evaluated for burn-in detection:**

| Approach | Pros | Cons | Recommendation |
|----------|------|------|----------------|
| Fixed fraction skip (skip_initial=0.2) | Simple, deterministic, no false positives | May skip too much or too little; not adaptive | Use as **default** |
| Effective rate stabilization | Theoretically optimal; identifies the exact moment single-exp kicks in | Requires computing log-derivative which amplifies noise; for DM trace distance (noiseless) this is clean | Use as **validation** diagnostic, not as primary method |
| Two-exponential fit (A1*exp(-g1*t) + A2*exp(-g2*t) + C) | Absorbs burn-in contamination; g1 gives true gap | Sensitive to initialization; 5-parameter fit needs more data | **Defer** to future milestone (already deferred from v1.4) |
| Piecewise: fit only tail where d(t) < 0.1*d(0) | Robustly avoids transient; simple threshold | May lose too much data for short runs | Use as **fallback** |

**Recommendation:** Use `skip_initial=0.2` as default. Add an optional `detect_burnin=false` flag that, when true, uses the effective rate method to automatically determine `t_min`. The effective rate method works well for noiseless DM data:

```julia
function _detect_burnin(times, trace_dists; lag=5, stability_threshold=0.1)
    # Compute effective rate: lambda_eff(t) = -log(d(t+tau)/d(t)) / tau
    # Find first t where lambda_eff stabilizes (relative variation < threshold
    # over a window)
    # Return that t as the burn-in cutoff
end
```

**Extrapolation:**

When `extrapolate=true` and `target_epsilon` is provided:

```
t_mix(epsilon) = -ln(epsilon / A) / gap
```

where A and gap come from the fit. This is valid when:
- The fit has R^2 > 0.95 (good single-exponential fit)
- The offset C is small relative to epsilon (C << epsilon)
- The fitted gap is physically reasonable (gap > 0, gap_se / gap < 0.5)

If any validity check fails, return `t_mix_extrapolated = nothing` with a warning.

**Relationship to existing staging/fitting.jl:**

The existing `fit_exponential_decay` does exactly what's needed. The mixing time estimation function wraps it with:
1. Pre-processing: skip burn-in, select fit window, validate data quality
2. Post-processing: compute mixing time from fitted gap, extrapolate if requested
3. Quality gates: R^2 threshold, CI width check, offset sanity check

The `FitResult` struct from staging/fitting.jl should be promoted from staging to active code and embedded within `MixingTimeEstimate`, or its fields replicated.

---

### TS-05: Multi-threaded BLAS for DM Thermalization

| Aspect | Detail |
|--------|--------|
| **Why Expected** | The trajectory simulator already manages BLAS threads (BLAS.set_num_threads(1) during threaded execution). DM thermalization is single-threaded on Julia side, so it should use ALL BLAS threads for the dense matrix operations. |
| **Complexity** | LOW |
| **Notes** | May already be the default behavior, but explicit management is good practice |

**What to build:**

Ensure that during `run_thermalize`, BLAS has access to all available threads. The trajectory simulator explicitly sets `BLAS.set_num_threads(1)` to avoid oversubscription. The DM simulator should do the opposite:

```julia
function run_thermalize(...)
    # Ensure BLAS uses all available threads for dense operations
    blas_threads = BLAS.get_num_threads()
    # No change needed if already at max; but if a prior trajectory run
    # set it to 1, restore it
    BLAS.set_num_threads(Sys.CPU_THREADS)
    try
        # ... main loop ...
    finally
        BLAS.set_num_threads(blas_threads)
    end
end
```

The main beneficiaries are the `mul!` calls in `_apply_precomputed_channel!` and the rho_jump sandwich. For dim=256 (n=8), multi-threaded BLAS can give 4-8x speedup on the matrix multiplications.

**Caution:** Do NOT set BLAS threads inside the hot loop. Set once at entry, restore at exit.

---

## Differentiators

Features that set the product apart. Not expected, but valued.

### DIFF-01: Multi-Threaded Frequency Summation

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| Thread-parallel w-loop in rho_jump accumulation | For large frequency grids (num_energy_bits >= 8), the w-loop is the dominant cost after precomputation | HIGH | Requires thread-local accumulators |

**Analysis:**

After per-jump precomputation (TS-01), the remaining per-step work is the rho_jump accumulation:

```julia
for w in energy_labels
    oft!(scratch.jump_oft, jump.in_eigenbasis, bohr_freqs, w, inv_4sigma2)
    rate2 = base_prefactor * transition(w)
    mul!(scratch.sandwich_tmp, evolving_dm, scratch.jump_oft')  # rho * L_w'
    mul!(scratch.rho_jump, scratch.jump_oft, scratch.sandwich_tmp, delta*rate2, 1.0)  # += L_w * (rho*L_w')
end
```

Each iteration involves 2 BLAS gemm calls (O(d^3) each). For num_energy_bits=8, there are ~256 energy labels (128 for half-grid with hermitian jumps). This is 256 * 2 = 512 gemm calls per step.

**Threading approach:**

Partition the energy labels across threads. Each thread has its own `rho_jump_local` accumulator. After the loop, sum the thread-local accumulators:

```julia
# Thread-local accumulators
rho_jump_locals = [zeros(CT, dim, dim) for _ in 1:nthreads()]
scratch_locals = [ThermalizeScratch(CT, dim) for _ in 1:nthreads()]

Threads.@threads for w_idx in 1:length(energy_labels_to_process)
    tid = Threads.threadid()
    w = energy_labels_to_process[w_idx]
    sc = scratch_locals[tid]
    # ... compute OFT, sandwich, accumulate into rho_jump_locals[tid] ...
end

# Reduce
for t in 1:nthreads()
    scratch.rho_jump .+= rho_jump_locals[t]
end
```

**Risk assessment:** This requires careful BLAS thread management (set BLAS threads to 1 during Julia-threaded w-loop to avoid oversubscription). The overhead of thread spawning and reduction may negate the benefit for small grids (num_energy_bits <= 6). Recommend gating this behind a `threaded_omega::Bool=false` flag and only enabling for num_energy_bits >= 8.

**Recommendation:** Defer to later in the milestone or make optional. The per-jump precomputation (TS-01/TS-02) provides the biggest speedup. Multi-threaded w-loop is a secondary optimization for large grids only.

---

### DIFF-02: Effective Rate Diagnostic for Burn-In Validation

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| Compute and return lambda_eff(t) alongside the fit | Provides a model-free diagnostic that validates the exponential fit quality | LOW | Leverages noiseless DM trace distance data |

**What to build:**

From the error catalogue (Section 6, Recommended Estimator, step (b)):

```julia
function effective_rate(trace_dists::Vector{Float64}, time_steps::Vector{Float64}; lag::Int=3)
    n = length(trace_dists)
    t_vals = time_steps[1:n-lag]
    lambda_eff = Vector{Float64}(undef, n-lag)
    for i in 1:n-lag
        d_now = trace_dists[i]
        d_later = trace_dists[i+lag]
        if d_now > 0 && d_later > 0
            tau = time_steps[i+lag] - time_steps[i]
            lambda_eff[i] = -log(d_later / d_now) / tau
        else
            lambda_eff[i] = NaN
        end
    end
    return t_vals, lambda_eff
end
```

For DM trace distance (noiseless), the effective rate plot shows:
- **Plateau at lambda_gap** in the exponential regime
- **Higher values** during the transient (where fast modes contribute)
- **Instability** near convergence (d(t) ~ C, log(d/d) ~ 0/0)

This is the "Rosetta Stone" diagnostic from the spectral gap refinement instructions. Including it in the `MixingTimeEstimate` output makes the mixing time result self-validating.

**Add to MixingTimeEstimate:**

```julia
struct MixingTimeEstimate
    # ... existing fields ...
    effective_rate_times::Vector{Float64}
    effective_rate_values::Vector{Float64}
end
```

---

### DIFF-03: Fit Quality Gates with Actionable Warnings

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| Structured warnings when fit quality is poor | Prevents users from trusting bad mixing time estimates | LOW | Critical for thesis-quality results |

**Quality gates to implement:**

1. **R-squared < 0.95:** Warning: "Single-exponential fit explains only X% of variance. Consider running longer (more data in exponential regime) or increasing skip_initial (transient contamination)."

2. **Fitted offset C > 0.1 * fitted amplitude A:** Warning: "Large offset C=X suggests the trace distance has not entered the exponential regime, or the fixed point accuracy is limiting. The offset should be approximately 0 for trace distance to the exact Gibbs state."

3. **Fitted gap negative or zero:** Error: "Non-physical negative decay rate. The trace distance may be increasing, suggesting the CPTP channel does not have the Gibbs state as its fixed point."

4. **CI width > 100% of gap:** Warning: "Spectral gap uncertainty (X%) is very large. The exponential regime may be too short for reliable estimation."

5. **Extrapolation sanity check:** When extrapolating, if `t_mix_extrapolated > 100 * t_data_max`, warning: "Extrapolation extends X times beyond observed data. Confidence in extrapolated mixing time is low."

---

## Anti-Features

Features to explicitly NOT build.

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| Two-exponential fit for trace distance | Adds 5-parameter fit complexity; for DM trace distance (noiseless), skip_initial + single-exp is sufficient. Two-exp was needed for noisy trajectory observables. | Use single-exponential with burn-in skip. If R^2 < 0.95, increase skip_initial. |
| Bootstrap uncertainty on DM trace distance | DM evolution is deterministic -- there is no statistical noise to bootstrap over. Bootstrap is for trajectory-averaged quantities. | Report LsqFit confidence intervals from the Jacobian-based covariance estimate. |
| Richardson extrapolation in delta for mixing time | Requires running thermalization at multiple delta values. This is a diagnostic tool, not a production feature. Mixing time estimation should work at a single delta. | Document that delta-convergence can be checked by running at delta/2 and comparing. |
| Observable-based mixing time estimation | The DM simulator tracks trace distance directly -- this is the gold standard. Observable decay is a proxy used when trace distance is unavailable (trajectory simulations). | Use trace distance directly. Observable fitting stays in staging/ for trajectory gap estimation. |
| Automatic delta selection for run_thermalize | Adds complexity without clear benefit. The user should choose delta based on domain knowledge (delta << 1 for CPTP validity). | Document recommended delta ranges in the function docstring. |
| GPU acceleration for DM operations | Not needed for current system sizes (n <= 12, dim <= 4096). BLAS already handles dim=4096 matrices efficiently. | Multi-threaded BLAS is sufficient. |

---

## Feature Dependencies

```
TS-01 (Per-jump precomputation)
  |
  +---> TS-02 (Precomputed CPTP application)
  |       |
  |       +---> DIFF-01 (Multi-threaded w-loop) -- only rho_jump remains in loop
  |
  +---> TS-05 (BLAS thread management) -- independent but test together

TS-03 (save_every) -- fully independent, can be done first

TS-04 (Mixing time estimation)
  |
  +---> requires staging/fitting.jl to be promoted to active code
  +---> uses ThermalizeResults output (trace_distances, time_steps)
  +---> DIFF-02 (Effective rate diagnostic) -- enhances TS-04
  +---> DIFF-03 (Quality gates) -- enhances TS-04
```

**Critical path:** TS-03 (easiest, do first) -> TS-01 -> TS-02 -> TS-05 -> TS-04

**Parallelizable:** TS-03, TS-05 are independent of each other and of TS-01.

---

## MVP Recommendation

### Phase 1: Quick Wins (Low Risk, High Visibility)

1. **TS-03: save_every** -- 30 min. Copy the pattern from trajectory simulator. Immediate usability improvement.

### Phase 2: Performance Core

2. **TS-01: Per-jump precomputation** -- 2 hrs. Follow the exact pattern from `_build_trajectory_workspace`. Add precomputed Rs, K0s, U_residuals to the thermalization setup.

3. **TS-02: Precomputed CPTP application** -- 1 hr. Replace `_finalize_kraus_step!` with `_apply_precomputed_channel!` using precomputed K0/U_residual. This eliminates the per-step eigendecomposition.

4. **TS-05: BLAS thread management** -- 30 min. Explicit BLAS thread control in run_thermalize.

### Phase 3: Mixing Time Estimation

5. **TS-04: Mixing time estimation** -- 3 hrs. Promote fitting.jl from staging, write `estimate_mixing_time`, add `MixingTimeEstimate` struct, implement extrapolation logic with quality gates.

6. **DIFF-02: Effective rate diagnostic** -- 1 hr. Add lambda_eff computation to MixingTimeEstimate output.

7. **DIFF-03: Quality gates** -- 1 hr. Add structured warnings for poor fit quality.

### Defer:

- **DIFF-01: Multi-threaded w-loop** -- Defer until benchmarking shows the w-loop is the bottleneck after precomputation. For typical parameter regimes (num_energy_bits=6), the loop has ~64 iterations which may not benefit from threading.

---

## Expected Behavior Specification: Mixing Time Estimation

### Input

`ThermalizeResults` from `run_thermalize` containing:
- `trace_distances::Vector{Float64}` -- d(t_i) = ||rho(t_i) - rho_beta||_1
- `time_steps::Vector{Float64}` -- t_i values

### Model

```
d(t) = A * exp(-gap * t) + C
```

where:
- `A > 0` is the initial amplitude (related to initial state distance from Gibbs)
- `gap > 0` is the effective spectral gap (the fitted parameter of interest)
- `C >= 0` is the asymptotic offset (should be ~0 for exact fixed point, ~1e-8 for TrotterDomain)

### Fit Procedure

1. **Skip burn-in:** Drop the first `skip_initial` fraction of data points. Default: 0.2 (20%).

2. **Initial guess:** Use `_log_linear_initial_guess` from fitting.jl:
   - Estimate C from last 20% of remaining data
   - Log-linear regression on log(d(t) - C) vs t gives A and gap

3. **Levenberg-Marquardt fit:** Via LsqFit.jl with bounds `gap >= 0`.

4. **Compute mixing time from data:** `t_mix_fitted = time at which fitted curve crosses convergence_cutoff`

5. **Extrapolate (optional):** `t_mix_extrapolated = -ln(target_epsilon / A) / gap`

### Output

`MixingTimeEstimate` struct containing all fitted parameters, quality metrics, and optionally the extrapolated mixing time.

### Validation Criteria

For a known system (e.g., n=4 EnergyDomain with exact gap from `run_lindblad`):
- Fitted gap should agree with exact gap within 10% (single-exponential approximation)
- R-squared should be > 0.99 for clean DM trace distance data
- Extrapolated t_mix should agree with actually observed convergence time within 20%

### Edge Cases

| Case | Expected Behavior |
|------|------------------|
| Insufficient data (< 10 points after skip) | Error: "Not enough data points for fitting" |
| Trace distance not monotonically decreasing | Warning: possible non-convergence or oscillatory modes |
| All trace distances below convergence_cutoff | Return gap=NaN, t_mix=0 (already converged) |
| Very short run (< 3 relaxation times) | Warning: "Run may be too short for reliable fitting" |
| Offset C > epsilon target | Warning: "Fixed point accuracy (C=X) exceeds target epsilon. Extrapolation not meaningful." |

---

## Sources

- **Codebase:** Direct reading of all source files in `/Users/bence/code/QuantumFurnace.jl/src/`, especially `furnace.jl`, `jump_workers.jl`, `furnace_utensils.jl`, `trajectories.jl`, `staging/fitting.jl`, `staging/gap_estimation.jl`
- **Chen et al. 2023:** Theorem III.1 (weak measurement CPTP channel), Proposition II.3 (mixing time from Hermitian gap), Corollary II.1 (fixed point accuracy)
- **Ramkumar & Soleimanifar 2024:** Theorems 2.1-2.3 (spectral gap bounds for graph-local and random jumps)
- **Error catalogue:** `supplementary-informations/error_catalogue_spectral_gap_estimation.md` -- Errors 1, 6 directly relevant to burn-in and fitting window selection
- **Spectral gap refinement instructions:** `supplementary-informations/spectral-gap-refinements-instructions.md` -- effective rate diagnostic, fitting window optimization
- **Project structure:** `supplementary-informations/quantumfurnace-structure.md` -- current architecture constraints and future plans
- **Confidence level:** HIGH -- all recommendations verified against existing code patterns and theoretical foundations
