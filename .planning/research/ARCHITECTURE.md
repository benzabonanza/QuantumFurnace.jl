# Architecture Research: v2.1 Speedup & Mixing Time Integration

**Domain:** QuantumFurnace.jl performance optimization and feature extension
**Researched:** 2026-03-01
**Confidence:** HIGH (analysis based on full source code reading, not external docs)

## Current Architecture Summary

The v2.0 architecture has five core layers for the `run_thermalize` path:

```
run_thermalize (furnace.jl:143-223)
    |
    |-- precompute phase (one-time)
    |       _precompute_data (furnace_utensils.jl:30-141)
    |       _precompute_coherent_unitary (coherent.jl:58-100)
    |
    |-- main loop (per-step)
    |       _jump_contribution! (jump_workers.jl:194-440)
    |           domain-dispatched: BohrDomain, EnergyDomain, TimeDomain/TrotterDomain
    |           accumulates R and rho_jump over omega-loop
    |           calls _finalize_kraus_step! at end
    |       _finalize_kraus_step! (jump_workers.jl:172-192)
    |           calls _build_cptp_channel (furnace_utensils.jl:183-200)
    |           K0, U_residual from R via eigen decomposition
    |           rho_next = K0 * rho * K0' + rho_jump + U_res * rho * U_res'
    |
    |-- return ThermalizeResults (structs.jl:201-207)
```

### Key Hot-Path Bottleneck: _build_cptp_channel Called Every Step

The critical performance issue is in `_finalize_kraus_step!` (jump_workers.jl:178):

```julia
(; K0, U_residual) = _build_cptp_channel(scratch.R, delta)
```

This calls `_build_cptp_channel` (furnace_utensils.jl:183-200) which:
1. Allocates K0 = I - alpha * R (dim x dim allocation)
2. Computes R^2 via matrix multiply
3. Computes S = (2*alpha - delta)*R - alpha^2 * R^2
4. Calls `eigen(Hermitian(S))` -- O(dim^3) eigendecomposition
5. Computes U_residual via Diagonal * eigvecs'

**This eigendecomposition happens every single thermalization step.** For a single randomly-selected jump, R is computed fresh each step from the omega-loop, then K0 and U_residual are recomputed from R. But when only one jump is picked per step, R^a is deterministic for a given jump `a` -- it depends only on the jump operator and precomputed domain data, not on the evolving state.

### Existing Precedent: Trajectory Workspace Already Has Per-Jump Precomputation

The trajectory code (trajectories.jl:60-152) already precomputes per-jump K0^a, U_residual^a:

```julia
# trajectories.jl:104-122
@inbounds for a in 1:n_jumps
    _precompute_R([jumps_for_diss[a]], ham_or_trott, config, precomputed_data, builder_scratch)
    R_a = copy(builder_scratch.R)
    R_a .*= (1.0 / p_jump)

    (; K0, U_residual) = _build_cptp_channel(R_a, delta)

    Rs[a] = R_a
    K0s[a] = K0
    U_residuals[a] = U_residual
end
```

These are stored in `Workspace.Rs`, `Workspace.K0s`, `Workspace.U_residuals`, and `Workspace.U_Bs` (per-jump coherent unitaries). This is exactly the pattern needed for `run_thermalize` speedup.

## Integration Analysis: New Features

### Feature 1: Per-Jump K0^a, U_residual^a Precomputation

**Current flow (furnace.jl:191-209, jump_workers.jl:194-440):**

```
for step in 1:num_steps
    idx = rand(rng, 1:length(jumps))
    jump = jumps[idx]
    _jump_contribution!(evolving_dm, jump, ..., scratch)
        # Inside: accumulates R from omega-loop
        # Then: _finalize_kraus_step!(evolving_dm, delta, scratch)
        #   -> _build_cptp_channel(scratch.R, delta)  # EXPENSIVE: eigen() per step
```

**Proposed flow:**

```
# Precompute phase (new, one-time):
for a in 1:n_jumps
    _precompute_R([jumps[a]], ham_or_trott, config, precomputed_data, scratch)
    R_a = copy(scratch.R) .* (1.0 / p_jump)
    (; K0, U_residual) = _build_cptp_channel(R_a, delta)
    K0s[a], U_residuals[a] = K0, U_residual
end

# Hot loop (modified):
for step in 1:num_steps
    idx = rand(rng, 1:length(jumps))
    _apply_coherent_unitary!(evolving_dm, U_Bs[idx], scratch)
    fill!(scratch.rho_jump, 0)
    _accumulate_rho_jump_only!(scratch.rho_jump, evolving_dm, jumps[idx], ...)
    _apply_precomputed_channel!(evolving_dm, K0s[idx], U_residuals[idx], scratch)
```

**Why the omega-loop for rho_jump still runs per step:** The rho_jump term is `delta * sum_w rate^2 * A_w * rho * A_w'` -- it depends on the current evolving_dm (the rho in the sandwich). Only the R accumulation (`R += rate^2 * A_w' * A_w`) and the resulting K0/U_residual are rho-independent and precomputable. The per-step omega-loop is cut roughly in half: only rho_jump sandwiches, no R accumulation.

**Where to store precomputed data:** Do NOT extend the Workspace struct. The `Workspace{S,D,C,T}` struct already has 38 fields and is shared across 4 simulation types. Instead, store K0s/U_residuals as local variables inside `run_thermalize`, matching the existing pattern where `coherent_unitaries` (furnace.jl:180) is already a local variable, not a struct field.

Specifically:
- `K0s::Vector{Matrix{CT}}` -- local to run_thermalize
- `U_residuals::Vector{Matrix{CT}}` -- local to run_thermalize
- `coherent_unitaries` -- already exists as local (furnace.jl:180)

**Modified components:**

| File | Location | Change |
|------|----------|--------|
| `furnace.jl:143-223` | `run_thermalize` | Add precompute loop before main loop; change main loop to use precomputed channel |
| `jump_workers.jl` | new function | Add `_apply_precomputed_channel!` |
| `jump_workers.jl:194-440` | `_jump_contribution!` (all 3 Thermalize variants) | Split into rho_jump-only accumulation + precomputed channel application; or add new `_accumulate_rho_jump_only!` functions |
| `furnace_utensils.jl` | existing file | Move `_precompute_R` here from trajectories.jl; add BohrDomain variant |

**New functions:**

| Function | Purpose | Location |
|----------|---------|----------|
| `_apply_precomputed_channel!(evolving_dm, K0, U_residual, scratch)` | Apply K0*rho*K0' + rho_jump + U_res*rho*U_res' using scratch buffers | `jump_workers.jl` |
| `_accumulate_rho_jump_only!(rho_jump, evolving_dm, jump, ...)` | Domain-dispatched omega-loop for rho_jump sandwich only (no R accumulation) | `jump_workers.jl` (3 variants: BohrDomain, EnergyDomain, TimeDomain/TrotterDomain) |

**Code reuse:** The `_precompute_R` function already exists in trajectories.jl:193-300 for EnergyDomain and TimeDomain/TrotterDomain. It should be moved to `furnace_utensils.jl` as a shared utility alongside `_precompute_data` and `_build_cptp_channel`.

**BohrDomain gap:** `_precompute_R` in trajectories.jl has no BohrDomain variant (trajectories never use BohrDomain). For `run_thermalize`, BohrDomain IS used. Need to add a BohrDomain `_precompute_R` method. The R accumulation logic is embedded in `_jump_contribution!` (jump_workers.jl:226-282):

```julia
@inbounds for (k, nu_2) in pairs(bohr_keys)
    @. scratch.jump_oft = alpha(hamiltonian.bohr_freqs, nu_2) * jump.in_eigenbasis
    # ... accumulate both R and rho_jump interleaved
end
```

For precomputation, extract only the R accumulation part into a new `_precompute_R` method dispatched on BohrDomain.

### Feature 2: save_every for Trace Distance Computation

**Current structure (furnace.jl:191-209):**

```julia
for step in 1:num_steps
    # ... jump contribution ...
    dist = trace_distance_h(Hermitian(evolving_dm), gibbs)
    push!(trace_distances, dist)
    @printf("Dist to Gibbs: %s\n", dist)
    if dist < convergence_cutoff
        break
    end
end
```

`trace_distance_h` involves an eigendecomposition (O(dim^3)), computed EVERY step. For long runs, this is wasteful.

**Existing precedent:** The trajectory code uses `save_every` extensively (trajectories.jl:396-464):

```julia
if step % save_every == 0
    save_idx = div(step, save_every) + 1
    # ... save observable measurements ...
end
```

**Proposed integration:** Add `save_every::Int = 1` keyword to `run_thermalize`:

```julia
function run_thermalize(
    jumps, config, hamiltonian, trotter=nothing;
    initial_dm=nothing, rng=Random.default_rng(),
    rescale_by_inv_prob=true,
    save_every::Int = 1,  # NEW
)
```

**Modified loop:**

```julia
num_saves = div(num_steps, save_every) + 1
trace_distances = Vector{T}(undef, num_saves)
trace_distances[1] = trace_distance_h(Hermitian(evolving_dm), gibbs)
save_idx = 1

for step in 1:num_steps
    # ... jump contribution ...
    if step % save_every == 0
        save_idx += 1
        dist = trace_distance_h(Hermitian(evolving_dm), gibbs)
        trace_distances[save_idx] = dist
        if dist < convergence_cutoff
            break
        end
    end
end
trace_distances = trace_distances[1:save_idx]
time_steps = [(s - 1) * save_every * config.delta for s in 1:save_idx]
```

**ThermalizeResults impact:** No struct change needed. `trace_distances` and `time_steps` already have matching lengths by construction. Callers get fewer data points when save_every > 1, which is the intended behavior.

**Backward compatibility:** Default `save_every=1` preserves current behavior exactly.

### Feature 3: Multi-threaded BLAS for Channel Application

**After Feature 1:** The per-step `_apply_precomputed_channel!` consists of 4 `mul!` calls (BLAS.gemm) on dim x dim matrices. For dim >= 64 (n >= 6 qubits), multi-threaded BLAS provides speedup.

**Current trajectory behavior:** The trajectory path sets `BLAS.set_num_threads(1)` during threaded execution (trajectories.jl:494-496) to avoid thread oversubscription. But `run_thermalize` is single-threaded at the Julia level, so multi-threaded BLAS on the matrix multiplies is beneficial and safe.

**Implementation:** Ensure `BLAS.set_num_threads()` is at its default (number of CPU threads) during `run_thermalize`. Since the trajectory path restores BLAS threads in try/finally, this should already be the case. Verify and document.

**No code changes needed** unless `run_thermalize` is called from within a trajectory batch (which it is not -- they are separate entry points).

### Feature 4: Multi-threaded omega-loops (Precomputation Only)

**After Feature 1:** The omega-loop only runs during precomputation (once per jump), not in the hot loop. Multi-threading the omega-loop in `_precompute_R` gives a one-time speedup on construction.

**Implementation:** Thread the energy_labels loop with per-thread R accumulators:

```julia
function _precompute_R_threaded(jump, ham_or_trott, config, precomputed_data)
    n_threads = Threads.nthreads()
    R_threads = [zeros(CT, dim, dim) for _ in 1:n_threads]
    scratch_threads = [ThermalizeScratch(CT, dim) for _ in 1:n_threads]

    Threads.@threads for w_idx in eachindex(energy_labels)
        tid = Threads.threadid()
        # accumulate into R_threads[tid] using scratch_threads[tid]
    end
    R = sum(R_threads)
    hermitianize!(R)
    return R
end
```

**The per-step rho_jump loop could also be threaded**, but this is riskier: the evolving_dm is shared read-only (safe), rho_jump accumulation needs atomic or per-thread reduction (manageable). However, the per-step loop is already smaller after Feature 1 (no R accumulation), and the BLAS calls within it may already use multi-threaded BLAS (Feature 3). Threading at the Julia level on top of multi-threaded BLAS risks oversubscription. **Defer per-step omega-loop threading** until profiling shows it is the bottleneck.

### Feature 5: Mixing Time Estimation via Exponential Fit

**Existing code:** `fit_exponential_decay` in `staging/fitting.jl` fits `A * exp(-gap * t) + C`. `FitResult` struct stores gap, amplitude, offset, CI, R-squared.

**Integration approach:** Post-processing function, NOT integrated into `run_thermalize`. Reasons:
1. Couples simulation logic with analysis logic if embedded
2. Fit parameters (skip_initial, threshold) are analysis choices, not simulation parameters
3. User may want to re-fit with different parameters without re-running

```julia
function estimate_mixing_time(result::ThermalizeResults;
    skip_initial::Float64 = 0.0,
    threshold::Float64 = 1e-4,
    extrapolate::Bool = false,
)
    fit = fit_exponential_decay(result.time_steps, result.trace_distances;
        skip_initial=skip_initial)

    if fit.gap > 0 && fit.amplitude > 0 && (threshold - fit.offset) / fit.amplitude > 0
        t_mix = -log((threshold - fit.offset) / fit.amplitude) / fit.gap
    else
        t_mix = NaN  # fit doesn't support meaningful extrapolation
    end

    return (; mixing_time=t_mix, fit=fit, extrapolated=extrapolate, threshold=threshold)
end
```

**File location:** Move `fitting.jl` from `staging/` to active `src/`. Add `estimate_mixing_time` alongside it.

**Return type:** NamedTuple, not a new struct. This is a lightweight post-processing step.

## Component Modification Summary

### Files Modified

| File | Lines | Change | Risk |
|------|-------|--------|------|
| `furnace.jl` | 143-223 | Add precompute loop, modify main loop, add save_every kwarg | MEDIUM |
| `jump_workers.jl` | new code | Add `_apply_precomputed_channel!`, add `_accumulate_rho_jump_only!` (3 domain variants) | MEDIUM |
| `furnace_utensils.jl` | append | Move `_precompute_R` here, add BohrDomain variant | LOW |
| `QuantumFurnace.jl` | exports | Add `estimate_mixing_time`, `FitResult`, `fit_exponential_decay`, include fitting.jl | LOW |

### Files Moved from Staging

| From | To | Purpose |
|------|----|---------|
| `staging/fitting.jl` | `src/fitting.jl` | Activate exponential fitting + add mixing time estimation |

### Files NOT Modified

| File | Reason |
|------|--------|
| `structs.jl` | No new structs needed; ThermalizeResults unchanged; K0s/U_residuals are local variables in run_thermalize |
| `results.jl` | ThermalizeResults serialization unchanged (same fields, fewer trace_distance entries with save_every) |
| `coherent.jl` | `_precompute_coherent_unitary` already produces per-jump U_B -- unchanged |
| `trajectories.jl` | Keeps its own `_precompute_R` usage internally; shared function moves to furnace_utensils.jl |
| `krylov_workspace.jl` | Krylov path unaffected |
| `krylov_matvec.jl` | Krylov path unaffected |
| `krylov_eigsolve.jl` | Krylov path unaffected |

## Data Flow Diagrams

### Current run_thermalize Data Flow (v2.0)

```
run_thermalize()
    |
    +-- _precompute_data(config, ham_or_trott)
    |       -> precomputed_data (transition, gamma_norm_factor, energy_labels, ...)
    |
    +-- _precompute_coherent_unitary(jumps, ...)
    |       -> coherent_unitaries: Vector{Matrix{CT}} (per-jump U_B)
    |
    +-- ThermalizeScratch(CT, dim)
    |       -> scratch (R, rho_jump, sandwich_tmp, rho_next, jump_oft, LdagL, ...)
    |
    +-- for step in 1:num_steps  -------- HOT LOOP --------
    |       |
    |       +-- idx = rand(rng, 1:n_jumps)
    |       |
    |       +-- _jump_contribution!(evolving_dm, jumps[idx], ...)
    |       |       |
    |       |       +-- _apply_coherent_unitary!(evolving_dm, U_B[idx], scratch)
    |       |       +-- fill!(scratch.R, 0); fill!(scratch.rho_jump, 0)
    |       |       +-- for w in energy_labels  ---- OMEGA LOOP ----
    |       |       |       scratch.R += rate^2 * A_w' * A_w       <-- rho-INDEPENDENT
    |       |       |       scratch.rho_jump += delta * rate^2 * A_w * rho * A_w'  <-- rho-DEPENDENT
    |       |       |
    |       |       +-- hermitianize!(scratch.R)
    |       |       +-- _finalize_kraus_step!(evolving_dm, delta, scratch)
    |       |               |
    |       |               +-- _build_cptp_channel(scratch.R, delta)
    |       |               |       eigen(Hermitian(S))    <-- O(dim^3) PER STEP
    |       |               |       -> K0, U_residual       (ALLOCATES EVERY STEP)
    |       |               |
    |       |               +-- rho_next = K0*rho*K0' + rho_jump + U_res*rho*U_res'
    |       |
    |       +-- dist = trace_distance_h(evolving_dm, gibbs)  <-- EVERY STEP
    |
    +-- return ThermalizeResults(evolving_dm, trace_distances, time_steps, metadata)
```

### Proposed run_thermalize Data Flow (v2.1)

```
run_thermalize(; save_every=1)
    |
    +-- _precompute_data(config, ham_or_trott)
    |       -> precomputed_data
    |
    +-- _precompute_coherent_unitary(jumps, ...)
    |       -> coherent_unitaries (per-jump U_B, unchanged)
    |
    +-- PER-JUMP CHANNEL PRECOMPUTATION (NEW, one-time)
    |       for a in 1:n_jumps:
    |           _precompute_R([jumps[a]], ..., scratch)  <-- omega-loop once per jump
    |           R_a = copy(scratch.R) .* (1/p_jump)
    |           (; K0, U_residual) = _build_cptp_channel(R_a, delta)
    |           K0s[a] = K0; U_residuals[a] = U_residual
    |       -> K0s, U_residuals: Vector{Matrix{CT}} (local vars)
    |
    +-- ThermalizeScratch(CT, dim)
    |       -> scratch (rho_jump, sandwich_tmp, rho_next, jump_oft, LdagL)
    |          (scratch.R no longer needed in hot loop)
    |
    +-- for step in 1:num_steps  -------- HOT LOOP (FASTER) --------
    |       |
    |       +-- idx = rand(rng, 1:n_jumps)
    |       |
    |       +-- _apply_coherent_unitary!(evolving_dm, U_Bs[idx], scratch)
    |       |
    |       +-- fill!(scratch.rho_jump, 0)
    |       +-- _accumulate_rho_jump_only!(scratch, evolving_dm, jumps[idx],
    |       |       config, precomputed_data)
    |       |       for w in energy_labels:      <-- still per-step (rho-dependent)
    |       |           rho_jump += delta * rate^2 * A_w * rho * A_w'
    |       |       (NO R accumulation -- eliminated)
    |       |
    |       +-- _apply_precomputed_channel!(evolving_dm, K0s[idx],
    |       |       U_residuals[idx], scratch)
    |       |       rho_next = K0*rho*K0' + rho_jump + U_res*rho*U_res'
    |       |       (NO eigendecomposition -- 4 mul! calls only)
    |       |
    |       +-- if step % save_every == 0:       <-- NOT every step
    |               dist = trace_distance_h(evolving_dm, gibbs)
    |
    +-- return ThermalizeResults(evolving_dm, trace_distances, time_steps, metadata)
    |
    |   (Post-processing, separate call):
    +-- estimate_mixing_time(result; threshold, extrapolate)
            -> (; mixing_time, fit, extrapolated)
```

### Key Differences (v2.0 vs v2.1)

| Aspect | v2.0 (current) | v2.1 (proposed) |
|--------|----------------|-----------------|
| R accumulation | Every step, in omega-loop | Once per jump, during precomputation |
| K0/U_residual computation | `eigen()` every step | Precomputed, lookup by jump index |
| Omega-loop per step | Full (R + rho_jump) | Half (rho_jump only) |
| trace_distance_h calls | Every step | Every save_every steps |
| Channel application | `_finalize_kraus_step!` (allocates K0, U_res) | `_apply_precomputed_channel!` (zero allocation) |
| Mixing time | Not available | Post-processing via exponential fit |

## Performance Impact Analysis

### Per-step savings from Feature 1

For each step, the eliminated operations are:
1. **R accumulation over omega-loop:** |energy_labels| iterations, each doing a `mul!(LdagL, A', A)` and element-wise add. For EnergyDomain with n=6: ~128 frequencies, each with 64x64 BLAS.gemm + elementwise add.
2. **eigendecomposition of S:** `eigen(Hermitian(S))` where S is dim x dim. For dim=64 (n=6), this is O(dim^3) and dominates.
3. **K0 and U_residual construction:** I - alpha*R (elementwise), R^2 (gemm), S construction, Diagonal * eigvecs' (gemm).
4. **All allocations in `_build_cptp_channel`:** K0, R2, S, eig.vectors -- 4 dim x dim matrix allocations per step eliminated.

The remaining per-step cost is:
1. **rho_jump accumulation:** Still |energy_labels| iterations with A_w*rho*A_w' sandwiches (two gemms per frequency).
2. **Channel application:** 4 mul! calls (K0*rho*K0' + rho_jump + U_res*rho*U_res').

### Estimated speedup by system size

| System Size | dim | eigen() cost | Omega-loop (R half) | Estimated total speedup |
|-------------|-----|-------------|---------------------|------------------------|
| n=4 | 16 | Small | Small | ~1.5-2x |
| n=6 | 64 | Significant | Moderate | ~2-4x |
| n=8 | 256 | Dominant | Large | ~5-10x |
| n=10 | 1024 | Massive | Massive | Essential for feasibility |

### Per-step savings from save_every

`trace_distance_h` calls `eigen(Hermitian(rho - gibbs))` -- another O(dim^3) eigendecomposition per evaluation. With `save_every=100` on a 10,000-step run, this drops from 10,000 to 100 evaluations.

## Build Order (Dependency-Driven)

### Phase A: Per-Jump Precomputation (Core Speedup)

**Depends on:** Nothing (start immediately)
**Enables:** Phase C (BLAS threading), Phase D (mixing time benefits from faster runs)

**Steps:**
1. Move `_precompute_R` from trajectories.jl to furnace_utensils.jl (EnergyDomain and Time/TrotterDomain variants). Keep trajectories.jl calling the shared function.
2. Add BohrDomain `_precompute_R` variant (extract R accumulation from jump_workers.jl:226-282).
3. Add `_apply_precomputed_channel!` to jump_workers.jl.
4. Add `_accumulate_rho_jump_only!` for each domain (3 variants: BohrDomain, EnergyDomain, Time/TrotterDomain). These are extracted from existing `_jump_contribution!` by removing the R accumulation lines.
5. Modify `run_thermalize` to: (a) precompute K0s/U_residuals per jump, (b) change main loop to call `_accumulate_rho_jump_only!` + `_apply_precomputed_channel!`.
6. Tests: verify numerical equivalence of new vs old path (trace distances match to machine precision for all domains).

**Risk:** MEDIUM -- restructuring 3 domain-dispatched functions. Mitigated by: testing equivalence against old code before removing it.

### Phase B: save_every Integration

**Depends on:** Phase A (both modify the loop, do them together or A then B)
**Enables:** Phase D (mixing time estimation benefits from save_every for practical use)

**Steps:**
1. Add `save_every::Int = 1` kwarg to `run_thermalize`.
2. Modify loop to check `step % save_every == 0`.
3. Pre-allocate trace_distances with correct size.
4. Adjust time_steps to match save points only.
5. Tests: save_every=1 matches old behavior; save_every=N produces correct subset; convergence_cutoff still works.

**Risk:** LOW -- additive change, backward compatible with default save_every=1.

### Phase C: Multi-threaded BLAS + Precomputation Threading

**Depends on:** Phase A (omega-loop threading only helps during precomputation after Phase A)
**Enables:** Nothing (independent performance feature)

**Steps:**
1. Verify `BLAS.set_num_threads()` is at default (not 1) during run_thermalize. Document.
2. Optionally: thread the omega-loop in `_precompute_R` with per-thread accumulators.
3. Optionally: thread the rho_jump omega-loop with per-thread accumulators (requires BLAS single-thread to avoid oversubscription, profile first).
4. Tests: numerical equivalence; verify no thread safety issues.

**Risk:** LOW for BLAS threading (just verify default). MEDIUM for omega-loop threading (reduction pattern).

### Phase D: Mixing Time Estimation

**Depends on:** Phase B (save_every needed for practical use), fitting code activation
**Enables:** Nothing (end-user feature)

**Steps:**
1. Move fitting.jl from staging/ to src/.
2. Add `estimate_mixing_time` function.
3. Add module includes and exports.
4. Tests: known exponential decay produces correct mixing time; extrapolation handles edge cases (NaN for non-convergent fits).

**Risk:** LOW -- post-processing function, no impact on simulation core.

## Anti-Patterns to Avoid

### Anti-Pattern 1: Extending Workspace{S,D,C,T} with Thermalize-Only Fields

**What people do:** Add K0s, U_residuals as new fields to Workspace.
**Why it's wrong:** Workspace is shared across 4 simulation types and already has 38 fields. Adding Thermalize-specific fields bloats it further. The trajectory path stores these as trajectory-specific fields, which makes sense for Workspace{Trajectory} but not for a shared struct.
**Do this instead:** Store K0s/U_residuals as local variables inside run_thermalize. The precomputation functions can be shared; the storage location differs.

### Anti-Pattern 2: Merging Thermalize and Trajectory Precomputation Paths

**What people do:** Make run_thermalize call `_build_trajectory_workspace` to reuse its precomputation.
**Why it's wrong:** The trajectory workspace has completely different structure (per-thread scratch, state-vector vs density-matrix, different step function). Coupling them creates maintenance risk and conceptual confusion.
**Do this instead:** Share `_precompute_R` and `_build_cptp_channel` functions (already shared), but keep workspace construction and hot loops separate.

### Anti-Pattern 3: Eliminating the Per-Step omega-loop Entirely via Superoperator

**What people do:** Precompute the full Lindblad channel E(rho) = sum_w rate^2 A_w rho A_w' as a dim^2 x dim^2 superoperator matrix and apply it per step.
**Why it's wrong:** The superoperator is dim^2 x dim^2 (4096 x 4096 for n=6), far more expensive to store and apply than the sum of dim x dim sandwiches. Storage is O(dim^4), application is O(dim^4) matmul vs O(|energy_labels| * dim^3) for the sandwich loop.
**Do this instead:** Keep the per-step omega-loop for rho_jump. Only precompute the rho-independent parts (R, K0, U_residual).

### Anti-Pattern 4: Adding Mixing Time Estimation as a run_thermalize Keyword

**What people do:** Add `estimate_mixing_time::Bool=false` or `extrapolate::Bool=false` to run_thermalize.
**Why it's wrong:** Couples simulation logic with analysis logic. The fit parameters (skip_initial, threshold) are analysis choices, not simulation parameters. Makes the function signature unwieldy and the function responsible for too many things.
**Do this instead:** Keep `estimate_mixing_time` as a separate post-processing function that takes `ThermalizeResults` as input.

### Anti-Pattern 5: Modifying ThermalizeScratch to Remove R Field

**What people do:** Since R is no longer accumulated per-step, remove `scratch.R` from ThermalizeScratch.
**Why it's wrong:** ThermalizeScratch is used by `_precompute_R` during the precomputation phase. Also, removing a struct field is a breaking change for any code that constructs ThermalizeScratch.
**Do this instead:** Keep `scratch.R` in ThermalizeScratch. It is used during precomputation even though it is not used in the per-step hot loop. The unused per-step R accumulation becomes a scratch buffer for precomputation.

## Integration Points

### Internal Boundaries

| Boundary | Communication | Notes |
|----------|---------------|-------|
| `run_thermalize` <-> `_precompute_R` | Calls with jumps/config/precomputed_data, gets R matrix back | Same pattern as trajectory path |
| `run_thermalize` <-> `_accumulate_rho_jump_only!` | Passes scratch (mutable), evolving_dm (read), jump, config, precomputed_data | Domain-dispatched, 3 variants |
| `run_thermalize` <-> `_apply_precomputed_channel!` | Passes evolving_dm (mutated), K0/U_residual (read), scratch (mutable) | Simple, no domain dispatch needed |
| `estimate_mixing_time` <-> `fit_exponential_decay` | Passes time_steps + trace_distances from ThermalizeResults | Clean function call, no shared state |
| `_precompute_R` shared between furnace_utensils.jl and trajectories.jl | Same function, called from different contexts | trajectories.jl changes to call shared function |

### Function Extraction Detail: _accumulate_rho_jump_only!

For EnergyDomain (jump_workers.jl:318-359), the current `_jump_contribution!` interleaves R and rho_jump accumulation:

```julia
# Lines 326-334 (positive frequency):
rate2_pos = base_prefactor * transition(w)
mul!(scratch.LdagL, scratch.jump_oft', scratch.jump_oft)
@. scratch.R += rate2_pos * scratch.LdagL                    # R part (REMOVE)
mul!(scratch.sandwich_tmp, evolving_dm, scratch.jump_oft')
mul!(scratch.rho_jump, scratch.jump_oft, scratch.sandwich_tmp, config.delta * rate2_pos, 1.0)  # rho_jump part (KEEP)
```

The new `_accumulate_rho_jump_only!` keeps only the rho_jump lines. The `scratch.LdagL` multiply for R is removed. The `scratch.R` fill is removed.

For BohrDomain (jump_workers.jl:218-282), the structure is different -- R accumulation is interleaved differently:

```julia
# Lines 256-257 (rho_jump):
mul!(scratch.rho_jump, scratch.jump_oft, scratch.sandwich_tmp, scaled_delta, 1.0)
# Lines 260-281 (R):
scratch.R[j, q] += v * scratch.jump_oft[i, q]
```

Again, `_accumulate_rho_jump_only!` drops the R lines and the jump_weight_scaling computation for R.

## Sources

All analysis based on direct source code reading:
- `furnace.jl` lines 143-223 (run_thermalize)
- `jump_workers.jl` lines 172-440 (_finalize_kraus_step!, _jump_contribution! x3 Thermalize domains)
- `furnace_utensils.jl` lines 30-200 (_precompute_data, _build_cptp_channel)
- `trajectories.jl` lines 60-300 (_build_trajectory_workspace, _precompute_R, _copy_workspace_for_thread)
- `coherent.jl` lines 58-100 (_precompute_coherent_unitary)
- `structs.jl` lines 291-417 (ThermalizeScratch, Workspace)
- `staging/fitting.jl` full file (FitResult, fit_exponential_decay)
- `staging/gap_estimation.jl` full file (estimate_spectral_gap, save_every usage)

---
*Architecture research for: v2.1 Speedup & Mixing Time*
*Researched: 2026-03-01*
