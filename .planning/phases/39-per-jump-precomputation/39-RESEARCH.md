# Phase 39: Per-Jump Precomputation - Research

**Researched:** 2026-03-01
**Domain:** Per-jump CPTP channel precomputation for `run_thermalize` hot loop optimization
**Confidence:** HIGH

## Summary

The `run_thermalize` hot loop currently recomputes R, K0, and U_residual (including an eigendecomposition inside `_build_cptp_channel`) on every step for Energy, Time, and Trotter domains. Since R depends only on the jump operator, the Hamiltonian eigendata, and the config (not on the evolving density matrix rho), it can be precomputed once per jump at simulation start.

The trajectory engine (`_build_trajectory_workspace` in `trajectories.jl`) already implements exactly this pattern: it calls `_precompute_R` per jump, scales by `1/p_jump`, and stores the resulting `Rs`, `K0s`, `U_residuals` vectors in the `Workspace{Trajectory}` struct. The per-jump coherent unitaries are also precomputed via `_precompute_coherent_unitary`. This phase ports the same precomputation pattern to `run_thermalize`, eliminating the per-step eigendecomposition that currently lives inside `_finalize_kraus_step! -> _build_cptp_channel`.

**Primary recommendation:** Refactor `run_thermalize` to precompute per-jump K0s, U_residuals (and coherent unitaries, already done) before the hot loop, then replace `_jump_contribution!` + `_finalize_kraus_step!` with a simpler `_apply_precomputed_channel!` that uses pre-stored matrices. Reuse the existing `_precompute_R` functions from `trajectories.jl`.

## Standard Stack

### Core

No new libraries required. This is a pure algorithmic refactoring within the existing Julia codebase.

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| LinearAlgebra (stdlib) | Julia 1.x | `eigen()`, `mul!`, `Hermitian`, `Diagonal` | Core matrix operations for CPTP channel |
| FINUFFT.jl | existing | NUFFT prefactors for Time/Trotter domains | Already used in `_precompute_data` |

### Supporting

No new supporting libraries needed.

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Per-jump precompute (chosen) | Sum-of-R precompute only | Per-jump gives exact same CPTP channel per step (no mixing), but sum-of-R is what Krylov does and loses per-jump structure |
| Eigendecomposition for sqrt_psd | Cholesky-based sqrt | Eigendecomposition is already validated, handles non-PSD noise robustly via clamping |

## Architecture Patterns

### Current Architecture (Before This Phase)

```
run_thermalize (furnace.jl)
  |-- _precompute_data(config, ham_or_trott)       # Once at start
  |-- _precompute_coherent_unitary(...)             # Once at start (already precomputed!)
  |-- for step in 1:num_steps                       # HOT LOOP
  |     |-- idx = rand(rng, 1:length(jumps))
  |     |-- _jump_contribution!(evolving_dm, jump, ...)
  |     |     |-- fill!(scratch.R, 0)
  |     |     |-- for w in energy_labels             # omega loop: accumulate R and rho_jump
  |     |     |     |-- oft! or nufft_prefactor_view
  |     |     |     |-- mul!(scratch.LdagL, ...)
  |     |     |     |-- scratch.R += rate2 * LdagL    # R recomputed every step (SAME for same jump)
  |     |     |     |-- rho_jump += delta * rate2 * Aw * rho * Aw'
  |     |     |-- hermitianize!(scratch.R)
  |     |     |-- _finalize_kraus_step!(evolving_dm, delta, scratch)
  |     |           |-- _build_cptp_channel(scratch.R, delta)
  |     |           |     |-- eigen(Hermitian(S))      # EIGENDECOMPOSITION EVERY STEP!
  |     |           |-- K0 * rho * K0' + rho_jump + U_res * rho * U_res'
```

### Target Architecture (After This Phase)

```
run_thermalize (furnace.jl)
  |-- _precompute_data(config, ham_or_trott)       # Once at start
  |-- _precompute_coherent_unitary(...)             # Once at start (already done)
  |-- _precompute_per_jump_channels(...)            # NEW: Once at start
  |     |-- for each jump a:
  |     |     |-- _precompute_R([jump_a], ...)       # Reuse existing from trajectories.jl
  |     |     |-- R_a *= (1/p_jump)                  # Scale by 1/p_jump to match DM convention
  |     |     |-- _build_cptp_channel(R_a, delta)    # One eigen() per jump, not per step
  |     |     |-- Store K0s[a], U_residuals[a]
  |-- for step in 1:num_steps                       # HOT LOOP
  |     |-- idx = rand(rng, 1:length(jumps))
  |     |-- _apply_coherent_unitary!(...)             # Use precomputed coherent_unitaries[idx]
  |     |-- _accumulate_rho_jump!(rho_jump, jump, ...) # NEW: omega loop for rho_jump ONLY
  |     |-- _apply_precomputed_channel!(evolving_dm, K0s[idx], U_residuals[idx], rho_jump)
  |           |-- rho_next = K0 * rho * K0' + rho_jump + U_res * rho * U_res'
  |           |-- NO eigen()!
```

### Key Architectural Insight: R is rho-independent

The critical insight enabling this optimization is that R^a is independent of the evolving density matrix rho:

```
R^a = sum_w rate^2(w) * L_{a,w}^dagger * L_{a,w}
```

where `L_{a,w}` depends only on `jump.in_eigenbasis`, `bohr_freqs`, and `w` (config/Hamiltonian data), NOT on rho. The only rho-dependent part of `_jump_contribution!` is `rho_jump`:

```
rho_jump = delta * sum_w rate^2(w) * L_{a,w} * rho * L_{a,w}^dagger
```

Therefore:
1. R^a can be precomputed once per jump
2. K0^a = I - alpha * R^a can be precomputed once per jump
3. U_residual^a = sqrt_psd(S^a) can be precomputed once per jump (eliminates per-step eigen()!)
4. rho_jump must still be computed per step (depends on rho)

### Pattern: Existing Trajectory Precomputation (Reference Implementation)

The trajectory engine in `trajectories.jl` (lines 61-152 of `_build_trajectory_workspace`) already implements this exact pattern:

```julia
# From _build_trajectory_workspace (trajectories.jl):
builder_scratch = ThermalizeScratch(CT, dim)
for a in 1:n_jumps
    _precompute_R([jumps_for_diss[a]], ham_or_trott, config, precomputed_data, builder_scratch)
    R_a = copy(builder_scratch.R)
    R_a .*= (1.0 / p_jump)
    (; K0, U_residual) = _build_cptp_channel(R_a, delta)
    Rs[a] = R_a
    K0s[a] = K0
    U_residuals[a] = U_residual
end
```

This is directly reusable for `run_thermalize`.

### Pattern: Per-Step Channel Application

After precomputation, the per-step code simplifies to:

```julia
# Simplified hot loop body (pseudocode):
_apply_coherent_unitary!(evolving_dm, coherent_unitaries[idx], scratch)
_accumulate_rho_jump!(scratch.rho_jump, jump, ...)  # omega loop (rho-dependent)
# Apply precomputed channel:
mul!(scratch.sandwich_tmp, K0s[idx], evolving_dm)
mul!(scratch.rho_next, scratch.sandwich_tmp, K0s[idx]')
scratch.rho_next .+= scratch.rho_jump
mul!(scratch.sandwich_tmp, U_residuals[idx], evolving_dm)
mul!(scratch.rho_next, scratch.sandwich_tmp, U_residuals[idx]', 1.0, 1.0)
hermitianize!(scratch.rho_next)
copyto!(evolving_dm, scratch.rho_next)
```

### BohrDomain: Why No Per-Jump Precomputation

For BohrDomain, R^a has the form:

```
R^a = sum_{nu_2} A_{nu_2}^dagger * B_{nu_2}^a
```

where `B_{nu_2}^a = sum_{nu_1} alpha(nu_1, nu_2) * A_{nu_1}^a`. This involves summing over ALL unique Bohr frequencies for each jump operator. For a system with `dim` eigenvalues, there are O(dim^2) unique Bohr frequencies. At n=12 qubits (dim=4096), this is ~16 million unique frequencies -- precomputing per-Bohr-frequency data would be prohibitive.

However, the BohrDomain `_jump_contribution!` for Thermalize does NOT do eigendecomposition either (it calls `_finalize_kraus_step!` which does `_build_cptp_channel` with eigen()). So BohrDomain ALSO benefits from precomputation, but through a different _precompute_R function specific to the Bohr structure.

**Wait -- there is no `_precompute_R` for BohrDomain yet.** The existing `_precompute_R` functions in `trajectories.jl` only cover EnergyDomain and TimeDomain/TrotterDomain. A new `_precompute_R` for BohrDomain would need to be written, extracting the R-accumulation logic from the Thermalize BohrDomain `_jump_contribution!`.

### Anti-Patterns to Avoid

- **Modifying `_jump_contribution!` signatures in a way that breaks the Lindbladian path:** The Lindbladian versions of `_jump_contribution!` accumulate into a vectorized Liouvillian and do NOT call `_finalize_kraus_step!`. These must remain unchanged.
- **Storing R^a in the ThermalizeScratch:** R^a is per-jump, not per-step scratch. Store in a separate precomputed structure, not in the mutable scratch buffers.
- **Forgetting the `1/p_jump` scaling:** The trajectory engine scales R_a by `1/p_jump = n_jumps` because each operator channel is applied with probability `p_jump = 1/n_jumps`. The DM code uses `jump_weight_scaling = rescale_by_inv_prob ? (gamma_norm_factor / jump_prob) : gamma_norm_factor`. Must match.
- **Breaking the rho_jump accumulation:** The rho_jump calculation in the omega loop depends on rho and MUST remain in the per-step loop. Only R-accumulation moves out.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| R precomputation for Energy/Time/Trotter | New R computation code | `_precompute_R` from `trajectories.jl` | Already validated, tested via CPTP completeness tests |
| CPTP channel from R | New K0/U_residual derivation | `_build_cptp_channel` from `furnace_utensils.jl` | Already validated, handles PSD guard correctly |
| Per-jump coherent unitaries | New coherent precomputation | `_precompute_coherent_unitary` from `coherent.jl` | Already called by `run_thermalize` -- already precomputed! |

**Key insight:** The coherent unitaries are ALREADY precomputed in `run_thermalize` (line 180-181 of furnace.jl). The only missing precomputation is R^a -> K0^a, U_residual^a.

## Common Pitfalls

### Pitfall 1: Delta Scaling Mismatch Between DM and Trajectory

**What goes wrong:** The trajectory engine scales R_a by `1/p_jump` and uses bare `delta` in `_build_cptp_channel`. The DM code uses `jump_weight_scaling` which combines `gamma_norm_factor / jump_prob`. If the scaling is not matched exactly, the CPTP channel will be different and results will diverge.

**Why it happens:** The DM code applies `scaled_delta = config.delta * jump_weight_scaling` where `jump_weight_scaling = gamma_norm_factor / jump_prob`. The R accumulation in `_jump_contribution!` uses `base_prefactor = oft_domain_prefactor * gamma_norm_factor`. The trajectory's `_precompute_R` uses the same `base_prefactor`. Then R_a is scaled by `1/p_jump`.

**How to avoid:** Trace the exact numerical path for a single step and verify that `K0 * rho * K0' + rho_jump + U_res * rho * U_res'` matches the current `_finalize_kraus_step!` output to floating-point precision. The existing `test_cptp.jl` CPTP completeness test verifies `K0'K0 + delta*R + U'U = I` which is the algebraic identity that must hold.

**Warning signs:** Trace distances differ by more than 1e-12 compared to baseline (not just 1e-14 eigendecomposition noise).

### Pitfall 2: BohrDomain _precompute_R Requires Careful Extraction

**What goes wrong:** The BohrDomain `_jump_contribution!` for Thermalize interleaves R accumulation and rho_jump accumulation in the same inner loop over Bohr frequencies. Extracting just the R part requires understanding which operations contribute to R vs rho_jump.

**Why it happens:** In the BohrDomain inner loop (jump_workers.jl lines 226-282):
- `scratch.jump_oft = alpha(bohr_freqs, nu_2) * jump.in_eigenbasis` -- this is B_{nu_2}
- `scratch.sandwich_tmp` accumulates `rho * A_{nu_2}^dagger` -- this is rho-dependent (for rho_jump)
- `scratch.rho_jump += scaled_delta * B_{nu_2} * (rho * A_{nu_2}^dagger)` -- rho-dependent
- `scratch.R += jump_weight_scaling * A_{nu_2}^dagger * B_{nu_2}` -- rho-independent!

So the R accumulation is: `R^a = sum_{nu_2} A_{nu_2}^dagger * B_{nu_2}^a` scaled by `gamma_norm_factor / jump_prob`.

**How to avoid:** Write a dedicated `_precompute_R` for BohrDomain that mirrors the Krylov's `_accumulate_R_total!` for BohrDomain (krylov_workspace.jl lines 188-216), but operates per-jump rather than summing across all jumps.

**Warning signs:** BohrDomain trace distances drift from baseline by more than eigendecomposition noise.

### Pitfall 3: Eigendecomposition FP Variation

**What goes wrong:** When R^a is precomputed once (instead of recomputed identically each step), the eigendecomposition in `_build_cptp_channel` happens once. The previous code recomputed an IDENTICAL R each step and got an identical eigen() result. The precomputed path does eigen() once, so the result is used for all steps. Any FP noise from the single eigendecomposition propagates consistently rather than being "re-randomized" each step.

**Why it happens:** `eigen(Hermitian(S))` is deterministic for identical input. Since R was recomputed identically each step, the old code got the same eigen() result each time. The new code also gets one deterministic result. The results SHOULD match exactly.

**How to avoid:** The key guarantee is that the R computed by `_precompute_R` is bit-identical to the R computed inside `_jump_contribution!`. If the same accumulation loop and `hermitianize!` are used, results should be identical. Verify with `@test isapprox(R_precomputed, R_from_jump_contribution; atol=0)`.

**Warning signs:** If R differs at all between precomputed and per-step paths, there is a bug.

### Pitfall 4: rho_jump Still Needs the omega Loop

**What goes wrong:** After precomputing R/K0/U_residual, someone might think the entire omega loop can be eliminated. It cannot -- rho_jump depends on rho and must be recomputed each step.

**Why it happens:** `rho_jump = delta * sum_w rate^2(w) * L_{a,w} * rho * L_{a,w}^dagger` contains rho in every term.

**How to avoid:** The refactored hot loop MUST retain the omega-loop for rho_jump accumulation. Only R accumulation and `_finalize_kraus_step!` eigendecomposition are eliminated.

## Code Examples

Verified patterns from the existing codebase:

### Existing _precompute_R for EnergyDomain (trajectories.jl)

```julia
# Source: trajectories.jl lines 193-246
function _precompute_R(
    jumps::AbstractVector{<:JumpOp},
    hamiltonian::HamHam,
    config::Config{Thermalize, EnergyDomain},
    precomputed_data,
    scratch::ThermalizeScratch{<:Complex},
)
    dim = size(hamiltonian.data, 1)
    (; transition, gamma_norm_factor, energy_labels) = precomputed_data
    base_prefactor = precomputed_data.oft_domain_prefactor * gamma_norm_factor
    inv_4sigma2 = 1.0 / (4 * config.sigma^2)
    fill!(scratch.R, 0)

    @inbounds for jump in jumps
        if jump.hermitian
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)
                oft!(scratch.jump_oft, jump.in_eigenbasis, hamiltonian.bohr_freqs, w, inv_4sigma2)
                rate2_pos = base_prefactor * transition(w)
                mul!(scratch.LdagL, scratch.jump_oft', scratch.jump_oft)
                @. scratch.R += rate2_pos * scratch.LdagL
                if w > 1e-12
                    rate2_neg = base_prefactor * transition(-w)
                    mul!(scratch.LdagL, scratch.jump_oft, scratch.jump_oft')
                    @. scratch.R += rate2_neg * scratch.LdagL
                end
            end
        else
            for w in energy_labels
                oft!(scratch.jump_oft, jump.in_eigenbasis, hamiltonian.bohr_freqs, w, inv_4sigma2)
                rate2 = base_prefactor * transition(w)
                mul!(scratch.LdagL, scratch.jump_oft', scratch.jump_oft)
                @. scratch.R += rate2 * scratch.LdagL
            end
        end
    end
    hermitianize!(scratch.R)
    return scratch.R
end
```

### Existing _build_cptp_channel (furnace_utensils.jl)

```julia
# Source: furnace_utensils.jl lines 183-201
function _build_cptp_channel(R::Matrix{T}, delta::Real) where {T<:Complex}
    dim = size(R, 1)
    alpha = 1 - sqrt(1 - delta)
    K0 = Matrix{T}(I, dim, dim) .- alpha .* R
    R2 = R * R
    S = (2 * alpha - delta) .* R .- (alpha^2) .* R2
    hermitianize!(S)
    eig = eigen(Hermitian(S))           # <-- THIS EIGEN() IS ELIMINATED BY PRECOMPUTATION
    eig.values .= max.(eig.values, 0.0)
    U_residual = Matrix{T}(Diagonal(sqrt.(eig.values)) * eig.vectors')
    return (; K0, U_residual, alpha)
end
```

### Existing Trajectory Workspace Build Pattern (trajectories.jl)

```julia
# Source: trajectories.jl lines 105-122
# Per-operator Kraus data precomputation (reference for DM precomputation):
Rs = Vector{Matrix{CT}}(undef, n_jumps)
K0s = Vector{Matrix{CT}}(undef, n_jumps)
U_residuals = Vector{Matrix{CT}}(undef, n_jumps)
builder_scratch = ThermalizeScratch(CT, dim)

for a in 1:n_jumps
    _precompute_R([jumps_for_diss[a]], ham_or_trott, config, precomputed_data, builder_scratch)
    R_a = copy(builder_scratch.R)
    R_a .*= (1.0 / p_jump)
    (; K0, U_residual) = _build_cptp_channel(R_a, delta)
    Rs[a] = R_a
    K0s[a] = K0
    U_residuals[a] = U_residual
end
```

### New BohrDomain _precompute_R (Needs to be Written)

Based on the Krylov `_accumulate_R_total!` for BohrDomain (krylov_workspace.jl lines 188-216), adapted for per-jump:

```julia
# Reference pattern from krylov_workspace.jl, adapted for single-jump R:
function _precompute_R(
    jumps::AbstractVector{<:JumpOp},
    hamiltonian::HamHam,
    config::Config{Thermalize, BohrDomain},
    precomputed_data,
    scratch::ThermalizeScratch{<:Complex},
)
    dim = size(hamiltonian.data, 1)
    (; alpha, gamma_norm_factor) = precomputed_data

    fill!(scratch.R, 0)

    for jump in jumps
        for nu_2 in keys(hamiltonian.bohr_dict)
            # B_{nu_2} = sum_{nu_1} alpha(nu_1, nu_2) * A^a
            @. scratch.jump_oft = alpha(hamiltonian.bohr_freqs, nu_2) * jump.in_eigenbasis

            # R += gamma_norm_factor * A_{nu_2}^dagger * B_{nu_2}
            indices = hamiltonian.bohr_dict[nu_2]
            @inbounds for idx in indices
                i = idx[1]; j = idx[2]
                v = conj(jump.in_eigenbasis[i, j]) * gamma_norm_factor
                @inbounds for q in 1:dim
                    scratch.R[j, q] += v * scratch.jump_oft[i, q]
                end
            end
        end
    end

    hermitianize!(scratch.R)
    return scratch.R
end
```

### Refactored rho_jump-only Accumulation (New)

After extracting R precomputation, the per-step omega loop only accumulates rho_jump:

```julia
# For EnergyDomain (pseudocode for the rho_jump-only inner loop):
fill!(scratch.rho_jump, 0)
@inbounds for w_raw in energy_labels
    # ... same omega loop but ONLY the rho_jump accumulation lines:
    oft!(scratch.jump_oft, jump.in_eigenbasis, hamiltonian.bohr_freqs, w, inv_4sigma2)
    rate2_pos = base_prefactor * transition(w)
    mul!(scratch.sandwich_tmp, evolving_dm, scratch.jump_oft')  # rho * Aw^dagger
    mul!(scratch.rho_jump, scratch.jump_oft, scratch.sandwich_tmp, config.delta * rate2_pos, 1.0)
    # ... mirror for hermitian jumps
end
# Then apply precomputed channel:
# rho_next = K0 * rho * K0' + rho_jump + U_res * rho * U_res'
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `_build_cptp_channel` per step | Per-jump precomputation at startup | This phase (39) | Eliminates `n_steps * n_jumps` eigendecompositions |
| R + rho_jump interleaved in omega loop | R precomputed, only rho_jump in omega loop | This phase (39) | Removes redundant L'L accumulation from hot path |

**Already done (no change needed):**
- Coherent unitaries: Already precomputed in `run_thermalize` via `_precompute_coherent_unitary`
- `_precompute_R` for Energy/Time/Trotter: Already exists in `trajectories.jl`
- `_build_cptp_channel`: Already exists in `furnace_utensils.jl`

**Newly needed:**
- `_precompute_R` for BohrDomain: Must be written (no existing function)
- Refactored `_jump_contribution!` or new `_accumulate_rho_jump!` for all 4 domains: Extracts rho_jump-only omega loop
- Precomputation orchestration in `run_thermalize`: Calls `_precompute_R` + `_build_cptp_channel` per jump
- Updated channel application: Uses precomputed K0/U_residual instead of `_finalize_kraus_step!`

## Memory Analysis

Precomputing K0 and U_residual per jump adds memory proportional to `n_jumps * dim^2`:

| n_qubits | dim | n_jumps (3*n) | Memory per matrix (ComplexF64) | Total added (K0+U_res) | System RAM context |
|----------|-----|---------------|-------------------------------|------------------------|-------------------|
| 4 | 16 | 12 | 2 KB | 48 KB | Negligible |
| 6 | 64 | 18 | 32 KB | 1.2 MB | Negligible |
| 8 | 256 | 24 | 512 KB | 24 MB | Small |
| 10 | 1024 | 30 | 8 MB | 480 MB | Moderate |
| 12 | 4096 | 36 | 128 MB | 9.2 GB | Significant |

At n=10-12, memory may be a concern. However, the density matrix itself is already 8 MB (n=10) or 128 MB (n=12), and the NUFFT prefactors for Time/Trotter are O(dim^2 * n_energies) which is already much larger. The precomputed K0/U_residual adds a constant factor per jump, not per energy label.

**Mitigation for large systems:** The precomputed matrices have the same dimension as the density matrix -- if the density matrix fits in memory, the precomputed channels will too.

## Implementation Strategy

### Phase Decomposition

**Plan 39-01: Precomputation Infrastructure**
1. Write `_precompute_R` for BohrDomain (new function in `trajectories.jl`)
2. Create `_precompute_per_jump_channels` helper function that:
   - Creates temporary ThermalizeScratch for building R
   - Calls `_precompute_R` per jump
   - Scales R by `1/p_jump`
   - Calls `_build_cptp_channel` per jump
   - Returns `(K0s, U_residuals)` vectors
3. Write `_accumulate_rho_jump!` functions for Energy, Time/Trotter, and Bohr domains (extract rho_jump-only logic from `_jump_contribution!`)
4. Write `_apply_precomputed_channel!` (replaces `_finalize_kraus_step!` for precomputed case)

**Plan 39-02: Integration and Validation**
1. Refactor `run_thermalize` to call precomputation before loop
2. Replace per-step `_jump_contribution!` with `_accumulate_rho_jump!` + `_apply_precomputed_channel!`
3. Validate numerical equivalence (trace distance to Gibbs matches pre-change baselines)
4. Verify no `eigen()` or `_build_cptp_channel` calls in hot loop
5. Update/regenerate regression baselines if needed (expect O(1e-14) shifts)
6. Run full test suite

### DM vs Trajectory Scaling Convention

The existing `_precompute_R` functions already use the correct `base_prefactor = oft_domain_prefactor * gamma_norm_factor`. The trajectory engine then scales by `1/p_jump`. The DM `_jump_contribution!` uses `jump_weight_scaling = rescale_by_inv_prob ? (gamma_norm_factor / jump_prob) : gamma_norm_factor`.

**Critical detail:** In the DM code, the R-accumulation weights are `base_prefactor = oft_domain_prefactor * jump_weight_scaling` where `jump_weight_scaling = gamma_norm_factor / jump_prob` (when `rescale_by_inv_prob=true`, which is the default). So the total scaling is `oft_domain_prefactor * gamma_norm_factor / p_jump`.

In `_precompute_R`, the weights are `oft_domain_prefactor * gamma_norm_factor`, then the trajectory engine separately multiplies by `1/p_jump`. This matches.

For rho_jump, the DM code uses `scaled_delta = config.delta * jump_weight_scaling = config.delta * gamma_norm_factor / p_jump`. In the precomputed path, the rho_jump accumulation should use the same `scaled_delta`.

The `_build_cptp_channel` receives R already scaled by `1/p_jump`, and uses bare `config.delta`. This is consistent between DM and trajectory:
- `alpha = 1 - sqrt(1 - delta)` (bare delta)
- `K0 = I - alpha * R` (R already scaled by 1/p_jump)
- `S = (2*alpha - delta) * R - alpha^2 * R^2` (bare delta, scaled R)

### BohrDomain: Precomputation Is Still Beneficial

PRECOMP-04 says "BohrDomain receives general speedups and threading but no per-Bohr-frequency vector precomputation." This means we should NOT precompute per-Bohr-frequency data. But we CAN and SHOULD precompute R^a (which is a single dim x dim matrix per jump), K0^a, and U_residual^a. The R^a computation iterates over Bohr frequencies, but the result is a single matrix. This eliminates the per-step eigendecomposition for BohrDomain too.

The "Bohr frequency count grows too fast" concern is about storing per-frequency data (O(dim^2) unique frequencies, each with associated matrices). Precomputing R^a just sums over those frequencies into one matrix -- same computation as today, just done once instead of every step.

## Open Questions

1. **Should `_precompute_R` move to a shared file?**
   - What we know: Currently lives in `trajectories.jl` but will be used by `run_thermalize` too.
   - What's unclear: Best file organization -- `furnace_utensils.jl` (where `_build_cptp_channel` lives) or keep in `trajectories.jl`?
   - Recommendation: Move `_precompute_R` functions to `furnace_utensils.jl` alongside `_build_cptp_channel`, since both are now shared between DM and trajectory paths.

2. **Should the refactored `_jump_contribution!` replace the existing one or be a new function?**
   - What we know: The existing `_jump_contribution!` signatures serve both Lindbladian (vectorized Liouvillian) and Thermalize (DM evolution) paths. Only the Thermalize path benefits from precomputation.
   - What's unclear: Whether to add new dispatch signatures or modify existing ones.
   - Recommendation: Create new `_accumulate_rho_jump!` functions (simpler than modifying `_jump_contribution!`), keeping the existing `_jump_contribution!` intact for the Lindbladian path and as fallback.

3. **How to handle `rescale_by_inv_prob=false` case?**
   - What we know: The default is `true` (delta rescaled by 1/p_jump). When `false`, `jump_weight_scaling = gamma_norm_factor` (no 1/p_jump factor).
   - What's unclear: Whether anyone uses `rescale_by_inv_prob=false` in practice.
   - Recommendation: Precompute for the default `rescale_by_inv_prob=true` case. If `false`, could either precompute differently or fall back to non-precomputed path. Simplest: precompute R without the 1/p_jump scaling, then apply the scaling in `_build_cptp_channel` based on the flag.

## Sources

### Primary (HIGH confidence)
- `src/furnace.jl` -- `run_thermalize` implementation (lines 143-223)
- `src/jump_workers.jl` -- `_jump_contribution!` for all domains (Thermalize variants at lines 194-440)
- `src/furnace_utensils.jl` -- `_build_cptp_channel` (lines 183-201), `_precompute_data` (lines 30-141)
- `src/trajectories.jl` -- `_build_trajectory_workspace` (lines 61-152), `_precompute_R` (lines 193-300)
- `src/krylov_workspace.jl` -- `_accumulate_R_total!` for BohrDomain (lines 188-216)
- `src/coherent.jl` -- `_precompute_coherent_unitary` (lines 58-100)
- `src/structs.jl` -- `ThermalizeScratch`, `Workspace` definitions
- `test/test_cptp.jl` -- CPTP completeness verification
- `test/test_regression.jl` -- Regression test baselines

### Secondary (MEDIUM confidence)
- `.planning/ROADMAP.md` -- Phase 39 requirements and success criteria
- `.planning/REQUIREMENTS.md` -- PRECOMP-01 through PRECOMP-04 specifications
- `.planning/PROJECT.md` -- Architectural context and key decisions

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- no new libraries, pure refactoring of existing validated code
- Architecture: HIGH -- trajectory engine already implements the exact same pattern; DM path is a direct port
- Pitfalls: HIGH -- all identified from direct code reading; scaling conventions traced through actual code paths
- BohrDomain _precompute_R: MEDIUM -- new function needed, but pattern clearly derivable from existing `_accumulate_R_total!` and current `_jump_contribution!` Bohr Thermalize code

**Research date:** 2026-03-01
**Valid until:** 2026-04-01 (stable -- internal refactoring, no external dependencies changing)
