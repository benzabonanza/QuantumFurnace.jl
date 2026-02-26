# Phase 34: Code Deduplication - Research

**Researched:** 2026-02-26
**Domain:** Julia multiple dispatch refactoring, cross-simulation-type code deduplication in quantum simulation codebase
**Confidence:** HIGH

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

#### domain_prefactor
- **Name:** `domain_prefactor` (purely scalar constant, no NUFFT/OFT involvement)
- **Dispatch:** Domain-only, no construction type needed
  - `domain_prefactor(::EnergyDomain, w0, sigma)` -> `w0 / (sigma * sqrt(2*pi))`
  - `domain_prefactor(::TimeDomain, w0, sigma, t0)` / `domain_prefactor(::TrotterDomain, w0, sigma, t0)` -> `w0 * t0^2 * sigma * sqrt(2/pi) / (2*pi)`
  - BohrDomain: **no method** -- callers use `gamma_norm_factor` directly, don't call `domain_prefactor`
- **Arguments:** Explicit fields (w0, sigma, t0), NOT Config
- **Called in:** `_precompute_data`, result stored in `precomputed_data` NamedTuple
- **Callers multiply:** `precomputed_data.domain_prefactor * gamma_norm_factor * [their own scaling]`

#### Frequency iteration (SKIPPED)
- **Decision:** Keep the hermitian half-grid for-loop pattern as-is in each call site
- **Rationale:** User prefers the explicit for-loop structure; deduplication effort should target higher-level functions, not micro-patterns

#### OFT unification
- **Merge:** `oft!` + `_krylov_oft!` -> single `oft!` with optimized signature
- **Signature:** `oft!(out, eigenbasis::Matrix, bohr_freqs::Matrix, energy, inv_4sigma2)` -- concrete types, precomputed `inv_4sigma2`
- **Old JumpOp-based `oft!` signature:** Removed
- **Scope:** EnergyDomain only. Time/TrotterDomain production path uses NUFFT prefactors as-is
- **`inv_4sigma2`:** Computed at call site (`1/(4*sigma^2)`), NOT stored in precomputed_data
- **Test utilities:** `time_oft!` and `trotter_oft!` remain as clearly-marked test/debug utilities

#### Cross-simulation-type deduplication (NEW -- identified during discussion)

Three shared helpers to extract into `furnace_utensils.jl`:

1. **R accumulation** (`_accumulate_R!`)
   - Same formula in 3 places: `R += sum_w rate^2(w) * (L'L)`
   - Thermalize: per-operator R^a (per step, hot path)
   - Krylov: summed R_total across all operators (once at construction)
   - Trajectories: per-operator R^a (once at framework build)
   - Difference: Thermalize per-step vs others precomputed; Trajectories rescales by n_jumps
   - Domain variants stay inline (Bohr/Energy/Time each have their own OFT)

2. **CPTP channel construction** (`_build_kraus_matrices!`)
   - Same formula in 3 places: `K0 = I - alpha*R`, `S = (2*alpha-delta)*R - alpha^2*R^2`, `U_res = sqrt_psd(S)`
   - Thermalize: jump_workers.jl
   - Krylov: krylov_workspace.jl
   - Trajectories: trajectories.jl

3. **Physics-convention sandwich** (shared helper, after convention bug fix)
   - After fixing krylov_matvec to use `L * rho * L'` everywhere
   - Thermalize and Krylov eigsolve already use physics convention
   - Extract shared `_accumulate_sandwich!(out, L, rho, scalar, ws)` using `L * rho * L'`

#### Sandwich convention bug (prerequisite quick task)
- krylov_matvec.jl `_accumulate_sandwich!` NOW uses `L * rho * L'` (fixed in quick-39)
- krylov_eigsolve.jl `_accumulate_jump_sandwich!` and Thermalize already used `L * rho * L'`
- **Status:** FIXED in commits 9c49f40..b4b1cec (quick-39 branch)

### Claude's Discretion
- Exact internal structure of `_accumulate_R!` (how to parameterize per-operator vs all-operator accumulation)
- Whether sandwich helper needs separate forward/adjoint variants or can be parameterized
- Order of extraction (prefactor -> OFT -> R -> CPTP -> sandwich, or different sequencing)

### Deferred Ideas (OUT OF SCOPE)
- Sandwich convention consolidation for BohrDomain 2-operator sandwich (`_accumulate_sandwich_2op!` / `_accumulate_adjoint_sandwich_2op!`) -- different enough to warrant separate treatment
- `_run_chunk_*` variant consolidation (3 variants differing only in measurement flags) -- minor, can be done in Phase 38 test cleanup
- `B_time`/`B_trotter` single/vector deduplication in coherent.jl -- low priority
- `_pick_transition_kms`/`_pick_transition_gns` consolidation -- belongs in Phase 35 or later when construction type dispatch is more mature
</user_constraints>

## Summary

Phase 34 targets five concrete deduplication tasks across three simulation types (Thermalize, Krylov/Lindbladian, Trajectory). The codebase contains 16 copy-pasted prefactor formula instances across 5 source files, 3 near-identical CPTP channel construction blocks, 3 R-accumulation implementations (2 domain variants each), and 2 functionally identical OFT implementations (`oft!` vs `_krylov_oft!`). In addition, the sandwich helpers in krylov_matvec.jl now have 4 functions that reduce to only 2 distinct computations (`L*rho*L'` and `L'*rho*L`), which can be consolidated.

The frequency iteration half-grid pattern is explicitly **skipped** per user decision. BohrDomain 2-operator sandwiches are explicitly **deferred**.

**Primary recommendation:** Extract helpers in dependency order: (1) `domain_prefactor` -> (2) unified `oft!` -> (3) `_accumulate_R!` -> (4) `_build_kraus_matrices!` -> (5) consolidate sandwich helpers. Each step enables the next while keeping existing tests passing throughout.

## Standard Stack

Not applicable -- this phase is pure refactoring within an existing Julia project. No new libraries needed.

### Relevant Julia Patterns
| Pattern | Purpose | Why Standard |
|---------|---------|--------------|
| Multiple dispatch on domain singletons | `domain_prefactor(::EnergyDomain, ...)` | Already used throughout codebase (Config{S,D,C,T}) |
| `@inline` annotation on hot-path helpers | Avoid call overhead in inner loops | Used on existing `_krylov_oft!`, `_accumulate_sandwich!` |
| Named tuples for precomputed data | Zero-overhead struct-like containers | Already used for `precomputed_data` throughout |

## Architecture Patterns

### Pattern 1: domain_prefactor as a Pure Function

**What:** A domain-dispatched function returning a scalar constant.
**When to use:** Every place that currently computes `w0 / (sigma * sqrt(2*pi))` (EnergyDomain) or `w0 * t0^2 * sigma * sqrt(2/pi) / (2*pi)` (Time/TrotterDomain).

**Current locations (16 instances):**

| File | Line(s) | Domain | Context |
|------|---------|--------|---------|
| `jump_workers.jl` | 65 | Energy | Lindbladian `_jump_contribution!` |
| `jump_workers.jl` | 109 | Time/Trotter | Lindbladian `_jump_contribution!` |
| `jump_workers.jl` | 334 | Energy | Thermalize `_jump_contribution!` |
| `jump_workers.jl` | 410 | Time/Trotter | Thermalize `_jump_contribution!` |
| `krylov_workspace.jl` | 189 | Energy | `_accumulate_R_total!` |
| `krylov_workspace.jl` | 234 | Time/Trotter | `_accumulate_R_total!` |
| `krylov_matvec.jl` | 166 | Energy | `apply_lindbladian!` |
| `krylov_matvec.jl` | 234 | Energy | `apply_adjoint_lindbladian!` |
| `krylov_matvec.jl` | 477 | Time/Trotter | `apply_lindbladian!` |
| `krylov_matvec.jl` | 547 | Time/Trotter | `apply_adjoint_lindbladian!` |
| `krylov_eigsolve.jl` | 231 | Energy | `_accumulate_jump_sandwich!` |
| `krylov_eigsolve.jl` | 277 | Time/Trotter | `_accumulate_jump_sandwich!` |
| `trajectories.jl` | 186 | Energy | `build_trajectoryframework` |
| `trajectories.jl` | 189 | Time/Trotter | `build_trajectoryframework` |
| `trajectories.jl` | 245 | Energy | `_precompute_R` |
| `trajectories.jl` | 301 | Time/Trotter | `_precompute_R` |

**Implementation:**
```julia
# In furnace_utensils.jl
domain_prefactor(::EnergyDomain, w0, sigma) = w0 / (sigma * sqrt(2 * pi))
domain_prefactor(::TimeDomain, w0, sigma, t0) = w0 * t0^2 * (sigma * sqrt(2 / pi)) / (2 * pi)
domain_prefactor(::TrotterDomain, w0, sigma, t0) = w0 * t0^2 * (sigma * sqrt(2 / pi)) / (2 * pi)
# No method for BohrDomain (callers use gamma_norm_factor directly)
```

**Call in `_precompute_data`:**
```julia
# Store in precomputed_data NamedTuple for zero-overhead access:
function _precompute_data(config::Config{<:Any, EnergyDomain}, ...)
    ...
    dp = domain_prefactor(config.domain, config.w0, config.sigma)
    return (
        ...,
        domain_prefactor = dp,
    )
end
```

**Callers then use:**
```julia
prefactor = precomputed_data.domain_prefactor * gamma_norm_factor
# instead of:
# prefactor = config.w0 / (config.sigma * sqrt(2 * pi)) * gamma_norm_factor
```

**Important nuance:** Some call sites multiply by additional factors:
- `jump_workers.jl` Thermalize variants multiply by `jump_weight_scaling` (not `gamma_norm_factor`)
- `trajectories.jl` `build_trajectoryframework` additionally divides by `(1.0 / n_jumps)` i.e. multiplies by `n_jumps`
- These must be preserved as caller-side multiplications.

### Pattern 2: Unified oft! Signature

**What:** Merge `oft!` (ofts.jl) and `_krylov_oft!` (krylov_matvec.jl) into one function.
**Current state:**

```julia
# ofts.jl -- current (uses JumpOp, no inv_4sigma2 precomputation)
function oft!(out_matrix, jump::JumpOp, energy, hamiltonian::HamHam, sigma)
    @. out_matrix = jump.in_eigenbasis * exp(-(energy - hamiltonian.bohr_freqs)^2 / (4 * sigma^2))
end

# krylov_matvec.jl -- current (uses concrete types, precomputed inv_4sigma2)
@inline function _krylov_oft!(out, eigenbasis::Matrix{T}, bohr_freqs::Matrix{<:Real}, energy, inv_4sigma2) where {T<:Complex}
    @. out = eigenbasis * exp(-(energy - bohr_freqs)^2 * inv_4sigma2)
end
```

**Target:** Single `oft!` with the optimized signature:
```julia
@inline function oft!(out::Matrix{T}, eigenbasis::Matrix{T}, bohr_freqs::Matrix{<:Real}, energy::Real, inv_4sigma2::Real) where {T<:Complex}
    @. out = eigenbasis * exp(-(energy - bohr_freqs)^2 * inv_4sigma2)
    return nothing
end
```

**Callers to update:**
- `jump_workers.jl` EnergyDomain Lindbladian and Thermalize (4 call sites) -- compute `inv_4sigma2 = 1/(4*sigma^2)` at call site, pass `jump.in_eigenbasis`, `hamiltonian.bohr_freqs`
- `trajectories.jl` EnergyDomain `_precompute_R` and `step_along_trajectory!` (3 call sites)
- `krylov_workspace.jl` EnergyDomain `_accumulate_R_total!` (already uses `_krylov_oft!` pattern)
- `krylov_matvec.jl` (6 call sites, already uses `_krylov_oft!`)
- `krylov_eigsolve.jl` (3 call sites, already uses `_krylov_oft!`)

**Callers that keep old interface (test/debug):** `time_oft!` and `trotter_oft!` in ofts.jl remain unchanged (they use OFTCaches, distinct purpose).

### Pattern 3: _accumulate_R! Shared Helper

**What:** Extract the R-accumulation loop (sum of `rate^2 * L'L` over frequencies) into a shared function.
**Current 3 implementations:**
1. `krylov_workspace.jl` `_accumulate_R_total!` -- sums over ALL jumps, used at workspace construction
2. `trajectories.jl` `_precompute_R` -- operates on a subset of jumps (typically one), returns `scratch.R`
3. `jump_workers.jl` `_jump_contribution!` Thermalize variants -- inline R accumulation per step

**Key differences:**
| Aspect | krylov_workspace | trajectories | jump_workers (Thermalize) |
|--------|-----------------|--------------|---------------------------|
| Jump iteration | All jumps | Subset (typically 1) | Single jump |
| Buffer | Dedicated `R` matrix + local `jump_oft`/`LdagL` | `scratch.R`, `scratch.jump_oft`, `scratch.LdagL` | `scratch.R` (same KrausScratch) |
| Domain dispatch | EnergyDomain / Time-TrotterDomain / BohrDomain | EnergyDomain / Time-TrotterDomain (no Bohr) | All 3 domains |
| Scaling | `* gamma_norm_factor` | `* gamma_norm_factor` | `* jump_weight_scaling` (may differ from gnf) |
| Called when | Once at construction | Once per operator at framework build | Every step (hot path) |

**Recommendation for discretion area:**

Extract a core helper for EnergyDomain and Time/TrotterDomain only (BohrDomain has fundamentally different iteration structure with bohr_dict):

```julia
function _accumulate_R_energy!(
    R::Matrix{T},
    eigenbases::AbstractVector{<:Matrix{T}},
    hermitian_flags::AbstractVector{Bool},
    energy_labels::Vector{Float64},
    transition,
    prefactor::Real,
    bohr_freqs::Matrix{<:Real},
    inv_4sigma2::Real,
    jump_oft::Matrix{T},
    LdagL::Matrix{T},
) where {T<:Complex}
```

And a NUFFT variant:
```julia
function _accumulate_R_nufft!(
    R::Matrix{T},
    eigenbases::AbstractVector{<:Matrix{T}},
    hermitian_flags::AbstractVector{Bool},
    energy_labels::Vector{Float64},
    transition,
    prefactor::Real,
    oft_nufft_prefactors,
    jump_oft::Matrix{T},
    LdagL::Matrix{T},
) where {T<:Complex}
```

Both accept a vector of eigenbases (length 1 for per-operator, length N for all-operator). The prefactor already includes any scaling (gamma_norm_factor or jump_weight_scaling) -- callers compose.

The `jump_workers.jl` Thermalize hot path also accumulates `rho_jump` interleaved with R. This means the Thermalize per-step code CANNOT use the R-only helper directly -- it must keep its interleaved loop. The helper extracts the R-only precomputation case (krylov_workspace and trajectories).

### Pattern 4: _build_kraus_matrices! Shared CPTP Channel Construction

**What:** Extract the common K0/S/U_residual construction from the R matrix.
**Current 3 implementations:**
1. `jump_workers.jl` `_finalize_kraus_step!` (lines 176-211) -- operates on `scratch.R`
2. `krylov_workspace.jl` Config{Thermalize} constructor (lines 386-397) -- operates on `R_total`
3. `trajectories.jl` `build_trajectoryframework` (lines 159-177) -- operates on per-operator `R_a`

**Common formula (Chen Eq. 3.2):**
```
alpha = 1 - sqrt(1 - delta)
K0 = I - alpha * R
S = (2*alpha - delta)*R - alpha^2 * R^2
U_residual = sqrt_psd(S)   # via eigen + clamp negatives + sqrt
```

**Key differences:**
- `_finalize_kraus_step!` also applies the channel to `evolving_dm` in the same function
- `krylov_workspace` and `trajectories` only precompute the matrices (K0, U_residual) and store them
- All three use the same formula for K0, S, U_residual

**Recommended extraction:**
```julia
function _build_cptp_channel(R::Matrix{T}, delta::Real) where {T<:Complex}
    dim = size(R, 1)
    alpha = 1 - sqrt(1 - delta)
    K0 = Matrix{T}(I, dim, dim) .- alpha .* R
    R2 = R * R
    S = (2 * alpha - delta) .* R .- (alpha^2) .* R2
    hermitianize!(S)
    eig = eigen(Hermitian(S))
    eig.values .= max.(eig.values, 0.0)
    U_residual = Matrix{T}(Diagonal(sqrt.(eig.values)) * eig.vectors')
    return (; K0, U_residual, alpha)
end
```

Then:
- `krylov_workspace.jl` calls `(; K0, U_residual, alpha) = _build_cptp_channel(R_total, delta)`
- `trajectories.jl` calls `(; K0, U_residual, _) = _build_cptp_channel(R_a, delta)` per operator
- `_finalize_kraus_step!` can call `(; K0, U_residual, _) = _build_cptp_channel(scratch.R, delta)` then apply the channel

**Concern:** `_finalize_kraus_step!` is called every step in the Thermalize hot path. Moving the eigendecomposition into a separate function call should not affect performance (Julia inlines small functions). However, the function currently uses `scratch.K0`, `scratch.tmp1`, `scratch.tmp2`, `scratch.LdagL` buffers. The extracted function would allocate K0 and U_residual -- but `_finalize_kraus_step!` already allocates via `eigen()`, so the channel construction is not on the zero-allocation path. This is acceptable.

### Pattern 5: Sandwich Helper Consolidation

**What:** After quick-39 fix, the 4 sandwich functions in `krylov_matvec.jl` reduce to 2 distinct computations:
- `_accumulate_sandwich!` and `_accumulate_adjoint_sandwich_adj_L!` both compute `L * rho * L'`
- `_accumulate_sandwich_adj_L!` and `_accumulate_adjoint_sandwich!` both compute `L' * rho * L`

**Post-fix state (from krylov_matvec.jl):**
```julia
# _accumulate_sandwich!:              L * rho * L'  (gemm 'N','N' then 'N','C')
# _accumulate_adjoint_sandwich_adj_L!: L * rho * L'  (gemm 'N','N' then 'N','C')  <- IDENTICAL
# _accumulate_sandwich_adj_L!:         L' * rho * L  (gemm 'C','N' then 'N','N')
# _accumulate_adjoint_sandwich!:       L' * rho * L  (gemm 'C','N' then 'N','N')  <- IDENTICAL
```

**Recommendation:** Consolidate to 2 functions:
```julia
_accumulate_sandwich!(out, L, rho, scalar, ws)        # L * rho * L'
_accumulate_sandwich_adj!(out, L, rho, scalar, ws)     # L' * rho * L
```

All 4 current callers redirect to the appropriate one. No separate forward/adjoint variant needed because the adjoint of `L*rho*L'` is `L'*rho*L` -- the same computation as the negative-frequency partner.

**BohrDomain 2-operator sandwiches** (`_accumulate_sandwich_2op!`, `_accumulate_adjoint_sandwich_2op!`) are explicitly deferred per user decision.

### Recommended Execution Order

1. **domain_prefactor** -- pure addition, no callers change behavior. Add function, compute in `_precompute_data`, store in NamedTuple. Update all 16 call sites.
2. **oft! unification** -- replace old `oft!` signature with new concrete-typed one. Update callers in jump_workers.jl and trajectories.jl to pass `inv_4sigma2`. Delete `_krylov_oft!`.
3. **Sandwich consolidation** -- merge 4 helpers into 2 in krylov_matvec.jl. Update caller sites.
4. **_build_cptp_channel** -- extract CPTP construction, used by 3 files. Update krylov_workspace.jl, trajectories.jl, jump_workers.jl.
5. **_accumulate_R!** -- extract R-accumulation helpers. Unify krylov_workspace.jl `_accumulate_R_total!` and trajectories.jl `_precompute_R`.

This order minimizes risk: each step is independently testable, and later steps build on earlier ones (e.g., `_accumulate_R!` uses the unified `oft!`).

### Anti-Patterns to Avoid

- **Over-abstracting BohrDomain:** BohrDomain R-accumulation iterates over `bohr_dict` keys with index scatter, fundamentally different from the energy-grid loop. Do not force it into the same helper.
- **Breaking allocation invariants:** The Krylov hot path (apply_lindbladian!, step_along_trajectory!) must remain zero-allocation. Extracted helpers must accept pre-allocated scratch buffers, never allocate internally.
- **Changing NamedTuple fields without updating all consumers:** Adding `domain_prefactor` to `precomputed_data` requires updating every destructuring `(; ...)` that consumes it.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| PSD matrix square root | Custom sqrt via SVD | `eigen(Hermitian(S))` + clamp + sqrt | Already used, handles numerical negatives correctly |
| In-place BLAS operations | Manual loops for matrix multiply | `BLAS.gemm!` / `mul!` | Already used, critical for zero-allocation hot paths |

**Key insight:** The deduplication targets are purely internal patterns -- no external libraries are relevant. The risk is in breaking numerical equivalence, not in choosing wrong tools.

## Common Pitfalls

### Pitfall 1: NamedTuple Field Addition Breaking Destructuring
**What goes wrong:** Adding `domain_prefactor` to the `precomputed_data` NamedTuple is safe (extra fields are ignored by destructuring), BUT the NamedTuple type parameter `PD` in `KrylovWorkspace{T, PD}` changes, potentially affecting precompilation.
**Why it happens:** Julia NamedTuples are parametric on their field names and types.
**How to avoid:** The type parameter `PD` is already `<:NamedTuple` so it will accommodate any shape. Just ensure tests pass after the change.
**Warning signs:** Type inference failures or unexpected recompilation.

### Pitfall 2: Thermalize Hot Path Cannot Use R-Only Helper
**What goes wrong:** Attempting to extract R accumulation from `_jump_contribution!` Thermalize variants breaks because R and rho_jump are accumulated in the same loop.
**Why it happens:** The Thermalize EnergyDomain and Time/TrotterDomain `_jump_contribution!` functions interleave `R += rate2 * L'L` with `rho_jump += delta * rate2 * L * rho * L'` in the same frequency loop. Splitting them would double the OFT computations.
**How to avoid:** Only extract R-only precomputation (krylov_workspace, trajectories). Leave Thermalize per-step code with its interleaved loop.
**Warning signs:** Duplicated OFT computation in the extracted helper path.

### Pitfall 3: Prefactor Scaling Differences Across Simulation Types
**What goes wrong:** Naively replacing all prefactor computations with `domain_prefactor * gamma_norm_factor` misses that some callers use `jump_weight_scaling` instead of `gamma_norm_factor`.
**Why it happens:** Thermalize `_jump_contribution!` computes `jump_weight_scaling = rescale_by_inv_prob ? (gamma_norm_factor / jump_prob) : gamma_norm_factor`. Trajectories multiply by `n_jumps`.
**How to avoid:** `domain_prefactor` returns ONLY the domain-dependent part (no `gamma_norm_factor`). Callers compose: `domain_prefactor * their_scaling`.
**Warning signs:** Numerical differences in regression tests.

### Pitfall 4: OFT Signature Change Breaking External Callers
**What goes wrong:** The old `oft!(out, jump::JumpOp, energy, hamiltonian, sigma)` is exported and may be used in user scripts/notebooks.
**Why it happens:** `oft!` is listed in the module exports.
**How to avoid:** Keep the old signature as a thin wrapper that computes `inv_4sigma2 = 1/(4*sigma^2)` and delegates to the new signature. Or add the new method as an overload (multiple dispatch handles both).
**Warning signs:** `MethodError` in test files or notebooks.

### Pitfall 5: Sandwich Consolidation Missing a Call Site
**What goes wrong:** After renaming/removing `_accumulate_adjoint_sandwich_adj_L!` and `_accumulate_adjoint_sandwich!`, a call site is missed.
**Why it happens:** These functions are called from 6 methods in krylov_matvec.jl (3 forward + 3 adjoint, for Energy/Time-Trotter/Bohr domains) plus krylov_eigsolve.jl sandwich helpers.
**How to avoid:** Use grep to find all call sites before renaming. Run full test suite after each consolidation step.
**Warning signs:** `MethodError: no method matching` at runtime.

## Code Examples

### Example 1: domain_prefactor Implementation and Integration

```julia
# furnace_utensils.jl -- add at top of file
domain_prefactor(::EnergyDomain, w0::Real, sigma::Real) = w0 / (sigma * sqrt(2 * pi))
domain_prefactor(::TimeDomain, w0::Real, sigma::Real, t0::Real) = w0 * t0^2 * (sigma * sqrt(2 / pi)) / (2 * pi)
domain_prefactor(::TrotterDomain, w0::Real, sigma::Real, t0::Real) = w0 * t0^2 * (sigma * sqrt(2 / pi)) / (2 * pi)

# In _precompute_data for EnergyDomain:
function _precompute_data(config::Config{<:Any, EnergyDomain}, ham_or_trott)
    energy_labels, = _precompute_labels(config)
    transition = pick_transition(config)
    gamma_norm_factor = 1.0 / maximum(transition.(energy_labels))
    dp = domain_prefactor(config.domain, config.w0, config.sigma)
    return (
        transition = transition,
        gamma_norm_factor = gamma_norm_factor,
        energy_labels = energy_labels,
        domain_prefactor = dp,
    )
end

# Caller sites change from:
#   prefactor = (config.w0 / (config.sigma * sqrt(2 * pi))) * gamma_norm_factor
# to:
#   prefactor = precomputed_data.domain_prefactor * gamma_norm_factor
```

### Example 2: Unified oft! with Backward-Compatible Wrapper

```julia
# ofts.jl -- new primary signature
@inline function oft!(
    out::Matrix{T},
    eigenbasis::Matrix{T},
    bohr_freqs::Matrix{<:Real},
    energy::Real,
    inv_4sigma2::Real,
) where {T<:Complex}
    @. out = eigenbasis * exp(-(energy - bohr_freqs)^2 * inv_4sigma2)
    return nothing
end

# Backward-compatible wrapper (for any remaining callers using JumpOp)
function oft!(out_matrix::Matrix{<:Complex}, jump::JumpOp, energy::Real, hamiltonian::HamHam, sigma::Real)
    inv_4sigma2 = 1.0 / (4 * sigma^2)
    @. out_matrix = jump.in_eigenbasis * exp(-(energy - hamiltonian.bohr_freqs)^2 * inv_4sigma2)
    return out_matrix
end
```

### Example 3: _build_cptp_channel Extraction

```julia
# furnace_utensils.jl
function _build_cptp_channel(R::Matrix{T}, delta::Real) where {T<:Complex}
    dim = size(R, 1)
    alpha = 1 - sqrt(1 - delta)

    # K0 = I - alpha * R
    K0 = Matrix{T}(I, dim, dim) .- alpha .* R

    # S = (2*alpha - delta)*R - alpha^2 * R^2
    R2 = R * R
    S = (2 * alpha - delta) .* R .- (alpha^2) .* R2
    hermitianize!(S)

    # PSD guard + sqrt
    eig = eigen(Hermitian(S))
    eig.values .= max.(eig.values, 0.0)
    U_residual = Matrix{T}(Diagonal(sqrt.(eig.values)) * eig.vectors')

    return (; K0, U_residual, alpha)
end
```

### Example 4: Consolidated Sandwich Helpers

```julia
# krylov_matvec.jl -- two functions instead of four
@inline function _accumulate_sandwich!(
    out::Matrix{T}, L_op::Matrix{T}, rho::Matrix{T}, scalar::Real, ws,
) where {T<:Complex}
    CT = one(T); ZT = zero(T)
    BLAS.gemm!('N', 'N', CT, L_op, rho, ZT, ws.tmp1)       # tmp1 = L * rho
    BLAS.gemm!('N', 'C', CT, ws.tmp1, L_op, ZT, ws.LdagL)  # LdagL = L * rho * L'
    BLAS.axpy!(T(scalar), ws.LdagL, out)
    return nothing
end

@inline function _accumulate_sandwich_adj!(
    out::Matrix{T}, L_op::Matrix{T}, rho::Matrix{T}, scalar::Real, ws,
) where {T<:Complex}
    CT = one(T); ZT = zero(T)
    BLAS.gemm!('C', 'N', CT, L_op, rho, ZT, ws.tmp1)       # tmp1 = L' * rho
    BLAS.gemm!('N', 'N', CT, ws.tmp1, L_op, ZT, ws.LdagL)  # LdagL = L' * rho * L
    BLAS.axpy!(T(scalar), ws.LdagL, out)
    return nothing
end
```

Call site mapping:
- `_accumulate_sandwich!` (forward, positive freq) -> `_accumulate_sandwich!` (unchanged)
- `_accumulate_sandwich_adj_L!` (forward, negative freq) -> `_accumulate_sandwich_adj!`
- `_accumulate_adjoint_sandwich!` (adjoint, positive freq) -> `_accumulate_sandwich_adj!`
- `_accumulate_adjoint_sandwich_adj_L!` (adjoint, negative freq) -> `_accumulate_sandwich!`

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `conj(L)*rho*L^T` sandwich (Krylov) | `L*rho*L'` sandwich (physics convention) | quick-39 (just before Phase 34) | Sandwich helpers now consistent with Thermalize/eigsolve |
| Separate `oft!` + `_krylov_oft!` | Both exist, functionally identical | Phase 27-28 introduced `_krylov_oft!` | Duplication to eliminate |
| Per-domain inline prefactor | To be extracted as `domain_prefactor()` | Phase 34 (this phase) | 16 instances to consolidate |
| 3 independent CPTP constructions | To be extracted as `_build_cptp_channel` | Phase 34 (this phase) | 3 instances to consolidate |

## Open Questions

1. **Should `_finalize_kraus_step!` call `_build_cptp_channel` internally?**
   - What we know: `_finalize_kraus_step!` runs every Thermalize step and already allocates via `eigen()`. The CPTP construction is a subset of its work.
   - What's unclear: Whether splitting the function reduces readability vs the dedup benefit.
   - Recommendation: YES, extract and call. The allocation cost is dominated by `eigen()` which stays regardless. The code clarity gain from having one CPTP formula is worth it.

2. **Should `_accumulate_R!` handle BohrDomain?**
   - What we know: BohrDomain R-accumulation iterates over `bohr_dict` keys with sparse index scatter, fundamentally different from the energy-grid pattern.
   - What's unclear: Whether a sufficiently abstract interface could cover both without harming readability.
   - Recommendation: NO. Keep BohrDomain R-accumulation inline. The code paths are structurally different enough that forcing them into one helper would be more confusing than helpful.

3. **What about `_precompute_data` NamedTuple growing with `domain_prefactor`?**
   - What we know: The NamedTuple currently varies by domain (BohrDomain has `alpha`/`gamma_norm_factor`, EnergyDomain has `transition`/`gamma_norm_factor`/`energy_labels`, etc.).
   - What's unclear: Whether adding `domain_prefactor` to only Energy/Time/Trotter precomputed_data (not Bohr) creates confusion.
   - Recommendation: Only add to domains that have it. BohrDomain has no `domain_prefactor` method, so its precomputed_data naturally lacks the field. Callers for BohrDomain never access it.

## Sources

### Primary (HIGH confidence)
- Direct codebase analysis of all 28 source files in `src/`
- Git history analysis of quick-39 convention fix (commits 709f006..b4b1cec)
- Full test suite structure (32 test files)

### Confidence Assessment
All findings are based on direct source code reading with exact line-number references. No external library research was needed.

## Metadata

**Confidence breakdown:**
- Architecture patterns: HIGH -- based on exact code analysis of all 16 prefactor instances, 3 CPTP constructions, 3 R-accumulation implementations, and 4 sandwich helpers
- Pitfalls: HIGH -- derived from actual code structure (interleaved loops, scaling differences)
- Execution order: HIGH -- dependency analysis verified against actual call graphs

**Research date:** 2026-02-26
**Valid until:** No expiration (internal codebase refactoring, no external dependencies)
