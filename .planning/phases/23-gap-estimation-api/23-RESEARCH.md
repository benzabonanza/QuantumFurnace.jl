# Phase 23: Gap Estimation API - Research

**Researched:** 2026-02-17
**Domain:** Orchestration API composing trajectory simulation, multi-observable fitting, and best-estimate selection for Lindbladian spectral gap estimation
**Confidence:** HIGH

## Summary

Phase 23 is a pure composition phase. All building blocks are implemented and tested: `run_observable_trajectories` (Phase 22) runs trajectory simulations with time-resolved observable measurement, `fit_exponential_decay` (Phase 21) fits `A * exp(-gap * t) + C` with LsqFit.jl and returns `FitResult` with gap, CI, R-squared, and convergence status, and `build_gap_estimation_observables` (Phase 20) constructs the H + M_z observable bundle in the correct basis. Phase 23 wires these together into a single `estimate_spectral_gap` function that: (1) builds default observables if none provided, (2) runs observable-only trajectories, (3) fits each observable's time series, and (4) selects the best estimate by fit quality. The result is packaged in a new `SpectralGapResult` struct.

The main technical challenge is the best-observable selection logic (GAP-03). The selection must handle edge cases: all fits may fail to converge, some may produce gap = 0 (hitting the lower bound), and R-squared can be negative (model worse than mean). The selection criteria should be: (1) converged, (2) gap > 0, (3) highest R-squared among valid fits. If no valid fit exists, the function should still return a result with diagnostic information rather than throwing an error.

No new external dependencies are needed. The function lives in a new `src/gap_estimation.jl` file (following the project's pattern of separating distinct concerns into their own files), with `SpectralGapResult` defined in the same file (co-located with its producer, matching the architecture decision that this struct is module-internal). Cross-validation against exact Liouvillian eigenvalues is deferred to Phase 24 -- Phase 23 does NOT include `cross_validate_gap`.

**Primary recommendation:** Create `estimate_spectral_gap` in `src/gap_estimation.jl` as a thin orchestration function that calls `build_gap_estimation_observables` -> `run_observable_trajectories` -> `fit_exponential_decay` per observable -> select best by R-squared. Define `SpectralGapResult` in the same file. Export both from `QuantumFurnace.jl`. Test with synthetic trajectories and the 3-qubit SMALL system.

## Standard Stack

### Core (already in project)
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| LsqFit.jl | 0.15 | Curve fitting (via `fit_exponential_decay`) | Already added in Phase 21 |
| LinearAlgebra | stdlib | Matrix operations | Already used throughout |
| Random | stdlib | Seed handling (pass-through to trajectory runner) | Already used |

### Supporting (already in project)
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| All trajectory infrastructure | Phase 13+ | `run_observable_trajectories`, `_build_framework_and_seed` | Always -- core dependency |
| Fitting infrastructure | Phase 21 | `fit_exponential_decay`, `FitResult` | Always -- core dependency |
| Observable builders | Phase 20 | `build_gap_estimation_observables` | Default observable construction |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| New `gap_estimation.jl` file | Add to existing `fitting.jl` | `fitting.jl` is pure curve fitting (no trajectory dependency). Gap estimation orchestrates trajectories + fitting + selection -- distinct concern. Separate file keeps each file focused. |
| `SpectralGapResult` in `gap_estimation.jl` | Put in `structs.jl` | Architecture decision: `SpectralGapResult` is produced and consumed only by gap estimation code. Co-locating with producer keeps `structs.jl` focused on core types. |

**Installation:** No new packages needed.

## Architecture Patterns

### Recommended Project Structure
```
src/
  convergence.jl       # Phase 20: build_gap_estimation_observables, build_total_magnetization
  fitting.jl           # Phase 21: FitResult, fit_exponential_decay
  trajectories.jl      # Phase 22: ObservableTrajectoryResult, run_observable_trajectories
  gap_estimation.jl    # Phase 23 NEW: SpectralGapResult, estimate_spectral_gap
  QuantumFurnace.jl    # Add include + exports

test/
  test_gap_estimation.jl  # Phase 23 NEW: tests for estimate_spectral_gap
  runtests.jl             # Add include
```

### Pattern 1: Thin Orchestration Function
**What:** `estimate_spectral_gap` orchestrates existing building blocks without reimplementing any logic.
**When to use:** Always. This is the entire point of the phase.
**Flow:**
```
estimate_spectral_gap(jumps, config, psi0, ham; ...)
  |
  v
build_gap_estimation_observables(ham, num_qubits; trotter)  [if observables=nothing]
  -> [H, Mz], ["H", "Mz"]
  |
  v
run_observable_trajectories(jumps, config, psi0, ham; observables, ...)
  -> ObservableTrajectoryResult (times, measurements_mean)
  |
  v
for each observable i:
  fit_exponential_decay(times, measurements_mean[i, :]; skip_initial)
  -> FitResult (gap, r_squared, converged, gap_ci, gap_se, ...)
  |
  v
Select best: converged AND gap > 0 AND highest R-squared
  |
  v
SpectralGapResult(...)
```

### Pattern 2: SpectralGapResult Struct (Co-located with Producer)
**What:** A result struct containing the gap estimate, CI, per-observable fit details, and metadata.
**Design rationale:** Phase 24 (Cross-Validation) needs the per-observable data to analyze which observables agree with the exact gap. Users need the confidence interval and best-observable name for reporting.

```julia
struct SpectralGapResult
    # Best estimate
    gap::Float64                           # gap from best-fit observable
    gap_ci::Tuple{Float64, Float64}        # 95% CI on gap
    gap_se::Float64                        # standard error on gap

    # Best observable identity
    best_observable::String                # name of observable used for gap
    best_r_squared::Float64                # R-squared of best fit

    # Per-observable results (GAP-02, GAP-03 explainability)
    per_observable::Vector{FitResult}      # all fits, one per observable
    observable_names::Vector{String}       # names corresponding to per_observable

    # Fit metadata (GAP-02)
    ntraj::Int
    total_time::Float64
    save_every::Int
    seed::Int
    skip_initial::Float64
end
```

**Key design choices:**
- Store full `FitResult` per observable (not just gap values) -- downstream Phase 24 needs R-squared, convergence status, residuals for analysis.
- No `exact_gap` or `relative_error` fields -- cross-validation is Phase 24's concern.
- `skip_initial` stored for reproducibility (user can see what fraction was skipped).
- `best_r_squared` stored separately for easy access without indexing into per_observable.

### Pattern 3: Observable Selection Logic (GAP-03)
**What:** Automatically select the best observable by fit quality.
**Algorithm:**
```
valid_fits = [(i, fit) for (i, fit) in enumerate(fits)
              if fit.converged && fit.gap > 0.0]

if isempty(valid_fits):
    # Fallback: use the fit with highest R-squared regardless of convergence
    best_idx = argmax(r_squared for fit in fits)
else:
    # Among valid fits, select highest R-squared
    best_idx = argmax(fit.r_squared for (i, fit) in valid_fits)
```

**Why R-squared over residual norm:** R-squared is scale-independent. An observable with large amplitude and large residuals may have a higher R-squared than one with small amplitude and small residuals. R-squared measures how well the model explains the variance, which is the right criterion for "did the exponential decay model capture this observable's dynamics?"

### Pattern 4: Follow Existing API Conventions
**What:** Match the calling convention of existing trajectory runners.
**Details:**
- Same positional arguments: `jumps, config, psi0, hamiltonian`
- Same keyword patterns: `trotter=nothing`, `seed=nothing`, `total_time=config.mixing_time`
- Same `delta=config.delta` pass-through
- New keywords specific to gap estimation: `observables`, `observable_names`, `ntraj`, `save_every`, `skip_initial`

```julia
function estimate_spectral_gap(
    jumps::Vector{JumpOp},
    config::AbstractThermalizeConfig,
    psi0::Vector{<:Complex},
    hamiltonian::HamHam;
    # Observable control
    observables::Union{Nothing, Vector{<:Matrix{<:Complex}}} = nothing,
    observable_names::Union{Nothing, Vector{String}} = nothing,
    # Trajectory parameters
    ntraj::Int = 1000,
    save_every::Int = 10,
    total_time::Real = config.mixing_time,
    delta::Real = config.delta,
    seed::Union{Int, Nothing} = nothing,
    trotter::Union{TrottTrott, Nothing} = nothing,
    # Fitting parameters
    skip_initial::Float64 = 0.0,
)
```

### Pattern 5: Inner/Outer Constructor for Aqua Compliance
**What:** If `SpectralGapResult` needs type parameters (it does not in the recommended design since all fields are concrete types), follow the Phase 22 pattern with explicit inner constructor.
**When to use:** Only if `SpectralGapResult` becomes parameterized. The recommended design uses concrete `Float64` types throughout, so a simple struct definition suffices.

### Anti-Patterns to Avoid
- **Hand-rolling trajectory execution:** Do NOT copy step loops or threading logic. Call `run_observable_trajectories`.
- **Hand-rolling fitting:** Do NOT access LsqFit directly. Call `fit_exponential_decay`.
- **Hand-rolling observable construction:** Do NOT build H or M_z manually. Call `build_gap_estimation_observables`.
- **Including cross-validation in Phase 23:** Cross-validation (`cross_validate_gap`) is Phase 24. Do NOT add `exact_result` parameter or cross-validation fields to `SpectralGapResult`.
- **Storing `ObservableTrajectoryResult` in `SpectralGapResult`:** The raw trajectory data is potentially large (`measurements_mean` is `n_obs x n_saves`). The user already has the trajectory result if needed. Store only metadata (ntraj, total_time, save_every, seed) in `SpectralGapResult`.
- **Modifying existing files unnecessarily:** The only existing files that need changes are `QuantumFurnace.jl` (include + exports) and `runtests.jl` (include test file). No changes to `fitting.jl`, `trajectories.jl`, `convergence.jl`, or `structs.jl`.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Observable construction | Manual H/Mz matrix building | `build_gap_estimation_observables(ham, n; trotter)` | Basis transforms are error-prone; existing function is tested |
| Trajectory simulation | Custom step loop / threading | `run_observable_trajectories(jumps, config, psi0, ham; ...)` | Thread-safe, deterministic, cross-validated against `run_trajectories` |
| Exponential fitting | Manual LsqFit calls | `fit_exponential_decay(times, values; skip_initial)` | Auto-initialization, bounds, CI, R-squared all handled |
| R-squared computation | Manual RSS/TSS | Already in `FitResult.r_squared` | Computed by `fit_exponential_decay` |
| Parameter bounds | Manual clamping | Already in `fit_exponential_decay` via `lower=[..., 0.0, ...]` | LM-integrated bounds |

**Key insight:** Phase 23 should contain approximately 100-150 lines of new code: the `SpectralGapResult` struct, `estimate_spectral_gap` function, selection logic, and a few validation helpers. Everything else is reuse.

## Common Pitfalls

### Pitfall 1: Mismatched Observable Names and Observable Matrices
**What goes wrong:** User provides `observables` but not `observable_names`, or provides them in different order, causing the "best_observable" field to point to the wrong observable.
**Why it happens:** The function accepts both as optional arguments. If only one is provided, the mapping breaks.
**How to avoid:** Validate that when `observables` is provided, `observable_names` must also be provided (or vice versa). When both are `nothing`, use `build_gap_estimation_observables` which returns matched pairs.
**Warning signs:** `best_observable` name does not match the actual observable that was fitted.

### Pitfall 2: All Fits Fail or Produce gap = 0
**What goes wrong:** The function has no valid fit to select as "best." Accessing `best_fit.gap` throws a null reference error.
**Why it happens:** Too few trajectories (noisy data), too short total_time (no decay visible), or observables have zero overlap with the gap mode.
**How to avoid:** Handle the "no valid fit" case gracefully. Return a `SpectralGapResult` with `gap = NaN`, `best_observable = "none"`, and the full `per_observable` vector so the user can diagnose the issue.
**Warning signs:** All `FitResult.converged == false` or all `FitResult.gap == 0.0`.

### Pitfall 3: skip_initial Mismatch Between FitResult and SpectralGapResult
**What goes wrong:** The `skip_initial` parameter is passed to `fit_exponential_decay` as a fraction, but stored in `SpectralGapResult` as... what? If units are inconsistent, reproducibility breaks.
**Why it happens:** `fit_exponential_decay` accepts `skip_initial::Float64` as a fraction in [0, 1). The `SpectralGapResult` should store the same value.
**How to avoid:** Store `skip_initial` as the same Float64 fraction that was passed to all `fit_exponential_decay` calls. Document clearly that it is a fraction, not a number of points.

### Pitfall 4: Observable Time Series Extraction from measurements_mean
**What goes wrong:** `measurements_mean` is `n_obs x n_saves` (rows = observables, columns = time points). Extracting `measurements_mean[i, :]` gives the time series for observable `i`. Getting the indexing wrong (column vs row) corrupts the fit.
**Why it happens:** Julia uses column-major order, and the `_accumulate_measurements!` function writes `acc[i, save_idx]` where `i` is the observable index and `save_idx` is the time index. So row = observable, column = time.
**How to avoid:** Use `measurements_mean[i, :]` (row slice) to get the time series for observable `i`. This matches the existing convention in `_accumulate_measurements!`.
**Warning signs:** Fit converges but produces nonsensical gap values.

### Pitfall 5: Default ntraj Too Low for Meaningful Fits
**What goes wrong:** With `ntraj=100`, observable means have ~10% noise per time point. The exponential fit is dominated by noise, producing wide CIs and low R-squared.
**Why it happens:** Default parameter chosen without considering the noise floor.
**How to avoid:** Default `ntraj=1000` provides ~3% noise for 3-qubit systems, sufficient for clean exponential fits. Document that larger systems need more trajectories (noise scales as 1/sqrt(ntraj)).
**Warning signs:** R-squared < 0.5 for all observables.

### Pitfall 6: Trotter Keyword Not Forwarded to Both Observable Builder and Trajectory Runner
**What goes wrong:** Observables are built in Hamiltonian eigenbasis but trajectories run in Trotter basis (or vice versa), causing basis mismatch.
**Why it happens:** The `trotter` keyword is needed by both `build_gap_estimation_observables(...; trotter)` and `run_observable_trajectories(...; trotter)`. Forgetting one or the other creates a silent basis mismatch.
**How to avoid:** Forward the `trotter` keyword to both calls. This is the exact pattern from Phase 20's `build_gap_estimation_observables(ham, n; trotter=trotter)`.
**Warning signs:** Observable expectation values at t=0 are incorrect (wrong basis).

### Pitfall 7: Not Exporting SpectralGapResult
**What goes wrong:** User cannot construct or type-check `SpectralGapResult` from outside the module.
**Why it happens:** Forgetting to add `SpectralGapResult` to the export list in `QuantumFurnace.jl`.
**How to avoid:** Add both `estimate_spectral_gap` and `SpectralGapResult` to the exports. Also add `include("gap_estimation.jl")` to the module file after `fitting.jl`.

### Pitfall 8: File Include Order in QuantumFurnace.jl
**What goes wrong:** `gap_estimation.jl` references `FitResult`, `fit_exponential_decay`, `run_observable_trajectories`, `build_gap_estimation_observables` -- all defined in other files. If `gap_estimation.jl` is included before those files, compilation fails.
**Why it happens:** Julia processes includes in order. Dependencies must be included first.
**How to avoid:** Include `gap_estimation.jl` after `fitting.jl` (which is after `convergence.jl` and `trajectories.jl`). Current order: `trajectories.jl` -> `convergence.jl` -> `fitting.jl` -> `results.jl`. Insert `gap_estimation.jl` between `fitting.jl` and `results.jl`.

## Code Examples

Verified patterns from the existing codebase, showing how Phase 23 composes them:

### Calling build_gap_estimation_observables (Phase 20)
```julia
# Source: src/convergence.jl lines 127-145
# Returns [H, Mz] in eigenbasis (or trotter basis), with names ["H", "Mz"]
obs, names = build_gap_estimation_observables(hamiltonian, num_qubits)
# With trotter basis:
obs, names = build_gap_estimation_observables(hamiltonian, num_qubits; trotter=trotter)
```

### Calling run_observable_trajectories (Phase 22)
```julia
# Source: src/trajectories.jl lines 736-749
traj_result = run_observable_trajectories(
    jumps, config, psi0, hamiltonian;
    observables=obs, save_every=10, ntraj=1000, seed=42,
    trotter=trotter, total_time=config.mixing_time,
)
# Returns ObservableTrajectoryResult with:
#   .times::Vector{Float64}           # time points
#   .measurements_mean::Matrix{Float64}  # n_obs x n_saves
#   .n_trajectories::Int
#   .seed::Int
#   .rho_mean::Nothing (when reconstruct_dm=false)
```

### Calling fit_exponential_decay (Phase 21)
```julia
# Source: src/fitting.jl lines 153-211
times = traj_result.times
obs_series = traj_result.measurements_mean[i, :]  # row i = observable i

fit = fit_exponential_decay(times, Float64.(obs_series); skip_initial=0.1)
# Returns FitResult with:
#   .gap::Float64           # decay rate
#   .amplitude::Float64     # A
#   .offset::Float64        # C
#   .gap_ci::Tuple{Float64, Float64}  # confidence interval
#   .gap_se::Float64        # standard error
#   .r_squared::Float64     # goodness of fit (can be negative)
#   .converged::Bool        # LM convergence
#   .residuals::Vector{Float64}
#   .times_used::Vector{Float64}
#   .values_used::Vector{Float64}
```

### Best Observable Selection Logic
```julia
# Source: Phase 23 design (from architecture research)
function _select_best_observable(fits::Vector{FitResult}, names::Vector{String})
    best_idx = 0
    best_r2 = -Inf

    for (i, fit) in enumerate(fits)
        if fit.converged && fit.gap > 0.0 && fit.r_squared > best_r2
            best_idx = i
            best_r2 = fit.r_squared
        end
    end

    # Fallback: if no valid fit, pick highest R-squared regardless
    if best_idx == 0
        best_idx = argmax(fit.r_squared for fit in fits)
        best_r2 = fits[best_idx].r_squared
    end

    return best_idx, names[best_idx], best_r2
end
```

### Complete estimate_spectral_gap Skeleton
```julia
function estimate_spectral_gap(
    jumps::Vector{JumpOp},
    config::AbstractThermalizeConfig,
    psi0::Vector{<:Complex},
    hamiltonian::HamHam;
    observables::Union{Nothing, Vector{<:Matrix{<:Complex}}} = nothing,
    observable_names::Union{Nothing, Vector{String}} = nothing,
    ntraj::Int = 1000,
    save_every::Int = 10,
    total_time::Real = config.mixing_time,
    delta::Real = config.delta,
    seed::Union{Int, Nothing} = nothing,
    trotter::Union{TrottTrott, Nothing} = nothing,
    skip_initial::Float64 = 0.0,
)
    # 1. Build default observables if not provided
    if observables === nothing
        observables, observable_names = build_gap_estimation_observables(
            hamiltonian, config.num_qubits; trotter=trotter)
    end
    # Validate: names must be provided with custom observables
    observable_names === nothing && throw(ArgumentError(
        "observable_names must be provided when custom observables are given"))
    length(observables) == length(observable_names) || throw(ArgumentError(
        "observables and observable_names must have same length"))

    # 2. Run observable-only trajectories
    traj_result = run_observable_trajectories(
        jumps, config, psi0, hamiltonian;
        observables=observables, save_every=save_every,
        ntraj=ntraj, total_time=total_time, delta=delta,
        seed=seed, trotter=trotter,
    )

    # 3. Fit each observable
    fits = FitResult[]
    for i in eachindex(observables)
        obs_series = Float64.(traj_result.measurements_mean[i, :])
        fit = fit_exponential_decay(traj_result.times, obs_series;
                                     skip_initial=skip_initial)
        push!(fits, fit)
    end

    # 4. Select best observable (GAP-03)
    best_idx, best_name, best_r2 = _select_best_observable(fits, observable_names)
    best_fit = fits[best_idx]

    # 5. Build result
    return SpectralGapResult(
        best_fit.gap,
        best_fit.gap_ci,
        best_fit.gap_se,
        best_name,
        best_r2,
        fits,
        observable_names,
        traj_result.n_trajectories,
        Float64(total_time),
        save_every,
        traj_result.seed,
        skip_initial,
    )
end
```

### Test Pattern: 3-Qubit System with Known Gap
```julia
# Source: test pattern from test_observable_trajectories.jl
# Use SMALL system (3-qubit, dim=8) for fast tests
psi0 = zeros(ComplexF64, SMALL_DIM); psi0[1] = 1.0
config = make_small_thermalize_config(TimeDomain();
    delta=0.01, mixing_time=2.0, with_coherent=false)

result = estimate_spectral_gap(
    SMALL_JUMPS, config, psi0, SMALL_HAM;
    ntraj=500, save_every=5, seed=42,
)

@test result isa SpectralGapResult
@test result.gap > 0.0
@test result.gap_ci[1] < result.gap_ci[2]
@test result.best_observable in ["H", "Mz"]
@test length(result.per_observable) == 2
@test length(result.observable_names) == 2
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Manual script: run trajectories, extract data, fit in separate notebook | Single `estimate_spectral_gap` call | Phase 23 | One-call API for gap estimation |
| Select observable by visual inspection | Automatic selection by R-squared | Phase 23 | Reproducible, non-subjective |
| Store results in ad-hoc named tuples | `SpectralGapResult` struct with typed fields | Phase 23 | Type-safe, introspectable |

**Important context:** The earlier architecture research (`.planning/research/ARCHITECTURE.md`) proposed a `SpectralGapResult` with many fields including `exact_gap`, `relative_error`, `fit_amplitude`, `fit_offset`, etc. The actual implementation should be simpler because: (1) cross-validation is Phase 24, not Phase 23; (2) storing per-observable `FitResult` objects provides full access to amplitude, offset, residuals without duplicating fields; (3) the `FitResult` struct already contains everything needed.

## Open Questions

1. **Should `estimate_spectral_gap` accept custom `skip_initial` per observable?**
   - What we know: Different observables may have different transient behaviors (e.g., energy decays faster than magnetization).
   - What's unclear: Whether per-observable skip_initial adds meaningful value for the current use case.
   - Recommendation: Use a single `skip_initial` for all observables in Phase 23. Per-observable fitting parameters can be added in a future extension.

2. **Should `SpectralGapResult` store the raw `ObservableTrajectoryResult`?**
   - What we know: The raw trajectory result contains the full `measurements_mean` matrix (potentially large for many observables and long simulations).
   - What's unclear: Whether downstream consumers (Phase 24, experiments) need access to the raw time series.
   - Recommendation: Do NOT store it. The user can always re-run with the same seed to reproduce. Store only metadata (ntraj, total_time, save_every, seed). If needed later, add an optional `keep_trajectory_result` flag.

3. **Should the function use `build_gap_estimation_observables` or the user-provided observables differently for TrotterDomain?**
   - What we know: `build_gap_estimation_observables(ham, n; trotter=trotter)` already handles Trotter basis transforms. `run_observable_trajectories(...; trotter=trotter)` already handles Trotter stepping.
   - Recommendation: Just forward `trotter` to both calls. No special TrotterDomain logic needed in `estimate_spectral_gap`.

4. **What file name: `gap_estimation.jl` or `spectral_gap.jl`?**
   - What we know: The architecture research proposed `spectral_gap.jl`. The existing naming convention uses descriptive names (`fitting.jl`, `convergence.jl`, `trajectories.jl`).
   - Recommendation: Use `gap_estimation.jl`. It is more descriptive of what the file does (estimating the gap), while `spectral_gap.jl` could be confused with the exact spectral gap from Liouvillian diagonalization. Either would work; this is a minor naming preference.

## Sources

### Primary (HIGH confidence)
- `src/fitting.jl` -- `FitResult` struct (10 fields), `fit_exponential_decay` function signature and behavior, `_log_linear_initial_guess`, `_compute_r_squared`. Verified by reading source code directly.
- `src/trajectories.jl` -- `ObservableTrajectoryResult` struct (lines 37-66), `run_observable_trajectories` function (lines 736-849), `_run_chunk_obs_only!`. Verified by reading source code directly.
- `src/convergence.jl` -- `build_gap_estimation_observables` (lines 127-145, returns [H, Mz] with names ["H", "Mz"]), `build_total_magnetization`. Verified by reading source code directly.
- `src/QuantumFurnace.jl` -- Module structure, include order (trajectories -> convergence -> fitting -> results), export lists. Verified by reading source code directly.
- `test/test_fitting.jl` -- Existing test patterns for `fit_exponential_decay`. Verified by reading source code directly.
- `test/test_observable_trajectories.jl` -- Existing test patterns using SMALL system. Verified by reading source code directly.
- `test/test_helpers.jl` -- SMALL_HAM, SMALL_JUMPS, SMALL_DIM, make_small_thermalize_config factory. Verified by reading source code directly.

### Secondary (MEDIUM confidence)
- `.planning/research/ARCHITECTURE.md` -- SpectralGapResult design (original proposal with 16 fields; simplified for Phase 23 scope). Informs struct design but adapted based on actual Phase 20-22 implementations.
- `.planning/research/FEATURES.md` -- estimate_spectral_gap signature proposal (line 254-262). Informs API design but adapted.
- `.planning/ROADMAP.md` -- Phase 23 requirements (GAP-01, GAP-02, GAP-03) and success criteria. Authoritative for scope.
- `.planning/REQUIREMENTS.md` -- GAP-01, GAP-02, GAP-03 requirement definitions. Authoritative.

### Tertiary (LOW confidence)
- None needed -- all Phase 23 work is internal composition of verified building blocks.

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- no new dependencies, all building blocks exist and are tested
- Architecture: HIGH -- composition pattern is straightforward; struct design informed by actual FitResult/ObservableTrajectoryResult implementations
- Pitfalls: HIGH -- all pitfalls derived from actual code patterns (measurements_mean indexing, basis mismatch, include order)
- Code examples: HIGH -- all examples use verified function signatures from source code

**Research date:** 2026-02-17
**Valid until:** No expiry -- codebase-internal research, no external dependency version concerns
