# Architecture: Spectral Gap Refinement Diagnostics

**Domain:** Diagnostic and analysis layer for spectral gap estimation quality assessment
**Researched:** 2026-02-19
**Confidence:** HIGH (direct analysis of all 26 source files, all 19 test files, existing architecture patterns fully traced through data flow)

## Executive Summary

The spectral gap refinement diagnostics milestone adds six new capabilities to the existing `estimate_spectral_gap` pipeline: exact reference comparison, Trotter convergence analysis, effective rate plots, two-exponential fitting, symmetry sector analysis, and bootstrap error bars. The central architectural question is how these capabilities integrate with the existing flat-file module layout (no subdirectories in `src/`).

**Recommendation:** Add a single new file `src/diagnostics.jl` rather than a `src/diagnostics/` subdirectory. The existing codebase uses a flat layout with one file per concern. A subdirectory would break this convention and create an artificial boundary between tightly-coupled code. The diagnostics module is a pure analysis layer that consumes existing result structs and trajectory data -- it does not modify any simulation primitives.

## Existing Architecture Summary (Relevant to Diagnostics)

### Current Data Flow for Spectral Gap Estimation

```
estimate_spectral_gap(jumps, config, psi0, hamiltonian; ...)
  |
  +-> build_preset_trajectory_observables(hamiltonian, num_qubits)
  |     Returns: (observables::Vector{Matrix{ComplexF64}}, names::Vector{String})
  |
  +-> run_observable_trajectories(jumps, config, psi0, hamiltonian; ...)
  |     Internally: _build_framework_and_seed -> TrajectoryFramework
  |                 Multi-threaded _run_chunk_obs_only! / _run_chunk_with_obs!
  |     Returns: ObservableTrajectoryResult{T}
  |       .times::Vector{Float64}                    # time grid
  |       .measurements_mean::Matrix{Float64}        # n_obs x n_saves (AVERAGES ONLY)
  |       .n_trajectories::Int
  |       .seed::Int
  |       .rho_mean::Union{Nothing, Matrix{T}}       # nothing when reconstruct_dm=false
  |
  +-> fit_exponential_decay(times, obs_series; skip_initial)   [per observable]
  |     Returns: FitResult
  |       .gap, .amplitude, .offset, .gap_ci, .gap_se, .r_squared, .converged
  |       .residuals, .times_used, .values_used
  |
  +-> _select_best_observable(fits, names)
  |
  +-> SpectralGapResult
        .gap, .gap_ci, .gap_se, .best_observable, .best_r_squared
        .per_observable::Vector{FitResult}, .observable_names
        .ntraj, .total_time, .save_every, .seed, .skip_initial
```

### Critical Observation: Individual Trajectories Are Lost

The current pipeline accumulates measurements into `mean_data_local::Matrix{Float64}` via `_accumulate_measurements!`. Individual trajectory time series are summed in-place and divided by `ntraj` at the end. **There is no mechanism to retrieve per-trajectory data** -- the chunk functions (`_run_chunk_obs_only!`, `_run_chunk_with_obs!`) accept a shared `mean_data_local` accumulator and add to it.

This is the single most important architectural constraint for bootstrap error bars. The options are:
1. Store per-trajectory measurements (memory: `ntraj * n_obs * n_saves` Float64s)
2. Store per-trajectory measurements in batches and bootstrap at the batch level
3. Re-run trajectories with different seeds (resampling the seed space)

### Existing Integration Points

| Component | Location | Relevant for | Access Pattern |
|-----------|----------|-------------|----------------|
| `construct_lindbladian()` | `furnace.jl:42-84` | Exact reference, anti-Hermitian defect | Returns full `L::Matrix{ComplexF64}` (dim^2 x dim^2) |
| `run_lindbladian()` | `furnace.jl:1-40` | Exact spectral gap reference | Returns `LindbladianResult` with `.liouvillian`, `.spectral_gap` |
| `eigenbasis_overlap_analysis()` | `gap_estimation.jl:279-327` | Symmetry analysis, mode decomposition | Takes `L`, `observables`, `rho0`; returns `OverlapAnalysisResult` |
| `fit_exponential_decay()` | `fitting.jl:153-217` | Effective rate, two-exp fitting (extend) | Takes `times`, `values`; returns `FitResult` |
| `_select_best_observable()` | `gap_estimation.jl:71-99` | Gap estimation selection logic | Internal function, used by `estimate_spectral_gap` |
| `build_preset_trajectory_observables()` | `convergence.jl:35-106` | Observable construction for all diagnostics | Returns `(observables, names)` in eigenbasis |
| `run_observable_trajectories()` | `trajectories.jl:740-853` | Trajectory data source for all diagnostics | Returns `ObservableTrajectoryResult` |
| `_run_chunk_obs_only!()` | `trajectories.jl:478-516` | Bootstrap (needs modification for per-traj storage) | Currently accumulates into shared buffer |
| `TrajectoryFramework` | `trajectories.jl:79-101` | Thread-safe read-only simulation parameters | Immutable struct, concrete-typed fields |
| `HamHam` | `hamiltonian.jl:24-39` | Eigenbasis for symmetry analysis | `.eigvals`, `.eigvecs`, `.data`, `.gibbs` |
| `TrottTrott` | `trotter_domain.jl` | Trotter convergence analysis | `.eigvals`, `.eigvecs` |

---

## Recommended Architecture

### Design Principle: Analysis Layer, Not Simulation Layer

All six diagnostic capabilities are **post-hoc analysis** of data that the existing pipeline can produce. None of them require changes to the step-along-trajectory hot loop or the Lindbladian construction machinery. The sole exception is bootstrap, which needs per-trajectory data that the current pipeline discards.

### New Files

| File | Purpose | LOC Estimate |
|------|---------|-------------|
| `src/diagnostics.jl` | All diagnostic structs and functions | 500-700 |
| `src/bootstrap.jl` | Bootstrap resampling infrastructure (trajectory runner variant + analysis) | 300-400 |
| `test/test_diagnostics.jl` | Tests for diagnostics module | 300-400 |
| `test/test_bootstrap.jl` | Tests for bootstrap module | 200-300 |

### Modified Files

| File | Change | Risk |
|------|--------|------|
| `src/QuantumFurnace.jl` | Add `include()` for new files, new exports | Minimal -- append only |
| `src/fitting.jl` | Add `fit_two_exponential_decay()` function | Low -- new function, no changes to existing |
| `src/trajectories.jl` | Add `_run_chunk_obs_per_traj!()` and `run_observable_trajectories_per_traj()` | Moderate -- new function parallel to existing patterns |

### Unchanged Files

Everything else. The diagnostic layer is purely additive.

---

## Component Design

### 1. Exact Reference Comparison

**Purpose:** Compare trajectory-estimated spectral gap against the exact gap from dense Liouvillian eigendecomposition.

**Data flow:**
```
construct_lindbladian(jumps, config, hamiltonian)    [existing]
  -> L::Matrix{ComplexF64}
  -> eigen(L)                                        [dense eigendecomposition]
  -> exact_gap = |Re(lambda_2)|
  -> compare with SpectralGapResult.gap

OR:
run_lindbladian(jumps, config_liouv, hamiltonian)    [existing]
  -> LindbladianResult.spectral_gap
```

**Implementation:** A function `compare_gap_to_exact()` that takes a `SpectralGapResult` and either a `LindbladianResult` or raw `(jumps, config, hamiltonian)`. Returns a struct with the exact gap, trajectory gap, relative error, and whether the trajectory CI contains the exact gap.

**Integration point:** Uses `construct_lindbladian()` from `furnace.jl` to get the Liouvillian, then standard `LinearAlgebra.eigen()` for the full spectrum. The `LiouvConfig` can be constructed from a `ThermalizeConfig` by dropping `mixing_time` and `delta`.

**Key design decision:** The function should accept a pre-computed `LindbladianResult` to avoid redundant Liouvillian construction when the user already has one. Also accept raw inputs for the one-call convenience API.

```julia
struct GapComparisonResult
    exact_gap::Float64
    estimated_gap::Float64
    relative_error::Float64
    ci_contains_exact::Bool
    exact_eigenvalues::Vector{ComplexF64}  # first few, sorted by |Re|
end
```

### 2. Anti-Hermitian Defect Computation

**Purpose:** Quantify how far the Lindbladian is from being purely dissipative (Hermitian in the GNS inner product). The anti-Hermitian part indicates coherent contributions or numerical errors.

**Data flow:**
```
construct_lindbladian(jumps, config, hamiltonian)    [existing]
  -> L::Matrix{ComplexF64}

anti_hermitian_defect(L)
  -> L_H = (L + L')/2    (Hermitian part)
  -> L_A = (L - L')/2    (anti-Hermitian part)
  -> defect = ||L_A|| / ||L||    (relative norm)
```

**Integration point:** Direct access to `construct_lindbladian()`. No new dependencies.

**Note:** For GNS-detailed-balance Lindbladians (with_coherent=false), this defect should be zero up to numerical precision. For KMS with coherent terms, there will be a nonzero anti-Hermitian component from the `B` operator. This diagnostic helps validate whether the coherent term is correctly implemented.

### 3. Trotter Convergence Analysis

**Purpose:** Run spectral gap estimation at multiple Trotter step counts and plot convergence toward the exact (or high-fidelity) value.

**Data flow:**
```
For each num_trotter_steps in [1, 2, 5, 10, 20, ...]:
  trotter = TrottTrott(hamiltonian, t0, num_trotter_steps)
  jumps_trotter = make_jumps_in_trotter_basis(...)
  config_trotter = ThermalizeConfig(..., num_trotter_steps_per_t0=num_trotter_steps)

  gap_result = estimate_spectral_gap(jumps_trotter, config_trotter, psi0, hamiltonian;
      trotter=trotter, ...)

  OR (for exact reference):
  liouv_config = LiouvConfig(..., num_trotter_steps_per_t0=num_trotter_steps)
  liouv_result = run_lindbladian(jumps_trotter, liouv_config, hamiltonian; trotter=trotter)

Collect: [(num_steps, gap_estimate, gap_ci, exact_gap)] for each
```

**Implementation:** A function `trotter_convergence_sweep()` that takes the base parameters and a vector of step counts, runs the estimation at each, and returns a struct collecting all results.

**Key design decision:** This is an orchestration function that calls existing entry points in a loop. It should NOT duplicate any simulation logic. The user provides the base config, and the function clones it with different `num_trotter_steps_per_t0` values.

```julia
struct TrotterConvergenceResult
    num_trotter_steps::Vector{Int}
    estimated_gaps::Vector{Float64}
    estimated_gap_cis::Vector{Tuple{Float64, Float64}}
    exact_gaps::Union{Nothing, Vector{Float64}}  # if run_exact=true
    reference_gap::Union{Nothing, Float64}        # exact gap at highest step count or from BohrDomain
end
```

### 4. Effective Rate Plots

**Purpose:** Compute the "instantaneous" or "effective" spectral gap as a function of time window, revealing multi-exponential behavior.

**Data flow:**
```
Given ObservableTrajectoryResult (times, measurements_mean):

For each window [t_start, t_end] sliding or expanding:
  fit = fit_exponential_decay(times[window], obs_series[window])
  effective_rate[t_start] = fit.gap

Plot: effective_rate vs t_start
  - Constant: clean single-exponential
  - Decreasing: multi-exponential (fast transients polluting early times)
  - Increasing: something wrong (noise dominating late times)
```

**Implementation:** A function `compute_effective_rates()` that takes a time series and computes the fitted gap over sliding or expanding windows.

**Key design decision:** Use expanding windows (fit from `t_start` to `t_end_fixed`) rather than sliding windows of fixed width. Expanding windows are more stable because they always include the long-time behavior, and the effective rate should converge to the true gap as `t_start` increases past the transient regime.

```julia
struct EffectiveRateResult
    t_starts::Vector{Float64}
    rates::Vector{Float64}
    r_squareds::Vector{Float64}
    converged::Vector{Bool}
    observable_name::String
end
```

### 5. Two-Exponential Fitting

**Purpose:** Fit `y(t) = A1*exp(-gap1*t) + A2*exp(-gap2*t) + C` to detect multi-exponential decay and separate the true spectral gap from faster-decaying transients.

**Data flow:**
```
Given (times, obs_series):
  fit_two_exponential_decay(times, obs_series)
  -> TwoExpFitResult
       .gap_slow (= spectral gap estimate)
       .gap_fast (= transient rate)
       .amplitude_slow, .amplitude_fast, .offset
       .r_squared, .converged
```

**Integration point:** Lives in `fitting.jl` alongside `fit_exponential_decay()`. Uses the same `LsqFit.curve_fit` infrastructure. The 5-parameter model (`A1, gap1, A2, gap2, C`) needs careful initialization to avoid the two exponentials collapsing onto the same rate.

**Initialization strategy:**
1. First fit single-exponential to get an initial rate
2. Use that as `gap_slow` initial guess
3. Set `gap_fast = 3 * gap_slow` as initial guess for the transient
4. Set `A1 = A_single * 0.7`, `A2 = A_single * 0.3` as initial amplitude split

**Key design decision:** The two-exp fit should be a separate function, not an overload of `fit_exponential_decay`. The interfaces differ (5 vs 3 parameters), and the initialization logic is fundamentally different. However, both should return types that `estimate_spectral_gap`-like functions can consume.

```julia
struct TwoExpFitResult
    gap_slow::Float64
    gap_fast::Float64
    amplitude_slow::Float64
    amplitude_fast::Float64
    offset::Float64
    gap_slow_ci::Tuple{Float64, Float64}
    gap_slow_se::Float64
    r_squared::Float64
    converged::Bool
    residuals::Vector{Float64}
    times_used::Vector{Float64}
    values_used::Vector{Float64}
end
```

### 6. Symmetry Sector Analysis

**Purpose:** Decompose the Lindbladian spectrum and observables by symmetry sectors (translation, spin-flip, total Sz) to understand which sectors contain the spectral gap and which observables couple to it.

**Data flow:**
```
Given L (Liouvillian), hamiltonian, observables:

eigenbasis_overlap_analysis(L, observables, names, rho0)    [EXISTING]
  -> OverlapAnalysisResult
       .eigenvalues, .exact_gap
       .overlap_coefficients (n_obs x n_modes)
       .gap_mode_overlap, .relative_gap_overlap

NEW: sector_decomposition(L, hamiltonian, sector_projectors)
  -> For each sector:
       Project L into sector: L_sector = P * L * P'
       Eigendecompose L_sector
       Find gap within sector

  -> SectorAnalysisResult
       .sector_names, .sector_gaps, .sector_eigenvalues
```

**Integration point:** Extends `eigenbasis_overlap_analysis()` from `gap_estimation.jl`. The existing function already does the core eigendecomposition and overlap computation. The new function adds symmetry-aware projectors.

**Key design decision:** Symmetry projectors should be built from the Hamiltonian's symmetries, not hardcoded. For the Heisenberg model specifically, the relevant symmetries are:
- Total Sz (conserved exactly)
- Translation (if periodic boundary conditions)
- Spin-flip Z2 (if no disordering term)

The projector construction is system-specific but the analysis machinery is generic.

```julia
struct SectorAnalysisResult
    sector_names::Vector{String}
    sector_dims::Vector{Int}
    sector_gaps::Vector{Float64}           # spectral gap within each sector
    sector_eigenvalues::Vector{Vector{ComplexF64}}  # first few eigenvalues per sector
    observable_sector_overlaps::Matrix{Float64}     # n_obs x n_sectors
    global_gap::Float64
    gap_sector::String                      # which sector contains the global gap
end
```

### 7. Bootstrap Error Bars

**Purpose:** Estimate confidence intervals on the spectral gap using bootstrap resampling of trajectory data, providing error bars that account for finite trajectory sampling noise.

**This is the only component that requires modifying the trajectory runner.**

#### The Trajectory Storage Problem

Currently, `_run_chunk_obs_only!` accumulates measurements into a shared `mean_data_local::Matrix{Float64}` (n_obs x n_saves). Individual trajectory measurements are lost after accumulation.

**Recommended approach: Batch-level bootstrap**

Instead of storing all `ntraj` individual trajectories (expensive), store measurements at the batch level:

```
Run N_batches of batch_size trajectories each.
For each batch b:
  batch_mean[b] = mean of batch_size trajectory measurements
  -> Store batch_mean[b] as one "sample" for bootstrap

Bootstrap:
  Resample batches with replacement B times
  For each bootstrap sample:
    Compute grand mean from resampled batch means
    Fit exponential decay to bootstrap grand mean
    Record fitted gap

CI = percentile of bootstrap gap distribution
```

This requires storing `N_batches x n_obs x n_saves` Float64s instead of `ntraj x n_obs x n_saves`. For typical parameters (100 batches, 8 obs, 200 saves), that is 160K Float64s (~1.3 MB) instead of (10000 traj, 8 obs, 200 saves) = 16M Float64s (~128 MB).

#### Implementation

**New trajectory runner variant:**

```julia
function run_observable_trajectories_batched(
    jumps, config, psi0, hamiltonian;
    observables, save_every, ntraj_per_batch, n_batches, seed, trotter
) -> BatchedTrajectoryResult
```

This function runs `n_batches` groups of `ntraj_per_batch` trajectories, storing per-batch mean measurements. It reuses `_build_framework_and_seed()` once and calls the existing `_run_chunk_obs_only!` machinery per batch.

**Bootstrap analysis function:**

```julia
function bootstrap_spectral_gap(
    batched_result::BatchedTrajectoryResult;
    n_bootstrap=1000, skip_initial=0.0, confidence_level=0.95
) -> BootstrapGapResult
```

This function resamples the batch-level measurements and re-fits the exponential decay for each bootstrap sample.

**New structs:**

```julia
struct BatchedTrajectoryResult
    times::Vector{Float64}
    batch_measurements::Array{Float64, 3}  # n_batches x n_obs x n_saves
    n_batches::Int
    ntraj_per_batch::Int
    total_ntraj::Int
    seed::Int
    observable_names::Vector{String}
end

struct BootstrapGapResult
    gap::Float64                    # point estimate (mean of bootstrap distribution)
    gap_median::Float64             # median of bootstrap distribution
    gap_ci::Tuple{Float64, Float64} # percentile CI
    gap_se::Float64                 # std of bootstrap distribution
    bootstrap_gaps::Vector{Float64} # full distribution
    n_bootstrap::Int
    per_observable::Vector{BootstrapObservableResult}
end

struct BootstrapObservableResult
    name::String
    gap::Float64
    gap_ci::Tuple{Float64, Float64}
    bootstrap_gaps::Vector{Float64}
end
```

---

## Integration with Existing Architecture

### Module Include Order in QuantumFurnace.jl

The include order in `QuantumFurnace.jl` follows the dependency DAG. The new files should be included at the end, after `gap_estimation.jl` and `results.jl`:

```julia
# Existing includes (unchanged)
include("fitting.jl")          # FitResult, fit_exponential_decay
include("gap_estimation.jl")   # SpectralGapResult, estimate_spectral_gap
include("results.jl")          # ExperimentResult, save/load

# New includes (added at end)
include("bootstrap.jl")        # BatchedTrajectoryResult, bootstrap infrastructure
include("diagnostics.jl")      # All diagnostic functions and result structs
```

**Rationale for order:**
- `bootstrap.jl` depends on `trajectories.jl` (TrajectoryFramework, _build_framework_and_seed) and `fitting.jl` (fit_exponential_decay)
- `diagnostics.jl` depends on `furnace.jl` (construct_lindbladian, run_lindbladian), `gap_estimation.jl` (SpectralGapResult, eigenbasis_overlap_analysis), `fitting.jl` (fit_exponential_decay, fit_two_exponential_decay), and `bootstrap.jl` (BatchedTrajectoryResult)

### Export Organization

Following the existing pattern of grouped exports:

```julia
# Diagnostics
export GapComparisonResult, compare_gap_to_exact,
       EffectiveRateResult, compute_effective_rates,
       TrotterConvergenceResult, trotter_convergence_sweep,
       SectorAnalysisResult, sector_gap_analysis,
       TwoExpFitResult, fit_two_exponential_decay

# Bootstrap
export BatchedTrajectoryResult, run_observable_trajectories_batched,
       BootstrapGapResult, bootstrap_spectral_gap
```

### Config Conversion: ThermalizeConfig -> LiouvConfig

Several diagnostics need a `LiouvConfig` to construct the Lindbladian for exact reference, but the user typically starts from a `ThermalizeConfig`. A helper is needed:

```julia
function _thermalize_to_liouv_config(config::ThermalizeConfig)
    LiouvConfig(
        num_qubits = config.num_qubits,
        with_coherent = config.with_coherent,
        with_linear_combination = config.with_linear_combination,
        domain = config.domain,
        beta = config.beta,
        sigma = config.sigma,
        a = config.a, b = config.b,
        num_energy_bits = config.num_energy_bits,
        t0 = config.t0, w0 = config.w0, eta = config.eta,
        num_trotter_steps_per_t0 = config.num_trotter_steps_per_t0,
    )
end
```

Similarly for GNS variants. This pattern already exists conceptually in the test helpers (`make_small_liouv_config` vs `make_small_thermalize_config`).

---

## Anti-Patterns to Avoid

### Anti-Pattern 1: Modifying the Trajectory Step Loop
**What:** Adding diagnostic hooks or counters inside `step_along_trajectory!`
**Why bad:** The step loop is the performance-critical hot path. It runs millions of times per gap estimation. Any added overhead there multiplies by `ntraj * num_steps`.
**Instead:** Keep all diagnostics as post-processing on the output data. The sole modification to the trajectory layer is adding a new runner variant (`run_observable_trajectories_batched`) that stores batch-level data, and this variant calls the existing step function unchanged.

### Anti-Pattern 2: Storing Full Per-Trajectory Time Series
**What:** Saving `ntraj x n_obs x n_saves` arrays for bootstrap
**Why bad:** For 10,000 trajectories with 8 observables and 200 time saves, that is 128 MB of Float64 data. The trajectory runner is already designed for accumulation, not storage.
**Instead:** Use batch-level storage (100 batches of 100 trajectories each). The batch means are sufficient for bootstrap resampling and require only ~1% of the memory.

### Anti-Pattern 3: Subclassing SpectralGapResult
**What:** Making diagnostics produce subtypes of `SpectralGapResult`
**Why bad:** Julia structs are final (no inheritance). Even if they were not, the diagnostic results contain fundamentally different information.
**Instead:** Each diagnostic produces its own result struct. A `DiagnosticReport` wrapper can bundle them with the original `SpectralGapResult`.

### Anti-Pattern 4: Making Diagnostics Depend on Each Other
**What:** Requiring `TrotterConvergenceResult` to call `compare_gap_to_exact` internally
**Why bad:** Forces the user to always compute both, even when they want only one. Creates coupling between independent analyses.
**Instead:** Each diagnostic function is standalone. A convenience `run_all_diagnostics()` can compose them, but they are independently callable.

---

## Patterns to Follow

### Pattern 1: Existing Result Struct Convention
All existing result structs (`SpectralGapResult`, `FitResult`, `LindbladianResult`, `OverlapAnalysisResult`) are flat immutable structs with no methods. Follow this pattern:
- Fields are plain types (Float64, Vector, Matrix, Tuple)
- No methods attached to structs
- Functions take structs as arguments
- Constructor is the default positional or @kwdef constructor

### Pattern 2: Existing Test Fixture Pattern
Tests use shared fixtures from `test_helpers.jl` with `SMALL_HAM`, `SMALL_JUMPS`, `SMALL_GIBBS` for fast tests (3-qubit, dim=8) and `TEST_HAM`, `TEST_JUMPS`, `TEST_GIBBS` for higher-fidelity tests (4-qubit, dim=16). New diagnostic tests should use the same fixtures.

### Pattern 3: Internal Functions Prefixed with _
Functions not exported are prefixed with `_`. The diagnostic module should follow this for all helpers.

### Pattern 4: Framework-Once, Run-Many
The trajectory infrastructure is designed around building `TrajectoryFramework` once and reusing it for many trajectory runs. The batched trajectory runner should use this same pattern: call `_build_framework_and_seed()` once, then loop over batches.

---

## Data Flow Diagrams

### Full Diagnostic Pipeline (One-Call API)

```
run_gap_diagnostics(jumps, config, psi0, hamiltonian; ...)
  |
  +-> estimate_spectral_gap(...)                      [existing, unchanged]
  |     -> SpectralGapResult
  |
  +-> compare_gap_to_exact(gap_result, jumps, config, hamiltonian)
  |     -> construct_lindbladian(...)                  [existing]
  |     -> eigen(L)
  |     -> GapComparisonResult
  |
  +-> run_observable_trajectories_batched(...)         [NEW runner]
  |     -> BatchedTrajectoryResult
  |
  +-> bootstrap_spectral_gap(batched_result)
  |     -> resample batch means B times
  |     -> fit_exponential_decay per resample
  |     -> BootstrapGapResult
  |
  +-> compute_effective_rates(times, obs_series)
  |     -> expanding window fits
  |     -> EffectiveRateResult
  |
  +-> fit_two_exponential_decay(times, obs_series)    [NEW in fitting.jl]
  |     -> TwoExpFitResult
  |
  +-> sector_gap_analysis(L, hamiltonian, observables)
  |     -> build symmetry projectors
  |     -> project L into sectors
  |     -> eigendecompose per sector
  |     -> SectorAnalysisResult
  |
  +-> DiagnosticReport (bundles everything)
```

### Bootstrap Data Flow (Detailed)

```
run_observable_trajectories_batched(jumps, config, psi0, hamiltonian;
    observables, ntraj_per_batch=100, n_batches=100, ...)
  |
  +-> _build_framework_and_seed(...)                  [existing, called once]
  |     -> fw::TrajectoryFramework, actual_seed::Int
  |
  +-> For batch_idx in 1:n_batches:
  |     batch_seed = actual_seed + (batch_idx - 1) * ntraj_per_batch
  |     mean_data_batch = zeros(n_obs, n_saves)
  |
  |     [Multi-threaded within batch: same pattern as run_observable_trajectories]
  |     _run_chunk_obs_only!(ws, fw, psi0, chunk, batch_seed, ...)
  |       -> accumulates into mean_data_batch
  |
  |     mean_data_batch ./= ntraj_per_batch
  |     batch_measurements[batch_idx, :, :] = mean_data_batch
  |
  +-> BatchedTrajectoryResult(times, batch_measurements, ...)

bootstrap_spectral_gap(batched_result; n_bootstrap=1000, ...)
  |
  +-> grand_mean = mean(batch_measurements, dims=1)   # for point estimate
  |
  +-> For b in 1:n_bootstrap:
  |     resampled_indices = rand(1:n_batches, n_batches)  # with replacement
  |     resampled_mean = mean(batch_measurements[resampled_indices, :, :], dims=1)
  |     For each observable:
  |       fit = fit_exponential_decay(times, resampled_mean[obs_idx, :])
  |       record fit.gap
  |
  +-> Compute percentile CI from bootstrap gap distribution
  +-> BootstrapGapResult
```

---

## Scalability Considerations

| Concern | 3-qubit (dim=8) | 4-qubit (dim=16) | 5-qubit (dim=32) |
|---------|-----------------|-------------------|-------------------|
| Lindbladian construction | 64x64 matrix, <1s | 256x256 matrix, ~2s | 1024x1024 matrix, ~30s |
| Dense eigendecomposition | 64x64, instant | 256x256, <1s | 1024x1024, ~10s |
| Trajectory simulation (1000 traj) | ~5s | ~30s | ~5min |
| Batched trajectories (100 batches x 100 traj) | ~50s | ~5min | ~50min |
| Bootstrap (1000 resamples) | <1s (fitting only) | <1s (fitting only) | <1s (fitting only) |
| Sector analysis (build projectors + project + eigen) | Trivial | ~2s | ~30s |
| Memory: batch measurements (100 x 8 x 200) | 1.3 MB | 1.3 MB | 1.3 MB |

**Key insight:** The bottleneck is always trajectory simulation, never the analysis. Bootstrap resampling and fitting is orders of magnitude cheaper than running trajectories. This validates the batch-level bootstrap approach: the dominant cost is running the batched trajectories, and the bootstrap analysis itself is negligible.

---

## Suggested Build Order

The build order follows the dependency chain and enables incremental testing:

### Phase 1: Two-Exponential Fitting (0 dependencies on new code)

**File:** `src/fitting.jl` (add `fit_two_exponential_decay`)
**Why first:** Pure numerical function with no dependencies on the rest of the diagnostic infrastructure. Can be tested immediately with synthetic data. Establishes the `TwoExpFitResult` struct that other diagnostics will reference.
**Tests:** Synthetic two-exponential data recovery, initialization from single-exp fit, degenerate case (two rates collapse to one).

### Phase 2: Effective Rate Analysis (depends on Phase 1)

**File:** `src/diagnostics.jl` (begin file)
**Why second:** Pure analysis of existing time series data. Uses `fit_exponential_decay` (existing) over sliding windows. No trajectory runner changes needed.
**Tests:** Synthetic multi-exponential data showing rate convergence, constant-rate data showing flat effective rate.

### Phase 3: Exact Reference Comparison (0 dependencies on new code)

**File:** `src/diagnostics.jl` (add to file)
**Why third:** Uses existing `construct_lindbladian()` and `run_lindbladian()`. Needs the `_thermalize_to_liouv_config()` helper. Independent of Phases 1-2.
**Tests:** Compare trajectory gap against exact gap for 3-qubit system, verify CI containment.

### Phase 4: Anti-Hermitian Defect (0 dependencies on new code)

**File:** `src/diagnostics.jl` (add to file)
**Why fourth:** Trivial function on a matrix. Can be done at any time. Grouped here because it pairs logically with exact reference.
**Tests:** Zero defect for Hermitian matrix, known defect for skew-Hermitian addition.

### Phase 5: Symmetry Sector Analysis (depends on Phase 3)

**File:** `src/diagnostics.jl` (add to file)
**Why fifth:** Needs Liouvillian from Phase 3's infrastructure. The projector construction is the complex part.
**Tests:** Verify total Sz conservation for Heisenberg model, verify gap sector identification.

### Phase 6: Batched Trajectory Runner (0 dependencies on diagnostics)

**File:** `src/bootstrap.jl` (begin file)
**Why sixth:** The new trajectory runner variant. Independent of all diagnostics. This is the highest-risk change because it touches the trajectory infrastructure (though only adding a new function, not modifying existing ones).
**Tests:** Verify batched results match non-batched grand mean, verify per-batch data storage.

### Phase 7: Bootstrap Analysis (depends on Phase 6)

**File:** `src/bootstrap.jl` (add to file)
**Why seventh:** Requires Phase 6's `BatchedTrajectoryResult`. The bootstrap resampling itself is straightforward.
**Tests:** Known-distribution bootstrap coverage, CI width decreases with ntraj.

### Phase 8: Trotter Convergence Sweep (0 dependencies on diagnostics)

**File:** `src/diagnostics.jl` (add to file)
**Why eighth:** Orchestration function that loops over existing `estimate_spectral_gap`. Independent of other diagnostics but computationally expensive to test. Placed last because it is the least likely to have subtle bugs (it is just a loop).
**Tests:** Verify convergence trend for 3-qubit system with 2-3 Trotter step counts.

### Phase 9: Integration (depends on all above)

**File:** `src/diagnostics.jl` (add convenience wrapper)
**Why last:** Bundle all diagnostics into a single-call `run_gap_diagnostics()` API.
**Tests:** End-to-end test on 3-qubit system.

---

## Result Struct Summary

| Struct | Fields | Produced By | Lives In |
|--------|--------|-------------|----------|
| `TwoExpFitResult` | gap_slow, gap_fast, amplitudes, offset, CI, R2, residuals | `fit_two_exponential_decay()` | `fitting.jl` |
| `EffectiveRateResult` | t_starts, rates, r_squareds, converged, observable_name | `compute_effective_rates()` | `diagnostics.jl` |
| `GapComparisonResult` | exact_gap, estimated_gap, relative_error, ci_contains_exact | `compare_gap_to_exact()` | `diagnostics.jl` |
| `TrotterConvergenceResult` | num_steps[], gaps[], CIs[], exact_gaps[] | `trotter_convergence_sweep()` | `diagnostics.jl` |
| `SectorAnalysisResult` | sector_names, sector_gaps, sector_eigenvalues, overlaps | `sector_gap_analysis()` | `diagnostics.jl` |
| `BatchedTrajectoryResult` | times, batch_measurements (3D array), metadata | `run_observable_trajectories_batched()` | `bootstrap.jl` |
| `BootstrapGapResult` | gap, gap_ci, gap_se, bootstrap_gaps[], per_observable | `bootstrap_spectral_gap()` | `bootstrap.jl` |
| `BootstrapObservableResult` | name, gap, gap_ci, bootstrap_gaps[] | (part of BootstrapGapResult) | `bootstrap.jl` |
| `DiagnosticReport` | comparison, effective_rates, two_exp_fits, sector, bootstrap | `run_gap_diagnostics()` | `diagnostics.jl` |

---

## Sources

- Direct analysis of all source files in `src/` (26 files, ~6,274 LOC)
- Direct analysis of all test files in `test/` (19 files, ~4,366 LOC)
- Existing architecture documentation: `.planning/codebase/ARCHITECTURE.md`
- Previous milestone research: `.planning/research/ARCHITECTURE.md` (v1.3 spectral gap estimation)
- Julia LsqFit.jl documentation for multi-parameter fitting patterns
- Bootstrap methodology: Efron & Tibshirani (1993) batch bootstrap for dependent data
