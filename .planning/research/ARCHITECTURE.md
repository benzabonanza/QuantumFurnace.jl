# Architecture: Spectral Gap Estimation from Trajectory Observable Decay

**Domain:** Observable-only trajectory runner + exponential fitting for Lindbladian spectral gap estimation
**Researched:** 2026-02-16
**Confidence:** HIGH (direct codebase analysis of all source files, existing architecture patterns well understood from v1.0-v1.2)

## System Overview: New Components in Existing Architecture

```
EXISTING (read-only or lightly modified)         NEW (this milestone)
========================================         ====================

trajectories.jl                                  spectral_gap.jl (NEW)
  TrajectoryFramework (read-only)                  SpectralGapResult struct
  TrajectoryWorkspace (per-thread)                 estimate_spectral_gap()
  step_along_trajectory!(psi, fw, ws, rng)           -> calls run_observable_trajectories
  _accumulate_measurements!                          -> calls fit_exponential_decay
  _accumulate_density_matrix!                        -> optionally cross-validates vs exact
  _run_chunk_with_obs!  (TEMPLATE)
  _run_batch_no_obs!    (reused for final DM)      run_observable_trajectories() (NEW)
  run_trajectories(...; observables)                  -> _run_chunk_obs_only! (NEW)
                                                      -> observable-only hot loop
convergence.jl                                       -> DM once at end via _run_batch_no_obs!
  run_trajectories_convergence
  run_trajectories_adaptive                        fit_exponential_decay() (NEW)
  build_convergence_observables                      -> model: f(t) = A*exp(-lambda*t) + C
                                                     -> uses LsqFit.jl curve_fit with bounds
furnace.jl
  run_lindbladian -> LindbladianResult             convergence.jl (MODIFIED)
    .spectral_gap (exact reference)                  build_total_magnetization(ham, n_qubits)
                                                     build_gap_estimation_observables(ham, n)

structs.jl
  LindbladianResult{T} (.spectral_gap)
  TrajectoryResult (measurements_mean, times)
  ConvergenceData
```

## Recommended Architecture

### Principle: Compose Existing Primitives, Add Minimal New Code

The existing codebase provides all the trajectory simulation machinery. The new milestone needs:

1. A trajectory runner variant that measures observables during simulation but skips per-trajectory DM reconstruction (reconstructing the DM only once from a separate final run or from the same trajectories at the end)
2. An exponential fitting function to extract decay rate from observable time series
3. A result struct to package the estimated gap with metadata
4. A total magnetization observable builder
5. A cross-validation helper that compares trajectory-estimated gap against exact eigenvalue gap

### Component 1: Observable-Only Trajectory Runner

**What exists:** `_run_chunk_with_obs!` measures observables at `save_every` intervals AND accumulates DM at the end of each trajectory. The DM accumulation (`_accumulate_density_matrix!`) is a rank-1 update (O(dim^2)) at the end of each trajectory -- cheap for dim=16 (n=4) but nontrivial for dim=256 (n=8).

**What is needed:** A trajectory runner that:
- Measures observables at regular intervals (same as `_run_chunk_with_obs!`)
- Does NOT reconstruct the DM per trajectory during the observable measurement pass
- Optionally reconstructs the DM once at the end (by running `_run_batch_no_obs!` separately or by accumulating DM from the same trajectories at the final step only)

**Key insight:** The existing `_run_chunk_with_obs!` already does both observable measurement and DM accumulation. The "observable-only" variant is a simplification: remove the `_accumulate_density_matrix!` call from the per-trajectory loop, and add it back only at the very end (or as a separate pass). This is a small modification.

**Recommended approach -- new function `_run_chunk_obs_only!`:**

```julia
function _run_chunk_obs_only!(
    ws::TrajectoryWorkspace{<:Complex},
    fw::TrajectoryFramework{<:Complex},
    psi0::Vector{<:Complex},
    chunk::UnitRange{Int},
    master_seed::Int,
    total_time::Real,
    observables::Vector{<:Matrix{<:Complex}},
    save_every::Int,
    num_steps::Int,
    num_saves::Int,
    mean_data_local::Matrix{Float64},
    # NEW: optional final-step DM accumulation
    accumulate_final_dm::Bool,
)
    psi = copy(psi0)
    tmp_meas = ws.psi_tmp

    for traj_id in chunk
        rng = Random.Xoshiro(master_seed + traj_id)
        copyto!(psi, psi0)
        n2 = real(dot(psi, psi))
        rmul!(psi, 1.0 / sqrt(max(n2, eps(Float64))))

        # Measure at t=0
        _accumulate_measurements!(mean_data_local, 1, psi, observables, tmp_meas)

        save_idx = 1
        for step in 1:num_steps
            step_along_trajectory!(psi, fw, ws, rng)
            if step % save_every == 0
                save_idx += 1
                _accumulate_measurements!(mean_data_local, save_idx, psi, observables, tmp_meas)
            end
        end

        # DM accumulation: ONLY if requested (for final-step DM)
        if accumulate_final_dm
            _accumulate_density_matrix!(ws.rho_acc, psi)
        end
    end
    return nothing
end
```

**Why not just use `_run_chunk_with_obs!` as-is:** The existing function always accumulates DM. For the spectral gap estimation use case, the user wants many short trajectories with frequent observable snapshots to fit the decay curve, and DM reconstruction is secondary. Skipping per-trajectory DM accumulation (or making it optional) clarifies intent and slightly reduces memory pressure on the workspace `rho_acc` buffer (though the cost is negligible for current system sizes).

**Alternative considered and rejected -- reuse existing `run_trajectories` with `accumulate_final_dm=true`:** Adding a boolean flag to the existing `_run_chunk_with_obs!` would work, but pollutes the existing function with a parameter that 95% of callers do not need. A separate function is cleaner. The new function can share the inner loop structure via copy (both functions are ~25 lines; DRY is not a concern at this size).

### Component 2: Public API -- `run_observable_trajectories`

This is the public entry point for observable-only trajectory simulation:

```julia
"""
    run_observable_trajectories(jumps, config, psi0, hamiltonian;
        observables, save_every, ntraj, total_time, delta, seed,
        trotter, reconstruct_dm)

Run trajectories measuring observables at regular intervals.
Unlike run_trajectories, DM reconstruction is optional (controlled
by `reconstruct_dm`; default false). When false, rho_mean in the
returned TrajectoryResult is zeros (not computed).

Returns TrajectoryResult with times and measurements_mean populated.
"""
function run_observable_trajectories(
    jumps::Vector{JumpOp},
    config::AbstractThermalizeConfig,
    psi0::Vector{<:Complex},
    hamiltonian::HamHam;
    observables::Vector{<:Matrix{<:Complex}},
    save_every::Int = 1,
    ntraj::Int = 1000,
    total_time::Real = config.mixing_time,
    delta::Real = config.delta,
    seed::Union{Int,Nothing} = nothing,
    trotter::Union{TrottTrott,Nothing} = nothing,
    reconstruct_dm::Bool = false,
)
```

**Integration with existing architecture:** This function follows the exact same pattern as `run_trajectories`:
1. Call `_build_framework_and_seed` (existing) to get `(fw, actual_seed)`
2. Compute `num_steps`, `num_saves`, `times` (same logic as `run_trajectories`)
3. Dispatch to multi-threaded or serial `_run_chunk_obs_only!` (new)
4. Return `TrajectoryResult` with `measurements_mean` and `times` populated, `rho_mean` zeros or computed depending on `reconstruct_dm`

**Where this lives:** `src/spectral_gap.jl` (new file). Keeping it separate from `trajectories.jl` prevents the already-large trajectories file (927 lines) from growing further, and groups all spectral-gap-related code together.

### Component 3: Exponential Fitting

**The physics:** For a Lindbladian with spectral gap lambda (the real part of the second-smallest eigenvalue of the Liouvillian), observables decay exponentially toward their steady-state value:

```
<O>(t) - <O>_ss  ~  A * exp(-lambda * t)  as t -> infinity
```

where `<O>_ss = tr(O * rho_ss)` is the steady-state expectation value. The decay rate lambda equals the real part of the Lindbladian spectral gap (|Re(eigenvalue_2)|).

In practice, the observable time series has:
- An initial transient (first few steps) where the system has not yet "felt" the gap mode
- An exponential decay regime (most of the time series)
- Statistical noise from finite trajectory averaging
- A plateau at the steady-state value

**Fitting model:**

```
f(t, p) = p[1] * exp(-p[2] * t) + p[3]
```

where:
- `p[1]` = amplitude A (can be positive or negative)
- `p[2]` = decay rate lambda (this is the spectral gap estimate)
- `p[3]` = steady-state offset C = <O>_ss

**Implementation using LsqFit.jl:**

LsqFit.jl is the standard Julia package for nonlinear curve fitting. It provides Levenberg-Marquardt optimization with parameter bounds, automatic Jacobian via ForwardDiff, and built-in confidence intervals / standard errors. This is preferred over raw Optim.jl because:

1. `curve_fit` with `lower`/`upper` bounds constrains `gap > 0` natively
2. `confidence_interval(fit, alpha)` gives t-distribution-based CIs on the gap
3. `standard_error(fit)` gives Jacobian-based parameter uncertainties
4. The transitive deps (ForwardDiff, NLSolversBase) are already in the Manifest via Optim

```julia
using LsqFit

function fit_exponential_decay(
    times::Vector{Float64},
    values::Vector{Float64};
    p0::Union{Nothing, Vector{Float64}} = nothing,
    skip_initial::Int = 0,
    gibbs_value::Union{Nothing, Float64} = nothing,
)
    # Trim initial transient
    t = times[skip_initial+1:end]
    y = values[skip_initial+1:end]

    # Model
    model(t, p) = p[1] .* exp.(-p[2] .* t) .+ p[3]

    # Auto-initialize if no guess provided
    if p0 === nothing
        C_guess = gibbs_value !== nothing ? gibbs_value : mean(y[max(1,end-div(length(y),5)):end])
        A_guess = y[1] - C_guess
        # Log-linear estimate for gap
        shifted = abs.(y .- C_guess) .+ 1e-15
        log_y = log.(shifted)
        gap_guess = -(log_y[end] - log_y[1]) / (t[end] - t[1])
        gap_guess = max(gap_guess, 1e-6)
        p0 = [A_guess, gap_guess, C_guess]
    end

    # Bounds: gap must be positive
    lower = [-Inf, 0.0, -Inf]
    upper = [Inf, Inf, Inf]

    fit = curve_fit(model, t, y, p0; lower=lower, upper=upper)

    # Extract results
    params = fit.param
    gap_ci = confidence_interval(fit, 0.05)[2]  # 95% CI on gap
    gap_se = standard_error(fit)[2]

    return (
        amplitude = params[1],
        decay_rate = params[2],   # This is the spectral gap estimate
        offset = params[3],
        converged = fit.converged,
        residual_norm = sqrt(sum(fit.resid.^2) / length(t)),
        confidence_interval = gap_ci,
        standard_error = gap_se,
        fit_result = fit,  # Raw LsqFit result for further analysis
    )
end
```

**Why LsqFit.jl over Optim.jl:** While Optim.jl is already a dependency and can do nonlinear optimization, LsqFit.jl wraps Levenberg-Marquardt with purpose-built curve fitting infrastructure: automatic Jacobian computation, covariance estimation, confidence intervals, and parameter bounds. Reimplementing these with raw Optim.jl would be error-prone and redundant. LsqFit.jl adds one direct dependency; its transitive deps (ForwardDiff, NLSolversBase, StatsAPI) are already resolved in the Manifest.toml via Optim. The only truly new transitive package is Distributions.jl (used internally for t-distribution quantiles in confidence intervals).

**Why Levenberg-Marquardt over NelderMead:** LM uses gradient information (via the Jacobian) for faster convergence and provides the Jacobian-based covariance matrix needed for confidence intervals. NelderMead is derivative-free but cannot produce uncertainty estimates without additional work.

**Alternative considered:** Linear regression on log-transformed data `log(|<O>(t) - C|) = log(A) - lambda*t`. This avoids nonlinear optimization but requires knowing `C` (the steady-state value) a priori, and breaks when `<O>(t) - C` changes sign (which happens with noisy data). The nonlinear fit is more robust. However, the log-linear approach is used as the initial guess strategy (see STACK.md for details).

### Component 4: Result Struct

```julia
"""
Spectral gap estimation result from trajectory observable decay.
"""
struct SpectralGapResult{T<:AbstractFloat}
    # Estimated gap
    estimated_gap::T                    # |lambda| from exponential fit
    gap_confidence_interval::Tuple{T,T} # 95% CI from LsqFit
    gap_standard_error::T               # SE from Jacobian covariance
    fit_amplitude::T                    # A from fit
    fit_offset::T                       # C from fit (steady-state estimate)
    fit_converged::Bool                 # Did LM converge?
    fit_residual_norm::T                # RMS residual of fit

    # Observable used
    observable_name::String             # Which observable was fitted
    observable_index::Int               # Index in the observables vector

    # Per-observable fits (for multi-observable consistency check)
    all_gap_estimates::Vector{T}        # Gap from each observable
    all_observable_names::Vector{String}

    # Cross-validation (optional, only for small systems)
    exact_gap::Union{Nothing, Complex{T}}  # From run_lindbladian, if available
    relative_error::Union{Nothing, T}      # |estimated - |Re(exact)|| / |Re(exact)|

    # Trajectory metadata
    n_trajectories::Int
    total_time::T
    save_every::Int
    seed::Int
end
```

**Where this lives:** `src/spectral_gap.jl` (with other gap estimation code). Unlike `LindbladianResult` and `TrajectoryResult` which are used across multiple files, `SpectralGapResult` is produced and consumed only by the gap estimation module. Co-locating the struct with its producer keeps `structs.jl` focused on core types and makes the gap estimation module self-contained.

### Component 5: Total Magnetization Observable

**What exists:** `build_convergence_observables` in `convergence.jl` builds nearest-neighbor `Z_iZ_{i+1}` correlations and energy `<H>` in the Hamiltonian eigenbasis.

**What is needed:** Total magnetization `M_z = sum_i Z_i` in the Hamiltonian eigenbasis.

```julia
"""
    build_total_magnetization(hamiltonian::HamHam, num_qubits::Int) -> Matrix{ComplexF64}

Build the total magnetization operator M_z = sum_i Z_i in the Hamiltonian eigenbasis.
"""
function build_total_magnetization(hamiltonian::HamHam, num_qubits::Int)
    V = hamiltonian.eigvecs
    dim = 2^num_qubits
    Mz_comp = zeros(ComplexF64, dim, dim)
    for i in 1:num_qubits
        Mz_comp .+= pad_term([Z], num_qubits, i)
    end
    return V' * Mz_comp * V  # Transform to eigenbasis
end
```

**Where this lives:** `src/convergence.jl` alongside the existing `build_convergence_observables`. Add a convenience function that returns the full set needed for gap estimation:

```julia
"""
    build_gap_estimation_observables(hamiltonian::HamHam, num_qubits::Int)

Build observables for spectral gap estimation: energy + total magnetization + ZZ correlations.
Returns (observables, names).
"""
function build_gap_estimation_observables(hamiltonian::HamHam, num_qubits::Int)
    observables = Matrix{ComplexF64}[]
    names = String[]

    V = hamiltonian.eigvecs

    # Energy (diagonal in eigenbasis) -- strongest gap coupling
    H_eigen = Matrix{ComplexF64}(diagm(ComplexF64.(hamiltonian.eigvals)))
    push!(observables, H_eigen)
    push!(names, "H")

    # Total magnetization
    Mz = build_total_magnetization(hamiltonian, num_qubits)
    push!(observables, Mz)
    push!(names, "Mz_total")

    # ZZ correlations (following build_convergence_observables pattern)
    for i in 1:num_qubits-1
        ZZ_comp = Matrix{ComplexF64}(pad_term([Z, Z], num_qubits, i))
        ZZ_eigen = V' * ZZ_comp * V
        push!(observables, ZZ_eigen)
        push!(names, "ZZ_$(i)$(i+1)")
    end

    return observables, names
end
```

### Component 6: Cross-Validation Against Exact Eigenvalues

For n=4 (dim=16) and n=6 (dim=64), the full Liouvillian can be constructed and diagonalized. The exact spectral gap from `run_lindbladian` serves as ground truth.

```julia
"""
    cross_validate_gap(estimated::SpectralGapResult, exact_result::LindbladianResult)

Compare trajectory-estimated spectral gap against exact Liouvillian eigenvalue.
Returns the relative error |estimated - |Re(exact)|| / |Re(exact)|.

IMPORTANT: Uses -real(spectral_gap), NOT abs(spectral_gap). The imaginary
part of the eigenvalue produces oscillations, not decay.
"""
function cross_validate_gap(
    estimated::SpectralGapResult,
    exact_result::LindbladianResult,
)
    exact_gap_real = abs(real(exact_result.spectral_gap))
    imag_fraction = abs(imag(exact_result.spectral_gap)) / exact_gap_real

    if imag_fraction > 0.1
        @warn "Liouvillian gap eigenvalue has significant imaginary part " *
              "(|Im/Re| = $(round(imag_fraction, digits=3))). " *
              "Observable decay may show oscillations; pure exponential fit may be poor."
    end

    rel_error = abs(estimated.estimated_gap - exact_gap_real) / exact_gap_real
    return rel_error
end
```

**Where this lives:** `src/spectral_gap.jl` alongside the other gap estimation code.

### Component 7: Top-Level API -- `estimate_spectral_gap`

The main entry point that orchestrates everything:

```julia
"""
    estimate_spectral_gap(jumps, config, psi0, hamiltonian;
        observables=nothing, observable_names=nothing,
        ntraj=1000, save_every=10, total_time=config.mixing_time,
        seed=nothing, trotter=nothing, skip_initial=0,
        exact_result=nothing, reconstruct_dm=false)

Estimate the Lindbladian spectral gap from trajectory observable decay.

Runs observable-only trajectories, fits exponential decay to each observable,
and returns the best gap estimate (lowest fit residual with valid decay rate).

If `exact_result::LindbladianResult` is provided, cross-validates against
the exact spectral gap.

Returns SpectralGapResult.
"""
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
    seed::Union{Int,Nothing} = nothing,
    trotter::Union{TrottTrott,Nothing} = nothing,
    skip_initial::Int = 0,
    exact_result::Union{Nothing, LindbladianResult} = nothing,
    reconstruct_dm::Bool = false,
)
    # Default observables: energy + total magnetization + ZZ correlations
    if observables === nothing
        observables, observable_names = build_gap_estimation_observables(
            hamiltonian, config.num_qubits)
    end

    # Run observable-only trajectories
    traj_result = run_observable_trajectories(
        jumps, config, psi0, hamiltonian;
        observables=observables, save_every=save_every,
        ntraj=ntraj, total_time=total_time, seed=seed,
        trotter=trotter, reconstruct_dm=reconstruct_dm,
    )

    # Fit exponential decay to each observable
    all_gap_estimates = Float64[]
    best_fit = nothing
    best_residual = Inf

    for (i, name) in enumerate(observable_names)
        obs_series = traj_result.measurements_mean[i, :]
        fit = fit_exponential_decay(
            traj_result.times, obs_series;
            skip_initial=skip_initial,
        )
        push!(all_gap_estimates, fit.converged && fit.decay_rate > 0 ? fit.decay_rate : NaN)

        if fit.converged && fit.decay_rate > 0 && fit.residual_norm < best_residual
            best_fit = fit
            best_residual = fit.residual_norm
            best_obs_name = name
            best_obs_idx = i
        end
    end

    # Cross-validate if exact result available
    exact_gap = exact_result !== nothing ? exact_result.spectral_gap : nothing
    rel_error = if exact_gap !== nothing && best_fit !== nothing
        abs(best_fit.decay_rate - abs(real(exact_gap))) / abs(real(exact_gap))
    else
        nothing
    end

    actual_seed = traj_result.seed

    return SpectralGapResult(
        best_fit.decay_rate,
        best_fit.confidence_interval,
        best_fit.standard_error,
        best_fit.amplitude,
        best_fit.offset,
        best_fit.converged,
        best_fit.residual_norm,
        best_obs_name,
        best_obs_idx,
        all_gap_estimates,
        observable_names,
        exact_gap,
        rel_error,
        ntraj,
        Float64(total_time),
        save_every,
        actual_seed,
    )
end
```

## Component Boundaries

### New Files

| File | Purpose | Key Types/Functions |
|------|---------|-------------------|
| `src/spectral_gap.jl` | Spectral gap estimation from observable decay | `SpectralGapResult`, `estimate_spectral_gap`, `run_observable_trajectories`, `fit_exponential_decay`, `_run_chunk_obs_only!`, `cross_validate_gap` |

### Modified Files

| File | What Changes | Why |
|------|-------------|-----|
| `src/convergence.jl` | Add `build_total_magnetization`, `build_gap_estimation_observables` | Observable builders grouped with existing `build_convergence_observables` |
| `src/QuantumFurnace.jl` | Add `include("spectral_gap.jl")`, `using LsqFit`, export new public API | Module registration + new dependency |
| `Project.toml` | Add `LsqFit = "2fda8390-95c7-5789-9bda-21331edee243"` to `[deps]` and `LsqFit = "0.15"` to `[compat]` | New dependency |

### Unchanged Files

| File | Why Unchanged |
|------|--------------|
| `src/trajectories.jl` | All primitives (`step_along_trajectory!`, `_accumulate_measurements!`, `_build_framework_and_seed`, `_run_batch_no_obs!`, `_partition_trajectories`) are reused as-is |
| `src/furnace.jl` | `run_lindbladian` already returns `LindbladianResult` with `.spectral_gap` -- used as-is for cross-validation |
| `src/structs.jl` | No changes needed; `SpectralGapResult` lives in `spectral_gap.jl` |
| `src/results.jl` | `ExperimentResult` serialization not needed for this milestone (can be added later) |

## Data Flow

### Gap Estimation Flow (Complete)

```
User calls: estimate_spectral_gap(jumps, config, psi0, ham; ntraj=1000)
    |
    v
build_gap_estimation_observables(ham, n_qubits)
    -> [H (eigenbasis), Mz_total (eigenbasis), ZZ_12, ZZ_23, ...]
    |
    v
run_observable_trajectories(jumps, config, psi0, ham; observables, save_every, ntraj)
    |
    |  _build_framework_and_seed(jumps, config, psi0, ham)  [EXISTING]
    |    -> (fw::TrajectoryFramework, seed)
    |
    |  Compute num_steps, num_saves, times  [same logic as run_trajectories]
    |
    |  Multi-threaded dispatch:
    |    For each thread chunk:
    |      _run_chunk_obs_only!(ws, fw, psi0, chunk, seed, total_time,
    |        observables, save_every, num_steps, num_saves, mean_data_local,
    |        accumulate_final_dm=false)
    |
    |  Reduce: mean_data = sum(per_thread) / ntraj
    |
    |  Return TrajectoryResult(rho_mean=zeros, ..., times=times,
    |    measurements_mean=mean_data)
    |
    v
For each observable i:
    fit_exponential_decay(times, measurements_mean[i, :])
        -> model: f(t) = A*exp(-lambda*t) + C
        -> LsqFit.curve_fit with lower=[..., 0.0, ...] bounds
        -> returns (amplitude, decay_rate, offset, converged,
                    residual_norm, confidence_interval, standard_error)
    |
    v
Select best fit (lowest residual_norm, decay_rate > 0, converged)
    |
    v
Optional: cross_validate_gap(estimated, exact_result::LindbladianResult)
    -> rel_error = |estimated_gap - |Re(exact_gap)|| / |Re(exact_gap)|
    -> Warns if Im(exact_gap) is significant
    |
    v
Return SpectralGapResult{Float64}(
    estimated_gap, gap_ci, gap_se, amplitude, offset, converged,
    residual_norm, observable_name, observable_index,
    all_gap_estimates, all_observable_names,
    exact_gap, relative_error,
    n_trajectories, total_time, save_every, seed
)
```

### Cross-Validation Flow (n=4,6 only)

```
# Step 1: Compute exact gap via full Liouvillian
liouv_config = LiouvConfig(... same params as ThermalizeConfig ...)
exact_result = run_lindbladian(jumps, liouv_config, ham)
# exact_result.spectral_gap is the exact gap (Complex, Re part is decay rate)

# Step 2: Estimate gap from trajectories
gap_result = estimate_spectral_gap(jumps, config, psi0, ham;
    ntraj=5000, exact_result=exact_result)

# Step 3: Compare
@info "Gap estimation" exact=abs(real(exact_result.spectral_gap))
    estimated=gap_result.estimated_gap
    relative_error=gap_result.relative_error
    ci_95=gap_result.gap_confidence_interval
```

## Key Design Decisions

### Decision 1: Separate File vs Extending trajectories.jl

**Choice:** New file `src/spectral_gap.jl`.

**Rationale:** `trajectories.jl` is already 927 lines. Adding ~200 lines of gap estimation code would push it past 1100 lines. The gap estimation logic is a distinct concern (fitting, cross-validation) that happens to use trajectory primitives. Separation keeps both files focused. The new file imports from trajectories via the module (no circular dependencies).

### Decision 2: `_run_chunk_obs_only!` vs Modifying `_run_chunk_with_obs!`

**Choice:** New function `_run_chunk_obs_only!`.

**Rationale:** The existing `_run_chunk_with_obs!` works correctly for its use case (observable measurement + DM accumulation). Adding a boolean flag `accumulate_dm` to it would be a minor change, but it introduces a branch into a per-trajectory hot loop and complicates the function's contract. A separate function is 25 lines and crystal clear in purpose. The DRY cost is minimal -- the two functions share the same structure but differ in one line (the `_accumulate_density_matrix!` call).

**Note:** The new `_run_chunk_obs_only!` DOES accept an `accumulate_final_dm` flag, but this controls whether DM is accumulated at the END of each trajectory (not at every step). This is distinct from the mid-simulation DM reconstruction that the convergence pipeline does. The flag defaults to `false` for pure observable runs, and `true` when the caller wants a DM estimate alongside the observable time series.

### Decision 3: LsqFit.jl for Exponential Fitting

**Choice:** LsqFit.jl (new dependency).

**Rationale:** LsqFit.jl provides the complete curve fitting pipeline needed for spectral gap estimation:

1. **Parameter bounds** -- `curve_fit(model, t, y, p0; lower=[...], upper=[...])` constrains gap > 0 natively. Essential because unconstrained optimization can converge to negative decay rates on noisy data.
2. **Confidence intervals** -- `confidence_interval(fit, alpha)` computes t-distribution-based CIs from the Jacobian covariance matrix. No separate statistics package needed.
3. **Standard errors** -- `standard_error(fit)` gives Jacobian-based parameter uncertainties.
4. **Weighted fitting** -- `wt = 1 ./ variance` for non-uniform noise across the time series.
5. **Automatic Jacobian** -- ForwardDiff computes the exact Jacobian, enabling proper Levenberg-Marquardt and covariance estimation.

The alternative (Optim.jl, already a dependency) would require manually implementing all of the above. While Optim.jl's NelderMead can minimize the residual, it cannot produce confidence intervals or standard errors without additional work. For a curve fitting task, the purpose-built tool is the right choice.

The dependency cost is minimal: LsqFit.jl depends on ForwardDiff, NLSolversBase, StatsAPI (all already in Manifest via Optim) plus Distributions.jl (the only truly new transitive package, used internally for t-distribution quantiles).

### Decision 4: Observable Selection for Gap Estimation

**Choice:** Energy (H) + Total magnetization (M_z) + ZZ correlations as defaults, with user-override.

**Rationale:** The spectral gap governs the slowest-decaying mode. Energy `<H>` has the strongest coupling to the gap mode (it is diagonal in the eigenbasis where the gap mode describes population redistribution). ZZ correlations `<Z_iZ_{i+1}>` probe two-body terms that make up the Hamiltonian. Total magnetization `M_z` is included as a consistency check but may have weak gap coupling for near-isotropic systems (see PITFALLS.md Pitfall 10). Fitting multiple observables and selecting the best provides robustness.

### Decision 5: SpectralGapResult Location

**Choice:** Define `SpectralGapResult` in `src/spectral_gap.jl` (not `src/structs.jl`).

**Rationale:** Unlike `LindbladianResult` and `TrajectoryResult` which are used across multiple files, `SpectralGapResult` is produced and consumed only by the gap estimation module. Co-locating the struct with its producer keeps `structs.jl` focused on core types and makes the gap estimation module self-contained. If BSON serialization is needed later, add a converter in `results.jl`.

## Patterns to Follow

### Pattern 1: Follow the `run_trajectories_convergence` Template

`run_trajectories_convergence` in `convergence.jl` demonstrates the correct pattern for building a higher-level trajectory runner on top of existing primitives:

1. Call `_build_framework_and_seed` once
2. Use `_run_batch_no_obs!` for trajectory execution
3. Accumulate results across batches
4. Return a result struct

The new `run_observable_trajectories` follows the same pattern but uses `_run_chunk_obs_only!` (observable variant) instead of `_run_batch_no_obs!`.

### Pattern 2: Multi-Threading via Existing Chunk Pattern

The existing codebase uses `_partition_trajectories` to split work across threads, and `@sync for (idx, chunk) in enumerate(chunks); Threads.@spawn ...` for parallel execution. The new observable-only runner reuses this exact pattern, including the `BLAS.set_num_threads(1)` guard.

### Pattern 3: Observable Builders in convergence.jl

The existing `build_convergence_observables` builds observables in the eigenbasis by: (1) constructing the operator in computational basis, (2) transforming to eigenbasis via `V' * O * V`. The new `build_total_magnetization` follows the same pattern.

## Anti-Patterns to Avoid

### Anti-Pattern 1: Fitting to Noisy Plateau Region

**What goes wrong:** If `total_time` is much longer than 1/gap, the observable time series is mostly flat (at the steady-state value) with noise. Fitting an exponential to mostly-flat data gives unreliable lambda estimates.

**Prevention:** The `skip_initial` parameter trims the early transient. More importantly, `total_time` should be chosen to be ~3-5x the expected mixing time (3-5/gap), so the exponential decay is well-sampled. The cross-validation against exact eigenvalues will catch bad fits.

### Anti-Pattern 2: Selecting the Wrong Observable

**What goes wrong:** If the observable chosen for fitting is orthogonal to the gap mode, its projection onto the slowest-decaying mode is zero. The fit then captures the second-slowest mode (or noise), giving the wrong gap estimate.

**Prevention:** Fit multiple observables and select the one with the best fit quality (lowest residual). Energy and ZZ correlations are both global observables that generically have nonzero overlap with the gap mode for Heisenberg chains.

### Anti-Pattern 3: Ignoring the Sign of the Decay Rate

**What goes wrong:** The optimizer finds a negative `lambda` (exponentially growing solution), which is unphysical.

**Prevention:** Use LsqFit.jl's `lower` bounds to constrain `p[2] >= 0`. The `estimate_spectral_gap` function additionally filters for `fit.decay_rate > 0`.

### Anti-Pattern 4: Using Too Few Trajectories

**What goes wrong:** With N=100 trajectories, the observable time series has O(1/sqrt(100)) = 10% noise per time point. The exponential fit is dominated by noise rather than signal.

**Prevention:** Use at least N=1000 trajectories for n=4, and N=5000+ for n=8. The noise level scales as 1/sqrt(N), so 10000 trajectories give ~1% noise.

### Anti-Pattern 5: Comparing abs(spectral_gap) Instead of abs(real(spectral_gap))

**What goes wrong:** The Liouvillian eigenvalue is complex. `abs(lambda) = sqrt(Re^2 + Im^2) >= abs(Re(lambda))`. Cross-validation using `abs` systematically overestimates agreement.

**Prevention:** Always use `abs(real(exact_result.spectral_gap))` for comparison. Warn when `|Im/Re| > 0.1`.

## Scalability Considerations

| Concern | n=4 (dim=16) | n=6 (dim=64) | n=8 (dim=256) |
|---------|--------------|--------------|---------------|
| Observable measurement cost per step | ~1us (2+ obs, dim=16 gemv) | ~10us (2+ obs, dim=64 gemv) | ~200us (2+ obs, dim=256 gemv) |
| DM accumulation per trajectory | ~0.5us (rank-1 update, 16x16) | ~8us (64x64) | ~130us (256x256) |
| Memory per thread (obs-only) | ~10 KB | ~100 KB | ~1.5 MB |
| Exact Liouvillian (cross-validation) | ~0.1s (256x256 eigs) | ~30s (4096x4096 eigs) | INFEASIBLE (65536x65536) |
| Recommended N_traj for gap fit | 1000-5000 | 2000-5000 | 5000-10000 |
| LsqFit curve_fit time | <1ms (3-param fit, ~200 points) | <1ms | <1ms |

**Key scaling insight:** For n=8, the exact Liouvillian is infeasible (dim^2 = 65536, so the Liouvillian is 65536x65536). This is precisely why trajectory-based gap estimation is needed. Cross-validation is limited to n<=6. The fitting itself is negligible cost; all computation time is in trajectory sampling.

## Suggested Build Order

### Phase 1: Observable Builders (No Dependencies)

**Add** `build_total_magnetization` and `build_gap_estimation_observables` to `convergence.jl`. Export from module.

**Test:** Verify M_z is Hermitian in eigenbasis, `tr(gibbs * M_z)` gives expected magnetization.

**Complexity:** LOW (~20 lines).

### Phase 2: Observable-Only Trajectory Runner (Depends on Phase 1)

**Add** `_run_chunk_obs_only!` and `run_observable_trajectories` to new `src/spectral_gap.jl`. This is the core trajectory loop variant.

**Test:** Run with known observables, verify `measurements_mean` matches `run_trajectories` with same observables and seed. Verify DM is zeros when `reconstruct_dm=false`.

**Complexity:** MEDIUM (~120 lines, follows existing patterns exactly).

### Phase 3: Exponential Fitting (Can Parallel with Phase 2)

**Add** LsqFit.jl dependency. Add `fit_exponential_decay` to `src/spectral_gap.jl`.

**Test:** Fit synthetic exponential data `y = 2.0 * exp(-0.5 * t) + 1.0 + noise`, verify recovered parameters match within tolerance. Verify `confidence_interval` brackets the true value.

**Complexity:** LOW-MEDIUM (~50 lines).

### Phase 4: Gap Estimation API + Result Struct (Depends on Phases 2, 3)

**Add** `SpectralGapResult`, `estimate_spectral_gap`, `cross_validate_gap` to `src/spectral_gap.jl`. Wire up the full pipeline.

**Test:** For n=4 (where exact gap is cheap to compute), verify `relative_error < 0.2` (20% agreement). This is a loose threshold for noisy trajectory data; tighter thresholds are a tuning problem.

**Complexity:** MEDIUM (~80 lines).

### Phase 5: Cross-Validation Experiment (Depends on Phase 4)

**Add** a simulation script (e.g., `simulations/gap_estimation.jl`) that:
1. For n=4 and n=6: compute exact gap via `run_lindbladian`, estimate gap from trajectories, report relative error
2. For n=8: estimate gap from trajectories only (exact is infeasible)
3. Sweep over beta values to characterize how gap estimation quality varies with temperature

**Complexity:** LOW (script, not library code).

### Build Order Rationale

- **Phase 1 first:** Observable builders are pure functions with no dependencies. Fast to implement and test. Enable Phase 2.
- **Phase 2 and 3 can be parallel:** The trajectory runner and the fitting function are independent. However, Phase 2 depends on Phase 1 for observable construction.
- **Phase 4 integrates 2+3:** The top-level API wires the trajectory runner to the fitter. Must come after both.
- **Phase 5 is validation:** Uses everything, validates the approach, produces the deliverable.

## Sources

### Primary (HIGH confidence -- direct codebase analysis)
- `src/trajectories.jl` (927 lines) -- `_run_chunk_with_obs!`, `_run_chunk_no_obs!`, `_run_batch_no_obs!`, `run_trajectories`, `_build_framework_and_seed`, `step_along_trajectory!`, `_accumulate_measurements!`, `_accumulate_density_matrix!`, `_partition_trajectories`
- `src/convergence.jl` (387 lines) -- `run_trajectories_convergence`, `run_trajectories_adaptive`, `build_convergence_observables`, `build_convergence_observables_trotter`
- `src/furnace.jl` (163 lines) -- `run_lindbladian` returning `LindbladianResult` with `spectral_gap`
- `src/structs.jl` (358 lines) -- `LindbladianResult`, `TrajectoryResult`, `ConvergenceData`, domain types, config hierarchy
- `src/QuantumFurnace.jl` -- module structure, exports, `using Optim`
- `Project.toml` -- Optim.jl already a dependency; LsqFit.jl to be added

### Secondary (HIGH confidence -- verified via official documentation)
- [LsqFit.jl Documentation](https://julianlsolvers.github.io/LsqFit.jl/latest/) -- curve_fit API, confidence_interval, standard_error, parameter bounds
- [LsqFit.jl GitHub v0.15.1](https://github.com/JuliaNLSolvers/LsqFit.jl) -- Dependencies, Julia compat, release history

### Tertiary (MEDIUM confidence -- methodology references)
- [Excitation Gap from Optimized Correlation Functions in QMC](https://ar5iv.labs.arxiv.org/html/1112.2269) -- methodology for gap estimation from Monte Carlo data
- [Generalized Moment Method for Gap Estimation](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.115.080601) -- alternative gap estimation approaches
- Lindbladian spectral gap theory: observables decay exponentially at rate governed by the spectral gap
- Chen, Kastoryano, Gilyen (2025) "An efficient and exact noncommutative quantum Gibbs sampler" -- the KMS construction that QuantumFurnace implements

---
*Architecture research for: v1.3 Spectral gap estimation from trajectory observable decay*
*Researched: 2026-02-16*
