# Feature Research: v1.2 Multi-threaded Trajectory Engine and KMS-vs-GNS Experiments

**Domain:** Multi-threaded quantum trajectory sampling, GNS trajectory path, adaptive convergence, per-observable tracking, KMS-vs-GNS comparison experiments for quantum Gibbs samplers
**Researched:** 2026-02-15
**Confidence:** HIGH (codebase analysis + Chen 2023/2025 papers + established MCWF literature + condensed matter physics of Heisenberg chains)

## Scope

This research covers features needed for the **v1.2 Multi-threading milestone**: multi-threaded trajectory engine with shared precomputed data, GNS trajectory path (approximate, no B term), adaptive sampling until convergence, per-observable convergence tracking, data architecture for saving experiments, and KMS-vs-GNS comparison experiments. It builds on the validated trajectory infrastructure from v1.0 and the cleaned codebase from v1.1.

---

## Feature Landscape

### Table Stakes (Must Have for Milestone Completion)

Features without which the milestone goals cannot be met.

| Feature | Why Expected | Complexity | Dependencies |
|---------|--------------|------------|--------------|
| **Multi-threaded trajectory sampling** | Thousands of trajectories needed for statistics at n=8 (dim=256); serial execution impractical. Each trajectory is independent -- embarrassingly parallel. Must share precomputed data (R, K0, U_res, NUFFT prefactors) read-only across threads. | MEDIUM | Existing `run_trajectories`, `TrajectoryFramework`, `TrajectoryWorkspace` |
| **Per-thread TrajectoryWorkspace** | The `TrajectoryWorkspace` contains mutable scratch buffers (`jump_oft`, `psi_tmp`, `Rpsi`). Cannot be shared across threads. Each thread needs its own workspace; the `TrajectoryFramework` minus workspace can be shared read-only. | LOW | `TrajectoryWorkspace` struct refactor |
| **Per-thread RNG** | `rand()` in `step_along_trajectory!` must use thread-local RNG to avoid contention and ensure reproducibility. Julia's default `TaskLocalRNG` is already thread-safe since Julia 1.7, but explicit seeding per thread is needed for reproducibility. | LOW | Julia Base.Threads, StableRNGs or TaskLocalRNG |
| **GNS trajectory path** | GNS (approximate detailed balance, no B term) trajectories are the comparison baseline for the paper. The `ThermalizeConfigGNS` struct already exists and `pick_transition` dispatches correctly, but `build_trajectoryframework` and `step_along_trajectory!` currently only handle `ThermalizeConfig` (KMS). Must extend to `ThermalizeConfigGNS`. | MEDIUM | `ThermalizeConfigGNS`, existing trajectory machinery, `_pick_transition_gns` |
| **Density matrix accumulation from trajectories** | Average rho = (1/N) sum |psi><psi| across all trajectories. Already exists in `run_trajectories` (serial). Must be thread-safe: each thread accumulates locally, then reduce. | LOW | Multi-threaded trajectory sampling |
| **Trace distance to Gibbs tracking** | The primary convergence metric: `trace_distance_h(rho_avg, gibbs)` at periodic intervals during the sampling. Existing `trace_distance_h` works on Hermitian matrices. Need to compute the running average and measure periodically. | LOW | `trace_distance_h`, density matrix accumulation |
| **Per-observable convergence tracking** | Track expectation values `<O_i>` from trajectory averages. Essential for the paper: shows which physical quantities converge fastest, and where KMS advantage is most visible. `run_trajectories` already has `observables` parameter and `_accumulate_measurements!` -- needs multi-threaded extension and batch-level tracking. | MEDIUM | Observable definitions, multi-threaded trajectory, density matrix accumulation |
| **Adaptive sampling (run batches until converged)** | Run trajectory batches until a convergence criterion is met, rather than pre-specifying N_traj. Criterion: standard error of the mean trace distance (or observable) falls below a threshold. | MEDIUM | Trace distance tracking, per-observable tracking |
| **Data architecture for experiment results** | Save convergence curves, configs, final density matrices, per-observable time series. Must support parameter sweeps (beta, n, KMS vs GNS). Serialization format: BSON (already a dependency) or JLD2. | MEDIUM | All convergence tracking features |
| **KMS-vs-GNS experiment driver** | Script that runs matched experiments: same Hamiltonian, same beta, same delta, same N_traj -- one with KMS (with_coherent=true, ThermalizeConfig) and one with GNS (ThermalizeConfigGNS). Compares final trace distances, convergence rates, per-observable errors. | LOW-MEDIUM | All above features, GNS trajectory path |

### Differentiators (Valuable but Not Blocking)

Features that strengthen the paper results but are not strictly required for milestone sign-off.

| Feature | Value Proposition | Complexity | Dependencies |
|---------|-------------------|------------|--------------|
| **Convergence curve plotting** | Paper-ready convergence plots: trace distance vs. N_traj for KMS and GNS on same axes. Shows the advantage visually. | LOW-MEDIUM | Data architecture, Plots.jl or Makie.jl |
| **Bootstrap confidence intervals on trace distance** | Report error bars on convergence curves. Standard bootstrap from trajectory sub-batches. Makes paper results quantitative rather than visual. | LOW-MEDIUM | Multi-threaded trajectory, batch accumulation |
| **Hamiltonian simulation cost accounting** | Count total Hamiltonian simulation time (sum of all t0 * num_trotter_steps used) per Gibbs sample. The paper theory predicts O_tilde(beta) cost per unit Lindbladian evolution; empirical verification strengthens the paper. | MEDIUM | Instrumentation in step_along_trajectory! |
| **Spectral gap measurement for mixing time** | Compute the Liouvillian spectral gap for each (n, beta, KMS/GNS) configuration. Predicts mixing time = O(1/gap). For n=4,6 this is feasible via full Liouvillian; for n=8, only trajectory-based mixing time estimation works. | MEDIUM-HIGH | `run_lindbladian` (n<=6), trajectory convergence rate fitting (n=8) |
| **Multiple initial states** | Run from maximally mixed state (standard) and from random pure states. Convergence from different initial states strengthens claims about mixing. | LOW | Initial state selection in experiment driver |

### Anti-Features (Do Not Build for This Milestone)

| Feature | Why Requested | Why Problematic | Alternative |
|---------|---------------|-----------------|-------------|
| **GPU acceleration for trajectories** | "Matrix-vector products would be faster on GPU for dim=256" | dim=256 is far too small for GPU advantage; GPU kernel launch overhead dominates. GPU is only worthwhile at dim > ~4096 (n >= 12). All current targets are n <= 8. | Use multi-threaded CPU. Profile first; BLAS threading with MKL/OpenBLAS handles dim=256 matrix-vector products in microseconds. |
| **Distributed (MPI) trajectory sampling** | "Use cluster nodes for more parallelism" | Single-node multi-core is sufficient for N_traj ~ 10^4-10^5. MPI adds serialization overhead for precomputed data. The precomputed data (NUFFT prefactors, per-operator Kraus) is read-only and fits in single-node memory. | Use `Threads.@threads` or `Threads.@spawn` on a single node with 16-64 cores. For parameter sweeps, launch independent Julia processes per (n, beta) point. |
| **Automatic Hamiltonian generation for n=6,8** | "Generate Heisenberg Hamiltonians on the fly" | For publication-quality results, Hamiltonians should be pre-generated with `find_ideal_heisenberg` (optimized disorder) and saved via BSON. On-the-fly generation adds variability across runs and makes results non-reproducible. | Pre-generate and save n=4,6,8 Hamiltonians once. Load them for all experiments. Already done for n=3,4 in `hamiltonians/` directory. |
| **Continuous-time adaptive timestep** | "Adaptive delta for faster convergence" | QuantumFurnace implements Chen's discrete-time CPTP map unraveling, not continuous-time MCWF. Adaptive timestep would require a fundamentally different formulation. The convergence rate is determined by the Lindbladian gap, not the timestep. | Keep fixed delta. Use small enough delta (~1e-3) to keep discretization error negligible, then rely on number of steps for convergence. |
| **Float32 trajectories for speed** | "Half the memory, twice the throughput" | Float32 gives only ~7 decimal digits of precision. The PSD guard in TrajectoryFramework clamps eigenvalues near zero; at Float32 precision, floating-point noise could create negative probabilities. The trace distance to Gibbs at high beta is ~1e-6 for KMS -- below Float32 noise floor. | Keep Float64. The type parameterization from v1.1 enables Float32 paths in the future if someone needs them for very large systems, but the current paper targets do not benefit. |

---

## Feature Details

### Multi-threaded Trajectory Sampling

**Architecture:** The key insight is that trajectories are embarrassingly parallel. Each trajectory is a sequence of independent `step_along_trajectory!` calls on a private `psi` vector, using a private `TrajectoryWorkspace` for scratch buffers and a private RNG.

The shared read-only data is:
- `TrajectoryFramework` fields: `per_operator` (Vector of PerOperatorKraus -- R, K0, U_residual, U_B matrices), `jumps`, `precomputed_data` (NUFFT prefactors, transition function, energy labels), `delta`, `delta_eff`, `alpha`, `n_jumps`
- These are all immutable or read-only during trajectory stepping

The per-thread mutable data is:
- `psi`: Vector{ComplexF64} (dim,) -- the evolving state
- `TrajectoryWorkspace`: `jump_oft` (dim x dim), `psi_tmp` (dim,), `Rpsi` (dim,)
- RNG state

**Implementation pattern:**
```julia
# Build framework once (shared across threads)
fw = build_trajectoryframework(...)

# Each thread gets its own workspace and psi
Threads.@threads for traj_id in 1:ntraj
    ws = TrajectoryWorkspace(ComplexF64, dim)  # per-thread
    psi = copy(psi0)
    rng = StableRNG(seed + traj_id)  # per-thread deterministic

    for step in 1:num_steps
        step_along_trajectory!(psi, fw, ws, rng)  # pass workspace explicitly
    end

    # Thread-local accumulation
    lock(acc_lock) do
        rho_acc .+= psi * psi'
    end
    # OR: accumulate into thread-local buffer, reduce after
end
```

**Refactoring needed:** Currently `TrajectoryFramework` embeds a single `TrajectoryWorkspace` (`fw.ws`). For multi-threading, the workspace must be passed separately or cloned per thread. Two options:
1. **Separate workspace from framework:** `step_along_trajectory!(psi, fw, ws)` -- cleanest, requires signature change
2. **Clone framework per thread:** Each thread gets a full copy -- wasteful because the immutable data (per_operator matrices) would be duplicated

Recommendation: Option 1 (separate workspace). The workspace is tiny (3 buffers totaling ~3*dim^2 + 2*dim complex numbers = ~600 KB for n=8). The per-operator Kraus data is large (12 operators x 4 matrices x dim^2 = ~18 MB for n=8) and should not be duplicated.

**BLAS threading interaction:** When using multi-threaded trajectories, BLAS should use single-threaded mode (`BLAS.set_num_threads(1)`) to avoid oversubscription. Each trajectory step involves small matrix-vector products (dim=256 for n=8) -- single-threaded BLAS is faster for these sizes anyway.

**Confidence:** HIGH. This is the standard pattern used by QuantumOptics.jl and QuTiP for parallel MCWF trajectories. The pattern is well-established and the codebase is designed for it (workspace caching pattern from v1.1).

### GNS Trajectory Path

**What needs to change:** The GNS construction differs from KMS in two ways:
1. **Transition function:** `_pick_transition_gns` returns the unshifted gamma-tilde(omega) instead of the shifted gamma(omega) = gamma-tilde(omega + beta*sigma^2/2) used by KMS
2. **No coherent B term:** GNS configs have `with_coherent=false`, so `U_B = nothing` in PerOperatorKraus

The trajectory machinery (`build_trajectoryframework`, `step_along_trajectory!`, `_precompute_R`) already handles `with_coherent=false` correctly. The only missing piece is that `build_trajectoryframework` takes `AbstractThermalizeConfig` but the concrete dispatch for `_precompute_data` and `_precompute_R` may not handle `ThermalizeConfigGNS` correctly through the type hierarchy.

**Verification needed:** Check that `ThermalizeConfigGNS` dispatches correctly through:
- `pick_transition(config::ThermalizeConfigGNS)` -> already defined, returns `_pick_transition_gns`
- `_precompute_data(config::ThermalizeConfigGNS, ham_or_trott)` -> needs to verify: `ThermalizeConfigGNS <: AbstractThermalizeConfig <: AbstractConfig`, so the existing `_precompute_data(config::AbstractConfig{D}, ...)` methods should apply
- `_precompute_R(jumps, ham_or_trott, config::AbstractThermalizeConfig{D}, ...)` -> should work via supertype dispatch

**GNS sigma parameter space for the paper:** The GNS construction from Chen 2023 uses energy uncertainty sigma_E (the `sigma` parameter in configs). The standard choice is `sigma = 1/beta` (used in all existing tests and simulations). For GNS, the approximation error to the true Gibbs state depends on sigma:
- `sigma = 1/beta`: Standard. The energy resolution sigma_E = 1/beta means the Gaussian filter has width comparable to the thermal energy scale. The GNS fixed point deviates from Gibbs by an amount controlled by beta*sigma = 1.
- `sigma = c/beta` with `c < 1`: Better Gibbs approximation (narrower energy filter), but higher Hamiltonian simulation cost (longer t0 needed for the Fourier transform). Cost scales as O(1/sigma).
- `sigma = c/beta` with `c > 1`: Worse Gibbs approximation but cheaper. Interesting for showing the cost-accuracy tradeoff.

**Interesting sigma values for GNS in the paper:**
- `sigma = 1/beta` (standard, `c=1`): The baseline comparison point
- `sigma = 0.5/beta` (`c=0.5`): Twice the simulation cost, significantly better Gibbs approximation. Shows what GNS could achieve if willing to pay more.
- `sigma = 0.8/beta` (`c=0.8`): Used in the existing `main_thermalize.jl` script. Moderate improvement.

For the KMS construction, sigma affects only the Trotter/quadrature errors (not the detailed balance error, which is zero by construction). So KMS at `sigma = 1/beta` already achieves exact detailed balance, while GNS at `sigma = 1/beta` has an intrinsic approximation gap. This is the core message of the paper.

**Confidence:** HIGH for the trajectory path implementation. The GNS configs are already well-defined in the codebase. MEDIUM for the exact sigma values that make the most compelling comparison -- this depends on empirical results.

### Observable Selection for Convergence Tracking

**The question:** For 1D Heisenberg chains (n=4,6,8) at inverse temperatures beta=5,10,20, which observables are most informative for tracking convergence to the Gibbs state?

**Analysis of candidates:**

1. **Local magnetization Z_i (single-site Pauli Z on site i):**
   - The 1D isotropic Heisenberg Hamiltonian `H = sum_i (X_iX_{i+1} + Y_iY_{i+1} + Z_iZ_{i+1})` has SU(2) symmetry. For periodic boundary conditions and no disorder field, `<Z_i>_gibbs = 0` for all sites by symmetry. This makes Z_i uninformative as a convergence metric for the clean Heisenberg chain.
   - However, with the disordering Z-field (which QuantumFurnace uses -- see `find_ideal_heisenberg`), the SU(2) symmetry is broken. Local magnetization `<Z_i>` will be nonzero and site-dependent, making it a useful observable.
   - **Verdict: Use it.** The disordered Hamiltonian breaks symmetry, so `<Z_i>` is site-specific and nonzero.

2. **Total magnetization M_z = sum_i Z_i:**
   - For the disordered Heisenberg chain, `<M_z>_gibbs` is nonzero but small (the disorder breaks symmetry but doesn't create large net magnetization).
   - Low sensitivity: being a sum of N terms, the central limit theorem means trajectory fluctuations in M_z are sqrt(N) times larger than for single-site Z_i, but the signal is only sqrt(N) times larger. Signal-to-noise is the same.
   - **Verdict: Include as a global check, not primary.**

3. **Nearest-neighbor correlations Z_iZ_{i+1}:**
   - This is the most physically informative observable for the Heisenberg chain. The ground state has strong antiferromagnetic correlations (`<Z_iZ_{i+1}>` negative for antiferromagnetic coupling). At finite temperature, these correlations weaken.
   - For the isotropic Heisenberg model with J=1, the ground state correlation `<Z_iZ_{i+1}>` is approximately -0.44 for the infinite chain (Bethe ansatz result).
   - At high beta (low temperature), the Gibbs state is dominated by the ground state, so `<Z_iZ_{i+1}>` approaches the ground state value. This provides a nontrivial target that differs from the maximally mixed value (0).
   - Convergence of `<Z_iZ_{i+1}>` from 0 (maximally mixed initial state) to the Gibbs value is a clear signal.
   - **Verdict: Primary observable. Use it.**

4. **Energy <H> = sum_i J(X_iX_{i+1} + Y_iY_{i+1} + Z_iZ_{i+1}):**
   - The average energy `<H>_gibbs` has a known value from the partition function: `<H> = -d/d(beta) ln(Z)`. For the finite chain, this is computed exactly from the eigenvalues.
   - Energy convergence is a global check -- if the trajectory-averaged energy matches the Gibbs energy, the trajectory sampling is working.
   - The energy has a clear physical interpretation and connects to thermodynamic quantities (specific heat).
   - **Verdict: Include. Good global convergence check.**

5. **Staggered magnetization M_s = sum_i (-1)^i Z_i:**
   - The staggered magnetization is the order parameter for antiferromagnetic order. For the 1D spin-1/2 Heisenberg antiferromagnet, there is NO long-range order at any temperature (Mermin-Wagner theorem in 1D). However, for finite chains, `<M_s^2>` is nonzero and tracks short-range antiferromagnetic correlations.
   - For the disordered Heisenberg chain, `<M_s>` is generally nonzero due to symmetry breaking.
   - Less directly informative than Z_iZ_{i+1} for tracking convergence.
   - **Verdict: Skip. Correlations Z_iZ_{i+1} capture the same physics more directly.**

**Recommended observable set for the paper:**

| Observable | Formula | Why | Computational Cost |
|------------|---------|-----|-------------------|
| **Nearest-neighbor Z_iZ_{i+1}** (all pairs) | `Z_i tensor Z_{i+1}` padded to n qubits | Primary convergence signal; nontrivial Gibbs value; directly probes antiferromagnetic correlations | LOW -- single matrix-vector product per observable per step |
| **Local magnetization Z_i** (all sites) | `Z_i` padded to n qubits | Site-resolved convergence; nonzero due to disorder; catches localization effects | LOW |
| **Energy <H>** | The full Hamiltonian matrix | Global convergence check; connects to thermodynamics | LOW -- already have H as a matrix |
| **Trace distance to Gibbs** | `trace_distance_h(rho_avg, gibbs)` | The definitive convergence metric; requires full density matrix | MEDIUM -- requires eigendecomposition of rho_avg - gibbs every measurement interval |

Total observables for n=8: 8 (nearest-neighbor correlations, periodic) + 8 (local Z) + 1 (energy) = 17 observables. Each requires one matrix-vector product (O(dim^2) = O(65536) operations for n=8). Negligible compared to the trajectory step cost.

**Confidence:** HIGH for the observable selection. The physics is well-understood from decades of Heisenberg chain studies.

### Gibbs State Peakedness at High Beta

**The question:** For the 1D Heisenberg chain at n=8 and beta=5,10,20, how peaked is the Gibbs state? Is beta=30 reasonable?

**Analysis:**

The Gibbs state in the energy eigenbasis is diagonal: `rho_gibbs = diag(exp(-beta * E_i) / Z)` where `E_i` are the eigenvalues and `Z = sum exp(-beta * E_i)`.

The peakedness is determined by the ratio `exp(-beta * (E_1 - E_0))` where `E_0` is the ground state energy and `E_1` is the first excited state energy. If `beta * (E_1 - E_0)` is large, the Gibbs state is heavily concentrated on the ground state.

For the 1D isotropic Heisenberg chain with n sites and periodic boundary conditions:
- The spectrum is rescaled by QuantumFurnace to fit in [0, 0.45]. The eigenvalue spread is `E_max - E_min`, mapped to approximately 0.45.
- The rescaling factor for n=8 is approximately `(E_max - E_min) / 0.45`. For the isotropic Heisenberg chain with J=1, `E_max - E_min` scales linearly with n (extensive). For n=8, `E_max - E_min ~ 8J = 8`.
- The rescaling factor is approximately `8 / 0.45 ~ 17.8`.
- The effective beta in the rescaled spectrum is `beta_eff = beta / rescaling_factor`.
- The gap `Delta = E_1 - E_0` in the rescaled spectrum corresponds to `nu_min` (the smallest Bohr frequency, which QuantumFurnace optimizes via `find_ideal_heisenberg`).

For the disordered Heisenberg chain (which breaks degeneracies), typical values are:
- n=4: `nu_min ~ 0.01-0.05` (after rescaling). With beta=10, `beta * nu_min ~ 0.1-0.5`. The Gibbs state is not very peaked.
- n=6: `nu_min ~ 0.005-0.02`. With beta=10, `beta * nu_min ~ 0.05-0.2`. Moderately peaked.
- n=8: `nu_min ~ 0.002-0.01`. With beta=20, `beta * nu_min ~ 0.04-0.2`. The ground state population is `p_0 ~ 1/(1 + exp(-beta*nu_min) + ...) ~ 1/(1 + 0.96 + ...) ~ 0.03-0.1` depending on degeneracies.

**Key consideration:** The unrescaled Heisenberg chain Hamiltonian has extensive bandwidth. The rescaling maps the spectrum to [0, 0.45], making the effective `beta_eff = beta / rescaling_factor` much smaller. For n=8:
- `rescaling_factor ~ 18` (typical)
- `beta_eff(beta=5) ~ 0.28` -- warm, Gibbs state well spread
- `beta_eff(beta=10) ~ 0.56` -- moderately peaked
- `beta_eff(beta=20) ~ 1.1` -- ground state dominant but not degenerate
- `beta_eff(beta=30) ~ 1.67` -- ground state very dominant, higher levels negligible

**Practical assessment for beta=30, n=8:**
- The ground state population will be `p_0 ~ exp(-0) / Z ~ 1 - sum(exp(-30/18 * (E_i - E_0)))`. For the rescaled gap of ~0.005, `exp(-30*0.005) ~ exp(-0.15) ~ 0.86`, so the first excited state still has ~86% of the ground state weight. With the disorder lifting degeneracies, the effective number of significantly populated states might be ~10-30 (out of 256).
- beta=30 IS reasonable. The Gibbs state is peaked but not degenerate. There are enough populated states to make the convergence problem nontrivial.
- However, beta=30 will have SLOWER mixing (longer convergence time) due to the larger spectral gap being smaller relative to the thermal energy. More trajectories and more steps will be needed.

**Recommendation for the paper:**
- **beta=5**: Easy case. Fast mixing, well-spread Gibbs state. Shows that both KMS and GNS work.
- **beta=10**: Medium case. The standard test temperature used in existing code. Good data point.
- **beta=20**: Hard case. Ground state dominant. This is where KMS advantage should be most visible -- the exact detailed balance prevents trajectory drift away from the ground state sector.
- **beta=30**: Optional stretch. If beta=20 already shows a clear advantage, beta=30 adds a data point but may require very long runs. Include if computational budget allows.

**Confidence:** MEDIUM. The exact peakedness depends on the specific disordered Hamiltonian instance (which eigenvalues the disorder lifts, how large nu_min is). The rescaling factor is the key quantity -- the estimates above use typical values but the actual number depends on the pre-generated n=8 Hamiltonian. Must generate the Hamiltonian and inspect the spectrum empirically to get exact numbers.

### Adaptive Sampling Protocol

**The question:** How should adaptive sampling work -- run batches until convergence criterion is met?

**Protocol design:**

The standard approach for Monte Carlo convergence estimation is:

1. **Batch structure:** Run trajectory batches of fixed size B (e.g., B=1000). After each batch, compute the cumulative trace distance to Gibbs and per-observable estimates.

2. **Running estimates:** Maintain:
   - `rho_running = (1/N_total) * sum_{all trajectories} |psi><psi|`
   - `obs_running[i] = (1/N_total) * sum_{all trajectories} <psi|O_i|psi>`
   - `trace_dist_running = trace_distance_h(rho_running, gibbs)`

3. **Convergence criterion for trace distance:**
   - After each batch, record `trace_dist_running`. The convergence curve should plateau.
   - Criterion: the standard deviation of the last K batch-to-batch changes in trace distance is below a threshold.
   - Alternative: the relative change in trace distance between consecutive batch completions is below epsilon (e.g., 1%): `|d_new - d_old| / d_old < 0.01` for the last 3 consecutive batches.

4. **Standard error estimation for observables:**
   - Partition trajectories into M sub-batches of size B/M.
   - Compute the observable estimate from each sub-batch independently.
   - The standard error is `std(sub_batch_estimates) / sqrt(M)`.
   - Convergence criterion: `SE / |mean| < threshold` (e.g., 1% relative error) for all tracked observables.

5. **Practical convergence criterion:**
   ```
   CONVERGED when ALL of:
     (a) N_total >= N_min (e.g., 1000) -- minimum sample size
     (b) trace_dist stabilized: relative change < 1% for last 3 batches
     (c) SE(trace_dist) < epsilon_td (e.g., 0.001)
     (d) SE(<O_i>) / |<O_i>_gibbs| < epsilon_obs (e.g., 5%) for all observables
   ```

6. **Maximum budget:** Always set a hard cap `N_max` (e.g., 100,000) to prevent infinite loops if convergence is slow.

**Implementation pattern:**
```julia
function run_adaptive_trajectories(fw, psi0, gibbs, observables;
    batch_size=1000, N_min=1000, N_max=100_000,
    td_threshold=0.001, td_rel_change=0.01, obs_rel_se=0.05,
    num_stable_batches=3)

    N_total = 0
    rho_acc = zeros(ComplexF64, dim, dim)
    td_history = Float64[]

    while N_total < N_max
        # Run a batch
        rho_batch = run_batch(fw, psi0, batch_size)
        rho_acc .+= rho_batch .* batch_size
        N_total += batch_size

        # Update running average
        rho_avg = rho_acc ./ N_total
        td = trace_distance_h(Hermitian(rho_avg), gibbs)
        push!(td_history, td)

        # Check convergence
        if N_total >= N_min && is_converged(td_history, ...)
            break
        end
    end
end
```

**Confidence:** HIGH. This is standard Monte Carlo methodology. The specific thresholds (1%, 5%) are reasonable defaults that can be tuned based on the desired paper precision.

### KMS-vs-GNS Comparison Experiment Structure

**Experiment design:**

The paper's core claim is: "KMS (exact detailed balance with coherent B term) produces trajectories that converge closer to the Gibbs state than GNS (approximate detailed balance without B)."

**Experiment matrix:**

| Dimension | n=4 (dim=16) | n=6 (dim=64) | n=8 (dim=256) |
|-----------|--------------|--------------|---------------|
| beta=5    | KMS, GNS(sigma=1/beta) | KMS, GNS(sigma=1/beta) | KMS, GNS(sigma=1/beta) |
| beta=10   | KMS, GNS(sigma=1/beta) | KMS, GNS(sigma=1/beta) | KMS, GNS(sigma=1/beta) |
| beta=20   | KMS, GNS(sigma=1/beta) | KMS, GNS(sigma=1/beta) | KMS, GNS(sigma=1/beta) |

Total: 18 experiment runs (3 n values x 3 beta values x 2 methods). Each run requires adaptive sampling until convergence.

**Additional GNS sigma sweep (for one (n, beta) point):**

To show the cost-accuracy tradeoff of GNS, run GNS at sigma = {0.5, 0.8, 1.0, 1.5, 2.0}/beta for a fixed (n, beta) point, e.g., (n=6, beta=10). This produces a Pareto front: better Gibbs approximation costs more Hamiltonian simulation time.

| sigma | sigma * beta | GNS approximation quality | Relative Hamiltonian simulation cost |
|-------|-------------|--------------------------|--------------------------------------|
| 0.5/beta | 0.5 | Very good (narrow filter) | ~2x baseline |
| 0.8/beta | 0.8 | Good | ~1.25x baseline |
| 1.0/beta | 1.0 | Standard (baseline) | 1x baseline |
| 1.5/beta | 1.5 | Moderate | ~0.67x baseline |
| 2.0/beta | 2.0 | Poor (very broad filter) | ~0.5x baseline |

The KMS result serves as the "unbeatable" reference: exact detailed balance at sigma = 1/beta cost. If GNS at sigma = 0.5/beta (twice the cost) still does not match KMS at sigma = 1/beta, that is a strong argument for KMS.

**Per-experiment data to collect:**

1. **Convergence curve:** trace distance to Gibbs vs. N_traj (log-log scale)
2. **Final trace distance:** at convergence (or at N_max)
3. **Per-observable convergence:** `<Z_iZ_{i+1}>` and `<Z_i>` vs. N_traj
4. **Error of per-observable vs. Gibbs:** `|<O>_traj - <O>_gibbs|` at convergence
5. **Config snapshot:** all parameters (beta, sigma, delta, n, domain, with_coherent)
6. **Wall clock time:** total computation time per experiment

**Expected results (predictions based on theory):**

- **At all (n, beta):** KMS trace distance at convergence should be limited only by Trotter error (~1e-6 to 1e-8 depending on delta and Trotter parameters) plus statistical noise (~1e-3 at N=10^4).
- **GNS at sigma=1/beta:** The trace distance at convergence should plateau at a value determined by the GNS approximation gap, which grows with beta. At beta=5, the gap is small (~1e-2 to 1e-3). At beta=20, the gap is larger (~1e-1 to 1e-2).
- **KMS advantage:** Most visible at high beta (beta=20), where the GNS approximation gap is largest. At low beta (beta=5), both methods may perform similarly.
- **Per-observable:** The correlations `<Z_iZ_{i+1}>` should show the KMS advantage most clearly, because these are directly related to the energy structure that the detailed balance condition controls.

**Confidence:** HIGH for the experiment design. MEDIUM for the specific numerical predictions -- these depend on the Hamiltonian spectrum and the interaction between Trotter error, quadrature error, and the GNS approximation gap.

### Data Architecture

**What to save per experiment:**

```julia
struct ExperimentResult
    # Identity
    experiment_id::String      # e.g., "kms_n8_beta10_sigma0.1"
    method::Symbol             # :KMS or :GNS
    timestamp::DateTime

    # Configuration
    config::AbstractThermalizeConfig
    num_qubits::Int
    hamiltonian_file::String   # path to BSON

    # Convergence data
    n_traj_total::Int
    trace_distances::Vector{Float64}           # td at each batch checkpoint
    n_traj_at_checkpoint::Vector{Int}           # N_traj at each checkpoint
    converged::Bool
    convergence_n_traj::Int                    # N at convergence (or N_max)

    # Observable data
    observable_names::Vector{String}
    observable_gibbs_values::Vector{Float64}    # <O_i>_gibbs
    observable_estimates::Matrix{Float64}       # obs x checkpoints
    observable_std_errors::Matrix{Float64}      # obs x checkpoints

    # Final state
    rho_final::Matrix{ComplexF64}              # trajectory-averaged density matrix
    trace_distance_final::Float64

    # Timing
    wall_time_seconds::Float64
end
```

**Serialization:** Use BSON (already a dependency) for consistency with existing Hamiltonian storage. JLD2 is an alternative with better compression, but adding a dependency for this alone is not justified.

**Directory structure:**
```
results/
  v1.2_kms_vs_gns/
    hamiltonians/
      heis_n4.bson
      heis_n6.bson
      heis_n8.bson
    experiments/
      kms_n4_beta5.bson
      gns_n4_beta5.bson
      ...
    analysis/
      convergence_plots/
      comparison_tables/
```

**Confidence:** HIGH. Standard data management pattern.

---

## Feature Dependencies

```
[Pre-generate n=6,8 Hamiltonians]
    |
    v
[Refactor TrajectoryWorkspace out of TrajectoryFramework]
    |
    v
[Multi-threaded trajectory sampling]
    |               |
    v               v
[GNS trajectory path]    [Per-observable convergence tracking]
    |               |
    v               v
[Trace distance to Gibbs tracking]
    |
    v
[Adaptive sampling (batch until converged)]
    |
    v
[Data architecture for saving results]
    |
    v
[KMS-vs-GNS experiment driver]
    |
    v
[Convergence curve plotting] (differentiator)
```

### Dependency Notes

- **Workspace refactor is the gate:** Multi-threading requires per-thread workspaces. This must be done before any parallel trajectory work.
- **GNS path and observables are independent of each other:** Both depend on multi-threading but not on each other. Can be built in parallel.
- **Adaptive sampling requires convergence tracking:** Must have trace distance and observable tracking before implementing the adaptive loop.
- **Experiment driver requires everything:** It is the integration point that ties all features together. Build last.
- **Data architecture must precede experiments:** Need the serialization format defined before running experiments that produce data.
- **Hamiltonian generation is a prerequisite:** n=6 and n=8 Hamiltonians must be generated and saved before experiments can run at those sizes. Currently only n=3 and n=4 exist in `hamiltonians/`.

---

## MVP Definition (v1.2 Milestone)

### Must Complete

- [ ] **Pre-generate n=6,8 Hamiltonians** -- Cannot run paper experiments without them. Use `find_ideal_heisenberg` with batch_size=1000 or more for good nu_min.
- [ ] **Refactor workspace out of TrajectoryFramework** -- Gate for multi-threading. Change `step_along_trajectory!` signature to accept workspace explicitly.
- [ ] **Multi-threaded trajectory sampling** -- Core performance feature. Use `Threads.@threads` with per-thread workspace and RNG.
- [ ] **GNS trajectory path** -- Extend `build_trajectoryframework` to accept `ThermalizeConfigGNS`. Verify correct dispatch through precompute and stepping.
- [ ] **Per-observable convergence tracking** -- Track `<Z_iZ_{i+1}>` and `<Z_i>` across trajectory batches.
- [ ] **Trace distance to Gibbs tracking** -- Compute running trace distance at batch checkpoints.
- [ ] **Adaptive sampling** -- Run batches until convergence criterion met.
- [ ] **Data architecture** -- Define result struct and BSON serialization.
- [ ] **KMS-vs-GNS experiment driver** -- Script for running the full experiment matrix.

### Add After Core Implementation

- [ ] **Bootstrap confidence intervals** -- Error bars on convergence curves for paper quality.
- [ ] **Convergence curve plotting** -- Paper-ready plots with KMS and GNS on same axes.
- [ ] **GNS sigma sweep** -- Multiple sigma values for the Pareto front analysis.
- [ ] **Hamiltonian simulation cost accounting** -- Track total simulation time for cost comparison.
- [ ] **Multiple initial states** -- Verify convergence from different starting points.

### Defer

- [ ] **Spectral gap measurement for n=8** -- Requires O(dim^4) = O(4 billion) Liouvillian construction. Only feasible for n<=6. For n=8, estimate mixing time from trajectory convergence curves.
- [ ] **GPU acceleration** -- Not needed at current system sizes.
- [ ] **Distributed (MPI) trajectories** -- Not needed for single-node execution.

---

## Feature Prioritization Matrix

| Feature | Paper Value | Implementation Cost | Priority |
|---------|------------|---------------------|----------|
| Pre-generate n=6,8 Hamiltonians | CRITICAL | LOW | P0 |
| Workspace refactor | CRITICAL | LOW | P0 |
| Multi-threaded trajectories | CRITICAL | MEDIUM | P1 |
| GNS trajectory path | CRITICAL | LOW-MEDIUM | P1 |
| Trace distance tracking | CRITICAL | LOW | P1 |
| Per-observable tracking | HIGH | MEDIUM | P1 |
| Adaptive sampling | HIGH | MEDIUM | P1 |
| Data architecture | HIGH | MEDIUM | P1 |
| KMS-vs-GNS experiment driver | CRITICAL | LOW-MEDIUM | P2 |
| Bootstrap confidence intervals | MEDIUM | LOW-MEDIUM | P2 |
| Convergence plotting | MEDIUM | LOW-MEDIUM | P2 |
| GNS sigma sweep | MEDIUM | LOW | P2 |
| Hamiltonian simulation cost | LOW-MEDIUM | MEDIUM | P3 |
| Multiple initial states | LOW | LOW | P3 |
| Spectral gap (n<=6) | LOW | MEDIUM-HIGH | P3 |

**Priority key:**
- P0: Blocks all other work; do first
- P1: Core milestone features
- P2: Paper-quality enhancements; add once core works
- P3: Nice to have; defer if time-constrained

---

## Physics Reference: Heisenberg Chain Observables

### Observable Expectations at Different Temperatures

For the 1D spin-1/2 Heisenberg antiferromagnet (J=1, periodic) with n sites, key observable values:

| Observable | Maximally mixed (beta=0) | Low T (beta>>1) | Notes |
|------------|--------------------------|-----------------|-------|
| `<Z_i>` (clean) | 0 | 0 | SU(2) symmetry |
| `<Z_i>` (disordered) | 0 | ~disorder strength | Symmetry broken |
| `<Z_iZ_{i+1}>` | 0 | ~-0.44 (infinite chain) | Bethe ansatz result |
| `<H>/n` | 0 (rescaled) | E_0/n | Ground state energy density |
| `<M_s^2>/n` | 1/3 | depends on n | Finite-size AFM correlations |

### Why Nearest-Neighbor Correlations Win

The convergence of `<Z_iZ_{i+1}>` from 0 (maximally mixed) to its Gibbs value (~-0.3 to -0.4 at high beta) provides:
1. **Large signal:** The change from 0 to -0.4 is O(1), easily detectable above noise.
2. **Physical meaning:** Directly probes the antiferromagnetic ordering that the Gibbs state captures.
3. **Sensitivity to detailed balance:** The correlation is determined by the low-energy sector of the spectrum, which is exactly where KMS (exact detailed balance) differs from GNS (approximate).

---

## Sources

### Verified (HIGH confidence)
- QuantumFurnace.jl codebase -- Direct analysis of `trajectories.jl`, `furnace.jl`, `structs.jl`, `energy_domain.jl`, `coherent.jl`, `qi_tools.jl`, `test_helpers.jl`, `main_thermalize.jl`, `main_liouv.jl`
- Chen, Kastoryano, Brandao, Gilyen (2023) "Quantum Thermal State Preparation" [arXiv:2303.18224](https://arxiv.org/abs/2303.18224) -- GNS construction, sigma parameter, energy uncertainty
- Chen, Kastoryano, Gilyen (2025) "An efficient and exact noncommutative quantum Gibbs sampler" [arXiv:2311.09207](https://arxiv.org/abs/2311.09207) -- KMS construction, coherent B term, exact detailed balance
- Ding, Li, Lin (2024) "Efficient quantum Gibbs samplers with KMS detailed balance" [arXiv:2404.05998](https://arxiv.org/abs/2404.05998) -- Alternative KMS construction, comparison framework
- [Julia Multi-Threading Documentation](https://docs.julialang.org/en/v1/manual/multi-threading/) -- Threads.@threads, TaskLocalRNG, thread safety
- [QuantumToolbox.jl Monte Carlo Solver](https://qutip.org/QuantumToolbox.jl/stable/users_guide/time_evolution/mcsolve) -- Parallel trajectory methodology
- [QuantumOptics.jl Quantum Trajectories](https://docs.qojulia.org/timeevolution/mcwf/) -- MCWF ensemble averaging

### Partially verified (MEDIUM confidence)
- [Abdelhafez et al. (2019)](https://arxiv.org/abs/1803.08589) -- Adaptive MCWF convergence study
- 1D Heisenberg chain ground state correlations -- Bethe ansatz results (standard condensed matter textbook knowledge)
- [Heisenberg chain Wikipedia](https://en.wikipedia.org/wiki/Quantum_Heisenberg_model) -- General properties

### Domain knowledge (HIGH confidence, established physics)
- SU(2) symmetry of isotropic Heisenberg chain implies `<Z_i> = 0` without disorder
- Nearest-neighbor correlations `<Z_iZ_{i+1}>` are the natural order parameter for antiferromagnetic Heisenberg models
- 1D spin-1/2 Heisenberg antiferromagnet is gapless (Haldane conjecture, proven by Bethe ansatz)
- Monte Carlo ensemble averaging: standard error scales as 1/sqrt(N)
- Batch-based adaptive sampling is standard methodology in computational physics

---
*Feature research for: v1.2 Multi-threaded trajectory engine and KMS-vs-GNS experiments*
*Researched: 2026-02-15*
