# Project Research Summary

**Project:** QuantumFurnace.jl v1.2 Multi-Threaded Trajectory Engine & KMS-vs-GNS Experiments
**Domain:** High-performance quantum Monte Carlo trajectory sampling for Gibbs state preparation
**Researched:** 2026-02-15
**Confidence:** HIGH

## Executive Summary

QuantumFurnace.jl implements quantum trajectory sampling for Gibbs state preparation using the Monte Carlo wave function (MCWF) method with Chen et al.'s KMS and GNS detailed balance constructions. The v1.2 milestone adds multi-threaded trajectory parallelism, the GNS trajectory path (approximate detailed balance), adaptive convergence monitoring, and comparison experiments to validate that KMS (exact detailed balance with coherent correction) outperforms GNS.

The recommended approach is straightforward: use Julia's stdlib threading (`Threads.@spawn`) with per-thread workspaces and explicit RNG seeding for reproducibility. No new production dependencies are needed for threading — only JLD2 for experiment data serialization and optionally DataFrames/CSV for parameter sweep tables. The architecture separates the immutable read-only precomputed data (NUFFT prefactors, Kraus matrices) from mutable per-thread scratch buffers (TrajectoryWorkspace), enabling embarrassingly parallel execution across thousands of independent trajectories.

The critical risk is BLAS thread oversubscription: Julia threads calling OpenBLAS-threaded `mul!` creates nested parallelism that destroys performance. The fix is trivial (`BLAS.set_num_threads(1)` before spawning threads) but must be applied from the start. Secondary risks include false RNG seeding (breaking reproducibility), data races from shared mutable workspace, and conflating fixed-point approximation error with mixing rate when comparing KMS vs GNS.

## Key Findings

### Recommended Stack

Julia's standard library threading is sufficient for this milestone. The trajectory sampling is embarrassingly parallel: each trajectory is an independent sequence of CPTP map applications with its own state vector, workspace buffers, and RNG. Shared read-only data (precomputed NUFFT prefactors, per-operator Kraus matrices R/K0/U_residual/U_B) is safely accessed by all threads via a single `TrajectoryFramework` instance.

**Core technologies:**
- **Base.Threads (stdlib)** — `@spawn`/`@sync` pattern for parallel trajectories. TaskLocalRNG provides per-task thread-safe random streams. Already imported in the codebase.
- **JLD2.jl (new dependency)** — Experiment result serialization. Preserves Julia types exactly (nested NamedTuples, ComplexF64 matrices, parametric structs). HDF5-compatible for Python/MATLAB interop. Replaces BSON for new data.
- **DataFrames.jl + CSV.jl (optional)** — Parameter sweep result tables if sweep API lives in the package. Use JLD2 for full results + CSV for human-readable summary tables.
- **BLAS.set_num_threads(1)** — Critical for avoiding thread oversubscription. At dim=256 (n=8), single-threaded BLAS per trajectory + Julia thread-level parallelism across trajectories is optimal.

**No new packages needed for:**
- Adaptive convergence: Welford's online mean/variance algorithm is 10 lines of code. OnlineStats.jl is overkill.
- Thread-safe RNG: Julia 1.7+ TaskLocalRNG is built-in. Use explicit `Xoshiro(seed + trajectory_id)` for reproducibility.
- Observable tracking: Existing `_accumulate_measurements!` pattern extends naturally to per-thread buffers.

### Expected Features

**Must have (table stakes for v1.2):**
- **Multi-threaded trajectory sampling** — Thousands of trajectories at n=8 require parallelism. Each trajectory is independent. Share precomputed data (R, K0, NUFFT prefactors) read-only; clone workspace per thread.
- **Per-thread workspace** — TrajectoryWorkspace (jump_oft, psi_tmp, Rpsi buffers) is mutable. Cannot be shared. Each thread allocates its own.
- **GNS trajectory path** — ThermalizeConfigGNS already exists. Dispatches correctly through `pick_transition` to `_pick_transition_gns` (unshifted gamma). No coherent B term (`with_coherent=false` enforced). Verify correct integration.
- **Density matrix accumulation** — Thread-safe: each thread accumulates `rho_local`, final reduction merges all. Avoids contention.
- **Trace distance to Gibbs tracking** — Primary convergence metric. Compute `trace_distance_h(rho_avg, gibbs)` at batch intervals.
- **Per-observable convergence** — Track `<Z_iZ_{i+1}>` (nearest-neighbor correlations), `<Z_i>` (local magnetization), `<H>` (energy) per trajectory batch. Most informative for Heisenberg chains.
- **Adaptive sampling** — Run trajectory batches until convergence criterion met (trace distance stabilized, observable standard error below threshold). Maximum budget cap prevents infinite loops.
- **Experiment data architecture** — Serialize ExperimentResult (config + rho_mean + convergence curves + metadata) via JLD2. Enables parameter sweeps over (n, beta, KMS/GNS).
- **KMS-vs-GNS experiment driver** — Script running matched experiments: same Hamiltonian, same beta, same delta, different detailed balance construction. Compares final trace distances and convergence rates.

**Should have (differentiators):**
- **Convergence curve plotting** — Paper-ready plots of trace distance vs N_traj for KMS and GNS on same axes.
- **Bootstrap confidence intervals** — Error bars on trace distance from trajectory sub-batches.
- **Spectral gap measurement** — For n<=6, compute Liouvillian gap. For n=8, estimate from trajectory convergence rate.
- **Multiple initial states** — Verify convergence from maximally mixed and random pure states.

**Defer (v2+ or anti-features):**
- **GPU acceleration** — dim=256 is far too small for GPU advantage. GPU kernel launch overhead dominates.
- **Distributed (MPI) trajectories** — Single-node multi-core sufficient for N_traj ~ 10^4-10^5. MPI adds serialization overhead.
- **Float32 trajectories** — Trace distance to Gibbs at high beta is ~1e-6, below Float32 noise floor.
- **Adaptive timestep** — QuantumFurnace uses discrete-time CPTP unraveling, not continuous-time MCWF. Timestep is fixed.

### Architecture Approach

The architecture cleanly separates immutable shared data from mutable per-thread state. The current single-threaded design embeds a `TrajectoryWorkspace` in `TrajectoryFramework`, which blocks sharing. The fix: pass workspace as a separate argument to `step_along_trajectory!`, allowing one framework to serve many threads.

**Major components:**

1. **TrajectoryFramework (immutable, shared read-only)** — Contains per-operator Kraus matrices (R, K0, U_residual, U_B), precomputed NUFFT prefactors, transition function, config, jumps. Built once before spawning threads. All fields are read-only during trajectory stepping.

2. **ThreadTrajectoryState (mutable, per-thread)** — Contains TrajectoryWorkspace (scratch buffers), per-thread psi vector, per-thread rho accumulator, per-thread observable accumulator, per-thread RNG. Allocated per spawned thread/task. Merged after `@sync`.

3. **AdaptiveBatchManager (convergence control)** — Runs trajectory batches of fixed size, evaluates convergence criteria (trace distance stability, observable standard error), decides continue/stop. Maintains running mean and convergence history. Maximum budget prevents infinite loops.

4. **ExperimentResult (data model)** — Contains config, rho_mean, trace distance, observable time series, convergence history, metadata (timestamp, Julia version, thread count, seed). Serialized via JLD2 per experiment point.

5. **GNS dispatch (existing, no changes)** — ThermalizeConfigGNS <: AbstractThermalizeConfig. `pick_transition(config::ThermalizeConfigGNS)` returns unshifted gamma. `with_coherent=false` enforced, so B term is skipped. All trajectory machinery already handles this via polymorphic dispatch.

**Key patterns:**
- **Thread coordination:** `@spawn` per batch or chunk, each task allocates workspace in closure, reduces locally, returns result. Main thread sums batch results. No locks or atomics needed.
- **BLAS thread management:** Set `BLAS.set_num_threads(1)` before trajectory loop, restore after. Critical for n=4-8 (dim=16-256) where BLAS threading overhead exceeds computation time.
- **Reproducible RNG:** Seed per-trajectory RNGs deterministically: `rng_i = Xoshiro(master_seed + trajectory_id)`. Pass to `step_along_trajectory!(psi, fw, rng)`. Results reproducible within Julia version at fixed thread count.

### Critical Pitfalls

1. **BLAS thread oversubscription (CRITICAL)** — OpenBLAS uses multiple threads for `mul!`. Julia threads + BLAS threads = 64 OS threads on 8 cores. Performance collapse (73% slower than serial in documented cases). **Prevention:** `BLAS.set_num_threads(1)` at entry point of parallel trajectory engine. The codebase already has a TODO comment flagging this issue.

2. **Shared mutable TrajectoryWorkspace data race (CRITICAL)** — The current TrajectoryFramework embeds a single `ws::TrajectoryWorkspace{T}` with shared buffers (jump_oft, psi_tmp, Rpsi). Naive parallelization of the trajectory loop causes all threads to write to the same buffers concurrently. Silently wrong results. **Prevention:** Separate workspace from framework. Pass `ws` as explicit argument to `step_along_trajectory!`. Each thread allocates its own workspace.

3. **Global RNG breaks reproducibility (CRITICAL)** — `step_along_trajectory!` currently calls `rand()` (implicit global RNG). TaskLocalRNG is thread-safe but NOT reproducible across different thread counts (task seeding depends on spawn order). **Prevention:** Modify `step_along_trajectory!` to accept `rng::AbstractRNG` parameter. Seed per-trajectory RNGs independently: `Xoshiro(seed + i)`.

4. **KMS-vs-GNS comparison at same sigma conflates physics (MODERATE)** — KMS uses shifted gamma `(w + beta*sigma^2/2)`, GNS uses unshifted gamma `w`. Same sigma parameter produces different Lindbladians with different spectral gaps and fixed-point errors. **Prevention:** Define comparison protocol upfront (same sigma vs matched gap vs same compute budget). Document which protocol is used. Separate fixed-point error from mixing rate in analysis.

5. **Adaptive sampling stops prematurely from autocorrelation (MODERATE)** — Consecutive trajectories from same initial state produce correlated final states if mixing time is not much longer than total evolution time. Standard error estimate assumes i.i.d., underestimates uncertainty. **Prevention:** Use observable-based convergence (track `||rho(t) - rho(t-Delta)||`) rather than pure statistical criterion. Batch means with proper batch size (20-30 batches minimum).

## Implications for Roadmap

Based on research, the milestone naturally decomposes into sequential phases due to strong dependencies. The workspace refactor is a gate for all parallelism. GNS path verification and experiment data model can proceed in parallel after that. Adaptive sampling requires convergence tracking. Experiments integrate everything.

### Phase 1: Workspace Refactor (Gate for Parallelism)
**Rationale:** The embedded mutable workspace in TrajectoryFramework is the single blocker for thread safety. This refactor enables all subsequent parallel work.
**Delivers:** `step_along_trajectory!(psi, fw, ws, rng)` signature with explicit workspace and RNG arguments. Backward compatibility wrapper for existing code.
**Addresses:** Pitfall 2 (shared mutable workspace), Pitfall 3 (global RNG). Prerequisite for multi-threading.
**Complexity:** LOW — signature change to 2-3 functions, add workspace parameter.

### Phase 2: Multi-Threaded Trajectory Engine
**Rationale:** Thousands of trajectories at n=8 require parallelism. Embarrassingly parallel once workspace is separated.
**Delivers:** `run_trajectories_parallel(jumps, config, psi0, ham; ntraj, seed, nthreads)` using `@spawn`-per-batch pattern. Per-thread workspace allocation, per-thread partial rho accumulation, final reduction. BLAS thread management (`BLAS.set_num_threads(1)`).
**Addresses:** Core feature (multi-threaded sampling), Pitfall 1 (BLAS oversubscription), Pitfall 4 (false sharing via per-thread buffers), Pitfall 9 (GC pressure via allocation audit).
**Complexity:** MEDIUM — threading coordination, per-thread state management, RNG reproducibility testing.
**Depends on:** Phase 1 (workspace refactor).

### Phase 3: GNS Trajectory Path Verification
**Rationale:** GNS dispatch already exists but needs integration testing. Can proceed in parallel with Phase 2 (no code changes needed, only verification).
**Delivers:** Tests confirming GNS trajectory path produces correct fixed point. Verification that `per_op.U_B === nothing` for GNS configs. Documentation of sigma parameter space for paper experiments.
**Addresses:** Core feature (GNS path), prerequisite for KMS-vs-GNS experiments.
**Complexity:** LOW — test harness, no source changes.
**Depends on:** Nothing (GNS dispatch already works).

### Phase 4: Experiment Data Model + Serialization
**Rationale:** Need concrete data structures and persistence before building convergence tracking or experiments. JLD2 preserves Julia types, enabling clean round-trip of nested configs and matrices.
**Delivers:** `ExperimentResult` struct, `ExperimentSpec` for parameter grid, `save_experiment`/`load_experiment` via JLD2. Directory structure for experiment organization.
**Addresses:** Core feature (data architecture). Enables Phase 5 (adaptive sampling) and Phase 6 (experiments).
**Complexity:** MEDIUM — data model design, JLD2 integration, file organization.
**Depends on:** Nothing (pure data types).

### Phase 5: Convergence Tracking & Adaptive Sampling
**Rationale:** Adaptive sampling requires trace distance tracking, per-observable tracking, and batch-based convergence evaluation. Builds on multi-threading and data model.
**Delivers:** `AdaptiveBatchManager` with convergence criteria (trace distance stability, observable variance, Gibbs proximity). Batch execution with convergence evaluation. Welford online mean/variance for numerical stability.
**Addresses:** Core features (trace distance tracking, per-observable convergence, adaptive sampling). Pitfall 7 (premature stopping via observable-based convergence), Pitfall 11 (memory via scalar metrics), Pitfall 12 (batch size effects).
**Complexity:** MEDIUM — convergence criteria design, batch coordination, statistical methodology.
**Depends on:** Phase 2 (multi-threading), Phase 4 (data model).

### Phase 6: KMS-vs-GNS Experiment Driver
**Rationale:** Integration point tying all features together. Runs matched experiments across parameter grid (n=4,6,8; beta=5,10,20; KMS vs GNS). Produces paper-ready comparison data.
**Delivers:** `kms_vs_gns_grid()` parameter sweep specification. `run_experiment_grid()` driver script. Convergence curve plotting. Matched comparison protocol definition (same sigma, same initial state, separate fixed-point error from mixing rate).
**Addresses:** Core feature (KMS-vs-GNS experiments). Pitfall 5 (comparison protocol), Pitfall 6 (fixed-point vs mixing), Pitfall 15 (initial state control).
**Complexity:** MEDIUM — parameter sweep logic, experiment orchestration, result aggregation.
**Depends on:** All previous phases.

### Phase Ordering Rationale

- **Workspace refactor gates everything:** Cannot parallelize until workspace is separated from framework. Small change, large enabling effect. Do first.
- **GNS verification is independent:** Dispatches through existing polymorphic code. Can happen in parallel with Phase 2. No changes needed, only tests.
- **Data model before adaptive sampling:** Need concrete result types before building convergence manager. Data model is pure (no behavioral logic), can define early.
- **Adaptive sampling integrates threading + tracking:** Requires multi-threaded trajectory execution and data persistence. Must come after Phase 2 and Phase 4.
- **Experiments integrate everything:** Cannot run comparison experiments until all components (threading, GNS path, data model, adaptive sampling) are complete. Natural final phase.

### Research Flags

**Phases with standard patterns (skip research-phase):**
- **Phase 1 (Workspace refactor):** Straightforward signature change. No external research needed.
- **Phase 2 (Multi-threading):** Julia threading documentation is comprehensive. Pattern is well-established in QuantumOptics.jl, ITensors.jl. BLAS thread management is documented.
- **Phase 3 (GNS verification):** Physics is known (Chen 2023 paper). Implementation already exists. Testing only.
- **Phase 4 (Data model):** Standard serialization problem. JLD2 docs are sufficient.

**Phases needing validation during execution:**
- **Phase 5 (Adaptive sampling):** Convergence criteria tuning is empirical. Threshold values (rtol=0.01, batch_size=100) are reasonable defaults but may need adjustment based on pilot runs at n=4,6,8. Not a research issue — a parameter tuning issue.
- **Phase 6 (Experiments):** Sigma parameter space for fair KMS-vs-GNS comparison needs empirical validation. Initial proposal (sigma=1/beta for both) is defensible but may need refinement based on spectral gap measurements. This is experimental science, not software architecture research.

**No phases require `/gsd:research-phase`:** The architecture, stack, and pitfalls are well-understood from this research. Phase planning should focus on implementation details, not additional domain research.

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | Julia stdlib threading is verified sufficient. JLD2 is actively maintained, pure Julia. No speculative dependencies. |
| Features | HIGH | Codebase analysis shows GNS dispatch already works. Multi-threading pattern is standard MCWF methodology. Observable selection grounded in Heisenberg chain physics. |
| Architecture | HIGH | Direct source code analysis of all 23 src files. Workspace separation pattern is proven in Julia HPC community. Thread coordination matches established Monte Carlo patterns. |
| Pitfalls | HIGH | BLAS oversubscription, RNG reproducibility, and workspace data races are documented Julia issues with known solutions. KMS-vs-GNS comparison subtleties grounded in Chen et al. papers. |

**Overall confidence:** HIGH

### Gaps to Address

**Empirical validation needed:**
- **BLAS thread count crossover at n=12:** For dim=4096, benchmark `BLAS.set_num_threads(1)` with `julia -t 64` vs `BLAS.set_num_threads(4)` with `julia -t 16`. Optimal split is system-dependent. Current recommendation (BLAS=1 for n<=8) is safe but may be suboptimal at n=12.
- **Sigma parameter space for GNS:** The approximation error `||rho_fixedpoint_GNS - gibbs||` vs sigma needs empirical measurement for the specific Heisenberg Hamiltonians used. Chen 2023 theory gives asymptotic bounds but practical values depend on the spectrum.
- **Adaptive convergence threshold tuning:** The proposed rtol=0.01 for trace distance stability is reasonable but may need adjustment based on pilot runs. Too tight wastes computation; too loose gives false convergence.

**Implementation decisions deferred to planning:**
- **DataFrames in package vs experiment scripts:** If parameter sweep API (`sweep_kms_vs_gns(...)` returning DataFrame) lives in the package, add DataFrames to `[deps]`. If experiments are standalone scripts outside the package, keep it script-local. Decision depends on whether the package is a library (defer DataFrames) or includes experiment tooling.
- **Batch size formula:** Recommended `B = max(ntraj_min / 30, 100)` for 20-30 batches. Actual optimal value depends on convergence speed (fast convergence allows larger batches; slow convergence needs more checkpoints). Tune empirically.
- **Workspace allocation strategy:** Current recommendation is `@spawn`-per-batch with closure-captured workspace. Alternative: `@threads :static` with `threadid()` indexing into pre-allocated workspace pool. Both work; choose based on allocation profiling results.

**No major gaps:** The research is comprehensive for milestone planning. Remaining questions are tuning parameters and empirical validation during execution, not architecture uncertainties.

## Sources

### Primary (HIGH confidence)
- **QuantumFurnace.jl codebase** — Direct analysis of all 23 source files, all tests, and existing trajectory infrastructure. Verified workspace structure, GNS dispatch, and RNG usage.
- **Julia Multi-Threading Documentation** (https://docs.julialang.org/en/v1/manual/multi-threading/) — `@spawn`, `@sync`, TaskLocalRNG, thread safety.
- **Julia Random stdlib** (https://docs.julialang.org/en/v1/stdlib/Random/) — TaskLocalRNG per-task seeding, deterministic child RNG derivation.
- **JLD2.jl official documentation** (https://juliaio.github.io/JLD2.jl/stable/) — API, type preservation, HDF5 compatibility.
- **ITensors.jl Multithreading Guide** (https://itensor.github.io/ITensors.jl/dev/Multithreading.html) — BLAS.set_num_threads(1) pattern.
- **Chen, Kastoryano, Gilyen (2025)** "An efficient and exact noncommutative quantum Gibbs sampler" (arXiv:2311.09207) — KMS construction, coherent B term, exact detailed balance.
- **Chen, Kastoryano, Brandao, Gilyen (2023)** "Quantum Thermal State Preparation" (arXiv:2303.18224) — GNS construction, sigma parameter, energy uncertainty.

### Secondary (MEDIUM confidence)
- **Julia Discourse: BLAS vs Julia threads** (https://discourse.julialang.org/t/julia-threads-vs-blas-threads/8914) — Thread oversubscription problem and solutions.
- **Julia Discourse: Reproducible multithreaded Monte Carlo** (https://discourse.julialang.org/t/reproducible-multithreaded-monte-carlo-task-local-random/35269) — Per-task RNG patterns.
- **Julia issue #49455** (https://github.com/JuliaLang/julia/issues/49455) — Multi-threaded `mul!` 73% slower than serial due to oversubscription.
- **Julia issue #49522** (https://github.com/JuliaLang/julia/issues/49522) — `Random.seed!` reproducibility across thread counts.
- **QuantumToolbox.jl Monte Carlo Solver** (https://qutip.org/QuantumToolbox.jl/stable/users_guide/time_evolution/mcsolve) — EnsembleThreads() pattern for parallel trajectories.
- **1D Heisenberg chain physics** — Bethe ansatz ground state correlations (~-0.44 for infinite chain), SU(2) symmetry, antiferromagnetic order. Standard condensed matter textbook knowledge.

### Tertiary (domain knowledge, established methodology)
- **Welford's online algorithm** — Numerically stable streaming variance. Used by NumPy, Rust statistical crates, OnlineStats.jl internally.
- **Monte Carlo ensemble averaging** — Standard error scales as 1/sqrt(N). Batch means with proper batch size for variance estimation.
- **MCWF trajectory methodology** — Established in quantum optics since 1990s. QuantumOptics.jl and QuTiP use same patterns.

---
*Research completed: 2026-02-15*
*Ready for roadmap: YES*
