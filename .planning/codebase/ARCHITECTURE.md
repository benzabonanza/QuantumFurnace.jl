# Architecture

**Analysis Date:** 2026-02-25

## Pattern Overview

**Overall:** Layered scientific simulation library — domain-parameterized dispatch over a physics abstraction hierarchy

**Key Characteristics:**
- All public behavior dispatches on `AbstractDomain` subtypes (`BohrDomain`, `EnergyDomain`, `TimeDomain`, `TrotterDomain`), which represent approximation levels from exact quantum to hardware-realizable circuit
- Configs carry domain as a type parameter `AbstractConfig{D<:AbstractDomain, T<:AbstractFloat}`, so the compiler selects the right algorithm path at dispatch time with zero runtime overhead
- Workspaces are pre-allocated structs — all hot-path matrix computations write into pre-allocated scratch buffers with no allocations per step

## Layers

**Configuration Layer:**
- Purpose: Encode all physical and numerical parameters for a simulation run
- Location: `src/structs.jl`
- Contains: `LiouvConfig`, `LiouvConfigGNS`, `ThermalizeConfig`, `ThermalizeConfigGNS` — all parameterized `{D<:AbstractDomain, T<:AbstractFloat}`
- Depends on: `AbstractDomain` singletons (`BohrDomain()`, etc.), physical scalars
- Used by: Everything; configs flow through every computation

**Hamiltonian Layer:**
- Purpose: Store Hamiltonians, their spectral decompositions, Bohr frequencies, Gibbs states, and Trotter approximations
- Location: `src/hamiltonian.jl`, `src/trotter_domain.jl`
- Contains: `HamHam{T}` (Hamiltonian data + eigenbasis + Gibbs), `TrottTrott{T}` (Trotter unitary eigenbasis)
- Depends on: `LinearAlgebra`, `BSON` (for loading pre-computed files from `hamiltonians/`)
- Used by: All simulation layers; Hamiltonians are constructed once and passed everywhere

**Jump Operator Layer:**
- Purpose: Represent the physical jump operators (Lindblad operators) that drive thermalization
- Location: `src/structs.jl` (`JumpOp`), `src/bohr_domain.jl`, `src/energy_domain.jl`, `src/jump_workers.jl`
- Contains: `JumpOp` struct (data in computational basis + `in_eigenbasis` for fast domain computations)
- Depends on: `HamHam` or `TrottTrott` for basis; `AbstractConfig` for domain dispatch
- Used by: `furnace.jl`, `trajectories.jl`, `krylov_workspace.jl`

**Coherent Correction Layer:**
- Purpose: Compute and precompute the coherent term B (correction that makes fixed point exactly the Gibbs state)
- Location: `src/coherent.jl`, `src/furnace_utensils.jl`
- Contains: `B_bohr`, `B_time`, `B_trotter`; `_precompute_data` family (returns domain-specific NamedTuples)
- Depends on: `HamHam`/`TrottTrott`, domain-specific OFT functions from `src/ofts.jl` and `src/nufft.jl`
- Used by: `furnace.jl` (dense Liouvillian path), `krylov_workspace.jl` (matrix-free path)

**Simulation Layer — Dense Path:**
- Purpose: Construct the full `dim^2 × dim^2` Liouvillian matrix or run step-by-step density matrix thermalization
- Location: `src/furnace.jl`, `src/jump_workers.jl`
- Contains: `run_lindbladian`, `construct_lindbladian`, `run_thermalization`; `_jump_contribution!` (accumulates dissipator in-place)
- Depends on: All layers above; uses `LindbladianWorkspace` for scratch buffers
- Used by: Callers in `simulations/`, spectral analysis via Arpack

**Simulation Layer — Matrix-Free Krylov Path:**
- Purpose: Apply Lindbladian action as a function (matvec) without building the full matrix; enables scaling to larger n
- Location: `src/krylov_workspace.jl`, `src/krylov_matvec.jl`, `src/krylov_eigsolve.jl`
- Contains: `KrylovWorkspace` (pre-allocated scratch), `apply_lindbladian!`, `apply_adjoint_lindbladian!`, `krylov_spectral_gap`
- Depends on: `KrylovKit`, `JumpOp`, `HamHam`/`TrottTrott`, precomputed effective Hamiltonians `G_left`/`G_right`
- Used by: Large-system spectral gap estimation; cross-validated against dense path

**Trajectory Layer:**
- Purpose: Monte Carlo quantum trajectory simulation (stochastic unraveling of the Lindbladian)
- Location: `src/trajectories.jl`
- Contains: `TrajectoryFramework`, `PerOperatorKraus`, `step_along_trajectory!`, `run_observable_trajectories`, `run_trajectories`
- Depends on: `JumpOp`, `HamHam`/`TrottTrott`, `AbstractThermalizeConfig`; `TrajectoryWorkspace` for scratch
- Used by: Convergence tracking, gap estimation, experiments

**Analysis Layer:**
- Purpose: Extract spectral gap and validate convergence from trajectory data
- Location: `src/gap_estimation.jl`, `src/convergence.jl`, `src/fitting.jl`, `src/diagnostics.jl`
- Contains: `estimate_spectral_gap`, `run_trajectories_convergence`, `fit_exponential_decay`, `run_exact_diagnostics`
- Depends on: Trajectory layer results; `LsqFit` for curve fitting; Arpack for dense diagnostics
- Used by: Experiment scripts in `experiments/`, `simulations/`

**Persistence Layer:**
- Purpose: Save and load experiment results to/from BSON files with metadata and provenance
- Location: `src/results.jl`
- Contains: `ExperimentResult`, `save_experiment`, `load_experiment`; Dict-based BSON serialization to handle struct evolution
- Depends on: `BSON`, `LibGit2` (embeds git hash in metadata)
- Used by: Experiment scripts; results stored in `results/`

## Data Flow

**Dense Liouvillian Path:**

1. User constructs `HamHam` (from terms/coeffs) or loads from `hamiltonians/*.bson` via `load_hamiltonian`
2. User creates `LiouvConfig` (or `LiouvConfigGNS`) with domain, beta, sigma, etc.
3. User constructs `JumpOp` vector (one per site/interaction)
4. `construct_lindbladian(jumps, config, hamiltonian)` → calls `_precompute_data(config, ham_or_trott)` → domain-specific NamedTuple
5. `_jump_contribution!(liouv_matrix, jump, ...)` accumulates each jump's dissipator into the `dim^2 × dim^2` matrix in-place via `LindbladianWorkspace`
6. Arpack `eigs()` extracts near-zero eigenvalues → `LindbladianResult` (gap, fixed point, gap mode)

**Trajectory Thermalization Path:**

1. User constructs `HamHam`, `ThermalizeConfig`, `JumpOp` vector
2. `build_trajectoryframework(jumps, ham_or_trott, config, ...)` precomputes per-operator Kraus data (`PerOperatorKraus`: `R`, `K0`, `U_residual`, `U_B`)
3. `step_along_trajectory!(ws, fw, rng)` applies a random CPTP channel step using pre-allocated `TrajectoryWorkspace`
4. Multiple trajectories run (multi-threaded via `Base.Threads`), accumulating density matrix
5. Observable expectation values recorded at save intervals → `ObservableTrajectoryResult`
6. `fit_exponential_decay` fits observable time series → `FitResult` → `SpectralGapResult`

**Matrix-Free Krylov Path:**

1. `KrylovWorkspace(config, hamiltonian, jumps)` precomputes `G_left`, `G_right` (effective Hamiltonian terms), coherent `B_total`, and per-jump eigenbases
2. `apply_lindbladian!(out_vec, in_vec, ws)` computes `L(rho)` via BLAS GEMMs on scratch — zero allocations per call
3. `KrylovKit.eigsolve` invokes the linear map repeatedly → leading eigenvalues
4. Results packaged as `KrylovGapResult`

**State Management:**
- All mutable state lives in workspace structs (`LindbladianWorkspace`, `TrajectoryWorkspace`, `KrausScratch`, `KrylovWorkspace`) that are explicitly constructed before hot loops
- Config structs and Hamiltonian structs are fully immutable (all fields computed at construction)
- `ConvergenceData` accumulates scalar metrics only (not full density matrices) to keep memory O(n_batches)

## Key Abstractions

**AbstractDomain Hierarchy:**
- Purpose: Encodes the level of approximation from exact theory to quantum circuit
- Examples: `BohrDomain` (exact Bohr frequency decomposition), `EnergyDomain` (energy integral), `TimeDomain` (Fourier of time evolution), `TrotterDomain` (Trotterized circuit)
- Pattern: Singleton structs used purely for dispatch — `config.domain isa TrotterDomain` branch selection at function entry

**AbstractConfig Hierarchy:**
- Purpose: Unified carrier for all simulation parameters, parameterized to allow dispatch on domain and precision
- Examples: `src/structs.jl` — `LiouvConfig{D,T}`, `LiouvConfigGNS{D,T}`, `ThermalizeConfig{D,T}`, `ThermalizeConfigGNS{D,T}`
- Pattern: `@kwdef struct` with `Union{T, Nothing}` fields for domain-dependent params; `validate_config!` enforces consistency

**HamHam / TrottTrott Pair:**
- Purpose: Unified "physics object" carrying both the matrix and its spectral decomposition; avoids repeated diagonalization
- Examples: `src/hamiltonian.jl`, `src/trotter_domain.jl`
- Pattern: Constructed once with full spectral data; `ham_or_trott::Union{HamHam, TrottTrott}` appears throughout as a unified argument

**Workspace Pattern:**
- Purpose: Pre-allocated scratch matrices for zero-allocation hot paths
- Examples: `LindbladianWorkspace` (`src/structs.jl`), `TrajectoryWorkspace` (`src/trajectories.jl`), `KrausScratch` (`src/kraus.jl`), `KrylovWorkspace` (`src/krylov_workspace.jl`)
- Pattern: Constructed once before a loop; passed mutably into inner functions; never allocated inside the hot path

**OFT (Optimal Fourier Transform) family:**
- Purpose: Compute the frequency-filtered jump operator `A(ω)` — the core operation of the quantum bath coupling
- Examples: `src/ofts.jl` (`oft!`, `time_oft!`, `trotter_oft!`), `src/krylov_matvec.jl` (`_krylov_oft!`), NUFFT-accelerated version in `src/nufft.jl`
- Pattern: In-place mutation of output matrix; Gaussian envelope `exp(-(ω - ν)^2 / 4σ²)` applied per Bohr frequency

## Entry Points

**Library Entry Point:**
- Location: `src/QuantumFurnace.jl`
- Triggers: `using QuantumFurnace` in any script
- Responsibilities: Declares module, all imports, all exports, and `include`s all source files in dependency order

**Dense Simulation Entry Points:**
- Location: `src/furnace.jl`
- `run_lindbladian(jumps, config, hamiltonian)` → spectral analysis via Arpack, returns `LindbladianResult`
- `construct_lindbladian(jumps, config, hamiltonian)` → raw `dim^2 × dim^2` Liouvillian matrix
- `run_thermalization(jumps, config, evolving_dm, hamiltonian)` → step-by-step DM evolution, returns `DMSimulationResult`

**Trajectory Entry Points:**
- Location: `src/trajectories.jl`
- `build_trajectoryframework(jumps, ham_or_trott, config, ...)` → `TrajectoryFramework`
- `run_trajectories(jumps, config, hamiltonian; ...)` → `TrajectoryResult`
- `run_observable_trajectories(fw, observables, ...)` → `ObservableTrajectoryResult`

**Gap Estimation Entry Points:**
- Location: `src/gap_estimation.jl`, `src/krylov_eigsolve.jl`
- `estimate_spectral_gap(fw, hamiltonian, ...)` → `SpectralGapResult` (trajectory-based)
- `krylov_spectral_gap(config, hamiltonian, jumps; ...)` → `KrylovGapResult` (matrix-free)

**Simulation Scripts:**
- Location: `simulations/main_liouv.jl`, `simulations/main_thermalize.jl`, `simulations/main_krylov_benchmark.jl`
- Triggers: Run as standalone Julia scripts with `julia --project simulations/main_*.jl`

## Error Handling

**Strategy:** Eager validation at entry points; hard `error()` calls for type mismatches; `@assert` for preconditions in internal functions; `@warn` for advisory conditions (e.g., memory estimates)

**Patterns:**
- `validate_config!(config)` called at top of `run_lindbladian` and `run_thermalization` — enforces domain-specific parameter presence
- Type mismatch between `HamHam{Th}` and config `{Tc}` raises `error("Type mismatch: ...")` immediately
- `ThermalizeConfigGNS` outer constructor raises `error("GNS configs must have with_coherent=false")` at construction
- Krylov memory guard: `_check_krylov_memory(n_qubits, krylovdim)` issues `@warn` if estimated memory > 80% of free

## Cross-Cutting Concerns

**Logging:** `@printf` / `println` directly to stdout; no logging framework; progress tracked via `ProgressMeter` in multi-trajectory runs

**Validation:** `validate_config!` in `src/furnace_utensils.jl` (or similar); domain-specific required fields checked at simulation start

**Precision:** All structs parameterized on `T<:AbstractFloat`; default is `Float64`; mixed precision errors caught early via explicit type checks

**Multi-threading:** `Base.Threads.@threads` over trajectory batches; `SharedArrays` available for distributed runs; `Distributed` used in simulation scripts for multi-process parallelism

---

*Architecture analysis: 2026-02-25*
