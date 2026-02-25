# Codebase Structure

**Analysis Date:** 2026-02-25

## Directory Layout

```
QuantumFurnace.jl/
├── src/                        # Library source code (single Julia module)
│   ├── QuantumFurnace.jl       # Module entry point: all imports, exports, includes
│   ├── structs.jl              # Core types: domains, configs, result structs, workspaces
│   ├── hamiltonian.jl          # HamHam struct: construction, spectral decomposition, Gibbs state
│   ├── trotter_domain.jl       # TrottTrott struct: Trotter unitary, quasi-Bohr frequencies
│   ├── constants.jl            # Physical/numerical constants
│   ├── errors.jl               # Error utilities (minimal placeholder)
│   ├── qi_tools.jl             # Quantum information utilities: trace distance, fidelity, Paulis
│   ├── misc_tools.jl           # Matrix utilities: Kronecker, Liouvillian vectorization, Pauli tools
│   ├── energy_domain.jl        # Energy-domain transition functions (pick_transition dispatch)
│   ├── bohr_domain.jl          # Bohr-domain coherent B term: B_bohr
│   ├── coherent.jl             # All coherent B terms: B_time, B_trotter, precompute helpers
│   ├── furnace_utensils.jl     # _precompute_data dispatch family; validate_config!
│   ├── furnace.jl              # Top-level: run_lindbladian, construct_lindbladian, run_thermalization
│   ├── ofts.jl                 # OFT functions: oft!, time_oft!, trotter_oft!
│   ├── nufft.jl                # NUFFT-accelerated OFT prefactor computation (FINUFFT wrapper)
│   ├── time_domain.jl          # Time-domain helpers
│   ├── kraus.jl                # KrausScratch workspace struct
│   ├── jump_workers.jl         # _jump_contribution! dispatch family (all domains); vectorize helpers
│   ├── trajectories.jl         # Trajectory simulation: TrajectoryFramework, step_along_trajectory!, run_*
│   ├── krylov_workspace.jl     # KrylovWorkspace: pre-allocated matrix-free matvec workspace
│   ├── krylov_matvec.jl        # apply_lindbladian!, apply_adjoint_lindbladian!, sandwich helpers
│   ├── krylov_eigsolve.jl      # krylov_spectral_gap, apply_delta_channel!, KrylovGapResult
│   ├── linearmaps_liouv.jl     # LinearMaps wrapper stubs (mostly commented out / unused)
│   ├── log_sobolev.jl          # LSI alpha2 computation via Optim
│   ├── log_sobolev_manopt.jl   # (empty stub)
│   ├── convergence.jl          # run_trajectories_convergence, run_trajectories_adaptive, build_preset_trajectory_observables
│   ├── fitting.jl              # fit_exponential_decay, FitResult
│   ├── gap_estimation.jl       # estimate_spectral_gap, SpectralGapResult, eigenbasis_overlap_analysis
│   ├── diagnostics.jl          # run_exact_diagnostics and DIAG-01..06 functions
│   └── results.jl              # ExperimentResult, save_experiment, load_experiment (BSON)
├── test/                       # Test suite
│   ├── runtests.jl             # Test runner: single @testset includes all test files
│   ├── test_helpers.jl         # Shared fixtures: Hamiltonians, configs, jump operators
│   ├── test_*.jl               # Individual test files, one per feature area
│   ├── reference/              # Reference data files for regression tests
│   ├── trajectory_validation/  # Trajectory-specific validation data and scripts
│   └── old_tests/              # Archived tests (not run in CI)
├── hamiltonians/               # Pre-computed Hamiltonian BSON files
│   ├── generate_hamiltonians.jl  # Script to regenerate BSON files
│   └── heis_disordered_periodic_n{3..10}.bson  # Heisenberg chain Hamiltonians, n=3..10
├── simulations/                # Runnable simulation scripts (not part of the library)
│   ├── main_liouv.jl           # Dense Liouvillian spectral analysis
│   ├── main_thermalize.jl      # Step-by-step DM thermalization
│   ├── main_krylov_benchmark.jl  # Krylov matvec benchmarking
│   └── run_julia.sbatch        # SLURM job submission script
├── experiments/                # Analysis scripts and stored experiment output
│   ├── validate_spectral_gap.jl
│   ├── validate_gap_delta_scaling.jl
│   ├── diagnose_gap_momentum.jl
│   ├── investigate_delta_scaling_bug.jl
│   ├── *.bson                  # Stored experiment results (KMS and GNS runs)
│   └── *.txt                   # Human-readable summaries alongside BSON files
├── results/                    # Simulation result storage (BSON + derived outputs)
├── playground/                 # One-off Julia scripts and notebooks for exploration
│   ├── *.jl                    # Ad-hoc exploration scripts (not maintained)
│   └── jl_playground.ipynb    # Jupyter notebook
├── docs/                       # Documenter.jl documentation source
│   ├── src/                    # Documentation markdown and literate sources
│   └── build/                  # Built documentation (generated, not committed)
├── tools/                      # External tooling
│   ├── python/                 # Python analysis scripts
│   └── mathematica/            # Mathematica notebooks
├── supplementary-informations/ # Supporting materials
├── .planning/                  # GSD planning documents and phase tracking
│   ├── codebase/               # Codebase analysis documents (this directory)
│   ├── phases/                 # Phase plans and verification documents
│   ├── milestones/             # Milestone definitions
│   ├── quick/                  # Quick-fix task plans
│   └── todos/                  # Pending and done task lists
├── .github/workflows/          # CI configuration
├── Project.toml                # Julia package manifest: dependencies and compat bounds
└── Manifest.toml               # Exact dependency versions (lockfile equivalent)
```

## Directory Purposes

**`src/`:**
- Purpose: The entire library — a single Julia module `QuantumFurnace`
- Contains: All `.jl` source files; no subdirectories
- Key files: `QuantumFurnace.jl` (entry), `structs.jl` (types), `furnace.jl` (main API), `trajectories.jl` (trajectory engine), `krylov_workspace.jl` + `krylov_matvec.jl` + `krylov_eigsolve.jl` (matrix-free path)

**`test/`:**
- Purpose: Complete test suite run via `julia --project test/runtests.jl`
- Contains: One `test_*.jl` file per feature area; shared `test_helpers.jl` with fixtures
- Key files: `runtests.jl` (orchestrator), `test_helpers.jl` (shared setup), `test_krylov_crossvalidation.jl` (dense vs Krylov comparison)

**`hamiltonians/`:**
- Purpose: Pre-computed, serialized `HamHam` objects for standard test systems (disordered Heisenberg chains)
- Contains: `heis_disordered_periodic_n{3..10}.bson`; loaded via `load_hamiltonian("heis", n; beta=...)`
- Generated: Via `hamiltonians/generate_hamiltonians.jl`; committed to repo

**`simulations/`:**
- Purpose: Runnable scripts for full simulation runs (not unit tests; intended for cluster)
- Contains: `main_liouv.jl`, `main_thermalize.jl`, `main_krylov_benchmark.jl`, `run_julia.sbatch`
- Pattern: Each script loads the module with `includet` or `using QuantumFurnace`, defines a `main()` function, and calls it

**`experiments/`:**
- Purpose: Analysis and validation scripts with stored output from past runs
- Contains: `.jl` scripts for gap validation, delta scaling, momentum diagnostics; `.bson` + `.txt` result pairs

**`results/`:**
- Purpose: Output directory for `save_experiment` calls; stores `ExperimentResult` BSON files
- Generated: By simulation scripts; not part of library logic

**`playground/`:**
- Purpose: Exploratory scripts and notebooks — NOT maintained, NOT tested
- Contains: Ad-hoc one-off investigations

## Key File Locations

**Entry Points:**
- `src/QuantumFurnace.jl`: Module definition, all `export` declarations, all `include` calls in load order
- `simulations/main_liouv.jl`: Dense Liouvillian simulation script
- `simulations/main_thermalize.jl`: DM thermalization simulation script

**Core Types:**
- `src/structs.jl`: `AbstractDomain`, `LiouvConfig`, `ThermalizeConfig`, `JumpOp`, `LindbladianResult`, `DMSimulationResult`, `ConvergenceData`, `LindbladianWorkspace`, `LSIFramework`
- `src/hamiltonian.jl`: `HamHam{T}` struct and constructors
- `src/trotter_domain.jl`: `TrottTrott{T}` struct and Trotter builder
- `src/kraus.jl`: `KrausScratch` workspace
- `src/trajectories.jl`: `TrajectoryFramework`, `TrajectoryWorkspace`, `PerOperatorKraus`, `TrajectoryResult`, `ObservableTrajectoryResult`

**Simulation Logic:**
- `src/furnace.jl`: `run_lindbladian`, `construct_lindbladian`, `run_thermalization`
- `src/jump_workers.jl`: `_jump_contribution!` (all domain variants), Liouvillian vectorization helpers
- `src/furnace_utensils.jl`: `_precompute_data` dispatch family
- `src/coherent.jl`: `B_time`, `B_trotter`, coherent precomputation helpers

**Krylov / Matrix-Free:**
- `src/krylov_workspace.jl`: `KrylovWorkspace` struct and constructor
- `src/krylov_matvec.jl`: `apply_lindbladian!`, `apply_adjoint_lindbladian!`, sandwich helpers
- `src/krylov_eigsolve.jl`: `krylov_spectral_gap`, `apply_delta_channel!`, `KrylovGapResult`

**Analysis and Post-Processing:**
- `src/gap_estimation.jl`: `estimate_spectral_gap`, `SpectralGapResult`, `eigenbasis_overlap_analysis`
- `src/convergence.jl`: `run_trajectories_convergence`, `run_trajectories_adaptive`, `build_preset_trajectory_observables`
- `src/fitting.jl`: `fit_exponential_decay`, `FitResult`
- `src/diagnostics.jl`: `run_exact_diagnostics`, six DIAG-0x functions

**Persistence:**
- `src/results.jl`: `ExperimentResult`, `save_experiment`, `load_experiment`

**Configuration:**
- `Project.toml`: Package metadata, dependency list, compat bounds
- `Manifest.toml`: Exact resolved dependency versions (lockfile)

## Naming Conventions

**Files:**
- `snake_case.jl` for all source files
- Files named by their primary concept: `hamiltonian.jl`, `trajectories.jl`, `krylov_matvec.jl`
- Test files prefixed `test_`: `test_krylov_matvec.jl`, `test_gap_estimation.jl`
- Experiment/simulation files prefixed by role: `main_*.jl` for simulation scripts, `validate_*.jl` for validation

**Functions:**
- Public API: `snake_case` — `run_lindbladian`, `build_trajectoryframework`, `estimate_spectral_gap`
- Internal helpers: prefixed with `_` — `_jump_contribution!`, `_precompute_data`, `_krylov_oft!`
- In-place mutation: suffixed with `!` — `step_along_trajectory!`, `apply_lindbladian!`, `hermitianize!`

**Types/Structs:**
- `PascalCase` — `HamHam`, `TrottTrott`, `LiouvConfig`, `KrylovWorkspace`, `TrajectoryFramework`
- Abstract types prefixed `Abstract` — `AbstractDomain`, `AbstractConfig`, `AbstractLiouvConfig`
- Domain singletons suffixed `Domain` — `BohrDomain`, `EnergyDomain`, `TimeDomain`, `TrotterDomain`
- Result structs suffixed `Result` — `LindbladianResult`, `SpectralGapResult`, `KrylovGapResult`, `FitResult`

**Variables:**
- `snake_case` throughout: `num_qubits`, `beta`, `ham_or_trott`, `precomputed_data`
- Type parameters: single uppercase letters — `T` for float type, `D` for domain, `C` for config, `CT` for complex

## Where to Add New Code

**New Domain Variant:**
- Add domain type to `src/structs.jl` (new `struct XDomain <: AbstractDomain end`)
- Add `_precompute_data` dispatch method in `src/furnace_utensils.jl`
- Add `_jump_contribution!` dispatch method in `src/jump_workers.jl`
- Add coherent term function in `src/coherent.jl` if needed
- Add to `DOMAIN_LOOKUP` dict in `src/results.jl` for serialization

**New Config Type:**
- Add `@kwdef struct` in `src/structs.jl`, inheriting from `AbstractLiouvConfig` or `AbstractThermalizeConfig`
- Add `pick_transition` dispatch in `src/energy_domain.jl`
- Export from `src/QuantumFurnace.jl`

**New Result Struct:**
- Add to `src/structs.jl` (or the file most relevant to its consumer)
- Export from `src/QuantumFurnace.jl`

**New Analysis Function:**
- Trajectory-based analysis: add to `src/gap_estimation.jl` or `src/convergence.jl`
- Dense exact analysis: add to `src/diagnostics.jl`
- Post-processing / fitting: add to `src/fitting.jl`

**New Test:**
- Add to an existing `test/test_*.jl` file if it fits an existing area
- Create `test/test_<feature>.jl` for new feature areas; add `include("test_<feature>.jl")` to `test/runtests.jl`
- Add shared fixtures/helpers to `test/test_helpers.jl`

**New Experiment Script:**
- Place in `experiments/` as `validate_*.jl` or `diagnose_*.jl`
- Store outputs as `.bson` + `.txt` pairs in `experiments/` or `results/`

**New Simulation Script:**
- Place in `simulations/` as `main_*.jl`

**Utilities:**
- Quantum information math (distances, norms, Paulis): `src/qi_tools.jl`
- Low-level matrix utilities (Kronecker, in-place linear algebra): `src/misc_tools.jl`
- Physical constants: `src/constants.jl`

## Special Directories

**`hamiltonians/`:**
- Purpose: Pre-computed serialized Hamiltonians for use in tests and experiments
- Generated: By `hamiltonians/generate_hamiltonians.jl`
- Committed: Yes — enables tests to run without regenerating expensive eigendecompositions

**`test/reference/`:**
- Purpose: Reference numerical outputs for regression tests
- Generated: By first run of regression tests (golden values)
- Committed: Yes

**`test/old_tests/`:**
- Purpose: Archived tests no longer included in `runtests.jl`
- Generated: No
- Committed: Yes (for historical reference)

**`docs/build/`:**
- Purpose: Built Documenter.jl output (HTML)
- Generated: Yes, by `julia --project docs/make.jl`
- Committed: No (generated artifact)

**`.planning/`:**
- Purpose: GSD workflow planning documents — phase plans, verification, codebase analysis
- Generated: By GSD Claude commands
- Committed: Yes

---

*Structure analysis: 2026-02-25*
