# Technology Stack

**Analysis Date:** 2026-02-25

## Languages

**Primary:**
- Julia 1.12.5 - All core library code in `src/`, simulation scripts in `simulations/`, `experiments/`

**Secondary:**
- Python 3.x (3.11/3.12 seen in `__pycache__`) - Validation and prototype tooling in `tools/python/`
- Mathematica - Analytical reference notebooks in `tools/mathematica/`

## Runtime

**Environment:**
- Julia 1.12.5 (pinned in `Manifest.toml`, minimum 1.11 per `Project.toml` compat)

**Package Manager:**
- Julia built-in `Pkg`
- Lockfile: `Manifest.toml` present and committed (machine-generated)

## Frameworks

**Core:**
- No web framework - pure scientific computing library

**Testing:**
- `Test` (stdlib) v1.11+ - Test runner and `@testset`/`@test` macros
- `Aqua` v0.8 - Package quality assurance (ambiguity, stale deps, unbound type params)
- `BenchmarkTools` v1 - Performance regression benchmarking in tests
- `HypothesisTests` v0.11 - Statistical hypothesis testing for trajectory convergence
- `StableRNGs` v1 - Reproducible RNG seeding in tests
- `StatsBase` v0.34 - Statistical utilities

**Documentation:**
- `Documenter` v1 - API doc generation, deployed via GitHub Actions
- `Literate` v2 - Literate programming: `.jl` files in `docs/src/literate/` compiled to Markdown

**Build/Dev:**
- `Revise` v3 - Live code reloading during development (dev-only, in `[extras]`)
- `Debugger` v0.7 - Interactive debugging (dev-only, in `[extras]`)

## Key Dependencies

**Critical - Numerical Linear Algebra:**
- `KrylovKit` v0.8–0.10 - Matrix-free Krylov eigensolver (Arnoldi); used in `src/krylov_eigsolve.jl` for spectral gap computation without constructing full Liouvillian
- `Arpack` v0.5.4 - Dense eigensolver for full Liouvillian matrices; used in `src/gap_estimation.jl`
- `LinearAlgebra` (stdlib) - BLAS/LAPACK wrappers; `mul!`, `eigen`, `Hermitian`, `tr`, `I`
- `LinearMaps` v3 - Lazy linear map abstractions; used in `src/linearmaps_liouv.jl`
- `SparseArrays` (stdlib) - Sparse matrix support

**Critical - Physics/Signal Processing:**
- `FINUFFT` v3 - Non-Uniform Fast Fourier Transform; used in `src/nufft.jl` for operator Fourier transforms (OFTs) in `TimeDomain` and `TrotterDomain` — requires Float64 inputs
- `SpecialFunctions` v2 - `erfc` (complementary error function) for transition functions

**Optimization & Fitting:**
- `Optim` v1 - Nonlinear optimization; used in `src/log_sobolev.jl` for LSI alpha2 computation
- `LsqFit` v0.15 - Levenberg-Marquardt nonlinear least squares; used in `src/fitting.jl` for exponential decay fitting of spectral gaps
- `QuadGK` v2 - Adaptive Gauss-Kronrod integration
- `Roots` v2 - Root-finding algorithms

**Data & Persistence:**
- `BSON` v0.3 - Binary JSON serialization; used in `src/results.jl` for experiment persistence (`.bson` files in `results/` and `experiments/`)
- `DataStructures` v0.18–0.19 - Extended data structures (queues, etc.)

**Infrastructure:**
- `Distributed` (stdlib) - Multi-process parallelism; used in `src/QuantumFurnace.jl` and simulation scripts (`@everywhere` pattern)
- `SharedArrays` (stdlib) - Shared memory arrays for multi-process workflows
- `ProgressMeter` v1 - Progress bars for long-running trajectory simulations
- `Random` (stdlib) - Seeded RNG for reproducibility
- `Printf` (stdlib) - Formatted output
- `Dates` (stdlib) - Timestamp capture in experiment metadata
- `LibGit2` (stdlib) - Git hash capture for experiment provenance in `src/results.jl`
- `Pkg` (stdlib) - Project path resolution for results directory defaults
- `ClusterManagers` v0.4 - HPC cluster job management (in `[extras]`, for Slurm/SLURM HPC deployment)

## Python Tooling (non-library)

Located in `tools/python/`, these are standalone research scripts and notebooks, not part of the Julia package:

**Key Python libraries used:**
- `numpy` - Numerical arrays
- `scipy` (`linalg.logm`, `expm`) - Matrix functions
- `qutip` - Quantum optics toolbox (reference implementations)
- `qiskit` + `qiskit_aer` - Quantum circuit simulation (prototype quantum circuit implementations)

## Configuration

**Environment:**
- No `.env` files; no secrets or external service credentials required
- All parameters passed directly as Julia structs (`LiouvConfig`, `ThermalizeConfig`)
- BLAS thread count controlled explicitly in benchmark scripts: `BLAS.set_num_threads(4)`

**Build:**
- `Project.toml` - Package manifest with deps and compat bounds
- `Manifest.toml` - Exact resolved dependency tree (committed, pinned to Julia 1.12.5)
- `docs/Project.toml` - Separate docs environment (only `Documenter`, `Literate`, `QuantumFurnace`)

## Platform Requirements

**Development:**
- Julia 1.11+ (1.12.5 in lockfile)
- VSCode with Julia extension (`.vscode/settings.json`, `.vscode/launch.json` present)
- Optional: Mathematica for `tools/mathematica/` notebooks
- Optional: Python 3.11+ with `numpy`, `scipy`, `qutip`, `qiskit` for `tools/python/`

**Production/HPC:**
- Slurm HPC cluster supported via `ClusterManagers` (listed in `[extras]`)
- `simulations/` directory contains standalone scripts with `@everywhere` distributed patterns
- `experiments/` contains archived run data (`.bson` + `.txt` companion files)

---

*Stack analysis: 2026-02-25*
