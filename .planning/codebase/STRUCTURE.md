# Codebase Structure

**Analysis Date:** 2026-02-13

## Directory Layout

```
QuantumFurnace.jl/
├── src/                      # Core package source
│   ├── QuantumFurnace.jl     # Module definition & exports
│   ├── constants.jl          # Global constants (minimal)
│   ├── structs.jl            # Type definitions (configs, domains, results)
│   ├── hamiltonian.jl        # Hamiltonian construction & management
│   ├── qi_tools.jl           # Quantum info utilities (trace distance, fidelity, Gibbs state)
│   ├── misc_tools.jl         # Miscellaneous helpers (Pauli construction, padding, etc.)
│   ├── errors.jl             # Custom error types
│   ├── bohr_domain.jl        # Bohr domain specifics (coherent_bohr, pick_f)
│   ├── energy_domain.jl      # Energy domain & transition weights (Gaussian, Metro, Glauber)
│   ├── time_domain.jl        # Time domain OFT truncation
│   ├── trotter_domain.jl     # Trotter domain & Trotterization (TrottTrott struct)
│   ├── ofts.jl               # Oscillatory Fourier transforms (oft!, time_oft!, trotter_oft!)
│   ├── nufft.jl              # NUFFT prefactors for OFT acceleration
│   ├── coherent.jl           # Coherent term computation (B terms, unitary exponents)
│   ├── kraus.jl              # Kraus operator utilities (minimal; see jump_workers.jl)
│   ├── jump_workers.jl       # Jump contribution accumulation (jump_contribution! dispatcher)
│   ├── linearmaps_liouv.jl   # LinearMaps interface for sparse Liouvillians
│   ├── furnace_utensils.jl   # Precomputation & labels (precompute_data, truncation)
│   ├── furnace.jl            # Core entry points (run_lindbladian, run_thermalization)
│   ├── trajectories.jl       # Trajectory framework & stepping (TrajectoryFramework, step_along_trajectory!)
│   ├── log_sobolev.jl        # Log-Sobolev inequality computation
│   ├── log_sobolev_manopt.jl # Manifold optimization for LSI (empty; planned)
│   └── kossakowski.jl        # Kossakowski matrix tools (minimal)
│
├── test/                     # Test suite
│   ├── ham_test.jl           # Hamiltonian construction tests
│   ├── B_test.jl             # Coherent term (B) tests
│   ├── trajectory_test.jl    # Trajectory framework tests
│   ├── trott_test.jl         # Trotter domain tests
│   ├── time_tests.jl         # Time domain tests
│   ├── kossakowski_test.jl   # Kossakowski tests
│   └── log_sobolev_test.jl   # Log-Sobolev bound tests
│
├── playground/               # Experimental scripts (not part of package)
│   ├── integrals.jl          # Numerical integration experiments
│   ├── fourier.jl            # Fourier transform exploration
│   ├── trotter_match.jl      # Trotter approximation validation
│   ├── time_energy_match.jl  # Time-energy domain agreement checks
│   ├── spectral_gaps.jl      # Gap analysis scripts
│   └── [30+ exploratory files]
│
├── hamiltonians/             # Pre-saved Hamiltonian instances
│   └── [binary .jld files]
│
├── simulations/              # Result data & simulation outputs
│   └── [data files]
│
├── results/                  # Deprecated result storage
│   └── [data files]
│
├── tools/                    # External utilities
│   ├── python/               # Python helper scripts
│   └── mathematica/          # Mathematica notebooks
│
├── docs/                     # Sphinx/Documenter documentation
│   ├── src/                  # Documentation source
│   └── build/                # Built documentation
│
├── .github/workflows/        # CI/CD (GitHub Actions)
├── .planning/codebase/       # This codebase documentation
├── Project.toml              # Package manifest & dependencies
├── Manifest.toml             # Locked dependency versions
├── README.md                 # Quick start guide
└── LICENSE                   # MIT License
```

## Directory Purposes

**`src/`:** Core package implementation
- **Purpose:** Implement quantum Gibbs sampling algorithms across 4 approximation domains
- **Contains:** 23 Julia source files organized by concern (domains, operators, utilities, API)
- **Key files:** `QuantumFurnace.jl` (module), `furnace.jl` (entry points), `structs.jl` (types)

**`test/`:** Test suite for core functionality
- **Purpose:** Validate correctness of domain approximations and algorithm implementations
- **Contains:** 7 test files (one per major module)
- **Key files:** `ham_test.jl`, `trajectory_test.jl` for integration testing
- **Run:** `julia test/runtests.jl` or `using Pkg; Pkg.test("QuantumFurnace")`

**`playground/`:** Research & exploration (excluded from package)
- **Purpose:** Experiment with new algorithms, validate math, debug edge cases
- **Contains:** 30+ one-off scripts exploring Fourier transforms, Trotter errors, spectral gaps
- **Pattern:** Each script imports main module and runs standalone; not committed to public API
- **Status:** Pre-alpha; scripts may be outdated or incomplete

**`hamiltonians/`:** Precomputed Hamiltonian cache
- **Purpose:** Store computed HamHam objects for reuse across simulations
- **Format:** BSON binary files (Julia native serialization)
- **Pattern:** Load via `load_hamiltonian("name.jld")` to avoid recomputation

**`simulations/` & `results/`:** Output data
- **Purpose:** Store evolving density matrices, convergence histories, spectral data
- **Format:** BSON or JLD2 (Julia native)
- **Usage:** Save `HotAlgorithmResults` and `HotSpectralResults` from `run_thermalization` and `run_lindbladian`

**`tools/`:** Utilities outside Julia
- **Purpose:** Python helper scripts (data post-processing), Mathematica notebooks (symbolic verification)
- **Pattern:** Not imported by main package; optional accessories

**`docs/`:** User-facing documentation
- **Purpose:** Tutorials, API reference, theory background
- **Built with:** Documenter.jl + Literate.jl (converts Julia script comments to docs)
- **Access:** `https://tembence.github.io/QuantumFurnace.jl/`

## Key File Locations

**Entry Points:**

- `src/QuantumFurnace.jl` — Module boundary; includes all submodules; exports public API
- `src/furnace.jl` — `run_lindbladian()`, `run_thermalization()`, `construct_lindbladian()`

**Configuration & Types:**

- `src/structs.jl` — All type definitions:
  - Domains: `BohrDomain`, `EnergyDomain`, `TimeDomain`, `TrotterDomain`
  - Configs: `LiouvConfig`, `LiouvConfigGNS`, `ThermalizeConfig`, `ThermalizeConfigGNS`
  - Results: `HotSpectralResults`, `HotAlgorithmResults`
  - Operators: `JumpOp`, `HamHam`, `TrottTrott`
  - Workspaces: `LindbladianWorkspace`, `TrajectoryWorkspace`, `TrajectoryFramework`

**Core Logic by Domain:**

- `src/bohr_domain.jl` — Bohr-frequency domain jump contributions; `coherent_bohr()`
- `src/energy_domain.jl` — Energy integral approximation; `pick_transition()` (all variants)
- `src/time_domain.jl` — Time label truncation for OFT
- `src/trotter_domain.jl` — Trotterization; `TrottTrott()` constructor; `trotterize2()`, `compute_trotter_error()`

**Quantum Information & Utilities:**

- `src/qi_tools.jl` — `trace_distance_h()`, `fidelity()`, `gibbs_state()`, `kron!()`, `vectorize_liouv_diss_and_add!()`
- `src/hamiltonian.jl` — `HamHam()` constructor variants; `find_ideal_heisenberg()`; `finalize_hamham()`
- `src/misc_tools.jl` — `pad_term()`, `pick_transition()` overloads, Pauli constant defs (X, Y, Z, Had)

**Jump Operator Processing:**

- `src/jump_workers.jl` — `jump_contribution!()` (main dispatcher; 4 domain methods)
- `src/coherent.jl` — `precompute_coherent_total_B()`, `precompute_coherent_unitary_terms()`, per-domain B computations
- `src/kraus.jl` — `KrausScratch` (minimal wrapper for jump evolution)

**Precomputation & Numerical Methods:**

- `src/furnace_utensils.jl` — `precompute_data()` (4 domain methods); `create_energy_labels()`, `truncate_energy_labels()`
- `src/ofts.jl` — `oft!()`, `time_oft!()`, `trotter_oft!()` (in-place OFT computations)
- `src/nufft.jl` — `prepare_oft_nufft_prefactors()` (NUFFT acceleration of OFT)

**Trajectory Simulation:**

- `src/trajectories.jl` — `TrajectoryFramework` definition; `build_trajectoryframework()`; `step_along_trajectory!()`; trajectory event dispatch
- `src/trajectories.jl` — `precompute_R()` for dissipative skeleton (Riccati-like)

**Advanced Features (Incomplete):**

- `src/linearmaps_liouv.jl` — LinearMaps.jl integration for sparse Liouvillian matrix-vector products
- `src/log_sobolev.jl` — LSI α² computation (gap-dependent lower bound on mixing time)
- `src/log_sobolev_manopt.jl` — Empty; planned manifold optimization for LSI

## Naming Conventions

**Files:**
- `_domain.jl` — Domain-specific logic (e.g., `bohr_domain.jl`, `energy_domain.jl`)
- `_test.jl` — Test file for corresponding module (e.g., `ham_test.jl` tests `hamiltonian.jl`)
- `_tools.jl`, `_utensils.jl` — Utility/helper collections (e.g., `qi_tools.jl`, `furnace_utensils.jl`)
- No file naming prefix for core (e.g., `furnace.jl` not `core_furnace.jl`)

**Functions:**
- Dispatcher pattern: `jump_contribution!(L, ::BohrDomain, ...)` — method selected by domain type
- Domain-specific variants: `coherent_bohr()`, `B_time()`, `B_trotter()` — function overloads or suffixes
- Compute vs. accumulate: `compute_trotter_error()`, `oft!()` (in-place), `precompute_data()` (cached)
- Camel case: `LiouvConfig`, `TrajectoryFramework`, `JumpOp`; snake_case functions: `run_thermalization()`, `jump_contribution!()`

**Types:**
- Abstract base types capitalized + "Domain" suffix: `AbstractDomain`, `BohrDomain`
- Config variants: `LiouvConfig`, `LiouvConfigGNS` (GNS = Gelfand-Naimark-Segal approximation)
- Result containers: `HotAlgorithmResults{D}`, `HotSpectralResults{D}` (parameterized by domain)
- Operator containers: `JumpOp`, `HamHam` (two-word names for visual clarity)
- Workspace/cache: `LindbladianWorkspace`, `TrajectoryWorkspace`, `TrajectoryFramework` (clear purpose suffix)

**Variables:**
- Hamiltonian data: `H`, `h`, `hamiltonian` (HamHam instance); `eigvals`, `eigvecs`
- Jump operators: `A`, `jump` (JumpOp instance); `jumps` (vector)
- Domain markers: `domain` (AbstractDomain), config.domain
- Liouvillian: `L`, `liouv` (Liouvillian superoperator matrix)
- Density matrix: `rho`, `dm`, `sigma` (quantum states); `gibbs` (target state)
- Bohr frequencies: `nu` (ν), `bohr_freqs` (matrix ω_ij = E_i - E_j)
- Time/energy labels: `t0`, `w0`, `energy_labels`, `time_labels`
- Configuration params: `beta` (β, inverse temp), `sigma` (σ, Gaussian width)

## Where to Add New Code

**New Feature (e.g., new transition weight type):**
- **Primary code:** `src/energy_domain.jl` — add new `_pick_transition_*()` variant
- **Tests:** `test/energy_test.jl` (create if needed) — test new transition against known KMS/GNS conditions
- **Integration:** Edit `pick_transition()` overloads to dispatch to new function based on config flags
- **Example:** To add Fermi-Dirac transition, add `_pick_transition_fermi()` and new config field `fermi_temperature`

**New Domain (e.g., a hybrid approximation):**
- **Struct definition:** `src/structs.jl` — `struct HybridDomain <: AbstractDomain end` + config variant
- **Core dispatch:** `src/jump_workers.jl` — new `jump_contribution!(L, ::HybridDomain, ...)` method
- **Precomputation:** `src/furnace_utensils.jl` — new `precompute_data(::HybridDomain, ...)` method
- **Domain-specific logic:** Create `src/hybrid_domain.jl` if methods are substantial
- **Tests:** `test/hybrid_test.jl` — validate against reference domain (e.g., verify fallback to Bohr)

**New Quantum Algorithm (e.g., adiabatic evolution):**
- **Entry point:** Add to `src/furnace.jl` — new `run_adiabatic(config, ...)` function
- **Workspace:** Add specialized workspace type to `src/structs.jl` if complex caching needed
- **Implementation:** Extend precomputation in `src/furnace_utensils.jl` for adiabatic-specific labels
- **Tests:** `test/adiabatic_test.jl` — compare against thermalization for cross-validation

**New Diagnostic or Analysis Tool:**
- **Utilities:** `src/misc_tools.jl` if < 100 lines; otherwise create `src/diagnostics_*.jl`
- **Example:** `compute_mixing_time_bound(spectral_gap, dimension)` in `misc_tools.jl`
- **Integration:** Export from `src/QuantumFurnace.jl` if public; keep private if internal-only

**Bug Fix or Optimization:**
- **Locate:** Search for issue by function name (e.g., `jump_contribution!` in `src/jump_workers.jl`)
- **Test first:** Add failing test case to `test/*.jl` reproducing issue
- **Fix:** Modify function and verify test passes; check all overloads if applicable
- **Cross-reference:** Update ARCHITECTURE.md and STRUCTURE.md if logic flow changes

## Special Directories

**`.planning/codebase/`:**
- **Purpose:** GSD (Get Shit Done) command reference documents
- **Generated:** By `/gsd:map-codebase` commands; consumed by `/gsd:plan-phase` and `/gsd:execute-phase`
- **Committed:** Yes; version-controlled documentation
- **Manual edits:** Rare; regenerated when codebase architecture changes substantially

**`playground/`:**
- **Purpose:** Experimental code; not part of package public API
- **Generated:** Developers write one-off scripts for research/debugging
- **Committed:** Yes; tracked for reproducibility, but not imported elsewhere
- **Cleanup:** Old playground files can be archived or deleted without breaking main package

**`.github/workflows/`:**
- **Purpose:** CI/CD pipeline definitions
- **Format:** YAML GitHub Actions workflows
- **Current state:** Minimal; expand as testing coverage grows

---

*Structure analysis: 2026-02-13*
