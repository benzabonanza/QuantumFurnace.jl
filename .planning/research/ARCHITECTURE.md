# Architecture Research

**Domain:** Quantum Lindbladian simulation Julia package (QuantumFurnace.jl)
**Researched:** 2026-02-13
**Confidence:** HIGH (direct codebase analysis + ecosystem comparison)

## Standard Architecture

### System Overview

```
+-----------------------------------------------------------------------+
|                      USER INTERFACE LAYER                             |
|  run_lindbladian()   run_thermalization()   run_trajectories()        |
|  construct_lindbladian()   build_trajectoryframework()                |
+-------+------+--------+---------+--------+---+--------+--------------+
        |      |        |         |        |   |        |
+-------v------v--------v---------v--------+   |        |
|           DOMAIN DISPATCH LAYER              |        |
|  jump_contribution!(::BohrDomain, ...)       |        |
|  jump_contribution!(::EnergyDomain, ...)     |        |
|  jump_contribution!(::TimeDomain, ...)       |        |
|  jump_contribution!(::TrotterDomain, ...)    |        |
|  precompute_data(::Domain, ...)              |        |
|  precompute_R(::Domain, ...)                 |        |
+-------+-------+--------+--------+-----------+        |
        |       |        |        |                     |
+-------v-------v--------v--------v-----------+ +------v--------------+
|         PRECOMPUTATION LAYER                | |  COHERENT TERMS     |
|  furnace_utensils.jl:                       | |  coherent.jl:       |
|    precompute_labels()                      | |    coherent_bohr()  |
|    precompute_data()                        | |    B_time()         |
|    pick_transition()                        | |    B_trotter()      |
|    pick_alpha()                             | |    precompute_      |
|  nufft.jl:                                  | |    coherent_total_B |
|    prepare_oft_nufft_prefactors()           | |    precompute_      |
|  ofts.jl:                                   | |    coherent_unitary |
|    oft!(), time_oft!(), trotter_oft!()      | |    _terms()         |
+---------+--------+--------+-----------------+ +---------+-----------+
          |        |        |                             |
+---------v--------v--------v-----------------------------v-----------+
|                     DATA MODEL LAYER                                |
|  structs.jl:                                                        |
|    AbstractDomain -> BohrDomain, EnergyDomain, TimeDomain, Trotter  |
|    AbstractConfig -> LiouvConfig, ThermalizeConfig, *GNS variants   |
|    JumpOp, KrausScratch, LindbladWorkspace, TrajectoryWorkspace     |
|  hamiltonian.jl:                                                    |
|    HamHam (spectral decomposition + Bohr frequencies)               |
|  trotter_domain.jl:                                                 |
|    TrottTrott (Trotterized time evolution data)                     |
|  Result types: HotSpectralResults, HotAlgorithmResults              |
+-----+------+------+------+-----------------------------------------+
      |      |      |      |
+-----v------v------v------v-----------------------------------------+
|                   UTILITY LAYER                                     |
|  qi_tools.jl: kron!(), vectorize_liouv_diss_and_add!(),            |
|               trace_distance_h(), gibbs_state()                     |
|  misc_tools.jl: pad_term(), load_hamiltonian(), riemann_sum()       |
|  constants.jl: X, Y, Z, Had                                        |
|  errors.jl: validation logic                                        |
|  log_sobolev.jl: LSI alpha2 computation                            |
+--------------------------------------------------------------------+
```

### Component Responsibilities

| Component | Responsibility | Communicates With |
|-----------|----------------|-------------------|
| **furnace.jl** (Entry Points) | Top-level API: `run_lindbladian`, `run_thermalization`, `construct_lindbladian`. Orchestrates validation, precomputation, loop accumulation, spectral analysis. | All components via orchestration |
| **jump_workers.jl** (Domain Dispatch) | Per-domain `jump_contribution!` methods for both Liouvillian vectorization and density-matrix Kraus stepping. The central multiple-dispatch hub. | domains, precompute_data, qi_tools, coherent |
| **trajectories.jl** (Trajectory Framework) | `TrajectoryFramework` construction, `step_along_trajectory!`, `evolve_along_trajectory`, `run_trajectories`. State-vector trajectory simulation. | jump_workers (precompute_R), coherent, furnace_utensils |
| **furnace_utensils.jl** (Precomputation) | `precompute_data()` dispatcher, energy/time label creation, transition function selection (`pick_transition`, `pick_alpha`). | energy_domain, nufft, coherent, ofts |
| **structs.jl** (Type Definitions) | All type definitions: domains, configs, workspaces, results. The shared vocabulary across all layers. | Imported by every other file |
| **hamiltonian.jl** (Hamiltonian) | `HamHam` construction, spectral decomposition, Bohr frequency computation, disordering, `finalize_hamham()`. | structs, qi_tools, constants |
| **trotter_domain.jl** (Trotter) | `TrottTrott` construction, Trotterization of time evolution, Trotter error computation. | hamiltonian |
| **coherent.jl** (Coherent Terms) | B operator computation for exact KMS detailed balance. Per-domain methods: `coherent_bohr`, `B_time`, `B_trotter`. Unitary exponentiation for Kraus. | hamiltonian, bohr_domain, furnace_utensils |
| **bohr_domain.jl** (Bohr Specifics) | `coherent_bohr()`, `pick_f()`, `create_f()`. Bohr-frequency-specific rate functions. | hamiltonian, structs |
| **energy_domain.jl** (Energy Specifics) | `pick_transition()` variants (Gaussian, Metropolis, Glauber, GNS). Transition weight function factory. | structs |
| **ofts.jl** (OFT Engine) | `oft!()`, `time_oft!()`, `trotter_oft!()`. In-place oscillatory Fourier transforms on jump operators. | hamiltonian, trotter_domain |
| **nufft.jl** (NUFFT Acceleration) | `prepare_oft_nufft_prefactors()`. Precompute all OFT prefactors via FINUFFT for the full energy grid. | FINUFFT library |
| **qi_tools.jl** (QI Utilities) | `kron!()`, `vectorize_liouv_diss_and_add!()`, `trace_distance_h()`, `gibbs_state()`, `fidelity()`. Zero-allocation Liouvillian assembly. | LinearAlgebra |
| **log_sobolev.jl** (LSI Bounds) | `compute_LSI_alpha2()`. Optimization-based lower bound on log-Sobolev constant from Liouvillian spectrum. | furnace results, Optim |
| **linearmaps_liouv.jl** (Sparse Liouvillian) | `construct_lindbladian_map()`. LinearMaps.jl interface for matrix-free Liouvillian-vector products. Currently incomplete. | LinearMaps, jump_workers |

## Recommended Project Structure

```
src/
+-- QuantumFurnace.jl          # Module definition, imports, include order, exports
+-- constants.jl               # Pauli matrices X, Y, Z, Had
+-- structs.jl                 # ALL type definitions (domains, configs, workspaces, results)
+-- hamiltonian.jl             # HamHam construction, spectral decomposition, Bohr freqs
+-- trotter_domain.jl          # TrottTrott, trotterize2(), Trotter error
+-- qi_tools.jl                # Quantum info: trace distance, Gibbs state, kron!, vectorize
+-- misc_tools.jl              # Pauli padding, file I/O, Riemann sums, filename generation
+-- errors.jl                  # Validation logic, custom error types
+-- bohr_domain.jl             # Bohr-specific: coherent_bohr(), pick_f(), create_f()
+-- energy_domain.jl           # Energy-specific: pick_transition() variants
+-- time_domain.jl             # Time-domain label truncation
+-- ofts.jl                    # Oscillatory Fourier transforms (oft!, time_oft!, trotter_oft!)
+-- nufft.jl                   # NUFFT prefactors via FINUFFT
+-- kraus.jl                   # KrausScratch workspace type
+-- coherent.jl                # B-term computation, coherent unitaries
+-- jump_workers.jl            # jump_contribution!() -- all domain methods
+-- furnace_utensils.jl        # precompute_data(), label creation, transition selection
+-- furnace.jl                 # Entry points: run_lindbladian, run_thermalization
+-- trajectories.jl            # TrajectoryFramework, stepping, trajectory evolution
+-- log_sobolev.jl             # LSI alpha2 optimization
+-- log_sobolev_manopt.jl      # Manifold optimization (planned)
+-- linearmaps_liouv.jl        # LinearMaps interface (incomplete)
+-- kossakowski.jl             # Kossakowski matrix tools

test/
+-- runtests.jl                # MISSING: central test runner with @testset
+-- ham_test.jl                # Interactive; needs @test conversion
+-- B_test.jl                  # Interactive; needs @test conversion
+-- trajectory_test.jl         # Interactive; needs @test conversion
+-- trott_test.jl              # Interactive; needs @test conversion
+-- time_tests.jl              # Interactive; needs @test conversion
+-- kossakowski_test.jl        # Interactive; needs @test conversion
+-- log_sobolev_test.jl        # Interactive; needs @test conversion

docs/
+-- make.jl                    # Documenter.jl + Literate.jl build script
+-- src/
    +-- index.md               # Landing page
    +-- api.md                 # Auto-generated API reference
    +-- literate/              # Source .jl files for tutorials/theory
        +-- tutorial_*.jl      # Executable tutorial scripts
        +-- theory_*.jl        # Executable theory explanation scripts
```

### Structure Rationale

- **Flat src/:** Julia convention. No subdirectories needed -- 23 files is manageable. QuantumOptics.jl (19 files) and QuantumToolbox.jl (15 files + 3 subdirs) follow similar patterns. Subdirectories only warranted when file count exceeds ~30 and clear groupings emerge (HIGH confidence: verified against both QuantumOptics.jl and QuantumToolbox.jl GitHub repos).

- **Include order in QuantumFurnace.jl respects dependency DAG:** constants -> hamiltonian -> trotter_domain -> structs -> qi_tools -> misc_tools -> ... -> furnace. This is critical in Julia: each `include()` must see types it depends on.

- **test/ currently lacks runtests.jl:** The standard Julia package test runner expects `test/runtests.jl` with `@testset` blocks. Current tests are interactive scripts using `Revise.includet()` -- they cannot be run by `Pkg.test()`. This is the most urgent structural gap.

- **docs/ uses Documenter.jl + Literate.jl:** This is the gold-standard Julia documentation approach. Literate.jl source files in `docs/src/literate/` are executable Julia scripts whose comments become markdown. Documenter.jl builds HTML from these. CI deploys via GitHub Actions (HIGH confidence: verified Documentation.yml exists).

## Architectural Patterns

### Pattern 1: Domain Dispatch via Singleton Types

**What:** Four empty structs (`BohrDomain`, `EnergyDomain`, `TimeDomain`, `TrotterDomain`) subtyping `AbstractDomain`. Functions dispatch on domain type as a positional argument rather than using `if-else` chains.

**When to use:** Any time behavior varies by approximation level. This is the core extensibility mechanism.

**Trade-offs:**
- Pro: Adding a new domain means adding new method definitions, not modifying existing code (open/closed principle)
- Pro: Julia's compiler specializes per-domain, eliminating runtime dispatch overhead
- Con: All domain methods must have matching signatures; signature drift between domains becomes a bug source
- Con: When two domains share logic (TimeDomain and TrotterDomain often do), code duplication is tempting -- use `::Union{TimeDomain, TrotterDomain}` dispatch carefully

**Example:**
```julia
# Each domain gets its own method -- no if/else needed
function jump_contribution!(L, ::BohrDomain, jump, ham, config, precomp, ws; kw...)
    # Bohr-specific: iterate over Bohr frequency buckets
end

function jump_contribution!(L, ::EnergyDomain, jump, ham, config, precomp, ws; kw...)
    # Energy-specific: iterate over energy labels with transition weights
end

function jump_contribution!(L, ::Union{TimeDomain, TrotterDomain}, jump, ham, config, precomp, ws; kw...)
    # Shared Time/Trotter: iterate with NUFFT prefactors
end
```

**Ecosystem precedent:** QuantumToolbox.jl uses the same pattern with quantum object types (Ket, Bra, Operator, SuperOperator) and solver dispatch. QuantumOptics.jl uses basis types for dispatch. This is idiomatic Julia. (MEDIUM confidence: verified from architecture papers and GitHub source listings.)

### Pattern 2: Workspace/Precomputation Cache

**What:** Expensive calculations (NUFFT prefactors, transition function evaluations, coherent term matrices) are computed once and stored in named tuples or structs. Hot loops consume these cached values without allocation.

**When to use:** Any computation that is invariant across the inner loop (e.g., invariant across jumps or time steps).

**Trade-offs:**
- Pro: Eliminates allocation in tight loops -- critical for dim^2 x dim^2 Liouvillian assembly
- Pro: Named tuples make the cached data self-documenting (`precomputed_data.gamma_norm_factor`)
- Con: Memory footprint grows (NUFFT prefactors are dim x dim x num_energies arrays)
- Con: Named tuples are not type-stable when their contents vary by domain -- QuantumFurnace.jl works around this by having `precompute_data()` return domain-specific tuples

**Example:**
```julia
# Precompute once
precomputed_data = precompute_data(config.domain, config, hamiltonian)
ws = LindbladianWorkspace(dim)

# Use many times, zero allocation per iteration
for jump in jumps
    jump_contribution!(L, config.domain, jump, hamiltonian, config, precomputed_data, ws)
end
```

**Ecosystem precedent:** QuantumOptics.jl uses pre-allocated output arrays throughout. QuantumToolbox.jl caches quantum object metadata in the Qobj struct. This is standard Julia numerical computing practice. (HIGH confidence.)

### Pattern 3: Dual Simulation Modes (Density Matrix + Trajectory)

**What:** The same physical model (Lindbladian) is simulated two ways: (1) constructing the full dim^2 x dim^2 Liouvillian superoperator for spectral analysis (`run_lindbladian`), and (2) stochastic trajectory simulation of the corresponding quantum channel (`run_thermalization`, `run_trajectories`).

**When to use:** When users need both exact spectral properties (steady state, spectral gap) and scalable stochastic sampling (convergence curves, statistics over many runs).

**Trade-offs:**
- Pro: Cross-validation: trajectory steady state should match Liouvillian fixed point
- Pro: Liouvillian gives spectral gap directly; trajectories give sampling statistics
- Con: The two modes require different `jump_contribution!` signatures (Liouvillian accumulates into dim^2 x dim^2 matrix; Kraus evolves dim x dim density matrix), leading to parallel method trees in `jump_workers.jl`
- Con: Ensuring mathematical equivalence between the two modes is a significant validation burden

**Ecosystem precedent:** QuantumOptics.jl provides `master()` (density matrix master equation) alongside `mcwf()` (Monte Carlo wave function / quantum trajectories). QuantumToolbox.jl similarly offers `mesolve()` and `mcsolve()`. This dual-mode pattern is universal in open quantum system packages. (HIGH confidence: verified from QuantumOptics.jl GitHub source listings showing both `master.jl` and `mcwf.jl`.)

## Data Flow

### Liouvillian Construction Flow

```
[Config + Hamiltonian]
    |
    v
validate_config!(config)
    |
    v
finalize_hamham(ham, beta) -----> HamHam with bohr_freqs, bohr_dict, gibbs
    |
    v
precompute_data(domain, config, ham_or_trott) -----> NamedTuple(alpha, gamma_norm_factor, ...)
    |                                                     |
    v                                                     v
precompute_coherent_total_B(jumps, ...) -----> B_total (or nothing)
    |                                              |
    v                                              v
+-- for each jump in jumps: --+     vectorize_liouvillian_coherent!(L, B, ws)
|   jump_contribution!(       |                    |
|     L_target,               |                    v
|     domain,                 |              L_target (accumulated)
|     jump,                   |
|     ham_or_trott,           |
|     config,                 |
|     precomputed_data,       |
|     workspace               |
|   )                         |
+-----------------------------+
    |
    v
L_target: Full Liouvillian (dim^2 x dim^2)
    |
    v
eigs(L, nev=2, sigma=shift) -----> steady_state, spectral_gap
    |
    v
HotSpectralResults{D}
```

### Thermalization (Density Matrix) Flow

```
[Config + Hamiltonian + initial_dm]
    |
    v
precompute_data(domain, config, ham_or_trott)
    |
    v
precompute_coherent_unitary_terms(jumps, ...) -----> [U_B_1, U_B_2, ...] per jump
    |
    v
+-- for step = 1 to ceil(mixing_time/delta): --+
|   idx = rand(1:num_jumps)                     |
|   jump = jumps[idx]                           |
|   jump_contribution!(domain,                  |
|     evolving_dm,                              |
|     jump,                                     |
|     ham_or_trott,                             |
|     config,                                   |
|     precomputed_data,                         |
|     scratch;                                  |
|     coherent_unitary_cache=U_B[idx],          |
|     jump_prob=1/num_jumps)                    |
|                                               |
|   Inside jump_contribution!:                  |
|     1. Apply U_B: dm <- U_B dm U_B'          |
|     2. Compute R = sum L_k' L_k              |
|     3. Compute rho_jump = delta * sum_w ...   |
|     4. K0 = I - alpha*R                       |
|     5. S = (2alpha - delta)R - alpha^2 R^2    |
|     6. U_res = cholesky(S).U                  |
|     7. dm_next = K0 dm K0' + rho_jump         |
|                  + U_res dm U_res'             |
+-----------------------------------------------+
    |
    v
HotAlgorithmResults{D} (distances_to_gibbs, time_steps)
```

### Trajectory (State Vector) Flow

```
[Config + Hamiltonian + psi0]
    |
    v
build_trajectoryframework(jumps, ham, config, delta)
    |
    v  Precomputes once:
    |    R = sum_jump,w L_w^dag L_w
    |    K0 = I - alpha*R
    |    S = (2alpha-delta)R - alpha^2 R^2
    |    U_res = cholesky(S).U
    |    U_B = exp(-i delta B_total)  [if coherent]
    |
    v
TrajectoryFramework (immutable cache)
    |
    v
+-- for step = 1 to num_steps: ----------+
|   step_along_trajectory!(psi, fw, rng)  |
|                                         |
|   1. Apply U_B: psi <- U_B * psi       |
|   2. Apply K0: psi_next = K0 * psi     |
|   3. Sample jump k ~ probability       |
|   4. Branch:                            |
|      - no-jump: psi = K0*psi (norm)    |
|      - jump k:  psi = L_k*psi (norm)   |
|   5. Apply U_res: psi <- U_res * psi   |
+-----------------------------------------+
    |
    v
Final |psi> or accumulated rho = E[|psi><psi|]
```

### Key Data Flows

1. **Config -> Precomputation -> Accumulation:** Configuration drives precomputation (which labels, which transition function, which domain), and precomputed data feeds the tight inner loops. This is a compile-once-run-many pattern.

2. **Hamiltonian -> JumpOp Basis Transform:** Jump operators are created in the computational basis, then rotated into the Hamiltonian eigenbasis (or Trotter eigenbasis). This basis-transformed `in_eigenbasis` field is what all domain methods use.

3. **Domain -> Method Selection:** The domain type flows through as a dispatch tag. It never carries data -- it only selects which code path runs. All domain-varying data is in `precomputed_data`.

## Scaling Considerations

| Scale | Architecture Adjustments |
|-------|--------------------------|
| 3-6 qubits (dim 8-64) | Current architecture works well. Dense Liouvillian fits in RAM. All domains fast. |
| 7-10 qubits (dim 128-1024) | Liouvillian is dim^2 x dim^2 = 16K-1M entries. Still fits in memory. NUFFT prefactors (dim x dim x num_energies) become the dominant cost. |
| 11-12 qubits (dim 2048-4096) | Liouvillian is 4M-16M entries = 128MB-1GB dense. Full spectral analysis via eigs becomes expensive. Trajectory mode becomes preferred. |
| 13+ qubits (dim 8192+) | Dense Liouvillian infeasible (>64GB). Must use LinearMaps (matrix-free) + iterative eigensolvers, or trajectory-only mode. The incomplete `linearmaps_liouv.jl` would need to be completed. |

### Scaling Priorities

1. **First bottleneck: NUFFT prefactor memory.** The `NUFFTPrefactors` array is `dim x dim x num_energies`. For 12 qubits with `num_energy_bits=12`, this is `4096 x 4096 x 4096 = 64B ComplexF64 entries = ~512GB`. The current truncation (`truncate_energy_labels`) mitigates this, but it determines the practical qubit limit.

2. **Second bottleneck: Dense Liouvillian assembly.** The dim^2 x dim^2 matrix grows as 2^(4n). Beyond 10 qubits, switching to LinearMaps (matrix-free) or sparse representation is necessary for `run_lindbladian`.

3. **Third bottleneck: Trajectory parallelism.** For statistical convergence, thousands of trajectory samples are needed. The `@distributed` infrastructure is commented out but present. Multi-threaded trajectory sampling with shared `TrajectoryFramework` is the natural scaling path.

## Anti-Patterns

### Anti-Pattern 1: Interactive Tests Masquerading as Automated Tests

**What people do:** Write test scripts that use `includet()` and `display(norm(...))` instead of `@test` assertions and `@testset` blocks.

**Why it's wrong:** `Pkg.test()` cannot run these. CI cannot catch regressions. No `runtests.jl` entry point exists. This is the current state of QuantumFurnace.jl's test suite.

**Do this instead:** Create `test/runtests.jl` with `@testset` blocks that `include()` individual test files. Convert manual `norm(expected - actual)` checks to `@test norm(expected - actual) < tolerance`. Keep the interactive scripts in `playground/` for development; `test/` should be automated. Standard Julia practice (HIGH confidence: see Julia documentation and all major packages).

### Anti-Pattern 2: Duplicated Config Fields Across Struct Variants

**What people do:** Define 4+ config structs (`LiouvConfig`, `LiouvConfigGNS`, `ThermalizeConfig`, `ThermalizeConfigGNS`) with nearly identical field lists, differing only in a few defaults or constraints.

**Why it's wrong:** Adding a new field requires editing 4 structs. Forgetting one leads to silent bugs. The GNS variants differ only in `with_coherent = false` default and `pick_transition` dispatch.

**Do this instead:** Use a single parametric config type with a `detailed_balance_type::Symbol` field (`:KMS` or `:GNS`), or use a shared inner struct with an outer wrapper. Julia's `@kwdef` with default overrides can handle the variation. Alternatively, keep the hierarchy but extract shared fields into a composition struct. (MEDIUM confidence: this is a design judgment call; the current approach works but scales poorly if more variants are added.)

### Anti-Pattern 3: Named Tuples for Complex Precomputed Data

**What people do:** Return `precompute_data()` as a `NamedTuple` whose fields vary by domain, making the return type unstable and hard to document.

**Why it's wrong:** Different domains return tuples with different field names. Code that accesses `precomputed_data.oft_nufft_prefactors` will error silently if called with BohrDomain data. No static type checking is possible.

**Do this instead:** Define per-domain `PrecomputedData` structs (e.g., `BohrPrecomputed`, `EnergyPrecomputed`, `TimePrecomputed`) inheriting from an abstract type. This enables both documentation and type-safe dispatch. (MEDIUM confidence: named tuples are idiomatic for quick prototyping in Julia, but structs are better for stable APIs.)

## Integration Points

### External Services

| Service | Integration Pattern | Notes |
|---------|---------------------|-------|
| **FINUFFT** (NUFFT library) | Julia wrapper via `FINUFFT.jl`. Plan-based API: `finufft_makeplan -> finufft_setpts! -> finufft_exec! -> finufft_destroy!` | Most computationally critical external dependency. Used in `nufft.jl` for OFT acceleration. Version 3.4.2 in compat. |
| **Arpack** (Eigenvalue solver) | Julia wrapper via `Arpack.jl`. Used for shift-invert eigenvalue computation of the Liouvillian. | `eigs(L, nev=2, sigma=shift)` in `furnace.jl`. Critical for extracting steady state and spectral gap. |
| **BSON** (Serialization) | Julia package for binary serialization. Used for Hamiltonian caching and result persistence. | Simple load/save pattern. Files in `hamiltonians/` directory. |
| **Qiskit** (planned) | Via PythonCall.jl (recommended over PyCall.jl). Would generate quantum circuits from Trotter domain gate sequences. | Not yet implemented. PythonCall.jl is actively maintained as of 2025 and supports bidirectional Julia-Python calls. (MEDIUM confidence: based on ecosystem research.) |

### Internal Boundaries

| Boundary | Communication | Notes |
|----------|---------------|-------|
| **Config -> Precomputation** | Config struct fields read directly by `precompute_data()`. No intermediate protocol. | Tightly coupled by field names. Changes to config fields propagate to all precomputation methods. |
| **Precomputation -> Domain Methods** | `precomputed_data` named tuple destructured at method entry: `(; alpha, gamma_norm_factor) = precomputed_data`. | Each domain method expects specific fields. No shared interface type. |
| **Domain Methods -> Workspace** | Workspace structs passed as mutable arguments. Methods write into workspace buffers, caller reads back. | Classic Fortran-style buffer passing. Efficient but requires disciplined buffer management. |
| **Hamiltonian -> JumpOp** | JumpOps constructed externally using `hamiltonian.eigvecs` for basis rotation, then passed into domain methods. | JumpOps are immutable after creation. The `in_eigenbasis` field is the primary data consumed. |
| **Liouvillian mode <-> Trajectory mode** | Share: `precompute_data()`, `HamHam`, `JumpOp`, `config`. Differ: `jump_contribution!` signatures (matrix vs density matrix target). | Two parallel method trees in `jump_workers.jl`. Cross-validation between modes is a key correctness check. |
| **Julia <-> Python (planned)** | Via PythonCall.jl for Qiskit circuit generation. Julia exports gate sequences, Python constructs `QuantumCircuit` objects. | Boundary should be a simple data exchange: list of gates + parameters. Avoid passing Julia objects directly into Python. |

## Build Order Implications

The dependency structure of QuantumFurnace.jl components determines the optimal build (development) order for new features. Components higher in this graph must be stable before components below them can be developed or validated.

### Dependency Graph

```
                    constants.jl
                        |
                    hamiltonian.jl
                   /            \
          trotter_domain.jl    structs.jl
                   \          / |  \
                    \        /  |   \
                     \      /   |    \
                  qi_tools.jl  misc_tools.jl  errors.jl
                      |
              +-------+-------+
              |       |       |
         bohr_    energy_   time_
         domain   domain    domain
              \       |       /
               \      |      /
            ofts.jl  nufft.jl
                 \     |
                  \    |
               coherent.jl
                    |
               kraus.jl
                    |
            furnace_utensils.jl
                    |
            jump_workers.jl
              /           \
      furnace.jl     trajectories.jl
          |                |
    log_sobolev.jl   linearmaps_liouv.jl
```

### Suggested Build Order for New Features

1. **Foundation (must be stable first):**
   - `structs.jl` -- type definitions are consumed everywhere
   - `hamiltonian.jl` -- HamHam is passed to every computation
   - `constants.jl` -- trivial but must exist

2. **Domain infrastructure (parallel development possible):**
   - `trotter_domain.jl`, `bohr_domain.jl`, `energy_domain.jl`, `time_domain.jl`
   - `ofts.jl`, `nufft.jl` -- numerical transforms
   - These can be developed in parallel since each domain is independent

3. **Core computation (depends on 1 + 2):**
   - `furnace_utensils.jl` -- precomputation wiring
   - `coherent.jl` -- B-term computation
   - `jump_workers.jl` -- the central dispatcher

4. **Entry points (depends on 1 + 2 + 3):**
   - `furnace.jl` -- orchestration
   - `trajectories.jl` -- trajectory simulation

5. **Analysis tools (depends on 4):**
   - `log_sobolev.jl` -- requires complete Liouvillian
   - `linearmaps_liouv.jl` -- requires jump_workers interface

6. **External integration (depends on 4):**
   - Qiskit interop -- requires Trotter domain gate sequences
   - Test infrastructure -- requires all of the above to be testable
   - Documentation -- requires stable API

### Feature-Specific Build Dependencies

| Feature | Depends On | Can Start After |
|---------|-----------|-----------------|
| Trajectory validation | `trajectories.jl` + `furnace.jl` (for cross-validation DM reference) | Core computation layer complete |
| Test infrastructure (`runtests.jl`) | All testable components + stable APIs | Foundation layer complete (can be built incrementally) |
| Qiskit interop | `trotter_domain.jl` gate sequences + PythonCall.jl | Trotter domain validated |
| Documentation generation | Stable public API + Documenter.jl + Literate.jl | Entry points layer stable |
| Additional Hamiltonians (Ising, 2D) | `hamiltonian.jl` constructor extensions | Foundation layer complete |
| LinearMaps sparse Liouvillian | `jump_workers.jl` interface | Core computation layer complete |
| Multi-threaded trajectories | `trajectories.jl` + SharedArrays | Trajectory framework validated |

## How Mature Julia Packages Handle Key Concerns

### Testing

Mature Julia quantum packages use `test/runtests.jl` as the central test runner with nested `@testset` blocks. QuantumOptics.jl and QuantumToolbox.jl both follow this pattern. The standard approach:

```julia
# test/runtests.jl
using QuantumFurnace
using Test

@testset "QuantumFurnace.jl" begin
    @testset "Hamiltonian" begin include("test_hamiltonian.jl") end
    @testset "Bohr Domain" begin include("test_bohr.jl") end
    @testset "Energy Domain" begin include("test_energy.jl") end
    @testset "Trajectories" begin include("test_trajectories.jl") end
end
```

For numerical packages, tolerance-based testing is standard:

```julia
@test norm(B_bohr - B_time) < 1e-8  # Cross-domain validation
@test trace_distance_h(steady_state, gibbs) < 1e-6  # Convergence
@test isapprox(tr(dm), 1.0, atol=1e-14)  # Physical constraints
```

(HIGH confidence: verified from Julia docs and multiple package structures.)

### Documentation

The Documenter.jl + Literate.jl combination used by QuantumFurnace.jl is the ecosystem standard. Key practices from mature packages:

- Docstrings on all exported functions with `@ref` cross-links to related functions
- Literate.jl tutorials that are both executable Julia scripts and rendered documentation
- CI auto-deploys docs on push to main (already configured in `.github/workflows/Documentation.yml`)
- Theory pages with LaTeX math rendered via KaTeX in Documenter.jl
- API reference page auto-generated from docstrings via `@autodocs`

(HIGH confidence: QuantumFurnace.jl already uses this stack. The setup is correct; content needs expansion.)

### Python Interop

For Qiskit integration, PythonCall.jl is the recommended bridge over the older PyCall.jl because:

- Actively maintained (improvements noted in September 2025)
- Supports a wider range of type conversions
- Bidirectional: Julia can call Python and Python can call Julia (via JuliaCall)
- More robust garbage collection across language boundaries

The recommended pattern for quantum circuit interop: export gate parameters as plain Julia arrays/dicts, then construct Qiskit `QuantumCircuit` objects on the Python side. Avoid passing complex Julia types across the boundary.

(MEDIUM confidence: PythonCall.jl is well-documented but specific Qiskit integration patterns are not widely documented in the Julia ecosystem.)

## Sources

- [QuantumOptics.jl GitHub repository](https://github.com/qojulia/QuantumOptics.jl) -- src directory structure, master.jl + mcwf.jl dual mode
- [QuantumToolbox.jl GitHub repository](https://github.com/qutip/QuantumToolbox.jl) -- src directory structure, Qobj design, multiple dispatch
- [QuantumToolbox.jl paper (Quantum, 2025)](https://quantum-journal.org/papers/q-2025-09-29-1866/) -- architecture description, AD integration, multiple dispatch patterns
- [QuantumOptics.jl paper (arXiv:1707.01060)](https://arxiv.org/pdf/1707.01060) -- basis system, operator abstraction, open quantum system framework
- [Yao.jl (yaoquantum.org)](https://yaoquantum.org/) -- meta-package structure, QBIR intermediate representation, component packages
- [PythonCall.jl documentation](https://juliapy.github.io/PythonCall.jl/stable/) -- Julia-Python interop, comparison to PyCall.jl
- [Julia Test module documentation](https://docs.julialang.org/en/v1/stdlib/Test/) -- @testset, @test, runtests.jl standard
- [Documenter.jl documentation](https://documenter.juliadocs.org/stable/) -- documentation generation for Julia packages
- [Literate.jl documentation](https://fredrikekre.github.io/Literate.jl/stable/) -- executable documentation from Julia scripts
- [Julia Package Testing Best Practices (Great Lakes Consulting)](https://blog.glcs.io/package-testing) -- test organization patterns
- [Bridging Worlds: Julia-Python Interoperability (arXiv:2404.18170)](https://arxiv.org/html/2404.18170v1) -- PythonCall vs PyCall analysis

---
*Architecture research for: Quantum Lindbladian simulation Julia package (QuantumFurnace.jl)*
*Researched: 2026-02-13*
