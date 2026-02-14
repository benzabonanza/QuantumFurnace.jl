# Architecture

**Analysis Date:** 2026-02-13

## Pattern Overview

**Overall:** Domain-layered quantum simulation framework with progressive approximation levels

**Key Characteristics:**
- Four approximation domains (Bohr → Energy → Time → Trotter) representing different levels of physical approximation
- Central configuration-driven pattern: configuration objects determine behavior across all layers
- Workspace/precomputation pattern: expensive calculations cached upfront, consumed by operators
- Generic over domain types: single jump operator contribution logic dispatches by domain

## Layers

**Physical Abstraction Layers (Approximation Domains):**

The core innovation is a hierarchy of domains, each representing a different level of physical approximation. All domains implement the same computational interface but trade accuracy for performance.

- **BohrDomain** (`bohr_domain.jl`)
  - Purpose: Decompose jump operators by exact Bohr frequencies (ω_ij = E_i - E_j)
  - Location: `src/bohr_domain.jl`
  - Consumes: Hamiltonian eigendecomposition, bohr_dict (eigenindex mapping)
  - Produces: Direct spectral jump contributions
  - Used by: `jump_contribution!` for BohrDomain dispatch

- **EnergyDomain** (`energy_domain.jl`)
  - Purpose: Approximate Bohr frequencies as continuous energy integrals via Gaussian envelope
  - Location: `src/energy_domain.jl`
  - Consumes: Transition weight functions (Gaussian, Metropolis, Glauber), energy labels
  - Produces: Energy-discretized jump contributions
  - Used by: `jump_contribution!` for EnergyDomain dispatch

- **TimeDomain** (`time_domain.jl`, `ofts.jl`)
  - Purpose: Expand energy approximations via time-domain Fourier transforms (OFT: Oscillatory Fourier Transform)
  - Location: `src/time_domain.jl`, `src/ofts.jl`
  - Consumes: Time labels, time_oft_! caches with NUFFT prefactors
  - Produces: Time-discretized approximations of energy integrals
  - Used by: `jump_contribution!` for TimeDomain dispatch

- **TrotterDomain** (`trotter_domain.jl`)
  - Purpose: Implement time evolution via Trotterization for quantum circuit compilation
  - Location: `src/trotter_domain.jl`
  - Consumes: Trotter unitaries, quasi-Bohr frequencies from Trotter eigenvalues
  - Produces: Circuit-implementable jump operators
  - Used by: `jump_contribution!` for TrotterDomain dispatch

**Computation Layer:**

- **Configuration & Validation** (`structs.jl`)
  - Purpose: Hold all simulation parameters and validate consistency across domains
  - Location: `src/structs.jl`
  - Core types: `LiouvConfig`, `LiouvConfigGNS`, `ThermalizeConfig`, `ThermalizeConfigGNS`
  - Validation: `validate_config!` ensures coherence across domain parameters

- **Jump Operators** (`structs.jl`, `jump_workers.jl`)
  - Purpose: Encapsulate dissipative dynamics (Lindblad jump operators)
  - Location: `src/jump_workers.jl`
  - Core type: `JumpOp(data, in_eigenbasis, orthogonal, hermitian)`
  - Key function: `jump_contribution!(L_target, domain, jump, ...)` — main dispatcher
  - Pattern: Single function with domain-specific methods handles all jump accumulation

- **Hamiltonian Management** (`hamiltonian.jl`)
  - Purpose: Encapsulate system Hamiltonian with spectral decomposition and Bohr frequencies
  - Location: `src/hamiltonian.jl`
  - Core type: `HamHam(data, bohr_freqs, bohr_dict, eigvals, eigvecs, gibbs, ...)`
  - Functions: `HamHam()` constructor, `find_ideal_heisenberg()`, `finalize_hamham()`
  - Caches: Bohr frequency matrix (dim² × dim²), bohr_dict (frequency → CartesianIndex map)

- **Liouvillian Construction** (`furnace.jl`)
  - Purpose: Build vectorized Lindbladian superoperators ℒ (dim² × dim²)
  - Location: `src/furnace.jl`
  - Core functions:
    - `run_lindbladian(jumps, config, hamiltonian)` → `HotSpectralResults`
    - `construct_lindbladian(jumps, config, hamiltonian)` → sparse/dense Liouvillian matrix
  - Pattern: Accumulates per-jump contributions in-place without allocating full dim² × dim² per jump
  - Workspace: `LindbladianWorkspace` holds reusable buffers for jump ops and products

- **Thermalization Simulation** (`furnace.jl`, `trajectories.jl`)
  - Purpose: Implement step-by-step quantum algorithm emulation
  - Location: `src/furnace.jl`, `src/trajectories.jl`
  - Core functions:
    - `run_thermalization(jumps, config, dm_init, hamiltonian)` → `HotAlgorithmResults`
    - `build_trajectoryframework()` → `TrajectoryFramework` (precomputed Kraus operators)
    - `step_along_trajectory!(dm, framework, rng)`
  - Pattern: Precompute all Kraus operators ({K_0, K_res, U_coherent}) once, reuse per step
  - Measurement: Track trace distance to Gibbs state at each step

**Infrastructure Layer:**

- **Precomputation** (`furnace_utensils.jl`, `coherent.jl`)
  - Purpose: Cache expensive calculations (energy labels, transition functions, coherent terms)
  - Location: `src/furnace_utensils.jl`, `src/coherent.jl`
  - Core function: `precompute_data(domain, config, ham_or_trott)` → precomputed_data tuple
  - Outputs: Domain-specific cached data (e.g., alpha function, gamma_norm_factor, OFT prefactors)

- **Coherent Terms** (`coherent.jl`, `bohr_domain.jl`)
  - Purpose: Compute B = sum_k B_k (coherent correction for exact detailed balance)
  - Location: `src/coherent.jl`, `src/bohr_domain.jl`
  - Key functions:
    - `precompute_coherent_total_B()` — computes total coherent term
    - `precompute_coherent_unitary_terms()` — exp(-iδ B_k) per jump for Kraus
    - `coherent_bohr()` — domain-specific B computation
  - Pattern: Reused in both Liouvillian and trajectory (Kraus) contexts

- **Quantum Information Tools** (`qi_tools.jl`)
  - Purpose: Efficient quantum state operations and metrics
  - Location: `src/qi_tools.jl`
  - Functions: `kron!()` (in-place Kronecker), `vectorize_liouv_diss_and_add!()`, `trace_distance_h()`, `fidelity()`, `gibbs_state()`, `gibbs_state_in_eigen()`
  - Pattern: In-place operations to minimize allocations in large-dim Lindbladian assembly

- **Numerical Transforms** (`nufft.jl`, `ofts.jl`)
  - Purpose: Oscillatory Fourier transforms and NUFFT acceleration
  - Location: `src/nufft.jl`, `src/ofts.jl`
  - Functions: `oft!()`, `time_oft!()`, `trotter_oft!()`, `prepare_oft_nufft_prefactors()`
  - Pattern: Cache Gaussian prefactors; reuse across all jumps

- **Error Handling** (`errors.jl`)
  - Purpose: Custom exception types for validation failures
  - Location: `src/errors.jl`
  - Currently minimal; room for domain-specific error types

## Data Flow

**Liouvillian Construction Flow:**

```
1. User configures: LiouvConfig(domain, num_qubits, beta, sigma, ...)

2. validate_config!(config)
   → Ensures num_energy_bits, t0, w0, coherent flags are consistent

3. finalize_hamham(hamiltonian, beta)
   → Compute bohr_freqs, bohr_dict, gibbs state

4. construct_lindbladian(jumps, config, hamiltonian)
   → precompute_data(config.domain, config, hamiltonian)
      → Domain-specific: alpha functions, gamma_norm_factor, labels

   → For each jump in jumps:
      jump_contribution!(L, domain, jump, hamiltonian, config, precomputed_data, ws)
      → Domain-dispatched method accumulates L in-place
      → Uses precomputed_data (cached alpha, gamma, transition functions)

   → Result: Full Liouvillian L (dim² × dim²)

5. run_lindbladian(jumps, config, hamiltonian)
   → construct_lindbladian(...) → L
   → Spectral analysis via eigs(L)
   → Extract: steady_state, spectral_gap
   → Return: HotSpectralResults{D}
```

**Thermalization (Trajectory) Flow:**

```
1. build_trajectoryframework(jumps, hamiltonian, config, delta)
   → precompute_data(...) [same as Liouvillian]
   → precompute_R() — dissipative skeleton
   → K0 = (1 - √(1-δ)) I - αR  [Kraus K0 operator]
   → U_residual from Cholesky of S = 2αR - α²R²
   → precompute_coherent_unitaries: {U_{B,k}} per jump
   → Return: TrajectoryFramework (immutable precomputed structure)

2. run_thermalization(jumps, config, dm_init, hamiltonian)
   → framework = build_trajectoryframework(...)
   → For each step t = 1, 2, ..., ceil(mixing_time/delta):
      → Sample jump index k uniformly at random
      → step_along_trajectory!(dm, framework, k)
         → Apply K0 (coherent term optional)
         → Apply randomized jump (Kraus K_{k,i})
         → Apply U_residual (final Kraus)
      → Measure: trace_distance_h(dm, gibbs)
   → Return: HotAlgorithmResults{D}
```

**Energy Label & Transition Weight Flow:**

```
Energy/Time discretization:
   create_energy_labels(num_energy_bits, w0)
   → w_j = j * w0, j ∈ [-2^(num_energy_bits-1), +2^(num_energy_bits-1)]

   truncate_energy_labels(labels, config)
   → Discard labels where |transition(w)| < tolerance

   Time labels (for Time/Trotter domains):
   time_labels = energy_labels * (t0 / w0)
   truncate_time_labels_for_oft(time_labels, sigma)
   → Discard times where Gaussian decay negligible

Transition functions (domain-independent selection):
   pick_transition(config) → Returns ω ↦ γ(ω)

   For KMS (detailed balance):
   - Gaussian: exp(-(ω - ω_γ)² / (2σ_γ²))
   - Metropolis (a>0, b=0): exp(-2√((β/4)(1+4a)) |ω + βσ²/2| - βω/2 - ...)
   - Glauber (a>0, b>0): Metropolis × smoothing via erfc

   For GNS (approx. detailed balance):
   - Same forms but unshifted (ω instead of ω + βσ²/2)

Normalization:
   gamma_norm_factor = 1 / max_w γ(w)
   → All γ(w) scaled by this before use in jump_contribution!
```

## Key Abstractions

**AbstractDomain & Domain Types:**

- Purpose: Encode computational approach; no fields, only use for dispatch
- Definition: `abstract type AbstractDomain end`; `struct BohrDomain <: AbstractDomain end`, etc.
- Usage: `jump_contribution!(L, ::BohrDomain, ...)` has different implementation than `jump_contribution!(L, ::TimeDomain, ...)`
- Benefit: Single function name, type-safe per-domain logic

**AbstractConfig & Config Hierarchy:**

- Purpose: Unified parameter container validated at entry
- Hierarchy:
  ```
  AbstractConfig{D <: AbstractDomain}
    ├─ AbstractLiouvConfig{D}
    │   ├─ LiouvConfig{D}
    │   └─ LiouvConfigGNS{D}
    └─ AbstractThermalizeConfig{D}
        ├─ ThermalizeConfig{D}
        └─ ThermalizeConfigGNS{D}
  ```
- Shared fields: `num_qubits, with_coherent, with_linear_combination, domain, beta, sigma, a, b, num_energy_bits, t0, w0, eta, num_trotter_steps_per_t0`
- Thermalizing variants add: `mixing_time, delta`
- GNS variants: `with_coherent::Bool = false` (locked to false by design)

**JumpOp:**

- Purpose: Encapsulate one dissipative mode
- Fields:
  - `data::AbstractMatrix{ComplexF64}` — operator in computational basis
  - `in_eigenbasis::Matrix{ComplexF64}` — transformed to Hamiltonian eigenbasis
  - `orthogonal::Bool` — true iff A = A^T (X, Z Paulis)
  - `hermitian::Bool` — true iff A = A†
- Usage: `jump_contribution!` branches on orthogonal/hermitian flags to optimize computation

**TrajectoryFramework:**

- Purpose: Immutable cache of all precomputed Kraus operators for trajectory simulation
- Fields:
  - `jumps, ham_or_trott, config, precomputed_data` — input references
  - `B, U_B` — coherent term and its exponentiated unitary (optional)
  - `R, K0, U_residual` — Kraus skeleton operators
  - `delta, alpha` — time step and α = 1 - √(1-δ)
  - `ws::TrajectoryWorkspace` — mutable buffer for stepping
- Pattern: Build once, reuse across all trajectory steps

## Entry Points

**High-level API Functions:**

**`run_lindbladian(jumps, config, hamiltonian; trotter=nothing) → HotSpectralResults{D}`**
- Location: `src/furnace.jl` (lines 1–39)
- Triggers: User initiates spectral analysis of a single Liouvillian
- Responsibilities:
  1. Validate config
  2. Construct full Liouvillian
  3. Compute 2 smallest eigenvalues (shift-invert via Arpack)
  4. Extract steady state and spectral gap
  5. Return spectral analysis results

**`run_thermalization(jumps, config, dm_init, hamiltonian; trotter=None, rng, rescale_by_inv_prob) → HotAlgorithmResults{D}`**
- Location: `src/furnace.jl` (lines 80–145)
- Triggers: User wants step-by-step algorithm emulation
- Responsibilities:
  1. Validate config
  2. Precompute Kraus operators via `build_trajectoryframework()`
  3. Iterate: sample jump uniformly, apply Kraus, measure trace distance
  4. Track convergence; stop early if distance < 1e-5
  5. Return full evolution history

**`construct_lindbladian(jumps, config, hamiltonian; trotter=None) → Matrix{ComplexF64}`**
- Location: `src/furnace.jl` (lines 41–78)
- Triggers: Internal; called by `run_lindbladian`
- Responsibilities:
  1. Dispatch `precompute_data()` by domain
  2. Allocate workspace buffers
  3. For each jump, call `jump_contribution!` accumulating into L in-place
  4. Return full Liouvillian (no allocation of individual dim² × dim² per jump)

**Module Entry: `src/QuantumFurnace.jl`**
- Location: `src/QuantumFurnace.jl` (lines 1–67)
- Responsibilities:
  1. Define module boundary
  2. Import all dependencies
  3. Include all submodules in order (respecting DAG: constants → hamiltonian → ... → furnace)
  4. Export public API (configs, domains, entry points, utilities)

## Error Handling

**Strategy:** Configuration validation upfront; minimal runtime error recovery

**Patterns:**
- `validate_config!(config)` checks:
  - Coherence of approximation parameters (num_energy_bits, t0, w0 interplay)
  - Domain feasibility (TrotterDomain requires explicit trotter object)
  - Physical bounds (beta > 0, sigma > 0, etc.)
- Domain-specific assertions: TrotterDomain errors if `trotter === nothing`
- Linear algebra: Arpack eigs may fail to converge; wrapped with `tol=1e-12` and shift-invert for robustness
- Cholesky in trajectory: `cholesky!(..., check=false)` — numerical issues logged but proceed

**Current gaps:**
- No recovery from failed Liouvillian construction
- Limited messaging on why convergence fails
- Room for domain-specific exceptions (e.g., `BohrFrequencyResolutionError`)

## Cross-Cutting Concerns

**Logging:** Printf-based via `@printf("text")` to stdout; see `print_press(config)` for config summary

**Validation:** Central `validate_config!(config)` called at entry to `run_lindbladian` and `run_thermalization`

**Authentication:** Not applicable (no external services)

**Performance:**
- In-place operations throughout (Kronecker `kron!()`, OFT `oft!()`, jump accumulation)
- Workspace pattern: preallocate buffers, reuse across loop iterations
- NUFFT caching: `prepare_oft_nufft_prefactors()` avoids Fourier transforms in tight loops
- Thread pools: `@distributed (+)` loops in `construct_lindbladian` (commented; available)
- SharedArrays support for distributed jumps across nodes (config via `use_shared_array` flag)

---

*Architecture analysis: 2026-02-13*
