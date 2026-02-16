# Phase 18: KMS-vs-GNS Experiments - Research

**Researched:** 2026-02-16
**Domain:** Parameter sweep comparing KMS (exact) vs GNS (approximate) detailed balance convergence in QuantumFurnace.jl
**Confidence:** HIGH

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

#### Hamiltonian model
- 1D Heisenberg chain (XXX) with periodic boundaries, J=1.0 uniform coupling
- Experiment function accepts a pre-built Hamiltonian object (not constructed internally) -- supports future Hamiltonian models
- System sizes n=4,6,8 are qubit counts (Hilbert space dims 16, 64, 256)
- Trajectory-based simulation uses local jump operators (TrotterDomain), no full Lindbladian construction -- n=8 is feasible

#### Simulation parameters
- delta = 0.01 for all 27 experiments (fixed, same as Phase 14 baseline)
- Adaptive sampling with shared n_max = 10000 cap across all experiments
- Rationale for n_max: trajectory sampling error should be subdominant to the delta=0.01 systematic error (~0.03 trace distance)
- Mixing time, batch size, convergence thresholds: Claude's discretion (see below)

#### Data organization
- Output directory: `experiments/` in project root
- Descriptive file names: e.g., `kms_n4_beta5.bson`, `gns_sigma1_n6_beta10.bson`
- No summary/index file -- individual BSON files contain all metadata, glob the directory
- `experiments/` is gitignored -- data is regenerable

#### Execution strategy
- Per-experiment function: `run_experiment()` called in a loop over the parameter grid
- Sweep script lives in `experiments/run_sweep.jl` -- standalone script using the package, not part of package API
- Failure handling: skip failed experiments, log the error, save what completed, continue with remaining -- print warning at end listing failures
- No resume capability -- always rerun all experiments (overwrite existing files)

### Claude's Discretion
- Mixing time value (scaling with beta or fixed)
- Batch size for adaptive sampling
- Convergence thresholds (relative change threshold, window size)
- Exact file naming format details
- How `run_experiment()` is structured internally (kwargs, return type)
- Progress/logging output during sweep

### Deferred Ideas (OUT OF SCOPE)
None -- discussion stayed within phase scope
</user_constraints>

## Summary

Phase 18 runs a 27-experiment parameter sweep comparing KMS (exact detailed balance) and GNS (approximate detailed balance) convergence behavior using TrotterDomain trajectory simulation. The experiment grid is: n=4,6,8 qubits x beta=5,10,20 x {KMS, GNS@sigma=1/beta, GNS@sigma=0.5/beta}. Each experiment uses `run_trajectories_adaptive` (Phase 17) with a shared n_max=10000 cap and fixed delta=0.01.

The codebase is fully prepared for this phase. All required infrastructure exists:
- **Hamiltonian construction**: `HamHam(terms, coeffs, num_qubits, beta; periodic=true)` constructs the XXX Heisenberg chain directly from `terms = [[X,X], [Y,Y], [Z,Z]], coeffs = [1.0, 1.0, 1.0]`. No disorder term needed for uniform coupling.
- **TrotterDomain**: `TrottTrott(hamiltonian, t0, num_trotter_steps)` creates the Trotter object. Pre-built Hamiltonians for n=4,6,8 exist in `hamiltonians/`, but since the phase requires uniform XXX coupling (not disordered), Hamiltonians must be constructed at runtime.
- **KMS trajectories**: `ThermalizeConfig` with `with_coherent=true` and `run_trajectories_adaptive`.
- **GNS trajectories**: `ThermalizeConfigGNS` with `with_coherent=false` (enforced) and `run_trajectories_adaptive`.
- **Result persistence**: `ExperimentResult` struct with `save_experiment(result, path)` writes BSON + companion `.txt`. The `_capture_metadata` function auto-records Julia version, git hash, timestamp, thread count, and wall time.
- **Convergence tracking**: `build_convergence_observables_trotter` constructs ZZ correlations and H observable in the Trotter eigenbasis. `ConvergenceData` stores batch-level trace distances, observable values, and adaptive diagnostics.

The sweep script is a standalone Julia script (`experiments/run_sweep.jl`) that imports QuantumFurnace, constructs the parameter grid, builds Hamiltonians and jump operators, and calls `run_experiment()` for each configuration. The `run_experiment()` function is a thin wrapper that constructs the config, calls `run_trajectories_adaptive`, wraps the result into `ExperimentResult`, and saves to BSON.

**Primary recommendation:** Implement `run_experiment()` as a function in the sweep script (not in the package) that accepts a pre-built `HamHam`, `TrottTrott`, jump operators, and experiment parameters (beta, sigma, delta, db_type), then calls `run_trajectories_adaptive` and saves the result. The sweep loop constructs Hamiltonians per (n, beta) pair (9 combinations), creates Trotter objects, builds jump operators, then runs 3 experiments per pair (KMS, GNS@sigma=1/beta, GNS@sigma=0.5/beta).

## Standard Stack

### Core

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| `QuantumFurnace.jl` | current | All simulation infrastructure: configs, trajectories, convergence, results | The package under use |
| `BSON` | 0.3+ | Experiment result serialization | Already used by save_experiment/load_experiment |
| `LinearAlgebra` (stdlib) | Julia 1.11+ | Hermitian, tr, eigvals, I | Standard Julia scientific computing |
| `Dates` (stdlib) | Julia 1.11+ | Timestamps in metadata | Used by _capture_metadata |

### Supporting

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| `Printf` (stdlib) | Julia 1.11+ | Progress logging during sweep | Print experiment progress |
| `Pkg` (stdlib) | Julia 1.11+ | Project root detection for output paths | Used by load_hamiltonian pattern |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| TrotterDomain for all experiments | EnergyDomain | Locked decision: TrotterDomain. EnergyDomain skips Trotter approximation error but doesn't match quantum circuit implementation. |
| Standalone sweep script | Package API function | Locked decision: sweep script in experiments/, not part of package API |
| BSON serialization | JLD2 or Arrow | BSON is the existing pattern; ExperimentResult already has Dict-based BSON serialization |

## Architecture Patterns

### Pattern 1: Hamiltonian Construction for Uniform XXX Chain

**What:** Build the 1D Heisenberg XXX chain with J=1.0 uniform coupling and periodic boundaries. No disorder.

**When to use:** At the start of the sweep, once per (n, beta) pair.

**Implementation:**
```julia
function build_heisenberg_xxx(num_qubits::Int, beta::Float64)
    terms = [[X, X], [Y, Y], [Z, Z]]
    coeffs = [1.0, 1.0, 1.0]
    return HamHam(terms, coeffs, num_qubits, beta; periodic=true)
end
```

This uses the existing `HamHam(terms, coeffs, num_qubits, beta; periodic)` constructor in `hamiltonian.jl`, which:
1. Constructs the Hamiltonian matrix via `_construct_base_ham`
2. Rescales spectrum to [0, 0.45] via `_rescaling_and_shift_factors`
3. Computes eigendecomposition, Bohr frequencies, bohr_dict
4. Computes Gibbs state at the given beta

Note: The Gibbs state is computed from the **rescaled** eigenvalues. The rescaling factor and shift are stored in HamHam for provenance. Since all experiments at the same (n, beta) share the same Hamiltonian, the eigendecomposition is computed once and reused.

Important: A new HamHam must be constructed for each beta value because the Gibbs state depends on beta. The HamHam stores `gibbs` as a field computed at construction time.

### Pattern 2: TrotterDomain Setup

**What:** Create TrottTrott and jump operators in the correct basis for TrotterDomain experiments.

**When to use:** For each (n, beta) pair, after constructing the Hamiltonian.

**Implementation:**
```julia
function build_trotter_system(hamiltonian::HamHam, num_qubits::Int, t0::Float64, num_trotter_steps::Int)
    trotter = TrottTrott(hamiltonian, t0, num_trotter_steps)

    # Jump operators: single-site Paulis (X, Y, Z) on each site
    jump_paulis = [[X], [Y], [Z]]
    n_jumps = length(jump_paulis) * num_qubits
    jump_normalization = sqrt(n_jumps)

    jumps = JumpOp[]
    for pauli in jump_paulis
        for site in 1:num_qubits
            jump_op = Matrix(pad_term(pauli, num_qubits, site)) ./ jump_normalization
            jump_in_eigen = hamiltonian.eigvecs' * jump_op * hamiltonian.eigvecs
            orthogonal = (jump_op == transpose(jump_op))
            herm = (jump_op == jump_op')
            push!(jumps, JumpOp(jump_op, jump_in_eigen, orthogonal, herm))
        end
    end

    return trotter, jumps
end
```

**Critical TrotterDomain detail:** The jump operators' `in_eigenbasis` field stores the operator in the Hamiltonian eigenbasis. When `run_trajectories` is called with TrotterDomain, `build_trajectoryframework` internally transforms jumps to the Trotter eigenbasis via `transform_jumps_to_basis(jumps, trotter.eigvecs)`. The `rho_mean` returned by `run_trajectories` is in the **Trotter eigenbasis**. Therefore:
- Observables must be built via `build_convergence_observables_trotter(hamiltonian, trotter, num_qubits)` (Trotter basis)
- Gibbs state for trace distance must be `_gibbs_in_trotter_basis(hamiltonian, trotter)` (Trotter basis)
- This is the same pattern used in the Phase 16 convergence tests

### Pattern 3: Per-Experiment Function

**What:** A function that takes pre-built infrastructure and runs one adaptive trajectory experiment.

**When to use:** Called 27 times in the sweep loop.

**Implementation:**
```julia
function run_experiment(;
    jumps::Vector{JumpOp},
    hamiltonian::HamHam,
    trotter::TrottTrott,
    num_qubits::Int,
    beta::Float64,
    sigma::Float64,
    delta::Float64,
    db_type::Symbol,  # :KMS or :GNS
    mixing_time::Float64,
    batch_size::Int,
    n_max::Int,
    convergence_threshold::Float64,
    patience::Int,
    window_size::Int,
    seed::Int,
    output_path::String,
)
    # 1. Build config
    config = if db_type == :KMS
        ThermalizeConfig(
            num_qubits=num_qubits, with_coherent=true, with_linear_combination=true,
            domain=TrotterDomain(), beta=beta, sigma=sigma,
            a=beta/30.0, b=0.4,
            num_energy_bits=12, w0=0.05,
            t0=2pi / (2^12 * 0.05), num_trotter_steps_per_t0=10,
            mixing_time=mixing_time, delta=delta,
        )
    else  # :GNS
        ThermalizeConfigGNS(
            num_qubits=num_qubits, with_coherent=false, with_linear_combination=true,
            domain=TrotterDomain(), beta=beta, sigma=sigma,
            a=beta/30.0, b=0.4,
            num_energy_bits=12, w0=0.05,
            t0=2pi / (2^12 * 0.05), num_trotter_steps_per_t0=10,
            mixing_time=mixing_time, delta=delta,
        )
    end

    # 2. Gibbs and observables in Trotter basis
    gibbs_trotter = _gibbs_in_trotter_basis(hamiltonian, trotter)
    observables, obs_names = build_convergence_observables_trotter(hamiltonian, trotter, num_qubits)

    # 3. Initial state
    dim = 2^num_qubits
    psi0 = zeros(ComplexF64, dim)
    psi0[1] = 1.0  # |0...0> in Trotter eigenbasis

    # 4. Run adaptive trajectories
    wall_t0 = time()
    traj_result, conv_data = run_trajectories_adaptive(
        jumps, config, psi0, hamiltonian;
        gibbs=gibbs_trotter, observables=observables, observable_names=obs_names,
        batch_size=batch_size, n_max=n_max,
        convergence_threshold=convergence_threshold, patience=patience,
        window_size=window_size, seed=seed,
        trotter=trotter, total_time=mixing_time, delta=delta,
    )
    wall_time = time() - wall_t0

    # 5. Build ExperimentResult and save
    ham_params = QuantumFurnace._extract_hamiltonian_params(hamiltonian)
    metadata = QuantumFurnace._capture_metadata(; wall_time_seconds=wall_time,
        extra=Dict{Symbol,Any}(
            :convergence => QuantumFurnace._convergence_to_dict(conv_data),
            :sigma_rule => sigma == 1.0/beta ? "1/beta" : "0.5/beta",
        ))

    result = ExperimentResult(config, traj_result, ham_params, metadata)
    save_experiment(result, output_path)

    return result, conv_data
end
```

### Pattern 4: Sweep Script Structure

**What:** The top-level sweep script that iterates over the parameter grid.

```julia
# experiments/run_sweep.jl
using QuantumFurnace
using LinearAlgebra
using Printf
using Dates

# Parameter grid
system_sizes = [4, 6, 8]
betas = [5.0, 10.0, 20.0]
# For each (n, beta): KMS, GNS@sigma=1/beta, GNS@sigma=0.5/beta

# Grid parameters (shared across all experiments)
const DELTA = 0.01
const N_MAX = 10_000
const NUM_ENERGY_BITS = 12
const W0 = 0.05
const T0 = 2pi / (2^NUM_ENERGY_BITS * W0)
const NUM_TROTTER_STEPS_PER_T0 = 10

# Experiment output directory
output_dir = joinpath(dirname(Pkg.project().path), "experiments")
mkpath(output_dir)

failed_experiments = String[]

for n in system_sizes
    for beta in betas
        # Build Hamiltonian and Trotter system (shared across 3 experiments)
        hamiltonian = build_heisenberg_xxx(n, beta)
        trotter, jumps = build_trotter_system(hamiltonian, n, T0, NUM_TROTTER_STEPS_PER_T0)

        # 3 experiments per (n, beta)
        experiments = [
            (:KMS, 1.0/beta, "kms_n$(n)_beta$(Int(beta)).bson"),
            (:GNS, 1.0/beta, "gns_sigma1_n$(n)_beta$(Int(beta)).bson"),
            (:GNS, 0.5/beta, "gns_sigma05_n$(n)_beta$(Int(beta)).bson"),
        ]

        for (db_type, sigma, filename) in experiments
            label = "$(db_type) n=$(n) beta=$(Int(beta)) sigma=$(round(sigma; digits=4))"
            @printf("[%s] Starting: %s\n", Dates.format(now(), "HH:MM:SS"), label)

            try
                result, conv_data = run_experiment(;
                    jumps=jumps, hamiltonian=hamiltonian, trotter=trotter,
                    num_qubits=n, beta=beta, sigma=sigma, delta=DELTA,
                    db_type=db_type, mixing_time=..., batch_size=...,
                    n_max=N_MAX, convergence_threshold=..., patience=...,
                    window_size=..., seed=42,
                    output_path=joinpath(output_dir, filename),
                )

                status = conv_data.converged ? "CONVERGED" : "HIT CAP"
                @printf("  -> %s  td=%.4f  n_traj=%d  wall=%.1fs\n",
                    status, conv_data.trace_distances[end],
                    traj_result.n_trajectories, metadata[:wall_time_seconds])
            catch e
                @printf("  -> FAILED: %s\n", e)
                push!(failed_experiments, label)
            end
        end
    end
end

if !isempty(failed_experiments)
    @printf("\nWARNING: %d experiments failed:\n", length(failed_experiments))
    for label in failed_experiments
        @printf("  - %s\n", label)
    end
end
```

### Pattern 5: Saving Experiment Data to BSON

**What:** Each experiment result is saved as an individual BSON file containing the full config, trajectory result (density matrix + stats), Hamiltonian parameters, convergence data, and metadata.

**Key point:** The `ExperimentResult` struct + `save_experiment` function (Phase 15) handles all serialization. The `ConvergenceData` is stored in the metadata dict as a nested Dict via `_convergence_to_dict`. This means each BSON file contains:
- Full config (reconstructable via `_reconstruct_config`)
- `rho_mean` density matrix (final trajectory average)
- `n_trajectories`, `seed`
- Hamiltonian parameters (base_terms, base_coeffs, periodic, shift, rescaling_factor)
- Convergence data (trace distances, observable values, adaptive diagnostics)
- Metadata (timestamp, git hash, Julia version, wall time, thread count)

**Loading:** `load_experiment(path)` reconstructs the full `ExperimentResult` including typed config.

### Anti-Patterns to Avoid

- **Constructing Hamiltonians inside run_experiment:** The locked decision says the function accepts a pre-built Hamiltonian. This allows future reuse with different Hamiltonian models.
- **Using EnergyDomain:** The locked decision specifies TrotterDomain. All observables and Gibbs state must be in the Trotter eigenbasis.
- **Comparing trajectory rho_mean directly to hamiltonian.gibbs:** TrotterDomain `rho_mean` is in Trotter eigenbasis; `hamiltonian.gibbs` is in energy eigenbasis. Must use `_gibbs_in_trotter_basis(hamiltonian, trotter)` for the Gibbs reference.
- **Using `with_coherent=true` for GNS configs:** The GNS constructor enforces `with_coherent=false`. This is correct -- GNS omits the coherent B term by design.
- **Hardcoding sigma as a constant:** sigma varies across experiments (1/beta and 0.5/beta), and also varies with beta.
- **Pre-allocating all 27 Hamiltonians:** Each (n, beta) pair needs its own HamHam (because Gibbs depends on beta), but Hamiltonians can be constructed on-the-fly and garbage collected.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Adaptive trajectory sampling | Custom batch loop | `run_trajectories_adaptive(...)` | Phase 17 implementation handles convergence detection, hard cap, seed management |
| BSON experiment serialization | Custom BSON writing | `save_experiment(ExperimentResult, path)` | Phase 15 handles Dict conversion, companion .txt, type safety |
| Gibbs state in Trotter basis | Manual basis transform | `QuantumFurnace._gibbs_in_trotter_basis(ham, trotter)` | Already handles the V_T' V_H gibbs V_H' V_T transform |
| Trotter-basis observables | Manual eigenbasis transform | `build_convergence_observables_trotter(ham, trotter, n)` | Builds ZZ correlations + H in Trotter basis |
| Metadata capture | Manual git hash/timestamp | `QuantumFurnace._capture_metadata(; wall_time, extra)` | Auto-captures Julia version, git hash, timestamp, threads |
| Hamiltonian param extraction | Manual field copying | `QuantumFurnace._extract_hamiltonian_params(ham)` | Extracts minimal provenance dict |
| Heisenberg chain construction | Custom matrix building | `HamHam(terms, coeffs, n, beta; periodic=true)` | Handles rescaling, eigendecomposition, Gibbs |

**Key insight:** Phase 18 is an *orchestration* phase. All computational building blocks exist. The new code is: (1) a helper to build uniform XXX Hamiltonians, (2) a `run_experiment()` wrapper, and (3) a sweep loop. Total new code: ~150-200 lines in a standalone script.

## Common Pitfalls

### Pitfall 1: TrotterDomain Basis Mismatch for Gibbs State and Observables

**What goes wrong:** Comparing `rho_mean` (in Trotter eigenbasis) against `hamiltonian.gibbs` (in energy eigenbasis) gives spurious large trace distances (~0.83 for n=3 in Quick-20).

**Why it happens:** `run_trajectories` with TrotterDomain returns `rho_mean` in the Trotter eigenbasis because `build_trajectoryframework` transforms jumps to Trotter basis. The Gibbs state `hamiltonian.gibbs` is always in the energy eigenbasis.

**How to avoid:** Use `_gibbs_in_trotter_basis(hamiltonian, trotter)` for the Gibbs reference and `build_convergence_observables_trotter(hamiltonian, trotter, num_qubits)` for observables. This is the pattern established in Quick-20 and Phase 16.

**Warning signs:** Trace distance to Gibbs is >0.5 at the first checkpoint (should be ~0.3-0.5 for a random initial state, not >0.8).

### Pitfall 2: delta_eff Overflow for Larger Systems

**What goes wrong:** `build_trajectoryframework` asserts `delta_eff < 1.0`. With delta=0.01 and more jump operators, delta_eff could approach the limit.

**Why it happens:** `delta_eff = delta * n_jumps`. For the XXX Heisenberg chain with single-site Paulis:
- n=4: 12 jumps, delta_eff = 0.12 (safe)
- n=6: 18 jumps, delta_eff = 0.18 (safe)
- n=8: 24 jumps, delta_eff = 0.24 (safe)

**How to avoid:** With delta=0.01, all system sizes are safe. The constraint is delta < 1/n_jumps. For n=8: delta < 1/24 = 0.042. The locked delta=0.01 is well within bounds.

**Warning signs:** Assertion error from build_trajectoryframework on the first trajectory batch.

### Pitfall 3: KMS-vs-GNS Comparison Conflates Fixed-Point Error with Mixing Rate

**What goes wrong:** KMS converges to trace distance ~0 (exact Gibbs fixed point), GNS saturates at the approximation gap (~0.08 for sigma=0.1). Naive comparison says "KMS is better" but ignores that GNS may mix *faster* to its own fixed point.

**Why it happens:** KMS Lindbladian has the exact Gibbs state as fixed point (with coherent B term). GNS fixed point deviates from Gibbs by O(sigma). Trace distance to Gibbs conflates two effects: mixing rate and fixed-point quality.

**How to avoid:** The experiment design already handles this: both trace distance to Gibbs (via `run_trajectories_adaptive` gibbs reference) and convergence behavior (trace distance curve over batches) are recorded. The convergence data in each BSON file contains the full trace distance trajectory. Analysis (Phase 19+) can separate mixing rate from fixed-point quality by examining the plateau level.

**Warning signs:** GNS trace distance curve flattens at a nonzero value while KMS reaches near-zero.

### Pitfall 4: Insufficient Mixing Time for Large Beta

**What goes wrong:** At high beta (=20), the spectral gap shrinks and the system needs more time steps to mix. If mixing_time is too short, trajectory averages have not converged to the fixed point.

**Why it happens:** Mixing time scales roughly as 1/spectral_gap, which grows with beta (slower relaxation at low temperature). At beta=20, the system is colder and mixes more slowly.

**How to avoid:** Scale mixing_time with beta. A reasonable heuristic: `mixing_time = C * beta` where C is chosen to give sufficient steps. With delta=0.01, the number of steps is mixing_time/delta. For beta=20, mixing_time=20.0 gives 2000 steps. For beta=5, mixing_time=5.0 gives 500 steps. This matches the Phase 14 baseline of mixing_time=5.0 at beta=10.0 (50 steps/beta unit).

**Warning signs:** Adaptive sampling hits the n_max cap for all high-beta experiments; trace distances are still decreasing at termination.

### Pitfall 5: Forgetting to Add experiments/ to .gitignore

**What goes wrong:** Large BSON data files get committed to git, bloating the repository.

**Why it happens:** The locked decision says `experiments/` is gitignored, but the current `.gitignore` does not include it.

**How to avoid:** Add `/experiments/` to `.gitignore` as part of the implementation.

**Warning signs:** `git status` shows BSON files as untracked.

### Pitfall 6: Using `psi0 = zeros(...)` Instead of Computational Basis State

**What goes wrong:** If psi0 is the zero vector, all trajectories produce NaN density matrices.

**Why it happens:** psi0 must be a normalized state vector. The common pattern is `psi0[1] = 1.0` for the |0...0> computational basis state.

**How to avoid:** Always set `psi0 = zeros(ComplexF64, dim); psi0[1] = 1.0`.

**Warning signs:** NaN in trace distance or density matrix entries.

## Code Examples

### Uniform XXX Heisenberg Chain Construction

```julia
# Source: hamiltonian.jl HamHam(terms, coeffs, num_qubits, beta; periodic) constructor
function build_heisenberg_xxx(num_qubits::Int, beta::Float64)
    terms = Vector{Vector{Matrix{ComplexF64}}}([[X, X], [Y, Y], [Z, Z]])
    coeffs = [1.0, 1.0, 1.0]
    return HamHam(terms, coeffs, num_qubits, beta; periodic=true)
end
```

### TrotterDomain System Setup

```julia
# Source: test_helpers.jl make_test_trotter() and make_test_system() patterns
function build_trotter_system(hamiltonian::HamHam, num_qubits::Int)
    t0 = 2pi / (2^12 * 0.05)  # shared grid parameter
    trotter = TrottTrott(hamiltonian, t0, 10)  # 10 Trotter steps per t0

    jump_paulis = [[X], [Y], [Z]]
    n_jumps = length(jump_paulis) * num_qubits
    norm_factor = sqrt(n_jumps)

    jumps = JumpOp[]
    for pauli in jump_paulis
        for site in 1:num_qubits
            op = Matrix(pad_term(pauli, num_qubits, site)) ./ norm_factor
            in_eigen = hamiltonian.eigvecs' * op * hamiltonian.eigvecs
            push!(jumps, JumpOp(op, in_eigen, op == transpose(op), op == op'))
        end
    end

    return trotter, jumps
end
```

### TrotterDomain Gibbs State and Observables

```julia
# Source: convergence.jl _gibbs_in_trotter_basis and build_convergence_observables_trotter
gibbs_trotter = QuantumFurnace._gibbs_in_trotter_basis(hamiltonian, trotter)
observables, obs_names = build_convergence_observables_trotter(hamiltonian, trotter, num_qubits)
```

### Config Construction for KMS and GNS

```julia
# KMS config (with coherent B term -> exact Gibbs fixed point)
kms_config = ThermalizeConfig(
    num_qubits=n, with_coherent=true, with_linear_combination=true,
    domain=TrotterDomain(), beta=beta, sigma=1.0/beta,
    a=beta/30.0, b=0.4,
    num_energy_bits=12, w0=0.05, t0=2pi/(2^12*0.05), num_trotter_steps_per_t0=10,
    mixing_time=mixing_time, delta=0.01,
)

# GNS config at sigma=1/beta (approximate Gibbs fixed point, no B term)
gns_config = ThermalizeConfigGNS(
    num_qubits=n, with_coherent=false, with_linear_combination=true,
    domain=TrotterDomain(), beta=beta, sigma=1.0/beta,
    a=beta/30.0, b=0.4,
    num_energy_bits=12, w0=0.05, t0=2pi/(2^12*0.05), num_trotter_steps_per_t0=10,
    mixing_time=mixing_time, delta=0.01,
)

# GNS config at sigma=0.5/beta (tighter approximation, more expensive)
gns_config_half = ThermalizeConfigGNS(
    num_qubits=n, with_coherent=false, with_linear_combination=true,
    domain=TrotterDomain(), beta=beta, sigma=0.5/beta,
    a=beta/30.0, b=0.4,
    num_energy_bits=12, w0=0.05, t0=2pi/(2^12*0.05), num_trotter_steps_per_t0=10,
    mixing_time=mixing_time, delta=0.01,
)
```

### Calling run_trajectories_adaptive

```julia
# Source: convergence.jl run_trajectories_adaptive
psi0 = zeros(ComplexF64, 2^num_qubits)
psi0[1] = 1.0

traj_result, conv_data = run_trajectories_adaptive(
    jumps, config, psi0, hamiltonian;
    gibbs=gibbs_trotter,
    observables=observables,
    observable_names=obs_names,
    batch_size=200,
    n_max=10_000,
    convergence_threshold=0.01,
    patience=3,
    min_batches=5,
    window_size=3,
    seed=42,
    trotter=trotter,
    total_time=config.mixing_time,
    delta=config.delta,
)
```

### Saving Experiment Results

```julia
# Source: results.jl save_experiment
ham_params = QuantumFurnace._extract_hamiltonian_params(hamiltonian)
metadata = QuantumFurnace._capture_metadata(;
    wall_time_seconds=wall_time,
    extra=Dict{Symbol,Any}(
        :convergence => QuantumFurnace._convergence_to_dict(conv_data),
        :sigma_rule => "1/beta",
    )
)
result = ExperimentResult(config, traj_result, ham_params, metadata)
save_experiment(result, joinpath(output_dir, filename))
```

## Discretion Recommendations

### Mixing Time: Scale with beta

**Recommendation:** `mixing_time = 2.0 * beta`

**Reasoning:**
- Phase 14 used mixing_time=5.0 at beta=10.0, which gave 500 steps at delta=0.01 and achieved 0.029 trace distance to GNS fixed point. That is a ratio of 0.5 * beta.
- For a parameter sweep, we want enough mixing for all beta values. The spectral gap shrinks with increasing beta (colder system mixes slower).
- With `mixing_time = 2.0 * beta`:
  - beta=5: mixing_time=10.0, steps=1000 (generous)
  - beta=10: mixing_time=20.0, steps=2000 (generous)
  - beta=20: mixing_time=40.0, steps=4000 (should be sufficient)
- This is conservative. For n=8 at beta=20, each trajectory has 4000 steps with dim=256 -- still computationally feasible.
- 2.0 * beta ensures we have at least 2x the relaxation time scale, which is generally sufficient for convergence.

**Confidence:** MEDIUM -- the spectral gap depends on both n and beta. For n=8 at beta=20, the gap may be smaller than at n=3 at beta=10 (Phase 14 baseline). If experiments hit the n_max cap, mixing_time can be increased.

### Batch Size: 200

**Recommendation:** `batch_size = 200`

**Reasoning:**
- This is the default for `run_trajectories_adaptive` (Phase 17 locked decision).
- With n_max=10000 and batch_size=200, max_batches=50. This gives plenty of checkpoints for convergence detection.
- For n=8 (dim=256), each batch of 200 trajectories is computationally moderate.
- The windowed convergence detection with window_size=3 needs 2*3=6 minimum checkpoints, which is well under 50.

**Confidence:** HIGH -- matches Phase 17 defaults and provides good checkpoint granularity.

### Convergence Thresholds

**Recommendation:**
- `convergence_threshold = 0.01` (1% relative change in trace distance)
- `patience = 3` (3 consecutive stable checks)
- `window_size = 3` (compare mean of 3 recent vs 3 previous batches)
- `min_batches = 5` (burn-in before checking)

**Reasoning:**
- These are the Phase 17 defaults, validated in testing.
- The 1% threshold means trace distance has stabilized to within 1% relative change between windows of 3 batches each.
- With n_max=10000 and batch_size=200, we have 50 batches. After effective_min=max(5, 6)=6 burn-in batches (1200 trajectories), convergence checking begins. With patience=3, convergence requires 3 consecutive stable checks, so earliest convergence is at batch 9 (1800 trajectories).
- This is appropriate: with n_max=10000, the statistical error at 10000 trajectories is 1/sqrt(10000) = 0.01, matching the 0.01 delta systematic error.

**Confidence:** HIGH -- Phase 17 defaults are well-tested and the n_max=10000 cap ensures subdominant statistical error.

### File Naming

**Recommendation:**
- KMS: `kms_n{n}_beta{beta}.bson` (e.g., `kms_n4_beta5.bson`)
- GNS at sigma=1/beta: `gns_sigma1_n{n}_beta{beta}.bson` (e.g., `gns_sigma1_n6_beta10.bson`)
- GNS at sigma=0.5/beta: `gns_sigma05_n{n}_beta{beta}.bson` (e.g., `gns_sigma05_n8_beta20.bson`)

This matches the CONTEXT.md examples and makes files glob-friendly. Beta is always integer-valued (5, 10, 20).

**Confidence:** HIGH -- follows CONTEXT.md examples exactly.

### Seed Strategy

**Recommendation:** Fixed seed=42 for all experiments (deterministic, reproducible).

**Reasoning:** The sweep should be fully deterministic for reproducibility. Using the same base seed across all experiments is fine because each experiment has its own independent trajectory execution. The adaptive batch loop uses `batch_seed = actual_seed + n_total` internally (Phase 16/17 pattern), so seed sequences are non-overlapping.

**Confidence:** HIGH -- follows existing test patterns.

### Progress Logging

**Recommendation:** Print one line per experiment start and one line per experiment completion, including convergence status, final trace distance, trajectory count, and wall time.

```
[10:15:23] Starting: KMS n=4 beta=5 sigma=0.2000
  -> CONVERGED  td=0.0312  n_traj=3600  wall=12.3s
[10:15:35] Starting: GNS n=4 beta=5 sigma=0.2000
  -> CONVERGED  td=0.0891  n_traj=4200  wall=14.1s
```

This gives sufficient visibility without overwhelming output.

## State of the Art

| What Phase 18 Needs | Already Exists (Phase) | Status |
|---|---|---|
| Adaptive trajectory sampling | `run_trajectories_adaptive` (Phase 17) | READY |
| Convergence data struct | `ConvergenceData` (Phase 16+17) | READY |
| BSON experiment serialization | `ExperimentResult`, `save_experiment` (Phase 15) | READY |
| GNS trajectory support | `ThermalizeConfigGNS`, `run_trajectories` (Phase 14) | READY |
| KMS trajectory support | `ThermalizeConfig`, `run_trajectories` (Phase 13) | READY |
| Trotter-basis Gibbs state | `_gibbs_in_trotter_basis` (Phase 16) | READY |
| Trotter-basis observables | `build_convergence_observables_trotter` (Phase 16) | READY |
| Hamiltonian construction (uniform) | `HamHam(terms, coeffs, n, beta)` (Phase 1) | READY |
| TrottTrott construction | `TrottTrott(ham, t0, steps)` (original code) | READY |
| Metadata capture | `_capture_metadata` (Phase 15) | READY |

**What's new in Phase 18:**
- `experiments/` directory (gitignored)
- `experiments/run_sweep.jl` standalone script (~150-200 lines)
- Helper functions in the script: `build_heisenberg_xxx`, `build_trotter_system`, `run_experiment`
- `.gitignore` update for `/experiments/`

## Key Technical Details

### delta_eff Constraints

| n | dim | n_jumps | delta_eff (delta=0.01) | Constraint (< 1.0) |
|---|-----|---------|------------------------|---------------------|
| 4 | 16 | 12 | 0.12 | OK |
| 6 | 64 | 18 | 0.18 | OK |
| 8 | 256 | 24 | 0.24 | OK |

### Trajectory Step Counts (mixing_time = 2*beta, delta=0.01)

| beta | mixing_time | steps/trajectory | n_max batches (200/batch) | max total trajectories |
|------|-------------|------------------|---------------------------|------------------------|
| 5 | 10.0 | 1000 | 50 | 10000 |
| 10 | 20.0 | 2000 | 50 | 10000 |
| 20 | 40.0 | 4000 | 50 | 10000 |

### Memory Estimates

| n | dim | rho_mean size | Per-operator Kraus | Total framework memory |
|---|-----|---------------|--------------------|------------------------|
| 4 | 16 | 16x16 = 2 KB | ~2 KB x 12 = 24 KB | ~100 KB |
| 6 | 64 | 64x64 = 32 KB | ~32 KB x 18 = 576 KB | ~2 MB |
| 8 | 256 | 256x256 = 512 KB | ~512 KB x 24 = 12 MB | ~50 MB |

The n=8 case uses ~50 MB per trajectory framework. With threading, this scales linearly with thread count but is still modest.

### Expected Wallclock Times (rough estimates)

Based on Phase 14 baseline (n=3, dim=8, 1000 trajectories, 500 steps, ~5s):
- n=4 (dim=16): ~4x slower per step -> ~20s for 1000 traj at 1000 steps
- n=6 (dim=64): ~64x slower per step -> ~5min for 1000 traj at 2000 steps
- n=8 (dim=256): ~1000x slower per step -> ~80min for 1000 traj at 4000 steps

For n=8 at beta=20, each batch of 200 trajectories takes ~16min. With 50 batches max, worst case ~13h. If convergence is reached early (e.g., 20 batches), ~5h.

Total sweep (27 experiments): most time is in n=8 experiments. Rough estimate: 1-2 hours for n=4,6 experiments, 5-20 hours for n=8 experiments. Total: ~6-24 hours depending on convergence speed.

### Accessing Internal Functions

The sweep script uses several internal (unexported) functions from QuantumFurnace:
- `QuantumFurnace._gibbs_in_trotter_basis`
- `QuantumFurnace._extract_hamiltonian_params`
- `QuantumFurnace._capture_metadata`
- `QuantumFurnace._convergence_to_dict`

These are accessed via module-qualified names (e.g., `QuantumFurnace._gibbs_in_trotter_basis`). This is acceptable for a standalone script that is not part of the public API.

### KMS sigma Parameter

For KMS experiments, sigma=1/beta is used (matching the GNS@sigma=1/beta experiments). The KMS transition function uses a shifted argument `w + beta*sigma^2/2`, while GNS uses the unshifted argument `w`. As documented in PITFALLS.md Pitfall 5, using the same sigma for both KMS and GNS means the two lines have different effective behavior. This is the "same sigma" comparison protocol (Option A from PITFALLS.md). The locked decision specifies sigma=1/beta for KMS and GNS at the same beta, which is the natural choice.

## Open Questions

1. **Exact wallclock time for n=8 experiments**
   - What we know: Rough estimates suggest 5-20 hours for each n=8 experiment depending on convergence.
   - What's unclear: The actual spectral gap for n=8 uniform XXX Heisenberg determines mixing speed. Disordered systems may have different gaps.
   - Recommendation: Start the sweep and observe. If n=8 experiments are too slow, the mixing_time can be reduced or n_max lowered for initial exploration.

2. **Whether mixing_time = 2*beta is sufficient for all (n, beta) combinations**
   - What we know: Phase 14 used 5.0 at beta=10 for n=3 and converged well. The spectral gap depends on both n and beta.
   - What's unclear: The spectral gap for n=8 at beta=20 with the uniform XXX chain.
   - Recommendation: Use 2*beta as the starting point. If experiments hit n_max without converging, increase to 4*beta for the problematic combinations. The adaptive sampling's `converged=false` flag will identify these cases.

3. **Whether `a = beta/30.0` is appropriate for all beta values**
   - What we know: Phase 14 used a=beta/30.0 (=0.333 at beta=10). This is the "Smooth Metro" parameter from the test fixtures.
   - What's unclear: Whether a=beta/30 is optimal for beta=5 (a=0.167) and beta=20 (a=0.667). The transition function shape depends on a.
   - Recommendation: Use a=beta/30 consistently across all experiments (matches existing test patterns). This is a research parameter, not an optimization parameter.

## Sources

### Primary (HIGH confidence)
- **QuantumFurnace.jl codebase** -- Direct analysis of:
  - `src/hamiltonian.jl`: HamHam constructors, find_ideal_heisenberg, Heisenberg chain pattern
  - `src/trotter_domain.jl`: TrottTrott constructor, Trotter eigendecomposition
  - `src/structs.jl`: Config type hierarchy, ThermalizeConfig/ThermalizeConfigGNS with all fields
  - `src/convergence.jl`: run_trajectories_adaptive signature, ConvergenceData, observable builders
  - `src/results.jl`: ExperimentResult, save_experiment, _convergence_to_dict, _capture_metadata
  - `src/trajectories.jl`: run_trajectories, delta_eff constraint, TrajectoryFramework
  - `src/furnace.jl`: construct_lindbladian TrotterDomain basis handling
  - `src/misc_tools.jl`: load_hamiltonian, validate_config!
  - `test/test_helpers.jl`: make_test_system, jump operator construction pattern
  - `hamiltonians/generate_hamiltonians.jl`: find_ideal_heisenberg usage pattern

- **Phase planning docs** -- Direct analysis of:
  - Phase 14 PLAN, SUMMARY, VERIFICATION: GNS trajectory validation, baseline gaps, convergence parameters
  - Phase 17 PLAN 01+02, SUMMARIES: Adaptive sampling implementation, API, defaults, test coverage
  - Phase 16 convergence tracking: Observable builders, Gibbs basis transforms
  - Quick-20 SUMMARY: TrotterDomain basis mismatch fix, GNS gap baselines

- **PITFALLS.md**: Pitfall 5 (KMS-vs-GNS sigma comparison), Pitfall 6 (fixed-point vs mixing rate), Pitfall 7 (adaptive premature termination)
- **STATE.md**: Prior decisions, Quick-20 baselines
- **REQUIREMENTS.md**: EXPT-01 through EXPT-04
- **ROADMAP.md**: Phase 18 success criteria

### Secondary (MEDIUM confidence)
- Wallclock time estimates: Based on Phase 14 baseline extrapolation (dim^3 scaling assumption for matrix operations)

### Tertiary (LOW confidence)
- None -- all findings verified against codebase analysis

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- No new dependencies; all infrastructure from Phases 13-17
- Architecture: HIGH -- Sweep script is pure orchestration of existing APIs
- Pitfalls: HIGH -- All identified from codebase analysis and documented prior issues (Quick-20, PITFALLS.md)
- Code examples: HIGH -- All patterns derived from existing working test/simulation code
- Discretion recommendations: HIGH for batch_size/convergence (matching defaults), MEDIUM for mixing_time (heuristic scaling)

**Research date:** 2026-02-16
**Valid until:** 2026-03-16 (stable codebase, no external dependencies)
