# ============================================================================
# KMS-vs-GNS Experiment Sweep
# ============================================================================
#
# Runs a 27-experiment parameter sweep comparing KMS (exact detailed balance)
# and GNS (approximate detailed balance) convergence across system sizes and
# inverse temperatures using TrotterDomain trajectory simulation.
#
# Usage:
#   julia --project=. experiments/run_sweep.jl          # all 27 experiments
#   julia --project=. experiments/run_sweep.jl 4        # n=4 only (9 experiments)
#   julia --project=. experiments/run_sweep.jl 4 6      # n=4,6 only (18 experiments)
#
# Output: BSON files in experiments/ (one per experiment).
# ============================================================================

using QuantumFurnace
using LinearAlgebra
using Printf
using Dates

# ============================================================================
# Shared constants
# ============================================================================

const DELTA = 0.01
const N_MAX = 10_000
const BATCH_SIZE = 200
const CONVERGENCE_THRESHOLD = 0.01
const PATIENCE = 3
const WINDOW_SIZE = 3
const MIN_BATCHES = 5
const NUM_ENERGY_BITS = 12
const W0 = 0.05
const T0 = 2pi / (2^NUM_ENERGY_BITS * W0)
const NUM_TROTTER_STEPS_PER_T0 = 10
const SEED = 42

# ============================================================================
# Helpers
# ============================================================================

"""
    build_heisenberg_xxx(num_qubits, beta) -> HamHam

Build a uniform 1D Heisenberg XXX chain with periodic boundaries and J=1.0.
"""
function build_heisenberg_xxx(num_qubits::Int, beta::Float64)
    terms = Vector{Vector{Matrix{ComplexF64}}}([[X, X], [Y, Y], [Z, Z]])
    coeffs = [1.0, 1.0, 1.0]
    return HamHam(terms, coeffs, num_qubits, beta; periodic=true)
end

"""
    build_trotter_system(hamiltonian, num_qubits) -> (TrottTrott, Vector{JumpOp})

Create TrottTrott and single-site Pauli jump operators for TrotterDomain experiments.
"""
function build_trotter_system(hamiltonian::HamHam, num_qubits::Int)
    trotter = TrottTrott(hamiltonian, T0, NUM_TROTTER_STEPS_PER_T0)

    jump_paulis = [[X], [Y], [Z]]
    n_jumps = length(jump_paulis) * num_qubits
    norm_factor = sqrt(n_jumps)

    jumps = JumpOp{Matrix{ComplexF64}}[]
    for pauli in jump_paulis
        for site in 1:num_qubits
            op = Matrix(pad_term(pauli, num_qubits, site)) ./ norm_factor
            in_eigen = hamiltonian.eigvecs' * op * hamiltonian.eigvecs
            push!(jumps, JumpOp(op, in_eigen, op == transpose(op), op == op'))
        end
    end

    return trotter, jumps
end

"""
    run_experiment(; kwargs...) -> (ExperimentResult, ConvergenceData, Float64)

Run a single adaptive trajectory experiment and save the result to BSON.
Returns (result, convergence_data, wall_time_seconds).
"""
function run_experiment(;
    jumps::Vector{JumpOp{Matrix{ComplexF64}}},
    hamiltonian::HamHam,
    trotter::TrottTrott,
    num_qubits::Int,
    beta::Float64,
    sigma::Float64,
    delta::Float64,
    db_type::Symbol,
    mixing_time::Float64,
    n_max::Int,
    seed::Int,
    output_path::String,
    label::String,
)
    # 1. Build config
    config = if db_type == :KMS
        ThermalizeConfig(
            num_qubits=num_qubits, with_coherent=true, with_linear_combination=true,
            domain=TrotterDomain(), beta=beta, sigma=sigma,
            a=beta/30.0, b=0.4,
            num_energy_bits=NUM_ENERGY_BITS, w0=W0, t0=T0,
            num_trotter_steps_per_t0=NUM_TROTTER_STEPS_PER_T0,
            mixing_time=mixing_time, delta=delta,
        )
    else  # :GNS
        ThermalizeConfigGNS(
            num_qubits=num_qubits, with_coherent=false, with_linear_combination=true,
            domain=TrotterDomain(), beta=beta, sigma=sigma,
            a=beta/30.0, b=0.4,
            num_energy_bits=NUM_ENERGY_BITS, w0=W0, t0=T0,
            num_trotter_steps_per_t0=NUM_TROTTER_STEPS_PER_T0,
            mixing_time=mixing_time, delta=delta,
        )
    end

    # 2. Gibbs and observables in Trotter basis (CRITICAL: not energy basis)
    gibbs_trotter = QuantumFurnace._gibbs_in_trotter_basis(hamiltonian, trotter)
    observables, obs_names = build_convergence_observables_trotter(hamiltonian, trotter, num_qubits)

    # 3. Initial state: |0...0> computational basis
    dim = 2^num_qubits
    psi0 = zeros(ComplexF64, dim)
    psi0[1] = 1.0

    # 4. Run adaptive trajectories
    wall_t0 = time()
    traj_result, conv_data = run_trajectories_adaptive(
        jumps, config, psi0, hamiltonian;
        gibbs=gibbs_trotter,
        observables=observables,
        observable_names=obs_names,
        batch_size=BATCH_SIZE,
        n_max=n_max,
        convergence_threshold=CONVERGENCE_THRESHOLD,
        patience=PATIENCE,
        min_batches=MIN_BATCHES,
        window_size=WINDOW_SIZE,
        seed=seed,
        trotter=trotter,
        total_time=mixing_time,
        delta=delta,
    )
    wall_time = time() - wall_t0

    # 5. Build ExperimentResult and save
    sigma_rule = sigma == 1.0/beta ? "1/beta" : "0.5/beta"
    ham_params = QuantumFurnace._extract_hamiltonian_params(hamiltonian)
    metadata = QuantumFurnace._capture_metadata(;
        wall_time_seconds=wall_time,
        extra=Dict{Symbol,Any}(
            :convergence => QuantumFurnace._convergence_to_dict(conv_data),
            :sigma_rule  => sigma_rule,
            :db_type     => string(db_type),
            :label       => label,
        ),
    )

    result = ExperimentResult(config, traj_result, ham_params, metadata)
    save_experiment(result, output_path)

    return result, conv_data, wall_time
end

# ============================================================================
# Main sweep
# ============================================================================

# Parameter grid
system_sizes = [4, 6, 8]
betas = [5.0, 10.0, 20.0]

# Command-line filter: optionally restrict system sizes
filter_sizes = isempty(ARGS) ? system_sizes : [parse(Int, a) for a in ARGS]

# Output directory
project_root = dirname(dirname(@__FILE__))
output_dir = joinpath(project_root, "experiments")
mkpath(output_dir)

sweep_start = time()
failed_experiments = String[]
total_experiments = length(filter_sizes) * length(betas) * 3
experiment_count = 0

@printf("\n=== KMS-vs-GNS Sweep ===\n")
@printf("System sizes: %s\n", string(filter_sizes))
@printf("Betas: %s\n", string(betas))
@printf("Total experiments: %d\n\n", total_experiments)

for n in filter_sizes
    for beta in betas
        mixing_time = 2.0 * beta

        @printf("[%s] Building Hamiltonian: n=%d, beta=%.0f\n",
            Dates.format(now(), "HH:MM:SS"), n, beta)
        hamiltonian = build_heisenberg_xxx(n, beta)
        trotter, jumps = build_trotter_system(hamiltonian, n)

        # 3 experiments per (n, beta) pair
        experiments = [
            (:KMS, 1.0/beta, "kms_n$(n)_beta$(Int(beta)).bson"),
            (:GNS, 1.0/beta, "gns_sigma1_n$(n)_beta$(Int(beta)).bson"),
            (:GNS, 0.5/beta, "gns_sigma05_n$(n)_beta$(Int(beta)).bson"),
        ]

        for (db_type, sigma, filename) in experiments
            experiment_count += 1
            label = "$(db_type) n=$(n) beta=$(Int(beta)) sigma=$(round(sigma; digits=4))"
            @printf("[%s] (%d/%d) Starting: %s\n",
                Dates.format(now(), "HH:MM:SS"), experiment_count, total_experiments, label)

            try
                result, conv_data, wall_time = run_experiment(;
                    jumps=jumps,
                    hamiltonian=hamiltonian,
                    trotter=trotter,
                    num_qubits=n,
                    beta=beta,
                    sigma=sigma,
                    delta=DELTA,
                    db_type=db_type,
                    mixing_time=mixing_time,
                    n_max=N_MAX,
                    seed=SEED,
                    output_path=joinpath(output_dir, filename),
                    label=label,
                )

                status = conv_data.converged ? "CONVERGED" : "HIT CAP"
                td = conv_data.trace_distances[end]
                n_traj = result.trajectory_result.n_trajectories
                @printf("  -> %s  td=%.4f  n_traj=%d  wall=%.1fs\n",
                    status, td, n_traj, wall_time)
            catch e
                @printf("  -> FAILED: %s\n", sprint(showerror, e))
                push!(failed_experiments, label)
            end
        end
    end
end

sweep_wall = time() - sweep_start
@printf("\n=== Sweep Complete ===\n")
@printf("Total wall time: %.1fs (%.1f min)\n", sweep_wall, sweep_wall / 60.0)
@printf("Experiments: %d/%d succeeded\n",
    total_experiments - length(failed_experiments), total_experiments)

if !isempty(failed_experiments)
    @printf("\nWARNING: %d experiments failed:\n", length(failed_experiments))
    for label in failed_experiments
        @printf("  - %s\n", label)
    end
end
