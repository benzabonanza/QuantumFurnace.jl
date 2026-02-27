# ============================================================================
# Trajectory simulation with observable tracking
# ============================================================================
#
# Demonstrates the run_trajectory API with Config{Trajectory} and
# observable-mode trajectory averaging.
#
# Usage:
#   julia --project=@. simulations/main_trajectory.jl

using Revise
includet("../src/QuantumFurnace.jl")
using .QuantumFurnace

using LinearAlgebra, Printf

function main()
    #* Config
    num_qubits = 4
    dim = 2^num_qubits
    beta = 10.0
    sigma = 0.8 / beta
    w_gamma = 1 / beta
    sigma_gamma = sqrt(2 * w_gamma / beta - sigma^2)

    # Smooth Metro
    a = 1 / 10
    b = 0.4
    eta = 0.0

    construction = KMS()
    with_linear_combination = true
    domain = EnergyDomain()
    num_energy_bits = 12
    w0 = 0.05
    t0 = 2pi / (2^num_energy_bits * w0)
    num_trotter_steps_per_t0 = 10

    mixing_time = 50.0
    delta = 0.005

    config = Config(
        sim = Trajectory(),
        domain = domain,
        construction = construction,
        num_qubits = num_qubits,
        with_linear_combination = with_linear_combination,
        beta = beta,
        sigma = sigma,
        gaussian_parameters = (w_gamma, sigma_gamma),
        a = a,
        b = b,
        num_energy_bits = num_energy_bits,
        w0 = w0,
        t0 = t0,
        eta = eta,
        num_trotter_steps_per_t0 = num_trotter_steps_per_t0,
        mixing_time = mixing_time,
        delta = delta,
    )

    #* Hamiltonian
    hamiltonian = load_hamiltonian("heis", num_qubits; beta=beta)

    #* Jumps
    jump_paulis = [[X], [Y], [Z]]
    num_of_jumps = length(jump_paulis) * num_qubits
    jump_normalization = sqrt(num_of_jumps)
    jumps::Vector{JumpOp} = []
    for pauli in jump_paulis
        for site in 1:num_qubits
            jump_op = pad_term(pauli, num_qubits, site) / jump_normalization
            basis_unitary = hamiltonian.eigvecs
            jump_op_in_eigenbasis = basis_unitary' * jump_op * basis_unitary
            orthogonal = (jump_op == transpose(jump_op))
            hermitian = (jump_op == jump_op')
            push!(jumps, JumpOp(jump_op, jump_op_in_eigenbasis, orthogonal, hermitian))
        end
    end

    #* Initial state: first computational basis state
    psi0 = zeros(ComplexF64, dim)
    psi0[1] = 1.0

    #* Build observables
    observables, observable_names = build_preset_trajectory_observables(hamiltonian, num_qubits)

    #* Run trajectory with observable tracking
    result = @time run_trajectory(
        jumps, config, hamiltonian, nothing;
        psi0 = psi0,
        ntraj = 1000,
        seed = 42,
        save_every = 10,
        observables = observables,
        observable_names = observable_names,
    )

    @printf("\nTrajectory results:\n")
    @printf("  N trajectories: %d\n", result.n_trajectories)
    @printf("  Seed:           %d\n", result.seed)
    @printf("  Wall time:      %.2f s\n", result.metadata[:wall_time_seconds])
    @printf("  rho_mean dim:   %dx%d\n", size(result.rho_mean)...)

    # Save
    # save_result(result, "results/trajectory_result.bson")
end

main()
