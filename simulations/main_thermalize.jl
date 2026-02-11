using Distributed
# addprocs(4, exeflags="--project=@.")

using Revise
includet("../src/QuantumFurnace.jl")
using .QuantumFurnace

@everywhere begin
    using Revise
    includet("../src/QuantumFurnace.jl")
    using .QuantumFurnace
end

using Pkg, LinearAlgebra, Random, Printf, SparseArrays, BSON, Arpack

function main()
    #* Config
    num_qubits = 4
    dim = 2^num_qubits
    beta = 10.  # 5, 10, 30
    sigma = 1 / beta
    w_gamma = 1 / beta
    sigma_gamma = 1 / beta

    #! change from 1 / beta
    sigma = 0.2 / beta  # w0 = 0.005, for broad enough time integrals in OFTs
    w_gamma = 1 / beta
    sigma_gamma = sqrt(2 * w_gamma / beta - sigma^2)

    # Smooth Metro
    a = 1 / 10
    b = 0.4
    eta = 0.0 

    # Kinky Metro 
    a = 0.0
    b = 0.0
    # eta = 0.002

    with_coherent = true
    with_linear_combination = false
    domain = TrotterDomain()
    num_energy_bits = 12 # 11
    w0 = 0.05
    max_E = w0 * 2^num_energy_bits / 2
    t0 = 2pi / (2^num_energy_bits * w0)  # Max time evolution pi / w0
    num_trotter_steps_per_t0 = 10

    # Thermalizing configs:
    mixing_time = 100.0 * 3 * num_qubits
    delta = 0.005

    # config = ThermalizeConfig(
    #     num_qubits = num_qubits, 
    #     with_coherent = with_coherent,
    #     with_linear_combination = with_linear_combination, 
    #     domain = domain,
    #     beta = beta,
    #     sigma = sigma,
    #     gaussian_parameters = (w_gamma, sigma_gamma),
    #     a = a,
    #     b = b,
    #     num_energy_bits = num_energy_bits,
    #     w0 = w0,
    #     t0 = t0,
    #     eta = eta,
    #     num_trotter_steps_per_t0 = num_trotter_steps_per_t0, 
    #     mixing_time = mixing_time,
    #     delta = delta,
    # )

    #* Approx GNS Config
    sigma_gamma = sqrt(2 * w_gamma / beta)
    config = ThermalizeConfigGNS(
        num_qubits = num_qubits,
        with_linear_combination = with_linear_combination,
        domain = domain,
        beta = beta,
        sigma = sigma,
        gaussian_parameters = (w_gamma, sigma_gamma),
        a = a,
        b = b,
        num_energy_bits = num_energy_bits,
        w0 = w0,
        t0 = t0,
        num_trotter_steps_per_t0 = num_trotter_steps_per_t0,
        mixing_time = mixing_time,
        delta = delta,
    )

    #* Hamiltonian    
    # Hamiltonian path
    project_root = Pkg.project().path |> dirname
    data_dir = joinpath(project_root, "hamiltonians")
    output_filename = join(["heis", "disordered", "periodic", "n$num_qubits"], "_") * ".bson"
    ham_path = joinpath(data_dir, output_filename)
    bson_hamiltonian_data = BSON.load(ham_path)
    hamiltonian = bson_hamiltonian_data[:hamiltonian]
    hamiltonian = finalize_hamham(hamiltonian, beta)
    
    initial_dm = Matrix{ComplexF64}(I(dim) / dim)
    @assert norm(real(tr(initial_dm)) - 1.) < 1e-15 "Trace is not 1.0"
    @assert norm(initial_dm - initial_dm') < 1e-15 "Not Hermitian"

    #* Trotter
    if domain == TrotterDomain()
        trotter = TrottTrott(hamiltonian, t0, num_trotter_steps_per_t0)
        trotter_error_T = compute_trotter_error(hamiltonian, trotter, 2^num_energy_bits * t0 / 2)
        gibbs_in_trotter = Hermitian(trotter.eigvecs' * gibbs_state(hamiltonian, beta) * trotter.eigvecs)
        @printf("Trotter is created.\n")
    else
        trotter = nothing
    end

    #* Jumps
    jump_paulis = [[X], [Y], [Z]]

    num_of_jumps = length(jump_paulis) * num_qubits
    jump_normalization = sqrt(num_of_jumps)
    jumps::Vector{JumpOp} = []
    for pauli in jump_paulis
            for site in 1:num_qubits
                jump_op = pad_term(pauli, num_qubits, site) / jump_normalization

                basis_unitary = (domain isa TrotterDomain) ? trotter.eigvecs : hamiltonian.eigvecs
                jump_op_in_eigenbasis = basis_unitary' * jump_op * basis_unitary

                orthogonal = (jump_op == transpose(jump_op))
                hermitian = (jump_op == jump_op')
                jump = JumpOp(jump_op,
                        jump_op_in_eigenbasis,
                        orthogonal,
                        hermitian)
                push!(jumps, jump)
            end
    end

    #* Thermalization
    alg_results = @time run_thermalization(jumps, config, initial_dm, hamiltonian; trotter=trotter)

    @printf("\n Last distance to Gibbs: %s\n", alg_results.distances_to_gibbs[end])
    @printf("Number of steps taken: %s\n", length(alg_results.time_steps))
    # plot(alg_results.time_steps, alg_results.distances_to_gibbs, label="Distance to Gibbs", xlabel="Time", ylabel="Distance", title="Distance to Gibbs over time")

    # Save
    # project_root = Pkg.project().path |> dirname
    # project_root = joinpath(project_root, "julia")  #! Omit this on cluster
    # results_dir = joinpath(project_root, "results")
    # output_filename = generate_filename(config)
    # full_path = joinpath(results_dir, output_filename)

    # println("Saving results to: ", full_path)
    # BSON.bson(full_path, Dict("results" => alg_results)) # Save as a dictionary
    # println("Save complete.")
end

if myid() == 1
    main()
end

# Load
# bson_data = BSON.load(full_path)
# loaded_results = bson_data["results"]