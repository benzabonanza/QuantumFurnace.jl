using Revise
using BenchmarkTools

if !isdefined(Main, :QuantumFurnace)
    includet("../src/QuantumFurnace.jl")
end

using .QuantumFurnace
using LinearAlgebra, Random

#* Config
num_qubits = 4
dim = 2^num_qubits
beta = 10.

# Smooth Metro
a = beta / 30. # a = beta / 50.
b = 0.4  # b = 0.5
eta = 0.0  # eta = 0.2

# Kinky Metro 
# a = 0.0
# b = 0.0
# eta = 0.002

with_coherent = true
with_linear_combination = true
domain = EnergyDomain()
num_energy_bits = 12  # 11
w0 = 0.05
max_E = w0 * 2^num_energy_bits / 2
t0 = 2pi / (2^num_energy_bits * w0)  # Max time evolution pi / w0
num_trotter_steps_per_t0 = 10

delta = 0.1

config = LiouvConfig(
                num_qubits = num_qubits, 
                with_coherent = with_coherent,
                with_linear_combination = with_linear_combination, 
                domain = domain,
                beta = beta,
                a = a,
                b = b,
                num_energy_bits = num_energy_bits,
                w0 = w0,
                t0 = t0,
                eta = eta,
                num_trotter_steps_per_t0 = num_trotter_steps_per_t0
        )

precomputed_data = precompute_data(config.domain, config)
time_oft_caches = OFTCaches(dim)

#* Hamiltonian
# hamiltonian_terms = [[X, X], [Y, Y], [Z, Z]]
# hamiltonian_coeffs = fill(1.0, length(hamiltonian_terms))
# hamiltonian = create_hamham(hamiltonian_terms, hamiltonian_coeffs, num_qubits)
hamiltonian = load_hamiltonian("heis", num_qubits)
hamiltonian = finalize_hamham(hamiltonian, beta)

#* Trotter
trotter = nothing
# trotter = create_trotter(hamiltonian, t0, num_trotter_steps_per_t0)
# trotter_error_T = compute_trotter_error(hamiltonian, trotter, 2^num_energy_bits * t0 / 2)
# gibbs_in_trotter = Hermitian(trotter.eigvecs' * gibbs_state(hamiltonian, beta) * trotter.eigvecs)

#* Jumps
jump_paulis = [[X], [Y], [Z]]
jump_sites = 1:num_qubits
num_of_jumps = length(jump_paulis) * length(jump_sites)
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

kraus_jumps = @time precompute_kraus_jumps(config.domain, jumps, hamiltonian, config, precomputed_data, time_oft_caches)
# verify_completeness(kraus_jumps)

# num_sparse_jumps = 0
# for kraus_jump in kraus_jumps
#         replace!(x -> abs(x) <= 1e-12 ? zero(x) : x, kraus_jump)
#         sparsity = count(iszero, kraus_jump) / length(kraus_jump)
#         if sparsity >= 0.9
#                 num_sparse_jumps += 1
#         end
# end
# num_sparse_jumps

# kraus_jumps = precompute_kraus_jumps(config.domain, jumps, trotter, config, precomputed_data, time_oft_caches)
R = @time precompute_R(kraus_jumps)

if with_coherent
        B = @time B_time(jumps, hamiltonian, precomputed_data.b_minus, precomputed_data.b_plus, t0, beta)
        # B = B_trotter(jumps, trotter, precomputed_data.b_minus, precomputed_data.b_plus, beta)
else
        B = zeros(ComplexF64, dim, dim)
end

fw = build_krausframework(B, kraus_jumps, R, delta)
# Base.summarysize(kraus_jumps)  # Size

psi0 = ones(ComplexF64, dim)
psi0 ./= norm(psi0)

total_time = 20.0
num_trajectories = 1
rho_combined = zeros(ComplexF64, dim, dim)

psi = @btime evolve_along_trajectory(psi0, fw, total_time)

# @time for trajectory in 1:num_trajectories
#         evolve_along_trajectory(psi0, fw, total_time)
#         BLAS.herk!('U', 'N', 1.0/num_trajectories, psi, 1.0, rho_combined)
# end
# println("done")
# LinearAlgebra.copytri!(rho_combined, 'U', true)

# norm(hamiltonian.gibbs - rho_combined)


# lindbladian_from_kraus = construct_gksl_lindbladian(B, kraus_jumps)
# result = run_liouvillian(jumps, config, hamiltonian; trotter=trotter)
# norm(lindbladian_from_kraus - result.data)
# norm(result.fixed_point - hamiltonian.gibbs)

# eigvals_near_zero, eigvecs_near_zero = eigen(lindbladian_from_kraus)
# sorted_permutation_eigen = sortperm(abs.(real.(eigvals_near_zero)))

# ss_index = sorted_permutation_eigen[1]   # Smallest
# gap_index = sorted_permutation_eigen[2]  # Second smallest
# spectral_gap = eigvals_near_zero[gap_index] # Spectral gap

# steady_state_vec = eigvecs_near_zero[:, ss_index]
# steady_state_dm = reshape(steady_state_vec, size(hamiltonian.data))
# steady_state_dm = (steady_state_dm + steady_state_dm') / 2
# steady_state_dm ./= tr(steady_state_dm) # Normalize

# norm(hamiltonian.gibbs - steady_state_dm)
# norm(gibbs_in_trotter - steady_state_dm)
