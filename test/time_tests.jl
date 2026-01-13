using Revise

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
domain = TimeDomain()
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

#* Jumps
jump_paulis = [[X], [Y], [Z]]
jump_sites = 1:num_qubits
num_of_jumps = length(jump_paulis) * length(jump_sites)
jump_normalization = sqrt(num_of_jumps)
jumps::Vector{JumpOp} = []
for pauli in jump_paulis
        for site in jump_sites
                jump_op = Matrix(pad_term(pauli, num_qubits, site)) / jump_normalization
                jump_op_in_eigenbasis = hamiltonian.eigvecs' * jump_op * hamiltonian.eigvecs
                orthogonal = (jump_op == transpose(jump_op))
                jump = JumpOp(jump_op,
                        jump_op_in_eigenbasis,
                        orthogonal)
                push!(jumps, jump)
        end
end

dim = 2^num_qubits
A = jumps[2]
w = -3 * w0
A_oft = oft(A, w, hamiltonian, beta) * sqrt(beta / sqrt(2 * pi))
A_oft_time = time_oft(A, w, hamiltonian, precomputed_data.oft_time_labels, beta) * t0 * sqrt((sqrt(2 / pi)/beta) / (2 * pi))
norm(A_oft - A_oft_time)

A_oft_time_fast = zeros(ComplexF64, dim, dim)
time_oft_caches = OFTCaches(dim)
time_oft_fast!(A_oft_time_fast, time_oft_caches, A, w, hamiltonian, precomputed_data.oft_time_labels, beta)
A_oft_time_fast *= (t0 * sqrt((sqrt(2 / pi)/beta) / (2 * pi)))
norm(A_oft - A_oft_time_fast)