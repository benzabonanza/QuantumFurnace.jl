using Revise

if !isdefined(Main, :QuantumFurnace)
    includet("../src/QuantumFurnace.jl")
end

using .QuantumFurnace
using LinearAlgebra, Random

#* Config
num_qubits = 8
dim = 2^num_qubits
beta = 10.
sigma = 0.1 / beta
w_gamma = 1 / beta
sigma_gamma = sqrt(2 * w_gamma / beta - sigma^2)

sigma = 1 / beta
w_gamma = 1 / beta
sigma_gamma = sqrt(2 * w_gamma / beta - sigma^2)

# Smooth Metro
a = 1 / 30. # a = beta / 50.
b = 0.4  # b = 0.5
eta = 0.0  # eta = 0.2

# Kinky Metro 
# a = 0.0
# b = 0.0
# eta = 0.002

with_coherent = true
with_linear_combination = true
domain = TimeDomain()
num_energy_bits = 16  # 11
w0 = 0.0005
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
                sigma = sigma,
                gaussian_parameters = (w_gamma, sigma_gamma),
                a = a,
                b = b,
                num_energy_bits = num_energy_bits,
                w0 = w0,
                t0 = t0,
                eta = eta,
                num_trotter_steps_per_t0 = num_trotter_steps_per_t0
        )

precomputed_data = precompute_data(config)
time_oft_caches = OFTCaches(dim)

#* Hamiltonian
# hamiltonian_terms = [[X, X], [Y, Y], [Z, Z]]
# hamiltonian_coeffs = fill(1.0, length(hamiltonian_terms))
# hamiltonian = HamHam(hamiltonian_terms, hamiltonian_coeffs, num_qubits)
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

energy_oft_prefactor = 1 / sqrt(config.sigma * sqrt(2 * pi))
time_oft_prefactor = config.t0 * sqrt(config.sigma * sqrt(2 / pi) / (2 * pi))

dim = 2^num_qubits
A = jumps[5]
w = -3 * w0
A_oft = Matrix{ComplexF64}(undef, dim, dim)
oft!(A_oft, A, w, hamiltonian, sigma)
A_oft .*= energy_oft_prefactor
A_oft_time = Matrix{ComplexF64}(undef, dim, dim)
@time begin
        time_oft_caches = OFTCaches(dim)
        time_oft!(A_oft_time, time_oft_caches, A, w, hamiltonian, precomputed_data.oft_time_labels, sigma)
end
A_oft_time *= time_oft_prefactor
norm(A_oft - A_oft_time)

A_nufft = Matrix{ComplexF64}(undef, dim, dim)
@time begin
        nufft_caches = NUFFTCaches(hamiltonian, precomputed_data.oft_time_labels, sigma)  #!
        nufft_prefactor_matrix!(nufft_caches, w, precomputed_data.oft_time_labels)  #!
        oft_nufft!(A_nufft, A, w, hamiltonian, precomputed_data.oft_time_labels, nufft_caches)  #!
end
A_nufft *= time_oft_prefactor

norm(A_oft - A_nufft)
norm(A_oft_time - A_nufft)