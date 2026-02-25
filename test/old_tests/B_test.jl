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
sigma = 0.1 / beta
w_gamma = 1 / beta
sigma_gamma = sqrt(2 * w_gamma / beta - sigma^2)

# sigma = 1 / beta
# w_gamma = 1 / beta
# sigma_gamma = sqrt(2 * w_gamma / beta - sigma^2)

# Smooth Metro
a = 1 / 30. # a = beta / 50.
b = 0.4  # b = 0.5
eta = 0.0  # eta = 0.2

# Kinky Metro 
# a = 0.0
# b = 0.0
# eta = 0.002

with_linear_combination = true
domain = TimeDomain()
num_energy_bits = 18  # 11
w0 = 0.005
max_E = w0 * 2^num_energy_bits / 2
t0 = 2pi / (2^num_energy_bits * w0)  # Max time evolution pi / w0
num_trotter_steps_per_t0 = 10

config = Config(
                sim = Lindbladian(),
                domain = domain,
                construction = KMS(),
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

jump = jumps[1]
B_bohr = B_bohr(hamiltonian, jump, config) 
rmul!(B_bohr, precomputed_data.gamma_norm_factor)

B_t = B_time(jump, hamiltonian, precomputed_data.b_minus, precomputed_data.b_plus, config.t0, config.beta, config.sigma)
rmul!(B_t, precomputed_data.gamma_norm_factor)

norm(B_bohr - B_t)



