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
sigma = 1 / beta

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

config = Config(
                sim = Lindbladian(),
                domain = domain,
                construction = with_coherent ? KMS() : GNS(),
                num_qubits = num_qubits,
                with_linear_combination = with_linear_combination,
                beta = beta,
                sigma = sigma,
                a = a,
                b = b,
                num_energy_bits = num_energy_bits,
                w0 = w0,
                t0 = t0,
                eta = eta,
                num_trotter_steps_per_t0 = num_trotter_steps_per_t0
        )

#* Hamiltonian
hamiltonian = load_hamiltonian("heis", num_qubits)
hamiltonian = finalize_hamham(hamiltonian, beta)

unique_freqs = keys(hamiltonian.bohr_dict)

alpha_A_nu1_old = zeros(ComplexF64, dim, dim)
alpha_A_nu1 = zeros(ComplexF64, dim, dim)

for nu_2 in unique_freqs
    @. alpha_A_nu1_old = create_alpha_old(hamiltonian.bohr_freqs, nu_2, config.beta, config.a, config.b)
    @. alpha_A_nu1 = create_alpha(hamiltonian.bohr_freqs, nu_2, config.beta, config.sigma, config.a, config.b)
    display(norm(alpha_A_nu1_old - alpha_A_nu1))
end
