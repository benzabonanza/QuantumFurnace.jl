using Revise

if !isdefined(Main, :QuantumFurnace)
    includet("../src/QuantumFurnace.jl")
end

using .QuantumFurnace
using LinearAlgebra, Random

num_qubits = 3

hamiltonian_terms = [[X, X], [Y, Y], [Z, Z]]
hamiltonian_coeffs = fill(1.0, length(hamiltonian_terms))
hamiltonian = create_hamham(hamiltonian_terms, hamiltonian_coeffs, num_qubits)

jump_rate = 0.5
jump_op = (X + 1im * Y) / 2
padded_jump = pad_term([jump_op], num_qubits, 1)
L = sqrt(jump_rate) * padded_jump
L_in_eigenbasis = hamiltonian.eigvecs' * L * hamiltonian.eigvecs
orthogonal = (L == transpose(L))
jump = JumpOp(L, L_in_eigenbasis, orthogonal)

[jump.in_eigenbasis]

delta = 0.01
fw = build_krausframework(hamiltonian.data, [jump.in_eigenbasis], delta)