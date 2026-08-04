"""Single-qubit Pauli `X` matrix."""
const X = ComplexF64[0 1; 1 0]

"""Single-qubit Pauli `Y` matrix."""
const Y = ComplexF64[0 -im; im 0]

"""Single-qubit Pauli `Z` matrix."""
const Z = ComplexF64[1 0; 0 -1]

"""Single-qubit Hadamard matrix."""
const Had = ComplexF64[1 1; 1 -1] / sqrt(2)
