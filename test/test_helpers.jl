"""
Shared test fixtures, tolerance constants, and factory functions for QuantumFurnace test suite.

Constants here are computed once at include time and reused across all test files.
"""

using QuantumFurnace
using LinearAlgebra
using BSON

# ---------------------------------------------------------------------------
# Physical parameters (LOCKED decisions)
# ---------------------------------------------------------------------------
const NUM_QUBITS = 4
const DIM = 2^NUM_QUBITS  # 16
const BETA = 10.0
const SIGMA = 1.0 / BETA  # 0.1

# ---------------------------------------------------------------------------
# Tolerance tiers (LOCKED decisions)
# ---------------------------------------------------------------------------
const TOL_EXACT = 1e-12          # machine precision identities
const TOL_QUADRATURE = 1e-6      # quadrature / discretization errors
TOL_DELTA(delta) = 5.0 * delta   # unraveling error, C = 5.0

# ---------------------------------------------------------------------------
# Test step size
# ---------------------------------------------------------------------------
const TEST_DELTA = 0.01

# ---------------------------------------------------------------------------
# Grid parameters (matching trajectory_test.jl conventions)
# ---------------------------------------------------------------------------
const NUM_ENERGY_BITS = 12
const W0 = 0.05
const T0 = 2pi / (2^NUM_ENERGY_BITS * W0)
const NUM_TROTTER_STEPS_PER_T0 = 10

# ---------------------------------------------------------------------------
# Test system: Hamiltonian, jump operators, Gibbs state
# ---------------------------------------------------------------------------
"""
    make_test_system() -> (; hamiltonian, jumps, gibbs)

Loads the 4-qubit disordered Heisenberg Hamiltonian, finalizes it at inverse
temperature BETA, and creates 12 single-site Pauli jump operators (X, Y, Z
on each of 4 sites), normalized by sqrt(3 * NUM_QUBITS).

Returns a named tuple with:
- `hamiltonian`: finalized HamHam with bohr_dict and gibbs populated
- `jumps`: Vector{JumpOp} of 12 jump operators
- `gibbs`: the Gibbs state matrix (Hermitian, trace 1)
"""
function make_test_system()
    # Load Hamiltonian directly using the source tree path
    # (load_hamiltonian uses Pkg.project().path which points to a temp dir during Pkg.test())
    source_root = dirname(@__DIR__)
    ham_path = joinpath(source_root, "hamiltonians", "heis_disordered_periodic_n$(NUM_QUBITS).bson")
    bson_data = BSON.load(ham_path)
    hamiltonian = bson_data[:hamiltonian]
    hamiltonian = finalize_hamham(hamiltonian, BETA)

    # Create jump operators: single-site Paulis (X, Y, Z) on each site
    jump_paulis = [[X], [Y], [Z]]
    jump_sites = 1:NUM_QUBITS
    num_of_jumps = length(jump_paulis) * length(jump_sites)
    jump_normalization = sqrt(num_of_jumps)

    jumps = JumpOp[]
    for pauli in jump_paulis
        for site in jump_sites
            jump_op = Matrix(pad_term(pauli, NUM_QUBITS, site)) ./ jump_normalization

            jump_in_eigen = hamiltonian.eigvecs' * jump_op * hamiltonian.eigvecs
            orthogonal = (jump_op == transpose(jump_op))
            herm = (jump_op == jump_op')
            push!(jumps, JumpOp(jump_op, jump_in_eigen, orthogonal, herm))
        end
    end

    gibbs = hamiltonian.gibbs

    return (; hamiltonian, jumps, gibbs)
end

# Compute once at include time
const TEST_SYSTEM = make_test_system()
const TEST_HAM = TEST_SYSTEM.hamiltonian
const TEST_JUMPS = TEST_SYSTEM.jumps
const TEST_GIBBS = TEST_SYSTEM.gibbs

# ---------------------------------------------------------------------------
# Small test system: 3-qubit Hamiltonian, jump operators, Gibbs state
# ---------------------------------------------------------------------------
"""
    make_small_test_system() -> (; hamiltonian, jumps, gibbs)

Loads the 3-qubit disordered Heisenberg Hamiltonian, finalizes it at inverse
temperature BETA, and creates 9 single-site Pauli jump operators (X, Y, Z
on each of 3 sites), normalized by sqrt(9).

Returns a named tuple with:
- `hamiltonian`: finalized HamHam with bohr_dict and gibbs populated
- `jumps`: Vector{JumpOp} of 9 jump operators
- `gibbs`: the Gibbs state matrix (Hermitian, trace 1)
"""
function make_small_test_system()
    small_num_qubits = 3
    source_root = dirname(@__DIR__)
    ham_path = joinpath(source_root, "hamiltonians", "heis_disordered_periodic_n$(small_num_qubits).bson")
    bson_data = BSON.load(ham_path)
    hamiltonian = bson_data[:hamiltonian]
    hamiltonian = finalize_hamham(hamiltonian, BETA)

    # Create jump operators: single-site Paulis (X, Y, Z) on each site
    jump_paulis = [[X], [Y], [Z]]
    jump_sites = 1:small_num_qubits
    num_of_jumps = length(jump_paulis) * length(jump_sites)
    jump_normalization = sqrt(num_of_jumps)

    jumps = JumpOp[]
    for pauli in jump_paulis
        for site in jump_sites
            jump_op = Matrix(pad_term(pauli, small_num_qubits, site)) ./ jump_normalization

            jump_in_eigen = hamiltonian.eigvecs' * jump_op * hamiltonian.eigvecs
            orthogonal = (jump_op == transpose(jump_op))
            herm = (jump_op == jump_op')
            push!(jumps, JumpOp(jump_op, jump_in_eigen, orthogonal, herm))
        end
    end

    gibbs = hamiltonian.gibbs

    return (; hamiltonian, jumps, gibbs)
end

const SMALL_SYSTEM = make_small_test_system()
const SMALL_HAM = SMALL_SYSTEM.hamiltonian
const SMALL_JUMPS = SMALL_SYSTEM.jumps
const SMALL_GIBBS = SMALL_SYSTEM.gibbs
const SMALL_DIM = 2^3  # 8

# ---------------------------------------------------------------------------
# Trotter helper
# ---------------------------------------------------------------------------
"""
    make_test_trotter() -> TrottTrott

Create a TrottTrott object for TrotterDomain tests using the shared test Hamiltonian.
"""
function make_test_trotter()
    TrottTrott(TEST_HAM, T0, NUM_TROTTER_STEPS_PER_T0)
end

const TEST_TROTTER = make_test_trotter()

# ---------------------------------------------------------------------------
# Factory functions for configs
# ---------------------------------------------------------------------------
"""
    make_liouv_config(domain; with_coherent=true) -> LiouvConfig

Create a LiouvConfig with locked test parameters.
"""
function make_liouv_config(domain; with_coherent::Bool=true)
    LiouvConfig(
        num_qubits = NUM_QUBITS,
        with_coherent = with_coherent,
        with_linear_combination = true,
        domain = domain,
        beta = BETA,
        sigma = SIGMA,
        a = BETA / 30.0,
        b = 0.4,
        num_energy_bits = NUM_ENERGY_BITS,
        w0 = W0,
        t0 = T0,
        num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
    )
end

"""
    make_thermalize_config(domain; with_coherent=true, delta=TEST_DELTA, mixing_time=1.0) -> ThermalizeConfig

Create a ThermalizeConfig with locked test parameters.
"""
function make_thermalize_config(domain;
    with_coherent::Bool=true,
    delta::Float64=TEST_DELTA,
    mixing_time::Float64=1.0,
)
    ThermalizeConfig(
        num_qubits = NUM_QUBITS,
        with_coherent = with_coherent,
        with_linear_combination = true,
        domain = domain,
        beta = BETA,
        sigma = SIGMA,
        a = BETA / 30.0,
        b = 0.4,
        num_energy_bits = NUM_ENERGY_BITS,
        w0 = W0,
        t0 = T0,
        num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
        mixing_time = mixing_time,
        delta = delta,
    )
end

# ---------------------------------------------------------------------------
# Small Trotter helper (3-qubit)
# ---------------------------------------------------------------------------
function make_small_test_trotter()
    TrottTrott(SMALL_HAM, T0, NUM_TROTTER_STEPS_PER_T0)
end

const SMALL_TROTTER = make_small_test_trotter()

# ---------------------------------------------------------------------------
# Small config factories (3-qubit)
# ---------------------------------------------------------------------------
"""
    make_small_thermalize_config(domain; with_coherent=false, delta=TEST_DELTA, mixing_time=1.0) -> ThermalizeConfig

Create a ThermalizeConfig for the 3-qubit SMALL system.
"""
function make_small_thermalize_config(domain;
    with_coherent::Bool=false,
    delta::Float64=TEST_DELTA,
    mixing_time::Float64=1.0,
)
    ThermalizeConfig(
        num_qubits = 3,
        with_coherent = with_coherent,
        with_linear_combination = true,
        domain = domain,
        beta = BETA,
        sigma = SIGMA,
        a = BETA / 30.0,
        b = 0.4,
        num_energy_bits = NUM_ENERGY_BITS,
        w0 = W0,
        t0 = T0,
        num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
        mixing_time = mixing_time,
        delta = delta,
    )
end

"""
    make_small_liouv_config(domain; with_coherent=false) -> LiouvConfig

Create a LiouvConfig for the 3-qubit SMALL system.
"""
function make_small_liouv_config(domain; with_coherent::Bool=false)
    LiouvConfig(
        num_qubits = 3,
        with_coherent = with_coherent,
        with_linear_combination = true,
        domain = domain,
        beta = BETA,
        sigma = SIGMA,
        a = BETA / 30.0,
        b = 0.4,
        num_energy_bits = NUM_ENERGY_BITS,
        w0 = W0,
        t0 = T0,
        num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
    )
end
