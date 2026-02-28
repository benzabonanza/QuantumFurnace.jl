"""
Shared test fixtures, tolerance constants, and factory functions for QuantumFurnace test suite.

Constants here are computed once at include time and reused across all test files.
"""

using QuantumFurnace
using LinearAlgebra
using BSON

# ---------------------------------------------------------------------------
# BSON legacy loader (test-only)
# ---------------------------------------------------------------------------
"""
    _load_test_hamiltonian(ham_path, beta) -> HamHam

Load a legacy BSON-serialized HamHam and reconstruct it with the new fully-initialized
struct definition. Uses BSON.parse to avoid deserialization failure from the changed
HamHam field types (bohr_freqs, bohr_dict, gibbs are no longer Union{..., Nothing}).
"""
function _load_test_hamiltonian(ham_path::String, beta::Float64)
    raw = open(ham_path) do io
        BSON.parse(io)
    end
    ham_raw = raw[:hamiltonian]
    fields = ham_raw[:data]

    # Legacy HamHam field order (14 fields):
    #   1:data, 2:bohr_freqs(nothing), 3:bohr_dict(nothing), 4:base_terms,
    #   5:base_coeffs, 6:disordering_term, 7:disordering_coeffs,
    #   8:eigvals, 9:eigvecs, 10:nu_min, 11:shift, 12:rescaling_factor,
    #   13:periodic, 14:gibbs
    cache = IdDict()
    init = QuantumFurnace

    data_matrix = BSON.raise_recursive(fields[1], cache, init)::Matrix{ComplexF64}
    base_terms = Vector{Vector{Matrix{ComplexF64}}}(BSON.raise_recursive(fields[4], cache, init))
    base_coeffs = BSON.raise_recursive(fields[5], cache, init)::Vector{Float64}
    disordering_term = let dt = BSON.raise_recursive(fields[6], cache, init)
        dt === nothing ? nothing : Vector{Matrix{ComplexF64}}(dt)
    end
    disordering_coeffs = let dc = BSON.raise_recursive(fields[7], cache, init)
        dc === nothing ? nothing : Vector{Float64}(dc)
    end
    eigvals_vec = BSON.raise_recursive(fields[8], cache, init)::Vector{Float64}
    eigvecs_mat = BSON.raise_recursive(fields[9], cache, init)::Matrix{ComplexF64}
    nu_min = Float64(fields[10])
    shift = Float64(fields[11])
    rescaling_factor = Float64(fields[12])
    periodic = Bool(fields[13])

    raw_nt = (
        matrix = data_matrix,
        terms = base_terms,
        base_coeffs = base_coeffs,
        disordering_term = disordering_term,
        disordering_coeffs = disordering_coeffs,
        eigvals = eigvals_vec,
        eigvecs = eigvecs_mat,
        nu_min = nu_min,
        shift = shift,
        rescaling_factor = rescaling_factor,
        periodic = periodic,
    )

    return HamHam(raw_nt, beta)
end

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
# All domain singletons
# ---------------------------------------------------------------------------
const ALL_DOMAINS = [EnergyDomain(), TimeDomain(), TrotterDomain(), BohrDomain()]

# ---------------------------------------------------------------------------
# Unified test system factory
# ---------------------------------------------------------------------------
"""
    make_test_system(; num_qubits=NUM_QUBITS, trotter=nothing) -> (; hamiltonian, jumps, gibbs)

Load the `num_qubits`-qubit disordered Heisenberg Hamiltonian at inverse temperature BETA
and create `3 * num_qubits` single-site Pauli jump operators (X, Y, Z on each site),
normalized by sqrt(3 * num_qubits).

Returns a named tuple with:
- `hamiltonian`: fully-initialized HamHam with bohr_dict and gibbs populated
- `jumps`: Vector{JumpOp} of jump operators
- `gibbs`: the Gibbs state matrix (Hermitian, trace 1)
"""
function make_test_system(; num_qubits::Int=NUM_QUBITS, trotter::Union{Nothing, TrottTrott}=nothing)
    # Load Hamiltonian directly using the source tree path
    # (load_hamiltonian uses Pkg.project().path which points to a temp dir during Pkg.test())
    source_root = dirname(@__DIR__)
    ham_path = joinpath(source_root, "hamiltonians", "heis_disordered_periodic_n$(num_qubits).bson")
    hamiltonian = _load_test_hamiltonian(ham_path, BETA)

    # Create jump operators: single-site Paulis (X, Y, Z) on each site
    jump_paulis = [[X], [Y], [Z]]
    jump_sites = 1:num_qubits
    num_of_jumps = length(jump_paulis) * length(jump_sites)
    jump_normalization = sqrt(num_of_jumps)

    # Select basis: trotter.eigvecs for TrotterDomain, hamiltonian.eigvecs otherwise
    basis_unitary = trotter !== nothing ? trotter.eigvecs : hamiltonian.eigvecs

    jumps = JumpOp[]
    for pauli in jump_paulis
        for site in jump_sites
            jump_op = Matrix(pad_term(pauli, num_qubits, site)) ./ jump_normalization

            jump_in_eigen = basis_unitary' * jump_op * basis_unitary
            orthogonal = (jump_op == transpose(jump_op))
            herm = (jump_op == jump_op')
            push!(jumps, JumpOp(jump_op, jump_in_eigen, orthogonal, herm))
        end
    end

    gibbs = hamiltonian.gibbs

    return (; hamiltonian, jumps, gibbs)
end

# Compute once at include time (4-qubit)
const TEST_SYSTEM = make_test_system()
const TEST_HAM = TEST_SYSTEM.hamiltonian
const TEST_JUMPS = TEST_SYSTEM.jumps
const TEST_GIBBS = TEST_SYSTEM.gibbs

# ---------------------------------------------------------------------------
# 3-qubit test system (N3_* globals)
# ---------------------------------------------------------------------------
const N3_SYSTEM = make_test_system(; num_qubits=3)
const N3_HAM = N3_SYSTEM.hamiltonian
const N3_JUMPS = N3_SYSTEM.jumps
const N3_GIBBS = N3_SYSTEM.gibbs
const N3_DIM = 2^3  # 8

# ---------------------------------------------------------------------------
# Trotter helpers
# ---------------------------------------------------------------------------
const TEST_TROTTER = TrottTrott(TEST_HAM, T0, NUM_TROTTER_STEPS_PER_T0)
const TEST_TROTTER_JUMPS = make_test_system(; trotter=TEST_TROTTER).jumps

const N3_TROTTER = TrottTrott(N3_HAM, T0, NUM_TROTTER_STEPS_PER_T0)
const N3_TROTTER_JUMPS = make_test_system(; num_qubits=3, trotter=N3_TROTTER).jumps

# ---------------------------------------------------------------------------
# Unified config factory
# ---------------------------------------------------------------------------
"""
    make_config(sim, domain; num_qubits=NUM_QUBITS, construction=KMS(), delta=TEST_DELTA, mixing_time=1.0) -> Config

Create a Config with locked test parameters.

For `Lindbladian()` sim: creates `Config{Lindbladian}` (no mixing_time/delta).
For `Thermalize()` sim: creates `Config{Thermalize}` with mixing_time and delta.
"""
function make_config(sim, domain;
    num_qubits::Int=NUM_QUBITS,
    construction=KMS(),
    delta::Float64=TEST_DELTA,
    mixing_time::Float64=1.0,
)
    therm_kw = sim isa Thermalize ? (; mixing_time=mixing_time, delta=delta) : (;)
    Config(;
        sim = sim,
        domain = domain,
        construction = construction,
        num_qubits = num_qubits,
        with_linear_combination = true,
        beta = BETA,
        sigma = SIGMA,
        a = BETA / 30.0,
        b = 0.4,
        num_energy_bits = NUM_ENERGY_BITS,
        w0 = W0,
        t0 = T0,
        num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
        therm_kw...,
    )
end
