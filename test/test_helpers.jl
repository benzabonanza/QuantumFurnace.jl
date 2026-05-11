"""
Shared test fixtures, tolerance constants, and factory functions for QuantumFurnace test suite.

Constants here are computed once at include time and reused across all test files.
"""

using QuantumFurnace
using LinearAlgebra
using Test: @test
using BSON

# ---------------------------------------------------------------------------
# BSON legacy loader (test-only) — wraps the package-internal loader so that
# fixtures can pass an explicit path (Pkg.test() runs from a temp project dir
# where load_hamiltonian's Pkg.project() lookup misses the source-tree files).
# ---------------------------------------------------------------------------
const _load_test_hamiltonian = QuantumFurnace._load_hamiltonian_bson

# ---------------------------------------------------------------------------
# Shared physics fixtures (test-only)
# ---------------------------------------------------------------------------
"""
    make_dll_n3_system(beta) -> (; ham, jumps, gibbs)

Build the n=3 disordered Heisenberg fixture at a specified `beta` for DLL
tests (test_dll_coherent.jl, test_dll_dissipator.jl, test_dll_kms_db.jl all
share this fixture). Differs from `make_test_system` only in that `beta`
is a per-call argument rather than the fixed `BETA` constant — DLL tests
need to sweep β ∈ {1, 5, 10}.

Returns the same shape as `make_test_system`.
"""
function make_dll_n3_system(beta::Real)
    source_root = dirname(@__DIR__)
    ham_path = joinpath(source_root, "hamiltonians", "heis_disordered_periodic_n3.bson")
    ham = _load_test_hamiltonian(ham_path, Float64(beta))
    jump_paulis = [[X], [Y], [Z]]
    num_jumps = length(jump_paulis) * 3
    jump_norm = sqrt(num_jumps)
    jumps = JumpOp[]
    for pauli in jump_paulis
        for site in 1:3
            op = Matrix(pad_term(pauli, 3, site)) ./ jump_norm
            op_eb = ham.eigvecs' * op * ham.eigvecs
            push!(jumps, JumpOp(op, op_eb, op == transpose(op), op == op'))
        end
    end
    return (; ham, jumps, gibbs = ham.gibbs)
end

"""
    assert_kms_skew_symmetric(α, ν_grid, β; atol=1e-12)

Structural witness of KMS detailed balance at the Kossakowski level
(Ding–Li–Lin 2024, Eq. 4.7):

    α(ν, ν') = α(-ν', -ν) · exp(-β(ν+ν')/2)

`ν_grid` must be symmetric around 0 (so each ν has its negation in the
grid) and `α` must be a `K × K` matrix with `K = length(ν_grid)`. Used
as a per-filter assertion across the DLL Kossakowski tests.
"""
function assert_kms_skew_symmetric(α::AbstractMatrix, ν_grid::AbstractVector,
                                   β::Real; atol::Real = 1e-12)
    K = length(ν_grid)
    @assert size(α) == (K, K)
    neg_idx = [findfirst(==(-ν_grid[k]), ν_grid) for k in 1:K]
    @assert all(!isnothing, neg_idx) "ν_grid must be symmetric around 0"
    for q in 1:K, p in 1:K
        lhs = α[p, q]
        rhs = α[neg_idx[q], neg_idx[p]] * exp(-β * (ν_grid[p] + ν_grid[q]) / 2)
        @test abs(lhs - rhs) <= atol
    end
end

# ---------------------------------------------------------------------------
# Physical parameters (LOCKED decisions)
# ---------------------------------------------------------------------------
# `BETA` is the algorithm-side inverse temperature β_alg (against the
# rescaled spectrum stored in `ham.eigvals`). Equal to `cfg.beta` in every
# Config constructed below. `SIGMA = 1/BETA` lives on the same scale. The
# qf-6vr refactor (Phase qf-bphys) keeps these semantics unchanged — see
# `Config.beta` docstring in `src/structs.jl`.
const NUM_QUBITS = 4
const DIM = 2^NUM_QUBITS  # 16
const BETA = 10.0
const BETA_ALG = BETA      # qf-6vr explicit alias for self-documenting tests
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
# qf-6vr: physical inverse temperature for the n=4 fixture; satisfies
# `BETA == BETA_PHYS · TEST_HAM.rescaling_factor` (= BETA_ALG by construction).
const BETA_PHYS = BETA / TEST_HAM.rescaling_factor

# ---------------------------------------------------------------------------
# 3-qubit test system (N3_* globals)
# ---------------------------------------------------------------------------
const N3_SYSTEM = make_test_system(; num_qubits=3)
const N3_HAM = N3_SYSTEM.hamiltonian
const N3_JUMPS = N3_SYSTEM.jumps
const N3_GIBBS = N3_SYSTEM.gibbs
const N3_DIM = 2^3  # 8
# qf-6vr: physical inverse temperature for the n=3 fixture; satisfies
# `BETA == N3_BETA_PHYS · N3_HAM.rescaling_factor`.
const N3_BETA_PHYS = BETA / N3_HAM.rescaling_factor

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
        s = 0.4,
        num_energy_bits = NUM_ENERGY_BITS,
        w0 = W0,
        t0 = T0,
        num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
        therm_kw...,
    )
end
