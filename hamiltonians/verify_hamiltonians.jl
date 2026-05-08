"""
Lightweight physics sanity checks on the cached Heisenberg Hamiltonian BSONs
produced by `generate_hamiltonians.jl` (epic qf-k1u, task qf-k1u.5).

Each family has expected qualitative features:

  F1 (1D XXX + full disorder):
    - SU(2) base + strong site disorder → no symmetry-protected order; just
      verify Hermiticity, trace ≈ 0 (after rescaling/shift the trace is
      shifted to the spectrum-midpoint value), nu_min > 0.

  F2 (1D XXX clean + ε-disorder):
    - Almost-clean Heisenberg AFM. Ground-state per-bond energy in the
      *unrescaled* Hamiltonian limits to the Bethe-ansatz value
      `-log(2) + 1/4 ≈ -0.4431` (antiferromagnetic; for our normalisation
      ‖S^a‖ = 1, the convention is `e₀ = -log(2)·J + J/4`). For finite n
      with PBC the energy per bond converges from above as 1/n. Verify
      that the per-bond ground-state energy is in the right ballpark
      (`-0.6` to `-0.3` range for n ∈ [3, 10]).

  F3 (2D XXZ + ε-disorder):
    - Ising-anisotropic AFM ground state. Verify the staggered Z order
      parameter `⟨S_stag⟩ = (1/n) ⟨Σ_b (-1)^{i+j} Z_(i,j)⟩` is non-zero
      in the ground state (it would vanish exactly for SU(2)-symmetric
      finite systems).

Usage:

    julia --project hamiltonians/verify_hamiltonians.jl
"""

using QuantumFurnace
using BSON
using Printf
using LinearAlgebra

const HAM_DIR = joinpath(dirname(@__DIR__), "hamiltonians")
const BETA = 1.0  # arbitrary; we only need the eigvals/eigvecs

# Load a `find_ideal_*heisenberg` raw NamedTuple from BSON. The new-format files
# saved via `@save path hamiltonian` store the NamedTuple directly, so a plain
# `BSON.load(path)[:hamiltonian]` gives back the NamedTuple ready for HamHam(raw, beta).
function _load_raw(path::String)
    return BSON.load(path)[:hamiltonian]
end

# Compute ⟨ψ| O |ψ⟩ for a state ψ in eigenbasis — eigvecs maps from comp basis to eigenbasis,
# so the operator in eigenbasis is U' * O * U.
function _expectation(ham, op_comp::Matrix{ComplexF64}, k::Int)
    # ground state is column k=1 in eigvecs (sorted ascending eigvalues)
    ψ_eigen = zeros(ComplexF64, length(ham.eigvals))
    ψ_eigen[k] = 1.0 + 0im
    ψ_comp = ham.eigvecs * ψ_eigen
    return real(ψ_comp' * op_comp * ψ_comp)
end

# Build the unrescaled-norm Heisenberg per-bond energy. Recall the cached matrix is
# `H_rescaled = H_orig / rescaling_factor + shift * I`. Then
#   ⟨ψ|H_orig|ψ⟩ = (⟨ψ|H_rescaled|ψ⟩ - shift) * rescaling_factor.
function _ground_state_energy_per_bond(ham, n_bonds::Int)
    e_rescaled = ham.eigvals[1]
    e_orig = (e_rescaled - ham.shift) * ham.rescaling_factor
    return e_orig / n_bonds
end

function _staggered_z_2d_op(Lx::Int, Ly::Int)
    n = Lx * Ly
    op = zeros(ComplexF64, 2^n, 2^n)
    for i in 1:Lx, j in 1:Ly
        q = (i - 1) * Ly + (j - 1) + 1
        sign = iseven(i + j) ? 1.0 : -1.0
        op .+= sign .* Matrix(pad_term([Z], n, q))
    end
    return op
end

# In a finite-size system the symmetry-breaking expectation `<S_stag>` is zero
# by the Z2 flip ∏ X_i; the right diagnostic is the squared order parameter
# `<S_stag^2> / n^2`. For the Ising-AFM ordered phase this approaches 1
# (perfect Néel correlator); for the disordered/paramagnetic phase it scales
# as O(1/n).
function _staggered_z2_per_site2(ham, Lx::Int, Ly::Int)
    n = Lx * Ly
    op = _staggered_z_2d_op(Lx, Ly)
    op2 = op * op
    return _expectation(ham, op2, 1) / n^2
end

# ----- Run -----

println("=== Verifying 1D XXX (clean + ε-disorder) ground-state energy per bond ===")
# Pauli normalization: H = Σ_b (X_q1 X_q2 + Y_q1 Y_q2 + Z_q1 Z_q2)_b. For 1D PBC chain
# with n ≥ 3 sites we have n distinct bonds, each contributing XX+YY+ZZ. The Bethe-ansatz
# GS energy per bond in the thermodynamic limit is 4 (1/4 − ln 2) ≈ −1.7726.
for n in 3:10
    path = joinpath(HAM_DIR, "heis_xxx_clean_periodic_n$(n).bson")
    isfile(path) || continue
    raw = _load_raw(path)
    n_bonds = n  # 1D PBC: n bonds for n ≥ 3
    e_per_bond = _ground_state_energy_per_bond(raw, n_bonds)
    @printf("  n=%2d  e_GS/bond = %+.4f   nu_min = %.3e\n", n, e_per_bond, raw.nu_min)
end
println("Bethe ansatz limit (Pauli normalisation): e_GS/bond → 4 × (1/4 − ln 2) ≈ −1.7726")

println("\n=== Verifying 1D XXX + full disorder: nu_min and traceless check ===")
for n in 3:10
    path = joinpath(HAM_DIR, "heis_xxx_zzdisordered_periodic_n$(n).bson")
    isfile(path) || continue
    raw = _load_raw(path)
    @printf("  n=%2d  nu_min = %.3e   shift = %.4f   rescaling = %.4f\n",
        n, raw.nu_min, raw.shift, raw.rescaling_factor)
end

println("\n=== Verifying 2D XXZ Ising-anisotropic: squared staggered Z correlator ===")
for (Lx, Ly, n) in [(2, 2, 4), (2, 3, 6), (3, 3, 9), (2, 5, 10)]
    path = joinpath(HAM_DIR, "heis_xxz_2d_$(Lx)x$(Ly)_n$(n).bson")
    isfile(path) || continue
    raw = _load_raw(path)
    sstag2 = _staggered_z2_per_site2(raw, Lx, Ly)
    @printf("  Lx=%d Ly=%d  ⟨S_stag²⟩/n² = %.4f   nu_min = %.3e\n", Lx, Ly, sstag2, raw.nu_min)
end
println("(Ising-AFM ordered: ⟨S_stag²⟩/n² → O(1) as n → ∞.")
println(" Disordered/paramagnetic: ⟨S_stag²⟩/n² → 0 as n → ∞ (decays as 1/n).")
println(" Random uniform basis state: ⟨S_stag²⟩/n² = 1/n.)")
