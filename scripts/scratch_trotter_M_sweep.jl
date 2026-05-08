#!/usr/bin/env julia
# Sweep `num_trotter_steps_per_t0` (M) for the TrotterDomain coherent
# Lindbladian and check whether KMS-DB rel_norm converges as M → ∞.
# Compares two Gibbs reference conventions:
#   1. gibbs of Trotter Hamiltonian H_T (diagonal in trotter.eigvecs).
#   2. gibbs of original H, expressed in trotter.eigvecs basis.
# If neither converges below ~1e-3 by M=200, there's a real bug in the
# TrotterDomain coherent path.

using QuantumFurnace
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "scratch_non_hermitian_kms_db.jl"))

ham = _load_n3_ham(5.0)

println("M_user sweep at β=5, n=3, Hermitian X single-site jump on site 1.")
println("--------------------------------------------------------------------------------")
for M in (10, 50, 100, 200, 500)
    cfg = Config(;
        sim = Lindbladian(),
        domain = TrotterDomain(),
        construction = KMS(),
        num_qubits = 3,
        with_linear_combination = true,
        beta = 5.0,
        sigma = 1.0/5.0,
        a = 5.0/30.0,
        s = 0.4,
        num_energy_bits = 14,
        w0 = 0.1,
        t0 = 2pi/(2^14 * 0.1),
        num_trotter_steps_per_t0 = M,
    )
    trotter = make_trotter_for_config(ham, cfg)
    jumps = build_hermitian_baseline(ham; basis=trotter.eigvecs)
    L = Matrix{ComplexF64}(construct_lindbladian(jumps, cfg, ham; trotter=trotter))

    # Convention 1: gibbs of H_T, diagonal in Trotter eigenbasis.
    gibbs_T = trotter_gibbs(trotter, 5.0)
    res_T = verify_detailed_balance(L, gibbs_T)

    # Convention 2: gibbs of H expressed in Trotter eigenbasis.
    # ham.gibbs is diagonal in H eigenbasis. Transform comp ← H eigen, then T eigen ← comp.
    gibbs_H_comp = ham.eigvecs * ham.gibbs * ham.eigvecs'
    gibbs_H_in_T = Hermitian(trotter.eigvecs' * gibbs_H_comp * trotter.eigvecs)
    res_H = verify_detailed_balance(L, gibbs_H_in_T)

    @printf("  M=%4d:  gibbs_T rel_norm = %.3e   gibbs_H_in_T rel_norm = %.3e\n",
            M, res_T.relative_norm, res_H.relative_norm)
end
