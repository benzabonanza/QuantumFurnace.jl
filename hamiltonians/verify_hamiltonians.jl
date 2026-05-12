"""
Lightweight physics sanity checks on the cached Heisenberg Hamiltonian BSONs
produced by `generate_hamiltonians.jl` (qf-2kd).

Expected qualitative features for the 1D XXX + [[Z], [Z, Z]] typical family:

  - Strong site + bond disorder lifts SU(2) and ∏ X_i symmetries; the spectrum
    has no degeneracies, `nu_min > 0`.
  - The W₂-typical sample is the L²-spectral-median of `BATCH_SIZE` random
    disorder realisations, not the extremum, so the ground-state per-bond
    energy is closer to a typical disordered chain than to the clean XXX
    Bethe-ansatz value.
  - Per-bond ground-state energy (unrescaled Hamiltonian) is bounded; the
    exact clean XXX Bethe-ansatz limit `4(1/4 − ln 2) ≈ −1.7726` (Pauli
    normalisation) is NOT recovered because strong [[Z],[Z,Z]] disorder
    shifts the local fields. We just verify the magnitude is in a sane
    range (no `Inf` / `NaN`, |e_GS/bond| ≤ a few units).

Usage:

    julia --project hamiltonians/verify_hamiltonians.jl
"""

using QuantumFurnace
using BSON
using Printf
using LinearAlgebra

const HAM_DIR = joinpath(dirname(@__DIR__), "hamiltonians")

function _load_raw(path::String)
    return BSON.load(path)[:hamiltonian]
end

# Per-bond unrescaled-norm ground-state energy. The cached matrix is
# `H_rescaled = H_orig / rescaling_factor + shift * I`, so
#   ⟨ψ|H_orig|ψ⟩ = (⟨ψ|H_rescaled|ψ⟩ - shift) * rescaling_factor.
function _ground_state_energy_per_bond(ham, n_bonds::Int)
    e_rescaled = ham.eigvals[1]
    e_orig = (e_rescaled - ham.shift) * ham.rescaling_factor
    return e_orig / n_bonds
end

println("=== 1D XXX + [[Z], [Z, Z]] typical: spectrum + GS energy diagnostics ===")
for n in 3:9
    path = joinpath(HAM_DIR, "heis_xxx_zzdisordered_periodic_n$(n).bson")
    isfile(path) || (println("  n=$n: file missing, skipped"); continue)
    raw = _load_raw(path)
    n_bonds = n  # 1D PBC: n bonds for n ≥ 3
    e_per_bond = _ground_state_energy_per_bond(raw, n_bonds)
    typicality = hasproperty(raw, :typicality_distance) ?
                 (@sprintf("  typicality=%.3e", raw.typicality_distance)) : ""
    @printf("  n=%2d  e_GS/bond = %+.4f  nu_min = %.3e  rescale = %.3f%s\n",
        n, e_per_bond, raw.nu_min, raw.rescaling_factor, typicality)
end
println("(typical W₂-spectral-median selection; not equal to clean XXX Bethe-ansatz limit −1.7726.)")
