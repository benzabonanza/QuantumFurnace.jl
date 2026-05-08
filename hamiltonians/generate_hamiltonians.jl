"""
Generate cached Heisenberg Hamiltonian BSON files for the Phase 48 thesis-numerics
chapter (epic qf-k1u, prerequisite issue qf-k1u.5).

Three families are produced:

  (F1) 1D XXX with full disorder ``[[Z], [Z, Z]]`` (strength 1.0). Files:
       `heis_xxx_zzdisordered_periodic_n{3..10}.bson`. Replaces the legacy
       `heis_disordered_periodic_n*.bson` which used [[Z]] only and exhibited
       bipartite-pairing Bohr collisions for even n.

  (F2) 1D XXX with weak ε-disorder ``[[Z]]`` (strength 1e-2) — "clean
       Heisenberg with degeneracy-lifting noise". Files:
       `heis_xxx_clean_periodic_n{3..10}.bson`.

  (F3) 2D XXZ on a square lattice with periodic BCs in both directions, with
       Ising-like anisotropy ``[J_x, J_y, J_z] = [1, 1, 1.5]``. The bulk model
       has a finite-temperature Néel phase transition into a Z₂-symmetry-broken
       Ising-AFM ordered phase — cleaner physics than the SU(2)-symmetric XXX
       case, which has no finite-T transition by Mermin–Wagner. Weak ε-disorder
       ``[[X], [Z]]`` (strength 1e-2) breaks both Sz conservation and lattice
       symmetries, giving a finite minimum Bohr gap. Files:
       `heis_xxz_2d_{Lx}x{Ly}_n{n}.bson` for `(Lx, Ly, n) ∈
       {(2,2,4), (2,3,6), (3,3,9), (2,5,10)}`.

Usage:

    julia --project hamiltonians/generate_hamiltonians.jl                # laptop, n ≤ 10
    LAPTOP=false julia --project hamiltonians/generate_hamiltonians.jl   # cluster, allow n ≥ 11

The script prints the achieved `nu_min` and warns when below the soft target
1e-4. The simulator does not require `nu_min ≥ 1e-4`; smaller values just need
finer simulator parameters at consumption time.
"""

using QuantumFurnace
using Random
using BSON
using Printf

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

const NU_MIN_TARGET = 1e-4
const SEED = 42                    # fixed seed for reproducibility
const LAPTOP = get(ENV, "LAPTOP", "true") == "true"
const LAPTOP_MAX_N = 10

# Approximate per-qubit BSON file size (bytes), measured empirically:
# n=3 → 6.6 KB, n=8 → 3.1 MB, n=9 → 13 MB, n=10 → 49 MB → ~ 4× per qubit.
function _projected_bson_bytes(n::Int)
    return Int(round(6.6e3 * 4.0^(n - 3)))
end

const HAM_DIR = joinpath(dirname(@__DIR__), "hamiltonians")
isdir(HAM_DIR) || mkpath(HAM_DIR)

# -----------------------------------------------------------------------------
# Common helpers
# -----------------------------------------------------------------------------

function _save_and_report(path::String, hamiltonian::NamedTuple, label::String)
    # Use BSON.bson(path, Dict) instead of @save: the macro form has a known
    # truncation issue on BSON.jl 0.3.9 for large NamedTuples (≳ 32 MB on
    # virtio-fs in this Docker sandbox). The dict form serialises correctly.
    BSON.bson(path, Dict(:hamiltonian => hamiltonian))
    file_size_mb = filesize(path) / 1e6
    is_warning = hamiltonian.nu_min < NU_MIN_TARGET
    flag = is_warning ? "  ⚠️ below 1e-4 soft target" : ""
    @printf("  saved %s  (%.2f MB)  nu_min=%.3e%s\n", basename(path),
        file_size_mb, hamiltonian.nu_min, flag)
end

function _abort_if_too_big(n::Int, family::String)
    if LAPTOP && n > LAPTOP_MAX_N
        projected_mb = _projected_bson_bytes(n) / 1e6
        @printf("\n[%s] n=%d would produce a ≈ %.0f MB BSON file; ", family, n, projected_mb)
        println("skipping on laptop. Set LAPTOP=false to override.")
        return true
    end
    return false
end

# -----------------------------------------------------------------------------
# Family 1: 1D XXX with [[Z], [Z, Z]] disorder
# -----------------------------------------------------------------------------

function generate_family_1d_zz_disordered(; n_range = 3:10)
    println("\n=== Family 1: 1D XXX + [[Z], [Z, Z]] disorder ===")
    coeffs = [1.0, 1.0, 1.0]
    dis_terms = [[Z], [Z, Z]]
    for n in n_range
        _abort_if_too_big(n, "F1") && continue
        Random.seed!(SEED + 1000 + n)
        @printf("  n=%2d (dim=%4d): ", n, 2^n)
        raw = find_ideal_heisenberg(n, coeffs;
            batch_size=200, periodic=true,
            disordering_terms=dis_terms, disorder_strength=1.0)
        path = joinpath(HAM_DIR, "heis_xxx_zzdisordered_periodic_n$(n).bson")
        _save_and_report(path, raw, "F1")
    end
end

# -----------------------------------------------------------------------------
# Family 2: 1D XXX clean + ε-disorder [[Z]]
# -----------------------------------------------------------------------------

function generate_family_1d_clean(; n_range = 3:10)
    println("\n=== Family 2: 1D XXX clean + ε-disorder [[Z]] (strength 1e-2) ===")
    coeffs = [1.0, 1.0, 1.0]
    dis_terms = [[Z]]
    for n in n_range
        _abort_if_too_big(n, "F2") && continue
        Random.seed!(SEED + 2000 + n)
        @printf("  n=%2d (dim=%4d): ", n, 2^n)
        raw = find_ideal_heisenberg(n, coeffs;
            batch_size=500, periodic=true,
            disordering_terms=dis_terms, disorder_strength=1e-2)
        path = joinpath(HAM_DIR, "heis_xxx_clean_periodic_n$(n).bson")
        _save_and_report(path, raw, "F2")
    end
end

# -----------------------------------------------------------------------------
# Family 3: 2D XXZ Ising-anisotropic + ε-disorder [[X], [Z]]
# -----------------------------------------------------------------------------

# 2D lattices to sweep, paired with their qubit count.
const LATTICES_2D = [(2, 2, 4), (2, 3, 6), (3, 3, 9), (2, 5, 10)]
# Ising-like anisotropy: J_z slightly larger than J_xy. The bulk model has a
# finite-T Néel transition; the small finite-size systems retain that physics
# in their low-T Gibbs state structure (long-range AFM correlations).
const COEFFS_2D = [1.0, 1.0, 1.5]
# Mixed transverse + longitudinal field disorder: random ε X-field + ε Z-field
# breaks Sz conservation, lattice translation, point-group symmetries, and
# the spin-flip symmetry ∏ X_i. Without the X part, Sz remains a good quantum
# number and cross-Sz accidental gaps stay astronomically small.
const DIS_TERMS_2D = [[X], [Z]]

function generate_family_2d(; lattices = LATTICES_2D)
    println("\n=== Family 3: 2D XXZ (Jz=1.5) + ε-disorder [[X], [Z]] (strength 1e-2) ===")
    for (Lx, Ly, n) in lattices
        _abort_if_too_big(n, "F3") && continue
        Random.seed!(SEED + 3000 + n)
        @printf("  Lx=%d Ly=%d n=%2d (dim=%4d): ", Lx, Ly, n, 2^n)
        raw = find_ideal_2d_heisenberg(Lx, Ly, COEFFS_2D;
            batch_size=500, periodic_x=true, periodic_y=true,
            disordering_terms=DIS_TERMS_2D, disorder_strength=1e-2)
        path = joinpath(HAM_DIR, "heis_xxz_2d_$(Lx)x$(Ly)_n$(n).bson")
        _save_and_report(path, raw, "F3")
    end
end

# -----------------------------------------------------------------------------
# Top level
# -----------------------------------------------------------------------------

println("Hamiltonian generation for Phase 48 thesis numerics (epic qf-k1u, task qf-k1u.5)")
println("Output dir: ", HAM_DIR)
@printf("LAPTOP=%s, max n on laptop = %d. Set LAPTOP=false to override.\n",
    LAPTOP, LAPTOP_MAX_N)

generate_family_1d_zz_disordered()
generate_family_1d_clean()
generate_family_2d()

println("\nDone.")
