"""
Multi-seed disordered fixture generator (qf-yi4 replacement for find_typical_*).

For each (n, seed) pair we generate ONE ε-disorder realisation of the chosen
Hamiltonian and save it as a BSON in the same NamedTuple schema as the legacy
find_typical_* builders so existing `HamHam(raw; beta_phys=β)` loaders work
unchanged.

Run:
    julia --project scripts/scratch_multiseed_disordered_fixtures.jl

Output:
    scripts/output/multiseed_fixtures/heis_xxx_disordered_periodic_n{n}_seed{seed}.bson
    scripts/output/multiseed_fixtures/tfim_2d_{phase}_{Lx}x{Ly}_seed{seed}.bson

Memory: 1D up to n=8 stays well below 200 MB total disk (5 seeds × 6 n ≈ 150 MB);
2D up to 3×3 stays below 600 MB. RAM peaks at one dense eigendecomp at a time —
n=9 dense Hermitian eigen is ~25 MB transient, far below the sandbox cap.

The "right" protocol per the qf-yi4 follow-up: one fixture per (n, seed), no
selector. Downstream analysis takes the median over seeds and shows the IQR
as a band. No magic spectral-typicality criterion.
"""

using Random
using LinearAlgebra
using Printf
using BSON
using QuantumFurnace
using QuantumFurnace: X, Y, Z,
    _construct_base_ham,
    _construct_disordering_terms,
    _construct_disordering_terms_2d,
    _construct_2d_heisenberg_base,
    _rescaling_and_shift_factors

const OUTPUT_DIR = joinpath(@__DIR__, "output", "multiseed_fixtures")
mkpath(OUTPUT_DIR)

"""
    build_heis_1d(num_qubits, coeffs; seed, periodic=true,
        disordering_terms=[[Z], [Z,Z]], disorder_strength=0.1)

Build ONE disordered 1D Heisenberg fixture at a given seed. No batch, no
selector — the random draw is fully reproducible from `seed`.

Schema matches the legacy find_typical_heisenberg NamedTuple plus extra
`seed` and `disorder_strength` diagnostic fields. `HamHam(raw, β)` ignores
the extras.
"""
function build_heis_1d(num_qubits::Int, coeffs::Vector{Float64};
        seed::Int,
        periodic::Bool=true,
        disordering_terms::Vector{Vector{Matrix{ComplexF64}}}=Vector{Matrix{ComplexF64}}[[Z], [Z, Z]],
        disorder_strength::Float64=0.1)

    base_terms = Vector{Matrix{ComplexF64}}[[X, X], [Y, Y], [Z, Z]]
    base_hamiltonian = _construct_base_ham(base_terms, coeffs, num_qubits; periodic=periodic)

    rng = MersenneTwister(seed)
    sample_coeffs = [zeros(Float64, num_qubits) for _ in disordering_terms]
    for dc in sample_coeffs
        rand!(rng, dc)
        dc .*= disorder_strength
    end
    disordering_ham = _construct_disordering_terms(disordering_terms, sample_coeffs, num_qubits;
        periodic=periodic)

    total_ham = Hermitian(Matrix(base_hamiltonian) + Matrix(disordering_ham))
    rescaling_factor, shift = _rescaling_and_shift_factors(total_ham)
    rescaled_ham = (Matrix(total_ham) ./ rescaling_factor) + shift * I(2^num_qubits)
    rescaled_eigvals, rescaled_eigvecs = eigen(Hermitian(rescaled_ham))
    nu_min = minimum(diff(rescaled_eigvals))

    return (
        matrix = rescaled_ham,
        terms = base_terms,
        base_coeffs = coeffs ./ rescaling_factor,
        disordering_terms = disordering_terms,
        disordering_coeffs = [dc ./ rescaling_factor for dc in sample_coeffs],
        eigvals = rescaled_eigvals,
        eigvecs = rescaled_eigvecs,
        nu_min = nu_min,
        shift = shift,
        rescaling_factor = rescaling_factor,
        periodic = periodic,
        seed = seed,
        disorder_strength = disorder_strength,
    )
end

"""
    build_tfim_2d(Lx, Ly; J=1.0, h=1.0, seed, periodic_x=true, periodic_y=true,
        disordering_terms=[[Z], [Z,Z]], disorder_strength=1e-3)

Build ONE 2D transverse-field Ising fixture
    H = -J Σ_{<i,j>} Z_i Z_j − h Σ_i X_i + ε-disorder
at the given seed. The disorder uses the new `_construct_disordering_terms_2d`
so two-site disorder terms ride actual lattice bonds (right + up neighbour
per site), not 1D-chain bonds on the linearised index.

Schema matches find_typical_2d_heisenberg.
"""
function build_tfim_2d(Lx::Int, Ly::Int;
        J::Float64=1.0, h::Float64=1.0,
        seed::Int,
        periodic_x::Bool=true, periodic_y::Bool=true,
        disordering_terms::Vector{Vector{Matrix{ComplexF64}}}=Vector{Matrix{ComplexF64}}[[Z], [Z, Z]],
        disorder_strength::Float64=1e-3)

    num_qubits = Lx * Ly
    H_bond = _construct_2d_heisenberg_base(Lx, Ly,
        Vector{Matrix{ComplexF64}}[[Z, Z]], [-J];
        periodic_x=periodic_x, periodic_y=periodic_y)
    # Transverse field as a uniform per-site X coefficient.
    field_coeffs = fill(-h, num_qubits)
    H_field = _construct_disordering_terms(
        Vector{Matrix{ComplexF64}}[[X]], [field_coeffs], num_qubits)
    base_clean = Hermitian(Matrix(H_bond) + Matrix(H_field))

    rng = MersenneTwister(seed)
    sample_coeffs = [zeros(Float64, num_qubits) for _ in disordering_terms]
    for dc in sample_coeffs
        rand!(rng, dc)
        dc .*= disorder_strength
    end
    disordering_ham = _construct_disordering_terms_2d(Lx, Ly,
        disordering_terms, sample_coeffs;
        periodic_x=periodic_x, periodic_y=periodic_y)

    total_ham = Hermitian(Matrix(base_clean) + Matrix(disordering_ham))
    rescaling_factor, shift = _rescaling_and_shift_factors(total_ham)
    rescaled_ham = (Matrix(total_ham) ./ rescaling_factor) + shift * I(2^num_qubits)
    rescaled_eigvals, rescaled_eigvecs = eigen(Hermitian(rescaled_ham))
    nu_min = minimum(diff(rescaled_eigvals))

    return (
        matrix = rescaled_ham,
        terms = Vector{Matrix{ComplexF64}}[[Z, Z], [X]],
        base_coeffs = [-J / rescaling_factor, -h / rescaling_factor],
        disordering_terms = disordering_terms,
        disordering_coeffs = [dc ./ rescaling_factor for dc in sample_coeffs],
        eigvals = rescaled_eigvals,
        eigvecs = rescaled_eigvecs,
        nu_min = nu_min,
        shift = shift,
        rescaling_factor = rescaling_factor,
        periodic = periodic_x && periodic_y,
        seed = seed,
        disorder_strength = disorder_strength,
        Lx = Lx, Ly = Ly,
        J = J, h = h,
    )
end

# -------------------------------------------------------------------------
# Generation driver
# -------------------------------------------------------------------------

const SEEDS = (42, 43, 44, 45, 46)  # five independent realisations per cell

function run_1d_heisenberg()
    @printf("\n=== 1D Heisenberg (XXX with [Z]+[Z,Z] disorder, ε=0.1) ===\n")
    coeffs = [1.0, 1.0, 1.0]
    for n in 3:8
        for seed in SEEDS
            raw = build_heis_1d(n, coeffs; seed=seed, periodic=true,
                disorder_strength=0.1)
            path = joinpath(OUTPUT_DIR,
                @sprintf("heis_xxx_disordered_periodic_n%d_seed%d.bson", n, seed))
            BSON.bson(path, hamiltonian=raw)
            @printf("  n=%d seed=%d  nu_min=%.4e  R=%.3f  → %s\n",
                n, seed, raw.nu_min, raw.rescaling_factor, basename(path))
        end
    end
end

function run_2d_tfim()
    # qf-1jj: ordered (h=1, β_phys=2 → T/T_c ≈ 0.24) vs disordered (h=1, β_phys=0.25
    # → T/T_c ≈ 1.93). We save the fixture once per (Lx, Ly, seed); the operating
    # point (β_phys) is applied at HamHam construction time downstream so a single
    # raw fixture serves both phases. We do however choose h to put each fixture
    # in its target phase at the canonical β_phys:
    #   ordered:    h = 1.0  (T_c(h=1) ≈ 2.07 ⇒ β_phys = 2 → T ≈ 0.5 ≪ T_c)
    #   disordered: h = 3.5  (h_c ≈ 3.044 ⇒ ground state already disordered at any T)
    #
    # 5 seeds × 3 lattices × 2 phases = 30 fixtures, ≤ 600 MB on disk.
    @printf("\n=== 2D TFIM (-J ZZ + -h X with [Z]+[Z,Z] ε-disorder, ε=1e-3) ===\n")
    for (Lx, Ly) in ((2, 2), (2, 3), (3, 3))
        for (phase, h) in (("ordered", 1.0), ("disordered", 3.5))
            for seed in SEEDS
                raw = build_tfim_2d(Lx, Ly; J=1.0, h=h, seed=seed,
                    periodic_x=true, periodic_y=true,
                    disorder_strength=1e-3)
                path = joinpath(OUTPUT_DIR,
                    @sprintf("tfim_2d_%s_%dx%d_seed%d.bson", phase, Lx, Ly, seed))
                BSON.bson(path, hamiltonian=raw, phase=phase, h=h)
                @printf("  %s %dx%d seed=%d  nu_min=%.4e  R=%.3f  → %s\n",
                    phase, Lx, Ly, seed, raw.nu_min, raw.rescaling_factor,
                    basename(path))
            end
        end
    end
end

function main()
    @printf("Output dir: %s\n", OUTPUT_DIR)
    run_1d_heisenberg()
    run_2d_tfim()
    @printf("\nDone.\n")
end

main()
