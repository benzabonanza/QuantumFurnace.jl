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
using QuantumFurnace: X, Y, Z

const OUTPUT_DIR = joinpath(@__DIR__, "output", "multiseed_fixtures")
mkpath(OUTPUT_DIR)

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
