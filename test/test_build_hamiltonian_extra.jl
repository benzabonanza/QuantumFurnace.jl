"""
Extra coverage for `build_heis_1d` / `build_tfim_2d` (qf-yi4 verifier).

The existing `test_hamiltonian.jl` covers schema, determinism (1D), BC
propagation, and argument validation. This file adds:

- 2D determinism (byte-identical at same seed).
- 2D non-determinism across seeds.
- TFIM negative-sign convention (-J Σ ZZ - h Σ X) via physical-Hamiltonian
  reconstruction.
- TFIM clean-limit (h=0, disorder=0) spectrum integer-spaced.
- HamHam wrap consumes the extra TFIM fields silently.
"""

using Test
using LinearAlgebra
using Random
using QuantumFurnace
using QuantumFurnace: _construct_2d_heisenberg_base, _construct_disordering_terms,
    X, Y, Z

@testset "qf-yi4: build_heis_1d / build_tfim_2d extra coverage" begin

    # ----------------------------------------------------------------------
    # 2D determinism
    # ----------------------------------------------------------------------
    @testset "build_tfim_2d: same seed gives byte-identical fixture" begin
        raw_a = build_tfim_2d(2, 3; J=1.0, h=1.0, seed=20260516,
            disordering_terms=Vector{Matrix{ComplexF64}}[[Z], [Z, Z]],
            disorder_strength=1e-3)
        raw_b = build_tfim_2d(2, 3; J=1.0, h=1.0, seed=20260516,
            disordering_terms=Vector{Matrix{ComplexF64}}[[Z], [Z, Z]],
            disorder_strength=1e-3)
        @test raw_a.matrix == raw_b.matrix
        @test raw_a.eigvals == raw_b.eigvals
        @test raw_a.disordering_coeffs == raw_b.disordering_coeffs
    end

    @testset "build_tfim_2d: different seeds give different fixtures" begin
        raw_a = build_tfim_2d(2, 3; J=1.0, h=1.0, seed=1,
            disorder_strength=1e-3)
        raw_b = build_tfim_2d(2, 3; J=1.0, h=1.0, seed=2,
            disorder_strength=1e-3)
        @test !isapprox(raw_a.matrix, raw_b.matrix; atol=1e-10)
    end

    # ----------------------------------------------------------------------
    # TFIM sign convention: H_phys = -J Σ ZZ - h Σ X
    # ----------------------------------------------------------------------
    @testset "build_tfim_2d: negative-sign convention recovered exactly" begin
        Lx, Ly = 2, 2
        J, h = 1.0, 1.0
        n = Lx * Ly
        raw = build_tfim_2d(Lx, Ly; J=J, h=h, seed=1, disorder_strength=0.0,
            disordering_terms=Vector{Matrix{ComplexF64}}[[Z]])
        # Recover physical from rescaled storage
        H_recovered = raw.rescaling_factor * (raw.matrix - raw.shift * I(2^n))
        # Construct expected directly
        H_bond_expected = _construct_2d_heisenberg_base(Lx, Ly,
            Vector{Matrix{ComplexF64}}[[Z, Z]], [-J];
            periodic_x=true, periodic_y=true)
        H_field_expected = _construct_disordering_terms(
            Vector{Matrix{ComplexF64}}[[X]], [fill(-h, n)], n)
        H_expected = Matrix(H_bond_expected) + Matrix(H_field_expected)
        @test isapprox(H_recovered, H_expected; atol=1e-10)
    end

    @testset "build_tfim_2d: clean (h=0) limit is diagonal in Z basis" begin
        raw = build_tfim_2d(2, 2; J=1.0, h=0.0, seed=1, disorder_strength=0.0,
            disordering_terms=Vector{Matrix{ComplexF64}}[[Z]])
        off = raw.matrix - Diagonal(diag(raw.matrix))
        @test sum(abs2, off) < 1e-20
    end

    @testset "build_tfim_2d: clean Ising h=0 J=1 spectrum integer-spaced after un-rescale" begin
        # 2x2 PBC: each bond doubled by Ly=2 wrap → effective J=2 over 4 distinct bonds.
        # Spectrum should be (-2J·m) for m ∈ {-4, ..., +4}, i.e. multiples of 4.
        n = 4
        raw = build_tfim_2d(2, 2; J=1.0, h=0.0, seed=1, disorder_strength=0.0,
            disordering_terms=Vector{Matrix{ComplexF64}}[[Z]])
        H_phys = raw.rescaling_factor * (raw.matrix - raw.shift * I(2^n))
        eigs_phys = sort(real.(eigvals(Hermitian(H_phys))))
        # All eigenvalues are integers (multiples of J)
        @test all(abs.(eigs_phys .- round.(eigs_phys)) .< 1e-10)
        # Min/max are ±8 (4 doubled bonds at J=1 ⇒ ZZ eigvals ∈ ±(0,2,4,6,8))
        @test isapprox(minimum(eigs_phys), -8.0; atol=1e-10)
        @test isapprox(maximum(eigs_phys), +8.0; atol=1e-10)
    end

    # ----------------------------------------------------------------------
    # HamHam wrap consumes ALL the extra fields silently
    # ----------------------------------------------------------------------
    @testset "HamHam wrap: build_tfim_2d Lx/Ly/J/h/seed/disorder_strength silently ignored" begin
        raw = build_tfim_2d(2, 3; J=1.5, h=0.5, seed=99,
            periodic_x=true, periodic_y=false,
            disordering_terms=Vector{Matrix{ComplexF64}}[[Z], [Z, Z]],
            disorder_strength=1e-2)
        # All extra fields present in raw
        @test :seed in keys(raw)
        @test :disorder_strength in keys(raw)
        @test :Lx in keys(raw) && :Ly in keys(raw)
        @test :J in keys(raw) && :h in keys(raw)
        # HamHam construction succeeds (silently ignored extras)
        ham = HamHam(raw, 2.0)
        @test ham isa HamHam{Float64}
        @test size(ham.data) == (64, 64)
        # Periodic flag composite: PBC_x AND OBC_y → false
        @test ham.periodic === false
        # β_phys keyword path consumes the same NamedTuple
        ham2 = HamHam(raw; beta_phys=0.5)
        @test isapprox(ham2.gibbs, HamHam(raw, 0.5 * ham2.rescaling_factor).gibbs;
            atol=1e-12)
    end

    # ----------------------------------------------------------------------
    # build_heis_1d: same expectations hold for 1D
    # ----------------------------------------------------------------------
    @testset "build_heis_1d: extra fields present and HamHam silently ignores them" begin
        raw = build_heis_1d(4, [1.0, 1.0, 1.0]; seed=42, disorder_strength=0.1)
        @test :seed in keys(raw)
        @test :disorder_strength in keys(raw)
        @test raw.seed == 42
        @test raw.disorder_strength == 0.1
        ham = HamHam(raw; beta_phys=0.25)
        @test ham isa HamHam{Float64}
        @test isapprox(tr(ham.gibbs), 1.0; atol=1e-12)
    end

    # ----------------------------------------------------------------------
    # Stress: many seeds, n=6 — every realisation must be Hermitian and live in [0, 0.45]
    # ----------------------------------------------------------------------
    @testset "build_heis_1d: 10 seeds × n=6 all give valid fixtures" begin
        for seed in 1:10
            raw = build_heis_1d(6, [1.0, 1.0, 1.0]; seed=seed,
                disordering_terms=Vector{Matrix{ComplexF64}}[[Z], [Z, Z]],
                disorder_strength=0.1)
            @test isapprox(raw.matrix, adjoint(raw.matrix); atol=1e-12)
            @test minimum(raw.eigvals) >= -1e-10
            @test maximum(raw.eigvals) <= 0.45 + 1e-10
            @test raw.nu_min > 0
        end
    end

    # ----------------------------------------------------------------------
    # Disorder semantics: build_tfim_2d uses 2D lattice bonds, NOT 1D-chain
    # ----------------------------------------------------------------------
    @testset "build_tfim_2d: 2-site disorder uses 2D lattice (bond (3,4) absent in 2x3)" begin
        using QuantumFurnace: _construct_disordering_terms_2d
        Lx, Ly = 2, 3
        n = Lx * Ly
        # Construct disorder via the 2D builder vs 1D builder for the SAME coefficients.
        sample_coeffs = [Float64.(1:n)]
        H_2d = Matrix(_construct_disordering_terms_2d(Lx, Ly,
            Vector{Matrix{ComplexF64}}[[Z, Z]], sample_coeffs;
            periodic_x=true, periodic_y=true))
        H_1d_chain = Matrix(_construct_disordering_terms(
            Vector{Matrix{ComplexF64}}[[Z, Z]], sample_coeffs, n; periodic=true))
        @test !isapprox(H_2d, H_1d_chain; atol=1e-10)
        # Bond (3,4) is a 1D-chain bond but NOT a 2D nearest-neighbour bond.
        ZZ_34 = Matrix(QuantumFurnace._pad_two_site_op([Z, Z], n, 3, 4))
        ip_1d = real(tr(H_1d_chain * ZZ_34)) / 2^n
        ip_2d = real(tr(H_2d * ZZ_34)) / 2^n
        @test abs(ip_1d) > 0.1     # 1D chain includes (3,4)
        @test abs(ip_2d) < 1e-12   # 2D lattice does NOT include (3,4)
    end
end
