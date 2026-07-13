"""
Tests for Heisenberg Hamiltonian builders in `src/hamiltonian.jl`.

Covers the underlying primitives (`_pad_two_site_op`,
`_construct_2d_heisenberg_base`, `_construct_disordering_terms*`) plus
the public seed-driven builders [`build_heis_1d`] and [`build_tfim_2d`]
that replaced the find_typical_* / find_ideal_* spectral-selector path
(qf-yi4, 2026-05-15). HamHam constructor coverage (direct, single-term,
multi-term with periodic) closes out the file.
"""

using LinearAlgebra
using SparseArrays
using Random
using Statistics: median

@testset "Heisenberg Hamiltonian builders" begin

    # _pad_two_site_op smoke tests --------------------------------------------------
    @testset "_pad_two_site_op: adjacent and non-adjacent placements" begin
        # n=3 adjacent at the chain edge: place X at q=1, q=2
        op = QuantumFurnace._pad_two_site_op([X, X], 3, 1, 2)
        @test size(op) == (8, 8)
        @test ishermitian(Matrix(op))
        @test Matrix(op) ≈ kron(X, X, I(2))

        # n=4 adjacent: place X at q=2, q=3
        op = QuantumFurnace._pad_two_site_op([X, X], 4, 2, 3)
        @test Matrix(op) ≈ kron(I(2), X, X, I(2))

        # n=4 non-adjacent: place Z at q=1, q=4 (separation 3)
        op = QuantumFurnace._pad_two_site_op([Z, Z], 4, 1, 4)
        @test Matrix(op) ≈ kron(Z, I(2), I(2), Z)

        # Order-independence for symmetric terms: q1<q2 vs q1>q2
        op_a = QuantumFurnace._pad_two_site_op([Z, Z], 4, 1, 4)
        op_b = QuantumFurnace._pad_two_site_op([Z, Z], 4, 4, 1)
        @test Matrix(op_a) ≈ Matrix(op_b)
    end

    @testset "_pad_two_site_op: argument validation" begin
        @test_throws ArgumentError QuantumFurnace._pad_two_site_op([X], 4, 1, 2)         # 1-site term
        @test_throws ArgumentError QuantumFurnace._pad_two_site_op([X, X], 4, 2, 2)      # q1 == q2
        @test_throws ArgumentError QuantumFurnace._pad_two_site_op([X, X], 4, 0, 2)      # q < 1
        @test_throws ArgumentError QuantumFurnace._pad_two_site_op([X, X], 4, 2, 5)      # q > num_qubits
    end

    # _construct_2d_heisenberg_base ----------------------------------------------------
    @testset "_construct_2d_heisenberg_base: 3x1 and 1x3 lattices match 1D periodic n=3" begin
        terms_xyz = [[X, X], [Y, Y], [Z, Z]]
        coeffs = [1.0, 1.0, 1.0]

        ham_1d = QuantumFurnace._construct_base_ham(terms_xyz, coeffs, 3; periodic=true)

        # 3x1 lattice with full periodic BC: x-direction wraps, y-direction has no bonds (Ly=1)
        ham_3x1 = QuantumFurnace._construct_2d_heisenberg_base(3, 1, terms_xyz, coeffs;
            periodic_x=true, periodic_y=true)
        @test Matrix(ham_3x1) ≈ Matrix(ham_1d)

        # 1x3 lattice with full periodic BC: y-direction wraps, x-direction has no bonds (Lx=1)
        ham_1x3 = QuantumFurnace._construct_2d_heisenberg_base(1, 3, terms_xyz, coeffs;
            periodic_x=true, periodic_y=true)
        @test Matrix(ham_1x3) ≈ Matrix(ham_1d)
    end

    @testset "_construct_2d_heisenberg_base: dimension and Hermiticity for several lattices" begin
        for (Lx, Ly) in [(2, 3), (3, 3), (2, 5)]
            n = Lx * Ly
            ham = QuantumFurnace._construct_2d_heisenberg_base(Lx, Ly,
                [[X, X], [Y, Y], [Z, Z]], [1.0, 1.0, 1.5];
                periodic_x=true, periodic_y=true)
            @test size(ham) == (2^n, 2^n)
            @test ishermitian(Matrix(ham))
            # Heisenberg model is traceless
            @test abs(tr(Matrix(ham))) < 1e-10
        end
    end

    @testset "_construct_2d_heisenberg_base: bond counting via Frobenius norm" begin
        # Distinct ZZ-bond Pauli strings are pairwise Hilbert–Schmidt-orthogonal
        # (every product of two distinct ZZ-strings is a higher-weight Pauli, traceless),
        # so for H = Σ_b c_b · Z_{q1_b} Z_{q2_b}, we have ‖H‖_F² = (Σ_b |c_b|²) · 2^n.
        #
        # For Lx == 2 (or Ly == 2) with periodic BC, the wrap-around bond coincides
        # with the original bond, so each distinct bond is added twice with coefficient 1
        # → effective coefficient 2 per distinct bond → contribution 4 × #distinct.
        function expected_sumsq(Lx, Ly; periodic_x=true, periodic_y=true)
            x_contrib = if Lx == 1
                0
            elseif Lx == 2 && periodic_x
                4 * Ly         # Ly distinct bonds, coeff 2 each → 4·Ly
            elseif periodic_x
                Lx * Ly        # Lx·Ly distinct bonds, coeff 1
            else
                (Lx - 1) * Ly  # OBC: (Lx−1)·Ly distinct bonds, coeff 1
            end
            y_contrib = if Ly == 1
                0
            elseif Ly == 2 && periodic_y
                4 * Lx
            elseif periodic_y
                Lx * Ly
            else
                Lx * (Ly - 1)
            end
            return x_contrib + y_contrib
        end

        for (Lx, Ly) in [(2, 3), (3, 3), (2, 5)]
            n = Lx * Ly
            ham_zz = QuantumFurnace._construct_2d_heisenberg_base(Lx, Ly, [[Z, Z]], [1.0];
                periodic_x=true, periodic_y=true)
            fro2 = sum(abs2, Matrix(ham_zz))
            @test fro2 ≈ expected_sumsq(Lx, Ly) * 2^n  rtol=1e-10
        end
    end

    @testset "_construct_2d_heisenberg_base: open boundary disables wrap" begin
        # 3x3 OBC: (Lx-1)*Ly + Lx*(Ly-1) = 6 + 6 = 12 distinct ZZ-bonds, coeff 1 each
        # 3x3 PBC: 9 + 9 = 18 (periodic_x adds 3 wrap bonds, periodic_y adds 3)
        ham_pbc_3 = QuantumFurnace._construct_2d_heisenberg_base(3, 3, [[Z, Z]], [1.0];
            periodic_x=true, periodic_y=true)
        ham_obc_3 = QuantumFurnace._construct_2d_heisenberg_base(3, 3, [[Z, Z]], [1.0];
            periodic_x=false, periodic_y=false)
        @test sum(abs2, Matrix(ham_pbc_3)) > sum(abs2, Matrix(ham_obc_3))
        @test sum(abs2, Matrix(ham_obc_3)) ≈ 12 * 2^9  rtol=1e-10
        @test sum(abs2, Matrix(ham_pbc_3)) ≈ 18 * 2^9  rtol=1e-10
    end

    # build_heis_1d / build_tfim_2d (qf-yi4 replacements for find_*) -----------------
    @testset "build_heis_1d: returns a valid raw NamedTuple" begin
        raw = build_heis_1d(3, [1.0, 1.0, 1.0]; seed=20260515,
            disordering_terms=Vector{Matrix{ComplexF64}}[[Z]], disorder_strength=1.0)
        @test raw.nu_min > 0
        @test size(raw.matrix) == (8, 8)
        @test size(raw.eigvecs) == (8, 8)
        @test length(raw.eigvals) == 8
        @test raw.periodic === true
        @test length(raw.disordering_terms) == 1
        @test length(raw.disordering_coeffs) == 1
        @test length(raw.disordering_coeffs[1]) == 3
        @test raw.seed === 20260515
        @test raw.disorder_strength == 1.0
        @test minimum(raw.eigvals) ≥ -1e-10
        @test maximum(raw.eigvals) ≤ 0.45 + 1e-10
        @test isapprox(raw.matrix, raw.matrix'; atol=1e-12)
    end

    @testset "build_heis_1d: HamHam wrap end-to-end (extra fields ignored)" begin
        raw = build_heis_1d(3, [1.0, 1.0, 1.0]; seed=20260516, disorder_strength=1.0)
        ham = HamHam(raw, 1.0)
        @test ham isa HamHam{Float64}
        @test size(ham.data) == (8, 8)
        @test ham.nu_min > 0
        @test isapprox(tr(ham.gibbs), 1.0; atol=1e-10)
    end

    @testset "build_heis_1d: same seed gives identical fixture" begin
        raw_a = build_heis_1d(4, [1.0, 1.0, 1.0]; seed=20260517,
            disordering_terms=Vector{Matrix{ComplexF64}}[[Z], [Z, Z]],
            disorder_strength=0.5)
        raw_b = build_heis_1d(4, [1.0, 1.0, 1.0]; seed=20260517,
            disordering_terms=Vector{Matrix{ComplexF64}}[[Z], [Z, Z]],
            disorder_strength=0.5)
        @test raw_a.eigvals ≈ raw_b.eigvals
        @test isapprox(raw_a.matrix, raw_b.matrix; atol=1e-14)
    end

    @testset "build_heis_1d: different seeds give different fixtures" begin
        raw_a = build_heis_1d(4, [1.0, 1.0, 1.0]; seed=1, disorder_strength=0.5)
        raw_b = build_heis_1d(4, [1.0, 1.0, 1.0]; seed=2, disorder_strength=0.5)
        @test !isapprox(raw_a.matrix, raw_b.matrix; atol=1e-8)
    end

    @testset "build_heis_1d: disorder_strength bounds the per-coefficient magnitude" begin
        raw = build_heis_1d(3, [1.0, 1.0, 1.0]; seed=20260518,
            disordering_terms=Vector{Matrix{ComplexF64}}[[Z]], disorder_strength=1e-2)
        # raw.disordering_coeffs is the *rescaled* version; each rescaled entry is at
        # most `disorder_strength / rescaling_factor` < 1e-2.
        @test all(abs.(raw.disordering_coeffs[1]) .≤ 1e-2)
    end

    @testset "build_heis_1d: OBC vs PBC produce different fixtures" begin
        raw_pbc = build_heis_1d(4, [1.0, 1.0, 1.0]; seed=20260519,
            disordering_terms=Vector{Matrix{ComplexF64}}[[Z], [Z, Z]],
            disorder_strength=1.0, periodic=true)
        raw_obc = build_heis_1d(4, [1.0, 1.0, 1.0]; seed=20260519,
            disordering_terms=Vector{Matrix{ComplexF64}}[[Z], [Z, Z]],
            disorder_strength=1.0, periodic=false)
        @test raw_pbc.periodic === true
        @test raw_obc.periodic === false
        @test !isapprox(raw_pbc.matrix, raw_obc.matrix; atol=1e-8)
    end

    @testset "build_tfim_2d: returns a valid raw NamedTuple" begin
        raw = build_tfim_2d(2, 2; J=1.0, h=1.0, seed=20260520,
            periodic_x=true, periodic_y=true,
            disordering_terms=Vector{Matrix{ComplexF64}}[[Z]],
            disorder_strength=1e-2)
        @test raw.nu_min > 0
        @test size(raw.matrix) == (16, 16)
        @test length(raw.eigvals) == 16
        @test raw.periodic === true
        @test raw.Lx == 2
        @test raw.Ly == 2
        @test raw.J == 1.0
        @test raw.h == 1.0
        @test length(raw.disordering_coeffs[1]) == 4
        @test minimum(raw.eigvals) ≥ -1e-10
        @test maximum(raw.eigvals) ≤ 0.45 + 1e-10
    end

    @testset "build_tfim_2d: HamHam wrap end-to-end" begin
        raw = build_tfim_2d(2, 2; J=1.0, h=1.0, seed=20260521, disorder_strength=1e-3)
        ham = HamHam(raw, 2.0)
        @test ham isa HamHam{Float64}
        @test size(ham.data) == (16, 16)
        @test isapprox(tr(ham.gibbs), 1.0; atol=1e-10)
    end

    @testset "build_tfim_2d: mixed periodic_x / periodic_y all distinct" begin
        rawpp = build_tfim_2d(2, 3; J=1.0, h=1.5, seed=20260522,
            periodic_x=true,  periodic_y=true,  disorder_strength=1e-3)
        rawpf = build_tfim_2d(2, 3; J=1.0, h=1.5, seed=20260522,
            periodic_x=true,  periodic_y=false, disorder_strength=1e-3)
        rawfp = build_tfim_2d(2, 3; J=1.0, h=1.5, seed=20260522,
            periodic_x=false, periodic_y=true,  disorder_strength=1e-3)
        rawff = build_tfim_2d(2, 3; J=1.0, h=1.5, seed=20260522,
            periodic_x=false, periodic_y=false, disorder_strength=1e-3)
        @test !isapprox(rawpp.matrix, rawpf.matrix; atol=1e-8)
        @test !isapprox(rawpp.matrix, rawfp.matrix; atol=1e-8)
        @test !isapprox(rawpp.matrix, rawff.matrix; atol=1e-8)
        @test !isapprox(rawpf.matrix, rawfp.matrix; atol=1e-8)
    end

    @testset "build_tfim_2d: argument validation" begin
        @test_throws ArgumentError build_tfim_2d(0, 2; seed=1)
        @test_throws ArgumentError build_tfim_2d(2, 0; seed=1)
    end


    @testset "HamHam ctor (no disorder): direct coverage" begin
        # Ctor (1) was reachable only transitively via the NamedTuple ctor (3)
        # before this test. Smoke-test the direct path.
        n = 3
        terms = Vector{Matrix{ComplexF64}}[[X, X], [Y, Y], [Z, Z]]
        coeffs = [1.0, 1.0, 1.0]
        h = HamHam(terms, coeffs, n, 1.0)
        @test h.periodic === true
        @test h.disordering_terms === nothing
        @test h.disordering_coeffs === nothing
        @test size(h.data) == (2^n, 2^n)
        @test isapprox(tr(h.gibbs), 1.0; atol=1e-12)
    end

    @testset "HamHam ctor (single-term convenience): wraps to multi-term" begin
        # Ctor (2b): single-term sugar over (2). Verify the wrapped multi-term
        # storage (disordering_terms is a 1-element vector).
        n = 3
        terms = Vector{Matrix{ComplexF64}}[[X, X], [Y, Y], [Z, Z]]
        coeffs = [1.0, 1.0, 1.0]
        dis_term = Matrix{ComplexF64}[Z]           # singular term: Vector{Matrix}
        dis_coeffs = [0.1, -0.1, 0.05]             # singular coeff vector
        h = HamHam(terms, coeffs, dis_term, dis_coeffs, n, 1.0)
        @test h.disordering_terms isa Vector
        @test length(h.disordering_terms) == 1
        @test length(h.disordering_coeffs) == 1
        @test isapprox(tr(h.gibbs), 1.0; atol=1e-12)
    end

    @testset "HamHam ctor (multi-term disorder) forwards `periodic` kwarg" begin
        # Regression for qf-fzj.3 audit §5b: the multi-term ctor previously dropped
        # `periodic` before calling _construct_base_ham, silently building a periodic
        # base even when periodic=false was requested. With zero disorder, the spectrum
        # of an OBC chain differs from a PBC chain (no wrap bond) — after rescaling
        # both fit [0, 0.5] but the rescaled matrices are not equal.
        n = 3
        terms = Vector{Matrix{ComplexF64}}[[X, X], [Y, Y], [Z, Z]]
        coeffs = [1.0, 1.0, 1.0]
        dis_terms = Vector{Matrix{ComplexF64}}[[Z]]
        dis_coeffs = [zeros(n)]
        H_per = HamHam(terms, coeffs, dis_terms, dis_coeffs, n, 1.0; periodic=true)
        H_obc = HamHam(terms, coeffs, dis_terms, dis_coeffs, n, 1.0; periodic=false)
        @test H_per.periodic === true
        @test H_obc.periodic === false
        @test !isapprox(H_per.data, H_obc.data; atol=1e-8)
    end
end
