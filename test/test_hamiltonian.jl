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

    @testset "Public Pauli conversion and grouping helpers" begin
        @test isapprox(Had' * Had, I(2); atol=1e-15)
        converted = pauli_string_to_matrix(["X", "Y", "Z", "I"])
        @test converted[1] == X
        @test converted[2] == Y
        @test converted[3] == Z
        @test converted[4] == Matrix{ComplexF64}(I, 2, 2)
        @test_throws KeyError pauli_string_to_matrix(["not-a-Pauli"])

        terms = Vector{Matrix{ComplexF64}}[[X, X], [Y, Y], [X, Y], [Z]]
        grouped = group_hamiltonian_terms(HamHam(terms, [1.0, 2.0, 3.0, 4.0], 3, 1.0))
        @test grouped.commuting[1] == terms[1:2]
        @test grouped.noncommuting[1] == terms[3:3]
        @test grouped.one_sites[1] == terms[4:4]
    end

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

    @testset "build_tfim_2d: deterministic disorder and physical conventions" begin
        kwargs = (;
            J=1.0, h=1.0, disordering_terms=Vector{Matrix{ComplexF64}}[[Z], [Z, Z]],
            disorder_strength=1e-3,
        )
        raw_a = build_tfim_2d(2, 3; kwargs..., seed=20260516)
        raw_b = build_tfim_2d(2, 3; kwargs..., seed=20260516)
        raw_c = build_tfim_2d(2, 3; kwargs..., seed=20260517)
        @test raw_a.matrix == raw_b.matrix
        @test raw_a.eigvals == raw_b.eigvals
        @test raw_a.disordering_coeffs == raw_b.disordering_coeffs
        @test !isapprox(raw_a.matrix, raw_c.matrix; atol=1e-10)

        Lx = Ly = 2
        n = Lx * Ly
        raw = build_tfim_2d(Lx, Ly; J=1.0, h=1.0, seed=1,
            disordering_terms=Vector{Matrix{ComplexF64}}[[Z]], disorder_strength=0.0)
        H_physical = raw.rescaling_factor * (raw.matrix - raw.shift * I(2^n))
        H_bonds = QuantumFurnace._construct_2d_heisenberg_base(
            Lx, Ly, Vector{Matrix{ComplexF64}}[[Z, Z]], [-1.0];
            periodic_x=true, periodic_y=true)
        H_field = QuantumFurnace._construct_disordering_terms(
            Vector{Matrix{ComplexF64}}[[X]], [fill(-1.0, n)], n)
        @test isapprox(H_physical, Matrix(H_bonds) + Matrix(H_field); atol=1e-10)

        ham_phys = HamHam(raw; beta_phys=0.5)
        ham_alg = HamHam(raw, 0.5 * ham_phys.rescaling_factor)
        @test ham_phys.gibbs == ham_alg.gibbs
    end

    @testset "build_tfim_2d: clean Ising structure and spectrum" begin
        n = 4
        raw = build_tfim_2d(2, 2; J=1.0, h=0.0, seed=1,
            disordering_terms=Vector{Matrix{ComplexF64}}[[Z]], disorder_strength=0.0)
        @test raw.matrix == Diagonal(diag(raw.matrix))

        H_physical = raw.rescaling_factor * (raw.matrix - raw.shift * I(2^n))
        eigs = sort(real.(eigvals(Hermitian(H_physical))))
        @test all(abs.(eigs ./ 4 .- round.(eigs ./ 4)) .< 1e-10)
        @test minimum(eigs) ≈ -8.0 atol=1e-10
        @test maximum(eigs) ≈ 8.0 atol=1e-10
    end

    @testset "build_tfim_2d: two-site disorder follows lattice bonds" begin
        Lx, Ly = 2, 3
        n = Lx * Ly
        coeffs = [Float64.(1:n)]
        H_2d = Matrix(QuantumFurnace._construct_disordering_terms_2d(
            Lx, Ly, Vector{Matrix{ComplexF64}}[[Z, Z]], coeffs;
            periodic_x=true, periodic_y=true))
        H_chain = Matrix(QuantumFurnace._construct_disordering_terms(
            Vector{Matrix{ComplexF64}}[[Z, Z]], coeffs, n; periodic=true))
        @test H_2d != H_chain

        ZZ_34 = Matrix(QuantumFurnace._pad_two_site_op([Z, Z], n, 3, 4))
        @test abs(real(tr(H_chain * ZZ_34)) / 2^n) > 0.1
        @test abs(real(tr(H_2d * ZZ_34)) / 2^n) < 1e-12
    end
end
