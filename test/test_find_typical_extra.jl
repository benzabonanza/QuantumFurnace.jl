"""
Additional coverage for `find_typical_heisenberg` / `find_typical_2d_heisenberg`
(W₂-spectral-median selector). Written by the code-verifier agent to close
gaps in the existing test_hamiltonian.jl coverage:

  (a) multi-term disorder (`[[Z], [Z, Z]]`) end-to-end + HamHam wrap
  (b) OBC vs PBC propagation: `periodic=false` must produce a different
      Hamiltonian from `periodic=true` (matches PBC/OBC distinction for the
      ideal selector).
  (c) `disorder_strength` scaling propagates: ε=1e-3 → stored coeffs all
      ≪ ε=1.0 (modulo rescaling).
  (d) `recompute eigen` step matches the spectrum that scored: rebuild the
      selected Hamiltonian from `disordering_coeffs` and confirm eigvals
      agree to floating-point precision.
  (e) Independent seeds → different but valid spectra (anti-determinism
      complement to the existing seeded determinism test).
  (f) 2D `periodic = periodic_x && periodic_y` (mixed BC).

All assertions run on n=3 sandbox-fit fixtures.
"""

using LinearAlgebra
using Random
using Statistics: median

@testset "find_typical extras" begin

    @testset "find_typical_heisenberg: multi-term disorder propagation" begin
        # [[Z], [Z, Z]] disorder: two separate coefficient vectors, both length=n.
        # The W₂ kernel should accept this and produce a valid HamHam wrap.
        Random.seed!(20260520)
        n = 3
        dis_terms = Vector{Matrix{ComplexF64}}[[Z], [Z, Z]]
        raw = find_typical_heisenberg(n, [1.0, 1.0, 1.0]; batch_size=30,
            disordering_terms=dis_terms, disorder_strength=1.0)
        @test length(raw.disordering_terms) == 2
        @test length(raw.disordering_coeffs) == 2
        @test length(raw.disordering_coeffs[1]) == n
        @test length(raw.disordering_coeffs[2]) == n
        ham = HamHam(raw, 1.0)
        @test ham isa HamHam{Float64}
        @test isapprox(tr(ham.gibbs), 1.0; atol=1e-10)
        @test ham.disordering_terms !== nothing
        @test length(ham.disordering_terms) == 2
    end

    @testset "find_typical_heisenberg: PBC vs OBC produce distinct matrices" begin
        Random.seed!(20260521)
        raw_pbc = find_typical_heisenberg(3, [1.0, 1.0, 1.0]; batch_size=20,
            periodic=true,
            disordering_terms=Vector{Matrix{ComplexF64}}[[Z]],
            disorder_strength=1.0)
        Random.seed!(20260521)
        raw_obc = find_typical_heisenberg(3, [1.0, 1.0, 1.0]; batch_size=20,
            periodic=false,
            disordering_terms=Vector{Matrix{ComplexF64}}[[Z]],
            disorder_strength=1.0)
        @test raw_pbc.periodic === true
        @test raw_obc.periodic === false
        # Matrices must differ because the base Hamiltonian differs (PBC adds
        # a wrap bond at n=3).
        @test !isapprox(raw_pbc.matrix, raw_obc.matrix; atol=1e-8)
    end

    @testset "find_typical_heisenberg: disorder_strength scales coefficients" begin
        Random.seed!(20260522)
        raw_low = find_typical_heisenberg(3, [1.0, 1.0, 1.0]; batch_size=20,
            disordering_terms=Vector{Matrix{ComplexF64}}[[Z]],
            disorder_strength=1e-3)
        Random.seed!(20260522)
        raw_high = find_typical_heisenberg(3, [1.0, 1.0, 1.0]; batch_size=20,
            disordering_terms=Vector{Matrix{ComplexF64}}[[Z]],
            disorder_strength=1.0)
        # The stored coefficients are post-rescaling, so the upper bound is
        # disorder_strength / rescaling_factor. The rescaling factor depends
        # on disorder, but for the same Heisenberg base the rescaling is
        # O(1)-comparable across the two ε runs → the ratio of stored
        # magnitudes should be ≫ 100.
        mag_low = maximum(abs.(raw_low.disordering_coeffs[1]))
        mag_high = maximum(abs.(raw_high.disordering_coeffs[1]))
        @test mag_high > 100 * mag_low
    end

    @testset "find_typical_heisenberg: recomputed eigvals match stored eigvals" begin
        # The final eigen() pass on the chosen sample must reconstruct the
        # same Hamiltonian that was scored during the sweep — otherwise the
        # returned eigvals/eigvecs would disagree with the typicality_distance.
        Random.seed!(20260523)
        n = 3
        raw = find_typical_heisenberg(n, [1.0, 1.0, 1.0]; batch_size=20,
            disordering_terms=Vector{Matrix{ComplexF64}}[[Z]],
            disorder_strength=1.0)
        # Rebuild manually using the returned (pre-rescaled) coefficients
        base = QuantumFurnace._construct_base_ham(
            Vector{Matrix{ComplexF64}}[[X, X], [Y, Y], [Z, Z]],
            [1.0, 1.0, 1.0], n; periodic=true)
        # The returned disordering_coeffs are POST-rescaling → undo by ×rsf
        pre_dc = [c .* raw.rescaling_factor for c in raw.disordering_coeffs]
        dh = QuantumFurnace._construct_disordering_terms(
            Vector{Matrix{ComplexF64}}[[Z]], pre_dc, n)
        total = base + dh
        rsf, sh = QuantumFurnace._rescaling_and_shift_factors(total)
        rescaled = (total ./ rsf) + sh * I
        ev_manual = eigvals(Hermitian(rescaled))
        @test isapprox(ev_manual, raw.eigvals; atol=1e-12)
        @test isapprox(rsf, raw.rescaling_factor; atol=1e-14)
        @test isapprox(sh, raw.shift; atol=1e-14)
    end

    @testset "find_typical_heisenberg: independent seeds → different valid results" begin
        Random.seed!(20260524)
        r1 = find_typical_heisenberg(3, [1.0, 1.0, 1.0]; batch_size=20,
            disordering_terms=Vector{Matrix{ComplexF64}}[[Z]],
            disorder_strength=1.0)
        Random.seed!(20260525)
        r2 = find_typical_heisenberg(3, [1.0, 1.0, 1.0]; batch_size=20,
            disordering_terms=Vector{Matrix{ComplexF64}}[[Z]],
            disorder_strength=1.0)
        # Different seeds → different eigvals
        @test !isapprox(r1.eigvals, r2.eigvals; atol=1e-8)
        # Both spectra valid in [0, 0.45]
        @test minimum(r1.eigvals) ≥ -1e-10 && maximum(r1.eigvals) ≤ 0.45 + 1e-10
        @test minimum(r2.eigvals) ≥ -1e-10 && maximum(r2.eigvals) ≤ 0.45 + 1e-10
        @test r1.typicality_distance >= 0
        @test r2.typicality_distance >= 0
    end

    @testset "find_typical_2d_heisenberg: mixed BC (periodic field = AND)" begin
        Random.seed!(20260526)
        # 2x3 mixed: x periodic, y open
        raw = find_typical_2d_heisenberg(2, 3, [1.0, 1.0, 1.5]; batch_size=15,
            periodic_x=true, periodic_y=false,
            disordering_terms=Vector{Matrix{ComplexF64}}[[Z]],
            disorder_strength=1e-2)
        @test raw.periodic === false  # AND of mixed BC
        n = 6
        @test size(raw.matrix) == (2^n, 2^n)
        @test length(raw.disordering_coeffs[1]) == n
        @test raw.typicality_distance >= 0
    end

    @testset "find_typical_heisenberg: typicality_distance is a lower bound across batch" begin
        # Stronger version of the existing argmin test: verify the returned
        # distance is min over *all* valid samples by computing every sample's
        # distance independently and confirming the returned one is the global min.
        Random.seed!(20260527)
        n = 3
        d = 2^n
        B = 25
        dis_terms = Vector{Matrix{ComplexF64}}[[Z]]
        coeffs = [1.0, 1.0, 1.0]
        strength = 1.0

        # Independent ref loop tracking every sample's distance
        Random.seed!(20260527)
        base_ref = QuantumFurnace._construct_base_ham(
            Vector{Matrix{ComplexF64}}[[X, X], [Y, Y], [Z, Z]], coeffs, n; periodic=true)
        specs = Matrix{Float64}(undef, d, B)
        for k in 1:B
            sc = [zeros(Float64, n) for _ in dis_terms]
            for dc in sc; rand!(dc); dc .*= strength; end
            dh = QuantumFurnace._construct_disordering_terms(dis_terms, sc, n)
            total = base_ref + dh
            rsf, sh = QuantumFurnace._rescaling_and_shift_factors(total)
            rescaled = (total ./ rsf) + sh * I
            ev = eigvals(Hermitian(rescaled))
            bw = ev[end] - ev[1]
            specs[:, k] = (ev .- ev[1]) ./ bw
        end
        med = [median(@view specs[i, :]) for i in 1:d]
        all_dists = [sqrt(sum((specs[:, k] .- med).^2)) for k in 1:B]
        min_dist = minimum(all_dists)

        # Now run find_typical with the same seed and confirm
        Random.seed!(20260527)
        raw = find_typical_heisenberg(n, coeffs; batch_size=B,
            disordering_terms=dis_terms, disorder_strength=strength)

        # The returned distance must equal the global min
        @test isapprox(raw.typicality_distance, min_dist; rtol=1e-12)
        # And must be ≤ every other sample's distance
        @test all(raw.typicality_distance .≤ all_dists .+ 1e-14)
    end

end
