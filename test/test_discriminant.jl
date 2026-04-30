@testset "Discriminant (Phase 1: primitives)" begin

    # -----------------------------------------------------------------------
    # DiscriminantBuffers
    # -----------------------------------------------------------------------
    @testset "DiscriminantBuffers" begin
        bufs = DiscriminantBuffers(N3_DIM)

        @test bufs isa DiscriminantBuffers{ComplexF64}
        @test size(bufs.work1) == (N3_DIM, N3_DIM)
        @test size(bufs.work2) == (N3_DIM, N3_DIM)
        @test size(bufs.work3) == (N3_DIM, N3_DIM)

        # Parametric construction
        bufs32 = DiscriminantBuffers{ComplexF32}(4)
        @test bufs32 isa DiscriminantBuffers{ComplexF32}
        @test eltype(bufs32.work1) === ComplexF32
    end

    # -----------------------------------------------------------------------
    # gibbs_fractional_powers: identities
    # -----------------------------------------------------------------------
    @testset "gibbs_fractional_powers identities (N3 fixture)" begin
        powers = gibbs_fractional_powers(N3_GIBBS)

        @test powers isa NamedTuple
        @test propertynames(powers) == (:sigma_quarter, :sigma_inv_quarter, :sigma_half)
        @test length(powers.sigma_quarter) == N3_DIM
        @test length(powers.sigma_inv_quarter) == N3_DIM
        @test length(powers.sigma_half) == N3_DIM

        diag_real = real.(diag(N3_GIBBS))

        # sigma_quarter^4 = sigma  (after eps_trunc floor; for N3 fixture
        # the floor is not active, so equality should be exact-up-to-fp).
        @test isapprox(powers.sigma_quarter .^ 4, diag_real; atol=TOL_EXACT)

        # sigma_quarter * sigma_inv_quarter = 1
        @test isapprox(powers.sigma_quarter .* powers.sigma_inv_quarter,
                       ones(Float64, N3_DIM); atol=TOL_EXACT)

        # sigma_half = sigma_quarter^2
        @test isapprox(powers.sigma_half, powers.sigma_quarter .^ 2; atol=TOL_EXACT)

        # sigma_half^2 = sigma
        @test isapprox(powers.sigma_half .^ 2, diag_real; atol=TOL_EXACT)

        # All powers strictly positive (Gibbs has no zero eigenvalues at finite beta)
        @test all(powers.sigma_quarter .> 0)
        @test all(powers.sigma_inv_quarter .> 0)
        @test all(powers.sigma_half .> 0)
    end

    # -----------------------------------------------------------------------
    # eps_trunc floor: defensive behaviour at synthetic low temperature
    # -----------------------------------------------------------------------
    @testset "gibbs_fractional_powers eps_trunc floor" begin
        # Synthetic 4x4 diagonal Gibbs with one entry below eps_trunc=1e-12.
        diag_vals = [0.5, 0.49, 1e-20, 0.01]
        rho = Hermitian(diagm(ComplexF64.(diag_vals)))

        # Default eps_trunc = 1e-12: third entry is floored.
        powers = gibbs_fractional_powers(rho)

        @test powers.sigma_inv_quarter[3] ≈ (1e-12)^(-0.25)
        @test isfinite(powers.sigma_inv_quarter[3])

        # Other entries untouched.
        for i in (1, 2, 4)
            @test powers.sigma_quarter[i] ≈ diag_vals[i]^0.25
            @test powers.sigma_inv_quarter[i] ≈ diag_vals[i]^(-0.25)
        end

        # Lower eps_trunc lets the small entry through.
        powers_loose = gibbs_fractional_powers(rho; eps_trunc=1e-30)
        @test powers_loose.sigma_inv_quarter[3] ≈ (1e-20)^(-0.25)
        @test powers_loose.sigma_inv_quarter[3] > powers.sigma_inv_quarter[3]
    end
end
