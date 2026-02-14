using Test
using LinearAlgebra

# CPTP completeness verification (TVAL-01)
# Verifies: K0'*K0 + delta*R + U_res'*U_res = I
#
# Derivation (Chen 2023, Theorem III.1):
#   K0 = I - alpha*R, alpha = 1 - sqrt(1-delta)
#   S = (2*alpha - delta)*R - alpha^2*R^2
#   U_res'*U_res = S (by construction)
#   K0'*K0 = I - 2*alpha*R + alpha^2*R^2
#   K0'*K0 + delta*R + S = I - 2*alpha*R + alpha^2*R^2 + delta*R + (2*alpha - delta)*R - alpha^2*R^2 = I
#
# Tolerance: 1e-10 (per user decision -- allows small numerical accumulation)

@testset "CPTP Completeness (TVAL-01)" begin

    @testset "EnergyDomain" begin
        config = make_thermalize_config(EnergyDomain(); delta=TEST_DELTA)
        precomputed = precompute_data(config.domain, config, TEST_HAM)
        scratch = KrausScratch(ComplexF64, DIM)
        fw = build_trajectoryframework(
            TEST_JUMPS, TEST_HAM, config, precomputed, scratch, TEST_DELTA
        )

        completeness = fw.K0' * fw.K0 + TEST_DELTA * fw.R + fw.U_residual' * fw.U_residual
        identity = Matrix{ComplexF64}(I, DIM, DIM)
        @test isapprox(completeness, identity; atol=1e-10)
    end

    @testset "TimeDomain" begin
        config = make_thermalize_config(TimeDomain(); delta=TEST_DELTA)
        precomputed = precompute_data(config.domain, config, TEST_HAM)
        scratch = KrausScratch(ComplexF64, DIM)
        fw = build_trajectoryframework(
            TEST_JUMPS, TEST_HAM, config, precomputed, scratch, TEST_DELTA
        )

        completeness = fw.K0' * fw.K0 + TEST_DELTA * fw.R + fw.U_residual' * fw.U_residual
        identity = Matrix{ComplexF64}(I, DIM, DIM)
        @test isapprox(completeness, identity; atol=1e-10)
    end

    @testset "TrotterDomain" begin
        config = make_thermalize_config(TrotterDomain(); delta=TEST_DELTA)
        precomputed = precompute_data(config.domain, config, TEST_TROTTER)
        scratch = KrausScratch(ComplexF64, DIM)
        fw = build_trajectoryframework(
            TEST_JUMPS, TEST_TROTTER, config, precomputed, scratch, TEST_DELTA
        )

        completeness = fw.K0' * fw.K0 + TEST_DELTA * fw.R + fw.U_residual' * fw.U_residual
        identity = Matrix{ComplexF64}(I, DIM, DIM)
        @test isapprox(completeness, identity; atol=1e-10)
    end

end
