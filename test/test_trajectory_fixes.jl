using Test
using LinearAlgebra
using Random

@testset "Trajectory Bug Fixes" begin

    # TFIX-02: U_B ordering
    @testset "TFIX-02: U_B applied before branch selection" begin
        # Test with coherent term enabled
        # EnergyDomain had the bug (U_B was applied after probability computation);
        # TimeDomain already had correct ordering but we test both to ensure consistency.
        for domain in [EnergyDomain(), TimeDomain()]
            @testset "$(typeof(domain))" begin
                config = make_thermalize_config(domain; with_coherent=true, delta=TEST_DELTA)
                ham_or_trott = TEST_HAM
                precomputed = precompute_data(config.domain, config, ham_or_trott)
                scratch = KrausScratch(ComplexF64, DIM)
                fw = build_trajectoryframework(
                    TEST_JUMPS, ham_or_trott, config, precomputed, scratch, TEST_DELTA
                )

                # Run a single step with a fixed RNG seed
                psi = zeros(ComplexF64, DIM)
                psi[1] = 1.0
                Random.seed!(42)
                step_along_trajectory!(psi, fw)

                # Output must be normalized (U_B ordering bug produced unnormalized output)
                @test isapprox(norm(psi), 1.0; atol=1e-12)
                # Output must not contain NaN
                @test all(isfinite, psi)
            end
        end
    end

    # TFIX-03: Normalization warning
    @testset "TFIX-03: Normalization check present" begin
        # The normalization check is a @warn that triggers when total_weight deviates from 1.0.
        # With a properly-built framework and normalized input, total_weight should be ~1.0
        # and NO warning should fire. We verify the code path runs without error.
        config = make_thermalize_config(EnergyDomain(); delta=TEST_DELTA)
        precomputed = precompute_data(config.domain, config, TEST_HAM)
        scratch = KrausScratch(ComplexF64, DIM)
        fw = build_trajectoryframework(
            TEST_JUMPS, TEST_HAM, config, precomputed, scratch, TEST_DELTA
        )

        psi = zeros(ComplexF64, DIM)
        psi[1] = 1.0
        Random.seed!(123)

        # Should complete without error; with properly built framework, no warning expected
        step_along_trajectory!(psi, fw)
        @test isapprox(norm(psi), 1.0; atol=1e-12)

        # Verify total_weight is approximately 1.0 by checking framework CPTP property
        # (if K0'K0 + delta*R + U_res'U_res = I, then total_weight must be ~1.0 for normalized psi)
        completeness = fw.K0' * fw.K0 + TEST_DELTA * fw.R + fw.U_residual' * fw.U_residual
        @test isapprox(completeness, Matrix{ComplexF64}(I, DIM, DIM); atol=1e-10)
    end

    # TFIX-04: PSD guard
    @testset "TFIX-04: PSD guard prevents Cholesky crash" begin
        # Build framework -- if PSD guard works, no NaN/Inf in U_residual
        for domain in [EnergyDomain(), TimeDomain()]
            @testset "$(typeof(domain))" begin
                config = make_thermalize_config(domain; delta=TEST_DELTA)
                ham_or_trott = TEST_HAM
                precomputed = precompute_data(config.domain, config, ham_or_trott)
                scratch = KrausScratch(ComplexF64, DIM)
                fw = build_trajectoryframework(
                    TEST_JUMPS, ham_or_trott, config, precomputed, scratch, TEST_DELTA
                )

                # U_residual must be all-finite (no NaN from failed Cholesky)
                @test all(isfinite, fw.U_residual)
                # U_residual' * U_residual must be PSD (all eigenvalues >= 0)
                UtU = fw.U_residual' * fw.U_residual
                eigenvalues = eigvals(Hermitian(UtU))
                @test all(v -> v >= -1e-14, eigenvalues)
            end
        end
    end

    # TFIX-05: Jump sampling faithfulness
    @testset "TFIX-05: Jump sampling matches paper construction" begin
        # Verify channel structure matches Chen 2023 Theorem III.1:
        # p_nojump + p_jump_total + p_res should sum to ~1.0 for any normalized psi
        # This verifies the trajectory code implements the correct probability decomposition
        config = make_thermalize_config(EnergyDomain(); delta=TEST_DELTA)
        precomputed = precompute_data(config.domain, config, TEST_HAM)
        scratch = KrausScratch(ComplexF64, DIM)
        fw = build_trajectoryframework(
            TEST_JUMPS, TEST_HAM, config, precomputed, scratch, TEST_DELTA
        )

        # Test with several random initial states
        Random.seed!(999)
        for _ in 1:5
            psi = randn(ComplexF64, DIM)
            psi ./= norm(psi)

            # Compute probabilities manually from the Kraus operators
            K0_psi = fw.K0 * psi
            p_nojump = real(dot(K0_psi, K0_psi))
            p_res = real(dot(fw.U_residual * psi, fw.U_residual * psi))
            p_jump_total = TEST_DELTA * real(dot(psi, fw.R * psi))

            total = p_nojump + p_res + p_jump_total
            @test isapprox(total, 1.0; atol=1e-10)

            # Run a single step and verify output is still normalized
            psi_copy = copy(psi)
            Random.seed!(42)
            step_along_trajectory!(psi_copy, fw)
            @test isapprox(norm(psi_copy), 1.0; atol=1e-12)
        end
    end

end
