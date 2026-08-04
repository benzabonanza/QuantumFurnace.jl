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
                config = make_config(Thermalize(),domain; construction=KMS(), delta=TEST_DELTA)
                ws = QuantumFurnace._build_trajectory_workspace(config, TEST_HAM, TEST_JUMPS; delta=TEST_DELTA)

                # Run a single step with a fixed RNG seed
                psi = zeros(ComplexF64, DIM)
                psi[1] = 1.0
                rng = Random.Xoshiro(42)
                step_along_trajectory!(psi, ws, rng)

                # Output must be normalized (U_B ordering bug produced unnormalized output)
                psi_norm = norm(psi)
                @test isapprox(psi_norm, 1.0; atol=1e-12)  # Single-step normalization: unitary + projection preserves norm to machine precision
                @info "TFIX-02: post-step norm ($(typeof(domain)))" norm=psi_norm threshold_atol=1e-12
                # Output must not contain NaN
                @test all(isfinite, psi)
            end
        end
    end

    # TFIX-03: Normalization warning
    @testset "TFIX-03: Normalization check present" begin
        # The normalization check is a @warn that triggers when total_weight deviates from 1.0.
        # With a properly-built workspace and normalized input, total_weight should be ~1.0
        # and NO warning should fire. We verify the code path runs without error.
        config = make_config(Thermalize(),EnergyDomain(); delta=TEST_DELTA)
        ws = QuantumFurnace._build_trajectory_workspace(config, TEST_HAM, TEST_JUMPS; delta=TEST_DELTA)

        psi = zeros(ComplexF64, DIM)
        psi[1] = 1.0
        rng = Random.Xoshiro(123)

        # Should complete without error; with properly built workspace, no warning expected
        step_along_trajectory!(psi, ws, rng)
        psi_norm = norm(psi)
        @test isapprox(psi_norm, 1.0; atol=1e-12)  # Post-step normalization: same as TFIX-02
        @info "TFIX-03: post-step norm" norm=psi_norm threshold_atol=1e-12

        # Verify total_weight is approximately 1.0 by checking per-operator CPTP property
        # (if K0_a'K0_a + delta*R_a + U_res_a'U_res_a = I, then total_weight must be ~1.0)
        max_err = 0.0
        for a in 1:ws.n_jumps
            completeness = ws.K0s[a]' * ws.K0s[a] + ws.delta * ws.Rs[a] + ws.U_residuals[a]' * ws.U_residuals[a]
            err = norm(completeness - Matrix{ComplexF64}(I, DIM, DIM))
            max_err = max(max_err, err)
            @test isapprox(completeness, Matrix{ComplexF64}(I, DIM, DIM); atol=1e-10)  # CPTP algebraic identity (same as TVAL-01)
        end
        @info "TFIX-03: CPTP completeness" n_jumps=ws.n_jumps max_error=max_err threshold_atol=1e-10
    end

    # TFIX-04: PSD guard
    @testset "TFIX-04: PSD guard prevents Cholesky crash" begin
        # Build workspace -- if PSD guard works, no NaN/Inf in U_residual
        for domain in [EnergyDomain(), TimeDomain()]
            @testset "$(typeof(domain))" begin
                config = make_config(Thermalize(),domain; delta=TEST_DELTA)
                ws = QuantumFurnace._build_trajectory_workspace(config, TEST_HAM, TEST_JUMPS; delta=TEST_DELTA)

                # Each per-operator U_residual must be all-finite (no NaN from failed decomposition)
                min_eigenvalue = Inf
                for a in 1:ws.n_jumps
                    @test all(isfinite, ws.U_residuals[a])
                    # U_residual' * U_residual must be PSD (all eigenvalues >= 0)
                    UtU = ws.U_residuals[a]' * ws.U_residuals[a]
                    eigenvalues = eigvals(Hermitian(UtU))
                    min_eigenvalue = min(min_eigenvalue, minimum(eigenvalues))
                    @test all(v -> v >= -1e-14, eigenvalues)  # PSD: eigenvalues >= 0, allow -1e-14 for FP rounding in Hermitian eigvals
                end
                @info "TFIX-04: PSD guard ($(typeof(domain)))" n_jumps=ws.n_jumps min_eigenvalue=min_eigenvalue threshold=-1e-14
            end
        end
    end

    # TFIX-05: Jump sampling faithfulness
    @testset "TFIX-05: Jump sampling matches paper construction" begin
        # Verify per-operator channel structure matches Chen 2023 Theorem III.1:
        # p_nojump + p_jump_total + p_res should sum to ~1.0 for any normalized psi
        config = make_config(Thermalize(),EnergyDomain(); delta=TEST_DELTA)
        ws = QuantumFurnace._build_trajectory_workspace(config, TEST_HAM, TEST_JUMPS; delta=TEST_DELTA)

        # Test with several random initial states
        Random.seed!(999)
        max_prob_err = 0.0
        max_norm_err = 0.0
        for _ in 1:5
            psi = randn(ComplexF64, DIM)
            psi ./= norm(psi)

            # Compute probabilities manually for each per-operator channel
            for a in 1:ws.n_jumps
                K0_psi = ws.K0s[a] * psi
                p_nojump = real(dot(K0_psi, K0_psi))
                p_res = real(dot(ws.U_residuals[a] * psi, ws.U_residuals[a] * psi))
                p_jump_total = ws.delta * real(dot(psi, ws.Rs[a] * psi))

                total = p_nojump + p_res + p_jump_total
                max_prob_err = max(max_prob_err, abs(total - 1.0))
                @test isapprox(total, 1.0; atol=1e-10)  # Probability conservation: algebraic identity from CPTP (same as TVAL-01)
            end

            # Run a single step and verify output is still normalized
            psi_copy = copy(psi)
            ws_check = QuantumFurnace._copy_workspace_for_thread(ws)
            rng_check = Random.Xoshiro(42)
            step_along_trajectory!(psi_copy, ws_check, rng_check)
            norm_err = abs(norm(psi_copy) - 1.0)
            max_norm_err = max(max_norm_err, norm_err)
            @test isapprox(norm(psi_copy), 1.0; atol=1e-12)  # Post-step normalization: unitary channels preserve norm
        end
        @info "TFIX-05: probability conservation" n_states=5 n_jumps=ws.n_jumps max_prob_error=max_prob_err threshold_prob=1e-10 max_norm_error=max_norm_err threshold_norm=1e-12
    end

end
