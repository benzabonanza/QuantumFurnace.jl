using Test
using LinearAlgebra
using Random
using QuantumFurnace

# test_helpers.jl is already included by runtests.jl

# ============================================================================
# Comprehensive tests for krylov_eigsolve.jl: spectral gap, channel path,
# domain coverage, guard rails, and eigenvalue properties.
# Phase 29: Eigensolver Integration
# Quick-36: Faithful Chen channel for apply_delta_channel!
# ============================================================================

@testset "Krylov Eigsolve" begin

    # ========================================================================
    # Testset 1: apply_delta_channel! faithful Chen channel
    # ========================================================================
    @testset "apply_delta_channel! faithful Chen channel" begin
        config_therm = make_thermalize_config(EnergyDomain(); with_coherent=true, delta=0.01)
        config_liouv = make_liouv_config(EnergyDomain(); with_coherent=true)
        ws = KrylovWorkspace(config_therm, TEST_HAM, TEST_JUMPS)
        delta = config_therm.delta

        # Dense Lindbladian for Euler comparison
        L_dense = construct_lindbladian(TEST_JUMPS, config_liouv, TEST_HAM)
        I_d2 = Matrix{ComplexF64}(LinearAlgebra.I(DIM^2))

        for _ in 1:5
            rho = Matrix(random_density_matrix(NUM_QUBITS))

            # Faithful Chen channel
            apply_delta_channel!(ws, rho, config_liouv, TEST_HAM)
            rho_chen = copy(ws.rho_out)

            # Trace preservation: tr(E(rho)) == tr(rho)
            @test isapprox(real(tr(rho_chen)), real(tr(rho)); atol=1e-10)

            # Positivity: eigenvalues of E(rho) >= -eps for valid density matrix input
            eigs = eigvals(Hermitian(rho_chen))
            @test all(eigs .> -1e-10)

            # O(delta^2) close to Euler: |E_chen(rho) - E_euler(rho)| < C * delta^2
            v_euler = (I_d2 + delta * L_dense) * vec(rho)
            @test norm(vec(rho_chen) - v_euler) < 50 * delta^2
        end
    end

    # ========================================================================
    # Testset 2: KrylovGapResult struct fields
    # ========================================================================
    @testset "KrylovGapResult struct fields" begin
        config_kms = make_liouv_config(EnergyDomain(); with_coherent=true)
        result = krylov_spectral_gap(config_kms, TEST_HAM, TEST_JUMPS;
            krylovdim=30, howmany=4)

        @test result isa KrylovGapResult
        @test length(result.eigenvalues) >= 2
        @test result.spectral_gap > 0
        @test size(result.fixed_point) == (DIM, DIM)
        @test size(result.gap_mode) == (DIM, DIM)
        @test result.converged >= 2
        @test result.matvec_count > 0
        @test result.channel_eigenvalues === nothing
        @test result.delta_used === nothing
    end

    # ========================================================================
    # Testset 3: Lindbladian eigsolve accuracy (EnergyDomain KMS)
    # ========================================================================
    @testset "Lindbladian eigsolve accuracy (EnergyDomain KMS)" begin
        config = make_liouv_config(EnergyDomain(); with_coherent=true)
        L_dense = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)
        dense_result = extract_leading_eigendata(L_dense; n_modes=4)

        result = krylov_spectral_gap(config, TEST_HAM, TEST_JUMPS;
            krylovdim=30, howmany=4)

        # Spectral gap matches dense reference
        @test isapprox(result.spectral_gap, dense_result.spectral_gap; rtol=1e-6)

        # Fixed point is close to Gibbs state
        @test trace_distance_h(Hermitian(result.fixed_point), TEST_GIBBS) < 1e-4

        # Steady-state eigenvalue near zero
        @test abs(real(result.eigenvalues[1])) < 1e-8
    end

    # ========================================================================
    # Testset 4: Lindbladian eigsolve accuracy (EnergyDomain GNS)
    # ========================================================================
    @testset "Lindbladian eigsolve accuracy (EnergyDomain GNS)" begin
        config = make_liouv_config_gns(EnergyDomain())
        L_dense = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)
        dense_result = extract_leading_eigendata(L_dense; n_modes=4)

        result = krylov_spectral_gap(config, TEST_HAM, TEST_JUMPS;
            krylovdim=30, howmany=4)

        # Spectral gap matches dense reference
        @test isapprox(result.spectral_gap, dense_result.spectral_gap; rtol=1e-6)

        # GNS fixed point differs from exact Gibbs -- skip trace distance check
        # Just verify the fixed point is a valid density matrix
        @test isapprox(tr(result.fixed_point), 1.0; atol=1e-6)
        @test abs(real(result.eigenvalues[1])) < 1e-8
    end

    # ========================================================================
    # Testset 5: Channel eigsolve accuracy (ThermalizeConfig)
    # ========================================================================
    @testset "Channel eigsolve accuracy (ThermalizeConfig)" begin
        config_therm = make_thermalize_config(EnergyDomain();
            with_coherent=true, delta=0.01)
        # Dense reference from Lindbladian path
        config_liouv = make_liouv_config(EnergyDomain(); with_coherent=true)
        L_dense = construct_lindbladian(TEST_JUMPS, config_liouv, TEST_HAM)
        dense_result = extract_leading_eigendata(L_dense; n_modes=4)

        result = krylov_spectral_gap(config_therm, TEST_HAM, TEST_JUMPS;
            krylovdim=30, howmany=4)

        # Channel path has O(delta^2) error from eigenvalue conversion formula
        # lambda_L = (mu-1)/delta is exact for linear channel but approximate for
        # faithful Chen channel whose eigenvalues differ from 1+delta*lambda by O(delta^2)
        @test isapprox(result.spectral_gap, dense_result.spectral_gap; rtol=2e-3)

        # Channel-specific fields are populated
        @test result.channel_eigenvalues !== nothing
        @test result.delta_used == 0.01

        # Channel eigenvalue near 1 for steady state
        @test isapprox(abs(result.channel_eigenvalues[1]), 1.0; atol=0.01)
    end

    # ========================================================================
    # Testset 6: All domains work (Lindbladian path)
    # ========================================================================
    @testset "All domains work (Lindbladian path)" begin
        @testset "EnergyDomain" begin
            config = make_liouv_config(EnergyDomain(); with_coherent=true)
            result = krylov_spectral_gap(config, TEST_HAM, TEST_JUMPS;
                krylovdim=30, howmany=2)
            @test result.spectral_gap > 0
        end

        @testset "TimeDomain" begin
            config = make_liouv_config(TimeDomain(); with_coherent=true)
            result = krylov_spectral_gap(config, TEST_HAM, TEST_JUMPS;
                krylovdim=30, howmany=2)
            @test result.spectral_gap > 0
        end

        @testset "TrotterDomain" begin
            config = make_liouv_config(TrotterDomain(); with_coherent=true)
            result = krylov_spectral_gap(config, TEST_HAM, TEST_TROTTER_JUMPS;
                trotter=TEST_TROTTER, krylovdim=30, howmany=2)
            @test result.spectral_gap > 0
        end

        @testset "BohrDomain" begin
            config = make_liouv_config(BohrDomain(); with_coherent=true)
            result = krylov_spectral_gap(config, TEST_HAM, TEST_JUMPS;
                krylovdim=30, howmany=2)
            @test result.spectral_gap > 0
        end
    end

    # ========================================================================
    # Testset 7: Guard rails
    # ========================================================================
    @testset "Guard rails" begin
        config = make_liouv_config(EnergyDomain(); with_coherent=true)

        # krylovdim <= howmany must error
        @test_throws ErrorException krylov_spectral_gap(config, TEST_HAM, TEST_JUMPS;
            krylovdim=2, howmany=4)

        # Memory guard should not throw for small n=4 system
        result = krylov_spectral_gap(config, TEST_HAM, TEST_JUMPS;
            krylovdim=50, howmany=2)
        @test result isa KrylovGapResult
    end

    # ========================================================================
    # Testset 8: Eigenvalue sorting and conversion
    # ========================================================================
    @testset "Eigenvalue sorting and conversion" begin
        # Lindbladian path: eigenvalues sorted by |Re(lambda)| ascending
        config = make_liouv_config(EnergyDomain(); with_coherent=true)
        result_liouv = krylov_spectral_gap(config, TEST_HAM, TEST_JUMPS;
            krylovdim=30, howmany=4)

        @test abs(real(result_liouv.eigenvalues[1])) <= abs(real(result_liouv.eigenvalues[2]))

        # Channel path: verify eigenvalue conversion formula
        config_therm = make_thermalize_config(EnergyDomain();
            with_coherent=true, delta=0.01)
        result_chan = krylov_spectral_gap(config_therm, TEST_HAM, TEST_JUMPS;
            krylovdim=30, howmany=4)

        # Converted Lindbladian eigenvalues should satisfy: lambda_L = (mu - 1) / delta
        # The result.eigenvalues are already the converted values.
        # Verify by back-computing: mu = 1 + delta * lambda_L
        reconstructed_mu = 1.0 .+ result_chan.delta_used .* result_chan.eigenvalues
        @test isapprox(reconstructed_mu, result_chan.channel_eigenvalues; atol=1e-10)
    end

end  # @testset "Krylov Eigsolve"
