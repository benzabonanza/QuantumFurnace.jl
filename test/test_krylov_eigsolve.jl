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
        config_therm = make_config(Thermalize(),EnergyDomain(); construction=KMS(), delta=0.01)
        config_liouv = make_config(Lindbladian(),EnergyDomain(); construction=KMS())
        ws = Workspace(config_therm, TEST_HAM, TEST_JUMPS)
        delta = config_therm.delta

        # Dense Lindbladian for Euler comparison
        L_dense = construct_lindbladian(TEST_JUMPS, config_liouv, TEST_HAM)
        I_d2 = Matrix{ComplexF64}(LinearAlgebra.I(DIM^2))

        # Threshold rationale (trace preservation): tr(E(rho)) should equal tr(rho) exactly
        # for a CPTP map. Real FP drift from the per-jump Lie–Trotter sweep is at the
        # rounding floor (~5e-16 measured); `hermitianize!` averages off-diagonals
        # without touching the diagonal, so it adds no trace drift. Threshold 1e-10
        # gives >1e5× margin.
        #
        # Threshold rationale (positivity): eigenvalues of E(rho) should be >= 0 for CPTP.
        # FP rounding can produce tiny negatives: O(eps * ||rho||) ~ 1e-16. Threshold -1e-10 is generous.
        #
        # Threshold rationale (Euler closeness): faithful per-jump Lie–Trotter Φ_δ differs
        # from Euler by O(delta^2). The constant absorbs ‖𝓛‖·DIM and the Lie–Trotter
        # splitting constant on the n_jumps substeps. Empirically ~1.6e-6 at δ=1e-2,
        # so the 50·δ² = 5e-3 threshold has >3000x margin.
        max_trace_err = 0.0
        min_eigenvalue = Inf
        max_euler_err = 0.0
        euler_threshold = 50 * delta^2
        for _ in 1:5
            rho = Matrix(random_density_matrix(NUM_QUBITS))

            # Faithful Chen channel
            apply_delta_channel!(ws, rho, config_therm, TEST_HAM)
            rho_chen = copy(ws.scratch.rho_next)

            # Trace preservation: tr(E(rho)) == tr(rho)
            trace_err = abs(real(tr(rho_chen)) - real(tr(rho)))
            @test isapprox(real(tr(rho_chen)), real(tr(rho)); atol=1e-10)
            max_trace_err = max(max_trace_err, trace_err)

            # Positivity: eigenvalues of E(rho) >= -eps for valid density matrix input
            eigs = eigvals(Hermitian(rho_chen))
            @test all(eigs .> -1e-10)
            min_eigenvalue = min(min_eigenvalue, minimum(eigs))

            # O(delta^2) close to Euler: |E_chen(rho) - E_euler(rho)| < C * delta^2
            v_euler = (I_d2 + delta * L_dense) * vec(rho)
            euler_err = norm(vec(rho_chen) - v_euler)
            @test euler_err < euler_threshold
            max_euler_err = max(max_euler_err, euler_err)
        end
        @info "Chen channel trace preservation" max_trace_error=max_trace_err threshold=1e-8
        @info "Chen channel positivity" min_eigenvalue=min_eigenvalue threshold=-1e-10
        @info "Chen channel Euler closeness" max_euler_error=max_euler_err threshold=euler_threshold delta=delta
    end

    # ========================================================================
    # Testset 2: krylov_spectral_gap result fields
    # ========================================================================
    @testset "krylov_spectral_gap result fields" begin
        config_kms = make_config(Lindbladian(),EnergyDomain(); construction=KMS())
        result = krylov_spectral_gap(config_kms, TEST_HAM, TEST_JUMPS;
            krylovdim=30, howmany=4)

        @test result isa NamedTuple
        @test length(result.eigenvalues) >= 2
        @test result.spectral_gap > 0
        @test size(result.fixed_point) == (DIM, DIM)
        @test size(result.gap_mode) == (DIM, DIM)
        @test result.converged >= 2
        @test result.matvec_count > 0
        @test result.channel_eigenvalues === nothing
        @test result.delta_used === nothing
        @info "krylov_spectral_gap result fields" spectral_gap=result.spectral_gap n_eigenvalues=length(result.eigenvalues) converged=result.converged matvec_count=result.matvec_count
    end

    # ========================================================================
    # Testset 3: Lindbladian eigsolve accuracy (EnergyDomain KMS)
    # ========================================================================
    @testset "Lindbladian eigsolve accuracy (EnergyDomain KMS)" begin
        config = make_config(Lindbladian(),EnergyDomain(); construction=KMS())
        L_dense = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)
        dense_result = extract_leading_eigendata(L_dense; n_modes=4)

        result = krylov_spectral_gap(config, TEST_HAM, TEST_JUMPS;
            krylovdim=30, howmany=4)

        # Threshold rationale (gap rtol=1e-6): KrylovKit tol=1e-10 (default) bounds eigenvalue error.
        # rtol=1e-6 gives 10000x margin for iterative convergence variability.
        gap_err = abs(result.spectral_gap - dense_result.spectral_gap) / dense_result.spectral_gap
        @test isapprox(result.spectral_gap, dense_result.spectral_gap; rtol=1e-6)
        @info "Eigsolve gap accuracy (EnergyDomain KMS)" krylov_gap=result.spectral_gap dense_gap=dense_result.spectral_gap relative_error=gap_err rtol=1e-6

        # Threshold rationale (trace distance < 1e-4): KMS fixed point should converge to Gibbs.
        # KrylovKit eigenvector accuracy bounded by tol * condition_number. 1e-4 is conservative.
        td = trace_distance_h(Hermitian(result.fixed_point), TEST_GIBBS)
        @test td < 1e-4
        @info "Fixed point trace distance to Gibbs" trace_distance=td threshold=1e-4

        # Threshold rationale (|Re(lambda_1)| < 1e-8): steady-state eigenvalue is exactly 0.
        # KrylovKit residual bounded by tol. 1e-8 gives 100x margin over tol=1e-10.
        ss_err = abs(real(result.eigenvalues[1]))
        @test ss_err < 1e-8
        @info "Steady-state eigenvalue magnitude" abs_real_lambda1=ss_err threshold=1e-8
    end

    # ========================================================================
    # Testset 4: Lindbladian eigsolve accuracy (EnergyDomain GNS)
    # ========================================================================
    @testset "Lindbladian eigsolve accuracy (EnergyDomain GNS)" begin
        config = make_config(Lindbladian(), EnergyDomain(); construction=GNS())
        L_dense = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)
        dense_result = extract_leading_eigendata(L_dense; n_modes=4)

        result = krylov_spectral_gap(config, TEST_HAM, TEST_JUMPS;
            krylovdim=30, howmany=4)

        # Threshold rationale: same as KMS testset -- KrylovKit tol=1e-10, rtol=1e-6 gives 10000x margin.
        gap_err = abs(result.spectral_gap - dense_result.spectral_gap) / dense_result.spectral_gap
        @test isapprox(result.spectral_gap, dense_result.spectral_gap; rtol=1e-6)
        @info "Eigsolve gap accuracy (EnergyDomain GNS)" krylov_gap=result.spectral_gap dense_gap=dense_result.spectral_gap relative_error=gap_err rtol=1e-6

        # GNS fixed point differs from exact Gibbs -- skip trace distance check
        # Threshold rationale (trace=1 atol=1e-6): valid density matrix has tr(rho)=1.
        # KrylovKit eigenvector normalization should preserve this. 1e-6 is conservative.
        trace_err = abs(tr(result.fixed_point) - 1.0)
        @test isapprox(tr(result.fixed_point), 1.0; atol=1e-6)
        @info "GNS fixed point trace" trace=real(tr(result.fixed_point)) trace_error=trace_err atol=1e-6

        # Threshold rationale: same as KMS -- steady-state eigenvalue is exactly 0.
        ss_err = abs(real(result.eigenvalues[1]))
        @test ss_err < 1e-8
        @info "Steady-state eigenvalue magnitude (GNS)" abs_real_lambda1=ss_err threshold=1e-8
    end

    # ========================================================================
    # Testset 5: Channel eigsolve accuracy (Thermalize)
    # ========================================================================
    @testset "Channel eigsolve accuracy (Thermalize)" begin
        config_therm = make_config(Thermalize(),EnergyDomain();
            construction=KMS(), delta=0.01)
        # Dense reference from Lindbladian path
        config_liouv = make_config(Lindbladian(),EnergyDomain(); construction=KMS())
        L_dense = construct_lindbladian(TEST_JUMPS, config_liouv, TEST_HAM)
        dense_result = extract_leading_eigendata(L_dense; n_modes=4)

        result = krylov_spectral_gap(config_therm, TEST_HAM, TEST_JUMPS;
            krylovdim=30, howmany=4)

        # Threshold rationale (gap rtol=3e-3): channel path has O(delta^2) error from eigenvalue
        # conversion lambda_L = (mu-1)/delta. For delta=0.01, the absolute gap error is ~3e-4.
        # rtol=3e-3 gives ~10x margin at the n=4 KMS smooth-Metro fixture (post-qf-etx the gap
        # is its true continuum value ≈ 0.127; pre-qf-etx the gap was artificially inflated by
        # the grid-dependent 1/gamma_norm_factor sample-sup, which cancelled some of the relative
        # error). Looser than Lindbladian path due to channel approximation.
        gap_err = abs(result.spectral_gap - dense_result.spectral_gap) / dense_result.spectral_gap
        @test isapprox(result.spectral_gap, dense_result.spectral_gap; rtol=3e-3)
        @info "Channel eigsolve gap accuracy" krylov_gap=result.spectral_gap dense_gap=dense_result.spectral_gap relative_error=gap_err rtol=3e-3

        # Channel-specific fields are populated
        @test result.channel_eigenvalues !== nothing
        @test result.delta_used == 0.01

        # Threshold rationale (channel eigenvalue ~1, atol=0.01): steady-state channel eigenvalue
        # is exactly 1.0 for a CPTP map. KrylovKit convergence + delta discretization give O(delta) error.
        chan_ss_err = abs(abs(result.channel_eigenvalues[1]) - 1.0)
        @test isapprox(abs(result.channel_eigenvalues[1]), 1.0; atol=0.01)
        @info "Channel steady-state eigenvalue" abs_mu1=abs(result.channel_eigenvalues[1]) error=chan_ss_err atol=0.01
    end

    # ========================================================================
    # Testset 6: All domains work (Lindbladian path)
    # ========================================================================
    @testset "All domains work (Lindbladian path)" begin
        @testset "EnergyDomain" begin
            config = make_config(Lindbladian(),EnergyDomain(); construction=KMS())
            result = krylov_spectral_gap(config, TEST_HAM, TEST_JUMPS;
                krylovdim=30, howmany=2)
            @test result.spectral_gap > 0
            @info "Domain coverage (EnergyDomain)" spectral_gap=result.spectral_gap
        end

        @testset "TimeDomain" begin
            config = make_config(Lindbladian(),TimeDomain(); construction=KMS())
            result = krylov_spectral_gap(config, TEST_HAM, TEST_JUMPS;
                krylovdim=30, howmany=2)
            @test result.spectral_gap > 0
            @info "Domain coverage (TimeDomain)" spectral_gap=result.spectral_gap
        end

        @testset "TrotterDomain" begin
            config = make_config(Lindbladian(),TrotterDomain(); construction=KMS())
            result = krylov_spectral_gap(config, TEST_HAM, TEST_TROTTER_JUMPS;
                trotter=TEST_TROTTER, krylovdim=30, howmany=2)
            @test result.spectral_gap > 0
            @info "Domain coverage (TrotterDomain)" spectral_gap=result.spectral_gap
        end

        @testset "BohrDomain" begin
            config = make_config(Lindbladian(),BohrDomain(); construction=KMS())
            result = krylov_spectral_gap(config, TEST_HAM, TEST_JUMPS;
                krylovdim=30, howmany=2)
            @test result.spectral_gap > 0
            @info "Domain coverage (BohrDomain)" spectral_gap=result.spectral_gap
        end
    end

    # ========================================================================
    # Testset 7: Guard rails
    # ========================================================================
    @testset "Guard rails" begin
        config = make_config(Lindbladian(),EnergyDomain(); construction=KMS())

        # krylovdim <= howmany must error
        @test_throws ErrorException krylov_spectral_gap(config, TEST_HAM, TEST_JUMPS;
            krylovdim=2, howmany=4)

        # Memory guard should not throw for small n=4 system
        result = krylov_spectral_gap(config, TEST_HAM, TEST_JUMPS;
            krylovdim=50, howmany=2)
        @test result isa NamedTuple
    end

    # ========================================================================
    # Testset 8: Eigenvalue sorting and conversion
    # ========================================================================
    @testset "Eigenvalue sorting and conversion" begin
        # Lindbladian path: eigenvalues sorted by |Re(lambda)| ascending (structural check)
        config = make_config(Lindbladian(),EnergyDomain(); construction=KMS())
        result_liouv = krylov_spectral_gap(config, TEST_HAM, TEST_JUMPS;
            krylovdim=30, howmany=4)

        @test abs(real(result_liouv.eigenvalues[1])) <= abs(real(result_liouv.eigenvalues[2]))
        @info "Eigenvalue sorting" abs_Re_lambda1=abs(real(result_liouv.eigenvalues[1])) abs_Re_lambda2=abs(real(result_liouv.eigenvalues[2]))

        # Channel path: verify eigenvalue conversion formula
        config_therm = make_config(Thermalize(),EnergyDomain();
            construction=KMS(), delta=0.01)
        result_chan = krylov_spectral_gap(config_therm, TEST_HAM, TEST_JUMPS;
            krylovdim=30, howmany=4)

        # Threshold rationale (atol=1e-10): conversion mu = 1 + delta * lambda_L is algebraically
        # exact (just arithmetic). Error is FP rounding only: O(eps * |mu|) ~ 1e-16.
        # Threshold 1e-10 gives >1e6 margin.
        reconstructed_mu = 1.0 .+ result_chan.delta_used .* result_chan.eigenvalues
        conversion_err = maximum(abs.(reconstructed_mu .- result_chan.channel_eigenvalues))
        @test isapprox(reconstructed_mu, result_chan.channel_eigenvalues; atol=1e-10)
        @info "Eigenvalue conversion consistency" max_conversion_error=conversion_err atol=1e-10
    end

end  # @testset "Krylov Eigsolve"
