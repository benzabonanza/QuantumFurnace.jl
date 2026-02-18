using Test
using QuantumFurnace
using LinearAlgebra

# ============================================================================
# Tests for estimate_spectral_gap (Phase 23)
# ============================================================================
#
# Uses the 3-qubit SMALL system from test_helpers.jl for fast execution.

@testset "Gap Estimation" begin

    # Shared setup: 3-qubit initial state and config
    psi0 = zeros(ComplexF64, SMALL_DIM)
    psi0[1] = 1.0
    config = make_small_thermalize_config(TimeDomain();
        delta=0.01, mixing_time=10.0, with_coherent=false)

    # -----------------------------------------------------------------
    @testset "Basic estimate_spectral_gap returns SpectralGapResult" begin
        result = estimate_spectral_gap(
            SMALL_JUMPS, config, psi0, SMALL_HAM;
            ntraj=500, save_every=5, seed=42,
        )

        @test result isa SpectralGapResult
        @test result.gap > 0.0
        @test result.gap_ci[1] < result.gap < result.gap_ci[2]
        @test result.gap_se > 0.0
        @test result.best_observable in ["H", "Mz", "XX_avg", "YY_avg", "ZZ_avg"]
        @test result.best_r_squared > 0.0
        @test length(result.per_observable) == 5
        @test length(result.observable_names) == 5
        @test result.observable_names == ["H", "Mz", "XX_avg", "YY_avg", "ZZ_avg"]
        @test result.ntraj == 500
        @test result.save_every == 5
        @test result.seed == 42
        @test result.skip_initial == 0.0
    end

    # -----------------------------------------------------------------
    @testset "Per-observable results are accessible" begin
        result = estimate_spectral_gap(
            SMALL_JUMPS, config, psi0, SMALL_HAM;
            ntraj=500, save_every=5, seed=42,
        )

        for i in eachindex(result.per_observable)
            @test result.per_observable[i] isa FitResult
            @test result.per_observable[i].gap >= 0.0
        end

        # Best observable's fit matches best_r_squared
        best_idx = findfirst(==(result.best_observable), result.observable_names)
        @test best_idx !== nothing
        @test result.per_observable[best_idx].r_squared == result.best_r_squared
    end

    # -----------------------------------------------------------------
    @testset "Deterministic seeding produces identical results" begin
        r1 = estimate_spectral_gap(
            SMALL_JUMPS, config, psi0, SMALL_HAM;
            ntraj=500, save_every=5, seed=123,
        )
        r2 = estimate_spectral_gap(
            SMALL_JUMPS, config, psi0, SMALL_HAM;
            ntraj=500, save_every=5, seed=123,
        )

        @test r1.gap == r2.gap
        @test r1.best_observable == r2.best_observable
        for i in eachindex(r1.per_observable)
            @test r1.per_observable[i].gap == r2.per_observable[i].gap
        end
    end

    # -----------------------------------------------------------------
    @testset "skip_initial affects fit" begin
        r1 = estimate_spectral_gap(
            SMALL_JUMPS, config, psi0, SMALL_HAM;
            ntraj=500, save_every=5, seed=42, skip_initial=0.0,
        )
        r2 = estimate_spectral_gap(
            SMALL_JUMPS, config, psi0, SMALL_HAM;
            ntraj=500, save_every=5, seed=42, skip_initial=0.2,
        )

        # Gaps should differ because different data windows are fitted
        @test r1.gap != r2.gap
        @test r2.skip_initial == 0.2
    end

    # -----------------------------------------------------------------
    @testset "Custom observables with names" begin
        obs, names = build_preset_trajectory_observables(SMALL_HAM, 3)
        @test length(obs) == 5
        result = estimate_spectral_gap(
            SMALL_JUMPS, config, psi0, SMALL_HAM;
            observables=obs, observable_names=names,
            ntraj=200, save_every=5, seed=42,
        )

        @test result isa SpectralGapResult
        @test length(result.per_observable) == length(obs)
        @test result.observable_names == names
    end

    # -----------------------------------------------------------------
    @testset "ArgumentError on observables without names" begin
        obs, _ = build_preset_trajectory_observables(SMALL_HAM, 3)
        @test_throws ArgumentError estimate_spectral_gap(
            SMALL_JUMPS, config, psi0, SMALL_HAM;
            observables=[obs[1]],
        )
    end

    # -----------------------------------------------------------------
    @testset "Selection logic: _select_best_observable" begin
        # Create FitResult objects with known properties using synthetic data
        # fit_exponential_decay produces real FitResult instances

        # fit_good: converged=true, gap > 0, good R-squared
        times_good = collect(0.0:0.1:10.0)
        values_good = 2.0 .* exp.(-0.5 .* times_good) .+ 0.3
        fit_good = fit_exponential_decay(times_good, values_good)
        @test fit_good.converged == true
        @test fit_good.gap > 0.0
        @test fit_good.r_squared > 0.9

        # For the selection test, directly construct FitResult with controlled properties
        # FitResult is a plain struct, so we can construct it directly
        fit_converged_good = FitResult(
            0.5, 2.0, 0.3,        # gap, amplitude, offset
            (0.4, 0.6), 0.05,     # gap_ci, gap_se
            0.95, true,            # r_squared, converged
            zeros(10), collect(0.0:1.0:9.0), ones(10),  # residuals, times_used, values_used
        )
        fit_unconverged_high_r2 = FitResult(
            0.3, 1.5, 0.2,        # gap, amplitude, offset
            (0.1, 0.5), 0.1,      # gap_ci, gap_se
            0.99, false,           # r_squared=0.99 but NOT converged
            zeros(10), collect(0.0:1.0:9.0), ones(10),
        )
        fit_converged_zero_gap = FitResult(
            0.0, 1.0, 1.0,        # gap=0, amplitude, offset
            (-0.1, 0.1), 0.05,    # gap_ci, gap_se
            0.90, true,            # r_squared, converged
            zeros(10), collect(0.0:1.0:9.0), ones(10),
        )

        # Test selection: should pick fit_converged_good (index 1)
        # because it's the only valid fit (converged AND gap > 0), so trivially
        # the smallest gap among valid fits under the new selection criterion
        idx, name, r2 = QuantumFurnace._select_best_observable(
            [fit_converged_good, fit_unconverged_high_r2, fit_converged_zero_gap],
            ["A", "B", "C"],
        )
        @test idx == 1
        @test name == "A"
        @test r2 == 0.95

        # Test fallback: when no fit is converged with gap > 0, pick highest R-squared
        idx2, name2, r2_2 = QuantumFurnace._select_best_observable(
            [fit_unconverged_high_r2, fit_converged_zero_gap],
            ["B", "C"],
        )
        @test idx2 == 1   # fit_unconverged_high_r2 has R2=0.99 > 0.90
        @test name2 == "B"
        @test r2_2 == 0.99

        # Test smallest-gap selection: when two fits are both valid (converged,
        # gap > 0, R² > 0.8), the one with the smaller gap should be selected
        # (not the one with higher R-squared). The true spectral gap is the
        # slowest-decaying mode, so the smallest positive gap wins.
        fit_high_r2_high_gap = FitResult(
            0.8, 2.0, 0.3,        # gap=0.8, amplitude, offset
            (0.7, 0.9), 0.05,     # gap_ci, gap_se
            0.98, true,            # r_squared=0.98, converged
            zeros(10), collect(0.0:1.0:9.0), ones(10),
        )
        fit_lower_r2_small_gap = FitResult(
            0.3, 1.5, 0.2,        # gap=0.3, amplitude, offset
            (0.2, 0.4), 0.05,     # gap_ci, gap_se
            0.92, true,            # r_squared=0.92, converged
            zeros(10), collect(0.0:1.0:9.0), ones(10),
        )
        idx3, name3, r2_3 = QuantumFurnace._select_best_observable(
            [fit_high_r2_high_gap, fit_lower_r2_small_gap],
            ["high_r2", "small_gap"],
        )
        @test idx3 == 2          # smallest gap wins, not highest R²
        @test name3 == "small_gap"
        @test r2_3 == 0.92
    end

    # -----------------------------------------------------------------
    @testset "Eigenbasis Overlap Analysis" begin

        # Shared setup: construct Liouvillian for the 3-qubit SMALL system
        config_l = make_small_liouv_config(TimeDomain(); with_coherent=false)
        liouv_result = run_lindbladian(SMALL_JUMPS, config_l, SMALL_HAM)
        L = liouv_result.liouvillian

        # Build observables
        obs, names = build_preset_trajectory_observables(SMALL_HAM, 3)

        # Initial state: excited state (far from Gibbs at high beta)
        psi0 = zeros(ComplexF64, SMALL_DIM)
        psi0[end] = 1.0
        rho0 = psi0 * psi0'

        # Call the analysis function
        result = eigenbasis_overlap_analysis(L, obs, names, rho0)

        @testset "Basic return type and structure" begin
            @test result isa OverlapAnalysisResult
            @test length(result.eigenvalues) == SMALL_DIM^2
            @test result.exact_gap > 0.0
            @test length(result.observable_names) == 5
            @test result.observable_names == names
            @test size(result.overlap_coefficients) == (5, SMALL_DIM^2)
            @test length(result.gap_mode_overlap) == 5
            @test length(result.relative_gap_overlap) == 5
            @test all(x -> x >= 0.0, result.gap_mode_overlap)
            @test all(x -> 0.0 <= x <= 1.0, result.relative_gap_overlap)
        end

        @testset "Exact gap matches run_lindbladian" begin
            # ARPACK (shift-invert) and dense eigen agree to ~1e-10;
            # use atol=1e-8 to allow for numerical method differences
            @test isapprox(result.exact_gap, abs(real(liouv_result.spectral_gap)); atol=1e-8)
        end

        @testset "Overlap coefficients reconstruct time series at t=0" begin
            # At t=0: <O>(0) = sum_k c_k = tr(O * rho0)
            c_sum = sum(result.overlap_coefficients[1, :])
            tr_O_rho0 = tr(obs[1] * rho0)
            @test isapprox(real(c_sum), real(tr_O_rho0); atol=1e-10)
        end

        @testset "Steady state mode has zero decay" begin
            @test abs(real(result.eigenvalues[1])) < 1e-10
        end

    end

end
