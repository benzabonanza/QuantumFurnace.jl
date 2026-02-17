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
        @test result.best_observable in ["H", "Mz"]
        @test result.best_r_squared > 0.0
        @test length(result.per_observable) == 2
        @test length(result.observable_names) == 2
        @test result.observable_names == ["H", "Mz"]
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
        @test r1.per_observable[1].gap == r2.per_observable[1].gap
        @test r1.per_observable[2].gap == r2.per_observable[2].gap
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
        obs, names = build_gap_estimation_observables(SMALL_HAM, 3)
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
        obs, _ = build_gap_estimation_observables(SMALL_HAM, 3)
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
        # because it's the only one that is converged AND has gap > 0
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
    end

    # -----------------------------------------------------------------
    # Cross-Validation tests (Phase 24)
    # -----------------------------------------------------------------

    # Helper to create a SpectralGapResult with controlled gap and gap_ci
    function _make_test_spectral_gap_result(; gap=0.5, gap_ci=(0.3, 0.7), gap_se=0.1)
        fit = FitResult(gap, 1.0, 0.0, gap_ci, gap_se, 0.95, true, Float64[], Float64[], Float64[])
        SpectralGapResult(gap, gap_ci, gap_se, "H", 0.95, [fit], ["H"], 1000, 10.0, 10, 42, 0.0)
    end

    @testset "Cross-Validation" begin

        # -----------------------------------------------------------------
        @testset "CrossValidationResult fields are correct for known inputs" begin
            sgr = _make_test_spectral_gap_result(gap=0.48, gap_ci=(0.3, 0.7))
            cv = cross_validate_gap(sgr, -0.5 + 0.1im)

            @test cv.exact_gap == 0.5
            @test cv.fitted_gap == 0.48
            @test cv.absolute_error ≈ 0.02
            @test cv.relative_error ≈ 0.02 / 0.5
            @test cv.within_ci == true  # 0.5 in [0.3, 0.7]
            @test cv.imaginary_ratio ≈ 0.2   # |0.1| / |0.5|
            @test cv.imaginary_warning == true
        end

        # -----------------------------------------------------------------
        @testset "within_ci is true when exact gap falls inside CI" begin
            sgr = _make_test_spectral_gap_result(gap=0.5, gap_ci=(0.3, 0.7))
            cv = cross_validate_gap(sgr, -0.5 + 0.0im)

            @test cv.within_ci == true
        end

        # -----------------------------------------------------------------
        @testset "within_ci is false when exact gap falls outside CI" begin
            sgr = _make_test_spectral_gap_result(gap=0.7, gap_ci=(0.6, 0.8))
            cv = cross_validate_gap(sgr, -0.5 + 0.0im)

            @test cv.within_ci == false  # exact_gap=0.5 is outside [0.6, 0.8]
        end

        # -----------------------------------------------------------------
        @testset "imaginary_warning is false when |Im/Re| <= 0.1" begin
            sgr = _make_test_spectral_gap_result()
            cv = cross_validate_gap(sgr, -0.5 + 0.04im)  # ratio = 0.08

            @test cv.imaginary_ratio ≈ 0.08
            @test cv.imaginary_warning == false
        end

        # -----------------------------------------------------------------
        @testset "imaginary_warning is true when |Im/Re| > 0.1" begin
            sgr = _make_test_spectral_gap_result()
            cv = cross_validate_gap(sgr, -0.5 + 0.1im)  # ratio = 0.2

            @test cv.imaginary_ratio ≈ 0.2
            @test cv.imaginary_warning == true
        end

        # -----------------------------------------------------------------
        @testset "LindbladianResult dispatch extracts spectral_gap" begin
            sgr = _make_test_spectral_gap_result(gap=0.48, gap_ci=(0.3, 0.7))

            # Construct a minimal LindbladianResult with known spectral_gap
            dim = 2
            dummy_L = zeros(ComplexF64, dim^2, dim^2)
            dummy_dm = zeros(ComplexF64, dim, dim)
            dummy_dm[1, 1] = 1.0
            exact = LindbladianResult(
                liouvillian = dummy_L,
                fixed_point = dummy_dm,
                gap_mode = dummy_dm,
                spectral_gap = -0.5 + 0.1im,
            )

            cv = cross_validate_gap(sgr, exact)
            @test cv.exact_gap == 0.5
            @test cv.fitted_gap == 0.48
            @test cv.imaginary_warning == true
        end

        # -----------------------------------------------------------------
        @testset "@warn is emitted for significant imaginary part" begin
            sgr = _make_test_spectral_gap_result()

            @test_logs (:warn, r"significant imaginary part") begin
                cross_validate_gap(sgr, -0.5 + 0.1im)
            end
        end

        # -----------------------------------------------------------------
        @testset "No warning for small imaginary part" begin
            sgr = _make_test_spectral_gap_result()

            @test_logs cross_validate_gap(sgr, -0.5 + 0.04im)
        end

        # -----------------------------------------------------------------
        @testset "Edge case: Re(eigenvalue) == 0 gives Inf" begin
            sgr = _make_test_spectral_gap_result()
            cv = cross_validate_gap(sgr, 0.0 + 0.5im)

            @test cv.exact_gap == 0.0
            @test cv.relative_error == Inf
            @test cv.imaginary_ratio == Inf
            @test cv.imaginary_warning == true
        end

    end  # Cross-Validation

end
