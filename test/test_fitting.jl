using StableRNGs

@testset "Fitting" begin

    # -----------------------------------------------------------------------
    # FIT-01: basic exponential decay recovery (clean data)
    # -----------------------------------------------------------------------
    @testset "FIT-01: basic exponential decay recovery" begin
        A_true = 2.0
        gap_true = 0.5
        C_true = 0.3
        times = collect(0.0:0.1:20.0)
        values = A_true .* exp.(-gap_true .* times) .+ C_true

        result = fit_exponential_decay(times, values)

        @test result isa FitResult
        @test result.converged == true
        @test isapprox(result.gap, gap_true; atol=1e-6)
        @test isapprox(result.amplitude, A_true; atol=1e-6)
        @test isapprox(result.offset, C_true; atol=1e-6)
        @test result.r_squared > 0.999
    end

    # -----------------------------------------------------------------------
    # FIT-01: noisy data recovery within CI
    # -----------------------------------------------------------------------
    @testset "FIT-01: noisy data recovery within CI" begin
        rng = StableRNG(42)
        A_true = 2.0
        gap_true = 0.5
        C_true = 0.3
        times = collect(0.0:0.1:20.0)
        values = A_true .* exp.(-gap_true .* times) .+ C_true .+ 0.05 .* randn(rng, length(times))

        result = fit_exponential_decay(times, values)

        @test result.converged == true
        @test result.gap_ci[1] <= gap_true <= result.gap_ci[2]
        @test result.gap_se > 0.0
        @test result.r_squared > 0.9
        @test length(result.residuals) == length(times)
        @test result.times_used == times
        @test result.values_used == values
    end

    # -----------------------------------------------------------------------
    # FIT-02: auto-generated initial guess (log-linear)
    # -----------------------------------------------------------------------
    @testset "FIT-02: auto-generated initial guess (log-linear)" begin
        A_true = 5.0
        gap_true = 1.5
        C_true = -0.5
        times = collect(0.0:0.05:5.0)
        values = A_true .* exp.(-gap_true .* times) .+ C_true

        result = fit_exponential_decay(times, values)

        @test result.converged == true
        @test isapprox(result.gap, gap_true; atol=0.01)

        # Test internal helper directly
        p0 = QuantumFurnace._log_linear_initial_guess(times, values)
        @test length(p0) == 3
        @test p0[2] > 0.0  # gap guess positive
    end

    # -----------------------------------------------------------------------
    # FIT-02: log-linear fallback for difficult data
    # -----------------------------------------------------------------------
    @testset "FIT-02: log-linear fallback for difficult data" begin
        times = collect(0.0:0.1:4.9)
        values = fill(1.0, 50)

        p0 = QuantumFurnace._log_linear_initial_guess(times, values)
        @test length(p0) == 3
        @test all(isfinite, p0)
    end

    # -----------------------------------------------------------------------
    # FIT-03: skip_initial window selection
    # -----------------------------------------------------------------------
    @testset "FIT-03: skip_initial window selection" begin
        gap_true = 0.3
        times = collect(0.0:0.05:10.0)
        # Data with fast-decaying transient added to the slow exponential
        values = 2.0 .* exp.(-gap_true .* times) .+ 0.5 .+ 3.0 .* exp.(-2.0 .* times)

        r1 = fit_exponential_decay(times, values)
        r2 = fit_exponential_decay(times, values; skip_initial=0.3)

        @test abs(r2.gap - gap_true) < abs(r1.gap - gap_true)
        @test length(r2.times_used) < length(r1.times_used)
    end

    # -----------------------------------------------------------------------
    # FIT-04: quality metrics present and correct
    # -----------------------------------------------------------------------
    @testset "FIT-04: quality metrics present and correct" begin
        times = collect(0.0:0.1:10.0)
        values = 2.0 .* exp.(-0.5 .* times) .+ 0.3

        result = fit_exponential_decay(times, values)

        @test result.r_squared isa Float64
        @test result.gap_ci isa Tuple{Float64, Float64}
        @test result.gap_ci[1] < result.gap_ci[2]
        @test result.gap_se isa Float64
        @test result.gap_se >= 0.0
        @test result.converged isa Bool
        @test result.residuals isa Vector{Float64}
    end

    # -----------------------------------------------------------------------
    # FIT-05: gap > 0 enforced via bounds
    # -----------------------------------------------------------------------
    @testset "FIT-05: gap > 0 enforced via bounds" begin
        times = collect(0.0:0.1:10.0)
        # Negative amplitude: rising exponential that LM might explore negative gap for
        values = -1.5 .* exp.(-0.3 .* times) .+ 2.0

        result = fit_exponential_decay(times, values)

        @test result.gap >= 0.0
    end

    # -----------------------------------------------------------------------
    # FIT-04: R-squared not clamped for bad fits
    # -----------------------------------------------------------------------
    @testset "FIT-04: R-squared not clamped for bad fits" begin
        times = collect(0.0:0.1:10.0)
        values = sin.(times)

        result = fit_exponential_decay(times, values)

        @test result.r_squared < 0.5
    end

    # -----------------------------------------------------------------------
    # Custom p0 override
    # -----------------------------------------------------------------------
    @testset "custom p0 override" begin
        A_true = 2.0
        gap_true = 0.5
        C_true = 0.3
        times = collect(0.0:0.1:20.0)
        values = A_true .* exp.(-gap_true .* times) .+ C_true

        result = fit_exponential_decay(times, values; p0=[2.0, 0.5, 0.3])

        @test result.converged == true
        @test isapprox(result.gap, gap_true; atol=1e-6)
    end

    # ===================================================================
    # Bi-exponential fitting tests
    # ===================================================================

    # -----------------------------------------------------------------------
    # BIEXP-01: Clean bi-exponential data recovery
    # -----------------------------------------------------------------------
    @testset "BIEXP-01: clean bi-exp data recovery" begin
        A1_true = 1.0    # fast amplitude
        g1_true = 2.0    # fast gap
        A2_true = 0.5    # slow amplitude
        g2_true = 0.3    # slow gap (spectral gap estimate)
        C_true  = 0.001  # offset

        times = collect(0.0:0.05:30.0)
        values = A1_true .* exp.(-g1_true .* times) .+
                 A2_true .* exp.(-g2_true .* times) .+
                 C_true

        result = fit_biexponential_decay(times, values)

        @test result isa BiexpFitResult
        @test result.converged == true
        @test result.r_squared > 0.999

        # Slow mode (spectral gap)
        @test isapprox(result.gap, g2_true; rtol=0.05)
        @test isapprox(result.amplitude, A2_true; rtol=0.1)

        # Fast mode
        @test isapprox(result.gap_fast, g1_true; rtol=0.1)
        @test isapprox(result.amplitude_fast, A1_true; rtol=0.1)

        # Offset
        @test isapprox(result.offset, C_true; atol=1e-3)

        # Mode sorting: fast >= slow
        @test result.gap_fast >= result.gap

        @info "BIEXP-01" gap_slow=result.gap gap_fast=result.gap_fast offset=result.offset r2=result.r_squared
    end

    # -----------------------------------------------------------------------
    # BIEXP-02: Offset accuracy — bi-exp closer to true C than single-exp
    # -----------------------------------------------------------------------
    @testset "BIEXP-02: offset accuracy vs single-exp" begin
        # This is the key validation: bi-exp should give more accurate offset
        # when data has two timescales
        A1_true = 1.0    # fast
        g1_true = 2.0    # fast gap
        A2_true = 0.5    # slow
        g2_true = 0.3    # slow gap
        C_true  = 6.8e-5 # small offset (like floor from coherent unitary)

        times = collect(0.0:0.1:40.0)
        values = A1_true .* exp.(-g1_true .* times) .+
                 A2_true .* exp.(-g2_true .* times) .+
                 C_true

        single_fit = fit_exponential_decay(times, values; skip_initial=0.2)
        biexp_fit  = fit_biexponential_decay(times, values; skip_initial=0.2)

        single_err = abs(single_fit.offset - C_true)
        biexp_err  = abs(biexp_fit.offset - C_true)

        @test biexp_err < single_err
        @info "BIEXP-02 offset comparison" true_C=C_true single_C=single_fit.offset biexp_C=biexp_fit.offset single_err=single_err biexp_err=biexp_err
    end

    # -----------------------------------------------------------------------
    # BIEXP-03: skip_initial works with bi-exp
    # -----------------------------------------------------------------------
    @testset "BIEXP-03: skip_initial with bi-exp" begin
        A1_true = 1.0
        g1_true = 2.0
        A2_true = 0.5
        g2_true = 0.3
        C_true  = 0.001

        times = collect(0.0:0.05:30.0)
        values = A1_true .* exp.(-g1_true .* times) .+
                 A2_true .* exp.(-g2_true .* times) .+
                 C_true

        r1 = fit_biexponential_decay(times, values; skip_initial=0.0)
        r2 = fit_biexponential_decay(times, values; skip_initial=0.2)

        @test length(r2.times_used) < length(r1.times_used)
        # Both should recover reasonable gap
        @test isapprox(r1.gap, g2_true; rtol=0.1)
        @test isapprox(r2.gap, g2_true; rtol=0.1)
    end

    # -----------------------------------------------------------------------
    # BIEXP edge case: too few data points
    # -----------------------------------------------------------------------
    @testset "BIEXP: too few data points throws" begin
        times = collect(0.0:1.0:5.0)  # 6 points
        values = exp.(-0.3 .* times)
        @test_throws ArgumentError fit_biexponential_decay(times, values)
    end

end
