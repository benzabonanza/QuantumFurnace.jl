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

end
