using StableRNGs

@testset "Mixing Time Estimation" begin

    # -----------------------------------------------------------------------
    # Helper: build a synthetic ThermalizeResults for testing
    # -----------------------------------------------------------------------
    function _make_synthetic_result(times::Vector{Float64}, dists::Vector{Float64};
                                     num_qubits::Int=3, mixing_time::Float64=50.0)
        config = make_config(Thermalize(), EnergyDomain();
                             num_qubits=num_qubits, mixing_time=mixing_time)
        dim = 2^num_qubits
        final_dm = zeros(ComplexF64, dim, dim)
        metadata = Dict{Symbol, Any}(:test => true)
        return ThermalizeResults(config, final_dm, dists, times, metadata)
    end

    # -----------------------------------------------------------------------
    # MIX-01: Exponential fit on synthetic trace distance curve
    # -----------------------------------------------------------------------
    @testset "MIX-01: synthetic gap recovery" begin
        A_true = 1.5
        gap_true = 0.3
        C_true = 0.001
        times = collect(0.0:0.1:50.0)
        dists = A_true .* exp.(-gap_true .* times) .+ C_true

        result = _make_synthetic_result(times, dists; mixing_time=50.0)
        est = estimate_mixing_time(result; skip_initial=0.1)

        @test isapprox(est.fitted_gap, gap_true; atol=0.01)
        @test est.r_squared > 0.99
        @test est.converged == true
        @test est.amplitude > 0
        @test est.offset >= 0
        @test est.fit_result isa FitResult

        @info "MIX-01 gap recovery" fitted=est.fitted_gap true_gap=gap_true r2=est.r_squared
    end

    # -----------------------------------------------------------------------
    # MIX-02: Extrapolation mode
    # -----------------------------------------------------------------------
    @testset "MIX-02: extrapolation mode" begin
        A_true = 1.5
        gap_true = 0.3
        C_true = 0.001
        times = collect(0.0:0.1:50.0)
        dists = A_true .* exp.(-gap_true .* times) .+ C_true

        result = _make_synthetic_result(times, dists; mixing_time=50.0)
        est = estimate_mixing_time(result; skip_initial=0.1, extrapolate=true, target_epsilon=0.01)

        # Expected extrapolated mixing time from the model:
        # t = -log((epsilon - C) / A) / gap
        t_expected = -log((0.01 - C_true) / A_true) / gap_true

        @test isapprox(est.mixing_time_extrapolated, t_expected; rtol=0.05)
        @test est.mixing_time == est.mixing_time_extrapolated
        @test est.target_epsilon == 0.01

        @info "MIX-02 extrapolation" t_extrap=est.mixing_time_extrapolated t_expected=t_expected
    end

    # -----------------------------------------------------------------------
    # MIX-03: Actual mixing time from data (no extrapolation)
    # -----------------------------------------------------------------------
    @testset "MIX-03: actual mixing time from data" begin
        A_true = 1.5
        gap_true = 0.3
        C_true = 0.001
        target = 0.01
        times = collect(0.0:0.1:50.0)
        dists = A_true .* exp.(-gap_true .* times) .+ C_true

        # The trace distance crosses 0.01 somewhere in this data
        # Find expected crossing: first index where dists <= target
        expected_idx = findfirst(d -> d <= target, dists)
        @assert expected_idx !== nothing "synthetic data should cross target=$target"
        expected_time = times[expected_idx]

        result = _make_synthetic_result(times, dists; mixing_time=50.0)
        est = estimate_mixing_time(result; skip_initial=0.1, target_epsilon=target, extrapolate=false)

        @test est.mixing_time_actual !== nothing
        @test isapprox(est.mixing_time_actual, expected_time; atol=0.15)
        @test est.mixing_time == est.mixing_time_actual
        @test est.mixing_time_extrapolated === nothing

        @info "MIX-03 actual crossing" t_actual=est.mixing_time_actual t_expected=expected_time
    end

    # -----------------------------------------------------------------------
    # MIX-04: skip_initial keyword
    # -----------------------------------------------------------------------
    @testset "MIX-04: skip_initial improves fit" begin
        gap_true = 0.3
        times = collect(0.0:0.05:30.0)
        # Data with fast transient added
        dists = 2.0 .* exp.(-gap_true .* times) .+ 0.5 .+ 3.0 .* exp.(-2.0 .* times)

        result = _make_synthetic_result(times, dists; mixing_time=30.0)
        est_no_skip = estimate_mixing_time(result; skip_initial=0.0)
        est_skip = estimate_mixing_time(result; skip_initial=0.3)

        @test abs(est_skip.fitted_gap - gap_true) < abs(est_no_skip.fitted_gap - gap_true)

        @info "MIX-04 skip_initial" gap_no_skip=est_no_skip.fitted_gap gap_skip=est_skip.fitted_gap true_gap=gap_true
    end

    # -----------------------------------------------------------------------
    # MIX-05: MixingTimeEstimate struct fields
    # -----------------------------------------------------------------------
    @testset "MIX-05: struct fields and types" begin
        times = collect(0.0:0.1:50.0)
        dists = 1.5 .* exp.(-0.3 .* times) .+ 0.001

        result = _make_synthetic_result(times, dists; mixing_time=50.0)
        est = estimate_mixing_time(result; skip_initial=0.1, target_epsilon=0.01, extrapolate=true)

        # Check all expected fields exist
        expected_fields = [
            :fitted_gap, :amplitude, :offset, :gap_ci, :gap_se, :r_squared,
            :converged, :mixing_time, :mixing_time_extrapolated, :mixing_time_actual,
            :target_epsilon, :fit_result
        ]
        for f in expected_fields
            @test f in fieldnames(MixingTimeEstimate)
        end

        # Check types
        @test est.fitted_gap isa Float64
        @test est.amplitude isa Float64
        @test est.offset isa Float64
        @test est.gap_ci isa Tuple{Float64, Float64}
        @test est.gap_se isa Float64
        @test est.r_squared isa Float64
        @test est.converged isa Bool
        @test est.mixing_time isa Float64
        @test est.target_epsilon isa Float64
        @test est.fit_result isa FitResult
    end

    # -----------------------------------------------------------------------
    # MIX-07: Quality gate warnings
    # -----------------------------------------------------------------------
    @testset "MIX-07: R-squared warning for bad fit" begin
        times = collect(0.0:0.1:50.0)
        # Use absolute value of sin to keep values positive (trace distances must be >= 0)
        dists = abs.(sin.(times)) .+ 0.01

        result = _make_synthetic_result(times, dists; mixing_time=50.0)
        @test_warn "R-squared" estimate_mixing_time(result; skip_initial=0.0)
    end

    @testset "MIX-07: offset warning for large C" begin
        times = collect(0.0:0.1:50.0)
        # Large offset C=1.0 relative to target_epsilon=0.5
        dists = 1.0 .* exp.(-0.3 .* times) .+ 1.0

        result = _make_synthetic_result(times, dists; mixing_time=50.0)
        @test_warn "offset" estimate_mixing_time(result; skip_initial=0.1, target_epsilon=0.5)
    end

    # -----------------------------------------------------------------------
    # Edge cases
    # -----------------------------------------------------------------------
    @testset "Edge: extrapolate=true without target_epsilon throws" begin
        times = collect(0.0:0.1:50.0)
        dists = 1.5 .* exp.(-0.3 .* times) .+ 0.001

        result = _make_synthetic_result(times, dists; mixing_time=50.0)
        @test_throws ArgumentError estimate_mixing_time(result; extrapolate=true)
    end

    @testset "Edge: fewer than 10 data points throws" begin
        times = collect(0.0:1.0:5.0)  # 6 points
        dists = 1.5 .* exp.(-0.3 .* times) .+ 0.001

        result = _make_synthetic_result(times, dists; mixing_time=5.0)
        @test_throws ArgumentError estimate_mixing_time(result)
    end

    @testset "Edge: target not reached in data" begin
        times = collect(0.0:0.1:5.0)  # short simulation
        dists = 1.5 .* exp.(-0.3 .* times) .+ 0.001  # doesn't reach 0.0001

        result = _make_synthetic_result(times, dists; mixing_time=5.0)
        est = estimate_mixing_time(result; target_epsilon=0.0001, extrapolate=false)

        @test est.mixing_time_actual === nothing
        @test isnan(est.mixing_time)  # NaN because actual not reached

        @info "Edge: target not reached" mixing_time=est.mixing_time actual=est.mixing_time_actual
    end

    # -----------------------------------------------------------------------
    # Integration test: real run_thermalize output
    # -----------------------------------------------------------------------
    @testset "Integration: real run_thermalize" begin
        config = make_config(Thermalize(), EnergyDomain();
                             num_qubits=3, mixing_time=5.0)
        result = run_thermalize(N3_JUMPS, config, N3_HAM)

        # Should NOT throw
        est = estimate_mixing_time(result)

        @test est.fitted_gap > 0
        @test est.r_squared > 0
        @test est.converged isa Bool
        @test est.fit_result isa FitResult

        @info "Integration test" gap=est.fitted_gap r2=est.r_squared converged=est.converged mixing_time=est.mixing_time
    end

end
