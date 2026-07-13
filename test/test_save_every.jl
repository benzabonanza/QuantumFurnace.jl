using Test
using Random
using LinearAlgebra

@testset "save_every" begin

    # -----------------------------------------------------------------------
    # Test 1: Backward compatibility -- save_every=1 matches default behavior
    # -----------------------------------------------------------------------
    @testset "save_every=1 matches default" begin
        config = make_config(Thermalize(), EnergyDomain(); num_qubits=3, mixing_time=0.1)
        rng1 = Random.MersenneTwister(42)
        rng2 = Random.MersenneTwister(42)
        result_default = run_thermalize(N3_JUMPS, config, N3_HAM; rng=rng1)
        result_se1    = run_thermalize(N3_JUMPS, config, N3_HAM; rng=rng2, save_every=1)
        @test result_default.trace_distances == result_se1.trace_distances
        @test result_default.time_steps == result_se1.time_steps
    end

    # -----------------------------------------------------------------------
    # Test 2: save_every=10 produces correct array lengths
    # -----------------------------------------------------------------------
    @testset "save_every=10 array lengths" begin
        config = make_config(Thermalize(), EnergyDomain(); num_qubits=3, mixing_time=0.2)
        save_every = 10
        num_steps = Int(ceil(config.mixing_time / config.delta))
        expected_saves = 1 + div(num_steps, save_every)  # +1 for initial state
        result = run_thermalize(N3_JUMPS, config, N3_HAM; rng=Random.MersenneTwister(42), save_every=save_every)
        # Allow for early convergence (fewer saves)
        @test length(result.trace_distances) <= expected_saves
        @test length(result.trace_distances) >= 2  # at least initial + one save point
        @test length(result.time_steps) == length(result.trace_distances)
    end

    # -----------------------------------------------------------------------
    # Test 3: time_steps stride is correct
    # -----------------------------------------------------------------------
    @testset "time_steps stride" begin
        config = make_config(Thermalize(), EnergyDomain(); num_qubits=3, mixing_time=0.2)
        save_every = 5
        result = run_thermalize(N3_JUMPS, config, N3_HAM; rng=Random.MersenneTwister(42), save_every=save_every)
        @test result.time_steps[1] == 0.0  # initial state at t=0
        if length(result.time_steps) >= 3
            stride = result.time_steps[2] - result.time_steps[1]
            expected_stride = save_every * config.delta
            @test isapprox(stride, expected_stride; atol=1e-15)
            # All interior strides should be uniform
            for i in 2:(length(result.time_steps) - 1)
                @test isapprox(result.time_steps[i+1] - result.time_steps[i], expected_stride; atol=1e-15)
            end
        end
    end

    # -----------------------------------------------------------------------
    # Test 4: save_every stored in metadata
    # -----------------------------------------------------------------------
    @testset "metadata contains save_every" begin
        config = make_config(Thermalize(), EnergyDomain(); num_qubits=3, mixing_time=0.1)
        result = run_thermalize(N3_JUMPS, config, N3_HAM; rng=Random.MersenneTwister(42), save_every=7)
        @test haskey(result.metadata, :save_every)
        @test result.metadata[:save_every] == 7
    end

    # -----------------------------------------------------------------------
    # Test 5: save_every < 1 throws assertion error
    # -----------------------------------------------------------------------
    @testset "save_every validation" begin
        config = make_config(Thermalize(), EnergyDomain(); num_qubits=3, mixing_time=0.1)
        @test_throws AssertionError run_thermalize(N3_JUMPS, config, N3_HAM; rng=Random.MersenneTwister(42), save_every=0)
        @test_throws AssertionError run_thermalize(N3_JUMPS, config, N3_HAM; rng=Random.MersenneTwister(42), save_every=-1)
    end

    # -----------------------------------------------------------------------
    # Test 6: save_every with TimeDomain (cross-domain check)
    # -----------------------------------------------------------------------
    @testset "save_every with TimeDomain" begin
        config = make_config(Thermalize(), TimeDomain(); num_qubits=3, mixing_time=0.1)
        result = run_thermalize(N3_JUMPS, config, N3_HAM; rng=Random.MersenneTwister(42), save_every=5)
        @test length(result.time_steps) == length(result.trace_distances)
        @test result.time_steps[1] == 0.0
        @test haskey(result.metadata, :save_every)
    end

end
