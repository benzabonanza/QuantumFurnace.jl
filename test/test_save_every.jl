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
        @test result_default.final_dm == result_se1.final_dm
    end

    # -----------------------------------------------------------------------
    # Test 2: exact recorded schedule and array lengths
    # -----------------------------------------------------------------------
    @testset "exact saved schedule" begin
        config = make_config(Thermalize(), EnergyDomain(); num_qubits=3, mixing_time=0.2)
        save_every = 5
        result = run_thermalize(N3_JUMPS, config, N3_HAM; rng=Random.MersenneTwister(42), save_every=save_every)
        expected_times = collect(0.0:0.05:0.2)
        @test result.time_steps ≈ expected_times atol=1e-15
        @test length(result.trace_distances) == length(expected_times)
        @test length(result.time_steps) == length(result.trace_distances)
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

end
