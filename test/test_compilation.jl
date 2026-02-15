@testset "Compilation and Loading" begin
    @testset "Module loads without errors" begin
        @test true  # If we got here, `using QuantumFurnace` succeeded
    end

    @testset "build_trajectoryframework with coherent" begin
        config = make_thermalize_config(EnergyDomain(); with_coherent=true)
        precomputed = precompute_data(config, TEST_HAM)
        scratch = KrausScratch(ComplexF64, DIM)
        fw = build_trajectoryframework(
            TEST_JUMPS, TEST_HAM, config, precomputed, scratch, TEST_DELTA
        )
        @test fw isa TrajectoryFramework
        @test fw.delta == TEST_DELTA
        # Per-operator: each operator should have a coherent unitary
        @test all(per_op -> per_op.U_B !== nothing, fw.per_operator)
    end

    @testset "build_trajectoryframework without coherent" begin
        config = make_thermalize_config(EnergyDomain(); with_coherent=false)
        precomputed = precompute_data(config, TEST_HAM)
        scratch = KrausScratch(ComplexF64, DIM)
        fw = build_trajectoryframework(
            TEST_JUMPS, TEST_HAM, config, precomputed, scratch, TEST_DELTA
        )
        @test fw isa TrajectoryFramework
        # Per-operator: no coherent unitaries
        @test all(per_op -> per_op.U_B === nothing, fw.per_operator)
    end

    @testset "Fixtures available" begin
        @test size(TEST_HAM.data) == (DIM, DIM)
        @test length(TEST_JUMPS) == 3 * NUM_QUBITS  # 12 jumps
        @test size(TEST_GIBBS) == (DIM, DIM)
        @test isapprox(tr(TEST_GIBBS), 1.0; atol=TOL_EXACT)
    end

    @testset "Tolerance tiers defined" begin
        @test TOL_EXACT == 1e-12
        @test TOL_QUADRATURE == 1e-6
        @test TOL_DELTA(0.01) == 0.05
        @test TOL_DELTA(0.1) == 0.5
    end

    @testset "Config factories" begin
        lc = make_liouv_config(EnergyDomain())
        @test lc isa LiouvConfig
        @test lc.with_coherent == true
        @test lc.num_qubits == NUM_QUBITS

        tc = make_thermalize_config(EnergyDomain(); with_coherent=false, delta=0.05)
        @test tc isa ThermalizeConfig
        @test tc.with_coherent == false
        @test tc.delta == 0.05
    end
end
