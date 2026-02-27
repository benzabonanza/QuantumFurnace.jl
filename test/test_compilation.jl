@testset "Compilation and Loading" begin
    @testset "Module loads without errors" begin
        @test true  # If we got here, `using QuantumFurnace` succeeded
    end

    @testset "_build_trajectory_workspace with coherent (KMS)" begin
        config = make_thermalize_config(EnergyDomain(); construction=KMS())
        ws = QuantumFurnace._build_trajectory_workspace(config, TEST_HAM, TEST_JUMPS; delta=TEST_DELTA)
        @test ws isa Workspace{Trajectory}
        @test ws.delta == TEST_DELTA
        # Per-operator: each operator should have a coherent unitary
        @test all(u -> u !== nothing, ws.U_Bs)
    end

    @testset "_build_trajectory_workspace without coherent (GNS)" begin
        config = make_thermalize_config(EnergyDomain(); construction=GNS())
        ws = QuantumFurnace._build_trajectory_workspace(config, TEST_HAM, TEST_JUMPS; delta=TEST_DELTA)
        @test ws isa Workspace{Trajectory}
        # Per-operator: no coherent unitaries
        @test all(u -> u === nothing, ws.U_Bs)
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
        @test lc isa Config{Lindbladian}
        @test with_coherent(lc.construction) == true
        @test lc.num_qubits == NUM_QUBITS

        tc = make_thermalize_config(EnergyDomain(); construction=GNS(), delta=0.05)
        @test tc isa Config{Thermalize}
        @test with_coherent(tc.construction) == false
        @test tc.delta == 0.05
    end
end
