using Test
using Random
using LinearAlgebra

@testset "Workspace Independence" begin
    # Use 3-qubit test fixtures from test_helpers.jl
    # (SMALL_HAM, SMALL_JUMPS already loaded by runtests.jl include order)

    dim = size(SMALL_HAM.data, 1)
    CT = ComplexF64

    # Build framework (shared, read-only)
    therm_config = make_small_thermalize_config(TimeDomain();
        delta=0.01, mixing_time=1.0, construction=GNS())
    precomputed = QuantumFurnace._precompute_data(therm_config, SMALL_HAM)
    scratch = QuantumFurnace.KrausScratch(CT, dim)
    fw = build_trajectoryframework(SMALL_JUMPS, SMALL_HAM, therm_config, precomputed, scratch, 0.01)

    # Create two independent workspaces and RNGs
    ws1 = QuantumFurnace.TrajectoryWorkspace(fw)
    ws2 = QuantumFurnace.TrajectoryWorkspace(fw)
    rng1 = Random.Xoshiro(100)
    rng2 = Random.Xoshiro(200)

    # Same initial state
    psi0 = zeros(CT, dim)
    psi0[1] = 1.0
    psi1 = copy(psi0)
    psi2 = copy(psi0)

    # Step each workspace independently (10 steps)
    for _ in 1:10
        step_along_trajectory!(psi1, fw, ws1, rng1)
        step_along_trajectory!(psi2, fw, ws2, rng2)
    end

    # Different RNG seeds should produce different states (with very high probability)
    @test !isapprox(psi1, psi2; atol=1e-10)

    # Both states should be normalized
    @test isapprox(norm(psi1), 1.0; atol=1e-10)
    @test isapprox(norm(psi2), 1.0; atol=1e-10)

    # Workspace buffers should contain different data (they were used independently)
    @test !isapprox(ws1.psi_tmp, ws2.psi_tmp; atol=1e-10)

    # Verify determinism: re-run with same seed produces same result
    ws3 = QuantumFurnace.TrajectoryWorkspace(fw)
    rng3 = Random.Xoshiro(100)  # same seed as rng1
    psi3 = copy(psi0)
    for _ in 1:10
        step_along_trajectory!(psi3, fw, ws3, rng3)
    end
    @test isapprox(psi1, psi3; atol=1e-14)  # deterministic replay

    # Verify framework was not mutated (read-only)
    # Re-check that we can build a workspace from fw (fw still valid)
    ws4 = QuantumFurnace.TrajectoryWorkspace(fw)
    @test size(ws4.jump_oft) == size(ws1.jump_oft)
end

@testset "TrajectoryResult seed capture" begin
    dim = size(SMALL_HAM.data, 1)
    CT = ComplexF64
    psi0 = zeros(CT, dim)
    psi0[1] = 1.0

    therm_config = make_small_thermalize_config(TimeDomain();
        delta=0.01, mixing_time=0.1, construction=GNS())

    # With explicit seed
    result1 = run_trajectories(SMALL_JUMPS, therm_config, psi0, SMALL_HAM;
        delta=0.01, ntraj=10, seed=42)
    @test result1 isa QuantumFurnace.TrajectoryResult
    @test result1.seed == 42
    @test result1.n_trajectories == 10
    @test size(result1.rho_mean) == (dim, dim)
    @test isapprox(tr(result1.rho_mean), 1.0; atol=1e-10)

    # Deterministic: same seed -> same result
    result2 = run_trajectories(SMALL_JUMPS, therm_config, psi0, SMALL_HAM;
        delta=0.01, ntraj=10, seed=42)
    @test isapprox(result1.rho_mean, result2.rho_mean; atol=1e-14)

    # Without seed: auto-generated, stored in result
    result3 = run_trajectories(SMALL_JUMPS, therm_config, psi0, SMALL_HAM;
        delta=0.01, ntraj=10)
    @test result3.seed != 0  # seed was auto-generated
    @test result3.times === nothing  # no observables
    @test result3.measurements_mean === nothing
end
