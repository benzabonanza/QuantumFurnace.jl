using Test
using Random
using LinearAlgebra

@testset "Workspace Independence" begin
    # Use 3-qubit test fixtures from test_helpers.jl
    # (N3_HAM, N3_JUMPS already loaded by runtests.jl include order)

    dim = size(N3_HAM.data, 1)
    CT = ComplexF64

    # Build workspace (shared immutable data, independent scratch)
    therm_config = make_config(Thermalize(), TimeDomain(); num_qubits=3,
        delta=0.01, mixing_time=1.0, construction=GNS())
    ws_base = QuantumFurnace._build_trajectory_workspace(therm_config, N3_HAM, N3_JUMPS; delta=0.01)

    # Create two independent workspaces and RNGs
    ws1 = QuantumFurnace._copy_workspace_for_thread(ws_base)
    ws2 = QuantumFurnace._copy_workspace_for_thread(ws_base)
    rng1 = Random.Xoshiro(100)
    rng2 = Random.Xoshiro(200)

    # Same initial state
    psi0 = zeros(CT, dim)
    psi0[1] = 1.0
    psi1 = copy(psi0)
    psi2 = copy(psi0)

    # Step each workspace independently (10 steps)
    for _ in 1:10
        step_along_trajectory!(psi1, ws1, rng1)
        step_along_trajectory!(psi2, ws2, rng2)
    end

    # Different RNG seeds should produce different states (with very high probability)
    @test !isapprox(psi1, psi2; atol=1e-10)

    # Both states should be normalized
    norm1 = norm(psi1)
    norm2 = norm(psi2)
    @test isapprox(norm1, 1.0; atol=1e-10)  # Multi-step normalization: accumulated FP error over 10 unitary steps, N3_DIM * 10 * eps ~ 80 * 1e-16 ~ 8e-15
    @info "Workspace independence: psi1 norm" norm=norm1 threshold_atol=1e-10
    @test isapprox(norm2, 1.0; atol=1e-10)  # Same rationale as psi1
    @info "Workspace independence: psi2 norm" norm=norm2 threshold_atol=1e-10

    # Workspace scratch buffers should contain different data (they were used independently)
    @test !isapprox(ws1.scratch.psi_tmp, ws2.scratch.psi_tmp; atol=1e-10)

    # Verify determinism: re-run with same seed produces same result
    ws3 = QuantumFurnace._copy_workspace_for_thread(ws_base)
    rng3 = Random.Xoshiro(100)  # same seed as rng1
    psi3 = copy(psi0)
    for _ in 1:10
        step_along_trajectory!(psi3, ws3, rng3)
    end
    replay_err = norm(psi1 - psi3)
    @test isapprox(psi1, psi3; atol=1e-14)  # Deterministic replay: same seed + same workspace -> bitwise identical (1e-14 allows tiny platform FP differences)
    @info "Workspace independence: deterministic replay" error=replay_err threshold_atol=1e-14

    # Verify immutable data is shared (read-only)
    @test ws1.Rs === ws2.Rs  # same object reference
    @test ws1.K0s === ws2.K0s
    @test ws1.jumps === ws2.jumps
end

@testset "TrajectoryResult seed capture" begin
    dim = size(N3_HAM.data, 1)
    CT = ComplexF64
    psi0 = zeros(CT, dim)
    psi0[1] = 1.0

    therm_config = make_config(Thermalize(), TimeDomain(); num_qubits=3,
        delta=0.01, mixing_time=0.1, construction=GNS())

    # With explicit seed
    result1 = run_trajectories(N3_JUMPS, therm_config, psi0, N3_HAM;
        delta=0.01, ntraj=10, seed=42)
    @test result1 isa QuantumFurnace.TrajectoryResult
    @test result1.seed == 42
    @test result1.n_trajectories == 10
    @test size(result1.rho_mean) == (dim, dim)
    trace_val = real(tr(result1.rho_mean))
    @test isapprox(trace_val, 1.0; atol=1e-10)  # Density matrix trace: normalization from trajectory average (10 trajectories, accumulated FP rounding)
    @info "TrajectoryResult trace" value=trace_val threshold_atol=1e-10

    # Deterministic: same seed -> same result
    result2 = run_trajectories(N3_JUMPS, therm_config, psi0, N3_HAM;
        delta=0.01, ntraj=10, seed=42)
    determ_err = maximum(abs.(result1.rho_mean - result2.rho_mean))
    @test isapprox(result1.rho_mean, result2.rho_mean; atol=1e-14)  # Deterministic: same seed + same code -> bitwise identical across runs
    @info "TrajectoryResult determinism" max_element_error=determ_err threshold_atol=1e-14

    # Without seed: auto-generated, stored in result
    result3 = run_trajectories(N3_JUMPS, therm_config, psi0, N3_HAM;
        delta=0.01, ntraj=10)
    @test result3.seed != 0  # seed was auto-generated
    @test result3.times === nothing  # no observables
    @test result3.measurements_mean === nothing
end
