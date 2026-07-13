using Test
using Random
using LinearAlgebra

# All tests in this file are exact bitwise cross-validations (same seed, same results)
# or structural checks (type, size, field existence). No numerical threshold comparisons,
# so no @info additions needed per test output policy.

@testset "Observable-Only Trajectory Runner" begin
    # Setup: 3-qubit test system (dim=8)
    dim = N3_DIM  # 8
    CT = ComplexF64
    psi0 = zeros(CT, dim); psi0[1] = 1.0
    therm_config = make_config(Thermalize(), TimeDomain(); num_qubits=3, delta=0.01, mixing_time=0.5, construction=GNS())

    # Z on first qubit in eigenbasis
    Z1_comp = kron(Z, Matrix{Float64}(I, div(dim, 2), div(dim, 2)))
    Z1_eigen = Matrix{CT}(N3_HAM.eigvecs' * Z1_comp * N3_HAM.eigvecs)
    obs = [Z1_eigen]

    # Precompute expected num_saves
    delta_step = 0.01
    num_steps = ceil(Int, 0.5 / delta_step)
    num_saves = div(num_steps, 5) + 1

    @testset "Basic observable-only run" begin
        result = run_observable_trajectories(N3_JUMPS, therm_config, psi0, N3_HAM;
            observables=obs, save_every=5, ntraj=10, seed=42)

        @test result isa ObservableTrajectoryResult
        @test result.rho_mean === nothing
        @test result.n_trajectories == 10
        @test result.seed == 42
        @test size(result.measurements_mean) == (1, num_saves)
        @test length(result.times) == num_saves
        @test result.times[1] == 0.0
        @test all(isfinite, result.measurements_mean)
    end

    @testset "Cross-validation: bitwise match with run_trajectories" begin
        result_old = run_trajectories(N3_JUMPS, therm_config, psi0, N3_HAM;
            observables=obs, save_every=5, ntraj=20, seed=42)
        result_new = run_observable_trajectories(N3_JUMPS, therm_config, psi0, N3_HAM;
            observables=obs, save_every=5, ntraj=20, seed=42)

        @test result_new.measurements_mean == result_old.measurements_mean
        @test result_new.times == result_old.times
    end

    @testset "Cross-validation with reconstruct_dm=true" begin
        result_old = run_trajectories(N3_JUMPS, therm_config, psi0, N3_HAM;
            observables=obs, save_every=5, ntraj=20, seed=42)
        result_dm = run_observable_trajectories(N3_JUMPS, therm_config, psi0, N3_HAM;
            observables=obs, save_every=5, ntraj=20, seed=42, reconstruct_dm=true)

        @test result_dm.rho_mean !== nothing
        @test result_dm.measurements_mean == result_old.measurements_mean
        @test result_dm.rho_mean == result_old.rho_mean
        @test result_dm.times == result_old.times
    end

    @testset "Deterministic: same seed gives same results" begin
        r1 = run_observable_trajectories(N3_JUMPS, therm_config, psi0, N3_HAM;
            observables=obs, save_every=5, ntraj=20, seed=42)
        r2 = run_observable_trajectories(N3_JUMPS, therm_config, psi0, N3_HAM;
            observables=obs, save_every=5, ntraj=20, seed=42)

        @test r1.measurements_mean == r2.measurements_mean
        @test r1.times == r2.times
    end

    @testset "Different seed gives different results" begin
        r1 = run_observable_trajectories(N3_JUMPS, therm_config, psi0, N3_HAM;
            observables=obs, save_every=5, ntraj=20, seed=42)
        r2 = run_observable_trajectories(N3_JUMPS, therm_config, psi0, N3_HAM;
            observables=obs, save_every=5, ntraj=20, seed=99)

        @test !(r1.measurements_mean == r2.measurements_mean)
    end

    @testset "Single trajectory (serial path)" begin
        result_old = run_trajectories(N3_JUMPS, therm_config, psi0, N3_HAM;
            observables=obs, save_every=5, ntraj=1, seed=42)
        result_new = run_observable_trajectories(N3_JUMPS, therm_config, psi0, N3_HAM;
            observables=obs, save_every=5, ntraj=1, seed=42)

        @test result_new.measurements_mean == result_old.measurements_mean
    end
end
