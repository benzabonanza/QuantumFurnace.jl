using Test
using Random
using LinearAlgebra

@testset "Deterministic multi-threaded results" begin
    if Threads.nthreads() > 1
        dim = size(SMALL_HAM.data, 1)
        CT = ComplexF64
        psi0 = zeros(CT, dim)
        psi0[1] = 1.0

        therm_config = make_small_thermalize_config(TimeDomain();
            delta=0.01, mixing_time=0.5, with_coherent=false)

        # Run twice with same seed => bitwise identical
        result1 = run_trajectories(SMALL_JUMPS, therm_config, psi0, SMALL_HAM;
            delta=0.01, ntraj=20, seed=42)
        result2 = run_trajectories(SMALL_JUMPS, therm_config, psi0, SMALL_HAM;
            delta=0.01, ntraj=20, seed=42)

        @test result1.rho_mean == result2.rho_mean  # bitwise identical, NOT isapprox
        @test result1.seed == result2.seed == 42
        @test result1.n_trajectories == 20

        # Different seed => different result
        result3 = run_trajectories(SMALL_JUMPS, therm_config, psi0, SMALL_HAM;
            delta=0.01, ntraj=20, seed=99)
        @test !(result1.rho_mean == result3.rho_mean)
    else
        @info "Skipping multi-thread determinism test (nthreads=$(Threads.nthreads()))"
        @test true  # placeholder so testset is not empty
    end
end

@testset "BLAS thread restoration" begin
    dim = size(SMALL_HAM.data, 1)
    CT = ComplexF64
    psi0 = zeros(CT, dim)
    psi0[1] = 1.0

    therm_config = make_small_thermalize_config(TimeDomain();
        delta=0.01, mixing_time=0.1, with_coherent=false)

    old_blas = BLAS.get_num_threads()
    result = run_trajectories(SMALL_JUMPS, therm_config, psi0, SMALL_HAM;
        delta=0.01, ntraj=10, seed=42)
    @test BLAS.get_num_threads() == old_blas
end

@testset "Serial-threaded agreement" begin
    if Threads.nthreads() > 1
        dim = size(SMALL_HAM.data, 1)
        CT = ComplexF64
        psi0 = zeros(CT, dim)
        psi0[1] = 1.0

        therm_config = make_small_thermalize_config(TimeDomain();
            delta=0.01, mixing_time=0.5, with_coherent=false)

        # Force serial: use julia -t 1 semantics by running ntraj=1
        # Instead, run with ntraj=20 and compare threaded vs explicit serial loop.
        # The threaded path uses Xoshiro(seed + traj_id) per trajectory.
        # We manually replicate the serial accumulation with the same per-trajectory seeds.
        seed = 42
        ntraj = 20
        result_threaded = run_trajectories(SMALL_JUMPS, therm_config, psi0, SMALL_HAM;
            delta=0.01, ntraj=ntraj, seed=seed)

        # Manual serial reference: accumulate density matrices with same per-trajectory seeds
        precomputed = QuantumFurnace._precompute_data(therm_config, SMALL_HAM)
        scratch = QuantumFurnace.KrausScratch(CT, dim)
        fw = build_trajectoryframework(SMALL_JUMPS, SMALL_HAM, therm_config, precomputed, scratch, 0.01)
        rho_ref = zeros(CT, dim, dim)
        for traj_id in 1:ntraj
            rng = Random.Xoshiro(seed + traj_id)
            ws_ref = QuantumFurnace.TrajectoryWorkspace(CT, dim)
            psi = copy(psi0)
            QuantumFurnace._evolve_along_trajectory!(psi, fw, ws_ref, rng, 0.5)
            QuantumFurnace._accumulate_density_matrix!(rho_ref, psi)
        end
        rho_ref ./= ntraj
        QuantumFurnace.hermitianize!(rho_ref)

        # Threaded accumulation order may differ, so use isapprox with tight tolerance
        @test isapprox(result_threaded.rho_mean, rho_ref; atol=1e-13)
    else
        @info "Skipping serial-threaded agreement test (nthreads=$(Threads.nthreads()))"
        @test true
    end
end

@testset "Deterministic observable path" begin
    if Threads.nthreads() > 1
        dim = size(SMALL_HAM.data, 1)
        CT = ComplexF64
        psi0 = zeros(CT, dim)
        psi0[1] = 1.0

        therm_config = make_small_thermalize_config(TimeDomain();
            delta=0.01, mixing_time=0.1, with_coherent=false)

        # Create a simple observable (Z on first qubit)
        obs = [Matrix{CT}(kron(Z, Matrix{Float64}(I, div(dim, 2), div(dim, 2))))]

        result1 = run_trajectories(SMALL_JUMPS, therm_config, psi0, SMALL_HAM;
            delta=0.01, ntraj=20, seed=42, observables=obs, save_every=5)
        result2 = run_trajectories(SMALL_JUMPS, therm_config, psi0, SMALL_HAM;
            delta=0.01, ntraj=20, seed=42, observables=obs, save_every=5)

        @test result1.rho_mean == result2.rho_mean  # bitwise identical
        @test result1.measurements_mean == result2.measurements_mean  # bitwise identical
        @test result1.times == result2.times
    else
        @info "Skipping observable determinism test (nthreads=$(Threads.nthreads()))"
        @test true
    end
end
