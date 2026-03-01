using Test
using Random
using LinearAlgebra

@testset "Deterministic multi-threaded results" begin
    if Threads.nthreads() > 1
        dim = size(N3_HAM.data, 1)
        CT = ComplexF64
        psi0 = zeros(CT, dim)
        psi0[1] = 1.0

        therm_config = make_config(Thermalize(), TimeDomain(); num_qubits=3,
            delta=0.01, mixing_time=0.5, construction=GNS())

        # Run twice with same seed => bitwise identical
        result1 = run_trajectories(N3_JUMPS, therm_config, psi0, N3_HAM;
            delta=0.01, ntraj=20, seed=42)
        result2 = run_trajectories(N3_JUMPS, therm_config, psi0, N3_HAM;
            delta=0.01, ntraj=20, seed=42)

        @test result1.rho_mean == result2.rho_mean  # bitwise identical, NOT isapprox
        @test result1.seed == result2.seed == 42
        @test result1.n_trajectories == 20

        # Different seed => different result
        result3 = run_trajectories(N3_JUMPS, therm_config, psi0, N3_HAM;
            delta=0.01, ntraj=20, seed=99)
        @test !(result1.rho_mean == result3.rho_mean)
    else
        @info "Skipping multi-thread determinism test (nthreads=$(Threads.nthreads()))"
        @test true  # placeholder so testset is not empty
    end
end

@testset "BLAS thread restoration" begin
    dim = size(N3_HAM.data, 1)
    CT = ComplexF64
    psi0 = zeros(CT, dim)
    psi0[1] = 1.0

    therm_config = make_config(Thermalize(), TimeDomain(); num_qubits=3,
        delta=0.01, mixing_time=0.1, construction=GNS())

    old_blas = BLAS.get_num_threads()
    result = run_trajectories(N3_JUMPS, therm_config, psi0, N3_HAM;
        delta=0.01, ntraj=10, seed=42)
    @test BLAS.get_num_threads() == old_blas
end

@testset "Serial-threaded agreement" begin
    if Threads.nthreads() > 1
        dim = size(N3_HAM.data, 1)
        CT = ComplexF64
        psi0 = zeros(CT, dim)
        psi0[1] = 1.0

        therm_config = make_config(Thermalize(), TimeDomain(); num_qubits=3,
            delta=0.01, mixing_time=0.5, construction=GNS())

        # Force serial: use julia -t 1 semantics by running ntraj=1
        # Instead, run with ntraj=20 and compare threaded vs explicit serial loop.
        # The threaded path uses Xoshiro(seed + traj_id) per trajectory.
        # We manually replicate the serial accumulation with the same per-trajectory seeds.
        seed = 42
        ntraj = 20
        result_threaded = run_trajectories(N3_JUMPS, therm_config, psi0, N3_HAM;
            delta=0.01, ntraj=ntraj, seed=seed)

        # Manual serial reference: accumulate density matrices with same per-trajectory seeds
        ws = QuantumFurnace._build_trajectory_workspace(therm_config, N3_HAM, N3_JUMPS; delta=0.01)
        num_steps = ceil(Int, 0.5 / ws.delta)
        rho_ref = zeros(CT, dim, dim)
        for traj_id in 1:ntraj
            rng = Random.Xoshiro(seed + traj_id)
            ws_ref = QuantumFurnace._copy_workspace_for_thread(ws)
            psi = copy(psi0)
            psi_norm2 = real(dot(psi, psi))
            rmul!(psi, 1.0 / sqrt(max(psi_norm2, eps(Float64))))
            for _ in 1:num_steps
                step_along_trajectory!(psi, ws_ref, rng)
            end
            QuantumFurnace._accumulate_density_matrix!(rho_ref, psi)
        end
        rho_ref ./= ntraj
        QuantumFurnace.hermitianize!(rho_ref)

        # Threaded accumulation order may differ, so use isapprox with tight tolerance
        agreement_err = maximum(abs.(result_threaded.rho_mean - rho_ref))
        @test isapprox(result_threaded.rho_mean, rho_ref; atol=1e-13)  # FP accumulation order difference between threaded and serial: O(ntraj * DIM^2 * eps) ~ 20 * 64 * 1e-16 ~ 1e-13
        @info "Serial-threaded agreement" max_element_error=agreement_err threshold_atol=1e-13 ntraj=ntraj
    else
        @info "Skipping serial-threaded agreement test (nthreads=$(Threads.nthreads()))"
        @test true
    end
end

@testset "Deterministic observable path" begin
    if Threads.nthreads() > 1
        dim = size(N3_HAM.data, 1)
        CT = ComplexF64
        psi0 = zeros(CT, dim)
        psi0[1] = 1.0

        therm_config = make_config(Thermalize(), TimeDomain(); num_qubits=3,
            delta=0.01, mixing_time=0.1, construction=GNS())

        # Create a simple observable (Z on first qubit)
        obs = [Matrix{CT}(kron(Z, Matrix{Float64}(I, div(dim, 2), div(dim, 2))))]

        result1 = run_trajectories(N3_JUMPS, therm_config, psi0, N3_HAM;
            delta=0.01, ntraj=20, seed=42, observables=obs, save_every=5)
        result2 = run_trajectories(N3_JUMPS, therm_config, psi0, N3_HAM;
            delta=0.01, ntraj=20, seed=42, observables=obs, save_every=5)

        @test result1.rho_mean == result2.rho_mean  # bitwise identical
        @test result1.measurements_mean == result2.measurements_mean  # bitwise identical
        @test result1.times == result2.times
    else
        @info "Skipping observable determinism test (nthreads=$(Threads.nthreads()))"
        @test true
    end
end

@testset "Threading speedup" begin
    if Threads.nthreads() >= 2
        dim = size(N3_HAM.data, 1)
        CT = ComplexF64
        psi0 = zeros(CT, dim)
        psi0[1] = 1.0

        # Use longer mixing_time and more trajectories to amortize threading overhead.
        # At dim=8 each step is very fast, so we need enough total work to see speedup.
        therm_config = make_config(Thermalize(), TimeDomain(); num_qubits=3,
            delta=0.01, mixing_time=10.0, construction=GNS())

        ntraj = 2000

        # Warmup both paths
        run_trajectories(N3_JUMPS, therm_config, psi0, N3_HAM;
            delta=0.01, ntraj=ntraj, seed=1)

        # Measure threaded execution time (automatic with nthreads > 1)
        t_threaded = @elapsed begin
            for _ in 1:3  # average over 3 runs
                run_trajectories(N3_JUMPS, therm_config, psi0, N3_HAM;
                    delta=0.01, ntraj=ntraj, seed=42)
            end
        end
        t_threaded /= 3

        # Measure serial execution time: force single-thread behavior by running
        # trajectories one at a time in a loop, using inline step loop directly
        # This simulates serial execution regardless of thread count
        ws = QuantumFurnace._build_trajectory_workspace(therm_config, N3_HAM, N3_JUMPS; delta=0.01)
        num_steps_perf = ceil(Int, 5.0 / ws.delta)

        # Warmup serial path
        ws_s = QuantumFurnace._copy_workspace_for_thread(ws)
        psi_s = copy(psi0)
        rng_s = Random.Xoshiro(1)
        psi_s_norm2 = real(dot(psi_s, psi_s))
        rmul!(psi_s, 1.0 / sqrt(max(psi_s_norm2, eps(Float64))))
        for _ in 1:num_steps_perf
            step_along_trajectory!(psi_s, ws_s, rng_s)
        end

        t_serial = @elapsed begin
            for _ in 1:3
                rho_acc = zeros(CT, dim, dim)
                for traj_id in 1:ntraj
                    ws_loop = QuantumFurnace._copy_workspace_for_thread(ws)
                    rng_loop = Random.Xoshiro(42 + traj_id)
                    psi_loop = copy(psi0)
                    psi_loop_norm2 = real(dot(psi_loop, psi_loop))
                    rmul!(psi_loop, 1.0 / sqrt(max(psi_loop_norm2, eps(Float64))))
                    for _ in 1:num_steps_perf
                        step_along_trajectory!(psi_loop, ws_loop, rng_loop)
                    end
                    QuantumFurnace._accumulate_density_matrix!(rho_acc, psi_loop)
                end
            end
        end
        t_serial /= 3

        @info "Threading performance" nthreads=Threads.nthreads() t_serial t_threaded speedup=t_serial/t_threaded

        # Threaded should be faster than serial (no regression from threading overhead)
        @test t_threaded < t_serial  # Timing comparison: system-dependent, no fixed threshold -- just verifies threading doesn't regress
        @info "Threading speedup test" t_threaded t_serial passed=(t_threaded < t_serial)
    else
        @info "Skipping threading speedup test (nthreads=$(Threads.nthreads()))"
        @test true
    end
end

# ============================================================================
# DM thermalization BLAS threading tests (THREAD-03, THREAD-05)
# ============================================================================

@testset "DM BLAS thread restoration" begin
    old_blas = BLAS.get_num_threads()

    # Normal completion (EnergyDomain)
    therm_config = make_config(Thermalize(), EnergyDomain(); num_qubits=3,
        delta=0.01, mixing_time=0.1)
    result = run_thermalize(N3_JUMPS, therm_config, N3_HAM)
    @test BLAS.get_num_threads() == old_blas
    @info "DM BLAS restoration (Energy)" blas_before=old_blas blas_after=BLAS.get_num_threads()

    # Test across multiple domains to ensure all paths restore
    for (domain, jumps, trott, name) in [
        (TimeDomain(), N3_JUMPS, nothing, "Time"),
        (TrotterDomain(), N3_TROTTER_JUMPS, N3_TROTTER, "Trotter"),
    ]
        cfg = make_config(Thermalize(), domain; num_qubits=3, delta=0.01, mixing_time=0.1)
        run_thermalize(jumps, cfg, N3_HAM, trott)
        @test BLAS.get_num_threads() == old_blas
        @info "DM BLAS restoration ($name)" blas_after=BLAS.get_num_threads()
    end
end

@testset "DM serial-threaded BLAS agreement" begin
    # Run with BLAS threads = 1 (serial BLAS) vs default (multi-threaded BLAS)
    therm_config = make_config(Thermalize(), EnergyDomain(); num_qubits=3,
        delta=0.01, mixing_time=0.5)
    rng_seed = 42

    # Serial BLAS reference
    old_blas = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    result_serial = run_thermalize(N3_JUMPS, therm_config, N3_HAM;
        rng=Random.Xoshiro(rng_seed))
    BLAS.set_num_threads(old_blas)

    # Multi-threaded BLAS
    result_multi = run_thermalize(N3_JUMPS, therm_config, N3_HAM;
        rng=Random.Xoshiro(rng_seed))

    # BLAS threading does NOT change FP results for BLAS gemm (BLAS is deterministic
    # at fixed thread count, but may differ between 1 and N threads due to different
    # reduction order in gemm). At dim=8 (3 qubits), the difference is negligible.
    # For larger systems the FP accumulation order in gemm may differ.
    td_diff = maximum(abs.(result_serial.trace_distances .- result_multi.trace_distances))
    dm_diff = maximum(abs.(result_serial.final_dm .- result_multi.final_dm))
    @test isapprox(result_serial.trace_distances, result_multi.trace_distances; atol=1e-10)
    @test isapprox(result_serial.final_dm, result_multi.final_dm; atol=1e-10)
    @info "DM serial-threaded BLAS agreement" trace_dist_diff=td_diff dm_diff=dm_diff threshold=1e-10
end

@testset "DM-trajectory BLAS isolation" begin
    # Verify that calling run_thermalize then run_trajectories does not leak BLAS state
    old_blas = BLAS.get_num_threads()

    therm_config_dm = make_config(Thermalize(), EnergyDomain(); num_qubits=3,
        delta=0.01, mixing_time=0.1)
    run_thermalize(N3_JUMPS, therm_config_dm, N3_HAM)
    @test BLAS.get_num_threads() == old_blas

    therm_config_traj = make_config(Thermalize(), TimeDomain(); num_qubits=3,
        delta=0.01, mixing_time=0.1, construction=GNS())
    dim = size(N3_HAM.data, 1)
    psi0 = zeros(ComplexF64, dim); psi0[1] = 1.0
    run_trajectories(N3_JUMPS, therm_config_traj, psi0, N3_HAM;
        delta=0.01, ntraj=10, seed=42)
    @test BLAS.get_num_threads() == old_blas

    @info "DM-trajectory BLAS isolation" blas_threads=old_blas verified=true
end
