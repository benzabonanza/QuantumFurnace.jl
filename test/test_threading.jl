using Test
using Random
using LinearAlgebra

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
        delta=0.01, mixing_time=0.05)
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

# ============================================================================
# Omega-loop threading tests (THREAD-01, THREAD-02, THREAD-04)
# ============================================================================

@testset "Omega-loop threading determinism" begin
    if Threads.nthreads() > 1
        for (domain, jumps, trott, name) in [
            (EnergyDomain(), N3_JUMPS, nothing, "Energy"),
            (TimeDomain(), N3_JUMPS, nothing, "Time"),
            (TrotterDomain(), N3_TROTTER_JUMPS, N3_TROTTER, "Trotter"),
            (BohrDomain(), N3_JUMPS, nothing, "Bohr"),
        ]
            cfg = make_config(Thermalize(), domain; num_qubits=3,
                delta=0.01, mixing_time=0.05)
            rng_seed = 42

            result1 = run_thermalize(jumps, cfg, N3_HAM, trott;
                rng=Random.Xoshiro(rng_seed))
            result2 = run_thermalize(jumps, cfg, N3_HAM, trott;
                rng=Random.Xoshiro(rng_seed))

            # Deterministic: same seed, same thread count => identical results
            @test result1.trace_distances == result2.trace_distances
            @test result1.final_dm == result2.final_dm
            @info "Omega-loop threading determinism ($name)" domain=name passed=true
        end
    else
        @info "Skipping omega-loop threading determinism tests (nthreads=$(Threads.nthreads()))"
        @test true
    end
end

@testset "Serial vs threaded omega-loop agreement" begin
    if Threads.nthreads() > 1
        cfg = make_config(Thermalize(), EnergyDomain(); num_qubits=3,
            delta=0.01, mixing_time=0.05)
        rng_seed = 42

        # Run with threading enabled (default path on -t 4)
        result_threaded = run_thermalize(N3_JUMPS, cfg, N3_HAM;
            rng=Random.Xoshiro(rng_seed))

        # Compare against run with BLAS threads forced to 1
        # (omega-loop threading sets BLAS=1 internally, so BLAS thread count
        # only affects the outer run_thermalize try/finally block)
        old_blas = BLAS.get_num_threads()
        BLAS.set_num_threads(1)
        result_serial_blas = run_thermalize(N3_JUMPS, cfg, N3_HAM;
            rng=Random.Xoshiro(rng_seed))
        BLAS.set_num_threads(old_blas)

        # Threaded omega-loop + multi-BLAS vs threaded omega-loop + serial BLAS
        # should agree within FP tolerance (BLAS thread count affects gemm
        # accumulation order but not correctness)
        @test isapprox(result_threaded.trace_distances,
            result_serial_blas.trace_distances; atol=1e-10)
        @test isapprox(result_threaded.final_dm,
            result_serial_blas.final_dm; atol=1e-10)

        td_diff = maximum(abs.(result_threaded.trace_distances .- result_serial_blas.trace_distances))
        dm_diff = maximum(abs.(result_threaded.final_dm .- result_serial_blas.final_dm))
        @info "Serial vs threaded omega-loop agreement" td_diff=td_diff dm_diff=dm_diff atol=1e-10
    else
        @info "Skipping serial vs threaded omega-loop test (nthreads=$(Threads.nthreads()))"
        @test true
    end
end

@testset "Omega-loop BLAS restoration" begin
    old_blas = BLAS.get_num_threads()
    cfg = make_config(Thermalize(), EnergyDomain(); num_qubits=3,
        delta=0.01, mixing_time=0.05)
    run_thermalize(N3_JUMPS, cfg, N3_HAM)
    @test BLAS.get_num_threads() == old_blas
    @info "Omega-loop BLAS restoration" blas_before=old_blas blas_after=BLAS.get_num_threads()
end

# ============================================================================
# Lindbladian-matvec ω-loop threading tests (qf-in3)
# These exercise `_apply_lindbladian_threaded_energy!` and
# `_apply_lindbladian_threaded_timetrot!` directly, comparing their output to
# the serial public path. Threading is correctness-preserving up to floating-
# point summation order, so the threshold for the bit-match check is set
# relative to ‖L(rho)‖.
# ============================================================================

# Helper: run threaded variant directly and return the resulting rho_out.
# Mirrors what `apply_lindbladian!` does internally when the threshold dispatch
# fires, so bypasses the threshold gate.
function _run_threaded_lindbladian!(ws, rho, config, ham; adjoint::Bool=false)
    sc = ws.scratch::QuantumFurnace.KrylovScratch{ComplexF64}
    if adjoint
        G_left, G_right = ws.G_left_adj, ws.G_right_adj
    else
        G_left, G_right = ws.G_left, ws.G_right
    end
    fill!(sc.rho_out, 0)
    mul!(sc.rho_out, G_left, rho)
    mul!(sc.rho_out, rho, G_right, 1.0, 1.0)

    prefactor = ws.oft_domain_prefactor * ws.gamma_norm_factor
    if config.domain isa EnergyDomain
        inv_4sigma2 = 1.0 / (4 * config.sigma^2)
        QuantumFurnace._apply_lindbladian_threaded_energy!(
            sc, rho, ws.jump_eigenbases, ws.jump_hermitian,
            ham.bohr_freqs, ws.energy_labels, config, prefactor, inv_4sigma2;
            adjoint=adjoint)
    else
        nufft = ws.oft_nufft_prefactors
        QuantumFurnace._apply_lindbladian_threaded_timetrot!(
            sc, rho, ws.jump_eigenbases, ws.jump_hermitian,
            nufft.data, nufft.energy_to_index, ws.energy_labels, config, prefactor;
            adjoint=adjoint)
    end
    return copy(sc.rho_out)
end

@testset "Lindbladian threaded matvec: serial ≡ threaded" begin
    if Threads.nthreads() > 1
        rng = MersenneTwister(123)

        # Build a non-Hermitian complex jump for breadth (covers the
        # `is_herm == false` branch of the work-list builder). KMS detailed
        # balance in production needs the conjugate-paired partner, but a
        # single non-Hermitian jump is fine for a serial vs threaded unit
        # test — both paths run the same (possibly non-physical) physics.
        raw_complex = randn(rng, ComplexF64, DIM, DIM) ./ sqrt(DIM)
        complex_jump = JumpOp(raw_complex,
                              TEST_HAM.eigvecs' * raw_complex * TEST_HAM.eigvecs,
                              false, false)

        for (domain, jumps, trott, name) in [
            (EnergyDomain(),  TEST_JUMPS,                 nothing,        "Energy / Hermitian jumps"),
            (EnergyDomain(),  JumpOp[complex_jump],       nothing,        "Energy / non-Hermitian jump"),
            (TimeDomain(),    TEST_JUMPS,                 nothing,        "Time / Hermitian jumps"),
            (TimeDomain(),    JumpOp[complex_jump],       nothing,        "Time / non-Hermitian jump"),
            (TrotterDomain(), TEST_TROTTER_JUMPS,         TEST_TROTTER,   "Trotter / Hermitian jumps"),
        ]
            config = make_config(Lindbladian(), domain; construction=KMS())
            ws_a = Workspace(config, TEST_HAM, jumps; trotter=trott)
            ws_b = Workspace(config, TEST_HAM, jumps; trotter=trott)

            for adjoint in (false, true)
                Random.seed!(rng, 7)
                rho = Matrix(random_density_matrix(NUM_QUBITS))

                serial = if adjoint
                    copy(apply_adjoint_lindbladian!(ws_a, rho, config, TEST_HAM))
                else
                    copy(apply_lindbladian!(ws_a, rho, config, TEST_HAM))
                end
                threaded = _run_threaded_lindbladian!(ws_b, rho, config, TEST_HAM; adjoint=adjoint)

                # FP-accumulation tolerance scaled to ‖L(rho)‖. Per-jump
                # sandwich GEMMs keep the error in the 1e-13/‖L(rho)‖ regime
                # for any reduction order; the helpers above empirically land
                # at ~1e-16, ten orders below the gate.
                tol = max(norm(serial), 1.0) * 1e-12
                err = norm(threaded - serial)
                @test err < tol
                @info "Lindbladian threaded match" path=name adjoint=adjoint err=err tol=tol
            end
        end
    else
        @info "Skipping Lindbladian threaded matvec test (nthreads=$(Threads.nthreads()))"
        @test true
    end
end

@testset "Lindbladian threaded matvec: BLAS thread restoration" begin
    if Threads.nthreads() > 1
        config = make_config(Lindbladian(), EnergyDomain(); construction=KMS())
        ws = Workspace(config, TEST_HAM, TEST_JUMPS)
        rho = Matrix(random_density_matrix(NUM_QUBITS))

        old_blas = BLAS.get_num_threads()
        _run_threaded_lindbladian!(ws, rho, config, TEST_HAM; adjoint=false)
        @test BLAS.get_num_threads() == old_blas
        _run_threaded_lindbladian!(ws, rho, config, TEST_HAM; adjoint=true)
        @test BLAS.get_num_threads() == old_blas
        @info "Lindbladian threaded BLAS restoration" blas_before=old_blas blas_after=BLAS.get_num_threads()
    else
        @info "Skipping Lindbladian threaded BLAS restoration test (nthreads=$(Threads.nthreads()))"
        @test true
    end
end

@testset "Lindbladian threaded matvec: empty work list short-circuit" begin
    if Threads.nthreads() > 1
        # An all-Hermitian jump set with energy_labels that are all > 1e-12
        # produces an empty work list; the helper must return without
        # touching `sc.rho_out` beyond the coherent terms already there.
        config = make_config(Lindbladian(), EnergyDomain(); construction=KMS())
        ws = Workspace(config, TEST_HAM, TEST_JUMPS)
        rho = Matrix(random_density_matrix(NUM_QUBITS))

        sc = ws.scratch
        prefactor = ws.oft_domain_prefactor * ws.gamma_norm_factor
        inv_4sigma2 = 1.0 / (4 * config.sigma^2)

        # Pre-fill rho_out with a sentinel value
        fill!(sc.rho_out, ComplexF64(7.0))
        # Empty energy labels -> empty work list
        QuantumFurnace._apply_lindbladian_threaded_energy!(
            sc, rho, ws.jump_eigenbases, ws.jump_hermitian,
            TEST_HAM.bohr_freqs, Float64[], config, prefactor, inv_4sigma2;
            adjoint=false)
        @test all(sc.rho_out .== ComplexF64(7.0))
    else
        @info "Skipping empty work-list test (nthreads=$(Threads.nthreads()))"
        @test true
    end
end

# ============================================================================
# Channel-Krylov ω-loop threading tests (qf-in3 follow-up)
# Mirrors the Lindbladian threading tests above for `apply_delta_channel!`.
#
# qf-po5 Commit 2 deleted the `_run_threaded_channel!` helper and its
# "Channel threaded matvec: serial ≡ threaded" testset. The new faithful
# `apply_delta_channel!` consumes `_accumulate_rho_jump_threaded_*!` (whose
# serial ≡ threaded equivalence is regressioned by the per-step run_thermalize
# tests in `test_thermalization.jl` and the byte-identity check in
# `test_predict_channel.jl` (a)/(b1)) — driving the threaded variants directly
# from a re-implementation of the deleted summed-channel matvec is no longer
# meaningful. The BLAS thread-restoration testset below stays.
# ============================================================================

@testset "Channel threaded matvec: BLAS thread restoration" begin
    if Threads.nthreads() > 1
        config = make_config(Thermalize(), EnergyDomain(); construction=KMS())
        ws = Workspace(config, TEST_HAM, TEST_JUMPS)
        rho = Matrix(random_density_matrix(NUM_QUBITS))

        old_blas = BLAS.get_num_threads()
        apply_delta_channel!(ws, rho, config, TEST_HAM)
        @test BLAS.get_num_threads() == old_blas
        @info "Channel threaded BLAS restoration" blas_before=old_blas blas_after=BLAS.get_num_threads()
    else
        @info "Skipping channel threaded BLAS restoration test (nthreads=$(Threads.nthreads()))"
        @test true
    end
end
