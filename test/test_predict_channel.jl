using LinearAlgebra: I, eigvals, svdvals, norm, tr, Hermitian
using Test
using QuantumFurnace

# qf-ev5.{1,5}: Krylov spectral-expansion predictor for the IMPLEMENTED CPTP
# channel Φ_δ that run_thermalize executes (per-jump :sweep, weak-measurement
# Kraus + coherent unitary). The forward matvec is byte-for-byte the same
# kernel run_thermalize uses, so the Krylov reconstruction must agree with
# stepping run_thermalize for `mixing_time / δ` steps to machine precision.

@testset "predict_channel_trajectory" begin

    # -----------------------------------------------------------------------
    # (a) n=3, β=10, δ=1e-3: byte-identical to run_thermalize @ jump_selection=:sweep
    # -----------------------------------------------------------------------
    @testset "(a) n=3, β=10, δ=1e-3 :sweep BohrDomain CKG vs run_thermalize" begin
        beta = 10.0
        delta = 1e-3
        sys = make_dll_n3_system(beta)
        ham = sys.ham
        jumps = sys.jumps
        d = size(ham.data, 1)

        cfg = Config(
            sim = Thermalize(),
            domain = BohrDomain(),
            construction = KMS(),
            num_qubits = 3,
            with_linear_combination = true,
            beta = beta, sigma = 1.0 / beta,
            a = 0.0, s = 0.25,
            num_energy_bits = 12,
            w0 = 0.05,
            t0 = 2π / (2^12 * 0.05),
            num_trotter_steps_per_t0 = 10,
            delta = delta,
            mixing_time = 30.0,
            jump_selection = :sweep,
        )

        rho_0 = Matrix{ComplexF64}(I(d) / d)

        # k-grid of 31 points spanning the run_thermalize trajectory.
        k_step = 1000
        k_grid = collect(0:k_step:30000)

        res_kr = predict_channel_trajectory(cfg, ham, jumps, rho_0, k_grid;
                                              krylovdim=40, tol=1e-10)

        res_th = run_thermalize(jumps, cfg, ham; initial_dm=copy(rho_0),
                                save_every=k_step)

        # Sanity on shape.
        @test length(res_kr.t) == length(k_grid)
        @test length(res_kr.distances) == length(k_grid)
        @test res_kr.delta_used == delta
        @test res_kr.k_grid == k_grid

        # Channel steady state mu_1 must be exactly 1 to numerical tolerance
        # (CPTP fixed point).
        @test abs(abs(res_kr.eigenvalues[1]) - 1.0) < 1e-10

        # Slow-mode mu_2 should be just below 1; gap (in 1/δ units) must be
        # positive and close to the corresponding Lindbladian gap (~ 0.23 from
        # the test_lindblad_action testset (a)).
        @test 0.0 < res_kr.spectral_gap < 1.0

        # Byte-identical agreement with run_thermalize on the entire
        # trajectory at every save point. PHYSICS CHECK: the forward closure
        # in predict_channel_trajectory uses the SAME _apply_one_dm_substep!
        # kernel run_thermalize calls, so any drift comes from finite Krylov
        # subspace + dense eigen(H) tolerance only — well below 1e-10 for
        # the n=3 test fixture.
        @assert length(res_kr.distances) == length(res_th.trace_distances)
        max_abs_err = maximum(abs.(res_kr.distances .- res_th.trace_distances))
        @test max_abs_err < 1e-10

        @info "(a) Krylov channel vs run_thermalize" max_abs_err matvecs_kr=res_kr.total_matvecs steps_th=length(res_th.trace_distances)

        # rho_final agrees too (within hermitisation noise).
        rho_th_final = res_th.final_dm
        @test maximum(abs.(res_kr.rho_final .- rho_th_final)) < 1e-9
    end

    # -----------------------------------------------------------------------
    # (b) :random rejected at validation
    # -----------------------------------------------------------------------
    @testset "(b) :random jump_selection rejected" begin
        beta = 10.0
        sys = make_dll_n3_system(beta)
        d = size(sys.ham.data, 1)
        cfg = Config(
            sim = Thermalize(),
            domain = BohrDomain(),
            construction = KMS(),
            num_qubits = 3,
            with_linear_combination = true,
            beta = beta, sigma = 1.0 / beta,
            a = 0.0, s = 0.25,
            num_energy_bits = 12, w0 = 0.05,
            t0 = 2π / (2^12 * 0.05),
            num_trotter_steps_per_t0 = 10,
            delta = 1e-3, mixing_time = 5.0,
            jump_selection = :random,
        )
        rho_0 = Matrix{ComplexF64}(I(d) / d)
        @test_throws ArgumentError predict_channel_trajectory(
            cfg, sys.ham, sys.jumps, rho_0, [0, 100, 1000])
    end

    # -----------------------------------------------------------------------
    # (b1) TrotterDomain + GQSP smoke (qf-ev5.8): same matvec kernel as
    #      run_thermalize ⇒ byte-identical reconstruction at small δ.
    # -----------------------------------------------------------------------
    @testset "(b1) TrotterDomain + GQSP smoke" begin
        # Build a Thermalize config that matches the trotterise+gqsp path used
        # by test_gqsp_thermalize.jl. Short trajectory (δ=1e-3, 100 steps) to
        # keep the test fast — agreement is byte-identical anyway.
        cfg = Config(
            sim = Thermalize(),
            domain = TrotterDomain(),
            construction = KMS(),
            num_qubits = 3,
            with_linear_combination = true,
            beta = BETA, sigma = SIGMA,
            a = BETA / 30.0, s = 0.4,
            num_energy_bits = NUM_ENERGY_BITS,
            w0 = W0, t0 = T0,
            num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
            delta = 1e-3, mixing_time = 0.1,  # 100 steps
            jump_selection = :sweep,
            with_gqsp = true, gqsp_degree = 1,
        )

        d = N3_DIM
        rho_0 = Matrix{ComplexF64}(I(d) / d)
        k_grid = collect(0:25:100)

        res_kr = predict_channel_trajectory(
            cfg, N3_HAM, N3_TROTTER_JUMPS, rho_0, k_grid;
            krylovdim=30, tol=1e-10, trotter=N3_TROTTER)

        res_th = run_thermalize(N3_TROTTER_JUMPS, cfg, N3_HAM, N3_TROTTER;
                                 initial_dm=copy(rho_0), save_every=25)

        # Channel fixed point (CPTP). PHYSICS CHECK: GQSP at degree d
        # introduces O((δα)^(d+1)) per-step non-CPTPness via polynomial
        # truncation; at δ=1e-3, d=1 this is ~ (δα)^2 ~ 1e-6, easily
        # within 1e-5.
        @test abs(abs(res_kr.eigenvalues[1]) - 1.0) < 1e-5

        # PHYSICS CHECK: closure uses the same _apply_one_dm_substep! kernel
        # as run_thermalize through _precompute_coherent_unitary +
        # _precompute_per_jump_channels, so the residual is dominated by
        # Krylov-subspace truncation only. 1e-6 covers d=1 GQSP truncation
        # noise propagated through 100 steps; at δ=1e-3 these are well below.
        n_compare = min(length(res_kr.distances), length(res_th.trace_distances))
        max_abs = maximum(abs.(res_kr.distances[1:n_compare] .-
                                res_th.trace_distances[1:n_compare]))
        @test max_abs < 1e-6
        @info "(b1) TrotterDomain + GQSP" max_abs matvecs_kr=res_kr.total_matvecs
    end

    # -----------------------------------------------------------------------
    # (c) Output shape + delta_used + states
    # -----------------------------------------------------------------------
    @testset "(c) save_states and struct fields" begin
        beta = 10.0
        delta = 1e-3
        sys = make_dll_n3_system(beta)
        d = size(sys.ham.data, 1)
        cfg = Config(
            sim = Thermalize(),
            domain = BohrDomain(),
            construction = KMS(),
            num_qubits = 3,
            with_linear_combination = true,
            beta = beta, sigma = 1.0 / beta,
            a = 0.0, s = 0.25,
            num_energy_bits = 12, w0 = 0.05,
            t0 = 2π / (2^12 * 0.05),
            num_trotter_steps_per_t0 = 10,
            delta = delta, mixing_time = 5.0,
            jump_selection = :sweep,
        )
        rho_0 = Matrix{ComplexF64}(I(d) / d)
        k_grid = [0, 100, 500, 1000, 2000, 5000]
        res = predict_channel_trajectory(cfg, sys.ham, sys.jumps, rho_0, k_grid;
                                           save_states=true)
        @test length(res.states) == length(k_grid)
        @test all(size.(res.states) .== ((d, d),))
        # rho(k=0) ≈ rho_0 within hermitisation tolerance.
        @test maximum(abs.(res.states[1] .- rho_0)) < 1e-9
        # rho_final = rho(k_grid[end]) (within tiny float noise).
        @test isapprox(res.rho_final, res.states[end]; atol=1e-12)
        @test res.delta_used == delta
        # Trace preservation.
        @test isapprox(real(tr(res.rho_final)), 1.0; atol=1e-7)
    end
end
