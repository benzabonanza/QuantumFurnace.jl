# test/test_predict_sandbox.jl
#
# Sandbox shadow of test_predict_lindbladian.jl + test_predict_channel.jl
# (qf-x56.3). Two headline NO_SANDBOX invariants from qf-ev5:
#
#   1. predict_lindbladian_trajectory matches the dense reference
#      exp(t·L)·vec(ρ₀) along the trajectory — the spectral-expansion
#      accuracy invariant (Krylov subspace captures the slow modes).
#   2. predict_channel_trajectory matches run_thermalize byte-for-byte at
#      small δ — the channel byte-identity invariant (closure uses the same
#      _apply_one_dm_substep! kernel as run_thermalize).
#
# Both at n = 3, β = 10. The Lindbladian check uses 5 time points (range
# [0, 5/gap]) at krylovdim = 30 — captures the diagonal-sector dynamics on
# d² = 64 modes ⇒ 1e-7 absolute error matches the heavy test threshold.
# The channel check uses ~50 outer steps at δ = 1e-3 (mixing_time = 0.05)
# at krylovdim = 30 — matches run_thermalize to 1e-10 in trace distance
# along the full save_every-spaced grid.

using LinearAlgebra: I, eigvals, svdvals, norm, tr, Hermitian
using Test
using QuantumFurnace


@testset "Predictor sandbox shadows (qf-x56.3)" begin

    # -----------------------------------------------------------------------
    # (a) predict_lindbladian_trajectory accuracy vs dense reference.
    # PHYSICS CHECK: at n=3 with the 3n=9 single-Pauli jump set + KMS-DB
    # the slow spectrum is well-separated; krylovdim=30 captures the entire
    # diagonal-sector dynamics on d²=64 modes. Same threshold as the heavy
    # test (1e-7) — the regime where bi-exp τ_mix extraction is accurate.
    # -----------------------------------------------------------------------
    @testset "(a) predict_lindbladian vs dense @ n=3, β=10" begin
        beta = 10.0
        sys = make_dll_n3_system(beta)
        ham = sys.ham; jumps = sys.jumps
        d = size(ham.data, 1)

        cfg = Config(
            sim = Lindbladian(),
            domain = EnergyDomain(),
            construction = KMS(),
            num_qubits = 3,
            with_linear_combination = true,
            beta = beta, sigma = 1.0 / beta,
            a = 0.0, s = 0.25,
            num_energy_bits = 12, w0 = 0.05,
            t0 = 2π / (2^12 * 0.05),
            num_trotter_steps_per_t0 = 10,
        )

        rho_0 = Matrix{ComplexF64}(I(d) / d)
        sigma_beta = Matrix{ComplexF64}(ham.gibbs)

        # Use Krylov gap to pick t_max (same as heavy test).
        gap_res = krylov_spectral_gap(cfg, ham, jumps;
                                       krylovdim = 30, howmany = 4, tol = 1e-10)
        gap = gap_res.spectral_gap
        @test gap > 0

        # qf-6yw: Pass-2 (krylov_spectral_gap) carries operator-side diagnostics
        # (no seeded ρ₀ ⇒ c-side is NaN; R-side is the ρ₀-independent picture).
        gsm = gap_res.spectral_modes
        @test gsm isa SpectralModeDiagnostics
        @test all(isnan, gsm.c_abs2)
        @test all(0.0 .<= gsm.off_diag_weight .<= 1.0)
        @test gsm.off_diag_weight[1] < 1e-6   # fixed_point ≈ σ_β: diagonal

        # 5 time points across [0, 5/gap] — shrunk from heavy test's 21 to
        # keep dense exp(t·L) workload trivial (5 exponentiations of a
        # 64-dim generator).
        t_grid = collect(range(0.0, 5.0 / gap, length = 5))

        # Dense reference trajectory.
        L_dense = construct_lindbladian(jumps, cfg, ham)
        v0 = vec(rho_0)
        distances_dense = Vector{Float64}(undef, length(t_grid))
        for (k, t) in enumerate(t_grid)
            rho_t = reshape(exp(t * L_dense) * v0, d, d)
            rho_t .= (rho_t + rho_t') ./ 2
            distances_dense[k] = sum(svdvals(rho_t .- sigma_beta)) / 2
        end

        # Krylov spectral expansion.
        res_kr = predict_lindbladian_trajectory(cfg, ham, jumps, rho_0, t_grid;
                                                 krylovdim = 30, tol = 1e-10)
        @test length(res_kr.t) == length(t_grid)
        @test size(res_kr.rho_final) == (d, d)
        @test res_kr.total_matvecs <= length(t_grid) * 30

        # qf-6yw: per-mode spectral diagnostics are attached to every result.
        sm = res_kr.spectral_modes
        @test sm isa SpectralModeDiagnostics
        @test length(sm.off_diag_weight) == length(res_kr.eigenvalues)
        @test all(0.0 .<= sm.off_diag_weight .<= 1.0)
        @test sm.off_diag_weight[1] < 1e-6   # steady mode ≈ σ_β: diagonal in the energy eigenbasis

        max_abs_err = maximum(abs.(res_kr.distances .- distances_dense))
        @test max_abs_err < 1e-7
        @info "(a) predict_lindbladian sandbox" max_abs_err matvecs=res_kr.total_matvecs

        # Spectral gap from Krylov must match dense (KrylovKit tolerance).
        eigs_dense = sort(eigvals(L_dense); by = v -> abs(real(v)))
        gap_dense = abs(real(eigs_dense[2]))
        @test isapprox(res_kr.spectral_gap, gap_dense; rtol = 1e-9)
    end

    # -----------------------------------------------------------------------
    # (b) predict_channel_trajectory byte-identity vs run_thermalize.
    # 50 outer steps at δ=1e-3 (mixing_time=0.05). The forward closure uses
    # the SAME _apply_one_dm_substep! kernel run_thermalize calls, so the
    # residual is dominated by finite Krylov subspace + dense eigen(H)
    # tolerance only — must be below 1e-10 for the n=3 fixture.
    # -----------------------------------------------------------------------
    @testset "(b) predict_channel byte-identity vs run_thermalize @ n=3, δ=1e-3" begin
        beta = 10.0
        delta = 1e-3
        sys = make_dll_n3_system(beta)
        ham = sys.ham; jumps = sys.jumps
        d = size(ham.data, 1)

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
            delta = delta,
            mixing_time = 0.05,        # 50 outer steps
            jump_selection = :sweep,
        )
        rho_0 = Matrix{ComplexF64}(I(d) / d)

        k_step = 10
        k_grid = collect(0:k_step:50)   # [0, 10, …, 50] — 6 save points

        res_kr = predict_channel_trajectory(cfg, ham, jumps, rho_0, k_grid;
                                              krylovdim = 30, tol = 1e-10)
        res_th = run_thermalize(jumps, cfg, ham; initial_dm = copy(rho_0),
                                save_every = k_step)

        @test length(res_kr.t) == length(k_grid)
        @test res_kr.delta_used == delta
        @test res_kr.k_grid == k_grid

        # CPTP fixed point at mu_1.
        @test abs(abs(res_kr.eigenvalues[1]) - 1.0) < 1e-10

        @assert length(res_kr.distances) == length(res_th.trace_distances)
        max_abs_err = maximum(abs.(res_kr.distances .- res_th.trace_distances))
        @test max_abs_err < 1e-10
        @info "(b) predict_channel sandbox byte-identity" max_abs_err matvecs=res_kr.total_matvecs

        # Final density matrix agrees too (within hermitisation noise).
        rho_th_final = res_th.final_dm
        @test maximum(abs.(res_kr.rho_final .- rho_th_final)) < 1e-9
    end
end
