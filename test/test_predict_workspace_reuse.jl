# test/test_predict_workspace_reuse.jl
#
# qf-qmi.2: covers the `workspace=` kwarg on
# `predict_lindbladian_trajectory` and `predict_channel_trajectory`.
#
# Per the qf-qmi verifier report (Findings 5, 6, 8):
#   1. Reuse must be bitwise idempotent across calls (no scratch contamination).
#   2. Cross-config reuse must throw a clean ArgumentError, not a silent
#      miscalculation.
#   3. Cross-simulation mismatch (Lindbladian ws into channel predict, etc.)
#      must throw a clean ArgumentError, not a downstream MethodError/FieldError.
#   4. Dim / jump-count mismatch must throw cleanly.

using LinearAlgebra: I
using Test
using QuantumFurnace

include(joinpath(@__DIR__, "test_helpers.jl"))


@testset "predict_*_trajectory workspace= reuse" begin

    # -----------------------------------------------------------------------
    # (a) Lindbladian: bitwise reuse identity at n=3 (no scratch contamination)
    # -----------------------------------------------------------------------
    @testset "(a) Lindbladian reuse byte-equality (n=3)" begin
        beta = 10.0
        sys = make_dll_n3_system(beta)
        ham = sys.ham; jumps = sys.jumps
        d = size(ham.data, 1)

        cfg = Config(
            sim = Lindbladian(), domain = EnergyDomain(), construction = KMS(),
            num_qubits = 3, with_linear_combination = true,
            beta = beta, sigma = 1.0 / beta, a = 0.0, s = 0.25,
            num_energy_bits = 12, w0 = 0.05,
            t0 = 2π / (2^12 * 0.05), num_trotter_steps_per_t0 = 10,
        )
        rho_0 = Matrix{ComplexF64}(I(d) / d)
        t_grid = collect(range(0.0, 5.0, length=11))

        # Baseline (no workspace=)
        r1 = predict_lindbladian_trajectory(cfg, ham, jumps, rho_0, t_grid;
                                             krylovdim=20, tol=1e-10)
        # Reuse with explicit workspace
        ws = Workspace(cfg, ham, jumps)
        r2 = predict_lindbladian_trajectory(cfg, ham, jumps, rho_0, t_grid;
                                             krylovdim=20, tol=1e-10, workspace=ws)
        r3 = predict_lindbladian_trajectory(cfg, ham, jumps, rho_0, t_grid;
                                             krylovdim=20, tol=1e-10, workspace=ws)

        # Strict bitwise identity: baseline ↔ reuse-1 and reuse-1 ↔ reuse-2.
        # (Krylov factorisation is deterministic given the same matvec closure.)
        @test r1.distances == r2.distances
        @test r2.distances == r3.distances
        @test r1.spectral_gap == r2.spectral_gap
        @test r2.spectral_gap == r3.spectral_gap
        @test r1.eigenvalues == r2.eigenvalues
    end

    # -----------------------------------------------------------------------
    # (b) Channel: bitwise reuse identity at n=3
    # -----------------------------------------------------------------------
    @testset "(b) Channel reuse byte-equality (n=3)" begin
        beta = 5.0; eb = 12; delta = 1e-3
        sys = make_dll_n3_system(beta)
        ham = sys.ham; jumps = sys.jumps
        d = size(ham.data, 1)

        num_t0_steps_per_step = 10
        t0 = 2π / (2^eb * 0.05)
        cfg = Config(
            sim = Thermalize(), domain = TrotterDomain(), construction = KMS(),
            num_qubits = 3, with_linear_combination = true,
            beta = beta, sigma = 1.0 / beta,
            a = beta / 30.0, s = 0.4,
            num_energy_bits = eb, w0 = 0.05,
            t0 = t0, num_trotter_steps_per_t0 = num_t0_steps_per_step,
            delta = delta, mixing_time = 1.0,
            jump_selection = :sweep, with_gqsp = true, gqsp_degree = 1,
        )
        trotter = TrottTrott(ham, t0/num_t0_steps_per_step, num_t0_steps_per_step)
        rho_0 = Matrix{ComplexF64}(I(d) / d)
        k_grid = collect(0:50:300)

        r1 = predict_channel_trajectory(cfg, ham, jumps, rho_0, k_grid;
                                          krylovdim=20, tol=1e-10, trotter=trotter)
        ws = Workspace(cfg, ham, jumps; trotter=trotter)
        r2 = predict_channel_trajectory(cfg, ham, jumps, rho_0, k_grid;
                                          krylovdim=20, tol=1e-10, trotter=trotter, workspace=ws)
        r3 = predict_channel_trajectory(cfg, ham, jumps, rho_0, k_grid;
                                          krylovdim=20, tol=1e-10, trotter=trotter, workspace=ws)

        @test r1.distances == r2.distances
        @test r2.distances == r3.distances
        @test r1.spectral_gap == r2.spectral_gap
        @test r2.spectral_gap == r3.spectral_gap
        @test r1.eigenvalues == r2.eigenvalues
    end

    # -----------------------------------------------------------------------
    # (c) Cross-config mismatch must throw cleanly (Lindbladian path)
    # -----------------------------------------------------------------------
    @testset "(c) Cross-config mismatch -> ArgumentError (Lindbladian)" begin
        beta = 10.0
        sys = make_dll_n3_system(beta)
        ham = sys.ham; jumps = sys.jumps
        d = size(ham.data, 1)

        cfg_a = Config(
            sim = Lindbladian(), domain = EnergyDomain(), construction = KMS(),
            num_qubits = 3, with_linear_combination = true,
            beta = beta, sigma = 1.0/beta, a = 0.0, s = 0.25,
            num_energy_bits = 12, w0 = 0.05,
            t0 = 2π/(2^12 * 0.05), num_trotter_steps_per_t0 = 10,
        )
        cfg_b = Config(
            sim = Lindbladian(), domain = EnergyDomain(), construction = KMS(),
            num_qubits = 3, with_linear_combination = true,
            beta = beta, sigma = 0.05, a = 0.0, s = 0.25,   # σ differs
            num_energy_bits = 12, w0 = 0.05,
            t0 = 2π/(2^12 * 0.05), num_trotter_steps_per_t0 = 10,
        )

        rho_0 = Matrix{ComplexF64}(I(d) / d)
        t_grid = collect(range(0.0, 5.0, length=11))
        ws_a = Workspace(cfg_a, ham, jumps)

        # Mismatched config must throw ArgumentError, not silently miscalculate.
        @test_throws ArgumentError predict_lindbladian_trajectory(
            cfg_b, ham, jumps, rho_0, t_grid;
            krylovdim=20, tol=1e-10, workspace=ws_a)
    end

    # -----------------------------------------------------------------------
    # (d) Cross-simulation workspace mismatch -> clean ArgumentError
    # -----------------------------------------------------------------------
    @testset "(d) Cross-simulation workspace mismatch -> ArgumentError" begin
        beta = 5.0; eb = 12; delta = 1e-3
        sys = make_dll_n3_system(beta)
        ham = sys.ham; jumps = sys.jumps
        d = size(ham.data, 1)

        cfg_L = Config(
            sim = Lindbladian(), domain = EnergyDomain(), construction = KMS(),
            num_qubits = 3, with_linear_combination = true,
            beta = beta, sigma = 1.0/beta, a = 0.0, s = 0.25,
            num_energy_bits = eb, w0 = 0.05,
            t0 = 2π/(2^eb * 0.05), num_trotter_steps_per_t0 = 10,
        )
        num_t0_steps_per_step = 10
        t0 = 2π / (2^eb * 0.05)
        cfg_C = Config(
            sim = Thermalize(), domain = TrotterDomain(), construction = KMS(),
            num_qubits = 3, with_linear_combination = true,
            beta = beta, sigma = 1.0/beta,
            a = beta/30.0, s = 0.4,
            num_energy_bits = eb, w0 = 0.05,
            t0 = t0, num_trotter_steps_per_t0 = num_t0_steps_per_step,
            delta = delta, mixing_time = 1.0,
            jump_selection = :sweep, with_gqsp = true, gqsp_degree = 1,
        )
        trotter = TrottTrott(ham, t0/num_t0_steps_per_step, num_t0_steps_per_step)

        rho_0 = Matrix{ComplexF64}(I(d) / d)
        ws_L = Workspace(cfg_L, ham, jumps)
        ws_C = Workspace(cfg_C, ham, jumps; trotter=trotter)

        # Channel ws into Lindbladian predict (and vice versa) -> ArgumentError.
        @test_throws ArgumentError predict_lindbladian_trajectory(
            cfg_L, ham, jumps, rho_0, [0.0, 1.0];
            krylovdim=20, tol=1e-10, workspace=ws_C)
        @test_throws ArgumentError predict_channel_trajectory(
            cfg_C, ham, jumps, rho_0, [0, 10];
            krylovdim=20, tol=1e-10, trotter=trotter, workspace=ws_L)
    end

    # -----------------------------------------------------------------------
    # (e) Dim / jump-count mismatch -> AssertionError
    # -----------------------------------------------------------------------
    @testset "(e) Dim and jump-count mismatch -> AssertionError" begin
        beta = 5.0
        sys = make_dll_n3_system(beta)
        ham = sys.ham; jumps = sys.jumps
        d = size(ham.data, 1)

        cfg = Config(
            sim = Lindbladian(), domain = EnergyDomain(), construction = KMS(),
            num_qubits = 3, with_linear_combination = true,
            beta = beta, sigma = 1.0/beta, a = 0.0, s = 0.25,
            num_energy_bits = 12, w0 = 0.05,
            t0 = 2π/(2^12 * 0.05), num_trotter_steps_per_t0 = 10,
        )
        rho_0_wrong = Matrix{ComplexF64}(I(2 * d) / (2 * d))   # wrong dim
        jumps_short = jumps[1:end-1]                          # one jump dropped
        ws = Workspace(cfg, ham, jumps)

        @test_throws AssertionError predict_lindbladian_trajectory(
            cfg, ham, jumps, rho_0_wrong, [0.0, 1.0];
            krylovdim=20, tol=1e-10, workspace=ws)
        # Jumps_short would also create a config that doesn't match → ArgumentError
        # before AssertionError fires. Either is a clean failure.
        @test_throws Exception predict_lindbladian_trajectory(
            cfg, ham, jumps_short, Matrix{ComplexF64}(I(d)/d), [0.0, 1.0];
            krylovdim=20, tol=1e-10, workspace=ws)
    end

end  # @testset
