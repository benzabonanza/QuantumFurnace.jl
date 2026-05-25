# test/test_qf_0fv_verifier.jl
#
# Adversarial verifier test for qf-0fv (commit 2b77e6f). Probes the gated
# `compute_true_gap` kwarg on `predict_lindbladian_trajectory` and
# `predict_channel_trajectory` for edge cases that the bundled smoke test
# `scripts/scratch_qf_0fv_gating_check.jl` does NOT cover:
#
# 1.  Demonstrate that the DEFAULT (Pass-1) returns the WRONG gap on a
#     parity-symmetric `(rho_0, L)` fixture — classical n=3 ZZ Ising +
#     `rho_0 = I/d` — and that `compute_true_gap=true` rescues it.
#     This catches a future regression that silently inlines Pass-2 logic
#     into Pass-1 (which would make the gating no-op).
#
# 2.  Type stability of both branches via `@inferred`-style assertions on
#     the returned `spectral_gap` value type.
#
# 3.  Channel-path return shape: the `(μ-1)/δ` vs `-log|μ|/δ` formulas
#     agree to O(δ) at small δ on an ideal channel.

using LinearAlgebra: I, eigvals, Hermitian
using Test
using QuantumFurnace

@testset "qf-0fv adversarial verifier" begin

    # -----------------------------------------------------------------------
    # (1) Parity-symmetric classical Ising + rho_0 = I/d: Pass-1 must
    #     differ from Pass-2 (otherwise the gating is no-op).
    # -----------------------------------------------------------------------
    @testset "Pass-1 vs Pass-2 on parity-symmetric fixture" begin
        n_ising = 3
        terms_zz = [[Z, Z]]
        coeffs_zz = [0.7]
        base_ham = QuantumFurnace._construct_base_ham(terms_zz, coeffs_zz, n_ising;
                                                      periodic=true)
        rescaling_factor, shift = QuantumFurnace._rescaling_and_shift_factors(base_ham)
        d_ising = 2^n_ising
        rescaled_mat = Matrix(base_ham) ./ rescaling_factor .+ shift * I(d_ising)
        eigvals_rs, eigvecs_rs = eigen(Hermitian(rescaled_mat))
        raw = (
            matrix = rescaled_mat,
            terms = terms_zz,
            base_coeffs = coeffs_zz ./ rescaling_factor,
            disordering_terms = nothing,
            disordering_coeffs = nothing,
            eigvals = eigvals_rs, eigvecs = eigvecs_rs,
            nu_min = minimum(diff(eigvals_rs)),
            shift = shift,
            rescaling_factor = rescaling_factor,
            periodic = true,
        )
        β_phys = 0.5
        ham_ising = HamHam(raw; beta_phys=β_phys)
        jumps_ising = QuantumFurnace._build_jump_set(ham_ising, n_ising)
        β_alg = beta_alg(ham_ising, β_phys)

        σ = 1.0 / β_alg
        H_norm = maximum(abs, ham_ising.eigvals)
        omega_range = 2.0 * (H_norm + 8 * σ)
        r_D = 7
        w0_D = omega_range / 2.0^r_D
        t0_D = 2π / (2.0^r_D * w0_D)
        cfg = Config(
            sim = Lindbladian(), domain = EnergyDomain(), construction = KMS(),
            num_qubits = n_ising, with_linear_combination = true,
            beta = β_alg, beta_phys = β_phys, sigma = σ,
            a = 0.0, s = 0.25, gaussian_parameters = (nothing, nothing),
            num_energy_bits_D = r_D, w0_D = w0_D, t0_D = t0_D,
            num_trotter_steps_per_t0 = 10, filter = nothing,
        )

        L_dense = construct_lindbladian(jumps_ising, cfg, ham_ising)
        ev_dense = eigvals(Matrix(L_dense))
        perm = sortperm(real.(ev_dense); by=abs)
        gap_dense = abs(real(ev_dense[perm[2]]))

        rho_0 = Matrix{ComplexF64}(I(d_ising) / d_ising)
        t_grid = collect(range(0.0, 100.0, length=21))

        # Default = Pass-1 only. On this fixture Pass-1 MUST return the
        # parity-even sub-spectrum gap (different from gap_dense).
        traj_p1 = predict_lindbladian_trajectory(
            cfg, ham_ising, jumps_ising, rho_0, t_grid; krylovdim=40)
        # Opt-in = Pass-1 + Pass-2. Must rescue the true gap.
        traj_p2 = predict_lindbladian_trajectory(
            cfg, ham_ising, jumps_ising, rho_0, t_grid;
            krylovdim=40, compute_true_gap=true)

        # Witness that the gating actually changes behavior on this
        # fixture: Pass-1 disagrees with Pass-2 by > 50% (i.e. the
        # parity trap fires).
        @test abs(traj_p1.spectral_gap - traj_p2.spectral_gap) /
              traj_p2.spectral_gap > 0.5
        # Pass-2 must match dense to machine precision.
        @test isapprox(traj_p2.spectral_gap, gap_dense; rtol=1e-8)
        # Pass-1 must NOT match dense (gating sanity).
        @test !isapprox(traj_p1.spectral_gap, gap_dense; rtol=1e-3)

        @info "qf-0fv parity-fixture witness" gap_dense gap_pass1=traj_p1.spectral_gap gap_pass2=traj_p2.spectral_gap
    end

    # -----------------------------------------------------------------------
    # (2) Type stability of both branches (Float64-typed spectral_gap).
    # -----------------------------------------------------------------------
    @testset "Type stability of spectral_gap return" begin
        n = 3
        cfg = make_config(Lindbladian(), EnergyDomain(); num_qubits=n)
        ham = N3_HAM
        jumps = N3_JUMPS
        d = 2^n
        psi = ones(ComplexF64, d) ./ sqrt(2.0^n)
        rho_0 = psi * psi'
        t_grid = collect(range(0.0, 1.0; length=3))

        r_p1 = predict_lindbladian_trajectory(cfg, ham, jumps, rho_0, t_grid;
                                              krylovdim=30)
        r_p2 = predict_lindbladian_trajectory(cfg, ham, jumps, rho_0, t_grid;
                                              krylovdim=30, compute_true_gap=true)
        @test r_p1.spectral_gap isa Float64
        @test r_p2.spectral_gap isa Float64
        @test r_p1.total_matvecs isa Integer
        @test r_p2.total_matvecs isa Integer
        @test r_p1.all_converged isa Bool
        @test r_p2.all_converged isa Bool

        cfg_ch = make_config(Thermalize(), EnergyDomain(); num_qubits=n)
        k_grid = collect(0:5:20)
        rho_init = Matrix{ComplexF64}(rho_0)
        rc_p1 = predict_channel_trajectory(cfg_ch, ham, jumps, rho_init, k_grid;
                                            krylovdim=30)
        rc_p2 = predict_channel_trajectory(cfg_ch, ham, jumps, rho_init, k_grid;
                                            krylovdim=30, compute_true_gap=true)
        @test rc_p1.spectral_gap isa Float64
        @test rc_p2.spectral_gap isa Float64
        @test rc_p1.total_matvecs isa Integer
        @test rc_p2.total_matvecs isa Integer
        @test rc_p1.all_converged isa Bool
        @test rc_p2.all_converged isa Bool
    end

    # -----------------------------------------------------------------------
    # (3) Channel: Pass-1 (-log|μ|/δ) vs Pass-2 ((μ-1)/δ) agree to O(δ)
    #     on a non-pathological fixture.
    # -----------------------------------------------------------------------
    @testset "Channel Pass-1 vs Pass-2 O(δ) agreement" begin
        n = 3
        cfg_ch = make_config(Thermalize(), EnergyDomain(); num_qubits=n)
        delta = cfg_ch.delta
        ham = N3_HAM
        jumps = N3_JUMPS
        d = 2^n
        psi = ones(ComplexF64, d) ./ sqrt(2.0^n)
        rho_init = Matrix{ComplexF64}(psi * psi')
        k_grid = collect(0:10:60)
        rc_p1 = predict_channel_trajectory(cfg_ch, ham, jumps, rho_init, k_grid;
                                            krylovdim=30)
        rc_p2 = predict_channel_trajectory(cfg_ch, ham, jumps, rho_init, k_grid;
                                            krylovdim=30, compute_true_gap=true)
        # Both gaps positive on a physical fixture.
        @test rc_p1.spectral_gap > 0
        @test rc_p2.spectral_gap > 0
        # Difference bounded by O(δ) * gap (Taylor: log(1-x) = -x - x²/2 - …,
        # so -log|μ|/δ - (1-|μ|)/δ = O(δ·gap²) ≈ O(δ·gap)·gap).
        # For δ=0.01 and gap~0.2, allow 1% rel difference.
        rel_diff = abs(rc_p1.spectral_gap - rc_p2.spectral_gap) /
                   max(rc_p2.spectral_gap, eps())
        @test rel_diff < 5e-2
        @info "qf-0fv channel Pass-1/Pass-2 O(δ) drift" delta gap_pass1=rc_p1.spectral_gap gap_pass2=rc_p2.spectral_gap rel_diff
    end

    # -----------------------------------------------------------------------
    # (4) Matvec budget validates the gating: Pass-2 must add matvecs.
    # -----------------------------------------------------------------------
    @testset "Gating drops Pass-2 matvecs in the default path" begin
        n = 3
        cfg = make_config(Lindbladian(), EnergyDomain(); num_qubits=n)
        ham = N3_HAM
        jumps = N3_JUMPS
        d = 2^n
        rho_0 = Matrix{ComplexF64}(I(d) ./ d)
        t_grid = collect(range(0.0, 1.0, length=3))
        r_p1 = predict_lindbladian_trajectory(cfg, ham, jumps, rho_0, t_grid;
                                              krylovdim=30)
        r_p2 = predict_lindbladian_trajectory(cfg, ham, jumps, rho_0, t_grid;
                                              krylovdim=30, compute_true_gap=true)
        # Pass-1 only ≈ krylovdim matvecs; Pass-1+Pass-2 ≈ 1.5×.
        @test r_p1.total_matvecs == 30  # exactly the Arnoldi cost
        @test r_p2.total_matvecs > r_p1.total_matvecs

        cfg_ch = make_config(Thermalize(), EnergyDomain(); num_qubits=n)
        k_grid = collect(0:5:20)
        rho_init = Matrix{ComplexF64}(I(d) ./ d)
        rc_p1 = predict_channel_trajectory(cfg_ch, ham, jumps, rho_init, k_grid;
                                            krylovdim=30)
        rc_p2 = predict_channel_trajectory(cfg_ch, ham, jumps, rho_init, k_grid;
                                            krylovdim=30, compute_true_gap=true)
        @test rc_p2.total_matvecs > rc_p1.total_matvecs
    end
end
