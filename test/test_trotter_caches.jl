"""
Tests for the qf-d0w shared-δt₀ TrottTrott scheme.

Coverage:
- Legacy single-cache TrottTrott has `nothing` per-register fields and yields
  byte-identical behaviour through B_trotter (regression).
- Shared-δt₀ constructor picks the right elementary step and computes
  per-register eigvals as vector powers of `λ_S` (no extra Trotterizations).
- Shared eigenbasis: the three caches' eigvals_t0_X are diagonal in the SAME
  `eigvecs` (no basis-alignment needed).
- B_trotter with shared caches has slope -2 in M_user vs ‖B_time‖ reference,
  confirming the dissipative-step rounding bug is removed.
- `make_trotter_for_config` picks shared-δt₀ for KMS coherent and legacy
  single-cache for GNS.
- Integer-M failure path is rejected with a clear error.
"""

using Test
using LinearAlgebra
using QuantumFurnace
using QuantumFurnace: B_trotter, B_time, _compute_b_minus, _compute_b_plus,
    _compute_b_plus_metro, _compute_truncated_func, _trotterize2

# test_helpers.jl is included once at the top of runtests.jl; the fixtures
# (BETA, SIGMA, N3_HAM, N3_TROTTER, etc.) are already in scope.

# ---------------------------------------------------------------------------
# Helpers (test-local).
# ---------------------------------------------------------------------------

function _make_b_dicts(beta::Real, sigma::Real; r_b::Int = 10,
                       T_minus::Real = 10.0, T_plus::Real = 5.0,
                       eta::Real = 1e-3)
    N = 2^r_b
    t0_minus = 2 * T_minus / N
    t0_plus  = 2 * T_plus  / N
    grid_minus = collect(-(N÷2):((N÷2)-1)) .* t0_minus
    grid_plus  = collect(-(N÷2):((N÷2)-1)) .* t0_plus
    b_minus = _compute_truncated_func(_compute_b_minus, grid_minus, beta, sigma)
    b_plus  = _compute_truncated_func(_compute_b_plus_metro, grid_plus, beta, sigma, eta, 0.25)
    return (; b_minus, b_plus, t0_minus, t0_plus)
end

function _build_jumps_in_basis(ham, basis_eigvecs, n)
    jumps = JumpOp[]
    for pauli in [[X], [Y], [Z]], site in 1:n
        op = Matrix(pad_term(pauli, n, site)) ./ sqrt(3 * n)
        op_eb = basis_eigvecs' * op * basis_eigvecs
        push!(jumps, JumpOp(op, op_eb, op == transpose(op), op == op'))
    end
    return jumps
end

@testset "qf-d0w shared-δt₀ TrottTrott" begin
    ham = N3_HAM
    n = 3
    beta = BETA  # from test_helpers.jl, locked at 10.0
    sigma = SIGMA

    @testset "Legacy constructor: per-register fields are nothing" begin
        trotter = TrottTrott(ham, 0.5, 4)
        @test trotter.eigvals_t0_b_minus === nothing
        @test trotter.eigvals_t0_b_plus  === nothing
        @test trotter.t0_b_minus === nothing
        @test trotter.t0_b_plus  === nothing
        @test trotter.t0 ≈ 0.5
        @test trotter.num_trotter_steps_per_t0 == 4
        @test length(trotter.eigvals_t0) == 8
    end

    @testset "Shared-δt₀ constructor: picks elementary step and powers eigvals" begin
        # Default thesis grid: t0_D = 10β/2^r_D, t0_bm = 20/2^r_bm, t0_bp = 10/2^r_bp.
        # With r_D = 8, r_bm = r_bp = 10, β=10:
        #   t0_D = 10·10/256 ≈ 0.391
        #   t0_b_minus_evol = β · 20/1024 ≈ 0.195
        #   t0_b_plus_evol  = β · 10/1024 ≈ 0.0977
        # binding = t0_b_plus_evol → δt₀ = 0.0977 / M_user.
        # M_D = t0_D / δt₀ = 4·M_user, M_bm = 2·M_user, M_bp = M_user.
        r_D = 8; r_bm = 10; r_bp = 10
        t0_D = 2pi / (2^r_D * pi / (5 * beta))
        t0_bm_evol = beta * 20.0 / 2^r_bm
        t0_bp_evol = beta *  10.0 / 2^r_bp

        for M_user in (1, 2, 4)
            trotter = TrottTrott(ham, t0_D, t0_bm_evol, t0_bp_evol, M_user)
            @test trotter.t0 ≈ t0_D
            @test trotter.t0_b_minus ≈ t0_bm_evol
            @test trotter.t0_b_plus  ≈ t0_bp_evol
            @test trotter.num_trotter_steps_per_t0 == 4 * M_user  # M_D
            @test length(trotter.eigvals_t0_b_minus) == 8
            @test length(trotter.eigvals_t0_b_plus)  == 8
            # eigvals_t0_X are POWERS of the elementary λ_S — verify that the
            # ratio of phases is the geometric M ratio.
            angles_D  = angle.(trotter.eigvals_t0)
            angles_bm = angle.(trotter.eigvals_t0_b_minus)
            angles_bp = angle.(trotter.eigvals_t0_b_plus)
            # angle ratios are commensurate up to mod 2π wrapping.
            @test all(abs.(exp.(im .* angles_D)  .- (exp.(im .* angles_bp)) .^ 4) .< 1e-10)
            @test all(abs.(exp.(im .* angles_bm) .- (exp.(im .* angles_bp)) .^ 2) .< 1e-10)
        end
    end

    @testset "Shared eigenbasis across registers" begin
        r_D = 8; r_bm = 10; r_bp = 10
        t0_D = 2pi / (2^r_D * pi / (5 * beta))
        t0_bm_evol = beta * 20.0 / 2^r_bm
        t0_bp_evol = beta * 10.0 / 2^r_bp
        trotter = TrottTrott(ham, t0_D, t0_bm_evol, t0_bp_evol, 2)
        # Each cache's evolution operator U_X = eigvecs · diag(eigvals_t0_X) · eigvecs'
        # should be Hermitian-conjugate of itself when combined as
        # U_X · U_X^† = I (unitary). This is a sanity check that eigvals are on the
        # unit circle and eigvecs are orthonormal.
        for evals in (trotter.eigvals_t0, trotter.eigvals_t0_b_minus, trotter.eigvals_t0_b_plus)
            @test all(abs.(abs.(evals) .- 1) .< 1e-10)
        end
        @test isapprox(trotter.eigvecs * trotter.eigvecs', I; atol = 1e-10)
    end

    @testset "B_trotter slope -2 in M_user (split caches)" begin
        # Thesis fixture: smooth Metropolis at β=10, n=3, r_D=8, r_b=10.
        r_D = 8; r_b = 10
        T_minus = 10.0; T_plus = 5.0
        eta = 1e-3
        w0_D = pi / (5 * beta)
        t0_D = 2pi / (2^r_D * w0_D)
        bd = _make_b_dicts(beta, sigma; r_b, T_minus, T_plus, eta)
        t0_bm = bd.t0_minus
        t0_bp = bd.t0_plus
        b_minus = bd.b_minus; b_plus = bd.b_plus

        # Reference: B_time at the same grid (no Trotterization).
        jumps_h = _build_jumps_in_basis(ham, ham.eigvecs, n)
        B_ref = B_time(jumps_h, ham, b_minus, b_plus, t0_bm, t0_bp, beta, sigma)
        norm_ref = opnorm(B_ref)

        errs = Float64[]
        for M_user in (1, 2, 4, 8)
            # Per-leg natural Trotter steps (qf-d0w.6): outer b_-(t/σ) → 1/σ;
            # inner b_+(τβ) → β. Coincide for σ = 1/β; written generally so the
            # formula carries to the σ ≠ 1/β regression below.
            trotter = TrottTrott(ham, t0_D, t0_bm / sigma, beta * t0_bp, M_user)
            jumps_t = _build_jumps_in_basis(ham, trotter.eigvecs, n)
            B_t = B_trotter(jumps_t, trotter, b_minus, b_plus, t0_bm, t0_bp, beta, sigma)
            # Lift back to H eigenbasis for comparison.
            U = ham.eigvecs' * trotter.eigvecs
            B_lifted = U * B_t * U'
            push!(errs, opnorm(B_ref - B_lifted))
        end
        # Errors should monotonically decrease.
        @test issorted(errs; rev = true)
        # Slope -2 in M_user: each doubling of M_user should decrease error by
        # roughly ×4. Check that ratio of consecutive errors > 3 (slack for
        # higher-order Strang corrections).
        @test all(errs[k] / errs[k+1] > 3 for k in 1:length(errs)-1)
        # Absolute scale: at M_user=8, ‖ΔB‖/‖B‖ should be < 1e-5
        # (legacy single-cache saturated at ~5e-5/2e-3 = 2.5%).
        @test errs[end] / norm_ref < 1e-5
    end

    @testset "Legacy single-cache produces M-saturated error (regression for the bug)" begin
        # With single cache at trotter.t0 = t0_D, error is dominated by τ·β / t0_D
        # rounding which does NOT decrease with M. Verify the saturation explicitly.
        r_D = 8; r_b = 10
        T_minus = 10.0; T_plus = 5.0
        eta = 1e-3
        w0_D = pi / (5 * beta)
        t0_D = 2pi / (2^r_D * w0_D)
        bd = _make_b_dicts(beta, sigma; r_b, T_minus, T_plus, eta)
        t0_bm = bd.t0_minus; t0_bp = bd.t0_plus
        b_minus = bd.b_minus; b_plus = bd.b_plus

        jumps_h = _build_jumps_in_basis(ham, ham.eigvecs, n)
        B_ref = B_time(jumps_h, ham, b_minus, b_plus, t0_bm, t0_bp, beta, sigma)

        errs_legacy = Float64[]
        for M in (4, 8, 16)
            trotter = TrottTrott(ham, t0_D, M)  # Legacy single-cache
            jumps_t = _build_jumps_in_basis(ham, trotter.eigvecs, n)
            B_t = B_trotter(jumps_t, trotter, b_minus, b_plus, t0_bm, t0_bp, beta, sigma)
            U = ham.eigvecs' * trotter.eigvecs
            B_lifted = U * B_t * U'
            push!(errs_legacy, opnorm(B_ref - B_lifted))
        end
        # Saturated: ratio of M=4 to M=16 error should be ~1 (within 50%) since
        # the t0_D-rounding error dominates.
        ratio = errs_legacy[1] / errs_legacy[end]
        @test 0.5 < ratio < 2.0
    end

    @testset "make_trotter_for_config dispatch" begin
        cfg_kms = make_config(Lindbladian(), TrotterDomain(); num_qubits=3, construction=KMS())
        trot_kms = make_trotter_for_config(N3_HAM, cfg_kms)
        @test trot_kms.eigvals_t0_b_minus !== nothing
        @test trot_kms.eigvals_t0_b_plus  !== nothing
        @test trot_kms.t0_b_minus !== nothing
        @test trot_kms.t0_b_plus  !== nothing

        cfg_gns = make_config(Lindbladian(), TrotterDomain(); num_qubits=3, construction=GNS())
        trot_gns = make_trotter_for_config(N3_HAM, cfg_gns)
        @test trot_gns.eigvals_t0_b_minus === nothing
        @test trot_gns.eigvals_t0_b_plus  === nothing
        @test trot_gns.t0_b_minus === nothing
        @test trot_gns.t0_b_plus  === nothing
    end

    @testset "Integer-M failure path is rejected" begin
        # Pick three natural steps that are NOT pairwise commensurate (rationals
        # with non-power-of-2 ratios).
        @test_throws ArgumentError TrottTrott(ham, 0.5, 0.7, 0.3, 2)
        # Negative / zero arguments are also rejected.
        @test_throws ArgumentError TrottTrott(ham, -0.1, 0.3, 0.2, 2)
        @test_throws ArgumentError TrottTrott(ham,  0.1, 0.3, 0.2, 0)
    end

    @testset "Construction: shared scheme uses ONE Strang Trotterization" begin
        # If we manually compute S_2(δt₀) and exponentiate eigvals, we should get
        # exactly the values stored in the shared trotter. This guards against
        # accidental separate-Trotterization regressions.
        r_D = 8; r_bm = 10; r_bp = 10
        t0_D = 2pi / (2^r_D * pi / (5 * beta))
        t0_bm_evol = beta * 20.0 / 2^r_bm
        t0_bp_evol = beta * 10.0 / 2^r_bp
        M_user = 2
        trotter = TrottTrott(ham, t0_D, t0_bm_evol, t0_bp_evol, M_user)

        # Reconstruct elementary δt₀ from the binding register.
        delta_t0 = min(t0_D, t0_bm_evol, t0_bp_evol) / M_user
        S = _trotterize2(ham, delta_t0, 1)
        eigvals_S, _ = eigen(S)

        # Verify that each per-register eigvals is exactly λ_S^M_X (up to the
        # eigen() ordering — we check that the multisets match).
        M_D  = trotter.num_trotter_steps_per_t0
        M_bm = round(Int, t0_bm_evol / delta_t0)
        M_bp = round(Int, t0_bp_evol / delta_t0)
        # Because eigen() may permute, compare sorted sets of phases.
        sorted_phases(v) = sort(angle.(v) .% (2π))
        @test sorted_phases(trotter.eigvals_t0)         ≈ sorted_phases(eigvals_S .^ M_D ) atol = 1e-10
        @test sorted_phases(trotter.eigvals_t0_b_minus) ≈ sorted_phases(eigvals_S .^ M_bm) atol = 1e-10
        @test sorted_phases(trotter.eigvals_t0_b_plus)  ≈ sorted_phases(eigvals_S .^ M_bp) atol = 1e-10
    end

    @testset "B_trotter slope -2 holds for σ ≠ 1/β" begin
        # Verifier note 6: the natural Trotter step for b_- is
        # `t0_grid_b_minus / σ`, NOT `β · t0_grid_b_minus`. They coincide iff
        # σ = 1/β. With σ ≠ 1/β, the formula must keep `num_steps` integer in
        # B_trotter, otherwise the dissipative-step rounding error reappears.
        # Use σ = 0.05 = 1/(2β) at β=10 (σβ ≠ 1) to stress this.
        sigma_off = 0.05  # ≠ 1/β = 0.1
        r_D = 8; r_b = 10
        T_minus = 10.0; T_plus = 5.0
        eta = 1e-3
        w0_D = pi / (5 * beta)
        t0_D = 2pi / (2^r_D * w0_D)
        bd = _make_b_dicts(beta, sigma_off; r_b, T_minus, T_plus, eta)
        t0_bm = bd.t0_minus
        t0_bp = bd.t0_plus
        b_minus = bd.b_minus; b_plus = bd.b_plus

        jumps_h = _build_jumps_in_basis(ham, ham.eigvecs, n)
        B_ref = B_time(jumps_h, ham, b_minus, b_plus, t0_bm, t0_bp, beta, sigma_off)

        errs = Float64[]
        for M_user in (1, 2, 4)
            # Note: per-leg scaling — b_-(t/σ) → 1/σ, b_+(τβ) → β.
            trotter = TrottTrott(ham, t0_D, t0_bm / sigma_off, beta * t0_bp, M_user)
            jumps_t = _build_jumps_in_basis(ham, trotter.eigvecs, n)
            B_t = B_trotter(jumps_t, trotter, b_minus, b_plus, t0_bm, t0_bp, beta, sigma_off)
            U = ham.eigvecs' * trotter.eigvecs
            B_lifted = U * B_t * U'
            push!(errs, opnorm(B_ref - B_lifted))
        end
        @test issorted(errs; rev = true)
        @test all(errs[k] / errs[k+1] > 3 for k in 1:length(errs)-1)
    end

    @testset "Run thermalize with shared trotter (smoke)" begin
        cfg = make_config(Thermalize(), TrotterDomain(); num_qubits=3, construction=KMS(),
            delta=0.01, mixing_time=0.05)
        trotter = make_trotter_for_config(N3_HAM, cfg)
        jumps_t = _build_jumps_in_basis(N3_HAM, trotter.eigvecs, n)
        result = run_thermalize(jumps_t, cfg, N3_HAM, trotter)
        @test all(isfinite, result.final_dm)
        @test all(result.trace_distances .>= 0)
        @test result.trace_distances[end] <= result.trace_distances[1] + 1e-8
    end
end
