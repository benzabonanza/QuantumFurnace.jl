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

# ============================================================================
# qf-e4z.20 — independent per-leg Trotter caches (TrotterTriple).
# ============================================================================

@testset "qf-e4z.20 TrotterTriple — independent per-leg caches" begin
    ham = N3_HAM
    n = 3
    beta = BETA
    sigma = SIGMA

    @testset "Struct invariants — rotations are unitary; composition consistent" begin
        triple = TrotterTriple(ham, 0.5, 0.3, 0.1, 4, 2, 1)
        @test triple isa AbstractTrotter
        @test triple isa TrotterTriple
        d = size(ham.data, 1)
        # Each per-leg sub-cache is a legacy single-cache TrottTrott.
        for leg in (triple.D, triple.b_minus, triple.b_plus)
            @test leg isa TrottTrott
            @test leg.eigvals_t0_b_minus === nothing
            @test leg.eigvals_t0_b_plus  === nothing
            @test isapprox(leg.eigvecs * leg.eigvecs', I; atol = 1e-10)
        end
        # Per-leg t0's recover the input.
        @test triple.D.t0       ≈ 0.5
        @test triple.b_minus.t0 ≈ 0.3
        @test triple.b_plus.t0  ≈ 0.1
        # Per-leg M counts.
        @test triple.D.num_trotter_steps_per_t0       == 4
        @test triple.b_minus.num_trotter_steps_per_t0 == 2
        @test triple.b_plus.num_trotter_steps_per_t0  == 1

        # Inter-basis rotations are unitary.
        @test isapprox(triple.R_bm_in_D  * triple.R_bm_in_D',  I; atol = 1e-10)
        @test isapprox(triple.R_bp_in_D  * triple.R_bp_in_D',  I; atol = 1e-10)
        @test isapprox(triple.R_bm_in_bp * triple.R_bm_in_bp', I; atol = 1e-10)
        # Composition: R_bm_in_D = V_bm' V_D = (V_bm' V_bp)(V_bp' V_D) = R_bm_in_bp · R_bp_in_D.
        @test isapprox(triple.R_bm_in_D, triple.R_bm_in_bp * triple.R_bp_in_D; atol = 1e-10)
    end

    @testset "All-same legs → rotations are identity to machine precision" begin
        triple = TrotterTriple(ham, 0.5, 0.5, 0.5, 4, 4, 4)
        d = size(ham.data, 1)
        @test isapprox(triple.R_bm_in_D,  I(d); atol = 1e-10)
        @test isapprox(triple.R_bp_in_D,  I(d); atol = 1e-10)
        @test isapprox(triple.R_bm_in_bp, I(d); atol = 1e-10)
    end

    @testset "Field aliasing: D-leg fields exposed via getproperty" begin
        triple = TrotterTriple(ham, 0.5, 0.3, 0.1, 4, 2, 1)
        @test triple.t0                       == triple.D.t0
        @test triple.eigvecs                  === triple.D.eigvecs
        @test triple.bohr_freqs               === triple.D.bohr_freqs
        @test triple.eigvals_t0               === triple.D.eigvals_t0
        @test triple.num_trotter_steps_per_t0 == triple.D.num_trotter_steps_per_t0
        # Legacy qf-d0w per-leg accessors map to sub-cache fields.
        @test triple.t0_b_minus           == triple.b_minus.t0
        @test triple.t0_b_plus            == triple.b_plus.t0
        @test triple.eigvals_t0_b_minus   === triple.b_minus.eigvals_t0
        @test triple.eigvals_t0_b_plus    === triple.b_plus.eigvals_t0
    end

    @testset "Builder argument validation" begin
        @test_throws ArgumentError TrotterTriple(ham, -0.1, 0.3, 0.1, 1, 1, 1)
        @test_throws ArgumentError TrotterTriple(ham,  0.5, 0.0, 0.1, 1, 1, 1)
        @test_throws ArgumentError TrotterTriple(ham,  0.5, 0.3, 0.1, 0, 1, 1)
        @test_throws ArgumentError TrotterTriple(ham,  0.5, 0.3, 0.1, 1, 0, 1)
        @test_throws ArgumentError TrotterTriple(ham,  0.5, 0.3, 0.1, 1, 1, -2)
    end

    @testset "B_trotter byte-identity at matched per-leg M counts vs legacy shared-δt₀" begin
        # Reproduce the test fixture from the qf-d0w slope test: r_D = 8,
        # r_b = 10, smooth Metropolis, β=10, n=3.
        r_D = 8; r_b = 10
        T_minus = 10.0; T_plus = 5.0
        eta = 1e-3
        w0_D = pi / (5 * beta)
        t0_D = 2pi / (2^r_D * w0_D)
        bd = _make_b_dicts(beta, sigma; r_b, T_minus, T_plus, eta)
        t0_bm = bd.t0_minus; t0_bp = bd.t0_plus
        b_minus = bd.b_minus; b_plus = bd.b_plus
        t0_bm_evol = t0_bm / sigma
        t0_bp_evol = beta * t0_bp

        # Legacy shared-δt₀ at M_user = 2 picks δt₀ = min(t0_D, t0_bm_evol,
        # t0_bp_evol) / 2. Per-leg M_X = round(t0_X / δt₀). Build the
        # equivalent TrotterTriple at exactly the same per-leg counts.
        M_user = 2
        delta_t0 = min(t0_D, t0_bm_evol, t0_bp_evol) / M_user
        M_D  = round(Int, t0_D / delta_t0)
        M_bm = round(Int, t0_bm_evol / delta_t0)
        M_bp = round(Int, t0_bp_evol / delta_t0)

        trotter_legacy = TrottTrott(ham, t0_D, t0_bm_evol, t0_bp_evol, M_user)
        triple         = TrotterTriple(ham, t0_D, t0_bm_evol, t0_bp_evol, M_D, M_bm, M_bp)

        jumps_legacy = _build_jumps_in_basis(ham, trotter_legacy.eigvecs, n)
        jumps_triple = _build_jumps_in_basis(ham, triple.eigvecs, n)
        B_legacy = B_trotter(jumps_legacy, trotter_legacy, b_minus, b_plus, t0_bm, t0_bp, beta, sigma)
        B_triple = B_trotter(jumps_triple, triple,         b_minus, b_plus, t0_bm, t0_bp, beta, sigma)

        # Lift both to Hamiltonian eigenbasis for a basis-invariant comparison.
        U_legacy = ham.eigvecs' * trotter_legacy.eigvecs
        U_triple = ham.eigvecs' * triple.eigvecs
        B_legacy_h = U_legacy * B_legacy * U_legacy'
        B_triple_h = U_triple * B_triple * U_triple'
        @test isapprox(B_legacy_h, B_triple_h; atol = 1e-12)
    end

    @testset "Slope -2 in joint M tightening (Strang 2nd order)" begin
        # When all three per-leg M counts grow together, the Trotter error
        # should drop as M^{-2}. Reproduces the qf-d0w slope test inside the
        # new independent-cache scheme.
        r_D = 8; r_b = 10
        T_minus = 10.0; T_plus = 5.0
        eta = 1e-3
        w0_D = pi / (5 * beta)
        t0_D = 2pi / (2^r_D * w0_D)
        bd = _make_b_dicts(beta, sigma; r_b, T_minus, T_plus, eta)
        t0_bm = bd.t0_minus; t0_bp = bd.t0_plus
        b_minus = bd.b_minus; b_plus = bd.b_plus
        t0_bm_evol = t0_bm / sigma
        t0_bp_evol = beta * t0_bp

        jumps_h = _build_jumps_in_basis(ham, ham.eigvecs, n)
        B_ref   = B_time(jumps_h, ham, b_minus, b_plus, t0_bm, t0_bp, beta, sigma)

        errs = Float64[]
        for M in (1, 2, 4, 8)
            triple = TrotterTriple(ham, t0_D, t0_bm_evol, t0_bp_evol, M, M, M)
            jumps_t = _build_jumps_in_basis(ham, triple.eigvecs, n)
            B_t = B_trotter(jumps_t, triple, b_minus, b_plus, t0_bm, t0_bp, beta, sigma)
            U = ham.eigvecs' * triple.eigvecs
            push!(errs, opnorm(B_ref - U * B_t * U'))
        end
        @test issorted(errs; rev = true)
        # Each doubling of M should give roughly 4× improvement (Strang -2).
        @test all(errs[k] / errs[k+1] > 3 for k in 1:length(errs)-1)
    end

    @testset "Independent leg-knob tightening (M_b_minus only)" begin
        # Hold M_D = M_b_plus = 1 (coarse) and tighten M_b_minus. The error
        # should fall as the outer-leg substep shrinks, independently of the
        # other legs (no shared-δt₀ coupling).
        r_D = 8; r_b = 10
        T_minus = 10.0; T_plus = 5.0
        eta = 1e-3
        w0_D = pi / (5 * beta)
        t0_D = 2pi / (2^r_D * w0_D)
        bd = _make_b_dicts(beta, sigma; r_b, T_minus, T_plus, eta)
        t0_bm = bd.t0_minus; t0_bp = bd.t0_plus
        b_minus = bd.b_minus; b_plus = bd.b_plus
        t0_bm_evol = t0_bm / sigma
        t0_bp_evol = beta * t0_bp

        jumps_h = _build_jumps_in_basis(ham, ham.eigvecs, n)
        B_ref   = B_time(jumps_h, ham, b_minus, b_plus, t0_bm, t0_bp, beta, sigma)

        errs = Float64[]
        for M_bm in (1, 2, 4, 8, 16)
            triple = TrotterTriple(ham, t0_D, t0_bm_evol, t0_bp_evol, 1, M_bm, 8)
            jumps_t = _build_jumps_in_basis(ham, triple.eigvecs, n)
            B_t = B_trotter(jumps_t, triple, b_minus, b_plus, t0_bm, t0_bp, beta, sigma)
            U = ham.eigvecs' * triple.eigvecs
            push!(errs, opnorm(B_ref - U * B_t * U'))
        end
        # Errors must strictly decrease as M_b_minus grows.
        @test issorted(errs; rev = true)
        # Aggressive improvement at large M_b_minus.
        @test errs[end] / errs[1] < 0.01
    end

    @testset "make_trotter_for_config returns TrotterTriple for KMS+Trotter" begin
        cfg_kms = make_config(Lindbladian(), TrotterDomain(); num_qubits=3, construction=KMS())
        trotter_kms = make_trotter_for_config(N3_HAM, cfg_kms)
        @test trotter_kms isa TrotterTriple
        @test trotter_kms isa AbstractTrotter
        # Legacy num_trotter_steps_per_t0 = 10 → all three legs at M=10.
        @test trotter_kms.D.num_trotter_steps_per_t0       == 10
        @test trotter_kms.b_minus.num_trotter_steps_per_t0 == 10
        @test trotter_kms.b_plus.num_trotter_steps_per_t0  == 10

        # GNS branch still returns legacy single-cache TrottTrott.
        cfg_gns = make_config(Lindbladian(), TrotterDomain(); num_qubits=3, construction=GNS())
        trotter_gns = make_trotter_for_config(N3_HAM, cfg_gns)
        @test trotter_gns isa TrottTrott
        @test !(trotter_gns isa TrotterTriple)
    end

    @testset "Per-leg M_user fields in Config + accessors" begin
        cfg = Config(
            sim = Lindbladian(),
            domain = TrotterDomain(),
            construction = KMS(),
            num_qubits = 3,
            with_linear_combination = true,
            beta = BETA,
            sigma = SIGMA,
            s = 0.4,
            a = BETA / 30.0,
            num_energy_bits = NUM_ENERGY_BITS,
            w0 = W0,
            t0 = T0,
            num_trotter_steps_per_t0 = 4,
            num_trotter_steps_per_t0_b_minus = 12,
        )
        # Per-leg M_b_minus overrides legacy; the other two fall back to legacy.
        @test register_M_D(cfg)        == 4
        @test register_M_b_minus(cfg)  == 12
        @test register_M_b_plus(cfg)   == 4

        trotter = make_trotter_for_config(N3_HAM, cfg)
        @test trotter isa TrotterTriple
        @test trotter.D.num_trotter_steps_per_t0       == 4
        @test trotter.b_minus.num_trotter_steps_per_t0 == 12
        @test trotter.b_plus.num_trotter_steps_per_t0  == 4
    end

    @testset "construct_lindbladian with TrotterTriple (smoke)" begin
        cfg = make_config(Lindbladian(), TrotterDomain(); num_qubits=3, construction=KMS())
        trotter = make_trotter_for_config(N3_HAM, cfg)
        @test trotter isa TrotterTriple
        jumps = _build_jumps_in_basis(N3_HAM, trotter.eigvecs, n)
        L = construct_lindbladian(jumps, cfg, N3_HAM; trotter=trotter)
        @test size(L) == (64, 64)
        # σ_β in V_D
        sigma_beta = trotter.eigvecs' * N3_HAM.eigvecs * N3_HAM.gibbs *
                     N3_HAM.eigvecs' * trotter.eigvecs
        # ‖L · σ_β‖_HS should be small (residue of detailed balance); not
        # asserting an exact value here, just that it's bounded.
        @test norm(L * vec(sigma_beta)) < 1e-2
    end
end
