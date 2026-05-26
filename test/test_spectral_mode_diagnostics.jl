# test/test_spectral_mode_diagnostics.jl
#
# Unit + lightweight-integration tests for the per-mode spectral diagnostics
# helper (`spectral_mode_diagnostics` / `SpectralModeDiagnostics`, src/diagnostics.jl),
# the cheap always-on Pass-1 diagnostics of qf-6yw.
#
# These synthetic unit tests use hand-built decompositions (no fixtures) so they
# are fast and sandbox-safe. The predictor wiring (every predict_lindbladian_
# trajectory / predict_channel_trajectory result carries a `spectral_modes`
# field) is integration-tested in `test_predict_sandbox.jl::(a)`, which asserts
# the field is present and that the steady mode (R_modes[1] ≈ σ_β, diagonal in
# the energy eigenbasis) has off_diag_weight ≈ 0.

using LinearAlgebra: I, diagm
using Test
using QuantumFurnace

@testset "spectral_mode_diagnostics (qf-6yw)" begin

    # -----------------------------------------------------------------------
    # (a) Purely diagonal R ⇒ off_diag_weight ≈ 0.
    # -----------------------------------------------------------------------
    @testset "(a) diagonal R ⇒ off_diag_weight ≈ 0" begin
        R = diagm([1.0 + 0im, 2.0 + 0im, 3.0 + 0im])
        eig = ComplexF64[0.0, -1.0, -2.0]
        c = ComplexF64[0.0, 1.0, 1.0]
        diag_obj = spectral_mode_diagnostics(eig, [R, R, R], c)
        @test all(isapprox.(diag_obj.off_diag_weight, 0.0; atol = 1e-12))
    end

    # -----------------------------------------------------------------------
    # (b) Purely off-diagonal R ⇒ off_diag_weight ≈ 1.
    # -----------------------------------------------------------------------
    @testset "(b) off-diagonal R ⇒ off_diag_weight ≈ 1" begin
        R = ComplexF64[0 2 0; 0 0 3; 1 0 0]   # zero diagonal, nonzero off-diagonal
        eig = ComplexF64[0.0, -1.0, -2.0]
        c = ComplexF64[0.0, 1.0, 1.0]
        diag_obj = spectral_mode_diagnostics(eig, [R, R, R], c)
        @test all(isapprox.(diag_obj.off_diag_weight, 1.0; atol = 1e-12))
    end

    # -----------------------------------------------------------------------
    # (c) off_diag_weight ∈ [0,1] for a random complex R.
    # -----------------------------------------------------------------------
    @testset "(c) off_diag_weight ∈ [0,1] for random complex R" begin
        R = randn(ComplexF64, 5, 5)
        eig = ComplexF64[0.0]
        c = ComplexF64[1.0]
        diag_obj = spectral_mode_diagnostics(eig, [R], c)
        w = diag_obj.off_diag_weight[1]
        @test 0.0 <= w <= 1.0
        @test isfinite(w)
    end

    # -----------------------------------------------------------------------
    # (d) Phase/scale invariance: R_k → α R_k leaves off_diag_weight fixed.
    # -----------------------------------------------------------------------
    @testset "(d) off_diag_weight phase/scale invariance" begin
        R = randn(ComplexF64, 4, 4)
        α = 3im
        eig = ComplexF64[0.0]
        c = ComplexF64[1.0]
        w_base = spectral_mode_diagnostics(eig, [R], c).off_diag_weight[1]
        w_scal = spectral_mode_diagnostics(eig, [α .* R], c).off_diag_weight[1]
        @test isapprox(w_base, w_scal; atol = 1e-10)
    end

    # -----------------------------------------------------------------------
    # (e) modal_hs_weight invariance: R_k → α R_k AND c_k → c_k/α leaves
    #     modal_hs_weight = |c_k|²‖R_k‖²_HS unchanged.
    # -----------------------------------------------------------------------
    @testset "(e) modal_hs_weight invariance under (R→αR, c→c/α)" begin
        R = randn(ComplexF64, 4, 4)
        c0 = 0.7 - 0.3im
        α = 2.0 - 1.5im
        eig = ComplexF64[0.0]
        mhw_base = spectral_mode_diagnostics(eig, [R], ComplexF64[c0]).modal_hs_weight[1]
        mhw_scal = spectral_mode_diagnostics(eig, [α .* R], ComplexF64[c0 / α]).modal_hs_weight[1]
        @test isapprox(mhw_base, mhw_scal; rtol = 1e-10)
    end

    # -----------------------------------------------------------------------
    # (f) c_abs2 == abs2(c); mode_spacing == |eig[k]-eig[k+1]|; last == Inf.
    # -----------------------------------------------------------------------
    @testset "(f) c_abs2, mode_spacing definitions" begin
        R = diagm([1.0 + 0im, 1.0 + 0im])
        eig = ComplexF64[0.0, -1.5 + 0.5im]
        c = ComplexF64[0.0, 0.4 - 1.1im]
        diag_obj = spectral_mode_diagnostics(eig, [R, R], c)
        @test diag_obj.c_abs2[1] == abs2(c[1])
        @test diag_obj.c_abs2[2] == abs2(c[2])
        @test diag_obj.mode_spacing[1] == abs(eig[1] - eig[2])
        @test diag_obj.mode_spacing[2] == Inf
    end

    # -----------------------------------------------------------------------
    # (g) Complex eigenvalues (a conjugate pair) ⇒ no error, off_diag finite
    #     and in [0,1].
    # -----------------------------------------------------------------------
    @testset "(g) complex conjugate-pair eigenvalues" begin
        R1 = randn(ComplexF64, 3, 3)
        R2 = randn(ComplexF64, 3, 3)
        R3 = randn(ComplexF64, 3, 3)
        eig = ComplexF64[0.0, -0.5 + 1.2im, -0.5 - 1.2im]  # conjugate pair
        c = ComplexF64[0.0, 1.0 + 0im, 1.0 + 0im]
        diag_obj = spectral_mode_diagnostics(eig, [R1, R2, R3], c)
        @test all(isfinite, diag_obj.off_diag_weight)
        @test all(0.0 .<= diag_obj.off_diag_weight .<= 1.0)
        @test isapprox(diag_obj.mode_spacing[2], abs(eig[2] - eig[3]); atol = 1e-12)
    end

    # -----------------------------------------------------------------------
    # (h) ArgumentError on length mismatch.
    # -----------------------------------------------------------------------
    @testset "(h) length-mismatch ArgumentError" begin
        R = diagm([1.0 + 0im, 2.0 + 0im])
        @test_throws ArgumentError spectral_mode_diagnostics(
            ComplexF64[0.0], [R, R], ComplexF64[0.0, 1.0])
        @test_throws ArgumentError spectral_mode_diagnostics(
            ComplexF64[0.0, -1.0], [R, R], ComplexF64[0.0])
    end

    # -----------------------------------------------------------------------
    # (i) operator-only spectrum (c === nothing): c-side NaN, R-side computed.
    #     The Pass-2 (krylov_spectral_gap / run_krylov_spectrum) path — no ρ₀.
    # -----------------------------------------------------------------------
    @testset "(i) c === nothing ⇒ c-side NaN, R-side computed" begin
        Rd = diagm([1.0 + 0im, 2.0 + 0im])   # diagonal ⇒ off_diag ≈ 0
        Ro = ComplexF64[0 1; 1 0]            # off-diagonal ⇒ off_diag ≈ 1
        eig = ComplexF64[0.0, -1.0]
        dobj = spectral_mode_diagnostics(eig, [Rd, Ro])   # 2-arg: no c
        @test all(isnan, dobj.c_abs2)
        @test all(isnan, dobj.modal_hs_weight)
        @test isapprox(dobj.off_diag_weight[1], 0.0; atol = 1e-12)
        @test isapprox(dobj.off_diag_weight[2], 1.0; atol = 1e-12)
        @test dobj.mode_spacing[1] == abs(eig[1] - eig[2])
        @test dobj.mode_spacing[2] == Inf
    end
end
