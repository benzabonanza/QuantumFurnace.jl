# test/test_spectral_mode_diagnostics_extra.jl
#
# Verifier-added regression tests for spectral_mode_diagnostics (qf-6yw),
# closing gaps left by test_spectral_mode_diagnostics.jl:
#   (j) the CONJUGATE-PAIR EQUALITY that justifies "no hermitisation" — the
#       existing testset (g) only checks finiteness on independent random modes.
#   (k) zero-HS-norm / m=1 edge-case defaults pinned.
#   (l) algebraic equivalence of the two off-diagonal-mass formulas.
#
# Synthetic inputs only; sandbox-safe.

using LinearAlgebra: I, diagm
using Test
using QuantumFurnace
using QuantumFurnace: spectral_mode_diagnostics, SpectralModeDiagnostics

@testset "spectral_mode_diagnostics EXTRA (qf-6yw verifier)" begin

    # -----------------------------------------------------------------------
    # (j) CONJUGATE-PAIR EQUALITY — the load-bearing claim behind "no
    #     hermitisation". L is Hermiticity-preserving, so a complex eigenpair
    #     (λ, R) has partner (conj λ, R†). off_diag_weight MUST take the SAME
    #     value on R and on R† (the dagger maps the off-diagonal index set to
    #     itself and preserves |·| and ‖·‖_HS). test (g) only checks
    #     finiteness on *independent* random modes; this pins the equality
    #     that justifies computing on the raw R_k.
    # -----------------------------------------------------------------------
    @testset "(j) off_diag_weight equal on conjugate partners R and R†" begin
        for _ in 1:50
            R = randn(ComplexF64, 5, 5)
            Rdag = Matrix(R')
            eig = ComplexF64[0.0]
            wR = spectral_mode_diagnostics(eig, [R]).off_diag_weight[1]
            wRd = spectral_mode_diagnostics(eig, [Rdag]).off_diag_weight[1]
            @test isapprox(wR, wRd; atol = 1e-13)
        end
        # A case where hermitisation genuinely changes the answer, proving the
        # code does NOT hermitise: diag-heavy asymmetric Q vs (Q+Q†)/2.
        eig = ComplexF64[0.0]
        Q = ComplexF64[2 1; 0 2]
        odw_raw = spectral_mode_diagnostics(eig, [Q]).off_diag_weight[1]
        Qh = (Q .+ Q') ./ 2
        odw_herm = spectral_mode_diagnostics(eig, [Qh]).off_diag_weight[1]
        @test !isapprox(odw_raw, odw_herm; atol = 1e-3)
    end

    # -----------------------------------------------------------------------
    # (k) zero-HS-norm R ⇒ off_diag_weight defaults to 0.0 (documented safe
    #     default); modal_hs_weight = |c|²·0 = 0. m=1 ⇒ spacing[1] = Inf.
    # -----------------------------------------------------------------------
    @testset "(k) zero-HS-norm and m=1 edge cases" begin
        Zr = zeros(ComplexF64, 3, 3)
        d0 = spectral_mode_diagnostics(ComplexF64[0.0], [Zr], ComplexF64[2.0])
        @test d0.off_diag_weight[1] == 0.0
        @test d0.modal_hs_weight[1] == 0.0
        @test d0.c_abs2[1] == 4.0
        @test d0.mode_spacing[1] == Inf

        d1 = spectral_mode_diagnostics(ComplexF64[0.0], [diagm([1.0 + 0im, 2.0 + 0im])])
        @test d1.mode_spacing[1] == Inf
        @test isnan(d1.c_abs2[1])
        @test isnan(d1.modal_hs_weight[1])
    end

    # -----------------------------------------------------------------------
    # (l) 1 − diag2/hs2 is numerically identical to Σ_{i≠j}|R[i,j]|²/‖R‖².
    # -----------------------------------------------------------------------
    @testset "(l) off_diag_weight == explicit off-diagonal mass fraction" begin
        for _ in 1:50
            R = randn(ComplexF64, 6, 6)
            hs2 = sum(abs2, R)
            offsum = sum(abs2(R[i, j]) for i in 1:6, j in 1:6 if i != j)
            direct = offsum / hs2
            code_val = spectral_mode_diagnostics(ComplexF64[0.0], [R]).off_diag_weight[1]
            @test isapprox(code_val, direct; atol = 1e-12)
        end
    end
end
