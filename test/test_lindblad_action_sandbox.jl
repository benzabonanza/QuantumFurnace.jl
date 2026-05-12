# test/test_lindblad_action_sandbox.jl
#
# Sandbox shadow of test_lindblad_action.jl (qf-x56.2). The heavy test
# sweeps n ∈ {3, 4, 5}, β = 10 for the (q0) Energy ↔ Bohr dense-L
# 1e-9 cross-check plus several end-to-end integrator runs. This shadow
# keeps the canonical 1e-9 cross-domain invariant — the strongest
# regression check for CKG smooth-Metropolis register sizing — at the
# n = 3 cell only, plus a single integrator end-to-end run.
#
# Fixture: n = 3, β = 10, CKG smooth-Metropolis (thesis defaults
# a = 0, s = 0.25, r_D = 12).
# Dense L_super is 64 × 64 ⇒ ≪ 1 MiB. Sandbox-safe.

using LinearAlgebra: I, Hermitian, eigvals, tr, opnorm
using Test
using QuantumFurnace


@testset "Lindbladian-action integrator [sandbox shadow] (qf-x56.2)" begin

    # -----------------------------------------------------------------------
    # (q0) Energy ↔ Bohr dense Liouvillian agreement at 1e-9 (n=3 only).
    #
    # Canonical cross-domain controllability invariant per .claude/rules/
    # julia-code.md Test Suite section: at the recipe register size
    # (Eb=12, w0=0.05) the EnergyDomain Riemann sum has FP-accumulation-
    # floor error vs the closed-form BohrDomain. Threshold 1e-9 keeps the
    # NO_SANDBOX-tier invariant — loosening it would mask the kind of
    # register-sizing / index-map bug this regression is designed to catch.
    # -----------------------------------------------------------------------
    @testset "(q0) Bohr ≈ Energy dense L at 1e-9 (n=3, β=10)" begin
        beta = 10.0
        sys = make_dll_n3_system(beta)
        base_kw = (
            sim = Lindbladian(),
            construction = KMS(),
            num_qubits = 3,
            with_linear_combination = true,
            beta = beta,
            sigma = 1.0 / beta,
            a = 0.0,
            s = 0.25,
            num_energy_bits = 12,
            w0 = 0.05,
            t0 = 2π / (2^12 * 0.05),
            num_trotter_steps_per_t0 = 10,
        )
        cfg_b = Config(; domain = BohrDomain(),   base_kw...)
        cfg_e = Config(; domain = EnergyDomain(), base_kw...)
        L_b = Matrix{ComplexF64}(construct_lindbladian(sys.jumps, cfg_b, sys.ham))
        L_e = Matrix{ComplexF64}(construct_lindbladian(sys.jumps, cfg_e, sys.ham))
        rel = opnorm(L_b - L_e) / opnorm(L_b)
        @test rel < 1e-9
        @info "(q0) sandbox Bohr ≈ Energy" beta rel threshold=1e-9
    end

    # -----------------------------------------------------------------------
    # (e) Integrator end-to-end smoke: BohrDomain CKG, n=3, β=10.
    # Same horizon as the heavy testset (e) (t=120 ≈ 4·τ_mix) so the 1e-6
    # equilibrium-tail floor is genuinely resolved. Shrunk from 121 → 61
    # grid points to halve the number of Krylov-exp evaluations.
    # -----------------------------------------------------------------------
    @testset "(e) Integrator end-to-end @ n=3, β=10 (mode=:L)" begin
        beta = 10.0
        sys = make_dll_n3_system(beta)
        config = make_config(Lindbladian(), BohrDomain();
                              num_qubits = 3, construction = KMS())
        d = size(sys.ham.data, 1)
        rho_0 = Matrix{ComplexF64}(I(d) / d)
        t_grid = collect(range(0.0, 120.0, length = 61))

        res = integrate_to_gibbs(config, sys.ham, sys.jumps, rho_0, t_grid;
                                  mode = :L, krylovdim = 20, tol = 1e-10)

        @test res.all_converged
        @test res.distances[end] < 1e-6
        @test isapprox(real(tr(res.rho_final)), 1.0; atol = 1e-9)
        evs = eigvals(Hermitian((res.rho_final + res.rho_final') / 2))
        @test minimum(real.(evs)) > -1e-9

        est = estimate_mixing_time(res; model = :biexp,
                                    target_epsilon = 1e-3, extrapolate = true)
        @test isfinite(est.mixing_time)
        @test est.mixing_time > 0
        @info "(e) sandbox integrator" tau_mix=est.mixing_time dist_end=res.distances[end]
    end
end
