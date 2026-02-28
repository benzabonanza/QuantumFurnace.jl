using Test
using LinearAlgebra

# CPTP Per-Operator Completeness verification (TVAL-01)
# Verifies per-operator channel: K0_a'*K0_a + delta*R_a + U_res_a'*U_res_a = I
#
# Derivation (Chen 2023, Theorem III.1, adapted for per-operator Lie-Trotter splitting):
#   R_a scaled by 1/p_jump = n_jumps; CPTP channel uses bare delta
#   alpha = 1 - sqrt(1 - delta)
#   K0_a = I - alpha*R_a
#   S_a = (2*alpha - delta)*R_a - alpha^2*R_a^2
#   U_res_a'*U_res_a = S_a (by construction)
#   K0_a'*K0_a + delta*R_a + S_a = I (same algebraic identity, per operator)
#
# Tolerance: 1e-10 (algebraic identity; error scales as DIM^2 * eps ~ 16^2 * 1e-16 ~ 3e-13,
#   so 1e-10 gives ~300x margin for FP accumulation across DIM^2 matrix entries)

@testset "CPTP Per-Operator Completeness (TVAL-01)" begin

    @testset "EnergyDomain" begin
        config = make_config(Thermalize(),EnergyDomain(); delta=TEST_DELTA)
        ws = QuantumFurnace._build_trajectory_workspace(config, TEST_HAM, TEST_JUMPS; delta=TEST_DELTA)

        @test ws.n_jumps == length(TEST_JUMPS)
        identity = Matrix{ComplexF64}(I, DIM, DIM)
        max_err = 0.0
        for a in 1:ws.n_jumps
            completeness = ws.K0s[a]' * ws.K0s[a] + ws.delta * ws.Rs[a] + ws.U_residuals[a]' * ws.U_residuals[a]
            err = norm(completeness - identity)
            max_err = max(max_err, err)
            @test isapprox(completeness, identity; atol=1e-10)  # CPTP: K0'K0 + delta*R + U'U = I (algebraic identity)
        end
        @info "CPTP completeness (EnergyDomain)" n_jumps=ws.n_jumps max_error=max_err threshold_atol=1e-10
    end

    @testset "TimeDomain" begin
        config = make_config(Thermalize(),TimeDomain(); delta=TEST_DELTA)
        ws = QuantumFurnace._build_trajectory_workspace(config, TEST_HAM, TEST_JUMPS; delta=TEST_DELTA)

        @test ws.n_jumps == length(TEST_JUMPS)
        identity = Matrix{ComplexF64}(I, DIM, DIM)
        max_err = 0.0
        for a in 1:ws.n_jumps
            completeness = ws.K0s[a]' * ws.K0s[a] + ws.delta * ws.Rs[a] + ws.U_residuals[a]' * ws.U_residuals[a]
            err = norm(completeness - identity)
            max_err = max(max_err, err)
            @test isapprox(completeness, identity; atol=1e-10)  # CPTP: K0'K0 + delta*R + U'U = I (algebraic identity)
        end
        @info "CPTP completeness (TimeDomain)" n_jumps=ws.n_jumps max_error=max_err threshold_atol=1e-10
    end

    @testset "TrotterDomain" begin
        config = make_config(Thermalize(),TrotterDomain(); delta=TEST_DELTA)
        ws = QuantumFurnace._build_trajectory_workspace(config, TEST_HAM, TEST_TROTTER_JUMPS;
            trotter=TEST_TROTTER, delta=TEST_DELTA)

        @test ws.n_jumps == length(TEST_JUMPS)
        identity = Matrix{ComplexF64}(I, DIM, DIM)
        max_err = 0.0
        for a in 1:ws.n_jumps
            completeness = ws.K0s[a]' * ws.K0s[a] + ws.delta * ws.Rs[a] + ws.U_residuals[a]' * ws.U_residuals[a]
            err = norm(completeness - identity)
            max_err = max(max_err, err)
            @test isapprox(completeness, identity; atol=1e-10)  # CPTP: K0'K0 + delta*R + U'U = I (algebraic identity)
        end
        @info "CPTP completeness (TrotterDomain)" n_jumps=ws.n_jumps max_error=max_err threshold_atol=1e-10
    end

end
