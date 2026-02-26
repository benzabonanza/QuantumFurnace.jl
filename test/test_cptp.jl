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
# Tolerance: 1e-10 (per user decision -- allows small numerical accumulation)

@testset "CPTP Per-Operator Completeness (TVAL-01)" begin

    @testset "EnergyDomain" begin
        config = make_thermalize_config(EnergyDomain(); delta=TEST_DELTA)
        precomputed = QuantumFurnace._precompute_data(config, TEST_HAM)
        scratch = QuantumFurnace.ThermalizeScratch(ComplexF64, DIM)
        fw = build_trajectoryframework(
            TEST_JUMPS, TEST_HAM, config, precomputed, scratch, TEST_DELTA
        )

        @test fw.n_jumps == length(TEST_JUMPS)
        identity = Matrix{ComplexF64}(I, DIM, DIM)
        for (a, per_op) in enumerate(fw.per_operator)
            completeness = per_op.K0' * per_op.K0 + fw.delta * per_op.R + per_op.U_residual' * per_op.U_residual
            @test isapprox(completeness, identity; atol=1e-10)
        end
    end

    @testset "TimeDomain" begin
        config = make_thermalize_config(TimeDomain(); delta=TEST_DELTA)
        precomputed = QuantumFurnace._precompute_data(config, TEST_HAM)
        scratch = QuantumFurnace.ThermalizeScratch(ComplexF64, DIM)
        fw = build_trajectoryframework(
            TEST_JUMPS, TEST_HAM, config, precomputed, scratch, TEST_DELTA
        )

        @test fw.n_jumps == length(TEST_JUMPS)
        identity = Matrix{ComplexF64}(I, DIM, DIM)
        for (a, per_op) in enumerate(fw.per_operator)
            completeness = per_op.K0' * per_op.K0 + fw.delta * per_op.R + per_op.U_residual' * per_op.U_residual
            @test isapprox(completeness, identity; atol=1e-10)
        end
    end

    @testset "TrotterDomain" begin
        config = make_thermalize_config(TrotterDomain(); delta=TEST_DELTA)
        precomputed = QuantumFurnace._precompute_data(config, TEST_TROTTER)
        scratch = QuantumFurnace.ThermalizeScratch(ComplexF64, DIM)
        fw = build_trajectoryframework(
            TEST_TROTTER_JUMPS, TEST_TROTTER, config, precomputed, scratch, TEST_DELTA
        )

        @test fw.n_jumps == length(TEST_JUMPS)
        identity = Matrix{ComplexF64}(I, DIM, DIM)
        for (a, per_op) in enumerate(fw.per_operator)
            completeness = per_op.K0' * per_op.K0 + fw.delta * per_op.R + per_op.U_residual' * per_op.U_residual
            @test isapprox(completeness, identity; atol=1e-10)
        end
    end

end
