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

    @testset "BohrDomain (via _precompute_per_jump_channels)" begin
        config = make_config(Thermalize(), BohrDomain(); delta=TEST_DELTA)
        precomputed_data = QuantumFurnace._precompute_data(config, TEST_HAM)
        (; K0s, U_residuals) = QuantumFurnace._precompute_per_jump_channels(
            TEST_JUMPS, TEST_HAM, config, precomputed_data;
            rescale_by_inv_prob=true,
        )

        n_jumps = length(TEST_JUMPS)
        identity = Matrix{ComplexF64}(I, DIM, DIM)
        # Need Rs to verify CPTP completeness -- recompute them
        builder_scratch = QuantumFurnace.ThermalizeScratch(ComplexF64, DIM)
        p_jump = 1.0 / n_jumps
        max_err = 0.0
        for a in 1:n_jumps
            QuantumFurnace._precompute_R([TEST_JUMPS[a]], TEST_HAM, config, precomputed_data, builder_scratch)
            R_a = copy(builder_scratch.R)
            R_a .*= (1.0 / p_jump)
            completeness = K0s[a]' * K0s[a] + TEST_DELTA * R_a + U_residuals[a]' * U_residuals[a]
            err = norm(completeness - identity)
            max_err = max(max_err, err)
            @test isapprox(completeness, identity; atol=1e-10)
        end
        @info "CPTP completeness (BohrDomain)" n_jumps max_error=max_err threshold_atol=1e-10
    end

    @testset "DM precomputation matches trajectory workspace" begin
        for (domain, label, extra_kw) in [
            (EnergyDomain(), "Energy", (;)),
            (TimeDomain(), "Time", (;)),
            (TrotterDomain(), "Trotter", (; trotter=TEST_TROTTER)),
        ]
            config = make_config(Thermalize(), domain; delta=TEST_DELTA)
            jumps = domain isa TrotterDomain ? TEST_TROTTER_JUMPS : TEST_JUMPS
            ham = TEST_HAM

            # DM precomputation path
            ham_or_trott = domain isa TrotterDomain ? TEST_TROTTER : ham
            precomputed_data = QuantumFurnace._precompute_data(config, ham_or_trott)
            (; K0s, U_residuals) = QuantumFurnace._precompute_per_jump_channels(
                jumps, ham_or_trott, config, precomputed_data; rescale_by_inv_prob=true,
            )

            # Trajectory workspace path (reference). Pin rescale_by_inv_prob=true to match
            # the DM precomp path's `rescale_by_inv_prob=true` above; the default now follows
            # `config.jump_selection` (= :sweep ⇒ bare-rate channels) per qf-2vo.
            ws = QuantumFurnace._build_trajectory_workspace(config, ham, jumps; extra_kw...,
                delta=TEST_DELTA, rescale_by_inv_prob=true)

            for a in 1:length(jumps)
                @test isapprox(K0s[a], ws.K0s[a]; atol=1e-15)
                @test isapprox(U_residuals[a], ws.U_residuals[a]; atol=1e-15)
            end
            @info "DM precomp matches trajectory ($label)" n_jumps=length(jumps)
        end
    end

end
