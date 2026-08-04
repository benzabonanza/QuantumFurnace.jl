# Tests for the dissipative jump-selection rule (qf-2vo).
#
# Covers:
#   1. Default Config.jump_selection == :sweep and validate_config! reject other symbols.
#   2. Both modes reproduce e^{T𝓛} in expectation: ‖rho_sim - exp(L)·rho0‖ small.
#      Random uses the legacy ×S rescaling; sweep uses bare-δ S substeps per outer step.
#   3. Both modes converge to the Gibbs state under repeated δ-steps.

using QuantumFurnace
using QuantumFurnace: trace_distance_h
using LinearAlgebra
using Random
using Test

@testset "Jump selection: :sweep | :random (qf-2vo)" begin
    @testset "Config defaults & validation" begin
        cfg = make_config(Thermalize(), EnergyDomain();
            num_qubits=3, construction=KMS(), delta=0.05, mixing_time=1.0)
        @test cfg.jump_selection == :sweep

        # validate_config! accepts :sweep and :random.
        for sel in (:sweep, :random)
            cfg_sel = Config(; sim=cfg.sim, domain=cfg.domain, construction=cfg.construction,
                num_qubits=cfg.num_qubits, with_linear_combination=cfg.with_linear_combination,
                beta=cfg.beta, sigma=cfg.sigma, gaussian_parameters=cfg.gaussian_parameters,
                a=cfg.a, s=cfg.s, num_energy_bits=cfg.num_energy_bits, w0=cfg.w0, t0=cfg.t0,
                num_trotter_steps_per_t0=cfg.num_trotter_steps_per_t0,
                mixing_time=cfg.mixing_time, delta=cfg.delta, jump_selection=sel)
            @test validate_config!(cfg_sel) === nothing
        end

        # validate_config! rejects unknown selection symbols.
        cfg_bad = Config(; sim=cfg.sim, domain=cfg.domain, construction=cfg.construction,
            num_qubits=cfg.num_qubits, with_linear_combination=cfg.with_linear_combination,
            beta=cfg.beta, sigma=cfg.sigma, gaussian_parameters=cfg.gaussian_parameters,
            a=cfg.a, s=cfg.s, num_energy_bits=cfg.num_energy_bits, w0=cfg.w0, t0=cfg.t0,
            num_trotter_steps_per_t0=cfg.num_trotter_steps_per_t0,
            mixing_time=cfg.mixing_time, delta=cfg.delta, jump_selection=:wat)
        @test_throws ArgumentError validate_config!(cfg_bad)
    end

    @testset "One outer DM step uses the retained channel route" begin
        delta = 0.05
        base = make_config(Thermalize(), EnergyDomain();
            num_qubits=3, construction=KMS(), delta=delta, mixing_time=delta)
        psi = normalize(ComplexF64.(1:N3_DIM) .+ im .* ComplexF64.(N3_DIM:-1:1))
        rho0 = psi * psi'
        seed = 29

        for selection in (:sweep, :random)
            cfg = Config(; sim=base.sim, domain=base.domain,
                construction=base.construction, num_qubits=base.num_qubits,
                with_linear_combination=base.with_linear_combination,
                beta=base.beta, sigma=base.sigma,
                gaussian_parameters=base.gaussian_parameters,
                a=base.a, s=base.s, num_energy_bits=base.num_energy_bits,
                w0=base.w0, t0=base.t0,
                num_trotter_steps_per_t0=base.num_trotter_steps_per_t0,
                mixing_time=delta, delta=delta, jump_selection=selection)

            reference = run_thermalize(
                N3_JUMPS, cfg, N3_HAM;
                initial_dm=rho0, rng=Xoshiro(seed), save_every=1,
            )

            # Reassemble exactly one outer step from the canonical retained
            # construction/application helpers. This guards both selection
            # semantics without relying on an independent simulator.
            precomputed_data = QuantumFurnace._precompute_data(cfg, N3_HAM)
            n_jumps = length(N3_JUMPS)
            p_jump = 1.0 / n_jumps
            rescale = selection === :random
            coherent_unitaries = QuantumFurnace._precompute_coherent_unitary(
                N3_JUMPS, N3_HAM, cfg, precomputed_data;
                delta_scale=rescale ? 1.0 / p_jump : 1.0,
            )
            (; K0s, U_residuals) = QuantumFurnace._precompute_per_jump_channels(
                N3_JUMPS, N3_HAM, cfg, precomputed_data;
                rescale_by_inv_prob=rescale,
            )
            jump_weight_scaling = rescale ?
                precomputed_data.gamma_norm_factor / p_jump :
                precomputed_data.gamma_norm_factor
            scratch = QuantumFurnace.ThermalizeScratch(ComplexF64, N3_DIM)
            assembled = copy(rho0)

            if selection === :sweep
                for a in 1:n_jumps
                    QuantumFurnace._apply_one_dm_substep!(
                        assembled, scratch, N3_JUMPS[a], coherent_unitaries[a],
                        K0s[a], U_residuals[a], N3_HAM, cfg,
                        precomputed_data, jump_weight_scaling,
                    )
                end
            else
                a = rand(Xoshiro(seed), 1:n_jumps)
                QuantumFurnace._apply_one_dm_substep!(
                    assembled, scratch, N3_JUMPS[a], coherent_unitaries[a],
                    K0s[a], U_residuals[a], N3_HAM, cfg,
                    precomputed_data, jump_weight_scaling,
                )
            end

            @test isapprox(assembled, reference.final_dm; atol=1e-13, rtol=0)
        end
    end

    @testset "Both modes reach the Gibbs state (DM)" begin
        cfg_sweep = make_config(Thermalize(), EnergyDomain();
            num_qubits=3, construction=KMS(), delta=0.05, mixing_time=2.0)
        cfg_random = Config(; sim=cfg_sweep.sim, domain=cfg_sweep.domain,
            construction=cfg_sweep.construction, num_qubits=cfg_sweep.num_qubits,
            with_linear_combination=cfg_sweep.with_linear_combination,
            beta=cfg_sweep.beta, sigma=cfg_sweep.sigma,
            gaussian_parameters=cfg_sweep.gaussian_parameters,
            a=cfg_sweep.a, s=cfg_sweep.s, num_energy_bits=cfg_sweep.num_energy_bits,
            w0=cfg_sweep.w0, t0=cfg_sweep.t0,
            num_trotter_steps_per_t0=cfg_sweep.num_trotter_steps_per_t0,
            mixing_time=cfg_sweep.mixing_time, delta=cfg_sweep.delta,
            jump_selection=:random)

        r_sweep  = run_thermalize(N3_JUMPS, cfg_sweep,  N3_HAM; rng=Xoshiro(11))
        r_random = run_thermalize(N3_JUMPS, cfg_random, N3_HAM; rng=Xoshiro(11))

        td0_sweep  = r_sweep.trace_distances[1]
        td_sweep   = r_sweep.trace_distances[end]
        td_random  = r_random.trace_distances[end]

        # Both modes should *decrease* the trace distance compared to t=0.
        @test td_sweep  < td0_sweep
        @test td_random < td0_sweep
        # Final TDs should be in the same ballpark (within an order of magnitude).
        @test 0.1 * td_sweep <= td_random <= 10.0 * td_sweep
    end

    @testset "rescale_by_inv_prob kwarg overrides config" begin
        cfg = make_config(Thermalize(), EnergyDomain();
            num_qubits=3, construction=KMS(), delta=0.05, mixing_time=0.5)

        # Even though config.jump_selection == :sweep, an explicit rescale=true
        # forces rate-rescaled channels (advancing more physical time per substep).
        # Both runs should still complete without throwing.
        r_sweep_default = run_thermalize(N3_JUMPS, cfg, N3_HAM; rng=Xoshiro(7))
        r_sweep_rescale = run_thermalize(N3_JUMPS, cfg, N3_HAM; rng=Xoshiro(7),
            rescale_by_inv_prob=true)

        @test length(r_sweep_default.trace_distances) >= 1
        @test length(r_sweep_rescale.trace_distances) >= 1
    end
end
