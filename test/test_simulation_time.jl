@testset "Simulation Time (Phases 44-47)" begin

    # ---- Phase 44: QPE grid info ----
    @testset "QPE grid info" begin
        for r in [4, 8, 12, 16]
            grid = QuantumFurnace._qpe_grid_info(r, W0)
            @test grid.N == 2^r
            # Fourier relation: t0 * w0 * N = 2π
            @test grid.t0 * W0 * grid.N ≈ 2π atol=TOL_EXACT
            # Energy range endpoints
            @test grid.energy_range[1] == -grid.N÷2 * W0
            @test grid.energy_range[2] == (grid.N÷2 - 1) * W0
        end
    end

    # ---- Phase 44 (qf-e4z.18): SimulationTimeBudget struct ----
    @testset "SimulationTimeBudget struct" begin
        # New shape: per-term register triples + (with_gqsp, gqsp_degree).
        # Field order matches the struct definition in src/simulation_time.jl.
        budget = SimulationTimeBudget(
            # cost
            100.0, 5.0, 10.0, 210.0, 1000, 210000.0,
            # dissipative register
            12, 4096, 0.05, 2π/(4096*0.05), (-102.4, 102.35),
            # outer coherent (b_-)
            6, 64, 0.314, 2π/(64*0.314),
            # inner coherent (b_+)
            14, 16384, 0.628, 2π/(16384*0.628),
            # GQSP
            true, 2,
            # physics
            10.0, 0.1, 0.01, :KMS, 4, 1.5, 10.0,
            # filter
            :smooth_metropolis, Dict{Symbol,Float64}(:a => 0.333, :s => 0.4),
        )
        @test budget.oft_time == 100.0
        @test budget.b_per_be == 5.0
        @test budget.b_time == 10.0
        @test budget.per_step_time == 210.0
        @test budget.n_steps == 1000
        @test budget.total_time == 210000.0
        @test budget.r_D == 12
        @test budget.N_D == 4096
        @test budget.r_bm == 6
        @test budget.N_bm == 64
        @test budget.r_bp == 14
        @test budget.N_bp == 16384
        @test budget.with_gqsp === true
        @test budget.gqsp_degree == 2
        @test budget.construction == :KMS
        @test budget.filter_type == :smooth_metropolis
        @test budget.T == 10.0

        # Compact show
        compact = sprint(show, budget)
        @test occursin("SimulationTimeBudget", compact)
        @test occursin("r_D=12", compact)
        @test occursin("gqsp d=2", compact)

        # Verbose show
        verbose = sprint(show, MIME"text/plain"(), budget)
        @test occursin("Dissipative", verbose)
        @test occursin("Outer coh.", verbose)
        @test occursin("Inner coh.", verbose)
        @test occursin("OFT time", verbose)
        @test occursin("B per BE", verbose)
        @test occursin("B time", verbose)
        @test occursin("Per step", verbose)
        @test occursin("Total", verbose)
    end

    # ---- Phase 45: OFT time — closed-form validation ----
    @testset "OFT time — closed-form" begin
        for r in [4, 8, 10, 12]
            N = 2^r
            t0 = 2π / (N * W0)
            expected = t0 * N^2 / 4
            result = QuantumFurnace._oft_hamiltonian_time(r, W0, ones(N))
            @test result ≈ expected rtol=1e-10
        end
    end

    @testset "OFT time — weighted" begin
        r = 8
        N = 2^r
        t0 = 2π / (N * W0)
        unweighted = t0 * N^2 / 4

        # Zero weights → 0
        @test QuantumFurnace._oft_hamiltonian_time(r, W0, zeros(N)) == 0.0

        # Uniform 2.0 → 2× closed form
        @test QuantumFurnace._oft_hamiltonian_time(r, W0, fill(2.0, N)) ≈ 2.0 * unweighted atol=TOL_EXACT

        # Transition weights (≤ 1) → less than unweighted
        config = make_config(Thermalize(), TimeDomain())
        energy_labels = QuantumFurnace._create_energy_labels(r, W0)
        tw = Float64[pick_transition(config, w) for w in energy_labels]
        weighted = QuantumFurnace._oft_hamiltonian_time(r, W0, tw)
        @test 0.0 < weighted < unweighted
    end

    # ---- Phase 46: B time ----
    @testset "B time — GNS returns 0" begin
        @test QuantumFurnace._b_hamiltonian_time(nothing, nothing, 10.0, 0.1, 0.01) == 0.0
        @test QuantumFurnace._b_hamiltonian_time(Dict(), Dict(), 10.0, 0.1, 0.01) == 0.0
        @test QuantumFurnace._b_hamiltonian_time(nothing, Dict(1.0 => 0.5+0im), 10.0, 0.1, 0.01) == 0.0
    end

    @testset "B time — KMS positive" begin
        config = make_config(Thermalize(), TimeDomain())
        energy_labels = QuantumFurnace._create_energy_labels(NUM_ENERGY_BITS, W0)
        time_labels = energy_labels .* (T0 / W0)
        bm = QuantumFurnace._compute_truncated_func(
            QuantumFurnace._compute_b_minus, time_labels, BETA, SIGMA)
        bp_fn, bp_args = QuantumFurnace._select_b_plus_calculator(config)
        bp = QuantumFurnace._compute_truncated_func(bp_fn, time_labels, bp_args...)
        bt = QuantumFurnace._b_hamiltonian_time(bm, bp, BETA, SIGMA, T0)
        @test bt > 0.0
        @test isfinite(bt)
    end

    # ---- Phase 47: compute_simulation_time ----
    @testset "compute_simulation_time — KMS" begin
        config = make_config(Thermalize(), TimeDomain())
        budget = compute_simulation_time(config, TEST_HAM, 10.0)
        # Formula: per_step = 2×OFT + B
        @test budget.per_step_time ≈ 2.0 * budget.oft_time + budget.b_time
        @test budget.n_steps == ceil(Int, 10.0 / TEST_DELTA)
        @test budget.total_time ≈ budget.n_steps * budget.per_step_time
        @test budget.construction == :KMS
        @test budget.n_qubits == NUM_QUBITS
        @test budget.rescaling_factor == TEST_HAM.rescaling_factor
        @test budget.oft_time > 0.0
        @test budget.b_per_be > 0.0
        @test budget.b_time > 0.0
        @test isfinite(budget.total_time)
        @test budget.T == 10.0
        # Default config: with_gqsp=false → b_time === b_per_be (no multiplier).
        @test budget.with_gqsp === false
        @test budget.b_time ≈ budget.b_per_be rtol=TOL_EXACT
    end

    @testset "compute_simulation_time — GNS" begin
        config = make_config(Thermalize(), TimeDomain(); construction=GNS())
        budget = compute_simulation_time(config, TEST_HAM, 10.0)
        @test budget.b_per_be == 0.0
        @test budget.b_time == 0.0
        @test budget.per_step_time ≈ 2.0 * budget.oft_time
        @test budget.construction == :GNS
    end

    @testset "compute_simulation_time — TrotterDomain" begin
        config = make_config(Thermalize(), TrotterDomain())
        budget = compute_simulation_time(config, TEST_HAM, 10.0)
        @test budget.per_step_time ≈ 2.0 * budget.oft_time + budget.b_time
        @test budget.n_steps == ceil(Int, 10.0 / TEST_DELTA)
        @test budget.total_time > 0.0
        @test isfinite(budget.total_time)
    end

    @testset "compute_simulation_time — n=3 and n=4 validation" begin
        for (ham, nq) in [(TEST_HAM, NUM_QUBITS), (N3_HAM, 3)]
            config = make_config(Thermalize(), TimeDomain(); num_qubits=nq)
            budget = compute_simulation_time(config, ham, 100.0)
            @test budget.n_qubits == nq
            @test budget.total_time > 0.0
            @test isfinite(budget.total_time)
            @test budget.per_step_time > 0.0
        end
    end

    # ---- qf-e4z.18: GQSP cost-model multiplier (Form B = MW Eq. 46) ----
    # With Form B: 2d block-encoding queries per CoherentStep ⇒ b_time = 2·d·b_per_be.
    # Anticipates the Form-B circuit refactor tracked in qf-e4z.19.
    @testset "compute_simulation_time — GQSP multiplier (Form B)" begin
        # Build matched configs at d ∈ {1, 2, 3} with all other params identical.
        # `make_config` doesn't expose with_gqsp/gqsp_degree, so use Config(...) directly.
        function _gqsp_cfg(; with_gqsp::Bool, gqsp_degree::Int=1)
            Config(;
                sim = Thermalize(),
                domain = TimeDomain(),
                construction = KMS(),
                num_qubits = NUM_QUBITS,
                with_linear_combination = true,
                beta = BETA, sigma = SIGMA,
                a = BETA / 30.0, s = 0.4,
                num_energy_bits = NUM_ENERGY_BITS,
                w0 = W0, t0 = T0,
                num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
                mixing_time = 1.0, delta = TEST_DELTA,
                with_gqsp = with_gqsp, gqsp_degree = gqsp_degree,
            )
        end

        b0 = compute_simulation_time(_gqsp_cfg(with_gqsp=false),               TEST_HAM, 10.0)
        b1 = compute_simulation_time(_gqsp_cfg(with_gqsp=true,  gqsp_degree=1), TEST_HAM, 10.0)
        b2 = compute_simulation_time(_gqsp_cfg(with_gqsp=true,  gqsp_degree=2), TEST_HAM, 10.0)
        b3 = compute_simulation_time(_gqsp_cfg(with_gqsp=true,  gqsp_degree=3), TEST_HAM, 10.0)

        # `b_per_be` is the GQSP-blind per-BE cost — independent of the GQSP flag and degree
        # for matched configs (same registers, same b_± grids).
        @test b0.b_per_be ≈ b1.b_per_be rtol=TOL_EXACT
        @test b1.b_per_be ≈ b2.b_per_be rtol=TOL_EXACT
        @test b1.b_per_be ≈ b3.b_per_be rtol=TOL_EXACT

        # B-time multiplier (Form B): 2d × b_per_be.
        @test b0.b_time ≈ b0.b_per_be          rtol=TOL_EXACT  # no GQSP
        @test b1.b_time ≈ 2.0 * b1.b_per_be    rtol=TOL_EXACT  # 2·1
        @test b2.b_time ≈ 4.0 * b2.b_per_be    rtol=TOL_EXACT  # 2·2
        @test b3.b_time ≈ 6.0 * b3.b_per_be    rtol=TOL_EXACT  # 2·3

        # GQSP fields recorded in budget.
        @test b0.with_gqsp === false
        @test b1.with_gqsp === true && b1.gqsp_degree == 1
        @test b2.with_gqsp === true && b2.gqsp_degree == 2
        @test b3.with_gqsp === true && b3.gqsp_degree == 3

        # GNS path: with_gqsp is config-rejected (validate_config!), but if a budget were
        # computed for a GNS config, b_time stays 0 regardless of multiplier branch.
        cfg_gns = make_config(Thermalize(), TimeDomain(); construction=GNS())
        bg = compute_simulation_time(cfg_gns, TEST_HAM, 10.0)
        @test bg.b_per_be == 0.0
        @test bg.b_time == 0.0
    end

    # ---- qf-9z0 + qf-e4z.18: per-term registers honoured end-to-end ----
    # Vary the b_+ register independently and confirm (i) the budget records the
    # per-term values, (ii) the OFT (D-register) cost is unchanged, (iii) the
    # B-cost responds to the b_+ spacing change (factored-double-sum scaling).
    @testset "compute_simulation_time — per-term registers (qf-9z0)" begin
        function _per_term_cfg(; rbp::Int)
            # Use legacy single-register defaults for D and b_- (auto-promoted via
            # register_*_X fallback), and override only the b_+ triple.
            w0_bp = W0
            t0_bp = 2π / (2^rbp * w0_bp)
            Config(;
                sim = Thermalize(),
                domain = TimeDomain(),
                construction = KMS(),
                num_qubits = NUM_QUBITS,
                with_linear_combination = true,
                beta = BETA, sigma = SIGMA,
                a = BETA / 30.0, s = 0.4,
                num_energy_bits = NUM_ENERGY_BITS, w0 = W0, t0 = T0,
                num_energy_bits_b_plus = rbp, w0_b_plus = w0_bp, t0_b_plus = t0_bp,
                num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
                mixing_time = 1.0, delta = TEST_DELTA,
            )
        end

        ba = compute_simulation_time(_per_term_cfg(rbp=8),  TEST_HAM, 10.0)
        bb = compute_simulation_time(_per_term_cfg(rbp=10), TEST_HAM, 10.0)

        # Budget records the per-term r_b+.
        @test ba.r_bp == 8
        @test bb.r_bp == 10
        @test ba.N_bp == 256
        @test bb.N_bp == 1024

        # D-register quantities unchanged across the two configs.
        @test ba.r_D == bb.r_D == NUM_ENERGY_BITS
        @test ba.w0_D ≈ bb.w0_D rtol=TOL_EXACT
        @test ba.oft_time ≈ bb.oft_time rtol=TOL_EXACT

        # b_- register also unchanged (we only varied b_+).
        @test ba.r_bm == bb.r_bm

        # B per BE responds to the b_+ spacing change (factored sum has t0_outer · t0_inner).
        @test ba.b_per_be != bb.b_per_be
    end
end
