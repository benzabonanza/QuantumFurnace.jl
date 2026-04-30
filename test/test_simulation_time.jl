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

    # ---- Phase 44: SimulationTimeBudget struct ----
    @testset "SimulationTimeBudget struct" begin
        budget = SimulationTimeBudget(
            100.0, 10.0, 210.0, 1000, 210000.0,
            12, 4096, 0.05, 2π/(4096*0.05), (-102.4, 102.35),
            10.0, 0.1, 0.01, :KMS, 4, 1.5, 10.0,
            :smooth_metropolis, Dict{Symbol,Float64}(:a => 0.333, :s => 0.4),
        )
        @test budget.oft_time == 100.0
        @test budget.b_time == 10.0
        @test budget.per_step_time == 210.0
        @test budget.n_steps == 1000
        @test budget.total_time == 210000.0
        @test budget.construction == :KMS
        @test budget.filter_type == :smooth_metropolis
        @test budget.T == 10.0

        # Compact show
        compact = sprint(show, budget)
        @test occursin("SimulationTimeBudget", compact)
        @test occursin("r=12", compact)

        # Verbose show
        verbose = sprint(show, MIME"text/plain"(), budget)
        @test occursin("OFT time", verbose)
        @test occursin("Per step", verbose)
        @test occursin("2×", verbose)
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
        @test budget.b_time > 0.0
        @test isfinite(budget.total_time)
        @test budget.T == 10.0
    end

    @testset "compute_simulation_time — GNS" begin
        config = make_config(Thermalize(), TimeDomain(); construction=GNS())
        budget = compute_simulation_time(config, TEST_HAM, 10.0)
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
end
