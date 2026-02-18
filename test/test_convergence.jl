"""
Convergence tracking tests (Phase 16, Plan 02).

Tests cover:
- ConvergenceData struct construction and field access
- Observable builders (eigenbasis and Trotter basis)
- Gibbs state helpers (_gibbs_in_trotter_basis, _compute_gibbs_observable_values)
- ConvergenceData Dict and BSON serialization round-trips
- Integration test: run_trajectories_convergence with trace distance convergence (CONV-01),
  ZZ correlation tracking (CONV-02), and energy convergence (CONV-03)
- Determinism: same seed produces identical convergence data
- Programmatic access (CONV-04 success criterion 4)
"""

using Random
using BSON

@testset "Convergence tracking" begin

    # -----------------------------------------------------------------------
    # Testset 1: ConvergenceData struct
    # -----------------------------------------------------------------------
    @testset "ConvergenceData struct" begin
        n_obs = 3
        n_checkpoints = 4
        conv = ConvergenceData(
            [100, 100, 100, 100],
            [100, 200, 300, 400],
            [0.5, 0.4, 0.3, 0.2],
            ["ZZ_12", "ZZ_23", "H"],
            rand(n_obs, n_checkpoints),
            [0.1, -0.2, -1.5],
        )

        @test conv.batch_sizes == [100, 100, 100, 100]
        @test conv.cumulative_n_traj == [100, 200, 300, 400]
        @test conv.trace_distances == [0.5, 0.4, 0.3, 0.2]
        @test conv.observable_names == ["ZZ_12", "ZZ_23", "H"]
        @test conv.observable_gibbs_values == [0.1, -0.2, -1.5]

        @test conv.batch_sizes isa Vector{Int}
        @test conv.cumulative_n_traj isa Vector{Int}
        @test conv.trace_distances isa Vector{Float64}
        @test conv.observable_names isa Vector{String}
        @test conv.observable_values isa Matrix{Float64}
        @test conv.observable_gibbs_values isa Vector{Float64}

        @test size(conv.observable_values) == (n_obs, n_checkpoints)
    end

    # -----------------------------------------------------------------------
    # Testset 4: _gibbs_in_trotter_basis
    # -----------------------------------------------------------------------
    @testset "_gibbs_in_trotter_basis" begin
        gibbs_trotter = QuantumFurnace._gibbs_in_trotter_basis(TEST_HAM, TEST_TROTTER)

        # Hermitian
        @test isapprox(Matrix(gibbs_trotter), Matrix(gibbs_trotter)'; atol=1e-12)

        # Trace 1 (still a valid density matrix after basis change)
        @test isapprox(real(tr(gibbs_trotter)), 1.0; atol=1e-12)

        # All eigenvalues non-negative (valid density matrix)
        eigs = eigvals(Matrix(gibbs_trotter))
        @test all(real.(eigs) .>= -1e-12)
    end

    # -----------------------------------------------------------------------
    # Testset 5: _compute_gibbs_observable_values
    # -----------------------------------------------------------------------
    @testset "_compute_gibbs_observable_values" begin
        observables, _ = build_preset_trajectory_observables(TEST_HAM, NUM_QUBITS)
        obs_gibbs = QuantumFurnace._compute_gibbs_observable_values(TEST_GIBBS, observables)

        @test length(obs_gibbs) == length(observables)
        @test all(isfinite, obs_gibbs)
        @test all(x -> isa(x, Real), obs_gibbs)

        # Energy value (first entry, H) matches analytical Gibbs energy
        boltz = exp.(-BETA .* TEST_HAM.eigvals)
        gibbs_energy_analytical = sum(TEST_HAM.eigvals .* boltz) / sum(boltz)
        @test isapprox(obs_gibbs[1], gibbs_energy_analytical; atol=1e-10)
    end

    # -----------------------------------------------------------------------
    # Testset 6: ConvergenceData Dict serialization round-trip
    # -----------------------------------------------------------------------
    @testset "ConvergenceData Dict serialization round-trip" begin
        obs_vals = rand(2, 3)
        conv = ConvergenceData(
            [100, 100, 100],
            [100, 200, 300],
            [0.4, 0.25, 0.15],
            ["ZZ_12", "H"],
            obs_vals,
            [0.1, -1.5],
        )

        d = QuantumFurnace._convergence_to_dict(conv)

        # Verify Dict has all expected keys
        @test d isa Dict
        @test haskey(d, :batch_sizes)
        @test haskey(d, :cumulative_n_traj)
        @test haskey(d, :trace_distances)
        @test haskey(d, :observable_names)
        @test haskey(d, :observable_values)
        @test haskey(d, :observable_gibbs_values)

        # Reconstruct
        conv2 = QuantumFurnace._dict_to_convergence(d)

        @test conv2.batch_sizes == conv.batch_sizes
        @test conv2.cumulative_n_traj == conv.cumulative_n_traj
        @test conv2.trace_distances == conv.trace_distances
        @test conv2.observable_names == conv.observable_names
        @test isapprox(conv2.observable_values, conv.observable_values)
        @test conv2.observable_gibbs_values == conv.observable_gibbs_values

        # Forward compatibility: remove observable_gibbs_values, should default to Float64[]
        d_compat = copy(d)
        delete!(d_compat, :observable_gibbs_values)
        conv3 = QuantumFurnace._dict_to_convergence(d_compat)
        @test conv3.observable_gibbs_values == Float64[]
    end

    # -----------------------------------------------------------------------
    # Testset 7: ConvergenceData BSON round-trip
    # -----------------------------------------------------------------------
    @testset "ConvergenceData BSON round-trip" begin
        mktempdir() do tmpdir
            obs_vals = rand(3, 4)
            conv = ConvergenceData(
                [50, 50, 50, 50],
                [50, 100, 150, 200],
                [0.6, 0.45, 0.3, 0.2],
                ["ZZ_12", "ZZ_23", "H"],
                obs_vals,
                [0.05, -0.1, -2.0],
            )

            d = QuantumFurnace._convergence_to_dict(conv)
            bson_path = joinpath(tmpdir, "conv_test.bson")
            BSON.bson(bson_path, d)

            loaded_d = BSON.load(bson_path)
            conv2 = QuantumFurnace._dict_to_convergence(loaded_d)

            @test conv2.batch_sizes == conv.batch_sizes
            @test conv2.cumulative_n_traj == conv.cumulative_n_traj
            @test isapprox(conv2.trace_distances, conv.trace_distances)
            @test conv2.observable_names == conv.observable_names
            @test isapprox(conv2.observable_values, conv.observable_values)
            @test isapprox(conv2.observable_gibbs_values, conv.observable_gibbs_values)
        end
    end

    # -----------------------------------------------------------------------
    # Testset 8: run_trajectories_convergence integration (CONV-01, CONV-02, CONV-03)
    # -----------------------------------------------------------------------
    @testset "run_trajectories_convergence integration" begin
        # mixing_time=60.0 accounts for the corrected (slower) trajectory evolution rate
        # (the per-step CPTP channel uses bare delta, not delta*n_jumps).
        # Previous mixing_time=5.0 was sufficient when delta_eff=delta*n_jumps=0.12.
        config = make_thermalize_config(EnergyDomain(); with_coherent=true, delta=0.01, mixing_time=60.0)
        observables, names = build_preset_trajectory_observables(TEST_HAM, NUM_QUBITS)
        obs_gibbs = QuantumFurnace._compute_gibbs_observable_values(TEST_GIBBS, observables)

        # Ground state in eigenbasis: first basis vector
        psi0 = zeros(ComplexF64, DIM)
        psi0[1] = 1.0

        traj_result = run_trajectories_convergence(
            TEST_JUMPS, config, psi0, TEST_HAM;
            gibbs=TEST_GIBBS,
            observables=observables,
            observable_names=names,
            batch_size=200,
            n_batches=5,
            seed=42,
        )
        conv_data = traj_result.convergence

        # Basic structure
        @test conv_data isa ConvergenceData
        @test length(conv_data.batch_sizes) == 5
        @test all(conv_data.batch_sizes .== 200)
        @test conv_data.cumulative_n_traj == [200, 400, 600, 800, 1000]
        @test length(conv_data.trace_distances) == 5
        @test all(isfinite, conv_data.trace_distances)
        @test all(conv_data.trace_distances .>= 0)

        # CONV-01: Trace distance decreases (last batch < first batch)
        @test conv_data.trace_distances[end] < conv_data.trace_distances[1]

        # Observable values structure
        @test size(conv_data.observable_values) == (7, 5)
        @test conv_data.observable_names == names
        @test length(conv_data.observable_gibbs_values) == 7

        # CONV-02/CONV-03: Energy observable converges toward Gibbs value
        energy_idx = findfirst(==("H"), names)
        @test abs(conv_data.observable_values[energy_idx, end] - obs_gibbs[energy_idx]) <
              abs(conv_data.observable_values[energy_idx, 1] - obs_gibbs[energy_idx])

        # TrajectoryResult checks
        @test traj_result isa TrajectoryResult
        @test traj_result.n_trajectories == 1000
        @test traj_result.seed == 42
        @test size(traj_result.rho_mean) == (DIM, DIM)

        # rho_mean is approximately a valid density matrix
        @test isapprox(real(tr(traj_result.rho_mean)), 1.0; atol=1e-6)
        @test isapprox(traj_result.rho_mean, traj_result.rho_mean'; atol=1e-10)
    end

    # -----------------------------------------------------------------------
    # Testset 9: run_trajectories_convergence determinism
    # -----------------------------------------------------------------------
    @testset "run_trajectories_convergence determinism" begin
        config = make_thermalize_config(EnergyDomain(); with_coherent=true, delta=0.01, mixing_time=5.0)
        observables, names = build_preset_trajectory_observables(TEST_HAM, NUM_QUBITS)

        psi0 = zeros(ComplexF64, DIM)
        psi0[1] = 1.0

        conv1 = run_trajectories_convergence(
            TEST_JUMPS, config, psi0, TEST_HAM;
            gibbs=TEST_GIBBS,
            observables=observables,
            observable_names=names,
            batch_size=50,
            n_batches=2,
            seed=12345,
        ).convergence

        conv2 = run_trajectories_convergence(
            TEST_JUMPS, config, psi0, TEST_HAM;
            gibbs=TEST_GIBBS,
            observables=observables,
            observable_names=names,
            batch_size=50,
            n_batches=2,
            seed=12345,
        ).convergence

        # Bitwise identical
        @test conv1.trace_distances == conv2.trace_distances
        @test conv1.observable_values == conv2.observable_values
    end

    # -----------------------------------------------------------------------
    # Testset 10: Convergence data accessible programmatically (CONV-04)
    # -----------------------------------------------------------------------
    @testset "Convergence data accessible programmatically" begin
        config = make_thermalize_config(EnergyDomain(); with_coherent=true, delta=0.01, mixing_time=5.0)
        observables, names = build_preset_trajectory_observables(TEST_HAM, NUM_QUBITS)

        psi0 = zeros(ComplexF64, DIM)
        psi0[1] = 1.0

        conv_data = run_trajectories_convergence(
            TEST_JUMPS, config, psi0, TEST_HAM;
            gibbs=TEST_GIBBS,
            observables=observables,
            observable_names=names,
            batch_size=50,
            n_batches=2,
            seed=99,
        ).convergence

        # trace_distances is a Vector{Float64} that can be indexed and iterated
        @test conv_data.trace_distances isa Vector{Float64}
        @test length(conv_data.trace_distances) == 2
        td_first = conv_data.trace_distances[1]
        @test td_first isa Float64

        # observable_values is a Matrix{Float64} sliceable by row (per-observable) and column (per-batch)
        @test conv_data.observable_values isa Matrix{Float64}
        row_slice = conv_data.observable_values[1, :]  # first observable across batches
        @test row_slice isa Vector{Float64}
        @test length(row_slice) == 2
        col_slice = conv_data.observable_values[:, 1]  # all observables at first batch
        @test col_slice isa Vector{Float64}
        @test length(col_slice) == 7

        # observable_names can be used to look up specific observables
        idx = findfirst(==("H"), conv_data.observable_names)
        @test idx !== nothing
        h_curve = conv_data.observable_values[idx, :]
        @test h_curve isa Vector{Float64}
        @test length(h_curve) == 2
    end

    # -----------------------------------------------------------------------
    # Testset 11: ConvergenceData backward-compatible constructor
    # -----------------------------------------------------------------------
    @testset "ConvergenceData backward-compatible constructor" begin
        # 6-argument constructor (Phase 16 pattern)
        obs_vals = rand(2, 3)
        conv6 = ConvergenceData(
            [100, 100, 100],
            [100, 200, 300],
            [0.4, 0.3, 0.2],
            ["ZZ_12", "H"],
            obs_vals,
            [0.1, -1.5],
        )

        # Defaults for Phase 17 fields
        @test conv6.converged == false
        @test isnan(conv6.final_relative_change)
        @test conv6.consecutive_stable_batches == 0
        @test conv6.total_batches == length(conv6.batch_sizes)

        # 10-argument constructor (Phase 17 full constructor)
        conv10 = ConvergenceData(
            [200, 200],
            [200, 400],
            [0.3, 0.15],
            ["ZZ_12", "H"],
            rand(2, 2),
            [0.05, -2.0],
            true,
            0.02,
            3,
            2,
        )

        @test conv10.converged == true
        @test conv10.final_relative_change == 0.02
        @test conv10.consecutive_stable_batches == 3
        @test conv10.total_batches == 2
    end

    # -----------------------------------------------------------------------
    # Testset 12: _windowed_relative_change unit tests
    # -----------------------------------------------------------------------
    @testset "_windowed_relative_change unit tests" begin
        # Insufficient data: 2 points, window_size=3 needs 6
        @test QuantumFurnace._windowed_relative_change([0.5, 0.4], 3) == Inf

        # Exactly 2*window_size data points returns a finite positive value
        result = QuantumFurnace._windowed_relative_change([0.5, 0.45, 0.4, 0.38, 0.37, 0.36], 3)
        @test isfinite(result) && result > 0

        # Known values: [1.0, 1.0, 1.0, 0.5, 0.5, 0.5] with window_size=3
        # mean_previous = mean([1.0, 1.0, 1.0]) = 1.0
        # mean_recent   = mean([0.5, 0.5, 0.5]) = 0.5
        # relative_change = abs(0.5 - 1.0) / max(abs(1.0), eps()) = 0.5
        @test isapprox(QuantumFurnace._windowed_relative_change([1.0, 1.0, 1.0, 0.5, 0.5, 0.5], 3), 0.5; atol=1e-12)

        # Converged values: all identical -> relative change ~0.0
        @test QuantumFurnace._windowed_relative_change([0.1, 0.1, 0.1, 0.1, 0.1, 0.1], 3) < 1e-10

        # window_size=1: two-element data [0.4, 0.3]
        # mean_previous = 0.4, mean_recent = 0.3
        # relative_change = abs(0.3 - 0.4) / max(0.4, eps()) = 0.1/0.4 = 0.25
        @test isapprox(QuantumFurnace._windowed_relative_change([0.4, 0.3], 1), 0.25; atol=1e-12)

        # Edge: window_size=2 with exactly 4 elements [2.0, 2.0, 1.0, 1.0]
        # mean_previous = 2.0, mean_recent = 1.0, relative_change = 1.0/2.0 = 0.5
        @test isapprox(QuantumFurnace._windowed_relative_change([2.0, 2.0, 1.0, 1.0], 2), 0.5; atol=1e-12)
    end

    # -----------------------------------------------------------------------
    # Testset 13: run_trajectories_adaptive convergence (CONV-04)
    # -----------------------------------------------------------------------
    @testset "run_trajectories_adaptive convergence (CONV-04)" begin
        config = make_thermalize_config(EnergyDomain(); with_coherent=true, delta=0.01, mixing_time=5.0)
        observables, names = build_preset_trajectory_observables(TEST_HAM, NUM_QUBITS)

        # Ground state in eigenbasis
        psi0 = zeros(ComplexF64, DIM)
        psi0[1] = 1.0

        traj_result = run_trajectories_adaptive(
            TEST_JUMPS, config, psi0, TEST_HAM;
            gibbs=TEST_GIBBS,
            observables=observables,
            observable_names=names,
            batch_size=200,
            n_max=20_000,
            convergence_threshold=0.05,
            patience=3,
            min_batches=5,
            window_size=3,
            seed=42,
        )
        conv_data = traj_result.convergence

        # System should converge with generous threshold
        @test conv_data.converged == true

        # Converged before hitting hard cap: cld(20000, 200) = 100
        @test conv_data.total_batches < 100

        # Patience was met
        @test conv_data.consecutive_stable_batches >= 3

        # Final relative change is below threshold
        @test conv_data.final_relative_change < 0.05
        @test isfinite(conv_data.final_relative_change)

        # Trace distances are all positive and finite
        @test all(isfinite, conv_data.trace_distances)
        @test all(x -> x >= 0, conv_data.trace_distances)

        # Observable values matrix shape matches
        @test size(conv_data.observable_values) == (7, conv_data.total_batches)

        # TrajectoryResult has valid density matrix
        @test isapprox(real(tr(traj_result.rho_mean)), 1.0; atol=1e-6)
        @test isapprox(traj_result.rho_mean, traj_result.rho_mean'; atol=1e-10)

        # Total trajectories match batch accounting
        @test traj_result.n_trajectories == conv_data.total_batches * 200
    end

    # -----------------------------------------------------------------------
    # Testset 14: run_trajectories_adaptive hard cap (CONV-05)
    # -----------------------------------------------------------------------
    @testset "run_trajectories_adaptive hard cap (CONV-05)" begin
        config = make_thermalize_config(EnergyDomain(); with_coherent=true, delta=0.01, mixing_time=5.0)
        observables, names = build_preset_trajectory_observables(TEST_HAM, NUM_QUBITS)

        psi0 = zeros(ComplexF64, DIM)
        psi0[1] = 1.0

        # n_max=500, batch_size=50 -> max_batches=10, threshold=0.0001 is nearly impossible
        traj_result = run_trajectories_adaptive(
            TEST_JUMPS, config, psi0, TEST_HAM;
            gibbs=TEST_GIBBS,
            observables=observables,
            observable_names=names,
            batch_size=50,
            n_max=500,
            convergence_threshold=0.0001,
            patience=3,
            min_batches=5,
            window_size=3,
            seed=99,
        )
        conv_data = traj_result.convergence

        # Should NOT converge with 0.01% threshold and only 500 trajectories
        @test conv_data.converged == false

        # Ran all batches to the cap: cld(500, 50) = 10
        @test conv_data.total_batches == cld(500, 50)

        # Total trajectories match
        @test traj_result.n_trajectories == 500

        # Patience NOT met
        @test conv_data.consecutive_stable_batches < 3

        # All trace distances positive and finite
        @test all(isfinite, conv_data.trace_distances)
        @test all(x -> x >= 0, conv_data.trace_distances)
    end

    # -----------------------------------------------------------------------
    # Testset 15: run_trajectories_adaptive determinism
    # -----------------------------------------------------------------------
    @testset "run_trajectories_adaptive determinism" begin
        config = make_thermalize_config(EnergyDomain(); with_coherent=true, delta=0.01, mixing_time=5.0)
        observables, names = build_preset_trajectory_observables(TEST_HAM, NUM_QUBITS)

        psi0 = zeros(ComplexF64, DIM)
        psi0[1] = 1.0

        conv1 = run_trajectories_adaptive(
            TEST_JUMPS, config, psi0, TEST_HAM;
            gibbs=TEST_GIBBS,
            observables=observables,
            observable_names=names,
            batch_size=100,
            n_max=5000,
            convergence_threshold=0.05,
            patience=3,
            seed=12345,
        ).convergence

        conv2 = run_trajectories_adaptive(
            TEST_JUMPS, config, psi0, TEST_HAM;
            gibbs=TEST_GIBBS,
            observables=observables,
            observable_names=names,
            batch_size=100,
            n_max=5000,
            convergence_threshold=0.05,
            patience=3,
            seed=12345,
        ).convergence

        # Bitwise identical
        @test conv1.trace_distances == conv2.trace_distances
        @test conv1.observable_values == conv2.observable_values
        @test conv1.converged == conv2.converged
        @test conv1.total_batches == conv2.total_batches
    end

    # -----------------------------------------------------------------------
    # Testset 16: ConvergenceData serialization with adaptive fields
    # -----------------------------------------------------------------------
    @testset "ConvergenceData serialization with adaptive fields" begin
        # Create ConvergenceData with all 10 fields populated
        conv = ConvergenceData(
            [200, 200, 200],
            [200, 400, 600],
            [0.35, 0.20, 0.12],
            ["ZZ_12", "ZZ_23", "H"],
            rand(3, 3),
            [0.05, -0.1, -2.0],
            true,     # converged
            0.03,     # final_relative_change
            3,        # consecutive_stable_batches
            3,        # total_batches
        )

        d = QuantumFurnace._convergence_to_dict(conv)

        # Verify Dict has all 10 keys including Phase 17 fields
        @test haskey(d, :converged)
        @test haskey(d, :final_relative_change)
        @test haskey(d, :consecutive_stable_batches)
        @test haskey(d, :total_batches)

        # Round-trip through Dict
        conv2 = QuantumFurnace._dict_to_convergence(d)

        @test conv2.batch_sizes == conv.batch_sizes
        @test conv2.cumulative_n_traj == conv.cumulative_n_traj
        @test conv2.trace_distances == conv.trace_distances
        @test conv2.observable_names == conv.observable_names
        @test isapprox(conv2.observable_values, conv.observable_values)
        @test conv2.observable_gibbs_values == conv.observable_gibbs_values
        @test conv2.converged == true
        @test conv2.final_relative_change == 0.03
        @test conv2.consecutive_stable_batches == 3
        @test conv2.total_batches == 3

        # Forward compatibility: Dict WITHOUT the 4 new keys (simulating old Phase 16 data)
        d_old = Dict{Symbol, Any}(
            :batch_sizes => [100, 100],
            :cumulative_n_traj => [100, 200],
            :trace_distances => [0.5, 0.3],
            :observable_names => ["ZZ_12"],
            :observable_values => rand(1, 2),
            :observable_gibbs_values => [0.1],
        )
        conv_old = QuantumFurnace._dict_to_convergence(d_old)
        @test conv_old.converged == false
        @test isnan(conv_old.final_relative_change)
        @test conv_old.consecutive_stable_batches == 0
        @test conv_old.total_batches == 2
    end

    # -----------------------------------------------------------------------
    # Testset 17: ConvergenceData BSON round-trip with adaptive fields
    # -----------------------------------------------------------------------
    @testset "ConvergenceData BSON round-trip with adaptive fields" begin
        mktempdir() do tmpdir
            conv = ConvergenceData(
                [150, 150, 150, 150],
                [150, 300, 450, 600],
                [0.4, 0.25, 0.15, 0.10],
                ["ZZ_12", "ZZ_23", "H"],
                rand(3, 4),
                [0.05, -0.1, -2.0],
                true,     # converged
                0.018,    # final_relative_change
                3,        # consecutive_stable_batches
                4,        # total_batches
            )

            d = QuantumFurnace._convergence_to_dict(conv)
            bson_path = joinpath(tmpdir, "conv_adaptive_test.bson")
            BSON.bson(bson_path, d)

            loaded_d = BSON.load(bson_path)
            conv2 = QuantumFurnace._dict_to_convergence(loaded_d)

            @test conv2.batch_sizes == conv.batch_sizes
            @test conv2.cumulative_n_traj == conv.cumulative_n_traj
            @test isapprox(conv2.trace_distances, conv.trace_distances)
            @test conv2.observable_names == conv.observable_names
            @test isapprox(conv2.observable_values, conv.observable_values)
            @test isapprox(conv2.observable_gibbs_values, conv.observable_gibbs_values)
            @test conv2.converged == true
            @test isapprox(conv2.final_relative_change, 0.018; atol=1e-12)
            @test conv2.consecutive_stable_batches == 3
            @test conv2.total_batches == 4
        end
    end

    # -----------------------------------------------------------------------
    # Testset 18: Adaptive result programmatic access (CONV-04 extended)
    # -----------------------------------------------------------------------
    @testset "Adaptive result programmatic access (CONV-04 extended)" begin
        config = make_thermalize_config(EnergyDomain(); with_coherent=true, delta=0.01, mixing_time=5.0)
        observables, names = build_preset_trajectory_observables(TEST_HAM, NUM_QUBITS)

        psi0 = zeros(ComplexF64, DIM)
        psi0[1] = 1.0

        conv_data = run_trajectories_adaptive(
            TEST_JUMPS, config, psi0, TEST_HAM;
            gibbs=TEST_GIBBS,
            observables=observables,
            observable_names=names,
            batch_size=50,
            n_max=1000,
            convergence_threshold=0.1,
            seed=77,
        ).convergence

        # New diagnostic fields are accessible and correctly typed
        @test conv_data.converged isa Bool
        @test conv_data.final_relative_change isa Float64
        @test conv_data.consecutive_stable_batches isa Int
        @test conv_data.total_batches isa Int
        @test conv_data.total_batches > 0

        # Cumulative trajectory count matches batch accounting
        @test conv_data.cumulative_n_traj[end] == conv_data.total_batches * 50
    end

    # -----------------------------------------------------------------------
    # Testset 21: build_preset_trajectory_observables (eigenbasis)
    # -----------------------------------------------------------------------
    @testset "build_preset_trajectory_observables (eigenbasis)" begin
        observables, names = build_preset_trajectory_observables(TEST_HAM, NUM_QUBITS)

        # Returns exactly 7 observables and 7 names
        @test length(observables) == 7
        @test length(names) == 7
        @test names == ["H", "Mz", "XX_avg", "YY_avg", "ZZ_avg", "Mz_stagg", "Z1"]

        # All matrices are DIM x DIM ComplexF64
        for obs in observables
            @test size(obs) == (DIM, DIM)
            @test eltype(obs) == ComplexF64
        end

        # All matrices are Hermitian
        for obs in observables
            @test isapprox(obs, obs'; atol=1e-12)
        end

        # H observable (first) is diagonal in eigenbasis
        H = observables[1]
        @test maximum(abs.(H - diagm(diag(H)))) < 1e-12

        # H Gibbs trace matches analytical
        boltz = exp.(-BETA .* TEST_HAM.eigvals)
        gibbs_energy_analytical = sum(TEST_HAM.eigvals .* boltz) / sum(boltz)
        @test isapprox(real(tr(Matrix(TEST_GIBBS) * observables[1])), gibbs_energy_analytical; atol=1e-10)

        # M_z observable (second) matches inline reference construction
        V = TEST_HAM.eigvecs
        Mz_comp = sum(Matrix{ComplexF64}(pad_term([Z], NUM_QUBITS, i)) for i in 1:NUM_QUBITS) / NUM_QUBITS
        Mz_expected = Matrix{ComplexF64}(V' * Mz_comp * V)
        @test isapprox(observables[2], Mz_expected; atol=1e-14)

        # XX_avg, YY_avg, ZZ_avg (indices 3-5) are per-bond averaged correlations
        V = TEST_HAM.eigvecs
        for (idx, pauli_pair) in [(3, [X, X]), (4, [Y, Y]), (5, [Z, Z])]
            PP_sum = zeros(ComplexF64, DIM, DIM)
            for i in 1:NUM_QUBITS
                PP_sum .+= Matrix{ComplexF64}(pad_term(pauli_pair, NUM_QUBITS, i; periodic=true))
            end
            PP_sum ./= NUM_QUBITS
            PP_expected = Matrix{ComplexF64}(V' * PP_sum * V)
            @test isapprox(observables[idx], PP_expected; atol=1e-14)
        end

        # Mz_stagg observable (index 6) matches inline staggered construction
        V = TEST_HAM.eigvecs
        Mz_stagg_comp = sum((-1)^i .* Matrix{ComplexF64}(pad_term([Z], NUM_QUBITS, i)) for i in 1:NUM_QUBITS) / NUM_QUBITS
        Mz_stagg_expected = Matrix{ComplexF64}(V' * Mz_stagg_comp * V)
        @test isapprox(observables[6], Mz_stagg_expected; atol=1e-14)

        # Z1 observable (index 7) matches single-site construction
        Z1_comp = Matrix{ComplexF64}(pad_term([Z], NUM_QUBITS, 1))
        Z1_expected = Matrix{ComplexF64}(V' * Z1_comp * V)
        @test isapprox(observables[7], Z1_expected; atol=1e-14)
    end

    # -----------------------------------------------------------------------
    # Testset 22: build_preset_trajectory_observables (Trotter basis)
    # -----------------------------------------------------------------------
    @testset "build_preset_trajectory_observables (Trotter basis)" begin
        observables, names = build_preset_trajectory_observables(TEST_HAM, NUM_QUBITS; trotter=TEST_TROTTER)

        # Returns exactly 7 observables and 7 names
        @test length(observables) == 7
        @test length(names) == 7
        @test names == ["H", "Mz", "XX_avg", "YY_avg", "ZZ_avg", "Mz_stagg", "Z1"]

        # All matrices are DIM x DIM ComplexF64
        for obs in observables
            @test size(obs) == (DIM, DIM)
            @test eltype(obs) == ComplexF64
        end

        # All matrices are Hermitian
        for obs in observables
            @test isapprox(obs, obs'; atol=1e-12)
        end

        # H observable (first) is built via full basis transform V_T' * H_data * V_T,
        # NOT as diagm(eigvals) (which would only be valid in the Hamiltonian eigenbasis).
        # For this test system the Trotter error is small so off-diagonal elements may
        # be near zero, but the matrix should still match the explicit transform.
        H = observables[1]
        H_expected = Matrix{ComplexF64}(TEST_TROTTER.eigvecs' * TEST_HAM.data * TEST_TROTTER.eigvecs)
        @test isapprox(H, H_expected; atol=1e-12)

        # H Gibbs trace
        gibbs_trotter = QuantumFurnace._gibbs_in_trotter_basis(TEST_HAM, TEST_TROTTER)
        boltz = exp.(-BETA .* TEST_HAM.eigvals)
        gibbs_energy_analytical = sum(TEST_HAM.eigvals .* boltz) / sum(boltz)
        @test isapprox(real(tr(Matrix(gibbs_trotter) * observables[1])), gibbs_energy_analytical; atol=1e-10)

        # M_z observable (second) matches inline reference construction
        V_T = TEST_TROTTER.eigvecs
        Mz_comp = sum(Matrix{ComplexF64}(pad_term([Z], NUM_QUBITS, i)) for i in 1:NUM_QUBITS) / NUM_QUBITS
        Mz_expected = Matrix{ComplexF64}(V_T' * Mz_comp * V_T)
        @test isapprox(observables[2], Mz_expected; atol=1e-14)

        # Cross-basis consistency: M_z Gibbs trace matches analytical
        gibbs_comp = TEST_HAM.eigvecs * Matrix(TEST_GIBBS) * TEST_HAM.eigvecs'
        mz_analytical = sum(
            real(tr(gibbs_comp * Matrix{ComplexF64}(pad_term([Z], NUM_QUBITS, i))))
            for i in 1:NUM_QUBITS
        ) / NUM_QUBITS
        @test isapprox(real(tr(Matrix(gibbs_trotter) * observables[2])), mz_analytical; atol=1e-10)

        # XX_avg, YY_avg, ZZ_avg (indices 3-5) are per-bond averaged correlations in Trotter basis
        V_T = TEST_TROTTER.eigvecs
        for (idx, pauli_pair) in [(3, [X, X]), (4, [Y, Y]), (5, [Z, Z])]
            PP_sum = zeros(ComplexF64, DIM, DIM)
            for i in 1:NUM_QUBITS
                PP_sum .+= Matrix{ComplexF64}(pad_term(pauli_pair, NUM_QUBITS, i; periodic=true))
            end
            PP_sum ./= NUM_QUBITS
            PP_expected = Matrix{ComplexF64}(V_T' * PP_sum * V_T)
            @test isapprox(observables[idx], PP_expected; atol=1e-14)
        end
    end

end  # @testset "Convergence tracking"
