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

        # Hermitian: basis transform V'*G*V of a Hermitian matrix preserves Hermiticity
        # to machine precision. Error O(DIM^2 * eps) ~ 256 * 2.2e-16 ~ 5.6e-14.
        # Threshold 1e-12 gives ~18x margin.
        herm_err = maximum(abs.(Matrix(gibbs_trotter) - Matrix(gibbs_trotter)'))
        @test isapprox(Matrix(gibbs_trotter), Matrix(gibbs_trotter)'; atol=1e-12)
        @info "Gibbs in Trotter: Hermiticity" max_error=herm_err threshold_atol=1e-12

        # Trace 1 (still a valid density matrix after basis change).
        # Unitary basis transform preserves trace exactly. Error O(DIM * eps) ~ 16 * 2.2e-16 ~ 3.5e-15.
        # Threshold 1e-12 gives ~280x margin.
        gibbs_tr = real(tr(gibbs_trotter))
        @test isapprox(gibbs_tr, 1.0; atol=1e-12)
        @info "Gibbs in Trotter: trace" trace=gibbs_tr deviation=abs(gibbs_tr - 1.0) threshold_atol=1e-12

        # All eigenvalues non-negative (valid density matrix).
        # Eigenvalues of a PSD matrix; numerical error can make smallest eigenvalue slightly negative.
        # Threshold -1e-12 allows for FP rounding.
        eigs = eigvals(Matrix(gibbs_trotter))
        min_eig = minimum(real.(eigs))
        @test all(real.(eigs) .>= -1e-12)
        @info "Gibbs in Trotter: positivity" min_eigenvalue=min_eig threshold=-1e-12
    end

    # -----------------------------------------------------------------------
    # Testset 5: _compute_gibbs_observable_values
    # -----------------------------------------------------------------------
    @testset "_compute_gibbs_observable_values" begin
        observables, obs_names = build_preset_trajectory_observables(TEST_HAM, NUM_QUBITS)
        obs_gibbs = QuantumFurnace._compute_gibbs_observable_values(TEST_GIBBS, observables)

        @test length(obs_gibbs) == length(observables)
        @test all(isfinite, obs_gibbs)
        @test all(x -> isa(x, Real), obs_gibbs)

        # Energy value (H, index 4) matches analytical Gibbs energy.
        # Both computed as Tr[H * rho_gibbs]; one via matrix trace, other via Boltzmann sum.
        # Round-trip error O(DIM * eps) ~ 16 * 2.2e-16 ~ 3.5e-15.
        # Threshold 1e-10 gives ~28,000x margin.
        boltz = exp.(-BETA .* TEST_HAM.eigvals)
        gibbs_energy_analytical = sum(TEST_HAM.eigvals .* boltz) / sum(boltz)
        h_idx = findfirst(==("H"), obs_names)
        energy_err = abs(obs_gibbs[h_idx] - gibbs_energy_analytical)
        @test isapprox(obs_gibbs[h_idx], gibbs_energy_analytical; atol=1e-10)
        @info "Gibbs observable: energy" computed=obs_gibbs[h_idx] analytical=gibbs_energy_analytical error=energy_err threshold_atol=1e-10
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
        config = make_config(Thermalize(), EnergyDomain(); construction=KMS(), delta=0.01, mixing_time=60.0)
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

        # CONV-01: Trace distance decreases (last batch < first batch).
        # With 1000 trajectories and mixing_time=60, the averaged density matrix converges
        # toward the Gibbs state. The trace distance should monotonically decrease on average,
        # though individual batches may fluctuate. Comparing first vs last is a weak convergence check.
        @test conv_data.trace_distances[end] < conv_data.trace_distances[1]
        @info "CONV-01: trace distance convergence" first=conv_data.trace_distances[1] last=conv_data.trace_distances[end]

        # Observable values structure
        @test size(conv_data.observable_values) == (6, 5)
        @test conv_data.observable_names == names
        @test length(conv_data.observable_gibbs_values) == 6

        # CONV-02/CONV-03: Energy observable converges toward Gibbs value.
        # The energy observable <H> should approach the analytical Gibbs energy as
        # N_traj increases. Initial error ~ O(1), final error ~ O(1/sqrt(N_traj)).
        energy_idx = findfirst(==("H"), names)
        initial_energy_err = abs(conv_data.observable_values[energy_idx, 1] - obs_gibbs[energy_idx])
        final_energy_err = abs(conv_data.observable_values[energy_idx, end] - obs_gibbs[energy_idx])
        @test final_energy_err < initial_energy_err
        @info "CONV-02/03: energy convergence" initial_error=initial_energy_err final_error=final_energy_err gibbs_value=obs_gibbs[energy_idx]

        # TrajectoryResult checks
        @test traj_result isa TrajectoryResult
        @test traj_result.n_trajectories == 1000
        @test traj_result.seed == 42
        @test size(traj_result.rho_mean) == (DIM, DIM)

        # rho_mean is approximately a valid density matrix.
        # Trace should be 1.0; with 1000 trajectory average, trace error is O(1/sqrt(1000)) ~ 0.03.
        # Threshold 1e-6 is well within margin. In practice, trace is preserved exactly by the
        # averaging procedure (each trajectory preserves trace), so error is ~ eps * DIM.
        rho_tr = real(tr(traj_result.rho_mean))
        @test isapprox(rho_tr, 1.0; atol=1e-6)
        @info "CONV: rho_mean trace" trace=rho_tr deviation=abs(rho_tr - 1.0) threshold_atol=1e-6

        # Hermiticity of averaged density matrix.
        # Each trajectory produces |psi><psi| which is Hermitian; averaging preserves this.
        # Error O(DIM^2 * eps * N_traj) but dominated by averaging precision ~ eps.
        # Threshold 1e-10 gives large margin.
        herm_err = maximum(abs.(traj_result.rho_mean - traj_result.rho_mean'))
        @test isapprox(traj_result.rho_mean, traj_result.rho_mean'; atol=1e-10)
        @info "CONV: rho_mean Hermiticity" max_error=herm_err threshold_atol=1e-10
    end

    # -----------------------------------------------------------------------
    # Testset 9: run_trajectories_convergence determinism
    # -----------------------------------------------------------------------
    @testset "run_trajectories_convergence determinism" begin
        config = make_config(Thermalize(), EnergyDomain(); construction=KMS(), delta=0.01, mixing_time=5.0)
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
        config = make_config(Thermalize(), EnergyDomain(); construction=KMS(), delta=0.01, mixing_time=5.0)
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
        @test length(col_slice) == 6

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
        # Exact arithmetic on small integers: error O(eps). Threshold 1e-12 is generous.
        wrc_half = QuantumFurnace._windowed_relative_change([1.0, 1.0, 1.0, 0.5, 0.5, 0.5], 3)
        @test isapprox(wrc_half, 0.5; atol=1e-12)
        @info "Windowed relative change: half-step" computed=wrc_half expected=0.5 threshold_atol=1e-12

        # Converged values: all identical -> relative change ~0.0.
        # Exact arithmetic: identical values yield exactly 0 mean difference. Threshold 1e-10.
        wrc_converged = QuantumFurnace._windowed_relative_change([0.1, 0.1, 0.1, 0.1, 0.1, 0.1], 3)
        @test wrc_converged < 1e-10
        @info "Windowed relative change: converged" computed=wrc_converged threshold=1e-10

        # window_size=1: two-element data [0.4, 0.3]
        # mean_previous = 0.4, mean_recent = 0.3
        # relative_change = abs(0.3 - 0.4) / max(0.4, eps()) = 0.1/0.4 = 0.25
        # Exact arithmetic. Threshold 1e-12.
        wrc_quarter = QuantumFurnace._windowed_relative_change([0.4, 0.3], 1)
        @test isapprox(wrc_quarter, 0.25; atol=1e-12)
        @info "Windowed relative change: quarter" computed=wrc_quarter expected=0.25 threshold_atol=1e-12

        # Edge: window_size=2 with exactly 4 elements [2.0, 2.0, 1.0, 1.0]
        # mean_previous = 2.0, mean_recent = 1.0, relative_change = 1.0/2.0 = 0.5
        # Exact arithmetic. Threshold 1e-12.
        wrc_edge = QuantumFurnace._windowed_relative_change([2.0, 2.0, 1.0, 1.0], 2)
        @test isapprox(wrc_edge, 0.5; atol=1e-12)
        @info "Windowed relative change: edge case" computed=wrc_edge expected=0.5 threshold_atol=1e-12
    end

    # -----------------------------------------------------------------------
    # Testset 13: run_trajectories_adaptive convergence (CONV-04)
    # -----------------------------------------------------------------------
    @testset "run_trajectories_adaptive convergence (CONV-04)" begin
        config = make_config(Thermalize(), EnergyDomain(); construction=KMS(), delta=0.01, mixing_time=5.0)
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

        # System should converge with generous 5% threshold
        @test conv_data.converged == true

        # Converged before hitting hard cap: cld(20000, 200) = 100.
        # With 200-traj batches and 5% threshold on a 4-qubit system, convergence
        # typically takes 10-30 batches. Cap of 100 gives ample margin.
        @test conv_data.total_batches < 100
        @info "CONV-04: adaptive convergence" total_batches=conv_data.total_batches cap=100

        # Patience was met
        @test conv_data.consecutive_stable_batches >= 3

        # Final relative change is below threshold.
        # This is the windowed relative change that triggered convergence.
        @test conv_data.final_relative_change < 0.05
        @test isfinite(conv_data.final_relative_change)
        @info "CONV-04: final relative change" value=conv_data.final_relative_change threshold=0.05

        # Trace distances are all positive and finite
        @test all(isfinite, conv_data.trace_distances)
        @test all(x -> x >= 0, conv_data.trace_distances)

        # Observable values matrix shape matches
        @test size(conv_data.observable_values) == (6, conv_data.total_batches)

        # TrajectoryResult has valid density matrix.
        # Trace preserved by trajectory averaging; error ~ eps * DIM. Threshold 1e-6.
        rho_tr = real(tr(traj_result.rho_mean))
        @test isapprox(rho_tr, 1.0; atol=1e-6)
        @info "CONV-04: rho_mean trace" trace=rho_tr threshold_atol=1e-6

        # Hermiticity: |psi><psi| averaging preserves Hermiticity. Threshold 1e-10.
        herm_err = maximum(abs.(traj_result.rho_mean - traj_result.rho_mean'))
        @test isapprox(traj_result.rho_mean, traj_result.rho_mean'; atol=1e-10)
        @info "CONV-04: rho_mean Hermiticity" max_error=herm_err threshold_atol=1e-10

        # Total trajectories match batch accounting
        @test traj_result.n_trajectories == conv_data.total_batches * 200
    end

    # -----------------------------------------------------------------------
    # Testset 14: run_trajectories_adaptive hard cap (CONV-05)
    # -----------------------------------------------------------------------
    @testset "run_trajectories_adaptive hard cap (CONV-05)" begin
        config = make_config(Thermalize(), EnergyDomain(); construction=KMS(), delta=0.01, mixing_time=5.0)
        observables, names = build_preset_trajectory_observables(TEST_HAM, NUM_QUBITS)

        psi0 = zeros(ComplexF64, DIM)
        psi0[1] = 1.0

        # n_max=500, batch_size=50 -> max_batches=10, threshold=0.0001 is nearly impossible.
        # With only 500 trajectories and 0.01% convergence threshold, the system cannot converge.
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

        # Patience NOT met (couldn't sustain 3 consecutive stable batches at 0.01% threshold)
        @test conv_data.consecutive_stable_batches < 3
        @info "CONV-05: hard cap" total_batches=conv_data.total_batches consecutive_stable=conv_data.consecutive_stable_batches

        # All trace distances positive and finite
        @test all(isfinite, conv_data.trace_distances)
        @test all(x -> x >= 0, conv_data.trace_distances)
    end

    # -----------------------------------------------------------------------
    # Testset 15: run_trajectories_adaptive determinism
    # -----------------------------------------------------------------------
    @testset "run_trajectories_adaptive determinism" begin
        config = make_config(Thermalize(), EnergyDomain(); construction=KMS(), delta=0.01, mixing_time=5.0)
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
        config = make_config(Thermalize(), EnergyDomain(); construction=KMS(), delta=0.01, mixing_time=5.0)
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

        # Returns exactly 6 observables and 6 names
        @test length(observables) == 6
        @test length(names) == 6
        @test names == ["Z1", "X1", "Z1_Zhalf", "H", "Rand_traceless", "Mz_stagg"]

        # All matrices are DIM x DIM ComplexF64
        for obs in observables
            @test size(obs) == (DIM, DIM)
            @test eltype(obs) == ComplexF64
        end

        # All matrices are Hermitian.
        # Unitary basis transform V'*O*V of Hermitian O preserves Hermiticity to machine precision.
        # Error O(DIM^2 * eps) ~ 256 * 2.2e-16 ~ 5.6e-14. Threshold 1e-12 gives ~18x margin.
        max_herm_err = 0.0
        for obs in observables
            err = maximum(abs.(obs - obs'))
            max_herm_err = max(max_herm_err, err)
            @test isapprox(obs, obs'; atol=1e-12)
        end
        @info "Observable Hermiticity (eigenbasis)" max_error=max_herm_err threshold_atol=1e-12 n_observables=6

        V = TEST_HAM.eigvecs

        # Z1 observable (index 1) matches single-site construction.
        # V'*Z1_comp*V vs direct construction: round-trip through unitary transform.
        # Error O(DIM^2 * eps) ~ 5.6e-14. Threshold 1e-14 is tight but passes because
        # both sides compute exactly the same matrix product.
        Z1_comp = Matrix{ComplexF64}(pad_term([Z], NUM_QUBITS, 1))
        Z1_expected = Matrix{ComplexF64}(V' * Z1_comp * V)
        z1_err = maximum(abs.(observables[1] - Z1_expected))
        @test isapprox(observables[1], Z1_expected; atol=1e-14)
        @info "Observable Z1 match" max_error=z1_err threshold_atol=1e-14

        # X1 observable (index 2) matches single-site X construction.
        # Same error analysis as Z1.
        X1_comp = Matrix{ComplexF64}(pad_term([X], NUM_QUBITS, 1))
        X1_expected = Matrix{ComplexF64}(V' * X1_comp * V)
        x1_err = maximum(abs.(observables[2] - X1_expected))
        @test isapprox(observables[2], X1_expected; atol=1e-14)
        @info "Observable X1 match" max_error=x1_err threshold_atol=1e-14

        # Z1_Zhalf observable (index 3) matches two-point correlator.
        # Two matrix multiplications add one more O(DIM * eps) factor.
        # Error O(DIM^3 * eps) ~ 4096 * 2.2e-16 ~ 9e-13. Threshold 1e-14 is tight
        # but passes because pad_term product is exact on sparse Pauli kronecker products.
        half_site = floor(Int, NUM_QUBITS / 2)
        Z1_Zhalf_comp = Matrix{ComplexF64}(pad_term([Z], NUM_QUBITS, 1)) *
                        Matrix{ComplexF64}(pad_term([Z], NUM_QUBITS, half_site))
        Z1_Zhalf_expected = Matrix{ComplexF64}(V' * Z1_Zhalf_comp * V)
        zz_err = maximum(abs.(observables[3] - Z1_Zhalf_expected))
        @test isapprox(observables[3], Z1_Zhalf_expected; atol=1e-14)
        @info "Observable Z1_Zhalf match" max_error=zz_err threshold_atol=1e-14

        # H observable (index 4) is diagonal in eigenbasis.
        # H = V'*H_data*V = diag(eigvals) by construction.
        # Off-diagonal elements should be exactly zero up to FP accumulation.
        # Error O(DIM^2 * eps) ~ 5.6e-14. Threshold 1e-12 gives ~18x margin.
        H = observables[4]
        h_offdiag = maximum(abs.(H - diagm(diag(H))))
        @test h_offdiag < 1e-12
        @info "Observable H diagonality" max_offdiag=h_offdiag threshold=1e-12

        # H Gibbs trace matches analytical.
        # Tr[rho_gibbs * H] computed via matrix vs Boltzmann sum.
        # Round-trip error O(DIM * eps) ~ 3.5e-15. Threshold 1e-10 gives ~28,000x margin.
        boltz = exp.(-BETA .* TEST_HAM.eigvals)
        gibbs_energy_analytical = sum(TEST_HAM.eigvals .* boltz) / sum(boltz)
        h_gibbs_trace = real(tr(Matrix(TEST_GIBBS) * observables[4]))
        h_energy_err = abs(h_gibbs_trace - gibbs_energy_analytical)
        @test isapprox(h_gibbs_trace, gibbs_energy_analytical; atol=1e-10)
        @info "Observable H Gibbs energy" computed=h_gibbs_trace analytical=gibbs_energy_analytical error=h_energy_err threshold_atol=1e-10

        # Rand_traceless observable (index 5): Hermitian, traceless, operator-norm 1, reproducible
        Rand = observables[5]
        @test isapprox(Rand, Rand'; atol=1e-12)  # Hermitian (already checked above, but explicit)
        # Tracelessness in eigenbasis: tr(V' * R * V) = tr(R) since V is unitary.
        # Trace is O(eps) for a traceless matrix. Threshold 1e-10.
        rand_tr = abs(tr(Rand))
        @test rand_tr < 1e-10
        @info "Observable Rand_traceless: trace" abs_trace=rand_tr threshold=1e-10
        # Operator norm 1 (in computational basis, preserved by unitary transform).
        # opnorm computed via SVD; error O(DIM * eps) ~ 3.5e-15. Threshold 1e-10.
        rand_opnorm = opnorm(Rand)
        @test isapprox(rand_opnorm, 1.0; atol=1e-10)
        @info "Observable Rand_traceless: opnorm" opnorm=rand_opnorm threshold_atol=1e-10
        # Reproducibility: calling again gives the same matrix (seeded RNG).
        # Should be bitwise identical. Threshold 1e-14 allows for any platform variation.
        obs2, _ = build_preset_trajectory_observables(TEST_HAM, NUM_QUBITS)
        repro_err = maximum(abs.(observables[5] - obs2[5]))
        @test isapprox(observables[5], obs2[5]; atol=1e-14)
        @info "Observable Rand_traceless: reproducibility" max_error=repro_err threshold_atol=1e-14

        # Mz_stagg observable (index 6) matches inline staggered construction.
        # Sum of DIM single-site Pauli terms, each with unitary transform.
        # Error O(NUM_QUBITS * DIM^2 * eps) ~ 4 * 256 * 2.2e-16 ~ 2.3e-13.
        # Threshold 1e-14 is tight but passes because pad_term is exact on sparse Paulis.
        Mz_stagg_comp = sum((-1)^i .* Matrix{ComplexF64}(pad_term([Z], NUM_QUBITS, i)) for i in 1:NUM_QUBITS) / NUM_QUBITS
        Mz_stagg_expected = Matrix{ComplexF64}(V' * Mz_stagg_comp * V)
        mz_err = maximum(abs.(observables[6] - Mz_stagg_expected))
        @test isapprox(observables[6], Mz_stagg_expected; atol=1e-14)
        @info "Observable Mz_stagg match" max_error=mz_err threshold_atol=1e-14
    end

    # -----------------------------------------------------------------------
    # Testset 22: build_preset_trajectory_observables (Trotter basis)
    # -----------------------------------------------------------------------
    @testset "build_preset_trajectory_observables (Trotter basis)" begin
        observables, names = build_preset_trajectory_observables(TEST_HAM, NUM_QUBITS; trotter=TEST_TROTTER)

        # Returns exactly 6 observables and 6 names
        @test length(observables) == 6
        @test length(names) == 6
        @test names == ["Z1", "X1", "Z1_Zhalf", "H", "Rand_traceless", "Mz_stagg"]

        # All matrices are DIM x DIM ComplexF64
        for obs in observables
            @test size(obs) == (DIM, DIM)
            @test eltype(obs) == ComplexF64
        end

        # All matrices are Hermitian.
        # Same error analysis as eigenbasis: O(DIM^2 * eps) ~ 5.6e-14. Threshold 1e-12.
        max_herm_err = 0.0
        for obs in observables
            err = maximum(abs.(obs - obs'))
            max_herm_err = max(max_herm_err, err)
            @test isapprox(obs, obs'; atol=1e-12)
        end
        @info "Observable Hermiticity (Trotter)" max_error=max_herm_err threshold_atol=1e-12 n_observables=6

        V_T = TEST_TROTTER.eigvecs

        # Z1 observable (index 1) matches single-site Z in Trotter basis.
        # V_T'*Z1_comp*V_T: same error analysis as eigenbasis. Threshold 1e-14.
        Z1_comp = Matrix{ComplexF64}(pad_term([Z], NUM_QUBITS, 1))
        Z1_expected = Matrix{ComplexF64}(V_T' * Z1_comp * V_T)
        z1_err = maximum(abs.(observables[1] - Z1_expected))
        @test isapprox(observables[1], Z1_expected; atol=1e-14)
        @info "Observable Z1 Trotter match" max_error=z1_err threshold_atol=1e-14

        # X1 observable (index 2) matches single-site X in Trotter basis
        X1_comp = Matrix{ComplexF64}(pad_term([X], NUM_QUBITS, 1))
        X1_expected = Matrix{ComplexF64}(V_T' * X1_comp * V_T)
        x1_err = maximum(abs.(observables[2] - X1_expected))
        @test isapprox(observables[2], X1_expected; atol=1e-14)
        @info "Observable X1 Trotter match" max_error=x1_err threshold_atol=1e-14

        # Z1_Zhalf observable (index 3) matches two-point correlator in Trotter basis
        half_site = floor(Int, NUM_QUBITS / 2)
        Z1_Zhalf_comp = Matrix{ComplexF64}(pad_term([Z], NUM_QUBITS, 1)) *
                        Matrix{ComplexF64}(pad_term([Z], NUM_QUBITS, half_site))
        Z1_Zhalf_expected = Matrix{ComplexF64}(V_T' * Z1_Zhalf_comp * V_T)
        zz_err = maximum(abs.(observables[3] - Z1_Zhalf_expected))
        @test isapprox(observables[3], Z1_Zhalf_expected; atol=1e-14)
        @info "Observable Z1_Zhalf Trotter match" max_error=zz_err threshold_atol=1e-14

        # H observable (index 4) is built via full basis transform V_T' * H_data * V_T,
        # NOT as diagm(eigvals) (which would only be valid in the Hamiltonian eigenbasis).
        # Two matrix multiplications: error O(DIM^3 * eps) ~ 4096 * 2.2e-16 ~ 9e-13.
        # Threshold 1e-12 gives ~1.1x margin.
        H = observables[4]
        H_expected = Matrix{ComplexF64}(TEST_TROTTER.eigvecs' * TEST_HAM.data * TEST_TROTTER.eigvecs)
        h_err = maximum(abs.(H - H_expected))
        @test isapprox(H, H_expected; atol=1e-12)
        @info "Observable H Trotter match" max_error=h_err threshold_atol=1e-12

        # H Gibbs trace matches analytical (basis-independent quantity).
        # Tr[rho_gibbs * H] is independent of basis choice.
        gibbs_trotter = QuantumFurnace._gibbs_in_trotter_basis(TEST_HAM, TEST_TROTTER)
        boltz = exp.(-BETA .* TEST_HAM.eigvals)
        gibbs_energy_analytical = sum(TEST_HAM.eigvals .* boltz) / sum(boltz)
        h_gibbs_trace = real(tr(Matrix(gibbs_trotter) * observables[4]))
        h_energy_err = abs(h_gibbs_trace - gibbs_energy_analytical)
        @test isapprox(h_gibbs_trace, gibbs_energy_analytical; atol=1e-10)
        @info "Observable H Trotter Gibbs energy" computed=h_gibbs_trace analytical=gibbs_energy_analytical error=h_energy_err threshold_atol=1e-10

        # Rand_traceless observable (index 5): Hermitian, traceless, operator-norm 1
        Rand = observables[5]
        @test isapprox(Rand, Rand'; atol=1e-12)
        rand_tr = abs(tr(Rand))
        @test rand_tr < 1e-10
        @info "Observable Rand_traceless Trotter: trace" abs_trace=rand_tr threshold=1e-10
        rand_opnorm = opnorm(Rand)
        @test isapprox(rand_opnorm, 1.0; atol=1e-10)
        @info "Observable Rand_traceless Trotter: opnorm" opnorm=rand_opnorm threshold_atol=1e-10

        # Mz_stagg observable (index 6) matches inline staggered construction in Trotter basis
        Mz_stagg_comp = sum((-1)^i .* Matrix{ComplexF64}(pad_term([Z], NUM_QUBITS, i)) for i in 1:NUM_QUBITS) / NUM_QUBITS
        Mz_stagg_expected = Matrix{ComplexF64}(V_T' * Mz_stagg_comp * V_T)
        mz_err = maximum(abs.(observables[6] - Mz_stagg_expected))
        @test isapprox(observables[6], Mz_stagg_expected; atol=1e-14)
        @info "Observable Mz_stagg Trotter match" max_error=mz_err threshold_atol=1e-14
    end

end  # @testset "Convergence tracking"
