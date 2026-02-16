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
    # Testset 2: build_convergence_observables (eigenbasis)
    # -----------------------------------------------------------------------
    @testset "build_convergence_observables (eigenbasis)" begin
        observables, names = build_convergence_observables(TEST_HAM, NUM_QUBITS)

        # Correct count: NUM_QUBITS ZZ pairs + 1 H
        @test length(observables) == NUM_QUBITS + 1
        @test length(names) == length(observables)

        # Name pattern check
        @test names == ["ZZ_12", "ZZ_23", "ZZ_34", "ZZ_41", "H"]

        # Matrix dimensions and type
        for obs in observables
            @test size(obs) == (DIM, DIM)
            @test eltype(obs) == ComplexF64
        end

        # ZZ matrices are Hermitian
        for i in 1:NUM_QUBITS
            @test isapprox(observables[i], observables[i]'; atol=1e-12)
        end

        # Energy observable (last one) is diagonal in eigenbasis
        H_eigen = observables[end]
        H_offdiag = H_eigen - diagm(diag(H_eigen))
        @test maximum(abs.(H_offdiag)) < 1e-12

        # Gibbs energy: tr(gibbs * H_eigen) matches analytical
        boltz = exp.(-BETA .* TEST_HAM.eigvals)
        gibbs_energy_analytical = sum(TEST_HAM.eigvals .* boltz) / sum(boltz)
        gibbs_energy_from_obs = real(tr(Matrix(TEST_GIBBS) * H_eigen))
        @test isapprox(gibbs_energy_from_obs, gibbs_energy_analytical; atol=1e-10)
    end

    # -----------------------------------------------------------------------
    # Testset 3: build_convergence_observables_trotter
    # -----------------------------------------------------------------------
    @testset "build_convergence_observables_trotter" begin
        observables_t, names_t = build_convergence_observables_trotter(TEST_HAM, TEST_TROTTER, NUM_QUBITS)

        # Same structure
        @test length(observables_t) == NUM_QUBITS + 1
        @test length(names_t) == length(observables_t)
        @test names_t == ["ZZ_12", "ZZ_23", "ZZ_34", "ZZ_41", "H"]

        # Matrix dimensions and type
        for obs in observables_t
            @test size(obs) == (DIM, DIM)
            @test eltype(obs) == ComplexF64
        end

        # ZZ matrices are Hermitian in Trotter basis
        for i in 1:NUM_QUBITS
            @test isapprox(observables_t[i], observables_t[i]'; atol=1e-12)
        end
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
        observables, _ = build_convergence_observables(TEST_HAM, NUM_QUBITS)
        obs_gibbs = QuantumFurnace._compute_gibbs_observable_values(TEST_GIBBS, observables)

        @test length(obs_gibbs) == length(observables)
        @test all(isfinite, obs_gibbs)
        @test all(x -> isa(x, Real), obs_gibbs)

        # Energy value (last entry) matches analytical Gibbs energy
        boltz = exp.(-BETA .* TEST_HAM.eigvals)
        gibbs_energy_analytical = sum(TEST_HAM.eigvals .* boltz) / sum(boltz)
        @test isapprox(obs_gibbs[end], gibbs_energy_analytical; atol=1e-10)
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
        config = make_thermalize_config(EnergyDomain(); with_coherent=true, delta=0.01, mixing_time=5.0)
        observables, names = build_convergence_observables(TEST_HAM, NUM_QUBITS)
        obs_gibbs = QuantumFurnace._compute_gibbs_observable_values(TEST_GIBBS, observables)

        # Ground state in eigenbasis: first basis vector
        psi0 = zeros(ComplexF64, DIM)
        psi0[1] = 1.0

        traj_result, conv_data = run_trajectories_convergence(
            TEST_JUMPS, config, psi0, TEST_HAM;
            gibbs=TEST_GIBBS,
            observables=observables,
            observable_names=names,
            batch_size=200,
            n_batches=5,
            seed=42,
        )

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
        @test size(conv_data.observable_values) == (NUM_QUBITS + 1, 5)
        @test conv_data.observable_names == names
        @test length(conv_data.observable_gibbs_values) == NUM_QUBITS + 1

        # CONV-02/CONV-03: Energy observable converges toward Gibbs value
        energy_idx = length(names)  # H is the last observable
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
        observables, names = build_convergence_observables(TEST_HAM, NUM_QUBITS)

        psi0 = zeros(ComplexF64, DIM)
        psi0[1] = 1.0

        _, conv1 = run_trajectories_convergence(
            TEST_JUMPS, config, psi0, TEST_HAM;
            gibbs=TEST_GIBBS,
            observables=observables,
            observable_names=names,
            batch_size=50,
            n_batches=2,
            seed=12345,
        )

        _, conv2 = run_trajectories_convergence(
            TEST_JUMPS, config, psi0, TEST_HAM;
            gibbs=TEST_GIBBS,
            observables=observables,
            observable_names=names,
            batch_size=50,
            n_batches=2,
            seed=12345,
        )

        # Bitwise identical
        @test conv1.trace_distances == conv2.trace_distances
        @test conv1.observable_values == conv2.observable_values
    end

    # -----------------------------------------------------------------------
    # Testset 10: Convergence data accessible programmatically (CONV-04)
    # -----------------------------------------------------------------------
    @testset "Convergence data accessible programmatically" begin
        config = make_thermalize_config(EnergyDomain(); with_coherent=true, delta=0.01, mixing_time=5.0)
        observables, names = build_convergence_observables(TEST_HAM, NUM_QUBITS)

        psi0 = zeros(ComplexF64, DIM)
        psi0[1] = 1.0

        _, conv_data = run_trajectories_convergence(
            TEST_JUMPS, config, psi0, TEST_HAM;
            gibbs=TEST_GIBBS,
            observables=observables,
            observable_names=names,
            batch_size=50,
            n_batches=2,
            seed=99,
        )

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
        @test length(col_slice) == NUM_QUBITS + 1

        # observable_names can be used to look up specific observables
        idx = findfirst(==("H"), conv_data.observable_names)
        @test idx !== nothing
        h_curve = conv_data.observable_values[idx, :]
        @test h_curve isa Vector{Float64}
        @test length(h_curve) == 2
    end

end  # @testset "Convergence tracking"
