"""
ExperimentResult serialization round-trip tests (Phase 15, Plan 02).

Tests cover:
- Round-trip save/load for all 4 config variants (Config{Thermalize,KMS}, Config{Thermalize,GNS}, Config{Lindbladian,KMS}, Config{Lindbladian,GNS})
- TrajectoryResult with observables (times, measurements_mean)
- Forward compatibility (missing fields in BSON)
- Companion .txt file generation
- Filename generation
- Integration test with real trajectory run
- Metadata auto-capture
"""

using Random
using BSON

@testset "ExperimentResult serialization" begin

    # -----------------------------------------------------------------------
    # Helper: create a fake TrajectoryResult for unit tests
    # -----------------------------------------------------------------------
    function _make_fake_trajectory(dim; seed=42, n_trajectories=100)
        rho = Matrix(random_density_matrix(Int(log2(dim))))
        TrajectoryResult(rho, n_trajectories, seed, nothing, nothing, nothing)
    end

    # -----------------------------------------------------------------------
    # Test group 1: KMS Thermalize round-trip
    # -----------------------------------------------------------------------
    @testset "KMS Thermalize round-trip" begin
        mktempdir() do tmpdir
            config = make_small_thermalize_config(EnergyDomain())
            traj = _make_fake_trajectory(SMALL_DIM; seed=42, n_trajectories=100)
            ham_params = QuantumFurnace._extract_hamiltonian_params(SMALL_HAM)
            metadata = Dict{Symbol, Any}(
                :seed => 42,
                :n_threads => 4,
                :julia_version => "test",
                :timestamp => "2026-02-16",
                :git_hash => "abc123",
                :wall_time_seconds => 1.5,
            )

            result = ExperimentResult(config, traj, ham_params, metadata)
            path = joinpath(tmpdir, "test_kms.bson")
            save_experiment(result, path)

            loaded = load_experiment(path)

            # Config fields
            @test loaded.config.num_qubits == config.num_qubits
            @test loaded.config.beta == config.beta
            @test loaded.config.sigma == config.sigma
            @test with_coherent(loaded.config.construction) == with_coherent(config.construction)
            @test loaded.config.domain isa EnergyDomain
            @test loaded.config isa Config{Thermalize}
            @test loaded.config.mixing_time == config.mixing_time
            @test loaded.config.delta == config.delta

            # Trajectory result
            @test isapprox(loaded.trajectory_result.rho_mean, traj.rho_mean; atol=0)
            @test loaded.trajectory_result.n_trajectories == 100
            @test loaded.trajectory_result.seed == 42

            # Metadata
            @test loaded.metadata[:seed] == 42
            @test loaded.metadata[:n_threads] == 4
            @test loaded.metadata[:julia_version] == "test"
            @test loaded.metadata[:timestamp] == "2026-02-16"
            @test loaded.metadata[:git_hash] == "abc123"
            @test loaded.metadata[:wall_time_seconds] == 1.5

            # Hamiltonian params
            @test haskey(loaded.hamiltonian_params, :base_coeffs)
            @test haskey(loaded.hamiltonian_params, :periodic)
        end
    end

    # -----------------------------------------------------------------------
    # Test group 2: GNS Thermalize round-trip
    # -----------------------------------------------------------------------
    @testset "GNS Thermalize round-trip" begin
        mktempdir() do tmpdir
            config = make_small_thermalize_config_gns(EnergyDomain())
            traj = _make_fake_trajectory(SMALL_DIM)
            ham_params = QuantumFurnace._extract_hamiltonian_params(SMALL_HAM)
            metadata = Dict{Symbol, Any}(:seed => 1)

            result = ExperimentResult(config, traj, ham_params, metadata)
            path = joinpath(tmpdir, "test_gns.bson")
            save_experiment(result, path)

            loaded = load_experiment(path)

            # Key assertion: GNS type preserved
            @test loaded.config isa Config{Thermalize, <:Any, GNS}
            @test with_coherent(loaded.config.construction) == false
            @test loaded.config.domain isa EnergyDomain
            @test loaded.config.beta == config.beta
            @test loaded.config.mixing_time == config.mixing_time

            # Trajectory round-trip
            @test isapprox(loaded.trajectory_result.rho_mean, traj.rho_mean; atol=0)
        end
    end

    # -----------------------------------------------------------------------
    # Test group 3: Lindbladian round-trip
    # -----------------------------------------------------------------------
    @testset "Lindbladian round-trip" begin
        mktempdir() do tmpdir
            config = make_small_liouv_config(EnergyDomain())
            traj = _make_fake_trajectory(SMALL_DIM)
            ham_params = Dict{Symbol, Any}()
            metadata = Dict{Symbol, Any}()

            result = ExperimentResult(config, traj, ham_params, metadata)
            path = joinpath(tmpdir, "test_liouv.bson")
            save_experiment(result, path)

            loaded = load_experiment(path)

            # Key assertion: Lindbladian config (not thermalize)
            @test loaded.config isa Config{Lindbladian}
            @test !(loaded.config.sim isa Thermalize)
            @test loaded.config.num_qubits == config.num_qubits
            @test loaded.config.beta == config.beta
            @test loaded.config.domain isa EnergyDomain

            # Trajectory round-trip
            @test isapprox(loaded.trajectory_result.rho_mean, traj.rho_mean; atol=0)
        end
    end

    # -----------------------------------------------------------------------
    # Test group 4: Lindbladian GNS round-trip
    # -----------------------------------------------------------------------
    @testset "Lindbladian GNS round-trip" begin
        mktempdir() do tmpdir
            config = make_small_liouv_config_gns(EnergyDomain())
            traj = _make_fake_trajectory(SMALL_DIM)
            ham_params = Dict{Symbol, Any}()
            metadata = Dict{Symbol, Any}()

            result = ExperimentResult(config, traj, ham_params, metadata)
            path = joinpath(tmpdir, "test_liouv_gns.bson")
            save_experiment(result, path)

            loaded = load_experiment(path)

            # Key assertion: GNS type preserved for Lindbladian config
            @test loaded.config isa Config{Lindbladian, <:Any, GNS}
            @test with_coherent(loaded.config.construction) == false
            @test loaded.config.domain isa EnergyDomain

            # Trajectory round-trip
            @test isapprox(loaded.trajectory_result.rho_mean, traj.rho_mean; atol=0)
        end
    end

    # -----------------------------------------------------------------------
    # Test group 5: TrajectoryResult with observables
    # -----------------------------------------------------------------------
    @testset "TrajectoryResult with observables" begin
        mktempdir() do tmpdir
            config = make_small_thermalize_config(EnergyDomain())
            rho = Matrix(random_density_matrix(3))
            times = collect(0.0:0.1:1.0)
            measurements = rand(length(times), 2)
            traj = TrajectoryResult(rho, 50, 99, times, measurements, nothing)
            ham_params = Dict{Symbol, Any}()
            metadata = Dict{Symbol, Any}()

            result = ExperimentResult(config, traj, ham_params, metadata)
            path = joinpath(tmpdir, "test_obs.bson")
            save_experiment(result, path)

            loaded = load_experiment(path)

            @test loaded.trajectory_result.times !== nothing
            @test loaded.trajectory_result.measurements_mean !== nothing
            @test isapprox(loaded.trajectory_result.times, times; atol=0)
            @test isapprox(loaded.trajectory_result.measurements_mean, measurements; atol=0)
            @test loaded.trajectory_result.n_trajectories == 50
            @test loaded.trajectory_result.seed == 99
        end
    end

    # -----------------------------------------------------------------------
    # Test group 6: Forward compatibility - missing fields
    # -----------------------------------------------------------------------
    @testset "Forward compatibility - missing fields" begin
        mktempdir() do tmpdir
            # Manually create a minimal BSON Dict missing :metadata and :hamiltonian_params
            config_d = Dict{Symbol, Any}(
                :config_type => "KMS",
                :config_kind => "thermalize",
                :domain => "EnergyDomain",
                :num_qubits => 3,
                :with_coherent => false,
                :with_linear_combination => true,
                :beta => 10.0,
                :sigma => 0.1,
                :gaussian_parameters => [nothing, nothing],
                :a => nothing,
                :b => nothing,
                :num_energy_bits => nothing,
                :t0 => nothing,
                :w0 => nothing,
                :eta => nothing,
                :num_trotter_steps_per_t0 => nothing,
                :mixing_time => 1.0,
                :delta => 0.01,
            )

            rho = Matrix(random_density_matrix(3))
            traj_d = Dict{Symbol, Any}(
                :rho_mean => rho,
                :n_trajectories => 10,
                :seed => 1,
                :times => nothing,
                :measurements_mean => nothing,
            )

            d = Dict{Symbol, Any}(
                :config => config_d,
                :trajectory => traj_d,
                # NOTE: :metadata and :hamiltonian_params are intentionally missing
            )

            path = joinpath(tmpdir, "test_missing.bson")
            BSON.bson(path, d)

            loaded = load_experiment(path)

            # Should not error -- metadata and hamiltonian_params should be empty Dicts
            @test loaded.metadata isa Dict
            @test isempty(loaded.metadata)
            @test loaded.hamiltonian_params isa Dict
            @test isempty(loaded.hamiltonian_params)
            @test loaded.config isa Config{Thermalize}
            @test loaded.config.num_qubits == 3
        end
    end

    # -----------------------------------------------------------------------
    # Test group 7: Companion .txt file
    # -----------------------------------------------------------------------
    @testset "Companion .txt file" begin
        mktempdir() do tmpdir
            config = make_small_thermalize_config(EnergyDomain())
            traj = _make_fake_trajectory(SMALL_DIM; seed=42, n_trajectories=100)
            ham_params = QuantumFurnace._extract_hamiltonian_params(SMALL_HAM)
            metadata = Dict{Symbol, Any}(
                :seed => 42,
                :n_threads => 4,
                :julia_version => string(VERSION),
                :timestamp => "2026-02-16_12:00:00",
                :git_hash => "abc123",
                :wall_time_seconds => 1.5,
            )

            result = ExperimentResult(config, traj, ham_params, metadata)
            bson_path = joinpath(tmpdir, "test_companion.bson")
            save_experiment(result, bson_path)

            txt_path = replace(bson_path, ".bson" => ".txt")
            @test isfile(txt_path)

            content = read(txt_path, String)
            @test occursin("QuantumFurnace", content)
            @test occursin("EnergyDomain", content)
            @test occursin("3", content)  # num_qubits
        end
    end

    # -----------------------------------------------------------------------
    # Test group 8: Filename generation
    # -----------------------------------------------------------------------
    @testset "Filename generation" begin
        # KMS TrotterDomain config
        kms_config = make_small_thermalize_config(TrotterDomain(); construction=KMS())
        kms_filename = QuantumFurnace._generate_experiment_filename(kms_config)
        @test startswith(kms_filename, "kms_")
        @test occursin("n3", kms_filename)
        @test occursin("beta10", kms_filename)
        @test occursin("trotter", kms_filename)
        @test endswith(kms_filename, ".bson")

        # GNS config
        gns_config = make_small_thermalize_config_gns(EnergyDomain())
        gns_filename = QuantumFurnace._generate_experiment_filename(gns_config)
        @test startswith(gns_filename, "gns_")
        @test endswith(gns_filename, ".bson")
    end

    # -----------------------------------------------------------------------
    # Test group 9: Integration - save/load real trajectory result
    # -----------------------------------------------------------------------
    @testset "Integration: save/load real trajectory result" begin
        mktempdir() do tmpdir
            # Use SMALL test fixtures (3-qubit system)
            config = make_small_thermalize_config(TrotterDomain(); delta=0.01, mixing_time=1.0)
            psi0 = zeros(ComplexF64, SMALL_DIM)
            psi0[1] = 1.0  # computational basis |0>

            traj_result = run_trajectories(
                SMALL_JUMPS, config, psi0, SMALL_HAM;
                trotter=SMALL_TROTTER, ntraj=10, seed=12345,
            )

            ham_params = QuantumFurnace._extract_hamiltonian_params(SMALL_HAM)
            metadata = QuantumFurnace._capture_metadata(wall_time_seconds=0.5)

            result = ExperimentResult(config, traj_result, ham_params, metadata)
            path = joinpath(tmpdir, "test_integration.bson")
            save_experiment(result, path)

            loaded = load_experiment(path)

            # Exact density matrix match (serialization only, no float computation)
            @test isapprox(loaded.trajectory_result.rho_mean, traj_result.rho_mean; atol=0)
            @test loaded.trajectory_result.n_trajectories == 10
            @test loaded.trajectory_result.seed == 12345

            # Config round-trip
            @test loaded.config.beta == config.beta
            @test loaded.config.sigma == config.sigma
            @test loaded.config.domain isa TrotterDomain
            @test loaded.config.num_qubits == 3

            # Metadata has auto-captured fields (julia_version removed per Phase 36 decision)
            @test !haskey(loaded.metadata, :julia_version)
            @test haskey(loaded.metadata, :timestamp)
            @test haskey(loaded.metadata, :git_hash)

            # Hamiltonian params have reproduction-relevant keys
            @test haskey(loaded.hamiltonian_params, :base_coeffs)
            @test haskey(loaded.hamiltonian_params, :periodic)
        end
    end

    # -----------------------------------------------------------------------
    # Test group 10: Metadata auto-capture
    # -----------------------------------------------------------------------
    @testset "Metadata auto-capture" begin
        meta = QuantumFurnace._capture_metadata(n_threads=2, wall_time_seconds=3.14)

        # julia_version removed per Phase 36 locked decision
        @test !haskey(meta, :julia_version)
        @test haskey(meta, :timestamp)
        @test haskey(meta, :git_hash)
        @test haskey(meta, :n_threads)
        @test haskey(meta, :wall_time_seconds)

        @test meta[:n_threads] == 2
        @test meta[:wall_time_seconds] == 3.14
        @test meta[:git_hash] isa String
        @test !isempty(meta[:timestamp])
    end

end  # @testset "ExperimentResult serialization"
