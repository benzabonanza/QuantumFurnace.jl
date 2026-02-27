using Random
using BSON

@testset "New Result types serialization" begin

    # -----------------------------------------------------------------------
    # LindbladResults round-trip
    # -----------------------------------------------------------------------
    @testset "LindbladResults round-trip" begin
        mktempdir() do tmpdir
            config = make_small_liouv_config(EnergyDomain())
            dim = SMALL_DIM
            eigenvalues = [0.0 + 0.0im, -0.5 + 0.01im]
            fixed_point = Matrix(random_density_matrix(Int(log2(dim))))
            gap_mode = randn(ComplexF64, dim, dim)
            spectral_gap = -0.5 + 0.01im
            metadata = Dict{Symbol, Any}(:wall_time_seconds => 1.5, :n_threads => 2,
                :timestamp => "2026-02-27", :git_hash => "abc123")

            result = LindbladResults{Float64}(config, eigenvalues, fixed_point, gap_mode, spectral_gap, metadata)
            path = joinpath(tmpdir, "test_lindblad.bson")
            save_result(result, path)

            loaded = load_result(path)
            @test loaded isa LindbladResults
            @test loaded.config.num_qubits == config.num_qubits
            @test loaded.config.beta == config.beta
            @test loaded.config.domain isa EnergyDomain
            @test isapprox(loaded.eigenvalues, eigenvalues; atol=0)
            @test isapprox(loaded.fixed_point, fixed_point; atol=0)
            @test isapprox(loaded.gap_mode, gap_mode; atol=0)
            @test loaded.spectral_gap == spectral_gap
            @test loaded.metadata[:wall_time_seconds] == 1.5

            # Companion .txt exists
            txt_path = replace(path, ".bson" => ".txt")
            @test isfile(txt_path)
            content = read(txt_path, String)
            @test occursin("LindbladResults", content)
            @test occursin("EnergyDomain", content)
        end
    end

    # -----------------------------------------------------------------------
    # ThermalizeResults round-trip
    # -----------------------------------------------------------------------
    @testset "ThermalizeResults round-trip" begin
        mktempdir() do tmpdir
            config = make_small_thermalize_config(EnergyDomain())
            dim = SMALL_DIM
            final_dm = Matrix(random_density_matrix(Int(log2(dim))))
            trace_distances = [0.5, 0.3, 0.1, 0.05]
            time_steps = [0.0, 0.01, 0.02, 0.03]
            metadata = Dict{Symbol, Any}(:wall_time_seconds => 2.0, :n_threads => 1,
                :timestamp => "2026-02-27", :git_hash => "def456")

            result = ThermalizeResults{Float64}(config, final_dm, trace_distances, time_steps, metadata)
            path = joinpath(tmpdir, "test_thermalize.bson")
            save_result(result, path)

            loaded = load_result(path)
            @test loaded isa ThermalizeResults
            @test loaded.config isa Config{Thermalize}
            @test loaded.config.mixing_time == config.mixing_time
            @test loaded.config.delta == config.delta
            @test isapprox(loaded.final_dm, final_dm; atol=0)
            @test isapprox(loaded.trace_distances, trace_distances; atol=0)
            @test isapprox(loaded.time_steps, time_steps; atol=0)
            @test loaded.metadata[:wall_time_seconds] == 2.0

            txt_path = replace(path, ".bson" => ".txt")
            @test isfile(txt_path)
        end
    end

    # -----------------------------------------------------------------------
    # KrylovSpectrumResults round-trip
    # -----------------------------------------------------------------------
    @testset "KrylovSpectrumResults round-trip" begin
        mktempdir() do tmpdir
            config = make_small_liouv_config(EnergyDomain())
            dim = SMALL_DIM
            eigenvalues = [0.0+0.0im, -0.3+0.0im, -0.5+0.01im, -0.8+0.0im]
            spectral_gap = 0.3
            fixed_point = Matrix(random_density_matrix(Int(log2(dim))))
            gap_mode = randn(ComplexF64, dim, dim)
            metadata = Dict{Symbol, Any}(:wall_time_seconds => 5.0, :n_threads => 4,
                :timestamp => "2026-02-27", :git_hash => "ghi789")

            result = KrylovSpectrumResults{Float64}(
                config, eigenvalues, spectral_gap, fixed_point, gap_mode,
                4, 100, 2, [1e-11, 1e-10, 1e-9, 1e-8],
                nothing, nothing, metadata,
            )
            path = joinpath(tmpdir, "test_krylov.bson")
            save_result(result, path)

            loaded = load_result(path)
            @test loaded isa KrylovSpectrumResults
            @test isapprox(loaded.eigenvalues, eigenvalues; atol=0)
            @test loaded.spectral_gap == spectral_gap
            @test isapprox(loaded.fixed_point, fixed_point; atol=0)
            @test loaded.converged == 4
            @test loaded.matvec_count == 100
            @test loaded.num_restarts == 2
            @test loaded.channel_eigenvalues === nothing
            @test loaded.delta_used === nothing
            @test loaded.metadata[:wall_time_seconds] == 5.0

            txt_path = replace(path, ".bson" => ".txt")
            @test isfile(txt_path)
            content = read(txt_path, String)
            @test occursin("KrylovSpectrumResults", content)
        end
    end

    # -----------------------------------------------------------------------
    # KrylovSpectrumResults with channel eigenvalues round-trip
    # -----------------------------------------------------------------------
    @testset "KrylovSpectrumResults channel path round-trip" begin
        mktempdir() do tmpdir
            config = make_small_thermalize_config(EnergyDomain())
            dim = SMALL_DIM
            eigenvalues = [0.0+0.0im, -0.3+0.0im]
            channel_eigs = [1.0+0.0im, 0.997+0.0im]

            result = KrylovSpectrumResults{Float64}(
                config, eigenvalues, 0.3,
                randn(ComplexF64, dim, dim), randn(ComplexF64, dim, dim),
                2, 50, 1, [1e-11, 1e-10],
                channel_eigs, 0.01, Dict{Symbol, Any}(),
            )
            path = joinpath(tmpdir, "test_krylov_channel.bson")
            save_result(result, path)

            loaded = load_result(path)
            @test loaded isa KrylovSpectrumResults
            @test loaded.channel_eigenvalues !== nothing
            @test isapprox(loaded.channel_eigenvalues, channel_eigs; atol=0)
            @test loaded.delta_used == 0.01
        end
    end

    # -----------------------------------------------------------------------
    # TrajectoryResults round-trip (plain mode)
    # -----------------------------------------------------------------------
    @testset "TrajectoryResults plain round-trip" begin
        mktempdir() do tmpdir
            config = make_small_thermalize_config(EnergyDomain())
            dim = SMALL_DIM
            rho_mean = Matrix(random_density_matrix(Int(log2(dim))))
            metadata = Dict{Symbol, Any}(:wall_time_seconds => 3.0, :n_threads => 1,
                :timestamp => "2026-02-27", :git_hash => "jkl012")

            result = TrajectoryResults{Float64}(config, rho_mean, 100, 42,
                nothing, nothing, nothing, metadata)
            path = joinpath(tmpdir, "test_traj_plain.bson")
            save_result(result, path)

            loaded = load_result(path)
            @test loaded isa TrajectoryResults
            @test isapprox(loaded.rho_mean, rho_mean; atol=0)
            @test loaded.n_trajectories == 100
            @test loaded.seed == 42
            @test loaded.times === nothing
            @test loaded.measurements_mean === nothing
            @test loaded.convergence === nothing
            @test loaded.metadata[:wall_time_seconds] == 3.0

            txt_path = replace(path, ".bson" => ".txt")
            @test isfile(txt_path)
        end
    end

    # -----------------------------------------------------------------------
    # TrajectoryResults with convergence data round-trip
    # -----------------------------------------------------------------------
    @testset "TrajectoryResults convergence round-trip" begin
        mktempdir() do tmpdir
            config = make_small_thermalize_config(EnergyDomain())
            dim = SMALL_DIM
            rho_mean = Matrix(random_density_matrix(Int(log2(dim))))
            conv = ConvergenceData(
                [100, 100], [100, 200], [0.3, 0.1],
                ["Z1", "Z2"], [0.5 0.4; 0.6 0.3], [0.55, 0.45],
            )
            metadata = Dict{Symbol, Any}()

            result = TrajectoryResults{Float64}(config, rho_mean, 200, 99,
                nothing, nothing, conv, metadata)
            path = joinpath(tmpdir, "test_traj_conv.bson")
            save_result(result, path)

            loaded = load_result(path)
            @test loaded isa TrajectoryResults
            @test loaded.convergence !== nothing
            @test loaded.convergence isa ConvergenceData
            @test loaded.convergence.trace_distances == [0.3, 0.1]
            @test loaded.convergence.observable_names == ["Z1", "Z2"]
        end
    end

    # -----------------------------------------------------------------------
    # Metadata auto-capture excludes Julia version
    # -----------------------------------------------------------------------
    @testset "Metadata excludes Julia version" begin
        meta = QuantumFurnace._capture_metadata(n_threads=2, wall_time_seconds=1.0)
        @test !haskey(meta, :julia_version)
        @test haskey(meta, :timestamp)
        @test haskey(meta, :git_hash)
        @test haskey(meta, :n_threads)
        @test haskey(meta, :wall_time_seconds)
    end

end  # @testset "New Result types serialization"
