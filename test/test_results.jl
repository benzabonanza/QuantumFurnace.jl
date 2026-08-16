using Random
using BSON

# Every persisted field is compared exactly: a BSON round-trip must not alter
# numerical results or the physical configuration used to produce them.

function _test_config_roundtrip(loaded::Config, expected::Config)
    @test typeof(loaded) === typeof(expected)
    for field in fieldnames(typeof(expected))
        @test isequal(getfield(loaded, field), getfield(expected, field))
    end
end

@testset "New Result types serialization" begin

    # -----------------------------------------------------------------------
    # LindbladResults round-trip
    # -----------------------------------------------------------------------
    @testset "LindbladResults round-trip" begin
        mktempdir() do tmpdir
            beta = 5.0
            r_D = 8
            w0_D = 0.1
            config = Config(;
                sim=Lindbladian(), domain=TimeDomain(), construction=DLL(),
                num_qubits=3, with_linear_combination=true,
                beta=beta, beta_phys=2.0, sigma=1 / beta, a=0.0, s=0.25,
                num_energy_bits_D=r_D, w0_D=w0_D,
                t0_D=2π / (2^r_D * w0_D), num_trotter_steps_per_t0_D=7,
                filter=DLLGaussianFilter(beta),
            )
            dim = N3_DIM
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
            _test_config_roundtrip(loaded.config, config)
            @test loaded.eigenvalues == eigenvalues
            @test loaded.fixed_point == fixed_point
            @test loaded.gap_mode == gap_mode
            @test loaded.spectral_gap == spectral_gap
            @test loaded.metadata == metadata

            # Companion .txt exists
            txt_path = replace(path, ".bson" => ".txt")
            @test isfile(txt_path)
            content = read(txt_path, String)
            @test occursin("LindbladResults", content)
            @test occursin("TimeDomain", content)
        end
    end

    # -----------------------------------------------------------------------
    # ThermalizeResults round-trip
    # -----------------------------------------------------------------------
    @testset "ThermalizeResults round-trip" begin
        mktempdir() do tmpdir
            beta = 5.0
            config = Config(;
                sim=Thermalize(), domain=TimeDomain(), construction=KMS(),
                num_qubits=3, with_linear_combination=true,
                beta=beta, beta_phys=2.0, sigma=1 / beta, a=0.0, s=0.25,
                num_energy_bits_D=8, w0_D=0.1, t0_D=2π / (2^8 * 0.1),
                num_energy_bits_b_minus=7, w0_b_minus=0.2,
                t0_b_minus=2π / (2^7 * 0.2),
                num_energy_bits_b_plus=6, w0_b_plus=0.25,
                t0_b_plus=2π / (2^6 * 0.25),
                num_trotter_steps_per_t0=9,
                num_trotter_steps_per_t0_D=7,
                num_trotter_steps_per_t0_b_minus=5,
                num_trotter_steps_per_t0_b_plus=3,
                mixing_time=0.4, delta=0.01,
                with_gqsp=true, gqsp_degree=2, jump_selection=:random,
            )
            dim = N3_DIM
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
            _test_config_roundtrip(loaded.config, config)
            @test loaded.final_dm == final_dm
            @test loaded.trace_distances == trace_distances
            @test loaded.time_steps == time_steps
            @test loaded.metadata == metadata

            txt_path = replace(path, ".bson" => ".txt")
            @test isfile(txt_path)
        end
    end

    # -----------------------------------------------------------------------
    # KrylovSpectrumResults round-trip
    # -----------------------------------------------------------------------
    @testset "KrylovSpectrumResults round-trip" begin
        mktempdir() do tmpdir
            config = Config(;
                sim=KrylovSpectrum(), domain=EnergyDomain(), construction=GNS(),
                num_qubits=3, with_linear_combination=true,
                beta=10.0, sigma=0.1, a=0.0, s=0.25,
                num_energy_bits=8, w0=0.05,
            )
            dim = N3_DIM
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
            _test_config_roundtrip(loaded.config, config)
            @test loaded.eigenvalues == eigenvalues
            @test loaded.spectral_gap == spectral_gap
            @test loaded.fixed_point == fixed_point
            @test loaded.gap_mode == gap_mode
            @test loaded.converged == 4
            @test loaded.matvec_count == 100
            @test loaded.num_restarts == 2
            @test loaded.normres == [1e-11, 1e-10, 1e-9, 1e-8]
            @test loaded.channel_eigenvalues === nothing
            @test loaded.delta_used === nothing
            @test loaded.metadata == metadata

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
            config = make_config(Thermalize(), EnergyDomain(); num_qubits=3, construction=GNS())
            dim = N3_DIM
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
            @test loaded.channel_eigenvalues == channel_eigs
            @test loaded.delta_used == 0.01
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

    @testset "Serialization path and tag validation" begin
        config = Config(;
            sim=Lindbladian(), domain=EnergyDomain(), construction=GNS(),
            num_qubits=3, with_linear_combination=false,
            beta=1.25, sigma=0.8,
        )
        result = LindbladResults{Float64}(
            config,
            ComplexF64[0.0, -0.5],
            Matrix{ComplexF64}(I, N3_DIM, N3_DIM) / N3_DIM,
            zeros(ComplexF64, N3_DIM, N3_DIM),
            -0.5 + 0.0im,
            Dict{Symbol, Any}(),
        )

        mktempdir() do tmpdir
            requested_path = joinpath(tmpdir, "result.data")
            write(requested_path, "unrelated payload")

            saved_path = save_result(result, requested_path)

            @test saved_path == joinpath(tmpdir, "result.bson")
            @test isfile(saved_path)
            @test isfile(joinpath(tmpdir, "result.txt"))
            @test read(requested_path, String) == "unrelated payload"
            @test load_result(saved_path) isa LindbladResults
        end

        serialized = QuantumFurnace._config_to_dict(config)
        legacy = copy(serialized)
        legacy[:config_kind] = "liouv"
        @test QuantumFurnace._reconstruct_config(legacy).sim isa Lindbladian

        bad_construction = copy(serialized)
        bad_construction[:config_type] = "KMZ"
        @test_throws ArgumentError QuantumFurnace._reconstruct_config(bad_construction)

        bad_simulation = copy(serialized)
        bad_simulation[:config_kind] = "thermalise"
        @test_throws ArgumentError QuantumFurnace._reconstruct_config(bad_simulation)
    end

    @testset "Generated filenames preserve beta" begin
        function filename_result(beta)
            config = Config(;
                sim=Lindbladian(), domain=EnergyDomain(), construction=GNS(),
                num_qubits=3, with_linear_combination=false,
                beta=beta, sigma=inv(beta),
            )
            return LindbladResults{Float64}(
                config,
                ComplexF64[0.0, -0.5],
                Matrix{ComplexF64}(I, N3_DIM, N3_DIM) / N3_DIM,
                zeros(ComplexF64, N3_DIM, N3_DIM),
                -0.5 + 0.0im,
                Dict{Symbol, Any}(),
            )
        end

        name_1 = QuantumFurnace._generate_result_filename(filename_result(1.25))
        name_2 = QuantumFurnace._generate_result_filename(filename_result(1.4))
        @test name_1 != name_2
        @test occursin("beta1.25", name_1)
        @test occursin("beta1.4", name_2)
    end

end  # @testset "New Result types serialization"
