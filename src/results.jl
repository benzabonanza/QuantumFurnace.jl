# ============================================================================
# ExperimentResult: struct + Dict-based BSON serialization
# ============================================================================

"""
    ExperimentResult{C<:Config, T<:AbstractFloat}

Complete experiment result with config, trajectory data, and metadata.
Parameterized on config type `C` (preserves KMS/GNS distinction via type parameter) and
element type `T` (typically Float64).

# Fields
- `config::C`: Full configuration used (everything needed to re-run).
- `trajectory_result::TrajectoryResult{Complex{T}}`: Trajectory-averaged density matrix + stats.
- `hamiltonian_params::Dict{Symbol, Any}`: Hamiltonian provenance (base_terms, base_coeffs, etc.).
- `metadata::Dict{Symbol, Any}`: Run metadata (timestamp, git hash, Julia version, wall time, etc.).
"""
struct ExperimentResult{C<:Config, T<:AbstractFloat}
    config::C
    trajectory_result::TrajectoryResult{Complex{T}}
    hamiltonian_params::Dict{Symbol, Any}
    metadata::Dict{Symbol, Any}
end

# ---------------------------------------------------------------------------
# Domain string <-> singleton lookup
# ---------------------------------------------------------------------------

const DOMAIN_LOOKUP = Dict(
    "TrotterDomain" => TrotterDomain(),
    "TimeDomain"    => TimeDomain(),
    "EnergyDomain"  => EnergyDomain(),
    "BohrDomain"    => BohrDomain(),
)

_string_to_domain(s::AbstractString) = DOMAIN_LOOKUP[s]

# ---------------------------------------------------------------------------
# ExperimentResult <-> Dict conversion (for safe BSON serialization)
# ---------------------------------------------------------------------------

"""
    _experiment_to_dict(result::ExperimentResult) -> Dict{Symbol, Any}

Convert an ExperimentResult to a plain Dict for BSON serialization.
All fields are stored as primitive types (String, Float64, Int, Matrix, Vector, Nothing).
"""
function _experiment_to_dict(result::ExperimentResult)
    return Dict{Symbol, Any}(
        :config             => _config_to_dict(result.config),
        :trajectory         => _trajectory_to_dict(result.trajectory_result),
        :hamiltonian_params => result.hamiltonian_params,
        :metadata           => result.metadata,
    )
end

"""
    _config_to_dict(config::Config) -> Dict{Symbol, Any}

Serialize a config struct to a Dict with string-tagged type info.
"""
function _config_to_dict(config::Config)
    d = Dict{Symbol, Any}()

    # Type tags
    d[:config_type] = config.construction isa GNS ? "GNS" : config.construction isa KMS ? "KMS" : "DLL"
    d[:config_kind] = config.sim isa Thermalize ? "thermalize" : "lindbladian"
    d[:domain] = string(typeof(config.domain))

    # Shared fields (all config types have these)
    d[:num_qubits]              = config.num_qubits
    d[:with_coherent]           = with_coherent(config.construction)
    d[:with_linear_combination] = config.with_linear_combination
    d[:beta]                    = config.beta
    d[:sigma]                   = config.sigma
    d[:gaussian_parameters]     = config.gaussian_parameters
    d[:a]                       = config.a
    d[:b]                       = config.b
    d[:num_energy_bits]         = config.num_energy_bits
    d[:t0]                      = config.t0
    d[:w0]                      = config.w0
    d[:eta]                     = config.eta
    d[:num_trotter_steps_per_t0] = config.num_trotter_steps_per_t0

    # Thermalize-specific fields
    if config.sim isa Thermalize
        d[:mixing_time] = config.mixing_time
        d[:delta]       = config.delta
    end

    return d
end

"""
    _trajectory_to_dict(traj::TrajectoryResult) -> Dict{Symbol, Any}

Serialize TrajectoryResult to a Dict. Unwraps Hermitian to Matrix for BSON safety.
"""
function _trajectory_to_dict(traj::TrajectoryResult)
    d = Dict{Symbol, Any}(
        :rho_mean          => Matrix(traj.rho_mean),
        :n_trajectories    => traj.n_trajectories,
        :seed              => traj.seed,
        :times             => traj.times,
        :measurements_mean => traj.measurements_mean,
    )
    if traj.convergence !== nothing
        d[:convergence] = _convergence_to_dict(traj.convergence)
    end
    return d
end

# ---------------------------------------------------------------------------
# Dict -> ExperimentResult reconstruction
# ---------------------------------------------------------------------------

"""
    _dict_to_experiment(d::Dict) -> ExperimentResult

Reconstruct an ExperimentResult from a Dict loaded from BSON.
Uses `get()` with defaults for forward compatibility (missing fields = nothing/empty).
"""
function _dict_to_experiment(d::Dict)
    # Reconstruct config
    config = _reconstruct_config(d[:config])

    # Reconstruct TrajectoryResult
    traj_d = d[:trajectory]
    conv = if haskey(traj_d, :convergence) && traj_d[:convergence] !== nothing
        _dict_to_convergence(traj_d[:convergence])
    else
        nothing
    end
    traj = TrajectoryResult(
        traj_d[:rho_mean],
        traj_d[:n_trajectories],
        traj_d[:seed],
        get(traj_d, :times, nothing),
        get(traj_d, :measurements_mean, nothing),
        conv,
    )

    # Forward-compatible metadata and hamiltonian_params
    metadata    = get(d, :metadata, Dict{Symbol, Any}())
    ham_params  = get(d, :hamiltonian_params, Dict{Symbol, Any}())

    # Infer type parameters
    T = eltype(real(traj.rho_mean))
    C = typeof(config)

    return ExperimentResult{C, T}(config, traj, ham_params, metadata)
end

"""
    _reconstruct_config(d::Dict) -> Config

Reconstruct the correct Config struct from a serialized Dict.
Uses config_type ("KMS"/"GNS"/"DLL") and config_kind ("liouv"/"thermalize") to pick the singletons.
"""
function _reconstruct_config(d::Dict)
    domain = _string_to_domain(d[:domain])
    config_type = d[:config_type]   # "KMS", "GNS", or "DLL"

    # Determine liouv vs thermalize: prefer config_kind tag, fall back to presence of mixing_time
    # "liouv" (old) and "lindbladian" (new) both map to Lindbladian() (the else branch)
    config_kind = get(d, :config_kind, nothing)
    is_thermalize = if config_kind !== nothing
        config_kind == "thermalize"
    else
        haskey(d, :mixing_time) && d[:mixing_time] !== nothing
    end

    kwargs = _dict_to_config_kwargs(d, domain)

    construction = config_type == "GNS" ? GNS() : config_type == "DLL" ? DLL() : KMS()
    sim = is_thermalize ? Thermalize() : Lindbladian()

    Config(; sim=sim, domain=domain, construction=construction, kwargs...)
end

"""
    _dict_to_config_kwargs(d::Dict, domain) -> Dict{Symbol, Any}

Build a kwargs Dict from serialized config fields, suitable for @kwdef constructors.
Filters out nothing values for optional fields to let defaults apply.
"""
function _dict_to_config_kwargs(d::Dict, domain)
    kwargs = Dict{Symbol, Any}()

    # Required fields (with_coherent derived from construction type, not stored)
    kwargs[:num_qubits]              = d[:num_qubits]
    kwargs[:with_linear_combination] = d[:with_linear_combination]
    kwargs[:beta]                    = d[:beta]
    kwargs[:sigma]                   = d[:sigma]

    # Optional fields (only set if present and non-nothing)
    for key in (:a, :b, :num_energy_bits, :t0, :w0, :eta, :num_trotter_steps_per_t0)
        val = get(d, key, nothing)
        if val !== nothing
            kwargs[key] = val
        end
    end

    # gaussian_parameters: BSON stores tuples as arrays, so convert back
    gp = get(d, :gaussian_parameters, nothing)
    if gp !== nothing
        if gp isa AbstractVector
            kwargs[:gaussian_parameters] = (gp[1], gp[2])
        else
            kwargs[:gaussian_parameters] = gp
        end
    end

    # Thermalize-specific
    for key in (:mixing_time, :delta)
        val = get(d, key, nothing)
        if val !== nothing
            kwargs[key] = val
        end
    end

    return kwargs
end

# ---------------------------------------------------------------------------
# ConvergenceData <-> Dict conversion (for safe BSON serialization)
# ---------------------------------------------------------------------------

"""
    _convergence_to_dict(conv::ConvergenceData) -> Dict{Symbol, Any}

Convert a ConvergenceData to a plain Dict for BSON serialization.
"""
function _convergence_to_dict(conv::ConvergenceData)
    return Dict{Symbol, Any}(
        :batch_sizes              => conv.batch_sizes,
        :cumulative_n_traj        => conv.cumulative_n_traj,
        :trace_distances          => conv.trace_distances,
        :observable_names         => conv.observable_names,
        :observable_values        => conv.observable_values,
        :observable_gibbs_values  => conv.observable_gibbs_values,
        # Phase 17: adaptive diagnostics
        :converged                => conv.converged,
        :final_relative_change    => conv.final_relative_change,
        :consecutive_stable_batches => conv.consecutive_stable_batches,
        :total_batches            => conv.total_batches,
    )
end

"""
    _dict_to_convergence(d::Dict) -> ConvergenceData

Reconstruct a ConvergenceData from a Dict loaded from BSON.
Uses `get` with default for forward compatibility.
"""
function _dict_to_convergence(d::Dict)
    return ConvergenceData(
        d[:batch_sizes],
        d[:cumulative_n_traj],
        d[:trace_distances],
        d[:observable_names],
        d[:observable_values],
        get(d, :observable_gibbs_values, Float64[]),
        # Phase 17: adaptive diagnostics (backward-compatible defaults for pre-Phase 17 data)
        get(d, :converged, false),
        get(d, :final_relative_change, NaN),
        get(d, :consecutive_stable_batches, 0),
        get(d, :total_batches, length(d[:batch_sizes])),
    )
end

# ============================================================================
# Metadata auto-capture
# ============================================================================

"""
    _capture_metadata(; n_threads, wall_time_seconds, extra) -> Dict{Symbol, Any}

Auto-capture run metadata: timestamp, git hash, thread count, wall time.
Merges any extra key-value pairs from the `extra` Dict.
"""
function _capture_metadata(;
    n_threads::Int = Threads.nthreads(),
    wall_time_seconds::Union{Float64, Nothing} = nothing,
    extra::Dict{Symbol, Any} = Dict{Symbol, Any}(),
)
    meta = Dict{Symbol, Any}(
        :timestamp         => Dates.format(Dates.now(), dateformat"yyyy-mm-dd_HH:MM:SS"),
        :git_hash          => _capture_git_hash(),
        :n_threads         => n_threads,
        :wall_time_seconds => wall_time_seconds,
    )
    merge!(meta, extra)
    return meta
end

"""
    _capture_git_hash() -> String

Capture the current HEAD commit hash via LibGit2. Returns "unknown" on failure.
"""
function _capture_git_hash()
    try
        project_root = dirname(Pkg.project().path)
        repo = LibGit2.GitRepo(project_root)
        hash = string(LibGit2.head_oid(repo))
        close(repo)
        return hash
    catch
        return "unknown"
    end
end

# ============================================================================
# Hamiltonian parameter extraction
# ============================================================================

"""
    _extract_hamiltonian_params(ham::HamHam) -> Dict{Symbol, Any}

Extract the minimal set of Hamiltonian parameters needed for provenance/reconstruction.
Does NOT store derived quantities (eigendecomposition, bohr_freqs, bohr_dict, gibbs).
"""
function _extract_hamiltonian_params(ham::HamHam)
    return Dict{Symbol, Any}(
        :num_qubits        => Int(log2(size(ham.data, 1))),
        :base_coeffs       => ham.base_coeffs,
        :base_terms        => [Matrix.(term_group) for term_group in ham.base_terms],
        :disordering_term  => ham.disordering_term === nothing ? nothing : Matrix.(ham.disordering_term),
        :disordering_coeffs => ham.disordering_coeffs,
        :periodic          => ham.periodic,
        :shift             => ham.shift,
        :rescaling_factor  => ham.rescaling_factor,
    )
end

# ============================================================================
# Save / Load wrappers
# ============================================================================

"""
    save_experiment(result::ExperimentResult, path::String) -> String

Save an ExperimentResult to a BSON file at `path`, plus a companion `.txt` file.
Creates parent directories as needed. Returns the path.
"""
function save_experiment(result::ExperimentResult, path::String)
    d = _experiment_to_dict(result)
    mkpath(dirname(path))
    BSON.bson(path, d)
    _write_companion_txt(result, replace(path, ".bson" => ".txt"))
    return path
end

"""
    save_experiment(result::ExperimentResult) -> String

Save an ExperimentResult to the default results directory with an auto-generated filename.
"""
function save_experiment(result::ExperimentResult)
    dir = _default_results_dir(result.config)
    filename = _generate_experiment_filename(result.config)
    path = joinpath(dir, filename)
    return save_experiment(result, path)
end

"""
    load_experiment(path::String) -> ExperimentResult

Load an ExperimentResult from a BSON file.
"""
function load_experiment(path::String)
    d = BSON.load(path)
    return _dict_to_experiment(d)
end

# ============================================================================
# Companion text file
# ============================================================================

"""
    _write_companion_txt(result::ExperimentResult, path::String)

Write a human-readable summary alongside the BSON file.
Simple key-value format for quick browsing without Julia.
"""
function _write_companion_txt(result::ExperimentResult, path::String)
    open(path, "w") do io
        cfg  = result.config
        meta = result.metadata

        println(io, "=== QuantumFurnace Experiment Result ===")
        println(io)
        println(io, "Date:       ", get(meta, :timestamp, "unknown"))
        println(io, "Git:        ", get(meta, :git_hash, "unknown"))
        println(io, "Julia:      ", get(meta, :julia_version, "unknown"))
        println(io)
        println(io, "--- Config ---")
        println(io, "Type:       ", cfg.construction isa GNS ? "GNS" : cfg.construction isa KMS ? "KMS" : "DLL")
        println(io, "Kind:       ", cfg.sim isa Thermalize ? "thermalize" : "liouv")
        println(io, "Domain:     ", typeof(cfg.domain))
        println(io, "n_qubits:   ", cfg.num_qubits)
        println(io, "beta:       ", cfg.beta)
        println(io, "sigma:      ", cfg.sigma)
        if cfg.sim isa Thermalize
            println(io, "mix_time:   ", cfg.mixing_time)
            println(io, "delta:      ", cfg.delta)
        end
        println(io)
        println(io, "--- Results ---")
        traj = result.trajectory_result
        println(io, "N_traj:     ", traj.n_trajectories)
        println(io, "Seed:       ", traj.seed)
        println(io, "Threads:    ", get(meta, :n_threads, "unknown"))
        println(io, "Wall time:  ", get(meta, :wall_time_seconds, "unknown"), " s")
        println(io, "rho dim:    ", size(traj.rho_mean, 1), "x", size(traj.rho_mean, 2))
    end
end

# ============================================================================
# Filename generation and default paths
# ============================================================================

"""
    _generate_experiment_filename(config::Config) -> String

Generate a descriptive filename: `{db}_{n}_{beta}_{domain}_{date}.bson`.
"""
function _generate_experiment_filename(config::Config)
    db_str     = config.construction isa GNS ? "gns" : "kms"
    domain_str = lowercase(replace(string(typeof(config.domain)), "Domain" => ""))
    n_str      = "n$(config.num_qubits)"
    beta_str   = "beta$(round(Int, config.beta))"
    date_str   = Dates.format(Dates.now(), dateformat"yyyymmdd")
    return "$(db_str)_$(n_str)_$(beta_str)_$(domain_str)_$(date_str).bson"
end

"""
    _default_results_dir(config::Config) -> String

Return the default results subdirectory for the given config type.
"""
function _default_results_dir(config::Config)
    subdir = config.construction isa GNS ? "approx_gns" : "kms"
    return joinpath(dirname(Pkg.project().path), "results", subdir)
end

# ============================================================================
# New Result serialization (Phase 36)
# ============================================================================

# ---------------------------------------------------------------------------
# Result type tag helper
# ---------------------------------------------------------------------------

_result_type_tag(::LindbladResults) = "lindblad"
_result_type_tag(::ThermalizeResults) = "thermalize"
_result_type_tag(::KrylovSpectrumResults) = "krylov_spectrum"
_result_type_tag(::TrajectoryResults) = "trajectory"

# ---------------------------------------------------------------------------
# Result -> Dict conversion (for BSON serialization)
# ---------------------------------------------------------------------------

function _lindblad_to_dict(r::LindbladResults)
    return Dict{Symbol, Any}(
        :result_type  => "lindblad",
        :config       => _config_to_dict(r.config),
        :eigenvalues  => r.eigenvalues,
        :fixed_point  => Matrix(r.fixed_point),
        :gap_mode     => Matrix(r.gap_mode),
        :spectral_gap => r.spectral_gap,
        :metadata     => r.metadata,
    )
end

function _thermalize_to_dict(r::ThermalizeResults)
    return Dict{Symbol, Any}(
        :result_type      => "thermalize",
        :config           => _config_to_dict(r.config),
        :final_dm         => Matrix(r.final_dm),
        :trace_distances  => r.trace_distances,
        :time_steps       => r.time_steps,
        :metadata         => r.metadata,
    )
end

function _krylov_spectrum_to_dict(r::KrylovSpectrumResults)
    return Dict{Symbol, Any}(
        :result_type          => "krylov_spectrum",
        :config               => _config_to_dict(r.config),
        :eigenvalues          => r.eigenvalues,
        :spectral_gap         => r.spectral_gap,
        :fixed_point          => Matrix(r.fixed_point),
        :gap_mode             => Matrix(r.gap_mode),
        :converged            => r.converged,
        :matvec_count         => r.matvec_count,
        :num_restarts         => r.num_restarts,
        :normres              => r.normres,
        :channel_eigenvalues  => r.channel_eigenvalues,
        :delta_used           => r.delta_used,
        :metadata             => r.metadata,
    )
end

function _trajectory_to_dict_new(r::TrajectoryResults)
    d = Dict{Symbol, Any}(
        :result_type        => "trajectory",
        :config             => _config_to_dict(r.config),
        :rho_mean           => Matrix(r.rho_mean),
        :n_trajectories     => r.n_trajectories,
        :seed               => r.seed,
        :times              => r.times,
        :measurements_mean  => r.measurements_mean,
        :metadata           => r.metadata,
    )
    if r.convergence !== nothing
        d[:convergence] = _convergence_to_dict(r.convergence)
    end
    return d
end

"""
    _result_to_dict(r::AbstractResults) -> Dict{Symbol, Any}

Convert any typed Result to a Dict for BSON serialization.
Dispatches to the appropriate type-specific conversion function.
"""
function _result_to_dict(r::LindbladResults)
    return _lindblad_to_dict(r)
end
function _result_to_dict(r::ThermalizeResults)
    return _thermalize_to_dict(r)
end
function _result_to_dict(r::KrylovSpectrumResults)
    return _krylov_spectrum_to_dict(r)
end
function _result_to_dict(r::TrajectoryResults)
    return _trajectory_to_dict_new(r)
end

# ---------------------------------------------------------------------------
# Dict -> Result reconstruction
# ---------------------------------------------------------------------------

function _dict_to_lindblad_results(d::Dict)
    config = _reconstruct_config(d[:config])
    T = real(eltype(d[:eigenvalues]))
    return LindbladResults{T}(
        config,
        d[:eigenvalues],
        d[:fixed_point],
        d[:gap_mode],
        d[:spectral_gap],
        d[:metadata],
    )
end

function _dict_to_thermalize_results(d::Dict)
    config = _reconstruct_config(d[:config])
    T = eltype(d[:trace_distances])
    return ThermalizeResults{T}(
        config,
        d[:final_dm],
        d[:trace_distances],
        d[:time_steps],
        d[:metadata],
    )
end

function _dict_to_krylov_spectrum_results(d::Dict)
    config = _reconstruct_config(d[:config])
    T = eltype(d[:normres])
    return KrylovSpectrumResults{T}(
        config,
        d[:eigenvalues],
        d[:spectral_gap],
        d[:fixed_point],
        d[:gap_mode],
        d[:converged],
        d[:matvec_count],
        d[:num_restarts],
        d[:normres],
        get(d, :channel_eigenvalues, nothing),
        get(d, :delta_used, nothing),
        d[:metadata],
    )
end

function _dict_to_trajectory_results(d::Dict)
    config = _reconstruct_config(d[:config])
    conv = if haskey(d, :convergence) && d[:convergence] !== nothing
        _dict_to_convergence(d[:convergence])
    else
        nothing
    end
    T = real(eltype(d[:rho_mean]))
    return TrajectoryResults{T}(
        config,
        d[:rho_mean],
        d[:n_trajectories],
        d[:seed],
        get(d, :times, nothing),
        get(d, :measurements_mean, nothing),
        conv,
        d[:metadata],
    )
end

# ---------------------------------------------------------------------------
# save_result / load_result
# ---------------------------------------------------------------------------

"""
    save_result(result::AbstractResults, path::String) -> String

Save a typed Result to a BSON file at `path`, plus a companion `.txt` file.
Creates parent directories as needed. Returns the path.
"""
function save_result(result::AbstractResults, path::String)
    d = _result_to_dict(result)
    mkpath(dirname(path))
    BSON.bson(path, d)
    _write_result_companion_txt(result, replace(path, ".bson" => ".txt"))
    return path
end

"""
    load_result(path::String) -> AbstractResults

Load a typed Result from a BSON file. Auto-detects the result type via the
`:result_type` tag stored in the BSON Dict.
"""
function load_result(path::String)
    d = BSON.load(path)
    tag = d[:result_type]
    if tag == "lindblad"
        return _dict_to_lindblad_results(d)
    elseif tag == "thermalize"
        return _dict_to_thermalize_results(d)
    elseif tag == "krylov_spectrum"
        return _dict_to_krylov_spectrum_results(d)
    elseif tag == "trajectory"
        return _dict_to_trajectory_results(d)
    else
        error("Unknown result type: $tag")
    end
end

# ---------------------------------------------------------------------------
# Companion .txt file (per result type)
# ---------------------------------------------------------------------------

"""
    _write_result_companion_txt(result::AbstractResults, path::String)

Write a human-readable summary alongside the BSON file.
Format varies by result type, showing key metrics for quick browsing.
"""
function _write_result_companion_txt(result::AbstractResults, path::String)
    open(path, "w") do io
        cfg  = result.config
        meta = result.metadata
        type_name = typeof(result).name.name  # e.g. :LindbladResults

        println(io, "=== QuantumFurnace [$type_name] ===")
        println(io)
        println(io, "Date:       ", get(meta, :timestamp, "unknown"))
        println(io, "Git:        ", get(meta, :git_hash, "unknown"))
        println(io, "Threads:    ", get(meta, :n_threads, "unknown"))
        wt = get(meta, :wall_time_seconds, nothing)
        println(io, "Wall time:  ", wt === nothing ? "unknown" : "$wt s")
        println(io)
        println(io, "--- Config ---")
        println(io, "Construction: ", cfg.construction isa GNS ? "GNS" : cfg.construction isa KMS ? "KMS" : "DLL")
        println(io, "Domain:     ", typeof(cfg.domain))
        println(io, "n_qubits:   ", cfg.num_qubits)
        println(io, "beta:       ", cfg.beta)
        println(io)
        println(io, "--- Results ---")

        if result isa LindbladResults
            println(io, "Spectral gap (real): ", real(result.spectral_gap))
            println(io, "Spectral gap (imag): ", imag(result.spectral_gap))
            println(io, "Fixed point dim:     ", size(result.fixed_point, 1), "x", size(result.fixed_point, 2))
            println(io, "N eigenvalues:       ", length(result.eigenvalues))

        elseif result isa ThermalizeResults
            td = result.trace_distances
            println(io, "Final trace dist:    ", isempty(td) ? "N/A" : last(td))
            println(io, "N time steps:        ", length(result.time_steps))
            println(io, "Final DM dim:        ", size(result.final_dm, 1), "x", size(result.final_dm, 2))

        elseif result isa KrylovSpectrumResults
            println(io, "Spectral gap:        ", result.spectral_gap)
            println(io, "Matvec count:        ", result.matvec_count)
            println(io, "Converged:           ", result.converged)
            println(io, "Num restarts:        ", result.num_restarts)
            if result.delta_used !== nothing
                println(io, "Delta used:          ", result.delta_used)
            end
            if result.channel_eigenvalues !== nothing
                println(io, "Channel eigenvalues: ", length(result.channel_eigenvalues))
            end

        elseif result isa TrajectoryResults
            println(io, "N trajectories:      ", result.n_trajectories)
            println(io, "Seed:                ", result.seed)
            println(io, "rho_mean dim:        ", size(result.rho_mean, 1), "x", size(result.rho_mean, 2))
            if result.convergence !== nothing
                println(io, "Converged:           ", result.convergence.converged)
                println(io, "Total batches:       ", result.convergence.total_batches)
            end
        end
    end
end

# ---------------------------------------------------------------------------
# Auto-filename generation for new Results
# ---------------------------------------------------------------------------

"""
    _generate_result_filename(result::AbstractResults) -> String

Generate a descriptive filename: `{type}_{construction}_{n}_{beta}_{domain}_{date}.bson`.
"""
function _generate_result_filename(result::AbstractResults)
    cfg = result.config
    type_str = _result_type_tag(result)
    db_str = cfg.construction isa GNS ? "gns" : cfg.construction isa KMS ? "kms" : "dll"
    domain_str = lowercase(replace(string(typeof(cfg.domain)), "Domain" => ""))
    n_str = "n$(cfg.num_qubits)"
    beta_str = "beta$(round(Int, cfg.beta))"
    date_str = Dates.format(Dates.now(), dateformat"yyyymmdd")
    return "$(type_str)_$(db_str)_$(n_str)_$(beta_str)_$(domain_str)_$(date_str).bson"
end
