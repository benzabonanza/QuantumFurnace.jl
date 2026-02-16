# ============================================================================
# ExperimentResult: struct + Dict-based BSON serialization
# ============================================================================

"""
    ExperimentResult{C<:AbstractConfig, T<:AbstractFloat}

Complete experiment result with config, trajectory data, and metadata.
Parameterized on config type `C` (preserves KMS/GNS distinction) and
element type `T` (typically Float64).

# Fields
- `config::C`: Full configuration used (everything needed to re-run).
- `trajectory_result::TrajectoryResult{Complex{T}}`: Trajectory-averaged density matrix + stats.
- `hamiltonian_params::Dict{Symbol, Any}`: Hamiltonian provenance (base_terms, base_coeffs, etc.).
- `metadata::Dict{Symbol, Any}`: Run metadata (timestamp, git hash, Julia version, wall time, etc.).
"""
struct ExperimentResult{C<:AbstractConfig, T<:AbstractFloat}
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
    _config_to_dict(config::AbstractConfig) -> Dict{Symbol, Any}

Serialize a config struct to a Dict with string-tagged type info.
"""
function _config_to_dict(config::AbstractConfig)
    d = Dict{Symbol, Any}()

    # Type tags
    d[:config_type] = (config isa Union{LiouvConfigGNS, ThermalizeConfigGNS}) ? "GNS" : "KMS"
    d[:config_kind] = (config isa AbstractThermalizeConfig) ? "thermalize" : "liouv"
    d[:domain] = string(typeof(config.domain))

    # Shared fields (all config types have these)
    d[:num_qubits]              = config.num_qubits
    d[:with_coherent]           = config.with_coherent
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
    if config isa AbstractThermalizeConfig
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
    return Dict{Symbol, Any}(
        :rho_mean          => Matrix(traj.rho_mean),
        :n_trajectories    => traj.n_trajectories,
        :seed              => traj.seed,
        :times             => traj.times,
        :measurements_mean => traj.measurements_mean,
    )
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
    traj = TrajectoryResult(
        traj_d[:rho_mean],
        traj_d[:n_trajectories],
        traj_d[:seed],
        get(traj_d, :times, nothing),
        get(traj_d, :measurements_mean, nothing),
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
    _reconstruct_config(d::Dict) -> AbstractConfig

Reconstruct the correct config struct from a serialized Dict.
Uses config_type ("KMS"/"GNS") and config_kind ("liouv"/"thermalize") to pick the constructor.
"""
function _reconstruct_config(d::Dict)
    domain = _string_to_domain(d[:domain])
    config_type = d[:config_type]   # "KMS" or "GNS"

    # Determine liouv vs thermalize: prefer config_kind tag, fall back to presence of mixing_time
    config_kind = get(d, :config_kind, nothing)
    is_thermalize = if config_kind !== nothing
        config_kind == "thermalize"
    else
        haskey(d, :mixing_time) && d[:mixing_time] !== nothing
    end

    kwargs = _dict_to_config_kwargs(d, domain)

    if config_type == "GNS" && is_thermalize
        ThermalizeConfigGNS(; kwargs...)
    elseif config_type == "GNS"
        LiouvConfigGNS(; kwargs...)
    elseif is_thermalize
        ThermalizeConfig(; kwargs...)
    else
        LiouvConfig(; kwargs...)
    end
end

"""
    _dict_to_config_kwargs(d::Dict, domain) -> Dict{Symbol, Any}

Build a kwargs Dict from serialized config fields, suitable for @kwdef constructors.
Filters out nothing values for optional fields to let defaults apply.
"""
function _dict_to_config_kwargs(d::Dict, domain)
    kwargs = Dict{Symbol, Any}()

    # Required fields
    kwargs[:num_qubits]              = d[:num_qubits]
    kwargs[:with_coherent]           = d[:with_coherent]
    kwargs[:with_linear_combination] = d[:with_linear_combination]
    kwargs[:domain]                  = domain
    kwargs[:beta]                    = d[:beta]
    kwargs[:sigma]                   = d[:sigma]

    # Optional fields (only set if present and non-nothing)
    for key in (:gaussian_parameters, :a, :b, :num_energy_bits, :t0, :w0, :eta, :num_trotter_steps_per_t0)
        val = get(d, key, nothing)
        if val !== nothing
            kwargs[key] = val
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
