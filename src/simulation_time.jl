# ============================================================================
# Hamiltonian Simulation Time Counting (Phases 44-47)
# ============================================================================
#
# Resource estimation for Chen's quantum Gibbs sampling algorithm.
# Computes total Hamiltonian simulation time broken down by component:
#   - OFT (Operator Fourier Transform): dominant per-step cost
#   - B (coherent correction): KMS-only correction term
#   - Per-step: 2×OFT + B (OFT appears twice: forward + backward)
#   - Total: n_steps × per_step

# ---------------------------------------------------------------------------
# Phase 44: SimulationTimeBudget struct and QPE grid utilities
# ---------------------------------------------------------------------------

"""
    SimulationTimeBudget

Immutable result container for Hamiltonian simulation time resource estimation.
Stores per-component time costs, QPE grid parameters, physical parameters,
and filter configuration — enough to reproduce a paper table row.

# Cost fields
- `oft_time`: OFT Hamiltonian simulation time per step
- `b_time`: B coherent correction time per step (0.0 for GNS)
- `per_step_time`: `2 × oft_time + b_time`
- `n_steps`: `ceil(T / delta)`
- `total_time`: `n_steps × per_step_time`

# Grid info
- `r`, `N`, `w0`, `t0`, `energy_range`

# Physical parameters
- `beta`, `sigma`, `delta`, `construction`, `n_qubits`, `rescaling_factor`
- `T`: target simulation time (e.g. mixing time)

# Filter info
- `filter_type`: `:gaussian`, `:metropolis`, `:smooth_metropolis`, or `:kinky_metropolis`
- `filter_params`: Dict with filter-specific parameters
"""
struct SimulationTimeBudget
    # Cost fields
    oft_time::Float64
    b_time::Float64
    per_step_time::Float64
    n_steps::Int
    total_time::Float64
    # Grid info
    r::Int
    N::Int
    w0::Float64
    t0::Float64
    energy_range::Tuple{Float64, Float64}
    # Physical parameters
    beta::Float64
    sigma::Float64
    delta::Float64
    construction::Symbol
    n_qubits::Int
    rescaling_factor::Float64
    T::Float64
    # Filter info
    filter_type::Symbol
    filter_params::Dict{Symbol, Float64}
end

"""
    _qpe_grid_info(r, w0) -> (; N, t0, energy_range)

Compute QPE grid parameters from resolution `r` and energy spacing `w0`.

- `N = 2^r`
- `t0 = 2π / (N × w0)` (Fourier relation)
- `energy_range = (-N/2 × w0, (N/2 - 1) × w0)`
"""
function _qpe_grid_info(r::Int, w0::Real)
    N = 2^r
    t0 = 2π / (N * w0)
    energy_range = (Float64(-N÷2 * w0), Float64((N÷2 - 1) * w0))
    return (; N, t0, energy_range)
end

function Base.show(io::IO, b::SimulationTimeBudget)
    print(io, "SimulationTimeBudget(r=$(b.r), total=$(b.total_time), n_steps=$(b.n_steps), $(b.construction))")
end

function Base.show(io::IO, ::MIME"text/plain", b::SimulationTimeBudget)
    println(io, "SimulationTimeBudget")
    println(io, "  Grid: r=$(b.r), N=$(b.N), w0=$(b.w0), t0=$(b.t0)")
    println(io, "  Energy range: [$(b.energy_range[1]), $(b.energy_range[2])]")
    println(io, "  Physics: β=$(b.beta), σ=$(b.sigma), δ=$(b.delta), n=$(b.n_qubits)")
    println(io, "  Construction: $(b.construction), rescaling=$(b.rescaling_factor)")
    fp = isempty(b.filter_params) ? "" : "($(join(["$k=$v" for (k,v) in b.filter_params], ", ")))"
    println(io, "  Filter: $(b.filter_type)$fp")
    println(io, "  T: $(b.T)")
    println(io, "  ─────────────────────────────────")
    println(io, "  OFT time:      $(b.oft_time)")
    println(io, "  B time:        $(b.b_time)")
    println(io, "  Per step:      $(b.per_step_time)  (2×$(b.oft_time) + $(b.b_time))")
    println(io, "  Steps:         $(b.n_steps)")
    print(io,   "  Total:         $(b.total_time)")
end

# ---------------------------------------------------------------------------
# Phase 45: OFT Hamiltonian simulation time
# ---------------------------------------------------------------------------

"""
    _oft_hamiltonian_time(r, w0, transition_weights) -> Float64

Compute the OFT Hamiltonian simulation time: sum of `|t_k| × γ(w_k)`
over the full 2^r QPE grid.

When all weights are 1.0, equals `t0 × N² / 4` (closed-form validation).
"""
function _oft_hamiltonian_time(r::Int, w0::Real, transition_weights::AbstractVector{<:Real})
    N = 2^r
    t0 = 2π / (N * w0)
    @assert length(transition_weights) == N "transition_weights length must be 2^r=$N, got $(length(transition_weights))"
    total = 0.0
    @inbounds for i in 1:N
        k = -N÷2 + i - 1
        total += abs(k * t0) * transition_weights[i]
    end
    return total
end

# ---------------------------------------------------------------------------
# Phase 46: B coherent term Hamiltonian simulation time
# ---------------------------------------------------------------------------

"""
    _b_hamiltonian_time(b_minus, b_plus, beta, sigma, t0) -> Float64

Compute the B coherent correction Hamiltonian simulation time as a
double sum over truncated b-dict entries.

Returns `0.0` when either dictionary is `nothing` or empty (GNS case).

For KMS, each (t,s) pair in the B operator involves:
- Inner: 3 time evolutions of total duration 4|sβ|
- Outer: 2 time evolutions of total duration 2|t/σ|

Total: `t0² × Σ_t Σ_s |b_minus[t]| × |b_plus[s]| × (4|sβ| + 2|t/σ|)`
"""
function _b_hamiltonian_time(
    b_minus::Union{Nothing, Dict},
    b_plus::Union{Nothing, Dict},
    beta::Real,
    sigma::Real,
    t0::Real,
)
    (b_minus === nothing || b_plus === nothing) && return 0.0
    (isempty(b_minus) || isempty(b_plus)) && return 0.0

    # Inner contribution (b_plus): Σ_s |b_plus[s]| × |s × β|
    inner_weighted = sum(abs(v) * abs(s * beta) for (s, v) in b_plus)
    inner_norm = sum(abs(v) for (_, v) in b_plus)

    # Outer contribution (b_minus): Σ_t |b_minus[t]| × |t / σ|
    outer_weighted = sum(abs(v) * abs(t / sigma) for (t, v) in b_minus)
    outer_norm = sum(abs(v) for (_, v) in b_minus)

    # Factored double sum: t0² × [4 × ‖b_minus‖₁ × inner_weighted + 2 × outer_weighted × ‖b_plus‖₁]
    return t0^2 * (4.0 * outer_norm * inner_weighted + 2.0 * outer_weighted * inner_norm)
end

# ---------------------------------------------------------------------------
# Phase 47: Public API
# ---------------------------------------------------------------------------

function _determine_filter_info(config::Config)
    if !config.with_linear_combination
        params = Dict{Symbol, Float64}()
        config.gaussian_parameters[1] !== nothing && (params[:w_gamma] = Float64(config.gaussian_parameters[1]))
        config.gaussian_parameters[2] !== nothing && (params[:sigma_gamma] = Float64(config.gaussian_parameters[2]))
        return :gaussian, params
    end
    a_val = something(config.a, 0.0)
    s_val = something(config.s, 0.0)
    if a_val > 0 && s_val > 0
        return :smooth_metropolis, Dict{Symbol, Float64}(:a => a_val, :s => s_val)
    elseif a_val > 0
        return :metropolis, Dict{Symbol, Float64}(:a => a_val)
    else
        return :kinky_metropolis, Dict{Symbol, Float64}()
    end
end

"""
    compute_simulation_time(config, ham, T) -> SimulationTimeBudget

Compute a complete Hamiltonian simulation time budget for Chen's quantum
Gibbs sampling algorithm.

# Arguments
- `config::Config{Thermalize, <:Union{TimeDomain, TrotterDomain}}`: algorithm configuration
- `ham::HamHam`: Hamiltonian (extracts `rescaling_factor`)
- `T::Real`: target simulation time (e.g. mixing time)

# Returns
[`SimulationTimeBudget`](@ref) with `per_step = 2×OFT + B`, `total = n_steps × per_step`.
"""
function compute_simulation_time(
    config::Config{Thermalize, D},
    ham::HamHam,
    T::Real,
) where {D <: Union{TimeDomain, TrotterDomain}}
    delta = config.delta
    r = config.num_energy_bits
    w0 = config.w0

    delta !== nothing && delta > 0 || throw(ArgumentError("config.delta must be set and positive"))
    r !== nothing && r > 0 || throw(ArgumentError("config.num_energy_bits must be set and positive"))
    w0 !== nothing || throw(ArgumentError("config.w0 must be set"))
    T > 0 || throw(ArgumentError("T must be positive"))

    # Grid
    grid = _qpe_grid_info(r, w0)
    N = grid.N
    t0 = grid.t0

    # Transition weights on full QPE grid
    energy_labels = _create_energy_labels(r, w0)
    transition_fn = pick_transition(config)
    # NaN can occur at extreme energies (Inf * 0 in erfc branch); clamp to 0
    transition_weights = Float64[let v = transition_fn(w); isnan(v) ? 0.0 : v end for w in energy_labels]

    # OFT time
    oft_time = _oft_hamiltonian_time(r, w0, transition_weights)

    # B coherent correction time
    construction = config.construction
    b_coh_time = if with_coherent(construction)
        time_labels = energy_labels .* (t0 / w0)
        bm = _compute_truncated_func(_compute_b_minus, time_labels, config.beta, config.sigma)
        bp_fn, bp_args = _select_b_plus_calculator(config)
        bp = _compute_truncated_func(bp_fn, time_labels, bp_args...)
        _b_hamiltonian_time(bm, bp, config.beta, config.sigma, t0)
    else
        0.0
    end

    # Assembly
    per_step_time = 2.0 * oft_time + b_coh_time
    n_steps = ceil(Int, T / delta)
    total_time = Float64(n_steps) * per_step_time

    # Metadata
    n_qubits = config.num_qubits
    construction_sym = construction isa KMS ? :KMS : construction isa GNS ? :GNS : :DLL
    filter_type, filter_params = _determine_filter_info(config)

    return SimulationTimeBudget(
        oft_time, b_coh_time, per_step_time, n_steps, total_time,
        r, N, w0, t0, grid.energy_range,
        config.beta, config.sigma, delta, construction_sym, n_qubits, ham.rescaling_factor, Float64(T),
        filter_type, filter_params,
    )
end
