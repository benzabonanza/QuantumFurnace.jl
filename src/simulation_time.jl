# Resource estimates for the OFT and coherent parts of the sampling channel.
# Math: $t_total = ceil(T / delta) (2 t_OFT + t_B)$. GQSP-labelled counts are
# formal until the current unscaled polynomial is made contractive and its
# complementary polynomial/circuit is synthesized.

"""
    SimulationTimeBudget

Hamiltonian-simulation time budget for a thermalisation run.

# Fields
- `oft_time`, `b_per_be`, `b_time`: per-step OFT and coherent costs. GQSP
  uses `2 * gqsp_degree` coherent block-encoding queries (Motlagh--Wiebe,
  Theorem 6, Form B).
- `per_step_time`, `n_steps`, `total_time`: combined cost, step count, and total.
- `r_D`/`r_bm`/`r_bp` groups: independent dissipative and coherent registers.
- `with_gqsp`, `gqsp_degree`: coherent implementation metadata.
- `cost_interpretation`: `:physical_channel` for the implemented direct
  coherent step, or `:formal_unscaled_gqsp_surrogate` when the query count
  prices the current unsynthesised polynomial surrogate.
- `beta`, `sigma`, `delta`, `construction`, `n_qubits`, `rescaling_factor`,
  `T`: physical and algorithmic metadata.
- `filter_type`, `filter_params`: transition-filter metadata.
"""
struct SimulationTimeBudget
    # Cost fields
    oft_time::Float64
    b_per_be::Float64
    b_time::Float64
    per_step_time::Float64
    n_steps::Int
    total_time::Float64
    # Dissipative register
    r_D::Int
    N_D::Int
    w0_D::Float64
    t0_D::Float64
    energy_range::Tuple{Float64, Float64}
    # Outer coherent register (b_-)
    r_bm::Int
    N_bm::Int
    w0_bm::Float64
    t0_bm::Float64
    # Inner coherent register (b_+)
    r_bp::Int
    N_bp::Int
    w0_bp::Float64
    t0_bp::Float64
    # GQSP info
    with_gqsp::Bool
    gqsp_degree::Int
    cost_interpretation::Symbol
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

Return QPE grid size, Fourier time step, and represented energy range.
"""
function _qpe_grid_info(r::Int, w0::Real)
    N = 2^r
    t0 = 2π / (N * w0)
    energy_range = (Float64(-N÷2 * w0), Float64((N÷2 - 1) * w0))
    return (; N, t0, energy_range)
end

_channel_cost_interpretation(with_gqsp::Bool) = with_gqsp ?
    :formal_unscaled_gqsp_surrogate : :physical_channel

function Base.show(io::IO, b::SimulationTimeBudget)
    gqsp = b.with_gqsp ? ", gqsp d=$(b.gqsp_degree), formal surrogate" : ""
    print(io, "SimulationTimeBudget(r_D=$(b.r_D), total=$(b.total_time), n_steps=$(b.n_steps), $(b.construction)$gqsp)")
end

function Base.show(io::IO, ::MIME"text/plain", b::SimulationTimeBudget)
    println(io, "SimulationTimeBudget")
    println(io, "  Dissipative   (D):  r=$(b.r_D), N=$(b.N_D), w0=$(b.w0_D), t0=$(b.t0_D)")
    println(io, "  Outer coh.   (b_-): r=$(b.r_bm), N=$(b.N_bm), w0=$(b.w0_bm), t0=$(b.t0_bm)")
    println(io, "  Inner coh.   (b_+): r=$(b.r_bp), N=$(b.N_bp), w0=$(b.w0_bp), t0=$(b.t0_bp)")
    println(io, "  Energy range:       [$(b.energy_range[1]), $(b.energy_range[2])]")
    println(io, "  Physics: β=$(b.beta), σ=$(b.sigma), δ=$(b.delta), n=$(b.n_qubits)")
    println(io, "  Construction: $(b.construction), rescaling=$(b.rescaling_factor)")
    fp = isempty(b.filter_params) ? "" : "($(join(["$k=$v" for (k,v) in b.filter_params], ", ")))"
    println(io, "  Filter: $(b.filter_type)$fp")
    println(io, "  GQSP: with_gqsp=$(b.with_gqsp), gqsp_degree=$(b.gqsp_degree)")
    println(io, "  Cost interpretation: $(b.cost_interpretation)")
    b.with_gqsp && println(io,
        "  Warning: formal query count for the unscaled polynomial surrogate; not a synthesised GQSP circuit")
    println(io, "  T: $(b.T)")
    println(io, "  ─────────────────────────────────")
    println(io, "  OFT time:      $(b.oft_time)")
    println(io, "  B per BE:      $(b.b_per_be)")
    mult_note = b.with_gqsp ? "= 2·$(b.gqsp_degree)·b_per_be" : "= b_per_be (no GQSP)"
    println(io, "  B time:        $(b.b_time)  ($mult_note)")
    println(io, "  Per step:      $(b.per_step_time)  (2×$(b.oft_time) + $(b.b_time))")
    println(io, "  Steps:         $(b.n_steps)")
    print(io,   "  Total:         $(b.total_time)")
end

"""
    _oft_hamiltonian_time(r, w0, transition_weights) -> Float64

Sum weighted absolute evolution times over the full QPE grid.
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

"""
    _b_hamiltonian_time(b_minus, b_plus, beta, sigma, t0_outer, t0_inner) -> Float64

Compute the coherent-correction evolution time using independent outer and
inner grids. Returns zero when either kernel is absent.
"""
function _b_hamiltonian_time(
    b_minus::Union{Nothing, Dict},
    b_plus::Union{Nothing, Dict},
    beta::Real,
    sigma::Real,
    t0_outer::Real,
    t0_inner::Real,
)
    (b_minus === nothing || b_plus === nothing) && return 0.0
    (isempty(b_minus) || isempty(b_plus)) && return 0.0

    # Inner contribution: $sum_s abs(b_+(s)) abs(s beta)$.
    inner_weighted = sum(abs(v) * abs(s * beta) for (s, v) in b_plus)
    inner_norm = sum(abs(v) for (_, v) in b_plus)

    # Outer contribution: $sum_t abs(b_-(t)) abs(t / sigma)$.
    outer_weighted = sum(abs(v) * abs(t / sigma) for (t, v) in b_minus)
    outer_norm = sum(abs(v) for (_, v) in b_minus)

    # Math: $t_- t_+ (4 norm(b_-, 1) I_+ + 2 I_- norm(b_+, 1))$.
    return t0_outer * t0_inner * (4.0 * outer_norm * inner_weighted + 2.0 * outer_weighted * inner_norm)
end

# Compatibility overload for a shared outer/inner grid spacing.
_b_hamiltonian_time(b_minus, b_plus, beta::Real, sigma::Real, t0::Real) =
    _b_hamiltonian_time(b_minus, b_plus, beta, sigma, t0, t0)

function _determine_filter_info(config::Config)
    if !config.with_linear_combination
        params = Dict{Symbol, Float64}()
        config.gaussian_parameters[1] !== nothing && (params[:w_gamma] = Float64(config.gaussian_parameters[1]))
        config.gaussian_parameters[2] !== nothing && (params[:sigma_gamma] = Float64(config.gaussian_parameters[2]))
        return :gaussian, params
    end
    a_val = something(config.a, 0.0)
    s_val = something(config.s, 0.0)
    # `validate_config!` rejects the unsupported `s == 0 && a > 0` case.
    if s_val > 0
        return :smooth_metropolis, Dict{Symbol, Float64}(:a => a_val, :s => s_val)
    else
        return :kinky_metropolis, Dict{Symbol, Float64}()
    end
end

"""
    compute_simulation_time(config, ham, T; n_steps=nothing) -> SimulationTimeBudget

Compute the Hamiltonian-simulation time budget up to target time `T`.

# Arguments
- `config`: thermalisation configuration with independent register triples.
- `ham`: Hamiltonian supplying the physical rescaling factor.
- `T`: target evolution time. It must be finite and positive when `n_steps` is
  omitted, and must equal `n_steps * config.delta` when an exact count is
  supplied (including `T = 0` for zero steps).

# Keywords
- `n_steps`: Exact nonnegative channel-step count. When omitted, use
  `ceil(Int, T / config.delta)` for a general continuous target time. When
  supplied, the budget records the exact effective time `n_steps * delta`.

# Returns
A [`SimulationTimeBudget`](@ref). OFT uses the dissipative register; the
coherent term uses its outer and inner registers. With GQSP, the coherent
cost is formally `2 * gqsp_degree * b_per_be`. The current `with_gqsp=true`
simulator applies an unscaled polynomial without complementary synthesis, so
this is labelled `:formal_unscaled_gqsp_surrogate`, not presented as the cost
of a realizable postselected GQSP circuit.
"""
function compute_simulation_time(
    config::Config{Thermalize, D},
    ham::HamHam,
    T::Real,
    ;
    n_steps::Union{Nothing, Integer}=nothing,
) where {D <: Union{TimeDomain, TrotterDomain}}
    delta = config.delta
    # OFT and the two coherent kernels use independent register triples.
    r_D  = register_r_D(config)
    w0_D = register_w0_D(config)

    delta !== nothing && delta > 0 || throw(ArgumentError("config.delta must be set and positive"))
    r_D !== nothing && r_D > 0 || throw(ArgumentError("dissipative register r_D must be set and positive"))
    w0_D !== nothing || throw(ArgumentError("dissipative register w0_D must be set"))
    T_float = Float64(T)
    isfinite(T_float) || throw(ArgumentError("T must be finite"))
    if isnothing(n_steps)
        T_float > 0 || throw(ArgumentError(
            "T must be positive when n_steps is omitted"))
    else
        n_steps >= 0 || throw(ArgumentError(
            "n_steps must be nonnegative when provided"))
    end

    # Dissipative grid
    grid_D = _qpe_grid_info(r_D, w0_D)
    N_D = grid_D.N
    t0_D = grid_D.t0

    # Transition weights on full QPE grid
    energy_labels = _create_energy_labels(r_D, w0_D)
    transition_fn = pick_transition(config)
    # NaN can occur at extreme energies (Inf * 0 in erfc branch); clamp to 0
    transition_weights = Float64[let v = transition_fn(w); isnan(v) ? 0.0 : v end for w in energy_labels]

    # OFT time
    oft_time = _oft_hamiltonian_time(r_D, w0_D, transition_weights)

    # Record coherent registers even when the construction disables their cost.
    r_bm  = register_r_b_minus(config)
    w0_bm = register_w0_b_minus(config)
    t0_bm = register_t0_b_minus(config)
    r_bp  = register_r_b_plus(config)
    w0_bp = register_w0_b_plus(config)
    t0_bp = register_t0_b_plus(config)
    grid_bm = _qpe_grid_info(r_bm, w0_bm)
    grid_bp = _qpe_grid_info(r_bp, w0_bp)
    N_bm = grid_bm.N
    N_bp = grid_bp.N

    # Coherent cost for one block-encoding query, before the GQSP multiplier.
    construction = config.construction
    b_per_be = if with_coherent(construction)
        time_labels_bm = _create_energy_labels(r_bm, w0_bm) .* (t0_bm / w0_bm)
        time_labels_bp = _create_energy_labels(r_bp, w0_bp) .* (t0_bp / w0_bp)
        bm = _compute_truncated_func(_compute_b_minus, time_labels_bm, config.beta, config.sigma)
        bp_fn, bp_args = _select_b_plus_calculator(config)
        bp = _compute_truncated_func(bp_fn, time_labels_bp, bp_args...)
        _b_hamiltonian_time(bm, bp, config.beta, config.sigma, t0_bm, t0_bp)
    else
        0.0
    end

    # Motlagh--Wiebe Form B uses `d` controlled-W and `d` controlled-W† slots.
    # Math: $N_BE = 2 d$ with GQSP, and $N_BE = 1$ for direct exponentiation.
    with_gqsp = config.with_gqsp
    gqsp_degree = config.gqsp_degree
    b_time = (with_gqsp && with_coherent(construction)) ? 2.0 * gqsp_degree * b_per_be : b_per_be

    # Assembly
    per_step_time = 2.0 * oft_time + b_time
    step_count = isnothing(n_steps) ? ceil(Int, T_float / delta) : Int(n_steps)
    effective_T = isnothing(n_steps) ? T_float : Float64(step_count) * Float64(delta)
    isfinite(effective_T) || throw(ArgumentError(
        "n_steps * config.delta must be finite"))
    if n_steps !== nothing
        time_tol = 8eps(Float64) * max(1.0, abs(effective_T))
        isapprox(T_float, effective_T; atol=time_tol, rtol=8eps(Float64)) ||
            throw(ArgumentError(
                "T must agree with n_steps * config.delta when n_steps is provided " *
                "(got T=$T_float, n_steps=$step_count, delta=$delta, expected T=$effective_T)"))
    end
    total_time = Float64(step_count) * per_step_time

    # Metadata
    n_qubits = config.num_qubits
    construction_sym = construction isa KMS ? :KMS : construction isa GNS ? :GNS : :DLL
    filter_type, filter_params = _determine_filter_info(config)

    return SimulationTimeBudget(
        oft_time, b_per_be, b_time, per_step_time, step_count, total_time,
        Int(r_D), Int(N_D), Float64(w0_D), Float64(t0_D), grid_D.energy_range,
        Int(r_bm), Int(N_bm), Float64(w0_bm), Float64(t0_bm),
        Int(r_bp), Int(N_bp), Float64(w0_bp), Float64(t0_bp),
        Bool(with_gqsp), Int(gqsp_degree), _channel_cost_interpretation(with_gqsp),
        Float64(config.beta), Float64(config.sigma), Float64(delta),
        construction_sym, n_qubits, Float64(ham.rescaling_factor), effective_T,
        filter_type, filter_params,
    )
end

"""
    TrotterStepBudget

Gate-level second-order Strang-substep budget for a channel run.

Unlike [`SimulationTimeBudget`](@ref), this counts every controlled QPE
ladder operation without transition-amplitude weighting. Full coherent
register ladders are counted even when truncated kernels have smaller support.
Plain substep counts are a lower bound until control overheads are applied.

# Fields
- `oft_*`, `b_*`, `n_be_queries`: per-pass and per-step substep counts.
- `substeps_per_step`, `n_steps`, `total_substeps`: aggregate counts.
- `blocks_per_step`, `total_blocks`: contiguous ladder blocks used for
  transpiler-intercept corrections.
- `r_*`, `N_*`, `t0_*`, `M_*`: independent register and Strang parameters.
- `cost_interpretation`: `:physical_channel` or
  `:formal_unscaled_gqsp_surrogate`, with the latter marking counts that assume
  a GQSP query pattern but not a synthesized contractive polynomial.
- Remaining fields record GQSP, construction, temperature, and system metadata.
"""
struct TrotterStepBudget
    # Substep counts
    oft_substeps_per_pass::Int
    oft_substeps_per_step::Int
    b_outer_substeps_per_be::Int
    b_inner_substeps_per_be::Int
    b_substeps_per_be::Int
    n_be_queries::Int
    b_substeps_per_step::Int
    substeps_per_step::Int
    n_steps::Int
    total_substeps::Int
    # Contiguous controlled-block (ladder-rung) counts
    blocks_per_step::Int
    total_blocks::Int
    # Dissipative register
    r_D::Int
    N_D::Int
    t0_D::Float64
    M_D::Int
    # Outer coherent register (b_-)
    r_bm::Int
    N_bm::Int
    t0_bm::Float64
    M_bm::Int
    # Inner coherent register (b_+)
    r_bp::Int
    N_bp::Int
    t0_bp::Float64
    M_bp::Int
    # GQSP info
    with_gqsp::Bool
    gqsp_degree::Int
    cost_interpretation::Symbol
    # Physical parameters
    beta::Float64
    sigma::Float64
    delta::Float64
    construction::Symbol
    n_qubits::Int
    rescaling_factor::Float64
    T::Float64
end

function Base.show(io::IO, b::TrotterStepBudget)
    gqsp = b.with_gqsp ? ", gqsp d=$(b.gqsp_degree), formal surrogate" : ""
    print(io, "TrotterStepBudget(total=$(b.total_substeps), per_step=$(b.substeps_per_step), n_steps=$(b.n_steps), $(b.construction)$gqsp)")
end

function Base.show(io::IO, ::MIME"text/plain", b::TrotterStepBudget)
    println(io, "TrotterStepBudget (Strang substeps)")
    println(io, "  Dissipative   (D):  r=$(b.r_D), N=$(b.N_D), t0=$(b.t0_D), M=$(b.M_D)")
    println(io, "  Outer coh.   (b_-): r=$(b.r_bm), N=$(b.N_bm), t0=$(b.t0_bm), M=$(b.M_bm)")
    println(io, "  Inner coh.   (b_+): r=$(b.r_bp), N=$(b.N_bp), t0=$(b.t0_bp), M=$(b.M_bp)")
    println(io, "  Physics: β=$(b.beta), σ=$(b.sigma), δ=$(b.delta), n=$(b.n_qubits)")
    println(io, "  Construction: $(b.construction), GQSP: with_gqsp=$(b.with_gqsp), d=$(b.gqsp_degree)")
    println(io, "  Cost interpretation: $(b.cost_interpretation)")
    b.with_gqsp && println(io,
        "  Warning: formal query count for the unscaled polynomial surrogate; not a synthesized GQSP circuit")
    println(io, "  T: $(b.T)")
    println(io, "  ─────────────────────────────────")
    println(io, "  OFT substeps/pass:   $(b.oft_substeps_per_pass)  (= (N_D−1)·M_D)")
    println(io, "  OFT substeps/step:   $(b.oft_substeps_per_step)  (×2 fwd+bwd)")
    println(io, "  B substeps/BE:       $(b.b_substeps_per_be)  (outer $(b.b_outer_substeps_per_be) + inner $(b.b_inner_substeps_per_be))")
    println(io, "  BE queries/step:     $(b.n_be_queries)")
    println(io, "  Substeps/step:       $(b.substeps_per_step)")
    println(io, "  Blocks/step:         $(b.blocks_per_step)  (ladder rungs)")
    println(io, "  Steps:               $(b.n_steps)")
    print(io,   "  Total substeps:      $(b.total_substeps)")
end

"""
    count_trotter_steps(config, ham, T) -> TrotterStepBudget

Count second-order Strang substeps and contiguous ladder blocks up to time `T`.

The configuration must provide each active leg's register and `M` value. The
count uses full QPE ladders without transition-amplitude weighting or kernel
truncation and returns a [`TrotterStepBudget`](@ref). For `with_gqsp=true`,
the count prices the formal `2 * gqsp_degree` query pattern but is labelled as
an unscaled-polynomial surrogate until contraction and complementary synthesis
are implemented.
"""
function count_trotter_steps(
    config::Config{Thermalize, D},
    ham::HamHam,
    T::Real,
) where {D <: Union{TimeDomain, TrotterDomain}}
    delta = config.delta
    r_D = register_r_D(config)
    M_D = register_M_D(config)

    delta !== nothing && delta > 0 || throw(ArgumentError("config.delta must be set and positive"))
    r_D !== nothing && r_D > 0 || throw(ArgumentError("dissipative register r_D must be set and positive"))
    M_D !== nothing && M_D > 0 || throw(ArgumentError("dissipative Strang substep count M_D must be set and positive"))
    T > 0 || throw(ArgumentError("T must be positive"))

    N_D = 2^r_D
    t0_D = register_t0_D(config)

    # Math: $S_OFT = 2 (N_D - 1) M_D$ per channel step.
    oft_substeps_per_pass = (N_D - 1) * M_D
    oft_substeps_per_step = 2 * oft_substeps_per_pass

    construction = config.construction
    coherent = with_coherent(construction)

    if coherent
        r_bm = register_r_b_minus(config)
        M_bm = register_M_b_minus(config)
        r_bp = register_r_b_plus(config)
        M_bp = register_M_b_plus(config)
        r_bm !== nothing && r_bm > 0 || throw(ArgumentError("outer coherent register r_b_minus must be set and positive"))
        M_bm !== nothing && M_bm > 0 || throw(ArgumentError("outer coherent Strang substep count M_b_minus must be set and positive"))
        r_bp !== nothing && r_bp > 0 || throw(ArgumentError("inner coherent register r_b_plus must be set and positive"))
        M_bp !== nothing && M_bp > 0 || throw(ArgumentError("inner coherent Strang substep count M_b_plus must be set and positive"))
        t0_bm = Float64(register_t0_b_minus(config))
        t0_bp = Float64(register_t0_b_plus(config))
    else
        # No coherent term: registers recorded as zero, b-substeps vanish.
        r_bm = 0; M_bm = 0; t0_bm = 0.0
        r_bp = 0; M_bp = 0; t0_bp = 0.0
    end
    N_bm = coherent ? 2^r_bm : 0
    N_bp = coherent ? 2^r_bp : 0

    # Coherent block encoding uses two outer and four inner ladder passes.
    b_outer_substeps_per_be = coherent ? 2 * (N_bm - 1) * M_bm : 0
    b_inner_substeps_per_be = coherent ? 4 * (N_bp - 1) * M_bp : 0
    b_substeps_per_be = b_outer_substeps_per_be + b_inner_substeps_per_be

    # Formal query model: $N_BE = 2 d$ for the current unscaled GQSP
    # polynomial surrogate, one for a direct exponential, or zero.
    n_be_queries = !coherent ? 0 : (config.with_gqsp ? 2 * config.gqsp_degree : 1)
    b_substeps_per_step = n_be_queries * b_substeps_per_be

    substeps_per_step = oft_substeps_per_step + b_substeps_per_step
    n_steps = ceil(Int, T / delta)
    total_substeps = n_steps * substeps_per_step

    # Contiguous controlled blocks (ladder rungs): 2 OFT passes of r_D rungs;
    # per BE query 2 outer evolutions × r_bm rungs + 3 inner evolutions × r_bp.
    blocks_per_step = 2 * r_D + n_be_queries * (2 * r_bm + 3 * r_bp)
    total_blocks = n_steps * blocks_per_step

    construction_sym = construction isa KMS ? :KMS : construction isa GNS ? :GNS : :DLL

    return TrotterStepBudget(
        oft_substeps_per_pass, oft_substeps_per_step,
        b_outer_substeps_per_be, b_inner_substeps_per_be, b_substeps_per_be,
        n_be_queries, b_substeps_per_step,
        substeps_per_step, n_steps, total_substeps,
        blocks_per_step, total_blocks,
        Int(r_D), Int(N_D), Float64(something(t0_D, 0.0)), Int(M_D),
        Int(r_bm), Int(N_bm), t0_bm, Int(M_bm),
        Int(r_bp), Int(N_bp), t0_bp, Int(M_bp),
        Bool(config.with_gqsp), Int(config.gqsp_degree),
        _channel_cost_interpretation(config.with_gqsp),
        Float64(config.beta), Float64(config.sigma), Float64(delta),
        construction_sym, config.num_qubits, Float64(ham.rescaling_factor), Float64(T),
    )
end

"""
    RxxBudget

Native RXX gate-count estimate from an affine transpiler model.

# Fields
- `rxx_total`: plain-evolution lower bound.
- `rxx_total_controlled`: estimate with measured one- and two-control factors.
- `rxx_substep_part`, `rxx_boundary_part`: slope and block-intercept parts.
- `rxx_per_substep`, `rxx_intercept`, `f_ctrl1`, `f_ctrl2`: fitted costs.
- `hamiltonian`, `qiskit_version`: measurement provenance.
- `steps`: underlying [`TrotterStepBudget`](@ref).

The estimate excludes state preparation, transition-weight rotations, and LCU
PREP/SELECT beyond the counted block encodings. Generic-angle transpiler fits
do not credit accidental angle simplifications. If `steps.with_gqsp` is true,
the result is a formal count for the current unscaled polynomial surrogate,
not a resource estimate for a synthesized postselected GQSP circuit.
"""
struct RxxBudget
    rxx_total::Float64               # plain lower bound
    rxx_total_controlled::Float64    # per-pass f1/f2 controlled (defendable)
    rxx_substep_part::Float64        # plain: slope · total_substeps
    rxx_boundary_part::Float64       # plain: intercept · total_blocks
    rxx_per_substep::Float64
    rxx_intercept::Float64
    f_ctrl1::Float64
    f_ctrl2::Float64
    hamiltonian::String
    qiskit_version::String
    steps::TrotterStepBudget
end

function Base.show(io::IO, b::RxxBudget)
    suffix = b.steps.with_gqsp ? ", formal unscaled-GQSP surrogate" : ""
    print(io, "RxxBudget(rxx_total=$(b.rxx_total), controlled=$(b.rxx_total_controlled), $(b.hamiltonian) n=$(b.steps.n_qubits)$suffix)")
end

function Base.show(io::IO, ::MIME"text/plain", b::RxxBudget)
    println(io, "RxxBudget — $(b.hamiltonian), n=$(b.steps.n_qubits) (qiskit $(b.qiskit_version))")
    println(io, "  RXX/substep (slope): $(b.rxx_per_substep), block intercept: $(b.rxx_intercept)")
    println(io, "  Substeps total:      $(b.steps.total_substeps)  (× slope = $(b.rxx_substep_part))")
    println(io, "  Blocks total:        $(b.steps.total_blocks)  (× intercept = $(b.rxx_boundary_part))")
    println(io, "  Control factors:     f1=$(b.f_ctrl1) (OFT fwd), f2=$(b.f_ctrl2) (OFT bwd + GQSP)")
    println(io, "  Cost interpretation: $(b.steps.cost_interpretation)")
    println(io, "  ─────────────────────────────────")
    println(io, "  RXX total (plain):       $(b.rxx_total)   [lower bound]")
    controlled_label = b.steps.with_gqsp ?
        "formal unscaled-GQSP surrogate" : "implemented Ham-sim model"
    print(io,   "  RXX total (controlled):  $(b.rxx_total_controlled)   [$controlled_label]")
end

"""
    load_rxx_table(path=default) -> Dict{Tuple{String, Int}, NamedTuple}

Load a tab-separated RXX transpiler measurement table keyed by Hamiltonian and
qubit count. Entries contain the affine fit, control overheads, geometry, and
Qiskit version. Nine-column tables map their single control factor to
`f_ctrl1` and leave `f_ctrl2=NaN`.
"""
function load_rxx_table(path::AbstractString = joinpath(
        dirname(@__DIR__), "scripts", "output", "qf_5hg", "rxx_per_step.tsv"))
    isfile(path) || throw(ArgumentError("RXX table not found at $path — run scripts/qf_5hg_rxx_per_step.py"))
    table = Dict{Tuple{String, Int}, NamedTuple}()
    for (i, line) in enumerate(eachline(path))
        i == 1 && continue  # header
        isempty(strip(line)) && continue
        f = split(line, '\t')
        key = (String(f[1]), parse(Int, f[2]))
        base = (;
            geometry = String(f[3]),
            rxx_slope_per_substep = parse(Float64, f[4]),
            rxx_intercept = parse(Float64, f[5]),
            rxx_L1 = parse(Float64, f[6]),
            rxx_fit_max_abs_dev = parse(Float64, f[7]),
        )
        if length(f) == 10        # current: f_ctrl1, f_ctrl2, qiskit_version
            table[key] = (; base...,
                f_ctrl1 = parse(Float64, f[8]),
                f_ctrl2 = parse(Float64, f[9]),
                qiskit_version = String(f[10]))
        elseif length(f) == 9     # compatibility schema with one control factor
            table[key] = (; base...,
                f_ctrl1 = parse(Float64, f[8]),
                f_ctrl2 = NaN,
                qiskit_version = String(f[9]))
        else
            throw(ArgumentError("malformed RXX table row ($(length(f)) cols): $line"))
        end
    end
    return table
end

"""
    estimate_rxx_count(config, ham, T; rxx_table, hamiltonian) -> RxxBudget

Estimate native RXX two-qubit gates up to time `T` from a measured cost table.

`rxx_table` comes from [`load_rxx_table`](@ref); `hamiltonian` selects its
family key. Returns an [`RxxBudget`](@ref) containing both the plain lower
bound and the control-adjusted estimate. A `with_gqsp=true` configuration is
explicitly labelled as a formal unscaled-polynomial surrogate through its
underlying step budget.
"""
function estimate_rxx_count(
    config::Config{Thermalize, D},
    ham::HamHam,
    T::Real;
    rxx_table::Dict{Tuple{String, Int}, <:NamedTuple},
    hamiltonian::AbstractString,
) where {D <: Union{TimeDomain, TrotterDomain}}
    steps = count_trotter_steps(config, ham, T)
    key = (String(hamiltonian), config.num_qubits)
    haskey(rxx_table, key) || throw(ArgumentError(
        "no RXX measurement for $key — available: $(sort(collect(keys(rxx_table))))"))
    e = rxx_table[key]
    slope = e.rxx_slope_per_substep
    intercept = e.rxx_intercept
    f1 = e.f_ctrl1
    f2 = e.f_ctrl2

    # Math: $N_RXX = slope S_total + intercept B_total$.
    substep_part = slope * steps.total_substeps
    boundary_part = intercept * steps.total_blocks
    rxx_total = substep_part + boundary_part

    # Forward OFT has one outer control; backward OFT and GQSP have two.
    s_fwd = steps.oft_substeps_per_pass               # forward U  (1 control)
    s_bwd = steps.oft_substeps_per_pass               # backward U† (2 controls)
    s_coh = steps.b_substeps_per_step                 # GQSP        (2 controls)
    b_fwd = steps.r_D                                 # forward ladder rungs
    b_bwd = steps.r_D                                 # backward ladder rungs
    b_coh = steps.blocks_per_step - 2 * steps.r_D     # coherent block rungs
    rxx_ctrl_per_step =
        slope * (f1 * s_fwd + f2 * s_bwd + f2 * s_coh) +
        intercept * (f1 * b_fwd + f2 * b_bwd + f2 * b_coh)
    rxx_total_controlled = steps.n_steps * rxx_ctrl_per_step

    return RxxBudget(
        rxx_total, rxx_total_controlled, substep_part, boundary_part,
        slope, intercept, f1, f2,
        String(hamiltonian), e.qiskit_version, steps,
    )
end
