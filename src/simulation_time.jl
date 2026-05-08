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
Stores per-component time costs, **per-term** QPE grid parameters (qf-9z0),
GQSP cost-model flags, physical parameters, and filter configuration — enough
to reproduce a paper table row from a saved BSON.

# Cost fields
- `oft_time`: OFT Hamiltonian simulation time per step
- `b_per_be`: B Hamiltonian-simulation time **per block-encoding query** (legacy
  formula, GQSP-blind). Useful as the audit-side scalar from which the
  GQSP-aware `b_time` derives.
- `b_time`: B coherent correction time per step. Equal to `b_per_be` when
  `with_gqsp = false` (direct `exp(-iδB)`) and to `2 · gqsp_degree · b_per_be`
  when `with_gqsp = true` — Motlagh & Wiebe 2024 Thm. 6 / **Eq. 46** (Form B,
  mixed slot pattern) with the Jacobi-Anger truncation `d = gqsp_degree`:
  `d` controlled-`W` slots + `d` closed-controlled-`W†` slots (the `A'` of MW
  Eq. 45, fires on `|1⟩`) = `2d` block-encoding queries to the walk operator
  `W = R_T · U_{B_a}`. Each `W` invokes one block encoding of `B_a`. The
  earlier Form C realisation (`2d` controlled-`W` + uncontrolled `W^{-d}` tail
  = `3d` queries) is mathematically equivalent (MW Eq. 49→53) but uses 1.5×
  more block-encoding queries; it is being phased out (qf-e4z.19). `0.0` for
  GNS (no coherent term).
- `per_step_time`: `2 × oft_time + b_time`
- `n_steps`: `ceil(T / delta)`
- `total_time`: `n_steps × per_step_time`

# Per-term register info (qf-9z0)
- Dissipative: `r_D`, `N_D`, `w0_D`, `t0_D`, `energy_range`
- Outer coherent (`b_-`): `r_bm`, `N_bm`, `w0_bm`, `t0_bm`
- Inner coherent (`b_+`): `r_bp`, `N_bp`, `w0_bp`, `t0_bp`

# GQSP info
- `with_gqsp`: cost-model flag (matches `config.with_gqsp`)
- `gqsp_degree`: bilateral Jacobi-Anger truncation index `d`. Default `1`.

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
    print(io, "SimulationTimeBudget(r_D=$(b.r_D), total=$(b.total_time), n_steps=$(b.n_steps), $(b.construction)$(b.with_gqsp ? ", gqsp d=$(b.gqsp_degree)" : ""))")
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
    _b_hamiltonian_time(b_minus, b_plus, beta, sigma, t0_outer, t0_inner) -> Float64

Compute the B coherent correction Hamiltonian simulation time as a nested
Riemann sum over truncated b-dict entries with **independent outer/inner**
grid spacings (qf-9z0).

Returns `0.0` when either dictionary is `nothing` or empty (GNS case).

For KMS, each (t,s) pair in the B operator involves:
- Inner: 3 time evolutions of total duration 4|sβ|
- Outer: 2 time evolutions of total duration 2|t/σ|

Total: `t0_outer × t0_inner × Σ_t Σ_s |b_minus[t]| × |b_plus[s]| × (4|sβ| + 2|t/σ|)`
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

    # Inner contribution (b_plus): Σ_s |b_plus[s]| × |s × β|
    inner_weighted = sum(abs(v) * abs(s * beta) for (s, v) in b_plus)
    inner_norm = sum(abs(v) for (_, v) in b_plus)

    # Outer contribution (b_minus): Σ_t |b_minus[t]| × |t / σ|
    outer_weighted = sum(abs(v) * abs(t / sigma) for (t, v) in b_minus)
    outer_norm = sum(abs(v) for (_, v) in b_minus)

    # Factored double sum: t0_outer · t0_inner ×
    #   [4 × ‖b_minus‖₁ × inner_weighted + 2 × outer_weighted × ‖b_plus‖₁]
    return t0_outer * t0_inner * (4.0 * outer_norm * inner_weighted + 2.0 * outer_weighted * inner_norm)
end

# Legacy 5-argument form used by call sites that haven't migrated yet
# (and by historical tests). Forwards to the explicit-outer-inner form with
# `t0_outer = t0_inner = t0`.
_b_hamiltonian_time(b_minus, b_plus, beta::Real, sigma::Real, t0::Real) =
    _b_hamiltonian_time(b_minus, b_plus, beta, sigma, t0, t0)

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

Compute a complete Hamiltonian simulation time budget for the CKG quantum
Gibbs sampling algorithm.

# Arguments
- `config::Config{Thermalize, <:Union{TimeDomain, TrotterDomain}}`: algorithm configuration
- `ham::HamHam`: Hamiltonian (extracts `rescaling_factor`)
- `T::Real`: target simulation time (e.g. mixing time)

# Returns
[`SimulationTimeBudget`](@ref) with `per_step = 2 · oft_time + b_time` and
`total = n_steps · per_step`. The B-time honours `config.with_gqsp`:
GQSP-on multiplies the per-block-encoding cost by `3 · gqsp_degree`
(MW2024 Thm. 6 / Eq. 52, see the struct docstring).

The B-time integrand reads its outer/inner registers separately
(`register_*_b_minus` / `register_*_b_plus` — qf-9z0); the OFT integrand
reads the dissipative `D` register (`register_*_D`). All three triples are
recorded in the returned budget so that downstream sweep BSONs can be
re-priced without re-running anything.
"""
function compute_simulation_time(
    config::Config{Thermalize, D},
    ham::HamHam,
    T::Real,
) where {D <: Union{TimeDomain, TrotterDomain}}
    delta = config.delta
    # Per-register triple plumbing (qf-9z0): the dissipative OFT grid is built
    # from the `D` triple; the B coherent budget below pulls in the `b_minus`
    # and `b_plus` triples separately so the per-term register design is
    # honoured end-to-end (and surfaced in the returned budget).
    r_D  = register_r_D(config)
    w0_D = register_w0_D(config)

    delta !== nothing && delta > 0 || throw(ArgumentError("config.delta must be set and positive"))
    r_D !== nothing && r_D > 0 || throw(ArgumentError("dissipative register r_D must be set and positive"))
    w0_D !== nothing || throw(ArgumentError("dissipative register w0_D must be set"))
    T > 0 || throw(ArgumentError("T must be positive"))

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

    # Per-term coherent grids (recorded regardless of construction; with
    # `with_coherent(construction) = false` the cost is zero but the registers
    # are still part of the config, so we keep them in the budget for audit).
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

    # B coherent correction — per-block-encoding cost. Outer integration uses
    # the `b_minus` register, inner uses the `b_plus` register — each with its
    # own grid (qf-9z0). This is the GQSP-blind formula (one direct
    # `exp(-iδB)` call) and equals one application of the block encoding `U_B`
    # in Hamiltonian-simulation time.
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

    # GQSP cost-model branch — MW2024 Thm. 6 / Eq. 46 (Form B, mixed slots).
    # `with_gqsp = true`: with `d = gqsp_degree`, the symmetric Laurent target
    # `L_d(z) = z^{-d} P(z)` (P of ordinary degree 2d) is realised by
    # `d` controlled-`W` slots + `d` closed-controlled-`W†` slots
    # (the `A' = |0⟩⟨0|⊗I + |1⟩⟨1|⊗U†` of MW Eq. 45; fires on `|1⟩`,
    # *not* the open-controlled `A† = |0⟩⟨0|⊗U† + |1⟩⟨1|⊗I` — the
    # distinction matters because `A' = X_anc · A† · X_anc`) — for a total of
    # `2 · gqsp_degree` block-encoding queries to `W = R_T · U_{B_a}`. Same
    # angles from BS+MW Algorithm 1 transfer directly. Each `W` invokes one
    # block encoding of `B_a` (the joint reflection `R_T` is a Clifford and
    # costs nothing in Ham-sim time). The `2d+1` rotation triples on the QSP
    # ancilla are constant-cost single-qubit gates and contribute nothing to
    # Hamiltonian-simulation time either.
    #
    # NOTE (qf-e4z.19): the Python POC and the v9 thesis figure currently
    # implement the equivalent Form C (Eq. 52: all controlled-`W` + an
    # uncontrolled `W^{-d}` tail), which costs `3d` block-encoding queries —
    # 1.5× more than Form B. This cost model anticipates the Form-B refactor
    # tracked in `qf-e4z.19`; flip the multiplier back to `3.0` if the
    # refactor stalls.
    #
    # `with_gqsp = false`: direct `exp(-iδB)` matrix exponential — one BE.
    with_gqsp = config.with_gqsp
    gqsp_degree = config.gqsp_degree
    b_time = (with_gqsp && with_coherent(construction)) ? 2.0 * gqsp_degree * b_per_be : b_per_be

    # Assembly
    per_step_time = 2.0 * oft_time + b_time
    n_steps = ceil(Int, T / delta)
    total_time = Float64(n_steps) * per_step_time

    # Metadata
    n_qubits = config.num_qubits
    construction_sym = construction isa KMS ? :KMS : construction isa GNS ? :GNS : :DLL
    filter_type, filter_params = _determine_filter_info(config)

    return SimulationTimeBudget(
        oft_time, b_per_be, b_time, per_step_time, n_steps, total_time,
        Int(r_D), Int(N_D), Float64(w0_D), Float64(t0_D), grid_D.energy_range,
        Int(r_bm), Int(N_bm), Float64(w0_bm), Float64(t0_bm),
        Int(r_bp), Int(N_bp), Float64(w0_bp), Float64(t0_bp),
        Bool(with_gqsp), Int(gqsp_degree),
        Float64(config.beta), Float64(config.sigma), Float64(delta),
        construction_sym, n_qubits, Float64(ham.rescaling_factor), Float64(T),
        filter_type, filter_params,
    )
end
