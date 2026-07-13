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
    # (a, s) taxonomy (qf-nq5): kinky Metropolis is exactly (s = 0, a = 0);
    # smooth Metropolis is (s > 0, any a ≥ 0). The (s = 0, a > 0) case is
    # rejected by validate_config!, so we don't need a third branch here.
    if s_val > 0
        return :smooth_metropolis, Dict{Symbol, Float64}(:a => a_val, :s => s_val)
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

# ---------------------------------------------------------------------------
# qf-5hg.2: Trotter-step (gate-level) accounting
# ---------------------------------------------------------------------------

"""
    TrotterStepBudget

Immutable result container for the **gate-level Strang-substep count** of the
CKG channel circuit (qf-5hg.2) — the integer sibling of
[`SimulationTimeBudget`](@ref). Where `compute_simulation_time` sums weighted
Hamiltonian-simulation *time*, [`count_trotter_steps`](@ref) counts the
2nd-order Strang substeps `S_2(δt₀_X)` the circuit physically executes; the
total RXX gate count is then `total_substeps × RXX-per-substep` (qf-5hg.4).

# Accounting model (THE qf-5hg.2 decision, documented here)

The γ(w_k)-weighted OFT time `Σ_k |k·t0_D|·γ(w_k)` in
`compute_simulation_time` is an **expected-time** model: it amortises the
controlled evolutions over the transition-weight amplitudes of the
superposition. Hardware gate counts cannot amortise — every controlled
evolution in the QPE ladder is physically implemented regardless of its
amplitude. One OFT pass is the binary ladder of `r_D` controlled
`e^{iH·2^j·t0_D}` blocks (two's-complement signed grid: the MSB rung carries
weight `−N_D/2`, the others `+2^j`; total evolution duration
`(N_D/2 − 1 + N_D/2)·t0_D = (N_D − 1)·t0_D` either way), i.e.
`(N_D − 1)·M_D` Strang substeps per pass — **not** `Σ_k |k|·M_D` and **not**
γ-weighted. The same ladder rule prices the coherent `B_a` block encoding:
per query, the outer leg runs 2 evolutions `e^{∓iH·t/σ}` (duration weights
1+1 = 2 ladder passes on the `b_-` register) and the inner leg 3 evolutions
`e^{iHτβ}·…·e^{−2iHτβ}·…·e^{iHτβ}` (duration weights 1+2+1 = 4 ladder passes
on the `b_+` register), matching the integrand structure of
`_b_hamiltonian_time` with the amplitude weights `|b_±|` replaced by 1.
Truncation of the `b_±` dicts (`_compute_truncated_func`) is deliberately
**ignored**: the registers are sized at `r_bm`/`r_bp`, so the full ladders
are built; capping at the largest retained label would under-count.

# Strang substep / layer convention

One Strang substep is `S_2(δt₀_X) = A/2·B·A/2` at `δt₀_X = t0_X/M_X`
(per-leg `TrotterTriple` convention, `src/trotter_domain.jl`). `M` composed
substeps telescope, `S_2(t0/M)^M = A/2·(B·A)^{M−1}·B·A/2`, so interior
half-layers merge; the same merge continues across the `2^j` repetitions
inside one ladder rung. Counting therefore works per substep: the
RXX-per-substep constant measured in qf-5hg.3 (asymptotic slope of the
L = 1, 2, 4, 8 linearity fit) already reflects the interior merge, and the
residual boundary correction is O(1) per **contiguous controlled block** =
per ladder rung (the control qubit changes between rungs, blocking further
merging). The rung counts are exposed as `blocks_per_step`/`total_blocks`
so qf-5hg.4 can apply the intercept correction `N_RXX = slope·total_substeps
+ intercept·total_blocks` if desired.

# Controlled-vs-plain caveat

Ancilla-**controlled** `e^{iHt}` costs more than the plain Trotter step on
`{RX, RY, RXX}` (each RXX inside a controlled block decomposes further).
This counter prices **plain** Strang substeps; combined with a plain-step
RXX constant the result is an explicit **lower bound** on the controlled
circuit. qf-5hg.3 may measure the controlled overhead factor separately.

# Fields
- `oft_substeps_per_pass = (N_D − 1)·M_D`; `oft_substeps_per_step = 2×` that
  (forward + backward OFT per δ-step, matching `per_step = 2·oft + b`).
- `b_outer_substeps_per_be = 2·(N_bm − 1)·M_bm`,
  `b_inner_substeps_per_be = 4·(N_bp − 1)·M_bp`, summed in `b_substeps_per_be`.
- `n_be_queries`: `2·gqsp_degree` when `with_gqsp` (MW2024 Thm. 6 / Eq. 46
  Form B — same multiplier as `SimulationTimeBudget.b_time`), `1` for the
  direct `exp(-iδB)`, `0` when `with_coherent(construction) = false`.
- `substeps_per_step = oft_substeps_per_step + n_be_queries·b_substeps_per_be`
- `n_steps = ceil(T/δ)`, `total_substeps = n_steps·substeps_per_step`.
- `blocks_per_step = 2·r_D + n_be_queries·(2·r_bm + 3·r_bp)`, `total_blocks`.
- Register triples `(r_X, N_X, t0_X, M_X)` for X ∈ {D, b_minus, b_plus};
  cross-check identity: `oft_substeps_per_pass·(t0_D/M_D) = (N_D−1)·t0_D`,
  `b_outer_substeps_per_be·(t0_bm/σ)/M_bm = 2·(N_bm−1)·t0_bm/σ`,
  `b_inner_substeps_per_be·(β·t0_bp)/M_bp = 4·(N_bp−1)·β·t0_bp` — the
  unweighted ladder durations of the corresponding budget components.
- GQSP flags, physics (`beta`, `sigma`, `delta`), `construction`, `n_qubits`,
  `rescaling_factor`, `T` — as in `SimulationTimeBudget`.
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
    print(io, "TrotterStepBudget(total=$(b.total_substeps), per_step=$(b.substeps_per_step), n_steps=$(b.n_steps), $(b.construction)$(b.with_gqsp ? ", gqsp d=$(b.gqsp_degree)" : ""))")
end

function Base.show(io::IO, ::MIME"text/plain", b::TrotterStepBudget)
    println(io, "TrotterStepBudget (Strang substeps)")
    println(io, "  Dissipative   (D):  r=$(b.r_D), N=$(b.N_D), t0=$(b.t0_D), M=$(b.M_D)")
    println(io, "  Outer coh.   (b_-): r=$(b.r_bm), N=$(b.N_bm), t0=$(b.t0_bm), M=$(b.M_bm)")
    println(io, "  Inner coh.   (b_+): r=$(b.r_bp), N=$(b.N_bp), t0=$(b.t0_bp), M=$(b.M_bp)")
    println(io, "  Physics: β=$(b.beta), σ=$(b.sigma), δ=$(b.delta), n=$(b.n_qubits)")
    println(io, "  Construction: $(b.construction), GQSP: with_gqsp=$(b.with_gqsp), d=$(b.gqsp_degree)")
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

Count the 2nd-order Strang substeps the CKG channel circuit executes to
simulate up to time `T` (e.g. a mixing time) at δ-step `config.delta` —
the gate-level sibling of [`compute_simulation_time`](@ref) (qf-5hg.2).
See the [`TrotterStepBudget`](@ref) docstring for the full accounting model
(QPE-ladder counting, no γ-amortisation, no `b_±` truncation, Strang merge
and controlled-vs-plain conventions).

Requires the per-leg Strang substep counts (`register_M_D`, and for coherent
constructions `register_M_b_minus`/`register_M_b_plus` — `TrotterTriple`
convention) to be set on the config, alongside the usual register triples.

Cross-check contract (qf-5hg.5 sanity gate): each substep count × its leg's
substep duration `t0_X/M_X` reproduces the **unweighted ladder duration** of
the corresponding `SimulationTimeBudget` component — by construction here,
and asserted against an independently-built budget in the tests.
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

    # OFT: one QPE-ladder pass = (N_D − 1)·M_D substeps; ×2 per δ-step.
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

    # B block encoding: outer 2 ladder passes (e^{∓iH t/σ}), inner 4
    # ladder-pass equivalents (duration weights 1 + 2 + 1 in |τβ|).
    b_outer_substeps_per_be = coherent ? 2 * (N_bm - 1) * M_bm : 0
    b_inner_substeps_per_be = coherent ? 4 * (N_bp - 1) * M_bp : 0
    b_substeps_per_be = b_outer_substeps_per_be + b_inner_substeps_per_be

    # BE queries per δ-step: GQSP Form B (2d) | direct exponential (1) | none.
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
        Float64(config.beta), Float64(config.sigma), Float64(delta),
        construction_sym, config.num_qubits, Float64(ham.rescaling_factor), Float64(T),
    )
end

# ---------------------------------------------------------------------------
# qf-5hg.4: RXX gate-count estimator (QVLS-Q1 native two-qubit gates)
# ---------------------------------------------------------------------------

"""
    RxxBudget

Immutable result of [`estimate_rxx_count`](@ref) — the total native RXX
two-qubit gate count of the CKG channel circuit on the QVLS-Q1 ion-trap gate
set `{RX, RY, RXX}` (Schmale et al., IEEE QSW 2022, arXiv:2206.00544), with
the full defendable chain of factors.

# Composition (qf-5hg.4)

The qf-5hg.3 Qiskit measurement (transpile at `optimization_level = 3` to
`basis_gates = [rx, ry, rxx]`, exact affine fits over L = 1, 2, 4, 8 Strang
substeps) gives, per Hamiltonian and `n`,

    RXX(L contiguous substeps) = slope · L + intercept.

Every contiguous controlled-evolution block of the channel circuit (one QPE
ladder rung — interior substeps consolidate under the transpiler, block
boundaries don't, because the control qubit changes) therefore costs
`slope · L_block + intercept`, and summing over the circuit:

    rxx_total = slope · total_substeps + intercept · total_blocks,

with `total_substeps`/`total_blocks` from [`count_trotter_steps`](@ref)
(embedded as `steps` for the full per-component breakdown: `n_steps` δ-steps,
OFT substeps ×2/step, B block-encoding substeps × `2·gqsp_degree` GQSP
Form-B queries, per-leg `M_D`/`M_bm`/`M_bp` Strang substep counts).

# Two reported numbers (qf-5hg.4 + the control-layer refinement)

`rxx_total` — **plain-evolution lower bound**: every Strang substep priced as
an uncontrolled `exp(-iH·δt)`. Unphysical (the OFT is intrinsically
clock-controlled) but a clean, assumption-free floor.

`rxx_total_controlled` — the **defendable Hamiltonian-simulation count** that
prices the control structure of the algorithm's block-encoding applications,
applying a measured per-control-layer overhead `f_ctrl1`/`f_ctrl2` (one/two
outer controls; from the qf-5hg.3 measurement) per pass:

- **OFT forward** block `U` — the QPE binary ladder of clock-register
  -controlled `e^{iH·2^j t0}`: **one** control ⇒ `f_ctrl1` on
  `oft_substeps_per_pass` substeps (+ `r_D` block intercepts).
- **OFT backward** `controlled-U†` — the weak-measurement uncompute,
  controlled by the weak ancilla *on top of* the clock control: **two**
  controls ⇒ `f_ctrl2` on the second `oft_substeps_per_pass`.
- **GQSP coherent** `2d · controlled-W`, `W = R_T·U_{B_a}` — the `b_±`
  register control *and* the QSP-ancilla control on the whole block:
  **two** controls ⇒ `f_ctrl2` on `b_substeps_per_step`.

Hence
`rxx_total_controlled = n_steps · [ slope·(f1·S_fwd + f2·S_bwd + f2·S_coh)
                                  + intercept·(f1·B_fwd + f2·B_bwd + f2·B_coh) ]`
with `S_fwd = S_bwd = oft_substeps_per_pass`, `S_coh = b_substeps_per_step`,
`B_fwd = B_bwd = r_D`, `B_coh = blocks_per_step − 2·r_D` (all per-δ-step
quantities; the whole bracket is then scaled by `n_steps`).

This deliberately **excludes** (user scope 2026-06-08): initial state
preparation and the Boltzmann / `γ(ω)`-weight ancilla rotation (the
single-qubit "weak rotation" — sub-dominant, scheme-dependent), and the
LCU PREP/SELECT machinery beyond the counted block-encoding queries
(`R_T` is Clifford → free). It captures the *majority* of the GQSP +
weak-measurement Hamiltonian-simulation complexity in RXX. `NaN` when the
control factors were not measured for this `n`.

# Caveats carried from qf-5hg.3
- `slope`/`intercept` are structural generic-angle counts (DT* protocol;
  angle-coincidental Weyl-tolerance savings not credited).
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
    print(io, "RxxBudget(rxx_total=$(b.rxx_total), controlled=$(b.rxx_total_controlled), $(b.hamiltonian) n=$(b.steps.n_qubits))")
end

function Base.show(io::IO, ::MIME"text/plain", b::RxxBudget)
    println(io, "RxxBudget — $(b.hamiltonian), n=$(b.steps.n_qubits) (qiskit $(b.qiskit_version))")
    println(io, "  RXX/substep (slope): $(b.rxx_per_substep), block intercept: $(b.rxx_intercept)")
    println(io, "  Substeps total:      $(b.steps.total_substeps)  (× slope = $(b.rxx_substep_part))")
    println(io, "  Blocks total:        $(b.steps.total_blocks)  (× intercept = $(b.rxx_boundary_part))")
    println(io, "  Control factors:     f1=$(b.f_ctrl1) (OFT fwd), f2=$(b.f_ctrl2) (OFT bwd + GQSP)")
    println(io, "  ─────────────────────────────────")
    println(io, "  RXX total (plain):       $(b.rxx_total)   [lower bound]")
    print(io,   "  RXX total (controlled):  $(b.rxx_total_controlled)   [defendable Ham-sim]")
end

"""
    load_rxx_table(path = scripts/output/qf_5hg/rxx_per_step.tsv)
        -> Dict{Tuple{String, Int}, NamedTuple}

Load the qf-5hg.3 Qiskit measurement table, keyed by `(hamiltonian_name, n)`.
Each entry carries `rxx_slope_per_substep`, `rxx_intercept`, `rxx_L1`,
`rxx_fit_max_abs_dev`, the control-layer overheads `f_ctrl1` (one outer
control) and `f_ctrl2` (two outer controls — both `NaN` where not measured),
`geometry`, and `qiskit_version`. The default path is the committed
measurement produced by `scripts/qf_5hg_rxx_per_step.py`; pass an explicit
path to re-price against a different measurement run.

Accepts both the current 10-column schema (`…, f_ctrl1, f_ctrl2,
qiskit_version`) and the legacy 9-column schema (`…, f_ctrl, qiskit_version`),
mapping the legacy single `f_ctrl` to `f_ctrl1` and setting `f_ctrl2 = NaN`.
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
        elseif length(f) == 9     # legacy: single f_ctrl, qiskit_version
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

Estimate the total native RXX two-qubit gate count for running the CKG
channel up to time `T` (e.g. a mixing time) at the δ-step / register /
Strang parameters of `config` — the qf-5hg epic deliverable. Same signature
family as [`compute_simulation_time`](@ref); like it, this prices saved sweep
configurations without re-running anything.

`rxx_table` is the [`load_rxx_table`](@ref) Dict; `hamiltonian` is the table
key name (`"heis1d_xxx_disordered_periodic_seed46"`). The lookup is
`(hamiltonian, config.num_qubits)`.

Returns both the plain lower bound `rxx_total = slope · total_substeps +
intercept · total_blocks` and the defendable `rxx_total_controlled` that
applies the measured per-control-layer overheads — `f_ctrl1` on the OFT
forward block, `f_ctrl2` on the weak-measurement backward controlled-U† and
the GQSP controlled-W coherent block. See [`RxxBudget`](@ref) for the full
per-pass formula, the documented exclusions (state prep, Boltzmann rotation,
LCU PREP/SELECT), and the carried caveats.
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

    # Plain lower bound.
    substep_part = slope * steps.total_substeps
    boundary_part = intercept * steps.total_blocks
    rxx_total = substep_part + boundary_part

    # Per-pass controlled count. The OFT appears twice per δ-step
    # (forward U + backward controlled-U†); split the lumped
    # oft_substeps_per_step back into the two equal passes. The coherent leg
    # (GQSP controlled-W) is f2 throughout.
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
