"""
    load_hamiltonian(type, num_qubits; beta) -> HamHam{Float64}

Load a pre-computed Hamiltonian from BSON and construct a fully-initialized HamHam{Float64}.

The `beta` keyword is required -- it is used to compute bohr_freqs, bohr_dict, and gibbs
via the `HamHam(NamedTuple, beta)` constructor.

BSON files store legacy HamHam structs with `nothing` for bohr_freqs/bohr_dict/gibbs.
This function uses `BSON.parse` to load the raw field data and reconstructs the HamHam
with the new fully-initialized struct definition.
"""
function load_hamiltonian(type::String, num_qubits::Int; beta::Float64)
    project_root = Pkg.project().path |> dirname
    data_dir = joinpath(project_root, "hamiltonians")
    output_filename = join([type, "disordered", "periodic", "n$num_qubits"], "_") * ".bson"
    ham_path = joinpath(data_dir, output_filename)
    return _load_hamiltonian_bson(ham_path, beta)
end

"""
    _load_hamiltonian_bson(path, beta) -> HamHam{Float64}

Low-level BSON loader that handles two on-disk schemas:

  1. Legacy `HamHam`-typed BSON files where `raw[:hamiltonian]` is a serialised
     struct with 14 positional fields (single-term `disordering_term`).
     Used by the legacy `heis_disordered_periodic_n*.bson` family.

  2. NamedTuple-typed BSON files where `raw[:hamiltonian]` is a NamedTuple
     with multi-term `disordering_terms`. Used by the newer
     `heis_xxx_zzdisordered_periodic_n*` and `heis_xxx_clean_periodic_n*`
     families produced by `hamiltonians/generate_hamiltonians.jl`.

Both paths reconstruct via `HamHam(raw_nt, beta)`, which infers `T` from
`eltype(raw.eigvals)` and recomputes `bohr_freqs`, `bohr_dict`, and `gibbs`.
"""
function _load_hamiltonian_bson(path::String, beta::Float64)
    raw = open(path) do io
        BSON.parse(io)
    end

    ham_raw = raw[:hamiltonian]
    fields = ham_raw[:data]

    # NamedTuple-typed schema (new families) has 11 positional fields and the
    # `:type` blob declares a NamedTuple. Legacy `HamHam`-typed schema has 14
    # positional fields (bohr_freqs/bohr_dict/gibbs slots are `nothing`).
    if length(fields) != 14
        type_name = ham_raw[:type][:name]  # Vector{Any} like ["Core", "NamedTuple"]
        is_namedtuple = type_name isa AbstractVector && !isempty(type_name) &&
                        last(type_name) == "NamedTuple"
        if is_namedtuple
            return HamHam(_namedtuple_schema_to_raw(ham_raw), beta)
        end
        error("Unrecognised Hamiltonian BSON schema (length=$(length(fields)), type=$type_name)")
    end

    # Legacy HamHam field order (14 fields):
    #   1:data, 2:bohr_freqs(nothing), 3:bohr_dict(nothing), 4:base_terms,
    #   5:base_coeffs, 6:disordering_term, 7:disordering_coeffs,
    #   8:eigvals, 9:eigvecs, 10:nu_min, 11:shift, 12:rescaling_factor,
    #   13:periodic, 14:gibbs
    cache = IdDict()
    init = @__MODULE__

    data_matrix = BSON.raise_recursive(fields[1], cache, init)::Matrix{ComplexF64}
    base_terms = Vector{Vector{Matrix{ComplexF64}}}(BSON.raise_recursive(fields[4], cache, init))
    base_coeffs = BSON.raise_recursive(fields[5], cache, init)::Vector{Float64}
    disordering_term_raw = let dt = BSON.raise_recursive(fields[6], cache, init)
        dt === nothing ? nothing : Vector{Matrix{ComplexF64}}(dt)
    end
    disordering_coeffs_raw = let dc = BSON.raise_recursive(fields[7], cache, init)
        dc === nothing ? nothing : Vector{Float64}(dc)
    end
    eigvals_vec = BSON.raise_recursive(fields[8], cache, init)::Vector{Float64}
    eigvecs_mat = BSON.raise_recursive(fields[9], cache, init)::Matrix{ComplexF64}
    nu_min = Float64(fields[10])
    shift = Float64(fields[11])
    rescaling_factor = Float64(fields[12])
    periodic = Bool(fields[13])

    raw_nt = (
        matrix = data_matrix,
        terms = base_terms,
        base_coeffs = base_coeffs,
        disordering_term = disordering_term_raw,
        disordering_coeffs = disordering_coeffs_raw,
        eigvals = eigvals_vec,
        eigvecs = eigvecs_mat,
        nu_min = nu_min,
        shift = shift,
        rescaling_factor = rescaling_factor,
        periodic = periodic,
    )

    return HamHam(raw_nt, beta)
end

"""
    _namedtuple_schema_to_raw(ham_raw::Dict) -> NamedTuple

Convert a `BSON.parse`-d NamedTuple-typed Hamiltonian (new families) into the
NamedTuple shape expected by `HamHam(::NamedTuple, beta)`.

`ham_raw` is the parsed `raw[:hamiltonian]` Dict — a `:tag => :struct` blob
with `:data` a 11-element vector matching the NamedTuple field order:
`(matrix, terms, base_coeffs, disordering_terms, disordering_coeffs,
eigvals, eigvecs, nu_min, shift, rescaling_factor, periodic)`.
"""
function _namedtuple_schema_to_raw(ham_raw::Dict)
    cache = IdDict()
    init = @__MODULE__
    fields = ham_raw[:data]

    matrix = Matrix{ComplexF64}(BSON.raise_recursive(fields[1], cache, init))
    base_terms = [Vector{Matrix{ComplexF64}}(t) for t in BSON.raise_recursive(fields[2], cache, init)]
    base_coeffs = Vector{Float64}(BSON.raise_recursive(fields[3], cache, init))
    disordering_terms = let dt = BSON.raise_recursive(fields[4], cache, init)
        dt === nothing ? nothing : [Vector{Matrix{ComplexF64}}(t) for t in dt]
    end
    disordering_coeffs = let dc = BSON.raise_recursive(fields[5], cache, init)
        dc === nothing ? nothing : [Vector{Float64}(c) for c in dc]
    end
    eigvals_vec = Vector{Float64}(BSON.raise_recursive(fields[6], cache, init))
    eigvecs_mat = Matrix{ComplexF64}(BSON.raise_recursive(fields[7], cache, init))
    nu_min = Float64(fields[8])
    shift = Float64(fields[9])
    rescaling_factor = Float64(fields[10])
    periodic = Bool(fields[11])

    return (
        matrix = matrix,
        terms = base_terms,
        base_coeffs = base_coeffs,
        disordering_terms = disordering_terms,
        disordering_coeffs = disordering_coeffs,
        eigvals = eigvals_vec,
        eigvecs = eigvecs_mat,
        nu_min = nu_min,
        shift = shift,
        rescaling_factor = rescaling_factor,
        periodic = periodic,
    )
end

function _generate_filename(config::Config{Lindbladian})
    pic_str = string(typeof(config.domain))
    db_str = config.construction isa GNS ? "GNS" : "KMS"

    beta_str = "beta=$(config.beta)"
    a_str = "a=$(config.a)"
    s_str = "s=$(config.s)"
    nqb_str = "n=$(config.num_qubits)"
    B = with_coherent(config.construction) ? "B" : "noB"

    return join(["liouv", db_str, pic_str, nqb_str, beta_str, B, a_str, s_str], "_") * ".bson"
end

function _generate_filename(config::Config{Thermalize})
    pic_str = string(typeof(config.domain))
    db_str = config.construction isa GNS ? "GNS" : "KMS"

    beta_str = "beta=$(config.beta)"
    a_str = "a=$(config.a)"
    s_str = "s=$(config.s)"
    nqb_str = "n=$(config.num_qubits)"
    B = with_coherent(config.construction) ? "B" : "noB"
    mix = "mix=$(config.mixing_time)"

    return join(["alg", db_str, pic_str, nqb_str, beta_str, B, a_str, s_str, mix], "_") * ".bson"
end

function _riemann_sum(f::Function, grid::Vector{Float64})
    """Uniform grid, rectangle method"""
    d0 = grid[2] - grid[1]
    return d0 * sum(f, grid)
end

function _riemann_sum(fvals::Vector{Float64}, d0::Float64)
    return d0 * sum(fvals)
end

function _riemann_sum(fvals::Vector{ComplexF64}, d0::Float64)
    return d0 * sum(fvals)
end

"""
    _resolve_filter(config::Config) -> AbstractFilter

Returns `config.filter` if set, otherwise constructs a default
`GaussianFilter(config.sigma)`. Resolves the `Union{Nothing, AbstractFilter}`
field at the precompute boundary so all hot paths see a concrete filter
type and can specialise on it.
"""
@inline _resolve_filter(config::Config) =
    isnothing(config.filter) ? GaussianFilter(config.sigma) : config.filter

"""
    make_trotter_for_config(hamiltonian, config) -> TrottTrott

Build a `TrottTrott` consistent with the per-register grids of `config`.
For `config.domain isa TrotterDomain`:

- **KMS coherent (`with_coherent(config.construction)`)** — calls the
  qf-d0w shared-δt₀ constructor with the natural per-register Trotter steps
  derived from the integration variables in the coherent integrals:

  - `b_-(t/σ)` evolves under `exp(-iH·t/σ)` over the outer grid
    `t = k · register_t0_b_minus(config)` — natural step
    `t0_b_minus_evol = register_t0_b_minus(config) / config.sigma`.
  - `b_+(τβ)` evolves under `exp(-iH·τβ)` over the inner grid
    `τ = k · register_t0_b_plus(config)`  — natural step
    `t0_b_plus_evol  = config.beta · register_t0_b_plus(config)`.

  (For the project convention `σ = 1/β` the two scalings coincide; written in
  general form so that any future runs with `σ ≠ 1/β` still yield exact
  integer step counts in `B_trotter`.)

- **GNS / no-coherent** — calls the legacy single-cache constructor
  `TrottTrott(ham, t0_D, M_user)`; there is no coherent term so the
  per-register caches would be unused.

Throws `ArgumentError` if `config.num_trotter_steps_per_t0` or
`register_t0_D(config)` is `nothing`, or if the integer-M condition fails
(forwarded from the shared-δt₀ constructor).
"""
function make_trotter_for_config(hamiltonian::HamHam, config::Config)
    config.domain isa TrotterDomain ||
        throw(ArgumentError("make_trotter_for_config: config.domain must be TrotterDomain (got $(typeof(config.domain)))."))
    M_user = config.num_trotter_steps_per_t0
    M_user === nothing &&
        throw(ArgumentError("make_trotter_for_config: config.num_trotter_steps_per_t0 must be set."))
    t0_D = register_t0_D(config)
    t0_D === nothing &&
        throw(ArgumentError("make_trotter_for_config: register_t0_D(config) must be set."))
    if with_coherent(config.construction)
        # b_-(t/σ): outer evolution per grid step = t0_grid_b_minus / σ.
        # b_+(τβ): inner evolution per grid step = β · t0_grid_b_plus.
        t0_bm_evol = register_t0_b_minus(config) / config.sigma
        t0_bp_evol = config.beta * register_t0_b_plus(config)
        return TrottTrott(hamiltonian, t0_D, t0_bm_evol, t0_bp_evol, M_user)
    else
        return TrottTrott(hamiltonian, t0_D, M_user)
    end
end

function validate_config!(config::Config)
    errors = String[]

    # --- Domain-Specific Validation ---
    _collect_config_errors!(errors, config)

    # --- Common Validation Logic ---
    # GNS coherent check removed: type system enforces with_coherent(::GNS) = false via trait.

    if !(config.with_linear_combination) && config.gaussian_parameters == (nothing, nothing)
        push!(errors, "If with_linear_combination is false, gaussian_parameters must be set.")
    end

    if !(config.with_linear_combination)
        w_gamma, sigma_gamma = config.gaussian_parameters
        if w_gamma === nothing || sigma_gamma === nothing
            push!(errors, "For Gaussian transitions gaussian_parameters=(ω_γ, σ_γ) must be set.")
        else
            rhs = if config.construction isa GNS
                2 * w_gamma / (sigma_gamma^2)
            else
                2 * w_gamma / (config.sigma^2 + sigma_gamma^2)
            end
            parameter_relation_holds = isapprox(config.beta, rhs)
            if !(parameter_relation_holds)
                if config.construction isa GNS
                    push!(errors, "For Gaussian transitions (GNS line) require beta ≈ 2*ω_γ/σ_γ^2")
                else
                    push!(errors, "For Gaussian transitions (KMS line) require beta ≈ 2*ω_γ/(σ^2+σ_γ^2)")
                end
            end
        end
    end

    if config.with_linear_combination
        a_val = something(config.a, 0.0)
        s_val = something(config.s, 0.0)
        # (a, s) taxonomy: kinky Metropolis is exactly (s = 0, a = 0); smooth
        # Metropolis is (s > 0, any a ≥ 0). The (s = 0, a > 0) combination is
        # an a-regularised but unsmoothed rate that the thesis numerics never
        # use — reject it so we don't silently dispatch into an out-of-scope
        # rate function.
        if s_val == 0.0 && a_val != 0.0
            push!(errors, "For linear combinations require (s = 0, a = 0) for kinky Metropolis or (s > 0) for smooth Metropolis; got (s=0, a=$(a_val)).")
        end
        # Note: the (a=0, s != 0) case is supported in eta-regularized smooth Metro (Task 8).
        if a_val == 0.0 && config.domain isa Union{TimeDomain, TrotterDomain} && with_coherent(config.construction) && (isnothing(config.eta) || config.eta <= 0.0)
            push!(errors, "For linear combinations in the KMS DB case with a=0 in TIME or TROTTER domain, eta must be > 0.")
        end
    end

    # --- GQSP coherent-step validation ---
    if config.with_gqsp
        if !with_coherent(config.construction)
            push!(errors, "with_gqsp requires a construction with coherent term (currently KMS only).")
        end
        # DLL+GQSP is not implemented: _precompute_data for (TimeDomain, DLL) does not
        # produce b_minus/b_plus/gamma_norm_factor that _gqsp_block_encoding_alpha needs.
        if config.construction isa DLL
            push!(errors, "with_gqsp is not supported with DLL construction (no DLL block-encoding norm yet).")
        end
        if !(config.domain isa Union{TimeDomain, TrotterDomain})
            push!(errors, "with_gqsp is only supported for TimeDomain or TrotterDomain.")
        end
        if config.gqsp_degree < 1
            push!(errors, "gqsp_degree must be ≥ 1.")
        end
        if config.gqsp_degree > 100
            push!(errors, "gqsp_degree must be ≤ 100.")
        end
    end

    # --- Jump-selection validation (qf-2vo) ---
    if !(config.jump_selection in (:sweep, :random))
        push!(errors, "jump_selection must be :sweep or :random (got $(config.jump_selection)).")
    end

    # --- DLL filter validation (DLL-1, qf-wmg) ---
    # All DLL filters carry a `beta` that must agree with Config.beta —
    # the filter's KMS factor `e^{-βν/4}` is locked to the simulator's β.
    if config.filter isa DLLGaussianFilter
        beta_tol = 10 * eps(typeof(config.beta))
        if !isapprox(config.filter.beta, config.beta; atol=beta_tol)
            push!(errors, "DLLGaussianFilter.beta must match Config.beta.")
        end
    end
    if config.filter isa DLLMetropolisFilter
        beta_tol = 10 * eps(typeof(config.beta))
        if !isapprox(config.filter.beta, config.beta; atol=beta_tol)
            push!(errors, "DLLMetropolisFilter.beta must match Config.beta.")
        end
        # Bump radius S must be positive — guards against accidental zero/negative.
        if config.filter.S <= 0
            push!(errors, "DLLMetropolisFilter.S must be > 0 (got $(config.filter.S)).")
        end
    end
    # Multi-channel DLL filter (qf-7go.1): every channel must agree with
    # Config.beta and pass its own per-type validation.
    if config.filter isa DLLMultiChannelFilter
        beta_tol = 10 * eps(typeof(config.beta))
        if !isapprox(config.filter.beta, config.beta; atol=beta_tol)
            push!(errors, "DLLMultiChannelFilter.beta must match Config.beta.")
        end
        if isempty(config.filter.channels)
            push!(errors, "DLLMultiChannelFilter must have at least one channel.")
        end
        for (ℓ, ch) in enumerate(config.filter.channels)
            if !hasproperty(ch, :beta)
                push!(errors, "DLLMultiChannelFilter channel $ℓ ($(typeof(ch))) " *
                              "lacks a `beta` field — only DLL-style channels supported.")
                continue
            end
            if !isapprox(ch.beta, config.beta; atol=beta_tol)
                push!(errors, "DLLMultiChannelFilter channel $ℓ.beta=$(ch.beta) " *
                              "does not match Config.beta=$(config.beta).")
            end
            if ch isa DLLMetropolisFilter && ch.S <= 0
                push!(errors, "DLLMultiChannelFilter channel $ℓ.S must be > 0 " *
                              "(got $(ch.S)).")
            end
        end
    end

    # --- DLL construction validation (DLL-2) ---
    if config.construction isa DLL
        # DLL needs an explicit DLL filter at the OFT stage (Eq. 3.4 weighting).
        if config.filter === nothing
            push!(errors, "DLL construction requires an explicit AbstractFilter " *
                          "(e.g. DLLGaussianFilter(beta) or DLLMetropolisFilter(beta)).")
        end
        # EnergyDomain DLL is not in the DLL-2 scope; the paper's EnergyDomain
        # analogue would re-introduce an outer ω-grid and is deferred.
        if config.domain isa EnergyDomain
            push!(errors, "DLL construction is not supported in EnergyDomain (out of scope for DLL-2).")
        end
        # TrotterDomain DLL is deferred until the time-grid quadrature is
        # ported to the Trotter eigvals; see beads epic qf-3i8 notes 2026-04-30.
        if config.domain isa TrotterDomain
            push!(errors, "DLL construction in TrotterDomain is deferred — not yet supported.")
        end
    end

    # --- Error Throwing ---
    if !isempty(errors)
        error_message = "Invalid configuration found:\n" * join(["  - " * err for err in errors], "\n")
        throw(ArgumentError(error_message))
    end

    return nothing
end

"""
    validate_config!(config::Config, ham::HamHam; atol=1e-12, rtol=1e-10)

Two-argument validation: runs the 1-arg `validate_config!(config)` checks
**and**, when `config.beta_phys` is set, enforces the relation

    config.beta == config.beta_phys · ham.rescaling_factor   (β_alg = β_phys · rescale)

within the supplied tolerances. Throws `ArgumentError` on mismatch.

Drivers that author at the *physical* temperature scale (the qf-6vr /
β_phys-first contract) should set both `beta_phys = β_phys` and
`beta = β_phys * ham.rescaling_factor` at construction and call this
2-arg form once the HamHam is in hand, so the pair cannot drift apart.
Callers that do not set `beta_phys` skip the consistency check; this
matches the legacy contract where `cfg.beta` alone is the algorithm-side
β_alg.
"""
function validate_config!(config::Config, ham::HamHam; atol::Real = 1e-12, rtol::Real = 1e-10)
    validate_config!(config)
    if config.beta_phys !== nothing
        expected_beta_alg = config.beta_phys * ham.rescaling_factor
        if !isapprox(config.beta, expected_beta_alg; atol=atol, rtol=rtol)
            throw(ArgumentError(
                "Inconsistent (β_phys, β_alg) pair: config.beta_phys=$(config.beta_phys) and " *
                "ham.rescaling_factor=$(ham.rescaling_factor) imply β_alg=$(expected_beta_alg), " *
                "but config.beta=$(config.beta). Set them at construction so " *
                "`beta == beta_phys * ham.rescaling_factor`."))
        end
    end
    return nothing
end

function _collect_config_errors!(errors::Vector{String}, config::Config{<:Any, BohrDomain})
    return # No specific checks
end

# ---------------------------------------------------------------------------
# Per-register Fourier-relation helpers (qf-9z0).
#
# For a register `X ∈ {D, b_minus, b_plus}` we enforce
#   `t0_X · w0_X ≈ 2π / 2^{r_X}`
# whenever all three are required and non-`nothing` (after legacy fallback).
# `_check_register_fourier!` writes one human-readable error per missing or
# inconsistent field — keeping the message specific to the offending register.
# ---------------------------------------------------------------------------

function _check_register_fourier!(
    errors::Vector{String}, name::AbstractString, r, t0, w0;
    require_t0::Bool = true, require_w0::Bool = true,
)
    if isnothing(r) || r <= 0
        push!(errors, "register '$name': num_energy_bits_$name must be > 0.")
    end
    if require_t0 && (isnothing(t0) || t0 <= 0.0)
        push!(errors, "register '$name': t0_$name must be > 0.")
    end
    if require_w0 && (isnothing(w0) || w0 <= 0.0)
        push!(errors, "register '$name': w0_$name must be > 0.")
    end
    if require_t0 && require_w0 &&
       !isnothing(t0) && !isnothing(w0) && !isnothing(r) &&
       !isapprox(t0 * w0, 2pi / 2^r)
        push!(errors,
              "register '$name': Fourier relation t0_$name * w0_$name ≈ 2π / 2^r_$name must hold (got " *
              "t0=$t0, w0=$w0, r=$r).")
    end
    return errors
end

function _collect_config_errors!(errors::Vector{String}, config::Config{<:Any, EnergyDomain})
    # EnergyDomain dissipator uses analytical A(ω) — only (r_D, w0_D) needed,
    # no t0_D. Coherent term is built in BohrDomain so b_minus/b_plus registers
    # are not consulted here (validate_config! does not require them).
    _check_register_fourier!(
        errors, "D", register_r_D(config), register_t0_D(config), register_w0_D(config);
        require_t0 = false, require_w0 = true,
    )
end

function _collect_config_errors!(errors::Vector{String}, config::Config{<:Any, TimeDomain})
    # CKG/GNS TimeDomain dissipator: full (r_D, t0_D, w0_D) Fourier triple.
    _check_register_fourier!(
        errors, "D", register_r_D(config), register_t0_D(config), register_w0_D(config);
        require_t0 = true, require_w0 = true,
    )
    # Coherent term (KMS only; GNS short-circuited via with_coherent trait).
    if with_coherent(config.construction)
        _check_register_fourier!(
            errors, "b_minus",
            register_r_b_minus(config), register_t0_b_minus(config), register_w0_b_minus(config);
            require_t0 = true, require_w0 = true,
        )
        _check_register_fourier!(
            errors, "b_plus",
            register_r_b_plus(config), register_t0_b_plus(config), register_w0_b_plus(config);
            require_t0 = true, require_w0 = true,
        )
    end
end

# DLL TimeDomain has no ω-grid for the dissipator — `w0_D` is not part of
# the construction (Ding–Li–Lin 2024, Eq. 3.4). Only `r_D` and `t0_D` are
# required. DLL's coherent operator G is built directly on the same dissipative
# time grid via the DLL filter (Eq. 3.7 second equality) — there is no clean
# outer/inner split as in CKG B, so DLL never consumes `b_minus / b_plus`
# registers.
function _collect_config_errors!(errors::Vector{String}, config::Config{<:Any, TimeDomain, DLL})
    _check_register_fourier!(
        errors, "D", register_r_D(config), register_t0_D(config), register_w0_D(config);
        require_t0 = true, require_w0 = false,
    )
end

# DLL TrotterDomain is rejected later in `validate_config!` ("TrotterDomain DLL
# is deferred — not yet supported"). Specialise here only to skip the spurious
# `b_minus / b_plus` register checks that the generic TrotterDomain branch
# would otherwise emit for DLL constructions.
function _collect_config_errors!(errors::Vector{String}, config::Config{<:Any, TrotterDomain, DLL})
    _check_register_fourier!(
        errors, "D", register_r_D(config), register_t0_D(config), register_w0_D(config);
        require_t0 = true, require_w0 = false,
    )
    if isnothing(config.num_trotter_steps_per_t0) || config.num_trotter_steps_per_t0 <= 0
        push!(errors, "For TrotterDomain, num_trotter_steps_per_t0 must be > 0.")
    end
end

function _collect_config_errors!(errors::Vector{String}, config::Config{<:Any, TrotterDomain})
    _check_register_fourier!(
        errors, "D", register_r_D(config), register_t0_D(config), register_w0_D(config);
        require_t0 = true, require_w0 = true,
    )
    if isnothing(config.num_trotter_steps_per_t0) || config.num_trotter_steps_per_t0 <= 0
        push!(errors, "For TrotterDomain, num_trotter_steps_per_t0 must be > 0.")
    end
    if with_coherent(config.construction)
        _check_register_fourier!(
            errors, "b_minus",
            register_r_b_minus(config), register_t0_b_minus(config), register_w0_b_minus(config);
            require_t0 = true, require_w0 = true,
        )
        _check_register_fourier!(
            errors, "b_plus",
            register_r_b_plus(config), register_t0_b_plus(config), register_w0_b_plus(config);
            require_t0 = true, require_w0 = true,
        )
    end
end


function _print_press(config::Config{Lindbladian})
    params = [
        ("db", config.construction isa GNS ? :GNS : (config.construction isa DLL ? :DLL : :KMS)),
        ("domain", config.domain),
        ("num_qubits", config.num_qubits),
        ("r_D", register_r_D(config)),
        ("t0_D", register_t0_D(config)),
        ("w0_D", register_w0_D(config)),
        ("r_b_minus", register_r_b_minus(config)),
        ("t0_b_minus", register_t0_b_minus(config)),
        ("w0_b_minus", register_w0_b_minus(config)),
        ("r_b_plus", register_r_b_plus(config)),
        ("t0_b_plus", register_t0_b_plus(config)),
        ("w0_b_plus", register_w0_b_plus(config)),
        ("beta", config.beta),
        ("sigma", config.sigma),
        ("gaussian_parameters", config.gaussian_parameters),
        ("a", config.a),
        ("s", config.s),
        ("eta", config.eta),
        ("with_coherent", with_coherent(config.construction)),
        ("with_linear_combination", config.with_linear_combination),
        ("num_trotter_steps_per_t0", config.num_trotter_steps_per_t0)
    ]
    provided = filter(p -> p[2] !== nothing, params)
    if isempty(provided)
        return
    end

    println("--- The Press ---")
    for (name, value) in provided
        println("$name: $value")
    end
    println("-----------------")
end

function _print_press(config::Config{Thermalize})
    params = [
        ("db", config.construction isa GNS ? :GNS : (config.construction isa DLL ? :DLL : :KMS)),
        ("domain", config.domain),
        ("num_qubits", config.num_qubits),
        ("r_D", register_r_D(config)),
        ("t0_D", register_t0_D(config)),
        ("w0_D", register_w0_D(config)),
        ("r_b_minus", register_r_b_minus(config)),
        ("t0_b_minus", register_t0_b_minus(config)),
        ("w0_b_minus", register_w0_b_minus(config)),
        ("r_b_plus", register_r_b_plus(config)),
        ("t0_b_plus", register_t0_b_plus(config)),
        ("w0_b_plus", register_w0_b_plus(config)),
        ("beta", config.beta),
        ("sigma", config.sigma),
        ("gaussian_parameters", config.gaussian_parameters),
        ("a", config.a),
        ("s", config.s),
        ("eta", config.eta),
        ("with_coherent", with_coherent(config.construction)),
        ("with_linear_combination", config.with_linear_combination),
        ("num_trotter_steps_per_t0", config.num_trotter_steps_per_t0),
        ("mixing time", config.mixing_time),
        ("delta", config.delta),
    ]
    provided = filter(p -> p[2] !== nothing, params)
    if isempty(provided)
        return
    end

    println("--- The Press ---")
    for (name, value) in provided
        println("$name: $value")
    end
    println("-----------------")
end

function pauli_string_to_matrix(paulistring::Vector{String})
    sigmax::Matrix{ComplexF64} = [0 1; 1 0]
    sigmay::Matrix{ComplexF64} = [0.0 -im; im 0.0]
    sigmaz::Matrix{ComplexF64} = [1 0; 0 -1]

    pauli_matrices::Vector{Matrix{ComplexF64}} = []
    pauli_dict = Dict("X" => sigmax, "Y" => sigmay, "Z" => sigmaz, "I" => Matrix{ComplexF64}(I(2)))
    for pauli_str in paulistring
        push!(pauli_matrices, pauli_dict[pauli_str])
    end
    return pauli_matrices
end

function expm_pauli_padded(pauli_list::Vector{Matrix{ComplexF64}}, coeff::Float64, num_qubits::Int64, position::Int64)
    """Arg e.g. NN terms: [X, X], and it pads it with identities in the rest of the sites. Then creates the expm."""

    padded_term = pad_term(pauli_list, num_qubits, position)
    expm = cos(coeff) * I(2^num_qubits) + 1im * sin(coeff) * padded_term
    return expm
end

function pad_term(terms::Vector{Matrix{ComplexF64}}, num_qubits::Int64, position::Int; periodic::Bool = true)
    
    term_length = length(terms)
    terms = [sparse(term) for term in terms]
    last_position = position + term_length - 1
    # Drop boundary overstepping terms for aperiodic boundary condition 
    if (!(periodic) && last_position > num_qubits)
        return zeros(2^num_qubits, 2^num_qubits)
    end

    if last_position <= num_qubits
        id_before = sparse(I, 2^(position - 1), 2^(position - 1))
        id_after = sparse(I, 2^(num_qubits - last_position), 2^(num_qubits - last_position))
        padded_tensor_list = [id_before, terms..., id_after]
    else
        id_between = sparse(I, 2^(num_qubits - term_length), 2^(num_qubits - term_length))
        not_overflown_terms = terms[1:num_qubits - position + 1]
        overflown_terms = terms[num_qubits - position + 2:end]
        padded_tensor_list = [overflown_terms..., id_between, not_overflown_terms...]
    end

    padded_term::SparseMatrixCSC{ComplexF64} = kron(padded_tensor_list...)
    return padded_term
end