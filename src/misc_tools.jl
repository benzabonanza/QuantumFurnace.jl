"""
    load_hamiltonian(type, num_qubits; beta) -> HamHam{Float64}

Load a supported precomputed Hamiltonian and initialise its temperature-dependent
Bohr data and Gibbs state at algorithmic inverse temperature `beta`.
"""
function load_hamiltonian(type::String, num_qubits::Int; beta::Float64)
    type == "heis" || error("load_hamiltonian: only type=\"heis\" is supported " *
                            "(maps to heis_xxx_disordered_periodic_n*_seed46.bson). Got: $type")
    project_root = Pkg.project().path |> dirname
    data_dir = joinpath(project_root, "hamiltonians")
    output_filename = "heis_xxx_disordered_periodic_n$(num_qubits)_seed46.bson"
    ham_path = joinpath(data_dir, output_filename)
    return _load_hamiltonian_bson(ham_path, beta)
end

"""
    _load_hamiltonian_bson(path, beta) -> HamHam{Float64}

Load a NamedTuple-schema Hamiltonian BSON and construct `HamHam` at `beta`.
"""
function _load_hamiltonian_bson(path::String, beta::Float64)
    return HamHam(_parse_hamiltonian_bson(path), beta)
end

"""
    _parse_hamiltonian_bson(path) -> NamedTuple

Parse a supported Hamiltonian BSON into the canonical raw NamedTuple without
constructing `HamHam`. The result retains `rescaling_factor` for physical-to-
algorithmic temperature conversion and drops trailing metadata fields.
"""
function _parse_hamiltonian_bson(path::String)
    raw = open(path) do io
        BSON.parse(io)
    end

    ham_raw = raw[:hamiltonian]
    type_name = ham_raw[:type][:name]
    is_namedtuple = type_name isa AbstractVector && !isempty(type_name) &&
                    last(type_name) == "NamedTuple"
    is_namedtuple ||
        error("Unrecognised Hamiltonian BSON schema (expected NamedTuple, got " *
              "type=$type_name). The older `HamHam`-typed schema is no longer " *
              "supported; regenerate via `hamiltonians/generate_hamiltonians.jl`.")
    return _namedtuple_schema_to_raw(ham_raw)
end

"""
    _namedtuple_schema_to_raw(ham_raw::Dict) -> NamedTuple

Convert parsed BSON data to the canonical eleven-field Hamiltonian NamedTuple.
Trailing metadata fields are ignored.
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
    # Uniform-grid rectangle rule.
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
    make_trotter_for_config(hamiltonian, config) -> AbstractTrotter

Build the Trotter cache required by a `TrotterDomain` configuration.

Coherent constructions return a [`TrotterTriple`](@ref) with independent
dissipative, outer-coherent, and inner-coherent Strang steps. Non-coherent
constructions return a single [`TrottTrott`](@ref). Missing register times or
substep counts raise `ArgumentError`.
"""
function make_trotter_for_config(hamiltonian::HamHam, config::Config)
    config.domain isa TrotterDomain ||
        throw(ArgumentError("make_trotter_for_config: config.domain must be TrotterDomain (got $(typeof(config.domain)))."))
    M_user_legacy = config.num_trotter_steps_per_t0
    t0_D = register_t0_D(config)
    t0_D === nothing &&
        throw(ArgumentError("make_trotter_for_config: register_t0_D(config) must be set."))
    if with_coherent(config.construction)
        # Math: $t_(b-) = t0_(b-) / sigma$ and $t_(b+) = beta t0_(b+)$.
        t0_bm_evol = register_t0_b_minus(config) / config.sigma
        t0_bp_evol = config.beta * register_t0_b_plus(config)

        # Accessors resolve per-leg values through the shared compatibility field.
        M_D  = register_M_D(config)
        M_bm = register_M_b_minus(config)
        M_bp = register_M_b_plus(config)
        if M_D === nothing || M_bm === nothing || M_bp === nothing
            throw(ArgumentError(
                "make_trotter_for_config: per-leg M counts must be resolvable; got " *
                "M_D=$M_D, M_b_minus=$M_bm, M_b_plus=$M_bp (set " *
                "`num_trotter_steps_per_t0` for a common default or " *
                "`num_trotter_steps_per_t0_X` for per-leg control)."))
        end
        return TrotterTriple(hamiltonian, t0_D, t0_bm_evol, t0_bp_evol, M_D, M_bm, M_bp)
    else
        M_user_legacy === nothing &&
            throw(ArgumentError("make_trotter_for_config: config.num_trotter_steps_per_t0 must be set (GNS branch)."))
        return TrottTrott(hamiltonian, t0_D, M_user_legacy)
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
        # Smooth Metropolis with `a == 0` requires positive `eta` in time domains.
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

    # Jump-selection validation.
    if !(config.jump_selection in (:sweep, :random))
        push!(errors, "jump_selection must be :sweep or :random (got $(config.jump_selection)).")
    end

    # DLL filter validation.
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
    # Every DLL channel must use the configuration temperature.
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
        # Trotter-domain DLL needs a quadrature defined on Trotter eigenvalues.
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

Validate a configuration and its physical/algorithmic temperature pair.

When `config.beta_phys` is set, require `config.beta` to equal
`config.beta_phys * ham.rescaling_factor` within `atol` and `rtol`. Otherwise
only the one-argument validation runs. Throws `ArgumentError` on failure.
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

# Validate each required register's Fourier relation independently.
# Math: $t0_X w0_X approx 2 pi / 2^(r_X)$.

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

"""
    pauli_string_to_matrix(paulistring) -> Vector{Matrix{ComplexF64}}

Convert labels from `"I"`, `"X"`, `"Y"`, and `"Z"` to single-qubit matrices.

# Arguments
- `paulistring`: Ordered Pauli labels.

# Returns
One matrix per label; the Kronecker product is not formed.
"""
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

"""
    expm_pauli_padded(pauli_list, coeff, num_qubits, position; periodic=true) -> Matrix

Embed a Pauli string and return its exponential.

# Arguments
- `pauli_list`: Consecutive single-qubit Pauli factors.
- `coeff`: Rotation coefficient.
- `num_qubits`: Register size.
- `position`: First site, using one-based indexing.

# Keywords
- `periodic`: Wrap the support across the boundary. For an absent open-boundary
  term, return the identity.

# Returns
The dense unitary on the full register.
"""
function expm_pauli_padded(pauli_list::Vector{Matrix{ComplexF64}}, coeff::Float64, num_qubits::Int64, position::Int64; periodic::Bool=true)
    term_length = length(pauli_list)
    last_position = position + term_length - 1
    if !periodic && last_position > num_qubits
        return Matrix{ComplexF64}(I, 2^num_qubits, 2^num_qubits)
    end

    padded_term = pad_term(pauli_list, num_qubits, position; periodic=periodic)
    # Math: $exp(i theta P) = cos(theta) I + i sin(theta) P$ because $P^2 = I$.
    expm = cos(coeff) * I(2^num_qubits) + 1im * sin(coeff) * padded_term
    return expm
end

"""
    pad_term(terms, num_qubits, position; periodic=true) -> SparseMatrixCSC

Embed consecutive local operators into a qubit register.

# Arguments
- `terms`: Ordered local factors.
- `num_qubits`: Register size.
- `position`: First site, using one-based indexing.

# Keywords
- `periodic`: Wrap the support; an absent open-boundary term returns zero.

# Returns
The sparse full-register operator.
"""
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
