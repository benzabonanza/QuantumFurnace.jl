"""
    load_hamiltonian(type, num_qubits; beta) -> HamHam

Load a pre-computed Hamiltonian from BSON and construct a fully-initialized HamHam.

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
    _load_hamiltonian_bson(path, beta) -> HamHam

Low-level BSON loader that handles legacy HamHam serialization format.
Extracts raw fields from BSON and constructs a new HamHam via HamHam(NamedTuple, beta).
"""
function _load_hamiltonian_bson(path::String, beta::Float64)
    raw = open(path) do io
        BSON.parse(io)
    end

    ham_raw = raw[:hamiltonian]
    fields = ham_raw[:data]

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
    disordering_term = let dt = BSON.raise_recursive(fields[6], cache, init)
        dt === nothing ? nothing : Vector{Matrix{ComplexF64}}(dt)
    end
    disordering_coeffs = let dc = BSON.raise_recursive(fields[7], cache, init)
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
        disordering_term = disordering_term,
        disordering_coeffs = disordering_coeffs,
        eigvals = eigvals_vec,
        eigvecs = eigvecs_mat,
        nu_min = nu_min,
        shift = shift,
        rescaling_factor = rescaling_factor,
        periodic = periodic,
    )

    return HamHam(raw_nt, beta)
end

function generate_filename(config::AbstractLiouvConfig)
    pic_str = string(typeof(config.domain))
    db_str = (config isa LiouvConfigGNS) ? "GNS" : "KMS"
    
    beta_str = "beta=$(config.beta)"
    a_str = "a=$(config.a)"
    b_str = "b=$(config.b)"
    nqb_str = "n=$(config.num_qubits)"
    B = config.with_coherent ? "B" : "noB"

    return join(["liouv", db_str, pic_str, nqb_str, beta_str, B, a_str, b_str], "_") * ".bson"
end

function generate_filename(config::AbstractThermalizeConfig)
    pic_str = string(typeof(config.domain))
    db_str = (config isa ThermalizeConfigGNS) ? "GNS" : "KMS"

    beta_str = "beta=$(config.beta)"
    a_str = "a=$(config.a)"
    b_str = "b=$(config.b)"
    nqb_str = "n=$(config.num_qubits)"
    B = config.with_coherent ? "B" : "noB"
    mix = "mix=$(config.mixing_time)"

    return join(["alg", db_str, pic_str, nqb_str, beta_str, B, a_str, b_str, mix], "_") * ".bson"
end

function riemann_sum(f::Function, grid::Vector{Float64})
    """Uniform grid, rectangle method"""
    d0 = grid[2] - grid[1]
    return d0 * sum(f, grid)
end

function riemann_sum(fvals::Vector{Float64}, d0::Float64)
    return d0 * sum(fvals)
end

function riemann_sum(fvals::Vector{ComplexF64}, d0::Float64)
    return d0 * sum(fvals)
end

function validate_config!(config::AbstractConfig)
    errors = String[]

    # --- Domain-Specific Validation ---
    _collect_config_errors!(errors, config.domain, config)

    # --- Common Validation Logic ---
    # GNS configs are defined without the coherent correction term B.
    if (config isa Union{LiouvConfigGNS, ThermalizeConfigGNS}) && config.with_coherent
        push!(errors, "GNS configs must have with_coherent=false (no coherent B term in this line).")
    end

    if !(config.with_linear_combination) && config.gaussian_parameters == (nothing, nothing)
        push!(errors, "If with_linear_combination is false, gaussian_parameters must be set.")
    end

    if !(config.with_linear_combination)
        w_gamma, sigma_gamma = config.gaussian_parameters
        if w_gamma === nothing || sigma_gamma === nothing
            push!(errors, "For Gaussian transitions gaussian_parameters=(ω_γ, σ_γ) must be set.")
        else
            rhs = if config isa Union{LiouvConfigGNS, ThermalizeConfigGNS}
                2 * w_gamma / (sigma_gamma^2)
            else
                2 * w_gamma / (config.sigma^2 + sigma_gamma^2)
            end
            parameter_relation_holds = isapprox(config.beta, rhs)
            if !(parameter_relation_holds)
                if config isa Union{LiouvConfigGNS, ThermalizeConfigGNS}
                    push!(errors, "For Gaussian transitions (GNS line) require beta ≈ 2*ω_γ/σ_γ^2")
                else
                    push!(errors, "For Gaussian transitions (KMS line) require beta ≈ 2*ω_γ/(σ^2+σ_γ^2)")
                end
            end
        end
    end

    if config.with_linear_combination && config.a == 0.0
        if config.b != 0.0
            push!(errors, "For linear combinations with b != 0, a must also be non-zero.")
        end
        if config.domain isa Union{TimeDomain, TrotterDomain} && config.with_coherent && (isnothing(config.eta) || config.eta <= 0.0)
            push!(errors, "For linear combinations in the KMS DB case with a=0 in TIME or TROTTER domain, eta must be > 0.")
        end
    end

    # --- Error Throwing ---
    if !isempty(errors)
        error_message = "Invalid configuration found:\n" * join(["  - " * err for err in errors], "\n")
        throw(ArgumentError(error_message))
    end

    return nothing
end

function _collect_config_errors!(errors::Vector{String}, ::BohrDomain, config)
    return # No specific checks
end

function _collect_config_errors!(errors::Vector{String}, ::EnergyDomain, config)
    if isnothing(config.num_energy_bits) || config.num_energy_bits <= 0
        push!(errors, "For EnergyDomain, num_energy_bits must be > 0.")
    end
    if isnothing(config.w0) || config.w0 <= 0.0
        push!(errors, "For EnergyDomain, w0 must be > 0.")
    end
end

function _collect_config_errors!(errors::Vector{String}, ::TimeDomain, config)
    if isnothing(config.num_energy_bits) || config.num_energy_bits <= 0
        push!(errors, "For TimeDomain, num_energy_bits must be > 0.")
    end
    if isnothing(config.t0) || config.t0 <= 0.0
        push!(errors, "For TimeDomain, t0 must be > 0.")
    end
    if isnothing(config.w0) || config.w0 <= 0.0
        push!(errors, "For TimeDomain, w0 must be > 0.")
    end
    if !isnothing(config.t0) && !isnothing(config.w0) && !isnothing(config.num_energy_bits) &&
       !isapprox(config.t0 * config.w0, 2pi / 2^config.num_energy_bits)
        push!(errors, "For TimeDomain, the relation t0 * w0 ≈ 2π / 2^N must hold.")
    end
end

function _collect_config_errors!(errors::Vector{String}, ::TrotterDomain, config)
    if isnothing(config.num_energy_bits) || config.num_energy_bits <= 0
        push!(errors, "For TrotterDomain, num_energy_bits must be > 0.")
    end
    if isnothing(config.t0) || config.t0 <= 0.0
        push!(errors, "For TrotterDomain, t0 must be > 0.")
    end
    if isnothing(config.w0) || config.w0 <= 0.0
        push!(errors, "For TrotterDomain, w0 must be > 0.")
    end
    if isnothing(config.num_trotter_steps_per_t0) || config.num_trotter_steps_per_t0 <= 0
        push!(errors, "For TrotterDomain, num_trotter_steps_per_t0 must be > 0.")
    end
    if !isnothing(config.t0) && !isnothing(config.w0) && !isnothing(config.num_energy_bits) &&
       !isapprox(config.t0 * config.w0, 2pi / 2^config.num_energy_bits)
        push!(errors, "For TrotterDomain, the relation t0 * w0 ≈ 2π / 2^N must hold.")
    end
end


function print_press(config::AbstractLiouvConfig)
    params = [
        ("db", (config isa LiouvConfigGNS) ? :GNS : :KMS),
        ("domain", config.domain),
        ("num_qubits", config.num_qubits),
        ("num_energy_bits", config.num_energy_bits),
        ("beta", config.beta),
        ("sigma", config.sigma),
        ("gaussian_parameters", config.gaussian_parameters),
        ("a", config.a),
        ("b", config.b),
        ("eta", config.eta),
        ("t0", config.t0),
        ("w0", config.w0),
        ("with_coherent", config.with_coherent),
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

function print_press(config::AbstractThermalizeConfig)
    params = [
        ("db", (config isa ThermalizeConfigGNS) ? :GNS : :KMS),
        ("domain", config.domain),
        ("num_qubits", config.num_qubits),
        ("num_energy_bits", config.num_energy_bits),
        ("beta", config.beta),
        ("sigma", config.sigma),
        ("gaussian_parameters", config.gaussian_parameters),
        ("a", config.a),
        ("b", config.b),
        ("eta", config.eta),
        ("t0", config.t0),
        ("w0", config.w0),
        ("with_coherent", config.with_coherent),
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