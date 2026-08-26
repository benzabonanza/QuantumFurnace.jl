"""
    load_hamiltonian(type, num_qubits; beta) -> HamHam{Float64}

Load a packaged Hamiltonian fixture and initialise its temperature-dependent
Bohr data and Gibbs state at algorithmic inverse temperature `beta`.
"""
const _PACKAGE_DATA_DIR = normpath(joinpath(@__DIR__, "..", "data"))
const _PACKAGED_HAMILTONIAN_DIR = joinpath(_PACKAGE_DATA_DIR, "hamiltonians")

@inline _package_data_path(parts::AbstractString...) = joinpath(_PACKAGE_DATA_DIR, parts...)

function _packaged_hamiltonian_path(type::AbstractString, num_qubits::Integer)
    type == "heis" || throw(ArgumentError(
        "load_hamiltonian supports type=\"heis\" (the disordered periodic XXX fixture), got \"$type\"."))
    num_qubits > 0 || throw(ArgumentError("num_qubits must be > 0."))
    filename = "heis_xxx_disordered_periodic_n$(num_qubits)_seed46.bson"
    path = joinpath(_PACKAGED_HAMILTONIAN_DIR, filename)
    isfile(path) || throw(ArgumentError(
        "No packaged heis fixture is available for num_qubits=$num_qubits. " *
        "Use build_heis_1d to construct another system."))
    return path
end

function load_hamiltonian(type::AbstractString, num_qubits::Integer; beta::Real)
    return _load_hamiltonian_bson(
        _packaged_hamiltonian_path(type, num_qubits), Float64(beta))
end

"""
    _load_hamiltonian_bson(path, beta) -> HamHam{Float64}

Load a NamedTuple-schema Hamiltonian BSON and construct `HamHam` at `beta`.
"""
function _load_hamiltonian_bson(path::AbstractString, beta::Real)
    return HamHam(_parse_hamiltonian_bson(path), beta)
end

"""
    _parse_hamiltonian_bson(path) -> NamedTuple

Parse a supported Hamiltonian BSON into the canonical raw NamedTuple without
constructing `HamHam`. The result retains `rescaling_factor` for physical-to-
algorithmic temperature conversion and drops trailing metadata fields.
"""
function _parse_hamiltonian_bson(path::AbstractString)
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
              "supported; regenerate it with `build_heis_1d`.")
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
