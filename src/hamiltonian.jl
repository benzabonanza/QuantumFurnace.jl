"""
    HamHam{T<:AbstractFloat}

    Container for Hamiltonian data, spectral decompositions, Bohr frequencies, and Gibbs state.
    Parameterized on element type `T`: real fields use `T`, complex fields use `Complex{T}`.
    Default construction produces `HamHam{Float64}`.

    All fields are fully initialized at construction time -- there are no `Nothing`-typed fields
    for `bohr_freqs`, `bohr_dict`, or `gibbs`. The constructors take a `beta` (inverse temperature)
    parameter and compute these derived quantities directly.

    # Fields
    - `data`: The full Hamiltonian matrix in the computational basis.
    - `bohr_freqs`: Precomputed Bohr frequency matrix (eigvals[i] - eigvals[j]).
    - `bohr_dict`: Mapping from Bohr frequencies to their matrix indices.
    - `base_terms`, `base_coeffs`: The 1, 2 or more site terms that constitute the Hamiltonians, and their uniform coefficients.
    - `disordering_terms`, `disordering_coeffs`: Disordering field terms with per-site coefficients (optional, may be `nothing`).
      Each entry is a Pauli term (e.g. `[Z]` for single-site, `[Z,Z]` for two-site) with a corresponding
      vector of per-site coefficients. Multiple disordering terms are supported (e.g. `[[Z], [Z,Z]]`).
    - `eigvals`, `eigvecs`: Spectral decomposition of the Hamiltonian.
    - `nu_min`: Smallest Bohr frequency in the spectrum, which has to be resolved by all approximations in the algorithm.
    - `shift`, `rescaling_factor`: Values to rescale the spectrum to [0; 0.45].
    - `periodic`: Sets the boundary conditions periodic if `true`.
    - `gibbs`: The theoretical Gibbs state with respect to the Hamiltonian ``\\rho \\propto e^{-\\beta H}``, in the eigenbasis.
"""
struct HamHam{T<:AbstractFloat}
    data::Matrix{Complex{T}}
    bohr_freqs::Matrix{T}
    bohr_dict::Dict{T, Vector{CartesianIndex{2}}}
    base_terms::Vector{Vector{Matrix{Complex{T}}}}
    base_coeffs::Vector{T}
    disordering_terms::Union{Vector{Vector{Matrix{Complex{T}}}}, Nothing}
    disordering_coeffs::Union{Vector{Vector{T}}, Nothing}
    eigvals::Vector{T}
    eigvecs::Matrix{Complex{T}}
    nu_min::T  # Smallest bohr frequency
    shift::T
    rescaling_factor::T
    periodic::Bool
    gibbs::Hermitian{Complex{T}, Matrix{Complex{T}}}
end

"""
    _gibbs_in_eigen(eigvals::Vector{T}, beta::T) where {T<:AbstractFloat} -> Matrix{Complex{T}}

Compute the Gibbs state in the eigenbasis: diagonal matrix with entries exp(-beta * E_i) / Z.
Internal helper used by HamHam constructors before the HamHam object exists.
"""
function _gibbs_in_eigen(eigvals::Vector{T}, beta::T) where {T<:AbstractFloat}
    dim = length(eigvals)
    CT = Complex{T}
    Z = sum(exp.(-beta .* eigvals))
    rho = zeros(CT, dim, dim)
    for i in 1:dim
        rho[i, i] = CT(exp(-beta * eigvals[i]) / Z)
    end
    return rho
end

function HamHam(terms::Vector{Vector{Matrix{ComplexF64}}}, coeffs::Vector{Float64},
    num_qubits::Int64, beta::Float64;
    periodic::Bool = true, hermitian_check = false,
    precision::Type{T} = Float64) where {T<:AbstractFloat}
    """Creates a HamHam{T} object from terms and coefficients, fully initialized with bohr_freqs, bohr_dict, and gibbs."""

    # Mixed-precision policy: downward mismatch errors, upward promotion allowed
    if T !== Float64 && T <: Union{Float16, Float32}
        throw(ArgumentError(
            "Expected $(Complex{T}) term data, got ComplexF64. " *
            "Reconstruct with $(Complex{T}) inputs or use default Float64 precision."))
    end

    hamiltonian_matrix = _construct_base_ham(terms, coeffs, num_qubits; periodic=periodic)

    rescaling_factor, shift = _rescaling_and_shift_factors(hamiltonian_matrix)
    rescaled_hamiltonian::Hermitian{ComplexF64, Matrix{ComplexF64}} = hamiltonian_matrix / rescaling_factor +
                                                                                    shift * I(2^num_qubits)

    rescaled_eigvals, rescaled_eigvecs = eigen(rescaled_hamiltonian)
    rescaled_base_coeffs = coeffs / rescaling_factor
    smallest_bohr_freq = minimum(diff(rescaled_eigvals))

    if hermitian_check
        @assert ishermitian(rescaled_hamiltonian) "The resulting matrix is not Hermitian!"
    end

    bohr_freqs = rescaled_eigvals .- transpose(rescaled_eigvals)
    bohr_dict = create_bohr_dict(bohr_freqs)
    gibbs = Hermitian(_gibbs_in_eigen(rescaled_eigvals, beta))

    return HamHam{T}(
        Matrix(rescaled_hamiltonian),
        bohr_freqs,
        bohr_dict,
        terms,
        rescaled_base_coeffs,
        nothing,  # disordering_terms absent
        nothing,  # disordering_coeffs absent
        rescaled_eigvals,
        rescaled_eigvecs,
        smallest_bohr_freq,
        shift,
        rescaling_factor,
        periodic,
        gibbs,
    )
end

function HamHam(terms::Vector{Vector{Matrix{ComplexF64}}}, coeffs::Vector{Float64},
    disordering_terms::Vector{Vector{Matrix{ComplexF64}}}, disordering_coeffs::Vector{Vector{Float64}},
    num_qubits::Int64, beta::Float64;
    periodic::Bool = true, hermitian_check = false,
    precision::Type{T} = Float64) where {T<:AbstractFloat}
    """Creates a HamHam{T} object from terms, coefficients, and multiple disordering terms, fully initialized."""

    # Mixed-precision policy: downward mismatch errors, upward promotion allowed
    if T !== Float64 && T <: Union{Float16, Float32}
        throw(ArgumentError(
            "Expected $(Complex{T}) term data, got ComplexF64. " *
            "Reconstruct with $(Complex{T}) inputs or use default Float64 precision."))
    end

    if length(disordering_terms) != length(disordering_coeffs)
        throw(ArgumentError("Number of disordering terms must match number of coefficient vectors"))
    end

    base_hamiltonian = _construct_base_ham(terms, coeffs, num_qubits)
    disordering_hamiltonian = _construct_disordering_terms(disordering_terms, disordering_coeffs, num_qubits)
    disordered_ham = base_hamiltonian + disordering_hamiltonian

    rescaling_factor, shift = _rescaling_and_shift_factors(disordered_ham)
    rescaled_hamiltonian::Hermitian{ComplexF64, Matrix{ComplexF64}} = disordered_ham / rescaling_factor +
                                                                                    shift * I(2^num_qubits)

    rescaled_eigvals, rescaled_eigvecs = eigen(rescaled_hamiltonian)
    rescaled_base_coeffs = coeffs / rescaling_factor
    rescaled_disordering_coeffs = [dc / rescaling_factor for dc in disordering_coeffs]
    smallest_bohr_freq = minimum(diff(rescaled_eigvals))

    if hermitian_check
        @assert ishermitian(rescaled_hamiltonian) "The resulting matrix is not Hermitian!"
    end

    bohr_freqs = rescaled_eigvals .- transpose(rescaled_eigvals)
    bohr_dict = create_bohr_dict(bohr_freqs)
    gibbs = Hermitian(_gibbs_in_eigen(rescaled_eigvals, beta))

    return HamHam{T}(
        Matrix(rescaled_hamiltonian),
        bohr_freqs,
        bohr_dict,
        terms,
        rescaled_base_coeffs,
        disordering_terms,
        rescaled_disordering_coeffs,
        rescaled_eigvals,
        rescaled_eigvecs,
        smallest_bohr_freq,
        shift,
        rescaling_factor,
        periodic,
        gibbs,
    )
end

# Single-term convenience: wraps a single disordering term into the multi-term format
function HamHam(terms::Vector{Vector{Matrix{ComplexF64}}}, coeffs::Vector{Float64},
    disordering_term::Vector{Matrix{ComplexF64}}, disordering_coeffs::Vector{Float64},
    num_qubits::Int64, beta::Float64;
    periodic::Bool = true, hermitian_check = false,
    precision::Type{T} = Float64) where {T<:AbstractFloat}

    return HamHam(terms, coeffs, [disordering_term], [disordering_coeffs],
        num_qubits, beta; periodic=periodic, hermitian_check=hermitian_check, precision=precision)
end

"""
    HamHam(raw::NamedTuple, beta) -> HamHam{T}

Construct a fully-initialized HamHam from a NamedTuple of raw data (as returned by
`find_ideal_heisenberg`) plus inverse temperature `beta`.

Infers T from `eltype(raw.eigvals)`. Computes `bohr_freqs`, `bohr_dict`, and `gibbs`
from the raw eigvals.

Supports both legacy single-term format (with `disordering_term` field) and the current
multi-term format (with `disordering_terms` field).
"""
function HamHam(raw::NamedTuple, beta::Real)
    T = eltype(raw.eigvals)
    beta_T = T(beta)
    bohr_freqs = raw.eigvals .- transpose(raw.eigvals)
    bohr_dict = create_bohr_dict(bohr_freqs)
    gibbs = Hermitian(_gibbs_in_eigen(raw.eigvals, beta_T))

    # Handle both legacy single-term and new multi-term format
    dis_terms, dis_coeffs = _unpack_disordering_fields(raw, T)

    return HamHam{T}(
        Matrix{Complex{T}}(raw.matrix),
        bohr_freqs,
        bohr_dict,
        Vector{Vector{Matrix{Complex{T}}}}(raw.terms),
        Vector{T}(raw.base_coeffs),
        dis_terms,
        dis_coeffs,
        Vector{T}(raw.eigvals),
        Matrix{Complex{T}}(raw.eigvecs),
        T(raw.nu_min),
        T(raw.shift),
        T(raw.rescaling_factor),
        raw.periodic,
        gibbs,
    )
end

"""
    _unpack_disordering_fields(raw::NamedTuple, T) -> (terms, coeffs)

Convert disordering fields from a NamedTuple into the multi-term format.
Handles legacy NamedTuples with singular `disordering_term`/`disordering_coeffs` fields
by wrapping them into vectors.
"""
function _unpack_disordering_fields(raw::NamedTuple, ::Type{T}) where {T}
    # New multi-term format
    if haskey(raw, :disordering_terms)
        if raw.disordering_terms === nothing
            return (nothing, nothing)
        end
        terms = [Vector{Matrix{Complex{T}}}(t) for t in raw.disordering_terms]
        coeffs = [Vector{T}(c) for c in raw.disordering_coeffs]
        return (terms, coeffs)
    end

    # Legacy single-term format
    if haskey(raw, :disordering_term)
        if raw.disordering_term === nothing
            return (nothing, nothing)
        end
        terms = [Vector{Matrix{Complex{T}}}(raw.disordering_term)]
        coeffs = [Vector{T}(raw.disordering_coeffs)]
        return (terms, coeffs)
    end

    return (nothing, nothing)
end

"""find_ideal_heisenberg(num_qubits::Int, coeffs::Vector{Float64};
    batch_size::Int=1, periodic::Bool=true,
    disordering_terms=[[Z]]) -> NamedTuple

    Constructs and optimizes a disordered 1D Heisenberg Hamiltonian to maximize the minimum level spacing (smallest Bohr frequency).

    The function generates `batch_size` random realizations of the disordering fields. For each realization, it constructs the Hamiltonian:
    ```math
    H = H_{base} + H_{disorder}
    ```
    where ``H_{base}`` is the Heisenberg chain defined by `coeffs` (XX, YY, ZZ interaction strengths) and ``H_{disorder}`` is the sum of all disordering terms with random per-site coefficients.

    The Hamiltonian is rescaled and shifted to ensure the spectrum fits within specific bounds.

    # Arguments
    - `num_qubits`: The number of sites on the spin chain.
    - `coeffs`: A vector of the uniform interaction strengths for ``\\sigma_x \\sigma_x``, ``\\sigma_y \\sigma_y``, and ``\\sigma_z \\sigma_z`` terms respectively.

    # Keywords
    - `batch_size`: The number of random disorder configurations to sample (default: 1).
    - `periodic`: If `true`, applies periodic boundary conditions to the chain.
    - `disordering_terms`: A vector of Pauli terms for disorder, e.g. `[[Z]]` (default) or `[[Z], [Z,Z]]`.
      Each term gets its own random per-site coefficients.

    # Returns
    - `NamedTuple`: Raw Hamiltonian data with fields: `matrix`, `terms`, `base_coeffs`,
      `disordering_terms`, `disordering_coeffs`, `eigvals`, `eigvecs`, `nu_min`, `shift`,
      `rescaling_factor`, `periodic`. Use `HamHam(raw, beta)` to construct a fully-initialized HamHam.
"""
function find_ideal_heisenberg(num_qubits::Int64,
    coeffs::Vector{Float64}; batch_size::Int64 = 1, periodic::Bool = true,
    disordering_terms::Vector{Vector{Matrix{ComplexF64}}} = [[Z]])

    dim = 2^num_qubits
    terms = [[X, X], [Y, Y], [Z, Z]]

    base_hamiltonian = _construct_base_ham(terms, coeffs, num_qubits; periodic=periodic)

    # Find best config for smallest bohr frequency
    best_nu_min = -1.0
    best_ham_matrix = Matrix{ComplexF64}(undef, 0, 0)
    best_eigvals = Float64[]
    best_eigvecs = Matrix{ComplexF64}(undef, 0, 0)
    best_shift = 0.0
    best_rescaling_factor = 1.0
    best_disordering_coeffs = [Float64[] for _ in disordering_terms]
    all_disordering_coeffs = [zeros(Float64, num_qubits) for _ in disordering_terms]

    p = Progress(batch_size; desc="Optimizing Heisenberg Hamiltonian...")
    for _ in 1:batch_size
        for dc in all_disordering_coeffs
            rand!(dc)
        end
        disordering_ham = _construct_disordering_terms(disordering_terms, all_disordering_coeffs, num_qubits)

        total_ham = base_hamiltonian + disordering_ham
        rescaling_factor, shift = _rescaling_and_shift_factors(total_ham)

        rescaled_ham = (total_ham ./ rescaling_factor) + shift * I

        rescaled_eigvals, rescaled_eigvecs = eigen(Hermitian(rescaled_ham))
        # Check all differences between consecutive eigenvalues
        nu_min = minimum(diff(rescaled_eigvals))
        if nu_min > best_nu_min
            best_nu_min = nu_min
            best_ham_matrix = copy(rescaled_ham)
            best_disordering_coeffs = [copy(dc) for dc in all_disordering_coeffs]
            best_eigvals = rescaled_eigvals
            best_eigvecs = rescaled_eigvecs
            best_shift = shift
            best_rescaling_factor = rescaling_factor

            next!(p, showvalues = [(:nu_min, best_nu_min)])
        else
            next!(p)
        end
    end

    if best_nu_min < 0
        error("Optimization failed to find a valid Hamiltonian")
    end

    return (
        matrix = best_ham_matrix,
        terms = terms,
        base_coeffs = coeffs ./ best_rescaling_factor,
        disordering_terms = disordering_terms,
        disordering_coeffs = [dc ./ best_rescaling_factor for dc in best_disordering_coeffs],
        eigvals = best_eigvals,
        eigvecs = best_eigvecs,
        nu_min = best_nu_min,
        shift = best_shift,
        rescaling_factor = best_rescaling_factor,
        periodic = periodic,
    )
end

function _construct_base_ham(terms::Vector{Vector{Matrix{ComplexF64}}}, coeffs::Vector{Float64},
    num_qubits::Int64; periodic::Bool = true)

    if length(terms) != length(coeffs)
        throw(ArgumentError("The number of terms and coefficients must be equal"))
    end

    hamiltonian::SparseMatrixCSC{ComplexF64} = spzeros(2^num_qubits, 2^num_qubits)
    for (i, term) in enumerate(terms)
        for q in 1:num_qubits
            padded_term = pad_term(term, num_qubits, q; periodic=periodic)  # e.g. term = XX
            hamiltonian += coeffs[i] * padded_term
        end
    end

    return Hermitian(Matrix(hamiltonian))
end

function _construct_disordering_terms(term::Vector{Matrix{ComplexF64}},
    coeffs::Vector{Float64}, num_qubits::Int64)

    if length(coeffs) != num_qubits
        throw(ArgumentError("The number of disordering coeffs must be equal to the number of qubits"))
    end

    disordering_hamiltonian::SparseMatrixCSC{ComplexF64} = spzeros(2^num_qubits, 2^num_qubits)
    for q in 1:num_qubits
        disordering_hamiltonian += coeffs[q] * pad_term(term, num_qubits, q)
    end

    return Hermitian(Matrix(disordering_hamiltonian))
end

function _construct_disordering_terms(terms::Vector{Vector{Matrix{ComplexF64}}},
    coeffs::Vector{Vector{Float64}}, num_qubits::Int64)

    disordering_hamiltonian::SparseMatrixCSC{ComplexF64} = spzeros(2^num_qubits, 2^num_qubits)
    for (term, term_coeffs) in zip(terms, coeffs)
        if length(term_coeffs) != num_qubits
            throw(ArgumentError("Each disordering coefficient vector must have length num_qubits ($num_qubits), got $(length(term_coeffs))"))
        end
        for q in 1:num_qubits
            disordering_hamiltonian += term_coeffs[q] * pad_term(term, num_qubits, q)
        end
    end

    return Hermitian(Matrix(disordering_hamiltonian))
end

function _rescaling_and_shift_factors(hamiltonian::Hermitian)
    """Computes rescaling and shifting factors for a Hamiltonian, s.t. the spectrum is in ``[0, 0.5*(1-ϵ)]`` """

    eps = 0.1  # to avoid 0.5 ~ 0.0 in algorithm
    eigenergies = eigvals(hamiltonian)
    smallest_eigval = minimum(eigenergies)
    largest_eigval = maximum(eigenergies)

    rescaling_factor = (largest_eigval - smallest_eigval) * (2 / (1 - eps))
    shift = - (largest_eigval - smallest_eigval * eps) / (2 * (largest_eigval - smallest_eigval)) + 0.5
    return rescaling_factor, shift
end
