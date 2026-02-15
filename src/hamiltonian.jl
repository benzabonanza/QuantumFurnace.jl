"""
    HamHam

    Container for Hamiltonian data, spectral decompositions, Bohr frequencies, and Gibbs state.

    All fields are fully initialized at construction time -- there are no `Nothing`-typed fields
    for `bohr_freqs`, `bohr_dict`, or `gibbs`. The constructors take a `beta` (inverse temperature)
    parameter and compute these derived quantities directly.

    # Fields
    - `data`: The full Hamiltonian matrix in the computational basis.
    - `bohr_freqs`: Precomputed Bohr frequency matrix (eigvals[i] - eigvals[j]).
    - `bohr_dict`: Mapping from Bohr frequencies to their matrix indices.
    - `base_terms`, `base_coeffs`: The 1, 2 or more site terms that constitute the Hamiltonians, and their uniform coefficients.
    - `disordering_term`, `disordering_coeffs`: Some external field term, that can have different coeffs. on each site (optional, may be `nothing`).
    - `eigvals`, `eigvecs`: Spectral decomposition of the Hamiltonian.
    - `nu_min`: Smallest Bohr frequency in the spectrum, which has to be resolved by all approximations in the algorithm.
    - `shift`, `rescaling_factor`: Values to rescale the spectrum to [0; 0.45].
    - `periodic`: Sets the boundary conditions periodic if `true`.
    - `gibbs`: The theoretical Gibbs state with respect to the Hamiltonian ``\\rho \\propto e^{-\\beta H}``, in the eigenbasis.
"""
struct HamHam
    data::Matrix{ComplexF64}
    bohr_freqs::Matrix{Float64}
    bohr_dict::Dict{Float64, Vector{CartesianIndex{2}}}
    base_terms::Vector{Vector{Matrix{ComplexF64}}}
    base_coeffs::Vector{Float64}
    disordering_term::Union{Vector{Matrix{ComplexF64}}, Nothing}
    disordering_coeffs::Union{Vector{Float64}, Nothing}
    eigvals::Vector{Float64}
    eigvecs::Matrix{ComplexF64}
    nu_min::Float64  # Smallest bohr frequency
    shift::Float64
    rescaling_factor::Float64
    periodic::Bool
    gibbs::Hermitian{ComplexF64, Matrix{ComplexF64}}
end

"""
    _gibbs_in_eigen(eigvals::Vector{Float64}, beta::Float64) -> Matrix{ComplexF64}

Compute the Gibbs state in the eigenbasis: diagonal matrix with entries exp(-beta * E_i) / Z.
Internal helper used by HamHam constructors before the HamHam object exists.
"""
function _gibbs_in_eigen(eigvals::Vector{Float64}, beta::Float64)
    dim = length(eigvals)
    Z = sum(exp.(-beta .* eigvals))
    rho = zeros(ComplexF64, dim, dim)
    for i in 1:dim
        rho[i, i] = exp(-beta * eigvals[i]) / Z
    end
    return rho
end

function HamHam(terms::Vector{Vector{Matrix{ComplexF64}}}, coeffs::Vector{Float64},
    num_qubits::Int64, beta::Float64;
    periodic::Bool = true, hermitian_check = false)
    """Creates a HamHam object from terms and coefficients, fully initialized with bohr_freqs, bohr_dict, and gibbs."""

    hamiltonian_matrix = construct_base_ham(terms, coeffs, num_qubits; periodic=periodic)

    rescaling_factor, shift = rescaling_and_shift_factors(hamiltonian_matrix)
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

    return HamHam(
        rescaled_hamiltonian,
        bohr_freqs,
        bohr_dict,
        terms,
        rescaled_base_coeffs,
        nothing,  # disordering_term absent
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
    disordering_terms::Vector{Matrix{ComplexF64}}, disordering_coeffs::Vector{Float64},
    num_qubits::Int64, beta::Float64;
    periodic::Bool = true, hermitian_check = false)
    """Creates a HamHam object from terms, coefficients, and disordering terms, fully initialized."""

    base_hamiltonian = construct_base_ham(terms, coeffs, num_qubits)
    disordering_hamiltonian = construct_disordering_terms(disordering_terms, disordering_coeffs, num_qubits)
    disordered_ham = base_hamiltonian + disordering_hamiltonian

    rescaling_factor, shift = rescaling_and_shift_factors(disordered_ham)
    rescaled_hamiltonian::Hermitian{ComplexF64, Matrix{ComplexF64}} = disordered_ham / rescaling_factor +
                                                                                    shift * I(2^num_qubits)

    rescaled_eigvals, rescaled_eigvecs = eigen(rescaled_hamiltonian)
    rescaled_base_coeffs = coeffs / rescaling_factor
    rescaled_disordering_coeffs = disordering_coeffs / rescaling_factor
    smallest_bohr_freq = minimum(diff(rescaled_eigvals))

    if hermitian_check
        @assert ishermitian(rescaled_hamiltonian) "The resulting matrix is not Hermitian!"
    end

    bohr_freqs = rescaled_eigvals .- transpose(rescaled_eigvals)
    bohr_dict = create_bohr_dict(bohr_freqs)
    gibbs = Hermitian(_gibbs_in_eigen(rescaled_eigvals, beta))

    return HamHam(
        rescaled_hamiltonian,
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

"""
    HamHam(raw::NamedTuple, beta::Float64) -> HamHam

Construct a fully-initialized HamHam from a NamedTuple of raw data (as returned by
`find_ideal_heisenberg`) plus inverse temperature `beta`.

Computes `bohr_freqs`, `bohr_dict`, and `gibbs` from the raw eigvals.
"""
function HamHam(raw::NamedTuple, beta::Float64)
    bohr_freqs = raw.eigvals .- transpose(raw.eigvals)
    bohr_dict = create_bohr_dict(bohr_freqs)
    gibbs = Hermitian(_gibbs_in_eigen(raw.eigvals, beta))

    return HamHam(
        Hermitian(raw.matrix),
        bohr_freqs,
        bohr_dict,
        raw.terms,
        raw.base_coeffs,
        raw.disordering_term,
        raw.disordering_coeffs,
        raw.eigvals,
        raw.eigvecs,
        raw.nu_min,
        raw.shift,
        raw.rescaling_factor,
        raw.periodic,
        gibbs,
    )
end

"""find_ideal_heisenberg(num_qubits::Int, coeffs::Vector{Float64};
    batch_size::Int=1, periodic::Bool=true) -> NamedTuple

    Constructs and optimizes a disordered 1D Heisenberg Hamiltonian to maximize the minimum level spacing (smallest Bohr frequency).

    The function generates `batch_size` random realizations of a disordering ``Z``-field. For each realization, it constructs the Hamiltonian:
    ```math
    H = H_{base} + H_{disorder}
    ```
    where ``H_{base}`` is the Heisenberg chain defined by `coeffs` (XX, YY, ZZ interaction strengths) and ``H_{disorder}`` is a site-dependent ``Z`` term with random coefficients.

    The Hamiltonian is rescaled and shifted to ensure the spectrum fits within specific bounds.

    # Arguments
    - `num_qubits`: The number of sites on the spin chain.
    - `coeffs`: A vector of the uniform interaction strengths for ``\\sigma_x \\sigma_x``, ``\\sigma_y \\sigma_y``, and ``\\sigma_z \\sigma_z`` terms respectively.

    # Keywords
    - `batch_size`: The number of random disorder configurations to sample (default: 1).
    - `periodic`: If `true`, applies periodic boundary conditions to the chain.

    # Returns
    - `NamedTuple`: Raw Hamiltonian data with fields: `matrix`, `terms`, `base_coeffs`,
      `disordering_term`, `disordering_coeffs`, `eigvals`, `eigvecs`, `nu_min`, `shift`,
      `rescaling_factor`, `periodic`. Use `HamHam(raw, beta)` to construct a fully-initialized HamHam.
"""
function find_ideal_heisenberg(num_qubits::Int64,
    coeffs::Vector{Float64}; batch_size::Int64 = 1, periodic::Bool = true)


    dim = 2^num_qubits
    terms = [[X, X], [Y, Y], [Z, Z]]
    disordering_term = [Z]

    base_hamiltonian = construct_base_ham(terms, coeffs, num_qubits; periodic=periodic)

    # Find best config for smallest bohr frequency
    best_nu_min = -1.0
    best_ham_matrix = Matrix{ComplexF64}(undef, 0, 0)
    best_eigvals = Float64[]
    best_eigvecs = Matrix{ComplexF64}(undef, 0, 0)
    best_shift = 0.0
    best_rescaling_factor = 1.0
    best_disordering_coeffs = Float64[]
    disordering_coeffs = zeros(Float64, num_qubits)

    p = Progress(batch_size; desc="Optimizing Heisenberg Hamiltonian...")
    for _ in 1:batch_size
        rand!(disordering_coeffs)
        disordering_ham = construct_disordering_terms(disordering_term, disordering_coeffs, num_qubits)

        total_ham = base_hamiltonian + disordering_ham
        rescaling_factor, shift = rescaling_and_shift_factors(total_ham)

        rescaled_ham = (total_ham ./ rescaling_factor) + shift * I

        rescaled_eigvals, rescaled_eigvecs = eigen(Hermitian(rescaled_ham))
        # Check all differences between consecutive eigenvalues
        nu_min = minimum(diff(rescaled_eigvals))
        if nu_min > best_nu_min
            best_nu_min = nu_min
            best_ham_matrix = copy(rescaled_ham)
            best_disordering_coeffs = copy(disordering_coeffs)
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
        disordering_term = disordering_term,
        disordering_coeffs = best_disordering_coeffs ./ best_rescaling_factor,
        eigvals = best_eigvals,
        eigvecs = best_eigvecs,
        nu_min = best_nu_min,
        shift = best_shift,
        rescaling_factor = best_rescaling_factor,
        periodic = periodic,
    )
end

function construct_base_ham(terms::Vector{Vector{Matrix{ComplexF64}}}, coeffs::Vector{Float64},
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

function construct_disordering_terms(term::Vector{Matrix{ComplexF64}},
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

function rescaling_and_shift_factors(hamiltonian::Hermitian{ComplexF64, Matrix{ComplexF64}})
    """Computes rescaling and shifting factors for a Hamiltonian, s.t. the spectrum is in ``[0, 0.5*(1-ϵ)]`` """

    eps = 0.1  # to avoid 0.5 ~ 0.0 in algorithm
    eigenergies = eigvals(hamiltonian)
    smallest_eigval = minimum(eigenergies)
    largest_eigval = maximum(eigenergies)

    rescaling_factor = (largest_eigval - smallest_eigval) * (2 / (1 - eps))
    shift = - (largest_eigval - smallest_eigval * eps) / (2 * (largest_eigval - smallest_eigval)) + 0.5
    return rescaling_factor, shift
end
