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

    base_hamiltonian = _construct_base_ham(terms, coeffs, num_qubits; periodic=periodic)
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
    HamHam(raw::NamedTuple; beta_phys::Real) -> HamHam{T}

Keyword constructor that takes the **physical** inverse temperature `β_phys`
(against the un-rescaled Hamiltonian) and resolves the algorithm-side
`β_alg = β_phys · raw.rescaling_factor` internally. The stored Gibbs state is
`ρ ∝ exp(-β_alg · H_rescaled)` (== `exp(-β_phys · H_phys)` in physical units).

Use this form when scripts want to type a physical temperature; the positional
`HamHam(raw, beta)` form keeps the legacy interpretation `beta = β_alg`
(against the rescaled spectrum stored in `data` / `eigvals`).
"""
function HamHam(raw::NamedTuple; beta_phys::Real)
    rescale = raw.rescaling_factor
    return HamHam(raw, beta_phys * rescale)
end

"""
    beta_alg(ham::HamHam{T}, beta_phys::Real)   where T<:AbstractFloat -> T
    beta_phys(ham::HamHam{T}, beta_alg::Real)   where T<:AbstractFloat -> T

Convert between physical inverse temperature `β_phys` (against the
un-rescaled Hamiltonian) and algorithm-side `β_alg` (against the rescaled
spectrum stored in `ham.data` / `ham.eigvals`). The relation is

    β_alg = β_phys · ham.rescaling_factor.

`ham.rescaling_factor ≈ 2(λ_max − λ_min)/(1 − ε)` is set at construction by
`find_ideal_heisenberg` to map the physical spectrum into `[0, 0.45]`. For
extensive Hamiltonians (Heisenberg, TFIM) it grows roughly linearly with the
qubit count `n`, so a sweep "at fixed `β_alg` across n" silently varies
`β_phys` by the same factor — and vice versa. Use these helpers whenever
crossing the boundary between user-facing physical temperature and the
algorithm's rescaled units.
"""
beta_alg(ham::HamHam{T}, beta_phys::Real) where {T<:AbstractFloat} = T(beta_phys * ham.rescaling_factor)
beta_phys(ham::HamHam{T}, beta_alg::Real) where {T<:AbstractFloat} = T(beta_alg / ham.rescaling_factor)

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
    disordering_terms=[[Z]], disorder_strength=1.0) -> NamedTuple

    Constructs and optimizes a disordered 1D Heisenberg Hamiltonian to maximize the minimum level spacing (smallest Bohr frequency).

    The function generates `batch_size` random realizations of the disordering fields. For each realization, it constructs the Hamiltonian:
    ```math
    H = H_{base} + H_{disorder}
    ```
    where ``H_{base}`` is the Heisenberg chain defined by `coeffs` (XX, YY, ZZ interaction strengths) and ``H_{disorder}`` is the sum of all disordering terms with random per-site coefficients drawn from `[0, disorder_strength)`.

    The Hamiltonian is rescaled and shifted to ensure the spectrum fits within specific bounds.

    # Arguments
    - `num_qubits`: The number of sites on the spin chain.
    - `coeffs`: A vector of the uniform interaction strengths for ``\\sigma_x \\sigma_x``, ``\\sigma_y \\sigma_y``, and ``\\sigma_z \\sigma_z`` terms respectively.

    # Keywords
    - `batch_size`: The number of random disorder configurations to sample (default: 1).
    - `periodic`: If `true`, applies periodic boundary conditions to the chain.
    - `disordering_terms`: A vector of Pauli terms for disorder, e.g. `[[Z]]` (default) or `[[Z], [Z,Z]]`.
      Each term gets its own random per-site coefficients.
    - `disorder_strength`: Scale factor for the random coefficients (default 1.0). Use ε ≈ 1e-3
      for "clean + ε-disorder" Hamiltonians where disorder only lifts exact degeneracies.

    # Returns
    - `NamedTuple`: Raw Hamiltonian data with fields: `matrix`, `terms`, `base_coeffs`,
      `disordering_terms`, `disordering_coeffs`, `eigvals`, `eigvecs`, `nu_min`, `shift`,
      `rescaling_factor`, `periodic`. Use `HamHam(raw, beta)` to construct a fully-initialized HamHam.
"""
function find_ideal_heisenberg(num_qubits::Int64,
    coeffs::Vector{Float64}; batch_size::Int64 = 1, periodic::Bool = true,
    disordering_terms::Vector{Vector{Matrix{ComplexF64}}} = [[Z]],
    disorder_strength::Float64 = 1.0)

    terms = [[X, X], [Y, Y], [Z, Z]]
    base_hamiltonian = _construct_base_ham(terms, coeffs, num_qubits; periodic=periodic)

    return _optimize_disordered_heisenberg(base_hamiltonian, terms, coeffs, num_qubits,
        disordering_terms; batch_size=batch_size, periodic=periodic,
        disorder_strength=disorder_strength)
end

"""find_ideal_2d_heisenberg(Lx::Int, Ly::Int, coeffs::Vector{Float64};
    batch_size::Int=1, periodic_x::Bool=true, periodic_y::Bool=true,
    disordering_terms=[[Z]], disorder_strength=1e-3) -> NamedTuple

    Constructs and optimizes a 2D anisotropic Heisenberg Hamiltonian on an `Lx × Ly`
    square lattice. Bonds are placed along the right neighbour (x-direction) and the
    up neighbour (y-direction); periodic boundary conditions wrap each direction.

    The Hamiltonian is
    ```math
    H = \\sum_{\\langle i,j\\rangle} \\big(J_x X_i X_j + J_y Y_i Y_j + J_z Z_i Z_j\\big)
        + \\sum_{q,k} c_{k,q} P^{(k)}_q,
    ```
    where ``\\langle i,j\\rangle`` runs over nearest-neighbour bonds in the 2D lattice.
    For Ising-anisotropic couplings (``J_z > J_x = J_y``) the bulk model has a finite-
    temperature thermal phase transition into a Néel-ordered ground state. The 2D
    isotropic model (XXX) does not exhibit a finite-T transition due to Mermin–Wagner.

    Site-to-qubit ordering: site ``(i, j)`` with ``1 \\le i \\le L_x``, ``1 \\le j \\le L_y``
    is mapped to qubit index ``q = (i-1) \\, L_y + (j-1) + 1``. Adjacent ``j`` values
    correspond to consecutive qubit indices; adjacent ``i`` values are separated by ``L_y``.

    # Arguments
    - `Lx`, `Ly`: Lattice dimensions; total qubit count is `num_qubits = Lx * Ly`.
    - `coeffs`: Uniform exchange couplings ``[J_x, J_y, J_z]``.

    # Keywords
    - `batch_size`: Random-disorder realisations to sample (default 1). For very weak
      disorder the optimisation surface is flat — use `batch_size ≥ 200`.
    - `periodic_x`, `periodic_y`: Periodic boundary conditions per direction.
    - `disordering_terms`: Pauli terms used as on-site/two-site disorder (default `[[Z]]`).
    - `disorder_strength`: Per-coefficient scale (default 1e-3 for the "clean + ε" use case).

    # Returns
    - `NamedTuple` with the same schema as `find_ideal_heisenberg` so it is
      compatible with `HamHam(raw, beta)`. The boolean `periodic` field is set to
      `periodic_x && periodic_y`.
"""
function find_ideal_2d_heisenberg(Lx::Int64, Ly::Int64,
    coeffs::Vector{Float64}; batch_size::Int64 = 1,
    periodic_x::Bool = true, periodic_y::Bool = true,
    disordering_terms::Vector{Vector{Matrix{ComplexF64}}} = [[Z]],
    disorder_strength::Float64 = 1e-3)

    if Lx < 1 || Ly < 1
        throw(ArgumentError("Lx and Ly must be at least 1; got Lx=$Lx, Ly=$Ly"))
    end

    num_qubits = Lx * Ly
    terms = [[X, X], [Y, Y], [Z, Z]]
    base_hamiltonian = _construct_2d_heisenberg_base(Lx, Ly, terms, coeffs;
        periodic_x=periodic_x, periodic_y=periodic_y)

    return _optimize_disordered_heisenberg(base_hamiltonian, terms, coeffs, num_qubits,
        disordering_terms; batch_size=batch_size, periodic=(periodic_x && periodic_y),
        disorder_strength=disorder_strength)
end

# Shared inner kernel for find_ideal_heisenberg and find_ideal_2d_heisenberg:
# generate `batch_size` random realisations of the disorder field, keep the one
# with the largest minimum Bohr gap after rescaling/shifting.
function _optimize_disordered_heisenberg(base_hamiltonian::Hermitian{ComplexF64},
    terms::Vector{Vector{Matrix{ComplexF64}}}, coeffs::Vector{Float64},
    num_qubits::Int64, disordering_terms::Vector{Vector{Matrix{ComplexF64}}};
    batch_size::Int64, periodic::Bool, disorder_strength::Float64)

    best_nu_min = -1.0
    best_ham_matrix = Matrix{ComplexF64}(undef, 0, 0)
    best_eigvals = Float64[]
    best_eigvecs = Matrix{ComplexF64}(undef, 0, 0)
    best_shift = 0.0
    best_rescaling_factor = 1.0
    best_disordering_coeffs = [Float64[] for _ in disordering_terms]
    all_disordering_coeffs = [zeros(Float64, num_qubits) for _ in disordering_terms]

    p = Progress(batch_size; desc="Optimizing disordered Heisenberg Hamiltonian...")
    for _ in 1:batch_size
        for dc in all_disordering_coeffs
            rand!(dc)
            dc .*= disorder_strength
        end
        disordering_ham = _construct_disordering_terms(disordering_terms, all_disordering_coeffs, num_qubits)

        total_ham = base_hamiltonian + disordering_ham
        rescaling_factor, shift = _rescaling_and_shift_factors(total_ham)

        rescaled_ham = (total_ham ./ rescaling_factor) + shift * I

        rescaled_eigvals, rescaled_eigvecs = eigen(Hermitian(rescaled_ham))
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

"""find_typical_heisenberg(num_qubits::Int, coeffs::Vector{Float64};
    batch_size::Int=1000, periodic::Bool=true,
    disordering_terms=[[Z]], disorder_strength=1.0) -> NamedTuple

    Construct a 1D disordered Heisenberg Hamiltonian whose normalised sorted spectrum
    is **closest** to the elementwise median sorted spectrum across `batch_size`
    random disorder realisations.

    This selects a "typical" realisation rather than an extremal one. Whereas
    `find_ideal_heisenberg` maximises the minimum Bohr gap (and so can return outlier
    spectra with anomalously large *or* small first gaps at small `n`),
    `find_typical_heisenberg` picks the realisation whose entire sorted spectrum is
    closest in L² to the ensemble median — i.e. the L²-Wasserstein-spectral-median.

    Per-sample work is identical to `find_ideal_heisenberg`: build `H_base + H_disorder`,
    rescale and shift the spectrum to the canonical algorithm range, then diagonalise.
    After the sweep, the elementwise median of the sorted, bandwidth-normalised
    spectra defines the consensus; the winning realisation is the one whose
    spectrum minimises Euclidean distance to that median.

    The returned NamedTuple has the same schema as `find_ideal_heisenberg`, plus an
    extra diagnostic field `typicality_distance` (the L² distance of the chosen
    sample's normalised spectrum to the ensemble-median spectrum). Use
    `HamHam(raw, beta)` to wrap as a `HamHam` — extra fields are ignored by the
    constructor.

    # Arguments
    - `num_qubits`: Spin-chain length.
    - `coeffs`: Uniform exchange couplings ``[J_x, J_y, J_z]``.

    # Keywords
    - `batch_size`: Number of disorder realisations to sample (default 1000). Must
      be ≥ 2 for the W₂ criterion to be meaningful; values ≥ a few hundred are
      recommended for a stable median.
    - `periodic`, `disordering_terms`, `disorder_strength`: same as `find_ideal_heisenberg`.

    # Returns
    NamedTuple with fields `matrix`, `terms`, `base_coeffs`, `disordering_terms`,
    `disordering_coeffs`, `eigvals`, `eigvecs`, `nu_min`, `shift`, `rescaling_factor`,
    `periodic`, `typicality_distance`.
"""
function find_typical_heisenberg(num_qubits::Int64,
    coeffs::Vector{Float64}; batch_size::Int64 = 1000, periodic::Bool = true,
    disordering_terms::Vector{Vector{Matrix{ComplexF64}}} = [[Z]],
    disorder_strength::Float64 = 1.0)

    terms = [[X, X], [Y, Y], [Z, Z]]
    base_hamiltonian = _construct_base_ham(terms, coeffs, num_qubits; periodic=periodic)

    return _optimize_typical_heisenberg(base_hamiltonian, terms, coeffs, num_qubits,
        disordering_terms; batch_size=batch_size, periodic=periodic,
        disorder_strength=disorder_strength)
end

"""find_typical_2d_heisenberg(Lx::Int, Ly::Int, coeffs::Vector{Float64};
    batch_size::Int=1000, periodic_x::Bool=true, periodic_y::Bool=true,
    disordering_terms=[[Z]], disorder_strength=1e-3) -> NamedTuple

    2D analogue of [`find_typical_heisenberg`](@ref) on an `Lx × Ly` lattice with
    nearest-neighbour Heisenberg bonds in both directions. The W₂-spectral-median
    selection criterion is identical; see `find_typical_heisenberg` for details.

    Returns the same NamedTuple schema as `find_ideal_2d_heisenberg`, plus
    `typicality_distance`.
"""
function find_typical_2d_heisenberg(Lx::Int64, Ly::Int64,
    coeffs::Vector{Float64}; batch_size::Int64 = 1000,
    periodic_x::Bool = true, periodic_y::Bool = true,
    disordering_terms::Vector{Vector{Matrix{ComplexF64}}} = [[Z]],
    disorder_strength::Float64 = 1e-3)

    if Lx < 1 || Ly < 1
        throw(ArgumentError("Lx and Ly must be at least 1; got Lx=$Lx, Ly=$Ly"))
    end

    num_qubits = Lx * Ly
    terms = [[X, X], [Y, Y], [Z, Z]]
    base_hamiltonian = _construct_2d_heisenberg_base(Lx, Ly, terms, coeffs;
        periodic_x=periodic_x, periodic_y=periodic_y)

    return _optimize_typical_heisenberg(base_hamiltonian, terms, coeffs, num_qubits,
        disordering_terms; batch_size=batch_size, periodic=(periodic_x && periodic_y),
        disorder_strength=disorder_strength)
end

# Shared inner kernel for find_typical_heisenberg and find_typical_2d_heisenberg.
# Mirrors _optimize_disordered_heisenberg but with a W₂-spectral-median selection
# criterion in place of the greedy max-min-gap criterion. Per-sample work is
# identical (build → rescale → diagonalise); only the across-sample selection
# differs.
function _optimize_typical_heisenberg(base_hamiltonian::Hermitian{ComplexF64},
    terms::Vector{Vector{Matrix{ComplexF64}}}, coeffs::Vector{Float64},
    num_qubits::Int64, disordering_terms::Vector{Vector{Matrix{ComplexF64}}};
    batch_size::Int64, periodic::Bool, disorder_strength::Float64)

    if batch_size < 2
        throw(ArgumentError(
            "find_typical_* requires batch_size >= 2 (got $batch_size); the W₂ selector " *
            "compares against an ensemble median and needs at least two samples"))
    end

    d = 2^num_qubits

    # Per-sample storage: normalised sorted spectra (d × batch_size) and disorder
    # coefficients. We only need eigvals during the sweep — full eigen()/eigvecs
    # is computed once at the end on the chosen sample.
    normalized_specs = Matrix{Float64}(undef, d, batch_size)
    disordering_coeffs_by_sample = Vector{Vector{Vector{Float64}}}(undef, batch_size)
    valid_mask = falses(batch_size)

    p = Progress(batch_size; desc="Sampling Hamiltonians for W2-typical selection...")
    for k in 1:batch_size
        sample_coeffs = [zeros(Float64, num_qubits) for _ in disordering_terms]
        for dc in sample_coeffs
            rand!(dc)
            dc .*= disorder_strength
        end
        disordering_ham = _construct_disordering_terms(disordering_terms, sample_coeffs, num_qubits)
        total_ham = base_hamiltonian + disordering_ham
        rescaling_factor, shift = _rescaling_and_shift_factors(total_ham)
        rescaled_ham = (total_ham ./ rescaling_factor) + shift * I

        # eigvals returns the spectrum already sorted ascending
        eigvals_k = eigvals(Hermitian(rescaled_ham))

        bw = eigvals_k[end] - eigvals_k[1]
        if bw > 1e-14
            @inbounds for i in 1:d
                normalized_specs[i, k] = (eigvals_k[i] - eigvals_k[1]) / bw
            end
            valid_mask[k] = true
        else
            # Degenerate (e.g. disorder_strength == 0 and a symmetric base):
            # exclude from the median / argmin
            @inbounds for i in 1:d
                normalized_specs[i, k] = NaN
            end
        end

        disordering_coeffs_by_sample[k] = sample_coeffs
        next!(p)
    end

    valid_count = sum(valid_mask)
    if valid_count < 2
        error("Only $valid_count valid (non-degenerate) samples among $batch_size; " *
              "cannot compute W₂-spectral-median (increase batch_size or disorder_strength)")
    end

    # Elementwise median of the normalised sorted spectra across valid samples.
    # By bandwidth-normalisation, median_spec[1] == 0 and median_spec[end] == 1,
    # so the selection is driven by the interior d-2 levels.
    median_spec = Vector{Float64}(undef, d)
    buf = Vector{Float64}(undef, valid_count)
    @inbounds for i in 1:d
        j = 0
        for k in 1:batch_size
            if valid_mask[k]
                j += 1
                buf[j] = normalized_specs[i, k]
            end
        end
        median_spec[i] = median(buf)
    end

    # argmin L² distance to median_spec
    best_idx = 0
    best_dist = Inf
    @inbounds for k in 1:batch_size
        valid_mask[k] || continue
        dist_sq = 0.0
        for i in 1:d
            diff_i = normalized_specs[i, k] - median_spec[i]
            dist_sq += diff_i * diff_i
        end
        dist_k = sqrt(dist_sq)
        if dist_k < best_dist
            best_dist = dist_k
            best_idx = k
        end
    end
    @assert best_idx > 0  # guaranteed by valid_count >= 2

    # Reconstruct the chosen sample to obtain eigvecs (avoids storing d^2 complex
    # entries per sample during the sweep).
    best_disordering_coeffs = disordering_coeffs_by_sample[best_idx]
    disordering_ham = _construct_disordering_terms(disordering_terms, best_disordering_coeffs, num_qubits)
    total_ham = base_hamiltonian + disordering_ham
    rescaling_factor, shift = _rescaling_and_shift_factors(total_ham)
    rescaled_ham = (total_ham ./ rescaling_factor) + shift * I
    rescaled_eigvals, rescaled_eigvecs = eigen(Hermitian(rescaled_ham))
    nu_min = minimum(diff(rescaled_eigvals))

    return (
        matrix = rescaled_ham,
        terms = terms,
        base_coeffs = coeffs ./ rescaling_factor,
        disordering_terms = disordering_terms,
        disordering_coeffs = [dc ./ rescaling_factor for dc in best_disordering_coeffs],
        eigvals = rescaled_eigvals,
        eigvecs = rescaled_eigvecs,
        nu_min = nu_min,
        shift = shift,
        rescaling_factor = rescaling_factor,
        periodic = periodic,
        typicality_distance = best_dist,
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

"""
    _pad_two_site_op(term, num_qubits, q1, q2) -> SparseMatrixCSC{ComplexF64}

Place a two-site Pauli term `[op_a, op_b]` at distinct qubit indices `q1, q2`,
padding the remaining sites with identities. Used by the 2D builder where
nearest-neighbour qubit indices may not be consecutive (right-neighbour bonds
have index step `Ly`, not 1).

The returned operator is `op at q1` ⊗ `op at q2` in the natural left-to-right
qubit ordering of `kron`. For symmetric terms (`op_a == op_b`, e.g. `[X, X]`)
the placement order does not matter.
"""
function _pad_two_site_op(term::Vector{Matrix{ComplexF64}}, num_qubits::Int, q1::Int, q2::Int)
    if length(term) != 2
        throw(ArgumentError("_pad_two_site_op expects a 2-site term, got length $(length(term))"))
    end
    if q1 == q2
        throw(ArgumentError("q1 and q2 must be distinct, got q1=q2=$q1"))
    end
    if !(1 <= q1 <= num_qubits) || !(1 <= q2 <= num_qubits)
        throw(ArgumentError("q1 and q2 must be in 1:$num_qubits, got q1=$q1, q2=$q2"))
    end

    # WLOG a < b in tensor order; assign the matching operator to each side
    if q1 < q2
        a, b = q1, q2
        op_a, op_b = term[1], term[2]
    else
        a, b = q2, q1
        op_a, op_b = term[2], term[1]
    end

    id_before  = sparse(I, 2^(a - 1), 2^(a - 1))
    id_between = sparse(I, 2^(b - a - 1), 2^(b - a - 1))
    id_after   = sparse(I, 2^(num_qubits - b), 2^(num_qubits - b))

    return kron(id_before, sparse(op_a), id_between, sparse(op_b), id_after)
end

"""
    _construct_2d_heisenberg_base(Lx, Ly, terms, coeffs; periodic_x=true, periodic_y=true)
        -> Hermitian{ComplexF64, Matrix{ComplexF64}}

Build the dense base Hamiltonian for a 2D Heisenberg-family model on an
`Lx × Ly` square lattice. `terms` is `[[X,X], [Y,Y], [Z,Z]]` (or any subset)
and `coeffs` carries the uniform exchange strengths `[J_x, J_y, J_z]`.

Bond placement: for each site `(i, j)`, add bonds to its right neighbour
`(i+1, j)` (along x) and up neighbour `(i, j+1)` (along y). Periodic flags
wrap each direction. `Lx == 1` (resp. `Ly == 1`) skips x-direction (resp.
y-direction) bonds even when the periodic flag is true — wrapping a
single-cell direction would create self-bonds.

When the lattice has `Lx == 2` or `Ly == 2` with periodic BC, the wrap-around
bond coincides with the original bond and is added twice, matching the
double-counting convention of `_construct_base_ham` for 1D `n=2` chains.
"""
function _construct_2d_heisenberg_base(Lx::Int64, Ly::Int64,
    terms::Vector{Vector{Matrix{ComplexF64}}}, coeffs::Vector{Float64};
    periodic_x::Bool = true, periodic_y::Bool = true)

    if length(terms) != length(coeffs)
        throw(ArgumentError("The number of terms and coefficients must be equal"))
    end

    num_qubits = Lx * Ly
    site_index(i, j) = (i - 1) * Ly + (j - 1) + 1

    hamiltonian::SparseMatrixCSC{ComplexF64} = spzeros(2^num_qubits, 2^num_qubits)
    for (k, term) in enumerate(terms)
        for i in 1:Lx, j in 1:Ly
            # Right neighbour (x-direction): (i, j) -> (i+1, j)
            if i < Lx
                hamiltonian += coeffs[k] * _pad_two_site_op(term, num_qubits,
                    site_index(i, j), site_index(i + 1, j))
            elseif periodic_x && Lx > 1
                hamiltonian += coeffs[k] * _pad_two_site_op(term, num_qubits,
                    site_index(Lx, j), site_index(1, j))
            end

            # Up neighbour (y-direction): (i, j) -> (i, j+1)
            if j < Ly
                hamiltonian += coeffs[k] * _pad_two_site_op(term, num_qubits,
                    site_index(i, j), site_index(i, j + 1))
            elseif periodic_y && Ly > 1
                hamiltonian += coeffs[k] * _pad_two_site_op(term, num_qubits,
                    site_index(i, Ly), site_index(i, 1))
            end
        end
    end

    return Hermitian(Matrix(hamiltonian))
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
