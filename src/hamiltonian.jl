"""
    HamHam{T<:AbstractFloat}

Hamiltonian, spectral data, Bohr-frequency groups, and Gibbs state.

Real fields use `T`; matrix fields use `Complex{T}`. `data` and `eigvals`
contain the spectrum rescaled to `[0, 0.45]`, while `rescaling_factor` converts
between physical and algorithm units. `gibbs` is stored in the eigenbasis.
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
    _gibbs_in_eigen(eigvals, beta) -> Matrix

Return the diagonal eigenbasis state `\$rho_i = exp(-beta E_i) / Z\$`.
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
    # Downward precision conversion is rejected; upward promotion is allowed.
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
    # Downward precision conversion is rejected; upward promotion is allowed.
    if T !== Float64 && T <: Union{Float16, Float32}
        throw(ArgumentError(
            "Expected $(Complex{T}) term data, got ComplexF64. " *
            "Reconstruct with $(Complex{T}) inputs or use default Float64 precision."))
    end

    if length(disordering_terms) != length(disordering_coeffs)
        throw(ArgumentError("Number of disordering terms must match number of coefficient vectors"))
    end

    base_hamiltonian = _construct_base_ham(terms, coeffs, num_qubits; periodic=periodic)
    disordering_hamiltonian = _construct_disordering_terms(disordering_terms, disordering_coeffs, num_qubits;
        periodic=periodic)
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

# Wrap a single disorder term in the multi-term representation.
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

Construct a fully initialised Hamiltonian from builder output.

# Arguments
- `raw`: Named tuple returned by `build_heis_1d` or `build_tfim_2d`.
- `beta`: Algorithm-side inverse temperature for the rescaled spectrum.

# Returns
A `HamHam` with derived Bohr data and Gibbs state.
"""
function HamHam(raw::NamedTuple, beta::Real)
    T = eltype(raw.eigvals)
    beta_T = T(beta)
    bohr_freqs = raw.eigvals .- transpose(raw.eigvals)
    bohr_dict = create_bohr_dict(bohr_freqs)
    gibbs = Hermitian(_gibbs_in_eigen(raw.eigvals, beta_T))

    # Accept scalar and multi-term disorder schemas.
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

Construct a Hamiltonian at physical inverse temperature `beta_phys`.

The constructor sets `\$beta_alg = beta_phys dot raw.rescaling_factor\$` before
forming the Gibbs state. The positional constructor instead accepts `beta_alg`.
"""
function HamHam(raw::NamedTuple; beta_phys::Real)
    rescale = raw.rescaling_factor
    return HamHam(raw, beta_phys * rescale)
end

"""
    beta_alg(ham::HamHam{T}, beta_phys::Real)   where T<:AbstractFloat -> T
    beta_phys(ham::HamHam{T}, beta_alg::Real)   where T<:AbstractFloat -> T

Convert inverse temperature between the physical and rescaled Hamiltonians.

The convention is `\$beta_alg = beta_phys dot ham.rescaling_factor\$`. Use these
helpers whenever a sweep crosses the physical/algorithm boundary.
"""
beta_alg(ham::HamHam{T}, beta_phys::Real) where {T<:AbstractFloat} = T(beta_phys * ham.rescaling_factor)
beta_phys(ham::HamHam{T}, beta_alg::Real) where {T<:AbstractFloat} = T(beta_alg / ham.rescaling_factor)

"""
    _unpack_disordering_fields(raw::NamedTuple, T) -> (terms, coeffs)

Return typed disorder terms and coefficients, or `(nothing, nothing)`.
"""
function _unpack_disordering_fields(raw::NamedTuple, ::Type{T}) where {T}
    if !haskey(raw, :disordering_terms) || raw.disordering_terms === nothing
        return (nothing, nothing)
    end
    terms = [Vector{Matrix{Complex{T}}}(t) for t in raw.disordering_terms]
    coeffs = [Vector{T}(c) for c in raw.disordering_coeffs]
    return (terms, coeffs)
end


"""
    build_heis_1d(num_qubits, coeffs; seed, periodic=true,
        disordering_terms=[[Z], [Z,Z]], disorder_strength=0.1) -> NamedTuple

Build one reproducible disordered 1D Heisenberg Hamiltonian.

# Arguments
- `num_qubits`: chain length.
- `coeffs`: `[J_x, J_y, J_z]` uniform exchange couplings.

# Keywords
- `seed`: MersenneTwister seed for the per-site disorder draw.
- `periodic`: Apply periodic boundaries to the base and two-site disorder.
- `disordering_terms`: Pauli terms for disorder, e.g. `[[Z], [Z, Z]]` (default).
- `disorder_strength`: per-coefficient amplitude `c ∈ [0, disorder_strength)`.

# Returns
A named tuple accepted by either `HamHam` constructor.
"""
function build_heis_1d(num_qubits::Int, coeffs::Vector{Float64};
        seed::Int,
        periodic::Bool=true,
        disordering_terms::Vector{Vector{Matrix{ComplexF64}}}=Vector{Matrix{ComplexF64}}[[Z], [Z, Z]],
        disorder_strength::Float64=0.1)

    base_terms = Vector{Matrix{ComplexF64}}[[X, X], [Y, Y], [Z, Z]]
    base_hamiltonian = _construct_base_ham(base_terms, coeffs, num_qubits; periodic=periodic)

    rng = MersenneTwister(seed)
    sample_coeffs = [zeros(Float64, num_qubits) for _ in disordering_terms]
    for dc in sample_coeffs
        rand!(rng, dc)
        dc .*= disorder_strength
    end
    disordering_ham = _construct_disordering_terms(disordering_terms, sample_coeffs, num_qubits;
        periodic=periodic)

    total_ham = Hermitian(Matrix(base_hamiltonian) + Matrix(disordering_ham))
    rescaling_factor, shift = _rescaling_and_shift_factors(total_ham)
    rescaled_ham = (Matrix(total_ham) ./ rescaling_factor) + shift * I(2^num_qubits)
    rescaled_eigvals, rescaled_eigvecs = eigen(Hermitian(rescaled_ham))
    nu_min = minimum(diff(rescaled_eigvals))

    return (
        matrix = rescaled_ham,
        terms = base_terms,
        base_coeffs = coeffs ./ rescaling_factor,
        disordering_terms = disordering_terms,
        disordering_coeffs = [dc ./ rescaling_factor for dc in sample_coeffs],
        eigvals = rescaled_eigvals,
        eigvecs = rescaled_eigvecs,
        nu_min = nu_min,
        shift = shift,
        rescaling_factor = rescaling_factor,
        periodic = periodic,
        seed = seed,
        disorder_strength = disorder_strength,
    )
end

"""
    build_tfim_2d(Lx, Ly; J=1.0, h=1.0, seed, periodic_x=true, periodic_y=true,
        disordering_terms=[[Z], [Z,Z]], disorder_strength=1e-3) -> NamedTuple

Build one reproducible disordered 2D transverse-field Ising Hamiltonian.

The model is `\$H = -J sum_(angle.l i,j angle.r) Z_i Z_j - h sum_i X_i + H_dis\$`.
Two-site disorder follows nearest-neighbour lattice bonds.

# Arguments
- `Lx`, `Ly`: Lattice dimensions.

# Keywords
- `J`: Nearest-neighbour `ZZ` coupling.
- `h`: Transverse-field magnitude.
- `seed`: MersenneTwister seed.
- `periodic_x`, `periodic_y`: per-direction BCs.
- `disordering_terms`, `disorder_strength`: as in [`build_heis_1d`], placed on
   the 2D lattice via the 2D builder.

# Returns
A named tuple accepted by either `HamHam` constructor.

Zero disorder can retain exact symmetries and make gap computations
ill-conditioned; the default weak disorder breaks those degeneracies.
"""
function build_tfim_2d(Lx::Int, Ly::Int;
        J::Float64=1.0, h::Float64=1.0,
        seed::Int,
        periodic_x::Bool=true, periodic_y::Bool=true,
        disordering_terms::Vector{Vector{Matrix{ComplexF64}}}=Vector{Matrix{ComplexF64}}[[Z], [Z, Z]],
        disorder_strength::Float64=1e-3)

    if Lx < 1 || Ly < 1
        throw(ArgumentError("Lx and Ly must be at least 1; got Lx=$Lx, Ly=$Ly"))
    end

    num_qubits = Lx * Ly
    H_bond = _construct_2d_heisenberg_base(Lx, Ly,
        Vector{Matrix{ComplexF64}}[[Z, Z]], [-J];
        periodic_x=periodic_x, periodic_y=periodic_y)
    # Math: the transverse-field term is $-h sum_i X_i$.
    field_coeffs = fill(-h, num_qubits)
    H_field = _construct_disordering_terms(
        Vector{Matrix{ComplexF64}}[[X]], [field_coeffs], num_qubits)
    base_clean = Hermitian(Matrix(H_bond) + Matrix(H_field))

    rng = MersenneTwister(seed)
    sample_coeffs = [zeros(Float64, num_qubits) for _ in disordering_terms]
    for dc in sample_coeffs
        rand!(rng, dc)
        dc .*= disorder_strength
    end
    disordering_ham = _construct_disordering_terms_2d(Lx, Ly,
        disordering_terms, sample_coeffs;
        periodic_x=periodic_x, periodic_y=periodic_y)

    total_ham = Hermitian(Matrix(base_clean) + Matrix(disordering_ham))
    rescaling_factor, shift = _rescaling_and_shift_factors(total_ham)
    rescaled_ham = (Matrix(total_ham) ./ rescaling_factor) + shift * I(2^num_qubits)
    rescaled_eigvals, rescaled_eigvecs = eigen(Hermitian(rescaled_ham))
    nu_min = minimum(diff(rescaled_eigvals))

    return (
        matrix = rescaled_ham,
        terms = Vector{Matrix{ComplexF64}}[[Z, Z], [X]],
        base_coeffs = [-J / rescaling_factor, -h / rescaling_factor],
        disordering_terms = disordering_terms,
        disordering_coeffs = [dc ./ rescaling_factor for dc in sample_coeffs],
        eigvals = rescaled_eigvals,
        eigvecs = rescaled_eigvecs,
        nu_min = nu_min,
        shift = shift,
        rescaling_factor = rescaling_factor,
        periodic = periodic_x && periodic_y,
        seed = seed,
        disorder_strength = disorder_strength,
        Lx = Lx, Ly = Ly,
        J = J, h = h,
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

Place a two-site operator at distinct qubit indices and pad with identities.

The returned sparse matrix follows the left-to-right tensor order of `kron`.
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

    # Preserve operator order after sorting sites into tensor order.
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

Build a dense nearest-neighbour Hamiltonian on an `Lx × Ly` square lattice.

Each site contributes right and upward bonds. Periodic directions wrap unless
their length is one; length-two periodic directions retain the package's
double-bond convention.
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

"""
    _construct_disordering_terms(terms, coeffs, num_qubits; periodic=true)
        -> Hermitian{ComplexF64, Matrix{ComplexF64}}

Build site and nearest-neighbour disorder on a 1D chain.

`periodic` controls the final two-site wrap bond. Use
`_construct_disordering_terms_2d` for lattice bonds in two dimensions.
"""
function _construct_disordering_terms(terms::Vector{Vector{Matrix{ComplexF64}}},
    coeffs::Vector{Vector{Float64}}, num_qubits::Int64; periodic::Bool=true)

    disordering_hamiltonian::SparseMatrixCSC{ComplexF64} = spzeros(2^num_qubits, 2^num_qubits)
    for (term, term_coeffs) in zip(terms, coeffs)
        if length(term_coeffs) != num_qubits
            throw(ArgumentError("Each disordering coefficient vector must have length num_qubits ($num_qubits), got $(length(term_coeffs))"))
        end
        for q in 1:num_qubits
            disordering_hamiltonian += term_coeffs[q] * pad_term(term, num_qubits, q; periodic=periodic)
        end
    end

    return Hermitian(Matrix(disordering_hamiltonian))
end

"""
    _construct_disordering_terms_2d(Lx, Ly, terms, coeffs; periodic_x=true, periodic_y=true)
        -> Hermitian{ComplexF64, Matrix{ComplexF64}}

Build site and nearest-neighbour disorder on an `Lx × Ly` lattice.

Each site's two-site coefficient is shared by its right and upward bonds, so
the two bond directions are correlated. Terms must act on one or two sites.
"""
function _construct_disordering_terms_2d(Lx::Int64, Ly::Int64,
    terms::Vector{Vector{Matrix{ComplexF64}}},
    coeffs::Vector{Vector{Float64}};
    periodic_x::Bool=true, periodic_y::Bool=true)

    num_qubits = Lx * Ly
    site_index(i, j) = (i - 1) * Ly + (j - 1) + 1

    disordering_hamiltonian::SparseMatrixCSC{ComplexF64} = spzeros(2^num_qubits, 2^num_qubits)
    for (term, term_coeffs) in zip(terms, coeffs)
        if length(term_coeffs) != num_qubits
            throw(ArgumentError("Each disordering coefficient vector must have length num_qubits ($num_qubits), got $(length(term_coeffs))"))
        end
        if length(term) == 1
            # Site-local term.
            for i in 1:Lx, j in 1:Ly
                q = site_index(i, j)
                disordering_hamiltonian += term_coeffs[q] * pad_term(term, num_qubits, q; periodic=true)
            end
        elseif length(term) == 2
            # Share each site coefficient between its right and upward bonds.
            for i in 1:Lx, j in 1:Ly
                c = term_coeffs[site_index(i, j)]
                # Right neighbour (x-direction): (i, j) -> (i+1, j)
                if i < Lx
                    disordering_hamiltonian += c * _pad_two_site_op(term, num_qubits,
                        site_index(i, j), site_index(i + 1, j))
                elseif periodic_x && Lx > 1
                    disordering_hamiltonian += c * _pad_two_site_op(term, num_qubits,
                        site_index(Lx, j), site_index(1, j))
                end
                # Up neighbour (y-direction): (i, j) -> (i, j+1)
                if j < Ly
                    disordering_hamiltonian += c * _pad_two_site_op(term, num_qubits,
                        site_index(i, j), site_index(i, j + 1))
                elseif periodic_y && Ly > 1
                    disordering_hamiltonian += c * _pad_two_site_op(term, num_qubits,
                        site_index(i, Ly), site_index(i, 1))
                end
            end
        else
            throw(ArgumentError(
                "_construct_disordering_terms_2d only supports 1- or 2-site terms, " *
                "got length-$(length(term))"))
        end
    end

    return Hermitian(Matrix(disordering_hamiltonian))
end

"""Return affine factors that map a Hamiltonian spectrum to `[0, 0.45]`."""
function _rescaling_and_shift_factors(hamiltonian::Hermitian)

    eps = 0.1  # Keep the upper endpoint below the algorithmic wrap at 0.5.
    eigenergies = eigvals(hamiltonian)
    smallest_eigval = minimum(eigenergies)
    largest_eigval = maximum(eigenergies)

    rescaling_factor = (largest_eigval - smallest_eigval) * (2 / (1 - eps))
    shift = - (largest_eigval - smallest_eigval * eps) / (2 * (largest_eigval - smallest_eigval)) + 0.5
    return rescaling_factor, shift
end
