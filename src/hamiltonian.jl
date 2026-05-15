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
[`build_heis_1d`] / [`build_tfim_2d`]) plus inverse temperature `beta`.

Infers T from `eltype(raw.eigvals)`. Computes `bohr_freqs`, `bohr_dict`, and `gibbs`
from the raw eigvals. Extra NamedTuple fields beyond the canonical 11 (e.g.
`seed`, `disorder_strength`, `Lx`, `Ly`, `J`, `h` from the builders) are silently ignored.
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
[`build_heis_1d`] / [`build_tfim_2d`] to map the physical spectrum into `[0, 0.45]`. For
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

Unpack the multi-term `disordering_terms` / `disordering_coeffs` NamedTuple
fields produced by [`build_heis_1d`] / [`build_tfim_2d`] into the
`Vector{Vector{Matrix{Complex{T}}}}` / `Vector{Vector{T}}` shape that the
HamHam struct stores. Returns `(nothing, nothing)` when the NamedTuple has no
`disordering_terms` key or when the field is `nothing`.
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

Build ONE disordered 1D Heisenberg fixture at the given random seed. No batch,
no spectral-typicality selector — the disorder draw is fully reproducible
from `seed`.

Produces a NamedTuple with fields
`matrix, terms, base_coeffs, disordering_terms, disordering_coeffs, eigvals,
eigvecs, nu_min, shift, rescaling_factor, periodic, seed, disorder_strength`
so that `HamHam(raw, β)` / `HamHam(raw; beta_phys=β)` accept it unchanged
(extra fields beyond the canonical 11 are silently ignored).

# Methodology
Downstream analyses that want a disorder average should generate N
realisations across N seeds, run the observable of interest at each, and
report median + IQR / min-max band. The qf-yi4 decision (2026-05-15)
deprecated the prior find_typical / find_ideal spectral-selector
approach because its "typicality" criterion (L² to bandwidth-normalised
median spectrum) has no physical justification for any specific
observable.

# Arguments
- `num_qubits`: chain length.
- `coeffs`: `[J_x, J_y, J_z]` uniform exchange couplings.

# Keywords
- `seed`: MersenneTwister seed for the per-site disorder draw.
- `periodic`: 1D-chain BCs for both base and 2-site disorder.
- `disordering_terms`: Pauli terms for disorder, e.g. `[[Z], [Z, Z]]` (default).
- `disorder_strength`: per-coefficient amplitude `c ∈ [0, disorder_strength)`.
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

Build ONE 2D transverse-field Ising fixture
    H = −J Σ_{<i,j>} Z_i Z_j − h Σ_i X_i + ε·(per-site Z + per-bond ZZ)
at the given seed. Two-site disorder rides actual nearest-neighbour
lattice bonds via [`_construct_disordering_terms_2d`] (right + up
neighbour per site), respecting `periodic_x` / `periodic_y` independently.

Schema matches [`build_heis_1d`] plus `Lx, Ly, J, h` diagnostic fields.
`HamHam(raw, β)` / `HamHam(raw; beta_phys=β)` accept it unchanged.

# Operating-point references (see memory `tc_2d_tfim_phase_diagram.md`)
- Disordered: `h > h_c ≈ 3.044` puts the ground state in the
  paramagnetic phase at any temperature (e.g. `h = 3.5`).
- Ordered:    `h < h_c`, `β_phys > 1/T_c(h)` puts the system in the
  symmetry-broken phase (e.g. `h = 1.0`, `β_phys ≥ 2.0`, `T_c(1) ≈ 2.07`).

# Keywords
- `J`: ZZ exchange coupling (negative-sign convention: −J Σ ZZ).
- `h`: transverse field magnitude (negative-sign convention: −h Σ X).
- `seed`: MersenneTwister seed.
- `periodic_x`, `periodic_y`: per-direction BCs.
- `disordering_terms`, `disorder_strength`: as in [`build_heis_1d`], placed on
   the 2D lattice via the 2D builder.
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
    # Transverse field as a uniform per-site X coefficient (single-site, no BC issue).
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

"""
    _construct_disordering_terms(terms, coeffs, num_qubits; periodic=true)
        -> Hermitian{ComplexF64, Matrix{ComplexF64}}

Build the per-site disordering Hamiltonian on a 1D chain of `num_qubits` sites
from a list of Pauli terms `terms` and per-site coefficient vectors `coeffs`.
The `periodic` flag controls whether two-site terms wrap around the chain
boundary; `pad_term(...; periodic=false)` returns zero for wrap-around
placements, so OBC fixtures get an OBC disordering Hamiltonian.

# Semantics
- Single-site terms (e.g. `[Z]`): one Pauli per site, no boundary issue.
- Two-site terms (e.g. `[Z, Z]`): bond `(q, q+1)` along a 1D chain.
  - `periodic=true` (default): adds the wrap-around bond `(num_qubits, 1)`.
  - `periodic=false`: drops the wrap-around bond.

# 2D caveat
This is a **1D-chain builder**. Calling it directly with a 2D HamHam and
two-site disorder places the disorder on 1D-chain bonds `(q, q+1)` on
the linearised site index — NOT on the actual 2D lattice bonds. The 2D
fixture builder [`build_tfim_2d`] uses [`_construct_disordering_terms_2d`]
instead. Single-site disorder is unaffected (site-local).
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

2D-lattice version of [`_construct_disordering_terms`]. Places single-site
terms at each of the `num_qubits = Lx * Ly` sites and two-site terms on the
nearest-neighbour bonds of the `Lx × Ly` square lattice (right neighbour
along x, up neighbour along y; site-to-qubit map matches
`_construct_2d_heisenberg_base`).

Per-term `coeffs` is a length-`num_qubits` vector: site `(i, j)` carries
`coeffs[k][site_index(i, j)]`. For two-site terms, that coefficient is
applied to BOTH the right-neighbour bond (i, j)→(i+1, j) and the
up-neighbour bond (i, j)→(i, j+1) that emanate from `(i, j)`. Wrap-around
bonds in each direction are included only when the corresponding periodic
flag is true and the lattice has length > 1 in that direction.

**Disorder correlation note**: the x- and y-bond disorder coefficients
sharing the same per-site random number means the disorder is *correlated*
across the two bond directions, not i.i.d. per-bond. This is appropriate
for ε ≈ 1e-3 symmetry-breaking (the only requirement is that the
coefficients are generically non-zero so every relevant symmetry is
broken). For strong-disorder MBL-style fixtures where bond-level
independence is expected, use the per-bond builder
[`_construct_disordering_terms_2d_per_bond`] (not implemented — file a
beads issue when the use case arises).

Throws `ArgumentError` for term lengths other than 1 or 2.
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
            # Single-site: site-local Pauli at each lattice site
            for i in 1:Lx, j in 1:Ly
                q = site_index(i, j)
                disordering_hamiltonian += term_coeffs[q] * pad_term(term, num_qubits, q; periodic=true)
            end
        elseif length(term) == 2
            # Two-site: place on actual nearest-neighbour lattice bonds, mirroring
            # _construct_2d_heisenberg_base. Per-site coeff is applied to both the
            # site's right-bond and up-bond.
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
