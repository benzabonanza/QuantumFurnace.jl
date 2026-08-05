"""
    AbstractTrotter{T<:AbstractFloat}

Supertype for Trotter-domain caches. [`TrottTrott`](@ref) stores one Strang
cache; [`TrotterTriple`](@ref) stores independent dissipative and coherent
caches with explicit basis rotations.
"""
abstract type AbstractTrotter{T<:AbstractFloat} end

"""
    TrottTrott{T<:AbstractFloat}

Single-cache Strang Trotter data at duration `t0`.

# Fields
- `num_trotter_steps_per_t0`: elementary Strang substeps per cached step.
- `eigvals_t0`, `eigvecs`: cached step eigendecomposition.
- `bohr_freqs`: quasi-Bohr frequencies at scale `t0`.
- `source_hamiltonian`: exact Hamiltonian object used to build the cache.
"""
struct TrottTrott{T<:AbstractFloat} <: AbstractTrotter{T}
    t0::T
    num_trotter_steps_per_t0::Int
    eigvals_t0::Vector{Complex{T}}
    eigvecs::Matrix{Complex{T}}
    bohr_freqs::Matrix{T}
    source_hamiltonian::HamHam{T}
end

"""
    TrottTrott(hamiltonian, t::Real, num_trotter_steps::Int) -> TrottTrott

Single-cache constructor. Builds `S_2(t/M)^M` via `_trotterize2`, diagonalizes
it, and stores the result.
"""
function TrottTrott(hamiltonian::HamHam{T}, t::Real, num_trotter_steps::Int64) where {T<:AbstractFloat}
    _check_1d_trotter_compatible(hamiltonian)
    t_f64 = Float64(t)
    # Trotter computation always in Float64 (Pauli matrices are ComplexF64 constants).
    # Convert results to T at the end.
    trottU = _trotterize2(hamiltonian, t_f64, num_trotter_steps)
    trottU_eigvals, trottU_eigvecs = eigen(trottU)
    bfreqs = _trotter_bohr_freqs(trottU_eigvals, t_f64)  # quasi Bohr frequencies due to Trotterization.
    return TrottTrott{T}(
        T(t),
        num_trotter_steps,
        Vector{Complex{T}}(trottU_eigvals),
        Matrix{Complex{T}}(trottU_eigvecs),
        Matrix{T}(bfreqs),
        hamiltonian,
        )
end

function _trotter_bohr_freqs(trottU_T_eigvals::Vector{ComplexF64}, t::Float64)
    bohr_freqs = angle.(trottU_T_eigvals) ./ t  # quasi Bohr frequencies due to Trotterization.
    return bohr_freqs .- bohr_freqs'  # dim×dim
end

# Each channel leg has an independent Strang step and eigenbasis.

"""
    TrotterTriple{T<:AbstractFloat}  <:  AbstractTrotter{T}

Independent dissipative and coherent Strang caches with basis rotations.

# Fields
- `D`, `b_minus`, `b_plus`: per-leg [`TrottTrott`](@ref) caches.
- `R_bm_in_D`, `R_bp_in_D`, `R_bm_in_bp`: unitary coordinate rotations.

Single-cache properties such as `eigvecs` and `bohr_freqs` delegate to `D`,
which is the channel input/output basis.
"""
struct TrotterTriple{T<:AbstractFloat} <: AbstractTrotter{T}
    D::TrottTrott{T}
    b_minus::TrottTrott{T}
    b_plus::TrottTrott{T}
    R_bm_in_D::Matrix{Complex{T}}
    R_bp_in_D::Matrix{Complex{T}}
    R_bm_in_bp::Matrix{Complex{T}}
end

"""
    TrotterTriple(ham::HamHam{T}, t0_D, t0_b_minus_evol, t0_b_plus_evol,
                  M_D, M_b_minus, M_b_plus) -> TrotterTriple{T}

Build three independent Strang caches and their basis rotations.

The three `t0` values are physical evolution durations for their respective
integrals and need not be commensurate. Each `M` must be positive.
"""
function TrotterTriple(
    hamiltonian::HamHam{T},
    t0_D::Real,
    t0_b_minus_evol::Real,
    t0_b_plus_evol::Real,
    M_D::Int,
    M_b_minus::Int,
    M_b_plus::Int,
) where {T<:AbstractFloat}
    M_D       > 0 || throw(ArgumentError("TrotterTriple: M_D must be > 0 (got $M_D)."))
    M_b_minus > 0 || throw(ArgumentError("TrotterTriple: M_b_minus must be > 0 (got $M_b_minus)."))
    M_b_plus  > 0 || throw(ArgumentError("TrotterTriple: M_b_plus must be > 0 (got $M_b_plus)."))
    t0_D            > 0 || throw(ArgumentError("TrotterTriple: t0_D must be > 0 (got $t0_D)."))
    t0_b_minus_evol > 0 || throw(ArgumentError("TrotterTriple: t0_b_minus_evol must be > 0 (got $t0_b_minus_evol)."))
    t0_b_plus_evol  > 0 || throw(ArgumentError("TrotterTriple: t0_b_plus_evol must be > 0 (got $t0_b_plus_evol)."))

    D       = TrottTrott(hamiltonian, t0_D,            M_D)
    bminus  = TrottTrott(hamiltonian, t0_b_minus_evol, M_b_minus)
    bplus   = TrottTrott(hamiltonian, t0_b_plus_evol,  M_b_plus)

    # Inter-basis rotations. Convention: R_{Y←X} := V_Y' · V_X. Then for an
    # operator M expressed in V_X-coords, M_Y = R_{Y←X} · M_X · R_{Y←X}'.
    CT = Complex{T}
    R_bm_in_D   = Matrix{CT}(bminus.eigvecs' * D.eigvecs)
    R_bp_in_D   = Matrix{CT}(bplus.eigvecs'  * D.eigvecs)
    R_bm_in_bp  = Matrix{CT}(bminus.eigvecs' * bplus.eigvecs)

    return TrotterTriple{T}(D, bminus, bplus, R_bm_in_D, R_bp_in_D, R_bm_in_bp)
end

# Field-access aliasing: trotter.eigvecs / .bohr_freqs / .t0 / etc. all route
# to the D leg (the channel I/O basis). Consumers that read these single-cache
# names work unchanged on a TrotterTriple.
function Base.getproperty(t::TrotterTriple, s::Symbol)
    if s === :D || s === :b_minus || s === :b_plus ||
       s === :R_bm_in_D || s === :R_bp_in_D || s === :R_bm_in_bp
        return getfield(t, s)
    elseif s === :t0 || s === :eigvecs || s === :bohr_freqs ||
           s === :eigvals_t0 || s === :num_trotter_steps_per_t0 ||
           s === :source_hamiltonian
        return getproperty(getfield(t, :D), s)
    else
        return getfield(t, s)  # falls through and Julia raises on unknown.
    end
end

function Base.propertynames(::TrotterTriple, private::Bool=false)
    return (:D, :b_minus, :b_plus,
            :R_bm_in_D, :R_bp_in_D, :R_bm_in_bp,
            :t0, :eigvecs, :bohr_freqs, :eigvals_t0,
            :num_trotter_steps_per_t0, :source_hamiltonian)
end

"""
    compute_trotter_error(hamiltonian, trotter, t) -> Real

Return the operator-norm error between exact and cached Trotter evolution at
time `t`. `t` is interpreted as an integer number of `trotter.t0` steps.
"""
function compute_trotter_error(hamiltonian::HamHam, trotter::TrottTrott, t::Float64)

    num_t0_steps = Int(t / trotter.t0)
    exact_time_evolution = Diagonal(exp.(1im * hamiltonian.eigvals * t))  # In energy eigenbasis
    trotter_time_evolution = Diagonal(trotter.eigvals_t0.^num_t0_steps)
    trotter_time_evolution = (hamiltonian.eigvecs' * trotter.eigvecs
                                * trotter_time_evolution * trotter.eigvecs' * hamiltonian.eigvecs)
    return norm(exact_time_evolution - trotter_time_evolution)
end

"""
    _check_1d_trotter_compatible(ham; tol=1e-10)

Check that `_trotterize2`'s 1D-chain reconstruction matches `ham.data`.
Returns the operator-norm deviation and throws when it exceeds `tol`.
"""
function _check_1d_trotter_compatible(ham::HamHam{T}; tol::Real=1e-10) where {T<:AbstractFloat}
    n = Int(log2(size(ham.data, 1)))
    rescale = Float64(ham.rescaling_factor)
    # Reconstruct the 1D-chain Hamiltonian the way _trotterize2 sees it.
    H_phys = Matrix{ComplexF64}(_construct_base_ham(
        Vector{Vector{Matrix{ComplexF64}}}(ham.base_terms),
        Vector{Float64}(ham.base_coeffs) .* rescale,
        n;
        periodic=ham.periodic,
    ))
    if ham.disordering_terms !== nothing
        H_phys .+= Matrix(_construct_disordering_terms(
            Vector{Vector{Matrix{ComplexF64}}}(ham.disordering_terms),
            [Vector{Float64}(c) .* rescale for c in ham.disordering_coeffs],
            n;
            periodic=ham.periodic,
        ))
    end
    H_alg = H_phys ./ rescale .+ Float64(ham.shift) * Matrix{ComplexF64}(I, 2^n, 2^n)
    err = opnorm(H_alg .- Matrix{ComplexF64}(ham.data))
    if err > tol
        throw(ArgumentError(
            "_trotterize2 / TrottTrott expects a 1D-chain HamHam. The stored " *
            "`ham.data` deviates from `_construct_base_ham(...) + " *
            "_construct_disordering_terms(...)` (with `periodic = ham.periodic`) " *
            "by ‖ΔH‖_op = $err > tol = $tol. This usually means a 2D HamHam " *
            "from `find_*_2d_heisenberg` was passed; `_trotterize2` does not " *
            "yet model 2D lattice bond structure."))
    end
    return err
end

function _trotterize2(hamiltonian::HamHam, t::Float64, num_trotter_steps::Int64)
    # Second-order Strang Trotterization for 1D one- and two-site terms.
    # Math: $S_2(dt) = exp(A dt / 2) exp(B dt) exp(A dt / 2)$.
    # Open boundaries omit wrapping terms; constructors reject 2D bond layouts.
    timestep::Float64 = t / num_trotter_steps
    num_qubits::Int64 = Int(log2(size(hamiltonian.data)[1]))
    dim = 2^num_qubits
    odd_system::Bool = (num_qubits % 2 == 1)
    periodic::Bool = hamiltonian.periodic
    is_bdr_strange::Bool = (odd_system && periodic)

    U::Matrix{ComplexF64} = exp(im * t * Float64(hamiltonian.shift)) * I(2^num_qubits)  # Shift

    groups = group_hamiltonian_terms(hamiltonian)

    # Base terms
    odd_sites = collect(1:2:(num_qubits - 1))
    U_odd = _compute_U_group(groups.commuting[1], groups.commuting[2], odd_sites, num_qubits, timestep; periodic=periodic)

    even_sites = collect(2:2:num_qubits)
    U_even = _compute_U_group(groups.commuting[1], groups.commuting[2], even_sites, num_qubits, timestep; periodic=periodic)

    U_odd_bdr = Matrix{ComplexF64}(I, dim, dim)
    if is_bdr_strange  # Strange odd boundary — only present for PBC
        odd_bdr_site = [num_qubits]
        U_odd_bdr *= _compute_U_group(groups.commuting[1], groups.commuting[2], odd_bdr_site, num_qubits, timestep; periodic=periodic)
    end

    # 1-site terms in the Hamiltonian (with same coeffs on all sites). No BC issue.
    U_1site_terms = I(2^num_qubits)
    if length(groups.one_sites[1]) != 0
        all_sites = collect(1:num_qubits)
        U_1site_terms = _compute_U_group(groups.one_sites[1], groups.one_sites[2], all_sites, num_qubits, timestep; periodic=periodic)
    end

    # disordering part (per-site terms with different coeffs on each site, i.e. disordered)
    U_disordering = Matrix{ComplexF64}(I, dim, dim)
    if hamiltonian.disordering_terms !== nothing
        for (term, term_coeffs) in zip(hamiltonian.disordering_terms, hamiltonian.disordering_coeffs)
            term_f64 = Vector{Matrix{ComplexF64}}(term)
            for q in 1:num_qubits
                coeff_f64 = Float64(term_coeffs[q])
                expm_disordering_pauli_term = expm_pauli_padded(term_f64,
                        timestep * coeff_f64 / 2, num_qubits, q; periodic=periodic)
                U_disordering *= expm_disordering_pauli_term
            end
        end
    end

    # 2-site terms in the Hamiltonian that do not commute with e.g. XX on site (1, 2)
    sequence_2site_not_commuting = Matrix{ComplexF64}[]
    if length(groups.noncommuting[1]) != 0
        for (term, coupling) in zip(groups.noncommuting[1], groups.noncommuting[2])
            for q in 1:num_qubits
                term_f64 = Vector{Matrix{ComplexF64}}(term)
                expm_pauli_term = expm_pauli_padded(term_f64, timestep * Float64(coupling) / 2, num_qubits, q; periodic=periodic)
                push!(sequence_2site_not_commuting, expm_pauli_term)
            end
        end
    end

    # Assemble a delta step
    left_unitary_sequence = Matrix{ComplexF64}[]
    append!(left_unitary_sequence, [U_odd, U_even, U_odd_bdr, U_1site_terms, U_disordering], sequence_2site_not_commuting)
    U_step = foldl(*, left_unitary_sequence) * foldl(*, reverse(left_unitary_sequence))
    for step in 1:num_trotter_steps
        U *= U_step
    end
    return U
end

function _does_term_differ_at_both_sites(term, list_to_compare_with)::Bool

    if isempty(list_to_compare_with)
        return true
    else
        ref_term = list_to_compare_with[1]
        first_site_good::Bool = (term[1] != ref_term[1])
        second_site_good::Bool = (term[2] != ref_term[2])
        # return true if (diff, diff) or (same, same) for commutation
        return !(xor(first_site_good, second_site_good))
    end
end

"""
    group_hamiltonian_terms(hamiltonian) -> NamedTuple

Partition one- and two-site base terms for product-formula construction.

# Arguments
- `hamiltonian`: Hamiltonian with local base terms and coefficients.

# Returns
The `commuting`, `noncommuting`, and `one_sites` term/coefficient groups.
"""
function group_hamiltonian_terms(hamiltonian::HamHam{T}) where {T<:AbstractFloat}
    CT = Complex{T}
    list_of_kinda_commuting_2site_terms::Vector{Vector{Matrix{CT}}} = []
    coeffs_kinda_commuting_2site::Vector{T} = []

    list_of_not_commuting_2site_terms::Vector{Vector{Matrix{CT}}} = []
    coeffs_not_commuting_2site::Vector{T} = []

    list_of_1site_terms::Vector{Vector{Matrix{CT}}} = []
    coeffs_1site::Vector{T} = []

    for (i, term) in enumerate(hamiltonian.base_terms)
        if length(term) == 1
            push!(list_of_1site_terms, term)
            push!(coeffs_1site, hamiltonian.base_coeffs[i])
        elseif length(term) == 2
            if _does_term_differ_at_both_sites(term, list_of_kinda_commuting_2site_terms)
                push!(list_of_kinda_commuting_2site_terms, term)
                push!(coeffs_kinda_commuting_2site, hamiltonian.base_coeffs[i])
            else
                push!(list_of_not_commuting_2site_terms, term)
                push!(coeffs_not_commuting_2site, hamiltonian.base_coeffs[i])
            end
        else
            throw(ErrorException("Can only handle 1- or 2-site terms atm."))
        end
    end
    return (
        commuting = (list_of_kinda_commuting_2site_terms, coeffs_kinda_commuting_2site),
        noncommuting = (list_of_not_commuting_2site_terms, coeffs_not_commuting_2site),
        one_sites = (list_of_1site_terms, coeffs_1site)
    )
end

function _compute_U_group(terms, couplings, sites::Vector{Int64},
    num_qubits::Int64, timestep::Float64; periodic::Bool=true)::Matrix{ComplexF64}

    U_group = Matrix{ComplexF64}(I, 2^num_qubits, 2^num_qubits)
    for q in sites
        for (term, coupling) in zip(terms, couplings)
            term_f64 = Vector{Matrix{ComplexF64}}(term)
            expm_pauli_term = expm_pauli_padded(term_f64, timestep * Float64(coupling) / 2, num_qubits, q; periodic=periodic)
            U_group *= expm_pauli_term
        end
    end
    return U_group
end

"""
    trotterize(hamiltonian, T, num_trotter_steps) -> Matrix

Approximate the time-`T` propagator with a first-order product formula.

# Arguments
- `hamiltonian`: Local Hamiltonian, including boundary-condition metadata.
- `T`: Total evolution time.
- `num_trotter_steps`: Number of equal product-formula steps.

# Returns
The dense approximate propagator.
"""
function trotterize(hamiltonian::HamHam, T::Float64, num_trotter_steps::Int64)
    timestep::Float64 = T / num_trotter_steps
    num_qubits::Int64 = Int(log2(size(hamiltonian.data)[1]))
    periodic::Bool = hamiltonian.periodic

    U::Matrix{ComplexF64} = exp(im * T * Float64(hamiltonian.shift)) * I(2^num_qubits)  # Shift
    p = Progress(num_trotter_steps)
    @showprogress dt=1 desc="Trotterizing (1st order)..." for step in 1:num_trotter_steps
        # Base Hamiltonian
        for q in 1:num_qubits
            for (i, term) in enumerate(hamiltonian.base_terms)
                    term_f64 = Vector{Matrix{ComplexF64}}(term)
                    expm_pauli_term = expm_pauli_padded(term_f64, timestep * Float64(hamiltonian.base_coeffs[i]), num_qubits, q; periodic=periodic)
                    U *= expm_pauli_term
            end

        # disordering
            if hamiltonian.disordering_terms !== nothing
                for (dis_term, dis_coeffs) in zip(hamiltonian.disordering_terms, hamiltonian.disordering_coeffs)
                    dis_term_f64 = Vector{Matrix{ComplexF64}}(dis_term)
                    expm_disordering_pauli_term = expm_pauli_padded(dis_term_f64,
                                                                timestep * Float64(dis_coeffs[q]),
                                                                num_qubits, q; periodic=periodic)
                    U *= expm_disordering_pauli_term
                end
            end
        end
    end
    return U
end
