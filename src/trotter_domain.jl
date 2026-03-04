"""
    TrottTrott{T<:AbstractFloat}

    Stores precomputed data for Trotterized time evolution.
    Parameterized on element type `T`, inferred from the HamHam used for construction.

    # Fields
    - `t0`: The time unit for the Trotter step.
    - `num_trotter_steps_per_t0`: Self-explanatory. Usually `t0` is small enough to just use 1 Trotter step for it.
    - `eigvals_t0`, `eigvecs`: Eigenvalues of the evolution operator for one time unit `t0`, and corresponding eigenvectors.
"""
struct TrottTrott{T<:AbstractFloat}
    t0::T
    num_trotter_steps_per_t0::Int
    eigvals_t0::Vector{Complex{T}}
    eigvecs::Matrix{Complex{T}}
    bohr_freqs::Matrix{T}
end

function TrottTrott(hamiltonian::HamHam{T}, t::Real, num_trotter_steps::Int64) where {T<:AbstractFloat}
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
        )
end

function _trotter_bohr_freqs(trottU_T_eigvals::Vector{ComplexF64}, t::Float64)
    bohr_freqs = angle.(trottU_T_eigvals) ./ t  # quasi Bohr frequencies due to Trotterization.
    return bohr_freqs .- bohr_freqs'  # dim×dim
end

function compute_trotter_error(hamiltonian::HamHam, trotter::TrottTrott, t::Float64)

    num_t0_steps = Int(t / trotter.t0)
    exact_time_evolution = Diagonal(exp.(1im * hamiltonian.eigvals * t))  # In energy eigenbasis
    trotter_time_evolution = Diagonal(trotter.eigvals_t0.^num_t0_steps)
    trotter_time_evolution = (hamiltonian.eigvecs' * trotter.eigvecs
                                * trotter_time_evolution * trotter.eigvecs' * hamiltonian.eigvecs)
    return norm(exact_time_evolution - trotter_time_evolution)
end

function _trotterize2(hamiltonian::HamHam, t::Float64, num_trotter_steps::Int64)
    """For 1 and 2 site Hamiltonians"""
    timestep::Float64 = t / num_trotter_steps
    num_qubits::Int64 = Int(log2(size(hamiltonian.data)[1]))
    dim = 2^num_qubits
    odd_system::Bool = (num_qubits % 2 == 1)
    is_bdr_strange::Bool = (odd_system && hamiltonian.periodic)

    U::Matrix{ComplexF64} = exp(im * t * Float64(hamiltonian.shift)) * I(2^num_qubits)  # Shift

    groups = group_hamiltonian_terms(hamiltonian)

    # Base terms
    odd_sites = collect(1:2:(num_qubits - 1))
    U_odd = _compute_U_group(groups.commuting[1], groups.commuting[2], odd_sites, num_qubits, timestep)

    even_sites = collect(2:2:num_qubits)
    U_even = _compute_U_group(groups.commuting[1], groups.commuting[2], even_sites, num_qubits, timestep)

    U_odd_bdr = Matrix{ComplexF64}(I, dim, dim)
    if is_bdr_strange  # Strange odd boundary
        odd_bdr_site = [num_qubits]
        U_odd_bdr *= _compute_U_group(groups.commuting[1], groups.commuting[2], odd_bdr_site, num_qubits, timestep)
    end

    # 1-site terms in the Hamiltonian (with same coeffs on all sites)
    U_1site_terms = I(2^num_qubits)
    if length(groups.one_sites[1]) != 0
        all_sites = collect(1:num_qubits)
        U_1site_terms = _compute_U_group(groups.one_sites[1], groups.one_sites[2], all_sites, num_qubits, timestep)
    end

    # disordering part (per-site terms with different coeffs on each site, i.e. disordered)
    U_disordering = Matrix{ComplexF64}(I, dim, dim)
    if hamiltonian.disordering_terms !== nothing
        for (term, term_coeffs) in zip(hamiltonian.disordering_terms, hamiltonian.disordering_coeffs)
            term_f64 = Vector{Matrix{ComplexF64}}(term)
            for q in 1:num_qubits
                coeff_f64 = Float64(term_coeffs[q])
                expm_disordering_pauli_term = expm_pauli_padded(term_f64,
                        timestep * coeff_f64 / 2, num_qubits, q)
                U_disordering *= expm_disordering_pauli_term
            end
        end
    end

    # 2-site terms in the Hamiltonian that do not commute with e.g. XX on site (1, 2)
    sequence_2site_not_commuting = []
    if length(groups.noncommuting[1]) != 0
        for (term, coupling) in (groups.noncommuting[1], groups.noncommuting[2])
            for q in 1:num_qubits
                term_f64 = Vector{Matrix{ComplexF64}}(term)
                expm_pauli_term = expm_pauli_padded(term_f64, timestep * Float64(coupling) / 2, num_qubits, q)
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
    num_qubits::Int64, timestep::Float64)::Matrix{ComplexF64}

    U_group = Matrix{ComplexF64}(I, 2^num_qubits, 2^num_qubits)
    for q in sites
        for (term, coupling) in zip(terms, couplings)
            term_f64 = Vector{Matrix{ComplexF64}}(term)
            expm_pauli_term = expm_pauli_padded(term_f64, timestep * Float64(coupling) / 2, num_qubits, q)
            U_group *= expm_pauli_term
        end
    end
    return U_group
end

function trotterize(hamiltonian::HamHam, T::Float64, num_trotter_steps::Int64)
    """1st order Trotter, periodic"""

    timestep::Float64 = T / num_trotter_steps
    num_qubits::Int64 = Int(log2(size(hamiltonian.data)[1]))

    U::Matrix{ComplexF64} = exp(im * T * Float64(hamiltonian.shift)) * I(2^num_qubits)  # Shift
    p = Progress(num_trotter_steps)
    @showprogress dt=1 desc="Trotterizing (1st order)..." for step in 1:num_trotter_steps
        # Base Hamiltonian
        for q in 1:num_qubits
            for (i, term) in enumerate(hamiltonian.base_terms)
                    term_f64 = Vector{Matrix{ComplexF64}}(term)
                    expm_pauli_term = expm_pauli_padded(term_f64, timestep * Float64(hamiltonian.base_coeffs[i]), num_qubits, q)
                    U *= expm_pauli_term
            end

        # disordering
            if hamiltonian.disordering_terms !== nothing
                for (dis_term, dis_coeffs) in zip(hamiltonian.disordering_terms, hamiltonian.disordering_coeffs)
                    dis_term_f64 = Vector{Matrix{ComplexF64}}(dis_term)
                    expm_disordering_pauli_term = expm_pauli_padded(dis_term_f64,
                                                                timestep * Float64(dis_coeffs[q]),
                                                                num_qubits, q)
                    U *= expm_disordering_pauli_term
                end
            end
        end
    end
    return U
end
