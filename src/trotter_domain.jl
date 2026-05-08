"""
    TrottTrott{T<:AbstractFloat}

    Stores precomputed data for Trotterized time evolution.
    Parameterized on element type `T`, inferred from the HamHam used for construction.

    # Legacy fields (single-cache mode — also set in shared-δt₀ mode where they alias the dissipative cache)
    - `t0`: Trotter-step duration of the **dissipative** register (= `t0_D`).
    - `num_trotter_steps_per_t0`: Number of elementary Strang substeps composing one `t0` step
      (= `M_D` in shared-δt₀ mode; = the user's `M` in legacy mode).
    - `eigvals_t0`: Eigenvalues of `S_2(t0/M)^M`, the Strang one-step operator at duration `t0`.
    - `eigvecs`: Strang eigenbasis. In shared-δt₀ mode this is the eigenbasis of the elementary
      `S_2(δt₀)` — the same set of eigenvectors as `eigvecs` of any power, hence shared across
      all three per-register caches.
    - `bohr_freqs`: Quasi-Bohr frequencies extracted from `eigvals_t0` at scale `t0`.

    # Per-register fields (qf-d0w shared-δt₀ scheme; `nothing` in legacy single-cache mode)
    - `eigvals_t0_b_minus`: Eigenvalues at the **outer coherent** Trotter step
      `t0_b_minus = β · register_t0_b_minus(config)`. Used by the `b_-(t)` outer loop in
      `B_trotter`.
    - `eigvals_t0_b_plus`: Eigenvalues at the **inner coherent** Trotter step
      `t0_b_plus = β · register_t0_b_plus(config)`. Used by the `b_+(τ)` inner loop in
      `B_trotter`.
    - `t0_b_minus`, `t0_b_plus`: The corresponding Trotter-step durations (β times the
      config grid spacings). Stored so the consumer can compute `num_steps =
      round(grid_index · β · t0_grid / t0_X)` and recover an integer step count
      independent of the dissipative `t0_D`.

    All three `eigvals_t0_X` derive from the same elementary `S_2(δt₀)` eigenvalues
    `λ_S`: `eigvals_t0_X = λ_S .^ M_X`. This guarantees a single shared eigenbasis
    (`eigvecs` of `λ_S`) and avoids any runtime basis alignment between the three
    coherent-term loops.
"""
struct TrottTrott{T<:AbstractFloat}
    t0::T
    num_trotter_steps_per_t0::Int
    eigvals_t0::Vector{Complex{T}}
    eigvecs::Matrix{Complex{T}}
    bohr_freqs::Matrix{T}
    # qf-d0w per-register caches (nothing in legacy single-cache mode).
    eigvals_t0_b_minus::Union{Nothing, Vector{Complex{T}}}
    eigvals_t0_b_plus::Union{Nothing, Vector{Complex{T}}}
    t0_b_minus::Union{Nothing, T}
    t0_b_plus::Union{Nothing, T}
end

"""
    TrottTrott(hamiltonian, t::Real, num_trotter_steps::Int) -> TrottTrott

Legacy single-cache constructor. Builds `S_2(t/M)^M` via `_trotterize2`, diagonalizes
it, and stores the result. The per-register fields (`eigvals_t0_b_minus`,
`eigvals_t0_b_plus`, `t0_b_minus`, `t0_b_plus`) are set to `nothing` — `B_trotter`
falls back to `eigvals_t0` for both the outer `b_-(t)` and inner `b_+(τ)` loops,
which reproduces the pre-qf-d0w behaviour byte-for-byte.
"""
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
        nothing, nothing, nothing, nothing,
        )
end

"""
    TrottTrott(hamiltonian, t0_D::Real, t0_b_minus::Real, t0_b_plus::Real, M_user::Int)
        -> TrottTrott

qf-d0w shared-δt₀ constructor for KMS coherent in TrotterDomain. Picks an elementary
Strang substep `δt₀ = min(t0_D, t0_b_minus, t0_b_plus) / M_user`, asserts that all
three natural steps are integer multiples of `δt₀`, and builds a single eigenbasis
`U_S` from `S_2(δt₀)`. The per-register `eigvals_t0_X = λ_S .^ M_X` are vector
powers — no extra Trotterizations or diagonalizations.

`t0_b_minus` and `t0_b_plus` here are the **Trotter-step durations** of the
outer and inner coherent integral evolutions:
- `b_-(t/σ)` over the outer grid → `t0_b_minus = register_t0_b_minus(config) / σ`.
- `b_+(τβ)` over the inner grid → `t0_b_plus  = β · register_t0_b_plus(config)`.

For the standard project convention `σ = 1/β` both scalings coincide. The
per-leg formulas are written separately so that `B_trotter` can recover an
**integer** step count for any σ. The dissipative step `t0_D =
register_t0_D(config)` is in raw Hamiltonian-time units and carries no
σ-rescaling.

Throws `ArgumentError` if any natural step is not an integer multiple of `δt₀`
(the integer-M condition fails for non-default `w0_D` or non-power-of-two grid
ratios; in that case the user should adjust the grid or — when implemented —
opt into the independent-cache fallback).
"""
function TrottTrott(
    hamiltonian::HamHam{T},
    t0_D::Real,
    t0_b_minus::Real,
    t0_b_plus::Real,
    M_user::Int,
) where {T<:AbstractFloat}
    M_user > 0 || throw(ArgumentError("TrottTrott shared-δt₀: M_user must be > 0 (got $M_user)."))
    t0_D > 0    || throw(ArgumentError("TrottTrott shared-δt₀: t0_D must be > 0 (got $t0_D)."))
    t0_b_minus > 0 || throw(ArgumentError("TrottTrott shared-δt₀: t0_b_minus must be > 0 (got $t0_b_minus)."))
    t0_b_plus > 0  || throw(ArgumentError("TrottTrott shared-δt₀: t0_b_plus must be > 0 (got $t0_b_plus)."))

    natural = (Float64(t0_D), Float64(t0_b_minus), Float64(t0_b_plus))
    delta_t0 = minimum(natural) / M_user
    M_D, M_bm, M_bp = _shared_delta_t0_steps(natural, delta_t0)

    # One elementary Strang one-step at δt₀ — single Trotterization, single diagonalization.
    trottU_S = _trotterize2(hamiltonian, delta_t0, 1)
    eigvals_S, eigvecs_S = eigen(trottU_S)

    eigvals_t0_D  = eigvals_S .^ M_D
    eigvals_t0_bm = eigvals_S .^ M_bm
    eigvals_t0_bp = eigvals_S .^ M_bp

    # bohr_freqs is a property of the dissipative t0 (mirrors legacy single-cache
    # semantics so existing consumers — NUFFT prefactors, OFT lifts — see the
    # same scale they always have).
    bfreqs = _trotter_bohr_freqs(eigvals_t0_D, Float64(t0_D))

    return TrottTrott{T}(
        T(t0_D),
        M_D,
        Vector{Complex{T}}(eigvals_t0_D),
        Matrix{Complex{T}}(eigvecs_S),
        Matrix{T}(bfreqs),
        Vector{Complex{T}}(eigvals_t0_bm),
        Vector{Complex{T}}(eigvals_t0_bp),
        T(t0_b_minus),
        T(t0_b_plus),
    )
end

# Helper: enforces the integer-M condition and returns (M_D, M_bm, M_bp) as Ints.
# The condition fails for non-default `w0_D` or non-power-of-two grid ratios —
# in that case the caller should adjust parameters or use a different scheme.
function _shared_delta_t0_steps(natural::NTuple{3, Float64}, delta_t0::Float64;
                                tol::Float64 = 1e-9)
    M = ntuple(i -> natural[i] / delta_t0, 3)
    M_int = ntuple(i -> round(Int, M[i]), 3)
    @inbounds for i in 1:3
        if abs(M[i] - M_int[i]) > tol * max(1.0, M[i])
            label = (i == 1) ? "t0_D" : (i == 2 ? "t0_b_minus" : "t0_b_plus")
            throw(ArgumentError(
                "TrottTrott shared-δt₀: $(label) / δt₀ = $(M[i]) is not integer " *
                "(δt₀ = $(delta_t0), tol = $tol). The shared-δt₀ scheme requires all " *
                "three natural Trotter steps to be commensurate; commonly this fails " *
                "for non-default w0_D or grid ratios that are not powers of two."))
        end
    end
    return M_int
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
    sequence_2site_not_commuting = Matrix{ComplexF64}[]
    if length(groups.noncommuting[1]) != 0
        for (term, coupling) in zip(groups.noncommuting[1], groups.noncommuting[2])
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
