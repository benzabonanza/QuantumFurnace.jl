"""
    AbstractTrotter{T<:AbstractFloat}

Common supertype for Trotter-cache objects passed to TrotterDomain code paths
(`construct_lindbladian`, `run_thermalize`, `B_trotter`, ...).

Two concrete subtypes:

- [`TrottTrott`](@ref): single-cache mode — one Strang Trotterization at a
  single step `t0`, one eigenbasis. Used for the dissipator-only (GNS) path
  and as the per-leg building block of `TrotterTriple`.
- [`TrotterTriple`](@ref) (qf-e4z.20): **three independent** Strang substeps and
  eigenbases — one each for the dissipative (`D`), outer coherent (`b_-`), and
  inner coherent (`b_+`) legs. The canonical KMS coherent scheme in
  TrotterDomain. The dissipator runs in `V_D` and `B_trotter` performs
  explicit inter-basis rotations.
"""
abstract type AbstractTrotter{T<:AbstractFloat} end

"""
    TrottTrott{T<:AbstractFloat}

    Single-cache Strang Trotter data for time evolution `exp(-iHt)` at a fixed
    Trotter step `t0`. Parameterized on element type `T`, inferred from the
    HamHam used for construction.

    # Fields
    - `t0`: Trotter-step duration.
    - `num_trotter_steps_per_t0`: Number of elementary Strang substeps composing one `t0` step.
    - `eigvals_t0`: Eigenvalues of `S_2(t0/M)^M`, the Strang one-step operator at duration `t0`.
    - `eigvecs`: Strang eigenbasis.
    - `bohr_freqs`: Quasi-Bohr frequencies extracted from `eigvals_t0` at scale `t0`.
"""
struct TrottTrott{T<:AbstractFloat} <: AbstractTrotter{T}
    t0::T
    num_trotter_steps_per_t0::Int
    eigvals_t0::Vector{Complex{T}}
    eigvecs::Matrix{Complex{T}}
    bohr_freqs::Matrix{T}
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
        )
end

function _trotter_bohr_freqs(trottU_T_eigvals::Vector{ComplexF64}, t::Float64)
    bohr_freqs = angle.(trottU_T_eigvals) ./ t  # quasi Bohr frequencies due to Trotterization.
    return bohr_freqs .- bohr_freqs'  # dim×dim
end

# ---------------------------------------------------------------------------
# qf-e4z.20 — independent per-leg Trotter caches (canonical KMS coherent scheme).
#
# `TrotterTriple` gives each coherent-term leg its own Trotterization,
# eigenbasis, and substep count `M_X`. `B_trotter` performs explicit
# inter-basis rotations to glue the three legs together. The three substeps
# `δt₀_X = t0_X / M_X` are fully independent — tightening one leg never
# inflates another.
# ---------------------------------------------------------------------------

"""
    TrotterTriple{T<:AbstractFloat}  <:  AbstractTrotter{T}

Three independent Strang Trotter caches (one per coherent-term leg) plus the
inter-basis rotation matrices.

# Fields
- `D::TrottTrott{T}`: dissipative leg, canonical basis `V_D`. The final
  Lindbladian / channel state `ρ` lives in `V_D`; `σ_β` rotates into `V_D`.
- `b_minus::TrottTrott{T}`: outer coherent leg (`b_-(t)` loop).
- `b_plus::TrottTrott{T}`:  inner coherent leg (`b_+(τ)` loop).
- `R_bm_in_D = V_bm' · V_D`: maps `V_D` operator coords to `V_bm` via
  `M_bm = R_bm_in_D · M_D · R_bm_in_D'`.
- `R_bp_in_D = V_bp' · V_D`: same pattern, `D → b_+`.
- `R_bm_in_bp = V_bm' · V_bp`: same pattern, `b_+ → b_-`.

# `getproperty` aliasing
Field access for the single-cache names — `:t0`, `:eigvecs`, `:bohr_freqs`,
`:eigvals_t0`, `:num_trotter_steps_per_t0` — transparently delegate to the
`D` leg. This lets every TrotterDomain consumer that reads
`trotter.eigvecs` / `trotter.bohr_freqs` / etc. work without modification
(the channel I/O basis is `V_D`; `predict_channel_trajectory` and the
qf-72g `superop_distance` tool rely on this delegation).
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

Build three independent Strang Trotter caches at substeps
`δt₀_X = t0_X / M_X` for `X ∈ {D, b_minus, b_plus}`.

`t0_b_minus_evol` and `t0_b_plus_evol` are the **Trotter-step durations**
of the outer (`b_-(t/σ)` → grid step `register_t0_b_minus / σ`) and inner
(`b_+(τβ)` → grid step `β · register_t0_b_plus`) coherent-integral
evolutions. Each leg is Trotterized independently with its OWN `M_X`, so the
three substeps need not be commensurate.

Each sub-cache is built via the single-cache `TrottTrott(ham, t0_X, M_X)`
constructor. The inter-basis rotations are computed once and stored.
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
           s === :eigvals_t0 || s === :num_trotter_steps_per_t0
        return getproperty(getfield(t, :D), s)
    else
        return getfield(t, s)  # falls through and Julia raises on unknown.
    end
end

function Base.propertynames(::TrotterTriple, private::Bool=false)
    return (:D, :b_minus, :b_plus,
            :R_bm_in_D, :R_bp_in_D, :R_bm_in_bp,
            :t0, :eigvecs, :bohr_freqs, :eigvals_t0, :num_trotter_steps_per_t0)
end

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

Verify that the 1D-chain decomposition assumed by `_trotterize2` matches the
stored Hamiltonian `ham.data`. The 1D Trotter places base and two-site
disorder bonds on consecutive qubits `(q, q+1)` (with optional wrap when
`ham.periodic`). A 2D HamHam from `find_*_2d_heisenberg` stores a different
bond structure that `_trotterize2` cannot represent.

Returns the operator-norm deviation between the 1D reconstruction and
`ham.data`. Callers that require a 1D-chain HamHam (currently the
`TrottTrott` constructors) should throw on `err > tol`.

The check costs one extra `_construct_base_ham` + `_construct_disordering_terms`
plus an `opnorm` on a `2^n × 2^n` matrix — negligible against the rest of the
Trotter construction at the sandbox n ≤ 8 envelope.
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
            "yet model 2D lattice bond structure (see qf-91g.3)."))
    end
    return err
end

function _trotterize2(hamiltonian::HamHam, t::Float64, num_trotter_steps::Int64)
    """
    2nd-order Strang Trotter for 1D chain Hamiltonians with 1- and 2-site terms.

    Honors `hamiltonian.periodic`: with `periodic=false`, every 2-site term whose
    placement would wrap past site `num_qubits` is dropped (via the `periodic`
    kwarg of `expm_pauli_padded`). 1-site terms have no boundary issue.

    **Limitation**: this is a 1D-chain Trotterizer. Passing a 2D HamHam yields
    an operator that does not reproduce `exp(-i δ H_2D)` — the bond structure of
    the 2D lattice (right + up neighbours) is not modeled here. The
    `_check_1d_trotter_compatible` guard (invoked from `TrottTrott` constructors)
    throws on this mismatch so 2D HamHams cannot reach this function silently.
    2D HamHams are currently only used in EnergyDomain/BohrDomain pipelines.
    """
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

function trotterize(hamiltonian::HamHam, T::Float64, num_trotter_steps::Int64)
    """1st order Trotter. Honors `hamiltonian.periodic`."""

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
