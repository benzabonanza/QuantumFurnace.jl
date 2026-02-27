"""
    Workspace(config::Config{Lindbladian}, hamiltonian, jumps; trotter=nothing)

Construct a `Workspace{Krylov}` pre-allocating all scratch matrices for the given
(config, hamiltonian) pair. Mirrors `construct_lindbladian` setup in `furnace.jl`.

Returns `Workspace{Krylov,D,C,T,KrylovScratch{Complex{T}}}`.
"""
function Workspace(
    config::Config{Lindbladian},
    hamiltonian::HamHam,
    jumps::Vector{JumpOp};
    trotter::Union{TrottTrott, Nothing}=nothing,
)
    # Determine ham_or_trott (mirrors construct_lindbladian in furnace.jl)
    ham_or_trott = if config.domain isa TrotterDomain
        trotter === nothing && error("A Trotter object must be provided for the TrotterDomain")
        trotter
    else
        hamiltonian
    end

    # Precompute domain-specific data (transition, energy_labels, gamma_norm_factor, ...)
    precomputed_data = _precompute_data(config, ham_or_trott)

    # Precompute coherent B_total (returns nothing for GNS / with_coherent=false)
    B_total = _precompute_coherent_B(jumps, ham_or_trott, config, precomputed_data)

    # Determine dimensions and element type
    dim = size(hamiltonian.data, 1)
    T = eltype(hamiltonian.eigvals)
    CT = Complex{T}

    # Extract concrete-typed eigenbasis matrices and hermitian flags
    jump_eigenbases = [Matrix{CT}(j.in_eigenbasis) for j in jumps]
    jump_hermitian  = [j.hermitian for j in jumps]

    # Allocate KrylovScratch (no channel_rho_jump for Lindbladian)
    sc = KrylovScratch(CT, dim; with_channel_rho_jump=false)

    # Precompute G_left/G_right for optimized Lindbladian matvec (Phase 32)
    R_total = zeros(CT, dim, dim)
    _accumulate_R_total!(R_total, jump_eigenbases, jump_hermitian,
                         precomputed_data, config, ham_or_trott)
    hermitianize!(R_total)

    if B_total !== nothing
        B_T = Matrix{CT}(transpose(B_total))
        G_left  = Matrix{CT}(1im .* B_T .- 0.5 .* R_total)
        G_right = Matrix{CT}(-1im .* B_T .- 0.5 .* R_total)
    else
        G_left  = Matrix{CT}(-0.5 .* R_total)
        G_right = Matrix{CT}(-0.5 .* R_total)
    end

    G_left_adj  = G_right
    G_right_adj = G_left

    # Absorb precomputed_data fields into flat workspace fields
    pd_transition = hasproperty(precomputed_data, :transition) ? precomputed_data.transition : nothing
    pd_gnf = hasproperty(precomputed_data, :gamma_norm_factor) ? precomputed_data.gamma_norm_factor : nothing
    pd_el = hasproperty(precomputed_data, :energy_labels) ? precomputed_data.energy_labels : nothing
    pd_odp = hasproperty(precomputed_data, :oft_domain_prefactor) ? precomputed_data.oft_domain_prefactor : nothing
    pd_nufft = hasproperty(precomputed_data, :oft_nufft_prefactors) ? precomputed_data.oft_nufft_prefactors : nothing
    pd_alpha = hasproperty(precomputed_data, :alpha) ? precomputed_data.alpha : nothing
    pd_bkeys = hasproperty(precomputed_data, :bohr_keys) ? precomputed_data.bohr_keys : nothing
    pd_bis = hasproperty(precomputed_data, :bohr_is) ? precomputed_data.bohr_is : nothing
    pd_bjs = hasproperty(precomputed_data, :bohr_js) ? precomputed_data.bohr_js : nothing
    pd_bminus = hasproperty(precomputed_data, :b_minus) ? precomputed_data.b_minus : nothing
    pd_bplus = hasproperty(precomputed_data, :b_plus) ? precomputed_data.b_plus : nothing

    D = typeof(config.domain)
    C = typeof(config.construction)

    return Workspace{Krylov, D, C, T, typeof(sc)}(
        jump_eigenbases, jump_hermitian, jumps, B_total,
        G_left, G_right, G_left_adj, G_right_adj,
        nothing, nothing, nothing, nothing, nothing,  # channel fields
        pd_transition, pd_gnf, pd_el, pd_odp, pd_nufft,
        pd_alpha, pd_bkeys, pd_bis, pd_bjs, pd_bminus, pd_bplus,
        nothing,  # coherent_unitaries
        nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing,  # trajectory fields
        sc,
    )
end

# ---------------------------------------------------------------------------
# R_total accumulation helpers (physics convention: R = sum rate^2 * L' * L)
# ---------------------------------------------------------------------------

"""
    _accumulate_R_total!(R, ws, config, hamiltonian) -> nothing

Accumulate R_total = sum over all jumps and frequencies of rate^2 * (L' * L)
in physics convention. Used at workspace construction time (not per-matvec).
"""
function _accumulate_R_total!(
    R::Matrix{T},
    ws_eigenbases::Vector{Matrix{T}},
    ws_hermitian::Vector{Bool},
    precomputed_data,
    config::Config{<:Any, EnergyDomain},
    hamiltonian::HamHam,
) where {T<:Complex}
    (; transition, gamma_norm_factor, energy_labels) = precomputed_data
    bohr_freqs = hamiltonian.bohr_freqs
    inv_4sigma2 = 1.0 / (4 * config.sigma^2)
    prefactor = precomputed_data.oft_domain_prefactor * gamma_norm_factor

    dim = size(R, 1)
    jump_oft = zeros(T, dim, dim)
    LdagL = zeros(T, dim, dim)

    for (k, eigenbasis) in enumerate(ws_eigenbases)
        is_herm = ws_hermitian[k]
        if is_herm
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)
                oft!(jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)
                rate2 = prefactor * transition(w)
                mul!(LdagL, jump_oft', jump_oft)
                @. R += rate2 * LdagL
                if w > 1e-12
                    rate2_neg = prefactor * transition(-w)
                    mul!(LdagL, jump_oft, jump_oft')
                    @. R += rate2_neg * LdagL
                end
            end
        else
            for w in energy_labels
                oft!(jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)
                rate2 = prefactor * transition(w)
                mul!(LdagL, jump_oft', jump_oft)
                @. R += rate2 * LdagL
            end
        end
    end
    return nothing
end

function _accumulate_R_total!(
    R::Matrix{T},
    ws_eigenbases::Vector{Matrix{T}},
    ws_hermitian::Vector{Bool},
    precomputed_data,
    config::Config{<:Any, D},
    ham_or_trott::Union{HamHam, TrottTrott},
) where {T<:Complex, D<:Union{TimeDomain, TrotterDomain}}
    (; transition, gamma_norm_factor, energy_labels, oft_nufft_prefactors) = precomputed_data
    prefactor = precomputed_data.oft_domain_prefactor * gamma_norm_factor

    dim = size(R, 1)
    jump_oft = zeros(T, dim, dim)
    LdagL = zeros(T, dim, dim)

    for (k, eigenbasis) in enumerate(ws_eigenbases)
        is_herm = ws_hermitian[k]
        if is_herm
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)
                nufft_pf = _prefactor_view(oft_nufft_prefactors, w)
                @. jump_oft = eigenbasis * nufft_pf
                rate2 = prefactor * transition(w)
                mul!(LdagL, jump_oft', jump_oft)
                @. R += rate2 * LdagL
                if w > 1e-12
                    rate2_neg = prefactor * transition(-w)
                    mul!(LdagL, jump_oft, jump_oft')
                    @. R += rate2_neg * LdagL
                end
            end
        else
            for w in energy_labels
                nufft_pf = _prefactor_view(oft_nufft_prefactors, w)
                @. jump_oft = eigenbasis * nufft_pf
                rate2 = prefactor * transition(w)
                mul!(LdagL, jump_oft', jump_oft)
                @. R += rate2 * LdagL
            end
        end
    end
    return nothing
end

function _accumulate_R_total!(
    R::Matrix{T},
    ws_eigenbases::Vector{Matrix{T}},
    ws_hermitian::Vector{Bool},
    precomputed_data,
    config::Config{<:Any, BohrDomain},
    hamiltonian::HamHam,
) where {T<:Complex}
    (; alpha, gamma_norm_factor) = precomputed_data
    dim = size(R, 1)
    jump_oft = zeros(T, dim, dim)
    A_nu2_dag = zeros(T, dim, dim)

    for (k, eigenbasis) in enumerate(ws_eigenbases)
        for nu_2 in keys(hamiltonian.bohr_dict)
            @. jump_oft = alpha(hamiltonian.bohr_freqs, nu_2) * eigenbasis

            fill!(A_nu2_dag, 0)
            indices = hamiltonian.bohr_dict[nu_2]
            @inbounds for idx in indices
                i = idx[1]; j = idx[2]
                A_nu2_dag[j, i] = conj(eigenbasis[i, j])
            end

            mul!(R, A_nu2_dag, jump_oft, gamma_norm_factor, 1.0)
        end
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Config{Thermalize} workspace constructor
# ---------------------------------------------------------------------------

"""
    Workspace(config::Config{Thermalize}, hamiltonian, jumps; trotter=nothing)

Construct a `Workspace{Krylov}` with precomputed CPTP channel matrices for
the faithful Chen channel (Eq. 3.2).

Returns `Workspace{Krylov,D,C,T,KrylovScratch{Complex{T}}}`.
"""
function Workspace(
    config::Config{Thermalize},
    hamiltonian::HamHam,
    jumps::Vector{JumpOp};
    trotter::Union{TrottTrott, Nothing}=nothing,
)
    # Determine ham_or_trott
    ham_or_trott = if config.domain isa TrotterDomain
        trotter === nothing && error("A Trotter object must be provided for the TrotterDomain")
        trotter
    else
        hamiltonian
    end

    # Precompute domain-specific data
    precomputed_data = _precompute_data(config, ham_or_trott)

    # Precompute coherent B_total (returns nothing for GNS / with_coherent=false)
    B_total = _precompute_coherent_B(jumps, ham_or_trott, config, precomputed_data)

    # Determine dimensions and element type
    dim = size(hamiltonian.data, 1)
    T = eltype(hamiltonian.eigvals)
    CT = Complex{T}

    # Extract concrete-typed eigenbasis matrices and hermitian flags
    jump_eigenbases = [Matrix{CT}(j.in_eigenbasis) for j in jumps]
    jump_hermitian  = [j.hermitian for j in jumps]

    # Allocate KrylovScratch (with channel_rho_jump for Thermalize)
    sc = KrylovScratch(CT, dim; with_channel_rho_jump=true)

    # --- Precompute channel matrices ---
    delta = config.delta

    # 1. Compute R_total (physics convention)
    R_total = zeros(CT, dim, dim)
    _accumulate_R_total!(R_total, jump_eigenbases, jump_hermitian,
                         precomputed_data, config, ham_or_trott)
    hermitianize!(R_total)

    # Precompute G_left/G_right
    if B_total !== nothing
        B_T = Matrix{CT}(transpose(B_total))
        G_left  = Matrix{CT}(1im .* B_T .- 0.5 .* R_total)
        G_right = Matrix{CT}(-1im .* B_T .- 0.5 .* R_total)
    else
        G_left  = Matrix{CT}(-0.5 .* R_total)
        G_right = Matrix{CT}(-0.5 .* R_total)
    end

    G_left_adj  = G_right
    G_right_adj = G_left

    # 2. Compute K0, U_residual (Chen Eq. 3.2)
    channel = _build_cptp_channel(R_total, delta)

    # 3. Compute U_coherent = exp(-i*delta*B_total) if coherent
    U_coherent = if B_total !== nothing
        Matrix{CT}(exp(-1im * delta * Hermitian(B_total)))
    else
        nothing
    end

    # Absorb precomputed_data fields
    pd_transition = hasproperty(precomputed_data, :transition) ? precomputed_data.transition : nothing
    pd_gnf = hasproperty(precomputed_data, :gamma_norm_factor) ? precomputed_data.gamma_norm_factor : nothing
    pd_el = hasproperty(precomputed_data, :energy_labels) ? precomputed_data.energy_labels : nothing
    pd_odp = hasproperty(precomputed_data, :oft_domain_prefactor) ? precomputed_data.oft_domain_prefactor : nothing
    pd_nufft = hasproperty(precomputed_data, :oft_nufft_prefactors) ? precomputed_data.oft_nufft_prefactors : nothing
    pd_alpha = hasproperty(precomputed_data, :alpha) ? precomputed_data.alpha : nothing
    pd_bkeys = hasproperty(precomputed_data, :bohr_keys) ? precomputed_data.bohr_keys : nothing
    pd_bis = hasproperty(precomputed_data, :bohr_is) ? precomputed_data.bohr_is : nothing
    pd_bjs = hasproperty(precomputed_data, :bohr_js) ? precomputed_data.bohr_js : nothing
    pd_bminus = hasproperty(precomputed_data, :b_minus) ? precomputed_data.b_minus : nothing
    pd_bplus = hasproperty(precomputed_data, :b_plus) ? precomputed_data.b_plus : nothing

    D = typeof(config.domain)
    C = typeof(config.construction)

    return Workspace{Krylov, D, C, T, typeof(sc)}(
        jump_eigenbases, jump_hermitian, jumps, B_total,
        G_left, G_right, G_left_adj, G_right_adj,
        channel.K0, channel.U_residual, U_coherent, nothing, Float64(delta),
        pd_transition, pd_gnf, pd_el, pd_odp, pd_nufft,
        pd_alpha, pd_bkeys, pd_bis, pd_bjs, pd_bminus, pd_bplus,
        nothing,  # coherent_unitaries
        nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing,  # trajectory fields
        sc,
    )
end

# Backward-compatible alias
const KrylovWorkspace = Workspace
