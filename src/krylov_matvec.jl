# ---------------------------------------------------------------------------
# Transition wrapper (Phase 35: zero-allocation function barrier)
# ---------------------------------------------------------------------------

"""
    _TransitionWrap{F}

Immutable wrapper around the transition function. Wrapping `transition` in a
parameterized struct forces Julia to dispatch on the concrete closure type `F`
at the function barrier call boundary, avoiding the 80-byte Union boxing that
occurs when passing `Union{Nothing, Function}` fields directly.

Created on-the-fly in the `apply_lindbladian!` / `apply_adjoint_lindbladian!`
barrier functions; never stored persistently.
"""
struct _TransitionWrap{F}
    f::F
end

# ---------------------------------------------------------------------------
# Sandwich-only helpers (Phase 32 optimization)
# These are stripped-down versions of the dissipator functions with the
# L'L and anticommutator terms removed (absorbed into precomputed G_left/G_right).
# Each performs only 2 GEMMs per call (down from 5 in the full dissipator).
# ---------------------------------------------------------------------------

"""
    _accumulate_sandwich!(out, L_op, rho, scalar, ws) -> nothing

Accumulate `scalar * L * rho * L'` into `out`. Used for:
- Forward Lindbladian positive-frequency sandwich (L * rho * L')
- Adjoint Lindbladian negative-frequency Hermitian partner (HS adjoint of L'*rho*L is L*rho*L')

Uses `ws.scratch.sandwich_tmp`, `ws.scratch.sandwich_out` as scratch.
"""
@inline function _accumulate_sandwich!(
    out::Matrix{T},
    L_op::Matrix{T},
    rho::Matrix{T},
    scalar::Real,
    ws::Workspace{Krylov},
) where {T<:Complex}
    _accumulate_sandwich_scratch!(out, L_op, rho, scalar, ws.scratch.sandwich_tmp, ws.scratch.sandwich_out)
end

"""
    _accumulate_sandwich_adj!(out, L_op, rho, scalar, ws) -> nothing

Accumulate `scalar * L' * rho * L` into `out`. Used for:
- Forward Lindbladian negative-frequency Hermitian partner (L_neg = L')
- Adjoint Lindbladian positive-frequency sandwich (HS adjoint of L*rho*L' is L'*rho*L)

Uses `ws.scratch.sandwich_tmp`, `ws.scratch.sandwich_out` as scratch.
"""
@inline function _accumulate_sandwich_adj!(
    out::Matrix{T},
    L_op::Matrix{T},
    rho::Matrix{T},
    scalar::Real,
    ws::Workspace{Krylov},
) where {T<:Complex}
    _accumulate_sandwich_adj_scratch!(out, L_op, rho, scalar, ws.scratch.sandwich_tmp, ws.scratch.sandwich_out)
end

# --- Scratch-only sandwich helpers (no Workspace access, zero allocation on hot path) ---

@inline function _accumulate_sandwich_scratch!(
    out::Matrix{T},
    L_op::Matrix{T},
    rho::Matrix{T},
    scalar::Real,
    sandwich_tmp::Matrix{T},
    sandwich_out::Matrix{T},
) where {T<:Complex}
    CT = one(T)
    ZT = zero(T)
    BLAS.gemm!('N', 'N', CT, L_op, rho, ZT, sandwich_tmp)           # sandwich_tmp = L * rho
    BLAS.gemm!('N', 'C', CT, sandwich_tmp, L_op, ZT, sandwich_out)  # sandwich_out = L * rho * L'
    BLAS.axpy!(T(scalar), sandwich_out, out)
    return nothing
end

@inline function _accumulate_sandwich_adj_scratch!(
    out::Matrix{T},
    L_op::Matrix{T},
    rho::Matrix{T},
    scalar::Real,
    sandwich_tmp::Matrix{T},
    sandwich_out::Matrix{T},
) where {T<:Complex}
    CT = one(T)
    ZT = zero(T)
    BLAS.gemm!('C', 'N', CT, L_op, rho, ZT, sandwich_tmp)           # sandwich_tmp = L' * rho
    BLAS.gemm!('N', 'N', CT, sandwich_tmp, L_op, ZT, sandwich_out)  # sandwich_out = L' * rho * L
    BLAS.axpy!(T(scalar), sandwich_out, out)
    return nothing
end

# ---------------------------------------------------------------------------
# EnergyDomain forward and adjoint Lindbladian
# ---------------------------------------------------------------------------

"""
    apply_lindbladian!(ws, rho, config, hamiltonian) -> ws.scratch.rho_out

Compute L(rho) for `EnergyDomain` configs, storing the result in `ws.scratch.rho_out`.
"""
function apply_lindbladian!(
    ws::Workspace{Krylov},
    rho::Matrix{T},
    config::Config{Lindbladian, EnergyDomain},
    hamiltonian::HamHam,
) where {T<:Complex}
    # Function barrier: wrap transition in _TransitionWrap{F} so Julia dispatches on
    # the concrete closure type (avoids Union{Nothing,Function} boxing).
    # All other Union-typed fields are narrowed via type assertions.
    _apply_lindbladian_energy!(
        _TransitionWrap(ws.transition), ws.scratch, rho,
        ws.G_left::Matrix{T}, ws.G_right::Matrix{T},
        ws.jump_eigenbases::Vector{Matrix{T}},
        ws.jump_hermitian::Vector{Bool},
        (ws.oft_domain_prefactor::Float64) * (ws.gamma_norm_factor::Float64),
        ws.energy_labels::Vector{Float64},
        hamiltonian.bohr_freqs,
        1.0 / (4 * config.sigma^2),
    )
end

function _apply_lindbladian_energy!(
    tw::_TransitionWrap{F}, scratch::SC, rho::Matrix{T},
    G_left::Matrix{T}, G_right::Matrix{T},
    jump_eigenbases::Vector{Matrix{T}},
    jump_hermitian::Vector{Bool},
    prefactor::Float64,
    energy_labels::Vector{Float64},
    bohr_freqs::Matrix{Float64},
    inv_4sigma2::Float64,
) where {T<:Complex, F, SC}
    transition = tw.f

    CT = one(T)
    ZT = zero(T)

    BLAS.gemm!('N', 'N', CT, G_left, rho, ZT, scratch.rho_out)
    BLAS.gemm!('N', 'N', CT, rho, G_right, CT, scratch.rho_out)

    for (k, eigenbasis) in enumerate(jump_eigenbases)
        is_herm = jump_hermitian[k]
        if is_herm
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)

                oft!(scratch.jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)

                scalar_w = prefactor * transition(w)
                _accumulate_sandwich_scratch!(scratch.rho_out, scratch.jump_oft, rho, scalar_w, scratch.sandwich_tmp, scratch.sandwich_out)

                if w > 1e-12
                    scalar_neg = prefactor * transition(-w)
                    _accumulate_sandwich_adj_scratch!(scratch.rho_out, scratch.jump_oft, rho, scalar_neg, scratch.sandwich_tmp, scratch.sandwich_out)
                end
            end
        else
            for w in energy_labels
                oft!(scratch.jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)
                scalar_w = prefactor * transition(w)
                _accumulate_sandwich_scratch!(scratch.rho_out, scratch.jump_oft, rho, scalar_w, scratch.sandwich_tmp, scratch.sandwich_out)
            end
        end
    end

    return scratch.rho_out
end

"""
    apply_adjoint_lindbladian!(ws, rho, config, hamiltonian) -> ws.scratch.rho_out

Compute L*(rho) (Hilbert-Schmidt adjoint) for `EnergyDomain` configs.
"""
function apply_adjoint_lindbladian!(
    ws::Workspace{Krylov},
    rho::Matrix{T},
    config::Config{Lindbladian, EnergyDomain},
    hamiltonian::HamHam,
) where {T<:Complex}
    _apply_adjoint_lindbladian_energy!(
        _TransitionWrap(ws.transition), ws.scratch, rho,
        ws.G_left_adj::Matrix{T}, ws.G_right_adj::Matrix{T},
        ws.jump_eigenbases::Vector{Matrix{T}},
        ws.jump_hermitian::Vector{Bool},
        (ws.oft_domain_prefactor::Float64) * (ws.gamma_norm_factor::Float64),
        ws.energy_labels::Vector{Float64},
        hamiltonian.bohr_freqs,
        1.0 / (4 * config.sigma^2),
    )
end

function _apply_adjoint_lindbladian_energy!(
    tw::_TransitionWrap{F}, scratch::SC, rho::Matrix{T},
    G_left_adj::Matrix{T}, G_right_adj::Matrix{T},
    jump_eigenbases::Vector{Matrix{T}},
    jump_hermitian::Vector{Bool},
    prefactor::Float64,
    energy_labels::Vector{Float64},
    bohr_freqs::Matrix{Float64},
    inv_4sigma2::Float64,
) where {T<:Complex, F, SC}
    transition = tw.f

    CT = one(T)
    ZT = zero(T)

    BLAS.gemm!('N', 'N', CT, G_left_adj, rho, ZT, scratch.rho_out)
    BLAS.gemm!('N', 'N', CT, rho, G_right_adj, CT, scratch.rho_out)

    for (k, eigenbasis) in enumerate(jump_eigenbases)
        is_herm = jump_hermitian[k]
        if is_herm
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)

                oft!(scratch.jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)

                scalar_w = prefactor * transition(w)
                _accumulate_sandwich_adj_scratch!(scratch.rho_out, scratch.jump_oft, rho, scalar_w, scratch.sandwich_tmp, scratch.sandwich_out)

                if w > 1e-12
                    scalar_neg = prefactor * transition(-w)
                    _accumulate_sandwich_scratch!(scratch.rho_out, scratch.jump_oft, rho, scalar_neg, scratch.sandwich_tmp, scratch.sandwich_out)
                end
            end
        else
            for w in energy_labels
                oft!(scratch.jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)
                scalar_w = prefactor * transition(w)
                _accumulate_sandwich_adj_scratch!(scratch.rho_out, scratch.jump_oft, rho, scalar_w, scratch.sandwich_tmp, scratch.sandwich_out)
            end
        end
    end

    return scratch.rho_out
end

# ---------------------------------------------------------------------------
# BohrDomain forward and adjoint Lindbladian
# ---------------------------------------------------------------------------

"""
    _accumulate_sandwich_2op!(out, A, B_dag, rho, scalar, ws) -> nothing

Accumulate `scalar * A * rho * B_dag` into `out`. BohrDomain two-operator sandwich.
"""
function _accumulate_sandwich_2op!(
    out::Matrix{T},
    A::Matrix{T},
    B_dag::Matrix{T},
    rho::Matrix{T},
    scalar::Real,
    ws::Workspace{Krylov},
) where {T<:Complex}
    CT = one(T)
    ZT = zero(T)
    BLAS.gemm!('N', 'N', CT, A, rho, ZT, ws.scratch.sandwich_tmp)
    BLAS.gemm!('N', 'N', CT, ws.scratch.sandwich_tmp, B_dag, ZT, ws.scratch.sandwich_out)
    BLAS.axpy!(T(scalar), ws.scratch.sandwich_out, out)
    return nothing
end

"""
    _accumulate_adjoint_sandwich_2op!(out, A, B_dag, rho, scalar, ws) -> nothing

Accumulate `scalar * A' * rho * B_dag'` into `out`. BohrDomain adjoint two-operator sandwich.
"""
function _accumulate_adjoint_sandwich_2op!(
    out::Matrix{T},
    A::Matrix{T},
    B_dag::Matrix{T},
    rho::Matrix{T},
    scalar::Real,
    ws::Workspace{Krylov},
) where {T<:Complex}
    CT = one(T)
    ZT = zero(T)
    BLAS.gemm!('C', 'N', CT, A, rho, ZT, ws.scratch.sandwich_tmp)
    BLAS.gemm!('N', 'C', CT, ws.scratch.sandwich_tmp, B_dag, ZT, ws.scratch.sandwich_out)
    BLAS.axpy!(T(scalar), ws.scratch.sandwich_out, out)
    return nothing
end

"""
    apply_lindbladian!(ws, rho, config, hamiltonian) -> ws.scratch.rho_out

Compute L(rho) for `BohrDomain` configs.
"""
function apply_lindbladian!(
    ws::Workspace{Krylov},
    rho::Matrix{T},
    config::Config{Lindbladian, BohrDomain},
    hamiltonian::HamHam,
) where {T<:Complex}
    bohr_alpha = ws.bohr_alpha
    gamma_norm_factor = ws.gamma_norm_factor
    dim = size(rho, 1)

    CT = one(T)
    ZT = zero(T)

    BLAS.gemm!('N', 'N', CT, ws.G_left, rho, ZT, ws.scratch.rho_out)
    BLAS.gemm!('N', 'N', CT, rho, ws.G_right, CT, ws.scratch.rho_out)

    A_nu2_dag = zeros(T, dim, dim)

    for (k, eigenbasis) in enumerate(ws.jump_eigenbases)
        for nu_2 in keys(hamiltonian.bohr_dict)
            @. ws.scratch.jump_oft = bohr_alpha(hamiltonian.bohr_freqs, nu_2) * eigenbasis

            fill!(A_nu2_dag, 0)
            indices = hamiltonian.bohr_dict[nu_2]
            @inbounds for idx in indices
                i = idx[1]; j = idx[2]
                A_nu2_dag[j, i] = conj(eigenbasis[i, j])
            end

            _accumulate_sandwich_2op!(ws.scratch.rho_out, ws.scratch.jump_oft, A_nu2_dag, rho, gamma_norm_factor, ws)
        end
    end

    return ws.scratch.rho_out
end

"""
    apply_adjoint_lindbladian!(ws, rho, config, hamiltonian) -> ws.scratch.rho_out

Compute L*(rho) (Hilbert-Schmidt adjoint) for `BohrDomain` configs.
"""
function apply_adjoint_lindbladian!(
    ws::Workspace{Krylov},
    rho::Matrix{T},
    config::Config{Lindbladian, BohrDomain},
    hamiltonian::HamHam,
) where {T<:Complex}
    bohr_alpha = ws.bohr_alpha
    gamma_norm_factor = ws.gamma_norm_factor
    dim = size(rho, 1)

    CT = one(T)
    ZT = zero(T)

    BLAS.gemm!('N', 'N', CT, ws.G_left_adj, rho, ZT, ws.scratch.rho_out)
    BLAS.gemm!('N', 'N', CT, rho, ws.G_right_adj, CT, ws.scratch.rho_out)

    A_nu2_dag = zeros(T, dim, dim)

    for (k, eigenbasis) in enumerate(ws.jump_eigenbases)
        for nu_2 in keys(hamiltonian.bohr_dict)
            @. ws.scratch.jump_oft = bohr_alpha(hamiltonian.bohr_freqs, nu_2) * eigenbasis

            fill!(A_nu2_dag, 0)
            indices = hamiltonian.bohr_dict[nu_2]
            @inbounds for idx in indices
                i = idx[1]; j = idx[2]
                A_nu2_dag[j, i] = conj(eigenbasis[i, j])
            end

            _accumulate_adjoint_sandwich_2op!(ws.scratch.rho_out, ws.scratch.jump_oft, A_nu2_dag, rho, gamma_norm_factor, ws)
        end
    end

    return ws.scratch.rho_out
end

# ---------------------------------------------------------------------------
# TimeDomain / TrotterDomain forward and adjoint Lindbladian
# ---------------------------------------------------------------------------

"""
    apply_lindbladian!(ws, rho, config, hamiltonian) -> ws.scratch.rho_out

Compute L(rho) for `TimeDomain` and `TrotterDomain` configs.
"""
function apply_lindbladian!(
    ws::Workspace{Krylov},
    rho::Matrix{T},
    config::Config{Lindbladian, D},
    hamiltonian::HamHam,
) where {T<:Complex, D<:Union{TimeDomain, TrotterDomain}}
    _apply_lindbladian_time!(
        _TransitionWrap(ws.transition), ws.oft_nufft_prefactors,
        ws.scratch, rho,
        ws.G_left::Matrix{T}, ws.G_right::Matrix{T},
        ws.jump_eigenbases::Vector{Matrix{T}},
        ws.jump_hermitian::Vector{Bool},
        (ws.oft_domain_prefactor::Float64) * (ws.gamma_norm_factor::Float64),
        ws.energy_labels::Vector{Float64},
    )
end

function _apply_lindbladian_time!(
    tw::_TransitionWrap{F}, oft_nufft_prefactors::P,
    scratch::SC, rho::Matrix{T},
    G_left::Matrix{T}, G_right::Matrix{T},
    jump_eigenbases::Vector{Matrix{T}},
    jump_hermitian::Vector{Bool},
    prefactor::Float64,
    energy_labels::Vector{Float64},
) where {T<:Complex, F, P, SC}
    transition = tw.f

    CT = one(T)
    ZT = zero(T)

    BLAS.gemm!('N', 'N', CT, G_left, rho, ZT, scratch.rho_out)
    BLAS.gemm!('N', 'N', CT, rho, G_right, CT, scratch.rho_out)

    for (k, eigenbasis) in enumerate(jump_eigenbases)
        is_herm = jump_hermitian[k]
        if is_herm
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)

                nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)
                @. scratch.jump_oft = eigenbasis * nufft_prefactor_matrix

                scalar_w = prefactor * transition(w)
                _accumulate_sandwich_scratch!(scratch.rho_out, scratch.jump_oft, rho, scalar_w, scratch.sandwich_tmp, scratch.sandwich_out)

                if w > 1e-12
                    scalar_neg = prefactor * transition(-w)
                    _accumulate_sandwich_adj_scratch!(scratch.rho_out, scratch.jump_oft, rho, scalar_neg, scratch.sandwich_tmp, scratch.sandwich_out)
                end
            end
        else
            for w in energy_labels
                nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)
                @. scratch.jump_oft = eigenbasis * nufft_prefactor_matrix
                scalar_w = prefactor * transition(w)
                _accumulate_sandwich_scratch!(scratch.rho_out, scratch.jump_oft, rho, scalar_w, scratch.sandwich_tmp, scratch.sandwich_out)
            end
        end
    end

    return scratch.rho_out
end

"""
    apply_adjoint_lindbladian!(ws, rho, config, hamiltonian) -> ws.scratch.rho_out

Compute L*(rho) (Hilbert-Schmidt adjoint) for `TimeDomain` and `TrotterDomain` configs.
"""
function apply_adjoint_lindbladian!(
    ws::Workspace{Krylov},
    rho::Matrix{T},
    config::Config{Lindbladian, D},
    hamiltonian::HamHam,
) where {T<:Complex, D<:Union{TimeDomain, TrotterDomain}}
    _apply_adjoint_lindbladian_time!(
        _TransitionWrap(ws.transition), ws.oft_nufft_prefactors,
        ws.scratch, rho,
        ws.G_left_adj::Matrix{T}, ws.G_right_adj::Matrix{T},
        ws.jump_eigenbases::Vector{Matrix{T}},
        ws.jump_hermitian::Vector{Bool},
        (ws.oft_domain_prefactor::Float64) * (ws.gamma_norm_factor::Float64),
        ws.energy_labels::Vector{Float64},
    )
end

function _apply_adjoint_lindbladian_time!(
    tw::_TransitionWrap{F}, oft_nufft_prefactors::P,
    scratch::SC, rho::Matrix{T},
    G_left_adj::Matrix{T}, G_right_adj::Matrix{T},
    jump_eigenbases::Vector{Matrix{T}},
    jump_hermitian::Vector{Bool},
    prefactor::Float64,
    energy_labels::Vector{Float64},
) where {T<:Complex, F, P, SC}
    transition = tw.f

    CT = one(T)
    ZT = zero(T)

    BLAS.gemm!('N', 'N', CT, G_left_adj, rho, ZT, scratch.rho_out)
    BLAS.gemm!('N', 'N', CT, rho, G_right_adj, CT, scratch.rho_out)

    for (k, eigenbasis) in enumerate(jump_eigenbases)
        is_herm = jump_hermitian[k]
        if is_herm
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)

                nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)
                @. scratch.jump_oft = eigenbasis * nufft_prefactor_matrix

                scalar_w = prefactor * transition(w)
                _accumulate_sandwich_adj_scratch!(scratch.rho_out, scratch.jump_oft, rho, scalar_w, scratch.sandwich_tmp, scratch.sandwich_out)

                if w > 1e-12
                    scalar_neg = prefactor * transition(-w)
                    _accumulate_sandwich_scratch!(scratch.rho_out, scratch.jump_oft, rho, scalar_neg, scratch.sandwich_tmp, scratch.sandwich_out)
                end
            end
        else
            for w in energy_labels
                nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)
                @. scratch.jump_oft = eigenbasis * nufft_prefactor_matrix
                scalar_w = prefactor * transition(w)
                _accumulate_sandwich_adj_scratch!(scratch.rho_out, scratch.jump_oft, rho, scalar_w, scratch.sandwich_tmp, scratch.sandwich_out)
            end
        end
    end

    return scratch.rho_out
end
