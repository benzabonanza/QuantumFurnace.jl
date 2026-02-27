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

Uses `sc.sandwich_tmp`, `sc.sandwich_out` as scratch (sc obtained via type assertion).
"""
@inline function _accumulate_sandwich!(
    out::Matrix{T},
    L_op::Matrix{T},
    rho::Matrix{T},
    scalar::Real,
    ws::Workspace{KrylovSpectrum},
) where {T<:Complex}
    sc = ws.scratch::KrylovScratch{T}
    _accumulate_sandwich_scratch!(out, L_op, rho, scalar, sc.sandwich_tmp, sc.sandwich_out)
end

"""
    _accumulate_sandwich_adj!(out, L_op, rho, scalar, ws) -> nothing

Accumulate `scalar * L' * rho * L` into `out`. Used for:
- Forward Lindbladian negative-frequency Hermitian partner (L_neg = L')
- Adjoint Lindbladian positive-frequency sandwich (HS adjoint of L*rho*L' is L'*rho*L)

Uses `sc.sandwich_tmp`, `sc.sandwich_out` as scratch (sc obtained via type assertion).
"""
@inline function _accumulate_sandwich_adj!(
    out::Matrix{T},
    L_op::Matrix{T},
    rho::Matrix{T},
    scalar::Real,
    ws::Workspace{KrylovSpectrum},
) where {T<:Complex}
    sc = ws.scratch::KrylovScratch{T}
    _accumulate_sandwich_adj_scratch!(out, L_op, rho, scalar, sc.sandwich_tmp, sc.sandwich_out)
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
    apply_lindbladian!(ws, rho, config, hamiltonian) -> sc.rho_out

Compute L(rho) for `EnergyDomain` configs, storing the result in `sc.rho_out`.
"""
function apply_lindbladian!(
    ws::Workspace{KrylovSpectrum},
    rho::Matrix{T},
    config::Config{Lindbladian, EnergyDomain},
    hamiltonian::HamHam,
) where {T<:Complex}
    sc = ws.scratch::KrylovScratch{T}
    transition = ws.transition
    G_left = ws.G_left::Matrix{T}
    G_right = ws.G_right::Matrix{T}
    jump_eigenbases = ws.jump_eigenbases::Vector{Matrix{T}}
    jump_hermitian = ws.jump_hermitian::Vector{Bool}
    prefactor = (ws.oft_domain_prefactor::Float64) * (ws.gamma_norm_factor::Float64)
    energy_labels = ws.energy_labels::Vector{Float64}
    bohr_freqs = hamiltonian.bohr_freqs
    inv_4sigma2 = 1.0 / (4 * config.sigma^2)

    CT = one(T)
    ZT = zero(T)

    BLAS.gemm!('N', 'N', CT, G_left, rho, ZT, sc.rho_out)
    BLAS.gemm!('N', 'N', CT, rho, G_right, CT, sc.rho_out)

    for (k, eigenbasis) in enumerate(jump_eigenbases)
        is_herm = jump_hermitian[k]
        if is_herm
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)

                oft!(sc.jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)

                scalar_w = prefactor * transition(w)
                _accumulate_sandwich_scratch!(sc.rho_out, sc.jump_oft, rho, scalar_w, sc.sandwich_tmp, sc.sandwich_out)

                if w > 1e-12
                    scalar_neg = prefactor * transition(-w)
                    _accumulate_sandwich_adj_scratch!(sc.rho_out, sc.jump_oft, rho, scalar_neg, sc.sandwich_tmp, sc.sandwich_out)
                end
            end
        else
            for w in energy_labels
                oft!(sc.jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)
                scalar_w = prefactor * transition(w)
                _accumulate_sandwich_scratch!(sc.rho_out, sc.jump_oft, rho, scalar_w, sc.sandwich_tmp, sc.sandwich_out)
            end
        end
    end

    return sc.rho_out
end

"""
    apply_adjoint_lindbladian!(ws, rho, config, hamiltonian) -> sc.rho_out

Compute L*(rho) (Hilbert-Schmidt adjoint) for `EnergyDomain` configs.
"""
function apply_adjoint_lindbladian!(
    ws::Workspace{KrylovSpectrum},
    rho::Matrix{T},
    config::Config{Lindbladian, EnergyDomain},
    hamiltonian::HamHam,
) where {T<:Complex}
    sc = ws.scratch::KrylovScratch{T}
    transition = ws.transition
    G_left_adj = ws.G_left_adj::Matrix{T}
    G_right_adj = ws.G_right_adj::Matrix{T}
    jump_eigenbases = ws.jump_eigenbases::Vector{Matrix{T}}
    jump_hermitian = ws.jump_hermitian::Vector{Bool}
    prefactor = (ws.oft_domain_prefactor::Float64) * (ws.gamma_norm_factor::Float64)
    energy_labels = ws.energy_labels::Vector{Float64}
    bohr_freqs = hamiltonian.bohr_freqs
    inv_4sigma2 = 1.0 / (4 * config.sigma^2)

    CT = one(T)
    ZT = zero(T)

    BLAS.gemm!('N', 'N', CT, G_left_adj, rho, ZT, sc.rho_out)
    BLAS.gemm!('N', 'N', CT, rho, G_right_adj, CT, sc.rho_out)

    for (k, eigenbasis) in enumerate(jump_eigenbases)
        is_herm = jump_hermitian[k]
        if is_herm
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)

                oft!(sc.jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)

                scalar_w = prefactor * transition(w)
                _accumulate_sandwich_adj_scratch!(sc.rho_out, sc.jump_oft, rho, scalar_w, sc.sandwich_tmp, sc.sandwich_out)

                if w > 1e-12
                    scalar_neg = prefactor * transition(-w)
                    _accumulate_sandwich_scratch!(sc.rho_out, sc.jump_oft, rho, scalar_neg, sc.sandwich_tmp, sc.sandwich_out)
                end
            end
        else
            for w in energy_labels
                oft!(sc.jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)
                scalar_w = prefactor * transition(w)
                _accumulate_sandwich_adj_scratch!(sc.rho_out, sc.jump_oft, rho, scalar_w, sc.sandwich_tmp, sc.sandwich_out)
            end
        end
    end

    return sc.rho_out
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
    ws::Workspace{KrylovSpectrum},
) where {T<:Complex}
    sc = ws.scratch::KrylovScratch{T}
    CT = one(T)
    ZT = zero(T)
    BLAS.gemm!('N', 'N', CT, A, rho, ZT, sc.sandwich_tmp)
    BLAS.gemm!('N', 'N', CT, sc.sandwich_tmp, B_dag, ZT, sc.sandwich_out)
    BLAS.axpy!(T(scalar), sc.sandwich_out, out)
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
    ws::Workspace{KrylovSpectrum},
) where {T<:Complex}
    sc = ws.scratch::KrylovScratch{T}
    CT = one(T)
    ZT = zero(T)
    BLAS.gemm!('C', 'N', CT, A, rho, ZT, sc.sandwich_tmp)
    BLAS.gemm!('N', 'C', CT, sc.sandwich_tmp, B_dag, ZT, sc.sandwich_out)
    BLAS.axpy!(T(scalar), sc.sandwich_out, out)
    return nothing
end

"""
    apply_lindbladian!(ws, rho, config, hamiltonian) -> sc.rho_out

Compute L(rho) for `BohrDomain` configs.
"""
function apply_lindbladian!(
    ws::Workspace{KrylovSpectrum},
    rho::Matrix{T},
    config::Config{Lindbladian, BohrDomain},
    hamiltonian::HamHam,
) where {T<:Complex}
    sc = ws.scratch::KrylovScratch{T}
    bohr_alpha = ws.bohr_alpha
    gamma_norm_factor = ws.gamma_norm_factor
    dim = size(rho, 1)

    CT = one(T)
    ZT = zero(T)

    BLAS.gemm!('N', 'N', CT, ws.G_left, rho, ZT, sc.rho_out)
    BLAS.gemm!('N', 'N', CT, rho, ws.G_right, CT, sc.rho_out)

    A_nu2_dag = zeros(T, dim, dim)

    for (k, eigenbasis) in enumerate(ws.jump_eigenbases)
        for nu_2 in keys(hamiltonian.bohr_dict)
            @. sc.jump_oft = bohr_alpha(hamiltonian.bohr_freqs, nu_2) * eigenbasis

            fill!(A_nu2_dag, 0)
            indices = hamiltonian.bohr_dict[nu_2]
            @inbounds for idx in indices
                i = idx[1]; j = idx[2]
                A_nu2_dag[j, i] = conj(eigenbasis[i, j])
            end

            _accumulate_sandwich_2op!(sc.rho_out, sc.jump_oft, A_nu2_dag, rho, gamma_norm_factor, ws)
        end
    end

    return sc.rho_out
end

"""
    apply_adjoint_lindbladian!(ws, rho, config, hamiltonian) -> sc.rho_out

Compute L*(rho) (Hilbert-Schmidt adjoint) for `BohrDomain` configs.
"""
function apply_adjoint_lindbladian!(
    ws::Workspace{KrylovSpectrum},
    rho::Matrix{T},
    config::Config{Lindbladian, BohrDomain},
    hamiltonian::HamHam,
) where {T<:Complex}
    sc = ws.scratch::KrylovScratch{T}
    bohr_alpha = ws.bohr_alpha
    gamma_norm_factor = ws.gamma_norm_factor
    dim = size(rho, 1)

    CT = one(T)
    ZT = zero(T)

    BLAS.gemm!('N', 'N', CT, ws.G_left_adj, rho, ZT, sc.rho_out)
    BLAS.gemm!('N', 'N', CT, rho, ws.G_right_adj, CT, sc.rho_out)

    A_nu2_dag = zeros(T, dim, dim)

    for (k, eigenbasis) in enumerate(ws.jump_eigenbases)
        for nu_2 in keys(hamiltonian.bohr_dict)
            @. sc.jump_oft = bohr_alpha(hamiltonian.bohr_freqs, nu_2) * eigenbasis

            fill!(A_nu2_dag, 0)
            indices = hamiltonian.bohr_dict[nu_2]
            @inbounds for idx in indices
                i = idx[1]; j = idx[2]
                A_nu2_dag[j, i] = conj(eigenbasis[i, j])
            end

            _accumulate_adjoint_sandwich_2op!(sc.rho_out, sc.jump_oft, A_nu2_dag, rho, gamma_norm_factor, ws)
        end
    end

    return sc.rho_out
end

# ---------------------------------------------------------------------------
# TimeDomain / TrotterDomain forward and adjoint Lindbladian
# ---------------------------------------------------------------------------

"""
    apply_lindbladian!(ws, rho, config, hamiltonian) -> sc.rho_out

Compute L(rho) for `TimeDomain` and `TrotterDomain` configs.
"""
function apply_lindbladian!(
    ws::Workspace{KrylovSpectrum},
    rho::Matrix{T},
    config::Config{Lindbladian, D},
    hamiltonian::HamHam,
) where {T<:Complex, D<:Union{TimeDomain, TrotterDomain}}
    sc = ws.scratch::KrylovScratch{T}
    transition = ws.transition
    oft_nufft_prefactors = ws.oft_nufft_prefactors
    G_left = ws.G_left::Matrix{T}
    G_right = ws.G_right::Matrix{T}
    jump_eigenbases = ws.jump_eigenbases::Vector{Matrix{T}}
    jump_hermitian = ws.jump_hermitian::Vector{Bool}
    prefactor = (ws.oft_domain_prefactor::Float64) * (ws.gamma_norm_factor::Float64)
    energy_labels = ws.energy_labels::Vector{Float64}

    CT = one(T)
    ZT = zero(T)

    BLAS.gemm!('N', 'N', CT, G_left, rho, ZT, sc.rho_out)
    BLAS.gemm!('N', 'N', CT, rho, G_right, CT, sc.rho_out)

    for (k, eigenbasis) in enumerate(jump_eigenbases)
        is_herm = jump_hermitian[k]
        if is_herm
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)

                nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)
                @. sc.jump_oft = eigenbasis * nufft_prefactor_matrix

                scalar_w = prefactor * transition(w)
                _accumulate_sandwich_scratch!(sc.rho_out, sc.jump_oft, rho, scalar_w, sc.sandwich_tmp, sc.sandwich_out)

                if w > 1e-12
                    scalar_neg = prefactor * transition(-w)
                    _accumulate_sandwich_adj_scratch!(sc.rho_out, sc.jump_oft, rho, scalar_neg, sc.sandwich_tmp, sc.sandwich_out)
                end
            end
        else
            for w in energy_labels
                nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)
                @. sc.jump_oft = eigenbasis * nufft_prefactor_matrix
                scalar_w = prefactor * transition(w)
                _accumulate_sandwich_scratch!(sc.rho_out, sc.jump_oft, rho, scalar_w, sc.sandwich_tmp, sc.sandwich_out)
            end
        end
    end

    return sc.rho_out
end

"""
    apply_adjoint_lindbladian!(ws, rho, config, hamiltonian) -> sc.rho_out

Compute L*(rho) (Hilbert-Schmidt adjoint) for `TimeDomain` and `TrotterDomain` configs.
"""
function apply_adjoint_lindbladian!(
    ws::Workspace{KrylovSpectrum},
    rho::Matrix{T},
    config::Config{Lindbladian, D},
    hamiltonian::HamHam,
) where {T<:Complex, D<:Union{TimeDomain, TrotterDomain}}
    sc = ws.scratch::KrylovScratch{T}
    transition = ws.transition
    oft_nufft_prefactors = ws.oft_nufft_prefactors
    G_left_adj = ws.G_left_adj::Matrix{T}
    G_right_adj = ws.G_right_adj::Matrix{T}
    jump_eigenbases = ws.jump_eigenbases::Vector{Matrix{T}}
    jump_hermitian = ws.jump_hermitian::Vector{Bool}
    prefactor = (ws.oft_domain_prefactor::Float64) * (ws.gamma_norm_factor::Float64)
    energy_labels = ws.energy_labels::Vector{Float64}

    CT = one(T)
    ZT = zero(T)

    BLAS.gemm!('N', 'N', CT, G_left_adj, rho, ZT, sc.rho_out)
    BLAS.gemm!('N', 'N', CT, rho, G_right_adj, CT, sc.rho_out)

    for (k, eigenbasis) in enumerate(jump_eigenbases)
        is_herm = jump_hermitian[k]
        if is_herm
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)

                nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)
                @. sc.jump_oft = eigenbasis * nufft_prefactor_matrix

                scalar_w = prefactor * transition(w)
                _accumulate_sandwich_adj_scratch!(sc.rho_out, sc.jump_oft, rho, scalar_w, sc.sandwich_tmp, sc.sandwich_out)

                if w > 1e-12
                    scalar_neg = prefactor * transition(-w)
                    _accumulate_sandwich_scratch!(sc.rho_out, sc.jump_oft, rho, scalar_neg, sc.sandwich_tmp, sc.sandwich_out)
                end
            end
        else
            for w in energy_labels
                nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)
                @. sc.jump_oft = eigenbasis * nufft_prefactor_matrix
                scalar_w = prefactor * transition(w)
                _accumulate_sandwich_adj_scratch!(sc.rho_out, sc.jump_oft, rho, scalar_w, sc.sandwich_tmp, sc.sandwich_out)
            end
        end
    end

    return sc.rho_out
end
