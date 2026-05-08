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
    hamiltonian::HamHam;
    include_coherent::Bool = true,
) where {T<:Complex}
    sc = ws.scratch::KrylovScratch{T}
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

    if include_coherent
        BLAS.gemm!('N', 'N', CT, G_left, rho, ZT, sc.rho_out)
        BLAS.gemm!('N', 'N', CT, rho, G_right, CT, sc.rho_out)
    else
        neg_R = sc.sandwich_tmp
        @. neg_R = G_left + G_right
        half = T(0.5)
        BLAS.gemm!('N', 'N', half, neg_R, rho, ZT, sc.rho_out)
        BLAS.gemm!('N', 'N', half, rho, neg_R, CT, sc.rho_out)
    end

    if Threads.nthreads() > 1 && length(energy_labels) >= OMEGA_THREAD_THRESHOLD
        return _apply_lindbladian_threaded_energy!(
            sc, rho, jump_eigenbases, jump_hermitian, bohr_freqs,
            energy_labels, config, prefactor, inv_4sigma2; adjoint=false)
    end

    for (k, eigenbasis) in enumerate(jump_eigenbases)
        is_herm = jump_hermitian[k]
        if is_herm
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)

                oft!(sc.jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)

                scalar_w = prefactor * pick_transition(config, w)
                _accumulate_sandwich_scratch!(sc.rho_out, sc.jump_oft, rho, scalar_w, sc.sandwich_tmp, sc.sandwich_out)

                if w > 1e-12
                    scalar_neg = prefactor * pick_transition(config, -w)
                    _accumulate_sandwich_adj_scratch!(sc.rho_out, sc.jump_oft, rho, scalar_neg, sc.sandwich_tmp, sc.sandwich_out)
                end
            end
        else
            for w in energy_labels
                oft!(sc.jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)
                scalar_w = prefactor * pick_transition(config, w)
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
    hamiltonian::HamHam;
    include_coherent::Bool = true,
) where {T<:Complex}
    sc = ws.scratch::KrylovScratch{T}
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

    if include_coherent
        BLAS.gemm!('N', 'N', CT, G_left_adj, rho, ZT, sc.rho_out)
        BLAS.gemm!('N', 'N', CT, rho, G_right_adj, CT, sc.rho_out)
    else
        neg_R = sc.sandwich_tmp
        @. neg_R = G_left_adj + G_right_adj
        half = T(0.5)
        BLAS.gemm!('N', 'N', half, neg_R, rho, ZT, sc.rho_out)
        BLAS.gemm!('N', 'N', half, rho, neg_R, CT, sc.rho_out)
    end

    if Threads.nthreads() > 1 && length(energy_labels) >= OMEGA_THREAD_THRESHOLD
        return _apply_lindbladian_threaded_energy!(
            sc, rho, jump_eigenbases, jump_hermitian, bohr_freqs,
            energy_labels, config, prefactor, inv_4sigma2; adjoint=true)
    end

    for (k, eigenbasis) in enumerate(jump_eigenbases)
        is_herm = jump_hermitian[k]
        if is_herm
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)

                oft!(sc.jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)

                scalar_w = prefactor * pick_transition(config, w)
                _accumulate_sandwich_adj_scratch!(sc.rho_out, sc.jump_oft, rho, scalar_w, sc.sandwich_tmp, sc.sandwich_out)

                if w > 1e-12
                    scalar_neg = prefactor * pick_transition(config, -w)
                    _accumulate_sandwich_scratch!(sc.rho_out, sc.jump_oft, rho, scalar_neg, sc.sandwich_tmp, sc.sandwich_out)
                end
            end
        else
            for w in energy_labels
                oft!(sc.jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)
                scalar_w = prefactor * pick_transition(config, w)
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
    hamiltonian::HamHam;
    include_coherent::Bool = true,
) where {T<:Complex}
    sc = ws.scratch::KrylovScratch{T}
    gamma_norm_factor = ws.gamma_norm_factor
    dim = size(rho, 1)

    CT = one(T)
    ZT = zero(T)

    if include_coherent
        BLAS.gemm!('N', 'N', CT, ws.G_left, rho, ZT, sc.rho_out)
        BLAS.gemm!('N', 'N', CT, rho, ws.G_right, CT, sc.rho_out)
    else
        neg_R = sc.sandwich_tmp
        @. neg_R = ws.G_left + ws.G_right
        half = T(0.5)
        BLAS.gemm!('N', 'N', half, neg_R, rho, ZT, sc.rho_out)
        BLAS.gemm!('N', 'N', half, rho, neg_R, CT, sc.rho_out)
    end

    A_nu2_dag = zeros(T, dim, dim)

    for (k, eigenbasis) in enumerate(ws.jump_eigenbases)
        for nu_2 in keys(hamiltonian.bohr_dict)
            @. sc.jump_oft = _pick_alpha($Ref(config), hamiltonian.bohr_freqs, nu_2) * eigenbasis

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
    hamiltonian::HamHam;
    include_coherent::Bool = true,
) where {T<:Complex}
    sc = ws.scratch::KrylovScratch{T}
    gamma_norm_factor = ws.gamma_norm_factor
    dim = size(rho, 1)

    CT = one(T)
    ZT = zero(T)

    if include_coherent
        BLAS.gemm!('N', 'N', CT, ws.G_left_adj, rho, ZT, sc.rho_out)
        BLAS.gemm!('N', 'N', CT, rho, ws.G_right_adj, CT, sc.rho_out)
    else
        neg_R = sc.sandwich_tmp
        @. neg_R = ws.G_left_adj + ws.G_right_adj
        half = T(0.5)
        BLAS.gemm!('N', 'N', half, neg_R, rho, ZT, sc.rho_out)
        BLAS.gemm!('N', 'N', half, rho, neg_R, CT, sc.rho_out)
    end

    A_nu2_dag = zeros(T, dim, dim)

    for (k, eigenbasis) in enumerate(ws.jump_eigenbases)
        for nu_2 in keys(hamiltonian.bohr_dict)
            @. sc.jump_oft = _pick_alpha($Ref(config), hamiltonian.bohr_freqs, nu_2) * eigenbasis

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
# BohrDomain DLL forward and adjoint Lindbladian (qf-lkb.9)
# ---------------------------------------------------------------------------
#
# DLL collapses the CKG outer ω-loop into a single Lindblad operator per
# coupling: `L_a = dll_lindblad_op_bohr(jump, hamiltonian, filter)` (Ding-Li-Lin
# 2024 Eq. 3.4 first form). The matrix-free hot path is:
#
#   L(ρ)  = G_left · ρ + ρ · G_right + Σ_a L_a · ρ · L_a†
#   L*(ρ) = G_left_adj · ρ + ρ · G_right_adj + Σ_a L_a† · ρ · L_a
#
# with `G_left = +1im · transpose(G) − 0.5 · R_total`, `G_right = −1im ·
# transpose(G) − 0.5 · R_total`, `G_left_adj = G_right`, `G_right_adj = G_left`
# (precomputed in the DLL specialised `Workspace` constructor). The sandwich
# helpers are reused from the EnergyDomain/TimeDomain path so allocation
# profile matches CKG.

"""
    apply_lindbladian!(ws, rho, config, hamiltonian) -> sc.rho_out

Compute L(rho) for `Config{Lindbladian, BohrDomain, DLL}`. Reuses
`_accumulate_sandwich_scratch!` for the per-jump `L_a · ρ · L_a†` term.
"""
function apply_lindbladian!(
    ws::Workspace{KrylovSpectrum},
    rho::Matrix{T},
    config::Config{Lindbladian, BohrDomain, DLL},
    hamiltonian::HamHam;
    include_coherent::Bool = true,
) where {T<:Complex}
    sc = ws.scratch::KrylovScratch{T}
    G_left  = ws.G_left::Matrix{T}
    G_right = ws.G_right::Matrix{T}
    dll_lindblads = ws.dll_lindblads::Vector{Matrix{T}}

    CT = one(T)
    ZT = zero(T)

    if include_coherent
        BLAS.gemm!('N', 'N', CT, G_left, rho, ZT, sc.rho_out)
        BLAS.gemm!('N', 'N', CT, rho, G_right, CT, sc.rho_out)
    else
        neg_R = sc.sandwich_tmp
        @. neg_R = G_left + G_right
        half = T(0.5)
        BLAS.gemm!('N', 'N', half, neg_R, rho, ZT, sc.rho_out)
        BLAS.gemm!('N', 'N', half, rho, neg_R, CT, sc.rho_out)
    end

    for L_a in dll_lindblads
        _accumulate_sandwich_scratch!(sc.rho_out, L_a, rho, 1.0,
                                      sc.sandwich_tmp, sc.sandwich_out)
    end

    return sc.rho_out
end

"""
    apply_adjoint_lindbladian!(ws, rho, config, hamiltonian) -> sc.rho_out

Compute L*(rho) (Hilbert-Schmidt adjoint) for `Config{Lindbladian, BohrDomain,
DLL}`. The HS-adjoint of `L_a · ρ · L_a†` is `L_a† · ρ · L_a`; the coherent
contribution sign-flips through `G_left_adj = G_right`, `G_right_adj = G_left`.
"""
function apply_adjoint_lindbladian!(
    ws::Workspace{KrylovSpectrum},
    rho::Matrix{T},
    config::Config{Lindbladian, BohrDomain, DLL},
    hamiltonian::HamHam;
    include_coherent::Bool = true,
) where {T<:Complex}
    sc = ws.scratch::KrylovScratch{T}
    G_left_adj  = ws.G_left_adj::Matrix{T}
    G_right_adj = ws.G_right_adj::Matrix{T}
    dll_lindblads = ws.dll_lindblads::Vector{Matrix{T}}

    CT = one(T)
    ZT = zero(T)

    if include_coherent
        BLAS.gemm!('N', 'N', CT, G_left_adj, rho, ZT, sc.rho_out)
        BLAS.gemm!('N', 'N', CT, rho, G_right_adj, CT, sc.rho_out)
    else
        neg_R = sc.sandwich_tmp
        @. neg_R = G_left_adj + G_right_adj
        half = T(0.5)
        BLAS.gemm!('N', 'N', half, neg_R, rho, ZT, sc.rho_out)
        BLAS.gemm!('N', 'N', half, rho, neg_R, CT, sc.rho_out)
    end

    for L_a in dll_lindblads
        _accumulate_sandwich_adj_scratch!(sc.rho_out, L_a, rho, 1.0,
                                          sc.sandwich_tmp, sc.sandwich_out)
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
    hamiltonian::HamHam;
    include_coherent::Bool = true,
) where {T<:Complex, D<:Union{TimeDomain, TrotterDomain}}
    sc = ws.scratch::KrylovScratch{T}
    _nufft = ws.oft_nufft_prefactors::NUFFTPrefactors{real(T), Array{T, 3}}
    nufft_data = _nufft.data
    nufft_idx = _nufft.energy_to_index
    G_left = ws.G_left::Matrix{T}
    G_right = ws.G_right::Matrix{T}
    jump_eigenbases = ws.jump_eigenbases::Vector{Matrix{T}}
    jump_hermitian = ws.jump_hermitian::Vector{Bool}
    prefactor = (ws.oft_domain_prefactor::Float64) * (ws.gamma_norm_factor::Float64)
    energy_labels = ws.energy_labels::Vector{Float64}

    CT = one(T)
    ZT = zero(T)

    if include_coherent
        BLAS.gemm!('N', 'N', CT, G_left, rho, ZT, sc.rho_out)
        BLAS.gemm!('N', 'N', CT, rho, G_right, CT, sc.rho_out)
    else
        neg_R = sc.sandwich_tmp
        @. neg_R = G_left + G_right
        half = T(0.5)
        BLAS.gemm!('N', 'N', half, neg_R, rho, ZT, sc.rho_out)
        BLAS.gemm!('N', 'N', half, rho, neg_R, CT, sc.rho_out)
    end

    if Threads.nthreads() > 1 && length(energy_labels) >= OMEGA_THREAD_THRESHOLD
        return _apply_lindbladian_threaded_timetrot!(
            sc, rho, jump_eigenbases, jump_hermitian,
            nufft_data, nufft_idx, energy_labels, config, prefactor; adjoint=false)
    end

    for (k, eigenbasis) in enumerate(jump_eigenbases)
        is_herm = jump_hermitian[k]
        if is_herm
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)

                nufft_prefactor_matrix = @view nufft_data[:, :, nufft_idx[w]]
                @. sc.jump_oft = eigenbasis * nufft_prefactor_matrix

                scalar_w = prefactor * pick_transition(config, w)
                _accumulate_sandwich_scratch!(sc.rho_out, sc.jump_oft, rho, scalar_w, sc.sandwich_tmp, sc.sandwich_out)

                if w > 1e-12
                    scalar_neg = prefactor * pick_transition(config, -w)
                    _accumulate_sandwich_adj_scratch!(sc.rho_out, sc.jump_oft, rho, scalar_neg, sc.sandwich_tmp, sc.sandwich_out)
                end
            end
        else
            for (i, w) in enumerate(energy_labels)
                nufft_prefactor_matrix = @view nufft_data[:, :, i]
                @. sc.jump_oft = eigenbasis * nufft_prefactor_matrix
                scalar_w = prefactor * pick_transition(config, w)
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
    hamiltonian::HamHam;
    include_coherent::Bool = true,
) where {T<:Complex, D<:Union{TimeDomain, TrotterDomain}}
    sc = ws.scratch::KrylovScratch{T}
    _nufft = ws.oft_nufft_prefactors::NUFFTPrefactors{real(T), Array{T, 3}}
    nufft_data = _nufft.data
    nufft_idx = _nufft.energy_to_index
    G_left_adj = ws.G_left_adj::Matrix{T}
    G_right_adj = ws.G_right_adj::Matrix{T}
    jump_eigenbases = ws.jump_eigenbases::Vector{Matrix{T}}
    jump_hermitian = ws.jump_hermitian::Vector{Bool}
    prefactor = (ws.oft_domain_prefactor::Float64) * (ws.gamma_norm_factor::Float64)
    energy_labels = ws.energy_labels::Vector{Float64}

    CT = one(T)
    ZT = zero(T)

    if include_coherent
        BLAS.gemm!('N', 'N', CT, G_left_adj, rho, ZT, sc.rho_out)
        BLAS.gemm!('N', 'N', CT, rho, G_right_adj, CT, sc.rho_out)
    else
        neg_R = sc.sandwich_tmp
        @. neg_R = G_left_adj + G_right_adj
        half = T(0.5)
        BLAS.gemm!('N', 'N', half, neg_R, rho, ZT, sc.rho_out)
        BLAS.gemm!('N', 'N', half, rho, neg_R, CT, sc.rho_out)
    end

    if Threads.nthreads() > 1 && length(energy_labels) >= OMEGA_THREAD_THRESHOLD
        return _apply_lindbladian_threaded_timetrot!(
            sc, rho, jump_eigenbases, jump_hermitian,
            nufft_data, nufft_idx, energy_labels, config, prefactor; adjoint=true)
    end

    for (k, eigenbasis) in enumerate(jump_eigenbases)
        is_herm = jump_hermitian[k]
        if is_herm
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)

                nufft_prefactor_matrix = @view nufft_data[:, :, nufft_idx[w]]
                @. sc.jump_oft = eigenbasis * nufft_prefactor_matrix

                scalar_w = prefactor * pick_transition(config, w)
                _accumulate_sandwich_adj_scratch!(sc.rho_out, sc.jump_oft, rho, scalar_w, sc.sandwich_tmp, sc.sandwich_out)

                if w > 1e-12
                    scalar_neg = prefactor * pick_transition(config, -w)
                    _accumulate_sandwich_scratch!(sc.rho_out, sc.jump_oft, rho, scalar_neg, sc.sandwich_tmp, sc.sandwich_out)
                end
            end
        else
            for (i, w) in enumerate(energy_labels)
                nufft_prefactor_matrix = @view nufft_data[:, :, i]
                @. sc.jump_oft = eigenbasis * nufft_prefactor_matrix
                scalar_w = prefactor * pick_transition(config, w)
                _accumulate_sandwich_adj_scratch!(sc.rho_out, sc.jump_oft, rho, scalar_w, sc.sandwich_tmp, sc.sandwich_out)
            end
        end
    end

    return sc.rho_out
end

# ---------------------------------------------------------------------------
# Threaded ω-loop variants (qf-in3) — mirror channel pattern in jump_workers.jl
# ---------------------------------------------------------------------------
#
# Build a flat work-list of (jump_idx, label_idx) pairs, partition across
# threads, and accumulate per-thread `rho_out` chunks into a final reduction
# that is added to the (already-coherent-populated) `sc.rho_out`.
#
# Reuses `OMEGA_THREAD_THRESHOLD` and `_partition_range` from `jump_workers.jl`
# (loaded earlier; both are package-internal).

# Build flat (k, li) work list honoring the hermitian-fold convention:
# Hermitian: only li with w_raw <= 1e-12 are queued (non-positive labels).
# Non-Hermitian: all li are queued.
function _populate_lindblad_work_list!(
    work::Vector{Tuple{Int, Int}},
    jump_hermitian::Vector{Bool},
    energy_labels::AbstractVector{<:Real},
)
    n_jumps = length(jump_hermitian)
    n_labels = length(energy_labels)
    empty!(work)
    sizehint!(work, n_jumps * n_labels)
    @inbounds for k in 1:n_jumps
        is_herm = jump_hermitian[k]
        for li in 1:n_labels
            if is_herm && energy_labels[li] > 1e-12
                continue
            end
            push!(work, (k, li))
        end
    end
    return work
end

# --- EnergyDomain threaded variant ---

function _apply_lindbladian_threaded_energy!(
    sc::KrylovScratch{T},
    rho::Matrix{T},
    jump_eigenbases::Vector{Matrix{T}},
    jump_hermitian::Vector{Bool},
    bohr_freqs::AbstractMatrix{<:Real},
    energy_labels::Vector{Float64},
    config::Config{Lindbladian, EnergyDomain},
    prefactor::Float64,
    inv_4sigma2::Float64;
    adjoint::Bool,
) where {T<:Complex}
    # `work` is the scratch's pre-allocated buffer; `_populate_…!` does
    # `empty!` + `push!` which is zero-alloc when the buffer is large enough
    # (the Workspace constructor sized it for the production label set).
    work = sc.work_list
    _populate_lindblad_work_list!(work, jump_hermitian, energy_labels)
    n_work = length(work)
    n_work == 0 && return sc.rho_out

    pool = sc.task_scratches
    nt = min(Threads.nthreads(), n_work, length(pool))
    chunks = _partition_range(1:n_work, nt)

    old_blas = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        @sync for (idx, chunk) in enumerate(chunks)
            Threads.@spawn _apply_lindbladian_chunk_energy!(
                pool[idx], rho, jump_eigenbases, jump_hermitian,
                bohr_freqs, energy_labels, work, chunk, config,
                prefactor, inv_4sigma2; adjoint=adjoint)
        end
    finally
        BLAS.set_num_threads(old_blas)
    end

    @inbounds for idx in 1:length(chunks)
        sc.rho_out .+= pool[idx].rho_out
    end

    return sc.rho_out
end

function _apply_lindbladian_chunk_energy!(
    task_sc::KrylovScratch{T},
    rho::Matrix{T},
    jump_eigenbases::Vector{Matrix{T}},
    jump_hermitian::Vector{Bool},
    bohr_freqs::AbstractMatrix{<:Real},
    energy_labels::Vector{Float64},
    work::Vector{Tuple{Int, Int}},
    chunk::UnitRange{Int},
    config::Config{Lindbladian, EnergyDomain},
    prefactor::Float64,
    inv_4sigma2::Float64;
    adjoint::Bool,
) where {T<:Complex}
    fill!(task_sc.rho_out, 0)

    @inbounds for w_idx in chunk
        (k, li) = work[w_idx]
        eigenbasis = jump_eigenbases[k]
        is_herm = jump_hermitian[k]

        w_raw = energy_labels[li]
        # Hermitian fold: only `w_raw <= 1e-12` is queued, OFT and rate use
        # `w = |w_raw|` (matches serial). Non-Hermitian: `w = w_raw` directly,
        # OFT and rate take the signed value.
        w = is_herm ? abs(w_raw) : w_raw

        oft!(task_sc.jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)

        scalar_w = prefactor * pick_transition(config, w)
        if adjoint
            _accumulate_sandwich_adj_scratch!(task_sc.rho_out, task_sc.jump_oft, rho, scalar_w,
                                              task_sc.sandwich_tmp, task_sc.sandwich_out)
        else
            _accumulate_sandwich_scratch!(task_sc.rho_out, task_sc.jump_oft, rho, scalar_w,
                                          task_sc.sandwich_tmp, task_sc.sandwich_out)
        end

        if is_herm && w > 1e-12
            scalar_neg = prefactor * pick_transition(config, -w)
            if adjoint
                _accumulate_sandwich_scratch!(task_sc.rho_out, task_sc.jump_oft, rho, scalar_neg,
                                              task_sc.sandwich_tmp, task_sc.sandwich_out)
            else
                _accumulate_sandwich_adj_scratch!(task_sc.rho_out, task_sc.jump_oft, rho, scalar_neg,
                                                  task_sc.sandwich_tmp, task_sc.sandwich_out)
            end
        end
    end

    return nothing
end

# --- TimeDomain / TrotterDomain threaded variant ---

function _apply_lindbladian_threaded_timetrot!(
    sc::KrylovScratch{T},
    rho::Matrix{T},
    jump_eigenbases::Vector{Matrix{T}},
    jump_hermitian::Vector{Bool},
    nufft_data::AbstractArray{T, 3},
    nufft_idx::AbstractDict,
    energy_labels::Vector{Float64},
    config::Config{Lindbladian, D},
    prefactor::Float64;
    adjoint::Bool,
) where {T<:Complex, D<:Union{TimeDomain, TrotterDomain}}
    work = sc.work_list
    _populate_lindblad_work_list!(work, jump_hermitian, energy_labels)
    n_work = length(work)
    n_work == 0 && return sc.rho_out

    pool = sc.task_scratches
    nt = min(Threads.nthreads(), n_work, length(pool))
    chunks = _partition_range(1:n_work, nt)

    old_blas = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        @sync for (idx, chunk) in enumerate(chunks)
            Threads.@spawn _apply_lindbladian_chunk_timetrot!(
                pool[idx], rho, jump_eigenbases, jump_hermitian,
                nufft_data, nufft_idx, energy_labels, work, chunk, config,
                prefactor; adjoint=adjoint)
        end
    finally
        BLAS.set_num_threads(old_blas)
    end

    @inbounds for idx in 1:length(chunks)
        sc.rho_out .+= pool[idx].rho_out
    end

    return sc.rho_out
end

function _apply_lindbladian_chunk_timetrot!(
    task_sc::KrylovScratch{T},
    rho::Matrix{T},
    jump_eigenbases::Vector{Matrix{T}},
    jump_hermitian::Vector{Bool},
    nufft_data::AbstractArray{T, 3},
    nufft_idx::AbstractDict,
    energy_labels::Vector{Float64},
    work::Vector{Tuple{Int, Int}},
    chunk::UnitRange{Int},
    config::Config{Lindbladian, D},
    prefactor::Float64;
    adjoint::Bool,
) where {T<:Complex, D<:Union{TimeDomain, TrotterDomain}}
    fill!(task_sc.rho_out, 0)

    @inbounds for w_idx in chunk
        (k, li) = work[w_idx]
        eigenbasis = jump_eigenbases[k]
        is_herm = jump_hermitian[k]

        w_raw = energy_labels[li]
        # Hermitian fold: only `w_raw <= 1e-12` queued; rate uses `|w_raw|`,
        # NUFFT prefactor index found via `nufft_idx[|w_raw|]`. Non-Hermitian:
        # rate uses signed `w_raw`; prefactor index is the label index `li`.
        w = is_herm ? abs(w_raw) : w_raw
        prefactor_idx = is_herm ? nufft_idx[w] : li
        nufft_prefactor_matrix = @view nufft_data[:, :, prefactor_idx]
        @. task_sc.jump_oft = eigenbasis * nufft_prefactor_matrix

        scalar_w = prefactor * pick_transition(config, w)
        if adjoint
            _accumulate_sandwich_adj_scratch!(task_sc.rho_out, task_sc.jump_oft, rho, scalar_w,
                                              task_sc.sandwich_tmp, task_sc.sandwich_out)
        else
            _accumulate_sandwich_scratch!(task_sc.rho_out, task_sc.jump_oft, rho, scalar_w,
                                          task_sc.sandwich_tmp, task_sc.sandwich_out)
        end

        if is_herm && w > 1e-12
            scalar_neg = prefactor * pick_transition(config, -w)
            if adjoint
                _accumulate_sandwich_scratch!(task_sc.rho_out, task_sc.jump_oft, rho, scalar_neg,
                                              task_sc.sandwich_tmp, task_sc.sandwich_out)
            else
                _accumulate_sandwich_adj_scratch!(task_sc.rho_out, task_sc.jump_oft, rho, scalar_neg,
                                                  task_sc.sandwich_tmp, task_sc.sandwich_out)
            end
        end
    end

    return nothing
end
