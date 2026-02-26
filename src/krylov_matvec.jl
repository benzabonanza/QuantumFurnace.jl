"""
    _krylov_oft!(out, eigenbasis, bohr_freqs, energy, inv_4sigma2) -> nothing

Compute the Optimal Fourier Transform A(omega) in-place, without accessing
any abstractly-typed `JumpOp` fields.

    out[i,j] = eigenbasis[i,j] * exp(-(energy - bohr_freqs[i,j])^2 * inv_4sigma2)

This is the zero-allocation replacement for `oft!` in the Krylov hot path.
"""
@inline function _krylov_oft!(
    out::Matrix{T},
    eigenbasis::Matrix{T},
    bohr_freqs::Matrix{<:Real},
    energy::Real,
    inv_4sigma2::Real,
) where {T<:Complex}
    @. out = eigenbasis * exp(-(energy - bohr_freqs)^2 * inv_4sigma2)
    return nothing
end

# ---------------------------------------------------------------------------
# Sandwich-only helpers (Phase 32 optimization)
# These are stripped-down versions of the dissipator functions with the
# L'L and anticommutator terms removed (absorbed into precomputed G_left/G_right).
# Each performs only 2 GEMMs per call (down from 5 in the full dissipator).
# ---------------------------------------------------------------------------

"""
    _accumulate_sandwich!(out, L_op, rho, scalar, ws) -> nothing

Accumulate `scalar * L * rho * L'` into `out`. This is the sandwich-only
part of the dissipator, used after the anticommutator has been absorbed into
precomputed G_left/G_right (Phase 32 optimization).

Uses `ws.tmp1`, `ws.LdagL` as scratch.
"""
function _accumulate_sandwich!(
    out::Matrix{T},
    L_op::Matrix{T},
    rho::Matrix{T},
    scalar::Real,
    ws,
) where {T<:Complex}
    CT = one(T)
    ZT = zero(T)
    BLAS.gemm!('N', 'N', CT, L_op, rho, ZT, ws.tmp1)           # tmp1 = L * rho
    BLAS.gemm!('N', 'C', CT, ws.tmp1, L_op, ZT, ws.LdagL)      # LdagL = L * rho * L'
    BLAS.axpy!(T(scalar), ws.LdagL, out)
    return nothing
end

"""
    _accumulate_sandwich_adj_L!(out, L_op, rho, scalar, ws) -> nothing

Accumulate `scalar * L' * rho * L` into `out`. This is the sandwich-only
part for the negative-frequency Hermitian partner (M = L'), used after the
anticommutator has been absorbed into G_left/G_right.

Uses `ws.tmp1`, `ws.LdagL` as scratch.
"""
function _accumulate_sandwich_adj_L!(
    out::Matrix{T},
    L_op::Matrix{T},
    rho::Matrix{T},
    scalar::Real,
    ws,
) where {T<:Complex}
    CT = one(T)
    ZT = zero(T)
    BLAS.gemm!('C', 'N', CT, L_op, rho, ZT, ws.tmp1)           # tmp1 = L' * rho
    BLAS.gemm!('N', 'N', CT, ws.tmp1, L_op, ZT, ws.LdagL)      # LdagL = L' * rho * L
    BLAS.axpy!(T(scalar), ws.LdagL, out)
    return nothing
end

"""
    _accumulate_adjoint_sandwich!(out, L_op, rho, scalar, ws) -> nothing

Accumulate `scalar * L' * rho * L` into `out`. This is the HS adjoint
of the forward sandwich `L * rho * L'`, used in apply_adjoint_lindbladian!.

Note: identical computation to _accumulate_sandwich_adj_L! (the HS adjoint of
L*rho*L' equals L'*rho*L, which is the same as the negative-freq partner).
Uses `ws.tmp1`, `ws.LdagL` as scratch.
"""
function _accumulate_adjoint_sandwich!(
    out::Matrix{T},
    L_op::Matrix{T},
    rho::Matrix{T},
    scalar::Real,
    ws,
) where {T<:Complex}
    CT = one(T)
    ZT = zero(T)
    BLAS.gemm!('C', 'N', CT, L_op, rho, ZT, ws.tmp1)           # tmp1 = L' * rho
    BLAS.gemm!('N', 'N', CT, ws.tmp1, L_op, ZT, ws.LdagL)      # LdagL = L' * rho * L
    BLAS.axpy!(T(scalar), ws.LdagL, out)
    return nothing
end

"""
    _accumulate_adjoint_sandwich_adj_L!(out, L_op, rho, scalar, ws) -> nothing

Accumulate `scalar * L * rho * L'` into `out`. This is the HS adjoint of
the negative-frequency sandwich `L' * rho * L`.

Note: identical computation to _accumulate_sandwich! (the HS adjoint of
L'*rho*L equals L*rho*L', which is the same as the positive-freq forward).
Uses `ws.tmp1`, `ws.LdagL` as scratch.
"""
function _accumulate_adjoint_sandwich_adj_L!(
    out::Matrix{T},
    L_op::Matrix{T},
    rho::Matrix{T},
    scalar::Real,
    ws,
) where {T<:Complex}
    CT = one(T)
    ZT = zero(T)
    BLAS.gemm!('N', 'N', CT, L_op, rho, ZT, ws.tmp1)           # tmp1 = L * rho
    BLAS.gemm!('N', 'C', CT, ws.tmp1, L_op, ZT, ws.LdagL)      # LdagL = L * rho * L'
    BLAS.axpy!(T(scalar), ws.LdagL, out)
    return nothing
end

"""
    apply_lindbladian!(ws, rho, config, hamiltonian) -> ws.rho_out

Compute L(rho) for `EnergyDomain` configs, storing the result in `ws.rho_out`.

Uses precomputed G_left/G_right (Phase 32) to absorb the coherent term and
anticommutator into 2 GEMMs, then iterates sandwich-only terms (2 GEMMs each).
Total: 2 + 2N GEMMs (down from 2 + 5N before Phase 32).

Uses concrete-typed `ws.jump_eigenbases` and `ws.jump_hermitian` to avoid
boxing allocations from `JumpOp`'s abstract `in_eigenbasis::Matrix{<:Complex}` field.

# Arguments
- `ws::KrylovWorkspace`: pre-allocated workspace (scratch matrices + precomputed data)
- `rho::Matrix{<:Complex}`: input density matrix (dim x dim)
- `config::Config{Lindbladian, EnergyDomain}`: Lindbladian configuration
- `hamiltonian::HamHam`: Hamiltonian with eigenbasis data

# Returns
`ws.rho_out` -- the output matrix (dim x dim), overwritten in-place each call.
"""
function apply_lindbladian!(
    ws::KrylovWorkspace{T},
    rho::Matrix{T},
    config::Config{Lindbladian, EnergyDomain},
    hamiltonian::HamHam,
) where {T<:Complex}
    (; transition, gamma_norm_factor, energy_labels) = ws.precomputed_data
    bohr_freqs = hamiltonian.bohr_freqs
    inv_4sigma2 = 1.0 / (4 * config.sigma^2)

    CT = one(T)
    ZT = zero(T)

    # Effective Hamiltonian: rho_out = G_left * rho + rho * G_right
    BLAS.gemm!('N', 'N', CT, ws.G_left, rho, ZT, ws.rho_out)     # G_left * rho -> rho_out
    BLAS.gemm!('N', 'N', CT, rho, ws.G_right, CT, ws.rho_out)     # + rho * G_right

    # Sandwich-only loop: sum_i scalar_i * L_i * rho * L_i'
    prefactor = ws.precomputed_data.domain_prefactor * gamma_norm_factor

    for (k, eigenbasis) in enumerate(ws.jump_eigenbases)
        is_herm = ws.jump_hermitian[k]
        if is_herm
            for w_raw in energy_labels
                # Iterate only half-grid (w <= 0) and mirror manually
                w_raw > 1e-12 && continue
                w = abs(w_raw)

                _krylov_oft!(ws.jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)

                scalar_w = prefactor * transition(w)
                _accumulate_sandwich!(ws.rho_out, ws.jump_oft, rho, scalar_w, ws)

                if w > 1e-12
                    # Negative-frequency partner: L = A(omega)'
                    scalar_neg = prefactor * transition(-w)
                    _accumulate_sandwich_adj_L!(ws.rho_out, ws.jump_oft, rho, scalar_neg, ws)
                end
            end
        else
            for w in energy_labels
                _krylov_oft!(ws.jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)
                scalar_w = prefactor * transition(w)
                _accumulate_sandwich!(ws.rho_out, ws.jump_oft, rho, scalar_w, ws)
            end
        end
    end

    return ws.rho_out
end

"""
    apply_adjoint_lindbladian!(ws, rho, config, hamiltonian) -> ws.rho_out

Compute L*(rho) (Hilbert-Schmidt adjoint) for `EnergyDomain` configs.

Uses precomputed G_left_adj/G_right_adj (Phase 32) to absorb the adjoint coherent
term and anticommutator into 2 GEMMs, then iterates adjoint sandwich-only terms.
Total: 2 + 2N GEMMs (down from 2 + 5N before Phase 32).

Shares the same `KrylovWorkspace` as the forward function.

# Arguments
Same as `apply_lindbladian!`.

# Returns
`ws.rho_out` -- the adjoint Lindbladian action on rho.
"""
function apply_adjoint_lindbladian!(
    ws::KrylovWorkspace{T},
    rho::Matrix{T},
    config::Config{Lindbladian, EnergyDomain},
    hamiltonian::HamHam,
) where {T<:Complex}
    (; transition, gamma_norm_factor, energy_labels) = ws.precomputed_data
    bohr_freqs = hamiltonian.bohr_freqs
    inv_4sigma2 = 1.0 / (4 * config.sigma^2)

    CT = one(T)
    ZT = zero(T)

    # Adjoint effective Hamiltonian: uses G_left_adj/G_right_adj (swapped G for adjoint)
    BLAS.gemm!('N', 'N', CT, ws.G_left_adj, rho, ZT, ws.rho_out)
    BLAS.gemm!('N', 'N', CT, rho, ws.G_right_adj, CT, ws.rho_out)

    # Adjoint sandwich-only loop: sum_i scalar_i * L_i' * rho * L_i
    prefactor = ws.precomputed_data.domain_prefactor * gamma_norm_factor

    for (k, eigenbasis) in enumerate(ws.jump_eigenbases)
        is_herm = ws.jump_hermitian[k]
        if is_herm
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)

                _krylov_oft!(ws.jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)

                scalar_w = prefactor * transition(w)
                _accumulate_adjoint_sandwich!(ws.rho_out, ws.jump_oft, rho, scalar_w, ws)

                if w > 1e-12
                    scalar_neg = prefactor * transition(-w)
                    # Negative-frequency partner: original L = A(omega)'
                    _accumulate_adjoint_sandwich_adj_L!(ws.rho_out, ws.jump_oft, rho, scalar_neg, ws)
                end
            end
        else
            for w in energy_labels
                _krylov_oft!(ws.jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)
                scalar_w = prefactor * transition(w)
                _accumulate_adjoint_sandwich!(ws.rho_out, ws.jump_oft, rho, scalar_w, ws)
            end
        end
    end

    return ws.rho_out
end

"""
    _accumulate_sandwich_2op!(out, A, B_dag, rho, scalar, ws) -> nothing

Accumulate `scalar * A * rho * B_dag` into `out`. This is the sandwich-only
part of the two-operator dissipator for BohrDomain, used after the anticommutator
has been absorbed into G_left/G_right.

Uses `ws.tmp1`, `ws.tmp2` as scratch.
"""
function _accumulate_sandwich_2op!(
    out::Matrix{T},
    A::Matrix{T},
    B_dag::Matrix{T},
    rho::Matrix{T},
    scalar::Real,
    ws,
) where {T<:Complex}
    CT = one(T)
    ZT = zero(T)
    BLAS.gemm!('N', 'N', CT, A, rho, ZT, ws.tmp1)           # tmp1 = A * rho
    BLAS.gemm!('N', 'N', CT, ws.tmp1, B_dag, ZT, ws.tmp2)   # tmp2 = A * rho * B_dag
    BLAS.axpy!(T(scalar), ws.tmp2, out)
    return nothing
end

"""
    _accumulate_adjoint_sandwich_2op!(out, A, B_dag, rho, scalar, ws) -> nothing

Accumulate `scalar * A' * rho * B_dag'` into `out`. This is the HS adjoint
of the two-operator sandwich `A * rho * B_dag`, used in apply_adjoint_lindbladian!
for BohrDomain.

Uses `ws.tmp1`, `ws.tmp2` as scratch.
"""
function _accumulate_adjoint_sandwich_2op!(
    out::Matrix{T},
    A::Matrix{T},
    B_dag::Matrix{T},
    rho::Matrix{T},
    scalar::Real,
    ws,
) where {T<:Complex}
    CT = one(T)
    ZT = zero(T)
    BLAS.gemm!('C', 'N', CT, A, rho, ZT, ws.tmp1)           # tmp1 = A' * rho
    BLAS.gemm!('N', 'C', CT, ws.tmp1, B_dag, ZT, ws.tmp2)   # tmp2 = A' * rho * B_dag'
    BLAS.axpy!(T(scalar), ws.tmp2, out)
    return nothing
end

"""
    apply_lindbladian!(ws, rho, config, hamiltonian) -> ws.rho_out

Compute L(rho) for `BohrDomain` configs, storing the result in `ws.rho_out`.

Uses precomputed G_left/G_right (Phase 32) to absorb the coherent term and
anticommutator into 2 GEMMs, then iterates sandwich-only 2op terms (2 GEMMs each).
Total: 2 + 2N GEMMs (down from 2 + 5N before Phase 32).

BohrDomain iterates over Bohr frequency buckets (`hamiltonian.bohr_dict` keys)
using a two-operator sandwich where A and B_dag differ. One dense buffer is
allocated per matvec call for A_nu2_dag (acceptable for Bohr).

# Arguments
- `ws::KrylovWorkspace`: pre-allocated workspace (scratch matrices + precomputed data)
- `rho::Matrix{<:Complex}`: input density matrix (dim x dim)
- `config::Config{Lindbladian, BohrDomain}`: Lindbladian configuration
- `hamiltonian::HamHam`: Hamiltonian with bohr_dict and bohr_freqs

# Returns
`ws.rho_out` -- the output matrix (dim x dim), overwritten in-place each call.
"""
function apply_lindbladian!(
    ws::KrylovWorkspace{T},
    rho::Matrix{T},
    config::Config{Lindbladian, BohrDomain},
    hamiltonian::HamHam,
) where {T<:Complex}
    (; alpha, gamma_norm_factor) = ws.precomputed_data
    dim = size(rho, 1)

    CT = one(T)
    ZT = zero(T)

    # Effective Hamiltonian: rho_out = G_left * rho + rho * G_right
    BLAS.gemm!('N', 'N', CT, ws.G_left, rho, ZT, ws.rho_out)
    BLAS.gemm!('N', 'N', CT, rho, ws.G_right, CT, ws.rho_out)

    # Allocate A_nu2_dag buffer (one allocation per matvec -- acceptable for Bohr)
    A_nu2_dag = zeros(T, dim, dim)

    for (k, eigenbasis) in enumerate(ws.jump_eigenbases)
        for nu_2 in keys(hamiltonian.bohr_dict)
            # alpha_A: dense alpha-weighted eigenbasis
            @. ws.jump_oft = alpha(hamiltonian.bohr_freqs, nu_2) * eigenbasis

            # Build A_nu2_dag via index scatter (avoid sparse() allocation)
            fill!(A_nu2_dag, 0)
            indices = hamiltonian.bohr_dict[nu_2]
            @inbounds for idx in indices
                i = idx[1]; j = idx[2]
                A_nu2_dag[j, i] = conj(eigenbasis[i, j])
            end

            _accumulate_sandwich_2op!(ws.rho_out, ws.jump_oft, A_nu2_dag, rho, gamma_norm_factor, ws)
        end
    end

    return ws.rho_out
end

"""
    apply_adjoint_lindbladian!(ws, rho, config, hamiltonian) -> ws.rho_out

Compute L*(rho) (Hilbert-Schmidt adjoint) for `BohrDomain` configs.

Uses precomputed G_left_adj/G_right_adj (Phase 32) to absorb the adjoint coherent
term and anticommutator into 2 GEMMs, then iterates adjoint sandwich-only 2op terms.
Total: 2 + 2N GEMMs (down from 2 + 5N before Phase 32).

For BohrDomain, G_left_adj/G_right_adj differ from G_right/G_left because
R_total is not Hermitian (adjoint anticommutator uses conj(R_total), not R_total^T).

# Arguments
Same as `apply_lindbladian!` for `BohrDomain`.

# Returns
`ws.rho_out` -- the adjoint Lindbladian action on rho.
"""
function apply_adjoint_lindbladian!(
    ws::KrylovWorkspace{T},
    rho::Matrix{T},
    config::Config{Lindbladian, BohrDomain},
    hamiltonian::HamHam,
) where {T<:Complex}
    (; alpha, gamma_norm_factor) = ws.precomputed_data
    dim = size(rho, 1)

    CT = one(T)
    ZT = zero(T)

    # Adjoint effective Hamiltonian: uses G_left_adj/G_right_adj
    BLAS.gemm!('N', 'N', CT, ws.G_left_adj, rho, ZT, ws.rho_out)
    BLAS.gemm!('N', 'N', CT, rho, ws.G_right_adj, CT, ws.rho_out)

    # Allocate A_nu2_dag buffer (one allocation per matvec -- acceptable for Bohr)
    A_nu2_dag = zeros(T, dim, dim)

    for (k, eigenbasis) in enumerate(ws.jump_eigenbases)
        for nu_2 in keys(hamiltonian.bohr_dict)
            # alpha_A: dense alpha-weighted eigenbasis (same as forward)
            @. ws.jump_oft = alpha(hamiltonian.bohr_freqs, nu_2) * eigenbasis

            # Build A_nu2_dag via index scatter (same as forward)
            fill!(A_nu2_dag, 0)
            indices = hamiltonian.bohr_dict[nu_2]
            @inbounds for idx in indices
                i = idx[1]; j = idx[2]
                A_nu2_dag[j, i] = conj(eigenbasis[i, j])
            end

            # Use dedicated adjoint 2op sandwich helper
            _accumulate_adjoint_sandwich_2op!(ws.rho_out, ws.jump_oft, A_nu2_dag, rho, gamma_norm_factor, ws)
        end
    end

    return ws.rho_out
end

"""
    apply_lindbladian!(ws, rho, config, hamiltonian) -> ws.rho_out

Compute L(rho) for `TimeDomain` and `TrotterDomain` configs, storing the result in `ws.rho_out`.

Uses precomputed G_left/G_right (Phase 32) to absorb the coherent term and
anticommutator into 2 GEMMs, then iterates sandwich-only terms (2 GEMMs each).
Total: 2 + 2N GEMMs (down from 2 + 5N before Phase 32).

Key differences from `EnergyDomain`:
- OFT computation uses NUFFT prefactors (`_prefactor_view`) instead of Gaussian filter (`_krylov_oft!`)
- Scalar prefactor: `w0 * t0^2 * sigma * sqrt(2/pi) / (2pi) * gamma_norm_factor`

Uses concrete-typed `ws.jump_eigenbases` and `ws.jump_hermitian` to avoid
boxing allocations from `JumpOp`'s abstract `in_eigenbasis::Matrix{<:Complex}` field.

# Arguments
- `ws::KrylovWorkspace`: pre-allocated workspace (scratch matrices + precomputed data)
- `rho::Matrix{<:Complex}`: input density matrix (dim x dim)
- `config::Config{Lindbladian, D}`: Lindbladian configuration (TimeDomain or TrotterDomain)
- `hamiltonian::HamHam`: Hamiltonian with eigenbasis data

# Returns
`ws.rho_out` -- the output matrix (dim x dim), overwritten in-place each call.
"""
function apply_lindbladian!(
    ws::KrylovWorkspace{T},
    rho::Matrix{T},
    config::Config{Lindbladian, D},
    hamiltonian::HamHam,
) where {T<:Complex, D<:Union{TimeDomain, TrotterDomain}}
    (; transition, gamma_norm_factor, energy_labels, oft_nufft_prefactors) = ws.precomputed_data

    CT = one(T)
    ZT = zero(T)

    # Effective Hamiltonian: rho_out = G_left * rho + rho * G_right
    BLAS.gemm!('N', 'N', CT, ws.G_left, rho, ZT, ws.rho_out)
    BLAS.gemm!('N', 'N', CT, rho, ws.G_right, CT, ws.rho_out)

    # Sandwich-only loop: sum_i scalar_i * L_i * rho * L_i'
    prefactor = ws.precomputed_data.domain_prefactor * gamma_norm_factor

    for (k, eigenbasis) in enumerate(ws.jump_eigenbases)
        is_herm = ws.jump_hermitian[k]
        if is_herm
            for w_raw in energy_labels
                # Iterate only half-grid (w <= 0) and mirror manually
                w_raw > 1e-12 && continue
                w = abs(w_raw)

                # NUFFT OFT: replace _krylov_oft! with prefactor view multiply
                nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)
                @. ws.jump_oft = eigenbasis * nufft_prefactor_matrix

                scalar_w = prefactor * transition(w)
                _accumulate_sandwich!(ws.rho_out, ws.jump_oft, rho, scalar_w, ws)

                if w > 1e-12
                    # Negative-frequency partner: L = A(omega)'
                    scalar_neg = prefactor * transition(-w)
                    _accumulate_sandwich_adj_L!(ws.rho_out, ws.jump_oft, rho, scalar_neg, ws)
                end
            end
        else
            for w in energy_labels
                nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)
                @. ws.jump_oft = eigenbasis * nufft_prefactor_matrix
                scalar_w = prefactor * transition(w)
                _accumulate_sandwich!(ws.rho_out, ws.jump_oft, rho, scalar_w, ws)
            end
        end
    end

    return ws.rho_out
end

"""
    apply_adjoint_lindbladian!(ws, rho, config, hamiltonian) -> ws.rho_out

Compute L*(rho) (Hilbert-Schmidt adjoint) for `TimeDomain` and `TrotterDomain` configs.

Uses precomputed G_left_adj/G_right_adj (Phase 32) to absorb the adjoint coherent
term and anticommutator into 2 GEMMs, then iterates adjoint sandwich-only terms.
Total: 2 + 2N GEMMs (down from 2 + 5N before Phase 32).

Same NUFFT prefactor computation and scalar prefactor as the forward method.

# Arguments
Same as `apply_lindbladian!` for `Union{TimeDomain, TrotterDomain}`.

# Returns
`ws.rho_out` -- the adjoint Lindbladian action on rho.
"""
function apply_adjoint_lindbladian!(
    ws::KrylovWorkspace{T},
    rho::Matrix{T},
    config::Config{Lindbladian, D},
    hamiltonian::HamHam,
) where {T<:Complex, D<:Union{TimeDomain, TrotterDomain}}
    (; transition, gamma_norm_factor, energy_labels, oft_nufft_prefactors) = ws.precomputed_data

    CT = one(T)
    ZT = zero(T)

    # Adjoint effective Hamiltonian: uses G_left_adj/G_right_adj
    BLAS.gemm!('N', 'N', CT, ws.G_left_adj, rho, ZT, ws.rho_out)
    BLAS.gemm!('N', 'N', CT, rho, ws.G_right_adj, CT, ws.rho_out)

    # Adjoint sandwich-only loop
    prefactor = ws.precomputed_data.domain_prefactor * gamma_norm_factor

    for (k, eigenbasis) in enumerate(ws.jump_eigenbases)
        is_herm = ws.jump_hermitian[k]
        if is_herm
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)

                # NUFFT OFT: replace _krylov_oft! with prefactor view multiply
                nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)
                @. ws.jump_oft = eigenbasis * nufft_prefactor_matrix

                scalar_w = prefactor * transition(w)
                _accumulate_adjoint_sandwich!(ws.rho_out, ws.jump_oft, rho, scalar_w, ws)

                if w > 1e-12
                    scalar_neg = prefactor * transition(-w)
                    # Negative-frequency partner: original L = A(omega)'
                    _accumulate_adjoint_sandwich_adj_L!(ws.rho_out, ws.jump_oft, rho, scalar_neg, ws)
                end
            end
        else
            for w in energy_labels
                nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)
                @. ws.jump_oft = eigenbasis * nufft_prefactor_matrix
                scalar_w = prefactor * transition(w)
                _accumulate_adjoint_sandwich!(ws.rho_out, ws.jump_oft, rho, scalar_w, ws)
            end
        end
    end

    return ws.rho_out
end
