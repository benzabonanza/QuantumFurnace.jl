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

"""
    _accumulate_dissipator!(out, L_op, rho, scalar, ws) -> nothing

Accumulate `scalar * D_L(rho)` into `out` in-place:

    D_L(rho) = L * rho * L' - 0.5 * (L'L * rho + rho * L'L)

Uses `ws.tmp1`, `ws.tmp2`, `ws.LdagL` as scratch. Uses BLAS.gemm! for
Adjoint arguments to avoid boxing allocations from `mul!` with Adjoint wrappers.
"""
function _accumulate_dissipator!(
    out::Matrix{T},
    L_op::Matrix{T},
    rho::Matrix{T},
    scalar::Real,
    ws,  # KrylovWorkspace -- uses tmp1, tmp2, LdagL
) where {T<:Complex}
    CT = one(T)
    ZT = zero(T)

    # L'L -> ws.LdagL
    BLAS.gemm!('C', 'N', CT, L_op, L_op, ZT, ws.LdagL)

    # Term 1: scalar * L * rho * L'
    BLAS.gemm!('N', 'N', CT, L_op, rho, ZT, ws.tmp1)      # tmp1 = L * rho
    BLAS.gemm!('N', 'C', CT, ws.tmp1, L_op, ZT, ws.tmp2)   # tmp2 = L * rho * L'
    BLAS.axpy!(T(scalar), ws.tmp2, out)

    # Term 2: -0.5 * scalar * L'L * rho
    BLAS.gemm!('N', 'N', CT, ws.LdagL, rho, ZT, ws.tmp1)   # tmp1 = L'L * rho
    BLAS.axpy!(T(-0.5 * scalar), ws.tmp1, out)

    # Term 3: -0.5 * scalar * rho * L'L
    BLAS.gemm!('N', 'N', CT, rho, ws.LdagL, ZT, ws.tmp1)   # tmp1 = rho * L'L
    BLAS.axpy!(T(-0.5 * scalar), ws.tmp1, out)

    return nothing
end

"""
    _accumulate_dissipator_adj_L!(out, L_op, rho, scalar, ws) -> nothing

Accumulate `scalar * D_{L'}(rho)` into `out` in-place, where the operator is L' (adjoint of L):

    D_{L'}(rho) = L' * rho * L - 0.5 * ((L')' L' * rho + rho * (L')' L')
                = L' * rho * L - 0.5 * (L L' * rho + rho * L L')

Used for the negative-frequency partner in the Hermitian half-grid optimization.
"""
function _accumulate_dissipator_adj_L!(
    out::Matrix{T},
    L_op::Matrix{T},
    rho::Matrix{T},
    scalar::Real,
    ws,
) where {T<:Complex}
    CT = one(T)
    ZT = zero(T)

    # (L')' L' = L L' -> ws.LdagL
    BLAS.gemm!('N', 'C', CT, L_op, L_op, ZT, ws.LdagL)

    # Term 1: scalar * L' * rho * L
    BLAS.gemm!('C', 'N', CT, L_op, rho, ZT, ws.tmp1)       # tmp1 = L' * rho
    BLAS.gemm!('N', 'N', CT, ws.tmp1, L_op, ZT, ws.tmp2)   # tmp2 = L' * rho * L
    BLAS.axpy!(T(scalar), ws.tmp2, out)

    # Term 2: -0.5 * scalar * L L' * rho
    BLAS.gemm!('N', 'N', CT, ws.LdagL, rho, ZT, ws.tmp1)
    BLAS.axpy!(T(-0.5 * scalar), ws.tmp1, out)

    # Term 3: -0.5 * scalar * rho * L L'
    BLAS.gemm!('N', 'N', CT, rho, ws.LdagL, ZT, ws.tmp1)
    BLAS.axpy!(T(-0.5 * scalar), ws.tmp1, out)

    return nothing
end

"""
    _accumulate_adjoint_dissipator!(out, L_op, rho, scalar, ws) -> nothing

Accumulate `scalar * D_L*(rho)` (Hilbert-Schmidt adjoint dissipator) into `out` in-place:

    D_L*(rho) = L' * rho * L - 0.5 * (L'L * rho + rho * L'L)

The sandwich term changes from `L rho L'` (forward) to `L' rho L` (adjoint),
but the anticommutator `{L'L, rho}` stays the same.
"""
function _accumulate_adjoint_dissipator!(
    out::Matrix{T},
    L_op::Matrix{T},
    rho::Matrix{T},
    scalar::Real,
    ws,
) where {T<:Complex}
    CT = one(T)
    ZT = zero(T)

    # L'L -> ws.LdagL (same as forward -- anticommutator is unchanged)
    BLAS.gemm!('C', 'N', CT, L_op, L_op, ZT, ws.LdagL)

    # Term 1: scalar * L' * rho * L  (adjoint sandwich: L' rho L instead of L rho L')
    BLAS.gemm!('C', 'N', CT, L_op, rho, ZT, ws.tmp1)       # tmp1 = L' * rho
    BLAS.gemm!('N', 'N', CT, ws.tmp1, L_op, ZT, ws.tmp2)   # tmp2 = L' * rho * L
    BLAS.axpy!(T(scalar), ws.tmp2, out)

    # Term 2: -0.5 * scalar * L'L * rho  (same as forward)
    BLAS.gemm!('N', 'N', CT, ws.LdagL, rho, ZT, ws.tmp1)
    BLAS.axpy!(T(-0.5 * scalar), ws.tmp1, out)

    # Term 3: -0.5 * scalar * rho * L'L  (same as forward)
    BLAS.gemm!('N', 'N', CT, rho, ws.LdagL, ZT, ws.tmp1)
    BLAS.axpy!(T(-0.5 * scalar), ws.tmp1, out)

    return nothing
end

"""
    _accumulate_adjoint_dissipator_adj_L!(out, L_op, rho, scalar, ws) -> nothing

Accumulate `scalar * D_{L'}*(rho)` into `out` in-place, where the original operator is L'
and we take the Hilbert-Schmidt adjoint of the resulting superoperator:

    D_{L'}*(rho) = (L')' * rho * L' - 0.5 * ((L')' L' * rho + rho * (L')' L')
                 = L * rho * L' - 0.5 * (L L' * rho + rho * L L')

Used for the negative-frequency partner in the adjoint Hermitian half-grid.
"""
function _accumulate_adjoint_dissipator_adj_L!(
    out::Matrix{T},
    L_op::Matrix{T},
    rho::Matrix{T},
    scalar::Real,
    ws,
) where {T<:Complex}
    CT = one(T)
    ZT = zero(T)

    # (L')' L' = L L' -> ws.LdagL
    BLAS.gemm!('N', 'C', CT, L_op, L_op, ZT, ws.LdagL)

    # Term 1: scalar * L * rho * L'  (adjoint of D_{L'} sandwich)
    BLAS.gemm!('N', 'N', CT, L_op, rho, ZT, ws.tmp1)       # tmp1 = L * rho
    BLAS.gemm!('N', 'C', CT, ws.tmp1, L_op, ZT, ws.tmp2)   # tmp2 = L * rho * L'
    BLAS.axpy!(T(scalar), ws.tmp2, out)

    # Term 2: -0.5 * scalar * L L' * rho
    BLAS.gemm!('N', 'N', CT, ws.LdagL, rho, ZT, ws.tmp1)
    BLAS.axpy!(T(-0.5 * scalar), ws.tmp1, out)

    # Term 3: -0.5 * scalar * rho * L L'
    BLAS.gemm!('N', 'N', CT, rho, ws.LdagL, ZT, ws.tmp1)
    BLAS.axpy!(T(-0.5 * scalar), ws.tmp1, out)

    return nothing
end

"""
    apply_lindbladian!(ws, rho, config, hamiltonian) -> ws.rho_out

Compute L(rho) for `EnergyDomain` configs, storing the result in `ws.rho_out`.

Mirrors `_jump_contribution!` for `AbstractLiouvConfig{EnergyDomain}` in
`jump_workers.jl` line-by-line, replacing vectorized superoperator accumulation
with direct matrix-on-matrix operations via `_accumulate_dissipator!`.

Uses concrete-typed `ws.jump_eigenbases` and `ws.jump_hermitian` to avoid
boxing allocations from `JumpOp`'s abstract `in_eigenbasis::Matrix{<:Complex}` field.

# Arguments
- `ws::KrylovWorkspace`: pre-allocated workspace (scratch matrices + precomputed data)
- `rho::Matrix{<:Complex}`: input density matrix (dim x dim)
- `config::AbstractLiouvConfig{EnergyDomain}`: Lindbladian configuration
- `hamiltonian::HamHam`: Hamiltonian with eigenbasis data

# Returns
`ws.rho_out` -- the output matrix (dim x dim), overwritten in-place each call.
"""
function apply_lindbladian!(
    ws::KrylovWorkspace{T},
    rho::Matrix{T},
    config::AbstractLiouvConfig{EnergyDomain},
    hamiltonian::HamHam,
) where {T<:Complex}
    (; transition, gamma_norm_factor, energy_labels) = ws.precomputed_data
    bohr_freqs = hamiltonian.bohr_freqs
    inv_4sigma2 = 1.0 / (4 * config.sigma^2)

    # Zero output accumulator (CRITICAL: must be zeroed each call)
    fill!(ws.rho_out, 0)

    # Coherent term: -i[B, rho] = -i*B*rho + i*rho*B
    B = ws.B_total
    if B !== nothing
        mul!(ws.tmp1, B, rho)
        @. ws.rho_out += -1im * ws.tmp1
        mul!(ws.tmp1, rho, B)
        @. ws.rho_out += 1im * ws.tmp1
    end

    # Dissipator: sum over jumps, sum over energy labels
    prefactor = (config.w0 / (config.sigma * sqrt(2 * pi))) * gamma_norm_factor

    for (k, eigenbasis) in enumerate(ws.jump_eigenbases)
        is_herm = ws.jump_hermitian[k]
        if is_herm
            for w_raw in energy_labels
                # Iterate only half-grid (w <= 0) and mirror manually
                w_raw > 1e-12 && continue
                w = abs(w_raw)

                _krylov_oft!(ws.jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)

                scalar_w = prefactor * transition(w)
                _accumulate_dissipator!(ws.rho_out, ws.jump_oft, rho, scalar_w, ws)

                if w > 1e-12
                    # Negative-frequency partner: L = A(omega)'
                    scalar_neg = prefactor * transition(-w)
                    _accumulate_dissipator_adj_L!(ws.rho_out, ws.jump_oft, rho, scalar_neg, ws)
                end
            end
        else
            for w in energy_labels
                _krylov_oft!(ws.jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)
                scalar_w = prefactor * transition(w)
                _accumulate_dissipator!(ws.rho_out, ws.jump_oft, rho, scalar_w, ws)
            end
        end
    end

    return ws.rho_out
end

"""
    apply_adjoint_lindbladian!(ws, rho, config, hamiltonian) -> ws.rho_out

Compute L*(rho) (Hilbert-Schmidt adjoint) for `EnergyDomain` configs.

Differences from forward `apply_lindbladian!`:
1. Coherent term sign flip: `+i[B, rho]` instead of `-i[B, rho]`
2. Dissipator sandwich swap: `L rho L'` becomes `L' rho L` (anticommutator unchanged)

Uses `_accumulate_adjoint_dissipator!` which correctly preserves the `{L'L, rho}`
anticommutator while swapping only the sandwich term.

Shares the same `KrylovWorkspace` as the forward function.

# Arguments
Same as `apply_lindbladian!`.

# Returns
`ws.rho_out` -- the adjoint Lindbladian action on rho.
"""
function apply_adjoint_lindbladian!(
    ws::KrylovWorkspace{T},
    rho::Matrix{T},
    config::AbstractLiouvConfig{EnergyDomain},
    hamiltonian::HamHam,
) where {T<:Complex}
    (; transition, gamma_norm_factor, energy_labels) = ws.precomputed_data
    bohr_freqs = hamiltonian.bohr_freqs
    inv_4sigma2 = 1.0 / (4 * config.sigma^2)

    # Zero output accumulator
    fill!(ws.rho_out, 0)

    # Adjoint coherent term: +i[B, rho] = +i*B*rho - i*rho*B (sign flip vs forward)
    B = ws.B_total
    if B !== nothing
        mul!(ws.tmp1, B, rho)
        @. ws.rho_out += 1im * ws.tmp1      # +i instead of -i
        mul!(ws.tmp1, rho, B)
        @. ws.rho_out -= 1im * ws.tmp1      # -i instead of +i
    end

    # Adjoint dissipator: D_L*(rho) = L' rho L - 0.5 {L'L, rho}
    # Same L_op as forward, but uses _accumulate_adjoint_dissipator! which
    # swaps the sandwich (L rho L' -> L' rho L) while keeping anticommutator {L'L, rho}.
    prefactor = (config.w0 / (config.sigma * sqrt(2 * pi))) * gamma_norm_factor

    for (k, eigenbasis) in enumerate(ws.jump_eigenbases)
        is_herm = ws.jump_hermitian[k]
        if is_herm
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)

                _krylov_oft!(ws.jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)

                scalar_w = prefactor * transition(w)
                _accumulate_adjoint_dissipator!(ws.rho_out, ws.jump_oft, rho, scalar_w, ws)

                if w > 1e-12
                    scalar_neg = prefactor * transition(-w)
                    # Negative-frequency partner: original L = A(omega)'
                    _accumulate_adjoint_dissipator_adj_L!(ws.rho_out, ws.jump_oft, rho, scalar_neg, ws)
                end
            end
        else
            for w in energy_labels
                _krylov_oft!(ws.jump_oft, eigenbasis, bohr_freqs, w, inv_4sigma2)
                scalar_w = prefactor * transition(w)
                _accumulate_adjoint_dissipator!(ws.rho_out, ws.jump_oft, rho, scalar_w, ws)
            end
        end
    end

    return ws.rho_out
end

"""
    apply_lindbladian!(ws, rho, config, hamiltonian) -> ws.rho_out

Compute L(rho) for `TimeDomain` and `TrotterDomain` configs, storing the result in `ws.rho_out`.

Mirrors `_jump_contribution!` for `AbstractLiouvConfig{Union{TimeDomain, TrotterDomain}}` in
`jump_workers.jl`, replacing vectorized superoperator accumulation with direct matrix-on-matrix
operations via `_accumulate_dissipator!`.

Key differences from `EnergyDomain`:
- OFT computation uses NUFFT prefactors (`_prefactor_view`) instead of Gaussian filter (`_krylov_oft!`)
- Scalar prefactor: `w0 * t0^2 * sigma * sqrt(2/pi) / (2pi) * gamma_norm_factor`
  (from `jump_workers.jl:109`)

Uses concrete-typed `ws.jump_eigenbases` and `ws.jump_hermitian` to avoid
boxing allocations from `JumpOp`'s abstract `in_eigenbasis::Matrix{<:Complex}` field.

# Arguments
- `ws::KrylovWorkspace`: pre-allocated workspace (scratch matrices + precomputed data)
- `rho::Matrix{<:Complex}`: input density matrix (dim x dim)
- `config::AbstractLiouvConfig{D}`: Lindbladian configuration (TimeDomain or TrotterDomain)
- `hamiltonian::HamHam`: Hamiltonian with eigenbasis data

# Returns
`ws.rho_out` -- the output matrix (dim x dim), overwritten in-place each call.
"""
function apply_lindbladian!(
    ws::KrylovWorkspace{T},
    rho::Matrix{T},
    config::AbstractLiouvConfig{D},
    hamiltonian::HamHam,
) where {T<:Complex, D<:Union{TimeDomain, TrotterDomain}}
    (; transition, gamma_norm_factor, energy_labels, oft_nufft_prefactors) = ws.precomputed_data

    # Zero output accumulator (CRITICAL: must be zeroed each call)
    fill!(ws.rho_out, 0)

    # Coherent term: -i[B, rho] = -i*B*rho + i*rho*B (identical to EnergyDomain)
    B = ws.B_total
    if B !== nothing
        mul!(ws.tmp1, B, rho)
        @. ws.rho_out += -1im * ws.tmp1
        mul!(ws.tmp1, rho, B)
        @. ws.rho_out += 1im * ws.tmp1
    end

    # Dissipator: sum over jumps, sum over energy labels
    # Time/Trotter prefactor from jump_workers.jl:109
    prefactor = config.w0 * config.t0^2 * (config.sigma * sqrt(2 / pi)) / (2 * pi) * gamma_norm_factor

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
                _accumulate_dissipator!(ws.rho_out, ws.jump_oft, rho, scalar_w, ws)

                if w > 1e-12
                    # Negative-frequency partner: L = A(omega)'
                    scalar_neg = prefactor * transition(-w)
                    _accumulate_dissipator_adj_L!(ws.rho_out, ws.jump_oft, rho, scalar_neg, ws)
                end
            end
        else
            for w in energy_labels
                nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)
                @. ws.jump_oft = eigenbasis * nufft_prefactor_matrix
                scalar_w = prefactor * transition(w)
                _accumulate_dissipator!(ws.rho_out, ws.jump_oft, rho, scalar_w, ws)
            end
        end
    end

    return ws.rho_out
end

"""
    apply_adjoint_lindbladian!(ws, rho, config, hamiltonian) -> ws.rho_out

Compute L*(rho) (Hilbert-Schmidt adjoint) for `TimeDomain` and `TrotterDomain` configs.

Mirrors `_jump_contribution!` for `AbstractLiouvConfig{Union{TimeDomain, TrotterDomain}}`
with the adjoint modifications:
1. Coherent term sign flip: `+i[B, rho]` instead of `-i[B, rho]`
2. Dissipator sandwich swap: `L rho L'` becomes `L' rho L` (anticommutator unchanged)

Same NUFFT prefactor computation and scalar prefactor as the forward method
(no prefactor changes for adjoint, per CONTEXT.md).

# Arguments
Same as `apply_lindbladian!` for `Union{TimeDomain, TrotterDomain}`.

# Returns
`ws.rho_out` -- the adjoint Lindbladian action on rho.
"""
function apply_adjoint_lindbladian!(
    ws::KrylovWorkspace{T},
    rho::Matrix{T},
    config::AbstractLiouvConfig{D},
    hamiltonian::HamHam,
) where {T<:Complex, D<:Union{TimeDomain, TrotterDomain}}
    (; transition, gamma_norm_factor, energy_labels, oft_nufft_prefactors) = ws.precomputed_data

    # Zero output accumulator
    fill!(ws.rho_out, 0)

    # Adjoint coherent term: +i[B, rho] = +i*B*rho - i*rho*B (sign flip vs forward)
    B = ws.B_total
    if B !== nothing
        mul!(ws.tmp1, B, rho)
        @. ws.rho_out += 1im * ws.tmp1      # +i instead of -i
        mul!(ws.tmp1, rho, B)
        @. ws.rho_out -= 1im * ws.tmp1      # -i instead of +i
    end

    # Adjoint dissipator: D_L*(rho) = L' rho L - 0.5 {L'L, rho}
    # Same prefactor as forward (no changes for adjoint)
    prefactor = config.w0 * config.t0^2 * (config.sigma * sqrt(2 / pi)) / (2 * pi) * gamma_norm_factor

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
                _accumulate_adjoint_dissipator!(ws.rho_out, ws.jump_oft, rho, scalar_w, ws)

                if w > 1e-12
                    scalar_neg = prefactor * transition(-w)
                    # Negative-frequency partner: original L = A(omega)'
                    _accumulate_adjoint_dissipator_adj_L!(ws.rho_out, ws.jump_oft, rho, scalar_neg, ws)
                end
            end
        else
            for w in energy_labels
                nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)
                @. ws.jump_oft = eigenbasis * nufft_prefactor_matrix
                scalar_w = prefactor * transition(w)
                _accumulate_adjoint_dissipator!(ws.rho_out, ws.jump_oft, rho, scalar_w, ws)
            end
        end
    end

    return ws.rho_out
end
