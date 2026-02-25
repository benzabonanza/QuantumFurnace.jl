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

    D_L(rho) = conj(L) * rho * L^T - 0.5 * ((L'L)^T * rho + rho * (L'L)^T)

Matches the dense kron convention exactly:
- Sandwich: `kron(L, conj(L))` un-vectorizes to `conj(L) rho L^T`
- Anticommutator: `kron(L'L, I)` un-vectorizes to `rho * (L'L)^T`, and
  `kron(I, (L'L)^T)` un-vectorizes to `(L'L)^T * rho`

Uses `ws.tmp1`, `ws.tmp2`, `ws.LdagL` as scratch. Uses BLAS.gemm! for
Adjoint/Transpose arguments to avoid boxing allocations from `mul!` with wrappers.
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

    # Term 2: -0.5 * scalar * (L'L)^T * rho
    BLAS.gemm!('T', 'N', CT, ws.LdagL, rho, ZT, ws.tmp1)
    BLAS.axpy!(T(-0.5 * scalar), ws.tmp1, out)

    # Term 3: -0.5 * scalar * rho * (L'L)^T
    BLAS.gemm!('N', 'T', CT, rho, ws.LdagL, ZT, ws.tmp1)
    BLAS.axpy!(T(-0.5 * scalar), ws.tmp1, out)

    # Term 1: scalar * conj(L) * rho * L^T  (LdagL now free as scratch)
    @. ws.tmp2 = conj(L_op)                                    # tmp2 = conj(L)
    BLAS.gemm!('N', 'N', CT, ws.tmp2, rho, ZT, ws.tmp1)        # tmp1 = conj(L) * rho
    BLAS.gemm!('N', 'T', CT, ws.tmp1, L_op, ZT, ws.LdagL)     # LdagL = conj(L) * rho * L^T
    BLAS.axpy!(T(scalar), ws.LdagL, out)

    return nothing
end

"""
    _accumulate_dissipator_adj_L!(out, L_op, rho, scalar, ws) -> nothing

Accumulate `scalar * D_{L'}(rho)` into `out` in-place, where the operator is L' (adjoint of L):

    D_{L'}(rho) = L^T * rho * conj(L) - 0.5 * ((LL')^T * rho + rho * (LL')^T)

Matches kron convention with M = L^H: `conj(M) rho M^T = L^T rho conj(L)`.
Anticommutator uses `(M'M)^T = (LL')^T` matching `kron(LL', I)` and `kron(I, (LL')^T)`.
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

    # Term 2: -0.5 * scalar * (LL')^T * rho
    BLAS.gemm!('T', 'N', CT, ws.LdagL, rho, ZT, ws.tmp1)
    BLAS.axpy!(T(-0.5 * scalar), ws.tmp1, out)

    # Term 3: -0.5 * scalar * rho * (LL')^T
    BLAS.gemm!('N', 'T', CT, rho, ws.LdagL, ZT, ws.tmp1)
    BLAS.axpy!(T(-0.5 * scalar), ws.tmp1, out)

    # Term 1: scalar * L^T * rho * conj(L)  (LdagL now free)
    @. ws.tmp2 = conj(L_op)                                    # tmp2 = conj(L)
    BLAS.gemm!('T', 'N', CT, L_op, rho, ZT, ws.tmp1)           # tmp1 = L^T * rho
    BLAS.gemm!('N', 'N', CT, ws.tmp1, ws.tmp2, ZT, ws.LdagL)  # LdagL = L^T * rho * conj(L)
    BLAS.axpy!(T(scalar), ws.LdagL, out)

    return nothing
end

"""
    _accumulate_adjoint_dissipator!(out, L_op, rho, scalar, ws) -> nothing

Accumulate `scalar * D_L*(rho)` (Hilbert-Schmidt adjoint dissipator) into `out` in-place:

    D_L*(rho) = L^T * rho * conj(L) - 0.5 * ((L'L)^T * rho + rho * (L'L)^T)

HS adjoint of `conj(L) rho L^T` is `L^T rho conj(L)` (X -> X^H, Y -> Y^H).
The anticommutator `{(L'L)^T, rho}` is the same as the forward (HS adjoint of
`M^T * rho` is `conj(M) * rho = M^T * rho` for Hermitian M).
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

    # L'L -> ws.LdagL (same as forward)
    BLAS.gemm!('C', 'N', CT, L_op, L_op, ZT, ws.LdagL)

    # Term 2: -0.5 * scalar * (L'L)^T * rho
    BLAS.gemm!('T', 'N', CT, ws.LdagL, rho, ZT, ws.tmp1)
    BLAS.axpy!(T(-0.5 * scalar), ws.tmp1, out)

    # Term 3: -0.5 * scalar * rho * (L'L)^T
    BLAS.gemm!('N', 'T', CT, rho, ws.LdagL, ZT, ws.tmp1)
    BLAS.axpy!(T(-0.5 * scalar), ws.tmp1, out)

    # Term 1: scalar * L^T * rho * conj(L)  (LdagL now free)
    @. ws.tmp2 = conj(L_op)                                    # tmp2 = conj(L)
    BLAS.gemm!('T', 'N', CT, L_op, rho, ZT, ws.tmp1)           # tmp1 = L^T * rho
    BLAS.gemm!('N', 'N', CT, ws.tmp1, ws.tmp2, ZT, ws.LdagL)  # LdagL = L^T * rho * conj(L)
    BLAS.axpy!(T(scalar), ws.LdagL, out)

    return nothing
end

"""
    _accumulate_adjoint_dissipator_adj_L!(out, L_op, rho, scalar, ws) -> nothing

Accumulate `scalar * D_{L'}*(rho)` into `out` in-place, where the original operator is L'
and we take the Hilbert-Schmidt adjoint of the resulting superoperator:

    D_{L'}*(rho) = conj(L) * rho * L^T - 0.5 * ((LL')^T * rho + rho * (LL')^T)

Forward for M=L^H gives sandwich `L^T rho conj(L)`. HS adjoint: `conj(L) rho L^T`.
Anticommutator uses `(LL')^T` matching the kron convention.
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

    # Term 2: -0.5 * scalar * (LL')^T * rho
    BLAS.gemm!('T', 'N', CT, ws.LdagL, rho, ZT, ws.tmp1)
    BLAS.axpy!(T(-0.5 * scalar), ws.tmp1, out)

    # Term 3: -0.5 * scalar * rho * (LL')^T
    BLAS.gemm!('N', 'T', CT, rho, ws.LdagL, ZT, ws.tmp1)
    BLAS.axpy!(T(-0.5 * scalar), ws.tmp1, out)

    # Term 1: scalar * conj(L) * rho * L^T  (LdagL now free)
    @. ws.tmp2 = conj(L_op)                                    # tmp2 = conj(L)
    BLAS.gemm!('N', 'N', CT, ws.tmp2, rho, ZT, ws.tmp1)        # tmp1 = conj(L) * rho
    BLAS.gemm!('N', 'T', CT, ws.tmp1, L_op, ZT, ws.LdagL)     # LdagL = conj(L) * rho * L^T
    BLAS.axpy!(T(scalar), ws.LdagL, out)

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

Accumulate `scalar * conj(L) * rho * L^T` into `out`. This is the sandwich-only
part of the dissipator, used after the anticommutator has been absorbed into
precomputed G_left/G_right (Phase 32 optimization).

Uses `ws.tmp1`, `ws.tmp2`, `ws.LdagL` as scratch.
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
    @. ws.tmp2 = conj(L_op)                                    # tmp2 = conj(L)
    BLAS.gemm!('N', 'N', CT, ws.tmp2, rho, ZT, ws.tmp1)        # tmp1 = conj(L) * rho
    BLAS.gemm!('N', 'T', CT, ws.tmp1, L_op, ZT, ws.LdagL)     # LdagL = conj(L) * rho * L^T
    BLAS.axpy!(T(scalar), ws.LdagL, out)
    return nothing
end

"""
    _accumulate_sandwich_adj_L!(out, L_op, rho, scalar, ws) -> nothing

Accumulate `scalar * L^T * rho * conj(L)` into `out`. This is the sandwich-only
part for the negative-frequency Hermitian partner (M = L'), used after the
anticommutator has been absorbed into G_left/G_right.

Uses `ws.tmp1`, `ws.tmp2`, `ws.LdagL` as scratch.
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
    @. ws.tmp2 = conj(L_op)                                    # tmp2 = conj(L)
    BLAS.gemm!('T', 'N', CT, L_op, rho, ZT, ws.tmp1)           # tmp1 = L^T * rho
    BLAS.gemm!('N', 'N', CT, ws.tmp1, ws.tmp2, ZT, ws.LdagL)  # LdagL = L^T * rho * conj(L)
    BLAS.axpy!(T(scalar), ws.LdagL, out)
    return nothing
end

"""
    _accumulate_adjoint_sandwich!(out, L_op, rho, scalar, ws) -> nothing

Accumulate `scalar * L^T * rho * conj(L)` into `out`. This is the HS adjoint
of the forward sandwich `conj(L) * rho * L^T`, used in apply_adjoint_lindbladian!.

Note: identical computation to _accumulate_sandwich_adj_L! (the HS adjoint of
conj(L)*rho*L^T equals L^T*rho*conj(L), which is the same as the negative-freq partner).
Uses `ws.tmp1`, `ws.tmp2`, `ws.LdagL` as scratch.
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
    @. ws.tmp2 = conj(L_op)                                    # tmp2 = conj(L)
    BLAS.gemm!('T', 'N', CT, L_op, rho, ZT, ws.tmp1)           # tmp1 = L^T * rho
    BLAS.gemm!('N', 'N', CT, ws.tmp1, ws.tmp2, ZT, ws.LdagL)  # LdagL = L^T * rho * conj(L)
    BLAS.axpy!(T(scalar), ws.LdagL, out)
    return nothing
end

"""
    _accumulate_adjoint_sandwich_adj_L!(out, L_op, rho, scalar, ws) -> nothing

Accumulate `scalar * conj(L) * rho * L^T` into `out`. This is the HS adjoint of
the negative-frequency sandwich `L^T * rho * conj(L)`.

Note: identical computation to _accumulate_sandwich! (the HS adjoint of
L^T*rho*conj(L) equals conj(L)*rho*L^T, which is the same as the positive-freq forward).
Uses `ws.tmp1`, `ws.tmp2`, `ws.LdagL` as scratch.
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
    @. ws.tmp2 = conj(L_op)                                    # tmp2 = conj(L)
    BLAS.gemm!('N', 'N', CT, ws.tmp2, rho, ZT, ws.tmp1)        # tmp1 = conj(L) * rho
    BLAS.gemm!('N', 'T', CT, ws.tmp1, L_op, ZT, ws.LdagL)     # LdagL = conj(L) * rho * L^T
    BLAS.axpy!(T(scalar), ws.LdagL, out)
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

    # Coherent term: i[B^T, rho] = i*B^T*rho - i*rho*B^T
    # Dense convention: kron(B, I) -> rho*B^T and kron(I, B^T) -> B^T*rho
    B = ws.B_total
    if B !== nothing
        CT = one(T)
        ZT = zero(T)
        BLAS.gemm!('T', 'N', CT, B, rho, ZT, ws.tmp1)   # tmp1 = B^T * rho
        @. ws.rho_out += 1im * ws.tmp1
        BLAS.gemm!('N', 'T', CT, rho, B, ZT, ws.tmp1)    # tmp1 = rho * B^T
        @. ws.rho_out -= 1im * ws.tmp1
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
2. Dissipator sandwich swap: `conj(L) rho L^T` becomes `L^T rho conj(L)` (anticommutator unchanged)

Uses `_accumulate_adjoint_dissipator!` which correctly preserves the `{L'L, rho}`
anticommutator while computing the HS-adjoint sandwich term.

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

    # Adjoint coherent term: -i[B^T, rho] = -i*B^T*rho + i*rho*B^T (sign flip vs forward)
    B = ws.B_total
    if B !== nothing
        CT = one(T)
        ZT = zero(T)
        BLAS.gemm!('T', 'N', CT, B, rho, ZT, ws.tmp1)   # tmp1 = B^T * rho
        @. ws.rho_out -= 1im * ws.tmp1                    # -i instead of +i
        BLAS.gemm!('N', 'T', CT, rho, B, ZT, ws.tmp1)    # tmp1 = rho * B^T
        @. ws.rho_out += 1im * ws.tmp1                    # +i instead of -i
    end

    # Adjoint dissipator: D_L*(rho) = L^T rho conj(L) - 0.5 {L'L, rho}
    # Same L_op as forward, but uses _accumulate_adjoint_dissipator! which
    # computes the HS-adjoint sandwich while keeping anticommutator {L'L, rho}.
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
    _accumulate_dissipator_2op!(out, A, B_dag, rho, scalar, ws) -> nothing

Accumulate `scalar * D_kron(A, B_dag)(rho)` into `out` in-place for the BohrDomain
two-operator dissipator. The formula is derived from the dense code's kron vectorization
convention (see `_vectorize_liouv_diss_and_add!` in qi_tools.jl):

    D_kron(A, B_dag)(rho) = B_dag^T * rho * A^T - 0.5*((B_dag*A)^T * rho + rho * (B_dag*A)^T)

Where `A` = alpha_A_nu1 (jump_1 in the dense code) and `B_dag` = A_nu2_dag (jump_2).
The kron identity `kron(j1, j2^T) vec(rho) = vec(j2^T * rho * j1^T)` with j1=A, j2=B_dag.

Uses `ws.LdagL`, `ws.tmp1`, `ws.tmp2` as scratch. The caller must ensure `A` and `B_dag`
are NOT stored in any of these scratch buffers.
"""
function _accumulate_dissipator_2op!(
    out::Matrix{T},
    A::Matrix{T},
    B_dag::Matrix{T},
    rho::Matrix{T},
    scalar::Real,
    ws,
) where {T<:Complex}
    CT = one(T)
    ZT = zero(T)

    # B_dag * A -> ws.LdagL (anticommutator base product)
    BLAS.gemm!('N', 'N', CT, B_dag, A, ZT, ws.LdagL)

    # Term 2: -0.5 * scalar * (B_dag*A)^T * rho  (anticommutator)
    BLAS.gemm!('T', 'N', CT, ws.LdagL, rho, ZT, ws.tmp1)
    BLAS.axpy!(T(-0.5 * scalar), ws.tmp1, out)

    # Term 3: -0.5 * scalar * rho * (B_dag*A)^T  (anticommutator)
    BLAS.gemm!('N', 'T', CT, rho, ws.LdagL, ZT, ws.tmp1)
    BLAS.axpy!(T(-0.5 * scalar), ws.tmp1, out)

    # Term 1: scalar * B_dag^T * rho * A^T
    BLAS.gemm!('T', 'N', CT, B_dag, rho, ZT, ws.tmp1)      # tmp1 = B_dag^T * rho
    BLAS.gemm!('N', 'T', CT, ws.tmp1, A, ZT, ws.tmp2)       # tmp2 = B_dag^T * rho * A^T
    BLAS.axpy!(T(scalar), ws.tmp2, out)

    return nothing
end

"""
    _accumulate_adjoint_dissipator_2op!(out, A, B_dag, rho, scalar, ws) -> nothing

Accumulate `scalar * D*_kron(A, B_dag)(rho)` (Hilbert-Schmidt adjoint of the two-operator
dissipator) into `out` in-place:

    D*_kron(A, B_dag)(rho) = conj(B_dag) * rho * conj(A) - 0.5*(conj(B_dag*A) * rho + rho * conj(B_dag*A))

Derived from the HS adjoint of the kron-based forward superoperator. The HS adjoint of the
map `rho -> X * rho * Y` is `rho -> X^H * rho * Y^H`, so:
- Sandwich `B_dag^T rho A^T` -> adjoint `conj(B_dag) rho conj(A)`
- Anticomm `(B_dag*A)^T rho` -> adjoint `conj(B_dag*A) rho`

For real-valued operators (BohrDomain's eigenbasis and alpha are real), this simplifies to:
    D*(A, B_dag)(rho) = B_dag * rho * A - 0.5*(B_dag*A * rho + rho * B_dag*A)

Uses `ws.LdagL`, `ws.tmp1`, `ws.tmp2` as scratch.
"""
function _accumulate_adjoint_dissipator_2op!(
    out::Matrix{T},
    A::Matrix{T},
    B_dag::Matrix{T},
    rho::Matrix{T},
    scalar::Real,
    ws,
) where {T<:Complex}
    CT = one(T)
    ZT = zero(T)

    # B_dag * A -> ws.LdagL
    BLAS.gemm!('N', 'N', CT, B_dag, A, ZT, ws.LdagL)

    # conj(B_dag*A) for anticommutator terms
    @. ws.tmp2 = conj(ws.LdagL)                                # tmp2 = conj(B_dag*A)

    # Term 2: -0.5 * scalar * conj(B_dag*A) * rho
    BLAS.gemm!('N', 'N', CT, ws.tmp2, rho, ZT, ws.tmp1)
    BLAS.axpy!(T(-0.5 * scalar), ws.tmp1, out)

    # Term 3: -0.5 * scalar * rho * conj(B_dag*A)
    BLAS.gemm!('N', 'N', CT, rho, ws.tmp2, ZT, ws.tmp1)
    BLAS.axpy!(T(-0.5 * scalar), ws.tmp1, out)

    # Term 1: scalar * conj(B_dag) * rho * conj(A)
    @. ws.tmp2 = conj(B_dag)                                   # tmp2 = conj(B_dag)
    BLAS.gemm!('N', 'N', CT, ws.tmp2, rho, ZT, ws.tmp1)        # tmp1 = conj(B_dag) * rho
    @. ws.LdagL = conj(A)                                      # LdagL = conj(A), free to reuse
    BLAS.gemm!('N', 'N', CT, ws.tmp1, ws.LdagL, ZT, ws.tmp2)  # tmp2 = conj(B_dag) * rho * conj(A)
    BLAS.axpy!(T(scalar), ws.tmp2, out)

    return nothing
end

"""
    _accumulate_sandwich_2op!(out, A, B_dag, rho, scalar, ws) -> nothing

Accumulate `scalar * B_dag^T * rho * A^T` into `out`. This is the sandwich-only
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
    BLAS.gemm!('T', 'N', CT, B_dag, rho, ZT, ws.tmp1)      # tmp1 = B_dag^T * rho
    BLAS.gemm!('N', 'T', CT, ws.tmp1, A, ZT, ws.tmp2)       # tmp2 = B_dag^T * rho * A^T
    BLAS.axpy!(T(scalar), ws.tmp2, out)
    return nothing
end

"""
    _accumulate_adjoint_sandwich_2op!(out, A, B_dag, rho, scalar, ws) -> nothing

Accumulate `scalar * conj(B_dag) * rho * conj(A)` into `out`. This is the HS adjoint
of the two-operator sandwich `B_dag^T * rho * A^T`, used in apply_adjoint_lindbladian!
for BohrDomain.

Uses `ws.tmp1`, `ws.tmp2`, `ws.LdagL` as scratch.
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
    @. ws.tmp2 = conj(B_dag)                                   # tmp2 = conj(B_dag)
    BLAS.gemm!('N', 'N', CT, ws.tmp2, rho, ZT, ws.tmp1)        # tmp1 = conj(B_dag) * rho
    @. ws.LdagL = conj(A)                                      # LdagL = conj(A), reuse as scratch
    BLAS.gemm!('N', 'N', CT, ws.tmp1, ws.LdagL, ZT, ws.tmp2)  # tmp2 = conj(B_dag) * rho * conj(A)
    BLAS.axpy!(T(scalar), ws.tmp2, out)
    return nothing
end

"""
    apply_lindbladian!(ws, rho, config, hamiltonian) -> ws.rho_out

Compute L(rho) for `BohrDomain` configs, storing the result in `ws.rho_out`.

BohrDomain is fundamentally different from Energy/Time/Trotter:
- Iterates over Bohr frequency buckets (`hamiltonian.bohr_dict` keys) instead of energy labels
- Uses a generalized two-operator dissipator D(A, B_dag) where A and B_dag differ
- Entrywise alpha computation: `alpha(bohr_freqs, nu_2) * eigenbasis` per bucket
- No Hermitian half-grid optimization (no energy_labels loop, no transition function)
- Scalar is just `gamma_norm_factor` (no w0/sigma prefactor formula)

A_nu2_dag is built by scattering `conj(eigenbasis[i,j])` to position `[j,i]`
(the dagger operation: swap indices and conjugate). One dense buffer is allocated
per matvec call and reused across bucket iterations (zeroed + scattered each iteration).

# Arguments
- `ws::KrylovWorkspace`: pre-allocated workspace (scratch matrices + precomputed data)
- `rho::Matrix{<:Complex}`: input density matrix (dim x dim)
- `config::AbstractLiouvConfig{BohrDomain}`: Lindbladian configuration
- `hamiltonian::HamHam`: Hamiltonian with bohr_dict and bohr_freqs

# Returns
`ws.rho_out` -- the output matrix (dim x dim), overwritten in-place each call.
"""
function apply_lindbladian!(
    ws::KrylovWorkspace{T},
    rho::Matrix{T},
    config::AbstractLiouvConfig{BohrDomain},
    hamiltonian::HamHam,
) where {T<:Complex}
    (; alpha, gamma_norm_factor) = ws.precomputed_data
    dim = size(rho, 1)
    fill!(ws.rho_out, 0)

    # Coherent term: i[B^T, rho] = i*B^T*rho - i*rho*B^T (identical to other domains)
    # Dense convention: kron(B, I) -> rho*B^T and kron(I, B^T) -> B^T*rho
    B = ws.B_total
    if B !== nothing
        CT = one(T)
        ZT = zero(T)
        BLAS.gemm!('T', 'N', CT, B, rho, ZT, ws.tmp1)   # tmp1 = B^T * rho
        @. ws.rho_out += 1im * ws.tmp1
        BLAS.gemm!('N', 'T', CT, rho, B, ZT, ws.tmp1)    # tmp1 = rho * B^T
        @. ws.rho_out -= 1im * ws.tmp1
    end

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

            _accumulate_dissipator_2op!(ws.rho_out, ws.jump_oft, A_nu2_dag, rho, gamma_norm_factor, ws)
        end
    end

    return ws.rho_out
end

"""
    apply_adjoint_lindbladian!(ws, rho, config, hamiltonian) -> ws.rho_out

Compute L*(rho) (Hilbert-Schmidt adjoint) for `BohrDomain` configs.

Same bucket iteration structure as the forward method but with:
1. Coherent term sign flip: `+i[B, rho]` instead of `-i[B, rho]`
2. Dissipator: uses `_accumulate_adjoint_dissipator_2op!` (dedicated adjoint helper,
   NOT a simple argument swap -- see RESEARCH.md Pitfall 5)

The adjoint two-operator dissipator computes:
    D*(A, B_dag)(rho) = A' * rho * B_dag - 0.5*(A' * B_dag * rho + rho * A' * B_dag)

# Arguments
Same as `apply_lindbladian!` for `BohrDomain`.

# Returns
`ws.rho_out` -- the adjoint Lindbladian action on rho.
"""
function apply_adjoint_lindbladian!(
    ws::KrylovWorkspace{T},
    rho::Matrix{T},
    config::AbstractLiouvConfig{BohrDomain},
    hamiltonian::HamHam,
) where {T<:Complex}
    (; alpha, gamma_norm_factor) = ws.precomputed_data
    dim = size(rho, 1)
    fill!(ws.rho_out, 0)

    # Adjoint coherent term: -i[B^T, rho] = -i*B^T*rho + i*rho*B^T (sign flip vs forward)
    B = ws.B_total
    if B !== nothing
        CT = one(T)
        ZT = zero(T)
        BLAS.gemm!('T', 'N', CT, B, rho, ZT, ws.tmp1)   # tmp1 = B^T * rho
        @. ws.rho_out -= 1im * ws.tmp1                    # -i instead of +i
        BLAS.gemm!('N', 'T', CT, rho, B, ZT, ws.tmp1)    # tmp1 = rho * B^T
        @. ws.rho_out += 1im * ws.tmp1                    # +i instead of -i
    end

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

            # Use dedicated adjoint 2op helper (NOT simple argument swap)
            _accumulate_adjoint_dissipator_2op!(ws.rho_out, ws.jump_oft, A_nu2_dag, rho, gamma_norm_factor, ws)
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

    # Coherent term: i[B^T, rho] = i*B^T*rho - i*rho*B^T (identical to EnergyDomain)
    # Dense convention: kron(B, I) -> rho*B^T and kron(I, B^T) -> B^T*rho
    B = ws.B_total
    if B !== nothing
        CT = one(T)
        ZT = zero(T)
        BLAS.gemm!('T', 'N', CT, B, rho, ZT, ws.tmp1)   # tmp1 = B^T * rho
        @. ws.rho_out += 1im * ws.tmp1
        BLAS.gemm!('N', 'T', CT, rho, B, ZT, ws.tmp1)    # tmp1 = rho * B^T
        @. ws.rho_out -= 1im * ws.tmp1
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
2. Dissipator sandwich swap: `conj(L) rho L^T` becomes `L^T rho conj(L)` (anticommutator unchanged)

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

    # Adjoint coherent term: -i[B^T, rho] = -i*B^T*rho + i*rho*B^T (sign flip vs forward)
    B = ws.B_total
    if B !== nothing
        CT = one(T)
        ZT = zero(T)
        BLAS.gemm!('T', 'N', CT, B, rho, ZT, ws.tmp1)   # tmp1 = B^T * rho
        @. ws.rho_out -= 1im * ws.tmp1                    # -i instead of +i
        BLAS.gemm!('N', 'T', CT, rho, B, ZT, ws.tmp1)    # tmp1 = rho * B^T
        @. ws.rho_out += 1im * ws.tmp1                    # +i instead of -i
    end

    # Adjoint dissipator: D_L*(rho) = L^T rho conj(L) - 0.5 {L'L, rho}
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
