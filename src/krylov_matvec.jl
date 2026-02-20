"""
    _accumulate_dissipator!(out, L_op, rho, scalar, ws) -> nothing

Accumulate `scalar * D_L(rho)` into `out` in-place:

    D_L(rho) = L * rho * L' - 0.5 * (L'L * rho + rho * L'L)

Uses `ws.tmp1`, `ws.tmp2`, `ws.LdagL` as scratch. All operations are
allocation-free for `Matrix` and `Adjoint{Matrix}` arguments.
"""
function _accumulate_dissipator!(
    out::Matrix{<:Complex},
    L_op::AbstractMatrix{<:Complex},
    rho::AbstractMatrix{<:Complex},
    scalar::Real,
    ws,  # KrylovWorkspace -- uses tmp1, tmp2, LdagL
)
    # L'L -> ws.LdagL
    mul!(ws.LdagL, L_op', L_op)

    # Term 1: scalar * L * rho * L'
    mul!(ws.tmp1, L_op, rho)         # tmp1 = L * rho
    mul!(ws.tmp2, ws.tmp1, L_op')    # tmp2 = L * rho * L'
    @. out += scalar * ws.tmp2

    # Term 2: -0.5 * scalar * L'L * rho
    mul!(ws.tmp1, ws.LdagL, rho)     # tmp1 = L'L * rho
    @. out -= (0.5 * scalar) * ws.tmp1

    # Term 3: -0.5 * scalar * rho * L'L
    mul!(ws.tmp1, rho, ws.LdagL)     # tmp1 = rho * L'L
    @. out -= (0.5 * scalar) * ws.tmp1

    return nothing
end

"""
    apply_lindbladian!(ws, rho, config, hamiltonian) -> ws.rho_out

Compute L(rho) for `EnergyDomain` configs, storing the result in `ws.rho_out`.

Mirrors `_jump_contribution!` for `AbstractLiouvConfig{EnergyDomain}` in
`jump_workers.jl` line-by-line, replacing vectorized superoperator accumulation
with direct matrix-on-matrix operations via `_accumulate_dissipator!`.

# Arguments
- `ws::KrylovWorkspace`: pre-allocated workspace (scratch matrices + precomputed data)
- `rho::AbstractMatrix{<:Complex}`: input density matrix (dim x dim)
- `config::AbstractLiouvConfig{EnergyDomain}`: Lindbladian configuration
- `hamiltonian::HamHam`: Hamiltonian with eigenbasis data

# Returns
`ws.rho_out` -- the output matrix (dim x dim), overwritten in-place each call.
"""
function apply_lindbladian!(
    ws::KrylovWorkspace,
    rho::AbstractMatrix{<:Complex},
    config::AbstractLiouvConfig{EnergyDomain},
    hamiltonian::HamHam,
)
    (; transition, gamma_norm_factor, energy_labels) = ws.precomputed_data

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

    for jump in ws.jumps
        if jump.hermitian
            for w_raw in energy_labels
                # Iterate only half-grid (w <= 0) and mirror manually
                w_raw > 1e-12 && continue
                w = abs(w_raw)

                oft!(ws.jump_oft, jump, w, hamiltonian, config.sigma)

                scalar_w = prefactor * transition(w)
                _accumulate_dissipator!(ws.rho_out, ws.jump_oft, rho, scalar_w, ws)

                if w > 1e-12
                    # Negative-frequency partner: L = A(omega)'
                    scalar_neg = prefactor * transition(-w)
                    _accumulate_dissipator!(ws.rho_out, ws.jump_oft', rho, scalar_neg, ws)
                end
            end
        else
            for w in energy_labels
                oft!(ws.jump_oft, jump, w, hamiltonian, config.sigma)
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
2. Dissipator L <-> L' swap: where forward passes `L_op`, adjoint passes `L_op'`

Shares the same `KrylovWorkspace` as the forward function.

# Arguments
Same as `apply_lindbladian!`.

# Returns
`ws.rho_out` -- the adjoint Lindbladian action on rho.
"""
function apply_adjoint_lindbladian!(
    ws::KrylovWorkspace,
    rho::AbstractMatrix{<:Complex},
    config::AbstractLiouvConfig{EnergyDomain},
    hamiltonian::HamHam,
)
    (; transition, gamma_norm_factor, energy_labels) = ws.precomputed_data

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

    # Adjoint dissipator: swap L <-> L'
    # D*(rho) = L' * rho * L - 0.5*(L*L' * rho + rho * L*L')
    # Achieved by passing L_op' where forward passes L_op (and vice versa)
    prefactor = (config.w0 / (config.sigma * sqrt(2 * pi))) * gamma_norm_factor

    for jump in ws.jumps
        if jump.hermitian
            for w_raw in energy_labels
                w_raw > 1e-12 && continue
                w = abs(w_raw)

                oft!(ws.jump_oft, jump, w, hamiltonian, config.sigma)

                scalar_w = prefactor * transition(w)
                # Forward uses L=A(w), adjoint uses L=A(w)' (swap)
                _accumulate_dissipator!(ws.rho_out, ws.jump_oft', rho, scalar_w, ws)

                if w > 1e-12
                    scalar_neg = prefactor * transition(-w)
                    # Forward uses L=A(w)', adjoint uses L=(A(w)')' = A(w) (swap)
                    _accumulate_dissipator!(ws.rho_out, ws.jump_oft, rho, scalar_neg, ws)
                end
            end
        else
            for w in energy_labels
                oft!(ws.jump_oft, jump, w, hamiltonian, config.sigma)
                scalar_w = prefactor * transition(w)
                # Forward uses L=A(w), adjoint uses L=A(w)' (swap)
                _accumulate_dissipator!(ws.rho_out, ws.jump_oft', rho, scalar_w, ws)
            end
        end
    end

    return ws.rho_out
end
