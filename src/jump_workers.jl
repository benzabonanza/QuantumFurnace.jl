#* Liouvillian (vectorized) jump contributions ------------------------------------------------------------------------------
"""
    Accumulate the Liouvillian contribution of a single jump operator in-place.

    This avoids allocating a full `dim^2 × dim^2` matrix per jump. Call with a
    preallocated `L_target` (dense) and a `LindbladianWorkspace`.

    If `with_coherent(config.construction)==true`, pass `coherent_term` already scaled by
    `gamma_norm_factor` to avoid modifying cached matrices.
"""
function _jump_contribution!(
    L_target::AbstractMatrix{<:Complex},
    jump::JumpOp,
    hamiltonian::HamHam,
    config::Config{Lindbladian, BohrDomain},
    precomputed_data,
    ws::LindbladianWorkspace;
    coherent_term::Union{Nothing, AbstractMatrix{<:Complex}} = nothing,
    )
    dim = size(hamiltonian.data, 1)
    unique_freqs = keys(hamiltonian.bohr_dict)
    (; alpha, gamma_norm_factor) = precomputed_data

    B = coherent_term
    if B !== nothing
        _vectorize_liouvillian_coherent!(L_target, B, ws)
    end

    alpha_A_nu1 = ws.jump_tmp
    for nu_2 in unique_freqs
        @. alpha_A_nu1 = alpha(hamiltonian.bohr_freqs, nu_2) * jump.in_eigenbasis

        indices = hamiltonian.bohr_dict[nu_2]
        A_nu_2_vals = view(jump.in_eigenbasis, indices)

        # swapped for dagger
        rows_dag = getindex.(indices, 2)
        cols_dag = getindex.(indices, 1)

        A_nu_2_dag = sparse(rows_dag, cols_dag, conj.(A_nu_2_vals), dim, dim)

        _vectorize_liouv_diss_and_add!(L_target, alpha_A_nu1, A_nu_2_dag, gamma_norm_factor, ws)
    end
    return L_target
end

function _jump_contribution!(
    L_target::AbstractMatrix{<:Complex},
    jump::JumpOp,
    hamiltonian::HamHam,
    config::Config{Lindbladian, EnergyDomain},
    precomputed_data,
    ws::LindbladianWorkspace;
    coherent_term::Union{Nothing, AbstractMatrix{<:Complex}} = nothing,
    )

    (; transition, gamma_norm_factor, energy_labels) = precomputed_data

    B = coherent_term
    if B !== nothing
        _vectorize_liouvillian_coherent!(L_target, B, ws)
    end

    jump_oft = ws.jump_tmp
    prefactor = (config.w0 / (config.sigma * sqrt(2 * pi))) * gamma_norm_factor

    if jump.hermitian
        for w_raw in energy_labels
            # iterate only half-grid (w<=0) and mirror manually
            w_raw > 1e-12 && continue
            w = abs(w_raw)
            oft!(jump_oft, jump, w, hamiltonian, config.sigma)
            scalar_w = prefactor * transition(w)
            _vectorize_liouv_diss_and_add!(L_target, jump_oft, scalar_w, ws)
            if w > 1e-12
                scalar_negative_w = prefactor * transition(-w)
                _vectorize_liouv_diss_and_add!(L_target, jump_oft', scalar_negative_w, ws)
            end
        end
    else
        for w in energy_labels
            oft!(jump_oft, jump, w, hamiltonian, config.sigma)
            scalar_w = prefactor * transition(w)
            _vectorize_liouv_diss_and_add!(L_target, jump_oft, scalar_w, ws)
        end
    end

    return L_target
end

function _jump_contribution!(
    L_target::AbstractMatrix{<:Complex},
    jump::JumpOp,
    ham_or_trott::Union{HamHam, TrottTrott},
    config::Config{Lindbladian, D},
    precomputed_data,
    ws::LindbladianWorkspace;
    coherent_term::Union{Nothing, AbstractMatrix{<:Complex}} = nothing,
    ) where {D<:Union{TimeDomain, TrotterDomain}}

    (; transition, gamma_norm_factor, energy_labels, oft_nufft_prefactors, b_minus, b_plus) = precomputed_data

    B = coherent_term
    if B !== nothing
        _vectorize_liouvillian_coherent!(L_target, B, ws)
    end

    jump_oft = ws.jump_tmp
    prefactor = config.w0 * config.t0^2 * (config.sigma * sqrt(2 / pi)) / (2 * pi) * gamma_norm_factor

    if jump.hermitian
        for w_raw in energy_labels
            w_raw > 1e-12 && continue
            w = abs(w_raw)
            nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)
            @. jump_oft = jump.in_eigenbasis * nufft_prefactor_matrix

            scalar_w = prefactor * transition(w)
            _vectorize_liouv_diss_and_add!(L_target, jump_oft, scalar_w, ws)
            if w > 1e-12
                scalar_negative_w = prefactor * transition(-w)
                _vectorize_liouv_diss_and_add!(L_target, jump_oft', scalar_negative_w, ws)
            end
        end
    else
        for w in energy_labels
            nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)
            @. jump_oft = jump.in_eigenbasis * nufft_prefactor_matrix
            scalar_w = prefactor * transition(w)
            _vectorize_liouv_diss_and_add!(L_target, jump_oft, scalar_w, ws)
        end
    end

    return L_target
end

#* Algorithmic jump contributions -------------------------------------------------------------------------------------------

"""
    _apply_coherent_unitary!(evolving_dm, U_B, scratch) -> nothing

Apply coherent unitary evolution: rho -> U_B * rho * U_B'.
No-op if U_B is nothing.
"""
@inline function _apply_coherent_unitary!(
    evolving_dm::Matrix{<:Complex},
    U_B::Union{Nothing,Matrix{<:Complex}},
    scratch::KrausScratch{<:Complex},
)
    U_B === nothing && return nothing
    mul!(scratch.tmp1, U_B, evolving_dm)
    mul!(scratch.rho_next, scratch.tmp1, U_B')
    copyto!(evolving_dm, scratch.rho_next)
    return nothing
end

"""
    _finalize_kraus_step!(evolving_dm, delta, scratch) -> evolving_dm

Apply the CPTP weak-measurement channel after R and rho_jump have been accumulated.

Implements Chen Eq. 3.2:
  K0 = I - alpha * R,  alpha = 1 - sqrt(1 - delta)
  S  = (2*alpha - delta)*R - alpha^2 * R^2  (residual, O(delta^2))
  U_residual = sqrt_psd(S)  (eigendecomposition with clamped eigenvalues)
  rho_next = K0 * rho * K0' + rho_jump + U_res * rho * U_res'

Expects scratch.R (Hermitianized) and scratch.rho_jump to be pre-filled by the
domain-specific dissipative accumulation loop.
"""
function _finalize_kraus_step!(
    evolving_dm::Matrix{<:Complex},
    delta::Real,
    scratch::KrausScratch{<:Complex},
)
    dim = size(evolving_dm, 1)

    # Build K0 = I - alpha * R   (Chen Eq. 3.2)
    delta_factor_for_K0 = 1 - sqrt(1 - delta)
    copyto!(scratch.K0, scratch.R)
    scratch.K0 .*= (-delta_factor_for_K0)
    @inbounds for i in 1:dim
        scratch.K0[i,i] += 1
    end

    # Residual TP fix: S := I - K0'K0 - delta*R = (2*alpha - delta)*R - alpha^2 * R^2
    mul!(scratch.LdagL, scratch.R, scratch.R)  # reuse as R^2
    s1 = 2 * delta_factor_for_K0 - delta
    s2 = delta_factor_for_K0 * delta_factor_for_K0
    @. scratch.tmp2 = s1 * scratch.R - s2 * scratch.LdagL

    # PSD guard: clamp negative eigenvalues to zero (more robust than Cholesky + eps shift)
    hermitianize!(scratch.tmp2)
    S_herm = Hermitian(scratch.tmp2)
    eig = eigen(S_herm)
    eig.values .= max.(eig.values, 0.0)
    CT = eltype(evolving_dm)
    U_residual = Matrix{CT}(Diagonal(sqrt.(eig.values)) * eig.vectors')

    # rho_next = K0 * rho * K0' + rho_jump + U_res * rho * U_res'
    mul!(scratch.tmp1, scratch.K0, evolving_dm)
    mul!(scratch.rho_next, scratch.tmp1, scratch.K0')
    scratch.rho_next .+= scratch.rho_jump

    mul!(scratch.tmp1, U_residual, evolving_dm)
    mul!(scratch.rho_next, scratch.tmp1, U_residual', 1.0, 1.0)

    # Keep it a density matrix numerically
    hermitianize!(scratch.rho_next)
    copyto!(evolving_dm, scratch.rho_next)
    return evolving_dm
end

function _jump_contribution!(
    evolving_dm::Matrix{<:Complex},
    jump::JumpOp,
    hamiltonian::HamHam,
    config::Config{Thermalize, BohrDomain},
    precomputed_data,
    scratch::KrausScratch{<:Complex};
    coherent_unitary_cache::Union{Nothing,Matrix{<:Complex}} = nothing,
    jump_prob::Real = 1.0,
    rescale_by_inv_prob::Bool = false
    )

    dim = size(evolving_dm, 1)
    (; alpha, gamma_norm_factor) = precomputed_data

    bohr_keys = hasproperty(precomputed_data, :bohr_keys) ? precomputed_data.bohr_keys : collect(keys(hamiltonian.bohr_dict))
    bohr_is   = hasproperty(precomputed_data, :bohr_is)   ? precomputed_data.bohr_is   : nothing
    bohr_js   = hasproperty(precomputed_data, :bohr_js)   ? precomputed_data.bohr_js   : nothing

    jump_weight_scaling = rescale_by_inv_prob ? (gamma_norm_factor / jump_prob) : gamma_norm_factor
    scaled_delta = config.delta * jump_weight_scaling

    _apply_coherent_unitary!(evolving_dm, coherent_unitary_cache, scratch)

    fill!(scratch.R, 0)
    fill!(scratch.rho_jump, 0)

    # For each fixed "right" Bohr label ν2 build the composite
    #   B_{ν2} = \sum_{ν1} α_{ν1,ν2} A_{ν1}
    # and accumulate
    #   rho_jump += δ * B_{ν2} ρ A_{ν2}†
    #   R        +=     A_{ν2}† B_{ν2}
    @inbounds for (k, nu_2) in pairs(bohr_keys)
        # B_{ν2}
        @. scratch.jump_oft = alpha(hamiltonian.bohr_freqs, nu_2) * jump.in_eigenbasis

        # tmp1 := ρ A_{ν2}† without explicitly building A_{ν2}†.
        # If (i,j) is in the ν2 bucket, then (A_{ν2}†)_{j,i} = conj(A_{i,j}).
        fill!(scratch.tmp1, 0)
        if bohr_is !== nothing
            is = bohr_is[k]
            js = bohr_js[k]
            @inbounds for t in eachindex(is)
                i = is[t]
                j = js[t]
                v = conj(jump.in_eigenbasis[i, j])
                @inbounds for p in 1:dim
                    scratch.tmp1[p, i] += evolving_dm[p, j] * v
                end
            end
        else
            indices = hamiltonian.bohr_dict[nu_2]
            @inbounds for idx in indices
                i = idx[1]
                j = idx[2]
                v = conj(jump.in_eigenbasis[i, j])
                @inbounds for p in 1:dim
                    scratch.tmp1[p, i] += evolving_dm[p, j] * v
                end
            end
        end

        # rho_jump += δ * B_{ν2} * (ρ A_{ν2}†)
        mul!(scratch.rho_jump, scratch.jump_oft, scratch.tmp1, scaled_delta, 1.0)

        # R += A_{ν2}† * B_{ν2}  (no δ factor)
        if bohr_is !== nothing
            is = bohr_is[k]
            js = bohr_js[k]
            @inbounds for t in eachindex(is)
                i = is[t]
                j = js[t]
                v = conj(jump.in_eigenbasis[i, j]) * jump_weight_scaling
                @inbounds for q in 1:dim
                    scratch.R[j, q] += v * scratch.jump_oft[i, q]
                end
            end
        else
            indices = hamiltonian.bohr_dict[nu_2]
            @inbounds for idx in indices
                i = idx[1]
                j = idx[2]
                v = conj(jump.in_eigenbasis[i, j]) * jump_weight_scaling
                @inbounds for q in 1:dim
                    scratch.R[j, q] += v * scratch.jump_oft[i, q]
                end
            end
        end
    end

    # Hermitianize R (numerical)
    hermitianize!(scratch.R)

    # Apply R, K0, U_residual
    _finalize_kraus_step!(evolving_dm, config.delta, scratch)
    return evolving_dm
end

function _jump_contribution!(
    evolving_dm::Matrix{<:Complex},
    jump::JumpOp,
    hamiltonian::HamHam,
    config::Config{Thermalize, EnergyDomain},
    precomputed_data,
    scratch::KrausScratch{<:Complex};
    coherent_unitary_cache::Union{Nothing,Matrix{<:Complex}} = nothing,
    jump_prob::Real = 1.0,
    rescale_by_inv_prob::Bool = false
    )

    dim = size(evolving_dm, 1)
    (; transition, gamma_norm_factor, energy_labels) = precomputed_data

    jump_weight_scaling = rescale_by_inv_prob ? (gamma_norm_factor / jump_prob) : gamma_norm_factor

    _apply_coherent_unitary!(evolving_dm, coherent_unitary_cache, scratch)

    # --- Dissipative part ---
    # Matches the Euler prefactor in jump_contribution(::EnergyDomain, ...):
    # prefactor = (delta * gamma_norm_factor) * w0 / (sigma*sqrt(2π))
    base_prefactor = config.w0 / (config.sigma * sqrt(2 * pi)) * jump_weight_scaling

    fill!(scratch.R, 0)
    fill!(scratch.rho_jump, 0)

    if jump.hermitian
        @inbounds for w_raw in energy_labels
            w_raw > 1e-12 && continue
            w = abs(w_raw)

            # Aω
            oft!(scratch.jump_oft, jump, w, hamiltonian, config.sigma)

            rate2_pos = base_prefactor * transition(w)

            # R += rate^2 * (Aω† Aω)
            mul!(scratch.LdagL, scratch.jump_oft', scratch.jump_oft)
            @. scratch.R += rate2_pos * scratch.LdagL

            # rho_jump += delta * rate^2 * (Aω ρ Aω†)
            mul!(scratch.tmp1, evolving_dm, scratch.jump_oft')  # ρ Aω†
            mul!(scratch.rho_jump, scratch.jump_oft, scratch.tmp1, config.delta * rate2_pos, 1.0)

            if w > 1e-12
                rate2_neg = base_prefactor * transition(-w)

                # Negative-frequency partner uses (Aω)† as Lindblad operator.
                # Then L†L = Aω Aω† and jump term is Aω† ρ Aω.
                mul!(scratch.LdagL, scratch.jump_oft, scratch.jump_oft')
                @. scratch.R += rate2_neg * scratch.LdagL

                mul!(scratch.tmp1, evolving_dm, scratch.jump_oft)  # ρ Aω
                mul!(scratch.rho_jump, scratch.jump_oft', scratch.tmp1, config.delta * rate2_neg, 1.0)
            end
        end
    else
        @inbounds for w in energy_labels
            oft!(scratch.jump_oft, jump, w, hamiltonian, config.sigma)

            rate2 = base_prefactor * transition(w)

            mul!(scratch.LdagL, scratch.jump_oft', scratch.jump_oft)
            @. scratch.R += rate2 * scratch.LdagL

            mul!(scratch.tmp1, evolving_dm, scratch.jump_oft')
            mul!(scratch.rho_jump, scratch.jump_oft, scratch.tmp1, config.delta * rate2, 1.0)
        end
    end

    # Hermitianize R (numerical)
    hermitianize!(scratch.R)

    # Apply R, K0, U_residual
    _finalize_kraus_step!(evolving_dm, config.delta, scratch)
    return evolving_dm
end

function _jump_contribution!(
    evolving_dm::Matrix{<:Complex},
    jump::JumpOp,
    ham_or_trott,              # HamHam or TrottTrott depending on domain
    config::Config{Thermalize, D},
    precomputed_data,
    scratch::KrausScratch{<:Complex};
    coherent_unitary_cache::Union{Nothing,Matrix{<:Complex}} = nothing,
    jump_prob::Real = 1.0,
    rescale_by_inv_prob::Bool = false
    ) where {D<:Union{TimeDomain, TrotterDomain}}

    dim = size(evolving_dm, 1)
    (; transition, gamma_norm_factor, energy_labels, oft_nufft_prefactors, b_minus, b_plus) = precomputed_data

    jump_weight_scaling = rescale_by_inv_prob ? (gamma_norm_factor / jump_prob) : gamma_norm_factor

    _apply_coherent_unitary!(evolving_dm, coherent_unitary_cache, scratch)

    base_prefactor = config.w0 * config.t0^2 * (config.sigma * sqrt(2 / pi)) / (2 * pi) * jump_weight_scaling

    fill!(scratch.R, 0)
    fill!(scratch.rho_jump, 0)

    if jump.hermitian
        @inbounds for w_raw in energy_labels
            w_raw > 1e-12 && continue
            w = abs(w_raw)

            nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)
            @. scratch.jump_oft = jump.in_eigenbasis * nufft_prefactor_matrix

            rate2_pos = base_prefactor * transition(w)

            mul!(scratch.LdagL, scratch.jump_oft', scratch.jump_oft)
            @. scratch.R += rate2_pos * scratch.LdagL

            mul!(scratch.tmp1, evolving_dm, scratch.jump_oft')
            mul!(scratch.rho_jump, scratch.jump_oft, scratch.tmp1, config.delta*rate2_pos, 1.0)

            if w > 1e-12
                rate2_neg = base_prefactor * transition(-w)

                mul!(scratch.LdagL, scratch.jump_oft, scratch.jump_oft')
                @. scratch.R += rate2_neg * scratch.LdagL

                mul!(scratch.tmp1, evolving_dm, scratch.jump_oft)
                mul!(scratch.rho_jump, scratch.jump_oft', scratch.tmp1, config.delta*rate2_neg, 1.0)
            end
        end
    else
        for w in energy_labels
            nufft_prefactor_matrix = _prefactor_view(oft_nufft_prefactors, w)
            @. scratch.jump_oft = jump.in_eigenbasis * nufft_prefactor_matrix

            rate2_pos = base_prefactor * transition(w)

            mul!(scratch.LdagL, scratch.jump_oft', scratch.jump_oft)
            @. scratch.R += rate2_pos * scratch.LdagL

            mul!(scratch.tmp1, evolving_dm, scratch.jump_oft')
            mul!(scratch.rho_jump, scratch.jump_oft, scratch.tmp1, config.delta*rate2_pos, 1.0)
        end
    end

    # Hermitianize R
    hermitianize!(scratch.R)

    # Apply R, K0, U_residual
    _finalize_kraus_step!(evolving_dm, config.delta, scratch)
    return evolving_dm
end