#* Liouvillian (vectorized) jump contributions ------------------------------------------------------------------------------
"""
    Accumulate the Liouvillian contribution of a single jump operator in-place.

    This avoids allocating a full `dim^2 × dim^2` matrix per jump. Call with a
    preallocated `L_target` (dense) and a `LindbladianWorkspace`.

    If `config.with_coherent==true`, pass `coherent_term` already scaled by
    `gamma_norm_factor` to avoid modifying cached matrices.
"""
function jump_contribution!(
    L_target::AbstractMatrix{ComplexF64},
    ::BohrDomain,
    jump::JumpOp,
    hamiltonian::HamHam,
    config::AbstractLiouvConfig,
    precomputed_data,
    ws::LindbladianWorkspace;
    coherent_term::Union{Nothing, AbstractMatrix{ComplexF64}} = nothing,
    )
    dim = size(hamiltonian.data, 1)
    unique_freqs = keys(hamiltonian.bohr_dict)
    (; alpha, gamma_norm_factor) = precomputed_data

    B = coherent_term
    if B !== nothing
        vectorize_liouvillian_coherent!(L_target, B, ws)
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

        vectorize_liouv_diss_and_add!(L_target, alpha_A_nu1, A_nu_2_dag, gamma_norm_factor, ws)
    end
    return L_target
end

function jump_contribution!(
    L_target::AbstractMatrix{ComplexF64},
    ::EnergyDomain,
    jump::JumpOp,
    hamiltonian::HamHam,
    config::AbstractLiouvConfig,
    precomputed_data,
    ws::LindbladianWorkspace;
    coherent_term::Union{Nothing, AbstractMatrix{ComplexF64}} = nothing,
    )

    (; transition, gamma_norm_factor, energy_labels) = precomputed_data

    B = coherent_term
    if B !== nothing
        vectorize_liouvillian_coherent!(L_target, B, ws)
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
            vectorize_liouv_diss_and_add!(L_target, jump_oft, scalar_w, ws)
            if w > 1e-12
                scalar_negative_w = prefactor * transition(-w)
                vectorize_liouv_diss_and_add!(L_target, jump_oft', scalar_negative_w, ws)
            end
        end
    else
        for w in energy_labels
            oft!(jump_oft, jump, w, hamiltonian, config.sigma)
            scalar_w = prefactor * transition(w)
            vectorize_liouv_diss_and_add!(L_target, jump_oft, scalar_w, ws)
        end
    end

    return L_target
end

function jump_contribution!(
    L_target::AbstractMatrix{ComplexF64},
    ::Union{TimeDomain, TrotterDomain},
    jump::JumpOp,
    ham_or_trott::Union{HamHam, TrottTrott},
    config::AbstractLiouvConfig,
    precomputed_data,
    ws::LindbladianWorkspace;
    coherent_term::Union{Nothing, AbstractMatrix{ComplexF64}} = nothing,
    )

    (; transition, gamma_norm_factor, energy_labels, oft_nufft_prefactors, b_minus, b_plus) = precomputed_data

    B = coherent_term
    if B !== nothing
        vectorize_liouvillian_coherent!(L_target, B, ws)
    end

    jump_oft = ws.jump_tmp
    prefactor = config.w0 * config.t0^2 * (config.sigma * sqrt(2 / pi)) / (2 * pi) * gamma_norm_factor

    if jump.hermitian
        for w_raw in energy_labels
            w_raw > 1e-12 && continue
            w = abs(w_raw)
            nufft_prefactor_matrix = prefactor_view(oft_nufft_prefactors, w)
            @. jump_oft = jump.in_eigenbasis * nufft_prefactor_matrix

            scalar_w = prefactor * transition(w)
            vectorize_liouv_diss_and_add!(L_target, jump_oft, scalar_w, ws)
            if w > 1e-12
                scalar_negative_w = prefactor * transition(-w)
                vectorize_liouv_diss_and_add!(L_target, jump_oft', scalar_negative_w, ws)
            end
        end
    else
        for w in energy_labels
            nufft_prefactor_matrix = prefactor_view(oft_nufft_prefactors, w)
            @. jump_oft = jump.in_eigenbasis * nufft_prefactor_matrix
            scalar_w = prefactor * transition(w)
            vectorize_liouv_diss_and_add!(L_target, jump_oft, scalar_w, ws)
        end
    end

    return L_target
end

#* Algorithmic jump contributions -------------------------------------------------------------------------------------------

"""
    apply_coherent_unitary!(evolving_dm, U_B, scratch) -> nothing

Apply coherent unitary evolution: rho -> U_B * rho * U_B'.
No-op if U_B is nothing.
"""
@inline function apply_coherent_unitary!(
    evolving_dm::Matrix{ComplexF64},
    U_B::Union{Nothing,Matrix{ComplexF64}},
    scratch::KrausScratch{ComplexF64},
)
    U_B === nothing && return nothing
    mul!(scratch.tmp1, U_B, evolving_dm)
    mul!(scratch.rho_next, scratch.tmp1, U_B')
    copyto!(evolving_dm, scratch.rho_next)
    return nothing
end

function jump_contribution!(::BohrDomain,
    evolving_dm::Matrix{ComplexF64},
    jump::JumpOp,
    hamiltonian::HamHam,
    config::AbstractThermalizeConfig,
    precomputed_data,
    scratch::KrausScratch{ComplexF64};
    coherent_unitary_cache::Union{Nothing,Matrix{ComplexF64}} = nothing,
    jump_prob::Float64 = 1.0,
    rescale_by_inv_prob::Bool = false
    )

    dim = size(evolving_dm, 1)
    (; alpha, gamma_norm_factor) = precomputed_data

    bohr_keys = hasproperty(precomputed_data, :bohr_keys) ? precomputed_data.bohr_keys : collect(keys(hamiltonian.bohr_dict))
    bohr_is   = hasproperty(precomputed_data, :bohr_is)   ? precomputed_data.bohr_is   : nothing
    bohr_js   = hasproperty(precomputed_data, :bohr_js)   ? precomputed_data.bohr_js   : nothing

    jump_weight_scaling = rescale_by_inv_prob ? (gamma_norm_factor / jump_prob) : gamma_norm_factor
    scaled_delta = config.delta * jump_weight_scaling

    apply_coherent_unitary!(evolving_dm, coherent_unitary_cache, scratch)

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

    # Build K0 = I - (1 - sqrt(1 - delta)) R   (Chen Eq. 3.2)
    delta_factor_for_K0 = 1 - sqrt(1 - config.delta)
    copyto!(scratch.K0, scratch.R)
    scratch.K0 .*= (-delta_factor_for_K0)
    @inbounds for i in 1:dim
        scratch.K0[i,i] += 1
    end

    # Residual TP fix: S := I - K0†K0 - δR = (2α-δ)R - α² R²  (O(δ²))
    mul!(scratch.LdagL, scratch.R, scratch.R)  # reuse as R^2
    s1 = 2 * delta_factor_for_K0 - config.delta
    s2 = delta_factor_for_K0 * delta_factor_for_K0
    @. scratch.tmp2 = s1 * scratch.R - s2 * scratch.LdagL

    eps_shift = 10 * eps(Float64)
    @inbounds for i in 1:dim
        scratch.tmp2[i,i] += eps_shift
    end

    cholesky_S = cholesky!(Hermitian(scratch.tmp2), check=false)
    U_residual = cholesky_S.U

    # ρ_next = K0 ρ K0† + ρ_jump + Ures ρ Ures†
    mul!(scratch.tmp1, scratch.K0, evolving_dm)
    mul!(scratch.rho_next, scratch.tmp1, scratch.K0')
    scratch.rho_next .+= scratch.rho_jump

    mul!(scratch.tmp1, U_residual, evolving_dm)
    mul!(scratch.rho_next, scratch.tmp1, U_residual', 1.0, 1.0)

    hermitianize!(scratch.rho_next)
    copyto!(evolving_dm, scratch.rho_next)
    return evolving_dm
end

function jump_contribution!(::EnergyDomain,
    evolving_dm::Matrix{ComplexF64},
    jump::JumpOp,
    hamiltonian::HamHam,
    config::AbstractThermalizeConfig,
    precomputed_data,
    scratch::KrausScratch{ComplexF64};
    coherent_unitary_cache::Union{Nothing,Matrix{ComplexF64}} = nothing,
    jump_prob::Float64 = 1.0,
    rescale_by_inv_prob::Bool = false
    )

    dim = size(evolving_dm, 1)
    (; transition, gamma_norm_factor, energy_labels) = precomputed_data

    jump_weight_scaling = rescale_by_inv_prob ? (gamma_norm_factor / jump_prob) : gamma_norm_factor

    apply_coherent_unitary!(evolving_dm, coherent_unitary_cache, scratch)

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

    # Build K0 = I - (1 - sqrt(1 - delta)) R   (Chen Eq. 3.2)
    delta_factor_for_K0 = 1 - sqrt(1 - config.delta)
    copyto!(scratch.K0, scratch.R)
    scratch.K0 .*= (-delta_factor_for_K0)

    @inbounds for i in 1:dim
        scratch.K0[i,i] += 1
    end

    # Residual TP fix: S := I - K0†K0 - δR = (2α-δ)R - α² R²  (O(δ²))
    mul!(scratch.LdagL, scratch.R, scratch.R)  # reuse as R^2
    s1 = 2 * delta_factor_for_K0 - config.delta
    s2 = delta_factor_for_K0 * delta_factor_for_K0
    @. scratch.tmp2 = s1 * scratch.R - s2 * scratch.LdagL

    eps_shift = 10 * eps(Float64)
    @inbounds for i in 1:dim
        scratch.tmp2[i,i] += eps_shift
    end

    cholesky_S = cholesky!(Hermitian(scratch.tmp2), check=false)
    U_residual = cholesky_S.U

    # ρ_next = K0 ρ K0† + ρ_jump + Ures ρ Ures†
    mul!(scratch.tmp1, scratch.K0, evolving_dm)
    mul!(scratch.rho_next, scratch.tmp1, scratch.K0')
    scratch.rho_next .+= scratch.rho_jump

    mul!(scratch.tmp1, U_residual, evolving_dm)
    mul!(scratch.rho_next, scratch.tmp1, U_residual', 1.0, 1.0)

    hermitianize!(scratch.rho_next)
    copyto!(evolving_dm, scratch.rho_next)

    return evolving_dm
end

function jump_contribution!(::Union{TimeDomain, TrotterDomain},
    evolving_dm::Matrix{ComplexF64},
    jump::JumpOp,
    ham_or_trott,              # HamHam or TrottTrott depending on domain
    config::AbstractThermalizeConfig,
    precomputed_data,
    scratch::KrausScratch{ComplexF64};
    coherent_unitary_cache::Union{Nothing,Matrix{ComplexF64}} = nothing,
    jump_prob::Float64 = 1.0,
    rescale_by_inv_prob::Bool = false
    )

    dim = size(evolving_dm, 1)
    (; transition, gamma_norm_factor, energy_labels, oft_nufft_prefactors, b_minus, b_plus) = precomputed_data

    jump_weight_scaling = rescale_by_inv_prob ? (gamma_norm_factor / jump_prob) : gamma_norm_factor

    apply_coherent_unitary!(evolving_dm, coherent_unitary_cache, scratch)

    base_prefactor = config.w0 * config.t0^2 * (config.sigma * sqrt(2 / pi)) / (2 * pi) * jump_weight_scaling

    fill!(scratch.R, 0)
    fill!(scratch.rho_jump, 0)

    energies = jump.hermitian ? abs.(filter(w -> w < 1e-12, energy_labels)) : energy_labels
    for w in energies
        nufft_prefactor_matrix = prefactor_view(oft_nufft_prefactors, w)

        # Aω := A .* nufft_prefactor_matrix(ω)
        @. scratch.jump_oft = jump.in_eigenbasis * nufft_prefactor_matrix

        # rate^2(ω)
        rate2_pos = base_prefactor * transition(w)

        # R += rate^2 * (Aω† Aω)
        mul!(scratch.LdagL, scratch.jump_oft', scratch.jump_oft)
        @. scratch.R += rate2_pos * scratch.LdagL

        # rho_jump += delta * rate^2 * (Aω ρ Aω†)
        mul!(scratch.tmp1, evolving_dm, scratch.jump_oft')  # ρ Aω†
        mul!(scratch.rho_jump, scratch.jump_oft, scratch.tmp1, config.delta*rate2_pos, 1.0)

        if jump.hermitian && w > 1e-12
            rate2_neg = base_prefactor * transition(-w)

            # For the negative-frequency partner, the Lindblad operator is (Aω)†.
            # Then L†L = Aω Aω†, and jump term is Aω† ρ Aω.
            mul!(scratch.LdagL, scratch.jump_oft, scratch.jump_oft')
            @. scratch.R += rate2_neg * scratch.LdagL

            mul!(scratch.tmp1, evolving_dm, scratch.jump_oft)              # ρ Aω
            mul!(scratch.rho_jump, scratch.jump_oft', scratch.tmp1, config.delta*rate2_neg, 1.0)
        end
    end

    # Hermitianize R
    hermitianize!(scratch.R)

    # Build K0 = I - (1 - sqrt(1 - delta)) R
    delta_factor_for_K0 = 1 - sqrt(1 - config.delta)  # Chen: Eq. 3.2
    copyto!(scratch.K0, scratch.R)
    # scratch.K0 .*= (0.5 * config.delta)
    scratch.K0 .*= (- delta_factor_for_K0)

    # add identity on the diagonal
    @inbounds for i in 1:dim
        scratch.K0[i,i] += 1
    end

    # (Necessary) residual delta^2 trash to have TP channel
    # S := I - K0†K0 - δ R  = (2α-δ)R - α² R²
    #     Choose Kres s.t. Kres†Kres = S. Use Cholesky: S = U' U => Kres = U.

    # Reuse LdagL as R^2
    mul!(scratch.LdagL, scratch.R, scratch.R)
    s1 = 2 * delta_factor_for_K0 - config.delta
    s2 = delta_factor_for_K0 * delta_factor_for_K0

    # tmp2 = S
    @. scratch.tmp2 = s1 * scratch.R - s2 * scratch.LdagL

    # Guard against tiny negative eigenvalues from roundoff (S is O(δ²)):
    # add a microscopic diagonal shift before Cholesky.
    eps_shift = 10 * eps(Float64)
    @inbounds for i in 1:dim
        scratch.tmp2[i,i] += eps_shift
    end

    cholesky_S = cholesky!(Hermitian(scratch.tmp2), check=false)  # S ≈ U' U
    U_residual = cholesky_S.U  

    # Combine into CPTP channel, delta step of weak measurement scheme
    # ρ_next = K0 ρ K0† + ρ_jump + Ures ρ Ures†
    mul!(scratch.tmp1, scratch.K0, evolving_dm)
    mul!(scratch.rho_next, scratch.tmp1, scratch.K0')
    scratch.rho_next .+= scratch.rho_jump

    # residual CP term
    mul!(scratch.tmp1, U_residual, evolving_dm)
    mul!(scratch.rho_next, scratch.tmp1, U_residual', 1.0, 1.0)

    # Keep it a density matrix numerically
    hermitianize!(scratch.rho_next)
    copyto!(evolving_dm, scratch.rho_next)

    return evolving_dm
end

#TODO: test it; set BLAS threads to 1, let julia threads be more.