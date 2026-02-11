struct KrausScratch{T}
    jump_oft::Matrix{T}
    LdagL::Matrix{T}
    R::Matrix{T}
    rho_jump::Matrix{T}
    K0::Matrix{T}
    tmp1::Matrix{T}
    tmp2::Matrix{T}
    rho_next::Matrix{T}
end

function KrausScratch(T::Type{ComplexF64}, dim::Int)
    Zm() = zeros(T, dim, dim)
    return KrausScratch(Zm(), Zm(), Zm(), Zm(), Zm(), Zm(), Zm(), Zm())
end

function apply_kraus_step!(::EnergyDomain,
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

    # Evolve via exp(-i delta B^a)
    U_B = coherent_unitary_cache
    if U_B !== nothing
        mul!(scratch.tmp1, U_B, evolving_dm)
        mul!(scratch.rho_next, scratch.tmp1, U_B')
        copyto!(evolving_dm, scratch.rho_next)
    end

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
    scratch.R .= 0.5 .* (scratch.R .+ scratch.R')

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

    evolving_dm .= 0.5 .* (scratch.rho_next .+ scratch.rho_next')

    return evolving_dm
end

function apply_kraus_step!(::Union{TimeDomain, TrotterDomain},
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

    # Evolve via exp(-i delta B^a)
    U_B = coherent_unitary_cache
    if U_B !== nothing
        mul!(scratch.tmp1, U_B, evolving_dm)
        mul!(scratch.rho_next, scratch.tmp1, U_B')
        copyto!(evolving_dm, scratch.rho_next)
    end

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
    scratch.R .= 0.5 .* (scratch.R .+ scratch.R')

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
    evolving_dm .= 0.5 .* (scratch.rho_next .+ scratch.rho_next')

    return evolving_dm
end