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

function apply_kraus_step!(::Union{TimeDomain, TrotterDomain},
    evolving_dm::Matrix{ComplexF64},
    jump::JumpOp,
    ham_or_trott,              # HamHam or TrottTrott depending on domain
    config::ThermalizeConfig,
    precomputed_data,
    scratch::KrausScratch{ComplexF64};
    coherent_unitary_cache::Union{Nothing,Matrix{ComplexF64}} = nothing,
    # If you randomize jump choice with prob p, and you want the expected generator to match
    # the *sum over all jumps* version, set rescale_by_inv_prob=true and pass p.
    jump_prob::Float64 = 1.0,
    rescale_by_inv_prob::Bool = false
)

    dim = size(evolving_dm, 1)
    (; transition, gamma_norm_factor, energy_labels, oft_nufft_prefactors, b_minus, b_plus) = precomputed_data

    jump_weight_scaling = rescale_by_inv_prob ? (gamma_norm_factor / jump_prob) : gamma_norm_factor

    # Coherent term H := jump_weight_scaling * B(jump)  (so dt*iH matches Euler's -i*dt*gamma_norm_factor*[B,ρ])
    if config.with_coherent
        U_B = coherent_unitary_cache
        if U_B === nothing  # Fallback: compute now (allocates unless your B_* has an in-place version)
            B = (config.domain isa TrotterDomain) ?  #TODO: Check again how Hermitian B's are, if not so much, Hermitianize it here or there
                B_trotter(jump, ham_or_trott, b_minus, b_plus, config.beta, config.sigma) :
                B_time(jump, ham_or_trott, b_minus, b_plus, config.t0, config.beta, config.sigma)
            copyto!(scratch.tmp2, B)
            scratch.tmp2 .*= jump_weight_scaling
            scratch.tmp2 .= 0.5 .* (scratch.tmp2 .+ scratch.tmp2')
            U_B = exp(-1im * config.delta * Hermitian(scratch.tmp2))
        end

        # Apply unitary evolution via U_B
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
    scratch.R .= 0.5 .* (scratch.R .+ scratch.R)

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

function run_thermalization_kraus(
    jumps::Vector{JumpOp},
    config::ThermalizeConfig,
    evolving_dm::Matrix{ComplexF64},
    hamiltonian::HamHam;
    trotter::Union{TrottTrott, Nothing}=nothing,
    rng::AbstractRNG = Random.default_rng(),
    rescale_by_inv_prob::Bool = false,
)

    dim = size(hamiltonian.data, 1)
    validate_config!(config)
    print_press(config)

    if config.domain isa TrotterDomain
        @assert trotter !== nothing
        ham_or_trott = trotter
        gibbs = Hermitian(trotter.eigvecs' * hamiltonian.eigvecs * hamiltonian.gibbs *
                          hamiltonian.eigvecs' * trotter.eigvecs)
    else
        ham_or_trott = hamiltonian
        gibbs = hamiltonian.gibbs
    end

    precomputed_data = precompute_data(config.domain, config, ham_or_trott)

    # precompute coherent U_B = exp(-i delta B(jump)) per jump to avoid allocations
    jump_weight_scaling = rescale_by_inv_prob ? (precomputed_data.gamma_norm_factor / jump_prob) : precomputed_data.gamma_norm_factor
    coherent_unitary_cache = if config.with_coherent
    [begin
        B = (config.domain isa TrotterDomain) ?
            B_trotter(j, ham_or_trott, precomputed_data.b_minus, precomputed_data.b_plus, config.beta, config.sigma) :
            B_time(j, ham_or_trott, precomputed_data.b_minus, precomputed_data.b_plus, config.t0, config.beta, config.sigma)
        B .= 0.5 .* (B .+ B')                  # hermitianize in-place
        B .*= jump_weight_scaling
        exp(-1im * config.delta * Hermitian(B))
    end for j in jumps]
else
    nothing
end

    scratch = KrausScratch(ComplexF64, dim)

    num_steps = Int(ceil(config.mixing_time / config.delta))

    convergence_cutoff = 1e-5
    distances_to_gibbs = [trace_distance_h(Hermitian(evolving_dm), gibbs)]

    p_jump = 1.0 / length(jumps)

    for step in 1:num_steps
        idx = rand(rng, 1:length(jumps))
        jump = jumps[idx]

        apply_kraus_step!(config.domain,
            evolving_dm,
            jump,
            ham_or_trott,
            config,
            precomputed_data,
            scratch;
            coherent_unitary_cache = (coherent_unitary_cache === nothing ? nothing : coherent_unitary_cache[idx]),
            jump_prob = p_jump,
            rescale_by_inv_prob = rescale_by_inv_prob,
        )

        dist = trace_distance_h(Hermitian(evolving_dm), gibbs)
        push!(distances_to_gibbs, dist)
        @printf("Dist to Gibbs: %s\n", dist)
        if dist < convergence_cutoff
            num_steps = step
            break
        end
    end

    time_steps = [0.0:config.delta:(num_steps * config.delta);]
    return HotAlgorithmResults(evolving_dm, distances_to_gibbs, time_steps, hamiltonian, trotter, config)
end

#TODO: Either fix this Kraus way to not trace drift, or incorporate explicit delta^2 terms from Chen's paper.