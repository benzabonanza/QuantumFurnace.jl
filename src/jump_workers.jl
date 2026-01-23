#* Liouvillian (vectorized) jump contributions
#TODO: use thread caches
function jump_contribution(::BohrDomain, 
    jump::JumpOp, 
    hamiltonian::HamHam, 
    config::LiouvConfig,
    precomputed_data)

    dim = size(hamiltonian.data, 1)
    unique_freqs = keys(hamiltonian.bohr_dict)
    (; alpha, gamma_norm_factor) = precomputed_data

    liouv_for_jump = zeros(ComplexF64, dim^2, dim^2)
    if config.with_coherent 
        coherent_term = coherent_bohr(hamiltonian, jump, config)
        rmul!(coherent_term, gamma_norm_factor) 
        vectorize_liouvillian_coherent!(liouv_for_jump, coherent_term)
    end

    alpha_A_nu1 = zeros(ComplexF64, dim, dim)
    for nu_2 in unique_freqs
        @. alpha_A_nu1 = alpha(hamiltonian.bohr_freqs, nu_2, config.beta, config.a, config.b) * jump.in_eigenbasis

        indices = hamiltonian.bohr_dict[nu_2]
        A_nu_2_vals = view(jump.in_eigenbasis, indices)

        # swapped for dagger
        rows_dag = getindex.(indices, 2) 
        cols_dag = getindex.(indices, 1)

        A_nu_2_dag = sparse(rows_dag, cols_dag, conj.(A_nu_2_vals), dim, dim)

        vectorize_liouv_diss_and_add!(liouv_for_jump, alpha_A_nu1, A_nu_2_dag, gamma_norm_factor)
    end

    return liouv_for_jump
end

function jump_contribution(::EnergyDomain, 
    jump::JumpOp, 
    hamiltonian::HamHam, 
    config::LiouvConfig, 
    precomputed_data)

    dim = size(hamiltonian.data, 1)
    (;transition, gamma_norm_factor, energy_labels) = precomputed_data

    liouv_for_jump = zeros(ComplexF64, dim^2, dim^2)
    if config.with_coherent 
        coherent_term = coherent_bohr(hamiltonian, jump, config) 
        rmul!(coherent_term, gamma_norm_factor)
        vectorize_liouvillian_coherent!(liouv_for_jump, coherent_term)
    end

    jump_oft = zeros(ComplexF64, dim, dim)
    prefactor = (config.w0 * config.beta / sqrt(2 * pi)) * gamma_norm_factor # w0, OFT norm^2, Fourier, norm factor
    for w in energy_labels
        oft!(jump_oft, jump, w, hamiltonian, config.beta) # subnorm = t0 * sqrt((sqrt(2 / pi)/beta) / (2 * pi))
        loop_scalar = prefactor * transition(w)

        vectorize_liouv_diss_and_add!(liouv_for_jump, jump_oft, loop_scalar)
    end
    return liouv_for_jump
end

function jump_contribution(::TimeDomain, 
    jump::JumpOp, 
    hamiltonian::HamHam, 
    config::LiouvConfig, 
    precomputed_data)

    dim = size(hamiltonian.data, 1)
    (; transition, gamma_norm_factor, energy_labels, oft_time_labels, b_minus, b_plus) = precomputed_data

    liouv_for_jump = zeros(ComplexF64, dim^2, dim^2)
    if config.with_coherent 
        coherent_term = B_time(jump, hamiltonian, b_minus, b_plus, config.t0, config.beta)
        rmul!(coherent_term, gamma_norm_factor) 
        vectorize_liouvillian_coherent!(liouv_for_jump, coherent_term)
    end

    jump_oft = zeros(ComplexF64, dim, dim)
    time_oft_caches = OFTCaches(dim)
    prefactor = config.w0 * config.t0^2 * (sqrt(2 / pi) / config.beta) / (2 * pi) * gamma_norm_factor
    for w in energy_labels
        time_oft_fast!(jump_oft, time_oft_caches, jump, w, hamiltonian, oft_time_labels, config.beta) # subnorm = t0 * sqrt((sqrt(2 / pi)/beta) / (2 * pi))
        loop_scalar = prefactor * transition(w)

        vectorize_liouv_diss_and_add!(liouv_for_jump, jump_oft, loop_scalar)
    end
    return liouv_for_jump
end

function jump_contribution(::TrotterDomain, 
    jump::JumpOp, 
    trotter::TrottTrott, 
    config::LiouvConfig, 
    precomputed_data)

    dim = size(trotter.eigvecs, 1)
    (;transition, gamma_norm_factor, energy_labels, oft_time_labels, b_minus, b_plus) = precomputed_data

    liouv_for_jump = zeros(ComplexF64, dim^2, dim^2)
    if config.with_coherent 
        coherent_term = B_trotter(jump, trotter, b_minus, b_plus, config.beta)
        rmul!(coherent_term, gamma_norm_factor) 
        vectorize_liouvillian_coherent!(liouv_for_jump, coherent_term)
    end

    jump_oft = zeros(ComplexF64, dim, dim)
    time_oft_caches = OFTCaches(dim)
    prefactor = config.w0 * config.t0^2 * (sqrt(2 / pi) / config.beta) / (2 * pi) * gamma_norm_factor # time ints t0^2, energy int w0, OFT time norm^2, Fourier
    for w in energy_labels
        trotter_oft_fast!(jump_oft, time_oft_caches, jump, w, trotter, oft_time_labels, config.beta) # subnorm = t0 * sqrt((sqrt(2 / pi)/beta) / (2 * pi))
        loop_scalar = prefactor * transition(w)

        vectorize_liouv_diss_and_add!(liouv_for_jump, jump_oft, loop_scalar)
    end
    return liouv_for_jump
end

#* Algorithmic jump contributions -----
function jump_contribution(::BohrDomain, 
    evolving_dm::Matrix{ComplexF64}, 
    jump::JumpOp, 
    hamiltonian::HamHam, 
    config::ThermalizeConfig,
    precomputed_data)

    (; alpha, gamma_norm_factor) = precomputed_data

    dim = size(evolving_dm, 1)
    unique_freqs = keys(hamiltonian.bohr_dict)

    jump_dm_contribution = zeros(ComplexF64, dim, dim)
    scaled_delta = config.delta * gamma_norm_factor  # easiest way to include gamma_norm_factor
    # Coherent part
    if config.with_coherent
        coherent_term = coherent_bohr(hamiltonian, jump, config)
        mul!(jump_dm_contribution, coherent_term, evolving_dm, -1im * scaled_delta, 1.0)
        mul!(jump_dm_contribution, evolving_dm, coherent_term, 1im * scaled_delta, 1.0)
    end

    # Dissipative part
    alpha_A_nu1 = zeros(ComplexF64, dim, dim)
    temp1 = similar(alpha_A_nu1)
    A_nu_2_dag_A_nu1 = similar(alpha_A_nu1)
    # mul!(C, A, B, α, β) computes C = α*A*B + β*C
    for nu_2 in unique_freqs
        @. alpha_A_nu1 = alpha(hamiltonian.bohr_freqs, nu_2, config.beta, config.a, config.b) * jump.in_eigenbasis
        
        indices = hamiltonian.bohr_dict[nu_2]
        A_nu_2_vals = view(jump.in_eigenbasis, indices)

        # swapped for dagger
        rows_dag = getindex.(indices, 2) 
        cols_dag = getindex.(indices, 1)

        A_nu_2_dag = sparse(rows_dag, cols_dag, conj.(A_nu_2_vals), dim, dim)

        mul!(A_nu_2_dag_A_nu1, A_nu_2_dag, alpha_A_nu1)

        # Term 1
        mul!(temp1, evolving_dm, A_nu_2_dag)
        mul!(jump_dm_contribution, alpha_A_nu1, temp1, scaled_delta, 1.0)

        # Term 2
        mul!(jump_dm_contribution, A_nu_2_dag_A_nu1, evolving_dm, -0.5 * scaled_delta, 1.0)

        # Term 3
        mul!(jump_dm_contribution, evolving_dm, A_nu_2_dag_A_nu1, -0.5 * scaled_delta, 1.0)
    end
        
    return jump_dm_contribution
end

function jump_contribution(::EnergyDomain, 
    evolving_dm::Matrix{ComplexF64}, 
    jump::JumpOp, 
    hamiltonian::HamHam, 
    config::ThermalizeConfig, 
    precomputed_data)

    (; transition, gamma_norm_factor, energy_labels) = precomputed_data

    dim = size(evolving_dm, 1)
    scaled_delta = config.delta * gamma_norm_factor

    jump_dm_contribution = zeros(ComplexF64, dim, dim)
    # Coherent part
    if config.with_coherent
        coherent_term = coherent_bohr(hamiltonian, jump, config)
        mul!(jump_dm_contribution, coherent_term, evolving_dm, -1im * scaled_delta, 1.0)
        mul!(jump_dm_contribution, evolving_dm, coherent_term, 1im * scaled_delta, 1.0)
    end

    # Dissipative part
    prefactor = scaled_delta * config.w0 * config.beta / sqrt(2 * pi)
    jump_oft = zeros(ComplexF64, dim, dim)
    jump_dag_jump = similar(jump_oft)
    temp1 = similar(jump_oft)
    # mul!(C, A, B, α, β) computes C = α*A*B + β*C
    for w in energy_labels
        oft!(jump_oft, jump, w, hamiltonian, config.beta)
        mul!(jump_dag_jump, jump_oft', jump_oft)

        loop_factor = transition(w) * prefactor
        # Term 1
        mul!(temp1, evolving_dm, jump_oft')  # rho * A'
        mul!(jump_dm_contribution, jump_oft, temp1, loop_factor, 1.0)  # L += prefactor * A * rho * A'

        # Term 2
        mul!(jump_dm_contribution, jump_dag_jump, evolving_dm, -0.5 * loop_factor, 1.0)
        
        # Term 3
        mul!(jump_dm_contribution, evolving_dm, jump_dag_jump, -0.5 * loop_factor, 1.0)
    end
    
    return jump_dm_contribution
end

function jump_contribution(
    ::TimeDomain, 
    evolving_dm::Matrix{ComplexF64}, 
    jump::JumpOp, 
    hamiltonian::HamHam, 
    config::ThermalizeConfig,
    precomputed_data)

    dim = size(evolving_dm, 1)
    (; transition, gamma_norm_factor, energy_labels, oft_time_labels, b_minus, b_plus) = precomputed_data

    scaled_delta = config.delta * gamma_norm_factor
    jump_dm_contribution = zeros(ComplexF64, dim, dim)
    # Coherent part
    if config.with_coherent
        # coherent_term_f = coherent_term_time(jump, hamiltonian, f_minus, f_plus, t0)
        coherent_term = B_time(jump, hamiltonian, b_minus, b_plus, config.t0, config.beta)
        mul!(jump_dm_contribution, coherent_term, evolving_dm, -1im * scaled_delta, 1.0)
        mul!(jump_dm_contribution, evolving_dm, coherent_term, 1im * scaled_delta, 1.0)
    end

    jump_oft = zeros(ComplexF64, dim, dim)
    jump_dag_jump = similar(jump_oft)
    temp1 = similar(jump_oft)

    # Pre-allocate caches for the time_oft function as well
    oft_caches = OFTCaches(dim)
    prefactor = scaled_delta * config.w0 * config.t0^2 * (sqrt(2 / pi) / config.beta) / (2 * pi)
    for w in energy_labels
        time_oft_fast!(jump_oft, oft_caches, jump, w, hamiltonian, oft_time_labels, config.beta)
        
        # jump_dag_jump = jump_oft' * jump_oft
        mul!(jump_dag_jump, jump_oft', jump_oft)

        loop_factor = transition(w) * prefactor
        # Term 1
        mul!(temp1, evolving_dm, jump_oft')  # rho * A'
        mul!(jump_dm_contribution, jump_oft, temp1, loop_factor, 1.0)  # L += prefactor * A * rho * A'

        # Term 2
        mul!(jump_dm_contribution, jump_dag_jump, evolving_dm, -0.5 * loop_factor, 1.0)

        # Term 3
        mul!(jump_dm_contribution, evolving_dm, jump_dag_jump, -0.5 * loop_factor, 1.0)
    end
    return jump_dm_contribution
end

function jump_contribution(::TrotterDomain, 
    evolving_dm::Matrix{ComplexF64}, 
    jump::JumpOp, 
    trotter::TrottTrott, 
    config::ThermalizeConfig,
    precomputed_data)

    dim = size(evolving_dm, 1)
    (; transition, gamma_norm_factor, energy_labels, oft_time_labels, b_minus, b_plus) = precomputed_data

    scaled_delta = config.delta * gamma_norm_factor

    jump_dm_contribution = zeros(ComplexF64, dim, dim)
    # Coherent part
    if config.with_coherent
        # coherent_term_f = coherent_term_trotter(jump, trotter, f_minus, f_plus)
        coherent_term = B_trotter(jump, trotter, b_minus, b_plus, beta)
        mul!(jump_dm_contribution, coherent_term, evolving_dm, -1im * scaled_delta, 1.0)
        mul!(jump_dm_contribution, evolving_dm, coherent_term, 1im * scaled_delta, 1.0)
    end

    # Pre-allocate caches
    jump_oft = zeros(ComplexF64, dim, dim)
    jump_dag_jump = similar(jump_oft)
    temp1 = similar(jump_oft)

    oft_caches = OFTCaches(dim)
    prefactor = scaled_delta * config.w0 * config.t0^2 * (sqrt(2 / pi) / config.beta) / (2 * pi)
    for w in energy_labels
        trotter_oft_fast!(jump_oft, oft_caches, jump, w, trotter, oft_time_labels, config.beta)
        
        # jump_dag_jump = jump_oft' * jump_oft
        mul!(jump_dag_jump, jump_oft', jump_oft)

        loop_factor = transition(w) * prefactor
        # Term 1
        mul!(temp1, evolving_dm, jump_oft')  # rho * A'
        mul!(jump_dm_contribution, jump_oft, temp1, loop_factor, 1.0)  # L += prefactor * A * rho * A'

        # Term 2
        mul!(jump_dm_contribution, jump_dag_jump, evolving_dm, -0.5 * loop_factor, 1.0)

        # Term 3
        mul!(jump_dm_contribution, evolving_dm, jump_dag_jump, -0.5 * loop_factor, 1.0)
    end

    return jump_dm_contribution
end

#* Precompute B, R for Quantum Trajector & efficient action of L(X) -----
# No Bohr version, because we can only diagonalize the Kossakowski matrix up to very low number of qubits.
# At n = 10, we have 10^6 x 10^6 Kossakowski.

function precompute_kraus_jumps(
    ::EnergyDomain,
    jumps::Vector{JumpOp}, 
    hamiltonian::HamHam,
    config::LiouvConfig,
    precomputed_data,
    _
    )

    dim = size(hamiltonian.data, 1)
    (; transition, gamma_norm_factor, energy_labels) = precomputed_data


    misc_factors = config.w0 * config.beta / sqrt(2 * pi)
    transition_rates = transition.(energy_labels) .* gamma_norm_factor
    
    kraus_jump_rates = @. sqrt(misc_factors * transition_rates)

    num_kraus_jumps = length(energy_labels) * length(jumps)
    kraus_jumps = Vector{Matrix{ComplexF64}}(undef, num_kraus_jumps)
    k = 1
    for jump in jumps
        for (i, w) in enumerate(energy_labels)
            kraus_jump = Matrix{ComplexF64}(undef, dim, dim)
            oft!(kraus_jump, jump, w, hamiltonian, config.beta)
            rmul!(kraus_jump, kraus_jump_rates[i])
            kraus_jumps[k] = kraus_jump
            k += 1
        end
    end
    return kraus_jumps
end

function precompute_kraus_jumps(
    ::TimeDomain,
    jumps::Vector{JumpOp}, 
    hamiltonian::HamHam,
    config::LiouvConfig,
    precomputed_data,
    oft_caches
    )

    dim = size(hamiltonian.data, 1)
    (; transition, gamma_norm_factor, b_minus, b_plus, energy_labels, oft_time_labels) = precomputed_data

    misc_factors = config.w0 * config.t0^2 * (sqrt(2 / pi) / config.beta) / (2 * pi)
    transition_rates = transition.(energy_labels) .* gamma_norm_factor
    println("MAXIMUM OF GAMMA")
    println(maximum(transition_rates))
    
    kraus_jump_rates = @. sqrt(misc_factors * transition_rates)

    num_kraus_jumps = length(energy_labels) * length(jumps)
    kraus_jumps = Vector{Matrix{ComplexF64}}(undef, num_kraus_jumps)
    k = 1
    for jump in jumps
        for (i, w) in enumerate(energy_labels)
            kraus_jump = Matrix{ComplexF64}(undef, dim, dim)
            time_oft_fast!(kraus_jump, oft_caches, jump, w, hamiltonian, oft_time_labels, config.beta)
            rmul!(kraus_jump, kraus_jump_rates[i])
            kraus_jumps[k] = kraus_jump
            k += 1
        end
    end
    return kraus_jumps
end

function precompute_kraus_jumps(
    ::TrotterDomain,
    jumps::Vector{JumpOp}, 
    trotter::TrottTrott,
    config::LiouvConfig,
    precomputed_data,
    oft_caches
    )

    dim = size(trotter.eigvecs, 1)
    (; transition, gamma_norm_factor, b_minus, b_plus, energy_labels, oft_time_labels) = precomputed_data

    misc_factors = config.w0 * config.t0^2 * (sqrt(2 / pi) / config.beta) / (2 * pi)
    transition_rates = transition.(energy_labels) .* gamma_norm_factor
    
    kraus_jump_rates = @. sqrt(misc_factors * transition_rates)

    num_kraus_jumps = length(energy_labels) * length(jumps)
    kraus_jumps = Vector{Matrix{ComplexF64}}(undef, num_kraus_jumps)
    # temp_buffer = Matrix{ComplexF64}(undef, dim, dim)

    k = 1
    for jump in jumps
        for (i, w) in enumerate(energy_labels)
            kraus_jump = Matrix{ComplexF64}(undef, dim, dim)

            trotter_oft_fast!(kraus_jump, oft_caches, jump, w, trotter, oft_time_labels, config.beta)
            # temp_buffer = OFT in Trotter basis:
            # mul!(oft_caches.temp_op, temp_buffer, trotter.trafo_from_eigen_to_trotter)  # temp_op = OFT U
            # mul!(kraus_jump, trotter.trafo_from_eigen_to_trotter', oft_caches.temp_op)  # kraus_jump = OFT in eigenbasis H

            rmul!(kraus_jump, kraus_jump_rates[i])
            kraus_jumps[k] = kraus_jump
            k += 1
        end
    end
    return kraus_jumps  # in Trotter basis
end

"""
R = \\sum_k L_k^\\dagger L_k, in eigenbasis of H.
"""
function precompute_R(kraus_jumps::Vector{Matrix{ComplexF64}})

    dim = size(kraus_jumps[1], 1)
    R = zeros(ComplexF64, dim, dim)
    # herk! only updates the upper triangle
    for kraus_jump in kraus_jumps  
        BLAS.herk!('U', 'C', 1.0, kraus_jump, 1.0, R)  # R += L' L
    end

    # Symmetrize it to full Hermitian matrix
    LinearAlgebra.copytri!(R, 'U', true)
    return R
end

# Should match up but obviously, not "faithful" since the QC have to go from the time domain with a Fourier transform.
function precompute_B(R::Matrix{ComplexF64}, hamiltonian::HamHam, beta::Float64)
    return @. R * 1im * tanh(beta * hamiltonian.bohr_freqs / 4) / 2
end

struct LindbladWorkspace{T}
    temp_buffers::Vector{Matrix{T}} 
    accumulators::Vector{Matrix{T}} 
    
    function LindbladWorkspace(dim::Int, T=ComplexF64)
        temp = [Matrix{T}(undef, dim, dim) for _ in 1:nthreads()]
        acc  = [Matrix{T}(undef, dim, dim) for _ in 1:nthreads()]
        new{T}(temp, acc)
    end
end

#TODO: test it; set BLAS threads to 1, let julia threads be more. 
#* In-place Lindbladian action at once
"""
Applies L(X) = -i[B, X] - 0.5{R, X} + sum(L X L')
"""
function apply_lindbladian!(
    target::Matrix{ComplexF64},
    rho::Matrix{ComplexF64},
    kraus_jumps::Vector{Matrix{ComplexF64}}, 
    B::Matrix{ComplexF64},
    R::Matrix{ComplexF64},
    ws::LindbladWorkspace{ComplexF64})
    
    #  mul!(C, A, B, alpha, beta) = alpha A B + beta C
    # Coherent: -i [B, rho]
    mul!(target, B, rho, -1im, 0.0) # target <- -i * B * rho
    mul!(target, rho, B, 1im, 1.0)  # target += i * rho * B
    
    # Reflection: -0.5 {R, rho}
    mul!(target, R, rho, -0.5, 1.0)  # target -= 0.5 * R * rho
    mul!(target, rho, R, -0.5, 1.0)  # target -= 0.5 * rho * R

    # Reset thread accumulators to zero
    @threads for i in 1:nthreads()
        fill!(ws.accumulators[i], 0.0)
    end

    # Translation: sum(L rho L') (Multi-threaded)
    @threads for k in 1:length(kraus_jumps)
        id = threadid()
        L = kraus_jumps[k]
        
        buf = ws.temp_buffers[id]
        acc = ws.accumulators[id]
        
        # Compute L * rho * L' 
        mul!(buf, rho, L')  # buf = rho * L'
        mul!(acc, L, buf, 1.0, 1.0)  # acc += L * rho * L'
    end

    # Sum all thread accumulators
    for i in 1:nthreads()
        @. target += ws.accumulators[i]
    end
    
    return target
end

"""
Applies L^dagger(X) = +i[B, X] - 0.5{R, X} + sum(L' X L),
i.e. the Hilbert-Schmidt adjoint.
"""
function apply_lindbladian_dagger!(
    target::Matrix{ComplexF64},
    rho::Matrix{ComplexF64},
    kraus_jumps::Vector{Matrix{ComplexF64}}, 
    B::Matrix{ComplexF64},
    R::Matrix{ComplexF64},
    ws::LindbladWorkspace{ComplexF64})
    
    # Coherent: -i [B, rho]
    mul!(target, B, rho, 1im, 0.0)  # target <- i * B * rho
    mul!(target, rho, B, -1im, 1.0)  # target += -i * rho * B
    
    # Reflection: -0.5 {R, rho}
    mul!(target, R, rho, -0.5, 1.0)  # target -= 0.5 * R * rho
    mul!(target, rho, R, -0.5, 1.0)  # target -= 0.5 * rho * R

    # Reset thread accumulators to zero
    @threads for i in 1:nthreads()
        fill!(ws.accumulators[i], 0.0)
    end

    # Translation: sum(L rho L') (Multi-threaded)
    @threads for k in 1:length(kraus_jumps)
        id = threadid()
        L = kraus_jumps[k]
        
        buf = ws.temp_buffers[id]
        acc = ws.accumulators[id]
        
        # Compute L * rho * L' 
        mul!(buf, rho, L)  # buf = rho * L
        mul!(acc, L', buf, 1.0, 1.0)  # acc += L' * rho * L
    end

    # Sum all thread accumulators
    for i in 1:nthreads()
        @. target += ws.accumulators[i]
    end
    
    return target
end

#TODO: Modernize these, but only after thread caches implementation above
#* In-place Lindbladian action jump-by-jump
function jump_contribution!(
    target::AbstractMatrix{ComplexF64}, 
    ::BohrDomain, 
    rho::AbstractMatrix{ComplexF64}, 
    jump::JumpOp, 
    hamiltonian::HamHam,
    config::LiouvConfig,
    precomputed_data,
    caches  # (jump_caches, oft_caches)
    )

    (; jump_caches, oft_caches) = caches

    alpha = pick_alpha(config)
    unique_freqs = keys(hamiltonian.bohr_dict)

    # Coherent part
    if config.with_coherent
        coherent_term = coherent_bohr(hamiltonian, jump, config)
        mul!(target, coherent_term, rho, -1im, 1.0)
        mul!(target, rho, coherent_term, 1im, 1.0)
    end

    # Dissipative part
    # mul!(C, A, B, α, β) computes C = α*A*B + β*C
    for nu_2 in unique_freqs
        @. jump_caches.jump_1 = alpha(hamiltonian.bohr_freqs, nu_2, config.beta, config.a, config.b) * jump.in_eigenbasis

        A_nu_2::SparseMatrixCSC{ComplexF64} = spzeros(dim, dim)
        indices = hamiltonian.bohr_dict[nu_2]
        A_nu_2[indices] .= jump.in_eigenbasis[indices]

        mul!(jump_caches.jump_2_dag_jump_1, A_nu_2', jump_caches.jump_1)

        # Term 1
        mul!(jump_caches.temp1, rho, A_nu_2')
        mul!(target, jump_caches.jump_1, jump_caches.temp1, 1.0, 1.0)

        # Term 2
        mul!(target, jump_caches.jump_2_dag_jump_1, rho, -0.5, 1.0)

        # Term 3
        mul!(target, rho, jump_caches.jump_2_dag_jump_1, -0.5, 1.0)
    end
        
    return target
end

function jump_contribution!(
    target::AbstractMatrix{ComplexF64}, 
    ::EnergyDomain, 
    input::AbstractMatrix{ComplexF64}, 
    jump::JumpOp, 
    hamiltonian::HamHam,
    config::LiouvConfig,
    precomputed_data,
    caches
    )

    (; w0, t0, transition, f_minus, f_plus, energy_labels, oft_time_labels) = precomputed_data
    (; jump_caches, oft_caches) = caches

    # Coherent part
    if config.with_coherent
        coherent_term = coherent_bohr(hamiltonian, jump, config)
        mul!(target, coherent_term, input, -1im, 1.0)
        mul!(target, input, coherent_term, 1im, 1.0)
    end

    # Dissipative part
    prefactor = w0 * config.beta / sqrt(2 * pi)
    # mul!(C, A, B, α, β) computes C = α*A*B + β*C
    for w in energy_labels
        oft!(jump_caches.jump_1, jump, w, hamiltonian, config.beta)
        mul!(jump_caches.jump_2_dag_jump_1, jump_caches.jump_1', jump_caches.jump_1)

        loop_factor = transition(w) * prefactor
        # Term 1
        mul!(jump.temp1, input, jump_caches.jump_1')  # rho * A'
        mul!(target, jump_caches.jump_1, jump.temp1, loop_factor, 1.0)  # L += prefactor * A * rho * A'

        # Term 2
        mul!(target, jump_caches.jump_2_dag_jump_1, input, -0.5 * loop_factor, 1.0)
        
        # Term 3
        mul!(target, input, jump_caches.jump_2_dag_jump_1, -0.5 * loop_factor, 1.0)
    end
    
    return target
end

function jump_contribution!(
    target::AbstractMatrix{ComplexF64}, 
    ::TimeDomain, 
    input::AbstractMatrix{ComplexF64}, 
    jump::JumpOp, 
    hamiltonian::HamHam,
    config::LiouvConfig,
    precomputed_data,
    caches
    )

    (; w0, t0, transition, f_minus, f_plus, energy_labels, oft_time_labels) = precomputed_data
    (; jump_caches, oft_caches) = caches

    # Coherent part
    if config.with_coherent
        coherent_term = coherent_term_time(jump, hamiltonian, f_minus, f_plus, t0)
        mul!(target, coherent_term, input, -1im, 1.0)
        mul!(target, input, coherent_term, 1im, 1.0)
    end

    prefactor = w0 * t0^2 * (sqrt(2 / pi) / config.beta) / (2 * pi)
    for w in energy_labels
        time_oft_fast!(jump_caches.jump_1, oft_caches, jump, w, hamiltonian, oft_time_labels, config.beta)
        
        # jump_dag_jump = jump_oft' * jump_oft
        mul!(jump_caches.jump_2_dag_jump_1, jump_caches.jump_1', jump_caches.jump_1)

        loop_factor = transition(w) * prefactor
        # Term 1
        mul!(jump_caches.temp1, input, jump_caches.jump_1')  # rho * A'
        mul!(target, jump_caches.jump_1, jump_caches.temp1, loop_factor, 1.0)  # L += prefactor * A * rho * A'

        # Term 2
        mul!(target, jump_caches.jump_2_dag_jump_1, input, -0.5 * loop_factor, 1.0)

        # Term 3
        mul!(target, input, jump_caches.jump_2_dag_jump_1, -0.5 * loop_factor, 1.0)
    end
    return target
end

function jump_contribution!(
    target::AbstractMatrix{ComplexF64}, 
    ::TrotterDomain, 
    input::AbstractMatrix{ComplexF64}, 
    jump::JumpOp, 
    trotter::TrottTrott,
    config::LiouvConfig,
    precomputed_data,
    caches
    )

    (; w0, t0, transition, f_minus, f_plus, energy_labels, oft_time_labels) = precomputed_data
    (; jump_caches, oft_caches) = caches

    # Coherent part
    if config.with_coherent
        coherent_term = coherent_term_trotter(jump, trotter, f_minus, f_plus)
        mul!(target, coherent_term, input, -1im, 1.0)
        mul!(target, input, coherent_term, 1im, 1.0)
    end

    prefactor = w0 * t0^2 * (sqrt(2 / pi) / config.beta) / (2 * pi)
    for w in energy_labels
        trotter_oft_fast!(jump_caches.jump_1, oft_caches, jump, w, trotter, oft_time_labels, config.beta)
        
        # jump_dag_jump = jump_oft' * jump_oft
        mul!(jump_caches.jump_2_dag_jump_1, jump_caches.jump_1', jump_caches.jump_1)

        loop_factor = transition(w) * prefactor
        # Term 1
        mul!(jump_caches.temp1, input, jump_caches.jump_1')  # rho * A'
        mul!(target, jump_caches.jump_1, jump_caches.temp1, loop_factor, 1.0)  # L += prefactor * A * rho * A'

        # Term 2
        mul!(target, jump_caches.jump_2_dag_jump_1, input, -0.5 * loop_factor, 1.0)

        # Term 3
        mul!(target, input, jump_caches.jump_2_dag_jump_1, -0.5 * loop_factor, 1.0)
    end

    return target
end



#* Slow and old
# function jump_contribution_slow(::TrotterDomain, jump::JumpOp, trotter::TrottTrott, config::LiouvConfig, 
#     energy_labels::Vector{Float64}, time_labels::Vector{Float64})

#     dim = size(trotter.eigvecs, 1)
#     w0 = abs(energy_labels[2] - energy_labels[1])
#     oft_time_labels = truncate_time_labels_for_oft(time_labels, config.beta)

#     transition = pick_transition(config.beta, config.a, config.b, config.with_linear_combination)

#     if config.with_coherent
#         f_minus = compute_truncated_f(compute_f_minus, time_labels, config.beta)
#         if config.with_linear_combination  
#             if config.a != 0.0  # Improved Metro / Glauber
#                 f_plus = compute_truncated_f(compute_f_plus_eh, time_labels, config.beta, config.a, config.b)
#             else  # Metro
#                 f_plus = compute_truncated_f(compute_f_plus_metro, time_labels, config.beta, config.eta)
#             end
#         else  # Gaussian
#             f_plus = compute_truncated_f(compute_f_plus, time_labels, config.beta)
#         end
#     end

#     liouv_coherent_part_for_jump = zeros(ComplexF64, dim^2, dim^2)
#     liouv_diss_part_for_jump = zeros(ComplexF64, dim^2, dim^2)
#     if config.with_coherent  # There is no energy formulation of the coherent term, only Bohr and time.
#         coherent_term = coherent_term_trotter(jump, trotter, f_minus, f_plus)
#         # coherent_term = trotter.trafo_from_eigen_to_trotter' * coherent_term * trotter.trafo_from_eigen_to_trotter
#         liouv_coherent_part_for_jump .+= vectorize_liouvillian_coherent(coherent_term)
#     end

#     for w in energy_labels
#         jump_oft = trotter_oft(jump, w, trotter, oft_time_labels, config.beta) # t0 * sqrt((sqrt(2 / pi)/beta) / (2 * pi))
#         # jump_oft = trotter.trafo_from_eigen_to_trotter' * jump_oft * trotter.trafo_from_eigen_to_trotter
#         liouv_diss_part_for_jump .+= transition(w) * vectorize_liouvillian_diss(jump_oft)
#     end
    
#     prefactor = w0 * trotter.t0^2 * (sqrt(2 / pi) / config.beta) / (2 * pi)
#     return liouv_coherent_part_for_jump .+ prefactor * liouv_diss_part_for_jump  #! L in trotter basis
# end

# function jump_contribution_slow(::TimeDomain, jump::JumpOp, hamiltonian::HamHam, config::LiouvConfig, 
#     energy_labels::Vector{Float64}, time_labels::Vector{Float64})

#     dim = size(hamiltonian.data, 1)
#     w0 = abs(energy_labels[2] - energy_labels[1])
#     t0 = time_labels[2] - time_labels[1]
#     oft_time_labels = truncate_time_labels_for_oft(time_labels, config.beta)

#     transition = pick_transition(config.beta, config.a, config.b, config.with_linear_combination)

#     if config.with_coherent
#         f_minus = compute_truncated_f(compute_f_minus, time_labels, config.beta)
#         if config.with_linear_combination  
#             if config.a != 0.0  # Improved Metro / Glauber
#                 f_plus = compute_truncated_f(compute_f_plus_eh, time_labels, config.beta, config.a, config.b)
#             else  # Metro
#                 f_plus = compute_truncated_f(compute_f_plus_metro, time_labels, config.beta, config.eta)
#             end
#         else  # Gaussian
#             f_plus = compute_truncated_f(compute_f_plus, time_labels, config.beta)
#         end
#     end

#     liouv_coherent_part_for_jump = zeros(ComplexF64, dim^2, dim^2)
#     liouv_diss_part_for_jump = zeros(ComplexF64, dim^2, dim^2)
#     if config.with_coherent 
#         coherent_term = coherent_term_time(jump, hamiltonian, f_minus, f_plus, t0)  
#         liouv_coherent_part_for_jump .+= vectorize_liouvillian_coherent(coherent_term)
#     end

#     for w in energy_labels
#         jump_oft = time_oft(jump, w, hamiltonian, oft_time_labels, config.beta) # subnorm = t0 * sqrt((sqrt(2 / pi)/beta) / (2 * pi))
#         liouv_diss_part_for_jump .+= transition(w) * vectorize_liouvillian_diss(jump_oft)
#     end
#     prefactor = w0 * t0^2 * (sqrt(2 / pi) / config.beta) / (2 * pi)  # time ints t0^2, energy int w0, OFT time norm^2, Fourier
#     return liouv_coherent_part_for_jump .+ prefactor * liouv_diss_part_for_jump
# end

# function jump_contribution_slow(::EnergyDomain, jump::JumpOp, hamiltonian::HamHam, config::LiouvConfig, 
#     energy_labels::Vector{Float64})

#     dim = size(hamiltonian.data, 1)
#     w0 = abs(energy_labels[2] - energy_labels[1])

#     transition = pick_transition(config.beta, config.a, config.b, config.with_linear_combination)

#     liouv_coherent_part_for_jump = zeros(ComplexF64, dim^2, dim^2)
#     liouv_diss_part_for_jump = zeros(ComplexF64, dim^2, dim^2)
#     if config.with_coherent
#         coherent_term = coherent_bohr(hamiltonian, jump, config)
#         liouv_coherent_part_for_jump .+= vectorize_liouvillian_coherent(coherent_term)
#     end

#     for w in energy_labels
#         jump_oft = oft(jump, w, hamiltonian, config.beta)
#         liouv_diss_part_for_jump .+= transition(w) * vectorize_liouvillian_diss(jump_oft)
#     end
#     oft_norm_squared = config.beta / sqrt(2 * pi)
#     return liouv_coherent_part_for_jump .+ w0 * oft_norm_squared * liouv_diss_part_for_jump
# end

# function jump_contribution_slow(::BohrDomain, jump::JumpOp, hamiltonian::HamHam, config::LiouvConfig)
#     dim = size(hamiltonian.data, 1)
#     unique_freqs = keys(hamiltonian.bohr_dict)

#     alpha = pick_alpha(config)

#     liouv_for_jump = zeros(ComplexF64, dim^2, dim^2)  
#     # Coherent part
#     if config.with_coherent
#         coherent_term = coherent_bohr(hamiltonian, jump, config)
#         liouv_for_jump .+= vectorize_liouvillian_coherent(coherent_term)
#     end

#     # Dissipative part
#     for nu_2 in unique_freqs
#         alpha_nu1_matrix = alpha.(hamiltonian.bohr_freqs, nu_2, config.beta, config.a, config.b)

#         A_nu_2::SparseMatrixCSC{ComplexF64} = spzeros(dim, dim)
#         A_nu_2[hamiltonian.bohr_dict[nu_2]] .= jump.in_eigenbasis[hamiltonian.bohr_dict[nu_2]]

#         liouv_for_jump .+= vectorize_liouvillian_diss(alpha_nu1_matrix .* jump.in_eigenbasis, A_nu_2')
#     end
#     return liouv_for_jump
# end

# function jump_contribution_slow(::BohrDomain, evolving_dm::Matrix{ComplexF64}, jump::JumpOp, hamiltonian::HamHam, 
#     config::ThermalizeConfig)

#     dim = size(evolving_dm, 1)
#     alpha = pick_alpha(config)
#     unique_freqs = keys(hamiltonian.bohr_dict)

#     jump_coherent = zeros(ComplexF64, dim, dim)
#     jump_dissipative = zeros(ComplexF64, dim, dim)
#     # Coherent part
#     if config.with_coherent
#         coherent_term = coherent_bohr(hamiltonian, jump, config)
#         jump_coherent .+= - 1im * (coherent_term * evolving_dm - evolving_dm * coherent_term)
#     end

#     # Dissipative part
#     for nu_2 in unique_freqs
#         A_nu_2::SparseMatrixCSC{ComplexF64} = spzeros(dim, dim)
#         indices_nu_2 = hamiltonian.bohr_dict[nu_2]
#         A_nu_2[indices_nu_2] .= jump.in_eigenbasis[indices_nu_2]

#         alpha_A_nu1 = alpha.(hamiltonian.bohr_freqs, nu_2, config.beta, config.a, config.b) .* jump.in_eigenbasis

#         jump_dissipative .+= (alpha_A_nu1 * evolving_dm * A_nu_2' - 0.5 * (A_nu_2' * alpha_A_nu1 * evolving_dm 
#                                         + evolving_dm * A_nu_2' * alpha_A_nu1)
#                                         )
#     end
#     return config.delta * (jump_coherent + jump_dissipative)
# end

# function jump_contribution_slow(::EnergyDomain, evolving_dm::Matrix{ComplexF64}, jump::JumpOp, hamiltonian::HamHam, 
#     config::ThermalizeConfig)

#     dim = size(evolving_dm, 1)
#     w0 = abs(energy_labels[2] - energy_labels[1])

#     jump_coherent = zeros(ComplexF64, dim, dim)
#     jump_dissipative = zeros(ComplexF64, dim, dim)
#     # Coherent part
#     if config.with_coherent
#         coherent_term = coherent_bohr(hamiltonian, jump, config)
#         jump_coherent .+= - 1im * (coherent_term * evolving_dm - evolving_dm * coherent_term)
#     end

#     # Dissipative part
#     for w in energy_labels
#         jump_oft = oft(jump, w, hamiltonian, config.beta)
#         jump_dag_jump = jump_oft' * jump_oft
#         jump_dissipative .+= transition(w) * (
#             jump_oft * evolving_dm * jump_oft' - 0.5 * (jump_dag_jump * evolving_dm + evolving_dm * jump_dag_jump)
#             )
#     end

#     oft_prefactor = config.beta / sqrt(2 * pi)
#     return config.delta * (jump_coherent + w0 * oft_prefactor * jump_dissipative)
# end

# function jump_contribution_slow(::TimeDomain, evolving_dm::Matrix{ComplexF64}, jump::JumpOp, hamiltonian::HamHam, 
#     config::ThermalizeConfig)

#     dim = size(evolving_dm, 1)
#     w0 = abs(energy_labels[2] - energy_labels[1])
#     t0 = time_labels[2] - time_labels[1]
#     oft_time_labels = truncate_time_labels_for_oft(time_labels, config.beta)

#     jump_coherent = zeros(ComplexF64, dim, dim)
#     jump_dissipative = zeros(ComplexF64, dim, dim)

#     # Coherent part
#     if config.with_coherent
#         coherent_term = coherent_term_time(jump, hamiltonian, f_minus, f_plus, t0)
#         jump_coherent .+= - 1im * (coherent_term * evolving_dm - evolving_dm * coherent_term)
#     end

#     # Dissipative part
#     for w in energy_labels
#         jump_oft = time_oft(jump, w, hamiltonian, oft_time_labels, config.beta)
#         jump_dag_jump = jump_oft' * jump_oft
#         jump_dissipative .+= transition(w) * (
#             jump_oft * evolving_dm * jump_oft' - 0.5 * (jump_dag_jump * evolving_dm + evolving_dm * jump_dag_jump)
#             )
#     end

#     oft_prefactor = (sqrt(2 / pi) / config.beta) / (2 * pi)
#     return config.delta * (jump_coherent + w0 * t0^2 * oft_prefactor * jump_dissipative)
# end

# function jump_contribution_slow(::TrotterDomain, evolving_dm::Matrix{ComplexF64}, jump::JumpOp, trotter::TrottTrott, 
#     config::ThermalizeConfig)

#     dim = size(evolving_dm, 1)
#     w0 = abs(energy_labels[2] - energy_labels[1])
#     oft_time_labels = truncate_time_labels_for_oft(time_labels, config.beta)

#     jump_coherent = zeros(ComplexF64, dim, dim)
#     jump_dissipative = zeros(ComplexF64, dim, dim)

#     # Coherent part
#     if config.with_coherent
#         coherent_term = coherent_term_trotter(jump, trotter, f_minus, f_plus)
#         jump_coherent .+= - 1im * (coherent_term * evolving_dm - evolving_dm * coherent_term)
#     end

#     # Dissipative part
#     for w in energy_labels
#         jump_oft = trotter_oft(jump, w, trotter, oft_time_labels, config.beta) # t0 * sqrt((sqrt(2 / pi)/beta) / (2 * pi))
#         jump_dag_jump = jump_oft' * jump_oft
#         jump_dissipative .+= transition(w) * (
#             jump_oft * evolving_dm * jump_oft' - 0.5 * (jump_dag_jump * evolving_dm + evolving_dm * jump_dag_jump)
#             )
#     end

#     oft_prefactor = (sqrt(2 / pi) / config.beta) / (2 * pi)
#     return config.delta * (jump_coherent + w0 * trotter.t0^2 * oft_prefactor * jump_dissipative)
# end

