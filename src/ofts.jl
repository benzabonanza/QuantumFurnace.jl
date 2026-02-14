function oft(jump::JumpOp, energy::Float64, hamiltonian::HamHam, sigma::Float64)
    """Subnormalized, multiply by sqrt(1 / sigma sqrt(2 * pi))"""
    return @. jump.in_eigenbasis * exp(-(energy - hamiltonian.bohr_freqs)^2 / (4 * sigma^2))
end

function oft!(out_matrix::Matrix{ComplexF64}, jump::JumpOp, energy::Float64, hamiltonian::HamHam, sigma::Float64)
    """Subnormalized, multiply by sqrt(1 / sigma sqrt(2 * pi))"""
    @. out_matrix = jump.in_eigenbasis * exp(-(energy - hamiltonian.bohr_freqs)^2 / (4 * sigma^2)) 
    return out_matrix
end

# Depricated but used for tests. We use precomputed NUFFT prefactors now in jump_contributions!()
function time_oft!(
    out_matrix::Matrix{ComplexF64},
    caches::OFTCaches,
    jump::JumpOp,
    energy::Float64,
    hamiltonian::HamHam, 
    time_labels::Vector{Float64},
    sigma::Float64
    )
    
    # Ensure the prefactor cache is the right size
    if length(caches.prefactors) != length(time_labels)
        resize!(caches.prefactors, length(time_labels))
    end
    
    # In-place calculation of prefactors
    @fastmath @. caches.prefactors = exp(-time_labels^2 * sigma^2 - 1im * energy * time_labels)
    
    zero_index = findfirst(t -> t >= -1e-12, time_labels)

    # Zero out the output matrix before we start accumulating
    fill!(out_matrix, 0.0)
    
    # --- Re-use the cache matrices U and temp_op inside the loops ---
    if jump.orthogonal  # Orthogonal (X, Z)
        # t = 0.0 case: U = I
        @. out_matrix += caches.prefactors[zero_index] * jump.in_eigenbasis

        for i in (zero_index + 1):length(time_labels)
            t = time_labels[i]
            @fastmath caches.U.diag .= exp.(1im .* hamiltonian.eigvals .* t)
            
            copyto!(caches.temp_op, jump.in_eigenbasis)
            # temp_op = U*jump*U', for diagonal U's:
            caches.temp_op .*= (caches.U.diag * caches.U.diag')
            
            @. out_matrix += caches.prefactors[i] * caches.temp_op
            @. out_matrix += conj(caches.prefactors[i]) * $(transpose(caches.temp_op))  # We learnt: @. makes transpose.()
        end
    else  # Non-orthogonal (Y)
        for i in eachindex(time_labels)
            t = time_labels[i]
            @fastmath caches.U.diag .= exp.(1im .* hamiltonian.eigvals .* t)
            
            mul!(caches.temp_op, caches.U, jump.in_eigenbasis)  # temp = U * jump
            mul!(out_matrix, caches.temp_op, caches.U', caches.prefactors[i], 1.0) # out += prefactor * U*jump*U'
        end
    end
    
    return out_matrix
end

function trotter_oft!(
    out_matrix::Matrix{ComplexF64},
    caches::OFTCaches,
    jump::JumpOp,
    energy::Float64,
    trotter::TrottTrott, 
    time_labels::Vector{Float64},
    sigma::Float64
    )

    if length(caches.prefactors) != length(time_labels)
        resize!(caches.prefactors, length(time_labels))
    end
    
    @fastmath @. caches.prefactors = exp(-time_labels^2 * sigma^2 - 1im * energy * time_labels)
    
    zero_index = findfirst(t -> t >= -1e-12, time_labels)
    
    fill!(out_matrix, 0.0)

    if jump.orthogonal
        # t = 0.0 case: U = I
        @fastmath out_matrix .+= caches.prefactors[zero_index] .* jump.in_eigenbasis

        for i in (zero_index + 1):length(time_labels)
            num_t0_steps = i - zero_index
            @fastmath caches.U.diag .= trotter.eigvals_t0 .^ num_t0_steps

            copyto!(caches.temp_op, jump.in_eigenbasis)
            # temp_op = U*jump*U', for diagonal U's:
            caches.temp_op .*= (caches.U.diag * caches.U.diag')  
            
            # Accumulate both terms in-place
            LinearAlgebra.axpby!(caches.prefactors[i], caches.temp_op, 1.0, out_matrix)
            LinearAlgebra.axpby!(conj(caches.prefactors[i]), transpose(caches.temp_op), 1.0, out_matrix)
        end
    else # Non-orthogonal jumps
        for i in eachindex(time_labels)
            num_t0_steps = i - zero_index
            @fastmath caches.U.diag .= trotter.eigvals_t0 .^ num_t0_steps
            
            # temp_op = U * jump.in_eigenbasis
            mul!(caches.temp_op, caches.U, jump.in_eigenbasis)
            
            # out_matrix += prefactor * (temp_op * U')
            mul!(out_matrix, caches.temp_op, caches.U', caches.prefactors[i], 1.0)
        end
    end
    
    return out_matrix
end


# function time_oft_integrated(energy::Float64, jump::JumpOp, hamiltonian::HamHam, beta::Float64)

#     diag_exponentiate(t) = Diagonal(exp.(1im * hamiltonian.eigvals * t))
#     integrand(t) = exp(-t^2 / beta^2 - 1im * energy * t) * diag_exponentiate(t) * jump.in_eigenbasis * diag_exponentiate(-t)
#     jump_oft = quadgk(integrand, -Inf, Inf)[1] / sqrt(2 * pi) * sqrt(sqrt(2 / pi) / beta)
#     return jump_oft
# end


#* Trotter OFT check
# energy_labels, time_labels = precompute_labels(config.domain, config)
# truncated_energy_labels = truncate_energy_labels(energy_labels, beta,
# a, b, with_linear_combination)
# time_labels = energy_labels .* (t0 / w0)
# w = -0.12
# oft_time_labels = truncate_time_labels_for_oft(time_labels, beta)
# jump = jumps[6]
# oft_trott = trotter_oft(jump, w, trotter, oft_time_labels, beta) * t0 * sqrt((sqrt(2 / pi)/beta) / (2 * pi))
# To compare trotter OFT with energy OFT: transform via trotter.eigvecs and hamiltonian.eigvecs
# oft_w = oft(jump, w, hamiltonian, beta) * sqrt(beta / sqrt(2 * pi))
# norm(oft_w - oft_trott)
