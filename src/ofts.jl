using LinearAlgebra
using SparseArrays
using Random
using Printf
using Plots
using QuadGK
using Pkg

function oft(jump::JumpOp, energy::Float64, hamiltonian::HamHam, sigma::Float64)
    """Subnormalized, multiply by sqrt(1 / sigma sqrt(2 * pi))"""
    return @. jump.in_eigenbasis * exp(-(energy - hamiltonian.bohr_freqs)^2 / (4 * sigma^2))
end

function oft!(out_matrix::Matrix{ComplexF64}, jump::JumpOp, energy::Float64, hamiltonian::HamHam, sigma::Float64)
    """Subnormalized, multiply by sqrt(1 / sigma sqrt(2 * pi))"""
    @. out_matrix = jump.in_eigenbasis * exp(-(energy - hamiltonian.bohr_freqs)^2 / (4 * sigma^2)) 
    return out_matrix
end


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

# TODO: Do NUFFT for trotter too.

mutable struct NUFFTCaches{T<:AbstractVector{Float64}}
    nufft_prefactor_temp::Matrix{ComplexF64}
    gaussian_weights::Vector{ComplexF64}
    bohr_freqs_flat::T
    combined_freqs::Vector{Float64}  # omega - nu 
    out_vector::Vector{ComplexF64}
    eps::Float64
end

function NUFFTCaches(hamiltonian::HamHam, time_labels::Vector{Float64}, sigma::Float64; eps=1e-12)

    dim = size(hamiltonian.data, 1)
    nufft_prefactor_temp = Matrix{ComplexF64}(undef, dim, dim)
    gaussian_weights = ComplexF64.(exp.(-(sigma^2) .* (time_labels .^ 2)))
    bohr_freqs_flat = vec(hamiltonian.bohr_freqs)
    combined_freqs = Vector{Float64}(undef, length(bohr_freqs_flat))
    out_vector = Vector{ComplexF64}(undef, length(bohr_freqs_flat))
    return NUFFTCaches{typeof(bohr_freqs_flat)}(
        nufft_prefactor_temp, gaussian_weights, bohr_freqs_flat, combined_freqs, out_vector, eps
        )
end

function time_oft_nufft!(
    out_matrix::Matrix{ComplexF64},
    jump::JumpOp,
    energy::Float64, 
    time_labels::Vector{Float64},
    cache::NUFFTCaches
    )

    nufft_prefactor_matrix!(cache, energy, time_labels)
    @. out_matrix = jump.in_eigenbasis * cache.nufft_prefactor_temp
    return out_matrix
end

function nufft_prefactor_matrix!(
    cache::NUFFTCaches,
    energy::Float64,
    time_labels::Vector{Float64},
    )

    # Combine freqs
    @fastmath @. cache.combined_freqs = energy - cache.bohr_freqs_flat

    # Type-3 NUFFT: out_vector[k] = Σ_j gaussian_weights[j] * exp(-i * combined_freqs[k] * time_labels[j])
    FINUFFT.nufft1d3!(
        time_labels,
        cache.gaussian_weights,
        -1,
        cache.eps,
        cache.combined_freqs,
        cache.out_vector
    )

    copyto!(cache.nufft_prefactor_temp, reshape(cache.out_vector, size(cache.nufft_prefactor_temp)))
    return cache.nufft_prefactor_temp
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
# oft_trott = trotter.trafo_from_eigen_to_trotter' * oft_trott * trotter.trafo_from_eigen_to_trotter
# oft_w = oft(jump, w, hamiltonian, beta) * sqrt(beta / sqrt(2 * pi))
# norm(oft_w - oft_trott)
