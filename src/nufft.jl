mutable struct NUFFTCaches{T<:AbstractVector{Float64}}
    nufft_prefactor_temp::Matrix{ComplexF64}
    gaussian_weights::Vector{ComplexF64}
    bohr_freqs_flat::T
    combined_freqs::Vector{Float64}  # omega - nu 
    out_vector::Vector{ComplexF64}
    eps::Float64
end
"""
Creates caches for OFT computation via NUFFT.
bohr_freqs: either the exact Bohr frequencies of the Hamiltonian, or the quasi Bohr frequencies via Trotter approximation.
"""
function NUFFTCaches(bohr_freqs::AbstractMatrix{Float64}, time_labels::Vector{Float64}, sigma::Float64; eps=1e-12)

    dim1, dim2 = size(bohr_freqs)
    @assert dim1 == dim2

    nufft_prefactor_temp = Matrix{ComplexF64}(undef, dim1, dim1)
    gaussian_weights = ComplexF64.(exp.(-(sigma^2) .* (time_labels .^ 2)))

    bohr_freqs_flat = vec(bohr_freqs)
    combined_freqs = Vector{Float64}(undef, length(bohr_freqs_flat))
    out_vector = Vector{ComplexF64}(undef, length(bohr_freqs_flat))
    return NUFFTCaches{typeof(bohr_freqs_flat)}(
        nufft_prefactor_temp, gaussian_weights, bohr_freqs_flat, combined_freqs, out_vector, eps
        )
end

"""
Precomputes NUFFT prefactor matrix (i.e. everything in the OFT that is not the jump operator). This will be then applied
entry-wise to the jump operator for fast computation of time / trotter OFTs.
"""
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
