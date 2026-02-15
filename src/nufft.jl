struct NUFFTPrefactors{T<:AbstractFloat, A<:AbstractArray{Complex{T}, 3}}
    data::A
    energy_labels::Vector{T}
    energy_to_index::Dict{T, Int}
end

"""Exact deduplication (no approximation) plus inverse map to reconstruct the full array."""
function _unique_with_invmap(v::AbstractVector{<:AbstractFloat})
    uniq = unique(v)
    T = eltype(v)
    idx = Dict{T,Int}(u => i for (i,u) in enumerate(uniq))
    invmap = Vector{Int32}(undef, length(v))
    @inbounds for k in eachindex(v)
        invmap[k] = Int32(idx[v[k]])
    end
    return uniq, invmap
end

function prepare_oft_nufft_prefactors(
    bohr_freqs::AbstractMatrix{<:AbstractFloat},
    time_labels::Vector{<:AbstractFloat},
    energy_labels::Vector{<:AbstractFloat},
    sigma::AbstractFloat;
    eps::Float64 = 1e-12,
    nthreads::Int = 1,
    use_shared_array::Bool = (nprocs() > 1),
)
    dim1, dim2 = size(bohr_freqs)
    @assert dim1 == dim2
    dim = dim1

    # Determine the element type T from energy_labels (the stored output type).
    T = eltype(energy_labels)

    # Promote inputs to Float64 for FINUFFT computation (FINUFFT requires Float64).
    bohr_freqs_f64 = Float64.(bohr_freqs)
    time_labels_f64 = Float64.(time_labels)
    energy_labels_f64 = Float64.(energy_labels)
    sigma_f64 = Float64(sigma)

    # Flatten Bohr frequencies but only retain unique ones.
    bohr_flat = vec(bohr_freqs_f64)
    unique_bohr_flat, invmap = _unique_with_invmap(bohr_flat)

    base_weights = ComplexF64.(exp.(-(sigma_f64^2) .* (time_labels_f64 .^ 2)))
    input_weights = Matrix{ComplexF64}(undef, length(time_labels_f64), 1)
    out_nufft = Matrix{ComplexF64}(undef, length(unique_bohr_flat), 1)

    # Allocate prefactor stack in Complex{T}.
    CT = Complex{T}
    prefactors = if use_shared_array && (nprocs() > 1)
        # Use all processes by default (SharedArray requires co-located processes).
        SharedArray{CT}(dim, dim, length(energy_labels))
    else
        Array{CT}(undef, dim, dim, length(energy_labels))
    end

    # FINUFFT plan - type 3, 1D, + sign: exp(+i s x)
    plan = FINUFFT.finufft_makeplan(3, 1, +1, 1, eps; dtype=Float64, nthreads=nthreads)
    # For dim=1: xj are sources; s are targets; yj,zj,t,u unused.
    empty = Float64[]
    FINUFFT.finufft_setpts!(plan,
        time_labels_f64,    # xj
        empty,              # yj (unused)
        empty,              # zj (unused)
        unique_bohr_flat,   # s (targets)
        empty,              # t (unused)
        empty               # u (unused)
    )

    # Fill prefactors[:,:,k] for each w=energy_labels[k].
    @inbounds for (k, omega) in enumerate(energy_labels_f64)
        # input_weights[j] = exp(-sigma^2 t_j^2) * exp(-i w t_j)
        @fastmath @. input_weights[:, 1] = base_weights * cis(-omega * time_labels_f64)

        # out_nufft[u] = sum_j input_weights[j] * exp(+i unique_targets[u] * time_labels[j])
        FINUFFT.finufft_exec!(plan, input_weights, out_nufft)

        # Scatter back to the full (dim x dim) ordering, converting to Complex{T}.
        @views full_bohr_prefac_omega = prefactors[:, :, k]
        @inbounds for p in eachindex(invmap)
            full_bohr_prefac_omega[p] = CT(out_nufft[Int(invmap[p]), 1])
        end
    end

    FINUFFT.finufft_destroy!(plan)

    energy_to_index = Dict{T,Int}(omega => i for (i, omega) in enumerate(energy_labels))
    return NUFFTPrefactors(prefactors, energy_labels, energy_to_index)
end

"""Convenience view: returns prefactor matrix for energy w without allocating."""
@inline function prefactor_view(nufft_prefactors::NUFFTPrefactors, omega)
    k = nufft_prefactors.energy_to_index[omega]
    return @view nufft_prefactors.data[:, :, k]
end
