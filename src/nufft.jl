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

function _prepare_oft_nufft_prefactors(
    bohr_freqs::AbstractMatrix{<:AbstractFloat},
    time_labels::Vector{<:AbstractFloat},
    energy_labels::Vector{<:AbstractFloat},
    filter::AbstractFilter;
    eps::Float64 = 1e-12,
    nthreads::Int = 1,
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

    # Flatten Bohr frequencies but only retain unique ones.
    bohr_flat = vec(bohr_freqs_f64)
    unique_bohr_flat, invmap = _unique_with_invmap(bohr_flat)

    # Filter time-domain weights. Gaussian: real exp(-σ² t²); DLL: complex
    # closed form. Both promote to ComplexF64 for FINUFFT.
    base_weights = ComplexF64.(time_kernel.(Ref(filter), time_labels_f64))
    input_weights = Matrix{ComplexF64}(undef, length(time_labels_f64), 1)
    out_nufft = Matrix{ComplexF64}(undef, length(unique_bohr_flat), 1)

    # Allocate prefactor stack in Complex{T}.
    CT = Complex{T}
    prefactors = Array{CT}(undef, dim, dim, length(energy_labels))

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
@inline function _prefactor_view(nufft_prefactors::NUFFTPrefactors, omega)
    k = nufft_prefactors.energy_to_index[omega]
    return @view nufft_prefactors.data[:, :, k]
end

# qf-e60.1: EnergyDomain Gaussian prefactor cache. Real-valued analogue of
# NUFFTPrefactors. The TimeDomain/TrotterDomain `oft_nufft_prefactors` carries
# a complex time-domain phase from FINUFFT; the EnergyDomain Gaussian
# `exp(-(ω - bohr)^2 / 4σ^2)` has no phase, so storage is `Float64` instead
# of `ComplexF64` — half the memory of the NUFFT cache at the same `(d, d, N_w)`
# shape.
struct EnergyDomainPrefactors{T<:Real, A<:AbstractArray{T, 3}}
    data::A
    energy_labels::Vector{T}
    energy_to_index::Dict{T, Int}
end

"""
    _prepare_oft_prefactors_energy(bohr_freqs, energy_labels, sigma)
        -> EnergyDomainPrefactors

Build the closed-form EnergyDomain Gaussian OFT prefactor table:

    data[i, j, k] = exp(-(bohr_freqs[i, j] - energy_labels[k])^2 / (4σ^2))

Bit-equivalent to `oft!(out, eigenbasis=I, bohr_freqs, energy_labels[k], 1/(4σ^2))`
at FP64 precision (same `exp` instruction, same rounding). Mainline EnergyDomain
matvec / construction paths consume this table via `_prefactor_view`; `oft!`
remains exported as the gold-standard cross-check oracle.
"""
function _prepare_oft_prefactors_energy(
    bohr_freqs::AbstractMatrix{<:Real},
    energy_labels::Vector{T},
    sigma::Real,
) where {T<:AbstractFloat}
    inv_4sigma2 = T(1) / (4 * T(sigma)^2)
    d1, d2 = size(bohr_freqs)
    @assert d1 == d2 "bohr_freqs must be square"
    d = d1
    N_w = length(energy_labels)
    bf = T.(bohr_freqs)

    data = Array{T, 3}(undef, d, d, N_w)
    @inbounds for k in 1:N_w
        w = energy_labels[k]
        @inbounds for j in 1:d, i in 1:d
            Δ = bf[i, j] - w
            data[i, j, k] = exp(-Δ^2 * inv_4sigma2)
        end
    end

    energy_to_index = Dict{T, Int}(omega => i for (i, omega) in enumerate(energy_labels))
    return EnergyDomainPrefactors(data, energy_labels, energy_to_index)
end

"""Convenience view: returns the real-valued Gaussian prefactor matrix for ω."""
@inline function _prefactor_view(p::EnergyDomainPrefactors, omega)
    k = p.energy_to_index[omega]
    return @view p.data[:, :, k]
end
