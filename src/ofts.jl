"""
    oft!(out, eigenbasis, bohr_freqs, energy, inv_4sigma2) -> nothing

Compute the Operator Fourier Transform A(omega) in-place:
    out[i,j] = eigenbasis[i,j] * exp(-(energy - bohr_freqs[i,j])^2 * inv_4sigma2)

This is the optimized concrete-typed signature for all simulation paths.
EnergyDomain only (Time/TrotterDomain uses NUFFT prefactors).
"""
@inline function oft!(
    out::Matrix{T},
    eigenbasis::Matrix{T},
    bohr_freqs::Matrix{<:Real},
    energy::Real,
    inv_4sigma2::Real,
) where {T<:Complex}
    @. out = eigenbasis * exp(-(energy - bohr_freqs)^2 * inv_4sigma2)
    return nothing
end

# Note: the legacy direct-summation Time/Trotter OFT routines (`time_oft!`,
# `trotter_oft!`) and their `OFTCaches` struct have been retired. The mainline
# code paths use NUFFT prefactors (`_prepare_oft_nufft_prefactors`) instead.
# The original implementations are preserved for personal reference at
# `src/staging/ofts.jl` (commented-out, not included in the build).
