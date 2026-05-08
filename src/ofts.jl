"""
    oft!(out, eigenbasis, bohr_freqs, energy, inv_4sigma2) -> nothing

Compute the Operator Fourier Transform A(omega) in-place:
    out[i,j] = eigenbasis[i,j] * exp(-(energy - bohr_freqs[i,j])^2 * inv_4sigma2)

# Role after qf-e60

Mainline EnergyDomain paths (`apply_lindbladian!`, `_jump_contribution!`,
`_accumulate_jump_sandwich!`, `_accumulate_R_total!`, `_precompute_R!`,
`step_along_trajectory!`, `_accumulate_rho_jump!`) no longer call `oft!` —
they read the precomputed `EnergyDomainPrefactors` table built once at
`Workspace` construction time (`ws.oft_prefactors_energy`,
`_prepare_oft_prefactors_energy` in `src/energy_domain.jl`) and apply it
via `_prefactor_view` plus an element-wise broadcast. The cache stores
the closed-form `exp(-(ω - bohr_freqs)^2 / 4σ²)` envelope as
`Array{Float64, 3}`.

`oft!` is retained as the **gold-standard reference implementation** for
the per-(jump, ω) bit-equivalence harness in `test_oft_prefactors_energy.jl`.
Use the cache (`_prefactor_view(ws.oft_prefactors_energy, ω)`) in new code.
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
