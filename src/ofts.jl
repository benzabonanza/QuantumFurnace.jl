"""
    oft!(out, eigenbasis, bohr_freqs, energy, inv_4sigma2) -> nothing

Compute one energy-domain operator Fourier component in place.

# Arguments
- `out`: Destination matrix.
- `eigenbasis`: Jump operator in the Hamiltonian eigenbasis.
- `bohr_freqs`: Matrix of Bohr frequencies.
- `energy`: Target energy label.
- `inv_4sigma2`: Gaussian exponent factor `1 / (4 sigma^2)`.

# Returns
`nothing`.
"""
@inline function oft!(
    out::Matrix{T},
    eigenbasis::Matrix{T},
    bohr_freqs::Matrix{<:Real},
    energy::Real,
    inv_4sigma2::Real,
) where {T<:Complex}
    # Math: $A(omega)_(i j) = A_(i j) exp(-(omega - Delta_(i j))^2 / (4 sigma^2))$.
    @. out = eigenbasis * exp(-(energy - bohr_freqs)^2 * inv_4sigma2)
    return nothing
end
