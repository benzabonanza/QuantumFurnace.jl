# ============================================================================
# Quantum discriminant of a KMS detailed-balance Lindbladian
# ============================================================================
#
# Computes the discriminant
#
#     D(sigma, L) := sigma^{-1/4} L( sigma^{1/4} . sigma^{1/4} ) sigma^{-1/4}
#
# (thesis 1_preliminaries.tex Eq. eq:discriminant; identical to Chen et al.
# 2025 "Purifying Lindbladians" Sec. I).  KMS detailed balance is equivalent
# to D being Hermitian under the Hilbert-Schmidt inner product.
#
# Phase 1 (this file): foundational primitives
#   - DiscriminantBuffers{T}        : preallocated dxd scratch matrices
#   - gibbs_fractional_powers       : sigma^{1/4}, sigma^{-1/4}, sigma^{1/2}
#                                     extracted from a diagonal Gibbs state
#                                     (Bohr / Energy / Time domains).

"""
    DiscriminantBuffers{T<:Complex}

Preallocated `dxd` scratch matrices for the closure-style discriminant
action `apply_discriminant!` (Phase 2).  Three buffers are needed: one to
hold `sigma^{1/4} X sigma^{1/4}`, one to receive the Lindbladian action,
and one for the final left/right multiplication by `sigma^{-1/4}`.

# Fields
- `work1`, `work2`, `work3`: `Matrix{T}` of size `dim x dim`.
"""
struct DiscriminantBuffers{T<:Complex}
    work1::Matrix{T}
    work2::Matrix{T}
    work3::Matrix{T}
end

"""
    DiscriminantBuffers{T}(dim::Int) where {T<:Complex}

Construct empty `DiscriminantBuffers` with three `Matrix{T}(undef, dim, dim)`
scratch buffers.
"""
function DiscriminantBuffers{T}(dim::Int) where {T<:Complex}
    return DiscriminantBuffers{T}(
        Matrix{T}(undef, dim, dim),
        Matrix{T}(undef, dim, dim),
        Matrix{T}(undef, dim, dim),
    )
end

"""
    DiscriminantBuffers(dim::Int)

Convenience constructor with `T = ComplexF64`.
"""
DiscriminantBuffers(dim::Int) = DiscriminantBuffers{ComplexF64}(dim)

"""
    gibbs_fractional_powers(gibbs::Hermitian{Complex{T}}; eps_trunc=1e-12)

Extract the diagonal of `gibbs` (assumed diagonal in the working basis --
true for `BohrDomain`, `EnergyDomain`, `TimeDomain`) and return three
length-`d` vectors

    sigma_quarter      = sigma^{1/4}
    sigma_inv_quarter  = sigma^{-1/4}
    sigma_half         = sigma^{1/2}

as a `NamedTuple{(:sigma_quarter, :sigma_inv_quarter, :sigma_half)}`.

Diagonal entries below `eps_trunc` are floored before exponentiation to
prevent `sigma^{-1/4}` from blowing up at very low temperature.  In the
QuantumFurnace.jl regime (rescaled `H in [0, 0.45]`, `beta <= 20`) the
floor is never triggered in practice; it is kept as a defensive guard.

`TrotterDomain` is not yet supported -- there `gibbs` is generally
non-diagonal in the working basis and an eigendecomposition would be
required.
"""
function gibbs_fractional_powers(
    gibbs::Hermitian{Complex{T}, Matrix{Complex{T}}};
    eps_trunc::Real = 1e-12,
) where {T<:AbstractFloat}
    diag_real = real.(diag(gibbs))
    diag_safe = max.(diag_real, T(eps_trunc))
    return (
        sigma_quarter     = diag_safe .^ T(0.25),
        sigma_inv_quarter = diag_safe .^ T(-0.25),
        sigma_half        = diag_safe .^ T(0.5),
    )
end
