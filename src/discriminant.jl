# KMS quantum discriminant.
# Math: $D(X) = sigma^(-1/4) L(sigma^(1/4) X sigma^(1/4)) sigma^(-1/4)$.

"""
    DiscriminantBuffers{T<:Complex}

Three reusable matrices for `apply_discriminant!`.

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

Allocate three uninitialised `dim × dim` buffers of element type `T`.
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

Allocate `ComplexF64` discriminant buffers.
"""
DiscriminantBuffers(dim::Int) = DiscriminantBuffers{ComplexF64}(dim)

"""
    gibbs_fractional_powers(gibbs::Hermitian{Complex{T}}; eps_trunc=1e-12)

Return diagonal vectors for `sigma^(1/4)`, `sigma^(-1/4)`, and `sigma^(1/2)`.

# Keywords
- `eps_trunc`: Floor applied before the negative fractional power.

The input must be diagonal in the working basis; Trotter-domain Gibbs states
generally do not meet this requirement.
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

"""
    apply_discriminant!(out, X, lindblad_action!, sigma_quarter, sigma_inv_quarter, buffers)

Apply the KMS quantum discriminant without materialising a superoperator.

# Arguments
- `out`: Destination matrix.
- `X::AbstractMatrix{T}`: input operator.
- `lindblad_action!`: In-place callable `L(out, X)`.
- `sigma_quarter`, `sigma_inv_quarter`: Diagonal fractional powers.
- `buffers`: Reusable scratch matrices.

# Returns
`out`.
"""
function apply_discriminant!(
    out::AbstractMatrix{T},
    X::AbstractMatrix{T},
    lindblad_action!::F,
    sigma_quarter::AbstractVector{<:Real},
    sigma_inv_quarter::AbstractVector{<:Real},
    buffers::DiscriminantBuffers{T},
) where {T<:Complex, F}
    work1 = buffers.work1
    work2 = buffers.work2
    d = size(X, 1)

    # Math: $work1 = sigma^(1/4) X sigma^(1/4)$.
    @inbounds for j in 1:d, i in 1:d
        work1[i, j] = sigma_quarter[i] * X[i, j] * sigma_quarter[j]
    end

    lindblad_action!(work2, work1)

    # Math: $out = sigma^(-1/4) L(work1) sigma^(-1/4)$.
    @inbounds for j in 1:d, i in 1:d
        out[i, j] = sigma_inv_quarter[i] * work2[i, j] * sigma_inv_quarter[j]
    end

    return out
end

"""
    materialize_discriminant!(D, L, sigma_quarter, sigma_inv_quarter)
    materialize_discriminant!(D, L, gibbs::Hermitian; eps_trunc=1e-12)

Materialise the KMS quantum discriminant in `D`.

Column stacking gives the diagonal similarity transform
`D = (sigma^(-1/4) tensor sigma^(-1/4)) L
(sigma^(1/4) tensor sigma^(1/4))`; Kronecker factors are not allocated.

# Returns
`D`.
"""
function materialize_discriminant!(
    D::AbstractMatrix{<:Complex},
    L::AbstractMatrix{<:Complex},
    sigma_quarter::AbstractVector{<:Real},
    sigma_inv_quarter::AbstractVector{<:Real},
)
    d = length(sigma_quarter)
    d2 = d * d

    @inbounds for l in 1:d2
        c2 = ((l - 1) ÷ d) + 1
        r2 = ((l - 1) % d) + 1
        scale_right = sigma_quarter[r2] * sigma_quarter[c2]
        for k in 1:d2
            c1 = ((k - 1) ÷ d) + 1
            r1 = ((k - 1) % d) + 1
            scale_left = sigma_inv_quarter[r1] * sigma_inv_quarter[c1]
            D[k, l] = scale_left * L[k, l] * scale_right
        end
    end
    return D
end

function materialize_discriminant!(
    D::AbstractMatrix{<:Complex},
    L::AbstractMatrix{<:Complex},
    gibbs::Hermitian;
    eps_trunc::Real = 1e-12,
)
    powers = gibbs_fractional_powers(gibbs; eps_trunc = eps_trunc)
    return materialize_discriminant!(D, L, powers.sigma_quarter, powers.sigma_inv_quarter)
end

"""
    materialize_discriminant(L, gibbs::Hermitian; eps_trunc=1e-12) -> Matrix

Allocate and return the dense KMS quantum discriminant.
"""
function materialize_discriminant(
    L::AbstractMatrix{T},
    gibbs::Hermitian;
    eps_trunc::Real = 1e-12,
) where {T<:Complex}
    D = similar(L)
    return materialize_discriminant!(D, L, gibbs; eps_trunc = eps_trunc)
end

"""
    hermitian_antihermitian_split(D::AbstractMatrix) -> (H, A)
    hermitian_antihermitian_split!(H, A, D)

Split `D` into Hermitian and anti-Hermitian parts.

Math: `\$H = (D + D^dagger)/2\$` and `\$A = (D-D^dagger)/2\$`. KMS detailed
balance implies `A = 0`.
"""
function hermitian_antihermitian_split(D::AbstractMatrix)
    H = (D + D') / 2
    A = (D - D') / 2
    return (H, A)
end

"""
    hermitian_antihermitian_split!(H, A, D) -> (H, A)

Write the Hermitian and anti-Hermitian parts of `D` into preallocated matrices.

# Arguments
- `H`, `A`: Destination matrices.
- `D`: Input square matrix.

# Returns
The tuple `(H, A)`.
"""
function hermitian_antihermitian_split!(
    H::AbstractMatrix{T},
    A::AbstractMatrix{T},
    D::AbstractMatrix{T},
) where {T<:Complex}
    @inbounds for j in axes(D, 2), i in axes(D, 1)
        d_ij = D[i, j]
        d_ji_conj = conj(D[j, i])
        H[i, j] = (d_ij + d_ji_conj) / 2
        A[i, j] = (d_ij - d_ji_conj) / 2
    end
    return (H, A)
end

"""
    DiscriminantSpectrum

Leading eigenvalues of the Hermitian discriminant part.

# Fields
- `H_eigenvalues::Vector{Float64}`: leading `n_modes` eigenvalues of the
  Hermitian part, sorted ascending by `|λ|` (steady state at index 1).
- `H_gap::Float64`: `|H_eigenvalues[2]|`, the parent-Hamiltonian gap.
- `n_modes::Int`: number of leading eigenvalues retained.
"""
struct DiscriminantSpectrum
    H_eigenvalues::Vector{Float64}
    H_gap::Float64
    n_modes::Int
end

"""
    discriminant_spectrum(L, gibbs::Hermitian; n_modes=20, eps_trunc=1e-12) -> DiscriminantSpectrum

Compute the dense Hermitian-discriminant spectrum.

# Keywords
- `n_modes`: Number of eigenvalues nearest zero to retain.
- `eps_trunc`: Gibbs fractional-power floor.

# Returns
A `DiscriminantSpectrum`. Dense eigendecomposition costs `O(d^6)`.
"""
function discriminant_spectrum(
    L::AbstractMatrix{<:Complex},
    gibbs::Hermitian;
    n_modes::Int = 20,
    eps_trunc::Real = 1e-12,
)
    D = materialize_discriminant(L, gibbs; eps_trunc = eps_trunc)
    H, _ = hermitian_antihermitian_split(D)

    # Remove roundoff asymmetry before the Hermitian eigensolve.
    H_clean = Hermitian((H + H') / 2)
    H_eigs = eigvals(H_clean)

    # Sort by distance to the stationary eigenvalue.
    perm = sortperm(abs.(H_eigs))
    n = min(n_modes, length(H_eigs))
    leading = H_eigs[perm[1:n]]

    gap = n >= 2 ? abs(leading[2]) : 0.0

    return DiscriminantSpectrum(leading, gap, n)
end

"""
    DBVerificationResult

Diagnostics returned by `verify_detailed_balance`.

# Fields
- `antihermitian_norm::Float64`: `||A||_{2→2}`, the operator-2 norm of the
  anti-Hermitian part of the discriminant.  KMS-DB iff this is zero.
- `discriminant_norm::Float64`: `||D||_{2→2}`, for normalisation.
- `relative_norm::Float64`: `antihermitian_norm / discriminant_norm`.
- `fixed_point_residual::Float64`: `||D · vec(σ^{1/2})||_2`.  Should be
  ≈ 0 because `σ^{1/2}` is the zero eigenvector of D for any Lindbladian
  with σ as its steady state (independent of detailed balance).
- `hermitian_part_gap::Float64`: smallest |λ| of the Hermitian part above
  the steady-state zero -- the parent-Hamiltonian gap.
- `spectral_gap_L::Float64`: `min |Re(λ)|` over nonzero eigenvalues of L.
  For a KMS-DB Lindbladian this equals `hermitian_part_gap` (similarity
  transform preserves spectrum).
- `is_kms_db::Bool`: `relative_norm < atol`.
- `atol::Float64`: the threshold used for `is_kms_db`.
"""
struct DBVerificationResult
    antihermitian_norm::Float64
    discriminant_norm::Float64
    relative_norm::Float64
    fixed_point_residual::Float64
    hermitian_part_gap::Float64
    spectral_gap_L::Float64
    is_kms_db::Bool
    atol::Float64
end

"""
    verify_detailed_balance(L, gibbs::Hermitian; atol=1e-10, eps_trunc=1e-12) -> DBVerificationResult

Verify KMS detailed balance from the dense quantum discriminant.

# Keywords
- `atol`: Threshold for the relative anti-Hermitian norm.
- `eps_trunc`: Gibbs fractional-power floor.

# Returns
A `DBVerificationResult` containing the defect, fixed-point residual, and
dense gap cross-check. The calculation costs `O(d^6)`.
"""
function verify_detailed_balance(
    L::AbstractMatrix{<:Complex},
    gibbs::Hermitian;
    atol::Float64 = 1e-10,
    eps_trunc::Real = 1e-12,
)
    D = materialize_discriminant(L, gibbs; eps_trunc = eps_trunc)
    H_part, A_part = hermitian_antihermitian_split(D)

    A_norm   = opnorm(A_part)
    D_norm   = opnorm(D)
    rel_norm = A_norm / max(D_norm, 1e-30)

    # Math: the stationary discriminant vector is $vec(sigma^(1/2))$.
    powers = gibbs_fractional_powers(gibbs; eps_trunc = eps_trunc)
    sigma_half = powers.sigma_half
    d = length(sigma_half)
    vec_sh = zeros(eltype(D), d * d)
    @inbounds for i in 1:d
        vec_sh[(i - 1) * d + i] = sigma_half[i]
    end
    fp_residual = norm(D * vec_sh)

    H_clean = Hermitian((H_part + H_part') / 2)
    H_eigs  = eigvals(H_clean)
    H_gap   = sort(abs.(H_eigs))[2]

    L_eigs = eigvals(Matrix(L))
    L_gap  = sort(abs.(real.(L_eigs)))[2]

    is_kms_db = rel_norm < atol

    return DBVerificationResult(
        A_norm, D_norm, rel_norm, fp_residual,
        H_gap, L_gap, is_kms_db, atol,
    )
end
