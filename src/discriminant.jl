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

"""
    apply_discriminant!(out, X, lindblad_action!, sigma_quarter, sigma_inv_quarter, buffers)

Apply the KMS quantum discriminant

    D(X) = sigma^{-1/4} * L(sigma^{1/4} X sigma^{1/4}) * sigma^{-1/4}

to a `dxd` operator `X`, writing the result into `out`.  The Lindbladian
appears only through the closure `lindblad_action!(out_mat, in_mat)`, so
this routine never materialises the `d^2 x d^2` superoperator.

# Arguments
- `out::AbstractMatrix{T}`: output buffer of size `dxd`.
- `X::AbstractMatrix{T}`: input operator.
- `lindblad_action!::F`: a callable that, given `(out_mat, in_mat)`, writes
  `L(in_mat)` to `out_mat` in place.  Typical wrappers:
    - dense Lindbladian `L::Matrix` of size `d^2 x d^2`:
      `f!(o, x) = (mul!(vec(o), L, vec(x)); o)`
    - QuantumFurnace `apply_lindbladian!` (closes over a `Workspace`):
      `f!(o, x) = (apply_lindbladian!(ws, x, config, ham); copyto!(o, ws.scratch.rho_out))`
- `sigma_quarter`, `sigma_inv_quarter`: real length-`d` vectors from
  `gibbs_fractional_powers` (diagonal of `sigma^{1/4}`, `sigma^{-1/4}`).
- `buffers::DiscriminantBuffers{T}`: pre-allocated `dxd` scratch.

The routine performs zero allocations after warmup.
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

    # work1 = sigma^{1/4} X sigma^{1/4}  (fused row + column scale)
    @inbounds for j in 1:d, i in 1:d
        work1[i, j] = sigma_quarter[i] * X[i, j] * sigma_quarter[j]
    end

    # work2 = L(work1)
    lindblad_action!(work2, work1)

    # out = sigma^{-1/4} * work2 * sigma^{-1/4}
    @inbounds for j in 1:d, i in 1:d
        out[i, j] = sigma_inv_quarter[i] * work2[i, j] * sigma_inv_quarter[j]
    end

    return out
end

"""
    materialize_discriminant!(D, L, sigma_quarter, sigma_inv_quarter)
    materialize_discriminant!(D, L, gibbs::Hermitian; eps_trunc=1e-12)

In-place materialisation of the KMS quantum discriminant as a `d^2 x d^2`
matrix.  Writes into the preallocated buffer `D`.

In column-stacking vec convention,

    D[k, l] = sigma_inv_quarter[r1] * sigma_inv_quarter[c1]
             * L[k, l]
             * sigma_quarter[r2]    * sigma_quarter[c2]

where `k = (c1-1)*d + r1` and `l = (c2-1)*d + r2`.  This is the diagonal
similarity transform `D = (sigma^{-1/4} kron sigma^{-1/4}) * L *
(sigma^{1/4} kron sigma^{1/4})` written without ever forming the
`d^2 x d^2` Kronecker products explicitly -- only the two length-`d`
fractional-power vectors are needed.
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

Allocating wrapper around `materialize_discriminant!`.  Returns a freshly
allocated `d^2 x d^2` matrix of the same element type as `L`.
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

Split a square matrix into its Hermitian and anti-Hermitian parts under
the conjugate-transpose involution:

    H = (D + D') / 2   (Hermitian)
    A = (D - D') / 2   (anti-Hermitian)

`H + A = D` exactly.  For the discriminant of a KMS detailed-balance
Lindbladian, `A = 0` (the discriminant is HS-self-adjoint), so the
operator norm of `A` is the diagnostic for KMS-DB violation.
"""
function hermitian_antihermitian_split(D::AbstractMatrix)
    H = (D + D') / 2
    A = (D - D') / 2
    return (H, A)
end

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

Leading eigenvalues of the Hermitian part of the KMS quantum discriminant.

When the underlying Lindbladian is KMS detailed-balanced, the discriminant
is HS-self-adjoint (`A_part = 0`), so `H_part = D` and
`eigvals(H_part) == eigvals(L)` -- all real, with one near-zero eigenvalue
corresponding to the steady state (purification of the Gibbs state) and
the rest negative.  The "parent Hamiltonian" is `-H_part`, and `H_gap`
is its smallest positive eigenvalue (the spectral gap).

For non-KMS-DB Lindbladians, `H_part` is still Hermitian; its spectrum
no longer matches `eigvals(L)`, but it remains a useful diagnostic.

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

Compute leading eigenvalues of the Hermitian part of the discriminant of
the Lindbladian `L`.  Materialises the discriminant once and uses dense
`eigvals(Hermitian(...))`, which is `O(d^6)` and feasible up to `n=6`
(`d^2 = 4096`).
"""
function discriminant_spectrum(
    L::AbstractMatrix{<:Complex},
    gibbs::Hermitian;
    n_modes::Int = 20,
    eps_trunc::Real = 1e-12,
)
    D = materialize_discriminant(L, gibbs; eps_trunc = eps_trunc)
    H, _ = hermitian_antihermitian_split(D)

    # Numerical Hermitization (cleanup floating-point asymmetry).
    H_clean = Hermitian((H + H') / 2)
    H_eigs = eigvals(H_clean)

    # Sort by ascending |λ| (steady state -- nearest zero -- first).
    perm = sortperm(abs.(H_eigs))
    n = min(n_modes, length(H_eigs))
    leading = H_eigs[perm[1:n]]

    # Parent-Hamiltonian gap = second-smallest |λ| (first ≈ 0).
    gap = n >= 2 ? abs(leading[2]) : 0.0

    return DiscriminantSpectrum(leading, gap, n)
end

"""
    DBVerificationResult

Result of `verify_detailed_balance` -- a complete diagnostic bundle for
KMS detailed balance of a Lindbladian.

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

End-to-end KMS detailed-balance verification.  Materialises the
discriminant once, splits it into Hermitian / anti-Hermitian parts,
and computes:
- the operator-2 norm of the anti-Hermitian part (KMS-DB diagnostic),
- the discriminant norm and their ratio,
- the fixed-point residual `||D · vec(σ^{1/2})||`,
- the Hermitian-part gap (parent-Hamiltonian gap), and
- the spectral gap of `L` for cross-validation.

Uses dense `opnorm` and `eigvals`, both `O(d^6)` -- feasible up to `n=6`
(`d^2 = 4096`).  For larger systems an action-only Krylov path will be
needed; that is deferred until needed.
"""
function verify_detailed_balance(
    L::AbstractMatrix{<:Complex},
    gibbs::Hermitian;
    atol::Float64 = 1e-10,
    eps_trunc::Real = 1e-12,
)
    # Materialise the discriminant once and split.
    D = materialize_discriminant(L, gibbs; eps_trunc = eps_trunc)
    H_part, A_part = hermitian_antihermitian_split(D)

    # Operator-2 norms (Schatten-∞), per thesis Def. def:approx-db.
    A_norm   = opnorm(A_part)
    D_norm   = opnorm(D)
    rel_norm = A_norm / max(D_norm, 1e-30)

    # Fixed-point residual: D · vec(σ^{1/2}) should be 0.
    powers = gibbs_fractional_powers(gibbs; eps_trunc = eps_trunc)
    sigma_half = powers.sigma_half
    d = length(sigma_half)
    vec_sh = zeros(eltype(D), d * d)
    @inbounds for i in 1:d
        vec_sh[(i - 1) * d + i] = sigma_half[i]
    end
    fp_residual = norm(D * vec_sh)

    # Parent-Hamiltonian gap.
    H_clean = Hermitian((H_part + H_part') / 2)
    H_eigs  = eigvals(H_clean)
    H_gap   = sort(abs.(H_eigs))[2]

    # L spectral gap (smallest |Re(λ)| above zero).
    L_eigs = eigvals(Matrix(L))
    L_gap  = sort(abs.(real.(L_eigs)))[2]

    is_kms_db = rel_norm < atol

    return DBVerificationResult(
        A_norm, D_norm, rel_norm, fp_residual,
        H_gap, L_gap, is_kms_db, atol,
    )
end
