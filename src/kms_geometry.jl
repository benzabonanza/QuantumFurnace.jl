# ============================================================================
# KMS-geometry utilities for Lindbladian comparison (epic qf-mto.{1,2,3})
# ============================================================================
#
# Status (2026-05-04): diagnostic only. None of the helpers below are called
# from the simulator hot path. They were used to confirm that DLL and CKG
# Lindbladians have comparable natural normalisations, so mainline mixing-time
# comparisons can be run without rescaling either generator. Kept available
# (and tested in `test/test_kms_geometry.jl`) so the qf-mto sweep results
# remain reproducible; new comparison work that doesn't need the KMS-geometry
# decomposition can ignore this module.
#
# Compute the KMS-Poincaré gap λ(L), maximum KMS Dirichlet rate Λ_max(L), and
# the *intrinsic* (scale-invariant) ratio ρ = λ / Λ_max for any (Lindbladian,
# Gibbs state) pair, plus a 1→1 norm upper bound surrogate for the dissipator.
# These are the numbers needed to compare CKG smooth-Metro vs DLL Metro vs DLL
# Gauss on a fair footing (see `.claude-memory/fair_comparison_dirichlet_qf_mto.md`).
#
# Mathematical setup.
#
#   Operators X, Y on H_d carry the (β, σ)-KMS inner product
#       ⟨X, Y⟩_{KMS, σ} = tr(X† σ^{1/2} Y σ^{1/2}).
#   KMS detailed balance for the Schrödinger generator L_S is exactly the
#   HS-self-adjointness of the *quantum discriminant*
#       D_S(X) = σ^{-1/4} L_S(σ^{1/4} X σ^{1/4}) σ^{-1/4}.
#   In the column-stacking vec convention this is the matrix
#       D_S = (σ^{-1/4} ⊗ σ^{-1/4}) · L_S · (σ^{1/4} ⊗ σ^{1/4})
#           = Φ⁻¹ L_S Φ,                    Φ := kron(σ^{1/4}, σ^{1/4}).
#   See `src/discriminant.jl::materialize_discriminant!` and
#   `verify_detailed_balance` for the same convention used by the production code.
#
#   The KMS Dirichlet form is
#       E_L(X) = -Re ⟨X, L^H(X)⟩_{KMS, σ}
#              = -Re tr(X† σ^{1/2} L^H(X) σ^{1/2}).
#   For a KMS-DBC L: E_L(X) ≥ 0 for all Hermitian X, with equality on X ∝ I.
#   The KMS-Poincaré gap and maximum Dirichlet rate are
#       λ(L)   = inf  E_L(X) / Var_{KMS}(X)    over X ⊥ I
#       Λ_max  = sup  E_L(X) / Var_{KMS}(X)    over X ⊥ I.
#   Both equal the smallest / largest absolute eigenvalues of -D_S restricted
#   to the subspace orthogonal to the constant direction.  D_S has σ^{1/2} as
#   its zero eigenvector — see `verify_detailed_balance`'s fixed-point check.
#
#   Scale-invariant ratio:  ρ_intrinsic(L) := λ(L) / Λ_max(L) ∈ (0, 1].
#
# Coherent term: the Lamb-shift -i[B_total, ρ] in QuantumFurnace's KMS
# construction is anti-Hermitian on operator space and contributes an
# *HS-anti-Hermitian* piece to D_S.  Symmetrising D_S → (D_S + D_S†)/2
# annihilates that piece, so λ and Λ_max do not change when the coherent
# contribution is restored.
#
# For the dissipator 1→1 norm we use the upper bound 4·sup_a ‖L_a† L_a‖_∞
# (Wolf & Pérez-García style); for DLL the L_a are obtained directly via
# `dll_lindblad_op_bohr`, while CKG would require an SDP eigendecomposition
# of α (deferred).
#
# Precondition: `materialize_discriminant` requires a *diagonal* Gibbs state in
# the working basis (true for BohrDomain / EnergyDomain / TimeDomain). The
# TrotterDomain Gibbs is generally non-diagonal; consult the discriminant module.

# ---------------------------------------------------------------------------
# Section 1 — KMS inner-product algebra
# ---------------------------------------------------------------------------

"""
    kms_inner_product(X, Y, sigma; sigma_sqrt=nothing) -> Complex

KMS inner product `⟨X, Y⟩_{KMS, σ} = tr(X† σ^{1/2} Y σ^{1/2})`.

Pass an explicit `sigma_sqrt = sqrt(Hermitian(σ))` to amortise the matrix
square root over many calls.
"""
function kms_inner_product(
    X::AbstractMatrix,
    Y::AbstractMatrix,
    sigma::AbstractMatrix;
    sigma_sqrt::Union{Nothing, AbstractMatrix} = nothing,
)
    s12 = sigma_sqrt === nothing ? sqrt(Hermitian(sigma)) : sigma_sqrt
    return tr(X' * s12 * Y * s12)
end

"""
    kms_norm(X, sigma; sigma_sqrt=nothing) -> Real

KMS norm `√⟨X, X⟩_{KMS}`.  Always returns a real number ≥ 0.
"""
function kms_norm(
    X::AbstractMatrix,
    sigma::AbstractMatrix;
    sigma_sqrt::Union{Nothing, AbstractMatrix} = nothing,
)
    val = kms_inner_product(X, X, sigma; sigma_sqrt = sigma_sqrt)
    return sqrt(max(real(val), 0.0))
end

"""
    kms_variance(X, sigma; sigma_sqrt=nothing) -> Real

KMS variance `Var_{KMS}(X) = ‖X − ⟨I, X⟩_{KMS} · I‖_{KMS}²`.

Note `⟨I, X⟩_{KMS} = tr(σ^{1/2} X σ^{1/2}) = tr(σ X)` is the Gibbs expectation
of `X`.  Variance vanishes iff `X` is a multiple of the identity.
"""
function kms_variance(
    X::AbstractMatrix,
    sigma::AbstractMatrix;
    sigma_sqrt::Union{Nothing, AbstractMatrix} = nothing,
)
    s12 = sigma_sqrt === nothing ? sqrt(Hermitian(sigma)) : sigma_sqrt
    mean_X = real(tr(sigma * X))   # ⟨I, X⟩_{KMS} = tr(σ X) for Hermitian X
    Y = X .- mean_X .* I(size(X, 1))
    return real(tr(Y' * s12 * Y * s12))
end


# ---------------------------------------------------------------------------
# Section 2 — Dense superoperator + KMS quantum discriminant (reuse production)
# ---------------------------------------------------------------------------

"""
    build_dense_superoperator(L_apply!, d; T=ComplexF64) -> Matrix{T}

Materialise a `d² × d²` dense superoperator from a column-major matvec
`L_apply!(out_buf, X)` that maps a `d×d` matrix to a `d×d` matrix.  Each
column is built via a unit basis matrix and a copy.

Use as
```julia
L_super = build_dense_superoperator((out, X) -> begin
    apply_lindbladian!(ws, X, cfg, ham); copyto!(out, ws.scratch.rho_out); out
end, d)
```
"""
function build_dense_superoperator(L_apply!::F, d::Integer; T = ComplexF64) where {F}
    out = zeros(T, d * d, d * d)
    in_buf  = zeros(T, d, d)
    out_buf = zeros(T, d, d)
    for col in 1:(d * d)
        fill!(in_buf, 0)
        in_buf[col] = one(T)              # column-major linear index
        L_apply!(out_buf, in_buf)
        @views out[:, col] .= vec(out_buf)
    end
    return out
end

"""
    _kms_discriminant_with_sym(L_super, sigma) -> (D_super, D_sym)

Build the quantum discriminant `D_S = Φ⁻¹ L_S Φ` and its Hermitian part
`D_sym = (D_S + D_S†)/2` via the production
`materialize_discriminant(L_super, Hermitian(σ))` from `src/discriminant.jl`
(no explicit Kronecker products).  For a KMS-DBC Schrödinger generator L_S,
D_S itself is HS-self-adjoint and `D_sym ≈ D_S` (anti-Hermitian residual at
machine precision).  When the Lamb-shift coherent piece is included, D_S has
a non-trivial anti-Hermitian part whose contribution to the Dirichlet form
is identically zero, so `D_sym` carries all spectral data relevant to the
gap and Λ_max.

Convention matches `src/discriminant.jl::materialize_discriminant!`.
"""
function _kms_discriminant_with_sym(
    L_super::AbstractMatrix{<:Complex},
    sigma::AbstractMatrix,
)
    D_super = materialize_discriminant(L_super, Hermitian(Matrix(sigma)))
    D_sym   = Hermitian((D_super .+ D_super') ./ 2)
    return D_super, D_sym
end


# ---------------------------------------------------------------------------
# Section 3 — KMS Dirichlet form
# ---------------------------------------------------------------------------

"""
    kms_dirichlet_form(L_super_H, X, sigma; sigma_sqrt=nothing) -> Real

KMS Dirichlet form `E_L(X) = -Re ⟨X, L^H(X)⟩_{KMS, σ}` where `L^H` is the
**Heisenberg-picture** generator on observables (HS-adjoint of the Schrödinger
generator).  For KMS-DBC L:
  • `E_L(X) ≥ 0` for all Hermitian X,
  • `E_L(X) = 0` for `X = c·I` (since `L^H(I) = 0` — every Lindbladian
    preserves trace, so the Heisenberg generator is unital).

In the column-stacking vec convention, if `L_super_S` is the Schrödinger
superoperator built by `build_dense_superoperator(...)`, then
`L_super_H = L_super_S'` (conjugate-transpose) is its HS-adjoint — pass the
Heisenberg superoperator explicitly to keep the algebra unambiguous.
"""
function kms_dirichlet_form(
    L_super_H::AbstractMatrix,
    X::AbstractMatrix,
    sigma::AbstractMatrix;
    sigma_sqrt::Union{Nothing, AbstractMatrix} = nothing,
)
    d = size(X, 1)
    LX_vec = L_super_H * vec(X)
    LX     = reshape(LX_vec, d, d)
    return -real(kms_inner_product(X, LX, sigma; sigma_sqrt = sigma_sqrt))
end


# ---------------------------------------------------------------------------
# Section 4 — Spectral gap λ(L) and Section 5 — Λ_max + intrinsic ratio
# ---------------------------------------------------------------------------

"""
    _kms_dirichlet_eigvals(L_super, sigma; project_constants=true) -> Vector{Float64}

Return the *non-spurious* real eigenvalues of `-D_sym`, where D_sym is the
Hermitian part of the quantum discriminant `D_S = Φ⁻¹ L_S Φ`.  For a
KMS-DBC L_S these are the Dirichlet rates; the smallest is the
KMS-Poincaré gap, the largest is Λ_max.

`project_constants=true` projects out the constant direction (eigenvector of
L_S corresponding to the steady state, which under Φ⁻¹ becomes σ^{1/2}) and
drops the spurious zero that the projection introduces.
"""
function _kms_dirichlet_eigvals(
    L_super::AbstractMatrix{<:Complex},
    sigma::AbstractMatrix;
    project_constants::Bool = true,
)
    _, D_sym = _kms_discriminant_with_sym(L_super, sigma)
    M = -Matrix(D_sym)                           # eigenvalues = Dirichlet rates ≥ 0
    if project_constants
        # The discriminant zero direction is vec(σ^{1/2}) (verified by
        # `verify_detailed_balance`'s fp_residual check).  Normalise: since
        # ‖σ^{1/2}‖_F^2 = tr(σ) = 1, vec(σ^{1/2}) is already unit-norm.
        s12 = sqrt(Hermitian(Matrix(sigma)))
        v = vec(s12) ./ norm(vec(s12))
        # Symmetric projector P = I - v v†; M' := P M P
        Mp = M .- v * (v' * M) .- (M * v) * v' .+ v * (v' * (M * v)) * v'
        Mp = Hermitian((Mp .+ Mp') ./ 2)
        ev = eigvals(Mp)
        ev_sorted = sort(real.(ev); by = abs)
        # First eigenvalue is the spurious zero introduced by the projector.
        return ev_sorted[2:end]
    else
        return sort(real.(eigvals(Hermitian((M .+ M') ./ 2))); by = abs)
    end
end

"""
    spectral_gap_kms(L_super, sigma; project_constants=true) -> NamedTuple

KMS-Poincaré gap

    λ(L) = inf_{X ⊥ I, X ≠ 0}  -⟨X, L(X)⟩^KMS_σ / Var^KMS_σ(X)

= smallest non-zero eigenvalue of the symmetrised KMS discriminant.
Controls the τ_mix UPPER BOUND via the quantum Poincaré inequality
`τ_mix(ε) ≤ (1/λ) · log(‖σ^{-½}‖ / ε)` (Kochanowski–Alhambra–Capel–Rouzé
2024 Eq. (1)).

Returns `(; gap, eigvals)` where `eigvals` are the d²−1 Dirichlet rates
sorted by absolute value (smallest first); `gap = abs(eigvals[1])`.

Dense path (O(d⁶)).  For matrix-free large-n problems use
`krylov_spectral_gap` (which produces the same gap on KMS-DBC Lindbladians).
"""
function spectral_gap_kms(
    L_super::AbstractMatrix{<:Complex},
    sigma::AbstractMatrix;
    project_constants::Bool = true,
)
    ev = _kms_dirichlet_eigvals(L_super, sigma; project_constants = project_constants)
    return (; gap = abs(ev[1]), eigvals = ev)
end

"""
    max_dirichlet_rate_kms(L_super, sigma; project_constants=true) -> Real

Largest Dirichlet rate `Λ_max(L)` — the fastest-decaying mode of the KMS
Dirichlet form.  Together with the gap `λ(L)` it characterises the
spectral spread of the dissipative part.

Note: `Λ_max` does NOT directly bound τ_mix from below in any standard
Poincaré-style sense — see `intrinsic_mixing_ratio` for the role this
quantity plays as the "max" half of a conditioning ratio.
"""
function max_dirichlet_rate_kms(
    L_super::AbstractMatrix{<:Complex},
    sigma::AbstractMatrix;
    project_constants::Bool = true,
)
    ev = _kms_dirichlet_eigvals(L_super, sigma; project_constants = project_constants)
    return maximum(abs.(ev))
end

"""
    intrinsic_mixing_ratio(L_super, sigma) -> Real

Scale-invariant **conditioning ratio** of the KMS Dirichlet form,

    ρ(L) = λ(L) / Λ_max(L)  ∈ (0, 1]

where `λ` is the spectral gap (`spectral_gap_kms`) and `Λ_max` is the
largest Dirichlet rate (`max_dirichlet_rate_kms`).  Both numerator and
denominator scale linearly under `L → c·L`, so ρ is invariant under any
uniform time rescaling — exposing the "geometric efficiency" of the
sampler independent of its absolute generator strength.

NOTE this is NOT a τ_mix bound (neither upper nor lower).  It is a
conditioning number / functional-inequality ratio analogous to the
Poincaré-vs-Sobolev gap ratio used in classical Markov chain mixing
analysis [Kastoryano–Temme 2013, *J. Math. Phys.* 54, 052202;
Diaconis–Saloff-Coste 1996].  Two samplers with the same ρ have
equivalent rate-per-unit-norm Dirichlet structure; the absolute mixing
time still depends on the generator's overall scale.
"""
function intrinsic_mixing_ratio(
    L_super::AbstractMatrix{<:Complex},
    sigma::AbstractMatrix,
)
    ev = _kms_dirichlet_eigvals(L_super, sigma)
    return abs(ev[1]) / maximum(abs.(ev))
end


# ---------------------------------------------------------------------------
# Section 6 — 1→1 norm bound, HS-induced norm, and Tr(α) diagnostic
# ---------------------------------------------------------------------------

"""
    hs_operator_norm(L_super) -> Real

Hilbert–Schmidt-induced operator norm `‖L‖_{2→2}` = largest singular value of
the dense superoperator. Cheap to compute; bounded relative to the trace-norm-
induced 1→1 norm by

    ‖L‖_{2→2} ≤ ‖L‖_{1→1} ≤ d · ‖L‖_{2→2}

(Watrous, *The Theory of Quantum Information* §3.3.2). I.e., `hs_operator_norm`
is a *lower* bound on `‖L‖_{1→1}`, loose by up to a factor of d.
"""
hs_operator_norm(L_super::AbstractMatrix) = opnorm(L_super)

"""
    hs_operator_norm_krylov(L_apply!, L_apply_adj!, d; ...) -> Real

Matrix-free Hilbert–Schmidt-induced operator norm via Golub–Kahan–Lanczos
bidiagonalization (`KrylovKit.svdsolve`).  Sibling of `hs_operator_norm` for
problems where building the dense `d² × d²` superoperator is infeasible
(memory: `d⁴ · 16 B`; at d=128 / n=7 already 4 GB).

Inputs are two closures that act on `d × d` matrices in place:
- `L_apply!(out::Matrix, X::Matrix)` writes the forward action `L(X)` into `out`.
- `L_apply_adj!(out::Matrix, X::Matrix)` writes the HS-adjoint action `L*(X)`.

Returns the largest singular value of the implicit superoperator (the HS-induced
operator norm).  Internally uses a single ComplexF64 in/out buffer pair shared
between forward and adjoint matvecs since `svdsolve` calls them sequentially
within the GKL iteration.

# Keyword arguments
- `T::Type=ComplexF64`: working scalar type
- `tol::Real=1e-12`: GKL convergence tolerance
- `krylovdim::Int=30`: maximum Krylov subspace dimension
- `maxiter::Int=100`: maximum GKL restarts
- `max_retries::Int=3`: retries with `krylovdim *= 1.5` on partial convergence

# Example
```julia
ws_ref  = Workspace(cfg_ref, ham, jumps)
ws_test = Workspace(cfg_test, ham, jumps)
buf = zeros(ComplexF64, d, d)
diff!(out, X) = begin
    apply_lindbladian!(ws_ref,  X, cfg_ref,  ham); copyto!(buf, ws_ref.scratch.rho_out)
    apply_lindbladian!(ws_test, X, cfg_test, ham); axpby!(-1, ws_test.scratch.rho_out, 1, buf)
    copyto!(out, buf)
end
adj!(out, X) = begin
    apply_adjoint_lindbladian!(ws_ref,  X, cfg_ref,  ham); copyto!(buf, ws_ref.scratch.rho_out)
    apply_adjoint_lindbladian!(ws_test, X, cfg_test, ham); axpby!(-1, ws_test.scratch.rho_out, 1, buf)
    copyto!(out, buf)
end
err = hs_operator_norm_krylov(diff!, adj!, d)
```
"""
function hs_operator_norm_krylov(
    L_apply!::F1,
    L_apply_adj!::F2,
    d::Integer;
    T::Type = ComplexF64,
    tol::Real = 1e-12,
    krylovdim::Int = 30,
    maxiter::Int = 100,
    max_retries::Int = 3,
) where {F1, F2}
    in_buf  = Matrix{T}(undef, d, d)
    out_buf = Matrix{T}(undef, d, d)
    fwd = function (v::AbstractVector)
        copyto!(in_buf, v)
        L_apply!(out_buf, in_buf)
        return copy(vec(out_buf))
    end
    adj = function (v::AbstractVector)
        copyto!(in_buf, v)
        L_apply_adj!(out_buf, in_buf)
        return copy(vec(out_buf))
    end
    x0 = randn(T, d * d)

    # Pre-flight: estimate the operator scale via one fwd/adj pair.  When the
    # operator is below the noise floor of the matvec (i.e.  ‖A* u₀‖ / ‖u₀‖ at
    # the order of the per-call roundoff), KrylovKit's GKL initialiser throws
    # `ArgumentError("operator and its adjoint are not compatible")` because
    # its self-consistency check `α² ≈ α*α` fails at sqrt(eps()) tolerance —
    # see KrylovKit/factorizations/gkl.jl.  In that regime the answer is
    # effectively zero and we return the single-shot lower-bound estimate.
    β₀ = norm(x0)
    v₀ = adj(x0)
    α_lb = norm(v₀) / β₀

    current_kdim = krylovdim
    local vals, info
    for attempt in 1:(max_retries + 1)
        try
            vals, _, _, info = svdsolve(
                (fwd, adj), copy(x0), 1, :LR;
                krylovdim = current_kdim, tol = tol, maxiter = maxiter,
            )
        catch e
            if e isa ArgumentError && occursin("not compatible", e.msg)
                @warn "hs_operator_norm_krylov: GKL self-consistency failed " *
                      "(operator at noise floor: ‖A*u₀‖/‖u₀‖ = $α_lb). " *
                      "Returning the lower-bound estimate."
                return α_lb
            else
                rethrow(e)
            end
        end
        if info.converged >= 1
            return real(vals[1])
        end
        if attempt <= max_retries
            new_kdim = ceil(Int, current_kdim * 1.5)
            @warn "hs_operator_norm_krylov: $(info.converged)/1 converged. " *
                  "Retrying with krylovdim=$new_kdim (attempt $(attempt+1)/$(max_retries+1))"
            current_kdim = new_kdim
        end
    end
    error("hs_operator_norm_krylov failed to converge after $(max_retries + 1) attempts")
end

"""
    dissipator_one_to_one_norm_bound(L_a_list) -> Real

Upper bound on the dissipator 1→1 induced norm `‖L_diss‖_{1→1}` for
`L_diss(ρ) = Σ_a (L_a ρ L_a† - ½{L_a† L_a, ρ})`.

Term-by-term in trace norm:
- `‖L_a ρ L_a†‖_1 ≤ ‖L_a† L_a‖_∞ · ‖ρ‖_1`
- `‖{L_a† L_a, ρ}‖_1 ≤ 2 · ‖L_a† L_a‖_∞ · ‖ρ‖_1`

Summing in absolute value bounds the difference of CP and anti-commutator
parts by `2·‖Σ_a L_a L_a†‖_∞ + 2·‖Σ_a L_a† L_a‖_∞`. For paired symmetric
jump sets (CKG / DLL convention) `Σ_a L_a L_a† = Σ_a L_a† L_a`, so this
collapses to `4·‖Σ_a L_a† L_a‖_∞`.

NOTE the inner quantity is the operator norm of the *sum* `Σ_a L_a† L_a`,
NOT `max_a ‖L_a† L_a‖_∞`. For K rank-1 channels the latter can be up to K×
too small and is not a valid upper bound.
"""
function dissipator_one_to_one_norm_bound(L_a_list::AbstractVector{<:AbstractMatrix})
    @assert !isempty(L_a_list) "L_a_list must contain at least one operator"
    sum_LdagL = sum(L_a' * L_a for L_a in L_a_list)
    return 4.0 * opnorm(sum_LdagL)
end

"""
    dissipator_trace_alpha(α) -> Real

Trace of the Kossakowski matrix, `Tr(α) = sum_i α[i, i]`.  A cheap diagnostic
of the total dissipative weight: for a rank-1 channel built from `L_a` the
diagonal entries of `α` equal `tr(L_a† L_a)`.
"""
dissipator_trace_alpha(α::AbstractMatrix) = real(sum(diag(α)))
