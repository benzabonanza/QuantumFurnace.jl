# KMS geometry diagnostics for comparing generator scales and decay rates.
#
# Math: $inner(X,Y)_KMS = tr(X^dagger sigma^(1/2) Y sigma^(1/2))$.
# The Hermitian part of the quantum discriminant determines the Dirichlet
# rates; the coherent commutator contributes only to its anti-Hermitian part.
# These dense diagnostics require a Gibbs state diagonal in the working basis.

"""
    kms_inner_product(X, Y, sigma; sigma_sqrt=nothing) -> Complex

Return the KMS inner product of `X` and `Y` with respect to `sigma`.

# Keywords
- `sigma_sqrt`: Optional cached square root of `sigma`.
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

Return the nonnegative KMS norm of `X`.
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

Return the KMS variance of `X` about its Gibbs expectation.
"""
function kms_variance(
    X::AbstractMatrix,
    sigma::AbstractMatrix;
    sigma_sqrt::Union{Nothing, AbstractMatrix} = nothing,
)
    s12 = sigma_sqrt === nothing ? sqrt(Hermitian(sigma)) : sigma_sqrt
    # Math: $mean(X) = tr(sigma X)$ for Hermitian `X`.
    mean_X = real(tr(sigma * X))
    Y = X .- mean_X .* I(size(X, 1))
    return real(tr(Y' * s12 * Y * s12))
end


"""
    build_dense_superoperator(L_apply!, d; T=ComplexF64) -> Matrix{T}

Materialise a column-stacked `d^2 × d^2` superoperator.

# Arguments
- `L_apply!`: In-place matrix action `L_apply!(out, X)`.
- `d`: Hilbert-space dimension.

# Keywords
- `T`: Superoperator element type.

# Returns
The dense column-stacked superoperator.
"""
function build_dense_superoperator(L_apply!::F, d::Integer; T = ComplexF64) where {F}
    out = zeros(T, d * d, d * d)
    in_buf  = zeros(T, d, d)
    out_buf = zeros(T, d, d)
    for col in 1:(d * d)
        fill!(in_buf, 0)
        in_buf[col] = one(T)
        L_apply!(out_buf, in_buf)
        @views out[:, col] .= vec(out_buf)
    end
    return out
end

"""
    _kms_discriminant_with_sym(L_super, sigma) -> (D_super, D_sym)

Return the KMS discriminant and its Hermitian part.
"""
function _kms_discriminant_with_sym(
    L_super::AbstractMatrix{<:Complex},
    sigma::AbstractMatrix,
)
    D_super = materialize_discriminant(L_super, Hermitian(Matrix(sigma)))
    D_sym   = Hermitian((D_super .+ D_super') ./ 2)
    return D_super, D_sym
end


"""
    kms_dirichlet_form(L_super_H, X, sigma; sigma_sqrt=nothing) -> Real

Return the KMS Dirichlet form for a Heisenberg-picture superoperator.

Pass the conjugate transpose of a column-stacked Schrödinger generator as
`L_super_H`.
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


"""
    _kms_dirichlet_eigvals(L_super, sigma; project_constants=true) -> Vector{Float64}

Return the real eigenvalues of the negative Hermitian discriminant part.

`project_constants=true` removes the stationary `sigma^(1/2)` direction.
"""
function _kms_dirichlet_eigvals(
    L_super::AbstractMatrix{<:Complex},
    sigma::AbstractMatrix;
    project_constants::Bool = true,
)
    _, D_sym = _kms_discriminant_with_sym(L_super, sigma)
    M = -Matrix(D_sym)
    if project_constants
        # Math: the zero direction is $vec(sigma^(1/2))$.
        s12 = sqrt(Hermitian(Matrix(sigma)))
        v = vec(s12) ./ norm(vec(s12))
        # Math: project with $P = I - v v^dagger$.
        Mp = M .- v * (v' * M) .- (M * v) * v' .+ v * (v' * (M * v)) * v'
        Mp = Hermitian((Mp .+ Mp') ./ 2)
        ev = eigvals(Mp)
        ev_sorted = sort(real.(ev); by = abs)
        return ev_sorted[2:end]
    else
        return sort(real.(eigvals(Hermitian((M .+ M') ./ 2))); by = abs)
    end
end

"""
    spectral_gap_kms(L_super, sigma; project_constants=true) -> NamedTuple

Return the dense KMS–Poincaré gap and all nonconstant Dirichlet rates.

# Keywords
- `project_constants`: Remove the stationary direction before diagonalisation.

# Returns
`(; gap, eigvals)`, sorted from slowest to fastest rate. The dense calculation
costs `O(d^6)`; use `krylov_spectral_gap` at larger dimension.
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

Return the largest nonconstant KMS Dirichlet rate.
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

Return the scale-invariant diagnostic `gap / maximum_dirichlet_rate`.

This is an internal spectral-spread diagnostic, not a standard sampler
comparison metric.
"""
function intrinsic_mixing_ratio(
    L_super::AbstractMatrix{<:Complex},
    sigma::AbstractMatrix,
)
    ev = _kms_dirichlet_eigvals(L_super, sigma)
    return abs(ev[1]) / maximum(abs.(ev))
end


"""
    hs_operator_norm(L_super) -> Real

Return the Hilbert–Schmidt-induced norm of a dense superoperator.
"""
hs_operator_norm(L_super::AbstractMatrix) = opnorm(L_super)

"""
    hs_operator_norm_krylov(L_apply!, L_apply_adj!, d; ...) -> Real

Estimate the Hilbert–Schmidt-induced norm by matrix-free GKL iteration.

# Arguments
- `L_apply!`, `L_apply_adj!`: Forward and Hilbert–Schmidt-adjoint matrix actions.
- `d`: Hilbert-space dimension.

# Keywords
- `T`: Working scalar type.
- `tol`, `krylovdim`, `maxiter`: GKL solver controls.
- `max_retries`: Retries with an enlarged Krylov subspace.

# Returns
The largest singular value of the implicit superoperator.
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

    # Near the matvec noise floor, return a one-shot lower-bound estimate
    # because GKL's forward/adjoint compatibility check loses significance.
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

Return the paired-jump bound `4 * opnorm(sum(L_a' * L_a))` on the dissipator's
trace-norm-induced operator norm.
"""
function dissipator_one_to_one_norm_bound(L_a_list::AbstractVector{<:AbstractMatrix})
    @assert !isempty(L_a_list) "L_a_list must contain at least one operator"
    sum_LdagL = sum(L_a' * L_a for L_a in L_a_list)
    return 4.0 * opnorm(sum_LdagL)
end

"""
    dissipator_trace_alpha(α) -> Real

Return the real trace of a Kossakowski matrix.
"""
dissipator_trace_alpha(α::AbstractMatrix) = real(sum(diag(α)))
