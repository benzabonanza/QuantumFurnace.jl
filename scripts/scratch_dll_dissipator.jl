#!/usr/bin/env julia
# scripts/scratch_dll_dissipator.jl
#
# Phase 51: DLL-2 — DLL dissipator (standalone prototype).
# Beads: qf-3i8.2 (parent epic qf-3i8).
#
# Goal: build the DLL dissipator superoperator (no coherent term — that is
# DLL-3) in BohrDomain and TimeDomain for a small toy system, and verify it
# against the equations of Ding–Li–Lin 2024, Section 3.
#
# Run with:  julia --project scripts/scratch_dll_dissipator.jl
#
# ============================================================================
# DLL PAPER EQUATIONS (Ding–Li–Lin 2024, Sec. 3)
# ----------------------------------------------------------------------------
# Eq. 3.2  (KMS condition on q^a):
#     q^a(-ν) = conj(q^a(ν)),   f̂^a(ν) := q^a(ν) e^{-β ν / 4} ∈ L¹(ℝ)
#
# Eq. 3.3  (time-domain filter via inverse Fourier transform):
#     f^a(t) = (1 / 2π) ∫ q^a(ν) e^{-β ν / 4} e^{-i t ν} dν
#
# Eq. 3.4  (THE key Lindblad-operator formula — three equivalent expressions):
#     L_a = Σ_{ν ∈ B_H} q^a(ν) e^{-β ν / 4} A^a_ν
#         = Σ_ν f̂^a(ν) A^a_ν
#         = ∫_{-∞}^∞ f^a(t) A^a(t) dt,   A^a(t) = e^{i H t} A^a e^{-i H t}
# Remark 12 reads L_a as Â^a_{f^a}(ω = 0) — the OFT of A^a with filter f^a
# evaluated at ω = 0. NO outer ω-loop, unlike CKG.
#
# Eq. 3.8 / 3.9  (master equation, dissipative part only — coherent G is DLL-3):
#     ∂_t ρ = -i[G, ρ] + Σ_a (L_a ρ L_a† − (1/2) {L_a† L_a, ρ})
#
# Eq. 3.13  (uniform quadrature grid for time-domain L_a and G):
#     t_m = -M τ + m τ,   0 ≤ m < 2M,   M = 2^{m-1}
#
# Eq. 3.15  (Proposition 16, quadrature error bound for the Gaussian-type
# filter from Eq. 3.21):
#     ‖L_a − Σ_m f^a(t_m) A^a(t_m) τ‖
#       ≤ C_f · exp(− s ((M-1)τ)^{1/s} / (2 (β A_u + A_w)^{1/s} · e))
# with s = 1/2 for the Gaussian filter (Gevrey order s_u = 1/2; w ≡ 1 makes
# s_w irrelevant). The constants A_u, A_w, C_f are O(1); we choose (M-1)τ
# so the exponent dominates and verify the residual numerically.
#
# Eq. 3.21 / 3.22  (Gaussian-type weighting and resulting filter):
#     q(ν)  = e^{-(β ν)² / 8},   w ≡ 1
#     f̂(ν) = q(ν) e^{-β ν / 4} = e^{1/8} exp(-(β ν + 1)² / 8)
#     f(t)  = (e^{1/8} √(2/π) / β) exp(-2 t² / β² + i t / β)
#
# ============================================================================
# COLUMN-STACKING VECTORISATION CONVENTION (matches src/jump_workers.jl):
#     vec(ρ) is column-stacked (Julia default).
#     L ρ R     ↔  (Rᵀ ⊗ L)            · vec(ρ)
#     {A, ρ}    ↔  (I ⊗ A + Aᵀ ⊗ I)    · vec(ρ)
#     L_a ρ L_a† ↔ ((L_a†)ᵀ ⊗ L_a) · vec(ρ) = (conj(L_a) ⊗ L_a) · vec(ρ)
# Hence the dissipator superoperator (per coupling a) is
#     D_a = conj(L_a) ⊗ L_a − (1/2) [I ⊗ (L_a† L_a) + (L_a† L_a)ᵀ ⊗ I]
# This matches `src/jump_workers.jl::_vectorize_liouv_diss_and_add!`.
# ============================================================================

import Pkg
Pkg.activate(dirname(@__DIR__))

using QuantumFurnace
using LinearAlgebra
using Printf
using Random

# ============================================================================
# Part A: BohrDomain DLL Lindblad operator and dissipator
# ============================================================================

"""
    lindblad_operator_bohr(filter, eigvals, A_eb) -> Matrix{ComplexF64}

DLL Lindblad operator from Eq. 3.4 (Bohr decomposition), in the eigenbasis:

    L^Bohr_a[i, j] = freq_kernel(filter, ν_ij) · A_eb[i, j],  ν_ij = λ_i − λ_j.

For `DLLGaussianFilter`, `freq_kernel(filter, ν) = q(ν) e^{-βν/4}` is the FULL
DLL weighting (Eq. 3.22) — no extra `e^{-βν/4}` factor needed.

Inputs are already in the Hamiltonian eigenbasis: `eigvals` (length n) and
`A_eb` (n×n Hermitian coupling rotated to the eigenbasis).
"""
function lindblad_operator_bohr(
    filter::AbstractFilter,
    eigvals::AbstractVector{<:Real},
    A_eb::AbstractMatrix{<:Complex},
)
    n = length(eigvals)
    @assert size(A_eb) == (n, n)
    L = zeros(ComplexF64, n, n)
    @inbounds for j in 1:n, i in 1:n
        ν_ij = eigvals[i] - eigvals[j]
        L[i, j] = freq_kernel(filter, ν_ij) * A_eb[i, j]
    end
    return L
end

"""
    dissipator_superoperator(L::Matrix) -> Matrix

Vectorised single-coupling dissipator:

    D_a = conj(L) ⊗ L − (1/2) [I ⊗ (L† L) + (L† L)ᵀ ⊗ I]

Acts on `vec(ρ)` (column-stacked). Matches `src/jump_workers.jl`'s convention
byte-for-byte.
"""
function dissipator_superoperator(L::AbstractMatrix{<:Complex})
    n = size(L, 1)
    @assert size(L, 2) == n
    Lc = conj(L)
    LdL = L' * L
    Id = Matrix{ComplexF64}(I, n, n)
    return kron(Lc, L) - 0.5 * (kron(Id, LdL) + kron(transpose(LdL), Id))
end

"""
    build_dll_dissipator_bohr(filter, eigvals, A_list_eb) -> Matrix

Sum the per-coupling dissipators built via `lindblad_operator_bohr`.
"""
function build_dll_dissipator_bohr(
    filter::AbstractFilter,
    eigvals::AbstractVector{<:Real},
    A_list_eb::AbstractVector{<:AbstractMatrix{<:Complex}},
)
    n = length(eigvals)
    D = zeros(ComplexF64, n^2, n^2)
    for A_eb in A_list_eb
        L_a = lindblad_operator_bohr(filter, eigvals, A_eb)
        D .+= dissipator_superoperator(L_a)
    end
    return D
end

# ============================================================================
# Part B: TimeDomain DLL Lindblad operator (OFT @ ω = 0) and dissipator
# ============================================================================

"""
    choose_quadrature_params(filter::DLLGaussianFilter; H_norm, target_err) -> (M, τ)

Pick `(M, τ)` for the trapezoidal quadrature of `L_a^time = ∫ f(t) A(t) dt`
on the uniform grid `t_m = -Mτ + mτ, m=0,…,2M-1` (Eq. 3.13).

Strategy:
  τ chosen via Nyquist-like rule  τ = π / (‖H‖ + S_eff)  where S_eff is a
  Gaussian-tail proxy for the spread of `q(ν) e^{-βν/4}` (≈ 4/β at the
  1e-6 level). This keeps the integrand `f(t) A(t)` well-resolved.

  M chosen so that the Eq. 3.15 exponent drives the bound below
  `target_err`. For s = 1/2, A_u = A_w = 1, the exponent is
      Φ(M, τ) = 0.5 · ((M-1)τ)^2 / (2 (β + 1)^2 · e)
  Solve Φ ≥ -log(target_err) for `(M-1) τ`, then round M up to the next
  power of 2.

# PHYSICS CHECK: A_u = A_w = 1 is a paper-O(1) placeholder; the prefactor
# C_f is also treated as O(1). We verify the residual numerically afterwards
# rather than relying on the bound; if the empirical residual exceeds 1e-4
# we double M (handled by the caller).
"""
function choose_quadrature_params(
    filter::DLLGaussianFilter{T};
    H_norm::Real,
    target_err::Real = 1e-6,
) where {T<:AbstractFloat}
    β = filter.beta
    # τ: Nyquist-like; q(ν) e^{-βν/4} is centred at ν_* = -1/β with σ_ν = 2/β,
    # so its support is roughly [-1/β - 4·2/β, -1/β + 4·2/β] = [-9/β, 7/β].
    # Take S_eff = 9/β (slight over-estimate); add 2·H_norm for the A(t) factor.
    S_eff = 9.0 / β
    τ = π / (2 * H_norm + S_eff)

    # M from Eq. 3.15 exponent (s = 1/2, A_u = A_w = 1):
    #   Φ = 0.5 · ((M-1)τ)^2 / (2 (β + 1)^2 · e)  ≥  -log(target_err)
    target_log = -log(target_err)
    Mtau_min = sqrt(4 * (β + 1)^2 * exp(1) * target_log)  # = (M-1) τ from above
    Mreq = ceil(Int, Mtau_min / τ) + 1  # +1 for the (M-1) shift
    # Round up to next power of 2 (paper M = 2^{m-1})
    M = 2^ceil(Int, log2(max(Mreq, 2)))
    return M, τ
end

"""
    lindblad_operator_time(filter, eigvals, A_eb; M, τ) -> Matrix

DLL Lindblad operator via discrete OFT @ ω = 0 (Eq. 3.4 third equality):

    L^time_a = Σ_{m=0}^{2M-1} f(t_m) · D(t_m) · A_eb · D(t_m)† · τ,

with `t_m = -M τ + m τ` and `D(t) = Diagonal(exp(i λ_k t))`. This is exactly
the trapezoidal quadrature on the uniform grid from Eq. 3.13.

`time_kernel(::DLLGaussianFilter, t)` returns the FULL `f(t)` of Eq. 3.22
including the `e^{1/8} √(2/π)/β` prefactor — no extra normalisation needed.
"""
function lindblad_operator_time(
    filter::AbstractFilter,
    eigvals::AbstractVector{<:Real},
    A_eb::AbstractMatrix{<:Complex};
    M::Int,
    τ::Real,
)
    n = length(eigvals)
    @assert size(A_eb) == (n, n)
    L = zeros(ComplexF64, n, n)
    # Pre-compute t_m grid; we accumulate f(t_m) · D(t_m) A_eb D(t_m)† · τ.
    @inbounds for m in 0:(2 * M - 1)
        t = (m - M) * τ  # = -M τ + m τ
        ft = time_kernel(filter, t)
        # D(t) A_eb D(t)†: element [i, j] = exp(i (λ_i − λ_j) t) · A_eb[i, j].
        for j in 1:n, i in 1:n
            phase = cis((eigvals[i] - eigvals[j]) * t)
            L[i, j] += ft * phase * A_eb[i, j] * τ
        end
    end
    return L
end

"""
    build_dll_dissipator_time(filter, eigvals, A_list_eb; M, τ) -> Matrix

Time-domain analogue of `build_dll_dissipator_bohr`.
"""
function build_dll_dissipator_time(
    filter::AbstractFilter,
    eigvals::AbstractVector{<:Real},
    A_list_eb::AbstractVector{<:AbstractMatrix{<:Complex}};
    M::Int,
    τ::Real,
)
    n = length(eigvals)
    D = zeros(ComplexF64, n^2, n^2)
    for A_eb in A_list_eb
        L_a = lindblad_operator_time(filter, eigvals, A_eb; M=M, τ=τ)
        D .+= dissipator_superoperator(L_a)
    end
    return D
end

# ============================================================================
# Symmetric toy filter for Verification 4
# ============================================================================

"""
    SymmetricToyFilter{T}(width::T)

Hand-crafted KMS-VIOLATING filter with `freq_kernel(f, ν) = exp(-(ν/width)²)`
— real, even in ν. Useful for the hermiticity sanity check (Task 11): with
`A_a` Hermitian and a real even kernel, the resulting `L_a` is Hermitian.
"""
struct SymmetricToyFilter{T<:AbstractFloat} <: AbstractFilter
    width::T
end
@inline QuantumFurnace.freq_kernel(f::SymmetricToyFilter{T}, nu::Real) where {T} =
    exp(-(nu / f.width)^2)
# time_kernel/filter_time_cutoff not implemented — only used at Bohr level here.

# ============================================================================
# DLL Kossakowski preview (Part D)
# ============================================================================

"""
    dll_kossakowski(filter::DLLGaussianFilter, bohr_freqs) -> Matrix

For unique Bohr frequencies `{ν_k}`, compute the DLL Kossakowski matrix per
single coupling A^a:

    α^DLL[k, l] = e^{-β(ν_k+ν_l)/4} · q(ν_k) · conj(q(ν_l)).

By construction this is a rank-1 outer product `v vᵀ*`, with
`v_k = e^{-βν_k/4} q(ν_k) = freq_kernel(DLL, ν_k)`.
"""
function dll_kossakowski(filter::DLLGaussianFilter, bohr_freqs::AbstractVector{<:Real})
    v = ComplexF64.(freq_kernel.(Ref(filter), bohr_freqs))  # length-K vector
    return v * v'  # K × K Hermitian rank-1
end

"""
    ckg_kossakowski(σ::Real, gp, bohr_freqs) -> Matrix

CKG Gaussian Kossakowski via `create_alpha_gauss(ν_k, ν_l, σ, gp)` for each
pair (Bohr-domain helper from `src/bohr_domain.jl`). `gp = (w_γ, σ_γ)`.
"""
function ckg_kossakowski(
    σ::Real,
    gp::Tuple{<:Real, <:Real},
    bohr_freqs::AbstractVector{<:Real},
)
    K = length(bohr_freqs)
    α = zeros(ComplexF64, K, K)
    for l in 1:K, k in 1:K
        α[k, l] = create_alpha_gauss(bohr_freqs[k], bohr_freqs[l], σ, gp)
    end
    return α
end

# ============================================================================
# Helpers
# ============================================================================

"""
    eigen_setup(H) -> (λ, U, bohr_dict)

Eigendecomposition of `H`, plus the Bohr frequency dictionary keyed by
unique pairwise differences λ_i − λ_j.
"""
function eigen_setup(H::AbstractMatrix{<:Complex})
    F = eigen(Hermitian(Matrix(H)))
    λ, U = F.values, F.vectors
    bohr_freqs_mat = λ .- transpose(λ)
    bohr_dict = create_bohr_dict(bohr_freqs_mat)
    return λ, U, bohr_dict
end

"""
    gibbs_diag(λ, β) -> Vector

Diagonal Gibbs state in the eigenbasis: `exp(-β λ_i)/Z`.
"""
function gibbs_diag(λ::AbstractVector{<:Real}, β::Real)
    weights = exp.(-β .* λ)
    return weights / sum(weights)
end

"""
    fro(M) -> Float64

Convenience: Frobenius norm of a matrix as a plain Float64.
"""
fro(M::AbstractMatrix) = norm(M, 2)  # Julia's `norm(matrix)` is Frobenius

"""
    op_norm(M) -> Float64

Operator (spectral) norm. Robust against tiny non-Hermitian noise.
"""
op_norm(M::AbstractMatrix) = opnorm(Matrix{ComplexF64}(M), 2)

# ============================================================================
# MAIN
# ============================================================================

function run_for_beta(β::Float64; verbose::Bool=true)
    verbose && println("="^76)
    verbose && println("Running DLL-2 prototype at β = $β")
    verbose && println("="^76)

    # ----- Toy system ------------------------------------------------------
    Random.seed!(20260430)
    n_qubits = 2
    dim = 4

    # 4×4 Hermitian Hamiltonian with non-degenerate Bohr spectrum.
    # Use the same XXZ + h_a Z + h_b Z form as scratch_dll_filter.jl so
    # the eigenstructure is well-understood.
    X = ComplexF64[0 1; 1 0]
    Y = ComplexF64[0 -im; im 0]
    Z = ComplexF64[1 0; 0 -1]
    Id2 = ComplexF64[1 0; 0 1]
    H = kron(X, X) + kron(Y, Y) + 1.5 * kron(Z, Z) +
        0.7 * kron(Z, Id2) + 0.3 * kron(Id2, Z)
    H = Matrix(Hermitian(H))
    H_norm = opnorm(H, 2)

    # Couplings: three Hermitian Pauli-style operators padded to dim 4.
    # PHYSICS CHECK: Hermitian A^a is required for the DLL DBC sanity in
    # Verification 4 (symmetric-filter hermiticity); the dissipator at Eq. 3.4
    # only requires bounded A^a, but Hermitian is the natural physical choice.
    A_list = Matrix{ComplexF64}[
        Matrix(kron(X, Id2)),
        Matrix(kron(Id2, Y)),
        Matrix(kron(Z, Z)),
    ]
    A_labels = ["X⊗I", "I⊗Y", "Z⊗Z"]

    # Eigendecomposition + Bohr dictionary.
    λ, U, bohr_dict = eigen_setup(H)
    verbose && @printf("  H_norm = %.6f, ‖λ‖∞ = %.6f, dim = %d\n", H_norm, maximum(abs, λ), dim)

    # Rotate couplings into the eigenbasis once.
    A_list_eb = [U' * A * U for A in A_list]

    # Gibbs state in eigenbasis (diagonal).
    σdiag = gibbs_diag(λ, β)
    σβ = Diagonal(ComplexF64.(σdiag))   # n×n
    σβ_mat = Matrix(σβ)
    vec_σβ = vec(σβ_mat)

    # Identity (for trace preservation in dual).
    Id_mat = Matrix{ComplexF64}(I, dim, dim)
    vec_Id = vec(Id_mat)

    # ----- Filter ----------------------------------------------------------
    filter = DLLGaussianFilter(β)

    # ----- Part A: Bohr domain --------------------------------------------
    verbose && println()
    verbose && println("[Part A] BohrDomain L_a (Eq. 3.4) — kernel q(ν)e^{-βν/4} = freq_kernel(DLL, ν).")
    L_bohr = [lindblad_operator_bohr(filter, λ, A_eb) for A_eb in A_list_eb]
    if verbose
        for (a, lbl) in enumerate(A_labels)
            @printf("    L^Bohr_a [%s] : ‖L‖_F = %.6e, ‖L‖_op = %.6e\n",
                    lbl, fro(L_bohr[a]), op_norm(L_bohr[a]))
        end
    end
    D_bohr = build_dll_dissipator_bohr(filter, λ, A_list_eb)
    verbose && @printf("    ‖D_DLL^Bohr‖_F = %.6e\n", fro(D_bohr))

    # ----- Part B: Time domain --------------------------------------------
    verbose && println()
    verbose && println("[Part B] TimeDomain L_a via OFT @ ω=0 (Eq. 3.4 third equality, Eq. 3.13).")
    M, τ = choose_quadrature_params(filter; H_norm=H_norm, target_err=1e-6)
    verbose && @printf("    initial (M, τ) from Eq. 3.15: M = %d, τ = %.6f, span = ±%.6f\n",
                       M, τ, M * τ)

    # Refinement loop: doubling M until per-coupling Bohr↔Time op-norm error
    # is ≤ target_residual or we hit the cap.
    target_residual = 1e-6
    cap_M = 8192  # safety cap for the prototype
    L_time = Vector{Matrix{ComplexF64}}(undef, length(A_list))
    while true
        L_time = [lindblad_operator_time(filter, λ, A_eb; M=M, τ=τ) for A_eb in A_list_eb]
        residual_op = maximum(op_norm(L_time[a] - L_bohr[a]) for a in eachindex(A_list))
        if residual_op ≤ target_residual || M ≥ cap_M
            verbose && @printf("    settled at M = %d, max op-norm residual = %.3e\n",
                               M, residual_op)
            break
        end
        verbose && @printf("    doubling M: %d -> %d (residual was %.3e)\n", M, 2M, residual_op)
        M *= 2
    end
    D_time = build_dll_dissipator_time(filter, λ, A_list_eb; M=M, τ=τ)
    verbose && @printf("    ‖D_DLL^time‖_F = %.6e\n", fro(D_time))

    # ----- Part C verifications -------------------------------------------
    verbose && println()
    verbose && println("[Part C] Verifications (DLL-paper-tracked tolerances)")

    # (1) Bohr ↔ Time per-coupling consistency (Eq. 3.4 third equality)
    bohr_time_op = Float64[]
    bohr_time_fro = Float64[]
    for a in eachindex(A_list)
        push!(bohr_time_op, op_norm(L_time[a] - L_bohr[a]))
        push!(bohr_time_fro, fro(L_time[a] - L_bohr[a]))
    end
    bohr_time_op_max = maximum(bohr_time_op)
    bohr_time_fro_max = maximum(bohr_time_fro)
    if verbose
        for (a, lbl) in enumerate(A_labels)
            @printf("    [1] Bohr↔Time [%s]: ‖ΔL‖_op = %.3e, ‖ΔL‖_F = %.3e\n",
                    lbl, bohr_time_op[a], bohr_time_fro[a])
        end
    end
    pass_1 = bohr_time_op_max ≤ 1e-4
    verbose && @printf("    [1] max op-norm = %.3e  =>  %s (tol 1e-4)\n",
                       bohr_time_op_max, pass_1 ? "PASS" : "FAIL")

    # (2) Hermitian Gibbs fixed point: ‖D[σ_β]‖_F
    img_bohr = D_bohr * vec_σβ
    img_time = D_time * vec_σβ
    fix_bohr = norm(img_bohr)
    fix_time = norm(img_time)
    pass_2_bohr = fix_bohr ≤ 1e-10
    pass_2_time = fix_time ≤ 1e-4
    verbose && @printf("    [2] Bohr  ‖D[σ_β]‖_F = %.3e  =>  %s (tol 1e-10)\n",
                       fix_bohr, pass_2_bohr ? "PASS" : "FAIL")
    verbose && @printf("    [2] Time  ‖D[σ_β]‖_F = %.3e  =>  %s (tol 1e-4)\n",
                       fix_time, pass_2_time ? "PASS" : "FAIL")

    # (3) Trace preservation in dual: ‖D†[I]‖_F
    # The vectorised adjoint of the superoperator is `D'` (Hermitian conjugate
    # of the dim²×dim² matrix). For the dissipator `D_a = conj(L) ⊗ L − ...`,
    # the adjoint reads `D_a† = transpose(L) ⊗ L† − (1/2) [I ⊗ (L†L) + (L†L)ᵀ ⊗ I]`,
    # which on `vec(I)` returns vec(L†IL − (1/2)(L†L · I + I · L†L)) = 0 by
    # construction. We use the Hermitian conjugate of the matrix directly.
    dual_bohr = D_bohr' * vec_Id
    dual_time = D_time' * vec_Id
    tp_bohr = norm(dual_bohr)
    tp_time = norm(dual_time)
    pass_3_bohr = tp_bohr ≤ 1e-10
    pass_3_time = tp_time ≤ 1e-4
    verbose && @printf("    [3] Bohr  ‖D†[I]‖_F = %.3e  =>  %s (tol 1e-10)\n",
                       tp_bohr, pass_3_bohr ? "PASS" : "FAIL")
    verbose && @printf("    [3] Time  ‖D†[I]‖_F = %.3e  =>  %s (tol 1e-4)\n",
                       tp_time, pass_3_time ? "PASS" : "FAIL")

    # (4) Symmetric-filter sanity: with f̂(ν) = exp(-(ν/w)²) (real, even),
    # L_a should be Hermitian for Hermitian A_a.
    sym_filter = SymmetricToyFilter(2.0 / β)  # similar width to DLL filter
    herm_residuals = Float64[]
    for (a, A_eb) in enumerate(A_list_eb)
        Ls = lindblad_operator_bohr(sym_filter, λ, A_eb)
        push!(herm_residuals, fro(Ls - Ls'))
    end
    herm_residual_max = maximum(herm_residuals)
    pass_4 = herm_residual_max ≤ 1e-12
    if verbose
        for (a, lbl) in enumerate(A_labels)
            @printf("    [4] symmetric filter [%s]: ‖L − L†‖_F = %.3e\n",
                    lbl, herm_residuals[a])
        end
        @printf("    [4] max ‖L − L†‖_F = %.3e  =>  %s (tol 1e-12)\n",
                herm_residual_max, pass_4 ? "PASS" : "FAIL")
    end

    # ----- Part D: Kossakowski preview ------------------------------------
    verbose && println()
    verbose && println("[Part D] DLL Kossakowski preview vs CKG (per coupling)")

    # Unique Bohr frequencies (sorted, as a deterministic ordering).
    unique_bohr = sort(collect(keys(bohr_dict)))
    K = length(unique_bohr)
    verbose && @printf("    K unique Bohr frequencies = %d\n", K)

    # DLL: rank-1 outer product per coupling.
    α_DLL = dll_kossakowski(filter, unique_bohr)
    sv_DLL = svdvals(α_DLL)
    rank_DLL = count(sv -> sv > 1e-12 * sv_DLL[1], sv_DLL)
    fro_DLL = norm(α_DLL)

    # CKG: full Gaussian Kossakowski with σ = 2/β (matched width to f̂(ν))
    # and (w_γ, σ_γ) = (0, σ).
    σ_match = 2.0 / β
    α_CKG = ckg_kossakowski(σ_match, (0.0, σ_match), unique_bohr)
    sv_CKG = svdvals(α_CKG)
    rank_CKG = count(sv -> sv > 1e-12 * sv_CKG[1], sv_CKG)
    fro_CKG = norm(α_CKG)

    if verbose
        @printf("    DLL  α: rank(ε=1e-12) = %d, top 5 sv = [%s], ‖α‖_F = %.6e\n",
                rank_DLL,
                join(map(s -> @sprintf("%.3e", s), sv_DLL[1:min(5, K)]), ", "),
                fro_DLL)
        @printf("    CKG  α: rank(ε=1e-12) = %d, top 5 sv = [%s], ‖α‖_F = %.6e\n",
                rank_CKG,
                join(map(s -> @sprintf("%.3e", s), sv_CKG[1:min(5, K)]), ", "),
                fro_CKG)
        @printf("    Spread (sv_max / sv_min) DLL: %.3e (rank-1 ⇒ effectively ∞)\n",
                sv_DLL[1] / max(sv_DLL[end], eps()))
        @printf("    Spread (sv_max / sv_min) CKG: %.3e\n",
                sv_CKG[1] / max(sv_CKG[end], eps()))
    end
    pass_D = (rank_DLL == 1)
    verbose && @printf("    [D] rank(α^DLL) = %d (expect 1)  =>  %s\n",
                       rank_DLL, pass_D ? "PASS" : "FAIL")

    return (
        β = β,
        H_norm = H_norm,
        A_labels = A_labels,
        L_bohr = L_bohr,
        L_time = L_time,
        D_bohr = D_bohr,
        D_time = D_time,
        M = M,
        τ = τ,
        bohr_time_op_max = bohr_time_op_max,
        bohr_time_fro_max = bohr_time_fro_max,
        fix_bohr = fix_bohr,
        fix_time = fix_time,
        tp_bohr = tp_bohr,
        tp_time = tp_time,
        herm_residual_max = herm_residual_max,
        rank_DLL = rank_DLL,
        rank_CKG = rank_CKG,
        fro_DLL = fro_DLL,
        fro_CKG = fro_CKG,
        sv_DLL = sv_DLL,
        sv_CKG = sv_CKG,
        K = K,
        passes = (pass_1, pass_2_bohr, pass_2_time, pass_3_bohr, pass_3_time, pass_4, pass_D),
    )
end

function main()
    println("="^76)
    println("Phase 51 / DLL-2: dissipator superoperator prototype")
    println("scripts/scratch_dll_dissipator.jl")
    println("="^76)
    println()

    results = Dict{Float64, Any}()
    for β in (1.0, 5.0)
        results[β] = run_for_beta(β; verbose=true)
        println()
    end

    # ---- Final SUMMARY ----------------------------------------------------
    println("="^76)
    println("DLL-2 DISSIPATOR PROTOTYPE — SUMMARY")
    println("="^76)
    for β in (1.0, 5.0)
        r = results[β]
        @printf("β = %.1f\n", β)
        # Part A summary
        fro_str = join([@sprintf("%.3e", fro(L)) for L in r.L_bohr], ", ")
        @printf("  Part A (Bohr): L_a built for {%s}, dim 4. Frobenius norms: [%s].\n",
                join(r.A_labels, ", "), fro_str)
        @printf("  Part B (Time): M = %d, τ = %.6f, quadrature span = ±%.6f. Per-jump Bohr↔Time op-norm error ≤ %.3e.\n",
                r.M, r.τ, r.M * r.τ, r.bohr_time_op_max)
        @printf("  Part C verifications:\n")
        @printf("     [1] Bohr↔Time   max op-norm = %.3e (tol 1e-4)         %s\n",
                r.bohr_time_op_max, r.passes[1] ? "PASS" : "FAIL")
        @printf("     [2] Gibbs fix   Bohr ‖D[σ]‖_F = %.3e (tol 1e-10) %s   |  Time = %.3e (tol 1e-4) %s\n",
                r.fix_bohr, r.passes[2] ? "PASS" : "FAIL",
                r.fix_time, r.passes[3] ? "PASS" : "FAIL")
        @printf("     [3] Trace pres. Bohr ‖D†[I]‖_F = %.3e (tol 1e-10) %s | Time = %.3e (tol 1e-4) %s\n",
                r.tp_bohr, r.passes[4] ? "PASS" : "FAIL",
                r.tp_time, r.passes[5] ? "PASS" : "FAIL")
        @printf("     [4] Sym-filter herm  ‖L − L†‖_F = %.3e (tol 1e-12)    %s\n",
                r.herm_residual_max, r.passes[6] ? "PASS" : "FAIL")
        @printf("  Part D Kossakowski preview: rank(α^DLL) = %d per coupling (expect 1), ‖α^DLL‖_F = %.6e ; rank(α^CKG) = %d, ‖α^CKG‖_F = %.6e\n",
                r.rank_DLL, r.fro_DLL, r.rank_CKG, r.fro_CKG)
        all_pass = all(r.passes)
        @printf("  ALL CHECKS β=%.1f: %s\n", β, all_pass ? "PASS" : "FAIL")
        println()
    end
    println("="^76)
    overall = all(all(results[β].passes) for β in (1.0, 5.0))
    println(overall ? "ALL VERIFICATIONS PASSED." : "SOME VERIFICATIONS FAILED.")
    println("="^76)

    return overall ? 0 : 1
end

# Exit with appropriate status code.
exit_code = main()
exit(exit_code)
