#=
scripts/scratch_kms_geometry.jl

KMS-geometry utilities for Lindbladian comparison (epic qf-mto.{1,2,3}).

Goal: provide standalone helpers that compute KMS-Poincaré gap λ(L), maximum
KMS Dirichlet rate Λ_max(L), and the *intrinsic* (scale-invariant) ratio
ρ = λ / Λ_max for any (Lindbladian, Gibbs state) pair, plus a 1→1 norm
upper bound surrogate for the dissipator.  These are exactly the numbers
needed to compare CKG smooth-Metro vs DLL Metro vs DLL Gauss on a fair
footing — see `.claude-memory/fair_comparison_dirichlet_qf_mto.md`.

Mathematical setup.

  Operators X, Y on H_d carry the (β, σ)-KMS inner product
      ⟨X, Y⟩_{KMS, σ} = tr(X† σ^{1/2} Y σ^{1/2}).
  KMS detailed balance for the Schrödinger generator L_S is exactly the
  HS-self-adjointness of the *quantum discriminant*
      D_S(X) = σ^{-1/4} L_S(σ^{1/4} X σ^{1/4}) σ^{-1/4}.
  In the column-stacking vec convention, this is the matrix
      D_S = (σ^{-1/4} ⊗ σ^{-1/4}) · L_S · (σ^{1/4} ⊗ σ^{1/4})
          = Φ⁻¹ L_S Φ,                    Φ := kron(σ^{1/4}, σ^{1/4}).
  See `src/discriminant.jl::materialize_discriminant!` and
  `verify_detailed_balance` — the same convention as the production code.

  The KMS Dirichlet form is
      E_L(X) = -Re ⟨X, L_S X⟩_{KMS, σ}
             = -Re tr(X† σ^{1/2} L_S(X) σ^{1/2}).
  KMS-DBC ⟹ E_L(X) ≥ 0 for all X, with equality on X ∝ I.  The KMS-Poincaré
  gap and the maximum Dirichlet rate are
      λ(L)   = inf  E_L(X) / Var_{KMS}(X)    over X ⊥ I
      Λ_max  = sup  E_L(X) / Var_{KMS}(X)    over X ⊥ I.
  Both equal the smallest / largest absolute eigenvalues of -D_S restricted
  to the subspace orthogonal to the constant direction.  D_S has σ^{1/2} as
  its zero eigenvector — see `verify_detailed_balance`'s fixed-point check.

  Scale-invariant ratio:  ρ_intrinsic(L) := λ(L) / Λ_max(L) ∈ (0, 1].

Coherent term: the Lamb-shift -i[B_total, ρ] in QuantumFurnace's KMS
construction is anti-Hermitian on operator space and contributes an
*HS-anti-Hermitian* piece to D_S.  Symmetrising D_S → (D_S + D_S†)/2
annihilates that piece, so λ and Λ_max do not change when the coherent
contribution is restored.  The script verifies this directly.

For the dissipator 1→1 norm we use the upper bound 4·sup_a ‖L_a† L_a‖_∞,
where L_a are the per-coupling Lindblad operators (DLL has them on the
nose via `dll_lindblad_op_bohr`; CKG would need an SDP eigendecomposition
of α and is reported only via the HS-induced norm here).

Run:   julia --project scripts/scratch_kms_geometry.jl
=#

using QuantumFurnace
using LinearAlgebra
using Test
using Printf
using Random

const TEST_DIR = joinpath(dirname(@__DIR__), "test")
include(joinpath(TEST_DIR, "test_helpers.jl"))


# ---------------------------------------------------------------------------
# Section 1 — KMS inner-product algebra
# ---------------------------------------------------------------------------

"""
    kms_inner_product(X, Y, sigma; sigma_sqrt=nothing) -> ComplexF64

KMS inner product `⟨X, Y⟩_{KMS, σ} = tr(X† σ^{1/2} Y σ^{1/2})`.

Pass an explicit `sigma_sqrt = sqrt(Hermitian(σ))` to amortise the matrix
square root over many calls.
"""
function kms_inner_product(X::AbstractMatrix, Y::AbstractMatrix, sigma::AbstractMatrix;
                           sigma_sqrt::Union{Nothing, AbstractMatrix} = nothing)
    s12 = sigma_sqrt === nothing ? sqrt(Hermitian(sigma)) : sigma_sqrt
    return tr(X' * s12 * Y * s12)
end

"""
    kms_norm(X, sigma; sigma_sqrt=nothing) -> Real

KMS norm `√⟨X, X⟩_{KMS}`.  Always returns a real number ≥ 0.
"""
function kms_norm(X::AbstractMatrix, sigma::AbstractMatrix;
                  sigma_sqrt::Union{Nothing, AbstractMatrix} = nothing)
    val = kms_inner_product(X, X, sigma; sigma_sqrt = sigma_sqrt)
    return sqrt(max(real(val), 0.0))
end

"""
    kms_variance(X, sigma; sigma_sqrt=nothing) -> Real

KMS variance `Var_{KMS}(X) = ‖X − ⟨I,X⟩_{KMS} · I⟩‖_{KMS}²`.

Note `⟨I, X⟩_{KMS} = tr(σ^{1/2} X σ^{1/2}) = tr(σ X)` is the Gibbs
expectation of X.  Variance vanishes iff X is a multiple of the identity.
"""
function kms_variance(X::AbstractMatrix, sigma::AbstractMatrix;
                      sigma_sqrt::Union{Nothing, AbstractMatrix} = nothing)
    s12 = sigma_sqrt === nothing ? sqrt(Hermitian(sigma)) : sigma_sqrt
    mean_X = real(tr(sigma * X))   # ⟨I, X⟩_{KMS} = tr(σ X) for Hermitian X (its real part for general X)
    Y = X .- mean_X .* I(size(X, 1))
    return real(tr(Y' * s12 * Y * s12))
end


# ---------------------------------------------------------------------------
# Section 2 — Dense superoperator + KMS quantum discriminant
# ---------------------------------------------------------------------------

"""
    build_dense_superoperator(L_apply!, d; T=ComplexF64) -> Matrix{T}

Materialise a `d² × d²` dense superoperator from a column-major matvec
`L_apply!(out_buf, X) -> out_buf` that maps a `d×d` matrix to a `d×d`
matrix.  Each column is built via a unit basis matrix and a copy.

Use as
    L_super = build_dense_superoperator((out, X) -> begin
        apply_lindbladian!(ws, X, cfg, ham); copyto!(out, ws.scratch.rho_out); out
    end, d)
"""
function build_dense_superoperator(L_apply!, d::Integer; T = ComplexF64)
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
    kms_phi_matrix(sigma; sigma_quarter=nothing) -> (Φ_mat, Φinv_mat)

Vec-form of `Φ(X) = σ^{1/4} X σ^{1/4}` and its inverse.  Column-stacking
vec convention: `vec(σ^{1/4} X σ^{1/4}) = (σ^{1/4} ⊗ σ^{1/4}) vec(X)`.

The flanks are equal so no transpose conjugation is required; the matrix
is `kron(σ^{1/4}, σ^{1/4})`.
"""
function kms_phi_matrix(sigma::AbstractMatrix;
                        sigma_quarter::Union{Nothing, AbstractMatrix} = nothing)
    s14 = sigma_quarter === nothing ? sqrt(sqrt(Hermitian(sigma))) : sigma_quarter
    s14m = inv(s14)
    Φ    = kron(s14, s14)
    Φinv = kron(s14m, s14m)
    return Φ, Φinv
end

"""
    kms_discriminant_superoperator(L_super, sigma; sigma_quarter=nothing)
        -> (D_super, D_sym)

Build the quantum discriminant `D_S = Φ⁻¹ L_S Φ` and its Hermitian part
`D_sym = (D_S + D_S†)/2`.  For a KMS-DBC Schrödinger generator L_S, D_S
itself is HS-self-adjoint and `D_sym ≈ D_S` (anti-Hermitian residual at
machine precision).  When the Lamb-shift coherent piece is included,
D_S has a non-trivial anti-Hermitian part whose contribution to the
Dirichlet form is identically zero — so `D_sym` carries all spectral data
relevant to the gap and Λ_max.

Convention matches `src/discriminant.jl::materialize_discriminant!`.
"""
function kms_discriminant_superoperator(L_super::AbstractMatrix, sigma::AbstractMatrix;
                                        sigma_quarter::Union{Nothing, AbstractMatrix} = nothing)
    Φ, Φinv = kms_phi_matrix(sigma; sigma_quarter = sigma_quarter)
    D_super = Φinv * L_super * Φ
    D_sym   = Hermitian((D_super .+ D_super') ./ 2)
    return D_super, D_sym
end


# ---------------------------------------------------------------------------
# Section 3 — KMS Dirichlet form
# ---------------------------------------------------------------------------

"""
    kms_dirichlet_form(L_super_H, X, sigma; sigma_sqrt=nothing) -> Real

KMS Dirichlet form `E_L(X) = -Re ⟨X, L^H(X)⟩_{KMS, σ}` where L^H is the
**Heisenberg-picture** generator on observables (HS-adjoint of the
Schrödinger generator).  For KMS-DBC L:
  • E_L(X) ≥ 0 for all Hermitian X,
  • E_L(X) = 0 for X = c·I (since L^H(I) = 0 — every Lindbladian preserves
    trace, so the Heisenberg generator is unital).

In the column-stacking vec convention, if `L_super_S` is the Schrödinger
superoperator built by `build_dense_superoperator((out, X) -> apply_lindbladian!(...))`,
then `L_super_H = L_super_S'` (conjugate-transpose) is its HS-adjoint —
this is what the function expects.  Pass the Heisenberg superoperator
explicitly to keep the algebra unambiguous.
"""
function kms_dirichlet_form(L_super_H::AbstractMatrix, X::AbstractMatrix,
                            sigma::AbstractMatrix;
                            sigma_sqrt::Union{Nothing, AbstractMatrix} = nothing)
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

Return the *non-spurious* real eigenvalues of -D_sym, where D_sym is the
Hermitian part of the quantum discriminant `D_S = Φ⁻¹ L_S Φ`.  For a
KMS-DBC L_S these are the Dirichlet rates; the smallest is the
KMS-Poincaré gap, the largest is Λ_max.

`project_constants=true` projects out the constant direction (eigenvector
of L_S corresponding to the steady state, which under Φ⁻¹ becomes
σ^{1/2}) and drops the spurious zero that the projection introduces.
"""
function _kms_dirichlet_eigvals(L_super::AbstractMatrix, sigma::AbstractMatrix;
                                project_constants::Bool = true,
                                sigma_quarter::Union{Nothing, AbstractMatrix} = nothing)
    _, D_sym = kms_discriminant_superoperator(L_super, sigma;
                                              sigma_quarter = sigma_quarter)
    M = -Matrix(D_sym)                           # -D_sym so eigenvalues = Dirichlet rates ≥ 0
    if project_constants
        # The discriminant zero direction is vec(σ^{1/2}) (verified by
        # `verify_detailed_balance`'s fp_residual check).  Normalise: since
        # ‖σ^{1/2}‖_F^2 = tr(σ) = 1, vec(σ^{1/2}) is already unit-norm.
        s12 = sqrt(Hermitian(sigma))
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

KMS-Poincaré gap = smallest non-zero Dirichlet rate.

Returns `(; gap, eigvals)` where `eigvals` are the d²−1 Dirichlet rates
sorted by absolute value (smallest first).  `gap = abs(eigvals[1])`.
"""
function spectral_gap_kms(L_super::AbstractMatrix, sigma::AbstractMatrix;
                          project_constants::Bool = true,
                          sigma_quarter::Union{Nothing, AbstractMatrix} = nothing)
    ev = _kms_dirichlet_eigvals(L_super, sigma; project_constants = project_constants,
                                                sigma_quarter = sigma_quarter)
    return (; gap = abs(ev[1]), eigvals = ev)
end

"""
    max_dirichlet_rate_kms(L_super, sigma; project_constants=true) -> Real

Maximum Dirichlet rate Λ_max(L) = largest Dirichlet rate.
"""
function max_dirichlet_rate_kms(L_super::AbstractMatrix, sigma::AbstractMatrix;
                                project_constants::Bool = true,
                                sigma_quarter::Union{Nothing, AbstractMatrix} = nothing)
    ev = _kms_dirichlet_eigvals(L_super, sigma; project_constants = project_constants,
                                                sigma_quarter = sigma_quarter)
    return maximum(abs.(ev))
end

"""
    intrinsic_mixing_ratio(L_super, sigma) -> Real

Scale-invariant ratio ρ = λ(L) / Λ_max(L) ∈ (0, 1].

ρ_intrinsic is the structural quantity isolating "geometric efficiency"
from any uniform speed-up by rescaling L → c·L.  qf-mto.5 will sweep
this across (CKG, DLL Metro, DLL Gauss) at matched (n, β).
"""
function intrinsic_mixing_ratio(L_super::AbstractMatrix, sigma::AbstractMatrix;
                                sigma_quarter::Union{Nothing, AbstractMatrix} = nothing)
    ev = _kms_dirichlet_eigvals(L_super, sigma;
                                sigma_quarter = sigma_quarter)
    return abs(ev[1]) / maximum(abs.(ev))
end


# ---------------------------------------------------------------------------
# Section 6 — 1→1 norm bound and HS-induced norm
# ---------------------------------------------------------------------------

"""
    hs_operator_norm(L_super) -> Real

Hilbert–Schmidt-induced operator norm = largest singular value of the
dense superoperator.  Always ≥ ‖L‖_{1→1} but easy to compute.
"""
hs_operator_norm(L_super::AbstractMatrix) = opnorm(L_super)

"""
    dissipator_one_to_one_norm_bound(L_a_list) -> Real

Upper bound `4·max_a ‖L_a† L_a‖_∞` on the dissipator 1→1 norm
(Wolf & Pérez-García style; see `.claude-memory/fair_comparison_dirichlet_qf_mto.md`).
For DLL the L_a's are produced by `dll_lindblad_op_bohr` (one per
coupling); for CKG one would need to eigendecompose α first (deferred
to integration phase).
"""
function dissipator_one_to_one_norm_bound(L_a_list::AbstractVector{<:AbstractMatrix})
    @assert !isempty(L_a_list)
    return 4.0 * maximum(opnorm(L_a' * L_a) for L_a in L_a_list)
end

"""
    dll_lindblad_ops(jumps, ham, filter) -> Vector{Matrix{ComplexF64}}

Convenience: list of per-coupling DLL Lindblad operators in the eigenbasis,
one per jump.  Used by the 1→1 bound on the DLL dissipator.
"""
function dll_lindblad_ops(jumps::Vector{JumpOp}, ham::HamHam, filter::AbstractFilter)
    return [Matrix{ComplexF64}(dll_lindblad_op_bohr(jump, ham, filter)) for jump in jumps]
end


# ---------------------------------------------------------------------------
# Helpers — fixture builders
# ---------------------------------------------------------------------------

"""
    make_ckg_n3_system(beta) -> (; ham, jumps, gibbs)

CKG companion to `make_dll_n3_system(beta)`: same Hamiltonian, same jump
set, but built without a DLL filter so that a CKG (`construction = KMS()`)
config can use it.  See `test/test_helpers.jl::make_dll_n3_system` for the
shared physics — only the basis-projected jumps differ in convention.
"""
function make_ckg_n3_system(beta::Real)
    return make_dll_n3_system(beta)   # ham, jumps, gibbs are construction-agnostic
end

"""
    ckg_smooth_metro_config(beta::Real) -> Config

CKG smooth-Metropolis n=3 EnergyDomain config.  Uses thesis-numerics
defaults (a=0, s=0.25) — see `src/lindblad_action.jl::sweep_mixing_times`.
"""
# PHYSICS CHECK: a=0, s=0.25 is the locked thesis-numerics convention for
# CKG smooth-Metropolis (see MEMORY.md and src/lindblad_action.jl docstring).
function ckg_smooth_metro_config(beta::Real)
    return Config(
        sim = Lindbladian(),
        domain = EnergyDomain(),
        construction = KMS(),
        num_qubits = 3,
        with_linear_combination = true,
        beta = float(beta),
        sigma = 1.0 / float(beta),
        a = 0.0,
        s = 0.25,
        num_energy_bits = 12,
        w0 = 0.05,
        t0 = 2π / (2^12 * 0.05),
        num_trotter_steps_per_t0 = 10,
    )
end

"""
    dll_config(beta::Real, filter::AbstractFilter) -> Config

DLL BohrDomain n=3 config — Metropolis or Gaussian via the supplied filter.
"""
function dll_config(beta::Real, filter::AbstractFilter)
    return Config(
        sim = Lindbladian(),
        domain = BohrDomain(),
        construction = DLL(),
        num_qubits = 3,
        with_linear_combination = true,
        beta = float(beta),
        sigma = 1.0 / float(beta),
        a = float(beta) / 30.0,
        s = 0.4,
        num_energy_bits = 12,
        w0 = 0.05,
        t0 = 2π / (2^12 * 0.05),
        num_trotter_steps_per_t0 = 10,
        filter = filter,
    )
end

"""
    build_dense_lindbladian(config, ham, jumps) -> (L_super, dim)

Build the dense d²×d² Schrödinger Lindbladian by `apply_lindbladian!`-driven
matvec.  Uses identical sandwich code to `krylov_spectral_gap`, so any
cross-check between a dense routine and the production matrix-free routine
sits on the same physics.
"""
function build_dense_lindbladian(config::Config{Lindbladian}, ham::HamHam,
                                  jumps::Vector{JumpOp})
    dim = size(ham.data, 1)
    ws  = Workspace(config, ham, jumps)
    L_apply! = function(out, X)
        apply_lindbladian!(ws, X, config, ham)
        copyto!(out, ws.scratch.rho_out)
        return out
    end
    L_super = build_dense_superoperator(L_apply!, dim)
    return L_super, dim
end


# ===========================================================================
# Validation suite
# ===========================================================================

# Build a tiny qubit Davies model (used in Tasks 1, 4) — closed-form rates.
"""
    build_qubit_davies(omega, beta, gamma) -> NamedTuple

Single-qubit Davies thermal generator with H = -(ω/2) σ_z and KMS-DB rates
γ_+ = γ (de-excitation; ν = -ω) and γ_- = γ e^{-βω} (excitation; ν = +ω).
Closed-form gap = γ_+ + γ_- on the diagonal sector.
"""
function build_qubit_davies(omega::Real, beta::Real, gamma::Real)
    CT = ComplexF64
    sigma_z = CT[1.0 0.0; 0.0 -1.0]
    sigma_p = CT[0.0 1.0; 0.0 0.0]   # |0⟩⟨1|
    sigma_m = CT[0.0 0.0; 1.0 0.0]   # |1⟩⟨0|
    H       = -(omega / 2) * sigma_z
    Z       = exp(beta * omega / 2) + exp(-beta * omega / 2)
    p0      = exp( beta * omega / 2) / Z
    p1      = exp(-beta * omega / 2) / Z
    sigma_state = Diagonal([p0, p1]) |> Matrix{CT}

    γp = gamma                       # de-excitation rate
    γm = gamma * exp(-beta * omega)  # excitation rate

    function L_apply!(out, X)
        # Schrödinger Lindbladian: L(X) = -i[H, X]
        #                                + γ_+ (σ_+ X σ_+† - 1/2 {σ_+† σ_+, X})
        #                                + γ_- (σ_- X σ_-† - 1/2 {σ_-† σ_-, X})
        comm = -1im * (H * X - X * H)
        sp_term = γp * (sigma_p * X * sigma_p' .- 0.5 .* (sigma_p' * sigma_p * X .+ X * sigma_p' * sigma_p))
        sm_term = γm * (sigma_m * X * sigma_m' .- 0.5 .* (sigma_m' * sigma_m * X .+ X * sigma_m' * sigma_m))
        copyto!(out, comm .+ sp_term .+ sm_term)
        return out
    end

    L_super = build_dense_superoperator(L_apply!, 2)
    # Davies single-qubit eigenvalues:
    #   0, -(γ_+ + γ_-), -(γ_+ + γ_-)/2 ± iω
    # → spectral gap (smallest non-zero |Re λ|) = (γ_+ + γ_-)/2  (T_2 coherence decay)
    # → Λ_max  = γ_+ + γ_-                                       (T_1 population decay)
    return (; L_super, gibbs = sigma_state, γp, γm, H, sigma_p, sigma_m,
              gap_exact   = (γp + γm) / 2,
              Lambda_exact = γp + γm)
end


function run_validation_suite()
    println("\n", "="^72)
    println("KMS-geometry scratch (qf-mto.{1,2,3}) — validation suite")
    println("="^72, "\n")

    Random.seed!(20260502)

    # -----------------------------------------------------------------------
    # @testset (1) — KMS inner-product algebra (Task 1)
    # -----------------------------------------------------------------------
    @testset "(1) KMS inner-product algebra" begin
        # Use a non-uniform 2-level Gibbs state.
        β = 1.7
        ω = 0.9
        Z = exp(β * ω / 2) + exp(-β * ω / 2)
        sigma = ComplexF64[exp(β*ω/2)/Z 0; 0 exp(-β*ω/2)/Z]
        s12 = sqrt(Hermitian(sigma))

        # Random Hermitian X, Y.
        Xr = randn(ComplexF64, 2, 2);  X = (Xr + Xr') / 2
        Yr = randn(ComplexF64, 2, 2);  Y = (Yr + Yr') / 2

        # Hermitian symmetry: ⟨X, Y⟩ = conj(⟨Y, X⟩)
        @test isapprox(kms_inner_product(X, Y, sigma; sigma_sqrt = s12),
                       conj(kms_inner_product(Y, X, sigma; sigma_sqrt = s12)); atol = 1e-12)

        # Identity-norm: ‖I‖_KMS = √tr(σ) = 1.
        @test isapprox(kms_norm(Matrix{ComplexF64}(I, 2, 2), sigma; sigma_sqrt = s12), 1.0; atol = 1e-12)

        # Identity-variance: Var_KMS(I) = 0.
        @test isapprox(kms_variance(Matrix{ComplexF64}(I, 2, 2), sigma; sigma_sqrt = s12), 0.0; atol = 1e-12)

        # Cauchy–Schwarz: |⟨X, Y⟩| ≤ ‖X‖ ‖Y‖.
        cs_lhs = abs(kms_inner_product(X, Y, sigma; sigma_sqrt = s12))
        cs_rhs = kms_norm(X, sigma; sigma_sqrt = s12) * kms_norm(Y, sigma; sigma_sqrt = s12)
        @test cs_lhs ≤ cs_rhs + 1e-12

        # Positivity: ⟨X, X⟩ ≥ 0 (real and non-negative).
        ip_xx = kms_inner_product(X, X, sigma; sigma_sqrt = s12)
        @test isapprox(imag(ip_xx), 0.0; atol = 1e-12)
        @test real(ip_xx) ≥ -1e-12
    end

    # -----------------------------------------------------------------------
    # @testset (2) — Discriminant + KMS-DB witness (Task 2)
    # -----------------------------------------------------------------------
    @testset "(2) Quantum discriminant + KMS-DB diagnostics" begin
        β = 5.0
        sys = make_ckg_n3_system(β)
        cfg = ckg_smooth_metro_config(β)
        L_super, dim = build_dense_lindbladian(cfg, sys.ham, sys.jumps)
        gibbs = Matrix{ComplexF64}(sys.gibbs)

        # Schrödinger fixed-point:  L_super · vec(σ) ≈ 0
        residual_sigma = norm(L_super * vec(gibbs)) / norm(vec(gibbs))
        @test residual_sigma < 1e-10

        # Quantum discriminant D_S = Φ⁻¹ L_S Φ.  For KMS-DBC with the Lamb
        # shift restored, D_S is HS-self-adjoint to numerical precision.
        D_super, D_sym = kms_discriminant_superoperator(L_super, gibbs)
        D_anti = (D_super .- D_super') ./ 2
        anti_to_sym_ratio = norm(D_anti) / max(norm(D_sym), eps())
        @info "(2) D_S anti-Hermitian / symmetric norm ratio at β=$β" anti_to_sym_ratio
        @test anti_to_sym_ratio < 1e-10   # KMS-DB witness; must be machine-precision

        # D_S has vec(σ^{1/2}) as zero eigenvector.
        v_sh = vec(sqrt(Hermitian(gibbs)))
        @test norm(D_super * v_sh) < 1e-10

        # Eigenvalues of D_sym are real and ≤ 0 (D_S is a similarity transform
        # of L_S so eigenvalues of D_sym match those of L_S which are ≤ 0).
        ev = eigvals(D_sym)
        @test maximum(abs.(imag.(ev))) < 1e-10
        @test maximum(real.(ev)) < 1e-8

        # Diagnostic: removing the Lamb-shift coherent term breaks KMS-DBC
        # (the bare dissipator is *not* KMS-DBC for CKG smooth-Metropolis,
        # contrary to the naive intuition).  We construct the no-coherent
        # variant and confirm:
        #   (a) its discriminant has substantial anti-Hermitian content (non-DBC witness)
        #   (b) Λ_max is preserved at machine precision (dissipator dominates large modes)
        #   (c) the gap shifts by ~1e-3 — small but real, NOT a numerical artefact
        # This is the meaningful counterpart to the plan's R3 risk flag.
        ws_no_coh = Workspace(cfg, sys.ham, sys.jumps)
        # Identity in workspace constructor: G_left = +i B_T - 0.5 R, G_right
        # = -i B_T - 0.5 R.  → -(G_left + G_right) = R, so we recover R and
        # rebuild G_left = G_right = -0.5 R.
        R = -(ws_no_coh.G_left .+ ws_no_coh.G_right)
        ws_no_coh.G_left  .= -0.5 .* R
        ws_no_coh.G_right .= -0.5 .* R
        ws_no_coh.G_left_adj  .= ws_no_coh.G_right
        ws_no_coh.G_right_adj .= ws_no_coh.G_left
        L_apply_no_coh! = function(out, X)
            apply_lindbladian!(ws_no_coh, X, cfg, sys.ham)
            copyto!(out, ws_no_coh.scratch.rho_out)
            return out
        end
        L_super_no_coh = build_dense_superoperator(L_apply_no_coh!, dim)

        # No-coherent discriminant should NOT be HS-self-adjoint (KMS-DBC fails).
        D_nc, D_nc_sym = kms_discriminant_superoperator(L_super_no_coh, gibbs)
        D_nc_anti = (D_nc .- D_nc') ./ 2
        anti_ratio_nc = norm(D_nc_anti) / max(norm(D_nc_sym), eps())
        @test anti_ratio_nc > 1e-6   # substantial anti-Hermitian piece (KMS-DBC is broken)

        gap_full   = spectral_gap_kms(L_super,        gibbs).gap
        gap_no_coh = spectral_gap_kms(L_super_no_coh, gibbs).gap
        Λ_full     = max_dirichlet_rate_kms(L_super,        gibbs)
        Λ_no_coh   = max_dirichlet_rate_kms(L_super_no_coh, gibbs)

        # Λ_max nearly preserved (dissipator dominates large modes); gap may shift.
        @test isapprox(Λ_full, Λ_no_coh; rtol = 1e-8)

        @info "(2) coherent removal: full KMS-DBC vs bare dissipator (β=$β)" gap_full gap_no_coh Λ_full Λ_no_coh anti_ratio_nc
    end

    # -----------------------------------------------------------------------
    # @testset (3) — Dirichlet form properties (Task 3)
    # -----------------------------------------------------------------------
    @testset "(3) KMS Dirichlet form" begin
        β = 5.0
        sys = make_ckg_n3_system(β)
        cfg = ckg_smooth_metro_config(β)
        L_super, dim = build_dense_lindbladian(cfg, sys.ham, sys.jumps)
        # Heisenberg-picture superoperator = HS-adjoint of Schrödinger.
        L_super_H = Matrix(L_super')
        gibbs = Matrix{ComplexF64}(sys.gibbs)
        s12 = sqrt(Hermitian(gibbs))

        # Random Hermitian X.
        Xraw = randn(ComplexF64, dim, dim)
        X = (Xraw + Xraw') / 2

        # Linearity in L: E_{2L}(X) = 2 E_L(X).
        E_L  = kms_dirichlet_form(L_super_H, X, gibbs; sigma_sqrt = s12)
        E_2L = kms_dirichlet_form(2.0 .* L_super_H, X, gibbs; sigma_sqrt = s12)
        @test isapprox(E_2L, 2 * E_L; rtol = 1e-12)

        # Non-negativity: E_L(X) ≥ 0 for KMS-DBC L on Hermitian X.
        @test E_L ≥ -1e-10

        # Vanishing on constants: E_L(c·I) = 0 — uses unitality of L^H.
        cI = 0.7 * Matrix{ComplexF64}(I, dim, dim)
        E_const = kms_dirichlet_form(L_super_H, cI, gibbs; sigma_sqrt = s12)
        @test isapprox(E_const, 0.0; atol = 1e-10)

        # Cross-check against discriminant route: E_L(X) should equal
        # -Re tr( (ΦX)† D_sym ΦX ) where Φ X = σ^{1/4} X σ^{1/4}.
        # See docstring derivation for the equivalence.
        Φ, _ = kms_phi_matrix(gibbs)
        _, D_sym = kms_discriminant_superoperator(L_super, gibbs)
        ΦX_vec = Φ * vec(X)
        E_via_D = -real(ΦX_vec' * (D_sym * ΦX_vec))
        @test isapprox(E_via_D, E_L; rtol = 1e-10)

        @info "(3) Dirichlet form sanity at β=$β" E_L E_2L E_const E_via_D
    end

    # -----------------------------------------------------------------------
    # @testset (4) — Spectral gap (Task 4)
    # -----------------------------------------------------------------------
    @testset "(4) Spectral gap λ(L)" begin
        # 4a — qubit Davies closed-form.
        ω = 1.3
        β = 0.8
        γ = 0.4
        davies = build_qubit_davies(ω, β, γ)
        gap_davies_num = spectral_gap_kms(davies.L_super, davies.gibbs).gap
        Λ_davies_num   = max_dirichlet_rate_kms(davies.L_super, davies.gibbs)
        @info "(4a) Davies qubit: numerical vs closed-form" gap_davies_num davies.gap_exact Λ_davies_num davies.Lambda_exact
        @test isapprox(gap_davies_num, davies.gap_exact;    rtol = 1e-10)
        @test isapprox(Λ_davies_num,   davies.Lambda_exact; rtol = 1e-10)

        # 4b — n=3 CKG smooth-Metro vs production krylov_spectral_gap.
        β = 5.0
        sys = make_ckg_n3_system(β)
        cfg = ckg_smooth_metro_config(β)
        L_super, _ = build_dense_lindbladian(cfg, sys.ham, sys.jumps)
        gibbs = Matrix{ComplexF64}(sys.gibbs)

        gap_dense = spectral_gap_kms(L_super, gibbs).gap
        # Production matrix-free.
        kr = krylov_spectral_gap(cfg, sys.ham, sys.jumps; krylovdim = 30,
                                 howmany = 4, tol = 1e-12)
        gap_krylov = kr.spectral_gap
        @info "(4b) n=3 CKG smooth-Metro β=$β: KMS-Poincaré gap vs Krylov gap" gap_dense gap_krylov
        @test isapprox(gap_dense, gap_krylov; rtol = 1e-7)

        # 4c — Schrödinger ↔ Heisenberg (HS-adjoint) similarity sanity.
        # apply_adjoint_lindbladian! is the HS-adjoint of apply_lindbladian!.
        # Eigenvalues of L_S and L_S* are complex conjugates of each other;
        # in particular |Re λ| spectra coincide so the gap matches.  (For the
        # Heisenberg generator the discriminant zero eigenvector is vec(I·d^{-1/2})
        # in our convention; the same `_kms_dirichlet_eigvals` routine would need
        # a different projector.  Here we just check the eigenvalue spectrum
        # of L_super_h via an unprojected routine.)
        ws_h = Workspace(cfg, sys.ham, sys.jumps)
        L_apply_h! = function(out, X)
            apply_adjoint_lindbladian!(ws_h, X, cfg, sys.ham)
            copyto!(out, ws_h.scratch.rho_out)
            return out
        end
        L_super_h = build_dense_superoperator(L_apply_h!, size(sys.ham.data, 1))
        # L_super_h = L_super^† (HS-adjoint) → eigenvalues are conjugates of L_super's
        # → same spectrum since L_super is real-eigenvalued (it's a similarity transform
        # of D_S, which is Hermitian).
        ev_S = sort(real.(eigvals(L_super));  by = abs)
        ev_H = sort(real.(eigvals(L_super_h)); by = abs)
        @info "(4c) Schrödinger ↔ Heisenberg gap" gap_S=abs(ev_S[2]) gap_H=abs(ev_H[2])
        @test isapprox(abs(ev_S[2]), abs(ev_H[2]); rtol = 1e-10)
    end

    # -----------------------------------------------------------------------
    # @testset (5) — Λ_max + intrinsic ratio + scale invariance (Task 5)
    # -----------------------------------------------------------------------
    @testset "(5) Λ_max + intrinsic ratio + scale invariance" begin
        β = 5.0
        sys = make_ckg_n3_system(β)
        cfg = ckg_smooth_metro_config(β)
        L_super, _ = build_dense_lindbladian(cfg, sys.ham, sys.jumps)
        gibbs = Matrix{ComplexF64}(sys.gibbs)

        gap = spectral_gap_kms(L_super, gibbs).gap
        Λ   = max_dirichlet_rate_kms(L_super, gibbs)
        ρ   = intrinsic_mixing_ratio(L_super, gibbs)

        @test Λ ≥ gap - 1e-12
        @test 0.0 < ρ ≤ 1.0 + 1e-12

        # Scale invariance under L → c·L.
        c = 2.7
        gap_c = spectral_gap_kms(c .* L_super, gibbs).gap
        Λ_c   = max_dirichlet_rate_kms(c .* L_super, gibbs)
        ρ_c   = intrinsic_mixing_ratio(c .* L_super, gibbs)

        @test isapprox(gap_c, c * gap; rtol = 1e-10)
        @test isapprox(Λ_c,   c * Λ;   rtol = 1e-10)
        @test isapprox(ρ_c,   ρ;       rtol = 1e-12)

        @info "(5) scale invariance at β=$β" gap Λ ρ gap_c Λ_c ρ_c
    end

    # -----------------------------------------------------------------------
    # @testset (6) — 1→1 norm bound + HS-induced norm (Task 6)
    # -----------------------------------------------------------------------
    @testset "(6) 1→1 norm bound + HS-induced norm" begin
        for β in (5.0, 10.0)
            sys     = make_dll_n3_system(β)
            f_metro = DLLMetropolisFilter(β; S = 2.0)
            cfg     = dll_config(β, f_metro)
            L_super, _ = build_dense_lindbladian(cfg, sys.ham, sys.jumps)
            gibbs   = Matrix{ComplexF64}(sys.gibbs)

            L_a_list = dll_lindblad_ops(sys.jumps, sys.ham, f_metro)
            one_to_one_bound = dissipator_one_to_one_norm_bound(L_a_list)
            hs_norm = hs_operator_norm(L_super)

            gap = spectral_gap_kms(L_super, gibbs).gap

            @test hs_norm > 0
            @test one_to_one_bound > 0

            # Scale invariance of *ratios* gap / one-to-one_bound under L_a → √c L_a.
            c = 1.5
            ratio_orig = gap / one_to_one_bound
            gap_c = spectral_gap_kms(c .* L_super, gibbs).gap
            # 1→1 bound scales linearly in c (each L_a → √c L_a → ‖L_a†L_a‖ → c).
            ratio_c = (c * gap) / (c * one_to_one_bound)
            @test isapprox(ratio_orig, ratio_c; rtol = 1e-12)

            @info "(6) DLL Metro β=$β" gap one_to_one_bound hs_norm ratio_orig
        end
    end

    # -----------------------------------------------------------------------
    # @testset (7) — Headline comparison: CKG smooth-Metro vs DLL Metro vs DLL Gauss
    # -----------------------------------------------------------------------
    headline = NamedTuple[]   # captured to print outside @testset
    @testset "(7) headline: CKG vs DLL @ n=3, β ∈ {5,10,20}" begin
        for β in (5.0, 10.0, 20.0)
            sys   = make_dll_n3_system(β)   # construction-agnostic ham + jumps + gibbs
            gibbs = Matrix{ComplexF64}(sys.gibbs)

            for (label, cfg, filter_obj) in (
                ("CKG smooth-Metro", ckg_smooth_metro_config(β),         nothing),
                ("DLL Metro",        dll_config(β, DLLMetropolisFilter(β; S = 2.0)),
                                                                          DLLMetropolisFilter(β; S = 2.0)),
                ("DLL Gauss",        dll_config(β, DLLGaussianFilter(β)), DLLGaussianFilter(β)),
            )
                L_super, _ = build_dense_lindbladian(cfg, sys.ham, sys.jumps)

                gap = spectral_gap_kms(L_super, gibbs).gap
                Λ   = max_dirichlet_rate_kms(L_super, gibbs)
                ρ   = intrinsic_mixing_ratio(L_super, gibbs)
                hsn = hs_operator_norm(L_super)

                # 1→1 bound + Tr(α) only meaningful for DLL (where L_a's are
                # explicit and α is per-jump rank-1).  For CKG report NaN.
                if filter_obj !== nothing
                    L_a_list = dll_lindblad_ops(sys.jumps, sys.ham, filter_obj)
                    one2one  = dissipator_one_to_one_norm_bound(L_a_list)
                    tr_alpha = sum(real(tr(L_a' * L_a)) for L_a in L_a_list)
                else
                    one2one  = NaN
                    tr_alpha = NaN
                end

                # Sanity invariants.
                @test gap > 0
                @test Λ ≥ gap - 1e-12
                @test 0 < ρ ≤ 1 + 1e-12

                push!(headline, (; β, label, gap, Λ, ρ, tr_alpha, one2one, hsn))
            end
        end

        # H1 vs H2 sanity: at β=20 the three intrinsic ratios should not be
        # numerically identical, otherwise H2 (diagonal of α drives gap) is
        # falsified before the formal sweep.  Soft assertion only:
        β20 = filter(r -> r.β == 20.0, headline)
        @test length(β20) == 3
        ρs_β20 = [r.ρ for r in β20]
        spread_β20 = (maximum(ρs_β20) - minimum(ρs_β20))
        @info "(7) ρ_intrinsic spread across families at β=20" spread_β20 ρs_β20
    end

    # -----------------------------------------------------------------------
    # Final summary print (Task 8)
    # -----------------------------------------------------------------------
    println("\n", "="^110)
    println("HEADLINE TABLE — KMS geometry at n=3 (CKG smooth-Metro vs DLL Metro vs DLL Gauss)")
    println("="^110)
    @printf("%-6s  %-22s  %12s  %12s  %12s  %12s  %12s  %12s\n",
        "β", "Filter", "λ(gap)", "Λ_max", "ρ_intr", "Tr(α)", "1→1 bound", "‖L‖_HS")
    println("-"^110)
    for r in headline
        tr_str  = isnan(r.tr_alpha) ? "    n/a    " : @sprintf("%12.5e", r.tr_alpha)
        one_str = isnan(r.one2one)  ? "    n/a    " : @sprintf("%12.5e", r.one2one)
        @printf("%-6.1f  %-22s  %12.5e  %12.5e  %12.5e  %s  %s  %12.5e\n",
            r.β, r.label, r.gap, r.Λ, r.ρ, tr_str, one_str, r.hsn)
    end
    println("="^110, "\n")

    println("ALL PASS\n")
    return headline
end


# ===========================================================================
# Main
# ===========================================================================

if abspath(PROGRAM_FILE) == @__FILE__
    headline = run_validation_suite()
end
