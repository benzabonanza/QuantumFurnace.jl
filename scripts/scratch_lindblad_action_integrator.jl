# scratch_lindblad_action_integrator.jl
#
# Matrix-free Krylov-subspace integrators for two equivalent flows on a
# finite-dimensional state space:
#
#   L-mode :  d/dt rho(t)  = L(rho(t))                    (Lindbladian on density operator)
#   K-mode :  d/dt psi(t)  = K(psi(t)),  K = D(sigma, L_D) (KMS quantum discriminant action)
#
# In Heisenberg-picture (column-stacking vec convention),
# rho_vec(t) = exp(t L_super) rho_vec(0) and  psi(t) = sigma^{-1/4} rho(t) sigma^{-1/4},
# so the two flows are exact reparametrisations of each other on the diagonal
# (Bohr-zero) sector of operator space, where the H commutator vanishes.
#
# Equation references:
# - eq:discriminant     :  D(X) = sigma^{-1/4} L(sigma^{1/4} X sigma^{1/4}) sigma^{-1/4}
#                         (thesis Eq. eq:discriminant; Chen-Kastoryano-Gilyen 2023; Chen 2025 Sec. III)
# - chi^2 distance      :  chi^2(rho, sigma) = || sigma^{-1/4}(rho - sigma) sigma^{-1/4} ||_F^2
# - Davies generator    :  L(rho) = -i[H, rho]
#                                 + sum_a gamma_a [ J_a rho J_a^dag
#                                                   - 1/2 { J_a^dag J_a, rho } ]
#                         with KMS-DB ratio  gamma(-omega) / gamma(omega) = e^{-beta omega}.
# - K HS-self-adjoint   :  Under KMS-DB, the *dissipative* part L_D yields
#                         K_D = D(sigma, L_D) self-adjoint w.r.t.
#                         <X, Y>_HS = tr(X^dag Y).  The bare-H commutator
#                         contributes an HS-anti-Hermitian piece to the full
#                         discriminant; with a Lamb-shift restoration of the
#                         coherent term the full L is KMS-DB and the full D is
#                         HS-self-adjoint (Chen 2025; QuantumFurnace
#                         src/discriminant.jl).  This prototype uses bare H,
#                         so the L-mode trajectory uses the full L = -i[H,.] + L_D
#                         while the K-mode trajectory builds K from the
#                         dissipator L_D only.  For the diagonal initial state
#                         the two trajectories coincide because [H, rho_diag] = 0.
#
# Numerics:
#   KrylovKit v0.10.2   `exponentiate(A, t, v; ishermitian, krylovdim, tol)` returns (y, info)
#                       with `info.numops` matvec count and `info.converged in {0, 1}`.
#                       `ishermitian=true` short-circuits to Lanczos; `false` to Arnoldi.
#
# Run as:    julia --project scripts/scratch_lindblad_action_integrator.jl

using KrylovKit
using LinearAlgebra
using Printf


# -----------------------------------------------------------------------------
# Single-qubit Davies thermal Lindbladian — analytic-rate toy
# -----------------------------------------------------------------------------

"""
    build_davies_2x2(omega::T, beta::T, gamma::T) where {T<:Real}

Construct the single-qubit Davies thermal generator with Hamiltonian
`H = -(omega/2) sigma_z` (so |0> is the energetic ground state) and
KMS-DB jump rates  `gamma_+ = gamma`  for de-excitation (jump `sigma_+ = |0><1|`)
and  `gamma_- = gamma * exp(-beta * omega)`  for excitation (jump
`sigma_- = |1><0|`).

Returns a NamedTuple containing:
- `L_super::Matrix{Complex{T}}`        :  4x4 full Lindblad superoperator H + L_D (column-stacking vec)
- `L_super_D::Matrix{Complex{T}}`      :  4x4 dissipator-only superoperator L_D (used to define K)
- `L_apply!`                           :  in-place closure for the FULL L (used by L-mode)
- `K_apply!`                           :  in-place closure for K = D(sigma, L_D)  (used by K-mode)
- `sigma_beta::Matrix{Complex{T}}`     :  Gibbs state (diagonal in computational basis)
- `sigma_quarter::Vector{T}`           :  diag(sigma^{1/4})
- `sigma_inv_quarter::Vector{T}`       :  diag(sigma^{-1/4})
- `sigma_half::Vector{T}`              :  diag(sigma^{1/2})
- `gamma_plus::T`, `gamma_minus::T`    :  the two KMS-DB rates
- `decay_rate::T`                      :  gamma_+ + gamma_- (diagonal/T_1 relaxation rate)
"""
function build_davies_2x2(omega::T, beta::T, gamma::T) where {T<:Real}
    CT = Complex{T}

    # Hamiltonian and Pauli operators (CT-typed)
    sigma_z = CT[T(1)  T(0); T(0)  T(-1)]
    sigma_p = CT[T(0)  T(1); T(0)   T(0)]   # |0><1|
    sigma_m = CT[T(0)  T(0); T(1)   T(0)]   # |1><0|
    Id2     = Matrix{CT}(I, 2, 2)

    # H = -(omega/2) sigma_z  =>  E_0 = -omega/2 (ground), E_1 = +omega/2 (excited)
    H       = -(omega / 2) * sigma_z

    # Gibbs state (Boltzmann weights). Ground state |0> dominant.
    Z       = exp( beta * omega / 2) + exp(-beta * omega / 2)
    p0      = exp( beta * omega / 2) / Z
    p1      = exp(-beta * omega / 2) / Z
    sigma_beta        = Matrix{CT}(undef, 2, 2)
    sigma_beta       .= 0
    sigma_beta[1, 1]  = p0
    sigma_beta[2, 2]  = p1
    sigma_quarter     = T[p0^T(0.25), p1^T(0.25)]
    sigma_inv_quarter = T[p0^T(-0.25), p1^T(-0.25)]
    sigma_half        = T[p0^T(0.5),  p1^T(0.5)]

    # KMS-DB jumps. With H = -(omega/2) sigma_z and |0> ground:
    #   sigma_+ = |0><1|  takes |1> -> |0|: Bohr frequency E_0 - E_1 = -omega.
    #     "de-excitation" jump, rate gamma_+ = gamma.
    #   sigma_- = |1><0|  takes |0> -> |1|: Bohr frequency E_1 - E_0 = +omega.
    #     "excitation" jump, rate gamma_- = gamma * e^{-beta*omega} (KMS).
    gamma_plus  = gamma
    gamma_minus = gamma * exp(-beta * omega)
    decay_rate  = gamma_plus + gamma_minus

    # ---- Build dense superoperators (column-stacking vec convention). ------
    # For column-stacking, vec(A B C) = (C^T \otimes A) vec(B). Hence:
    #   vec([H, rho]) = (I \otimes H - H^T \otimes I) vec(rho)
    #   vec(J rho J^dag) = (conj(J) \otimes J) vec(rho)
    #   vec({J^dag J, rho}/2) = 1/2 (I \otimes J^dag J + (J^dag J)^T \otimes I) vec(rho)
    H_super  = -1im * (kron(Id2, H) - kron(transpose(H), Id2))

    # De-excitation jump: J = sigma_+ at rate gamma_+.  J^dag J = sigma_- sigma_+ = |1><1|.
    JdJp = sigma_m * sigma_p   # = |1><1|
    D_plus = gamma_plus * (
                  kron(conj(sigma_p), sigma_p)
                - T(0.5) * (kron(Id2, JdJp) + kron(transpose(JdJp), Id2))
              )
    # Excitation jump: J = sigma_- at rate gamma_-.  J^dag J = sigma_+ sigma_- = |0><0|.
    JdJm = sigma_p * sigma_m   # = |0><0|
    D_minus = gamma_minus * (
                  kron(conj(sigma_m), sigma_m)
                - T(0.5) * (kron(Id2, JdJm) + kron(transpose(JdJm), Id2))
              )

    L_super_D = Matrix{CT}(D_plus + D_minus)             # dissipator only
    L_super   = Matrix{CT}(H_super + L_super_D)          # full L = -i[H,.] + L_D

    # ---- L-apply closure (matrix-free, FULL L). ----------------------------
    # Using the dense 4x4 L_super here is fine for a 2x2 toy; the integrator's
    # matrix-free property is about not forming a (d^2 x d^2) operator inside
    # the Krylov body for arbitrary d, which is guaranteed by feeding this
    # closure to `exponentiate` (it never sees the matrix).
    function L_apply!(out::AbstractMatrix, x::AbstractMatrix)
        mul!(vec(out), L_super, vec(x))
        return out
    end

    # ---- L_D-apply closure (dissipator only) used inside K_apply!. ---------
    function L_D_apply!(out::AbstractMatrix, x::AbstractMatrix)
        mul!(vec(out), L_super_D, vec(x))
        return out
    end

    # ---- K-apply closure: K(X) = D(sigma, L_D)(X)
    #   = sigma^{-1/4} L_D(sigma^{1/4} X sigma^{1/4}) sigma^{-1/4}.
    # Two persistent 2x2 buffers captured by the closure.
    work1 = Matrix{CT}(undef, 2, 2)
    work2 = Matrix{CT}(undef, 2, 2)
    function K_apply!(out::AbstractMatrix, x::AbstractMatrix)
        @inbounds for j in 1:2, i in 1:2
            work1[i, j] = sigma_quarter[i] * x[i, j] * sigma_quarter[j]
        end
        L_D_apply!(work2, work1)
        @inbounds for j in 1:2, i in 1:2
            out[i, j] = sigma_inv_quarter[i] * work2[i, j] * sigma_inv_quarter[j]
        end
        return out
    end

    return (
        L_super           = L_super,
        L_super_D         = L_super_D,
        L_apply!          = L_apply!,
        K_apply!          = K_apply!,
        sigma_beta        = sigma_beta,
        sigma_quarter     = sigma_quarter,
        sigma_inv_quarter = sigma_inv_quarter,
        sigma_half        = sigma_half,
        gamma_plus        = gamma_plus,
        gamma_minus       = gamma_minus,
        decay_rate        = decay_rate,
    )
end


# -----------------------------------------------------------------------------
# Construction-time sanity asserts (printed)
# -----------------------------------------------------------------------------

function _print_construction_checks(toy)
    println("\n--- Construction checks ---")
    p0, p1 = real(toy.sigma_beta[1, 1]), real(toy.sigma_beta[2, 2])
    @printf "  sigma_beta diag         : (%.6f, %.6f), sum = %.16f\n" p0 p1 (p0 + p1)
    @assert p0 > 0 && p1 > 0
    @assert isapprox(p0 + p1, 1.0; atol=1e-14)

    # KMS ratio: gamma_-/gamma_+ should equal e^{-beta omega}.
    # Derive beta*omega from p_0/p_1 = e^{beta omega}.
    bw         = log(p0 / p1)                # = beta * omega
    kms_ratio  = toy.gamma_minus / toy.gamma_plus
    kms_target = exp(-bw)
    @printf "  KMS-DB rate ratio       : gamma_-/gamma_+ = %.16f, target e^{-beta*omega} = %.16f\n" kms_ratio kms_target
    @assert isapprox(kms_ratio, kms_target; rtol=1e-14)

    # K HS-self-adjoint check.
    # In column-stacking, D_super = (S^{-1} \otimes S^{-1}) * L * (S \otimes S),
    # so the full discriminant uses M_left = kron(Si, Si), M_right = kron(S, S).
    # The bare-H commutator contributes an HS-anti-Hermitian piece, so the
    # full-L discriminant is NOT HS-self-adjoint; the dissipator-only
    # discriminant IS HS-self-adjoint (KMS-DB witness for L_D).
    Sq      = Diagonal(toy.sigma_quarter)
    Sqi     = Diagonal(toy.sigma_inv_quarter)
    M_left  = kron(Sqi, Sqi)
    M_right = kron(Sq,  Sq)
    D_super_full = M_left * toy.L_super   * M_right
    D_super_diss = M_left * toy.L_super_D * M_right
    skew_full = opnorm(D_super_full - D_super_full') / opnorm(D_super_full)
    skew_diss = opnorm(D_super_diss - D_super_diss') / opnorm(D_super_diss)
    @printf "  K HS-self-adjoint check (dissipator-only) : ||K_D - K_D'||/||K_D|| = %.3e\n" skew_diss
    @printf "  K HS-self-adjoint check (full L = H + L_D): ||K - K'||/||K||       = %.3e   (expected non-zero: bare-H breaks DB)\n" skew_full
    @assert skew_diss < 1e-12

    # L_super spectrum: 0 (steady state), -(g+ + g-) (T_1, real),
    #                   -(g+ + g-)/2 ± i*omega (T_2 with Hamiltonian rotation).
    evs = eigvals(toy.L_super)
    perm = sortperm(evs; by = z -> abs(real(z)))
    evs = evs[perm]
    println("  L_super eigenvalues     :")
    for (k, lam) in enumerate(evs)
        @printf "      lambda_%d = %+.10f %+.10fi\n" k real(lam) imag(lam)
    end
    @assert abs(real(evs[1])) < 1e-12 && abs(imag(evs[1])) < 1e-12

    others = evs[2:end]
    real_ones = filter(z -> abs(imag(z)) < 1e-10, others)
    imag_ones = filter(z -> abs(imag(z)) >= 1e-10, others)
    @assert length(real_ones) == 1
    @assert length(imag_ones) == 2
    @assert isapprox(real(real_ones[1]), -toy.decay_rate; rtol=1e-12)
    for z in imag_ones
        @assert isapprox(real(z), -toy.decay_rate / 2; rtol=1e-12)
    end
    println("  L_super spectrum matches Davies T_1 (-(g+ + g-)) and T_2 (-(g+ + g-)/2 +- iw) modes.")
    return nothing
end


# -----------------------------------------------------------------------------
# L-mode integrator: Arnoldi on L
# -----------------------------------------------------------------------------

"""
    lindblad_action_integrate(L_apply!, rho_0, sigma_beta, t_grid;
                              krylovdim=30, tol=1e-10, save_states=false)

Integrate `d/dt rho = L(rho)` along `t_grid` using `KrylovKit.exponentiate`
on the matrix-free closure `L_apply!(out, in)`.

# Arguments
- `L_apply!::F`                     :  in-place closure writing `L(in)` into `out`.
- `rho_0::Matrix{T}`                :  initial dxd density operator (T<:Complex).
- `sigma_beta::Matrix{T}`           :  reference state for the trace-distance metric.
- `t_grid::AbstractVector{<:Real}`  :  ordered time grid (t_grid[1] is the start).

# Keywords
- `krylovdim::Int=30`               :  Arnoldi subspace size.
- `tol::Real=1e-10`                 :  KrylovKit per-step tolerance.
- `save_states::Bool=false`         :  also return the full trajectory of `rho` at each grid point.

# Returns
NamedTuple with `t`, `distances`, `rho_final`, `total_matvecs`, `all_converged`,
and `states::Vector{Matrix{T}}` if `save_states=true`.
"""
function lindblad_action_integrate(
    L_apply!::F,
    rho_0::Matrix{T},
    sigma_beta::Matrix{T},
    t_grid::AbstractVector{<:Real};
    krylovdim::Int = 30,
    tol::Real = 1e-10,
    save_states::Bool = false,
)::NamedTuple where {T<:Complex, F}
    d = size(rho_0, 1)
    @assert size(rho_0, 2) == d  "rho_0 must be square"

    # Persistent buffers captured by the closure (no allocations in the hot loop).
    rho_buf = Matrix{T}(undef, d, d)
    out_buf = Matrix{T}(undef, d, d)

    # KrylovKit-side closure on flat vectors.
    # Pitfall 1 (per krylov_eigsolve.jl:373-378): the buffer `out_buf` is overwritten
    # on the next call, so we MUST `copy(vec(out_buf))` before returning. We also
    # `copyto!` the input into our private buffer in case `v` aliases internal state.
    function L_vec_apply(v::AbstractVector)
        copyto!(rho_buf, reshape(v, d, d))
        L_apply!(out_buf, rho_buf)
        return copy(vec(out_buf))
    end

    n_steps   = length(t_grid)
    distances = Vector{Float64}(undef, n_steps)
    states    = save_states ? Vector{Matrix{T}}(undef, n_steps) : Matrix{T}[]

    # Initial state: defensive copy + trace-distance.
    rho     = copy(rho_0)
    distances[1] = sum(svdvals(rho - sigma_beta)) / 2
    save_states && (states[1] = copy(rho))

    # Working flat vector for Krylov.
    v_rho   = copy(vec(rho))

    total_matvecs = 0
    all_converged = true

    @inbounds for i in 1:(n_steps - 1)
        dt = float(t_grid[i + 1] - t_grid[i])

        v_next, info = exponentiate(L_vec_apply, dt, v_rho;
                                    krylovdim = krylovdim,
                                    tol = tol,
                                    ishermitian = false)
        total_matvecs += info.numops
        if info.converged == 0
            all_converged = false
            @warn "L-mode exponentiate did not converge at step" i numops=info.numops
        end
        copyto!(v_rho, v_next)

        # Reshape into rho, then defensively re-Hermitise + re-trace-normalise:
        # Davies preserves Hermiticity and trace exactly in the continuum, but
        # Krylov truncation introduces O(tol)-level violations.
        copyto!(rho, reshape(v_rho, d, d))
        @inbounds for j in 1:d, k in 1:d
            rho[k, j] = (rho[k, j] + conj(rho[j, k])) / 2
        end
        # rho is now Hermitian; re-vec and renormalise trace (only the real diagonal).
        tr_now = real(tr(rho))
        if tr_now != 0
            rho ./= tr_now
        end
        copyto!(v_rho, vec(rho))

        distances[i + 1] = sum(svdvals(rho - sigma_beta)) / 2
        save_states && (states[i + 1] = copy(rho))
    end

    return (
        t              = collect(t_grid),
        distances      = distances,
        rho_final      = copy(rho),
        total_matvecs  = total_matvecs,
        all_converged  = all_converged,
        states         = states,
    )
end


# -----------------------------------------------------------------------------
# K-mode integrator: Lanczos on K = D(sigma, L_D) (HS-self-adjoint under KMS-DB)
# -----------------------------------------------------------------------------

"""
    discriminant_action_integrate(K_apply!, psi_0, psi_eq, t_grid;
                                  krylovdim=30, tol=1e-10,
                                  is_hermitian=true, save_states=false)

Integrate `d/dt psi = K(psi)` along `t_grid` using `KrylovKit.exponentiate`
on the matrix-free closure `K_apply!(out, in)`.

When `K = D(sigma, L_D)` and `L_D` satisfies KMS-DB w.r.t. `sigma`, `K` is
HS-self-adjoint, so Lanczos (`is_hermitian=true`) is correct and faster.
For non-DB Lindbladians (or when including the bare-H commutator without
a Lamb-shift correction), set `is_hermitian=false` to force Arnoldi.

# Returns
NamedTuple with `t`, `distances`, `psi_final`, `total_matvecs`, `all_converged`,
and `states::Vector{Matrix{T}}` if `save_states=true`.
"""
function discriminant_action_integrate(
    K_apply!::F,
    psi_0::Matrix{T},
    psi_eq::Matrix{T},
    t_grid::AbstractVector{<:Real};
    krylovdim::Int = 30,
    tol::Real = 1e-10,
    is_hermitian::Bool = true,
    save_states::Bool = false,
)::NamedTuple where {T<:Complex, F}
    d = size(psi_0, 1)
    @assert size(psi_0, 2) == d  "psi_0 must be square"

    in_buf  = Matrix{T}(undef, d, d)
    out_buf = Matrix{T}(undef, d, d)

    function K_vec_apply(v::AbstractVector)
        copyto!(in_buf, reshape(v, d, d))
        K_apply!(out_buf, in_buf)
        return copy(vec(out_buf))
    end

    n_steps   = length(t_grid)
    distances = Vector{Float64}(undef, n_steps)
    states    = save_states ? Vector{Matrix{T}}(undef, n_steps) : Matrix{T}[]

    psi = copy(psi_0)
    distances[1] = norm(psi - psi_eq)             # Frobenius distance (chi metric)
    save_states && (states[1] = copy(psi))

    v_psi = copy(vec(psi))
    total_matvecs = 0
    all_converged = true

    @inbounds for i in 1:(n_steps - 1)
        dt = float(t_grid[i + 1] - t_grid[i])

        v_next, info = exponentiate(K_vec_apply, dt, v_psi;
                                    krylovdim = krylovdim,
                                    tol = tol,
                                    ishermitian = is_hermitian)
        total_matvecs += info.numops
        if info.converged == 0
            all_converged = false
            @warn "K-mode exponentiate did not converge at step" i numops=info.numops
        end
        copyto!(v_psi, v_next)
        copyto!(psi, reshape(v_psi, d, d))

        # No re-Hermitisation/normalisation here: psi is in the chi-metric
        # representation and need not be Hermitian or trace-1; rounding errors
        # stay bounded in F-norm by the per-step tol.

        distances[i + 1] = norm(psi - psi_eq)
        save_states && (states[i + 1] = copy(psi))
    end

    return (
        t              = collect(t_grid),
        distances      = distances,
        psi_final      = copy(psi),
        total_matvecs  = total_matvecs,
        all_converged  = all_converged,
        states         = states,
    )
end


# -----------------------------------------------------------------------------
# Asymptotic-slope tau_mix estimator
# -----------------------------------------------------------------------------

"""
    estimate_tau_mix_slope(t, dist, t_lo, t_hi)
        -> (tau_mix, slope, intercept, ix_lo, ix_hi)

Closed-form linear regression of `log(dist) vs t` on the window
`t in [t_lo, t_hi]`. Returns `tau_mix = -1 / slope` and the fit details.
"""
function estimate_tau_mix_slope(t::AbstractVector{<:Real}, dist::AbstractVector{<:Real},
                                t_lo::Real, t_hi::Real)
    ix = findall(τ -> t_lo <= τ <= t_hi, t)
    @assert length(ix) >= 5  "fit window too small (need at least 5 points)"
    x = Float64.(@view t[ix])
    y = log.(Float64.(@view dist[ix]))
    @assert all(isfinite, y)  "log(dist) is non-finite in fit window — distance hit floor"
    x_bar = sum(x) / length(x)
    y_bar = sum(y) / length(y)
    cov_xy = sum((x .- x_bar) .* (y .- y_bar))
    var_x  = sum((x .- x_bar) .^ 2)
    slope     = cov_xy / var_x
    intercept = y_bar - slope * x_bar
    return (tau_mix = -1.0 / slope, slope = slope, intercept = intercept,
            ix_lo = first(ix), ix_hi = last(ix))
end


# -----------------------------------------------------------------------------
# Main validation block
# -----------------------------------------------------------------------------

function main()
    # ---- Parameters --------------------------------------------------------
    omega = 1.0    # PHYSICS CHECK: Hamiltonian gap; sets time scale (omega = 1 is the natural unit).
    beta  = 2.0    # PHYSICS CHECK: moderate temperature; KMS rate ratio e^{-2} ≈ 0.135 is well-resolved
                    #                 (avoids beta -> 0 trivial limit and beta -> infinity sigma^{-1/4} blow-up).
    gamma = 0.5    # PHYSICS CHECK: jump rate; sets decay scale via gamma_+ + gamma_- = 0.5(1 + e^{-2}) ≈ 0.5677.

    toy = build_davies_2x2(omega, beta, gamma)
    println("=== Single-qubit Davies thermal Lindbladian validation ===")
    @printf "  omega = %.4f, beta = %.4f, gamma = %.4f\n" omega beta gamma
    @printf "  gamma_+ = %.10f, gamma_- = %.10f\n" toy.gamma_plus toy.gamma_minus
    @printf "  decay rate (gamma_+ + gamma_-) = %.10f\n" toy.decay_rate
    tau_mix_analytic = 1.0 / toy.decay_rate
    @printf "  analytic tau_mix = 1/(gamma_+ + gamma_-) = %.10f\n" tau_mix_analytic

    _print_construction_checks(toy)

    # ---- Time grid and initial state --------------------------------------
    t_max  = 10.0 / toy.decay_rate
    # PHYSICS CHECK: t_max = 10*tau_mix gives ~10 e-foldings — sufficient for an
    # asymptotic-slope fit while keeping signal above the Krylov tol floor.
    n_grid = 121
    t_grid = collect(range(0.0, t_max, length = n_grid))
    @printf "\n--- Integration setup ---\n  t_max = %.6f, n_grid = %d, dt = %.6f\n" t_max n_grid (t_grid[2]-t_grid[1])

    # PHYSICS CHECK: |1><1| is diagonal (and is the *excited* state when |0> is ground).
    # Lindbladian dynamics live entirely in the population (T_1) sector — the
    # leading decay rate is exactly gamma_+ + gamma_-. Off-diagonal initial state
    # would bring in the slower T_2 = (gamma_+ + gamma_-)/2 mode plus iw rotation.
    rho_0 = Matrix{ComplexF64}(undef, 2, 2)
    rho_0 .= 0
    rho_0[2, 2] = 1

    # psi_0 = sigma^{-1/4} rho_0 sigma^{-1/4};  psi_eq = sigma^{1/2}.
    psi_0  = Matrix{ComplexF64}(undef, 2, 2)
    @inbounds for j in 1:2, i in 1:2
        psi_0[i, j] = toy.sigma_inv_quarter[i] * rho_0[i, j] * toy.sigma_inv_quarter[j]
    end
    psi_eq = Matrix{ComplexF64}(undef, 2, 2)
    psi_eq .= 0
    psi_eq[1, 1] = toy.sigma_half[1]
    psi_eq[2, 2] = toy.sigma_half[2]

    # Sanity: chi^2 == || psi_0 - psi_eq ||_F^2.
    chi2_direct = 0.0
    let
        diff = rho_0 - toy.sigma_beta
        for j in 1:2, i in 1:2
            chi2_direct += real(conj(diff[i, j]) * diff[i, j]) *
                           toy.sigma_inv_quarter[i]^2 * toy.sigma_inv_quarter[j]^2
        end
    end
    chi2_via_psi = norm(psi_0 - psi_eq)^2
    @printf "  chi^2(rho_0, sigma) (direct)   = %.10f\n" chi2_direct
    @printf "  chi^2(rho_0, sigma) (via psi)  = %.10f\n" chi2_via_psi
    @assert isapprox(chi2_direct, chi2_via_psi; rtol=1e-12)

    # ---- Run integrators ---------------------------------------------------
    krylovdim = 8     # PHYSICS CHECK: 4-dim space; krylovdim=8 is overkill but cheap and
                       #                guarantees machine-precision per-step expv.
    tol_step  = 1e-12

    println("\n--- Running integrators ---")
    t0 = time()
    res_L = lindblad_action_integrate(toy.L_apply!, rho_0, toy.sigma_beta, t_grid;
                                      krylovdim = krylovdim, tol = tol_step,
                                      save_states = true)
    res_K = discriminant_action_integrate(toy.K_apply!, psi_0, psi_eq, t_grid;
                                          krylovdim = krylovdim, tol = tol_step,
                                          is_hermitian = true, save_states = true)
    wall_time = time() - t0

    @printf "  L-mode total matvecs : %d   (all converged = %s)\n" res_L.total_matvecs string(res_L.all_converged)
    @printf "  K-mode total matvecs : %d   (all converged = %s)\n" res_K.total_matvecs string(res_K.all_converged)
    @printf "  wall time            : %.4f s\n" wall_time

    # ---- Slope fit on [t_max/3, 2 t_max/3] window -------------------------
    # PHYSICS CHECK: window past initial transient and well above Krylov tol floor;
    # for diagonal init on Davies the decay is essentially mono-exponential, so any
    # internal sub-window of the decay tail gives the same slope.
    fit_lo = t_max / 3
    fit_hi = 2 * t_max / 3
    fit_L = estimate_tau_mix_slope(res_L.t, res_L.distances, fit_lo, fit_hi)
    fit_K = estimate_tau_mix_slope(res_K.t, res_K.distances, fit_lo, fit_hi)

    @printf "\n--- Slope fits on t in [%.4f, %.4f] ---\n" fit_lo fit_hi
    @printf "  L-mode tau_mix = %.10f  (slope = %+.10f)\n" fit_L.tau_mix fit_L.slope
    @printf "  K-mode tau_mix = %.10f  (slope = %+.10f)\n" fit_K.tau_mix fit_K.slope

    # ---- Trajectory equivalence: rho^L vs rho-from-psi^K ------------------
    # rho_from_psi[i,j] = sigma^{1/4}_i * psi[i,j] * sigma^{1/4}_j.
    # NOTE: equivalent only on the diagonal (Bohr-zero) sector here, since
    # L-mode uses the full L (with bare-H rotation in the off-diagonal block)
    # while K-mode uses the dissipator-only L_D. For diagonal rho_0 the
    # off-diagonal block is identically zero in both modes -> equivalence holds.
    max_traj_mismatch = 0.0
    for k in 1:length(t_grid)
        psi_k = res_K.states[k]
        rho_from_psi = Matrix{ComplexF64}(undef, 2, 2)
        @inbounds for j in 1:2, i in 1:2
            rho_from_psi[i, j] = toy.sigma_quarter[i] * psi_k[i, j] * toy.sigma_quarter[j]
        end
        diff = norm(res_L.states[k] - rho_from_psi)
        max_traj_mismatch = max(max_traj_mismatch, diff)
    end
    @printf "\n--- Trajectory equivalence (rho^L vs sigma^{1/4} psi^K sigma^{1/4}) ---\n"
    @printf "  max ||rho^L_i - rho_from_psi^K_i||_F = %.3e\n" max_traj_mismatch

    # ---- Density-matrix invariants along L-mode trajectory ----------------
    min_trace = Inf;  max_trace = -Inf
    min_eig   = Inf;  max_eig   = -Inf
    for rho_i in res_L.states
        tr_i = real(tr(rho_i))
        min_trace = min(min_trace, tr_i)
        max_trace = max(max_trace, tr_i)
        # rho_i is already (rho+rho')/2 from the integrator.
        evs_i = eigvals(Hermitian((rho_i + rho_i') / 2))
        min_eig = min(min_eig, minimum(evs_i))
        max_eig = max(max_eig, maximum(evs_i))
    end
    @printf "\n--- L-mode density-matrix invariants ---\n"
    @printf "  trace : min = %.16f, max = %.16f\n" min_trace max_trace
    @printf "  eigs  : min = %.3e, max = %.6f\n" min_eig max_eig

    # ---- Mid-time direct dense-expv cross-check ---------------------------
    t_mid = t_max / 2
    ix_mid = findfirst(τ -> τ >= t_mid, t_grid)
    @assert ix_mid !== nothing
    t_mid_actual = t_grid[ix_mid]
    v_ref = exp(t_mid_actual * toy.L_super) * vec(rho_0)
    rho_ref = reshape(v_ref, 2, 2)
    rho_at_t = res_L.states[ix_mid]
    mid_mismatch = norm(rho_at_t - rho_ref)
    @printf "\n--- Mid-time dense-expv cross-check at t = %.6f (grid index %d) ---\n" t_mid_actual ix_mid
    @printf "  ||rho^Krylov - rho^dense||_F = %.3e\n" mid_mismatch

    # ---- Acceptance checks -----------------------------------------------
    rel_err_L     = abs(fit_L.tau_mix - tau_mix_analytic) / tau_mix_analytic
    rel_err_K     = abs(fit_K.tau_mix - tau_mix_analytic) / tau_mix_analytic
    pass_L        = rel_err_L < 0.01
    pass_K        = rel_err_K < 0.01
    pass_traj     = max_traj_mismatch < 1e-10
    pass_invariants = (abs(min_trace - 1) < 1e-10) && (abs(max_trace - 1) < 1e-10) &&
                      (min_eig > -1e-10)
    pass_mid      = mid_mismatch < 1e-9
    all_pass      = pass_L && pass_K && pass_traj && pass_invariants && pass_mid &&
                    res_L.all_converged && res_K.all_converged

    println("\n=== SUMMARY ===")
    @printf "Analytic tau_mix         : %.6f  (= 1/(gamma_+ + gamma_-) at omega=%.2f, beta=%.2f, gamma=%.2f)\n" tau_mix_analytic omega beta gamma
    @printf "L-mode   tau_mix         : %.6f  (rel err = %.3e %s)\n" fit_L.tau_mix rel_err_L (pass_L ? "PASS" : "FAIL")
    @printf "K-mode   tau_mix         : %.6f  (rel err = %.3e %s)\n" fit_K.tau_mix rel_err_K (pass_K ? "PASS" : "FAIL")
    @printf "Max trajectory mismatch  : %.3e  (%s)\n" max_traj_mismatch (pass_traj ? "PASS" : "FAIL")
    @printf "Mid-time expv mismatch   : %.3e  (%s)\n" mid_mismatch (pass_mid ? "PASS" : "FAIL")
    @printf "Density-matrix invariants: trace in [%.3e, %.3e], min eig = %.3e  (%s)\n" (min_trace - 1) (max_trace - 1) min_eig (pass_invariants ? "PASS" : "FAIL")
    @printf "Total matvecs (L-mode)   : %d\n" res_L.total_matvecs
    @printf "Total matvecs (K-mode)   : %d\n" res_K.total_matvecs
    @printf "Wall time                : %.4f s\n" wall_time
    @printf "All checks               : %s\n" (all_pass ? "PASS" : "FAIL")

    return all_pass ? 0 : 1
end


exit(main())
