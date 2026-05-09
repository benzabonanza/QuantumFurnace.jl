# ============================================================================
# Matrix-free Krylov-subspace integrators for Lindbladian and discriminant flows
# ============================================================================
#
# Two equivalent flows on a finite-dimensional state space:
#
#   L-mode :  d/dt rho(t)  = L(rho(t))                    (Lindbladian on density operator)
#   K-mode :  d/dt psi(t)  = K(psi(t)),  K = D(sigma, L_D) (KMS quantum discriminant action)
#
# In Heisenberg-picture (column-stacking vec convention),
# rho_vec(t) = exp(t L_super) rho_vec(0)  and
# psi(t) = sigma^{-1/4} rho(t) sigma^{-1/4}, so the two flows are exact
# reparametrisations of each other on the diagonal (Bohr-zero) sector of
# operator space, where the H commutator vanishes.
#
# Equation references:
# - eq:discriminant : D(X) = sigma^{-1/4} L(sigma^{1/4} X sigma^{1/4}) sigma^{-1/4}
#                     (thesis Eq. eq:discriminant; Chen-Kastoryano-Gilyen 2023; Chen 2025 Sec. III)
# - chi^2 distance  : chi^2(rho, sigma) = || sigma^{-1/4}(rho - sigma) sigma^{-1/4} ||_F^2
# - KrylovKit       : `exponentiate(A, t, v; ishermitian, krylovdim, tol)` returns (y, info)
#                     with `info.numops` matvec count and `info.converged in {0, 1}`.
#                     `ishermitian=true` short-circuits to Lanczos; `false` to Arnoldi.

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

    # KrylovKit-side closure on flat vectors. KrylovKit may overwrite the
    # buffer `out_buf` on the next call, so we MUST `copy(vec(out_buf))`
    # before returning; we also `copyto!` the input into our private buffer
    # in case `v` aliases internal state.
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

    # KrylovKit may overwrite `out_buf` on the next call, so we
    # `copy(vec(out_buf))` before returning; we also `copyto!` the input
    # into our private buffer in case `v` aliases internal state.
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

        # Defensive re-Hermitisation + invariant renormalisation, mirroring the
        # L-mode integrator. KMS-DB ⇒ K preserves ψ = ψ^† and ⟨ψ_eq, ψ⟩_F = tr(ρ)
        # exactly in the continuum (= 1 for a normalised initial state, since
        # ⟨σ^{1/2}, σ^{-1/4} ρ σ^{-1/4}⟩_F = tr(ρ)). Krylov truncation introduces
        # O(tol) violations that accumulate over many small steps and matter
        # for tight target_epsilon (qf-lkb.10).
        @inbounds for j in 1:d, k in 1:d
            psi[k, j] = (psi[k, j] + conj(psi[j, k])) / 2
        end
        c_now = real(dot(psi_eq, psi))
        if c_now != 0
            psi ./= c_now
        end
        copyto!(v_psi, vec(psi))

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


"""
    integrate_to_gibbs(config, hamiltonian, jumps, rho_0, t_grid;
                       mode=:L, krylovdim=30, tol=1e-10, save_states=false)

End-to-end wrapper around `lindblad_action_integrate` (mode=:L) and
`discriminant_action_integrate` (mode=:K). Supports KMS in `BohrDomain`
or `EnergyDomain` (matrix-free `apply_lindbladian!`), and DLL in
`BohrDomain` (qf-lkb.9 specialised workspace). Both paths reach the same
Krylov-subspace integrator and return the same NamedTuple shape.

EnergyDomain CKG (qf-lkb.11) is the production path for n ≥ 5 sweeps:
its matvec scales as `O(n_jumps · N · d²)` with `N = 2^num_energy_bits`
fixed by the OFT grid, vs BohrDomain's `O(d⁴)` over Bohr-pair indices.

# Arguments
- `config::Config{Lindbladian, <:Union{BohrDomain, EnergyDomain}}`:
  simulator config. DLL is restricted to BohrDomain by `validate_config!`.
- `hamiltonian::HamHam`                     :  Hamiltonian + cached Gibbs state.
- `jumps::Vector{JumpOp}`                   :  jump operators (computational + eigenbasis).
- `rho_0::Matrix{T}`                        :  initial dxd density operator (T<:Complex).
- `t_grid::AbstractVector{<:Real}`          :  ordered time grid.

# Keywords
- `mode::Symbol = :L`                       :  `:L` for Lindbladian flow on rho;
                                                `:K` for discriminant flow on
                                                psi = sigma^{-1/4} rho sigma^{-1/4}.
- `krylovdim::Int = 30`                     :  forwarded to KrylovKit.exponentiate.
- `tol::Real = 1e-10`                       :  forwarded to KrylovKit.exponentiate.
- `save_states::Bool = false`               :  forwarded.

# Returns
NamedTuple from the underlying integrator: `t`, `distances`,
`rho_final` (mode=:L) or `psi_final` (mode=:K), `total_matvecs`,
`all_converged`, and `states::Vector{Matrix{T}}` if `save_states=true`.

# Caveats
- BohrDomain or EnergyDomain only by signature dispatch; `validate_config!`
  is called defensively at the start so DLL EnergyDomain/TrotterDomain
  and TimeDomain/TrotterDomain are rejected upstream.
- For DLL, `config.filter` must be a `DLLGaussianFilter(beta)` or
  `DLLMetropolisFilter(beta)` matching `config.beta`; this is also enforced
  by `validate_config!`.
- mode=:K hardcodes `is_hermitian=true` because both BohrDomain and
  EnergyDomain + KMS / DLL satisfy KMS detailed balance by construction
  (the discriminant is HS-self-adjoint).
- The integrator internals defensively rebalance Hermiticity and trace at the
  Krylov-tol level (mode=:L); the K-mode does not (psi is in the chi-metric
  representation and need not be Hermitian or trace-1).
"""
function integrate_to_gibbs(
    config::Config{Lindbladian, <:Union{BohrDomain, EnergyDomain}},
    hamiltonian::HamHam,
    jumps::Vector{JumpOp},
    rho_0::Matrix{T},
    t_grid::AbstractVector{<:Real};
    mode::Symbol = :L,
    krylovdim::Int = 30,
    tol::Real = 1e-10,
    save_states::Bool = false,
    allow_unpaired_nonhermitian::Bool = false,
)::NamedTuple where {T<:Complex}
    mode in (:L, :K) || throw(ArgumentError("mode must be :L or :K (got :$mode)"))
    d = size(rho_0, 1)
    @assert size(rho_0, 2) == d  "rho_0 must be square"

    # validate_config! is invoked by run_lindblad/run_thermalize; we call it
    # explicitly here since this entry point bypasses those.
    validate_config!(config)
    validate_jump_pairing(jumps; allow_unpaired_nonhermitian=allow_unpaired_nonhermitian)

    # Build the L_apply!(out, in) closure: matrix-free for both KMS and DLL.
    # `let` scope binds the captured state directly to dodge Box wrapping under
    # Julia 1.11+ closure capture rules.
    ws = Workspace(config, hamiltonian, jumps)
    L_apply! = let ws = ws, config = config, ham = hamiltonian
        (out::AbstractMatrix, x::AbstractMatrix) -> begin
            apply_lindbladian!(ws, x, config, ham)
            copyto!(out, ws.scratch.rho_out)
            return out
        end
    end

    if mode == :L
        sigma_beta = Matrix{T}(hamiltonian.gibbs)
        return lindblad_action_integrate(
            L_apply!, rho_0, sigma_beta, t_grid;
            krylovdim = krylovdim, tol = tol, save_states = save_states,
        )
    else  # mode == :K
        powers = gibbs_fractional_powers(hamiltonian.gibbs)
        sq, sq_inv, sh = powers.sigma_quarter, powers.sigma_inv_quarter, powers.sigma_half

        # psi_0 = sigma^{-1/4} rho_0 sigma^{-1/4} (diagonal multiply, BohrDomain).
        psi_0 = Matrix{T}(undef, d, d)
        @inbounds for j in 1:d, i in 1:d
            psi_0[i, j] = sq_inv[i] * rho_0[i, j] * sq_inv[j]
        end
        # psi_eq = sigma^{1/2} as a full Matrix (the integrator wants a Matrix, not Diagonal).
        psi_eq = Matrix{T}(Diagonal(complex.(sh)))

        bufs = DiscriminantBuffers{T}(d)
        K_apply! = let L = L_apply!, sq = sq, sq_inv = sq_inv, bufs = bufs
            (out::AbstractMatrix, x::AbstractMatrix) -> begin
                apply_discriminant!(out, x, L, sq, sq_inv, bufs)
                return out
            end
        end

        return discriminant_action_integrate(
            K_apply!, psi_0, psi_eq, t_grid;
            krylovdim = krylovdim, tol = tol,
            is_hermitian = true,  # KMS-DB ⇒ K is HS-self-adjoint (Lanczos OK)
            save_states = save_states,
        )
    end
end


# ---------------------------------------------------------------------------
# qf-ev5.{3,4}: Krylov spectral-expansion trajectory
# ---------------------------------------------------------------------------
#
# Idea: instead of integrating d/dt rho = L(rho) step-by-step, build the
# leading h eigenpairs of L matrix-free (right + adjoint Krylov), biorthogonalise,
# project rho_0 - rho_inf onto the slow-mode subspace, and evaluate
#
#   rho(t) = rho_inf + sum_{i ne 0}  c_i * exp(lambda_i * t) * R_i
#
# at any t in O(h * d^2). Cost is dominated by the two Krylov eigsolves
# (~ h * krylov_iter matvecs, typically ~30-100 for tight convergence), to be
# compared with `lindblad_action_integrate`'s ~ length(t_grid) *
# krylov_exp_subspace ~ 1000-3000 matvecs. Both routes return the same
# NamedTuple shape so the bi-exp + mixing pipeline plugs in unchanged.
#
# Adjoint subtlety: L preserves trace ⇒ L*(I) = 0 exactly, so vec(I/d) is an
# exact eigenvector of L* and Arnoldi terminates on the first matvec. The
# adjoint eigsolve must therefore start from a *generic* Hermitian seed.

"""
    _arnoldi_factorize(f, x0, m) -> (Q, H, broke)

Plain Arnoldi factorization: build orthonormal basis `Q` and upper-Hessenberg
`H` such that `f(Q[:, j]) ≈ Q[:, 1:j+1] * H[1:j+1, j]` for each j ≤ m.
Returns the full `m`-step factorization, or a smaller one if the basis breaks
down before step `m` (in which case `broke=true`). Uses modified
Gram–Schmidt with one reorthogonalisation pass for stability.

Output: `Q::Matrix{T}` (size N × m_eff), `H::Matrix{T}` (size m_eff × m_eff).
"""
function _arnoldi_factorize(f, x0::AbstractVector{T}, m::Int) where {T}
    N = length(x0)
    Q = zeros(T, N, m + 1)
    H = zeros(T, m + 1, m)
    Q[:, 1] .= x0 ./ norm(x0)
    broke_at = m
    @inbounds for j in 1:m
        w = f(Q[:, j])
        # Modified Gram–Schmidt
        for i in 1:j
            H[i, j] = dot(view(Q, :, i), w)
            w .-= H[i, j] .* view(Q, :, i)
        end
        # Reorthogonalisation pass for numerical stability
        for i in 1:j
            corr = dot(view(Q, :, i), w)
            H[i, j] += corr
            w .-= corr .* view(Q, :, i)
        end
        h_jp1 = norm(w)
        H[j + 1, j] = h_jp1
        if h_jp1 < eps(real(T)) * sqrt(N)
            broke_at = j
            break
        end
        Q[:, j + 1] .= w ./ h_jp1
    end
    return Q[:, 1:broke_at], H[1:broke_at, 1:broke_at], broke_at < m
end


"""
    _krylov_spectral_decomposition(forward_apply!, rho_0, dim; krylovdim,
                                   tol, fwd_init)

Compute a biorthogonal eigendecomposition of a generic linear operator F on
operator space (matvec via `forward_apply!(out, in)`) restricted to a
Krylov subspace of dimension `krylovdim`, and project `rho_0 - rho_inf`
onto the slow-mode subspace.

Algorithm (single forward Arnoldi, no adjoint Krylov needed):
1. Run plain Arnoldi: `f(Q[:, j]) ≈ Q[:, 1:j+1] * H[1:j+1, j]`, with `Q`
   orthonormal in the standard ℓ² (HS) inner product.
2. Diagonalise the small Hessenberg `H = W * Λ * W^{-1}` densely.
3. Right eigvecs in operator space: `R_i = reshape(Q * W[:, i], dim, dim)`.
4. Left eigvecs (biorthogonal duals): `L_i = reshape(Q * V[:, i], dim, dim)`
   where `V = (W^{-1})'`. By construction `<L_i, R_j>_HS = delta_{ij}` to
   the orthonormality of `Q` and the inverse relation `V' W = I`.
5. Coefficients: `c_i = <L_i, rho_0 - rho_inf>_HS = (V' Q' vec(rho_0 - rho_inf))[i]`.

The fixed point `rho_inf` is taken as the eigenvector with smallest
`|Re(lambda)|`, hermitised + trace-normalised.

# Returns
NamedTuple with
- `eigenvalues::Vector{ComplexF64}`  — eigenvalues sorted by |Re| ascending
- `R_modes::Vector{Matrix}`          — right eigvecs (R_modes[1] = rho_inf)
- `L_modes::Vector{Matrix}`          — biorthogonal left eigvecs
- `c::Vector{ComplexF64}`            — biorthogonal projection coefficients
- `rho_inf::Matrix`                  — hermitised, trace-normalised steady state
- `matvec_count::Int`                — total matvecs of `forward_apply!`
- `converged::Bool`                  — true if Arnoldi did not break down
"""
function _krylov_spectral_decomposition(
    forward_apply!::F1,
    rho_0::Matrix{T},
    dim::Integer;
    krylovdim::Integer = 40,
    tol::Real = 1e-10,
    fwd_init::Union{Nothing, AbstractVector} = nothing,
    sort_mode::Symbol = :lindbladian,
) where {T<:Complex, F1}

    dim2 = dim * dim

    # Closure on flat vectors with private scratch.
    rho_buf = Matrix{T}(undef, dim, dim)
    out_buf = Matrix{T}(undef, dim, dim)
    function fwd_vec(v::AbstractVector)
        copyto!(rho_buf, reshape(v, dim, dim))
        forward_apply!(out_buf, rho_buf)
        return copy(vec(out_buf))
    end

    # Arnoldi seed: rho_0. Generically has nonzero overlap with all relevant
    # slow modes so Arnoldi locks onto them. Pure I/d would also work for
    # the forward problem (only L^*(I) = 0 makes I the trivial adjoint
    # eigenvector; for the forward L applied to I/d, the slow non-steady
    # components are typically O(1) for non-trivial Hamiltonians).
    x0 = if fwd_init === nothing
        vec(rho_0)
    else
        Vector{T}(fwd_init)
    end
    m_target = min(Int(krylovdim), dim2)

    # Single forward Arnoldi factorisation. No adjoint Krylov needed: the
    # biorthogonal left eigvecs are derived from the small dense
    # eigendecomposition of H, which is unitarily related to a sub-block
    # of L's Schur form. This is tighter than running an independent
    # adjoint Krylov and matching by eigenvalue (which fails when forward
    # and adjoint Arnoldi capture different fast-mode subsets).
    Q, H, broke = _arnoldi_factorize(fwd_vec, x0, m_target)
    m = size(H, 1)
    m >= 2 || error("Arnoldi broke down with only $m basis vectors; cannot decompose")

    # Diagonalise H densely. eigen returns columns of W as right eigvecs;
    # we form V = (W^{-1})' so V' W = I (the biorthogonality condition).
    F = eigen(H)
    Λ = F.values
    W = F.vectors
    W_inv = inv(W)
    V = Matrix{T}(W_inv')  # so V' = W_inv ⇒ V' W = I

    # Sort eigenvalues so the steady state lands at index 1.
    # :lindbladian -> Re(lambda) ~ 0 is steady ⇒ sort by |Re(lambda)| ascending
    # :channel     -> |mu| ~ 1 is steady ⇒ sort by |1 - mu| ascending
    sort_key = if sort_mode === :lindbladian
        v -> abs(real(v))
    elseif sort_mode === :channel
        v -> abs(abs(v) - 1.0)
    else
        throw(ArgumentError("sort_mode must be :lindbladian or :channel (got :$sort_mode)"))
    end
    perm = sortperm(Λ; by = sort_key)
    Λ = Λ[perm]
    W = W[:, perm]
    V = V[:, perm]

    # Lift small-space eigvecs to operator space.
    QW = Q * W            # right eigvecs as flat columns
    QV = Q * V            # left eigvecs as flat columns
    R_modes = [reshape(QW[:, i], dim, dim) for i in 1:m]
    L_modes = [reshape(QV[:, i], dim, dim) for i in 1:m]

    # Steady state: hermitise R_modes[1], trace-normalise.
    rho_inf = (R_modes[1] .+ R_modes[1]') ./ 2
    tr_inf = real(tr(rho_inf))
    if tr_inf == 0
        error("steady-state mode has zero trace; cannot normalise")
    end
    rho_inf ./= tr_inf
    R_modes[1] = rho_inf

    # Biorthogonal projection: c_i = <L_i, rho_0 - rho_inf>_HS.
    # Using the lift QV: c = V' Q' vec(rho_0 - rho_inf).
    delta_rho_vec = vec(rho_0 .- rho_inf)
    c = (V') * (Q' * delta_rho_vec)

    # Diagnostic: c[1] should be ~0 for trace-preserving F (since L^*(I) = 0
    # ⇒ vec(I)' is the steady-state left eigvec, and tr(rho_0 - rho_inf) = 0
    # for both Lindbladian L and CPTP channel Phi). Force it to 0 to
    # suppress numerical noise that would otherwise contaminate every
    # time step with a constant offset (lambda_1 ~ 0 ⇒ exp ~ 1; or
    # mu_1 ~ 1 ⇒ mu_1^k ~ 1).
    c[1] = zero(T)

    return (
        eigenvalues  = Complex{Float64}.(Λ),
        R_modes      = R_modes,
        L_modes      = L_modes,
        c            = Complex{Float64}.(c),
        rho_inf      = rho_inf,
        matvec_count = m,  # Arnoldi performs exactly m matvecs
        converged    = !broke,
    )
end


"""
    predict_lindbladian_trajectory(config::Config{Lindbladian, <:Union{BohrDomain, EnergyDomain}},
                                   hamiltonian, jumps, rho_0, t_grid;
                                   krylovdim=40, tol=1e-10, save_states=false)

Reconstruct `rho(t) = e^{tL} rho_0` on `t_grid` via Krylov spectral
expansion of the Lindbladian, returning the same NamedTuple shape as
`lindblad_action_integrate` (so the bi-exp + mixing pipeline plugs in
unchanged).

The cost is `~ krylovdim` matvecs of `apply_lindbladian!` (one Arnoldi
factorisation), vs `~ length(t_grid) * krylov_subspace` matvecs for
`lindblad_action_integrate`. Typical cost reduction: 20-50x at n in 4-8
with `length(t_grid) ~ 80`.

# Arguments
Same as `integrate_to_gibbs(::Config{Lindbladian}, ..., mode=:L)`.

# Keywords
- `krylovdim::Int = 40`: Arnoldi subspace size = number of eigenmodes
  retained. The reconstruction is biorthogonally exact in the captured
  subspace; truncated fast modes contribute O(exp(-|Re(lambda_{m+1})| t))
  to the trace distance, where m+1 is the slowest excluded mode. For
  mixing-time fits on the slow tail this is benign; for early-time
  accuracy increase `krylovdim`.
- `tol::Real = 1e-10`: not used directly (the dense eigendecomposition of
  the small Hessenberg has no tolerance knob); kept for API symmetry with
  `lindblad_action_integrate` and to allow a future Krylov-shift-invert
  variant. Currently a no-op.
- `save_states::Bool = false`: also return the reconstructed rho(t_k).

# Returns
NamedTuple with `t`, `distances`, `rho_final`, `total_matvecs`,
`all_converged`, `states`, plus extra fields `eigenvalues`, `c`,
`spectral_gap`, `rho_inf`, `R_modes` (biorthogonal right-eigenvector
matrices, one `d×d` per Krylov mode), and `sigma_beta` (the trace-
distance reference used in `distances`, equal to `hamiltonian.gibbs`).
`distances[1]` is the trace distance of `rho_0` to `sigma_beta`;
subsequent entries are the spectral-expansion predictions. The
`R_modes` + `sigma_beta` fields enable closed-form τ_mix(ε)
extrapolation off the grid via `eigenmode_mixing_time` (`src/mixing.jl`).

# Notes
- Defensive Hermitisation is applied at every grid point.
- For DLL via the matrix-free path (`qf-lkb.9`), the same `apply_lindbladian!`
  machinery applies; this entry point dispatches through
  `Workspace(::Config{Lindbladian}, ...)` which selects KMS vs DLL.
- The biorthogonal projection uses the small dense eigendecomposition of
  the Arnoldi Hessenberg (no separate adjoint Krylov solve), which is
  essential: independent forward / adjoint Arnoldi runs on a non-normal
  generator generically converge to *different* fast-mode subsets, which
  breaks pair-by-pair biorthogonalisation. Lifting `eigen(H)` back via
  `Q` gives a guaranteed-biorthogonal pair on the captured subspace.
"""
function predict_lindbladian_trajectory(
    config::Config{Lindbladian, <:Union{BohrDomain, EnergyDomain}},
    hamiltonian::HamHam,
    jumps::Vector{JumpOp},
    rho_0::Matrix{T},
    t_grid::AbstractVector{<:Real};
    krylovdim::Integer = 40,
    tol::Real = 1e-10,
    save_states::Bool = false,
    allow_unpaired_nonhermitian::Bool = false,
    workspace::Union{Nothing, Workspace{KrylovSpectrum}} = nothing,
)::NamedTuple where {T<:Complex}
    d = size(rho_0, 1)
    @assert size(rho_0, 2) == d  "rho_0 must be square"

    validate_config!(config)
    validate_jump_pairing(jumps; allow_unpaired_nonhermitian=allow_unpaired_nonhermitian)

    # qf-qmi.2: optional pre-built Workspace reuse for parameter sweeps. The
    # ctor cost is O(d^3.x) for EnergyDomain Lindbladian + KMS (B_bohr loop),
    # vastly larger than a 30-matvec Krylov factorisation at n>=6. A caller
    # sweeping (β, σ, ε) at fixed (n, ham, jumps, domain, construction) can
    # build the workspace once and pass it here.
    ws = if workspace === nothing
        Workspace(config, hamiltonian, jumps)
    else
        @assert size(workspace.G_left, 1) == d  "workspace dim must match rho_0"
        @assert length(workspace.jumps) == length(jumps)  "workspace jump count mismatch"
        workspace
    end
    fwd! = let ws = ws, config = config, ham = hamiltonian
        (out::AbstractMatrix, x::AbstractMatrix) -> begin
            apply_lindbladian!(ws, x, config, ham)
            copyto!(out, ws.scratch.rho_out)
            return out
        end
    end

    decomp = _krylov_spectral_decomposition(
        fwd!, rho_0, d;
        krylovdim=krylovdim, tol=tol,
    )

    sigma_beta = Matrix{T}(hamiltonian.gibbs)
    n_t = length(t_grid)
    distances = Vector{Float64}(undef, n_t)
    states = save_states ? Vector{Matrix{T}}(undef, n_t) : Matrix{T}[]
    rho_t = Matrix{T}(undef, d, d)
    h = length(decomp.eigenvalues)

    @inbounds for k in 1:n_t
        t = float(t_grid[k])
        copyto!(rho_t, decomp.rho_inf)
        for i in 1:h
            # Skip the steady-state mode (its c is ~0 by trace preservation).
            abs(decomp.eigenvalues[i]) < 1e-10 && continue
            phase = exp(decomp.eigenvalues[i] * t)
            rho_t .+= (decomp.c[i] * phase) .* decomp.R_modes[i]
        end
        # Defensive Hermitisation (mirrors lindblad_action_integrate).
        @inbounds for j in 1:d, kc in 1:d
            rho_t[kc, j] = (rho_t[kc, j] + conj(rho_t[j, kc])) / 2
        end
        distances[k] = sum(svdvals(rho_t .- sigma_beta)) / 2
        save_states && (states[k] = copy(rho_t))
    end

    spectral_gap = h >= 2 ? abs(real(decomp.eigenvalues[2])) : NaN

    return (
        t              = collect(t_grid),
        distances      = distances,
        rho_final      = copy(rho_t),
        total_matvecs  = decomp.matvec_count,
        all_converged  = decomp.converged,
        states         = states,
        eigenvalues    = decomp.eigenvalues,
        c              = decomp.c,
        spectral_gap   = spectral_gap,
        rho_inf        = decomp.rho_inf,
        R_modes        = decomp.R_modes,
        sigma_beta     = sigma_beta,
    )
end


# ---------------------------------------------------------------------------
# qf-ev5.{1,5}: Faithful Φ_δ channel matvec + spectral trajectory predictor
# ---------------------------------------------------------------------------
#
# Sister of `predict_lindbladian_trajectory` for the *implemented* CPTP channel
# Φ_δ that `run_thermalize` actually executes. The forward matvec composes
# the per-jump weak-measurement substeps in :sweep order, exactly matching
# what one full δ-step of `run_thermalize` does. Single Arnoldi factorisation
# of Φ_δ then powers the same biorthogonal closed-form expansion used for
# the Lindbladian path: rho_k = rho_inf + Σ c_i  μ_i^k  R_i, with μ_i the
# channel eigenvalues (clustered near 1 for slow modes).

"""
    predict_channel_trajectory(config::Config{Thermalize, <:Union{BohrDomain,
                                EnergyDomain, TimeDomain}}, hamiltonian, jumps,
                                rho_0, k_grid; krylovdim=40, tol=1e-10,
                                save_states=false)

Reconstruct `rho_k = (Φ_δ)^k(rho_0)` on integer step grid `k_grid` via Krylov
spectral expansion of the implemented CPTP channel `Φ_δ`. The forward matvec
matches `run_thermalize` byte-for-byte in `:sweep` mode: per-jump coherent
unitary, per-jump weak-measurement Kraus, applied in jump order. The
reconstructed trajectory captures all algorithmic errors of the implemented
channel (δ-step weak measurement, coherent splitting, OFT/Trotter), unlike
the ideal-L `predict_lindbladian_trajectory` which targets `e^{tL} rho_0`.

# Arguments
- `config::Config{Thermalize, ...}`: must have `jump_selection = :sweep`
  (the deterministic Φ_δ requires a fixed jump order; `:random` would
  describe an ensemble average over jump pickings, not a single channel).
- `hamiltonian::HamHam`: with `gibbs` cached. The Gibbs state is used as
  the trace-distance reference; the channel's *actual* fixed point is
  separately recovered from the leading eigenvector and is generally
  O(δ²) close to but not equal to `hamiltonian.gibbs`.
- `jumps::Vector{JumpOp}`: jump operators (computational + eigenbasis).
- `rho_0::Matrix{T}`: initial dxd density operator.
- `k_grid::AbstractVector{<:Integer}`: integer step counts at which to
  evaluate `(Φ_δ)^k(rho_0)`. The associated physical time is `k * δ`.

# Keywords
- `krylovdim::Integer = 40`: Arnoldi subspace size.
- `tol::Real = 1e-10`: reserved (no-op currently — the dense small-space
  eigendecomposition has no tolerance knob).
- `save_states::Bool = false`: also return the reconstructed rho_k.
- `trotter::Union{Nothing, TrottTrott} = nothing`: required for TrotterDomain.

# Returns
NamedTuple with `t = collect(k_grid) .* δ`, `distances` (trace distance to
the basis-aligned Gibbs reference at each step), `rho_final`,
`total_matvecs`, `all_converged`, `states`, plus `eigenvalues` (channel
μ_i, sorted by |μ| descending so steady-state ~1 is first), `c`,
`spectral_gap` (in Lindbladian units, `(1 - |μ_2|) / δ`), `rho_inf`,
`R_modes` (biorthogonal right-eigenvector matrices, one `d×d` per
captured Krylov mode), `sigma_beta` (the basis-aligned Gibbs reference;
on TrotterDomain it has been transformed into the Trotter eigenbasis),
`delta_used`, `k_grid`. The `R_modes` + `sigma_beta` fields enable
closed-form τ_mix(ε) extrapolation via `eigenmode_mixing_time` after
the channel-to-Lindbladian eigenvalue conversion `λ_eff = log(μ) / δ`.

# Notes
- Limited to `jump_selection = :sweep`. `:random` runs a stochastic
  process whose deterministic mean is `e^{δ𝓛}` only in expectation; the
  step-by-step matvec is not a CPTP channel in the deterministic sense.
- The leading channel eigenvalue μ_1 should be exactly 1 (CPTP fixed
  point). Numerical deviation is reported but not enforced — small
  drifts from 1 indicate tolerance-level non-CPTPness from the weak-
  measurement assembly (`_build_cptp_channel` PSD clamp).
- The biorthogonal expansion via single Arnoldi + dense `eigen(H)` works
  unchanged for non-Hermitian Φ_δ: complex eigenvalues come in conjugate
  pairs whose contributions to the real reconstruction cancel imaginary
  parts up to machine precision.
"""
function predict_channel_trajectory(
    config::Config{Thermalize, <:Union{BohrDomain, EnergyDomain, TimeDomain, TrotterDomain}},
    hamiltonian::HamHam,
    jumps::Vector{JumpOp},
    rho_0::Matrix{T},
    k_grid::AbstractVector{<:Integer};
    krylovdim::Integer = 40,
    tol::Real = 1e-10,
    save_states::Bool = false,
    trotter::Union{Nothing, TrottTrott} = nothing,
    allow_unpaired_nonhermitian::Bool = false,
    workspace::Union{Nothing, Workspace{KrylovSpectrum}} = nothing,
)::NamedTuple where {T<:Complex}
    d = size(rho_0, 1)
    @assert size(rho_0, 2) == d  "rho_0 must be square"

    validate_config!(config)
    validate_jump_pairing(jumps; allow_unpaired_nonhermitian=allow_unpaired_nonhermitian)
    config.jump_selection === :sweep || throw(ArgumentError(
        "predict_channel_trajectory requires config.jump_selection = :sweep " *
        "(got :$(config.jump_selection)). The :random selection runs a stochastic " *
        "process whose deterministic Φ_δ matvec is e^{δ𝓛} only in expectation."))

    # Domain dispatch for ham_or_trott (mirrors run_thermalize). Computed
    # locally for the Gibbs-reference basis transform below; the channel
    # matvec consumes the same `ham_or_trott` via `ws.ham_or_trott`.
    ham_or_trott = config.domain isa TrotterDomain ? begin
        trotter === nothing && error("TrotterDomain requires a trotter object")
        trotter
    end : hamiltonian

    # qf-po5: route through `Workspace + apply_delta_channel!`. The Workspace
    # constructor builds per-jump K0s/U_residuals/U_coherents via the SAME
    # helpers run_thermalize calls; the matvec then runs the per-jump Lie–Trotter
    # sweep. Single source of truth — bit-identical to run_thermalize :sweep
    # modulo Krylov truncation.
    # qf-qmi.2: optional pre-built Workspace for sweep reuse (see
    # predict_lindbladian_trajectory rationale). Matvec scales as d^3 here,
    # so reuse savings are smaller than for the Lindbladian path but still
    # sizeable when iterating across (β, σ, δ) at fixed (n, ham, jumps).
    ws = if workspace === nothing
        Workspace(config, hamiltonian, jumps; trotter=trotter)
    else
        @assert length(workspace.jumps) == length(jumps)  "workspace jump count mismatch"
        @assert size(workspace.scratch.rho_next, 1) == d  "workspace dim must match rho_0"
        workspace
    end

    CT = Complex{eltype(hamiltonian.eigvals)}

    # Forward channel matvec: applies one full Φ_δ step on `x`, writes to `out`.
    fwd! = let ws = ws, config = config, ham = hamiltonian
        (out::AbstractMatrix, x::AbstractMatrix) -> begin
            rho_in = Matrix{CT}(x)
            apply_delta_channel!(ws, rho_in, config, ham)
            copyto!(out, ws.scratch.rho_next)
            return out
        end
    end

    decomp = _krylov_spectral_decomposition(
        fwd!, Matrix{CT}(rho_0), d;
        krylovdim=krylovdim, tol=tol,
        sort_mode=:channel,
    )

    # rho_inf for the IMPLEMENTED channel is whatever the leading-mu eigenvector
    # collapses to (NOT exactly hamiltonian.gibbs, but O(δ²) close). The
    # decomposition already hermitised + trace-normalised it.
    #
    # On TrotterDomain the dm evolves in the Trotter eigenbasis (Kraus
    # operators + coherent unitaries are built there), so the Gibbs reference
    # for the trace distance must be transformed into that basis. Mirrors
    # the basis-transform `run_thermalize` does at the top of its loop.
    sigma_beta = if config.domain isa TrotterDomain
        @assert trotter !== nothing
        Matrix{CT}(Hermitian(trotter.eigvecs' * hamiltonian.eigvecs *
                              hamiltonian.gibbs *
                              hamiltonian.eigvecs' * trotter.eigvecs))
    else
        Matrix{CT}(hamiltonian.gibbs)
    end

    delta = float(config.delta)
    n_k = length(k_grid)
    t_grid = collect(k_grid) .* delta
    distances = Vector{Float64}(undef, n_k)
    states = save_states ? Vector{Matrix{CT}}(undef, n_k) : Matrix{CT}[]
    rho_k = Matrix{CT}(undef, d, d)
    h = length(decomp.eigenvalues)

    @inbounds for j in 1:n_k
        k = Int(k_grid[j])
        copyto!(rho_k, decomp.rho_inf)
        for i in 1:h
            # Steady state (i=1) is suppressed via c[1] = 0 in the engine,
            # so we can include all modes uniformly.
            mu_pow = decomp.eigenvalues[i]^k
            rho_k .+= (decomp.c[i] * mu_pow) .* decomp.R_modes[i]
        end
        # Defensive Hermitisation.
        @inbounds for jj in 1:d, kc in 1:d
            rho_k[kc, jj] = (rho_k[kc, jj] + conj(rho_k[jj, kc])) / 2
        end
        distances[j] = sum(svdvals(rho_k .- sigma_beta)) / 2
        save_states && (states[j] = copy(rho_k))
    end

    # Convert leading channel eigenvalue gap to Lindbladian units:
    # mu_2 = exp(δ λ_2) ⇒ |1 - μ_2| ≈ δ |λ_2| at small δ.
    spectral_gap = h >= 2 ? abs(1.0 - abs(decomp.eigenvalues[2])) / delta : NaN

    return (
        t              = t_grid,
        distances      = distances,
        rho_final      = copy(rho_k),
        total_matvecs  = decomp.matvec_count,
        all_converged  = decomp.converged,
        states         = states,
        eigenvalues    = decomp.eigenvalues,
        c              = decomp.c,
        spectral_gap   = spectral_gap,
        rho_inf        = decomp.rho_inf,
        R_modes        = decomp.R_modes,
        sigma_beta     = sigma_beta,
        delta_used     = delta,
        k_grid         = collect(k_grid),
    )
end


# ---------------------------------------------------------------------------
# qf-lkb.3 helpers + sweep harness
# ---------------------------------------------------------------------------

"""
    _build_jump_set(ham::HamHam, num_qubits::Integer) -> Vector{JumpOp}

Standard 3n single-site Pauli jump set (X/Y/Z on each site), normalised by
`sqrt(3n)`, paired with their eigenbasis transforms. Mirrors the recipe in
`make_dll_n3_system` (`test/test_helpers.jl`, test-only) so that the library
sweep harness is self-contained without depending on test infrastructure.
"""
function _build_jump_set(ham::HamHam, num_qubits::Integer)
    paulis = ([X], [Y], [Z])
    num_jumps = length(paulis) * num_qubits
    jump_norm = sqrt(num_jumps)
    jumps = JumpOp[]
    for pauli in paulis
        for site in 1:num_qubits
            op = Matrix(pad_term(pauli, num_qubits, site)) ./ jump_norm
            op_eb = ham.eigvecs' * op * ham.eigvecs
            push!(jumps, JumpOp(op, op_eb,
                                op == transpose(op), op == op'))
        end
    end
    return jumps
end

"""
    _make_init_state(init_state::Symbol, d::Integer,
                     sigma_beta::AbstractMatrix, seed::Integer) -> Matrix{ComplexF64}

Build the initial density matrix per `init_state`:

- `:maximally_mixed`  → `I(d) / d` (deterministic; seed ignored).
- `:random_pure`      → Haar-random pure state via `MersenneTwister(seed)`.
- `:thermal_perturbed`→ `σ_β + ε · H_GUE` rescaled to `ε = 1e-3 / d`,
  re-Hermitised + trace-renormalised. Falls back to `:maximally_mixed` if
  the renormalised state has non-positive trace (numerically pathological).
"""
function _make_init_state(init_state::Symbol, d::Integer,
                          sigma_beta::AbstractMatrix, seed::Integer)
    if init_state === :maximally_mixed
        return Matrix{ComplexF64}(I(d) / d)
    elseif init_state === :random_pure
        rng = MersenneTwister(seed)
        psi = randn(rng, ComplexF64, d)
        psi ./= norm(psi)
        return psi * psi'
    elseif init_state === :thermal_perturbed
        rng = MersenneTwister(seed)
        # GUE perturbation: H = (G + G') / 2, G ~ ComplexNormal.
        G = randn(rng, ComplexF64, d, d)
        H = (G + G') / 2
        H_norm = max(norm(H), 1e-30)
        H .*= (1e-3 / d) / H_norm
        rho = Matrix{ComplexF64}(sigma_beta) .+ H
        # Re-Hermitise + renormalise trace.
        rho .= (rho + rho') / 2
        tr_now = real(tr(rho))
        tr_now <= 0 && return Matrix{ComplexF64}(I(d) / d)  # fallback
        rho ./= tr_now
        return rho
    else
        throw(ArgumentError(
            "init_state must be :maximally_mixed, :random_pure, or :thermal_perturbed (got :$init_state)"))
    end
end

"""
    _sweep_sidecar_path(output_dir, n, beta, seed, mode, construction_tag, domain_tag) -> String

Canonical sweep sidecar filename:
`sweep_n<n>_beta<beta>_seed<seed>_<mode>_<construction>_<domain>.bson`. β is
formatted to up to 6 decimals with trailing zeros and dangling decimal points
stripped so integer or simple-decimal values render compactly. The `domain`
suffix prevents BohrDomain and EnergyDomain caches from colliding when both
are exercised against the same `output_dir`.
"""
function _sweep_sidecar_path(output_dir::AbstractString, n::Integer, beta::Real,
                             seed::Integer, mode::Symbol,
                             construction_tag::AbstractString,
                             domain_tag::AbstractString)
    beta_str = let s = @sprintf("%.6f", float(beta))
        s = rstrip(s, '0')
        s = rstrip(s, '.')
        isempty(s) ? "0" : s
    end
    fname = "sweep_n$(n)_beta$(beta_str)_seed$(seed)_$(mode)_$(construction_tag)_$(domain_tag).bson"
    return joinpath(output_dir, fname)
end

"""
    _channel_sweep_sidecar_path(output_dir, n, beta, seed, eps, filter_kind,
                                construction_tag, domain_tag) -> String

Channel-sweep sidecar filename (qf-e4z.2):
`channel_n<n>_beta<beta>_seed<seed>_eps<eps>_<filter>_<construction>_<domain>.bson`.
The extra `eps` and `filter` segments prevent collisions across the (n, β, ε,
filter) grid that `sweep_channel_mixing` walks.
"""
function _channel_sweep_sidecar_path(output_dir::AbstractString, n::Integer,
                                     beta::Real, seed::Integer, eps::Real,
                                     filter_kind::Symbol,
                                     construction_tag::AbstractString,
                                     domain_tag::AbstractString)
    beta_str = let s = @sprintf("%.6f", float(beta))
        s = rstrip(s, '0'); s = rstrip(s, '.')
        isempty(s) ? "0" : s
    end
    eps_str = @sprintf("%.0e", float(eps))            # e.g. "1e-03"
    fname = "channel_n$(n)_beta$(beta_str)_seed$(seed)_eps$(eps_str)_$(filter_kind)_$(construction_tag)_$(domain_tag).bson"
    return joinpath(output_dir, fname)
end

"""
    sweep_mixing_times(n_values, beta_values; kwargs...) -> Vector{NamedTuple}

Drive the matrix-free ODE-integrator + bi-exp extrapolation pipeline
(`integrate_to_gibbs` + `estimate_mixing_time(...; model=:biexp)`) over a
flat product of `(n, β, seed)` triples. Optionally multithreaded across
triples; optionally persisted as BSON sidecars (one per triple) for
cluster-job resumability.

# Arguments
- `n_values`: vector of qubit counts. Each `n` requires
  `<hamiltonian_dir>/heis_disordered_periodic_n<n>.bson` to exist.
- `beta_values`: vector of inverse temperatures.

# Keywords
- `construction = KMS()`: `KMS()` (default) or `DLL()`. For DLL a `filter`
  is required (validated by `validate_config!` per (n, β) point).
- `domain::AbstractDomain = BohrDomain()`: `BohrDomain()` (DLL default)
  or `EnergyDomain()` (CKG production path; qf-lkb.11). DLL is restricted
  to BohrDomain by `validate_config!`.
- `filter`: `nothing` (CKG default) or a `DLLGaussianFilter(β)` /
  `DLLMetropolisFilter(β)` instance. For DLL sweeps, the filter's β must
  match the sweep's β values; we rebuild the filter per `(n, β)` point
  internally to satisfy this — the user-supplied `filter` instance is only
  used as a *type tag* to choose between Gaussian and Metropolis.
- `mode::Symbol = :L`: `:L` for Lindbladian ρ-flow, `:K` for discriminant
  ψ-flow.
- `method::Symbol = :ode` (qf-ev5.7): integration backend.
  `:ode` is the legacy `lindblad_action_integrate` Krylov-exponentiate path
  (~ `t_grid_length` × `krylovdim` matvecs of `apply_lindbladian!`). `:krylov`
  switches to `predict_lindbladian_trajectory` which builds the spectral
  decomposition once (~ `spectral_krylovdim` matvecs total) and evaluates
  `rho(t_k)` in closed form on the grid. The two paths return identical
  NamedTuple shape so the bi-exp pipeline is unchanged. `:krylov` requires
  `mode = :L` (it spectrally decomposes L itself; the discriminant K-flow
  has no analogous closed-form expansion in this codebase).
- `spectral_krylovdim::Integer = 60`: Arnoldi subspace size for `method=:krylov`.
  Larger captures more eigenmodes ⇒ better early-time accuracy; the slow-tail
  regime relevant for τ_mix needs only ~ 30 modes at n ≤ 6.
- `seeds = [42]`: vector of integer seeds. The flat sweep is `(n, β, seed)`.
- `init_state::Symbol = :maximally_mixed`: `:maximally_mixed`,
  `:random_pure`, or `:thermal_perturbed`.
- `a::Real = 0.0`, `s::Real = 0.25`: smooth-Metropolis parameters of the
  CKG transition weight γ(ω). Default a=0, s=0.25 matches the thesis-numerics
  convention (qf-lkb.11). Forwarded to the per-cell `Config`.
- `target_epsilon::Real = 1e-3`: target distance for the bi-exp extrapolation.
- `t_max_factor::Union{Real, Symbol} = :auto`: scales the integration horizon
  via `t_max = t_max_factor / gap_estimate`. With `:auto` (qf-lkb.10), the
  factor is derived from `target_epsilon` as `max(5.0, 1.5 * log10(1/eps))`,
  giving 5 for `eps=1e-3`, 9 for `eps=1e-6`, 14 for `eps=1e-9`. Pass a
  numeric value to override (legacy default was `3.0`).
- `t_grid_length::Integer = 81`: number of grid points in `[0, t_max]`.
- `krylovdim::Integer = 30`, `tol::Real = 1e-10`: forwarded to
  `KrylovKit.exponentiate`.
- `output_dir = nothing`: if non-`nothing`, save BSON sidecars there.
- `hamiltonian_dir`: defaults to `joinpath(dirname(@__DIR__), "hamiltonians")`
  relative to the package source.
- `use_threads::Bool = true`: outer loop via `Threads.@threads`. Each thread
  builds its own `Workspace` internally via `integrate_to_gibbs`, so
  per-thread isolation is automatic; no shared mutable state.
- `skip_existing::Bool = true`: if `output_dir` is set and a sidecar already
  exists for `(n, β, seed)`, load it instead of re-running.

# Returns

`Vector{NamedTuple}` with one entry per `(n, β, seed)` point and fields:
`(n, beta, seed, init_state, mode, construction, domain, filter_name,
gap_est, t_max, n_grid, total_matvecs, all_converged, fitted_gap,
mixing_time, r_squared, converged_fit, wall_time)`.

# Memory model

For CKG and DLL (both matrix-free via `apply_lindbladian!` after qf-lkb.9),
per-thread footprint is dominated by two Krylov subspaces: one for the
spectral-gap pre-pass (`krylovdim ≈ 30`) and one for the time evolution
(`krylov_dim` kwarg). Each subspace is `~k · d²` complex floats. At n=10,
`k = 30`, d²=10⁶ → ~480 MB/thread per subspace, ~960 MB/thread total.
No dense d²×d² Liouvillian is ever constructed, so the cell is feasible
at all n where the per-thread Krylov budget fits. The DLL specialised
`Workspace` additionally caches `|jumps|` per-jump `d × d` Lindblad
operators (`16 · d² · |jumps|` bytes for ComplexF64), cheap relative to
the Krylov subspaces.

# Threading note

BLAS is multi-threaded by default, so a `Threads.@threads`-driven sweep
with `Threads.nthreads()` Julia workers may saturate the CPU twice if
BLAS is not constrained. For cluster runs, set `JULIA_NUM_THREADS=N` and
`BLAS.set_num_threads(1)` externally. The harness does **not** alter
BLAS threading itself.

# Example

```julia
# CKG production sweep using EnergyDomain (qf-lkb.11):
results = sweep_mixing_times([3, 4, 5, 6], [1.0, 5.0, 10.0];
    mode = :L, seeds = [42, 43, 44],
    domain = EnergyDomain(),
    output_dir = "results/ckg_sweep")
```
"""
function sweep_mixing_times(
    n_values::AbstractVector{<:Integer},
    beta_values::AbstractVector{<:Real};
    construction::AbstractConstruction = KMS(),
    domain::AbstractDomain = BohrDomain(),
    filter::Union{Nothing, AbstractFilter} = nothing,
    mode::Symbol = :L,
    method::Symbol = :ode,
    seeds::AbstractVector{<:Integer} = [42],
    init_state::Symbol = :maximally_mixed,
    a::Real = 0.0,
    s::Real = 0.25,
    target_epsilon::Real = 1e-3,
    t_max_factor::Union{Real, Symbol} = :auto,
    t_grid_length::Integer = 81,
    krylovdim::Integer = 30,
    spectral_krylovdim::Integer = 60,
    tol::Real = 1e-10,
    output_dir::Union{Nothing, AbstractString} = nothing,
    hamiltonian_dir::AbstractString = joinpath(dirname(@__DIR__), "hamiltonians"),
    hamiltonian_filename::Function = n -> "heis_disordered_periodic_n$(n).bson",
    use_threads::Bool = true,
    skip_existing::Bool = true,
    # qf-e4z: thread per-cell ideal-Lindbladian recipe from the parameter table
    # (`scripts/output/ideal_lindbladian_param_table.bson`). When `param_table_bson`
    # is non-`nothing` AND `construction isa KMS` (CKG) AND `domain isa EnergyDomain`,
    # the per-cell `(r_D, w0_D, t0_D)` is read from the table row matching
    # `(n, β, target_epsilon, filter_kind)` instead of using the legacy
    # over-provisioned `num_energy_bits = 12, w0 = 0.05` defaults. DLL on
    # BohrDomain is analytical so the table is unused even when supplied.
    #
    # The recipe is calibrated on cells at `n ≤ 6` (qf-7xt). For `n ≥ 8` the
    # K-prefactors are unexplored — spot-check `‖L_E − L_B‖_op / ‖L_B‖_op` at
    # one cell per new (n, β) and bump `_r_D_ideal` in
    # `scripts/numerics_param_table.jl` if the relative error exceeds 10⁻⁹.
    # See `drafts/error-analysis/channel-param-table.md` §"Caveats".
    param_table_bson::Union{Nothing, AbstractString} = nothing,
    filter_kind::Symbol = :smooth_metro,
)::Vector{NamedTuple}
    mode in (:L, :K) || throw(ArgumentError("mode must be :L or :K (got :$mode)"))
    method in (:ode, :krylov) || throw(ArgumentError(
        "method must be :ode or :krylov (got :$method)"))
    method === :krylov && mode === :K && throw(ArgumentError(
        "method=:krylov requires mode=:L (the K-mode discriminant flow needs the K matvec, " *
        "not the L spectral expansion). Run mode=:L if τ_mix is the goal."))
    init_state in (:maximally_mixed, :random_pure, :thermal_perturbed) || throw(
        ArgumentError("init_state ∈ {:maximally_mixed, :random_pure, :thermal_perturbed} (got :$init_state)"))
    domain isa Union{BohrDomain, EnergyDomain} || throw(
        ArgumentError("sweep_mixing_times supports BohrDomain or EnergyDomain only (got $(typeof(domain)))"))
    domain isa EnergyDomain && construction isa DLL && throw(
        ArgumentError("DLL construction is not supported in EnergyDomain (out of scope for DLL-2)"))

    # qf-lkb.10: derive t_max_factor from target_epsilon when :auto.
    # Heuristic: factor = max(5.0, 1.5 * log10(1 / target_eps)) — gives 5 at
    # 1e-3 (preserving legacy behaviour), 9 at 1e-6, 14 at 1e-9. Bi-exp's
    # offset C must drop below ~target_eps/100 for the extrapolation to be
    # well-conditioned; the trajectory needs to reach distance ~target_eps,
    # which takes ~log10(initial_dist/target_eps) decay times.
    t_max_factor_resolved::Float64 = if t_max_factor isa Symbol
        t_max_factor === :auto || throw(ArgumentError(
            "t_max_factor must be :auto or a positive number (got :$t_max_factor)"))
        max(5.0, 1.5 * log10(1.0 / target_epsilon))
    else
        Float64(t_max_factor) > 0 || throw(ArgumentError(
            "t_max_factor must be positive (got $t_max_factor)"))
        Float64(t_max_factor)
    end

    construction_tag = construction isa DLL ? "DLL" : "KMS"
    domain_tag = domain isa EnergyDomain ? "Energy" : "Bohr"
    filter_name = if filter === nothing
        "default"
    elseif filter isa DLLGaussianFilter
        "DLLGaussian"
    elseif filter isa DLLMetropolisFilter
        "DLLMetropolis"
    else
        string(typeof(filter).name.name)
    end

    # Load the ideal-Lindbladian parameter table once when threading is requested.
    # Only consumed by the CKG / EnergyDomain branch below; ignored otherwise.
    use_param_table = (param_table_bson !== nothing
                        && construction isa KMS
                        && domain isa EnergyDomain)
    table_rows = use_param_table ? _load_channel_param_table(param_table_bson) : NamedTuple[]
    if use_param_table
        filter_kind in (:gaussian, :smooth_metro, :kinky_metro) || throw(ArgumentError(
            "param_table_bson requires filter_kind ∈ {:gaussian, :smooth_metro, :kinky_metro} (got :$filter_kind)"))
    end

    if output_dir !== nothing && !isdir(output_dir)
        mkpath(output_dir)
    end

    # Flat product of (n, β, seed). Materialised so @threads can index it.
    points = [(Int(n), Float64(β), Int(s)) for n in n_values
                                            for β in beta_values
                                            for s in seeds]
    n_points = length(points)
    results = Vector{NamedTuple}(undef, n_points)
    skipped = falses(n_points)

    # Pre-pass: load and stamp sidecars for skip_existing. Done serially —
    # strictly before the threaded launcher — to dodge concurrent BSON.load
    # contention and any race with per-thread sidecar writes downstream.
    if output_dir !== nothing && skip_existing
        for i in 1:n_points
            n_i, β_i, seed_i = points[i]
            sidecar = _sweep_sidecar_path(output_dir, n_i, β_i, seed_i,
                                          mode, construction_tag, domain_tag)
            if isfile(sidecar)
                try
                    d_loaded = BSON.load(sidecar, @__MODULE__)
                    results[i] = NamedTuple(d_loaded[:result])
                    skipped[i] = true
                catch err
                    @warn "skip_existing: failed to load sidecar; will recompute" sidecar err
                end
            end
        end
    end

    runner = function (i)
        skipped[i] && return
        n_i, β_i, seed_i = points[i]
        ham_path = joinpath(hamiltonian_dir, hamiltonian_filename(n_i))
        if !isfile(ham_path)
            @warn "Hamiltonian file missing; skipping point" n=n_i β=β_i ham_path
            results[i] = (n=n_i, beta=β_i, seed=seed_i, init_state=init_state,
                          mode=mode, method=method, construction=construction_tag,
                          domain=domain_tag,
                          filter_name=filter_name,
                          target_epsilon=float(target_epsilon),
                          filter_kind=filter_kind,
                          r_D=0, w0_D=NaN, t0_D=NaN,
                          gap_est=NaN, t_max=NaN,
                          t_max_factor=t_max_factor_resolved,
                          tau_mix_bound=NaN, n_grid=0,
                          total_matvecs=0, all_converged=false,
                          fitted_gap=NaN, mixing_time=NaN,
                          mixing_time_source=:nan, r_squared=NaN,
                          converged_fit=false, wall_time=0.0)
            return
        end

        t0_run = time()
        ham = _load_hamiltonian_bson(ham_path, β_i)
        jumps = _build_jump_set(ham, n_i)

        # Per-(n, β) filter: DLL filter must match β; rebuild from the user's
        # filter type tag if construction is DLL. CKG passes filter through
        # unchanged (typically `nothing`).
        local_filter = if construction isa DLL
            if filter isa DLLGaussianFilter
                DLLGaussianFilter(β_i)
            elseif filter isa DLLMetropolisFilter
                DLLMetropolisFilter(β_i)
            else
                throw(ArgumentError(
                    "DLL sweeps require filter::Union{DLLGaussianFilter, DLLMetropolisFilter}"))
            end
        else
            filter
        end

        # Per-cell parameter resolution. Default to the legacy over-provisioned
        # (r_D=12, w0=0.05) — fine for BohrDomain (analytical) and DLL. For CKG
        # on EnergyDomain with `param_table_bson` provided, look up the recipe
        # from `ideal_lindbladian_param_table.bson` and use (r_D, w0_D, t0_D)
        # from the row matching (n_i, β_i, target_epsilon, filter_kind).
        cell_r_D, cell_w0_D, cell_t0_D = 12, 0.05, 2π / (2^12 * 0.05)
        cell_a, cell_s = float(a), float(s)
        cell_with_lc = true
        cell_gauss_params = (nothing, nothing)
        if use_param_table
            row = _lookup_channel_params(table_rows, n_i, β_i, target_epsilon, filter_kind)
            cell_r_D, cell_w0_D, cell_t0_D = row.r_D, row.w0_D, row.t0_D
            cell_a, cell_s = row.a, row.s
            cell_with_lc = row.with_linear_combination
            cell_gauss_params = row.with_linear_combination ? (nothing, nothing) :
                (row.gaussian_omega, row.gaussian_sigma)
        end

        config = Config(
            sim = Lindbladian(),
            domain = domain,
            construction = construction,
            num_qubits = n_i,
            with_linear_combination = cell_with_lc,
            beta = β_i,
            sigma = 1.0 / β_i,
            a = cell_a,
            s = cell_s,
            gaussian_parameters = cell_gauss_params,
            num_energy_bits = cell_r_D,
            w0 = cell_w0_D,
            t0 = cell_t0_D,
            num_trotter_steps_per_t0 = 10,
            filter = local_filter,
        )

        d = size(ham.data, 1)
        rho_0 = _make_init_state(init_state, d,
                                  Matrix{ComplexF64}(ham.gibbs), seed_i)

        # Two paths diverge here (qf-e4y.3):
        # • :krylov — `predict_lindbladian_trajectory` builds a single
        #   bi-orthogonal eigendecomposition; spectral gap and τ_mix(ε)
        #   are read off the eigendecomposition analytically via
        #   `eigenmode_mixing_time`. No second Arnoldi pre-pass, no
        #   bi-exp curve fit. Sidecar emits {gap_est, mixing_time,
        #   mixing_time_source ∈ {:extrapolated, :floor, :nan},
        #   floor_distance}.
        # • :ode — matrix-free ODE integrator over `t_grid`; spectral
        #   gap pre-pass via `krylov_spectral_gap` sets `t_max`, and
        #   `estimate_mixing_time(:biexp)` post-fits the trajectory.
        #   Sidecar emits {gap_est, fitted_gap, r_squared,
        #   converged_fit, mixing_time, mixing_time_source ∈
        #   {:extrapolated, :observed, :nan}}.
        # The two NamedTuple shapes are intentionally different:
        # `:ode` is the legacy / debug path and exposes biexp
        # diagnostics that benchmark harnesses depend on; `:krylov`
        # is the production path the thesis numerics use.

        result = if method === :krylov
            # qf-ev5 / qf-e4y: spectral expansion path. Single Arnoldi
            # factorisation of size `spectral_krylovdim` ⇒ bi-orthogonal
            # eigendecomposition ⇒ closed-form rho(t_k) =
            # rho_inf + Σ c_i exp(lambda_i t_k) R_i. Gap and τ_mix(ε)
            # both come from this single pass.
            #
            # Bracket on the t-axis: pick `t_max = factor / β` as the
            # initial cap; `eigenmode_mixing_time` widens by 3× on
            # demand if the residual hasn't reached `target_epsilon`.
            t_max_seed = t_max_factor_resolved * float(β_i)
            t_grid = collect(range(0.0, t_max_seed, length=t_grid_length))

            predict_res = predict_lindbladian_trajectory(
                config, ham, jumps, rho_0, t_grid;
                krylovdim=spectral_krylovdim, tol=tol,
            )

            # gap_est and t_max are reconstructed from the predictor's own
            # eigendecomposition — `predict_res.spectral_gap` is the smallest
            # |Re(λ_i)| over non-steady modes (already excludes λ_1 ≈ 0).
            gap_est = predict_res.spectral_gap > 0 ?
                       predict_res.spectral_gap : 1.0 / β_i
            t_max = t_max_factor_resolved / max(gap_est, 1e-12)

            res_eig = eigenmode_mixing_time(
                predict_res.eigenvalues, predict_res.c, predict_res.R_modes,
                predict_res.rho_inf, predict_res.sigma_beta,
                target_epsilon;
                t_upper = max(t_max, t_max_seed),
            )
            mixing_time        = res_eig.mixing_time
            mixing_time_source = res_eig.source

            wall = time() - t0_run
            (
                n                   = n_i,
                beta                = β_i,
                seed                = seed_i,
                init_state          = init_state,
                mode                = mode,
                method              = method,
                construction        = construction_tag,
                domain              = domain_tag,
                filter_name         = filter_name,
                target_epsilon      = float(target_epsilon),
                filter_kind         = filter_kind,
                r_D                 = cell_r_D,
                w0_D                = cell_w0_D,
                t0_D                = cell_t0_D,
                gap_est             = gap_est,
                t_max               = t_max,
                t_max_factor        = t_max_factor_resolved,
                tau_mix_bound       = log(size(ham.data, 1) / float(target_epsilon)) /
                                        max(gap_est, 1e-12),
                n_grid              = t_grid_length,
                total_matvecs       = predict_res.total_matvecs,
                all_converged       = predict_res.all_converged,
                mixing_time         = mixing_time,
                mixing_time_source  = mixing_time_source,
                floor_distance      = res_eig.floor_distance,
                wall_time           = wall,
            )
        else  # method === :ode — legacy biexp path
            # Spectral gap estimate via matrix-free Krylov (Arnoldi on
            # apply_lindbladian!). Memory-safe at all n: never materialises
            # the d²×d² Liouvillian (which would OOM beyond n=7); only
            # allocates `krylovdim · d²` complex floats for the subspace.
            # For KMS-DB the Lindbladian gap = parent-Hamiltonian gap H_gap
            # (similarity transform preserves spectrum). On Krylov failure
            # (rare), fall back to the 1/β heuristic so the sweep cell
            # still records data.
            gap_est = try
                krylov_result = krylov_spectral_gap(config, ham, jumps;
                    krylovdim = 30, howmany = 2, tol = 1e-8)
                krylov_result.spectral_gap > 0 ?
                    krylov_result.spectral_gap : 1.0 / β_i
            catch err
                @warn "krylov_spectral_gap failed; falling back to 1/β" n=n_i β=β_i err
                1.0 / β_i
            end
            t_max = t_max_factor_resolved / max(gap_est, 1e-12)
            t_grid = collect(range(0.0, t_max, length=t_grid_length))

            res_int = integrate_to_gibbs(config, ham, jumps, rho_0, t_grid;
                                          mode=mode, krylovdim=krylovdim, tol=tol)
            est = estimate_mixing_time(res_int; model=:biexp,
                                        target_epsilon=target_epsilon,
                                        extrapolate=true)

            # qf-lkb.10: direct-observation fallback. When bi-exp extrapolation
            # fails (offset C close to or above target_eps blocks the inverse
            # solve) but the trajectory itself reached target_eps, return that
            # as τ_mix. The `:observed` source flag distinguishes this from a
            # successful extrapolation.
            mixing_time_source = if isfinite(est.mixing_time)
                :extrapolated
            elseif est.mixing_time_actual !== nothing && isfinite(est.mixing_time_actual)
                :observed
            else
                :nan
            end
            mixing_time = if mixing_time_source === :extrapolated
                est.mixing_time
            elseif mixing_time_source === :observed
                est.mixing_time_actual::Float64
            else
                NaN
            end

            wall = time() - t0_run
            (
                n                   = n_i,
                beta                = β_i,
                seed                = seed_i,
                init_state          = init_state,
                mode                = mode,
                method              = method,
                construction        = construction_tag,
                domain              = domain_tag,
                filter_name         = filter_name,
                target_epsilon      = float(target_epsilon),
                filter_kind         = filter_kind,
                r_D                 = cell_r_D,
                w0_D                = cell_w0_D,
                t0_D                = cell_t0_D,
                gap_est             = gap_est,
                t_max               = t_max,
                t_max_factor        = t_max_factor_resolved,
                tau_mix_bound       = log(size(ham.data, 1) / float(target_epsilon)) /
                                        max(gap_est, 1e-12),
                n_grid              = t_grid_length,
                total_matvecs       = res_int.total_matvecs,
                all_converged       = res_int.all_converged,
                fitted_gap          = est.fitted_gap,
                mixing_time         = mixing_time,
                mixing_time_source  = mixing_time_source,
                r_squared           = est.r_squared,
                converged_fit       = est.converged,
                wall_time           = wall,
            )
        end
        results[i] = result

        if output_dir !== nothing
            sidecar = _sweep_sidecar_path(output_dir, n_i, β_i, seed_i,
                                          mode, construction_tag, domain_tag)
            try
                BSON.bson(sidecar, Dict(:result => Dict(pairs(result))))
            catch err
                @warn "Sidecar write failed (continuing)" sidecar err
            end
        end
        return
    end

    if use_threads
        Threads.@threads for i in 1:n_points
            runner(i)
        end
    else
        for i in 1:n_points
            runner(i)
        end
    end

    return results
end

# ---------------------------------------------------------------------------
# qf-e4z.2 (P0b): sweep_channel_mixing — per-cell channel τ_mix + Ham-sim time
# ---------------------------------------------------------------------------

"""
    _load_channel_param_table(path) -> Vector{NamedTuple}

Load the parameter-table BSON produced by `scripts/numerics_param_table.jl`.
Returns the `:rows` field as a `Vector{NamedTuple}`; throws if the file
is missing or malformed.
"""
function _load_channel_param_table(path::AbstractString)
    isfile(path) || throw(ArgumentError(
        "channel parameter table not found at $path — run `julia --project scripts/numerics_param_table.jl` first"))
    raw = BSON.load(path, @__MODULE__)
    haskey(raw, :rows) || throw(ArgumentError(
        "channel parameter table $path missing :rows entry"))
    return Vector{NamedTuple}(raw[:rows])
end

"""
    _lookup_channel_params(rows, n, β, ε, filter_kind) -> NamedTuple

Find the unique row of the parameter table matching `(n, β, ε, filter_kind)`.
Throws with a clear message if no match (the harness can then skip the cell
or fail loudly depending on the caller's preference).
"""
function _lookup_channel_params(rows::AbstractVector{<:NamedTuple},
                                n::Integer, β::Real, ε::Real, filter_kind::Symbol)
    for r in rows
        r.n == n && r.beta ≈ β && r.eps ≈ ε && r.filter === filter_kind && return r
    end
    throw(ArgumentError(
        "no parameter-table row for (n=$n, β=$β, ε=$ε, filter=$filter_kind). " *
        "Either add the cell to scripts/numerics_param_table.jl or pick a covered cell."))
end

"""
    _build_channel_config(row, n, β, domain, construction; mixing_time_target)
        -> Config{Thermalize, ...}

Construct a `Config{Thermalize, <:Union{TimeDomain, TrotterDomain}, KMS}` from a
parameter-table row. The `mixing_time` field is a placeholder
(`predict_channel_trajectory` only uses `k_grid`); pass any positive number.
"""
function _build_channel_config(row::NamedTuple, n::Integer, β::Real,
                                domain::AbstractDomain,
                                construction::AbstractConstruction;
                                mixing_time_target::Real = 5.0)
    return Config(
        sim = Thermalize(),
        domain = domain,
        construction = construction,
        num_qubits = n,
        beta = β, sigma = row.sigma,
        with_linear_combination = row.with_linear_combination,
        a = row.a, s = row.s,
        gaussian_parameters = row.with_linear_combination ? (nothing, nothing) :
            (row.gaussian_omega, row.gaussian_sigma),
        eta = row.eta,
        num_energy_bits_D = row.r_D, t0_D = row.t0_D, w0_D = row.w0_D,
        num_energy_bits_b_minus = row.r_bm, t0_b_minus = row.t0_bm, w0_b_minus = row.w0_bm,
        num_energy_bits_b_plus  = row.r_bp, t0_b_plus  = row.t0_bp,  w0_b_plus  = row.w0_bp,
        num_trotter_steps_per_t0 = row.M_D,
        delta = row.delta,
        mixing_time = mixing_time_target,
        with_gqsp = row.with_gqsp, gqsp_degree = row.gqsp_degree,
        jump_selection = row.jump_selection,
    )
end

"""
    sweep_channel_mixing(n_values, beta_values; kwargs...) -> Vector{NamedTuple}

Channel-side analogue of [`sweep_mixing_times`](@ref). Drives
[`predict_channel_trajectory`](@ref) over the (n, β, ε, filter, seed) grid,
estimates `τ_mix` via biexp extrapolation, and records the Hamiltonian-
simulation time per cell via [`compute_simulation_time`](@ref). One BSON
sidecar per cell, with `skip_existing` for cluster-job resumability.

Per-cell parameters come from a parameter-table BSON
(`scripts/output/channel_param_table.bson`, generated by P0a) — never
hardcoded inside this harness. To extend the supported (n, β, ε, filter)
grid, edit and re-run `scripts/numerics_param_table.jl`.

# Arguments
- `n_values::AbstractVector{<:Integer}`: qubit counts.
- `beta_values::AbstractVector{<:Real}`: inverse temperatures.

# Keywords
- `target_epsilons = [1e-3]`: list of ε targets per cell.
- `filter_kinds = [:smooth_metro]`: subset of `{:gaussian, :smooth_metro, :kinky_metro}`.
- `domain::AbstractDomain = TrotterDomain()`: faithful δ-channel
  (`TrotterDomain` + GQSP) or analytical ω-OFT (`TimeDomain`).
- `construction::AbstractConstruction = KMS()`: KMS (default) or GNS (no
  coherent term — GQSP rejected by `validate_config!`).
- `param_table_bson::String`: path to the table BSON.
- `family::Symbol = :xxx_zzdis`: Hamiltonian family tag — selects the file
  `<hamiltonian_dir>/heis_<family>_periodic_n<n>.bson` (matches
  `scripts/numerics_param_table.jl`).
- `seeds = [42]`, `init_state = :maximally_mixed`.
- `krylovdim = 30`: Arnoldi subspace size for `predict_channel_trajectory`.
- `k_grid_max_log = 5`: trajectory's `k_max = 10^k_grid_max_log` steps.
- `k_grid_length = 50`: log-spaced k-grid resolution.
- `output_dir = nothing`: if non-`nothing`, write per-cell sidecars there.
- `skip_existing = true`: skip cells whose sidecar already exists.
- `hamiltonian_dir = ".../hamiltonians"`: parent of the family-specific
  Hamiltonian BSONs.

# Per-cell record schema (single BSON layout for all channel sweeps)
```
(
    # identifier
    n, beta, seed, eps, filter, family, construction, domain,
    # parameters (channel only; sourced from the param table BSON)
    r_D, w0_D, t0_D, r_bm, w0_bm, t0_bm, r_bp, w0_bp, t0_bp,
    M_D, M_bm, M_bp, delta, eta,
    with_gqsp, gqsp_degree,
    # results (qf-e4y.4: eigenmode schema)
    tau_mix, tau_mix_source, lambda_gap_channel,
    floor_distance,
    n_steps_to_target, k_max, t_max,
    achieved_dist_at_kmax,
    total_matvecs, all_converged_predict,
    # Hamiltonian simulation time (per `compute_simulation_time`, qf-e4z.18 cost model)
    oft_time_per_step, b_per_be_per_step, b_time_per_step,
    per_step_time, n_steps_total, total_ham_sim_time,
    # bookkeeping
    wall_time_seconds, init_state, family_tag,
)
```

`tau_mix_source ∈ {:extrapolated, :floor, :nan}`:
- `:extrapolated` — eigenmode bisection found a crossing of the trace
  distance with `eps`; `tau_mix` is the continuous-time bisection result.
- `:floor` — `eps` is below the asymptotic floor
  ``\\| \\rho_\\infty - \\sigma_\\beta \\|_1 / 2`` (typically the O(δ)
  channel shift); `tau_mix` is populated with the conservative
  `log(d / eps) / lambda_gap_channel` bound (finite, for plotting).
- `:nan` — degenerate predictor (no slow modes captured); `tau_mix = NaN`.
"""
function sweep_channel_mixing(
    n_values::AbstractVector{<:Integer},
    beta_values::AbstractVector{<:Real};
    target_epsilons::AbstractVector{<:Real} = [1e-3],
    filter_kinds::AbstractVector{Symbol} = [:smooth_metro],
    domain::AbstractDomain = TrotterDomain(),
    construction::AbstractConstruction = KMS(),
    param_table_bson::AbstractString = joinpath(dirname(@__DIR__),
        "scripts", "output", "channel_param_table.bson"),
    family::Symbol = :xxx_zzdisordered,
    seeds::AbstractVector{<:Integer} = [42],
    init_state::Symbol = :maximally_mixed,
    krylovdim::Integer = 30,
    k_grid_max_log::Real = 5,
    k_grid_length::Integer = 50,
    output_dir::Union{Nothing, AbstractString} = nothing,
    skip_existing::Bool = true,
    hamiltonian_dir::AbstractString = joinpath(dirname(@__DIR__), "hamiltonians"),
)::Vector{NamedTuple}
    domain isa Union{TimeDomain, TrotterDomain} || throw(ArgumentError(
        "sweep_channel_mixing supports TimeDomain or TrotterDomain (got $(typeof(domain)))"))
    init_state in (:maximally_mixed, :random_pure, :thermal_perturbed) || throw(
        ArgumentError("init_state ∈ {:maximally_mixed, :random_pure, :thermal_perturbed} (got :$init_state)"))
    all(f -> f in (:gaussian, :smooth_metro, :kinky_metro), filter_kinds) || throw(ArgumentError(
        "filter_kinds must be ⊆ {:gaussian, :smooth_metro, :kinky_metro} (got $filter_kinds)"))
    Float64(k_grid_max_log) > 0 || throw(ArgumentError(
        "k_grid_max_log must be > 0 (got $k_grid_max_log)"))
    k_grid_length >= 10 || throw(ArgumentError(
        "k_grid_length must be ≥ 10 (got $k_grid_length)"))

    construction_tag = construction isa KMS ? "KMS" : construction isa GNS ? "GNS" : "DLL"
    domain_tag = domain isa TrotterDomain ? "Trotter" : "Time"
    family_str = string(family)
    ham_filename = (n) -> "heis_$(family_str)_periodic_n$(n).bson"

    output_dir !== nothing && !isdir(output_dir) && mkpath(output_dir)

    rows_table = _load_channel_param_table(param_table_bson)

    # Flat product of (n, β, ε, filter, seed). Materialised so loops are stable.
    points = [(Int(n), Float64(β), Float64(ε), f, Int(s))
              for n in n_values for β in beta_values
              for ε in target_epsilons for f in filter_kinds for s in seeds]
    n_points = length(points)
    results = Vector{NamedTuple}(undef, n_points)
    skipped = falses(n_points)

    # Pre-pass: load existing sidecars under skip_existing, before the main loop.
    if output_dir !== nothing && skip_existing
        for i in 1:n_points
            n_i, β_i, ε_i, f_i, seed_i = points[i]
            sidecar = _channel_sweep_sidecar_path(output_dir, n_i, β_i, seed_i, ε_i,
                                                   f_i, construction_tag, domain_tag)
            if isfile(sidecar)
                try
                    d_loaded = BSON.load(sidecar, @__MODULE__)
                    results[i] = NamedTuple(d_loaded[:result])
                    skipped[i] = true
                catch err
                    @warn "skip_existing: failed to load channel sidecar; will recompute" sidecar err
                end
            end
        end
    end

    k_max = round(Int, exp10(Float64(k_grid_max_log)))

    for i in 1:n_points
        skipped[i] && continue
        n_i, β_i, ε_i, f_i, seed_i = points[i]
        ham_path = joinpath(hamiltonian_dir, ham_filename(n_i))
        if !isfile(ham_path)
            @warn "Hamiltonian file missing; skipping channel cell" n=n_i β=β_i ham_path
            continue
        end

        local row, ham, jumps, cfg, trotter, rho_0, rho_init, predict_res, res_eig, sim_budget
        t0_run = time()
        try
            row = _lookup_channel_params(rows_table, n_i, β_i, ε_i, f_i)
            ham = _load_hamiltonian_bson(ham_path, β_i)
            cfg = _build_channel_config(row, n_i, β_i, domain, construction)

            # Trotter cache (TrotterDomain only). Use make_trotter_for_config which
            # picks the qf-d0w shared-δt₀ scheme for KMS coherent. May fail at
            # ε=1e-6 if the per-term grids are not commensurate; surfaced below.
            trotter = domain isa TrotterDomain ? make_trotter_for_config(ham, cfg) : nothing

            # Build jump set in the right basis. TimeDomain uses Hamiltonian
            # eigenbasis; TrotterDomain uses the Trotter eigenbasis (Kraus path).
            jumps = if domain isa TrotterDomain
                _jumps_in_basis(ham, n_i, trotter.eigvecs)
            else
                _jumps_in_basis(ham, n_i, ham.eigvecs)
            end

            # Initial state in the basis the predictor expects.
            d = size(ham.data, 1)
            rho_0 = _make_init_state(init_state, d, Matrix{ComplexF64}(ham.gibbs), seed_i)
            rho_init = if domain isa TrotterDomain
                Matrix{ComplexF64}(trotter.eigvecs' * rho_0 * trotter.eigvecs)
            else
                Matrix{ComplexF64}(rho_0)
            end

            k_grid = unique(round.(Int, exp10.(range(0, Float64(k_grid_max_log),
                                                       length=k_grid_length))))

            predict_res = predict_channel_trajectory(cfg, ham, jumps, rho_init, k_grid;
                krylovdim=krylovdim, trotter=trotter)

            # qf-e4y.4: eigenmode τ_mix on the continuous-time axis. Convert
            # channel μ_i → λ_eff_i = log(μ_i) / δ. Steady mode μ_1 ≈ 1 maps
            # to λ_eff_1 ≈ 0 (within numerical drift), which the helper
            # treats as steady via `eigenvalue_zero_tol`. Complex μ near 1
            # have small principal-branch log; far from the branch cut at
            # arg = π for any reasonable δ.
            delta_used = predict_res.delta_used
            lambda_eff = log.(predict_res.eigenvalues) ./ delta_used
            # Bisection upper bracket: take the larger of the trajectory's
            # k_max·δ horizon and a generous gap-based estimate
            # `5 · log(d/ε) / λ_gap`. The trajectory's k_max is set to
            # observe the decay; the bisection needs to bracket τ_mix(ε)
            # which can lie past the observation window when ε is near
            # the channel's asymptotic floor.
            d_dim = size(ham.data, 1)
            gap_ch = predict_res.spectral_gap
            t_upper_ch = if isfinite(gap_ch) && gap_ch > 0
                max(predict_res.t[end], 5.0 * log(d_dim / ε_i) / gap_ch)
            else
                predict_res.t[end]
            end
            res_eig = eigenmode_mixing_time(
                lambda_eff, predict_res.c, predict_res.R_modes,
                predict_res.rho_inf, predict_res.sigma_beta, ε_i;
                t_upper = t_upper_ch,
            )
            # Pass a FINITE τ to compute_simulation_time. The eigenmode helper
            # returns Inf on :floor; substitute the conservative log(d/ε)/λ
            # bound (which is what the prior `:gap` branch used).
            tau_for_budget = let
                if res_eig.source === :extrapolated &&
                   isfinite(res_eig.mixing_time) && res_eig.mixing_time > 0
                    res_eig.mixing_time
                elseif isfinite(gap_ch) && gap_ch > 0
                    log(d_dim / ε_i) / gap_ch
                elseif isfinite(predict_res.t[end]) && predict_res.t[end] > 0
                    predict_res.t[end]
                else
                    1.0
                end
            end
            sim_budget = compute_simulation_time(cfg, ham, tau_for_budget)
        catch err
            @warn "channel cell failed; recording NaN row" n=n_i β=β_i eps=ε_i filter=f_i err
            results[i] = (
                n=n_i, beta=β_i, seed=seed_i, eps=ε_i, filter=f_i,
                family=family_str, construction=construction_tag, domain=domain_tag,
                r_D=NaN, w0_D=NaN, t0_D=NaN, r_bm=NaN, w0_bm=NaN, t0_bm=NaN,
                r_bp=NaN, w0_bp=NaN, t0_bp=NaN,
                M_D=NaN, M_bm=NaN, M_bp=NaN, delta=NaN, eta=NaN,
                with_gqsp=false, gqsp_degree=0,
                tau_mix=NaN, tau_mix_source=:nan, lambda_gap_channel=NaN,
                floor_distance=NaN,
                n_steps_to_target=0, k_max=k_max, t_max=NaN,
                achieved_dist_at_kmax=NaN,
                total_matvecs=0, all_converged_predict=false,
                oft_time_per_step=NaN, b_per_be_per_step=NaN, b_time_per_step=NaN,
                per_step_time=NaN, n_steps_total=0, total_ham_sim_time=NaN,
                wall_time_seconds=time() - t0_run,
                init_state=init_state, family_tag=family_str,
            )
            continue
        end

        # Resolve τ_mix from the eigenmode helper output:
        # - :extrapolated → bisection found a crossing; tau_mix is its result.
        # - :floor → ε below the channel shift `‖ρ_∞ - σ_β‖_1 / 2`; tau_mix
        #   is the conservative `log(d/ε) / λ` bound (matches the prior
        #   `:gap` branch's value), so plotting tooling has a finite number
        #   even when no crossing exists.
        # - :nan → degenerate input (no slow mode captured); tau_mix = NaN.
        tau_mix_source = res_eig.source
        tau_mix = if tau_mix_source === :extrapolated
            res_eig.mixing_time
        elseif tau_mix_source === :floor
            (isfinite(predict_res.spectral_gap) && predict_res.spectral_gap > 0) ?
                log(size(ham.data, 1) / ε_i) / predict_res.spectral_gap : NaN
        else
            NaN
        end
        # Actual k count to reach ε from the trajectory (for sanity-check tooling).
        n_steps_to_target = let dists = predict_res.distances
            idx = findfirst(d -> d <= ε_i, dists)
            idx === nothing ? 0 : Int(idx)
        end

        wall = time() - t0_run
        result = (
            n                       = n_i,
            beta                    = β_i,
            seed                    = seed_i,
            eps                     = ε_i,
            filter                  = f_i,
            family                  = family_str,
            construction            = construction_tag,
            domain                  = domain_tag,
            # parameters (echoed from the param-table row)
            r_D = row.r_D,  w0_D = row.w0_D,  t0_D = row.t0_D,
            r_bm = row.r_bm, w0_bm = row.w0_bm, t0_bm = row.t0_bm,
            r_bp = row.r_bp, w0_bp = row.w0_bp, t0_bp = row.t0_bp,
            M_D = row.M_D, M_bm = row.M_bm, M_bp = row.M_bp,
            delta = row.delta, eta = row.eta,
            with_gqsp = row.with_gqsp, gqsp_degree = row.gqsp_degree,
            # results (qf-e4y.4: eigenmode schema)
            tau_mix                 = tau_mix,
            tau_mix_source          = tau_mix_source,
            lambda_gap_channel      = predict_res.spectral_gap,
            floor_distance          = res_eig.floor_distance,
            n_steps_to_target       = n_steps_to_target,
            k_max                   = k_max,
            t_max                   = predict_res.t[end],
            achieved_dist_at_kmax   = predict_res.distances[end],
            total_matvecs           = predict_res.total_matvecs,
            all_converged_predict   = predict_res.all_converged,
            # Hamiltonian-simulation time
            oft_time_per_step       = sim_budget.oft_time,
            b_per_be_per_step       = sim_budget.b_per_be,
            b_time_per_step         = sim_budget.b_time,
            per_step_time           = sim_budget.per_step_time,
            n_steps_total           = sim_budget.n_steps,
            total_ham_sim_time      = sim_budget.total_time,
            # bookkeeping
            wall_time_seconds       = wall,
            init_state              = init_state,
            family_tag              = family_str,
        )
        results[i] = result

        if output_dir !== nothing
            sidecar = _channel_sweep_sidecar_path(output_dir, n_i, β_i, seed_i, ε_i,
                                                   f_i, construction_tag, domain_tag)
            try
                BSON.bson(sidecar, Dict(:result => Dict(pairs(result))))
            catch err
                @warn "Channel sidecar write failed (continuing)" sidecar err
            end
        end
    end

    return results
end

# Helper: build the standard CKG Pauli jump set in a given basis.
function _jumps_in_basis(ham::HamHam, num_qubits::Integer,
                          basis_eigvecs::AbstractMatrix)
    jumps = JumpOp[]
    jump_norm = sqrt(3 * num_qubits)
    for pauli in (X, Y, Z), site in 1:num_qubits
        op = Matrix(pad_term([pauli], num_qubits, site)) ./ jump_norm
        op_eb = basis_eigvecs' * op * basis_eigvecs
        push!(jumps, JumpOp(op, op_eb, op == transpose(op), op == op'))
    end
    return jumps
end
