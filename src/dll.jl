# DLL Lindblad operators (Ding–Li–Lin 2024, Sec. 3).
#
# The DLL construction collapses the CKG outer ω-loop into a single Lindblad
# operator per coupling A^a (Eq. 3.4 / Remark 12). Two equivalent expressions:
#
#   Bohr-domain (Eq. 3.4, first form):
#     L_a = Σ_{ν ∈ B_H} q^a(ν) e^{-βν/4} A^a_ν
#         = Σ_ν f̂^a(ν) A^a_ν
#
#   Time-domain (Eq. 3.4, third form, OFT at ω = 0):
#     L_a = ∫_{-∞}^∞ f^a(t) e^{i H t} A^a e^{-i H t} dt
#         ≈ Σ_{m=0}^{2M-1} f^a(t_m) e^{i H t_m} A^a e^{-i H t_m} τ
#       on the uniform grid t_m = -Mτ + mτ (Eq. 3.13).
#
# `freq_kernel(::DLLGaussianFilter, ν)` already returns the FULL DLL weighting
# `q(ν) e^{-βν/4}` (Eq. 3.22) — no extra prefactor is applied here. Likewise,
# `time_kernel(::DLLGaussianFilter, t)` already includes the
# `e^{1/8} √(2/π) / β` prefactor of the inverse FT (Eq. 3.3 / 3.22) — no extra
# normalisation either.
#
# Validated end-to-end by `scripts/scratch_dll_dissipator.jl`: Bohr ↔ Time
# operator-norm residual ≤ 4.4e-16, Bohr Gibbs fixed point and trace
# preservation at machine precision.

"""
    dll_lindblad_op_bohr(jump, hamiltonian, filter) -> Matrix

DLL Bohr-domain Lindblad operator (Ding–Li–Lin 2024, Eq. 3.4, first form):

    L_a[i, j] = freq_kernel(filter, λ_i − λ_j) · A_eb[i, j]

with `A_eb = jump.in_eigenbasis` and `λ` the Hamiltonian eigenvalues. Returned
matrix is in the Hamiltonian eigenbasis (the basis in which `jump.in_eigenbasis`
is expressed).

For `DLLGaussianFilter`, `freq_kernel(filter, ν) = q(ν) e^{-βν/4}` is the FULL
DLL weighting (Eq. 3.22) — no extra `e^{-βν/4}` factor is applied here.
"""
function dll_lindblad_op_bohr(
    jump::JumpOp,
    hamiltonian::HamHam{T},
    filter::AbstractFilter,
) where {T<:AbstractFloat}
    eigvals = hamiltonian.eigvals
    A_eb = jump.in_eigenbasis
    n = length(eigvals)
    @assert size(A_eb) == (n, n)
    CT = Complex{T}
    L = zeros(CT, n, n)
    @inbounds for j in 1:n, i in 1:n
        ν_ij = eigvals[i] - eigvals[j]
        L[i, j] = freq_kernel(filter, ν_ij) * A_eb[i, j]
    end
    return L
end

"""
    dll_lindblad_op_time(jump, hamiltonian, time_labels, filter, t0) -> Matrix

DLL Time-domain Lindblad operator via discrete OFT at ω = 0 (Ding–Li–Lin 2024,
Eq. 3.4 third form, Eq. 3.13):

    L_a = Σ_m f(t_m) · D(t_m) · A_eb · D(t_m)† · τ

with `D(t) = Diagonal(exp(i λ_k t))`, `t_m = time_labels[m]`, and `τ = t0` the
uniform grid spacing (the `t0` field of the simulator config matches the
trapezoidal weight here).

`time_kernel(::DLLGaussianFilter, t)` returns the FULL `f(t)` of Eq. 3.22 —
including the `e^{1/8} √(2/π) / β` prefactor — so no extra normalisation is
applied. Returned matrix is in the Hamiltonian eigenbasis.

# Production path note

This explicit Riemann-sum form is the canonical Eq. 3.4 third-form computation
and is used by tests for direct verification. The simulator's production path
(`Lindbladian` + `TimeDomain` + `DLL`) amortises the same computation via the
precomputed `oft_nufft_at_zero` slice in `_precompute_data`, hence
`_jump_contribution!` does not call this function on the hot path.
"""
function dll_lindblad_op_time(
    jump::JumpOp,
    hamiltonian::HamHam{T},
    time_labels::AbstractVector{<:Real},
    filter::AbstractFilter,
    t0::Real,
) where {T<:AbstractFloat}
    eigvals = hamiltonian.eigvals
    A_eb = jump.in_eigenbasis
    n = length(eigvals)
    @assert size(A_eb) == (n, n)
    CT = Complex{T}
    L = zeros(CT, n, n)
    @inbounds for t in time_labels
        ft = time_kernel(filter, t)
        weight = ft * t0
        for j in 1:n, i in 1:n
            phase = cis((eigvals[i] - eigvals[j]) * t)
            L[i, j] += weight * phase * A_eb[i, j]
        end
    end
    return L
end

# ---------------------------------------------------------------------------
# DLL coherent operator G (Ding–Li–Lin 2024, Eqs. 3.5–3.7)
# ---------------------------------------------------------------------------
#
# Frequency-domain kernel (Eq. 3.5, factored through `freq_kernel`):
#     ĝ^a(ν, ν') = (1/2i) · tanh(β(ν'-ν)/4) · freq_kernel(filter, ν) · conj(freq_kernel(filter, ν'))
#
# Bohr form (Eq. 3.7, first equality):
#     G = Σ_a Σ_{ν, ν' ∈ B_H} ĝ^a(ν, ν') (A^a_{ν'})† A^a_ν
#
# Time form (Eq. 3.7, second equality, with the operator-ordering correction):
#     G = Σ_a ∬ g^a(t, t') · A^a(t') · A^a(t) dt dt'
# where `g^a(t, t')` is the inverse 2D FT of `ĝ^a(ν, ν')` (Eq. 3.6) with the
# asymmetric sign convention `e^{-iνt + iν't'}`. **NOTE on Eq. 3.7**: the
# paper's printed third equality reads `A^a(t) A^a(t')` (left/right) but
# substituting Eq. 3.6 into the integral and using the Bohr decomposition
# `A^a(t) = Σ_ν e^{iνt} A^a_ν` yields the order `A^a(t') A^a(t)` (paired-
# with-`t'` on the LEFT). Verified numerically to ~1e-10 op-norm vs the Bohr
# form for β ∈ {1, 5, 10}; see `scripts/scratch_dll_coherent.jl`.
#
# Unlike the CKG coherent term, the DLL `g^a(t, t')` does NOT factorise as
# `b_-(t) · b_+(t')`: the `tanh(β(ν'-ν)/4)` term entangles ν and ν', so the
# 2D inverse FT is not a tensor product. We therefore tabulate `g^a` on a
# uniform `(ν, ν')` grid and 2D-DFT to the time grid directly.

"""
    dll_coherent_kernel_bohr(filter, ν, νp) -> Complex

DLL frequency-domain coherent kernel `ĝ^a(ν, ν')` from Ding–Li–Lin 2024
Eq. 3.5, factored through `freq_kernel`:

    ĝ^a(ν, ν') = (1/2i) · tanh(β(ν'-ν)/4) · freq_kernel(filter, ν) · conj(freq_kernel(filter, ν'))

Defined on any DLL filter (`DLLGaussianFilter` or `DLLMetropolisFilter`)
that exposes a `.beta` field and a `freq_kernel`. The kernel always
factors as a `freq_kernel(ν) · conj(freq_kernel(ν'))` outer product times
a `tanh` weight that mixes (ν, ν'); the difference between filters is
the shape of `freq_kernel`.
"""
@inline function dll_coherent_kernel_bohr(
    filter::Union{DLLGaussianFilter{T}, DLLMetropolisFilter{T}},
    ν::Real,
    νp::Real,
) where {T<:AbstractFloat}
    β = filter.beta
    pref = one(T) / (2im)
    th = tanh(β * (νp - ν) / 4)
    fkν = freq_kernel(filter, ν)
    fkνp = freq_kernel(filter, νp)
    return Complex{T}(pref) * Complex{T}(th) * fkν * conj(fkνp)
end

"""
    dll_coherent_op_bohr(jumps, hamiltonian, filter, beta) -> Matrix

DLL Bohr-domain coherent operator `G` (Ding–Li–Lin 2024, Eqs. 3.5 + 3.7 first
form), evaluated via the closed-form `O(n³)` expression that is algebraically
equivalent to the `O(|B_H|² n³)` Bohr double sum:

    G = (1/2i) · T ⊙ Σ_a (M^a)† M^a
        T[i, j]    = tanh(β(λ_j − λ_i)/4)                    (real, anti-Hermitian)
        M^a[k, j]  = freq_kernel(filter, λ_k − λ_j) · A^a_eb[k, j]

where `⊙` is the Hadamard product. Derivation: in the eigenbasis,
`A^a_ν[k, j]` is supported on `(k, j)` with `ν = λ_k − λ_j`, so the (ν, ν')
sum in Eq. 3.7 collapses on each `k`, and the `tanh(β(ν'−ν)/4)` factor in
Eq. 3.5 evaluates to `tanh(β(λ_j − λ_i)/4)` independent of `k`. Verified to
machine precision against the original brute-force form in
`scripts/scratch_dll_coherent_closedform.jl` (residual `4.5e-17`, relative
error `8.2e-16` at `n=5, β=5`).

Self-adjointness (Theorem 10): `T` is real anti-Hermitian and `Σ_a (M^a)† M^a`
is Hermitian, so the Hadamard product is anti-Hermitian and the `1/(2i)` factor
makes `G` Hermitian to roundoff.

The `beta` argument is accepted for symmetry with `dll_coherent_op_time`
but the kernel reads `β` from the filter; mismatch is the caller's
responsibility (validated upstream by `validate_config!`).
"""
function dll_coherent_op_bohr(
    jumps::AbstractVector{<:JumpOp},
    hamiltonian::HamHam{T},
    filter::AbstractFilter,
    beta::Real,
) where {T<:AbstractFloat}
    eigvals = hamiltonian.eigvals
    n = length(eigvals)
    CT = Complex{T}

    # f̂(λ_k − λ_j) on a square grid; both supported filters return real-valued
    # `freq_kernel` (Gaussian and DLL Gevrey), so this is `T`-typed.
    f_mat = Matrix{T}(undef, n, n)
    @inbounds for j in 1:n, k in 1:n
        f_mat[k, j] = T(freq_kernel(filter, eigvals[k] - eigvals[j]))
    end

    # tanh(β(λ_j − λ_i)/4): real, anti-symmetric (T[j, i] = −T[i, j]).
    tanh_mat = Matrix{T}(undef, n, n)
    βT = T(beta)
    @inbounds for j in 1:n, i in 1:n
        tanh_mat[i, j] = tanh(βT * (eigvals[j] - eigvals[i]) / 4)
    end

    G = zeros(CT, n, n)
    M = Matrix{CT}(undef, n, n)
    MdM = Matrix{CT}(undef, n, n)
    @inbounds for jump in jumps
        @. M = f_mat * jump.in_eigenbasis
        mul!(MdM, M', M)
        @. G += MdM
    end
    pref = CT(1) / (2im)
    @. G = pref * tanh_mat * G
    return G
end

# ---------------------------------------------------------------------------
# Closed-form `g(t, t')` for the DLL Gaussian-type filter (Phase C / qf-hur.3)
# ---------------------------------------------------------------------------
#
# For `f̂(ν) = exp(-(βν+1)²/8) · exp(1/8) = exp(-β²ν²/8 - βν/4)` (Eq. 3.22),
# substituting (u, v) = (ν' − ν, (ν + ν')/2) into the 2D inverse FT of
# `ĝ(ν, ν') = (1/2i) tanh(β(ν'−ν)/4) f̂(ν) f̂(ν')` yields a separable form:
#
#     g(t, t') = (1 / (8π²)) · I_v(t' − t) · J(t + t')
#         I_v(δ) = (2√π / β) · exp(1/4 − iδ/β − δ²/β²)              (closed form)
#         J(s)   = 2 ∫_0^∞ tanh(βu/4) · exp(-β²u²/16) · sin(us/2) du   (1D quad)
#
# Derivation: f̂(v − u/2) f̂(v + u/2) = exp(-β²v²/4 − β²u²/16 − βv/2), so the
# `(u, v)` integrals decouple; the v-integral is a pure Gaussian (closed form),
# the u-integral has tanh times a Gaussian (1D quadrature).
#
# This eliminates the internal `(ν, ν')`-tabulation grid in the previous
# `dll_coherent_op_time`. The Eq. 3.17 quadrature error structure becomes
# exactly the time-grid Riemann sum on `g(t_m, t_n)` with no frequency-domain
# discretisation. Verified to machine precision against an Nν=1024 reference
# in `scripts/scratch_dll_coherent_v2.jl` (`‖g_cf − g_finegrid‖_op = 5e-14`).

# 1D quadrature for `J(s)`. The integrand decays as `exp(-(βu/4)²)` so
# truncate at `u_max = 24/β` (`exp(-36) ≈ 2e-16`). At large `s`, `J(s)` is
# tiny but oscillatory; the absolute tolerance prevents QuadGK from overworking
# the trivially-small values.
function _dll_J_quadrature(beta::T, s::Real; rtol::Real=1e-12, atol::Real=1e-14) where {T<:AbstractFloat}
    integrand(u) = tanh(beta * u / 4) * exp(-(beta * u)^2 / 16) * sin(u * s / 2)
    cutoff = T(24) / beta
    val, _ = quadgk(integrand, T(0), cutoff; rtol=rtol, atol=atol)
    return T(2) * val
end

# Closed-form `g(t, t') = (1 / (8π²)) · I_v(t' − t) · J(t + t')`.
@inline function _dll_g_closed_form(beta::T, t::Real, tp::Real, J_val::Real) where {T<:AbstractFloat}
    δ = T(tp) - T(t)
    factor_diff = (T(2) * sqrt(T(π)) / beta) * exp(T(1) / 4 - im * δ / beta - δ^2 / beta^2)
    # Combined prefactor: 1/(2π)² · 1/2 = 1/(8π²).
    pref = inv(T(8) * T(π)^2)
    return Complex{T}(pref) * factor_diff * J_val
end

# Tabulate `J(s)` on the unique `t_m + t_n` values in `time_labels`. Returns
# the J value vector and a `Dict` mapping each unique sum to its index.
function _dll_J_table(beta::T, time_labels::AbstractVector{<:Real}; rtol::Real=1e-12, atol::Real=1e-14) where {T<:AbstractFloat}
    Nt = length(time_labels)
    sums_set = Set{T}()
    @inbounds for nidx in 1:Nt, m in 1:Nt
        push!(sums_set, T(time_labels[m]) + T(time_labels[nidx]))
    end
    sums_vec = sort!(collect(sums_set))
    J_vals = Vector{T}(undef, length(sums_vec))
    @inbounds for (i, s) in enumerate(sums_vec)
        J_vals[i] = _dll_J_quadrature(beta, s; rtol=rtol, atol=atol)
    end
    sum_to_index = Dict{T, Int}(s => i for (i, s) in enumerate(sums_vec))
    return J_vals, sum_to_index
end

"""
    dll_coherent_op_time(jumps, hamiltonian, time_labels, filter, beta, τ) -> Matrix

DLL Time-domain coherent operator `G` (Ding–Li–Lin 2024, Eqs. 3.6 + 3.7,
second equality) on the simulator's uniform time grid:

    G ≈ Σ_a Σ_{m, n} g(t_m, t_n) · A^a(t_n) · A^a(t_m) · τ²

with `A^a(t) = D(t) A^a D(t)†`, `D(t) = Diagonal(exp(i λ_k t))`, `t_m =
time_labels[m]`, and `τ` the trapezoidal weight (typically `config.t0`).

**Operator-ordering correction**: the paper's printed Eq. 3.7 third equality
reads `A^a(t) A^a(t')` (left/right) but re-deriving from Eq. 3.6 + the Bohr
decomposition gives `A^a(t') A^a(t)` — the operator paired with the `t'`
argument of `g` multiplies on the LEFT. Verified numerically against the Bohr
form (see `scripts/scratch_dll_coherent.jl`).

Implementation (Phase C / qf-hur.3):
1. **Closed-form `g(t_m, t_n)`** via the `(u, v)` substitution that makes the
   DLL filter separable: `g(t, t') = (1/(8π²)) · I_v(t' − t) · J(t + t')`,
   where `I_v` is closed-form Gaussian and `J(s)` is a 1D QuadGK integral
   tabulated at the `2Nt − 1` unique `t_m + t_n` values. Eliminates the
   internal `(ν, ν')` grid (a Julia-numerics artefact, not in the paper).
2. **NUFFT-factored Riemann sum**: the per-jump cost
   `Σ_{m, n} g(t_m, t_n) · cis((λ_i − λ_k) t_n + (λ_k − λ_j) t_m) · τ²`
   is a 2D type-3 NUFFT evaluated at the `n³` Bohr-target tuples
   `(λ_k − λ_j, λ_i − λ_k)` for `(i, j, k) ∈ [1, n]³`. The result `Q_ijk` is
   jump-independent — built once per Liouvillian, reused across all jumps.

The Eq. 3.17 error structure now matches the paper: it is the time-grid
Riemann sum on `g(t_m, t_n)` with no frequency-domain discretisation.

The `beta` argument must equal `filter.beta` (validated upstream by
`validate_config!`).
"""
function dll_coherent_op_time(
    jumps::AbstractVector{<:JumpOp},
    hamiltonian::HamHam{T},
    time_labels::AbstractVector{<:Real},
    filter::DLLGaussianFilter,
    beta::Real,
    τ::Real,
) where {T<:AbstractFloat}
    Nt = length(time_labels)
    CT = Complex{T}
    βT = T(beta)

    # ------------------------------------------------------------------
    # Step 1 + 2: closed-form g_tt[m, n] (Gaussian-specific).
    # ------------------------------------------------------------------
    J_vals, sum_to_index = _dll_J_table(βT, time_labels)
    g_tt = Matrix{CT}(undef, Nt, Nt)
    @inbounds for nidx in 1:Nt, m in 1:Nt
        s = T(time_labels[m]) + T(time_labels[nidx])
        J_val = J_vals[sum_to_index[s]]
        g_tt[m, nidx] = _dll_g_closed_form(βT, time_labels[m], time_labels[nidx], J_val)
    end

    # ------------------------------------------------------------------
    # Steps 3 + 4: filter-agnostic NUFFT contraction g_tt → Q_ijk → G.
    # ------------------------------------------------------------------
    return _dll_coherent_from_g_tt(jumps, hamiltonian, g_tt, time_labels, τ)
end

"""
    dll_coherent_op_time(jumps, hamiltonian, time_labels, filter::DLLMetropolisFilter,
                         beta, τ; nu_grid_size=256) -> Matrix

DLL Time-domain coherent operator `G` (Eq. 3.7) for the Metropolis-type
filter. Mirrors the Gaussian path's structure but computes `g(t, t')`
numerically because the Metropolis `f̂(ν) = q(ν) e^{-βν/4}` (with
`q(ν) = exp(-√(1+(βν)²)/4) · w(ν/S)`) does **not** separate under the
`(u, v)` substitution that yields the Gaussian closed form.

Implementation:
1. **Tabulate `ĝ(ν, ν')` on a uniform `(ν, ν')` grid on `[-S, S]²`** —
   the bump's compact support makes this an *exact* truncation of an
   otherwise-infinite integral. Smoothness of the integrand at the
   boundaries (Hörmander bump) means uniform-trapezoidal quadrature
   converges super-polynomially in `nu_grid_size`; the default
   `nu_grid_size = 256` gives <1e-12 quadrature error.
2. **2D type-3 NUFFT to the time grid** — sources are the
   `(ν_p, ν'_q)` tensor grid (`Nν²` source points) with strengths
   `ĝ(ν_p, ν'_q) · (Δν)² / (2π)²`; targets are the
   `(t_m, t_n)` tensor grid (`Nt²` target points). The sign convention
   `e^{-iνt + iν't'}` (Eq. 3.6) is realised by feeding the source-x
   coordinates as `-ν_p` (FINUFFT's `isign = +1` then gives `e^{+i s_x t_x}`,
   so `e^{-i ν_p t_m}`).
3. **Filter-agnostic contraction** `g_tt → Q_ijk → G` via the shared
   `_dll_coherent_from_g_tt` helper (qf-wmg.6 refactor).

Performance: `Nν² + Nt²` for the NUFFT, plus the Step 3+4 helper cost.
At `n = 3`, `Nt = 4096`, `Nν = 256`, expected ~50 ms total.

The `beta` argument must equal `filter.beta` (validated upstream by
`validate_config!`).
"""
function dll_coherent_op_time(
    jumps::AbstractVector{<:JumpOp},
    hamiltonian::HamHam{T},
    time_labels::AbstractVector{<:Real},
    filter::DLLMetropolisFilter{T},
    beta::Real,
    τ::Real;
    nu_grid_size::Int = 256,
) where {T<:AbstractFloat}
    Nt = length(time_labels)
    CT = Complex{T}
    βT = T(beta)
    S = filter.S

    # ------------------------------------------------------------------
    # Step 1: tabulate ĝ(ν, ν') on a uniform grid over [-S, S]².
    # ĝ(ν, ν') = (1/2i) · tanh(β(ν'-ν)/4) · f̂(ν) · conj(f̂(ν'))
    # f̂ is real for the Metropolis filter (q is real, e^{-βν/4} is real).
    # ------------------------------------------------------------------
    Nν = nu_grid_size
    nu = collect(range(-S, S; length = Nν))
    Δν = T(2) * S / T(Nν - 1)
    f_vec = T[T(freq_kernel(filter, ν)) for ν in nu]
    g_hat = Matrix{CT}(undef, Nν, Nν)
    pref_g = CT(one(T) / (2im))
    @inbounds for q in 1:Nν, p in 1:Nν
        th = tanh(βT * (nu[q] - nu[p]) / 4)
        # f̂ is real → conj is identity in T; keep the structure for clarity.
        g_hat[p, q] = pref_g * Complex{T}(th * f_vec[p] * f_vec[q])
    end

    # ------------------------------------------------------------------
    # Step 2: 2D type-3 NUFFT to the time grid.
    #   g(t_m, t_n) = (Δν)² / (2π)² · Σ_{p, q} ĝ[p, q] · exp(-iν_p t_m + iν'_q t_n)
    #
    # FINUFFT type-3 (isign = +1) computes Σ c[j] · exp(+i (sx[j] tx[k] + sy[j] ty[k])).
    # Feed sx = -ν_p so that exp(+i sx tx) = exp(-i ν_p t_m), and sy = +ν'_q.
    # ------------------------------------------------------------------
    Nsrc = Nν * Nν
    src_x = Vector{Float64}(undef, Nsrc)
    src_y = Vector{Float64}(undef, Nsrc)
    src_c = Vector{ComplexF64}(undef, Nsrc)
    idx = 0
    @inbounds for q in 1:Nν, p in 1:Nν
        idx += 1
        src_x[idx] = -Float64(nu[p])  # negate to flip sign on first dimension
        src_y[idx] = Float64(nu[q])
        src_c[idx] = ComplexF64(g_hat[p, q])
    end

    Ntgt = Nt * Nt
    tgt_x = Vector{Float64}(undef, Ntgt)
    tgt_y = Vector{Float64}(undef, Ntgt)
    @inbounds for nn in 1:Nt, m in 1:Nt
        idx_t = (nn - 1) * Nt + m
        tgt_x[idx_t] = Float64(time_labels[m])
        tgt_y[idx_t] = Float64(time_labels[nn])
    end

    plan = FINUFFT.finufft_makeplan(3, 2, +1, 1, 1e-12; dtype = Float64, nthreads = 1)
    FINUFFT.finufft_setpts!(plan, src_x, src_y, Float64[], tgt_x, tgt_y, Float64[])
    out_g = Vector{ComplexF64}(undef, Ntgt)
    FINUFFT.finufft_exec!(plan, src_c, out_g)
    FINUFFT.finufft_destroy!(plan)

    norm_factor = (Float64(Δν) / (2π))^2
    g_tt = Matrix{CT}(undef, Nt, Nt)
    @inbounds for nn in 1:Nt, m in 1:Nt
        idx_t = (nn - 1) * Nt + m
        g_tt[m, nn] = CT(out_g[idx_t] * norm_factor)
    end

    # ------------------------------------------------------------------
    # Step 3 + 4: filter-agnostic NUFFT contraction g_tt → Q_ijk → G.
    # ------------------------------------------------------------------
    return _dll_coherent_from_g_tt(jumps, hamiltonian, g_tt, time_labels, τ)
end

"""
    _dll_coherent_from_g_tt(jumps, hamiltonian, g_tt, time_labels, τ) -> Matrix

Filter-agnostic helper for `dll_coherent_op_time`. Given the time-domain
kernel matrix `g_tt[m, n] = g(t_m, t_n)` (computed by whatever closed-form
or numerical path the concrete filter supports), evaluate

    G[i, j] = Σ_a Σ_k A^a[i, k] · A^a[k, j]
              · (τ² · Σ_{m, n} g_tt[m, n] · cis((λ_i − λ_k) t_n + (λ_k − λ_j) t_m))

The inner double sum is a 2D type-3 NUFFT (sources on the uniform
`(t_m, t_n)` tensor grid, `n³` targets at the Bohr triples
`(λ_k − λ_j, λ_i − λ_k)`). The result `Q_ijk` is jump-independent and
contracted across jumps via the standard `Σ_a A^a A^a Q` reduction.

This is the "fast path" — the only path used in production. The legacy
`O(Nt² · n³ · |𝓐|)` Riemann sum is preserved as `dll_coherent_op_time_legacy`
for reference, never used here.
"""
function _dll_coherent_from_g_tt(
    jumps::AbstractVector{<:JumpOp},
    hamiltonian::HamHam{T},
    g_tt::AbstractMatrix{<:Complex},
    time_labels::AbstractVector{<:Real},
    τ::Real,
) where {T<:AbstractFloat}
    eigvals = hamiltonian.eigvals
    n = length(eigvals)
    Nt = length(time_labels)
    CT = Complex{T}
    @assert size(g_tt) == (Nt, Nt)

    # ------------------------------------------------------------------
    # Step 3: Q_ijk[i, j, k] = τ² · Σ_{m, n} g_tt[m, n] ·
    #                          cis((λ_i − λ_k) t_n + (λ_k − λ_j) t_m)
    # via 2D type-3 NUFFT (sources on uniform tensor grid, n³ targets).
    # ------------------------------------------------------------------
    src_x = Vector{Float64}(undef, Nt * Nt)  # paired with γ = λ_k − λ_j
    src_y = Vector{Float64}(undef, Nt * Nt)  # paired with α = λ_i − λ_k
    src_c = Vector{ComplexF64}(undef, Nt * Nt)
    idx = 0
    @inbounds for nidx in 1:Nt, m in 1:Nt
        idx += 1
        src_x[idx] = Float64(time_labels[m])
        src_y[idx] = Float64(time_labels[nidx])
        src_c[idx] = ComplexF64(g_tt[m, nidx])
    end

    n_targets = n^3
    tgt_s = Vector{Float64}(undef, n_targets)
    tgt_t = Vector{Float64}(undef, n_targets)
    @inbounds for k in 1:n, j in 1:n, i in 1:n
        idx_t = ((i - 1) + (j - 1) * n + (k - 1) * n^2) + 1
        tgt_s[idx_t] = Float64(eigvals[k] - eigvals[j])
        tgt_t[idx_t] = Float64(eigvals[i] - eigvals[k])
    end

    plan = FINUFFT.finufft_makeplan(3, 2, +1, 1, 1e-12; dtype=Float64, nthreads=1)
    FINUFFT.finufft_setpts!(plan, src_x, src_y, Float64[], tgt_s, tgt_t, Float64[])
    out_q = Vector{ComplexF64}(undef, n_targets)
    FINUFFT.finufft_exec!(plan, src_c, out_q)
    FINUFFT.finufft_destroy!(plan)

    Q_ijk = reshape(out_q, n, n, n) .* (τ^2)

    # ------------------------------------------------------------------
    # Step 4: G[i, j] = Σ_a Σ_k A^a[i, k] · A^a[k, j] · Q_ijk[i, j, k]
    # ------------------------------------------------------------------
    G = zeros(CT, n, n)
    @inbounds for jump in jumps
        A_eb = jump.in_eigenbasis
        for j in 1:n, i in 1:n
            acc = CT(0)
            for k in 1:n
                acc += A_eb[i, k] * A_eb[k, j] * Q_ijk[i, j, k]
            end
            G[i, j] += acc
        end
    end
    return G
end

"""
    dll_coherent_op_time_legacy(jumps, hamiltonian, time_labels, filter, beta, τ; nu_grid=nothing)

Reference implementation of `dll_coherent_op_time` using a 2D internal
`(ν, ν')`-tabulation grid for `ĝ^a(ν, ν')` (Eq. 3.5) followed by a 2D
inverse FT to the time grid and an explicit `O(Nt² · n³)` Riemann sum.
Superseded by the closed-form + NUFFT path in `dll_coherent_op_time`
(Phase C / qf-hur.3); kept here as a numerical reference for tests.

The default `(ν, ν')` grid is centered at `-1/β` with half-width `12/β`
(`Nν ≥ 64`) and resolves `ĝ^a` to ≤ 1e-12 in the tail. Pass `nu_grid` to
override.
"""
function dll_coherent_op_time_legacy(
    jumps::AbstractVector{<:JumpOp},
    hamiltonian::HamHam{T},
    time_labels::AbstractVector{<:Real},
    filter::AbstractFilter,
    beta::Real,
    τ::Real;
    nu_grid::Union{Nothing, AbstractVector{<:Real}} = nothing,
) where {T<:AbstractFloat}
    eigvals = hamiltonian.eigvals
    n = length(eigvals)
    CT = Complex{T}

    νs = if nu_grid === nothing
        ν_centre = -one(T) / T(beta)
        ν_half = T(12) / T(beta)
        Δν_target = T(beta) / 16
        Nν_pre = 2 * ceil(Int, ν_half / Δν_target)
        Nν = max(Nν_pre, 64)
        collect(range(ν_centre - ν_half, ν_centre + ν_half; length=Nν))
    else
        collect(nu_grid)
    end
    Nν = length(νs)
    Δν = νs[2] - νs[1]

    G_hat = Matrix{CT}(undef, Nν, Nν)
    @inbounds for q in 1:Nν, p in 1:Nν
        G_hat[p, q] = dll_coherent_kernel_bohr(filter, νs[p], νs[q])
    end

    Nt = length(time_labels)
    pref_g = CT(Δν * Δν / (2 * T(π))^2)

    Φ = Matrix{CT}(undef, Nt, Nν)
    Ψ = Matrix{CT}(undef, Nt, Nν)
    @inbounds for p in 1:Nν, m in 1:Nt
        Φ[m, p] = cis(-νs[p] * time_labels[m])
    end
    @inbounds for q in 1:Nν, n in 1:Nt
        Ψ[n, q] = cis(νs[q] * time_labels[n])
    end

    tmp_FH = Matrix{CT}(undef, Nt, Nν)
    g_tt = Matrix{CT}(undef, Nt, Nt)
    mul!(tmp_FH, Φ, G_hat)
    mul!(g_tt, tmp_FH, transpose(Ψ))
    rmul!(g_tt, pref_g)

    G = zeros(CT, n, n)
    weight_outer = τ^2

    phases_t = Matrix{CT}(undef, Nt, n)
    @inbounds for k in 1:n, m in 1:Nt
        phases_t[m, k] = cis(eigvals[k] * time_labels[m])
    end

    Atm = Matrix{CT}(undef, n, n)
    Atn = Matrix{CT}(undef, n, n)
    prod_buf = Matrix{CT}(undef, n, n)

    @inbounds for jump in jumps
        A_eb = jump.in_eigenbasis
        for nidx in 1:Nt
            for j in 1:n, i in 1:n
                Atn[i, j] = phases_t[nidx, i] * conj(phases_t[nidx, j]) * A_eb[i, j]
            end
            for m in 1:Nt
                for j in 1:n, i in 1:n
                    Atm[i, j] = phases_t[m, i] * conj(phases_t[m, j]) * A_eb[i, j]
                end
                mul!(prod_buf, Atn, Atm)
                w = g_tt[m, nidx] * weight_outer
                for j in 1:n, i in 1:n
                    G[i, j] += w * prod_buf[i, j]
                end
            end
        end
    end
    return G
end

# ---------------------------------------------------------------------------
# DLL Kossakowski matrix α^{DLL}_{ν, ν'} (Ding–Li–Lin 2024, Sec. 4)
# ---------------------------------------------------------------------------
#
# Per coupling A^a, the DLL Kossakowski matrix on the Bohr-frequency grid is
# (DLL paper Sec. 4 / Remark 23):
#
#   α^{DLL}_{ν, ν'} = e^{-β(ν+ν')/4} q^a(ν) conj(q^a(ν'))
#                   = freq_kernel(filter, ν) · conj(freq_kernel(filter, ν'))
#
# This is RANK-1 by construction (outer product `v · v†` with
# v_k = freq_kernel(filter, ν_k)). In contrast the CKG Kossakowski (Eq. 4.6,
# `α^{CKG}_{ν, ν'} = (2π)² ∫ γ(ω) f̂(ω-ν) conj(f̂(ω-ν')) dω`) is generally
# full-rank (Fig. 3 of the paper). The rank gap is the structural distinction
# between DLL (one Lindblad operator per coupling) and CKG (a `K × K`
# Kossakowski per coupling, where `K = |B_H|`).

"""
    dll_kossakowski_bohr(filter, bohr_freqs::AbstractVector{<:Real}) -> Matrix

DLL Kossakowski matrix α^{DLL}_{ν, ν'} per coupling on a given Bohr-frequency
grid (Ding–Li–Lin 2024, Sec. 4):

    α^{DLL}_{ν, ν'} = freq_kernel(filter, ν) · conj(freq_kernel(filter, ν'))

This is a rank-1 outer product `v · v†` with `v_k = freq_kernel(filter, ν_k)`.
The result is the same for every coupling A^a because the DLL construction
factors the dissipator through a single weighting function `q^a` that is
coupling-independent in the `DLLGaussianFilter` case.

# Arguments
- `filter`: any `AbstractFilter`. For non-DLL filters, the formula above
  reduces to whatever `freq_kernel` returns (which may not be the proper
  KMS-DB Kossakowski — caller's responsibility).
- `bohr_freqs`: vector of (typically unique) Bohr frequencies `ν ∈ B_H`.

# Returns
A `length(bohr_freqs) × length(bohr_freqs)` complex matrix.
"""
function dll_kossakowski_bohr(
    filter::AbstractFilter,
    bohr_freqs::AbstractVector{<:Real},
)
    K = length(bohr_freqs)
    v = [freq_kernel(filter, ν) for ν in bohr_freqs]
    α = Matrix{eltype(v)}(undef, K, K)
    @inbounds for q in 1:K, p in 1:K
        α[p, q] = v[p] * conj(v[q])
    end
    return α
end

"""
    dll_kossakowski_bohr(filter, hamiltonian::HamHam) -> (alpha, bohr_freqs)

Convenience overload that uses the unique Bohr frequencies from
`hamiltonian.bohr_dict`. Returns the Kossakowski matrix and the matching
frequency vector (so the caller can label rows/columns).
"""
function dll_kossakowski_bohr(
    filter::AbstractFilter,
    hamiltonian::HamHam{T},
) where {T<:AbstractFloat}
    bohr_freqs = sort!(collect(keys(hamiltonian.bohr_dict)))
    α = dll_kossakowski_bohr(filter, bohr_freqs)
    return α, bohr_freqs
end
