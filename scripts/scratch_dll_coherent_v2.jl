#=
Phase C prototype: closed-form g^a(t, t') for `DLLGaussianFilter` + NUFFT-factored
Riemann sum for `dll_coherent_op_time`.

DLL paper (Ding-Li-Lin 2024) Eq. 3.5 + 3.6:
    ĝ(ν, ν') = (1/2i) tanh(β(ν'-ν)/4) f̂(ν) f̂(ν')      (f̂ real for DLL)
    g(t, t') = (1/(2π)²) ∫∫ ĝ(ν, ν') e^{-iνt + iν't'} dν dν'

with f̂(ν) = e^{-β²ν²/8 - βν/4} for `DLLGaussianFilter`.

Substitute (u, v) = (ν' - ν, (ν + ν')/2), so ν = v - u/2, ν' = v + u/2:
    f̂(v - u/2) f̂(v + u/2) = exp(-β²v²/4 - β²u²/16 - βv/2)

The (ν, ν') exponent rewrites as:
    -iνt + iν't' = iv(t' - t) + iu(t + t')/2

So g(t, t') = (1/(2π)² · 2i) · I_v(t' - t) · I_u(t + t')

where
    I_v(δ) := ∫ exp(-β²v²/4 - βv/2 + ivδ) dv          (closed-form Gaussian)
    I_u(s) := ∫ tanh(βu/4) · exp(-β²u²/16 + ius/2) du   (1D quadrature)

Closed form for I_v:
    I_v(δ) = (2√π/β) · exp(1/4 - iδ/β - δ²/β²)

I_u(s) is purely imaginary (tanh is odd), so write I_u(s) = i · J(s) with
    J(s) = -2 ∫_0^∞ tanh(βu/4) · exp(-β²u²/16) · sin(us/2) du     (real)

Final structure:
    g(t, t') = factor_diff(t - t') · J(t + t') · constant

where factor_diff is closed-form (no quadrature) and J is a 1D function that
only needs 2Nt - 1 distinct s = t_m + t_n samples (instead of Nt² 2D ones).

This script:
  (a) Verifies the closed-form decomposition against the current
      internal-(ν,ν')-grid implementation in `dll_coherent_op_time`.
  (b) Verifies that a 2D type-3 NUFFT factors the Riemann sum exactly.

Run: julia --project scripts/scratch_dll_coherent_v2.jl
=#

using QuantumFurnace
using LinearAlgebra
using Printf
using QuadGK
using FINUFFT

# ----- Closed-form ingredients ------------------------------------------

# I_v(δ) = (2√π/β) · exp(1/4 - iδ/β - δ²/β²)
@inline function _factor_diff(beta::Real, δ::Real)
    return (2 * sqrt(π) / beta) * exp(1//4 - im * δ / beta - δ^2 / beta^2)
end

# I_u(s) = ∫_{-∞}^∞ tanh(β u / 4) · exp(-β² u² / 16) · e^{i u s / 2} du
# Tanh is odd, Gaussian is even, so the cosine part vanishes and
#     I_u(s) = i · J(s),     J(s) = 2 ∫_0^∞ tanh(β u / 4) · exp(-β² u² / 16) · sin(u s / 2) du.
# Use QuadGK on [0, 12/β] — Gaussian width ~4/β, tail beyond 3σ is below 1e-15.
function _J_quadrature(beta::Real, s::Real; rtol::Real=1e-12, atol::Real=1e-14)
    integrand(u) = tanh(beta * u / 4) * exp(-(beta * u)^2 / 16) * sin(u * s / 2)
    # Gaussian width 4/β; cut at 24/β where exp(-(βu/4)²) ≈ exp(-36) ≈ 2e-16.
    cutoff = 24 / beta
    val, _ = quadgk(integrand, 0.0, cutoff; rtol=rtol, atol=atol)
    return 2 * val
end

# Combined: g(t, t') = (1/(2π)²) · (1/2i) · I_v(t-t') · I_u(t+t')
# I_u(s) = i · J(s); so 1/(2i) · I_u(s) = J(s) / 2.
# Hence g(t, t') = (1/(2π)²) · (1/2) · I_v(t' - t) · J(t + t')   (note I_v argument sign)
# Wait let me re-check the sign convention.

# Going back to derivation: integrand has e^{iv(t'-t)}, so I_v(δ) is evaluated at δ = t' - t.
@inline function _g_closed_form(beta::Real, t::Real, tp::Real, J_val::Number)
    return (1 / (2π)^2) * (1//2) * _factor_diff(beta, tp - t) * J_val
end

# Tabulated J(s) on a 1D grid of unique s = t_m + t_n values.
function _build_J_table(beta::Real, time_labels::AbstractVector{<:Real}; rtol::Real=1e-12)
    Nt = length(time_labels)
    # All unique values of t_m + t_n.
    unique_sums = sort!(unique([time_labels[m] + time_labels[n] for m in 1:Nt for n in 1:Nt]))
    J_vals = Vector{Float64}(undef, length(unique_sums))
    for (i, s) in enumerate(unique_sums)
        J_vals[i] = _J_quadrature(beta, s; rtol=rtol)
    end
    sum_to_index = Dict{Float64, Int}(s => i for (i, s) in enumerate(unique_sums))
    return unique_sums, J_vals, sum_to_index
end

# ----- Setup -------------------------------------------------------------
# Toggle PROBLEM_SIZE to switch between small (n=3, fast diagnose) and full (n=5, benchmark).
const PROBLEM_SIZE = get(ENV, "DLL_PROBLEM_SIZE", "small")
n = PROBLEM_SIZE == "full" ? 5 : 3
beta = 5.0
println("Loading n=$n Heisenberg fixture, β=$beta (PROBLEM_SIZE=$PROBLEM_SIZE)...")
ham = load_hamiltonian("heis", n; beta=beta)
dim = size(ham.data, 1)

# Single jump for easy timing
A1 = Matrix(pad_term([X], n, 1)) ./ sqrt(3*n)
B1 = ham.eigvecs' * A1 * ham.eigvecs
jumps = JumpOp[JumpOp(A1, B1, true, true)]

filter = DLLGaussianFilter(beta)

# Time grid: at n=5, the benchmark uses N_BITS=12 → Nt ≈ 1173.
N_BITS = PROBLEM_SIZE == "full" ? 12 : 10
W0 = 0.05
T0 = 2π / (2^N_BITS * W0)
N = 2^N_BITS
raw_time_labels = collect((-N÷2):(N÷2 - 1)) .* T0
time_labels = QuantumFurnace._truncate_time_labels_for_oft(raw_time_labels, 2/beta; filter=filter)
Nt = length(time_labels)
@printf "Truncated grid: Nt = %d (raw N = %d)\n" Nt N

# ----- (1) Closed-form g_tt vs internal-grid g_tt ------------------------
println("\n--- Step 1: Closed-form g(t, t') vs internal-grid -----")

# Closed-form g_tt
println("Building J(s) table (1D quadrature)...")
t_J = @elapsed (unique_sums, J_vals, sum_to_index) = _build_J_table(beta, time_labels)
@printf "  J table: %d unique sums, %.4f s\n" length(unique_sums) t_J

println("Building g_tt closed form...")
t_gtt_cf = @elapsed begin
    g_tt_cf = Matrix{ComplexF64}(undef, Nt, Nt)
    for nidx in 1:Nt, m in 1:Nt
        s = time_labels[m] + time_labels[nidx]
        J_val = J_vals[sum_to_index[s]]
        g_tt_cf[m, nidx] = _g_closed_form(beta, time_labels[m], time_labels[nidx], J_val)
    end
end
@printf "  g_tt closed form: %.4f s\n" t_gtt_cf

# Reference: rebuild g_tt the same way the existing dll_coherent_op_time does.
println("Building g_tt via internal-(ν,ν')-grid (Nν=64 default)...")
νs = let
    ν_centre = -1.0 / beta
    ν_half = 12.0 / beta
    Δν_target = beta / 16
    Nν_pre = 2 * ceil(Int, ν_half / Δν_target)
    Nν = max(Nν_pre, 64)
    collect(range(ν_centre - ν_half, ν_centre + ν_half; length=Nν))
end
Nν = length(νs)
Δν = νs[2] - νs[1]
G_hat = Matrix{ComplexF64}(undef, Nν, Nν)
for q in 1:Nν, p in 1:Nν
    G_hat[p, q] = QuantumFurnace.dll_coherent_kernel_bohr(filter, νs[p], νs[q])
end
pref_g = ComplexF64(Δν * Δν / (2π)^2)
Φ = Matrix{ComplexF64}(undef, Nt, Nν)
Ψ = Matrix{ComplexF64}(undef, Nt, Nν)
for p in 1:Nν, m in 1:Nt
    Φ[m, p] = cis(-νs[p] * time_labels[m])
end
for q in 1:Nν, nidx in 1:Nt
    Ψ[nidx, q] = cis(νs[q] * time_labels[nidx])
end
g_tt_grid = Φ * G_hat * transpose(Ψ) .* pref_g

# Compare
err = opnorm(g_tt_cf - g_tt_grid)
rel = err / opnorm(g_tt_grid)
@printf "  ‖g_cf − g_grid‖_op = %.3e   relative = %.3e\n" err rel

# Increase ν-grid to confirm the closed-form is actually the more accurate one
νs_fine = let
    Nν_fine = 1024
    ν_centre = -1.0 / beta
    ν_half = 18.0 / beta  # wider tail
    collect(range(ν_centre - ν_half, ν_centre + ν_half; length=Nν_fine))
end
Nν_fine = length(νs_fine)
Δν_fine = νs_fine[2] - νs_fine[1]
G_hat_fine = Matrix{ComplexF64}(undef, Nν_fine, Nν_fine)
for q in 1:Nν_fine, p in 1:Nν_fine
    G_hat_fine[p, q] = QuantumFurnace.dll_coherent_kernel_bohr(filter, νs_fine[p], νs_fine[q])
end
pref_g_fine = ComplexF64(Δν_fine * Δν_fine / (2π)^2)
Φ_fine = Matrix{ComplexF64}(undef, Nt, Nν_fine)
Ψ_fine = Matrix{ComplexF64}(undef, Nt, Nν_fine)
for p in 1:Nν_fine, m in 1:Nt
    Φ_fine[m, p] = cis(-νs_fine[p] * time_labels[m])
end
for q in 1:Nν_fine, nidx in 1:Nt
    Ψ_fine[nidx, q] = cis(νs_fine[q] * time_labels[nidx])
end
g_tt_fine = Φ_fine * G_hat_fine * transpose(Ψ_fine) .* pref_g_fine

err_cf_fine = opnorm(g_tt_cf - g_tt_fine)
rel_cf_fine = err_cf_fine / opnorm(g_tt_fine)
@printf "  ‖g_cf − g_finegrid‖_op = %.3e   relative = %.3e   (Nν=1024)\n" err_cf_fine rel_cf_fine

# ----- (2) NUFFT-factored Riemann sum vs explicit Riemann sum -----------
println("\n--- Step 2: NUFFT-factored Riemann sum vs explicit -----")

# Build G via the existing (g_tt-based) Riemann sum directly. We use g_tt_cf
# as the source so the only difference between paths is the DFT method
# (explicit Σ_{m,n} vs 2D type-3 NUFFT).
function _G_explicit(g_tt::AbstractMatrix, jumps, ham, time_labels::AbstractVector{<:Real}, t0::Real)
    n = size(ham.data, 1)
    eigvals = ham.eigvals
    Nt = length(time_labels)
    weight_outer = t0^2

    phases_t = Matrix{ComplexF64}(undef, Nt, n)
    @inbounds for k in 1:n, m in 1:Nt
        phases_t[m, k] = cis(eigvals[k] * time_labels[m])
    end

    G = zeros(ComplexF64, n, n)
    Atm = Matrix{ComplexF64}(undef, n, n)
    Atn = Matrix{ComplexF64}(undef, n, n)
    prod_buf = Matrix{ComplexF64}(undef, n, n)
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
                # Operator order: A^a(t_n) · A^a(t_m)  (paper has typo in 3.7 third)
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

println("Building G via explicit Riemann sum (g_tt_cf as source)...")
t_explicit = @elapsed G_explicit = _G_explicit(g_tt_cf, jumps, ham, time_labels, T0)
@printf "  G explicit: %.3f s\n" t_explicit

# NUFFT-factored Riemann sum:
#   Q(α, γ) = τ² · Σ_{m,n} g_tt[m, n] · cis(α · t_n + γ · t_m)
#   G[i, j] = Σ_a Σ_k A^a[i, k] · A^a[k, j] · Q(λ_i − λ_k, λ_k − λ_j)
function _G_nufft(g_tt::AbstractMatrix, jumps, ham, time_labels::AbstractVector{<:Real}, t0::Real; eps::Real=1e-12)
    n = size(ham.data, 1)
    eigvals = ham.eigvals
    Nt = length(time_labels)

    # Build flat source coords (xj, yj) for FINUFFT (uniform tensor grid).
    src_x = Vector{Float64}(undef, Nt * Nt)  # xj = t_m  (paired with γ in target)
    src_y = Vector{Float64}(undef, Nt * Nt)  # yj = t_n  (paired with α in target)
    src_c = Vector{ComplexF64}(undef, Nt * Nt)  # source values g_tt[m, n]
    idx = 0
    @inbounds for nidx in 1:Nt, m in 1:Nt
        idx += 1
        src_x[idx] = time_labels[m]
        src_y[idx] = time_labels[nidx]
        src_c[idx] = g_tt[m, nidx]
    end

    # Targets: (γ, α) = (λ_k − λ_j, λ_i − λ_k) for all (i, j, k) ∈ [1, n]³.
    # Use n³ targets (with potential duplicates) for simplicity.
    n_targets = n^3
    tgt_s = Vector{Float64}(undef, n_targets)  # s = γ
    tgt_t = Vector{Float64}(undef, n_targets)  # t = α
    @inbounds for k in 1:n, j in 1:n, i in 1:n
        idx_t = ((i - 1) + (j - 1) * n + (k - 1) * n^2) + 1
        tgt_s[idx_t] = eigvals[k] - eigvals[j]
        tgt_t[idx_t] = eigvals[i] - eigvals[k]
    end

    # 2D type-3 NUFFT, sign +1: out_t = Σ_j c_j · cis(+1 · (s_t · x_j + t_t · y_j))
    plan = FINUFFT.finufft_makeplan(3, 2, +1, 1, eps; dtype=Float64, nthreads=1)
    FINUFFT.finufft_setpts!(plan, src_x, src_y, Float64[], tgt_s, tgt_t, Float64[])
    out_q = Vector{ComplexF64}(undef, n_targets)
    FINUFFT.finufft_exec!(plan, src_c, out_q)
    FINUFFT.finufft_destroy!(plan)

    # Apply τ² weight and reshape into (i, j, k) array.
    Q_ijk = reshape(out_q, n, n, n) .* (t0^2)  # Q_ijk[i, j, k] = Q(λ_k − λ_j, λ_i − λ_k)

    # Final assembly: G[i, j] = Σ_a Σ_k A^a[i, k] · A^a[k, j] · Q_ijk[i, j, k]
    G = zeros(ComplexF64, n, n)
    @inbounds for jump in jumps
        A_eb = jump.in_eigenbasis
        for j in 1:n, i in 1:n
            acc = ComplexF64(0)
            for k in 1:n
                acc += A_eb[i, k] * A_eb[k, j] * Q_ijk[i, j, k]
            end
            G[i, j] += acc
        end
    end
    return G
end

println("Building G via NUFFT-factored Riemann sum (warmup + timed)...")
_ = _G_nufft(g_tt_cf, jumps, ham, time_labels, T0)
t_nufft = @elapsed G_nufft = _G_nufft(g_tt_cf, jumps, ham, time_labels, T0)
@printf "  G NUFFT:    %.4f s\n" t_nufft

err_g = opnorm(G_explicit - G_nufft)
rel_g = err_g / opnorm(G_explicit)
@printf "  ‖G_explicit − G_nufft‖_op = %.3e   relative = %.3e\n" err_g rel_g

# Also compare against the public dll_coherent_op_time.
println("\n--- Step 3: Cross-check vs current dll_coherent_op_time ----")
G_public = QuantumFurnace.dll_coherent_op_time(jumps, ham, time_labels, filter, beta, T0)
err_public_explicit = opnorm(G_explicit - G_public)
err_public_nufft = opnorm(G_nufft - G_public)
@printf "  ‖G_explicit − G_public‖_op (closed-form vs Nν=64): %.3e\n" err_public_explicit
@printf "  ‖G_nufft − G_public‖_op:    %.3e\n" err_public_nufft

@printf "\nSpeedup (G_explicit time / G_nufft time): %.1fx\n" t_explicit / t_nufft

println("\nDONE")
