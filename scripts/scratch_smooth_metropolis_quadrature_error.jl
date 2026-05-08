#!/usr/bin/env julia
#
# Smooth-Metropolis ω-quadrature error analysis (qf-3il).
#
# Headline question
# -----------------
# How many ω-grid points (= qubits in the energy/time register) are
# needed to evaluate the CKG ω-integral
#
#     α_{ν₁,ν₂} = ∫ γ_M^{(s)}(ω) f̂(ω − ν₁) f̂(ω − ν₂) dω
#
# to a target precision ε, as a function of s ∈ {0, 0.1, 0.25, 0.5, 1.0}
# and β ∈ {5, 10, 20}? The kinky case (s = 0) has a discontinuous γ' at
# ω = −β σ²/2 → only polynomial trapezoidal convergence; smooth s > 0 is
# Gevrey-1/2 (Prop. prop:smooth-metro-gevrey, thesis §2_methods.tex) and
# converges super-polynomially.
#
# This analysis depends on β, σ = 1/β, s only — *not* on n. So the qubit
# argument extends uniformly to large n where full mixing-time sweeps are
# out of reach.
#
# Strategy
# --------
# 1. Define γ_M^{(s)}(ω) via QuantumFurnace.pick_transition (same code path
#    the simulator uses).
# 2. For (β, s), compute the reference integral I_ref via QuadGK with a
#    tight tolerance and a wide window.
# 3. For N_ω ∈ {8, 16, …, 8192}, compute I_trap(N_ω) on [−W, W] using the
#    midpoint Riemann sum (the convention `_create_energy_labels` uses;
#    trapezoidal differs only by O(1/N) at endpoints which contribute
#    Gaussian-tail-negligible mass for our window choice).
# 4. Plot rel_err vs N_ω on log-log axes. Tag the smallest N_ω achieving
#    rel_err < ε for ε ∈ {1e-3, 1e-4, 1e-6}; convert to bit count
#    log₂(N_ω).
# 5. Repeat for representative ν-pairs: (ν₁, ν₂) = (0, 0) (centred), and
#    (0.5, 0.5), (1.0, 1.0) — to test sensitivity to where the integrand
#    mass sits relative to the kink (the kink is at ω = −β σ²/2 = −1/(2β),
#    deep inside the support for ν ≈ 0).
#
# Window choice
# -------------
# W = 10 / β. The Gaussian filters f̂(ω − ν) have stddev σ = 1/β, so
# 10 σ tail bound is ~exp(−50) ≪ any ε we care about, ensuring the
# truncation error is dominated by the discretisation error.
#
# PHYSICS CHECK: the Riemann sum with midpoint nodes matches the energy
# register convention used in `src/energy_domain.jl` (centred, equally
# spaced, Δω = 2W/N).
#
# Output
# ------
# - drafts/figures/numerics/smooth_metro_quadrature_error.bson  (raw data)
# - drafts/figures/numerics/smooth_metro_quadrature_error.{png,pdf}
# - drafts/figures/numerics/smooth_metro_qubits_vs_s.{png,pdf}
# - Summary table to stdout: N_ω → bits at each (β, s, ε).

using Printf
using LinearAlgebra
using BSON
using QuadGK
using QuantumFurnace
using QuantumFurnace: pick_transition
ENV["GKSwstype"] = "100"
using Plots

# ── Sweep grid ────────────────────────────────────────────────────────────────
const BETA_VALUES = [5.0, 10.0, 20.0]
const S_VALUES    = [0.0, 0.1, 0.25, 0.5, 1.0]
const NU_PAIRS    = [(0.0, 0.0), (0.5, 0.5), (1.0, 1.0)]      # representative
const N_GRID      = 2 .^ (3:13)                                 # 8 .. 8192
const EPS_TARGETS = [1e-3, 1e-4, 1e-6]

# ── Output paths ──────────────────────────────────────────────────────────────
const OUT_DIR = joinpath(@__DIR__, "..", "drafts", "figures", "numerics")
mkpath(OUT_DIR)
const BSON_OUT  = joinpath(OUT_DIR, "smooth_metro_quadrature_error.bson")
const FIG_ERR_PNG = joinpath(OUT_DIR, "smooth_metro_quadrature_error.png")
const FIG_ERR_PDF = joinpath(OUT_DIR, "smooth_metro_quadrature_error.pdf")
const FIG_QUB_PNG = joinpath(OUT_DIR, "smooth_metro_qubits_vs_s.png")
const FIG_QUB_PDF = joinpath(OUT_DIR, "smooth_metro_qubits_vs_s.pdf")

# ── Thesis colour palette ─────────────────────────────────────────────────────
const COLOR_S = Dict(
    0.0  => "#222222",   # black (kinky baseline)
    0.1  => "#5C7794",   # slateblue
    0.25 => "#5F8B8E",   # dustyteal
    0.5  => "#B5654A",   # terracotta
    1.0  => "#7A2E39",   # bordeaux
)
const STYLE_S = Dict(
    0.0  => :dash,
    0.1  => :solid,
    0.25 => :solid,
    0.5  => :solid,
    1.0  => :solid,
)

# ── Building blocks ───────────────────────────────────────────────────────────
"""Build the γ_M^{(s)}(ω) closure via the simulator's pick_transition."""
function gamma_smooth(beta::Real, s::Real)
    cfg = Config(;
        sim                       = Lindbladian(),
        domain                    = EnergyDomain(),
        construction              = KMS(),
        num_qubits                = 1,        # immaterial for γ
        with_linear_combination   = true,
        beta                      = float(beta),
        sigma                     = 1.0 / float(beta),
        a                         = 0.0,
        s                         = float(s),
        num_energy_bits           = 12,
        w0                        = 0.05,
        t0                        = 2π / (2^12 * 0.05),
        num_trotter_steps_per_t0  = 10,
    )
    return pick_transition(cfg)
end

"""Gaussian filter f̂(ω) = exp(−ω²/(4σ²)) / (σ √(2π))^(1/2). Mirrors the
energy_domain integrand normalisation so the absolute α value matches the
simulator; only the *relative* error against the reference matters here."""
function f_hat(omega::Real, sigma::Real)
    return exp(- omega^2 / (4 * sigma^2)) * sqrt(1.0 / (sigma * sqrt(2π)))
end

"""Window that always contains the integrand mass.
The Gaussian filters f̂(ω−ν) have stddev σ = 1/β; γ has decay scale 1/β.
We pick the larger of (10/β) and |ν|+(10/β) on each side so the truncation
error is ≪ any ε we care about even at non-zero ν."""
window_W(beta::Real, nu1::Real, nu2::Real) = max(abs(nu1), abs(nu2)) + 10.0 / beta

"""Reference integral via QuadGK, very tight tolerance, wide window."""
function alpha_ref(beta::Real, s::Real, nu1::Real, nu2::Real;
                   W::Real = window_W(beta, nu1, nu2) + 20.0 / beta,
                   rtol::Real = 1e-12)
    sigma = 1.0 / beta
    g     = gamma_smooth(beta, s)
    integrand(omega) = g(omega) * f_hat(omega - nu1, sigma) * f_hat(omega - nu2, sigma)
    val, _ = quadgk(integrand, -W, W; rtol = rtol, atol = 1e-16)
    return val
end

"""Midpoint Riemann sum on [−W, W] with N nodes (matches the centred energy
register convention used in `_create_energy_labels`)."""
function alpha_riemann(beta::Real, s::Real, nu1::Real, nu2::Real,
                        N::Integer; W::Real = window_W(beta, nu1, nu2))
    sigma = 1.0 / beta
    g     = gamma_smooth(beta, s)
    dw    = 2W / N
    # Centred nodes ω_k = -W + (k + 1/2) dw, k = 0 .. N-1
    s_acc = 0.0
    @inbounds for k in 0:(N - 1)
        omega = -W + (k + 0.5) * dw
        s_acc += g(omega) * f_hat(omega - nu1, sigma) * f_hat(omega - nu2, sigma)
    end
    return s_acc * dw
end

# ── Driver ────────────────────────────────────────────────────────────────────
println("="^72)
println("Smooth-Metropolis quadrature error (qf-3il)")
println("="^72)
@printf("β values  : %s\n", string(BETA_VALUES))
@printf("s values  : %s\n", string(S_VALUES))
@printf("ν pairs   : %s\n", string(NU_PAIRS))
@printf("N_ω grid  : %s\n", string(N_GRID))
@printf("ε targets : %s\n", string(EPS_TARGETS))
println("="^72)
flush(stdout)

# Triple loop: (β, s, ν-pair) × N_ω. Ref integral once per (β, s, ν-pair).
# Result struct: relative error array, plus per-(β,s,ν-pair,ε) min N_ω.
Result = Dict{Tuple{Float64, Float64, Tuple{Float64, Float64}}, NamedTuple}()

for beta in BETA_VALUES
    for s in S_VALUES
        for (nu1, nu2) in NU_PAIRS
            ref = alpha_ref(beta, s, nu1, nu2)
            rel_errs = Float64[]
            for N in N_GRID
                approx = alpha_riemann(beta, s, nu1, nu2, N)
                push!(rel_errs, abs(approx - ref) / max(abs(ref), 1e-30))
            end
            min_N_for_eps = Dict{Float64, Int}()
            for eps in EPS_TARGETS
                idx = findfirst(<=(eps), rel_errs)
                min_N_for_eps[eps] = idx === nothing ? 0 : N_GRID[idx]
            end
            Result[(beta, s, (nu1, nu2))] = (
                ref          = ref,
                rel_errs     = rel_errs,
                min_N_for_eps = min_N_for_eps,
            )
            @printf("β=%-5.1f s=%-5.2f ν=%-12s   I_ref=%-10.3e   min N (1e-3, 1e-4, 1e-6) = (%5d, %5d, %5d)\n",
                    beta, s, string((nu1, nu2)), ref,
                    min_N_for_eps[1e-3], min_N_for_eps[1e-4], min_N_for_eps[1e-6])
            flush(stdout)
        end
    end
end

# ── Persist ───────────────────────────────────────────────────────────────────
bson_payload = Dict{Symbol, Any}(
    :result       => Result,
    :beta_values  => BETA_VALUES,
    :s_values     => S_VALUES,
    :nu_pairs     => NU_PAIRS,
    :N_grid       => collect(N_GRID),
    :eps_targets  => EPS_TARGETS,
    :git_sha      => try
        read(`git -C $(dirname(@__DIR__)) rev-parse HEAD`, String) |> strip
    catch
        "unknown"
    end,
)
BSON.bson(BSON_OUT, bson_payload)
@printf("\nWrote BSON: %s\n", BSON_OUT)

# ── Plot 1: rel_err vs N_ω, panel per β, ν=(0,0) ──────────────────────────────
plots_err = Plots.Plot[]
for beta in BETA_VALUES
    p = plot(;
        xlabel = "N_ω (grid points)",
        ylabel = "relative quadrature error",
        title  = "β = $(beta) — ν=(0,0)",
        xscale = :log10,
        yscale = :log10,
        legend = :bottomleft,
        xticks = ([8, 32, 128, 512, 2048, 8192], ["8","32","128","512","2048","8192"]),
        framestyle = :box,
        size = (520, 380),
    )
    for s in S_VALUES
        rec = Result[(beta, s, (0.0, 0.0))]
        # Filter out exact zeros (replace with floor for log plot)
        ys = [max(e, 1e-16) for e in rec.rel_errs]
        plot!(p, collect(N_GRID), ys;
            color = COLOR_S[s],
            linestyle = STYLE_S[s],
            linewidth = 2,
            marker = :circle, markersize = 3, markerstrokewidth = 0,
            label = "s = $(s)$(s == 0.0 ? " (kinky)" : "")",
        )
    end
    # Reference horizontal lines at ε targets
    for eps in EPS_TARGETS
        hline!(p, [eps]; color = :black, linestyle = :dot, linewidth = 0.7,
               alpha = 0.5, label = "")
    end
    push!(plots_err, p)
end
fig_err = plot(plots_err...; layout = (1, length(BETA_VALUES)),
                size = (520 * length(BETA_VALUES), 380))
savefig(fig_err, FIG_ERR_PNG)
savefig(fig_err, FIG_ERR_PDF)
@printf("Wrote: %s, %s\n", FIG_ERR_PNG, FIG_ERR_PDF)

# ── Plot 2: required qubits vs s, panel per β, ν=(0,0), curves per ε ─────────
plots_qub = Plots.Plot[]
COLOR_EPS = Dict(1e-3 => "#5C7794", 1e-4 => "#B5654A", 1e-6 => "#7A2E39")
for beta in BETA_VALUES
    p = plot(;
        xlabel = "s",
        ylabel = "required bits log₂(N_ω)",
        title  = "β = $(beta) — ν=(0,0)",
        legend = :topright,
        framestyle = :box,
        size = (520, 380),
    )
    for eps in EPS_TARGETS
        bits = Float64[]
        ss   = Float64[]
        for s in S_VALUES
            N = Result[(beta, s, (0.0, 0.0))].min_N_for_eps[eps]
            if N > 0
                push!(bits, log2(N))
                push!(ss, s)
            else
                push!(bits, NaN)   # never reached ε within N_GRID range
                push!(ss, s)
            end
        end
        plot!(p, ss, bits;
            color = COLOR_EPS[eps],
            linewidth = 2,
            marker = :diamond, markersize = 5, markerstrokewidth = 0,
            label = "ε = $(eps)",
        )
    end
    push!(plots_qub, p)
end
fig_qub = plot(plots_qub...; layout = (1, length(BETA_VALUES)),
                size = (520 * length(BETA_VALUES), 380))
savefig(fig_qub, FIG_QUB_PNG)
savefig(fig_qub, FIG_QUB_PDF)
@printf("Wrote: %s, %s\n", FIG_QUB_PNG, FIG_QUB_PDF)

# ── Summary table ─────────────────────────────────────────────────────────────
println("\n" * "="^88)
println("SUMMARY — required bits log₂(N_ω) at ν=(0,0), ε = 1e-4")
println("="^88)
@printf("%-8s", "β\\s")
for s in S_VALUES
    @printf("%-12s", "s=$(s)")
end
println()
for beta in BETA_VALUES
    @printf("%-8s", "β=$(beta)")
    for s in S_VALUES
        N = Result[(beta, s, (0.0, 0.0))].min_N_for_eps[1e-4]
        if N > 0
            @printf("%-12s", "$(Int(log2(N)))  ($(N))")
        else
            @printf("%-12s", ">$(Int(log2(maximum(N_GRID))))")
        end
    end
    println()
end
println()
println("="^88)
println("SUMMARY — bit savings smooth (s=0.25) vs kinky (s=0) at ε = 1e-4")
println("="^88)
@printf("%-8s%-12s%-12s%-10s\n", "β", "kinky bits", "s=0.25 bits", "savings")
for beta in BETA_VALUES
    Nk = Result[(beta, 0.0,  (0.0, 0.0))].min_N_for_eps[1e-4]
    Ns = Result[(beta, 0.25, (0.0, 0.0))].min_N_for_eps[1e-4]
    bk = Nk > 0 ? Int(log2(Nk)) : -1
    bs = Ns > 0 ? Int(log2(Ns)) : -1
    @printf("%-8.1f%-12s%-12s%-10s\n", beta,
            bk < 0 ? ">$(Int(log2(maximum(N_GRID))))" : "$(bk)",
            bs < 0 ? ">$(Int(log2(maximum(N_GRID))))" : "$(bs)",
            (bk < 0 || bs < 0) ? "—" : "$(bk - bs) bits")
end
println()
println("Done.")
