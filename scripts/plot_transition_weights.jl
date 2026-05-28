#!/usr/bin/env julia
#
# Plot 1: Transition weights γ(ω) for the four CKG variants used in the thesis.
#
#   • Gaussian γ_G(ω)            = exp( -(ω + ω_γ)² / (2 σ_γ²) )
#   • Kinky Metropolis γ_M(ω)    = exp( -β · max(ω + βσ²/2, 0) )
#   • Smooth Metropolis γ_M^{(s)}(ω) = γ_M(ω) · ½ [erfc(z₋) + e^{β|ω̃|} erfc(z₊)]
#       with ω̃ = ω + βσ²/2  and  z± = (βσ√s)/(2√2) ± |ω̃|/(σ√(2s))
#   • Glauber γ_Gl(ω)            = 1 / (1 + exp(β · (ω + βσ²/2)))
#       (Fermi-function form, shifted by the same βσ²/2 as γ_M so the
#        inflection point of γ_Gl coincides with the kink of γ_M)
#
# References (supplementary-informations/2_methods.tex):
#   eq:gaussian-transition (line 81)
#   eq:shifted-metro-func (line 228)
#   eq:smooth-metro       (line 296)
#   Glauber as the s = 2 limit of γ_M^{(s)} per Chen et al. \cite{chen2023quantum}
#
# PHYSICS CHECK: σ = ω_γ = σ_γ = 1/β̄ follows Chen et al., where β̄ = β_alg is
# the algorithm-space inverse temperature on the rescaled Bohr grid. With
# σ = 1/β̄ the Gaussian Kossakowski matrix is skew-symmetric (Eq. line 160).
# The Gaussian transition centre at ω = -1/β̄ and the Metropolis kink at
# -1/(2 β̄) sit in the same energy scale.
#
# This plot illustrates the *shapes* of γ(ω) on the rescaled-ω axis the
# algorithm sees natively. The physical β_phys is fixture-dependent
# (β_phys = β̄ / rescaling_factor); for these qualitative shape comparisons
# it is more distracting than helpful, so we just pick β̄ = 10 — the legacy
# value — and label the plot accordingly. Kossakowski heatmaps still
# parametrise by β_phys, where the physical temperature is the natural
# axis.
#
# PHYSICS CHECK: s = 0.25 is the chosen smoothing parameter for γ_M^{(s)}.
#   - s > 0 puts γ_M^{(s)} in the Gevrey-1/2 class (Prop. smooth-metro-gevrey),
#     giving quadrature errors exp(-c√N) instead of the kinky O(1/N).
#   - s = 2 is the Glauber-like extreme (worse mixing per the thesis).
#   - Production scripts use s = 0.4.
#   - 0.25 is small enough that the spectral gap stays close to the kinky
#     (Peskun-Tierney) optimum and large enough to materially smoothen the kink.
#
# Output: drafts/plots/transition_weights.pdf  (and .png for quick preview)
#
# Usage:  julia --project scripts/plot_transition_weights.jl

using Printf
using SpecialFunctions: erfc, erfcx

using Plots

# ── Parameters ────────────────────────────────────────────────────────────────
# β̄ = β_alg is the algorithm-space inverse temperature on the rescaled Bohr
# grid; we use the legacy value 10 here because the curves' qualitative shape
# is what matters for this plot (Kossakowski heatmaps stay on β_phys).
const β      = 10.0
const σ      = 1 / β                      # filter width on the rescaled Bohr grid
const ω_γ    = 1 / β                      # Gaussian transition centre  (Chen et al.)
const σ_γ    = 1 / β                      # Gaussian transition width   (Chen et al.)
const s      = 0.25                       # smoothing parameter for γ_M^{(s)}
const ω_kink = -β * σ^2 / 2               # kink position of γ_M (and γ_Gl inflection)

@printf("β̄ (= β_alg) = %.3f   σ = 1/β̄ = %.4f   ω_kink = %.4f\n", β, σ, ω_kink)

# ── Transition weights ────────────────────────────────────────────────────────
γ_gauss(ω) = exp(-(ω + ω_γ)^2 / (2 * σ_γ^2))

γ_kinky(ω) = exp(-β * max(ω + β * σ^2 / 2, 0.0))

# Glauber: Fermi-function form, shifted by βσ²/2 to align inflection with γ_M kink.
# Numerically stable for both signs of the exponent.
function γ_glauber(ω; β=β, σ=σ)
    x = β * (ω + β * σ^2 / 2)
    return x ≥ 0 ? exp(-x) / (1 + exp(-x)) : 1 / (1 + exp(x))
end

# Robust smooth Metropolis: the second term inside the bracket can overflow for
# ω̃ ≫ 0 because both e^{β|ω̃|} and erfc(z₊) are extreme. Use the scaled
# complementary error function erfcx(z) = e^{z²} erfc(z) and combine exponents.
function γ_smooth(ω; β=β, σ=σ, s=s)
    ω̃ = ω + β * σ^2 / 2
    aω̃ = abs(ω̃)
    zp = β * σ * sqrt(s) / (2 * sqrt(2)) + aω̃ / (σ * sqrt(2 * s))
    zm = β * σ * sqrt(s) / (2 * sqrt(2)) - aω̃ / (σ * sqrt(2 * s))
    # e^{β|ω̃|} erfc(z₊) = e^{β|ω̃| - z₊²} · erfcx(z₊)   (no overflow)
    log_factor = β * aω̃ - zp^2
    second = exp(log_factor) * erfcx(zp)
    bracket = (erfc(zm) + second) / 2
    return γ_kinky(ω) * bracket
end

# ── Grid ──────────────────────────────────────────────────────────────────────
# Rescaled-ω grid covering the n=5 Bohr range [-0.9, 0.9]; curves decay to ≈ 0
# by |ω| ≈ 0.3 at β_alg ≈ 22, so [-0.5, 0.5] is a tight visualisation window.
ω_grid = range(-0.5, 0.5; length=1601)

g_gauss   = γ_gauss.(ω_grid)
g_kinky   = γ_kinky.(ω_grid)
g_smooth  = γ_smooth.(ω_grid)
g_glauber = γ_glauber.(ω_grid)

# Sanity: smooth ≤ kinky pointwise (Prop 5 corollary)
@assert all(g_smooth .<= g_kinky .+ 1e-12) "smooth Metropolis exceeds kinky — formula bug"
# Sanity: γ_Gl(ω̃ = 0) = 1/2
@assert abs(γ_glauber(ω_kink) - 0.5) < 1e-12 "Glauber inflection should give 1/2 at ω = -βσ²/2"

# Print a few diagnostic values
@printf("Parameters: β = %.2f, σ = %.4f, s = %.2f, kink at ω = %.4f\n", β, σ, s, ω_kink)
@printf("γ values at ω = 0 (above kink):\n")
@printf("  γ_G       = %.4f\n", γ_gauss(0.0))
@printf("  γ_M       = %.4f\n", γ_kinky(0.0))
@printf("  γ_M^{(s)} = %.4f   (factor: %.3f)\n", γ_smooth(0.0), γ_smooth(0.0) / γ_kinky(0.0))
@printf("  γ_Gl      = %.4f\n", γ_glauber(0.0))
@printf("γ values at the kink ω = %.4f:\n", ω_kink)
@printf("  γ_M       = %.4f\n", γ_kinky(ω_kink))
@printf("  γ_M^{(s)} = %.4f   (factor: %.3f)\n", γ_smooth(ω_kink), γ_smooth(ω_kink) / γ_kinky(ω_kink))
@printf("  γ_Gl      = %.4f\n", γ_glauber(ω_kink))

# ── Plot ──────────────────────────────────────────────────────────────────────
default(
    fontfamily       = "Computer Modern",
    titlefontsize    = 12,
    guidefontsize    = 12,
    tickfontsize     = 10,
    legendfontsize   = 8,
    framestyle       = :box,
    grid             = true,
    gridalpha        = 0.25,
    linewidth        = 2.0,
    dpi              = 200,
    margin           = 4Plots.mm,
)

# Thesis colour palette (named).
# Primary curves: γ_M = bordeaux, γ_M^{(s)} = deepplum.
# Secondary curves use the "brighter variant" recipe (HSV V-boost to 0.80,
# preserves hue and saturation) of slateblue (Gaussian) and sage (Glauber).
c_kinky   = "#7A2E39"   # bordeaux
c_smooth  = "#735874"   # deepplum
c_gauss   = "#7FA4CC"   # slateblue, brighter variant of #5C7794
c_glauber = "#B2CCA1"   # sage, brighter variant of #8B9F7E

# Left panel — Gaussian only.
plt_gauss = plot(
    ω_grid, g_gauss;
    label  = "\$\\gamma_G\$ (Gaussian)",
    color  = c_gauss,
    xlabel = "\$\\omega\$",
    ylabel = "\$\\gamma(\\omega)\$",
    title  = "Gaussian  (\$\\bar{\\beta}=$(Int(β))\$)",
    legend                 = :topright,
    legendfontsize         = 8,
    foreground_color_legend = nothing,
    background_color_legend = RGBA(1, 1, 1, 0.7),
    xlims  = (-0.5, 0.5),
    ylims  = (-0.02, 1.08),
)
vline!(plt_gauss, [ω_kink];
    color     = :grey40,
    linestyle = :dash,
    linewidth = 1.2,
    label     = "")
annotate!(plt_gauss, ω_kink + 0.012, 1.04,
          text("\$\\omega = -\\beta\\sigma^{2}/2\$", :grey40, 8, :left, rotation=0))

# Right panel — Metropolis (kinky + smooth) and Glauber.
plt_metro = plot(
    ω_grid, g_kinky;
    label  = "\$\\gamma_M\$ (kinky Metropolis)",
    color  = c_kinky,
    xlabel = "\$\\omega\$",
    ylabel = "\$\\gamma(\\omega)\$",
    title  = "Metropolis & Glauber  (\$\\bar{\\beta}=$(Int(β))\$)",
    legend                 = :topright,
    legendfontsize         = 8,
    foreground_color_legend = nothing,
    background_color_legend = RGBA(1, 1, 1, 0.7),
    xlims  = (-0.5, 0.5),
    ylims  = (-0.02, 1.08),
)
plot!(plt_metro, ω_grid, g_smooth;
    label = "\$\\gamma_M^{(s)}\$ (smooth, \$s=$(s)\$)",
    color = c_smooth)
plot!(plt_metro, ω_grid, g_glauber;
    label = "\$\\gamma_{\\mathrm{Gl}}\$ (Glauber)",
    color = c_glauber)
vline!(plt_metro, [ω_kink];
    color     = :grey40,
    linestyle = :dash,
    linewidth = 1.2,
    label     = "")
annotate!(plt_metro, ω_kink + 0.012, 1.04,
          text("\$\\omega = -\\beta\\sigma^{2}/2\$", :grey40, 8, :left, rotation=0))

# Combined 1×2 figure.
plt = plot(plt_gauss, plt_metro;
           layout = (1, 2),
           size   = (1300, 460),
           link   = :all)

# ── Save ──────────────────────────────────────────────────────────────────────
out_pdf = joinpath(@__DIR__, "..", "drafts", "plots", "transition_weights.pdf")
out_png = joinpath(@__DIR__, "..", "drafts", "plots", "transition_weights.png")
mkpath(dirname(out_pdf))
savefig(plt, out_pdf)
savefig(plt, out_png)
@printf("Saved %s\n", relpath(out_pdf, joinpath(@__DIR__, "..")))
@printf("Saved %s\n", relpath(out_png, joinpath(@__DIR__, "..")))
