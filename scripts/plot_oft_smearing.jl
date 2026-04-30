#!/usr/bin/env julia
#
# Plot 4: OFT Gaussian smearing around two Bohr frequencies.
#
# The OFT \hat{A}(\omega) = \sum_{\nu \in B_H} A_\nu \hat{f}(\omega - \nu)
# replaces each delta-peaked Davies jump A_\nu by a Gaussian-smeared bump
# centred on the Bohr frequency \nu. The width of that bump is set by the
# filter's Fourier transform (2_methods.tex, eq. line 51):
#
#     \hat{f}(\omega) = (2 \pi \sigma^2)^{-1/4} \exp(-\omega^2 / (4 \sigma^2)).
#
# With \sigma = 1/\beta (Chen et al. canonical, eq. line 160), nearby Bohr
# frequencies smear into one another: that overlap is exactly the source of
# the Gaussian smearing error in \alpha_{\nu_1, \nu_2} (eq. line 116).
#
# This figure picks two arbitrary Bohr frequencies \nu_1 < \nu_2 separated by
# about 3 \sigma so that the L²-normalized filters \hat{f}(\omega - \nu_i)
# visibly overlap, then shades the overlap region to highlight the ambiguous
# energies.
#
# PHYSICS CHECK: \beta = 10, \sigma = 1/\beta = 0.1, identical to
# scripts/plot_transition_weights.jl. The factors in \hat{f} are chosen so
# that \sigma is the std dev of |\hat{f}|^2, equivalently |\hat{f}|^2 \propto
# \exp(-\omega^2 / (2\sigma^2)). For \hat{f} itself, \exp(-(\omega/(2\sigma))^2)
# means the half-width at 1/e of the peak is exactly 2\sigma — the "inverse
# width" of the time-domain Gaussian referenced in 2_methods.tex line 49–52.
# We pick \nu_1 = -0.10, \nu_2 = 0.30 (separation 0.40 = 4\sigma = the full
# width at 1/e of one filter), which produces a visible but mild overlap.
#
# PHYSICS CHECK: We plot the actual L²-normalized \hat{f} (no rescaling), so
# the peak value (2 \pi \sigma^2)^{-1/4} \approx 2.0 reflects the genuine
# amplitude that enters the OFT. The y-axis is labelled \hat{f}(\omega - \nu).
#
# Output: drafts/plots/oft_smearing.{pdf,png}
#
# Usage:  julia --project scripts/plot_oft_smearing.jl

using Printf
using Plots

# ── Parameters ────────────────────────────────────────────────────────────────
const β = 10.0                # inverse temperature (matches plot_transition_weights.jl)
const σ = 1 / β               # filter width (Chen et al. canonical, σ = 1/β)

# Two arbitrary Bohr frequencies. Separation ≈ 3σ ≈ 0.9·FWHM → visible overlap.
const ν1 = -0.10
const ν2 =  0.20

# Plot grid (zoomed to the relevant region around both peaks)
const ω_grid = range(-0.6, 0.6; length=2001)

# ── OFT Gaussian filter ──────────────────────────────────────────────────────
# \hat{f}(\omega) = (2\pi\sigma^2)^{-1/4} \exp(-\omega^2 / (4\sigma^2))
fhat(ω; σ=σ) = (2π * σ^2)^(-1/4) * exp(-ω^2 / (4 * σ^2))

g1 = fhat.(ω_grid .- ν1)
g2 = fhat.(ω_grid .- ν2)
g_overlap = min.(g1, g2)         # pointwise minimum: the visually overlapping mass

peak = (2π * σ^2)^(-1/4)
half_e_width = 2 * σ            # half-width at 1/e of the peak (thesis convention)
sep  = ν2 - ν1
@printf("Parameters: β = %.2f, σ = %.4f, peak \\hat{f}(0) = %.4f, half-width at 1/e = 2σ = %.4f\n",
        β, σ, peak, half_e_width)
@printf("Bohr frequencies: ν₁ = %+.3f, ν₂ = %+.3f, separation = %.3f (= %.2f·σ = %.2f·(2σ))\n",
        ν1, ν2, sep, sep/σ, sep/half_e_width)
@printf("Filter values at the midpoint ω = %.3f:\n", (ν1+ν2)/2)
@printf("  \\hat{f}(ω - ν₁) = %.4f\n", fhat((ν1+ν2)/2 - ν1))
@printf("  \\hat{f}(ω - ν₂) = %.4f   (equal by symmetry — peak ambiguity)\n",
        fhat((ν1+ν2)/2 - ν2))
@printf("Overlap area ∫ min(g₁, g₂) dω ≈ %.4f\n",
        sum(g_overlap) * step(ω_grid))

# ── Plot styling ──────────────────────────────────────────────────────────────
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
    size             = (700, 460),
    dpi              = 200,
    margin           = 4Plots.mm,
)

# Thesis colour palette (named).
# Two distinct hues — bordeaux (warm) and pinegreen (cool) — keep the two
# Bohr-frequency Gaussians easy to tell apart while staying within the thesis
# palette. The overlap is shaded with deepplum, a darker neutral that reads
# as "where the two coexist".
c_g1      = "#7A2E39"   # bordeaux
c_g2      = "#2D5A3D"   # pinegreen
c_overlap = "#735874"   # deepplum (overlap shading)
c_axis    = :grey40

ymax = 1.10 * peak

plt = plot(
    ω_grid, g1;
    label  = "\$\\hat{f}(\\omega - \\nu_1)\$",
    color  = c_g1,
    fillrange = 0,
    fillalpha = 0.18,
    fillcolor = c_g1,
    xlabel = "\$\\omega\$",
    ylabel = "\$\\hat{f}(\\omega - \\nu)\$",
    title  = "OFT Gaussian filtering around Bohr frequencies  (\$\\beta=$(Int(β)),\\; \\sigma = 1/\\beta\$)",
    legend                 = :topleft,
    legendfontsize         = 7,
    foreground_color_legend = nothing,
    background_color_legend = RGBA(1, 1, 1, 0.7),
    xlims  = (first(ω_grid), last(ω_grid)),
    ylims  = (-0.05, ymax),
)

plot!(plt, ω_grid, g2;
    label = "\$\\hat{f}(\\omega - \\nu_2)\$",
    color = c_g2,
    fillrange = 0,
    fillalpha = 0.18,
    fillcolor = c_g2)

# Overlap region: pointwise minimum of the two filters, shaded with deepplum
# to mark the energies that cannot be cleanly attributed to ν₁ or ν₂.
plot!(plt, ω_grid, g_overlap;
    label = "overlap",
    color = c_overlap,
    linewidth = 1.4,
    fillrange = 0,
    fillalpha = 0.45,
    fillcolor = c_overlap)

# Mark the two Bohr frequencies on the x-axis with small spikes + labels.
vline!(plt, [ν1];
    color     = c_g1,
    linestyle = :dash,
    linewidth = 1.2,
    label     = "")
vline!(plt, [ν2];
    color     = c_g2,
    linestyle = :dash,
    linewidth = 1.2,
    label     = "")

# GR's math renderer ignores \boldsymbol/\mathbf, so we get visual prominence
# from a larger font size instead of true bold (no ghosting).
annotate!(plt, ν1 + 0.012, 1.5, text("\$\\nu_1\$", c_g1, 16, :left, :bottom))
annotate!(plt, ν2 + 0.012, 1.5, text("\$\\nu_2\$", c_g2, 16, :left, :bottom))

# Width annotation: an arrow at height \hat{f}(0)/e showing the full width
# at 1/e of one Gaussian (= 4σ, the natural width in the thesis convention).
e_height = peak / ℯ
left_e  = ν1 - 2 * σ
right_e = ν1 + 2 * σ
plot!(plt, [left_e, right_e], [e_height, e_height];
    color = c_axis,
    linewidth = 1.0,
    arrow = true,
    label = "")
plot!(plt, [right_e, left_e], [e_height, e_height];
    color = c_axis,
    linewidth = 1.0,
    arrow = true,
    label = "")
# Place the 4/β label just right of ν₁'s dashed line so it doesn't sit on it.
annotate!(plt, ν1 + 0.04, e_height + 0.08, text("\$4/\\beta\$",
    c_axis, 10, :center))

# ── Save ──────────────────────────────────────────────────────────────────────
out_pdf = joinpath(@__DIR__, "..", "drafts", "plots", "oft_smearing.pdf")
out_png = joinpath(@__DIR__, "..", "drafts", "plots", "oft_smearing.png")
mkpath(dirname(out_pdf))
savefig(plt, out_pdf)
savefig(plt, out_png)
@printf("Saved %s\n", relpath(out_pdf, joinpath(@__DIR__, "..")))
@printf("Saved %s\n", relpath(out_png, joinpath(@__DIR__, "..")))
