#!/usr/bin/env julia
#
# Plot 2: time-domain coherent-term kernels b_+(t) and b_-(t) for the
# transition weights used in the thesis.
#
# Three figures are produced:
#   • drafts/plots/b_plus.pdf            -- |b_+(t)| for Gaussian + kinky/smooth η-reg Metro
#   • drafts/plots/b_minus.pdf           -- b_-(t)   (universal across γ choices)
#   • drafts/plots/b_plus_with_a_smooth.pdf -- adds the a-regularized "smooth Metro"
#                                          (iε-style; the failed-appendix attempt)
#
# b_+(t) is γ-dependent and complex-valued; we plot |b_+(t)|:
#   • Gaussian              b_+^{(G)}(t)     = (β σ_γ / π√π)·exp(-2 β ω_γ (2t² + i·t))
#                           (eq:b_plus-gauss, line 513 of 2_methods.tex)
#   • Kinky  Metro (η-reg)  b_+^{(0,η)}(t)   = eq:b_plus-s-eta with s = 0
#   • Smooth Metro (η-reg)  b_+^{(s,η)}(t)   = eq:b_plus-s-eta with s > 0
#                           (eq:b_plus-s-eta, line 524 of 2_methods.tex)
#   • Smooth Metro (a-reg)  b_+^{(a,s)}(t)   = √(4a+1)/(√2 π²) · e^{-a s/2}
#                                              · e^{-σ²β² t(2t+i)(1+s)} / (4t² + a + 2it)
#                           — the iε-style regularization from line 549 of 2_methods.tex.
#                           Both poles shift off the real axis (no t=η jump), at the
#                           price of an e^{-a/2}-like attenuation near the kink.
#                           Mirrors src/coherent.jl::_compute_b_plus_smooth.
#
# Both η-Metropolis variants use the η-regularization that replaces the 1/t pole
# on |t| ≤ η with a smooth bounded kernel (Chen et al., Prop. B.1). The kernel
# is therefore bounded everywhere but has a finite jump at |t| = η. The
# a-regularized variant has no jump anywhere.
#
# b_-(t) (universal, real-valued) is identical for all γ choices:
#   b_-(t) = (2√π / (βσ))·e^{β²σ²/8} ·[(1/cosh(2π t/(βσ))) ∗_t (sin(-βσ t) e^{-2t²})]
#   (eq:b_minus, line 503 of 2_methods.tex)
#
# Convention: t is the natural time variable from the inner integral
# ∫ b_+(t')A^{a†}(βt')A^a(-βt') dt' (eq:coherent-time-domain). For the canonical
# Chen choice σ = ω_γ = σ_γ = 1/β we have σβ = 1 and both kernels live on the
# O(1) time scale.
#
# PHYSICS CHECK: η = 0.05 is the test-suite default (test_helpers.jl) and is
# small enough that the regularized window |t| ≤ η is sub-thesis-relevant for
# the kernel decay scale (~1), but large enough to keep the t = ±η jump
# visible on a thesis figure.
#
# PHYSICS CHECK: s = 0.25 matches plot_transition_weights.jl for visual
# consistency with Plot 1 of the epic.
#
# PHYSICS CHECK: a = β/30 ≈ 0.333 is the canonical test-helper default for
# the a-regularized variant (test_helpers.jl::make_config). With β = 10 and
# a = 1/3, the iε-shifted poles sit at t ≈ +i·0.135 and t ≈ -i·0.635 — both
# safely off the real axis.
#
# Output: drafts/plots/b_plus.{pdf,png}, drafts/plots/b_minus.{pdf,png}
#
# Usage: julia --project scripts/plot_b_plus_b_minus.jl

using Printf
using QuadGK: quadgk
using Plots

# ── Parameters ────────────────────────────────────────────────────────────────
const β        = 10.0
const σ        = 1 / β
const ω_γ      = 1 / β
const σ_γ      = 1 / β
const s_smooth = 0.25
const η        = 0.05
const a_reg    = β / 30      # canonical a-regularization value (test_helpers.jl)

# ── b_+ kernels (γ-dependent, complex-valued) ─────────────────────────────────
# Mirror of src/coherent.jl::_compute_b_plus
b_plus_gauss(t) = β * σ_γ * exp(-2 * β * ω_γ * (2 * t^2 + 1im * t)) / sqrt(pi^3)

# Mirror of src/coherent.jl::_compute_b_plus_metro (s = 0 ↔ kinky, s > 0 ↔ smooth)
function b_plus_metro(t::Real; s::Real=0.0)
    if abs(t) < 1e-12
        # L'Hôpital limit at t = 0; reduces to 1/(2√2 π²) for σβ = 1, s = 0
        return complex((2 - σ^2 * β^2 * (1 + s)) / (2 * sqrt(2) * pi^2))
    elseif abs(t) <= η
        numerator = exp(-σ^2 * β^2 * (2 * t^2 + 1im * t) * (1 + s)) + 1im * (2 * t + 1im)
    else
        numerator = exp(-σ^2 * β^2 * (2 * t^2 + 1im * t) * (1 + s))
    end
    denominator = t * (2 * t + 1im)
    return (1 / (2 * sqrt(2) * pi^2)) * numerator / denominator
end

b_plus_metro_s0(t) = b_plus_metro(t; s=0.0)        # Metropolis,        η-reg
b_plus_smooth(t)   = b_plus_metro(t; s=s_smooth)   # smooth Metropolis, η-reg

# Mirror of src/coherent.jl::_compute_b_plus_smooth — a-regularized "iε-style"
# smooth Metro. Both poles of 1/(4t² + a + 2it) sit on the imaginary axis for
# any a > 0, so the function is C^∞ on the real line (no |t|=η jump).
function b_plus_a_smooth(t::Real; a::Real=a_reg, s::Real=s_smooth)
    num = exp(-a * s / 2) * exp(-σ^2 * β^2 * t * (2 * t + 1im) * (1 + s))
    den = 4 * t^2 + a + 2im * t
    return sqrt(4 * a + 1) * num / (sqrt(2) * pi^2 * den)
end

# ── b_- kernel (universal, real-valued) ───────────────────────────────────────
# Mirror of src/coherent.jl::_compute_b_minus + ::_convolute
function b_minus(t::Real; atol=1e-12, rtol=1e-12)
    f1(s) = 1 / cosh(2 * pi * s / (β * σ))
    f2(u) = sin(-u * β * σ) * exp(-2 * u^2)
    integrand(s) = f1(s) * f2(t - s)
    conv, _ = quadgk(integrand, -Inf, Inf; atol=atol, rtol=rtol)
    return 2 * sqrt(pi) * exp(β^2 * σ^2 / 8) * conv / (β * σ)
end

# ── Grids ─────────────────────────────────────────────────────────────────────
# Dense base grid for b_+, plus extra resolution around |t| = η for the jump.
t_grid_plus = sort(unique(vcat(
    range(-2.0, 2.0; length=2001),
    range(-2η, 2η; length=401),
    [0.0],
)))

# b_- has effective width O(βσ) = 1 in our config; [-3, 3] captures decay.
t_grid_minus = collect(range(-3.0, 3.0; length=801))

# ── Compute kernel values ─────────────────────────────────────────────────────
@info "Computing |b_+(t)| on a $(length(t_grid_plus))-point grid…"
g_p_gauss     = b_plus_gauss.(t_grid_plus)
g_p_metro_s0  = b_plus_metro_s0.(t_grid_plus)
g_p_smooth    = b_plus_smooth.(t_grid_plus)
g_p_a_smooth  = b_plus_a_smooth.(t_grid_plus)

@info "Computing b_-(t) on a $(length(t_grid_minus))-point grid (quadgk)…"
g_m = b_minus.(t_grid_minus)

# Sanity: b_- is real (real integrand → real result, up to QuadGK roundoff)
@assert all(abs(imag(z)) < 1e-10 for z in g_m) "b_- should be real-valued"
g_m_real = real.(g_m)

# Sanity: |b_+|(0) matches the analytic L'Hôpital limit
zero_idx_p = findfirst(==(0.0), t_grid_plus)
@assert zero_idx_p !== nothing "t = 0 missing from b_+ grid"
analytic_metro_at_0    = abs((2 - σ^2 * β^2) / (2 * sqrt(2) * pi^2))
analytic_smooth_at_0   = abs((2 - σ^2 * β^2 * (1 + s_smooth)) / (2 * sqrt(2) * pi^2))
analytic_a_smooth_at_0 = abs(sqrt(4 * a_reg + 1) * exp(-a_reg * s_smooth / 2) /
                              (sqrt(2) * pi^2 * a_reg))
@assert abs(abs(g_p_metro_s0[zero_idx_p]) - analytic_metro_at_0)    < 1e-12
@assert abs(abs(g_p_smooth[zero_idx_p])   - analytic_smooth_at_0)   < 1e-12
@assert abs(abs(g_p_a_smooth[zero_idx_p]) - analytic_a_smooth_at_0) < 1e-12

# Sanity: b_-(0) = 0 by parity (even ∗ odd = odd)
@assert abs(g_m_real[findfirst(==(0.0), t_grid_minus)]) < 1e-9 "b_-(0) should vanish"

# ── L¹ norms via trapezoidal rule on the (possibly non-uniform) grid ──────────
function trapz(ts::AbstractVector, fs::AbstractVector)
    @assert length(ts) == length(fs)
    s = zero(promote_type(eltype(ts), eltype(fs)))
    @inbounds for i in 1:length(ts) - 1
        s += 0.5 * (fs[i] + fs[i + 1]) * (ts[i + 1] - ts[i])
    end
    return s
end

norm_b_plus_gauss    = trapz(t_grid_plus,  abs.(g_p_gauss))
norm_b_plus_metro_s0 = trapz(t_grid_plus,  abs.(g_p_metro_s0))
norm_b_plus_smooth   = trapz(t_grid_plus,  abs.(g_p_smooth))
norm_b_minus         = trapz(t_grid_minus, abs.(g_m_real))

# ── Diagnostics ───────────────────────────────────────────────────────────────
@printf("Parameters: β = %.2f, σ = %.4f, σβ = %.2f, s = %.2f, η = %.3f\n",
        β, σ, β * σ, s_smooth, η)
@printf("|b_+|(0):\n")
@printf("  Gaussian:                     %.6f   (analytic β σ_γ / π√π             = %.6f)\n",
        abs(g_p_gauss[zero_idx_p]), β * σ_γ / (pi * sqrt(pi)))
@printf("  Metropolis        (η-reg):    %.6f   (analytic 1/(2√2 π²)              = %.6f)\n",
        abs(g_p_metro_s0[zero_idx_p]), analytic_metro_at_0)
@printf("  smooth Metropolis (η-reg):    %.6f   (analytic (2-σ²β²(1+s))/(2√2 π²)  = %.6f)\n",
        abs(g_p_smooth[zero_idx_p]), analytic_smooth_at_0)
@printf("  smooth Metropolis (a-reg, a=%.3f): %.6f   (analytic √(4a+1) e^{-as/2}/(√2 π² a) = %.6f)\n",
        a_reg, abs(g_p_a_smooth[zero_idx_p]), analytic_a_smooth_at_0)
@printf("|b_+| just outside |t|=η (jump):\n")
i_just_out = findfirst(t -> t >= η + 1e-4, t_grid_plus)
@printf("  Metropolis        at t=%.4f:  %.4f\n",
        t_grid_plus[i_just_out], abs(g_p_metro_s0[i_just_out]))
@printf("  smooth Metropolis at t=%.4f:  %.4f\n",
        t_grid_plus[i_just_out], abs(g_p_smooth[i_just_out]))
@printf("L¹ norms (trapezoidal):\n")
@printf("  ‖b_+^{(G)}‖_1     = %.4f\n", norm_b_plus_gauss)
@printf("  ‖b_+^{(0,η)}‖_1   = %.4f   (Chen bound: 2/(5π^{3/2}) + log(1/η)/(√2 π²) = %.4f)\n",
        norm_b_plus_metro_s0,
        2 / (5 * pi^1.5) + log(1 / η) / (sqrt(2) * pi^2))
@printf("  ‖b_+^{(s,η)}‖_1   = %.4f   (Chen bound: previous - log(1+s)/(2√2 π²)   = %.4f)\n",
        norm_b_plus_smooth,
        2 / (5 * pi^1.5) + log(1 / η) / (sqrt(2) * pi^2) -
            log(1 + s_smooth) / (2 * sqrt(2) * pi^2))
@printf("  ‖b_-‖_1           = %.4f   (Chen Hölder bound: ‖b_-‖_1 ≤ 1; no closed form)\n",
        norm_b_minus)

# ── Plot defaults (match plot_transition_weights.jl) ──────────────────────────
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

# Thesis colour palette (matches plot_transition_weights.jl).
c_kinky    = "#7A2E39"   # bordeaux
c_smooth   = "#735874"   # deepplum
c_gauss    = "#7FA4CC"   # slateblue, brighter variant of #5C7794
c_a_smooth = "#B5654A"   # terracotta — sets the a-regularized variant apart
c_minus    = "#2D5A3D"   # pinegreen — single neutral colour for the universal b_-

# Helper: split a curve at the |t| = η discontinuity so the discontinuous Metro
# kernels do not visually connect across the jump.
function split_at_eta(ts::AbstractVector, ys::AbstractVector)
    idx_neg   = findall(t -> t <  -η, ts)
    idx_inner = findall(t -> abs(t) <= η, ts)
    idx_pos   = findall(t -> t >   η, ts)
    return (idx_neg, idx_inner, idx_pos)
end

# Helper: draw a styled white-fill / grey-border box at (xb1..xb2, yb1..yb2)
# and place each `lines[i]` (LaTeX-rendered text) on its own row, centred.
function add_norm_box!(plt, xb1, xb2, yb1, yb2, lines::Vector{<:AbstractString};
                       fontsize::Int=8)
    plot!(plt, Plots.Shape([xb1, xb2, xb2, xb1], [yb1, yb1, yb2, yb2]);
          fillcolor = RGBA(1, 1, 1, 0.85),
          linecolor = :grey60,
          linewidth = 0.6,
          label     = "")
    n = length(lines)
    height = yb2 - yb1
    for (i, line) in enumerate(lines)
        # Vertical centre of the i-th row (top to bottom).
        y = yb2 - (i - 0.5) * height / n
        annotate!(plt, (xb1 + xb2) / 2, y, text(line, fontsize, :center))
    end
    return plt
end

# ── Plot 1: |b_+(t)| ──────────────────────────────────────────────────────────
plt_p = plot(
    t_grid_plus, abs.(g_p_gauss);
    label  = "\$|b_+^{(G)}|\$ (Gaussian)",
    color  = c_gauss,
    xlabel = "\$t\$",
    ylabel = "\$|b_+(t)|\$",
    title  = "Coherent-term kernel  \$|b_+(t)|\$  (\$\\bar{\\beta}=$(Int(β)),\\; \\sigma = 1/\\bar{\\beta},\\; \\eta=$(η)\$)",
    legend                  = :topright,
    foreground_color_legend = nothing,
    background_color_legend = RGBA(1, 1, 1, 0.7),
    xlims  = (-2.0, 2.0),
    ylims  = (-0.02, 0.85),
)

# Metropolis (s = 0, η-regularized) — three segments split at ±η so the line
# breaks at the jump.
i_neg_k, i_in_k, i_pos_k = split_at_eta(t_grid_plus, abs.(g_p_metro_s0))
plot!(plt_p, t_grid_plus[i_neg_k], abs.(g_p_metro_s0[i_neg_k]);
      label = "\$|b_+^{(0,\\eta)}|\$ (Metropolis)",
      color = c_kinky)
plot!(plt_p, t_grid_plus[i_in_k],  abs.(g_p_metro_s0[i_in_k]);  label = "", color = c_kinky)
plot!(plt_p, t_grid_plus[i_pos_k], abs.(g_p_metro_s0[i_pos_k]); label = "", color = c_kinky)

# Smooth Metropolis (s > 0, η-regularized) — three segments split at ±η.
i_neg_s, i_in_s, i_pos_s = split_at_eta(t_grid_plus, abs.(g_p_smooth))
plot!(plt_p, t_grid_plus[i_neg_s], abs.(g_p_smooth[i_neg_s]);
      label = "\$|b_+^{(s,\\eta)}|\$ (smooth Metropolis, \$s=$(s_smooth)\$)",
      color = c_smooth)
plot!(plt_p, t_grid_plus[i_in_s],  abs.(g_p_smooth[i_in_s]);  label = "", color = c_smooth)
plot!(plt_p, t_grid_plus[i_pos_s], abs.(g_p_smooth[i_pos_s]); label = "", color = c_smooth)

# Mark the η boundaries
vline!(plt_p, [-η, η];
       color     = :grey40,
       linestyle = :dash,
       linewidth = 1.0,
       label     = "")
annotate!(plt_p, η + 0.04, 0.81, text("\$|t|=\\eta\$", :grey40, 8, :left))

# L¹-norm side panel, placed below the topright legend.
add_norm_box!(plt_p,
    0.78, 1.95,        # x extent
    0.36, 0.58,        # y extent (just below the 3-line legend)
    [
        "\$\\Vert b_+^{(G)}\\Vert_1 = $(@sprintf("%.3f", norm_b_plus_gauss))\$",
        "\$\\Vert b_+^{(0,\\eta)}\\Vert_1 = $(@sprintf("%.3f", norm_b_plus_metro_s0))\$",
        "\$\\Vert b_+^{(s,\\eta)}\\Vert_1 = $(@sprintf("%.3f", norm_b_plus_smooth))\$",
    ];
    fontsize = 8,
)

out_pdf_p = joinpath(@__DIR__, "..", "drafts", "plots", "b_plus.pdf")
out_png_p = joinpath(@__DIR__, "..", "drafts", "plots", "b_plus.png")
mkpath(dirname(out_pdf_p))
savefig(plt_p, out_pdf_p)
savefig(plt_p, out_png_p)
@printf("Saved %s\n", relpath(out_pdf_p, joinpath(@__DIR__, "..")))
@printf("Saved %s\n", relpath(out_png_p, joinpath(@__DIR__, "..")))

# ── Plot 1b: |b_+(t)| including the a-regularized smooth Metropolis ──────────
# The a-regularization (a > 0) corresponds to the iε-prescription regularization
# from the failed-attempt appendix of the thesis (line 549 of 2_methods.tex):
# both poles of 1/(4t² + a + 2it) sit on the imaginary axis, so the kernel is
# C^∞ on the real line. Visually: no |t|=η jump, a smooth single-peak shape
# that bridges what the η-regularized kernel achieves only piecewise.
plt_p2 = plot(
    t_grid_plus, abs.(g_p_gauss);
    label  = "\$|b_+^{(G)}|\$ (Gaussian)",
    color  = c_gauss,
    xlabel = "\$t\$",
    ylabel = "\$|b_+(t)|\$",
    title  = "\$|b_+(t)|\$  with a-regularized variant  (\$\\bar{\\beta}=$(Int(β)),\\; \\sigma=1/\\bar{\\beta},\\; s=$(s_smooth)\$)",
    legend                  = :topright,
    foreground_color_legend = nothing,
    background_color_legend = RGBA(1, 1, 1, 0.7),
    xlims  = (-2.0, 2.0),
    ylims  = (-0.02, 0.85),
)

# Metropolis / smooth Metropolis (η-reg) — split at ±η to render the discontinuity faithfully.
i_neg_k, i_in_k, i_pos_k = split_at_eta(t_grid_plus, abs.(g_p_metro_s0))
plot!(plt_p2, t_grid_plus[i_neg_k], abs.(g_p_metro_s0[i_neg_k]);
      label = "\$|b_+^{(0,\\eta)}|\$ (Metropolis, η-reg)",
      color = c_kinky)
plot!(plt_p2, t_grid_plus[i_in_k],  abs.(g_p_metro_s0[i_in_k]);  label = "", color = c_kinky)
plot!(plt_p2, t_grid_plus[i_pos_k], abs.(g_p_metro_s0[i_pos_k]); label = "", color = c_kinky)

i_neg_s, i_in_s, i_pos_s = split_at_eta(t_grid_plus, abs.(g_p_smooth))
plot!(plt_p2, t_grid_plus[i_neg_s], abs.(g_p_smooth[i_neg_s]);
      label = "\$|b_+^{(s,\\eta)}|\$ (smooth Metropolis, η-reg, \$s=$(s_smooth)\$)",
      color = c_smooth)
plot!(plt_p2, t_grid_plus[i_in_s],  abs.(g_p_smooth[i_in_s]);  label = "", color = c_smooth)
plot!(plt_p2, t_grid_plus[i_pos_s], abs.(g_p_smooth[i_pos_s]); label = "", color = c_smooth)

# a-regularized smooth Metropolis — single continuous curve.
plot!(plt_p2, t_grid_plus, abs.(g_p_a_smooth);
      label = "\$|b_+^{(a,s)}|\$ (smooth Metropolis, a-reg, \$a=$(round(a_reg, digits=3))\$)",
      color = c_a_smooth)

# Mark the η boundaries (still meaningful for the η-reg curves)
vline!(plt_p2, [-η, η];
       color     = :grey40,
       linestyle = :dash,
       linewidth = 1.0,
       label     = "")
annotate!(plt_p2, η + 0.04, 0.81, text("\$|t|=\\eta\$", :grey40, 8, :left))

out_pdf_p2 = joinpath(@__DIR__, "..", "drafts", "plots", "b_plus_with_a_smooth.pdf")
out_png_p2 = joinpath(@__DIR__, "..", "drafts", "plots", "b_plus_with_a_smooth.png")
savefig(plt_p2, out_pdf_p2)
savefig(plt_p2, out_png_p2)
@printf("Saved %s\n", relpath(out_pdf_p2, joinpath(@__DIR__, "..")))
@printf("Saved %s\n", relpath(out_png_p2, joinpath(@__DIR__, "..")))

# ── Plot 2: b_-(t) ────────────────────────────────────────────────────────────
plt_m = plot(
    t_grid_minus, g_m_real;
    label  = "\$b_-(t)\$",
    color  = c_minus,
    xlabel = "\$t\$",
    ylabel = "\$b_-(t)\$",
    title  = "Outer kernel  \$b_-(t)\$  (\$\\bar{\\beta}=$(Int(β)),\\; \\sigma = 1/\\bar{\\beta}\$)",
    legend                  = :topright,
    foreground_color_legend = nothing,
    background_color_legend = RGBA(1, 1, 1, 0.7),
    xlims  = (-3.0, 3.0),
    ylims  = (-0.55, 0.55),
)
hline!(plt_m, [0.0]; color = :grey60, linestyle = :dot, linewidth = 1.0, label = "")

# L¹-norm side panel below the topright legend (single line).
add_norm_box!(plt_m,
    1.65, 2.95,
    0.30, 0.42,
    [
        "\$\\Vert b_-\\Vert_1 = $(@sprintf("%.3f", norm_b_minus))\$",
    ];
    fontsize = 9,
)

out_pdf_m = joinpath(@__DIR__, "..", "drafts", "plots", "b_minus.pdf")
out_png_m = joinpath(@__DIR__, "..", "drafts", "plots", "b_minus.png")
savefig(plt_m, out_pdf_m)
savefig(plt_m, out_png_m)
@printf("Saved %s\n", relpath(out_pdf_m, joinpath(@__DIR__, "..")))
@printf("Saved %s\n", relpath(out_png_m, joinpath(@__DIR__, "..")))
