#!/usr/bin/env julia
#
# Plot: ‖b_+^{(s,η)}‖_1 as a function of η for kinky (s=0) and smooth (s=0.25)
# Metropolis. Marks the threshold η_crit where ‖b_+‖_1 = 1, the practical limit
# beyond which the LCU block encoding can no longer be made deterministic.
#
# Sources:
#   • Thesis 2_methods.tex eq:b_plus-s-eta (line 524) — kernel definition
#   • Thesis 2_methods.tex line 528          — analytic upper bound,
#       ‖b_+^{(s,η)}‖_1 < 2/(5π^{3/2}) − ln(1+s)/(2√2π²) + ln(1/η)/(√2π²),
#     which corrects Chen et al.'s prefactor by 1/π^{3/2} (line 532).
#
# Strategy: integrate |b_+^{(s,η)}(t)| via quadgk on a wide window for a sweep
# of η on a log grid. Bisect to locate the crossing ‖b_+‖_1 = 1 for each s.
# Overlay the thesis upper bound to make the looseness explicit (factor ~7
# between numerical norm and bound at L¹ = 1).
#
# PHYSICS CHECK: σβ = 1 (canonical Chen choice σ = 1/β). The kernel shape is
# β-independent in this scaling, so a single sweep characterises the kinky /
# smooth Metropolis kernels universally.
#
# PHYSICS CHECK: s = 0.25 matches the choice in plot_b_plus_b_minus.jl and
# plot_transition_weights.jl for visual consistency across the thesis figures.
#
# Output: drafts/plots/b_plus_l1_vs_eta.{pdf,png}
#
# Usage: julia --project scripts/plot_b_plus_l1_vs_eta.jl

using Printf
using QuadGK: quadgk
using Plots

# ── Parameters (canonical Chen choice σβ = 1) ────────────────────────────────
const β        = 10.0
const σ        = 1 / β
const s_smooth = 0.25

# ── b_+^{(s,η)}(t) — mirror of thesis eq:b_plus-s-eta ─────────────────────────
function b_plus_metro(t::Real, η::Real; s::Real=0.0)
    if abs(t) < 1e-14
        return complex((2 - σ^2 * β^2 * (1 + s)) / (2 * sqrt(2) * pi^2))
    elseif abs(t) <= η
        numerator = exp(-σ^2 * β^2 * (2 * t^2 + 1im * t) * (1 + s)) +
                    1im * (2 * t + 1im)
    else
        numerator = exp(-σ^2 * β^2 * (2 * t^2 + 1im * t) * (1 + s))
    end
    return (1 / (2 * sqrt(2) * pi^2)) * numerator / (t * (2 * t + 1im))
end

l1_norm(η::Real, s::Real) =
    quadgk(t -> abs(b_plus_metro(t, η; s=s)), -50.0, 50.0; rtol=1e-11)[1]

# Thesis upper bound (line 528): increases as log(1/η) with η ↘ 0.
thesis_bound(η::Real, s::Real) = 2 / (5 * pi^1.5) -
                                  log(1 + s) / (2 * sqrt(2) * pi^2) +
                                  log(1 / η) / (sqrt(2) * pi^2)

# Bisect log10(η) for ‖b_+^{(s,η)}‖_1 = 1.
function η_crit_numeric(s::Real; lo::Real=-12.0, hi::Real=0.0, tol::Real=1e-7)
    while hi - lo > tol
        m = (lo + hi) / 2
        (l1_norm(10.0^m, s) > 1.0) ? (lo = m) : (hi = m)
    end
    return 10.0^hi
end

# Closed-form thesis-bound η_crit: solve bound(η, s) = 1 for log(1/η).
η_crit_bound(s::Real) = exp(-(1 - 2 / (5 * pi^1.5) +
                              log(1 + s) / (2 * sqrt(2) * pi^2)) * sqrt(2) * pi^2)

# ── Compute crossings ─────────────────────────────────────────────────────────
@info "Locating ‖b_+‖_1 = 1 crossings via bisection…"
η_crit_kinky_num   = η_crit_numeric(0.0)
η_crit_smooth_num  = η_crit_numeric(s_smooth)
η_crit_kinky_bnd   = η_crit_bound(0.0)
η_crit_smooth_bnd  = η_crit_bound(s_smooth)

@printf("Kinky  Metropolis (s=0):     η_crit_num   = %.3e   ‖b_+‖_1 = %.4f\n",
        η_crit_kinky_num, l1_norm(η_crit_kinky_num, 0.0))
@printf("Kinky  Metropolis (s=0):     η_crit_bound = %.3e\n", η_crit_kinky_bnd)
@printf("Smooth Metropolis (s=0.25):  η_crit_num   = %.3e   ‖b_+‖_1 = %.4f\n",
        η_crit_smooth_num, l1_norm(η_crit_smooth_num, s_smooth))
@printf("Smooth Metropolis (s=0.25):  η_crit_bound = %.3e\n", η_crit_smooth_bnd)
@printf("Bound looseness at L¹=1: kinky  factor %.2f, smooth factor %.2f\n",
        η_crit_kinky_bnd  / η_crit_kinky_num,
        η_crit_smooth_bnd / η_crit_smooth_num)

# ── Sweep η on a log grid ─────────────────────────────────────────────────────
log10_η_grid = collect(range(0.0, -10.0; length=121))   # η from 1 down to 1e-10
η_grid       = 10.0 .^ log10_η_grid

@info "Computing L¹ norm on a $(length(η_grid))-point log grid (kinky + smooth)…"
norm_kinky  = [l1_norm(η, 0.0)      for η in η_grid]
norm_smooth = [l1_norm(η, s_smooth) for η in η_grid]
bnd_kinky   = [thesis_bound(η, 0.0)      for η in η_grid]
bnd_smooth  = [thesis_bound(η, s_smooth) for η in η_grid]

# ── Plot defaults (match plot_b_plus_b_minus.jl) ──────────────────────────────
default(
    fontfamily       = "Computer Modern",
    titlefontsize    = 12,
    guidefontsize    = 12,
    tickfontsize     = 10,
    legendfontsize   = 9,
    framestyle       = :box,
    grid             = true,
    gridalpha        = 0.25,
    linewidth        = 2.0,
    size             = (720, 480),
    dpi              = 200,
    margin           = 4Plots.mm,
)

# Thesis colour palette (matches plot_b_plus_b_minus.jl).
c_kinky  = "#7A2E39"   # bordeaux
c_smooth = "#735874"   # deepplum
c_bound  = "#C29A3A"   # ochre — analytic upper bound
c_thresh = "#2D5A3D"   # pinegreen — L¹ = 1 threshold

# ── Figure ────────────────────────────────────────────────────────────────────
plt = plot(;
    xlabel = "\$-\\log_{10}\\eta\$",
    ylabel = "\$\\Vert b_+^{(s,\\eta)}\\Vert_1\$",
    title  = "\$\\ell_1\$ norm of the regularized Metropolis kernel  (\$\\sigma = 1/\\beta\$)",
    legend                  = :topleft,
    foreground_color_legend = nothing,
    background_color_legend = RGBA(1, 1, 1, 0.85),
    xlims  = (0.0, 10.0),
    ylims  = (0.0, 1.7),
)

# Numerical curve (smooth Metropolis).
plot!(plt, -log10_η_grid, norm_smooth;
      label = "\$\\Vert b_+^{(s,\\eta)}\\Vert_1\$  (smooth, \$s=$(s_smooth)\$, numeric)",
      color = c_smooth)

# Upper bound.
plot!(plt, -log10_η_grid, bnd_smooth;
      label     = "upper bound (\$s=$(s_smooth)\$)",
      color     = c_bound,
      linestyle = :dash,
      linewidth = 1.6)

# Block-encoding threshold ‖b_+‖_1 = 1.
hline!(plt, [1.0]; color = c_thresh, linestyle = :dashdot, linewidth = 1.4,
       label = "\$\\Vert b_+\\Vert_1 = 1\$")

# Mark the numerical crossing.
scatter!(plt, [-log10(η_crit_smooth_num)], [1.0];
         color = c_smooth, markershape = :diamond, markersize = 5, label = "")

# Annotate crossing position to the right of the marker, at the same height.
annotate!(plt, -log10(η_crit_smooth_num) + 0.20, 0.93,
          text(@sprintf("\$\\eta_{\\rm crit}^{(s)} = %.2f\\!\\times\\!10^{-7}\$",
                        η_crit_smooth_num * 1e7),
               c_smooth, 9, :left))

# ── Save ──────────────────────────────────────────────────────────────────────
out_pdf = joinpath(@__DIR__, "..", "drafts", "plots", "b_plus_l1_vs_eta.pdf")
out_png = joinpath(@__DIR__, "..", "drafts", "plots", "b_plus_l1_vs_eta.png")
mkpath(dirname(out_pdf))
savefig(plt, out_pdf)
savefig(plt, out_png)
@printf("Saved %s\n", relpath(out_pdf, joinpath(@__DIR__, "..")))
@printf("Saved %s\n", relpath(out_png, joinpath(@__DIR__, "..")))
