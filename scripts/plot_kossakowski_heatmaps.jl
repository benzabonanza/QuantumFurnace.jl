#!/usr/bin/env julia
#
# Plot 3: Kossakowski matrix α_{ν1,ν2} heat maps for the disordered
# Heisenberg Hamiltonian.
#
# Strategy: build a 1×3 panel that pedagogically isolates the
# off-diagonal contribution Chen et al. add over Davies.
#
#   Panel 1 — Davies (GNS DB) at β = 10:
#             α^GNS_{ν,ν} on the diagonal, zero elsewhere. Davies only
#             admits GNS detailed balance; per Corollary cor:kossakowski-shift
#             of 2_methods.tex (line 413), α^GNS_{ν1,ν2} = α^KMS_{ν1−βσ²/2,
#             ν2−βσ²/2}, so the diagonal entries are evaluated as
#             create_alpha(ν − βσ²/2, ν − βσ²/2, β, σ, 0, s).
#             The diagonal-only structure does not change qualitatively
#             with β, so a single panel suffices.
#   Panel 2 — KMS Chen at β = 10: full α^KMS_{ν1,ν2} = create_alpha(...).
#   Panel 3 — KMS Chen at β = 20: full α^KMS_{ν1,ν2} (colder).
#
# All three panels use the smooth Metropolis filter (a = 0, s = 0.25),
# matching scripts/plot_transition_weights.jl and plot_b_plus_b_minus.jl.
# Colour scale is shared across all three panels so the visual gap
# between Davies (just diagonal) and KMS Chen (diagonal + skirt) is
# directly the off-diagonal contribution.
#
# System: disordered Heisenberg, n = 3, periodic, ferromagnetic XXZ
# (coeffs = [-1, -1, -1]) with on-site Z disorder. n = 3 (odd) is not
# bipartite so Z-only disorder gives a unique spectrum with no
# Bohr-frequency collisions (memory: bohr_collision_sectors).
#
# PHYSICS CHECK: σ = 1/β follows the Chen et al. canonical setting and
# matches plots 1 and 2 in this epic. With σ = 1/β the KMS DB
# Kossakowski matrix is skew-symmetric (eq. line 160 of 2_methods.tex).
#
# PHYSICS CHECK: s = 0.25 matches plots 1, 2 of this epic.
#
# PHYSICS CHECK: β values 10, 20 (warm vs cold) — Bohr spectrum spans
# O(1), σ = 1/β. The Gaussian factor exp(-(ν1-ν2)²/(8σ²)) controls
# off-diagonal decay; at β = 20 the decay length is ν ≈ 0.14, at
# β = 10 it is ν ≈ 0.28, so the cold panel should show a tighter
# off-diagonal skirt than the warm one.
#
# Output: drafts/plots/kossakowski_heatmap.{pdf,png}
#
# Colour map: "ocean inverted" gradient (deep navy → mint), per
# memory: reference_gradient_palettes. Same dark-base / bright-peak feel as
# `:inferno` but anchored at navy (#003147) instead of pure black, which reads
# more elegantly against the white thesis page.
#
# Usage:  julia --project scripts/plot_kossakowski_heatmaps.jl

using Printf
using Random
using LinearAlgebra
using Plots
using QuantumFurnace

# ── Parameters ────────────────────────────────────────────────────────────────
const num_qubits = 3
const coeffs     = [-1.0, -1.0, -1.0]
const a_reg      = 0.0
const s_smooth   = 0.25
const β_warm     = 10.0
const β_cold     = 20.0
const seed       = 42
const batch_size = 200

# ── Build disordered Heisenberg ──────────────────────────────────────────────
Random.seed!(seed)
@info "Building disordered Heisenberg (n=$num_qubits, batch=$batch_size, periodic, Z disorder)…"
raw = find_ideal_heisenberg(num_qubits, coeffs; batch_size=batch_size, periodic=true)

# β passed to HamHam only initialises the (unused) Gibbs state; α uses β explicitly below.
ham = HamHam(raw, β_warm)
unique_freqs = sort(collect(keys(ham.bohr_dict)))
Nfreq = length(unique_freqs)
@printf("Unique Bohr frequencies: %d   (expected %d for n=%d)\n",
        Nfreq, 2^num_qubits * (2^num_qubits - 1) + 1, num_qubits)
@printf("ν_min (smallest positive Bohr): %.6f\n", raw.nu_min)
@printf("ν range: [%.4f, %.4f]\n", first(unique_freqs), last(unique_freqs))

# ── Compute Kossakowski matrices α_{ν1,ν2} ───────────────────────────────────
function compute_alpha_matrix(unique_freqs::Vector{Float64}, β::Real;
                              σ::Real=1/β, a::Real=a_reg, s::Real=s_smooth)
    N = length(unique_freqs)
    A = zeros(Float64, N, N)
    for j in 1:N, i in 1:N
        A[i, j] = create_alpha(unique_freqs[i], unique_freqs[j], β, σ, a, s)
    end
    return A
end

@info "Computing α at β=$(β_warm), σ=$(round(1/β_warm, digits=3))…"
α_warm = compute_alpha_matrix(unique_freqs, β_warm)
@info "Computing α at β=$(β_cold), σ=$(round(1/β_cold, digits=3))…"
α_cold = compute_alpha_matrix(unique_freqs, β_cold)

# Davies (GNS DB) panel: diagonal of α^GNS at β_warm, off-diagonals zeroed.
# Per Corollary cor:kossakowski-shift (2_methods.tex line 413):
#   α^GNS_{ν1,ν2} = α^KMS_{ν1 − βσ²/2, ν2 − βσ²/2},
# so on the diagonal α^GNS(ν, ν) = α^KMS(ν − βσ²/2, ν − βσ²/2).
# At σ = 1/β the shift βσ²/2 = 1/(2β) is small compared to the ν range
# (≈ 0.05 at β = 10), but it correctly identifies the Davies generator
# with GNS DB — the only kind of detailed balance Davies satisfies.
σ_warm = 1 / β_warm
shift_warm = β_warm * σ_warm^2 / 2   # = 1/(2β) for σ = 1/β
α_davies = zeros(size(α_warm))
for i in 1:Nfreq
    ν = unique_freqs[i]
    α_davies[i, i] = create_alpha(ν - shift_warm, ν - shift_warm,
                                   β_warm, σ_warm, a_reg, s_smooth)
end
@printf("Davies shift βσ²/2 at β=%d: %.4f   (ν range [%.2f, %.2f])\n",
        Int(β_warm), shift_warm, first(unique_freqs), last(unique_freqs))

# Sanity: skew symmetry of α at σ = 1/β (Eq. line 160 of 2_methods.tex):
#   α(ν1, ν2) = α(-ν2, -ν1) · exp(-β (ν1 + ν2) / 2)
function skew_residual(α::Matrix, freqs::Vector, β::Real)
    res = 0.0
    for j in 1:length(freqs), i in 1:length(freqs)
        ν1, ν2 = freqs[i], freqs[j]
        # find indices of (-ν2, -ν1) — since freqs is sorted symmetric set, mirror via reverse
        i_mirror = findfirst(≈(-ν2), freqs)
        j_mirror = findfirst(≈(-ν1), freqs)
        i_mirror === nothing && continue
        j_mirror === nothing && continue
        rhs = α[i_mirror, j_mirror] * exp(-β * (ν1 + ν2) / 2)
        res = max(res, abs(α[i, j] - rhs))
    end
    return res
end
@printf("Skew-symmetry residual α_warm: %.2e   (should be ≪ 1)\n",
        skew_residual(α_warm, unique_freqs, β_warm))
@printf("Skew-symmetry residual α_cold: %.2e   (should be ≪ 1)\n",
        skew_residual(α_cold, unique_freqs, β_cold))

# Shared colour scale across all three panels
α_max = max(maximum(α_warm), maximum(α_cold))
α_min = 0.0
@printf("α range across panels: [%.3e, %.3e]\n", α_min, α_max)

# Off-diagonal magnitude diagnostics
function offdiag_max(α::Matrix)
    m = 0.0
    for j in 1:size(α, 2), i in 1:size(α, 1)
        i == j && continue
        m = max(m, α[i, j])
    end
    return m
end
@printf("max diag(α_warm) = %.3e   max off-diag(α_warm) = %.3e   ratio = %.2f\n",
        maximum(diag(α_warm)), offdiag_max(α_warm),
        maximum(diag(α_warm)) / offdiag_max(α_warm))
@printf("max diag(α_cold) = %.3e   max off-diag(α_cold) = %.3e   ratio = %.2f\n",
        maximum(diag(α_cold)), offdiag_max(α_cold),
        maximum(diag(α_cold)) / offdiag_max(α_cold))

# ── Plotting defaults (match other thesis-plot scripts) ──────────────────────
default(
    fontfamily              = "Computer Modern",
    titlefontsize           = 11,
    guidefontsize           = 11,
    tickfontsize            = 8,
    legendfontsize          = 8,
    framestyle              = :axes,
    grid                    = false,
    margin                  = 3Plots.mm,
    # Hide axis lines and tick marks while keeping tick *labels* visible:
    # ν values render in the default text colour, but the L-shape of
    # bottom + left axis lines is painted white (i.e. invisible on a
    # white canvas).
    foreground_color_axis   = :white,
    foreground_color_border = :white,
    tick_direction          = :none,
)

# Index-based axes with ~7 sparse ν labels — uniform cells are easier to read
# than the non-uniformly-spaced raw Bohr frequencies.
tick_idx = round.(Int, range(1, Nfreq; length=7))
tick_lbl = [@sprintf("%.2f", unique_freqs[i]) for i in tick_idx]

# ── Colour map ────────────────────────────────────────────────────────────────
# "Ocean inverted": deep navy (#003147) zero-anchor → mint (#B7E6A5) peaks.
# Same dark-base / bright-peak feel as `:inferno`, but with a non-black base
# so the figure stays elegant on white paper.
const CMAP = cgrad(["#003147", "#045275", "#00718B", "#089099",
                    "#46AEA0", "#7CCBA2", "#B7E6A5"])

# All three panels: no individual colorbar, identical bottom margin so they
# render at the same physical size; ν₂ is labelled only on the middle panel
# (visually centred, serves as the shared x-axis label for all three).
function make_heatmap(α::Matrix, title::String;
                      show_y::Bool=true, show_xlabel::Bool=false,
                      left_margin=3Plots.mm)
    return heatmap(
        α;
        c             = CMAP,
        clims         = (α_min, α_max),
        xticks        = (tick_idx, tick_lbl),
        yticks        = show_y ? (tick_idx, tick_lbl) : nothing,
        xlabel        = show_xlabel ? "\$\\nu_2\$" : "",
        ylabel        = show_y ? "\$\\nu_1\$" : "",
        title         = title,
        colorbar      = false,
        aspect_ratio  = :equal,
        xrotation     = 45,
        framestyle    = :axes,
        left_margin   = left_margin,
        bottom_margin = 8Plots.mm,
        top_margin    = 3Plots.mm,
    )
end

p1 = make_heatmap(α_davies, "Davies (GNS DB), \$\\beta=$(Int(β_warm))\$";
                  show_y=true,  show_xlabel=false, left_margin=10Plots.mm)
p2 = make_heatmap(α_warm,   "CKG (KMS DB), \$\\beta=$(Int(β_warm))\$";
                  show_y=false, show_xlabel=true)
p3 = make_heatmap(α_cold,   "CKG (KMS DB), \$\\beta=$(Int(β_cold))\$";
                  show_y=false, show_xlabel=false)

# Shared vertical colorbar — thin gradient column as a 4th subplot.
n_cb     = 256
cb_y     = collect(range(α_min, α_max; length=n_cb))
cb_data  = repeat(reshape(cb_y, n_cb, 1), 1, 2)
cb_ticks = round.(range(α_min, α_max; length=5); digits=2)
cbar_panel = heatmap(
    [0.0, 1.0], cb_y, cb_data;
    c              = CMAP,
    clims          = (α_min, α_max),
    xticks         = nothing,
    yticks         = cb_ticks,
    ymirror        = true,
    xlabel         = "\$\\alpha_{\\nu_1,\\nu_2}\$",
    title          = "",
    colorbar       = false,
    framestyle     = :axes,
    grid           = false,
    aspect_ratio   = :auto,
    top_margin     = 4Plots.mm,
    bottom_margin  = 8Plots.mm,
    left_margin    = 0Plots.mm,
    right_margin   = 0Plots.mm,
)

l = @layout [grid(1,3) b{0.025w}]
plt_all = plot(p1, p2, p3, cbar_panel;
               layout = l,
               size   = (1300, 460),
               dpi    = 200,
               plot_title = "Kossakowski matrix \$\\alpha_{\\nu_1,\\nu_2}\$ (smooth Metropolis)",
               plot_titlefontsize = 15,
               plot_titlevspan   = 0.14)

out_pdf = joinpath(@__DIR__, "..", "drafts", "plots", "kossakowski_heatmap.pdf")
out_png = joinpath(@__DIR__, "..", "drafts", "plots", "kossakowski_heatmap.png")
mkpath(dirname(out_pdf))
savefig(plt_all, out_pdf)
savefig(plt_all, out_png)
@printf("Saved %s\n", relpath(out_pdf, joinpath(@__DIR__, "..")))
@printf("Saved %s\n", relpath(out_png, joinpath(@__DIR__, "..")))
