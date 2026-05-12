#!/usr/bin/env julia
#
# DLL Gaussian vs DLL Metropolis Kossakowski matrix α_{ν1,ν2} heat maps,
# at β ∈ {5, 10, 20}. Companion to scripts/plot_kossakowski_heatmaps.jl
# (which does the CKG/Davies side); same colour map, axes, and layout
# style for direct visual comparison in the thesis.
#
# Strategy: 2×3 panel grid + shared colorbar.
#   row 1: DLL Gaussian (Eq. 3.21–3.22) — narrow lobe centred at -1/β,
#          width 2/β. Lobe shrinks ∝ 1/β as β grows.
#   row 2: DLL Metropolis (Eq. 3.19–3.20, S=2) — fully-saturated (-,-)
#          quadrant, exponential decay into (+,+). Lobe SUPPORT stays
#          O(1) as β grows; the (+,+) corner shrinks (mandatory KMS-DB).
#   columns: β ∈ {5, 10, 20}.
#
# Both DLL Kossakowskis are rank-1 outer products α = v·vᵀ with
# v_k = freq_kernel(filter, ν_k); both non-negative (freq_kernel real).
# Plotted on a shared linear scale [0, max(α)] so panel sizes match
# the existing CKG figure exactly.
#
# System: same n=3 disordered Heisenberg fixture as
# `plot_kossakowski_heatmaps.jl` (apples-to-apples comparison with the
# CKG figure).
#
# PHYSICS CHECK: S = 2 default for the Metropolis filter — flat-top
# region [-1, 1] strictly contains the n=3 fixture's Bohr range
# [-0.45, 0.45], so the bump is invisible to the matrix on this grid.
# At β=20 the (-,-) quadrant saturates near the analytic value f̂(ν → -∞)
# = 1, while the (+,+) quadrant decays as exp(-β(ν+ν')) — KMS-DB forced.
#
# Output: drafts/plots/kossakowski_dll_heatmaps.{pdf,png}
# Usage:  julia --project scripts/plot_kossakowski_dll_heatmaps.jl

using Printf
using LinearAlgebra
using Plots
using QuantumFurnace

# Reuse the n=3 fixture loader from the test helpers (same fixture as
# every other DLL test/diagnostic in this repo).
include(joinpath(@__DIR__, "..", "test", "test_helpers.jl"))

# ── Parameters ────────────────────────────────────────────────────────────────
const β_values = (5.0, 10.0, 20.0)
const S_meta   = 2.0

# ── Build n=3 fixture ────────────────────────────────────────────────────────
ham_path = joinpath(@__DIR__, "..", "hamiltonians", "heis_xxx_zzdisordered_periodic_n3.bson")
ham = _load_test_hamiltonian(ham_path, first(β_values))
unique_freqs = sort(collect(keys(ham.bohr_dict)))
Nfreq = length(unique_freqs)
@printf("Unique Bohr frequencies: %d   ν range: [%.4f, %.4f]\n",
        Nfreq, first(unique_freqs), last(unique_freqs))

# ── Compute Kossakowski matrices (rank-1 outer products) ─────────────────────
function compute_alpha(filter::AbstractFilter, freqs::Vector{Float64})
    v = [Float64(freq_kernel(filter, ν)) for ν in freqs]
    return v * transpose(v)  # real, non-negative, rank-1
end

println("\nComputing α matrices…")
α_panels = Dict{Tuple{Symbol,Float64},Matrix{Float64}}()
for β in β_values
    α_panels[(:gauss, β)] = compute_alpha(DLLGaussianFilter(β), unique_freqs)
    α_panels[(:meta,  β)] = compute_alpha(DLLMetropolisFilter(β; S=S_meta), unique_freqs)
end

# Diagnostic summary.
@printf("\n%-6s %-10s %-12s %-12s %-12s\n", "β", "filter", "max(α)", "‖α‖_F", "Σ|α|")
println("-"^60)
for β in β_values
    for (sym, name) in ((:gauss, "Gaussian"), (:meta, "Metropolis"))
        α = α_panels[(sym, β)]
        @printf("%-6.1f %-10s %-12.4e %-12.4e %-12.4e\n",
                β, name, maximum(α), norm(α), sum(abs.(α)))
    end
end

# Shared colour scale across all 6 panels — matches the CKG-figure
# convention so the visual story is "same axes, different filter shape".
α_max = maximum(maximum(α) for α in values(α_panels))
α_min = 0.0
@printf("\nShared α range: [%.3e, %.3e]\n", α_min, α_max)

# ── Plotting defaults — match scripts/plot_kossakowski_heatmaps.jl ───────────
default(
    fontfamily              = "Computer Modern",
    titlefontsize           = 11,
    guidefontsize           = 11,
    tickfontsize            = 8,
    legendfontsize          = 8,
    framestyle              = :axes,
    grid                    = false,
    margin                  = 3Plots.mm,
    foreground_color_axis   = :white,
    foreground_color_border = :white,
    tick_direction          = :none,
)

# Sparse tick labels — 7 evenly-spaced ν values.
tick_idx = round.(Int, range(1, Nfreq; length=7))
tick_lbl = [@sprintf("%.2f", unique_freqs[i]) for i in tick_idx]

# Same "ocean inverted" colormap as the CKG figure (deep navy → mint).
const CMAP = cgrad(["#003147", "#045275", "#00718B", "#089099",
                    "#46AEA0", "#7CCBA2", "#B7E6A5"])

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

# Build the 6 heatmaps. Top row: Gaussian. Bottom row: Metropolis.
# Y labels on column 1 only; x labels on bottom-row middle column.
panels = Plots.Plot[]
for (i, β) in enumerate(β_values)
    push!(panels, make_heatmap(
        α_panels[(:gauss, β)],
        "DLL Gaussian, \$\\beta=$(Int(β))\$";
        show_y      = (i == 1),
        show_xlabel = false,
        left_margin = i == 1 ? 10Plots.mm : 3Plots.mm,
    ))
end
for (i, β) in enumerate(β_values)
    push!(panels, make_heatmap(
        α_panels[(:meta, β)],
        "DLL Metropolis (\$S=$(Int(S_meta))\$), \$\\beta=$(Int(β))\$";
        show_y      = (i == 1),
        show_xlabel = (i == 2),  # middle column carries the shared x-label
        left_margin = i == 1 ? 10Plots.mm : 3Plots.mm,
    ))
end

# Shared vertical colorbar — thin gradient column as a 7th subplot.
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

# Layout: 6 heatmap panels in a 2×3 grid + one thin colorbar column.
l = @layout [grid(2,3) b{0.025w}]
plt_all = plot(panels..., cbar_panel;
               layout = l,
               size   = (1300, 880),
               dpi    = 200,
               plot_title = "DLL Kossakowski matrix \$\\alpha_{\\nu_1,\\nu_2}\$ — Gaussian vs Metropolis",
               plot_titlefontsize = 15,
               plot_titlevspan    = 0.07)

out_pdf = joinpath(@__DIR__, "..", "drafts", "plots", "kossakowski_dll_heatmaps.pdf")
out_png = joinpath(@__DIR__, "..", "drafts", "plots", "kossakowski_dll_heatmaps.png")
mkpath(dirname(out_pdf))
savefig(plt_all, out_pdf)
savefig(plt_all, out_png)
@printf("Saved %s\n", relpath(out_pdf, joinpath(@__DIR__, "..")))
@printf("Saved %s\n", relpath(out_png, joinpath(@__DIR__, "..")))
