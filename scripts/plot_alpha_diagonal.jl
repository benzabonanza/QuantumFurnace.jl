#!/usr/bin/env julia
#
# Plot α(ν, ν) diagonal-of-Kossakowski curves for the three samplers
# (CKG smooth-Metro, DLL Metropolis, DLL Gaussian) across (n, β) cells.
#
# Companion to scripts/scratch_alpha_diagonal_analysis.jl (qf-mto.5).
# Reads scripts/output/alpha_diagonal_analysis.bson and produces
#   drafts/figures/numerics/alpha_diagonal_three_samplers.{png,pdf}
#
# Layout: 4 (β) × 3 (n) grid of panels.
#   Rows: β ∈ {1, 5, 10, 20} (top → bottom)
#   Cols: n ∈ {3, 4, 5}      (left → right)
#   Curves per panel: one per sampler, log-scale on y to show the
#   tails (CKG and DLL Metro share the Metropolis tail, DLL Gauss
#   collapses super-exponentially at large |βν|).
#
# Companion β-distance plot: distance vs β at fixed n (one row, three
# panels for n ∈ {3,4,5}; two lines per panel: DLL Metro and DLL Gauss
# baselined to CKG smooth-Metro). Saved as
#   drafts/figures/numerics/alpha_diagonal_distance_vs_beta.{png,pdf}
#
# Run: julia --project scripts/plot_alpha_diagonal.jl

using Printf
using BSON
using Plots

# ---------------------------------------------------------------------------
# Load BSON
# ---------------------------------------------------------------------------

const ANALYSIS_BSON = joinpath(@__DIR__, "output", "alpha_diagonal_analysis.bson")
@assert isfile(ANALYSIS_BSON) "Run scripts/scratch_alpha_diagonal_analysis.jl first ($(ANALYSIS_BSON) missing)."

raw = BSON.load(ANALYSIS_BSON)
const records = raw[:records]
const NS      = raw[:ns]
const BETAS   = raw[:betas]
const SAMPLERS = raw[:samplers]

@printf("Loaded %d records: %d (n, β) cells × %d samplers = %d records\n",
        length(records), length(NS) * length(BETAS), length(SAMPLERS),
        length(NS) * length(BETAS) * length(SAMPLERS))

# ---------------------------------------------------------------------------
# Thesis colour palette (matches numerics_ckg_vs_dll_taumix_comparison.jl
# and plot_kossakowski_dll_heatmaps.jl families)
# ---------------------------------------------------------------------------

const COLOR_CKG       = "#003147"  # deep navy (= CKG family root in heatmap palette)
const COLOR_DLL_METRO = "#089099"  # teal-cyan (DLL Metro)
const COLOR_DLL_GAUSS = "#7CCBA2"  # mint-green (DLL Gauss)

const SAMPLER_COLOR = Dict(
    "CKG smooth-Metro" => COLOR_CKG,
    "DLL Metropolis"   => COLOR_DLL_METRO,
    "DLL Gaussian"     => COLOR_DLL_GAUSS,
)

const SAMPLER_STYLE = Dict(
    "CKG smooth-Metro" => :solid,
    "DLL Metropolis"   => :dash,
    "DLL Gaussian"     => :dashdot,
)

const SAMPLER_MARKER = Dict(
    "CKG smooth-Metro" => :circle,
    "DLL Metropolis"   => :diamond,
    "DLL Gaussian"     => :star5,
)

# Floor for log-scale plots (clamp tiny α_diag values to avoid log(0)).
const Y_FLOOR = 1e-10

# ---------------------------------------------------------------------------
# Plotting defaults — match thesis-numerics conventions
# ---------------------------------------------------------------------------

default(
    fontfamily      = "Computer Modern",
    titlefontsize   = 11,
    guidefontsize   = 10,
    tickfontsize    = 7,
    legendfontsize  = 8,
    framestyle      = :box,
    grid            = true,
    gridlinewidth   = 0.4,
    gridstyle       = :dot,
    gridalpha       = 0.4,
    margin          = 2Plots.mm,
)

# ---------------------------------------------------------------------------
# Helper: fetch the three samplers for one (n, β) cell
# ---------------------------------------------------------------------------

function records_for(n::Int, beta::Real)
    matches = filter(r -> r.n == n && r.beta == beta, records)
    @assert length(matches) == length(SAMPLERS) "Expected $(length(SAMPLERS)) samplers, got $(length(matches)) at (n=$n, β=$beta)"
    return matches
end

# ---------------------------------------------------------------------------
# Build per-panel α_diag(ν) plots
# ---------------------------------------------------------------------------

function panel_alpha_diag(n::Int, beta::Real;
                          show_y::Bool, show_x::Bool, legend_in_panel::Bool)
    matches = records_for(n, beta)

    # Compute global y range across the three samplers in this panel
    all_y = Float64[]
    for r in matches
        append!(all_y, r.α_diag)
    end
    y_min = max(Y_FLOOR, minimum(filter(>(0), all_y)) / 2)
    y_max = maximum(all_y) * 2

    plt = plot(
        xlabel = show_x ? "\$\\nu\$" : "",
        ylabel = show_y ? "\$\\alpha_{\\mathrm{diag}}(\\nu)\$" : "",
        title  = "\$n=$n,\\ \\beta=$(Int(beta))\$",
        yscale = :log10,
        ylims  = (y_min, y_max),
        legend = legend_in_panel ? :bottomleft : false,
        legendfontsize = 7,
        size   = (320, 240),
    )

    for sampler in SAMPLERS
        r = matches[findfirst(rec -> rec.sampler == sampler, matches)]
        ν = r.ν_grid
        α = max.(r.α_diag, Y_FLOOR)  # clamp tiny values for log scale
        plot!(plt, ν, α;
            label     = sampler,
            color     = SAMPLER_COLOR[sampler],
            linestyle = SAMPLER_STYLE[sampler],
            linewidth = 1.6,
            marker    = SAMPLER_MARKER[sampler],
            markersize = 1.6,
            markerstrokewidth = 0,
        )
    end
    return plt
end

# ---------------------------------------------------------------------------
# Compose the 4×3 grid (β rows × n cols)
# ---------------------------------------------------------------------------

println("\nBuilding 4 (β) × 3 (n) panel grid…")

panels = Plots.Plot[]
for (i, β) in enumerate(BETAS)
    for (j, n) in enumerate(NS)
        legend_in_panel = (i == 1 && j == length(NS))  # legend in top-right cell only
        push!(panels, panel_alpha_diag(n, β;
            show_y = (j == 1),
            show_x = (i == length(BETAS)),
            legend_in_panel = legend_in_panel))
    end
end

plt_grid = plot(panels...;
    layout = (length(BETAS), length(NS)),
    size   = (1500, 1100),
    dpi    = 200,
    plot_title = "Kossakowski diagonal \$\\alpha_{\\mathrm{diag}}(\\nu) = \\sum_a \\alpha^a(\\nu, \\nu)\$ — three samplers",
    plot_titlefontsize = 14,
    plot_titlevspan    = 0.04,
)

# ---------------------------------------------------------------------------
# Save the main figure
# ---------------------------------------------------------------------------

const OUT_DIR = joinpath(@__DIR__, "..", "drafts", "figures", "numerics")
mkpath(OUT_DIR)
const OUT_PNG = joinpath(OUT_DIR, "alpha_diagonal_three_samplers.png")
const OUT_PDF = joinpath(OUT_DIR, "alpha_diagonal_three_samplers.pdf")

savefig(plt_grid, OUT_PNG)
savefig(plt_grid, OUT_PDF)
@printf("Saved %s\n", relpath(OUT_PNG, dirname(@__DIR__)))
@printf("Saved %s\n", relpath(OUT_PDF, dirname(@__DIR__)))

# ---------------------------------------------------------------------------
# Companion: distance vs β at fixed n (one panel per n)
# ---------------------------------------------------------------------------

println("\nBuilding distance-vs-β companion plot…")

dist_panels = Plots.Plot[]
for (j, n) in enumerate(NS)
    βs = sort(BETAS)
    dist_metro = Float64[]
    dist_gauss = Float64[]
    for β in βs
        matches = records_for(n, β)
        r_metro = matches[findfirst(r -> r.sampler == "DLL Metropolis", matches)]
        r_gauss = matches[findfirst(r -> r.sampler == "DLL Gaussian", matches)]
        push!(dist_metro, r_metro.dist_to_ckg)
        push!(dist_gauss, r_gauss.dist_to_ckg)
    end

    plt = plot(
        xlabel = "\$\\beta\$",
        ylabel = (j == 1) ? "\$\\|\\alpha_{\\mathrm{diag}}^{\\mathrm{DLL}} - \\alpha_{\\mathrm{diag}}^{\\mathrm{CKG}}\\|_2 \\,/\\, \\|\\alpha_{\\mathrm{diag}}^{\\mathrm{CKG}}\\|_2\$" : "",
        title  = "\$n=$n\$",
        ylims  = (0.0, max(maximum(dist_gauss), maximum(dist_metro)) * 1.15),
        legend = (j == length(NS)) ? :topleft : false,
        size   = (430, 320),
    )
    plot!(plt, βs, dist_metro;
        label     = "DLL Metropolis vs CKG sM",
        color     = COLOR_DLL_METRO,
        linestyle = :dash,
        linewidth = 1.8,
        marker    = :diamond,
        markersize = 5,
    )
    plot!(plt, βs, dist_gauss;
        label     = "DLL Gaussian vs CKG sM",
        color     = COLOR_DLL_GAUSS,
        linestyle = :dashdot,
        linewidth = 1.8,
        marker    = :star5,
        markersize = 6,
    )
    # Threshold guide
    hline!(plt, [0.30]; color = "#cc6677", linestyle = :dot, linewidth = 1.0,
           label = (j == length(NS)) ? "30% threshold" : "")
    push!(dist_panels, plt)
end

plt_dist = plot(dist_panels...;
    layout = (1, length(NS)),
    size   = (1300, 380),
    dpi    = 200,
    plot_title = "Normalised \$L^2(\\nu)\$ distance of \$\\alpha_{\\mathrm{diag}}\$ from CKG smooth-Metro reference",
    plot_titlefontsize = 13,
    plot_titlevspan    = 0.07,
)

const OUT_DIST_PNG = joinpath(OUT_DIR, "alpha_diagonal_distance_vs_beta.png")
const OUT_DIST_PDF = joinpath(OUT_DIR, "alpha_diagonal_distance_vs_beta.pdf")
savefig(plt_dist, OUT_DIST_PNG)
savefig(plt_dist, OUT_DIST_PDF)
@printf("Saved %s\n", relpath(OUT_DIST_PNG, dirname(@__DIR__)))
@printf("Saved %s\n", relpath(OUT_DIST_PDF, dirname(@__DIR__)))

println("\nDone.")
