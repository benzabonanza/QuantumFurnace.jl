#!/usr/bin/env julia
#
# FAIR.7 — synthesise the fair-comparison findings (qf-mto.4, .5, .6) into a
# single thesis-quality 3-panel figure.
#
#   Panel (a) — left:    ρ_intrinsic = λ/Λ_max  vs  β  at fixed n=3
#                        Three lines (one per sampler). Tests H1 (CKG ≈ DLL
#                        Metro) and shows the DLL Gauss low-T collapse.
#
#   Panel (b) — middle:  α_diagonal(ν) overlay at  n=3, β=20.
#                        Three log-y curves. Visualises the diagonal-of-α
#                        driver behind the τ_mix tie (H2) and the Gauss
#                        super-exponential collapse outside |βν| ≲ 2.
#
#   Panel (c) — right:   τ_pred (FAIR.6, λ-from-gap) vs τ_meas scatter.
#                        Pooled best-fit dashed line + slope-1 dotted reference.
#                        Three colours, 32 points total. Slope and intercept
#                        annotated on the plot.
#
# Inputs  (all produced by FAIR.4-6):
#   scripts/output/fair_comparison_kms_dirichlet.bson   (FAIR.4 — ρ_intrinsic)
#   scripts/output/alpha_diagonal_analysis.bson         (FAIR.5 — α_diag(ν))
#   scripts/output/tau_mix_bound_check.bson             (FAIR.6 — τ_pred/meas)
#
# Outputs:
#   drafts/figures/numerics/fair_comparison_summary.png
#   drafts/figures/numerics/fair_comparison_summary.pdf
#
# Run:  julia --project scripts/plot_fair_comparison_summary.jl
#
# PHYSICS CHECK: panel (b) fixed at β=20 because that's where the three-way
# separation is largest — DLL Gauss α drops below 1e-9 in the wings while CKG
# sM and DLL Metro stay > 1e-3. β=10 also separates but less dramatically.
#
# PHYSICS CHECK: panel (a) plots all four available β values {1, 5, 10, 20}.
# The CKG ≈ DLL Metro gap stays within ~10% across β; DLL Gauss drops by ~5×
# at β=20 only — quintessential super-exp filter behaviour.

using Printf
using BSON
ENV["GKSwstype"] = "100"  # headless GR — required for PNG savefig in sandbox
using Plots

# ──────────────────────────────────────────────────────────────────────────────
# Load inputs
# ──────────────────────────────────────────────────────────────────────────────

const FAIR4_BSON = joinpath(@__DIR__, "output", "fair_comparison_kms_dirichlet.bson")
const FAIR5_BSON = joinpath(@__DIR__, "output", "alpha_diagonal_analysis.bson")
const FAIR6_BSON = joinpath(@__DIR__, "output", "tau_mix_bound_check.bson")

@assert isfile(FAIR4_BSON) "Missing FAIR.4 BSON: $FAIR4_BSON — run scratch_fair_comparison_dirichlet_sweep.jl first."
@assert isfile(FAIR5_BSON) "Missing FAIR.5 BSON: $FAIR5_BSON — run scratch_alpha_diagonal_analysis.jl first."
@assert isfile(FAIR6_BSON) "Missing FAIR.6 BSON: $FAIR6_BSON — run scratch_tau_mix_bound_check.jl first."

raw4 = BSON.load(FAIR4_BSON)
raw5 = BSON.load(FAIR5_BSON)
raw6 = BSON.load(FAIR6_BSON)

const RESULTS_FAIR4 = raw4[:results]            # records with (:n, :β, :sampler_key, :ρ_intrinsic, …)
const BETAS_FAIR4   = raw4[:beta_values]        # [1.0, 5.0, 10.0, 20.0]

const RECORDS_FAIR5 = raw5[:records]            # records with (:n, :beta, :sampler, :ν_grid, :α_diag, …)

const ROWS_FAIR6    = raw6[:rows]               # records with (:sampler, :n, :β, :τ_pred, :τ_meas, …)
const FITS_FAIR6    = raw6[:fits]               # Dict with :pooled, :ckg_smooth_metro, :dll_metro, :dll_gauss

@printf("Loaded FAIR.4 (%d records), FAIR.5 (%d records), FAIR.6 (%d rows + %d fits).\n",
        length(RESULTS_FAIR4), length(RECORDS_FAIR5), length(ROWS_FAIR6), length(FITS_FAIR6))

# ──────────────────────────────────────────────────────────────────────────────
# Sampler labels, keys, palette — single source of truth across panels
# ──────────────────────────────────────────────────────────────────────────────
#
# The three input BSONs use slightly inconsistent sampler labels:
#   FAIR.4 (results.sampler_key) : :ckg_smooth_metro, :dll_metro, :dll_gauss
#   FAIR.5 (records.sampler)     : "CKG smooth-Metro", "DLL Metropolis", "DLL Gaussian"
#   FAIR.6 (rows.sampler)        : :ckg_smooth_metro, :dll_metro, :dll_gauss
# We canonicalise on the symbol keys and provide string labels for legends.

const SAMPLER_KEYS  = [:ckg_smooth_metro, :dll_metro, :dll_gauss]

const SAMPLER_LABEL = Dict(
    :ckg_smooth_metro => "CKG smooth-Metro",
    :dll_metro        => "DLL Metropolis",
    :dll_gauss        => "DLL Gaussian",
)

# FAIR.5 uses "DLL Metropolis"/"DLL Gaussian" as record.sampler strings;
# match those exactly when querying α_diagonal records.
const SAMPLER_FAIR5_STRING = Dict(
    :ckg_smooth_metro => "CKG smooth-Metro",
    :dll_metro        => "DLL Metropolis",
    :dll_gauss        => "DLL Gaussian",
)

# Thesis colours: same palette as scripts/plot_alpha_diagonal.jl and
# scripts/numerics_ckg_vs_dll_taumix_comparison.jl.
const COLOR_CKG       = "#003147"   # deep navy
const COLOR_DLL_METRO = "#089099"   # teal-cyan
const COLOR_DLL_GAUSS = "#7CCBA2"   # mint-green

const SAMPLER_COLOR = Dict(
    :ckg_smooth_metro => COLOR_CKG,
    :dll_metro        => COLOR_DLL_METRO,
    :dll_gauss        => COLOR_DLL_GAUSS,
)

const SAMPLER_STYLE = Dict(
    :ckg_smooth_metro => :solid,
    :dll_metro        => :dash,
    :dll_gauss        => :dashdot,
)

const SAMPLER_MARKER = Dict(
    :ckg_smooth_metro => :circle,
    :dll_metro        => :diamond,
    :dll_gauss        => :star5,
)

# ──────────────────────────────────────────────────────────────────────────────
# Plot defaults — match thesis-numerics conventions
# ──────────────────────────────────────────────────────────────────────────────

default(
    fontfamily      = "Computer Modern",
    titlefontsize   = 11,
    guidefontsize   = 10,
    tickfontsize    = 8,
    legendfontsize  = 8,
    framestyle      = :box,
    grid            = true,
    gridlinewidth   = 0.4,
    gridstyle       = :dot,
    gridalpha       = 0.4,
    margin          = 3Plots.mm,
)

# ──────────────────────────────────────────────────────────────────────────────
# Panel (a) — ρ_intrinsic vs β at n=3
# ──────────────────────────────────────────────────────────────────────────────

const PANEL_A_N = 3                                          # PHYSICS CHECK: n=3 is the smallest cell where all three samplers are well-resolved
const PANEL_A_BETAS = sort(BETAS_FAIR4)                      # {1, 5, 10, 20}

function rho_for(sampler_key::Symbol, n::Int, β::Real)
    idx = findfirst(r -> r.n == n && isapprox(r.β, β; atol=1e-12) && r.sampler_key == sampler_key,
                    RESULTS_FAIR4)
    idx === nothing && error("ρ_intrinsic missing for sampler=$sampler_key n=$n β=$β")
    return RESULTS_FAIR4[idx].ρ_intrinsic
end

function build_panel_a(; show_legend::Bool)
    plt = plot(
        xlabel = raw"$\beta$",
        ylabel = raw"$\rho_{\mathrm{intrinsic}} \;=\; \lambda \,/\, \Lambda_{\max}$",
        title  = "(a)  Intrinsic mixing ratio at \$n=$(PANEL_A_N)\$",
        xticks = (PANEL_A_BETAS, [string(Int(β)) for β in PANEL_A_BETAS]),
        legend = show_legend ? :bottomleft : false,
        xscale = :log10,                                     # 4 points span β=1→20; log-x makes the spacing readable
        ylims  = (0.0, 0.36),                                # all ρ values lie in [0, 0.31]; small headroom for label
    )
    for k in SAMPLER_KEYS
        ρs = [rho_for(k, PANEL_A_N, β) for β in PANEL_A_BETAS]
        plot!(plt, PANEL_A_BETAS, ρs;
            label     = SAMPLER_LABEL[k],
            color     = SAMPLER_COLOR[k],
            linestyle = SAMPLER_STYLE[k],
            linewidth = 1.8,
            marker    = SAMPLER_MARKER[k],
            markersize = 6,
            markerstrokewidth = 0,
        )
    end
    return plt
end

# ──────────────────────────────────────────────────────────────────────────────
# Panel (b) — α_diag(ν) overlay at n=3, β=20
# ──────────────────────────────────────────────────────────────────────────────

const PANEL_B_N = 3
const PANEL_B_BETA = 20.0                                    # PHYSICS CHECK: β=20 maximises the 3-way α_diag separation
const Y_FLOOR = 1e-10                                        # clamp tiny α values to avoid log(0)

function alpha_diag_record(sampler_key::Symbol, n::Int, β::Real)
    label = SAMPLER_FAIR5_STRING[sampler_key]
    idx = findfirst(r -> r.n == n && r.beta == β && r.sampler == label, RECORDS_FAIR5)
    idx === nothing && error("α_diagonal missing for sampler=$label n=$n β=$β")
    return RECORDS_FAIR5[idx]
end

function build_panel_b(; show_legend::Bool)
    # Pre-fetch all three records to size the y-axis sensibly
    records = [alpha_diag_record(k, PANEL_B_N, PANEL_B_BETA) for k in SAMPLER_KEYS]
    all_y = Float64[]
    for r in records
        append!(all_y, r.α_diag)
    end
    y_max = maximum(all_y) * 2
    y_min = max(Y_FLOOR, minimum(filter(>(0), all_y)) / 2)

    plt = plot(
        xlabel = raw"$\nu$",
        ylabel = raw"$\alpha_{\mathrm{diag}}(\nu) \;=\; \sum_a \alpha^a(\nu, \nu)$",
        title  = "(b)  Kossakowski diagonal at \$n=$(PANEL_B_N),\\ \\beta=$(Int(PANEL_B_BETA))\$",
        yscale = :log10,
        ylims  = (y_min, y_max),
        legend = show_legend ? :bottomleft : false,
    )
    for (k, r) in zip(SAMPLER_KEYS, records)
        ν = r.ν_grid
        α = max.(r.α_diag, Y_FLOOR)                          # clamp for log scale
        plot!(plt, ν, α;
            label     = SAMPLER_LABEL[k],
            color     = SAMPLER_COLOR[k],
            linestyle = SAMPLER_STYLE[k],
            linewidth = 1.8,
            marker    = SAMPLER_MARKER[k],
            markersize = 2.5,
            markerstrokewidth = 0,
        )
    end
    return plt
end

# ──────────────────────────────────────────────────────────────────────────────
# Panel (c) — τ_pred vs τ_meas scatter (FAIR.6)
# ──────────────────────────────────────────────────────────────────────────────

function build_panel_c(; show_legend::Bool)
    # Collect all valid points per sampler (32 total).
    pts_per_sampler = Dict(k => (Float64[], Float64[]) for k in SAMPLER_KEYS)
    for row in ROWS_FAIR6
        # Skip any row missing data (FAIR.6 already filters; defensive)
        if isfinite(row.τ_pred) && isfinite(row.τ_meas) && row.τ_pred > 0 && row.τ_meas > 0
            xs, ys = pts_per_sampler[row.sampler]
            push!(xs, row.τ_meas)                             # x: measured (truth)
            push!(ys, row.τ_pred)                             # y: predicted (bound)
        end
    end

    # Determine plot range from data, padded a little on each side in log space.
    all_x = Float64[]; all_y = Float64[]
    for (xs, ys) in values(pts_per_sampler)
        append!(all_x, xs); append!(all_y, ys)
    end
    @assert !isempty(all_x) "Panel (c) found 0 valid (τ_meas, τ_pred) points."
    lo = min(minimum(all_x), minimum(all_y)) * 0.7
    hi = max(maximum(all_x), maximum(all_y)) * 1.4

    plt = plot(
        xlabel = raw"$\tau_{\mathrm{measured}}$ from biexp fit",
        ylabel = raw"$\tau_{\mathrm{predicted}} = \log(1/(\varepsilon\,\sqrt{p_{\min}})) / \lambda$",
        title  = "(c)  Gap-bound validation",
        xscale = :log10,
        yscale = :log10,
        xlims  = (lo, hi),
        ylims  = (lo, hi),
        legend = show_legend ? :bottomright : false,
        aspect_ratio = :equal,
    )

    # y = x diagonal reference (dotted gray)
    plot!(plt, [lo, hi], [lo, hi];
        label = raw"$y=x$",
        color = "#7f7f7f",
        linestyle = :dot,
        linewidth = 1.4,
    )

    # Pooled best-fit line (in log–log space): log10(y) = slope * log10(x) + intercept
    pooled  = FITS_FAIR6[:pooled]
    slope_p = pooled[:slope]
    icpt_p  = pooled[:intercept]
    xs_fit  = exp10.(range(log10(lo), log10(hi); length=64))
    ys_fit  = @. exp10(slope_p * log10(xs_fit) + icpt_p)
    plot!(plt, xs_fit, ys_fit;
        label     = @sprintf("pooled fit: slope=%.2f, intercept=%.2f", slope_p, icpt_p),
        color     = "#222222",
        linestyle = :dash,
        linewidth = 1.6,
    )

    # Scatter: one colour per sampler. Restrict marker labels so the legend
    # has exactly {y=x, fit, sampler×3} = 5 entries.
    for k in SAMPLER_KEYS
        xs, ys = pts_per_sampler[k]
        scatter!(plt, xs, ys;
            label  = SAMPLER_LABEL[k],
            color  = SAMPLER_COLOR[k],
            marker = SAMPLER_MARKER[k],
            markersize = 6,
            markerstrokewidth = 0.6,
            markerstrokecolor = :black,
        )
    end

    return plt
end

# ──────────────────────────────────────────────────────────────────────────────
# Compose the 3-panel layout
# ──────────────────────────────────────────────────────────────────────────────
#
# Strategy: place the legend inside panel (a) so it doesn't compete with the
# fit annotation in panel (c). Panel (a)'s lower-left corner is empty (DLL
# Gauss collapses, but at low β all three start near 0.30), so :bottomleft
# fits cleanly.

println("\nBuilding panel (a) — ρ_intrinsic vs β …")
plt_a = build_panel_a(; show_legend = true)

println("Building panel (b) — α_diag(ν) at n=$(PANEL_B_N), β=$(Int(PANEL_B_BETA)) …")
plt_b = build_panel_b(; show_legend = false)

println("Building panel (c) — τ_pred vs τ_meas scatter …")
plt_c = build_panel_c(; show_legend = true)

# ──────────────────────────────────────────────────────────────────────────────
# Compose & save
# ──────────────────────────────────────────────────────────────────────────────

println("\nComposing 1×3 figure …")
fig = plot(plt_a, plt_b, plt_c;
    layout = (1, 3),
    size   = (1700, 540),
    dpi    = 200,
    plot_title = "Fair comparison of CKG vs DLL Lindbladians via the KMS Dirichlet form",
    plot_titlefontsize = 14,
    plot_titlevspan    = 0.06,
    bottom_margin = 5Plots.mm,
    left_margin   = 6Plots.mm,
    right_margin  = 4Plots.mm,
    top_margin    = 3Plots.mm,
)

const OUT_DIR = joinpath(@__DIR__, "..", "drafts", "figures", "numerics")
mkpath(OUT_DIR)
const OUT_PNG = joinpath(OUT_DIR, "fair_comparison_summary.png")
const OUT_PDF = joinpath(OUT_DIR, "fair_comparison_summary.pdf")

savefig(fig, OUT_PNG)
savefig(fig, OUT_PDF)
@printf("Saved %s\n", relpath(OUT_PNG, dirname(@__DIR__)))
@printf("Saved %s\n", relpath(OUT_PDF, dirname(@__DIR__)))

# ──────────────────────────────────────────────────────────────────────────────
# Headline diagnostics (printed to stdout — verifies the figure tells the story)
# ──────────────────────────────────────────────────────────────────────────────

println("\n" * "="^72)
println("Headline numbers (sanity check)")
println("="^72)

println("\n[Panel (a)] ρ_intrinsic at n=$(PANEL_A_N) across β:")
@printf("%-18s %-10s %-10s %-10s %-10s\n", "sampler", "β=1", "β=5", "β=10", "β=20")
for k in SAMPLER_KEYS
    ρs = [rho_for(k, PANEL_A_N, β) for β in PANEL_A_BETAS]
    @printf("%-18s %.4f     %.4f     %.4f     %.4f\n",
            SAMPLER_LABEL[k], ρs[1], ρs[2], ρs[3], ρs[4])
end
ρ_ckg_β20 = rho_for(:ckg_smooth_metro, PANEL_A_N, 20.0)
ρ_metro_β20 = rho_for(:dll_metro, PANEL_A_N, 20.0)
ρ_gauss_β20 = rho_for(:dll_gauss, PANEL_A_N, 20.0)
@printf("\n  → at β=20: CKG/Metro ratio = %.2f (should be ≈ 1 if H1 wins)\n",
        ρ_ckg_β20 / ρ_metro_β20)
@printf("  → at β=20: CKG/Gauss ratio = %.2f (Gauss collapse marker)\n",
        ρ_ckg_β20 / ρ_gauss_β20)

println("\n[Panel (b)] α_diag(ν) extremes at n=$(PANEL_B_N), β=$(Int(PANEL_B_BETA)):")
for k in SAMPLER_KEYS
    r = alpha_diag_record(k, PANEL_B_N, PANEL_B_BETA)
    α_max = maximum(r.α_diag)
    α_min = minimum(r.α_diag)
    α_min_pos = minimum(filter(>(0), r.α_diag))
    @printf("  %-18s  α_max = %.3e   α_min(>0) = %.3e   ratio = %.1e\n",
            SAMPLER_LABEL[k], α_max, α_min_pos, α_max / α_min_pos)
end

println("\n[Panel (c)] Pooled gap-bound regression (FAIR.6):")
pooled = FITS_FAIR6[:pooled]
@printf("  n_points = %d\n", pooled[:n])
@printf("  slope     = %.3f   (1.0 means τ_pred ∝ τ_meas)\n", pooled[:slope])
@printf("  intercept = %.3f   (10^intercept = %.2f× constant-factor loose)\n",
        pooled[:intercept], 10^pooled[:intercept])
@printf("  R²        = %.3f\n", pooled[:r2])

println("\nDone.")
