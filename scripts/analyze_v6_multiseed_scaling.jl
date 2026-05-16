#!/usr/bin/env julia
# analyze_v6_multiseed_scaling.jl  (qf-e4z.23)
#
# Aggregates the v6 multi-seed CKG smooth-Metro Heisenberg sidecars at
#   scripts/output/sweep_S1_v6_ckg_ideal_multiseed/smooth_metro_eps1e-03/
# and runs:
#   (a) median + IQR (and min/max) of gap_phys and τ_mix across 5 seeds per
#       (n, β_phys) cell.
#   (b) `fit_scaling` for M0 (separable power) vs M1 (Arrhenius) on the
#       per-cell MEDIAN τ_mix(n, β_phys). AICc weights, residuals, formula.
#   (c) Even/odd-n diagnostic: per-β_phys, is the inter-parity gap_phys spread
#       at matched n LARGER than the intra-parity seed scatter? If not, no
#       even/odd structure can be claimed (per
#       [[feedback_more_data_points_for_scaling_claims]]).
#
# Produces:
#   drafts/figures/numerics/v6_multiseed_gap_phys.{png,pdf}
#   drafts/figures/numerics/v6_multiseed_taumix.{png,pdf}
#   drafts/figures/numerics/v6_multiseed_even_odd.{png,pdf}
# and prints a markdown-ready summary on stdout.
#
# Usage:
#   julia --project scripts/analyze_v6_multiseed_scaling.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using QuantumFurnace
using LinearAlgebra
using BSON
using Printf
using Dates
using Statistics

println("[init] $(now())")

const S1_DIR = joinpath(@__DIR__, "output", "sweep_S1_v6_ckg_ideal_multiseed", "smooth_metro_eps1e-03")
const FIG_DIR = joinpath(@__DIR__, "..", "drafts", "figures", "numerics")
mkpath(FIG_DIR)

# --- Sidecar loader ------------------------------------------------------

function load_sidecars(dir)
    rows = NamedTuple[]
    for fname in sort(readdir(dir))
        endswith(fname, ".bson") || continue
        path = joinpath(dir, fname)
        try
            d = BSON.load(path, QuantumFurnace)
            r = d[:result]
            push!(rows, (; (Symbol(k) => v for (k, v) in pairs(r))...))
        catch err
            @warn "Failed to load sidecar" path err
        end
    end
    sort!(rows; by = r -> (r.n, r.beta_phys, r.seed))
    return rows
end

rows = load_sidecars(S1_DIR)
@printf("[load] %d sidecars\n", length(rows))

# --- Helpers -------------------------------------------------------------

gap_phys(r) = r.gap_arnoldi * r.rescaling_factor

# Group by (n, beta_phys) — each group should have up to 5 seeds.
function group_by_cell(rows)
    cells = Dict{Tuple{Int, Float64}, Vector{NamedTuple}}()
    for r in rows
        key = (r.n, r.beta_phys)
        push!(get!(cells, key, NamedTuple[]), r)
    end
    return cells
end

cells = group_by_cell(rows)

# --- Per-cell aggregate table -------------------------------------------

struct CellStats
    n::Int
    beta_phys::Float64
    n_seeds::Int
    median_gap_phys::Float64
    iqr_gap_phys::Float64
    min_gap_phys::Float64
    max_gap_phys::Float64
    median_tau_mix::Float64
    iqr_tau_mix::Float64
    min_tau_mix::Float64
    max_tau_mix::Float64
    median_rescale::Float64
end

function _iqr(xs)
    isempty(xs) && return NaN
    length(xs) == 1 && return 0.0
    return quantile(xs, 0.75) - quantile(xs, 0.25)
end

function aggregate_cells(cells)
    stats = CellStats[]
    for (k, cell) in cells
        n, β_phys = k
        gps = [gap_phys(r) for r in cell]
        τs  = [r.mixing_time for r in cell]
        rs  = [r.rescaling_factor for r in cell]
        push!(stats, CellStats(n, β_phys, length(cell),
            median(gps), _iqr(gps), minimum(gps), maximum(gps),
            median(τs),  _iqr(τs),  minimum(τs),  maximum(τs),
            median(rs)))
    end
    sort!(stats; by = s -> (s.n, s.beta_phys))
    return stats
end

stats = aggregate_cells(cells)

function print_cell_table(stats)
    println("\n", "="^140)
    println("Per-cell aggregate (5-seed median + IQR)")
    println("="^140)
    @printf("%-3s %-7s %-3s %-9s %-12s %-12s %-12s %-12s %-12s %-12s\n",
            "n", "β_phys", "#s", "med R",
            "med gap_phys", "IQR gap", "min..max gap",
            "med τ_mix", "IQR τ_mix", "min..max τ_mix")
    println("-"^140)
    for s in stats
        @printf("%-3d %-7.2f %-3d %-9.3f %-12.4g %-12.3g %-12s %-12.4g %-12.3g %-12s\n",
                s.n, s.beta_phys, s.n_seeds, s.median_rescale,
                s.median_gap_phys, s.iqr_gap_phys,
                @sprintf("%.3g..%.3g", s.min_gap_phys, s.max_gap_phys),
                s.median_tau_mix, s.iqr_tau_mix,
                @sprintf("%.3g..%.3g", s.min_tau_mix, s.max_tau_mix))
    end
end

print_cell_table(stats)

# --- Markdown-ready table -----------------------------------------------

function print_md_table(stats)
    println("\n## Markdown table (paste into draft):\n")
    println("| n | β_phys | #seeds | med R | med gap_phys | IQR | med τ_mix | IQR |")
    println("|---|---|---|---|---|---|---|---|")
    for s in stats
        @printf("| %d | %.2f | %d | %.3f | %.4g | %.3g | %.4g | %.3g |\n",
                s.n, s.beta_phys, s.n_seeds, s.median_rescale,
                s.median_gap_phys, s.iqr_gap_phys,
                s.median_tau_mix, s.iqr_tau_mix)
    end
end

print_md_table(stats)

# --- Scaling fit on per-cell median --------------------------------------

# Synthesize a NamedTuple row per (n, β_phys) cell using MEDIAN τ_mix and β_phys,
# so fit_scaling sees one effective datapoint per cell.

function rows_from_medians(stats)
    return NamedTuple[
        (
            n = s.n,
            beta_phys = s.beta_phys,
            mixing_time = s.median_tau_mix,
            mixing_time_source = :extrapolated,  # synthesized from medians, treat as extrapolated
            beta = s.beta_phys,                  # so legacy fallback picks the right column
        )
        for s in stats
    ]
end

println("\n", "="^88)
println("Scaling fit on per-cell MEDIAN τ_mix(n, β_phys)")
println("="^88)

median_rows = rows_from_medians(stats)
fits = nothing
ranking = nothing
try
    global fits = fit_scaling(median_rows; beta_kind = :phys)
    global ranking = compare_models(fits)
    println("Model ranking (best AICc first):")
    for (i, m) in enumerate(ranking.ranked)
        @printf("  %d. %-3s  AICc=%-8.3f  Δ=%-7.3f  weight=%.3f\n",
                i, string(m), ranking.aicc[i], ranking.delta_aicc[i], ranking.weights[i])
    end
    for m in (:M0, :M1)
        haskey(fits, m) || continue
        f = fits[m]
        println("  $(m): ", formula_string(f), "  σ_residual=", round(f.sigma_residual; digits=4))
    end
catch err
    @warn "fit_scaling on median failed" err
end

# --- Even/odd-n diagnostic -----------------------------------------------

println("\n", "="^88)
println("Even/odd-n diagnostic (per β_phys)")
println("="^88)
println("If max{intra-parity seed scatter} ≥ |even − odd median spread at matched n|,")
println("the even/odd distinction is consistent with noise.\n")

function even_odd_diagnostic(stats, rows)
    by_b = Dict{Float64, Vector{CellStats}}()
    for s in stats
        push!(get!(by_b, s.beta_phys, CellStats[]), s)
    end
    for β in sort!(collect(keys(by_b)))
        cells = sort(by_b[β]; by = s -> s.n)
        # For each n: median gap_phys, plus the per-seed max scatter (max - min) /
        # 2 as a "half-range" proxy for intra-parity spread.
        evens = filter(s -> iseven(s.n), cells)
        odds  = filter(s ->  isodd(s.n), cells)

        # Intra-parity half-range across all even-n cells: max gap_phys / min gap_phys
        # at each parity class (across n, but at fixed β_phys).
        function summarize(group, label)
            isempty(group) && return
            meds = [s.median_gap_phys for s in group]
            scat = [(s.max_gap_phys - s.min_gap_phys) / max(2 * s.median_gap_phys, 1e-30) for s in group]
            @printf("    %s n=%s  median gap_phys=[%s]  per-cell half-range/median: %s\n",
                    label,
                    [s.n for s in group],
                    join((@sprintf("%.3g", g) for g in meds), ", "),
                    join((@sprintf("%.2f", x) for x in scat), ", "))
        end

        @printf("  β_phys = %.2f\n", β)
        summarize(evens, "even")
        summarize(odds,  "odd ")

        # Geometric mean of medians per parity, ratio between parities.
        if !isempty(evens) && !isempty(odds)
            geo_even = exp(mean(log.([s.median_gap_phys for s in evens])))
            geo_odd  = exp(mean(log.([s.median_gap_phys for s in odds])))
            ratio = geo_odd / geo_even
            # Typical intra-parity scatter (geomean of half-ranges/median over both parities)
            half_ranges = vcat(
                [(s.max_gap_phys - s.min_gap_phys) / max(2 * s.median_gap_phys, 1e-30) for s in evens],
                [(s.max_gap_phys - s.min_gap_phys) / max(2 * s.median_gap_phys, 1e-30) for s in odds],
            )
            scatter_typ = isempty(half_ranges) ? NaN : median(half_ranges)
            @printf("    → geomean ratio odd/even = %.3f  (intra-parity median half-range/median = %.3f)\n",
                    ratio, scatter_typ)
            log_ratio = abs(log(ratio))
            verdict = log_ratio > 2 * scatter_typ ? "REAL" : "consistent with noise"
            @printf("    → |log(ratio)| = %.3f vs 2·intra-scatter = %.3f  →  %s\n\n",
                    log_ratio, 2 * scatter_typ, verdict)
        end
    end
end

even_odd_diagnostic(stats, rows)

# --- Plots --------------------------------------------------------------

println("\n[plot] preparing figures…")
ENV["GKSwstype"] = "100"
using Plots

function fig_gap(stats, fig_basename, title_)
    by_n = Dict{Int, Vector{CellStats}}()
    by_b = Dict{Float64, Vector{CellStats}}()
    for s in stats
        push!(get!(by_n, s.n, CellStats[]), s)
        push!(get!(by_b, s.beta_phys, CellStats[]), s)
    end

    # Panel 1: gap_phys vs n (one curve per β_phys), even/odd shape distinction.
    p1 = Plots.plot(xlabel = "n", ylabel = "gap_phys",
        yscale = :log10,
        title = "$(title_) — vs n  (5-seed median + min..max band)",
        legend = :outerright)
    for β in sort(collect(keys(by_b)))
        cells = sort(by_b[β]; by = s -> s.n)
        ns = [s.n for s in cells]
        meds = [s.median_gap_phys for s in cells]
        lows = [s.min_gap_phys for s in cells]
        highs = [s.max_gap_phys for s in cells]
        Plots.plot!(p1, ns, meds;
            label = @sprintf("β_phys=%.2f", β),
            ribbon = (meds .- lows, highs .- meds),
            fillalpha = 0.15,
            marker = :auto, markersize = 5)
        # Mark even-n with squares and odd-n with diamonds for parity diagnosis.
        evens = [(s.n, s.median_gap_phys) for s in cells if iseven(s.n)]
        odds  = [(s.n, s.median_gap_phys) for s in cells if isodd(s.n)]
        if !isempty(evens)
            Plots.scatter!(p1, first.(evens), last.(evens);
                label = "", marker = :square, color = :black, markersize = 6, markeralpha = 0.4)
        end
        if !isempty(odds)
            Plots.scatter!(p1, first.(odds), last.(odds);
                label = "", marker = :diamond, color = :black, markersize = 6, markeralpha = 0.4)
        end
    end

    # Panel 2: gap_phys vs β_phys (one curve per n).
    p2 = Plots.plot(xlabel = "β_phys", ylabel = "gap_phys",
        xscale = :log10, yscale = :log10,
        title = "vs β_phys",
        legend = :outerright)
    for n in sort(collect(keys(by_n)))
        cells = sort(by_n[n]; by = s -> s.beta_phys)
        βs = [s.beta_phys for s in cells]
        meds = [s.median_gap_phys for s in cells]
        lows = [s.min_gap_phys for s in cells]
        highs = [s.max_gap_phys for s in cells]
        Plots.plot!(p2, βs, meds;
            label = "n=$(n)",
            ribbon = (meds .- lows, highs .- meds),
            fillalpha = 0.15,
            marker = :circle, markersize = 4)
    end

    fig = Plots.plot(p1, p2; layout = (1, 2), size = (1500, 550))
    Plots.savefig(fig, joinpath(FIG_DIR, "$(fig_basename).png"))
    Plots.savefig(fig, joinpath(FIG_DIR, "$(fig_basename).pdf"))
    @info "Wrote gap_phys plot" base=fig_basename
end

fig_gap(stats, "v6_multiseed_gap_phys", "gap_phys — 1D Heisenberg PBC multi-seed (CKG smooth-Metro)")

function fig_taumix(stats, fits, fig_basename, title_)
    by_n = Dict{Int, Vector{CellStats}}()
    for s in stats
        push!(get!(by_n, s.n, CellStats[]), s)
    end
    p1 = Plots.plot(xscale = :log10, yscale = :log10,
        xlabel = "β_phys", ylabel = "τ_mix",
        title = "$(title_)  (5-seed median + min..max band)",
        legend = :outerright, size = (1000, 550))
    ns_sorted = sort(collect(keys(by_n)))
    for n in ns_sorted
        cells = sort(by_n[n]; by = s -> s.beta_phys)
        βs = [s.beta_phys for s in cells]
        meds = [s.median_tau_mix for s in cells]
        lows = [s.min_tau_mix for s in cells]
        highs = [s.max_tau_mix for s in cells]
        Plots.plot!(p1, βs, meds;
            label = "n=$(n)",
            ribbon = (meds .- lows, highs .- meds),
            fillalpha = 0.15,
            marker = :circle, markersize = 5)
    end
    if fits !== nothing && haskey(fits, :M0)
        f = fits[:M0]
        β_grid = exp.(range(log(minimum(s.beta_phys for s in stats)),
                            log(maximum(s.beta_phys for s in stats));
                            length = 50))
        for n in ns_sorted
            τ_pred = [predict_scaling(f, n, β) for β in β_grid]
            Plots.plot!(p1, β_grid, τ_pred;
                label = "", linestyle = :dash, alpha = 0.5)
        end
    end
    Plots.savefig(p1, joinpath(FIG_DIR, "$(fig_basename).png"))
    Plots.savefig(p1, joinpath(FIG_DIR, "$(fig_basename).pdf"))
    @info "Wrote τ_mix plot" base=fig_basename
end

fig_taumix(stats, fits, "v6_multiseed_taumix",
           "τ_mix vs β_phys — 1D Heisenberg PBC multi-seed (CKG smooth-Metro)")

# Even/odd specific plot: gap_phys vs n with even/odd colored separately, at
# every β_phys (compact small-multiple).
function fig_even_odd(stats, fig_basename)
    by_b = Dict{Float64, Vector{CellStats}}()
    for s in stats
        push!(get!(by_b, s.beta_phys, CellStats[]), s)
    end
    βs_sorted = sort(collect(keys(by_b)))
    panels = []
    for β in βs_sorted
        cells = sort(by_b[β]; by = s -> s.n)
        ns = [s.n for s in cells]
        meds = [s.median_gap_phys for s in cells]
        lows = [s.min_gap_phys for s in cells]
        highs = [s.max_gap_phys for s in cells]
        even_idx = findall(iseven, ns)
        odd_idx  = findall(isodd, ns)
        p = Plots.plot(xlabel = "n", ylabel = "gap_phys",
            title = @sprintf("β_phys=%.2f", β),
            yscale = :log10,
            legend = (β == βs_sorted[1] ? :bottomleft : false))
        # Bands first so dots overlay.
        Plots.plot!(p, ns, meds;
            label = "", color = :gray, alpha = 0.4,
            ribbon = (meds .- lows, highs .- meds),
            fillalpha = 0.15)
        if !isempty(even_idx)
            Plots.scatter!(p, ns[even_idx], meds[even_idx];
                label = "even n", marker = :square, color = :blue, markersize = 6)
        end
        if !isempty(odd_idx)
            Plots.scatter!(p, ns[odd_idx], meds[odd_idx];
                label = "odd n", marker = :diamond, color = :red, markersize = 6)
        end
        push!(panels, p)
    end
    fig = Plots.plot(panels...; layout = (2, 3), size = (1400, 700))
    Plots.savefig(fig, joinpath(FIG_DIR, "$(fig_basename).png"))
    Plots.savefig(fig, joinpath(FIG_DIR, "$(fig_basename).pdf"))
    @info "Wrote even/odd plot" base=fig_basename
end

fig_even_odd(stats, "v6_multiseed_even_odd")

println("\n[done] $(now())")
