#!/usr/bin/env julia
# analyze_betaphys_extended_scaling.jl  (qf-e4z.22)
#
# Consumes the v4 + v5 CKG smooth-Metro Heisenberg sidecars at
#   scripts/output/sweep_S1_v4_ckg_ideal/smooth_metro_eps1e-03/
#   scripts/output/sweep_S2_v4_dll_ideal/smooth_metro_eps1e-03/
# and runs:
#   (a) `fit_scaling` for M0 (separable power law) vs M1 (power × Arrhenius)
#       on τ_mix(n, β_phys); AICc weights, formula strings, residuals.
#   (b) gap_phys = gap_arnoldi · rescaling_factor analysis:
#       — gap_phys vs n at each β_phys (table + log fit slope, qualitative)
#       — gap_phys vs β_phys at each n   (table + sanity check vs O(1) prior)
#       gap_arnoldi (= gap_alg) is the rescaled-spectrum gap and IS NOT the
#       thesis-relevant gap; CLAUDE.md / [[krylov_x0_symmetric_bug_qf_8fr]].
#
# Produces:
#   drafts/figures/numerics/heis1d_betaphys_extended_taumix.{png,pdf}
#   drafts/figures/numerics/heis1d_betaphys_extended_gap.{png,pdf}
# and prints a markdown-ready table on stdout for paste-into-summary.
#
# Usage:
#   julia --project scripts/analyze_betaphys_extended_scaling.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using QuantumFurnace
using LinearAlgebra
using BSON
using Printf
using Dates
using Statistics

println("[init] $(now())")

const S1_DIR = joinpath(@__DIR__, "output", "sweep_S1_v4_ckg_ideal", "smooth_metro_eps1e-03")
const S2_DIR = joinpath(@__DIR__, "output", "sweep_S2_v4_dll_ideal", "smooth_metro_eps1e-03")
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
    sort!(rows; by = r -> (r.n, r.beta_phys))
    return rows
end

s1 = load_sidecars(S1_DIR)
s2 = load_sidecars(S2_DIR)
@printf("[load] S1 (CKG / EnergyDomain): %d sidecars\n", length(s1))
@printf("[load] S2 (DLL / BohrDomain):   %d sidecars\n", length(s2))

# --- Per-cell table ------------------------------------------------------

function gap_phys(r)
    haskey(r, :gap_arnoldi) || return NaN
    haskey(r, :rescaling_factor) || return NaN
    return r.gap_arnoldi * r.rescaling_factor
end

function print_table(rows, label)
    println("\n", "="^115)
    @printf("%s  (%d cells)\n", label, length(rows))
    println("="^115)
    @printf("%-3s %-7s %-7s %-9s %-3s %-11s %-11s %-11s %-9s %-8s %-7s %-6s\n",
            "n", "β_phys", "β_alg", "R(n)", "r_D", "gap_alg", "gap_phys", "τ_mix",
            "floor", "wall_s", "conv", "ver")
    println("-"^115)
    for r in rows
        ver = string(get(r, :sweep_version, :v4))
        @printf("%-3d %-7.2f %-7.2f %-9.3f %-3d %-11.5g %-11.4g %-11.4g %-9.3g %-8.2f %-7s %-6s\n",
                r.n, r.beta_phys, r.beta_alg, r.rescaling_factor,
                get(r, :r_D, -1), r.gap_arnoldi, gap_phys(r), r.mixing_time,
                r.floor_distance, r.wall_time, string(r.all_converged), ver)
    end
end

print_table(s1, "S1 — CKG smooth-Metropolis (EnergyDomain, Krylov)")
print_table(s2, "S2 — DLL Metropolis (BohrDomain, analytical)")

# --- Markdown-ready table for the draft ---------------------------------

function print_md_table(rows)
    println("\n## Markdown table (paste into draft):\n")
    println("| n | β_phys | β_alg | R(n) | r_D | gap_alg | gap_phys | τ_mix | floor | source |")
    println("|---|---|---|---|---|---|---|---|---|---|")
    for r in rows
        @printf("| %d | %.2f | %.2f | %.3f | %d | %.4g | %.3g | %.3g | %.2g | %s |\n",
                r.n, r.beta_phys, r.beta_alg, r.rescaling_factor,
                get(r, :r_D, -1), r.gap_arnoldi, gap_phys(r), r.mixing_time,
                r.floor_distance, string(r.mixing_time_source))
    end
end

# --- Scaling fit τ_mix ~ (n, β_phys) ------------------------------------

println("\n", "="^88)
println("Scaling fits  τ_mix(n, β_phys)")
println("="^88)

function run_scaling_fit(rows, label)
    println("\n--- $label ($(length(rows)) cells) ---")
    # Use beta_kind=:phys explicitly to make the convention unambiguous in the output.
    try
        fits = fit_scaling(rows; beta_kind = :phys)
        ranking = compare_models(fits)
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
        return fits, ranking
    catch err
        @warn "fit_scaling failed for $label" err
        return nothing, nothing
    end
end

fits_s1, rank_s1 = run_scaling_fit(s1, "S1 (CKG)")
fits_s2, rank_s2 = run_scaling_fit(s2, "S2 (DLL)")

# --- Gap_phys analyses ---------------------------------------------------

function gap_analysis(rows, label)
    println("\n", "="^88)
    @printf("gap_phys analysis — %s\n", label)
    println("="^88)

    by_n = Dict{Int, Vector{NamedTuple}}()
    by_b = Dict{Float64, Vector{NamedTuple}}()
    for r in rows
        push!(get!(by_n, r.n, NamedTuple[]), r)
        push!(get!(by_b, r.beta_phys, NamedTuple[]), r)
    end

    println("\n  gap_phys vs n at each β_phys:")
    for β in sort!(collect(keys(by_b)))
        cells = sort(by_b[β]; by = r -> r.n)
        gps = [gap_phys(r) for r in cells]
        ns  = [r.n for r in cells]
        @printf("    β_phys=%.2f  n=%s  gap_phys=[%s]  ratio max/min=%.2f\n",
                β, ns, join((@sprintf("%.3g", g) for g in gps), ", "),
                maximum(gps) / max(minimum(gps), 1e-15))
        # Crude log-slope across n (only meaningful when ≥3 n points).
        if length(ns) ≥ 3
            xs = log.(Float64.(ns))
            ys = log.(gps)
            x̄ = mean(xs); ȳ = mean(ys)
            slope = sum((xs .- x̄) .* (ys .- ȳ)) / sum((xs .- x̄).^2)
            @printf("      → log-log slope (gap_phys ~ n^slope) ≈ %+.3f\n", slope)
        end
    end

    println("\n  gap_phys vs β_phys at each n:")
    for n in sort!(collect(keys(by_n)))
        cells = sort(by_n[n]; by = r -> r.beta_phys)
        gps = [gap_phys(r) for r in cells]
        βs  = [r.beta_phys for r in cells]
        @printf("    n=%d  β_phys=%s  gap_phys=[%s]  ratio max/min=%.2f\n",
                n, βs, join((@sprintf("%.3g", g) for g in gps), ", "),
                maximum(gps) / max(minimum(gps), 1e-15))
    end
end

gap_analysis(s1, "S1 (CKG)")
gap_analysis(s2, "S2 (DLL)")

# --- Print markdown summary for the draft -------------------------------

println("\n", "="^88)
println("Markdown summary (S1 — CKG)")
println("="^88)
print_md_table(s1)

println("\n", "="^88)
println("Markdown summary (S2 — DLL)")
println("="^88)
print_md_table(s2)

# --- Diagnostic figures ------------------------------------------------

println("\n[plot] preparing figures…")
ENV["GKSwstype"] = "100"
using Plots

# Two-panel τ_mix vs β_phys (color = n) and best-fit overlay.
function fig_taumix(rows, fits, ranking, fig_basename, title_)
    by_n = Dict{Int, Vector{NamedTuple}}()
    for r in rows
        push!(get!(by_n, r.n, NamedTuple[]), r)
    end
    p1 = Plots.plot(xscale = :log10, yscale = :log10,
        xlabel = "β_phys", ylabel = "τ_mix", title = title_,
        legend = :outerright, size = (900, 450))
    ns_sorted = sort(collect(keys(by_n)))
    for n in ns_sorted
        cells = sort(by_n[n]; by = r -> r.beta_phys)
        Plots.scatter!(p1, [r.beta_phys for r in cells], [r.mixing_time for r in cells];
            label = "n=$(n)", markersize = 5)
    end
    # Overlay best-fit lines (M0) at each n.
    if fits !== nothing && haskey(fits, :M0)
        f = fits[:M0]
        β_grid = exp.(range(log(minimum(r.beta_phys for r in rows)),
                            log(maximum(r.beta_phys for r in rows));
                            length = 50))
        for n in ns_sorted
            τ_pred = [predict_scaling(f, n, β) for β in β_grid]
            Plots.plot!(p1, β_grid, τ_pred; label = "", linestyle = :dash, alpha = 0.5)
        end
    end
    Plots.savefig(p1, joinpath(FIG_DIR, "$(fig_basename).png"))
    Plots.savefig(p1, joinpath(FIG_DIR, "$(fig_basename).pdf"))
    @info "Wrote τ_mix plot" base=fig_basename
end

fig_taumix(s1, fits_s1, rank_s1, "heis1d_betaphys_extended_taumix_s1_ckg",
           "τ_mix vs β_phys — 1D Heisenberg disordered, CKG smooth-Metro")
fig_taumix(s2, fits_s2, rank_s2, "heis1d_betaphys_extended_taumix_s2_dll",
           "τ_mix vs β_phys — 1D Heisenberg disordered, DLL Metropolis")

# Gap_phys: 2-panel (gap_phys vs n / gap_phys vs β_phys).
function fig_gap(rows, fig_basename, title_)
    by_n = Dict{Int, Vector{NamedTuple}}()
    by_b = Dict{Float64, Vector{NamedTuple}}()
    for r in rows
        push!(get!(by_n, r.n, NamedTuple[]), r)
        push!(get!(by_b, r.beta_phys, NamedTuple[]), r)
    end
    p1 = Plots.plot(xlabel = "n", ylabel = "gap_phys",
        title = "$(title_) — vs n",
        legend = :outerright)
    for β in sort(collect(keys(by_b)))
        cells = sort(by_b[β]; by = r -> r.n)
        Plots.plot!(p1, [r.n for r in cells], [gap_phys(r) for r in cells];
            label = "β_phys=$(β)", marker = :circle, markersize = 4)
    end
    Plots.hline!(p1, [0]; color = :black, linestyle = :dot, label = "")

    p2 = Plots.plot(xlabel = "β_phys", ylabel = "gap_phys",
        xscale = :log10,
        title = "$(title_) — vs β_phys",
        legend = :outerright)
    for n in sort(collect(keys(by_n)))
        cells = sort(by_n[n]; by = r -> r.beta_phys)
        Plots.plot!(p2, [r.beta_phys for r in cells], [gap_phys(r) for r in cells];
            label = "n=$(n)", marker = :circle, markersize = 4)
    end
    fig = Plots.plot(p1, p2; layout = (1, 2), size = (1300, 500))
    Plots.savefig(fig, joinpath(FIG_DIR, "$(fig_basename).png"))
    Plots.savefig(fig, joinpath(FIG_DIR, "$(fig_basename).pdf"))
    @info "Wrote gap_phys plot" base=fig_basename
end

fig_gap(s1, "heis1d_betaphys_extended_gap_s1_ckg",
        "gap_phys (CKG smooth-Metro)")
fig_gap(s2, "heis1d_betaphys_extended_gap_s2_dll",
        "gap_phys (DLL Metropolis)")

println("\n[done] $(now())")
