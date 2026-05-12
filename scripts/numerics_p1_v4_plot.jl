#!/usr/bin/env julia
#
# Plot P1 v4 — ideal-Lindbladian τ_mix: DLL Metropolis vs CKG smooth-Metropolis.
# Single panel at ε=1e-3 with the spectral-gap upper-bound overlay
# τ_bound = log(d/ε) / λ_2 (per the 2026-05-08 spec in qf-e4z.11).
#
# v4 vs v3: β_phys-first sidecars (file naming `betaphys*` instead of `beta*`),
# legend labels in β_phys, fresh post-fix data after qf-6vr / qf-9z0 / qf-2kd
# refactors. CKG filter is fixed s = 0.25 (decided 2026-05-12, see CLAUDE.md
# and quadrature-convergence-summary-v2.md). Skips cells whose sidecar isn't
# present yet — so the plot auto-grows as the n=8 / n=9 sweep cells land.
#
# Source data (no new simulations — pure post-processing):
#   S1 (qf-e4z.21): scripts/output/sweep_S1_v4_ckg_ideal/smooth_metro_eps1e-03/sweep_n*_betaphys*_seed42_L_KMS_Energy.bson
#   S2 (qf-e4z.21): scripts/output/sweep_S2_v4_dll_ideal/smooth_metro_eps1e-03/sweep_n*_betaphys*_seed42_L_DLL_Bohr.bson
#
# Output:
#   drafts/figures/numerics/p1_taumix_v4.{png, pdf}
#
# Usage: JULIA_NUM_THREADS=1 julia --project scripts/numerics_p1_v4_plot.jl

using Printf
using BSON
using QuantumFurnace                                # for BSON symbol resolution
ENV["GKSwstype"] = "100"
using Plots

# ── Paths ─────────────────────────────────────────────────────────────────────
const REPO_ROOT = joinpath(@__DIR__, "..")
const S1_DIR    = joinpath(REPO_ROOT, "scripts", "output", "sweep_S1_v4_ckg_ideal",
                            "smooth_metro_eps1e-03")
const S2_DIR    = joinpath(REPO_ROOT, "scripts", "output", "sweep_S2_v4_dll_ideal",
                            "smooth_metro_eps1e-03")
const OUT_DIR   = joinpath(REPO_ROOT, "drafts", "figures", "numerics")
const FIG_PNG   = joinpath(OUT_DIR, "p1_taumix_v4.png")
const FIG_PDF   = joinpath(OUT_DIR, "p1_taumix_v4.pdf")
mkpath(OUT_DIR)

# ── Sweep grid (matches v4 sweep) ────────────────────────────────────────────
const N_VALUES    = 3:8
const BETA_PHYS   = (0.25, 0.5, 1.0)
const TARGET_EPS  = 1e-3
const SEED        = 42

# ── Thesis colour palette (β_phys-gradient, warm → dark) ─────────────────────
# β_phys low (high T) = terracotta; β_phys high (low T) = aubergine.
const COLOR_BETA = Dict(
    0.25 => "#B5654A",   # terracotta
    0.5  => "#7A2E39",   # bordeaux
    1.0  => "#4F3B5C",   # aubergine
)

const MARKER_CKG = :circle
const MARKER_DLL = :rect

# ── Schema-aware BSON loader (β_phys filename tag) ────────────────────────────
function _betaphys_str(β::Real)
    s = @sprintf("%.6f", float(β))
    s = rstrip(s, '0'); s = rstrip(s, '.')
    return isempty(s) ? "0" : s
end

function _load_cell(dir::AbstractString, n::Integer, β_phys::Real,
                    construction_tag::AbstractString, domain_tag::AbstractString)
    path = joinpath(dir, "sweep_n$(n)_betaphys$(_betaphys_str(β_phys))_seed$(SEED)_L_$(construction_tag)_$(domain_tag).bson")
    isfile(path) || return nothing
    data = BSON.load(path, QuantumFurnace)
    return data[:result]
end

function _tau_and_bound(cell)
    cell === nothing && return (NaN, NaN)
    τ = get(cell, :mixing_time, NaN)
    τ_b = get(cell, :tau_mix_bound, NaN)
    return (Float64(τ), Float64(τ_b))
end

function _beta_alg(cell)
    cell === nothing && return NaN
    return Float64(get(cell, :beta_alg, NaN))
end

# ── Collect data ──────────────────────────────────────────────────────────────
println("="^72)
println("P1 v4 plot — DLL Metropolis vs CKG smooth-Metro (ideal Lindbladian)")
println("S1 dir: ", S1_DIR)
println("S2 dir: ", S2_DIR)
println("="^72)

results = Dict{Tuple{Symbol,Float64},Vector{NamedTuple{(:n,:tau,:bound,:beta_alg),Tuple{Int,Float64,Float64,Float64}}}}()
for β in BETA_PHYS
    results[(:CKG, β)] = NamedTuple{(:n,:tau,:bound,:beta_alg),Tuple{Int,Float64,Float64,Float64}}[]
    results[(:DLL, β)] = NamedTuple{(:n,:tau,:bound,:beta_alg),Tuple{Int,Float64,Float64,Float64}}[]
end

for n in N_VALUES, β in BETA_PHYS
    ckg_cell = _load_cell(S1_DIR, n, β, "KMS", "Energy")
    dll_cell = _load_cell(S2_DIR, n, β, "DLL", "Bohr")
    τ_ckg, b_ckg = _tau_and_bound(ckg_cell)
    τ_dll, b_dll = _tau_and_bound(dll_cell)
    β_alg_ckg = _beta_alg(ckg_cell)
    β_alg_dll = _beta_alg(dll_cell)
    push!(results[(:CKG, β)], (n=n, tau=τ_ckg, bound=b_ckg, beta_alg=β_alg_ckg))
    push!(results[(:DLL, β)], (n=n, tau=τ_dll, bound=b_dll, beta_alg=β_alg_dll))
    if isfinite(τ_ckg) || isfinite(τ_dll)
        @printf("n=%d  β_phys=%-5.2f β_alg=%.2f | CKG τ=%.3e  bound=%.3e | DLL τ=%.3e | ratio DLL/CKG=%.2f\n",
                n, β, isfinite(β_alg_ckg) ? β_alg_ckg : β_alg_dll,
                τ_ckg, b_ckg, τ_dll,
                (isfinite(τ_ckg) && τ_ckg > 0) ? τ_dll/τ_ckg : NaN)
    end
end

# ── Build the figure ──────────────────────────────────────────────────────────
plt = Plots.plot(
    xlabel = "n",
    ylabel = raw"$\tau_{\mathrm{mix}}$",
    yscale = :log10,
    legend = :outerright,
    xticks = collect(N_VALUES),
    framestyle = :box,
    guidefontsize = 11,
    tickfontsize = 9,
    legendfontsize = 8,
    size = (820, 520),
    left_margin = 4Plots.mm,
    bottom_margin = 4Plots.mm,
    right_margin = 2Plots.mm,
    top_margin = 2Plots.mm,
)

# Plot the simulated τ_mix curves (solid, filled markers).
for β in BETA_PHYS
    col = COLOR_BETA[β]

    rows_ckg = results[(:CKG, β)]
    ns_ckg = [r.n for r in rows_ckg if isfinite(r.tau) && r.tau > 0]
    τs_ckg = [r.tau for r in rows_ckg if isfinite(r.tau) && r.tau > 0]
    if !isempty(ns_ckg)
        Plots.plot!(plt, ns_ckg, τs_ckg;
            color = col,
            linestyle = :solid,
            linewidth = 2,
            marker = MARKER_CKG,
            markersize = 6,
            markerstrokecolor = col,
            label = @sprintf("CKG sM, β_phys=%g", β))
    end

    rows_dll = results[(:DLL, β)]
    ns_dll = [r.n for r in rows_dll if isfinite(r.tau) && r.tau > 0]
    τs_dll = [r.tau for r in rows_dll if isfinite(r.tau) && r.tau > 0]
    if !isempty(ns_dll)
        Plots.plot!(plt, ns_dll, τs_dll;
            color = col,
            linestyle = :solid,
            linewidth = 2,
            marker = MARKER_DLL,
            markersize = 6,
            markerstrokecolor = col,
            label = @sprintf("DLL M, β_phys=%g", β))
    end
end

# Plot the gap-upper-bound overlay (dashed, hollow markers, no legend entry).
for β in BETA_PHYS
    col = COLOR_BETA[β]

    rows_ckg = results[(:CKG, β)]
    ns_ckg_b = [r.n for r in rows_ckg if isfinite(r.bound) && r.bound > 0]
    bs_ckg = [r.bound for r in rows_ckg if isfinite(r.bound) && r.bound > 0]
    if !isempty(ns_ckg_b)
        Plots.plot!(plt, ns_ckg_b, bs_ckg;
            color = col,
            linestyle = :dash,
            linewidth = 1.5,
            marker = MARKER_CKG,
            markersize = 5,
            markercolor = :white,
            markerstrokecolor = col,
            markerstrokewidth = 1.5,
            label = "")
    end

    rows_dll = results[(:DLL, β)]
    ns_dll_b = [r.n for r in rows_dll if isfinite(r.bound) && r.bound > 0]
    bs_dll = [r.bound for r in rows_dll if isfinite(r.bound) && r.bound > 0]
    if !isempty(ns_dll_b)
        Plots.plot!(plt, ns_dll_b, bs_dll;
            color = col,
            linestyle = :dash,
            linewidth = 1.5,
            marker = MARKER_DLL,
            markersize = 5,
            markercolor = :white,
            markerstrokecolor = col,
            markerstrokewidth = 1.5,
            label = "")
    end
end

# Proxy legend entry for the dashed gap-bound line style.
Plots.plot!(plt, Float64[], Float64[];
    color = :black, linestyle = :dash, linewidth = 1.2,
    marker = :circle, markersize = 5,
    markercolor = :white, markerstrokecolor = :black,
    label = raw"$\log(d/\varepsilon)/\lambda_2$")

Plots.title!(plt,
    @sprintf("Ideal-Lindbladian τ_mix · 1D XXX (zz-disorder, find_typical) · ε = %.0e · s = 0.25", TARGET_EPS);
    titlefontsize = 10)

Plots.savefig(plt, FIG_PNG)
Plots.savefig(plt, FIG_PDF)
@info "Saved figure" FIG_PNG FIG_PDF

println("\n[done] wrote ", FIG_PNG)
println("[done] wrote ", FIG_PDF)
