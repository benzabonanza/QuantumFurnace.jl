#!/usr/bin/env julia
#
# Plot P1 — ideal-Lindbladian τ_mix: DLL Metropolis vs CKG smooth-Metropolis.
# Single panel at ε=1e-3 with the spectral-gap upper-bound overlay
# τ_bound = log(d/ε) / λ_2 (per the 2026-05-08 spec correction in qf-e4z.11).
#
# Source data (no new simulations — pure post-processing):
#   S1 (qf-e4z.3, re-run v3 on 2026-05-10):
#     scripts/output/sweep_S1_v3_ckg_ideal/smooth_metro_eps1e-03/sweep_n*_beta*_seed42_L_KMS_Energy.bson
#   S2 (qf-e4z.4, re-run v3 on 2026-05-10):
#     scripts/output/sweep_S2_v3_dll_ideal/smooth_metro_eps1e-03/sweep_n*_beta*_seed42_L_DLL_Bohr.bson
#
# Output:
#   drafts/figures/numerics/p1_taumix.{png, pdf}
#
# Layout (per qf-e4z.11 corrected spec):
#   - x = n (3..8 linear)
#   - y = τ_mix (log10)
#   - 6 simulated curves: {CKG smooth-Metro, DLL Metro} × {β = 5, 10, 20}
#   - 6 gap-bound overlays τ_bound = log(d/ε)/λ_2 (dashed, hollow markers)
#   - Color = β (terracotta=5 → bordeaux=10 → aubergine=20)
#   - Marker = filter (○ = CKG smooth-Metro, □ = DLL Metropolis)
#
# Usage: JULIA_NUM_THREADS=1 julia --project scripts/numerics_p1_plot.jl

using Printf
using BSON
using QuantumFurnace                                # needed for BSON symbol resolution
ENV["GKSwstype"] = "100"
using Plots

# ── Paths ─────────────────────────────────────────────────────────────────────
const REPO_ROOT = joinpath(@__DIR__, "..")
const S1_DIR    = joinpath(REPO_ROOT, "scripts", "output", "sweep_S1_v3_ckg_ideal",
                            "smooth_metro_eps1e-03")
const S2_DIR    = joinpath(REPO_ROOT, "scripts", "output", "sweep_S2_v3_dll_ideal",
                            "smooth_metro_eps1e-03")
const OUT_DIR   = joinpath(REPO_ROOT, "drafts", "figures", "numerics")
const FIG_PNG   = joinpath(OUT_DIR, "p1_taumix.png")
const FIG_PDF   = joinpath(OUT_DIR, "p1_taumix.pdf")
mkpath(OUT_DIR)

# ── Sweep grid (matches v3 sweeps) ────────────────────────────────────────────
const N_VALUES    = 3:8
const BETA_VALUES = (5.0, 10.0, 20.0)
const TARGET_EPS  = 1e-3
const SEED        = 42

# ── Thesis colour palette (memory: reference_thesis_colors.md) ───────────────
# β-gradient (warm→dark): low β = warm, high β = aubergine.
const COLOR_BETA = Dict(
    5.0  => "#B5654A",   # terracotta
    10.0 => "#7A2E39",   # bordeaux
    20.0 => "#4F3B5C",   # aubergine
)

# Filter glyphs (one shape per construction so colour can encode β cleanly).
const MARKER_CKG = :circle
const MARKER_DLL = :rect

# ── Schema-aware BSON loader ──────────────────────────────────────────────────
# v3 sidecars store the per-cell record under :result (a Dict).
function _beta_str(β::Real)
    s = @sprintf("%.6f", float(β))
    s = rstrip(s, '0'); s = rstrip(s, '.')
    return isempty(s) ? "0" : s
end

function _load_cell(dir::AbstractString, n::Integer, β::Real,
                    construction_tag::AbstractString, domain_tag::AbstractString)
    path = joinpath(dir, "sweep_n$(n)_beta$(_beta_str(β))_seed$(SEED)_L_$(construction_tag)_$(domain_tag).bson")
    isfile(path) || return nothing
    data = BSON.load(path, QuantumFurnace)
    return data[:result]
end

# Pull τ_mix and τ_bound from a cell dict (returns NaN if either is missing).
function _tau_and_bound(cell)
    cell === nothing && return (NaN, NaN)
    τ = get(cell, :mixing_time, NaN)
    τ_b = get(cell, :tau_mix_bound, NaN)
    return (Float64(τ), Float64(τ_b))
end

# ── Collect data ──────────────────────────────────────────────────────────────
println("="^72)
println("P1 plot — DLL Metropolis vs CKG smooth-Metro (ideal Lindbladian)")
println("S1 dir: ", S1_DIR)
println("S2 dir: ", S2_DIR)
println("="^72)

results = Dict{Tuple{Symbol,Float64},Vector{NamedTuple{(:n,:tau,:bound),Tuple{Int,Float64,Float64}}}}()
for β in BETA_VALUES
    results[(:CKG, β)] = NamedTuple{(:n,:tau,:bound),Tuple{Int,Float64,Float64}}[]
    results[(:DLL, β)] = NamedTuple{(:n,:tau,:bound),Tuple{Int,Float64,Float64}}[]
end

for n in N_VALUES, β in BETA_VALUES
    ckg_cell = _load_cell(S1_DIR, n, β, "KMS", "Energy")
    dll_cell = _load_cell(S2_DIR, n, β, "DLL", "Bohr")
    τ_ckg, b_ckg = _tau_and_bound(ckg_cell)
    τ_dll, b_dll = _tau_and_bound(dll_cell)
    push!(results[(:CKG, β)], (n=n, tau=τ_ckg, bound=b_ckg))
    push!(results[(:DLL, β)], (n=n, tau=τ_dll, bound=b_dll))
    @printf("n=%d  β=%-4g | CKG τ=%.3e  bound=%.3e | DLL τ=%.3e  bound=%.3e | ratio τ DLL/CKG = %.2f\n",
            n, β, τ_ckg, b_ckg, τ_dll, b_dll,
            (isfinite(τ_ckg) && τ_ckg > 0) ? τ_dll/τ_ckg : NaN)
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
for β in BETA_VALUES
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
            label = @sprintf("CKG sM, β=%g", β))
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
            label = @sprintf("DLL M, β=%g", β))
    end
end

# Plot the gap-upper-bound overlay (dashed, hollow markers, no legend entry).
for β in BETA_VALUES
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

# Add a single proxy legend entry for the dashed gap-bound line style.
Plots.plot!(plt, Float64[], Float64[];
    color = :black, linestyle = :dash, linewidth = 1.2,
    marker = :circle, markersize = 5,
    markercolor = :white, markerstrokecolor = :black,
    label = raw"$\log(d/\varepsilon)/\lambda_2$")

Plots.title!(plt,
    @sprintf("Ideal-Lindbladian τ_mix · 1D XXX (zz-disorder) · ε = %.0e", TARGET_EPS);
    titlefontsize = 11)

Plots.savefig(plt, FIG_PNG)
Plots.savefig(plt, FIG_PDF)
@info "Saved figure" FIG_PNG FIG_PDF

println("\n[done] wrote ", FIG_PNG)
println("[done] wrote ", FIG_PDF)
