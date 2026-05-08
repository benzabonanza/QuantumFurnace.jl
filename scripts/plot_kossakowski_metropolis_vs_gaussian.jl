#!/usr/bin/env julia
#
# Side-by-side DLL Metropolis vs DLL Gaussian Kossakowski matrix comparison
# (Ding–Li–Lin 2024, Sec. 3.2 Fig. 1 + Sec. 4 motivation).
#
# Strategy: build a 2×3 panel:
#   row 1: DLL Gaussian filter, β ∈ {1, 5, 10}
#   row 2: DLL Metropolis filter (S=2), β ∈ {1, 5, 10}
#
# Each panel shows |α(ν, ν')| with the thesis WARM gradient (cream → mulberry).
# Both Kossakowski matrices are rank-1 outer products of `freq_kernel`, so the
# panel structure mirrors the *shape* of `f̂` along the diagonal.
#
# Headline qualitative finding (the whole point of the Metropolis filter):
#   • Gaussian (Eq. 3.21–3.22): support of f̂ shrinks ∝ 1/β around -1/β.
#     At β=10, |α| is concentrated in a tiny ~0.2 × 0.2 region near (-0.1, -0.1).
#   • Metropolis (Eq. 3.19–3.20): support of f̂ stays O(1) — the flat-top
#     region [-S/2, S/2] gives |α| ≈ 1 for ν, ν' ≤ 0 even at β=10.
#
# This is THE motivating qualitative contrast for swapping Gaussian → Metropolis.
#
# System: n=3 disordered Heisenberg fixture (heis_disordered_periodic_n3.bson),
# rescaled to ‖H‖ ≤ 0.45 so all Bohr frequencies fit inside the flat top
# |ν| ≤ S/2 = 1 (S=2 default).
#
# PHYSICS CHECK: S = 2 default. The fixture has max|ν_BH| ≤ 0.9, so the
# Hörmander bump is invisible to the Lindbladian — the Metropolis filter
# acts as the bare smoothed-Metropolis weight `min(1, e^{-βν/2})` on the
# entire Bohr grid. Verified by `validate_config!` (qf-wmg.4).
#
# Output: stdout summary + (if Plots.jl is loaded) a PNG file
# `drafts/plots/kossakowski_metropolis_vs_gaussian.png`.
#
# Usage: julia --project scripts/plot_kossakowski_metropolis_vs_gaussian.jl

using Printf
using LinearAlgebra
using BSON
using QuantumFurnace

# Reuse the test-helper loader which handles the legacy BSON struct schema.
include(joinpath(@__DIR__, "..", "test", "test_helpers.jl"))

# ── Parameters ────────────────────────────────────────────────────────────────
const β_values = (1.0, 5.0, 10.0)
const S_meta   = 2.0  # Metropolis bump radius

# ── Build n=3 fixture ────────────────────────────────────────────────────────
script_dir = @__DIR__
ham_path = joinpath(script_dir, "..", "hamiltonians", "heis_disordered_periodic_n3.bson")
@info "Loading n=3 disordered Heisenberg fixture" ham_path
ham = _load_test_hamiltonian(ham_path, first(β_values))
@info "Hamiltonian (rescaled): max|eigval| = $(maximum(abs.(ham.eigvals)))"

unique_freqs = sort(collect(keys(ham.bohr_dict)))
K = length(unique_freqs)
@printf("Unique Bohr frequencies: K = %d, range [%.4f, %.4f]\n",
        K, first(unique_freqs), last(unique_freqs))

# ── Kossakowski matrix loop ──────────────────────────────────────────────────
println("\n", "="^78)
println("Kossakowski matrix sizes (Frobenius norm and entrywise extrema)")
println("="^78)
@printf("%-6s %-12s %-14s %-14s %-14s %-14s\n",
        "β", "filter", "‖α‖_F", "max|α|", "min|α(diag)|", "Σ|α|")
println("-"^78)

results = []  # cache for plotting
for β in β_values
    # Gaussian (Eq. 3.21–3.22)
    f_gauss = DLLGaussianFilter(β)
    α_gauss = dll_kossakowski_bohr(f_gauss, unique_freqs)

    # Metropolis (Eq. 3.19–3.20)
    f_meta = DLLMetropolisFilter(β; S = S_meta)
    α_meta = dll_kossakowski_bohr(f_meta, unique_freqs)

    @printf("%-6.1f %-12s %-14.4e %-14.4e %-14.4e %-14.4e\n",
            β, "Gaussian",
            norm(α_gauss), maximum(abs.(α_gauss)),
            minimum(abs.(diag(α_gauss))), sum(abs.(α_gauss)))
    @printf("%-6.1f %-12s %-14.4e %-14.4e %-14.4e %-14.4e\n",
            β, "Metropolis",
            norm(α_meta), maximum(abs.(α_meta)),
            minimum(abs.(diag(α_meta))), sum(abs.(α_meta)))
    println("-"^78)

    push!(results, (; β, α_gauss, α_meta))
end

# ── Headline qualitative summary at β = 10 ───────────────────────────────────
println("\n", "="^78)
println("Headline qualitative contrast at β = 10:")
println("="^78)
r10 = results[findfirst(r -> r.β == 10.0, results)]
α_g, α_m = r10.α_gauss, r10.α_meta

# Diagonal entries: |α(ν, ν)| = |f̂(ν)|², the cleanest probe of filter shape.
println("Diagonal |α(ν, ν)| = |f̂(ν)|² at β = 10:")
@printf("%-10s %-18s %-18s %-12s\n", "ν", "Gaussian", "Metropolis", "ratio M/G")
for k in 1:K
    ν = unique_freqs[k]
    g_diag = abs(α_g[k, k])
    m_diag = abs(α_m[k, k])
    rt = g_diag > 0 ? m_diag / g_diag : Inf
    @printf("%-10.3f %-18.4e %-18.4e %-12.2e\n", ν, g_diag, m_diag, rt)
end
println("-"^78)
println("Gaussian filter centred at -1/β = -0.1 with width ~ 2/β = 0.2:")
println("  → |α_gauss(ν, ν)| collapses for |ν| ≫ 1/β.")
println("Metropolis filter saturates at f̂(ν) → 1 for ν ≪ 0 on the flat top:")
println("  → |α_meta(ν, ν)| ≈ 1 across the entire negative half of the Bohr grid.")
println("  → THIS is the qualitative win: low-T Metropolis preserves O(1) coupling")
println("    weights, while Gaussian drives them to numerical zero.")
println("="^78)

# ── Optional plot output ─────────────────────────────────────────────────────
global plot_output = nothing
try
    @eval using Plots
    global plot_output = joinpath(@__DIR__, "..", "drafts", "plots",
                                  "kossakowski_metropolis_vs_gaussian.png")
    mkpath(dirname(plot_output))
catch
    @info "Plots.jl not available; skipping figure generation."
end

if plot_output !== nothing
    # Thesis warm gradient: cream → deep mulberry.
    WARM_GRAD = Plots.cgrad(["#FCE1A4", "#FABF7B", "#F08F6E", "#E05C5C",
                             "#D12959", "#AB1866", "#6E005F"])

    plots = []
    # Common log-scale color limits across all panels for fair comparison.
    all_vals = vcat([abs.(r.α_gauss) for r in results]...,
                    [abs.(r.α_meta) for r in results]...)
    floor_v = max(minimum(all_vals[all_vals .> 0]), 1e-10)
    ceil_v  = maximum(all_vals)

    # Layout: row 1 Gaussian, row 2 Metropolis, columns are β.
    for (i, r) in enumerate(results)
        # Gaussian panel
        z_g = max.(abs.(r.α_gauss), floor_v)  # clip for log10
        p_g = Plots.heatmap(unique_freqs, unique_freqs, log10.(z_g);
            title = "Gaussian, β=$(r.β)",
            xlabel = "ν'", ylabel = i == 1 ? "ν" : "",
            c = WARM_GRAD,
            clim = (log10(floor_v), log10(ceil_v)),
            aspect_ratio = :equal,
            colorbar = i == 3,
            colorbar_title = i == 3 ? "log₁₀ |α|" : "")
        push!(plots, p_g)
    end
    for (i, r) in enumerate(results)
        # Metropolis panel
        z_m = max.(abs.(r.α_meta), floor_v)
        p_m = Plots.heatmap(unique_freqs, unique_freqs, log10.(z_m);
            title = "Metropolis (S=$(S_meta)), β=$(r.β)",
            xlabel = "ν'", ylabel = i == 1 ? "ν" : "",
            c = WARM_GRAD,
            clim = (log10(floor_v), log10(ceil_v)),
            aspect_ratio = :equal,
            colorbar = i == 3,
            colorbar_title = i == 3 ? "log₁₀ |α|" : "")
        push!(plots, p_m)
    end

    fig = Plots.plot(plots...; layout = (2, 3), size = (1500, 900))
    Plots.savefig(fig, plot_output)
    @info "Saved Kossakowski heatmap" plot_output
end
