#!/usr/bin/env julia
#
# qf-7go.7 — Plot τ_mix vs k for the multi-rank DLL sweep.
#
# 2-panel figure: left for DLL Metropolis base, right for DLL Gaussian base.
# Each panel shows τ_mix(k) at β ∈ {1, 5, 10, 20} on a log-y axis. The
# CKG smooth-Metropolis baseline (loaded from the qf-mto / fair-comparison
# BSON when available) is drawn as a horizontal reference per β.
#
# Inputs:
#   scripts/output/dll_multirank/taumix_vs_k.bson   (qf-7go.6 sweep)
#   scripts/output/fair_comparison_kms_dirichlet.bson (CKG baseline; optional)
#
# Outputs:
#   drafts/figures/numerics/dll_multirank_taumix.png
#   drafts/figures/numerics/dll_multirank_taumix.pdf
#
# Run:
#   julia --project scripts/plot_dll_multirank_taumix.jl
#

using Printf
using BSON
ENV["GKSwstype"] = "100"
using Plots

const ROOT = dirname(@__DIR__)
const SWEEP_BSON   = joinpath(ROOT, "scripts", "output", "dll_multirank", "taumix_vs_k.bson")
const FAIR4_BSON   = joinpath(ROOT, "scripts", "output", "fair_comparison_kms_dirichlet.bson")
const OUT_DIR      = joinpath(ROOT, "drafts", "figures", "numerics")
const OUT_PNG      = joinpath(OUT_DIR, "dll_multirank_taumix.png")
const OUT_PDF      = joinpath(OUT_DIR, "dll_multirank_taumix.pdf")

isfile(SWEEP_BSON) || error("Missing sweep BSON: $SWEEP_BSON")
mkpath(OUT_DIR)

raw = BSON.load(SWEEP_BSON)
const RESULTS    = raw[:results]
const BETAS      = Float64.(raw[:beta_values])
const KS         = Int.(raw[:k_values])
const SHAPES     = Symbol.(raw[:base_shapes])
const N_QUBITS   = Int(raw[:n])

@printf("Loaded multi-rank sweep: n=%d, %d cells.\n", N_QUBITS, length(RESULTS))

# ── Optional CKG smooth-Metropolis τ_mix reference (from fair-comparison) ────
ckg_τmix = Dict{Float64, Float64}()
if isfile(FAIR4_BSON)
    raw4 = BSON.load(FAIR4_BSON)
    fair4_results = raw4[:results]
    for r in fair4_results
        if Symbol(r.sampler_key) === :ckg_smooth_metro && Int(r.n) == N_QUBITS
            β = Float64(r.β)
            # qf-mto records the dimensionless intrinsic mixing ratio
            # ρ_intrinsic = λ/Λ_max; convert via τ ≈ 2/λ. For a CKG-vs-DLL
            # τ_mix overlay we want absolute τ_mix; if the BSON also stores
            # mixing_time use that, otherwise leave the dict empty.
            if hasproperty(r, :mixing_time) && isfinite(r.mixing_time)
                ckg_τmix[β] = Float64(r.mixing_time)
            end
        end
    end
    @printf("Loaded %d CKG sM τ_mix references from %s.\n",
            length(ckg_τmix), FAIR4_BSON)
end

# ── Helpers ──────────────────────────────────────────────────────────────────
function _τmix_table(shape::Symbol)
    # Returns a Dict β → Vector{τmix} of length length(KS).
    out = Dict{Float64, Vector{Float64}}()
    for β in BETAS
        row = Float64[]
        for k in KS
            idx = findfirst(r -> r.beta == β && r.k == k && r.base_shape == shape,
                            RESULTS)
            push!(row, idx === nothing ? NaN : Float64(RESULTS[idx].mixing_time))
        end
        out[β] = row
    end
    return out
end

# Thesis palette (β-coloured lines).
const COLOR_BETA = Dict(
    1.0  => "#bdd7e7",   # pale blue
    5.0  => "#6baed6",   # blue
    10.0 => "#3182bd",   # deep blue
    20.0 => "#08519c",   # navy
)
const COLOR_CKG = "#d62728"   # tomato red (reference baselines)

plot_metro = let
    table = _τmix_table(:metropolis)
    plt = plot(;
        xlabel = "k (number of channels)",
        ylabel = "τ_mix",
        title = "DLL Metropolis (n = $N_QUBITS)",
        yscale = :log10,
        xticks = (KS, string.(KS)),
        legend = :topright,
        framestyle = :box,
        gridalpha = 0.35,
    )
    for β in BETAS
        ys = table[β]
        plot!(plt, KS, ys;
              label = "β = $(Int(β))",
              color = COLOR_BETA[β],
              lw = 2.4,
              marker = :circle,
              markersize = 5,
              markerstrokewidth = 0)
        if haskey(ckg_τmix, β)
            hline!(plt, [ckg_τmix[β]];
                   color = COLOR_BETA[β],
                   linestyle = :dot,
                   lw = 1.4,
                   label = "")
        end
    end
    plt
end

plot_gauss = let
    table = _τmix_table(:gaussian)
    plt = plot(;
        xlabel = "k (number of channels)",
        ylabel = "τ_mix",
        title = "DLL Gaussian (n = $N_QUBITS)",
        yscale = :log10,
        xticks = (KS, string.(KS)),
        legend = :topright,
        framestyle = :box,
        gridalpha = 0.35,
    )
    for β in BETAS
        ys = table[β]
        plot!(plt, KS, ys;
              label = "β = $(Int(β))",
              color = COLOR_BETA[β],
              lw = 2.4,
              marker = :circle,
              markersize = 5,
              markerstrokewidth = 0)
        if haskey(ckg_τmix, β)
            hline!(plt, [ckg_τmix[β]];
                   color = COLOR_BETA[β],
                   linestyle = :dot,
                   lw = 1.4,
                   label = "")
        end
    end
    plt
end

fig = plot(plot_metro, plot_gauss;
           layout = (1, 2),
           size = (1000, 420),
           left_margin = 4Plots.mm,
           bottom_margin = 4Plots.mm,
           top_margin = 2Plots.mm,
           right_margin = 2Plots.mm)

savefig(fig, OUT_PNG)
savefig(fig, OUT_PDF)
@printf("Wrote %s\n", OUT_PNG)
@printf("Wrote %s\n", OUT_PDF)

# ── Numerical summary ───────────────────────────────────────────────────────
println("\nMulti-rank speedup (τ_mix at k=1) / (τ_mix at k=k_max):")
for shape in SHAPES
    table = _τmix_table(shape)
    @printf("  [%-10s] β  k=1→k=%d ratio\n", string(shape), maximum(KS))
    for β in BETAS
        ys = table[β]
        ratio = ys[1] / ys[end]
        @printf("    β=%-5.1f  %5.2f×  (%.3g → %.3g)\n",
                β, ratio, ys[1], ys[end])
    end
end
