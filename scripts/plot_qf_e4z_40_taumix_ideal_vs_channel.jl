#!/usr/bin/env julia
# plot_qf_e4z_40_taumix_ideal_vs_channel.jl  (qf-e4z.40 wrap-up)
#
# τ_mix^phys vs n, for each β_phys ∈ {0.25, 0.5, 1.0}, comparing:
#   • IDEAL Lindbladian L_Energy  (qf-e4z.34 sidecars)
#   • CHANNEL Φ_δ Trotter         (qf-e4z.39 v2 Pass-1-patched sidecars)
#
# Single-seed = 46, ρ_0 = |+⟩⟨+|^⊗N, ε_target = 1e-3 — the canonical setup
# for this thesis comparison. All 18 cells (n=3..8 × β∈{0.25,0.5,1.0}) are
# :extrapolated (no floor-bisection artefacts).
#
# Plot convention:
#   • x-axis: n (system size)
#   • y-axis: τ_mix^phys (β_phys frame)
#   • Color = β_phys (uses the COLOR_BETA palette shared across thesis figs)
#   • Marker / line:  ideal = filled circle + solid line
#                     channel = open square + dashed line
#
# Output:
#   drafts/figures/numerics/qf_e4z_40_taumix_ideal_vs_channel.{png,pdf}

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BSON
using Printf
ENV["GKSwstype"] = "100"
using Plots

const REPO_ROOT = joinpath(@__DIR__, "..")
const ID_DIR    = joinpath(REPO_ROOT, "scripts", "output", "sweep_qf_e4z_34_ckg_vs_dll_plot/ckg")
const CH_DIR    = joinpath(REPO_ROOT, "scripts", "output", "sweep_qf_e4z_39_channel_ckg_seed46_v2")
const OUT_DIR   = joinpath(REPO_ROOT, "drafts", "figures", "numerics")
mkpath(OUT_DIR)

const SEED         = 46
const N_RANGE      = 3:8
const BETA_PHYS_ALL = (0.25, 0.5, 1.0)

# Shared thesis palette (matches plot_qf_e4z_34_ckg_spectral_bound.jl)
const COLOR_BETA = Dict(
    0.25 => "#B5654A",
    0.5  => "#7A2E39",
    1.0  => "#4F3B5C",
)

function _bp_str(β)
    return β == 0.25 ? "0.25" : (β == 0.5 ? "0.5" : "1")
end

function _load_ideal(n, β)
    path = joinpath(ID_DIR, "sweep_n$(n)_betaphys$(_bp_str(β))_seed$(SEED)_L_KMS_Energy.bson")
    isfile(path) || return nothing
    r = BSON.load(path)[:result]
    src = get(r, :tau_mix_source, get(r, :mixing_time_source, :unknown))
    τ   = get(r, :mixing_time_phys, get(r, :tau_mix_phys, NaN))
    return (τ = τ, src = src)
end

function _load_channel(n, β)
    path = joinpath(CH_DIR,
        "channel_n$(n)_betaphys$(_bp_str(β))_seed$(SEED)_eps1e-03_smooth_metro_KMS_Trotter.bson")
    isfile(path) || return nothing
    r = BSON.load(path)[:result]
    src = get(r, :tau_mix_channel_source, :unknown)
    τ   = get(r, :mixing_time_phys_channel, NaN)
    return (τ = τ, src = src)
end

function main()
    # Collect: ns[β] -> Vector{Int} ; ideal[β] -> Vector{Float64}; channel[β] -> Vector{Float64}
    ns       = Dict{Float64, Vector{Int}}()
    τ_id_per = Dict{Float64, Vector{Float64}}()
    τ_ch_per = Dict{Float64, Vector{Float64}}()

    for β in BETA_PHYS_ALL
        ns[β] = Int[]; τ_id_per[β] = Float64[]; τ_ch_per[β] = Float64[]
        for n in N_RANGE
            id = _load_ideal(n, β)
            ch = _load_channel(n, β)
            (id === nothing || ch === nothing) && continue
            id.src === :extrapolated || @warn "ideal not :extrapolated" n β src=id.src
            ch.src === :extrapolated || @warn "channel not :extrapolated" n β src=ch.src
            push!(ns[β], n); push!(τ_id_per[β], id.τ); push!(τ_ch_per[β], ch.τ)
        end
    end

    # ------ Plot --------------------------------------------------------
    plt = plot(
        size            = (700, 460),
        xlabel          = "system size n",
        ylabel          = "τ_mix^phys",
        title           = "τ_mix: ideal Lindbladian vs channel Φ_δ (β_phys frame, seed=$SEED, ρ_0 = |+⟩^⊗N)",
        legend          = :topleft,
        xticks          = collect(N_RANGE),
        gridstyle       = :dash,
        gridalpha       = 0.25,
        framestyle      = :box,
        titlefontsize   = 9,
        guidefontsize   = 10,
        tickfontsize    = 9,
        legendfontsize  = 8,
        bottom_margin   = 3Plots.mm,
        left_margin     = 3Plots.mm,
    )

    for β in BETA_PHYS_ALL
        c = COLOR_BETA[β]
        plot!(plt, ns[β], τ_id_per[β];
            seriestype = :path,
            color      = c, lw = 2.0, linestyle = :solid,
            marker     = :circle, ms = 6, msw = 0.5, mc = c,
            label      = "ideal L,    β_phys = $β",
        )
        plot!(plt, ns[β], τ_ch_per[β];
            seriestype = :path,
            color      = c, lw = 1.8, linestyle = :dash,
            marker     = :square, ms = 6, msw = 1.5, mc = :white, msc = c,
            label      = "channel Φ_δ, β_phys = $β",
        )
    end

    out_png = joinpath(OUT_DIR, "qf_e4z_40_taumix_ideal_vs_channel.png")
    out_pdf = joinpath(OUT_DIR, "qf_e4z_40_taumix_ideal_vs_channel.pdf")
    savefig(plt, out_png)
    savefig(plt, out_pdf)
    @printf("[saved] %s\n", out_png)
    @printf("[saved] %s\n", out_pdf)

    # ------ Per-cell summary table -------------------------------------
    println("\n" * "="^72)
    println("τ_mix^phys summary (β_phys frame)")
    println("="^72)
    @printf("%-3s %-5s | %-9s %-9s %-7s\n", "n", "β", "ideal", "channel", "ratio")
    println("-"^72)
    for β in BETA_PHYS_ALL
        for (i, n) in enumerate(ns[β])
            τi, τc = τ_id_per[β][i], τ_ch_per[β][i]
            r = τc / τi
            mark = abs(r - 1) > 0.15 ? " ⚠" : ""
            @printf("%-3d %-5g | %-9.4f %-9.4f %-7.4f%s\n", n, β, τi, τc, r, mark)
        end
        println()
    end
    println("[done]")
end

main()
