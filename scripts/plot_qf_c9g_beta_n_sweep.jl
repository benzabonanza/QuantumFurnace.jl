#!/usr/bin/env julia
# plot_qf_c9g_beta_n_sweep.jl  (qf-c9g)
#
# Render the n × β gap-closing diagnostic from
#   scripts/output/qf_c9g_ordered_gap_mechanism/qf_c9g_beta_n_sweep.bson
# Produces two figures:
#   drafts/figures/numerics/qf_c9g_gap_phys_vs_n.{png,pdf}
#       gap_phys vs n, one curve per β_phys (log y-axis). Companion:
#       (a) doublet split ΔE_1^phys(n), Hamiltonian-only reference;
#       (b) λ_3/λ_2 ratio (doublet-mode capture).
#   drafts/figures/numerics/qf_c9g_phase_indicators.{png,pdf}
#       Binder U_4 + ⟨m²⟩ + eff_rank vs β at each n.
#
# Run:  julia --project scripts/plot_qf_c9g_beta_n_sweep.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BSON
using Plots
using Printf

const BSON_PATH = joinpath(@__DIR__, "output", "qf_c9g_ordered_gap_mechanism",
                            "qf_c9g_beta_n_sweep.bson")
const FIG_DIR = joinpath(@__DIR__, "..", "drafts", "figures", "numerics")
mkpath(FIG_DIR)

# T_c at h=1
const T_C = 2.07

function load_data()
    raw = BSON.load(BSON_PATH)
    sweep_a = [NamedTuple{Tuple(Symbol(k) for k in keys(d))}(values(d)) for d in raw[:sweep_a_rows]]
    sweep_b = [NamedTuple{Tuple(Symbol(k) for k in keys(d))}(values(d)) for d in raw[:sweep_b_rows]]
    return (
        sweep_a = sweep_a,
        sweep_b = sweep_b,
        n_list = raw[:n_list],
        beta_grid = raw[:beta_grid],
    )
end

function plot_gap_vs_n(data)
    sweep_a = data.sweep_a
    n_list = data.n_list
    beta_grid = data.beta_grid

    # Three panels: λ_L vs n (log y), ΔE_1^phys vs n, λ_3/λ_2 vs n
    p1 = plot(xlabel = "n (sites)",
              ylabel = "Lindbladian gap_phys",
              yscale = :log10,
              legend = :topright,
              size = (700, 420),
              title = "(a) λ_L^phys(n) at each β_phys — closing only below T_c")
    p2 = plot(xlabel = "n (sites)",
              ylabel = "ΔE_1^phys (Hamiltonian doublet split)",
              yscale = :log10,
              legend = false,
              size = (700, 420),
              title = "(b) Hamiltonian doublet split (β-independent)")
    p3 = plot(xlabel = "n (sites)",
              ylabel = "|λ_3 / λ_2|",
              yscale = :log10,
              legend = :topleft,
              size = (700, 420),
              title = "(c) doublet vs bulk: |λ_3 / λ_2|")

    # Color per β (cold→hot)
    sorted_betas = sort(beta_grid)
    colors = palette(:viridis, length(sorted_betas))

    for (i, β) in enumerate(sorted_betas)
        ord = (1/β) < T_C ? "*" : ""   # mark ordered cells
        β_rows = sort([r for r in sweep_a if isapprox(r.beta_phys, β; atol=1e-12)], by = r->r.n)
        ns = [r.n for r in β_rows]
        gaps = [r.gap_phys for r in β_rows]
        ratios = [r.λ3_over_λ2 for r in β_rows]
        label = @sprintf("β=%.2g (T/T_c=%.2g)%s", β, (1/β)/T_C, ord)
        plot!(p1, ns, gaps; lw=2, marker=:circle, color=colors[i], label=label)
        plot!(p3, ns, ratios; lw=2, marker=:circle, color=colors[i], label=label)
    end

    # ΔE_1^phys is β-independent; plot from sweep_b only
    sb = sort(data.sweep_b, by = r->r.n)
    plot!(p2, [r.n for r in sb], [r.doublet_split for r in sb];
        lw=2, marker=:circle, color=:black, label="ΔE_1^phys")
    plot!(p2, [r.n for r in sb], [r.bulk_gap for r in sb];
        lw=2, marker=:square, color=:darkred, label="bulk gap (E_4 − E_1)")

    fig = plot(p1, p2, p3; layout = (3, 1), size = (700, 1260))
    out_png = joinpath(FIG_DIR, "qf_c9g_gap_phys_vs_n.png")
    out_pdf = joinpath(FIG_DIR, "qf_c9g_gap_phys_vs_n.pdf")
    savefig(fig, out_png); savefig(fig, out_pdf)
    println("[saved] $out_png")
    println("[saved] $out_pdf")
    return fig
end

function plot_phase_indicators(data)
    sweep_a = data.sweep_a
    n_list = sort(data.n_list)
    beta_grid = sort(data.beta_grid)

    p1 = plot(xlabel = "β_phys",
              ylabel = "<m_z²>",
              xscale = :log10,
              ylim = (0, 1.05),
              legend = :bottomright,
              title = "(a) magnetization-squared")
    p2 = plot(xlabel = "β_phys",
              ylabel = "Binder U_4",
              xscale = :log10,
              ylim = (0, 0.75),
              legend = :bottomright,
              title = "(b) Binder cumulant")
    p3 = plot(xlabel = "β_phys",
              ylabel = "eff_rank",
              xscale = :log10,
              yscale = :log10,
              legend = :topright,
              title = "(c) Gibbs eff_rank")

    colors = palette(:plasma, length(n_list))
    for (i, n) in enumerate(n_list)
        n_rows = sort([r for r in sweep_a if r.n == n], by = r->r.beta_phys)
        βs = [r.beta_phys for r in n_rows]
        m2s = [r.m2 for r in n_rows]
        u4s = [r.U4 for r in n_rows]
        eff_ranks = [r.eff_rank for r in n_rows]
        plot!(p1, βs, m2s; lw=2, marker=:circle, color=colors[i], label="n=$n")
        plot!(p2, βs, u4s; lw=2, marker=:circle, color=colors[i], label="n=$n")
        plot!(p3, βs, eff_ranks; lw=2, marker=:circle, color=colors[i], label="n=$n")
    end

    # mark T_c
    vline!(p1, [1/T_C]; ls=:dash, color=:gray, label="β_c = 1/T_c")
    vline!(p2, [1/T_C]; ls=:dash, color=:gray, label="β_c = 1/T_c")
    vline!(p3, [1/T_C]; ls=:dash, color=:gray, label="β_c = 1/T_c")

    # Binder reference lines
    hline!(p2, [2/3]; ls=:dot, color=:green, label="U_4 → 2/3 (ordered)")
    hline!(p2, [0.0]; ls=:dot, color=:red, label="U_4 → 0 (paramagnetic)")

    # eff_rank reference: rank 2 (doublet)
    hline!(p3, [2.0]; ls=:dot, color=:green, label="eff_rank = 2 (doublet)")

    fig = plot(p1, p2, p3; layout = (3, 1), size = (700, 1260))
    out_png = joinpath(FIG_DIR, "qf_c9g_phase_indicators.png")
    out_pdf = joinpath(FIG_DIR, "qf_c9g_phase_indicators.pdf")
    savefig(fig, out_png); savefig(fig, out_pdf)
    println("[saved] $out_png")
    println("[saved] $out_pdf")
    return fig
end

function print_table(data)
    println("\n" * "="^120)
    println("Summary table — gap_phys (Lindbladian), ΔE_1^phys (Hamiltonian), λ3/λ2")
    println("="^120)
    @printf("%-3s %-7s %-7s %-9s %-12s %-12s %-9s %-9s %-9s\n",
        "n", "β_phys", "T/T_c", "phase", "λ_L^phys", "ΔE_1^phys",
        "λ_L/ΔE_1", "λ3/λ2", "eff_rank")
    println("-"^120)
    for n in sort(data.n_list)
        for β in sort(data.beta_grid)
            row = first(r for r in data.sweep_a if r.n == n && isapprox(r.beta_phys, β; atol=1e-12))
            phase = (1/β) < T_C ? "ORD" : "DIS"
            @printf("%-3d %-7.3g %-7.3g %-9s %-12.4e %-12.4e %-9.3g %-9.3g %-9.3g\n",
                n, β, (1/β)/T_C, phase,
                row.gap_phys, row.Δ1_phys, row.ratio_gap_dE1, row.λ3_over_λ2, row.eff_rank)
        end
        println()
    end
end

function main()
    data = load_data()
    print_table(data)
    plot_gap_vs_n(data)
    plot_phase_indicators(data)
end

main()
