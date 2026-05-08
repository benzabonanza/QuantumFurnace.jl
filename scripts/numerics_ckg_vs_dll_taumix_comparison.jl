#!/usr/bin/env julia
#
# CKG vs DLL τ_mix(n, β) comparison plot for the thesis numerics chapter (qf-lkb.6).
#
# qf-lkb.11 update (2026-05-02):
#   - CKG path now runs in EnergyDomain instead of BohrDomain (matvec scales
#     as O(N · d²) with N = 2^num_energy_bits = 4096 fixed, vs the BohrDomain
#     O(d⁴) over Bohr-pair indices). This unblocks n ≥ 5 sweeps on the laptop.
#   - Smooth-Metropolis defaults switched to a=0, s=0.25 (the thesis-numerics
#     convention). Previous a=β/30, s=0.4 only differs in transition-weight
#     shape, but the new value is the locked thesis convention.
#   - sweep_cache/ckg sidecar names now embed the domain tag, so old
#     BohrDomain caches won't be reused (they were under different params
#     anyway). DLL caches keep their BohrDomain names — DLL path is unchanged.
#
# Strategy: drive the matrix-free ODE-integrator + bi-exp extrapolation pipeline
# (`integrate_to_gibbs` + `estimate_mixing_time(...; model=:biexp,
# extrapolate=true)`) over the (n, β) product for three constructions on the
# same disordered Heisenberg fixture family
# (`hamiltonians/heis_disordered_periodic_n<n>.bson`):
#
#   1) CKG (KMS, EnergyDomain) — matrix-free `apply_lindbladian!` (Workspace path)
#   2) DLL Gaussian filter      — matrix-free DLL apply (qf-lkb.9, BohrDomain)
#   3) DLL Metropolis filter    — matrix-free DLL apply (qf-lkb.9, BohrDomain)
#
# Sweep grid (3 × 4 × 3 = 36 cells):
#   n ∈ {3, 4, 5},  β ∈ {1.0, 2.0, 5.0, 10.0},  construction × filter as above.
#
# Output:
#   drafts/figures/numerics/ckg_vs_dll_taumix.{png, pdf}     (figure)
#   drafts/figures/numerics/ckg_vs_dll_taumix.bson           (raw sweep data)
#   drafts/figures/numerics/sweep_cache/{ckg, dll_gauss, dll_metro}/  (per-cell sidecars; resumable)
#
# Plot layout (2×3 subplots):
#   row 1: τ_mix vs n  with one line per β  (CKG | DLL Gaussian | DLL Metropolis)
#   row 2: τ_mix vs β  with one line per n  (CKG | DLL Gaussian | DLL Metropolis)
#   y-axis: log-scale (τ_mix spans orders of magnitude across β)
#   x-axis: linear (small range)
#   colour palettes: thesis colours (memory: reference_thesis_colors.md)
#     β series: slateblue (β=1) → dustyteal (β=2) → terracotta (β=5) → bordeaux (β=10)
#     n series: pinegreen (n=3) → sage (n=4) → ochre (n=5)
#
# PHYSICS CHECK: β min bumped from the qf-lkb.6 audit's 0.5 → 1.0 because at
# β = 0.5 the system is essentially infinite-temperature; τ_mix is dominated
# by H rotation rather than thermal relaxation, so the CKG-vs-DLL contrast
# is uninteresting at β < 1.
#
# PHYSICS CHECK: n cap at 5 for thesis-scale visualisation (n=6,7 are feasible
# on this laptop in EnergyDomain CKG and matrix-free DLL after qf-lkb.9, but
# require fresh Hamiltonian fixtures). Cluster sweeps can extend to n ≥ 7.
#
# PHYSICS CHECK: target_epsilon = 1e-3, t_max_factor = 5.0. The qf-lkb.3 testset
# settled on 5.0 / gap_estimate as a robust horizon for bi-exp extrapolation
# (single-exp regime starts ~3/gap, bi-exp wants ~3× past crossing).
#
# Usage:  JULIA_NUM_THREADS=4 julia --project scripts/numerics_ckg_vs_dll_taumix_comparison.jl
#
# Notes:
# - skip_existing = true ⇒ re-runs of this script reuse already-finished cells
#   from the per-construction sweep_cache subdirectories.
# - Plots.jl is in the project's [extras]/test target and the Project.toml
#   `[compat]` block. It is *not* a runtime [deps] entry, but `using Plots`
#   succeeds when the user has it precompiled into the project environment.
#   If the script aborts on `using Plots`, run
#       julia --project -e 'using Pkg; Pkg.add("Plots")'
#   first, or rely on Plots being inherited from the default global env.

using Printf
using Random
using LinearAlgebra
using BSON
using QuantumFurnace
# Force the GR backend's offscreen PNG-rendering driver before `using Plots`,
# otherwise PNG `savefig` fails with `GKS: cannot open display` and writes a
# 0-byte file in headless environments. The PDF path uses Cairo and is unaffected.
ENV["GKSwstype"] = "100"
using Plots

# ── Sweep grid ────────────────────────────────────────────────────────────────
# n_values = {3, 4, 5}: n=5 is now affordable in EnergyDomain CKG (qf-lkb.11)
# and matrix-free DLL (qf-lkb.9). The legacy heis_disordered_periodic_n*.bson
# fixture family covers exactly n ∈ {3, 4, 5}; n=6,7 extension requires
# generating fresh fixtures or switching to heis_xxx_zzdisordered_periodic_n*.
const N_VALUES    = [3, 4, 5]
const BETA_VALUES = [1.0, 2.0, 5.0, 10.0]
const TARGET_EPS  = 1e-3
const T_MAX_FACTOR = 5.0      # PHYSICS CHECK: per qf-lkb.3 testset finding
const T_GRID_LEN   = 81       # bi-exp fitting wants ≥ 50 well-spaced samples
const KRYLOV_DIM   = 30
const TOL          = 1e-10
const SEEDS        = [42]      # single seed; this script targets the figure, not seed-stats

# ── Thesis colour palette (memory: reference_thesis_colors.md) ───────────────
const COLOR_BETA = Dict(
    1.0  => "#5C7794",   # slateblue
    2.0  => "#5F8B8E",   # dustyteal
    5.0  => "#B5654A",   # terracotta
    10.0 => "#7A2E39",   # bordeaux
)
const COLOR_N = Dict(
    3 => "#2D5A3D",      # pinegreen
    4 => "#8B9F7E",      # sage
    5 => "#B89143",      # ochre
)

# ── Output paths ──────────────────────────────────────────────────────────────
const OUT_DIR     = joinpath(@__DIR__, "..", "drafts", "figures", "numerics")
const CACHE_DIR   = joinpath(OUT_DIR, "sweep_cache")
const CACHE_CKG   = joinpath(CACHE_DIR, "ckg")
const CACHE_GAUSS = joinpath(CACHE_DIR, "dll_gauss")
const CACHE_METRO = joinpath(CACHE_DIR, "dll_metro")
const FIG_PNG     = joinpath(OUT_DIR, "ckg_vs_dll_taumix.png")
const FIG_PDF     = joinpath(OUT_DIR, "ckg_vs_dll_taumix.pdf")
const BSON_OUT    = joinpath(OUT_DIR, "ckg_vs_dll_taumix.bson")

mkpath(OUT_DIR)
mkpath(CACHE_CKG)
mkpath(CACHE_GAUSS)
mkpath(CACHE_METRO)

# ── Threading info ────────────────────────────────────────────────────────────
println("="^72)
println("CKG vs DLL τ_mix comparison sweep (qf-lkb.6)")
println("="^72)
@printf("Julia threads: %d, BLAS threads: %d\n",
        Threads.nthreads(), BLAS.get_num_threads())
@printf("n_values    : %s\n", string(N_VALUES))
@printf("beta_values : %s\n", string(BETA_VALUES))
@printf("seeds       : %s\n", string(SEEDS))
@printf("target ε    : %.1e\n", TARGET_EPS)
@printf("t_max_factor: %.2f / gap_est\n", T_MAX_FACTOR)
@printf("t_grid_len  : %d\n", T_GRID_LEN)
println("="^72)
flush(stdout)

# ── Run the three sweeps ──────────────────────────────────────────────────────
# `sweep_mixing_times` accepts `filter::Union{Nothing, AbstractFilter}` and
# rebuilds DLL filter instances internally per (n, β). Using a *type tag*
# (a representative instance at any β) is sufficient.

t_total_start = time()

println("\n[1/3] CKG (KMS, EnergyDomain) sweep…");  flush(stdout)
t_ckg_start = time()
results_ckg = sweep_mixing_times(
    N_VALUES, BETA_VALUES;
    construction   = KMS(),
    domain         = EnergyDomain(),         # qf-lkb.11: production CKG path
    filter         = nothing,
    mode           = :L,
    method         = :krylov,                # qf-e4y.7: eigenmode τ_mix
    seeds          = SEEDS,
    a              = 0.0,                    # qf-lkb.11: thesis smooth-Metropolis
    s              = 0.25,                   # qf-lkb.11: thesis smooth-Metropolis
    target_epsilon = TARGET_EPS,
    t_max_factor   = T_MAX_FACTOR,
    t_grid_length  = T_GRID_LEN,
    krylovdim      = KRYLOV_DIM,
    tol            = TOL,
    output_dir     = CACHE_CKG,
    use_threads    = false,  # BLAS handles parallelism; Julia-level threading causes oversubscription at n=5
    skip_existing  = true,
)
@printf("    CKG sweep wall time : %.1f s\n", time() - t_ckg_start)
flush(stdout)

println("\n[2/3] DLL Gaussian sweep…"); flush(stdout)
t_g_start = time()
results_dll_gauss = sweep_mixing_times(
    N_VALUES, BETA_VALUES;
    construction   = DLL(),
    domain         = BohrDomain(),               # DLL is BohrDomain only
    filter         = DLLGaussianFilter(1.0),     # type tag; harness rebuilds per β
    mode           = :L,
    method         = :krylov,                    # qf-e4y.7: eigenmode τ_mix
    seeds          = SEEDS,
    a              = 0.0,                        # DLL ignores a, s but keep
    s              = 0.25,                       # config consistent across runs
    target_epsilon = TARGET_EPS,
    t_max_factor   = T_MAX_FACTOR,
    t_grid_length  = T_GRID_LEN,
    krylovdim      = KRYLOV_DIM,
    tol            = TOL,
    output_dir     = CACHE_GAUSS,
    use_threads    = false,  # BLAS handles parallelism; Julia-level threading causes oversubscription at n=5
    skip_existing  = true,
)
@printf("    DLL Gaussian sweep wall time : %.1f s\n", time() - t_g_start)
flush(stdout)

println("\n[3/3] DLL Metropolis sweep…"); flush(stdout)
t_m_start = time()
results_dll_metro = sweep_mixing_times(
    N_VALUES, BETA_VALUES;
    construction   = DLL(),
    domain         = BohrDomain(),               # DLL is BohrDomain only
    filter         = DLLMetropolisFilter(1.0),   # type tag; harness rebuilds per β
    mode           = :L,
    method         = :krylov,                    # qf-e4y.7: eigenmode τ_mix
    seeds          = SEEDS,
    a              = 0.0,
    s              = 0.25,
    target_epsilon = TARGET_EPS,
    t_max_factor   = T_MAX_FACTOR,
    t_grid_length  = T_GRID_LEN,
    krylovdim      = KRYLOV_DIM,
    tol            = TOL,
    output_dir     = CACHE_METRO,
    use_threads    = false,  # BLAS handles parallelism; Julia-level threading causes oversubscription at n=5
    skip_existing  = true,
)
@printf("    DLL Metropolis sweep wall time : %.1f s\n", time() - t_m_start)
flush(stdout)

t_total = time() - t_total_start

# ── Summary table ─────────────────────────────────────────────────────────────
function print_summary(label::String, results::Vector{<:NamedTuple})
    println("\n", "="^88)
    println("Summary: ", label)
    println("="^88)
    # qf-e4y.7: switched from biexp diagnostics (fitted_gap / r²) to the
    # eigenmode schema. gap_est comes from the Arnoldi spectral pass,
    # floor_distance is the captured ‖ρ_inf - σ_β‖_1 / 2, source ∈
    # {:extrapolated, :floor, :nan}.
    @printf("%-3s %-6s %-12s %-12s %-12s %-13s %-10s\n",
            "n", "β", "τ_mix", "gap_est", "floor_dist", "source", "all_conv")
    println("-"^88)
    sorted = sort(collect(results); by = r -> (r.n, r.beta))
    for r in sorted
        @printf("%-3d %-6.1f %-12.4e %-12.4e %-12.4e %-13s %-10s\n",
                r.n, r.beta, r.mixing_time, r.gap_est, r.floor_distance,
                string(r.mixing_time_source), string(r.all_converged))
    end
end

print_summary("CKG (KMS)",        results_ckg)
print_summary("DLL Gaussian",     results_dll_gauss)
print_summary("DLL Metropolis",   results_dll_metro)
@printf("\nTotal sweep wall time: %.1f s\n", t_total)
flush(stdout)

# ── Persist raw results ───────────────────────────────────────────────────────
# BSON.@save accepts namedtuple vectors directly when their concrete element
# type is round-trippable; cast to Vector{Any} for the BSON layer's safety.
let
    rc_any = Vector{Any}(results_ckg)
    rg_any = Vector{Any}(results_dll_gauss)
    rm_any = Vector{Any}(results_dll_metro)
    BSON.bson(BSON_OUT, Dict(
        :results_ckg        => rc_any,
        :results_dll_gauss  => rg_any,
        :results_dll_metro  => rm_any,
        :n_values           => N_VALUES,
        :beta_values        => BETA_VALUES,
        :target_epsilon     => TARGET_EPS,
        :t_max_factor       => T_MAX_FACTOR,
        :t_grid_length      => T_GRID_LEN,
    ))
    @info "Saved raw sweep data" BSON_OUT
end

# ── Build the figure ──────────────────────────────────────────────────────────
# Helper: pluck τ_mix for given (n, β) from a sweep result vector.
function get_taumix(results::Vector{<:NamedTuple}, n::Integer, β::Real)
    idx = findfirst(r -> r.n == n && isapprox(r.beta, β; atol=1e-12), results)
    idx === nothing && return NaN
    return results[idx].mixing_time
end

function build_panel_vs_n(results::Vector{<:NamedTuple}, title::String;
                          ylabel::String = "")
    plt = Plots.plot(
        title = title,
        xlabel = "n",
        ylabel = ylabel,
        yscale = :log10,
        legend = :outerright,
        xticks = N_VALUES,
        guidefontsize = 10,
        titlefontsize = 11,
    )
    for β in BETA_VALUES
        τs = [get_taumix(results, n, β) for n in N_VALUES]
        # Drop NaNs to keep the line connected through valid cells
        valid_idx = findall(!isnan, τs)
        if isempty(valid_idx)
            continue
        end
        Plots.plot!(plt, N_VALUES[valid_idx], τs[valid_idx];
            color = COLOR_BETA[β],
            marker = :circle,
            markersize = 5,
            linewidth = 2,
            label = @sprintf("β=%.0f", β))
    end
    return plt
end

function build_panel_vs_beta(results::Vector{<:NamedTuple}, title::String;
                             ylabel::String = "", xlabel::String = "β")
    plt = Plots.plot(
        title = title,
        xlabel = xlabel,
        ylabel = ylabel,
        yscale = :log10,
        legend = :outerright,
        xticks = BETA_VALUES,
        guidefontsize = 10,
        titlefontsize = 11,
    )
    for n in N_VALUES
        τs = [get_taumix(results, n, β) for β in BETA_VALUES]
        valid_idx = findall(!isnan, τs)
        if isempty(valid_idx)
            continue
        end
        Plots.plot!(plt, BETA_VALUES[valid_idx], τs[valid_idx];
            color = COLOR_N[n],
            marker = :diamond,
            markersize = 5,
            linewidth = 2,
            label = "n=$n")
    end
    return plt
end

# Y-axis labels: only on the leftmost panel of each row.
ylab_top = raw"$\tau_{\mathrm{mix}}$"
ylab_bot = raw"$\tau_{\mathrm{mix}}$"

# Top row: τ_mix vs n
p1 = build_panel_vs_n(results_ckg,        "CKG (KMS)";       ylabel = ylab_top)
p2 = build_panel_vs_n(results_dll_gauss,  "DLL Gaussian";    ylabel = "")
p3 = build_panel_vs_n(results_dll_metro,  "DLL Metropolis";  ylabel = "")

# Bottom row: τ_mix vs β
p4 = build_panel_vs_beta(results_ckg,       "CKG (KMS)";       ylabel = ylab_bot)
p5 = build_panel_vs_beta(results_dll_gauss, "DLL Gaussian";    ylabel = "")
p6 = build_panel_vs_beta(results_dll_metro, "DLL Metropolis";  ylabel = "")

# Top-row panels: blank x-axis label (only bottom row should have one).
for p in (p1, p2, p3)
    Plots.plot!(p; xlabel = "")
end

fig = Plots.plot(p1, p2, p3, p4, p5, p6;
    layout = (2, 3),
    size = (1500, 800),
    plot_title = "Mixing time τ_mix from ODE-integrator + bi-exp extrapolation, target ε=1e-3",
    plot_titlefontsize = 12,
    bottom_margin = 4Plots.mm,
    left_margin   = 5Plots.mm,
    right_margin  = 2Plots.mm,
    top_margin    = 2Plots.mm,
)

Plots.savefig(fig, FIG_PNG)
Plots.savefig(fig, FIG_PDF)
@info "Saved figure" FIG_PNG FIG_PDF

println("\nDone.")
