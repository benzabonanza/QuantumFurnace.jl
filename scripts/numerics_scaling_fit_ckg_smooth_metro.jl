#!/usr/bin/env julia
# numerics_scaling_fit_ckg_smooth_metro.jl  (qf-6vr.8 / Task 8)
#
# Promoted from the gitignored `scratch_scaling_fit_ckg_smooth_metro.jl`
# placeholder named in the qf-bphys plan — the β_phys-first re-run is the
# canonical Task 8 driver, so it lives under `numerics_*` with the rest of
# the production sweep scripts.
#
# Re-run the CKG smooth-Metropolis (n, β_phys) sweep + scaling-law fit under
# the qf-6vr β_phys-first contract. Replaces the legacy β_alg-first run that
# produced `drafts/figures/numerics/scaling_fit_ckg_smooth_metro.{png,pdf}`
# (with x≈0.82, y≈0.42 fitted exponents now suspected of being a β_phys/n
# coupling artifact — see drafts/scaling-fit-physics-check.md).
#
# Sweep grid (Phase A: smoke / laptop budget):
#   n_values       = [3, 4, 5, 6, 7, 8]        — fits in a few hours on the laptop
#   β_phys_values  = [0.25, 0.5, 1.0]          — qf-6vr canonical grid (replaces
#                                                legacy β_alg ∈ {5, 10, 20})
#   construction   = KMS, domain = EnergyDomain, filter = nothing
#   smooth-Metro: a = 0, s = 0.25 (legacy fixed value — see CAVEAT below)
#   target_eps     = 1e-3, method = :krylov, seeds = [42] (single seed)
#
# CAVEAT (qf-6vr / `s` blowup at β_phys=1, n≥9): the σ=1/β_alg convention
# combined with the qf-96o `default_smooth_s(β,σ) = (0.05/σ)²` rule would
# force s ≈ 25 at (β_phys=1, n=11) to preserve the absolute kink width
# σ·√s = 0.05. We do not know yet how a smooth-Metro kernel with s=O(10)
# affects the γ-rate suppression vs the optimal kinky Metropolis, nor what
# it does to the τ_mix. This script holds s = 0.25 fixed (legacy thesis
# convention) until that question is resolved. If the β_phys=1 cells
# scale terribly, revisit either (a) the s-vs-β rule, (b) the σ-vs-β
# convention (σ = c/β with c < 1 per `sigma_sweep_findings_qf_bw1`), or
# (c) reinstate fixed σ = 0.1 as a control sweep.
#
# Outputs:
#   scripts/output/sweep_S1_ckg_ideal_betaphys/                     per-cell BSON sidecars
#   drafts/figures/numerics/scaling_fit_ckg_smooth_metro_betaphys.{png,pdf}   diagnostic figure
#   drafts/scaling-fit-bphys-rerun.md                                results note
#
# PHYSICS CHECK (qf-6vr): the sweep harness reads
# `ham.rescaling_factor` for the family-specific n-fixture and derives
# `β_alg = β_phys · rescaling_factor` per (n) cell. The Hamiltonian family
# `heis_xxx_zzdisordered_periodic_n*` has `rescaling_factor` ≈ 20–30 at n=3
# growing roughly linearly with n. β_phys=3 at n=8 ⇒ β_alg ≈ 30 × 3 ≈ 90 —
# beyond the calibration envelope of `quadrature_register_recipe_qf_7xt.md`.
# This script uses BohrDomain / EnergyDomain (analytical α-form, no per-term
# register table required) so the register-sizing audit deferred for the
# channel-side sweep (S3) does not block S1.
#
# Usage:
#   JULIA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 \
#       julia --project scripts/scratch_scaling_fit_ckg_smooth_metro.jl
#
# Smoke mode (single cell, n=3 / β_phys=1): set the env var SMOKE=1
#   SMOKE=1 julia --project scripts/scratch_scaling_fit_ckg_smooth_metro.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using QuantumFurnace
using LinearAlgebra
using Printf
using BSON
using Dates

BLAS.set_num_threads(1)
println("[init] $(now())  Julia threads = $(Threads.nthreads()), BLAS = $(BLAS.get_num_threads())")
println("[init] hostname = $(gethostname())")

# --- Sweep grid -----------------------------------------------------------

const SMOKE_MODE = get(ENV, "SMOKE", "0") == "1"

const N_VALUES         = SMOKE_MODE ? [3]              : [3, 4, 5, 6, 7, 8]
# qf-6vr canonical β_phys grid (decided 2026-05-11 from the Gibbs-state
# entropy analysis): three log-spaced values factor-2 apart, factor-4 total
# span. 0.25 is the lower bound (below is essentially uniform/infinite-T:
# S/log(d) > 0.96 across all n); 1.0 is the practical upper bound (σ at
# n=10 is already 0.011, and σ·√s_default would force s ≈ 25, which is
# physically defined but pushes the smooth-Metro kernel into a regime
# we haven't tested — `default_smooth_s(β,σ)` may need to be capped or
# the legacy s=0.25 reinstated if the τ_mix at β_phys=1 looks bad).
const BETA_PHYS_VALUES = SMOKE_MODE ? [0.25]           : [0.25, 0.5, 1.0]
const TARGET_EPS       = 1e-3
const T_MAX_FACTOR     = 5.0
const T_GRID_LEN       = 81
const KRYLOV_DIM       = 30
const SPECTRAL_KDIM    = 60
const TOL              = 1e-10
const SEEDS            = [42]

const HAM_FILENAME     = (n) -> "heis_xxx_zzdisordered_periodic_n$(n).bson"
const HAM_DIR          = joinpath(@__DIR__, "..", "hamiltonians")

const OUT_DIR          = joinpath(@__DIR__, "output", "sweep_S1_ckg_ideal_betaphys")
const FIG_DIR          = joinpath(@__DIR__, "..", "drafts", "figures", "numerics")
const FIG_PNG          = joinpath(FIG_DIR, "scaling_fit_ckg_smooth_metro_betaphys.png")
const FIG_PDF          = joinpath(FIG_DIR, "scaling_fit_ckg_smooth_metro_betaphys.pdf")
const NOTE_PATH        = joinpath(@__DIR__, "..", "drafts", "scaling-fit-bphys-rerun.md")

mkpath(OUT_DIR)
mkpath(FIG_DIR)

# --- Run the sweep --------------------------------------------------------

println("\n=== CKG smooth-Metropolis β_phys-first sweep ===")
@printf("n_values         : %s\n", string(N_VALUES))
@printf("beta_phys_values : %s\n", string(BETA_PHYS_VALUES))
@printf("smoke mode       : %s\n", SMOKE_MODE)
@printf("output dir       : %s\n", OUT_DIR)
flush(stdout)

t_sweep_start = time()
results = sweep_mixing_times(
    N_VALUES;
    beta_phys_values     = BETA_PHYS_VALUES,
    construction         = KMS(),
    domain               = EnergyDomain(),
    filter               = nothing,
    mode                 = :L,
    method               = :krylov,
    seeds                = SEEDS,
    a                    = 0.0,
    s                    = 0.25,
    target_epsilon       = TARGET_EPS,
    t_max_factor         = T_MAX_FACTOR,
    t_grid_length        = T_GRID_LEN,
    krylovdim            = KRYLOV_DIM,
    spectral_krylovdim   = SPECTRAL_KDIM,
    tol                  = TOL,
    output_dir           = OUT_DIR,
    hamiltonian_dir      = HAM_DIR,
    hamiltonian_filename = HAM_FILENAME,
    use_threads          = false,  # BLAS handles parallelism for EnergyDomain matvecs
    skip_existing        = true,
)
t_sweep_total = time() - t_sweep_start
@printf("[sweep] done in %.1f s — %d cells\n", t_sweep_total, length(results))
flush(stdout)

# --- Summary table --------------------------------------------------------

println("\n", "="^88)
println("Per-cell summary (β_phys / β_alg / τ_mix / source)")
println("="^88)
@printf("%-3s %-8s %-8s %-10s %-12s %-12s %-13s\n",
        "n", "β_phys", "β_alg", "rescale", "τ_mix", "gap_est", "source")
println("-"^88)
sorted = sort(collect(results); by = r -> (r.n, r.beta_phys))
for r in sorted
    rescale = haskey(r, :rescaling_factor) ? r.rescaling_factor : NaN
    @printf("%-3d %-8.2f %-8.2f %-10.2f %-12.4e %-12.4e %-13s\n",
            r.n, r.beta_phys, r.beta_alg, rescale,
            r.mixing_time, r.gap_est,
            string(r.mixing_time_source))
end

if SMOKE_MODE
    println("\n[SMOKE] single-cell smoke test complete — full sweep skipped.")
    println("[SMOKE] Inspect $(OUT_DIR) for the sidecar, then unset SMOKE to run the full grid.")
    exit(0)
end

# --- Scaling-fit on the full grid ----------------------------------------

println("\n=== Scaling-law fit (qf-now / qf-6vr) ===")
fits = fit_scaling(results)  # :auto → reads :beta_phys preferentially, kind = :phys
ranking = compare_models(fits)
println("Model ranking (best AICc first):")
for (i, m) in enumerate(ranking.ranked)
    @printf("  %d. %s   AICc=%.3f   Δ=%.3f   weight=%.3f\n",
            i, m, ranking.aicc[i], ranking.delta_aicc[i], ranking.weights[i])
end
for (m, f) in fits
    println("\n--- $m ---")
    println("  ", formula_string(f))
    println("  beta_kind = $(f.beta_kind)")
    println("  σ_residual = $(round(f.sigma_residual; digits=4))")
end

# --- Persist the figure + results note -----------------------------------

ENV["GKSwstype"] = "100"
using Plots

best = fits[ranking.ranked[1]]
grid = scaling_fit_grid(best)

# Two-panel diagnostic: data vs fit + residuals.
function _make_diagnostic(grid_data, fit::ScalingFit)
    p1 = Plots.scatter(grid_data.beta_obs, grid_data.tau_obs;
        xscale = :log10, yscale = :log10,
        xlabel = "β_phys", ylabel = "τ_mix",
        title = "Data vs fit ($(fit.model))",
        marker_z = grid_data.n_obs,
        legend = false, markersize = 5)
    for nv in unique(grid_data.n_obs)
        ix = findall(==(nv), grid_data.n_obs)
        Plots.plot!(p1, grid_data.beta_obs[ix],
                    grid_data.tau_pred_at_obs[ix];
                    linewidth = 2, alpha = 0.6,
                    label = "")
    end
    p2 = Plots.scatter(grid_data.beta_obs, grid_data.residuals_log;
        xscale = :log10,
        xlabel = "β_phys", ylabel = "log(τ_obs) − log(τ_pred)",
        title = "Log residuals", legend = false, markersize = 5)
    Plots.hline!(p2, [0]; color = :black, linewidth = 1, label = "")
    return Plots.plot(p1, p2; layout = (1, 2), size = (1200, 450))
end

fig = _make_diagnostic(grid, best)
Plots.savefig(fig, FIG_PNG)
Plots.savefig(fig, FIG_PDF)
@info "Saved diagnostic figure" FIG_PNG FIG_PDF

# Write results note
open(NOTE_PATH, "w") do io
    println(io, "# Scaling-fit re-run under β_phys convention (qf-6vr / Task 8)")
    println(io, "")
    println(io, "Sweep date: ", today())
    println(io, "")
    println(io, "## Grid")
    println(io, "")
    println(io, "- `n_values` = ", N_VALUES)
    println(io, "- `β_phys_values` = ", BETA_PHYS_VALUES)
    println(io, "- Construction: CKG (KMS, EnergyDomain), smooth-Metropolis (a=0, s=0.25)")
    println(io, "- Fixture: heis_xxx_zzdisordered_periodic_n* (Z+ZZ disorder)")
    println(io, "- target_ε = ", TARGET_EPS, ", seed = ", SEEDS[1])
    println(io, "")
    println(io, "## Per-cell summary")
    println(io, "")
    println(io, "| n | β_phys | β_alg | rescale | τ_mix | gap_est | source |")
    println(io, "|---|---|---|---|---|---|---|")
    for r in sorted
        rescale = haskey(r, :rescaling_factor) ? r.rescaling_factor : NaN
        println(io, @sprintf("| %d | %.2f | %.2f | %.2f | %.4e | %.4e | %s |",
                             r.n, r.beta_phys, r.beta_alg, rescale,
                             r.mixing_time, r.gap_est, r.mixing_time_source))
    end
    println(io, "")
    println(io, "## Scaling-law fit")
    println(io, "")
    println(io, "Model ranking:")
    println(io, "")
    println(io, "| model | AICc | ΔAICc | weight | formula |")
    println(io, "|---|---|---|---|---|")
    for (i, m) in enumerate(ranking.ranked)
        f = fits[m]
        println(io, @sprintf("| %s | %.3f | %.3f | %.3f | %s |",
                             m, ranking.aicc[i], ranking.delta_aicc[i],
                             ranking.weights[i], formula_string(f)))
    end
    println(io, "")
    println(io, "## Comparison to the legacy β_alg sweep")
    println(io, "")
    println(io, "Legacy fit (`drafts/scaling-fit-physics-check.md`):")
    println(io, "  `τ_mix ≈ 3.48 · n^0.82 · β_alg^0.42` (single seed, β_alg ∈ {5, 10, 20}, n ∈ {3..8})")
    println(io, "")
    println(io, "New fit (this run, β_phys ∈ {1, 2, 3}):")
    println(io, "  ", formula_string(fits[:M0]))
    println(io, "")
    println(io, "## Files")
    println(io, "")
    println(io, "- Sidecars: `", OUT_DIR, "`")
    println(io, "- Figure: `", FIG_PNG, "`, `", FIG_PDF, "`")
end
@info "Saved results note" NOTE_PATH

println("\n=== Done ", now(), " ===")
