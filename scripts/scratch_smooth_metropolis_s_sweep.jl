#!/usr/bin/env julia
#
# Smooth-Metropolis s-sweep mixing-time experiment (qf-3il).
#
# Headline question
# -----------------
# How does the CKG smooth-Metropolis mixing time τ_mix change with the
# smoothing parameter s ∈ {0, 0.1, 0.25, 0.5, 1.0} at β ∈ {5, 10, 20}
# for n ∈ {3, 4, 5}? Combined with the quadrature-cost analysis (a 4-bit
# saving at ε = 1e-4 going from kinky → smooth, uniformly across β; see
# scripts/scratch_smooth_metropolis_quadrature_error.jl), this fixes the
# Pareto-optimal s.
#
# Prediction
# ----------
# - prop:optim-metro: γ_M (= s = 0, kinky) is optimal among Metropolis-like
#   transitions → s > 0 gives a strictly smaller spectral gap.
# - Per the thesis (eq:smooth-metro), the smooth γ_M^{(s)} ≤ γ_M^{(0)}
#   pointwise, so τ_mix(s) ≥ τ_mix(0). The question is *how much*.
# - "Small s" should give negligible mixing-time degradation while
#   recovering ε-quadrature ancillas (4 bits at ε = 1e-4 per the Gevrey-1/2
#   prop:smooth-metro-gevrey).
#
# Sweep grid (5 × 3 × 3 = 45 cells)
# --------------------------------
#   s     ∈ {0.0, 0.1, 0.25, 0.5, 1.0}
#   β     ∈ {5.0, 10.0, 20.0}
#   n     ∈ {3, 4, 5}                 — laptop n cap (n=6 deferred)
#   a = 0.0  (thesis-numerics convention; locked at qf-lkb.11)
#   σ = 1/β  (CKG, locked thesis convention)
#
# Fixture: heis_xxx_zzdisordered_periodic_n*.bson — covers n=3..10 and
# eliminates bipartite collisions (memory: Bohr Frequency Collision Root
# Cause). Same family as the qf-lkb.11 EnergyDomain CKG path uses.
#
# Domain: EnergyDomain (qf-lkb.11). Matrix-free `apply_lindbladian!`.
#
# PHYSICS CHECK: target_ε = 1e-3, t_max_factor = :auto (qf-lkb.10
# adaptive horizon for bi-exp extrapolation).
#
# Output
# ------
# - drafts/figures/numerics/sweep_cache/smooth_metro_s/<s_tag>/<sidecar>.bson
#   per-cell BSON via `sweep_mixing_times(output_dir=...)` so re-runs are
#   resumable.
# - drafts/figures/numerics/smooth_metro_s_sweep.bson  (aggregated rows)
# - drafts/figures/numerics/smooth_metro_s_taumix.{png,pdf}
# - Summary table to stdout.

using Printf
using LinearAlgebra
using BSON
using QuantumFurnace
ENV["GKSwstype"] = "100"
using Plots

# ── Sweep grid ────────────────────────────────────────────────────────────────
const N_VALUES    = [3, 4, 5]
const BETA_VALUES = [5.0, 10.0, 20.0]
const S_VALUES    = [0.0, 0.1, 0.25, 0.5, 1.0]
const TARGET_EPS  = 1e-3
const T_GRID_LEN  = 81
const KRYLOV_DIM  = 30
const TOL         = 1e-10
const SEEDS       = [42]

# ── Fixture family (zz-disordered → no bipartite collisions, covers n=3..10) ──
const HAM_FAMILY = n -> "heis_xxx_zzdisordered_periodic_n$(n).bson"

# ── Output paths ──────────────────────────────────────────────────────────────
const OUT_DIR    = joinpath(@__DIR__, "..", "drafts", "figures", "numerics")
const CACHE_ROOT = joinpath(OUT_DIR, "sweep_cache", "smooth_metro_s")
mkpath(OUT_DIR)
mkpath(CACHE_ROOT)
const FIG_PNG    = joinpath(OUT_DIR, "smooth_metro_s_taumix.png")
const FIG_PDF    = joinpath(OUT_DIR, "smooth_metro_s_taumix.pdf")
const BSON_OUT   = joinpath(OUT_DIR, "smooth_metro_s_sweep.bson")

# ── Thesis colour palette ─────────────────────────────────────────────────────
const COLOR_S = Dict(
    0.0  => "#222222",   # black (kinky baseline)
    0.1  => "#5C7794",   # slateblue
    0.25 => "#5F8B8E",   # dustyteal
    0.5  => "#B5654A",   # terracotta
    1.0  => "#7A2E39",   # bordeaux
)
const STYLE_S = Dict(
    0.0  => :dash,
    0.1  => :solid,
    0.25 => :solid,
    0.5  => :solid,
    1.0  => :solid,
)
const MARKER_N = Dict(
    3 => :circle,
    4 => :diamond,
    5 => :utriangle,
)

# Subdir tag per s (filesystem-safe).
s_tag(s::Real) = replace(@sprintf("%g", s), "." => "p")

# ── Threading info ────────────────────────────────────────────────────────────
println("="^72)
println("Smooth-Metropolis s-sweep mixing-time experiment (qf-3il)")
println("="^72)
@printf("Julia threads : %d, BLAS threads: %d\n",
        Threads.nthreads(), BLAS.get_num_threads())
@printf("n_values      : %s\n", string(N_VALUES))
@printf("beta_values   : %s\n", string(BETA_VALUES))
@printf("s_values      : %s\n", string(S_VALUES))
@printf("target ε      : %.1e\n", TARGET_EPS)
@printf("t_max_factor  : :auto (qf-lkb.10 adaptive)\n")
@printf("hamiltonian   : heis_xxx_zzdisordered_periodic_n*.bson\n")
println("="^72)
flush(stdout)

# ── Run the s-sweep ───────────────────────────────────────────────────────────
all_results = NamedTuple[]
t_total_start = time()

for s in S_VALUES
    cache_dir = joinpath(CACHE_ROOT, s_tag(s))
    mkpath(cache_dir)
    @printf("\n>> s = %.2f  (cache: %s)\n", s, cache_dir);  flush(stdout)
    t_s_start = time()
    res = sweep_mixing_times(
        N_VALUES, BETA_VALUES;
        construction        = KMS(),
        domain              = EnergyDomain(),
        filter              = nothing,
        mode                = :L,
        seeds               = SEEDS,
        a                   = 0.0,
        s                   = s,
        target_epsilon      = TARGET_EPS,
        t_max_factor        = :auto,
        t_grid_length       = T_GRID_LEN,
        krylovdim           = KRYLOV_DIM,
        tol                 = TOL,
        output_dir          = cache_dir,
        hamiltonian_filename = HAM_FAMILY,
        use_threads         = false,    # BLAS handles parallelism; avoid oversubscription at n=5
        skip_existing       = true,
    )
    @printf("   wall : %.1f s  (cells: %d)\n", time() - t_s_start, length(res))
    flush(stdout)
    for row in res
        push!(all_results, merge(row, (s = s,)))
    end
end

t_total = time() - t_total_start
@printf("\nTotal wall time: %.1f s\n", t_total)

# ── Persist aggregated BSON ───────────────────────────────────────────────────
bson_payload = Dict{Symbol, Any}(
    :results       => all_results,
    :n_values      => N_VALUES,
    :beta_values   => BETA_VALUES,
    :s_values      => S_VALUES,
    :seeds         => SEEDS,
    :target_epsilon => TARGET_EPS,
    :git_sha       => try
        read(`git -C $(dirname(@__DIR__)) rev-parse HEAD`, String) |> strip
    catch
        "unknown"
    end,
)
BSON.bson(BSON_OUT, bson_payload)
@printf("Wrote BSON: %s\n", BSON_OUT)

# ── Plot: τ_mix vs s, panel per β, line per n ────────────────────────────────
plots_panels = Plots.Plot[]
for beta in BETA_VALUES
    p = plot(;
        xlabel = "s (smoothing parameter)",
        ylabel = "τ_mix",
        title  = "β = $(beta)",
        legend = :topleft,
        framestyle = :box,
        size = (520, 380),
    )
    for n in N_VALUES
        ss = Float64[]
        ts = Float64[]
        for s in S_VALUES
            row = findfirst(r -> r.n == n && r.beta == beta && r.s == s,
                            all_results)
            if row !== nothing
                mt = all_results[row].mixing_time
                if isfinite(mt)
                    push!(ss, s); push!(ts, mt)
                end
            end
        end
        plot!(p, ss, ts;
            marker = MARKER_N[n], markersize = 5, markerstrokewidth = 0,
            linewidth = 2,
            label = "n = $(n)",
        )
    end
    push!(plots_panels, p)
end
fig = plot(plots_panels...; layout = (1, length(BETA_VALUES)),
           size = (520 * length(BETA_VALUES), 380))
savefig(fig, FIG_PNG)
savefig(fig, FIG_PDF)
@printf("Wrote: %s, %s\n", FIG_PNG, FIG_PDF)

# ── Summary table ─────────────────────────────────────────────────────────────
println("\n" * "="^88)
println("SUMMARY — τ_mix vs s")
println("="^88)
for beta in BETA_VALUES
    @printf("\n[β = %.1f]\n", beta)
    @printf("  n        ")
    for s in S_VALUES
        @printf("s=%-7.2f  ", s)
    end
    println()
    @printf("  ---------")
    for _ in S_VALUES
        @printf("-----------")
    end
    println()
    for n in N_VALUES
        @printf("  n=%-5d  ", n)
        for s in S_VALUES
            row = findfirst(r -> r.n == n && r.beta == beta && r.s == s,
                            all_results)
            mt = row === nothing ? NaN : all_results[row].mixing_time
            @printf("%-10.3g ", mt)
        end
        println()
    end
end

# ── Relative-to-kinky table ───────────────────────────────────────────────────
println("\n" * "="^88)
println("τ_mix(s) / τ_mix(s = 0) — fractional degradation from smoothing")
println("="^88)
for beta in BETA_VALUES
    @printf("\n[β = %.1f]\n", beta)
    @printf("  n        ")
    for s in S_VALUES
        @printf("s=%-7.2f  ", s)
    end
    println()
    @printf("  ---------")
    for _ in S_VALUES
        @printf("-----------")
    end
    println()
    for n in N_VALUES
        @printf("  n=%-5d  ", n)
        kinky_row = findfirst(r -> r.n == n && r.beta == beta && r.s == 0.0,
                              all_results)
        kinky = kinky_row === nothing ? NaN : all_results[kinky_row].mixing_time
        for s in S_VALUES
            row = findfirst(r -> r.n == n && r.beta == beta && r.s == s,
                            all_results)
            mt = row === nothing ? NaN : all_results[row].mixing_time
            @printf("%-10.3g ", mt / kinky)
        end
        println()
    end
end
println()
println("Done.")
