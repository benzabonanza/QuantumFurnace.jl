#!/usr/bin/env julia
#
# β=20 (low-temperature) extension of the CKG-vs-DLL τ_mix comparison.
#
# Motivation
# ----------
# The qf-lkb.6 figure (β ∈ {1,2,5,10}, n ∈ {3,4,5}) shows DLL Gaussian
# beating CKG/DLL Metropolis by 30–50% at β ≤ 5, but the gap *narrows*
# at β=10 — at n=3, β=10 DLL Gaussian is actually slower than CKG.
#
# Folklore prediction: DLL has β-independent Kossakowski rank, while
# CKG's Kossakowski matrix has spectral norm shrinking with β. So at
# very low temperature, DLL should pull ahead.
#
# This script tests that prediction by extending to β ∈ {15, 20} on the
# same Heisenberg n ∈ {3,4,5} fixture family, reusing the existing
# adaptive-t_max bi-exp pipeline. β=10 cells are kept so cache hits
# verify reproducibility against the production figure.
#
# Strategy: same matrix-free `sweep_mixing_times` path the production
# figure uses (qf-lkb.10 adaptive horizon, qf-lkb.11 EnergyDomain CKG,
# qf-lkb.9 matrix-free DLL). The only differences vs the production
# script are the β grid and a separate cache directory (so the canonical
# figure data doesn't get mixed up with this experiment).
#
# PHYSICS CHECK: target_epsilon = 1e-3, t_max_factor = :auto (qf-lkb.10
# adaptive horizon — required because gap shrinks fast with β; a fixed
# 5/gap_est would still work but :auto is the tested default).
#
# PHYSICS CHECK: smooth-Metropolis a=0, s=0.25 (thesis convention; locked
# in qf-lkb.11). DLL filters have no a/s knobs but kwargs are ignored.

using Printf
using Random
using LinearAlgebra
using BSON
using QuantumFurnace

# ── Sweep grid ────────────────────────────────────────────────────────────────
const N_VALUES    = [3, 4, 5]
const BETA_VALUES = [10.0, 15.0, 20.0]   # β=10 stays so we can sanity-check
const TARGET_EPS  = 1e-3
const T_MAX_FACTOR = :auto                 # adaptive (qf-lkb.10)
const T_GRID_LEN   = 81
const KRYLOV_DIM   = 30
const TOL          = 1e-10
const SEEDS        = [42]

# ── Output paths (separate cache so production figure stays clean) ────────────
const OUT_DIR     = joinpath(@__DIR__, "output", "taumix_beta20")
const CACHE_DIR   = joinpath(OUT_DIR, "sweep_cache")
const CACHE_CKG   = joinpath(CACHE_DIR, "ckg")
const CACHE_GAUSS = joinpath(CACHE_DIR, "dll_gauss")
const CACHE_METRO = joinpath(CACHE_DIR, "dll_metro")
const BSON_OUT    = joinpath(OUT_DIR, "taumix_beta20.bson")

mkpath(OUT_DIR)
mkpath(CACHE_CKG)
mkpath(CACHE_GAUSS)
mkpath(CACHE_METRO)

println("="^72)
println("CKG vs DLL τ_mix — low-temperature extension β ∈ {10,15,20}")
println("="^72)
@printf("Julia threads: %d, BLAS threads: %d\n",
        Threads.nthreads(), BLAS.get_num_threads())
@printf("n_values    : %s\n", string(N_VALUES))
@printf("beta_values : %s\n", string(BETA_VALUES))
@printf("target ε    : %.1e\n", TARGET_EPS)
@printf("t_max_factor: :auto (gap-adaptive, qf-lkb.10)\n")
println("="^72)
flush(stdout)

t_total_start = time()

println("\n[1/3] CKG (KMS, EnergyDomain) sweep…");  flush(stdout)
t0 = time()
results_ckg = sweep_mixing_times(
    N_VALUES, BETA_VALUES;
    construction   = KMS(),
    domain         = EnergyDomain(),
    filter         = nothing,
    mode           = :L,
    method         = :krylov,                # qf-e4y.7: eigenmode τ_mix
    seeds          = SEEDS,
    a              = 0.0,
    s              = 0.25,
    target_epsilon = TARGET_EPS,
    t_max_factor   = T_MAX_FACTOR,
    t_grid_length  = T_GRID_LEN,
    krylovdim      = KRYLOV_DIM,
    tol            = TOL,
    output_dir     = CACHE_CKG,
    use_threads    = false,
    skip_existing  = true,
)
@printf("    CKG sweep wall time : %.1f s\n", time() - t0); flush(stdout)

println("\n[2/3] DLL Gaussian sweep…"); flush(stdout)
t0 = time()
results_dll_gauss = sweep_mixing_times(
    N_VALUES, BETA_VALUES;
    construction   = DLL(),
    domain         = BohrDomain(),
    filter         = DLLGaussianFilter(1.0),
    mode           = :L,
    method         = :krylov,                # qf-e4y.7: eigenmode τ_mix
    seeds          = SEEDS,
    a              = 0.0,
    s              = 0.25,
    target_epsilon = TARGET_EPS,
    t_max_factor   = T_MAX_FACTOR,
    t_grid_length  = T_GRID_LEN,
    krylovdim      = KRYLOV_DIM,
    tol            = TOL,
    output_dir     = CACHE_GAUSS,
    use_threads    = false,
    skip_existing  = true,
)
@printf("    DLL Gaussian sweep wall time : %.1f s\n", time() - t0); flush(stdout)

println("\n[3/3] DLL Metropolis sweep…"); flush(stdout)
t0 = time()
results_dll_metro = sweep_mixing_times(
    N_VALUES, BETA_VALUES;
    construction   = DLL(),
    domain         = BohrDomain(),
    filter         = DLLMetropolisFilter(1.0),
    mode           = :L,
    method         = :krylov,                # qf-e4y.7: eigenmode τ_mix
    seeds          = SEEDS,
    a              = 0.0,
    s              = 0.25,
    target_epsilon = TARGET_EPS,
    t_max_factor   = T_MAX_FACTOR,
    t_grid_length  = T_GRID_LEN,
    krylovdim      = KRYLOV_DIM,
    tol            = TOL,
    output_dir     = CACHE_METRO,
    use_threads    = false,
    skip_existing  = true,
)
@printf("    DLL Metropolis sweep wall time : %.1f s\n", time() - t0); flush(stdout)

t_total = time() - t_total_start

# ── Results table ─────────────────────────────────────────────────────────────
function print_summary(label::String, results::Vector{<:NamedTuple})
    println("\n", "="^88)
    println("Summary: ", label)
    println("="^88)
    # qf-e4y.7: eigenmode schema (gap_est / floor_distance / source).
    # Fall back to legacy biexp diagnostics if a row is from an :ode-route
    # cache (handles mixed-vintage BSON sweeps).
    has_eigenmode = !isempty(results) && haskey(results[1], :gap_est) &&
                     haskey(results[1], :floor_distance)
    if has_eigenmode
        @printf("%-3s %-6s %-12s %-12s %-12s %-13s %-10s\n",
                "n", "β", "τ_mix", "gap_est", "floor_dist", "source", "all_conv")
    else
        @printf("%-3s %-6s %-12s %-12s %-9s %-10s %-10s\n",
                "n", "β", "τ_mix", "fitted_gap", "R²", "fit_conv", "all_conv")
    end
    println("-"^88)
    sorted = sort(collect(results); by = r -> (r.n, r.beta))
    for r in sorted
        if has_eigenmode
            @printf("%-3d %-6.1f %-12.4e %-12.4e %-12.4e %-13s %-10s\n",
                    r.n, r.beta, r.mixing_time, r.gap_est, r.floor_distance,
                    string(r.mixing_time_source), string(r.all_converged))
        else
            @printf("%-3d %-6.1f %-12.4e %-12.4e %-9.4f %-10s %-10s\n",
                    r.n, r.beta, r.mixing_time, r.fitted_gap, r.r_squared,
                    string(r.converged_fit), string(r.all_converged))
        end
    end
end

print_summary("CKG (KMS, smooth-Metropolis a=0,s=0.25)", results_ckg)
print_summary("DLL Gaussian",   results_dll_gauss)
print_summary("DLL Metropolis", results_dll_metro)

# ── Side-by-side comparison ───────────────────────────────────────────────────
function pluck(res::Vector{<:NamedTuple}, n::Int, β::Float64)
    idx = findfirst(r -> r.n == n && isapprox(r.beta, β; atol=1e-12), res)
    idx === nothing ? NaN : res[idx].mixing_time
end

println("\n", "="^88)
println("Head-to-head τ_mix (lower is faster)")
println("="^88)
@printf("%-3s %-6s | %-12s %-12s %-12s | %-12s %-12s\n",
        "n", "β", "CKG (Metro)", "DLL Gauss", "DLL Metro", "Gauss/CKG", "DMetro/CKG")
println("-"^88)
for n in N_VALUES, β in BETA_VALUES
    τ_ckg   = pluck(results_ckg,         n, β)
    τ_gauss = pluck(results_dll_gauss,   n, β)
    τ_metro = pluck(results_dll_metro,   n, β)
    @printf("%-3d %-6.1f | %-12.4e %-12.4e %-12.4e | %-12.4f %-12.4f\n",
            n, β, τ_ckg, τ_gauss, τ_metro,
            τ_gauss / τ_ckg, τ_metro / τ_ckg)
end

@printf("\nTotal sweep wall time: %.1f s\n", t_total)
flush(stdout)

# ── Persist raw results ───────────────────────────────────────────────────────
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
        :t_grid_length      => T_GRID_LEN,
    ))
    @info "Saved raw sweep data" BSON_OUT
end

println("\nDone.")
