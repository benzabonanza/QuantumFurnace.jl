#!/usr/bin/env julia
# analyze_qf_e4z_35_sigma_sweep_plot.jl  (qf-e4z.35)
#
# Consumes the qf-e4z.35 plot-grade σ-sweep sidecars (see
# scratch_qf_e4z_35_sigma_sweep_plot.jl) and emits a plot-ready summary
# BSON + a stdout table.
#
# Per (n, β_phys) we report τ_mix, gap, ‖L‖_HS, d_{1→1} along the c-grid
# and three diagnostic ratios that disentangle "structural" vs "rate-scale"
# components of the σ-trend (cf. scratch_qf_e4z_34_norm_diagnostic.jl
# and scratch_sigma_sweep_norm_proxies.jl):
#
#   • τ_mix · ‖L‖_HS  ("mixing in L-time-units") — flat in c ⇒ σ-trend in
#     τ_mix is pure rate-scaling of L; rising in c ⇒ structural slow-down.
#   • τ_mix · gap     (gap-resolved mixing) — flat ⇒ τ_mix is gap-driven
#     (Poincaré bound saturated); drifting ⇒ slowest-mode coefficient |c_2|
#     changes with σ.
#   • gap / ‖L‖_HS    ("relative spectral gap" of L) — measures where in
#     the spectrum the slow mode sits.
#
# Output:
#   scripts/output/sweep_qf_e4z_35_sigma_sweep_plot_summary_stats.bson
#     :rows        => Vector{Dict} (one row per (n, β_phys, c)) — all metadata
#     :provenance  => sidecar dir
#     :timestamp   => UTC
#
# Usage:
#   julia --project scripts/analyze_qf_e4z_35_sigma_sweep_plot.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using QuantumFurnace
using BSON
using Printf
using Statistics
using Dates

const OUTPUT_ROOT  = joinpath(@__DIR__, "output")
const OUT_ROOT     = joinpath(OUTPUT_ROOT, "sweep_qf_e4z_35_sigma_sweep_plot")
const OUT_CKG      = joinpath(OUT_ROOT, "ckg")
const STATS_BSON   = joinpath(OUTPUT_ROOT, "sweep_qf_e4z_35_sigma_sweep_plot_summary_stats.bson")

const BETA_PHYS_ALL = (0.25, 0.5, 1.0)
const SIGMA_FACTORS = (0.25, 0.5, 0.75, 1.0, 1.5, 2.0)
const N_MIN, N_MAX  = 3, 8
const SEED          = 46

function _bp_str(beta_phys::Real)
    s = @sprintf("%.6f", float(beta_phys))
    s = rstrip(s, '0'); s = rstrip(s, '.')
    return isempty(s) ? "0" : s
end

function _c_str(sigma_factor::Real)
    s = @sprintf("%.4f", float(sigma_factor))
    s = rstrip(s, '0'); s = rstrip(s, '.')
    return isempty(s) ? "0" : s
end

function _sidecar_path(n::Integer, beta_phys::Real, sigma_factor::Real, seed::Integer)
    return joinpath(OUT_CKG,
        "sweep_n$(n)_betaphys$(_bp_str(beta_phys))_sigma$(_c_str(sigma_factor))_seed$(seed)_L_KMS_Energy.bson")
end

function _load_all()
    rows = NamedTuple[]
    for n in N_MIN:N_MAX, β in BETA_PHYS_ALL, c in SIGMA_FACTORS
        path = _sidecar_path(n, float(β), float(c), SEED)
        if isfile(path)
            try
                d = BSON.load(path, QuantumFurnace)
                r = d[:result]
                push!(rows, (; (Symbol(k) => v for (k, v) in pairs(r))...))
            catch err
                @warn "failed to load $path" err
            end
        else
            @warn "missing sidecar" n β_phys=β c seed=SEED path
        end
    end
    return rows
end

function _row_for(rows, n::Integer, β::Real, c::Real)
    i = findfirst(r -> r.n == n &&
                       isapprox(r.beta_phys, float(β); atol=1e-12) &&
                       isapprox(r.sigma_factor, float(c); atol=1e-12),
                  rows)
    return i === nothing ? nothing : rows[i]
end

function main()
    println("="^110)
    println("qf-e4z.35 σ-sweep summary (single seed=$SEED, n=$N_MIN..$N_MAX, β_phys=$BETA_PHYS_ALL)")
    println("="^110)
    rows = _load_all()
    @printf("[load] %d / %d expected cells loaded from %s\n",
            length(rows),
            (N_MAX - N_MIN + 1) * length(BETA_PHYS_ALL) * length(SIGMA_FACTORS),
            OUT_CKG)
    flush(stdout)

    # Header
    @printf("\n%-3s %-7s %-6s | %-11s %-11s | %-10s %-10s | %-11s %-11s %-11s | %-8s\n",
            "n", "β_phys", "c",
            "τ_mix_phys", "gap_phys",
            "‖L‖_HS_phys", "d_{1→1}_phys",
            "τ_mix·‖L‖", "τ_mix·gap", "gap/‖L‖",
            "src")
    println("-"^110)

    # Per (n, β_phys, c) row
    for n in N_MIN:N_MAX, β in BETA_PHYS_ALL
        first_in_block = true
        for c in SIGMA_FACTORS
            r = _row_for(rows, n, β, c)
            r === nothing && continue
            τ = r.mixing_time_phys
            g = r.gap_phys
            N = r.hs_norm_phys
            d = r.d_1to1_phys
            τN  = isfinite(τ) ? τ * N : NaN
            τg  = isfinite(τ) ? τ * g : NaN
            gN  = g / N
            sep = first_in_block ? "" : ""
            @printf("%-3d %-7g %-6g | %-11.4g %-11.4g | %-10.4g %-10.4g | %-11.4g %-11.4g %-11.4g | %-8s\n",
                    n, β, c,
                    τ, g, N, d,
                    τN, τg, gN,
                    string(r.mixing_time_source))
            first_in_block = false
        end
        println("-"^110)
    end

    # Diagnostic call-outs:
    # - flag cells with floor (τ_mix = Inf)
    floors = filter(r -> r.mixing_time_source == :floor, rows)
    if !isempty(floors)
        println("\n[diag] τ_mix = Inf (floor) cells:")
        for r in sort(floors, by = r -> (r.n, r.beta_phys, r.sigma_factor))
            @printf("    n=%d β_phys=%g c=%g  σ_alg=%.4g  floor_dist=%.4g  kdim=%d\n",
                    r.n, r.beta_phys, r.sigma_factor,
                    r.sigma_alg, r.floor_distance, r.krylovdim)
        end
    end

    # - flag non-converged Arnoldi
    not_conv = filter(r -> !r.all_converged, rows)
    if !isempty(not_conv)
        println("\n[diag] Arnoldi not all-converged cells:")
        for r in sort(not_conv, by = r -> (r.n, r.beta_phys, r.sigma_factor))
            @printf("    n=%d β_phys=%g c=%g  mv=%d  kdim=%d\n",
                    r.n, r.beta_phys, r.sigma_factor, r.total_matvecs, r.krylovdim)
        end
    end

    # Save the stats BSON
    rows_dict = Dict[]
    for r in rows
        push!(rows_dict, Dict(pairs(r)...))
    end
    BSON.bson(STATS_BSON, Dict(
        :rows        => rows_dict,
        :n_range     => collect(N_MIN:N_MAX),
        :beta_phys_grid => collect(BETA_PHYS_ALL),
        :sigma_factor_grid => collect(SIGMA_FACTORS),
        :seed        => SEED,
        :provenance  => OUT_CKG,
        :timestamp   => string(now()),
    ))
    @printf("\n[save] %d rows → %s\n", length(rows_dict), STATS_BSON)
    println("[done] $(now())")
end

main()
