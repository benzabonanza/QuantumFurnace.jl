#!/usr/bin/env julia
# scratch_qf_e4z_35_structural_pattern.jl  (qf-e4z.35)
#
# Post-hoc: at which σ_factor does the smooth-Metro CKG sampler reach its
# BEST "structural" mixing, after normalising by an L-norm? Two norms:
#
#   • τ·‖L‖_HS — mixing in HS-time-units. Removes pure rate-scale.
#   • τ·d_{1→1} — mixing in Kossakowski opnorm units (dissipator-only bound).
#
# Reading: the minimum c* of τ·N picks where the GENERATOR L is most
# efficient per unit norm. If τ·N is monotonic decreasing in c → wider σ
# is structurally better; monotonic increasing → smaller σ wins; a
# turning point → there is an optimal σ for structural mixing.
#
# Also reports τ·gap (rate-resolved, ≈ const if the slow mode coefficient
# |c_2| is σ-independent) and gap/‖L‖ (relative spectral gap).
#
# Usage: julia --project scripts/scratch_qf_e4z_35_structural_pattern.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BSON
using QuantumFurnace
using Printf
using Statistics

const STATS_BSON = joinpath(@__DIR__, "output",
    "sweep_qf_e4z_35_sigma_sweep_plot_summary_stats.bson")
const N_RANGE       = 3:7                          # complete-n only
const BETA_PHYS_ALL = (0.25, 0.5, 1.0)
const SIGMA_FACTORS = (0.25, 0.5, 0.75, 1.0, 1.5, 2.0)

function load_rows()
    d = BSON.load(STATS_BSON, QuantumFurnace)
    return d[:rows]
end

function get_metric(rows, n, β, c, key)
    for r in rows
        if r[:n] == n && isapprox(r[:beta_phys], float(β); atol=1e-12) &&
           isapprox(r[:sigma_factor], float(c); atol=1e-12)
            return float(r[key])
        end
    end
    return NaN
end

function summarize_axis(rows, axis_label::String, axis_key::Symbol)
    println("\n" * "="^110)
    println("STRUCTURAL METRIC = $axis_label   (lower = better mixing per unit norm)")
    println("="^110)
    @printf("%-3s %-6s | %-9s %-9s %-9s %-9s %-9s %-9s | %-6s %-8s %-9s\n",
            "n", "β_phys",
            "c=0.25", "c=0.5", "c=0.75", "c=1.0", "c=1.5", "c=2.0",
            "c*", "min", "max/min")
    println("-"^110)

    for n in N_RANGE, β in BETA_PHYS_ALL
        vals = Float64[get_metric(rows, n, β, c, axis_key) for c in SIGMA_FACTORS]
        any(isnan, vals) && continue
        imin = argmin(vals)
        c_star = SIGMA_FACTORS[imin]
        spread = maximum(vals) / minimum(vals)
        @printf("%-3d %-6g | %-9.4g %-9.4g %-9.4g %-9.4g %-9.4g %-9.4g | %-6g %-8.4g %-9.3f\n",
                n, β, vals..., c_star, vals[imin], spread)
    end
end

function summarize_argmin_pattern(rows, axis_label::String, axis_key::Symbol)
    println("\n[$(axis_label)] c* (argmin) lookup table:")
    @printf("%-12s | %-6s %-6s %-6s %-6s %-6s\n",
            "β_phys ↓ / n →", "n=3", "n=4", "n=5", "n=6", "n=7")
    println("-"^60)
    for β in BETA_PHYS_ALL
        cs = String[]
        for n in N_RANGE
            vals = Float64[get_metric(rows, n, β, c, axis_key) for c in SIGMA_FACTORS]
            if any(isnan, vals)
                push!(cs, "--")
            else
                push!(cs, @sprintf("%.2f", SIGMA_FACTORS[argmin(vals)]))
            end
        end
        @printf("%-12g | %-6s %-6s %-6s %-6s %-6s\n", β, cs...)
    end
end

function main()
    println("="^110)
    println("qf-e4z.35 σ-sweep: structural mixing pattern across (n, β_phys)")
    println("Single seed=46;  n ∈ ", collect(N_RANGE), ";  β_phys ∈ ", BETA_PHYS_ALL,
            ";  c ∈ ", SIGMA_FACTORS)
    println("="^110)
    rows = load_rows()
    println("[load] $(length(rows)) cells from $STATS_BSON")

    # 1. τ_mix · ‖L‖_HS (HS-norm normalised mixing)
    summarize_axis(rows, "τ_mix · ‖L‖_HS  (HS-norm normalised)", :hs_norm_phys)
    summarize_argmin_pattern(rows, "τ_mix · ‖L‖_HS", :hs_norm_phys)
    # ⚠ summarize_axis uses :hs_norm_phys as the column — we need τ·N, not N alone.
    # Redo properly below.

    # Helper: τ_mix · X
    function product_metric(rows, n, β, c, key)
        τ = get_metric(rows, n, β, c, :mixing_time_phys)
        X = get_metric(rows, n, β, c, key)
        return τ * X
    end

    function summarize_product(label, key)
        println("\n" * "="^110)
        println("STRUCTURAL METRIC = $label   (lower = better)")
        println("="^110)
        @printf("%-3s %-6s | %-9s %-9s %-9s %-9s %-9s %-9s | %-6s %-8s %-9s\n",
                "n", "β_phys",
                "c=0.25", "c=0.5", "c=0.75", "c=1.0", "c=1.5", "c=2.0",
                "c*", "min", "max/min")
        println("-"^110)
        for n in N_RANGE, β in BETA_PHYS_ALL
            vals = Float64[product_metric(rows, n, β, c, key) for c in SIGMA_FACTORS]
            any(isnan, vals) && continue
            imin = argmin(vals)
            spread = maximum(vals) / minimum(vals)
            @printf("%-3d %-6g | %-9.4g %-9.4g %-9.4g %-9.4g %-9.4g %-9.4g | %-6g %-8.4g %-9.3f\n",
                    n, β, vals..., SIGMA_FACTORS[imin], vals[imin], spread)
        end
        println("\n[$label] c* (argmin) table:")
        @printf("%-12s | %-6s %-6s %-6s %-6s %-6s\n",
                "β_phys ↓ / n →", "n=3", "n=4", "n=5", "n=6", "n=7")
        println("-"^60)
        for β in BETA_PHYS_ALL
            cs = String[]
            for n in N_RANGE
                vals = Float64[product_metric(rows, n, β, c, key) for c in SIGMA_FACTORS]
                if any(isnan, vals)
                    push!(cs, "--")
                else
                    push!(cs, @sprintf("%.2f", SIGMA_FACTORS[argmin(vals)]))
                end
            end
            @printf("%-12g | %-6s %-6s %-6s %-6s %-6s\n", β, cs...)
        end
    end

    summarize_product("τ_mix · ‖L‖_HS",   :hs_norm_phys)
    summarize_product("τ_mix · d_{1→1}",  :d_1to1_phys)
    summarize_product("τ_mix · gap_phys", :gap_phys)
    summarize_product("τ_mix",            :rescaling_factor)  # rescaling_factor=R; τ_mix·R = τ_mix_alg

    # Plain gap/‖L‖ (relative spectral gap of L)
    println("\n" * "="^110)
    println("RELATIVE SPECTRAL GAP  gap/‖L‖_HS   (higher = better structural — slow mode is well separated)")
    println("="^110)
    @printf("%-3s %-6s | %-9s %-9s %-9s %-9s %-9s %-9s | %-6s %-8s %-9s\n",
            "n", "β_phys",
            "c=0.25", "c=0.5", "c=0.75", "c=1.0", "c=1.5", "c=2.0",
            "c*", "max", "max/min")
    println("-"^110)
    for n in N_RANGE, β in BETA_PHYS_ALL
        vals = Float64[get_metric(rows, n, β, c, :gap_phys) / get_metric(rows, n, β, c, :hs_norm_phys) for c in SIGMA_FACTORS]
        any(isnan, vals) && continue
        imax = argmax(vals)
        spread = maximum(vals) / minimum(vals)
        @printf("%-3d %-6g | %-9.4g %-9.4g %-9.4g %-9.4g %-9.4g %-9.4g | %-6g %-8.4g %-9.3f\n",
                n, β, vals..., SIGMA_FACTORS[imax], vals[imax], spread)
    end

    println("\nReading the c* tables:")
    println("  small c*    → SHARPER kernels give better structural mixing")
    println("  large c*    → WIDER kernels give better structural mixing")
    println("  c* in the middle → there is an interior optimum (tradeoff)")
end

main()
