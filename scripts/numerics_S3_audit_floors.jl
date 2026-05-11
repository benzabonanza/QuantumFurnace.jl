#!/usr/bin/env julia
# numerics_S3_audit_floors.jl  (qf-e4z.5)
#
# Audit S3 sidecars (scripts/output/sweep_S3_channel/*.bson):
#   - Summary count of tau_mix_source ∈ {:extrapolated, :floor, :nan}
#   - Per-cell table of :extrapolated cells (USABLE for thesis plots)
#   - Per-cell retune table for :floor / :nan cells with δ_suggest
#   - rm commands to delete only the offending sidecars before a rerun
#
# Usage (run from repo root):
#   julia --project scripts/numerics_S3_audit_floors.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using QuantumFurnace
using BSON
using Printf

const OUTPUT_DIR = joinpath(@__DIR__, "output", "sweep_S3_channel")
const HEADROOM   = 0.3   # target ε / floor on retune (3× headroom below the new floor)

function load_results(dir::AbstractString)
    isdir(dir) || error("Sweep output dir not found: $dir  (run scripts/numerics_sweep_S3.jl first)")
    files = filter(f -> startswith(f, "channel_n") && endswith(f, ".bson"), readdir(dir))
    isempty(files) && error("No channel_n*.bson sidecars in $dir  (run scripts/numerics_sweep_S3.jl first)")
    files = sort(files)
    rows = NamedTuple[]
    for f in files
        d = BSON.load(joinpath(dir, f), QuantumFurnace)
        push!(rows, NamedTuple(d[:result]))
    end
    return rows, files
end

function summarise(rows)
    n_total  = length(rows)
    n_extrap = count(r -> r.tau_mix_source === :extrapolated, rows)
    n_floor  = count(r -> r.tau_mix_source === :floor,        rows)
    n_nan    = count(r -> r.tau_mix_source === :nan,          rows)
    println("\n=== S3 audit summary ($(n_total) cells) ===")
    @printf "  :extrapolated  %3d  (%.1f%%)\n"  n_extrap  100*n_extrap/n_total
    @printf "  :floor         %3d  (%.1f%%)\n"  n_floor   100*n_floor/n_total
    @printf "  :nan           %3d  (%.1f%%)\n"  n_nan     100*n_nan/n_total
    return (n_extrap, n_floor, n_nan)
end

function print_extrapolated_table(rows)
    e = filter(r -> r.tau_mix_source === :extrapolated, rows)
    isempty(e) && return
    sort!(e, by = r -> (r.n, r.beta))
    println("\n=== :extrapolated cells — usable for plots ===")
    println("   n   β     ε       δ         τ_mix         λ_gap        ham_sim_T     wall (s)")
    println("  ────────────────────────────────────────────────────────────────────────────────")
    for r in e
        @printf "  %2d  %4.1f  %.0e  %.2e  %12.3e  %12.3e  %12.3e  %8.1f\n" r.n r.beta r.eps r.delta r.tau_mix r.lambda_gap_channel r.total_ham_sim_time r.wall_time_seconds
    end
end

function print_retune_table(rows, files)
    bad_indices = findall(r -> r.tau_mix_source !== :extrapolated, rows)
    if isempty(bad_indices)
        println("\n✓ All cells :extrapolated — no retune needed. Plots P5/P6/P7 can read this BSON pool directly.")
        return
    end
    println("\n=== Cells requiring retune ($(length(bad_indices))) ===")
    println("  Channel floor exceeds ε for these cells; their tau_mix and total_ham_sim_time")
    println("  encode the gap-bound proxy `log(d/ε)/λ_gap`, NOT an achievable mixing time.")
    println("  Suggested:  δ_new = δ_now · (ε / floor) · $(HEADROOM)   (3× headroom).")
    println()
    println("   n   β     ε       δ_now      floor       ε/floor    λ_gap        δ_suggest    source")
    println("  ──────────────────────────────────────────────────────────────────────────────────────")
    for i in bad_indices
        r = rows[i]
        ratio = r.eps / r.floor_distance
        δ_sug = r.delta * ratio * HEADROOM
        @printf "  %2d  %4.1f  %.0e  %.2e   %.3e   %.3e   %.3e   %.2e     %s\n" r.n r.beta r.eps r.delta r.floor_distance ratio r.lambda_gap_channel δ_sug String(r.tau_mix_source)
    end
    println("\nTo retune: edit the δ entry for the flagged (n, β, ε, filter) rows in")
    println("scripts/numerics_param_table.jl, regenerate channel_param_table.bson, then")
    println("delete ONLY the flagged sidecars and rerun the sweep:")
    println()
    for i in bad_indices
        println("  rm '$(joinpath(OUTPUT_DIR, files[i]))'")
    end
    println()
    println("  julia --project scripts/numerics_sweep_S3.jl")
    println("  julia --project scripts/numerics_S3_audit_floors.jl    # confirm all :extrapolated")
end

# --- Main -----------------------------------------------------------------

rows, files = load_results(OUTPUT_DIR)
summarise(rows)
print_extrapolated_table(rows)
print_retune_table(rows, files)
