#!/usr/bin/env julia
# analyze_qf_e4z_39_channel_vs_ideal.jl  (qf-e4z.39 + qf-e4z.37 fix)
#
# Compare the qf-e4z.39 v2 channel CKG sweep (recipe2: per-cell adaptive M_D)
# to the qf-e4z.34 ideal CKG arm at seed=46. Mirrors
# analyze_qf_e4z_36_channel_vs_ideal.jl but reads from the v2 directory.
# Also reports floor_distance and tau_mix_channel_source per cell — every
# cell MUST report :extrapolated for the thesis P5 plot (qf-e4z.38).
#
# qf-e4z.37 fix (2026-05-19): legacy v2 sidecars stored Pass-2 spectral gap as
# `gap_alg_channel`, which is unreliable on this Heisenberg family (7/18 cells
# inflated 11-47% vs ideal). The Pass-1 trajectory gap (matches ideal
# driver's convention) is the comparable quantity. This analyzer:
#   (1) Prefers `gap_alg_pass1_channel` from new sidecars (qf-e4z.37 patch).
#   (2) Falls back to a separately-computed Pass-1 BSON
#       (`qf_e4z_37_pass1_extraction.bson`) if the sidecar lacks Pass-1.
#   (3) Else falls back to `gap_alg_channel` (Pass-2) with a warning column.
# It also reports the dense-LAPACK ground truth gap where available (n ≤ 6).
#
# Usage:
#   julia --project scripts/analyze_qf_e4z_39_channel_vs_ideal.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BSON
using Printf
using QuantumFurnace

const OUT_CH = joinpath(@__DIR__, "output", "sweep_qf_e4z_39_channel_ckg_seed46_v2")
const OUT_ID = joinpath(@__DIR__, "output", "sweep_qf_e4z_34_ckg_vs_dll_plot", "ckg")
const N_RANGE  = 3:8
const BETA_ALL = (0.25, 0.5, 1.0)
const SEED     = 46
const EPS      = 1e-3

function _bp_str(βphys::Real)
    s = let t = @sprintf("%.6f", float(βphys))
        t = rstrip(t, '0'); t = rstrip(t, '.')
        isempty(t) ? "0" : t
    end
    return s
end

_channel_path(n, βphys, seed) = joinpath(OUT_CH,
    "channel_n$(n)_betaphys$(_bp_str(βphys))_seed$(seed)_eps1e-03_smooth_metro_KMS_Trotter.bson")
_ideal_path(n, βphys, seed) = joinpath(OUT_ID,
    "sweep_n$(n)_betaphys$(_bp_str(βphys))_seed$(seed)_L_KMS_Energy.bson")

"""
Load the qf-e4z.37 Pass-1 extraction BSON (if present) and index by (n, β_phys).
Returns Dict((n, β) -> (pass1_gap_phys, dense_gap_phys)).
"""
function _load_pass1_extraction()
    path = joinpath(@__DIR__, "output", "qf_e4z_37_pass1_extraction.bson")
    isfile(path) || return Dict{Tuple{Int, Float64}, NamedTuple}()
    d = BSON.load(path)
    rows = d[:rows]
    out = Dict{Tuple{Int, Float64}, NamedTuple}()
    for r in rows
        key = (Int(r[:n]), float(r[:β_phys]))
        out[key] = (pass1_gap_phys = r[:pass1_gap_phys],
                    pass2_gap_phys = r[:pass2_gap_phys],
                    dense_gap_phys = r[:dense_gap_phys])
    end
    return out
end

function main()
    println("\nqf-e4z.39 v2 channel (recipe2) vs qf-e4z.34 ideal CKG arm @ seed=$SEED")
    println("qf-e4z.37 fix: gap_phys reported is Pass-1 (matches ideal); Pass-2 (legacy) and dense (ground truth, n≤6) shown for audit.")
    println("=" ^ 145)
    @printf("%3s  %6s   %4s   %10s  %10s  %6s   %10s  %10s  %10s   %10s  %10s  %6s   %10s   %10s\n",
            "n", "β_phys",
            "M_D",
            "ch_gap_p1", "id_gap_p", "g_rat",
            "ch_gap_p2", "ch_dense ", "src_used ",
            "ch_τ_p  ", "id_τ_p  ", "τ_rat",
            "ch_floor ", "ch_src ")
    println("=" ^ 145)

    extraction = _load_pass1_extraction()

    rows = NamedTuple[]
    for n in N_RANGE, βphys in BETA_ALL
        ch_p = _channel_path(n, βphys, SEED)
        id_p = _ideal_path(n, βphys, SEED)
        if !isfile(ch_p)
            @printf("%3d  %6.2f   (channel sidecar missing)\n", n, βphys); continue
        end
        if !isfile(id_p)
            @printf("%3d  %6.2f   (ideal sidecar missing)\n", n, βphys); continue
        end
        ch = BSON.load(ch_p, QuantumFurnace)[:result]
        id = BSON.load(id_p, QuantumFurnace)[:result]

        # qf-e4z.37 gap resolution: prefer Pass-1 from updated sidecar; fall back
        # to extraction BSON; finally fall back to Pass-2 (legacy) with warning.
        ext = get(extraction, (Int(n), float(βphys)), nothing)
        ch_gap_p1, src_used = if haskey(ch, :gap_phys_pass1_channel)
            (ch[:gap_phys_pass1_channel], :pass1_sidecar)
        elseif ext !== nothing && isfinite(ext.pass1_gap_phys)
            (ext.pass1_gap_phys, :pass1_extraction)
        else
            (ch[:gap_phys_channel], :pass2_fallback)
        end
        ch_gap_p2 = haskey(ch, :gap_phys_pass2_channel) ? ch[:gap_phys_pass2_channel] : ch[:gap_phys_channel]
        ch_gap_dense = (ext === nothing || !isfinite(ext.dense_gap_phys)) ? NaN : ext.dense_gap_phys

        id_gap_phys = id[:gap_phys]
        ch_tau_phys = ch[:mixing_time_phys_channel]; id_tau_phys = id[:mixing_time_phys]
        g_ratio   = ch_gap_p1 / id_gap_phys
        tau_ratio = ch_tau_phys / id_tau_phys
        M_D       = ch[:M_D]
        @printf("%3d  %6.2f   %4d   %10.4f  %10.4f  %6.3f   %10.4f  %10s   %-10s  %10.3f  %10.3f  %6.3f   %10.3e   %10s\n",
                n, βphys, M_D,
                ch_gap_p1, id_gap_phys, g_ratio,
                ch_gap_p2, isnan(ch_gap_dense) ? "n/a" : @sprintf("%.4f", ch_gap_dense),
                string(src_used),
                ch_tau_phys, id_tau_phys, tau_ratio,
                ch[:floor_distance], string(ch[:tau_mix_channel_source]))
        push!(rows, (n = n, beta_phys = βphys, M_D = M_D,
                     ch_gap_phys = ch_gap_p1, ch_gap_phys_pass2 = ch_gap_p2,
                     ch_gap_phys_dense = ch_gap_dense,
                     id_gap_phys = id_gap_phys, g_ratio = g_ratio,
                     ch_tau_phys = ch_tau_phys, id_tau_phys = id_tau_phys, tau_ratio = tau_ratio,
                     ch_floor = ch[:floor_distance],
                     ch_source = ch[:tau_mix_channel_source],
                     gap_source_used = src_used))
    end

    println("=" ^ 145)
    n_cells = length(rows)
    if n_cells > 0
        # KEY CHECK: every cell must be :extrapolated.
        n_extrap = count(r -> r.ch_source === :extrapolated, rows)
        n_floor  = count(r -> r.ch_source === :floor, rows)
        @printf("\n[verdict] :extrapolated cells: %d / %d  (target: ALL %d)\n",
                n_extrap, n_cells, n_cells)
        @printf("[verdict] :floor cells       : %d / %d  (must be 0 for thesis plot)\n",
                n_floor, n_cells)

        g_ratios = [r.g_ratio for r in rows]
        tau_ratios = [r.tau_ratio for r in rows if isfinite(r.tau_ratio)]
        floors = [r.ch_floor for r in rows]
        @printf("\nSummary (%d cells):\n", n_cells)
        @printf("  gap_ratio  range = [%.3f, %.3f]   median = %.3f\n",
                minimum(g_ratios), maximum(g_ratios), sort(g_ratios)[end÷2 + 1])
        if !isempty(tau_ratios)
            @printf("  tau_ratio  range = [%.3f, %.3f]   median = %.3f\n",
                    minimum(tau_ratios), maximum(tau_ratios), sort(tau_ratios)[end÷2 + 1])
        end
        @printf("  channel-floor range = [%.3e, %.3e]   median = %.3e   (target ε/2 = 5e-4)\n",
                minimum(floors), maximum(floors), sort(floors)[end÷2 + 1])

        # Flag cells with floor above ε/2 (won't bisect cleanly).
        flagged = [r for r in rows if r.ch_floor > EPS / 2]
        if !isempty(flagged)
            println("\nCells with floor > ε/2 (bisection headroom risk):")
            for r in flagged
                @printf("    n=%d β=%g  M_D=%d  floor=%.3e  source=%s\n",
                        r.n, r.beta_phys, r.M_D, r.ch_floor, string(r.ch_source))
            end
        end

        # Flag cells with gap_ratio outside ±10 %.
        flagged_gap = [r for r in rows if r.g_ratio > 1.10 || r.g_ratio < 0.90]
        if !isempty(flagged_gap)
            println("\nCells with gap_ratio outside ±10 %:")
            for r in flagged_gap
                @printf("    n=%d β=%g  g_ratio=%.3f\n", r.n, r.beta_phys, r.g_ratio)
            end
        end
    end

    summary_path = joinpath(@__DIR__, "output", "sweep_qf_e4z_39_channel_vs_ideal_analysis.bson")
    BSON.bson(summary_path, Dict(:rows => Dict.(pairs.(rows))))
    @printf("\nWrote %s\n", summary_path)
    return rows
end

main()
