#!/usr/bin/env julia
# analyze_qf_e4z_33_kdim_sweep.jl  (qf-e4z.33)
#
# Reads the per-cell sidecars from `scratch_qf_e4z_33_kdim_sweep.jl` and
# tabulates:
#   • Per-cell Pass-1 self-saturation:  |gap_p1(k=100) − gap_p1(k=80)| / gap
#   • Per-cell Pass-2 self-saturation:  |gap_p2(k=100) − gap_p2(k=80)| / gap
#   • Cross-agreement at largest kdims: |gap_p1(100) − gap_p2(100)| / gap
#   • Dense ground-truth match (n ≤ 6): |gap_pass − gap_dense| / gap_dense
# Acceptance threshold: rel_err < 1e-6 ("PASS"), else "FAIL".
#
# Also compares against the canonical v6_plus sidecars (Pass 1 kdim=60 +
# Pass 2 kdim=30) — see scripts/output/sweep_S1_v6_plus_ckg_ideal_multiseed/
#
# Usage:  julia --project scripts/analyze_qf_e4z_33_kdim_sweep.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using QuantumFurnace
using BSON
using Printf

const SIDECAR_DIR = joinpath(@__DIR__, "output", "sweep_qf_e4z_33_kdim_sweep", "seed42")
const V6PLUS_DIR  = joinpath(@__DIR__, "output", "sweep_S1_v6_plus_ckg_ideal_multiseed", "smooth_metro_eps1e-03")

const SEED = 42
const BETA_PHYS_ALL = (0.25, 0.5, 1.0, 1.5, 2.0, 2.5)
const N_RANGE = 3:7
const ACCEPT_REL = 1e-6

function _bp_str(β_phys::Real)
    s = @sprintf("%.6f", float(β_phys))
    s = rstrip(s, '0'); s = rstrip(s, '.')
    return isempty(s) ? "0" : s
end

_kdim_sidecar(n, β_phys) =
    joinpath(SIDECAR_DIR,
        "kdim_n$(n)_betaphys$(_bp_str(β_phys))_seed$(SEED)_L_KMS_Energy.bson")

_v6plus_sidecar(n, β_phys) =
    joinpath(V6PLUS_DIR,
        "sweep_n$(n)_betaphys$(_bp_str(β_phys))_seed$(SEED)_L_KMS_Energy.bson")

function _load_sidecar(path::AbstractString)
    isfile(path) || return nothing
    try
        return BSON.load(path, QuantumFurnace)[:result]
    catch err
        @warn "Failed to load $path" err
        return nothing
    end
end

verdict_str(rel) = isnan(rel) ? "  N/A " : (rel < ACCEPT_REL ? " PASS " : "*FAIL*")

function tabulate_dense_ref()
    println("=" ^ 100)
    println("DENSE GROUND-TRUTH MATCH  (n ≤ 6: gap_dense vs gap_p1[k=100], gap_p2[k=100])")
    println("Acceptance: rel_err < $(ACCEPT_REL)")
    println("=" ^ 100)
    @printf("%-3s %-6s %-14s %-14s %-12s %-8s %-14s %-12s %-8s\n",
            "n", "β_phys", "gap_dense", "gap_p1(100)", "rel_p1", "p1?",
            "gap_p2(100)", "rel_p2", "p2?")
    println("-" ^ 100)
    n_fail_p1 = 0; n_fail_p2 = 0; n_cells = 0
    for n in N_RANGE
        for β in BETA_PHYS_ALL
            r = _load_sidecar(_kdim_sidecar(n, β))
            r === nothing && continue
            isnan(r[:gap_dense]) && continue
            gd = r[:gap_dense]
            gp1 = r[:gap_p1_grid][end]
            gp2 = r[:gap_p2_grid][end]
            rel1 = abs(gp1 - gd) / gd
            rel2 = abs(gp2 - gd) / gd
            v1 = verdict_str(rel1); v2 = verdict_str(rel2)
            (v1 == "*FAIL*") && (n_fail_p1 += 1)
            (v2 == "*FAIL*") && (n_fail_p2 += 1)
            n_cells += 1
            @printf("%-3d %-6.2f %-14.6e %-14.6e %-12.2e %-8s %-14.6e %-12.2e %-8s\n",
                    n, β, gd, gp1, rel1, v1, gp2, rel2, v2)
        end
    end
    println("-" ^ 100)
    @printf("Dense cells: %d   P1 failures: %d   P2 failures: %d\n", n_cells, n_fail_p1, n_fail_p2)
    println()
end

function tabulate_self_saturation()
    println("=" ^ 100)
    println("SELF-SATURATION  (relative diff between two largest kdims per pass)")
    println("Acceptance: rel_diff < $(ACCEPT_REL)  → pass is saturated")
    println("=" ^ 100)
    @printf("%-3s %-6s | %-14s %-14s %-12s %-8s | %-14s %-14s %-12s %-8s\n",
            "n", "β_phys",
            "p1[k=80]", "p1[k=100]", "rel", "ok?",
            "p2[k=80]", "p2[k=100]", "rel", "ok?")
    println("-" ^ 110)
    for n in N_RANGE
        for β in BETA_PHYS_ALL
            r = _load_sidecar(_kdim_sidecar(n, β))
            r === nothing && continue
            # Find kdim entries.
            p1g = r[:gap_p1_grid]
            p1k = r[:kdim_p1_grid]
            p2g = r[:gap_p2_grid]
            p2k = r[:kdim_p2_grid]
            i80_1 = findfirst(==(80), p1k); i100_1 = findfirst(==(100), p1k)
            i80_2 = findfirst(==(80), p2k); i100_2 = findfirst(==(100), p2k)
            (i80_1 === nothing || i100_1 === nothing) && continue
            (i80_2 === nothing || i100_2 === nothing) && continue
            rel1 = abs(p1g[i100_1] - p1g[i80_1]) / max(abs(p1g[i100_1]), 1e-30)
            rel2 = abs(p2g[i100_2] - p2g[i80_2]) / max(abs(p2g[i100_2]), 1e-30)
            @printf("%-3d %-6.2f | %-14.6e %-14.6e %-12.2e %-8s | %-14.6e %-14.6e %-12.2e %-8s\n",
                    n, β,
                    p1g[i80_1], p1g[i100_1], rel1, verdict_str(rel1),
                    p2g[i80_2], p2g[i100_2], rel2, verdict_str(rel2))
        end
    end
    println()
end

function tabulate_cross_agreement()
    println("=" ^ 100)
    println("CROSS-AGREEMENT  (Pass 1 vs Pass 2 at largest kdim)")
    println("Acceptance: rel_diff < $(ACCEPT_REL)")
    println("=" ^ 100)
    @printf("%-3s %-6s | %-14s %-14s %-12s %-8s | dense=%s\n",
            "n", "β_phys", "p1[k=100]", "p2[k=100]", "rel_p1_p2", "ok?", "if avail")
    println("-" ^ 100)
    for n in N_RANGE
        for β in BETA_PHYS_ALL
            r = _load_sidecar(_kdim_sidecar(n, β))
            r === nothing && continue
            gp1 = r[:gap_p1_grid][end]
            gp2 = r[:gap_p2_grid][end]
            rel = abs(gp1 - gp2) / max(abs(gp2), 1e-30)
            dense_str = isnan(r[:gap_dense]) ? "(n=7, none)" :
                        @sprintf("%.6e", r[:gap_dense])
            @printf("%-3d %-6.2f | %-14.6e %-14.6e %-12.2e %-8s | %s\n",
                    n, β, gp1, gp2, rel, verdict_str(rel), dense_str)
        end
    end
    println()
end

function tabulate_v6plus_vs_kdim_sweep()
    println("=" ^ 100)
    println("CANONICAL v6_plus (Pass 1 k=60, Pass 2 k=30)  vs  kdim_sweep (P1 k=100, P2 k=100)")
    println("Highlights how much the previous publication-target setup over/under-estimated")
    println("=" ^ 100)
    @printf("%-3s %-6s | v6+ p1[60]    v6+ p2[30]    | kdim p1[100]   kdim p2[100]  | Δp1     Δp2     dense?\n",
            "n", "β_phys")
    println("-" ^ 100)
    for n in N_RANGE
        for β in BETA_PHYS_ALL
            rk = _load_sidecar(_kdim_sidecar(n, β))
            rv = _load_sidecar(_v6plus_sidecar(n, β))
            rk === nothing && continue
            rv === nothing && continue
            gp1_v6 = rv[:gap_arnoldi_pass1]
            gp2_v6 = rv[:gap_arnoldi_pass2]
            gp1_k  = rk[:gap_p1_grid][end]
            gp2_k  = rk[:gap_p2_grid][end]
            d1 = abs(gp1_k - gp1_v6) / max(abs(gp1_k), 1e-30)
            d2 = abs(gp2_k - gp2_v6) / max(abs(gp2_k), 1e-30)
            dense_str = isnan(rk[:gap_dense]) ? "(n7)" :
                        @sprintf("%.4e", rk[:gap_dense])
            @printf("%-3d %-6.2f | %-13.6e %-13.6e | %-13.6e %-13.6e | %-7.1e %-7.1e %s\n",
                    n, β, gp1_v6, gp2_v6, gp1_k, gp2_k, d1, d2, dense_str)
        end
    end
    println()
end

function tabulate_kdim_convergence_n7()
    println("=" ^ 100)
    println("n=7 KDIM CONVERGENCE TRACES  (no dense ref available)")
    println("=" ^ 100)
    for β in BETA_PHYS_ALL
        r = _load_sidecar(_kdim_sidecar(7, β))
        r === nothing && continue
        println("\n--- n=7  β_phys=$β  (β_alg=$(round(r[:beta_alg]; digits=2))  r_D=$(r[:r_D])) ---")
        @printf("  Pass 1 (Arnoldi from |+⟩⟨+|^⊗N):\n")
        for (i, k) in enumerate(r[:kdim_p1_grid])
            @printf("    kdim=%-3d  gap=%.10e  matvecs=%-4d conv=%s wall=%.2fs\n",
                    k, r[:gap_p1_grid][i], r[:matvecs_p1_grid][i],
                    r[:converged_p1_grid][i], r[:wall_p1_each][i])
        end
        @printf("  Pass 2 (krylov_spectral_gap + thick-restart):\n")
        for (i, k) in enumerate(r[:kdim_p2_grid])
            @printf("    kdim=%-3d  gap=%.10e  matvecs=%-4d conv=%d   wall=%.2fs\n",
                    k, r[:gap_p2_grid][i], r[:matvecs_p2_grid][i],
                    r[:converged_p2_grid][i], r[:wall_p2_each][i])
        end
        # Cross diff at max kdim.
        gp1 = r[:gap_p1_grid][end]; gp2 = r[:gap_p2_grid][end]
        @printf("  Cross at max kdim: rel_diff = %.3e\n", abs(gp1 - gp2)/max(abs(gp2), 1e-30))
    end
    println()
end

function main()
    println("\n[analysis] $(now()) — sidecar dir: $SIDECAR_DIR")
    n_files = isdir(SIDECAR_DIR) ? length(readdir(SIDECAR_DIR)) : 0
    println("[analysis] $n_files sidecars present\n")

    tabulate_dense_ref()
    tabulate_self_saturation()
    tabulate_cross_agreement()
    tabulate_kdim_convergence_n7()
    tabulate_v6plus_vs_kdim_sweep()

    println("[analysis] done $(now())")
end

using Dates
main()
