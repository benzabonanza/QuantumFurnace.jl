#!/usr/bin/env julia
# analyze_qf_c9g_verdict.jl
#
# After the qf-c9g main sweep + slow-mode characterization complete, this
# script reads the BSON sidecars and produces the verdict-checklist table
# that maps onto §3 of `drafts/qf-c9g-2d-tfim-ordered-mechanism.md`.
#
# Runs:
#   1. n × β gap-closing table — does gap_phys close with n at every β below T_c?
#   2. Per-β n-slope (log gap_phys vs n) — quantitative closing rate.
#   3. Paramagnetic control β ∈ {0.10, 0.25} — gap should be Ω(1) flat.
#   4. Magnetization-indicator pivot (U_4, eff_rank) — ordered/paramagnet
#      character at each (n, β).
#   5. Slow-mode signature (from `qf_c9g_slow_mode_n4.bson`) — doublet × bulk
#      mass at each β.
# Then a final "verdict checklist" table mapping observations → diagnostic
# passes/fails.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BSON
using Printf
using Statistics

const BSON_DIR = joinpath(@__DIR__, "output", "qf_c9g_ordered_gap_mechanism")
const T_C = 2.07
const BETA_C = 1.0 / T_C   # ≈ 0.483

function load_main_sweep()
    path = joinpath(BSON_DIR, "qf_c9g_beta_n_sweep.bson")
    isfile(path) || error("Main sweep BSON not found: $path")
    raw = BSON.load(path)
    rows = [NamedTuple{Tuple(Symbol(k) for k in keys(d))}(values(d)) for d in raw[:sweep_a_rows]]
    sweep_b = [NamedTuple{Tuple(Symbol(k) for k in keys(d))}(values(d)) for d in raw[:sweep_b_rows]]
    return (rows = rows, sweep_b = sweep_b,
            n_list = raw[:n_list], beta_grid = raw[:beta_grid])
end

function load_slow_mode()
    path = joinpath(BSON_DIR, "qf_c9g_slow_mode_n4.bson")
    isfile(path) || return nothing
    raw = BSON.load(path)
    rows = [NamedTuple{Tuple(Symbol(k) for k in keys(d))}(values(d)) for d in raw[:rows]]
    return rows
end

function print_main_table(data)
    println("\n" * "="^130)
    println("VERDICT TABLE 1 — gap_phys(n, β) on the full 15-cell grid")
    println("="^130)
    @printf("%-3s %-7s %-7s %-9s %-12s %-12s %-9s %-9s %-9s %-9s %-9s\n",
        "n", "β_phys", "T/T_c", "phase", "λ_L^phys", "ΔE_1^phys",
        "λ_L/ΔE_1", "λ3/λ2", "<m_z²>", "U_4", "eff_rank")
    println("-"^130)
    for β in sort(data.beta_grid)
        for n in sort(data.n_list)
            row = first(r for r in data.rows if r.n == n && isapprox(r.beta_phys, β; atol=1e-12))
            phase = (1/β) < T_C ? "ORDERED" : "paramagn"
            @printf("%-3d %-7.3g %-7.3g %-9s %-12.4e %-12.4e %-9.3g %-9.3g %-9.4f %-9.4f %-9.3g\n",
                n, β, (1/β)/T_C, phase,
                row.gap_phys, row.Δ1_phys, row.ratio_gap_dE1, row.λ3_over_λ2,
                row.m2, row.U4, row.eff_rank)
        end
        println()
    end
end

"Linear regression of log(gap_phys) vs n. Returns (slope, intercept, R²)."
function log_n_slope(ns::Vector{Int}, gaps::Vector{Float64})
    xs = Float64.(ns)
    ys = log.(gaps)
    n = length(xs)
    x̄ = mean(xs); ȳ = mean(ys)
    sxx = sum((xs .- x̄).^2)
    sxy = sum((xs .- x̄) .* (ys .- ȳ))
    syy = sum((ys .- ȳ).^2)
    slope = sxy / sxx
    intercept = ȳ - slope * x̄
    r2 = sxy^2 / (sxx * syy)
    return (slope = slope, intercept = intercept, r2 = r2)
end

function print_slopes(data)
    println("="^110)
    println("VERDICT TABLE 2 — log(λ_L^phys) vs n slope at each β (negative = closing; flat ≈ no closing)")
    println("="^110)
    @printf("%-7s %-7s %-9s %-10s %-10s %-8s %-15s\n",
        "β_phys", "T/T_c", "phase", "slope", "exp(slope)", "R²", "verdict")
    println("-"^110)
    for β in sort(data.beta_grid)
        ord_rows = sort([r for r in data.rows if isapprox(r.beta_phys, β; atol=1e-12)], by = r->r.n)
        ns = [Int(r.n) for r in ord_rows]
        gaps = [r.gap_phys for r in ord_rows]
        fit = log_n_slope(ns, gaps)
        phase = (1/β) < T_C ? "ORDERED" : "paramagn"
        verdict = if abs(fit.slope) < 0.1
            "FLAT — no closing"
        elseif fit.slope < -0.5
            "STRONG closing"
        elseif fit.slope < -0.1
            "moderate closing"
        else
            "(growing??)"
        end
        @printf("%-7.3g %-7.3g %-9s %-+10.4f %-10.4f %-8.4f %-15s\n",
            β, (1/β)/T_C, phase,
            fit.slope, exp(fit.slope),
            fit.r2, verdict)
    end
end

function print_slow_mode(slow_rows)
    if slow_rows === nothing
        println("\n[SWEEP C slow-mode characterisation not yet run — skipping]")
        return
    end
    println("\n" * "="^120)
    println("VERDICT TABLE 3 — slow-mode operator-space character at n=4 (dense L eigendecomp)")
    println("="^120)
    @printf("%-7s %-7s %-9s %-12s %-12s %-12s %-12s %-12s %-12s\n",
        "β_phys", "T/T_c", "phase", "|λ_2|^phys",
        "M(d×d)", "M(d×bulk)", "M(b×b)", "Z₂ parity", "⟨R₂,M_z⟩")
    println("-"^120)
    for r in sort(slow_rows, by = x->x.beta_phys)
        phase = r.T_over_Tc < 1.0 ? "ORDERED" : "paramagn"
        @printf("%-7.3g %-7.3g %-9s %-12.4e %-12.4f %-12.4f %-12.4f %-12.3g %-12.3g\n",
            r.beta_phys, r.T_over_Tc, phase,
            abs(r.λ_2_phys),
            r.mass_doublet_doublet, r.mass_doublet_bulk, r.mass_bulk_bulk,
            r.z2_parity, r.overlap_Mz)
    end
end

function print_verdict_checklist(data, slow_rows)
    println("\n" * "="^110)
    println("FINAL VERDICT CHECKLIST")
    println("="^110)

    # Test (i) — β-persistence: gap_phys must close at every ordered β
    ordered_betas = [β for β in data.beta_grid if (1/β) < T_C]
    closing_at = String[]
    flat_at = String[]
    for β in ordered_betas
        ord_rows = sort([r for r in data.rows if isapprox(r.beta_phys, β; atol=1e-12)], by = r->r.n)
        ns = [Int(r.n) for r in ord_rows]
        gaps = [r.gap_phys for r in ord_rows]
        fit = log_n_slope(ns, gaps)
        if fit.slope < -0.2
            push!(closing_at, @sprintf("β=%.2g (slope=%+.2f)", β, fit.slope))
        else
            push!(flat_at, @sprintf("β=%.2g (slope=%+.2f)", β, fit.slope))
        end
    end
    test_i = isempty(flat_at) ? "✓ PASS" : "✗ FAIL"
    println("(i) Closing persists across all β below T_c:        $test_i")
    println("    closing at: " * join(closing_at, ", "))
    isempty(flat_at) || println("    FLAT at:    " * join(flat_at, ", "))

    # Test (ii) — Binder U_4 saturates near 2/3 across ordered cells
    ord_u4 = Float64[]
    for β in ordered_betas, n in data.n_list
        row = first(r for r in data.rows if r.n == n && isapprox(r.beta_phys, β; atol=1e-12))
        push!(ord_u4, row.U4)
    end
    u4_min, u4_max = extrema(ord_u4)
    test_ii = (u4_min > 0.5) ? "✓ PASS" : "△ partial"
    @printf("(ii) U_4 in ordered window in [%.3f, %.3f]:           %s\n",
        u4_min, u4_max, test_ii)

    # Test (iii) — slow mode in doublet × bulk block
    if slow_rows !== nothing
        ord_slow = [r for r in slow_rows if r.T_over_Tc < 1.0]
        if isempty(ord_slow)
            println("(iii) [no ordered β in slow-mode rows]")
        else
            min_db_bulk = minimum(r.mass_doublet_bulk for r in ord_slow)
            test_iii = (min_db_bulk > 0.5) ? "✓ PASS" : "△ partial"
            @printf("(iii) doublet × bulk mass ≥ 0.5 in all ord β (min %.3f): %s\n",
                min_db_bulk, test_iii)
        end
    else
        println("(iii) [SWEEP C slow-mode not yet run — skipped]")
    end

    # Test (iv) — paramagnetic control: gap flat
    para_betas = [β for β in data.beta_grid if (1/β) >= T_C]
    flat_para = String[]; closing_para = String[]
    for β in para_betas
        ord_rows = sort([r for r in data.rows if isapprox(r.beta_phys, β; atol=1e-12)], by = r->r.n)
        ns = [Int(r.n) for r in ord_rows]
        gaps = [r.gap_phys for r in ord_rows]
        fit = log_n_slope(ns, gaps)
        if abs(fit.slope) < 0.2
            push!(flat_para, @sprintf("β=%.2g (slope=%+.2f)", β, fit.slope))
        else
            push!(closing_para, @sprintf("β=%.2g (slope=%+.2f)", β, fit.slope))
        end
    end
    test_iv = isempty(closing_para) ? "✓ PASS" : "✗ FAIL"
    println("(iv) Paramagnetic control: gap flat in n:           $test_iv")
    println("    flat at: " * join(flat_para, ", "))
    isempty(closing_para) || println("    CLOSING at: " * join(closing_para, ", "))

    println()
    println("="^110)
end

function main()
    data = load_main_sweep()
    slow_rows = load_slow_mode()

    print_main_table(data)
    print_slopes(data)
    print_slow_mode(slow_rows)
    print_verdict_checklist(data, slow_rows)
end

main()
