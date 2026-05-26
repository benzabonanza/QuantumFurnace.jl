#!/usr/bin/env julia
#
# qf-2dn: σ-baseline-normalised gap ratio WITHIN CKG smooth-Metropolis.
# Within-sampler analogue of qf-79h (CKG vs DLL) and qf-cwg (KMS vs GNS).
# Same machinery, different question: vary the smoothing factor σ_f at
# fixed (n, β) and ask whether the d_{1→1} norm captures the σ-induced
# rate variation, or whether a STRUCTURAL piece (kink-window width and
# similar effects) remains.
#
# Per-cell ratio (σ_f=1 is the baseline ≡ qf-e4z.34 canonical operating
# point):
#
#     R_σ(n, β; σ_f) =  (gap_phys / d11_phys) |_σ_f
#                     ─────────────────────────────────
#                       (gap_phys / d11_phys) |_{σ_f=1}
#
# At σ_f = 1, R_σ ≡ 1 by construction. The CENTRAL QUESTION is whether
# R_σ(n; σ_f) DRIFTS WITH n at fixed (β, σ_f):
#   • FLAT in n   → σ-dependence is rate-only as captured by d11 within ε;
#                   d11 is a clean rate canceller for σ.
#   • DRIFTS in n → σ-dependence has a STRUCTURAL piece d11 misses
#                   (kink-window width, broadening of the Metropolis
#                   acceptance shape …). Quantify its size vs the rate
#                   piece.
#
# # PHYSICS CHECK: d_{1→1} is the operator norm of the Kossakowski matrix
# M[k,j] = Σ_i conj(A[i,k]) A[i,j] α(ν_{ij}, ν_{ik}) (see
# scripts/scratch_qf_e4z_35_sigma_sweep_plot.jl:294-298). For fixed jump
# operators A_i and a single sampler, varying σ_f changes ONLY the α
# kernel through the kink width σ_alg = c · σ_β and its smoothing factor
# s. Under a uniform time-scale L → kL the ratio gap/d11 is invariant, so
# any residual σ-dependence of gap/d11 is necessarily NON-rate (i.e.
# structural — γ(ω) shape, not amplitude). The same-cell baseline ratio
# isolates that structural piece.
#
# Data (no new sweeps; all on disk):
#   scripts/output/sweep_qf_e4z_35_sigma_sweep_plot/ckg/
#     sweep_n{n}_betaphys{β}_sigma{σ_f}_seed46_L_KMS_Energy.bson
#   Grid: n ∈ {3..7} × β_phys ∈ {0.25, 0.5, 1.0} × σ_f ∈ {0.25, 0.5,
#         0.75, 1.0, 1.5, 2.0}, seed = 46 only (qf-e4z.35 is single-seed).
#   Plus 4 partial cells at (n=8, β_phys=0.25, σ_f ∈ {0.25, 0.5, 0.75, 1})
#   used as a small-N=8 footnote (NOT in the main fit).
#
# Output (mirror qf-79h / qf-cwg layout):
#   drafts/figures/numerics/qf_2dn_ckg_sigma_ratio_main.{png,pdf}
#     3-panel β=0.5: (a) raw gap_phys vs n at each σ_f,
#                    (b) raw τ_mix_phys vs n at each σ_f,
#                    (c) R_σ(n) vs n, one curve per σ_f≠1, hline at 1.
#   drafts/figures/numerics/qf_2dn_ckg_sigma_ratio_checks.{png,pdf}
#     3-panel: (a) β-stability of R_σ at σ_f=2,
#              (b) cross-norm — R_σ(d11) vs R_σ(HS) at β=0.5,
#              (c) all-β collapse at σ_f=0.25 — does R_σ(n) coincide
#                  across β_phys ∈ {0.25, 0.5, 1.0}?
#   drafts/figures/numerics/qf_2dn_ckg_sigma_ratio_data.bson
#
# NB: per qf-2dn issue notes, the only seed=42 cell with d11/hs/gap on
# disk is qf_e4z_34_norm_diagnostic.bson at σ_f=1 (the baseline itself,
# R_σ ≡ 1 trivially) — so a seed cross-check is NOT available without a
# new ~90-cell sweep. The issue explicitly drops it from the deliverable;
# seed=46 throughout.

using Printf, Statistics, BSON, Plots
ENV["GKSwstype"] = "100"

const REPO_ROOT = abspath(joinpath(@__DIR__, ".."))
const SIGMA_DIR = joinpath(REPO_ROOT, "scripts", "output",
                           "sweep_qf_e4z_35_sigma_sweep_plot", "ckg")
const OUT_DIR   = joinpath(REPO_ROOT, "drafts", "figures", "numerics")

const SEED          = 46
const BASELINE_σF   = 1.0
const HEADLINE_β    = 0.5             # main panels use β_phys = 0.5
const SIGMA_GRID    = (0.25, 0.5, 0.75, 1.0, 1.5, 2.0)
const BETAS         = (0.25, 0.5, 1.0)
const NS_MAIN       = 3:7             # fit window for the verdict is n=4..7
const N_FOOTNOTE    = 8               # n=8 β=0.25 partial coverage footnote

# Colours.  σ_f curves use a cool→neutral→warm sequence (purples below 1,
# warms above), matching the thesis palette (reference_thesis_colors.md).
const COL_BASELINE = RGB(0.55, 0.55, 0.55)
const SIGMA_COL = Dict(
    0.25 => RGB(0.310, 0.231, 0.361),  # aubergine    — deepest cool
    0.5  => RGB(0.451, 0.345, 0.455),  # deepplum
    0.75 => RGB(0.557, 0.435, 0.549),  # dustyplum
    1.0  => COL_BASELINE,
    1.5  => RGB(0.722, 0.569, 0.263),  # ochre        — warm
    2.0  => RGB(0.710, 0.396, 0.290),  # terracotta   — warmest
)
# β colours match qf-79h / qf-cwg conventions
const BETA_COL = Dict(
    0.25 => RGB(0.373, 0.545, 0.557),  # dustyteal
    0.5  => RGB(0.710, 0.396, 0.290),  # terracotta
    1.0  => RGB(0.310, 0.231, 0.361),  # aubergine
)

# Σ_f marker shapes — keeps the curves distinguishable in B/W too.
const SIGMA_MARKER = Dict(
    0.25 => :diamond, 0.5 => :utriangle, 0.75 => :circle,
    1.0  => :star5,   1.5 => :pentagon,  2.0 => :square,
)

# ── Loader ──────────────────────────────────────────────────────────────────

# Same filename convention as qf-cwg/qf-e4z.35 driver: %f formatted then
# trailing zeros + dot stripped.
function _sf_str(σf::Real)
    s = @sprintf("%.6f", float(σf))
    s = rstrip(s, '0'); s = rstrip(s, '.')
    return isempty(s) ? "0" : s
end

function load_cell(n::Integer, β_phys::Real, σf::Real)
    β_str = β_phys == 1.0 ? "1" : (β_phys == 0.5 ? "0.5" : string(β_phys))
    f = joinpath(SIGMA_DIR,
        "sweep_n$(n)_betaphys$(β_str)_sigma$(_sf_str(σf))_seed$(SEED)_L_KMS_Energy.bson")
    isfile(f) || return nothing
    r = BSON.load(f)[:result]
    return (
        n = r[:n], beta_phys = r[:beta_phys], beta_alg = r[:beta_alg],
        sigma_factor = r[:sigma_factor], sigma_alg = r[:sigma_alg],
        sigma_phys = r[:sigma_phys], seed = r[:seed],
        gap_phys = r[:gap_phys], gap_alg = r[:gap_alg],
        tau_phys = r[:mixing_time_phys], tau_src = r[:mixing_time_source],
        d11_phys = r[:d_1to1_phys], d11_alg = r[:d_1to1_alg],
        hs_phys  = r[:hs_norm_phys], hs_alg  = r[:hs_norm_alg],
        all_converged = r[:all_converged], floor_distance = r[:floor_distance],
        rescaling_factor = r[:rescaling_factor],
    )
end

# Returns (cells = Dict{(n, β, σ_f) => row}, missing = Vector{(n, β, σ_f)})
function load_all()
    cells = Dict{Tuple{Int, Float64, Float64}, NamedTuple}()
    miss  = Tuple{Int, Float64, Float64}[]
    for n in NS_MAIN, β in BETAS, σf in SIGMA_GRID
        r = load_cell(n, β, σf)
        if r === nothing
            push!(miss, (n, β, σf))
        else
            cells[(n, Float64(β), Float64(σf))] = r
        end
    end
    # n=8 footnote: only β=0.25, σ_f ∈ {0.25, 0.5, 0.75, 1.0} on disk
    for σf in (0.25, 0.5, 0.75, 1.0)
        r = load_cell(N_FOOTNOTE, 0.25, σf)
        r === nothing || (cells[(N_FOOTNOTE, 0.25, Float64(σf))] = r)
    end
    return (cells = cells, missing = miss)
end

# ── Ratios ──────────────────────────────────────────────────────────────────

# Compute R_σ(n, β; σ_f) and the HS-norm twin from baseline at σ_f = 1.
# Returns nothing if either the σ_f cell or the baseline cell is missing.
function ratio_at(cells, n::Integer, β::Real, σf::Real;
                  baseline_σf::Real = BASELINE_σF)
    cur = get(cells, (n, Float64(β), Float64(σf)),    nothing)
    bse = get(cells, (n, Float64(β), Float64(baseline_σf)), nothing)
    (cur === nothing || bse === nothing) && return nothing
    R_d11 = (cur.gap_phys / cur.d11_phys) / (bse.gap_phys / bse.d11_phys)
    R_hs  = (cur.gap_phys / cur.hs_phys)  / (bse.gap_phys / bse.hs_phys)
    R_tau = bse.tau_phys / cur.tau_phys   # >1 means current σ_f mixes faster
    return (n = n, beta_phys = Float64(β), sigma_factor = Float64(σf),
            R_d11 = R_d11, R_hs = R_hs, R_tau = R_tau,
            gap_phys = cur.gap_phys, tau_phys = cur.tau_phys,
            d11_phys = cur.d11_phys, hs_phys = cur.hs_phys,
            sigma_alg = cur.sigma_alg, sigma_phys = cur.sigma_phys)
end

# ── Pretty-print ────────────────────────────────────────────────────────────

function print_main_table(cells)
    println()
    println("=== Per-σ_f gap, τ, d11 at β_phys=$HEADLINE_β (seed=$SEED) ===")
    for σf in SIGMA_GRID
        println("\nσ_f = $σf  (",
                σf == BASELINE_σF ? "BASELINE" :
                @sprintf("R_σ vs baseline σ_f=%.1f", BASELINE_σF), ")")
        @printf("%4s | %8s %8s %8s %8s | %8s %8s %8s\n",
            "n", "gap_phys", "τ_phys", "d11_phys", "hs_phys",
            "R_d11", "R_hs", "R_tau")
        println(repeat("-", 84))
        for n in NS_MAIN
            r = ratio_at(cells, n, HEADLINE_β, σf)
            r === nothing && continue
            @printf("%4d | %8.4f %8.4f %8.3f %8.3f | %8.4f %8.4f %8.4f\n",
                r.n, r.gap_phys, r.tau_phys, r.d11_phys, r.hs_phys,
                r.R_d11, r.R_hs, r.R_tau)
        end
    end
end

function print_drift_summary(cells)
    println()
    println("=== Drift of R_σ over n=4..7 at β_phys=$HEADLINE_β ===")
    println("(rel. range = (max - min)/mean across n; FLAT ≲ 3 % ⇒ d11 cancels σ)")
    @printf("%-8s | %8s %8s %8s %8s %8s | %10s\n",
        "σ_f", "n=4", "n=5", "n=6", "n=7", "mean", "rel.range")
    println(repeat("-", 78))
    for σf in SIGMA_GRID
        σf == BASELINE_σF && continue
        Rs = Float64[]
        for n in 4:7
            r = ratio_at(cells, n, HEADLINE_β, σf)
            r === nothing || push!(Rs, r.R_d11)
        end
        isempty(Rs) && continue
        rng = (maximum(Rs) - minimum(Rs)) / mean(Rs)
        @printf("%-8.2f | ", σf)
        for v in Rs;  @printf("%8.4f ", v); end
        @printf("| %10.4f\n", rng)
    end
end

function print_n3_anomaly(cells)
    println()
    println("=== n=3 anomaly check at β_phys=$HEADLINE_β (compare R_σ(n=3) vs mean R_σ at n=4..7) ===")
    @printf("%-8s | %10s %10s %10s\n", "σ_f", "R_σ(n=3)", "mean(n=4..7)", "rel.diff")
    println(repeat("-", 50))
    for σf in SIGMA_GRID
        σf == BASELINE_σF && continue
        r3 = ratio_at(cells, 3, HEADLINE_β, σf)
        r3 === nothing && continue
        means = Float64[]
        for n in 4:7
            r = ratio_at(cells, n, HEADLINE_β, σf)
            r === nothing || push!(means, r.R_d11)
        end
        isempty(means) && continue
        m = mean(means)
        @printf("%-8.2f | %10.4f %10.4f %10.4f\n", σf, r3.R_d11, m, abs(r3.R_d11 - m) / m)
    end
end

function print_beta_stability(cells; σf::Real = 2.0)
    println()
    println("=== β-stability of R_σ at σ_f=$σf (does R_σ(n) collapse across β?) ===")
    @printf("%4s | %12s %12s %12s\n", "n",
        "β=0.25", "β=0.50", "β=1.00")
    println(repeat("-", 50))
    for n in NS_MAIN
        vals = String[]
        for β in BETAS
            r = ratio_at(cells, n, β, σf)
            push!(vals, r === nothing ? "    —    " : @sprintf("%10.4f  ", r.R_d11))
        end
        @printf("%4d | %s\n", n, join(vals, " "))
    end
end

function print_cross_norm(cells; β::Real = HEADLINE_β)
    println()
    println("=== Cross-norm: R_σ(d11) vs R_σ(HS) at β_phys=$β  ===")
    @printf("%-8s | %4s %12s %12s %10s\n", "σ_f", "n", "R_d11", "R_hs",
            "|ΔR|/R_d11")
    println(repeat("-", 60))
    for σf in SIGMA_GRID
        σf == BASELINE_σF && continue
        for n in NS_MAIN
            r = ratio_at(cells, n, β, σf)
            r === nothing && continue
            δ = abs(r.R_d11 - r.R_hs) / abs(r.R_d11)
            @printf("%-8.2f | %4d %12.5f %12.5f %10.5f\n", σf, n, r.R_d11, r.R_hs, δ)
        end
    end
end

function print_n8_footnote(cells)
    println()
    println("=== n=8 footnote (partial coverage at β_phys=0.25 only) ===")
    @printf("%-8s | %10s %10s %10s\n", "σ_f", "R_σ(n=8)",
            "mean(n=4..7,β=0.25)", "rel.diff")
    println(repeat("-", 50))
    for σf in (0.25, 0.5, 0.75)
        r8 = ratio_at(cells, N_FOOTNOTE, 0.25, σf)
        r8 === nothing && continue
        means = Float64[]
        for n in 4:7
            r = ratio_at(cells, n, 0.25, σf)
            r === nothing || push!(means, r.R_d11)
        end
        isempty(means) && continue
        m = mean(means)
        @printf("%-8.2f | %10.4f %10.4f %10.4f\n", σf, r8.R_d11, m, (r8.R_d11 - m) / m)
    end
end

# ── Plots ───────────────────────────────────────────────────────────────────

# Panel-level helpers.

# Sort (n, σ_f) cells into per-σ_f vectors (skipping missing cells).
function _series(cells, β::Real, σf::Real, getter::Function, ns = NS_MAIN)
    xs = Int[]; ys = Float64[]
    for n in ns
        r = get(cells, (n, Float64(β), Float64(σf)), nothing)
        r === nothing && continue
        push!(xs, n); push!(ys, getter(r))
    end
    return xs, ys
end

# Per-σ_f R-vector (excludes baseline).
function _R_series(cells, β::Real, σf::Real, which::Symbol = :R_d11,
                   ns = NS_MAIN)
    xs = Int[]; ys = Float64[]
    for n in ns
        r = ratio_at(cells, n, β, σf)
        r === nothing && continue
        push!(xs, n); push!(ys, getproperty(r, which))
    end
    return xs, ys
end

function plot_main(cells; out_path::String)
    β = HEADLINE_β

    # ── Panel A: raw gap_phys vs n, one curve per σ_f
    pA = plot(xlabel="n", ylabel="gap_phys", legend=:topright,
              title="(a) Raw spectral gap (β_phys=$β)", titlefontsize=10,
              size=(380, 320), framestyle=:box, legendfontsize=7,
              legend_background_color=RGBA(1,1,1,0.7))
    for σf in SIGMA_GRID
        xs, ys = _series(cells, β, σf, r -> r.gap_phys)
        isempty(xs) && continue
        lbl = σf == BASELINE_σF ? "σ_f=$σf (baseline)" : "σ_f=$σf"
        ls  = σf == BASELINE_σF ? :solid : :solid
        lw  = σf == BASELINE_σF ? 2.5     : 1.6
        plot!(pA, xs, ys, color=SIGMA_COL[σf], lw=lw, ls=ls,
              marker=SIGMA_MARKER[σf], ms=4, label=lbl)
    end

    # ── Panel B: raw τ_mix vs n, one curve per σ_f
    pB = plot(xlabel="n", ylabel="τ_mix_phys", legend=:topleft,
              title="(b) Raw mixing time (β_phys=$β)", titlefontsize=10,
              size=(380, 320), framestyle=:box, legendfontsize=7,
              legend_background_color=RGBA(1,1,1,0.7))
    for σf in SIGMA_GRID
        xs, ys = _series(cells, β, σf, r -> r.tau_phys)
        isempty(xs) && continue
        lw  = σf == BASELINE_σF ? 2.5 : 1.6
        plot!(pB, xs, ys, color=SIGMA_COL[σf], lw=lw,
              marker=SIGMA_MARKER[σf], ms=4, label="σ_f=$σf")
    end

    # ── Panel C: R_σ(n) at each σ_f ≠ 1, with hline at 1 for baseline.
    # The visual punchline: if all curves stay near 1 across n, d11 cancels
    # the σ-induced rate variation cleanly; if they drift, structural piece.
    pC = plot(xlabel="n",
              ylabel="R_σ(n; σ_f)  =  (gap/d11)_σ_f / (gap/d11)_{σ_f=1}",
              legend=:topright,
              title="(c) σ-baseline ratio at β_phys=$β", titlefontsize=10,
              size=(380, 320), framestyle=:box, legendfontsize=7,
              legend_background_color=RGBA(1,1,1,0.7),
              ylabelfontsize=7)
    hline!(pC, [1.0], ls=:dot, color=:gray, label="ratio = 1 (baseline)")
    for σf in SIGMA_GRID
        σf == BASELINE_σF && continue
        xs, ys = _R_series(cells, β, σf, :R_d11)
        isempty(xs) && continue
        plot!(pC, xs, ys, color=SIGMA_COL[σf], lw=2.2,
              marker=SIGMA_MARKER[σf], ms=4.5, label="σ_f=$σf")
    end

    fig = plot(pA, pB, pC, layout=(1, 3), size=(1140, 320),
               left_margin=5Plots.mm, bottom_margin=5Plots.mm,
               top_margin=3Plots.mm)
    savefig(fig, out_path * ".png")
    savefig(fig, out_path * ".pdf")
    println("\nMain figure saved: ", out_path, ".{png,pdf}")
    return fig
end

function plot_checks(cells; out_path::String, check_σf_a::Real = 2.0,
                     check_σf_c::Real = 0.25, check_β_b::Real = HEADLINE_β)

    # ── Check A: β-stability of R_σ at one σ_f. If R_σ(n) at fixed σ_f
    # depends on β, the σ-effect interacts with temperature (each β has
    # its own σ_alg = c/β_alg at the baseline, so β-collapse is non-trivial
    # evidence that d11 cleanly factors out the rate part of σ).
    pA = plot(xlabel="n", ylabel="R_σ at σ_f=$check_σf_a",
              legend=:topright,
              title="(a) β-stability of R_σ", titlefontsize=10,
              size=(380, 320), framestyle=:box, legendfontsize=7,
              legend_background_color=RGBA(1,1,1,0.7))
    hline!(pA, [1.0], ls=:dot, color=:gray, label=nothing)
    for β in BETAS
        xs, ys = _R_series(cells, β, check_σf_a, :R_d11)
        isempty(xs) && continue
        plot!(pA, xs, ys, color=BETA_COL[β], lw=2, marker=:circle, ms=4,
              label="β_phys=$β")
    end

    # ── Check B: cross-norm — R_σ(d11) vs R_σ(HS) across σ_f at β_phys=0.5,
    # n on the x-axis grouped by σ_f. If R(d11) and R(HS) coincide, the
    # norm choice is robust (qf-79h / qf-cwg precedent).
    pB = plot(xlabel="n", ylabel="R_σ at β_phys=$check_β_b",
              legend=:topright,
              title="(b) Cross-norm  d_{1→1} vs HS",
              titlefontsize=10, size=(380, 320), framestyle=:box,
              legendfontsize=6,
              legend_background_color=RGBA(1,1,1,0.7))
    hline!(pB, [1.0], ls=:dot, color=:gray, label=nothing)
    for σf in SIGMA_GRID
        σf == BASELINE_σF && continue
        xs_d, ys_d = _R_series(cells, check_β_b, σf, :R_d11)
        xs_h, ys_h = _R_series(cells, check_β_b, σf, :R_hs)
        isempty(xs_d) && continue
        plot!(pB, xs_d, ys_d, color=SIGMA_COL[σf], lw=2,
              marker=SIGMA_MARKER[σf], ms=4,
              label="d11 σ_f=$σf")
        plot!(pB, xs_h, ys_h, color=SIGMA_COL[σf], lw=1.2, ls=:dash,
              marker=SIGMA_MARKER[σf], ms=3, alpha=0.7,
              label="HS  σ_f=$σf")
    end

    # ── Check C: all-β collapse at one σ_f (small-σ extreme: σ_f=0.25).
    # Same data as panel A but at the OPPOSITE end of the σ-range — if the
    # collapse holds at both σ_f extremes, β-independence of R_σ(n) is
    # robust across the smoothing-window range.
    pC = plot(xlabel="n", ylabel="R_σ at σ_f=$check_σf_c",
              legend=:topright,
              title="(c) β-collapse at small σ_f", titlefontsize=10,
              size=(380, 320), framestyle=:box, legendfontsize=7,
              legend_background_color=RGBA(1,1,1,0.7))
    hline!(pC, [1.0], ls=:dot, color=:gray, label=nothing)
    for β in BETAS
        xs, ys = _R_series(cells, β, check_σf_c, :R_d11)
        isempty(xs) && continue
        plot!(pC, xs, ys, color=BETA_COL[β], lw=2, marker=:circle, ms=4,
              label="β_phys=$β")
    end

    fig = plot(pA, pB, pC, layout=(1, 3), size=(1140, 320),
               left_margin=5Plots.mm, bottom_margin=5Plots.mm,
               top_margin=3Plots.mm)
    savefig(fig, out_path * ".png")
    savefig(fig, out_path * ".pdf")
    println("Checks figure saved: ", out_path, ".{png,pdf}")
    return fig
end

# ── Main ────────────────────────────────────────────────────────────────────

function main()
    isdir(OUT_DIR) || mkpath(OUT_DIR)

    println("Loading qf-e4z.35 CKG σ-sweep (seed=$SEED) …")
    loaded = load_all()
    cells = loaded.cells
    @assert !isempty(cells) "no cells loaded — check $SIGMA_DIR"
    println("  → loaded $(length(cells)) cells; missing $(length(loaded.missing)) of the n∈$(collect(NS_MAIN))×β×σ_f grid.")
    if !isempty(loaded.missing)
        for m in loaded.missing
            println("    missing: n=$(m[1]) β=$(m[2]) σ_f=$(m[3])")
        end
    end

    # Convergence sanity: every cell on the headline grid converged.
    bad = NamedTuple[]
    for ((n, β, σf), r) in cells
        r.all_converged || push!(bad, (n=n, β=β, σf=σf, floor=r.floor_distance))
    end
    if !isempty(bad)
        println("\nWARNING — non-converged cells:")
        for b in bad
            @printf("    n=%d β=%.2f σ_f=%.2f floor_distance=%.3e\n",
                    b.n, b.β, b.σf, b.floor)
        end
    end

    print_main_table(cells)
    print_drift_summary(cells)
    print_n3_anomaly(cells)
    print_beta_stability(cells; σf = 2.0)
    print_beta_stability(cells; σf = 0.25)
    print_cross_norm(cells; β = HEADLINE_β)
    print_n8_footnote(cells)

    # ── Build flat row tables for downstream consumers and save ──────────────
    raw_rows = NamedTuple[]
    for n in NS_MAIN, β in BETAS, σf in SIGMA_GRID
        r = get(cells, (n, Float64(β), Float64(σf)), nothing)
        r === nothing && continue
        push!(raw_rows, r)
    end
    # n=8 footnote rows (only β=0.25)
    n8_rows = NamedTuple[]
    for σf in (0.25, 0.5, 0.75, 1.0)
        r = get(cells, (N_FOOTNOTE, 0.25, Float64(σf)), nothing)
        r === nothing || push!(n8_rows, r)
    end

    ratio_rows = NamedTuple[]
    for n in NS_MAIN, β in BETAS, σf in SIGMA_GRID
        σf == BASELINE_σF && continue
        r = ratio_at(cells, n, β, σf)
        r === nothing || push!(ratio_rows, r)
    end
    n8_ratio_rows = NamedTuple[]
    for σf in (0.25, 0.5, 0.75)
        r = ratio_at(cells, N_FOOTNOTE, 0.25, σf)
        r === nothing || push!(n8_ratio_rows, r)
    end

    # ── Figures ──────────────────────────────────────────────────────────────
    plot_main(cells,
              out_path=joinpath(OUT_DIR, "qf_2dn_ckg_sigma_ratio_main"))
    plot_checks(cells,
                out_path=joinpath(OUT_DIR, "qf_2dn_ckg_sigma_ratio_checks"))

    bson_out = joinpath(OUT_DIR, "qf_2dn_ckg_sigma_ratio_data.bson")
    BSON.@save bson_out raw_rows n8_rows ratio_rows n8_ratio_rows
    println("\nData table saved: ", bson_out)

    # ── Verdict ──────────────────────────────────────────────────────────────
    # Aggregate drift of R_σ over n=4..7 across all (β, σ_f≠1) cells.
    drift_by = Dict{Tuple{Float64, Float64}, Float64}()
    for β in BETAS, σf in SIGMA_GRID
        σf == BASELINE_σF && continue
        Rs = Float64[]
        for n in 4:7
            r = ratio_at(cells, n, β, σf)
            r === nothing || push!(Rs, r.R_d11)
        end
        length(Rs) >= 2 || continue
        drift_by[(Float64(β), Float64(σf))] =
            (maximum(Rs) - minimum(Rs)) / mean(Rs)
    end
    cross_norm_max_diff = maximum(abs(r.R_d11 - r.R_hs) / abs(r.R_d11)
                                  for r in ratio_rows if r.n >= 4)
    drift_vals = collect(values(drift_by))

    println()
    println("==================== VERDICT ====================")
    @printf("R_σ drift (rel. range over n=4..7) per (β, σ_f):\n")
    βs = sort(unique(k[1] for k in keys(drift_by)))
    σfs = sort(unique(k[2] for k in keys(drift_by)))
    @printf("%-8s |", "β \\ σ_f")
    for σf in σfs; @printf(" %8.2f", σf); end; println()
    println(repeat("-", 10 + 9 * length(σfs)))
    for β in βs
        @printf("%-8.2f |", β)
        for σf in σfs
            v = get(drift_by, (β, σf), NaN)
            isnan(v) ? @printf("    —    ") : @printf(" %8.4f", v)
        end
        println()
    end
    @printf("\n    aggregate over %d (β, σ_f) cells:\n", length(drift_vals))
    @printf("        median drift = %.4f\n", median(drift_vals))
    @printf("        max    drift = %.4f\n", maximum(drift_vals))
    @printf("    cross-norm |R(d11) − R(HS)|/R(d11), max over n≥4 cells = %.5f\n",
            cross_norm_max_diff)
    println()
    println("READING:")
    println("  • Each drift entry is the n-range of R_σ(n; σ_f) across n=4..7")
    println("    at fixed (β, σ_f). FLAT (≲ 3 %) means d_{1→1} normalisation")
    println("    cancels the σ-induced rate change to within ε; DRIFTS means")
    println("    there is a STRUCTURAL piece of the σ-dependence that d11")
    println("    misses (kink-window width, broadening of γ(ω), etc.).")
    println("  • Cross-norm tracks whether HS-norm tells the same story.")
    println("  • n=3 is reported separately (small-Hilbert edge effect, c.f.")
    println("    qf-cwg). n=8 β=0.25 is a footnote (partial coverage).")
    println("=================================================\n")
    println("Done.")
end

main()
