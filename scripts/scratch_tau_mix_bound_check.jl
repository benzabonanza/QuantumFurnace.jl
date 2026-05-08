#!/usr/bin/env julia
#
# τ_mix bound check — predicted (gap-based) vs measured (bi-exp + extrapolation).
# Beads issue: qf-mto.6  (epic qf-mto: KMS Dirichlet fair comparison).
#
# Goal:
#   Validate the spectral-gap prediction
#       τ_pred(ε) = (1/λ) · log( 1 / (p_min · ε) )                    (1)
#   against the empirical τ_mix from the existing biexp-fit pipeline. If (1) is
#   tight (slope ≈ 1, small intercept on log–log axes), the FAIR.4 ρ_intrinsic
#   comparison built on λ alone is meaningful. A loose bound (large positive
#   intercept) means the comparison still works *up to a constant factor*.
#
# Inputs:
#   1) scripts/output/fair_comparison_kms_dirichlet.bson   (FAIR.4 — has λ,
#      tau_mix_predicted, tau_mix_measured (often), p_min per (n, β, sampler)).
#   2) drafts/figures/numerics/ckg_vs_dll_taumix.bson      (qf-lkb.6 figure
#      data — measured τ_mix at β ∈ {1, 2, 5, 10}, n ∈ {3, 4, 5}).
#   3) scripts/output/taumix_beta20/taumix_beta20.bson     (β extension —
#      measured τ_mix at β ∈ {10, 15, 20}, n ∈ {3, 4, 5}).
#
# FAIR.4 grid:  n ∈ {3, 4, 5},  β ∈ {1, 5, 10, 20},  3 samplers ⇒ 36 cells.
# Measured grid (union of #2 and #3): β ∈ {1, 2, 5, 10, 15, 20}.
# Cross-check intersection: β ∈ {1, 5, 10, 20}, n ∈ {3, 4, 5}, 3 samplers ⇒
#   ≤ 36 cells; some have NaN measured τ_mix (bi-exp fit failed) — skipped.
#
# Output:
#   scripts/output/tau_mix_bound_check.bson
#   drafts/figures/numerics/tau_mix_predicted_vs_measured.{png, pdf}
#
# PHYSICS CHECK: ε = target_epsilon = 1e-3 is the τ_mix target threshold used
# throughout the qf-lkb pipeline. Verified by reading FAIR.4's :epsilon field.
#
# PHYSICS CHECK: p_min = min eigenvalue of the Gibbs state σ(β,H). At low T
# (β = 20), several p_min values are ~1e-4 or numerically zero (n=5, β=20),
# which sends τ_pred → ∞ via log(1/(p_min·ε)). Such cells are dropped from
# the regression because the predicted side becomes degenerate.
#
# Usage: julia --project scripts/scratch_tau_mix_bound_check.jl

using Printf
using Statistics
using LinearAlgebra
using BSON
ENV["GKSwstype"] = "100"   # headless GR backend (matches the qf-lkb.6 plot script)
using Plots

# ── Paths ─────────────────────────────────────────────────────────────────────
const REPO_ROOT     = normpath(joinpath(@__DIR__, ".."))
const FAIR4_BSON    = joinpath(REPO_ROOT, "scripts", "output",
                               "fair_comparison_kms_dirichlet.bson")
const TAU_BSON_LO   = joinpath(REPO_ROOT, "drafts", "figures", "numerics",
                               "ckg_vs_dll_taumix.bson")
const TAU_BSON_HI   = joinpath(REPO_ROOT, "scripts", "output", "taumix_beta20",
                               "taumix_beta20.bson")
const OUT_BSON      = joinpath(REPO_ROOT, "scripts", "output",
                               "tau_mix_bound_check.bson")
const FIG_PNG       = joinpath(REPO_ROOT, "drafts", "figures", "numerics",
                               "tau_mix_predicted_vs_measured.png")
const FIG_PDF       = joinpath(REPO_ROOT, "drafts", "figures", "numerics",
                               "tau_mix_predicted_vs_measured.pdf")

mkpath(dirname(OUT_BSON))
mkpath(dirname(FIG_PNG))

# ── Sampler key mapping ──────────────────────────────────────────────────────
# FAIR.4 uses Symbol keys (:ckg_smooth_metro, :dll_metro, :dll_gauss).
# Measured BSONs identify the sampler by (construction::String, filter_name::String):
#   ckg_smooth_metro  ↔ ("KMS", "default")
#   dll_metro         ↔ ("DLL", "DLLMetropolis")
#   dll_gauss         ↔ ("DLL", "DLLGaussian")
const SAMPLER_KEYS = (:ckg_smooth_metro, :dll_metro, :dll_gauss)

const SAMPLER_LABEL = Dict(
    :ckg_smooth_metro => "CKG smooth-Metro",
    :dll_metro        => "DLL Metropolis",
    :dll_gauss        => "DLL Gaussian",
)

# Thesis colour palette — same family as the qf-lkb.6 comparison plot.
# Pinegreen (CKG), terracotta (DLL Metro), bordeaux (DLL Gauss) read as a
# distinct, perceptually well-separated triple on log–log axes.
const SAMPLER_COLOR = Dict(
    :ckg_smooth_metro => "#2D5A3D",   # pinegreen
    :dll_metro        => "#B5654A",   # terracotta
    :dll_gauss        => "#7A2E39",   # bordeaux
)

const SAMPLER_MARKER = Dict(
    :ckg_smooth_metro => :circle,
    :dll_metro        => :diamond,
    :dll_gauss        => :utriangle,
)

function classify_measured(r::NamedTuple)::Union{Symbol, Nothing}
    # `r` is a row from one of the measured BSONs; (construction, filter_name)
    # determine the sampler key. Return `nothing` if the combo isn't recognised.
    if r.construction == "KMS" && r.filter_name == "default"
        return :ckg_smooth_metro
    elseif r.construction == "DLL" && r.filter_name == "DLLMetropolis"
        return :dll_metro
    elseif r.construction == "DLL" && r.filter_name == "DLLGaussian"
        return :dll_gauss
    else
        return nothing
    end
end

# ── Build measured-τ_mix lookup ───────────────────────────────────────────────
# Key: (sampler_key::Symbol, n::Int, β::Float64) → mixing_time::Float64
# Priority: FAIR.4 :tau_mix_measured (since it's already populated for many cells
# and is consistent with the other FAIR.4 fields), then beta20 BSON, then the
# legacy qf-lkb.6 BSON. NaN values are treated as "missing" and let lower-priority
# sources fill in.
const MeasuredKey = Tuple{Symbol, Int, Float64}

"""
    build_measured_lookup(fair4, tau_lo, tau_hi)

Return `(measured::Dict, counts::NamedTuple)` where `counts` reports how many
entries each source contributed (only counts new keys — the priority order
inside the function controls precedence).
"""
function build_measured_lookup(fair4, tau_lo, tau_hi)
    measured = Dict{MeasuredKey, Float64}()

    # Priority 1: FAIR.4's own tau_mix_measured field.
    fair4_count = 0
    for r in fair4[:results]
        skey = r.sampler_key
        n    = Int(r.n)
        β    = float(r.β)
        if isfinite(r.tau_mix_measured) && !haskey(measured, (skey, n, β))
            measured[(skey, n, β)] = r.tau_mix_measured
            fair4_count += 1
        end
    end

    # Priority 2: taumix_beta20 BSON (β ∈ {10, 15, 20}).
    hi_count = 0
    for key in (:results_ckg, :results_dll_metro, :results_dll_gauss)
        for r in tau_hi[key]
            skey = classify_measured(r)
            skey === nothing && continue
            n = Int(r.n)
            β = float(r.beta)
            if isfinite(r.mixing_time) && !haskey(measured, (skey, n, β))
                measured[(skey, n, β)] = r.mixing_time
                hi_count += 1
            end
        end
    end

    # Priority 3: legacy qf-lkb.6 BSON (β ∈ {1, 2, 5, 10}).
    lo_count = 0
    for key in (:results_ckg, :results_dll_metro, :results_dll_gauss)
        for r in tau_lo[key]
            skey = classify_measured(r)
            skey === nothing && continue
            n = Int(r.n)
            β = float(r.beta)
            if isfinite(r.mixing_time) && !haskey(measured, (skey, n, β))
                measured[(skey, n, β)] = r.mixing_time
                lo_count += 1
            end
        end
    end

    return measured, (fair4=fair4_count, beta20=hi_count, lkb6=lo_count)
end

struct CellRow
    sampler::Symbol
    n::Int
    β::Float64
    λ::Float64
    p_min::Float64
    τ_pred::Float64
    τ_meas::Float64
    ratio::Float64    # τ_pred / τ_meas
end

"""
    build_rows(fair4, measured)

Walk the FAIR.4 rows; pair the predicted side with `measured` lookup. Return
`(rows::Vector{CellRow}, skipped::Vector{NamedTuple})`.
"""
function build_rows(fair4, measured::Dict{MeasuredKey, Float64})
    rows = CellRow[]
    skipped = NamedTuple{(:sampler, :n, :β, :reason),
                         Tuple{Symbol, Int, Float64, String}}[]
    for r in fair4[:results]
        skey = r.sampler_key
        n    = Int(r.n)
        β    = float(r.β)
        τp   = r.tau_mix_predicted
        pmin = r.p_min
        λ    = r.λ
        τm   = get(measured, (skey, n, β), NaN)

        if !isfinite(τp)
            push!(skipped, (sampler=skey, n=n, β=β, reason="τ_pred non-finite"))
            continue
        end
        if !isfinite(τm)
            push!(skipped, (sampler=skey, n=n, β=β, reason="τ_meas missing"))
            continue
        end
        if τm <= 0 || τp <= 0
            push!(skipped, (sampler=skey, n=n, β=β, reason="non-positive τ"))
            continue
        end
        push!(rows, CellRow(skey, n, β, λ, pmin, τp, τm, τp/τm))
    end
    return rows, skipped
end

# ── Statistics ────────────────────────────────────────────────────────────────
function log_ratio_stats(rows, skey)
    sub = filter(r -> r.sampler == skey, rows)
    isempty(sub) && return (mean=NaN, std=NaN, n=0)
    lr = log10.([r.ratio for r in sub])
    return (mean=Statistics.mean(lr),
            std=length(lr) > 1 ? Statistics.std(lr) : 0.0,
            n=length(lr))
end

# Linear regression log10(τ_pred) = α · log10(τ_meas) + β  per sampler.
# Returns (slope, intercept, R²). Closed-form OLS.
function regress_logpred_vs_logmeas(rows, skey::Union{Symbol, Nothing}=nothing)
    sub = isnothing(skey) ? rows : filter(r -> r.sampler == skey, rows)
    length(sub) < 2 && return (slope=NaN, intercept=NaN, r2=NaN, n=length(sub))
    x = log10.([r.τ_meas for r in sub])
    y = log10.([r.τ_pred for r in sub])
    n = length(x)
    x̄, ȳ = mean(x), mean(y)
    Sxx = sum((x .- x̄).^2)
    Sxy = sum((x .- x̄) .* (y .- ȳ))
    slope = Sxy / Sxx
    intercept = ȳ - slope * x̄
    yhat = slope .* x .+ intercept
    ss_res = sum((y .- yhat).^2)
    ss_tot = sum((y .- ȳ).^2)
    r2 = ss_tot > 0 ? 1.0 - ss_res / ss_tot : 1.0
    return (slope=slope, intercept=intercept, r2=r2, n=n)
end

function print_sampler_table(rows, skey)
    sub = filter(r -> r.sampler == skey, rows)
    isempty(sub) && (println("(no rows for $skey)"); return)
    sort!(sub, by = r -> (r.n, r.β))
    println("\n", "="^88)
    println("Sampler: ", SAMPLER_LABEL[skey], "  (key=:", skey, ")")
    println("="^88)
    @printf("%-3s %-6s %-10s %-10s %-12s %-12s %-10s\n",
            "n", "β", "λ", "p_min", "τ_pred", "τ_meas", "ratio")
    println("-"^88)
    for r in sub
        @printf("%-3d %-6.1f %-10.4f %-10.3e %-12.4e %-12.4e %-10.3f\n",
                r.n, r.β, r.λ, r.p_min, r.τ_pred, r.τ_meas, r.ratio)
    end
end

function find_anomalies(rows; z_threshold=2.0)
    anomalies = NamedTuple{(:sampler, :n, :β, :ratio, :z),
                           Tuple{Symbol, Int, Float64, Float64, Float64}}[]
    for skey in SAMPLER_KEYS
        sub = filter(r -> r.sampler == skey, rows)
        isempty(sub) && continue
        lr = log10.([r.ratio for r in sub])
        μ, σ = mean(lr), length(lr) > 1 ? Statistics.std(lr) : 0.0
        σ <= 0 && continue
        for r in sub
            dev = log10(r.ratio) - μ
            if abs(dev) > z_threshold * σ
                push!(anomalies, (sampler=skey, n=r.n, β=r.β,
                                  ratio=r.ratio, z=dev/σ))
            end
        end
    end
    return anomalies
end

# ── Main ──────────────────────────────────────────────────────────────────────
function main()
    println("="^78)
    println("τ_mix bound check (qf-mto.6)")
    println("="^78)
    @printf("FAIR.4 BSON           : %s\n", FAIR4_BSON)
    @printf("Measured BSON (β≤10)  : %s\n", TAU_BSON_LO)
    @printf("Measured BSON (β=10..20): %s\n", TAU_BSON_HI)
    println("="^78)
    flush(stdout)

    isfile(FAIR4_BSON)  || error("Missing FAIR.4 BSON: $FAIR4_BSON")
    isfile(TAU_BSON_LO) || error("Missing measured BSON: $TAU_BSON_LO")
    isfile(TAU_BSON_HI) || error("Missing measured BSON: $TAU_BSON_HI")

    fair4  = BSON.load(FAIR4_BSON)
    tau_lo = BSON.load(TAU_BSON_LO)
    tau_hi = BSON.load(TAU_BSON_HI)

    epsilon = fair4[:epsilon]
    @printf("ε (target_epsilon, FAIR.4) : %.1e\n", epsilon)
    # Sanity: the measured BSONs must use the same ε.
    let ε_lo = tau_lo[:target_epsilon], ε_hi = tau_hi[:target_epsilon]
        if !(isapprox(ε_lo, epsilon) && isapprox(ε_hi, epsilon))
            @warn "ε mismatch across BSONs" epsilon ε_lo ε_hi
        end
    end

    measured, counts = build_measured_lookup(fair4, tau_lo, tau_hi)
    @printf("Measured-τ_mix entries: %d (fair4: +%d, beta20: +%d, lkb6: +%d)\n",
            length(measured), counts.fair4, counts.beta20, counts.lkb6)
    flush(stdout)

    rows, skipped = build_rows(fair4, measured)
    @printf("Comparison cells included: %d\n", length(rows))
    @printf("Comparison cells skipped : %d\n", length(skipped))
    if !isempty(skipped)
        println("Skip details:")
        for s in skipped
            @printf("  sampler=%-16s n=%d β=%-5.1f reason=%s\n",
                    string(s.sampler), s.n, s.β, s.reason)
        end
    end
    flush(stdout)

    # Per-sampler tables
    for skey in SAMPLER_KEYS
        print_sampler_table(rows, skey)
    end

    # Summary stats
    println("\n", "="^88)
    println("Summary statistics: log10(τ_pred / τ_meas) and OLS fit")
    println("="^88)
    @printf("%-22s %5s %12s %12s %10s %12s %8s\n",
            "sampler", "n", "mean log10R", "std log10R", "slope", "intercept", "R²")
    println("-"^88)
    for skey in SAMPLER_KEYS
        stats = log_ratio_stats(rows, skey)
        fit   = regress_logpred_vs_logmeas(rows, skey)
        @printf("%-22s %5d %12.4f %12.4f %10.4f %12.4f %8.4f\n",
                SAMPLER_LABEL[skey], stats.n, stats.mean, stats.std,
                fit.slope, fit.intercept, fit.r2)
    end
    pooled = regress_logpred_vs_logmeas(rows, nothing)
    pooled_lr = log10.([r.ratio for r in rows])
    pooled_mean = isempty(pooled_lr) ? NaN : mean(pooled_lr)
    pooled_std  = length(pooled_lr) > 1 ? Statistics.std(pooled_lr) : 0.0
    @printf("%-22s %5d %12.4f %12.4f %10.4f %12.4f %8.4f\n",
            "POOLED (all)", pooled.n, pooled_mean, pooled_std,
            pooled.slope, pooled.intercept, pooled.r2)
    flush(stdout)

    # Acceptance checks
    println("\n", "="^88)
    println("Acceptance checks (per sampler)")
    println("="^88)
    for skey in SAMPLER_KEYS
        fit = regress_logpred_vs_logmeas(rows, skey)
        slope_ok     = isfinite(fit.slope) && 0.8 ≤ fit.slope ≤ 1.2
        intercept_ok = isfinite(fit.intercept) && abs(fit.intercept) < 1.0
        @printf("  %-22s  slope=%6.3f (in [0.8, 1.2]? %s)  intercept=%6.3f (|·|<1? %s)\n",
                SAMPLER_LABEL[skey], fit.slope, slope_ok ? "yes" : "NO ",
                fit.intercept, intercept_ok ? "yes" : "NO ")
    end

    # Headline ratio
    if !isempty(rows)
        lr = log10.([r.ratio for r in rows])
        headline_factor = 10.0 ^ mean(lr)
        headline_low    = 10.0 ^ (mean(lr) - Statistics.std(lr))
        headline_high   = 10.0 ^ (mean(lr) + Statistics.std(lr))
        println("\n", "="^88)
        @printf("HEADLINE: τ_pred = (1/λ) log(1/(p_min·ε)) overestimates measured τ_mix\n")
        @printf("          by a factor of ~%.2fx (1σ band: %.2fx-%.2fx) across %d cells.\n",
                headline_factor, headline_low, headline_high, length(rows))
        println("="^88)
    end

    # Anomalies
    println("\n", "="^88)
    println("Anomalous cells (|log10 ratio - sampler mean| > 2σ)")
    println("="^88)
    anomalies = find_anomalies(rows; z_threshold=2.0)
    if isempty(anomalies)
        println("  none — all per-sampler residuals within 2σ.")
    else
        for a in anomalies
            @printf("  %-22s n=%d β=%-5.1f ratio=%6.3f  z=%+.2f\n",
                    SAMPLER_LABEL[a.sampler], a.n, a.β, a.ratio, a.z)
        end
    end
    flush(stdout)

    # ── Save BSON ────────────────────────────────────────────────────────────
    rows_any = Vector{Any}(undef, length(rows))
    for (i, r) in enumerate(rows)
        rows_any[i] = (sampler=r.sampler, n=r.n, β=r.β, λ=r.λ, p_min=r.p_min,
                       τ_pred=r.τ_pred, τ_meas=r.τ_meas, ratio=r.ratio)
    end
    fits = Dict{Symbol, Any}()
    for skey in SAMPLER_KEYS
        f = regress_logpred_vs_logmeas(rows, skey)
        fits[skey] = Dict(:slope => f.slope, :intercept => f.intercept,
                          :r2 => f.r2, :n => f.n)
    end
    fits[:pooled] = Dict(:slope => pooled.slope, :intercept => pooled.intercept,
                         :r2 => pooled.r2, :n => pooled.n)
    BSON.bson(OUT_BSON, Dict(
        :rows           => rows_any,
        :skipped        => Vector{Any}(skipped),
        :fits           => fits,
        :epsilon        => epsilon,
        :sources        => Dict(
            :fair4   => FAIR4_BSON,
            :beta20  => TAU_BSON_HI,
            :legacy  => TAU_BSON_LO,
        ),
        :samplers       => collect(SAMPLER_KEYS),
        :sampler_labels => Dict(skey => SAMPLER_LABEL[skey] for skey in SAMPLER_KEYS),
    ))
    @info "Saved comparison BSON" OUT_BSON

    # ── Plot ─────────────────────────────────────────────────────────────────
    if isempty(rows)
        println("\nNo rows to plot.")
        return
    end

    τ_meas_all = [r.τ_meas for r in rows]
    τ_pred_all = [r.τ_pred for r in rows]
    xmin = minimum(τ_meas_all) * 0.7
    xmax = maximum(τ_meas_all) * 1.4
    ymin = min(minimum(τ_pred_all), xmin) * 0.7
    ymax = max(maximum(τ_pred_all), xmax) * 1.4
    ref_lo = min(xmin, ymin)
    ref_hi = max(xmax, ymax)

    plt = Plots.plot(
        xlabel = raw"measured $\tau_{\mathrm{mix}}$",
        ylabel = raw"predicted $\tau_{\mathrm{mix}} = \lambda^{-1}\,\log\frac{1}{p_{\min}\varepsilon}$",
        xscale = :log10, yscale = :log10,
        xlim = (xmin, xmax), ylim = (ymin, ymax),
        legend = :topleft,
        size   = (820, 720),
        title  = "Spectral-gap bound vs measured τ_mix (qf-mto.6)",
        titlefontsize = 11,
        guidefontsize = 11,
        legendfontsize = 9,
        dpi = 200,
    )

    Plots.plot!(plt, [ref_lo, ref_hi], [ref_lo, ref_hi];
        color = :black, linestyle = :dot, linewidth = 1.2,
        label = "y = x (tight bound)")

    for skey in SAMPLER_KEYS
        sub = filter(r -> r.sampler == skey, rows)
        isempty(sub) && continue
        xs = [r.τ_meas for r in sub]
        ys = [r.τ_pred for r in sub]
        Plots.scatter!(plt, xs, ys;
            color = SAMPLER_COLOR[skey],
            marker = SAMPLER_MARKER[skey],
            markersize = 7, markerstrokewidth = 0.3,
            markeralpha = 0.85,
            label = SAMPLER_LABEL[skey])
        fit = regress_logpred_vs_logmeas(rows, skey)
        if isfinite(fit.slope)
            xfit_lo, xfit_hi = extrema(xs)
            xgrid = exp10.(range(log10(xfit_lo*0.9), log10(xfit_hi*1.1), length=64))
            ygrid = exp10.(fit.slope .* log10.(xgrid) .+ fit.intercept)
            Plots.plot!(plt, xgrid, ygrid;
                color = SAMPLER_COLOR[skey],
                linestyle = :dash, linewidth = 1.6,
                label = @sprintf("%s fit: slope=%.2f, intercept=%.2f",
                                 SAMPLER_LABEL[skey], fit.slope, fit.intercept))
        end
    end

    Plots.savefig(plt, FIG_PNG)
    Plots.savefig(plt, FIG_PDF)
    @info "Saved figure" FIG_PNG FIG_PDF

    println("\nDone.")
    return nothing
end

main()
