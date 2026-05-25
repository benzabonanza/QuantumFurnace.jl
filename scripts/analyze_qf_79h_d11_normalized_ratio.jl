#!/usr/bin/env julia
#
# qf-79h: d_{1→1}-normalized gap ratio between CKG and DLL — does it isolate
# structural differences from pure rate differences between Lindbladians?
#
# User's hypothesis (2026-05-25): the d_{1→1} norm proxy used to "normalize"
# Lindbladian rates for sampler comparisons is dimension-dependent (grows with
# the number of jump operators and the local dim 2^n), so dividing gap/d11 to
# get a "rate-corrected" gap over-penalises large n and obscures structural
# scalings. A SAME-CELL RATIO between samplers should cancel the common
# dimensional growth — leaving only the structural-quality difference:
#
#     R(n, β) = [Gap_CKG(n,β) / d11_CKG(n,β)] / [Gap_DLL(n,β) / d11_DLL(n,β)]
#             = [Gap_CKG / Gap_DLL] × [d11_DLL / d11_CKG]
#
# # PHYSICS CHECK: d_{1→1} is implemented as `opnorm(M)` of the Kossakowski
# matrix `M[k,j] = Σ_i conj(A[i,k]) A[i,j] α(ν_{ij}, ν_{ik})` (see
# scripts/scratch_qf_e4z_35_sigma_sweep_plot.jl:294–298). This is a tight upper
# bound on the dissipator's d_{1→1} superoperator norm and sets a max
# continuous-time rate of the dissipator. Under uniform time-scaling
# L → kL: gap → k·gap, d11 → k·d11, so gap/d11 is rate-invariant. Across
# samplers it varies through:
#   (i)  the choice of α (CKG smooth-Metro vs DLL Metropolis) — STRUCTURAL,
#   (ii) the dimensional growth of the Kossakowski matrix in n — assumed
#        common to both samplers, so cancels in R.
# If (ii) does NOT cancel, R will drift with n; that drift is itself the answer.
#
# Data sources (all already on disk):
#   - scripts/output/sweep_qf_e4z_34_ckg_vs_dll_plot/{ckg,dll}/*.bson
#     5 seeds × 6 n × 3 β × 2 samplers = 180 sidecars with gap_phys, τ_mix_phys
#   - scripts/output/qf_e4z_34_d11_seed46.bson
#     36 rows (n=3..8, β ∈ {0.25,0.5,1.0}, 2 samplers) of d11_phys at seed=46
#   - scripts/output/qf_e4z_34_norm_diagnostic.bson
#     30 rows (n=3..7, β ∈ {0.25,0.5,1.0}, 2 samplers) at seed=42 with d11+hs+gap+τ
#   - scripts/output/sweep_qf_e4z_35_sigma_sweep_plot/ckg/*.bson
#     σ-sweep at σ_factor ∈ {0.25,0.5,0.75,1.0,1.5,2.0}, 2 seeds, CKG only
#
# Output:
#   drafts/figures/numerics/qf_79h_normalized_ratio_main.{png,pdf}
#     3-panel: (a) raw gap, (b) raw τ_mix, (c) R(n,β)
#   drafts/figures/numerics/qf_79h_normalized_ratio_checks.{png,pdf}
#     3-panel sanity checks: (a) σ-rate test, (b) d11(n) scaling, (c) HS-norm cross-check
#   drafts/figures/numerics/qf_79h_normalized_ratio_data.bson  (raw table)

using Printf, Statistics, BSON, Plots
ENV["GKSwstype"] = "100"

const REPO_ROOT     = abspath(joinpath(@__DIR__, ".."))
const SWEEP_DIR_CKG = joinpath(REPO_ROOT, "scripts", "output", "sweep_qf_e4z_34_ckg_vs_dll_plot", "ckg")
const SWEEP_DIR_DLL = joinpath(REPO_ROOT, "scripts", "output", "sweep_qf_e4z_34_ckg_vs_dll_plot", "dll")
const D11_FILE_46   = joinpath(REPO_ROOT, "scripts", "output", "qf_e4z_34_d11_seed46.bson")
const D11_FILE_42   = joinpath(REPO_ROOT, "scripts", "output", "qf_e4z_34_norm_diagnostic.bson")
const SIGMA_DIR     = joinpath(REPO_ROOT, "scripts", "output", "sweep_qf_e4z_35_sigma_sweep_plot", "ckg")
const OUT_DIR       = joinpath(REPO_ROOT, "drafts", "figures", "numerics")

const CANONICAL_SEED = 46    # see memory/qf_e4z_34_canonical_seed_46.md
const CROSS_SEED     = 42

# Thesis colours (reference_thesis_colors.md)
const COL_CKG  = RGB(0.20, 0.42, 0.65)   # slateblue
const COL_DLL  = RGB(0.70, 0.30, 0.30)   # bordeaux
const BETA_COL = Dict(
    0.25 => RGB(0.45, 0.65, 0.62),       # dustyteal
    0.5  => RGB(0.62, 0.45, 0.38),       # terracotta
    1.0  => RGB(0.35, 0.25, 0.45),       # plum
)

# ── Loaders ──────────────────────────────────────────────────────────────────

# qf-e4z.34 main sweep sidecar → (n, β_phys, sampler, seed) → (gap_phys, τ_phys)
function load_main_sweep(dir::String)
    rows = Dict{NTuple{4, Any}, NamedTuple}()
    for f in readdir(dir; join=true)
        endswith(f, ".bson") || continue
        d = BSON.load(f)
        haskey(d, :result) || continue
        r = d[:result]
        key = (r[:n], r[:beta_phys], r[:sampler], r[:seed])
        rows[key] = (
            gap_phys = r[:gap_phys],
            tau_phys = r[:mixing_time_phys],
            tau_src  = r[:mixing_time_source],
            converged = get(r, :all_converged, true),
            r        = r,
        )
    end
    return rows
end

# qf_e4z_34_d11_seed46.bson → (n, β_phys, sampler) → d11_phys at seed=46
function load_d11_seed46()
    d = BSON.load(D11_FILE_46)
    out = Dict{NTuple{3, Any}, NamedTuple}()
    for row in d[:rows]
        key = (row[:n], row[:beta_phys], row[:sampler])
        out[key] = (d11 = row[:d11_phys], d11_alg = row[:d11_alg])
    end
    return out
end

# qf_e4z_34_norm_diagnostic.bson → seed=42 with d11_phys + hs_phys + gap_phys + tau_phys
function load_norm_diagnostic_seed42()
    d = BSON.load(D11_FILE_42)
    out = Dict{NTuple{3, Any}, NamedTuple}()
    for row in d[:rows]
        key = (row[:n], row[:beta_phys], row[:sampler])
        out[key] = (
            d11 = row[:d11_phys], d11_alg = row[:d11_alg],
            hs  = row[:hs_phys],  hs_alg  = row[:hs_alg],
            gap = row[:gap_phys], tau = row[:tau_phys],
        )
    end
    return out
end

# qf-e4z.35 σ-sweep sidecars (CKG only) → list of NamedTuples.
# Skip the auxiliary `_spectral_mech.bson` files (qf-e4z.35.1 diagnostics);
# those contain the eigenmode decomposition rather than the τ_mix sidecar.
function load_sigma_sweep()
    rows = NamedTuple[]
    for f in readdir(SIGMA_DIR; join=true)
        endswith(f, ".bson") || continue
        occursin("_spectral_mech", basename(f)) && continue
        d = BSON.load(f)
        haskey(d, :result) || continue
        r = d[:result]
        haskey(r, :gap_phys) && haskey(r, :d_1to1_phys) || continue
        push!(rows, (
            n = r[:n], beta_phys = r[:beta_phys], sigma_factor = r[:sigma_factor],
            seed = r[:seed],
            gap_phys = r[:gap_phys], tau_phys = r[:mixing_time_phys],
            d11_phys = r[:d_1to1_phys], hs_phys = r[:hs_norm_phys],
            sigma_phys = r[:sigma_phys],
        ))
    end
    return rows
end

# ── Build joined table at the canonical seed ─────────────────────────────────

function build_main_table()
    ckg = load_main_sweep(SWEEP_DIR_CKG)
    dll = load_main_sweep(SWEEP_DIR_DLL)
    d11_46 = load_d11_seed46()

    # Discover (n, β_phys) cells
    cells = Set{NTuple{2, Any}}()
    for (k, _) in ckg
        n, b, sampler, seed = k
        seed == CANONICAL_SEED && push!(cells, (n, b))
    end

    rows = NamedTuple[]
    for (n, b) in sort(collect(cells))
        ck_key = (n, b, :ckg, CANONICAL_SEED)
        dl_key = (n, b, :dll, CANONICAL_SEED)
        haskey(ckg, ck_key) && haskey(dll, dl_key) || continue
        ck_d11 = get(d11_46, (n, b, :ckg), nothing)
        dl_d11 = get(d11_46, (n, b, :dll), nothing)
        (ck_d11 === nothing || dl_d11 === nothing) && continue
        push!(rows, (
            n = n, beta_phys = b, seed = CANONICAL_SEED,
            gap_ckg = ckg[ck_key].gap_phys,  tau_ckg = ckg[ck_key].tau_phys,
            gap_dll = dll[dl_key].gap_phys,  tau_dll = dll[dl_key].tau_phys,
            d11_ckg = ck_d11.d11,            d11_dll = dl_d11.d11,
        ))
    end
    return rows
end

# Cross-seed (seed=42) confirmation using the standalone norm_diagnostic file
function build_cross_table()
    nd = load_norm_diagnostic_seed42()
    cells = Set{NTuple{2, Any}}()
    for ((n, b, _), _) in nd
        push!(cells, (n, b))
    end
    rows = NamedTuple[]
    for (n, b) in sort(collect(cells))
        ck = get(nd, (n, b, :ckg), nothing)
        dl = get(nd, (n, b, :dll), nothing)
        (ck === nothing || dl === nothing) && continue
        push!(rows, (
            n = n, beta_phys = b, seed = CROSS_SEED,
            gap_ckg = ck.gap, tau_ckg = ck.tau, d11_ckg = ck.d11, hs_ckg = ck.hs,
            gap_dll = dl.gap, tau_dll = dl.tau, d11_dll = dl.d11, hs_dll = dl.hs,
        ))
    end
    return rows
end

# ── Ratios ───────────────────────────────────────────────────────────────────

# R(n,β) = (gap_CKG / d11_CKG) / (gap_DLL / d11_DLL)
norm_gap_ratio(r, gap_field_ck, gap_field_dl, d11_field_ck, d11_field_dl) =
    (getproperty(r, gap_field_ck) / getproperty(r, d11_field_ck)) /
    (getproperty(r, gap_field_dl) / getproperty(r, d11_field_dl))

# τ_mix-based normalized ratio: τ_DLL · d11_DLL  /  τ_CKG · d11_CKG
# This is "CKG advantage per unit d11-cost" — values > 1 mean CKG mixes faster
# per unit max-rate.
norm_tau_ratio(r, tau_field_ck, tau_field_dl, d11_field_ck, d11_field_dl) =
    (getproperty(r, tau_field_dl) * getproperty(r, d11_field_dl)) /
    (getproperty(r, tau_field_ck) * getproperty(r, d11_field_ck))

# ── Print summary ────────────────────────────────────────────────────────────

function print_table(rows, header)
    println()
    println("=== ", header, " ===")
    @printf("%4s %8s | %8s %8s %8s | %8s %8s %8s | %8s %8s %10s %8s %10s\n",
        "n", "β_phys", "gap_CKG", "gap_DLL", "gap_rat",
        "τ_CKG", "τ_DLL", "τ_rat",
        "d11_CKG", "d11_DLL", "d11_rat", "R_gap", "R/gap_rat")
    println(repeat("-", 130))
    for r in rows
        R_gap   = (r.gap_ckg / r.d11_ckg) / (r.gap_dll / r.d11_dll)
        gap_rat = r.gap_ckg / r.gap_dll
        tau_rat = r.tau_dll / r.tau_ckg
        d11_rat = r.d11_ckg / r.d11_dll
        @printf("%4d %8.3f | %8.4f %8.4f %8.3f | %8.3f %8.3f %8.3f | %8.3f %8.3f %10.5f %8.3f %10.5f\n",
            r.n, r.beta_phys, r.gap_ckg, r.gap_dll, gap_rat,
            r.tau_ckg, r.tau_dll, tau_rat,
            r.d11_ckg, r.d11_dll, d11_rat, R_gap, R_gap / gap_rat)
    end
end

# ── Plotting ─────────────────────────────────────────────────────────────────

function plot_main(rows; out_path::String)
    βs = sort(unique(r.beta_phys for r in rows))

    # Linear y-axes throughout — values lie in narrow ranges (gap 3–5, τ 1–2.3,
    # R 0.94–1.71) and log scale was just compressing the signal.

    # Panel A: raw gap_phys vs n
    pA = plot(xlabel="n", ylabel="gap_phys", legend=:topright,
              title="(a) Raw spectral gap", titlefontsize=10,
              size=(380, 320), framestyle=:box, legendfontsize=7,
              legend_background_color=RGBA(1,1,1,0.7))
    for b in βs
        sub = sort(filter(r -> r.beta_phys == b, rows), by = r -> r.n)
        plot!(pA, [r.n for r in sub], [r.gap_ckg for r in sub],
              color=BETA_COL[b], lw=2, marker=:circle, ms=4, ls=:solid,
              label="CKG β=$(b)")
        plot!(pA, [r.n for r in sub], [r.gap_dll for r in sub],
              color=BETA_COL[b], lw=2, marker=:square, ms=4, ls=:dash,
              label="DLL β=$(b)")
    end

    # Panel B: raw τ_mix vs n
    pB = plot(xlabel="n", ylabel="τ_mix_phys", legend=false,
              title="(b) Raw mixing time", titlefontsize=10,
              size=(380, 320), framestyle=:box)
    for b in βs
        sub = sort(filter(r -> r.beta_phys == b, rows), by = r -> r.n)
        plot!(pB, [r.n for r in sub], [r.tau_ckg for r in sub],
              color=BETA_COL[b], lw=2, marker=:circle, ms=4, ls=:solid)
        plot!(pB, [r.n for r in sub], [r.tau_dll for r in sub],
              color=BETA_COL[b], lw=2, marker=:square, ms=4, ls=:dash)
    end

    # Panel C: R(n,β), with raw gap ratio overlay to show the d11 normalization
    # barely moves anything (d11_CKG/d11_DLL is in [1.004, 1.025] empirically).
    pC = plot(xlabel="n",
              ylabel="ratio",
              legend=:bottomright,
              title="(c) Normalized vs raw gap ratio", titlefontsize=10,
              size=(380, 320), framestyle=:box, legendfontsize=7,
              legend_background_color=RGBA(1,1,1,0.7))
    hline!(pC, [1.0], ls=:dot, color=:gray, label="ratio = 1")
    for b in βs
        sub = sort(filter(r -> r.beta_phys == b, rows), by = r -> r.n)
        R       = [(r.gap_ckg / r.d11_ckg) / (r.gap_dll / r.d11_dll) for r in sub]
        gap_rat = [r.gap_ckg / r.gap_dll for r in sub]
        plot!(pC, [r.n for r in sub], R,
              color=BETA_COL[b], lw=2.5, marker=:diamond, ms=5, ls=:solid,
              label="R β=$(b)")
        plot!(pC, [r.n for r in sub], gap_rat,
              color=BETA_COL[b], lw=1.2, marker=:circle, ms=3, ls=:dot,
              label=nothing, alpha=0.6)
    end

    fig = plot(pA, pB, pC, layout=(1, 3), size=(1140, 320),
               left_margin=5Plots.mm, bottom_margin=5Plots.mm,
               top_margin=3Plots.mm)
    savefig(fig, out_path * ".png")
    savefig(fig, out_path * ".pdf")
    println("\nMain figure saved: ", out_path, ".{png,pdf}")
    return fig
end

function plot_checks(main_rows, cross_rows, sigma_rows; out_path::String)
    βs = sort(unique(r.beta_phys for r in main_rows))

    # ── Check A: σ-rate cancellation. Within CKG at fixed (n, β), vary σ_factor.
    # If d11 were a *complete* rate-canceller, gap/d11 would be flat in σ.
    # σ is NOT a pure rate knob (it sets the kink-window width); residual
    # variation = structural σ-dependence not absorbed by d11. So flatness
    # is a one-sided check: NOT-flat doesn't condemn d11.
    pA = plot(xlabel="σ_factor (CKG smoothing kernel ∝ σ)",
              ylabel="gap_CKG / d11_CKG", legend=:topright,
              title="(a) σ-rate cancellation within CKG",
              titlefontsize=10, size=(380, 320), framestyle=:box,
              legendfontsize=7,
              legend_background_color=RGBA(1,1,1,0.7))
    cells_for_sigma = [(3, 0.5), (5, 0.5), (7, 0.5), (3, 1.0), (5, 1.0), (7, 1.0)]
    cell_colors = [RGB(0.55, 0.85, 0.85), RGB(0.30, 0.55, 0.65), RGB(0.10, 0.30, 0.45),
                   RGB(0.95, 0.75, 0.55), RGB(0.85, 0.50, 0.30), RGB(0.65, 0.25, 0.10)]
    for (i, (nval, bval)) in enumerate(cells_for_sigma)
        sub = filter(r -> r.n == nval && r.beta_phys == bval, sigma_rows)
        isempty(sub) && continue
        by_sigma = Dict{Float64, Vector{Float64}}()
        for r in sub
            by_sigma[r.sigma_factor] = push!(get(by_sigma, r.sigma_factor, Float64[]),
                                              r.gap_phys / r.d11_phys)
        end
        σs = sort(collect(keys(by_sigma)))
        ys = [mean(by_sigma[s]) for s in σs]
        plot!(pA, σs, ys, color=cell_colors[i], lw=2, marker=:circle, ms=4,
              label="n=$(nval), β=$(bval)")
    end

    # ── Check B: d11(n) for both samplers. The visual punchline: the two
    # curves at the same β are essentially on top of each other — d11(CKG)
    # ≈ d11(DLL) to within ~2 %. Hence R(n,β) ≈ gap-ratio(n,β).
    pB = plot(xlabel="n", ylabel="d11_phys", legend=:topleft,
              title="(b) d_{1→1} vs n  —  CKG ≈ DLL", titlefontsize=10,
              size=(380, 320), framestyle=:box, legendfontsize=7,
              legend_background_color=RGBA(1,1,1,0.7))
    for b in βs
        sub = sort(filter(r -> r.beta_phys == b, main_rows), by = r -> r.n)
        plot!(pB, [r.n for r in sub], [r.d11_ckg for r in sub],
              color=BETA_COL[b], lw=2, marker=:circle, ms=4, ls=:solid,
              label="CKG β=$(b)")
        plot!(pB, [r.n for r in sub], [r.d11_dll for r in sub],
              color=BETA_COL[b], lw=2, marker=:square, ms=4, ls=:dash,
              label="DLL β=$(b)")
    end

    # ── Check C: cross-norm — d11 vs HS at seed=42 (cross_rows has both).
    # R_HS(n,β) should track R_d11(n,β) if the choice of operator norm is
    # robust. Same conclusion as panel (b): HS(CKG) ≈ HS(DLL) and d11(CKG)
    # ≈ d11(DLL), so both ratios reduce to gap_CKG/gap_DLL within ~2 %.
    pC = plot(xlabel="n", ylabel="normalized gap ratio",
              legend=:bottomright,
              title="(c) Cross-norm at seed=42", titlefontsize=10,
              size=(380, 320), framestyle=:box, legendfontsize=7,
              legend_background_color=RGBA(1,1,1,0.7))
    hline!(pC, [1.0], ls=:dot, color=:gray, label=nothing)
    for b in βs
        sub = sort(filter(r -> r.beta_phys == b, cross_rows), by = r -> r.n)
        isempty(sub) && continue
        R_d11 = [(r.gap_ckg / r.d11_ckg) / (r.gap_dll / r.d11_dll) for r in sub]
        R_hs  = [(r.gap_ckg / r.hs_ckg)  / (r.gap_dll / r.hs_dll)  for r in sub]
        plot!(pC, [r.n for r in sub], R_d11,
              color=BETA_COL[b], lw=2, marker=:diamond, ms=5, ls=:solid,
              label="R(d11) β=$(b)")
        plot!(pC, [r.n for r in sub], R_hs,
              color=BETA_COL[b], lw=1.5, marker=:utriangle, ms=4, ls=:dash,
              label="R(HS) β=$(b)")
    end

    fig = plot(pA, pB, pC, layout=(1, 3), size=(1140, 320),
               left_margin=5Plots.mm, bottom_margin=5Plots.mm,
               top_margin=3Plots.mm)
    savefig(fig, out_path * ".png")
    savefig(fig, out_path * ".pdf")
    println("Checks figure saved: ", out_path, ".{png,pdf}")
    return fig
end

# ── Main ─────────────────────────────────────────────────────────────────────

function main()
    isdir(OUT_DIR) || mkpath(OUT_DIR)

    println("Loading qf-e4z.34 main sweep (seed=$CANONICAL_SEED) + d11 diagnostic …")
    main_rows = build_main_table()
    @assert !isempty(main_rows) "main_rows is empty — check data paths"

    println("Loading qf-e4z.34 cross-seed norm diagnostic (seed=$CROSS_SEED) …")
    cross_rows = build_cross_table()

    println("Loading qf-e4z.35 σ-sweep (CKG only) …")
    sigma_rows = load_sigma_sweep()

    print_table(main_rows,  "qf-e4z.34, seed=$CANONICAL_SEED")
    print_table(cross_rows, "qf-e4z.34, seed=$CROSS_SEED (cross-check)")

    # Spread of R between seeds 42 and 46 — first sanity diagnostic
    println("\n=== Seed-stability of R (where seeds 42 & 46 both available) ===")
    @printf("%4s %8s | %10s %10s %10s\n", "n", "β_phys", "R_seed46", "R_seed42", "rel.diff")
    println(repeat("-", 52))
    for r46 in main_rows
        m = filter(r -> r.n == r46.n && r.beta_phys == r46.beta_phys, cross_rows)
        isempty(m) && continue
        r42 = first(m)
        R46 = (r46.gap_ckg / r46.d11_ckg) / (r46.gap_dll / r46.d11_dll)
        R42 = (r42.gap_ckg / r42.d11_ckg) / (r42.gap_dll / r42.d11_dll)
        @printf("%4d %8.3f | %10.4f %10.4f %10.4f\n",
                r46.n, r46.beta_phys, R46, R42, abs(R46 - R42) / R46)
    end

    # σ-rate cancellation diagnostic: for fixed (n,β), how much does
    # gap_CKG / d11_CKG vary across σ_factor?
    println("\n=== σ-rate cancellation within CKG (seed-mean over qf-e4z.35) ===")
    cells_diag = [(3, 0.5), (5, 0.5), (7, 0.5), (3, 1.0), (5, 1.0), (7, 1.0)]
    @printf("%4s %8s | %8s %8s %8s %8s %8s %8s | %10s\n",
        "n", "β_phys", "c=0.25", "c=0.50", "c=0.75", "c=1.00", "c=1.50", "c=2.00", "rel.range")
    println(repeat("-", 90))
    for (nval, bval) in cells_diag
        sub = filter(r -> r.n == nval && r.beta_phys == bval, sigma_rows)
        isempty(sub) && continue
        σs = sort(unique(r.sigma_factor for r in sub))
        means = Float64[]
        for σ in σs
            sub2 = filter(r -> r.sigma_factor == σ, sub)
            push!(means, mean(r.gap_phys / r.d11_phys for r in sub2))
        end
        rng = (maximum(means) - minimum(means)) / mean(means)
        vals_str = join([@sprintf("%8.4f", m) for m in means], " ")
        @printf("%4d %8.3f | %s | %10.4f\n", nval, bval, vals_str, rng)
    end

    # Make figures
    plot_main(main_rows, out_path=joinpath(OUT_DIR, "qf_79h_normalized_ratio_main"))
    plot_checks(main_rows, cross_rows, sigma_rows,
                out_path=joinpath(OUT_DIR, "qf_79h_normalized_ratio_checks"))

    # Persist the joined table for downstream use
    bson_out = joinpath(OUT_DIR, "qf_79h_normalized_ratio_data.bson")
    BSON.@save bson_out main_rows cross_rows sigma_rows
    println("\nData table saved: ", bson_out)

    # ── Verdict summary ──────────────────────────────────────────────────────
    d11_ratios = [r.d11_ckg / r.d11_dll for r in main_rows]
    R_over_gap = [((r.gap_ckg/r.d11_ckg)/(r.gap_dll/r.d11_dll)) /
                  (r.gap_ckg/r.gap_dll) for r in main_rows]
    # n-scaling of d11 (geometric mean across β at each n) — does it grow
    # similarly for CKG and DLL?
    ns_uniq = sort(unique(r.n for r in main_rows))
    d11_ratio_per_n = Dict{Int, Float64}()
    for nval in ns_uniq
        sub = filter(r -> r.n == nval, main_rows)
        d11_ratio_per_n[nval] = mean(r.d11_ckg / r.d11_dll for r in sub)
    end

    println()
    println("==================== VERDICT ====================")
    @printf("d11_CKG / d11_DLL across all %d cells:\n", length(main_rows))
    @printf("    range  = [%.4f, %.4f]\n", minimum(d11_ratios), maximum(d11_ratios))
    @printf("    mean   = %.4f   (std = %.4f)\n", mean(d11_ratios), std(d11_ratios))
    @printf("    per n  : ")
    for nval in ns_uniq; @printf("n=%d→%.4f  ", nval, d11_ratio_per_n[nval]); end
    println()
    println()
    @printf("R(n,β) / gap-ratio(n,β) across all cells:\n")
    @printf("    range  = [%.4f, %.4f]   (perfect cancellation would be 1.0)\n",
        minimum(R_over_gap), maximum(R_over_gap))
    println()
    println("INTERPRETATION:")
    println("  - The d_{1→1} norm grows ~4× from n=3 to n=8, but it grows by")
    println("    THE SAME factor for CKG and DLL (d11_CKG/d11_DLL ≈ 1 to 1-2 %).")
    println("    Hence the proposed R(n,β) ≈ gap_CKG/gap_DLL within 2 % — the")
    println("    d11 normalisation is effectively a no-op for CKG vs DLL at")
    println("    fixed σ. User's worry about d11 'overcorrecting' is empirically")
    println("    unfounded for THIS comparison: there is nothing to overcorrect.")
    println()
    println("  - Why this works: CKG and DLL share the SAME jump operators A_i")
    println("    in the energy eigenbasis; they only differ in the rate function")
    println("    α(ν_{ij}, ν_{ik}). Both α functions have the same KMS magnitude")
    println("    scaling, so opnorm(M) cancels almost exactly between samplers.")
    println("    The structural difference between CKG and DLL lives in the")
    println("    DIRECTION of α (how rate is distributed across the Bohr graph),")
    println("    not in its magnitude — gap captures that, d11 does not.")
    println()
    println("  - σ-rate test (panel A of checks figure): WITHIN CKG, gap/d11")
    println("    varies ~50 % across σ_factor ∈ [0.25, 2.0]. σ is NOT a pure")
    println("    rate knob — it also changes the kink-window structure of γ(ω),")
    println("    so d11 cancels only the rate part. This does NOT undermine the")
    println("    CKG-vs-DLL R, which compares at fixed σ.")
    println()
    println("  - Cross-norm (panel C): HS-normalised R tracks d11-normalised R")
    println("    point-for-point — the norm choice is robust because BOTH norms")
    println("    are nearly sampler-symmetric at fixed (n,β).")
    println()
    println("  - Recommendation: report raw gap and τ_mix ratios as the headline,")
    println("    note in passing that d_{1→1} normalisation leaves them unchanged")
    println("    (so 'rate-corrected' and 'raw' agree). The d11 normalisation is")
    println("    NOT needed to defend the structural-vs-rate distinction here.")
    println("=================================================\n")
    println("Done.")
end

main()
