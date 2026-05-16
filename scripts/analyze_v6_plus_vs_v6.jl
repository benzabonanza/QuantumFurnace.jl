#!/usr/bin/env julia
# analyze_v6_plus_vs_v6.jl  (qf-e4z.31)
#
# Compares the new rho_0 = |+⟩⟨+|^⊗N + krylovdim=60 multiseed sweep (v6_plus)
# to the original rho_0 = I/d + krylovdim=40 multiseed sweep (v6, qf-e4z.23).
#
# The goal is to settle the even/odd-n splitting interpretation:
#   • If v6_plus gap_phys per (n, β_phys) median MATCHES v6 → splitting is
#     real physics (case a).
#   • If v6_plus gap_phys per (n, β_phys) median is SYSTEMATICALLY LARGER
#     than v6 on even-n (and equal on odd-n) → v6 was reporting the
#     parity-EVEN sub-spectrum eigenvalue, the true gap doesn't collapse,
#     and the even/odd-n splitting is an I/d-convention artefact (case b).
#
# Also reports:
#   • Pass1↔Pass2 agreement per cell (sanity check that |+⟩⟨+|^⊗N is in
#     fact parity-broken: the two passes should agree to 1e-9 or better).
#   • τ_mix change v6 → v6_plus (expected ≈ 0 since trajectory accuracy is
#     symmetry-invariant per qf-e4z.27 — c_i for P̂-odd modes are zero for
#     a P̂-even rho_0, so the v6 trajectory was correct even with the
#     parity-trapped Krylov decomp).
#
# Usage:
#   julia --project scripts/analyze_v6_plus_vs_v6.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using QuantumFurnace
using LinearAlgebra
using BSON
using Printf
using Dates
using Statistics

println("[init] $(now())")

const V6_DIR      = joinpath(@__DIR__, "output",
                             "sweep_S1_v6_ckg_ideal_multiseed", "smooth_metro_eps1e-03")
const V6_PLUS_DIR = joinpath(@__DIR__, "output",
                             "sweep_S1_v6_plus_ckg_ideal_multiseed", "smooth_metro_eps1e-03")
const FIG_DIR = joinpath(@__DIR__, "..", "drafts", "figures", "numerics")
mkpath(FIG_DIR)

# --- Loaders -------------------------------------------------------------

function load_sidecars(dir)
    rows = NamedTuple[]
    isdir(dir) || return rows
    for fname in sort(readdir(dir))
        endswith(fname, ".bson") || continue
        path = joinpath(dir, fname)
        try
            d = BSON.load(path, QuantumFurnace)
            r = d[:result]
            push!(rows, (; (Symbol(k) => v for (k, v) in pairs(r))...))
        catch err
            @warn "Failed to load sidecar" path err
        end
    end
    sort!(rows; by = r -> (r.n, r.beta_phys, r.seed))
    return rows
end

v6_rows      = load_sidecars(V6_DIR)
v6_plus_rows = load_sidecars(V6_PLUS_DIR)
@printf("[load] v6      : %d sidecars\n", length(v6_rows))
@printf("[load] v6_plus : %d sidecars\n", length(v6_plus_rows))

gap_phys(r) = r.gap_arnoldi * r.rescaling_factor

function group_by_cell(rows)
    cells = Dict{Tuple{Int, Float64}, Vector{NamedTuple}}()
    for r in rows
        push!(get!(cells, (r.n, r.beta_phys), NamedTuple[]), r)
    end
    return cells
end

v6_cells      = group_by_cell(v6_rows)
v6_plus_cells = group_by_cell(v6_plus_rows)

# Only compare cells present in BOTH sweeps (v6_plus is currently n=3..7
# while v6 has n=3..8; this excludes the n=8 row which is qf-e4z.32 scope).
common_keys = sort!(collect(intersect(keys(v6_cells), keys(v6_plus_cells))))
@printf("[load] common cells: %d  (n ∈ %s, β_phys ∈ %s)\n",
        length(common_keys),
        sort!(unique([k[1] for k in common_keys])),
        sort!(unique([k[2] for k in common_keys])))

# --- Per-cell comparison table --------------------------------------------

println("\n" * "="^130)
println("Per-cell comparison: median gap_phys v6 vs v6_plus, ratio = v6_plus / v6")
println("="^130)
@printf("%-3s %-7s %-3s | %-12s %-12s %-12s | %-12s %-12s %-7s\n",
        "n", "β_phys", "#s", "v6 med gap", "v6+ med gap", "v6+/v6",
        "v6 med τ", "v6+ med τ", "τ ratio")
println("-"^130)

ratios = NamedTuple[]
for k in common_keys
    n, β = k
    c_v6  = v6_cells[k]
    c_v6p = v6_plus_cells[k]
    nseeds = min(length(c_v6), length(c_v6p))

    gap_v6_med  = median(gap_phys.(c_v6))
    gap_v6p_med = median(gap_phys.(c_v6p))
    gap_ratio   = gap_v6p_med / max(gap_v6_med, 1e-30)

    tau_v6_med  = median([r.mixing_time for r in c_v6])
    tau_v6p_med = median([r.mixing_time for r in c_v6p])
    tau_ratio   = tau_v6p_med / max(tau_v6_med, 1e-30)

    pass12_diffs = [r.gap_rel_diff_pass12 for r in c_v6p if haskey(r, :gap_rel_diff_pass12)]
    pass12_max   = isempty(pass12_diffs) ? NaN : maximum(pass12_diffs)

    push!(ratios, (n=n, beta_phys=β, gap_ratio=gap_ratio, tau_ratio=tau_ratio,
                   gap_v6=gap_v6_med, gap_v6p=gap_v6p_med,
                   tau_v6=tau_v6_med, tau_v6p=tau_v6p_med,
                   pass12_max=pass12_max, n_seeds=nseeds))

    @printf("%-3d %-7.2f %-3d | %-12.4g %-12.4g %-12.3f | %-12.4g %-12.4g %-7.3f\n",
            n, β, nseeds,
            gap_v6_med, gap_v6p_med, gap_ratio,
            tau_v6_med, tau_v6p_med, tau_ratio)
end

# --- Pass1 ↔ Pass2 sanity ------------------------------------------------

println("\n" * "="^88)
println("Pass1 ↔ Pass2 gap agreement on the v6_plus sweep (sanity check)")
println("="^88)
println("With rho_0 = |+⟩⟨+|^⊗N parity-broken, Pass 1 (eigenvalues[2]) should agree")
println("with Pass 2 (krylov_spectral_gap) to ~1e-9 or better (qf-e4z.30 finding).\n")
@printf("%-3s %-7s | %-12s %-12s\n", "n", "β_phys", "max rel_diff", "max |Δ|")
println("-"^60)
for k in common_keys
    n, β = k
    c_v6p = v6_plus_cells[k]
    diffs    = Float64[]
    abs_diff = Float64[]
    for r in c_v6p
        haskey(r, :gap_rel_diff_pass12) || continue
        push!(diffs, r.gap_rel_diff_pass12)
        push!(abs_diff, abs(r.gap_arnoldi_pass1 - r.gap_arnoldi_pass2))
    end
    isempty(diffs) && continue
    @printf("%-3d %-7.2f | %-12.3e %-12.3e\n",
            n, β, maximum(diffs), maximum(abs_diff))
end

# --- Even/odd-n diagnostic on v6_plus -------------------------------------

println("\n" * "="^88)
println("Even/odd-n diagnostic on v6_plus (gap_phys medians vs β_phys)")
println("="^88)
println("If v6_plus still shows monotone collapse on even-n at high β, the splitting is")
println("real physics. If v6_plus shows even-n gap_phys flat across β, v6 was tracking")
println("the parity-EVEN sub-spectrum eigenvalue.\n")

by_n_v6p = Dict{Int, Vector{NamedTuple}}()
for k in common_keys
    n, β = k
    c_v6p = v6_plus_cells[k]
    gaps = gap_phys.(c_v6p)
    push!(get!(by_n_v6p, n, NamedTuple[]), (beta_phys=β, gap_med=median(gaps),
        gap_min=minimum(gaps), gap_max=maximum(gaps)))
end

for n in sort!(collect(keys(by_n_v6p)))
    parity = iseven(n) ? "even" : "odd "
    rows = sort(by_n_v6p[n]; by = r -> r.beta_phys)
    βs   = [r.beta_phys for r in rows]
    meds = [r.gap_med   for r in rows]
    # Trend across β: ratio of largest to smallest β
    ratio = meds[end] / max(meds[1], 1e-30)
    @printf("  n=%d (%s)  gap_phys[β_phys=%s]: %s   →  ratio min/max = %.3g\n",
            n, parity,
            join(@sprintf("%.2f", β) for β in βs),
            join(@sprintf(" %.3g", g) for g in meds),
            ratio)
end

# --- Side-by-side: v6 vs v6_plus per-n gap_phys-vs-β --------------------

println("\n" * "="^110)
println("Side-by-side gap_phys medians per n (v6 vs v6_plus)")
println("="^110)
@printf("%-3s %-7s | %-14s %-14s | %-7s\n", "n", "β_phys",
        "v6 (I/d)", "v6+ (|+⟩^⊗N)", "ratio")
println("-"^110)
for n in sort!(collect(intersect(keys(by_n_v6p),
                                  Set(k[1] for k in common_keys))))
    rows = sort(by_n_v6p[n]; by = r -> r.beta_phys)
    for r in rows
        k = (n, r.beta_phys)
        c_v6 = v6_cells[k]
        gap_v6_med = median(gap_phys.(c_v6))
        @printf("%-3d %-7.2f | %-14.4g %-14.4g | %-7.3f\n",
                n, r.beta_phys, gap_v6_med, r.gap_med, r.gap_med / max(gap_v6_med, 1e-30))
    end
    println()
end

# --- Plots ----------------------------------------------------------------

println("[plot] preparing comparison figure…")
ENV["GKSwstype"] = "100"
using Plots

function build_v6_byb()
    by_b = Dict{Float64, Vector{NamedTuple}}()
    for k in common_keys
        c_v6 = v6_cells[k]
        gp = median(gap_phys.(c_v6))
        push!(get!(by_b, k[2], NamedTuple[]), (n=k[1], gap=gp))
    end
    return Dict(b => sort(rs; by = r -> r.n) for (b, rs) in by_b)
end

function build_v6p_byb()
    by_b = Dict{Float64, Vector{NamedTuple}}()
    for k in common_keys
        c_v6p = v6_plus_cells[k]
        gp = median(gap_phys.(c_v6p))
        push!(get!(by_b, k[2], NamedTuple[]), (n=k[1], gap=gp))
    end
    return Dict(b => sort(rs; by = r -> r.n) for (b, rs) in by_b)
end

v6_byb  = build_v6_byb()
v6p_byb = build_v6p_byb()

p = Plots.plot(layout = (1, 2), size = (1300, 540),
    legend = :outerright)

# Left: v6 (I/d, parity-trapped) — same data shown in v6_multiseed_gap_phys.pdf
p1 = p[1]
Plots.plot!(p1, title = "v6: ρ₀ = I/d + krylovdim=40 (parity-trapped)",
    xlabel = "n", ylabel = "median gap_phys",
    yscale = :log10)
for β in sort(collect(keys(v6_byb)))
    rs = v6_byb[β]
    Plots.plot!(p1, [r.n for r in rs], [r.gap for r in rs];
        label = @sprintf("β_phys=%.2f", β), marker = :auto, markersize = 5)
end

# Right: v6_plus (|+⟩⟨+|^⊗N, parity-broken) — true L gap per qf-e4z.30
p2 = p[2]
Plots.plot!(p2, title = "v6_plus: ρ₀ = |+⟩⟨+|^⊗N + krylovdim=60 (true L gap)",
    xlabel = "n", ylabel = "median gap_phys",
    yscale = :log10)
for β in sort(collect(keys(v6p_byb)))
    rs = v6p_byb[β]
    Plots.plot!(p2, [r.n for r in rs], [r.gap for r in rs];
        label = @sprintf("β_phys=%.2f", β), marker = :auto, markersize = 5)
end

fig_path_png = joinpath(FIG_DIR, "v6_plus_vs_v6_gap_phys.png")
fig_path_pdf = joinpath(FIG_DIR, "v6_plus_vs_v6_gap_phys.pdf")
Plots.savefig(p, fig_path_png)
Plots.savefig(p, fig_path_pdf)
@printf("[plot] wrote %s\n", fig_path_png)

println("\n[done] $(now())")
