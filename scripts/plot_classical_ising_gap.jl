#!/usr/bin/env julia
# plot_classical_ising_gap.jl  (qf-8fr)
#
# Plot the classical Ising sweep results: λ_gap vs n at β_phys = 0.5
# Two panels:
#   (a) gap_alg vs n          — what the code outputs (rescaled-frame gap)
#   (b) gap_phys = gap_alg × rescaling_factor vs n   — physical-frame gap
# A "constant gap" reference line is drawn on the gap_phys panel to anchor
# the Ω(1) prediction from KB16 / Bardet23 / BC25.
#
# Disordered Heisenberg β_phys=0.5 (qf-e4z.21) is overlaid in dashed lines
# to highlight that the "n^{-1.2}" gap_alg decay reported for Heisenberg
# is the SAME rescaling-factor artefact: in gap_phys both curves are roughly
# flat in n.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BSON
using Printf
using Statistics
using QuantumFurnace
ENV["GKSwstype"] = "100"
using Plots

const REPO_ROOT = joinpath(@__DIR__, "..")
const ISING_DIR = joinpath(REPO_ROOT, "scripts", "output", "sweep_classical_ising")
const HEIS_DIR  = joinpath(REPO_ROOT, "scripts", "output", "sweep_S1_v4_ckg_ideal",
                           "smooth_metro_eps1e-03")
const OUT_DIR   = joinpath(REPO_ROOT, "drafts", "plots")
mkpath(OUT_DIR)

function loadrows(dir::AbstractString)
    rows = NamedTuple[]
    isdir(dir) || return rows
    for p in sort(readdir(dir; join=true))
        endswith(p, ".bson") || continue
        startswith(basename(p), "sweep_classical_ising_summary") && continue
        try
            d = BSON.load(p, QuantumFurnace)
            haskey(d, :result) || continue
            r = d[:result]
            if r isa Dict
                push!(rows, (; (Symbol(k) => v for (k, v) in pairs(r))...))
            elseif r isa NamedTuple
                push!(rows, r)
            end
        catch err
            @warn "BSON load failed; skipping" path=p err=err
        end
    end
    return rows
end

# Filter EnergyDomain rows, optionally by β_phys.  Heisenberg P1 v4 sidecars
# write `domain = "Energy"`; the classical Ising sweep writes `"EnergyDomain"`.
function filter_rows(rows; domain_alts = ("EnergyDomain", "Energy"),
                     β_phys::Union{Real, Nothing} = nothing)
    out = filter(r -> haskey(r, :domain) ? r.domain in domain_alts : true, rows)
    if β_phys !== nothing
        out = filter(r -> haskey(r, :beta_phys) && r.beta_phys ≈ β_phys, out)
    end
    sort!(out, by = r -> r.n)
    return out
end

function loglog_slope(xs, ys)
    length(xs) >= 2 || return NaN
    lx = log.(xs); ly = log.(ys)
    xm = mean(lx); ym = mean(ly)
    cov = sum((lx .- xm) .* (ly .- ym))
    var = sum((lx .- xm) .^ 2)
    return cov / var
end

# --- Load data ---
println("Loading classical Ising sweep from $ISING_DIR ...")
rows_ising = loadrows(ISING_DIR)
ising = filter_rows(rows_ising)
println("  found $(length(ising)) cells")

println("Loading Heisenberg P1 v4 β_phys=0.5 from $HEIS_DIR ...")
rows_heis = loadrows(HEIS_DIR)
heis = filter_rows(rows_heis; β_phys = 0.5)
println("  found $(length(heis)) cells")

# Extract series
function series(rows)
    n  = [r.n for r in rows]
    rescale = [haskey(r, :rescaling_factor) ? r.rescaling_factor : NaN for r in rows]
    g_alg = [haskey(r, :gap) ? r.gap : (haskey(r, :gap_arnoldi) ? r.gap_arnoldi : NaN) for r in rows]
    g_phys = g_alg .* rescale
    return (; n, gap_alg = g_alg, gap_phys = g_phys, rescale)
end
s_ising = series(ising)
s_heis  = series(heis)

# --- Plot ---
default(size=(960, 420), titlefontsize=11, legendfontsize=8, guidefontsize=10,
        tickfontsize=8, framestyle=:box, dpi=150)

# Panel (a): gap_alg vs n (log-log)
p_alg = plot(xscale=:log10, yscale=:log10, xlabel="n", ylabel="λ_gap (algorithm frame)",
    title="(a) gap_alg vs n  (β_phys = 0.5)", legend=:bottomleft)
if !isempty(s_ising.n)
    sl = loglog_slope(s_ising.n, s_ising.gap_alg)
    plot!(p_alg, s_ising.n, s_ising.gap_alg; marker=:circle, color="#7A2E39",
        label=@sprintf("Classical Ising  (slope %.2f)", sl), linewidth=2)
end
if !isempty(s_heis.n)
    sl = loglog_slope(s_heis.n, s_heis.gap_alg)
    plot!(p_alg, s_heis.n, s_heis.gap_alg; marker=:square, color="#4F3B5C",
        linestyle=:dash, label=@sprintf("Heisenberg disordered  (slope %.2f)", sl), linewidth=2)
end
# Reference n^{-1} line for visual
if !isempty(s_ising.n)
    nref = float.(s_ising.n)
    base = s_ising.gap_alg[1] * float(s_ising.n[1])
    plot!(p_alg, nref, base ./ nref; linestyle=:dot, color=:gray, label="n⁻¹ reference")
end

# Panel (b): gap_phys vs n  (log scale on y for clarity but linear could also work)
p_phys = plot(xscale=:log10, yscale=:log10, xlabel="n", ylabel="λ_gap × rescaling (physical frame)",
    title="(b) gap_phys vs n  (β_phys = 0.5)", legend=:bottomleft)
if !isempty(s_ising.n)
    sl_full = loglog_slope(s_ising.n, s_ising.gap_phys)
    sl_no3  = length(s_ising.n) >= 3 ? loglog_slope(s_ising.n[2:end], s_ising.gap_phys[2:end]) : NaN
    plot!(p_phys, s_ising.n, s_ising.gap_phys; marker=:circle, color="#7A2E39",
        label=@sprintf("Classical Ising  (slope %.2f / n≥4: %.2f)", sl_full, sl_no3), linewidth=2)
end
if !isempty(s_heis.n)
    sl_full = loglog_slope(s_heis.n, s_heis.gap_phys)
    sl_no3  = length(s_heis.n) >= 3 ? loglog_slope(s_heis.n[2:end], s_heis.gap_phys[2:end]) : NaN
    plot!(p_phys, s_heis.n, s_heis.gap_phys; marker=:square, color="#4F3B5C",
        linestyle=:dash, label=@sprintf("Heisenberg disordered  (slope %.2f / n≥4: %.2f)", sl_full, sl_no3),
        linewidth=2)
end
# Constant-gap reference at the mean of the n≥4 Ising data
if length(s_ising.n) >= 3
    const_g = mean(s_ising.gap_phys[2:end])
    hline!(p_phys, [const_g]; linestyle=:dot, color=:gray,
        label=@sprintf("constant (mean Ising n≥4: %.2f)", const_g))
end

final = plot(p_alg, p_phys; layout=(1, 2))

savefig(final, joinpath(OUT_DIR, "classical_ising_gap.png"))
savefig(final, joinpath(OUT_DIR, "classical_ising_gap.pdf"))
@printf("\nWrote %s and .pdf\n", joinpath(OUT_DIR, "classical_ising_gap.png"))
