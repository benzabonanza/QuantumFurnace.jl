#!/usr/bin/env julia
# Test: does τ_mix(1e-6) = τ_mix(1e-3) + (3 ln10)/λ across the existing
# CKG / DLL ideal-Lindbladian sweeps? (Answers user's question.)

using BSON
using Printf
using Statistics

const LN10 = log(10.0)
const OFFSET_FACTOR = 3 * LN10  # τ_mix(1e-6) - τ_mix(1e-3) = 3 ln10 / λ

const SWEEPS = [
    ("CKG sM",  "scripts/output/sweep_S1_ckg_ideal/smooth_metro_eps1e-03",
                "scripts/output/sweep_S1_ckg_ideal/smooth_metro_eps1e-06"),
    ("DLL Metro", "scripts/output/sweep_S2_dll_ideal/smooth_metro_eps1e-03",
                  "scripts/output/sweep_S2_dll_ideal/smooth_metro_eps1e-06"),
]

function load_cell(path::AbstractString)
    haskey(BSON.load(path), :result) || return nothing
    return BSON.load(path)[:result]
end

function pair_files(dir3, dir6)
    files3 = sort(filter(s -> endswith(s, ".bson"), readdir(dir3)))
    files6 = sort(filter(s -> endswith(s, ".bson"), readdir(dir6)))
    pairs = NamedTuple[]
    for f3 in files3
        f6 = f3  # same naming convention
        if f6 in files6
            push!(pairs, (f3 = joinpath(dir3, f3), f6 = joinpath(dir6, f6), name = f3))
        end
    end
    return pairs
end

function relerr(pred, actual)
    return abs(pred - actual) / actual
end

println(repeat("=", 92))
@printf("%-9s %-3s %-5s   %-10s %-10s %-10s %-10s %-10s   %-10s %-7s\n",
    "sweep", "n", "β",
    "τ(1e-3)", "τ(1e-6)", "λ_fit (3)", "λ_fit (6)",
    "Δτ_meas", "Δτ_pred", "rel.err")
println(repeat("=", 92))

all_results = NamedTuple[]

for (label, dir3, dir6) in SWEEPS
    isdir(dir3) || continue
    isdir(dir6) || continue
    for p in pair_files(dir3, dir6)
        c3 = load_cell(p.f3)
        c6 = load_cell(p.f6)
        (c3 === nothing || c6 === nothing) && continue
        n = c3[:n]; β = c3[:beta]
        τ3 = c3[:mixing_time]; τ6 = c6[:mixing_time]
        # qf-e4y.7: :krylov route sidecars no longer contain fitted_gap.
        # Fall back to gap_est (which agrees with the eigenmode helper to
        # ~1e-12) for any sidecar that lacks the biexp diagnostics.
        λ3 = haskey(c3, :fitted_gap) ? c3[:fitted_gap] : c3[:gap_est]
        λ6 = haskey(c6, :fitted_gap) ? c6[:fitted_gap] : c6[:gap_est]
        gap_est3 = c3[:gap_est]; gap_est6 = c6[:gap_est]
        # Use the spectral-gap pre-pass (gap_est) for the prediction — it is the
        # cleanest λ; on the legacy :ode route fitted_gap is a curve fit
        # and may be noisier at small Δ.
        Δmeas = τ6 - τ3
        Δpred_fitted = OFFSET_FACTOR / λ3            # use the τ(1e-3)-cell's fit
        Δpred_gap    = OFFSET_FACTOR / gap_est3      # use the spectral gap
        rel_fit  = relerr(Δpred_fitted, Δmeas)
        rel_gap  = relerr(Δpred_gap,    Δmeas)
        push!(all_results, (
            sweep=label, n=n, β=β,
            τ3=τ3, τ6=τ6, λ_fit3=λ3, λ_fit6=λ6,
            gap_est3=gap_est3, gap_est6=gap_est6,
            Δmeas=Δmeas, Δpred_fit=Δpred_fitted, Δpred_gap=Δpred_gap,
            rel_fit=rel_fit, rel_gap=rel_gap,
            src3=c3[:mixing_time_source], src6=c6[:mixing_time_source],
        ))
        @printf("%-9s %-3d %-5g   %-10.3f %-10.3f %-10.5f %-10.5f %-10.3f   %-10.3f %-7.2f%%\n",
            label, n, β,
            τ3, τ6, λ3, λ6, Δmeas, Δpred_fitted, 100 * rel_fit)
    end
end

println(repeat("=", 92))
println()
println("Summary — relative error of the prediction Δτ ≈ 3 ln(10) / λ:")
println()
@printf("%-12s %-12s %-12s %-12s %-12s\n",
    "predictor", "median", "mean", "p90", "max")
for (key, label) in [(:rel_fit, "fitted_gap (3)"), (:rel_gap, "gap_est (3)")]
    vals = [r[key] for r in all_results if isfinite(r[key])]
    @printf("%-12s %-12.3f %-12.3f %-12.3f %-12.3f\n",
        label, median(vals), mean(vals), quantile(vals, 0.9), maximum(vals))
end

println()
println("Per-sweep breakdown (median rel.err):")
for sw in unique(r.sweep for r in all_results)
    sub = [r for r in all_results if r.sweep == sw]
    medfit = median([r.rel_fit for r in sub])
    medgap = median([r.rel_gap for r in sub])
    @printf("  %-10s  fit-based: %.2f%%   gap-based: %.2f%%   (n cells = %d)\n",
        sw, 100*medfit, 100*medgap, length(sub))
end

# Check if the offset is uniform across (n, β) — that's the real test of the
# "single-gap exponential tail" assumption.
println()
println("Δτ_meas · λ_fit (should be ≈ 3 ln 10 ≈ 6.908 if gap-limited):")
println()
@printf("%-9s %-3s %-5s   %-10s %-10s %-10s\n",
    "sweep", "n", "β", "Δτ·λ_fit", "Δτ·gap_est", "vs 6.908")
println(repeat("-", 60))
for r in all_results
    prod_fit = r.Δmeas * r.λ_fit3
    prod_gap = r.Δmeas * r.gap_est3
    @printf("%-9s %-3d %-5g   %-10.4f %-10.4f %+8.2f%%\n",
        r.sweep, r.n, r.β, prod_fit, prod_gap, 100*(prod_fit/OFFSET_FACTOR - 1))
end
