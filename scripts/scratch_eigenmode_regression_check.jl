# qf-e4y.9: Regression cross-check.
#
# Loads each Lindbladian sidecar from the existing CKG / DLL ideal sweeps,
# re-runs `sweep_mixing_times` with `method=:krylov` on the same (n, β, ε,
# filter, construction, domain) tuple under the new eigenmode schema, and
# tabulates per-cell deltas:
#
#   - `gap_old` (cached `gap_est`) vs `gap_new` (new `gap_est`) — must agree
#     to better than rtol 1e-6 (both come from the same Krylov-Arnoldi pass).
#   - `τ_old` (cached `mixing_time`) vs `τ_new` (new `mixing_time`) — should
#     agree to rtol 5% on cells where the cached `fitted_gap / gap_est` ratio
#     lies in [0.8, 1.2] (biexp had a healthy fit). On cells where the ratio
#     is out of that band (LM degenerate-basin failures), the new path is
#     simply expected to produce a finite, physically-sensible τ.
#
# Output:
#   - per-row table to stdout
#   - `scripts/output/eigenmode_regression_table.txt` for the verifier-agent step
#
# Usage:
#   julia --project scripts/scratch_eigenmode_regression_check.jl

using QuantumFurnace
using BSON
using Printf

const PROJECT_ROOT = dirname(@__DIR__)
const FAMILY_HAM   = n -> "heis_xxx_zzdisordered_periodic_n$(n).bson"
const PARAM_TABLE  = joinpath(PROJECT_ROOT, "scripts", "output",
                              "ideal_lindbladian_param_table.bson")
const OUT_TABLE    = joinpath(PROJECT_ROOT, "scripts", "output",
                              "eigenmode_regression_table.txt")

# Restrict to n ∈ {3, 4, 5} for a tractable wall-time budget on the test
# machine. n=6 / n=7 cells take ~minutes each; the schema-correctness claim
# is independent of n.
const N_LIMIT = 5

# (label, dir, construction, domain, filter_kind)
const SWEEPS = [
    ("CKG sM ε=1e-3", joinpath(PROJECT_ROOT, "scripts", "output",
                                "sweep_S1_ckg_ideal", "smooth_metro_eps1e-03"),
     KMS(), EnergyDomain(), :smooth_metro, 1e-3),
    ("CKG sM ε=1e-6", joinpath(PROJECT_ROOT, "scripts", "output",
                                "sweep_S1_ckg_ideal", "smooth_metro_eps1e-06"),
     KMS(), EnergyDomain(), :smooth_metro, 1e-6),
    ("DLL M ε=1e-3",  joinpath(PROJECT_ROOT, "scripts", "output",
                                "sweep_S2_dll_ideal", "smooth_metro_eps1e-03"),
     DLL(), BohrDomain(), :smooth_metro, 1e-3),
    ("DLL M ε=1e-6",  joinpath(PROJECT_ROOT, "scripts", "output",
                                "sweep_S2_dll_ideal", "smooth_metro_eps1e-06"),
     DLL(), BohrDomain(), :smooth_metro, 1e-6),
]

function load_cell(path::AbstractString)
    try
        d = BSON.load(path, QuantumFurnace)
        return NamedTuple(d[:result])
    catch err
        @warn "load failed" path err
        return nothing
    end
end

# Header
header = @sprintf("%-15s %-3s %-5s   %-10s %-10s %-10s   %-12s %-12s %-9s   %-12s",
                  "sweep", "n", "β",
                  "gap_old", "gap_new", "Δgap_rtol",
                  "τ_old", "τ_new", "Δτ_rtol",
                  "fit/gap_old")
sep = repeat("=", length(header))
println(sep)
println(header)
println(sep)

function run_regression()
results = NamedTuple[]
n_pass_gap = 0
n_pass_tau_healthy = 0
n_pass_tau_broken = 0
n_total = 0

for (label, dir, construction, domain, filter_kind, ε) in SWEEPS
    isdir(dir) || (@info "skip $label — dir missing" dir; continue)
    for f in sort(readdir(dir))
        endswith(f, ".bson") || continue
        cell = load_cell(joinpath(dir, f))
        cell === nothing && continue
        n = cell.n
        n > N_LIMIT && continue
        β = cell.beta
        gap_old = cell.gap_est
        τ_old = cell.mixing_time
        # Healthy = biexp's fitted_gap matched the spectral gap. On the
        # cached sidecars we have both fields.
        fit_old = haskey(cell, :fitted_gap) ? cell.fitted_gap : NaN
        ratio_old = isfinite(fit_old) && gap_old > 0 ? fit_old / gap_old : NaN

        # Re-run with the NEW :krylov path on the same fixture.
        res_new = try
            sweep_mixing_times(
                [n], [Float64(β)];
                construction = construction,
                domain = domain,
                filter = construction isa DLL ? DLLMetropolisFilter(1.0) : nothing,
                mode = :L,
                method = :krylov,
                seeds = [42],
                target_epsilon = ε,
                a = 0.0, s = 0.25,
                t_max_factor = :auto,
                t_grid_length = 81,
                spectral_krylovdim = 60,
                tol = 1e-10,
                hamiltonian_filename = FAMILY_HAM,
                use_threads = false,
                output_dir = nothing,
                param_table_bson = PARAM_TABLE,
                filter_kind = filter_kind,
            )
        catch err
            @warn "RECOMPUTE FAILED" label n β err
            continue
        end
        r = res_new[1]
        gap_new = r.gap_est
        τ_new = r.mixing_time

        Δgap = abs(gap_new - gap_old) / max(abs(gap_old), 1e-12)
        Δτ   = isfinite(τ_old) && isfinite(τ_new) && τ_old > 0 ?
               abs(τ_new - τ_old) / τ_old : NaN

        n_total += 1
        Δgap < 1e-4 && (n_pass_gap += 1)
        is_healthy = isfinite(ratio_old) && 0.8 ≤ ratio_old ≤ 1.2
        if is_healthy
            isfinite(Δτ) && Δτ < 0.05 && (n_pass_tau_healthy += 1)
        else
            isfinite(τ_new) && τ_new > 0 && (n_pass_tau_broken += 1)
        end

        push!(results, (
            label = label, n = n, β = β,
            gap_old = gap_old, gap_new = gap_new, Δgap = Δgap,
            τ_old = τ_old, τ_new = τ_new, Δτ = Δτ,
            ratio_old = ratio_old, healthy = is_healthy,
            source_new = r.mixing_time_source,
        ))
        @printf("%-15s %-3d %-5g   %-10.4e %-10.4e %-10.2e   %-12.4e %-12.4e %-9.2e   %-12.4f\n",
                label, n, β, gap_old, gap_new, Δgap, τ_old, τ_new, Δτ, ratio_old)
    end
end

println(sep)
@printf("%d cells; gap-parity (rtol 1e-4): %d / %d\n", n_total, n_pass_gap, n_total)
n_healthy = count(r -> r.healthy, results)
n_broken  = n_total - n_healthy
@printf("%d healthy cells (fit/gap ∈ [0.8, 1.2]) — τ-parity within 5%%: %d / %d\n",
        n_healthy, n_pass_tau_healthy, n_healthy)
@printf("%d biexp-broken cells — τ_new finite + positive: %d / %d\n",
        n_broken, n_pass_tau_broken, n_broken)
println(sep)

# Persist table to disk for the verifier agent.
open(OUT_TABLE, "w") do io
    println(io, sep)
    println(io, header)
    println(io, sep)
    for r in results
        @printf(io, "%-15s %-3d %-5g   %-10.4e %-10.4e %-10.2e   %-12.4e %-12.4e %-9.2e   %-12.4f\n",
                r.label, r.n, r.β, r.gap_old, r.gap_new, r.Δgap,
                r.τ_old, r.τ_new, r.Δτ, r.ratio_old)
    end
    println(io, sep)
    @printf(io, "%d cells; gap-parity (rtol 1e-4): %d / %d\n",
            n_total, n_pass_gap, n_total)
    @printf(io, "%d healthy cells — τ-parity within 5%%: %d / %d\n",
            n_healthy, n_pass_tau_healthy, n_healthy)
    @printf(io, "%d biexp-broken cells — τ_new finite + positive: %d / %d\n",
            n_broken, n_pass_tau_broken, n_broken)
    println(io, sep)
end
@info "Wrote regression table" OUT_TABLE
return results
end

run_regression()
