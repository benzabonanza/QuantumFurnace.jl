#!/usr/bin/env julia
# scratch_dll_vs_ckg_overnight_sweep.jl
#
# Overnight driver for the DLL-vs-CKG comparison plot (qf-e4z.11 / Plot P1).
# Drives sweeps S1 (CKG smooth-Metro EnergyDomain) and S2 (DLL Metropolis
# BohrDomain) cell-by-cell, printing per-cell wall time and process RSS so
# the operator can extrapolate cost for larger n.
#
# Strategy:
#   - For each `n`, run all (β, ε) cells for BOTH S1 and S2 before moving to
#     the next n. That way the comparison plot has matching data points.
#   - `skip_existing=true` so re-runs only fill missing sidecars.
#   - Per ε / per filter, output_dir is unique to avoid sidecar collisions
#     (the Lindbladian sidecar name does NOT include ε / filter_kind).
#   - Walks n = 3, 4, … until budget runs low or scaling makes the next cell
#     too expensive (predicted wall > 0.6 * remaining budget).
#
# Threading:
#   - JULIA_NUM_THREADS=8 to enable Threads.@threads on the omega-loop in
#     `_accumulate_jump_sandwich!` (threshold N_ω ≥ 10, see jump_workers.jl:563).
#   - OPENBLAS_NUM_THREADS=1 to avoid 8 × 8 thread oversubscription on the
#     dense `mul!` calls inside the threaded omega loop.
#
# Usage:
#   JULIA_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 julia --project \
#       scripts/scratch_dll_vs_ckg_overnight_sweep.jl
#
# Wall-time + memory budget: 10h, 4 GB RAM, 8 cores.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using QuantumFurnace
using LinearAlgebra
using Printf
using Dates
using BSON

# --- Process-RSS helper (Linux). Returns RSS in MB. ---
function rss_mb()
    try
        for line in eachline("/proc/self/status")
            if startswith(line, "VmRSS:")
                kb = parse(Int, split(line)[2])
                return kb / 1024.0
            end
        end
        return NaN
    catch
        return NaN
    end
end

function rss_peak_mb()
    try
        for line in eachline("/proc/self/status")
            if startswith(line, "VmHWM:")
                kb = parse(Int, split(line)[2])
                return kb / 1024.0
            end
        end
        return NaN
    catch
        return NaN
    end
end

# Force BLAS to single-thread; outer omega-loop is multi-threaded (Julia threads).
BLAS.set_num_threads(1)
println("[init] Julia threads = ", Threads.nthreads(), ", BLAS threads = ", BLAS.get_num_threads())
println("[init] hostname = ", gethostname(), "   ", now())

const OUTPUT_ROOT = joinpath(@__DIR__, "output")
const PARAM_TABLE  = joinpath(OUTPUT_ROOT, "ideal_lindbladian_param_table.bson")

# Time budget (seconds). Leave 30 min headroom for cleanup.
const T_BUDGET_SEC = 9.5 * 3600.0
const T_START = time()

remaining_budget() = T_BUDGET_SEC - (time() - T_START)

# Per-cell output discriminated by (filter, ε) since the Lindbladian sidecar
# filename does NOT include those tags.
function _cell_output_dir(sweep_tag::AbstractString, filter_kind::Symbol, ε::Real)
    eps_str = @sprintf("%.0e", ε)
    return joinpath(OUTPUT_ROOT, sweep_tag, "$(filter_kind)_eps$(eps_str)")
end

# Bigger t_max_factor than the :auto heuristic — the bi-exp offset C
# only stays below ε if the trajectory itself reaches near ε within
# t_max. We need at least gap·t_max ≳ log(d/ε) to clear the offset
# gate at line 166 of mixing.jl. The :auto heuristic (1.5·log10(1/ε))
# underestimates this for tighter ε.
_t_max_factor(ε::Real) = ε >= 1e-3 ? 12.0 : 22.0
_t_grid_length(ε::Real) = ε >= 1e-3 ? 81 : 121

# Run one cell, print stats, return the per-cell NamedTuple (or `nothing` on err).
function run_one_cell(; sweep_tag, n, β, ε, filter_kind, construction, domain,
                       filter, krylovdim, label)
    out_dir = _cell_output_dir(sweep_tag, filter_kind, ε)
    mkpath(out_dir)

    rss_before = rss_mb()
    GC.gc()
    t0 = time()
    local result
    try
        results = sweep_mixing_times(
            [n], [Float64(β)];
            construction      = construction,
            domain            = domain,
            filter            = filter,
            mode              = :L,
            method            = :krylov,
            seeds             = [42],
            init_state        = :maximally_mixed,
            target_epsilon    = float(ε),
            t_max_factor      = _t_max_factor(ε),
            t_grid_length     = _t_grid_length(ε),
            spectral_krylovdim = krylovdim,
            tol               = 1e-10,
            output_dir        = out_dir,
            hamiltonian_filename = n -> "heis_xxx_zzdisordered_periodic_n$(n).bson",
            use_threads       = false,
            skip_existing     = true,
            param_table_bson  = PARAM_TABLE,
            filter_kind       = filter_kind,
        )
        result = first(results)
    catch err
        elapsed = time() - t0
        rss_now = rss_mb()
        msg = sprint(showerror, err)
        @error "[$label] CELL FAILED" n=n β=β ε=ε filter_kind=filter_kind elapsed_s=elapsed rss_mb=rss_now msg
        return nothing
    end
    elapsed = time() - t0
    rss_after = rss_mb()
    rss_peak = rss_peak_mb()

    # qf-e4y.7: :krylov route emits eigenmode schema (gap_est, floor_distance,
    # mixing_time_source ∈ {:extrapolated, :floor, :nan}) — no fitted_gap /
    # r_squared / converged_fit on this path.
    @printf("[%s] n=%d β=%-4g ε=%.0e filt=%s | τ_mix=%.4g (%s) | gap=%.4g floor=%.4g | matvecs=%d conv=%s | wall=%.2fs ΔRSS=%+.1fMB peak=%.0fMB\n",
            label, n, β, ε, String(filter_kind),
            result.mixing_time, result.mixing_time_source,
            result.gap_est, result.floor_distance,
            result.total_matvecs, result.all_converged,
            elapsed, rss_after - rss_before, rss_peak)
    flush(stdout)
    return result
end

# Run S1 (CKG smooth-Metro Energy) for one n across (β, ε).
function run_S1_at_n(n; krylovdim=60)
    rs = NamedTuple[]
    for β in [5.0, 10.0, 20.0], ε in [1e-3, 1e-6]
        if remaining_budget() < 60.0
            @warn "Budget < 1 min; skipping S1 cell" n=n β=β ε=ε
            return rs
        end
        r = run_one_cell(; sweep_tag="sweep_S1_ckg_ideal",
                          n=n, β=β, ε=ε,
                          filter_kind=:smooth_metro,
                          construction=KMS(), domain=EnergyDomain(),
                          filter=nothing, krylovdim=krylovdim,
                          label="S1·n$(n)")
        r === nothing || push!(rs, r)
    end
    return rs
end

# Run S2 (DLL Metropolis Bohr) for one n across (β, ε).
function run_S2_at_n(n; krylovdim=60)
    rs = NamedTuple[]
    for β in [5.0, 10.0, 20.0], ε in [1e-3, 1e-6]
        if remaining_budget() < 60.0
            @warn "Budget < 1 min; skipping S2 cell" n=n β=β ε=ε
            return rs
        end
        r = run_one_cell(; sweep_tag="sweep_S2_dll_ideal",
                          n=n, β=β, ε=ε,
                          filter_kind=:smooth_metro,    # ignored for DLL
                          construction=DLL(), domain=BohrDomain(),
                          filter=DLLMetropolisFilter(Float64(β)), krylovdim=krylovdim,
                          label="S2·n$(n)")
        r === nothing || push!(rs, r)
    end
    return rs
end

# Compose: do BOTH S1 and S2 at one n.
function run_both_at_n(n; krylovdim=60)
    println("\n" * "="^72)
    println("=== n = $n  (S1 then S2)  remaining budget = $(round(remaining_budget()/3600, digits=2)) h")
    println("="^72)
    rs1 = run_S1_at_n(n; krylovdim=krylovdim)
    rs2 = run_S2_at_n(n; krylovdim=krylovdim)
    return rs1, rs2
end

function summary_log(all_results::Vector{NamedTuple}; tag::AbstractString)
    println("\n" * "─"^72)
    println("Summary [$tag] — $(length(all_results)) cells")
    println("─"^72)
    for r in all_results
        @printf("  %s n=%d β=%-4g ε=%.0e | τ=%-10.4g (%-12s) | gap=%-9.4g | wall=%6.1fs\n",
                r.construction * "/" * r.domain,
                r.n, r.beta, r.target_epsilon,
                r.mixing_time, String(r.mixing_time_source),
                r.gap_est, r.wall_time)
    end
    if !isempty(all_results)
        total = sum(r.wall_time for r in all_results)
        @printf("  TOTAL wall (sum of cells, includes cached): %.1fs (%.2f h)\n", total, total/3600)
    end
    println("─"^72)
end

# ---------------------------------------------------------------------------
# Main: walk n upward, push as far as the budget allows.
# ---------------------------------------------------------------------------
function main()
    println("\n[main] T_BUDGET = $(T_BUDGET_SEC/3600) h, starting at $(now())")

    rs_s1 = NamedTuple[]
    rs_s2 = NamedTuple[]

    # Phase A · n=3..6 — fast, walks both sweeps.
    for n in 3:6
        rs1, rs2 = run_both_at_n(n)
        append!(rs_s1, rs1)
        append!(rs_s2, rs2)
    end
    summary_log(rs_s1; tag="A · S1 (n=3..6)")
    summary_log(rs_s2; tag="A · S2 (n=3..6)")

    # Live wall-time scaling: average of fresh n=k cells (filter cached using
    # whether wall_time was the result of THIS run). We approximate by using
    # all cells; cached ones have small wall_time so they only bias DOWN —
    # safe for triggering early-exit decisions.
    wall_at_n(rs, n) = sum(r.wall_time for r in rs if r.n == n; init=0.0)

    # Walk upward starting at n=7. The per-cell budget guard inside
    # `run_S1_at_n` / `run_S2_at_n` is the real safety net: it skips
    # cells once the budget falls below 60s, so partial-n results are
    # always saved as sidecars. The outer projection just logs an estimate
    # and bails only if the *single-cell* predicted wall would exceed the
    # remaining budget (i.e. we'd never finish even one n=k cell).
    n = 7
    while n <= 12
        if remaining_budget() < 60.0
            @info "Budget < 1 min; stopping push at n=$n" n remaining=remaining_budget()
            break
        end
        prev = wall_at_n(rs_s1, n-1) + wall_at_n(rs_s2, n-1)
        prev_prev = wall_at_n(rs_s1, n-2) + wall_at_n(rs_s2, n-2)
        ratio = if prev_prev > 1.0 && prev > 1.0
            max(prev / prev_prev, 4.0)
        else
            10.0
        end
        # Per-cell estimate: prev / 12 cells × ratio per n-step.
        est_per_cell = prev / 12.0 * ratio
        @printf("[scale] n=%d ratio_est=%.1fx, prev wall (S1+S2)=%.0fs, est per-cell=%.0fs, remaining=%.0fs (%.1fh)\n",
                n, ratio, prev, est_per_cell, remaining_budget(), remaining_budget()/3600)
        flush(stdout)
        if est_per_cell > remaining_budget()
            @info "Even one n=$n cell exceeds remaining budget — stopping" est_per_cell remaining=remaining_budget()
            break
        end
        rs1, rs2 = run_both_at_n(n)
        append!(rs_s1, rs1)
        append!(rs_s2, rs2)
        n += 1
    end

    summary_log(rs_s1; tag="Final S1")
    summary_log(rs_s2; tag="Final S2")

    println("\n[main] All phases complete at $(now())")
    @printf("[main] Total elapsed: %.2f h, peak RSS: %.0f MB\n",
            (time() - T_START) / 3600, rss_peak_mb())
end

main()
