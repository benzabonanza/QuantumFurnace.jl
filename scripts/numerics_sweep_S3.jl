#!/usr/bin/env julia
# numerics_sweep_S3.jl  (qf-e4z.5 / qf-6vr)
#
# Sweep S3 — implemented δ-channel, smooth-Metro KMS, TrotterDomain + GQSP,
# jump-sweep splitting. Feeds plots P5 (channel-vs-ideal), P6 (Ham-sim-time
# scaling), P7 (filter-family smooth-Metro arm).
#
# Cells: family = 1D-XXX-zzdis, β_phys ∈ {1, 2, 3}, ε = 1e-3, filter = smooth_metro.
#   Phase A: n ∈ {3, 4, 5, 6}        — laptop smoke test.
#   Phase B: n ∈ {7, 8, 9} incremental — bail at first cell exceeding 1 h.
#
# PHYSICS CHECK (qf-6vr / Phase qf-bphys): the β grid is now in *physical*
# inverse-temperature units (against the un-rescaled Hamiltonian). The sweep
# harness loads each fixture's `rescaling_factor` and derives the algorithm-
# side `β_alg = β_phys · rescaling_factor` per cell. Sidecar filenames carry
# the `betaphys<β_phys>` prefix and the sidecar dict records `:beta_phys`,
# `:beta_alg`, `:rescaling_factor`. Required param-table cells must cover the
# derived β_alg ranges — re-run `numerics_param_table.jl` if you extend the
# `n` range outside the current calibration envelope.
#
# Per-cell parameters come from scripts/output/channel_param_table.bson; the
# δ entry there is what gets retuned after the floor audit.
#
# Threading recipe per .claude/rules/julia-code.md "Threading Defaults":
#   JULIA_NUM_THREADS = max (Julia=8 typical), BLAS=1.
#
# Usage (run from repo root):
#   JULIA_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 \
#       julia --project scripts/numerics_sweep_S3.jl
#
# After the run, audit floors with:
#   julia --project scripts/numerics_S3_audit_floors.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using QuantumFurnace
using LinearAlgebra
using Dates

BLAS.set_num_threads(1)
@assert BLAS.get_num_threads() == 1
println("[init] $(now())  Julia threads = $(Threads.nthreads()), BLAS = $(BLAS.get_num_threads())")
println("[init] hostname = $(gethostname())")

# --- Constants ------------------------------------------------------------

const OUTPUT_DIR    = joinpath(@__DIR__, "output", "sweep_S3_channel_betaphys")
const BETAS_PHYS    = [1.0, 2.0, 3.0]     # qf-6vr: β_phys grid (replaces legacy [5, 10, 20] β_alg)
const EPSILONS      = [1e-3]              # ε=1e-3 only per qf-e4z.5 decision 2026-05-08
const FILTER_KINDS  = [:smooth_metro]
const PHASE_A_N     = collect(3:6)
const PHASE_B_N     = [7, 8, 9]           # tried incrementally; bail at first >1h cell
const CELL_BUDGET_S = 3600.0              # per-cell wall-time cap for Phase B continuation

mkpath(OUTPUT_DIR)

# --- Helpers --------------------------------------------------------------

function _run_block(n_values::AbstractVector{Int}; tag::AbstractString)
    println("\n[block] $(tag)  n ∈ $(collect(n_values))  $(now())")
    t0 = time()
    results = sweep_channel_mixing(
        n_values;
        beta_phys_values = BETAS_PHYS,    # qf-6vr β_phys-first sweep
        target_epsilons  = EPSILONS,
        filter_kinds     = FILTER_KINDS,
        domain           = TrotterDomain(),
        construction     = KMS(),
        output_dir       = OUTPUT_DIR,
        skip_existing    = true,
    )
    elapsed = time() - t0
    println("[block] $(tag) done in $(round(elapsed; digits=1)) s — $(length(results)) cells")
    return results
end

function _max_fresh_wall(results)
    # wall_time_seconds reflects the most recent run; for skip_existing-loaded
    # cells it carries the ORIGINAL wall time. That's still useful for the
    # Phase-B-continuation decision: if a cell ever took > 1 h, the higher-n
    # extrapolation isn't going to be cheaper.
    walls = Float64[]
    for r in results
        w = r.wall_time_seconds
        (isnan(w) || w <= 0) && continue
        push!(walls, w)
    end
    isempty(walls) && return 0.0
    return maximum(walls)
end

function _summarise(results)
    n_ext   = count(r -> r.tau_mix_source === :extrapolated, results)
    n_floor = count(r -> r.tau_mix_source === :floor,        results)
    n_nan   = count(r -> r.tau_mix_source === :nan,          results)
    println("[block] tau_mix_source: $(n_ext) :extrapolated, $(n_floor) :floor, $(n_nan) :nan")
end

# --- Phase A --------------------------------------------------------------

println("\n=== Phase A: n ∈ $(PHASE_A_N) (β_phys = $(BETAS_PHYS); laptop ≤ 30 min budget) ===")
results_a = _run_block(PHASE_A_N; tag="A")
_summarise(results_a)
println("[block] Phase A max cell wall = $(round(_max_fresh_wall(results_a); digits=1)) s")

# --- Phase B (incremental) ------------------------------------------------

println("\n=== Phase B: incremental n ∈ $(PHASE_B_N), bail at first cell > $(CELL_BUDGET_S) s ===")
for n in PHASE_B_N
    results_n = _run_block([n]; tag="B-n$(n)")
    _summarise(results_n)
    w_max = _max_fresh_wall(results_n)
    println("[block] n=$(n) max cell wall = $(round(w_max; digits=1)) s (budget $(CELL_BUDGET_S) s)")
    if w_max > CELL_BUDGET_S
        @warn "Cell wall exceeded budget; halting Phase B before higher n" n max_wall_seconds=w_max budget_seconds=CELL_BUDGET_S
        break
    end
end

println("\n=== S3 sweep done $(now()) ===")
println("Next: julia --project scripts/numerics_S3_audit_floors.jl")
