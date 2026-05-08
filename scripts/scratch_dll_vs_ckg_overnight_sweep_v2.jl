#!/usr/bin/env julia
# scratch_dll_vs_ckg_overnight_sweep_v2.jl
#
# Faster overnight driver for the DLL-vs-CKG comparison plot. v1 used
# `sweep_mixing_times` which calls `krylov_spectral_gap` separately before
# `predict_lindbladian_trajectory` — both build their own Arnoldi
# factorisation, doubling the per-cell wall. v2 calls
# `predict_lindbladian_trajectory` directly: it already returns
# `spectral_gap` from its dense Hessenberg eigendecomposition, so we
# skip the redundant gap call.
#
# Saves to the same sidecar tree as v1 with the same record schema, so
# previously-saved cells are skip_existing-loaded and analysis scripts
# don't need to know which version produced each cell.
#
# Usage:
#   JULIA_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 julia --project \
#       scripts/scratch_dll_vs_ckg_overnight_sweep_v2.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using QuantumFurnace
using LinearAlgebra
using BSON
using Printf
using Dates

# --- Process-RSS helper (Linux). ---
function rss_mb()
    try
        for line in eachline("/proc/self/status")
            if startswith(line, "VmRSS:")
                return parse(Int, split(line)[2]) / 1024.0
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
                return parse(Int, split(line)[2]) / 1024.0
            end
        end
        return NaN
    catch
        return NaN
    end
end

BLAS.set_num_threads(1)
println("[init] Julia threads = ", Threads.nthreads(), ", BLAS threads = ", BLAS.get_num_threads())
println("[init] hostname = ", gethostname(), "   ", now())

const OUTPUT_ROOT = joinpath(@__DIR__, "output")
const PARAM_TABLE = joinpath(OUTPUT_ROOT, "ideal_lindbladian_param_table.bson")
const T_BUDGET_SEC = 9.5 * 3600.0
const T_START = time()

remaining_budget() = T_BUDGET_SEC - (time() - T_START)

function _cell_output_dir(sweep_tag::AbstractString, filter_kind::Symbol, ε::Real)
    eps_str = @sprintf("%.0e", ε)
    return joinpath(OUTPUT_ROOT, sweep_tag, "$(filter_kind)_eps$(eps_str)")
end

# Same sidecar filename as sweep_mixing_times so existing files load via skip_existing.
function _sidecar_path(out_dir::AbstractString, n::Integer, beta::Real, seed::Integer,
                       construction_tag::AbstractString, domain_tag::AbstractString)
    beta_str = let s = @sprintf("%.6f", float(beta))
        s = rstrip(s, '0'); s = rstrip(s, '.')
        isempty(s) ? "0" : s
    end
    return joinpath(out_dir, "sweep_n$(n)_beta$(beta_str)_seed$(seed)_L_$(construction_tag)_$(domain_tag).bson")
end

# Parameter-table lookup (CKG/EnergyDomain only).
function _ckg_param_row(n::Integer, beta::Real, eps::Real)
    rows = BSON.load(PARAM_TABLE, QuantumFurnace)[:rows]
    for r in rows
        r[:n] == n && r[:beta] ≈ beta && r[:eps] ≈ eps && r[:filter] === :smooth_metro && return r
    end
    error("no param row for n=$n β=$beta ε=$eps smooth_metro")
end

# Build a Lindbladian Config for CKG/EnergyDomain, smooth-Metro.
function _build_ckg_config(n::Integer, beta::Real, row)
    return Config(
        sim = Lindbladian(),
        domain = EnergyDomain(),
        construction = KMS(),
        num_qubits = n,
        with_linear_combination = true,
        beta = float(beta),
        sigma = 1.0 / float(beta),
        a = 0.0, s = 0.25,
        gaussian_parameters = (nothing, nothing),
        num_energy_bits = row[:r_D],
        w0 = row[:w0_D], t0 = row[:t0_D],
        num_trotter_steps_per_t0 = 10,
        filter = nothing,
    )
end

# Build a Lindbladian Config for DLL/BohrDomain.
function _build_dll_config(n::Integer, beta::Real)
    return Config(
        sim = Lindbladian(),
        domain = BohrDomain(),
        construction = DLL(),
        num_qubits = n,
        with_linear_combination = true,
        beta = float(beta),
        sigma = 1.0 / float(beta),
        a = 0.0, s = 0.25,
        gaussian_parameters = (nothing, nothing),
        num_energy_bits = 12, w0 = 0.05, t0 = 2π / (2^12 * 0.05),
        num_trotter_steps_per_t0 = 10,
        filter = DLLMetropolisFilter(float(beta)),
    )
end

# Run one cell: predict_lindbladian_trajectory → bi-exp fit → τ_mix.
function run_one_cell(; sweep_tag, n, β, ε, filter_kind, construction_tag, domain_tag,
                       config, ham, jumps, krylovdim, label)
    out_dir = _cell_output_dir(sweep_tag, filter_kind, ε)
    mkpath(out_dir)
    sidecar = _sidecar_path(out_dir, n, β, 42, construction_tag, domain_tag)

    if isfile(sidecar)
        # Skip-existing: load record and report.
        try
            d = BSON.load(sidecar, QuantumFurnace)
            r = NamedTuple(d[:result])
            @printf("[%s] n=%d β=%-4g ε=%.0e filt=%s | CACHED τ_mix=%.4g (%s) | gap=%.4g | wall=%.1fs\n",
                    label, n, β, ε, String(filter_kind),
                    r.mixing_time, r.mixing_time_source, r.gap_est, r.wall_time)
            flush(stdout)
            return r
        catch err
            @warn "Failed to load cached sidecar; recomputing" sidecar err
        end
    end

    rss_before = rss_mb()
    GC.gc()
    t0_run = time()

    # Generous t_max: gap is unknown a priori. From observed n=3..7 data,
    # gap ≥ 0.05 across all cells. For ε=1e-6 we need gap·t_max ≥ 22 ⇒
    # t_max ≥ 440. Use 600 (ε=1e-6) / 300 (ε=1e-3). Trajectory plateaus
    # are absorbed by the bi-exp offset.
    t_max = ε >= 1e-3 ? 300.0 : 600.0
    t_grid_length = ε >= 1e-3 ? 81 : 121
    t_grid = collect(range(0.0, t_max, length=t_grid_length))

    d_dim = size(ham.data, 1)
    rho_0 = Matrix{ComplexF64}(I, d_dim, d_dim) ./ d_dim    # maximally mixed

    local result
    try
        traj = predict_lindbladian_trajectory(
            config, ham, jumps, rho_0, t_grid;
            krylovdim = krylovdim, tol = 1e-10)
        gap_est = traj.spectral_gap > 0 ? traj.spectral_gap : 1.0 / β

        # qf-e4y.7: eigenmode τ_mix in place of biexp curve fit. Same
        # schema as sweep_mixing_times :krylov route now emits.
        res_eig = eigenmode_mixing_time(
            traj.eigenvalues, traj.c, traj.R_modes,
            traj.rho_inf, traj.sigma_beta, float(ε);
            t_upper = t_max,
        )
        mixing_time = res_eig.mixing_time
        mixing_time_source = res_eig.source

        wall = time() - t0_run
        result = (
            n = n, beta = float(β), seed = 42, init_state = :maximally_mixed,
            mode = :L, method = :krylov,
            construction = construction_tag, domain = domain_tag,
            filter_name = construction_tag == "DLL" ? "DLLMetropolis" : "default",
            target_epsilon = float(ε), filter_kind = filter_kind,
            r_D = config.num_energy_bits, w0_D = config.w0, t0_D = config.t0,
            gap_est = gap_est,
            t_max = t_max,
            t_max_factor = t_max * gap_est,    # post-hoc factor in e-folds
            tau_mix_bound = log(d_dim / float(ε)) / max(gap_est, 1e-12),
            n_grid = t_grid_length,
            total_matvecs = traj.total_matvecs,
            all_converged = traj.all_converged,
            mixing_time = mixing_time,
            mixing_time_source = mixing_time_source,
            floor_distance = res_eig.floor_distance,
            wall_time = wall,
        )

        try
            BSON.bson(sidecar, Dict(:result => Dict(pairs(result))))
        catch err
            @warn "Sidecar write failed (continuing)" sidecar err
        end
    catch err
        wall = time() - t0_run
        rss_now = rss_mb()
        @error "[$label] CELL FAILED" n=n β=β ε=ε wall_s=wall rss_mb=rss_now exception=err
        return nothing
    end

    rss_after = rss_mb()
    rss_peak = rss_peak_mb()

    @printf("[%s] n=%d β=%-4g ε=%.0e filt=%s | τ_mix=%.4g (%s) | gap=%.4g floor=%.4g | matvecs=%d conv=%s | wall=%.2fs ΔRSS=%+.1fMB peak=%.0fMB\n",
            label, n, β, ε, String(filter_kind),
            result.mixing_time, result.mixing_time_source,
            result.gap_est, result.floor_distance,
            result.total_matvecs, result.all_converged,
            result.wall_time, rss_after - rss_before, rss_peak)
    flush(stdout)
    return result
end

function run_S1_at_n(n; krylovdim=40)
    rs = NamedTuple[]
    for β in [5.0, 10.0, 20.0], ε in [1e-3, 1e-6]
        if remaining_budget() < 60.0
            @warn "Budget < 1 min; skipping S1 cell" n=n β=β ε=ε
            return rs
        end
        ham = QuantumFurnace._load_hamiltonian_bson(
            "hamiltonians/heis_xxx_zzdisordered_periodic_n$(n).bson", float(β))
        jumps = QuantumFurnace._build_jump_set(ham, n)
        row = _ckg_param_row(n, float(β), float(ε))
        cfg = _build_ckg_config(n, float(β), row)
        r = run_one_cell(; sweep_tag="sweep_S1_ckg_ideal",
                          n=n, β=float(β), ε=float(ε),
                          filter_kind=:smooth_metro,
                          construction_tag="KMS", domain_tag="Energy",
                          config=cfg, ham=ham, jumps=jumps,
                          krylovdim=krylovdim, label="S1·n$(n)")
        r === nothing || push!(rs, r)
    end
    return rs
end

function run_S2_at_n(n; krylovdim=40)
    rs = NamedTuple[]
    for β in [5.0, 10.0, 20.0], ε in [1e-3, 1e-6]
        if remaining_budget() < 60.0
            @warn "Budget < 1 min; skipping S2 cell" n=n β=β ε=ε
            return rs
        end
        ham = QuantumFurnace._load_hamiltonian_bson(
            "hamiltonians/heis_xxx_zzdisordered_periodic_n$(n).bson", float(β))
        jumps = QuantumFurnace._build_jump_set(ham, n)
        cfg = _build_dll_config(n, float(β))
        r = run_one_cell(; sweep_tag="sweep_S2_dll_ideal",
                          n=n, β=float(β), ε=float(ε),
                          filter_kind=:smooth_metro,
                          construction_tag="DLL", domain_tag="Bohr",
                          config=cfg, ham=ham, jumps=jumps,
                          krylovdim=krylovdim, label="S2·n$(n)")
        r === nothing || push!(rs, r)
    end
    return rs
end

function run_both_at_n(n; krylovdim=40)
    println("\n" * "="^72)
    println("=== n = $n  (S1 then S2)  remaining = $(round(remaining_budget()/3600, digits=2))h")
    println("="^72)
    rs1 = run_S1_at_n(n; krylovdim=krylovdim)
    rs2 = run_S2_at_n(n; krylovdim=krylovdim)
    return rs1, rs2
end

function main()
    println("\n[main] T_BUDGET = $(T_BUDGET_SEC/3600) h, starting at $(now())")

    rs_s1 = NamedTuple[]
    rs_s2 = NamedTuple[]

    for n in 3:6
        rs1, rs2 = run_both_at_n(n)
        append!(rs_s1, rs1); append!(rs_s2, rs2)
    end

    wall_at_n(rs, n) = sum(r.wall_time for r in rs if r.n == n; init=0.0)

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
        est_per_cell = prev / 12.0 * ratio
        @printf("[scale] n=%d ratio=%.1fx prev=%.0fs est_per_cell=%.0fs remaining=%.1fh\n",
                n, ratio, prev, est_per_cell, remaining_budget()/3600)
        flush(stdout)
        if est_per_cell > remaining_budget()
            @info "Single n=$n cell exceeds budget — stopping" est_per_cell remaining=remaining_budget()
            break
        end
        rs1, rs2 = run_both_at_n(n)
        append!(rs_s1, rs1); append!(rs_s2, rs2)
        n += 1
    end

    println("\n[main] Done at $(now())  elapsed=$(round((time()-T_START)/3600, digits=2))h  peak_rss=$(round(rss_peak_mb()))MB")
end

main()
