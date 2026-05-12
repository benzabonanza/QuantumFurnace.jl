#!/usr/bin/env julia
#
# qf-w2u — How much harder is the dissipator quadrature at β_phys=2 vs β_phys=1?
#
# Smooth Metropolis (default_smooth_s), n ∈ {3, 4, 5, 6}, matrix-free both
# sides:
#
#   ‖ apply_lindbladian(EnergyDomain, r_test)  -  apply_lindbladian(EnergyDomain, r_ref) ‖_op
#
# The reference at r_ref = R_REF takes the place of the dense BohrDomain
# reference used in `scratch_quad_S1_krylov_campaign.jl`. This is sound when
# the EnergyDomain ↔ BohrDomain operator-norm error at r_ref is well below
# the smallest target ε (1e-9 here): then L_E(r_ref) is indistinguishable
# from L_B at our resolution. β_phys=1 n=6 hits the matrix-element floor
# (≈ 3e-16) already at r=7 [see scripts/output/quad_redo_S1_krylov/
# n6_bphys1.0_smooth_def_s.bson], so R_REF = 10 buys ~3 extra bits of
# headroom over the β_phys=1 saturation point — generous, and still fast.
#
# Avoiding the dense L_b lets us run n=6 in ~15s of per-cell sweep wall
# (vs 113s on the dense-reference path, where 98s was the one-time L_b
# build). All matvecs go through `apply_lindbladian!` and benefit from the
# OMEGA_THREAD_THRESHOLD ω-loop threading.
#
# Run with:
#   JULIA_NUM_THREADS=8 julia --project scripts/scratch_quad_S1_bphys2_krylov.jl
#
# Env overrides:
#   N_LIST        comma-separated n list   (default 3,4,5,6)
#   BETA_PHYS     β_phys                   (default 2.0)
#   R_MIN, R_MAX  test r_D sweep range     (default 3..11)
#   R_REF         reference r_D            (default 13)
#   FILTER        smooth_def_s|smooth_fixed_s|gaussian|kinky  (default smooth_def_s)

using Printf
using LinearAlgebra
using Random
using BSON
using QuantumFurnace
using QuantumFurnace: apply_lindbladian!, apply_adjoint_lindbladian!, default_smooth_s

# ── Campaign grid ───────────────────────────────────────────────────────────
const N_LIST    = parse.(Int, split(get(ENV, "N_LIST",   "3,4,5,6"), ','))
const BETA_PHYS = parse(Float64, get(ENV, "BETA_PHYS", "2.0"))
const FILTER    = Symbol(get(ENV, "FILTER", "smooth_def_s"))
const R_MIN     = parse(Int, get(ENV, "R_MIN", "3"))
const R_MAX     = parse(Int, get(ENV, "R_MAX", "11"))
const R_REF     = parse(Int, get(ENV, "R_REF", "13"))
const TAIL_C    = parse(Float64, get(ENV, "TAIL_C", "8.0"))
const ETA       = parse(Float64, get(ENV, "ETA", "1e-3"))
const CELL_TIME_BUDGET = parse(Float64, get(ENV, "CELL_TIME_BUDGET", "2400.0"))
const MAX_PI_ITER = parse(Int,     get(ENV, "MAX_PI_ITER", "60"))
const REL_TOL     = parse(Float64, get(ENV, "REL_TOL",     "1e-6"))
const MIN_PI_ITER = parse(Int,     get(ENV, "MIN_PI_ITER", "6"))
const ABS_FLOOR   = parse(Float64, get(ENV, "ABS_FLOOR",   "1e-14"))
const OUTDIR      = get(ENV, "OUTDIR",
                        joinpath(@__DIR__, "output", "quad_redo_S1_bphys2_krylov"))

Random.seed!(20260512)
mkpath(OUTDIR)

# Matrix-free path: keep BLAS single-threaded so it doesn't fight the Julia
# ω-loop threads inside `apply_lindbladian!`. The ω-loop threshold is
# OMEGA_THREAD_THRESHOLD = 10 (src/krylov_matvec.jl:122), so as long as
# 2^r ≥ 10 the per-jump ω accumulation is parallelised across Julia threads.
BLAS.set_num_threads(1)

@info "Campaign grid" N_LIST BETA_PHYS FILTER R_MIN R_MAX R_REF TAIL_C ETA CELL_TIME_BUDGET MAX_PI_ITER REL_TOL ABS_FLOOR Threads.nthreads() BLAS.get_num_threads()
flush(stderr); flush(stdout)

# ── Per-cell EnergyDomain config builder (per-register API, qf-9z0) ─────────
function build_cfg_e(n_qubits, beta_alg, sigma, omega_range, filter::Symbol,
                     r_D::Int, beta_phys::Float64)
    w0_D = omega_range / 2^r_D
    t0_D = 2π / (2^r_D * w0_D)
    common = (
        sim = Lindbladian(), domain = EnergyDomain(), construction = KMS(),
        num_qubits = n_qubits, beta = beta_alg, beta_phys = beta_phys, sigma = sigma,
        num_energy_bits_D = r_D, w0_D = w0_D, t0_D = t0_D,
    )
    if filter === :gaussian
        sigma_gamma = sigma
        w_gamma = beta_alg * (sigma^2 + sigma_gamma^2) / 2
        return Config(; common...,
            with_linear_combination = false,
            gaussian_parameters = (w_gamma, sigma_gamma),
        )
    elseif filter === :kinky
        return Config(; common..., with_linear_combination = true,
            a = 0.0, s = 0.0, eta = ETA)
    elseif filter === :smooth_def_s
        return Config(; common..., with_linear_combination = true,
            a = 0.0, s = default_smooth_s(beta_alg, sigma), eta = ETA)
    elseif filter === :smooth_fixed_s
        return Config(; common..., with_linear_combination = true,
            a = 0.0, s = 0.25, eta = ETA)
    else
        error("Unknown filter: $filter")
    end
end

function jumps_for(ham, n_qubits)
    jp = [[X], [Y], [Z]]
    norm = sqrt(length(jp) * n_qubits)
    out = JumpOp[]
    for pauli in jp, site in 1:n_qubits
        op = Matrix(pad_term(pauli, n_qubits, site)) ./ norm
        op_eb = ham.eigvecs' * op * ham.eigvecs
        push!(out, JumpOp(op, op_eb, op == transpose(op), op == op'))
    end
    return out
end

# ── Two-EnergyDomain matrix-free op-norm (no dense reference) ───────────────
# Custom power iteration on A^† A with A = L_E(r_test) − L_E(r_ref). We read
# the norm as ‖A x_final‖ rather than √λ_max(A^†A) to preserve absolute
# matrix-element precision through the cancellation between two close
# Lindbladians. Same scheme as scratch_quad_S1_krylov_campaign.jl, but with
# the dense BohrDomain matvec replaced by a second EnergyDomain matvec.
function op_norm_matrix_free_2e(ws_test, cfg_test, ws_ref, cfg_ref, ham, d::Int;
                                max_iter::Int = 60, rel_tol::Real = 1e-6,
                                abs_floor::Real = 1e-14, min_iter::Int = 6)
    function fwd!(out_buf, rho)
        L_t = apply_lindbladian!(ws_test, rho, cfg_test, ham; include_coherent = false)
        L_r = apply_lindbladian!(ws_ref,  rho, cfg_ref,  ham; include_coherent = false)
        @inbounds for i in eachindex(out_buf)
            out_buf[i] = L_t[i] - L_r[i]
        end
        return out_buf
    end
    function adj!(out_buf, rho)
        L_t = apply_adjoint_lindbladian!(ws_test, rho, cfg_test, ham; include_coherent = false)
        L_r = apply_adjoint_lindbladian!(ws_ref,  rho, cfg_ref,  ham; include_coherent = false)
        @inbounds for i in eachindex(out_buf)
            out_buf[i] = L_t[i] - L_r[i]
        end
        return out_buf
    end

    x = randn(ComplexF64, d, d); x ./= norm(x)
    y = similar(x); z = similar(x)
    sigma_prev = 0.0
    converged_iter = 0
    for k in 1:max_iter
        fwd!(y, x)
        sigma = norm(y)
        if sigma < abs_floor
            converged_iter = -k
            break
        end
        adj!(z, y)
        nz = norm(z)
        nz == 0.0 && break
        x .= z ./ nz
        if k >= min_iter && abs(sigma - sigma_prev) <= max(rel_tol * sigma, abs_floor)
            converged_iter = k
            break
        end
        sigma_prev = sigma
    end
    fwd!(y, x)
    return norm(y), converged_iter
end

# ── Per-cell runner ─────────────────────────────────────────────────────────
function run_cell(n_qubits::Int, beta_phys::Float64, filter::Symbol)
    ham_path = joinpath(@__DIR__, "..", "hamiltonians",
                       "heis_xxx_zzdisordered_periodic_n$(n_qubits).bson")
    ham_raw  = BSON.load(ham_path)[:hamiltonian]
    ham      = HamHam(ham_raw; beta_phys = beta_phys)
    rescale  = ham.rescaling_factor
    beta_alg = beta_phys * rescale
    sigma    = 1.0 / beta_alg
    H_norm   = opnorm(ham.data)
    omega_range = 2 * (H_norm + TAIL_C * sigma)
    s_used   = filter === :smooth_def_s  ? default_smooth_s(beta_alg, sigma) :
              filter === :smooth_fixed_s ? 0.25 : NaN
    jumps    = jumps_for(ham, n_qubits)
    d        = size(ham.data, 1)

    @info "Cell" n_qubits beta_phys beta_alg sigma s_used omega_range filter R_REF R_MIN R_MAX d
    flush(stdout); flush(stderr)

    cfg_ref = build_cfg_e(n_qubits, beta_alg, sigma, omega_range, filter,
                          R_REF, beta_phys)
    t_build_ref = @elapsed ws_ref = Workspace(cfg_ref, ham, jumps)

    rows = NamedTuple[]
    cell_wall = t_build_ref
    @printf "  build ref ws (EnergyDomain, r=%d) wall = %.2fs   d=%d   2^r·n_jumps = %d\n" R_REF t_build_ref d (2^R_REF * length(jumps))
    @printf "  %-4s %12s %14s %10s %10s %8s\n" "r_D" "w0_D" "‖ΔL‖_op" "build_ws_s" "krylov_s" "conv"
    flush(stdout)
    for r in R_MIN:R_MAX
        cell_wall > CELL_TIME_BUDGET && (@warn "  ⏰ cell budget exceeded at r=$r"; break)
        cfg_test = build_cfg_e(n_qubits, beta_alg, sigma, omega_range, filter,
                                r, beta_phys)
        w0_D = omega_range / 2^r
        t_build_ws_test = @elapsed ws_test = Workspace(cfg_test, ham, jumps)
        t_norm = @elapsed (op_norm, converged) = op_norm_matrix_free_2e(
            ws_test, cfg_test, ws_ref, cfg_ref, ham, d;
            max_iter = MAX_PI_ITER, rel_tol = REL_TOL,
            min_iter = MIN_PI_ITER, abs_floor = ABS_FLOOR,
        )
        push!(rows, (; r, w0_D, err = op_norm,
                       t_build_ws_test, t_norm, converged))
        @printf "  %-4d %12.6e %14.6e %10.3f %10.3f %8d\n" r w0_D op_norm t_build_ws_test t_norm converged
        flush(stdout)
        cell_wall += t_build_ws_test + t_norm
        ws_test = nothing
        GC.gc()
        if op_norm < 1e-13 && r >= R_MIN + 4
            @info "  → matrix-element floor at r=$r; stopping sweep early."
            break
        end
    end

    cutoffs = Dict(
        :eps_1e_3 => let i = findfirst(row -> row.err <= 1e-3, rows); i === nothing ? missing : rows[i].r end,
        :eps_1e_6 => let i = findfirst(row -> row.err <= 1e-6, rows); i === nothing ? missing : rows[i].r end,
        :eps_1e_9 => let i = findfirst(row -> row.err <= 1e-9, rows); i === nothing ? missing : rows[i].r end,
    )
    println("  cutoffs: ", cutoffs)

    sidecar_path = joinpath(OUTDIR,
                            "n$(n_qubits)_bphys$(beta_phys)_$(filter).bson")
    BSON.bson(sidecar_path, Dict(
        :cell => (
            n_qubits = n_qubits, beta_phys = beta_phys, beta_alg = beta_alg,
            rescaling_factor = rescale, sigma = sigma, h_norm = H_norm,
            omega_range = omega_range, filter = filter, s_used = s_used,
            tail_c = TAIL_C, eta = ETA, d = d,
            r_ref = R_REF,
            max_iter = MAX_PI_ITER, rel_tol = REL_TOL,
            min_iter = MIN_PI_ITER, abs_floor = ABS_FLOOR,
            julia_threads = Threads.nthreads(),
        ),
        :sweep => rows,
        :t_build_ref => t_build_ref,
        :total_wall => cell_wall,
        :cutoffs => cutoffs,
    ))
    @printf "  → saved %s   (cell wall %.1fs)\n\n" sidecar_path cell_wall

    ws_ref = nothing
    GC.gc()
    return (; n_qubits, beta_phys, filter, cell_wall, cutoffs, sidecar_path)
end

# ── Campaign loop ───────────────────────────────────────────────────────────
campaign_summary = NamedTuple[]
for n_qubits in N_LIST
    try
        rec = run_cell(n_qubits, BETA_PHYS, FILTER)
        push!(campaign_summary, rec)
    catch err
        @error "  Cell failed" n_qubits BETA_PHYS FILTER exception=(err, catch_backtrace())
    end
end

println("\n=== Campaign summary ===")
@printf "%-3s %-6s %-18s %10s %8s %8s %8s\n" "n" "β_phys" "filter" "wall_s" "r@1e-3" "r@1e-6" "r@1e-9"
for r in campaign_summary
    @printf "%-3d %-6.2f %-18s %10.1f %8s %8s %8s\n" r.n_qubits r.beta_phys string(r.filter) r.cell_wall (
        ismissing(r.cutoffs[:eps_1e_3]) ? "—" : string(r.cutoffs[:eps_1e_3])
    ) (
        ismissing(r.cutoffs[:eps_1e_6]) ? "—" : string(r.cutoffs[:eps_1e_6])
    ) (
        ismissing(r.cutoffs[:eps_1e_9]) ? "—" : string(r.cutoffs[:eps_1e_9])
    )
end
summary_path = joinpath(OUTDIR, "_campaign_summary.bson")
BSON.bson(summary_path, Dict(:summary => campaign_summary))
@printf "→ saved campaign summary to %s\n" summary_path
