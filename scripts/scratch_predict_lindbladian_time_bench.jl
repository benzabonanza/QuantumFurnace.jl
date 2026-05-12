#!/usr/bin/env julia
#
# qf-1au — speed benchmark for predict_lindbladian_trajectory on TimeDomain
# with the FULL Lindbladian (include_coherent ON by default), matrix-free,
# multithreaded.
#
# Background
# ----------
# `predict_lindbladian_trajectory` (src/lindblad_action.jl:601) is currently
# dispatched only on `Config{Lindbladian, <:Union{BohrDomain, EnergyDomain}}`.
# Its body, however, is domain-agnostic: it builds an Arnoldi factorisation
# via repeated `apply_lindbladian!` matvecs through a closure. To benchmark
# the TimeDomain path without touching src/ in this scratch pass, we call
# the same inner helper `_krylov_spectral_decomposition(fwd!, rho_0, d;
# krylovdim, tol)` directly with a TimeDomain workspace. If TimeDomain is
# competitive here, the right follow-up is to relax the dispatch guard in
# src/lindblad_action.jl to `Union{BohrDomain, EnergyDomain, TimeDomain,
# TrotterDomain}`.
#
# Registers (canonical v2 recipe, smooth_def_s)
# ---------------------------------------------
# r_D       = 7      (dissipator, qf-w2u — 10⁻⁹ at β_phys ≤ 2)
# r_b_minus = 7      (coherent outer leg, v2 recipe — 10⁻⁶ floor)
# r_b_plus  = 17     (coherent inner leg, v2 recipe — slope-(−1) cap)
#
# At β_phys = 2 the coherent term may need r_b_plus = 18 to keep its 10⁻⁶
# envelope; per-cell ε is not the goal here, speed is. Holding r_b_plus =
# 17 gives a fair speed comparison across (n, β_phys).
#
# Per cell we report:
#   wall_workspace   — Workspace ctor wall (includes B_time 2D NUFFT,
#                      OFT NUFFT prefactor cache)
#   wall_krylov      — _krylov_spectral_decomposition wall
#   total_wall       — sum
#   matvecs          — number of `apply_lindbladian!` calls
#   ms_per_matvec    — wall_krylov / matvecs * 1000
#
# Run with:
#   JULIA_NUM_THREADS=8 julia --project scripts/scratch_predict_lindbladian_time_bench.jl

using Printf
using LinearAlgebra
using Random
using BSON
using QuantumFurnace
using QuantumFurnace: apply_lindbladian!, _krylov_spectral_decomposition,
                      default_smooth_s

# ── Campaign grid ───────────────────────────────────────────────────────────
const N_LIST          = parse.(Int,     split(get(ENV, "N_LIST",   "3,4,5,6"), ','))
const BETA_PHYS_LIST  = parse.(Float64, split(get(ENV, "BETA_PHYS_LIST", "1.0,2.0"), ','))
const FILTER          = Symbol(get(ENV, "FILTER", "smooth_def_s"))
const R_D             = parse(Int, get(ENV, "R_D",        "7"))
const R_BMINUS        = parse(Int, get(ENV, "R_BMINUS",   "7"))
const R_BPLUS         = parse(Int, get(ENV, "R_BPLUS",    "17"))
const TAIL_C          = parse(Float64, get(ENV, "TAIL_C", "8.0"))
const ETA             = parse(Float64, get(ENV, "ETA",    "1e-3"))
const KRYLOVDIM       = parse(Int, get(ENV, "KRYLOVDIM",  "30"))
const TOL             = parse(Float64, get(ENV, "TOL",    "1e-10"))
const T_RANGE_MINUS   = parse(Float64, get(ENV, "T_RANGE_MINUS", "18.0"))
const T_RANGE_PLUS    = parse(Float64, get(ENV, "T_RANGE_PLUS",  "12.0"))
const OUTDIR          = get(ENV, "OUTDIR",
                            joinpath(@__DIR__, "output", "predict_lindbladian_time_bench"))

Random.seed!(20260512)
mkpath(OUTDIR)
BLAS.set_num_threads(1)

@info "Bench grid" N_LIST BETA_PHYS_LIST FILTER R_D R_BMINUS R_BPLUS KRYLOVDIM TOL TAIL_C T_RANGE_MINUS T_RANGE_PLUS Threads.nthreads() BLAS.get_num_threads()
flush(stderr); flush(stdout)

# ── Config builder with full per-register triples ───────────────────────────
function build_cfg(n_qubits, beta_alg, sigma, omega_range,
                   filter::Symbol, r_D::Int, r_bm::Int, r_bp::Int,
                   beta_phys::Float64)
    # Dissipator register
    w0_D = omega_range / 2^r_D
    t0_D = 2π / (2^r_D * w0_D)
    # Coherent registers — `w0_X = T_RANGE_X / 2^r_X · w0_D / t0_X` is not
    # the structure; the per-leg triple (`r_X`, `t0_X`, `w0_X`) just has to
    # satisfy `t0_X = 2π / (2^r_X · w0_X)`. Set `t0_X = T_RANGE_X / 2^r_X`
    # (so the coherent grid covers [-T/2, T/2]) and derive `w0_X`.
    t0_bm = T_RANGE_MINUS / 2^r_bm
    w0_bm = 2π / (2^r_bm * t0_bm)
    t0_bp = T_RANGE_PLUS / 2^r_bp
    w0_bp = 2π / (2^r_bp * t0_bp)
    common = (
        sim = Lindbladian(), domain = TimeDomain(), construction = KMS(),
        num_qubits = n_qubits, beta = beta_alg, beta_phys = beta_phys, sigma = sigma,
        num_energy_bits_D = r_D, w0_D = w0_D, t0_D = t0_D,
        num_energy_bits_b_minus = r_bm, w0_b_minus = w0_bm, t0_b_minus = t0_bm,
        num_energy_bits_b_plus  = r_bp, w0_b_plus  = w0_bp, t0_b_plus  = t0_bp,
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
    nrm = sqrt(length(jp) * n_qubits)
    out = JumpOp[]
    for pauli in jp, site in 1:n_qubits
        op = Matrix(pad_term(pauli, n_qubits, site)) ./ nrm
        op_eb = ham.eigvecs' * op * ham.eigvecs
        push!(out, JumpOp(op, op_eb, op == transpose(op), op == op'))
    end
    return out
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

    @info "Cell" n_qubits beta_phys beta_alg sigma s_used omega_range filter R_D R_BMINUS R_BPLUS d
    flush(stdout); flush(stderr)

    cfg = build_cfg(n_qubits, beta_alg, sigma, omega_range, filter,
                    R_D, R_BMINUS, R_BPLUS, beta_phys)

    # Workspace builds B_time (full coherent term) + OFT NUFFT prefactors.
    t_ws = @elapsed ws = Workspace(cfg, ham, jumps)

    # Closure: full Lindbladian matvec, include_coherent default = true.
    fwd! = let ws = ws, cfg = cfg, ham = ham
        (out::AbstractMatrix, x::AbstractMatrix) -> begin
            apply_lindbladian!(ws, x, cfg, ham)
            copyto!(out, ws.scratch.rho_out)
            return out
        end
    end

    # Initial state: Gibbs of an unrelated Hamiltonian (here just I/d for
    # generality; rho_0 has nonzero overlap with all slow modes).
    rho_0 = Matrix{ComplexF64}(I, d, d) ./ d

    # Spectral decomposition — the timed core of predict_lindbladian_trajectory.
    GC.gc()
    rss_before = Sys.maxrss() / (1024^2)
    t_kry = @elapsed decomp = _krylov_spectral_decomposition(
        fwd!, rho_0, d;
        krylovdim = KRYLOVDIM, tol = TOL,
        sort_mode = :lindbladian,
    )
    rss_after = Sys.maxrss() / (1024^2)

    matvecs = decomp.matvec_count
    gap_estimate = length(decomp.eigenvalues) >= 2 ? abs(real(decomp.eigenvalues[2])) : NaN
    total_wall = t_ws + t_kry

    @printf "  d=%-4d  ws=%6.2fs  krylov=%6.2fs  total=%6.2fs  matvecs=%3d  ms/mv=%6.1f  gap≈%.3e  RSS=%.0fMB\n" d t_ws t_kry total_wall matvecs (t_kry / matvecs * 1000) gap_estimate rss_after
    flush(stdout)

    sidecar_path = joinpath(OUTDIR,
                            "n$(n_qubits)_bphys$(beta_phys)_$(filter)_time.bson")
    BSON.bson(sidecar_path, Dict(
        :cell => (
            n_qubits = n_qubits, beta_phys = beta_phys, beta_alg = beta_alg,
            rescaling_factor = rescale, sigma = sigma, h_norm = H_norm,
            omega_range = omega_range, filter = filter, s_used = s_used,
            tail_c = TAIL_C, eta = ETA, d = d,
            r_D = R_D, r_bminus = R_BMINUS, r_bplus = R_BPLUS,
            t_range_minus = T_RANGE_MINUS, t_range_plus = T_RANGE_PLUS,
            krylovdim = KRYLOVDIM, tol = TOL,
            include_coherent = true,
            domain = :TimeDomain,
            julia_threads = Threads.nthreads(),
        ),
        :timing => (
            t_workspace = t_ws,
            t_krylov    = t_kry,
            total_wall  = total_wall,
            matvecs     = matvecs,
            ms_per_matvec = t_kry / matvecs * 1000,
            rss_before_mb = rss_before,
            rss_after_mb  = rss_after,
        ),
        :spectrum => (
            gap_estimate = gap_estimate,
            n_eigenmodes = length(decomp.eigenvalues),
            converged    = decomp.converged,
        ),
    ))

    ws = nothing; decomp = nothing
    GC.gc()
    return (; n_qubits, beta_phys, d,
             t_ws, t_kry, total_wall, matvecs,
             ms_per_matvec = t_kry / matvecs * 1000,
             gap = gap_estimate, rss_mb = rss_after,
             sidecar_path)
end

# ── Campaign loop ───────────────────────────────────────────────────────────
println("\n══════════════════════════════════════════════════════════════════════════")
println("Bench: predict_lindbladian_trajectory @ TimeDomain, coherent ON, krylovdim=$KRYLOVDIM")
println("Registers: r_D=$R_D  r_b_minus=$R_BMINUS  r_b_plus=$R_BPLUS")
println("══════════════════════════════════════════════════════════════════════════\n")

campaign_summary = NamedTuple[]
for n_qubits in N_LIST
    for beta_phys in BETA_PHYS_LIST
        try
            rec = run_cell(n_qubits, beta_phys, FILTER)
            push!(campaign_summary, rec)
        catch err
            @error "  Cell failed" n_qubits beta_phys FILTER exception=(err, catch_backtrace())
        end
    end
end

println("\n══════════════════════════════════════════════════════════════════════════")
println("Bench summary")
println("══════════════════════════════════════════════════════════════════════════")
@printf "%-3s %-7s %-5s %10s %10s %10s %8s %10s %8s\n" "n" "β_phys" "d" "ws_s" "krylov_s" "total_s" "matvecs" "ms/mv" "RSS_MB"
for r in campaign_summary
    @printf "%-3d %-7.2f %-5d %10.2f %10.2f %10.2f %8d %10.1f %8.0f\n" r.n_qubits r.beta_phys r.d r.t_ws r.t_kry r.total_wall r.matvecs r.ms_per_matvec r.rss_mb
end
summary_path = joinpath(OUTDIR, "_campaign_summary.bson")
BSON.bson(summary_path, Dict(:summary => campaign_summary))
@printf "\n→ saved campaign summary to %s\n" summary_path
