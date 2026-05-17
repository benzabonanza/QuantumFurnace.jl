#!/usr/bin/env julia
# scratch_qf_e4z_33_kdim_sweep.jl  (qf-e4z.33)
#
# Pin gap accuracy: dense reference at n ≤ 6 + dual krylovdim sweeps for
# Pass 1 / Pass 2 at every cell of the canonical v6_plus grid, on a SINGLE
# seed (42) per n. Goal: see whether the n=7 low-β Pass1↔Pass2 disagreement
# of 1e-4 in qf-e4z.31 is a Krylov-truncation residual (saturates with kdim
# bump), an algorithm bias (cross-disagreement persists), or genuinely
# Krylov-resistant (neither pass self-saturates).
#
# Methodology matches scratch_p1_v6_multiseed_plus.jl, with two differences:
#   1. Single seed (42) per n instead of 5 seeds — user's call for a kdim
#      diagnostic, not a statistics-of-disorder run.
#   2. Per cell:
#        • Pass 1 sweep:  kdim_p1 ∈ {40, 60, 80, 100}     (one
#          `_krylov_spectral_decomposition` call per kdim, workspace reused)
#        • Pass 2 sweep:  kdim_p2 ∈ {30, 60, 80, 100}     (one
#          `krylov_spectral_gap` call per kdim, workspace reused; uses the
#          qf-e4z.33 `krylovdim_gap_pass` kwarg semantics)
#        • Dense ref:     n ≤ 6 only — `construct_lindbladian + eigvals`,
#          smallest |Re(λ)| above 1e-10 threshold.
#
# Fixtures: heis_xxx_disordered_periodic_n{3..7}_seed42.bson (same as v6_plus).
# r_D map: β_phys ≤ 1.0 → r_D=7, β_phys ≥ 1.5 → r_D=8 (same as v6_plus).
# Smooth Metro filter: s = 0.25, a = 0 (qf-yt9 canonical, same as v6_plus).
# Initial state for Pass 1: rho_0 = |+⟩⟨+|^⊗N (qf-e4z.30 canonical).
# Coherent term: ON (script-default `include_coherent=true`, the physical L).
#
# Cells: 5 n × 6 β_phys × 1 seed = 30
# Pass 1 calls per cell: 4 → 120 total
# Pass 2 calls per cell: 4 → 120 total
# Dense ref calls: 4 n × 6 β_phys = 24 (n=3..6 only)
#
# Sidecar dir: scripts/output/sweep_qf_e4z_33_kdim_sweep/seed42/
# Filename:    kdim_n{n}_betaphys{β_phys}_seed42_L_KMS_Energy.bson
#
# Usage (run from repo root):
#   JULIA_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 \
#     julia --project scripts/scratch_qf_e4z_33_kdim_sweep.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using QuantumFurnace
using QuantumFurnace: _parse_hamiltonian_bson, _build_jump_set,
                     _krylov_spectral_decomposition,
                     Workspace, apply_lindbladian!,
                     construct_lindbladian
using LinearAlgebra
using BSON
using Printf
using Dates

BLAS.set_num_threads(1)
@assert BLAS.get_num_threads() == 1
println("[init] Julia threads = ", Threads.nthreads(), ", BLAS threads = ", BLAS.get_num_threads())
println("[init] hostname = ", gethostname(), "   ", now())

# --- Constants -------------------------------------------------------------

const OUTPUT_ROOT = joinpath(@__DIR__, "output")
const FIXTURE_DIR = joinpath(OUTPUT_ROOT, "multiseed_fixtures")

const BETA_PHYS_ALL = (0.25, 0.5, 1.0, 1.5, 2.0, 2.5)
const N_MIN, N_MAX  = 3, 7
const SEED          = 42
const TAIL_C        = 8.0

const KDIM_P1 = (40, 60, 80, 100)
const KDIM_P2 = (30, 60, 80, 100)

# Dense gap is feasible up to d²=4096 (n=6) at ~256 MB. n=7 = 4 GB and
# infeasible in the sandbox.
const DENSE_MAX_N = 6
const ZERO_MODE_TOL = 1e-10            # |λ| below this counted as the steady-state zero

_r_D_for_beta_phys(β_phys::Real) = β_phys ≥ 1.5 ? 8 : 7

const OUT_DIR = joinpath(OUTPUT_ROOT, "sweep_qf_e4z_33_kdim_sweep", "seed42")
mkpath(OUT_DIR)

# --- |+⟩⟨+|^⊗N density matrix (parity-broken seed) -----------------------

function rho_plus_tensor(n::Integer)
    plus = ComplexF64[0.5 0.5; 0.5 0.5]
    rho = plus
    for _ in 2:n
        rho = kron(rho, plus)
    end
    return rho
end

# --- Sidecar / fixture path helpers --------------------------------------

function _sidecar_path(n::Integer, beta_phys::Real)
    bp_str = let s = @sprintf("%.6f", float(beta_phys))
        s = rstrip(s, '0'); s = rstrip(s, '.')
        isempty(s) ? "0" : s
    end
    return joinpath(OUT_DIR,
        "kdim_n$(n)_betaphys$(bp_str)_seed$(SEED)_L_KMS_Energy.bson")
end

function _load_ham(n::Integer, beta_phys::Real)
    path = joinpath(FIXTURE_DIR,
        "heis_xxx_disordered_periodic_n$(n)_seed$(SEED).bson")
    isfile(path) || error("Missing fixture: $path")
    raw = _parse_hamiltonian_bson(path)
    return HamHam(raw; beta_phys = float(beta_phys))
end

# --- Per-cell config builder (identical to v6_plus) -----------------------

function _build_ckg_cfg(n::Integer, beta_phys::Real, ham; r_D::Integer)
    β_alg = beta_alg(ham, float(beta_phys))
    σ     = 1.0 / β_alg
    H_norm = maximum(abs, ham.eigvals)
    omega_range = 2.0 * (H_norm + TAIL_C * σ)
    w0_D = omega_range / 2.0^r_D
    t0_D = 2π / (2.0^r_D * w0_D)
    return Config(
        sim = Lindbladian(),
        domain = EnergyDomain(),
        construction = KMS(),
        num_qubits = n,
        with_linear_combination = true,
        beta = β_alg,
        beta_phys = float(beta_phys),
        sigma = σ,
        a = 0.0, s = 0.25,
        gaussian_parameters = (nothing, nothing),
        num_energy_bits_D = r_D, w0_D = w0_D, t0_D = t0_D,
        num_trotter_steps_per_t0 = 10,
        filter = nothing,
    )
end

# --- Dense L gap (n ≤ 6) --------------------------------------------------

"""
    dense_gap(cfg, ham, jumps; zero_tol=1e-10) -> (gap, top5_re)

Build the d²×d² dense Lindbladian via `construct_lindbladian`, compute its
full eigvals, and return the smallest |Re(λ)| above `zero_tol` (the
nonzero leading mode = canonical "spectral gap"). Coherent term is ON (the
default of `construct_lindbladian`). Also returns the 5 smallest |Re(λ)|
values for diagnostic (e.g. parity sub-spectrum structure).
"""
function dense_gap(cfg, ham, jumps; zero_tol::Float64 = ZERO_MODE_TOL)
    L = construct_lindbladian(jumps, cfg, ham)          # dense d²×d² complex
    λs = eigvals(L)
    re_abs = abs.(real.(λs))
    # Sort ascending; smallest is zero (trace preservation); take first
    # nonzero entry above tolerance.
    perm = sortperm(re_abs)
    sorted = re_abs[perm]
    nonzero = filter(x -> x > zero_tol, sorted)
    gap = isempty(nonzero) ? NaN : nonzero[1]
    top5 = first(nonzero, min(5, length(nonzero)))
    return gap, top5
end

# --- Pass 1: single-seed Arnoldi from vec(|+⟩⟨+|^⊗N) at chosen kdim ------

function pass1_gap(ws, cfg, ham, rho_0::Matrix{ComplexF64}, kdim::Int)
    d = size(rho_0, 1)
    fwd! = let ws = ws, cfg = cfg, ham = ham
        (out::AbstractMatrix, x::AbstractMatrix) -> begin
            apply_lindbladian!(ws, x, cfg, ham)
            copyto!(out, ws.scratch.rho_out)
            return out
        end
    end
    decomp = _krylov_spectral_decomposition(fwd!, rho_0, d;
                                            krylovdim=kdim, tol=1e-10)
    gap = abs(real(decomp.eigenvalues[2]))
    return gap, decomp.matvec_count, decomp.converged
end

# --- Pass 2: krylov_spectral_gap from _krylov_default_x0 + thick restart -

function pass2_gap(cfg, ham, jumps, kdim::Int)
    res = krylov_spectral_gap(cfg, ham, jumps;
                              krylovdim = kdim, howmany = 4, tol = 1e-10)
    return res.spectral_gap, res.matvec_count, Int(res.converged)
end

# --- Per-cell runner ------------------------------------------------------

function run_cell(; n, beta_phys, skip_existing::Bool = true)
    sidecar = _sidecar_path(n, beta_phys)

    if skip_existing && isfile(sidecar)
        try
            d = BSON.load(sidecar, QuantumFurnace)
            r = d[:result]
            @printf("[kdim] n=%d β_phys=%-4g | SKIP (cached) dense=%.6e p1@100=%.6e p2@100=%.6e\n",
                    n, beta_phys, r[:gap_dense],
                    r[:gap_p1_grid][end], r[:gap_p2_grid][end])
            return (; (Symbol(k) => v for (k, v) in pairs(r))...)
        catch err
            @warn "Cached sidecar failed to load; re-running" sidecar err
        end
    end

    t0_run = time()
    ham = _load_ham(n, beta_phys)
    jumps = _build_jump_set(ham, n)
    r_D = _r_D_for_beta_phys(beta_phys)
    cfg = _build_ckg_cfg(n, beta_phys, ham; r_D = r_D)
    rho_0 = rho_plus_tensor(n)

    # Build workspace once, reuse across all Pass 1 kdims (Pass 2 builds
    # its own inside krylov_spectral_gap — that's fine, the matvec cost
    # dominates the workspace cost in the Pass 2 path).
    ws = Workspace(cfg, ham, jumps)

    # --- Dense ref (n ≤ 6 only) ---
    gap_dense = NaN
    dense_top5 = Float64[]
    dense_wall = 0.0
    if n ≤ DENSE_MAX_N
        t_dense = time()
        gap_dense, dense_top5 = dense_gap(cfg, ham, jumps)
        dense_wall = time() - t_dense
    end

    # --- Pass 1 sweep ---
    kdim_p1_grid = collect(KDIM_P1)
    gap_p1 = Float64[]
    mv_p1  = Int[]
    conv_p1 = Bool[]
    p1_wall_each = Float64[]
    for kdim in kdim_p1_grid
        t_p1 = time()
        gap, mv, conv = pass1_gap(ws, cfg, ham, rho_0, kdim)
        push!(gap_p1, gap); push!(mv_p1, mv); push!(conv_p1, conv)
        push!(p1_wall_each, time() - t_p1)
    end

    # --- Pass 2 sweep ---
    kdim_p2_grid = collect(KDIM_P2)
    gap_p2 = Float64[]
    mv_p2  = Int[]
    conv_p2 = Int[]
    p2_wall_each = Float64[]
    for kdim in kdim_p2_grid
        t_p2 = time()
        gap, mv, conv = pass2_gap(cfg, ham, jumps, kdim)
        push!(gap_p2, gap); push!(mv_p2, mv); push!(conv_p2, conv)
        push!(p2_wall_each, time() - t_p2)
    end

    wall = time() - t0_run

    result = (
        n = n, beta_phys = float(beta_phys), beta_alg = cfg.beta,
        rescaling_factor = ham.rescaling_factor,
        nu_min = ham.nu_min, disorder_strength = 0.1,
        seed = SEED,
        r_D = register_r_D(cfg), w0_D = register_w0_D(cfg), t0_D = register_t0_D(cfg),
        s = cfg.s, a = cfg.a, sigma = cfg.sigma,

        # Dense reference (NaN if n > DENSE_MAX_N).
        gap_dense = gap_dense,
        dense_top5 = dense_top5,
        dense_wall = dense_wall,

        # Pass 1 sweep.
        kdim_p1_grid = kdim_p1_grid,
        gap_p1_grid = gap_p1,
        matvecs_p1_grid = mv_p1,
        converged_p1_grid = conv_p1,
        wall_p1_each = p1_wall_each,

        # Pass 2 sweep.
        kdim_p2_grid = kdim_p2_grid,
        gap_p2_grid = gap_p2,
        matvecs_p2_grid = mv_p2,
        converged_p2_grid = conv_p2,
        wall_p2_each = p2_wall_each,

        wall_time = wall,
        sweep_version = :qf_e4z_33_kdim_single_seed,
    )

    try
        BSON.bson(sidecar, Dict(:result => Dict(pairs(result))))
    catch err
        @warn "Sidecar write failed (continuing)" sidecar err
    end

    # Human-readable summary line.
    p1_end = gap_p1[end]; p2_end = gap_p2[end]
    cross_diff = abs(p1_end - p2_end) / max(abs(p2_end), 1e-30)
    if !isnan(gap_dense)
        rel_p1 = abs(p1_end - gap_dense) / gap_dense
        rel_p2 = abs(p2_end - gap_dense) / gap_dense
        @printf("[kdim] n=%d β_phys=%-4g (β_alg=%.2f, r_D=%d) | dense=%.6e  p1@%d=%.6e (rel=%.2e)  p2@%d=%.6e (rel=%.2e)  cross=%.2e  wall=%.1fs\n",
                n, beta_phys, cfg.beta, r_D,
                gap_dense,
                kdim_p1_grid[end], p1_end, rel_p1,
                kdim_p2_grid[end], p2_end, rel_p2,
                cross_diff, wall)
    else
        @printf("[kdim] n=%d β_phys=%-4g (β_alg=%.2f, r_D=%d) | dense=  N/A   p1@%d=%.6e  p2@%d=%.6e  cross=%.2e  wall=%.1fs\n",
                n, beta_phys, cfg.beta, r_D,
                kdim_p1_grid[end], p1_end,
                kdim_p2_grid[end], p2_end,
                cross_diff, wall)
    end
    # Per-kdim line (Pass 1 self-saturation).
    @printf("    P1 ")
    for (i, kdim) in enumerate(kdim_p1_grid)
        @printf(" k=%d:%.5e", kdim, gap_p1[i])
    end
    println()
    @printf("    P2 ")
    for (i, kdim) in enumerate(kdim_p2_grid)
        @printf(" k=%d:%.5e", kdim, gap_p2[i])
    end
    println()
    flush(stdout)
    return result
end

# --- Main ------------------------------------------------------------------

function main()
    println("\n[main] start $(now())  seed=$SEED  β_phys=$(BETA_PHYS_ALL)  n=$N_MIN..$N_MAX")
    println("[main] kdim_p1 grid = $KDIM_P1  | kdim_p2 grid = $KDIM_P2")
    println("[main] r_D map: β_phys ≤ 1.0 → r_D=7, β_phys ≥ 1.5 → r_D=8")
    println("[main] init state (Pass 1): rho_0 = |+⟩⟨+|^⊗N (qf-e4z.30 canonical)")
    println("[main] dense ref: enabled for n ≤ $DENSE_MAX_N (skipped at n=7)")
    println("[main] sidecar dir: $OUT_DIR")
    println("[main] expected cells: $(N_MAX - N_MIN + 1) × $(length(BETA_PHYS_ALL)) = $((N_MAX - N_MIN + 1) * length(BETA_PHYS_ALL))")

    rs = NamedTuple[]
    for n in N_MIN:N_MAX
        println("\n" * "="^72)
        println("=== n=$n  ($(length(BETA_PHYS_ALL)) β_phys cells)")
        println("="^72)
        for β_phys in BETA_PHYS_ALL
            try
                r = run_cell(; n = n, beta_phys = float(β_phys))
                push!(rs, r)
            catch err
                @error "[kdim n=$n β_phys=$β_phys] CELL CRASHED" exception=(err, catch_backtrace())
            end
            GC.gc(true)
        end
    end

    # Summary BSON.
    summary_path = joinpath(OUTPUT_ROOT, "sweep_qf_e4z_33_kdim_sweep_summary.bson")
    rows = Dict[]
    for r in rs
        push!(rows, Dict(pairs(r)...))
    end
    BSON.bson(summary_path, Dict(:rows => rows))
    @printf("\n[summary] wrote %d rows to %s\n", length(rows), summary_path)

    dir_size = sum(filesize(p) for p in readdir(OUT_DIR; join=true); init=0)
    @printf("[storage] sidecar dir: %d files, %.1f KB\n", length(readdir(OUT_DIR)), dir_size/1024)

    println("\n[main] done $(now())")
end

main()
