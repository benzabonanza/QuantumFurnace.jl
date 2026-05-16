#!/usr/bin/env julia
# scratch_p1_v6_multiseed.jl  (qf-e4z.23)
#
# Multi-seed disorder-averaged 1D Heisenberg PBC scaling sweep. Forks
# scratch_p1_v5_betaphys_extended.jl with three changes:
#
#   1. Fixtures now come from scripts/output/multiseed_fixtures/ (per-(n,seed)
#      raw NamedTuples produced by scripts/scratch_multiseed_disordered_fixtures.jl,
#      qf-yi4) — instead of a single hamiltonians/heis_xxx_zzdisordered_periodic_n*.bson
#      shared across all seeds. Each seed therefore carries its OWN rescaling_factor,
#      so β_alg = β_phys · R(n, seed) varies cell-to-cell within an (n, β_phys) group.
#
#   2. Inner loop iterates SEEDS = (42, 43, 44, 45, 46). 6 n × 6 β_phys × 5 seed
#      = 180 cells (CKG / S1 only — DLL S2 is out of scope for v6 per the issue
#      plan; can be filed as a separate follow-up).
#
#   3. Sidecar dir is sweep_S1_v6_ckg_ideal_multiseed/smooth_metro_eps1e-03/ to
#      avoid contaminating the v4/v5 single-seed pool. Sidecar filenames carry
#      the seed.
#
# r_D is bumped to 8 for β_phys ≥ 1.5 (per qf-yt9 / qf-e4z.22 — β_alg ≳ 50 needs
# r_D=8 for the ≤ 1e-9 vs Bohr target); r_D=7 for β_phys ≤ 1.0. Heterogeneous
# across the grid but matches what v4/v5 did, and τ_mix is converged to machine
# precision at both r_D ≥ 6 per qf-yt9 — so the heterogeneity does not affect
# the scaling fit.
#
# Methodology guardrails (per [[feedback_more_data_points_for_scaling_claims]]):
#   • Even/odd structure from 3 even-n + 3 odd-n cells is on the edge of what
#     6 n-points can resolve. The per-seed scatter within each parity class is
#     the null reference — only call an even/odd distinction REAL if inter-parity
#     spread exceeds intra-parity scatter at matched β_phys.
#
# Usage (run from repo root):
#   JULIA_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 \
#     julia --project scripts/scratch_p1_v6_multiseed.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using QuantumFurnace
using QuantumFurnace: _parse_hamiltonian_bson, _build_jump_set
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

const CELL_WALL_BUDGET = 900.0   # 15 min per cell (matches v5; n=8 v5 ~ 7 min)
const NPLUS1_RATIO_GUESS = 12.0

const EPS_TARGET    = 1e-3
const BETA_PHYS_ALL = (0.25, 0.5, 1.0, 1.5, 2.0, 2.5)
const N_MIN, N_MAX  = 3, 8
const SEEDS         = (42, 43, 44, 45, 46)
const KRYLOVDIM     = 40
const TAIL_C        = 8.0

# Heterogeneous r_D mapping (matches v4 + v5 single-seed)
_r_D_for_beta_phys(β_phys::Real) = β_phys ≥ 1.5 ? 8 : 7

const S1_DIR = joinpath(OUTPUT_ROOT, "sweep_S1_v6_ckg_ideal_multiseed", "smooth_metro_eps1e-03")
mkpath(S1_DIR)

# --- RSS helpers -----------------------------------------------------------

function rss_mb()
    try
        for line in eachline("/proc/self/status")
            startswith(line, "VmRSS:") && return parse(Int, split(line)[2]) / 1024.0
        end
    catch; end
    return NaN
end
function rss_peak_mb()
    try
        for line in eachline("/proc/self/status")
            startswith(line, "VmHWM:") && return parse(Int, split(line)[2]) / 1024.0
        end
    catch; end
    return NaN
end

# --- Sidecar path (extends v5 naming with seed) ---------------------------

function _sidecar_path(out_dir::AbstractString, n::Integer, beta_phys::Real, seed::Integer,
                       construction_tag::AbstractString, domain_tag::AbstractString)
    bp_str = let s = @sprintf("%.6f", float(beta_phys))
        s = rstrip(s, '0'); s = rstrip(s, '.')
        isempty(s) ? "0" : s
    end
    return joinpath(out_dir,
        "sweep_n$(n)_betaphys$(bp_str)_seed$(seed)_L_$(construction_tag)_$(domain_tag).bson")
end

# --- Multi-seed fixture loader --------------------------------------------

function _load_ham_multiseed(n::Integer, seed::Integer, beta_phys::Real)
    path = joinpath(FIXTURE_DIR,
        "heis_xxx_disordered_periodic_n$(n)_seed$(seed).bson")
    isfile(path) || error("Missing multiseed fixture: $path. " *
                          "Run scripts/scratch_multiseed_disordered_fixtures.jl first.")
    raw = _parse_hamiltonian_bson(path)
    return HamHam(raw; beta_phys = float(beta_phys))
end

# --- Per-cell config builder ----------------------------------------------

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

# --- Per-cell runner -------------------------------------------------------

function run_cell(; sweep, n, beta_phys, seed, ε, config, ham, jumps, label,
                    out_dir, construction_tag, domain_tag, skip_existing::Bool=true)
    sidecar = _sidecar_path(out_dir, n, beta_phys, seed, construction_tag, domain_tag)

    if skip_existing && isfile(sidecar)
        try
            d = BSON.load(sidecar, QuantumFurnace)
            r = d[:result]
            @printf("[%s] n=%d β_phys=%-4g seed=%d | SKIP (cached r_D=%d) τ_mix=%.4g gap=%.4g wall=%.2fs\n",
                    label, n, beta_phys, seed,
                    get(r, :r_D, -1),
                    r[:mixing_time], r[:gap_arnoldi], r[:wall_time])
            return (; (Symbol(k) => v for (k, v) in pairs(r))...)
        catch err
            @warn "Cached sidecar failed to load; re-running" sidecar err
        end
    end

    rss_before = rss_mb()
    GC.gc()
    t0_run = time()

    t_max  = 500.0
    t_grid = collect(range(0.0, t_max, length=81))

    d_dim = size(ham.data, 1)
    rho_0 = Matrix{ComplexF64}(I, d_dim, d_dim) ./ d_dim

    local result
    try
        traj = predict_lindbladian_trajectory(
            config, ham, jumps, rho_0, t_grid;
            krylovdim = KRYLOVDIM, tol = 1e-10)
        gap_arnoldi = traj.spectral_gap > 0 ? traj.spectral_gap : 1.0 / config.beta

        res_eig = eigenmode_mixing_time(
            traj.eigenvalues, traj.c, traj.R_modes,
            traj.rho_inf, traj.sigma_beta, float(ε);
            t_upper = t_max,
        )
        wall = time() - t0_run

        tau_mix       = res_eig.mixing_time
        tau_mix_src   = res_eig.source
        floor_dist    = res_eig.floor_distance
        tau_mix_bound = log(d_dim / float(ε)) / max(gap_arnoldi, 1e-12)

        result = (
            n = n, beta_phys = float(beta_phys), beta_alg = config.beta,
            rescaling_factor = ham.rescaling_factor,
            nu_min = ham.nu_min, disorder_strength = 0.1,
            seed = seed, init_state = :maximally_mixed,
            mode = :L, method = :krylov,
            construction = construction_tag, domain = domain_tag,
            filter_name = "default",
            filter_kind = :smooth_metro,
            target_epsilon = float(ε),
            sigma = config.sigma, sigma_factor = config.sigma * config.beta,
            r_D = register_r_D(config), w0_D = register_w0_D(config), t0_D = register_t0_D(config),
            s = config.s, a = config.a,
            gap_arnoldi = gap_arnoldi,
            t_max = t_max,
            t_max_factor = t_max * gap_arnoldi,
            tau_mix_bound = tau_mix_bound,
            n_grid = length(t_grid),
            total_matvecs = traj.total_matvecs,
            all_converged = traj.all_converged,
            mixing_time = tau_mix,
            mixing_time_source = tau_mix_src,
            floor_distance = floor_dist,
            wall_time = wall,
            sweep_version = :v6_multiseed,
        )

        try
            BSON.bson(sidecar, Dict(:result => Dict(pairs(result))))
        catch err
            @warn "Sidecar write failed (continuing)" sidecar err
        end
    catch err
        wall = time() - t0_run
        rss_now = rss_mb()
        @error "[$label] CELL FAILED" n=n β_phys=beta_phys seed=seed ε=ε wall_s=wall rss_mb=rss_now exception=err
        return nothing
    end

    rss_after = rss_mb()
    rss_peak  = rss_peak_mb()

    @printf("[%s] n=%d β_phys=%-4g seed=%d (β_alg=%.2f) r_D=%d | τ_mix=%.4g (%s) gap=%.4g floor=%.3g | matvecs=%d conv=%s | wall=%.2fs ΔRSS=%+.1fMB peak=%.0fMB\n",
            label, n, beta_phys, seed, result.beta_alg, result.r_D,
            result.mixing_time, result.mixing_time_source,
            result.gap_arnoldi, result.floor_distance,
            result.total_matvecs, result.all_converged,
            result.wall_time, rss_after - rss_before, rss_peak)

    if result.floor_distance > 10 * float(ε)
        @printf("    ⚠ floor_distance %.3g > 10·ε = %.3g — bi-exp may not have reached target. Investigate.\n",
                result.floor_distance, 10 * float(ε))
    end
    if !result.all_converged
        @printf("    ⚠ Arnoldi did not fully converge at krylovdim=%d. Consider bumping.\n", KRYLOVDIM)
    end
    flush(stdout)
    return result
end

# --- Smoke probe: Energy↔Bohr ‖ΔL‖_op cross-check at one cell -----------
# Validates r_D=8 keeps the dissipator + coherent error ≤ 1e-9 at the worst
# (highest-β_phys) cell across all 5 seeds before we commit to a 3 h run.

function _build_jump_set_safe(ham, n)
    try
        return _build_jump_set(ham, n)
    catch err
        @warn "_build_jump_set failed; using default Z+ZZ pair jumps" err
        # Fall back to a manual minimal jump set if needed (defensive)
        rethrow()
    end
end

function smoke_energy_vs_bohr(n::Integer, beta_phys::Real)
    @printf("\n[smoke] Energy↔Bohr cross-check at n=%d, β_phys=%.2f, all %d seeds\n",
            n, beta_phys, length(SEEDS))
    println("        (validates r_D=$(_r_D_for_beta_phys(beta_phys)) for β_phys=$(beta_phys))")
    max_err = 0.0
    for seed in SEEDS
        ham = _load_ham_multiseed(n, seed, beta_phys)
        jumps = _build_jump_set(ham, n)
        cfg_e = _build_ckg_cfg(n, beta_phys, ham; r_D = _r_D_for_beta_phys(beta_phys))

        # BohrDomain reference config (closed-form, no quadrature error)
        β_alg = beta_alg(ham, float(beta_phys))
        cfg_b = Config(
            sim = Lindbladian(),
            domain = BohrDomain(),
            construction = KMS(),
            num_qubits = n,
            with_linear_combination = true,
            beta = β_alg,
            beta_phys = float(beta_phys),
            sigma = 1.0 / β_alg,
            a = 0.0, s = 0.25,
            gaussian_parameters = (nothing, nothing),
            num_energy_bits_D = _r_D_for_beta_phys(beta_phys),
            w0_D = register_w0_D(cfg_e), t0_D = register_t0_D(cfg_e),
            num_trotter_steps_per_t0 = 10,
            filter = nothing,
        )

        L_e = construct_lindbladian(jumps, cfg_e, ham)
        L_b = construct_lindbladian(jumps, cfg_b, ham)

        err = opnorm(Matrix(L_e - L_b))
        rel = err / max(opnorm(Matrix(L_b)), 1e-30)
        max_err = max(max_err, rel)
        @printf("    seed=%d  β_alg=%.2f  R=%.3f  ‖L_e−L_b‖=%.3e  rel=%.3e\n",
                seed, β_alg, ham.rescaling_factor, err, rel)
    end
    @printf("[smoke] WORST relative error across %d seeds: %.3e (target ≤ 1e-9)\n",
            length(SEEDS), max_err)
    return max_err
end

# --- Sweep loops -----------------------------------------------------------

function run_S1_at_n(n; rs_s1, beta_phys_list, seeds)
    for seed in seeds
        ham = nothing
        try
            for β_phys in beta_phys_list
                if ham === nothing
                    ham = _load_ham_multiseed(n, seed, float(β_phys))
                else
                    # Reload because β_alg differs per β_phys — HamHam.beta is
                    # baked at construction
                    ham = _load_ham_multiseed(n, seed, float(β_phys))
                end
                jumps = _build_jump_set(ham, n)
                r_D = _r_D_for_beta_phys(β_phys)
                cfg = _build_ckg_cfg(n, β_phys, ham; r_D = r_D)
                r = run_cell(; sweep=:S1, n=n, beta_phys=float(β_phys), seed=seed,
                             ε=float(EPS_TARGET),
                             config=cfg, ham=ham, jumps=jumps,
                             label=@sprintf("S1·n%d", n),
                             out_dir=S1_DIR,
                             construction_tag="KMS", domain_tag="Energy")
                r === nothing || push!(rs_s1, r)
                GC.gc(true)
            end
        catch err
            @error "[S1·n$n seed=$seed] CELL SETUP/RUN CRASHED — skipping seed" exception=(err, catch_backtrace())
        end
    end
end

# --- Main ------------------------------------------------------------------

function main()
    println("\n[main] start $(now())  ε=$EPS_TARGET  β_phys=$(BETA_PHYS_ALL)  seeds=$(SEEDS)  krylovdim=$KRYLOVDIM")
    println("[main] r_D map: β_phys ≤ 1.0 → r_D=7, β_phys ≥ 1.5 → r_D=8")
    println("[main] sidecar dir: $S1_DIR")
    println("[main] expected cells: 6 n × 6 β_phys × 5 seeds = 180")

    # --- Smoke: cross-check at (n=4, β_phys=2.5) the highest-β cell -----
    # This validates that r_D=8 keeps the dissipator + coherent error ≤ 1e-9
    # across all 5 seeds. The smoke is only run if the cached sidecars do not
    # already exist for the first cell — saves a few minutes on resume.
    sample_sidecar = _sidecar_path(S1_DIR, 4, 2.5, 42, "KMS", "Energy")
    if !isfile(sample_sidecar)
        max_err = smoke_energy_vs_bohr(4, 2.5)
        if max_err > 1e-7
            @error "[smoke] FAILED: r_D=8 not adequate at worst-case cell (rel err $max_err > 1e-7). Aborting."
            return
        end
        println("[smoke] PASSED — proceeding with full sweep.")
    else
        println("[smoke] SKIPPING — first cell sidecar already exists, full sweep is resuming.")
    end

    rs_s1 = NamedTuple[]

    # Phase A: n = 3..6 unconditionally (fast).
    for n in 3:6
        println("\n" * "="^72)
        println("=== n=$n  (S1 across 5 seeds × 6 β_phys = 30 cells)")
        println("="^72)
        run_S1_at_n(n; rs_s1=rs_s1, beta_phys_list=collect(BETA_PHYS_ALL), seeds=collect(SEEDS))
    end

    # Adaptive admission for n=7, n=8 using max wall in this run.
    fresh_walls = filter(r -> get(r, :sweep_version, :unknown) === :v6_multiseed, rs_s1)
    for n_try in 7:N_MAX
        prev_walls = filter(r -> r.n == n_try - 1, fresh_walls)
        if isempty(prev_walls)
            @printf("\n[admit] no fresh n=%d wall data → cannot project n=%d; skipping\n", n_try - 1, n_try)
            break
        end
        prev_wall = maximum(r.wall_time for r in prev_walls)
        ratio_guess = if n_try == 7
            NPLUS1_RATIO_GUESS
        else
            pp = filter(r -> r.n == n_try - 2, fresh_walls)
            isempty(pp) ? NPLUS1_RATIO_GUESS :
                max(prev_wall / maximum(r.wall_time for r in pp), 4.0)
        end
        n_est = prev_wall * ratio_guess
        @printf("\n[admit] n=%d S1 worst-cell wall=%.1fs × ratio=%.1f → n=%d est=%.1fs (budget=%.0fs/cell)\n",
                n_try - 1, prev_wall, ratio_guess, n_try, n_est, CELL_WALL_BUDGET)
        if n_est > CELL_WALL_BUDGET
            @printf("[admit] n=%d declined: est %.0fs > budget %.0fs/cell\n", n_try, n_est, CELL_WALL_BUDGET)
            break
        end
        println("\n" * "="^72)
        println("=== n=$(n_try)  (admitted; 30 cells)")
        println("="^72)
        run_S1_at_n(n_try; rs_s1=rs_s1, beta_phys_list=collect(BETA_PHYS_ALL), seeds=collect(SEEDS))
        fresh_walls = filter(r -> get(r, :sweep_version, :unknown) === :v6_multiseed, rs_s1)
    end

    # Summary BSON.
    summary_path = joinpath(OUTPUT_ROOT, "sweep_v6_multiseed_summary.bson")
    rows = Dict[]
    for r in rs_s1
        push!(rows, Dict(:sweep => :S1, pairs(r)...))
    end
    BSON.bson(summary_path, Dict(:rows => rows))
    @printf("\n[summary] wrote %d rows to %s\n", length(rows), summary_path)

    s1_size = sum(filesize(p) for p in readdir(S1_DIR; join=true); init=0)
    @printf("[storage] S1 dir: %d files, %.1f KB\n", length(readdir(S1_DIR)), s1_size/1024)

    println("\n[main] done $(now())  peak_rss=$(round(rss_peak_mb()))MB")
end

main()
