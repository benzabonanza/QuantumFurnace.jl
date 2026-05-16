#!/usr/bin/env julia
# scratch_p1_v6_multiseed_plus.jl  (qf-e4z.31, PARTIAL: n=3..7)
#
# Re-runs the qf-e4z.23 multi-seed Heisenberg PBC sweep with the canonical
# rho_0 = |+⟩⟨+|^⊗N initial state (parity-broken; qf-e4z.30) and
# krylovdim = 60. v6 used rho_0 = I/d which gets parity-trapped on
# Z+ZZ-disordered fixtures (P = Z^⊗N commutes with both H and rho_0): the
# Arnoldi factorisation stays inside the P̂-even sub-spectrum and the
# reported gap_arnoldi is the parity-EVEN sub-spectrum eigenvalue, NOT the
# true L gap. The supposed "even/odd-n parity split" in v6 at high β may
# be an I/d-convention artefact rather than real physics.
#
# This sweep settles the interpretation by saving, per cell:
#   • gap_arnoldi_pass1 = abs(real(traj.eigenvalues[2])) — Pass 1 single-seed
#     Arnoldi from vec(|+⟩⟨+|^⊗N). With a parity-broken seed this captures
#     the true L gap (qf-e4z.30: rel_err 2.4e-9 at n=7 krylovdim=60).
#   • gap_arnoldi_pass2 = traj.spectral_gap — the post-qf-e4z.27 patch
#     pass through krylov_spectral_gap with _krylov_default_x0
#     (I/d + 1e-10·H_GUE-traceless) + KrylovKit thick restart. Independent
#     check: should agree with pass1 to machine precision now.
#
# PARTIAL SCOPE: n ∈ {3, 4, 5, 6, 7}. n=8 split out to qf-e4z.32.
#
# Methodology is otherwise identical to scratch_p1_v6_multiseed.jl:
#   • Same fixtures (scripts/output/multiseed_fixtures/heis_xxx_disordered_periodic_n{3..7}_seed{42..46}.bson)
#   • Same r_D mapping (β_phys ≤ 1.0 → r_D=7, β_phys ≥ 1.5 → r_D=8)
#   • Same Smooth Metropolis filter at s = 0.25, a = 0 (qf-yt9 canonical)
#   • Same ε target 1e-3
#   • Same t_grid (0..500, 81 points)
#
# Cells: 5 n × 6 β_phys × 5 seeds = 150
#
# Sidecar dir: scripts/output/sweep_S1_v6_plus_ckg_ideal_multiseed/smooth_metro_eps1e-03/
# Filename: same as v6 (sweep_n{n}_betaphys{β_phys}_seed{seed}_L_KMS_Energy.bson)
#
# Usage (run from repo root):
#   JULIA_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 \
#     julia --project scripts/scratch_p1_v6_multiseed_plus.jl

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

const EPS_TARGET    = 1e-3
const BETA_PHYS_ALL = (0.25, 0.5, 1.0, 1.5, 2.0, 2.5)
const N_MIN, N_MAX  = 3, 7        # PARTIAL: n=8 split out to qf-e4z.32
const SEEDS         = (42, 43, 44, 45, 46)
const KRYLOVDIM     = 60          # qf-e4z.30 canonical
const TAIL_C        = 8.0

# Heterogeneous r_D mapping (matches v6)
_r_D_for_beta_phys(β_phys::Real) = β_phys ≥ 1.5 ? 8 : 7

const S1_DIR = joinpath(OUTPUT_ROOT, "sweep_S1_v6_plus_ckg_ideal_multiseed", "smooth_metro_eps1e-03")
mkpath(S1_DIR)

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

function _sidecar_path(out_dir::AbstractString, n::Integer, beta_phys::Real, seed::Integer,
                       construction_tag::AbstractString, domain_tag::AbstractString)
    bp_str = let s = @sprintf("%.6f", float(beta_phys))
        s = rstrip(s, '0'); s = rstrip(s, '.')
        isempty(s) ? "0" : s
    end
    return joinpath(out_dir,
        "sweep_n$(n)_betaphys$(bp_str)_seed$(seed)_L_$(construction_tag)_$(domain_tag).bson")
end

function _load_ham_multiseed(n::Integer, seed::Integer, beta_phys::Real)
    path = joinpath(FIXTURE_DIR,
        "heis_xxx_disordered_periodic_n$(n)_seed$(seed).bson")
    isfile(path) || error("Missing multiseed fixture: $path")
    raw = _parse_hamiltonian_bson(path)
    return HamHam(raw; beta_phys = float(beta_phys))
end

# --- Per-cell config builder (identical to v6) ----------------------------

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

# --- Per-cell runner ------------------------------------------------------

function run_cell(; sweep, n, beta_phys, seed, ε, config, ham, jumps, label,
                    out_dir, construction_tag, domain_tag, skip_existing::Bool=true)
    sidecar = _sidecar_path(out_dir, n, beta_phys, seed, construction_tag, domain_tag)

    if skip_existing && isfile(sidecar)
        try
            d = BSON.load(sidecar, QuantumFurnace)
            r = d[:result]
            @printf("[%s] n=%d β_phys=%-4g seed=%d | SKIP (cached) τ_mix=%.4g gap1=%.4g gap2=%.4g wall=%.2fs\n",
                    label, n, beta_phys, seed,
                    r[:mixing_time], r[:gap_arnoldi_pass1], r[:gap_arnoldi_pass2], r[:wall_time])
            return (; (Symbol(k) => v for (k, v) in pairs(r))...)
        catch err
            @warn "Cached sidecar failed to load; re-running" sidecar err
        end
    end

    t0_run = time()
    t_max  = 500.0
    t_grid = collect(range(0.0, t_max, length=81))

    rho_0 = rho_plus_tensor(n)             # qf-e4z.30 canonical seed

    local result
    try
        traj = predict_lindbladian_trajectory(
            config, ham, jumps, rho_0, t_grid;
            krylovdim = KRYLOVDIM, tol = 1e-10)

        # Pass 1: single-seed Arnoldi from vec(|+⟩⟨+|^⊗N). With a
        # parity-broken seed this captures the true L gap (qf-e4z.30).
        gap_pass1 = abs(real(traj.eigenvalues[2]))
        # Pass 2: krylov_spectral_gap with _krylov_default_x0 — the
        # post-qf-e4z.27 patch. Independent reference.
        gap_pass2 = traj.spectral_gap > 0 ? traj.spectral_gap : 1.0 / config.beta

        # Use Pass 1 for the canonical gap_arnoldi field (matches qf-e4z.30
        # canonical recipe), but save Pass 2 alongside for cross-check.
        gap_arnoldi = gap_pass1

        res_eig = eigenmode_mixing_time(
            traj.eigenvalues, traj.c, traj.R_modes,
            traj.rho_inf, traj.sigma_beta, float(ε);
            t_upper = t_max,
        )
        wall = time() - t0_run
        d_dim = size(ham.data, 1)

        tau_mix       = res_eig.mixing_time
        tau_mix_src   = res_eig.source
        floor_dist    = res_eig.floor_distance
        tau_mix_bound = log(d_dim / float(ε)) / max(gap_arnoldi, 1e-12)

        # Relative agreement between the two passes (sanity check; should be
        # ~1e-9 or better on parity-broken seed).
        gap_rel_diff = abs(gap_pass1 - gap_pass2) / max(abs(gap_pass2), 1e-30)

        result = (
            n = n, beta_phys = float(beta_phys), beta_alg = config.beta,
            rescaling_factor = ham.rescaling_factor,
            nu_min = ham.nu_min, disorder_strength = 0.1,
            seed = seed, init_state = :plus_tensor,
            mode = :L, method = :krylov,
            construction = construction_tag, domain = domain_tag,
            filter_name = "default",
            filter_kind = :smooth_metro,
            target_epsilon = float(ε),
            sigma = config.sigma, sigma_factor = config.sigma * config.beta,
            r_D = register_r_D(config), w0_D = register_w0_D(config), t0_D = register_t0_D(config),
            s = config.s, a = config.a,
            gap_arnoldi = gap_arnoldi,
            gap_arnoldi_pass1 = gap_pass1,
            gap_arnoldi_pass2 = gap_pass2,
            gap_rel_diff_pass12 = gap_rel_diff,
            krylovdim = KRYLOVDIM,
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
            sweep_version = :v6_plus_multiseed,
        )

        try
            BSON.bson(sidecar, Dict(:result => Dict(pairs(result))))
        catch err
            @warn "Sidecar write failed (continuing)" sidecar err
        end
    catch err
        wall = time() - t0_run
        @error "[$label] CELL FAILED" n=n β_phys=beta_phys seed=seed ε=ε wall_s=wall exception=err
        return nothing
    end

    @printf("[%s] n=%d β_phys=%-4g seed=%d (β_alg=%.2f) r_D=%d | τ_mix=%.4g (%s) gap1=%.4g gap2=%.4g rel_diff=%.2e | matvecs=%d conv=%s | wall=%.2fs\n",
            label, n, beta_phys, seed, result.beta_alg, result.r_D,
            result.mixing_time, result.mixing_time_source,
            result.gap_arnoldi_pass1, result.gap_arnoldi_pass2,
            result.gap_rel_diff_pass12,
            result.total_matvecs, result.all_converged,
            result.wall_time)

    if result.gap_rel_diff_pass12 > 1e-5
        @printf("    ⚠ Pass1↔Pass2 gap disagreement %.2e > 1e-5 — Krylov truncation residual at this n?\n",
                result.gap_rel_diff_pass12)
    end
    if result.floor_distance > 10 * float(ε)
        @printf("    ⚠ floor_distance %.3g > 10·ε = %.3g — bi-exp may not have reached target.\n",
                result.floor_distance, 10 * float(ε))
    end
    if !result.all_converged
        @printf("    ⚠ Arnoldi did not fully converge at krylovdim=%d. Consider bumping.\n", KRYLOVDIM)
    end
    flush(stdout)
    return result
end

# --- Sweep loops -----------------------------------------------------------

function run_S1_at_n(n; rs_s1, beta_phys_list, seeds)
    for seed in seeds
        try
            for β_phys in beta_phys_list
                ham = _load_ham_multiseed(n, seed, float(β_phys))
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
            @error "[S1·n$n seed=$seed] CELL SETUP CRASHED — skipping seed" exception=(err, catch_backtrace())
        end
    end
end

# --- Main ------------------------------------------------------------------

function main()
    println("\n[main] start $(now())  ε=$EPS_TARGET  β_phys=$(BETA_PHYS_ALL)  seeds=$(SEEDS)  krylovdim=$KRYLOVDIM")
    println("[main] n range: $N_MIN..$N_MAX (PARTIAL — n=8 deferred to qf-e4z.32)")
    println("[main] r_D map: β_phys ≤ 1.0 → r_D=7, β_phys ≥ 1.5 → r_D=8")
    println("[main] init state: rho_0 = |+⟩⟨+|^⊗N (qf-e4z.30 canonical)")
    println("[main] sidecar dir: $S1_DIR")
    println("[main] expected cells: $(N_MAX - N_MIN + 1) n × $(length(BETA_PHYS_ALL)) β_phys × $(length(SEEDS)) seeds = $((N_MAX - N_MIN + 1) * length(BETA_PHYS_ALL) * length(SEEDS))")

    rs_s1 = NamedTuple[]
    for n in N_MIN:N_MAX
        println("\n" * "="^72)
        println("=== n=$n  (S1 across $(length(SEEDS)) seeds × $(length(BETA_PHYS_ALL)) β_phys = $(length(SEEDS)*length(BETA_PHYS_ALL)) cells)")
        println("="^72)
        run_S1_at_n(n; rs_s1=rs_s1, beta_phys_list=collect(BETA_PHYS_ALL), seeds=collect(SEEDS))
    end

    # Summary BSON.
    summary_path = joinpath(OUTPUT_ROOT, "sweep_v6_plus_multiseed_summary.bson")
    rows = Dict[]
    for r in rs_s1
        push!(rows, Dict(:sweep => :S1, pairs(r)...))
    end
    BSON.bson(summary_path, Dict(:rows => rows))
    @printf("\n[summary] wrote %d rows to %s\n", length(rows), summary_path)

    s1_size = sum(filesize(p) for p in readdir(S1_DIR; join=true); init=0)
    @printf("[storage] S1 dir: %d files, %.1f KB\n", length(readdir(S1_DIR)), s1_size/1024)

    println("\n[main] done $(now())")
end

main()
