#!/usr/bin/env julia
# scratch_qf_e4z_35_sigma_sweep_plot.jl  (qf-e4z.35)
#
# Plot-grade σ-sweep for the numerics chapter:
#   τ_mix (and gap) of the smooth-Metro CKG sampler as a function of σ at
#   fixed β_phys. Same canonical pipeline as qf-e4z.34 (CKG-vs-DLL plot
#   sweep); only the σ-axis changes and we bake HS-norm + d_{1→1}
#   diagnostics into every sidecar.
#
# Pipeline (locked, same as qf-e4z.34)
#   • Sampler: CKG = KMS + EnergyDomain + smooth-Metro (s=0.25, a=0)
#              filter = nothing
#   • Fixtures: heis_xxx_disordered_periodic_n{n}_seed{seed}.bson (Z+ZZ)
#              (legacy disorder; canonical per CLAUDE.md / [[canonical-taumix-setup-qf-e4z-30]])
#   • rho_0 = |+⟩⟨+|^⊗N  (qf-e4z.30 canonical, parity-broken)
#   • Single Arnoldi pass via predict_lindbladian_trajectory:
#       gap_alg = abs(real(traj.eigenvalues[2]))
#       τ_mix  via eigenmode_mixing_time at ε=1e-3
#   • t_max=500, n_grid=81 (matches qf-e4z.34)
#   • r_D = 8 baseline (raised from issue spec's r_D=7 after the floor-cell
#     diagnostic in scratch_qf_e4z_35_diagnostic.jl showed r_D=7 yields
#     floor_distance ≈ 2e-3 at (n=5, β_phys=1, c=0.25) — quadrature bias
#     of the EnergyDomain fixed point exceeds ε. r_D=8 drops the floor to
#     ≤ 1e-5 at the same cell; r_D=9 to ≤ 1e-11. r_D=8 is the cheapest
#     setting that uniformly satisfies floor < ε across the c∈[0.25, 2]
#     and β_phys ∈ {0.25, 0.5, 1.0} grid.
#   • w0_D sized from (H_norm + 8·σ_max) per (n, β_phys) so the SAME
#     energy grid covers every c value in the sweep — fixes the quadrature
#     reference frame as σ varies.
#   • krylovdim recipe (same as qf-e4z.34):
#       n ≤ 6 → 60,   n = 7 → 80,   n = 8 → 80
#
# Diagnostics baked into every sidecar (so the analyzer doesn't redo ‖L‖)
#   • hs_norm_alg, hs_norm_phys  — matrix-free GKL via hs_operator_norm_krylov
#     on (apply_lindbladian!, apply_adjoint_lindbladian!), krylovdim=20
#   • d_1to1_alg, d_1to1_phys    — opnorm of Kossakowski M built from
#     dissipator_M_from_alpha (same recipe as scratch_qf_e4z_34_norm_diagnostic.jl).
#     Bound on the DISSIPATOR only; coherent (i·B) part lives in ‖L‖_HS.
#
# Grid (revised 2026-05-17, single-seed scope)
#   • n ∈ {3, 4, 5, 6, 7, 8}                   6 system sizes
#   • β_phys ∈ {0.25, 0.5, 1.0}                canonical CLAUDE.md grid
#   • σ_factor (c) ∈ {0.25, 0.5, 0.75, 1.0, 1.5, 2.0}    σ_alg = c · (1/β_alg)
#   • seed = 46 only                           (canonical "median-proxy" — qf-e4z.34)
#   • Total cells: 6 × 3 × 6 × 1 = 108
#
# σ-cap policy
#   For the canonical β_phys ≤ 1 grid and the per-n β_alg ≈ R(n)·β_phys, we
#   have β_alg in [0.7 (n=3 β=0.25), ~25 (n=8 β=1)]. σ_alg = c/β_alg, so the
#   max σ_alg is at c=2, β_alg=0.7 → ~2.86, well above H_norm_alg ≈ 0.45 —
#   but the kink is still well-resolved inside the omega_range. We DO NOT
#   cap the c-grid: the spec is six values everywhere.
#
# Physical units
#   • gap_phys         = gap_alg × ham.rescaling_factor
#   • mixing_time_phys = mixing_time_alg / ham.rescaling_factor
#   • hs_norm_phys     = hs_norm_alg × R
#   • d_1to1_phys      = d_1to1_alg  × R
#   • sigma_phys       = sigma_alg / R  (β_phys = β_alg / R ⇒ same scaling)
#
# Output
#   scripts/output/sweep_qf_e4z_35_sigma_sweep_plot/ckg/
#     sweep_n{n}_betaphys{β}_sigma{c}_seed{seed}_L_KMS_Energy.bson
#
# Note on qf-e4z.34 reuse
#   The issue draft mentions "re-using qf-e4z.34 sidecars at c=1". We
#   intentionally re-run that column fresh here because qf-e4z.34's w0_D
#   was sized for σ = 1/β_alg only (TAIL_C·σ = 8/β_alg), whereas qf-e4z.35
#   sizes w0_D from σ_max = 2/β_alg (TAIL_C·σ_max = 16/β_alg). The two
#   quadratures differ at the per-mille level; mixing them in the same
#   plot is not clean. Cost of fresh c=1 column: 18 cells, ~3 minutes.
#
# Usage (run from repo root):
#   JULIA_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 \
#     julia --project scripts/scratch_qf_e4z_35_sigma_sweep_plot.jl

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
println("[init] Julia threads = ", Threads.nthreads(),
        ", BLAS threads = ", BLAS.get_num_threads())
println("[init] hostname = ", gethostname(), "   ", now())

# --- Constants -----------------------------------------------------------

const OUTPUT_ROOT = joinpath(@__DIR__, "output")
const FIXTURE_DIR = joinpath(OUTPUT_ROOT, "multiseed_fixtures")

const EPS_TARGET    = 1e-3
const BETA_PHYS_ALL = (0.25, 0.5, 1.0)
const SIGMA_FACTORS = (0.25, 0.5, 0.75, 1.0, 1.5, 2.0)
const SIGMA_FACTOR_MAX = maximum(SIGMA_FACTORS)
const N_MIN, N_MAX  = 3, 8
const SEED          = 46
const TAIL_C        = 8.0
const R_D           = 8                # raised from 7 after the qf-e4z.35 floor-cell diagnostic

const T_MAX         = 500.0
const T_GRID_LEN    = 81
const KRYLOV_DIM_HS = 20
const KRYLOV_TOL_HS = 1e-10

const OUT_ROOT = joinpath(OUTPUT_ROOT, "sweep_qf_e4z_35_sigma_sweep_plot")
const OUT_CKG  = joinpath(OUT_ROOT, "ckg")
mkpath(OUT_CKG)

# Flat krylovdim per n; verified saturating at qf-e4z.33.
function _krylovdim_for(n::Integer)
    if n ≤ 6
        return 60
    elseif n == 7 || n == 8
        return 80
    else
        error("_krylovdim_for: n=$n out of declared range 3..8")
    end
end

# --- |+⟩⟨+|^⊗N density matrix (parity-broken seed; qf-e4z.30 canonical) ---

function rho_plus_tensor(n::Integer)
    plus = ComplexF64[0.5 0.5; 0.5 0.5]
    rho = plus
    for _ in 2:n
        rho = kron(rho, plus)
    end
    return rho
end

# --- Sidecar / fixture paths ----------------------------------------------

function _bp_str(beta_phys::Real)
    s = @sprintf("%.6f", float(beta_phys))
    s = rstrip(s, '0'); s = rstrip(s, '.')
    return isempty(s) ? "0" : s
end

function _c_str(sigma_factor::Real)
    s = @sprintf("%.4f", float(sigma_factor))
    s = rstrip(s, '0'); s = rstrip(s, '.')
    return isempty(s) ? "0" : s
end

function _sidecar_path(n::Integer, beta_phys::Real, sigma_factor::Real, seed::Integer)
    return joinpath(OUT_CKG,
        "sweep_n$(n)_betaphys$(_bp_str(beta_phys))_sigma$(_c_str(sigma_factor))_seed$(seed)_L_KMS_Energy.bson")
end

function _load_ham_multiseed(n::Integer, seed::Integer, beta_phys::Real)
    path = joinpath(FIXTURE_DIR,
        "heis_xxx_disordered_periodic_n$(n)_seed$(seed).bson")
    isfile(path) || error("Missing multiseed fixture: $path")
    raw = _parse_hamiltonian_bson(path)
    return HamHam(raw; beta_phys = float(beta_phys))
end

# --- Quadrature register for the (n, β_phys) cell -------------------------
#
# omega_range is sized from σ_max so the same r_D, w0_D, t0_D triple is
# used for every c in the c-grid at this (n, β_phys). This fixes the
# energy discretization while sweeping σ — only α(ν, ν′) (the smooth-Metro
# kernel) changes between cells.

function _shared_omega_range(ham, beta_phys::Real)
    β_alg = beta_alg(ham, float(beta_phys))
    σ_max = SIGMA_FACTOR_MAX / β_alg
    H_norm = maximum(abs, ham.eigvals)
    return 2.0 * (H_norm + TAIL_C * σ_max)
end

function _build_ckg_cfg(n::Integer, beta_phys::Real, ham, sigma_factor::Real,
                        omega_range::Real; r_D::Integer = R_D)
    β_alg = beta_alg(ham, float(beta_phys))
    σ     = float(sigma_factor) / β_alg
    w0_D  = omega_range / 2.0^r_D
    t0_D  = 2π / (2.0^r_D * w0_D)
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
        filter = nothing,                  # CKG smooth-Metro
    )
end

# --- Kossakowski M for d_{1→1} bound --------------------------------------
# Ported from scripts/scratch_qf_e4z_34_norm_diagnostic.jl.

function dissipator_M_from_alpha(jumps::Vector{JumpOp},
                                 eigvals::AbstractVector{<:Real},
                                 alpha_func::F) where {F}
    d = length(eigvals)
    M = zeros(ComplexF64, d, d)
    for (a_idx, jump) in enumerate(jumps)
        A = jump.in_eigenbasis
        for j in 1:d, k in 1:d
            s = zero(ComplexF64)
            @inbounds for i in 1:d
                ν  = eigvals[i] - eigvals[j]
                ν′ = eigvals[i] - eigvals[k]
                s += conj(A[i, k]) * A[i, j] * alpha_func(ν, ν′, a_idx)
            end
            M[k, j] += s
        end
    end
    return M
end

ckg_alpha_func(β::Real, σ::Real) = (ν, ν′, _a) -> create_alpha(ν, ν′, β, σ, 0.0, 0.25)

# --- HS norm via matrix-free GKL ------------------------------------------

function hs_norm_cell(cfg::Config{Lindbladian}, ham::HamHam, jumps::Vector{JumpOp})
    dim = size(ham.data, 1)
    ws_fwd = Workspace(cfg, ham, jumps)
    ws_adj = Workspace(cfg, ham, jumps)
    L_apply!     = (out, X) -> (apply_lindbladian!(ws_fwd, X, cfg, ham);
                                copyto!(out, ws_fwd.scratch.rho_out); out)
    L_apply_adj! = (out, X) -> (apply_adjoint_lindbladian!(ws_adj, X, cfg, ham);
                                copyto!(out, ws_adj.scratch.rho_out); out)
    return hs_operator_norm_krylov(L_apply!, L_apply_adj!, dim;
                                   krylovdim = KRYLOV_DIM_HS, tol = KRYLOV_TOL_HS,
                                   maxiter = 100)
end

# --- Per-cell runner ------------------------------------------------------

function run_cell(; n::Integer, beta_phys::Real, sigma_factor::Real,
                    seed::Integer = SEED, ham, jumps, omega_range::Real,
                    skip_existing::Bool = true)
    sidecar = _sidecar_path(n, beta_phys, sigma_factor, seed)

    if skip_existing && isfile(sidecar)
        try
            d = BSON.load(sidecar, QuantumFurnace)
            r = d[:result]
            @printf("[ckg σ-sweep] n=%d β_phys=%-4g c=%-4g seed=%d | SKIP (cached)  τ_mix_phys=%.4g  gap_phys=%.4g  ‖L‖_HS_phys=%.4g  d11_phys=%.4g  wall=%.2fs\n",
                    n, beta_phys, sigma_factor, seed,
                    r[:mixing_time_phys], r[:gap_phys], r[:hs_norm_phys], r[:d_1to1_phys],
                    r[:wall_time])
            return (; (Symbol(k) => v for (k, v) in pairs(r))...)
        catch err
            @warn "Cached sidecar failed to load; re-running" sidecar err
        end
    end

    t0_run = time()
    cfg = _build_ckg_cfg(n, beta_phys, ham, sigma_factor, omega_range)
    rho_0 = rho_plus_tensor(n)
    kdim  = _krylovdim_for(n)
    t_grid = collect(range(0.0, T_MAX, length = T_GRID_LEN))
    R      = ham.rescaling_factor

    local result
    try
        # 1. Trajectory + gap + τ_mix (qf-e4z.30 single-pass recipe).
        traj = predict_lindbladian_trajectory(
            cfg, ham, jumps, rho_0, t_grid;
            krylovdim = kdim, tol = 1e-10)

        gap_alg  = abs(real(traj.eigenvalues[2]))
        gap_phys = gap_alg * R

        res_eig = eigenmode_mixing_time(
            traj.eigenvalues, traj.c, traj.R_modes,
            traj.rho_inf, traj.sigma_beta, float(EPS_TARGET);
            t_upper = T_MAX,
        )
        tau_mix_alg     = res_eig.mixing_time
        tau_mix_alg_src = res_eig.source
        floor_dist      = res_eig.floor_distance
        tau_mix_phys    = isfinite(tau_mix_alg) ? tau_mix_alg / R : tau_mix_alg

        # 2. HS norm of L (matrix-free; alg frame, units 1/time_alg).
        hs_alg  = hs_norm_cell(cfg, ham, jumps)
        hs_phys = hs_alg * R

        # 3. d_{1→1} dissipator bound via Kossakowski M.
        α = ckg_alpha_func(cfg.beta, cfg.sigma)
        M = dissipator_M_from_alpha(jumps, ham.eigvals, α)
        d11_alg  = opnorm(M)
        d11_phys = d11_alg * R

        wall = time() - t0_run

        result = (
            # Cell coordinates
            n = n, beta_phys = float(beta_phys), beta_alg = cfg.beta,
            sigma_alg = cfg.sigma, sigma_phys = cfg.sigma / R,
            sigma_factor = float(sigma_factor),
            rescaling_factor = R, nu_min = ham.nu_min,
            disorder_strength = 0.1, seed = seed,

            # Pipeline metadata
            sampler = :ckg,
            init_state = :plus_tensor,
            mode = :L, method = :krylov,
            construction = "KMS", domain = "Energy",
            filter_name = "default_smooth_metro",
            filter_kind = :smooth_metro,
            target_epsilon = float(EPS_TARGET),
            r_D = register_r_D(cfg), w0_D = register_w0_D(cfg), t0_D = register_t0_D(cfg),
            omega_range_shared = omega_range,
            s = cfg.s, a = cfg.a,
            krylovdim = kdim,
            t_max = T_MAX, n_grid = T_GRID_LEN,

            # Spectral gap (alg + phys)
            gap_alg  = gap_alg,
            gap_phys = gap_phys,
            t_max_factor = T_MAX * gap_alg,

            # Mixing time (alg + phys)
            mixing_time_alg = tau_mix_alg,
            mixing_time_phys = tau_mix_phys,
            mixing_time_source = tau_mix_alg_src,
            floor_distance = floor_dist,

            # Norm diagnostics (alg + phys)
            hs_norm_alg = hs_alg,
            hs_norm_phys = hs_phys,
            d_1to1_alg = d11_alg,
            d_1to1_phys = d11_phys,
            hs_krylovdim = KRYLOV_DIM_HS,

            # Trajectory diagnostics
            total_matvecs = traj.total_matvecs,
            all_converged = traj.all_converged,

            wall_time = wall,
            sweep_version = :qf_e4z_35_sigma_sweep_plot,
        )

        try
            BSON.bson(sidecar, Dict(:result => Dict(pairs(result))))
        catch err
            @warn "Sidecar write failed (continuing)" sidecar err
        end
    catch err
        wall = time() - t0_run
        @error "[ckg σ-sweep] CELL FAILED" n=n β_phys=beta_phys c=sigma_factor seed=seed wall_s=wall exception=(err, catch_backtrace())
        return nothing
    end

    @printf("[ckg σ-sweep] n=%d β_phys=%-4g c=%-4g seed=%d (β_alg=%.2f σ_alg=%.4g) kdim=%d | τ_mix_phys=%.4g (%s)  gap_phys=%.4g  ‖L‖_HS_phys=%.4g  d11_phys=%.4g  mv=%d  conv=%s  wall=%.2fs\n",
            n, beta_phys, sigma_factor, seed, result.beta_alg, result.sigma_alg, result.krylovdim,
            result.mixing_time_phys, result.mixing_time_source, result.gap_phys,
            result.hs_norm_phys, result.d_1to1_phys,
            result.total_matvecs, result.all_converged,
            result.wall_time)

    if result.floor_distance > 10 * float(EPS_TARGET)
        @printf("    ⚠ floor_distance %.3g > 10·ε = %.3g — bi-exp may not have reached target.\n",
                result.floor_distance, 10 * float(EPS_TARGET))
    end
    if !result.all_converged
        @printf("    ⚠ Arnoldi did not fully converge at krylovdim=%d. Bump on rerun.\n", result.krylovdim)
    end
    flush(stdout)
    return result
end

# --- Sweep loops -----------------------------------------------------------

function run_at_n(n; rows::Vector{NamedTuple}, beta_phys_list, sigma_factors, seed::Integer)
    for β_phys in beta_phys_list
        # Load the (n, β_phys, seed) fixture once; sweep c at constant grid.
        ham   = _load_ham_multiseed(n, seed, β_phys)
        jumps = _build_jump_set(ham, n)
        omega_range = _shared_omega_range(ham, β_phys)
        β_alg = beta_alg(ham, float(β_phys))
        H_norm_alg = maximum(abs, ham.eigvals)
        σ_min_alg  = minimum(sigma_factors) / β_alg
        σ_max_alg  = maximum(sigma_factors) / β_alg
        w0_D_grid  = omega_range / 2.0^R_D
        @printf("  [cell-prep] n=%d β_phys=%g  β_alg=%.3f  H_norm_alg=%.3g  σ_alg ∈ [%.3g, %.3g]  ω_range=%.3g  w0_D=%.3g  resolution_min(kink/grid)=%.3g\n",
                n, β_phys, β_alg, H_norm_alg, σ_min_alg, σ_max_alg, omega_range, w0_D_grid,
                σ_min_alg * sqrt(0.25) / w0_D_grid)
        flush(stdout)

        for c in sigma_factors
            try
                r = run_cell(; n = n, beta_phys = float(β_phys),
                             sigma_factor = float(c), seed = seed,
                             ham = ham, jumps = jumps,
                             omega_range = omega_range)
                r === nothing || push!(rows, r)
            catch err
                @error "[ckg·n$n β=$β_phys c=$c seed=$seed] CELL SETUP CRASHED" exception=(err, catch_backtrace())
            end
            GC.gc(true)
        end
    end
end

# --- Main ------------------------------------------------------------------

function main()
    println("\n[main] start $(now())  ε=$EPS_TARGET  β_phys=$(BETA_PHYS_ALL)  c=$(SIGMA_FACTORS)  seed=$SEED  n=$N_MIN..$N_MAX")
    println("[main] sampler: CKG = KMS + EnergyDomain + smooth-Metro (s=0.25, a=0), r_D=$R_D, 8σ_max tail")
    println("[main] init state: rho_0 = |+⟩⟨+|^⊗N (qf-e4z.30 canonical)")
    println("[main] t_max=$T_MAX, n_grid=$T_GRID_LEN")
    println("[main] krylovdim recipe: n≤6 → 60 (saturating); n=7,8 → 80")
    println("[main] HS norm via GKL (krylovdim=$KRYLOV_DIM_HS, tol=$KRYLOV_TOL_HS)")
    println("[main] output: $OUT_ROOT")
    n_cells = (N_MAX - N_MIN + 1) * length(BETA_PHYS_ALL) * length(SIGMA_FACTORS)
    println("[main] expected cells: $(N_MAX - N_MIN + 1) n × $(length(BETA_PHYS_ALL)) β × $(length(SIGMA_FACTORS)) c = $n_cells")

    rows = NamedTuple[]
    for n in N_MIN:N_MAX
        println("\n" * "="^72)
        println("=== n=$n  ($(length(SIGMA_FACTORS) * length(BETA_PHYS_ALL)) cells)")
        println("="^72)
        run_at_n(n; rows = rows,
                 beta_phys_list = collect(BETA_PHYS_ALL),
                 sigma_factors  = collect(SIGMA_FACTORS),
                 seed = SEED)
    end

    # Summary BSON.
    summary_path = joinpath(OUTPUT_ROOT, "sweep_qf_e4z_35_sigma_sweep_plot_summary.bson")
    rows_dict = Dict[]
    for r in rows
        push!(rows_dict, Dict(pairs(r)...))
    end
    BSON.bson(summary_path, Dict(:rows => rows_dict))
    @printf("\n[summary] wrote %d rows to %s\n", length(rows_dict), summary_path)

    ndir_size = sum(filesize(p) for p in readdir(OUT_CKG; join = true); init = 0)
    @printf("[storage] ckg dir: %d files, %.1f KB\n",
            length(readdir(OUT_CKG)), ndir_size / 1024)

    println("\n[main] done $(now())")
end

main()
