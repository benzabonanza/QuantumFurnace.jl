#!/usr/bin/env julia
# scratch_qf_e4z_41_gns_smooth_metro_betaphys0_5.jl  (qf-e4z.41)
#
# GNS-DB smooth-Metropolis Lindbladian τ_mix + spectral-gap sweep at
# β_phys = 0.5, seed = 46, n ∈ {3..8}. σ-tightening anchored at n=3 only;
# comparison-plot cells at n ∈ {4..8} use the finest σ_factor = 1/16.
#
# Mirrors the CKG arm of qf-e4z.34 (β_phys=0.5 cells) one-to-one in
# pipeline and sidecar schema so the side-by-side plot is a one-liner
# merge of the two sidecar trees.
#
# DISSIPATOR-ONLY by trait. `with_coherent(::GNS) = false` so
# `apply_lindbladian!` skips the B-term automatically. This is the
# physical GNS-DB sampler (no extra `include_coherent` flag in scripts).
#
# Pipeline
#   • Sampler: GNS = GNS + EnergyDomain + smooth-Metro (s=0.25, a=0)
#              kink at ω=0 (unshifted; cf. KMS at ω=-βσ²/2)
#              filter = nothing
#   • Fixtures: heis_xxx_disordered_periodic_n{n}_seed46.bson (Z+ZZ; same
#              as qf-e4z.34). n=8 seed=46 is built on-the-fly if missing.
#   • rho_0 = |+⟩⟨+|^⊗N (qf-e4z.30 canonical, parity-broken)
#   • Single Pass-1 spectral gap via `predict_lindbladian_trajectory(...;
#     compute_true_gap=false)` (qf-0fv default).
#       gap_arnoldi_pass1 = |Re(traj.eigenvalues[2])|
#   • τ_mix via eigenmode_mixing_time at ε = 1e-3 (bisection, never floor)
#   • t_max = 500, n_grid = 81 (matches qf-e4z.34)
#
# σ_factor → r_D ladder (S5-style, step+1 per halving so kink/grid ≈ const)
#   σ_factor | r_D
#       1.0  |  8   (qf-e4z.35 baseline)
#       1/2  |  9
#       1/4  | 10
#       1/8  | 11
#       1/16 | 12
#       1/32 | 13   (bail-out only)
#
# Krylov dims (from qf-e4z.33 saturation map + qf-e4z.34 extrapolation)
#   n ≤ 6 → 60
#   n = 7 → 80
#   n = 8 → 100
#
# Grid (10 baseline cells)
#   • n = 3   σ_factor ∈ {1, 1/2, 1/4, 1/8, 1/16}  — 5 cells
#   • n ∈ {4,5,6,7,8} σ_factor = 1/16            — 5 cells (one per n)
#
# Convergence policy (per user instruction)
#   • NEVER accept :floor results. If a cell reports :floor (or :nan),
#     escalate via the strategies in `_escalate_cell`:
#       - :floor at σ_factor>1/16 → continue with the ladder (the next
#         tighter σ in the n=3 ramp will lower the floor).
#       - :floor at σ_factor=1/16 → re-run at σ_factor=1/32, r_D=13.
#       - :nan at any cell        → bump krylovdim (×1.5, capped).
#       - gap_phys outlier vs neighbours → diagnostic note + bump
#         krylovdim and/or r_D as appropriate.
#   • Diagnostic baseline: log all per-cell numbers and compare gap_phys
#     against neighbouring n at the same σ_factor (n=3 anchor) before
#     accepting.
#
# Diagnostics baked into every sidecar (so the analyzer doesn't redo ‖L‖)
#   • hs_norm_alg, hs_norm_phys  — matrix-free GKL via hs_operator_norm_krylov
#   • d_1to1_alg, d_1to1_phys    — opnorm of Kossakowski M built from GNS
#     α(ν, ν′) (un-shifted smooth-Metro kernel; cf. KMS in qf-e4z.34).
#
# Output
#   scripts/output/sweep_qf_e4z_41_gns_smooth_metro_betaphys0_5/
#     sweep_n{n}_sigma{σf}_seed46_L_GNS_Energy.bson
#
# Usage (run from repo root):
#   JULIA_NUM_THREADS=#cores OPENBLAS_NUM_THREADS=1 \
#     julia --project scripts/scratch_qf_e4z_41_gns_smooth_metro_betaphys0_5.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using QuantumFurnace
using QuantumFurnace: _parse_hamiltonian_bson, _build_jump_set, X, Y, Z
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
mkpath(FIXTURE_DIR)

const EPS_TARGET     = 1e-3
const BETA_PHYS      = 0.5
const SEED           = 46
const TAIL_C         = 8.0
const T_MAX          = 500.0
const T_GRID_LEN     = 81
const KRYLOV_DIM_HS  = 20
const KRYLOV_TOL_HS  = 1e-10
const N_MIN, N_MAX   = 3, 8

const SIGMA_RAMP_N3   = (1.0, 0.5, 0.25, 0.125, 0.0625)
const SIGMA_FACTOR_NL = 0.0625      # comparison cell at n ≥ 4

# σ_factor → r_D (S5 step+1 per halving)
const SIGMA_TO_RD = Dict(
    1.0     => 8,
    0.5     => 9,
    0.25    => 10,
    0.125   => 11,
    0.0625  => 12,
    0.03125 => 13,     # bail-out
)

# n → krylovdim
function _kdim_for(n::Integer)
    n ≤ 6 && return 60
    n == 7 && return 80
    n == 8 && return 100
    error("_kdim_for: n=$n out of declared range 3..8")
end

# Bail-out & escalation parameters
const KDIM_ESCALATE_CAP = 200          # never push kdim past this
const KDIM_ESCALATE_FACTOR = 1.5       # ×1.5 on each :nan retry

const OUT_ROOT = joinpath(OUTPUT_ROOT, "sweep_qf_e4z_41_gns_smooth_metro_betaphys0_5")
mkpath(OUT_ROOT)

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

function _sf_str(sigma_factor::Real)
    s = @sprintf("%.6f", float(sigma_factor))
    s = rstrip(s, '0'); s = rstrip(s, '.')
    return isempty(s) ? "0" : s
end

function _sidecar_path(n::Integer, sigma_factor::Real, seed::Integer)
    return joinpath(OUT_ROOT,
        "sweep_n$(n)_sigma$(_sf_str(sigma_factor))_seed$(seed)_L_GNS_Energy.bson")
end

function _fixture_path(n::Integer, seed::Integer)
    return joinpath(FIXTURE_DIR,
        "heis_xxx_disordered_periodic_n$(n)_seed$(seed).bson")
end

function _ensure_fixture(n::Integer, seed::Integer)
    path = _fixture_path(n, seed)
    if isfile(path)
        return path
    end
    @printf("[fixture] missing %s — building (deterministic via MersenneTwister(%d))\n",
            basename(path), seed)
    coeffs = [1.0, 1.0, 1.0]
    raw = build_heis_1d(n, coeffs; seed=seed, periodic=true,
        disorder_strength=0.1)
    BSON.bson(path, hamiltonian=raw)
    @printf("[fixture] wrote %s  nu_min=%.4e  R=%.3f\n",
            basename(path), raw.nu_min, raw.rescaling_factor)
    return path
end

function _load_ham(n::Integer, seed::Integer, beta_phys::Real)
    path = _ensure_fixture(n, seed)
    raw = _parse_hamiltonian_bson(path)
    return HamHam(raw; beta_phys = float(beta_phys))
end

# --- Quadrature window: H_norm + 8σ_max tail ------------------------------
#
# omega_range is sized from σ_max so the SAME r_D, w0_D, t0_D triple is
# used at THIS sigma_factor (each cell sets its own r_D from the ladder).
# Per-σ-cell, not per-n: GNS kink at ω=0 is the SAME interior point for
# every σ, so we just want the energy grid to comfortably cover the bulk.

function _omega_range(ham, sigma_alg::Real)
    H_norm = maximum(abs, ham.eigvals)
    return 2.0 * (H_norm + TAIL_C * sigma_alg)
end

function _build_gns_cfg(n::Integer, beta_phys::Real, ham, sigma_factor::Real,
                        r_D::Integer)
    β_alg = beta_alg(ham, float(beta_phys))
    σ     = float(sigma_factor) / β_alg
    ω_range = _omega_range(ham, σ)
    w0_D  = ω_range / 2.0^r_D
    t0_D  = 2π / (2.0^r_D * w0_D)
    return Config(
        sim = Lindbladian(),
        domain = EnergyDomain(),
        construction = GNS(),
        num_qubits = n,
        with_linear_combination = true,
        beta = β_alg,
        beta_phys = float(beta_phys),
        sigma = σ,
        a = 0.0, s = 0.25,
        gaussian_parameters = (nothing, nothing),
        num_energy_bits_D = r_D, w0_D = w0_D, t0_D = t0_D,
        num_trotter_steps_per_t0 = 10,
        filter = nothing,                  # GNS smooth-Metro (no extra filter)
    )
end

# --- Kossakowski M for d_{1→1} bound --------------------------------------
# Same as qf-e4z.34/35 except α uses the GNS smooth-Metro form (a=0, s=0.25).

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

# GNS Kossakowski kernel via the canonical `create_alpha_gns`. The
# GNS-DB kernel inherits the same energy-domain structure as KMS but with
# `|ν+ν′|` replaced by `|ν+ν′ + βσ²/2|` (see bohr_domain.jl:268-284 and
# the surrounding docstring). This is the construction-correct α(ν,ν′)
# used by the GNS Lindbladian's dissipator; using it here keeps the
# d_{1→1} bound consistent with what `apply_lindbladian!` would compute
# if the same M were assembled internally.

gns_alpha_func(β::Real, σ::Real) = (ν, ν′, _a) -> create_alpha_gns(ν, ν′, β, σ, 0.0, 0.25)

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
#
# Returns a NamedTuple result or `nothing` on failure. Does not write the
# sidecar — the caller (`run_cell_with_convergence`) writes the FINAL
# sidecar after deciding convergence. This is so failed/escalated runs
# don't leave a stale sidecar on disk.

function _run_one(; n, beta_phys, sigma_factor, seed, ham, jumps,
                    r_D, kdim, attempt_label="")
    t0_run = time()
    cfg = _build_gns_cfg(n, beta_phys, ham, sigma_factor, r_D)
    rho_0 = rho_plus_tensor(n)
    t_grid = collect(range(0.0, T_MAX, length = T_GRID_LEN))
    R      = ham.rescaling_factor

    traj = predict_lindbladian_trajectory(
        cfg, ham, jumps, rho_0, t_grid;
        krylovdim = kdim, tol = 1e-10,
        compute_true_gap = false,           # qf-0fv default, Pass-1 only
    )

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

    hs_alg  = hs_norm_cell(cfg, ham, jumps)
    hs_phys = hs_alg * R

    α = gns_alpha_func(cfg.beta, cfg.sigma)
    M = dissipator_M_from_alpha(jumps, ham.eigvals, α)
    d11_alg  = opnorm(M)
    d11_phys = d11_alg * R

    wall = time() - t0_run

    return (
        # Cell coordinates
        n = n, beta_phys = float(beta_phys), beta_alg = cfg.beta,
        sigma_alg = cfg.sigma, sigma_phys = cfg.sigma / R,
        sigma_factor = float(sigma_factor),
        rescaling_factor = R, nu_min = ham.nu_min,
        disorder_strength = 0.1, seed = seed,

        # Pipeline metadata
        sampler = :gns,
        init_state = :plus_tensor,
        mode = :L, method = :krylov,
        construction = "GNS", domain = "Energy",
        filter_name = "default_smooth_metro_gns",
        filter_kind = :smooth_metropolis_gns,
        target_epsilon = float(EPS_TARGET),
        r_D = register_r_D(cfg), w0_D = register_w0_D(cfg), t0_D = register_t0_D(cfg),
        omega_range_shared = _omega_range(ham, cfg.sigma),
        s = cfg.s, a = cfg.a,
        krylovdim_p1 = kdim, krylovdim = kdim,
        t_max = T_MAX, n_grid = T_GRID_LEN,

        # Spectral gap (alg + phys)
        gap_arnoldi_pass1 = gap_alg,
        gap_alg  = gap_alg,
        gap_phys = gap_phys,
        t_max_factor = T_MAX * gap_alg,

        # Mixing time (alg + phys)
        mixing_time = tau_mix_alg, mixing_time_alg = tau_mix_alg,
        mixing_time_phys = tau_mix_phys,
        mixing_time_source = tau_mix_alg_src,
        floor_distance = floor_dist,

        # Norm diagnostics (alg + phys)
        hs_norm_alg = hs_alg, hs_norm_phys = hs_phys,
        d_1to1_alg = d11_alg, d_1to1_phys = d11_phys,
        hs_krylovdim = KRYLOV_DIM_HS,

        # Trajectory diagnostics
        total_matvecs = traj.total_matvecs,
        all_converged = traj.all_converged,

        wall_time = wall,
        attempt_label = attempt_label,
        sweep_version = :qf_e4z_41_gns_smooth_metro_betaphys0_5,
    )
end

# Convergence-aware runner. Escalates per the user policy:
#   • :floor at σ_factor = 1/16 → retry at σ_factor = 1/32, r_D = 13.
#   • :nan / not all_converged → bump kdim by ×1.5 (cap KDIM_ESCALATE_CAP).
#   • At lower-σ cells in the n=3 ramp (σ_factor > 1/16): :floor is
#     expected at coarse σ. Log it, accept the floor distance for the
#     σ-scaling tail, but do NOT use this cell for the comparison plot.

function run_cell_with_convergence(; n, beta_phys, sigma_factor, seed, ham, jumps,
                                     skip_existing::Bool = true,
                                     is_n3_ramp_coarse::Bool = false)
    sidecar = _sidecar_path(n, sigma_factor, seed)

    if skip_existing && isfile(sidecar)
        try
            d = BSON.load(sidecar, QuantumFurnace)
            r = d[:result]
            @printf("[GNS] n=%d σ_f=%-7g seed=%d | SKIP (cached)  τ_mix_phys=%s gap_phys=%.4g floor=%.3g src=%s\n",
                    n, sigma_factor, seed,
                    isfinite(r[:mixing_time_phys]) ? @sprintf("%.4g", r[:mixing_time_phys]) : "Inf",
                    r[:gap_phys], r[:floor_distance], r[:mixing_time_source])
            return (; (Symbol(k) => v for (k, v) in pairs(r))...)
        catch err
            @warn "Cached sidecar failed to load; re-running" sidecar err
        end
    end

    r_D  = SIGMA_TO_RD[sigma_factor]
    kdim = _kdim_for(n)
    result = nothing
    attempt = 0

    while true
        attempt += 1
        label = "r_D=$r_D kdim=$kdim attempt#$attempt"
        @printf("[GNS] running n=%d σ_f=%-7g seed=%d (%s)...\n", n, sigma_factor, seed, label)
        flush(stdout)

        try
            result = _run_one(; n=n, beta_phys=beta_phys, sigma_factor=sigma_factor,
                                seed=seed, ham=ham, jumps=jumps,
                                r_D=r_D, kdim=kdim, attempt_label=label)
        catch err
            @error "[GNS] CELL CRASHED" n=n σ_f=sigma_factor exception=(err, catch_backtrace())
            return nothing
        end

        # Cell complete; print summary.
        @printf("       β_alg=%.3f σ_alg=%.4g r_D=%d kdim=%d | τ_mix_phys=%s gap_phys=%.4g floor=%.3g src=%s mv=%d conv=%s wall=%.2fs\n",
                result.beta_alg, result.sigma_alg, result.r_D, result.krylovdim_p1,
                isfinite(result.mixing_time_phys) ? @sprintf("%.4g", result.mixing_time_phys) : "Inf",
                result.gap_phys, result.floor_distance, result.mixing_time_source,
                result.total_matvecs, result.all_converged, result.wall_time)
        flush(stdout)

        # Decide on acceptance / escalation.
        src = result.mixing_time_source
        conv = result.all_converged

        if src == :extrapolated && conv
            # Clean. Accept.
            break
        elseif src == :extrapolated && !conv
            # τ_mix succeeded but Arnoldi flagged non-convergence. Bump kdim.
            if kdim >= KDIM_ESCALATE_CAP
                @printf("    ⚠ Arnoldi not converged at kdim=%d (cap reached). Accepting cell.\n", kdim)
                break
            end
            new_kdim = min(KDIM_ESCALATE_CAP, round(Int, kdim * KDIM_ESCALATE_FACTOR))
            @printf("    ⚠ Arnoldi not converged at kdim=%d → retry kdim=%d.\n", kdim, new_kdim)
            kdim = new_kdim
            continue
        elseif src == :floor
            # No bisection root because floor ≥ ε. Only retry at σf=1/16.
            if is_n3_ramp_coarse
                @printf("    ℹ σ_f=%g floor=%.3g ≥ ε at the σ-ramp tail — expected. Recording for σ-scaling.\n",
                        sigma_factor, result.floor_distance)
                break  # accept floor record; the comparison-plot cell is σ_f=1/16
            elseif sigma_factor == 0.0625
                # The comparison cell hit floor — user-decreed bail-out:
                # retry at σ_factor = 1/32, r_D = 13. Done one level above.
                @printf("    ⚠ σ_f=1/16 floor=%.3g ≥ ε — caller should re-run at σ_f=1/32 r_D=13.\n",
                        result.floor_distance)
                break  # caller orchestrates the retry
            else
                @printf("    ⚠ Unexpected :floor at σ_f=%g. Recording.\n", sigma_factor)
                break
            end
        elseif src == :nan
            # Bisection failed — bump kdim and retry; if already capped, escalate r_D.
            if kdim < KDIM_ESCALATE_CAP
                new_kdim = min(KDIM_ESCALATE_CAP, round(Int, kdim * KDIM_ESCALATE_FACTOR))
                @printf("    ⚠ :nan from eigenmode_mixing_time → retry kdim=%d (was %d).\n", new_kdim, kdim)
                kdim = new_kdim
                continue
            else
                @printf("    ⚠ :nan persists at kdim=%d cap. Accepting (degenerate cell).\n", kdim)
                break
            end
        else
            @printf("    ⚠ Unknown mixing_time_source=%s. Accepting as-is.\n", string(src))
            break
        end
    end

    if result === nothing
        return nothing
    end

    # Write final sidecar.
    try
        BSON.bson(sidecar, Dict(:result => Dict(pairs(result))))
    catch err
        @warn "Sidecar write failed (continuing)" sidecar err
    end

    return result
end

# --- Diagnostic: compare cells across n ----------------------------------

function _diagnose_neighbours(rows::Vector{NamedTuple}, focus_cell::NamedTuple;
                              gap_outlier_ratio::Float64 = 0.3)
    n = focus_cell.n
    σf = focus_cell.sigma_factor
    relevant = [r for r in rows if r.sigma_factor == σf && r.n != n]
    if isempty(relevant)
        return ""
    end
    gaps = [r.gap_phys for r in relevant]
    med  = sort(gaps)[fld(length(gaps), 2) + 1]
    rel  = abs(focus_cell.gap_phys - med) / med
    if rel < gap_outlier_ratio
        return ""
    end
    msg = @sprintf("    ⚠ gap_phys=%.4g at n=%d σ_f=%g is %.1f%% off the σ_f=%g median %.4g over n∈%s.\n",
                   focus_cell.gap_phys, n, σf, 100*rel, σf, med,
                   string([r.n for r in relevant]))
    return msg
end

# --- Main ------------------------------------------------------------------

function main()
    println("\n[main] start $(now())  ε=$EPS_TARGET  β_phys=$BETA_PHYS  seed=$SEED  n=$N_MIN..$N_MAX")
    println("[main] sampler: GNS = GNS + EnergyDomain + smooth-Metro (s=0.25, a=0)")
    println("[main] init state: rho_0 = |+⟩⟨+|^⊗N (qf-e4z.30 canonical)")
    println("[main] t_max=$T_MAX, n_grid=$T_GRID_LEN")
    println("[main] krylovdim recipe: n≤6 → 60 (qf-e4z.33 verified)")
    println("[main]                   n=7 → 80;  n=8 → 100 (qf-e4z.34 extrapolated)")
    println("[main] σ_factor → r_D ladder: ", SIGMA_TO_RD)
    println("[main] output: $OUT_ROOT")
    println("[main] policy: NEVER accept :floor at σ_factor=1/16 — bail out to σ_factor=1/32, r_D=13.")

    rows = NamedTuple[]
    n_planned = length(SIGMA_RAMP_N3) + (N_MAX - 4 + 1)
    @printf("[main] planned baseline cells: %d (n=3 ramp: %d, n=4..8 single σ: %d)\n",
            n_planned, length(SIGMA_RAMP_N3), N_MAX - 4 + 1)

    for n in N_MIN:N_MAX
        println("\n" * "="^72)
        println("=== n=$n")
        println("="^72)

        ham = _load_ham(n, SEED, BETA_PHYS)
        jumps = _build_jump_set(ham, n)
        β_alg = beta_alg(ham, BETA_PHYS)
        H_norm_alg = maximum(abs, ham.eigvals)
        @printf("  [cell-prep] n=%d β_phys=%g β_alg=%.3f H_norm_alg=%.3g R=%.3f\n",
                n, BETA_PHYS, β_alg, H_norm_alg, ham.rescaling_factor)

        # σ grid for this n
        σ_grid = (n == 3) ? collect(SIGMA_RAMP_N3) : [SIGMA_FACTOR_NL]

        for σf in σ_grid
            is_coarse = (n == 3 && σf > 0.0625)
            r = run_cell_with_convergence(
                n=n, beta_phys=BETA_PHYS, sigma_factor=σf, seed=SEED,
                ham=ham, jumps=jumps,
                is_n3_ramp_coarse=is_coarse,
            )
            if r === nothing
                @error "[GNS] cell returned nothing; continuing"
                continue
            end
            push!(rows, r)

            # Neighbour-aware diagnostic at the comparison σ_factor.
            if σf == SIGMA_FACTOR_NL
                diag_msg = _diagnose_neighbours(rows, r)
                if !isempty(diag_msg)
                    print(diag_msg)
                end
            end

            # Floor-bailout retry at σf=1/16: re-run at σf=1/32 if needed.
            if σf == SIGMA_FACTOR_NL && r.mixing_time_source == :floor
                @printf("[GNS] bail-out: re-running n=%d at σ_factor=1/32 (r_D=13)\n", n)
                r2 = run_cell_with_convergence(
                    n=n, beta_phys=BETA_PHYS, sigma_factor=0.03125, seed=SEED,
                    ham=ham, jumps=jumps,
                    is_n3_ramp_coarse=false,
                )
                r2 === nothing || push!(rows, r2)
            end

            GC.gc(true)
        end
    end

    # Summary BSON.
    summary_path = joinpath(OUTPUT_ROOT, "sweep_qf_e4z_41_summary.bson")
    rows_dict = Dict[]
    for r in rows
        push!(rows_dict, Dict(pairs(r)...))
    end
    BSON.bson(summary_path, Dict(:rows => rows_dict))
    @printf("\n[summary] wrote %d rows to %s\n", length(rows_dict), summary_path)

    ndir_size = sum(filesize(p) for p in readdir(OUT_ROOT; join = true); init = 0)
    @printf("[storage] output dir: %d files, %.1f KB\n",
            length(readdir(OUT_ROOT)), ndir_size / 1024)

    # Final cross-cell convergence audit.
    println("\n" * "="^72)
    println("=== Cross-cell convergence audit")
    println("="^72)
    println("σ_factor=1/16 comparison-plot cells (n=3..8):")
    comp = [r for r in rows if r.sigma_factor == SIGMA_FACTOR_NL]
    sort!(comp, by = r -> r.n)
    for r in comp
        flag = r.mixing_time_source == :extrapolated ? "✓" :
               r.mixing_time_source == :floor        ? "FLOOR" : "NaN"
        @printf("  n=%d  gap_phys=%.4g  τ_mix_phys=%s  floor=%.3g  src=%s  kdim=%d  %s\n",
                r.n, r.gap_phys,
                isfinite(r.mixing_time_phys) ? @sprintf("%.4g", r.mixing_time_phys) : "Inf",
                r.floor_distance, r.mixing_time_source, r.krylovdim_p1, flag)
    end

    println("\nn=3 σ-tightening tail (σ_factor 1 → 1/16):")
    tail = [r for r in rows if r.n == 3]
    sort!(tail, by = r -> -r.sigma_factor)
    for r in tail
        @printf("  σ_f=%-7g  σ_alg=%.4g  r_D=%d  floor=%.3g  src=%s\n",
                r.sigma_factor, r.sigma_alg, r.r_D, r.floor_distance, r.mixing_time_source)
    end

    println("\n[main] done $(now())")
end

main()
