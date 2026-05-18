#!/usr/bin/env julia
# scratch_qf_e4z_35_diagnostic.jl  (qf-e4z.35)
#
# Focused diagnostic for the σ-sweep "floor cell" anomaly at small σ:
# at n=5, β_phys=1, c=0.25 (σ_alg ≈ 0.0086) the single-pass Arnoldi from
# |+⟩⟨+|^⊗N reported floor_distance > ε ⇒ τ_mix = Inf with source=:floor.
#
# Two competing explanations:
#   (1) Krylov truncation: kdim=60 cannot capture the Lindbladian's
#       fixed point to within ε=1e-3 of sigma_beta.  Bumping kdim should
#       drop floor_distance to ~1e-15.
#   (2) Quadrature error: at small σ the smooth-Metropolis kink kernel
#       in α(ν, ν′) becomes sub-cell on the r_D=7 energy grid, so the
#       computed L has a different fixed point than the true one.
#       Bumping r_D should align the EnergyDomain fixed point with Bohr.
#
# What we check at the cell:
#   (a) Single-pass kdim sweep with r_D = 7: kdim ∈ {60, 80, 100, 140, 200}.
#       Track floor_distance, gap, τ_mix per kdim.
#   (b) r_D sweep at kdim = 100: r_D ∈ {7, 8, 9, 10}. Track the same
#       quantities; gap should converge as r_D grows.
#   (c) BohrDomain reference at the same (n, β, σ): dense `construct_lindbladian`
#       + eigvals (n=5, d=32 → d² = 1024, fits easily). Gap is exact;
#       floor distance is 0 by construction.
#   (d) Repeat (b) at the SHARED quadrature (omega_range from σ_max=2/β_alg,
#       same as the σ-sweep driver) to mirror exactly what the driver does.
#
# Decision rule:
#   • If (a) collapses floor_distance with kdim alone → bump kdim in driver.
#   • If (a) plateaus above ε regardless of kdim, but (b) reduces floor on
#     bigger r_D → bump r_D in driver.
#   • If neither helps and Bohr also produces a non-trivial Krylov-floor-like
#     behaviour ⇒ the small-σ Lindbladian is genuinely near-degenerate.
#
# Usage:
#   JULIA_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 \
#     julia --project scripts/scratch_qf_e4z_35_diagnostic.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using QuantumFurnace
using QuantumFurnace: _parse_hamiltonian_bson, _build_jump_set
using LinearAlgebra
using BSON
using Printf
using Dates

BLAS.set_num_threads(1)
println("[init] Julia threads = ", Threads.nthreads(),
        ", BLAS threads = ", BLAS.get_num_threads())

const OUTPUT_ROOT = joinpath(@__DIR__, "output")
const FIXTURE_DIR = joinpath(OUTPUT_ROOT, "multiseed_fixtures")

# Pick the worst observed cell.
const N         = 5
const BETA_PHYS = 1.0
const C         = 0.25
const SEED      = 46
const TAIL_C    = 8.0
const SIGMA_FACTOR_MAX = 2.0           # mirror driver's σ-sweep upper bound
const T_MAX        = 500.0
const T_GRID_LEN   = 81
const EPS_TARGET   = 1e-3

function rho_plus_tensor(n::Integer)
    plus = ComplexF64[0.5 0.5; 0.5 0.5]
    rho = plus
    for _ in 2:n
        rho = kron(rho, plus)
    end
    return rho
end

function build_cfg(; n::Integer, beta_phys::Real, sigma_factor::Real, r_D::Integer,
                     omega_range_sigma_factor::Real = SIGMA_FACTOR_MAX, ham)
    β_alg = beta_alg(ham, float(beta_phys))
    σ     = sigma_factor / β_alg
    σ_max = omega_range_sigma_factor / β_alg
    H_norm = maximum(abs, ham.eigvals)
    omega_range = 2.0 * (H_norm + TAIL_C * σ_max)
    w0_D = omega_range / 2.0^r_D
    t0_D = 2π / (2.0^r_D * w0_D)
    cfg = Config(
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
    return cfg, σ
end

function build_bohr_cfg(; n::Integer, beta_phys::Real, sigma_factor::Real, ham)
    β_alg = beta_alg(ham, float(beta_phys))
    σ     = sigma_factor / β_alg
    return Config(
        sim = Lindbladian(),
        domain = BohrDomain(),
        construction = KMS(),
        num_qubits = n,
        with_linear_combination = true,
        beta = β_alg,
        beta_phys = float(beta_phys),
        sigma = σ,
        a = 0.0, s = 0.25,
        gaussian_parameters = (nothing, nothing),
        num_energy_bits = 12, w0 = 0.05, t0 = 2π / (2.0^12 * 0.05),  # ignored by Bohr
        num_trotter_steps_per_t0 = 10,
        filter = nothing,
    )
end

function load_ham(n::Integer, seed::Integer, beta_phys::Real)
    path = joinpath(FIXTURE_DIR,
        "heis_xxx_disordered_periodic_n$(n)_seed$(seed).bson")
    return HamHam(_parse_hamiltonian_bson(path); beta_phys = float(beta_phys))
end

function run_single(cfg, ham, jumps, rho_0, kdim::Integer)
    t_grid = collect(range(0.0, T_MAX, length = T_GRID_LEN))
    traj = predict_lindbladian_trajectory(cfg, ham, jumps, rho_0, t_grid;
                                          krylovdim = kdim, tol = 1e-10)
    gap_alg = abs(real(traj.eigenvalues[2]))
    res = eigenmode_mixing_time(traj.eigenvalues, traj.c, traj.R_modes,
                                 traj.rho_inf, traj.sigma_beta, EPS_TARGET;
                                 t_upper = T_MAX)
    return (gap_alg = gap_alg, tau_mix_alg = res.mixing_time,
            source = res.source, floor_distance = res.floor_distance,
            all_converged = traj.all_converged,
            total_matvecs = traj.total_matvecs)
end

function bohr_dense_gap(cfg_bohr, ham, jumps)
    # n = 5, d = 32, d² = 1024, dense Lindbladian = 1024×1024 complex ⇒ 16 MB.
    L = construct_lindbladian(jumps, cfg_bohr, ham)
    λ = eigvals(L)
    # Sort by abs, pick first nonzero
    sort!(λ; by = abs)
    gap_alg = abs(real(λ[2]))                # 2nd smallest = slowest non-stationary mode
    fixed_eval = abs(λ[1])
    return (gap_alg = gap_alg, fixed_eval = fixed_eval)
end

function main()
    println("="^104)
    println("qf-e4z.35 floor-cell diagnostic   n=$N  β_phys=$BETA_PHYS  c=$C  seed=$SEED  ε=$EPS_TARGET")
    println("="^104)
    ham = load_ham(N, SEED, BETA_PHYS)
    jumps = _build_jump_set(ham, N)
    β_alg = beta_alg(ham, BETA_PHYS)
    σ_target = C / β_alg
    H_norm = maximum(abs, ham.eigvals)
    @printf("[setup] β_alg = %.4f   σ_alg = %.6g   σ_alg·√s = %.6g   H_norm_alg = %.4f   R = %.4f\n",
            β_alg, σ_target, σ_target * sqrt(0.25), H_norm, ham.rescaling_factor)
    rho_0 = rho_plus_tensor(N)

    # ============== (a) kdim sweep at r_D = 7 (matches initial driver) ====
    println("\n--- (a) kdim sweep, r_D = 7 (driver baseline), σ-sweep ω_range from σ_max=2/β_alg ---")
    cfg7, _ = build_cfg(; n=N, beta_phys=BETA_PHYS, sigma_factor=C, r_D=7,
                          ham=ham, omega_range_sigma_factor=SIGMA_FACTOR_MAX)
    @printf("[cfg] r_D=7  w0_D=%.4g  grid pts per kink width = %.3g\n",
            register_w0_D(cfg7), σ_target * sqrt(0.25) / register_w0_D(cfg7))
    @printf("%-6s %-12s %-12s %-12s %-8s %-6s\n",
            "kdim", "gap_alg", "tau_mix_alg", "floor_dist", "source", "conv")
    for kdim in (60, 80, 100, 140, 200)
        r = run_single(cfg7, ham, jumps, rho_0, kdim)
        @printf("%-6d %-12.6g %-12.6g %-12.6g %-8s %-6s\n",
                kdim, r.gap_alg, r.tau_mix_alg, r.floor_distance,
                string(r.source), string(r.all_converged))
    end

    # ============== (b) r_D sweep at kdim = 100 ============================
    println("\n--- (b) r_D sweep at kdim=100, σ-sweep ω_range from σ_max=2/β_alg ---")
    @printf("%-6s %-12s %-12s %-12s %-12s %-8s %-6s\n",
            "r_D", "w0_D", "gap_alg", "tau_mix_alg", "floor_dist", "source", "conv")
    for r_D in (7, 8, 9, 10)
        cfg, _ = build_cfg(; n=N, beta_phys=BETA_PHYS, sigma_factor=C, r_D=r_D,
                            ham=ham, omega_range_sigma_factor=SIGMA_FACTOR_MAX)
        r = run_single(cfg, ham, jumps, rho_0, 100)
        @printf("%-6d %-12.4g %-12.6g %-12.6g %-12.6g %-8s %-6s\n",
                r_D, register_w0_D(cfg), r.gap_alg, r.tau_mix_alg,
                r.floor_distance, string(r.source), string(r.all_converged))
    end

    # ============== (c) Bohr dense reference ================================
    println("\n--- (c) BohrDomain dense reference (exact quadrature; n=5 d=32) ---")
    cfg_bohr = build_bohr_cfg(; n=N, beta_phys=BETA_PHYS, sigma_factor=C, ham=ham)
    bohr = bohr_dense_gap(cfg_bohr, ham, jumps)
    @printf("[Bohr dense] gap_alg = %.6g   fixed_eval = %.3g (should be ~0)\n",
            bohr.gap_alg, bohr.fixed_eval)

    # ============== (d) r_D sweep at PER-σ ω_range (tighter) ===============
    # If we don't share ω_range with σ_max, each r_D gets a tighter grid
    # around the c=0.25 kink.  This is NOT what the driver does, but it
    # tells us whether the shared-ω_range choice is causing the under-resolution.
    println("\n--- (d) r_D sweep with PER-σ ω_range (no shared ω_range) — diagnostic only ---")
    @printf("%-6s %-12s %-12s %-12s %-12s %-8s %-6s\n",
            "r_D", "w0_D", "gap_alg", "tau_mix_alg", "floor_dist", "source", "conv")
    for r_D in (7, 8, 9)
        cfg, _ = build_cfg(; n=N, beta_phys=BETA_PHYS, sigma_factor=C, r_D=r_D,
                             ham=ham, omega_range_sigma_factor=C)
        r = run_single(cfg, ham, jumps, rho_0, 100)
        @printf("%-6d %-12.4g %-12.6g %-12.6g %-12.6g %-8s %-6s\n",
                r_D, register_w0_D(cfg), r.gap_alg, r.tau_mix_alg,
                r.floor_distance, string(r.source), string(r.all_converged))
    end

    println("\n[done] $(now())")
end

main()
