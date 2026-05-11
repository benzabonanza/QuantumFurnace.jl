#!/usr/bin/env julia
#
# qf-e4z.20.6 — Independent per-leg Trotter caches controllability proof.
#
# Goal: demonstrate that with the qf-e4z.20 `TrotterTriple` scheme, the
# Lindbladian KMS-DBC residue `‖L_TrotterDomain · σ_β‖_HS` is controllable to
# ≤ 1e-6 at the canonical fixture (n=3, β=10, σ=1/β, smooth_metro, s=0.25) by
# INDEPENDENTLY sizing the per-leg knobs (r_D, r_b_minus, r_b_plus, M_D,
# M_b_minus, M_b_plus) — no commensurability constraint forcing simultaneous
# over-resolution across legs (as in the qf-e4z.5.3 Option A m_D = 1280 recipe).
#
# Cross-check: β = 10 is the algorithm-level inverse temperature (no
# beta_phys / beta_alg split exists in src/ yet, so the bare `Config.beta = 10`
# coincides with the old code's β convention).

using QuantumFurnace
using LinearAlgebra
using Printf

const QF = QuantumFurnace

# ---------------------------------------------------------------------------
# Fixture: n=3 disordered Heisenberg, β=10, smooth Metro (s=0.25, a=β/30).
# ---------------------------------------------------------------------------
const N_QUBITS = 3
const BETA = 10.0
const SIGMA = 1.0 / BETA
const SMOOTH_S = 0.25
const SMOOTH_A = BETA / 30.0
const ETA = 1e-3

const HAM = QF._load_hamiltonian_bson(
    joinpath(@__DIR__, "..", "hamiltonians", "heis_disordered_periodic_n3.bson"), BETA)
const DIM = size(HAM.data, 1)

# ω_range for the dissipative grid: thesis fixture uses ω_range ≈ 2.5.
const OMEGA_RANGE = 2.5
const T_MINUS = 18.0
const T_PLUS = 12.0

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
function jumps_in_basis(V::AbstractMatrix)
    jumps = JumpOp[]
    for pauli in [[X], [Y], [Z]], site in 1:N_QUBITS
        op = Matrix(pad_term(pauli, N_QUBITS, site)) ./ sqrt(3 * N_QUBITS)
        op_eb = V' * op * V
        push!(jumps, JumpOp(op, op_eb, op == transpose(op), op == op'))
    end
    return jumps
end

function make_cfg(; r_D, r_bm, r_bp, M_D, M_bm, M_bp)
    w0_D = OMEGA_RANGE / 2^r_D
    t0_D = 2pi / (2^r_D * w0_D)
    t0_bm = 2 * T_MINUS / 2^r_bm
    w0_bm = 2pi / (2^r_bm * t0_bm)
    t0_bp = 2 * T_PLUS / 2^r_bp
    w0_bp = 2pi / (2^r_bp * t0_bp)
    return Config(
        sim = Lindbladian(),
        domain = TrotterDomain(),
        construction = KMS(),
        num_qubits = N_QUBITS,
        with_linear_combination = true,
        beta = BETA,
        sigma = SIGMA,
        s = SMOOTH_S,
        a = SMOOTH_A,
        eta = ETA,
        num_energy_bits_D = r_D, t0_D = t0_D, w0_D = w0_D,
        num_energy_bits_b_minus = r_bm, t0_b_minus = t0_bm, w0_b_minus = w0_bm,
        num_energy_bits_b_plus = r_bp,  t0_b_plus  = t0_bp,  w0_b_plus = w0_bp,
        num_trotter_steps_per_t0 = M_D,
        num_trotter_steps_per_t0_b_minus = M_bm,
        num_trotter_steps_per_t0_b_plus  = M_bp,
    )
end

function residue_at(; r_D, r_bm, r_bp, M_D, M_bm, M_bp)
    cfg = make_cfg(; r_D, r_bm, r_bp, M_D, M_bm, M_bp)
    validate_config!(cfg)
    trotter = make_trotter_for_config(HAM, cfg)
    jumps = jumps_in_basis(trotter.eigvecs)
    L = construct_lindbladian(jumps, cfg, HAM; trotter=trotter)
    sigma_beta = Hermitian(trotter.eigvecs' * HAM.eigvecs * HAM.gibbs *
                            HAM.eigvecs' * trotter.eigvecs)
    return norm(L * vec(Matrix(sigma_beta)))
end

println("="^72)
println("qf-e4z.20.6 — Independent per-leg Trotter caches controllability proof")
println("="^72)
@printf("Fixture: n=%d, β=%.1f, σ=%.4f, smooth Metropolis (s=%.2f, a=%.4f)\n",
        N_QUBITS, BETA, SIGMA, SMOOTH_S, SMOOTH_A)
@printf("Hamiltonian: %s\n",
        joinpath("hamiltonians", "heis_disordered_periodic_n3.bson"))
println()

# ---------------------------------------------------------------------------
# Sweep 1: joint M tightening to identify Strang slope and find where it
# stops dominating. At low M the Strang error swamps any quadrature error.
# ---------------------------------------------------------------------------
println("=== Sweep 1: joint M ramp (r_D=11, r_bm=14, r_bp=18 — tight grid) ===")
println("Identifies Strang slope and where it drops below quadrature.")
println()
@printf("  %-8s  %-12s  %-12s\n", "M_user", "‖L · σ_β‖_HS", "ratio")
let prev = NaN
    for M in (1, 2, 4, 8, 16, 32, 64, 128, 256)
        res = residue_at(r_D=11, r_bm=14, r_bp=18, M_D=M, M_bm=M, M_bp=M)
        ratio = isnan(prev) ? "—" : @sprintf("%.2f", prev/res)
        @printf("  %-8d  %.3e   %s\n", M, res, ratio)
        prev = res
    end
end
println()

# ---------------------------------------------------------------------------
# Sweep 2: ASYMMETRIC M — find which leg drives the Strang error.
# Hold two legs at large M and sweep the third.
# ---------------------------------------------------------------------------
println("=== Sweep 2: asymmetric M (hold two legs at M=128, sweep the third) ===")
println()
println("Tighten M_D only (M_bm = M_bp = 128):")
@printf("  %-8s  %-12s\n", "M_D", "‖L · σ_β‖_HS")
for M_D in (1, 4, 16, 64, 128, 256)
    res = residue_at(r_D=11, r_bm=14, r_bp=18, M_D=M_D, M_bm=128, M_bp=128)
    @printf("  %-8d  %.3e\n", M_D, res)
end
println()
println("Tighten M_b_minus only (M_D = M_bp = 128):")
@printf("  %-8s  %-12s\n", "M_bm", "‖L · σ_β‖_HS")
for M_bm in (1, 4, 16, 64, 128, 256)
    res = residue_at(r_D=11, r_bm=14, r_bp=18, M_D=128, M_bm=M_bm, M_bp=128)
    @printf("  %-8d  %.3e\n", M_bm, res)
end
println()
println("Tighten M_b_plus only (M_D = M_bm = 128):")
@printf("  %-8s  %-12s\n", "M_bp", "‖L · σ_β‖_HS")
for M_bp in (1, 4, 16, 64, 128, 256)
    res = residue_at(r_D=11, r_bm=14, r_bp=18, M_D=128, M_bm=128, M_bp=M_bp)
    @printf("  %-8d  %.3e\n", M_bp, res)
end
println()

# ---------------------------------------------------------------------------
# Sweep 3: at very large M (Strang error << 1e-6), sweep r_D, r_bm, r_bp to
# expose the QUADRATURE floor in each leg.
# ---------------------------------------------------------------------------
println("=== Sweep 3: quadrature floor (M=128 per leg) ===")
println("Joint M=128 puts Strang error ≈ 0; remaining error is quadrature.")
println()
println("r_b_plus sweep at r_D=11, r_b_minus=14, M=128:")
@printf("  %-8s  %-12s\n", "r_bp", "‖L · σ_β‖_HS")
for r_bp in (10, 12, 14, 16, 18, 20)
    res = residue_at(r_D=11, r_bm=14, r_bp=r_bp, M_D=128, M_bm=128, M_bp=128)
    @printf("  %-8d  %.3e\n", r_bp, res)
end
println()
println("r_D sweep at r_b_minus=14, r_b_plus=18, M=128:")
@printf("  %-8s  %-12s\n", "r_D", "‖L · σ_β‖_HS")
for r_D in (7, 9, 11, 13, 15)
    res = residue_at(r_D=r_D, r_bm=14, r_bp=18, M_D=128, M_bm=128, M_bp=128)
    @printf("  %-8d  %.3e\n", r_D, res)
end
println()

# ---------------------------------------------------------------------------
# Final tight recipe demonstration.
# ---------------------------------------------------------------------------
println("=== Final tight recipe at ≤ 1e-6 ===")
recipes = [
    (label="(11, 14, 18) M=(128, 128, 128)", r_D=11, r_bm=14, r_bp=18, M_D=128, M_bm=128, M_bp=128),
    (label="(11, 14, 18) M=(256, 64, 64)",   r_D=11, r_bm=14, r_bp=18, M_D=256, M_bm=64,  M_bp=64),
    (label="(11, 14, 18) M=(64, 256, 256)",  r_D=11, r_bm=14, r_bp=18, M_D=64,  M_bm=256, M_bp=256),
    (label="(13, 14, 18) M=(256, 128, 64)",  r_D=13, r_bm=14, r_bp=18, M_D=256, M_bm=128, M_bp=64),
]
println()
println("Each recipe sized to drive Strang << 1e-6 per leg.")
@printf("  %-40s  %-12s\n", "Recipe", "‖L · σ_β‖_HS")
for r in recipes
    res = residue_at(; r_D=r.r_D, r_bm=r.r_bm, r_bp=r.r_bp,
                       M_D=r.M_D, M_bm=r.M_bm, M_bp=r.M_bp)
    @printf("  %-40s  %.3e\n", r.label, res)
end
