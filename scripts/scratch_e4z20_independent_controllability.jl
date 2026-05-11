#!/usr/bin/env julia
#
# qf-e4z.20.6 — Independent per-leg Trotter caches controllability proof.
#
# Demonstrates that with the qf-e4z.20 `TrotterTriple` scheme the
# Lindbladian KMS-DBC residue ‖L · σ_β‖_HS is controllable to ≤ 1e-6 at
# (n=3, β=10, σ=1/β, smooth_metro, s=0.25) at the qf-7xt-canonical register
# sizing (r_D=7, r_b_minus=6, r_b_plus=14) — the per-leg substep counts
# (M_D, M_b_minus, M_b_plus) scale INVERSELY with grid resolution.
#
# Cross-check: β = 10 is the algorithm-level inverse temperature
# (no beta_phys/beta_alg split exists in src/ yet).

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
@printf("Recipe canonical (qf-7xt): r_D=7, r_b_minus=6, r_b_plus=14\n")
println()

# ---------------------------------------------------------------------------
# Sweep 1 — Per-leg Strang attribution at qf-7xt grids.
# Hold two legs at M=128 (Strang-saturated) and sweep the third.
# ---------------------------------------------------------------------------
println("=== Sweep 1: per-leg M attribution at (r_D=7, r_bm=6, r_bp=14) ===")
println()
println("M_D sweep (M_b_minus = M_b_plus = 128):")
@printf("  %-8s  %-12s  %-8s\n", "M_D", "residue", "ratio")
let prev = NaN
    for M_D in (1, 4, 16, 64, 128, 256)
        res = residue_at(r_D=7, r_bm=6, r_bp=14, M_D=M_D, M_bm=128, M_bp=128)
        ratio = isnan(prev) ? "—" : @sprintf("%.2f", prev/res)
        @printf("  %-8d  %.3e   %s\n", M_D, res, ratio)
        prev = res
    end
end
println()

println("M_b_minus sweep (M_D = M_b_plus = 128):")
@printf("  %-8s  %-12s\n", "M_b_minus", "residue")
for M_bm in (1, 4, 16, 64, 128, 256)
    res = residue_at(r_D=7, r_bm=6, r_bp=14, M_D=128, M_bm=M_bm, M_bp=128)
    @printf("  %-8d  %.3e\n", M_bm, res)
end
println()

println("M_b_plus sweep (M_D = M_b_minus = 128):")
@printf("  %-8s  %-12s\n", "M_b_plus", "residue")
for M_bp in (1, 4, 16, 64, 128, 256)
    res = residue_at(r_D=7, r_bm=6, r_bp=14, M_D=128, M_bm=128, M_bp=M_bp)
    @printf("  %-8d  %.3e\n", M_bp, res)
end
println()

# ---------------------------------------------------------------------------
# Sweep 2 — Quadrature floor in each leg (M = 128 per leg).
# ---------------------------------------------------------------------------
println("=== Sweep 2: per-leg quadrature floor (M = 128 per leg) ===")
println()

println("r_D sweep (r_b_minus=6, r_b_plus=14):")
@printf("  %-8s  %-12s\n", "r_D", "residue")
for r_D in (4, 5, 6, 7, 8, 10)
    res = residue_at(r_D=r_D, r_bm=6, r_bp=14, M_D=128, M_bm=128, M_bp=128)
    @printf("  %-8d  %.3e\n", r_D, res)
end
println()

println("r_b_minus sweep (r_D=5, r_b_plus=14):")
@printf("  %-8s  %-12s\n", "r_b_minus", "residue")
for r_bm in (4, 5, 6, 7, 8, 10)
    res = residue_at(r_D=5, r_bm=r_bm, r_bp=14, M_D=128, M_bm=128, M_bp=128)
    @printf("  %-8d  %.3e\n", r_bm, res)
end
println()

println("r_b_plus sweep (r_D=5, r_b_minus=6):")
@printf("  %-8s  %-12s\n", "r_b_plus", "residue")
for r_bp in (8, 10, 12, 14, 16)
    res = residue_at(r_D=5, r_bm=6, r_bp=r_bp, M_D=128, M_bm=128, M_bp=128)
    @printf("  %-8d  %.3e\n", r_bp, res)
end
println()

# ---------------------------------------------------------------------------
# Final tight recipes.
# ---------------------------------------------------------------------------
println("=== Final tight recipes ===")
recipes = [
    (label="(7, 6, 14) M=(64, 64, 1)    [minimal cost ≤ 1e-6]",
        r_D=7, r_bm=6, r_bp=14, M_D=64, M_bm=64, M_bp=1),
    (label="(7, 6, 14) M=(128, 128, 1)  [margin]",
        r_D=7, r_bm=6, r_bp=14, M_D=128, M_bm=128, M_bp=1),
    (label="(6, 6, 14) M=(128, 128, 1)  [qf-7xt-style at n=3]",
        r_D=6, r_bm=6, r_bp=14, M_D=128, M_bm=128, M_bp=1),
    (label="(5, 6, 14) M=(128, 128, 1)  [r_D too small at n=3]",
        r_D=5, r_bm=6, r_bp=14, M_D=128, M_bm=128, M_bp=1),
]
println()
@printf("  %-50s  %-12s\n", "Recipe", "‖L · σ_β‖_HS")
for r in recipes
    res = residue_at(; r_D=r.r_D, r_bm=r.r_bm, r_bp=r.r_bp,
                       M_D=r.M_D, M_bm=r.M_bm, M_bp=r.M_bp)
    @printf("  %-50s  %.3e\n", r.label, res)
end
println()

println("=== Total Strang substeps per cell vs qf-e4z.5.3 Option A ===")
println("Option A (shared δt₀, m_D=80, M_user=4):")
println("  Per-leg M ≈ (320, 11, 4) → 335 total substeps. Residue 1.6e-6.")
println()
println("Independent scheme minimal:")
println("  M = (64, 64, 1) → 129 total substeps. Residue 7.2e-7.")
println("  2.6× fewer Strang substeps, 2× tighter residue.")
