#!/usr/bin/env julia
#
# qf-f45 spot-check — compute the "true" gap of the CKG smooth_def_s
# Lindbladian at n=6, β_phys=2 by building it in EnergyDomain (where
# B_energy ≡ B_bohr exactly, no coherent quadrature error). Use r_D=8
# to keep dissipator error well below 1e-10. Compare against the
# previously-measured TimeDomain gaps at ε=10⁻³ and "ε=10⁻⁶" register
# sizings to verify the 2.4e-5 TimeDomain coherent floor explains the
# gap shifts.

using Printf
using LinearAlgebra
using BSON
using QuantumFurnace
using QuantumFurnace: apply_lindbladian!, _krylov_spectral_decomposition,
                      default_smooth_s

BLAS.set_num_threads(1)

const N        = 6
const BPHYS    = 2.0
const TAIL_C   = 8.0
const ETA      = 1e-3
const R_D_REF  = 8           # well past the r_D@1e-9 = 7 cutoff
const KRYLOVDIM = 30

ham_path = joinpath(@__DIR__, "..", "hamiltonians",
                    "heis_xxx_zzdisordered_periodic_n$(N).bson")
ham_raw  = BSON.load(ham_path)[:hamiltonian]
ham      = HamHam(ham_raw; beta_phys = BPHYS)
beta_alg = BPHYS * ham.rescaling_factor
sigma    = 1.0 / beta_alg
H_norm   = opnorm(ham.data)
omega_range = 2 * (H_norm + TAIL_C * sigma)
s_used   = default_smooth_s(beta_alg, sigma)

w0_D = omega_range / 2^R_D_REF
t0_D = 2π / (2^R_D_REF * w0_D)

cfg = Config(
    sim = Lindbladian(), domain = EnergyDomain(), construction = KMS(),
    num_qubits = N, beta = beta_alg, beta_phys = BPHYS, sigma = sigma,
    num_energy_bits_D = R_D_REF, w0_D = w0_D, t0_D = t0_D,
    with_linear_combination = true,
    a = 0.0, s = s_used, eta = ETA,
)

jp = [[X], [Y], [Z]]
nrm = sqrt(length(jp) * N)
jumps = JumpOp[]
for pauli in jp, site in 1:N
    op = Matrix(pad_term(pauli, N, site)) ./ nrm
    op_eb = ham.eigvecs' * op * ham.eigvecs
    push!(jumps, JumpOp(op, op_eb, op == transpose(op), op == op'))
end

d = size(ham.data, 1)
println("n=$N  β_phys=$BPHYS  β_alg=$(round(beta_alg, digits=2))  σ=$(round(sigma, digits=5))  s=$(round(s_used, digits=2))  d=$d")
println("EnergyDomain reference (B_energy ≡ B_bohr exact; dissipator at r_D=$R_D_REF, ε_diss < 1e-10)")
flush(stdout)

t_ws = @elapsed ws = Workspace(cfg, ham, jumps)
fwd! = (out, x) -> begin
    apply_lindbladian!(ws, x, cfg, ham)
    copyto!(out, ws.scratch.rho_out)
    return out
end
rho_0 = Matrix{ComplexF64}(I, d, d) ./ d
t_kry = @elapsed decomp = _krylov_spectral_decomposition(
    fwd!, rho_0, d; krylovdim = KRYLOVDIM, tol = 1e-10, sort_mode = :lindbladian)

gap_true = abs(real(decomp.eigenvalues[2]))

println()
@printf "EnergyDomain (TRUE):  gap = %.6e   τ_mix ≈ %.2f  wall = %.2fs (ws=%.2fs, krylov=%.2fs)\n" gap_true (1/gap_true) (t_ws+t_kry) t_ws t_kry
println()
println("Comparison with TimeDomain benches:")
@printf "  TimeDomain ε=1e-6 reg (r_b_plus=17):  gap = 4.0600e-03   shift vs TRUE = %.2f%% (Δgap = %.2e)\n" (abs(4.060e-3 - gap_true)/gap_true * 100) (4.060e-3 - gap_true)
@printf "  TimeDomain ε=1e-3 reg (r_b_plus=7):   gap = 3.3983e-03   shift vs TRUE = %.2f%% (Δgap = %.2e)\n" (abs(3.3983e-3 - gap_true)/gap_true * 100) (3.3983e-3 - gap_true)
println()
println("Coherent-term ‖ΔB‖_op vs B_bohr from the S34 sweep (this cell):")
println("  TimeDomain r_b_plus=17: ‖ΔB‖ ≈ 2.36e-5  (TimeDomain coherent floor at β_phys=2, n=6)")
println("  TimeDomain r_b_plus=7 (extrap): ‖ΔB‖ ≈ 1.15e-3  (well above floor; r_+ too small)")
