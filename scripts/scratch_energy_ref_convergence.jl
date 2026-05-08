#!/usr/bin/env julia
#
# qf-7xt.3: EnergyDomain → BohrDomain convergence at fixed ω-range.
#
# Question
# --------
# Fix the ω-range to the integrand support 2·(‖H‖ + TAIL_C·σ).  Sweep r — and
# implicitly w0 = range/2^r — and measure how well the EnergyDomain Riemann
# sum approximates the analytical BohrDomain integral, separately for the
# three filters: Gaussian, kinky Metropolis (s=0), smooth Metropolis (s=0.25).
#
# Method
# ------
# At each r ∈ R_GRID:
#   w0(r)        = omega_range / 2^r          (window-conjugate spacing)
#   cfg_bohr(r)  = Config(BohrDomain,   r, w0(r), filter)
#   cfg_energy(r)= Config(EnergyDomain, r, w0(r), filter)
#   L_b(r) = construct_lindbladian(cfg_bohr,   ham, jumps; include_coherent=false)
#   L_e(r) = construct_lindbladian(cfg_energy, ham, jumps; include_coherent=false)
#   err(r) = opnorm(L_e(r) − L_b(r))                        (raw difference)
#
# Both Lindbladians are built dissipator-only (`include_coherent=false`) and
# at the SAME (r, w0).  Post-qf-etx.2, `gamma_norm_factor` is grid-independent
# (= 1.0 for every standard family) so no `/gnf_*` rescaling is needed.
#
# Unified ω-range principle
# -------------------------
# The integrand `f(ω−ν)·γ(ω)·f(ω−ν')` decays as Gaussian f² (kinky γ saturates
# at 1, doesn't help), so it falls below ε at |ω − ν| ≥ σ·√(−2 ln ε).  For
# (ν, ν') ∈ [-‖H‖, ‖H‖]², the support spans 2·(‖H‖ + σ·√(−2 ln ε)).  ONE knob
# — a tail constant TAIL_C — fixes the grid extent (and, equivalently, the
# integrand cutoff):
#     omega_range = 2·(‖H‖ + TAIL_C·σ)         ⇔   ε = exp(−TAIL_C²/2).
# At TAIL_C = 8 → ε ≈ 1.3e-14 (machine precision); the grid IS the integrand
# support and the library's `_truncate_energy_labels` (cutoff 1e-12, c≈7.43)
# becomes a near-no-op.
#
# Predictions
# -----------
# Riemann-sum convergence in 1/N (here N = 2^r at fixed range):
#   - Gaussian γ_G:        super-algebraic / exponential in r — machine-precision floor.
#   - Smooth Metro γ_M^s:  super-algebraic, slower than Gaussian (Gevrey-1/2).
#   - Kinky Metro γ_M^0:   algebraic, slope -2 in 1/N (γ has a kink at ω = -βσ²/2).

using Printf
using LinearAlgebra
using Random
using QuantumFurnace
using QuantumFurnace: construct_lindbladian, _load_hamiltonian_bson

Random.seed!(20260506)

# ── Setup ────────────────────────────────────────────────────────────────────
const N_QUBITS = 4
const BETA     = 10.0
const SIGMA    = 1.0 / BETA
const TAIL_C   = 8.0           # f² < exp(-TAIL_C²/2) ≈ 1.3e-14 outside
const ETA      = 1e-3
const FILTERS  = (:gaussian, :smooth, :kinky)
const R_GRID   = 3:1:14        # 2^r ∈ {8, ..., 16384}

ham_path = "/Users/bence/code/QuantumFurnace.jl/hamiltonians/heis_xxx_zzdisordered_periodic_n$(N_QUBITS).bson"
ham = _load_hamiltonian_bson(ham_path, BETA)

jump_paulis = [[X], [Y], [Z]]
num_jumps_total = length(jump_paulis) * N_QUBITS
jump_norm = sqrt(num_jumps_total)
jumps = JumpOp[]
for pauli in jump_paulis
    for site in 1:N_QUBITS
        op = Matrix(pad_term(pauli, N_QUBITS, site)) ./ jump_norm
        op_eb = ham.eigvecs' * op * ham.eigvecs
        push!(jumps, JumpOp(op, op_eb, op == transpose(op), op == op'))
    end
end
const D = size(ham.data, 1)

H_norm = opnorm(ham.data)
omega_range = 2 * (H_norm + TAIL_C * SIGMA)
@info "Fixture" N_QUBITS D BETA SIGMA H_norm TAIL_C omega_range cutoff_ε = exp(-TAIL_C^2 / 2)

# ── Config builder ───────────────────────────────────────────────────────────
function build_cfg(domain, filter::Symbol, r::Int, w0::Real)
    t0 = 2π / (2^r * w0)
    common = (
        sim = Lindbladian(), domain = domain, construction = KMS(),
        num_qubits = N_QUBITS, beta = BETA, sigma = SIGMA,
        num_energy_bits = r, w0 = w0, t0 = t0,
    )
    if filter === :gaussian
        sigma_gamma = SIGMA
        w_gamma = BETA * (SIGMA^2 + sigma_gamma^2) / 2
        return Config(; common...,
            with_linear_combination = false,
            gaussian_parameters = (w_gamma, sigma_gamma),
        )
    elseif filter === :kinky
        return Config(; common...,
            with_linear_combination = true,
            a = 0.0, s = 0.0, eta = ETA,
        )
    elseif filter === :smooth
        return Config(; common...,
            with_linear_combination = true,
            a = 0.0, s = 0.25, eta = ETA,
        )
    end
end

# ── Convergence sweep ────────────────────────────────────────────────────────
results = Dict{Tuple{Symbol, Int}, NamedTuple}()
println("\n=== EnergyDomain → BohrDomain convergence ===")
@printf "%-9s %4s %12s %14s %14s\n" "filter" "r" "w0" "‖ΔL‖" "‖L_b‖"
for filt in FILTERS
    for r in R_GRID
        w0 = omega_range / 2^r
        cfg_b = build_cfg(BohrDomain(),   filt, r, w0)
        cfg_e = build_cfg(EnergyDomain(), filt, r, w0)
        L_b = construct_lindbladian(jumps, cfg_b, ham; include_coherent = false)
        L_e = construct_lindbladian(jumps, cfg_e, ham; include_coherent = false)
        err = opnorm(L_e .- L_b)
        nrm_b = opnorm(L_b)
        results[(filt, r)] = (; r, w0, err, nrm_b)
        @printf "%-9s %4d %12.6e %14.6e %14.6e\n" string(filt) r w0 err nrm_b
    end
    println()
end

# ── Slope analysis (log10 err vs r) ─────────────────────────────────────────
println("\n=== Slope analysis (log₁₀ ‖ΔL‖ per r, range above noise floor) ===")
for filt in FILTERS
    rs = collect(R_GRID)
    errs = [results[(filt, r)].err for r in rs]
    use = findall(e -> e > 1e-13 && isfinite(e), errs)
    length(use) < 3 && (println("  $filt — saturated at noise floor; skipping."); continue)
    rs_used = rs[use]
    log_e = log10.(errs[use])
    n_pts = length(use)
    mean_x = sum(rs_used) / n_pts
    mean_y = sum(log_e) / n_pts
    cov = sum((rs_used .- mean_x) .* (log_e .- mean_y)) / n_pts
    var = sum((rs_used .- mean_x) .^ 2) / n_pts
    slope = cov / var
    @printf "  %-9s slope log₁₀(‖ΔL‖)/r = %+7.3f   (errs at r∈%s: %s)\n" string(filt) slope rs_used [round(e, sigdigits = 3) for e in errs[use]]
end
