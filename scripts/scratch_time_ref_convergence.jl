#!/usr/bin/env julia
#
# qf-7xt: TimeDomain → BohrDomain convergence at fixed ω-range.
#
# Mirror of `scratch_energy_ref_convergence.jl` for the TimeDomain Lindbladian.
# In TimeDomain the dissipator is built via NUFFT — the t-grid spacing t0 and
# the ω-grid spacing w0 are paired by  t0·w0·2^r = 2π.
#
# Two physical constraints, both must hold:
#   (T1) ω-range = 2^r · w0 large enough to cover the integrand
#        f(ω-ν)·γ(ω)·f(ω-ν') for (ν, ν') ∈ [-‖H‖, ‖H‖]².  Gaussian f² hits ε_ω
#        at |ω-ν| ≥ σ·√(-2 ln ε_ω), so ω-range ≥ 2·(‖H‖ + σ·√(-2 ln ε_ω)).
#   (T2) t-range = 2π/w0 large enough to capture the time kernel f(t):
#        |f(t)| < ε_t at |t| ≥ filter_time_cutoff(σ, ε_t).
#
# Picking ω-range from (T1) sets the smallest allowed `2^r · w0`.  Then for
# any r ≥ R_MIN we have BOTH constraints satisfied IF w0 = ω-range / 2^r
# (which gives t-range = 2π·2^r/ω-range, growing with r above (T2)'s bound).
#
# This sweep mirrors the EnergyDomain sweep exactly:
#   ω-range fixed = 2·(‖H‖ + TAIL_C·σ),    w0(r) = ω-range / 2^r,
#   t0   = 2π / ω-range  (CONSTANT in r),  t-range(r) = 2π · 2^r / ω-range.
# At fixed t0, increasing r adds time samples by extending the time window
# (the time grid NEVER GETS FINER).  The ω-Riemann sum gets finer (w0 → 0).
#
# A separate, complementary sweep — fix w0, sweep r → t0 shrinks — tests how
# fine the t-grid needs to be at fixed ω-resolution; at fixed w0 the err is
# determined by w0 alone (the ω-Riemann sum) and is flat in r once
# t-truncation is below ε_t.  We do that sweep too, second.
#
# Method
# ------
# At each r ∈ R_GRID:
#   w0(r)    = ω-range / 2^r       (varies)
#   t0       = 2π / ω-range         (CONSTANT)
#   cfg_bohr = Config(BohrDomain,   r, w0(r), filter)
#   cfg_time = Config(TimeDomain,   r, w0(r), filter)
#   L_b = construct_lindbladian(cfg_bohr, ham, jumps; include_coherent=false)
#   L_t = construct_lindbladian(cfg_time, ham, jumps; include_coherent=false)
#   err(r) = opnorm(L_t − L_b)         (post-qf-etx.2: gnf is grid-independent)
#
# Predictions
# -----------
# Same as EnergyDomain (Riemann sum in ω):
#   - Gaussian γ_G:        super-algebraic — machine-precision floor.
#   - Smooth Metro γ_M^s:  super-algebraic, slower (Gevrey-1/2).
#   - Kinky Metro γ_M^0:   slope -2 in 1/N (γ has a kink at ω = -βσ²/2).
# Plus a NUFFT precision floor (FINUFFT eps = 1e-12) for the OFT computation.

using Printf
using LinearAlgebra
using Random
using QuantumFurnace
using QuantumFurnace: construct_lindbladian, _load_hamiltonian_bson, filter_time_cutoff

Random.seed!(20260506)

# ── Setup ────────────────────────────────────────────────────────────────────
const N_QUBITS  = 4
const BETA      = 10.0
const SIGMA     = 1.0 / BETA
const TAIL_C    = 8.0          # ω-side: f² < exp(-TAIL_C²/2) outside ‖H‖+TAIL_C·σ
const ETA       = 1e-3
const FILTERS   = (:gaussian, :smooth, :kinky)
const R_GRID    = 3:1:14

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

# ── Fix ω-range from the unified principle ──────────────────────────────────
# Strategy (per user's physics intuition): the harder kernel to resolve here
# is the ENERGY kernel f(ω-ν)·γ(ω)·f(ω-ν') (kinky γ has slope -2 in 1/N).
# The time-OFT kernel is just a Gaussian f(t), super-algebraically resolved
# by a much coarser t-grid.  So pick (r, w0) such that the ω-Riemann sum is
# accurate, and the resulting t0 = 2π/ω-range is automatically fine enough
# for the t-OFT.
const H_NORM     = opnorm(ham.data)
const OMEGA_RANGE = 2 * (H_NORM + TAIL_C * SIGMA)
const T0_FIXED    = 2π / OMEGA_RANGE        # CONSTANT t0 across the sweep
const T_CUTOFF_GAUSS_E14 = filter_time_cutoff(GaussianFilter(SIGMA), 1e-14)

@info "Fixture / window choice" N_QUBITS D BETA SIGMA H_NORM TAIL_C OMEGA_RANGE T0_FIXED T_CUTOFF_GAUSS_E14

# ── Config builder ───────────────────────────────────────────────────────────
function build_cfg(domain, filter::Symbol, r::Int)
    w0 = OMEGA_RANGE / 2^r       # varies with r so 2^r·w0 = OMEGA_RANGE constant
    t0 = 2π / (2^r * w0)         # = T0_FIXED, constant in r
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

# ── Sweep ───────────────────────────────────────────────────────────────────
results = Dict{Tuple{Symbol, Int}, NamedTuple}()
println("\n=== TimeDomain → BohrDomain convergence (fixed ω-range, sweep r → w0) ===")
@printf "%-9s %4s %12s %12s %12s %14s\n" "filter" "r" "w0" "t0" "t-range" "‖ΔL‖"
for filt in FILTERS
    for r in R_GRID
        w0 = OMEGA_RANGE / 2^r
        t0 = 2π / (2^r * w0)        # = T0_FIXED
        t_range_at_r = 2π / w0
        cfg_b = build_cfg(BohrDomain(),  filt, r)
        cfg_t = build_cfg(TimeDomain(),  filt, r)
        L_b = construct_lindbladian(jumps, cfg_b, ham; include_coherent = false)
        L_t = construct_lindbladian(jumps, cfg_t, ham; include_coherent = false)
        # Post-qf-etx.2: gnf is grid-independent, raw L_t-L_b is the Riemann error.
        err = opnorm(L_t .- L_b)
        results[(filt, r)] = (; r, w0, t0, t_range = t_range_at_r, err)
        @printf "%-9s %4d %12.6e %12.6e %12.6e %14.6e\n" string(filt) r w0 t0 t_range_at_r err
    end
    println()
end

# ── Slope analysis (log10 err vs r, above noise floor) ──────────────────────
println("\n=== Slope analysis (log₁₀ ‖ΔL‖ per r) ===")
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
