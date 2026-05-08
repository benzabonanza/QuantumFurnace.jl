#!/usr/bin/env julia
#
# qf-7xt: Coherent-term B convergence — B_time → B_bohr at fixed b_± grid range.
#
# Mirror of `scratch_energy_ref_convergence.jl` and `scratch_time_ref_convergence.jl`,
# but for the COHERENT term B (Lamb shift).  In TimeDomain B is built as a
# nested Riemann sum with two INDEPENDENT registers (qf-9z0):
#
#     B = t0_-·t0_+  Σ_t Σ_τ  b_-(t) · b_+(τ) · U(jumps; t, τ)
#
# Each register has its own (r, w0, t0) triple.  The truncation
# `_compute_truncated_func(_compute_b_*, time_labels, ...; atol=1e-12)` chops
# the time labels where the function magnitude is below 1e-12 — so the grid
# range only matters insofar as it must comfortably cover the SUPPORT of the
# function (where |b(t)| ≥ atol).
#
# Empirical supports at β=10, σ=0.1 (filter-independent for b_-, filter-
# dependent for b_+):
#   b_-(t)         : |t| ≤ 5.55       (set by 1/cosh(2π t / (βσ)) decay)
#   b_+(t) gaussian: |t| ≤ 2.50       (exp(-4t²) at our params)
#   b_+(t) kinky   : |t| ≤ 3.25       (exp(-2t²) / (t·(2t+i)))
#   b_+(t) smooth  : |t| ≤ 2.90       (smoothed Metro)
#
# We pick w0_b_minus, w0_b_plus such that the GRID range = 2π/w0 is ≥ ~3×
# support (generous safety margin so the truncation has nothing to do beyond
# the actual support).  No new "time-range" constants are introduced —
# `_compute_truncated_func` sets the effective Riemann set.
#
# Method
# ------
# Two SPLIT sweeps to separate r_- and r_+ contributions (per thesis 2_methods
# eq:r_- vs eq:r_+_gaussian / eq:r_+_metro — they're independent knobs):
#   (A) Sweep r_b_minus at fixed r_b_plus = R_LARGE  →  isolates outer slope.
#   (B) Sweep r_b_plus  at fixed r_b_minus = R_LARGE →  isolates inner slope.
#
# At each r:
#   t0_-(r) = T_RANGE_MINUS / 2^r        (only varies in sweep A)
#   t0_+(r) = T_RANGE_PLUS  / 2^r        (only varies in sweep B)
#   B_time  = nested Riemann sum (calls _compute_truncated_func internally)
#   B_bohr  = analytical reference
#   err(r)  = opnorm(B_time(r, R_LARGE) − B_bohr)   (or with r_-, r_+ swapped)
#
# Predictions (thesis eq:r_-, eq:r_+_*):
#   b_-(t) is smooth (cosh-decay convolved with exp·sin) → super-algebraic
#     for ANY filter (the outer kernel is filter-independent).
#   b_+(t) is filter-dependent:
#     - gaussian: smooth Gaussian → super-algebraic in t0_+.
#     - kinky / smooth Metro: 1/t · (2t+i) envelope + η-step regularisation
#       → algebraic slope -1 in 1/N_+ (per thesis eq:r_+_metro and the t=0
#       Cauchy P.V. anomaly noted in `.claude-memory/trap_rule_t0_lhopital_origin.md`).

using Printf
using LinearAlgebra
using Random
using QuantumFurnace
using QuantumFurnace: construct_lindbladian, _load_hamiltonian_bson,
                      _compute_b_minus, _compute_b_plus, _compute_b_plus_metro,
                      _compute_truncated_func, B_time, B_bohr

Random.seed!(20260506)

# ── Setup ────────────────────────────────────────────────────────────────────
const N_QUBITS  = 4
const BETA      = 10.0
const SIGMA     = 1.0 / BETA
const ETA       = 1e-3
const FILTERS   = (:gaussian, :smooth, :kinky)
const R_GRID    = 4:1:15        # 2^r ∈ {16, …, 32768}
const R_LARGE   = 16            # "the other register is converged" baseline

# Grid ranges chosen ≥ 3× empirical b_± support so `_compute_truncated_func`
# (atol=1e-12) sets the effective range.  These yield w0_b_*  via 2π/range.
const T_RANGE_MINUS = 18.0      # 3.2× the |t| ≤ 5.55 support of b_-
const T_RANGE_PLUS  = 12.0      # 3.7× the worst-case |t| ≤ 3.25 (kinky) of b_+
const W0_MINUS  = 2π / T_RANGE_MINUS
const W0_PLUS   = 2π / T_RANGE_PLUS

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

@info "Fixture / window choice" N_QUBITS D BETA SIGMA T_RANGE_MINUS T_RANGE_PLUS W0_MINUS W0_PLUS

# ── Build B_bohr per filter (one-time analytical reference) ────────────────
function build_bohr_cfg(filter::Symbol)
    common = (
        sim = Lindbladian(), domain = BohrDomain(), construction = KMS(),
        num_qubits = N_QUBITS, beta = BETA, sigma = SIGMA,
    )
    if filter === :gaussian
        sigma_gamma = SIGMA
        w_gamma = BETA * (SIGMA^2 + sigma_gamma^2) / 2
        return Config(; common..., with_linear_combination = false,
            gaussian_parameters = (w_gamma, sigma_gamma))
    elseif filter === :kinky
        return Config(; common..., with_linear_combination = true,
            a = 0.0, s = 0.0, eta = ETA)
    elseif filter === :smooth
        return Config(; common..., with_linear_combination = true,
            a = 0.0, s = 0.25, eta = ETA)
    end
end

# Build B_time directly via the underlying API (avoids needing a TimeDomain
# Config with all its dissipative parameters; we only need the b_± dicts).
function build_time_b(filter::Symbol, r_minus::Int, r_plus::Int)
    N_minus = 2^r_minus;  t0_minus = T_RANGE_MINUS / N_minus
    N_plus  = 2^r_plus;   t0_plus  = T_RANGE_PLUS  / N_plus
    grid_minus = collect((-(N_minus ÷ 2)):((N_minus ÷ 2) - 1)) .* t0_minus
    grid_plus  = collect((-(N_plus  ÷ 2)):((N_plus  ÷ 2) - 1)) .* t0_plus
    b_minus = _compute_truncated_func(_compute_b_minus, grid_minus, BETA, SIGMA)
    b_plus = if filter === :gaussian
        sigma_gamma = SIGMA
        w_gamma = BETA * (SIGMA^2 + sigma_gamma^2) / 2
        _compute_truncated_func(_compute_b_plus, grid_plus, BETA, w_gamma, sigma_gamma)
    elseif filter === :kinky
        _compute_truncated_func(_compute_b_plus_metro, grid_plus, BETA, SIGMA, ETA, 0.0)
    elseif filter === :smooth
        _compute_truncated_func(_compute_b_plus_metro, grid_plus, BETA, SIGMA, ETA, 0.25)
    end
    return b_minus, b_plus, t0_minus, t0_plus
end

bohr_cache = Dict{Symbol, Matrix{ComplexF64}}()
println("\nBuilding B_bohr (analytical reference) per filter …")
for filt in FILTERS
    cfg_b = build_bohr_cfg(filt)
    bohr_cache[filt] = B_bohr(ham, jumps, cfg_b)
    @printf "  %s: ‖B_bohr‖_op = %.6e (#truncated b_-: …, #b_+: …)\n" string(filt) opnorm(bohr_cache[filt])
end

# ── Sweep helpers ───────────────────────────────────────────────────────────
function err_at(filt::Symbol, r_minus::Int, r_plus::Int, Bref)
    b_minus, b_plus, t0_minus, t0_plus = build_time_b(filt, r_minus, r_plus)
    Bt = B_time(jumps, ham, b_minus, b_plus, t0_minus, t0_plus, BETA, SIGMA)
    return opnorm(Bt .- Bref), length(b_minus), length(b_plus), t0_minus, t0_plus
end

results_minus = Dict{Tuple{Symbol, Int}, NamedTuple}()
results_plus  = Dict{Tuple{Symbol, Int}, NamedTuple}()

# (A) Sweep r_- at fixed large r_+ — isolates outer (b_-) slope.
println("\n=== Split A: sweep r_b_- at fixed r_b_+ = $R_LARGE ===")
@printf "%-9s %4s %12s %8s %14s\n" "filter" "r_-" "t0_-" "#b_-" "‖ΔB‖_op"
for filt in FILTERS
    Bref = bohr_cache[filt]
    for r in R_GRID
        err, n_m, n_p, t0_m, t0_p = err_at(filt, r, R_LARGE, Bref)
        results_minus[(filt, r)] = (; r, err, t0_minus = t0_m, nb_minus = n_m)
        @printf "%-9s %4d %12.6e %8d %14.6e\n" string(filt) r t0_m n_m err
    end
    println()
end

# (B) Sweep r_+ at fixed large r_- — isolates inner (b_+) slope.
println("\n=== Split B: sweep r_b_+ at fixed r_b_- = $R_LARGE ===")
@printf "%-9s %4s %12s %8s %14s\n" "filter" "r_+" "t0_+" "#b_+" "‖ΔB‖_op"
for filt in FILTERS
    Bref = bohr_cache[filt]
    for r in R_GRID
        err, n_m, n_p, t0_m, t0_p = err_at(filt, R_LARGE, r, Bref)
        results_plus[(filt, r)] = (; r, err, t0_plus = t0_p, nb_plus = n_p)
        @printf "%-9s %4d %12.6e %8d %14.6e\n" string(filt) r t0_p n_p err
    end
    println()
end

# ── Slope analysis ──────────────────────────────────────────────────────────
function slope_per_r(rs, errs)
    use = findall(e -> e > 1e-13 && isfinite(e), errs)
    length(use) < 3 && return nothing
    rs_u = rs[use]; le = log10.(errs[use])
    mx = sum(rs_u) / length(use); my = sum(le) / length(use)
    cov = sum((rs_u .- mx) .* (le .- my)) / length(use)
    var = sum((rs_u .- mx) .^ 2) / length(use)
    return (slope = cov / var, rs = rs_u, errs = errs[use])
end

println("\n=== Slope analysis ===")
for filt in FILTERS
    rs = collect(R_GRID)
    em = [results_minus[(filt, r)].err for r in rs]
    ep = [results_plus[(filt, r)].err for r in rs]
    sm = slope_per_r(rs, em)
    sp = slope_per_r(rs, ep)
    @printf "  %-9s  outer (r_-): %s   inner (r_+): %s\n" string(filt) (
        sm === nothing ? "saturated" : @sprintf("slope %+7.3f / r  errs %s", sm.slope, [round(e, sigdigits=3) for e in sm.errs])
    ) (
        sp === nothing ? "saturated" : @sprintf("slope %+7.3f / r  errs %s", sp.slope, [round(e, sigdigits=3) for e in sp.errs])
    )
end
