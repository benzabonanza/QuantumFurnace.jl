#!/usr/bin/env julia
# ============================================================================
# Delta-Scaling Validation & Richardson Extrapolation Script
# ============================================================================
#
# Tests whether spectral gap estimation error scales as O(delta) with Trotter
# step size, then applies Richardson extrapolation to cancel the leading error.
#
# Motivation (Quick Tasks 22-29):
#   - Quick-22: Fixed delta_eff double-counting; residual factor ~1.6x (n=4)
#     is discrete-step Kraus effect at delta=0.01.
#   - Quick-28: Disordered Heisenberg shows ~48.7% (n=4) gap error at delta=0.01.
#   - Quick-29: H alone gives gap/exact = 1.21x (n=4) -- Kraus overestimation.
#
# Hypothesis: If error = C*delta (linear in Trotter step size), then:
#   1. error/delta should be approximately constant across delta values.
#   2. Richardson extrapolation: gap_rich = (h2*gap(h1) - h1*gap(h2)) / (h2-h1)
#      should cancel the O(delta) term and yield much more accurate estimates.
#
# This script uses a pre-generated disordered Hamiltonian loaded via
# `load_hamiltonian("heis", n; beta=10.0)` from BSON files.
#
# Parameters (locked decisions):
#   - beta = 10.0
#   - System size: n = 4 only
#   - Trajectories: 20,000
#   - Delta values: [0.1, 0.01, 0.001]
#   - Domain: TimeDomain with with_coherent=false
#   - Initial state: psi0[end] = 1.0 (excited state)
#   - Seed: 42 (same for all runs)
#
# Usage:
#   cd QuantumFurnace.jl && julia --project experiments/validate_gap_delta_scaling.jl
# ============================================================================

using QuantumFurnace
using LinearAlgebra
using Printf

# ============================================================================
# Section 0: Constants and Setup
# ============================================================================

const BETA = 10.0
const NTRAJ = 20_000
const SAVE_EVERY = 10
const SEED = 42
const DELTAS = [0.1, 0.01, 0.001]

# Grid parameters matching test suite conventions
const NUM_ENERGY_BITS = 12
const W0 = 0.05
const T0 = 2pi / (2^NUM_ENERGY_BITS * W0)
const NUM_TROTTER_STEPS_PER_T0 = 10
const SIGMA = 1.0 / BETA  # 0.1

# Pauli matrices
const X = ComplexF64[0 1; 1 0]
const Y = ComplexF64[0 -im; im 0]
const Z = ComplexF64[1 0; 0 -1]

# ---------------------------------------------------------------------------
# Helper: create disordered Heisenberg system via load_hamiltonian
# ---------------------------------------------------------------------------
"""
    make_system(n, beta) -> (ham, jumps, dim)

Load a pre-generated disordered periodic Heisenberg Hamiltonian and build
single-site Pauli jump operators in the Hamiltonian eigenbasis.
"""
function make_system(n, beta)
    ham = load_hamiltonian("heis", n; beta=beta)
    dim = 2^n

    # Jump operators: single-site Paulis (X, Y, Z) on each site
    # Constructed via pad_term, then transformed to Hamiltonian eigenbasis
    jump_paulis = [[X], [Y], [Z]]
    num_of_jumps = 3 * n
    jump_normalization = sqrt(num_of_jumps)
    V = ham.eigvecs  # Hamiltonian eigenbasis (TimeDomain)

    jumps = JumpOp[]
    for pauli in jump_paulis
        for site in 1:n
            jump_op = Matrix(pad_term(pauli, n, site)) ./ jump_normalization
            jump_in_eigen = V' * jump_op * V
            orthogonal = (jump_op == transpose(jump_op))
            herm = (jump_op == jump_op')
            push!(jumps, JumpOp(jump_op, jump_in_eigen, orthogonal, herm))
        end
    end

    return ham, jumps, dim
end

# ---------------------------------------------------------------------------
# Helper: create LiouvConfig for TimeDomain
# ---------------------------------------------------------------------------
function make_liouv_config(n)
    LiouvConfig(;
        num_qubits = n,
        with_coherent = false,
        with_linear_combination = true,
        domain = TimeDomain(),
        beta = BETA,
        sigma = SIGMA,
        a = BETA / 30.0,
        b = 0.4,
        num_energy_bits = NUM_ENERGY_BITS,
        w0 = W0,
        t0 = T0,
        num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
    )
end

# ---------------------------------------------------------------------------
# Helper: create ThermalizeConfig for TimeDomain
# ---------------------------------------------------------------------------
function make_thermalize_config(n; mixing_time=50.0, delta=0.01)
    ThermalizeConfig(;
        num_qubits = n,
        with_coherent = false,
        with_linear_combination = true,
        domain = TimeDomain(),
        beta = BETA,
        sigma = SIGMA,
        a = BETA / 30.0,
        b = 0.4,
        num_energy_bits = NUM_ENERGY_BITS,
        w0 = W0,
        t0 = T0,
        num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
        mixing_time = mixing_time,
        delta = delta,
    )
end

println("=" ^ 70)
println("  TROTTER DELTA-SCALING VALIDATION & RICHARDSON EXTRAPOLATION")
println("  beta=$BETA, ntraj=$NTRAJ, seed=$SEED")
println("  Deltas: $DELTAS")
println("  Testing: Does gap estimation error scale as O(delta)?")
println("=" ^ 70)
println()

# ============================================================================
# Section 1: n=4 system setup and exact gap
# ============================================================================

const n = 4

@printf("Loading disordered Hamiltonian for n=%d...\n", n)
ham, jumps, dim = make_system(n, BETA)

config_l = make_liouv_config(n)
@printf("Building Lindbladian for n=%d...\n", n)
liouv_result = run_lindbladian(jumps, config_l, ham)

exact_gap = abs(real(liouv_result.spectral_gap))
@printf("\nExact gap (ARPACK): %.10f\n\n", exact_gap)

# Build observables ONCE, reuse for all 3 delta runs
obs, obs_names = build_preset_trajectory_observables(ham, n)
@printf("Observables built: %s\n", join(obs_names, ", "))

# Initial state: excited state (locked decision)
psi0 = zeros(ComplexF64, dim)
psi0[end] = 1.0

# Mixing time: generous value (5 decay times, minimum 10.0)
mixing_time = max(5.0 / exact_gap, 10.0)
@printf("Mixing time: %.1f (5/gap=%.1f)\n\n", mixing_time, 5.0 / exact_gap)

# Base config (delta=0.01 default, but we override per-run)
config_t = make_thermalize_config(n; mixing_time=mixing_time)

# ============================================================================
# Section 2: Run trajectory gap estimation for each delta
# ============================================================================

# Store results: delta => SpectralGapResult
gap_results = Dict{Float64, Any}()

for delta_val in DELTAS
    println("=" ^ 70)
    @printf("  Delta = %.4f -- Running %d trajectories\n", delta_val, NTRAJ)
    println("=" ^ 70)

    t_start = time()
    gap_result = estimate_spectral_gap(jumps, config_t, psi0, ham;
        observables=obs, observable_names=obs_names,
        ntraj=NTRAJ, save_every=SAVE_EVERY, seed=SEED,
        skip_initial=0.1, delta=delta_val)
    t_elapsed = time() - t_start
    @printf("Completed in %.1f seconds.\n\n", t_elapsed)

    gap_results[delta_val] = gap_result

    # Per-observable table
    @printf("Per-observable fits (delta=%.4f):\n", delta_val)
    @printf("%-12s  %12s  %8s  %10s  %s\n",
            "Observable", "Gap", "R2", "Converged", "CI")
    @printf("%-12s  %12s  %8s  %10s  %s\n",
            "-" ^ 12, "-" ^ 12, "-" ^ 8, "-" ^ 10, "-" ^ 20)
    for (i, name) in enumerate(gap_result.observable_names)
        fit = gap_result.per_observable[i]
        @printf("%-12s  %12.8f  %8.4f  %10s  [%.4f, %.4f]\n",
                name, fit.gap, fit.r_squared,
                fit.converged ? "yes" : "no",
                fit.gap_ci[1], fit.gap_ci[2])
    end
    println()

    error = gap_result.gap - exact_gap
    abs_error = abs(error)
    @printf("  Best observable: %s\n", gap_result.best_observable)
    @printf("  Estimated gap:   %.8f\n", gap_result.gap)
    @printf("  Exact gap:       %.8f\n", exact_gap)
    @printf("  Error:           %+.8f\n", error)
    @printf("  |Error|:         %.8f\n", abs_error)
    @printf("  Error/delta:     %+.6f\n", error / delta_val)
    @printf("  |Error|/delta:   %.6f\n", abs_error / delta_val)
    println()
end

# ============================================================================
# Section 3: Delta-scaling analysis
# ============================================================================

println()
println("=" ^ 70)
println("  SECTION 3: DELTA-SCALING ANALYSIS (best observable per run)")
println("=" ^ 70)
println()

@printf("%-10s  %14s  %14s  %14s  %14s\n",
        "delta", "est_gap", "exact_gap", "error", "error/delta")
@printf("%-10s  %14s  %14s  %14s  %14s\n",
        "-" ^ 10, "-" ^ 14, "-" ^ 14, "-" ^ 14, "-" ^ 14)

error_over_delta = Float64[]
for delta_val in DELTAS
    gr = gap_results[delta_val]
    err = gr.gap - exact_gap
    ratio = err / delta_val
    push!(error_over_delta, ratio)
    @printf("%-10.4f  %14.8f  %14.8f  %+14.8f  %+14.6f\n",
            delta_val, gr.gap, exact_gap, err, ratio)
end
println()

# Check if error/delta is approximately constant (within factor of 2)
abs_ratios = abs.(error_over_delta)
if length(abs_ratios) >= 2 && minimum(abs_ratios) > 0
    max_ratio = maximum(abs_ratios) / minimum(abs_ratios)
    o_delta_confirmed = max_ratio < 2.0
    @printf("Ratio spread: max/min = %.2f (threshold: < 2.0)\n", max_ratio)
    if o_delta_confirmed
        println("CONCLUSION: O(delta) scaling CONFIRMED")
    else
        println("CONCLUSION: O(delta) scaling NOT confirmed (ratio spread too large)")
    end
else
    o_delta_confirmed = false
    println("CONCLUSION: Cannot assess O(delta) scaling (insufficient data or zero ratios)")
end
println()

# --- Per-observable delta-scaling analysis ---
println("=" ^ 70)
println("  PER-OBSERVABLE DELTA-SCALING ANALYSIS")
println("=" ^ 70)
println()

# Find the "fixed" reference observable: best at delta=0.01
ref_result = gap_results[0.01]
fixed_obs_name = ref_result.best_observable
@printf("Fixed reference observable (best at delta=0.01): %s\n\n", fixed_obs_name)

# For each observable, check error/delta across deltas
for (obs_idx, name) in enumerate(obs_names)
    @printf("--- %s ---\n", name)
    all_converged = true
    for delta_val in DELTAS
        gr = gap_results[delta_val]
        fit = gr.per_observable[obs_idx]
        if fit.converged && fit.gap > 0
            err = fit.gap - exact_gap
            @printf("  delta=%-8.4f  gap=%-12.8f  error=%+12.8f  error/delta=%+12.6f\n",
                    delta_val, fit.gap, err, err / delta_val)
        else
            @printf("  delta=%-8.4f  (not converged or gap<=0)\n", delta_val)
            all_converged = false
        end
    end
    println()
end

# Fixed observable tracking across all 3 deltas
println("=" ^ 70)
@printf("  FIXED OBSERVABLE (%s) ACROSS ALL DELTAS\n", fixed_obs_name)
println("=" ^ 70)
println()

fixed_obs_idx = findfirst(==(fixed_obs_name), obs_names)
if fixed_obs_idx !== nothing
    @printf("%-10s  %14s  %14s  %14s\n",
            "delta", "gap", "error", "error/delta")
    @printf("%-10s  %14s  %14s  %14s\n",
            "-" ^ 10, "-" ^ 14, "-" ^ 14, "-" ^ 14)
    for delta_val in DELTAS
        gr = gap_results[delta_val]
        fit = gr.per_observable[fixed_obs_idx]
        if fit.converged && fit.gap > 0
            err = fit.gap - exact_gap
            @printf("%-10.4f  %14.8f  %+14.8f  %+14.6f\n",
                    delta_val, fit.gap, err, err / delta_val)
        else
            @printf("%-10.4f  (not converged or gap<=0)\n", delta_val)
        end
    end
    println()
end

# ============================================================================
# Section 4: Richardson extrapolation
# ============================================================================

println()
println("=" ^ 70)
println("  SECTION 4: RICHARDSON EXTRAPOLATION")
println("  Formula: gap_rich = (h2*gap(h1) - h1*gap(h2)) / (h2 - h1)")
println("  where h1 < h2 (h1 = finer delta)")
println("=" ^ 70)
println()

# Richardson pairs: (h1, h2) where h1 < h2
richardson_pairs = [(0.01, 0.1), (0.001, 0.01)]

for (h1, h2) in richardson_pairs
    println("-" ^ 50)
    @printf("Pair: h1=%.4f, h2=%.4f\n", h1, h2)
    println("-" ^ 50)

    gr1 = gap_results[h1]
    gr2 = gap_results[h2]

    # Best-observable Richardson
    gap1 = gr1.gap
    gap2 = gr2.gap
    gap_rich = (h2 * gap1 - h1 * gap2) / (h2 - h1)

    err_h1 = abs(gap1 - exact_gap)
    err_rich = abs(gap_rich - exact_gap)
    rel_err_h1 = err_h1 / exact_gap
    rel_err_rich = err_rich / exact_gap

    @printf("\nBest-observable Richardson:\n")
    @printf("  gap(h1=%.4f): %.8f (best: %s)\n", h1, gap1, gr1.best_observable)
    @printf("  gap(h2=%.4f): %.8f (best: %s)\n", h2, gap2, gr2.best_observable)
    @printf("  gap_rich:      %.8f\n", gap_rich)
    @printf("  exact_gap:     %.8f\n", exact_gap)
    @printf("  error(h1):     %.8f  (rel: %.4f%%)\n", err_h1, rel_err_h1 * 100)
    @printf("  error(rich):   %.8f  (rel: %.4f%%)\n", err_rich, rel_err_rich * 100)
    if err_rich > 0
        @printf("  improvement:   %.1fx\n", err_h1 / err_rich)
    else
        @printf("  improvement:   exact (zero error)\n")
    end
    println()

    # Per-observable Richardson
    @printf("Per-observable Richardson (h1=%.4f, h2=%.4f):\n", h1, h2)
    @printf("%-12s  %12s  %12s  %12s  %14s  %10s\n",
            "Observable", "gap(h1)", "gap(h2)", "gap_rich", "error_rich", "rel_err%")
    @printf("%-12s  %12s  %12s  %12s  %14s  %10s\n",
            "-" ^ 12, "-" ^ 12, "-" ^ 12, "-" ^ 12, "-" ^ 14, "-" ^ 10)

    for (obs_idx, name) in enumerate(obs_names)
        fit1 = gr1.per_observable[obs_idx]
        fit2 = gr2.per_observable[obs_idx]

        if fit1.converged && fit1.gap > 0 && fit2.converged && fit2.gap > 0
            g1 = fit1.gap
            g2 = fit2.gap
            g_rich = (h2 * g1 - h1 * g2) / (h2 - h1)
            err = abs(g_rich - exact_gap)
            rel = err / exact_gap * 100
            @printf("%-12s  %12.8f  %12.8f  %12.8f  %+14.8f  %10.4f\n",
                    name, g1, g2, g_rich, g_rich - exact_gap, rel)
        else
            @printf("%-12s  (skip -- not converged in both)\n", name)
        end
    end
    println()
end

# ============================================================================
# Section 5: Summary
# ============================================================================

println()
println("=" ^ 70)
println("  SECTION 5: FINAL SUMMARY")
println("=" ^ 70)
println()

@printf("System: n=%d disordered periodic Heisenberg, beta=%.1f\n", n, BETA)
@printf("Exact gap: %.10f\n\n", exact_gap)

# Main results table
@printf("%-10s  %14s  %14s  %10s  %14s\n",
        "delta", "est_gap", "error", "rel_err%", "error/delta")
@printf("%-10s  %14s  %14s  %10s  %14s\n",
        "-" ^ 10, "-" ^ 14, "-" ^ 14, "-" ^ 10, "-" ^ 14)
for delta_val in DELTAS
    gr = gap_results[delta_val]
    err = gr.gap - exact_gap
    rel = abs(err) / exact_gap * 100
    @printf("%-10.4f  %14.8f  %+14.8f  %10.4f  %+14.6f\n",
            delta_val, gr.gap, err, rel, err / delta_val)
end
println()

# Richardson results
println("Richardson extrapolation results:")
for (h1, h2) in richardson_pairs
    gr1 = gap_results[h1]
    gr2 = gap_results[h2]
    gap_rich = (h2 * gr1.gap - h1 * gr2.gap) / (h2 - h1)
    err_rich = abs(gap_rich - exact_gap)
    rel_rich = err_rich / exact_gap * 100
    err_fine = abs(gr1.gap - exact_gap)
    improvement = err_fine > 0 ? err_fine / err_rich : Inf

    @printf("  Pair (%.4f, %.4f): gap_rich=%.8f, error=%.8f (%.4f%%), improvement=%.1fx\n",
            h1, h2, gap_rich, err_rich, rel_rich, improvement)
end
println()

# O(delta) scaling verdict
println("O(delta) scaling assessment:")
if o_delta_confirmed
    println("  CONFIRMED -- error/delta is approximately constant across delta values")
else
    println("  NOT CONFIRMED -- error/delta varies significantly across delta values")
end
println("  error/delta values: ", [@sprintf("%.4f", r) for r in error_over_delta])
println()

# Richardson improvement summary
println("Richardson improvement summary:")
for (h1, h2) in richardson_pairs
    gr1 = gap_results[h1]
    gr2 = gap_results[h2]
    gap_rich = (h2 * gr1.gap - h1 * gr2.gap) / (h2 - h1)
    err_fine = abs(gr1.gap - exact_gap)
    err_rich = abs(gap_rich - exact_gap)
    @printf("  (%.4f, %.4f): fine_error=%.6f -> rich_error=%.6f (%.1fx improvement)\n",
            h1, h2, err_fine, err_rich, err_fine > 0 ? err_fine / err_rich : Inf)
end

println()
println("=" ^ 70)
println("  END OF VALIDATION")
println("=" ^ 70)
