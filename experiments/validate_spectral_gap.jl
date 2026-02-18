#!/usr/bin/env julia
# ============================================================================
# Unified Spectral Gap Validation Script
# ============================================================================
#
# Replaces 4 deleted experiment scripts with a single comprehensive validation:
#   1. ARPACK vs dense eigen verification (n=4)
#   2. Eigenbasis overlap analysis (n=4, n=6)
#   3. 20k-trajectory gap estimation (n=4, n=6)
#   4. Pass/fail with diagnostic evidence on failure
#
# Parameters (locked decisions):
#   - beta = 10.0, delta = 0.01
#   - System sizes: n = 4, n = 6
#   - Trajectories: 20,000
#   - Target: n=4 < 1% (k=0 gap, strong overlap); n=6 < 12% (k=pi gap, weak overlap)
#   - Domain: TimeDomain with with_coherent=false
#   - Initial state: psi0[end] = 1.0 (excited state)
#   - Periodic Heisenberg chain
#
# Usage:
#   cd QuantumFurnace.jl && julia --project experiments/validate_spectral_gap.jl
# ============================================================================

using QuantumFurnace
using LinearAlgebra
using Printf

# ============================================================================
# Section 0: Constants and Setup
# ============================================================================

const BETA = 10.0
const DELTA = 0.01
const NTRAJ = 20_000
const SAVE_EVERY = 10
const TARGET_REL_ERROR_N4 = 1e-2   # 1% for n=4 (k=0 gap mode, strong overlap)
const TARGET_REL_ERROR_N6 = 0.12   # 12% for n=6 (k=pi gap, SU(2)-protected; XZ_stagg also has zero overlap)
const SEED = 42

# Grid parameters matching test suite conventions
const NUM_ENERGY_BITS = 12
const W0 = 0.05
const T0 = 2pi / (2^NUM_ENERGY_BITS * W0)
const NUM_TROTTER_STEPS_PER_T0 = 10
const SIGMA = 1.0 / BETA  # 0.1

# ---------------------------------------------------------------------------
# Helper: create periodic Heisenberg chain system
# ---------------------------------------------------------------------------
"""
    make_system(n, beta) -> (ham, jumps, dim)

Create a periodic Heisenberg chain Hamiltonian and single-site Pauli jump operators.
Follows the exact pattern from test/test_helpers.jl (make_test_system).
"""
function make_system(n, beta)
    ham = HamHam([[X,X],[Y,Y],[Z,Z]], [1.0,1.0,1.0], n, beta; periodic=true)
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
function make_thermalize_config(n; mixing_time=50.0)
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
        delta = DELTA,
    )
end

# Tracking pass/fail results
results = Dict{String, Bool}()

println("=" ^ 70)
println("  SPECTRAL GAP VALIDATION")
println("  beta=$BETA, delta=$DELTA, ntraj=$NTRAJ, seed=$SEED")
println("=" ^ 70)
println()

# ============================================================================
# Section 1: ARPACK vs eigen verification (n=4 only)
# ============================================================================

println("=" ^ 70)
println("  Section 1: ARPACK vs eigen verification (n=4)")
println("=" ^ 70)
println()

ham4, jumps4, dim4 = make_system(4, BETA)
config_l4 = make_liouv_config(4)
liouv_result4 = run_lindbladian(jumps4, config_l4, ham4)
L4 = liouv_result4.liouvillian

# ARPACK gap (from run_lindbladian)
gap_arpack = abs(real(liouv_result4.spectral_gap))

# Dense eigen gap
F4 = eigen(Matrix(L4))
sorted_idx4 = sortperm(abs.(real.(F4.values)))
gap_eigen = abs(real(F4.values[sorted_idx4[2]]))

arpack_eigen_diff = abs(gap_arpack - gap_eigen)
arpack_pass = arpack_eigen_diff < 1e-8

@printf("ARPACK gap: %.10f\n", gap_arpack)
@printf("eigen  gap: %.10f\n", gap_eigen)
@printf("Difference: %.3e\n", arpack_eigen_diff)
@printf("Result: %s (threshold: 1e-8)\n", arpack_pass ? "PASS" : "FAIL")
println()

results["ARPACK vs eigen (n=4)"] = arpack_pass

# ============================================================================
# Section 2: Per-system-size analysis (n=4, n=6)
# ============================================================================

for n in [4, 6]
    println("=" ^ 70)
    @printf("  Section 2: System size n=%d (dim=%d, Liouvillian %dx%d)\n",
            n, 2^n, (2^n)^2, (2^n)^2)
    println("=" ^ 70)
    println()

    # --- 2a: System setup ---
    local ham, jumps, dim
    if n == 4
        ham, jumps, dim = ham4, jumps4, dim4
    else
        ham, jumps, dim = make_system(n, BETA)
    end

    # --- 2b: Exact gap computation ---
    local config_l, liouv_result, L_dense
    if n == 4
        config_l = config_l4
        liouv_result = liouv_result4
        L_dense = Matrix(L4)
    else
        config_l = make_liouv_config(n)
        @printf("Building Lindbladian for n=%d...\n", n)
        liouv_result = run_lindbladian(jumps, config_l, ham)
        L_dense = Matrix(liouv_result.liouvillian)
    end

    exact_gap = abs(real(liouv_result.spectral_gap))
    @printf("\nExact gap (ARPACK): %.10f\n\n", exact_gap)

    # --- 2c: Eigenbasis overlap analysis ---
    @printf("--- Eigenbasis Overlap Analysis (n=%d) ---\n", n)
    @printf("Computing dense eigendecomposition of %dx%d Liouvillian...\n",
            size(L_dense, 1), size(L_dense, 2))

    obs, obs_names = build_preset_trajectory_observables(ham, n)
    psi0 = zeros(ComplexF64, dim)
    psi0[end] = 1.0  # Excited state (locked decision)
    rho0 = psi0 * psi0'

    overlap = eigenbasis_overlap_analysis(L_dense, obs, obs_names, rho0)

    @printf("Exact gap (dense eigen): %.10f\n\n", overlap.exact_gap)
    @printf("%-12s  %12s  %16s\n", "Observable", "|c_gap|", "Relative overlap")
    @printf("%-12s  %12s  %16s\n", "-" ^ 12, "-" ^ 12, "-" ^ 16)
    for (i, name) in enumerate(obs_names)
        @printf("%-12s  %12.6f  %16.4f\n",
                name, overlap.gap_mode_overlap[i], overlap.relative_gap_overlap[i])
    end
    println()

    # --- 2d: Trajectory gap estimation ---
    @printf("--- Trajectory Gap Estimation (n=%d, ntraj=%d) ---\n", n, NTRAJ)

    # Mixing time: generous value (5 decay times, minimum 10.0)
    mixing_time = max(5.0 / exact_gap, 10.0)
    @printf("Using mixing_time=%.1f (5/gap=%.1f)\n", mixing_time, 5.0 / exact_gap)

    config_t = make_thermalize_config(n; mixing_time=mixing_time)

    @printf("Running %d trajectories with save_every=%d, skip_initial=0.1...\n",
            NTRAJ, SAVE_EVERY)
    t_start = time()
    gap_result = estimate_spectral_gap(jumps, config_t, psi0, ham;
        ntraj=NTRAJ, save_every=SAVE_EVERY, seed=SEED, skip_initial=0.1)
    t_elapsed = time() - t_start
    @printf("Completed in %.1f seconds.\n\n", t_elapsed)

    # --- 2e: Per-observable fit report ---
    rel_error = abs(gap_result.gap - exact_gap) / exact_gap
    target_rel_error = n == 4 ? TARGET_REL_ERROR_N4 : TARGET_REL_ERROR_N6
    estimation_pass = rel_error < target_rel_error

    @printf("--- Gap Estimation Results (n=%d) ---\n", n)
    @printf("Exact gap:     %.8f\n", exact_gap)
    @printf("Estimated gap: %.8f (best: %s)\n", gap_result.gap, gap_result.best_observable)
    @printf("Relative error: %.4f%% (%.2e)\n", rel_error * 100, rel_error)
    println()

    @printf("Per-observable fits:\n")
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

    @printf("Result: %s (target: relative error < %.1f%%)\n",
            estimation_pass ? "PASS" : "FAIL", target_rel_error * 100)
    println()

    results["n=$n gap estimation"] = estimation_pass

    # --- 2f: Diagnostic explanation if FAIL ---
    if !estimation_pass
        println("=" ^ 50)
        @printf("  DIAGNOSTIC: Why estimation failed for n=%d\n", n)
        println("=" ^ 50)

        # Find best observable's overlap
        best_idx = findfirst(==(gap_result.best_observable), obs_names)
        if best_idx !== nothing
            @printf("Best observable (%s) gap overlap: |c_gap| = %.6f, relative = %.4f\n",
                    gap_result.best_observable,
                    overlap.gap_mode_overlap[best_idx],
                    overlap.relative_gap_overlap[best_idx])
        end

        # Check if overlap is low
        max_overlap_idx = argmax(overlap.gap_mode_overlap)
        max_overlap = overlap.gap_mode_overlap[max_overlap_idx]
        @printf("Strongest gap-mode observable: %s (|c_gap| = %.6f)\n",
                obs_names[max_overlap_idx], max_overlap)

        if max_overlap < 0.01
            @printf("\nDiagnosis: Observable has insufficient overlap with gap mode.\n")
            @printf("All observables have very low gap-mode coupling.\n")
        else
            @printf("\nDiagnosis: Gap mode coupling is sufficient.\n")
            @printf("Possible cause: discrete-step Kraus effect (known systematic bias).\n")
            @printf("\nPer-observable gap/exact ratios:\n")
            for (i, name) in enumerate(gap_result.observable_names)
                fit = gap_result.per_observable[i]
                if fit.converged && fit.gap > 0
                    @printf("  %-12s  %.4f\n", name, fit.gap / exact_gap)
                else
                    @printf("  %-12s  (not converged or gap<=0)\n", name)
                end
            end
        end
        println()
    end
end

# ============================================================================
# Section 3: Final Summary
# ============================================================================

println("=" ^ 70)
println("  FINAL SUMMARY")
println("=" ^ 70)
println()

for (name, pass) in sort(collect(results); by=first)
    @printf("  %-30s  %s\n", name, pass ? "PASS" : "FAIL")
end
println()

all_pass = all(values(results))
@printf("Overall: %s\n", all_pass ? "ALL PASS" : "SOME FAILED")
println("=" ^ 70)
