#!/usr/bin/env julia
# ============================================================================
# Disordered Heisenberg Gap Validation Script
# ============================================================================
#
# Tests whether disordered Heisenberg Hamiltonians (random Z-field per site)
# break the symmetry protection on the Lindbladian gap mode for n=6.
#
# Motivation (Quick Tasks 25-27):
#   - Quick-25: n=6 gap mode is in k=pi momentum sector; all 5 original
#     observables are k=0 and have zero gap-mode overlap.
#   - Quick-26: Added Mz_stagg (k=pi) and Z1 (all k); still zero overlap
#     due to SU(2) spin-rotation symmetry of Heisenberg chain.
#   - Quick-27: Added XZ_stagg (breaks SU(2), has k=pi); still zero overlap
#     -- protection extends beyond SU(2) (discrete symmetries: parity, spin-flip).
#
# Hypothesis: Random Z-fields break ALL symmetries (translational, SU(2),
# parity, spin-flip). If the gap mode couples to observables after disorder,
# trajectory-based gap estimation becomes viable for n=6.
#
# This script uses pre-generated disordered Hamiltonians loaded via
# `load_hamiltonian("heis", n; beta=10.0)` from BSON files.
#
# Parameters (locked decisions):
#   - beta = 10.0, delta = 0.01
#   - System sizes: n = 4, n = 6
#   - Trajectories: 20,000
#   - Target: n=4 < 1%; n=6 < 20% (generous -- unknown territory)
#   - Domain: TimeDomain with with_coherent=false
#   - Initial state: psi0[end] = 1.0 (excited state)
#   - Disordered periodic Heisenberg chain
#
# Usage:
#   cd QuantumFurnace.jl && julia --project experiments/validate_gap_disordered.jl
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
const TARGET_REL_ERROR_N4 = 1e-2   # 1% for n=4
const TARGET_REL_ERROR_N6 = 0.20   # 20% for n=6 (generous -- disordered system, unknown territory)
const SEED = 42

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
println("  DISORDERED HEISENBERG GAP VALIDATION")
println("  beta=$BETA, delta=$DELTA, ntraj=$NTRAJ, seed=$SEED")
println("  Testing: Does disorder break n=6 symmetry protection?")
println("=" ^ 70)
println()

# ============================================================================
# Section 1: Per-system-size analysis (n=4, n=6)
# ============================================================================

for n in [4, 6]
    println("=" ^ 70)
    @printf("  System size n=%d (dim=%d, Liouvillian %dx%d)\n",
            n, 2^n, (2^n)^2, (2^n)^2)
    println("  Hamiltonian: disordered periodic Heisenberg (random Z-field)")
    println("=" ^ 70)
    println()

    # --- 1a: System setup ---
    @printf("Loading disordered Hamiltonian for n=%d...\n", n)
    ham, jumps, dim = make_system(n, BETA)

    # --- 1b: Exact gap computation ---
    config_l = make_liouv_config(n)
    @printf("Building Lindbladian for n=%d...\n", n)
    liouv_result = run_lindbladian(jumps, config_l, ham)
    L_dense = Matrix(liouv_result.liouvillian)

    exact_gap = abs(real(liouv_result.spectral_gap))
    @printf("\nExact gap (ARPACK): %.10f\n\n", exact_gap)

    # --- 1c: Eigenbasis overlap analysis ---
    @printf("--- Eigenbasis Overlap Analysis (n=%d, disordered) ---\n", n)
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

    # KEY RESULT: Did disorder break symmetry protection?
    max_cgap_idx = argmax(overlap.gap_mode_overlap)
    max_cgap = overlap.gap_mode_overlap[max_cgap_idx]
    symmetry_broken = max_cgap > 0.001

    println("=" ^ 50)
    if symmetry_broken
        @printf("KEY RESULT: Disorder DID break symmetry protection for n=%d\n", n)
    else
        @printf("KEY RESULT: Disorder DID NOT break symmetry protection for n=%d\n", n)
    end
    @printf("Max |c_gap| = %.6f (observable: %s)\n", max_cgap, obs_names[max_cgap_idx])
    @printf("Threshold: |c_gap| > 0.001\n")
    println("=" ^ 50)
    println()

    results["n=$n symmetry broken"] = symmetry_broken

    # --- 1d: Trajectory gap estimation ---
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

    # --- 1e: Per-observable fit report ---
    rel_error = abs(gap_result.gap - exact_gap) / exact_gap
    target_rel_error = n == 4 ? TARGET_REL_ERROR_N4 : TARGET_REL_ERROR_N6
    estimation_pass = rel_error < target_rel_error

    @printf("--- Gap Estimation Results (n=%d, disordered) ---\n", n)
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

    # --- 1f: Diagnostic explanation if FAIL ---
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
        @printf("Strongest gap-mode observable: %s (|c_gap| = %.6f)\n",
                obs_names[max_cgap_idx], max_cgap)

        if max_cgap < 0.01
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
# Section 2: Final Summary
# ============================================================================

println("=" ^ 70)
println("  FINAL SUMMARY")
println("=" ^ 70)
println()

for (name, pass) in sort(collect(results); by=first)
    @printf("  %-30s  %s\n", name, pass ? "PASS" : "FAIL")
end
println()

# Key question answer
n6_broken = get(results, "n=6 symmetry broken", false)
println("=" ^ 50)
println("  KEY QUESTION: Did disorder break n=6 symmetry protection?")
if n6_broken
    println("  ANSWER: YES -- Disordered Hamiltonians enable gap-mode coupling")
    println("  Trajectory-based gap estimation is viable for n=6 with disorder.")
else
    println("  ANSWER: NO -- Even with disorder, gap mode remains protected")
    println("  Further investigation needed for n=6 gap estimation.")
end
println("=" ^ 50)
println()

all_pass = all(values(results))
@printf("Overall: %s\n", all_pass ? "ALL PASS" : "SOME FAILED")
println("=" ^ 70)
