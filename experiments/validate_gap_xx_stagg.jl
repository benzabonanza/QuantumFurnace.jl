#!/usr/bin/env julia
# ============================================================================
# XX_stagg Disordered Heisenberg Gap Validation Script
# ============================================================================
#
# Tests XX_stagg = Sum((-1)^i * X_i * X_{i+1}) / n (staggered nearest-neighbor
# XX correlation, periodic) as a new observable for disordered Heisenberg gap
# estimation. Uses ONLY H + XX_stagg (no preset 8-observable bundle).
#
# Motivation (Quick Tasks 25-28):
#   - Quick-25: n=6 gap mode is in k=pi momentum sector
#   - Quick-26: Mz_stagg (k=pi, Z-type) has zero overlap due to SU(2) protection
#   - Quick-27: XZ_stagg (k=pi, SU(2)-breaking) still zero -- discrete symmetry
#   - Quick-28: Disorder breaks ALL symmetry protection; Mz_stagg achieves
#     |c_gap|=0.120 for n=6, but gap estimation biased at 34% (selection artifact)
#
# Hypothesis: XX_stagg combines k=pi momentum (staggered sign) with XX two-site
# correlation. Different structure from Mz_stagg may yield different overlap
# profile and potentially better gap estimation.
#
# Parameters (locked, same as Quick-28):
#   - beta = 10.0, delta = 0.01
#   - System sizes: n = 4, n = 6
#   - Trajectories: 20,000
#   - Domain: TimeDomain, with_coherent=false
#   - Initial state: psi0[end] = 1.0 (excited state)
#   - Disordered periodic Heisenberg chain
#
# Usage:
#   cd QuantumFurnace.jl && julia --threads=3 --project experiments/validate_gap_xx_stagg.jl
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

# Tracking results
results = Dict{String, Any}()

println("=" ^ 70)
println("  XX_STAGG DISORDERED HEISENBERG GAP VALIDATION")
println("  beta=$BETA, delta=$DELTA, ntraj=$NTRAJ, seed=$SEED")
println("  Observables: H + XX_stagg ONLY (no preset bundle)")
println("  XX_stagg = Sum((-1)^i * X_i * X_{i+1}) / n  (periodic)")
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

    num_qubits = n
    dim = 2^n

    # --- Section 1a: System setup ---
    @printf("Loading disordered Hamiltonian for n=%d...\n", n)
    ham, jumps, dim = make_system(n, BETA)
    V = ham.eigvecs

    # --- Section 1b: Exact gap via run_lindbladian + ARPACK ---
    config_l = make_liouv_config(n)
    @printf("Building Lindbladian for n=%d...\n", n)
    liouv_result = run_lindbladian(jumps, config_l, ham)
    L_dense = Matrix(liouv_result.liouvillian)

    exact_gap = abs(real(liouv_result.spectral_gap))
    @printf("\nExact gap (ARPACK): %.10f\n\n", exact_gap)

    # --- Section 1c: Build custom observables H_eigen and XX_stagg_eigen MANUALLY ---
    @printf("--- Building Custom Observables (n=%d) ---\n", n)

    # H in eigenbasis (diagonal)
    H_eigen = Matrix{ComplexF64}(diagm(ComplexF64.(ham.eigvals)))
    @printf("H_eigen: %dx%d diagonal matrix\n", size(H_eigen)...)

    # XX_stagg = Sum((-1)^i * X_i * X_{i+1}) / n  (staggered nearest-neighbor XX, periodic)
    XX_stagg_comp = zeros(ComplexF64, dim, dim)
    for i in 1:num_qubits
        sign = (-1)^i
        XX_stagg_comp .+= sign .* Matrix{ComplexF64}(pad_term([X, X], num_qubits, i; periodic=true))
    end
    XX_stagg_comp ./= num_qubits
    XX_stagg_eigen = Matrix{ComplexF64}(V' * XX_stagg_comp * V)
    @printf("XX_stagg_eigen: %dx%d matrix (staggered XX, periodic bonds)\n", size(XX_stagg_eigen)...)

    # Custom observable bundle: ONLY H + XX_stagg
    custom_obs = [H_eigen, XX_stagg_eigen]
    custom_names = ["H", "XX_stagg"]

    # Initial state: excited state (locked decision)
    psi0 = zeros(ComplexF64, dim)
    psi0[end] = 1.0
    rho0 = psi0 * psi0'
    println()

    # --- Section 1d: Eigenbasis overlap analysis ---
    @printf("--- Eigenbasis Overlap Analysis (n=%d, disordered, H + XX_stagg) ---\n", n)
    @printf("Computing dense eigendecomposition of %dx%d Liouvillian...\n",
            size(L_dense, 1), size(L_dense, 2))

    overlap = eigenbasis_overlap_analysis(L_dense, custom_obs, custom_names, rho0)

    @printf("Exact gap (dense eigen): %.10f\n\n", overlap.exact_gap)
    @printf("%-12s  %12s  %16s\n", "Observable", "|c_gap|", "Relative overlap")
    @printf("%-12s  %12s  %16s\n", "-" ^ 12, "-" ^ 12, "-" ^ 16)
    for (i, name) in enumerate(custom_names)
        @printf("%-12s  %12.6f  %16.4f\n",
                name, overlap.gap_mode_overlap[i], overlap.relative_gap_overlap[i])
    end
    println()

    # Symmetry broken check
    for (i, name) in enumerate(custom_names)
        broken = overlap.gap_mode_overlap[i] > 0.001
        @printf("%s gap-mode coupling: |c_gap| = %.6f -> %s\n",
                name, overlap.gap_mode_overlap[i],
                broken ? "COUPLED (|c_gap| > 0.001)" : "NOT COUPLED (|c_gap| <= 0.001)")
    end
    println()

    results["n=$n H |c_gap|"] = overlap.gap_mode_overlap[1]
    results["n=$n XX_stagg |c_gap|"] = overlap.gap_mode_overlap[2]

    # --- Section 1e: Full expansion coefficient spectrum for XX_stagg ---
    @printf("--- Full Expansion Coefficient Spectrum for XX_stagg (n=%d) ---\n", n)

    n_modes = size(overlap.overlap_coefficients, 2)
    n_show = min(20, n_modes)

    @printf("Showing |c_k| for modes k=1..%d (of %d total) for XX_stagg (row 2):\n\n", n_show, n_modes)
    @printf("%-6s  %14s  %14s  %14s\n", "Mode", "|c_k|", "Re(lambda_k)", "Im(lambda_k)")
    @printf("%-6s  %14s  %14s  %14s\n", "-" ^ 6, "-" ^ 14, "-" ^ 14, "-" ^ 14)
    for k in 1:n_show
        ck = abs(overlap.overlap_coefficients[2, k])
        lk = overlap.eigenvalues[k]
        @printf("k=%-4d  %14.8f  %14.8f  %14.8f\n", k, ck, real(lk), imag(lk))
    end
    println()

    # Rank top 5 modes by |c_k| for XX_stagg (excluding steady state k=1)
    all_ck = [abs(overlap.overlap_coefficients[2, k]) for k in 1:n_modes]
    # Sort by |c_k| descending, exclude k=1 (steady state)
    mode_indices = sortperm(all_ck[2:end]; rev=true) .+ 1  # offset by 1 to skip k=1
    top5 = mode_indices[1:min(5, length(mode_indices))]

    @printf("Top 5 modes by |c_k| for XX_stagg (excluding steady state):\n")
    @printf("%-6s  %14s  %14s  %14s\n", "Mode", "|c_k|", "Re(lambda_k)", "Is gap?")
    @printf("%-6s  %14s  %14s  %14s\n", "-" ^ 6, "-" ^ 14, "-" ^ 14, "-" ^ 14)
    for k in top5
        ck = abs(overlap.overlap_coefficients[2, k])
        lk = overlap.eigenvalues[k]
        is_gap = k == 2 ? "YES (gap mode)" : "no"
        @printf("k=%-4d  %14.8f  %14.8f  %s\n", k, ck, real(lk), is_gap)
    end
    println()

    results["n=$n XX_stagg dominant_mode"] = top5[1]
    results["n=$n XX_stagg dominant_ck"] = all_ck[top5[1]]

    # --- Section 1f: Eigenvalue spectrum near gap (modes 1-5) ---
    @printf("--- Eigenvalue Spectrum Near Gap (n=%d, modes 1-5) ---\n", n)
    n_eig_show = min(5, n_modes)
    @printf("%-6s  %16s  %16s\n", "Mode", "Re(lambda)", "Im(lambda)")
    @printf("%-6s  %16s  %16s\n", "-" ^ 6, "-" ^ 16, "-" ^ 16)
    for k in 1:n_eig_show
        lk = overlap.eigenvalues[k]
        @printf("k=%-4d  %16.10f  %16.10f\n", k, real(lk), imag(lk))
    end
    println()

    # Near-degenerate mode detection
    @printf("Near-degenerate mode check (|Re(lambda_i) - Re(lambda_j)| < 0.1 * exact_gap):\n")
    degen_threshold = 0.1 * exact_gap
    found_degen = false
    for i in 2:n_eig_show
        for j in (i+1):n_eig_show
            diff = abs(real(overlap.eigenvalues[i]) - real(overlap.eigenvalues[j]))
            if diff < degen_threshold
                @printf("  NEAR-DEGENERATE: modes %d and %d, |delta Re(lambda)| = %.2e (threshold %.2e)\n",
                        i, j, diff, degen_threshold)
                found_degen = true
            end
        end
    end
    if !found_degen
        @printf("  No near-degenerate modes found among modes 2-5 (threshold: %.2e)\n", degen_threshold)
    end
    println()

    results["n=$n near_degenerate"] = found_degen

    # --- Section 1g: Trajectory gap estimation ---
    @printf("--- Trajectory Gap Estimation (n=%d, ntraj=%d, H + XX_stagg only) ---\n", n, NTRAJ)

    mixing_time = max(5.0 / exact_gap, 10.0)
    @printf("Using mixing_time=%.1f (5/gap=%.1f)\n", mixing_time, 5.0 / exact_gap)

    config_t = make_thermalize_config(n; mixing_time=mixing_time)

    @printf("Running %d trajectories with save_every=%d, skip_initial=0.1...\n",
            NTRAJ, SAVE_EVERY)
    t_start = time()
    gap_result = estimate_spectral_gap(jumps, config_t, psi0, ham;
        observables=custom_obs, observable_names=custom_names,
        ntraj=NTRAJ, save_every=SAVE_EVERY, seed=SEED, skip_initial=0.1)
    t_elapsed = time() - t_start
    @printf("Completed in %.1f seconds.\n\n", t_elapsed)

    # --- Section 1h: Deep diagnostics (ALWAYS print) ---
    rel_error = abs(gap_result.gap - exact_gap) / exact_gap

    @printf("--- Gap Estimation Results (n=%d, disordered, H + XX_stagg) ---\n", n)
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

    results["n=$n exact_gap"] = exact_gap
    results["n=$n estimated_gap"] = gap_result.gap
    results["n=$n rel_error"] = rel_error
    results["n=$n best_obs"] = gap_result.best_observable

    # Per-observable gap/exact ratios (ALWAYS print)
    @printf("Per-observable gap/exact ratios:\n")
    for (i, name) in enumerate(gap_result.observable_names)
        fit = gap_result.per_observable[i]
        if fit.converged && fit.gap > 0
            ratio = fit.gap / exact_gap
            @printf("  %-12s  gap/exact = %.4f  (gap=%.6f)\n", name, ratio, fit.gap)
        else
            @printf("  %-12s  (not converged or gap<=0)\n", name)
        end
    end
    println()

    # Raw time series via run_observable_trajectories (SEPARATE call)
    @printf("--- Raw Time Series (n=%d, XX_stagg) ---\n", n)
    @printf("Running run_observable_trajectories for raw time series...\n")
    traj_result = run_observable_trajectories(jumps, config_t, psi0, ham;
        observables=custom_obs, save_every=SAVE_EVERY,
        ntraj=NTRAJ, total_time=config_t.mixing_time, delta=DELTA,
        seed=SEED)

    n_times = length(traj_result.times)
    sample_indices = [1, div(n_times, 4), div(n_times, 2), div(3 * n_times, 4), n_times]
    # Ensure valid unique indices
    sample_indices = unique(max.(1, min.(sample_indices, n_times)))

    @printf("Total time points: %d\n", n_times)
    @printf("Sample time points for XX_stagg (row 2 of measurements_mean):\n\n")
    @printf("%-8s  %14s  %14s\n", "Index", "Time", "XX_stagg value")
    @printf("%-8s  %14s  %14s\n", "-" ^ 8, "-" ^ 14, "-" ^ 14)
    for idx in sample_indices
        t = traj_result.times[idx]
        val = traj_result.measurements_mean[2, idx]
        @printf("%-8d  %14.6f  %14.8f\n", idx, t, val)
    end
    println()

    # Also show H (row 1) for comparison
    @printf("Sample time points for H (row 1 of measurements_mean):\n\n")
    @printf("%-8s  %14s  %14s\n", "Index", "Time", "H value")
    @printf("%-8s  %14s  %14s\n", "-" ^ 8, "-" ^ 14, "-" ^ 14)
    for idx in sample_indices
        t = traj_result.times[idx]
        val = traj_result.measurements_mean[1, idx]
        @printf("%-8d  %14.6f  %14.8f\n", idx, t, val)
    end
    println()

    # Investigation section if >10% error
    if rel_error > 0.10
        println("=" ^ 50)
        @printf("  INVESTIGATION: >10%% error for n=%d (%.2f%%)\n", n, rel_error * 100)
        println("=" ^ 50)
        println()

        # Is the dominant mode in XX_stagg's expansion actually the gap mode?
        dominant_k = top5[1]
        dominant_ck = all_ck[dominant_k]
        gap_ck = all_ck[2]

        @printf("XX_stagg dominant mode: k=%d (|c_k| = %.8f)\n", dominant_k, dominant_ck)
        @printf("XX_stagg gap mode (k=2): |c_2| = %.8f\n", gap_ck)

        if dominant_k == 2
            @printf("-> Dominant mode IS the gap mode. Bias is from fitting/Kraus effect.\n")
        else
            @printf("-> Dominant mode is NOT the gap mode (k=%d != k=2).\n", dominant_k)
            @printf("   The fitted decay rate may reflect mode k=%d (Re(lambda) = %.8f)\n",
                    dominant_k, real(overlap.eigenvalues[dominant_k]))
            @printf("   rather than the gap (Re(lambda_2) = %.8f)\n",
                    real(overlap.eigenvalues[2]))
            @printf("   Ratio: dominant_mode_rate / gap_rate = %.4f\n",
                    abs(real(overlap.eigenvalues[dominant_k])) / exact_gap)
        end
        println()

        # Compare H and XX_stagg gap/exact ratios
        @printf("Comparison of H vs XX_stagg gap estimation:\n")
        for (i, name) in enumerate(gap_result.observable_names)
            fit = gap_result.per_observable[i]
            if fit.converged && fit.gap > 0
                @printf("  %-12s: gap=%.6f, gap/exact=%.4f, R2=%.4f, |c_gap|=%.6f\n",
                        name, fit.gap, fit.gap / exact_gap, fit.r_squared,
                        overlap.gap_mode_overlap[i])
            end
        end
        println()

        # Check if time series shows clear exponential decay
        @printf("Time series behavior analysis (XX_stagg):\n")
        first_val = traj_result.measurements_mean[2, 1]
        last_val = traj_result.measurements_mean[2, end]
        mid_val = traj_result.measurements_mean[2, div(n_times, 2)]
        @printf("  Initial value: %.8f\n", first_val)
        @printf("  Mid-point:     %.8f\n", mid_val)
        @printf("  Final value:   %.8f\n", last_val)
        @printf("  Dynamic range: %.8f (|first - last|)\n", abs(first_val - last_val))
        if abs(first_val - last_val) < 1e-6
            @printf("  WARNING: Very small dynamic range -- XX_stagg may have minimal signal\n")
        end
        println()
    end

    println()
end

# ============================================================================
# Section 2: Final Summary
# ============================================================================

println("=" ^ 70)
println("  FINAL SUMMARY")
println("=" ^ 70)
println()

@printf("%-25s  %12s  %12s\n", "Metric", "n=4", "n=6")
@printf("%-25s  %12s  %12s\n", "-" ^ 25, "-" ^ 12, "-" ^ 12)
@printf("%-25s  %12.6f  %12.6f\n", "Exact gap",
        results["n=4 exact_gap"], results["n=6 exact_gap"])
@printf("%-25s  %12.6f  %12.6f\n", "Estimated gap",
        results["n=4 estimated_gap"], results["n=6 estimated_gap"])
@printf("%-25s  %12.4f%%  %11.4f%%\n", "Relative error",
        results["n=4 rel_error"] * 100, results["n=6 rel_error"] * 100)
@printf("%-25s  %12s  %12s\n", "Best observable",
        results["n=4 best_obs"], results["n=6 best_obs"])
@printf("%-25s  %12.6f  %12.6f\n", "H |c_gap|",
        results["n=4 H |c_gap|"], results["n=6 H |c_gap|"])
@printf("%-25s  %12.6f  %12.6f\n", "XX_stagg |c_gap|",
        results["n=4 XX_stagg |c_gap|"], results["n=6 XX_stagg |c_gap|"])
@printf("%-25s  %12d  %12d\n", "XX_stagg dominant mode",
        results["n=4 XX_stagg dominant_mode"], results["n=6 XX_stagg dominant_mode"])
@printf("%-25s  %12.6f  %12.6f\n", "XX_stagg dominant |c_k|",
        results["n=4 XX_stagg dominant_ck"], results["n=6 XX_stagg dominant_ck"])
@printf("%-25s  %12s  %12s\n", "Near-degenerate modes",
        results["n=4 near_degenerate"] ? "YES" : "no",
        results["n=6 near_degenerate"] ? "YES" : "no")
println()

# Key findings
println("=" ^ 50)
println("  KEY FINDINGS")
println("=" ^ 50)
println()

for n in [4, 6]
    cgap = results["n=$n XX_stagg |c_gap|"]
    coupled = cgap > 0.001
    @printf("n=%d: XX_stagg |c_gap| = %.6f -> %s\n",
            n, cgap, coupled ? "COUPLED to gap mode" : "NOT coupled to gap mode")
    if coupled
        @printf("     Gap estimation: %.4f%% relative error (best: %s)\n",
                results["n=$n rel_error"] * 100, results["n=$n best_obs"])
    end
end
println()

println("=" ^ 70)
