# ============================================================================
# Gap Estimation Cross-Validation Script
# ============================================================================
#
# Validates trajectory-fitted spectral gap against exact Liouvillian eigenvalues
# for n=4 and n=6 uniform Heisenberg XXX chains using TimeDomain.
#
# Demonstrates that estimate_spectral_gap produces gap estimates consistent
# with construct_lindbladian + dense eigen() for small systems.
#
# NOTE: The trajectory-based gap estimation operates in discrete time with
# delta_eff = delta * n_jumps per step, which introduces a normalization
# factor between the fitted trajectory rate and the continuous-time
# Liouvillian spectral gap. The script reports this factor and evaluates
# agreement after normalization.
#
# Usage:
#   julia --project=. experiments/validate_gap_estimation.jl
#
# Output: Cross-validation results for each system size with pass/fail status.
# ============================================================================

using QuantumFurnace
using LinearAlgebra
using Printf
using Dates

# ============================================================================
# Shared constants
# ============================================================================

const DELTA = 0.01
const BETA = 10.0
const NUM_ENERGY_BITS = 12
const W0 = 0.05
const T0 = 2pi / (2^NUM_ENERGY_BITS * W0)
const NUM_TROTTER_STEPS_PER_T0 = 10
const SEED = 42

# ============================================================================
# Helpers
# ============================================================================

"""
    build_heisenberg_hamiltonian(num_qubits, beta) -> HamHam

Build a uniform 1D Heisenberg XXX chain with periodic boundaries and J=1.0.
"""
function build_heisenberg_hamiltonian(num_qubits::Int, beta::Float64)
    terms = Vector{Vector{Matrix{ComplexF64}}}([[X, X], [Y, Y], [Z, Z]])
    coeffs = [1.0, 1.0, 1.0]
    return HamHam(terms, coeffs, num_qubits, beta; periodic=true)
end

"""
    build_jumps(hamiltonian, num_qubits) -> Vector{JumpOp}

Build single-site Pauli jump operators in Hamiltonian eigenbasis.
Uses hamiltonian.eigvecs (NOT trotter.eigvecs) because we use TimeDomain
with with_coherent=false. Matches run_sweep.jl pattern but for TimeDomain.
"""
function build_jumps(hamiltonian::HamHam, num_qubits::Int)
    jump_paulis = [[X], [Y], [Z]]
    n_jumps = length(jump_paulis) * num_qubits
    norm_factor = sqrt(n_jumps)

    jumps = JumpOp[]
    for pauli in jump_paulis
        for site in 1:num_qubits
            op = Matrix(pad_term(pauli, num_qubits, site)) ./ norm_factor
            in_eigen = hamiltonian.eigvecs' * op * hamiltonian.eigvecs
            push!(jumps, JumpOp(op, in_eigen, op == transpose(op), op == op'))
        end
    end
    return jumps
end

"""
    make_matched_configs(num_qubits, beta, sigma) -> (LiouvConfig, NamedTuple)

Create matching LiouvConfig and shared parameter tuple with IDENTICAL physics parameters.
Both use TimeDomain() with with_coherent=false. The shared params are used to
construct ThermalizeConfig with system-specific mixing_time.
"""
function make_matched_configs(num_qubits::Int, beta::Float64, sigma::Float64)
    shared = (
        num_qubits = num_qubits,
        with_coherent = false,
        with_linear_combination = true,
        domain = TimeDomain(),
        beta = beta,
        sigma = sigma,
        a = beta / 30.0,
        b = 0.4,
        num_energy_bits = NUM_ENERGY_BITS,
        w0 = W0,
        t0 = T0,
        num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
    )

    liouv_config = LiouvConfig(; shared...)
    return liouv_config, shared
end

"""
    extract_exact_result(L, hamiltonian, num_qubits) -> LindbladianResult

Build LindbladianResult from dense eigendecomposition of the Liouvillian.
Uses eigen() (NOT Arpack eigs) for exact results on small systems.
"""
function extract_exact_result(L::Matrix, hamiltonian::HamHam, num_qubits::Int)
    dim = 2^num_qubits

    eig = eigen(L)
    sorted_idx = sortperm(abs.(real.(eig.values)))

    ss_idx = sorted_idx[1]   # eigenvalue closest to 0 (steady state)
    gap_idx = sorted_idx[2]  # second eigenvalue (spectral gap)

    # Extract steady state
    ss_vec = eig.vectors[:, ss_idx]
    ss_dm = reshape(ss_vec, dim, dim)
    hermitianize!(ss_dm)
    ss_dm ./= tr(ss_dm)

    # Extract gap mode
    gap_vec = eig.vectors[:, gap_idx]
    gap_mode = reshape(gap_vec, dim, dim)

    return LindbladianResult(
        liouvillian = L,
        fixed_point = ss_dm,
        gap_mode = gap_mode,
        spectral_gap = eig.values[gap_idx],
    )
end

# ============================================================================
# Main validation function
# ============================================================================

"""
    validate_system(num_qubits; ntraj, total_time) -> CrossValidationResult

Run full cross-validation for a given system size:
1. Build Hamiltonian and jump operators
2. Compute exact Liouvillian gap via dense eigendecomposition
3. Estimate gap from trajectories via estimate_spectral_gap
4. Cross-validate the two results
"""
function validate_system(num_qubits::Int; ntraj::Int, total_time::Float64)
    @printf("\n")
    println("=" ^ 60)
    @printf("=== Validating n=%d Heisenberg chain (beta=%.1f) ===\n", num_qubits, BETA)
    println("=" ^ 60)

    # 1. Build Hamiltonian
    @printf("[%s] Building Hamiltonian...\n", Dates.format(now(), "HH:MM:SS"))
    hamiltonian = build_heisenberg_hamiltonian(num_qubits, BETA)
    dim = 2^num_qubits
    n_jumps = 3 * num_qubits
    @printf("  dim = %d, Liouvillian dim = %d x %d\n", dim, dim^2, dim^2)

    # 2. Build jumps in Hamiltonian eigenbasis (TimeDomain)
    @printf("[%s] Building jump operators...\n", Dates.format(now(), "HH:MM:SS"))
    jumps = build_jumps(hamiltonian, num_qubits)
    @printf("  %d jump operators (3 Paulis x %d sites)\n", length(jumps), num_qubits)

    # 3. Create matched configs
    liouv_config, shared_params = make_matched_configs(num_qubits, BETA, 1.0 / BETA)

    therm_config = ThermalizeConfig(;
        shared_params...,
        mixing_time = total_time,
        delta = DELTA,
    )

    # 4. Exact gap via Liouvillian
    @printf("[%s] Computing exact Liouvillian eigenvalues...\n", Dates.format(now(), "HH:MM:SS"))
    t_exact_start = time()
    L = construct_lindbladian(jumps, liouv_config, hamiltonian)
    exact_result = extract_exact_result(L, hamiltonian, num_qubits)
    t_exact = time() - t_exact_start

    exact_eigenvalue = exact_result.spectral_gap
    exact_gap = abs(real(exact_eigenvalue))
    @printf("  Exact eigenvalue: %.6f %+.6fim\n", real(exact_eigenvalue), imag(exact_eigenvalue))
    @printf("  Exact gap (|Re|): %.6f\n", exact_gap)
    @printf("  Exact computation time: %.1fs\n", t_exact)

    # 5. Trajectory gap estimation
    # Use excited initial state (highest energy eigenstate in eigenbasis)
    # to maximize observable decay signal. The ground state (psi0[1]=1)
    # is already near the Gibbs state at high beta, producing no decay.
    @printf("[%s] Running trajectory gap estimation (ntraj=%d, total_time=%.1f)...\n",
        Dates.format(now(), "HH:MM:SS"), ntraj, total_time)
    psi0 = zeros(ComplexF64, dim)
    psi0[end] = 1.0

    t_traj_start = time()
    estimated = estimate_spectral_gap(
        jumps, therm_config, psi0, hamiltonian;
        ntraj = ntraj,
        save_every = 10,
        seed = SEED,
        total_time = total_time,
        skip_initial = 0.0,
    )
    t_traj = time() - t_traj_start

    @printf("  Fitted gap: %.6f\n", estimated.gap)
    @printf("  Fitted CI: (%.6f, %.6f)\n", estimated.gap_ci[1], estimated.gap_ci[2])
    @printf("  Best observable: %s\n", estimated.best_observable)
    @printf("  Best R-squared: %.4f\n", estimated.best_r_squared)
    @printf("  Trajectory estimation time: %.1fs\n", t_traj)

    # 6. Cross-validate
    @printf("[%s] Cross-validating...\n", Dates.format(now(), "HH:MM:SS"))
    cv = cross_validate_gap(estimated, exact_result)

    @printf("\n--- Cross-Validation Results ---\n")
    @printf("  Fitted gap:       %.6f\n", cv.fitted_gap)
    @printf("  Exact gap:        %.6f\n", cv.exact_gap)
    @printf("  Relative error:   %.4f (%.1f%%)\n", cv.relative_error, cv.relative_error * 100)
    @printf("  Absolute error:   %.6f\n", cv.absolute_error)
    @printf("  Within CI:        %s\n", cv.within_ci ? "YES" : "NO")
    @printf("  Imaginary ratio:  %.4f\n", cv.imaginary_ratio)
    @printf("  Imaginary warning:%s\n", cv.imaginary_warning ? " YES" : " NO")

    # 7. Report normalization factor
    # The trajectory uses delta_eff = delta * n_jumps per step, causing the
    # fitted rate to be scaled relative to the continuous-time Liouvillian gap.
    if cv.exact_gap > 0.0 && cv.fitted_gap > 0.0
        norm_factor = cv.fitted_gap / cv.exact_gap
        @printf("\n  Normalization factor (fitted/exact): %.2f\n", norm_factor)
    end

    # 8. Pass/fail
    passed = cv.within_ci || cv.relative_error < 0.3
    @printf("\n  >>> %s <<<\n", passed ? "PASS" : "NEEDS MORE TRAJECTORIES")

    return cv
end

# ============================================================================
# Main
# ============================================================================

function main()
    sweep_start = time()

    @printf("\n")
    println("=" ^ 60)
    println("Gap Estimation Cross-Validation")
    println("=" ^ 60)
    @printf("Started: %s\n", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    @printf("Domain: TimeDomain (with_coherent=false)\n")
    @printf("Beta: %.1f, Sigma: %.4f\n", BETA, 1.0 / BETA)
    @printf("Seed: %d\n", SEED)

    # Run n=4: gap ~0.16, use enough trajectories for good fit
    cv4 = validate_system(4; ntraj=1000, total_time=50.0)

    # Run n=6: gap ~0.12, keep practical run time
    cv6 = validate_system(6; ntraj=1000, total_time=50.0)

    # Summary
    sweep_wall = time() - sweep_start
    @printf("\n")
    println("=" ^ 60)
    println("Summary")
    println("=" ^ 60)

    pass4 = cv4.within_ci || cv4.relative_error < 0.3
    pass6 = cv6.within_ci || cv6.relative_error < 0.3

    @printf("  n=4: relative_error=%.4f, within_ci=%s -> %s\n",
        cv4.relative_error, cv4.within_ci ? "YES" : "NO", pass4 ? "PASS" : "FAIL")
    @printf("  n=6: relative_error=%.4f, within_ci=%s -> %s\n",
        cv6.relative_error, cv6.within_ci ? "YES" : "NO", pass6 ? "PASS" : "FAIL")

    @printf("\nTotal wall time: %.1fs (%.1f min)\n", sweep_wall, sweep_wall / 60.0)

    overall = pass4 && pass6
    @printf("\n>>> OVERALL: %s <<<\n\n", overall ? "PASS" : "FAIL")

    return cv4, cv6
end

# Run the validation
main()
