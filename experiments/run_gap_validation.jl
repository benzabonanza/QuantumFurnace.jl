# ============================================================================
# Gap Estimation Validation: n=6 with Two-Site Correlations
# ============================================================================
#
# Validates spectral gap estimation improvement from adding XX_avg, YY_avg,
# ZZ_avg observables. Runs n=6 Heisenberg chain with high trajectory count
# and prints per-observable fit summary to evaluate which observables best
# capture the Lindbladian spectral gap.
#
# Usage:
#   julia --project=. experiments/run_gap_validation.jl
#
# Output: Gap comparison (fitted vs exact) and per-observable fit details.
# ============================================================================

using QuantumFurnace
using LinearAlgebra
using Printf

# Constants (same as run_sweep.jl)
const DELTA = 0.01
const NUM_ENERGY_BITS = 12
const W0 = 0.05
const T0 = 2pi / (2^NUM_ENERGY_BITS * W0)
const NUM_TROTTER_STEPS_PER_T0 = 10
const SEED = 42

function main()
    n = 6
    beta = 5.0
    ntraj = 10_000
    mixing_time = 20.0

    println("=== Gap Estimation Validation: n=$n, beta=$beta, ntraj=$ntraj ===")

    # Build system
    hamiltonian = HamHam([[X, X], [Y, Y], [Z, Z]], [1.0, 1.0, 1.0], n, beta; periodic=true)
    trotter = TrottTrott(hamiltonian, T0, NUM_TROTTER_STEPS_PER_T0)

    # Build jump operators (same as run_sweep.jl)
    jump_paulis = [[X], [Y], [Z]]
    n_jumps = length(jump_paulis) * n
    norm_factor = sqrt(n_jumps)
    jumps = JumpOp[]
    for pauli in jump_paulis
        for site in 1:n
            op = Matrix(pad_term(pauli, n, site)) ./ norm_factor
            in_eigen = trotter.eigvecs' * op * trotter.eigvecs
            push!(jumps, JumpOp(op, in_eigen, op == transpose(op), op == op'))
        end
    end

    # Trajectory config (ThermalizeConfig for estimate_spectral_gap)
    config = ThermalizeConfig(
        num_qubits=n, with_coherent=true, with_linear_combination=true,
        domain=TrotterDomain(), beta=beta, sigma=1.0 / beta,
        a=beta / 30.0, b=0.4,
        num_energy_bits=NUM_ENERGY_BITS, w0=W0, t0=T0,
        num_trotter_steps_per_t0=NUM_TROTTER_STEPS_PER_T0,
        mixing_time=mixing_time, delta=DELTA,
    )

    # Liouvillian config (LiouvConfig for exact gap computation)
    liouv_config = LiouvConfig(
        num_qubits=n, with_coherent=true, with_linear_combination=true,
        domain=TrotterDomain(), beta=beta, sigma=1.0 / beta,
        a=beta / 30.0, b=0.4,
        num_energy_bits=NUM_ENERGY_BITS, w0=W0, t0=T0,
        num_trotter_steps_per_t0=NUM_TROTTER_STEPS_PER_T0,
    )

    # Initial state -- use excited state (psi0[end]=1) for validation
    # Ground state at high beta is near Gibbs, giving no decay signal
    dim = 2^n
    psi0 = zeros(ComplexF64, dim)
    psi0[end] = 1.0

    # Estimate spectral gap (uses all 5 observables now)
    println("Running estimate_spectral_gap with ntraj=$ntraj...")
    t0_wall = time()
    result = estimate_spectral_gap(
        jumps, config, psi0, hamiltonian;
        ntraj=ntraj, save_every=10, seed=SEED,
        trotter=trotter, skip_initial=0.1,
    )
    wall = time() - t0_wall
    @printf("Wall time: %.1fs\n\n", wall)

    # Exact gap via Liouvillian eigendecomposition
    println("Computing exact Liouvillian gap...")
    exact_result = run_lindbladian(jumps, liouv_config, hamiltonian; trotter=trotter)

    # Cross-validate
    cv = cross_validate_gap(result, exact_result)

    # Results
    println("\n=== Results ===")
    @printf("Exact gap:   %.6f\n", cv.exact_gap)
    @printf("Fitted gap:  %.6f (best observable: %s)\n", cv.fitted_gap, result.best_observable)
    @printf("Factor:      %.4fx\n", cv.fitted_gap / cv.exact_gap)
    @printf("Rel error:   %.4f\n", cv.relative_error)
    @printf("Within CI:   %s\n", cv.within_ci)
    @printf("CI:          [%.6f, %.6f]\n", result.gap_ci...)

    println("\n=== Per-Observable Fits ===")
    @printf("%-10s %10s %10s %10s\n", "Name", "Gap", "R-squared", "Converged")
    for (name, fit) in zip(result.observable_names, result.per_observable)
        @printf("%-10s %10.6f %10.4f %10s\n", name, fit.gap, fit.r_squared, fit.converged)
    end
end

main()
