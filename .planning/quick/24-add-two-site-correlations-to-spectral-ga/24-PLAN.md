---
phase: 24-add-two-site-correlations
plan: 01
type: execute
wave: 1
depends_on: []
files_modified:
  - src/convergence.jl
  - test/test_gap_estimation.jl
  - experiments/run_gap_validation.jl
autonomous: true
must_haves:
  truths:
    - "build_gap_estimation_observables returns 5 observables: H, Mz, XX_avg, YY_avg, ZZ_avg"
    - "XX_avg, YY_avg, ZZ_avg are each an average over all nearest-neighbor bonds (periodic)"
    - "All observables are correctly transformed to the eigenbasis (Hamiltonian or Trotter)"
    - "estimate_spectral_gap uses all 5 observables for fitting and selection"
    - "n=6 gap estimation factor is measured with higher trajectory count"
  artifacts:
    - path: "src/convergence.jl"
      provides: "Updated build_gap_estimation_observables with XX_avg, YY_avg, ZZ_avg"
      contains: "XX_avg"
    - path: "test/test_gap_estimation.jl"
      provides: "Updated tests reflecting 5 observables"
      contains: "XX_avg"
    - path: "experiments/run_gap_validation.jl"
      provides: "Validation script for n=6 gap estimation with higher ntraj"
      contains: "estimate_spectral_gap"
  key_links:
    - from: "src/convergence.jl"
      to: "src/gap_estimation.jl"
      via: "build_gap_estimation_observables called by estimate_spectral_gap"
      pattern: "build_gap_estimation_observables"
---

<objective>
Add averaged two-site correlation observables (XX_avg, YY_avg, ZZ_avg) to the spectral gap estimation observable bundle, update tests, and create a validation script to measure n=6 gap estimation improvement.

Purpose: The current spectral gap estimation uses only H and Mz. Adding two-site correlations provides observables that may overlap better with the first excited Lindbladian mode, potentially improving the gap estimate (reducing the ~1.46x residual factor at n=6).

Output: Updated `build_gap_estimation_observables` returning 5 observables, updated tests, and a standalone validation script.
</objective>

<execution_context>
@./.claude/get-shit-done/workflows/execute-plan.md
@./.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@src/convergence.jl
@src/gap_estimation.jl
@src/constants.jl
@test/test_gap_estimation.jl
@experiments/run_sweep.jl
</context>

<tasks>

<task type="auto">
  <name>Task 1: Add XX_avg, YY_avg, ZZ_avg to build_gap_estimation_observables</name>
  <files>src/convergence.jl</files>
  <action>
In `build_gap_estimation_observables` (line 127-145 of convergence.jl), add three new averaged two-site correlation observables after the existing H and Mz construction:

For each Pauli pair PP in {[X,X], [Y,Y], [Z,Z]}:
1. Build the sum of nearest-neighbor PP correlations across all bonds: `PP_sum = sum_{i=1}^{n} PP_{i,(i%n)+1}` using `pad_term([P, P], num_qubits, i; periodic=true)` -- same pattern as `build_convergence_observables` lines 31-38.
2. Divide by `num_qubits` to get per-bond average (analogous to Mz per-site normalization).
3. Transform to eigenbasis: `V' * PP_avg_comp * V` where V is `trotter.eigvecs` if trotter is provided, else `hamiltonian.eigvecs`.
4. Push the matrix with name "XX_avg", "YY_avg", or "ZZ_avg" respectively.

The final returned observables should be: ["H", "Mz", "XX_avg", "YY_avg", "ZZ_avg"] (5 total).

Update the docstring to reflect the new observables and remove the "ZZ correlations are intentionally excluded" note on line 123.

Implementation pattern -- add after line 141 (`observables = vcat([H], mz_obs)`), before the return:

```julia
# Averaged two-site correlations (nearest-neighbor, periodic)
V = trotter !== nothing ? trotter.eigvecs : hamiltonian.eigvecs
for (pauli_pair, pair_name) in [([X, X], "XX_avg"), ([Y, Y], "YY_avg"), ([Z, Z], "ZZ_avg")]
    dim = size(hamiltonian.data, 1)
    PP_sum = zeros(ComplexF64, dim, dim)
    for i in 1:num_qubits
        PP_sum .+= Matrix{ComplexF64}(pad_term(pauli_pair, num_qubits, i; periodic=true))
    end
    PP_sum ./= num_qubits  # Per-bond average
    PP_eigen = Matrix{ComplexF64}(V' * PP_sum * V)
    push!(observables, PP_eigen)
    push!(names, pair_name)
end
```

Note: Change the observables/names construction from `vcat` to mutable `push!` or adjust accordingly so that pushing the correlation observables works. The simplest approach: initialize `observables = Matrix{ComplexF64}[H]` and `names = String["H"]`, then `append!(observables, mz_obs)` / `append!(names, mz_names)`, then the loop above pushes to them.
  </action>
  <verify>
Run: `cd /Users/bence/code/QuantumFurnace.jl && julia --project=. -e 'using QuantumFurnace, LinearAlgebra; h = HamHam([[X,X],[Y,Y],[Z,Z]], [1.0,1.0,1.0], 4, 5.0; periodic=true); obs, names = build_gap_estimation_observables(h, 4); println("n_obs=", length(obs), " names=", names); @assert length(obs) == 5; @assert names == ["H","Mz","XX_avg","YY_avg","ZZ_avg"]; println("OK")'`
  </verify>
  <done>build_gap_estimation_observables returns 5 observables: H, Mz, XX_avg, YY_avg, ZZ_avg. Each correlation observable is a per-bond average over all nearest-neighbor bonds, correctly transformed to eigenbasis.</done>
</task>

<task type="auto">
  <name>Task 2: Update tests for 5 observables and create validation script</name>
  <files>test/test_gap_estimation.jl, experiments/run_gap_validation.jl</files>
  <action>
**Part A: Update test/test_gap_estimation.jl**

1. In "Basic estimate_spectral_gap returns SpectralGapResult" testset (line 20-39):
   - Change `@test result.best_observable in ["H", "Mz"]` to `@test result.best_observable in ["H", "Mz", "XX_avg", "YY_avg", "ZZ_avg"]`
   - Change `@test length(result.per_observable) == 2` to `@test length(result.per_observable) == 5`
   - Change `@test length(result.observable_names) == 2` to `@test length(result.observable_names) == 5`
   - Change `@test result.observable_names == ["H", "Mz"]` to `@test result.observable_names == ["H", "Mz", "XX_avg", "YY_avg", "ZZ_avg"]`

2. In "Custom observables with names" testset (line 93-104):
   - Change `@test length(result.per_observable) == length(obs)` -- this is already generic, keep as-is.
   - Verify `length(obs)` is now 5 by adding: `@test length(obs) == 5`

3. In "Deterministic seeding produces identical results" testset (line 60-74):
   - The existing assertions for per_observable[1] and [2] are fine. Optionally add a check that all 5 per_observable gaps match: loop over all indices.

**Part B: Create experiments/run_gap_validation.jl**

Create a standalone script that:
1. Builds a Heisenberg XX chain for n=6, beta=5.0 (same as run_sweep.jl)
2. Builds TrotterDomain system (same constants: DELTA=0.01, T0, NUM_ENERGY_BITS=12, W0=0.05, etc.)
3. Runs `estimate_spectral_gap` with ntraj=10000 (high trajectory count), save_every=10, seed=42, skip_initial=0.1
4. Computes exact gap via `build_liouvillian` + `LindbladianResult`
5. Calls `cross_validate_gap` to get the factor
6. Prints results: fitted gap, exact gap, relative error, best observable, factor (fitted/exact)
7. Also prints per-observable fit summary (name, gap, R-squared, converged)

Use the same helper functions from run_sweep.jl (build_heisenberg_hamiltonian, build_trotter_system) -- duplicate them or call them directly.

Script structure:
```julia
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
    hamiltonian = HamHam([[X,X],[Y,Y],[Z,Z]], [1.0,1.0,1.0], n, beta; periodic=true)
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

    # Config
    config = ThermalizeConfig(
        num_qubits=n, with_coherent=true, with_linear_combination=true,
        domain=TrotterDomain(), beta=beta, sigma=1.0/beta,
        a=beta/30.0, b=0.4,
        num_energy_bits=NUM_ENERGY_BITS, w0=W0, t0=T0,
        num_trotter_steps_per_t0=NUM_TROTTER_STEPS_PER_T0,
        mixing_time=mixing_time, delta=DELTA,
    )

    # Initial state
    dim = 2^n
    psi0 = zeros(ComplexF64, dim)
    psi0[1] = 1.0

    # Estimate spectral gap (uses all 5 observables now)
    println("Running estimate_spectral_gap with ntraj=$ntraj...")
    t0 = time()
    result = estimate_spectral_gap(
        jumps, config, psi0, hamiltonian;
        ntraj=ntraj, save_every=10, seed=SEED,
        trotter=trotter, skip_initial=0.1,
    )
    wall = time() - t0
    @printf("Wall time: %.1fs\n\n", wall)

    # Exact gap
    println("Computing exact Liouvillian gap...")
    exact_result = build_liouvillian(jumps, config, hamiltonian; trotter=trotter)

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
```

This script is self-contained and can be run with:
`julia --project=. experiments/run_gap_validation.jl`
  </action>
  <verify>
1. Run tests: `cd /Users/bence/code/QuantumFurnace.jl && julia --project=. -e 'using Pkg; Pkg.test()'` -- all tests pass, including updated gap estimation tests.
2. Verify validation script parses: `cd /Users/bence/code/QuantumFurnace.jl && julia --project=. -e 'include("experiments/run_gap_validation.jl")' 2>&1 | head -5` -- prints header line and starts running (may take several minutes for n=6 with 10k trajectories; just verify it starts).
  </verify>
  <done>Tests updated to expect 5 observables. Validation script created and runnable. All existing tests pass with the updated observable count.</done>
</task>

</tasks>

<verification>
1. `build_gap_estimation_observables(ham, n)` returns 5 observables named ["H", "Mz", "XX_avg", "YY_avg", "ZZ_avg"]
2. All existing tests in test/test_gap_estimation.jl pass with updated assertions
3. The validation script `experiments/run_gap_validation.jl` runs and produces gap comparison output
4. The trotter kwarg path also produces 5 observables (same function, same code path)
</verification>

<success_criteria>
- build_gap_estimation_observables returns 5 observables in both Hamiltonian and Trotter eigenbases
- All tests pass
- Validation script is ready for user to run and evaluate n=6 improvement
</success_criteria>

<output>
After completion, create `.planning/quick/24-add-two-site-correlations-to-spectral-ga/24-SUMMARY.md`
</output>
