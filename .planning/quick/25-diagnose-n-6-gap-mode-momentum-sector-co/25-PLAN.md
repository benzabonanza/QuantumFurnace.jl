---
phase: 25-diagnose-n-6-gap-mode-momentum-sector
plan: 01
type: execute
wave: 1
depends_on: []
files_modified:
  - experiments/diagnose_gap_momentum.jl
autonomous: true
must_haves:
  truths:
    - "Script reports the momentum quantum number k of the gap mode (v_2) for n=4 and n=6"
    - "Script verifies T_L commutes with the Lindbladian (confirming translational invariance)"
    - "Script reports whether all 5 observables live in k=0 sector (confirming their translational invariance)"
    - "Output confirms or refutes the momentum sector hypothesis for n=6 zero overlap"
  artifacts:
    - path: "experiments/diagnose_gap_momentum.jl"
      provides: "Momentum sector diagnostic script"
      min_lines: 80
  key_links:
    - from: "experiments/diagnose_gap_momentum.jl"
      to: "src/hamiltonian.jl"
      via: "HamHam constructor for periodic Heisenberg chain"
      pattern: "HamHam"
    - from: "experiments/diagnose_gap_momentum.jl"
      to: "src/qi_tools.jl"
      via: "pad_term for constructing translation operator"
      pattern: "pad_term"
---

<objective>
Diagnose why n=6 periodic Heisenberg chain has zero overlap between all observables and the Lindbladian gap mode, while n=4 works fine.

Purpose: The hypothesis is that the Lindbladian inherits translational symmetry from the periodic chain, so eigenmodes have definite crystal momentum k = 2*pi*m/n. All 5 observables (H, Mz, XX_avg, YY_avg, ZZ_avg) are translationally invariant (k=0 sector). If the gap mode (second eigenvalue) lives in a nonzero momentum sector, all k=0 observables get exactly zero overlap, explaining the n=6 failure.

Output: A diagnostic script that constructs the translation operator in Liouville space, measures the momentum of the gap eigenvector, and reports which sector it belongs to for n=4 and n=6.
</objective>

<execution_context>
@./.claude/get-shit-done/workflows/execute-plan.md
@./.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@.planning/STATE.md
@experiments/validate_spectral_gap.jl
@src/hamiltonian.jl
@src/qi_tools.jl
@src/gap_estimation.jl
@src/convergence.jl
@src/constants.jl
@src/misc_tools.jl
</context>

<tasks>

<task type="auto">
  <name>Task 1: Create momentum sector diagnostic script</name>
  <files>experiments/diagnose_gap_momentum.jl</files>
  <action>
Create `experiments/diagnose_gap_momentum.jl` that performs the following analysis for n=4 and n=6:

**1. System setup (reuse pattern from validate_spectral_gap.jl):**
- Use `HamHam([[X,X],[Y,Y],[Z,Z]], [1.0,1.0,1.0], n, beta; periodic=true)` with beta=10.0
- Build jump operators and LiouvConfig exactly as in the existing validation script
- Build the full dense Lindbladian via `run_lindbladian` then `Matrix(liouv_result.liouvillian)`

**2. Construct the translation operator T in Hilbert space:**
- T cyclically permutes qubits: T|q_1, q_2, ..., q_n> = |q_2, q_3, ..., q_n, q_1>
- Build T as a permutation matrix of dimension 2^n x 2^n
- For each computational basis state |b_1 b_2 ... b_n> (integer j in 0..2^n-1), compute the permuted state index by cyclically shifting the binary representation left by one position: new_index = ((j << 1) | (j >> (n-1))) & (2^n - 1)
- Note: Julia uses 1-based indexing, so add 1 to both old and new indices
- Verify: T^n = I (T to the n-th power equals identity)
- Verify: T * H_comp * T' = H_comp (T commutes with computational-basis Hamiltonian)

**3. Construct T_L in Liouville space:**
- T_L = T tensor conj(T), acting on vectorized density matrices vec(rho) = (rho_11, rho_21, ..., rho_d1, rho_12, ...)
- Using the column-stacking (Watrous) convention that the codebase uses: T_L = kron(T, conj(T))
- Verify: T_L commutes with L (the Lindbladian matrix). Check ||[T_L, L]|| / ||L|| < 1e-10. Print this commutator norm.

**4. Eigendecompose L and measure gap mode momentum:**
- Compute `F = eigen(Matrix(L))`, sort by `abs.(real.(F.values))`
- The gap mode is the second eigenvector v_2 (sorted by smallest |Re(lambda)|)
- Compute `T_L * v_2` and extract the eigenvalue: `exp_ik = dot(v_2, T_L * v_2)` (since v_2 is an eigenvector of T_L if T_L commutes with L)
- Extract momentum: k = angle(exp_ik), momentum_index m = round(k * n / (2*pi))
- Print: eigenvalue of T_L on v_2, |eigenvalue| (should be ~1.0), momentum k, sector index m
- Also report T_L eigenvalues for the first ~10 eigenmodes to see the full momentum structure

**5. Verify observable momentum sectors:**
- Build the 5 observables via `build_preset_trajectory_observables(ham, n)` -- these are in the eigenbasis
- For each observable O, vectorize: O_vec = vec(O) (reshape to column vector)
- Compute T_L * O_vec vs O_vec: if O is translationally invariant, T_L * O_vec = O_vec (eigenvalue = 1, i.e. k=0)
- BUT: the observables from build_preset_trajectory_observables are in the Hamiltonian eigenbasis. The Lindbladian L is also in eigenbasis. So T must first be transformed: T_eigen = V' * T * V where V = ham.eigvecs. Then T_L_eigen = kron(T_eigen, conj(T_eigen)). Use this T_L_eigen for both the commutator check and the momentum measurement.
- Print for each observable: `<O_vec | T_L | O_vec> / <O_vec | O_vec>`, which should be 1.0 for k=0 observables
- Actually simpler: just check dot(O_vec, T_L_eigen * O_vec) / dot(O_vec, O_vec) for the Rayleigh quotient

**6. Summary:**
- Print a clear table: for each n, report:
  - Commutator norm ||[T_L, L]|| / ||L||
  - Gap eigenvalue lambda_2
  - Gap mode T_L eigenvalue (complex)
  - Gap mode momentum k and sector m
  - Whether gap mode is in k=0 sector
  - Observable momentum check results
- Print conclusion: "CONFIRMED: gap mode is in k={m} sector, not k=0 -- explains zero overlap" or "REFUTED: gap mode is in k=0 sector"

**Important implementation details:**
- Use `using QuantumFurnace` and `using LinearAlgebra`
- Copy the `make_system`, `make_liouv_config` helpers from validate_spectral_gap.jl (same parameters)
- The script goes in experiments/ (gitignored directory, but script can be force-added)
- For the binary permutation: Julia bit operations are `<<`, `>>`, `&`, `|` on integers. Computational basis state j (0-indexed) maps to ((j << 1) | (j >> (n-1))) & ((1 << n) - 1)
- Handle potential eigenvalue degeneracy near the gap: report eigenvalues for the first 10 modes and their momentum sectors, to see if there are near-degenerate modes in different sectors
  </action>
  <verify>
Run: `cd /Users/bence/code/QuantumFurnace.jl && julia --project experiments/diagnose_gap_momentum.jl`

Expected output includes:
- T^n = I verification passes for both n=4 and n=6
- Commutator norm ||[T_L, L]||/||L|| is < 1e-10 for both system sizes
- Gap mode momentum sector is reported for both n=4 and n=6
- All 5 observables confirmed as k=0 (translational eigenvalue ~1.0)
- Clear conclusion about whether momentum sector explains n=6 zero overlap
  </verify>
  <done>
Script runs successfully for n=4 and n=6, reports the momentum quantum number of the gap mode for each system size, confirms observables are in k=0 sector, and provides a definitive confirmation or refutation of the momentum sector hypothesis for the n=6 zero-overlap phenomenon.
  </done>
</task>

</tasks>

<verification>
- Script executes without error for both n=4 and n=6
- Translation operator T satisfies T^n = I
- T_L commutes with L (commutator norm < 1e-10)
- Momentum of gap mode is clearly reported
- Observable k=0 membership is verified
- Conclusion is stated explicitly
</verification>

<success_criteria>
1. The momentum quantum number of the gap mode is determined for both n=4 and n=6
2. The relationship between gap mode momentum and observable overlap is explained
3. The hypothesis (gap mode in non-zero momentum sector for n=6) is confirmed or refuted with numerical evidence
</success_criteria>

<output>
After completion, create `.planning/quick/25-diagnose-n-6-gap-mode-momentum-sector-co/25-SUMMARY.md`
</output>
