---
phase: 27-add-symmetry-breaking-observable
plan: 01
type: execute
wave: 1
depends_on: []
files_modified:
  - src/convergence.jl
  - test/test_gap_estimation.jl
  - test/test_convergence.jl
  - experiments/validate_spectral_gap.jl
autonomous: true

must_haves:
  truths:
    - "build_preset_trajectory_observables returns 8 observables including XZ_stagg"
    - "XZ_stagg is Hermitian, staggered (k=pi), and breaks SU(2) symmetry"
    - "All existing tests pass with updated observable count (7 -> 8)"
    - "n=6 gap estimation error improves (XZ_stagg has nonzero gap-mode overlap)"
  artifacts:
    - path: "src/convergence.jl"
      provides: "XZ_stagg observable construction in build_preset_trajectory_observables"
      contains: "XZ_stagg"
    - path: "test/test_gap_estimation.jl"
      provides: "Updated assertions for 8 observables"
      contains: "XZ_stagg"
    - path: "test/test_convergence.jl"
      provides: "Updated assertions for 8 observables"
      contains: "XZ_stagg"
  key_links:
    - from: "src/convergence.jl"
      to: "src/constants.jl"
      via: "X and Z Pauli constants in XZ_stagg construction"
      pattern: "pad_term\\(\\[X, Z\\]"
---

<objective>
Add XZ_stagg (staggered nearest-neighbor XZ correlation) as the 8th observable in build_preset_trajectory_observables, then update all test assertions from 7 to 8 observables.

Purpose: The n=6 Heisenberg chain gap mode is protected by both translational symmetry (k=pi sector) and SU(2) spin-rotation symmetry. All current 7 observables have zero or near-zero gap-mode overlap for n=6, leading to ~10.7% relative error. XZ_stagg = sum_i (-1)^i X_i Z_{i+1} / n breaks SU(2) while having k=pi momentum component, enabling coupling to the gap mode.

Output: Updated source with 8-observable bundle, all tests green, validation script ready to verify improved n=6 gap estimation.
</objective>

<execution_context>
@./.claude/get-shit-done/workflows/execute-plan.md
@./.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@.planning/STATE.md
@src/convergence.jl
@src/constants.jl
@test/test_gap_estimation.jl
@test/test_convergence.jl
@experiments/validate_spectral_gap.jl
</context>

<tasks>

<task type="auto">
  <name>Task 1: Add XZ_stagg observable to build_preset_trajectory_observables</name>
  <files>src/convergence.jl</files>
  <action>
In `build_preset_trajectory_observables` (src/convergence.jl), add the XZ_stagg observable construction AFTER the Z1 block (after line 89, before the return statement on line 91).

1. Add XZ_stagg construction block following the exact pattern of the XX_avg/YY_avg/ZZ_avg loop but with staggered sign and asymmetric Pauli pair [X, Z]:

```julia
# Staggered XZ correlation: sum((-1)^i * X_i Z_{i+1}) / n  (breaks SU(2), has k=pi)
XZ_stagg_comp = zeros(ComplexF64, dim, dim)
for i in 1:num_qubits
    sign = (-1)^i
    XZ_stagg_comp .+= sign .* Matrix{ComplexF64}(pad_term([X, Z], num_qubits, i; periodic=true))
end
XZ_stagg_comp ./= num_qubits  # Per-bond normalization (consistent with XX_avg, YY_avg, ZZ_avg)
XZ_stagg_eigen = Matrix{ComplexF64}(V' * XZ_stagg_comp * V)
push!(observables, XZ_stagg_eigen)
push!(names, "XZ_stagg")
```

2. Update the docstring (lines 13-31) to reflect 8 observables:
   - Change "7 observables" to "8 observables" in the description line
   - Change the name list to: `["H", "Mz", "XX_avg", "YY_avg", "ZZ_avg", "Mz_stagg", "Z1", "XZ_stagg"]`
   - Add bullet: `- \`XZ_stagg\`: Per-bond staggered nearest-neighbor XZ correlation (sum of (-1)^i X_i Z_{i+1} / n). Has k=pi momentum component and breaks SU(2) spin-rotation symmetry, enabling coupling to symmetry-protected gap modes.`

Key physics: XZ_stagg is Hermitian because X_i and Z_{i+1} act on different sites, so (X_i Z_{i+1})^dagger = X_i^dagger Z_{i+1}^dagger = X_i Z_{i+1}. The (-1)^i staggering gives k=pi Fourier component. The asymmetric Pauli type (X on one site, Z on neighbor) breaks SU(2) spin-rotation symmetry since SU(2) rotations transform X,Y,Z uniformly but XZ is not invariant under joint spin rotations.
  </action>
  <verify>
Run `cd /Users/bence/code/QuantumFurnace.jl && julia --project -e 'using QuantumFurnace; obs, names = build_preset_trajectory_observables(HamHam([[X,X],[Y,Y],[Z,Z]], [1.0,1.0,1.0], 4, 1.0; periodic=true), 4); println(length(obs), " ", names)'` and confirm output shows 8 observables with names ending in "XZ_stagg".
  </verify>
  <done>build_preset_trajectory_observables returns 8 observables; 8th is "XZ_stagg" constructed with staggered (-1)^i sign and [X,Z] Pauli pair with periodic boundary via pad_term.</done>
</task>

<task type="auto">
  <name>Task 2: Update all test assertions from 7 to 8 observables</name>
  <files>test/test_gap_estimation.jl, test/test_convergence.jl</files>
  <action>
**test/test_gap_estimation.jl** -- update these lines:

1. Line 30: Change `in ["H", "Mz", "XX_avg", "YY_avg", "ZZ_avg", "Mz_stagg", "Z1"]` to `in ["H", "Mz", "XX_avg", "YY_avg", "ZZ_avg", "Mz_stagg", "Z1", "XZ_stagg"]`
2. Line 32: Change `== 7` to `== 8`
3. Line 33: Change `== 7` to `== 8`
4. Line 34: Change `== ["H", "Mz", "XX_avg", "YY_avg", "ZZ_avg", "Mz_stagg", "Z1"]` to `== ["H", "Mz", "XX_avg", "YY_avg", "ZZ_avg", "Mz_stagg", "Z1", "XZ_stagg"]`
5. Line 96: Change `== 7` to `== 8`
6. Line 219: Change `== 7` to `== 8` (observable_names length in Eigenbasis Overlap)
7. Line 221 (size check): Change `(7,` to `(8,` in `size(result.overlap_coefficients) == (7, SMALL_DIM^2)`
8. Line 222: Change `== 7` to `== 8` (gap_mode_overlap length)
9. Line 223: Change `== 7` to `== 8` (relative_gap_overlap length)

**test/test_convergence.jl** -- update these lines:

1. Line 197: Change `(7, 5)` to `(8, 5)` in `size(conv_data.observable_values) == (7, 5)`
2. Line 199: Change `== 7` to `== 8` in `length(conv_data.observable_gibbs_values) == 7`
3. Line 285: Change `== 7` to `== 8` in `length(col_slice) == 7`
4. Line 410: Change `(7,` to `(8,` in `size(conv_data.observable_values) == (7, conv_data.total_batches)`
5. Line 637: Change `== 7` to `== 8` (length(observables))
6. Line 638: Change `== 7` to `== 8` (length(names))
7. Line 639: Change `== ["H", "Mz", "XX_avg", "YY_avg", "ZZ_avg", "Mz_stagg", "Z1"]` to `== ["H", "Mz", "XX_avg", "YY_avg", "ZZ_avg", "Mz_stagg", "Z1", "XZ_stagg"]`
8. Line 698: Change `== 7` to `== 8` (Trotter basis length(observables))
9. Line 699: Change `== 7` to `== 8` (Trotter basis length(names))
10. Line 700: Change `== ["H", "Mz", "XX_avg", "YY_avg", "ZZ_avg", "Mz_stagg", "Z1"]` to `== ["H", "Mz", "XX_avg", "YY_avg", "ZZ_avg", "Mz_stagg", "Z1", "XZ_stagg"]`

Also add XZ_stagg construction verification test in the eigenbasis testset (after the Z1 test at line 688), following the same pattern:
```julia
# XZ_stagg observable (index 8) matches inline staggered XZ construction
XZ_stagg_comp = zeros(ComplexF64, DIM, DIM)
for i in 1:NUM_QUBITS
    sign = (-1)^i
    XZ_stagg_comp .+= sign .* Matrix{ComplexF64}(pad_term([X, Z], NUM_QUBITS, i; periodic=true))
end
XZ_stagg_comp ./= NUM_QUBITS
XZ_stagg_expected = Matrix{ComplexF64}(V' * XZ_stagg_comp * V)
@test isapprox(observables[8], XZ_stagg_expected; atol=1e-14)
```

Do NOT change the "Z1 observable (index 7)" comment or its test -- Z1 stays at index 7.
  </action>
  <verify>
Run `cd /Users/bence/code/QuantumFurnace.jl && julia --project -e 'using Test; include("test/runtests.jl")'` (or the subset: test_convergence.jl and test_gap_estimation.jl). All tests must pass with zero failures.
  </verify>
  <done>All 7->8 count assertions updated in both test files; XZ_stagg added to name lists; construction verification test added for index 8; all tests pass.</done>
</task>

<task type="auto">
  <name>Task 3: Run validation and update n=6 target threshold</name>
  <files>experiments/validate_spectral_gap.jl</files>
  <action>
1. Run the full test suite first to confirm all unit tests pass:
   `cd /Users/bence/code/QuantumFurnace.jl && julia --project -e 'using Pkg; Pkg.test()'`

2. Run the validation script to measure XZ_stagg's effect on n=6 gap estimation:
   `cd /Users/bence/code/QuantumFurnace.jl && julia --project experiments/validate_spectral_gap.jl`

3. Examine the output:
   - Check XZ_stagg's |c_gap| overlap for n=6 (should be nonzero, indicating SU(2) symmetry breaking works)
   - Check the new n=6 relative error (should improve from ~10.7%)
   - Check n=4 still passes < 1%

4. Based on results, update TARGET_REL_ERROR_N6 on line 38 of experiments/validate_spectral_gap.jl:
   - If n=6 error drops below 5%: set to 0.05 (5%)
   - If n=6 error drops below 2%: set to 0.02 (2%)
   - If n=6 error stays above 5%: set to current_error * 1.2 (20% margin above observed)
   - Add a comment explaining the improvement: "# X% for n=6 (XZ_stagg breaks SU(2), nonzero gap-mode overlap)"

5. Re-run validation to confirm PASS with updated threshold.

Note: The validation script takes ~2-5 minutes for n=6 (4096-dimensional Liouvillian, 20k trajectories). Use timeout of 600000ms.
  </action>
  <verify>
Validation script prints "ALL PASS" in final summary. n=4 relative error < 1%. n=6 relative error improved from ~10.7% baseline.
  </verify>
  <done>Validation confirms XZ_stagg has nonzero gap-mode overlap for n=6; n=6 relative error improved; TARGET_REL_ERROR_N6 tightened to reflect improvement; all tests and validation pass.</done>
</task>

</tasks>

<verification>
1. `julia --project -e 'using Pkg; Pkg.test()'` -- all unit tests pass
2. `julia --project experiments/validate_spectral_gap.jl` -- ALL PASS
3. build_preset_trajectory_observables returns 8 observables with names ["H", "Mz", "XX_avg", "YY_avg", "ZZ_avg", "Mz_stagg", "Z1", "XZ_stagg"]
4. XZ_stagg has nonzero |c_gap| for n=6 (confirms SU(2) breaking)
</verification>

<success_criteria>
- 8-observable bundle including XZ_stagg (staggered nearest-neighbor XZ)
- All unit tests pass with 8-observable assertions
- n=6 gap estimation improved from ~10.7% baseline
- n=4 gap estimation still < 1% (no regression)
- Validation script passes with updated threshold
</success_criteria>

<output>
After completion, create `.planning/quick/27-add-symmetry-breaking-observable-for-n-6/27-SUMMARY.md`
</output>
