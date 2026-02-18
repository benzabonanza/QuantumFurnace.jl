---
phase: quick-29
plan: 01
type: execute
wave: 1
depends_on: []
files_modified: [experiments/validate_gap_xx_stagg.jl]
autonomous: true

must_haves:
  truths:
    - "XX_stagg observable is manually constructed as Sum((-1)^i * X_i * X_{i+1}) / n with periodic boundary"
    - "ONLY H and XX_stagg are passed as observables (not the full 8-observable bundle)"
    - "Script runs with julia --threads=3 --project for both n=4 and n=6"
    - "Eigenbasis overlap analysis shows |c_gap| for both H and XX_stagg"
    - "Full expansion coefficient spectrum printed for XX_stagg (all modes, not just gap mode)"
    - "Eigenvalue spectrum near gap (modes 1-5) printed"
    - "Near-degenerate mode detection with eigenvalues close to the gap"
    - "Per-observable gap/exact ratios printed when >10% error"
    - "Raw time series values at early/middle/late time points printed"
  artifacts:
    - path: "experiments/validate_gap_xx_stagg.jl"
      provides: "Rigorous XX_stagg disordered Heisenberg gap validation with deep diagnostics"
      min_lines: 250
  key_links:
    - from: "experiments/validate_gap_xx_stagg.jl"
      to: "src/gap_estimation.jl"
      via: "estimate_spectral_gap with custom observables=[H_eigen, XX_stagg_eigen]"
      pattern: "estimate_spectral_gap.*observables=\\[H"
    - from: "experiments/validate_gap_xx_stagg.jl"
      to: "src/gap_estimation.jl"
      via: "eigenbasis_overlap_analysis with custom observables"
      pattern: "eigenbasis_overlap_analysis.*L_dense.*\\[H"
---

<objective>
Create and run a focused validation script that tests XX_stagg (staggered nearest-neighbor XX correlation) as a new observable for disordered Heisenberg gap estimation.

Purpose: Quick-28 showed disorder breaks n=6 symmetry protection but gap estimation was biased at 34-49% error using the full 8-observable bundle. XX_stagg (Sum((-1)^i * X_i * X_{i+1}) / n) is a new observable not yet in the codebase that combines staggered sign (k=pi momentum) with two-site XX correlation. This script isolates just H + XX_stagg to measure overlap and gap estimation quality, with deep diagnostics if error exceeds 10%.

Output: experiments/validate_gap_xx_stagg.jl script + execution results captured in SUMMARY
</objective>

<execution_context>
@./.claude/get-shit-done/workflows/execute-plan.md
@./.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@.planning/STATE.md
@.planning/quick/28-test-gap-mode-coupling-with-random-field/28-SUMMARY.md
@experiments/validate_gap_disordered.jl
@src/convergence.jl
@src/gap_estimation.jl
</context>

<tasks>

<task type="auto">
  <name>Task 1: Create XX_stagg validation script with deep diagnostics</name>
  <files>experiments/validate_gap_xx_stagg.jl</files>
  <action>
Create experiments/validate_gap_xx_stagg.jl following the structure of experiments/validate_gap_disordered.jl but with these critical differences:

**Observable construction (MANUALLY build XX_stagg, do NOT use build_preset_trajectory_observables):**

```julia
# XX_stagg = Sum((-1)^i * X_i * X_{i+1}) / n  (staggered nearest-neighbor XX, periodic)
XX_stagg_comp = zeros(ComplexF64, dim, dim)
for i in 1:num_qubits
    sign = (-1)^i
    XX_stagg_comp .+= sign .* Matrix{ComplexF64}(pad_term([X, X], num_qubits, i; periodic=true))
end
XX_stagg_comp ./= num_qubits
XX_stagg_eigen = Matrix{ComplexF64}(V' * XX_stagg_comp * V)
```

And H in eigenbasis:
```julia
H_eigen = Matrix{ComplexF64}(diagm(ComplexF64.(ham.eigvals)))
```

**Pass ONLY these two observables to both eigenbasis_overlap_analysis and estimate_spectral_gap:**
```julia
custom_obs = [H_eigen, XX_stagg_eigen]
custom_names = ["H", "XX_stagg"]
```

**Parameters (same as Quick-28, locked):**
- beta=10.0, delta=0.01, ntraj=20_000, save_every=10, seed=42
- TimeDomain, with_coherent=false, excited initial state (psi0[end]=1.0)
- System sizes: n=4, n=6
- load_hamiltonian("heis", n; beta=10.0) for disordered Heisenberg

**Script sections (for each n in [4, 6]):**

Section 0: Constants (same as validate_gap_disordered.jl)

Section 1a: System setup via make_system(n, BETA) -- reuse helper from validate_gap_disordered.jl pattern

Section 1b: Exact gap via run_lindbladian + ARPACK

Section 1c: Build custom observables H_eigen and XX_stagg_eigen manually

Section 1d: Eigenbasis overlap analysis
- Call eigenbasis_overlap_analysis(L_dense, custom_obs, custom_names, rho0)
- Print |c_gap| and relative overlap for both H and XX_stagg
- Print whether symmetry is broken (|c_gap| > 0.001)

Section 1e: **Full expansion coefficient spectrum for XX_stagg**
- Access overlap.overlap_coefficients (n_obs x n_modes matrix)
- Print |c_k| for ALL modes k=1..min(20, n_modes) for XX_stagg (row 2)
- Identify which modes have the largest |c_k|, rank top 5

Section 1f: **Eigenvalue spectrum near gap (modes 1-5)**
- Print eigenvalues[1:5] from overlap.eigenvalues (real and imaginary parts)
- Flag near-degenerate modes: any pair with |Re(lambda_i) - Re(lambda_j)| < 0.1 * exact_gap

Section 1g: Trajectory gap estimation
- Call estimate_spectral_gap with custom_obs, custom_names
- Print per-observable fit details (gap, R2, converged, CI)

Section 1h: **Deep diagnostics (ALWAYS print, not just on >10% error)**
- Per-observable gap/exact ratios for both H and XX_stagg
- Raw time series values: access the ObservableTrajectoryResult to get measurements_mean at time indices [1, div(end,4), div(end,2), div(3*end,4), end] for XX_stagg
  NOTE: estimate_spectral_gap does NOT return the raw time series. To get raw time series, call run_observable_trajectories SEPARATELY before or after estimate_spectral_gap:
  ```julia
  traj_result = run_observable_trajectories(jumps, config_t, psi0, ham;
      observables=custom_obs, save_every=SAVE_EVERY,
      ntraj=NTRAJ, total_time=config_t.mixing_time, delta=DELTA,
      seed=SEED)
  ```
  Then access traj_result.measurements_mean (n_obs x n_times matrix) and traj_result.times.
  Print 5 sample time points with their XX_stagg values.
- If >10% error: print additional "INVESTIGATION" section with analysis of whether the dominant mode in XX_stagg's expansion is actually the gap mode or a different mode

Section 2: Final summary table

**Usage line in header:**
```
# Usage:
#   cd QuantumFurnace.jl && julia --threads=3 --project experiments/validate_gap_xx_stagg.jl
```
  </action>
  <verify>
Verify file exists and has the right structure:
- `grep -c "XX_stagg" experiments/validate_gap_xx_stagg.jl` shows 20+ occurrences
- `grep "pad_term\(\[X, X\]" experiments/validate_gap_xx_stagg.jl` finds the manual construction
- `grep "(-1)\^i" experiments/validate_gap_xx_stagg.jl` finds staggered sign
- `grep "observables=\[H_eigen" experiments/validate_gap_xx_stagg.jl` OR `grep "observables=custom_obs" experiments/validate_gap_xx_stagg.jl` confirms custom 2-observable pass
- `grep "build_preset_trajectory_observables" experiments/validate_gap_xx_stagg.jl` returns NOTHING (must NOT use preset bundle)
- `grep "overlap_coefficients" experiments/validate_gap_xx_stagg.jl` finds expansion spectrum analysis
- `grep "run_observable_trajectories" experiments/validate_gap_xx_stagg.jl` finds raw time series call
- File is 250+ lines
  </verify>
  <done>
experiments/validate_gap_xx_stagg.jl exists with:
- Manual XX_stagg construction using staggered (-1)^i * X_i * X_{i+1} / n
- Only H + XX_stagg passed as observables (no preset bundle)
- Eigenbasis overlap analysis with |c_gap| output
- Full expansion coefficient spectrum for XX_stagg (modes 1-20)
- Eigenvalue spectrum near gap (modes 1-5) with degeneracy check
- Gap estimation via estimate_spectral_gap with custom observables
- Raw time series sampling via run_observable_trajectories
- Per-observable gap/exact ratios
- Deep investigation section for >10% error cases
  </done>
</task>

<task type="auto">
  <name>Task 2: Run validation script and capture results</name>
  <files></files>
  <action>
Run the script with 3 threads:
```bash
cd /Users/bence/code/QuantumFurnace.jl && julia --threads=3 --project experiments/validate_gap_xx_stagg.jl
```

Use a 10-minute timeout (600000ms). The script should complete in ~3-5 minutes based on Quick-28 timing (22min for 8 observables, this has only 2).

Capture the FULL output. Key things to extract for the SUMMARY:
1. XX_stagg |c_gap| values for n=4 and n=6
2. Whether XX_stagg has non-zero gap mode overlap (|c_gap| > 0.001)
3. Full expansion coefficient spectrum -- which modes dominate?
4. Eigenvalue near-degeneracy findings
5. Gap estimation results: estimated vs exact gap, relative error
6. Per-observable gap/exact ratios
7. Raw time series behavior at sample points
8. If >10% error: what the deep diagnostics reveal about WHY

If the script fails with an error:
- Read the error carefully
- Fix the script (likely typo or API mismatch)
- Re-run

Record all findings for the SUMMARY.
  </action>
  <verify>
Script ran to completion (printed "FINAL SUMMARY" section). Both n=4 and n=6 results captured.
  </verify>
  <done>
Script output captured with: XX_stagg overlap values, expansion spectrum, eigenvalue analysis, gap estimates, raw time series, and diagnostic findings for both n=4 and n=6.
  </done>
</task>

</tasks>

<verification>
- experiments/validate_gap_xx_stagg.jl exists and is 250+ lines
- Script uses ONLY H + XX_stagg (no preset bundle)
- XX_stagg is manually constructed with staggered sign and periodic XX bonds
- Script ran successfully for both n=4 and n=6 with --threads=3
- Deep diagnostics printed (expansion spectrum, eigenvalue analysis, raw time series)
- Results captured and analyzed
</verification>

<success_criteria>
- XX_stagg |c_gap| values measured for both n=4 and n=6 disordered Heisenberg
- If |c_gap| > 0.001: gap estimation results with relative error quantified
- If >10% error despite overlap: diagnostic output explains WHY (dominant mode analysis, degeneracy, time series behavior)
- All findings recorded in SUMMARY for scientific understanding
</success_criteria>

<output>
After completion, create `.planning/quick/29-rigorous-disordered-heisenberg-gap-estim/29-SUMMARY.md`
</output>
