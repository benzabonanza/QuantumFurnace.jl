---
phase: quick-31
plan: 01
type: execute
wave: 1
depends_on: []
files_modified:
  - experiments/validate_gap_delta_scaling_v2.jl
autonomous: true
must_haves:
  truths:
    - "Script runs with 4 threads and completes delta-scaling analysis for all 3 deltas"
    - "mixing_time is fixed at 20 (not computed from gap)"
    - "skip_initial is 0.3 (not 0.1)"
    - "Initial state psi0 is uniform superposition ones(ComplexF64, dim)/sqrt(dim)"
    - "All analysis sections (delta-scaling, per-observable, Richardson) produce output"
  artifacts:
    - path: "experiments/validate_gap_delta_scaling_v2.jl"
      provides: "Revalidation script with longer mixing, higher skip, uniform psi0"
      contains: "mixing_time = 20"
  key_links:
    - from: "experiments/validate_gap_delta_scaling_v2.jl"
      to: "QuantumFurnace module"
      via: "using QuantumFurnace"
      pattern: "estimate_spectral_gap"
---

<objective>
Create validate_gap_delta_scaling_v2.jl -- a variant of the delta-scaling validation script with longer mixing time (20), higher skip_initial (0.3), and uniform superposition initial state. This tests whether the systematic observable bias identified in Quick-30 is reduced when the system has more time to mix and the initial state covers all Hilbert space subspaces.

Purpose: Quick-30 showed gap error is NOT O(delta) and Richardson extrapolation is ineffective. The hypothesis is that the previous script's mixing_time (computed from gap, ~5/gap) and excited-state psi0 may not have been in the correct asymptotic regime for exponential gap convergence.

Output: New script at experiments/validate_gap_delta_scaling_v2.jl
</objective>

<execution_context>
@./.claude/get-shit-done/workflows/execute-plan.md
@./.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@.planning/STATE.md
@experiments/validate_gap_delta_scaling.jl
</context>

<tasks>

<task type="auto">
  <name>Task 1: Create v2 delta-scaling validation script with modified parameters</name>
  <files>experiments/validate_gap_delta_scaling_v2.jl</files>
  <action>
Copy experiments/validate_gap_delta_scaling.jl to experiments/validate_gap_delta_scaling_v2.jl, then make these specific changes:

1. **Header comment updates:**
   - Update title to "Delta-Scaling Revalidation (v2: longer mixing, uniform psi0)"
   - Update Parameters section:
     - `mixing_time = 20` (fixed, not computed from gap)
     - `skip_initial = 0.3` (was 0.1)
     - `Initial state: psi0 = ones(ComplexF64, dim) / sqrt(dim) (uniform superposition)`
   - Update Usage line: `julia -t 4 --project experiments/validate_gap_delta_scaling_v2.jl`

2. **Initial state (around line 167-168):**
   Replace:
   ```julia
   psi0 = zeros(ComplexF64, dim)
   psi0[end] = 1.0
   ```
   With:
   ```julia
   psi0 = ones(ComplexF64, dim) / sqrt(dim)
   ```
   Update the comment from "excited state (locked decision)" to "uniform superposition (all subspaces covered)"

3. **Mixing time (around line 171-172):**
   Replace:
   ```julia
   mixing_time = max(5.0 / exact_gap, 10.0)
   @printf("Mixing time: %.1f (5/gap=%.1f)\n\n", mixing_time, 5.0 / exact_gap)
   ```
   With:
   ```julia
   mixing_time = 20.0
   @printf("Mixing time: %.1f (fixed)\n\n", mixing_time)
   ```

4. **skip_initial (around line 193):**
   Replace `skip_initial=0.1` with `skip_initial=0.3` in the `estimate_spectral_gap` call.

5. **Banner update (around line 139-143):**
   Update the banner to mention "v2" and the changed parameters.

Everything else (DELTAS, NTRAJ, SEED, SAVE_EVERY, n=4, observables, all analysis sections 3-5) stays IDENTICAL to the v1 script.
  </action>
  <verify>
Run: `diff experiments/validate_gap_delta_scaling.jl experiments/validate_gap_delta_scaling_v2.jl` to confirm only the expected lines changed:
- Header comments updated
- psi0 changed to uniform superposition
- mixing_time changed to fixed 20.0
- skip_initial changed from 0.1 to 0.3
- Banner updated
- Everything else identical
  </verify>
  <done>
experiments/validate_gap_delta_scaling_v2.jl exists with exactly 4 parameter changes (mixing_time=20, skip_initial=0.3, uniform psi0, updated comments). Diff shows no other logic changes.
  </done>
</task>

<task type="auto">
  <name>Task 2: Run the v2 script with 4 threads and capture output</name>
  <files>experiments/validate_gap_delta_scaling_v2.jl</files>
  <action>
Run the script:
```bash
cd /Users/bence/code/QuantumFurnace.jl && julia -t 4 --project experiments/validate_gap_delta_scaling_v2.jl
```

Capture full output. The script will:
1. Load the n=4 disordered Heisenberg Hamiltonian
2. Compute exact gap via ARPACK
3. Run 20k trajectories for each of delta=0.1, 0.01, 0.001
4. Print per-observable fits, delta-scaling analysis, Richardson extrapolation, and summary

Key metrics to record in the SUMMARY:
- Exact gap value
- Estimated gap at each delta (and which observable was best)
- error/delta ratios (are they more constant than Quick-30's 96x spread?)
- Richardson extrapolation improvement factor
- Whether O(delta) scaling is now confirmed or still not
- Comparison with Quick-30 results (excited state, mixing_time~5/gap, skip_initial=0.1)
  </action>
  <verify>
Script completes without errors. Output contains all 5 sections:
1. System setup and exact gap
2. Per-delta trajectory results
3. Delta-scaling analysis with error/delta table
4. Richardson extrapolation results
5. Final summary with O(delta) verdict
  </verify>
  <done>
Script runs to completion with 4 threads. All analysis sections produce numeric results. SUMMARY.md documents the key findings and comparison with Quick-30.
  </done>
</task>

</tasks>

<verification>
- experiments/validate_gap_delta_scaling_v2.jl exists and differs from v1 only in the 4 specified parameters
- Script output contains complete delta-scaling analysis for all 3 deltas
- SUMMARY.md documents whether longer mixing + uniform psi0 changes the O(delta) scaling conclusion
</verification>

<success_criteria>
- v2 script created with exactly the specified parameter changes
- Script runs successfully with `julia -t 4`
- Results documented with clear comparison to Quick-30 findings
- Verdict on whether longer mixing time regime shows O(delta) scaling
</success_criteria>

<output>
After completion, create `.planning/quick/31-delta-scaling-revalidation-with-longer-m/31-SUMMARY.md`
</output>
