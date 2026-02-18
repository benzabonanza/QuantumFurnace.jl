---
phase: quick-30
plan: 01
type: execute
wave: 1
depends_on: []
files_modified:
  - experiments/validate_gap_delta_scaling.jl
autonomous: true
must_haves:
  truths:
    - "Script runs n=4 disordered Heisenberg for delta = 0.1, 0.01, 0.001 using same observables"
    - "Error/delta ratio is reported for each delta value to assess O(delta) scaling"
    - "Richardson extrapolation is computed for two pairs and compared to exact gap"
  artifacts:
    - path: "experiments/validate_gap_delta_scaling.jl"
      provides: "Delta-scaling validation and Richardson extrapolation script"
      min_lines: 150
  key_links:
    - from: "experiments/validate_gap_delta_scaling.jl"
      to: "estimate_spectral_gap"
      via: "delta keyword override per run"
      pattern: "estimate_spectral_gap.*delta="
---

<objective>
Validate that spectral gap estimation error scales as O(delta) with Trotter step size, then test Richardson extrapolation to improve estimates.

Purpose: The discrete-step Kraus effect introduces systematic bias in trajectory gap estimates. If error = C*delta (linear in delta), Richardson extrapolation can cancel the leading error term and yield much more accurate gap estimates without needing tiny delta values.

Output: experiments/validate_gap_delta_scaling.jl script with results printed to stdout
</objective>

<execution_context>
@./.claude/get-shit-done/workflows/execute-plan.md
@./.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@.planning/STATE.md
@experiments/validate_gap_disordered.jl
</context>

<tasks>

<task type="auto">
  <name>Task 1: Create delta-scaling validation script with Richardson extrapolation</name>
  <files>experiments/validate_gap_delta_scaling.jl</files>
  <action>
Create `experiments/validate_gap_delta_scaling.jl` following the exact pattern of `experiments/validate_gap_disordered.jl`.

**Section 0: Constants and Setup**
- BETA = 10.0, NTRAJ = 20_000, SAVE_EVERY = 10, SEED = 42
- DELTAS = [0.1, 0.01, 0.001] -- three delta values to test
- Same grid parameters: NUM_ENERGY_BITS=12, W0=0.05, T0=2pi/(2^12*0.05), NUM_TROTTER_STEPS_PER_T0=10, SIGMA=1/BETA
- Reuse `make_system(n, beta)` helper exactly from validate_gap_disordered.jl (load_hamiltonian, pad_term jump construction)
- Reuse `make_liouv_config(n)` helper exactly
- Reuse `make_thermalize_config(n; mixing_time)` helper exactly (use delta=0.01 as config default -- actual delta overridden per-run)

**Section 1: n=4 system setup and exact gap**
- n = 4 only (locked decision)
- ham, jumps, dim = make_system(4, BETA)
- Build Lindbladian via run_lindbladian for exact gap
- exact_gap = abs(real(liouv_result.spectral_gap))
- Print exact gap
- Build observables via build_preset_trajectory_observables(ham, 4) -- ONCE, reuse for all 3 delta runs
- psi0 = zeros(ComplexF64, dim); psi0[end] = 1.0 (excited state, locked decision)

**Section 2: Run trajectory gap estimation for each delta**
- mixing_time = max(5.0/exact_gap, 10.0) -- computed once
- config_t = make_thermalize_config(4; mixing_time=mixing_time) -- one config
- For each delta in DELTAS:
  - Call estimate_spectral_gap(jumps, config_t, psi0, ham; ntraj=NTRAJ, save_every=SAVE_EVERY, seed=SEED, skip_initial=0.1, delta=delta)
    NOTE: The `delta` keyword overrides config.delta per-run
  - Store gap_result in a Dict or Vector indexed by delta
  - Print per-observable table: name, gap, R2, converged, CI
  - Print: estimated_gap (best), exact_gap, error = estimated - exact, abs_error, error/delta ratio

**Section 3: Delta-scaling analysis**
- Print table: delta | estimated_gap | exact_gap | error | error/delta
- Check if error/delta is approximately constant (within factor of 2 across the 3 delta values)
- Print conclusion: "O(delta) scaling CONFIRMED" or "O(delta) scaling NOT confirmed"

Also analyze per-observable: For each of the 8 observables that converged across all 3 delta values, print error/delta. This lets us see if individual observables also show O(delta) scaling.

**Section 4: Richardson extrapolation**
- Use the formula for O(h) error: gap_rich = (h2 * gap(h1) - h1 * gap(h2)) / (h2 - h1)
  where h1 < h2 (h1 is the finer delta)
- Pair 1: h1=0.01, h2=0.1 (using best-observable gap from each run)
- Pair 2: h1=0.001, h2=0.01
- For each pair print: gap_rich, exact_gap, error, relative_error
- Also do per-observable Richardson for observables converged in both runs of the pair

**Section 5: Summary**
- Print final summary table with all results
- Print whether O(delta) scaling was observed
- Print Richardson improvement factor (error_richardson / error_fine_delta for each pair)

IMPORTANT: Use the SAME seed=42 for all runs so noise is comparable. Use the SAME observables object for all 3 runs (built once from build_preset_trajectory_observables).

IMPORTANT: For the "best observable" selection across delta values, use the best_observable from each individual gap_result (it may differ per delta). But also track a FIXED observable (the one selected as best at delta=0.01) and report its error/delta across all 3 delta values for a clean comparison.
  </action>
  <verify>
Run: `cd /Users/bence/code/QuantumFurnace.jl && julia --project experiments/validate_gap_delta_scaling.jl`

Expected:
- Script completes without error
- Three gap estimates produced (one per delta)
- error/delta ratio printed for each delta
- Richardson extrapolation results printed for 2 pairs
- All numerical values are reasonable (gaps positive, errors decrease with delta)
  </verify>
  <done>
- Script runs successfully for n=4 disordered Heisenberg with delta = 0.1, 0.01, 0.001
- error/delta ratio is reported and assessed for O(delta) scaling
- Richardson extrapolation computed for both pairs (0.1/0.01 and 0.01/0.001)
- Richardson error vs exact gap reported with improvement factor
  </done>
</task>

<task type="auto">
  <name>Task 2: Write summary documenting delta-scaling and Richardson results</name>
  <files>.planning/quick/30-validate-trotter-delta-order-gap-error-a/30-SUMMARY.md</files>
  <action>
After running the script, create the summary at `.planning/quick/30-validate-trotter-delta-order-gap-error-a/30-SUMMARY.md` following the standard summary template.

Include:
- Whether O(delta) error scaling was confirmed (error/delta approximately constant)
- The error/delta ratio values for each delta
- Richardson extrapolation results (extrapolated gap, error, improvement factor)
- Per-observable analysis highlights
- Key numerical results table
- Any surprises or deviations from expected behavior
  </action>
  <verify>File exists and contains numerical results from the script run.</verify>
  <done>Summary written with delta-scaling confirmation/denial and Richardson extrapolation results.</done>
</task>

</tasks>

<verification>
- experiments/validate_gap_delta_scaling.jl exists and runs without error
- Output shows gap estimates for 3 delta values with error/delta analysis
- Richardson extrapolation results shown for 2 pairs
- Summary documents findings
</verification>

<success_criteria>
- Script produces numerical evidence for/against O(delta) gap estimation error scaling
- Richardson extrapolation tested and improvement factor quantified
- Results documented in summary
</success_criteria>

<output>
After completion, create `.planning/quick/30-validate-trotter-delta-order-gap-error-a/30-SUMMARY.md`
</output>
