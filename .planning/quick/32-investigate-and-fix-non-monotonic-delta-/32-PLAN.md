---
phase: quick-32
plan: 01
type: execute
wave: 1
depends_on: []
files_modified:
  - experiments/investigate_delta_scaling_bug.jl
  - src/trajectories.jl  # potential fix if bug found
  - src/gap_estimation.jl  # potential fix if bug found
autonomous: true
must_haves:
  truths:
    - "Root cause of non-monotonic delta scaling identified and documented"
    - "Multi-step trajectory observable averages compared against exact exp(t*L) evolution for delta=0.1, 0.01, 0.001"
    - "If a code bug is found, it is fixed and existing tests still pass"
  artifacts:
    - path: "experiments/investigate_delta_scaling_bug.jl"
      provides: "Diagnostic script isolating DM-vs-trajectory divergence at multiple deltas"
  key_links:
    - from: "experiments/investigate_delta_scaling_bug.jl"
      to: "src/trajectories.jl"
      via: "step_along_trajectory! and run_observable_trajectories"
      pattern: "step_along_trajectory!|run_observable_trajectories"
---

<objective>
Investigate and diagnose the non-monotonic delta scaling bug where smaller Trotter step
sizes (delta) produce WORSE spectral gap estimation errors -- opposite to expected O(delta^2)
Trotter convergence. Specifically: Mz_stagg error goes +0.6% -> -6.5% -> -19.7% as delta
shrinks 0.1 -> 0.01 -> 0.001 (Quick-31 data).

Purpose: Determine whether this is (a) a code bug in the trajectory or gap estimation
pipeline, (b) a fundamental property of the discrete Kraus channel interacting with
the fitting procedure, or (c) an error in how multi-step trajectories accumulate
Lie-Trotter splitting bias.

Output: Diagnostic script with results, identification of root cause, and fix if applicable.
</objective>

<context>
@.planning/STATE.md
@src/trajectories.jl
@src/gap_estimation.jl
@src/fitting.jl
@src/jump_workers.jl
@src/kraus.jl
@src/convergence.jl
@test/test_regression.jl
@test/trajectory_validation/run_trajectory_validation.jl
@experiments/validate_gap_delta_scaling_v2.jl
</context>

<tasks>

<task type="auto">
  <name>Task 1: Create diagnostic script comparing multi-step trajectory rho against exact exp(t*L) at 3 deltas</name>
  <files>experiments/investigate_delta_scaling_bug.jl</files>
  <action>
Create a Julia script that isolates exactly where the DM-trajectory divergence happens
at different delta values. This is the KEY diagnostic: previous Quick-30/31 only looked
at gap estimation output, not at the underlying density matrix accuracy.

The script should:

1. **Setup:** Use the same 3-qubit SMALL system from test_helpers.jl (n=3, dim=8) to keep
   computation tractable. Build the Lindbladian L for TimeDomain with with_coherent=false.
   Use uniform superposition psi0 = ones(dim)/sqrt(dim).

2. **For each delta in [0.1, 0.01, 0.001]:**

   a. **Exact DM evolution via exp(t*L):** Compute rho_exact(t) = exp(t*L) * vec(rho0)
      for T_total = 5.0 (fixed total time). This is the ground truth -- continuous
      Lindbladian evolution with no Trotter splitting.

   b. **Multi-step trajectory average:** Run 50,000 trajectories with the SAME T_total
      and the given delta. Each trajectory does ceil(T_total/delta) steps. Average
      the final density matrices to get rho_traj(T_total).

   c. **DM Kraus evolution (deterministic):** Run the DM _jump_contribution! step
      (the density matrix version of the CPTP channel) for ceil(T_total/delta) steps,
      applying ALL operators each step (not random selection). This is the "expected"
      result of the trajectory channel -- the channel the trajectories are unraveling.
      Actually, the DM code applies a RANDOM jump each step just like trajectories.
      Instead, construct the CPTP superoperator for a single step:
        E(rho) = (1/n_jumps) * sum_a [K0_a rho K0_a' + U_res_a rho U_res_a' + delta * sum_w rate(w) * L_{a,w} rho L_{a,w}']
      Then E^N(rho0) = the exact N-step channel evolution.

      SIMPLIFICATION: Since computing the full superoperator E is complex, instead
      compare trajectory rho directly against exp(T*L) rho0. The key question is:
      does ||rho_traj - rho_exact|| scale as O(delta) or O(delta^2) with total time fixed?

      With T_total fixed and num_steps = T/delta:
      - Single-step splitting error: O(delta^2) per step (trajectory channel vs exp(delta*L))
      - N steps accumulate to: N * O(delta^2) = (T/delta) * O(delta^2) = O(delta) total error

      So we EXPECT trace_distance(rho_traj, rho_exact) to scale as O(delta).
      If it does NOT decrease monotonically with delta, that signals a bug.

   d. **Measure:** trace_distance(rho_traj, rho_exact) for each delta.
      Also measure observable expectation values: for each observable O,
      compute tr(rho_traj * O) vs tr(rho_exact * O). This directly tests
      whether the trajectory-averaged observable matches the exact Lindbladian
      result -- independent of the fitting procedure.

3. **Observable-level comparison:** For each of the 8 preset observables, compute:
   - exact_val = tr(rho_exact * O)
   - traj_val = tr(rho_traj * O)
   - error = traj_val - exact_val
   Do this for all 3 deltas. If |error| decreases monotonically with delta for each
   observable, the non-monotonic gap estimation error is in the FITTING procedure
   (exponential decay fitting), not in the trajectory simulation.

4. **Time-resolved observable comparison:** Additionally, for delta=0.01, compute
   the trajectory-averaged observable time series O(t) at save_every=10 and compare
   against the exact time series from iterated exp(delta*save_every*L). If these
   match well, the non-monotonic behavior is definitively in the fitting, not the
   simulation.

5. **Single-step sanity check:** As a control, verify that single-step trajectory
   rho matches exp(delta*L)*rho0 to O(delta^2) for all 3 deltas (replicating
   TVAL-02/03 but with explicit delta sweep including delta=0.001).

Print all results in clearly labeled sections. Use `@printf` for formatted output.

IMPORTANT NOTES:
- Use the 4-qubit system (n=4) from load_hamiltonian("heis", 4; beta=10.0) to match
  the exact same system from Quick-30/31. The 3-qubit system would be faster but
  wouldn't reproduce the reported bug.
- For the multi-step trajectory part, use run_observable_trajectories with
  reconstruct_dm=true to get both the density matrix AND observable time series.
- The script should be self-contained with its own helper functions (like the v2 script).
- Use seed=42 for reproducibility.
- T_total=5.0 should be enough mixing time for n=4 (gap ~0.17, so 5/gap ~ 29 steps
  of size 1, which is about 5 mixing times).
- Actually use T_total = 20.0 to match Quick-31 mixing_time.

The script structure should be:
  Section 0: Constants and system setup
  Section 1: Single-step DM-vs-trajectory sanity check (3 deltas)
  Section 2: Multi-step rho comparison: rho_traj vs exp(T*L)*rho0 (3 deltas)
  Section 3: Observable expectation values from final rho (3 deltas, 8 observables)
  Section 4: Time-resolved observable comparison for delta=0.01
  Section 5: Summary and diagnosis
  </action>
  <verify>
Run: `cd /Users/bence/code/QuantumFurnace.jl && julia -t 4 --project experiments/investigate_delta_scaling_bug.jl`
Script completes without errors and prints results for all 5 sections.
  </verify>
  <done>
- Script produces data showing whether trace_distance(rho_traj, rho_exact) decreases
  monotonically with delta (expected O(delta) for multi-step)
- Observable expectation value errors are computed for all 3 deltas and 8 observables
- Root cause is identified: either (a) trajectory rho is wrong at small delta (code bug),
  or (b) trajectory rho is correct but gap fitting produces non-monotonic errors
  (fitting/estimation issue, not a bug)
  </done>
</task>

<task type="auto">
  <name>Task 2: Analyze results and fix any bugs found</name>
  <files>src/trajectories.jl, src/gap_estimation.jl, src/fitting.jl</files>
  <action>
Based on the diagnostic results from Task 1:

**If trajectory rho is correct (trace_distance decreases monotonically with delta):**
- The non-monotonic gap estimation is NOT a code bug but a property of the
  exponential decay fitting procedure interacting with the discrete Kraus channel
  spectral properties at different step sizes.
- Document this finding clearly in the summary.
- No code changes needed. The "bug" is actually expected behavior: the fitting
  model (A*exp(-gap*t) + C) may not capture the true multi-exponential decay
  structure of the discrete Kraus channel, and different delta values excite
  different spectral modes.

**If trajectory rho is wrong (trace_distance does NOT decrease monotonically):**
- There is a genuine code bug. Investigate the most likely suspects:

  a. **num_steps rounding:** `ceil(Int, total_time / delta)` -- for delta=0.001 and
     total_time=20.0, this gives 20000 steps. Verify this is correct and that the
     trajectory actually runs this many steps.

  b. **Normalization drift:** After 20000 steps at delta=0.001, numerical errors in
     the normalization `rmul!(psi, 1/sqrt(norm2))` could accumulate. Check if the
     final psi norm deviates from 1.0 significantly.

  c. **scaled_prefactor computation:** The `scaled_prefactor` in TrajectoryFramework
     involves `1/p_jump` rescaling. Verify this is independent of delta (it should be
     -- it depends only on n_jumps and physical constants).

  d. **p_jump_total = delta * expR:** This uses bare delta. With R already scaled
     by n_jumps, verify that the probability branching is correct. At very small
     delta, p_jump_total becomes very small, making almost all steps "no-jump".
     This is correct behavior but could amplify numerical issues.

  e. **U_residual PSD clamping:** At small delta, S = (2*alpha - delta)*R - alpha^2*R^2
     with alpha = 1-sqrt(1-delta) ~ delta/2. For small delta, S ~ delta^2 * (R - R^2/4).
     The eigenvalue clamping could introduce O(delta^2) artifacts that accumulate over
     T/delta = O(1/delta) steps.

- Apply the fix to the relevant source file.
- Run `cd /Users/bence/code/QuantumFurnace.jl && julia --project -e 'using Pkg; Pkg.test()'`
  to verify all existing tests still pass.

**In either case:** Update the diagnosis in the summary with the root cause determination.
  </action>
  <verify>
If fix applied: `julia --project -e 'using Pkg; Pkg.test()'` passes all tests.
If no fix needed: Document that trajectory rho scales correctly and the non-monotonic
behavior is in the gap estimation fitting, not the simulation.
  </verify>
  <done>
Root cause of non-monotonic delta scaling is definitively identified as either:
(a) A code bug that has been fixed and all tests pass, OR
(b) An inherent property of the fitting procedure documented with evidence
  </done>
</task>

</tasks>

<verification>
- experiments/investigate_delta_scaling_bug.jl exists and runs to completion
- Single-step trajectory matches exp(delta*L) to O(delta^2) for all 3 deltas
- Multi-step comparison data available for all 3 deltas
- Root cause identified with supporting numerical evidence
- If fix applied, all existing tests pass
</verification>

<success_criteria>
- The diagnostic script produces clear, quantitative evidence about whether the
  non-monotonic delta scaling originates in the trajectory simulation or the fitting
- Root cause is identified with a confidence level supported by numerical data
- If a bug is found, it is fixed; if not, the finding is documented
</success_criteria>

<output>
After completion, create `.planning/quick/32-investigate-and-fix-non-monotonic-delta-/32-SUMMARY.md`
</output>
