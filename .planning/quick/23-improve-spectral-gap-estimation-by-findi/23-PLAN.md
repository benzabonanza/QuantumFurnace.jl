---
phase: quick-23
plan: 01
type: execute
wave: 1
depends_on: []
files_modified:
  - src/gap_estimation.jl
  - test/test_gap_estimation.jl
  - experiments/validate_gap_estimation.jl
  - experiments/diagnose_observable_overlap.jl
autonomous: true
must_haves:
  truths:
    - "_select_best_observable picks the smallest gap among good fits, not the highest R-squared"
    - "Validation script passes with tighter factor range [0.8, 1.5] for both n=4 and n=6"
    - "All existing tests still pass; new selection test verifies smallest-gap logic"
  artifacts:
    - path: "src/gap_estimation.jl"
      provides: "Corrected _select_best_observable with smallest-gap-among-good-fits criterion"
      contains: "smallest gap"
    - path: "test/test_gap_estimation.jl"
      provides: "New test case for smallest-gap selection with two valid fits"
      contains: "smallest gap"
    - path: "experiments/validate_gap_estimation.jl"
      provides: "Updated pass criterion and comments reflecting smallest-gap selection"
      contains: "0.8 <= residual_factor <= 1.5"
  key_links:
    - from: "src/gap_estimation.jl"
      to: "test/test_gap_estimation.jl"
      via: "_select_best_observable selection logic"
      pattern: "_select_best_observable"
---

<objective>
Fix spectral gap estimation by changing observable selection from highest-R-squared to smallest-gap-among-good-fits.

Purpose: The current selection picks observables with the highest R-squared, but M_z has zero overlap with the gap mode for the Heisenberg XXX chain (it commutes with H). M_z's high R-squared fit is to noise, giving a meaningless gap (~1.7x overestimate). The physically correct approach is to pick the smallest positive gap among fits with acceptable R-squared, since the true spectral gap is the slowest-decaying mode.

Output: Corrected `_select_best_observable`, updated tests, updated validation script, removed diagnostic script.
</objective>

<execution_context>
@./.claude/get-shit-done/workflows/execute-plan.md
@./.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@.planning/STATE.md
@src/gap_estimation.jl
@test/test_gap_estimation.jl
@experiments/validate_gap_estimation.jl
</context>

<tasks>

<task type="auto">
  <name>Task 1: Fix _select_best_observable and add test</name>
  <files>src/gap_estimation.jl, test/test_gap_estimation.jl</files>
  <action>
  **In `src/gap_estimation.jl`:**

  1. Remove the TODO comment at line 56.

  2. Replace the `_select_best_observable` function (lines 57-86) with new logic:
     - Update the docstring to describe the new selection criterion: among fits where (converged AND gap > 0 AND R-squared > 0.8), select the one with the **smallest gap** (not highest R-squared). The rationale is that the spectral gap is the slowest-decaying mode, so the smallest positive fitted gap from a good fit is closest to the true gap.
     - Primary criterion: filter fits where `fit.converged && fit.gap > 0.0 && fit.r_squared > 0.8`, then pick the one with the smallest `fit.gap`.
     - Fallback 1: if no fits pass the R-squared > 0.8 filter, among fits where `fit.converged && fit.gap > 0.0`, pick the one with the smallest `fit.gap`.
     - Fallback 2: if still none, pick the fit with the highest `fit.r_squared` overall (existing fallback behavior, for diagnostic purposes).
     - Return `(best_idx, names[best_idx], fits[best_idx].r_squared)` as before.

  **In `test/test_gap_estimation.jl`:**

  3. In the `"Selection logic: _select_best_observable"` testset (starting line 116), add a new test case AFTER the existing tests (do not remove existing tests). The new test creates two FitResult objects that are BOTH valid (converged=true, gap > 0, R-squared > 0.8) but with different gaps and R-squared values:
     - `fit_high_r2_high_gap`: converged=true, gap=0.8, R-squared=0.98
     - `fit_lower_r2_small_gap`: converged=true, gap=0.3, R-squared=0.92
     - Call `_select_best_observable` with both fits and verify it picks `fit_lower_r2_small_gap` (the one with the smaller gap), NOT `fit_high_r2_high_gap` (the one with higher R-squared).
     - This directly tests that the new "smallest gap among good fits" logic works.

  4. Update the existing test at line 149-157 ("Test selection: should pick fit_converged_good") -- the comment says "should pick fit_converged_good ... because it's the only one that is converged AND has gap > 0". This test still works with the new logic because fit_converged_good is the ONLY valid fit, so it's trivially the smallest gap among valid fits. Just update the comment to mention the new criterion. The assertion (`idx == 1`) stays correct.

  5. Run `julia --project=. -e 'using Pkg; Pkg.test()'` or the specific test to confirm all tests pass.
  </action>
  <verify>
  Run: `cd /Users/bence/code/QuantumFurnace.jl && julia --project=. -e 'using Test; include("test/runtests.jl")'`
  All tests pass, including existing selection logic tests and the new smallest-gap test.
  </verify>
  <done>
  - `_select_best_observable` selects smallest gap among good fits (converged, gap > 0, R-squared > 0.8)
  - TODO comment at line 56 is removed
  - New test case verifies smallest-gap selection when two valid fits exist
  - All existing tests still pass
  </done>
</task>

<task type="auto">
  <name>Task 2: Update validation script and delete diagnostic</name>
  <files>experiments/validate_gap_estimation.jl, experiments/diagnose_observable_overlap.jl</files>
  <action>
  **In `experiments/validate_gap_estimation.jl`:**

  1. Update the header comment block (lines 1-20):
     - Remove the explanation about "discrete-step Kraus decomposition effects" causing ~1.5-1.7x factor.
     - Replace with: the smallest-gap selection criterion picks the observable whose fitted decay rate best approximates the true Lindbladian spectral gap, yielding a residual factor of ~1.0-1.1x.
     - Update pass criterion comment from `[1.0, 3.0]` to `[0.8, 1.5]`.

  2. Update the residual factor analysis comments (around lines 233-237):
     - Remove: "The discrete-step Kraus decomposition introduces a factor between the fitted trajectory rate and the continuous Liouvillian gap."
     - Replace with explanation that smallest-gap selection picks the observable with the best physical overlap with the gap mode, so the residual factor should be close to 1.0.
     - Change the expected factor from "~1.5-1.7" to "~1.0-1.1" in the printf on line 241.

  3. Update the pass criterion (line 248):
     - Change `factor_in_range = 1.0 <= residual_factor <= 3.0` to `factor_in_range = 0.8 <= residual_factor <= 1.5`
     - Update the printf on line 253 from `[1.0, 3.0]` to `[0.8, 1.5]`.

  **Delete diagnostic script:**

  4. Delete `experiments/diagnose_observable_overlap.jl` using `rm` (it was a temporary diagnostic for this investigation).
  </action>
  <verify>
  Run: `cd /Users/bence/code/QuantumFurnace.jl && julia --project=. experiments/validate_gap_estimation.jl`
  Both n=4 and n=6 should show PASS with residual factors near 1.0 (in the [0.8, 1.5] range).
  The diagnostic file should no longer exist: `ls experiments/diagnose_observable_overlap.jl` returns not found.
  </verify>
  <done>
  - Validation script pass criterion is [0.8, 1.5]
  - Comments accurately describe smallest-gap selection (no more "Kraus decomposition effects")
  - Expected factor comments say ~1.0-1.1
  - Both n=4 and n=6 pass the validation
  - `experiments/diagnose_observable_overlap.jl` is deleted
  </done>
</task>

</tasks>

<verification>
1. All unit tests pass: `julia --project=. -e 'using Test; include("test/runtests.jl")'`
2. Validation script passes for both system sizes: `julia --project=. experiments/validate_gap_estimation.jl` shows OVERALL: PASS
3. No reference to deleted diagnostic file remains in tracked code
</verification>

<success_criteria>
- `_select_best_observable` picks smallest gap among good fits (R-squared > 0.8 threshold)
- New test case proves smallest-gap selection beats highest-R-squared selection
- Validation residual factors are ~1.0-1.1x (down from ~1.6x)
- Both n=4 and n=6 pass the tightened [0.8, 1.5] criterion
- All existing tests continue to pass
- Diagnostic script removed
</success_criteria>

<output>
After completion, create `.planning/quick/23-improve-spectral-gap-estimation-by-findi/23-SUMMARY.md`
</output>
