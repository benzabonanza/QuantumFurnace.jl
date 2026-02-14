---
phase: quick-5
plan: 01
type: execute
wave: 1
depends_on: []
files_modified:
  - test/test_dm_scaling.jl
autonomous: true
must_haves:
  truths:
    - "TrotterDomain error thresholds are tight, not loose 0.1 placeholders"
    - "All three Trotter sanity checks (DMTST-05, DMTST-06, DMTST-06b) still pass"
  artifacts:
    - path: "test/test_dm_scaling.jl"
      provides: "Tightened Trotter error thresholds"
      contains: "1e-5"
  key_links: []
---

<objective>
Tighten the three loose `< 0.1` TrotterDomain error thresholds in test/test_dm_scaling.jl
to values that are tight but still safe given measured errors.

Purpose: The 0.1 thresholds were temporary safety margins introduced during quick fixes 1-4.
Actual measured errors are far smaller, so these thresholds should be tightened for safety
(catching regressions that increase Trotter error).

Output: Updated test/test_dm_scaling.jl with tight thresholds.
</objective>

<execution_context>
@./.claude/get-shit-done/workflows/execute-plan.md
@./.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@test/test_dm_scaling.jl
@test/test_helpers.jl
</context>

<tasks>

<task type="auto">
  <name>Task 1: Tighten TrotterDomain error thresholds in test_dm_scaling.jl</name>
  <files>test/test_dm_scaling.jl</files>
  <action>
Replace three `< 0.1` Trotter sanity thresholds in test/test_dm_scaling.jl with tight values
based on measured errors:

1. **Line 132** (DMTST-05: coherent term B consistency):
   - Current: `@test dist_bohr_trott < 0.1`
   - Measured error: ~0.0106
   - Change to: `@test dist_bohr_trott < 0.02`
   - Update comment to: `# Trotter error on B term (measured ~0.011, threshold 2x margin)`
   - NOTE: Cannot use 1e-5 here because the actual B-term Trotter error is ~0.011
     (this is a real physical approximation error from the Trotter decomposition, not a bug).

2. **Line 184** (DMTST-06: OFT consistency):
   - Current: `@test dist_energy_trott < 0.1`
   - Measured error: ~1.47e-8
   - Change to: `@test dist_energy_trott < 1e-5`
   - Update comment to: `# Trotter OFT error (measured ~1.5e-8, tight threshold)`

3. **Line 264** (DMTST-06b: NUFFT OFT consistency):
   - Current: `@test dist_nufft_trott_vs_energy < 0.1`
   - Measured error: ~1.47e-8
   - Change to: `@test dist_nufft_trott_vs_energy < 1e-5`
   - Update comment to: `# NUFFT Trotter OFT error (measured ~1.5e-8, tight threshold)`

Do NOT change any other thresholds (TOL_QUADRATURE, 1e-10, 1e-12 etc. are already tight).
Do NOT change the delta sweep values `[0.1, 0.05, 0.025, 0.0125]` in DMTST-03/04 --
those are step sizes for scaling tests, not error thresholds.
  </action>
  <verify>
Run the full test suite: `cd /Users/bence/code/QuantumFurnace.jl && julia --project -e 'include("test/test_helpers.jl"); include("test/test_dm_scaling.jl")'`

All tests must pass. Specifically verify:
- DMTST-05 passes with the tighter 0.02 threshold
- DMTST-06 passes with 1e-5 threshold
- DMTST-06b passes with 1e-5 threshold
  </verify>
  <done>
All three TrotterDomain `< 0.1` thresholds replaced:
- DMTST-05: 0.1 -> 0.02 (tight but safe for ~0.011 measured error)
- DMTST-06: 0.1 -> 1e-5 (tight but safe for ~1.5e-8 measured error)
- DMTST-06b: 0.1 -> 1e-5 (tight but safe for ~1.5e-8 measured error)
All tests pass.
  </done>
</task>

</tasks>

<verification>
Run the full DM scaling test suite and confirm all 12 tests pass (DMTST-03 through DMTST-06b).
No threshold should be looser than 10x the measured error.
</verification>

<success_criteria>
- Zero occurrences of `< 0.1` as an error threshold in test_dm_scaling.jl
  (the delta sweep `[0.1, ...]` values are step sizes, not thresholds)
- All tests pass with tightened thresholds
- DMTST-05 threshold is 0.02 (cannot go to 1e-5 due to real Trotter approximation error)
- DMTST-06 and DMTST-06b thresholds are 1e-5
</success_criteria>

<output>
After completion, create `.planning/quick/5-tighten-trotterdomain-error-thresholds-f/5-SUMMARY.md`
</output>
