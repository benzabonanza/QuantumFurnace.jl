---
phase: 38-fix-failing-diagnostics-and-filename-tes
plan: 01
type: execute
wave: 1
depends_on: []
files_modified:
  - test/test_results.jl
  - test/test_diagnostics.jl
autonomous: true

must_haves:
  truths:
    - "test_results.jl filename test passes with kms_ prefix"
    - "test_diagnostics.jl fixed point distance < 0.01 for BohrDomain"
    - "test_diagnostics.jl Trotter/Bohr ratio test passes"
    - "test_diagnostics.jl backward compat test passes"
  artifacts:
    - path: "test/test_results.jl"
      provides: "Filename generation test with KMS construction"
    - path: "test/test_diagnostics.jl"
      provides: "Diagnostics tests with KMS construction at all 4 call sites"
  key_links:
    - from: "make_small_thermalize_config(TrotterDomain())"
      to: "construction=KMS()"
      via: "keyword argument"
    - from: "make_small_liouv_config(BohrDomain()) / make_small_liouv_config(TrotterDomain())"
      to: "construction=KMS()"
      via: "keyword argument"
---

<objective>
Fix 5 failing tests by passing construction=KMS() at the 5 call sites that currently default to GNS().

Purpose: These tests verify the diagnostics module behavior, not the construction method. They require the Gibbs state to be the exact fixed point, which only KMS (exact detailed balance) provides. GNS (approximate) yields trace_distance≈0.035, failing the < 0.01 threshold.
Output: All 5 tests pass without changing thresholds or test logic.
</objective>

<execution_context>
@./.claude/get-shit-done/workflows/execute-plan.md
@./.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@.planning/PROJECT.md
@.planning/ROADMAP.md
@.planning/STATE.md
</context>

<tasks>

<task type="auto">
  <name>Task 1: Fix test_results.jl filename test (issue #1)</name>
  <files>test/test_results.jl</files>
  <action>
At line 282, change:
  `kms_config = make_small_thermalize_config(TrotterDomain())`
to:
  `kms_config = make_small_thermalize_config(TrotterDomain(); construction=KMS())`

This is the only change needed in this file. The test at line 284 checks `startswith(kms_filename, "kms_")` which requires KMS construction to generate the kms_ prefix.
  </action>
  <verify>Run: `cd /Users/bence/code/QuantumFurnace.jl && julia --project=. -e 'using Pkg; Pkg.test(; test_args=["results"])'` — confirm the filename generation testset passes.</verify>
  <done>test_results.jl:284 passes: `@test startswith(kms_filename, "kms_")` is true.</done>
</task>

<task type="auto">
  <name>Task 2: Fix test_diagnostics.jl at all 4 call sites (issues #2-5)</name>
  <files>test/test_diagnostics.jl</files>
  <action>
Make the following 4 targeted changes:

1. Line 4 (top-level setup, fixes issues #2 and #3):
   `config = make_small_liouv_config(BohrDomain())`
   → `config = make_small_liouv_config(BohrDomain(); construction=KMS())`
   (L_dense built from this is reused throughout the outer @testset, fixing both DIAG-02 fixed-point test at line 63 and the bundle test at line 245)

2. Line 312 (TrotterDomain section setup, fixes issue #4 partially):
   `config_trott = make_small_liouv_config(TrotterDomain())`
   → `config_trott = make_small_liouv_config(TrotterDomain(); construction=KMS())`

3. Line 383 (Trotter/Bohr ratio comparison, fixes issue #4 fully):
   `config_bohr = make_small_liouv_config(BohrDomain())`
   → `config_bohr = make_small_liouv_config(BohrDomain(); construction=KMS())`

4. Line 513 (backward compat testset, fixes issue #5):
   `config_bohr = make_small_liouv_config(BohrDomain())`
   → `config_bohr = make_small_liouv_config(BohrDomain(); construction=KMS())`

Do NOT change any thresholds. Do NOT modify any other lines.
  </action>
  <verify>Run: `cd /Users/bence/code/QuantumFurnace.jl && julia --project=. -e 'using Pkg; Pkg.test(; test_args=["diagnostics"])'` — all 5 previously failing assertions must pass.</verify>
  <done>
- test_diagnostics.jl:63: `fp.trace_distance < 0.01` passes (was ~0.035 with GNS, now ≈0 with KMS)
- test_diagnostics.jl:245: bundle test passes
- test_diagnostics.jl:390: Trotter/Bohr ratio in [0.5, 2.0] passes
- test_diagnostics.jl:519: backward compat `trace_distance < 0.01` passes
  </done>
</task>

</tasks>

<verification>
Run the full test suite to confirm no regressions:
`cd /Users/bence/code/QuantumFurnace.jl && julia --project=. -e 'using Pkg; Pkg.test()'`

Expected: all tests pass except the known-deferred test_regression.jl:40.
</verification>

<success_criteria>
All 5 targeted test failures resolved. No new test failures introduced. test_regression.jl:40 remains deferred (not touched).
</success_criteria>

<output>
After completion, create `.planning/quick/38-fix-failing-diagnostics-and-filename-tes/38-SUMMARY.md`
</output>
