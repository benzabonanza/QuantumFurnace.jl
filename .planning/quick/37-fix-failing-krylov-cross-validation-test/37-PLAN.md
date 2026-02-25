---
phase: 37-fix-failing-krylov-cross-validation-test
plan: 01
type: execute
wave: 1
depends_on: []
files_modified:
  - src/krylov_eigsolve.jl
  - test/test_krylov_crossvalidation.jl
autonomous: true
must_haves:
  truths:
    - "All Krylov cross-validation tests pass (julia --project -e 'using Pkg; Pkg.test()' or targeted test run)"
    - "Docstrings accurately describe the O(delta) convergence behavior of the (mu-1)/delta conversion"
    - "Test thresholds match the mathematically expected convergence order for faithful Chen channel"
  artifacts:
    - path: "test/test_krylov_crossvalidation.jl"
      provides: "Corrected XVAL-03 convergence thresholds and comments"
      contains: "order >= 0.9"
    - path: "src/krylov_eigsolve.jl"
      provides: "Corrected docstring for channel-to-Lindbladian conversion"
  key_links:
    - from: "test/test_krylov_crossvalidation.jl"
      to: "src/krylov_eigsolve.jl"
      via: "krylov_spectral_gap called by run_le_convergence"
      pattern: "krylov_spectral_gap"
---

<objective>
Fix 8 failing XVAL-03 tests in test/test_krylov_crossvalidation.jl by correcting convergence
order thresholds and updating misleading docstrings.

Purpose: Phase 32 replaced the Euler channel (E = I + delta*L, where (mu-1)/delta was exact) with
the faithful Chen CPTP channel (mu = exp(delta*lambda_L) + O(delta^2)). The conversion formula
(mu-1)/delta now introduces O(delta) error from the Taylor expansion of exp, so convergence order
is ~1.0, not ~2.0. The test thresholds and docstrings need to reflect this correct mathematical behavior.

Output: All tests green, accurate documentation.
</objective>

<execution_context>
@./.claude/get-shit-done/workflows/execute-plan.md
@./.claude/get-shit-done/templates/summary.md
</execution_context>

<context>
@test/test_krylov_crossvalidation.jl
@src/krylov_eigsolve.jl
</context>

<tasks>

<task type="auto">
  <name>Task 1: Fix docstring in krylov_eigsolve.jl</name>
  <files>src/krylov_eigsolve.jl</files>
  <action>
In src/krylov_eigsolve.jl, fix the misleading docstring at lines 466-468 for
`krylov_spectral_gap(config::AbstractThermalizeConfig, ...)`.

Replace lines 466-468:
```
The channel eigenvalues mu are related to Lindbladian eigenvalues by the exact linear
formula: `lambda_L = (mu - 1) / delta`. The steady state has mu ~ 1 (largest magnitude),
and the gap is recovered from the second eigenvalue after conversion.
```

With:
```
The channel eigenvalues mu are related to Lindbladian eigenvalues by the first-order
approximation: `lambda_L = (mu - 1) / delta`. Since mu = exp(delta * lambda_L) + O(delta^2),
the conversion introduces O(delta) error. The steady state has mu ~ 1 (largest magnitude),
and the gap is recovered from the second eigenvalue after conversion.
```

This corrects "exact linear formula" to "first-order approximation" and explains the
error source. The formula on line 532 (`(vals .- 1) ./ delta`) stays unchanged since it
is the correct first-order conversion to use.
  </action>
  <verify>Read the modified docstring and confirm it no longer claims "exact".</verify>
  <done>Docstring accurately describes (mu-1)/delta as a first-order approximation with O(delta) error.</done>
</task>

<task type="auto">
  <name>Task 2: Fix XVAL-03 test thresholds and comments</name>
  <files>test/test_krylov_crossvalidation.jl</files>
  <action>
In test/test_krylov_crossvalidation.jl, make these changes:

1. Update the `run_le_convergence` docstring (lines 106-107). Replace:
```
The channel eigenvalue mapping error is O(delta^2), so the convergence order
should be >= 1.5 (with margin for sub-leading terms).
```
With:
```
The faithful Chen channel has mu = exp(delta*lambda_L) + O(delta^2), so the
first-order conversion (mu-1)/delta introduces O(delta) error. Expected
convergence order is ~1.0; threshold is >= 0.9 with margin for sub-leading terms.
```

2. Update the XVAL-03 section comments (lines 388-391). Replace:
```
# XVAL-03: L-vs-E convergence (KMS only, per locked decision)
# Tests that channel-to-Lindbladian gap mapping converges with O(delta^2)
# Deltas: [0.1, 0.01, 0.001]
# Hard assertion: convergence order >= 1.5 for each consecutive pair
```
With:
```
# XVAL-03: L-vs-E convergence (KMS only, per locked decision)
# Tests that channel-to-Lindbladian gap mapping converges with O(delta)
# The faithful Chen channel gives mu = exp(delta*lambda_L) + O(delta^2),
# so (mu-1)/delta has first-order error. Deltas: [0.1, 0.01, 0.001]
# Hard assertion: convergence order >= 0.9 for each consecutive pair
```

3. Change ALL FOUR convergence order thresholds from `1.5` to `0.9`:
   - Line 398: `@test order >= 0.9`  (EnergyDomain)
   - Line 405: `@test order >= 0.9`  (TimeDomain)
   - Line 413: `@test order >= 0.9`  (TrotterDomain)
   - Line 420: `@test order >= 0.9`  (BohrDomain)
  </action>
  <verify>Run the XVAL-03 tests specifically:
```bash
cd /Users/bence/code/QuantumFurnace.jl && julia --project -e '
  using Pkg; Pkg.instantiate()
  using Test
  include("test/test_krylov_crossvalidation.jl")
'
```
All 8 XVAL-03 tests (2 orders x 4 domains) should pass.
  </verify>
  <done>All 8 previously-failing XVAL-03 tests pass with the corrected threshold of >= 0.9. No other tests are broken.</done>
</task>

</tasks>

<verification>
Run the full Krylov cross-validation test file to confirm no regressions:
```bash
julia --project -e 'using Pkg; Pkg.instantiate(); using Test; include("test/test_krylov_crossvalidation.jl")'
```
All test sets (XVAL-01, XVAL-02, XVAL-03) should pass.
</verification>

<success_criteria>
- All 8 XVAL-03 tests pass (convergence order >= 0.9 threshold)
- No regressions in XVAL-01 or XVAL-02
- Docstrings in both files accurately describe O(delta) convergence behavior
- No code logic changes (only comments, docstrings, and test thresholds)
</success_criteria>

<output>
After completion, create `.planning/quick/37-fix-failing-krylov-cross-validation-test/37-01-SUMMARY.md`
</output>
