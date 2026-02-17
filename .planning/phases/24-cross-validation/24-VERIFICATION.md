---
phase: 24-cross-validation
verified: 2026-02-17T11:53:43Z
status: gaps_found
score: 2/3 must-haves verified
gaps:
  - truth: "Validation script demonstrates gap estimation agreement for n=4 and n=6 Heisenberg chains, with fitted gap within confidence interval of exact gap (or within documented tolerance)"
    status: failed
    reason: "The trajectory-fitted gap is ~20-28x larger than the exact Liouvillian spectral gap due to a time normalization factor (delta_eff = delta * n_jumps per step). The script pass criterion is `within_ci || relative_error < 0.3`, but relative_error is ~19-27 (1900-2700%) making both conditions false. The SUMMARY documents this as a 'scientific finding' but the script still outputs 'NEEDS MORE TRAJECTORIES', not 'PASS'. The 'documented tolerance' clause in the ROADMAP success criterion is not satisfied by merely naming it a scientific finding -- no corrected comparison or normalized agreement is demonstrated."
    artifacts:
      - path: "experiments/validate_gap_estimation.jl"
        issue: "Pass criterion `within_ci || relative_error < 0.3` will fail for both n=4 and n=6 because the normalization factor is ~20-28x. The script correctly reports the normalization factor but cannot evaluate agreement under normalized conditions."
    missing:
      - "Either: (a) apply the normalization factor to the fitted gap before cross-validation, OR (b) document the specific agreed-upon tolerance in the script and update the pass criterion to reflect it, OR (c) add a normalized_relative_error field to CrossValidationResult and check that instead"
      - "The ROADMAP criterion 'or within documented tolerance' requires a DEFINED tolerance that the script actually checks -- not just printing the factor and failing the comparison"
---

# Phase 24: Cross-Validation Verification Report

**Phase Goal:** Users can verify that trajectory-fitted spectral gap agrees with exact Liouvillian eigenvalues, establishing trust in the method for systems where exact diagonalization is infeasible
**Verified:** 2026-02-17T11:53:43Z
**Status:** gaps_found
**Re-verification:** No -- initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|---------|
| 1 | User can call `cross_validate_gap(estimated, exact_result)` comparing fitted gap against `abs(real(exact_result.spectral_gap))`, with relative error reported | VERIFIED | `cross_validate_gap(estimated::SpectralGapResult, exact_result::LindbladianResult)` exists at line 246 of `src/gap_estimation.jl`. It extracts `exact_result.spectral_gap` and delegates to the Complex method. The Complex method at line 275 computes `exact_gap = abs(real(exact_eigenvalue))` (locked decision enforced) and returns a `CrossValidationResult` with `relative_error`, `absolute_error`, `within_ci`, etc. Both methods are exported via `src/QuantumFurnace.jl` line 54. |
| 2 | Cross-validation warns when `|Im/Re| > 0.1`, indicating oscillatory decay | VERIFIED | Line 292-295 of `src/gap_estimation.jl` emits `@warn "Exact eigenvalue has significant imaginary part (|Im/Re| = ...)"` when `im_warning == true`. Threshold is `im_ratio > 0.1`. Unit test at test line 256 verifies with `@test_logs (:warn, r"significant imaginary part")`. |
| 3 | Validation script demonstrates gap estimation agreement for n=4 and n=6 Heisenberg chains, with fitted gap within CI of exact gap (or within documented tolerance) | FAILED | The script exists at `experiments/validate_gap_estimation.jl` and runs both n=4 and n=6, but the SUMMARY for plan 02 explicitly documents that the fitted trajectory gap is ~20x larger than the exact Liouvillian gap for n=4 and ~28x for n=6. This is due to a time normalization factor (`delta_eff = delta * n_jumps per step`). The script's own pass criterion `within_ci \|\| relative_error < 0.3` will fail: relative_error is ~19-27 (far above 0.3) and within_ci is NO. The script outputs "NEEDS MORE TRAJECTORIES", not "PASS". The ROADMAP says "or within documented tolerance" but the script does not define or check any alternative tolerance that accounts for the normalization factor. |

**Score:** 2/3 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `src/gap_estimation.jl` | CrossValidationResult struct and cross_validate_gap function | VERIFIED | Struct at line 217 with 7 fields (fitted_gap, exact_gap, relative_error, absolute_error, within_ci, imaginary_ratio, imaginary_warning). Two methods at lines 246 and 275. File is 301 lines -- substantive. |
| `test/test_gap_estimation.jl` | Unit tests for cross_validate_gap (Cross-Validation testset) | VERIFIED | "Cross-Validation" testset at line 180 with 8 test cases: field correctness, within_ci true, within_ci false, imaginary_warning false, imaginary_warning true, LindbladianResult dispatch, @warn emission, edge case Re==0. All substantive. |
| `experiments/validate_gap_estimation.jl` | Standalone validation script for n=4 and n=6 | PARTIAL | File exists (291 lines), is substantive, wired to `estimate_spectral_gap`, `cross_validate_gap`, `construct_lindbladian`. However, it will output "NEEDS MORE TRAJECTORIES" for both systems due to the ~20-28x normalization factor gap. The pass criterion in the script does not match the actual behavior. |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `cross_validate_gap` | `SpectralGapResult` | `estimated.gap`, `estimated.gap_ci` | WIRED | Lines 277-289 of `gap_estimation.jl` read `estimated.gap` and `estimated.gap_ci[1]`, `estimated.gap_ci[2]`. |
| `cross_validate_gap` | `LindbladianResult` | `exact_result.spectral_gap` | WIRED | Line 247 reads `exact_result.spectral_gap` and passes it to the Complex method. |
| `experiments/validate_gap_estimation.jl` | `estimate_spectral_gap` | function call | WIRED | Line 201 calls `estimate_spectral_gap(jumps, therm_config, psi0, hamiltonian; ...)`. |
| `experiments/validate_gap_estimation.jl` | `cross_validate_gap` | function call | WIRED | Line 219 calls `cv = cross_validate_gap(estimated, exact_result)`. |
| `experiments/validate_gap_estimation.jl` | `construct_lindbladian` | function call | WIRED | Line 181 calls `L = construct_lindbladian(jumps, liouv_config, hamiltonian)`. |

### Requirements Coverage

| Requirement | Status | Blocking Issue |
|-------------|--------|---------------|
| VAL-01: User can cross-validate fitted gap against exact Liouvillian spectral gap using `abs(real(spectral_gap))` | SATISFIED | `cross_validate_gap` correctly implements this. |
| VAL-02: Cross-validation warns when imaginary part of exact eigenvalue is significant (|Im/Re| > 0.1) | SATISFIED | `@warn` is emitted, tested with `@test_logs`. |
| VAL-03: Validation script demonstrates gap estimation agreement for n=4 and n=6 | BLOCKED | Script runs but reports ~19-27x relative error and outputs "NEEDS MORE TRAJECTORIES". Agreement is not demonstrated -- only the discrepancy is documented. |

### Anti-Patterns Found

| File | Pattern | Severity | Impact |
|------|---------|----------|--------|
| `experiments/validate_gap_estimation.jl` (line 239) | Pass criterion `within_ci \|\| relative_error < 0.3` will evaluate to false for both n=4 and n=6 given the known ~20x normalization factor. Script outputs "NEEDS MORE TRAJECTORIES". | Blocker | The validation script cannot demonstrate agreement under its own pass criterion. The ROADMAP success criterion 3 requires demonstration of agreement (within CI or within documented tolerance), not just observation of disagreement with a normalization note. |

### Human Verification Required

None -- all checks performed programmatically. The normalization factor issue is documented in the SUMMARY and is conclusive.

### Gaps Summary

The cross-validation API (plan 01) is fully implemented and correct:
- `CrossValidationResult` struct with 7 fields
- `cross_validate_gap` with two-method dispatch (LindbladianResult and Complex)
- `@warn` emission for imaginary ratio > 0.1
- 8 unit tests all passing

The validation script (plan 02) has a fundamental unresolved issue: the trajectory-based gap estimation produces rates ~20-28x larger than the continuous-time Liouvillian spectral gap, due to a time normalization factor (`delta_eff = delta * n_jumps per step`). The SUMMARY correctly documents this as a physics finding but does not resolve it before declaring the phase complete.

The ROADMAP success criterion 3 says "fitted gap within confidence interval of exact gap **(or within documented tolerance)**". The SUMMARY interprets "documented tolerance" as simply printing the normalization factor. But the script's own pass logic (`within_ci || relative_error < 0.3`) will print "NEEDS MORE TRAJECTORIES" for both systems -- the script itself declares failure.

For this criterion to be satisfied, one of the following must be true:
1. The script applies the normalization factor and shows normalized agreement within CI or < 30% error
2. The script defines a specific documented tolerance (e.g., "within 30x normalization factor") and checks it
3. The ROADMAP success criterion is updated to match the actual physics (the trajectory rate and Liouvillian gap differ by a known factor, and this is the intended result to demonstrate)

---

_Verified: 2026-02-17T11:53:43Z_
_Verifier: Claude (gsd-verifier)_
