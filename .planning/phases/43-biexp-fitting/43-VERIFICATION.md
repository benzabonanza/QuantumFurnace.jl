---
phase: 43-biexp-fitting
verified: 2026-03-04T09:44:14Z
status: passed
score: 5/5 must-haves verified
---

# Phase 43: Bi-Exponential Fitting Verification Report

**Phase Goal:** Reduce mixing time extrapolation error from 26% to <5% by adding bi-exponential fitting with explicit :biexp model keyword, preserving backward compatibility with default :single model.
**Verified:** 2026-03-04T09:44:14Z
**Status:** passed
**Re-verification:** No -- initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | BiexpFitResult struct exists with all required fields (gap, gap_fast, amplitude, amplitude_fast, offset, gap_ci, gap_se, r_squared, converged, residuals, times_used, values_used) | VERIFIED | `src/fitting.jl:269-282` -- struct with all 12 fields defined |
| 2 | fit_biexponential_decay function fits 5-parameter bi-exp model with mode sorting (g1 >= g2), residual-seeded initial guess, and >=8 data point requirement | VERIFIED | `src/fitting.jl:351-437` -- full implementation with bounds, mode sorting at L401-415, SE/CI pre-swap index tracking at L407/414, auto initial guess via `_biexp_initial_guess` at L298-326 |
| 3 | estimate_mixing_time accepts model=:biexp keyword, uses Roots.Bisection for numerical extrapolation, and stores BiexpFitResult in result | VERIFIED | `src/mixing.jl:288` default `:single`, L301 validation, L331-356 biexp path calls `fit_biexponential_decay`, `_extrapolate_mixing_time_biexp` (L157-203 uses `Roots.find_zero` with `Roots.Bisection()` at L199), synthetic FitResult at L348 |
| 4 | Bi-exp extrapolation achieves <5% error on synthetic bi-exp data (acceptance criterion) | VERIFIED | `test/test_mixing.jl:257` -- `@test biexp_err < 0.05` with synthetic data matching the problem scenario (C_true=6.8e-5, target=1e-4). SUMMARY reports actual error <0.001% |
| 5 | Backward compatibility preserved: default model=:single produces identical results, biexp_fit_result===nothing for single mode | VERIFIED | `src/mixing.jl:288` default `:single`, L328 returns `nothing` for biexp_fit_result. Test BIEXP-MIX-02 at `test/test_mixing.jl:276-299` verifies `est.model_used == :single`, `est.biexp_fit_result === nothing`, and identical gap/mixing_time/offset vs explicit `model=:single` |

**Score:** 5/5 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `src/fitting.jl` | BiexpFitResult struct + fit_biexponential_decay function + _biexp_initial_guess + _biexp_decay_model | VERIFIED | 437 lines total. Bi-exp code appended at L219-437 (~220 new lines). Struct at L269-282, model at L239-241, initial guess at L298-326, main function at L351-437 |
| `src/mixing.jl` | MixingTimeEstimate extended with model_used + biexp_fit_result fields; estimate_mixing_time with model keyword; _extrapolate_mixing_time_biexp; _biexp_to_single_fit_result | VERIFIED | 357 lines. Struct fields at L74-75, biexp extrapolation at L157-203, synthetic FitResult converter at L212-225, biexp branch in main API at L331-356 |
| `src/QuantumFurnace.jl` | export fit_biexponential_decay, BiexpFitResult | VERIFIED | Line 79: `export fit_biexponential_decay, BiexpFitResult` |
| `test/test_fitting.jl` | BIEXP-01, BIEXP-02, BIEXP-03, edge case tests | VERIFIED | Lines 162-256. BIEXP-01 (clean data recovery, L162-195), BIEXP-02 (offset accuracy comparison, L200-222), BIEXP-03 (skip_initial, L227-246), edge case too-few-points (L251-255) |
| `test/test_mixing.jl` | BIEXP-MIX-01 (<5% accuracy), BIEXP-MIX-02 (backward compat), edge cases, MIX-05 updated | VERIFIED | Lines 222-310. BIEXP-MIX-01 (L223-271, `@test biexp_err < 0.05`), BIEXP-MIX-02 (L276-299), invalid model edge case (L304-310), MIX-05 updated with model_used/biexp_fit_result fields (L122-126) |
| `scripts/mixing_time_extrapolate_verify.jl` | Side-by-side single vs biexp comparison | VERIFIED | Lines 157-206. Calls `estimate_mixing_time` with `model=:biexp` at L158-163, prints both predictions with error percentages at L199-206 |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `estimate_mixing_time` (mixing.jl:333) | `fit_biexponential_decay` (fitting.jl:351) | Direct function call in biexp branch | WIRED | L333: `bifit = fit_biexponential_decay(...)` |
| `estimate_mixing_time` (mixing.jl:337) | `_extrapolate_mixing_time_biexp` (mixing.jl:157) | Direct call with bifit result | WIRED | L337: `t_mix_extrap = extrapolate ? _extrapolate_mixing_time_biexp(bifit, target_epsilon) : nothing` |
| `_extrapolate_mixing_time_biexp` (mixing.jl:199) | `Roots.find_zero` / `Roots.Bisection` | Numerical root finding for multi-exp equation | WIRED | L199: `Roots.find_zero(biexp_residual, (0.0, t_upper), Roots.Bisection())` |
| `fit_biexponential_decay` (fitting.jl:379-380) | `_biexp_initial_guess` (fitting.jl:298) | Residual-seeded auto initial guess | WIRED | L379-380: `single_fit = fit_exponential_decay(...)` then `_biexp_initial_guess(times_fit, values_fit, single_fit)` |
| `estimate_mixing_time` (mixing.jl:348) | `_biexp_to_single_fit_result` (mixing.jl:212) | Backward-compat synthetic FitResult | WIRED | L348: `synthetic_fit = _biexp_to_single_fit_result(bifit)` |
| `QuantumFurnace.jl` (L79) | `BiexpFitResult` + `fit_biexponential_decay` | Module export | WIRED | L79: `export fit_biexponential_decay, BiexpFitResult` |
| `scripts/mixing_time_extrapolate_verify.jl` (L158-163) | `estimate_mixing_time` with `model=:biexp` | Script usage | WIRED | L159: `model = :biexp` in keyword arguments |

### Requirements Coverage

| Requirement | Status | Evidence |
|-------------|--------|----------|
| Bi-exp extrapolation <5% error on synthetic data | SATISFIED | Test BIEXP-MIX-01: `@test biexp_err < 0.05` |
| Backward compatibility with default model=:single | SATISFIED | Test BIEXP-MIX-02: default returns `:single`, `biexp_fit_result === nothing`, identical results |
| All existing tests pass unchanged | SATISFIED | SUMMARY reports 1273 tests pass (1246 existing + 27 new); MIX-05 updated to include new fields |

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| (none) | - | - | - | - |

No TODOs, FIXMEs, placeholders, empty implementations, or stub patterns found in any modified source files.

### Commits Verified

All 6 task commits from SUMMARY confirmed present in git log:

| Commit | Description | Verified |
|--------|-------------|----------|
| `094184f` | feat(43): add bi-exponential decay fitting to src/fitting.jl | Present |
| `13324e2` | feat(43): extend estimate_mixing_time with model=:biexp support | Present |
| `94e2c50` | feat(43): export fit_biexponential_decay and BiexpFitResult | Present |
| `c82c63b` | test(43): add bi-exponential fitting tests BIEXP-01/02/03 | Present |
| `3b2d690` | test(43): add bi-exp mixing time tests BIEXP-MIX-01/02 | Present |
| `6ccccf1` | feat(43): add bi-exp comparison to mixing time verification script | Present |

### Human Verification Required

### 1. Run full test suite

**Test:** `julia --project -e 'using Pkg; Pkg.test()'`
**Expected:** All 1273 tests pass (1246 existing + 27 new), no regressions
**Why human:** Cannot run Julia test suite in this verification environment

### 2. Run real-data verification script

**Test:** `OPENBLAS_NUM_THREADS=4 julia -t4 --project scripts/mixing_time_extrapolate_verify.jl`
**Expected:** Bi-exp prediction shows lower error % than single-exp prediction on actual Heisenberg chain data. Target epsilon should be reached in the verification run.
**Why human:** Requires full Julia environment with compiled QuantumFurnace module and Hamiltonian data files. This is the ultimate acceptance test on real (non-synthetic) data.

### 3. Confirm bi-exp improves real extrapolation error from 26% to <5%

**Test:** Compare single-exp vs bi-exp error percentages in the script output "VERIFICATION RESULTS" section
**Expected:** Bi-exp error significantly lower than single-exp error, ideally <5%
**Why human:** The 26% -> <5% claim comes from synthetic data tests; real-data validation requires running the full simulation

---

_Verified: 2026-03-04T09:44:14Z_
_Verifier: Claude (gsd-verifier)_
