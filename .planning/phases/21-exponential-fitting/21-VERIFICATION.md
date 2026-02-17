---
phase: 21-exponential-fitting
verified: 2026-02-17T08:46:57Z
status: passed
score: 5/5 must-haves verified
re_verification: false
---

# Phase 21: Exponential Fitting Verification Report

**Phase Goal:** Users can fit exponential decay curves to time-series data and extract decay rates with confidence intervals, validated on synthetic data before touching real trajectories
**Verified:** 2026-02-17T08:46:57Z
**Status:** passed
**Re-verification:** No -- initial verification

## Goal Achievement

### Observable Truths

| #  | Truth                                                                                                         | Status     | Evidence                                                                                                              |
|----|---------------------------------------------------------------------------------------------------------------|------------|-----------------------------------------------------------------------------------------------------------------------|
| 1  | User can call `fit_exponential_decay(times, values)` and recover a known decay rate within the returned CI    | VERIFIED   | `fitting.jl:153-211` full LM fit with bounds; testset "FIT-01: noisy data recovery within CI" asserts `gap_ci[1] <= gap_true <= gap_ci[2]` |
| 2  | Fitting auto-generates initial guess via log-linear estimate -- user need not provide initial parameters      | VERIFIED   | `_log_linear_initial_guess` at `fitting.jl:73-102`; `fit_exponential_decay` calls it when `p0 === nothing`; testset "FIT-02: auto-generated initial guess" asserts convergence + atol=0.01 |
| 3  | User can specify `skip_initial` fraction and the fitted gap changes measurably                                | VERIFIED   | `fitting.jl:169-174` slices data by `start_idx`; testset "FIT-03: skip_initial window selection" asserts `abs(r2.gap - gap_true) < abs(r1.gap - gap_true)` |
| 4  | Fit result includes R-squared, CI on gap, standard error; gap constrained positive via bounds                 | VERIFIED   | `FitResult` fields at `fitting.jl:44-55`; bounds `lower=[-Inf, 0.0, -Inf]` at `fitting.jl:184`; testsets "FIT-04" and "FIT-05" verify all fields and `result.gap >= 0.0` |
| 5  | LsqFit.jl is added to Project.toml with proper compat bounds                                                  | VERIFIED   | `Project.toml:16` has `LsqFit = "2fda8390..."` in `[deps]`; `Project.toml:46` has `LsqFit = "0.15"` in `[compat]` |

**Score:** 5/5 truths verified

### Required Artifacts

| Artifact              | Expected                                      | Status   | Details                                                                                                    |
|-----------------------|-----------------------------------------------|----------|------------------------------------------------------------------------------------------------------------|
| `src/fitting.jl`      | FitResult struct and fit_exponential_decay    | VERIFIED | 212 lines; complete struct (10 fields), full function implementation, two internal helpers, full docstrings |
| `test/test_fitting.jl`| Synthetic data tests for all FIT-* requirements | VERIFIED | 155 lines; 9 testsets covering FIT-01 through FIT-05 plus custom p0 and log-linear fallback               |
| `Project.toml`        | LsqFit dependency with compat bound           | VERIFIED | LsqFit in `[deps]` (line 16) and `[compat]` with `"0.15"` (line 46)                                       |

### Key Link Verification

| From                   | To                          | Via                                 | Status   | Details                                                                                           |
|------------------------|-----------------------------|-------------------------------------|----------|---------------------------------------------------------------------------------------------------|
| `src/fitting.jl`       | `LsqFit.curve_fit`          | `using LsqFit` in module             | WIRED    | `QuantumFurnace.jl:22` has `using LsqFit`; `fitting.jl:189` calls `curve_fit(...; lower=lower, upper=upper)` |
| `src/QuantumFurnace.jl`| `src/fitting.jl`            | `include("fitting.jl")`             | WIRED    | `QuantumFurnace.jl:108` has `include("fitting.jl")`; exports at line 51                           |
| `test/test_fitting.jl` | `fit_exponential_decay`     | synthetic data tests                | WIRED    | `fit_exponential_decay(...)` called in 7 of 9 testsets; `QuantumFurnace._log_linear_initial_guess` called in 2 |

### Requirements Coverage

| Requirement | Status    | Blocking Issue |
|-------------|-----------|----------------|
| FIT-01: recover known gap within CI from synthetic data           | SATISFIED | None |
| FIT-02: auto log-linear initial guess, no p0 required             | SATISFIED | None |
| FIT-03: skip_initial measurably changes fitted gap                | SATISFIED | None |
| FIT-04: R-squared, CI, SE fields; R-squared unclamped             | SATISFIED | None |
| FIT-05: gap >= 0 enforced via parameter bounds                    | SATISFIED | None |
| LsqFit.jl in [deps] and [compat] with "0.15" bound               | SATISFIED | None |

### Anti-Patterns Found

None. No TODO/FIXME/PLACEHOLDER comments, no stub return values, no empty handlers found in `src/fitting.jl` or `test/test_fitting.jl`.

### Human Verification Required

**1. Test suite execution**

**Test:** Run `julia --project=. -e 'using Pkg; Pkg.test()'` from the repo root
**Expected:** All 610+ tests pass (no regressions), all 9 fitting testsets pass including "FIT-01: noisy data recovery within CI" and "FIT-03: skip_initial window selection"
**Why human:** Julia compilation and LsqFit.jl numerical results cannot be verified via static analysis alone. The noisy-data CI coverage test in particular depends on runtime LM convergence behavior.

### Gaps Summary

No gaps found. All 5 must-have truths are verified at all three levels (exists, substantive, wired). All artifacts are present with correct content. All key links are connected. No anti-patterns detected. The two commits (`0eafd22`, `e50f29c`) referenced in SUMMARY.md are present in git history.

---

_Verified: 2026-02-17T08:46:57Z_
_Verifier: Claude (gsd-verifier)_
