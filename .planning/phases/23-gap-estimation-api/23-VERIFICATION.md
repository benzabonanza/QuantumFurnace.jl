---
phase: 23-gap-estimation-api
verified: 2026-02-17T09:46:48Z
status: passed
score: 3/3 must-haves verified
---

# Phase 23: Gap Estimation API Verification Report

**Phase Goal:** Users can estimate the Lindbladian spectral gap from a single function call that orchestrates trajectory simulation, multi-observable fitting, and best-estimate selection
**Verified:** 2026-02-17T09:46:48Z
**Status:** passed
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | User can call `estimate_spectral_gap(jumps, config, psi0, ham)` and receive a `SpectralGapResult` | VERIFIED | `src/gap_estimation.jl` line 127: full function definition; `src/QuantumFurnace.jl` line 54: exported; test/test_gap_estimation.jl line 21: called and asserted `result isa SpectralGapResult` |
| 2 | `SpectralGapResult` contains gap estimate, CI, per-observable FitResult vector, best observable name, and fit metadata (ntraj, total_time, save_every, seed, skip_initial) | VERIFIED | `src/gap_estimation.jl` lines 37-50: struct with exactly 12 fields (gap, gap_ci, gap_se, best_observable, best_r_squared, per_observable, observable_names, ntraj, total_time, save_every, seed, skip_initial); all fields populated in constructor at lines 177-190 |
| 3 | Best observable is automatically selected by fit quality (converged AND gap > 0 AND highest R-squared, with fallback to highest R-squared if no valid fit), and the selection is explainable from per_observable results | VERIFIED | `src/gap_estimation.jl` lines 67-85: `_select_best_observable` iterates fits, selects converged+gap>0+highest-r_squared, fallback to argmax(r_squared); test line 151-166: selection logic tested with controlled FitResult instances including fallback path |

**Score:** 3/3 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `src/gap_estimation.jl` | SpectralGapResult struct, estimate_spectral_gap function, _select_best_observable helper | VERIFIED | 192 lines; contains `struct SpectralGapResult` (line 37), `function _select_best_observable` (line 67), `function estimate_spectral_gap` (line 127); no stubs or TODOs |
| `test/test_gap_estimation.jl` | Integration tests for estimate_spectral_gap with SMALL 3-qubit system, min 50 lines | VERIFIED | 169 lines (well above 50-line minimum); 7 test sets covering GAP-01, GAP-02, GAP-03 |
| `src/QuantumFurnace.jl` | include and export of gap_estimation.jl | VERIFIED | Line 112: `include("gap_estimation.jl")`; Line 54: `export SpectralGapResult, estimate_spectral_gap`; include order correct (after fitting.jl line 111, before results.jl line 113) |
| `test/runtests.jl` | include("test_gap_estimation.jl") | VERIFIED | Line 24: `include("test_gap_estimation.jl")` present after test_observable_trajectories.jl |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `src/gap_estimation.jl` | `src/convergence.jl` | `build_gap_estimation_observables` call | WIRED | gap_estimation.jl line 144 calls `build_gap_estimation_observables`; function defined at convergence.jl line 127 |
| `src/gap_estimation.jl` | `src/trajectories.jl` | `run_observable_trajectories` call | WIRED | gap_estimation.jl line 156 calls `run_observable_trajectories`; function defined at trajectories.jl line 736 |
| `src/gap_estimation.jl` | `src/fitting.jl` | `fit_exponential_decay` call per observable | WIRED | gap_estimation.jl line 167 calls `fit_exponential_decay` inside loop over observables; function defined at fitting.jl line 153 |
| `src/QuantumFurnace.jl` | `src/gap_estimation.jl` | include and export | WIRED | QuantumFurnace.jl line 112: `include("gap_estimation.jl")`; line 54: `export SpectralGapResult, estimate_spectral_gap` |

### Requirements Coverage

| Requirement | Status | Notes |
|-------------|--------|-------|
| GAP-01: Single-call API returning SpectralGapResult | SATISFIED | `estimate_spectral_gap` callable with documented signature, returns `SpectralGapResult` |
| GAP-02: SpectralGapResult struct with all required fields | SATISFIED | All 12 fields present: gap, gap_ci, gap_se, best_observable, best_r_squared, per_observable, observable_names, ntraj, total_time, save_every, seed, skip_initial |
| GAP-03: Best observable selected automatically, explainable from per_observable | SATISFIED | `_select_best_observable` implements converged+gap>0+highest-R2 with fallback; `per_observable` vector enables full inspection |

### Anti-Patterns Found

No anti-patterns detected.

- Zero TODO/FIXME/HACK/PLACEHOLDER comments in `src/gap_estimation.jl`
- No stub implementations (return null, empty handlers, etc.)
- All functions have real implementations with data flow from inputs to return value
- `traj_result.measurements_mean[i, :]` correctly indexes row per observable, column per time step
- `Float64.(...)` cast applied before passing to fitter (correct type handling)
- `traj_result.seed` used for result (actual seed, not requested seed — correct for reproducibility)

### Human Verification Required

None. All goal-critical behaviors are verifiable programmatically:
- Struct existence and field presence: verified via source inspection
- Function signature and implementation: verified via source inspection
- Key link wiring: verified via grep of call sites against function definitions
- Test coverage: verified via test file inspection
- Commit existence: both `cf0c438` and `bebc1c1` confirmed in git log

The SUMMARY reports 666 total tests passing. This is consistent with 633 prior tests + approximately 33 new assertions across 7 test sets. The test runner includes `test_gap_estimation.jl` at line 24 of runtests.jl.

### Gaps Summary

No gaps. All three must-have truths are verified with full artifact and wiring evidence.

---

_Verified: 2026-02-17T09:46:48Z_
_Verifier: Claude (gsd-verifier)_
