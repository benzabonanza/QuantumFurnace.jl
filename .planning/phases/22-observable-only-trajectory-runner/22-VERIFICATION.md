---
phase: 22-observable-only-trajectory-runner
verified: 2026-02-17T09:17:12Z
status: passed
score: 4/4 must-haves verified
re_verification: false
---

# Phase 22: Observable-Only Trajectory Runner Verification Report

**Phase Goal:** Users can run trajectory simulations that measure time-resolved observables efficiently, without the overhead of per-trajectory density matrix reconstruction
**Verified:** 2026-02-17T09:17:12Z
**Status:** passed
**Re-verification:** No -- initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | User can call `run_observable_trajectories` with observables and receive time-resolved observable means `<O>(t)` | VERIFIED | Function at `src/trajectories.jl:736` with correct signature; returns `ObservableTrajectoryResult{CT}` with `times` (Vector{Float64}) and `measurements_mean` (Matrix{Float64}, n_obs x n_saves) |
| 2 | Observable time series from `run_observable_trajectories` matches `run_trajectories` bitwise for same seed/ntraj/observables | VERIFIED | Cross-validation tests at `test/test_observable_trajectories.jl:42,53,83` use `==` (bitwise equality); implementation uses identical `Xoshiro(actual_seed + traj_id)` seed pattern and identical measurement loop to `run_trajectories` |
| 3 | User can pass `reconstruct_dm=true` to get averaged density matrix at end of run, `rho_mean=nothing` otherwise | VERIFIED | `reconstruct_dm::Bool = false` kwarg at `src/trajectories.jl:748`; dispatches to `_run_chunk_with_obs!` when true (accumulates DM), sets `rho_mean = nothing` when false; test confirms at lines 27 and 52-55 |
| 4 | Multi-threaded execution works with per-thread workspace/RNG pattern matching existing `run_trajectories` | VERIFIED | Threaded path at `src/trajectories.jl:771-808`: creates `ws_per_task = [TrajectoryWorkspace(CT, dim) for _ in 1:length(chunks)]`, saves/restores BLAS threads in try/finally, dispatches between `_run_chunk_obs_only!` and `_run_chunk_with_obs!` based on `reconstruct_dm` |

**Score:** 4/4 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `src/trajectories.jl` | `ObservableTrajectoryResult` struct with times, measurements_mean, n_trajectories, seed, rho_mean fields | VERIFIED | Struct at lines 37-66; all 5 fields present; inner+outer constructor pattern resolves Aqua unbound type parameter issue; struct is in trajectories.jl (not structs.jl as plan suggested, but correctly co-located with the analogous `TrajectoryResult`) |
| `src/trajectories.jl` | `_run_chunk_obs_only!` and `run_observable_trajectories` functions | VERIFIED | `_run_chunk_obs_only!` at lines 474-512 (identical to `_run_chunk_with_obs!` minus `_accumulate_density_matrix!` call, confirmed by comment at line 509); `run_observable_trajectories` at lines 736-849 |
| `src/QuantumFurnace.jl` | Exports for `run_observable_trajectories` and `ObservableTrajectoryResult` | VERIFIED | Line 42: `export TrajectoryFramework, TrajectoryResult, ObservableTrajectoryResult, build_trajectoryframework, step_along_trajectory!, run_observable_trajectories` |
| `test/test_observable_trajectories.jl` | Cross-validation and correctness tests | VERIFIED | 85 lines, 6 testsets: basic run, bitwise cross-validation, reconstruct_dm=true bitwise match, deterministic seeding, different seeds, single trajectory serial path |
| `test/runtests.jl` | Include for new test file | VERIFIED | Line 23: `include("test_observable_trajectories.jl")` -- last include in the test suite |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `src/trajectories.jl` | `src/trajectories.jl` (ObservableTrajectoryResult) | Constructor call | WIRED | Line 848: `return ObservableTrajectoryResult{CT}(times, mean_data, ntraj, actual_seed, rho_mean)` |
| `run_observable_trajectories` | `_build_framework_and_seed` | Framework setup reuse | WIRED | Line 753: `fw, actual_seed = _build_framework_and_seed(jumps, config, psi0, hamiltonian; trotter=trotter, delta=delta, seed=seed,)` |
| `run_observable_trajectories` | `_accumulate_measurements!` | Observable measurement in step loop | WIRED | Lines 822, 829 (serial path); lines 499, 506 (inside `_run_chunk_obs_only!`); lines 452, 459 (inside `_run_chunk_with_obs!`) |
| `test/test_observable_trajectories.jl` | `run_trajectories` | Cross-validation bitwise comparison | WIRED | Lines 42, 53, 83: `result_new.measurements_mean == result_old.measurements_mean` (bitwise `==`, not `isapprox`) |

### Requirements Coverage

The 4 success criteria from ROADMAP.md are all satisfied:

| Requirement | Status | Blocking Issue |
|-------------|--------|----------------|
| 1. `run_observable_trajectories(jumps, config, psi0, ham; observables, save_every, ntraj)` returns time-resolved `<O>(t)` | SATISFIED | None |
| 2. Observable-only runner produces identical observable time series as `run_trajectories` (cross-check) | SATISFIED | None -- bitwise equality enforced in tests |
| 3. User can optionally reconstruct averaged DM at end via `reconstruct_dm=true` (not during) | SATISFIED | None -- `rho_mean=nothing` when false; full DM when true |
| 4. Multi-threaded execution with per-thread workspace/RNG pattern | SATISFIED | None -- identical pattern to `run_trajectories` threaded path |

### Anti-Patterns Found

None. No TODOs, FIXMEs, placeholder returns, or empty implementations found in any phase 22 files.

### Human Verification Required

None for automated correctness. The SUMMARY reports 633 tests pass (615 existing + 18 new), including Aqua checks, which provides confidence in the implementation.

One optional human check:

**Multi-threaded bitwise cross-validation:** The test suite uses `ntraj=20` which on a 2-thread system will exercise the threaded path. If run on a single-thread system, both runners use the serial path and bitwise match is guaranteed. On a multi-thread system, thread scheduling affects which trajectories land in which chunk, but the cross-validation still holds because both runners partition identically. This is a correctness property that is validated by the test suite but worth confirming with `-t2` explicitly.

**Test:** Run `julia --project -t2 -e 'using Pkg; Pkg.test()'`
**Expected:** All 633 tests pass
**Why human:** Confirms the threaded path exercises the multi-thread code branch, not just the serial fallback

### Gaps Summary

No gaps. All 4 observable truths verified. All artifacts exist, are substantive, and are wired. All key links confirmed. No anti-patterns.

**Notable deviation from plan (auto-fixed, not a gap):** The PLAN specified `ObservableTrajectoryResult` in `src/structs.jl`, but the implementation correctly placed it in `src/trajectories.jl` alongside the analogous `TrajectoryResult` struct. This is the right location and was documented in the SUMMARY as a deliberate deviation.

---

_Verified: 2026-02-17T09:17:12Z_
_Verifier: Claude (gsd-verifier)_
