---
phase: 40-save-every
verified: 2026-03-01T11:00:00Z
status: passed
score: 6/6 must-haves verified
re_verification: false
---

# Phase 40: Save Every Verification Report

**Phase Goal:** Users control how often trace distance to the Gibbs state is computed during run_thermalize, reducing observation overhead for long simulations while preserving backward compatibility
**Verified:** 2026-03-01T11:00:00Z
**Status:** PASSED
**Re-verification:** No -- initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | run_thermalize accepts save_every keyword that controls trace distance computation frequency | VERIFIED | `save_every::Int = 1` in function signature at line 152 of src/furnace.jl; observation block gated by `if step % save_every == 0` at line 226 |
| 2 | Default save_every=1 produces identical results to pre-change behavior | VERIFIED | Test "save_every=1 matches default" (test/test_save_every.jl:10-18) verifies bit-identical trace_distances and time_steps using same RNG seed |
| 3 | save_every=10 produces trace_distances and time_steps arrays with length 1 + div(actual_steps, save_every) | VERIFIED | Test "save_every=10 array lengths" (test/test_save_every.jl:23-33) checks `length <= expected_saves` and `length >= 2` with matching lengths. Array built from `recorded_steps` which only appends on `step % save_every == 0` (lines 202, 229) |
| 4 | time_steps and trace_distances always have matching lengths | VERIFIED | Both arrays are extended together inside the same `if` block (lines 228-229: `push!(trace_distances, dist)` then `push!(recorded_steps, step)`). `time_steps = T.(recorded_steps .* config.delta)` at line 237 ensures 1:1 mapping. Test at line 32 explicitly asserts `length(result.time_steps) == length(result.trace_distances)` |
| 5 | Convergence cutoff only triggers on save points (correct coarser detection) | VERIFIED | `if dist < convergence_cutoff; break; end` at lines 231-233 is inside the `if step % save_every == 0` block (line 226). No other convergence check exists outside this gate |
| 6 | save_every value is recorded in metadata for provenance | VERIFIED | `metadata[:save_every] = save_every` at line 241. Test "metadata contains save_every" (test/test_save_every.jl:57-62) verifies via `haskey` and value equality |

**Score:** 6/6 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `src/furnace.jl` | run_thermalize with save_every keyword argument | VERIFIED | 251 lines; contains `save_every` at signature (line 152), validation (line 167), gating (line 226), metadata (line 241) |
| `test/test_save_every.jl` | Behavioral tests for save_every feature (min 40 lines) | VERIFIED | 84 lines; 6 substantive test cases covering backward compatibility, array lengths, stride correctness, metadata, validation, cross-domain |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| src/furnace.jl | trace_distance_h | `step % save_every == 0` gating | WIRED | Line 226: `if step % save_every == 0` gates the `trace_distance_h` call at line 227. Physics blocks (lines 209-224) remain unconditionally outside the gate |
| src/furnace.jl | ThermalizeResults | `metadata[:save_every] = save_every` | WIRED | Line 241 stores save_every in metadata dict; metadata is passed to ThermalizeResults constructor at line 248 |
| test/test_save_every.jl | test/runtests.jl | include statement | WIRED | Line 22 of runtests.jl: `include("test_save_every.jl")`, positioned after test_convergence.jl and before test_observable_trajectories.jl |

### Requirements Coverage

| Requirement | Status | Blocking Issue |
|-------------|--------|----------------|
| SAVE-01: save_every keyword controls trace distance computation frequency | SATISFIED | None |
| SAVE-02: Default save_every=1 preserves backward compatibility | SATISFIED | None |

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| (none) | - | - | - | No anti-patterns detected |

No TODOs, FIXMEs, placeholders, empty implementations, or stub patterns found in modified files.

### Human Verification Required

### 1. Performance Improvement with Large save_every

**Test:** Run `run_thermalize` with mixing_time=10.0 and save_every=100 vs save_every=1 on a 4-qubit system; compare wall_time in metadata
**Expected:** save_every=100 should be noticeably faster due to skipping eigendecomposition in trace_distance_h on non-save steps
**Why human:** Performance characteristics depend on hardware and system size; cannot verify speedup programmatically without running the simulation

### Gaps Summary

No gaps found. All 6 observable truths are verified with concrete evidence in the codebase:

1. The `save_every` keyword is present in the function signature with default value 1 (backward compatible).
2. Input validation (`@assert save_every >= 1`) prevents invalid values.
3. The physics hot loop (coherent unitary, rho_jump accumulation, precomputed channel) runs unconditionally every step.
4. The observation block (trace_distance_h, push, printf, convergence check) is gated by `step % save_every == 0`.
5. `recorded_steps` tracks which steps produced observations, and `time_steps` is derived from it -- ensuring perfect alignment with `trace_distances`.
6. `metadata[:save_every]` provides provenance.
7. The test file contains 6 substantive test cases covering all key behaviors, is included in runtests.jl, and the SUMMARY reports 1158 tests passing.
8. Both commits (f50dddf, 796956a) exist in git history with correct content.

---

_Verified: 2026-03-01T11:00:00Z_
_Verifier: Claude (gsd-verifier)_
