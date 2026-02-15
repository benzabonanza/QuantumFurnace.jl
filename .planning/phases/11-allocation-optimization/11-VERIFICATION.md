---
phase: 11-allocation-optimization
verified: 2026-02-15T16:00:00Z
status: passed
score: 5/5 truths verified
re_verification: false
---

# Phase 11: Allocation Optimization Verification Report

**Phase Goal:** Core simulation hot paths avoid unnecessary heap allocations -- sparse matrices, Diagonal wrappers, filter intermediates, and redundant basis transforms are eliminated or precomputed

**Verified:** 2026-02-15T16:00:00Z
**Status:** passed
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | B_bohr inner loop does not allocate sparse matrices per iteration (A_nu matrices pre-allocated or precomputed) | ✓ VERIFIED | No `spzeros` calls in bohr_domain.jl; index-based accumulation loop at lines 17-23 (single-jump) and 44-50 (multi-jump); allocation test threshold set at num_freqs * dim^2 and passes |
| 2 | B_time and B_trotter closures compute phase rotations in-place without allocating Diagonal wrappers | ✓ VERIFIED | No `Diagonal(` calls in B_time/B_trotter functions (coherent.jl lines 163-348); pre-allocated diag_u/diag_u2 vector buffers used; allocation tests for both single/multi-jump variants pass with threshold 25 * d^2 |
| 3 | Time/Trotter thermalize hot path in jump_workers.jl avoids the abs.(filter(...)) intermediate allocation | ✓ VERIFIED | No `abs.(filter` pattern in jump_workers.jl; half-grid continue pattern at lines 416-440 (hermitian branch) with `w_raw > 1e-12 && continue`; allocation test passes with 50 * d^2 threshold |
| 4 | B_trotter multi-jump variant precomputes Trotter basis transforms instead of recomputing in the inner loop | ✓ VERIFIED | No `eigvecs' * jump` pattern in B_trotter (coherent.jl lines 253-348); all three callers (_precompute_coherent_total_B line 30, _precompute_coherent_unitary_terms line 85, _precompute_coherent_terms line 141) call transform_jumps_to_basis before B_trotter |
| 5 | All 224 existing tests pass with no regressions | ✓ VERIFIED | Test suite reports "231 tests passed" (224 original + 7 allocation tests) with no failures; numerical correctness preserved |

**Score:** 5/5 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `src/bohr_domain.jl` | Index-based accumulation B_bohr (single-jump and multi-jump) | ✓ VERIFIED | Lines 17-23 (single-jump) and 44-50 (multi-jump) contain `@inbounds for idx in indices` loops; no `spzeros` calls; pattern `for idx in indices` found |
| `src/jump_workers.jl` | Half-grid continue pattern replacing filter+abs | ✓ VERIFIED | Lines 416-440 implement `w_raw > 1e-12 && continue` pattern; no `abs.(filter` pattern found; hermitian and non-hermitian branches separated |
| `src/coherent.jl` | Diagonal-free B_time/B_trotter with pre-allocated vector buffers and redundant transform elimination | ✓ VERIFIED | Lines 171-172 (B_time single), 214-215 (B_time multi), 263-264 (B_trotter single), 307-308 (B_trotter multi) show `diag_u/diag_u2 = Vector{CT}(undef, d)`; no `Diagonal(` in function bodies; transform_jumps_to_basis calls at lines 30, 85, 141 |
| `test/test_allocation.jl` | Allocation regression tests for all four optimized hot paths | ✓ VERIFIED | 131 lines; contains @allocated assertions for B_bohr (lines 29, 42), B_time (lines 57, 69), B_trotter (lines 89, 96), _jump_contribution! (line 120) |
| `test/runtests.jl` | Updated test runner including allocation tests | ✓ VERIFIED | Line 16 contains `include("test_allocation.jl")` |

### Key Link Verification

| From | To | Via | Status | Details |
|------|-----|-----|--------|---------|
| src/bohr_domain.jl | src/coherent.jl | B_bohr called from _precompute_coherent_total_B and _precompute_coherent_unitary_terms | ✓ WIRED | Pattern `B_bohr(` found at coherent.jl lines 36, 96, 152 (total 3 call sites) |
| src/jump_workers.jl | test/test_thermalize.jl | _jump_contribution! called during thermalization steps | ✓ WIRED | Function used in thermalize.jl for both Time and Trotter domains; test suite passes |
| src/coherent.jl (_precompute_coherent_total_B) | src/coherent.jl (B_trotter) | Ensures jumps have Trotter-basis in_eigenbasis before B_trotter uses jump.in_eigenbasis directly | ✓ WIRED | Line 30: `trotter_jumps = transform_jumps_to_basis(jumps, ham_or_trott.eigvecs)` followed by `B_trotter(trotter_jumps, ...)` line 31 |
| src/coherent.jl (_precompute_coherent_unitary_terms) | src/coherent.jl (B_trotter) | Ensures jumps have Trotter-basis in_eigenbasis before B_trotter uses jump.in_eigenbasis directly | ✓ WIRED | Line 85: `trotter_jumps = transform_jumps_to_basis(jumps, trotter.eigvecs)` followed by loop calling `B_trotter(jump, trotter, ...)` line 87 |
| src/coherent.jl (_precompute_coherent_terms) | src/coherent.jl (B_trotter) | Ensures jumps have Trotter-basis in_eigenbasis before B_trotter uses jump.in_eigenbasis directly | ✓ WIRED | Line 141: `trotter_jumps = transform_jumps_to_basis(jumps, trotter.eigvecs)` followed by loop calling `B_trotter(jump, trotter, ...)` line 143 |
| test/test_allocation.jl | src/bohr_domain.jl | @allocated QuantumFurnace.B_bohr(..) | ✓ WIRED | Lines 28-29, 41-42 call B_bohr with @allocated measurement |
| test/test_allocation.jl | src/coherent.jl | @allocated QuantumFurnace.B_time(..) and @allocated QuantumFurnace.B_trotter(..) | ✓ WIRED | Lines 56-57, 68-69 (B_time); lines 88-89, 95-96 (B_trotter) |

### Requirements Coverage

No requirements explicitly mapped to Phase 11 in REQUIREMENTS.md. This is a performance optimization phase with observable success criteria defined in the phase goal.

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| src/coherent.jl | 351 | `#TODO: Reintroduce sigmas here` | ℹ️ Info | Pre-existing TODO, not related to Phase 11 work |
| src/log_sobolev.jl | 17 | `#TODO: Rewrite this with apply_lindbladian!()` | ℹ️ Info | Pre-existing TODO, not in modified files |
| src/jump_workers.jl | 464 | `#TODO: test it; set BLAS threads to 1` | ℹ️ Info | Pre-existing TODO, not in Phase 11 modified region |

**No blocker or warning-level anti-patterns found in Phase 11 modified code.**

### Human Verification Required

None. All verification criteria are programmatically verifiable through:
- Pattern absence checks (no spzeros, no Diagonal(), no abs.(filter))
- Pattern presence checks (index-based loops, diag_u buffers, transform_jumps_to_basis calls)
- Allocation regression tests with deterministic thresholds
- Full test suite pass (numerical correctness)

### Gaps Summary

None. All 5 observable truths verified, all artifacts substantive and wired, all key links connected, all tests pass, no blocker anti-patterns.

---

_Verified: 2026-02-15T16:00:00Z_
_Verifier: Claude (gsd-verifier)_
