---
status: testing
phase: 32-some-speedup-for-the-krylov-simulator
source: 32-01-SUMMARY.md, 32-02-SUMMARY.md
started: 2026-02-25T09:00:00Z
updated: 2026-02-25T09:00:00Z
---

## Current Test

number: 1
name: Matvec round-trip tests pass
expected: |
  `julia -e 'using Pkg; Pkg.test()' -- test/test_krylov_matvec.jl` runs all 198 matvec tests (round-trip, duality, zero-allocation) and they all pass. The optimized sandwich-only helpers produce identical results to the previous full-dissipator implementation.
awaiting: user response

## Tests

### 1. Matvec round-trip tests pass
expected: All 198 matvec round-trip, duality, and zero-allocation tests pass, confirming the 2+2N GEMM optimization is semantically identical to the previous 5N pattern
result: [pending]

### 2. Eigsolve tests pass
expected: All eigsolve tests pass (`test_krylov_eigsolve.jl`), confirming Krylov spectral gap computation still works correctly with the optimized matvec underneath
result: [pending]

### 3. KrylovWorkspace has precomputed G fields
expected: The KrylovWorkspace struct in `src/krylov_workspace.jl` contains fields G_left, G_right, G_left_adj, G_right_adj, and these are populated at construction time
result: [pending]

### 4. Legacy Euler apply_delta_channel! removed
expected: `src/krylov_eigsolve.jl` has NO 5-argument `apply_delta_channel!(ws, rho, delta, config, ham)` method. Only the faithful Chen 4-arg form remains.
result: [pending]

### 5. Dead _accumulate_dissipator! functions removed
expected: `src/krylov_matvec.jl` contains NO functions named `_accumulate_dissipator!`. Only `_accumulate_sandwich!` helpers remain.
result: [pending]

### 6. Cross-validation tests pass
expected: Cross-validation tests (`test_krylov_crossvalidation.jl`) for n=4 all-domain round-trip comparisons pass, confirming end-to-end Krylov gap accuracy is preserved
result: [pending]

## Summary

total: 6
passed: 0
issues: 0
pending: 6
skipped: 0

## Gaps

[none yet]
