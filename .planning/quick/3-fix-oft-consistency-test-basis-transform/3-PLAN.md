---
phase: quick
plan: 3
subsystem: testing
tags: [trotter, eigenbasis, oft, basis-transformation]
---

# Quick Task 3: Fix OFT Consistency Test Basis Transformation for TrotterDomain (DMTST-06)

<objective>
Fix the DMTST-06 OFT consistency test to correctly compare trotter_oft!() output against oft!() and time_oft!() by handling the Trotter eigenbasis transformation.
</objective>

## Problem

`trotter_oft!()` computes `A(w) = sum_n prefactor(n) * U^n * A * U^{-n}` where U is diagonal with Trotter eigenvalues. For this element-wise product to be correct, the jump operator must be in Trotter eigenbasis. The test was passing `jump.in_eigenbasis` (H-eigenbasis), producing a mixed-basis output that can't be meaningfully compared.

Also discovered: the Explore agent accidentally removed `trafo_from_eigen_to_trotter` from the TrottTrott struct, breaking DMTST-05. This was restored via `git checkout`.

## Tasks

### Task 1: Transform jump to Trotter eigenbasis, compute OFT, transform back
**File:** `test/test_dm_scaling.jl`

1. Transform jump to Trotter eigenbasis: `U * jump.in_eigenbasis * U'` where `U = trafo_from_eigen_to_trotter`
2. Call `trotter_oft!` with the Trotter-basis jump
3. Transform output back to H-eigenbasis: `U' * A_trott * U`
4. Add sanity assertion `dist < 0.1`

**Acceptance:** Trotter OFT error < 1e-6. All tests pass.
