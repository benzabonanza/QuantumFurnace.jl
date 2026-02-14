---
phase: quick
plan: 3
subsystem: testing
tags: [trotter, eigenbasis, oft, basis-transformation]
---

# Summary: Fix OFT Consistency Test Basis Transformation (DMTST-06)

## What Changed

**File:** `test/test_dm_scaling.jl` (DMTST-06 testset)

1. Jump operator is now transformed to Trotter eigenbasis before calling `trotter_oft!()`:
   `U * jump.in_eigenbasis * U'` where `U = trafo_from_eigen_to_trotter`
2. Output is transformed back to H-eigenbasis: `U' * A_trott * U`
3. Added sanity assertion `dist < 0.1`

**Also fixed:** Restored `trafo_from_eigen_to_trotter` field in TrottTrott struct that was accidentally removed by an exploration agent.

## Why

`trotter_oft!()` performs element-wise multiplication of the jump matrix with Trotter eigenvalue phases (`lambda_i^n * conj(lambda_j^n)`). For this to be correct, both the jump and the phases must use the same basis indices — i.e., the jump must be in Trotter eigenbasis. Previously the H-eigenbasis jump was passed, producing a mixed-basis output.

This is the same pattern as `B_trotter()` in DMTST-05 where the function internally uses Trotter eigenvalue phases with `jump.in_eigenbasis`, but `B_trotter` uses matrix products (not element-wise), making the basis handling implicit.

## Results

Trotter OFT error: **1.47e-8** (was 1.40 before fix)

All 70 tests pass (69 original + 1 new sanity check).

## Duration

~15 min (investigation of basis conventions + debugging accidental struct deletion)
