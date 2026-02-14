---
phase: quick
plan: 2
subsystem: testing
tags: [trotter, eigenbasis, basis-transformation, detailed-balance]
---

# Summary: Fix Gibbs Basis Transformation for TrotterDomain in DMTST-02

## What Changed

**File:** `test/test_dm_detailed_balance.jl` (lines 65-72)

Corrected the Gibbs state basis transformation for TrotterDomain comparison. The previous fix (quick task 1) directly applied `trotter.eigvecs' * TEST_GIBBS * trotter.eigvecs`, which mixed Hamiltonian eigenbasis with the computational-to-Trotter transformation. The correct path goes through computational basis as intermediate:

```julia
gibbs_comp = TEST_HAM.eigvecs * TEST_GIBBS * TEST_HAM.eigvecs'
gibbs_ref = Hermitian(TEST_TROTTER.eigvecs' * gibbs_comp * TEST_TROTTER.eigvecs)
```

## Results

Trotter steady-state distance to Gibbs: **9.3e-9** (was 0.688 with wrong basis)

Full hierarchy now verified:
- bohr: 3.4e-16
- energy: 2.3e-14
- time: 6.3e-14
- trotter: 9.3e-9

All 69 tests pass.

## Duration

< 1 min
