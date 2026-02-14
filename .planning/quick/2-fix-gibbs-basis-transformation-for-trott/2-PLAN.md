---
phase: quick
plan: 2
subsystem: testing
tags: [trotter, eigenbasis, basis-transformation, detailed-balance]
---

# Quick Task 2: Fix Gibbs Basis Transformation for TrotterDomain in DMTST-02

<objective>
Correct the basis transformation when comparing TrotterDomain steady state to Gibbs state in DMTST-02. The previous fix (quick task 1) used the wrong intermediate basis.
</objective>

## Problem

Quick task 1 applied `trotter.eigvecs' * TEST_GIBBS * trotter.eigvecs`, but:
- `TEST_GIBBS` is in **Hamiltonian eigenbasis** (diagonal)
- `trotter.eigvecs` transforms between **computational basis** and **Trotter eigenbasis**

These are three different bases. The correct path is: eigenbasis -> computational -> Trotter eigenbasis.

## Tasks

### Task 1: Two-step basis transformation
**File:** `test/test_dm_detailed_balance.jl`

Replace single transformation with two-step:
```julia
gibbs_comp = TEST_HAM.eigvecs * TEST_GIBBS * TEST_HAM.eigvecs'
gibbs_ref = Hermitian(TEST_TROTTER.eigvecs' * gibbs_comp * TEST_TROTTER.eigvecs)
```

**Acceptance:** Trotter distance drops from 0.688 to ~1e-9 or better. All 69 tests pass.
