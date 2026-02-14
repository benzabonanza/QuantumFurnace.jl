---
phase: quick
plan: 1
subsystem: testing
tags: [trotter, eigenbasis, detailed-balance, domain-hierarchy]
---

# Quick Task 1: Fix TrotterDomain Eigenbasis Mismatch in DMTST-02

<objective>
Fix the DMTST-02 domain error hierarchy test to correctly compare TrotterDomain steady state against Gibbs state in the Trotter eigenbasis.
</objective>

## Problem

The TrotterDomain Liouvillian operates in the Trotter eigenbasis (the eigenbasis of the Trotterized Hamiltonian), not the exact Hamiltonian eigenbasis. When comparing the steady state to the Gibbs state, the Gibbs state must be transformed into the Trotter eigenbasis first. Without this transformation, the trace distance is artificially large (~0.69), breaking the monotonic hierarchy check.

## Tasks

### Task 1: Transform Gibbs reference for TrotterDomain comparison
**File:** `test/test_dm_detailed_balance.jl`

In the DMTST-02 test loop, when the domain is `TrotterDomain`, transform `TEST_GIBBS` into the Trotter eigenbasis via `TEST_TROTTER.eigvecs' * TEST_GIBBS * TEST_TROTTER.eigvecs` before computing trace distance.

**Acceptance:** All 69 tests pass, including the hierarchy check `dist(Time) <= dist(Trotter)`.
