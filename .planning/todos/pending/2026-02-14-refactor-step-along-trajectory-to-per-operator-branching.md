---
created: 2026-02-14T12:02:38.471Z
title: Refactor step_along_trajectory to per-operator branching
area: trajectories
files:
  - src/trajectories.jl
  - src/QuantumFurnace.jl
---

## Problem

The current `step_along_trajectory!()` function logic does not exactly match the DM simulator's per-operator channel structure. It uses total/aggregated R, S, U_residual quantities rather than per-Lindblad-operator (per `a`) quantities. While the error scaling may match for some delta values, the mathematical structure should more faithfully reproduce Chen's quantum algorithm.

The correct conceptual flow per step should be:

1. **Precompute per-operator quantities:** R^a, S^a, U_residual^a for each Lindblad operator index `a`
2. **Coherent evolution:** If `with_coherent`, evolve deterministically with exp(-i delta B^a)
3. **Branching probabilities:** Compute probabilities with respect to per-operator R^a, S^a, U_residual^a (not the current total R)
4. **Branch selection:** No-jump, residual term, or jump — when the jump branch is selected for a given `a`, sample an omega-contribution of the jump
5. This should mathematically give the same channel as Chen's quantum algorithm or the DM simulator on average

## Solution

Refactor `step_along_trajectory!()` to:
- Store per-`a` precomputed R^a, S^a, U_residual^a in the TrajectoryFramework
- Implement two-level branching: first select operator `a`, then branch within that operator
- Ensure the resulting channel matches `jump_contribution!` in the DM code exactly
- Verify with CPTP test that the per-operator channel still sums to identity

This may need to happen before or during Phase 4 cross-validation to get meaningful trajectory-vs-DM agreement.
