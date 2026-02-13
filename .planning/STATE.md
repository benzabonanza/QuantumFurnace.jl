# State: QuantumFurnace.jl

## Current Position

Phase: Not started (defining requirements)
Plan: —
Status: Defining requirements
Last activity: 2026-02-13 — Milestone v1.0 Trajectories started

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-13)

**Core value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers
**Current focus:** Milestone v1.0 — Fix, validate, and test trajectory simulation

## Accumulated Context

- Trajectory sampling has a known algorithmic bug: jumps flatten all (a, omega) pairs instead of two-stage selection (pick A^a, then sample omega)
- Error hierarchy: Bohr (exact) → Energy (+quadrature) → Time (+time quadrature) → Trotter (+Trotter error)
- Even TrotterDomain with_coherent=true should achieve ≤ 1e-6 error to Gibbs with proper parameters
- Existing configs in main_liouv.jl and main_thermalize.jl are known-good parameter sets
- Test systems: 3-4 qubit Heisenberg (1-2 qubit may hit edge cases)
- Single-threaded trajectories for this milestone; parallelism deferred
