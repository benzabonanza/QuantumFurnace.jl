# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-13)

**Core value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers
**Current focus:** Phase 4 in progress - Trajectory Cross-Validation

## Current Position

Phase: 4 of 5 (Trajectory Cross-Validation)
Plan: 1 of 2 in current phase (Plan 1 COMPLETE)
Status: 04-01 complete, 04-02 next
Last activity: 2026-02-14 - Completed 04-01: Single-step trajectory-vs-DM cross-validation

Progress: [#######░░░] 70%

## Performance Metrics

**Velocity:**
- Total plans completed: 7
- Average duration: 5 min
- Total execution time: 0.63 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 01-foundation-and-compilation | 1 | 8 min | 8 min |
| 02-trajectory-bug-fixes | 2 | 7 min | 4 min |
| 03-dm-reference-test-suite | 3 | 17 min | 6 min |
| 04-trajectory-cross-validation | 1 | 6 min | 6 min |

**Recent Trend:**
- Last 5 plans: 2, 6, 5, 6, 6 min
- Trend: stable

*Updated after each plan completion*

## Accumulated Context

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- [Roadmap]: 5-phase structure: Foundation -> Bug Fixes -> DM Tests -> Trajectory Validation -> Regression
- [Roadmap]: TVAL-01 (CPTP check) grouped with bug fixes in Phase 2 -- it validates the fixes
- [Roadmap]: TINF-03 (Aqua) grouped with DM tests in Phase 3 -- natural first batch of real tests
- [Roadmap]: Phase 4 requires both Phase 2 and Phase 3 (trajectory code + DM reference)
- [01-01]: Relaxed stdlib compat bounds to support Julia 1.11+ (original were 1.12-only)
- [01-01]: Test Hamiltonian loaded via @__DIR__ path, not load_hamiltonian(), for Pkg.test() compatibility
- [01-01]: TOL_DELTA constant C=5.0, TEST_DELTA=0.01 for unraveling error tolerance
- [02-01]: EnergyDomain had U_B ordering bug; Time/Trotter variant was already correct
- [02-01]: Used TrottTrott() constructor directly -- create_trotter() is exported but never defined
- [02-01]: Replaced Cholesky with eigendecomposition for PSD guard (silent clamp to zero)
- [02-02]: CPTP test in its own file (test_cptp.jl), separate from bug fix tests
- [02-02]: Tolerance 1e-10 for CPTP completeness check (allows small numerical accumulation)
- [03-03]: Disabled Aqua ambiguities/piracies checks (legitimate multiple dispatch + kron! extension)
- [03-03]: Removed 4 undefined exports (create_trotter, create_hamham, evolve_along_trajectory, evolve_and_measure_along_trajectory)
- [03-03]: Added comprehensive compat bounds for all deps/extras/julia to satisfy Aqua quality gate
- [03-01]: Used eigen() instead of Arpack eigs for small dense Liouvillians (deterministic, no convergence issues)
- [03-01]: QuantumFurnace.trace_distance_h qualified access since function not exported
- [03-01]: Hierarchy tolerance 1e-12 for numerical noise in domain distance comparisons
- [03-02]: EnergyDomain Liouvillian for Euler scaling tests (simplest non-trivial domain)
- [03-02]: Maximally mixed initial state for clean scaling measurement
- [03-02]: Non-exported functions (time_oft!, trotter_oft!, etc.) accessed via QuantumFurnace.func()
- [04-01]: 50,000 trajectories per delta point (10k insufficient -- noise floor ~0.01 trace distance)
- [04-01]: Delta sweep [0.2, 0.1, 0.05] (delta=0.025 falls below statistical noise floor)
- [04-01]: Ratio bounds [2.0, 8.0] for O(delta^2) check (accommodates statistical fluctuations)

### Pending Todos

- [ ] Refactor step_along_trajectory to per-operator branching (trajectories) — [todo](./todos/pending/2026-02-14-refactor-step-along-trajectory-to-per-operator-branching.md)

### Quick Tasks Completed

| # | Description | Date | Commit | Directory |
|---|-------------|------|--------|-----------|
| 1 | Fix TrotterDomain eigenbasis mismatch in DMTST-02 | 2026-02-14 | 8ad104f | [1-fix-trotterdomain-eigenbasis-mismatch-in](./quick/1-fix-trotterdomain-eigenbasis-mismatch-in/) |
| 2 | Fix Gibbs basis transformation for TrotterDomain in DMTST-02 | 2026-02-14 | 7a117c7 | [2-fix-gibbs-basis-transformation-for-trott](./quick/2-fix-gibbs-basis-transformation-for-trott/) |
| 3 | Fix OFT consistency test basis transformation (DMTST-06) | 2026-02-14 | c3b94f3 | [3-fix-oft-consistency-test-basis-transform](./quick/3-fix-oft-consistency-test-basis-transform/) |
| 4 | Add NUFFT OFT consistency test (DMTST-06b) | 2026-02-14 | c4a2126 | [4-add-nufft-oft-consistency-test-alongside](./quick/4-add-nufft-oft-consistency-test-alongside/) |
| 5 | Tighten TrotterDomain error thresholds | 2026-02-14 | eb2569b | [5-tighten-trotterdomain-error-thresholds-f](./quick/5-tighten-trotterdomain-error-thresholds-f/) |
| 6 | Fix B_trotter() basis mismatch | 2026-02-14 | 2334a5e | [6-fix-b-time-vs-b-trott-basis-mismatch-in-](./quick/6-fix-b-time-vs-b-trott-basis-mismatch-in-/) |

### Blockers/Concerns

- TFIX-01 (compilation bug) -- RESOLVED in Phase 1
- TFIX-02/03/04/05 -- RESOLVED in Phase 2 Plan 1
- TFIX-05: Cross-check confirmed jump sampling is structurally correct; U_B ordering (TFIX-02) was the root cause
- TVAL-01: CPTP completeness verified at 1e-10 for all three domains -- RESOLVED

## Session Continuity

Last session: 2026-02-14
Stopped at: Completed 04-01-PLAN.md (single-step trajectory cross-validation)
Resume file: None
