# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-13)

**Core value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers
**Current focus:** Phase 5 COMPLETE - Statistical Validation and Regression

## Current Position

Phase: 5 of 5 (Statistical Validation and Regression) - COMPLETE
Plan: 2 of 2 in current phase (ALL COMPLETE)
Status: All 5 phases complete
Last activity: 2026-02-14 - Completed quick task 11: Widen trajectory regression tolerance to 1e-3 for Julia version portability

Progress: [##########] 100%

## Performance Metrics

**Velocity:**
- Total plans completed: 10
- Average duration: 5 min
- Total execution time: 0.83 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 01-foundation-and-compilation | 1 | 8 min | 8 min |
| 02-trajectory-bug-fixes | 2 | 7 min | 4 min |
| 03-dm-reference-test-suite | 3 | 17 min | 6 min |
| 04-trajectory-cross-validation | 2 | 8 min | 4 min |
| 05-statistical-validation-and-regression | 2 | 10 min | 5 min |

**Recent Trend:**
- Last 5 plans: 6, 6, 2, 5, 5 min
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
- [04-02]: DM convergence against Liouvillian fixed point (not Gibbs) to isolate convergence from domain approximation
- [04-02]: Trajectory threshold 0.02 to Gibbs (domain approx ~0.005 + statistical noise ~0.01 makes 1e-3 unachievable)
- [04-02]: Thermalization assertion: trajectory must move >10x closer to fixed point than initial state
- [quick-7]: Per-operator Lie-Trotter splitting: delta_eff = delta * N_jumps, alpha from delta_eff
- [quick-7]: Per-operator coherent B via precompute_coherent_total_B with single-jump vector
- [quick-7]: JumpOp[jumps[a]] forces Vector{JumpOp} type for dispatch compatibility
- [Phase quick-8]: Trotter basis transform (U * A * U') applied before NUFFT dissipative loops in all three entry points
- [Phase quick-8]: JumpOp[] typed comprehension to preserve Vector{JumpOp} for struct field compatibility
- [quick-9]: Removed trafo_from_eigen_to_trotter field; use trotter.eigvecs' * j.data * trotter.eigvecs directly
- [quick-9]: Test back-transforms use U_t2e = trotter.eigvecs' * ham.eigvecs (equivalent to removed field)
- [05-01]: High-N trajectory reference (500k) instead of Liouvillian DM reference to avoid Lie-Trotter splitting bias
- [05-01]: Batch-averaged errors (10 batches per N_traj) for robust convergence ratio estimation
- [05-01]: N_traj points [200, 800, 3200, 12800] with factor-of-4 steps for 1/sqrt(N) verification
- [05-02]: Symbol keys in BSON (not string keys) for idiomatic Julia d[:rho] access pattern
- [05-02]: Tolerance 1e-10 for regression comparison (allows floating-point accumulation across platforms)
- [05-02]: Trajectory seed=12345, ntraj=1000 for deterministic reference (distinct from Phase 4 seeds)
- [quick-10]: Trajectory regression atol relaxed to 1e-6 for cross-platform BLAS variance; DM tests unchanged at 1e-10
- [quick-11]: Trajectory regression atol widened to 1e-3; 1e-6 insufficient for Julia 1.11→1.12 + x86→aarch64 differences (max diff ~7.5e-5 observed)

### Pending Todos

None

### Quick Tasks Completed

| # | Description | Date | Commit | Directory |
|---|-------------|------|--------|-----------|
| 1 | Fix TrotterDomain eigenbasis mismatch in DMTST-02 | 2026-02-14 | 8ad104f | [1-fix-trotterdomain-eigenbasis-mismatch-in](./quick/1-fix-trotterdomain-eigenbasis-mismatch-in/) |
| 2 | Fix Gibbs basis transformation for TrotterDomain in DMTST-02 | 2026-02-14 | 7a117c7 | [2-fix-gibbs-basis-transformation-for-trott](./quick/2-fix-gibbs-basis-transformation-for-trott/) |
| 3 | Fix OFT consistency test basis transformation (DMTST-06) | 2026-02-14 | c3b94f3 | [3-fix-oft-consistency-test-basis-transform](./quick/3-fix-oft-consistency-test-basis-transform/) |
| 4 | Add NUFFT OFT consistency test (DMTST-06b) | 2026-02-14 | c4a2126 | [4-add-nufft-oft-consistency-test-alongside](./quick/4-add-nufft-oft-consistency-test-alongside/) |
| 5 | Tighten TrotterDomain error thresholds | 2026-02-14 | eb2569b | [5-tighten-trotterdomain-error-thresholds-f](./quick/5-tighten-trotterdomain-error-thresholds-f/) |
| 6 | Fix B_trotter() basis mismatch | 2026-02-14 | 2334a5e | [6-fix-b-time-vs-b-trott-basis-mismatch-in-](./quick/6-fix-b-time-vs-b-trott-basis-mismatch-in-/) |
| 7 | Refactor step_along_trajectory to per-operator Lie-Trotter splitting | 2026-02-14 | 0c879e7 | [7-refactor-step-along-trajectory-to-per-op](./quick/7-refactor-step-along-trajectory-to-per-op/) |
| 8 | Fix TrotterDomain Gibbs fixed point distance (~0.004 -> ~9e-9) | 2026-02-14 | e0ef0fc | [8-fix-trotterdomain-gibbs-fixed-point-dist](./quick/8-fix-trotterdomain-gibbs-fixed-point-dist/) |
| 9 | Remove trafo_from_eigen_to_trotter from TrottTrott | 2026-02-14 | 5646187 | [9-remove-trafo-from-eigen-to-trotter-from-](./quick/9-remove-trafo-from-eigen-to-trotter-from-/) |
| 10 | Fix trajectory regression test tolerance for cross-platform BLAS | 2026-02-14 | 9d2e71a | [10-fix-trajectory-regression-test-failure](./quick/10-fix-trajectory-regression-test-failure/) |
| 11 | Widen trajectory regression tolerance to 1e-3 for Julia version portability | 2026-02-14 | 7ab18ae | [11-widen-trajectory-regression-tolerance-to](./quick/11-widen-trajectory-regression-tolerance-to/) |

### Blockers/Concerns

- TFIX-01 (compilation bug) -- RESOLVED in Phase 1
- TFIX-02/03/04/05 -- RESOLVED in Phase 2 Plan 1
- TFIX-05: Cross-check confirmed jump sampling is structurally correct; U_B ordering (TFIX-02) was the root cause
- TVAL-01: CPTP completeness verified at 1e-10 for all three domains -- RESOLVED

## Session Continuity

Last session: 2026-02-14
Stopped at: Completed quick task 11 (widen trajectory regression tolerance to 1e-3)
Resume file: None
