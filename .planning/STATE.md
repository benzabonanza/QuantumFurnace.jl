# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-13)

**Core value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers
**Current focus:** Phase 3 in progress - DM Reference Test Suite

## Current Position

Phase: 3 of 5 (DM Reference Test Suite) -- IN PROGRESS
Plan: 3 of 3 in current phase (03-03 COMPLETE)
Status: Phase 03 plan 03 complete (Aqua.jl package quality)
Last activity: 2026-02-14 -- Executed 03-03-PLAN (Aqua.jl package quality: TINF-03)

Progress: [######░░░░] 60%

## Performance Metrics

**Velocity:**
- Total plans completed: 4
- Average duration: 5 min
- Total execution time: 0.35 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 01-foundation-and-compilation | 1 | 8 min | 8 min |
| 02-trajectory-bug-fixes | 2 | 7 min | 4 min |
| 03-dm-reference-test-suite | 1 | 6 min | 6 min |

**Recent Trend:**
- Last 5 plans: 8, 5, 2, 6 min
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

### Pending Todos

None yet.

### Blockers/Concerns

- TFIX-01 (compilation bug) -- RESOLVED in Phase 1
- TFIX-02/03/04/05 -- RESOLVED in Phase 2 Plan 1
- TFIX-05: Cross-check confirmed jump sampling is structurally correct; U_B ordering (TFIX-02) was the root cause
- TVAL-01: CPTP completeness verified at 1e-10 for all three domains -- RESOLVED

## Session Continuity

Last session: 2026-02-14
Stopped at: Completed 03-03-PLAN.md (Aqua.jl Package Quality)
Resume file: None
