# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-13)

**Core value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers
**Current focus:** Phase 2 complete - ready for Phase 3

## Current Position

Phase: 2 of 5 (Trajectory Bug Fixes) -- COMPLETE
Plan: 2 of 2 in current phase (COMPLETE)
Status: Phase 02 complete, Phase 03 next
Last activity: 2026-02-14 -- Executed 02-02-PLAN (CPTP completeness verification: TVAL-01)

Progress: [#####░░░░░] 50%

## Performance Metrics

**Velocity:**
- Total plans completed: 3
- Average duration: 5 min
- Total execution time: 0.25 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 01-foundation-and-compilation | 1 | 8 min | 8 min |
| 02-trajectory-bug-fixes | 2 | 7 min | 4 min |

**Recent Trend:**
- Last 5 plans: 8, 5, 2 min
- Trend: improving

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

### Pending Todos

None yet.

### Blockers/Concerns

- TFIX-01 (compilation bug) -- RESOLVED in Phase 1
- TFIX-02/03/04/05 -- RESOLVED in Phase 2 Plan 1
- TFIX-05: Cross-check confirmed jump sampling is structurally correct; U_B ordering (TFIX-02) was the root cause
- TVAL-01: CPTP completeness verified at 1e-10 for all three domains -- RESOLVED

## Session Continuity

Last session: 2026-02-14
Stopped at: Completed 02-02-PLAN.md (CPTP Completeness Verification)
Resume file: None
