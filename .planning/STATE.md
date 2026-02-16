# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-16)

**Core value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers
**Current focus:** v1.3 Mixing Time Estimation -- Phase 20 (Observable Infrastructure)

## Current Position

Phase: 20 of 24 (Observable Infrastructure)
Plan: 0 of ? in current phase
Status: Ready to plan
Last activity: 2026-02-16 -- Roadmap created for v1.3 Mixing Time Estimation (5 phases, 17 requirements)

Progress: [####################..........] 41/46 plans (v1.0-v1.2 complete, v1.3 not started)

## Performance Metrics

**Velocity:**
- Total plans completed: 46 (v1.0: 10, v1.1: 16, quick: 8, v1.2: 12, cleanup: 3)

**By Milestone:**

| Milestone | Phases | Plans | Timeline |
|-----------|--------|-------|----------|
| v1.0 Trajectories | 1-5 | 10 | 2026-02-13 to 2026-02-14 |
| v1.1 Reduce | 6-11 | 16 (+5 quick) | 2026-02-15 |
| v1.2 Multi-threading | 12-19 | 15 (+3 quick) | 2026-02-15 to 2026-02-16 |
| v1.3 Mixing Time | 20-24 | TBD | 2026-02-16 to ... |

## Accumulated Context

### Decisions

All decisions logged in PROJECT.md Key Decisions table.
Key context for v1.3:
- Use `abs(real(spectral_gap))` not `abs(spectral_gap)` for cross-validation (complex eigenvalues)
- Follow `build_convergence_observables` pattern for basis transforms (avoids pitfall from quick-task-20)
- LsqFit.jl is the single new dependency (Levenberg-Marquardt + CIs + bounds)
- Phases 20 and 21 can execute in parallel (independent work)

### Pending Todos

None

### Blockers/Concerns

None

### Quick Tasks Completed

| # | Description | Date | Commit | Directory |
|---|-------------|------|--------|-----------|
| 19 | Fix failing test after EnergyDomain to TrotterDomain rename in GNS trajectory | 2026-02-16 | 0683acc | [19-fix-failing-test-after-energydomain-to-t](./quick/19-fix-failing-test-after-energydomain-to-t/) |
| 20 | Fix basis mismatch in GNS TrotterDomain tests (0.83 -> 0.08 gap) | 2026-02-16 | 0161e01 | [20-debug-gns-trotterdomain-0-83-gap-suspect](./quick/20-debug-gns-trotterdomain-0-83-gap-suspect/) |
| 21 | Fix test errors after removing transform_jumps_to_basis | 2026-02-16 | b2e4123 | [21-fix-test-errors-after-removing-transform](./quick/21-fix-test-errors-after-removing-transform/) |

## Session Continuity

Last session: 2026-02-16
Stopped at: v1.3 roadmap created, ready to plan Phase 20 or Phase 21
Resume file: None
