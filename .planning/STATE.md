# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-25)

**Core value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers
**Current focus:** v2.0 Restructure

## Current Position

Phase: Not started (defining requirements)
Plan: —
Status: Defining requirements
Last activity: 2026-02-25 — Milestone v2.0 started

## Performance Metrics

**Velocity:**
- Total plans completed: 88 (v1.0: 10, v1.1: 16, quick: 24, v1.2: 12, cleanup: 3, v1.3: 10, v1.4: 2, v1.5: 12)

**By Milestone:**

| Milestone | Phases | Plans | Timeline |
|-----------|--------|-------|----------|
| v1.0 Trajectories | 1-5 | 10 | 2026-02-13 to 2026-02-14 |
| v1.1 Reduce | 6-11 | 16 (+5 quick) | 2026-02-15 |
| v1.2 Multi-threading | 12-19 | 15 (+3 quick) | 2026-02-15 to 2026-02-16 |
| v1.3 Mixing Time | 20-25 | 10 (+11 quick) | 2026-02-17 to 2026-02-18 |
| v1.4 Spectral Gap Refinement | 26 | 2 (+1 quick) | 2026-02-19 to 2026-02-20 |
| v1.5 Krylov Gap Estimation | 27-32 | 12 (+3 quick) | 2026-02-20 to 2026-02-25 |

## Accumulated Context

### Decisions

All decisions logged in PROJECT.md Key Decisions table.

### Deferred from v1.4

- FIT-01/02/03: Two-exponential fitting with Prony initialization
- RATE-01/02/03/04: Effective rate lambda_eff(t) and automatic window selection
- BOOT-01/02/03: Batched bootstrap uncertainty quantification
- RICH-01/02/03: Richardson extrapolation with monotonicity gate
- VAL-01/02/03/04: Diagnostic dashboard and final validation

### Deferred from v1.5

- BIEIG-01: bieigsolve for left+right eigenvectors simultaneously
- SECTOR-01: Sector-resolved gap computation
- SCALE-01: n=10, n=12 production runs on cluster
- ADAPT-01: Adaptive krylovdim auto-increase on partial convergence

### Pending Todos

None

### Blockers/Concerns

None

## Session Continuity

Last session: 2026-02-25
Stopped at: Completed v1.5 milestone archival
Resume file: None
