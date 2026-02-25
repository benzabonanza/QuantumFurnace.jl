---
phase: 32-some-speedup-for-the-krylov-simulator
plan: 02
subsystem: krylov-simulator
tags: [dead-code, cleanup, dissipator, matvec, eigsolve]

# Dependency graph
requires:
  - phase: 32-some-speedup-for-the-krylov-simulator
    plan: 01
    provides: "Sandwich-only helpers replacing full dissipator functions; G_left/G_right precomputation"
provides:
  - "Clean krylov_matvec.jl with only sandwich-only helpers (no legacy 5-GEMM dissipator functions)"
  - "Clean krylov_eigsolve.jl with only faithful Chen 4-arg apply_delta_channel! (no legacy Euler)"
  - "Net 309 lines of dead code removed"
affects: [krylov-benchmark, krylov-eigsolve]

# Tech tracking
tech-stack:
  added: []
  patterns: []

key-files:
  created: []
  modified:
    - "src/krylov_eigsolve.jl"
    - "src/krylov_matvec.jl"
    - "test/test_krylov_eigsolve.jl"

key-decisions:
  - "No new decisions required -- pure dead code removal following Plan 32-01 optimization"

patterns-established: []

# Metrics
duration: 7min
completed: 2026-02-25
---

# Phase 32 Plan 02: Dead Code Removal Summary

**Deleted legacy Euler apply_delta_channel!, its test, and 6 dead _accumulate_dissipator! functions (309 lines) superseded by Plan 32-01 sandwich-only optimization**

## Performance

- **Duration:** 7 min
- **Started:** 2026-02-25T08:26:31Z
- **Completed:** 2026-02-25T08:33:51Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments
- Deleted the legacy 5-argument Euler `apply_delta_channel!` and its testset from krylov_eigsolve.jl
- Deleted all 6 `_accumulate_dissipator!` family functions (each 5 GEMMs) from krylov_matvec.jl
- Net 309 lines removed with zero regressions (198 matvec + 42 eigsolve tests pass)

## Task Commits

Each task was committed atomically:

1. **Task 1: Delete legacy Euler apply_delta_channel! and its test** - `1473dea` (refactor)
2. **Task 2: Delete dead _accumulate_dissipator! family (6 functions)** - `5d5c3a8` (refactor)

## Files Created/Modified
- `src/krylov_eigsolve.jl` - Removed 5-arg legacy Euler apply_delta_channel! (22 lines); only faithful Chen 4-arg form remains
- `src/krylov_matvec.jl` - Removed 6 dead _accumulate_dissipator! functions (267 lines); only _accumulate_sandwich! helpers remain
- `test/test_krylov_eigsolve.jl` - Removed "apply_delta_channel! legacy Euler" testset (20 lines)

## Decisions Made
None - followed plan as specified. Pure dead code removal.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

Pre-existing test failures in `test_krylov_crossvalidation.jl` for TimeDomain/TrotterDomain L-vs-E convergence order (`order >= 1.5` assertion). These are unrelated to the dead code removal -- the crossvalidation file does not reference any deleted functions. All directly relevant tests (198 matvec + 42 eigsolve) pass cleanly.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- Phase 32 is complete: precomputed G_left/G_right optimization (Plan 01) + dead code cleanup (Plan 02)
- krylov_matvec.jl is clean: only _krylov_oft!, 6 sandwich helpers, 2 sandwich_2op helpers, and 8 apply_lindbladian!/apply_adjoint_lindbladian! functions
- krylov_eigsolve.jl is clean: only faithful Chen channel apply_delta_channel! and eigsolve API

## Self-Check: PASSED

All artifacts verified:
- src/krylov_eigsolve.jl: FOUND, no legacy Euler apply_delta_channel!
- src/krylov_matvec.jl: FOUND, no _accumulate_dissipator functions
- test/test_krylov_eigsolve.jl: FOUND, no legacy Euler testset
- Commit 1473dea: FOUND (Task 1)
- Commit 5d5c3a8: FOUND (Task 2)

---
*Phase: 32-some-speedup-for-the-krylov-simulator*
*Completed: 2026-02-25*
