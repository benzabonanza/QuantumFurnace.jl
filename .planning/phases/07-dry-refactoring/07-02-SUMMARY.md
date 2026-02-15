---
phase: 07-dry-refactoring
plan: 02
subsystem: core-simulation
tags: [julia, DRY, CPTP, coherent-unitary, refactoring, jump-workers]

# Dependency graph
requires:
  - phase: 07-dry-refactoring-01
    provides: hermitianize! helper in qi_tools.jl used by CPTP channel and R accumulation
provides:
  - apply_coherent_unitary! in-place helper in jump_workers.jl
  - apply_cptp_channel! in-place helper in jump_workers.jl
  - all 3 inline coherent unitary blocks replaced
  - all 3 inline CPTP channel application blocks replaced
affects: [08-naming-conventions]

# Tech tracking
tech-stack:
  added: []
  patterns: [apply_coherent_unitary! for U_B sandwich transforms, apply_cptp_channel! for Chen Eq. 3.2 weak-measurement CPTP step]

key-files:
  created: []
  modified:
    - src/jump_workers.jl

key-decisions:
  - "apply_cptp_channel! expects scratch.R pre-Hermitianized; hermitianize!(scratch.R) remains before call site"
  - "apply_coherent_unitary! marked @inline for zero-overhead dispatch on nothing vs Matrix"

patterns-established:
  - "apply_coherent_unitary!(evolving_dm, U_B, scratch): canonical coherent unitary sandwich transform"
  - "apply_cptp_channel!(evolving_dm, delta, scratch): canonical CPTP weak-measurement channel (K0 + residual Cholesky)"

# Metrics
duration: 4min
completed: 2026-02-15
---

# Phase 7 Plan 2: CPTP Channel and Coherent Unitary Helpers Summary

**Extracted apply_coherent_unitary! and apply_cptp_channel! helpers, replacing 3 identical 5-line unitary blocks and 3 identical 30-line CPTP blocks in jump_workers.jl**

## Performance

- **Duration:** 4 min
- **Started:** 2026-02-15T07:24:35Z
- **Completed:** 2026-02-15T07:29:01Z
- **Tasks:** 2
- **Files modified:** 1

## Accomplishments
- Created `apply_coherent_unitary!(evolving_dm, U_B, scratch)` helper replacing 3 identical inline U_B sandwich blocks
- Created `apply_cptp_channel!(evolving_dm, delta, scratch)` helper replacing 3 identical inline K0/residual/Cholesky CPTP blocks (~30 lines each)
- Net reduction of ~70 lines (22 insertions + 61 insertions - 21 deletions - 113 deletions)
- All 224 tests pass with zero regressions after each task

## Task Commits

Each task was committed atomically:

1. **Task 1: Extract apply_coherent_unitary! helper from 3 identical blocks** - `eaf5486` (feat)
2. **Task 2: Extract apply_cptp_channel! helper from 3 identical blocks** - `bd95b34` (feat)

## Files Created/Modified
- `src/jump_workers.jl` - Added apply_coherent_unitary! and apply_cptp_channel! helpers; replaced 3+3 inline blocks in BohrDomain, EnergyDomain, Time/TrotterDomain jump_contribution! methods

## Decisions Made
- apply_cptp_channel! expects scratch.R to be Hermitianized before the call; the hermitianize!(scratch.R) call remains at each call site (not inside the helper) to keep the helper focused on the CPTP step
- apply_coherent_unitary! is marked @inline since the nothing-check fast path should be zero-cost

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- DRY-02 (CPTP channel) and DRY-03 (coherent unitary) satisfied
- Phase 07 DRY refactoring complete (both plans done)
- Ready for Phase 08 (naming conventions)

## Self-Check: PASSED

All files verified present. Commits eaf5486 and bd95b34 verified in git log.

---
*Phase: 07-dry-refactoring*
*Completed: 2026-02-15*
