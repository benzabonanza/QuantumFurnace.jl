---
phase: 34-code-deduplication
plan: 02
subsystem: simulation-core
tags: [julia, BLAS, deduplication, CPTP, sandwich-helpers, Chen-channel]

# Dependency graph
requires:
  - phase: 34-01
    provides: "domain_prefactor and unified oft! (foundation for cross-simulation-type dedup)"
  - phase: quick-39
    provides: "Fixed krylov_matvec sandwich convention to L*rho*L'"
provides:
  - "Consolidated sandwich helpers: _accumulate_sandwich! (L*rho*L') and _accumulate_sandwich_adj! (L'*rho*L)"
  - "_build_cptp_channel(R, delta) shared helper for Chen Eq. 3.2 CPTP channel construction"
affects: [Phase 35, Phase 36]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Two canonical sandwich helpers replace four pairwise-identical functions"
    - "Single _build_cptp_channel function replaces three inline CPTP constructions"

key-files:
  created: []
  modified:
    - src/krylov_matvec.jl
    - src/furnace_utensils.jl
    - src/krylov_workspace.jl
    - src/trajectories.jl
    - src/jump_workers.jl

key-decisions:
  - "_accumulate_sandwich_adj! chosen as canonical name for L'*rho*L (clear, symmetric with _accumulate_sandwich!)"
  - "_build_cptp_channel returns NamedTuple (; K0, U_residual, alpha) -- callers destructure what they need"
  - "scratch.K0 field left in KrausScratch struct (dead after extraction) -- struct change deferred to Phase 35"

patterns-established:
  - "_accumulate_sandwich! (L*rho*L') and _accumulate_sandwich_adj! (L'*rho*L) are the only single-operator sandwich helpers"
  - "_build_cptp_channel(R, delta) is the single source for Chen Eq. 3.2 CPTP construction"

# Metrics
duration: 9min
completed: 2026-02-26
---

# Phase 34 Plan 02: Consolidate Sandwich Helpers and Extract CPTP Channel Summary

**4 sandwich helpers consolidated to 2, plus _build_cptp_channel extracting Chen Eq. 3.2 CPTP formula from 3 files into furnace_utensils.jl**

## Performance

- **Duration:** 9 min
- **Started:** 2026-02-26T13:12:51Z
- **Completed:** 2026-02-26T13:21:38Z
- **Tasks:** 2
- **Files modified:** 5

## Accomplishments
- Consolidated 4 pairwise-identical sandwich helpers into 2 canonical functions (`_accumulate_sandwich!` for L*rho*L', `_accumulate_sandwich_adj!` for L'*rho*L) with `@inline` for hot path
- Extracted `_build_cptp_channel(R, delta)` into furnace_utensils.jl as single source for Chen Eq. 3.2 CPTP weak-measurement channel construction
- Replaced 3 independent inline CPTP formulas (krylov_workspace.jl, trajectories.jl, jump_workers.jl) with calls to the shared helper
- All 1198 tests pass with identical numerical results

## Task Commits

Each task was committed atomically:

1. **Task 1: Consolidate 4 sandwich helpers into 2 in krylov_matvec.jl** - `a1b5363` (feat)
2. **Task 2: Extract _build_cptp_channel shared helper** - `aab8d8b` (feat)

## Files Created/Modified
- `src/krylov_matvec.jl` - Consolidated 4 sandwich helpers to 2; added @inline; updated all 8 call sites across Energy/Time/Trotter forward+adjoint methods
- `src/furnace_utensils.jl` - Added _build_cptp_channel(R, delta) returning (; K0, U_residual, alpha)
- `src/krylov_workspace.jl` - Replaced inline CPTP construction in Config{Thermalize} constructor with _build_cptp_channel call
- `src/trajectories.jl` - Replaced per-operator inline CPTP construction in build_trajectoryframework with _build_cptp_channel call
- `src/jump_workers.jl` - Replaced inline CPTP construction in _finalize_kraus_step! with _build_cptp_channel call

## Decisions Made
- `_accumulate_sandwich_adj!` chosen as canonical name (shorter than `_accumulate_sandwich_adj_L!`, symmetric with `_accumulate_sandwich!`)
- `_build_cptp_channel` returns a NamedTuple with `alpha` included for callers that need it (krylov_workspace.jl does not use it; trajectories.jl framework constructor's `alpha` field is still computed independently since it is also used for coherent terms)
- `scratch.K0` field in KrausScratch struct is now dead (no longer written in `_finalize_kraus_step!`) but struct modification deferred to Phase 35 per plan guidance
- BohrDomain 2-operator sandwiches (`_accumulate_sandwich_2op!`, `_accumulate_adjoint_sandwich_2op!`) left unchanged per user decision

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Phase 34 is now complete (both plans done)
- Sandwich helpers consolidated, CPTP channel extracted -- ready for Phase 35 (struct cleanup / dead field removal)
- `scratch.K0` is a known dead field ready for removal in Phase 35
- No blockers

## Self-Check: PASSED

All 5 modified files verified present. Both task commits (a1b5363, aab8d8b) verified in git log.

---
*Phase: 34-code-deduplication*
*Completed: 2026-02-26*
