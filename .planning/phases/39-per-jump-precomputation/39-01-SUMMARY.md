---
phase: 39-per-jump-precomputation
plan: 01
subsystem: thermalize
tags: [cptp-channel, precomputation, eigendecomposition, kraus-operators, bohr-domain]

# Dependency graph
requires:
  - phase: 38-trajectory-workspace
    provides: "_precompute_R for Energy/Time/Trotter, _build_cptp_channel, ThermalizeScratch"
provides:
  - "_precompute_R for BohrDomain (new dispatch method in trajectories.jl)"
  - "_precompute_per_jump_channels orchestrator (furnace_utensils.jl)"
  - "_accumulate_rho_jump! for Energy, Time/Trotter, Bohr domains (jump_workers.jl)"
  - "_apply_precomputed_channel! (furnace_utensils.jl)"
affects: [39-02-integration, run_thermalize, hot-loop-optimization]

# Tech tracking
tech-stack:
  added: []
  patterns: ["per-jump CPTP channel precomputation", "rho_jump-only omega loop extraction"]

key-files:
  created: []
  modified:
    - src/trajectories.jl
    - src/furnace_utensils.jl
    - src/jump_workers.jl

key-decisions:
  - "BohrDomain _precompute_R uses precomputed bohr_is/bohr_js from precomputed_data when available, with fallback to hamiltonian.bohr_dict"
  - "_precompute_per_jump_channels returns (K0s, U_residuals) only (no Rs) since DM path does not need stored Rs"
  - "_accumulate_rho_jump! takes jump_weight_scaling as keyword arg matching existing _jump_contribution! convention"

patterns-established:
  - "Per-jump R precomputation reuses existing _precompute_R + _build_cptp_channel pipeline"
  - "Omega-loop functions split into R-accumulation (_precompute_R) and rho_jump-only (_accumulate_rho_jump!)"

# Metrics
duration: 3min
completed: 2026-03-01
---

# Phase 39 Plan 01: Per-Jump Precomputation Infrastructure Summary

**Per-jump CPTP channel precomputation functions: BohrDomain _precompute_R, channel orchestrator, 3-domain rho_jump accumulators, and precomputed channel applicator**

## Performance

- **Duration:** 3 min
- **Started:** 2026-03-01T09:29:21Z
- **Completed:** 2026-03-01T09:32:29Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments
- Added `_precompute_R` for BohrDomain completing the 3-domain set (Energy, Time/Trotter, Bohr) needed for per-jump R precomputation
- Created `_precompute_per_jump_channels` orchestrator that calls `_precompute_R` + `_build_cptp_channel` per jump, returning `(; K0s, U_residuals)` vectors
- Implemented `_accumulate_rho_jump!` for all 3 domains extracting only rho_jump accumulation from `_jump_contribution!` (no R/LdagL/eigen computation)
- Added `_apply_precomputed_channel!` applying pre-stored K0/U_residual sandwich without eigendecomposition

## Task Commits

Each task was committed atomically:

1. **Task 1: BohrDomain _precompute_R and _precompute_per_jump_channels** - `5da6732` (feat)
2. **Task 2: _accumulate_rho_jump! and _apply_precomputed_channel!** - `bcc9793` (feat)

## Files Created/Modified
- `src/trajectories.jl` - Added `_precompute_R` for `Config{Thermalize, BohrDomain}` with bohr_is/bohr_js optimization
- `src/furnace_utensils.jl` - Added `_precompute_per_jump_channels` orchestrator and `_apply_precomputed_channel!`
- `src/jump_workers.jl` - Added 3 `_accumulate_rho_jump!` methods (EnergyDomain, TimeDomain/TrotterDomain, BohrDomain)

## Decisions Made
- BohrDomain `_precompute_R` uses precomputed `bohr_is`/`bohr_js` from `precomputed_data` when available, falling back to `hamiltonian.bohr_dict` -- matches the existing `_jump_contribution!` Bohr Thermalize pattern
- `_precompute_per_jump_channels` does NOT store per-jump Rs (only K0s and U_residuals) since the DM hot loop only needs the channel matrices, not the raw R
- `_accumulate_rho_jump!` uses `jump_weight_scaling` keyword parameter matching the existing `_jump_contribution!` convention for consistent scaling

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- All four function groups ready for integration into `run_thermalize` in Plan 39-02
- `_precompute_per_jump_channels` is the entry point for Plan 39-02 to call before the hot loop
- `_accumulate_rho_jump!` + `_apply_precomputed_channel!` replace `_jump_contribution!` + `_finalize_kraus_step!` in the hot loop
- Project compiles cleanly with all new functions accessible

## Self-Check: PASSED

All files and commits verified.

---
*Phase: 39-per-jump-precomputation*
*Completed: 2026-03-01*
