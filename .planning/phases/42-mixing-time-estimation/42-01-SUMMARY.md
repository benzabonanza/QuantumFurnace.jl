---
phase: 42-mixing-time-estimation
plan: 01
subsystem: simulation
tags: [lsqfit, exponential-decay, mixing-time, spectral-gap, curve-fitting, post-processing]

# Dependency graph
requires:
  - phase: 41-threading
    provides: "Completed BLAS + omega-loop threading for DM path"
provides:
  - "LsqFit.jl as active runtime dependency"
  - "fit_exponential_decay and FitResult exported from QuantumFurnace"
  - "estimate_mixing_time post-processing API with MixingTimeEstimate struct"
  - "Quality gate warnings for fit reliability"
  - "Extrapolation formula with guards for mixing time prediction"
affects: [42-02, testing, mixing-time]

# Tech tracking
tech-stack:
  added: [LsqFit.jl v0.15]
  patterns: [post-processing-wrapper, quality-gate-warnings, extrapolation-with-guards]

key-files:
  created:
    - src/fitting.jl
    - src/mixing.jl
  modified:
    - Project.toml
    - Manifest.toml
    - src/QuantumFurnace.jl

key-decisions:
  - "Pass skip_initial through to fit_exponential_decay (option b from Research pitfall 3)"
  - "MixingTimeEstimate is separate struct, not modification of ThermalizeResults"
  - "Primary mixing_time uses extrapolated value when extrapolate=true, actual crossing when target_epsilon provided, total sim time otherwise"

patterns-established:
  - "Post-processing wrapper: thin API over fit_exponential_decay with quality gates"
  - "Quality gates via @warn: non-blocking warnings for fit reliability issues"
  - "Guard-heavy extrapolation: multiple early returns for edge cases"

# Metrics
duration: 3min
completed: 2026-03-01
---

# Phase 42 Plan 01: Mixing Time Estimation Summary

**LsqFit.jl promoted from staging with estimate_mixing_time post-processing API, quality gates, and extrapolation formula**

## Performance

- **Duration:** 3 min
- **Started:** 2026-03-01T12:43:56Z
- **Completed:** 2026-03-01T12:46:57Z
- **Tasks:** 2
- **Files modified:** 5

## Accomplishments
- Promoted LsqFit.jl and fitting.jl from staging to active runtime dependency with full exports
- Implemented MixingTimeEstimate struct with 12 fields covering fit params, mixing times, and quality metrics
- Built estimate_mixing_time(::ThermalizeResults) API with skip_initial, target_epsilon, extrapolate, and level keywords
- Added quality gate warnings via @warn for R-squared, offset, convergence, and SE thresholds
- Implemented extrapolation formula with 5 guard clauses for edge cases

## Task Commits

Each task was committed atomically:

1. **Task 1: Promote LsqFit.jl and fitting.jl to active source** - `4f015cb` (feat)
2. **Task 2: Implement MixingTimeEstimate and estimate_mixing_time** - `0374cf5` (feat)

## Files Created/Modified
- `Project.toml` - Added LsqFit.jl to [deps] and [compat]
- `Manifest.toml` - Resolved LsqFit.jl and transitive dependencies
- `src/QuantumFurnace.jl` - Added using LsqFit, include fitting.jl/mixing.jl, exports for all new symbols
- `src/fitting.jl` - Promoted from staging: FitResult struct, fit_exponential_decay, log-linear initial guess
- `src/mixing.jl` - New: MixingTimeEstimate struct, estimate_mixing_time, quality gates, extrapolation helpers

## Decisions Made
- Pass skip_initial through to fit_exponential_decay without pre-processing (Research pitfall 3, option b) -- avoids double-truncation
- MixingTimeEstimate is a separate immutable struct, not a modification of ThermalizeResults -- preserves BSON compatibility
- Primary mixing_time field logic: extrapolated when extrapolate=true, actual crossing when target_epsilon provided, total simulation time otherwise

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Pkg.resolve() required before Pkg.instantiate()**
- **Found during:** Task 1 (LsqFit installation)
- **Issue:** `Pkg.instantiate()` failed because LsqFit was in Project.toml but not in Manifest.toml; Pkg required explicit resolve first
- **Fix:** Ran `Pkg.resolve()` before `Pkg.instantiate()` to populate the manifest
- **Files modified:** Manifest.toml
- **Verification:** Module loads with `using QuantumFurnace` and all exports work
- **Committed in:** 4f015cb (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (1 blocking)
**Impact on plan:** Minor package management step. No scope creep.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Core implementation complete: estimate_mixing_time API is exported and compiles
- Ready for Phase 42 Plan 02: testing the mixing time estimation
- All new symbols verified via Julia REPL: methods, fieldnames, compilation

## Self-Check: PASSED

All files verified present, all commit hashes found in git log.

---
*Phase: 42-mixing-time-estimation*
*Completed: 2026-03-01*
