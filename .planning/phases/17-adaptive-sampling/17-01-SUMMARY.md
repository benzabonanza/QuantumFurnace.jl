---
phase: 17-adaptive-sampling
plan: 01
subsystem: convergence
tags: [adaptive-sampling, convergence, trace-distance, windowed-average, batch-loop]

# Dependency graph
requires:
  - phase: 16-convergence-tracking
    provides: "ConvergenceData struct, run_trajectories_convergence batch loop pattern, observable builders, Dict serialization"
provides:
  - "Extended ConvergenceData with 4 adaptive diagnostic fields (converged, final_relative_change, consecutive_stable_batches, total_batches)"
  - "run_trajectories_adaptive function with automatic stopping based on windowed trace distance convergence"
  - "_windowed_relative_change helper for convergence detection"
  - "Backward-compatible 6-argument ConvergenceData constructor"
  - "Updated Dict serialization for 10-field ConvergenceData (forward and backward compatible)"
affects: [17-02-adaptive-sampling-tests, 18-experiments]

# Tech tracking
tech-stack:
  added: []
  patterns: ["push!-based dynamic storage for unknown-length batch loops", "windowed relative change for convergence detection"]

key-files:
  created: []
  modified:
    - src/convergence.jl
    - src/results.jl
    - src/QuantumFurnace.jl

key-decisions:
  - "Broad-typed outer constructor to accept BSON-deserialized Vector{Any} (fixes pre-existing BSON round-trip edge case)"
  - "Ceiling division (cld) for max_batches accepts slight n_max overshoot to maintain fixed batch size invariant"
  - "effective_min = max(min_batches, 2*window_size) silently clamps to ensure enough data for windowed comparison"

patterns-established:
  - "Push-based dynamic storage: use Float64[], Vector{Float64}[], etc. for unknown-length batch loops, convert at end with reduce(hcat, ...)"
  - "Windowed convergence detection: compare mean of last K vs previous K entries with eps-safe denominator"

# Metrics
duration: 9min
completed: 2026-02-16
---

# Phase 17 Plan 01: Adaptive Sampling Implementation Summary

**Adaptive trajectory batch loop with windowed trace distance convergence detection, extended ConvergenceData diagnostics, and backward-compatible serialization**

## Performance

- **Duration:** 9 min
- **Started:** 2026-02-16T10:56:01Z
- **Completed:** 2026-02-16T11:05:35Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments
- Extended ConvergenceData struct with 4 adaptive diagnostic fields while maintaining full backward compatibility with Phase 16 code (470 tests pass unchanged)
- Implemented run_trajectories_adaptive with configurable windowed convergence detection (threshold, patience, min_batches, window_size)
- Updated Dict serialization to handle both old 6-field and new 10-field ConvergenceData formats
- Added _windowed_relative_change helper with safe handling of insufficient data (returns Inf)

## Task Commits

Each task was committed atomically:

1. **Task 1: Extend ConvergenceData, implement _windowed_relative_change and run_trajectories_adaptive** - `bbdd32c` (feat)
2. **Task 2: Update serialization and module exports for adaptive sampling** - `3df92bd` (feat)

## Files Created/Modified
- `src/convergence.jl` - Extended ConvergenceData (10 fields), 6-arg backward-compat constructor, _windowed_relative_change helper, run_trajectories_adaptive function
- `src/results.jl` - Updated _convergence_to_dict (10 fields) and _dict_to_convergence (backward-compatible defaults)
- `src/QuantumFurnace.jl` - Added run_trajectories_adaptive to export list

## Decisions Made
- Used broad (untyped) parameter signatures in the 6-argument outer constructor to accept BSON-deserialized data (Vector{Any} for strings) -- fixes a latent type mismatch that was previously hidden by the auto-generated inner constructor
- Ceiling division (cld) for computing max_batches, accepting slight n_max overshoot rather than truncating the last batch (maintains fixed batch size invariant per locked decision)
- effective_min silently clamps min_batches to 2*window_size to ensure sufficient data for windowed comparison

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed BSON round-trip type mismatch in ConvergenceData outer constructor**
- **Found during:** Task 1 (extending ConvergenceData)
- **Issue:** BSON deserializes `Vector{String}` as `Vector{Any}`. The typed 6-argument outer constructor (`observable_names::Vector{String}`) rejected BSON-loaded data, causing the existing BSON round-trip test to error. Previously, the auto-generated inner constructor (which accepts Any) handled this implicitly.
- **Fix:** Changed outer constructor to use broad (untyped) parameter signatures, matching Julia's standard struct constructor pattern for interoperability with BSON.
- **Files modified:** src/convergence.jl
- **Verification:** BSON round-trip test passes, all 470 tests pass
- **Committed in:** bbdd32c (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (1 bug)
**Impact on plan:** Essential fix for backward compatibility. No scope creep.

## Issues Encountered
- Git index corruption after Task 1 commit required index rebuild (`rm .git/index && git reset`). No data lost; Task 2 committed successfully after rebuild.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Adaptive sampling function is implemented and exported, ready for comprehensive testing in Phase 17 Plan 02
- All 470 existing tests pass, confirming backward compatibility
- ConvergenceData serialization handles both old and new formats

## Self-Check: PASSED

All created/modified files verified present. All commit hashes verified in git log.

---
*Phase: 17-adaptive-sampling*
*Completed: 2026-02-16*
