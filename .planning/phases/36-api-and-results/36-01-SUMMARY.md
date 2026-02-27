---
phase: 36-api-and-results
plan: 01
subsystem: api
tags: [julia, structs, bson, serialization, results]

# Dependency graph
requires:
  - phase: 35-workspace-consolidation
    provides: "Unified Workspace{S,D,C,T,SC}, ConvergenceData, Config{S,D,C,T}"
provides:
  - "AbstractResults supertype for dispatch"
  - "LindbladResults, ThermalizeResults, KrylovSpectrumResults, TrajectoryResults structs"
  - "save_result/load_result with Dict-based BSON and auto type detection"
  - "_capture_metadata without julia_version (locked decision)"
  - "_generate_result_filename for auto-naming"
  - "Companion .txt with per-type result summaries"
affects: [36-02-PLAN, 36-03-PLAN, 36-04-PLAN, 37-removal]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Dict-based BSON serialization with :result_type tag for type dispatch on load"
    - "_result_to_dict / _dict_to_*_results pattern for each Result type"
    - "Companion .txt with type-specific sections"

key-files:
  created: []
  modified:
    - "src/structs.jl"
    - "src/results.jl"
    - "src/QuantumFurnace.jl"
    - "test/test_results.jl"

key-decisions:
  - "config_kind tag changed from 'liouv' to 'lindbladian' for new saves; _reconstruct_config accepts both for backward compat"
  - "_result_to_dict uses multiple dispatch (one method per concrete type) rather than if/elseif"
  - "_trajectory_to_dict_new suffix to avoid name clash with existing _trajectory_to_dict for TrajectoryResult"

patterns-established:
  - "Result type tag pattern: each Result stores :result_type in BSON Dict, load_result dispatches on it"
  - "Companion .txt pattern: type-specific sections via if/elseif in _write_result_companion_txt"

# Metrics
duration: 9min
completed: 2026-02-27
---

# Phase 36 Plan 01: Result Structs and Serialization Summary

**AbstractResults type hierarchy with 4 concrete Result structs and Dict-based BSON save_result/load_result with companion .txt**

## Performance

- **Duration:** 9 min
- **Started:** 2026-02-27T11:57:41Z
- **Completed:** 2026-02-27T12:06:13Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments
- Defined AbstractResults supertype and 4 concrete Result structs (LindbladResults, ThermalizeResults, KrylovSpectrumResults, TrajectoryResults) with correct fields matching locked decisions
- Implemented save_result/load_result with Dict-based BSON serialization and auto-detection of result type via :result_type tag
- Removed :julia_version from _capture_metadata per locked decision (git hash, timestamp, wall_time_seconds, n_threads only)
- Updated _config_to_dict to store "lindbladian" instead of "liouv" while maintaining backward compatibility in _reconstruct_config
- Added companion .txt generation with type-specific result summaries
- All 1198 existing tests pass with zero regressions

## Task Commits

Each task was committed atomically:

1. **Task 1: Define AbstractResults and 4 Result structs in structs.jl** - `8f9b890` (feat)
2. **Task 2: Update metadata capture, implement save_result/load_result, and companion .txt** - `3ae3d22` (feat)

## Files Created/Modified
- `src/structs.jl` - Added AbstractResults abstract type and 4 concrete Result struct definitions after ConvergenceData
- `src/results.jl` - Removed julia_version from _capture_metadata; added _result_type_tag, _result_to_dict, _dict_to_*_results, save_result, load_result, _write_result_companion_txt, _generate_result_filename
- `src/QuantumFurnace.jl` - Added exports for AbstractResults, LindbladResults, ThermalizeResults, KrylovSpectrumResults, TrajectoryResults, save_result, load_result
- `test/test_results.jl` - Updated metadata auto-capture tests to reflect julia_version removal

## Decisions Made
- Changed config_kind tag from "liouv" to "lindbladian" for new saves, with backward-compatible "liouv" acceptance in _reconstruct_config
- Used multiple dispatch for _result_to_dict (one method per concrete Result type) rather than a single function with if/elseif
- Named the TrajectoryResults serializer _trajectory_to_dict_new to avoid clash with existing _trajectory_to_dict for the old TrajectoryResult struct

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Updated existing tests for julia_version removal**
- **Found during:** Task 2 (metadata capture update)
- **Issue:** test_results.jl tests 9 and 10 asserted that _capture_metadata includes :julia_version, which we intentionally removed per locked decision
- **Fix:** Changed tests to assert !haskey(meta, :julia_version) instead of haskey
- **Files modified:** test/test_results.jl
- **Verification:** Full test suite passes (1198/1198)
- **Committed in:** 3ae3d22 (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (1 bug - test alignment with intentional behavior change)
**Impact on plan:** Essential fix for test correctness. No scope creep.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Result type hierarchy and serialization layer complete
- Plans 02-03 can now construct and return typed Results from run_* entry points
- Plan 04 can write comprehensive tests for save_result/load_result round-trips
- All existing code unchanged (old ExperimentResult, TrajectoryResult, etc. still functional)

## Self-Check: PASSED

All files verified present, all commits verified in git log.

---
*Phase: 36-api-and-results*
*Completed: 2026-02-27*
