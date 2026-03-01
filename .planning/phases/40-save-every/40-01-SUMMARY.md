---
phase: 40-save-every
plan: 01
subsystem: simulation
tags: [thermalization, save-every, trace-distance, density-matrix, convergence]

# Dependency graph
requires:
  - phase: 39-per-jump-precomputation
    provides: "Precomputed per-jump CPTP channels (K0s, U_residuals) in run_thermalize hot loop"
provides:
  - "save_every keyword argument for run_thermalize controlling trace distance computation frequency"
  - "recorded_steps-based time_steps construction (decoupled from num_steps)"
  - "metadata[:save_every] provenance tracking"
affects: [41-blas-threading, 42-mixing-time]

# Tech tracking
tech-stack:
  added: []
  patterns: ["save_every gating pattern: step % save_every == 0 for observation-only code"]

key-files:
  created:
    - test/test_save_every.jl
  modified:
    - src/furnace.jl
    - test/runtests.jl

key-decisions:
  - "Observation gating only: physics (coherent unitary, rho_jump, precomputed channel) runs every step; only trace_distance_h computation is gated by save_every"
  - "recorded_steps array tracks which steps were saved, replacing num_steps-based time_steps construction"
  - "Convergence cutoff checked only at save points (coarser detection is acceptable trade-off for reduced observation overhead)"

patterns-established:
  - "save_every gating: step % save_every == 0 guards observation-only code while physics runs unconditionally"

# Metrics
duration: 5min
completed: 2026-03-01
---

# Phase 40 Plan 01: Save Every Summary

**save_every keyword for run_thermalize gating trace distance computation frequency with full backward compatibility at save_every=1**

## Performance

- **Duration:** 5 min
- **Started:** 2026-03-01T10:34:36Z
- **Completed:** 2026-03-01T10:40:04Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments
- Added save_every::Int=1 keyword to run_thermalize that gates trace distance computation on step % save_every == 0
- Backward compatibility verified: save_every=1 produces bit-identical trace_distances and time_steps to pre-change behavior (same RNG seed)
- Built recorded_steps-based time_steps construction, ensuring time_steps and trace_distances always have matching lengths
- Stored save_every in metadata for provenance tracking
- Full test suite passes (1158 tests, 0 regressions)

## Task Commits

Each task was committed atomically:

1. **Task 1: Add save_every keyword to run_thermalize** - `f50dddf` (feat)
2. **Task 2: Add behavioral tests for save_every** - `796956a` (test)

## Files Created/Modified
- `src/furnace.jl` - Added save_every keyword, validation, recorded_steps tracking, observation gating, metadata storage
- `test/test_save_every.jl` - 6 behavioral test cases covering backward compatibility, array lengths, time stride, metadata, validation, and cross-domain
- `test/runtests.jl` - Added include for test_save_every.jl

## Decisions Made
- Observation gating only: physics blocks (coherent unitary, rho_jump accumulation, precomputed channel application) run every step unconditionally; only trace_distance_h computation, push, printf, and convergence check are gated
- recorded_steps array replaces num_steps-based time_steps construction -- time_steps = T.(recorded_steps .* config.delta)
- Convergence cutoff checked only at save points -- coarser detection is acceptable trade-off for reduced observation cost

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- save_every feature complete and tested, ready for Phase 41 (BLAS/omega threading)
- Phase 42 (mixing time estimation) can use save_every for efficient long-run simulations

## Self-Check: PASSED

All files and commits verified:
- FOUND: src/furnace.jl
- FOUND: test/test_save_every.jl
- FOUND: test/runtests.jl
- FOUND: f50dddf (Task 1 commit)
- FOUND: 796956a (Task 2 commit)

---
*Phase: 40-save-every*
*Completed: 2026-03-01*
