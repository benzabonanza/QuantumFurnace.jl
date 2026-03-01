---
phase: 41-threading
plan: 02
subsystem: threading
tags: [omega-loop, threading, density-matrix, blas, per-task-scratch, spawn]

# Dependency graph
requires:
  - phase: 41-threading
    plan: 01
    provides: "BLAS try/finally wrapping in run_thermalize hot loop"
  - phase: 39-precomputation
    provides: "per-jump precomputed channels (_accumulate_rho_jump!, _precompute_R)"
provides:
  - "Threaded _accumulate_rho_jump! for Energy, Time/Trotter, and Bohr domains (THREAD-01, THREAD-04)"
  - "Threaded _precompute_R for Energy, Time/Trotter, and Bohr domains (THREAD-02)"
  - "OMEGA_THREAD_THRESHOLD constant and _partition_range utility"
  - "Serial-threaded agreement tests for omega-loop parallelism"
affects: [mixing-time-estimation, performance-benchmarking]

# Tech tracking
tech-stack:
  added: []
  patterns: ["Per-task ThermalizeScratch for omega-loop parallelism with BLAS=1 during Julia-level threading", "Threading gate pattern: nthreads > 1 && work_count >= OMEGA_THREAD_THRESHOLD"]

key-files:
  created: []
  modified:
    - src/jump_workers.jl
    - src/trajectories.jl
    - test/test_threading.jl

key-decisions:
  - "OMEGA_THREAD_THRESHOLD = 50: balances task spawn overhead vs parallelism benefit"
  - "Per-task ThermalizeScratch allocation per @sync block (not pooled) for simplicity and correctness"
  - "Half-grid pre-filtering for Hermitian jumps ensures balanced partition across tasks"
  - "_partition_range utility duplicated from _partition_trajectories pattern for locality"

patterns-established:
  - "Omega-loop threading: gate -> partition -> per-task scratch -> BLAS=1 -> @sync/@spawn -> reduce"
  - "BLAS=1 during Julia-level parallelism (mutual exclusion with BLAS multi-threading)"

# Metrics
duration: 21min
completed: 2026-03-01
---

# Phase 41 Plan 02: Omega-Loop Threading Summary

**Omega-loop parallelism for _accumulate_rho_jump! and _precompute_R across all 4 domains with per-task ThermalizeScratch isolation and BLAS=1 mutual exclusion**

## Performance

- **Duration:** 21 min
- **Started:** 2026-03-01T11:54:06Z
- **Completed:** 2026-03-01T12:14:48Z
- **Tasks:** 3
- **Files modified:** 3

## Accomplishments
- _accumulate_rho_jump! for Energy, Time/Trotter, and Bohr domains dispatches to threaded variants when nthreads > 1 and frequency count >= 50
- _precompute_R for all 3 domain families (Energy, Time/Trotter, Bohr) has threaded variants
- Per-task ThermalizeScratch provides complete data race isolation (rho_jump, jump_oft, sandwich_tmp, R, LdagL all private per task)
- BLAS set to 1 during all omega-loop parallelism, restored via try/finally
- Serial fallback automatic when nthreads==1 or work count below threshold
- Full test suite: 1181 tests pass on -t 4 (11 new), 1168 on -t 1 (serial fallback)

## Task Commits

Each task was committed atomically:

1. **Task 1: Add threaded _accumulate_rho_jump! for Energy/Time/Trotter domains** - `d0b54e1` (feat)
2. **Task 2: Add BohrDomain threading and threaded _precompute_R for all domains** - `a51f9ce` (feat)
3. **Task 3: Add omega-loop serial-threaded agreement tests** - `4eed834` (test)

## Files Created/Modified
- `src/jump_workers.jl` - Added OMEGA_THREAD_THRESHOLD, _partition_range, threaded _accumulate_rho_jump! variants for all 4 domains with per-task scratch and BLAS=1
- `src/trajectories.jl` - Added threaded _precompute_R variants for Energy, Time/Trotter, and Bohr domains with per-task scratch and BLAS=1
- `test/test_threading.jl` - Added 3 omega-loop threading test blocks: determinism (all domains), serial vs threaded agreement, BLAS restoration

## Decisions Made
- OMEGA_THREAD_THRESHOLD set to 50 -- below this, serial execution is faster due to task spawn overhead
- Per-task ThermalizeScratch allocated fresh per @sync block rather than pooled -- simplicity wins at this allocation frequency (once per step, not per-frequency)
- Half-grid pre-filtering for Hermitian jumps ensures balanced partition across tasks (without it, some tasks would get only `continue` iterations)
- _partition_range duplicated as a utility function in jump_workers.jl rather than importing from trajectories.jl -- both files are in the same module so it works, but having the utility local to the file that defines the threading pattern improves readability

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
- Git index corruption during Task 2 commit (resolved by regenerating index with `rm .git/index && git reset`)

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- All threading optimizations for Phase 41 complete (BLAS threading + omega-loop parallelism)
- Ready for Phase 42 (mixing time estimation) or further phases in v2.1 milestone
- For large systems (N_w > 50), omega-loop parallelism provides complementary speedup to BLAS threading

## Self-Check: PASSED

All files exist, all commits verified.

---
*Phase: 41-threading*
*Completed: 2026-03-01*
