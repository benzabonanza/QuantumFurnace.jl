---
phase: 13-multi-threaded-trajectory-engine
plan: 02
subsystem: testing
tags: [threading, determinism, BLAS, performance, Xoshiro, bitwise-reproducibility]

# Dependency graph
requires:
  - phase: 13-multi-threaded-trajectory-engine
    plan: 01
    provides: "Multi-threaded trajectory engine with @sync/@spawn, BLAS control, per-trajectory Xoshiro seeding"
provides:
  - "Deterministic multi-threaded result verification (bitwise identical rho_mean for same seed+nthreads)"
  - "BLAS thread restoration test (get_num_threads before == after)"
  - "Serial-threaded agreement test (isapprox within 1e-13)"
  - "Observable path determinism test (bitwise identical measurements_mean)"
  - "Threading speedup regression test (threaded < serial for ntraj=200, dim=8)"
affects: [14-convergence-checking, 15-adaptive-batching]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Thread-count guarded testsets: if Threads.nthreads() > 1 with graceful skip"
    - "Manual serial baseline via _evolve_along_trajectory! loop for fair threading comparison"
    - "3-run averaging for performance measurements to smooth JIT/GC noise"

key-files:
  created:
    - "test/test_threading.jl"
  modified:
    - "test/runtests.jl"

key-decisions:
  - "Increased performance test workload from plan (ntraj=50, mixing_time=1.0) to (ntraj=200, mixing_time=5.0) to amortize threading overhead on dim=8 system"

patterns-established:
  - "Thread count guard pattern: wrap multi-threaded tests in if Threads.nthreads() > 1 with @info skip message and @test true placeholder"

# Metrics
duration: 10min
completed: 2026-02-16
---

# Phase 13 Plan 02: Threading Tests Summary

**5 threading testsets verifying bitwise determinism, BLAS safety, serial-threaded agreement, observable reproducibility, and 1.57x speedup with 2 threads**

## Performance

- **Duration:** 10 min
- **Started:** 2026-02-16T05:19:16Z
- **Completed:** 2026-02-16T05:29:12Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments
- Bitwise determinism verified: same seed + same nthreads produces identical rho_mean (strict ==, not isapprox)
- BLAS thread count restored to original value after run_trajectories (BLAS safety)
- Serial-threaded agreement within atol=1e-13 (validates accumulation order independence)
- Observable path determinism: measurements_mean bitwise identical across runs
- Threading speedup: 1.57x with 2 threads on 200 trajectories x 500 steps (dim=8)
- All tests skip gracefully with nthreads==1 (252 tests pass vs 257 with 2 threads)

## Task Commits

Each task was committed atomically:

1. **Task 1: Create threading correctness tests** - `57026b6` (test)
2. **Task 2: Add threading performance test** - `df99ef3` (test)

## Files Created/Modified
- `test/test_threading.jl` - 5 testsets: deterministic results, BLAS restoration, serial-threaded agreement, observable path determinism, threading speedup
- `test/runtests.jl` - Added include("test_threading.jl") after test_workspace_independence.jl

## Decisions Made
- Increased performance test workload from plan specification (ntraj=50, mixing_time=1.0) to (ntraj=200, mixing_time=5.0) because dim=8 per-step cost is too low for threading overhead amortization at smaller scales

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Adjusted performance test parameters for reliable speedup measurement**
- **Found during:** Task 2 (threading performance test)
- **Issue:** Plan specified ntraj=50 with mixing_time=1.0 (100 steps per trajectory), but dim=8 steps are so fast (~microseconds) that @sync/@spawn threading overhead dominates, producing 0.24x "speedup" (slower than serial)
- **Fix:** Increased to ntraj=200 with mixing_time=5.0 (500 steps per trajectory), giving ~100x more total work to amortize threading overhead. This reliably produces 1.4-1.6x speedup with 2 threads.
- **Files modified:** test/test_threading.jl
- **Verification:** Performance test passes: speedup=1.57x logged, @test t_threaded < t_serial succeeds
- **Committed in:** df99ef3 (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (1 bug - test parameters)
**Impact on plan:** Essential fix -- without larger workload the performance test would always fail on dim=8 systems. No scope creep; test still validates the same property (threaded faster than serial).

## Issues Encountered
None beyond the performance test parameter adjustment documented above.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Phase 13 fully complete: multi-threaded trajectory engine (Plan 01) + comprehensive threading tests (Plan 02)
- All 257 tests pass with 2 threads; 252 pass with 1 thread (threading tests skip gracefully)
- Ready for Phase 14 (convergence checking) which builds on the trajectory infrastructure

## Self-Check: PASSED

All files and commits verified:
- test/test_threading.jl: FOUND
- test/runtests.jl: FOUND
- 13-02-SUMMARY.md: FOUND
- Commit 57026b6: FOUND
- Commit df99ef3: FOUND

---
*Phase: 13-multi-threaded-trajectory-engine*
*Completed: 2026-02-16*
