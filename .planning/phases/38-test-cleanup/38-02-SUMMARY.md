---
phase: 38-test-cleanup
plan: 02
subsystem: testing
tags: [julia, test-output, info-logging, threshold-rationale, numerical-assertions]

# Dependency graph
requires:
  - phase: 38-01
    provides: "Unified make_config factory, N3_* globals, ALL_DOMAINS constant"
provides:
  - "@info after every numerical @test in 8 shorter test files"
  - "Loop-summary @info pattern (max_error tracking) for CPTP and probability conservation tests"
  - "Inline threshold rationale comments on all numerical assertions"
affects: [38-03, 38-04, 38-05]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "@info AFTER @test with keyword args: @info 'label' value=computed threshold=threshold"
    - "Loop-summary pattern: track max_err across iterations, emit ONE @info after loop"
    - "Threshold rationale as inline comments: # algebraic identity, FP accumulation scaling, statistical noise"

key-files:
  created: []
  modified:
    - test/test_compilation.jl
    - test/test_cptp.jl
    - test/test_dm_detailed_balance.jl
    - test/test_trajectory_fixes.jl
    - test/test_regression.jl
    - test/test_allocation.jl
    - test/test_workspace_independence.jl
    - test/test_threading.jl

key-decisions:
  - "CPTP threshold 1e-10 kept: algebraic identity error scales as DIM^2 * eps ~ 3e-13, giving ~300x margin"
  - "DMTST-01 @info moved from before @test to after @test (locked decision compliance)"
  - "TFIX-04 PSD eigenvalue threshold -1e-14: allows FP rounding in Hermitian eigvals while catching real PSD failures"
  - "Allocation @info includes both allocs_bytes and threshold for direct comparison in test output"
  - "Threading speedup @info kept alongside existing performance @info for double-check visibility"

patterns-established:
  - "@info 'label' key=value: structured logging after every numerical @test"
  - "max_err loop-summary: compute aggregate stats in loop, emit single @info after"
  - "Inline threshold rationale: every atol/rtol/< threshold has a # comment explaining derivation"

# Metrics
duration: 10min
completed: 2026-02-28
---

# Phase 38 Plan 02: @info Output and Threshold Rationale for 8 Test Files Summary

**Added @info after every numerical @test in 8 shorter test files with loop-summary patterns for CPTP/probability checks and inline threshold rationale for all numerical assertions**

## Performance

- **Duration:** ~10 min
- **Started:** 2026-02-28T09:03:17Z
- **Completed:** 2026-02-28T09:13:00Z
- **Tasks:** 2
- **Files modified:** 8

## Accomplishments
- Added 38 @info statements across 8 test files, covering all numerical @test assertions
- Applied loop-summary @info pattern for CPTP completeness (test_cptp, test_trajectory_fixes) and probability conservation (test_trajectory_fixes) to avoid per-iteration noise
- Documented threshold rationale for every numerical comparison: algebraic identity scaling (DIM^2 * eps), statistical noise (1/sqrt(N)), FP accumulation, PSD rounding, allocation budgets
- Moved DMTST-01 @info from before @test to after @test (locked decision compliance)
- Zero structural/type checks received @info (isa, ===, haskey, boolean equality all skipped correctly)

## Task Commits

Each task was committed atomically:

1. **Task 1: Add @info and threshold rationale to test_compilation, test_cptp, test_dm_detailed_balance, test_trajectory_fixes** - `4b25b1a` (feat)
2. **Task 2: Add @info and threshold rationale to test_regression, test_allocation, test_workspace_independence, test_threading** - `1323fcc` (feat)

## Files Created/Modified
- `test/test_compilation.jl` - @info for Gibbs trace normalization (1 @info)
- `test/test_cptp.jl` - Loop-summary @info for CPTP completeness across 3 domains (3 @info), threshold rationale in header comment
- `test/test_dm_detailed_balance.jl` - Moved DMTST-01 @info after @test, added per-comparison @info for DMTST-02 hierarchy (5 @info total)
- `test/test_trajectory_fixes.jl` - @info for TFIX-02..05 with summary stats for loop-based checks (5 @info)
- `test/test_regression.jl` - @info for DM and trajectory regression with max_element_error (4 @info)
- `test/test_allocation.jl` - @info for all 8 allocation checks with allocs_bytes and threshold (8 @info)
- `test/test_workspace_independence.jl` - @info for norm preservation, deterministic replay, trajectory trace (5 @info)
- `test/test_threading.jl` - @info for serial-threaded agreement and speedup timing (7 @info including existing skip messages)

## Decisions Made
- Kept CPTP threshold at 1e-10 (theory: DIM^2 * eps ~ 3e-13, so 1e-10 gives ~300x margin for FP accumulation)
- DMTST-01 @info was originally placed before @test; moved it after per locked decision
- PSD guard threshold of -1e-14 documented: allows FP rounding in Hermitian eigvals decomposition
- Allocation @info includes both measured bytes and computed threshold for easy comparison
- Skip @info for all negative checks (!isapprox for independence verification) -- these are structural assertions

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
- Full `Pkg.test()` was killed (OOM in sandbox environment) but individual file tests all pass
- Threading tests skip when nthreads=1 (expected behavior in single-threaded sandbox)

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- All 8 shorter test files now have self-documenting @info output
- Loop-summary pattern established for use in larger test files (plans 03-05)
- Threshold rationale convention consistent across files for future reference

## Self-Check: PASSED

All 8 modified files verified present on disk. Both task commits (4b25b1a, 1323fcc) verified in git log.

---
*Phase: 38-test-cleanup*
*Completed: 2026-02-28*
