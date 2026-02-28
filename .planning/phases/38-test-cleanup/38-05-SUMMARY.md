---
phase: 38-test-cleanup
plan: 05
subsystem: testing
tags: [julia, krylov, matvec, eigsolve, crossvalidation, test-observability, threshold-rationale]

# Dependency graph
requires:
  - phase: 38-01
    provides: "Unified make_config factory and test infrastructure consolidation"
provides:
  - "@info output after every numerical @test in 3 Krylov test files"
  - "Loop-summary max_error pattern for matvec round-trip tests"
  - "Inline threshold rationale comments with FP accumulation and KrylovKit tol theory"
  - "Allocation @info showing bytes and budget for all 3 domain hot paths"
affects: []

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Loop-summary @info: track max_err across iterations, emit single @info after loop with max_error, n_samples, threshold"
    - "Threshold rationale comments: document FP accumulation theory, KrylovKit tol relationship, and safety margins"
    - "Allocation @info: report allocs_bytes and threshold for both forward and adjoint hot paths"

key-files:
  created: []
  modified:
    - test/test_krylov_matvec.jl
    - test/test_krylov_eigsolve.jl
    - test/test_krylov_crossvalidation.jl

key-decisions:
  - "Round-trip threshold 1e-12 kept (theory: O(n_jumps*DIM^2*eps) ~ 3e-13, threshold gives 3x margin)"
  - "Duality threshold 1e-11 kept (2x matvec FP accumulation, 30x margin over expected error)"
  - "Cross-validation XVAL-01/04 atol=1e-8 documented as 100x margin over KrylovKit tol=1e-10"
  - "XVAL-02 n=6 atol=1e-6 documented as 100x looser than n=4 for larger Krylov subspace error"
  - "Printf tables in crossvalidation preserved per locked research pitfall 6 decision"

patterns-established:
  - "max_err loop-summary: track max_err=0.0 before loop, max_err=max(max_err,err) inside, @info after"
  - "Allocation @info pair: one @info for forward, one for adjoint, each showing allocs_bytes and threshold"

# Metrics
duration: 12min
completed: 2026-02-28
---

# Phase 38 Plan 05: Krylov Test @info and Threshold Rationale Summary

**Added 60+ @info outputs and inline threshold rationale to all numerical assertions in 3 Krylov test files (matvec round-trips, eigsolve accuracy, cross-validation gap comparisons)**

## Performance

- **Duration:** ~12 min
- **Started:** 2026-02-28T09:03:14Z
- **Completed:** 2026-02-28T09:15:00Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments
- Added 26 @info to test_krylov_matvec.jl covering 16 round-trip loops (max_error summary pattern), 6 allocation checks (forward + adjoint x 3 domains), and 3 adjoint duality checks
- Added 18 @info to test_krylov_eigsolve.jl covering Chen channel properties (trace preservation, positivity, Euler closeness), gap accuracy, domain coverage, eigenvalue sorting, and conversion consistency
- Added 16 @info to test_krylov_crossvalidation.jl covering XVAL-01 (n=4 KMS), XVAL-04 (n=4 GNS), XVAL-03 (convergence orders), and XVAL-02 (n=6 KMS) gap comparisons
- Added inline threshold rationale comments to every numerical assertion documenting FP accumulation theory, KrylovKit tolerance relationship, channel O(delta^2) error, and system dimension scaling
- All 1057 tests pass with zero regressions

## Task Commits

Each task was committed atomically:

1. **Task 1: Add @info and threshold rationale to test_krylov_matvec.jl** - `1f09d50` (feat)
2. **Task 2: Add @info and threshold rationale to test_krylov_eigsolve.jl and test_krylov_crossvalidation.jl** - `d17e7b2` (feat)

## Files Created/Modified
- `test/test_krylov_matvec.jl` - 26 @info additions: loop-summary max_error for 16 round-trip testsets, 6 allocation @info, 3 duality @info, plus threshold rationale comments
- `test/test_krylov_eigsolve.jl` - 18 @info additions: Chen channel properties, eigsolve gap accuracy (KMS + GNS), domain coverage, eigenvalue sorting/conversion, plus threshold rationale comments
- `test/test_krylov_crossvalidation.jl` - 16 @info additions: XVAL-01 through XVAL-04 gap errors and convergence orders, plus threshold rationale comments (Printf tables preserved)

## Decisions Made
- Kept all existing thresholds unchanged after reviewing theory-based rationale (all have adequate safety margins)
- Used loop-summary pattern (track max_err across iterations, single @info after loop) to avoid 10x noise per testset
- Documented allocation tests with budget thresholds (MATVEC_ALLOC_BUDGET=0 for EnergyDomain/BohrDomain, MATVEC_ALLOC_BUDGET_NUFFT=0 for TimeDomain/TrotterDomain)
- Preserved Printf-formatted diagnostic tables in crossvalidation per locked decision from research pitfall 6

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- All 3 Krylov test files now have full @info observability and threshold documentation
- This completes the Krylov-specific test cleanup work in Phase 38

## Self-Check: PASSED

All 3 modified test files verified present on disk. Both task commits (1f09d50, d17e7b2) verified in git log. @info counts: matvec=26 (>=15), eigsolve=18 (>=6), crossval=17 (>=4). Printf table count unchanged at 7.

---
*Phase: 38-test-cleanup*
*Completed: 2026-02-28*
