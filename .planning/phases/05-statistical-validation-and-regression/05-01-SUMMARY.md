---
phase: 05-statistical-validation-and-regression
plan: 01
subsystem: testing
tags: [monte-carlo, convergence, trajectory, statistics, trace-distance]

# Dependency graph
requires:
  - phase: 04-trajectory-cross-validation
    provides: trajectory framework, step_along_trajectory!, single-step cross-validation
  - phase: quick-7
    provides: per-operator Lie-Trotter splitting for trajectory steps
provides:
  - TVAL-05 gated convergence test verifying 1/sqrt(N_traj) scaling
  - batch-averaged convergence testing methodology for trajectory Monte Carlo
affects: [05-02, regression-tests]

# Tech tracking
tech-stack:
  added: [Statistics stdlib]
  patterns: [batch-averaged convergence testing, high-N trajectory reference]

key-files:
  created:
    - test/trajectory_validation/run_convergence_tests.jl

key-decisions:
  - "High-N trajectory reference (500k) instead of Liouvillian DM reference to avoid Lie-Trotter splitting bias"
  - "Batch-averaged errors (10 batches per N_traj) for robust ratio estimation vs single-realization"
  - "N_traj points [200, 800, 3200, 12800] with factor-of-4 steps (smaller N than planned due to batch approach)"

patterns-established:
  - "Convergence testing: compare against same-method high-N reference, not analytical reference, when systematic bias exists"
  - "Batch averaging: use independent seeds per batch for noise reduction in convergence ratio tests"

# Metrics
duration: 14min
completed: 2026-02-14
---

# Phase 5 Plan 1: Gated 1/sqrt(N) Convergence Test Summary

**Batch-averaged 1/sqrt(N_traj) convergence test using high-N trajectory reference for EnergyDomain and TrotterDomain (with_coherent=true)**

## Performance

- **Duration:** 14 min
- **Started:** 2026-02-14T16:26:05Z
- **Completed:** 2026-02-14T16:40:50Z
- **Tasks:** 1
- **Files created:** 1

## Accomplishments
- TVAL-05 requirement satisfied: 1/sqrt(N) scaling confirmed for both EnergyDomain and TrotterDomain with coherent
- All 6 ratio checks (3 per domain) within [1.5, 2.5], typical values ~1.9-2.0
- Test gated behind QUANTUMFURNACE_FULL_TESTS=true, does NOT run during Pkg.test()
- Total gated test runtime ~15 seconds

## Task Commits

Each task was committed atomically:

1. **Task 1: Create gated 1/sqrt(N) convergence test** - `7e71503` (feat)

**Plan metadata:** pending (docs: complete plan)

## Files Created/Modified
- `test/trajectory_validation/run_convergence_tests.jl` - Standalone gated convergence test verifying TVAL-05

## Decisions Made
- **High-N trajectory reference instead of DM reference:** The per-operator Lie-Trotter splitting (quick-7) introduces O(delta_eff) = O(delta * n_jumps) systematic error when comparing trajectory vs Liouvillian DM evolution. At delta=0.1 with 9 jumps (3-qubit system), this gives ~0.48 trace distance, completely dominating 1/sqrt(N) statistical error. Using 500k trajectory reference eliminates this bias.
- **Batch-averaged approach:** Single-realization convergence ratios are noisy (individual ratios can fall outside [1.5, 2.5] due to sampling variance). Averaging 10 independent batches per N_traj point gives robust mean errors with ratios consistently near 2.0.
- **N_traj progression [200, 800, 3200, 12800]:** Smaller base N than plan's [1000, 4000, 16000, 64000] because batch overhead (10 batches) compensates, and the reference noise floor at N_ref=500k is ~0.001 (negligible compared to test errors ~0.009-0.064).

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Plan assumed O(delta^2) systematic error for DM comparison**
- **Found during:** Task 1 (initial test implementation)
- **Issue:** Plan assumed trajectory vs DM trace distance is O(delta^2) systematic + 1/sqrt(N) statistical. After quick-7 per-operator Lie-Trotter splitting, the systematic error is O(delta * n_jumps) = O(delta_eff), which at delta=0.1 with 9 jumps gives 0.48 trace distance. This completely masks 1/sqrt(N) convergence.
- **Fix:** Replaced DM reference with high-N trajectory reference (N=500,000 with independent seed). The trajectory reference approximates the true mean of the trajectory channel (which includes the Lie-Trotter splitting), so the only remaining error is 1/sqrt(N) statistical noise.
- **Files modified:** test/trajectory_validation/run_convergence_tests.jl
- **Verification:** All 6 ratio checks pass, ratios consistently ~1.9-2.0
- **Committed in:** 7e71503

**2. [Rule 1 - Bug] Single-realization convergence ratios too noisy**
- **Found during:** Task 1 (testing single-realization approach)
- **Issue:** With a single batch per N_traj, individual trace distance realizations have high variance. Ratios can fall outside [1.5, 2.5] (observed 1.15 for TrotterDomain first ratio).
- **Fix:** Implemented batch-averaged methodology: 10 independent batches per N_traj point, each with unique seed, averaged to get robust mean error estimates.
- **Files modified:** test/trajectory_validation/run_convergence_tests.jl
- **Verification:** Batch-averaged ratios consistently 1.84-2.03 across both domains
- **Committed in:** 7e71503

---

**Total deviations:** 2 auto-fixed (2 bugs in plan assumptions)
**Impact on plan:** Both auto-fixes necessary due to interaction with quick-7 Lie-Trotter refactor. The test correctly validates 1/sqrt(N) convergence with a more robust methodology than originally planned.

## Issues Encountered
- Discovered that per-operator Lie-Trotter splitting (quick-7) introduces O(delta_eff) systematic bias between trajectory and Liouvillian DM evolution. This is expected mathematically (first-order product formula) but was not accounted for in the plan's delta=0.1 assumption. The fix (high-N trajectory reference) is more principled anyway since it tests the actual trajectory convergence without conflating splitting error.
- Phase 4 trajectory validation tests (run_trajectory_validation.jl) now fail for delta=0.2 due to delta_eff assertion (0.2*9=1.8 > 1.0). This is a pre-existing issue from quick-7 that should be addressed separately.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- TVAL-05 convergence test complete and passing
- Ready for plan 05-02 (next test in statistical validation phase)
- Note: Phase 4 cross-validation tests need updating for Lie-Trotter delta_eff constraint (separate task)

---
*Phase: 05-statistical-validation-and-regression*
*Completed: 2026-02-14*

## Self-Check: PASSED

- [x] test/trajectory_validation/run_convergence_tests.jl EXISTS
- [x] Commit 7e71503 EXISTS in git log
