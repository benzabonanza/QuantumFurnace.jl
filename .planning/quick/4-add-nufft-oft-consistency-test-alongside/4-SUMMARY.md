---
phase: quick-4
plan: 01
subsystem: testing
tags: [nufft, oft, consistency, quantum-furnace, julia]

# Dependency graph
requires:
  - phase: 03-dm-reference-test-suite
    provides: "DMTST-06 OFT consistency test framework and test helpers"
provides:
  - "DMTST-06b NUFFT OFT consistency testset verifying NUFFT acceleration faithfulness"
affects: [phase-04, nufft, oft]

# Tech tracking
tech-stack:
  added: []
  patterns: ["NUFFT prefactor extraction via precompute_data + prefactor_view for test verification"]

key-files:
  created: []
  modified: ["test/test_dm_scaling.jl"]

key-decisions:
  - "Used same test energy w=-3*W0 as DMTST-06 for direct comparability"
  - "Added haskey sanity checks to catch energy grid mismatches early"

patterns-established:
  - "NUFFT test pattern: precompute_data -> prefactor_view -> elementwise multiply -> compare with direct method"

# Metrics
duration: 1min
completed: 2026-02-14
---

# Quick Task 4: Add NUFFT OFT Consistency Test Summary

**DMTST-06b testset verifying NUFFT-accelerated OFT matches direct summation (time_oft!, trotter_oft!) and analytical OFT across Time and Trotter domains**

## Performance

- **Duration:** 1 min
- **Started:** 2026-02-14T11:07:12Z
- **Completed:** 2026-02-14T11:08:23Z
- **Tasks:** 1
- **Files modified:** 1

## Accomplishments
- Added DMTST-06b testset with 6 assertions (2 haskey + 4 numerical comparisons)
- NUFFT time OFT matches direct time_oft! at 3.3e-13 (threshold: 1e-10)
- NUFFT trotter OFT matches direct trotter_oft! at 4.4e-12 (threshold: 1e-10)
- NUFFT time OFT matches analytical oft! at 2.0e-12 (threshold: TOL_QUADRATURE = 1e-6)
- NUFFT trotter OFT matches analytical at 1.5e-8 (threshold: 0.1)

## Task Commits

Each task was committed atomically:

1. **Task 1: Add DMTST-06b NUFFT OFT consistency testset** - `c4a2126` (test)

## Files Created/Modified
- `test/test_dm_scaling.jl` - Added DMTST-06b testset (80 lines) after existing DMTST-06

## Decisions Made
- Used same test energy `w = -3 * W0` as DMTST-06 for direct comparability
- Added `haskey` sanity checks before accessing NUFFT prefactors to catch energy grid mismatches early
- Used manually reconstructed `oft_time_labels` (same as DMTST-06) for direct comparison with time_oft!/trotter_oft!

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- NUFFT OFT consistency is now verified alongside direct and analytical methods
- All existing tests (DMTST-03 through DMTST-06) remain passing

## Self-Check: PASSED

- FOUND: test/test_dm_scaling.jl
- FOUND: commit c4a2126
- FOUND: DMTST-06b in test file

---
*Quick Task: 4-add-nufft-oft-consistency-test-alongside*
*Completed: 2026-02-14*
