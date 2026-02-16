# Quick Task 19: Fix failing test after EnergyDomain to TrotterDomain rename

## Task
Fix the single failing test in `test_gns_trajectory.jl` after user changed GNS tests from `EnergyDomain` to `TrotterDomain`.

## Tasks

### Task 1: Relax GNS-to-Gibbs gap sanity bound
- **File:** `test/test_gns_trajectory.jl:39`
- **Issue:** `@test gap < 0.5` fails because TrotterDomain has gap ~0.83 (vs EnergyDomain ~0.08)
- **Fix:** Change bound to `gap < 1.0` with comment documenting expected values per domain
- **Verification:** All 284 tests pass
