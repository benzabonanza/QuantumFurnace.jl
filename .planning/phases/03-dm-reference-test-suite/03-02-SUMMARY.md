---
phase: 03-dm-reference-test-suite
plan: 02
subsystem: testing
tags: [lindbladian, euler-step, error-scaling, coherent-term, oft, density-matrix]

# Dependency graph
requires:
  - phase: 01-foundation-and-compilation
    provides: test_helpers.jl fixtures (TEST_HAM, TEST_JUMPS, TEST_TROTTER, make_liouv_config)
  - phase: 02-trajectory-bug-fixes
    provides: corrected construct_lindbladian, fixed EnergyDomain U_B ordering
provides:
  - DMTST-03: single-step Euler error O(delta^2) ratio test
  - DMTST-04: multi-step accumulated error O(delta) ratio test
  - DMTST-05: coherent term B cross-domain consistency verification
  - DMTST-06: OFT cross-domain consistency verification
affects: [04-trajectory-validation]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Liouvillian matrix exponential as exact DM evolution reference"
    - "Error ratio test for convergence order verification"
    - "Cross-domain B consistency with gamma_norm_factor normalization"
    - "OFT prefactor conventions: energy vs time domain"

key-files:
  created:
    - test/test_dm_scaling.jl
  modified:
    - test/runtests.jl

key-decisions:
  - "Used EnergyDomain Liouvillian for Euler scaling tests (simplest non-trivial domain)"
  - "Maximally mixed initial state for clean scaling measurement"
  - "Non-exported functions accessed via QuantumFurnace.func() for OFT tests"

patterns-established:
  - "Ratio test pattern: compute errors at halving delta values, verify ratios match expected convergence order"
  - "B comparison pattern: compute B_bohr, B_time, B_trotter with matching normalization, transform B_trotter to Hamiltonian eigenbasis"

# Metrics
duration: 6min
completed: 2026-02-14
---

# Phase 3 Plan 2: DM Scaling and Consistency Tests Summary

**Euler step error scaling O(delta^2)/O(delta) verified via ratio tests, coherent term B and OFT consistency confirmed across Bohr/Time/Trotter domains**

## Performance

- **Duration:** 6 min
- **Started:** 2026-02-14T09:32:19Z
- **Completed:** 2026-02-14T09:38:24Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments
- Single-step Euler error confirmed O(delta^2): ratios 3.97-3.99 across 4 delta values
- Multi-step accumulated error confirmed O(delta): ratios 2.01-2.03 across 4 delta values
- B_bohr matches B_time to 3.5e-14 (well within TOL_QUADRATURE=1e-6)
- B_trotter has 0.011 additional Trotter error, correctly larger than time quadrature error
- Energy OFT matches time OFT to 1.9e-12; trotter OFT has expected additional error

## Task Commits

Each task was committed atomically:

1. **Task 1: DM step error scaling tests (DMTST-03, DMTST-04)** - `c1d4aa3` (test)
2. **Task 2: Coherent term B and OFT consistency tests (DMTST-05, DMTST-06)** - `411adf2` (test)

## Files Created/Modified
- `test/test_dm_scaling.jl` - Four testsets: DMTST-03 (single-step O(delta^2)), DMTST-04 (multi-step O(delta)), DMTST-05 (coherent B consistency), DMTST-06 (OFT consistency)
- `test/runtests.jl` - Added include for test_dm_scaling.jl (committed by parallel Plan 01 execution)

## Decisions Made
- Used EnergyDomain (not BohrDomain) for Euler scaling tests: simplest non-trivial domain that exercises the full Liouvillian construction path
- Maximally mixed state (I/DIM) as initial condition: known nonzero distance to Gibbs, produces clean scaling without stochastic noise
- Accessed non-exported functions (time_oft!, trotter_oft!, create_energy_labels, truncate_time_labels_for_oft) via QuantumFurnace.func() qualification since they are not in the module's export list

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None - all tests passed on first execution with expected scaling ratios and error hierarchies.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- DM reference test suite (DMTST-01 through DMTST-06) complete
- Ready for Phase 4 trajectory validation which uses DM results as ground truth
- All error scaling and consistency properties verified quantitatively

## Self-Check: PASSED

- test/test_dm_scaling.jl: FOUND (4 testsets: DMTST-03, DMTST-04, DMTST-05, DMTST-06)
- test/runtests.jl: FOUND (includes test_dm_scaling.jl)
- Commit c1d4aa3: FOUND
- Commit 411adf2: FOUND

---
*Phase: 03-dm-reference-test-suite*
*Completed: 2026-02-14*
