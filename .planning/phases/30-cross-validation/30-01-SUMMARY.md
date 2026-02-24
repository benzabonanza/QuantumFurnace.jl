---
phase: 30-cross-validation
plan: 01
subsystem: testing
tags: [krylov, cross-validation, spectral-gap, dense-eigen, convergence]

# Dependency graph
requires:
  - phase: 29-eigensolver-integration
    provides: krylov_spectral_gap(), KrylovGapResult, extract_leading_eigendata()
  - phase: 27-krylov-matvec-workspace
    provides: KrylovWorkspace, apply_lindbladian!
provides:
  - Cross-validation test file with Krylov vs dense helpers
  - n=4 KMS/GNS gap comparison across all 4 domains at atol=1e-8
  - L-vs-E convergence analysis with O(delta^2) order assertion
  - Diagnostic output helpers (gap summary, eigenvalue table, failure diagnostics)
affects: [30-02, 31-benchmarks]

# Tech tracking
tech-stack:
  added: []
  patterns: [compare_krylov_dense helper pattern, L-vs-E convergence table format]

key-files:
  created:
    - test/test_krylov_crossvalidation.jl
  modified:
    - test/runtests.jl

key-decisions:
  - "atol=1e-8 for n=4 cross-validation (KrylovKit tol=1e-10 provides margin)"
  - "L-vs-E convergence order >= 1.5 hard assertion (O(delta^2) with sub-leading margin)"
  - "GNS via make_liouv_config_gns (with_coherent=false) per struct guard"
  - "TrotterDomain always uses TEST_TROTTER_JUMPS + trotter=TEST_TROTTER"

patterns-established:
  - "compare_krylov_dense: shared helper returning (krylov_result, dense_result, L_dense)"
  - "Diagnostic pattern: always print gap summary, print eigenvalue table only on failure"
  - "run_le_convergence: formatted convergence table with order computation from consecutive error pairs"

# Metrics
duration: 2min
completed: 2026-02-24
---

# Phase 30 Plan 01: Cross-Validation Test Infrastructure Summary

**Krylov vs dense eigen() cross-validation across 4 domains, KMS/GNS balance, and L-vs-E O(delta^2) convergence at n=4**

## Performance

- **Duration:** 2 min
- **Started:** 2026-02-24T14:27:05Z
- **Completed:** 2026-02-24T14:29:07Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments
- Created comprehensive cross-validation test file with 5 helper functions and 3 testset blocks
- XVAL-01: n=4 KMS gap comparison for EnergyDomain, TimeDomain, TrotterDomain, BohrDomain at atol=1e-8
- XVAL-04: n=4 GNS gap comparison for all 4 domains at atol=1e-8
- XVAL-03: L-vs-E convergence test across 3 deltas (0.1, 0.01, 0.001) with order >= 1.5 assertion for all 4 domains
- Diagnostic output: one-line gap summary always printed, top-6 eigenvalue table on failure

## Task Commits

Each task was committed atomically:

1. **Task 1: Create test_krylov_crossvalidation.jl with helpers and tests** - `06aabc2` (feat)
2. **Task 2: Add test_krylov_crossvalidation.jl to runtests.jl** - `3abcc0e` (chore)

## Files Created/Modified
- `test/test_krylov_crossvalidation.jl` - Cross-validation test file with compare_krylov_dense helper, diagnostic printing, L-vs-E convergence, n=4 KMS/GNS testsets
- `test/runtests.jl` - Added include for new test file as last entry in QuantumFurnace.jl testset

## Decisions Made
None - followed plan as specified. All implementation decisions (atol=1e-8, order >= 1.5, top_k=6, krylovdim=30) were locked in CONTEXT.md and PLAN.md.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
- Git index corruption after Task 2 commit; resolved by removing corrupt `.git/index` and running `git reset`. Both commits preserved intact.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Cross-validation test infrastructure complete for n=4
- Plan 30-02 can add n=6 tests (env-gated behind QUANTUMFURNACE_FULL_TESTS=true) using the same compare_krylov_dense and run_le_convergence helpers
- Placeholder comment in test file marks the n=6 insertion point

## Self-Check: PASSED

All files exist, all commits verified, all content checks pass.

---
*Phase: 30-cross-validation*
*Completed: 2026-02-24*
