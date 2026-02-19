---
phase: 26-exact-reference-and-structural-diagnostics
plan: 02
subsystem: diagnostics
tags: [observables, spectral-gap, eigenvectors, biorthogonal, lindbladian]

# Dependency graph
requires:
  - phase: 20-25 (v1.3 mixing time)
    provides: build_preset_trajectory_observables, eigenbasis_overlap_analysis, pad_term
provides:
  - Canonical 6-observable set [Z1, X1, Z1_Zhalf, H, Rand_traceless, Mz_stagg]
  - Biorthogonal left+right eigenvector overlap formula with rho_beta subtraction
affects: [27-reference-data-generation, 28-two-exponential-fitting, 29-batch-bootstrap, 30-validation]

# Tech tracking
tech-stack:
  added: []
  patterns: [biorthogonal-eigenvector-decomposition, seeded-random-observable]

key-files:
  created: []
  modified:
    - src/convergence.jl
    - src/gap_estimation.jl
    - test/test_convergence.jl
    - test/test_gap_estimation.jl

key-decisions:
  - "Canonical observable set: Z1, X1, Z1_Zhalf, H, Rand_traceless, Mz_stagg (replacing XX_avg, YY_avg, ZZ_avg, Mz, XZ_stagg)"
  - "Random traceless observable uses MersenneTwister(12345) for reproducibility, normalized by operator norm"
  - "Z1_Zhalf uses literal floor(n/2) formula: for n=2,3 this degenerates to Z1*Z1=I (valid but trivial)"
  - "rho_beta=nothing preserves v1.3 behavior for backward compatibility"
  - "Left eigenvectors computed via transpose(inv(V_right)), not conjugate transpose"

patterns-established:
  - "Biorthogonal overlap: c_k = Tr[O R_k] * Tr[L_k^dagger (rho_0 - rho_beta)]"
  - "Observable bundle versioning via canonical set replacement (not additive)"

# Metrics
duration: 6min
completed: 2026-02-19
---

# Phase 26 Plan 02: Observable Set and Overlap Formula Summary

**Canonical 6-observable set replacing v1.3 8-observable bundle, with biorthogonal left+right eigenvector overlap formula for proper steady-state subtraction**

## Performance

- **Duration:** 6 min
- **Started:** 2026-02-19T03:30:07Z
- **Completed:** 2026-02-19T03:35:51Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments
- Replaced 8-observable preset with physically motivated 6-observable canonical set from supplementary info
- Added X1 (Pauli-X site 1), Z1_Zhalf (two-point Z correlator), Rand_traceless (seeded random Hermitian)
- Updated eigenbasis_overlap_analysis with exact biorthogonal formula using left+right eigenvectors
- All 286 tests passing (200 convergence + 86 gap estimation)

## Task Commits

Each task was committed atomically:

1. **Task 1: Replace 8-observable set with 6-canonical set and update all tests** - `022a4f2` (feat)
2. **Task 2: Update eigenbasis_overlap_analysis with correct left+right eigenvector formula** - `77bceab` (feat)

## Files Created/Modified
- `src/convergence.jl` - Replaced build_preset_trajectory_observables with 6-observable canonical set
- `src/gap_estimation.jl` - Added rho_beta keyword to eigenbasis_overlap_analysis with biorthogonal formula
- `test/test_convergence.jl` - Updated all observable count/name references, added Trotter basis tests for new set
- `test/test_gap_estimation.jl` - Updated all references, added rho_beta subtraction and backward compatibility tests

## Decisions Made
- Canonical observable set from CONTEXT.md locked decisions: Z1, X1, Z1_Zhalf, H, Rand_traceless, Mz_stagg
- Random traceless observable: MersenneTwister(12345), Hermitianized, traceless, operator-norm normalized
- Z1_Zhalf: literal floor(n/2) formula as specified (degenerates to I for n=2,3)
- Left eigenvectors via transpose(inv(V_right)) -- transpose, not conjugate transpose, for proper biorthogonality
- Backward compatible: rho_beta=nothing preserves v1.3 alpha = V \ vec(rho0) formula

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Canonical 6-observable set ready for use by Phase 27 reference data generation
- eigenbasis_overlap_analysis ready with proper rho_beta subtraction for Phase 27 diagnostics
- All gap estimation tests updated and passing

## Self-Check: PASSED

All files verified present, all commit hashes found in git log.

---
*Phase: 26-exact-reference-and-structural-diagnostics*
*Completed: 2026-02-19*
