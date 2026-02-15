---
phase: 11-allocation-optimization
plan: 01
subsystem: simulation-core
tags: [allocation, sparse-matrix, index-accumulation, half-grid, performance]

# Dependency graph
requires:
  - phase: 09-type-parameterization
    provides: "Generic HamHam{T} and parameterized jump operators"
provides:
  - "Index-based accumulation B_bohr (single-jump and multi-jump)"
  - "Half-grid continue pattern in Time/Trotter thermalize _jump_contribution!"
affects: [11-02-allocation-optimization, 11-03-allocation-optimization]

# Tech tracking
tech-stack:
  added: []
  patterns: ["Index-based sparse accumulation via @inbounds CartesianIndex loops", "Half-grid continue pattern for hermitian symmetry exploitation"]

key-files:
  created: []
  modified:
    - "src/bohr_domain.jl"
    - "src/jump_workers.jl"

key-decisions:
  - "No SparseArrays import changes needed -- spzeros was only SparseArrays usage in bohr_domain.jl"

patterns-established:
  - "Index-based accumulation: loop over CartesianIndex entries from bohr_dict instead of constructing sparse matrices"
  - "Half-grid continue: w_raw > 1e-12 && continue with w = abs(w_raw) for hermitian symmetry exploitation"

# Metrics
duration: 6min
completed: 2026-02-15
---

# Phase 11 Plan 01: Hot-Loop Allocation Elimination Summary

**Index-based accumulation in B_bohr eliminating O(num_freqs) sparse matrix allocations, and half-grid continue pattern in Time/Trotter thermalize eliminating filter+abs vector allocations**

## Performance

- **Duration:** 6 min
- **Started:** 2026-02-15T15:21:00Z
- **Completed:** 2026-02-15T15:27:00Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments
- Replaced spzeros + scatter + transpose multiply in both B_bohr variants with direct @inbounds loops over CartesianIndex entries from bohr_dict
- Replaced abs.(filter(w -> w < 1e-12, energy_labels)) temporary vector allocation with inline half-grid continue pattern matching existing EnergyDomain and Liouvillian variants
- All 224 tests pass with numerically identical results (pure refactoring)

## Task Commits

Each task was committed atomically:

1. **Task 1: Replace sparse matrix allocation in B_bohr with index-based accumulation** - `4a42ed3` (feat)
2. **Task 2: Replace filter+abs allocation in Time/Trotter thermalize with half-grid continue** - `800f1a5` (feat)

## Files Created/Modified
- `src/bohr_domain.jl` - B_bohr single-jump and multi-jump variants now use index-based accumulation instead of spzeros
- `src/jump_workers.jl` - Time/Trotter thermalize _jump_contribution! uses half-grid continue pattern with explicit hermitian/non-hermitian branches

## Decisions Made
- No SparseArrays import changes needed: spzeros was the only SparseArrays usage in bohr_domain.jl; the main module import stays for other files (jump_workers.jl uses sparse())

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- Allocation hotspots in B_bohr and Time/Trotter thermalize eliminated
- Ready for plan 02 (further allocation optimization targets) and plan 03
- Patterns established here (index-based accumulation, half-grid continue) can be applied to remaining allocation hotspots

## Self-Check: PASSED

- FOUND: src/bohr_domain.jl
- FOUND: src/jump_workers.jl
- FOUND: 11-01-SUMMARY.md
- FOUND: 4a42ed3 (Task 1 commit)
- FOUND: 800f1a5 (Task 2 commit)

---
*Phase: 11-allocation-optimization*
*Completed: 2026-02-15*
