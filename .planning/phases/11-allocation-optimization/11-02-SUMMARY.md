---
phase: 11-allocation-optimization
plan: 02
subsystem: performance
tags: [diagonal-elimination, basis-transform, broadcasting, mul!, coherent-terms]

# Dependency graph
requires:
  - phase: 09-type-parameterization
    provides: "HamHam{T}, TrottTrott with typed eigvals/eigvecs"
provides:
  - "Diagonal-free B_time/B_trotter with pre-allocated vector buffers"
  - "B_trotter uses jump.in_eigenbasis directly (no redundant basis transform)"
  - "Callers ensure Trotter-basis in_eigenbasis via transform_jumps_to_basis"
affects: [11-03-PLAN, allocation-tests]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Vector buffer + element-wise .* broadcasting replaces Diagonal wrapper in loop bodies"
    - "transpose(vec) extracted outside @. macro to preserve row-vector semantics"
    - "Callers transform jumps to target eigenbasis before passing to B_trotter"

key-files:
  created: []
  modified:
    - src/coherent.jl
    - test/test_dm_scaling.jl

key-decisions:
  - "Used explicit .+= .* .* broadcasting instead of @. to avoid transpose-inside-@. pitfall"
  - "Test DMTST-05 updated to transform jump to Trotter basis before direct B_trotter call"

patterns-established:
  - "diag_u_row = transpose(diag_u) for column-scaling in fused broadcasts"
  - "transform_jumps_to_basis at caller site, not inside B_trotter"

# Metrics
duration: 7min
completed: 2026-02-15
---

# Phase 11 Plan 02: Diagonal Elimination and Redundant Transform Removal Summary

**Diagonal-free B_time/B_trotter via pre-allocated vector buffers and element-wise broadcasting, plus direct jump.in_eigenbasis usage in B_trotter per locked decision**

## Performance

- **Duration:** 7 min
- **Started:** 2026-02-15T15:21:06Z
- **Completed:** 2026-02-15T15:28:56Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments
- Eliminated all Diagonal() constructor calls from B_time and B_trotter loop bodies (4 variants)
- Replaced Diagonal wrappers with pre-allocated diag_u/diag_u2 vectors and element-wise .* broadcasting
- Removed redundant trotter.eigvecs' * jump.data * trotter.eigvecs recomputation in B_trotter
- Added transform_jumps_to_basis calls in all three TrotterDomain caller functions
- Fixed dim = size(matrix) tuple usage to d = size(matrix, 1) integer for correct vector buffer allocation
- All 224 tests pass with no regressions

## Task Commits

Each task was committed atomically:

1. **Task 1: Eliminate Diagonal wrappers in B_time and B_trotter** - `8e51827` (feat)
2. **Task 2: Use jump.in_eigenbasis directly in B_trotter** - `c497b05` (feat)

## Files Created/Modified
- `src/coherent.jl` - Replaced Diagonal closures with vector buffers + element-wise ops; B_trotter uses jump.in_eigenbasis directly; callers add transform_jumps_to_basis for TrotterDomain
- `test/test_dm_scaling.jl` - DMTST-05 updated to transform jump to Trotter basis before direct B_trotter call

## Decisions Made
- Used explicit `.+=` `.* .* .*` broadcasting instead of `@.` macro for the diagonal scaling operations, because `@.` converts `transpose(diag_u)` into element-wise `transpose.(diag_u)` which returns the vector unchanged (scalar transpose is identity), losing the row-vector shape needed for column-scaling
- Test DMTST-05 updated to match new B_trotter contract (expects jump.in_eigenbasis in Trotter basis)

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed @. macro transpose broadcasting semantics**
- **Found during:** Task 1
- **Issue:** Initial implementation used `@. b_plus_summand += b_s * diag_u * M * transpose(diag_u)` which converts `transpose(diag_u)` to `transpose.(diag_u)` (element-wise scalar transpose = identity), resulting in row-scale-only instead of row+column scaling
- **Fix:** Extracted `diag_u_row = transpose(diag_u)` outside the broadcast, used explicit `.+=` `.* .* .*` syntax
- **Files modified:** src/coherent.jl
- **Verification:** All 224 tests pass (4 had failed before fix)
- **Committed in:** 8e51827 (part of Task 1 commit)

**2. [Rule 3 - Blocking] Updated DMTST-05 test for new B_trotter contract**
- **Found during:** Task 2
- **Issue:** DMTST-05 calls B_trotter directly with a jump whose in_eigenbasis is in Hamiltonian basis, but B_trotter now expects Trotter-basis in_eigenbasis
- **Fix:** Added transform_jumps_to_basis call before direct B_trotter call in the test
- **Files modified:** test/test_dm_scaling.jl
- **Verification:** All 224 tests pass
- **Committed in:** c497b05 (part of Task 2 commit)

---

**Total deviations:** 2 auto-fixed (1 bug fix, 1 blocking)
**Impact on plan:** Both auto-fixes necessary for correctness. No scope creep.

## Issues Encountered
None beyond the deviations documented above.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- B_time and B_trotter are now Diagonal-free with pre-allocated buffers
- Ready for plan 11-03 (allocation tests and any remaining hotspot fixes)
- The jump_workers.jl filter intermediate fix (from plan 11-03 scope) is already present as an uncommitted change in the working tree

## Self-Check: PASSED

All files verified present. All commit hashes verified in git log.

---
*Phase: 11-allocation-optimization*
*Completed: 2026-02-15*
