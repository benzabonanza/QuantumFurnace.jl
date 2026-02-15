---
phase: 08-struct-simplification
plan: 02
subsystem: core-structs
tags: [julia, hamiltonian, struct-design, bson, initialization]

# Dependency graph
requires:
  - phase: 01-foundation-and-compilation
    provides: "Test infrastructure (test_helpers.jl), HamHam struct definition"
provides:
  - "Fully-initialized HamHam struct (no Nothing fields for bohr_freqs/bohr_dict/gibbs)"
  - "HamHam(NamedTuple, beta) constructor for BSON/raw data reconstruction"
  - "_gibbs_in_eigen helper for eigenbasis Gibbs state computation"
  - "find_ideal_heisenberg returns NamedTuple instead of partially-initialized HamHam"
  - "Legacy BSON loading via _load_hamiltonian_bson and _load_test_hamiltonian"
affects: [08-struct-simplification, simulation-scripts, test-infrastructure]

# Tech tracking
tech-stack:
  added: []
  patterns: ["Single-step initialization: all HamHam fields computed in constructor", "BSON.parse + raise_recursive for legacy deserialization"]

key-files:
  created: []
  modified:
    - "src/hamiltonian.jl"
    - "src/QuantumFurnace.jl"
    - "src/misc_tools.jl"
    - "test/test_helpers.jl"
    - "simulations/main_liouv.jl"
    - "simulations/main_thermalize.jl"

key-decisions:
  - "Inline Gibbs computation via _gibbs_in_eigen helper rather than modifying gibbs_state_in_eigen"
  - "Use BSON.parse + raise_recursive for legacy BSON deserialization instead of re-saving BSON files"
  - "load_hamiltonian now requires beta kwarg (breaking change, clean break per CONTEXT.md)"
  - "find_ideal_heisenberg returns NamedTuple for composability with HamHam(NamedTuple, beta)"

patterns-established:
  - "HamHam single-step init: all constructors take beta and compute bohr_freqs/bohr_dict/gibbs"
  - "Legacy BSON compat: BSON.parse -> extract fields -> HamHam(NamedTuple, beta)"

# Metrics
duration: 13min
completed: 2026-02-15
---

# Phase 8 Plan 2: HamHam Initialization Redesign Summary

**Eliminated two-step HamHam initialization: bohr_freqs/bohr_dict/gibbs computed directly in constructor via _gibbs_in_eigen helper, finalize_hamham deleted**

## Performance

- **Duration:** 13 min
- **Started:** 2026-02-15T08:58:22Z
- **Completed:** 2026-02-15T09:11:55Z
- **Tasks:** 2
- **Files modified:** 6

## Accomplishments
- HamHam struct fields bohr_freqs, bohr_dict, gibbs are no longer Union{..., Nothing} -- always fully initialized
- finalize_hamham function eliminated entirely (not exported, not defined)
- find_ideal_heisenberg returns NamedTuple of raw data instead of partially-initialized HamHam
- New HamHam(NamedTuple, beta) constructor for building from raw/serialized data
- Legacy BSON files handled via BSON.parse + manual field extraction (no re-serialization needed)
- All 224 tests pass with zero regressions, identical numerical results

## Task Commits

Each task was committed atomically:

1. **Task 1: Redesign HamHam struct and constructors** - `b612262` (feat)
2. **Task 2: Update all call sites and exports** - `97c0842` (feat)
3. **Task 2 cleanup: Remove finalize_hamham shim from parallel 08-01** - `fd83d2f` (fix)

**Plan metadata:** (pending)

## Files Created/Modified
- `src/hamiltonian.jl` - Redesigned HamHam struct (no Nothing fields), new constructors with beta, _gibbs_in_eigen helper, find_ideal_heisenberg -> NamedTuple, finalize_hamham deleted
- `src/QuantumFurnace.jl` - Removed finalize_hamham and load_hamiltonian_bson from exports
- `src/misc_tools.jl` - load_hamiltonian now requires beta kwarg, _load_hamiltonian_bson for legacy BSON
- `test/test_helpers.jl` - _load_test_hamiltonian for legacy BSON in test context, removed finalize_hamham calls
- `simulations/main_liouv.jl` - Updated to use load_hamiltonian with beta kwarg
- `simulations/main_thermalize.jl` - Updated to use load_hamiltonian with beta kwarg, simplified BSON loading

## Decisions Made
- **Inline Gibbs computation**: Created `_gibbs_in_eigen(eigvals, beta)` helper rather than modifying the existing `gibbs_state_in_eigen(HamHam, beta)` -- avoids circular dependency where constructor needs HamHam that doesn't exist yet
- **BSON.parse approach**: Used `BSON.parse` + `raise_recursive` for legacy deserialization rather than re-saving BSON files -- zero-touch approach to data migration
- **Breaking API change**: `load_hamiltonian` now requires `beta` kwarg -- clean break per CONTEXT.md decision "no backward compatibility shims"
- **NamedTuple return for find_ideal_heisenberg**: Enables composable construction pattern `HamHam(find_ideal_heisenberg(...), beta)`

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Legacy BSON deserialization handled via BSON.parse**
- **Found during:** Task 2 (test_helpers.jl update)
- **Issue:** BSON.load fails when struct definition changes (can't assign `nothing` to `Matrix{Float64}` field)
- **Fix:** Used `BSON.parse` + `raise_recursive` to load individual fields, then reconstruct via HamHam(NamedTuple, beta). Created `_load_hamiltonian_bson` in misc_tools.jl and `_load_test_hamiltonian` in test_helpers.jl
- **Files modified:** src/misc_tools.jl, test/test_helpers.jl
- **Verification:** All 224 tests pass with legacy BSON files
- **Committed in:** 97c0842 (part of Task 2 commit)

**2. [Rule 3 - Blocking] Parallel 08-01 execution overlap handled**
- **Found during:** Task 2 (commit phase)
- **Issue:** A parallel agent executing 08-01 had already committed changes to misc_tools.jl, test_helpers.jl, and simulation scripts that overlapped with 08-02's scope
- **Fix:** Identified the overlap, avoided double-committing the same changes, committed only the remaining export list cleanup
- **Files modified:** src/QuantumFurnace.jl
- **Verification:** git diff confirmed no content duplication
- **Committed in:** 97c0842

---

**Total deviations:** 2 auto-fixed (2 blocking)
**Impact on plan:** Both auto-fixes necessary for correctness. BSON deserialization fix was anticipated by the plan. Parallel execution overlap was handled gracefully.

## Issues Encountered
- BSON files store legacy HamHam with gibbs field containing a zero matrix (not nothing as expected) -- the stored gibbs was never the correct thermal state. This was a pre-existing data issue that is now irrelevant since gibbs is always recomputed from eigvals + beta in the new constructors.
- A parallel 08-01 plan execution committed some of the same files (test_helpers.jl, misc_tools.jl, simulation scripts) before this plan's Task 2 commit. The overlap was detected and handled by committing only the remaining export list change.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- HamHam is now fully initialized in one step -- no more two-step pattern
- find_ideal_heisenberg returns composable NamedTuple data
- All call sites updated for the new API
- Ready for further struct simplification (plans 08-03+)

## Self-Check: PASSED

- [x] All 6 source files exist on disk
- [x] Commit b612262 (Task 1) found in git log
- [x] Commit 97c0842 (Task 2) found in git log
- [x] Commit fd83d2f (Task 2 cleanup) found in git log
- [x] 08-02-SUMMARY.md exists
- [x] HamHam.bohr_freqs is Matrix{Float64} (not Union)
- [x] HamHam.gibbs is Hermitian{ComplexF64, Matrix{ComplexF64}} (not Union)
- [x] finalize_hamham not defined in QuantumFurnace module
- [x] find_ideal_heisenberg returns NamedTuple
- [x] All 224 tests pass

---
*Phase: 08-struct-simplification*
*Completed: 2026-02-15*
