---
phase: 34-code-deduplication
plan: 01
subsystem: simulation-core
tags: [julia, multiple-dispatch, deduplication, OFT, domain-prefactor]

# Dependency graph
requires:
  - phase: 33-type-foundation
    provides: "Unified Config{S,D,C,T} type with domain singletons"
  - phase: quick-39
    provides: "Fixed krylov_matvec sandwich convention to L*rho*L'"
provides:
  - "domain_prefactor() function with 3 methods (Energy, Time, Trotter) in furnace_utensils.jl"
  - "domain_prefactor field in precomputed_data NamedTuple for Energy/Time/Trotter domains"
  - "Unified oft!(out, eigenbasis, bohr_freqs, energy, inv_4sigma2) replacing old oft! and _krylov_oft!"
affects: [34-02-PLAN, Phase 35, Phase 36]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "domain_prefactor stored in precomputed_data NamedTuple, callers multiply by their own scaling"
    - "Single oft! with concrete-typed signature, inv_4sigma2 computed at call site"

key-files:
  created: []
  modified:
    - src/furnace_utensils.jl
    - src/ofts.jl
    - src/jump_workers.jl
    - src/krylov_workspace.jl
    - src/krylov_matvec.jl
    - src/krylov_eigsolve.jl
    - src/trajectories.jl
    - test/test_dm_scaling.jl
    - test/old_tests/time_tests.jl

key-decisions:
  - "domain_prefactor returns only domain-dependent scalar; callers compose with gamma_norm_factor or jump_weight_scaling"
  - "Old JumpOp-based oft! deleted entirely; test callers updated to new signature"
  - "inv_4sigma2 computed at call site, not stored in precomputed_data (per CONTEXT.md)"
  - "time_oft! and trotter_oft! preserved unchanged as test/debug utilities"

patterns-established:
  - "Precomputed domain_prefactor in NamedTuple: callers use precomputed_data.domain_prefactor * scaling"
  - "Unified oft! concrete-typed signature: oft!(out, eigenbasis, bohr_freqs, energy, inv_4sigma2)"

# Metrics
duration: 12min
completed: 2026-02-26
---

# Phase 34 Plan 01: Extract domain_prefactor and Unify OFT Summary

**domain_prefactor() centralizing 16 inline formulas across 5 files and unified oft! replacing both old JumpOp-based oft! and _krylov_oft!**

## Performance

- **Duration:** 12 min
- **Started:** 2026-02-26T12:58:23Z
- **Completed:** 2026-02-26T13:10:18Z
- **Tasks:** 2
- **Files modified:** 9

## Accomplishments
- Extracted `domain_prefactor()` with 3 methods (EnergyDomain, TimeDomain, TrotterDomain) as single source for domain-dependent scalar prefactor
- Stored `domain_prefactor` in `precomputed_data` NamedTuple for zero-overhead access at all 16 call sites
- Unified `oft!` and `_krylov_oft!` into single concrete-typed `oft!(out, eigenbasis, bohr_freqs, energy, inv_4sigma2)` function
- Deleted old JumpOp-based `oft!` from ofts.jl and `_krylov_oft!` from krylov_matvec.jl
- All 1198 tests pass with zero numerical regressions

## Task Commits

Each task was committed atomically:

1. **Task 1: Extract domain_prefactor and store in precomputed_data** - `02cb3ea` (feat)
2. **Task 2: Unify oft! and _krylov_oft! into single function** - `61871fb` (feat)

## Files Created/Modified
- `src/furnace_utensils.jl` - Added domain_prefactor() with 3 domain-dispatched methods; stored in _precompute_data NamedTuple
- `src/ofts.jl` - Replaced old JumpOp-based oft! with new concrete-typed signature
- `src/jump_workers.jl` - Updated 4 prefactor sites to precomputed_data.domain_prefactor; updated 4 oft! call sites to new signature
- `src/krylov_workspace.jl` - Updated 2 prefactor sites; replaced 2 inline OFT expressions with oft!()
- `src/krylov_matvec.jl` - Updated 4 prefactor sites; deleted _krylov_oft! function; replaced 8 _krylov_oft! calls with oft!
- `src/krylov_eigsolve.jl` - Updated 2 prefactor sites; replaced 3 _krylov_oft! calls with oft!
- `src/trajectories.jl` - Updated 4 prefactor sites (including build_trajectoryframework simplification); updated 4 oft! call sites to new signature
- `test/test_dm_scaling.jl` - Updated 2 oft! test calls to new signature
- `test/old_tests/time_tests.jl` - Updated 1 oft! test call to new signature

## Decisions Made
- `domain_prefactor` returns only the domain-dependent scalar (no gamma_norm_factor). Callers compose: `domain_prefactor * gamma_norm_factor` or `domain_prefactor * jump_weight_scaling`. This preserves the distinction between Thermalize (uses jump_weight_scaling) and other paths.
- `inv_4sigma2` computed at call site per CONTEXT.md decision (not stored in precomputed_data).
- Old JumpOp-based `oft!` fully deleted per user decision in CONTEXT.md, not kept as compatibility wrapper. Test callers updated to new signature.
- `build_trajectoryframework` simplified: EnergyDomain/TimeDomain/TrotterDomain branches collapsed into single `precomputed_data.domain_prefactor * gamma_norm_factor / (1.0 / n_jumps)` expression.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Updated test callers of old oft! signature**
- **Found during:** Task 2 (OFT unification)
- **Issue:** test/test_dm_scaling.jl (2 sites) and test/old_tests/time_tests.jl (1 site) called old JumpOp-based `oft!` signature which was deleted
- **Fix:** Updated to new signature: `oft!(out, jump.in_eigenbasis, hamiltonian.bohr_freqs, w, 1.0 / (4 * SIGMA^2))`
- **Files modified:** test/test_dm_scaling.jl, test/old_tests/time_tests.jl
- **Verification:** All 1198 tests pass
- **Committed in:** 61871fb (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (1 blocking)
**Impact on plan:** Essential fix -- tests would fail without updating callers of deleted function. No scope creep.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- domain_prefactor and unified oft! are in place, ready for Plan 02 (cross-simulation-type deduplication: _accumulate_R!, _build_cptp_channel, sandwich consolidation)
- precomputed_data NamedTuple now includes domain_prefactor field for Energy/Time/Trotter domains
- No blockers

## Self-Check: PASSED

All 9 modified files verified present. Both task commits (02cb3ea, 61871fb) verified in git log.

---
*Phase: 34-code-deduplication*
*Completed: 2026-02-26*
