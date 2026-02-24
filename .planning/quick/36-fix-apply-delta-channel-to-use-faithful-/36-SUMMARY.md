---
phase: quick-36
plan: 01
subsystem: krylov-eigsolve
tags: [chen-channel, cptp, spectral-gap, krylov, lindbladian]

# Dependency graph
requires:
  - phase: 29-eigensolver-integration
    provides: "KrylovWorkspace, apply_delta_channel!, krylov_spectral_gap channel path"
provides:
  - "Faithful Chen CPTP channel in apply_delta_channel! (Eq. 3.2)"
  - "KrylovWorkspace ThermalizeConfig constructor with precomputed K0, U_residual, U_coherent"
  - "Domain-dispatched _accumulate_jump_sandwich! for EnergyDomain, TimeDomain, TrotterDomain, BohrDomain"
  - "_accumulate_R_total! helpers for all domains"
affects: [30-cross-validation, channel-eigsolve, thermalization]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Physics-convention sandwich (L*rho*L') for channel, kron-convention for Lindbladian matvec"
    - "Precompute rho-independent channel matrices (K0, U_residual, U_coherent) at workspace construction"

key-files:
  created: []
  modified:
    - "src/krylov_workspace.jl"
    - "src/krylov_eigsolve.jl"
    - "test/test_krylov_eigsolve.jl"

key-decisions:
  - "alpha_chen variable name to avoid shadowing BohrDomain alpha function"
  - "Move _thermalize_to_liouv_config to krylov_workspace.jl to avoid circular dep"
  - "Physics convention for channel sandwiches (matching thermalization code), kron convention for Lindbladian matvec"
  - "Retain legacy 5-arg Euler apply_delta_channel! for backward compatibility"
  - "Relax channel eigsolve rtol from 1e-3 to 2e-3 (faithful channel O(delta^2) eigenvalue mapping error)"
  - "BohrDomain sandwich uses entrywise rho*A_nu2_dag scatter (matching jump_workers.jl thermalization code)"

# Metrics
duration: 12min
completed: 2026-02-24
---

# Quick-36: Faithful Chen CPTP Channel Summary

**Replace Euler-approximate apply_delta_channel! with Chen Eq. 3.2 faithful CPTP channel using precomputed K0, U_residual, U_coherent matrices**

## Performance

- **Duration:** 12 min
- **Started:** 2026-02-24T13:33:19Z
- **Completed:** 2026-02-24T13:45:13Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments
- KrylovWorkspace now has 5 optional channel fields (K0, U_residual, U_coherent, rho_jump scratch, delta)
- ThermalizeConfig constructor precomputes R_total, K0, S, U_residual at workspace creation time
- apply_delta_channel! implements E(rho) = K0*rho*K0' + rho_jump + U_res*rho*U_res' with coherent unitary rotation
- Channel is exactly trace-preserving and positive (CPTP), matching run_thermalization
- O(delta^2) agreement with Euler approximation verified across random density matrices
- All 1187 tests pass including new faithful channel tests

## Task Commits

Each task was committed atomically:

1. **Task 1: Add channel fields to KrylovWorkspace and ThermalizeConfig constructor** - `5f3c80b` (feat)
2. **Task 2: Rewrite apply_delta_channel! and update channel eigsolve path** - `f761e79` (feat)

## Files Created/Modified
- `src/krylov_workspace.jl` - Added 5 channel fields, ThermalizeConfig constructor, _accumulate_R_total! for all domains, moved _thermalize_to_liouv_config here
- `src/krylov_eigsolve.jl` - Rewrote apply_delta_channel! with Chen channel, added _accumulate_jump_sandwich! for all domains, updated channel eigsolve path
- `test/test_krylov_eigsolve.jl` - New faithful channel testset (trace preservation, positivity, O(delta^2)), kept legacy Euler test, relaxed channel eigsolve tolerance

## Decisions Made
- **alpha_chen naming:** Used `alpha_chen` for `1 - sqrt(1 - delta)` to avoid shadowing the `alpha` function in BohrDomain precomputed_data
- **Config conversion location:** Moved `_thermalize_to_liouv_config` to krylov_workspace.jl (included before krylov_eigsolve.jl) to avoid circular dependency
- **Physics convention for sandwiches:** Channel sandwich terms use `L * rho * L'` (physics convention) matching _jump_contribution! in jump_workers.jl, while Lindbladian matvec uses `conj(L) * rho * L^T` (kron convention)
- **Legacy compatibility:** Retained 5-arg `apply_delta_channel!(ws, rho, delta, config, ham)` Euler form for backward compatibility
- **Channel tolerance:** Relaxed Testset 5 rtol from 1e-3 to 2e-3 because the faithful channel's eigenvalues differ from 1+delta*lambda_L by O(delta^2), giving ~0.13% relative error in spectral gap

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Channel eigsolve tolerance too tight for faithful channel**
- **Found during:** Task 2 (test verification)
- **Issue:** rtol=1e-3 failed because faithful channel eigenvalue mapping introduces O(delta^2) error (~0.134%) vs exact linear Euler mapping
- **Fix:** Relaxed rtol to 2e-3 with explanatory comment documenting the error source
- **Files modified:** test/test_krylov_eigsolve.jl
- **Verification:** All tests pass
- **Committed in:** f761e79 (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (1 bug fix)
**Impact on plan:** Tolerance adjustment necessary due to physics of eigenvalue mapping. No scope creep.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Faithful channel ready for cross-validation (Phase 30)
- Channel eigsolve path now matches thermalization code exactly (same CPTP map)
- All domain dispatches implemented (EnergyDomain, TimeDomain, TrotterDomain, BohrDomain)

---
*Quick task: 36-fix-apply-delta-channel-to-use-faithful-*
*Completed: 2026-02-24*

## Self-Check: PASSED

All files exist, all commits verified.
