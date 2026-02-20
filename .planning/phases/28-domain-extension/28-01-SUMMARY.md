---
phase: 28-domain-extension
plan: 01
subsystem: krylov-matvec
tags: [lindbladian, matvec, nufft, time-domain, trotter-domain, zero-alloc]

# Dependency graph
requires:
  - phase: 27-core-matvec-infrastructure
    provides: "EnergyDomain apply_lindbladian!/apply_adjoint_lindbladian!, KrylovWorkspace, dissipator helpers"
provides:
  - "apply_lindbladian! for Union{TimeDomain, TrotterDomain}"
  - "apply_adjoint_lindbladian! for Union{TimeDomain, TrotterDomain}"
  - "Round-trip tests (KMS, GNS) and zero-allocation regression tests for both domains"
affects: [28-02, 29-krylov-integration, 30-gap-estimation]

# Tech tracking
tech-stack:
  added: []
  patterns: ["NUFFT prefactor multiply for Time/Trotter OFT", "_measure_matvec_allocs helper for soft-scope-safe allocation tests"]

key-files:
  created: []
  modified:
    - src/krylov_matvec.jl
    - test/test_krylov_matvec.jl

key-decisions:
  - "Function wrapper _measure_matvec_allocs to avoid @testset soft-scope variable boxing causing spurious 176-byte allocations"
  - "Same scalar prefactor for forward and adjoint (w0 * t0^2 * sigma * sqrt(2/pi) / (2pi) * gamma_norm_factor)"

patterns-established:
  - "_prefactor_view + broadcast multiply for zero-alloc NUFFT OFT in Krylov hot path"
  - "_measure_matvec_allocs helper pattern for allocation regression tests"

# Metrics
duration: 19min
completed: 2026-02-20
---

# Phase 28 Plan 01: Time/Trotter Domain Matvec Summary

**Zero-allocation matrix-free Lindbladian action for TimeDomain and TrotterDomain using NUFFT prefactors, with round-trip verification < 1e-12 against dense construction**

## Performance

- **Duration:** 19 min
- **Started:** 2026-02-20T14:00:56Z
- **Completed:** 2026-02-20T14:20:35Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments
- `apply_lindbladian!` and `apply_adjoint_lindbladian!` for `Union{TimeDomain, TrotterDomain}` -- mechanical translation from EnergyDomain with NUFFT prefactor OFT and t0^2 scalar prefactor
- 8 new testsets: KMS forward, GNS forward, adjoint, and zero-allocation for both TimeDomain and TrotterDomain
- All round-trip errors < 1e-12 across 10 random density matrices at n=4
- Zero heap allocations verified for both forward and adjoint hot paths in both domains

## Task Commits

Each task was committed atomically:

1. **Task 1: Implement Time/Trotter forward and adjoint matvec methods** - `e9560f8` (feat)
2. **Task 2: Add round-trip and allocation tests for TimeDomain and TrotterDomain** - `1bcbc81` (test)

## Files Created/Modified
- `src/krylov_matvec.jl` - Added 2 new method dispatches for `Union{TimeDomain, TrotterDomain}` (forward + adjoint), using `_prefactor_view` for NUFFT OFT and `w0 * t0^2 * sigma * sqrt(2/pi) / (2pi) * gamma_norm_factor` scalar prefactor
- `test/test_krylov_matvec.jl` - Added 8 new testsets (TimeDomain: KMS forward, GNS forward, adjoint, zero-alloc; TrotterDomain: KMS forward, GNS forward, adjoint, zero-alloc) plus `_measure_matvec_allocs` helper function

## Decisions Made
- Used `_measure_matvec_allocs` helper function wrapper instead of inline `@allocated` to avoid Julia's `@testset` soft-scope variable boxing. The `@testset` macro creates a soft scope where variables shared with `@allocated` can become boxed (176 bytes), producing spurious allocation counts. Wrapping the measurement in a function creates a hard scope, ensuring zero-allocation verification is accurate.
- Same scalar prefactor formula used for both forward and adjoint (no prefactor changes for adjoint, per CONTEXT.md)

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed @testset soft-scope boxing in allocation tests**
- **Found during:** Task 2 (allocation test verification)
- **Issue:** `@allocated apply_adjoint_lindbladian!(...)` inside nested `@testset` blocks reported 176 bytes due to Julia's soft-scope variable boxing. The actual function produces zero allocations (verified at global scope and in function contexts).
- **Fix:** Extracted allocation measurement into a `_measure_matvec_allocs` helper function that creates a hard scope via function boundary, preventing variable boxing.
- **Files modified:** test/test_krylov_matvec.jl
- **Verification:** All 1060 tests pass including zero-allocation assertions
- **Committed in:** 1bcbc81 (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (1 bug)
**Impact on plan:** Auto-fix necessary for test correctness. No scope creep. The underlying matvec code is genuinely zero-allocation.

## Issues Encountered
None beyond the @testset soft-scope issue documented above.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- TimeDomain and TrotterDomain matvec complete with full parity to EnergyDomain
- Ready for Plan 02 (BohrDomain matvec with 2-operator dissipator helper)
- All 1060 tests pass (996 existing + 64 new)

## Self-Check: PASSED

- [x] src/krylov_matvec.jl exists
- [x] test/test_krylov_matvec.jl exists
- [x] 28-01-SUMMARY.md exists
- [x] Commit e9560f8 found
- [x] Commit 1bcbc81 found

---
*Phase: 28-domain-extension*
*Completed: 2026-02-20*
