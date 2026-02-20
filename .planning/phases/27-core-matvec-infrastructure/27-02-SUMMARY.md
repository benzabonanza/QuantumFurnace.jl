---
phase: 27-core-matvec-infrastructure
plan: 02
subsystem: krylov
tags: [lindbladian, matvec, round-trip, adjoint, zero-allocation, BLAS, energy-domain]

# Dependency graph
requires:
  - phase: 27-01
    provides: "KrylovWorkspace struct, apply_lindbladian!, apply_adjoint_lindbladian!, _accumulate_dissipator!"
provides:
  - "Module integration: krylov_workspace.jl and krylov_matvec.jl included and exported"
  - "Round-trip validated: matvec matches dense construct_lindbladian() to < 1e-12 for KMS, GNS, coherent, adjoint"
  - "Zero-allocation hot path verified via @allocated"
  - "make_liouv_config_gns factory for 4-qubit GNS test configs"
  - "Correct adjoint dissipator with proper anticommutator preservation"
affects: [28 (other-domain dispatch), 29 (KrylovKit integration)]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "BLAS.gemm!/axpy! for zero-allocation matrix operations (replaces mul!/broadcast)"
    - "Concrete-typed jump_eigenbases in KrylovWorkspace to avoid JumpOp abstract field boxing"
    - "_krylov_oft! inline OFT avoiding JumpOp field access in hot path"
    - "Separate _accumulate_adjoint_dissipator! preserving {L'L, rho} anticommutator"

key-files:
  created:
    - "test/test_krylov_matvec.jl"
  modified:
    - "src/QuantumFurnace.jl"
    - "src/krylov_workspace.jl"
    - "src/krylov_matvec.jl"
    - "test/test_helpers.jl"
    - "test/runtests.jl"

key-decisions:
  - "BLAS.gemm!/axpy! instead of mul!/broadcast to achieve zero heap allocations in inner loop"
  - "Store concrete-typed jump_eigenbases::Vector{Matrix{T}} in KrylovWorkspace to avoid JumpOp abstract field boxing"
  - "Separate adjoint dissipator function preserving {L'L, rho} anticommutator instead of naive L<->L' swap"

patterns-established:
  - "Zero-alloc Krylov matvec pattern: concrete-typed workspace fields + BLAS primitives + inline OFT"
  - "Adjoint dissipator correctness: sandwich L rho L' -> L' rho L with unchanged {L'L, rho}"

# Metrics
duration: 17min
completed: 2026-02-20
---

# Phase 27 Plan 02: Round-trip Tests and Module Integration Summary

**Module-integrated Krylov matvec with BLAS-based zero-allocation hot path, validated against dense Lindbladian to <1e-12 for KMS/GNS/coherent/adjoint across 10 random density matrices at n=4**

## Performance

- **Duration:** 17 min
- **Started:** 2026-02-20T11:33:37Z
- **Completed:** 2026-02-20T11:51:31Z
- **Tasks:** 2
- **Files modified:** 6

## Accomplishments
- Integrated krylov_workspace.jl and krylov_matvec.jl into QuantumFurnace module with proper exports
- Created comprehensive test suite (7 testsets, 54 assertions) covering all Phase 27 success criteria
- Fixed adjoint dissipator bug: naive L<->L' swap incorrectly changed anticommutator from {L'L, rho} to {LL', rho}
- Eliminated all heap allocations (156KB -> 0 bytes) by replacing mul!/broadcast with BLAS.gemm!/axpy! and extracting concrete-typed jump eigenbasis matrices

## Task Commits

Each task was committed atomically:

1. **Task 1: Module integration and test helper additions** - `5b74abe` (feat)
2. **Task 2: Round-trip correctness and allocation regression tests** - `cb5665e` (feat)

## Files Created/Modified
- `src/QuantumFurnace.jl` - Added include statements for krylov files and export block for KrylovWorkspace, apply_lindbladian!, apply_adjoint_lindbladian!
- `src/krylov_workspace.jl` - Added jump_eigenbases and jump_hermitian fields with concrete types to avoid abstract field boxing
- `src/krylov_matvec.jl` - Rewrote dissipator helpers with BLAS primitives; added _krylov_oft!, _accumulate_adjoint_dissipator!, and separate adj_L variants; fixed adjoint dissipator math
- `test/test_helpers.jl` - Added make_liouv_config_gns factory for 4-qubit GNS test configs
- `test/test_krylov_matvec.jl` - 7 testsets: workspace construction, KMS no-coherent, KMS coherent, GNS, adjoint, duality, zero-allocation
- `test/runtests.jl` - Added test_krylov_matvec.jl include

## Decisions Made
- Used `BLAS.gemm!`/`BLAS.axpy!` instead of `mul!` with Adjoint wrappers and `@.` broadcast: `mul!(A, B', C)` allocates 16 bytes per call due to Adjoint wrapper boxing; BLAS.gemm! is zero-allocation. Over 564 dissipator calls per matvec, this eliminated 156KB of heap allocations.
- Stored `jump_eigenbases::Vector{Matrix{T}}` with concrete element type in KrylovWorkspace: `JumpOp.in_eigenbasis` has abstract type `Matrix{<:Complex}` which causes boxing allocation on every field access. Extracting to concrete-typed vectors eliminates this.
- Created separate `_accumulate_adjoint_dissipator!` instead of reusing forward dissipator with swapped L: The Hilbert-Schmidt adjoint D_L*(rho) = L' rho L - 0.5{L'L, rho} has the SAME anticommutator {L'L, rho} as the forward, but passing L' to the forward helper gives D_{L'}(rho) = L' rho L - 0.5{LL', rho} which has a WRONG anticommutator.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed incorrect adjoint dissipator implementation**
- **Found during:** Task 2 (round-trip tests)
- **Issue:** The plan and plan-01 implementation used `_accumulate_dissipator!(out, L_op', rho, ...)` for the adjoint, which computes D_{L'}(rho) = L' rho L - 0.5{LL', rho}. The correct Hilbert-Schmidt adjoint is D_L*(rho) = L' rho L - 0.5{L'L, rho} -- the anticommutator must stay {L'L, rho}, not change to {LL', rho}.
- **Fix:** Created separate `_accumulate_adjoint_dissipator!` that computes L'L (not LL') for the anticommutator while using the adjoint sandwich L' rho L. Added corresponding `_adj_L` variants for the negative-frequency partner case.
- **Files modified:** src/krylov_matvec.jl
- **Verification:** All 10 random density matrices pass round-trip adjoint test < 1e-12; duality tr(X' L(Y)) == tr(L*(X)' Y) passes for 5 pairs < 1e-11
- **Committed in:** cb5665e (Task 2 commit)

**2. [Rule 1 - Bug] Eliminated JumpOp abstract field boxing allocations**
- **Found during:** Task 2 (allocation regression test)
- **Issue:** `JumpOp.in_eigenbasis` has type `Matrix{<:Complex}` (abstract element parameter), causing Julia to box the field value on every access -- 832 bytes per `oft!` call, totaling 156KB per matvec invocation.
- **Fix:** Added `jump_eigenbases::Vector{Matrix{T}}` and `jump_hermitian::Vector{Bool}` with concrete types to KrylovWorkspace. Created `_krylov_oft!` that takes eigenbasis directly. Replaced all `mul!` with Adjoint wrappers and `@.` broadcasts with BLAS.gemm!/axpy!.
- **Files modified:** src/krylov_workspace.jl, src/krylov_matvec.jl
- **Verification:** `@allocated apply_lindbladian!(...) == 0` and `@allocated apply_adjoint_lindbladian!(...) == 0`
- **Committed in:** cb5665e (Task 2 commit)

---

**Total deviations:** 2 auto-fixed (2 bugs)
**Impact on plan:** Both fixes essential for correctness (adjoint formula) and performance requirement (zero allocations). No scope creep.

## Issues Encountered
- None beyond the deviations documented above.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- Phase 27 complete: KrylovWorkspace + apply_lindbladian! + apply_adjoint_lindbladian! integrated, exported, tested
- Phase 28 can add TimeDomain/TrotterDomain/BohrDomain dispatch methods without touching workspace struct
- Phase 29 can wrap apply_lindbladian! in KrylovKit eigsolve closure for gap estimation
- Full test suite passes (996 tests) with no regressions

## Self-Check: PASSED

- [x] `src/QuantumFurnace.jl` includes krylov_workspace.jl, krylov_matvec.jl and exports KrylovWorkspace
- [x] `src/krylov_workspace.jl` has jump_eigenbases field with concrete type
- [x] `src/krylov_matvec.jl` has _accumulate_adjoint_dissipator! and BLAS.gemm! usage
- [x] `test/test_helpers.jl` has make_liouv_config_gns factory
- [x] `test/test_krylov_matvec.jl` exists with 7 testsets
- [x] `test/runtests.jl` includes test_krylov_matvec.jl
- [x] Commit 5b74abe exists (Task 1)
- [x] Commit cb5665e exists (Task 2)
- [x] Full test suite passes (996 tests)

---
*Phase: 27-core-matvec-infrastructure*
*Completed: 2026-02-20*
