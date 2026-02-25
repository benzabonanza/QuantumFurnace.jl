---
phase: 32-some-speedup-for-the-krylov-simulator
plan: 01
subsystem: krylov-simulator
tags: [blas, gemm, lindbladian, matvec, performance, precompute]

# Dependency graph
requires:
  - phase: 27-krylov-matvec-infrastructure
    provides: "KrylovWorkspace struct, _accumulate_dissipator! functions, BLAS.gemm! hot path"
  - phase: 28-krylov-domain-expansion
    provides: "TimeDomain, TrotterDomain, BohrDomain matvec functions"
provides:
  - "Precomputed G_left/G_right/G_left_adj/G_right_adj effective Hamiltonian fields in KrylovWorkspace"
  - "6 sandwich-only helper functions (_accumulate_sandwich!, _accumulate_sandwich_adj_L!, etc.)"
  - "Optimized 2+2N GEMM pattern replacing 2+5N per-matvec GEMM pattern across all 8 apply_lindbladian!/apply_adjoint_lindbladian! functions"
affects: [32-02-PLAN, krylov-eigsolve, krylov-benchmark]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Precomputed effective Hamiltonian G_left/G_right absorbing coherent + anticommutator terms"
    - "Sandwich-only inner loop (2 GEMMs per term) replacing full dissipator (5 GEMMs per term)"
    - "Separate G_left_adj/G_right_adj for BohrDomain non-Hermitian R_total correctness"

key-files:
  created: []
  modified:
    - "src/krylov_workspace.jl"
    - "src/krylov_matvec.jl"

key-decisions:
  - "Store 4 separate G matrices (G_left, G_right, G_left_adj, G_right_adj) rather than deriving adjoint from forward"
  - "BohrDomain adjoint uses conj(R_total) not R_total^T because R_total is non-Hermitian for Bohr"
  - "Pre-transpose G matrices at construction so matvec uses gemm!('N','N',...) for zero-allocation"

patterns-established:
  - "G_left/G_right effective Hamiltonian pattern: L(rho) = G_left*rho + rho*G_right + sandwiches"
  - "Sandwich-only helpers reuse existing tmp1/tmp2/LdagL scratch without new allocations"

# Metrics
duration: 13min
completed: 2026-02-25
---

# Phase 32 Plan 01: Precomputed Effective Hamiltonian Summary

**Precomputed G_left/G_right effective Hamiltonian matrices reducing per-matvec GEMM count from 5N to 2+2N across all 4 domains**

## Performance

- **Duration:** 13 min
- **Started:** 2026-02-25T08:10:42Z
- **Completed:** 2026-02-25T08:24:21Z
- **Tasks:** 3
- **Files modified:** 2

## Accomplishments
- Extended KrylovWorkspace with 4 precomputed G matrices computed at construction time from R_total and B_total
- Created 6 sandwich-only helper functions (2 GEMMs each, down from 5 in full dissipator)
- Refactored all 8 apply_lindbladian!/apply_adjoint_lindbladian! functions to use G_left/G_right + sandwich-only inner loop
- All 198 matvec round-trip, duality, and zero-allocation tests pass unchanged
- All 47 eigsolve tests pass unchanged
- All 8 cross-validation round-trip tests pass unchanged

## Task Commits

Each task was committed atomically:

1. **Task 1: Extend KrylovWorkspace struct and constructors** - `ba6d5a3` (feat)
2. **Task 2: Create sandwich-only helper functions** - `1c62244` (feat)
3. **Task 3: Refactor all apply_lindbladian!/apply_adjoint_lindbladian!** - `0ef033d` (feat)

## Files Created/Modified
- `src/krylov_workspace.jl` - Added G_left/G_right/G_left_adj/G_right_adj fields to struct; both constructors compute and store them at construction time using existing _accumulate_R_total!
- `src/krylov_matvec.jl` - Added 6 sandwich-only helpers; refactored all 8 apply_lindbladian!/apply_adjoint_lindbladian! to use 2+2N GEMM pattern

## Decisions Made
- Store 4 separate G matrices rather than 2 + derivation: eliminates subtle conjugation errors for BohrDomain where R_total is non-Hermitian
- For Energy/Time/Trotter domains (Hermitian R_total): G_left_adj = G_right, G_right_adj = G_left (simple pointer sharing)
- For BohrDomain (non-Hermitian R_total): G_left_adj = -i*B^T - 0.5*conj(R_total), computed independently
- Pre-transpose all G matrices at construction time so hot path uses only gemm!('N','N',...) avoiding Adjoint/Transpose wrapper boxing

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- Optimized matvec foundation in place for Plan 32-02 (dead code removal of legacy Euler apply_delta_channel!)
- All existing tests verified passing -- optimization is semantically transparent

## Self-Check: PASSED

All artifacts verified:
- src/krylov_workspace.jl: FOUND, contains G_left field
- src/krylov_matvec.jl: FOUND, contains _accumulate_sandwich! and ws.G_left usage
- Commit ba6d5a3: FOUND (Task 1)
- Commit 1c62244: FOUND (Task 2)
- Commit 0ef033d: FOUND (Task 3)

---
*Phase: 32-some-speedup-for-the-krylov-simulator*
*Completed: 2026-02-25*
