---
phase: quick-39
plan: 01
subsystem: physics
tags: [lindbladian, krylov, convention, dissipator, sandwich]

# Dependency graph
requires:
  - phase: 32
    provides: Krylov sandwich-only helpers with precomputed G_left/G_right
provides:
  - Consistent L*rho*L' convention across dense Lindbladian, Krylov matvec, and Thermalize paths
  - Correct HS adjoint for all sandwich helpers including BohrDomain 2-operator case
affects: [krylov_eigsolve, krylov_crossvalidation, thermalization, spectral-gap]

# Tech tracking
tech-stack:
  added: []
  patterns: [L*rho*L' Lindblad convention via kron(conj(J),J)]

key-files:
  created: []
  modified:
    - src/qi_tools.jl
    - src/krylov_matvec.jl
    - src/krylov_workspace.jl
    - test/reference/energy_dm_reference.bson
    - test/reference/trotter_coherent_dm_reference.bson

key-decisions:
  - "Use kron(conj(J), J) for L*rho*L' vectorization (Watrous convention)"
  - "Simplify adjoint G matrices: G_left_adj=G_right for all domains since R is always Hermitianized"
  - "Regenerate frozen BSON references to match new convention"

patterns-established:
  - "Lindblad sandwich: L*rho*L' via BLAS gemm!('N','N') + gemm!('N','C') for forward, gemm!('C','N') + gemm!('N','N') for adjoint"
  - "G_left = iB^T - 0.5*R (not R^T) -- direct R usage with new kron ordering"

# Metrics
duration: 26min
completed: 2026-02-26
---

# Quick Task 39: Fix Krylov Lindblad Operator Convention Summary

**Unified Lindblad dissipator to L*rho*L' convention across dense vectorization, Krylov matvec, and Thermalize paths, fixing incorrect conj(L)*rho*L^T that differed from jump_workers.jl**

## Performance

- **Duration:** 26 min
- **Started:** 2026-02-26T11:40:24Z
- **Completed:** 2026-02-26T12:06:29Z
- **Tasks:** 3
- **Files modified:** 5

## Accomplishments
- Dense Lindbladian vectorization (qi_tools.jl) now uses kron(conj(J), J) for J*rho*J' sandwich
- All 6 Krylov sandwich helpers updated: 4 single-operator + 2 two-operator variants
- G_left/G_right precomputation corrected to use R directly instead of R^T
- All 1198 tests pass including round-trip, adjoint duality, cross-validation, and eigsolve

## Task Commits

Each task was committed atomically:

1. **Task 1: Fix dense Lindbladian vectorization convention in qi_tools.jl** - `2c9236e` (fix)
2. **Task 2: Fix Krylov sandwich helpers and R_total accumulation** - `9c49f40` (fix)
3. **Task 2b: Fix G_left/G_right and 2op adjoint (deviation fix)** - `b4b1cec` (fix)

## Files Created/Modified
- `src/qi_tools.jl` - Swapped kron argument order in single-op and 2-op dissipator vectorization
- `src/krylov_matvec.jl` - Updated all 6 sandwich helpers: forward uses gemm('N','C'), adjoint uses gemm('C','N')
- `src/krylov_workspace.jl` - G_left/G_right use R (not R^T), simplified adjoint G for all domains
- `test/reference/energy_dm_reference.bson` - Regenerated for new convention
- `test/reference/trotter_coherent_dm_reference.bson` - Regenerated for new convention

## Decisions Made
- kron(conj(J), J) chosen over kron(J, conj(J)) to match the Watrous convention identity: kron(conj(J), J)*vec(rho) = vec(J*rho*J')
- Anticommutator terms use kron(I, J'J) and kron((J'J)^T, I) for -0.5*(J'J*rho + rho*J'J)
- Adjoint G simplified to swap(G_left, G_right) for ALL domains since hermitianize!(R_total) ensures R'=R

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] G_left/G_right used R^T instead of R in workspace precomputation**
- **Found during:** Task 2 verification (BohrDomain adjoint round-trip failures)
- **Issue:** Old convention vectorized anticommutator as (J'J)^T*rho + rho*(J'J)^T, requiring R^T in G matrices. New convention uses J'J*rho + rho*J'J, requiring R directly. The workspace still used R^T.
- **Fix:** Changed G_left = iB^T - 0.5*R (was iB^T - 0.5*R^T). Simplified adjoint G to always swap G_left/G_right since R is Hermitianized.
- **Files modified:** src/krylov_workspace.jl (both Lindbladian and Thermalize constructors)
- **Verification:** All 198 Krylov matvec tests pass
- **Committed in:** b4b1cec

**2. [Rule 1 - Bug] 2-operator adjoint sandwich computed B_dag'*rho*A' instead of A'*rho*B_dag'**
- **Found during:** Task 2 verification (BohrDomain adjoint round-trip failures)
- **Issue:** HS adjoint of f(X)=A*X*B is f*(Y)=A'*Y*B'. The implementation had the factors swapped (B_dag'*rho*A' instead of A'*rho*B_dag'). For the old convention this was correct because the old forward was B_dag^T*rho*A^T, whose adjoint is conj(B_dag)*rho*conj(A).
- **Fix:** Swapped BLAS.gemm! arguments in _accumulate_adjoint_sandwich_2op! to compute A'*rho*B_dag'
- **Files modified:** src/krylov_matvec.jl
- **Verification:** BohrDomain adjoint round-trip and duality tests pass
- **Committed in:** b4b1cec

**3. [Rule 1 - Bug] Frozen BSON regression references outdated after convention change**
- **Found during:** Task 3 (full test suite)
- **Issue:** DM regression test compared against frozen reference computed with old conj(L)*rho*L^T convention
- **Fix:** Regenerated both BSON reference files using generate_references.jl
- **Files modified:** test/reference/energy_dm_reference.bson, test/reference/trotter_coherent_dm_reference.bson
- **Verification:** DM regression test passes
- **Committed in:** b4b1cec

---

**Total deviations:** 3 auto-fixed (3 bugs)
**Impact on plan:** All fixes necessary for correctness. The plan correctly identified the sandwich helper changes but missed the cascading effect on G_left/G_right precomputation and the adjoint sandwich argument order. No scope creep.

## Issues Encountered
None beyond the deviations documented above.

## Next Phase Readiness
- All Lindblad paths now consistently use L*rho*L' convention
- Krylov spectral gap estimates will match Thermalize fixed point for complex jump operators
- No blockers for Phase 34 (Code Deduplication)

## Self-Check: PASSED

All 5 modified files exist. All 3 commit hashes verified. Summary file exists.

---
*Quick Task: 39*
*Completed: 2026-02-26*
