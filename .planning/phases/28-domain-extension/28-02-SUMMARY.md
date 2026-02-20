---
phase: 28-domain-extension
plan: 02
subsystem: krylov-matvec
tags: [lindbladian, matvec, bohr-domain, two-operator-dissipator, kron-convention]

# Dependency graph
requires:
  - phase: 27-core-matvec-infrastructure
    provides: "EnergyDomain apply_lindbladian!/apply_adjoint_lindbladian!, KrylovWorkspace, dissipator helpers"
  - phase: 28-domain-extension
    plan: 01
    provides: "TimeDomain/TrotterDomain matvec methods (structural template)"
provides:
  - "apply_lindbladian! for BohrDomain with bucket iteration and 2-operator dissipator"
  - "apply_adjoint_lindbladian! for BohrDomain with dedicated adjoint 2op helper"
  - "_accumulate_dissipator_2op! and _accumulate_adjoint_dissipator_2op! helpers"
  - "Round-trip and duality tests for BohrDomain (KMS, GNS)"
affects: [29-krylov-integration, 30-gap-estimation]

# Tech tracking
tech-stack:
  added: []
  patterns: ["kron-derived 2op dissipator formula matching dense vectorization convention", "dense scatter for A_nu2_dag (one alloc per matvec, zero sparse())"]

key-files:
  created: []
  modified:
    - src/krylov_matvec.jl
    - test/test_krylov_matvec.jl

key-decisions:
  - "2op dissipator formula derived from kron vectorization convention (B_dag'*rho*A' sandwich), not physics convention (A*rho*B_dag')"
  - "Dedicated _accumulate_adjoint_dissipator_2op! with all-N BLAS flags (no transpose), derived from HS adjoint of kron superoperator"
  - "Single dense A_nu2_dag buffer allocated per matvec call (zeros+scatter), reused across bucket iterations"

patterns-established:
  - "kron(P,Q) vec(X) = vec(Q*X*P^T) un-vectorization for matching dense Lindbladian code"
  - "HS adjoint of kron-based superoperator: sandwich X*rho*Y -> X^H*rho*Y^H, anticomm M^T*rho -> conj(M)*rho"

# Metrics
duration: 19min
completed: 2026-02-20
---

# Phase 28 Plan 02: BohrDomain Matvec Summary

**Matrix-free BohrDomain Lindbladian action with kron-derived two-operator dissipator, bucket iteration over Bohr frequencies, and round-trip verification < 1e-12 against dense construction**

## Performance

- **Duration:** 19 min
- **Started:** 2026-02-20T14:23:32Z
- **Completed:** 2026-02-20T14:43:12Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments
- `apply_lindbladian!` and `apply_adjoint_lindbladian!` for `BohrDomain` -- bucket iteration over `hamiltonian.bohr_dict` keys with entrywise alpha computation
- `_accumulate_dissipator_2op!` and `_accumulate_adjoint_dissipator_2op!` helpers implementing the kron-derived two-operator dissipator formula
- 4 new testsets: BohrDomain KMS forward, GNS forward, adjoint, and duality check
- All round-trip errors < 1e-12 across 10 random density matrices at n=4
- Duality check |tr(X' * L(Y)) - tr(L*(X)' * Y)| < 1e-11 for 5 random pairs
- Full parity: all four domains (Energy, Time, Trotter, Bohr) now have both forward and adjoint matvec

## Task Commits

Each task was committed atomically:

1. **Task 1: Implement 2op dissipator helpers and BohrDomain forward/adjoint matvec** - `eb95618` (feat)
2. **Task 2: Add round-trip and duality tests for BohrDomain** - `5d370f7` (test + fix)

## Files Created/Modified
- `src/krylov_matvec.jl` - Added 4 new functions: `_accumulate_dissipator_2op!` (kron-derived forward 2op), `_accumulate_adjoint_dissipator_2op!` (HS adjoint of kron superoperator), `apply_lindbladian!` for BohrDomain (bucket iteration + dense scatter A_nu2_dag), `apply_adjoint_lindbladian!` for BohrDomain (sign-flipped coherent + adjoint 2op helper)
- `test/test_krylov_matvec.jl` - Added 4 new testsets (BohrDomain KMS forward, GNS forward, adjoint, duality check) totaling 35 new test assertions

## Decisions Made
- **Kron-derived formula instead of physics formula:** The plan specified the physics-convention dissipator `D(A,B)(rho) = A*rho*B' - 0.5*(B*A*rho + rho*B*A)`, but the dense code's kron vectorization convention produces `D_kron(j1,j2)(rho) = j2^T*rho*j1^T - 0.5*((j2*j1)^T*rho + rho*(j2*j1)^T)`. These differ because `kron(P,Q) vec(X) = vec(Q*X*P^T)` introduces transposes. The corrected forward formula is `B_dag'*rho*A' - 0.5*((B_dag*A)'*rho + rho*(B_dag*A)')`.
- **Dedicated adjoint helper with all-'N' BLAS flags:** The HS adjoint of the kron-based superoperator maps `X*rho*Y -> X^H*rho*Y^H` and `M^T*rho -> conj(M)*rho`. For real-valued operators (BohrDomain always produces real eigenbasis and real alpha values), this simplifies to `B_dag*rho*A - 0.5*(B_dag*A*rho + rho*B_dag*A)` with no transpose/adjoint flags.
- **Single dense buffer for A_nu2_dag:** Allocated once per matvec call (`zeros(T, dim, dim)`) and reused across bucket iterations via `fill!` + scatter. This avoids per-bucket `sparse()` allocation while keeping the code BLAS-friendly.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed 2op dissipator formula to match kron vectorization convention**
- **Found during:** Task 2 (round-trip test verification)
- **Issue:** The plan's formula `D(A,B)(rho) = A*rho*B' - 0.5*(B*A*rho + rho*B*A)` uses the physics convention, but the dense code vectorizes via kron with `(P kron Q) vec(X) = vec(Q*X*P^T)`, producing transposed operators. The forward round-trip error was ~0.07 (should be < 1e-12).
- **Fix:** Rewrote both `_accumulate_dissipator_2op!` and `_accumulate_adjoint_dissipator_2op!` to use the kron-derived formulas. Forward: `B_dag'*rho*A' - 0.5*((B_dag*A)'*rho + rho*(B_dag*A)')`. Adjoint (real): `B_dag*rho*A - 0.5*(B_dag*A*rho + rho*B_dag*A)`.
- **Files modified:** src/krylov_matvec.jl
- **Verification:** All 1095 tests pass including BohrDomain round-trip < 1e-12 and duality < 1e-11
- **Committed in:** 5d370f7 (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (1 bug)
**Impact on plan:** Formula correction necessary for mathematical correctness. The plan's physics-convention formula was wrong for this codebase's kron-based vectorization. No scope creep.

## Issues Encountered
None beyond the formula mismatch documented above. The root cause is that the kron vectorization identity `kron(P,Q) vec(X) = vec(Q*X*P^T)` introduces transposes that the physics-convention formula `D(A,B)(rho) = A*rho*B'` does not account for. This is a well-known subtlety when converting between vectorized (superoperator) and un-vectorized (matrix-on-matrix) Lindbladian representations.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Phase 28 complete: all four domains (Energy, Time, Trotter, Bohr) have `apply_lindbladian!` and `apply_adjoint_lindbladian!`
- Ready for Phase 29 (Krylov Integration) which will use these matvec methods as the `linmap` argument to `KrylovKit.eigsolve`
- All 1095 tests pass

## Self-Check: PASSED

- [x] src/krylov_matvec.jl exists
- [x] test/test_krylov_matvec.jl exists
- [x] 28-02-SUMMARY.md exists
- [x] Commit eb95618 found
- [x] Commit 5d370f7 found

---
*Phase: 28-domain-extension*
*Completed: 2026-02-20*
