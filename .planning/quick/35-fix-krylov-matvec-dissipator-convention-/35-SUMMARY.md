---
phase: quick-35
plan: 01
subsystem: krylov-matvec
tags: [BLAS, dissipator, kron-convention, complex-operators, Lindbladian]

# Dependency graph
requires:
  - phase: 27-krylov-matvec
    provides: "6 dissipator accumulation functions and apply_lindbladian!/apply_adjoint_lindbladian!"
  - phase: 28-domain-extension
    provides: "TimeDomain, TrotterDomain, BohrDomain matvec implementations"
provides:
  - "Krylov matvec correct for general complex non-Hermitian jump operators"
  - "All sandwich, anticommutator, and coherent terms match dense kron convention exactly"
  - "Complex jump operator round-trip tests for EnergyDomain and TimeDomain"
affects: [29-eigensolver-integration, 30-cross-validation, 31-production-pipeline]

# Tech tracking
tech-stack:
  added: []
  patterns: ["kron(A,B) vec(X) = vec(B*X*A^T) convention consistently applied"]

key-files:
  modified:
    - "src/krylov_matvec.jl"
    - "test/test_krylov_matvec.jl"

key-decisions:
  - "Dense kron convention kron(A,B)vec(X)=vec(B*X*A^T) is the ground truth for all Krylov matvec formulas"
  - "Anticommutator uses (L'L)^T not L'L -- matches kron(L'L, I) un-vectorization"
  - "Coherent term uses i[B^T, rho] not -i[B, rho] -- matches kron(B, I) + kron(I, B^T) un-vectorization"

patterns-established:
  - "All BLAS operations in dissipator helpers use transpose ('T') flags, never adjoint ('C'), for kron-convention correctness"
  - "Coherent term uses BLAS.gemm! with 'T' flag on B instead of mul!(tmp, B, rho)"

# Metrics
duration: 34min
completed: 2026-02-20
---

# Quick-35: Fix Krylov Matvec Dissipator Convention Summary

**Fixed all Krylov matvec sandwich, anticommutator, and coherent terms to match dense kron(A,B)vec(X)=vec(B*X*A^T) convention for complex non-Hermitian jump operators**

## Performance

- **Duration:** 34 min
- **Started:** 2026-02-20T15:09:38Z
- **Completed:** 2026-02-20T15:43:47Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments
- Fixed sandwich terms in all 6 dissipator functions: `conj(L)*rho*L^T` instead of `L*rho*L'` (and corresponding adjoint/2-op variants)
- Fixed anticommutator terms to use `(L'L)^T` via BLAS 'T' flag instead of `L'L` directly
- Fixed coherent terms in all 6 apply functions to use `i[B^T, rho]` matching the dense vectorization convention
- Added 4 new testsets (45 individual tests) with complex non-Hermitian jump operators validating EnergyDomain and TimeDomain round-trips
- All 1140 tests pass (1095 existing + 45 new)

## Task Commits

Each task was committed atomically:

1. **Task 1: Fix all 6 dissipator helper functions** - `7e8fb99` (fix) - sandwich terms
2. **Task 1 deviation: Fix anticommutator + coherent** - `d45b88b` (fix) - anticommutator transpose and coherent B^T
3. **Task 2: Add complex jump operator round-trip test** - `cd249c3` (test)

## Files Created/Modified
- `src/krylov_matvec.jl` - Fixed sandwich, anticommutator, and coherent terms in all 6 dissipator helpers and 6 apply functions
- `test/test_krylov_matvec.jl` - Added 4 testsets with complex non-Hermitian jump operator (EnergyDomain forward/adjoint/duality + TimeDomain forward+adjoint)

## Decisions Made
- Dense kron convention `kron(A,B) vec(X) = vec(B * X * A^T)` is the single source of truth. All Krylov formulas are derived from un-vectorizing the dense kron expressions.
- Anticommutator uses `(L'L)^T` not `L'L`: the dense code uses `kron(L'L, I)` which un-vectorizes to `rho * (L'L)^T` and `kron(I, (L'L)^T)` which un-vectorizes to `(L'L)^T * rho`.
- Coherent term uses `i[B^T, rho]` not `-i[B, rho]`: the dense code uses `kron(B, I, -i) + kron(I, B^T, +i)` which un-vectorizes to `i*B^T*rho - i*rho*B^T`.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Anticommutator terms used L'L instead of (L'L)^T**
- **Found during:** Task 2 (complex jump round-trip testing)
- **Issue:** Dense code uses `kron(L'L, I)` which un-vectorizes to `rho * (L'L)^T`, not `rho * L'L`. For complex-Hermitian `L'L`, `(L'L)^T = conj(L'L) != L'L`. Bug was masked by real-valued Pauli test operators.
- **Fix:** Changed all 4 single-op dissipator functions to use BLAS 'T' flag on ws.LdagL. Changed 2-op functions (already correct -- plan correctly specified 'T' flags).
- **Files modified:** src/krylov_matvec.jl
- **Verification:** Complex jump round-trip error dropped from ~0.19 to ~0.03 (remaining error was coherent term)
- **Committed in:** d45b88b

**2. [Rule 1 - Bug] Coherent term used -i[B, rho] instead of i[B^T, rho]**
- **Found during:** Task 2 (complex jump round-trip testing)
- **Issue:** Dense coherent vectorization uses `kron(B, I, -i) + kron(I, B^T, +i)` which un-vectorizes to `i*B^T*rho - i*rho*B^T = i[B^T, rho]`. Krylov code used `-i[B, rho]`. For Hermitian B with complex entries, `B^T = conj(B) != B`.
- **Fix:** Replaced all 6 coherent blocks (3 forward + 3 adjoint across EnergyDomain, BohrDomain, TimeDomain/TrotterDomain) to use `BLAS.gemm!('T', ...)` on B instead of `mul!(..., B, ...)`.
- **Files modified:** src/krylov_matvec.jl
- **Verification:** All complex jump round-trip errors dropped to machine precision (~1e-16)
- **Committed in:** d45b88b

---

**Total deviations:** 2 auto-fixed (2 bugs)
**Impact on plan:** Both auto-fixes were essential for correctness with complex operators. The plan correctly identified the sandwich term fix but did not identify the anticommutator transpose and coherent term transpose issues (which use the same kron un-vectorization convention). No scope creep.

## Issues Encountered
- `JumpOp{Matrix{ComplexF64}}` vs `JumpOp` type mismatch: `construct_lindbladian` expects `Vector{JumpOp}` but `[complex_jump]` infers parameterized type. Fixed by using `JumpOp[complex_jump]` syntax.
- Git index corruption required `rm .git/index && git reset` to recover.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Krylov matvec now correct for arbitrary complex jump operators
- Ready for Phase 29 (Eigensolver Integration) which will use the matvec as a linear operator

---
*Quick task: 35-fix-krylov-matvec-dissipator-convention*
*Completed: 2026-02-20*
