---
phase: 27-core-matvec-infrastructure
plan: 01
subsystem: krylov
tags: [lindbladian, matvec, workspace, energy-domain, krylov, zero-allocation]

# Dependency graph
requires:
  - phase: existing codebase
    provides: "_precompute_data, _precompute_coherent_total_B, oft!, JumpOp, AbstractLiouvConfig, HamHam"
provides:
  - "KrylovWorkspace struct with pre-allocated scratch matrices"
  - "apply_lindbladian!(ws, rho, config, ham) for EnergyDomain"
  - "apply_adjoint_lindbladian!(ws, rho, config, ham) for EnergyDomain"
  - "_accumulate_dissipator!(out, L_op, rho, scalar, ws) helper"
affects: [27-02-PLAN (tests), 28 (other-domain dispatch), 29 (KrylovKit integration)]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "KrylovWorkspace pre-allocation pattern (5 scratch matrices tied to config+ham)"
    - "Matrix-free Lindbladian matvec via _accumulate_dissipator! helper"
    - "Adjoint via L<->L' swap in dissipator + coherent sign flip"

key-files:
  created:
    - "src/krylov_workspace.jl"
    - "src/krylov_matvec.jl"
  modified: []

key-decisions:
  - "Two type parameters KrylovWorkspace{T,PD} for zero-overhead NamedTuple access"
  - "3-arg mul! + @. broadcasting pattern (not 5-arg mul! with alpha/beta) for clarity and guaranteed zero alloc"
  - "Lazy adjoint ws.jump_oft' passed directly to _accumulate_dissipator! for negative freq and adjoint swap"

patterns-established:
  - "KrylovWorkspace pattern: pre-allocate at construction, pass ws to all in-place functions"
  - "_accumulate_dissipator! as shared building block for forward and adjoint"
  - "Domain dispatch via config type parameter (Phase 28 adds more methods)"

# Metrics
duration: 3min
completed: 2026-02-20
---

# Phase 27 Plan 01: KrylovWorkspace and EnergyDomain Matvec Summary

**KrylovWorkspace struct with 5 scratch matrices and matrix-free apply_lindbladian!/apply_adjoint_lindbladian! for EnergyDomain, mirroring dense _jump_contribution! loop**

## Performance

- **Duration:** 3 min
- **Started:** 2026-02-20T11:28:26Z
- **Completed:** 2026-02-20T11:31:19Z
- **Tasks:** 2
- **Files created:** 2

## Accomplishments
- KrylovWorkspace struct with dual type parameters (T for Complex element type, PD for NamedTuple precomputed data) ensuring zero-overhead field access
- Constructor that mirrors construct_lindbladian setup: determines ham_or_trott, calls _precompute_data and _precompute_coherent_total_B, allocates 5 scratch matrices
- _accumulate_dissipator! helper factoring out the Lindblad dissipator formula D_L(rho) using 3-arg mul! chains
- apply_lindbladian! for EnergyDomain that mirrors _jump_contribution! line-by-line: coherent -i[B,rho], Hermitian half-grid optimization, oft! projection, prefactor computation
- apply_adjoint_lindbladian! with sign-flipped coherent +i[B,rho] and swapped L<->L' in dissipator calls

## Task Commits

Each task was committed atomically:

1. **Task 1: Create KrylovWorkspace struct and constructor** - `ff0a98f` (feat)
2. **Task 2: Implement apply_lindbladian!, apply_adjoint_lindbladian!, _accumulate_dissipator!** - `3206b1b` (feat)

## Files Created/Modified
- `src/krylov_workspace.jl` - KrylovWorkspace{T,PD} struct and constructor; pre-allocates 5 scratch matrices and stores precomputed data, B_total, and jump references
- `src/krylov_matvec.jl` - _accumulate_dissipator! helper, apply_lindbladian! for EnergyDomain (forward), apply_adjoint_lindbladian! for EnergyDomain (adjoint)

## Decisions Made
- Used two type parameters `KrylovWorkspace{T<:Complex, PD<:NamedTuple}` rather than `Any` for precomputed_data, ensuring Julia specializes on the concrete NamedTuple type for zero-overhead field access
- Used 3-arg `mul!` + `@.` broadcasting accumulation pattern rather than 5-arg `mul!` with complex alpha/beta for clarity and guaranteed allocation-free behavior with Adjoint wrappers
- Passed lazy `ws.jump_oft'` directly to `_accumulate_dissipator!` as L_op argument rather than materializing explicit adjoints, leveraging Julia's zero-cost Adjoint wrapper

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
- Git index corruption during Task 2 commit required index rebuild via `rm .git/index && git read-tree HEAD`. No work lost; commit succeeded on retry.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- Both source files ready for round-trip testing in plan 02 (matvec vs dense Liouvillian comparison)
- Files are NOT yet included in the QuantumFurnace.jl module (plan 02 or later will add the includes)
- Phase 28 can add Time/Trotter/Bohr domain dispatch methods without workspace restructuring
- Phase 29 can wrap apply_lindbladian! in KrylovKit eigsolve closure

## Self-Check: PASSED

- [x] `src/krylov_workspace.jl` exists with KrylovWorkspace struct (8 fields: precomputed_data, B_total, jumps, jump_oft, tmp1, tmp2, LdagL, rho_out)
- [x] `src/krylov_matvec.jl` exists with 3 functions: _accumulate_dissipator!, apply_lindbladian!, apply_adjoint_lindbladian!
- [x] Commit ff0a98f exists (Task 1)
- [x] Commit 3206b1b exists (Task 2)
- [x] Both files parse without errors (verified via Meta.parseall)

---
*Phase: 27-core-matvec-infrastructure*
*Completed: 2026-02-20*
