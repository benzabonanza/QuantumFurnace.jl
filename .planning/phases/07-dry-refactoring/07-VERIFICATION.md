---
phase: 07-dry-refactoring
verified: 2026-02-15T07:33:26Z
status: passed
score: 5/5 must-haves verified
re_verification: false
---

# Phase 7: DRY Refactoring Verification Report

**Phase Goal:** Repeated code patterns are extracted into single-source helpers -- Hermitianization, CPTP channel application, coherent unitary application, and Trotter basis transforms each have one canonical implementation

**Verified:** 2026-02-15T07:33:26Z
**Status:** passed
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| #   | Truth                                                                                                                       | Status     | Evidence                                                                                         |
| --- | --------------------------------------------------------------------------------------------------------------------------- | ---------- | ------------------------------------------------------------------------------------------------ |
| 1   | A single hermitianize! helper replaces all 8+ inline Hermitianization patterns                                             | ✓ VERIFIED | 13 call sites found; 0 inline `0.5 .* (X .+ X')` patterns remain                                |
| 2   | CPTP channel application (K0/residual/Cholesky sequence) exists as one shared function, called from 3 sites                | ✓ VERIFIED | apply_cptp_channel! defined once, 3 calls; delta_factor_for_K0 only in helper                   |
| 3   | Coherent unitary application exists as one shared function, called from 3 sites                                            | ✓ VERIFIED | apply_coherent_unitary! defined once, 3 calls                                                    |
| 4   | Trotter basis transform of jumps has one canonical implementation shared between furnace.jl and trajectories.jl            | ✓ VERIFIED | transform_jumps_to_basis defined once, 3 calls (2 in furnace.jl, 1 in trajectories.jl)          |
| 5   | All 224 existing tests pass with no regressions                                                                            | ✓ VERIFIED | Test suite passed: 224/224 tests pass                                                           |

**Score:** 5/5 truths verified

### Required Artifacts

| Artifact                                | Expected                                                               | Status     | Details                                                                         |
| --------------------------------------- | ---------------------------------------------------------------------- | ---------- | ------------------------------------------------------------------------------- |
| `src/qi_tools.jl`                       | hermitianize! in-place helper                                          | ✓ VERIFIED | Function defined at line 12, substantive (3 lines), called from 13 sites       |
| `src/qi_tools.jl`                       | transform_jumps_to_basis helper                                        | ✓ VERIFIED | Function defined at line 23, substantive (1 line comprehension), 3 call sites  |
| `src/jump_workers.jl`                   | apply_coherent_unitary! helper                                         | ✓ VERIFIED | Function defined at line 148, substantive (@inline, 9 lines), 3 call sites     |
| `src/jump_workers.jl`                   | apply_cptp_channel! helper                                             | ✓ VERIFIED | Function defined at line 174, substantive (42 lines), 3 call sites             |
| `src/jump_workers.jl`                   | jump_contribution! methods calling helpers instead of inline patterns  | ✓ VERIFIED | 3 algorithmic methods (Bohr/Energy/Time-Trotter) use both helpers              |
| `src/trajectories.jl`                   | precompute_R and run_trajectories calling hermitianize!               | ✓ VERIFIED | 7 hermitianize! calls found, 1 transform_jumps_to_basis call                   |
| `src/furnace.jl`                        | construct_lindbladian and run_thermalization calling helpers           | ✓ VERIFIED | 1 hermitianize! call, 2 transform_jumps_to_basis calls                         |
| `src/log_sobolev.jl`                    | calling hermitianize!                                                  | ✓ VERIFIED | 1 hermitianize! call at line 176                                               |

### Key Link Verification

| From                                                                  | To                                   | Via                          | Status  | Details                                                      |
| --------------------------------------------------------------------- | ------------------------------------ | ---------------------------- | ------- | ------------------------------------------------------------ |
| `src/jump_workers.jl`                                                 | `src/qi_tools.jl`                    | hermitianize! call           | ✓ WIRED | 4 call sites (lines 213, 309, 387, 449)                      |
| `src/furnace.jl`                                                      | `src/qi_tools.jl`                    | transform_jumps_to_basis     | ✓ WIRED | 2 call sites (lines 75, 118)                                 |
| `src/trajectories.jl`                                                 | `src/qi_tools.jl`                    | hermitianize! call           | ✓ WIRED | 7 call sites in precompute_R and run_trajectories           |
| `src/trajectories.jl`                                                 | `src/qi_tools.jl`                    | transform_jumps_to_basis     | ✓ WIRED | 1 call site (line 70)                                        |
| `jump_contribution!(::BohrDomain, ...)`                               | `apply_cptp_channel!`                | function call                | ✓ WIRED | Line 311 in jump_workers.jl                                  |
| `jump_contribution!(::EnergyDomain, ...)`                             | `apply_cptp_channel!`                | function call                | ✓ WIRED | Line 389 in jump_workers.jl                                  |
| `jump_contribution!(::Union{TimeDomain,TrotterDomain}, ...)`          | `apply_cptp_channel!`                | function call                | ✓ WIRED | Line 451 in jump_workers.jl                                  |
| `jump_contribution!(::BohrDomain, ...)`                               | `apply_coherent_unitary!`            | function call                | ✓ WIRED | Line 240 in jump_workers.jl                                  |
| `jump_contribution!(::EnergyDomain, ...)`                             | `apply_coherent_unitary!`            | function call                | ✓ WIRED | Line 332 in jump_workers.jl                                  |
| `jump_contribution!(::Union{TimeDomain,TrotterDomain}, ...)`          | `apply_coherent_unitary!`            | function call                | ✓ WIRED | Line 410 in jump_workers.jl                                  |

### Requirements Coverage

**From ROADMAP.md Phase 7 Requirements:** DRY-01, DRY-02, DRY-03, DRY-04

| Requirement | Description                                                      | Status      | Blocking Issue |
| ----------- | ---------------------------------------------------------------- | ----------- | -------------- |
| DRY-01      | Single hermitianize! helper replaces all inline patterns         | ✓ SATISFIED | None           |
| DRY-02      | CPTP channel application as one shared function                  | ✓ SATISFIED | None           |
| DRY-03      | Coherent unitary application as one shared function              | ✓ SATISFIED | None           |
| DRY-04      | Trotter basis transform has one canonical implementation         | ✓ SATISFIED | None           |

### Anti-Patterns Found

| File                   | Line | Pattern                              | Severity | Impact                                                 |
| ---------------------- | ---- | ------------------------------------ | -------- | ------------------------------------------------------ |
| jump_workers.jl        | 455  | TODO comment (test BLAS threading)   | ℹ️ Info  | Pre-existing TODO, not related to this phase           |
| log_sobolev.jl         | 17   | TODO (rewrite with apply_lindbladian)| ℹ️ Info  | Pre-existing TODO, outside scope of this phase         |

**No blocker anti-patterns found.** The two TODO comments are pre-existing and unrelated to the DRY refactoring work.

### Human Verification Required

None. All verification criteria are programmatically testable and have been verified.

### Verification Details

**Plan 07-01 (hermitianize! and transform_jumps_to_basis):**

✓ hermitianize! function exists in qi_tools.jl (line 12)
✓ transform_jumps_to_basis function exists in qi_tools.jl (line 23)
✓ 13 call sites to hermitianize! across 4 files:
  - jump_workers.jl: 4 sites (1 in apply_cptp_channel!, 3 pre-CPTP)
  - trajectories.jl: 7 sites
  - furnace.jl: 1 site
  - log_sobolev.jl: 1 site
✓ 3 call sites to transform_jumps_to_basis:
  - furnace.jl: 2 sites (construct_lindbladian, run_thermalization)
  - trajectories.jl: 1 site (build_trajectoryframework)
✓ No inline `0.5 .* (X .+ X')` patterns remain
✓ No inline `eigvecs' * j.data * eigvecs` comprehensions remain
✓ All 224 tests pass

**Plan 07-02 (apply_coherent_unitary! and apply_cptp_channel!):**

✓ apply_coherent_unitary! function exists in jump_workers.jl (line 148)
✓ apply_cptp_channel! function exists in jump_workers.jl (line 174)
✓ 3 call sites to apply_coherent_unitary! (Bohr, Energy, Time/Trotter domains)
✓ 3 call sites to apply_cptp_channel! (Bohr, Energy, Time/Trotter domains)
✓ delta_factor_for_K0 only appears in apply_cptp_channel! helper (4 times, lines 182, 184, 191, 192)
✓ cholesky! only appears in apply_cptp_channel! helper (1 time, line 201)
✓ No inline U_B application blocks remain
✓ No inline K0/Cholesky/residual blocks remain
✓ All 224 tests pass

**Commit verification:**

✓ Commit c8c4478: feat(07-01): extract hermitianize! and transform_jumps_to_basis helpers
✓ Commit eaf5486: feat(07-02): extract apply_coherent_unitary! helper from 3 identical blocks
✓ Commit bd95b34: feat(07-02): extract apply_cptp_channel! helper from 3 identical blocks
✓ All commits present in git log

---

_Verified: 2026-02-15T07:33:26Z_
_Verifier: Claude (gsd-verifier)_
