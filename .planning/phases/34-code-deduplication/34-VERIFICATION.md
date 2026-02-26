---
phase: 34-code-deduplication
verified: 2026-02-26T13:45:00Z
status: passed
score: 9/9 must-haves verified
---

# Phase 34: Code Deduplication Verification Report

**Phase Goal:** The 16+ copy-pasted prefactor formulas, hermitian half-grid branching patterns, and OFT variants are each single-source functions dispatched on domain type
**Verified:** 2026-02-26T13:45:00Z
**Status:** passed
**Re-verification:** No -- initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | domain_prefactor(::EnergyDomain, w0, sigma) returns w0/(sigma*sqrt(2*pi)) and is called from all 5 simulation files instead of inline formula | VERIFIED | Defined at furnace_utensils.jl:13. Used in jump_workers.jl (4 sites), krylov_workspace.jl (2 sites), krylov_matvec.jl (4 sites), krylov_eigsolve.jl (2 sites), trajectories.jl (3 sites). Zero inline formulas remain outside furnace_utensils.jl |
| 2 | domain_prefactor(::TimeDomain, w0, sigma, t0) and domain_prefactor(::TrotterDomain, w0, sigma, t0) return identical time-domain formula | VERIFIED | Both defined at furnace_utensils.jl:14-15 with identical body: `w0 * t0^2 * (sigma * sqrt(2 / pi)) / (2 * pi)` |
| 3 | No domain_prefactor method exists for BohrDomain -- BohrDomain callers use gamma_norm_factor directly | VERIFIED | Only 3 methods defined (Energy, Time, Trotter). BohrDomain _precompute_data (lines 31-80) does not include domain_prefactor field |
| 4 | A single oft!(out, eigenbasis, bohr_freqs, energy, inv_4sigma2) replaces both old JumpOp-based oft! (deleted) and _krylov_oft! (deleted) | VERIFIED | New signature at ofts.jl:10-19. grep for `_krylov_oft!` returns 0 results in src/. grep for `JumpOp` in ofts.jl returns only time_oft!/trotter_oft! parameter types (lines 25, 77) which are retained test utilities |
| 5 | time_oft! and trotter_oft! remain unchanged as test/debug utilities in ofts.jl | VERIFIED | time_oft! at ofts.jl:22-72, trotter_oft! at ofts.jl:74-124. Both intact with OFTCaches signatures |
| 6 | 4 sandwich helpers consolidated to 2: _accumulate_sandwich! (L*rho*L') and _accumulate_sandwich_adj! (L'*rho*L) | VERIFIED | Only 2 single-operator sandwich helpers exist in krylov_matvec.jl: _accumulate_sandwich! (lines 17-30) and _accumulate_sandwich_adj! (lines 41-54). grep for `_accumulate_adjoint_sandwich!` and `_accumulate_adjoint_sandwich_adj_L!` returns 0 results |
| 7 | _build_cptp_channel is called from krylov_workspace.jl, trajectories.jl, and jump_workers.jl instead of 3 independent CPTP constructions | VERIFIED | Defined in furnace_utensils.jl:184-201. Called at krylov_workspace.jl:386, trajectories.jl:160, jump_workers.jl:178. No inline `(2*alpha - delta)*R` formulas remain in caller files (only in comments/docstrings) |
| 8 | BohrDomain 2-operator sandwiches left unchanged | VERIFIED | `_accumulate_sandwich_2op!` (krylov_matvec.jl:204) and `_accumulate_adjoint_sandwich_2op!` (krylov_matvec.jl:229) both present and used in BohrDomain apply_lindbladian!/apply_adjoint_lindbladian! |
| 9 | trajectories.jl passes the 1/p_jump-rescaled per-operator R_a to _build_cptp_channel (not R_total), preserving per-operator Lie-Trotter semantics | VERIFIED | trajectories.jl:157 applies `R_a .*= (1.0 / p_jump)` then passes `R_a` to `_build_cptp_channel(R_a, delta)` at line 160 |

**Score:** 9/9 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `src/furnace_utensils.jl` | domain_prefactor definitions + _build_cptp_channel | VERIFIED | 3 domain_prefactor methods (lines 13-15), stored in precomputed_data for Energy (line 90) and Time/Trotter (line 131). _build_cptp_channel (lines 184-201) |
| `src/ofts.jl` | Unified oft! with concrete-typed signature | VERIFIED | Single `oft!(out, eigenbasis, bohr_freqs, energy, inv_4sigma2)` at lines 10-19. No JumpOp-based oft! remains. time_oft! and trotter_oft! preserved |
| `src/krylov_matvec.jl` | Consolidated sandwich helpers (2 instead of 4) | VERIFIED | Only `_accumulate_sandwich!` and `_accumulate_sandwich_adj!` defined, both @inline. 2op helpers preserved for BohrDomain |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| furnace_utensils.jl (domain_prefactor) | jump_workers.jl, krylov_workspace.jl, krylov_matvec.jl, krylov_eigsolve.jl, trajectories.jl | precomputed_data.domain_prefactor | WIRED | 15 usage sites across all 5 files confirmed via grep |
| ofts.jl (unified oft!) | jump_workers.jl, krylov_workspace.jl, krylov_matvec.jl, krylov_eigsolve.jl, trajectories.jl | oft!(out, eigenbasis, bohr_freqs, energy, inv_4sigma2) | WIRED | All callers use new concrete-typed signature. Zero _krylov_oft! references remain |
| krylov_matvec.jl (_accumulate_sandwich!, _accumulate_sandwich_adj!) | krylov_matvec.jl (Energy/Time/Trotter forward+adjoint) | direct call | WIRED | 16 call sites across all domain/direction combinations |
| furnace_utensils.jl (_build_cptp_channel) | krylov_workspace.jl, trajectories.jl, jump_workers.jl | function call | WIRED | 3 callers confirmed, no inline CPTP formulas remain |

### Requirements Coverage

| Requirement | Status | Blocking Issue |
|-------------|--------|----------------|
| DEDUP-01: Extract domain_prefactor() replacing 16 formula copies | SATISFIED | None |
| DEDUP-02: Extract foreach_frequency() iterator | DEFERRED | User decision: keep explicit for-loop patterns |
| DEDUP-03: Unify oft!() with _krylov_oft! | SATISFIED | None |

### Anti-Patterns Found

No anti-patterns found. All modified files are clean of TODO/FIXME/HACK/PLACEHOLDER markers. No empty implementations or stub patterns detected.

### Commit Verification

All 4 task commits verified in git history:
- `02cb3ea` feat(34-01): extract domain_prefactor and store in precomputed_data
- `61871fb` feat(34-01): unify oft! and _krylov_oft! into single function
- `a1b5363` feat(34-02): consolidate 4 sandwich helpers into 2 in krylov_matvec.jl
- `aab8d8b` feat(34-02): extract _build_cptp_channel shared helper

### Human Verification Required

### 1. Test Suite Passes

**Test:** Run `julia --project -e 'using Pkg; Pkg.test()'`
**Expected:** All 1198 tests pass with zero numerical regressions
**Why human:** Cannot run Julia test suite in verification environment; SUMMARY claims all tests pass but this needs runtime confirmation

### 2. Zero-Allocation Krylov Hot Path

**Test:** Run allocation tests targeting the Krylov matvec hot path
**Expected:** 0 bytes allocated in sandwich helper calls (the @inline annotation must be effective)
**Why human:** Allocation behavior depends on Julia compiler inlining decisions at runtime

---

_Verified: 2026-02-26T13:45:00Z_
_Verifier: Claude (gsd-verifier)_
