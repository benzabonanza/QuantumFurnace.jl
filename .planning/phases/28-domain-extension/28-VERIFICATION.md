---
phase: 28-domain-extension
verified: 2026-02-20T14:49:37Z
status: passed
score: 3/3 must-haves verified
---

# Phase 28: Domain Extension Verification Report

**Phase Goal:** Matrix-free Lindbladian action works for all four domains, with proper NUFFT prefactors (Time), Trotter eigenbasis (Trotter), and Bohr bucket iteration (Bohr)
**Verified:** 2026-02-20T14:49:37Z
**Status:** passed
**Re-verification:** No -- initial verification

## Goal Achievement

### Observable Truths

| #   | Truth                                                                                                          | Status     | Evidence                                                                                                                                 |
| --- | -------------------------------------------------------------------------------------------------------------- | ---------- | ---------------------------------------------------------------------------------------------------------------------------------------- |
| 1   | `apply_lindbladian!` for TimeDomain passes round-trip test against dense at n=4 (KMS + GNS), using NUFFT prefactors | VERIFIED | Test testsets 8-9 in test_krylov_matvec.jl; `_prefactor_view(oft_nufft_prefactors, w)` used in `Union{TimeDomain,TrotterDomain}` dispatch; `< 1e-12` tolerance; all 1095 tests pass |
| 2   | `apply_lindbladian!` for TrotterDomain passes round-trip test at n=4 (KMS + GNS), operating in Trotter eigenbasis (trotter.eigvecs) | VERIFIED | Test testsets 12-13; TEST_TROTTER_JUMPS built from `trotter.eigvecs` basis; `trotter=TEST_TROTTER` passed to both `construct_lindbladian` and `KrylovWorkspace`; `< 1e-12` tolerance; all 1095 tests pass |
| 3   | `apply_lindbladian!` for BohrDomain passes round-trip test at n=4 (KMS + GNS), using Bohr frequency bucket iteration with generalized two-operator dissipator | VERIFIED | Test testsets 16-19; `for nu_2 in keys(hamiltonian.bohr_dict)` bucket iteration; `_accumulate_dissipator_2op!` and `_accumulate_adjoint_dissipator_2op!` helpers; duality check `< 1e-11`; all 1095 tests pass |

**Score:** 3/3 truths verified

### Required Artifacts

| Artifact                      | Expected                                                         | Status    | Details                                                                                                                            |
| ----------------------------- | ---------------------------------------------------------------- | --------- | ---------------------------------------------------------------------------------------------------------------------------------- |
| `src/krylov_matvec.jl`        | apply_lindbladian! and apply_adjoint_lindbladian! for Time/Trotter/Bohr | VERIFIED | Contains `Union{TimeDomain, TrotterDomain}` dispatch (lines 585-717), `BohrDomain` dispatch (lines 455-557), and four 2op helper functions (lines 352-428) |
| `test/test_krylov_matvec.jl`  | Round-trip and allocation tests for TimeDomain, TrotterDomain, and BohrDomain | VERIFIED | 19 testsets total: testsets 8-11 (TimeDomain), 12-15 (TrotterDomain), 16-19 (BohrDomain) |

### Key Link Verification

| From                                          | To                          | Via                                            | Status  | Details                                                                          |
| --------------------------------------------- | --------------------------- | ---------------------------------------------- | ------- | -------------------------------------------------------------------------------- |
| krylov_matvec.jl (TimeDomain apply_lindbladian!) | src/nufft.jl               | `_prefactor_view(oft_nufft_prefactors, w)`      | WIRED   | Line 618, 632, 694, 708 -- called for each energy label in the loop             |
| krylov_matvec.jl (TimeDomain precomputed)     | ws.precomputed_data         | `oft_nufft_prefactors` destructured at lines 591, 668 | WIRED | `_precompute_data` for TimeDomain/TrotterDomain returns `oft_nufft_prefactors` |
| krylov_matvec.jl (BohrDomain apply_lindbladian!) | hamiltonian.bohr_dict      | `keys(hamiltonian.bohr_dict)` bucket iteration  | WIRED   | Lines 478, 539 -- bucket loop iterates over `hamiltonian.bohr_dict` keys        |
| krylov_matvec.jl (_accumulate_dissipator_2op!) | BLAS.gemm!                 | B_dag * A product for anticommutator           | WIRED   | Line 364: `BLAS.gemm!('N','N', CT, B_dag, A, ZT, ws.LdagL)`                   |
| krylov_matvec.jl (BohrDomain apply_lindbladian!) | hamiltonian.bohr_freqs     | `alpha(hamiltonian.bohr_freqs, nu_2)` entrywise | WIRED  | Lines 480, 541 -- `@. ws.jump_oft = alpha(hamiltonian.bohr_freqs, nu_2) * eigenbasis` |

### Requirements Coverage

| Requirement                                             | Status    | Notes                                               |
| ------------------------------------------------------- | --------- | --------------------------------------------------- |
| TimeDomain matvec with NUFFT prefactors (MATVEC-02)    | SATISFIED | Forward + adjoint implemented; zero-alloc verified  |
| TrotterDomain matvec in Trotter eigenbasis (MATVEC-03) | SATISFIED | Forward + adjoint implemented; zero-alloc verified  |
| BohrDomain matvec with 2op dissipator (MATVEC-04)      | SATISFIED | Forward + adjoint + duality check implemented       |

### Anti-Patterns Found

| File                   | Line | Pattern                                                        | Severity | Impact                                                                                    |
| ---------------------- | ---- | -------------------------------------------------------------- | -------- | ----------------------------------------------------------------------------------------- |
| src/krylov_matvec.jl   | 508  | Docstring for `apply_adjoint_lindbladian!` (BohrDomain) shows stale physics-convention formula `A' * rho * B_dag` instead of the corrected kron-convention `B_dag * rho * A` | Warning  | Documentation only -- actual implementation at lines 412-424 uses the correct formula and all 1095 tests pass |

No blocker anti-patterns found.

### Human Verification Required

None required. All correctness conditions are verified programmatically by the test suite (round-trip errors vs dense, duality check, zero-allocation counts).

### Gaps Summary

No gaps. All three observable truths are verified end-to-end:

1. TimeDomain: `_prefactor_view` wired into the hot loop; round-trip KMS + GNS < 1e-12; zero allocations.
2. TrotterDomain: Jump operators built in Trotter eigenbasis (`trotter.eigvecs`); `trotter=TEST_TROTTER` passed to `KrylovWorkspace`; round-trip < 1e-12; zero allocations.
3. BohrDomain: Bucket iteration over `hamiltonian.bohr_dict`; separate forward `_accumulate_dissipator_2op!` and adjoint `_accumulate_adjoint_dissipator_2op!` helpers (not a simple argument swap); round-trip < 1e-12 for KMS and GNS; duality check < 1e-11.

One warning-level documentation inconsistency exists: the outer docstring of `apply_adjoint_lindbladian!` for `BohrDomain` (line 508) still reflects the original plan's physics-convention formula rather than the corrected kron-convention formula. This does not affect correctness or tests.

---

_Verified: 2026-02-20T14:49:37Z_
_Verifier: Claude (gsd-verifier)_
