---
phase: 27-core-matvec-infrastructure
verified: 2026-02-20T11:57:02Z
status: passed
score: 11/11 must-haves verified
---

# Phase 27: Core Matvec Infrastructure Verification Report

**Phase Goal:** Users can apply the Lindbladian superoperator to a density matrix without forming the full dim^2 x dim^2 matrix, validated against dense construction at n=4
**Verified:** 2026-02-20T11:57:02Z
**Status:** passed
**Re-verification:** No -- initial verification

## Goal Achievement

### Observable Truths

Plan 01 truths:

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | KrylovWorkspace struct pre-allocates all scratch matrices at construction time, tied to (config, ham) | VERIFIED | `src/krylov_workspace.jl` lines 30-47: struct with 5 scratch matrices (jump_oft, tmp1, tmp2, LdagL, rho_out) + 2 concrete-typed data vectors; constructor lines 61-103 calls `_precompute_data`, `_precompute_coherent_total_B`, allocates all 5 with `zeros(CT, dim, dim)` |
| 2 | apply_lindbladian! computes L(rho) for EnergyDomain KMS configs with dissipator + coherent term | VERIFIED | `src/krylov_matvec.jl` lines 200-254: dispatches on `AbstractLiouvConfig{EnergyDomain}`; coherent branch lines 215-219 (`-i[B, rho]`); dissipator loop lines 225-251 |
| 3 | apply_lindbladian! computes L(rho) for EnergyDomain GNS configs (coherent term silently skipped) | VERIFIED | GNS config gives `B_total === nothing` (verified by test); `if B !== nothing` guard at line 215 skips coherent term |
| 4 | apply_adjoint_lindbladian! computes L*(rho) with swapped L<->L' and sign-flipped coherent term | VERIFIED | `src/krylov_matvec.jl` lines 276-331: coherent sign flip lines 293-295 (`+1im` / `-1im`); adjoint dissipator via `_accumulate_adjoint_dissipator!` (L' rho L sandwich, same {L'L, rho} anticommutator) |
| 5 | _accumulate_dissipator! helper factors out the Lindblad dissipator formula for reuse | VERIFIED | Lines 32-59: factors D_L(rho) = L rho L' - 0.5{L'L, rho} using BLAS.gemm!/axpy! |

Plan 02 truths:

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 6 | apply_lindbladian! for EnergyDomain KMS matches dense construct_lindbladian() to < 1e-12 for 10 random density matrices at n=4 | VERIFIED | Testset 2 (no coherent) and Testset 3 (with coherent) in test_krylov_matvec.jl; 996 tests all pass (confirmed via test run) |
| 7 | apply_lindbladian! for EnergyDomain GNS matches dense construct_lindbladian() to < 1e-12 for 10 random density matrices at n=4 | VERIFIED | Testset 4 in test_krylov_matvec.jl lines 67-78; 996 tests pass |
| 8 | Coherent term round-trip passes: KMS with_coherent=true matches dense to < 1e-12 | VERIFIED | Testset 3 (lines 50-62) explicitly uses `make_liouv_config(EnergyDomain(); with_coherent=true)`; 996 tests pass |
| 9 | Adjoint round-trip passes: apply_adjoint_lindbladian! matches dense adjoint to < 1e-12 for 10 random density matrices at n=4 | VERIFIED | Testset 5 (lines 83-94) uses `L_dense' * vec(rho)` as reference; 996 tests pass |
| 10 | KrylovWorkspace matvec hot path produces zero heap allocations verified by @allocated | VERIFIED | Testset 7 (lines 120-136): `@test allocs == 0` for both forward and adjoint; BLAS.gemm!/axpy! used throughout instead of mul!/broadcast; concrete-typed `jump_eigenbases::Vector{Matrix{T}}` avoids abstract field boxing |
| 11 | Module includes krylov_workspace.jl and krylov_matvec.jl and exports apply_lindbladian!, apply_adjoint_lindbladian!, KrylovWorkspace | VERIFIED | `src/QuantumFurnace.jl` lines 67-68 (exports), lines 117-118 (includes after furnace.jl) |

**Score:** 11/11 truths verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `src/krylov_workspace.jl` | KrylovWorkspace struct with constructor | VERIFIED | 103 lines; struct has 10 fields including 5 scratch matrices, precomputed_data, B_total, jumps, jump_eigenbases, jump_hermitian; constructor calls _precompute_data and _precompute_coherent_total_B |
| `src/krylov_matvec.jl` | apply_lindbladian!, apply_adjoint_lindbladian!, _accumulate_dissipator! | VERIFIED | 331 lines; contains _krylov_oft!, _accumulate_dissipator!, _accumulate_dissipator_adj_L!, _accumulate_adjoint_dissipator!, _accumulate_adjoint_dissipator_adj_L!, apply_lindbladian!, apply_adjoint_lindbladian! |
| `src/QuantumFurnace.jl` | Module integration -- includes new files and exports new symbols | VERIFIED | `include("krylov_workspace.jl")` at line 117, `include("krylov_matvec.jl")` at line 118; export block at lines 67-68 |
| `test/test_helpers.jl` | GNS config factory for 4-qubit system | VERIFIED | `make_liouv_config_gns` function present at line 250 |
| `test/test_krylov_matvec.jl` | Round-trip correctness tests and allocation regression tests | VERIFIED | 139 lines; 7 testsets covering construction, KMS no-coherent, KMS with-coherent, GNS, adjoint round-trip, duality check, zero-allocation |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `src/krylov_workspace.jl` | `src/furnace_utensils.jl` | `_precompute_data(config, ham_or_trott)` | WIRED | Line 76: `precomputed_data = _precompute_data(config, ham_or_trott)` |
| `src/krylov_workspace.jl` | `src/coherent.jl` | `_precompute_coherent_total_B(jumps, ham_or_trott, config, precomputed_data)` | WIRED | Line 79: `B_total = _precompute_coherent_total_B(jumps, ham_or_trott, config, precomputed_data)` |
| `src/krylov_matvec.jl` | `src/krylov_workspace.jl` | KrylovWorkspace scratch matrices used in hot path | WIRED | `ws.jump_eigenbases`, `ws.jump_hermitian`, `ws.jump_oft`, `ws.rho_out`, `ws.tmp1`, `ws.tmp2`, `ws.LdagL` all accessed throughout matvec functions |
| `test/test_krylov_matvec.jl` | `src/krylov_matvec.jl` | apply_lindbladian! and apply_adjoint_lindbladian! calls | WIRED | Lines 43, 59, 75, 91, 107, 110, 126-135 |
| `test/test_krylov_matvec.jl` | `src/furnace.jl` | construct_lindbladian for dense reference | WIRED | Lines 37, 53, 69, 85 |
| `src/QuantumFurnace.jl` | `src/krylov_workspace.jl` | include statement | WIRED | Line 117: `include("krylov_workspace.jl")` |

Note: `src/krylov_matvec.jl` uses a custom `_krylov_oft!` (defined in the file itself) instead of calling `oft!` from `src/ofts.jl` directly. This is an intentional optimization documented in 27-02-SUMMARY.md: avoids JumpOp abstract field boxing in the hot path. The OFT computation is functionally equivalent.

### Requirements Coverage

No requirements file mapping checked (REQUIREMENTS.md does not map individual requirements to phases at this granularity).

### Anti-Patterns Found

None. Scanned `src/krylov_workspace.jl`, `src/krylov_matvec.jl`, and `test/test_krylov_matvec.jl` for TODO, FIXME, XXX, HACK, PLACEHOLDER, empty returns, and console logging. All files are substantive with no stubs.

### Human Verification Required

None. All goal-critical behaviors are programmatically verifiable and were confirmed by running the full test suite (996 tests pass).

### Gaps Summary

No gaps. All 11 must-haves are verified at all three levels (exists, substantive, wired). The test suite confirms the phase goal: `apply_lindbladian!` and `apply_adjoint_lindbladian!` produce results matching the dense `construct_lindbladian()` to < 1e-12 for 10 random density matrices at n=4 (dim=16), for KMS, GNS, and adjoint variants, with zero heap allocations on the hot path.

Notable implementation decisions beyond the plan (auto-fixed bugs documented in 27-02-SUMMARY.md):
- Separate `_accumulate_adjoint_dissipator!` function to correctly preserve `{L'L, rho}` anticommutator in adjoint (naive L<->L' swap gives wrong anticommutator `{LL', rho}`)
- BLAS.gemm!/axpy! replacing mul!/broadcast to eliminate 156KB per-matvec heap allocations
- Concrete-typed `jump_eigenbases::Vector{Matrix{T}}` to avoid JumpOp abstract field boxing

---

_Verified: 2026-02-20T11:57:02Z_
_Verifier: Claude (gsd-verifier)_
