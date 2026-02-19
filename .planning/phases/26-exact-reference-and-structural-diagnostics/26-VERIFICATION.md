---
phase: 26-exact-reference-and-structural-diagnostics
verified: 2026-02-19T03:40:10Z
status: passed
score: 7/7 must-haves verified
re_verification: false
---

# Phase 26: Exact Reference and Structural Diagnostics Verification Report

**Phase Goal:** Researchers can extract exact Lindbladian spectral data and structural diagnostics (anti-Hermitian defect, symmetry sectors, observable overlaps) as ground truth for validating trajectory-based estimates

**Verified:** 2026-02-19T03:40:10Z
**Status:** passed
**Re-verification:** No — initial verification

---

## Goal Achievement

### Observable Truths

| #  | Truth | Status | Evidence |
|----|-------|--------|---------|
| 1  | `extract_leading_eigendata` returns leading eigenvalues sorted by |Re(lambda)| with biorthonormal left+right eigenvectors | VERIFIED | Lines 165-196 of `src/diagnostics.jl`: dense `eigen(L)`, sorted by `abs.(real.(F.values))`, left eigenvectors via `inv(V_full)[1:n_modes, :]'`. 140 tests pass including biorthonormality check `V_left' * V_right ≈ I(10)` at `atol=1e-8`. |
| 2  | `compute_fixed_point_distance` returns trace distance between Lindbladian fixed point and Gibbs state | VERIFIED | Lines 208-218 of `src/diagnostics.jl`: reshapes lambda_1 eigenvector, calls `hermitianize!` and normalizes, then `trace_distance_h`. Test confirms `fp.trace_distance < 0.01` for 3-qubit BohrDomain Lindbladian. |
| 3  | `compute_anti_hermitian_defect` returns defect ratio `||A||/lambda_gap(H_D)` with advisory warning above threshold | VERIFIED | Lines 236-271 of `src/diagnostics.jl`: diagonal KMS similarity transform, Hermitian/anti-Hermitian split, `opnorm(A_part)`, gap from `eigvals(Hermitian(H_part))`. Warning emitted via `@warn` when `defect_ratio > 0.1`. `DefectResult` fields `A_norm`, `H_gap`, `defect_ratio`, `warning`, `threshold` all populated. |
| 4  | Im/Re ratios `|Im(lambda_k)/Re(lambda_k)|` are reported for leading modes | VERIFIED | Lines 189-193 of `src/diagnostics.jl`: `im_re_ratios[k] = abs(imag(eigenvalues[k])) / max(abs(real(eigenvalues[k])), 1e-30)`. Field `im_re_ratios` in `EigenDecompositionResult`. Tested in DIAG-01 testset. |
| 5  | `compute_overlap_coefficients` returns `n_obs x n_modes` matrix using `c_k = Tr[O R_k] * Tr[L_k^dagger (rho_0 - rho_beta)]` | VERIFIED | Lines 290-323 of `src/diagnostics.jl`: explicit `rho_diff = rho0 - Matrix(rho_beta)`, factor `dot(vec(L_k), vec(rho_diff))` (conjugates first arg), factor `tr(O * R_k)`. Test confirms `c_1 < 1e-8` for all observables with any initial state. |
| 6  | `compute_sz_labels` assigns Delta_Sz quantum numbers with purity fractions and near-degeneracy multiplet detection | VERIFIED | Lines 339-391 and 407-441 of `src/diagnostics.jl`: builds Sz operator in eigenbasis, computes `|M_k[i,j]|^2` weight map per (Delta_Sz = sz_i - sz_j), reports dominant sector purity. `detect_multiplets` groups eigenvalues within `rel_tol=0.01`. Steady-state mode test confirms `delta_sz = 0.0` and `purity > 0.95`. |
| 7  | `run_exact_diagnostics` bundles all six DIAG functions into a single call returning `ExactDiagnosticsResult` | VERIFIED | Lines 466-548 of `src/diagnostics.jl`: calls all 6 DIAG functions in sequence, builds 3 default initial states with eigenbasis transforms, returns `ExactDiagnosticsResult`. Bundle test confirms 3 overlaps, 10 sz_labels, non-empty multiplets, correct sub-result types. |

**Score:** 7/7 truths verified

---

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `src/diagnostics.jl` | All DIAG-01 through DIAG-06 functions, result structs, run_exact_diagnostics bundle | VERIFIED | 549 lines. All 7 result structs defined (lines 33-147). DIAG-01 at line 165, DIAG-02 at 208, DIAG-03/04 at 236, DIAG-05 at 290, DIAG-06 at 339, detect_multiplets at 407, bundle at 466. No stubs or TODOs. |
| `test/test_diagnostics.jl` | Tests for all diagnostic functions at n=3 with correctness checks | VERIFIED | 304 lines. 9 testsets covering all DIAG functions, bundle, and custom inputs. 140 tests pass. |
| `src/QuantumFurnace.jl` | Updated includes and exports for diagnostics module | VERIFIED | Line 120: `include("diagnostics.jl")`. Lines 57-61: all 7 result structs and 7 functions exported under "Diagnostics (Phase 26)" comment. |
| `src/convergence.jl` | Updated build_preset_trajectory_observables with 6 canonical observables | VERIFIED | Lines 35-86: 6 observables [Z1, X1, Z1_Zhalf, H, Rand_traceless, Mz_stagg]. MersenneTwister(12345) seed for Rand_traceless. All 200 convergence tests pass. |
| `src/gap_estimation.jl` | Updated eigenbasis_overlap_analysis with correct left+right eigenvector overlap formula | VERIFIED | Lines 290-330+: `rho_beta` keyword, biorthogonal formula `c_k = Tr[O R_k] * Tr[L_k^dagger (rho_0 - rho_beta)]`, backward compatible `rho_beta=nothing` path. All 86 gap estimation tests pass. |
| `test/test_convergence.jl` | Updated tests reflecting 6-observable set | VERIFIED | Tests updated to count 6 observables, check `names == ["Z1", "X1", "Z1_Zhalf", "H", "Rand_traceless", "Mz_stagg"]`, and test individual observable properties at new indices. |
| `test/test_gap_estimation.jl` | Updated tests reflecting 6-observable set | VERIFIED | All old 8-observable references replaced with 6-observable names and counts. `rho_beta` subtraction test added. |

---

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|----|--------|---------|
| `src/diagnostics.jl` | `src/furnace.jl` | `construct_lindbladian` | WIRED | `construct_lindbladian` called in `test/test_diagnostics.jl` line 5 to build the Lindbladian matrix passed to all DIAG functions. Module includes `furnace.jl` before `diagnostics.jl`. |
| `src/diagnostics.jl` | `src/qi_tools.jl` | `trace_distance_h`, `hermitianize!`, `pad_term` | WIRED | `hermitianize!` at line 214, `trace_distance_h` at line 217, `pad_term` at line 348 and 493. All three functions actively used. |
| `test/test_diagnostics.jl` | `src/diagnostics.jl` | Direct function calls | WIRED | Test calls `extract_leading_eigendata`, `compute_fixed_point_distance`, `compute_anti_hermitian_defect`, `compute_overlap_coefficients`, `compute_sz_labels`, `detect_multiplets`, `run_exact_diagnostics` — all functions verified. |
| `src/convergence.jl` | `src/qi_tools.jl` | `pad_term` for Pauli construction | WIRED | `pad_term([Z], ...)` and `pad_term([X], ...)` used at lines 43, 47, 52, 53, 77. |
| `src/gap_estimation.jl` | `src/convergence.jl` | `build_preset_trajectory_observables` called by `estimate_spectral_gap` | WIRED | Confirmed via grep pattern `build_preset_trajectory_observables` present in `gap_estimation.jl`. |

---

### Requirements Coverage

| Requirement | Status | Notes |
|-------------|--------|-------|
| DIAG-01: Leading eigenvalues with left+right eigenvectors, sorted by proximity to zero | SATISFIED | `extract_leading_eigendata` verified. Dense `eigen(L)`, sorted by `abs.(real.(F.values))`, biorthonormal left eigenvectors via `inv(V_full)[1:n_modes,:]'`. Im/Re ratios reported. |
| DIAG-02: Fixed point distance to Gibbs state | SATISFIED | `compute_fixed_point_distance` verified. `trace_distance < 0.01` for 3-qubit test system (0.002 expected from Gaussian filter; plan acknowledges tolerance). |
| DIAG-03: Anti-Hermitian defect via KMS similarity transform | SATISFIED | `compute_anti_hermitian_defect` verified. KMS diagonal transform correct, Hermitian/anti-Hermitian split, `opnorm(A_part)`. |
| DIAG-04: Defect ratio determines appropriateness of real-exponential fitting | SATISFIED | `DefectResult.defect_ratio = A_norm / H_gap`, advisory `@warn` when `> 0.1`. Does not gate any computation. |
| DIAG-05: Observable overlap coefficients `c_k = Tr[O R_k] * Tr[L_k^dagger (rho_0 - rho_beta)]` | SATISFIED | `compute_overlap_coefficients` and `eigenbasis_overlap_analysis(rho_beta=...)` both implement the correct formula. Steady-state coefficient `c_1 < 1e-8` confirmed. |
| DIAG-06: Delta_Sz symmetry sector labels with near-degeneracy detection | SATISFIED | `compute_sz_labels` verified. `detect_multiplets` groups near-degenerate eigenvalues. Steady-state mode confirmed `delta_sz = 0.0`, `purity > 0.95`. |

---

### Anti-Patterns Found

None. The implementation shows:
- No TODO/FIXME/PLACEHOLDER comments in `src/diagnostics.jl`
- No empty function bodies or stub returns
- No `console.log`-only handlers
- No static/hardcoded return values where computation is expected
- `@warn` emission is real behavior (advisory defect warning), not a stub

---

### Human Verification Required

The following items are verified programmatically via test suite (140 diagnostics + 200 convergence + 86 gap estimation = 426 total tests). No human verification is required for goal achievement.

However, the following aspects would benefit from researcher review at n=4,6:

1. **Arpack vs dense eigen performance at n=6 (4096x4096 Lindbladian)**
   - Test: Run `extract_leading_eigendata` on a 4096x4096 matrix and measure wall time
   - Expected: Completes in reasonable time (seconds to low minutes)
   - Why human: Test suite uses 3-qubit (n=3, 64x64) system; n=6 timing not exercised in CI

2. **Physical correctness of n=6 zero-overlap mystery resolution**
   - Test: Run `compute_sz_labels` on n=6 system gap mode and check Delta_Sz label
   - Expected: Gap mode shows pure Delta_Sz sector different from zero, explaining zero overlap with Sz-conserving observables
   - Why human: Requires n=6 reference computation not part of automated test suite

---

### Gaps Summary

No gaps. All seven must-have truths verified across three levels (exists, substantive, wired). All 426 tests across diagnostics, convergence, and gap estimation pass with zero failures.

The one deviation from the original plan (fixed point trace distance tolerance relaxed from `< 1e-10` to `< 0.01`) reflects physical reality: the 3-qubit BohrDomain Lindbladian with Gaussian filter smoothing has approximately 0.002 trace distance from the exact Gibbs state. This is expected behavior, documented in the SUMMARY, and verified in the tests with an explanatory comment.

---

_Verified: 2026-02-19T03:40:10Z_
_Verifier: Claude (gsd-verifier)_
