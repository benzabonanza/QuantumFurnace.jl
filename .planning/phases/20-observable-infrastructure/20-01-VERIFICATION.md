---
phase: 20-observable-infrastructure
verified: 2026-02-17T08:22:17Z
status: passed
score: 6/6 must-haves verified
re_verification: false
---

# Phase 20: Observable Infrastructure Verification Report

**Phase Goal:** Users can construct all observables needed for spectral gap estimation, correctly transformed to the simulation basis
**Verified:** 2026-02-17T08:22:17Z
**Status:** passed
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| #  | Truth                                                                                   | Status     | Evidence                                                                                         |
|----|-----------------------------------------------------------------------------------------|------------|--------------------------------------------------------------------------------------------------|
| 1  | User can call `build_total_magnetization(ham, n)` and get M_z in Hamiltonian eigenbasis | VERIFIED   | `convergence.jl` line 93–111; `V = hamiltonian.eigvecs`, `Mz_eigen = V' * Mz_comp * V`         |
| 2  | User can call `build_total_magnetization(ham, n; trotter=trotter)` for Trotter basis    | VERIFIED   | Same function, `trotter` keyword at line 94; `V = trotter.eigvecs` branch at line 105            |
| 3  | User can call `build_gap_estimation_observables(ham, n)` and get H + M_z bundle         | VERIFIED   | `convergence.jl` line 127–145; returns `vcat([H], mz_obs)` with names `["H", "Mz"]`            |
| 4  | User can call `build_gap_estimation_observables(ham, n; trotter=trotter)` in Trotter basis | VERIFIED | Same function; Trotter branch at line 133 uses `V_T' * hamiltonian.data * V_T`                  |
| 5  | `tr(gibbs * M_z)` matches analytical `sum_i <Z_i>_gibbs / n` for both bases            | VERIFIED   | Testsets 19 (eigenbasis) and 20 (Trotter): `isapprox(..., atol=1e-10)` against analytical value |
| 6  | `tr(gibbs * H)` matches analytical Gibbs energy for both bases                          | VERIFIED   | Testsets 21 (eigenbasis) and 22 (Trotter): `isapprox(..., atol=1e-10)` against `sum(E*boltz)/sum(boltz)` |

**Score:** 6/6 truths verified

### Required Artifacts

| Artifact                       | Provides                                                           | Exists | Substantive | Wired   | Status   |
|--------------------------------|--------------------------------------------------------------------|--------|-------------|---------|----------|
| `src/convergence.jl`           | `build_total_magnetization` (lines 79–111)                        | Yes    | Yes         | Yes     | VERIFIED |
| `src/convergence.jl`           | `build_gap_estimation_observables` (lines 113–145)                | Yes    | Yes         | Yes     | VERIFIED |
| `src/QuantumFurnace.jl`        | Export of both functions (line 47)                                 | Yes    | Yes         | Yes     | VERIFIED |
| `test/test_convergence.jl`     | Regression testsets 19–22 (lines 688–830)                         | Yes    | Yes         | Yes     | VERIFIED |

### Key Link Verification

| From                  | To                             | Via                                          | Status  | Detail                                                                                   |
|-----------------------|--------------------------------|----------------------------------------------|---------|------------------------------------------------------------------------------------------|
| `src/convergence.jl`  | `pad_term([Z], ...)`           | Pauli Z tensor products for M_z construction | WIRED   | Line 100: `Mz_comp .+= Matrix{ComplexF64}(pad_term([Z], num_qubits, i))`                |
| `src/convergence.jl`  | `hamiltonian.eigvecs / trotter.eigvecs` | Basis transform `V' * Mz_comp * V`  | WIRED   | Line 105: `V = trotter !== nothing ? trotter.eigvecs : hamiltonian.eigvecs`; line 108: `V' * Mz_comp * V` |
| `src/convergence.jl`  | `build_total_magnetization`    | Called internally by `build_gap_estimation_observables` | WIRED | Line 130: `mz_obs, mz_names = build_total_magnetization(hamiltonian, num_qubits; trotter=trotter)` |

### Requirements Coverage

| Requirement | Status    | Notes                                                                                                         |
|-------------|-----------|---------------------------------------------------------------------------------------------------------------|
| OBS-01      | SATISFIED | `build_total_magnetization(ham, n)` returns M_z in Hamiltonian eigenbasis; testset 19 verifies Gibbs trace   |
| OBS-02      | SATISFIED | `build_total_magnetization(ham, n; trotter=...)` returns M_z in Trotter basis; testset 20 verifies Gibbs trace |
| OBS-03      | PARTIAL   | `build_gap_estimation_observables` returns H + M_z (no ZZ correlations). ZZ deliberately deferred per user decision. This is an accepted deviation; the deferred scope is documented in PLAN and SUMMARY. |

Note: REQUIREMENTS.md tracking table still shows OBS-01/02/03 as "Pending" — this is a documentation state issue only, not a code gap. The code satisfies OBS-01 and OBS-02 fully, and OBS-03 to the agreed reduced scope.

### Anti-Patterns Found

None. Scanned `src/convergence.jl`, `src/QuantumFurnace.jl`, and `test/test_convergence.jl` for TODO/FIXME/placeholder patterns, empty returns, and stub handlers. No issues found.

### Human Verification Required

None. All truths are programmatically verifiable through the test suite. The implementation uses pure linear algebra (matrix transforms, trace operations) against deterministic test fixtures.

### Deviations from ROADMAP

**ROADMAP success criterion #3** mentions "H, M_z, ZZ correlations" in the gap estimation bundle. Per explicit user decision documented in the PLAN (line 101: "IMPORTANT: Do NOT include ZZ correlations") and SUMMARY (key-decisions), ZZ correlations were deferred to a future phase. The implemented bundle contains H + M_z only. This is an accepted scope reduction, not a defect.

**ROADMAP criterion #4** (regression tests verifying Gibbs trace) is fully satisfied: testsets 19–22 each verify `tr(gibbs * O) ≈ analytical_value` at `atol=1e-10`.

### Commits Verified

Both task commits exist in git history:
- `d3b57bc` — feat(20-01): add build_total_magnetization and build_gap_estimation_observables
- `424ccb8` — test(20-01): add regression tests for magnetization and gap estimation observables

---

_Verified: 2026-02-17T08:22:17Z_
_Verifier: Claude (gsd-verifier)_
