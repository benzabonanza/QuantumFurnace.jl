---
phase: 29-eigensolver-integration
verified: 2026-02-24T12:52:00Z
status: human_needed
score: 10/10 must-haves verified (structural)
human_verification:
  - test: "Run julia --project -e 'using Pkg; Pkg.test()' and confirm 'Krylov Eigsolve' testset passes all 8 sub-testsets with zero failures"
    expected: "All 8 testsets pass: apply_delta_channel! round-trip atol=1e-12, Lindbladian gap rtol=1e-6, channel gap rtol=1e-3, all 4 domains positive gap, guard rail ErrorException, eigenvalue conversion round-trip"
    why_human: "Julia is not installed in the sandbox environment; tests could not be executed. The SUMMARY explicitly documents this limitation. All code is structurally correct and wired, but runtime correctness requires Julia execution."
  - test: "Run krylov_spectral_gap on a known n=4 EnergyDomain KMS system and compare result.spectral_gap to the value from extract_leading_eigendata"
    expected: "result.spectral_gap matches dense reference to rtol=1e-6; trace_distance_h(Hermitian(result.fixed_point), TEST_GIBBS) < 1e-4"
    why_human: "Numerical correctness of the KrylovKit Arnoldi convergence and eigenvalue sorting cannot be verified without executing Julia"
  - test: "Run krylov_spectral_gap with a ThermalizeConfig (channel path) and confirm result.channel_eigenvalues is non-nothing and delta_used == 0.01"
    expected: "Channel eigenvalue near 1 (|mu[1]| ~ 1.0 atol=0.01); spectral_gap matches Lindbladian path to rtol=1e-3"
    why_human: "Channel path mu -> lambda_L conversion (lambda_L = (mu-1)/delta) requires runtime verification"
---

# Phase 29: Eigensolver Integration Verification Report

**Phase Goal:** Users can compute spectral gaps via a single `krylov_spectral_gap()` call that wraps KrylovKit eigsolve, with both Lindbladian (`:LR`) and CPTP channel (`:LM`) targeting
**Verified:** 2026-02-24T12:52:00Z
**Status:** human_needed — all structural checks pass; functional verification requires Julia runtime
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|---------|
| 1 | KrylovKit in Project.toml [deps] and [compat], `using KrylovKit` in module | VERIFIED | `Project.toml` line 13: KrylovKit uuid present in [deps]; line 42: `KrylovKit = "0.8, 0.9, 0.10"` in [compat]; `src/QuantumFurnace.jl` line 23: `using KrylovKit` |
| 2 | `KrylovGapResult` struct has all 10 required fields | VERIFIED | `src/krylov_eigsolve.jl` lines 40-51: all 10 fields match spec: eigenvalues, spectral_gap, fixed_point, gap_mode, converged, matvec_count, num_restarts, normres, channel_eigenvalues, delta_used |
| 3 | `krylov_spectral_gap(config::AbstractLiouvConfig, ...)` wraps KrylovKit with `:LR` targeting | VERIFIED | `src/krylov_eigsolve.jl` lines 259-328: function exists, calls `_eigsolve_with_retry(..., :LR; ...)` at line 291 |
| 4 | `krylov_spectral_gap(config::AbstractThermalizeConfig, ...)` wraps KrylovKit with `:LM` targeting and converts `(mu-1)/delta` | VERIFIED | `src/krylov_eigsolve.jl` lines 356-439: function exists, calls `_eigsolve_with_retry(..., :LM; ...)` at line 395, conversion on line 404 |
| 5 | `apply_delta_channel!(ws, rho, delta, config_liouv, ham)` computes E(rho) = rho + delta*L(rho) | VERIFIED | `src/krylov_eigsolve.jl` lines 214-226: calls `apply_lindbladian!` then `@. ws.rho_out = rho + delta * ws.rho_out` |
| 6 | `_eigsolve_with_retry` retries with 1.5x krylovdim up to 3 times, errors on total failure | VERIFIED | `src/krylov_eigsolve.jl` lines 107-133: loop `1:(max_retries+1)`, `ceil(Int, current_krylovdim * 1.5)`, error on exhaustion |
| 7 | Pre-flight memory guard issues `@warn` when krylovdim * 4^n * 16 * 1.5 > 80% of `Sys.free_memory()` | VERIFIED | `src/krylov_eigsolve.jl` lines 68-78: formula matches spec, advisory `@warn`, returns nothing |
| 8 | `krylovdim > howmany` guard raises error | VERIFIED | `src/krylov_eigsolve.jl` line 271: `krylovdim > howmany \|\| error(...)` (Lindbladian path); line 368 (channel path) |
| 9 | Module exports `KrylovGapResult`, `krylov_spectral_gap`, `apply_delta_channel!` | VERIFIED | `src/QuantumFurnace.jl` line 72: `export KrylovGapResult, krylov_spectral_gap, apply_delta_channel!`; include at line 123 after dependency chain |
| 10 | `test/test_krylov_eigsolve.jl` has 8 testsets; included in `test/runtests.jl` | VERIFIED | `test/test_krylov_eigsolve.jl` is 194 lines with 8 named `@testset` blocks and 42 total test assertions; `test/runtests.jl` line 27: `include("test_krylov_eigsolve.jl")` |

**Score:** 10/10 truths structurally verified

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `src/krylov_eigsolve.jl` | KrylovGapResult, krylov_spectral_gap (2 dispatch), apply_delta_channel!, retry, memory guard | VERIFIED | 439 lines, all required functions present, no stubs/placeholders |
| `src/QuantumFurnace.jl` | using KrylovKit, include krylov_eigsolve.jl, exports | VERIFIED | All three changes present at lines 23, 72, 123 |
| `Project.toml` | KrylovKit in [deps] and [compat] | VERIFIED | [deps] line 13, [compat] line 42 |
| `test/test_krylov_eigsolve.jl` | 8 testsets covering all phase requirements | VERIFIED | 194 lines, 8 @testset blocks, 42 assertions, no TODOs/placeholders |
| `test/runtests.jl` | include test_krylov_eigsolve.jl | VERIFIED | Line 27 |

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `src/krylov_eigsolve.jl` | `src/krylov_matvec.jl` | `apply_lindbladian!` called inside closures | WIRED | Line 283 (Lindbladian closure), line 222 (apply_delta_channel!) |
| `src/krylov_eigsolve.jl` | `src/krylov_workspace.jl` | `KrylovWorkspace(config, ham, jumps)` allocation | WIRED | Line 275 (Lindbladian path), line 379 (channel path) |
| `src/krylov_eigsolve.jl` | `src/structs.jl` | Config type dispatch: `AbstractLiouvConfig` -> `:LR`, `AbstractThermalizeConfig` -> `:LM` | WIRED | Function signatures at lines 259 and 356 use the abstract type hierarchy from structs.jl |
| `test/test_krylov_eigsolve.jl` | `src/krylov_eigsolve.jl` | `krylov_spectral_gap` and `apply_delta_channel!` calls | WIRED | 12 calls to krylov_spectral_gap, 1 to apply_delta_channel! in the test file |
| `test/test_krylov_eigsolve.jl` | `src/diagnostics.jl` | `extract_leading_eigendata` for dense reference | WIRED | Called in testsets 3, 4, 5 for cross-validation |
| `test/test_krylov_eigsolve.jl` | `test/test_helpers.jl` | `make_liouv_config`, `make_thermalize_config`, `TEST_HAM`, `TEST_JUMPS`, `TEST_GIBBS` | WIRED | All fixtures used correctly; test_helpers.jl included by runtests.jl before test_krylov_eigsolve.jl |
| Channel path closure | `apply_delta_channel!` | Called within `channel_matvec` closure | WIRED | Line 387: `apply_delta_channel!(ws, rho, delta, config_liouv, hamiltonian)` |
| `copy(vec(ws.rho_out))` | Aliasing prevention | Critical copy in both closures | WIRED | Line 284 (Lindbladian), line 388 (channel): `return copy(vec(ws.rho_out))` |

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| None | — | — | — | No stubs, no TODOs, no empty handlers, no placeholder returns found in either implementation or test file |

### Human Verification Required

Both SUMMARYs explicitly document that Julia is not installed in the sandbox and no test suite run was possible. All code was reviewed structurally. The following runtime verifications are required:

#### 1. Full Test Suite Pass

**Test:** Run `julia --project -e 'using Pkg; Pkg.test()'` in the project directory
**Expected:** All 8 Krylov Eigsolve testsets pass with zero failures, including:
- apply_delta_channel! round-trip matches dense `(I + delta*L)*vec(rho)` to atol=1e-12 for 5 random density matrices
- Lindbladian spectral gap matches extract_leading_eigendata to rtol=1e-6 at n=4 (EnergyDomain KMS)
- GNS path spectral gap matches dense to rtol=1e-6; trace-normalized fixed point
- Channel path gap matches to rtol=1e-3; channel_eigenvalues populated; delta_used == 0.01
- All 4 domains (Energy, Time, Trotter, Bohr) return result.spectral_gap > 0
- `krylov_spectral_gap(...; krylovdim=2, howmany=4)` raises ErrorException
- Eigenvalue conversion round-trip: `1.0 .+ delta .* result.eigenvalues ≈ result.channel_eigenvalues` to atol=1e-10
**Why human:** Julia runtime not available in sandbox; no execution was possible

#### 2. Krylov vs Dense Spectral Gap Numerical Accuracy

**Test:** In a Julia REPL: call `krylov_spectral_gap` with `LiouvConfig(EnergyDomain())` and compare to `extract_leading_eigendata(construct_lindbladian(...))` result
**Expected:** `abs(result.spectral_gap - dense_result.spectral_gap) / dense_result.spectral_gap < 1e-6`
**Why human:** KrylovKit convergence depends on runtime numerical behavior; the specific n=4 problem may have near-degenerate eigenvalues affecting convergence

#### 3. Channel Path Eigenvalue Conversion

**Test:** Call `krylov_spectral_gap` with `ThermalizeConfig(EnergyDomain(); delta=0.01)` and verify `result.channel_eigenvalues` and conversion
**Expected:** `abs(result.channel_eigenvalues[1]) ≈ 1.0` (atol=0.01); `result.delta_used == 0.01`; `isapprox(1.0 .+ 0.01 .* result.eigenvalues, result.channel_eigenvalues; atol=1e-10)`
**Why human:** Requires Julia runtime to confirm the conversion formula produces correct values

### Gaps Summary

No gaps found. All phase artifacts are substantive and properly wired. The phase goal is structurally achieved: `krylov_spectral_gap()` exists as a single-call API with both dispatch paths (`:LR` for Lindbladian, `:LM` for CPTP channel), backed by comprehensive tests covering all requirements. The only outstanding item is runtime validation which requires Julia execution not available in this environment.

---

_Verified: 2026-02-24T12:52:00Z_
_Verifier: Claude (gsd-verifier)_
