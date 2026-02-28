---
phase: 38-test-cleanup
verified: 2026-02-28T10:30:00Z
status: human_needed
score: 3/4 must-haves verified
human_verification:
  - test: "Run full test suite with julia --project -e 'using Pkg; Pkg.test()'"
    expected: "All 1057 tests pass with zero failures"
    why_human: "Cannot run Julia test suite in this verification environment; need actual Julia runtime with project dependencies"
---

# Phase 38: Test Cleanup Verification Report

**Phase Goal:** Test infrastructure is consolidated with parametrized helpers, informative output, and validated thresholds so that tests are maintainable and trustworthy
**Verified:** 2026-02-28T10:30:00Z
**Status:** human_needed
**Re-verification:** No -- initial verification

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | Config factory functions in test_helpers.jl are parametrized by system size and construction type, eliminating duplicate setup patterns across test files | VERIFIED | Single `make_config(sim, domain; num_qubits, construction, delta, mixing_time)` at line 185-208 of test/test_helpers.jl; single `make_test_system(; num_qubits, trotter)` at line 117-148; zero references to old 8 factory names (grep returns no matches); 124 uses of `make_config(` across 22 test files; 145 uses of N3_* globals across 14 files |
| 2 | Every @testset block prints @info showing what is being tested and key numerical results (trace distances, gaps, allocation counts) | VERIFIED | 204 @info statements across 19 test files; automated script confirmed every numerical @test in all 15 test files is followed by @info with label, value, and threshold; loop-summary pattern used for multi-iteration tests (CPTP, matvec round-trips); test_results.jl and test_observable_trajectories.jl correctly skip @info with policy comment (all exact/structural checks) |
| 3 | Previously dubious test thresholds are reviewed and either tightened with documented rationale or relaxed with explanation of why the original threshold was wrong | VERIFIED | 163 inline threshold rationale comments across 18 test files; documents FP accumulation theory (O(DIM^2 * eps)), KrylovKit tolerance relationships, statistical noise (1/sqrt(N)), empirical bounds with margin factors; all existing thresholds confirmed appropriate after error analysis |
| 4 | All tests pass and the full test suite runs without regressions from the restructure | NEEDS HUMAN | Cannot run Julia test suite in verification environment; summaries claim 1057 tests pass across all 5 sub-plans; git history shows clean progression with no fix-up commits; test file structure is intact and self-consistent |

**Score:** 3/4 truths verified (1 needs human confirmation)

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `test/test_helpers.jl` | Unified make_config, make_test_system, N3_* globals, ALL_DOMAINS | VERIFIED | Single `make_config` (1 definition), single `make_test_system` (1 definition), N3_HAM/N3_JUMPS/N3_GIBBS/N3_DIM/N3_TROTTER/N3_TROTTER_JUMPS constants, ALL_DOMAINS with 4 domain singletons |
| `test/runtests.jl` | Updated test runner with trajectory_validation integration | VERIFIED | QUANTUMFURNACE_FULL_TESTS gate at lines 29-34; includes trajectory_validation/run_trajectory_validation.jl and run_convergence_tests.jl; @info skip message when gated |
| `test/old_tests/` | Deleted | VERIFIED | Directory does not exist (`test -d` returns false) |
| `test/test_cptp.jl` | Loop-summary @info pattern for CPTP completeness checks | VERIFIED | 3 @info statements with loop-summary pattern (max_err tracking), threshold rationale in header comment |
| `test/test_dm_detailed_balance.jl` | @info for trace distance and domain hierarchy tests | VERIFIED | 5 @info statements, @info placed after @test per locked decision |
| `test/test_allocation.jl` | @info after every allocation @test, threshold rationale comments | VERIFIED | 8 @info statements showing allocs_bytes and threshold, budget rationale documented |
| `test/test_dm_scaling.jl` | @info for Euler/OFT/NUFFT scaling tests with ratio diagnostics | VERIFIED | 18 @info statements; ratio tests show i, ratio, expected, lower, upper bounds |
| `test/test_gns_trajectory.jl` | @info for GNS trajectory trace distance and convergence tests | VERIFIED | 9 @info statements covering gap bounds, CPTP summary, trajectory convergence |
| `test/test_diagnostics.jl` | @info after numerical DIAG-01 through DIAG-06 tests | VERIFIED | 38 @info lines covering DIAG-01 through DIAG-06 plus multiplet and bundle tests |
| `test/test_convergence.jl` | @info after numerical convergence tests | VERIFIED | 36 @info lines covering Gibbs helpers, windowed relative change, CONV-01 through CONV-05 |
| `test/test_krylov_matvec.jl` | @info for matvec round-trips and allocation checks | VERIFIED | 26 @info: 16 loop-summary round-trips, 6 allocation checks, 3 duality checks |
| `test/test_krylov_eigsolve.jl` | @info for eigsolve accuracy and guard rail tests | VERIFIED | 18 @info covering Chen channel properties, gap accuracy, domain coverage |
| `test/test_krylov_crossvalidation.jl` | @info for cross-validation eigenvalue comparisons, Printf tables preserved | VERIFIED | 17 @info for XVAL-01 through XVAL-04; 7 Printf/printf references preserved |

### Key Link Verification

| From | To | Via | Status | Details |
|------|-----|-----|--------|---------|
| test/test_helpers.jl | all 22 test files | make_config replaces 8 old factory functions | WIRED | 124 uses of `make_config(` across 22 files; zero references to old factory names |
| test/test_helpers.jl | all test files | N3_ replaces SMALL_ globals | WIRED | 145 uses of N3_HAM/N3_JUMPS/N3_GIBBS/N3_DIM across 14 files; zero references to SMALL_* |
| test/test_helpers.jl | test files | TOL_EXACT, TOL_QUADRATURE, TOL_DELTA tolerance constants | WIRED | Constants defined in test_helpers.jl, referenced across test files |
| test/runtests.jl | trajectory_validation/ | QUANTUMFURNACE_FULL_TESTS gate | WIRED | Gate at lines 29-34, includes both run_trajectory_validation.jl and run_convergence_tests.jl |
| test/test_krylov_crossvalidation.jl | test/test_helpers.jl | make_config with num_qubits=6 for large system tests | WIRED | Local n=6 factories deleted, replaced with make_config(...; num_qubits=6) calls |

### Requirements Coverage

No specific REQUIREMENTS.md entries mapped to Phase 38.

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| test/test_threading.jl | 31 | `@test true  # placeholder so testset is not empty` | Info | Standard Julia pattern for gated testsets when nthreads=1; not a real placeholder |

No blockers or warnings found.

### Human Verification Required

### 1. Full Test Suite Execution

**Test:** Run `julia --project -e 'using Pkg; Pkg.test()'`
**Expected:** All 1057 tests pass with zero failures, @info output visible in test log showing computed values and thresholds
**Why human:** Cannot run Julia test suite in this verification environment; need Julia runtime with project dependencies compiled

### Gaps Summary

No gaps found. All automated checks pass. The only remaining item is human confirmation that the full test suite runs without regressions, which requires a Julia runtime environment with the project's dependencies available. The codebase artifacts, wiring, and instrumentation are all verified correct.

The consolidation is thorough:
- 8 factory functions reduced to 1 parametrized `make_config`
- 2 system factories reduced to 1 parametrized `make_test_system`
- SMALL_* renamed to N3_* consistently across all files
- 204 @info statements provide comprehensive test observability
- 163 threshold rationale comments document the scientific basis for every numerical bound
- 7 dead test files deleted, trajectory_validation properly integrated
- Printf diagnostic tables preserved where appropriate

---

_Verified: 2026-02-28T10:30:00Z_
_Verifier: Claude (gsd-verifier)_
