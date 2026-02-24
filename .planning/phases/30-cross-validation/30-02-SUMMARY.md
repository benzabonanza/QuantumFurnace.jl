---
phase: 30-cross-validation
plan: 02
subsystem: testing
tags: [krylov, cross-validation, spectral-gap, n6, env-gated, production-scale]

# Dependency graph
requires:
  - phase: 30-cross-validation
    plan: 01
    provides: compare_krylov_dense helper, print_gap_summary, on_failure_diagnostics, test file structure
  - phase: 29-eigensolver-integration
    provides: krylov_spectral_gap(), KrylovGapResult, extract_leading_eigendata()
provides:
  - n=6 env-gated KMS cross-validation tests for all 4 domains at atol=1e-6
  - make_n6_liouv_config, make_n6_test_system, make_n6_thermalize_config factory functions
  - Production-scale (dim=64, dim^2=4096) validation of Krylov spectral gap
affects: [31-benchmarks]

# Tech tracking
tech-stack:
  added: []
  patterns: [env-gated n=6 test block, n=6 config factories with num_qubits=6]

key-files:
  created: []
  modified:
    - test/test_krylov_crossvalidation.jl

key-decisions:
  - "atol=1e-6 for n=6 cross-validation (looser than n=4's 1e-8 due to larger Krylov subspace approximation error at dim^2=4096)"
  - "No GNS at n=6 per locked decision (n=4 sufficient for balance-type correctness)"
  - "Dedicated n6_trotter and n6_trotter_sys for TrotterDomain (separate eigenbasis from Hamiltonian)"

patterns-established:
  - "n=6 factory pattern: make_n6_liouv_config/make_n6_test_system/make_n6_thermalize_config with num_qubits=6"
  - "Env-gated production tests: QUANTUMFURNACE_FULL_TESTS=true skipping pattern"

# Metrics
duration: 2min
completed: 2026-02-24
---

# Phase 30 Plan 02: n=6 Env-Gated Cross-Validation Summary

**n=6 KMS cross-validation across 4 domains at atol=1e-6, gated behind QUANTUMFURNACE_FULL_TESTS=true, with dedicated 6-qubit test system factories**

## Performance

- **Duration:** 2 min
- **Started:** 2026-02-24T14:31:28Z
- **Completed:** 2026-02-24T14:33:14Z
- **Tasks:** 1
- **Files modified:** 1

## Accomplishments
- XVAL-02: n=6 KMS cross-validation tests for all 4 domains (EnergyDomain, TimeDomain, TrotterDomain, BohrDomain) at atol=1e-6
- Three n=6 factory functions (make_n6_liouv_config, make_n6_test_system, make_n6_thermalize_config) with num_qubits=6 and same physical parameters as n=4
- Env-gated behind QUANTUMFURNACE_FULL_TESTS=true, with @info skip message when not set
- TrotterDomain uses dedicated n6_trotter (TrottTrott) and n6_trotter_sys.jumps (separate eigenbasis from Hamiltonian eigenvectors)

## Task Commits

Each task was committed atomically:

1. **Task 1: Add n=6 test system construction and env-gated cross-validation block** - `74075da` (feat)

## Files Created/Modified
- `test/test_krylov_crossvalidation.jl` - Added make_n6_liouv_config, make_n6_test_system, make_n6_thermalize_config factories and env-gated n=6 KMS testset with 4 domain sub-testsets

## Decisions Made
None - followed plan as specified. All implementation decisions (atol=1e-6, KMS-only at n=6, env gating pattern, dedicated Trotter eigenbasis) were locked in CONTEXT.md and PLAN.md.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Phase 30 cross-validation complete: n=4 KMS/GNS + L-vs-E convergence (Plan 01) and n=6 KMS env-gated (Plan 02) all implemented
- Phase 31 (benchmarks) can proceed -- all cross-validation infrastructure is in place
- To run n=6 tests: `QUANTUMFURNACE_FULL_TESTS=true julia --project -e 'using Pkg; Pkg.test()'`

## Self-Check: PASSED

All files exist, all commits verified, all content checks pass.

---
*Phase: 30-cross-validation*
*Completed: 2026-02-24*
