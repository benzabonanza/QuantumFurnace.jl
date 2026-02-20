# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-19)

**Core value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers
**Current focus:** v1.4 Spectral Gap Refinement -- Phase 27 (Two-Exponential Fitting Infrastructure)

## Current Position

Phase: 27 of 30 (Two-Exponential Fitting Infrastructure)
Plan: 0 of TBD in current phase
Status: Ready to plan
Last activity: 2026-02-20 -- Completed quick task 34: Extend exact diagnostics to support TrotterDomain

Progress: [░░░░░░░░░░] 0%

## Performance Metrics

**Velocity:**
- Total plans completed: 71 (v1.0: 10, v1.1: 16, quick: 18, v1.2: 12, cleanup: 3, v1.3: 10, v1.4: 2)

**By Milestone:**

| Milestone | Phases | Plans | Timeline |
|-----------|--------|-------|----------|
| v1.0 Trajectories | 1-5 | 10 | 2026-02-13 to 2026-02-14 |
| v1.1 Reduce | 6-11 | 16 (+5 quick) | 2026-02-15 |
| v1.2 Multi-threading | 12-19 | 15 (+3 quick) | 2026-02-15 to 2026-02-16 |
| v1.3 Mixing Time | 20-25 | 10 (+11 quick) | 2026-02-17 to 2026-02-18 |
| v1.4 Spectral Gap Refinement | 26-30 | 2 | 2026-02-19 to - |

## Accumulated Context

### Decisions

All decisions logged in PROJECT.md Key Decisions table.
Recent context for v1.4:
- All trajectory runs use `with_coherent=true`, 4 threads, Trotter steps = 10, delta is the Trotter parameter
- Batch-level bootstrap (~3 MB) chosen over per-trajectory storage (~640 MB+) -- memory safe within 40 GB budget
- Two-exponential fitting is the root cause fix for non-monotonic delta-scaling (Quick-32 conclusion)
- Lindbladian fixed point (not Gibbs state) must be used as steady-state for TrotterDomain lambda_eff
- Dense eigen() used for DIAG-01 (Arpack cannot compute left eigenvectors needed for overlap formula)
- DIAG-03/04 defect_ratio > 0.1 is advisory-only warning, does not gate computation
- Canonical 6-observable set: Z1, X1, Z1_Zhalf, H, Rand_traceless, Mz_stagg (replaces v1.3 8-obs set)
- Biorthogonal overlap formula: c_k = Tr[O R_k] * Tr[L_k^dagger (rho_0 - rho_beta)] via explicit left+right eigenvectors
- run_exact_diagnostics accepts basis_eigvecs keyword for TrotterDomain; controls default observables, initial states, and Sz labels

### Pending Todos

None

### Blockers/Concerns

- Phase 27: `_thermalize_to_liouv_config` field mapping needs targeted audit before implementation
- Phase 29: Per-batch seeding arithmetic needs verification against existing `Xoshiro(seed + traj_id)` scheme
- Phase 26: n=4 has g2/g1 ~ 2x, right at identifiability boundary -- validate Prony initialization early with real data

### Quick Tasks Completed

| # | Description | Date | Commit | Directory |
|---|-------------|------|--------|-----------|
| 19 | Fix failing test after EnergyDomain to TrotterDomain rename in GNS trajectory | 2026-02-16 | 0683acc | [19-fix-failing-test-after-energydomain-to-t](./quick/19-fix-failing-test-after-energydomain-to-t/) |
| 20 | Fix basis mismatch in GNS TrotterDomain tests (0.83 -> 0.08 gap) | 2026-02-16 | 0161e01 | [20-debug-gns-trotterdomain-0-83-gap-suspect](./quick/20-debug-gns-trotterdomain-0-83-gap-suspect/) |
| 21 | Fix test errors after removing transform_jumps_to_basis | 2026-02-16 | b2e4123 | [21-fix-test-errors-after-removing-transform](./quick/21-fix-test-errors-after-removing-transform/) |
| 22 | Fix trajectory delta_eff double-counting (bare delta for CPTP channel) | 2026-02-17 | dc83bf0 | [22-fix-trajectory-delta-eff-double-counting](./quick/22-fix-trajectory-delta-eff-double-counting/) |
| 23 | Improve spectral gap estimation: smallest-gap selection criterion | 2026-02-17 | c5f6c68 | [23-improve-spectral-gap-estimation-by-findi](./quick/23-improve-spectral-gap-estimation-by-findi/) |
| 24 | Add XX_avg, YY_avg, ZZ_avg two-site correlations to gap estimation | 2026-02-17 | 595aba3 | [24-add-two-site-correlations-to-spectral-ga](./quick/24-add-two-site-correlations-to-spectral-ga/) |
| 25 | Diagnose n=6 zero gap-mode overlap: confirmed k=pi momentum sector | 2026-02-18 | cf61f88 | [25-diagnose-n-6-gap-mode-momentum-sector-co](./quick/25-diagnose-n-6-gap-mode-momentum-sector-co/) |
| 26 | Add Mz_stagg/Z1 observables; n=6 gap has SU(2) symmetry protection | 2026-02-18 | 9e5b0de | [26-add-staggered-non-symmetric-observables-](./quick/26-add-staggered-non-symmetric-observables-/) |
| 27 | Add XZ_stagg observable; n=6 gap protected beyond SU(2) | 2026-02-18 | a0f4c62 | [27-add-symmetry-breaking-observable-for-n-6](./quick/27-add-symmetry-breaking-observable-for-n-6/) |
| 28 | Disorder breaks n=6 gap-mode symmetry; estimation biased at 34% | 2026-02-18 | 82a8014 | [28-test-gap-mode-coupling-with-random-field](./quick/28-test-gap-mode-coupling-with-random-field/) |
| 29 | XX_stagg has zero gap-mode overlap; XX correlations don't couple to gap | 2026-02-18 | 4914103 | [29-rigorous-disordered-heisenberg-gap-estim](./quick/29-rigorous-disordered-heisenberg-gap-estim/) |
| 30 | Gap error NOT O(delta); Richardson extrapolation ineffective (1.0x) | 2026-02-18 | 502afe4 | [30-validate-trotter-delta-order-gap-error-a](./quick/30-validate-trotter-delta-order-gap-error-a/) |
| 31 | Longer mixing + uniform psi0 does NOT fix O(delta) scaling (199x ratio spread) | 2026-02-18 | 06fcff8 | [31-delta-scaling-revalidation-with-longer-m](./quick/31-delta-scaling-revalidation-with-longer-m/) |
| 32 | Trajectory simulation confirmed correct; non-monotonic gap estimation is fitting artifact | 2026-02-18 | dfd3baf | [32-investigate-and-fix-non-monotonic-delta-](./quick/32-investigate-and-fix-non-monotonic-delta-/) |
| 34 | Extend exact diagnostics to support TrotterDomain (basis_eigvecs keyword) | 2026-02-20 | 5b48b1e | [34-extend-exact-diagnostics-to-support-trot](./quick/34-extend-exact-diagnostics-to-support-trot/) |

## Session Continuity

Last session: 2026-02-20
Stopped at: Quick task 34 complete -- TrotterDomain diagnostics support added
Resume file: None
