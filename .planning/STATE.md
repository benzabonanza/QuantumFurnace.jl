# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-16)

**Core value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers
**Current focus:** Phase 25 (Spectral Gap Validation Overhaul) -- COMPLETE

## Current Position

Phase: 25 (Spectral Gap Validation Overhaul) -- COMPLETE
Plan: 3 of 3 in current phase (Plan 03 COMPLETE)
Status: Phase 25 complete. Quick-26 added Mz_stagg/Z1 observables (7 total); n=6 gap mode has additional SU(2) symmetry beyond translational.
Last activity: 2026-02-18 -- Quick task 26 executed (staggered observables)

Progress: [##############################] 48/48 plans (v1.0-v1.3) + 3/3 Phase 25

## Performance Metrics

**Velocity:**
- Total plans completed: 60 (v1.0: 10, v1.1: 16, quick: 12, v1.2: 12, cleanup: 3, v1.3: 7, Phase 25: 3)

**By Milestone:**

| Milestone | Phases | Plans | Timeline |
|-----------|--------|-------|----------|
| v1.0 Trajectories | 1-5 | 10 | 2026-02-13 to 2026-02-14 |
| v1.1 Reduce | 6-11 | 16 (+5 quick) | 2026-02-15 |
| v1.2 Multi-threading | 12-19 | 15 (+3 quick) | 2026-02-15 to 2026-02-16 |
| v1.3 Mixing Time | 20-24 | 7 | 2026-02-16 to 2026-02-17 |

## Accumulated Context

### Decisions

All decisions logged in PROJECT.md Key Decisions table.
Key context for v1.3:
- Use `abs(real(spectral_gap))` not `abs(spectral_gap)` for cross-validation (complex eigenvalues)
- Follow `build_convergence_observables` pattern for basis transforms (avoids pitfall from quick-task-20)
- LsqFit.jl is the single new dependency (Levenberg-Marquardt + CIs + bounds)
- Phases 20 and 21 can execute in parallel (independent work)
- Phase 20: Unified trotter keyword pattern (not separate _trotter suffix) for new observable builders
- Phase 20: M_z = sum(Z_i)/n per-site normalization, H + M_z bundle (no ZZ correlations)
- Phase 21: FitResult struct wraps LsqFit output; _IDX_A=1, _IDX_GAP=2, _IDX_C=3 parameter ordering
- Phase 21: R-squared not clamped (negative = valid diagnostic); no weight vector to curve_fit
- Phase 22: ObservableTrajectoryResult in trajectories.jl (not structs.jl); inner/outer constructor pattern for Aqua compliance
- Phase 22: reconstruct_dm=true reuses _run_chunk_with_obs!; reconstruct_dm=false uses new _run_chunk_obs_only!
- Phase 23: SpectralGapResult uses concrete types only (no type parameter) for Aqua compliance
- Phase 23: Best observable selected by smallest gap among converged + gap > 0 + R-squared > 0.8 fits (quick-23 fix)
- Phase 23: gap_estimation.jl included after fitting.jl and before results.jl for correct dependency order
- Phase 24: CrossValidationResult uses concrete types only (no type parameter) for Aqua compliance
- Phase 24: Two-method dispatch for cross_validate_gap (LindbladianResult delegates to Complex)
- Phase 24: abs(real(spectral_gap)) enforced in cross_validate_gap (locked decision)
- Phase 24: Excited initial state (psi0[end]=1) for validation -- ground state at high beta is near Gibbs, no decay signal
- Phase 24: Two-tier pass criterion (R-squared > 0.9 AND residual_factor in [0.8, 1.5]) tightened in quick-23
- Quick-22: Fixed delta_eff double-counting -- trajectory CPTP channel now uses bare delta (matching DM), R_a scaled by n_jumps is the single compensation
- Quick-22: Residual factor ~1.6x (n=4) and ~1.5x (n=6) between fitted and exact gap is discrete-step Kraus effect (was ~20x before fix)
- Quick-23: Smallest-gap selection reduces n=4 factor from ~1.6x to ~1.17x; n=6 at ~1.46x; both pass [0.8, 1.5]
- Quick-24: Added XX_avg, YY_avg, ZZ_avg per-bond averaged correlations to gap estimation (5 observables total)
- Quick-24: SingularException in LsqFit.stderror handled gracefully with Inf/(-Inf,Inf) fallback
- Phase 25-01: Single observable builder (build_preset_trajectory_observables) replaces 4 old builders
- Phase 25-01: Mz construction inlined into single builder (was delegated to deleted build_total_magnetization)
- Phase 25-01: CrossValidationResult and cross_validate_gap removed from source and exports
- Phase 25-02: ARPACK vs dense eigen gap agrees to ~1.8e-10 for 3-qubit system (test atol=1e-8)
- Phase 25-02: dot(O_vec, V[:,k]) computes vec(O)^H * v_k -- correct for Hermitian observable eigenbasis projection
- Phase 25-02: Relative gap overlap excludes steady-state mode (k=1) from denominator
- Phase 25-03: n=4 ARPACK vs eigen agrees to 1.2e-10 (well within 1e-8 threshold)
- Phase 25-03: n=4 gap estimation passes with 0.72% relative error (ZZ_avg best, 20k trajectories)
- Phase 25-03: n=6 all observables have zero gap-mode overlap -- estimation fails at 10.7% relative error
- Phase 25-03: /experiments/ is gitignored; validation script force-added to git
- Quick-25: n=6 gap mode in k=3 (k=pi) momentum sector; n=4 gap mode in k=0 -- translational symmetry momentum mismatch is root cause of zero overlap
- Quick-25: T_L = kron(T_eigen, conj(T_eigen)) commutes with L to ~1e-15 precision, confirming Lindbladian inherits translational symmetry
- Quick-25: n=6 gap eigenspace is 3-fold degenerate, all in k=3 sector; need k=pi observables or symmetry-breaking for viable estimation
- Quick-26: Added Mz_stagg (staggered Z) and Z1 (single-site) to 7-observable bundle; Mz_stagg has |c_gap|=0 for n=6 despite k=pi component
- Quick-26: n=6 gap mode protected by additional SU(2) spin-rotation symmetry -- staggered-Z alone insufficient
- Quick-26: Tiered validation thresholds: n=4 < 1%, n=6 < 12% (acknowledges symmetry-protected gap mode)

### Pending Todos

None

### Blockers/Concerns

None

### Roadmap Evolution

- Phase 25 added: Spectral Gap Validation Overhaul — consolidate validation, verify ARPACK vs eigen, eigenbasis overlap analysis, high-ntraj estimation

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

## Session Continuity

Last session: 2026-02-18
Stopped at: Completed quick task 26 (staggered observables). 7-observable bundle in place; n=6 needs symmetry-breaking beyond k=pi.
Resume file: None
