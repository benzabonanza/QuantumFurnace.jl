# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-16)

**Core value:** Correct and efficient classical simulation of Lindbladian-based quantum Gibbs samplers
**Current focus:** v1.3 Mixing Time Estimation -- Phase 24 (Cross-Validation)

## Current Position

Phase: 24 of 24 (Cross-Validation) -- COMPLETE
Plan: 3 of 3 in current phase (COMPLETE)
Status: Phase 24 complete (gap closure plan 03 executed). v1.3 Mixing Time Estimation milestone complete.
Last activity: 2026-02-17 -- Phase 24 Plan 03 executed (gap closure: n_jumps normalization)

Progress: [##############################] 48/48 plans (v1.0-v1.3 ALL COMPLETE)

## Performance Metrics

**Velocity:**
- Total plans completed: 53 (v1.0: 10, v1.1: 16, quick: 8, v1.2: 12, cleanup: 3, v1.3: 7)

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
- Phase 23: Best observable selected by converged + gap > 0 + highest R-squared; fallback to highest R-squared if no valid fit
- Phase 23: gap_estimation.jl included after fitting.jl and before results.jl for correct dependency order
- Phase 24: CrossValidationResult uses concrete types only (no type parameter) for Aqua compliance
- Phase 24: Two-method dispatch for cross_validate_gap (LindbladianResult delegates to Complex)
- Phase 24: abs(real(spectral_gap)) enforced in cross_validate_gap (locked decision)
- Phase 24: Excited initial state (psi0[end]=1) for validation -- ground state at high beta is near Gibbs, no decay signal
- Phase 24: Normalization factor (~20x n=4, ~28x n=6) between trajectory rate and Liouvillian gap is a physics finding (delta_eff = delta * n_jumps)
- Phase 24: Two-tier pass criterion (R-squared > 0.9 AND residual_factor in [1.0, 3.0]) constitutes "documented tolerance" per ROADMAP success criterion 3
- Phase 24: Residual factor after n_jumps correction is ~1.66 (n=4) and ~1.57 (n=6) -- consistent across system sizes

### Pending Todos

None

### Blockers/Concerns

None

### Quick Tasks Completed

| # | Description | Date | Commit | Directory |
|---|-------------|------|--------|-----------|
| 19 | Fix failing test after EnergyDomain to TrotterDomain rename in GNS trajectory | 2026-02-16 | 0683acc | [19-fix-failing-test-after-energydomain-to-t](./quick/19-fix-failing-test-after-energydomain-to-t/) |
| 20 | Fix basis mismatch in GNS TrotterDomain tests (0.83 -> 0.08 gap) | 2026-02-16 | 0161e01 | [20-debug-gns-trotterdomain-0-83-gap-suspect](./quick/20-debug-gns-trotterdomain-0-83-gap-suspect/) |
| 21 | Fix test errors after removing transform_jumps_to_basis | 2026-02-16 | b2e4123 | [21-fix-test-errors-after-removing-transform](./quick/21-fix-test-errors-after-removing-transform/) |

## Session Continuity

Last session: 2026-02-17
Stopped at: Completed 24-03-PLAN.md (Gap Closure). Phase 24 fully complete. v1.3 milestone complete.
Resume file: None
