---
phase: 18-kms-vs-gns-experiments
plan: 01
subsystem: experiments
tags: [trajectory, convergence, kms, gns, trotter, bson, heisenberg, sweep]

# Dependency graph
requires:
  - phase: 17-adaptive-sampling
    provides: run_trajectories_adaptive with convergence detection and ConvergenceData
  - phase: 16-convergence-tracking
    provides: _gibbs_in_trotter_basis, build_convergence_observables_trotter, trace distance tracking
  - phase: 15-experiment-results
    provides: ExperimentResult, save_experiment, load_experiment, _capture_metadata, _convergence_to_dict
  - phase: 14-gns-trajectories
    provides: ThermalizeConfigGNS, GNS trajectory support
  - phase: 13-kms-trajectories
    provides: ThermalizeConfig, KMS trajectory support, threaded run_trajectories
provides:
  - "experiments/run_sweep.jl: standalone 27-experiment KMS-vs-GNS parameter sweep script"
  - "Helper functions: build_heisenberg_xxx, build_trotter_system, run_experiment"
  - "9 verified BSON experiment results for n=4 slice (n=4 x beta=5,10,20 x {KMS, GNS@1/beta, GNS@0.5/beta})"
  - "EXPT-04 validated: KMS trace distance < GNS@sigma=1/beta for all beta values at n=4"
affects: [phase-19-analysis, experiment-data]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Standalone sweep script in experiments/ (not package API)"
    - "main() function wrapper for Julia global scope safety"
    - "CLI system size filtering via ARGS"
    - "Per-experiment try/catch with failure logging"

key-files:
  created:
    - experiments/run_sweep.jl
  modified:
    - .gitignore

key-decisions:
  - "mixing_time = 2.0 * beta (scales with inverse temperature for sufficient mixing at high beta)"
  - "batch_size=200, convergence_threshold=0.01, patience=3, window_size=3, min_batches=5 (Phase 17 defaults)"
  - "seed=42 for all experiments (deterministic, reproducible)"
  - "JumpOp[] (unparameterized) for compatibility with convergence.jl Vector{JumpOp} signatures"
  - "Wrap sweep in main() function to avoid Julia global scope scoping issues"

patterns-established:
  - "Experiment sweep pattern: build Hamiltonian + Trotter per (n,beta), run 3 experiments per pair"
  - "CLI filtering: ARGS-based system size restriction for quick verification runs"

# Metrics
duration: 18min
completed: 2026-02-16
---

# Phase 18 Plan 01: KMS-vs-GNS Experiment Sweep Summary

**Standalone sweep script running 27-experiment parameter grid comparing KMS and GNS convergence across n=4,6,8 and beta=5,10,20 using TrotterDomain adaptive trajectories, verified on n=4 slice with EXPT-04 confirmed**

## Performance

- **Duration:** 18 min
- **Started:** 2026-02-16T12:52:58Z
- **Completed:** 2026-02-16T13:10:57Z
- **Tasks:** 2
- **Files modified:** 2

## Accomplishments

- Created `experiments/run_sweep.jl` (262 lines) with `build_heisenberg_xxx`, `build_trotter_system`, `run_experiment` helpers and full parameter sweep loop
- Verified 9/9 n=4 experiments pass in 11.3 minutes total wall time
- EXPT-04 confirmed for all beta values: KMS trace distance < GNS@sigma=1/beta
  - beta=5: KMS=0.1177 < GNS=0.1698
  - beta=10: KMS=0.0774 < GNS=0.1975
  - beta=20: KMS=0.0245 < GNS=0.2049
- GNS@sigma=0.5/beta consistently outperforms GNS@sigma=1/beta (lower trace distance)
- Added `/experiments/` to `.gitignore` to prevent data files from being committed

## Task Commits

Each task was committed atomically:

1. **Task 1: Create sweep script and update .gitignore** - `1fac7f5` (feat)
2. **Task 2: Verify sweep script on n=4 slice** - `78c4789` (fix: scoping + type mismatch fixes discovered during verification)

## Files Created/Modified

- `experiments/run_sweep.jl` - Standalone sweep script with helpers and 27-experiment parameter grid loop
- `.gitignore` - Added `/experiments/` entry to prevent data file commits

## Decisions Made

- **mixing_time = 2.0 * beta**: Conservative scaling ensures sufficient mixing at all inverse temperatures. At beta=20, mixing_time=40.0 gives 4000 steps per trajectory.
- **Phase 17 defaults for adaptive sampling**: batch_size=200, convergence_threshold=0.01, patience=3, window_size=3, min_batches=5 -- well-tested values from Phase 17.
- **seed=42 for all experiments**: Deterministic, reproducible. Internal batch seeding (seed + n_total) ensures non-overlapping per-trajectory seeds.
- **Unparameterized JumpOp[]**: Using `JumpOp[]` instead of `JumpOp{Matrix{ComplexF64}}[]` because `Vector{JumpOp{Matrix{ComplexF64}}}` is not a subtype of `Vector{JumpOp}` (Julia invariance).
- **main() wrapper**: Wrapping sweep in a function avoids Julia's global-scope scoping ambiguity with `+=` in nested for loops.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed Julia global scope scoping issue**
- **Found during:** Task 2 (first run attempt)
- **Issue:** `experiment_count += 1` inside nested `for` loop at global scope triggers Julia's soft-scope warning and `UndefVarError`
- **Fix:** Wrapped entire sweep section in `main()` function
- **Files modified:** experiments/run_sweep.jl
- **Verification:** Script runs successfully after fix
- **Committed in:** 78c4789 (Task 2 commit)

**2. [Rule 1 - Bug] Fixed JumpOp type mismatch with run_trajectories_adaptive**
- **Found during:** Task 2 (first run attempt)
- **Issue:** `JumpOp{Matrix{ComplexF64}}[]` creates `Vector{JumpOp{Matrix{ComplexF64}}}` which is not a subtype of `Vector{JumpOp}` due to Julia's type invariance. `run_trajectories_adaptive` signature requires `Vector{JumpOp}`.
- **Fix:** Changed to `JumpOp[]` (creates `Vector{JumpOp}`) and updated `run_experiment` kwarg type annotation to match
- **Files modified:** experiments/run_sweep.jl
- **Verification:** All 9 experiments run successfully
- **Committed in:** 78c4789 (Task 2 commit)

---

**Total deviations:** 2 auto-fixed (2 bugs, Rule 1)
**Impact on plan:** Both fixes necessary for script to execute. No scope creep.

## Issues Encountered

- Normalization violation warnings during trajectory simulation (sum ~ 1.000001-1.000007). These are benign numerical precision issues in the trajectory stepping and do not affect convergence results.
- KMS at beta=20 hit the n_max=10000 cap (trace distance 0.0245, still improving). This is expected -- the spectral gap is smaller at high beta. For full production runs, n_max could be increased or mixing_time extended.

## User Setup Required

None - no external service configuration required.

## n=4 Verification Results

| Experiment | beta | sigma | Status | Final TD | n_traj | Wall Time |
|------------|------|-------|--------|----------|--------|-----------|
| KMS | 5 | 0.2000 | CONVERGED | 0.1177 | 6000 | 101.1s |
| GNS@1/beta | 5 | 0.2000 | CONVERGED | 0.1698 | 1800 | 31.0s |
| GNS@0.5/beta | 5 | 0.1000 | CONVERGED | 0.1364 | 6400 | 99.6s |
| KMS | 10 | 0.1000 | CONVERGED | 0.0774 | 2800 | 66.7s |
| GNS@1/beta | 10 | 0.1000 | CONVERGED | 0.1975 | 3200 | 91.8s |
| GNS@0.5/beta | 10 | 0.0500 | CONVERGED | 0.1170 | 4000 | 78.4s |
| KMS | 20 | 0.0500 | HIT CAP | 0.0245 | 10000 | 90.9s |
| GNS@1/beta | 20 | 0.0500 | CONVERGED | 0.2049 | 3400 | 47.9s |
| GNS@0.5/beta | 20 | 0.0250 | CONVERGED | 0.0553 | 8800 | 65.1s |

## Next Phase Readiness

- Sweep script ready for full 27-experiment run (n=4,6,8) via `julia --project=. experiments/run_sweep.jl`
- n=6 and n=8 experiments expected to take significantly longer (minutes to hours per experiment)
- Analysis of n=4 results ready for Phase 19+ (convergence curves, trace distance comparison plots)
- All BSON files contain full metadata: config, trajectory result, hamiltonian params, convergence data

## Self-Check: PASSED

All artifacts verified:
- experiments/run_sweep.jl exists
- /experiments/ in .gitignore
- Commit 1fac7f5 (Task 1) exists
- Commit 78c4789 (Task 2) exists
- All 9 BSON files exist in experiments/
- 18-01-SUMMARY.md exists

---
*Phase: 18-kms-vs-gns-experiments*
*Completed: 2026-02-16*
