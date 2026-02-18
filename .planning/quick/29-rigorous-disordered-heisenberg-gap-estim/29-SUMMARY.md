---
phase: quick-29
plan: 01
subsystem: experiments
tags: [spectral-gap, disorder, xx-stagg, heisenberg, lindbladian, eigenbasis-overlap, observable-design]

# Dependency graph
requires:
  - phase: quick-28
    provides: Disordered Heisenberg validation; Mz_stagg |c_gap|=0.120 for n=6; gap estimation biased at 34%
provides:
  - "XX_stagg observable validation: zero gap-mode overlap for disordered Heisenberg (both n=4 and n=6)"
  - "Evidence that XX two-site correlation structure does NOT couple to the Lindbladian gap mode"
  - "Focused 2-observable (H + XX_stagg) validation script with deep diagnostics"
affects: [gap-estimation, observable-design, symmetry-breaking]

# Tech tracking
tech-stack:
  added: []
  patterns: [manual observable construction with pad_term for targeted eigenbasis analysis]

key-files:
  created: [experiments/validate_gap_xx_stagg.jl]
  modified: []

key-decisions:
  - "XX_stagg (staggered nearest-neighbor XX) has zero gap-mode overlap even with disorder: |c_gap|=2.5e-5 (n=4), 3.5e-6 (n=6)"
  - "XX_stagg's dominant expansion modes (k~42-43) are ~2x the gap rate, meaning its decay is dominated by faster modes, not the gap"
  - "H alone at 1.21x (n=4) and 1.63x (n=6) gap/exact ratio -- discrete-step Kraus overestimation confirmed"
  - "XX two-site correlation structure is fundamentally misaligned with the gap mode even when symmetries are broken"

patterns-established:
  - "Focused 2-observable validation pattern: isolate H + single new observable for clean diagnostic"
  - "run_observable_trajectories called separately from estimate_spectral_gap for raw time series access"

# Metrics
duration: 19min
completed: 2026-02-18
---

# Quick Task 29: XX_stagg Disordered Heisenberg Gap Validation Summary

**XX_stagg (staggered nearest-neighbor XX correlation) has zero gap-mode overlap for disordered Heisenberg: |c_gap|=2.5e-5 (n=4), 3.5e-6 (n=6); XX two-site correlations do not couple to the Lindbladian gap mode**

## Performance

- **Duration:** 19 min
- **Started:** 2026-02-18T16:15:47Z
- **Completed:** 2026-02-18T16:35:00Z
- **Tasks:** 2
- **Files modified:** 1

## Accomplishments

- CONFIRMED: XX_stagg has essentially zero gap-mode overlap for both n=4 (|c_gap|=2.5e-5) and n=6 (|c_gap|=3.5e-6) disordered Heisenberg chains
- Full expansion coefficient spectrum reveals XX_stagg's dominant modes at k~42-43 (2x gap rate), explaining why it cannot estimate the gap
- H alone provides gap estimation at 1.21x (n=4) and 1.63x (n=6) -- consistent with discrete-step Kraus overestimation pattern
- Deep diagnostics confirm XX_stagg time series is pure noise (dynamic range < 0.002), no exponential decay signal

## Experimental Results

### Eigenbasis Overlap Analysis

**n=4 (disordered):**
| Observable | |c_gap| | Relative overlap |
|------------|---------|-----------------|
| H | 0.137050 | 0.2946 |
| XX_stagg | 0.000025 | 0.0004 |

**n=6 (disordered):**
| Observable | |c_gap| | Relative overlap |
|------------|---------|-----------------|
| H | 0.005729 | 0.0163 |
| XX_stagg | 0.000003 | 0.0001 |

### Gap Estimation Results

| System | Exact Gap | Estimated Gap | Best Obs | Rel Error | H gap/exact | XX_stagg gap/exact |
|--------|-----------|--------------|----------|-----------|-------------|-------------------|
| n=4 disordered | 0.17298 | 0.20995 | H | 21.4% | 1.2137 | 0.9723 (R2=0.02) |
| n=6 disordered | 0.11343 | 0.18493 | H | 63.0% | 1.6304 | 0.0283 (R2=0.38) |

### XX_stagg Expansion Coefficient Spectrum

**Top 5 modes by |c_k| for XX_stagg (n=4):**
| Mode | |c_k| | Re(lambda) | Rate/Gap |
|------|-------|------------|----------|
| k=42 | 0.01176 | -0.4015 | 2.32x |
| k=145 | 0.00938 | -0.5702 | 3.30x |
| k=39 | 0.00631 | -0.3996 | 2.31x |
| k=148 | 0.00539 | -0.5752 | 3.32x |
| k=139 | 0.00308 | -0.5643 | 3.26x |

**Top 5 modes by |c_k| for XX_stagg (n=6):**
| Mode | |c_k| | Re(lambda) | Rate/Gap |
|------|-------|------------|----------|
| k=43 | 0.00523 | -0.2435 | 2.15x |
| k=111 | 0.00385 | -0.2901 | 2.56x |
| k=105 | 0.00325 | -0.2882 | 2.54x |
| k=44 | 0.00319 | -0.2440 | 2.15x |
| k=100 | 0.00252 | -0.2875 | 2.53x |

### Eigenvalue Near-Degeneracy

**n=4:** Modes 2-4 are near-degenerate (modes 3,4 separated by 3.9e-16). Gap mode (k=2) at Re(lambda)=-0.1730 is close to modes 3,4 at Re(lambda)=-0.1815.

**n=6:** Modes 2-4 are near-degenerate (modes 3,4 separated by 6.5e-16). Gap mode (k=2) at Re(lambda)=-0.1134 is very close to modes 3,4 at Re(lambda)=-0.1136.

### Raw Time Series

**XX_stagg (n=4):** Values oscillate around zero with tiny amplitude (~0.001). No exponential decay pattern visible. Initial value: 0.0 (XX_stagg of computational basis excited state), confirming the observable starts from zero and has no signal.

**XX_stagg (n=6):** Similarly flat near zero with amplitude ~0.001. No decay signal.

**H (both sizes):** Clear exponential decay from excited state energy (~0.45) to Gibbs equilibrium value, confirming H has a real signal while XX_stagg does not.

### Key Scientific Findings

1. **XX_stagg does NOT couple to the gap mode.** Despite having k=pi momentum (staggered sign) and two-site correlation structure, the XX operator pair does not project onto the gap eigenmode of the Lindbladian. This is in contrast to Mz_stagg which achieves |c_gap|=0.120 for n=6 with disorder.

2. **The gap mode couples to Z-type operators, not XX-type.** Comparing Quick-28 results (Mz_stagg |c_gap|=0.120, Z1 |c_gap|=0.119) with this result (XX_stagg |c_gap|=3.5e-6), the gap mode has strong Z-character but essentially zero XX-character. This suggests the gap eigenmode has specific operator-type selectivity beyond momentum considerations.

3. **XX_stagg's expansion is dominated by fast modes.** The dominant modes (k~42-43) have decay rates 2-3x the gap rate. This means any fitted decay rate from XX_stagg would reflect these faster modes, not the spectral gap.

4. **H overestimates due to discrete-step Kraus effect.** With only H + XX_stagg, the selected observable is always H (since XX_stagg is pure noise). H gives gap/exact ratios of 1.21x (n=4) and 1.63x (n=6), consistent with the systematic overestimation seen in Quick-28.

5. **Near-degenerate eigenvalue clusters exist near the gap.** For both system sizes, modes 2-4 form tight clusters. This degeneracy structure may explain why single-observable fitting is biased -- the fitted exponential picks up contributions from the entire cluster, not just the gap mode.

## Task Commits

Each task was committed atomically:

1. **Task 1: Create XX_stagg validation script with deep diagnostics** - `4914103` (feat)
2. **Task 2: Run validation script and capture results** - no code changes (execution only)

## Files Created/Modified

- `experiments/validate_gap_xx_stagg.jl` - XX_stagg disordered Heisenberg gap validation experiment (499 lines)

## Decisions Made

- XX_stagg has zero gap-mode overlap even with disorder -- XX two-site correlation structure is fundamentally misaligned with the gap mode
- The gap mode has Z-character selectivity: Mz_stagg and Z1 couple (Quick-28), but XX_stagg does not
- H alone overestimates (1.21x-1.63x) due to discrete-step Kraus effect, consistent with prior observations

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None - script ran successfully on first attempt for both n=4 and n=6.

## User Setup Required

None - no external service configuration required.

## Next Steps

- XX_stagg is NOT a viable gap estimation observable for disordered Heisenberg chains
- The gap mode's Z-character selectivity suggests ZZ-type staggered observables (e.g., ZZ_stagg = Sum((-1)^i * Z_i * Z_{i+1})/n) might have stronger coupling than XX-type
- Focus should return to improving the selection algorithm (smallest-gap picks underestimating Mz_stagg) or finding a debiasing approach for H's overestimation

## Self-Check: PASSED

- FOUND: experiments/validate_gap_xx_stagg.jl (499 lines)
- FOUND: commit 4914103
- FOUND: 29-SUMMARY.md

---
*Quick Task: 29-rigorous-disordered-heisenberg-gap-estim*
*Completed: 2026-02-18*
