---
phase: 25-quick
plan: 25
subsystem: spectral-analysis
tags: [lindbladian, momentum-sectors, translational-symmetry, heisenberg-chain, spectral-gap]

# Dependency graph
requires:
  - phase: 25-03
    provides: "Unified validation script identifying n=6 zero gap-mode overlap problem"
provides:
  - "Definitive diagnosis: n=6 gap mode in k=3 (pi) momentum sector, not k=0"
  - "Translation operator construction in Hilbert and Liouville space"
  - "Momentum sector measurement via Rayleigh quotient of T_L eigenvectors"
affects: [spectral-gap-estimation, observable-design, symmetry-breaking]

# Tech tracking
tech-stack:
  added: []
  patterns: [translation-operator-construction, liouville-space-momentum-analysis, rayleigh-quotient-momentum-measurement]

key-files:
  created: [experiments/diagnose_gap_momentum.jl]
  modified: []

key-decisions:
  - "n=6 gap mode confirmed in k=3 (pi) sector -- translational symmetry momentum mismatch explains zero overlap"
  - "n=4 gap mode confirmed in k=0 sector -- same sector as observables, overlap is nonzero"
  - "T_L = kron(T_eigen, conj(T_eigen)) in eigenbasis correctly commutes with L (norm < 1e-15)"
  - "All 5 observables (H, Mz, XX_avg, YY_avg, ZZ_avg) confirmed k=0 for both n=4 and n=6"
  - "n=6 has 3-fold degenerate gap eigenspace, all in k=3 sector"

patterns-established:
  - "Translation operator: cyclic left-shift of n-bit computational basis index"
  - "Liouville space symmetry: T_L = kron(T_eigen, conj(T_eigen)) for Watrous convention"
  - "Momentum measurement: Rayleigh quotient dot(v, T_L*v)/dot(v,v) gives exp(ik)"

# Metrics
duration: 4min
completed: 2026-02-18
---

# Quick Task 25: Diagnose n=6 Gap-Mode Momentum Sector Summary

**Confirmed n=6 Lindbladian gap mode lives in k=pi (m=3) momentum sector while all observables are k=0, explaining zero overlap by translational symmetry orthogonality**

## Performance

- **Duration:** 4 min
- **Started:** 2026-02-18T09:09:36Z
- **Completed:** 2026-02-18T09:14:02Z
- **Tasks:** 1
- **Files created:** 1

## Accomplishments
- Definitively confirmed the momentum sector hypothesis for the n=6 zero-overlap phenomenon
- Built translation operator T in Hilbert space (cyclic qubit permutation), verified T^n = I and [T,H]=0
- Constructed T_L in Liouville eigenbasis via kron(T_eigen, conj(T_eigen)), verified [T_L, L]/||L|| < 1e-15
- For n=4: gap mode is in k=0 sector (same as observables) -- overlap nonzero, estimation works
- For n=6: gap mode is in k=3 sector (k=pi, T_L eigenvalue = -1) -- orthogonal to all k=0 observables
- All 5 observables (H, Mz, XX_avg, YY_avg, ZZ_avg) confirmed in k=0 sector for both system sizes
- n=6 gap eigenspace is 3-fold degenerate, all three modes in k=3 sector

## Key Results

| Property | n=4 | n=6 |
|----------|-----|-----|
| ||[T_L, L]|| / ||L|| | 1.2e-15 | 3.1e-15 |
| Gap eigenvalue | 0.1585339428 | 0.1156098797 |
| Gap mode T_L eigenvalue | +1.0 | -1.0 |
| Gap mode momentum sector | m=0 (k=0) | m=3 (k=pi) |
| Gap mode degeneracy | 1 | 3 |
| All observables k=0 | Yes | Yes |
| Overlap possible | Yes | No (orthogonal sectors) |

## Task Commits

1. **Task 1: Create momentum sector diagnostic script** - `cf61f88` (feat)

## Files Created
- `experiments/diagnose_gap_momentum.jl` - Momentum sector diagnostic script (~230 lines) that constructs translation operator, measures gap mode momentum via Rayleigh quotient, verifies observable sectors, and reports conclusion

## Decisions Made
- Used Rayleigh quotient dot(v, T_L*v)/dot(v,v) to extract T_L eigenvalue on each eigenvector (avoids full eigendecomposition of T_L)
- Reports momentum for first 10 eigenmodes and identifies near-degenerate modes within 10% of gap eigenvalue
- Confirmed that n=6 next gap in k=0 sector would need to be identified for viable estimation (modes 5-10 at Re(lam) ~ -0.136 have mixed momentum structure)

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None.

## Implications for Future Work
- To estimate spectral gap for n=6, need observables in the k=3 (k=pi) sector
- Possible approach: staggered observables like (-1)^i Z_i which have k=pi translational eigenvalue
- Alternative: break translational symmetry with disorder (disordered Heisenberg chain)
- The 3-fold degeneracy of the n=6 gap eigenspace in k=3 may relate to additional symmetries (spin rotation, parity)

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- Root cause of n=6 zero overlap definitively identified
- Clear path forward: either use k=pi observables or break translational symmetry
- No blockers

## Self-Check: PASSED

- [x] experiments/diagnose_gap_momentum.jl exists (308 lines, min 80 required)
- [x] Commit cf61f88 exists in git log
- [x] Script runs successfully for both n=4 and n=6
- [x] T^n = I verified for both system sizes
- [x] [T_L, L]/||L|| < 1e-10 for both system sizes
- [x] Gap mode momentum reported for both n=4 (m=0) and n=6 (m=3)
- [x] All 5 observables confirmed k=0 for both system sizes
- [x] Conclusion stated: CONFIRMED for n=6, REFUTED for n=4

---
*Quick Task: 25-diagnose-n-6-gap-mode-momentum-sector*
*Completed: 2026-02-18*
