---
phase: 06-dead-code-pruning
plan: 01
subsystem: codebase
tags: [dead-code, cleanup, julia, quantum-simulation]

# Dependency graph
requires: []
provides:
  - "Clean source files with ~987 lines of commented-out code removed across 10 files"
  - "All TODO/FIXME/HACK annotations preserved"
  - "All 224 existing tests pass unchanged"
affects: [07-unused-function-pruning, 08-dead-struct-pruning]

# Tech tracking
tech-stack:
  added: []
  patterns: []

key-files:
  created: []
  modified:
    - "src/qi_tools.jl"
    - "src/energy_domain.jl"
    - "src/trotter_domain.jl"
    - "src/coherent.jl"
    - "src/trajectories.jl"
    - "src/jump_workers.jl"
    - "src/errors.jl"
    - "src/misc_tools.jl"
    - "src/ofts.jl"
    - "src/bohr_domain.jl"

key-decisions:
  - "Preserved all inline single-line commented-out debug assertions within active functions (e.g., commented @assert, @printf)"
  - "Preserved 4-line commented alternative implementation in coherent.jl (non-Hermitian case) as design documentation"
  - "linearmaps_liouv.jl excluded from scope per CONTEXT.md explicit keep list"

patterns-established: []

# Metrics
duration: 8min
completed: 2026-02-15
---

# Phase 6 Plan 1: Dead Code Pruning - Commented Code Removal Summary

**Removed 987 lines of commented-out dead code from 10 Julia source files; all 224 tests pass with zero regressions**

## Performance

- **Duration:** 8 min
- **Started:** 2026-02-15T00:17:53Z
- **Completed:** 2026-02-15T00:26:00Z
- **Tasks:** 2 (1 execution + 1 verification)
- **Files modified:** 10

## Accomplishments
- Removed ~987 lines of commented-out code blocks across 10 source files (original total: 4087 lines, new total: 3100 lines)
- Preserved all TODO/FIXME/HACK comments (4 across codebase)
- Preserved all explanatory comments describing design decisions and rationale
- Verified all 224 tests pass with no regressions via `Pkg.test()`

## Task Commits

Each task was committed atomically:

1. **Task 1: Remove commented-out code blocks from all source files** - `3f23aa5` (chore)
2. **Task 2: Verify no commented code blocks remain and run full test suite** - (verification only, no file changes needed)

## Files Modified

| File | Lines Removed | What was removed |
|------|-------------|------------------|
| `src/qi_tools.jl` | 47 | Old QuantumOptics benchmark/test comparison code |
| `src/energy_domain.jl` | 140 | 8 old commented functions (transition_gauss_vectorized, create_alpha variants, integrate_gamma variants) |
| `src/trotter_domain.jl` | 165 | Pauli exponentiation tests, Hamiltonian time evolution tests, Trotter parameter sweeps, Fourier label tests |
| `src/coherent.jl` | 209 | Old compute_f_minus/f_plus/f_plus_metro/f_plus_eh functions + large test/benchmark block |
| `src/trajectories.jl` | 166 | Old evolve_and_measure_along_trajectory, measure!, KrausFramework, krausframework, construct_gksl_lindbladian, old evolve/step functions |
| `src/jump_workers.jl` | 111 | Old verify_completeness, precompute_B, apply_lindbladian!, apply_lindbladian_dagger! |
| `src/errors.jl` | 50 | Commented if-trotter block + quadrature error test code |
| `src/misc_tools.jl` | 38 | Old validate_config! implementation (pre-GNS) |
| `src/ofts.jl` | 23 | Old time_oft_integrated + Trotter OFT check script |
| `src/bohr_domain.jl` | 38 | Old exact_mask, approx_mask, transition_bohr_metro_gibbsed |

## Decisions Made
- Preserved all inline single-line commented-out debug assertions within active functions (e.g., commented `@assert`, `@printf` calls that serve as togglable debug aids)
- Preserved 4-line commented alternative implementation in `coherent.jl` (non-Hermitian B case at line 299) as design documentation for a known mathematical variant
- `linearmaps_liouv.jl` was excluded from this plan's scope entirely per CONTEXT.md explicit keep list ("future project dependency, preserve")

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
- `Aqua` test dependency was not installed in the environment; resolved by running `Pkg.test()` which handles test-only [extras] dependencies automatically. No changes to Project.toml/Manifest.toml were committed.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- Codebase is clean of commented-out dead code blocks
- Ready for 06-02 (unused function pruning) or subsequent phases
- No blockers or concerns

---
*Phase: 06-dead-code-pruning*
*Completed: 2026-02-15*
