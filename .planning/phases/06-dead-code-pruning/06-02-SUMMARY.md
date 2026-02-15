---
phase: 06-dead-code-pruning
plan: 02
subsystem: codebase
tags: [dead-code, cleanup, julia, quantum-simulation, structs, unused-functions]

# Dependency graph
requires:
  - "06-01: Commented code removal (cleared commented blocks that would have confused dead code analysis)"
provides:
  - "3 dead structs removed (LindbladianJumpCaches, LiouvLiouv, LindbladWorkspace)"
  - "~35 unreachable functions removed across 5 source files"
  - "Non-mutating oft() wrapper removed, in-place oft!() preserved"
  - "Clean export list with no references to deleted symbols"
  - "All 224 tests pass with zero regressions"
affects: [07-unused-function-pruning, 08-dead-struct-pruning, 10-api-cleanup]

# Tech tracking
tech-stack:
  added: []
  patterns: []

key-files:
  created: []
  modified:
    - "src/structs.jl"
    - "src/ofts.jl"
    - "src/bohr_domain.jl"
    - "src/coherent.jl"
    - "src/trotter_domain.jl"
    - "src/qi_tools.jl"
    - "src/errors.jl"
    - "src/QuantumFurnace.jl"

key-decisions:
  - "Preserved all qi_tools.jl functions except are_we_tp (pedagogical/user-facing value per locked decision)"
  - "Left errors.jl as placeholder file for Phase 10 API cleanup"
  - "Preserved compute_b_minus, compute_b_plus, and support functions in coherent.jl (used by furnace_utensils.jl)"

patterns-established: []

# Metrics
duration: 8min
completed: 2026-02-15
---

# Phase 6 Plan 2: Dead Code Pruning - Unused Functions and Dead Structs Summary

**Removed 3 dead struct types, non-mutating oft() wrapper, and ~35 unreachable functions totaling ~497 lines; all 224 tests pass**

## Performance

- **Duration:** 8 min
- **Started:** 2026-02-15T00:28:54Z
- **Completed:** 2026-02-15T00:37:13Z
- **Tasks:** 3
- **Files modified:** 8

## Accomplishments
- Removed 3 dead struct types: LindbladianJumpCaches, LiouvLiouv, LindbladWorkspace{T} (PRUNE-03)
- Removed non-mutating oft() wrapper; preserved in-place oft!(), time_oft!(), trotter_oft!()
- Removed ~35 unreachable functions across bohr_domain.jl (8), coherent.jl (9), trotter_domain.jl (4), qi_tools.jl (1), errors.jl (7) (PRUNE-02)
- Updated export list to remove LindbladWorkspace; all Aqua quality checks pass
- Preserved all functions on the keep list: qi_tools (minus are_we_tp), linearmaps_liouv, log_sobolev

## Task Commits

Each task was committed atomically:

1. **Task 1: Remove dead structs and non-mutating oft wrapper** - `b598002` (chore)
2. **Task 2: Remove unused active functions** - `44c09e9` (chore)
3. **Task 3: Update exports and run final verification** - `510f450` (chore)

## Files Modified

| File | Lines Removed | What was removed |
|------|-------------|------------------|
| `src/structs.jl` | 29 | LindbladianJumpCaches, LiouvLiouv, LindbladWorkspace{T} structs |
| `src/ofts.jl` | 4 | Non-mutating oft() wrapper |
| `src/bohr_domain.jl` | 162 | coherent_bohr_gauss, transition_bohr_gauss_vectorized, transition_bohr_gauss_gibbsed_vectorized, thermalize_bohr_gauss_vectorized, B_nu_gauss, R_nu_gauss, check_alpha_skew_symmetry, find_all_nu1s_to_nu2 |
| `src/coherent.jl` | 207 | coherent_term_time, coherent_term_trotter, coherent_term_time_metro_exact, coherent_term_timedomain_integrated_gauss, coherent_term_time_integrated_metro, coherent_term_time_integrated_eh, coherent_term_time_integrated_eh_b, check_B_gauss, check_B_metro |
| `src/trotter_domain.jl` | 52 | trotterize (1st order), trotter_diag, trotter2_diag, trotter2_t0_multiple |
| `src/qi_tools.jl` | 13 | are_we_tp (references undefined globals) |
| `src/errors.jl` | 41 | compute_errors, compute_quadrature_error, compute_energy_quadrature_error, 4 empty stubs; replaced with placeholder comment |
| `src/QuantumFurnace.jl` | 1 | LindbladWorkspace removed from export |

## Decisions Made
- Preserved all qi_tools.jl functions except are_we_tp per locked CONTEXT.md decision (pedagogical value)
- Left errors.jl as a 1-line placeholder rather than deleting the file, since Phase 10 API cleanup may add proper error types there
- Kept compute_b_minus, compute_b_plus, compute_b_plus_metro, compute_b_plus_smooth, compute_truncated_func, get_truncated_indices, and convolute in coherent.jl -- they are called from furnace_utensils.jl (active code), not dead

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
- Aqua "Undefined exports" test failed after Task 2 (before Task 3 export cleanup) because LindbladWorkspace was still exported but the struct was removed in Task 1. Resolved in Task 3 as planned.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- Phase 6 (Dead Code Pruning) is fully complete: PRUNE-01 (commented code), PRUNE-02 (unreachable functions), PRUNE-03 (dead structs) all satisfied
- Codebase is clean and ready for subsequent v1.1 phases (7-11)
- No blockers or concerns

## Self-Check: PASSED

- All 8 modified source files exist on disk
- All 3 task commits verified in git log (b598002, 44c09e9, 510f450)
- 224/224 tests pass

---
*Phase: 06-dead-code-pruning*
*Completed: 2026-02-15*
