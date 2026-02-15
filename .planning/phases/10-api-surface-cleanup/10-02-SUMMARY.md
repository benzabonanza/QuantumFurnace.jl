---
phase: 10-api-surface-cleanup
plan: 02
subsystem: api
tags: [internal-functions, underscore-prefix, naming-convention, julia-module]

# Dependency graph
requires:
  - phase: 10-api-surface-cleanup
    plan: 01
    provides: "Curated export block with labeled groups; internal functions identified for _ prefix"
provides:
  - "All ~45 internal function definitions prefixed with _ across 14 source files"
  - "All intra-file call sites updated to use _ prefixed names"
  - "Exported functions unchanged"
  - "Types unchanged"
affects: [10-03]

# Tech tracking
tech-stack:
  added: []
  patterns: ["Internal functions prefixed with _ to signal non-public API"]

key-files:
  created: []
  modified:
    - src/hamiltonian.jl
    - src/trotter_domain.jl
    - src/qi_tools.jl
    - src/time_domain.jl
    - src/nufft.jl
    - src/ofts.jl
    - src/energy_domain.jl
    - src/bohr_domain.jl
    - src/coherent.jl
    - src/log_sobolev.jl
    - src/jump_workers.jl
    - src/trajectories.jl
    - src/furnace_utensils.jl
    - src/misc_tools.jl

key-decisions:
  - "pick_transition kept without _ prefix (exported function, despite plan suggesting rename)"

patterns-established:
  - "All non-exported function definitions use _ prefix convention"
  - "Cross-file calls to old names remain temporarily broken until Plan 10-03"

# Metrics
duration: 7min
completed: 2026-02-15
---

# Phase 10 Plan 02: Internal Function Prefix Rename Summary

**Renamed ~45 internal function definitions to _ prefix across 14 source files, updated all intra-file call sites**

## Performance

- **Duration:** 7 min
- **Started:** 2026-02-15T13:20:00Z
- **Completed:** 2026-02-15T13:27:23Z
- **Tasks:** 2
- **Files modified:** 14

## Accomplishments
- All internal (non-exported) function definitions prefixed with `_` across 14 source files
- All intra-file call sites updated to match new names (definitions + calls within same file)
- All exported functions retain original names (trace_distance_h, create_f, trotterize, pick_transition, etc.)
- All types retain original names (HamHam, TrottTrott, NUFFTPrefactors, etc.)
- Module loads without errors (`julia --project -e 'using QuantumFurnace'` succeeds)

## Task Commits

Each task was committed atomically:

1. **Task 1: Rename in 10 self-contained files** - `f0739b7` (feat)
2. **Task 2: Rename in 4 worker/integration files** - `05734d6` (feat)

## Files Created/Modified
- `src/hamiltonian.jl` - _construct_base_ham, _construct_disordering_terms, _rescaling_and_shift_factors
- `src/trotter_domain.jl` - _trotter_bohr_freqs, _trotterize2, _does_term_differ_at_both_sites, _compute_U_group
- `src/qi_tools.jl` - _kron!, _vectorize_liouv_diss_and_add!, _vectorize_liouvillian_coherent!
- `src/time_domain.jl` - _truncate_time_labels_for_oft
- `src/nufft.jl` - _prepare_oft_nufft_prefactors, _prefactor_view
- `src/ofts.jl` - _time_oft!, _trotter_oft!
- `src/energy_domain.jl` - _create_energy_labels, _truncate_energy_labels
- `src/bohr_domain.jl` - _pick_f, _pick_alpha (4 dispatch variants)
- `src/coherent.jl` - _precompute_coherent_total_B, _precompute_coherent_unitary_terms, _precompute_coherent_terms, _compute_b_minus, _compute_b_plus, _compute_b_plus_metro, _compute_b_plus_smooth, _compute_truncated_func, _get_truncated_indices, _convolute
- `src/log_sobolev.jl` - _sandwich!
- `src/jump_workers.jl` - _jump_contribution! (6 dispatch variants), _apply_coherent_unitary!, _finalize_kraus_step!
- `src/trajectories.jl` - _precompute_R (2 dispatch variants)
- `src/furnace_utensils.jl` - _precompute_labels (2 variants), _precompute_data (4 variants), _select_b_plus_calculator
- `src/misc_tools.jl` - _generate_filename (2 variants), _riemann_sum (3 variants), _print_press (2 variants)

## Decisions Made
- `pick_transition` was NOT renamed despite the plan suggesting it, because it is in the export list (exported as a public physics building block function). The must_haves rule "All exported functions retain their original names" takes precedence.

## Deviations from Plan

None - plan executed as written (with the single correction of keeping `pick_transition` un-renamed since it is exported).

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- All internal function definitions now use _ prefix
- Cross-file calls still reference old names (e.g., `vectorize_liouv_diss_and_add!` called from jump_workers.jl, `precompute_data` called from furnace.jl/trajectories.jl)
- Plan 10-03 will update all cross-file call sites to use the new _ prefixed names

## Self-Check: PASSED

- SUMMARY.md exists
- Commit f0739b7 exists (Task 1)
- Commit 05734d6 exists (Task 2)
- Module loads without errors
- All non-exported function definitions have _ prefix
- All exported functions retain original names

---
*Phase: 10-api-surface-cleanup*
*Completed: 2026-02-15*
