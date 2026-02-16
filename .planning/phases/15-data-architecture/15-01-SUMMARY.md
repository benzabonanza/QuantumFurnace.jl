---
phase: 15-data-architecture
plan: 01
subsystem: data-persistence
tags: [bson, serialization, experiment-results, libgit2, metadata]

# Dependency graph
requires:
  - phase: 14-gns-trajectory
    provides: TrajectoryResult struct, run_trajectories returning trajectory-averaged density matrices
provides:
  - ExperimentResult{C,T} struct bundling config, trajectory result, hamiltonian params, metadata
  - save_experiment / load_experiment for BSON persistence with Dict-based serialization
  - Companion .txt file for human-readable browsing without Julia
  - Metadata auto-capture (git hash, Julia version, timestamp, thread count, wall time)
  - Hamiltonian parameter extraction for provenance
  - Auto-generated filenames and default results directory structure
affects: [16-convergence-tracking, 17-adaptive-sampling, 18-experiments]

# Tech tracking
tech-stack:
  added: [LibGit2 (stdlib), Dates (stdlib)]
  patterns: [Dict-based BSON serialization, string-tagged domain/config type reconstruction, Hermitian unwrapping]

key-files:
  created: [src/results.jl]
  modified: [src/QuantumFurnace.jl, Project.toml]

key-decisions:
  - "Dict-based BSON serialization avoids parametric struct and abstract type pitfalls"
  - "TrajectoryResult embedded as field (not flattened) to avoid name collisions"
  - "Domain singletons stored as strings with DOMAIN_LOOKUP for reconstruction"
  - "gaussian_parameters tuple converted back from BSON array on load"
  - "LibGit2 for git hash capture (no external dependency, works without git on PATH)"
  - "Hamiltonian params store only reproduction-relevant subset (no eigendecomposition)"

patterns-established:
  - "Dict-based BSON: convert to plain Dict before save, reconstruct from Dict on load"
  - "String-tagged types: config_type/config_kind/domain as strings for forward compatibility"
  - "Forward-compatible loading: get() with defaults for missing fields"

# Metrics
duration: 4min
completed: 2026-02-16
---

# Phase 15 Plan 01: ExperimentResult Summary

**ExperimentResult{C,T} struct with Dict-based BSON save/load, companion .txt, metadata auto-capture via LibGit2, and auto-generated filename/directory structure**

## Performance

- **Duration:** 4 min
- **Started:** 2026-02-16T09:10:36Z
- **Completed:** 2026-02-16T09:15:03Z
- **Tasks:** 2
- **Files modified:** 3

## Accomplishments
- ExperimentResult{C,T} struct with config, trajectory_result, hamiltonian_params, and metadata fields
- Full Dict-based BSON serialization avoiding all parametric struct and abstract type pitfalls
- save_experiment writes BSON + companion .txt; supports explicit path or auto-generated path
- load_experiment reconstructs ExperimentResult from BSON with forward-compatible field handling
- Metadata auto-captures git hash (LibGit2), Julia version, timestamp, thread count, wall time
- All four config types (LiouvConfig, LiouvConfigGNS, ThermalizeConfig, ThermalizeConfigGNS) serialize and reconstruct via string-tagged type dispatch
- Package compiles cleanly; all 284 existing tests pass

## Task Commits

Each task was committed atomically:

1. **Task 1: ExperimentResult struct and Dict conversion helpers** - `314fe72` (feat)
2. **Task 2: Save/load wrappers, companion txt, metadata capture, and exports** - `57480af` (feat)

## Files Created/Modified
- `src/results.jl` - ExperimentResult struct, Dict conversion, save/load, companion txt, metadata capture, filename generation
- `src/QuantumFurnace.jl` - Added using LibGit2/Dates, include("results.jl"), exports for ExperimentResult/save_experiment/load_experiment
- `Project.toml` - Added LibGit2 and Dates stdlib dependencies

## Decisions Made
- Dict-based BSON serialization: all fields converted to plain Dict of primitive types before BSON.bson, avoiding BSON.jl's issues with parametric structs and abstract type fields
- TrajectoryResult embedded as nested sub-dict under :trajectory key (not flattened) to avoid name collisions with metadata fields
- Domain singletons stored as strings ("TrotterDomain", etc.) with DOMAIN_LOOKUP constant for reconstruction
- Config type preserved via string tags: config_type ("KMS"/"GNS") + config_kind ("liouv"/"thermalize")
- gaussian_parameters tuple explicitly converted back from BSON array representation on load
- Hamiltonian params store only reproduction-relevant subset: base_coeffs, base_terms, disordering_term/coeffs, periodic, shift, rescaling_factor, num_qubits -- not eigendecomposition, bohr_freqs, bohr_dict, or gibbs

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Added LibGit2 and Dates to Project.toml**
- **Found during:** Task 2 (verification step - package compilation)
- **Issue:** LibGit2 and Dates are Julia stdlibs but still need explicit Project.toml entries for package dependencies
- **Fix:** Ran `Pkg.add("LibGit2")` and `Pkg.add("Dates")` to register them in Project.toml
- **Files modified:** Project.toml
- **Verification:** Package compiles cleanly after adding dependencies
- **Committed in:** 57480af (Task 2 commit)

**2. [Rule 1 - Bug] BSON tuple-to-array conversion for gaussian_parameters**
- **Found during:** Task 1 (Dict reconstruction implementation)
- **Issue:** BSON.jl stores Julia tuples as arrays; `(nothing, nothing)` or `(1.0, 2.0)` would come back as `[nothing, nothing]` or `[1.0, 2.0]`, causing constructor type mismatch
- **Fix:** Added explicit array-to-tuple conversion in `_dict_to_config_kwargs` for the gaussian_parameters field
- **Files modified:** src/results.jl
- **Verification:** Covered by forward-compatible loading pattern
- **Committed in:** 57480af (Task 2 commit, though logic was refined during Task 2)

---

**Total deviations:** 2 auto-fixed (1 blocking, 1 bug)
**Impact on plan:** Both auto-fixes necessary for correctness. No scope creep.

## Issues Encountered
None beyond the auto-fixed deviations above.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- ExperimentResult infrastructure is complete and ready for use by Phase 16 (convergence tracking) and Phase 17 (adaptive sampling)
- Phase 18 (experiments) can use save_experiment/load_experiment for persisting sweep results
- The results/ directory structure (kms/, approx_gns/) is created on first save

## Self-Check: PASSED

- FOUND: src/results.jl
- FOUND: .planning/phases/15-data-architecture/15-01-SUMMARY.md
- FOUND: 314fe72 (Task 1 commit)
- FOUND: 57480af (Task 2 commit)

---
*Phase: 15-data-architecture*
*Completed: 2026-02-16*
