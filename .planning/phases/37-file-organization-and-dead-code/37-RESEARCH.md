# Phase 37: File Organization and Dead Code - Research

**Researched:** 2026-02-27
**Domain:** Julia module organization, dead code removal, export list cleanup
**Confidence:** HIGH

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

#### Staging area design
- Move trajectory gap estimation files to `src/staging/` subdirectory
- Files to stage: convergence.jl, fitting.jl, gap_estimation.jl, log_sobolev.jl
- **Critical:** Verify each file contains ONLY trajectory gap estimation code. If any functions are used by Krylov or Lindbladian spectral gap paths, those functions must stay in active src/
- Staging code is **excluded from the module** -- not included in QuantumFurnace.jl
- Related tests move to `test/staging/` and are **not run** by the regular test suite (dormant)

#### Dead code removal
- Remove `@distributed` code and `using Distributed` import from furnace.jl; SharedArrays stays
- Delete old entry points outright: `run_lindbladian()`, `run_thermalization()` -- no deprecation warnings
- Delete all backward-compat config type aliases (LiouvConfig, ThermalizeConfig, LiouvConfigGNS, ThermalizeConfigGNS) if any remain
- Delete any dead or backward-compat structs -- target struct landscape is:
  - `Config{S,D,C,T}` -- unified config
  - `Workspace{S,D,C,T}` -- unified workspace (4 type params, no scratch type param)
  - 4 scratch types (one per simulation kind)
  - 4 Result structs (LindbladResults, ThermalizeResults, KrylovSpectrumResults, TrajectoryResults)
- Conservative approach overall: keep explicitly useful test utilities like `time_oft!`, `trotter_oft!` even if not called in production code
- Clean break: no deprecation shims, no "one more release" grace period

#### Export list organization
- Group by simulation type: `# Lindbladian`, `# Thermalize`, `# Krylov`, `# Trajectory` comment blocks
- Shared types/utilities (Config, domain types, qi_tools, hamiltonian utilities) go in `# Common` section **at the bottom**
- Staging code exports: keep as **commented-out block** with `# STAGING:` prefix for reference when code is reactivated
- Diagnostics types exported from main QuantumFurnace module directly (not a separate submodule)
- Old/deleted function names removed from export list entirely

### Claude's Discretion
- Exact ordering of exports within each simulation type section
- Which comment/docstring references to old types/functions to clean up vs leave
- Whether to clean up stale comments referencing deleted code in non-staging files

### Deferred Ideas (OUT OF SCOPE)
None -- discussion stayed within phase scope

</user_constraints>

## Summary

Phase 37 is a cleanup phase: remove dead code, move trajectory gap estimation code to a staging area, and reorganize the module export list. The codebase is a Julia physics simulation package (`QuantumFurnace.jl`) with 28 source files in `src/`, one main module file, and a comprehensive test suite.

The most critical finding of this research is that **convergence.jl cannot be moved wholesale to staging**. It contains functions (`_run_trajectory_convergence`, `_run_trajectory_adaptive`, `run_trajectories_convergence`, `run_trajectories_adaptive`, `build_preset_trajectory_observables`, `_windowed_relative_change`, `_compute_gibbs_observable_values`, `_gibbs_in_trotter_basis`) that are actively called by the new `run_trajectory` entry point in `trajectories.jl` (lines 1205, 1223). The file must be split: gap-estimation-specific code (the inner helpers `_run_trajectory_convergence` and `_run_trajectory_adaptive` that delegate to the main convergence functions) stays active, while the `estimate_spectral_gap` pipeline (gap_estimation.jl) and its fitting dependency (fitting.jl) can be staged. Alternatively, the convergence tracking functions themselves must remain in active `src/` with only gap_estimation.jl and fitting.jl moving to staging. log_sobolev.jl is already entirely commented out and can move cleanly.

**Primary recommendation:** Split convergence.jl -- keep all convergence tracking functions in active `src/` (they power `run_trajectory`'s convergence/adaptive modes), move only gap_estimation.jl and fitting.jl to `src/staging/`, and move log_sobolev.jl (entirely commented out) to staging. Then remove dead code (old entry points, `using Distributed`, dead structs) and reorganize exports.

## Architecture Patterns

### Current Source File Layout (28 files in src/)
```
src/
  QuantumFurnace.jl          -- module definition, using statements, exports, includes
  constants.jl               -- physical constants
  hamiltonian.jl             -- HamHam struct, Hamiltonian construction
  trotter_domain.jl          -- TrottTrott, Trotter decomposition
  structs.jl                 -- Config, JumpOp, Workspace, scratch types, result types, LSIFramework
  qi_tools.jl                -- quantum information utilities
  misc_tools.jl              -- validation, printing, precomputation
  time_domain.jl             -- TimeDomain support
  nufft.jl                   -- NUFFT wrappers (active, used by Krylov/trajectories/jump_workers)
  ofts.jl                    -- OFT functions (kept for debugging/pedagogy)
  errors.jl                  -- placeholder (2 lines, essentially empty)
  kraus.jl                   -- placeholder (2 lines, replaced by ThermalizeScratch)
  energy_domain.jl           -- EnergyDomain support
  bohr_domain.jl             -- BohrDomain support
  coherent.jl                -- coherent correction terms
  jump_workers.jl            -- jump operator processing
  trajectories.jl            -- trajectory simulation (run_trajectories, run_trajectory, step_along_trajectory!)
  furnace_utensils.jl        -- furnace helper utilities
  furnace.jl                 -- run_lindbladian, run_thermalization, run_lindblad, run_thermalize
  krylov_workspace.jl        -- Krylov workspace construction
  krylov_matvec.jl           -- Krylov matrix-vector products
  krylov_eigsolve.jl         -- Krylov eigensolver, run_krylov_spectrum
  log_sobolev.jl             -- entirely commented out (LSI framework)
  convergence.jl             -- convergence tracking, adaptive stopping (ACTIVELY USED)
  fitting.jl                 -- exponential decay fitting (FitResult)
  gap_estimation.jl          -- estimate_spectral_gap, eigenbasis_overlap_analysis
  diagnostics.jl             -- exact diagnostics analysis module
  results.jl                 -- result serialization (save_result/load_result)
```

### Target Export List Organization
```julia
# --- Lindbladian ---
export run_lindblad, construct_lindbladian
export LindbladResults
export apply_lindbladian!, apply_adjoint_lindbladian!
export krylov_spectral_gap, apply_delta_channel!

# --- Thermalize ---
export run_thermalize
export ThermalizeResults

# --- Krylov ---
export run_krylov_spectrum
export KrylovSpectrumResults

# --- Trajectory ---
export run_trajectory
export TrajectoryResults
export TrajectoryResult, ObservableTrajectoryResult  # old result types still used
export step_along_trajectory!, run_observable_trajectories
export run_trajectories  # old entry point still used in tests
export ConvergenceData, run_trajectories_convergence, run_trajectories_adaptive
export build_preset_trajectory_observables

# --- Diagnostics ---
export EigenDecompositionResult, FixedPointResult, DefectResult, OverlapResult, ...
export run_exact_diagnostics

# --- Common ---
export Config, AbstractSimulation, Lindbladian, Thermalize, KrylovSpectrum, Trajectory
export AbstractConstruction, KMS, GNS, DLL, with_coherent
export Workspace, LiouvillianScratch, ThermalizeScratch, KrylovScratch, TrajectoryScratch
export AbstractResults, save_result, load_result
export BohrDomain, EnergyDomain, TimeDomain, TrotterDomain
export HamHam, TrottTrott, JumpOp
export trace_distance_h, trace_distance_nh, ...  # qi_tools
export gibbs_state, gibbs_state_in_eigen, ...     # hamiltonian utils
export X, Y, Z, Had, pad_term, ...                # Pauli building blocks
export validate_config!
export oft!, time_oft!, trotter_oft!               # OFT (pedagogy/debugging)

# STAGING: estimate_spectral_gap, OverlapAnalysisResult, eigenbasis_overlap_analysis
# STAGING: fit_exponential_decay, FitResult
# STAGING: LSIFramework, compute_LSI_alpha2
```

### Anti-Patterns to Avoid
- **Moving functions to staging that are called by active code:** The `run_trajectory` entry point in trajectories.jl calls `_run_trajectory_convergence` and `_run_trajectory_adaptive` from convergence.jl. These delegate to `run_trajectories_convergence` and `run_trajectories_adaptive`. All these functions must stay.
- **Orphaned includes:** After moving files to staging, every `include("file.jl")` in QuantumFurnace.jl must be removed or the module will fail to load.
- **Breaking test imports:** Tests `using QuantumFurnace` will fail if staged exports are still referenced.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| File dependency analysis | Manual grep | Structured codebase analysis (already done here) | Easy to miss transitive dependencies |
| Test isolation for staging | Custom test runner | Simply remove `include()` from runtests.jl | Julia test infrastructure handles file-level inclusion |

## Common Pitfalls

### Pitfall 1: convergence.jl Has Active Dependencies
**What goes wrong:** Moving convergence.jl wholesale to staging breaks `run_trajectory` convergence and adaptive modes.
**Why it happens:** convergence.jl was originally created for trajectory gap estimation but its convergence-tracking functions were later integrated into the `run_trajectory` entry point (Phase 36-03).
**How to avoid:** convergence.jl must stay in active `src/`. Only gap_estimation.jl, fitting.jl, and log_sobolev.jl can move to staging.
**Warning signs:** Module load failure, test failures in test_trajectory_fixes.jl, test_threading.jl, test_workspace_independence.jl.

**Detailed dependency chain:**
- `trajectories.jl:1205` -> `_run_trajectory_adaptive()` (defined in convergence.jl:448)
- `trajectories.jl:1223` -> `_run_trajectory_convergence()` (defined in convergence.jl:408)
- Both call `run_trajectories_convergence` / `run_trajectories_adaptive` (convergence.jl)
- Both use `build_preset_trajectory_observables` (convergence.jl)
- Both use `_build_framework_and_seed` and `_run_batch_no_obs!` (trajectories.jl -- these are fine)
- Both use `_compute_gibbs_observable_values` and `_windowed_relative_change` (convergence.jl)

### Pitfall 2: ConvergenceData Struct in structs.jl
**What goes wrong:** ConvergenceData is defined in structs.jl, not in convergence.jl. Moving convergence.jl to staging still leaves ConvergenceData in active code.
**Why it happens:** ConvergenceData must be available before trajectories.jl (which uses it in TrajectoryResult).
**How to avoid:** ConvergenceData stays in structs.jl. It is used by TrajectoryResults (the new Phase 36 result type) and TrajectoryResult (the old result type), both active.

### Pitfall 3: LSIFramework in structs.jl
**What goes wrong:** LSIFramework struct is defined in structs.jl (lines 288-318) but is only used by the commented-out code in log_sobolev.jl.
**Why it happens:** Struct definitions were centralized in structs.jl.
**How to avoid:** When staging log_sobolev.jl, also remove LSIFramework from structs.jl. Or keep it (it is harmless dead code in structs.jl). Recommendation: remove it as part of dead code cleanup since it is unused.

### Pitfall 4: DMSimulationResult and LindbladianResult Are Dead
**What goes wrong:** These old result types in structs.jl are only used by the old entry points (`run_lindbladian`, `run_thermalization`).
**Why it happens:** Phase 36 introduced new result types (LindbladResults, ThermalizeResults) with metadata support. The old types are still exported.
**How to avoid:** When deleting `run_lindbladian` and `run_thermalization`, also delete `DMSimulationResult` and `LindbladianResult` from structs.jl, and remove their export lines.
**Verification:** `DMSimulationResult` -- only used in furnace.jl:156 (run_thermalization return). `LindbladianResult` -- only used in furnace.jl:29 (run_lindbladian return) and log_sobolev.jl:18 (commented out). Neither is used in tests.

### Pitfall 5: Old Entry Points Still Referenced in Active Code Comments
**What goes wrong:** Docstrings and comments reference `run_lindbladian` and `run_thermalization` after they are deleted.
**Why it happens:** Natural during code evolution.
**Where:** `src/furnace.jl:199` comment "same logic as existing run_lindbladian", `src/convergence.jl:96` comment "Follows the same transform as furnace.jl run_thermalization", `src/krylov_eigsolve.jl:438` comment "matching what run_thermalization", `src/structs.jl:360` comment "DM Kraus evolution (run_thermalization)".
**How to avoid:** Update or remove stale comments referencing deleted functions (Claude's discretion).

### Pitfall 6: OFTCaches Struct Stays
**What goes wrong:** OFTCaches might be confused for dead code since OFTs are "kept for debugging/pedagogy."
**Why it happens:** OFTCaches is used by ofts.jl which is actively included and exports `oft!`, `time_oft!`, `trotter_oft!`.
**How to avoid:** Keep OFTCaches. The user explicitly said to keep test utilities like `time_oft!`, `trotter_oft!`.

### Pitfall 7: Project.toml Distributed Entry
**What goes wrong:** Removing `using Distributed` from the module without removing it from Project.toml leaves a stale dependency.
**Why it happens:** Project.toml [deps] and QuantumFurnace.jl `using` statements must be kept in sync.
**How to avoid:** Remove `Distributed = "..."` from [deps] AND its [compat] entry. No `@distributed` usage exists anywhere in src/. The `using Distributed` on line 12 of QuantumFurnace.jl is the only reference.

### Pitfall 8: TrajectoryResult and ObservableTrajectoryResult Are Still Active
**What goes wrong:** These old result types might be confused with dead code since Phase 36 introduced TrajectoryResults.
**Why it happens:** `TrajectoryResult` is still the internal return type of `run_trajectories`, `run_trajectories_convergence`, `run_trajectories_adaptive`, and `run_observable_trajectories`. `ObservableTrajectoryResult` is returned by `run_observable_trajectories`.
**How to avoid:** Keep both. They are actively exported and used by internal code and tests. `run_trajectory` (new entry point) wraps these into `TrajectoryResults` at its boundary.

### Pitfall 9: errors.jl and kraus.jl Are Placeholder Files
**What goes wrong:** These files contain 1-2 lines (comments only) and are included by the module.
**Why it happens:** They were kept as placeholders after earlier refactoring.
**How to avoid:** Delete both files and remove their `include()` lines from QuantumFurnace.jl. They contribute zero code.

### Pitfall 10: run_trajectories Is Not Targeted for Deletion
**What goes wrong:** Overzealous cleanup deletes `run_trajectories` thinking it is superseded by `run_trajectory`.
**Why it happens:** The names are very similar and `run_trajectory` is the new Phase 36 entry point.
**How to avoid:** `run_trajectories` is still exported, actively used in tests (test_gns_trajectory.jl, test_observable_trajectories.jl, test_threading.jl, test_workspace_independence.jl, trajectory_validation/), and used internally by `run_trajectory`. It was NOT listed for deletion. Keep it.

## Code Examples

### Staging File Move Pattern
```bash
# Create staging directories
mkdir -p src/staging
mkdir -p test/staging

# Move files (git mv preserves history)
git mv src/gap_estimation.jl src/staging/gap_estimation.jl
git mv src/fitting.jl src/staging/fitting.jl
git mv src/log_sobolev.jl src/staging/log_sobolev.jl

# Move related tests
git mv test/test_gap_estimation.jl test/staging/test_gap_estimation.jl
git mv test/test_fitting.jl test/staging/test_fitting.jl
```

### Remove Include Lines from Module
```julia
# REMOVE these lines from QuantumFurnace.jl:
# include("log_sobolev.jl")    -- staged
# include("fitting.jl")        -- staged
# include("gap_estimation.jl") -- staged
# include("errors.jl")         -- empty placeholder, delete file
# include("kraus.jl")          -- empty placeholder, delete file
```

### Remove Dead Exports
```julia
# REMOVE from QuantumFurnace.jl exports:
# export LindbladianResult, DMSimulationResult          -- old result types
# export run_lindbladian, run_thermalization             -- old entry points
# export LSIFramework, compute_LSI_alpha2               -- already commented out
```

### Delete Dead Functions from furnace.jl
```julia
# DELETE function run_lindbladian(...) on lines 1-36 of furnace.jl
# DELETE function run_thermalization(...) on lines 91-161 of furnace.jl
# KEEP construct_lindbladian (used by run_lindblad)
# KEEP run_lindblad (new entry point)
# KEEP run_thermalize (new entry point)
```

### Remove Distributed from Module
```julia
# In QuantumFurnace.jl line 12, REMOVE:
# using Distributed

# In Project.toml [deps], REMOVE:
# Distributed = "8ba89e20-285c-5b6f-9357-94700520ee1b"

# In Project.toml [compat], REMOVE:
# Distributed = "1.11, 1.12"
```

## Detailed File-by-File Analysis

### Files to Move to src/staging/

| File | Size | Status | Reason | Dependencies |
|------|------|--------|--------|-------------|
| gap_estimation.jl | 12.9 KB | Can stage | Contains `estimate_spectral_gap`, `eigenbasis_overlap_analysis`, `OverlapAnalysisResult`. Only used by test_gap_estimation.jl. No active code depends on it. | Calls `build_preset_trajectory_observables` (convergence.jl, stays active), `run_observable_trajectories` (trajectories.jl, stays active), `fit_exponential_decay` (fitting.jl, also staged) |
| fitting.jl | 8.2 KB | Can stage | Contains `fit_exponential_decay`, `FitResult`, internal helpers. Only called by gap_estimation.jl (also staged) and test_fitting.jl. | Uses LsqFit (external package) |
| log_sobolev.jl | 8.4 KB | Can stage | Entirely commented out. Contains `compute_LSI_alpha2` function body, all commented. Zero active code. | None (commented out) |
| **convergence.jl** | **18.9 KB** | **CANNOT stage** | Contains functions actively called by `run_trajectory` in trajectories.jl. See Pitfall 1 for full dependency chain. | Called by trajectories.jl:1205, 1223 |

### Files/Structs to Delete (Dead Code)

| Item | Location | Reason |
|------|----------|--------|
| `run_lindbladian()` | furnace.jl:1-36 | Old entry point, replaced by `run_lindblad()` |
| `run_thermalization()` | furnace.jl:91-161 | Old entry point, replaced by `run_thermalize()` |
| `DMSimulationResult` struct | structs.jl:126-142 | Only returned by `run_thermalization()` (being deleted) |
| `LindbladianResult` struct | structs.jl:145-163 | Only returned by `run_lindbladian()` (being deleted), referenced in commented-out log_sobolev.jl |
| `LSIFramework` struct | structs.jl:288-318 | Only used by commented-out code in log_sobolev.jl |
| `using Distributed` | QuantumFurnace.jl:12 | No `@distributed` usage anywhere in src/ |
| Distributed in Project.toml | Project.toml:11, 37 | Follows from removing `using Distributed` |
| errors.jl | src/errors.jl | Placeholder file (1 comment line, no code) |
| kraus.jl | src/kraus.jl | Placeholder file (1 comment line, no code) |

### Files That Stay (Confirmed Active)

| File | Why Active |
|------|-----------|
| convergence.jl | Powers `run_trajectory` convergence/adaptive modes. Exports `ConvergenceData`, `run_trajectories_convergence`, `run_trajectories_adaptive`, `build_preset_trajectory_observables`. |
| nufft.jl | Used by krylov_matvec.jl, krylov_workspace.jl, jump_workers.jl, trajectories.jl, furnace_utensils.jl |
| ofts.jl | Exports `oft!`, `time_oft!`, `trotter_oft!` (user explicitly wants to keep) |

### Tests to Move to test/staging/

| Test File | For | Status |
|-----------|-----|--------|
| test_gap_estimation.jl | gap_estimation.jl | Move (gap_estimation.jl staged). Note: line 201 calls `run_lindbladian` (being deleted) -- this test line needs commenting out or removal. |
| test_fitting.jl | fitting.jl | Move (fitting.jl staged) |
| test_convergence.jl | convergence.jl | **KEEP ACTIVE** -- convergence.jl stays active |

### Backward-Compat Config Type Aliases

Research confirmed: **No config type aliases exist in the current codebase.** `LiouvConfig`, `ThermalizeConfig`, `LiouvConfigGNS`, `ThermalizeConfigGNS` appear only in old planning docs, README, and tutorial files -- not in any `src/` code. The unified `Config{S,D,C,T}` is already the only config type. No code changes needed for this item.

### Old Types Referenced Only in Non-Source Locations

- `LiouvConfig`, `ThermalizeConfig` appear in: README.md, docs/, old planning files
- These should be updated in README.md and docs/ as part of cleanup (Claude's discretion on whether to include in this phase)

## Simulation Scripts Status (ORG-07)

The phase requirements include ORG-07: Create/update 4 simulation scripts in `simulations/` matching the 4 `run_*` entry points.

Current state:
| Script | Entry Point | Status |
|--------|-------------|--------|
| main_liouv.jl | `run_lindblad` | Already uses new API (line 121) |
| main_thermalize.jl | `run_thermalize` | Already uses new API (line 121) |
| main_trajectory.jl | `run_trajectory` | Already uses new API (line 89) |
| main_krylov.jl | `run_krylov_spectrum` | **MISSING** -- only main_krylov_benchmark.jl exists (26KB, uses internal `krylov_spectral_gap` not the new `run_krylov_spectrum`) |

A simple `main_krylov.jl` script demonstrating `run_krylov_spectrum` needs to be created.

## Open Questions

1. **test_gap_estimation.jl line 201 calls `run_lindbladian`**
   - What we know: This test uses the old `run_lindbladian` to build a Liouvillian for overlap analysis tests
   - What's unclear: Whether to delete this test entirely (it moves to staging anyway) or update it to use `run_lindblad`
   - Recommendation: Since the entire test file moves to test/staging/ (dormant), it does not need updating. When the code is reactivated from staging, it can be updated then.

2. **Old `run_trajectories` and `run_observable_trajectories` exports**
   - What we know: These are still actively used in tests and internal code
   - What's unclear: Whether they should eventually be deprecated (not in this phase)
   - Recommendation: Keep them exported. They are NOT listed for deletion. Mark them clearly in the export list under `# Trajectory` section.

3. **ConvergenceData backward-compatible 6-arg constructor** (structs.jl:200-211)
   - What we know: It exists for "Phase 16 callers" who pass 6 args instead of 10
   - What's unclear: Whether this is still needed
   - Recommendation: Keep it. It is a convenience constructor, not a type alias. It does no harm and the convergence tracking code actively creates ConvergenceData with the 6-arg form (convergence.jl:231, 240).

## Sources

### Primary (HIGH confidence)
- Direct codebase analysis of all 28 source files in `src/`
- Direct analysis of `test/runtests.jl` and all 25 test files
- Direct analysis of 4 simulation scripts in `simulations/`
- grep/read of all cross-references between files

### Verification Protocol
- Every staging candidate file was checked for callers across the entire `src/` directory
- Every dead code candidate was checked for references in `src/`, `test/`, and `simulations/`
- The `run_trajectory` -> convergence.jl dependency chain was traced line by line

## Metadata

**Confidence breakdown:**
- Dead code identification: HIGH -- all references traced through src/, test/, simulations/
- Staging file analysis: HIGH -- convergence.jl dependency chain fully verified
- Export list design: HIGH -- based on exhaustive review of current exports and what each function does
- Simulation scripts: HIGH -- all 4 scripts read and analyzed

**Research date:** 2026-02-27
**Valid until:** 2026-03-27 (stable internal refactoring, no external dependencies changing)
