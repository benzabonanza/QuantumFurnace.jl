# Research Summary: QuantumFurnace.jl Codebase Restructure

**Domain:** Refactoring a Julia quantum Gibbs sampling simulation package (8,312 LOC src + 5,071 LOC test)
**Researched:** 2026-02-25
**Overall confidence:** HIGH

## Executive Summary

QuantumFurnace.jl has grown organically through 37+ development phases, accumulating significant structural debt: 4 config structs with 52 duplicated field lines, 5 workspace structs with overlapping scratch buffers, prefactor formulas copied verbatim 8+ times across 6 files, and 28 source files in a flat directory with no logical grouping. The restructure milestone addresses all four dimensions using Julia's type system and multiple dispatch as the primary tools -- no new dependencies required.

The central design decision is encoding the three orthogonal axes of variation (Simulation mode, Domain, Detailed Balance flavor) as type parameters on a single `SimConfig{S, D, DB, T}` struct, replacing the current 4 separate config types and their Union-type dispatch. This enables zero-cost compile-time dispatch on any axis independently: `pick_transition` dispatches on DB only, `compute_prefactor` dispatches on Domain only, `run_simulation` dispatches on SimMode only. The compiler specializes each combination, producing identical machine code to the current hand-duplicated methods but from a single source definition.

Workspace consolidation merges `KrylovWorkspace`, `KrausScratch`, and `LindbladianWorkspace` into a single `SimWorkspace` struct with `Union{Nothing, Matrix{T}}` fields for simulation-mode-specific buffers. `TrajectoryWorkspace` remains separate (state-vector vs density-matrix scratch). The `Union{Nothing, T}` pattern is safe here because these optional fields are accessed at setup boundaries, never in hot loops, and Julia's isbits union optimization stores them inline.

Code deduplication extracts domain-specific computations into small `@inline` dispatched helper functions (`compute_rate_prefactor`, `compute_jump_oft!`) and a unified frequency iteration pattern (`foreach_frequency`). The forward/adjoint sandwich variants collapse from 4 functions to 1 function with a `Val{:forward/:adjoint}` parameter. The 6 `apply_(adjoint_)lindbladian!` methods reduce to a single `_apply_lindbladian_impl!` with domain and direction dispatch. Estimated savings: ~600 LOC from deduplication alone.

## Key Findings

**Stack:** No new dependencies. Pure refactoring using Julia parametric types, multiple dispatch, and `@inline` annotations.
**Architecture:** `SimConfig{S, D, DB, T}` as the central type, with SimWorkspace and domain-dispatched helpers for zero-overhead deduplication.
**Critical pitfall:** Breaking type stability during incremental refactoring -- a single abstract-typed field access in a hot loop regresses from 0 allocations to thousands.

## Implications for Roadmap

Based on research, suggested phase structure:

1. **Phase 1: Type Foundation** - Define the new type hierarchy (SimConfig, DB flavors, SimMode types) with backward-compatible aliases
   - Addresses: Config struct duplication, future DLL extensibility
   - Avoids: Breaking all existing tests at once (aliases preserve old API)
   - Research confidence: HIGH (standard Julia parametric type patterns)

2. **Phase 2: Transition/Alpha Unification** - Consolidate `pick_transition` (4 methods to 2), `_pick_alpha` (4 methods to 2), and prefactor formulas (8 locations to 2 dispatched helpers)
   - Addresses: Most impactful code duplication (formulas that must stay in sync)
   - Avoids: Touching hot-path loop structures yet (lower risk)
   - Research confidence: HIGH

3. **Phase 3: Workspace Consolidation** - Merge KrylovWorkspace + KrausScratch + LindbladianWorkspace into SimWorkspace
   - Addresses: 5 workspace structs to 2, eliminates conversion helpers
   - Avoids: Premature optimization of loop structures (workspace is setup-time)
   - Research confidence: HIGH

4. **Phase 4: Hot-Path Deduplication** - Unify sandwich functions, frequency iteration, apply_lindbladian! methods
   - Addresses: The 600-line deduplication in performance-critical code
   - Avoids: Must come AFTER workspace consolidation (depends on unified ws type)
   - Research confidence: HIGH but needs careful allocation testing

5. **Phase 5: File Reorganization** - Move files into subdirectory structure, update includes
   - Addresses: Navigation and discoverability for 22+ files
   - Avoids: Functional changes (pure file moves)
   - Research confidence: HIGH (mechanical transformation)

6. **Phase 6: Test Cleanup** - Consolidate test helpers, eliminate duplicate setups
   - Addresses: Test maintainability
   - Avoids: Changing test logic (only restructuring helpers)
   - Research confidence: MEDIUM (test specifics need per-file analysis)

**Phase ordering rationale:**
- Phase 1 before all others: The new type hierarchy is the foundation everything else builds on
- Phase 2 before 3: Transition/alpha unification is purely function-level, no struct changes, validates the dispatch pattern
- Phase 3 before 4: Hot-path deduplication requires the unified workspace type to exist
- Phase 5 after 4: File reorganization is safest after the code structure stabilizes
- Phase 6 last: Tests should be cleaned up after the code they test is finalized

**Research flags for phases:**
- Phase 1: Standard patterns, unlikely to need further research
- Phase 4: Likely needs deeper allocation profiling (`@allocated`, `@code_warntype`) to verify zero-overhead guarantees for closure-based frequency iteration
- Phase 6: Likely needs per-file analysis of test helpers to identify consolidation targets

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | No new dependencies; Julia type system patterns well-documented |
| Config type hierarchy | HIGH | Standard parametric type pattern, used by DiffEqBase.jl, verified against Julia docs |
| Workspace consolidation | HIGH | Union{Nothing, T} for setup-time fields is explicitly documented as efficient in Julia |
| Hot-path deduplication | HIGH | @inline + dispatch on singleton types is the standard zero-overhead pattern |
| File organization | HIGH | Flat src/ with subdirectories is mechanical; standard Julia convention |
| Test cleanup | MEDIUM | Requires per-file analysis not yet performed |

## Gaps to Address

- **DLL variant**: The supplementary notes mention a future DLL detailed balance flavor. The `DLLDB <: AbstractDBFlavor` type is designed to accommodate this, but the specific fields/functions for DLL have not been analyzed yet.
- **Bohr domain hot-path**: The Bohr domain uses a fundamentally different loop structure (iterating over Bohr dictionary keys with 2-operator sandwiches). It may resist full unification with Energy/Time/Trotter and may need a separate `_apply_lindbladian_bohr_impl!`.
- **Trajectory framework**: `TrajectoryFramework` stores many concrete-typed hot-path fields extracted from config/precomputed_data. After the config restructure, these extractions may need updating. The research covers the pattern but not the specific field mapping.
- **BSON serialization**: Changing config struct layout may break BSON.load of existing saved experiments. A migration function or versioned config format may be needed.
