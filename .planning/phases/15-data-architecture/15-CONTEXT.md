# Phase 15: Data Architecture - Context

**Gathered:** 2026-02-16
**Status:** Ready for planning

<domain>
## Phase Boundary

Experiment result serialization and round-trip via BSON. An experiment result (config, convergence history, observables, density matrices) can be saved to disk and loaded back with all fields matching to machine precision. This phase builds the data structs and save/load functions — convergence tracking logic is Phase 16, adaptive sampling is Phase 17.

</domain>

<decisions>
## Implementation Decisions

### Result struct shape
- Single flat struct `ExperimentResult{C,T}` parameterized on config type `C<:AbstractConfig` and element type `T`
- Reuse existing result types (TrajectoryResult, LindbladianResult, etc.) where they already carry the right data — embed them rather than duplicating fields
- Full config object embedded in the result — everything needed to re-run the experiment is stored
- Config type parameter `C` preserves KMS vs GNS distinction at the type level (dispatch-friendly)

### File organization
- Descriptive file names: e.g., `kms_n4_beta10_trotter_20260216.bson`
- Default output directory `results/` in project root, with optional path override
- Subdirectories by construction type: `results/kms/`, `results/approx_gns/` (and later `results/ding/`)
- Overwriting existing files is allowed silently — re-running replaces the old result
- Companion `.txt` file saved alongside each `.bson` with key parameters for browsing without Julia

### Metadata scope
- Core: seed, thread count, Julia version, timestamp
- Provenance: git commit hash (auto-captured at save time via LibGit2 or shell)
- Hamiltonian: store parameters (n, J, h) plus the coefficient matrices and term vectors (small, not the full eigendecomposition)
- Performance: wall-clock time and total trajectory count
- No package version field — git hash is sufficient

### Forward compatibility
- No schema version number — keep it simple
- Missing fields on load should fill with `nothing`/defaults (best-effort partial load)
- Experiments are cheap enough to re-run if schema diverges significantly

### Claude's Discretion
- Whether to use dedicated save_experiment/load_experiment wrapper functions vs raw BSON.bson/BSON.load
- Whether ExperimentResult wraps TrajectoryResult as a field or flattens its contents
- Exact companion .txt format and content
- Auto-capture mechanism for git hash (LibGit2 vs shell command)

</decisions>

<specifics>
## Specific Ideas

- Directory structure mirrors the paper's comparison structure: `results/kms/` vs `results/approx_gns/` maps directly to KMS-vs-GNS comparison in Phase 18
- A `results/ding/` folder will be added in a future milestone for the Ding et al. construction
- Companion text file enables quick browsing of sweep results without starting Julia

</specifics>

<deferred>
## Deferred Ideas

- `results/ding/` directory for Ding et al. (2024) construction — future milestone
- Index/manifest file listing all experiments in a sweep — could be useful for Phase 18 but not needed for the data layer itself

</deferred>

---

*Phase: 15-data-architecture*
*Context gathered: 2026-02-16*
