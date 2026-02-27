# Phase 37: File Organization and Dead Code - Context

**Gathered:** 2026-02-27
**Status:** Ready for planning

<domain>
## Phase Boundary

Reorganize source files, remove dead code and backward-compat shims, move trajectory gap estimation code to staging, and clean up the module export list. File renaming/PRE-MID-POST reorganization of src/ is **excluded** (user will handle manually).

</domain>

<decisions>
## Implementation Decisions

### Staging area design
- Move trajectory gap estimation files to `src/staging/` subdirectory
- Files to stage: fitting.jl, gap_estimation.jl, log_sobolev.jl (convergence.jl stays active — research confirmed it powers run_trajectory's convergence/adaptive modes via _run_trajectory_convergence and _run_trajectory_adaptive)
- **Critical:** Verified: convergence.jl contains functions called by trajectories.jl:1205,1223 — must stay in active src/
- Staging code is **excluded from the module** — not included in QuantumFurnace.jl
- Related tests move to `test/staging/` and are **not run** by the regular test suite (dormant)

### Dead code removal
- Remove `@distributed` code and `using Distributed` import from furnace.jl; SharedArrays stays
- Delete old entry points outright: `run_lindbladian()`, `run_thermalization()` — no deprecation warnings
- Delete all backward-compat config type aliases (LiouvConfig, ThermalizeConfig, LiouvConfigGNS, ThermalizeConfigGNS) if any remain
- Delete any dead or backward-compat structs — target struct landscape is:
  - `Config{S,D,C,T}` — unified config
  - `Workspace{S,D,C,T}` — unified workspace (4 type params, no scratch type param)
  - 4 scratch types (one per simulation kind)
  - 4 Result structs (LindbladResults, ThermalizeResults, KrylovSpectrumResults, TrajectoryResults)
- Conservative approach overall: keep explicitly useful test utilities like `time_oft!`, `trotter_oft!` even if not called in production code
- Clean break: no deprecation shims, no "one more release" grace period

### Export list organization
- Group by simulation type: `# Lindbladian`, `# Thermalize`, `# Krylov`, `# Trajectory` comment blocks
- Shared types/utilities (Config, domain types, qi_tools, hamiltonian utilities) go in `# Common` section **at the bottom**
- Staging code exports: keep as **commented-out block** with `# STAGING:` prefix for reference when code is reactivated
- Diagnostics types exported from main QuantumFurnace module directly (not a separate submodule)
- Old/deleted function names removed from export list entirely

### Claude's Discretion
- Exact ordering of exports within each simulation type section
- Which comment/docstring references to old types/functions to clean up vs leave
- Whether to clean up stale comments referencing deleted code in non-staging files

</decisions>

<specifics>
## Specific Ideas

- "We basically want Workspace parametrized on config stuff, related 4 scratch structs, 4 main result structs for the 4 simulations"
- nufft.jl is active core code (used by Krylov, jump_workers, trajectories, furnace_utensils) — do NOT move to staging
- construct_lindbladian is still used by the new run_lindblad — only the old wrapper entry points are dead

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 37-file-organization-and-dead-code*
*Context gathered: 2026-02-27*
