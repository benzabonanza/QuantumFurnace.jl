# Phase 12: Workspace Refactor - Context

**Gathered:** 2026-02-15
**Status:** Ready for planning

<domain>
## Phase Boundary

Separate mutable workspace from shared framework to enable thread-safe trajectory stepping. `step_along_trajectory!` takes explicit workspace and RNG arguments. Two independent workspaces can step from the same framework without interference. This is internal restructuring — the high-level thermalization API signature does not change.

</domain>

<decisions>
## Implementation Decisions

### Workspace scope
- Density matrix accumulator lives inside the workspace (not separate)
- `TrajectoryWorkspace{T}` matches framework's type parameterization for type safety
- Claude's discretion on which scratch arrays go into workspace vs stay in framework (goal: framework becomes read-only during stepping)
- Claude's discretion on constructor design (from-framework vs independent)

### User-facing API
- Workspace is internal — not exported or exposed to users
- High-level API (thermalize functions) keeps the same external signature; workspace/RNG managed inside
- Trajectory run returns a `TrajectoryResult` struct containing averaged density matrix + step count (minimal — no forward-looking fields for convergence etc.)
- Result struct extended by later phases as needed

### Backward compatibility
- Clean break: old function signatures removed, all call sites updated to new signatures
- No deprecation wrappers — direct migration
- Test logic and assertions stay the same; call sites updated to new signatures
- `TrajectoryFramework` struct modified in place (mutable fields move to workspace) — no duplicate types

### RNG contract
- Low-level `step_along_trajectory!` takes `AbstractRNG` as explicit argument
- High-level API takes a seed integer; RNG created internally from seed (Xoshiro default)
- If no seed provided, auto-generate random seed from system entropy
- Seed stored in `TrajectoryResult` — every run is reproducible after the fact

</decisions>

<specifics>
## Specific Ideas

- "The averaged DM is the star of this simulation" — workspace is internal but the DM result must flow back to the user via the result struct
- Every run should be reproducible: random seeds are captured and stored even when not explicitly provided

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 12-workspace-refactor*
*Context gathered: 2026-02-15*
