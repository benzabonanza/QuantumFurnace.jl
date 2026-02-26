# Phase 35: Workspace and Channel Consolidation - Context

**Gathered:** 2026-02-26
**Status:** Ready for planning

<domain>
## Phase Boundary

Merge KrylovWorkspace, KrausScratch, and LindbladianWorkspace into a unified parameterized `Workspace{S,D,C,T}` struct. TrajectoryWorkspace also becomes `Workspace{Trajectory,D,C,T}` following the same pattern. CPTP channel computation uses shared helper functions with correct per-jump vs summed semantics.

</domain>

<decisions>
## Implementation Decisions

### Unified Workspace Struct Design
- Parameterized `Workspace{S,D,C,T}` matching Config's type parameters
- S expanded to 4 simulation singletons: `Lindbladian`, `Thermalize`, `Krylov`, `Trajectory`
- Each Workspace{S,...} has different fields appropriate to its simulation path
- Nested structure: precomputed/immutable physics data as flat fields, mutable scratch buffers bundled into a single nested `Scratch` struct field
- `Workspace{Krylov,...}` handles both Krylov-Lindbladian (G_left/G_right) and Krylov-Thermalize (CPTP channel) modes as a single type, using `Union{Nothing,...}` for mode-specific fields
- Constructor: `Workspace(config::Config{S,D,C,T}, hamiltonian, jumps)` — type constructor, dispatch on Config
- **Critical constraint:** Only ONE constructor per workspace variant (inner OR outer, never both) — multiple constructors caused bugs previously
- Domain-specific precomputed data (NUFFT for TimeDomain/TrotterDomain) absorbed into the workspace parameterization, replacing the generic `precomputed_data::NamedTuple`

### CPTP Channel Computation API
- `_build_cptp_channel(R, delta)` stays as a single pure function — callers pass R_total (summed) for Krylov or R^a (per-jump) for DM/Trajectory
- Function doesn't know or care whether input R is summed or per-jump
- Coherent B computation stays separate from `_build_cptp_channel` — clean physics separation
- Eliminate single-jump variants of `B_time`, `B_trotter`, etc. — keep only the `jumps` (vector) variant; callers wrap single jump as `[jump]`
- Same deduplication for `_precompute_coherent_B` and `_precompute_coherent_unitary` — vector-of-jumps only
- K0, U_residual, alpha stored as direct flat fields on the workspace (not wrapped in NamedTuple or sub-struct)
- G_left, G_right, G_left_adj, G_right_adj as direct workspace fields with `Union{Nothing,...}` for non-Krylov-Lindbladian workspaces

### TrajectoryWorkspace Integration
- `Workspace{Trajectory,D,C,T}` follows same pattern as other simulation types
- PerOperatorKraus data absorbed: per-jump R^a, K0^a, U_residual^a, U_B^a stored as vectors of matrices (e.g., `K0s::Vector{Matrix{T}}`) — flattened, not nested in PerOperatorKraus
- Only mutable scratch buffers (jump_oft, psi_tmp, Rpsi, rho_acc) go into the nested Scratch struct

### Field Naming
- Aim for physics-descriptive names for scratch buffers (not generic tmp1/tmp2)
- No `channel_` prefix on CPTP fields: `K0`, `U_residual` directly on workspace
- Claude to audit dead/redundant fields across all three workspaces during research (KrausScratch.K0 already known dead from Phase 34-02)

### Claude's Discretion
- Exact Scratch struct field names (descriptive but Claude determines the best names per simulation path)
- Whether domain-specific precomputed data gets its own type parameter or is handled via dispatch in constructors
- How to handle the Identity matrix (currently a field in LindbladianWorkspace — may become computed or shared)
- Optimal number of scratch matrices per simulation type

</decisions>

<specifics>
## Specific Ideas

- "Only one constructor, either inner or outer, but only one. Multiple constructors led to bugs previously" — hard constraint on workspace construction
- `[jump]` wrapping pattern for single-jump → vector-of-jumps unification at call sites
- Workspace{Krylov,...} is a single type handling both Lindbladian and Thermalize Krylov modes

</specifics>

<deferred>
## Deferred Ideas

- Adding `Krylov` and `Trajectory` as new Config simulation singleton types may overlap with Phase 36 (API and Results). The workspace phase defines the singletons; Phase 36 wires them to `run_*` entry points.

</deferred>

---

*Phase: 35-workspace-and-channel-consolidation*
*Context gathered: 2026-02-26*
