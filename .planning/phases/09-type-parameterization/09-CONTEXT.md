# Phase 9: Type Parameterization - Context

**Gathered:** 2026-02-15
**Status:** Ready for planning

<domain>
## Phase Boundary

Parameterize core structs on element type (`T<:AbstractFloat`) so the library supports Float64 (default) and Float32 precision without changing existing calling code. This phase adds the type parameter infrastructure; it does not add Float32-specific optimizations or benchmarks.

</domain>

<decisions>
## Implementation Decisions

### Default type behavior
- Float64 is the default everywhere — users never see `T` unless they want Float32
- Single type parameter `T` controls everything: `T` for real fields, `Complex{T}` for complex fields (e.g., `T=Float32` → `Float32` reals, `ComplexF32` density matrices)
- Users opt into Float32 via `precision=Float32` keyword argument on constructors

### Mixed-precision policy
- **Downward mismatch errors:** Passing Float64 data to a `precision=Float32` constructor throws an error ("Expected Float32 data, got Float64") — strict, no silent precision loss
- **Upward promotion allowed:** Passing Float32 data to a Float64 (default) constructor silently promotes — no precision loss, so it's safe
- **Cross-struct must match:** All structs in a simulation must share the same `T` — error if a Float32 workspace meets a Float64 config
- **No conversion utilities:** Users reconstruct with the right precision rather than converting existing structs

### Parameterization scope
- HamHam, LindbladianWorkspace, and Config structs (KMSConfig, GNSConfig, etc.) all get `T` parameter
- Return types from simulation functions match input `T` — Float32 in → Float32 out
- TrajectoryFramework stays unparameterized on `T` — infers from its components
- TrottTrott and scratch/cache arrays: Claude's discretion based on what fields actually need `T`

### API feel
- Keyword name: `precision=Float32` (not `eltype` or `T`)
- `precision=` kwarg appears only on top-level entry points (HamHam, Config) — Workspace and internals infer from HamHam
- Display always shows `T`: `HamHam{Float64}(...)` even for the default
- Fully backward compatible: all existing code that doesn't specify `precision=` continues to work unchanged

### Claude's Discretion
- Whether TrottTrott needs the T parameter (depends on its actual field types)
- How deep scratch/cache array parameterization goes
- Internal type propagation mechanics
- Exact constructor implementation for the precision kwarg

</decisions>

<specifics>
## Specific Ideas

- "The moment we opt into F32, density matrices etc. should all be ComplexF32" — single T controls both real and complex field types
- Config structs get the same precision kwarg — Float32 tolerances when running in Float32 mode

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 09-type-parameterization*
*Context gathered: 2026-02-15*
