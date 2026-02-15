# Phase 8: Struct Simplification - Context

**Gathered:** 2026-02-15
**Status:** Ready for planning

<domain>
## Phase Boundary

Core data structures are minimal and correct: config struct field duplication reduced via shared constructors, HamHam fully initialized in a single constructor call (no Nothing fields), TrottTrott immutable with correct field types, and unnecessary type parametrization simplified on TrajectoryFramework. All 224 tests must pass.

</domain>

<decisions>
## Implementation Decisions

### Config field deduplication
- Keep 4 separate flat structs (LiouvConfig, LiouvConfigGNS, ThermalizeConfig, ThermalizeConfigGNS) for KMS/GNS/future Ding extensibility
- Extract a shared `_build_common_fields()` helper function that each variant's constructor calls
- GNS variants keep the `with_coherent` field but the constructor enforces `false`
- Keep `@kwdef` on config structs for ergonomic keyword construction
- Replace -1 sentinel defaults with `Union{Float64, Nothing} = nothing` for t0, w0, eta, num_energy_bits, num_trotter_steps_per_t0
- Full cleanup: update ALL downstream checks from `== -1` to `isnothing()` in the same phase

### HamHam initialization redesign
- Eliminate the two-step initialization pattern entirely
- `find_ideal_heisenberg` returns a NamedTuple of raw data (matrix, eigvals, eigvecs, terms, etc.) instead of a partially-initialized HamHam
- HamHam constructor takes raw data + beta and computes bohr_freqs, bohr_dict, and gibbs state directly
- No more `Union{..., Nothing}` fields on HamHam
- `finalize_hamham` is eliminated entirely (not kept as wrapper)
- Keep base_terms, base_coeffs, disordering_term, disordering_coeffs fields for inspection/debugging

### TrottTrott immutability
- Make TrottTrott an immutable struct (researcher should verify no mutation sites exist first)
- Change `num_trotter_steps_per_t0` from `Float64` to `Int` (platform-native, not explicit Int64)
- Field naming: rename `bohr_freqs` to `quasi_bohr_freqs` IF the code can stay simple with the rename (polymorphic `ham_or_trott.bohr_freqs` access currently exists). If renaming breaks polymorphism elegance, keep as `bohr_freqs`

### Type parametrization cleanup
- TrajectoryFramework{T,C,H,PD,D}: reduce to only essential type parameters (1-2). Researcher should investigate which params enable dispatch or performance; remove the rest
- Config domain parametrization (Config{D}): KEEP the domain type param, but leverage it properly throughout the codebase
- Refactor dispatch signatures from `f(config, domain::TimeDomain, ...)` to `f(config::Config{TimeDomain}, ...)` so domain dispatches through the config type param (eliminating redundant domain arguments)

### API breakage tolerance
- Clean break: no backward compatibility shims or deprecation warnings
- Change all constructor signatures to be correct, update all call sites

### Claude's Discretion
- Whether to remove the `domain::D` field from configs (since the type param carries the same info) — depends on whether `config.domain` is accessed for non-dispatch purposes
- Exact shared helper function signature for config construction
- Which TrajectoryFramework type params are essential vs removable (after investigation)

</decisions>

<specifics>
## Specific Ideas

- User finds TrajectoryFramework{T,C,H,PD,D} "ugly" and wants it simplified to 1-2 essential params
- Config domain param should be the primary dispatch mechanism: `f(config::Config{TimeDomain}, ...)` instead of `f(config, ::TimeDomain, ...)`
- If TrottTrott.bohr_freqs is renamed, the preferred name is `quasi_bohr_freqs`

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 08-struct-simplification*
*Context gathered: 2026-02-15*
