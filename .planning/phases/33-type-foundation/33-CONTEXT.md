# Phase 33: Type Foundation - Context

**Gathered:** 2026-02-25
**Status:** Ready for planning

<domain>
## Phase Boundary

Replace 4 duplicate config types (`LiouvConfig`, `ThermalizeConfig`, `LiouvConfigGNS`, `ThermalizeConfigGNS`) with a single `Config{S,D,C,T}` struct. Define simulation, domain, and construction singleton type hierarchies. Derive `with_coherent` from construction type at compile time. Update all call sites — no backward-compatible aliases.

</domain>

<decisions>
## Implementation Decisions

### Type parameter naming
- Struct: `Config{S,D,C,T}` — single-letter type params (S=Simulation, D=Domain, C=Construction, T=Float)
- Simulation hierarchy: `AbstractSimulation` with subtypes `Lindbladian`, `Thermalize`, `KrylovSpectrum`, `Trajectory`
- Domain hierarchy: keep existing `AbstractDomain` with subtypes `BohrDomain`, `EnergyDomain`, `TimeDomain`, `TrotterDomain`
- Construction hierarchy: `AbstractConstruction` with subtypes `KMS`, `GNS`, `DLL`
- All three singletons stored as both fields AND type parameters: `sim::S`, `domain::D`, `construction::C`

### Migration strategy
- Hard swap — remove old type names immediately, no aliases, no deprecation warnings
- Update all call sites in src/ and test/ in this phase
- Update success criteria to remove the alias requirement
- All 4 simulation types defined in Phase 33 (not deferred to Phase 36)
- Dispatch approach: use full `Config{S,D,C,T}` where methods need config fields; dispatch on extracted singletons (e.g., `::KMS`, `::Lindbladian`) where only the type tag matters — Claude's discretion per call site

### with_coherent derivation
- Trait function: `with_coherent(::KMS) = true`, `with_coherent(::GNS) = false`, `with_coherent(::DLL) = true` (placeholder)
- Remove `with_coherent` field from Config — derived from construction type
- Researcher should audit all KMS/GNS branching patterns beyond `with_coherent` (e.g., `pick_transition`, `precompute_data` already dispatch on construction type)
- DLL-specific behavior (no `with_linear_combinations`, no energy labels) deferred to future phase

### Constructor API
- Single constructor with all keyword arguments: `Config(sim=Lindbladian(), domain=EnergyDomain(), construction=KMS(), N=4, beta=1.0, ...)`
- No convenience constructors, no defaults — explicit is better
- No shortcut aliases

### Claude's Discretion
- Per-method decision on dispatching full Config vs extracted singleton
- Exact field ordering in Config struct
- How to handle the Bohr domain's different loop structure during migration
- Whether to introduce helper accessor functions (e.g., `simulation(config)`, `construction(config)`)

</decisions>

<specifics>
## Specific Ideas

- `pick_transition` and `precompute_data` already dispatch on construction type — good patterns to follow
- DLL won't use `with_linear_combinations` or energy labels, but that's a future concern — just define the singleton now
- 4 simulation types match the 4 `run_*` entry points planned for Phase 36

</specifics>

<deferred>
## Deferred Ideas

- DLL-specific behavior (no linear combinations, no energy labels) — future phase when DLL construction is implemented
- KrylovSpectrum and Trajectory run_* entry points — Phase 36 (types defined here, methods added there)

</deferred>

---

*Phase: 33-type-foundation*
*Context gathered: 2026-02-25*
