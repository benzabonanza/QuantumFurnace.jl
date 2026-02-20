# Phase 27: Core Matvec Infrastructure - Context

**Gathered:** 2026-02-20
**Status:** Ready for planning

<domain>
## Phase Boundary

Matrix-free Lindbladian action for EnergyDomain: `apply_lindbladian!` computes L(rho) without forming the full 4^n x 4^n superoperator matrix, validated against dense `construct_lindbladian()` at n=4. Includes coherent term and adjoint. KrylovWorkspace pre-allocates all scratch so the hot path is zero-allocation.

</domain>

<decisions>
## Implementation Decisions

### Adjoint interface
- Separate function `apply_adjoint_lindbladian!` (not keyword flag) — explicit, no ambiguity
- Shared KrylovWorkspace between forward and adjoint — same scratch matrices, used differently
- Leverage single-site Hermitian Pauli structure: A^a are Hermitian, so A^a' is lazy/free. Kraus jumps A(omega) are not Hermitian but have symmetry with negative omegas — use existing patterns from DM simulator code
- Memory-efficient: use lazy `'` adjoint where possible, don't store redundant orderings
- Design signatures now for KrylovKit compatibility — `apply_lindbladian!` should be easily wrappable into a closure `rho -> L(rho)` for KrylovKit's `eigsolve`

### Workspace scope
- KrylovWorkspace designed for all domains from the start (Energy, Time, Trotter, Bohr) — Phase 28 adds dispatch methods, not workspace restructuring
- Tied to (config, ham) at construction — precomputes everything once, one workspace per problem instance
- Scratch matrix pattern: follow existing DM simulator workspace pattern (TrajectoryWorkspace or equivalent)
- Jump operator storage: match DM simulator — dense matrices for A(omega) in EnergyDomain; for Time/TrotterDomain store NUFFT matrices and create A(omega) on the fly via entrywise multiplication

### Coherent term control
- `with_coherent` read from config object (not keyword arg) — config already carries this information
- Coherent correction B precomputed at workspace construction, using existing functions that precompute full B and its action (or B^a and its action)
- GNS balance: silently skip coherent term (no error if GNS config encountered)
- Adjoint coherent: automatic sign flip — `apply_adjoint_lindbladian!` uses +i[B, rho] instead of -i[B, rho], same B matrix. B is approximately Hermitian (up to Trotter and quadrature errors)

### Dispatch pattern
- Dispatch on domain type via config parametrization: `apply_lindbladian!(... config::SamplerConfig{EnergyDomain} ...)` — matches existing DM simulator dispatch
- KMS vs GNS dispatch: follow DM simulator pattern with different config names (one with GNS, one without)
- Inner loop order: match existing `construct_lindbladian` pattern (sites outer, frequencies inner) — consistency, easier to validate
- Output convention: determined by KrylovKit eigsolve requirements — research phase should check what signature KrylovKit expects and design accordingly

### Claude's Discretion
- Exact number of scratch matrices in the workspace pool
- Cache-friendliness optimizations within the loop (as long as loop order matches existing pattern)
- Internal helper function decomposition
- Test matrix generation strategy for round-trip validation

</decisions>

<specifics>
## Specific Ideas

- "Follow whatever pattern we had for the DM simulators" — the DM simulator workspace and dispatch patterns are the reference architecture for this phase
- Lazy `'` adjoint (Julia's `adjoint()`) should be used where possible to avoid storing explicit adjoints of Hermitian operators
- A(omega) for negative omega has symmetry with positive omega — exploit this from existing code rather than storing both
- Existing functions for precomputing B and its action should be reused, not reimplemented
- KrylovKit integration is a first-class design concern — the matvec signature should be closure-wrappable for Phase 29

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 27-core-matvec-infrastructure*
*Context gathered: 2026-02-20*
