# Phase 20: Observable Infrastructure - Context

**Gathered:** 2026-02-17
**Status:** Ready for planning

<domain>
## Phase Boundary

Build total magnetization and combined gap estimation observables in both Hamiltonian eigenbasis and Trotter eigenbasis. Follows the existing `build_convergence_observables` pattern from v1.2. Exponential fitting, trajectory running, and gap estimation API are separate phases.

</domain>

<decisions>
## Implementation Decisions

### Observable selection
- Start with only H (energy) and M_z (total magnetization per site) for gap estimation
- No ZZ correlations in the gap estimation bundle — keep it lean
- Rationale: choose observables that overlap with the first excited mode (slowest-decaying) to capture the spectral gap
- If H + M_z don't give a good gap estimate (cross-checked against exact Liouvillian), reconsider adding more observables later

### M_z definition
- Per-site magnetization: M_z = sum(Z_i) / n (not total sum)
- Normalizes across system sizes, making amplitude comparison easier for n=4 vs n=6 cross-validation
- Decay rate (which determines spectral gap) is unaffected by normalization

### Trotter support
- Both `build_total_magnetization` and `build_gap_estimation_observables` must have Trotter variants
- The Trotter spectral gap (algorithmic evolution) is the primary quantity of interest for the paper
- Pattern: `build_total_magnetization(ham, n; trotter=trotter)` and `build_gap_estimation_observables(ham, n; trotter=trotter)`

### Return format
- Keep existing `(observables::Vector{Matrix{ComplexF64}}, names::Vector{String})` tuple pattern
- Consistent with `build_convergence_observables` — downstream code already handles this

### Claude's Discretion
- Whether `build_gap_estimation_observables` internally calls `build_total_magnetization` or constructs directly
- Test system sizes and tolerance values for regression tests
- Internal helper organization

</decisions>

<specifics>
## Specific Ideas

- "We are curious really what the spectral gap of the real algorithmic evolution is, i.e. the Trotter version, for the paper"
- The idea is to pick observables that for sure have the first excited mode in them to capture the slowest decaying part
- If the initial H + M_z pair doesn't produce good gap estimates, the plan is to reconsider and add more observables (not pre-optimize the bundle)

</specifics>

<deferred>
## Deferred Ideas

- ZZ correlations and other observables for gap estimation — revisit if H + M_z insufficient (after Phase 24 cross-validation)

</deferred>

---

*Phase: 20-observable-infrastructure*
*Context gathered: 2026-02-17*
