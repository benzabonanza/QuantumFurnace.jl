# Phase 14: GNS Trajectory Path - Context

**Gathered:** 2026-02-16
**Status:** Ready for planning

<domain>
## Phase Boundary

GNS (approximate detailed balance) trajectory simulation works end-to-end and produces physically correct results. Verify that `ThermalizeConfigGNS` dispatches through the full trajectory pipeline (workspace allocation, stepping, density matrix accumulation) and that the averaged trajectory DM converges to the GNS Lindbladian fixed point. This phase validates the existing GNS code path -- it does not build new GNS capabilities or compare KMS vs GNS (that's Phase 18).

</domain>

<decisions>
## Implementation Decisions

### GNS reference state
- Compute the GNS fixed point using the existing `run_lindbladian()` with a GNS config, which returns `HotSpectralResults` (fixed point, spectral gap, first excited mode)
- No separate null-space solver needed -- `run_lindbladian()` already provides this via `construct_lindbladian()` dispatching on config type
- The GNS Lindbladian is built by `construct_lindbladian()` when given a `LiouvConfigGNS` -- this is already implemented
- Also compute and document the trace distance between the GNS fixed point and the exact Gibbs state (the approximation gap) as a baseline for Phase 18

### Sigma parameter design
- Sigma already exists as a field on GNS config structs
- Phase 14 validation tests use sigma = 1/beta only
- The two-sigma comparison (1/beta and 0.5/beta) is Phase 18 scope
- Test at n=3 (dim=8) only -- small system for fast tests, matches existing test fixtures
- beta = 10.0 minimum for all simulations -- lower beta is high temperature and can mask convergence errors

### Approximation tolerance
- Use trace distance (`trace_distance_h()`) as the convergence metric
- Target: trace distance < 0.05 between averaged trajectory DM and GNS fixed point
- Adjust delta step size (try 0.1 or 0.01) and trajectory count to achieve this threshold
- Single delta value sufficient -- no delta sweep needed for Phase 14
- Validate DM properties (Hermitian, unit trace, PSD) on the final averaged result only
- If final result fails validation, escalate to batch checkpoint checks to isolate where it breaks

### B-term handling
- `with_coherent=false` is already the default for GNS configs -- no B-term construction in GNS runs
- `step_along_trajectory!` already respects `with_coherent=false` and skips the unitary step
- `pick_transition()` dispatches on GNS config type and uses GNS-specific transition weights gamma(omega)
- No explicit B-term absence tests needed -- the GNS Lindbladian converges to its own fixed point regardless of B presence
- The entire GNS code path (jump construction, B-term suppression, transition selection) is config-driven and already implemented

### Claude's Discretion
- Exact delta step size and trajectory count to achieve < 0.05 trace distance
- Test structure and organization (single test file vs integrated into existing test suite)
- Whether to test BohrDomain, EnergyDomain, or TrotterDomain for GNS (or multiple)

</decisions>

<specifics>
## Specific Ideas

- Use `run_lindbladian()` with GNS config to get the reference -- gives spectral gap and first excited mode for free, not just the fixed point
- Beta < 10 is considered too easy (high temperature) and can mask errors -- enforce beta >= 10 in all simulation tests
- "Everything is already working well for GNS case" -- the primary task is verification and testing, not building new functionality

</specifics>

<deferred>
## Deferred Ideas

- Two-sigma comparison (sigma=1/beta vs sigma=0.5/beta) demonstrating that reducing sigma improves GNS approximation -- Phase 18
- Larger system sizes (n=4,6) for GNS validation -- Phase 18 parameter grid
- Delta sweep showing convergence improvement with smaller steps -- future work if needed

</deferred>

---

*Phase: 14-gns-trajectory-path*
*Context gathered: 2026-02-16*
