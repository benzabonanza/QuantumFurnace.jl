# Phase 4: Trajectory Cross-Validation - Context

**Gathered:** 2026-02-14
**Status:** Ready for planning

<domain>
## Phase Boundary

Validate that trajectory-averaged density matrices match DM evolution for Energy, Time, and Trotter domains, and that TrotterDomain with coherent term converges toward the Gibbs state. This phase tests the trajectory simulation against the DM ground truth established in Phase 3.

</domain>

<decisions>
## Implementation Decisions

### Statistical matching criteria
- Fixed threshold comparison: trace distance between trajectory-averaged rho and DM rho must be < 0.01
- 10,000 trajectories per test point
- Fixed RNG seed via StableRNG for deterministic, reproducible results
- No empirical C/sqrt(N) bound fitting — that belongs in Phase 5

### System parameters
- 3-qubit Heisenberg system only (no 4-qubit)
- Reuse make_test_system() fixture from Phase 3 (same Hamiltonian, beta, sigma)
- Initial state: normalized all-ones pure state |psi> = (1,...,1)/sqrt(d), with DM starting from rho = |psi><psi|
- Single-step tests: all three domains (Energy, Time, Trotter) with a delta sweep to verify delta^2 scaling of trajectory-vs-DM error
- Multi-step convergence: TrotterDomain only, single delta value, with_coherent=true

### Coherent convergence test
- Both DM and trajectory-averaged rho must reach trace distance < 1e-3 to Gibbs state
- TrotterDomain with_coherent=true only (no Bohr baseline, no comparison with with_coherent=false)
- Run until 1e-3 threshold is reached for the given delta (not a fixed step count)
- Claude picks the delta value that achieves convergence practically

### Test runtime budget
- No delta sweep for multi-step tests — single delta value chosen by Claude
- Single-threaded execution (multi-threading deferred to a future milestone)
- Tests in a separate test group, not in main Pkg.test() suite
- Prioritize faster runtime over exhaustive parameter coverage

### Claude's Discretion
- Exact delta value for multi-step convergence test
- Delta values for the single-step delta sweep
- Number of multi-step thermalization steps (run until converged)
- Test file organization within the separate test group
- How to structure the delta^2 scaling check (number of delta values, regression method)

</decisions>

<specifics>
## Specific Ideas

- User noted that trajectory simulation starts from a pure state, so DM must start from |psi><psi| not I/d — the two simulators must use compatible initial states
- Multi-threading for trajectories is explicitly deferred to a future milestone rewrite
- 1e-6 Gibbs convergence target from roadmap was revised to 1e-3 due to delta-step error accumulation in both DM and trajectory modes

</specifics>

<deferred>
## Deferred Ideas

- Multi-threaded trajectory execution — future milestone
- 4-qubit system tests — not needed for current validation
- Empirical C/sqrt(N) trajectory scaling verification — Phase 5

</deferred>

---

*Phase: 04-trajectory-cross-validation*
*Context gathered: 2026-02-14*
