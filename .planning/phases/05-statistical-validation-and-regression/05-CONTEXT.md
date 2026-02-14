# Phase 5: Statistical Validation and Regression - Context

**Gathered:** 2026-02-14
**Status:** Ready for planning

<domain>
## Phase Boundary

Verify trajectory convergence properties (1/sqrt(N) scaling) and create frozen regression data for known-good DM and trajectory results. This phase delivers confidence that the validated code (Phases 1-4) stays correct over time.

</domain>

<decisions>
## Implementation Decisions

### Convergence test design
- 4-5 N_traj points in geometric progression (e.g., 1k, 4k, 16k, 64k)
- Ratio test: check error(N)/error(4N) is approximately 2.0 (consistent with Phase 4's delta scaling approach)
- Tight ratio bounds: [1.5, 2.5]
- Error metric: trace distance to DM result (not Gibbs state) — isolates trajectory sampling convergence from domain approximation
- Convergence tests are gated behind `QUANTUMFURNACE_FULL_TESTS=true` environment variable (expensive, skip by default)

### Regression data format
- Freeze both DM results and trajectory averages (fixed RNG seed via StableRNGs)
- File format: BSON (already a project dependency, used for Hamiltonian serialization)
- Storage location: `test/reference/` directory, committed to git
- Regression tolerance: ~1e-10 (allow for floating-point accumulation across Julia versions/platforms)
- Regression tests always run (fast — just load frozen file and compare against a fresh deterministic run)
- Trajectory regression: small fixed-seed average (~1000 trajectories) compared against frozen reference

### Domain coverage
- Convergence tests: EnergyDomain + TrotterDomain (with_coherent=true)
- Regression data: same scope — Energy + Trotter (with coherent)
- System size: 3-qubit Heisenberg only (reuse existing test fixtures)
- Rationale: Energy is simplest domain, Trotter with coherent is most complex — if both converge correctly, Time is implied

### Test runtime budget
- Total `Pkg.test()` target: under 5 minutes (including all phases 1-5)
- Convergence tests (expensive) behind env flag: `QUANTUMFURNACE_FULL_TESTS=true`
- Regression tests (cheap) always run — load BSON + compare, no trajectory averaging needed
- Pattern: `if get(ENV, "QUANTUMFURNACE_FULL_TESTS", "false") == "true"` to gate expensive tests

### Claude's Discretion
- Exact N_traj progression values (as long as 4-5 points, geometric progression)
- Number of DM steps for regression baseline
- Internal test organization (single file vs split)

</decisions>

<specifics>
## Specific Ideas

- Ratio test approach proven in Phase 4 (delta^2 scaling with ratio bounds [2.0, 8.0]) — adapt same pattern for 1/sqrt(N)
- StableRNGs already in test extras — use for reproducible trajectory regression data
- BSON serialization pattern already established for Hamiltonian objects — extend to density matrices

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 05-statistical-validation-and-regression*
*Context gathered: 2026-02-14*
