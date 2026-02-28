# Phase 38: Test Cleanup - Context

**Gathered:** 2026-02-28
**Status:** Ready for planning

<domain>
## Phase Boundary

Consolidate test infrastructure with parametrized helpers, informative output, and validated thresholds. Delete dead tests, integrate valuable orphaned tests. No new test coverage for new features — this phase cleans and strengthens what exists.

</domain>

<decisions>
## Implementation Decisions

### Helper consolidation
- Single `make_config(sim, domain; num_qubits=4, construction=KMS(), ...)` factory replacing all 8 current functions
- Single `make_test_system(; num_qubits=4, trotter=nothing)` replacing both `make_test_system()` and `make_small_test_system()`
- Keep precomputed globals (TEST_HAM, TEST_JUMPS, etc.) for both sizes — rename SMALL_* to N3_* for clarity
- Add `const ALL_DOMAINS = [EnergyDomain(), TimeDomain(), TrotterDomain(), BohrDomain()]` in test_helpers.jl

### @info output design
- Pattern: `@info "label" value=computed_value threshold=threshold_used` — placed AFTER the assertion
- Only for numerical comparisons (trace distance, eigenvalue, allocation bytes, convergence rate)
- Skip @info for structural/type checks (round-trip serialization, field existence, etc.)
- Every numerical @test shows: label, computed value, and threshold it was compared against

### Threshold review
- Audit ALL 239 threshold comparisons across all 22 test files
- Theory-based tightening first (O(delta), O(delta^2), machine epsilon, known error scaling)
- Empirical tightening acceptable but must account for scaling to larger systems (double qubit count) — small-system empirical values don't reliably bound larger systems
- Existing tiers (TOL_EXACT=1e-12, TOL_QUADRATURE=1e-6, TOL_DELTA(delta)=5*delta) are good — extend with new named constants where needed, but think through each case
- Rationale documented as inline comments next to each threshold check

### Old/staging test handling
- DELETE test/old_tests/ entirely (ham_test, kossakowski_test, trott_test, B_test, trajectory_test, time_tests)
- KEEP test/staging/ (test_fitting.jl, test_gap_estimation.jl) — matches dormant src/staging/ code
- INTEGRATE test/trajectory_validation/ into runtests.jl — behind a `FULL` test flag if tests are slow (few minutes)
- REVIEW test/reference/generate_references.jl — verify it still matches current types/API after v2.0 restructure

### Claude's Discretion
- Exact naming of new tolerance constants beyond the existing tiers
- Whether trajectory validation tests need the FULL flag (depends on measured runtime)
- How to implement the FULL test flag (ENV variable check)
- Grouping/ordering of @info output within test files

</decisions>

<specifics>
## Specific Ideas

- User emphasized: always show label + value + threshold for every numerical test, not just allocation tests
- User emphasized: theory-based tightening over empirical — "we will look at double amount of qubit numbers on the largest scale sim, so empirical tightening on smaller system tests can be off"
- Tolerance tiers are already well-designed; extension should be thoughtful, not mechanical

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 38-test-cleanup*
*Context gathered: 2026-02-28*
