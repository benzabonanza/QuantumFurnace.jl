# Phase 30: Cross-Validation - Context

**Gathered:** 2026-02-24
**Status:** Ready for planning

<domain>
## Phase Boundary

Validate Krylov spectral gap results against dense eigen() reference values across all domains and balance types. Establish trust that the Krylov method produces correct gaps for n>6 production use. No new capabilities — purely validation and cross-checking of existing Phase 29 eigensolver.

</domain>

<decisions>
## Implementation Decisions

### L-vs-E consistency testing
- Test multiple delta values (0.1, 0.01, 0.001) to demonstrate O(delta^2) convergence of channel-to-Lindbladian gap mapping
- Print a formatted convergence table: delta | gap_L | gap_from_E | error | order
- Hard assertion: convergence order must be >= 1.5 (test fails otherwise)
- Cover all 4 domains (EnergyDomain, TimeDomain, TrotterDomain, BohrDomain) with KMS balance
- L-vs-E convergence test runs for KMS only (channel math is balance-independent)

### GNS domain coverage
- GNS tested across all 4 domains at n=4
- GNS NOT tested at n=6 (n=4 is sufficient for balance-type correctness)
- KMS-vs-GNS comparison: each must match its own dense eigen() reference gap — the KMS/GNS relationship itself is not asserted
- L-vs-E consistency test runs KMS only, not GNS

### Diagnostic output
- On test failure: print top-k eigenvalue table from both Krylov and dense, with absolute/relative errors
- On test success: always print one-line summary per test (domain, gap_krylov, gap_dense, error)
- L-vs-E convergence: print formatted convergence table (delta, gap_L, gap_from_E, error, estimated order)
- No wall-clock timing in Phase 30 — correctness only; Phase 31 handles benchmarking

### Test organization
- New dedicated test file (e.g., test_krylov_crossvalidation.jl), separate from Phase 29 unit tests
- Single file with separate @testset blocks for n=4 and n=6
- Shared helper function (e.g., compare_krylov_dense) to run both methods and compare, reducing duplication across domain/balance combos
- n=6 tests gated behind `QUANTUMFURNACE_FULL_TESTS=true` environment variable — skipped in CI, run manually for validation

### Claude's Discretion
- Exact helper function signatures and internal structure
- How many top-k eigenvalues to print on failure
- @testset nesting structure within the file
- KrylovKit parameters (krylovdim, tol) for cross-validation runs

</decisions>

<specifics>
## Specific Ideas

- User wants `QUANTUMFURNACE_FULL_TESTS=true` env var to gate n=6 tests (not a custom test tag)
- Convergence table format should be human-readable, not just raw numbers
- Always-print summaries serve as a record of achieved accuracy for future reference

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 30-cross-validation*
*Context gathered: 2026-02-24*
