---
status: complete
phase: 28-domain-extension
source: 28-01-SUMMARY.md, 28-02-SUMMARY.md
started: 2026-02-24T00:00:00Z
updated: 2026-02-24T00:01:00Z
---

## Current Test

[testing complete]

## Tests

### 1. Full test suite passes
expected: All tests pass (1095+ assertions) with no errors or failures when running the full test suite.
result: pass

### 2. TimeDomain forward matvec matches dense construction
expected: `apply_lindbladian!(ws, rho, config, ham)` for a TimeDomain KMS config at n=4 produces output matching `construct_lindbladian(config, ham) * vec(rho)` (reshaped) to < 1e-12 norm error, for a random Hermitian PSD density matrix.
result: pass

### 3. TrotterDomain forward matvec matches dense construction
expected: `apply_lindbladian!(ws, rho, config, ham)` for a TrotterDomain KMS config at n=4 produces output matching dense Lindbladian to < 1e-12 norm error.
result: pass

### 4. BohrDomain forward matvec matches dense construction
expected: `apply_lindbladian!(ws, rho, config, ham)` for a BohrDomain KMS config at n=4 produces output matching dense Lindbladian to < 1e-12 norm error, using bucket iteration over Bohr frequencies with the two-operator dissipator.
result: pass

### 5. Adjoint Lindbladian matches for all new domains
expected: `apply_adjoint_lindbladian!` for TimeDomain, TrotterDomain, and BohrDomain each match their corresponding dense adjoint `construct_lindbladian(config, ham)'` to < 1e-12 norm error at n=4.
result: pass

### 6. BohrDomain forward/adjoint duality holds
expected: For random Hermitian matrices X, Y: `|tr(X' * L(Y)) - tr(L*(X)' * Y)| < 1e-11` confirming L and L* are Hilbert-Schmidt adjoints. This validates the BohrDomain two-operator dissipator adjoint formula.
result: pass

### 7. Zero heap allocations in hot paths
expected: `apply_lindbladian!` and `apply_adjoint_lindbladian!` for TimeDomain and TrotterDomain produce zero heap allocations (after warmup) when measured with `@allocated`. BohrDomain is allowed one buffer allocation per call (A_nu2_dag scatter buffer).
result: pass

## Summary

total: 7
passed: 7
issues: 0
pending: 0
skipped: 0

## Gaps

[none yet]
