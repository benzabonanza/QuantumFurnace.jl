---
status: testing
phase: 36-api-and-results
source: 36-01-SUMMARY.md, 36-02-SUMMARY.md, 36-03-SUMMARY.md, 36-04-SUMMARY.md
started: 2026-02-27T13:00:00Z
updated: 2026-02-27T13:00:00Z
---

## Current Test
<!-- OVERWRITE each test - shows where we are -->

number: 1
name: Full test suite passes
expected: |
  Running `julia --project -e 'using Pkg; Pkg.test()'` completes with all 1254 tests passing and 0 errors/failures.
awaiting: user response

## Tests

### 1. Full test suite passes
expected: Running `julia --project -e 'using Pkg; Pkg.test()'` completes with all 1254 tests passing and 0 errors/failures.
result: [pending]

### 2. New entry points are exported and callable
expected: In a Julia REPL with `using QuantumFurnace`, all 4 new entry points (`run_lindblad`, `run_thermalize`, `run_krylov_spectrum`, `run_trajectory`) are accessible as functions (not UndefVarError).
result: [pending]

### 3. New Result types are exported
expected: `LindbladResults`, `ThermalizeResults`, `KrylovSpectrumResults`, `TrajectoryResults`, and `AbstractResults` are all accessible types after `using QuantumFurnace`. `save_result` and `load_result` are also exported.
result: [pending]

### 4. Round-trip serialization for Result types
expected: Constructing a LindbladResults (or any Result type), calling `save_result(result, path)`, then `load_result(path)` returns an equivalent struct with matching fields (config, eigenvalues, metadata, etc.).
result: [pending]

### 5. Companion .txt file generated on save
expected: After `save_result(result, "test.bson")`, a companion file `test.txt` exists alongside it containing a human-readable summary with result type name and domain info.
result: [pending]

### 6. Simulation scripts use new API
expected: `simulations/main_liouv.jl` calls `run_lindblad` (not `run_lindbladian`), `simulations/main_thermalize.jl` calls `run_thermalize` (not `run_thermalization`), and `simulations/main_krylov_benchmark.jl` references `run_krylov_spectrum`.
result: [pending]

### 7. run_trajectory keyword-driven mode dispatch
expected: `run_trajectory` accepts keyword arguments to switch between default, observable, convergence, and adaptive modes. Calling with no mode keywords runs default batch; providing `observables=` runs observable mode; providing `convergence=true` runs convergence mode.
result: [pending]

## Summary

total: 7
passed: 0
issues: 0
pending: 7
skipped: 0

## Gaps

[none yet]
