---
status: testing
phase: 37-file-organization-and-dead-code
source: 37-01-SUMMARY.md, 37-02-SUMMARY.md
started: 2026-02-27T16:10:00Z
updated: 2026-02-27T16:10:00Z
---

## Current Test

number: 1
name: Module loads cleanly
expected: |
  `using QuantumFurnace` succeeds without errors. No warnings about missing dependencies (Distributed, LsqFit, Optim, SharedArrays).
awaiting: user response

## Tests

### 1. Module loads cleanly
expected: `using QuantumFurnace` succeeds without errors. No warnings about missing dependencies (Distributed, LsqFit, Optim, SharedArrays).
result: [pending]

### 2. Old entry points removed
expected: `run_lindbladian` and `run_thermalization` are NOT accessible via `QuantumFurnace.run_lindbladian` / `QuantumFurnace.run_thermalization` (should throw UndefVarError). The new entry points `run_lindblad`, `run_thermalize`, `run_krylov_spectrum`, `run_trajectory` ARE exported and accessible.
result: [pending]

### 3. Dead structs removed
expected: `DMSimulationResult`, `LindbladianResult`, `LSIFramework` are NOT accessible (should throw UndefVarError). Active structs like `LindbladResults`, `ThermalizeResults`, `KrylovSpectrumResults`, `TrajectoryResults` ARE exported.
result: [pending]

### 4. Staging area established
expected: `src/staging/` contains `gap_estimation.jl`, `fitting.jl`, `log_sobolev.jl`. `test/staging/` contains `test_gap_estimation.jl`, `test_fitting.jl`. These files are NOT included in the module (not in `include()` calls in QuantumFurnace.jl) and NOT run by `test/runtests.jl`.
result: [pending]

### 5. Exports organized by simulation type
expected: The export block in `src/QuantumFurnace.jl` is organized into labeled sections (Lindbladian, Thermalize, Krylov, Trajectory, Diagnostics, Common) with `# ---` comment separators. Dead exports like `fit_exponential_decay`, `FitResult`, `estimate_spectral_gap` are removed. Dormant exports preserved as `# STAGING:` commented block.
result: [pending]

### 6. Krylov simulation script exists
expected: `simulations/main_krylov.jl` exists and contains a call to `run_krylov_spectrum`. All 4 simulation scripts exist: `main_liouv.jl`, `main_thermalize.jl`, `main_trajectory.jl`, `main_krylov.jl`.
result: [pending]

### 7. Test suite passes
expected: `julia --project test/runtests.jl` completes with all tests passing. Known pre-existing allocation threshold failures in test_krylov_matvec.jl (Time/TrotterDomain) are acceptable. No NEW failures introduced by Phase 37 changes.
result: [pending]

## Summary

total: 7
passed: 0
issues: 0
pending: 7
skipped: 0

## Gaps

[none yet]
