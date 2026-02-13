# Phase 1: Foundation and Compilation - Context

**Gathered:** 2026-02-13
**Status:** Ready for planning

<domain>
## Phase Boundary

Trajectory code compiles and test infrastructure exists for all subsequent phases. This phase delivers: Project.toml cleanup (TINF-04), shared test helpers with fixtures and tolerances (TINF-01), and the `build_trajectoryframework` compilation fix (TFIX-01). No bug fixes beyond compilation, no validation tests, no DM tests.

</domain>

<decisions>
## Implementation Decisions

### Test system parameters
- 4-qubit Heisenberg chain (dim=16) as the primary test system
- Load pre-computed Hamiltonian from `hamiltonians/heis_disordered_periodic_n4.bson` (optimized disorder realization via `find_ideal_heisenberg`)
- beta = 10 (inverse temperature)
- `with_linear_combination=true` with `a=beta/30`, `b=0.4` (specific transition function for paper simulations)
- `with_coherent` controls exact detailed balance: true gives exact Gibbs in Bohr domain; Energy/Time/Trotter still have approximation errors regardless
- All 4 domain configs (Bohr, Energy, Time, Trotter) for Lindbladian construction and spectral analysis tests
- Only Energy, Time, Trotter configs for thermalization (DM and trajectory) — no Bohr thermalization (Bohr trajectories are out of scope)
- Single-site Pauli jump operators (X, Y, Z on each site), normalized by `sqrt(3 * num_qubits)`

### Tolerance tier design
- Three tiers matching the physics error sources:
  - `TOL_EXACT` ~ 1e-12: Machine precision identities (CPTP checks, Bohr detailed balance, algebraic identities)
  - `TOL_QUADRATURE` ~ 1e-6: Quadrature, discretization, and Trotter approximation errors (domain comparisons to exact/Bohr reference)
  - `TOL_DELTA(delta) = C * delta`: Unraveling error from finite step size in weak measurement scheme, scaling linearly with delta (C to be calibrated empirically)
- `TOL_EXACT` and `TOL_QUADRATURE` are hardcoded constants; `TOL_DELTA` is a function of the step size parameter

### Test suite organization
- Standard Julia pattern: single `runtests.jl` that includes test files with `@testset` blocks
- Shared helpers in `test/test_helpers.jl`, included by `runtests.jl` before test files
- New test files named `test_<feature>.jl` (e.g., `test_compilation.jl`, `test_cptp.jl`, `test_dm_evolution.jl`)
- Existing 7 test files kept as-is for reference (development scripts with Revise); new proper tests written fresh alongside them
- `Pkg.test()` runs the full suite through `runtests.jl`

### Fixture scope
- `make_test_system` computed once at include time, stored as top-level constants in `test_helpers.jl`
- Includes precomputed Gibbs state (`exp(-beta*H)/Z`) for trace distance comparisons
- Includes a standard `TEST_DELTA` step size value for thermalization tests (individual tests can override)
- Depth of precomputation: Claude's discretion based on minimizing boilerplate in common test patterns

### Claude's Discretion
- Exact `TEST_DELTA` value (calibrate for good test sensitivity)
- Calibration constant C for `TOL_DELTA(delta)`
- Exact depth of `make_test_system` precomputation (how far down the Ham -> jumps -> config -> precomputed_data chain to go)
- Grid parameters for domain configs (`num_energy_bits`, `w0`, `t0`, `num_trotter_steps_per_t0`) — should be consistent with existing test values
- The specific compilation fix for TFIX-01 (undefined `trotter` variable, uninitialized `B_total`)
- Project.toml cleanup details (which deps move to `[extras]`)

</decisions>

<specifics>
## Specific Ideas

- Use the same `a=beta/30`, `b=0.4` parameter values from existing `trajectory_test.jl` — these define the transition function for the paper's main simulations
- Existing pre-computed BSON Hamiltonians are already optimized for large minimum Bohr frequency gap — leverage this property for cleaner test signals
- The `with_coherent=true` vs `false` distinction is important to showcase (manifests around delta=0.01) — fixture should support both configurations

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 01-foundation-and-compilation*
*Context gathered: 2026-02-13*
