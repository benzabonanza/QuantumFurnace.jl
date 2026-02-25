# Testing Patterns

**Analysis Date:** 2026-02-25

## Test Framework

**Runner:**
- Julia's built-in `Test` standard library
- Aqua.jl (`test/test_aqua.jl`) for package quality checks
- StableRNGs (`using StableRNGs`) used in `test/test_fitting.jl` for reproducible noise
- No config file (Julia test runner requires none); tests run via `Pkg.test()`

**Assertion Library:**
- `@test`, `@testset` from `Test`
- `@test isapprox(x, y; atol=tol)` for numerical comparisons
- `@test x isa Type` for type membership
- `@test allocs == 0` with `@allocated` macro for allocation checks
- `@elapsed` for timing/performance assertions

**Run Commands:**
```bash
julia --project -e 'using Pkg; Pkg.test()'   # Run all tests
julia --project test/runtests.jl              # Run directly
julia --threads 4 --project -e 'using Pkg; Pkg.test()'  # With threading
```

## Test File Organization

**Location:**
- Separate `test/` directory; NOT co-located with source files
- All test files included from `test/runtests.jl` via `include()`

**Naming:**
- `test_<feature>.jl` pattern: `test_cptp.jl`, `test_regression.jl`, `test_krylov_matvec.jl`
- Shared helpers in `test/test_helpers.jl` (included first by `runtests.jl`)
- Reference BSON data in `test/reference/` with generation script `test/reference/generate_references.jl`
- Old/archived tests in `test/old_tests/` (not included in `runtests.jl`)
- Validation scripts in `test/trajectory_validation/` (standalone, not part of `Pkg.test()`)

**Structure:**
```
test/
├── runtests.jl                    # Top-level: includes all test files
├── test_helpers.jl                # Shared fixtures, constants, factory functions
├── test_aqua.jl                   # Package quality (Aqua.jl)
├── test_compilation.jl            # Module load and fixture availability
├── test_cptp.jl                   # CPTP channel completeness
├── test_dm_detailed_balance.jl    # Gibbs fixed point and domain hierarchy
├── test_dm_scaling.jl             # Scaling with system size
├── test_regression.jl             # BSON and DM-based regression
├── test_allocation.jl             # Zero-allocation hot path guards
├── test_workspace_independence.jl # Thread workspace isolation
├── test_threading.jl              # Multi-thread determinism and speedup
├── test_gns_trajectory.jl         # GNS variant trajectory tests
├── test_results.jl                # Result serialization
├── test_convergence.jl            # Convergence tracking
├── test_fitting.jl                # Exponential decay fitting
├── test_observable_trajectories.jl # Observable-based trajectories
├── test_gap_estimation.jl         # Spectral gap estimation API
├── test_diagnostics.jl            # Exact diagnostics (DIAG-01 through DIAG-06)
├── test_krylov_matvec.jl          # Krylov matvec round-trips + allocations
├── test_krylov_eigsolve.jl        # Krylov eigsolve accuracy and guard rails
├── test_krylov_crossvalidation.jl # Krylov vs dense eigen cross-validation
├── reference/
│   ├── generate_references.jl     # One-time BSON reference generator
│   ├── energy_dm_reference.bson   # Frozen DM reference (EnergyDomain)
│   └── trotter_coherent_dm_reference.bson  # Frozen DM reference (TrotterDomain)
├── trajectory_validation/         # Standalone validation scripts (not in Pkg.test())
└── old_tests/                     # Archived tests (not run)
```

## Test Structure

**Suite Organization:**
```julia
@testset "DMTST-01: Bohr detailed balance (3-qubit)" begin
    config = LiouvConfig(...)
    liouv = construct_lindbladian(...)
    dist = trace_distance_h(...)
    @test dist < 1e-10
end

@testset "Krylov Matvec" begin
    @testset "Round-trip: matvec vs dense (EnergyDomain KMS, no coherent)" begin
        config = make_liouv_config(EnergyDomain(); with_coherent=false)
        # ...
        @test norm(v_dense - vec(L_rho)) < 1e-12
    end

    @testset "Zero allocations in matvec hot path" begin
        # Warmup first, then measure
        apply_lindbladian!(ws, rho, config, TEST_HAM)
        allocs = @allocated apply_lindbladian!(ws, rho, config, TEST_HAM)
        @test allocs == 0
    end
end
```

**Naming convention for test IDs:**
- Functional tests use coded names: `TVAL-01`, `DMTST-01`, `DMTST-02`, `TINF-02`, `TINF-03`, `CONV-01`, `DIAG-01`
- These codes appear in both the `@testset` label and in test file comments/docstrings

**Patterns:**
- Setup pattern: build config via factory → build framework/workspace → run computation
- Warmup pattern: call function once (JIT), then `@allocated` to measure
- No `@setup` / `@teardown` blocks; setup is inline at top of `@testset`
- No test fixtures in the xUnit sense; shared state is module-level constants

## Shared Fixtures (test_helpers.jl)

The most important pattern in this codebase: all expensive setup is computed **once at include time** and stored as module-level constants. Every test file reads from these constants via the shared `test_helpers.jl`.

**Tolerance tiers (LOCKED):**
```julia
const TOL_EXACT = 1e-12          # machine precision identities
const TOL_QUADRATURE = 1e-6      # quadrature / discretization errors
TOL_DELTA(delta) = 5.0 * delta   # unraveling error bound (function, not const)
```

**Shared system constants:**
```julia
const NUM_QUBITS = 4
const DIM = 2^NUM_QUBITS  # 16
const BETA = 10.0
const SIGMA = 1.0 / BETA
const TEST_DELTA = 0.01
```

**Shared fixtures (4-qubit system):**
```julia
const TEST_SYSTEM = make_test_system()
const TEST_HAM = TEST_SYSTEM.hamiltonian    # HamHam
const TEST_JUMPS = TEST_SYSTEM.jumps        # Vector{JumpOp} (12 jumps)
const TEST_GIBBS = TEST_SYSTEM.gibbs        # Gibbs state matrix
const TEST_TROTTER = make_test_trotter()
const TEST_TROTTER_JUMPS = make_test_system(; trotter=TEST_TROTTER).jumps
```

**Small fixtures (3-qubit, for fast tests):**
```julia
const SMALL_SYSTEM = make_small_test_system()
const SMALL_HAM = SMALL_SYSTEM.hamiltonian
const SMALL_JUMPS = SMALL_SYSTEM.jumps      # 9 jumps
const SMALL_GIBBS = SMALL_SYSTEM.gibbs
const SMALL_DIM = 2^3  # 8
const SMALL_TROTTER = make_small_test_trotter()
const SMALL_TROTTER_JUMPS = make_small_test_system(; trotter=SMALL_TROTTER).jumps
```

**Factory functions (for config construction):**
```julia
# 4-qubit configs
make_liouv_config(domain; with_coherent=true) -> LiouvConfig
make_liouv_config_gns(domain) -> LiouvConfigGNS
make_thermalize_config(domain; with_coherent=true, delta=TEST_DELTA, mixing_time=1.0) -> ThermalizeConfig

# 3-qubit (SMALL) configs
make_small_liouv_config(domain; with_coherent=false) -> LiouvConfig
make_small_liouv_config_gns(domain) -> LiouvConfigGNS
make_small_thermalize_config(domain; with_coherent=false, delta=TEST_DELTA, mixing_time=1.0) -> ThermalizeConfig
make_small_thermalize_config_gns(domain; delta=TEST_DELTA, mixing_time=1.0) -> ThermalizeConfigGNS
```

## Mocking

**Framework:** None (no mocking library used)

**Patterns:**
- No mocking of dependencies; tests use real implementations throughout
- Internal functions tested directly via `QuantumFurnace._internal_function(...)` qualified access
- Deterministic RNG seeding used for reproducibility: `Random.Xoshiro(seed)`, `MersenneTwister(42)`, `StableRNG(42)`

**What to Mock:**
- Nothing — the codebase favors integration-style testing over unit-mocked testing

**What NOT to Mock:**
- Everything; all tests exercise real matrix computations against each other

## Reference Data (BSON)

**Frozen BSON references** (`test/reference/`):
- DM evolution reference stored once with `generate_references.jl`, checked into git
- Loaded in `test_regression.jl` with `BSON.load(joinpath(ref_dir, "*.bson"))`
- Platform portability note: DM evolution (via `exp(delta*L)`) is deterministic; trajectory BSON is NOT used (platform-dependent due to BLAS/RNG differences)

```julia
# Generating a reference (run once, commit result)
BSON.bson(joinpath(REF_DIR, filename), Dict(:rho => Matrix(rho_dm), :delta => delta, ...))

# Loading in tests
ref_data = BSON.load(joinpath(ref_dir, "energy_dm_reference.bson"))
rho_ref = ref_data[:rho]
```

## Coverage

**Requirements:** Not enforced (no `.github/workflows/ci.yml` for tests; only `Documentation.yml` exists)

**View Coverage:**
```bash
# No automated coverage reporting configured
julia --project -e 'using Pkg; Pkg.test(; coverage=true)'
```

## Test Types

**Correctness Tests (majority):**
- Round-trip tests: compare matrix-free Krylov matvec against dense Liouvillian: `norm(v_dense - vec(L_rho)) < 1e-12`
- Symmetry/duality tests: `tr(X' * L(Y)) == tr(L*(X)' * Y)` (adjoint duality)
- Fixed-point tests: Gibbs state is exact fixed point of BohrDomain Lindbladian
- Domain hierarchy tests: error monotonically increases Bohr → Energy → Time → Trotter

**Regression Tests (`test_regression.jl`):**
- DM regression: fresh computation vs frozen BSON reference, `atol=1e-10`
- Trajectory regression: trajectory average vs fresh DM evolution, `atol=0.05` (accounts for statistical noise at ntraj=1000)

**Allocation Tests (`test_allocation.jl`, `test_krylov_matvec.jl`):**
- Hot-path zero-allocation guard: `@test @allocated f(...) == 0`
- Pattern: always warm up once before measuring to force JIT
- Use function barriers to avoid soft-scope boxing in `@testset`:
```julia
function _measure_matvec_allocs(ws, rho, config, ham)
    apply_lindbladian!(ws, rho, config, ham)  # warmup
    return @allocated apply_lindbladian!(ws, rho, config, ham)
end
```

**Threading Tests (`test_threading.jl`):**
- Determinism: `result1.rho_mean == result2.rho_mean` (bitwise identical, NOT `isapprox`)
- Seed isolation: per-trajectory seeds via `Xoshiro(seed + traj_id)`
- BLAS thread restoration: `@test BLAS.get_num_threads() == old_blas`
- Speedup test: `@test t_threaded < t_serial`
- All thread tests skip gracefully when `Threads.nthreads() == 1`: `@test true` placeholder

**Structural Tests (`test_compilation.jl`):**
- Type membership: `@test fw isa TrajectoryFramework`
- Field access: `@test fw.delta == TEST_DELTA`
- Fixture sizes: `@test size(TEST_HAM.data) == (DIM, DIM)`

## Common Patterns

**Domain coverage pattern** (run same test for all 4 domains):
```julia
for (name, domain) in [(:bohr, BohrDomain()), (:energy, EnergyDomain()),
                        (:time, TimeDomain()), (:trotter, TrotterDomain())]
    config = make_liouv_config(domain)
    trotter_obj = (domain isa TrotterDomain) ? TEST_TROTTER : nothing
    domain_jumps = (domain isa TrotterDomain) ? TEST_TROTTER_JUMPS : TEST_JUMPS
    # ... test body
end
```

**Random input sampling pattern** (statistical robustness):
```julia
for _ in 1:10
    rho = Matrix(random_density_matrix(NUM_QUBITS))
    v_dense = L_dense * vec(rho)
    L_rho = apply_lindbladian!(ws, rho, config, TEST_HAM)
    @test norm(v_dense - vec(L_rho)) < 1e-12
end
```

**Async/RNG seeding pattern:**
```julia
rng = Random.Xoshiro(seed)
# per-trajectory seeds for deterministic threading
rng_traj = Random.Xoshiro(seed + traj_id)
```

**Error Testing:**
```julia
# Constructor invariant
@test_throws ErrorException LiouvConfigGNS(...; with_coherent=true)
# Assertion violation
@test_throws AssertionError build_trajectoryframework(..., delta=2.0)
```

**Diagnostic logging in tests:**
```julia
@info "DMTST-02: Domain distances to Gibbs" bohr=distances[:bohr] energy=distances[:energy]
```

**`let` blocks for local scope in `@testset` (prevents variable boxing):**
```julia
let rng = MersenneTwister(42)
    raw_jump = randn(rng, ComplexF64, DIM, DIM) ./ sqrt(DIM)
    # ... test body using complex_jumps
end
```

---

*Testing analysis: 2026-02-25*
