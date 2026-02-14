# Phase 5: Statistical Validation and Regression - Research

**Researched:** 2026-02-14
**Domain:** Statistical trajectory convergence (1/sqrt(N) scaling) and BSON-based regression testing
**Confidence:** HIGH

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

#### Convergence test design
- 4-5 N_traj points in geometric progression (e.g., 1k, 4k, 16k, 64k)
- Ratio test: check error(N)/error(4N) is approximately 2.0 (consistent with Phase 4's delta scaling approach)
- Tight ratio bounds: [1.5, 2.5]
- Error metric: trace distance to DM result (not Gibbs state) -- isolates trajectory sampling convergence from domain approximation
- Convergence tests are gated behind `QUANTUMFURNACE_FULL_TESTS=true` environment variable (expensive, skip by default)

#### Regression data format
- Freeze both DM results and trajectory averages (fixed RNG seed via StableRNGs)
- File format: BSON (already a project dependency, used for Hamiltonian serialization)
- Storage location: `test/reference/` directory, committed to git
- Regression tolerance: ~1e-10 (allow for floating-point accumulation across Julia versions/platforms)
- Regression tests always run (fast -- just load frozen file and compare against a fresh deterministic run)
- Trajectory regression: small fixed-seed average (~1000 trajectories) compared against frozen reference

#### Domain coverage
- Convergence tests: EnergyDomain + TrotterDomain (with_coherent=true)
- Regression data: same scope -- Energy + Trotter (with coherent)
- System size: 3-qubit Heisenberg only (reuse existing test fixtures)
- Rationale: Energy is simplest domain, Trotter with coherent is most complex -- if both converge correctly, Time is implied

#### Test runtime budget
- Total `Pkg.test()` target: under 5 minutes (including all phases 1-5)
- Convergence tests (expensive) behind env flag: `QUANTUMFURNACE_FULL_TESTS=true`
- Regression tests (cheap) always run -- load BSON + compare, no trajectory averaging needed
- Pattern: `if get(ENV, "QUANTUMFURNACE_FULL_TESTS", "false") == "true"` to gate expensive tests

### Claude's Discretion
- Exact N_traj progression values (as long as 4-5 points, geometric progression)
- Number of DM steps for regression baseline
- Internal test organization (single file vs split)

### Deferred Ideas (OUT OF SCOPE)
None -- discussion stayed within phase scope
</user_constraints>

## Summary

Phase 5 has two independent deliverables: (1) a statistical convergence test verifying that trajectory-averaged density matrix error to the DM reference decreases as 1/sqrt(N_traj), and (2) a regression test framework with BSON-frozen reference data for known-good DM and trajectory results.

The convergence test adapts the proven ratio-test approach from Phase 4 (where consecutive error ratios verified O(delta^2) scaling). Here, doubling sqrt(N) means quadrupling N_traj, so error(N)/error(4N) should be approximately 2.0. The test uses 4-5 N_traj points in geometric progression (factor-of-4 steps), computing trace distance between trajectory-averaged rho and exp(delta*L)*rho0 at each point. This test is expensive (the largest N_traj point may need 64k trajectories) and is gated behind the `QUANTUMFURNACE_FULL_TESTS=true` environment variable.

The regression test freezes known-good numerical results as BSON files in `test/reference/`. The DM regression data is fully deterministic (Liouvillian construction + matrix exponential). The trajectory regression requires a small fixed-seed trajectory average (~1000 trajectories with `Random.seed!(FIXED_SEED)`) that produces a deterministic result within a Julia session. The frozen reference is compared against a fresh run with the same seed, asserting elementwise agreement within 1e-10. Regression tests run in `Pkg.test()` (always-on, fast -- just load BSON and compare).

**Primary recommendation:** Implement Plan 05-01 (convergence test) as a gated test in `test/trajectory_validation/` alongside the existing Phase 4 tests. Implement Plan 05-02 (regression) as a new `test/test_regression.jl` included in `runtests.jl`, with a one-time generator script to create the frozen BSON files in `test/reference/`.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Test | stdlib | Julia test framework | Used throughout test suite |
| LinearAlgebra | stdlib | Matrix operations, eigen, exp, norm | DM reference computation |
| Random | stdlib | Global RNG seeding (`Random.seed!`) | Trajectory reproducibility |
| BSON | 0.3 | Serialize/deserialize reference data | Already a project dependency; used for Hamiltonian serialization; Dict-based save/load |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| StableRNGs | >= 1.0 | Version-stable RNG for regression seeds | In [extras]/test target; however, `step_along_trajectory!` uses global `rand()`, so StableRNG cannot be passed directly -- see RNG Pitfall below |

### Already Available (No Installation Needed)
All required libraries are already in Project.toml. BSON is in `[deps]`. StableRNGs is in `[extras]` and `[targets].test`. No new dependencies needed.

## Architecture Patterns

### Recommended Project Structure

```
test/
  runtests.jl                                  # Add include("test_regression.jl")
  test_helpers.jl                              # Existing fixtures (SMALL_*, make_small_*)
  test_regression.jl                           # NEW: regression tests (always-on, in Pkg.test())
  reference/                                   # NEW: frozen BSON reference data
    energy_dm_reference.bson                   # DM result for EnergyDomain (3-qubit)
    energy_traj_reference.bson                 # Trajectory average for EnergyDomain (3-qubit)
    trotter_coherent_dm_reference.bson         # DM result for TrotterDomain+coherent (3-qubit)
    trotter_coherent_traj_reference.bson       # Trajectory average for TrotterDomain+coherent (3-qubit)
  trajectory_validation/
    run_trajectory_validation.jl               # Existing Phase 4 tests
    run_convergence_tests.jl                   # NEW: 1/sqrt(N) convergence (gated behind env var)
```

### Pattern 1: 1/sqrt(N) Convergence Ratio Test

**What:** Run trajectory averages at increasing N_traj values, compute trace distance to DM reference, verify error(N)/error(4N) is approximately 2.0.

**When to use:** When `get(ENV, "QUANTUMFURNACE_FULL_TESTS", "false") == "true"`.

**Implementation:**

```julia
# DM reference (computed once)
liouv_config = make_small_liouv_config(EnergyDomain())
L = construct_lindbladian(SMALL_JUMPS, liouv_config, SMALL_HAM)
delta = 0.1  # single step, same as Phase 4
psi0 = fill(ComplexF64(1.0), SMALL_DIM) / sqrt(SMALL_DIM)
rho0 = psi0 * psi0'
rho_dm = reshape(exp(delta * L) * vec(rho0), SMALL_DIM, SMALL_DIM)
rho_dm = (rho_dm + rho_dm') / 2

# N_traj geometric progression (factor of 4)
ntraj_values = [1_000, 4_000, 16_000, 64_000]
errors = Float64[]

for ntraj in ntraj_values
    # Build framework (same as Phase 4 single_step_crossval)
    therm_config = make_small_thermalize_config(EnergyDomain(); delta=delta, mixing_time=Float64(delta))
    precomputed = precompute_data(EnergyDomain(), therm_config, SMALL_HAM)
    scratch = KrausScratch(ComplexF64, SMALL_DIM)
    fw = build_trajectoryframework(SMALL_JUMPS, SMALL_HAM, therm_config, precomputed, scratch, delta)

    rho_traj = zeros(ComplexF64, SMALL_DIM, SMALL_DIM)
    Random.seed!(42)  # Same seed each time for fair comparison
    for _ in 1:ntraj
        psi = copy(psi0)
        step_along_trajectory!(psi, fw)
        rho_traj .+= psi * psi'
    end
    rho_traj ./= ntraj
    rho_traj = (rho_traj + rho_traj') / 2

    dist = QuantumFurnace.trace_distance_h(Hermitian(rho_dm), Hermitian(rho_traj))
    push!(errors, dist)
end

# Ratio test: error(N)/error(4N) should be ~2.0 for 1/sqrt(N) scaling
ratios = [errors[i] / errors[i+1] for i in 1:length(errors)-1]
for ratio in ratios
    @test 1.5 <= ratio <= 2.5
end
```

**Key insight on RNG seeding for convergence tests:** Each N_traj point uses the same initial seed (`Random.seed!(42)`). The smaller-N runs use a prefix of the same random number sequence as the larger-N runs. This means the larger N_traj averages include the same trajectories as the smaller ones plus additional trajectories. This is actually ideal -- it reduces the variance of the ratio estimate because the first-N trajectories contribute identically to both points. The 1/sqrt(N) scaling test measures how additional trajectories reduce the remaining error.

### Pattern 2: BSON Reference Data Generation

**What:** Create frozen reference density matrices by running DM evolution and small trajectory averages, then save as BSON files.

**Implementation (one-time generator script):**

```julia
using BSON

# Example: generate EnergyDomain DM reference
rho_dm = ...  # compute as above
BSON.bson("test/reference/energy_dm_reference.bson", Dict(
    "rho" => rho_dm,
    "delta" => delta,
    "domain" => "EnergyDomain",
    "num_qubits" => 3,
    "beta" => BETA,
    "sigma" => SIGMA,
    "seed" => nothing,  # DM is deterministic
    "ntraj" => nothing,
    "description" => "Single-step DM reference via exp(delta*L)*rho0, EnergyDomain, 3-qubit Heisenberg"
))

# Example: generate EnergyDomain trajectory reference
Random.seed!(12345)  # Fixed seed for regression
ntraj_reg = 1000
rho_traj = ...  # compute as above
BSON.bson("test/reference/energy_traj_reference.bson", Dict(
    "rho" => rho_traj,
    "delta" => delta,
    "domain" => "EnergyDomain",
    "num_qubits" => 3,
    "beta" => BETA,
    "sigma" => SIGMA,
    "seed" => 12345,
    "ntraj" => ntraj_reg,
    "description" => "Trajectory average (1000 traj, seed=12345), EnergyDomain, 3-qubit Heisenberg"
))
```

### Pattern 3: BSON Reference Data Loading and Comparison

**What:** Load frozen reference, recompute fresh result, compare elementwise.

**Implementation (in test_regression.jl):**

```julia
@testset "TINF-02: DM regression (EnergyDomain)" begin
    # Load frozen reference
    source_root = dirname(@__DIR__)
    ref_path = joinpath(source_root, "test", "reference", "energy_dm_reference.bson")
    ref_data = BSON.load(ref_path)
    rho_ref = ref_data[:rho]

    # Recompute fresh
    liouv_config = make_small_liouv_config(EnergyDomain())
    L = construct_lindbladian(SMALL_JUMPS, liouv_config, SMALL_HAM)
    psi0 = fill(ComplexF64(1.0), SMALL_DIM) / sqrt(SMALL_DIM)
    rho0 = psi0 * psi0'
    rho_fresh = reshape(exp(ref_data[:delta] * L) * vec(rho0), SMALL_DIM, SMALL_DIM)
    rho_fresh = (rho_fresh + rho_fresh') / 2

    # Compare
    @test isapprox(rho_fresh, rho_ref; atol=1e-10)
end
```

### Pattern 4: Environment-Gated Expensive Tests

**What:** Skip expensive convergence tests unless explicitly enabled.

```julia
if get(ENV, "QUANTUMFURNACE_FULL_TESTS", "false") == "true"
    @testset "TVAL-05: 1/sqrt(N) convergence" begin
        # ... expensive convergence tests ...
    end
else
    @info "Skipping TVAL-05 convergence tests (set QUANTUMFURNACE_FULL_TESTS=true to run)"
end
```

### Anti-Patterns to Avoid

- **Comparing trajectory average against Gibbs state for convergence test:** This conflates domain approximation error with statistical sampling noise. Compare against the DM result (exp(delta*L)*rho0) which shares the same domain approximation and isolates the 1/sqrt(N) trajectory sampling noise.
- **Using different RNG seeds for different N_traj points:** Use the same seed so that the first N trajectories are identical across all points. This gives cleaner ratio estimates.
- **Storing reference data as JLD2 or HDF5:** BSON is already a project dependency and already used for Hamiltonian serialization. Adding another serialization library is unnecessary complexity.
- **Putting convergence tests in Pkg.test():** At 64k trajectories, the convergence test will take several minutes. It must be gated behind the environment variable.
- **Regenerating reference data in tests:** The reference data generation should be a separate one-time script. The regression test only loads and compares.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| BSON serialization | Custom binary format | `BSON.bson(path, dict)` / `BSON.load(path)` | Already used for Hamiltonians, handles ComplexF64 matrices natively |
| Trace distance computation | Manual eigenvalue SVD | `QuantumFurnace.trace_distance_h(Hermitian(A), Hermitian(B))` | Already verified in Phases 3-4 |
| Liouvillian construction | Custom matrix builder | `construct_lindbladian(jumps, config, ham; trotter=trotter)` | Verified in Phase 3, handles all domains correctly |
| Trajectory framework | Custom Kraus builder | `build_trajectoryframework(jumps, ham_or_trott, config, precomputed, scratch, delta)` | Verified CPTP in Phase 2, cross-validated in Phase 4 |
| Config factory (3-qubit) | Manual config construction | `make_small_liouv_config(domain)` / `make_small_thermalize_config(domain; ...)` | Added in Phase 4, already tested |

**Key insight:** Phase 5 builds entirely on verified infrastructure from Phases 1-4. No new mathematical operations are needed -- only new test patterns (ratio test for 1/sqrt(N), BSON save/load for regression).

## Common Pitfalls

### Pitfall 1: RNG Reproducibility for Trajectory Regression
**What goes wrong:** `step_along_trajectory!` uses `rand()` on the global RNG (lines 480, 514, 624, 658 in trajectories.jl). StableRNG cannot be passed directly.
**Why it happens:** The trajectory code was not designed with pluggable RNG. All random number generation goes through Julia's global/task-local default RNG.
**How to avoid:** Use `Random.seed!(FIXED_INTEGER_SEED)` before trajectory averaging. This seeds the global MersenneTwister. The trajectory output is deterministic within a Julia version + platform combination. Accept that results may differ between Julia versions (this is why regression tolerance is ~1e-10, not bitwise exact).
**Warning signs:** Regression test fails after Julia version upgrade. Resolution: regenerate reference data with the new Julia version.

**StableRNGs limitation:** Phase 2 discovered that `step_along_trajectory!` has no `rng` keyword argument. The Phase 4 tests used `Random.seed!(42)` successfully. For regression tests, use the same approach: `Random.seed!(FIXED_SEED)` gives within-version reproducibility. Cross-version reproducibility would require modifying `step_along_trajectory!` to accept an rng parameter (out of scope).

**Practical implication for the user's decision "fixed RNG seed via StableRNGs":** StableRNGs cannot be used directly with the trajectory code. The implementation should use `Random.seed!(FIXED_SEED)` for trajectory regression. For the DM regression, StableRNGs is not needed since DM evolution is fully deterministic (no RNG involved). The CONTEXT decision about StableRNGs should be interpreted as "use a fixed, documented seed for reproducibility" rather than literally requiring the StableRNGs library. The regression data is frozen as BSON files anyway, so cross-version RNG differences only matter when regenerating the reference data.

### Pitfall 2: Convergence Test Noise Floor
**What goes wrong:** At very large N_traj, the 1/sqrt(N) statistical error may approach numerical precision limits, causing the ratio to deviate from 2.0.
**Why it happens:** For N=64k trajectories on a dim=8 system, the trace distance to DM is very small. Machine epsilon effects and Hermitianization rounding can create a noise floor.
**How to avoid:** Choose N_traj progression so the largest point still has error well above machine precision. With delta=0.1 and 3-qubit system, Phase 4 showed that 50k trajectories give trace distance ~0.001-0.003. At 64k, the error should be ~0.001, well above 1e-15 machine precision.
**Warning signs:** The largest-N error is suspiciously small (< 1e-6) or the final ratio deviates strongly from 2.0.

### Pitfall 3: BSON Path Resolution During Pkg.test()
**What goes wrong:** `Pkg.test()` runs tests from a temporary directory, so relative paths to `test/reference/` fail.
**Why it happens:** Julia's `Pkg.test()` copies the test files to a temporary environment. Relative paths break because the working directory is not the project root.
**How to avoid:** Use `@__DIR__` to resolve paths relative to the test file location, then navigate to the reference directory:
```julia
source_root = dirname(@__DIR__)
ref_path = joinpath(source_root, "test", "reference", "energy_dm_reference.bson")
```
This pattern is already established in `test_helpers.jl` (line 57-58) for loading Hamiltonians:
```julia
source_root = dirname(@__DIR__)
ham_path = joinpath(source_root, "hamiltonians", "heis_disordered_periodic_n$(NUM_QUBITS).bson")
```
**Warning signs:** "File not found" errors during `Pkg.test()` that don't occur when running the test file directly.

### Pitfall 4: BSON Serialization of Hermitian Matrices
**What goes wrong:** BSON may not correctly round-trip `Hermitian{ComplexF64, Matrix{ComplexF64}}` wrapper types.
**Why it happens:** BSON serializes Julia types including their wrapper structure. Loading may reconstruct a different type.
**How to avoid:** Always save plain `Matrix{ComplexF64}`, not `Hermitian(...)`. Convert when loading:
```julia
# Save:
BSON.bson(path, Dict("rho" => Matrix(rho_dm)))

# Load:
ref_data = BSON.load(path)
rho_ref = ref_data[:rho]  # Plain Matrix{ComplexF64}
# Compare with:
@test isapprox(rho_fresh, rho_ref; atol=1e-10)
```
**Warning signs:** Type errors or unexpected data shapes when loading BSON files.

### Pitfall 5: Convergence Test Delta Value Selection
**What goes wrong:** Using a delta value where the systematic error (trajectory-vs-DM channel mismatch at order delta^2) is comparable to the statistical error, making the 1/sqrt(N) signal noisy.
**Why it happens:** The trajectory-vs-DM trace distance has two components: (i) the finite-N statistical error (scales as 1/sqrt(N)) and (ii) the channel mismatch at O(delta^2). If (ii) dominates, increasing N_traj won't improve the error, breaking the ratio test.
**How to avoid:** Use a delta value where the O(delta^2) systematic error is negligible compared to the 1/sqrt(N) statistical error at the smallest N_traj point. From Phase 4 data: at delta=0.1, the systematic channel mismatch (measured at N=50k where statistical noise is small) is of order 10^-4. The statistical error at N=1000 should be ~0.01 (sqrt(50) times larger than at 50k). So delta=0.1 gives a good separation: systematic ~10^-4 << statistical ~10^-2.
**Warning signs:** Ratios converge to a value significantly less than 2.0, indicating a systematic error floor.

### Pitfall 6: TrotterDomain Requires trotter Keyword
**What goes wrong:** `construct_lindbladian` and `build_trajectoryframework` for TrotterDomain require the `trotter` keyword argument.
**Why it happens:** TrotterDomain uses Trotter quasi-Bohr frequencies and Trotter eigenvectors, which live in the TrottTrott object.
**How to avoid:** Always pass `trotter=SMALL_TROTTER` when using TrotterDomain. Already well-established in Phase 4 code:
```julia
trotter_kw = domain isa TrotterDomain ? (; trotter=SMALL_TROTTER) : (;)
L = construct_lindbladian(SMALL_JUMPS, liouv_config, SMALL_HAM; trotter_kw...)
```
**Warning signs:** "TrotterDomain requires `trotter`" error.

## Code Examples

### Example 1: Complete 1/sqrt(N) Convergence Test

```julia
# Source: adapted from Phase 4 single_step_crossval + Phase 4 ratio test pattern

function convergence_test(domain; with_coherent::Bool = (domain isa TrotterDomain), delta::Float64 = 0.1)
    dim = SMALL_DIM  # 8

    # DM reference (computed once)
    liouv_config = make_small_liouv_config(domain; with_coherent=with_coherent)
    trotter_kw = domain isa TrotterDomain ? (; trotter=SMALL_TROTTER) : (;)
    L = construct_lindbladian(SMALL_JUMPS, liouv_config, SMALL_HAM; trotter_kw...)

    psi0 = fill(ComplexF64(1.0), dim) / sqrt(dim)
    rho0 = psi0 * psi0'
    rho_dm = reshape(exp(delta * L) * vec(rho0), dim, dim)
    rho_dm = (rho_dm + rho_dm') / 2

    # Trajectory framework (built once, reused across N_traj points)
    therm_config = make_small_thermalize_config(domain;
        with_coherent=with_coherent, delta=delta, mixing_time=Float64(delta))
    ham_or_trott = domain isa TrotterDomain ? SMALL_TROTTER : SMALL_HAM
    precomputed = precompute_data(domain, therm_config, ham_or_trott)
    scratch = KrausScratch(ComplexF64, dim)
    fw = build_trajectoryframework(SMALL_JUMPS, ham_or_trott, therm_config,
        precomputed, scratch, delta)

    # N_traj geometric progression (factor-of-4 steps)
    ntraj_values = [1_000, 4_000, 16_000, 64_000]
    errors = Float64[]

    for ntraj in ntraj_values
        rho_traj = zeros(ComplexF64, dim, dim)
        Random.seed!(42)
        for _ in 1:ntraj
            psi = copy(psi0)
            step_along_trajectory!(psi, fw)
            rho_traj .+= psi * psi'
        end
        rho_traj ./= ntraj
        rho_traj = (rho_traj + rho_traj') / 2

        dist = QuantumFurnace.trace_distance_h(Hermitian(rho_dm), Hermitian(rho_traj))
        push!(errors, dist)
        println("  N_traj=$ntraj  trace_dist=$(round(dist; sigdigits=4))")
    end

    # Ratio test
    ratios = [errors[i] / errors[i+1] for i in 1:length(errors)-1]
    return errors, ratios
end
```

### Example 2: BSON Reference Data Generation Script

```julia
# Source: codebase BSON.load pattern from test_helpers.jl + generate_hamiltonians.jl

"""
Generate frozen reference data for regression tests.
Run once: julia --project test/reference/generate_references.jl
"""

using QuantumFurnace
using LinearAlgebra
using Random
using BSON

# Include test helpers for fixtures
include(joinpath(@__DIR__, "..", "test_helpers.jl"))

const REF_DIR = @__DIR__
const REF_DELTA = 0.1
const REF_SEED = 12345
const REF_NTRAJ = 1000

function generate_dm_reference(domain; with_coherent::Bool, filename::String)
    dim = SMALL_DIM
    liouv_config = make_small_liouv_config(domain; with_coherent=with_coherent)
    trotter_kw = domain isa TrotterDomain ? (; trotter=SMALL_TROTTER) : (;)
    L = construct_lindbladian(SMALL_JUMPS, liouv_config, SMALL_HAM; trotter_kw...)

    psi0 = fill(ComplexF64(1.0), dim) / sqrt(dim)
    rho0 = psi0 * psi0'
    rho_dm = reshape(exp(REF_DELTA * L) * vec(rho0), dim, dim)
    rho_dm = (rho_dm + rho_dm') / 2

    BSON.bson(joinpath(REF_DIR, filename), Dict(
        "rho" => Matrix(rho_dm),
        "delta" => REF_DELTA,
        "domain" => string(typeof(domain)),
        "with_coherent" => with_coherent,
        "num_qubits" => 3,
    ))
    println("Saved: $filename")
end

function generate_traj_reference(domain; with_coherent::Bool, filename::String)
    dim = SMALL_DIM
    therm_config = make_small_thermalize_config(domain;
        with_coherent=with_coherent, delta=REF_DELTA, mixing_time=Float64(REF_DELTA))
    ham_or_trott = domain isa TrotterDomain ? SMALL_TROTTER : SMALL_HAM
    precomputed = precompute_data(domain, therm_config, ham_or_trott)
    scratch = KrausScratch(ComplexF64, dim)
    fw = build_trajectoryframework(SMALL_JUMPS, ham_or_trott, therm_config,
        precomputed, scratch, REF_DELTA)

    psi0 = fill(ComplexF64(1.0), dim) / sqrt(dim)
    rho_traj = zeros(ComplexF64, dim, dim)
    Random.seed!(REF_SEED)
    for _ in 1:REF_NTRAJ
        psi = copy(psi0)
        step_along_trajectory!(psi, fw)
        rho_traj .+= psi * psi'
    end
    rho_traj ./= REF_NTRAJ
    rho_traj = (rho_traj + rho_traj') / 2

    BSON.bson(joinpath(REF_DIR, filename), Dict(
        "rho" => Matrix(rho_traj),
        "delta" => REF_DELTA,
        "domain" => string(typeof(domain)),
        "with_coherent" => with_coherent,
        "num_qubits" => 3,
        "seed" => REF_SEED,
        "ntraj" => REF_NTRAJ,
    ))
    println("Saved: $filename")
end

# Generate all reference files
generate_dm_reference(EnergyDomain(); with_coherent=false, filename="energy_dm_reference.bson")
generate_traj_reference(EnergyDomain(); with_coherent=false, filename="energy_traj_reference.bson")
generate_dm_reference(TrotterDomain(); with_coherent=true, filename="trotter_coherent_dm_reference.bson")
generate_traj_reference(TrotterDomain(); with_coherent=true, filename="trotter_coherent_traj_reference.bson")
```

### Example 3: Regression Test (Always-On, in Pkg.test())

```julia
# Source: test_helpers.jl BSON.load pattern + Phase 3/4 trace_distance_h usage

@testset "TINF-02: Regression tests" begin
    source_root = dirname(@__DIR__)
    ref_dir = joinpath(source_root, "test", "reference")

    @testset "DM regression: EnergyDomain" begin
        ref_data = BSON.load(joinpath(ref_dir, "energy_dm_reference.bson"))
        rho_ref = ref_data[:rho]
        delta = ref_data[:delta]

        liouv_config = make_small_liouv_config(EnergyDomain())
        L = construct_lindbladian(SMALL_JUMPS, liouv_config, SMALL_HAM)
        psi0 = fill(ComplexF64(1.0), SMALL_DIM) / sqrt(SMALL_DIM)
        rho0 = psi0 * psi0'
        rho_fresh = reshape(exp(delta * L) * vec(rho0), SMALL_DIM, SMALL_DIM)
        rho_fresh = (rho_fresh + rho_fresh') / 2

        @test isapprox(rho_fresh, rho_ref; atol=1e-10)
    end

    @testset "Trajectory regression: EnergyDomain" begin
        ref_data = BSON.load(joinpath(ref_dir, "energy_traj_reference.bson"))
        rho_ref = ref_data[:rho]
        delta = ref_data[:delta]
        seed = ref_data[:seed]
        ntraj = ref_data[:ntraj]

        therm_config = make_small_thermalize_config(EnergyDomain();
            delta=delta, mixing_time=Float64(delta))
        precomputed = precompute_data(EnergyDomain(), therm_config, SMALL_HAM)
        scratch = KrausScratch(ComplexF64, SMALL_DIM)
        fw = build_trajectoryframework(SMALL_JUMPS, SMALL_HAM, therm_config,
            precomputed, scratch, delta)

        psi0 = fill(ComplexF64(1.0), SMALL_DIM) / sqrt(SMALL_DIM)
        rho_traj = zeros(ComplexF64, SMALL_DIM, SMALL_DIM)
        Random.seed!(seed)
        for _ in 1:ntraj
            psi = copy(psi0)
            step_along_trajectory!(psi, fw)
            rho_traj .+= psi * psi'
        end
        rho_traj ./= ntraj
        rho_traj = (rho_traj + rho_traj') / 2

        @test isapprox(rho_traj, rho_ref; atol=1e-10)
    end
end
```

## Discretion Recommendations

### N_traj Progression Values

**Recommendation:** `ntraj_values = [1_000, 4_000, 16_000, 64_000]`

Rationale:
- 4 points give 3 ratio checks, matching Phase 4's proven approach (3 delta values, 2 ratios)
- Factor-of-4 geometric progression means sqrt(4N)/sqrt(N) = 2, so each ratio should be ~2.0
- N=1,000 gives trace distance ~0.01 (well above machine precision noise floor)
- N=64,000 gives trace distance ~0.001 (still well above machine precision, matches Phase 4's 50k calibration)
- Total trajectory count across all 4 points: 85,000. At dim=8 and ~microseconds per step, this is ~10 seconds per domain. Two domains (Energy + Trotter): ~20 seconds total. Acceptable for gated test.

### Number of DM Steps for Regression Baseline

**Recommendation:** Single step with delta=0.1 (matching Phase 4 single-step cross-validation)

Rationale:
- Single-step DM via `exp(delta * L) * vec(rho0)` is fully deterministic and fast
- delta=0.1 gives non-trivial evolution (rho changes significantly from rho0)
- Already validated to be correct in Phase 4 TVAL-02/04
- Simple to verify: one matrix exponential, no iteration

### Internal Test Organization

**Recommendation:** Two files -- split by gating behavior:

1. `test/test_regression.jl` -- Regression tests, always-on, included in `runtests.jl`
   - Contains TINF-02 regression tests (DM + trajectory for both domains)
   - Fast: loads 4 BSON files, runs 2 DM computations + 2x1000-trajectory averages
   - Estimated runtime: ~5-10 seconds

2. `test/trajectory_validation/run_convergence_tests.jl` -- Convergence tests, gated
   - Contains TVAL-05 convergence test (1/sqrt(N) scaling for both domains)
   - Expensive: runs up to 64k trajectories per domain
   - Gated: only runs with `QUANTUMFURNACE_FULL_TESTS=true`
   - NOT included in runtests.jl
   - Run via: `QUANTUMFURNACE_FULL_TESTS=true julia --project test/trajectory_validation/run_convergence_tests.jl`

Rationale for split:
- Regression tests must always run (user decision) -> include in `Pkg.test()` via runtests.jl
- Convergence tests are expensive (user decision: gated) -> separate file, not in Pkg.test()
- Keeps the convergence tests alongside the Phase 4 trajectory validation tests (natural grouping)
- Alternative considered: single file with env-gated sections. Rejected because the convergence test is fundamentally a different test mode (gated) and should not be loaded during normal Pkg.test().

### Trajectory Regression Seed and Count

**Recommendation:** seed=12345, ntraj=1000

Rationale:
- 1000 trajectories at dim=8 takes ~0.01 seconds (negligible for regression test runtime)
- Seed 12345 is distinct from Phase 4 test seeds (42, 123, 999) to avoid any confusion
- 1000 trajectories produce a meaningful (non-trivial) density matrix that exercises all branches of step_along_trajectory!
- The resulting rho_traj is deterministic given `Random.seed!(12345)` on the same Julia version/platform

### Reference Data Delta Value

**Recommendation:** delta=0.1

Rationale:
- Matches Phase 4 single-step tests where delta=0.1 was used
- Large enough to produce non-trivial evolution (the DM reference is measurably different from rho0)
- Small enough that delta_eff = delta * n_jumps = 0.1 * 9 = 0.9 < 1.0 (required assertion in build_trajectoryframework)
- Single-step means the reference captures the core CPTP channel behavior without accumulated multi-step errors

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| No convergence rate verification | Ratio test for 1/sqrt(N) scaling | Phase 5 (this phase) | Proves trajectory averaging has correct statistical properties |
| No regression data | BSON-frozen reference files | Phase 5 (this phase) | Catches silent regressions in numerical output |
| Random.seed! for test reproducibility | Same (StableRNG not usable with trajectory code) | Phase 2 (discovered) | Regression data tied to Julia version's MersenneTwister implementation |

**Deprecated/outdated:**
- StableRNG for trajectory seeding: not possible since `step_along_trajectory!` uses global `rand()`. Use `Random.seed!(integer)` instead.
- `run_trajectories` for regression: use manual trajectory loop for fine-grained seed control (run_trajectories calls validate_config! and print_press which produce stdout noise).

## Open Questions

1. **BSON Round-Trip Fidelity for ComplexF64 Matrices**
   - What we know: BSON is used successfully for HamHam objects (which contain `Matrix{ComplexF64}` fields like `data` and `eigvecs`). The Hamiltonian BSON files in `hamiltonians/` load correctly.
   - What's unclear: Whether BSON preserves full Float64 bit-level precision for ComplexF64 matrices, or if there is any rounding during serialization.
   - Recommendation: Test BSON round-trip on a simple ComplexF64 matrix before committing reference data. If any precision loss occurs, the 1e-10 regression tolerance should absorb it. Verify with: `M = rand(ComplexF64, 8, 8); BSON.bson("/tmp/test.bson", Dict("m" => M)); M2 = BSON.load("/tmp/test.bson")[:m]; @assert M == M2`.
   - **Confidence:** HIGH that this works, since HamHam objects with ComplexF64 matrices round-trip successfully (the entire test suite depends on loading Hamiltonians from BSON).

2. **Julia Version Sensitivity of Random.seed! Output**
   - What we know: `Random.seed!(12345)` seeds the global MersenneTwister. Within a Julia version, the output is deterministic. Between Julia major versions (1.11 -> 1.12), the MersenneTwister implementation may change.
   - What's unclear: Whether the trajectory regression test will need regenerating when upgrading Julia.
   - Recommendation: Document the Julia version used to generate reference data in the BSON metadata. If regression tests fail after a Julia upgrade, regenerate reference data. The 1e-10 tolerance applies to within-version floating-point accumulation, not cross-version RNG differences. If cross-version RNG differences are found, the reference files must be regenerated.

3. **Convergence Test Runtime at N=64k**
   - What we know: At dim=8, each step_along_trajectory! involves 8x8 matrix-vector products. Phase 4 ran 50k trajectories in seconds.
   - What's unclear: Exact wall time for 64k trajectories x 2 domains.
   - Recommendation: Based on Phase 4 timing (~40 seconds for the entire suite including 50k-trajectory tests), 64k single-step trajectories at dim=8 should take ~2-5 seconds. Two domains: ~4-10 seconds. Well within the gated test budget.

## Sources

### Primary (HIGH confidence)
- Codebase direct analysis: `src/trajectories.jl` (rand() usage lines 480, 514, 624, 658), `src/furnace.jl` (construct_lindbladian, run_trajectories), `src/qi_tools.jl` (trace_distance_h), `src/furnace_utensils.jl` (precompute_data for all domains)
- Existing test infrastructure: `test/test_helpers.jl` (SMALL_*, make_small_*, BSON path resolution via dirname(@__DIR__)), `test/trajectory_validation/run_trajectory_validation.jl` (single_step_crossval pattern, Phase 4 ratio test)
- Phase 4 summaries: `.planning/phases/04-trajectory-cross-validation/04-01-SUMMARY.md` (50k trajectories, noise floor ~0.001-0.003), `.planning/phases/04-trajectory-cross-validation/04-02-SUMMARY.md` (thermalization validation)
- Phase 2 summary: `.planning/phases/02-trajectory-bug-fixes/02-01-SUMMARY.md` (StableRNG limitation discovered: "step_along_trajectory! uses global rand() with no rng keyword argument")
- Project.toml: BSON in [deps] (v0.3), StableRNGs in [extras] (v1)
- Hamiltonian BSON pattern: `hamiltonians/generate_hamiltonians.jl` (BSON.@save), `test/test_helpers.jl` (BSON.load)
- Phase 5 CONTEXT.md: locked decisions on convergence design, regression format, domain coverage, runtime budget

### Secondary (MEDIUM confidence)
- Phase 4 calibration data: 50k trajectories at delta=0.1 give trace distance ~0.001-0.003 (from 04-01-SUMMARY log-log slopes). Used to estimate 1/sqrt(N) scaling curve: at N=1000, expect ~0.007*sqrt(50) ~= 0.05; at N=64000, expect ~0.007*sqrt(50/64) ~= 0.006.
- BSON ComplexF64 fidelity: inferred from successful Hamiltonian serialization (eigvecs, data fields are Matrix{ComplexF64})

### Tertiary (LOW confidence)
- Convergence test runtime estimates: based on Phase 4 timing extrapolation, not direct measurement

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - all libraries already in Project.toml, no new deps, BSON patterns verified from existing usage
- Architecture: HIGH - all patterns adapted from Phase 4's proven single_step_crossval and ratio test approach; BSON save/load follows established Hamiltonian pattern
- Pitfalls: HIGH - RNG limitation confirmed from Phase 2 discovery (code inspection of rand() calls in trajectories.jl); BSON path resolution from existing test_helpers.jl pattern; all other pitfalls from direct codebase analysis
- Code examples: HIGH - adapted from verified Phase 4 code patterns with minimal modification

**Research date:** 2026-02-14
**Valid until:** 2026-03-14 (stable domain, no external dependency changes expected)
