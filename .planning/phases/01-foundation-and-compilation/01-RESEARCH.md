# Phase 1: Foundation and Compilation - Research

**Researched:** 2026-02-13
**Domain:** Julia package testing infrastructure, Project.toml management, compilation bug fixes
**Confidence:** HIGH

## Summary

Phase 1 has three concrete deliverables: (1) clean up Project.toml to separate dev deps from production deps, (2) create test infrastructure with shared fixtures and tiered tolerances, and (3) fix compilation bugs in `build_trajectoryframework`. All three are well-understood, concrete tasks with clear success criteria.

The compilation bugs in `build_trajectoryframework` are fully diagnosed: an undefined `trotter` variable on line 53, an uninitialized `B_total` on line 96 (when `with_coherent=false`), a dangling `where T` type parameter on line 48, and a spurious `trotter=trotter` keyword passed to `precompute_coherent_total_B` which no longer accepts it. Additionally, the export block in `QuantumFurnace.jl` lines 37-38 has a missing comma causing `precompute_data` and `verify_completeness` to not be exported. The same `trotter=trotter` kwarg bug also exists in `furnace.jl` line 66 (but is not blocking because `trotter` is defined in that scope as a function parameter -- it will error at runtime if called, not at compile time).

The test infrastructure follows standard Julia patterns: a `test/runtests.jl` entry point, a `test/test_helpers.jl` with shared fixtures, and test-specific dependencies in `[extras]`/`[targets]`. The 4-qubit Heisenberg chain test system is well-defined by user decisions, with pre-computed Hamiltonians already available in `hamiltonians/heis_disordered_periodic_n4.bson`.

**Primary recommendation:** Fix all identified bugs in `trajectories.jl` and `QuantumFurnace.jl` (export comma), create the test scaffolding from scratch (new `runtests.jl`, `test_helpers.jl`, `test_compilation.jl`), and restructure Project.toml. Keep existing 7 test files untouched as development scripts.

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- 4-qubit Heisenberg chain (dim=16) as the primary test system
- Load pre-computed Hamiltonian from `hamiltonians/heis_disordered_periodic_n4.bson` (optimized disorder realization via `find_ideal_heisenberg`)
- beta = 10 (inverse temperature)
- `with_linear_combination=true` with `a=beta/30`, `b=0.4` (specific transition function for paper simulations)
- `with_coherent` controls exact detailed balance: true gives exact Gibbs in Bohr domain; Energy/Time/Trotter still have approximation errors regardless
- All 4 domain configs (Bohr, Energy, Time, Trotter) for Lindbladian construction and spectral analysis tests
- Only Energy, Time, Trotter configs for thermalization (DM and trajectory) -- no Bohr thermalization (Bohr trajectories are out of scope)
- Single-site Pauli jump operators (X, Y, Z on each site), normalized by `sqrt(3 * num_qubits)`
- Three tolerance tiers: `TOL_EXACT` ~ 1e-12, `TOL_QUADRATURE` ~ 1e-6, `TOL_DELTA(delta) = C * delta`
- `TOL_EXACT` and `TOL_QUADRATURE` are hardcoded constants; `TOL_DELTA` is a function of step size
- Standard Julia pattern: single `runtests.jl` that includes test files with `@testset` blocks
- Shared helpers in `test/test_helpers.jl`, included by `runtests.jl` before test files
- New test files named `test_<feature>.jl` (e.g., `test_compilation.jl`, `test_cptp.jl`, `test_dm_evolution.jl`)
- Existing 7 test files kept as-is for reference; new proper tests written fresh alongside them
- `Pkg.test()` runs the full suite through `runtests.jl`
- `make_test_system` computed once at include time, stored as top-level constants in `test_helpers.jl`
- Includes precomputed Gibbs state (`exp(-beta*H)/Z`) for trace distance comparisons
- Includes a standard `TEST_DELTA` step size value for thermalization tests (individual tests can override)

### Claude's Discretion
- Exact `TEST_DELTA` value (calibrate for good test sensitivity)
- Calibration constant C for `TOL_DELTA(delta)`
- Exact depth of `make_test_system` precomputation (how far down the Ham -> jumps -> config -> precomputed_data chain to go)
- Grid parameters for domain configs (`num_energy_bits`, `w0`, `t0`, `num_trotter_steps_per_t0`) -- should be consistent with existing test values
- The specific compilation fix for TFIX-01 (undefined `trotter` variable, uninitialized `B_total`)
- Project.toml cleanup details (which deps move to `[extras]`)

### Deferred Ideas (OUT OF SCOPE)
None -- discussion stayed within phase scope
</user_constraints>

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Test (stdlib) | Julia stdlib | `@test`, `@testset` macros | Standard Julia testing framework |
| BSON | (already in deps) | Load pre-computed Hamiltonian fixtures | Already used by the project |
| LinearAlgebra (stdlib) | Julia stdlib | Matrix operations in test helpers | Already used throughout |

### Test-Only Dependencies (to add to [extras])
| Library | Purpose | When to Use |
|---------|---------|-------------|
| StableRNGs | Reproducible RNG streams across Julia versions | Any test requiring random numbers (trajectory sampling) |
| HypothesisTests | Statistical tests for convergence verification | Phase 5 trajectory convergence tests |
| StatsBase | `mean_and_std` for trajectory averages | Phase 4-5 trajectory statistics |
| Aqua | Package quality checks (ambiguities, stale deps) | Phase 3 quality gate |

### Dependencies to Move from [deps] to [extras]
| Library | Current Location | Reason |
|---------|-----------------|--------|
| BenchmarkTools | `[deps]` AND `[extras]` | Dev-only; already partially in extras |
| Debugger | `[deps]` AND `[extras]` | Dev-only; already partially in extras |
| Revise | `[deps]` | Dev-only; used only in interactive sessions |
| Plots | `[deps]` | Dev-only; not needed for package functionality |
| ClusterManagers | `[deps]` | Only for multi-node compute; not needed for core functionality |

**Note:** BenchmarkTools and Debugger are currently duplicated -- they appear in both `[deps]` (lines 9, 12) and `[extras]` (lines 48-49). They must be removed from `[deps]` while keeping them in `[extras]`.

**Installation:** No `npm install` equivalent -- Julia resolves dependencies from Project.toml. After editing Project.toml, run:
```julia
] resolve
] test  # to verify test deps resolve correctly
```

## Architecture Patterns

### Recommended Test Structure
```
test/
  runtests.jl            # Entry point for Pkg.test(); includes helpers then test files
  test_helpers.jl        # Shared fixtures, constants, utility functions
  test_compilation.jl    # Phase 1: compilation smoke tests (TFIX-01)
  # Future phases add:
  # test_cptp.jl         # Phase 2: CPTP verification
  # test_dm_evolution.jl # Phase 3: DM reference tests
  # test_trajectories.jl # Phase 4: trajectory cross-validation
  # Existing dev scripts (kept as-is):
  B_test.jl
  ham_test.jl
  kossakowski_test.jl
  log_sobolev_test.jl
  time_tests.jl
  trajectory_test.jl
  trott_test.jl
```

### Pattern 1: runtests.jl Entry Point
**What:** Single file that loads helpers and includes all test files
**When to use:** Always -- this is the Julia convention for `Pkg.test()`
**Example:**
```julia
using Test
using QuantumFurnace
using LinearAlgebra

# Shared fixtures and constants
include("test_helpers.jl")

@testset "QuantumFurnace.jl" begin
    include("test_compilation.jl")
    # Future phases add test includes here
end
```

### Pattern 2: test_helpers.jl Fixture Module
**What:** Computed-once constants providing the standard test system
**When to use:** Shared across all test files via top-level constants
**Example:**
```julia
using BSON

# === Physical parameters (locked decisions) ===
const NUM_QUBITS = 4
const DIM = 2^NUM_QUBITS  # 16
const BETA = 10.0
const SIGMA = 1.0 / BETA  # 0.1

# === Tolerance tiers ===
const TOL_EXACT = 1e-12        # Machine precision identities
const TOL_QUADRATURE = 1e-6    # Quadrature/discretization errors
TOL_DELTA(delta) = 5.0 * delta # Unraveling error (C=5.0, calibrate empirically)

# === Test step size ===
const TEST_DELTA = 0.01  # Small enough for good test sensitivity

# === Grid parameters (consistent with existing test values) ===
const NUM_ENERGY_BITS = 12
const W0 = 0.05
const T0 = 2pi / (2^NUM_ENERGY_BITS * W0)
const NUM_TROTTER_STEPS_PER_T0 = 10

# === Precomputed test system ===
function make_test_system()
    # Load pre-computed Hamiltonian
    hamiltonian = load_hamiltonian("heis", NUM_QUBITS)
    hamiltonian = finalize_hamham(hamiltonian, BETA)

    # Jump operators: single-site Paulis, normalized
    jump_normalization = sqrt(3 * NUM_QUBITS)
    jumps = JumpOp[]
    for pauli in [[X], [Y], [Z]]
        for site in 1:NUM_QUBITS
            jump_op = Matrix(pad_term(pauli, NUM_QUBITS, site)) / jump_normalization
            jump_in_eigen = hamiltonian.eigvecs' * jump_op * hamiltonian.eigvecs
            orthogonal = (jump_op == transpose(jump_op))
            hermitian = (jump_op == jump_op')
            push!(jumps, JumpOp(jump_op, jump_in_eigen, orthogonal, hermitian))
        end
    end

    # Gibbs state for trace distance comparisons
    gibbs = hamiltonian.gibbs  # Already computed by finalize_hamham

    return (; hamiltonian, jumps, gibbs)
end

# Compute once at include time
const TEST_SYSTEM = make_test_system()
const TEST_HAM = TEST_SYSTEM.hamiltonian
const TEST_JUMPS = TEST_SYSTEM.jumps
const TEST_GIBBS = TEST_SYSTEM.gibbs
```

### Pattern 3: Domain Config Factory Functions
**What:** Helper functions to create configs for each domain, reducing boilerplate
**When to use:** Whenever a test needs a domain-specific config
**Example:**
```julia
function make_liouv_config(domain; with_coherent=true)
    LiouvConfig(
        num_qubits = NUM_QUBITS,
        with_coherent = with_coherent,
        with_linear_combination = true,
        domain = domain,
        beta = BETA,
        sigma = SIGMA,
        a = BETA / 30.0,
        b = 0.4,
        num_energy_bits = NUM_ENERGY_BITS,
        w0 = W0,
        t0 = T0,
        eta = 0.0,
        num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
    )
end

function make_thermalize_config(domain; with_coherent=true, delta=TEST_DELTA, mixing_time=1.0)
    ThermalizeConfig(
        num_qubits = NUM_QUBITS,
        with_coherent = with_coherent,
        with_linear_combination = true,
        domain = domain,
        beta = BETA,
        sigma = SIGMA,
        a = BETA / 30.0,
        b = 0.4,
        num_energy_bits = NUM_ENERGY_BITS,
        w0 = W0,
        t0 = T0,
        eta = 0.0,
        num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
        mixing_time = mixing_time,
        delta = delta,
    )
end
```

### Anti-Patterns to Avoid
- **Including existing dev test files in runtests.jl:** These use `Revise` + `includet` and are not `Pkg.test()`-compatible. Leave them untouched and write new test files.
- **Using `using .QuantumFurnace` (relative import):** In `Pkg.test()`, the package is loaded as a full module with `using QuantumFurnace` (no dot). The existing dev scripts use `.QuantumFurnace` because they `includet` the source directly.
- **Defining constants inside `@testset` blocks:** Constants should be at module top-level in `test_helpers.jl` so all test files can access them.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Reproducible RNG | Custom seed management | `StableRNGs.StableRNG(seed)` | Streams are stable across Julia versions |
| Package quality checks | Manual export/ambiguity inspection | `Aqua.test_all(QuantumFurnace)` | Catches stale deps, undefined exports, ambiguities |
| Statistical convergence tests | Manual mean/std calculation | `StatsBase.mean_and_std`, `HypothesisTests` | Correct numerical stability and statistical rigor |
| Test dependency isolation | Manual path manipulation | `[extras]` + `[targets]` in Project.toml | Julia's built-in test dependency mechanism |

## Common Pitfalls

### Pitfall 1: Duplicate deps in [deps] and [extras]
**What goes wrong:** BenchmarkTools and Debugger currently appear in both `[deps]` and `[extras]`. This causes them to be loaded for production use AND inflates install time for users.
**Why it happens:** Incremental Project.toml editing without cross-checking.
**How to avoid:** Remove from `[deps]`, keep in `[extras]`. Verify with `] resolve`.
**Warning signs:** `] status` shows packages in both sections.

### Pitfall 2: Missing comma in export block silently fails
**What goes wrong:** In `QuantumFurnace.jl` lines 37-38:
```julia
export ..., run_trajectories, precompute_R
       precompute_data, verify_completeness
```
Line 38 is a no-op expression, not an export continuation. `precompute_data` and `verify_completeness` are never exported. Note: `verify_completeness` is actually commented out in `jump_workers.jl` line 611, so exporting it would cause an error (undefined).
**Why it happens:** Multi-line export statements without trailing comma.
**How to avoid:** Fix the comma on line 37. Remove `verify_completeness` from the export list (it's commented-out dead code). Ensure `precompute_data` is properly exported.
**Warning signs:** `QuantumFurnace.precompute_data` works but `precompute_data` doesn't (without module prefix).

### Pitfall 3: `where T` dangling type parameter
**What goes wrong:** `build_trajectoryframework` has `where T` at line 48 but no type parameter `T` is used in the function signature. Julia permits this but it's misleading and suggests incomplete refactoring.
**Why it happens:** The function was refactored to use concrete `ComplexF64` types but the `where T` was not removed.
**How to avoid:** Remove `where T` from the function signature.

### Pitfall 4: trotter keyword passed to function that doesn't accept it
**What goes wrong:** Both `trajectories.jl:53` and `furnace.jl:66` call `precompute_coherent_total_B(...; trotter=trotter)` but the function signature in `coherent.jl:14-19` does NOT have a `trotter` keyword parameter. Julia will error: "MethodError: no method matching precompute_coherent_total_B(... ; trotter=...)".
**Why it happens:** The function signature was refactored (the old version accepted `trotter` as a kwarg) but call sites were not updated.
**How to avoid:** Remove `;trotter=trotter` from both call sites. The function already receives `ham_or_trott` as a positional argument which can be either `HamHam` or `TrottTrott`.

### Pitfall 5: B_total undefined when with_coherent=false
**What goes wrong:** In `build_trajectoryframework`, when `config.with_coherent == false`, the else branch (line 56-58) sets `U_B = nothing` but does not set `B_total`. Then line 96 passes `B_total` to the constructor, causing `UndefVarError: B_total not defined`.
**Why it happens:** The `TrajectoryFramework` struct has a `B` field that needs a value regardless of whether coherent is enabled.
**How to avoid:** Set `B_total = nothing` in the else branch alongside `U_B = nothing`.

### Pitfall 6: Existing test files use Revise/includet pattern
**What goes wrong:** All 7 existing test files start with `using Revise` and `includet("../src/QuantumFurnace.jl")`. These are NOT compatible with `Pkg.test()` which loads the package via standard `using QuantumFurnace`.
**Why it happens:** These were written as interactive development scripts, not proper unit tests.
**How to avoid:** Do NOT include them in `runtests.jl`. Write new test files from scratch that use `using QuantumFurnace` (provided by the runtests.jl entry point).

### Pitfall 7: BohrDomain config needs precomputed bohr_dict
**What goes wrong:** `precompute_data(::BohrDomain, config::AbstractThermalizeConfig, hamiltonian)` accesses `hamiltonian.bohr_dict` (line 42 of `furnace_utensils.jl`). If `finalize_hamham` was not called, `bohr_dict` is `nothing` and the code crashes.
**Why it happens:** `HamHam` constructors set `bohr_dict = nothing`; only `finalize_hamham` populates it.
**How to avoid:** The `make_test_system` fixture calls `finalize_hamham(hamiltonian, BETA)` which populates both `bohr_dict` and `gibbs`. Ensure this is always done before creating BohrDomain configs.

## Code Examples

### Compilation Fix for build_trajectoryframework (TFIX-01)

The fix involves four changes in `src/trajectories.jl`:

```julia
# Line 48: Remove dangling `where T`
# BEFORE:
#   delta::Float64) where T
# AFTER:
    delta::Float64)

# Lines 52-58: Fix undefined trotter and uninitialized B_total
# BEFORE:
#   if config.with_coherent
#       B_total = precompute_coherent_total_B(jumps, ham_or_trott, config, precomputed_data; trotter=trotter)
#       B_total .= 0.5 .* (B_total .+ B_total')
#       U_B = exp(-1im * delta * Hermitian(B_total))
#   else
#       U_B = nothing
#   end
# AFTER:
    if config.with_coherent
        B_total = precompute_coherent_total_B(jumps, ham_or_trott, config, precomputed_data)
        B_total .= 0.5 .* (B_total .+ B_total')
        U_B = exp(-1im * delta * Hermitian(B_total))
    else
        B_total = nothing
        U_B = nothing
    end
```

And in `src/furnace.jl` line 66, the same keyword removal:
```julia
# BEFORE:
#   Btot = precompute_coherent_total_B(jumps, ham_or_trott, config, precomputed_data; trotter=trotter)
# AFTER:
    Btot = precompute_coherent_total_B(jumps, ham_or_trott, config, precomputed_data)
```

### Export Block Fix in QuantumFurnace.jl

```julia
# Lines 36-38: Fix missing comma and remove dead export
# BEFORE:
# export TrajectoryFramework, TrajectoryWorkspace, build_trajectoryframework, step_along_trajectory!, evolve_along_trajectory,
#        evolve_and_measure_along_trajectory, run_trajectories, precompute_R
#        precompute_data, verify_completeness
# AFTER:
export TrajectoryFramework, TrajectoryWorkspace, build_trajectoryframework, step_along_trajectory!, evolve_along_trajectory,
       evolve_and_measure_along_trajectory, run_trajectories, precompute_R,
       precompute_data
# Note: verify_completeness removed -- it's commented-out dead code in jump_workers.jl
```

### Project.toml Restructuring

```toml
[deps]
Arpack = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
BSON = "fbb218c0-5317-5bc6-957e-2ee96dd4b1f0"
DataStructures = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
Distributed = "8ba89e20-285c-5b6f-9357-94700520ee1b"
FINUFFT = "d8beea63-0952-562e-9c6a-8e8ef7364055"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
LinearMaps = "7a12625a-238d-50fd-b39a-03d52299707e"
Optim = "429524aa-4258-5aef-a3af-852621145aeb"
Pkg = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
ProgressMeter = "92933f4c-e287-5a05-a399-4b506db050ca"
QuadGK = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Roots = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
SharedArrays = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

# Removed from [deps]: BenchmarkTools, ClusterManagers, Debugger, Plots, Revise

[extras]
Aqua = "4c88cf16-eb10-579e-8560-4a9242c79595"
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
ClusterManagers = "34f1f09b-3a8b-5176-ab39-66d58a4d544e"
Debugger = "31a5f54b-26ea-5ae9-a837-f05ce5417438"
Documenter = "e30172f5-a6a5-5a46-863b-614d45cd2de4"
HypothesisTests = "09f84164-cd44-5f33-b23f-e6b0d136a0d5"
Literate = "98b081ad-f1c9-55d3-8b20-4c87d4299306"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
Profile = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"
Revise = "295af30f-e4ad-537b-8983-00126c2a3abe"
StableRNGs = "860ef19b-820b-49d6-a774-d7a799459cd3"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[targets]
docs = ["Documenter", "Literate"]
test = ["Test", "StableRNGs", "HypothesisTests", "StatsBase", "Aqua"]
```

**Key changes:**
1. Remove BenchmarkTools, Debugger, Revise, Plots, ClusterManagers from `[deps]`
2. Add them to `[extras]` (BenchmarkTools and Debugger already there -- just remove from deps)
3. Add StableRNGs, HypothesisTests, StatsBase, Aqua to `[extras]`
4. Update `[targets]` test to include all test-only deps
5. Add `[compat]` entries for new dependencies

### Compilation Smoke Test (test_compilation.jl)

```julia
@testset "Compilation and Loading" begin
    @testset "Module loads without errors" begin
        # If we got here, `using QuantumFurnace` succeeded (loaded in runtests.jl)
        @test true
    end

    @testset "build_trajectoryframework compiles" begin
        config = make_thermalize_config(EnergyDomain(); with_coherent=true)
        precomputed_data = precompute_data(config.domain, config, TEST_HAM)
        scratch = KrausScratch(ComplexF64, DIM)

        fw = build_trajectoryframework(
            TEST_JUMPS, TEST_HAM, config, precomputed_data, scratch, TEST_DELTA
        )
        @test fw isa TrajectoryFramework
        @test fw.delta == TEST_DELTA
    end

    @testset "build_trajectoryframework without coherent" begin
        config = make_thermalize_config(EnergyDomain(); with_coherent=false)
        precomputed_data = precompute_data(config.domain, config, TEST_HAM)
        scratch = KrausScratch(ComplexF64, DIM)

        fw = build_trajectoryframework(
            TEST_JUMPS, TEST_HAM, config, precomputed_data, scratch, TEST_DELTA
        )
        @test fw isa TrajectoryFramework
        @test fw.B === nothing
        @test fw.U_B === nothing
    end

    @testset "Fixtures available" begin
        @test size(TEST_HAM.data) == (DIM, DIM)
        @test length(TEST_JUMPS) == 3 * NUM_QUBITS  # 12 jumps
        @test size(TEST_GIBBS) == (DIM, DIM)
        @test isapprox(tr(TEST_GIBBS), 1.0; atol=TOL_EXACT)
    end
end
```

## Detailed Bug Analysis

### Bug 1: Undefined `trotter` in build_trajectoryframework (CRITICAL -- blocks compilation)
**Location:** `src/trajectories.jl` line 53
**Root cause:** `trotter` is not a parameter or local variable of `build_trajectoryframework`. It was likely a kwarg in an older version that got refactored out.
**Additional context:** The function already receives `ham_or_trott` which can be `HamHam` or `TrottTrott`. The called function `precompute_coherent_total_B` no longer accepts a `trotter` kwarg either -- it dispatches internally based on `config.domain`.
**Fix:** Remove `; trotter=trotter` from the call. The function works correctly without it.
**Confidence:** HIGH -- verified by reading both function signatures.

### Bug 2: Uninitialized `B_total` (CRITICAL -- crashes when with_coherent=false)
**Location:** `src/trajectories.jl` line 96
**Root cause:** When `config.with_coherent == false`, the else branch sets `U_B = nothing` but never defines `B_total`. The constructor on line 96 references `B_total` unconditionally.
**Fix:** Add `B_total = nothing` in the else branch.
**Confidence:** HIGH -- verified by reading the code paths.

### Bug 3: Dangling `where T` (NON-CRITICAL -- compiles but misleading)
**Location:** `src/trajectories.jl` line 48
**Root cause:** Leftover from a parametric version of the function.
**Fix:** Remove `where T`.
**Confidence:** HIGH.

### Bug 4: Missing comma in export block (NON-CRITICAL -- precompute_data not exported)
**Location:** `src/QuantumFurnace.jl` lines 37-38
**Root cause:** Missing comma after `precompute_R` causes line 38 to be a standalone expression.
**Fix:** Add comma, remove `verify_completeness` (dead code).
**Confidence:** HIGH -- verified that `verify_completeness` is commented out.

### Bug 5: Same trotter kwarg issue in furnace.jl (NON-CRITICAL -- runtime error, not compile)
**Location:** `src/furnace.jl` line 66
**Root cause:** Same stale kwarg as Bug 1. Here `trotter` IS defined (as a function parameter), so it compiles. But `precompute_coherent_total_B` doesn't accept the kwarg, so it will error at runtime.
**Fix:** Remove `; trotter=trotter` from the call.
**Confidence:** HIGH.

## Discretion Recommendations

### TEST_DELTA Value
**Recommendation:** `TEST_DELTA = 0.01`
**Reasoning:** The existing `trajectory_test.jl` uses `delta = 0.1`. For test sensitivity, a smaller delta is better because it reduces unraveling error (`TOL_DELTA` scales linearly). `0.01` gives `TOL_DELTA = 0.05` (with C=5), which separates well from `TOL_QUADRATURE = 1e-6`. Going smaller (e.g., 0.001) would make tests slower without adding much diagnostic value for Phase 1.

### C Constant for TOL_DELTA
**Recommendation:** `C = 5.0` initially, expect to calibrate empirically in Phase 4
**Reasoning:** The unraveling error scales as O(delta) per the paper's Theorem III.1. The constant C absorbs dimension-dependent and transition-function-dependent factors. For a 16-dimensional system, C=5.0 is a reasonable starting point. Phase 4 trajectory-vs-DM tests will reveal if this needs adjustment.

### Precomputation Depth in make_test_system
**Recommendation:** Precompute through `hamiltonian`, `jumps`, and `gibbs` only. Do NOT precompute domain configs or `precomputed_data` -- these should be constructed per-test because different tests need different domain/coherent settings.
**Reasoning:** Config creation is cheap (just struct construction). `precomputed_data` depends on the domain choice. Keeping these out of the fixture makes tests more explicit about what they're testing.

### Grid Parameters
**Recommendation:** Match `trajectory_test.jl` values exactly:
- `num_energy_bits = 12`
- `w0 = 0.05`
- `t0 = 2pi / (2^12 * 0.05)` (follows from the Fourier relation)
- `num_trotter_steps_per_t0 = 10`

**Reasoning:** These are the values used in the existing development test and in the paper simulations. Consistency avoids surprises.

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `[extras]`+`[targets]` for test deps | Pkg workspaces (`test/Project.toml`) | Julia 1.2+ (newer) | Legacy approach still fully supported; project currently uses it |
| `Revise` + `includet` for dev testing | Proper `Pkg.test()` with `runtests.jl` | Convention since Julia 1.0 | Both coexist; dev scripts stay for interactive use |

**Note:** The project uses the legacy `[extras]`+`[targets]` approach. This is fine and fully supported. No need to migrate to workspace-based test deps.

## Open Questions

1. **Which other deps should move out of [deps]?**
   - What we know: BenchmarkTools, Debugger, Revise, Plots are clearly dev-only. ClusterManagers is for multi-node execution.
   - What's unclear: Whether `Distributed` and `SharedArrays` are used in core package functionality or just for parallel execution that users may need.
   - Recommendation: Keep `Distributed` and `SharedArrays` in `[deps]` for now. They are stdlib modules with zero install cost and ARE used in production code (`jump_workers.jl`, `nufft.jl`).

2. **Aqua UUID**
   - What we know: Aqua.jl needs to be added to `[extras]` with its UUID.
   - Recommendation: The UUID is `4c88cf16-eb10-579e-8560-4a9242c79595`. Verify at install time.

3. **StableRNGs, HypothesisTests, StatsBase UUIDs**
   - These are needed for `[extras]`. UUIDs: StableRNGs = `860ef19b-820b-49d6-a774-d7a799459cd3`, HypothesisTests = `09f84164-cd44-5f33-b23f-e6b0d136a0d5`, StatsBase = `2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91`.
   - Recommendation: Verify UUIDs using `] add` in a temporary environment.

## Sources

### Primary (HIGH confidence)
- Direct source code analysis of `/Users/bence/code/QuantumFurnace.jl/src/trajectories.jl` -- bug identification
- Direct source code analysis of `/Users/bence/code/QuantumFurnace.jl/src/coherent.jl` -- function signature verification
- Direct source code analysis of `/Users/bence/code/QuantumFurnace.jl/src/QuantumFurnace.jl` -- export block bug
- Direct source code analysis of `/Users/bence/code/QuantumFurnace.jl/src/furnace.jl` -- secondary trotter kwarg bug
- Direct source code analysis of `/Users/bence/code/QuantumFurnace.jl/Project.toml` -- dep duplication analysis
- [Julia Pkg.jl docs - Creating Packages](https://pkgdocs.julialang.org/v1/creating-packages/) -- `[extras]`+`[targets]` conventions
- [Julia Pkg.jl docs - TOML files](https://pkgdocs.julialang.org/v1/toml-files/) -- Project.toml structure

### Secondary (MEDIUM confidence)
- [StableRNGs.jl](https://github.com/JuliaRandom/StableRNGs.jl) -- RNG with stable streams for testing
- [Aqua.jl](https://github.com/JuliaTesting/Aqua.jl) -- Auto quality assurance
- [HypothesisTests.jl](https://github.com/JuliaStats/HypothesisTests.jl) -- Statistical tests
- [StatsBase.jl](https://github.com/JuliaStats/StatsBase.jl) -- Basic statistics

### Tertiary (LOW confidence)
- Package UUIDs for StableRNGs, HypothesisTests, StatsBase, Aqua -- should be verified at install time

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- Julia testing patterns are well-established and verified against official docs
- Architecture: HIGH -- based on Julia conventions and direct codebase analysis
- Bug analysis: HIGH -- every bug identified by reading source code directly
- Pitfalls: HIGH -- based on observed patterns in the actual codebase
- Discretion recommendations: MEDIUM -- TEST_DELTA and C values are informed estimates, not empirically validated yet

**Research date:** 2026-02-13
**Valid until:** No expiration -- these findings are based on the actual codebase state, not external API versions
