# Stack Research: v1.0 Trajectories Milestone

**Domain:** Trajectory simulation validation, statistical comparison, and correctness testing for quantum Lindbladian evolution
**Researched:** 2026-02-13
**Confidence:** HIGH (most recommendations use Julia stdlib or verified JuliaStats packages; one MEDIUM item noted)

## Scope

This stack research covers ONLY the additions needed for the v1.0 Trajectories milestone:
1. Fixing two-stage jump sampling
2. Validating trajectory-averaged rho against DM simulation
3. Building a comprehensive correctness test suite
4. Statistical comparison tools for stochastic results

It does NOT re-research the existing stack (Arpack, FINUFFT, LinearAlgebra, etc.) or future features (Qiskit, multi-threading, >12 qubits). See the project-level STACK.md in git history for those.

## Current Stack (Relevant Subset)

These existing dependencies are directly used by trajectory code and need no changes:

| Dependency | Role in Trajectories | Status |
|------------|---------------------|--------|
| LinearAlgebra (stdlib) | BLAS.gerc! for rho accumulation, mul!, dot, cholesky, Hermitian | Keep as-is |
| Random (stdlib) | rand() for branch selection in step_along_trajectory! | Keep as-is |
| Test (stdlib) | @test, @testset macros | Keep -- extend usage |
| Statistics (stdlib) | mean, var, std for trajectory averaging | Keep -- already available but underutilized |

## New Dependencies to Add

### For Testing Infrastructure

| Technology | Version | Purpose | Why Recommended |
|------------|---------|---------|-----------------|
| StableRNGs.jl | >= 1.0.1 | Reproducible RNG for deterministic trajectory tests | Julia's default RNG streams are NOT stable across versions. For regression tests that depend on specific random sequences (e.g., verifying a particular jump selection path), you need `StableRNG(seed)` to get identical results on Julia 1.12, 1.13, etc. Passed Big Crush statistical tests. Used by QuantEcon.jl, Flux.jl, and many JuliaStats packages for exactly this purpose. Add to `[extras]` test target. |
| HypothesisTests.jl | >= 0.11 | Statistical significance tests for trajectory vs DM comparison | Provides `OneSampleTTest`, `UnequalVarianceTTest`, and `ChisqTest`. When comparing trajectory-averaged expectation values against DM-computed values, you need a principled way to assess whether discrepancies are within statistical noise or indicate a bug. A one-sample t-test of trajectory means against the DM value, with a chosen significance level, is the correct tool. Part of JuliaStats org, well-maintained. Add to `[extras]` test target. |
| Aqua.jl | >= 0.8 | Automated package quality checks | Catches method ambiguities (critical with QuantumFurnace's heavy multiple dispatch -- 4 domain types x 2 config types x multiple method signatures), stale dependencies, missing compat entries, undefined exports. A single `Aqua.test_all(QuantumFurnace)` call in CI prevents entire classes of subtle bugs. Standard practice for Julia packages. Add to `[extras]` test target. |

**Confidence: HIGH** -- StableRNGs verified via [GitHub repo](https://github.com/JuliaRandom/StableRNGs.jl) (v1.0.1, Jan 2024, stable API). HypothesisTests verified via [official docs](https://juliastats.org/HypothesisTests.jl/stable/parametric/) (OneSampleTTest, ChisqTest confirmed). Aqua verified via [official docs](https://juliatesting.github.io/Aqua.jl/dev/).

### For Statistical Analysis in Tests

| Technology | Version | Purpose | Why Recommended |
|------------|---------|---------|-----------------|
| Statistics (stdlib) | ships with Julia | mean(), var(), std() for trajectory ensemble statistics | Already available -- no install needed. Currently underutilized. Use `mean(samples)` and `std(samples) / sqrt(length(samples))` for standard error of trajectory-averaged observables. |
| StatsBase.jl | >= 0.34 | sem() (standard error of mean), weighted statistics | Provides `sem(x)` = `std(x)/sqrt(length(x))` as a convenience, plus `L2dist`, `kldivergence` for comparing probability distributions element-wise. Useful for comparing diagonal elements of trajectory-averaged rho against Gibbs state. Lightweight (few deps). Add to `[extras]` test target. |

**Confidence: HIGH** -- Statistics stdlib verified via [Julia 1.12 docs](https://docs.julialang.org/en/v1/stdlib/Statistics/). StatsBase verified via [official docs](https://juliastats.org/StatsBase.jl/stable/).

## What NOT to Add for This Milestone

| Avoid | Why | What to Do Instead |
|-------|-----|-------------------|
| OrdinaryDiffEqTsit5 / DifferentialEquations.jl | Overkill for this milestone. You do NOT need an independent ODE reference solver -- you already HAVE the DM simulation (`run_thermalization`) as the reference. Trajectory vs DM comparison uses the DM stepper output directly. ODE solvers are for continuous-time master equation integration, which is a different problem. Adds ~30+ transitive deps for zero value here. | Compare trajectory rho_mean against `run_thermalization` output directly. Both use the same delta-stepping Kraus map -- the DM version applies it to the full density matrix, trajectories sample from it. Their agreement validates the sampling, not the Kraus map itself. |
| ExponentialUtilities.jl | Not needed for trajectory validation. The trajectory code uses discrete Kraus-map steps (K0, jump operators, residual), not matrix exponentials. DM validation similarly uses discrete steps. expv/expm is irrelevant. | Keep using the existing discrete delta-step framework. |
| KrylovKit.jl | Arpack replacement is a separate concern. Eigendecomposition is done once during setup (finalize_hamham), not during trajectory evolution or testing. Swapping eigensolvers during the validation milestone adds risk with zero benefit. | Defer to a future cleanup milestone. |
| Distributions.jl | Full distribution library is heavyweight (~15 deps). For this milestone you need chi-squared p-values and normal quantiles -- HypothesisTests.jl already provides the tests you need without requiring you to manually construct distribution objects. | Use HypothesisTests.jl which handles the underlying distributions internally. |
| SafeTestsets.jl | Adds module isolation per testset. QuantumFurnace tests are not large enough to have state leakage issues. The Julia 1.12 `@testset` with `rng` keyword already handles the main isolation concern (RNG state). Module isolation adds complexity without solving a real problem at this scale. | Use plain `@testset` with explicit RNG seeding. Group tests into separate files included from `runtests.jl`. |
| MonteCarloMeasurements.jl | Propagates uncertainty through computations. Interesting but wrong abstraction -- you need to compare trajectory statistics against a known reference, not propagate uncertainty through a computation pipeline. | Direct statistical comparison via t-tests and confidence intervals. |
| QuantumOptics.jl / QuantumToolbox.jl | Do NOT import these as validation references. They have their own type systems and would force conversion layers. Your DM simulation IS the reference. Cross-validation against a third-party framework is a separate optional future task. | Use your own `run_thermalization` as the ground truth for trajectory validation. |

## Recommended Stack for This Milestone

### Production Dependencies (src/)

**No new production dependencies needed.** The trajectory fix is a code correction in `step_along_trajectory!` and `build_trajectoryframework`, not a dependency addition. All required numerical tools (LinearAlgebra, Random, BLAS) are already present.

### Test Dependencies (test/)

Add to `Project.toml` under `[extras]` and `[targets]`:

| Package | Purpose | Section |
|---------|---------|---------|
| StableRNGs | Reproducible trajectory seeds | [extras] + test target |
| HypothesisTests | Statistical correctness assertions | [extras] + test target |
| StatsBase | Standard error, distribution distances | [extras] + test target |
| Aqua | Package quality | [extras] + test target |

These are TEST-ONLY dependencies. They do not increase the dependency burden for users of the package.

## Installation

```julia
# From the package directory, add test-only dependencies:
using Pkg
Pkg.activate(".")

# These go in [extras] section of Project.toml, not [deps]
# Add UUIDs manually or via:
# Pkg.add(["StableRNGs", "HypothesisTests", "StatsBase", "Aqua"])
# then move entries from [deps] to [extras] and add to test target
```

Concrete `Project.toml` changes:

```toml
[extras]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
Debugger = "31a5f54b-26ea-5ae9-a837-f05ce5417438"
Documenter = "e30172f5-a6a5-5a46-863b-614d45cd2de4"
Literate = "98b081ad-f1c9-55d3-8b20-4c87d4299306"
Profile = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
StableRNGs = "860ef19b-820b-49d6-a774-d7a799459cd3"
HypothesisTests = "09f84164-cd44-5f33-b23f-e6b0d136a0d5"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
Aqua = "4c88cf16-eb10-579e-8560-4a9242c79595"

[targets]
docs = ["Documenter", "Literate"]
test = ["Test", "StableRNGs", "HypothesisTests", "StatsBase", "Aqua"]
```

## How Each Addition Integrates with Existing Code

### StableRNGs Integration

The existing `step_along_trajectory!` uses `rand()` (global RNG) for branch selection. For deterministic tests:

```julia
using StableRNGs
using Random

# Seed global RNG with StableRNG for reproducible trajectory
rng = StableRNG(42)
copy!(Random.default_rng(), rng)  # or use Julia 1.12 @testset rng= feature

# Now run_trajectories will produce identical results every time
result = run_trajectories(jumps, config, psi0, hamiltonian; ntraj=1)
```

Julia 1.12's `@testset` supports `rng` keyword directly:
```julia
@testset "trajectory regression" rng=StableRNG(12345) begin
    # Global RNG is seeded deterministically within this testset
    result = run_trajectories(...)
    @test isapprox(result.rho_mean, expected_rho, atol=1e-10)
end
```

**No changes to production code needed** -- the tests control the global RNG externally.

### HypothesisTests Integration

For trajectory vs DM statistical comparison:

```julia
using HypothesisTests

# Run many trajectories, collect per-trajectory observable values
observable_values = Float64[]
for _ in 1:ntraj
    psi = copy(psi0)
    _evolve_along_trajectory!(psi, fw, total_time)
    push!(observable_values, real(dot(psi, observable * psi)))
end

# DM gives the exact expected value
dm_expected = real(tr(observable * rho_dm))

# Statistical test: is the trajectory mean consistent with DM value?
test_result = OneSampleTTest(observable_values, dm_expected)
@test pvalue(test_result) > 0.001  # Fail only if extremely unlikely
```

For trace distance comparison with error bars:
```julia
# Collect per-trajectory rho contributions
rho_samples = [outer_product(psi_i) for psi_i in trajectory_finals]
trace_distances = [trace_distance_nh(rho_i, rho_dm) for rho_i in rho_samples]

# The MEAN trace distance should decrease with more trajectories
# Test that averaged rho is close to DM rho
@test trace_distance_nh(mean(rho_samples), rho_dm) < tolerance
```

### StatsBase Integration

```julia
using StatsBase

# Standard error of trajectory-averaged observable
obs_values = [measure(psi_i, O) for psi_i in trajectories]
stderr = sem(obs_values)  # std(obs_values) / sqrt(length(obs_values))

# Convergence rate validation: stderr should scale as 1/sqrt(N)
# Double trajectories -> stderr should halve (approximately)
```

### Aqua Integration

Single test file:
```julia
# test/aqua_test.jl
using Aqua
Aqua.test_all(QuantumFurnace;
    ambiguities=true,        # Critical with 4 domain types x multiple dispatch
    stale_deps=true,         # Catch deps in [deps] no longer used
    deps_compat=true,        # All deps need compat entries
    project_extras=true,     # Test extras must be consistent
    piracies=true,           # No type piracy
)
```

## Statistical Validation Strategy (Informing Test Design)

The stack choices above enable this validation hierarchy:

1. **Deterministic regression tests** (StableRNGs): Fixed seed, verify exact trajectory output matches known-good values. Catches any code change that alters the sampling logic.

2. **Single-trajectory property tests** (Test stdlib): For ANY trajectory, verify:
   - State is normalized after every step: `|dot(psi, psi) - 1| < eps`
   - Probabilities sum to ~1: `p_nojump + p_res + p_jump_total approx 1`
   - No NaN/Inf in state vector

3. **Statistical ensemble tests** (HypothesisTests + StatsBase): Run N trajectories, compare:
   - Trajectory-averaged rho vs DM rho (trace distance < f(N, delta))
   - Per-observable means vs DM values (t-test, p > threshold)
   - Convergence rate: error ~ 1/sqrt(N) verified by doubling N

4. **Cross-domain consistency tests** (Test stdlib): For small systems:
   - Energy domain trajectory avg should match Energy domain DM
   - Time domain trajectory avg should match Time domain DM
   - Error nesting: |Bohr - Energy| < |Bohr - Time| < |Bohr - Trotter|

## Alternatives Considered

| Recommended | Alternative | When to Use Alternative |
|-------------|-------------|-------------------------|
| StableRNGs for seed control | Julia 1.12 `@testset rng=MersenneTwister(seed)` | If you want zero extra deps and accept that streams may change between Julia patch versions. StableRNGs is safer for long-term regression tests. |
| HypothesisTests for statistical comparison | Manual z-test computation (`(mean - expected) / (std/sqrt(n))` compared to quantile) | If you want zero deps and are comfortable implementing the test statistic correctly. HypothesisTests handles edge cases (small n, unequal variance) that manual implementations often get wrong. |
| StatsBase for sem() | `std(x) / sqrt(length(x))` inline | Perfectly fine for this one function. StatsBase adds value if you also use `kldivergence` or `L2dist` for distribution comparison. If you only need sem, skip StatsBase. |
| Aqua for quality checks | Manual ambiguity checks | Never. Aqua is trivial to add and catches bugs you would never find manually. |
| No ODE solver | OrdinaryDiffEqTsit5 for independent reference | Only if you distrust the DM simulation AND the trajectory simulation simultaneously. In that case, solving drho/dt = L(rho) with an ODE solver gives a third independent reference. Not needed for this milestone -- the DM stepper is already validated. |

## Version Compatibility

| Package | Compatible With | Notes |
|---------|-----------------|-------|
| StableRNGs >= 1.0.1 | Julia >= 1.0 | Pure Julia, zero deps beyond Random stdlib. Trivially compatible. |
| HypothesisTests >= 0.11 | Julia >= 1.6 | Depends on Distributions.jl, StatsBase.jl, Combinatorics.jl. These are test-only deps so the transitive tree does not affect users. |
| StatsBase >= 0.34 | Julia >= 1.6 | Lightweight. Few transitive deps (SortingAlgorithms, Missings). |
| Aqua >= 0.8 | Julia >= 1.6 | Pure Julia, minimal deps. |
| All of the above | Julia 1.12.x | No known compatibility issues. All packages are actively maintained in the JuliaStats/JuliaTesting orgs. |

## Cleanup Needed in Current Project.toml

While adding test deps, also clean up the existing `[deps]` section. These packages are in `[deps]` but should NOT be (they are dev/doc tools, not runtime dependencies):

| Package | Current | Should Be | Impact |
|---------|---------|-----------|--------|
| Revise | [deps] | Remove entirely (load manually in REPL) | Users won't be forced to install Revise |
| Debugger | [deps] | Remove entirely (load manually) | Same |
| BenchmarkTools | [deps] AND [extras] | [extras] only | Remove from [deps] |
| Plots | [deps] | Remove from [deps]; use in scripts only | Heavy dep (50+ transitive). Not needed at runtime. |
| Pkg | [deps] | Remove | stdlib, auto-available |
| DocumenterTools | [deps] | [extras] with docs target | Only needed for doc generation |
| ClusterManagers | [deps] | Keep for now | Used by Distributed workflows |

This cleanup is not blocking but should happen during this milestone to get the package registry-ready.

## Sources

- [StableRNGs.jl GitHub](https://github.com/JuliaRandom/StableRNGs.jl) -- v1.0.1, LehmerRNG implementation, passed Big Crush. HIGH confidence.
- [HypothesisTests.jl parametric tests docs](https://juliastats.org/HypothesisTests.jl/stable/parametric/) -- OneSampleTTest, ChisqTest signatures verified. HIGH confidence.
- [StatsBase.jl scalar statistics docs](https://juliastats.org/StatsBase.jl/v0.20/scalarstats.html) -- sem() function verified. HIGH confidence.
- [Aqua.jl official docs](https://juliatesting.github.io/Aqua.jl/dev/) -- test_all() API verified. HIGH confidence.
- [Julia 1.12 Test stdlib docs](https://docs.julialang.org/en/v1/stdlib/Test/) -- @testset rng= keyword confirmed for Julia >= 1.12. HIGH confidence.
- [Julia Statistics stdlib docs](https://docs.julialang.org/en/v1/stdlib/Statistics/) -- mean, var, std confirmed. HIGH confidence.
- [QuantumOptics.jl MCWF docs](https://docs.qojulia.org/timeevolution/mcwf/) -- trajectory vs master equation validation pattern confirmed. HIGH confidence.
- [QuantumToolbox.jl Monte Carlo solver](https://qutip.org/QuantumToolbox.jl/stable/users_guide/time_evolution/mcsolve) -- mcsolve validation approach. HIGH confidence.
- [Monte Carlo wave-function convergence study (Daley, 2014)](https://arxiv.org/pdf/1803.08589) -- convergence rate analysis for MCWF. MEDIUM confidence (training data, not directly verified).

---
*Stack research for: QuantumFurnace.jl v1.0 Trajectories milestone*
*Researched: 2026-02-13*
