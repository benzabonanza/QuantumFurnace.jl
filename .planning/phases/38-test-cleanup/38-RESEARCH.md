# Phase 38: Test Cleanup - Research

**Researched:** 2026-02-28
**Domain:** Julia test infrastructure consolidation, threshold validation, informative test output
**Confidence:** HIGH

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

#### Helper consolidation
- Single `make_config(sim, domain; num_qubits=4, construction=KMS(), ...)` factory replacing all 8 current functions
- Single `make_test_system(; num_qubits=4, trotter=nothing)` replacing both `make_test_system()` and `make_small_test_system()`
- Keep precomputed globals (TEST_HAM, TEST_JUMPS, etc.) for both sizes -- rename SMALL_* to N3_* for clarity
- Add `const ALL_DOMAINS = [EnergyDomain(), TimeDomain(), TrotterDomain(), BohrDomain()]` in test_helpers.jl

#### @info output design
- Pattern: `@info "label" value=computed_value threshold=threshold_used` -- placed AFTER the assertion
- Only for numerical comparisons (trace distance, eigenvalue, allocation bytes, convergence rate)
- Skip @info for structural/type checks (round-trip serialization, field existence, etc.)
- Every numerical @test shows: label, computed value, and threshold it was compared against

#### Threshold review
- Audit ALL 239 threshold comparisons across all 22 test files
- Theory-based tightening first (O(delta), O(delta^2), machine epsilon, known error scaling)
- Empirical tightening acceptable but must account for scaling to larger systems (double qubit count) -- small-system empirical values don't reliably bound larger systems
- Existing tiers (TOL_EXACT=1e-12, TOL_QUADRATURE=1e-6, TOL_DELTA(delta)=5*delta) are good -- extend with new named constants where needed, but think through each case
- Rationale documented as inline comments next to each threshold check

#### Old/staging test handling
- DELETE test/old_tests/ entirely (ham_test, kossakowski_test, trott_test, B_test, trajectory_test, time_tests)
- KEEP test/staging/ (test_fitting.jl, test_gap_estimation.jl) -- matches dormant src/staging/ code
- INTEGRATE test/trajectory_validation/ into runtests.jl -- behind a `FULL` test flag if tests are slow (few minutes)
- REVIEW test/reference/generate_references.jl -- verify it still matches current types/API after v2.0 restructure

### Claude's Discretion
- Exact naming of new tolerance constants beyond the existing tiers
- Whether trajectory validation tests need the FULL flag (depends on measured runtime)
- How to implement the FULL test flag (ENV variable check)
- Grouping/ordering of @info output within test files

### Deferred Ideas (OUT OF SCOPE)
None -- discussion stayed within phase scope
</user_constraints>

## Summary

This phase consolidates test infrastructure for a Julia quantum physics package (QuantumFurnace.jl). The test suite has 664 `@test` assertions across 22 active test files (17 in runtests.jl + 2 in trajectory_validation + 2 in staging + test_aqua.jl), organized around a `test_helpers.jl` file that provides precomputed test systems and 8 config factory functions. The work involves three interrelated tasks: (1) consolidating the 8 config factories into a single parametrized `make_config` and the 2 system factories into a single `make_test_system`, (2) adding `@info` output after every numerical `@test` so test output is self-documenting, and (3) auditing all ~276 threshold comparisons to ensure they have theory-based or empirically-justified values with inline rationale.

The codebase also has 7 dead test files in `test/old_tests/` to delete, 2 standalone trajectory validation scripts to integrate into the main test suite (gated behind an environment variable), and a reference generator to verify still matches the current API. The existing `QUANTUMFURNACE_FULL_TESTS` environment variable pattern is already in use (in run_convergence_tests.jl and test_krylov_crossvalidation.jl) and should be reused for gating slow tests.

**Primary recommendation:** Work file-by-file through the test suite, consolidating factory calls to the new unified `make_config`/`make_test_system`, adding `@info` after numerical assertions, and tightening thresholds with inline rationale comments -- testing after each file to ensure no regressions.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Test (stdlib) | Julia 1.12 | @test, @testset macros | Julia standard test framework |
| LinearAlgebra (stdlib) | Julia 1.12 | Matrix operations, norms, eigen | Standard for numerical tests |
| Random (stdlib) | Julia 1.12 | Xoshiro RNG for deterministic tests | Standard Julia RNG |
| BSON | 0.3.9 | Serialization round-trip tests | Already used for reference data |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| Aqua.jl | (test dep) | Package quality tests | test_aqua.jl only |
| StableRNGs | (test dep) | Stable RNG for staging fitting tests | staging tests only |
| Printf (stdlib) | Julia 1.12 | Formatted diagnostic output | Cross-validation tables |

No new dependencies needed. This phase only reorganizes existing test code.

## Architecture Patterns

### Current Test File Layout
```
test/
  runtests.jl               # Main entry: includes test_helpers.jl then 17 test files
  test_helpers.jl            # Fixtures, constants, 8 factory functions, precomputed globals
  test_allocation.jl         # 7 allocation regression tests
  test_aqua.jl               # Package quality (Aqua.jl)
  test_compilation.jl        # Loading and fixture smoke tests
  test_convergence.jl        # 764 lines, largest file, convergence tracking
  test_cptp.jl               # CPTP completeness across 3 domains
  test_diagnostics.jl        # 525 lines, exact diagnostics (DIAG-01..06)
  test_dm_detailed_balance.jl # DMTST-01, DMTST-02
  test_dm_scaling.jl         # DMTST-03..06b, Euler/OFT/NUFFT scaling
  test_gns_trajectory.jl     # GNS-01, GNS-02
  test_krylov_crossvalidation.jl # XVAL-01..04 (+ n=6 env-gated)
  test_krylov_eigsolve.jl    # Krylov eigsolve accuracy
  test_krylov_matvec.jl      # Krylov matvec round-trip (23 testsets)
  test_observable_trajectories.jl # Observable-only trajectory runner
  test_regression.jl         # TINF-02 frozen BSON regression
  test_results.jl            # Result type serialization round-trips
  test_threading.jl          # Deterministic threading, speedup
  test_trajectory_fixes.jl   # TFIX-02..05 bug fix guards
  test_workspace_independence.jl # Workspace isolation
  old_tests/                 # 7 dead files (TO DELETE)
  staging/                   # 2 files matching src/staging/ (KEEP)
  trajectory_validation/     # 2 standalone scripts (TO INTEGRATE)
  reference/                 # BSON reference data + generator script
```

### Pattern 1: Unified Config Factory
**What:** Replace 8 separate factory functions with a single parametrized factory
**When to use:** Every test file that creates Config objects

The 8 existing factories are:
1. `make_liouv_config(domain; construction=KMS())` -- Lindbladian, 4-qubit
2. `make_liouv_config_gns(domain)` -- Lindbladian, 4-qubit, GNS hardcoded
3. `make_thermalize_config(domain; construction=KMS(), delta, mixing_time)` -- Thermalize, 4-qubit
4. `make_small_liouv_config(domain; construction=GNS())` -- Lindbladian, 3-qubit
5. `make_small_liouv_config_gns(domain)` -- Lindbladian, 3-qubit, GNS hardcoded
6. `make_small_thermalize_config(domain; construction=GNS(), delta, mixing_time)` -- Thermalize, 3-qubit
7. `make_small_thermalize_config_gns(domain; delta, mixing_time)` -- Thermalize, 3-qubit, GNS hardcoded
8. `make_small_liouv_config_gns(domain)` -- duplicate of #5

**Unified replacement:**
```julia
function make_config(sim, domain;
    num_qubits=NUM_QUBITS,
    construction=KMS(),
    delta=TEST_DELTA,
    mixing_time=1.0,
)
    Config(
        sim = sim,
        domain = domain,
        construction = construction,
        num_qubits = num_qubits,
        with_linear_combination = true,
        beta = BETA,
        sigma = SIGMA,
        a = BETA / 30.0,
        b = 0.4,
        num_energy_bits = NUM_ENERGY_BITS,
        w0 = W0,
        t0 = T0,
        num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
        # Only pass these for Thermalize:
        (sim isa Thermalize ? (mixing_time = mixing_time, delta = delta) : ())...,
    )
end
```

Note: The `Config` constructor likely uses keyword arguments that are ignored for `Lindbladian` sim type. Need to verify if `Config(sim=Lindbladian(), ..., mixing_time=1.0, delta=0.01)` errors or silently ignores extra kwargs. If it errors, use conditional splatting as shown above.

**Migration examples:**
```julia
# Old: make_liouv_config(EnergyDomain())
# New: make_config(Lindbladian(), EnergyDomain())

# Old: make_small_liouv_config_gns(TrotterDomain())
# New: make_config(Lindbladian(), TrotterDomain(); num_qubits=3, construction=GNS())

# Old: make_thermalize_config(EnergyDomain(); delta=0.01, mixing_time=60.0)
# New: make_config(Thermalize(), EnergyDomain(); delta=0.01, mixing_time=60.0)

# Old: make_small_thermalize_config_gns(TrotterDomain(); delta=0.01, mixing_time=100.0)
# New: make_config(Thermalize(), TrotterDomain(); num_qubits=3, construction=GNS(), delta=0.01, mixing_time=100.0)
```

### Pattern 2: Unified Test System Factory
**What:** Replace `make_test_system()` and `make_small_test_system()` with a single parametrized function
**When to use:** test_helpers.jl globals

```julia
function make_test_system(; num_qubits=NUM_QUBITS, trotter=nothing)
    source_root = dirname(@__DIR__)
    ham_path = joinpath(source_root, "hamiltonians", "heis_disordered_periodic_n$(num_qubits).bson")
    hamiltonian = _load_test_hamiltonian(ham_path, BETA)

    jump_paulis = [[X], [Y], [Z]]
    jump_sites = 1:num_qubits
    num_of_jumps = length(jump_paulis) * length(jump_sites)
    jump_normalization = sqrt(num_of_jumps)

    basis_unitary = trotter !== nothing ? trotter.eigvecs : hamiltonian.eigvecs

    jumps = JumpOp[]
    for pauli in jump_paulis
        for site in jump_sites
            jump_op = Matrix(pad_term(pauli, num_qubits, site)) ./ jump_normalization
            jump_in_eigen = basis_unitary' * jump_op * basis_unitary
            orthogonal = (jump_op == transpose(jump_op))
            herm = (jump_op == jump_op')
            push!(jumps, JumpOp(jump_op, jump_in_eigen, orthogonal, herm))
        end
    end

    gibbs = hamiltonian.gibbs
    return (; hamiltonian, jumps, gibbs)
end
```

### Pattern 3: Global Renaming (SMALL_* -> N3_*)
**What:** Rename precomputed globals for 3-qubit system from `SMALL_*` to `N3_*`
**Current globals to rename:**
```julia
# Current -> New
SMALL_SYSTEM    -> N3_SYSTEM
SMALL_HAM       -> N3_HAM
SMALL_JUMPS     -> N3_JUMPS
SMALL_GIBBS     -> N3_GIBBS
SMALL_DIM       -> N3_DIM
SMALL_TROTTER   -> N3_TROTTER
SMALL_TROTTER_JUMPS -> N3_TROTTER_JUMPS
```

These are referenced across many files. Use search-replace carefully. Count of references per name:
- `SMALL_HAM`: ~40 references across 10 files
- `SMALL_JUMPS`: ~20 references across 8 files
- `SMALL_GIBBS`: ~10 references across 5 files
- `SMALL_DIM`: ~25 references across 8 files
- `SMALL_TROTTER`: ~15 references across 6 files
- `SMALL_TROTTER_JUMPS`: ~12 references across 5 files

### Pattern 4: @info After Numerical Assertions
**What:** Add `@info` with label, computed value, and threshold after every numerical `@test`
**When to use:** Numerical comparisons only (trace distance, eigenvalue, allocation, convergence rate)
**Skip for:** Type checks, field existence, structural assertions, bitwise equality

```julia
# CORRECT: @info AFTER the assertion, with named keyword arguments
dist = trace_distance_h(Hermitian(ss_dm), SMALL_GIBBS)
@test dist < 1e-10
@info "DMTST-01: Bohr fixed point -> Gibbs" trace_distance=dist threshold=1e-10

# CORRECT: isapprox with atol
@test isapprox(tr(fp.fixed_point), 1.0; atol=1e-12)
@info "Fixed point trace" value=real(tr(fp.fixed_point)) threshold_atol=1e-12

# CORRECT: allocation test
@test allocs == 0
@info "step_along_trajectory! allocations" allocs_bytes=allocs threshold=0

# SKIP @info for these:
@test result isa ConvergenceData           # structural
@test conv_data.converged == true           # boolean
@test haskey(d, :batch_sizes)               # structural
@test result.rho_mean == result2.rho_mean   # bitwise equality (not numerical threshold)
```

### Pattern 5: FULL Test Flag
**What:** Gate slow tests behind `QUANTUMFURNACE_FULL_TESTS` environment variable
**Existing pattern already in use:**
```julia
# Already used in run_convergence_tests.jl and test_krylov_crossvalidation.jl
if get(ENV, "QUANTUMFURNACE_FULL_TESTS", "false") == "true"
    @testset "Expensive test" begin
        # ...
    end
else
    @info "Skipping expensive test (set QUANTUMFURNACE_FULL_TESTS=true to run)"
end
```

### Pattern 6: ALL_DOMAINS Constant
**What:** Centralize domain list for parametrized tests
```julia
const ALL_DOMAINS = [EnergyDomain(), TimeDomain(), TrotterDomain(), BohrDomain()]
```

### Anti-Patterns to Avoid
- **Inline Config construction**: DMTST-01 in test_dm_detailed_balance.jl constructs Config(...) inline with 11 kwargs instead of using a factory. After consolidation, all Config creation must go through `make_config`.
- **Duplicate factory functions**: The `_gns` variants are just hardcoded `construction=GNS()`. The unified factory eliminates this.
- **println for diagnostics**: Some files use `println()` instead of `@info`. The decision is to use `@info` with keyword args for structured output.
- **Threshold without rationale**: Bare `atol=1e-10` without a comment explaining why 1e-10. Every threshold needs an inline comment.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Test output | Custom println formatting | `@info "label" key=val` | Julia's logging is structured, filterable, consistent |
| Config creation | Per-file Config(...) construction | `make_config(sim, domain; ...)` | Prevents parameter drift, centralizes defaults |
| System creation | Per-file Hamiltonian loading | `make_test_system(; num_qubits=N)` | Prevents duplication, ensures consistent normalization |
| Env-gated tests | Custom boolean parsing | `get(ENV, "QUANTUMFURNACE_FULL_TESTS", "false") == "true"` | Already established pattern in codebase |

**Key insight:** The test suite's maintainability problem is fundamentally about duplication -- 8 factory functions doing nearly the same thing, Config() constructed inline in one file, thresholds scattered without rationale. Consolidation into single-point-of-truth functions solves the root cause.

## Common Pitfalls

### Pitfall 1: Breaking Precomputed Globals During Rename
**What goes wrong:** Renaming SMALL_* to N3_* breaks references in files not updated simultaneously
**Why it happens:** Globals are defined in test_helpers.jl but referenced across 10+ files
**How to avoid:** Do the rename atomically -- search-replace across ALL files in a single commit. Use Julia's compilation errors (undefined variable) as a safety net.
**Warning signs:** `UndefVarError: SMALL_HAM not defined` during test run

### Pitfall 2: Config Factory Signature Mismatch
**What goes wrong:** The unified `make_config` passes Thermalize-only kwargs (delta, mixing_time) to Lindbladian configs, causing an error
**Why it happens:** Julia's `Config` constructor may not accept unused keyword arguments
**How to avoid:** Test the unified factory with both `Lindbladian()` and `Thermalize()` sim types before mass-migration. Use conditional splatting if needed.
**Warning signs:** `MethodError: no method matching Config(...)` with extra kwargs

### Pitfall 3: @info Placement Before vs After @test
**What goes wrong:** Placing `@info` before `@test` means the info is printed even when the test would pass, cluttering output. Worse: if `@info` computation fails, it masks the test result.
**Why it happens:** Natural reading order suggests "announce then check"
**How to avoid:** Always place `@info` AFTER the `@test` assertion (locked decision). The `@test` macro does not short-circuit -- both pass and fail cases continue execution.
**Warning signs:** @info appearing before @test in code review

### Pitfall 4: Threshold Tightening Without Scale Consideration
**What goes wrong:** Threshold looks great on 3-qubit but fails on 4-qubit or larger systems
**Why it happens:** Numerical error scales with system dimension (DIM=8 for 3q, DIM=16 for 4q, DIM=64 for 6q)
**How to avoid:** Theory-based tightening accounts for dimension scaling. If empirical, test at both 3q and 4q. User explicitly warned about this.
**Warning signs:** Tests pass locally but fail in CI with different system sizes

### Pitfall 5: Missing @info for Tests Inside Loops
**What goes wrong:** Tests like `for a in 1:ws.n_jumps ... @test isapprox(...) ... end` run multiple numerical comparisons but the @info would fire 12 times
**Why it happens:** Loop structure makes per-iteration @info noisy
**How to avoid:** For loop-based numerical tests, compute max/summary statistics and emit a single @info after the loop with `max_error=...` or similar.
**Warning signs:** Test output flooded with 12 identical @info lines per testset

### Pitfall 6: Changing println to @info in Cross-Validation Tables
**What goes wrong:** The formatted tables in test_krylov_crossvalidation.jl use Printf formatting with alignment. @info doesn't support Printf-style formatting natively.
**Why it happens:** @info uses Julia's logging backend which formats key=value pairs differently than Printf
**How to avoid:** Keep Printf-formatted tables as println (they serve a different purpose -- tabular comparison data). The user decision is for @info on "numerical @test" assertions specifically, not replacing all diagnostic output.
**Warning signs:** Mangled table output, lost alignment

### Pitfall 7: Trajectory Validation Integration Timing
**What goes wrong:** Integrating trajectory_validation/ tests into runtests.jl without gating makes the test suite take minutes instead of seconds
**Why it happens:** run_trajectory_validation.jl runs 50,000 trajectories; run_convergence_tests.jl runs 500,000
**How to avoid:** Gate behind QUANTUMFURNACE_FULL_TESTS (already used by run_convergence_tests.jl). Measure actual runtime before deciding.
**Warning signs:** CI timeouts, test suite taking >5 minutes

## Code Examples

### Example 1: Unified make_config Usage
```julia
# Source: Codebase analysis of test_helpers.jl factory functions

# Lindbladian + KMS (default) for 4-qubit (default)
config = make_config(Lindbladian(), EnergyDomain())

# Lindbladian + GNS for 3-qubit
config = make_config(Lindbladian(), TrotterDomain(); num_qubits=3, construction=GNS())

# Thermalize + KMS for 4-qubit with custom delta
config = make_config(Thermalize(), EnergyDomain(); delta=0.01, mixing_time=60.0)

# Thermalize + GNS for 3-qubit
config = make_config(Thermalize(), TimeDomain(); num_qubits=3, construction=GNS(), delta=0.01, mixing_time=0.5)
```

### Example 2: @info After Numerical @test
```julia
# Source: Locked decision from CONTEXT.md

# Trace distance comparison
dist = trace_distance_h(Hermitian(ss_dm), TEST_GIBBS)
@test dist < 1e-10  # KMS detailed balance: Gibbs is exact fixed point (machine precision)
@info "DMTST-01: fixed point -> Gibbs" trace_distance=dist threshold=1e-10

# isapprox with atol
err = norm(v_dense - vec(L_rho))
@test err < 1e-12  # Matvec matches dense to machine precision (no approximation error)
@info "Matvec round-trip" error=err threshold=1e-12

# Allocation test
allocs = @allocated step_along_trajectory!(psi, ws, rng2)
@test allocs == 0  # Hot path must be allocation-free for parallel scaling
@info "step_along_trajectory! allocations" allocs_bytes=allocs threshold=0

# Ratio/scaling test
@test 3.0 <= ratio <= 5.0  # O(delta^2): halving delta -> 4x error reduction
@info "Euler error ratio" ratio expected=4.0 lower_bound=3.0 upper_bound=5.0
```

### Example 3: Summary @info for Loop-Based Tests
```julia
# Source: Pattern derived from test_cptp.jl and test_krylov_matvec.jl

# Instead of @info per iteration:
max_err = 0.0
for a in 1:ws.n_jumps
    completeness = ws.K0s[a]' * ws.K0s[a] + ws.delta * ws.Rs[a] + ws.U_residuals[a]' * ws.U_residuals[a]
    err = norm(completeness - identity)
    max_err = max(max_err, err)
    @test isapprox(completeness, identity; atol=1e-10)  # CPTP completeness (algebraic identity)
end
@info "CPTP completeness (EnergyDomain)" n_jumps=ws.n_jumps max_error=max_err threshold_atol=1e-10
```

### Example 4: Trajectory Validation Integration
```julia
# Source: Existing QUANTUMFURNACE_FULL_TESTS pattern in codebase

# In runtests.jl, add at the end:
if get(ENV, "QUANTUMFURNACE_FULL_TESTS", "false") == "true"
    include("trajectory_validation/run_trajectory_validation.jl")
    include("trajectory_validation/run_convergence_tests.jl")
else
    @info "Skipping trajectory validation tests (set QUANTUMFURNACE_FULL_TESTS=true to run)"
end
```

### Example 5: Threshold Rationale Comments
```julia
# Source: Locked decision from CONTEXT.md -- inline rationale for every threshold

# Machine precision identity (algebraic, no approximation):
@test isapprox(completeness, identity; atol=1e-10)  # CPTP: K0'K0 + delta*R + U'U = I (algebraic identity, 1e-10 allows FP accumulation across DIM^2 entries)

# Quadrature approximation:
@test dist_bohr_time < TOL_QUADRATURE  # Time quadrature error: O(1/N_time_points) with N=4096, well below 1e-6

# Statistical noise:
@test isapprox(rho_traj, rho_dm; atol=0.05)  # 1000 trajectories, expected error ~ 1/sqrt(1000) ~ 0.03, threshold 0.05 gives 1.6x margin

# Convergence rate:
@test 3.0 <= ratio <= 5.0  # O(delta^2): exact ratio is 4.0, bounds allow for sub-leading terms at finite delta
```

## Detailed Inventory for Planning

### Current Config Factory Functions (8 total in test_helpers.jl)

| Function | sim | num_qubits | default construction | Call count |
|----------|-----|------------|---------------------|------------|
| `make_liouv_config` | Lindbladian | 4 | KMS | ~25 |
| `make_liouv_config_gns` | Lindbladian | 4 | GNS (hardcoded) | ~8 |
| `make_thermalize_config` | Thermalize | 4 | KMS | ~12 |
| `make_small_liouv_config` | Lindbladian | 3 | GNS | ~10 |
| `make_small_liouv_config_gns` | Lindbladian | 3 | GNS (hardcoded) | ~4 |
| `make_small_thermalize_config` | Thermalize | 3 | GNS | ~15 |
| `make_small_thermalize_config_gns` | Thermalize | 3 | GNS (hardcoded) | ~2 |
| (inline in test_dm_detailed_balance.jl) | Lindbladian | 3 | KMS | 1 |

**Additional factories in test_krylov_crossvalidation.jl (n=6 system):**
- `make_n6_liouv_config` -- Lindbladian, 6-qubit
- `make_n6_thermalize_config` -- Thermalize, 6-qubit
- `make_n6_test_system` -- 6-qubit system factory

These n=6 factories should also migrate to the unified `make_config` and `make_test_system`.

### Files Requiring @info Addition

| File | Numerical @tests | @info present | @info needed |
|------|-----------------|---------------|--------------|
| test_allocation.jl | 7 | 0 | ~7 (allocation checks) |
| test_compilation.jl | 3 | 0 | ~3 |
| test_convergence.jl | ~54 | 0 | ~30 (many are structural) |
| test_cptp.jl | 3 | 0 | ~3 (loop-summary pattern) |
| test_dm_detailed_balance.jl | 4 | 2 (already correct) | ~2 more |
| test_dm_scaling.jl | 12 | 0 (uses println) | ~12 (replace println) |
| test_diagnostics.jl | ~47 | 0 | ~15 (many structural) |
| test_gns_trajectory.jl | 14 | 4 (already present) | ~6 more |
| test_krylov_crossvalidation.jl | 4 | 0 (uses Printf) | ~4 (keep Printf tables) |
| test_krylov_eigsolve.jl | 21 | 0 | ~10 |
| test_krylov_matvec.jl | 26 | 0 | ~6 (loop-summary for round-trips) |
| test_observable_trajectories.jl | 18 | 0 | ~2 (mostly structural) |
| test_regression.jl | 4 | 0 | ~4 |
| test_results.jl | 10 | 0 | ~0 (all exact round-trips, atol=0) |
| test_threading.jl | 3 | 5 (skip messages) | ~2 |
| test_trajectory_fixes.jl | 6 | 0 | ~4 |
| test_workspace_independence.jl | 7 | 0 | ~4 |

### Threshold Categories (for audit)

| Category | Count | Current tolerance | Theory-based tightening? |
|----------|-------|-------------------|--------------------------|
| Machine precision (algebraic identity) | ~40 | 1e-10 to 1e-14 | Yes: atol should reflect DIM * eps(Float64) |
| Quadrature approximation | ~8 | TOL_QUADRATURE (1e-6) | Already theory-based |
| Statistical noise (trajectory avg) | ~10 | 0.01 to 0.05 | Theory: 1/sqrt(N_traj) |
| Convergence ratios | ~12 | [1.5, 2.5] or [3.0, 5.0] | Theory: exact ratio +/- sub-leading |
| Krylov vs dense | ~12 | 1e-6 to 1e-8 | Theory: KrylovKit tol * safety factor |
| Round-trip serialization | ~15 | atol=0 (exact) | Already tight |
| Domain approximation error | ~8 | 1e-5 to 0.5 | Partly theory-based |
| Allocation budgets | ~8 | 0 to formula-based | Already well-reasoned |
| Trace/norm preservation | ~20 | 1e-6 to 1e-12 | Yes: machine eps * DIM |

### Files to Delete
```
test/old_tests/B_test.jl           (2862 bytes) - replaced by test_dm_scaling.jl
test/old_tests/ham_test.jl         (785 bytes)  - replaced by test_compilation.jl
test/old_tests/kossakowski_test.jl (1751 bytes) - replaced by test_dm_scaling.jl
test/old_tests/log_sobolev_test.jl (140 bytes)  - uses Revise+includet, broken
test/old_tests/time_tests.jl      (3641 bytes) - replaced by test_dm_scaling.jl
test/old_tests/trajectory_test.jl (4897 bytes) - replaced by trajectory_validation/
test/old_tests/trott_test.jl      (1041 bytes) - replaced by test_dm_scaling.jl
```

### Reference Generator Review
`test/reference/generate_references.jl` uses:
- `make_small_liouv_config` -- needs migration to `make_config`
- `SMALL_*` constants -- needs N3_* rename
- `construct_lindbladian` API -- verify still matches current signatures
- Generates 2 BSON files used by test_regression.jl

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| 8 separate factory functions | Single parametrized `make_config` | This phase | Eliminates duplicate code, centralizes defaults |
| `SMALL_*` naming | `N3_*` naming | This phase | Clearer semantics for system size |
| println diagnostics | @info with key=val | This phase | Structured, filterable output |
| Undocumented thresholds | Inline rationale comments | This phase | Maintainable, reviewable |
| Standalone trajectory validation | Integrated into runtests.jl (gated) | This phase | Single test entry point |
| Dead old_tests/ directory | Deleted | This phase | Less confusion |

## Open Questions

1. **Config constructor kwargs flexibility**
   - What we know: Each factory constructs Config with specific kwargs. Lindbladian configs don't pass mixing_time/delta.
   - What's unclear: Whether `Config(sim=Lindbladian(), ..., mixing_time=1.0, delta=0.01)` throws an error or silently ignores extra kwargs.
   - Recommendation: Test this first. If it errors, use conditional splatting. If it silently ignores, the unified factory is trivial.

2. **Trajectory validation runtime**
   - What we know: `run_trajectory_validation.jl` runs 50,000 trajectories per delta per domain (3 deltas x 3 domains + 1 multi-step = ~500k total). `run_convergence_tests.jl` runs 500,000 reference + 4 points x 10 batches = ~540k additional per domain (2 domains).
   - What's unclear: Actual wall-clock time on CI. Likely several minutes each.
   - Recommendation: Measure runtime, gate both behind QUANTUMFURNACE_FULL_TESTS. The convergence tests already use this gate.

3. **Whether n=6 factories should merge into make_config**
   - What we know: `make_n6_liouv_config`, `make_n6_thermalize_config`, `make_n6_test_system` are only used in test_krylov_crossvalidation.jl behind the FULL test gate.
   - What's unclear: Whether they should merge into the unified factories or remain local to that file.
   - Recommendation: Merge them. `make_config(...; num_qubits=6)` and `make_test_system(; num_qubits=6)` are natural extensions.

## Sources

### Primary (HIGH confidence)
- Direct codebase analysis of all 22 test files, test_helpers.jl, runtests.jl
- Julia stdlib Test documentation (Julia 1.12)
- Julia stdlib Logging documentation (@info syntax with keyword arguments)

### Secondary (MEDIUM confidence)
- [Julia Logging documentation](https://docs.julialang.org/en/v1/stdlib/Logging/) -- @info macro syntax and key=value pairs
- [JuliaLogging tutorials](https://julialogging.github.io/tutorials/logging-basics/) -- logging best practices

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- this is pure Julia stdlib (Test, Logging), no external deps needed
- Architecture: HIGH -- all patterns derived from reading the actual codebase, not hypothetical
- Pitfalls: HIGH -- identified from concrete code analysis (e.g., Config constructor kwargs, SMALL_* reference count, loop @info noise)

**Research date:** 2026-02-28
**Valid until:** 2026-03-28 (stable -- Julia test infrastructure changes slowly)
