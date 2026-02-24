# Phase 30: Cross-Validation - Research

**Researched:** 2026-02-24
**Domain:** Krylov spectral gap cross-validation against dense eigen() reference across domains and balance types
**Confidence:** HIGH

## Summary

Phase 30 is a pure validation phase: no new capabilities are added, only a comprehensive test file that cross-validates the Phase 29 Krylov eigensolver (`krylov_spectral_gap`) against the dense `eigen()`-based reference (`extract_leading_eigendata`) across all four domains (EnergyDomain, TimeDomain, TrotterDomain, BohrDomain), both balance types (KMS, GNS), and at two system sizes (n=4 for tight tolerance, n=6 for production-scale). Additionally, the L-vs-E consistency test validates that the faithful Chen CPTP channel eigenvalues converge to the Lindbladian eigenvalues with O(delta^2) error across multiple delta values.

The codebase already has all the building blocks: `krylov_spectral_gap()` dispatches on config type (LiouvConfig for Lindbladian `:LR` path, ThermalizeConfig for channel `:LM` path), `construct_lindbladian()` builds the full dense Liouvillian (feasible at n=4 with dim^2=256 and n=6 with dim^2=4096), and `extract_leading_eigendata()` extracts leading eigenvalues via dense `eigen()`. The test helpers provide factory functions (`make_liouv_config`, `make_liouv_config_gns`, `make_thermalize_config`) and pre-loaded test systems (TEST_HAM at n=4, TEST_JUMPS, TEST_GIBBS). An n=6 system will need to be constructed at test time using `_load_test_hamiltonian` or a new factory -- the Hamiltonian BSON file exists at `hamiltonians/heis_disordered_periodic_n6.bson`.

**Primary recommendation:** Create a single `test/test_krylov_crossvalidation.jl` with `@testset` blocks for n=4 (always-run) and n=6 (gated behind `QUANTUMFURNACE_FULL_TESTS=true`), a shared `compare_krylov_dense` helper, L-vs-E convergence table for KMS, and KMS-vs-GNS comparison. Print diagnostic summaries on every test run.

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

#### L-vs-E consistency testing
- Test multiple delta values (0.1, 0.01, 0.001) to demonstrate O(delta^2) convergence of channel-to-Lindbladian gap mapping
- Print a formatted convergence table: delta | gap_L | gap_from_E | error | order
- Hard assertion: convergence order must be >= 1.5 (test fails otherwise)
- Cover all 4 domains (EnergyDomain, TimeDomain, TrotterDomain, BohrDomain) with KMS balance
- L-vs-E convergence test runs for KMS only (channel math is balance-independent)

#### GNS domain coverage
- GNS tested across all 4 domains at n=4
- GNS NOT tested at n=6 (n=4 is sufficient for balance-type correctness)
- KMS-vs-GNS comparison: each must match its own dense eigen() reference gap -- the KMS/GNS relationship itself is not asserted
- L-vs-E consistency test runs KMS only, not GNS

#### Diagnostic output
- On test failure: print top-k eigenvalue table from both Krylov and dense, with absolute/relative errors
- On test success: always print one-line summary per test (domain, gap_krylov, gap_dense, error)
- L-vs-E convergence: print formatted convergence table (delta, gap_L, gap_from_E, error, estimated order)
- No wall-clock timing in Phase 30 -- correctness only; Phase 31 handles benchmarking

#### Test organization
- New dedicated test file (e.g., test_krylov_crossvalidation.jl), separate from Phase 29 unit tests
- Single file with separate @testset blocks for n=4 and n=6
- Shared helper function (e.g., compare_krylov_dense) to run both methods and compare, reducing duplication across domain/balance combos
- n=6 tests gated behind `QUANTUMFURNACE_FULL_TESTS=true` environment variable -- skipped in CI, run manually for validation

### Claude's Discretion
- Exact helper function signatures and internal structure
- How many top-k eigenvalues to print on failure
- @testset nesting structure within the file
- KrylovKit parameters (krylovdim, tol) for cross-validation runs

### Deferred Ideas (OUT OF SCOPE)
None -- discussion stayed within phase scope
</user_constraints>

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| KrylovKit.jl | 0.8/0.9/0.10 | Matrix-free eigsolve via Arnoldi | Already in Project.toml (Phase 29); provides `krylov_spectral_gap` |
| LinearAlgebra (stdlib) | -- | Dense `eigen()`, `norm()`, matrix ops | Already used throughout; `extract_leading_eigendata` uses it |
| Test (stdlib) | -- | `@testset`, `@test`, `@test_throws` | Standard Julia test framework |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| Printf (stdlib) | -- | `@printf` or `@sprintf` for formatted convergence tables | Diagnostic output printing |
| BSON.jl | (existing) | Load n=6 Hamiltonian from disk | n=6 test system construction |

### Alternatives Considered
None -- this phase uses only existing project infrastructure; no new libraries needed.

**Installation:**
No new dependencies. All required packages are already in Project.toml.

## Architecture Patterns

### Recommended File Structure
```
test/
├── test_helpers.jl                    # [EXISTS] Shared fixtures, factory functions
├── test_krylov_matvec.jl              # [EXISTS] Phase 27-28 matvec tests
├── test_krylov_eigsolve.jl            # [EXISTS] Phase 29 unit tests
├── test_krylov_crossvalidation.jl     # [NEW] Phase 30 cross-validation
└── runtests.jl                        # [MODIFY] Add include
```

### Pattern 1: compare_krylov_dense Helper
**What:** A shared helper that runs both Krylov and dense eigsolve on the same (config, hamiltonian, jumps) and returns structured comparison results.
**When to use:** Every individual domain/balance cross-validation test.
**Recommended signature:**
```julia
function compare_krylov_dense(config_liouv, hamiltonian, jumps;
    trotter=nothing,
    krylovdim=30,
    howmany=4,
    tol=1e-10,
    n_dense_modes=4,
)
    # 1. Dense reference
    L_dense = construct_lindbladian(jumps, config_liouv, hamiltonian; trotter=trotter)
    dense_result = extract_leading_eigendata(L_dense; n_modes=n_dense_modes)

    # 2. Krylov result
    krylov_result = krylov_spectral_gap(config_liouv, hamiltonian, jumps;
        trotter=trotter, krylovdim=krylovdim, howmany=howmany, tol=tol)

    return (; krylov_result, dense_result, L_dense)
end
```
**Key insight:** The helper returns all intermediate results so the calling testset can make custom assertions. The L_dense is returned for the L-vs-E convergence path to avoid recomputing it.

### Pattern 2: Diagnostic Printing on Success/Failure
**What:** Print a one-line summary on every test, and an eigenvalue comparison table on failure.
**When to use:** Every compare_krylov_dense invocation.
**Example:**
```julia
function print_gap_summary(domain_name, balance_name, krylov_gap, dense_gap)
    err = abs(krylov_gap - dense_gap)
    rel_err = err / dense_gap
    @printf("  %-15s %-5s | gap_krylov=%.8e  gap_dense=%.8e  err=%.2e  rel=%.2e\n",
            domain_name, balance_name, krylov_gap, dense_gap, err, rel_err)
end

function print_eigenvalue_table(krylov_eigs, dense_eigs; top_k=6)
    println("  Top-$top_k eigenvalue comparison:")
    println("  idx | Re(krylov)        | Re(dense)         | abs_err      | rel_err")
    println("  " * "-"^75)
    for k in 1:min(top_k, length(krylov_eigs), length(dense_eigs))
        kr = real(krylov_eigs[k])
        dr = real(dense_eigs[k])
        ae = abs(kr - dr)
        re = ae / max(abs(dr), 1e-30)
        @printf("  %3d | %+.10e | %+.10e | %.2e | %.2e\n", k, kr, dr, ae, re)
    end
end
```

### Pattern 3: L-vs-E Convergence Table
**What:** For multiple delta values, compute gap from Lindbladian path and channel path, show convergence order.
**When to use:** XVAL-03 L-vs-E consistency test.
**Example:**
```julia
function run_le_convergence(config_factory, hamiltonian, jumps; deltas=[0.1, 0.01, 0.001], trotter=nothing)
    # Get Lindbladian reference gap (delta-independent)
    config_liouv = config_factory(; with_coherent=true)
    gap_L = krylov_spectral_gap(config_liouv, hamiltonian, jumps;
        trotter=trotter, krylovdim=30, howmany=4).spectral_gap

    rows = []
    for delta in deltas
        config_therm = make_thermalize_config(domain; with_coherent=true, delta=delta)
        gap_from_E = krylov_spectral_gap(config_therm, hamiltonian, jumps;
            trotter=trotter, krylovdim=30, howmany=4).spectral_gap
        error = abs(gap_L - gap_from_E)
        push!(rows, (; delta, gap_L, gap_from_E, error))
    end

    # Compute convergence orders from consecutive pairs
    orders = Float64[]
    for i in 2:length(rows)
        if rows[i-1].error > 0 && rows[i].error > 0
            ratio = log(rows[i-1].error / rows[i].error) / log(rows[i-1].delta / rows[i].delta)
            push!(orders, ratio)
        end
    end

    return (; rows, orders)
end
```

### Pattern 4: n=6 Test System Construction
**What:** Load the n=6 Hamiltonian BSON and create jump operators, analogous to test_helpers.jl's `make_test_system()`.
**When to use:** Inside the n=6 `@testset` block.
**Example:**
```julia
function make_n6_test_system(; trotter=nothing)
    n_qubits = 6
    source_root = dirname(@__DIR__)
    ham_path = joinpath(source_root, "hamiltonians", "heis_disordered_periodic_n$(n_qubits).bson")
    hamiltonian = _load_test_hamiltonian(ham_path, BETA)

    jump_paulis = [[X], [Y], [Z]]
    jump_sites = 1:n_qubits
    num_of_jumps = length(jump_paulis) * length(jump_sites)
    jump_normalization = sqrt(num_of_jumps)

    basis_unitary = trotter !== nothing ? trotter.eigvecs : hamiltonian.eigvecs

    jumps = JumpOp[]
    for pauli in jump_paulis
        for site in jump_sites
            jump_op = Matrix(pad_term(pauli, n_qubits, site)) ./ jump_normalization
            jump_in_eigen = basis_unitary' * jump_op * basis_unitary
            orthogonal = (jump_op == transpose(jump_op))
            herm = (jump_op == jump_op')
            push!(jumps, JumpOp(jump_op, jump_in_eigen, orthogonal, herm))
        end
    end

    gibbs = hamiltonian.gibbs
    return (; hamiltonian, jumps, gibbs, n_qubits)
end
```
**Key insight:** The n=6 system has dim=64, dim^2=4096. Dense `eigen()` on a 4096x4096 complex matrix is feasible but takes seconds. The Krylov method at n=6 with krylovdim=30 requires 30 vectors of 4096 elements each = ~1MB, well within memory bounds. The dense Lindbladian construction (`construct_lindbladian`) allocates a 4096x4096 matrix = ~256MB of ComplexF64, which is fine.

### Pattern 5: Env-Gated n=6 Tests
**What:** Skip n=6 tests unless `QUANTUMFURNACE_FULL_TESTS=true`.
**When to use:** n=6 testset block.
**Example (matching existing pattern in run_convergence_tests.jl):**
```julia
if get(ENV, "QUANTUMFURNACE_FULL_TESTS", "false") == "true"
    @testset "n=6 cross-validation (KMS, all domains)" begin
        # ... n=6 tests
    end
else
    @info "Skipping n=6 cross-validation (set QUANTUMFURNACE_FULL_TESTS=true to run)"
end
```

### Anti-Patterns to Avoid
- **Computing dense Lindbladian at n=6 without gating:** A 4096x4096 dense eigendecomposition takes seconds per domain. Running it for all 4 domains on every CI run would add 20+ seconds. Always gate behind the env var.
- **Hardcoding tolerance without explanation:** Each tolerance should have a comment explaining its source (1e-8 for n=4 KMS comes from KrylovKit convergence at tol=1e-10 with room for numerical noise; 1e-6 for n=6 KMS allows for larger Krylov subspace approximation error at higher dimension).
- **Using `rtol` for L-vs-E convergence assertion:** The convergence order assertion should use the computed order, not a relative tolerance on the gap itself. The error between gap_L and gap_from_E is O(delta^2) by construction.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Dense eigendecomposition | Custom power iteration | `extract_leading_eigendata` (DIAG-01) | Already implemented in `diagnostics.jl`, handles sorting, left/right eigenvectors |
| Dense Lindbladian construction | Manual kron accumulation | `construct_lindbladian` from `furnace.jl` | Already implemented for all 4 domains, both KMS and GNS |
| n=6 Hamiltonian loading | Custom BSON parsing | `_load_test_hamiltonian` from `test_helpers.jl` | Already handles legacy BSON field order, bohr_dict construction |
| Jump operator creation | Manual Pauli kronecker products | `pad_term` + JumpOp constructor pattern from `make_test_system` | Established pattern that handles eigenbasis transformation |

**Key insight:** Phase 30 should not introduce any new computational functions -- only test infrastructure. Every computation needed already exists in the project.

## Common Pitfalls

### Pitfall 1: GNS with_coherent Mismatch
**What goes wrong:** Creating a ThermalizeConfigGNS or LiouvConfigGNS with `with_coherent=true` throws an error.
**Why it happens:** The GNS constructors have a hard `with_coherent && error("GNS configs must have with_coherent=false")` guard (see `structs.jl:124`).
**How to avoid:** Always use `with_coherent=false` for GNS configs. The existing `make_liouv_config_gns` already does this.
**Warning signs:** `ErrorException: GNS configs must have with_coherent=false` at config construction time.

### Pitfall 2: TrotterDomain Requires Separate Jumps
**What goes wrong:** Using TEST_JUMPS with TrotterDomain produces wrong results because the jump eigenbases are in the Hamiltonian eigenbasis, not the Trotter eigenbasis.
**Why it happens:** TrotterDomain uses `trotter.eigvecs` as its working basis, but TEST_JUMPS were created with `hamiltonian.eigvecs`. The eigenbases differ because Trotterization introduces approximation error.
**How to avoid:** Use `TEST_TROTTER_JUMPS` (and `TEST_TROTTER`) for TrotterDomain tests, as established in test_helpers.jl. For n=6, create a separate set of Trotter jumps using `make_n6_test_system(; trotter=trotter_n6)`.
**Warning signs:** Krylov gap does not match dense gap for TrotterDomain, even though other domains pass.

### Pitfall 3: L-vs-E Convergence Order Sensitive to Absolute Error Floor
**What goes wrong:** At very small delta (e.g., 0.001), the O(delta^2) error becomes comparable to KrylovKit convergence tolerance, and the computed convergence order drops below 2.
**Why it happens:** The channel eigenvalue mapping error is `O(delta^2)`, but KrylovKit's absolute error is `O(tol)` = O(1e-10). At delta=0.001, the O(delta^2) error is ~1e-6, which is well above the KrylovKit floor. But if delta were 1e-5, the O(delta^2) error would be ~1e-10, at the KrylovKit noise floor.
**How to avoid:** Use delta values that keep the O(delta^2) error well above the KrylovKit convergence floor. The user-specified range [0.1, 0.01, 0.001] is appropriate: O(delta^2) errors are [1e-2, 1e-4, 1e-6], all well above 1e-10. The hard assertion threshold of 1.5 (rather than 2.0) provides margin for sub-leading terms.
**Warning signs:** Computed convergence order < 1.5 at the smallest delta. Fix: remove the smallest delta from the sequence.

### Pitfall 4: Dense Lindbladian Construction Prints to stdout
**What goes wrong:** `construct_lindbladian` calls `println("Constructing Liouvillian ...")` and `@printf("Done.\n")` (see `furnace.jl:46,12`), cluttering test output.
**Why it happens:** The function was written for interactive use, not test infrastructure.
**How to avoid:** Accept the stdout noise -- it is harmless in test output and not worth silencing. Alternatively, redirect stdout in a `let` block, but this adds complexity for no functional benefit.
**Warning signs:** "Constructing Liouvillian (Energy)" lines interspersed with test output.

### Pitfall 5: n=6 Config Factory Needs num_qubits=6
**What goes wrong:** Using the existing `make_liouv_config(EnergyDomain())` for n=6 creates a config with `num_qubits=4`, producing wrong Krylov memory estimates and potentially wrong precomputation.
**Why it happens:** The test_helpers.jl factory functions hardcode `NUM_QUBITS = 4`.
**How to avoid:** Create n=6-specific factory functions (or pass num_qubits as a parameter) that use `num_qubits=6`. The n=6 config should use the same physical parameters (BETA=10, SIGMA=0.1, etc.) for consistency, just with `num_qubits=6`.
**Warning signs:** Krylov result at n=6 is suspiciously fast (dim=16 instead of 64) or dimension mismatch errors.

### Pitfall 6: Channel Path with ThermalizeConfigGNS Not Yet Tested
**What goes wrong:** The `krylov_spectral_gap(::AbstractThermalizeConfig, ...)` method dispatches on the abstract type, so ThermalizeConfigGNS should work. But `make_thermalize_config_gns` does not exist in test_helpers.jl for n=4 (only `make_small_thermalize_config_gns` for n=3).
**Why it happens:** Phase 29 only tested the channel path with ThermalizeConfig (KMS), not ThermalizeConfigGNS.
**How to avoid:** The L-vs-E test runs KMS only per user decision, so this is not needed for Phase 30. The GNS tests use the Lindbladian path (LiouvConfigGNS), not the channel path. If channel-path GNS testing is needed in the future, a `make_thermalize_config_gns` factory for n=4 would need to be added. For now, this is not in scope.
**Warning signs:** N/A -- intentionally out of scope per CONTEXT.md.

## Code Examples

### Existing Cross-Validation Pattern (from test_krylov_eigsolve.jl)
```julia
# Source: test/test_krylov_eigsolve.jl, Testset 3
config = make_liouv_config(EnergyDomain(); with_coherent=true)
L_dense = construct_lindbladian(TEST_JUMPS, config, TEST_HAM)
dense_result = extract_leading_eigendata(L_dense; n_modes=4)
result = krylov_spectral_gap(config, TEST_HAM, TEST_JUMPS; krylovdim=30, howmany=4)
@test isapprox(result.spectral_gap, dense_result.spectral_gap; rtol=1e-6)
```

### Channel Eigsolve with ThermalizeConfig (from test_krylov_eigsolve.jl)
```julia
# Source: test/test_krylov_eigsolve.jl, Testset 5
config_therm = make_thermalize_config(EnergyDomain(); with_coherent=true, delta=0.01)
config_liouv = make_liouv_config(EnergyDomain(); with_coherent=true)
L_dense = construct_lindbladian(TEST_JUMPS, config_liouv, TEST_HAM)
dense_result = extract_leading_eigendata(L_dense; n_modes=4)
result = krylov_spectral_gap(config_therm, TEST_HAM, TEST_JUMPS; krylovdim=30, howmany=4)
# O(delta^2) tolerance due to faithful Chen channel eigenvalue mapping
@test isapprox(result.spectral_gap, dense_result.spectral_gap; rtol=2e-3)
```

### Existing Env-Gated Test Pattern (from run_convergence_tests.jl)
```julia
# Source: test/trajectory_validation/run_convergence_tests.jl, line 110
if get(ENV, "QUANTUMFURNACE_FULL_TESTS", "false") == "true"
    @testset "TVAL-05: EnergyDomain 1/sqrt(N) convergence" begin
        # ... expensive tests
    end
else
    @info "Skipping TVAL-05 convergence tests (set QUANTUMFURNACE_FULL_TESTS=true to run)"
end
```

### n=4 Config Factories (from test_helpers.jl)
```julia
# Source: test/test_helpers.jl
# All factory functions use: NUM_QUBITS=4, BETA=10.0, SIGMA=0.1, W0=0.05, etc.
make_liouv_config(domain; with_coherent=true)      # LiouvConfig (KMS)
make_liouv_config_gns(domain)                       # LiouvConfigGNS (GNS, with_coherent=false)
make_thermalize_config(domain; with_coherent=true, delta=0.01)  # ThermalizeConfig (KMS)
```

### TrotterDomain Test System (from test_helpers.jl)
```julia
# Source: test/test_helpers.jl
const TEST_TROTTER = TrottTrott(TEST_HAM, T0, NUM_TROTTER_STEPS_PER_T0)
const TEST_TROTTER_JUMPS = make_test_system(; trotter=TEST_TROTTER).jumps
# Usage: krylov_spectral_gap(config, TEST_HAM, TEST_TROTTER_JUMPS; trotter=TEST_TROTTER)
```

## Specific Implementation Guidance

### n=4 Test Matrix (XVAL-01 + XVAL-04)

All combinations to test at n=4:

| Domain | Balance | Config Type | Jumps | Trotter | Tolerance | Req |
|--------|---------|-------------|-------|---------|-----------|-----|
| EnergyDomain | KMS | LiouvConfig | TEST_JUMPS | -- | 1e-8 | XVAL-01 |
| TimeDomain | KMS | LiouvConfig | TEST_JUMPS | -- | 1e-8 | XVAL-01 |
| TrotterDomain | KMS | LiouvConfig | TEST_TROTTER_JUMPS | TEST_TROTTER | 1e-8 | XVAL-01 |
| BohrDomain | KMS | LiouvConfig | TEST_JUMPS | -- | 1e-8 | XVAL-01 |
| EnergyDomain | GNS | LiouvConfigGNS | TEST_JUMPS | -- | 1e-8 | XVAL-04 |
| TimeDomain | GNS | LiouvConfigGNS | TEST_JUMPS | -- | 1e-8 | XVAL-04 |
| TrotterDomain | GNS | LiouvConfigGNS | TEST_TROTTER_JUMPS | TEST_TROTTER | 1e-8 | XVAL-04 |
| BohrDomain | GNS | LiouvConfigGNS | TEST_JUMPS | -- | 1e-8 | XVAL-04 |

### n=6 Test Matrix (XVAL-02)

KMS only, all 4 domains. Gated behind env var.

| Domain | Balance | Tolerance | Notes |
|--------|---------|-----------|-------|
| EnergyDomain | KMS | 1e-6 | Larger dim = slightly looser tolerance |
| TimeDomain | KMS | 1e-6 | |
| TrotterDomain | KMS | 1e-6 | Needs n=6 Trotter + Trotter jumps |
| BohrDomain | KMS | 1e-6 | |

**n=6 dimensions:** dim=64, dim^2=4096. Dense Liouvillian is 4096x4096 = ~256MB. Dense `eigen()` on this size takes a few seconds. `krylov_spectral_gap` with krylovdim=30 requires 30 vectors of 4096 elements = ~1MB.

**n=6 krylovdim guidance:** May need krylovdim=50 instead of 30 for n=6 due to larger spectral spread. Start with 30; if convergence fails, the retry mechanism (30->45->68) will handle it automatically.

### L-vs-E Convergence Test (XVAL-03)

For each of the 4 domains with KMS balance:
1. Compute `gap_L` via `krylov_spectral_gap` with `LiouvConfig` (Lindbladian path)
2. For each `delta` in [0.1, 0.01, 0.001]:
   - Compute `gap_from_E` via `krylov_spectral_gap` with `ThermalizeConfig` (channel path)
   - Compute `error = abs(gap_L - gap_from_E)`
3. Print formatted convergence table
4. Compute convergence order from consecutive error pairs: `order = log(e1/e2) / log(d1/d2)`
5. Assert order >= 1.5 for each consecutive pair

**Mathematical basis:** The faithful Chen channel eigenvalues satisfy `mu_E = 1 + delta*lambda_L + O(delta^2)`. The linear conversion `lambda_L_approx = (mu_E - 1)/delta` recovers `lambda_L + O(delta)` in the eigenvalue, giving `gap_from_E = gap_L + O(delta^2)` (the leading O(delta) cancels because both the numerator and denominator scale with delta, leaving a quadratic residual). The O(delta^2) scaling has been empirically confirmed in quick-36: at delta=0.01, the relative error was ~0.134%.

**Convergence table format:**
```
L-vs-E convergence (EnergyDomain, KMS):
  delta    | gap_L            | gap_from_E       | error        | order
  ---------|------------------|------------------|--------------|------
  1.00e-01 | 1.23456789e+00   | 1.23567890e+00   | 1.11e-03     |  --
  1.00e-02 | 1.23456789e+00   | 1.23458789e+00   | 2.00e-05     | 1.74
  1.00e-03 | 1.23456789e+00   | 1.23456809e+00   | 2.00e-07     | 2.00
```

### KrylovKit Parameters for Cross-Validation

**Recommendation:**
- n=4: `krylovdim=30, howmany=4, tol=1e-10` (same as Phase 29 tests; sufficient for 1e-8 accuracy)
- n=6: `krylovdim=30, howmany=4, tol=1e-10` (start here; retry mechanism handles convergence issues)
- L-vs-E channel path: `krylovdim=30, howmany=4, tol=1e-10` (same; channel eigsolve uses :LM targeting which typically converges faster than :LR)

### Failure Diagnostics: Top-k Eigenvalue Table

**Recommendation:** Print top 6 eigenvalues on failure. This covers the steady state + gap mode + 4 additional modes, which is enough to diagnose sorting issues, convergence failures, or systematic offsets.

```julia
function on_failure_diagnostics(krylov_result, dense_result)
    println("  FAILURE DIAGNOSTICS:")
    print_eigenvalue_table(krylov_result.eigenvalues, dense_result.eigenvalues; top_k=6)
    @printf("  Krylov converged: %d, matvec_count: %d, restarts: %d\n",
            krylov_result.converged, krylov_result.matvec_count, krylov_result.num_restarts)
end
```

### @testset Nesting Structure

**Recommendation:**
```
@testset "Krylov Cross-Validation" begin
    @testset "n=4 KMS (all domains)" begin
        @testset "EnergyDomain" ... end
        @testset "TimeDomain" ... end
        @testset "TrotterDomain" ... end
        @testset "BohrDomain" ... end
    end
    @testset "n=4 GNS (all domains)" begin
        @testset "EnergyDomain" ... end
        @testset "TimeDomain" ... end
        @testset "TrotterDomain" ... end
        @testset "BohrDomain" ... end
    end
    @testset "L-vs-E convergence (KMS)" begin
        @testset "EnergyDomain" ... end
        @testset "TimeDomain" ... end
        @testset "TrotterDomain" ... end
        @testset "BohrDomain" ... end
    end
    # n=6 block (gated)
    if get(ENV, "QUANTUMFURNACE_FULL_TESTS", "false") == "true"
        @testset "n=6 KMS (all domains)" begin
            @testset "EnergyDomain" ... end
            @testset "TimeDomain" ... end
            @testset "TrotterDomain" ... end
            @testset "BohrDomain" ... end
        end
    else
        @info "Skipping n=6 cross-validation (set QUANTUMFURNACE_FULL_TESTS=true to run)"
    end
end
```

### Config Factory Functions for n=6

The n=6 tests need config factories with `num_qubits=6`. These can be simple local functions in the test file:

```julia
function make_n6_liouv_config(domain; with_coherent=true)
    LiouvConfig(
        num_qubits = 6,
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
        num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
    )
end
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Phase 29 unit tests: per-domain spot checks at n=4 with rtol=1e-6 | Phase 30: systematic cross-validation across all domain/balance combos with tight tolerances (1e-8) and convergence analysis | Phase 30 (this) | Establishes quantitative trust for n>6 production use |
| No L-vs-E consistency testing | O(delta^2) convergence table with order assertion | Phase 30 (this) | Validates the channel-to-Lindbladian eigenvalue conversion for the faithful Chen channel |
| n=4 only | n=4 tight + n=6 production-scale | Phase 30 (this) | Validates Krylov at the scale boundary where dense becomes expensive |

**Deprecated/outdated:**
- The Phase 29 unit tests in `test_krylov_eigsolve.jl` remain as fast CI unit tests. Phase 30's cross-validation complements them, not replaces them.

## Open Questions

1. **n=6 TrotterDomain test system construction time**
   - What we know: `TrottTrott(ham, T0, num_trotter_steps_per_t0)` requires matrix exponentiation for n=6 (64x64), which should be fast. The n=6 jump operator creation requires 18 Pauli kronecker products (3 Paulis x 6 sites), each 64x64, plus eigenbasis transformation -- also fast.
   - What's unclear: Whether the total n=6 test system construction + dense Lindbladian build + Krylov eigsolve can complete within a reasonable timeout (say, 60 seconds per domain).
   - Recommendation: Profile during implementation. If TrotterDomain at n=6 is too slow, document it and consider deferring to a separate standalone validation script (not in the env-gated CI block).

2. **n=6 Krylov convergence at default krylovdim=30**
   - What we know: n=4 converges easily at krylovdim=30. n=6 has a 16x larger Liouville space (4096 vs 256). The Lindbladian spectrum may have more near-degenerate modes, making convergence harder.
   - What's unclear: Whether krylovdim=30 suffices for 4 eigenvalues at n=6, or whether the retry mechanism (30->45->68->102) will activate.
   - Recommendation: Start with krylovdim=30 and let the retry mechanism handle it. If all 4 domains consistently trigger retries, consider using krylovdim=50 as default for n=6 tests.

## Sources

### Primary (HIGH confidence)
- Codebase: `src/krylov_eigsolve.jl` -- `krylov_spectral_gap`, `apply_delta_channel!`, `KrylovGapResult`
- Codebase: `src/krylov_workspace.jl` -- `KrylovWorkspace`, ThermalizeConfig constructor with channel precomputation
- Codebase: `src/krylov_matvec.jl` -- `apply_lindbladian!` for all 4 domains
- Codebase: `src/diagnostics.jl` -- `extract_leading_eigendata` (dense reference)
- Codebase: `src/furnace.jl` -- `construct_lindbladian` (dense Liouvillian)
- Codebase: `src/structs.jl` -- Config type hierarchy, domain types
- Codebase: `test/test_helpers.jl` -- All factory functions, TEST_HAM, TEST_JUMPS, TEST_TROTTER, etc.
- Codebase: `test/test_krylov_eigsolve.jl` -- Phase 29 unit tests (pattern reference)
- Codebase: `test/trajectory_validation/run_convergence_tests.jl` -- Env-gated test pattern
- Codebase: `hamiltonians/heis_disordered_periodic_n6.bson` -- n=6 Hamiltonian data
- Phase 29 research: `.planning/phases/29-eigensolver-integration/29-RESEARCH.md`
- Quick-36 summary: `.planning/quick/36-fix-apply-delta-channel-to-use-faithful-/36-SUMMARY.md`

### Secondary (MEDIUM confidence)
- None needed -- all research is from primary codebase sources

### Tertiary (LOW confidence)
- None

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- no new libraries; all existing infrastructure is well-tested
- Architecture: HIGH -- test structure follows established patterns (test_helpers, @testset nesting, env gating)
- Dense reference correctness: HIGH -- `extract_leading_eigendata` uses stdlib `eigen()` which is LAPACK-backed
- L-vs-E convergence math: HIGH -- O(delta^2) error empirically confirmed in quick-36, mathematically derived from Chen channel structure
- n=6 feasibility: MEDIUM -- dense 4096x4096 eigen() should take seconds but has not been profiled in this specific test context
- Pitfalls: HIGH -- all identified from code inspection of existing infrastructure

**Research date:** 2026-02-24
**Valid until:** 2026-03-24 (internal codebase infrastructure; stable)
