# Phase 3: DM Reference Test Suite - Research

**Researched:** 2026-02-14
**Domain:** Density matrix Lindbladian correctness testing across approximation domains
**Confidence:** HIGH

## Summary

This phase builds a comprehensive test suite verifying the density matrix (DM) simulation as ground truth for all four approximation domains: BohrDomain, EnergyDomain, TimeDomain, and TrotterDomain. The tests establish three categories of verification: (1) detailed balance and domain error hierarchy, (2) DM step error scaling with coherent term and OFT consistency, and (3) Aqua.jl package quality checks.

The codebase already has all the necessary infrastructure for these tests: `construct_lindbladian` builds the full Liouvillian for any domain, `run_lindbladian` finds its fixed point via Arpack, and the coherent term B operators (`coherent_bohr`, `B_time`, `B_trotter`) plus OFT functions (`oft!`, `time_oft!`, `trotter_oft!`) are fully implemented. The test helpers from Phase 1/2 provide `make_liouv_config`, `make_thermalize_config`, `TEST_HAM`, `TEST_JUMPS`, `TEST_GIBBS`, and `TEST_TROTTER`. The key insight is that the Bohr domain with coherent term B produces the exact Gibbs state as its Lindbladian fixed point (this is the theoretically exact construction from Chen 2023). Each successive domain (Energy, Time, Trotter) introduces additional approximation errors that should form a monotonic hierarchy.

For the DM step error scaling tests (DMTST-03/04), the `run_thermalization` function evolves a density matrix via the weak-measurement channel (Kraus maps). Single-step error should scale as O(delta^2) and multi-step accumulated error as O(delta) -- these are empirical verification of Chen Theorem III.1. These tests require running thermalization with several delta values and fitting the observed error scaling.

**Primary recommendation:** Use the 4-qubit disordered Heisenberg system (already the test fixture) for all DM tests. For the detailed balance test (DMTST-01), also test with a 3-qubit system as specified in the success criteria. Use `run_lindbladian` for fixed-point tests and `run_thermalization` for step-error tests. Group Aqua.jl checks in a minimal separate test file.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Test | stdlib | Julia test framework | Used throughout existing test suite |
| LinearAlgebra | stdlib | Matrix norms, trace, eigen | Already in use for all DM operations |
| Aqua | test-only (extras) | Package quality checks | Already in Project.toml extras/test targets |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| BSON | deps | Load Hamiltonian fixtures | Loading n=3 and n=4 Hamiltonian .bson files |
| Random | stdlib | Reproducible initial states | Generating random initial DMs for step-error tests |

### Already Available (No Installation Needed)
All required libraries are already in Project.toml. Aqua is already listed in `[extras]` and `[targets] test`. No new dependencies needed.

## Architecture Patterns

### Test File Organization

Following the existing pattern from Phase 1/2:

```
test/
  runtests.jl                    # Add includes for new test files
  test_helpers.jl                # Add make_small_test_system() for 3-qubit
  test_dm_detailed_balance.jl    # DMTST-01, DMTST-02 (Plan 03-01)
  test_dm_scaling.jl             # DMTST-03, DMTST-04, DMTST-05, DMTST-06 (Plan 03-02)
  test_aqua.jl                   # TINF-03 (Plan 03-03)
```

### Pattern 1: Lindbladian Fixed Point Test (Detailed Balance)

**What:** Construct the Lindbladian for BohrDomain with `with_coherent=true`, find its fixed point via spectral analysis, verify it matches the Gibbs state.

**When to use:** DMTST-01 (detailed balance verification)

**How it works in the codebase:**
1. Create `LiouvConfig` with `BohrDomain()`, `with_coherent=true`
2. Call `construct_lindbladian(jumps, config, hamiltonian)` to get the dim^2 x dim^2 Liouvillian matrix
3. Find the zero eigenvector (fixed point) via `eigen` (small system) or `eigs` (larger system)
4. Reshape into DM, Hermitianize, normalize trace to 1
5. Compare with `hamiltonian.gibbs` using `trace_distance_h`

**Key code path:** `construct_lindbladian` -> `precompute_data(BohrDomain(), ...)` -> `precompute_coherent_total_B(...)` -> `vectorize_liouvillian_coherent!` + `jump_contribution!(BohrDomain(), ...)` for each jump

**Important:** The existing `run_lindbladian` uses Arpack's `eigs` which needs dim^2 >= 3 (works for n >= 2 qubits). For 3-qubit system (dim=8, dim^2=64) this works fine. But for a cleaner test, use `eigen` from LinearAlgebra directly since the matrix is small enough.

```julia
# Pattern for fixed-point extraction
liouv = construct_lindbladian(jumps, config, hamiltonian)
eig = eigen(liouv)
ss_idx = argmin(abs.(real.(eig.values)))
ss_vec = eig.vectors[:, ss_idx]
ss_dm = reshape(ss_vec, dim, dim)
ss_dm = (ss_dm + ss_dm') / 2
ss_dm ./= tr(ss_dm)
dist = trace_distance_h(Hermitian(ss_dm), hamiltonian.gibbs)
@test dist < 1e-10
```

### Pattern 2: Domain Error Hierarchy Test

**What:** Build Lindbladian fixed points for all four domains with matched parameters, verify distances to Gibbs form a monotonic hierarchy.

**When to use:** DMTST-02

**Hierarchy:** dist(Gibbs, Bohr) <= dist(Gibbs, Energy) <= dist(Gibbs, Time) <= dist(Gibbs, Trotter)

**Key insight from code analysis:** All four domains use `construct_lindbladian` with the same dispatch structure. The config factory `make_liouv_config(domain)` already sets matched parameters. The only difference per domain is:
- BohrDomain: exact Bohr frequency decomposition
- EnergyDomain: energy quadrature approximation (Riemann sum over w0 grid)
- TimeDomain: time-domain OFT via NUFFT prefactors (additional time quadrature)
- TrotterDomain: Trotterized time evolution (additional Trotter error)

**Trotter requires TrottTrott object:** `construct_lindbladian` for TrotterDomain dispatches to `jump_contribution!(::TrotterDomain, ...)` which requires the trotter object. The function signature: `construct_lindbladian(jumps, config, hamiltonian; trotter=trotter)`.

### Pattern 3: DM Step Error Scaling

**What:** Run thermalization for a fixed mixing_time with several delta values, measure trace distance to Gibbs at each step, verify O(delta^2) single-step and O(delta) accumulated error.

**When to use:** DMTST-03, DMTST-04

**How to implement:**
1. Use `run_thermalization(jumps, config, initial_dm, hamiltonian)` with different delta values
2. Single-step error: Run one step, measure ||rho_1 - rho_exact_1||. Use Gibbs state as initial (it should be a fixed point, so error = ||rho_1 - gibbs||)
3. Multi-step error: Run for fixed T with varying delta, measure final ||rho_T - gibbs||

**Delta sweep for scaling verification:**
```julia
deltas = [0.1, 0.05, 0.025, 0.0125]
# For single-step: error should halve when delta halves -> ratio ~ 4 (O(delta^2))
# For multi-step over fixed T: error should halve when delta halves -> ratio ~ 2 (O(delta))
```

**Alternative (simpler, more direct):** Use `jump_contribution!` on density matrices directly rather than the full `run_thermalization`, since the latter includes random jump selection. For DM channel step testing, we want the deterministic channel (apply all jumps, not random sampling). The DM stepping in `run_thermalization` randomly picks one jump per step. For a proper DM step test, we should apply the full Lindbladian channel E(rho) = sum_k channel_k(rho) deterministically.

**Best approach:** Use the Lindbladian matrix directly. Single DM step via Euler: rho(t+delta) = rho(t) + delta * L * vec(rho(t)), where L is the Liouvillian. This is exactly what the theory predicts: single-step error O(delta^2), accumulated O(delta).

### Pattern 4: Coherent Term B Consistency

**What:** Compute B via `coherent_bohr` (Bohr domain), `B_time` (time domain), and `B_trotter` (Trotter domain), verify they agree up to known approximation errors.

**When to use:** DMTST-05

**Code paths for B computation:**
```julia
# Bohr domain B (exact, analytical)
config_bohr = make_liouv_config(BohrDomain())
precomputed_bohr = precompute_data(BohrDomain(), config_bohr, TEST_HAM)
B_bohr = coherent_bohr(TEST_HAM, TEST_JUMPS[1], config_bohr)
rmul!(B_bohr, precomputed_bohr.gamma_norm_factor)

# Time domain B (time quadrature approximation)
config_time = make_liouv_config(TimeDomain())
precomputed_time = precompute_data(TimeDomain(), config_time, TEST_HAM)
B_time_val = B_time(TEST_JUMPS[1], TEST_HAM, precomputed_time.b_minus, precomputed_time.b_plus,
                    config_time.t0, config_time.beta, config_time.sigma)
rmul!(B_time_val, precomputed_time.gamma_norm_factor)

# Trotter domain B (time quadrature + Trotter approximation)
config_trott = make_liouv_config(TrotterDomain())
precomputed_trott = precompute_data(TrotterDomain(), config_trott, TEST_TROTTER)
B_trott = B_trotter(TEST_JUMPS[1], TEST_TROTTER, precomputed_trott.b_minus, precomputed_trott.b_plus,
                     config_trott.beta, config_trott.sigma)
rmul!(B_trott, precomputed_trott.gamma_norm_factor)
```

**Expected hierarchy:** `norm(B_bohr - B_time) <= TOL_QUADRATURE`, `norm(B_bohr - B_trott)` >= `norm(B_bohr - B_time)` (Trotter adds error on top of time quadrature).

**Note:** B_trotter is computed in the Trotter eigenbasis. To compare with B_bohr (in Hamiltonian eigenbasis), need basis transformation: `B_trott_in_eigen = TEST_TROTTER.trafo_from_eigen_to_trotter' * B_trott * TEST_TROTTER.trafo_from_eigen_to_trotter`.

Wait -- checking `bohr_domain.jl` more carefully, `coherent_bohr` takes the jump's `in_eigenbasis` field which is always in Hamiltonian eigenbasis. `B_time` also uses `hamiltonian.eigvals` for the time evolution, so B_time is also in Hamiltonian eigenbasis. But `B_trotter` uses `trotter.eigvals_t0` for evolution, so B_trotter is in the Trotter eigenbasis. For comparison, transform B_trotter to Hamiltonian eigenbasis via `trafo_from_eigen_to_trotter`.

### Pattern 5: OFT Consistency

**What:** For the same jump and energy value, compare `oft!` (analytical Bohr/Energy), `time_oft!` (time-domain numerical), and `trotter_oft!` (Trotter-domain numerical).

**When to use:** DMTST-06

**Code paths (from `src/ofts.jl`):**
```julia
dim = DIM
w = -0.12  # Some test energy
A = TEST_JUMPS[1]
caches = OFTCaches(dim)

# Energy domain OFT (analytical Gaussian filter)
A_energy = Matrix{ComplexF64}(undef, dim, dim)
oft!(A_energy, A, w, TEST_HAM, SIGMA)
energy_prefactor = 1 / sqrt(SIGMA * sqrt(2 * pi))
A_energy .*= energy_prefactor

# Time domain OFT (time quadrature)
config_time = make_liouv_config(TimeDomain())
precomputed_time = precompute_data(TimeDomain(), config_time, TEST_HAM)
# ... need to get time_labels from precomputed_time
A_time = Matrix{ComplexF64}(undef, dim, dim)
time_oft!(A_time, caches, A, w, TEST_HAM, time_labels, SIGMA)
time_prefactor = T0 * sqrt(SIGMA * sqrt(2 / pi) / (2 * pi))
A_time .*= time_prefactor

# Trotter domain OFT
A_trotter = Matrix{ComplexF64}(undef, dim, dim)
trotter_oft!(A_trotter, caches, A, w, TEST_TROTTER, time_labels, SIGMA)
A_trotter .*= time_prefactor
```

**Note on time_labels:** The time labels come from `precompute_labels(TimeDomain(), config)` which returns `(energy_labels, time_labels)`. Then they are truncated for OFT via `truncate_time_labels_for_oft(time_labels, sigma)`. The precomputed_data stores NUFFT prefactors, not the raw time labels. For the OFT test, we need to reconstruct the time labels: `energy_labels = create_energy_labels(NUM_ENERGY_BITS, W0)` then `time_labels = energy_labels .* (T0 / W0)` then truncate.

**Alternatively:** Use `precompute_labels(TimeDomain(), config)` directly and pass the time_labels to `time_oft!` and `trotter_oft!`.

**Trotter OFT basis:** `trotter_oft!` returns the result in the Trotter eigenbasis (it uses `trotter.eigvals_t0`). The jump's `in_eigenbasis` is in Hamiltonian eigenbasis when passed to the function. Wait -- examining the code more carefully, `trotter_oft!` uses `jump.in_eigenbasis` directly. If the jump was created with Hamiltonian eigenvectors (as in the test fixtures), then the OFT result mixes bases. Checking the existing test fixtures: `TEST_JUMPS` are created with `hamiltonian.eigvecs' * jump_op * hamiltonian.eigvecs`, so `in_eigenbasis` is in Hamiltonian eigenbasis. For TrotterDomain, the jumps should arguably be in Trotter basis. However, looking at how `precompute_R` for TrotterDomain works, it uses `jump.in_eigenbasis` with NUFFT prefactors from `ham_or_trott.bohr_freqs` where `ham_or_trott` is a TrottTrott. So the Trotter domain uses the Trotter's bohr_freqs but the jump in Hamiltonian eigenbasis -- this is the intended mixing.

**Important:** For Trotter OFT comparison, the output of `trotter_oft!` needs to be compared with `time_oft!` after accounting for the Trotter approximation error. They do NOT need basis transformation if both use the same `jump.in_eigenbasis`. The error comes from Trotter time evolution (discrete eigenvalues) vs exact time evolution (continuous eigenvalues).

### Pattern 6: Aqua.jl Quality Checks

**What:** Standard package quality tests.

**Implementation:**
```julia
using Aqua
using QuantumFurnace

@testset "Aqua.jl quality checks" begin
    Aqua.test_all(QuantumFurnace;
        ambiguities = false,  # May have legitimate ambiguities from multiple dispatch
        deps_compat = (ignore=[:Base, :Pkg, :Printf, :Random, :LinearAlgebra, :SparseArrays, :SharedArrays, :Distributed],),
        stale_deps = (ignore=[:Pkg],),
        piracies = false,  # May have legitimate extensions
    )
end
```

**Note on expected issues:**
- The package has many stdlib deps in `[deps]` (LinearAlgebra, SparseArrays, Random, etc.) that may not have compat entries -- this is normal for stdlibs
- Method ambiguities may exist due to heavy use of multiple dispatch
- Need to handle Aqua's expectations about stdlib compat entries (stdlibs don't need compat bounds in Julia 1.x)
- Stale deps: `Pkg` is used in the module (`using Pkg`) but might be flagged
- Piracy: The package defines methods on stdlib types (e.g., `kron!` on AbstractMatrix) which could be flagged

### Anti-Patterns to Avoid

- **Using Arpack `eigs` for small systems:** For dim^2 = 64 (3-qubit) or 256 (4-qubit), use `eigen` directly. Arpack is for large sparse problems and can be numerically unstable for small dense matrices.
- **Random initial states for deterministic tests:** Use fixed initial states (maximally mixed, specific eigenstates) rather than random for reproducibility of scaling tests.
- **Not accounting for gamma_norm_factor:** All B computations and Lindbladian contributions include a `gamma_norm_factor` normalization. When comparing B across domains, must apply the same normalization.
- **Comparing B_trotter directly with B_bohr without basis change:** B_trotter is in Trotter eigenbasis; B_bohr is in Hamiltonian eigenbasis. Transform B_trotter via `trafo_from_eigen_to_trotter`.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Lindbladian construction | Custom Liouvillian builder | `construct_lindbladian(jumps, config, hamiltonian)` | Already handles all domains, normalization, coherent terms |
| Fixed point extraction | Custom eigenvalue search | `eigen(liouv)` for small systems | Reliable, deterministic, no convergence issues |
| Trace distance | Manual eigenvalue computation | `trace_distance_h(Hermitian(rho), gibbs)` from qi_tools.jl | Already handles Hermitian input correctly |
| DM thermalization | Step-by-step channel application | `run_thermalization(jumps, config, dm, hamiltonian)` for multi-step | Full pipeline with distance tracking |
| Gibbs state | Manual Boltzmann weight computation | `hamiltonian.gibbs` (precomputed in `finalize_hamham`) | Already normalized, in eigenbasis |
| Package quality | Custom dep/compat checks | `Aqua.test_all(QuantumFurnace; ...)` | Standard, comprehensive, maintained |

**Key insight:** The codebase already has all building blocks; tests just need to call existing functions and verify outputs.

## Common Pitfalls

### Pitfall 1: 3-Qubit vs 4-Qubit System Choice
**What goes wrong:** Success criterion 1 specifies "3-qubit Heisenberg system" but the test fixtures use 4-qubit.
**Why it happens:** 3-qubit has dim=8 (dim^2=64) -- very fast but smaller error margins. 4-qubit has dim=16 (dim^2=256) -- slower but more representative.
**How to avoid:** Add a 3-qubit fixture to test_helpers.jl (the `heis_disordered_periodic_n3.bson` file exists in hamiltonians/). Use 3-qubit for DMTST-01 (exact detailed balance) and 4-qubit for all other tests (better error statistics).
**Warning signs:** If DMTST-01 passes at 1e-10 for 3-qubit but not 4-qubit, it may be a numerical precision issue rather than a bug.

### Pitfall 2: Lindbladian Eigenvalue Extraction
**What goes wrong:** The zero eigenvalue of the Liouvillian may not be exactly zero due to numerical roundoff.
**Why it happens:** `construct_lindbladian` accumulates many terms; floating point errors compound.
**How to avoid:** Find the eigenvalue with smallest |Re(lambda)|, extract its eigenvector, Hermitianize and normalize. The existing `run_lindbladian` does this correctly with `eigs`.
**Warning signs:** Multiple eigenvalues near zero, or the "fixed point" having negative eigenvalues.

### Pitfall 3: EnergyDomain Lindbladian Construction Prints to Stdout
**What goes wrong:** `_pick_transition_kms` calls `@printf("Smooth Metro\n")` during config creation.
**Why it happens:** Debugging prints left in production code.
**How to avoid:** Accept that test output will be noisy. Not a bug, just cosmetic. Don't suppress it (it confirms the correct transition type is selected).
**Warning signs:** Unexpected "Smooth Metro" or "Gaussian" print during tests.

### Pitfall 4: Normalization Factors Differ Between Domains
**What goes wrong:** B_bohr and B_time have different normalization conventions.
**Why it happens:** `precompute_coherent_total_B` handles this internally (multiplies by `gamma_norm_factor`), but if calling `coherent_bohr` or `B_time` directly, you must apply the factor manually.
**How to avoid:** Always use `precompute_coherent_total_B` for total B, or apply `rmul!(B, precomputed_data.gamma_norm_factor)` after direct calls (as seen in B_test.jl).
**Warning signs:** B values differ by orders of magnitude rather than small approximation errors.

### Pitfall 5: OFT Normalization Prefactors
**What goes wrong:** `oft!` returns subnormalized output ("multiply by sqrt(1 / sigma sqrt(2 * pi))"). `time_oft!` also returns subnormalized output (multiply by t0 * sqrt(sigma * sqrt(2/pi) / (2*pi))). These prefactors differ.
**Why it happens:** Different normalization conventions between energy-domain and time-domain OFT.
**How to avoid:** Apply the correct prefactor for each domain before comparing. See `time_tests.jl` for the exact prefactors:
- Energy: `energy_oft_prefactor = 1 / sqrt(sigma * sqrt(2 * pi))`
- Time/Trotter: `time_oft_prefactor = t0 * sqrt(sigma * sqrt(2 / pi) / (2 * pi))`
**Warning signs:** OFT outputs differ by a factor rather than small errors.

### Pitfall 6: Step Error Scaling Test Needs Deterministic DM Evolution
**What goes wrong:** `run_thermalization` uses random jump selection (picks one jump per step), introducing stochastic noise that obscures O(delta) scaling.
**Why it happens:** The thermalization function implements stochastic trajectory mixing, not deterministic DM channel application.
**How to avoid:** For DM step error scaling tests, use the Liouvillian matrix directly. Compute: `rho(t+delta) = rho(t) + delta * reshape(L * vec(rho(t)), dim, dim)` (Euler step). Or use the matrix exponential: `rho(t) = reshape(exp(t*L) * vec(rho(0)), dim, dim)` for exact reference.
**Warning signs:** Error scaling shows random fluctuations instead of clean power-law behavior.

### Pitfall 7: Aqua Stale Deps False Positive on Pkg
**What goes wrong:** Aqua flags `Pkg` as a stale dependency.
**Why it happens:** `Pkg` is used in `QuantumFurnace.jl` (`using Pkg`) but Aqua may not detect indirect usage via Pkg.project().
**How to avoid:** Use `stale_deps = (ignore=[:Pkg],)` in `Aqua.test_all`.
**Warning signs:** Aqua test failure on stale_deps with `Pkg` listed.

## Code Examples

### Example 1: Detailed Balance Test (DMTST-01)
```julia
# Source: codebase analysis of construct_lindbladian + run_lindbladian
@testset "DMTST-01: Bohr detailed balance" begin
    # Use 3-qubit system per success criteria
    ham3 = make_small_test_system()  # 3-qubit fixture
    jumps3 = ham3.jumps
    gibbs3 = ham3.gibbs

    config = LiouvConfig(
        num_qubits = 3,
        with_coherent = true,
        with_linear_combination = true,
        domain = BohrDomain(),
        beta = BETA,
        sigma = SIGMA,
        a = BETA / 30.0,
        b = 0.4,
        num_energy_bits = NUM_ENERGY_BITS,
        w0 = W0,
        t0 = T0,
        num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
    )

    liouv = construct_lindbladian(jumps3, config, ham3.hamiltonian)
    eig = eigen(liouv)
    ss_idx = argmin(abs.(real.(eig.values)))
    ss_vec = eig.vectors[:, ss_idx]
    dim3 = 2^3
    ss_dm = reshape(ss_vec, dim3, dim3)
    ss_dm = (ss_dm + ss_dm') / 2
    ss_dm ./= tr(ss_dm)

    dist = trace_distance_h(Hermitian(ss_dm), gibbs3)
    @test dist < 1e-10
end
```

### Example 2: Domain Hierarchy Test (DMTST-02)
```julia
# Source: codebase analysis of construct_lindbladian dispatch pattern
@testset "DMTST-02: Domain error hierarchy" begin
    distances = Dict{Symbol, Float64}()

    for (name, domain) in [(:bohr, BohrDomain()), (:energy, EnergyDomain()),
                            (:time, TimeDomain()), (:trotter, TrotterDomain())]
        config = make_liouv_config(domain)
        trotter_obj = (domain isa TrotterDomain) ? TEST_TROTTER : nothing
        liouv = construct_lindbladian(TEST_JUMPS, config, TEST_HAM; trotter=trotter_obj)

        eig = eigen(liouv)
        ss_idx = argmin(abs.(real.(eig.values)))
        ss_vec = eig.vectors[:, ss_idx]
        ss_dm = reshape(ss_vec, DIM, DIM)
        ss_dm = (ss_dm + ss_dm') / 2
        ss_dm ./= tr(ss_dm)

        distances[name] = trace_distance_h(Hermitian(ss_dm), TEST_GIBBS)
    end

    @test distances[:bohr] <= distances[:energy]
    @test distances[:energy] <= distances[:time]
    @test distances[:time] <= distances[:trotter]
end
```

### Example 3: OFT Consistency (DMTST-06)
```julia
# Source: time_tests.jl pattern + ofts.jl API
energy_oft_prefactor = 1 / sqrt(SIGMA * sqrt(2 * pi))
time_oft_prefactor = T0 * sqrt(SIGMA * sqrt(2 / pi) / (2 * pi))

w = -0.12
A = TEST_JUMPS[1]  # X on site 1
caches = OFTCaches(DIM)

# Energy OFT
A_energy = Matrix{ComplexF64}(undef, DIM, DIM)
oft!(A_energy, A, w, TEST_HAM, SIGMA)
A_energy .*= energy_oft_prefactor

# Time OFT (need time_labels)
energy_labels = create_energy_labels(NUM_ENERGY_BITS, W0)
time_labels_full = energy_labels .* (T0 / W0)
oft_time_labels = truncate_time_labels_for_oft(time_labels_full, SIGMA)

A_time = Matrix{ComplexF64}(undef, DIM, DIM)
time_oft!(A_time, caches, A, w, TEST_HAM, oft_time_labels, SIGMA)
A_time .*= time_oft_prefactor

@test norm(A_energy - A_time) < TOL_QUADRATURE

# Trotter OFT
A_trott = Matrix{ComplexF64}(undef, DIM, DIM)
trotter_oft!(A_trott, caches, A, w, TEST_TROTTER, oft_time_labels, SIGMA)
A_trott .*= time_oft_prefactor

# Trotter error should be larger than time quadrature error
@test norm(A_energy - A_trott) >= norm(A_energy - A_time) - 1e-10  # small tolerance for numerics
```

### Example 4: Aqua.jl Package Quality (TINF-03)
```julia
using Aqua
using QuantumFurnace

@testset "Aqua.jl package quality" begin
    Aqua.test_all(QuantumFurnace;
        ambiguities = false,     # Disable: multiple dispatch may create legitimate ambiguities
        deps_compat = (ignore=[:Pkg, :Printf, :Random, :LinearAlgebra, :SparseArrays,
                               :SharedArrays, :Distributed, :Base],),
        stale_deps = (ignore=[:Pkg],),
        piracies = false,        # Disable: kron! on AbstractMatrix may be flagged
    )
end
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Cholesky for PSD guard | Eigendecomposition (Phase 2) | 2026-02-14 | Silent clamp of negative eigenvalues |
| create_trotter() (exported, undefined) | TrottTrott() constructor directly | Phase 2 discovery | Use constructor, not the non-existent helper |
| B_test.jl informal playground test | Formal test_dm_scaling.jl with assertions | Phase 3 (planned) | Reproducible CI-compatible coherent term tests |

**Deprecated/outdated:**
- `create_trotter()`: exported but never defined; use `TrottTrott(hamiltonian, t0, num_trotter_steps)` directly
- Old `time_tests.jl`, `B_test.jl`, `trajectory_test.jl`: playground files (use Revise, hardcoded params), not proper tests

## Open Questions

1. **Performance of 4-qubit Lindbladian eigendecomposition**
   - What we know: dim^2 = 256 matrix; `eigen` should handle this in milliseconds
   - What's unclear: Whether constructing the Lindbladian for TimeDomain/TrotterDomain with NUM_ENERGY_BITS=12 is fast enough for CI (involves NUFFT precomputation)
   - Recommendation: Time it during implementation. If > 30s for one config, reduce NUM_ENERGY_BITS for DM tests or use only 3-qubit for slow domains.

2. **Exact tolerance for hierarchy test strict inequalities**
   - What we know: dist_bohr should be O(1e-13), dist_energy larger, etc.
   - What's unclear: Whether Energy/Time distances are close enough that numerical noise could violate strict ordering
   - Recommendation: Use `<=` (not `<`) for the hierarchy checks. If Energy and Time are numerically equal, that's still a valid hierarchy. Add a margin test: `@test distances[:energy] - distances[:bohr] >= -1e-12` (non-negative difference with tolerance).

3. **DM step error test initial conditions**
   - What we know: Starting from Gibbs state, single-step error should be zero (Gibbs is fixed point of exact channel). Starting from non-Gibbs state, error scaling with delta is what we want.
   - What's unclear: Best initial state for clean scaling measurement
   - Recommendation: Use maximally mixed state (rho = I/dim) as initial. It has known nonzero distance to Gibbs and produces clean scaling.

## Sources

### Primary (HIGH confidence)
- Codebase direct analysis: `src/bohr_domain.jl`, `src/coherent.jl`, `src/ofts.jl`, `src/energy_domain.jl`, `src/trotter_domain.jl`, `src/furnace.jl`, `src/furnace_utensils.jl`, `src/jump_workers.jl`, `src/trajectories.jl`, `src/qi_tools.jl`, `src/structs.jl`, `src/kraus.jl`
- Existing test infrastructure: `test/test_helpers.jl`, `test/test_cptp.jl`, `test/test_trajectory_fixes.jl`, `test/runtests.jl`
- Existing playground tests: `test/B_test.jl`, `test/time_tests.jl`, `test/trajectory_test.jl` (patterns for B comparison and OFT verification)
- Phase 2 research and plans: `.planning/phases/02-trajectory-bug-fixes/02-RESEARCH.md`

### Secondary (MEDIUM confidence)
- [Aqua.jl GitHub](https://github.com/JuliaTesting/Aqua.jl) - API for test_all and individual tests
- [Aqua.jl docs](https://juliatesting.github.io/Aqua.jl/dev/) - deps_compat options and ignore patterns

### Tertiary (LOW confidence)
- None

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - all libraries already in Project.toml, no new deps needed
- Architecture: HIGH - all test patterns derived from direct codebase analysis, existing playground files demonstrate exact calling patterns
- Pitfalls: HIGH - identified from analyzing Phase 2 issues, codebase normalization conventions, and existing test infrastructure
- Code examples: HIGH - patterns verified against existing working code in B_test.jl, time_tests.jl, and Phase 2 tests

**Research date:** 2026-02-14
**Valid until:** 2026-03-14 (stable domain, no external dependency changes expected)
