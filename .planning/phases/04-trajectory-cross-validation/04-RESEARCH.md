# Phase 4: Trajectory Cross-Validation - Research

**Researched:** 2026-02-14
**Domain:** Trajectory-averaged density matrix cross-validation against DM ground truth
**Confidence:** HIGH

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- Fixed threshold comparison: trace distance between trajectory-averaged rho and DM rho must be < 0.01
- 10,000 trajectories per test point
- Fixed RNG seed via StableRNG for deterministic, reproducible results
- No empirical C/sqrt(N) bound fitting -- that belongs in Phase 5
- 3-qubit Heisenberg system only (no 4-qubit)
- Reuse make_test_system() fixture from Phase 3 (same Hamiltonian, beta, sigma)
- Initial state: normalized all-ones pure state |psi> = (1,...,1)/sqrt(d), with DM starting from rho = |psi><psi|
- Single-step tests: all three domains (Energy, Time, Trotter) with a delta sweep to verify delta^2 scaling of trajectory-vs-DM error
- Multi-step convergence: TrotterDomain only, single delta value, with_coherent=true
- Both DM and trajectory-averaged rho must reach trace distance < 1e-3 to Gibbs state
- TrotterDomain with_coherent=true only (no Bohr baseline, no comparison with with_coherent=false)
- Run until 1e-3 threshold is reached for the given delta (not a fixed step count)
- Claude picks the delta value that achieves convergence practically
- No delta sweep for multi-step tests -- single delta value chosen by Claude
- Single-threaded execution (multi-threading deferred to a future milestone)
- Tests in a separate test group, not in main Pkg.test() suite
- Prioritize faster runtime over exhaustive parameter coverage

### Claude's Discretion
- Exact delta value for multi-step convergence test
- Delta values for the single-step delta sweep
- Number of multi-step thermalization steps (run until converged)
- Test file organization within the separate test group
- How to structure the delta^2 scaling check (number of delta values, regression method)

### Deferred Ideas (OUT OF SCOPE)
- Multi-threaded trajectory execution -- future milestone
- 4-qubit system tests -- not needed for current validation
- Empirical C/sqrt(N) trajectory scaling verification -- Phase 5
</user_constraints>

## Summary

This phase validates that trajectory-averaged density matrices match DM (density matrix) evolution across the three approximation domains (Energy, Time, Trotter). The core test is: run 10,000 trajectories each starting from the same pure state, average the outer products to form rho_traj, then compare against the DM rho produced by stepping the same initial density matrix through the CPTP channel. The comparison measures trace distance, which must be < 0.01 for statistical matching.

There are two distinct test types: (1) single-step delta-sweep tests that verify the trajectory-vs-DM error scales as O(delta^2) for all three domains, and (2) a multi-step convergence test for TrotterDomain with coherent term that verifies both DM and trajectory modes converge to within 1e-3 of the Gibbs state. All tests use a 3-qubit Heisenberg system with fixed StableRNG seeds.

A critical architectural difference in the codebase shapes the implementation: the trajectory mode (`run_trajectories` / `build_trajectoryframework`) sums ALL jump operators into a single R matrix and a single total B coherent term, applying the combined channel in each step. The DM mode (`run_thermalization`) randomly selects ONE jump per step and applies that single jump's channel (with rescaled delta). For cross-validation, we must NOT use `run_thermalization` as the DM reference, since it introduces additional randomness from jump selection. Instead, we must apply the all-jumps-at-once CPTP channel directly to the density matrix, which is what `jump_contribution!` does when applied with all jumps summed. Specifically, we construct the Liouvillian (L) and compute the exact DM step as rho -> rho + delta * L(rho), or equivalently apply the full CPTP map matching what the trajectory channel implements.

**Primary recommendation:** Build the DM reference by constructing the Liouvillian matrix L and applying it as a matrix exponential or Euler step to the initial rho = |psi><psi|. For single-step tests, use the exact matrix exponential exp(delta * L) as ground truth. For multi-step convergence, iterate the Euler step rho += delta * reshape(L * vec(rho), d, d). This avoids the random jump selection in `run_thermalization` and gives a deterministic DM reference that corresponds exactly to the all-jumps-combined Lindbladian.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Test | stdlib | Julia test framework | Used throughout test suite |
| LinearAlgebra | stdlib | Matrix operations, eigen, exp, norm | All DM and trajectory operations |
| Random | stdlib | Global RNG seeding for trajectory reproducibility | `step_along_trajectory!` uses `rand()` on global RNG |
| StableRNGs | >= 1.0 | Deterministic, version-stable RNG seeds | Already in [extras]/test target; `rand()` seeded via `Random.seed!(StableRNG(seed))` |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| BSON | deps | Load 3-qubit Hamiltonian | Via `make_small_test_system()` in test_helpers.jl |
| Printf | stdlib | Diagnostic output formatting | Already used in codebase for @printf |

### Already Available (No Installation Needed)
All required libraries are already in Project.toml. StableRNGs is in `[extras]` and `[targets].test`. No new dependencies needed.

## Architecture Patterns

### Critical: DM vs Trajectory Channel Mismatch

The most important architectural detail for this phase is that the DM and trajectory modes apply the CPTP channel differently:

**Trajectory mode** (`build_trajectoryframework` + `step_along_trajectory!`):
- Sums ALL jumps into a single R = sum_k sum_w gamma(w) * A_kw^dag A_kw
- Builds a single total B = sum_k B_k for the coherent term
- Each step applies: U_B first, then sample from {K0, K_res, K_{k,w}} for ALL (k,w) simultaneously
- This is the "combined Lindbladian" channel

**DM mode** (`run_thermalization`):
- Randomly picks ONE jump index k per step
- Applies that single jump k's coherent unitary U_{B_k} and dissipative channel
- Rescales delta by 1/p_jump = num_jumps to compensate for random selection
- Over many steps, this averages to the same Lindbladian, but each step differs

**Consequence for cross-validation:** We must NOT use `run_thermalization` as DM reference for single-step tests. Instead, use the Liouvillian matrix directly:

```julia
L = construct_lindbladian(jumps, liouv_config, hamiltonian; trotter=trotter)
rho_vec = vec(rho0)
# Exact single step:
rho_exact = reshape(exp(delta * L) * rho_vec, dim, dim)
# Or Euler step for multi-step:
rho_euler = rho0 + delta * reshape(L * rho_vec, dim, dim)
```

The Liouvillian L includes contributions from ALL jumps (including coherent terms B_k), so it matches what the trajectory channel implements.

### Pattern 1: Single-Step Trajectory-vs-DM Comparison

**What:** Run ntraj=10000 single trajectory steps from |psi>, accumulate rho_traj = (1/N) sum |psi_i><psi_i|. Compare against the exact DM step exp(delta * L) applied to rho0 = |psi><psi|.

**How to implement:**

```julia
# 1. Build the all-jumps Liouvillian
liouv_config = LiouvConfig(
    num_qubits = 3, with_coherent = false,  # or true for coherent tests
    with_linear_combination = true, domain = EnergyDomain(),
    beta = BETA, sigma = SIGMA, a = BETA/30.0, b = 0.4,
    num_energy_bits = NUM_ENERGY_BITS, w0 = W0, t0 = T0,
    num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
)
L = construct_lindbladian(SMALL_JUMPS, liouv_config, SMALL_HAM)

# 2. Compute exact DM reference
dim = SMALL_DIM  # 8
psi0 = fill(ComplexF64(1.0), dim) / sqrt(dim)
rho0 = psi0 * psi0'
rho_dm = reshape(exp(delta * L) * vec(rho0), dim, dim)
rho_dm = (rho_dm + rho_dm') / 2  # Hermitianize

# 3. Build trajectory framework and run N trajectories
therm_config = ThermalizeConfig(
    num_qubits = 3, with_coherent = false,
    with_linear_combination = true, domain = EnergyDomain(),
    beta = BETA, sigma = SIGMA, a = BETA/30.0, b = 0.4,
    num_energy_bits = NUM_ENERGY_BITS, w0 = W0, t0 = T0,
    num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
    mixing_time = delta, delta = delta,
)
precomputed = precompute_data(EnergyDomain(), therm_config, SMALL_HAM)
scratch = KrausScratch(ComplexF64, dim)
fw = build_trajectoryframework(SMALL_JUMPS, SMALL_HAM, therm_config, precomputed, scratch, delta)

rho_traj = zeros(ComplexF64, dim, dim)
using Random
Random.seed!(StableRNG(42))  # or: copy!(Random.default_rng(), StableRNG(42))
for _ in 1:10000
    psi = copy(psi0)
    step_along_trajectory!(psi, fw)
    rho_traj .+= psi * psi'
end
rho_traj ./= 10000
rho_traj = (rho_traj + rho_traj') / 2

dist = QuantumFurnace.trace_distance_h(Hermitian(rho_dm), Hermitian(rho_traj))
@test dist < 0.01
```

**Important:** The RNG seeding is tricky. `step_along_trajectory!` calls `rand()` on the global RNG. Phase 2 discovered that `StableRNG` cannot be directly passed to the function. The correct approach is to seed the global RNG task-local state. The cleanest way:

```julia
using StableRNGs
rng = StableRNG(42)
# Copy the StableRNG state into the global RNG
# Actually: Random.seed! only accepts integers, not StableRNG objects
# Must use: copy!(Random.default_rng(), rng)  -- but this may not work with StableRNG
# Safest: use Random.seed!(42) for reproducibility (not version-stable, but deterministic within a session)
```

**RNG resolution:** Given that `step_along_trajectory!` uses `rand()` (global RNG) without accepting an rng keyword, and `copy!(Random.default_rng(), StableRNG(...))` is not supported across RNG types, the practical approach is to use `Random.seed!(FIXED_SEED)` at the start of each test. This gives session-local reproducibility. For true cross-version reproducibility, the trajectory code would need to accept an rng parameter (out of scope for this phase).

### Pattern 2: Delta-Squared Scaling Verification

**What:** For multiple delta values, compute single-step trajectory-vs-DM trace distance and verify the error scales as O(delta^2).

**Recommended delta values:** [0.2, 0.1, 0.05, 0.025]
- These span one decade of delta
- At delta=0.2, the error should be detectable (large enough for 10000 trajectories to resolve)
- At delta=0.025, the error should be small (verifying scaling, not just noise floor)

**Scaling check method:** Compute log-log regression of (delta, error) pairs. For O(delta^2), the slope should be approximately 2.0. Accept slope in range [1.5, 2.5] to account for statistical noise from finite trajectory count.

```julia
deltas = [0.2, 0.1, 0.05, 0.025]
errors = Float64[]
for delta in deltas
    # ... compute rho_dm and rho_traj for this delta ...
    push!(errors, trace_distance_h(Hermitian(rho_dm), Hermitian(rho_traj)))
end

# Check consecutive ratios: halving delta should quarter the error
ratios = [errors[i] / errors[i+1] for i in 1:length(errors)-1]
for ratio in ratios
    @test 2.0 <= ratio <= 6.0  # Expect ~4 for O(delta^2), allow statistical noise
end
```

**Alternative (more robust):** Use log-log slope:
```julia
log_deltas = log.(deltas)
log_errors = log.(errors)
# Linear regression: log_error = slope * log_delta + intercept
n = length(deltas)
slope = (n * sum(log_deltas .* log_errors) - sum(log_deltas) * sum(log_errors)) /
        (n * sum(log_deltas.^2) - sum(log_deltas)^2)
@test 1.5 <= slope <= 2.5
```

### Pattern 3: Multi-Step Convergence Test (TrotterDomain + Coherent)

**What:** Run multiple DM Euler steps and trajectory batches until both converge to within 1e-3 of the Gibbs state.

**Recommended delta:** 0.01
- Small enough for DM Euler to converge (O(delta) accumulated error stays small)
- Large enough that each step makes meaningful progress toward Gibbs
- At 3 qubits (dim=8), the DM Euler step with construct_lindbladian is fast

**Implementation pattern:**

```julia
# Build Liouvillian for TrotterDomain + coherent
liouv_config = LiouvConfig(
    num_qubits = 3, with_coherent = true,
    ..., domain = TrotterDomain(), ...
)
L = construct_lindbladian(SMALL_JUMPS, liouv_config, SMALL_HAM; trotter=SMALL_TROTTER)

# DM evolution via Euler
dim = SMALL_DIM
psi0 = fill(ComplexF64(1.0), dim) / sqrt(dim)
rho_dm = psi0 * psi0'
delta = 0.01

# Gibbs state in Trotter eigenbasis (for TrotterDomain comparison)
gibbs_comp = SMALL_HAM.eigvecs * SMALL_GIBBS * SMALL_HAM.eigvecs'
gibbs_trott = Hermitian(SMALL_TROTTER.eigvecs' * gibbs_comp * SMALL_TROTTER.eigvecs)

max_steps = 10000  # Safety cap
for step in 1:max_steps
    rho_dm .+= delta .* reshape(L * vec(rho_dm), dim, dim)
    rho_dm = (rho_dm + rho_dm') / 2  # Hermitianize
    dist = trace_distance_h(Hermitian(rho_dm), gibbs_trott)
    if dist < 1e-3
        break
    end
end

# Also run trajectory average for the same number of steps
# Compare rho_traj to rho_dm at the end
```

**Important TrotterDomain basis considerations:**
- The Liouvillian operates in the Trotter eigenbasis when domain = TrotterDomain
- The initial state psi0 = (1,...,1)/sqrt(d) is in the COMPUTATIONAL basis
- For DM: rho0 must be transformed to Trotter eigenbasis: `rho0_trott = SMALL_TROTTER.eigvecs' * rho0_comp * SMALL_TROTTER.eigvecs`
- For trajectories: psi0 must be transformed to Trotter eigenbasis: `psi0_trott = SMALL_TROTTER.eigvecs' * psi0_comp`
- Wait -- looking at `build_trajectoryframework` more carefully: it takes `ham_or_trott` and internally uses `jump.in_eigenbasis`. For TrotterDomain, `ham_or_trott` is TrottTrott. The trajectory framework's R is built using NUFFT prefactors from `ham_or_trott.bohr_freqs` which for TrottTrott are the Trotter bohr frequencies. The jump's `in_eigenbasis` is in the HAMILTONIAN eigenbasis (as created in test_helpers).

**Critical basis detail:** The `SMALL_JUMPS` have `in_eigenbasis` in Hamiltonian eigenbasis. For TrotterDomain, both `construct_lindbladian` and `build_trajectoryframework` use these jumps with the Trotter's bohr_freqs. The NUFFT prefactors are computed from `ham_or_trott.bohr_freqs` (Trotter Bohr frequencies). The jump OFT `A_{a,w} = A_a .* prefactor(w)` uses `jump.in_eigenbasis` (H-eigenbasis) elementwise-multiplied by prefactors computed from Trotter Bohr frequencies. This mixing of bases is intentional and consistent between DM and trajectory modes.

However, the coherent B_trotter function does transform the jump to Trotter eigenbasis internally:
```julia
U = trotter.trafo_from_eigen_to_trotter
jump_in_trotter = U * jump.in_eigenbasis * U'
```

So the B_trotter result is in Trotter eigenbasis, and U_B = exp(-i delta B) is applied in whatever basis the state lives in. The entire framework works in the H-eigenbasis with Trotter Bohr frequencies providing the spectral information.

**Where does the Liouvillian L live for TrotterDomain?** Looking at `construct_lindbladian`: it creates a dim^2 x dim^2 matrix. The jumps are in H-eigenbasis, the NUFFT prefactors use Trotter Bohr freqs, and B_trotter is in Trotter eigenbasis. The vectorization via `vectorize_liouv_diss_and_add!` uses the jump OFT (H-eigenbasis with Trotter prefactors). The coherent part uses B in Trotter eigenbasis.

Actually wait -- re-reading `construct_lindbladian` more carefully:
```julia
Btot = precompute_coherent_total_B(jumps, ham_or_trott, config, precomputed_data)
if Btot !== nothing
    vectorize_liouvillian_coherent!(total_lindbladian, Btot, ws)
end
```
And `precompute_coherent_total_B` for TrotterDomain calls `B_trotter(jumps, ham_or_trott, ...)` which returns B in Trotter eigenbasis. Then the dissipative part: `jump_contribution!(total_lindbladian, TrotterDomain(), jump, ham_or_trott, ...)` uses NUFFT prefactors with `jump.in_eigenbasis` in H-eigenbasis.

This means the Liouvillian L mixes two bases. But this is consistent: the B term is added in Trotter eigenbasis, the dissipative part uses jumps in H-eigenbasis with Trotter-Bohr-frequency prefactors. The density matrix that this L acts on has a specific basis convention. Looking at `run_thermalization`:
```julia
if config.domain isa TrotterDomain
    ham_or_trott = trotter
    gibbs = Hermitian(trotter.eigvecs' * hamiltonian.eigvecs * hamiltonian.gibbs *
                      hamiltonian.eigvecs' * trotter.eigvecs)
```
It transforms the Gibbs state to Trotter eigenbasis. But the initial `evolving_dm` is passed in by the caller. The Liouvillian for TrotterDomain operates on density matrices in the "mixed" basis (H-eigenbasis for dissipative, Trotter-eigenbasis for coherent).

**Actually -- re-reading the jump_contribution! code more carefully:** For Time/TrotterDomain, the dissipative part computes `@. jump_oft = jump.in_eigenbasis * nufft_prefactor_matrix`. The NUFFT prefactors encode the OFT from Trotter/Time domain. The result `jump_oft` is an operator in whatever basis `jump.in_eigenbasis` is in (H-eigenbasis), modulated by Trotter Bohr-frequency-based prefactors. The `rho_jump += delta * jump_oft * rho * jump_oft'` operates in the same basis as rho.

So for TrotterDomain, BOTH the Liouvillian construction and the trajectory framework operate consistently in the H-eigenbasis for the state/operators, using Trotter Bohr frequencies for the spectral modulation. The B_trotter term transforms the jump to Trotter eigenbasis internally to compute B, but the result B is in Trotter eigenbasis, which is then used to rotate rho in Trotter eigenbasis...

Let me resolve this more carefully. In `construct_lindbladian`, the B_trotter is computed and added to L via `vectorize_liouvillian_coherent!(L, Btot, ws)`. This adds `-i * kron(B, I) + i * kron(I, B^T)` to L. So L acts on vec(rho) and adds `-i[B, rho]`. If B is in Trotter eigenbasis, then rho must also be in Trotter eigenbasis for this to be correct. BUT the dissipative part uses jump.in_eigenbasis (H-eigenbasis). This seems like a basis inconsistency.

**Resolution from Phase 3 DMTST-02 test:** Looking at the existing test_dm_detailed_balance.jl (DMTST-02), the Gibbs state for TrotterDomain is computed as:
```julia
gibbs_ref = if domain isa TrotterDomain
    gibbs_comp = TEST_HAM.eigvecs * TEST_GIBBS * TEST_HAM.eigvecs'
    Hermitian(TEST_TROTTER.eigvecs' * gibbs_comp * TEST_TROTTER.eigvecs)
```
This transforms Gibbs to Trotter eigenbasis. And the Liouvillian's fixed point is found and compared in this basis. So the Liouvillian for TrotterDomain operates in some effective basis, and the fixed point comparison uses Gibbs in Trotter eigenbasis. This suggests the Liouvillian L for TrotterDomain indeed operates on states in a basis that is consistent with Trotter eigenbasis for the coherent term.

**But wait -- in the DMTST-02 test, the fixed point ss_dm is found by eigendecomposing L. The comparison is:**
```julia
distances[name] = QuantumFurnace.trace_distance_h(Hermitian(ss_dm), gibbs_ref)
```
And the test PASSES (the hierarchy is verified). So the Liouvillian for TrotterDomain is self-consistent. The jumps are in H-eigenbasis, the B_trotter is in Trotter eigenbasis, and the resulting L operates on vec(rho) in a consistent way -- because the NUFFT prefactors effectively encode a change-of-basis through the spectral modulation.

**Practical impact for Phase 4:** For the DM reference, we construct L = construct_lindbladian(jumps, config, ham; trotter=trotter), and the initial state rho0 must be in the same basis that L operates on. Since we know from DMTST-02 that comparing with Gibbs in Trotter eigenbasis works, we should:
1. Transform the initial psi0 to the basis L operates on (this is the H-eigenbasis for dissipative part / mixed for coherent)
2. Actually, simplest approach: just start from psi0 in computational basis, transform to H-eigenbasis via `psi0_eigen = ham.eigvecs' * psi0_comp`, form rho0 = psi0_eigen * psi0_eigen', and evolve with L.

For comparing with Gibbs: use SMALL_GIBBS (already in H-eigenbasis) for Energy/Time domains, and the Trotter-transformed Gibbs for TrotterDomain.

**Even simpler approach:** For single-step comparison, we don't need to compare with Gibbs at all. We compare rho_traj with rho_dm. Both start from the same initial state and undergo the same channel. The basis is irrelevant as long as both use the same basis consistently.

For EnergyDomain and TimeDomain, use SMALL_HAM as ham_or_trott. psi0 in H-eigenbasis, L constructed with SMALL_HAM.
For TrotterDomain, use SMALL_TROTTER as ham_or_trott. The trajectory framework is built with SMALL_TROTTER. The L is constructed with SMALL_HAM + trotter=SMALL_TROTTER.

**For the trajectory framework:** `build_trajectoryframework(SMALL_JUMPS, ham_or_trott, config, precomputed, scratch, delta)` where `ham_or_trott = SMALL_TROTTER` for TrotterDomain. The framework operates on psi vectors. These must be in H-eigenbasis (since jump.in_eigenbasis is in H-eigenbasis and the framework applies jumps in that space).

So for all domains: start psi0 in H-eigenbasis, form rho0 in H-eigenbasis, construct L and apply to vec(rho0), run trajectories from psi0 in H-eigenbasis.

For TrotterDomain Gibbs comparison: transform Gibbs to same basis that L's fixed point lives in (Trotter eigenbasis per DMTST-02), OR simply don't compare against Gibbs for the DM reference -- compare rho_traj against rho_dm_stepped directly.

### Pattern 4: Trajectory Test File Organization

**Recommended structure:** Tests go in a separate directory and are NOT included in runtests.jl (per user decision: "separate test group, not in main Pkg.test() suite").

```
test/
  runtests.jl                          # Existing, unchanged
  test_helpers.jl                      # Add SMALL_TROTTER fixture
  trajectory_validation/
    run_trajectory_tests.jl            # Entry point: include helpers, run tests
    test_single_step_crossval.jl       # TVAL-02/03/04: single-step delta sweep
    test_coherent_convergence.jl       # TVAL-06: multi-step Gibbs convergence
```

**Running:** `julia --project test/trajectory_validation/run_trajectory_tests.jl`

### Anti-Patterns to Avoid

- **Using `run_thermalization` as DM reference:** This randomly selects one jump per step, introducing stochastic noise that obscures the trajectory-vs-DM comparison. Use Liouvillian matrix directly.
- **Using `run_trajectories` as-is without controlling RNG:** The function calls `_evolve_along_trajectory!` which calls `step_along_trajectory!` which calls `rand()`. Must seed global RNG before each trajectory batch.
- **Forgetting to Hermitianize rho_traj:** Trajectory averaging produces rho_traj = (1/N) sum |psi><psi|, which should be Hermitian by construction. But floating point accumulation can break this. Always Hermitianize before trace_distance_h.
- **Not creating SMALL_TROTTER fixture:** test_helpers.jl has SMALL_SYSTEM/SMALL_HAM/SMALL_JUMPS but no SMALL_TROTTER. Must add.
- **Using 4-qubit NUM_QUBITS in config:** make_thermalize_config uses NUM_QUBITS=4. Must create 3-qubit configs manually or add a make_small_thermalize_config helper.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Liouvillian construction | Custom L matrix builder | `construct_lindbladian(jumps, config, ham; trotter=trotter)` | Handles all domains, normalization, coherent terms consistently |
| Trajectory framework | Custom Kraus map builder | `build_trajectoryframework(jumps, ham_or_trott, config, precomputed, scratch, delta)` | Already implements Chen's channel correctly (verified Phase 2) |
| Trace distance | Manual eigenvalue computation | `QuantumFurnace.trace_distance_h(Hermitian(rho1), Hermitian(rho2))` | Non-exported but accessible, handles Hermitian inputs |
| Gibbs state | Manual Boltzmann weights | `SMALL_GIBBS` from test_helpers.jl (precomputed) | Already in H-eigenbasis, normalized, trace 1 |
| CPTP channel application | Custom Kraus iteration | `step_along_trajectory!(psi, fw)` for trajectory, `exp(delta * L) * vec(rho)` for DM | Both are verified correct implementations |

**Key insight:** The trajectory framework is verified CPTP (Phase 2), the Liouvillian is verified to have correct fixed points (Phase 3). Phase 4 just needs to run both and compare.

## Common Pitfalls

### Pitfall 1: RNG Seeding for Trajectory Reproducibility
**What goes wrong:** `step_along_trajectory!` uses `rand()` on the global RNG. Passing a StableRNG to the function is not supported.
**Why it happens:** The trajectory code was not designed with pluggable RNG. It uses the task-local default RNG.
**How to avoid:** Use `Random.seed!(FIXED_SEED)` before each test. This seeds the global MersenneTwister. For cross-session reproducibility, the exact RNG output may differ between Julia versions, but within a session, results are deterministic.
**Warning signs:** Running the same test twice gives different results. (Should not happen with seed.)

### Pitfall 2: DM Reference Must Match Trajectory Channel
**What goes wrong:** Using `run_thermalization` as DM reference produces rho_dm that differs from trajectory average even for N -> infinity.
**Why it happens:** `run_thermalization` applies one random jump per step (with rescaled delta), while trajectory applies all jumps simultaneously. Over many steps they converge to the same fixed point, but single-step DM output differs from single-step trajectory average.
**How to avoid:** Use the Liouvillian matrix L. `rho_dm = reshape(exp(delta * L) * vec(rho0), dim, dim)` gives the exact single-step DM that the trajectory channel should match.
**Warning signs:** Even with 100,000 trajectories, single-step trace distance remains > 0.1.

### Pitfall 3: TrotterDomain Requires SMALL_TROTTER Fixture
**What goes wrong:** TrotterDomain tests fail because there's no 3-qubit Trotter object.
**Why it happens:** test_helpers.jl only defines `TEST_TROTTER` for the 4-qubit system.
**How to avoid:** Add to test_helpers.jl:
```julia
function make_small_test_trotter()
    TrottTrott(SMALL_HAM, T0, NUM_TROTTER_STEPS_PER_T0)
end
const SMALL_TROTTER = make_small_test_trotter()
```
**Warning signs:** "TrotterDomain requires `trotter`" error in tests.

### Pitfall 4: Gibbs State Basis for TrotterDomain
**What goes wrong:** Comparing DM-stepped rho against Gibbs state uses wrong basis.
**Why it happens:** SMALL_GIBBS is in H-eigenbasis. TrotterDomain Liouvillian has its fixed point in a basis that corresponds to Trotter eigenbasis for comparison.
**How to avoid:** For the multi-step Gibbs convergence test, transform Gibbs to Trotter eigenbasis:
```julia
gibbs_comp = SMALL_HAM.eigvecs * SMALL_GIBBS * SMALL_HAM.eigvecs'
gibbs_trott = Hermitian(SMALL_TROTTER.eigvecs' * gibbs_comp * SMALL_TROTTER.eigvecs)
```
This pattern is already used in DMTST-02 and verified to work.
**Warning signs:** DM evolution "converges" but trace distance to Gibbs stays large (~0.1).

### Pitfall 5: Config Creation for 3-Qubit System
**What goes wrong:** Using `make_thermalize_config` creates configs with `num_qubits=4`.
**Why it happens:** The helper is hardcoded to `NUM_QUBITS=4`.
**How to avoid:** Create configs manually with `num_qubits=3`:
```julia
config = ThermalizeConfig(
    num_qubits = 3,
    with_coherent = with_coherent,
    with_linear_combination = true,
    domain = domain,
    beta = BETA, sigma = SIGMA,
    a = BETA / 30.0, b = 0.4,
    num_energy_bits = NUM_ENERGY_BITS,
    w0 = W0, t0 = T0,
    num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
    mixing_time = mixing_time,
    delta = delta,
)
```
Or add a `make_small_thermalize_config` helper.
**Warning signs:** Config validation fails or dimension mismatches.

### Pitfall 6: Matrix Exponential Performance for Multi-Step
**What goes wrong:** Using `exp(delta * L)` for every step in a multi-step test is expensive (dim^2 x dim^2 = 64 x 64 matrix exp at each step).
**Why it happens:** Matrix exponential is O(n^3) per call. For 10000 steps, that's 10000 * O(64^3).
**How to avoid:** For multi-step tests, compute `exp_L = exp(delta * L)` ONCE, then apply repeatedly: `rho_vec = exp_L * rho_vec`. This is just a matrix-vector multiply per step.
**Warning signs:** Multi-step test takes > 5 minutes.

### Pitfall 7: Trajectory Accumulation Numerical Precision
**What goes wrong:** After 10,000 trajectory accumulations, rho_traj may lose Hermiticity or have trace != 1.
**Why it happens:** Floating point accumulation errors from 10000 rank-1 additions.
**How to avoid:** Hermitianize and normalize after accumulation:
```julia
rho_traj = (rho_traj + rho_traj') / 2
rho_traj ./= tr(rho_traj)
```
**Warning signs:** trace_distance_h throws an error about non-Hermitian input.

## Code Examples

### Example 1: Single-Step Cross-Validation (EnergyDomain)

```julia
# Source: codebase analysis of construct_lindbladian + build_trajectoryframework
using StableRNGs

dim = SMALL_DIM  # 8
delta = 0.1

# 1. Exact DM reference via Liouvillian
liouv_config = LiouvConfig(
    num_qubits = 3, with_coherent = false,
    with_linear_combination = true, domain = EnergyDomain(),
    beta = BETA, sigma = SIGMA, a = BETA / 30.0, b = 0.4,
    num_energy_bits = NUM_ENERGY_BITS, w0 = W0, t0 = T0,
    num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
)
L = construct_lindbladian(SMALL_JUMPS, liouv_config, SMALL_HAM)

psi0 = fill(ComplexF64(1.0), dim) / sqrt(dim)
rho0 = psi0 * psi0'
rho_dm = reshape(exp(delta * L) * vec(rho0), dim, dim)
rho_dm = (rho_dm + rho_dm') / 2

# 2. Trajectory average
therm_config = ThermalizeConfig(
    num_qubits = 3, with_coherent = false,
    with_linear_combination = true, domain = EnergyDomain(),
    beta = BETA, sigma = SIGMA, a = BETA / 30.0, b = 0.4,
    num_energy_bits = NUM_ENERGY_BITS, w0 = W0, t0 = T0,
    num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
    mixing_time = delta, delta = delta,
)
precomputed = precompute_data(EnergyDomain(), therm_config, SMALL_HAM)
scratch = KrausScratch(ComplexF64, dim)
fw = build_trajectoryframework(SMALL_JUMPS, SMALL_HAM, therm_config, precomputed, scratch, delta)

ntraj = 10000
rho_traj = zeros(ComplexF64, dim, dim)
Random.seed!(42)
for _ in 1:ntraj
    psi = copy(psi0)
    step_along_trajectory!(psi, fw)
    BLAS.gerc!(one(ComplexF64), psi, psi, rho_traj)  # rho_traj += psi * psi'
end
rho_traj ./= ntraj
rho_traj = (rho_traj + rho_traj') / 2

dist = QuantumFurnace.trace_distance_h(Hermitian(rho_dm), Hermitian(rho_traj))
@test dist < 0.01
```

### Example 2: TrotterDomain with Coherent Term (Single Step)

```julia
dim = SMALL_DIM
delta = 0.1

# Liouvillian for TrotterDomain with coherent term
liouv_config = LiouvConfig(
    num_qubits = 3, with_coherent = true,
    with_linear_combination = true, domain = TrotterDomain(),
    beta = BETA, sigma = SIGMA, a = BETA / 30.0, b = 0.4,
    num_energy_bits = NUM_ENERGY_BITS, w0 = W0, t0 = T0,
    num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
)
L = construct_lindbladian(SMALL_JUMPS, liouv_config, SMALL_HAM; trotter = SMALL_TROTTER)

psi0 = fill(ComplexF64(1.0), dim) / sqrt(dim)
rho0 = psi0 * psi0'
rho_dm = reshape(exp(delta * L) * vec(rho0), dim, dim)
rho_dm = (rho_dm + rho_dm') / 2

# Trajectory framework for TrotterDomain
therm_config = ThermalizeConfig(
    num_qubits = 3, with_coherent = true,
    with_linear_combination = true, domain = TrotterDomain(),
    beta = BETA, sigma = SIGMA, a = BETA / 30.0, b = 0.4,
    num_energy_bits = NUM_ENERGY_BITS, w0 = W0, t0 = T0,
    num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
    mixing_time = delta, delta = delta,
)
precomputed = precompute_data(TrotterDomain(), therm_config, SMALL_TROTTER)
scratch = KrausScratch(ComplexF64, dim)
fw = build_trajectoryframework(SMALL_JUMPS, SMALL_TROTTER, therm_config, precomputed, scratch, delta)

# Note: for TrotterDomain, build_trajectoryframework takes SMALL_TROTTER as ham_or_trott
# The jumps' in_eigenbasis is in H-eigenbasis, which is correct (matches construct_lindbladian)

ntraj = 10000
rho_traj = zeros(ComplexF64, dim, dim)
Random.seed!(42)
for _ in 1:ntraj
    psi = copy(psi0)
    step_along_trajectory!(psi, fw)
    BLAS.gerc!(one(ComplexF64), psi, psi, rho_traj)
end
rho_traj ./= ntraj
rho_traj = (rho_traj + rho_traj') / 2

dist = QuantumFurnace.trace_distance_h(Hermitian(rho_dm), Hermitian(rho_traj))
@test dist < 0.01
```

### Example 3: Multi-Step Gibbs Convergence

```julia
dim = SMALL_DIM
delta = 0.01

# L for TrotterDomain + coherent
liouv_config = LiouvConfig(
    num_qubits = 3, with_coherent = true,
    with_linear_combination = true, domain = TrotterDomain(),
    beta = BETA, sigma = SIGMA, a = BETA / 30.0, b = 0.4,
    num_energy_bits = NUM_ENERGY_BITS, w0 = W0, t0 = T0,
    num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
)
L = construct_lindbladian(SMALL_JUMPS, liouv_config, SMALL_HAM; trotter = SMALL_TROTTER)
exp_L = exp(delta * L)  # Compute once, reuse

# Gibbs in appropriate basis (Trotter eigenbasis for TrotterDomain)
gibbs_comp = SMALL_HAM.eigvecs * SMALL_GIBBS * SMALL_HAM.eigvecs'
gibbs_trott = Hermitian(SMALL_TROTTER.eigvecs' * gibbs_comp * SMALL_TROTTER.eigvecs)

# DM evolution
psi0 = fill(ComplexF64(1.0), dim) / sqrt(dim)
rho_dm_vec = vec(psi0 * psi0')
max_steps = 5000
dm_converged_step = 0
for step in 1:max_steps
    rho_dm_vec = exp_L * rho_dm_vec
    rho_dm = reshape(rho_dm_vec, dim, dim)
    rho_dm = (rho_dm + rho_dm') / 2
    dist = QuantumFurnace.trace_distance_h(Hermitian(rho_dm), gibbs_trott)
    if dist < 1e-3
        dm_converged_step = step
        break
    end
end
@test dm_converged_step > 0  # DM must converge

# Trajectory evolution for same number of steps
therm_config = ThermalizeConfig(
    num_qubits = 3, with_coherent = true,
    with_linear_combination = true, domain = TrotterDomain(),
    beta = BETA, sigma = SIGMA, a = BETA / 30.0, b = 0.4,
    num_energy_bits = NUM_ENERGY_BITS, w0 = W0, t0 = T0,
    num_trotter_steps_per_t0 = NUM_TROTTER_STEPS_PER_T0,
    mixing_time = dm_converged_step * delta,
    delta = delta,
)
ntraj = 10000
result = run_trajectories(SMALL_JUMPS, therm_config, psi0, SMALL_HAM;
    trotter = SMALL_TROTTER, ntraj = ntraj,
    total_time = dm_converged_step * delta, delta = delta,
)
rho_traj = result.rho_mean

# Both must be near Gibbs
dist_dm = QuantumFurnace.trace_distance_h(Hermitian(rho_dm), gibbs_trott)
dist_traj = QuantumFurnace.trace_distance_h(Hermitian(rho_traj), gibbs_trott)
@test dist_dm < 1e-3
@test dist_traj < 1e-3
```

**Important note on `run_trajectories`:** This function already handles multi-step evolution and accumulates `rho_mean` (averaged over ntraj). Looking at the code:
```julia
rho_mean ./= ntraj
rho_mean .= 0.5 .* (rho_mean .+ rho_mean')
```
So it already Hermitianizes. This is the preferred API for multi-step trajectory tests.

However, `run_trajectories` calls `validate_config!` and `print_press` which may print to stdout. For test cleanliness, consider redirecting output or using `@test_nowarn` where appropriate.

## Discretion Recommendations

### Delta Values for Single-Step Sweep

**Recommendation:** `deltas = [0.2, 0.1, 0.05, 0.025]`

Rationale:
- 4 values give 3 ratio checks, enough to verify O(delta^2) scaling
- Delta=0.2 is large enough to produce measurable error above statistical noise floor
- Delta=0.025 is small enough that the error is clearly decreasing
- Halving pattern makes ratio analysis clean (expect ratio ~4 for O(delta^2))
- For dim=8 (3-qubit), CPTP channel is well-behaved at these delta values (verified in Phase 2 CPTP tests)

### Delta Value for Multi-Step Convergence

**Recommendation:** delta = 0.01

Rationale:
- Small enough for DM Euler convergence (O(delta) accumulated error stays below 1e-3)
- Large enough that convergence happens in reasonable number of steps (~100-500)
- Matches TEST_DELTA from test_helpers.jl (already validated in CPTP tests)
- For 3-qubit with TrotterDomain + coherent, spectral gap of Liouvillian determines convergence speed; delta=0.01 should converge in ~200-1000 steps

### Delta^2 Scaling Check Method

**Recommendation:** Consecutive ratio check with tolerance band [2.0, 6.0]

Rationale:
- For O(delta^2) with halving, expect ratio ~4
- Band [2.0, 6.0] allows for statistical noise from finite trajectory count
- Simpler than log-log regression, more interpretable
- 3 ratio checks (from 4 delta values) provide redundancy
- Additionally: compute log-log slope as diagnostic info (print but don't assert)

### Test File Organization

**Recommendation:** Single test file under `test/trajectory_validation/`

```
test/
  trajectory_validation/
    run_trajectory_validation.jl    # Self-contained entry point
```

Rationale:
- One file keeps it simple; two test types (single-step, multi-step) are small enough
- Self-contained: includes test_helpers.jl, adds SMALL_TROTTER, runs all tests
- Can be run via: `julia --project test/trajectory_validation/run_trajectory_validation.jl`
- Not included in runtests.jl (per user decision)

### Maximum Steps for Convergence Test

**Recommendation:** max_steps = 5000 (safety cap)

Rationale:
- At delta=0.01, 5000 steps = 50 time units of evolution
- If the Liouvillian gap is O(0.01) or larger (typical for 3-qubit), convergence should occur well within 5000 steps
- Safety cap prevents infinite loops if convergence fails
- Test should assert that convergence actually happened (dm_converged_step > 0)

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Random jump selection in DM (run_thermalization) | All-jumps Liouvillian as DM reference | Phase 4 (this phase) | Eliminates jump-selection randomness from DM reference |
| Cholesky for PSD guard | Eigendecomposition (Phase 2) | Phase 2 | Trajectory framework builds reliably |
| No coherent term in trajectory | precompute_coherent_total_B sums all B_k | Phase 2 TFIX-02 | U_B applied correctly before Kraus branching |

**Deprecated/outdated:**
- `run_thermalization` for cross-validation: not suitable as DM reference (random jump selection)
- `create_trotter()`: use `TrottTrott(ham, t0, num_steps)` directly

## Open Questions

1. **Basis consistency for TrotterDomain Liouvillian**
   - What we know: construct_lindbladian for TrotterDomain mixes H-eigenbasis (dissipative) and Trotter-eigenbasis (coherent B). The fixed point comparison works in DMTST-02 using Gibbs in Trotter eigenbasis.
   - What's unclear: The exact basis convention for rho when using this L. Is it H-eigenbasis or Trotter-eigenbasis?
   - Recommendation: For single-step cross-validation, basis does not matter -- both rho_dm and rho_traj use the same L/framework, so they are in the same basis. For multi-step Gibbs comparison, use Gibbs in Trotter eigenbasis (following DMTST-02 pattern). If trace distance is unexpectedly large, try Gibbs in H-eigenbasis as a diagnostic.

2. **Performance of 10,000 trajectories for 3-qubit system**
   - What we know: dim=8, so each step_along_trajectory! involves 8x8 matrix-vector products. Very fast.
   - What's unclear: Total wall time for 4 deltas x 3 domains x 10000 trajectories
   - Recommendation: Estimate ~0.1ms per trajectory step -> 4 * 3 * 10000 * 0.1ms = 12 seconds. Plus multi-step convergence: 5000 steps * 10000 traj * 0.1ms = 500 seconds. This might be slow.
   - Mitigation: Multi-step convergence test may need fewer trajectories or fewer steps. With run_trajectories, the 10000 trajectories are run sequentially (single-threaded). For 5000 steps * 10000 traj at dim=8, this is 50 million step calls. At ~microseconds per call, ~50 seconds. Should be acceptable.

3. **run_trajectories API compatibility**
   - What we know: run_trajectories builds the framework internally and runs ntraj trajectories. It returns a named tuple with rho_mean.
   - What's unclear: Whether run_trajectories re-seeds the RNG or uses whatever global state exists.
   - Recommendation: Seed RNG before calling run_trajectories. The function does not touch the RNG seed internally; it just calls step_along_trajectory! which calls rand().

## Sources

### Primary (HIGH confidence)
- Codebase direct analysis: `src/trajectories.jl`, `src/furnace.jl`, `src/jump_workers.jl`, `src/coherent.jl`, `src/structs.jl`, `src/kraus.jl`, `src/furnace_utensils.jl`, `src/energy_domain.jl`, `src/trotter_domain.jl`, `src/qi_tools.jl`
- Existing test infrastructure: `test/test_helpers.jl`, `test/test_cptp.jl`, `test/test_trajectory_fixes.jl`, `test/test_dm_detailed_balance.jl`, `test/test_dm_scaling.jl`, `test/runtests.jl`
- Phase 2 summaries: `.planning/phases/02-trajectory-bug-fixes/02-01-SUMMARY.md` (StableRNG limitation discovered)
- Phase 3 research: `.planning/phases/03-dm-reference-test-suite/03-RESEARCH.md` (DM patterns, basis conventions)
- Project.toml: StableRNGs in [extras] with compat "1", available in test target

### Secondary (MEDIUM confidence)
- Phase 2 TFIX-02 fix confirmed U_B ordering matches DM code
- Phase 3 DMTST-02 confirms TrotterDomain Liouvillian consistency using Gibbs in Trotter eigenbasis

### Tertiary (LOW confidence)
- Performance estimates for 10,000 trajectories at 3-qubit (not benchmarked, estimated from operation counts)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - all libraries already in Project.toml, no new deps needed
- Architecture: HIGH - all patterns derived from direct codebase analysis of trajectories.jl, furnace.jl, jump_workers.jl; critical DM-vs-trajectory channel difference identified from code
- Pitfalls: HIGH - identified from analyzing Phase 2 RNG issues, Phase 3 basis conventions, and DM channel randomness
- Code examples: HIGH - patterns verified against existing test patterns (Phase 2/3) and function signatures in source

**Research date:** 2026-02-14
**Valid until:** 2026-03-14 (stable domain, no external dependency changes expected)
