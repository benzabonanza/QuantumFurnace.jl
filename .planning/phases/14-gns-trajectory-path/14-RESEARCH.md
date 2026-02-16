# Phase 14: GNS Trajectory Path - Research

**Researched:** 2026-02-16
**Domain:** GNS (approximate detailed balance) trajectory simulation validation in QuantumFurnace.jl
**Confidence:** HIGH

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions

- **GNS reference state**: Compute the GNS fixed point using the existing `run_lindbladian()` with a GNS config, which returns `HotSpectralResults` (fixed point, spectral gap, first excited mode). No separate null-space solver needed -- `run_lindbladian()` already provides this via `construct_lindbladian()` dispatching on config type. The GNS Lindbladian is built by `construct_lindbladian()` when given a `LiouvConfigGNS` -- this is already implemented. Also compute and document the trace distance between the GNS fixed point and the exact Gibbs state (the approximation gap) as a baseline for Phase 18.

- **Sigma parameter design**: Sigma already exists as a field on GNS config structs. Phase 14 validation tests use sigma = 1/beta only. The two-sigma comparison (1/beta and 0.5/beta) is Phase 18 scope. Test at n=3 (dim=8) only -- small system for fast tests, matches existing test fixtures. beta = 10.0 minimum for all simulations -- lower beta is high temperature and can mask convergence errors.

- **Approximation tolerance**: Use trace distance (`trace_distance_h()`) as the convergence metric. Target: trace distance < 0.05 between averaged trajectory DM and GNS fixed point. Adjust delta step size (try 0.1 or 0.01) and trajectory count to achieve this threshold. Single delta value sufficient -- no delta sweep needed for Phase 14. Validate DM properties (Hermitian, unit trace, PSD) on the final averaged result only. If final result fails validation, escalate to batch checkpoint checks to isolate where it breaks.

- **B-term handling**: `with_coherent=false` is already the default for GNS configs -- no B-term construction in GNS runs. `step_along_trajectory!` already respects `with_coherent=false` and skips the unitary step. `pick_transition()` dispatches on GNS config type and uses GNS-specific transition weights gamma(omega). No explicit B-term absence tests needed -- the GNS Lindbladian converges to its own fixed point regardless of B presence. The entire GNS code path (jump construction, B-term suppression, transition selection) is config-driven and already implemented.

### Claude's Discretion

- Exact delta step size and trajectory count to achieve < 0.05 trace distance
- Test structure and organization (single test file vs integrated into existing test suite)
- Whether to test BohrDomain, EnergyDomain, or TrotterDomain for GNS (or multiple)

### Deferred Ideas (OUT OF SCOPE)

- Two-sigma comparison (sigma=1/beta vs sigma=0.5/beta) demonstrating that reducing sigma improves GNS approximation -- Phase 18
- Larger system sizes (n=4,6) for GNS validation -- Phase 18 parameter grid
- Delta sweep showing convergence improvement with smaller steps -- future work if needed
</user_constraints>

## Summary

Phase 14 validates the GNS (approximate detailed balance) trajectory code path end-to-end. The codebase already has complete GNS support implemented: `LiouvConfigGNS` and `ThermalizeConfigGNS` struct types, GNS-specific `pick_transition()` dispatch returning unshifted gamma(omega), `_pick_alpha_gns()` for BohrDomain Kossakowski matrix entries, and `with_coherent=false` enforcement in config constructors. The `construct_lindbladian()`, `run_lindbladian()`, `build_trajectoryframework()`, `step_along_trajectory!()`, and `run_trajectories()` functions all operate on the abstract `AbstractLiouvConfig` and `AbstractThermalizeConfig` types, meaning they already accept GNS configs without modification. However, **no tests exist anywhere in the test suite that use GNS config types**. Every existing test uses KMS configs (`LiouvConfig`, `ThermalizeConfig`).

The primary task is therefore verification testing, not new functionality. The tests need to: (1) construct a GNS Lindbladian via `construct_lindbladian()` with `LiouvConfigGNS` and verify its fixed point, (2) build a `TrajectoryFramework` from `ThermalizeConfigGNS` and verify CPTP completeness of per-operator channels, (3) run trajectory averaging via `run_trajectories()` and verify convergence to the GNS fixed point with trace distance < 0.05, and (4) document the GNS-to-Gibbs approximation gap as a Phase 18 baseline.

The test infrastructure is well-established. The 3-qubit test system (`SMALL_HAM`, `SMALL_JUMPS`, `SMALL_GIBBS`, `SMALL_DIM=8`) already exists in `test_helpers.jl` with beta=10.0 and is the correct fixture for Phase 14. Factory functions for GNS configs need to be added to `test_helpers.jl` following the existing `make_small_liouv_config()` and `make_small_thermalize_config()` patterns, but using `LiouvConfigGNS` and `ThermalizeConfigGNS` with `with_linear_combination=true`, `a=beta/30.0`, `b=0.4` (matching KMS test parameters).

**Primary recommendation:** Add `make_small_liouv_config_gns()` and `make_small_thermalize_config_gns()` to `test_helpers.jl`, then create a single new test file `test_gns_trajectory.jl` that validates the full GNS pipeline: Lindbladian fixed point, CPTP completeness, trajectory convergence, and DM validity. Use EnergyDomain as the primary domain (matching Phase 14's KMS counterpart tests, simplest dispatch, no NUFFT prefactors needed), with optional BohrDomain detailed balance test to mirror DMTST-01.

## Standard Stack

### Core

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| `QuantumFurnace.jl` | current | `LiouvConfigGNS`, `ThermalizeConfigGNS`, `construct_lindbladian`, `run_lindbladian`, `run_trajectories`, `build_trajectoryframework`, `step_along_trajectory!` | The library under test; all GNS dispatch paths already implemented |
| `LinearAlgebra` (stdlib) | Julia 1.11+ | `eigen()`, `Hermitian()`, `tr()`, `norm()`, `I`, `BLAS` | Dense eigendecomposition for Lindbladian fixed points, matrix operations |
| `Random` (stdlib) | Julia 1.11+ | `Xoshiro(seed)` | Deterministic RNG seeding for trajectory reproducibility |
| `Test` (stdlib) | Julia 1.11+ | `@testset`, `@test`, `isapprox` | Standard Julia test framework |
| `Arpack` | 0.5+ | Sparse eigendecomposition via `eigs` | Used internally by `run_lindbladian()` for fixed point extraction |

### Supporting

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| `BSON` | 0.3+ | Test Hamiltonian loading | Already used by `test_helpers.jl` for `_load_test_hamiltonian` |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| `construct_lindbladian` + `eigen` for reference | `run_lindbladian` (uses Arpack `eigs`) | `run_lindbladian` is the locked decision and gives spectral gap for free; full `eigen` is fine for dim=64 (3-qubit) but `run_lindbladian` follows the user decision |
| EnergyDomain only | Multiple domains (Bohr, Energy, Time, Trotter) | EnergyDomain is simplest for trajectory validation; BohrDomain adds value for detailed balance check; Time/Trotter add testing depth but may be overkill for Phase 14's scope |
| Single test file | Spread across existing test files | Single file is cleaner for a new code path validation; can always integrate later |

## Architecture Patterns

### GNS Config Type Hierarchy

```
AbstractConfig{D,T}
  AbstractLiouvConfig{D,T}
    LiouvConfig{D,T}        # KMS (with coherent B term option)
    LiouvConfigGNS{D,T}     # GNS (with_coherent=false enforced in constructor)
  AbstractThermalizeConfig{D,T}
    ThermalizeConfig{D,T}   # KMS
    ThermalizeConfigGNS{D,T} # GNS (with_coherent=false enforced in constructor)
```

All dispatch operates on abstract types (`AbstractLiouvConfig`, `AbstractThermalizeConfig`), except:
- `pick_transition()`: 4 methods dispatching on concrete types (LiouvConfig, LiouvConfigGNS, ThermalizeConfig, ThermalizeConfigGNS)
- `_pick_alpha()`: 4 methods dispatching on concrete types (same pattern)
- `_select_b_plus_calculator()`: only defined for KMS types (GNS never reaches it because `config.with_coherent == false`)
- `validate_config!()`: special GNS checks (with_coherent must be false, different beta constraint for Gaussian parameters)

### GNS Config Construction Pattern

For `with_linear_combination=true` (Smooth Metro variant, matching existing KMS test params):

```julia
LiouvConfigGNS(
    num_qubits = 3,
    with_coherent = false,          # enforced by constructor
    with_linear_combination = true,
    domain = EnergyDomain(),
    beta = 10.0,
    sigma = 0.1,                    # 1/beta
    a = 10.0 / 30.0,               # beta/30
    b = 0.4,
    num_energy_bits = 12,
    w0 = 0.05,
    t0 = 2pi / (2^12 * 0.05),      # only needed for Time/Trotter domains
    num_trotter_steps_per_t0 = 10,  # only needed for TrotterDomain
)
```

For `ThermalizeConfigGNS`, add `mixing_time` and `delta` fields.

**Critical parameter note for GNS Gaussian case:** If `with_linear_combination=false`, the constraint is `beta = 2*w_gamma/sigma_gamma^2` (NOT `beta = 2*w_gamma/(sigma^2 + sigma_gamma^2)` as in KMS). Since the locked decision uses `with_linear_combination=true`, this constraint is not relevant for Phase 14.

### GNS Transition Function (Unshifted)

The key difference from KMS: GNS transition weights are **unshifted**. For the Smooth Metro case (`a!=0, b!=0`):

```
KMS: sqrtB = sqrt(beta/4) * |w + beta*sigma^2/2|   (shifted by beta*sigma^2/2)
GNS: sqrtB = sqrt(beta/4) * |w|                     (unshifted)
```

This means the GNS Lindbladian has a **different fixed point** than the KMS Lindbladian. The GNS fixed point is an approximation of the Gibbs state, with the approximation quality depending on sigma. Phase 14 measures this gap; Phase 18 explores it.

### GNS Alpha (Kossakowski Matrix) for BohrDomain

Similarly unshifted in the `create_alpha_gns` function:
```
KMS: sqrtB = sqrt(beta/16) * |nu_1 + nu_2|
GNS: sqrtB = sqrt(beta/16) * |nu_1 + nu_2 + beta*sigma^2/2|
```

Note: the Kossakowski alpha for GNS is *partially shifted* (by `beta*sigma^2/2`), which is a consequence of integrating the unshifted gamma against Gaussian filters.

### Test Pattern: GNS Lindbladian Fixed Point

Follow the DMTST-01 pattern but with GNS config:

```julia
# Build GNS Lindbladian
liouv_config_gns = make_small_liouv_config_gns(EnergyDomain())
L_gns = construct_lindbladian(SMALL_JUMPS, liouv_config_gns, SMALL_HAM)

# Full eigendecomposition (64x64 dense matrix)
eig = eigen(L_gns)
ss_idx = argmin(abs.(real.(eig.values)))
ss_vec = eig.vectors[:, ss_idx]
ss_dm = reshape(ss_vec, SMALL_DIM, SMALL_DIM)
ss_dm = (ss_dm + ss_dm') / 2
ss_dm ./= tr(ss_dm)

# GNS fixed point is NOT the Gibbs state -- measure the gap
gap = trace_distance_h(Hermitian(ss_dm), SMALL_GIBBS)
@info "GNS-to-Gibbs approximation gap" gap
# Expect gap > 0 but bounded (this is the "approximation bound" for sigma=0.1)
```

Alternatively, use `run_lindbladian()` as per the locked decision:

```julia
result = run_lindbladian(SMALL_JUMPS, liouv_config_gns, SMALL_HAM)
gns_fixed_point = result.fixed_point
spectral_gap = result.spectral_gap
gap_to_gibbs = trace_distance_h(Hermitian(gns_fixed_point), SMALL_GIBBS)
```

### Test Pattern: GNS Trajectory Convergence

Follow the TVAL-06 pattern (TrotterDomain convergence to Gibbs) but targeting the GNS fixed point:

```julia
# Get GNS reference
result_lindblad = run_lindbladian(SMALL_JUMPS, liouv_config_gns, SMALL_HAM)
gns_fp = result_lindblad.fixed_point

# Run trajectories
therm_config_gns = make_small_thermalize_config_gns(EnergyDomain();
    delta=0.01, mixing_time=5.0)
psi0 = fill(ComplexF64(1.0), SMALL_DIM) / sqrt(SMALL_DIM)
result_traj = run_trajectories(SMALL_JUMPS, therm_config_gns, psi0, SMALL_HAM;
    ntraj=1000, seed=42)
rho_traj = result_traj.rho_mean

# Convergence to GNS fixed point (NOT Gibbs)
dist = trace_distance_h(Hermitian(rho_traj), Hermitian(gns_fp))
@test dist < 0.05
```

### Anti-Patterns to Avoid

- **Comparing GNS trajectories to the Gibbs state directly:** The GNS trajectory converges to the GNS Lindbladian fixed point, NOT the exact Gibbs state. Comparing to Gibbs gives a nonzero distance even with infinite trajectories. Always compare to the GNS Lindbladian fixed point for convergence testing.
- **Using `with_coherent=true` with GNS configs:** The constructor will throw an error. GNS configs enforce `with_coherent=false`.
- **Using `with_linear_combination=false` without setting gaussian_parameters:** The validator will throw. If using Gaussian transition variant, must set `gaussian_parameters=(w_gamma, sigma_gamma)` where `beta = 2*w_gamma/sigma_gamma^2` for GNS.
- **Testing with beta < 10:** High temperature can mask convergence errors -- the locked decision specifies beta >= 10.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| GNS fixed point computation | Custom null-space solver | `run_lindbladian(jumps, LiouvConfigGNS_config, ham)` | Returns `HotSpectralResults` with fixed point, spectral gap, gap mode -- complete diagnostic package |
| Trace distance computation | Manual eigenvalue computation | `trace_distance_h(Hermitian(A), Hermitian(B))` | Already handles Hermitian wrapper, uses eigvals for efficient computation |
| DM validity checking | Custom eigenvalue/trace checks | `is_density_matrix(Hermitian(rho))` | Checks Hermiticity, non-negative eigenvalues, unit trace with proper tolerances |
| Trajectory averaging | Manual psi*psi' accumulation loop | `run_trajectories(jumps, config, psi0, ham; ntraj=N, seed=S)` | Handles workspace allocation, threading, seeding, density matrix averaging, Hermitianization |
| GNS config factory | Inline config construction in tests | Factory function in `test_helpers.jl` | Centralizes parameter choices, matches existing pattern (`make_small_liouv_config`, `make_small_thermalize_config`) |

**Key insight:** The entire GNS code path is already implemented and config-driven. Phase 14's job is verification, not construction. Use the existing infrastructure directly.

## Common Pitfalls

### Pitfall 1: Wrong reference state for convergence testing
**What goes wrong:** Tests compare trajectory-averaged DM to the Gibbs state and get a persistent nonzero distance even with many trajectories.
**Why it happens:** The GNS Lindbladian has a different fixed point than the exact KMS Lindbladian. The GNS fixed point approximates the Gibbs state, but there is always an approximation gap for nonzero sigma.
**How to avoid:** Always compute the GNS Lindbladian fixed point via `run_lindbladian()` with `LiouvConfigGNS` and compare trajectories to THAT reference. Separately measure and document the GNS-to-Gibbs gap.
**Warning signs:** Trajectory distance plateaus at a nonzero value that doesn't decrease with more trajectories.

### Pitfall 2: GNS config validation failure
**What goes wrong:** Config construction or validation throws errors.
**Why it happens:** The GNS config validators enforce different constraints than KMS. Specifically: (1) `with_coherent` must be `false`, (2) for Gaussian transitions, the parameter relation is `beta = 2*w_gamma/sigma_gamma^2` (not `beta = 2*w_gamma/(sigma^2 + sigma_gamma^2)`).
**How to avoid:** Use `with_linear_combination=true` (Smooth Metro) to sidestep the Gaussian parameter constraint entirely. Set `with_coherent=false` explicitly (or rely on the default).
**Warning signs:** `ArgumentError` from `validate_config!`.

### Pitfall 3: delta_eff overflow for many jump operators
**What goes wrong:** `build_trajectoryframework` asserts `delta_eff < 1.0` and throws.
**Why it happens:** `delta_eff = delta * n_jumps`. For n=3 there are 9 jump operators, so delta=0.1 gives delta_eff=0.9 (fine) but delta=0.12 gives delta_eff=1.08 (fails).
**How to avoid:** With 9 jumps, delta must be < 1/9 ~ 0.111. Use delta=0.1 (delta_eff=0.9) or delta=0.01 (delta_eff=0.09). The smaller delta gives better Lie-Trotter accuracy but more steps.
**Warning signs:** Assertion error from `build_trajectoryframework`.

### Pitfall 4: Insufficient trajectories or mixing time
**What goes wrong:** Trace distance to GNS fixed point exceeds the 0.05 threshold.
**Why it happens:** Statistical noise from too few trajectories (1/sqrt(N) scaling) or insufficient mixing time (system hasn't converged).
**How to avoid:** Start with ntraj=1000 and mixing_time=5.0 (matching Phase 13 performance test parameters). If trace distance exceeds 0.05, increase ntraj to 2000-5000 or increase mixing_time. The spectral gap from `run_lindbladian()` gives the convergence time scale: mixing_time should be several multiples of 1/|spectral_gap|.
**Warning signs:** Trace distance decreasing but still above threshold; high variance across different seeds.

### Pitfall 5: BohrDomain _select_b_plus_calculator missing GNS dispatch
**What goes wrong:** If using BohrDomain with GNS, the coherent B-term code path should never be reached (since `with_coherent=false`), but `_select_b_plus_calculator` only has methods for `LiouvConfig`/`ThermalizeConfig`.
**Why it happens:** The function was written for KMS only. GNS never needs it because GNS never uses coherent terms.
**How to avoid:** This is not a problem as long as `with_coherent=false` (enforced by GNS constructors). The coherent path is gated by `config.with_coherent || return nothing` in `_precompute_coherent_total_B`.
**Warning signs:** None in practice -- the code path is never reached for GNS.

## Code Examples

### GNS Config Factory Functions (to add to test_helpers.jl)

```julia
"""
    make_small_liouv_config_gns(domain) -> LiouvConfigGNS

Create a LiouvConfigGNS for the 3-qubit SMALL system.
Uses Smooth Metro transition (with_linear_combination=true, a=beta/30, b=0.4)
matching the KMS test parameter choices.
"""
function make_small_liouv_config_gns(domain)
    LiouvConfigGNS(
        num_qubits = 3,
        with_coherent = false,
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

"""
    make_small_thermalize_config_gns(domain; delta=TEST_DELTA, mixing_time=1.0) -> ThermalizeConfigGNS

Create a ThermalizeConfigGNS for the 3-qubit SMALL system.
"""
function make_small_thermalize_config_gns(domain;
    delta::Float64=TEST_DELTA,
    mixing_time::Float64=1.0,
)
    ThermalizeConfigGNS(
        num_qubits = 3,
        with_coherent = false,
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
        mixing_time = mixing_time,
        delta = delta,
    )
end
```

### GNS Lindbladian Fixed Point with run_lindbladian

```julia
# Locked decision: use run_lindbladian() for GNS reference
liouv_config_gns = make_small_liouv_config_gns(EnergyDomain())
result = run_lindbladian(SMALL_JUMPS, liouv_config_gns, SMALL_HAM)

gns_fixed_point = result.fixed_point       # Matrix{ComplexF64}
spectral_gap = result.spectral_gap         # Complex{Float64}
gap_mode = result.gap_mode                 # Matrix{ComplexF64}

# Document the approximation gap
gap_to_gibbs = trace_distance_h(Hermitian(gns_fixed_point), SMALL_GIBBS)
@info "GNS-to-Gibbs gap" gap=gap_to_gibbs spectral_gap=abs(real(spectral_gap))
```

### GNS Trajectory CPTP Verification

```julia
therm_config_gns = make_small_thermalize_config_gns(EnergyDomain(); delta=0.01)
precomputed = QuantumFurnace._precompute_data(therm_config_gns, SMALL_HAM)
scratch = QuantumFurnace.KrausScratch(ComplexF64, SMALL_DIM)
fw = build_trajectoryframework(
    SMALL_JUMPS, SMALL_HAM, therm_config_gns, precomputed, scratch, 0.01
)

identity = Matrix{ComplexF64}(I, SMALL_DIM, SMALL_DIM)
for per_op in fw.per_operator
    completeness = per_op.K0' * per_op.K0 + fw.delta_eff * per_op.R + per_op.U_residual' * per_op.U_residual
    @test isapprox(completeness, identity; atol=1e-10)
    @test per_op.U_B === nothing  # No coherent term for GNS
end
```

### GNS Trajectory Convergence Test

```julia
psi0 = fill(ComplexF64(1.0), SMALL_DIM) / sqrt(SMALL_DIM)
therm_config_gns = make_small_thermalize_config_gns(EnergyDomain();
    delta=0.01, mixing_time=5.0)

result = run_trajectories(SMALL_JUMPS, therm_config_gns, psi0, SMALL_HAM;
    ntraj=1000, seed=42)
rho_traj = result.rho_mean

# Validate DM properties on final result
@test isapprox(tr(rho_traj), 1.0; atol=1e-10)
@test isapprox(rho_traj, rho_traj'; atol=1e-10)  # Hermitian
eigs = eigvals(Hermitian(rho_traj))
@test all(v -> v >= -1e-14, eigs)  # PSD

# Convergence to GNS fixed point
dist = trace_distance_h(Hermitian(rho_traj), Hermitian(gns_fixed_point))
@test dist < 0.05
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Shared mutable RNG | Per-trajectory `Xoshiro(seed + traj_id)` seeding | Phase 13 | Enables deterministic threaded trajectories; GNS tests inherit this |
| `Vector{JumpOp}` (abstract elements) | `Vector{JumpOp{Matrix{T}}}` (concrete elements) | Phase 13 | Zero-allocation hot path; GNS benefits without changes |
| Workspace mixed into TrajectoryFramework | Explicit ws/rng arguments | Phase 12 | Thread-safe stepping; GNS tests use same `step_along_trajectory!(psi, fw, ws, rng)` API |
| Abstract config field access in step loop | Concrete `scaled_prefactor`, `sigma`, `transition` fields on framework | Phase 13 | Hot-path performance; GNS transition function stored the same way |

**Deprecated/outdated:**
- `run_thermalization`: DM-based step-by-step evolution. Still works but `run_trajectories` is the preferred trajectory API for Phase 14.

## Discretion Recommendations

### Delta step size and trajectory count

**Recommendation: delta=0.01, ntraj=1000, mixing_time=5.0**

Reasoning:
- With 9 jumps (n=3), delta_eff = 0.01 * 9 = 0.09 (well within the < 1.0 constraint)
- The Lie-Trotter splitting error is O(delta), so delta=0.01 gives ~1% systematic bias
- 1/sqrt(1000) ~ 0.03 statistical noise, combining with ~0.01 splitting bias gives total error ~0.04 (under the 0.05 threshold)
- mixing_time=5.0 with delta=0.01 gives 500 steps per trajectory -- enough to mix for the dim=8 system at beta=10
- These parameters match the Phase 13 threading performance test (ntraj=200 x mixing_time=5.0), so we know they work
- If 0.05 threshold is tight, ntraj=2000 reduces statistical noise to 1/sqrt(2000) ~ 0.022

### Test structure

**Recommendation: Single new test file `test_gns_trajectory.jl` added to runtests.jl**

Reasoning:
- The GNS validation is a distinct code path that has never been tested
- A dedicated file makes it easy to identify, run in isolation, and reference
- Follow the existing naming pattern (`test_*.jl`)
- Add it to `runtests.jl` alongside the other `include()` calls
- Factory functions go in `test_helpers.jl` (shared infrastructure)

### Domain choice

**Recommendation: EnergyDomain as primary, with BohrDomain detailed balance as secondary**

Reasoning:
- **EnergyDomain**: Matches the existing trajectory regression tests (Phase 5), uses `oft!` for jump operators which is simpler than NUFFT prefactors, and exercises the most common GNS trajectory code path. The `_precompute_R` and `step_along_trajectory!` both have EnergyDomain-specific methods that handle the Gaussian filter, making it a good coverage target.
- **BohrDomain**: Valuable for a GNS detailed balance test (analogous to DMTST-01) -- verifying that the BohrDomain GNS Lindbladian has a well-defined fixed point and documenting its distance from Gibbs. This is purely a Lindbladian test (no trajectories), so it is fast and provides the approximation gap baseline for Phase 18.
- **TimeDomain/TrotterDomain**: Not recommended for Phase 14. These add NUFFT prefactor computation and (for Trotter) basis transforms, which increases test time and complexity without testing GNS-specific logic. The GNS dispatch difference is entirely in `pick_transition()` and `_pick_alpha()`, which are domain-independent.

## Open Questions

1. **Expected GNS-to-Gibbs approximation gap magnitude**
   - What we know: The gap depends on sigma. For sigma=1/beta=0.1, beta=10, n=3, the gap is nonzero but should be bounded.
   - What's unclear: The exact numerical value of the gap. It could be anywhere from 0.001 to 0.1.
   - Recommendation: The test should compute and log the gap, then assert it is bounded (e.g., < 0.5). The exact value becomes the Phase 18 baseline. If the gap is very small (< 0.01), it means sigma=0.1 is a good approximation, which is useful context.

2. **Spectral gap of GNS Lindbladian relative to KMS**
   - What we know: `run_lindbladian()` returns the spectral gap, which determines the mixing time.
   - What's unclear: Whether the GNS spectral gap is similar to KMS or very different, which affects how much mixing_time is needed.
   - Recommendation: Log the spectral gap in test output. If mixing_time=5.0 is insufficient (trajectory doesn't converge), increase it based on the spectral gap.

## Sources

### Primary (HIGH confidence)
- Source code analysis: `src/structs.jl` lines 105-216 -- GNS config struct definitions with constructor validation
- Source code analysis: `src/energy_domain.jl` lines 1-104 -- `pick_transition` dispatch and `_pick_transition_gns` implementation
- Source code analysis: `src/bohr_domain.jl` lines 82-149 -- `_pick_alpha` dispatch and `_pick_alpha_gns`/`create_alpha_gns` implementations
- Source code analysis: `src/trajectories.jl` -- Full trajectory pipeline (`TrajectoryFramework`, `TrajectoryWorkspace`, `step_along_trajectory!`, `run_trajectories`)
- Source code analysis: `src/furnace.jl` lines 1-93 -- `run_lindbladian` and `construct_lindbladian` using abstract type dispatch
- Source code analysis: `src/coherent.jl` lines 1-41 -- `_precompute_coherent_total_B` gated by `config.with_coherent`
- Source code analysis: `src/misc_tools.jl` lines 117-170 -- `validate_config!` with GNS-specific checks
- Source code analysis: `src/qi_tools.jl` lines 128-197 -- `trace_distance_h`, `is_density_matrix`, `gibbs_state`
- Source code analysis: `test/test_helpers.jl` -- All existing test infrastructure, fixtures, factory functions
- Source code analysis: `test/runtests.jl` -- Test suite structure
- Source code analysis: All existing test files -- Confirmed zero GNS config usage in test suite

### Secondary (MEDIUM confidence)
- Phase 13 research and decisions: Threading infrastructure, parameter choices (ntraj=200, mixing_time=5.0 for dim=8)
- Phase context decisions: Locked parameter choices for sigma, beta, tolerance thresholds

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Direct source code analysis, all functions verified to exist and accept GNS types
- Architecture: HIGH - Type hierarchy and dispatch paths traced through source code
- Pitfalls: HIGH - Derived from source code analysis of validation logic and known constraints
- Code examples: HIGH - Modeled on existing working test patterns in the test suite

**Research date:** 2026-02-16
**Valid until:** 2026-03-16 (stable codebase, no external dependencies)
