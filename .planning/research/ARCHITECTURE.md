# Architecture Research: Trajectory Validation and Test Suite Integration

**Domain:** Trajectory simulation validation + comprehensive test suite for QuantumFurnace.jl
**Researched:** 2026-02-13
**Confidence:** HIGH (direct codebase analysis -- every source file read in full)

## System Overview: How New Components Fit the Existing Architecture

```
EXISTING (working)                          NEW / MODIFIED (this milestone)
=================                          ================================

+-- furnace.jl --------------------------------+   +-- test/runtests.jl --------------------+
|  run_lindbladian() [DM spectral, KNOWN GOOD] |   |  @testset "QuantumFurnace" begin       |
|  run_thermalization() [DM stepping, KNOWN G.] |   |    include("test_helpers.jl")           |
|  construct_lindbladian()                      |   |    include("test_detailed_balance.jl")  |
+-------+--------------------------+-----------+   |    include("test_coherent_term.jl")     |
        |                          |               |    include("test_oft_consistency.jl")   |
        v                          v               |    include("test_dm_step_errors.jl")    |
+-- jump_workers.jl ---+  +-- trajectories.jl -+   |    include("test_trajectory.jl")        |
|  jump_contribution!  |  |  TrajectoryFW      |   |    include("test_traj_vs_dm.jl")        |
|  (Liouv + Kraus DM)  |  |  step_along_traj!  |   |  end                                   |
|  precompute_R()      |  |  run_trajectories() |   +------------------------------------------+
|  [KNOWN GOOD]        |  |  [NEEDS FIX+VALID.] |
+----------------------+  +----+----------------+   +-- test/test_helpers.jl -----------------+
                               |                    |  make_test_system(nq, beta, domain)     |
                               |                    |  Known-good parameter configs            |
                       FIX: two-stage sampling      |  Tolerance constants per error tier      |
                       in step_along_trajectory!     +------------------------------------------+
```

### What Exists vs What Needs to Change

| Component | File | Status | Action |
|-----------|------|--------|--------|
| **HamHam + finalize** | `hamiltonian.jl` | WORKING | No changes. Used by tests as-is. |
| **Domain types** | `structs.jl` | WORKING | No changes. |
| **Config structs** | `structs.jl` | WORKING | No changes. |
| **JumpOp construction** | `hamiltonian.jl` + scripts | WORKING | Extract into helper for tests. |
| **precompute_data()** | `furnace_utensils.jl` | WORKING | No changes. |
| **coherent_bohr()** | `bohr_domain.jl` | WORKING | Test target (reference for B comparisons). |
| **B_time(), B_trotter()** | `coherent.jl` | WORKING | Test target (compare against coherent_bohr). |
| **oft!()** | `ofts.jl` | WORKING | Test target (reference for OFT comparisons). |
| **time_oft!(), trotter_oft!()** | `ofts.jl` | WORKING | Test target (compare against oft!). |
| **jump_contribution!(Liouv)** | `jump_workers.jl` | WORKING | Test target (Liouvillian construction). |
| **jump_contribution!(Kraus DM)** | `jump_workers.jl` | WORKING | Test target (DM stepping). |
| **run_lindbladian()** | `furnace.jl` | WORKING | Used by tests as reference oracle. |
| **run_thermalization()** | `furnace.jl` | WORKING | Used by tests as DM reference. |
| **build_trajectoryframework()** | `trajectories.jl` | **BUG** | Fix: `trotter` var not in scope (line 53). |
| **precompute_R()** | `trajectories.jl` | WORKING | No changes needed. |
| **step_along_trajectory!()** | `trajectories.jl` | **NEEDS FIX** | Two-stage jump sampling restructure. |
| **run_trajectories()** | `trajectories.jl` | WORKING (modulo step bug) | No structural changes. |
| **test/*.jl** | `test/` | INTERACTIVE ONLY | Convert to @testset, add runtests.jl. |
| **~200 lines commented code** | `jump_workers.jl`, `trajectories.jl` | DEAD CODE | Clean up after validation. |

## The Two-Stage Jump Sampling Fix

### Current Architecture (Flat Sampling -- Incorrect)

The current `step_along_trajectory!` iterates over all `(jump, omega)` pairs in a flat scan:

```
for jump in fw.jumps           # outer: iterate over Lindblad operators A^a
    for w in energy_labels     # inner: iterate over frequencies omega
        # Compute A_{a,omega} and its probability
        p = delta * rate(omega) * ||A_{a,omega} psi||^2
        csum += p
        if csum >= target: pick this (a, omega)
```

This is a flat cumulative-probability scan over the Cartesian product `{A^a} x {omega}`. The probability of picking jump `a` at frequency `omega` is:

```
P(a, omega) = delta * rate(omega) * ||A_{a,omega} |psi>||^2 / p_jump_total
```

The problem: `p_jump_total = delta * <psi|R|psi>` where `R = sum_{a,omega} rate(omega) * A_{a,omega}^dag A_{a,omega}`. This means the probability of picking operator `a` depends on the *state-dependent* norms `||A_{a,omega} psi||^2` rather than on the mathematical structure of the Lindbladian. The theory (Chen et al.) requires a specific factorization.

### Required Architecture (Two-Stage Sampling -- Correct)

The Chen construction decomposes the CPTP map as a sum over Lindblad operator indices `a`. The quantum algorithm physically implements:

**Stage 1:** Pick which Lindblad operator `A^a` to apply. In the implementation, this means uniformly sampling a jump index `a` (since all jump operators are weighted equally via `1/sqrt(num_jumps)` normalization).

**Stage 2:** Given `a`, sample a frequency `omega` from the conditional distribution:

```
P(omega | a) proportional to rate(omega) * ||A_{a,omega} |psi>||^2
```

Then apply the corresponding Kraus operator `sqrt(rate(omega)) * A_{a,omega}` to `|psi>`.

### Implementation Plan

The fix modifies `step_along_trajectory!` in `trajectories.jl`. The key change is in the dissipative-jump branch (the `else` block after no-jump and residual):

```julia
# CURRENT: flat scan over all (a, omega) pairs
# PROPOSED: two-stage

# Stage 1: Pick jump operator index uniformly
a_idx = rand(1:length(fw.jumps))
jump = fw.jumps[a_idx]

# Stage 2: Sample omega from conditional distribution for this a
# (cumulative probability scan over energy_labels only, not over all jumps)
target_omega = rand() * p_jump_a  # p_jump_a = sum_omega p(a, omega)
csum = 0.0
for w in energy_labels
    # build A_{a,omega}, compute ||A_{a,omega} psi||^2
    p = delta * rate(omega) * n2  # contribution from this omega
    csum += p
    if csum >= target_omega: apply A_{a,omega} and break
end
```

**Critical subtlety:** The `p_jump_total` and `p_nojump` computations use the *total* `R = sum_a sum_omega ...`, which is already correct and does not change. What changes is only the *selection* of which jump operator to apply when a dissipative jump occurs.

**What does NOT change:**
- `precompute_R()` -- still computes the full R summed over all jumps and frequencies
- `K0 = I - alpha*R` -- unchanged
- `U_residual` -- unchanged
- `p_nojump`, `p_res`, `p_jump_total` -- unchanged
- The coherent unitary `U_B` -- unchanged

**What changes:**
- The dissipative-jump branch in `step_along_trajectory!` (both EnergyDomain and Time/TrotterDomain variants)
- Need per-jump `R_a` or per-jump `p_jump_a = delta * <psi|R_a|psi>` for Stage 1 weighting (not uniform -- proportional to `<psi|R_a|psi>`)

### Detailed Changes to `step_along_trajectory!`

```julia
# In the dissipative-jump branch:

# Stage 1: Sample jump index a, weighted by <psi|R_a|psi>
# Need: for each a, compute p_a = <psi| R_a |psi> where
#   R_a = sum_omega rate(omega) * A_{a,omega}^dag A_{a,omega}
# This can be done without precomputing per-jump R matrices by
# computing A_{a,omega}|psi> on the fly (we already do this in the flat scan).

# Option A (recommended): Precompute per-jump R_a matrices in TrajectoryFramework
# Then p_a = dot(psi, R_a * psi) is cheap (one matvec per jump).
# Total: num_jumps matvecs for Stage 1, then ~num_energies matvecs for Stage 2.
# vs current: num_jumps * num_energies matvecs worst case (same or better).

# Option B: Two-pass scan. First pass: scan all (a,omega) to compute cumulative
# per-jump probabilities. Second pass: scan only the selected a's frequencies.
# More complex, same cost.

# Recommendation: Option A. Add R_a::Vector{Matrix{ComplexF64}} to TrajectoryFramework.
```

### TrajectoryFramework Modifications

```julia
struct TrajectoryFramework{T,C,H,PD,D<:AbstractDomain}
    # ... existing fields ...

    # NEW: Per-jump R_a matrices for two-stage sampling
    R_per_jump::Vector{Matrix{T}}   # R_a = sum_omega rate(omega) A_{a,omega}^dag A_{a,omega}

    # ... existing fields ...
end
```

The `build_trajectoryframework` function needs a new loop that computes `R_a` for each jump, analogous to the existing `precompute_R` but without summing across jumps. This is straightforward -- extract the per-jump inner loop of `precompute_R` into a helper.

### Memory and Performance Impact

| Metric | Current | After Fix |
|--------|---------|-----------|
| TrajectoryFramework memory | R (dim x dim) + K0 + U_res + U_B | + num_jumps * (dim x dim) for R_per_jump |
| Per-step cost (expected) | O(num_jumps * num_energies * dim^2) worst case | O(num_jumps * dim^2) for Stage 1 + O(num_energies * dim^2) for Stage 2 |
| Per-step cost (3n Pauli jumps, 4 qubits) | ~12 * ~100 * 256 = 307K mul ops worst case | ~12 * 256 + ~100 * 256 = 28.7K mul ops |

For 4 qubits with 12 Pauli jumps and ~100 energy labels, the two-stage approach is ~10x fewer matrix-vector multiplications per step, because Stage 1 is just `num_jumps` matvecs (using precomputed R_a), not `num_jumps * num_energies`.

## Test Suite Architecture

### Test File Organization

```
test/
+-- runtests.jl                    # Central test runner
+-- test_helpers.jl                # Shared fixtures, reference system construction
+-- test_detailed_balance.jl       # Tier 1: Gibbs fixed point (BohrDomain + B)
+-- test_coherent_term.jl          # Tier 2: B operator consistency across domains
+-- test_oft_consistency.jl        # Tier 3: OFT consistency across domains
+-- test_dm_step_errors.jl         # Tier 4: DM step error scaling (delta)
+-- test_trajectory.jl             # Tier 5: Trajectory correctness (after fix)
+-- test_traj_vs_dm.jl             # Tier 6: Trajectory vs DM cross-validation
```

### Test Helpers: The Shared Fixture Pattern

All tests need the same setup: Hamiltonian, jumps, configs, precomputed data. Extract this into `test_helpers.jl`:

```julia
# test/test_helpers.jl

using QuantumFurnace
using LinearAlgebra, Random, Test

"""Standard 3-qubit Heisenberg test system with known parameters."""
function make_test_system(;
    num_qubits::Int = 3,
    beta::Float64 = 10.0,
    sigma_factor::Float64 = 0.8,  # sigma = sigma_factor / beta
    num_energy_bits::Int = 10,
    w0::Float64 = 0.05,
    a::Float64 = 1/10,
    b::Float64 = 0.4,
    with_coherent::Bool = true,
)
    dim = 2^num_qubits
    sigma = sigma_factor / beta
    w_gamma = 1 / beta
    sigma_gamma = sqrt(2 * w_gamma / beta - sigma^2)
    t0 = 2pi / (2^num_energy_bits * w0)
    num_trotter_steps_per_t0 = 10

    hamiltonian = load_hamiltonian("heis", num_qubits)
    hamiltonian = finalize_hamham(hamiltonian, beta)

    # Jump operators: X, Y, Z on each site
    jump_paulis = [[X], [Y], [Z]]
    num_of_jumps = length(jump_paulis) * num_qubits
    jump_normalization = sqrt(num_of_jumps)
    jumps = JumpOp[]
    for pauli in jump_paulis
        for site in 1:num_qubits
            jump_op = pad_term(pauli, num_qubits, site) / jump_normalization
            jump_op_in_eigenbasis = hamiltonian.eigvecs' * jump_op * hamiltonian.eigvecs
            orthogonal = (jump_op == transpose(jump_op))
            hermitian = (jump_op == jump_op')
            push!(jumps, JumpOp(jump_op, jump_op_in_eigenbasis, orthogonal, hermitian))
        end
    end

    trotter = TrottTrott(hamiltonian, t0, num_trotter_steps_per_t0)

    # Return everything tests might need
    return (; num_qubits, dim, beta, sigma, w_gamma, sigma_gamma,
             w0, t0, a, b, num_energy_bits, num_trotter_steps_per_t0,
             hamiltonian, trotter, jumps, with_coherent, jump_normalization)
end

"""Create a LiouvConfig from the test system parameters."""
function make_liouv_config(sys, domain::AbstractDomain)
    LiouvConfig(
        num_qubits = sys.num_qubits,
        with_coherent = sys.with_coherent,
        with_linear_combination = true,
        domain = domain,
        beta = sys.beta,
        sigma = sys.sigma,
        gaussian_parameters = (sys.w_gamma, sys.sigma_gamma),
        a = sys.a,
        b = sys.b,
        num_energy_bits = sys.num_energy_bits,
        w0 = sys.w0,
        t0 = sys.t0,
        eta = 0.0,
        num_trotter_steps_per_t0 = sys.num_trotter_steps_per_t0,
    )
end

"""Create a ThermalizeConfig from the test system parameters."""
function make_therm_config(sys, domain::AbstractDomain; delta=0.01, mixing_time=10.0)
    ThermalizeConfig(
        num_qubits = sys.num_qubits,
        with_coherent = sys.with_coherent,
        with_linear_combination = true,
        domain = domain,
        beta = sys.beta,
        sigma = sys.sigma,
        gaussian_parameters = (sys.w_gamma, sys.sigma_gamma),
        a = sys.a,
        b = sys.b,
        num_energy_bits = sys.num_energy_bits,
        w0 = sys.w0,
        t0 = sys.t0,
        eta = 0.0,
        num_trotter_steps_per_t0 = sys.num_trotter_steps_per_t0,
        mixing_time = mixing_time,
        delta = delta,
    )
end

# --- Tolerance tiers ---
# These map to the error hierarchy
const TOL_MACHINE = 1e-12    # Machine precision (exact operations)
const TOL_BOHR_EXACT = 1e-10 # Bohr+B should give Gibbs to near machine precision
const TOL_ENERGY_QUAD = 1e-4 # Energy domain introduces Gaussian quadrature error
const TOL_TIME_QUAD = 1e-3   # Time domain adds Riemann sum quadrature error
const TOL_TROTTER = 1e-2     # Trotter adds time discretization error
const TOL_DELTA_SINGLE = 0.1 # Single delta-step error (O(delta^2))
const TOL_TRAJECTORY_STAT = 0.05  # Trajectory statistical error (1/sqrt(ntraj))
```

### Test Tier 1: Detailed Balance (Gibbs Fixed Point)

**Tests:** BohrDomain with coherent term B gives Gibbs state as exact fixed point.

**What this validates:** The fundamental mathematical property -- the Lindbladian with exact KMS detailed balance has the Gibbs state as its kernel.

```julia
# test/test_detailed_balance.jl
@testset "Detailed Balance: Gibbs Fixed Point" begin
    sys = make_test_system(num_qubits=3)

    @testset "Bohr+B: Liouvillian kernel is Gibbs" begin
        config = make_liouv_config(sys, BohrDomain())
        result = run_lindbladian(sys.jumps, config, sys.hamiltonian)

        @test trace_distance_h(
            Hermitian(result.fixed_point),
            sys.hamiltonian.gibbs
        ) < TOL_BOHR_EXACT
    end

    @testset "Bohr+B: DM thermalization converges to Gibbs" begin
        config = make_therm_config(sys, BohrDomain(); delta=0.005, mixing_time=50.0)
        dm0 = Matrix{ComplexF64}(I(sys.dim) / sys.dim)
        result = run_thermalization(sys.jumps, config, dm0, sys.hamiltonian)

        @test result.distances_to_gibbs[end] < TOL_BOHR_EXACT
    end

    @testset "Energy domain: fixed point near Gibbs (quadrature error)" begin
        config = make_liouv_config(sys, EnergyDomain())
        result = run_lindbladian(sys.jumps, config, sys.hamiltonian)

        dist = trace_distance_h(Hermitian(result.fixed_point), sys.hamiltonian.gibbs)
        @test dist < TOL_ENERGY_QUAD
        @test dist > TOL_BOHR_EXACT  # Should NOT be exact -- verifies error is real
    end

    @testset "Error hierarchy: Bohr < Energy < Time < Trotter" begin
        distances = Float64[]
        for domain in [BohrDomain(), EnergyDomain(), TimeDomain(), TrotterDomain()]
            config = make_liouv_config(sys, domain)
            trotter_arg = domain isa TrotterDomain ? sys.trotter : nothing
            result = run_lindbladian(sys.jumps, config, sys.hamiltonian; trotter=trotter_arg)
            push!(distances, trace_distance_h(
                Hermitian(result.fixed_point), sys.hamiltonian.gibbs))
        end

        @test distances[1] < distances[2]  # Bohr < Energy
        @test distances[2] < distances[3]  # Energy < Time
        @test distances[3] < distances[4]  # Time < Trotter
    end
end
```

### Test Tier 2: Coherent Term Consistency

**Tests:** `coherent_bohr()` equals `B_time()` up to time quadrature error, equals `B_trotter()` up to Trotter error.

```julia
# test/test_coherent_term.jl
@testset "Coherent Term B Consistency" begin
    sys = make_test_system(num_qubits=3)
    jump = sys.jumps[1]

    @testset "B_bohr vs B_time (time quadrature error)" begin
        config_bohr = make_liouv_config(sys, BohrDomain())
        config_time = make_liouv_config(sys, TimeDomain())

        pd_bohr = precompute_data(BohrDomain(), config_bohr, sys.hamiltonian)
        pd_time = precompute_data(TimeDomain(), config_time, sys.hamiltonian)

        B_bohr = coherent_bohr(sys.hamiltonian, jump, config_bohr)
        rmul!(B_bohr, pd_bohr.gamma_norm_factor)

        B_t = B_time(jump, sys.hamiltonian, pd_time.b_minus, pd_time.b_plus,
                     sys.t0, sys.beta, sys.sigma)
        rmul!(B_t, pd_time.gamma_norm_factor)

        @test norm(B_bohr - B_t) < TOL_TIME_QUAD
        @test norm(B_bohr - B_t) > TOL_MACHINE  # Not exact (verifies error is real)
    end

    @testset "B_time vs B_trotter (Trotter error)" begin
        config_time = make_liouv_config(sys, TimeDomain())
        config_trotter = make_liouv_config(sys, TrotterDomain())

        pd_time = precompute_data(TimeDomain(), config_time, sys.hamiltonian)
        pd_trotter = precompute_data(TrotterDomain(), config_trotter, sys.trotter)

        B_t = B_time(jump, sys.hamiltonian, pd_time.b_minus, pd_time.b_plus,
                     sys.t0, sys.beta, sys.sigma)
        rmul!(B_t, pd_time.gamma_norm_factor)

        B_tr = B_trotter(jump, sys.trotter, pd_trotter.b_minus, pd_trotter.b_plus,
                         sys.beta, sys.sigma)
        rmul!(B_tr, pd_trotter.gamma_norm_factor)

        # Transform B_trotter to eigenbasis for comparison
        V = sys.hamiltonian.eigvecs' * sys.trotter.eigvecs
        B_tr_in_eigen = V * B_tr * V'

        @test norm(B_t - B_tr_in_eigen) < TOL_TROTTER
    end

    @testset "B is Hermitian" begin
        config = make_liouv_config(sys, BohrDomain())
        pd = precompute_data(BohrDomain(), config, sys.hamiltonian)
        B = coherent_bohr(sys.hamiltonian, jump, config)
        @test norm(B - B') < TOL_MACHINE
    end
end
```

### Test Tier 3: OFT Consistency

**Tests:** `oft!()` (Bohr-exact Gaussian filter) matches `time_oft!()` (Riemann sum) up to time quadrature, matches `trotter_oft!()` up to Trotter error.

```julia
# test/test_oft_consistency.jl
@testset "OFT Consistency Across Domains" begin
    sys = make_test_system(num_qubits=3)
    jump = sys.jumps[1]

    @testset "oft! vs NUFFT prefactors (Time domain)" begin
        config = make_liouv_config(sys, TimeDomain())
        pd = precompute_data(TimeDomain(), config, sys.hamiltonian)

        test_energies = [0.0, 0.1, -0.1]
        for w in test_energies
            if haskey(pd.oft_nufft_prefactors.energy_to_index, w)
                A_exact = oft(jump, w, sys.hamiltonian, sys.sigma)
                pref = prefactor_view(pd.oft_nufft_prefactors, w)
                A_nufft = jump.in_eigenbasis .* pref

                @test norm(A_exact - A_nufft) < TOL_TIME_QUAD
            end
        end
    end
end
```

### Test Tier 4: DM Step Error Scaling

**Tests:** Single delta-step error is O(delta^2), multi-step error is O(delta). This corresponds to Chen Theorem III.1.

```julia
# test/test_dm_step_errors.jl
@testset "DM Step Error Scaling (Chen Theorem III.1)" begin
    sys = make_test_system(num_qubits=3)

    @testset "Single step: error ~ O(delta^2)" begin
        domain = EnergyDomain()
        deltas = [0.1, 0.05, 0.025]
        errors = Float64[]

        for delta in deltas
            config = make_therm_config(sys, domain; delta=delta, mixing_time=delta)
            dm0 = Matrix{ComplexF64}(I(sys.dim) / sys.dim)
            result = run_thermalization(sys.jumps, config, copy(dm0), sys.hamiltonian)

            # Compare single-step DM output against exact channel
            # (exact channel = exp(delta * L) applied to dm0)
            liouv_config = make_liouv_config(sys, domain)
            liouv = construct_lindbladian(sys.jumps, liouv_config, sys.hamiltonian)
            exact_dm_vec = exp(delta * liouv) * vec(dm0)
            exact_dm = reshape(exact_dm_vec, sys.dim, sys.dim)

            push!(errors, norm(result.evolved_dm - exact_dm))
        end

        # Check quadratic scaling: error(delta/2) / error(delta) ~ 1/4
        for i in 1:(length(errors)-1)
            ratio = errors[i+1] / errors[i]
            @test ratio < 0.35  # Should be ~0.25 for O(delta^2)
        end
    end
end
```

### Test Tier 5: Trajectory Correctness (Post-Fix)

**Tests:** After the two-stage sampling fix, trajectory-averaged density matrix matches DM simulation.

```julia
# test/test_trajectory.jl
@testset "Trajectory Simulation" begin
    sys = make_test_system(num_qubits=3)

    @testset "Trajectory average matches DM (EnergyDomain)" begin
        domain = EnergyDomain()
        delta = 0.01
        mixing_time = 20.0
        ntraj = 500  # Enough for statistical convergence

        config = make_therm_config(sys, domain; delta=delta, mixing_time=mixing_time)

        # DM reference
        dm0 = Matrix{ComplexF64}(I(sys.dim) / sys.dim)
        dm_result = run_thermalization(sys.jumps, config, copy(dm0), sys.hamiltonian)

        # Trajectory average
        psi0 = ones(ComplexF64, sys.dim) / sqrt(sys.dim)
        traj_result = run_trajectories(sys.jumps, config, psi0, sys.hamiltonian;
                                       ntraj=ntraj, delta=delta, total_time=mixing_time)

        @test trace_distance_h(
            Hermitian(traj_result.rho_mean),
            Hermitian(dm_result.evolved_dm)
        ) < TOL_TRAJECTORY_STAT
    end

    @testset "Trajectory converges to Gibbs (EnergyDomain)" begin
        domain = EnergyDomain()
        delta = 0.005
        mixing_time = 100.0
        ntraj = 1000

        config = make_therm_config(sys, domain; delta=delta, mixing_time=mixing_time)
        psi0 = ones(ComplexF64, sys.dim) / sqrt(sys.dim)

        result = run_trajectories(sys.jumps, config, psi0, sys.hamiltonian;
                                  ntraj=ntraj, delta=delta, total_time=mixing_time)

        dist = trace_distance_h(Hermitian(result.rho_mean), sys.hamiltonian.gibbs)
        @test dist < TOL_ENERGY_QUAD + TOL_TRAJECTORY_STAT
    end
end
```

### Test Tier 6: Trajectory vs DM Cross-Validation

**Tests:** Same error hierarchy structure in trajectories as in DM mode.

```julia
# test/test_traj_vs_dm.jl
@testset "Trajectory vs DM Cross-Validation" begin
    sys = make_test_system(num_qubits=3)

    @testset "Error hierarchy preserved in trajectories" begin
        # Same test as Tier 1 error hierarchy, but using trajectories
        delta = 0.005
        mixing_time = 100.0
        ntraj = 500
        psi0 = ones(ComplexF64, sys.dim) / sqrt(sys.dim)

        distances = Float64[]
        for domain in [EnergyDomain(), TimeDomain(), TrotterDomain()]
            config = make_therm_config(sys, domain; delta=delta, mixing_time=mixing_time)
            trotter_arg = domain isa TrotterDomain ? sys.trotter : nothing
            result = run_trajectories(sys.jumps, config, psi0, sys.hamiltonian;
                                      ntraj=ntraj, delta=delta, total_time=mixing_time,
                                      trotter=trotter_arg)
            push!(distances, trace_distance_h(
                Hermitian(result.rho_mean), sys.hamiltonian.gibbs))
        end

        @test distances[1] < distances[2] + TOL_TRAJECTORY_STAT  # Energy < Time
        @test distances[2] < distances[3] + TOL_TRAJECTORY_STAT  # Time < Trotter
    end
end
```

## Data Flow: Trajectory Validation Pipeline

### Full Validation Data Flow

```
[make_test_system()]
    |
    +----> HamHam, JumpOps, TrottTrott, parameters
    |
    +----> [make_liouv_config(domain)]
    |         |
    |         v
    |     [run_lindbladian()] ----------> HotSpectralResults
    |         |                              |
    |         v                              v
    |     fixed_point (Liouv kernel)    spectral_gap
    |         |
    |         +----> REFERENCE: trace_distance(fixed_point, gibbs)
    |
    +----> [make_therm_config(domain, delta)]
              |
              +----> [run_thermalization()] --> HotAlgorithmResults
              |         |                          |
              |         v                          v
              |     evolved_dm               distances_to_gibbs[]
              |         |
              |         +----> DM REFERENCE for trajectory comparison
              |
              +----> [run_trajectories()] --> (rho_mean, ...)
                        |                        |
                        v                        v
                    rho_mean (traj avg)     [optional: measurements]
                        |
                        +----> CROSS-VALIDATE: trace_distance(rho_mean, evolved_dm)
                        +----> GIBBS CHECK:    trace_distance(rho_mean, gibbs)
```

### Error Budget Propagation

```
Gibbs state (exact target)
    |
    v
BohrDomain + B          error_1 ~ 0 (machine precision)
    |
    v  + Gaussian filter quadrature (w0 discretization)
EnergyDomain + B        error_2 = error_1 + O(w0 * exp(-sigma^2))
    |
    v  + Riemann sum time quadrature (t0 discretization)
TimeDomain + B           error_3 = error_2 + O(t0 * exp(-sigma^2 * t_max^2))
    |
    v  + Trotter time evolution approximation
TrotterDomain + B        error_4 = error_3 + O(t0 / num_trotter_steps)
    |
    v  + delta-stepping error (weak measurement discretization)
DM step simulation       error_5 = error_4 + O(delta) per mixing time
    |
    v  + statistical sampling error
Trajectory simulation    error_6 = error_5 + O(1/sqrt(ntraj))
```

Each test tier validates one layer of this error budget.

## Component Boundaries for New vs Existing Code

### New Files

| File | Purpose | Dependencies |
|------|---------|--------------|
| `test/runtests.jl` | Central test runner | All test files |
| `test/test_helpers.jl` | Shared fixtures, configs, tolerances | QuantumFurnace module |
| `test/test_detailed_balance.jl` | Tier 1: Gibbs fixed point tests | test_helpers, run_lindbladian, run_thermalization |
| `test/test_coherent_term.jl` | Tier 2: B operator cross-domain tests | test_helpers, coherent_bohr, B_time, B_trotter |
| `test/test_oft_consistency.jl` | Tier 3: OFT cross-domain tests | test_helpers, oft!, nufft prefactors |
| `test/test_dm_step_errors.jl` | Tier 4: DM step error scaling | test_helpers, run_thermalization, construct_lindbladian |
| `test/test_trajectory.jl` | Tier 5: Trajectory correctness | test_helpers, run_trajectories |
| `test/test_traj_vs_dm.jl` | Tier 6: Cross-validation | test_helpers, run_trajectories, run_thermalization |

### Modified Files

| File | What Changes | Why |
|------|-------------|-----|
| `trajectories.jl` | Fix `build_trajectoryframework` scope bug (line 53: `trotter` undefined). Add `R_per_jump` field to `TrajectoryFramework`. Restructure dissipative-jump branch in `step_along_trajectory!` for two-stage sampling. | Core sampling fix. |
| `structs.jl` | No changes needed (TrajectoryFramework is in trajectories.jl, not structs.jl) | -- |

### Unchanged Files (Used As-Is by Tests)

All other source files remain unchanged. Tests exercise them as black boxes through the existing public API.

## Build Order: Dependency-Aware Sequencing

The dependency between the trajectory fix and the test suite determines build order:

```
Phase 0: Test infrastructure (no dependency on trajectory fix)
    |
    +-- test/runtests.jl
    +-- test/test_helpers.jl
    +-- test/test_detailed_balance.jl   (uses only DM / Liouvillian)
    +-- test/test_coherent_term.jl      (uses only B computation)
    +-- test/test_oft_consistency.jl    (uses only OFT)
    +-- test/test_dm_step_errors.jl     (uses only DM stepping)
    |
    All of these can be written and run BEFORE fixing trajectories.
    They validate the known-good DM path and establish reference values.

Phase 1: Fix trajectory sampling (depends on understanding from Phase 0 tests)
    |
    +-- Fix build_trajectoryframework scope bug
    +-- Add R_per_jump precomputation
    +-- Restructure step_along_trajectory! for two-stage sampling
    +-- Quick smoke test: single trajectory doesn't crash
    |
    Run Phase 0 tests to verify nothing is broken.

Phase 2: Trajectory validation tests (depends on Phase 1)
    |
    +-- test/test_trajectory.jl         (trajectory average matches DM)
    +-- test/test_traj_vs_dm.jl         (error hierarchy in trajectories)
    |
    These can only be written after the sampling fix.

Phase 3: Cleanup
    |
    +-- Remove ~200 lines of commented-out code in jump_workers.jl, trajectories.jl
    +-- Remove old KrausFramework, precompute_kraus_jumps, etc.
```

**Why this order:**
1. Phase 0 establishes the DM reference values that Phase 2 will compare against. If DM tests fail, we know the bug is in DM code (not trajectories), which saves debugging time.
2. Phase 1 cannot be properly tested without the DM reference oracle from Phase 0.
3. Phase 2 depends on Phase 1 (fixed sampling) and Phase 0 (reference values).
4. Phase 3 is cosmetic and should only happen after validation confirms correctness.

## Known Bugs to Fix Before Validation

### Bug 1: `build_trajectoryframework` Scope Error (trajectories.jl:53)

```julia
# Line 48: function signature has `where T` but T is not used
function build_trajectoryframework(
    jumps::AbstractVector{<:JumpOp},
    ham_or_trott::Union{HamHam, TrottTrott},
    config::AbstractThermalizeConfig,
    precomputed_data,
    scratch::KrausScratch{ComplexF64},
    delta::Float64) where T    # <-- T is unused, remove `where T`

# Line 53: `trotter` variable is not defined in this scope
    B_total = precompute_coherent_total_B(jumps, ham_or_trott, config, precomputed_data; trotter=trotter)
    #                                                                                    ^^^^^^^
    # `trotter` was never passed in. Should be removed or derived from ham_or_trott.
```

**Fix:** Remove `where T` from the signature. Remove `; trotter=trotter` from the `precompute_coherent_total_B` call -- this function already receives `ham_or_trott` which is either `HamHam` or `TrottTrott`, and `precompute_coherent_total_B` already handles both cases via its domain dispatch (line 23-36 of coherent.jl).

### Bug 2: Coherent Unitary Applied AFTER Branching (trajectories.jl:472-477)

In `step_along_trajectory!`, the coherent unitary `U_B` is applied **after** computing `K0*psi` and the branch probabilities, but **before** the branch is taken. This means the branching probabilities are computed on the pre-U_B state, but the post-branch state includes U_B. The DM version (`jump_contribution!` in jump_workers.jl line 163-169) applies U_B **before** everything else, which is the correct order.

**Fix:** Move the U_B application block to the very beginning of `step_along_trajectory!`, before computing `Rpsi` and `K0*psi`. This matches the DM code.

```julia
# CORRECT ORDER:
# 1. Apply U_B (coherent unitary)
# 2. Compute R*psi, K0*psi, p_nojump, p_res, p_jump_total
# 3. Branch
```

### Bug 3: Coherent Term Variable Name Mismatch (trajectories.jl:96)

```julia
# Line 96: passes B_total but the variable was defined on line 54 as B_total
return TrajectoryFramework(
    ...
    B_total,    # This is the B matrix
    U_B,        # exp(-i delta B_total) -- correct
    ...
)
```

When `config.with_coherent == false`, `B_total` is never assigned (the `if` block is skipped). The `TrajectoryFramework` constructor would receive an undefined variable. Need to add `B_total = nothing` before the `if` block, or restructure.

**Fix:**
```julia
B_total = nothing
if config.with_coherent
    B_total = precompute_coherent_total_B(jumps, ham_or_trott, config, precomputed_data)
    B_total .= 0.5 .* (B_total .+ B_total')
    U_B = exp(-1im * delta * Hermitian(B_total))
else
    U_B = nothing
end
```

## Patterns for This Milestone

### Pattern: Reference Oracle Testing

**What:** Test correctness by comparing a new computation against an established reference implementation at a coarser approximation level.

**When to use:** Every layer of the error hierarchy.

**How it works in QuantumFurnace:**

```
Bohr (exact) --reference-for--> Energy (check: ||Bohr - Energy|| < tol_quadrature)
Energy       --reference-for--> Time   (check: ||Energy - Time||  < tol_time)
Time         --reference-for--> Trotter(check: ||Time - Trotter|| < tol_trotter)
DM mode      --reference-for--> Trajectory (check: ||DM_rho - Traj_rho|| < tol_stat)
```

Each test uses the next-higher level as an oracle. This is possible because the higher-level computations are already validated and working.

### Pattern: Statistical Tolerance with Confidence Bounds

**What:** Trajectory tests are inherently statistical. Use `ntraj` large enough that the standard error `~1/sqrt(ntraj)` is well below the tolerance.

**Rule of thumb:** For tolerance `tol`, use `ntraj >= 4 / tol^2`. For `tol = 0.05`, need `ntraj >= 1600`. For `tol = 0.1`, need `ntraj >= 400`.

**Trade-off:** More trajectories = slower tests. Use smaller `ntraj` (100-500) for quick CI, larger (1000-5000) for thorough validation runs.

### Pattern: Error Monotonicity Assertion

**What:** Assert that errors increase monotonically along the approximation hierarchy. If `error(Energy) > error(Time)`, something is wrong -- the extra approximation in Time should make things worse, not better.

**Implementation:**
```julia
@test dist_bohr < dist_energy < dist_time < dist_trotter
```

This is a powerful diagnostic: it catches both bugs in individual domains AND bugs in the error analysis.

## Anti-Patterns to Avoid

### Anti-Pattern: Testing Against Hardcoded Numerical Values

**What people do:** `@test norm(B_bohr) == 0.0023456789` (hardcoded from a previous run).

**Why it's wrong:** Numerical values depend on the specific Hamiltonian realization (which uses random disorder coefficients via `find_ideal_heisenberg`). If the BSON file changes or the random seed differs, all tests break even though the code is correct.

**Do this instead:** Test *relationships* (e.g., `||B_bohr - B_time|| < tolerance`) and *properties* (e.g., `B == B'` for Hermiticity). Use the deterministic `load_hamiltonian("heis", n)` which loads from a fixed BSON file, not random generation.

### Anti-Pattern: Single Tolerance for All Tests

**What people do:** `const TOL = 1e-8` used everywhere.

**Why it's wrong:** The error hierarchy spans 10+ orders of magnitude. A tolerance tight enough for Bohr (1e-12) will fail for Trotter (1e-2). A tolerance loose enough for Trotter (1e-2) won't catch bugs in Bohr.

**Do this instead:** Use tiered tolerances as shown in `test_helpers.jl` above. Each test tier has a tolerance that matches its expected error level.

### Anti-Pattern: Testing Only Convergence to Gibbs

**What people do:** Only test `trace_distance(rho, gibbs) < tolerance` at the end of a long simulation.

**Why it's wrong:** This is a weak test -- many bugs still converge to something near Gibbs. It doesn't catch: wrong convergence rate, wrong spectral gap, transient errors, or sampling bias that averages out over long time.

**Do this instead:** Test intermediate properties too: error scaling with delta, DM vs trajectory agreement at intermediate times, convergence rate consistency with spectral gap.

## Sources

- Direct codebase analysis of all 23 source files in `src/`
- Direct analysis of all 7 test files in `test/`
- Direct analysis of simulation scripts in `simulations/`
- Direct analysis of playground scripts (especially `unraveling.jl`)
- Chen, Kastoryano, Gilyen (2025) -- KMS detailed balance construction, Theorem III.1 on delta-step errors
- Chen, Kastoryano, Brandao, Gilyen (2023) -- approximate GNS construction, error hierarchy
- Julia Test module documentation -- `@testset`, `@test`, `runtests.jl` standard
- QuantumOptics.jl `test/` structure -- pattern for quantum simulation test suites

---
*Architecture research for: Trajectory validation and test suite integration (v1.0 Trajectories milestone)*
*Researched: 2026-02-13*
