# Feature Research: Trajectory Simulation Validation and Testing

**Domain:** Quantum trajectory (MCWF) validation for Gibbs sampling Lindbladians
**Researched:** 2026-02-13
**Confidence:** HIGH (codebase analysis + Chen et al. papers + established MCWF literature + QuantumOptics.jl/QuantumToolbox.jl reference implementations)

## Scope

This research covers features needed for the **v1.0 Trajectories milestone**: fixing the trajectory sampling bug, validating trajectory averages against density matrix evolution, and building a comprehensive correctness test suite. It does NOT cover general package features (CI, docs, Hamiltonians, circuit generation) already addressed in previous project-level research.

---

## Feature Landscape

### Table Stakes (Must Have for Milestone Completion)

Features without which the trajectory implementation cannot be trusted or published.

| Feature | Why Expected | Complexity | Dependencies |
|---------|--------------|------------|--------------|
| **Two-stage jump sampling fix** | Current code samples flat from (a, omega) Kraus pair pool. Physically, one must first pick Lindblad operator A^a with probability proportional to its total jump rate, then sample omega within that operator. Flat sampling produces incorrect jump statistics -- the density matrix average will not converge to the master equation solution. | MEDIUM | Existing `step_along_trajectory!` in `trajectories.jl` |
| **Single-step CPTP verification** | Each delta-step must be a valid CPTP map. The weak measurement scheme (Chen CKBG23 Theorem III.1) guarantees that sum(K_i^dag K_i) = I to O(delta^2). Must verify: `K0^dag K0 + delta*R + U_res^dag U_res = I` numerically. Without this, accumulating steps amplifies unitarity violations. | LOW | `build_trajectoryframework`, `precompute_R` |
| **Trajectory-averaged rho vs DM evolution** | The defining correctness criterion: rho_traj = (1/N) sum_n psi_n psi_n^dag must converge to rho_DM from `run_thermalize()` for the same parameters. Statistical error scales as 1/sqrt(N_traj). This is what "validation" means for trajectory simulation. | MEDIUM-HIGH | Two-stage sampling fix, `run_thermalize()`, `run_trajectories()` |
| **Trace distance and fidelity between trajectory rho and DM rho** | Standard metrics for comparing quantum states. Trace distance T(rho, sigma) = 0.5 * ||rho - sigma||_1 and fidelity F(rho, sigma) = (tr sqrt(sqrt(rho) sigma sqrt(rho)))^2. Both are needed to quantify agreement. `qi_tools.jl` already has `trace_distance_h` and `fidelity`. | LOW | Existing `qi_tools.jl` functions |
| **Detailed balance preservation test (DM)** | For the KMS construction (with_coherent=true), the Lindbladian satisfies exact detailed balance: L_beta^dag[.] = rho_beta^(1/4) L_beta[rho_beta^(-1/4) . rho_beta^(-1/4)] rho_beta^(1/4). The fixed point of the Liouvillian must be exactly the Gibbs state. Test: `trace_distance(fixed_point, gibbs) < tol` for BohrDomain. For approximate domains, the error should match the predicted error budget. | LOW | `run_lindbladian`, `gibbs_state` |
| **Domain error hierarchy test** | The core theoretical structure: Bohr (exact DB) -> Energy (+quadrature error) -> Time (+time quadrature error) -> Trotter (+Trotter error). Must verify that errors are monotonically non-decreasing along this chain: `dist_bohr <= dist_energy <= dist_time <= dist_trotter`. | MEDIUM | All four domains working, `run_lindbladian` or `run_thermalize` for each |
| **Single-step DM error scaling test** | Per Chen CKBG23 Theorem III.1: the single-step weak measurement scheme produces a channel O(delta^2)-close to (I + delta*L) in diamond norm. Must verify empirically that single-step DM error scales as delta^2 (not delta). | LOW | `jump_contribution!` for thermalize configs, single step application |
| **Multi-step DM error scaling test** | Over M = t/delta steps, the accumulated error is O(delta) = O(t^2/M) in diamond norm. Chen Theorem III.1 gives t^2/epsilon cost. Verify: doubling steps (halving delta) at fixed total_time halves the total DM error to Gibbs. | LOW-MEDIUM | `run_thermalize` with varying delta |
| **Statistical convergence test** | Verify that the standard error of the trajectory ensemble mean decreases as 1/sqrt(N_traj). Run N trajectories, compute trace distance to DM result, check that doubling N halves the standard error. | MEDIUM | Fixed trajectory sampling, multiple runs |
| **Coherent term (U_B) correctness for trajectories** | The deterministic coherent unitary U_B = exp(-i delta B_total) must be applied at each step. For with_coherent=true, omitting U_B breaks detailed balance. Test: with_coherent=true TrotterDomain should reach <= 1e-6 distance to Gibbs; without coherent term, it should be much worse. | LOW | `build_trajectoryframework`, coherent term computation |
| **build_trajectoryframework bug fixes** | The current `build_trajectoryframework` has a bug: the `trotter` variable is referenced but not in scope (line 53: `trotter=trotter`). Also, `B_total` is used before being assigned when `with_coherent=false`. These are compilation errors that prevent any trajectory simulation. | LOW | Direct code fixes in `trajectories.jl` |

### Differentiators (Valuable but Not Blocking)

Features that strengthen the validation but are not strictly required for milestone sign-off.

| Feature | Value Proposition | Complexity | Dependencies |
|---------|-------------------|------------|--------------|
| **OFT consistency check across domains** | Verify that the operator Fourier transform A_omega computed via Gaussian filter (EnergyDomain) and via NUFFT prefactors (Time/TrotterDomain) agree to quadrature precision for each jump operator. This catches NUFFT bugs early. | LOW | `oft!`, `prefactor_view`, OFT NUFFT machinery |
| **R matrix consistency check** | The reflection operator R = sum_k L_k^dag L_k must be the same whether computed by the Liouvillian pathway or the trajectory pathway (`precompute_R`). Cross-validate for each domain. | LOW | `precompute_R` in `trajectories.jl`, R computation in `jump_workers.jl` |
| **Per-observable trajectory convergence** | Beyond the full density matrix comparison: track individual observable expectations <O_i>(t) from trajectories vs DM, e.g. local magnetization <Z_i>. The `run_trajectories` function already supports `observables` parameter but it is untested. | LOW-MEDIUM | Observable measurement machinery in `run_trajectories` |
| **Jump statistics histogram** | Record which (a, omega) outcomes are selected across many trajectories. Verify that the empirical distribution matches the theoretical rates: P(a, omega) proportional to delta * rate^2(omega) * ||A_omega psi||^2. This directly tests the two-stage sampling correctness. | MEDIUM | Two-stage sampling fix, histogram analysis code |
| **Regression test with frozen reference data** | Compute trajectory-averaged rho for a fixed (seed, parameters, system) and store the result. Future code changes must reproduce this within tolerance. Catches silent regressions. | LOW | One successful trajectory validation run |
| **Normalization drift monitoring** | Track ||psi||^2 across trajectory steps. After each step, psi is re-normalized, but the pre-normalization norm reveals whether the CPTP map is well-conditioned. Large deviations (||psi||^2 far from 1) indicate numerical issues. | LOW | Minor instrumentation in `step_along_trajectory!` |
| **Multi-threaded trajectory validation** | Verify that multi-threaded trajectory accumulation produces the same rho_mean as serial execution (up to statistical noise). Catches race conditions in shared workspace. | MEDIUM | Thread-safe trajectory framework, `Base.Threads` |
| **Diamond norm channel comparison** | The strongest metric: compare the DM channel E_DM (one thermalization step) to the trajectory-averaged channel E_traj as superoperators using diamond norm. More rigorous than comparing output states for a single input. But diamond norm computation is expensive (SDP). | HIGH | Channel tomography (run DM step on basis states), SDP solver (SCS.jl or similar) |
| **Confidence interval for trace distance** | Report statistical uncertainty on the trajectory-vs-DM comparison. Use bootstrap resampling or the analytic 1/sqrt(N) bound on the ensemble mean. Enables quantitative statements like "trajectory rho agrees with DM rho to trace distance 0.01 +/- 0.003". | LOW-MEDIUM | Multiple trajectory runs, basic statistics |

### Anti-Features (Do Not Build for This Milestone)

| Feature | Why It Seems Useful | Why Problematic for This Milestone | What to Do Instead |
|---------|--------------------|------------------------------------|-------------------|
| **Continuous-time MCWF (adaptive timestep)** | Standard in QuantumOptics.jl and QuantumToolbox.jl. Detects jump times exactly instead of using fixed delta steps. | QuantumFurnace's algorithm is inherently discrete-time: the CPTP map is designed for fixed delta-steps per Chen's weak measurement scheme. Adaptive timestep would require a fundamentally different formulation (continuous H_eff propagation + jump detection), which contradicts the theoretical framework. | Keep the fixed-delta-step structure. It is the physically correct unraveling of the Chen channel, not an approximation of continuous dynamics. |
| **Bohr domain trajectories** | "Why not validate trajectories in the exact (Bohr) domain too?" | The BohrDomain Lindbladian uses a non-standard dissipative structure (the alpha(bohr_freqs, nu) weighting creates cross-frequency coupling). Diagonalizing the resulting Kossakowski matrix to extract independent Kraus operators is computationally prohibitive for >3 qubits. The project context explicitly excludes Bohr for trajectories. | Validate Bohr only via DM methods (`run_lindbladian`, `run_thermalize`). Use Energy/Time/Trotter for trajectory validation. |
| **Quantum state tomography from trajectories** | "Reconstruct rho from measurement outcomes, not from psi*psi^dag averaging." | This conflates two different things. Trajectory averaging of |psi><psi| is the exact unraveling -- it should give rho directly. Tomography from measurement outcomes is a much harder problem (needs many measurement bases, maximum likelihood estimation). It is a different validation target. | Use direct psi*psi^dag averaging. It is both correct and computationally cheaper. |
| **Large system trajectory benchmarks (>6 qubits)** | "Show it works at scale." | Validation requires comparing against DM results, which means computing the full density matrix anyway. For >6 qubits the DM is 64x64 to 4096x4096 -- still feasible but slow. The point of this milestone is correctness, not performance. Use small systems (2-4 qubits) where exact comparison is fast. | Validate on 2-4 qubit Heisenberg. Add 5-6 qubit tests as regression tests with looser tolerances. Performance benchmarks belong in a later milestone. |
| **Non-Hermitian jump operator support in trajectories** | The JumpOp struct has a `hermitian` flag. Non-Hermitian jumps follow a different code path. | The Heisenberg model with single-site Pauli jump operators (X, Y, Z) uses exclusively Hermitian jumps. Testing non-Hermitian paths requires constructing artificial systems. It is lower priority than validating the Hermitian path that all target systems use. | Test Hermitian path thoroughly first. Add non-Hermitian tests in a follow-up if needed for Ding et al. (2024) extensions. |

---

## Feature Details

### Two-Stage Jump Sampling (Critical Bug Fix)

**Current behavior (WRONG):** The `step_along_trajectory!` function iterates over all `(jump, omega)` pairs in a flat cumulative scan:
```
for jump in fw.jumps
    for w in energy_labels
        p = delta * rate2(w) * ||A_omega psi||^2
        csum += p
        if csum >= target: select this outcome
```
This samples from the joint distribution P(a, omega) = delta * rate2(omega) * ||A^a_omega psi||^2 with a single flat scan. While this may appear correct at first glance, it conflates the hierarchical structure.

**Required behavior (CORRECT):** The Chen algorithm defines Lindblad operators {sqrt(gamma(omega)) * A^a(omega)} indexed by both `a` and `omega`. The theoretically correct unraveling samples:
1. **Stage 1:** Pick which Lindblad operator class fires, weighted by total rate per operator:
   P(a) = sum_omega [delta * rate2(omega) * ||A^a_omega psi||^2] / p_jump_total
2. **Stage 2:** Given the selected operator A^a, pick omega with probability:
   P(omega | a) = rate2(omega) * ||A^a_omega psi||^2 / sum_omega' [rate2(omega') * ||A^a_omega' psi||^2]

**Why flat scan is actually equivalent:** After careful analysis, the flat cumulative scan IS mathematically correct for sampling from the joint distribution. The joint probability P(a, omega) = delta * rate2(omega) * ||A^a_omega psi||^2 can be sampled either hierarchically (two-stage) or via a single flat scan -- both produce the same distribution. The cumulative scan in the current code correctly weights each outcome.

**The ACTUAL bug:** The real sampling bug is more subtle. Looking at the code carefully:
- The code samples a single random number `r = rand() * total_weight` to decide between no-jump/residual/dissipative-jump
- For the dissipative-jump branch, it uses `target = r - p_nojump - p_res` and does a cumulative scan
- This is correct IF the probabilities sum correctly to `p_jump_total`

The bug may instead be in the probability normalization: `total_weight = p_nojump + p_res + p_jump_total` should equal 1.0 for a proper CPTP map (up to O(delta^2) from the residual). If it does not, the sampling is biased. **This must be verified numerically as the first validation step.**

**Confidence:** HIGH that a bug exists (project context states it explicitly), MEDIUM on the exact nature of the bug. The first task must be numerical CPTP verification.

### CPTP Verification Tests

The weak measurement scheme produces a CPTP map with Kraus operators:
- **K0** = I - alpha * R, where alpha = 1 - sqrt(1 - delta)
- **K_{a,omega}** = sqrt(delta * rate2(omega)) * A^a_omega (jump Kraus operators)
- **K_res** = U_residual (Cholesky factor of residual S = I - K0^dag K0 - delta*R)

The CPTP condition requires: K0^dag K0 + sum_{a,omega} K_{a,omega}^dag K_{a,omega} + K_res^dag K_res = I

This decomposes as:
- K0^dag K0 = I - 2*alpha*R + alpha^2 * R^2
- sum K_{a,omega}^dag K_{a,omega} = delta * R (by construction of R)
- K_res^dag K_res = S = (2*alpha - delta)*R - alpha^2 * R^2

Sum = I - 2*alpha*R + alpha^2*R^2 + delta*R + (2*alpha - delta)*R - alpha^2*R^2 = I

This is exact algebraically. Numerical deviations come from:
1. Cholesky factorization of S (S may have tiny negative eigenvalues from roundoff)
2. The epsilon shift applied to S before Cholesky (10*eps(Float64) per diagonal element)
3. Hermitianization of R

**Test:** Compute K0^dag K0 + delta*R + U_res^dag U_res and verify ||result - I|| < tolerance (expect < 1e-12 for well-conditioned systems).

### Trajectory vs DM Convergence Protocol

**Standard MCWF validation procedure** (verified across QuantumOptics.jl, QuantumToolbox.jl, QDYN):

1. **Fix parameters:** system (Heisenberg 3-4 qubit), domain (Energy or Time), delta, total_time, initial state (e.g., |0...0> or maximally mixed)
2. **Run DM evolution:** `run_thermalize()` to get rho_DM at final time
3. **Run N trajectories:** `run_trajectories(ntraj=N)` to get rho_traj = (1/N) sum |psi_n><psi_n|
4. **Compute discrepancy:** D = trace_distance(rho_traj, rho_DM)
5. **Verify scaling:** Repeat for N = 100, 500, 1000, 5000. Verify D decreases as approximately 1/sqrt(N).

**Expected tolerances** (based on domain and parameters):

| Domain | with_coherent | Expected trace_distance to Gibbs | Trajectory tolerance (N=1000) |
|--------|--------------|----------------------------------|-------------------------------|
| Energy | false (GNS) | O(quadrature error) ~ 1e-3 to 1e-5 | Should match DM to within ~0.03 |
| Energy | true (KMS) | O(quadrature error) ~ 1e-3 to 1e-5 | Should match DM to within ~0.03 |
| Time | true (KMS) | O(quad + time quad) ~ 1e-3 to 1e-4 | Should match DM to within ~0.03 |
| Trotter | true (KMS) | <= 1e-6 (with good Trotter params) | Should match DM to within ~0.03 |

The trajectory vs DM discrepancy should be dominated by statistical noise (1/sqrt(N)), NOT by systematic bias. If a systematic bias remains as N grows, the sampling or CPTP construction is wrong.

### Error Budget Decomposition

The total error from the Gibbs state has five independent contributions:

| Error Source | Scaling | Controlled By | How to Isolate |
|-------------|---------|---------------|----------------|
| **Detailed balance violation** (GNS only) | depends on sigma_E | sigma parameter | Compare KMS (with_coherent=true) vs GNS: only difference is B term |
| **Energy quadrature error** | O(1/N_omega) ~ O(w0) | w0, num_energy_bits | Compare Bohr vs Energy at same parameters |
| **Time quadrature error** | O(1/N_t) ~ O(t0) | t0, num_energy_bits | Compare Energy vs Time at same parameters |
| **Trotter error** | O(t0/num_trotter_steps)^2 | num_trotter_steps_per_t0 | Compare Time vs Trotter at same parameters |
| **Delta-step discretization** | O(delta) per total evolution | delta | Fix system, vary delta, measure distance to continuous-time Gibbs |
| **Trajectory sampling noise** | O(1/sqrt(N_traj)) | N_traj | Fix everything else, vary N_traj |

**Key insight for testing:** To isolate the trajectory sampling correctness, compare trajectory-averaged rho against the DM-evolved rho (not against the Gibbs state directly). Both DM and trajectory share the same quadrature/Trotter/delta errors. The only difference should be statistical noise from trajectory sampling.

### Statistical Methods for Validation

**Confidence:** HIGH -- these are standard methods from quantum state estimation literature.

| Method | What It Tests | When to Use | Implementation Cost |
|--------|-------------|-------------|---------------------|
| **Trace distance D(rho_traj, rho_DM)** | Overall state agreement | Primary metric for every validation test | LOW (already in `qi_tools.jl`) |
| **Fidelity F(rho_traj, rho_DM)** | Overlap between states; more sensitive near pure states | Complementary to trace distance | LOW (already in `qi_tools.jl`) |
| **Frobenius norm ||rho_traj - rho_DM||_F** | Cheap proxy for trace distance; faster to compute (no eigendecomposition) | Quick sanity checks, per-step monitoring | LOW (trivial with `norm()`) |
| **1/sqrt(N) convergence fit** | Verifies trajectory averaging is unbiased | After fixing sampling bug -- confirms fix is correct | LOW (run at several N values, fit power law) |
| **Bootstrap confidence intervals** | Uncertainty quantification on trace distance | When reporting quantitative agreement numbers | MEDIUM (resample subsets of trajectories, compute CI) |
| **Per-element comparison max_ij |rho_traj_ij - rho_DM_ij|** | Catches localized errors that global metrics might average out | Debugging when global metrics look fine but physics is wrong | LOW |
| **Eigenvalue spectrum comparison** | Verifies correct population distribution | Confirms Gibbs-like population ordering | LOW |
| **Chi-squared goodness of fit** | Tests if observed jump statistics match predicted rates | Validating two-stage sampling correctness | MEDIUM (requires collecting jump histograms) |

---

## Feature Dependencies

```
[build_trajectoryframework bug fix]
    |
    v
[CPTP verification test]
    |
    v
[Two-stage sampling fix (or verification that flat scan is correct)]
    |
    v
[Single-step DM error scaling test (delta^2)]
    |
    v
[Trajectory vs DM convergence] ----requires----> [run_thermalize() working for Energy/Time/Trotter]
    |                                                      |
    |                                                      v
    |                                             [Detailed balance test (DM only)]
    v                                                      |
[Statistical convergence test (1/sqrt(N))]                 v
    |                                             [Domain error hierarchy test]
    v                                                      |
[Per-observable trajectory convergence]                    v
    |                                             [Multi-step DM error scaling test]
    v
[Regression tests with frozen data]
    |
    v
[OFT consistency check]  [R matrix consistency check]  [Jump statistics histogram]
```

### Dependency Notes

- **build_trajectoryframework fix is the gate:** Nothing trajectory-related works until compilation errors in `build_trajectoryframework` are fixed (undefined `trotter` variable on line 53, `B_total` referenced when `with_coherent=false`).
- **CPTP verification before trajectory runs:** Must confirm the Kraus map is trace-preserving before trusting any trajectory output. Otherwise, non-unit total probability will bias all downstream results.
- **DM tests are independent of trajectory tests:** Detailed balance, domain hierarchy, and DM error scaling tests depend only on `run_lindbladian` and `run_thermalize`, which already work. These can proceed in parallel with trajectory bug fixes.
- **Trajectory vs DM is the milestone gate:** The entire milestone succeeds or fails on whether rho_traj matches rho_DM. Everything else supports this.
- **Statistical convergence proves the fix:** If the 1/sqrt(N) scaling holds, the sampling is provably unbiased. If it plateaus, there is a systematic error.

---

## MVP Definition (Milestone v1.0 Trajectories)

### Must Complete

- [ ] **Fix `build_trajectoryframework` compilation bugs** -- Cannot run any trajectory code without this.
- [ ] **CPTP verification test** -- Confirm K0^dag K0 + delta*R + U_res^dag U_res = I to machine precision.
- [ ] **Fix or verify two-stage jump sampling** -- Either implement proper two-stage sampling or prove the flat scan is equivalent and find the actual bug.
- [ ] **Trajectory vs DM validation for EnergyDomain** -- rho_traj matches rho_DM to within statistical noise for 3-qubit Heisenberg.
- [ ] **Trajectory vs DM validation for TimeDomain** -- Same test for TimeDomain.
- [ ] **Trajectory vs DM validation for TrotterDomain** -- Same test for TrotterDomain (with_coherent=true).
- [ ] **Domain error hierarchy test (DM)** -- Bohr <= Energy <= Time <= Trotter error chain verified.
- [ ] **Single-step delta^2 error scaling test** -- Empirical verification of Chen Theorem III.1.
- [ ] **Multi-step delta error scaling test** -- Verify error accumulation is O(delta) over full evolution.
- [ ] **Detailed balance test** -- BohrDomain fixed point is Gibbs to machine precision; TrotterDomain with_coherent=true reaches <= 1e-6.
- [ ] **Coherent term correctness test** -- With vs without coherent term: with_coherent=true must be dramatically closer to Gibbs.
- [ ] **Statistical 1/sqrt(N) convergence test** -- Trajectory averaging error decreases as expected with N.

### Add After Core Validation

- [ ] **Per-observable trajectory convergence** -- <Z_i>(t) from trajectories matches DM.
- [ ] **Jump statistics histogram** -- Empirical jump rates match theoretical predictions.
- [ ] **OFT consistency across domains** -- Energy OFT matches Time/Trotter NUFFT OFT.
- [ ] **R matrix cross-validation** -- Trajectory R matches Liouvillian R.
- [ ] **Regression tests with frozen data** -- Lock down a reference result.
- [ ] **Confidence interval reporting** -- Bootstrap CIs on trajectory-vs-DM trace distance.
- [ ] **Normalization drift monitoring** -- Track pre-normalization ||psi||^2 per step.

### Defer

- [ ] **Diamond norm channel comparison** -- Expensive (SDP), only needed for publication rigor.
- [ ] **Multi-threaded trajectory validation** -- Correctness first, performance later.
- [ ] **Non-Hermitian jump operator tests** -- All target systems use Hermitian jumps.
- [ ] **Large system benchmarks (>6 qubits)** -- Validation is about correctness, not scale.

---

## Feature Prioritization Matrix

| Feature | Correctness Value | Implementation Cost | Priority |
|---------|------------------|---------------------|----------|
| build_trajectoryframework fix | CRITICAL | LOW | P0 |
| CPTP verification test | CRITICAL | LOW | P0 |
| Two-stage sampling fix/verify | CRITICAL | MEDIUM | P0 |
| Trajectory vs DM (Energy) | CRITICAL | MEDIUM | P1 |
| Trajectory vs DM (Time) | CRITICAL | MEDIUM | P1 |
| Trajectory vs DM (Trotter) | CRITICAL | MEDIUM | P1 |
| Domain error hierarchy test | HIGH | LOW | P1 |
| Single-step delta^2 test | HIGH | LOW | P1 |
| Multi-step delta test | HIGH | LOW | P1 |
| Detailed balance test | HIGH | LOW | P1 |
| Coherent term test | HIGH | LOW | P1 |
| 1/sqrt(N) convergence test | HIGH | LOW-MEDIUM | P1 |
| Per-observable convergence | MEDIUM | LOW-MEDIUM | P2 |
| Jump statistics histogram | MEDIUM | MEDIUM | P2 |
| OFT consistency check | MEDIUM | LOW | P2 |
| R matrix cross-validation | MEDIUM | LOW | P2 |
| Regression tests | MEDIUM | LOW | P2 |
| Confidence intervals | LOW-MEDIUM | LOW-MEDIUM | P2 |
| Normalization drift | LOW | LOW | P2 |
| Diamond norm comparison | LOW | HIGH | P3 |
| Multi-threaded validation | LOW | MEDIUM | P3 |
| Non-Hermitian jump tests | LOW | MEDIUM | P3 |

**Priority key:**
- P0: Blocks all other work; fix immediately
- P1: Core validation; milestone success criteria
- P2: Strengthens validation; add once core passes
- P3: Nice to have; defer to later milestone

---

## Test System Specifications

### Recommended Test Systems

| System | Qubits | Hilbert Space dim | DM size | Why Use It |
|--------|--------|-------------------|---------|------------|
| **1D Heisenberg, n=2** | 2 | 4 | 4x4 | Fastest, good for debugging. Known exact spectrum. |
| **1D Heisenberg, n=3** | 3 | 8 | 8x8 | Primary validation target. Rich enough spectrum, fast enough for many trajectories. |
| **1D Heisenberg, n=4** | 4 | 16 | 16x16 | Secondary validation. Approaches "real" system complexity. 12 jump operators (3 Paulis x 4 sites). |
| **1D Heisenberg, n=4 + Z-disorder** | 4 | 16 | 16x16 | Breaks symmetries. Non-degenerate spectrum ensures cleaner Bohr frequencies. |

### Recommended Parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| beta | 1.0, 5.0, 10.0 | Low/medium/high temperature. High beta = harder thermalization, tests convergence. |
| sigma | 1/beta | Standard choice from Chen papers. |
| delta | 0.01, 0.05, 0.1 | Small enough for O(delta^2) per step to be negligible. 0.1 is aggressive; 0.01 is safe. |
| num_energy_bits | 10-12 | 2^10 = 1024 to 2^12 = 4096 energy grid points. Sufficient for 4-qubit systems. |
| w0 | 0.05 | Fine enough energy resolution for Heisenberg spectrum. |
| N_traj | 100, 500, 1000, 5000 | For convergence scaling tests. 1000 is the sweet spot for validation. |
| total_time | 5-20 * (1/spectral_gap) | Must be long enough to reach steady state. |
| with_coherent | true | Primary interest is KMS-DB. Test both true/false for comparison. |
| with_linear_combination | true | Metropolis-style transition weight for faster mixing. |
| a, b | beta/30, 0.4 | From existing test configs in codebase. |

---

## Comparison: Trajectory Validation Approaches

| Approach | QuantumOptics.jl | QuantumToolbox.jl | QuantumFurnace Needs |
|----------|-----------------|-------------------|---------------------|
| **Trajectory method** | Continuous-time MCWF with adaptive timestep | Continuous-time MCWF (mcsolve) | Fixed-delta discrete CPTP map unraveling |
| **Jump detection** | Detect when norm drops below random threshold | Same as QuantumOptics.jl | No detection needed -- fixed timestep with probability branching |
| **Jump selection** | P(L_n) ~ <psi|L_n^dag L_n|psi> | Same | P(a, omega) ~ delta * rate2(omega) * ||A^a_omega psi||^2 |
| **Ensemble averaging** | Direct psi*psi^dag sum | Direct psi*psi^dag sum | Same (already in `_accumulate_density_matrix!`) |
| **Validation target** | master equation rho(t) | master equation rho(t) | `run_thermalize` DM rho(t) |
| **Error metric** | Visual comparison, no automated tests | Visual comparison (plot convergence) | Automated trace distance + fidelity tests with tolerances |
| **Unique challenge** | None (standard Lindbladian) | None | Residual Kraus K_res for CPTP completion; coherent U_B unitary; error budget with four domain levels |

The key difference: QuantumOptics.jl and QuantumToolbox.jl validate continuous-time MCWF against continuous-time master equation. QuantumFurnace validates a discrete-CPTP-map unraveling against the same discrete-CPTP-map applied to the density matrix. This is simpler in some ways (no ODE solver needed) but more complex in others (the CPTP map itself has internal structure with K0, K_res, and many jump Kraus operators).

---

## Sources

### Verified (HIGH confidence)
- Chen, Kastoryano, Brandao, Gilyen (2023) "Quantum Thermal State Preparation" [arXiv:2303.18224](https://arxiv.org/abs/2303.18224) -- Theorem III.1 (weak measurement scheme, O(delta^2) per step, O(t^2/epsilon) cost), Corollary III.1 (random jump selection). Direct reading of paper text.
- Chen, Kastoryano, Gilyen (2025) "An efficient and exact noncommutative quantum Gibbs sampler" [arXiv:2311.09207](https://arxiv.org/abs/2311.09207) -- KMS detailed balance construction, coherent term B, Lindbladian fixed point stationarity. Direct reading.
- [QDYN documentation: Quantum Jump Method](https://ag-koch.gitpages.physik.fu-berlin.de/qdyn/main/concepts/mcwf.html) -- Ensemble averaging rho = (1/N) sum |psi_n><psi_n|, statistical error 1/sqrt(N), jump operator selection probabilities.
- [QuantumOptics.jl: Quantum trajectories](https://docs.qojulia.org/timeevolution/mcwf/) -- Standard MCWF validation methodology.
- [QuantumToolbox.jl: Monte Carlo Solver](https://qutip.org/QuantumToolbox.jl/stable/users_guide/time_evolution/mcsolve) -- Jump selection via cumulative probability, convergence with trajectory count.
- QuantumFurnace.jl codebase -- Direct analysis of `trajectories.jl`, `jump_workers.jl`, `furnace.jl`, `qi_tools.jl`, `structs.jl`, `errors.jl`.

### Partially verified (MEDIUM confidence)
- [Abdelhafez et al. (2019) "The Monte Carlo wave-function method: A robust adaptive algorithm and a study in convergence"](https://arxiv.org/abs/1803.08589) -- Convergence properties, discretization error analysis, adaptive timestep (not directly applicable but informs error understanding).
- [Rall et al. (2025) "A Randomized Method for Simulating Lindblad Equations and Thermal State Preparation"](https://arxiv.org/html/2407.06594) -- Per-step O(lambda^2 tau^2) error, multi-step O(T^2/M) average error, O(T/sqrt(M)) random channel error, exponential convergence to thermal state preservation.
- [Ding, Li, Lin (2024) "Efficient quantum Gibbs samplers with KMS detailed balance"](https://arxiv.org/abs/2404.05998) -- Alternative KMS construction with discrete jump operators, simplified error analysis.

### Domain knowledge (HIGH confidence, not web-verified)
- Trace distance and fidelity as standard quantum state comparison metrics: Fuchs-van de Graaf inequalities relating them.
- 1/sqrt(N) convergence of Monte Carlo ensemble averages: standard Central Limit Theorem result.
- CPTP map completeness relation sum K_i^dag K_i = I: fundamental requirement for physical quantum channels.
- Bootstrap confidence interval methodology: standard statistical technique.

---
*Feature research for: Trajectory simulation validation and testing (v1.0 Trajectories milestone)*
*Researched: 2026-02-13*
