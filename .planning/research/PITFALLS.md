# Pitfalls Research: Trajectory Simulation Validation and Testing

**Domain:** Quantum trajectory (MCWF) validation against density matrix evolution; correctness testing of numerical quantum simulation with multiple error sources
**Researched:** 2026-02-13
**Confidence:** HIGH (codebase-grounded + literature-verified)

---

## Critical Pitfalls

### Pitfall 1: Flat (a, omega) Sampling Instead of Two-Stage Jump Selection Produces Wrong Steady State

**What goes wrong:**
The current trajectory code in `step_along_trajectory!` (`trajectories.jl` lines 493-568 for Time/Trotter, 651-727 for Energy) iterates over all `(jump, omega)` pairs in a single flat cumulative-probability scan. This means the probability of selecting any particular jump outcome is `p(a,w) = delta * rate2(w) * ||A_{a,w} psi||^2`. For the Chen 2025 construction, the correct physical sampling should first select which Lindblad operator index `a` with probability proportional to `sum_w rate2(w) * ||A_{a,w} psi||^2`, then sample `omega` conditioned on `a`. When the jump operators have very different total rates (e.g., Pauli X on site 1 vs site 4 in a non-translationally-invariant chain), flat sampling biases toward operators with many high-rate frequency channels, regardless of their physical coupling to the current state.

**Why it happens:**
For the specific case of single-site Paulis on a Heisenberg chain where all jumps are normalized by `1/sqrt(num_jumps)`, flat and two-stage sampling give the same answer because the per-jump normalization is uniform. The bug only manifests when: (a) jump operators have different norms or couplings, (b) the Hamiltonian breaks translational symmetry, or (c) the gamma_norm_factor differs per jump. In the current codebase, all jumps share one `gamma_norm_factor` and are uniformly normalized, so the bug is latent.

**How to avoid:**
- Implement two-stage sampling: first accumulate `p_a = delta * sum_w rate2(w) * ||A_{a,w} psi||^2` per jump `a`, pick `a` from this distribution, then scan `omega` only within the selected `a`.
- Even if the current code gives correct results for uniform jumps, the fix should be done now because: (1) future Hamiltonians (Ising, 2D lattice) will break the symmetry, (2) the Ding et al. construction uses non-uniform discrete jumps, (3) it is the mathematically correct unraveling.
- Write a regression test with deliberately non-uniform jump norms and verify the trajectory steady state matches DM.

**Warning signs:**
- Trajectory convergence to wrong state when using non-uniform jump operator norms.
- Discrepancy between trajectory and DM results that does NOT shrink with more trajectories (systematic bias, not sampling noise).
- Works perfectly for symmetric Heisenberg + uniform Paulis but fails for any asymmetric case.

**Phase to address:**
Phase 1 (Trajectory Sampling Fix) -- this is the primary known bug.

---

### Pitfall 2: Comparing Trajectory Average to DM Without Accounting for the Error Budget

**What goes wrong:**
When cross-validating trajectories against density matrix evolution, the test checks `||rho_traj - rho_DM|| < tolerance`. If this fails, the developer assumes the trajectory code is buggy. But the actual error has multiple independent sources:

```
||rho_traj - gibbs|| <= (sampling noise: O(1/sqrt(ntraj)))
                      + (delta-step discretization: O(delta) per step, O(delta) total for Chen's scheme)
                      + (quadrature error from energy integral approximation)
                      + (Trotter error, if TrotterDomain)
                      + (coherent term approximation, if with_coherent=true on Time/Trotter)
```

The DM evolution (`run_thermalization`) has all the same errors except sampling noise. So `||rho_traj - rho_DM||` isolates only sampling noise + any trajectory-specific bugs, which is what you want. But `||rho_traj - gibbs||` includes the entire error budget, and a test using this metric can fail because of quadrature error, not a trajectory bug. Conversely, it can pass despite a trajectory bug if errors cancel.

**Why it happens:**
The natural validation instinct is "does the trajectory produce the Gibbs state?" but the trajectory's job is to reproduce the CPTP channel defined by the Lindbladian, not to produce the exact Gibbs state (which itself is only approximately the fixed point for non-Bohr domains). Confusing "matches Gibbs" with "matches DM" leads to either false failures or false confidence.

**How to avoid:**
- **Primary test:** Compare `rho_traj` vs `rho_DM` (both evolved with identical parameters). The only difference should be sampling noise. Use tolerance `C / sqrt(ntraj)` where `C` is estimated from the variance of a few trial runs.
- **Secondary test:** Compare both `rho_traj` and `rho_DM` to the Gibbs state separately. Document the expected error from each source so deviations can be attributed.
- **Error budget test:** For a 2-qubit system, compute the exact Lindbladian fixed point (via `run_lindbladian`) and verify `||rho_DM_evolved - fixed_point||` matches the expected delta-step error `O(delta)`.
- Never test trajectory correctness by comparing to the Gibbs state alone.

**Warning signs:**
- Test tolerance must be set very loose (e.g., 0.1) to pass -- indicates you are fighting approximation error, not isolating trajectory correctness.
- Test passes for BohrDomain (exact) but fails for TimeDomain (approximate) -- this is expected behavior, not a bug, if the tolerance is set against Gibbs.
- Increasing `ntraj` does not improve agreement with Gibbs beyond a floor -- the floor is the deterministic approximation error.

**Phase to address:**
Phase 2 (Cross-Validation) -- the central phase for trajectory validation.

---

### Pitfall 3: Wrong Delta-Step Error Scaling Expectation (O(delta^2) per Step vs O(delta) Total)

**What goes wrong:**
Chen 2025 Theorem III.1 establishes that each delta-step of the weak measurement scheme introduces `O(delta^2)` error in trace distance. Over `T/delta` steps, the total accumulated error is `O(T * delta)`, i.e., `O(delta)` for fixed total evolution time `T`. A developer testing the trajectory code might expect either: (a) per-step error to be `O(delta)` (too pessimistic), or (b) total error to be `O(delta^2)` (too optimistic). Testing with wrong scaling expectations leads to either dismissing real bugs as "expected error" or flagging correct behavior as bugs.

**Why it happens:**
The `O(delta^2)` per-step comes from the residual Kraus term `S = (2*alpha - delta)*R - alpha^2 * R^2` which is `O(delta^2)` (since `alpha = 1 - sqrt(1-delta) ~ delta/2 + O(delta^2)`). But errors accumulate over `T/delta ~ 1/delta` steps, giving `O(delta)` total. This is standard for first-order integrators but can be confused with the Trotter literature where "first order" means `O(delta)` per step.

**How to avoid:**
- Write a delta-scaling test: run DM thermalization for a fixed total time `T` at multiple delta values (0.1, 0.05, 0.02, 0.01). Measure `||rho(T) - rho_exact(T)||` where `rho_exact` comes from matrix exponentiation of the Lindbladian. Verify the error scales as `O(delta)`.
- For trajectory-DM comparison, the delta-step error is the same in both, so it cancels. Use this as the primary validation.
- Document the expected scaling in test comments: "Per step: O(delta^2). Total over T/delta steps: O(T*delta)."

**Warning signs:**
- Error does not decrease when halving delta -- indicates a bug, not just large constants in the error bound.
- Error decreases faster than O(delta) -- suspicious, may indicate error cancellation hiding a bug at different parameters.
- Error matches O(delta^2) total -- too good, likely only tested one parameter point where higher-order terms dominate.

**Phase to address:**
Phase 3 (Error Hierarchy Testing) -- dedicated to verifying error scaling.

---

### Pitfall 4: Stochastic Test Flakiness from Fixed Seeds Hiding Real Bugs

**What goes wrong:**
To make trajectory tests deterministic, a fixed random seed is used. The test passes because the specific sequence of random numbers avoids the code path where the bug lives (e.g., the fallback branch in the cumulative scan is never triggered with that seed). When the seed changes (CI upgrade, Julia version change, thread count change), the test suddenly fails, appearing as a flaky test when it is actually exposing a real bug that was hidden.

**Why it happens:**
Julia's `Random.seed!(42)` produces a deterministic sequence, but the trajectory branching logic means that the "no-jump" vs "jump" vs "residual" branch taken depends sensitively on the random draw. With seed 42 and delta=0.1, perhaps 95% of steps are no-jump and the dissipative scanning code is barely exercised. A different seed might trigger the dissipative branch 20% of the time, exposing the fallback-masking issue (Pitfall 1 from prior research).

**How to avoid:**
- **Two-tier testing strategy:**
  1. **Deterministic property tests (no randomness):** Verify channel CPTP properties, probability normalization, single-step DM evolution. These use `step_along_trajectory!` with a known state and check branch probabilities without actually sampling.
  2. **Statistical convergence tests (many trajectories, no fixed seed):** Run N_traj=1000+ trajectories, compute mean and standard error, assert that `||rho_traj_mean - rho_DM|| < 3 * sigma / sqrt(N_traj)`. This test is designed to pass at the 99.7% level; if it fails, it is a real signal.
- Never use a single fixed-seed trajectory test as the primary correctness gate.
- If using seeds for debugging reproducibility, run the statistical test first to establish correctness, then use seeds for debugging specific failures.

**Warning signs:**
- Test passes with `seed=42` but fails with `seed=43` -- the code has a bug that is path-dependent.
- Test passes with `ntraj=1` and fixed seed but fails with `ntraj=100` and no seed -- the single trajectory was not representative.
- Test that previously passed starts failing after unrelated Julia/BLAS version update -- the seed-dependent path changed.

**Phase to address:**
Phase 2 (Test Infrastructure) -- establish the two-tier pattern before writing trajectory validation tests.

---

### Pitfall 5: DM Thermalization's Random Jump Selection vs Trajectory's Deterministic Full Map

**What goes wrong:**
`run_thermalization` (`furnace.jl` line 119) randomly selects one jump operator per step (`idx = rand(rng, 1:length(jumps))`) and applies a rescaled channel (`gamma_norm_factor / jump_prob`). The trajectory code in `step_along_trajectory!` applies the full Kraus map (summing over all jumps in each step) without per-jump randomization. These are two different CPTP channels that converge to the same Lindbladian in the continuous-time limit, but for finite delta they are not identical. Naively comparing `run_thermalization(delta=0.1, 200 steps)` with `run_trajectories(delta=0.1, 200 steps)` will show discrepancy that is NOT a bug but a genuine mathematical difference in the channels.

**Why it happens:**
The DM code implements the randomized channel from Chen 2025 Eq. 3.2 (randomly pick one jump per step, compensate by 1/p_jump scaling). The trajectory code implements the full deterministic channel (all jumps contribute to R and the cumulative scan in each step). Both converge to `exp(L*t)` as `delta -> 0`, but at finite delta they differ by `O(delta)` terms.

**How to avoid:**
- **For trajectory-DM comparison:** either (a) modify the DM code to also apply the full map (loop over all jumps, no randomization), or (b) compare both to the exact Lindbladian fixed point and verify each is within its expected error.
- **Preferred approach for testing:** Create a `run_thermalization_deterministic` variant that applies all jumps per step (no `rand`), matching the trajectory channel. This isolates trajectory-specific issues from channel-variant differences.
- The quantity `rho_DM - rho_traj` should have two contributions: (1) sampling noise from trajectory stochasticity, and (2) channel difference from randomized vs full map. Contribution (2) can dominate at large delta and mask or mimic trajectory bugs.

**Warning signs:**
- Discrepancy between DM and trajectory is `O(delta)` and shrinks linearly when delta is halved -- this is the expected channel difference, not a bug.
- Discrepancy does NOT shrink with delta -- this IS a bug (either in DM or trajectory).
- Discrepancy varies with the DM seed but not with ntraj -- this is the DM randomization noise.

**Phase to address:**
Phase 2 (Cross-Validation) -- must be resolved before interpreting any trajectory-DM comparison results.

---

### Pitfall 6: Probability Normalization Sum Drift Masked by Fallback Logic

**What goes wrong:**
In `step_along_trajectory!`, the total probability branches are:
- `p_nojump = ||K0 psi||^2`
- `p_res = ||U_res psi||^2`
- `p_jump_total = delta * <psi|R|psi>`

These should sum to 1.0 (the CPTP condition for the state `psi`). However, `p_jump_total` is computed from the precomputed matrix `R`, while the per-jump probabilities in the dissipative scan are computed on-the-fly using `delta * base_prefactor * transition(w) * ||A_{a,w} psi||^2`. If the `base_prefactor` in `step_along_trajectory!` has even a slightly different formula from the one in `precompute_R()` (e.g., a missing factor, a sigma vs beta confusion), then `sum_of_per_jump_probs != p_jump_total`, and the cumulative scan will either: (a) fail to reach `target`, triggering the fallback at lines 564-568, or (b) overshoot before scanning all jumps. The fallback picks the last computed candidate, silently producing the wrong jump distribution.

**Why it happens:**
The `base_prefactor` formula is duplicated in three places: `precompute_R(::EnergyDomain)` (line 132), `precompute_R(::Time/Trotter)` (line 189), and each `step_along_trajectory!` variant (lines 441, 604). These must be identical. During refactoring, any change to one that is not propagated to the others creates a drift.

**How to avoid:**
- **Immediate:** Add a debug assertion in `step_along_trajectory!` that checks `|csum - p_jump_total| / max(p_jump_total, 1e-15) < 0.01` after the dissipative scan completes. This catches normalization drift at runtime. Disable in production via `@assert` or `@debug`.
- **Structural:** Extract the `base_prefactor` computation into a single function `dissipative_base_prefactor(domain, config, precomputed_data)` called by both `precompute_R` and `step_along_trajectory!`.
- **Test:** For a 1-qubit system, manually compute all branch probabilities and verify they sum to 1.0 to machine precision.
- **Monitoring:** Count how often the fallback fires. If `fallback_count / total_steps > 0.001`, there is a normalization bug.

**Warning signs:**
- Fallback fires frequently (> 0.1% of steps). Currently invisible because there is no counter.
- `csum` at the end of the scan is consistently less than or greater than `target` -- indicates systematic bias in per-jump rates vs `R`.
- Changing the order of jump operators changes the trajectory-averaged result (it should not, since the scan covers all of them).

**Phase to address:**
Phase 1 (Trajectory Sampling Fix) -- add the normalization assertion when fixing two-stage sampling.

---

### Pitfall 7: Coherent Unitary Applied to Wrong State in Trajectory Step

**What goes wrong:**
In `step_along_trajectory!` (lines 472-477 for Time/Trotter, 633-638 for Energy), the coherent unitary `U_B` is applied to `psi` BEFORE the branch selection, but AFTER the branch probabilities (`p_nojump`, `p_res`, `p_jump_total`) have been computed from the pre-U_B state. This means:
1. The probabilities are computed for state `psi`.
2. `psi` is overwritten to `U_B * psi`.
3. The branch is selected based on old probabilities.
4. If the no-jump branch is taken, `psi` is set to `ws.psi_tmp` which was computed from the pre-U_B `psi` (i.e., `K0 * old_psi`), NOT `K0 * U_B * psi` or `U_B * K0 * psi`.

This is actually the correct ordering for the Chen construction where the coherent unitary and the Kraus map are applied sequentially (the step is `U_B` followed by the Kraus map, or vice versa). BUT: the `psi` is mutated to `U_B * psi` on line 475, and then the no-jump branch uses `ws.psi_tmp` (which contains `K0 * original_psi`), effectively discarding the `U_B` application in the no-jump case. In the jump and residual branches, `ws.Rpsi` contains computations from the original `psi`, not the U_B-evolved one.

If the intent is "apply U_B then apply Kraus map to U_B*psi", the current code is wrong because K0 and the jump operators are applied to the pre-U_B state. If the intent is "apply Kraus map then apply U_B", the current code is also wrong because U_B is applied to `psi` but then overwritten by the branch result (which comes from pre-U_B `psi`).

**Why it happens:**
The DM code applies `U_B rho U_B'` first (lines 166-169 in `jump_workers.jl`), then the Kraus map. The trajectory code appears to attempt the same ordering but the mutation of `psi` in-place creates a state where U_B has been applied to the original state while all cached quantities (K0*psi, R*psi, U_res*psi) were computed from the original. The branch then selects from cached quantities, discarding U_B.

**How to avoid:**
- Clarify the mathematical ordering: In the DM code, the step is `E(rho) = U_B * (K0 rho K0' + sum_jumps + Ures rho Ures') * U_B'`, i.e., coherent rotation applied to the entire Kraus output.
- For trajectories, the correct implementation is: (1) apply Kraus step to get `psi_new`, (2) apply `psi_final = U_B * psi_new`, (3) normalize. Move the U_B application AFTER the branch selection, not before.
- Alternatively, if U_B is meant to be absorbed into each Kraus operator (i.e., `K0_eff = U_B * K0`, `K_jump_eff = U_B * K_jump`), precompute this.
- Write a test: single step on a known 2-qubit state, compare the trajectory step output to the DM step output.

**Warning signs:**
- Trajectory with `with_coherent=true` converges to a different state than DM with `with_coherent=true`.
- The coherent term has no observable effect on trajectory results (because U_B is applied then overwritten).
- This is currently masked by the `trotter` undefined variable bug (Pitfall 6 from prior research) which prevents `with_coherent=true` from running at all.

**Phase to address:**
Phase 1 (Trajectory Sampling Fix) -- must be fixed alongside the sampling fix, since both affect the step function.

---

### Pitfall 8: Testing Trace Distance Convergence with Insufficient Trajectories, Getting False Confidence

**What goes wrong:**
A test runs 100 trajectories and checks `trace_distance(rho_traj_avg, rho_DM) < 0.05`. This passes because with 100 trajectories the statistical error is `O(1/sqrt(100)) = 0.1` for unit-trace density matrices, and the tolerance is within noise. The developer concludes the trajectory code is correct. But the actual bias (from a sampling bug) is 0.03, well within the noise floor. Increasing to 10,000 trajectories would reveal the bias as `0.03 > 3 * 0.01 = 0.03` (borderline), while 100,000 trajectories would unambiguously detect it.

**Why it happens:**
The standard error of trajectory-averaged density matrix elements scales as `sigma / sqrt(N_traj)`, where `sigma` depends on the observable variance. For density matrix elements of an n-qubit system thermalized near the Gibbs state, `sigma ~ 1/dim` for diagonal elements but can be much larger for off-diagonal elements during early evolution. Developers often use `ntraj = 100` because 1000+ is slow, and the test "passes," but the statistical power is too low to detect `O(0.01)` biases.

**How to avoid:**
- **Power analysis before writing the test:** For a system of dimension `dim`, estimate the per-trajectory variance of `tr(O * |psi><psi|)` for a representative observable `O`. The standard error of the mean is `sigma / sqrt(N_traj)`. To detect a bias of size `epsilon`, need `N_traj > (3 * sigma / epsilon)^2`.
- **For 2-qubit systems (dim=4):** `sigma ~ 0.25` for diagonal elements, so `N_traj > (3 * 0.25 / 0.01)^2 = 5625` to detect 1% bias at 3-sigma.
- **For 4-qubit systems (dim=16):** `sigma ~ 0.06`, so `N_traj > (3 * 0.06 / 0.01)^2 = 324` -- more forgiving.
- Use the standard error from the trajectory ensemble itself: `se = std(observable_per_traj) / sqrt(N_traj)`. Assert that `|mean - expected| < 3 * se`.
- **Quick heuristic for CI:** Use 2-qubit tests with N_traj=500 for fast feedback, and 3-qubit tests with N_traj=2000 as a slower nightly gate.

**Warning signs:**
- Test passes at N_traj=100 but fails at N_traj=10000 -- the bias was always there.
- Computed standard error is large compared to tolerance -- the test has no power to detect biases.
- Test tolerance is suspiciously large (> 0.05 for trace distance) -- it is not testing anything useful.

**Phase to address:**
Phase 2 (Cross-Validation) -- determine trajectory counts and tolerances before writing validation tests.

---

### Pitfall 9: Cholesky Residual S Goes Non-PSD at Large Coupling, Silently Breaking CPTP

**What goes wrong:**
The residual matrix `S = (2*alpha - delta)*R - alpha^2 * R^2` computed in `build_trajectoryframework` (line 73 of `trajectories.jl`) and `jump_contribution!` (line 466-476 of `jump_workers.jl`) uses `cholesky!(Hermitian(S), check=false)`. When `check=false`, a failed Cholesky factorization does not throw an error; it silently produces an upper-triangular factor with NaN or incorrect values. The resulting `U_residual` breaks the CPTP property: the state can develop negative norm components, density matrices lose positive-semidefiniteness, and trace preservation fails.

**Why it happens:**
`S` is `O(delta^2)` with eigenvalues proportional to `delta^2 * eigenvalues_of_R * (1 - alpha * eigenvalues_of_R)`. When `R` has eigenvalues close to `1/alpha ~ 2/delta`, the factor `(1 - alpha * lambda_R)` goes negative. This happens when the total dissipation rate `||R||` is large relative to `1/delta`, i.e., when delta is too large for the given system. The `10 * eps(Float64)` diagonal shift is `~2.2e-15`, far too small to fix eigenvalues that are `O(delta^2) ~ 0.01` in magnitude.

**How to avoid:**
- After computing `S`, check `issuccess(cholesky_S)`. If it fails, either: (a) reduce delta, (b) use eigendecomposition `S = V * diag(max.(lambdas, 0)) * V'` and set `U_res = diag(sqrt.(max.(lambdas, 0))) * V'`, or (c) increase the diagonal shift adaptively.
- Add a pre-flight check in `build_trajectoryframework`: compute `eigmin(Hermitian(S))` and warn if it is negative.
- In tests, assert `all(eigvals(Hermitian(S)) .> -1e-12)` for every configuration tested.
- Provide guidance on delta selection: `delta < 1 / (2 * opnorm(R))` ensures `alpha * lambda_max(R) < 1` and thus `S` is PSD.

**Warning signs:**
- NaN appearing in trajectory state vectors.
- `norm(psi)` growing or shrinking despite explicit normalization (because `U_res * psi` has NaN components).
- Trace of trajectory-averaged rho is not close to 1.0.
- `issuccess(cholesky_S)` returns `false` (not currently checked).

**Phase to address:**
Phase 1 (Trajectory Sampling Fix) -- fix before any validation relies on trajectory results.

---

## Technical Debt Patterns

| Shortcut | Immediate Benefit | Long-term Cost | When Acceptable |
|----------|-------------------|----------------|-----------------|
| Duplicated `base_prefactor` formula across `precompute_R` and `step_along_trajectory!` | Quick to implement, avoids function call overhead | Normalization drift when one is updated but not the other; probability sum assertion catches this | Never -- extract to shared function |
| `check=false` on Cholesky of S matrix | Avoids crash when S is marginally non-PSD | Silent garbage in U_residual; breaks CPTP guarantee without warning | Only with a runtime PSD guard (check eigenvalues first) |
| Coherent unitary placement before branch selection in `step_along_trajectory!` | "Works" for `with_coherent=false` (the only tested path) | Wrong Kraus map when coherent term is enabled; blocks a major feature | Never -- fix ordering before enabling coherent trajectories |
| ~200 lines of commented-out old trajectory code | Preserves reference to old `KrausFramework` implementation | Obscures the actual API; makes code review harder; old patterns may be inadvertently copied | Delete after new trajectory code is validated (keep in git history) |
| Using `round(Int, total_time / delta)` vs `ceil(Int, ...)` for step count | Works when total_time is exact multiple of delta | Silently drops the last fractional step; trajectory evolves for slightly less time than DM | Never -- use `ceil` consistently (already done in new code, but old code uses `round`) |

## Integration Gotchas: Trajectory vs Density Matrix Cross-Validation

| Integration Point | Common Mistake | Correct Approach |
|-------------------|----------------|------------------|
| Channel convention | Assuming DM `run_thermalization` and `run_trajectories` implement the same CPTP map | They implement different maps (randomized vs deterministic). Compare both to `expm(L*delta)` or create a `run_thermalization_deterministic` that applies all jumps per step |
| Initial state | Using maximally mixed state `I/dim` for comparison (it is close to Gibbs at high temperature and masks convergence issues) | Use a pure state far from Gibbs (e.g., a computational basis state or the ground state). This maximizes the signal in convergence diagnostics |
| Time grid alignment | Comparing trajectory snapshot at step `k * save_every * delta` to DM state at same time, but DM uses `config.delta` which may differ from trajectory `delta` | Ensure identical delta and total evolution time between DM and trajectory runs |
| Basis convention | Trajectory operates in eigenbasis (for Energy/Time domains) but DM `evolving_dm` may be in computational basis depending on setup | Verify that `rho_traj_avg` and `rho_DM` are in the same basis before computing trace distance |
| Normalization of observables | Comparing `<O>_traj` vs `tr(O * rho_DM)` when `O` is not Hermitian | Always use Hermitian observables for trajectory expectation values. Non-Hermitian observables give complex expectation values whose imaginary part should vanish in the average but adds noise |

## Performance Traps in Trajectory Validation

| Trap | Symptoms | Prevention | When It Breaks |
|------|----------|------------|----------------|
| Running N_traj=10000 trajectories for a 4-qubit system in a single-threaded CI test | Test takes > 10 minutes; CI times out or developers skip it | Use 2-qubit system (dim=4) for large-N_traj statistical tests; use 4-qubit system with N_traj=200 for smoke tests only | At N_traj > 1000 for dim > 64 (6+ qubits) |
| Computing `eigvals(Hermitian(S))` in every step for PSD validation | 100x slowdown per trajectory step (eigvals is O(dim^3)) | Validate PSD of S only once during `build_trajectoryframework`, not per step. Per-step check is only needed during debugging | At dim > 64 with > 100 steps |
| Full energy grid scan in dissipative branch even when most `||A_w psi||^2` are negligible | Per-step cost is O(num_energy_labels * dim^2); with 2^12 labels and dim=256, each step costs ~10^8 FLOPS | For validation tests, use small `num_energy_bits` (8-10) so the grid is manageable. Production optimization (pruning, early exit) is a separate concern | At num_energy_bits >= 12 for dim >= 128 |
| Hermitianizing `rho_traj_avg` at every accumulation step | O(dim^2) copy per trajectory per accumulation | Only Hermitianize the final averaged rho, not during accumulation. The code already does this correctly (line 359: `rho_mean .= 0.5 .* (rho_mean .+ rho_mean')` only after the loop) | Not currently a problem; flagging to prevent regression |

## Numerical Validation Pitfalls

### Pitfall: Distinguishing Quadrature Error from Trajectory Bugs

**What goes wrong:**
The Energy, Time, and Trotter domains approximate the continuous Lindbladian with discrete sums. This introduces an `O(w0)` quadrature error in the Lindbladian itself. When testing trajectory convergence, this error is present in both the DM and trajectory paths, so it cancels in the DM-vs-trajectory comparison. BUT: if a test compares to the exact Gibbs state or to the BohrDomain Lindbladian fixed point, the quadrature error appears as a discrepancy. A developer seeing `||rho_Energy - gibbs|| = 0.02` may hunt for a trajectory bug when the cause is simply coarse energy discretization.

**Prevention:**
- Establish the "error floor" for each domain by running the DM code alone and measuring `||rho_DM(T) - gibbs||`. This is the best-achievable accuracy for that domain/config.
- Only expect trajectory accuracy to match this floor (plus sampling noise), never better.
- Use BohrDomain (no discretization error) for pure trajectory-correctness tests. Reserve Energy/Time/Trotter domains for integration tests that verify "trajectory matches DM within this domain's error budget."

### Pitfall: Hermiticity Enforcement Masking Non-Hermitian Errors

**What goes wrong:**
Both DM and trajectory code enforce Hermiticity at every step: `rho .= 0.5 .* (rho .+ rho')` (DM: line 270, 385, 498 in jump_workers.jl; trajectory: line 359, 407 in trajectories.jl). This makes the output look valid even if the underlying channel has a non-Hermitian error (e.g., from incorrect conjugate-transpose ordering in the jump operator). The Hermitianization absorbs the anti-Hermitian part, hiding bugs that produce `rho` with large `||rho - rho'||`.

**Prevention:**
- In tests, check `||rho_raw - rho_raw'|| / ||rho_raw||` BEFORE Hermitianization. This should be `O(epsilon_machine)` for a correct channel. If it is `O(delta)` or larger, there is a bug in the channel construction.
- For trajectories, the analogous check is that `|Im(<psi|O|psi>)| / |Re(<psi|O|psi>)|` is small for Hermitian observables O.
- Add a `@debug` check inside the step function that flags when the anti-Hermitian component exceeds `100 * eps(Float64) * norm(rho)`.

### Pitfall: Positive-Semidefiniteness Lost Through Accumulated Floating-Point Error

**What goes wrong:**
After many steps, the trajectory-averaged density matrix `rho_mean` can develop small negative eigenvalues (e.g., `-1e-14`) due to accumulated floating-point arithmetic in `_accumulate_density_matrix!` (BLAS `gerc!` at line 239). The `is_density_matrix` check in `qi_tools.jl` uses `round.(eigvals(rho), digits=13)` or `digits=15` (two different overloads!) and rejects matrices with negative eigenvalues after rounding. This means the validation function can either pass a matrix with eigenvalue `-1e-16` (rounded to zero) or fail a matrix with eigenvalue `-3e-14` (not rounded to zero), depending on which overload is called.

**Prevention:**
- Use a single tolerance-based PSD check: `minimum(eigvals(Hermitian(rho))) > -tol` where `tol = 100 * dim * eps(Float64)`.
- Never use `round` for numerical validation -- it creates sharp discontinuities where `1e-13.4` rounds to zero but `1e-13.6` does not.
- For trajectory tests, Hermitianize and then use `Hermitian(rho)` wrapper which guarantees real eigenvalues.

---

## "Looks Done But Isn't" Checklist

- [ ] **Two-stage jump sampling:** The linear scan in `step_along_trajectory!` looks like it samples jumps, but it uses flat `(a, omega)` scanning, not two-stage `a` then `omega|a`. Verify by checking that the sampling distribution is correct for non-uniform jump operators.
- [ ] **Coherent term in trajectories:** `build_trajectoryframework` calls `precompute_coherent_total_B` with undefined `trotter` kwarg (line 53). Currently crashes for `with_coherent=true`. Even after fixing the kwarg, the U_B placement in the step function (before branch selection) is mathematically wrong.
- [ ] **Probability normalization:** `p_nojump + p_res + p_jump_total` looks like it sums to 1.0, but the `p_jump_total = delta * expR` is computed from `R` while the individual scan uses `base_prefactor * transition(w) * n2`. Verify the two are consistent by checking `|csum - p_jump_total| < epsilon`.
- [ ] **Error scaling verification:** "Delta-step error is O(delta^2)" -- this is per-step, NOT total. The total error over `T/delta` steps is `O(T * delta)`. Tests must use the correct scaling or they validate nothing.
- [ ] **Test suite runnable via `Pkg.test()`:** Current test files use `includet` (Revise) and `include("../src/QuantumFurnace.jl")`. No `test/runtests.jl` exists. Tests are interactive scripts, not CI-compatible test suites.
- [ ] **S matrix PSD check:** `cholesky!(... check=false)` looks like it handles the residual, but `check=false` means silent failure. Verify `issuccess(cholesky_S)` or check eigenvalues.
- [ ] **Basis consistency:** `run_trajectories` returns `rho_mean` but the basis (eigenbasis vs computational) is implicit. Verify it matches the basis of the DM comparison target.

## Recovery Strategies

| Pitfall | Recovery Cost | Recovery Steps |
|---------|---------------|----------------|
| Flat sampling gives wrong steady state | MEDIUM | Implement two-stage sampling; rewrite the inner loop of `step_along_trajectory!` dissipative branch; add non-uniform jump test |
| Error budget confusion (trajectory vs Gibbs) | LOW | Document the error hierarchy; add the DM-only error floor as a constant in test comments; restructure tests to compare traj-vs-DM first |
| Wrong delta scaling expectation | LOW | Add delta-sweep test; document "O(delta^2) per step, O(delta) total" in code and tests |
| Flaky seed-dependent tests | MEDIUM | Restructure to two-tier pattern; convert existing trajectory_test.jl from interactive script to proper Test module with statistical assertions |
| DM vs trajectory channel mismatch | MEDIUM | Create `run_thermalization_deterministic`; or add a test that compares both to `expm(L*delta)` for a single step |
| Probability normalization drift | LOW | Extract `base_prefactor` to shared function; add debug assertion for probability sum |
| Coherent unitary ordering bug | MEDIUM | Move U_B application after branch selection; verify against DM single-step result with `with_coherent=true` |
| Insufficient trajectory count | LOW | Compute required N_traj from power analysis; use 2-qubit system for heavy statistical tests |
| Cholesky non-PSD silent failure | MEDIUM | Replace `check=false` with explicit PSD guard; add eigenvalue check in framework construction |
| PSD lost through accumulation | LOW | Use single tolerance-based check; fix dual `is_density_matrix` overloads |

## Pitfall-to-Phase Mapping

| Pitfall | Prevention Phase | Verification |
|---------|------------------|--------------|
| Flat sampling (wrong steady state) | Phase 1: Sampling Fix | Non-uniform jump test converges to correct Gibbs state |
| Error budget confusion | Phase 2: Cross-Validation | DM-vs-trajectory distance < C/sqrt(N_traj) for C estimated from variance |
| Wrong delta scaling | Phase 3: Error Hierarchy | Delta-sweep test shows O(delta) total error scaling |
| Flaky seed-dependent tests | Phase 2: Test Infrastructure | All tests pass without fixed seeds; statistical tests use 3-sigma bounds |
| DM vs trajectory channel mismatch | Phase 2: Cross-Validation | Single-step channel comparison shows agreement to O(eps_machine) when using same map |
| Probability normalization drift | Phase 1: Sampling Fix | `|csum - p_jump_total| / p_jump_total < 1e-10` assertion passes |
| Coherent unitary ordering | Phase 1: Sampling Fix | Coherent trajectory matches coherent DM for 2-qubit system |
| Insufficient trajectories | Phase 2: Cross-Validation | Power analysis determines N_traj; tests document required count |
| Cholesky non-PSD | Phase 1: Sampling Fix | `eigmin(S) > -1e-12` for all test configs; `issuccess(cholesky_S)` checked |
| Hermiticity masking | Phase 3: Error Hierarchy | Pre-Hermitianization anti-Hermitian norm < 100*eps for DM steps |
| PSD accumulation loss | Phase 2: Test Infrastructure | Single PSD tolerance function; both `is_density_matrix` overloads unified |
| Quadrature vs bug confusion | Phase 3: Error Hierarchy | DM-only error floor documented for each domain; trajectory tests use DM floor as baseline |

## Sources

- Codebase analysis: `/Users/bence/code/QuantumFurnace.jl/src/trajectories.jl`, `jump_workers.jl`, `coherent.jl`, `furnace.jl`, `qi_tools.jl`, `errors.jl`
- [Monte Carlo wave function method: robust adaptive algorithm and convergence study (Abdelhafez et al. 2019)](https://arxiv.org/abs/1803.08589)
- [Kraus is King: High-order CPTP low rank method for Lindblad master equation (Cai & Lu 2025)](https://arxiv.org/abs/2409.08898)
- [QuTiP Monte Carlo Solver documentation](https://qutip.readthedocs.io/en/latest/guide/dynamics/dynamics-monte.html)
- [QuantumOptics.jl MCWF documentation](https://docs.qojulia.org/timeevolution/mcwf/)
- [An efficient and exact noncommutative quantum Gibbs sampler (Chen, Kastoryano, Gilyen 2023)](https://arxiv.org/abs/2311.09207)
- [Efficient quantum Gibbs samplers with KMS detailed balance condition (Ding, Li, Lin 2024)](https://arxiv.org/abs/2404.05998)
- [To Seed or Not to Seed: Empirical Analysis of Seeds for Testing in ML Projects](http://misailo.cs.illinois.edu/papers/icst22.pdf)
- [FLEX: Fixing flaky tests by updating assertion bounds (Dutta et al. 2021)](https://dl.acm.org/doi/10.1145/3468264.3468615)
- [Quantum localization bounds Trotter errors in digital quantum simulation (Heyl et al. 2019)](https://www.science.org/doi/10.1126/sciadv.aau8342)
- [Quantum trajectory framework for general time-local master equations (Breur et al. 2022)](https://www.nature.com/articles/s41467-022-31533-8)

---
*Pitfalls research for: QuantumFurnace.jl -- Trajectory simulation validation and testing*
*Researched: 2026-02-13*
*Milestone: v1.0 Trajectories*
