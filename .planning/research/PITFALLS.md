# Domain Pitfalls: Spectral Gap Estimation from Trajectory Observables

**Domain:** Adding spectral gap estimation via exponential fitting to observable decay curves from quantum trajectory simulations, cross-validated against exact Liouvillian eigenvalues
**Researched:** 2026-02-16
**Confidence:** HIGH (codebase analysis + established numerical analysis + quantum open systems literature)

---

## Critical Pitfalls

Mistakes that cause wrong spectral gap estimates, silent numerical failures, or fundamentally flawed cross-validation.

---

### Pitfall 1: Confusing Complex Liouvillian Eigenvalue with Real Observable Decay Rate

**What goes wrong:**
The existing `run_lindbladian()` in `furnace.jl` (line 22) computes `spectral_gap = eigvals_near_zero[gap_index]`, which is a `Complex{T}` value -- the second eigenvalue of the Liouvillian superoperator. The Liouvillian is non-Hermitian, so its eigenvalues are generically complex: `lambda_2 = -gamma + i*omega` where `gamma > 0` is the decay rate and `omega` is an oscillation frequency. The observable decay rate from trajectory fitting will be a real number `gamma_fit`. The cross-validation must compare `gamma_fit` against `-Re(lambda_2)`, NOT against `abs(lambda_2)` or `real(lambda_2)`.

**Why it happens:**
The `LindbladianResult` struct stores `spectral_gap::Complex{T}` (structs.jl line 263). The existing simulation script `main_liouv.jl` (line 130) prints `abs(real(liouv_result.spectral_gap))`, which happens to be correct for the comparison but hides the sign convention. When writing the cross-validation, it is natural to compare `abs(spectral_gap)` (the modulus) against the fitted decay rate, which would be WRONG whenever `Im(lambda_2) != 0` because `|lambda_2| = sqrt(gamma^2 + omega^2) > gamma`.

**The physics:**
For a Lindbladian L, the semigroup `exp(Lt)` has eigenvalues `exp(lambda_k * t)`. For an observable `O`, the deviation from steady state decays as:

```
<O>(t) - <O>_ss = sum_k c_k * exp(lambda_k * t) = sum_k c_k * exp(-gamma_k * t) * exp(i * omega_k * t)
```

The REAL part of the eigenvalue controls the exponential envelope. If you measure `<O>(t)` from trajectories, the oscillations `exp(i*omega*t)` average out for a real-valued observable, and the dominant visible decay rate is `gamma_1 = -Re(lambda_2)`.

However, for KMS-detailed-balanced Lindbladians (which QuantumFurnace implements), the Liouvillian is self-adjoint with respect to the GNS inner product, meaning all eigenvalues are real and non-positive. So for the exact KMS construction (BohrDomain), `Im(lambda_2) = 0`. But for approximate domains (Energy, Time, Trotter), the detailed balance is approximate, and eigenvalues may have small imaginary parts. The cross-validation must handle both cases.

**Consequences:**
- Comparing `abs(spectral_gap)` to the fitted rate gives a systematic overestimate of agreement (the modulus is always >= the real part).
- For TrotterDomain at coarse Trotter steps, `Im(lambda_2)` can be non-negligible, making the discrepancy significant.
- False confidence in the cross-validation: you think trajectory-based and exact methods agree when they do not.

**Prevention:**
1. Define the comparison quantity explicitly: `exact_gap = -real(liouv_result.spectral_gap)`. Assert it is positive.
2. Also record `imag(liouv_result.spectral_gap)` and flag when `|Im/Re| > threshold` (e.g., 0.1), because a large imaginary part means observable oscillations that complicate exponential fitting.
3. For BohrDomain (exact KMS), assert `abs(imag(spectral_gap)) < 1e-10` as a sanity check.
4. When the imaginary part is significant, the observable time series will show damped oscillations, not pure exponential decay. The fitting model must account for this (see Pitfall 3).

**Detection:**
- Print both `real(spectral_gap)` and `imag(spectral_gap)` in every cross-validation.
- Test: for BohrDomain, verify `imag(spectral_gap) ~= 0`.
- Test: for TrotterDomain with fine Trotter steps, verify `|Im/Re| < 0.01`.

**Phase to address:** Cross-validation phase -- the very first thing to get right when defining the comparison metric.

---

### Pitfall 2: Observable Basis Mismatch -- Eigenbasis vs Computational vs Trotter Basis

**What goes wrong:**
The total magnetization observable `M_z = sum_i Z_i` is naturally defined in the computational basis. The trajectory simulation evolves `psi` in the Hamiltonian eigenbasis (for Energy/Time domains) or Trotter eigenbasis (for TrotterDomain). If the observable matrix is constructed in the computational basis but applied to `psi` in the eigenbasis, `<psi|M_z_comp|psi>` gives a wrong expectation value. The decay curve of this wrong observable decays at the wrong rate, and the fitted spectral gap is meaningless.

**Why it happens:**
This is the EXACT same class of bug that caused the v1.2 quick-task 20 crisis: the GNS TrotterDomain test compared a Trotter-basis fixed point against an energy-eigenbasis Gibbs state, producing a spurious gap of 0.83 instead of 0.0807. The existing `build_convergence_observables()` in `convergence.jl` (lines 25-46) correctly transforms ZZ operators into the eigenbasis via `V' * O_comp * V`. The `build_convergence_observables_trotter()` (lines 56-77) correctly uses `V_T' * O_comp * V_T`. But when building a NEW observable (total magnetization), there is high risk of forgetting this transform.

**The insidious part:**
If the observable is diagonal in the eigenbasis (like the Hamiltonian `H`), the computational-basis and eigenbasis representations differ only by the unitary transform, and the EXPECTATION VALUE is basis-independent for the correct `psi`. But the MATRIX REPRESENTATION used in `_accumulate_measurements!` (trajectories.jl line 325: `mul!(tmp, observables[i], psi)`) must be in the SAME basis as `psi`. Since `psi` is in the eigenbasis, `observables[i]` must also be in the eigenbasis.

For total magnetization `M_z = sum_i Z_i`:
- Computational basis: `M_z_comp = sum_i pad_term([Z], num_qubits, i)`
- Eigenbasis: `M_z_eigen = V' * M_z_comp * V` where `V = hamiltonian.eigvecs`
- Trotter basis: `M_z_trotter = V_T' * M_z_comp * V_T` where `V_T = trotter.eigvecs`

**Consequences:**
- The observable time series looks physically plausible (starts near zero, evolves to some value) but the values are wrong.
- The fitted decay rate does not correspond to any physical decay mode.
- Cross-validation against exact spectral gap fails, and the failure is attributed to "statistical noise" rather than a basis error.
- This bug is particularly hard to catch because `<M_z>` for a random initial state is close to zero in any basis for the Heisenberg chain (symmetry), so the wrong-basis result might pass superficial sanity checks.

**Prevention:**
1. Follow the EXISTING pattern from `build_convergence_observables()`: always transform `O_eigen = V' * O_comp * V`.
2. Build the total magnetization observable in the SAME function that builds ZZ observables, using the same basis transform.
3. Add a regression test: for the eigenbasis energy observable `H_eigen`, verify that `tr(gibbs * H_eigen)` matches the analytical Gibbs energy. Then do the same for `M_z_eigen`: verify `tr(gibbs * M_z_eigen) = sum_i <Z_i>_gibbs`. This catches basis mismatches.
4. For TrotterDomain: use `build_convergence_observables_trotter()` (which uses `V_T`) rather than `build_convergence_observables()` (which uses `V`). Do NOT mix them.

**Detection:**
- Compare `<M_z>` from trajectory vs `tr(gibbs * M_z)` analytically. If they disagree at convergence, suspect basis mismatch.
- For the Heisenberg chain WITHOUT disorder, `<M_z>_gibbs = 0` by SU(2) symmetry. If `<M_z>` from trajectories converges to a nonzero value for the clean chain, something is wrong.

**Phase to address:** Observable construction phase -- before ANY trajectory runs with the new observable.

---

### Pitfall 3: Fitting Single Exponential to Multi-Exponential Decay

**What goes wrong:**
The observable decay `<O>(t) - <O>_ss` is a superposition of ALL Liouvillian eigenmodes, not just the slowest:

```
<O>(t) - <O>_ss = c_1 * exp(-gamma_1 * t) + c_2 * exp(-gamma_2 * t) + ... + noise
```

where `gamma_1 < gamma_2 < ...` are the decay rates (negative real parts of Liouvillian eigenvalues). The spectral gap is `gamma_1`. Fitting a single exponential `A * exp(-gamma * t)` to this data gives a `gamma_fit` that is BIASED HIGH -- it is a weighted average of all decay rates, dominated by the fast-decaying components at short times and the slow component at long times.

**Why it happens:**
For a general initial state (like the maximally mixed state) and a general observable (like total magnetization), the coefficients `c_k` can be non-negligible for many eigenmodes. The single-exponential fit minimizes the sum of squared residuals over the entire time range, which means:
- At short times: the fast modes dominate, pulling `gamma_fit` up.
- At long times: only the slowest mode survives, but the signal is buried in noise.
- The least-squares optimum is a compromise that overestimates the spectral gap.

For the Heisenberg chain at n=4 (dim=16, Liouvillian dim=256), there are 255 non-zero eigenvalues. The spectral gap `gamma_1` might be ~0.08 while the second eigenvalue `gamma_2` might be ~0.15 -- a factor of 2 larger. The observable `M_z` will have overlap with both modes.

**Consequences:**
- The fitted spectral gap is systematically too large (optimistic).
- Cross-validation against exact eigenvalue `gamma_1` fails with the fit consistently above the true value.
- The overestimate gets WORSE with fewer trajectories (more noise at long times, fit dominated by short-time behavior).

**Prevention:**
1. **Use late-time fitting only.** Discard the initial transient (t < t_cutoff) and fit only the long-time tail where the slowest mode dominates. The cutoff should be `t_cutoff ~ 3 / gamma_2` (roughly 3 decay times of the second mode). For n=4, this is estimable from the exact Liouvillian. For larger systems, start with a conservative cutoff (e.g., discard the first 50% of the time series).
2. **Fit log-linear.** Plot `log|<O>(t) - <O>_ss|` vs `t`. At late times this should be linear with slope `-gamma_1`. A linear fit to the log-data in the tail region is more robust than nonlinear exponential fitting.
3. **Use the exact steady-state value.** Subtract `<O>_ss = tr(gibbs * O)` (known analytically from the Gibbs state) rather than fitting it as a free parameter. This eliminates one degree of freedom and prevents the fitter from absorbing baseline errors into the decay rate.
4. **Multi-exponential awareness.** If cross-validation fails, try a two-exponential fit `c_1 * exp(-gamma_1 * t) + c_2 * exp(-gamma_2 * t)` and check if `gamma_1` from the two-exponential fit matches the exact spectral gap. But beware: two-exponential fitting is ill-conditioned (see Pitfall 5).
5. **Observable choice matters.** Choose observables with maximal overlap with the slowest eigenmode. The energy observable `<H>` may or may not have good overlap. The optimal observable is the gap mode itself (`liouv_result.gap_mode`), but this is only known from the exact diagonalization, creating a circularity. In practice, try multiple observables and see which gives the most consistent gap estimate.

**Detection:**
- Fitted gap is consistently 20-100% larger than the exact gap.
- Residuals show systematic curvature (not random scatter) at early times.
- Changing the fit window (start time) changes the fitted gap significantly.

**Phase to address:** Exponential fitting phase -- core fitting methodology.

---

### Pitfall 4: Trajectory Noise Destroying the Late-Time Signal

**What goes wrong:**
At late times `t >> 1/gamma_1`, the observable deviation `<O>(t) - <O>_ss` is exponentially small. With `N_traj` trajectories, the statistical error on `<O>(t)` scales as `sigma(t) / sqrt(N_traj)`, where `sigma(t)` is the single-trajectory variance. At late times, the signal `|c_1 * exp(-gamma_1 * t)|` drops below the noise floor `sigma / sqrt(N_traj)`. Beyond this point, the data is pure noise, and including it in the fit corrupts the estimate.

For the QuantumFurnace setup:
- Each trajectory is a single pure-state evolution. The variance of `<psi|O|psi>` for a single trajectory is `<O^2> - <O>^2`, which is O(1) for bounded observables like `M_z` (eigenvalues in [-n, +n]).
- The signal decays as `~exp(-gamma_1 * t)` where `gamma_1 ~ 0.08` for n=4 Heisenberg.
- After time `t = 50`, the signal is `exp(-4) ~ 0.018`. With 1000 trajectories, noise is `1/sqrt(1000) ~ 0.032`. The signal is already buried.

**Why it happens:**
This is the fundamental signal-to-noise problem of extracting slow decay rates from stochastic simulations. It is well-known in quantum Monte Carlo for imaginary-time correlation functions (see Sandvik and Vidal, "Excitation Gap from Optimized Correlation Functions"). The problem is intrinsic to the method, not a bug.

**Consequences:**
- If fitting includes the noisy tail, the fitter either: (a) finds a spurious minimum far from the true gap, or (b) converges to a shallow minimum with enormous uncertainty.
- The fitted `gamma_fit` can be arbitrarily wrong if the noise-dominated region is large relative to the signal region.
- Bootstrap confidence intervals on the fit become huge, but only if computed correctly (see Pitfall 8).

**Prevention:**
1. **Determine the noise floor empirically.** Compute `std(<O>(t))` across trajectories at each time point. The noise floor is `std / sqrt(N_traj)`. Exclude time points where `|<O>(t) - <O>_ss| < k * noise_floor` (e.g., k=2 or k=3).
2. **Use weighted least squares.** Weight each data point by `1 / variance(t)`. This automatically downweights the noisy tail. LsqFit.jl supports the `wt` parameter for weighted fitting.
3. **Increase trajectory count.** The noise floor drops as `1/sqrt(N_traj)`. To extend the usable time range by a factor of 2, you need 4x more trajectories (the signal drops by `exp(-gamma_1 * Delta_t) ~ exp(-0.08 * 25) ~ 0.14` while noise drops by `1/2`).
4. **Use coarser time binning at late times.** Average the observable over wider time windows at late times to reduce noise. This is equivalent to a running average, which does NOT bias the decay rate.
5. **For cross-validation at n=4,6 (where exact gap is known):** compute the expected signal-to-noise ratio at the planned trajectory count, and verify it is sufficient to resolve the gap before running the full experiment.

**Detection:**
- Plot `|<O>(t) - <O>_ss|` on a log scale alongside the noise floor `std/sqrt(N_traj)`. The signal should be well above the noise for at least ~3 decay times (`3/gamma_1`).
- If the error bars on the fit exceed 50% of the fitted value, the signal-to-noise is inadequate.
- The fit residuals in the late-time region should be consistent with random noise, not systematic.

**Phase to address:** Trajectory runner and fitting phases -- must plan trajectory count and time range together.

---

## Moderate Pitfalls

Mistakes that cause significant debugging time, wrong intermediate results, or wasted computation, but are recoverable.

---

### Pitfall 5: Exponential Fitting Sensitivity to Initial Parameter Guesses

**What goes wrong:**
Nonlinear least squares fitting of `A * exp(-gamma * t) + offset` is highly sensitive to the initial guess for `gamma`. The Levenberg-Marquardt algorithm (used by LsqFit.jl) converges to local minima. If the initial `gamma_0` is far from the true value, the fit may:
- Converge to a wrong local minimum (e.g., fitting to a fast mode instead of the slow mode).
- Diverge entirely (LsqFit returns NaN or throws an exception).
- Converge to a physically nonsensical result (`gamma < 0` or `gamma ~ 0`).

**Why it happens:**
The exponential model is highly nonlinear in the rate parameter. Small changes in `gamma` produce large changes in the residual at late times but small changes at early times. The least-squares landscape has narrow valleys and shallow basins, especially when fitting to noisy data.

For the QuantumFurnace use case, the spectral gap `gamma_1` ranges from ~0.01 to ~0.2 depending on the Hamiltonian, domain, and temperature. An initial guess of `gamma_0 = 1.0` (off by 10x) is enough to cause convergence failure.

**Prevention:**
1. **Use the log-linear estimate as the initial guess.** Before nonlinear fitting, compute `gamma_init = -slope(linear_fit(log|signal|, t))` over the last 30% of the signal-above-noise time range. This is cheap, robust, and typically within 50% of the true value.
2. **Fix the offset parameter.** Set `offset = <O>_ss = tr(gibbs * O)` (known analytically) instead of fitting it. This reduces the problem from 3 parameters to 2, greatly improving convergence.
3. **Bound the parameters.** Use LsqFit.jl's `lower` and `upper` bounds: `gamma in [1e-4, 10.0]`, `A in [-max_signal, max_signal]`. This prevents runaway to unphysical values.
4. **Try multiple initial guesses.** Run the fit from 3-5 different `gamma_0` values spanning the expected range. Take the result with the lowest residual sum of squares.

**Detection:**
- LsqFit.jl throws a convergence warning or returns `converged = false`.
- The fitted `gamma` is negative or orders of magnitude from expected.
- The fit residual is comparable to the total signal variance (fit did not improve over a constant model).

**Phase to address:** Exponential fitting implementation.

---

### Pitfall 6: Using abs(spectral_gap) Instead of Sorting by Real Part

**What goes wrong:**
The existing `run_lindbladian()` (furnace.jl line 18) sorts eigenvalues by `abs.(real.(eigvals_near_zero))` to find the steady state (smallest) and gap mode (second smallest). This works when eigenvalues are real or nearly real. But the Arpack `eigs()` call with `sigma=shift` (shift-invert mode) returns the 2 eigenvalues closest to `shift = 1e-9 * (1 + 1im)`. The sorting by `abs(real(...))` could misidentify the gap eigenvalue if:
- Two eigenvalues have similar real parts but different imaginary parts.
- The shift-invert targets an eigenvalue with a small imaginary part but not the one with the smallest real part magnitude.

For KMS BohrDomain where all eigenvalues are real, this is fine. For TrotterDomain where eigenvalues have small imaginary parts, sorting by `abs(real(...))` could pick the wrong second eigenvalue, giving a wrong spectral gap.

**Consequences:**
- The "exact" spectral gap used as ground truth in cross-validation is itself wrong.
- The cross-validation comparison is meaningless because both sides are wrong.

**Prevention:**
1. After extracting the spectral gap, verify: compute `gamma_exact = -real(spectral_gap)` and `omega_exact = imag(spectral_gap)`. Assert `gamma_exact > 0` (the gap eigenvalue must have negative real part). Assert `gamma_exact < gamma_max` for some reasonable bound (e.g., `gamma_max < 1.0` for the rescaled Hamiltonian).
2. For the cross-validation study, also compute the full spectrum (not just 2 eigenvalues) for n=4 systems where this is feasible (Liouvillian is 256x256). Use `eigen()` instead of `eigs()` and verify that the gap identified by `eigs()` matches the second-smallest `abs(real(...))` in the full spectrum.
3. Consider requesting more eigenvalues from Arpack: `eigs(liouv, nev=5, sigma=shift)` and then sorting to find the true gap.

**Detection:**
- For BohrDomain: full eigenvalue check (all should be real and non-positive).
- For n=4: compare `eigs(nev=2)` result against `eigen()` full spectrum.
- If the fitted gap consistently disagrees with the "exact" gap in a systematic direction, suspect that the exact gap is wrong.

**Phase to address:** Cross-validation phase -- validate the ground truth before using it.

---

### Pitfall 7: Wrong Time Grid for Observable Sampling

**What goes wrong:**
The existing `run_trajectories()` with observables uses `save_every` to subsample the time grid: measurements are taken every `save_every * delta` time units. If `save_every` is too large, the Nyquist criterion for the oscillation frequency `omega` (imaginary part of the spectral gap) is violated, and the time series aliases. More commonly, if `save_every` is too small, the data points are highly autocorrelated (consecutive measurements on the same trajectory step are correlated), which inflates apparent signal-to-noise and leads to underestimated error bars on the fit.

For spectral gap fitting specifically:
- The TOTAL evolution time `total_time` must be long enough to see at least 3-5 decay times: `total_time > 5 / gamma_1`. For `gamma_1 ~ 0.08`, this means `total_time > 62`. With `delta = 0.01`, that is 6,200 steps.
- The sampling interval `save_every * delta` must be fine enough to resolve the decay: at least 10-20 points per decay time. For `gamma_1 ~ 0.08`, one decay time is `~12.5`, so `save_every * delta ~ 0.6-1.2`.
- But `delta = 0.01` with `save_every = 1` gives 6,200 data points, which is excessive and wastes memory. Use `save_every = 50-100`.

**Why it happens:**
The time grid parameters were designed for convergence tracking (where you want to see the full trajectory), not for spectral gap estimation (where you want a specific time range at specific resolution).

**Consequences:**
- Too coarse sampling: miss the initial transient, aliasing of oscillatory components.
- Too fine sampling: memory bloat (6200 points * N_traj * N_obs), slow fitting, autocorrelation bias.
- Too short total time: cannot resolve the spectral gap (the signal has not decayed enough).
- Too long total time: the late-time data is all noise (see Pitfall 4), wasting computation.

**Prevention:**
1. Choose `total_time` based on a rough estimate of the spectral gap. For cross-validation where the exact gap is known, use `total_time = 5 / gamma_exact`. Otherwise, use a conservative estimate based on the domain approximation level.
2. Choose `save_every` so that there are 100-500 usable data points across the signal region.
3. For memory efficiency, use `save_every >> 1`. The trajectory step loop runs at `delta` resolution for accuracy, but measurements can be coarser.
4. Run a short pilot trajectory (100 trajectories, full time range) to estimate the decay timescale before committing to the full run.

**Detection:**
- The observable time series is flat (no visible decay) -- total time too short.
- The observable time series is noisy everywhere -- too few trajectories or total time too long.
- Memory issues when storing `measurements_mean` -- `save_every` too small.

**Phase to address:** Trajectory runner phase -- parameter selection.

---

### Pitfall 8: Naive Error Estimation on Fitted Parameters

**What goes wrong:**
LsqFit.jl provides `stderror(fit)` which estimates parameter uncertainties from the covariance matrix of the fit. This estimate assumes: (a) the model is correct, (b) residuals are independent and identically distributed, and (c) the noise is Gaussian. For trajectory-averaged observables, condition (b) is violated: consecutive time points from the same trajectory are correlated (they share the same trajectory history), and the noise is non-Gaussian (it comes from projective measurements on quantum states, not Gaussian noise).

**Consequences:**
- `stderror(fit)` dramatically underestimates the true uncertainty on the fitted spectral gap.
- The cross-validation reports "agreement within error bars" when the actual uncertainty is much larger.
- Publication-quality error bars based on `stderror` are unreliable.

**Prevention:**
1. **Use bootstrap resampling over trajectories.** The correct error estimation procedure is:
   - Split the N_traj trajectories into M bootstrap samples (resample with replacement).
   - For each bootstrap sample, recompute the trajectory-averaged observable time series.
   - Fit each bootstrap time series independently.
   - The standard deviation of the M fitted gaps is the bootstrap standard error.
2. **Use block averaging.** Divide trajectories into K blocks. Compute the fitted gap from each block independently. The standard error is `std(block_gaps) / sqrt(K)`.
3. **Do NOT use jackknife on time points** (this would propagate autocorrelation). Resample over trajectories (which are independent), not over time points.
4. Bootstrap is computationally cheap: the expensive part is running trajectories. Refitting 100 bootstrap samples is trivial.

**Detection:**
- `stderror(fit)` gives uncertainties much smaller than the observed variation across different seeds or trajectory counts.
- Bootstrap uncertainty is 5-50x larger than `stderror(fit)`.

**Phase to address:** Statistical analysis phase -- error estimation on the spectral gap.

---

### Pitfall 9: Subtracting Wrong Steady-State Value in Observable Decay

**What goes wrong:**
The exponential decay analysis requires subtracting the steady-state expectation value: `signal(t) = <O>(t) - <O>_ss`. If `<O>_ss` is wrong, the entire decay analysis is corrupted. There are several ways to get `<O>_ss` wrong:

1. **Using the Gibbs state instead of the Liouvillian fixed point.** For exact KMS (BohrDomain), the fixed point IS the Gibbs state. But for approximate domains (Energy, Time, Trotter), the Liouvillian fixed point deviates from the Gibbs state. The observable decays toward the FIXED POINT, not toward the Gibbs state. Using `<O>_gibbs` instead of `<O>_fixed_point` introduces a constant offset that the exponential fitter will absorb, biasing the decay rate.

2. **Using `tr(gibbs_comp * O_comp)` in the computational basis when `O` is in the eigenbasis.** This gives the wrong trace. Must use `tr(gibbs_eigen * O_eigen)` or `tr(gibbs_trotter * O_trotter)`, consistent with the basis of the trajectory.

3. **Using the trajectory time-averaged value instead of the known analytical value.** The trajectory-averaged `<O>(t_final)` at the end of a long evolution approximates `<O>_ss` but with statistical noise. Using it as `<O>_ss` introduces noise into every data point of the decay signal, creating correlated errors.

**Consequences:**
- A systematic offset in the signal causes the exponential fit to converge to a wrong rate.
- For small offsets, the bias is approximately `delta_gamma ~ offset / (t_final * signal_amplitude)`, which can be significant.
- The cross-validation shows a consistent directional bias (fit always above or always below the exact gap).

**Prevention:**
1. For cross-validation where the Liouvillian is available: compute `<O>_ss = tr(fixed_point * O)` using `liouv_result.fixed_point` from `run_lindbladian()`. This is the correct steady-state value that the trajectories converge to.
2. When the Liouvillian is NOT available (large systems): use `<O>_gibbs` but acknowledge the domain approximation error. For KMS with coherent term, this error is controlled by Trotter/quadrature errors.
3. Always compute both `<O>_gibbs` and `<O>_fixed_point` and report the difference. If it exceeds the statistical precision, flag it.

**Detection:**
- The subtracted signal `<O>(t) - <O>_ss` does not approach zero at late times (it approaches a nonzero constant). This is a clear sign that `<O>_ss` is wrong.
- The fit residuals show a systematic constant offset.

**Phase to address:** Observable analysis phase -- computing the correct baseline.

---

### Pitfall 10: Total Magnetization Has Zero Overlap with Gap Mode (Symmetry Selection)

**What goes wrong:**
For the isotropic Heisenberg Hamiltonian `H = sum(XX + YY + ZZ)`, the total magnetization `M_z = sum_i Z_i` commutes with the Hamiltonian: `[H, M_z] = 0`. This means `M_z` is block-diagonal in the energy eigenbasis, with blocks corresponding to different total `S_z` sectors. If the gap mode of the Liouvillian connects states within the same `S_z` sector, then `M_z` may have ZERO overlap with the gap mode -- meaning `c_1 = 0` in the expansion `<M_z>(t) - <M_z>_ss = sum c_k exp(-gamma_k t)`, and the slowest visible decay is `gamma_2` (the second gap), not `gamma_1`.

For the DISORDERED Heisenberg chain (with the external Z-field that QuantumFurnace uses), `[H, M_z] != 0` in general, so this exact cancellation is broken. However, the overlap `c_1` may still be small if the disorder is weak, leading to a near-invisible slowest mode in the `M_z` time series.

**Why it happens:**
Symmetry-based selection rules are a fundamental feature of quantum systems. The Liouvillian preserves symmetry sectors if the jump operators respect the symmetry. For single-site Pauli jumps `{X_i, Y_i, Z_i}`, the `Z_i` jump preserves `M_z` but `X_i, Y_i` do not. So the full Lindbladian does NOT preserve `M_z`, and the gap mode generically has nonzero overlap with `M_z`. But the overlap may be small for systems close to the symmetric limit.

**Consequences:**
- Fitting `M_z` decay extracts the SECOND gap instead of the first.
- Cross-validation fails because the fitted rate is ~2x the true spectral gap.
- The failure mode is subtle: the fit quality (R^2) may be excellent, but the extracted gap is wrong.

**Prevention:**
1. **Use multiple observables and compare.** Fit the spectral gap from `M_z`, from `<H>`, from `<Z_1 Z_2>`, and from other observables. If they all give the same gap, it is likely the true spectral gap. If one gives a systematically different value, it may have selection-rule issues.
2. **Check overlap coefficients.** For small systems (n=4), compute the gap mode from the exact Liouvillian and calculate `c_1 = tr(gap_mode^dagger * O)` for each observable. Choose the observable with the largest `|c_1|`.
3. **Prefer observables that break symmetries.** Single-site `Z_i` (for site i with strong disorder) is better than `M_z` because it has overlap with all eigenmodes. Nearest-neighbor correlations `Z_i Z_{i+1}` are also generally safe.
4. **For the cross-validation: report the gap from each observable separately** and note which observables agree and which do not.

**Detection:**
- The fitted gap from `M_z` is approximately 2x the fitted gap from `Z_1 Z_2`.
- The `M_z` decay signal is much smaller in amplitude than expected from other observables.
- The overlap coefficient `c_1` (computed from exact diagonalization) is anomalously small.

**Phase to address:** Observable selection phase -- before running the full experiment.

---

## Minor Pitfalls

Mistakes that cause confusion, wasted time, or suboptimal results, but are easily fixed.

---

### Pitfall 11: Memory Blow-Up from Storing Full Observable Time Series per Trajectory

**What goes wrong:**
The current `run_trajectories()` with observables returns `measurements_mean` as an `n_obs x num_saves` matrix averaged over all trajectories. This is memory-efficient. But for bootstrap error estimation (Pitfall 8), you need per-trajectory or per-batch observable time series, not just the mean. Storing per-trajectory data requires `n_obs x num_saves x N_traj` memory, which for n=8 (17 observables, 500 time points, 10,000 trajectories) is 17 * 500 * 10000 * 8 bytes = 680 MB.

**Prevention:**
1. Store per-BATCH means, not per-trajectory. With batch_size=200 and 50 batches, you get 50 independent time series at a cost of 17 * 500 * 50 * 8 = 3.4 MB.
2. Compute bootstrap over batches (block bootstrap), not over individual trajectories.
3. The existing `run_trajectories_convergence` already uses a batch structure. Extend it to also record per-batch observable TIME SERIES (not just checkpoint values).

**Phase to address:** Data architecture phase.

---

### Pitfall 12: LsqFit.jl Convergence Failure with Default Parameters

**What goes wrong:**
LsqFit.jl's `curve_fit()` uses default tolerances that may be too tight for noisy data, or too few iterations for slow convergence. The default `maxIter=1000` is usually sufficient, but the default `x_tol=1e-8` and `g_tol=1e-12` can cause premature termination on noisy data where the gradient is dominated by noise.

**Prevention:**
1. Set explicit tolerances: `curve_fit(model, t, data, p0; maxIter=10000, x_tol=1e-6, g_tol=1e-8)`.
2. Check `fit.converged` after every fit. If false, increase iterations or relax tolerances.
3. Always provide parameter bounds via `lower` and `upper` to prevent physically nonsensical solutions.

**Phase to address:** Exponential fitting implementation.

---

### Pitfall 13: Initial State Contamination of Decay Rate

**What goes wrong:**
The initial state `psi0` determines the coefficients `c_k` in the eigenmode expansion. If `psi0` is a ground state (e.g., `psi0 = [1, 0, 0, ...]` in the eigenbasis), it may already be close to the Gibbs state, meaning `<O>(0) - <O>_ss` is small and the decay signal is weak. Alternatively, if `psi0` is an eigenstate of the observable, the initial transient may be dominated by a specific eigenmode that is not the gap mode.

**Prevention:**
1. Use the maximally mixed state as the initial DM (which corresponds to random pure states in the trajectory picture -- but NOTE: maximally mixed state is NOT a pure state and cannot be represented by a single trajectory starting from a specific `psi0`).
2. For trajectory simulations, use a fixed initial state like `psi0 = [1, 0, 0, ...]` (eigenbasis ground state) which is far from the Gibbs state at high temperature. This maximizes the decay signal.
3. Run from multiple initial states and verify that the fitted gap is consistent.

**Phase to address:** Experiment design phase.

---

### Pitfall 14: Confusing Observable Decay Time with Mixing Time

**What goes wrong:**
The "spectral gap" from the Liouvillian is the slowest decay rate, which sets the ASYMPTOTIC convergence rate. The actual mixing time (time to reach within epsilon of the steady state) depends on the initial state, the observable, and the target precision:

```
t_mix(epsilon) ~ (1/gamma_1) * log(c_max / epsilon)
```

where `c_max` is the largest eigenmode coefficient. The decay rate `gamma_1` from fitting gives the SLOPE of the log-distance curve at late times, not the total mixing time. Reporting `1/gamma_fit` as "the mixing time" without specifying the prefactor is misleading.

**Prevention:**
1. Report `gamma_fit` as the "spectral gap estimate" or "asymptotic decay rate", NOT as "the mixing time".
2. If estimating mixing time, also estimate the prefactor: `c_max ~ |<O>(0) - <O>_ss|`.
3. The mixing time bound is `t_mix ~ (1/gamma_1) * log(c_max / epsilon)` for target precision epsilon.

**Phase to address:** Results reporting -- terminology and interpretation.

---

## Phase-Specific Warnings

| Phase Topic | Likely Pitfall | Severity | Mitigation |
|-------------|---------------|----------|------------|
| **Observable construction (M_z)** | Basis mismatch (Pitfall 2) | CRITICAL | Follow `build_convergence_observables` pattern exactly |
| **Observable construction (M_z)** | Symmetry selection (Pitfall 10) | MODERATE | Use multiple observables, check overlap |
| **Trajectory runner for gap estimation** | Wrong time grid (Pitfall 7) | MODERATE | Compute time range from rough gap estimate |
| **Trajectory runner for gap estimation** | Noise floor at late times (Pitfall 4) | CRITICAL | Pre-compute signal-to-noise, choose N_traj accordingly |
| **Exponential fitting** | Single-exponential bias (Pitfall 3) | CRITICAL | Late-time fitting, log-linear pre-estimate |
| **Exponential fitting** | Initial guess sensitivity (Pitfall 5) | MODERATE | Log-linear pre-estimate, bounded parameters |
| **Exponential fitting** | LsqFit.jl convergence (Pitfall 12) | MINOR | Explicit tolerances, bounds, convergence check |
| **Cross-validation metric** | Complex vs real gap (Pitfall 1) | CRITICAL | Always use `-real(spectral_gap)`, check imaginary part |
| **Cross-validation ground truth** | Wrong eigenvalue from Arpack (Pitfall 6) | MODERATE | Full spectrum check for n=4, request more eigenvalues |
| **Cross-validation baseline** | Wrong steady-state value (Pitfall 9) | MODERATE | Use Liouvillian fixed point, not Gibbs state |
| **Error estimation** | Naive stderror (Pitfall 8) | MODERATE | Bootstrap over trajectory batches |
| **Memory management** | Per-trajectory storage (Pitfall 11) | MINOR | Per-batch storage, block bootstrap |
| **Results interpretation** | Decay rate vs mixing time (Pitfall 14) | MINOR | Correct terminology, include prefactor |

---

## Recommended Build Order to Minimize Risk

Based on the pitfall analysis, the safest implementation order is:

1. **Cross-validation metric first** (addresses Pitfalls 1, 6): Define exactly what you are comparing. Extract exact gap from Liouvillian. Validate with full spectrum for n=4.

2. **Observable construction second** (addresses Pitfalls 2, 10): Build M_z in the correct basis. Compute overlap coefficients with gap mode for n=4 to verify observability.

3. **Trajectory runner with correct time grid third** (addresses Pitfalls 4, 7): Compute required total_time and save_every from the exact gap. Run pilot to verify signal-to-noise.

4. **Exponential fitting fourth** (addresses Pitfalls 3, 5, 12): Implement log-linear pre-estimate, late-time fitting, fixed baseline subtraction.

5. **Error estimation fifth** (addresses Pitfall 8): Implement block bootstrap over trajectory batches.

6. **Full cross-validation last** (addresses Pitfall 9): Compare trajectory-derived gap against exact gap with proper error bars.

---

## Sources

### Verified (HIGH confidence)
- QuantumFurnace.jl codebase -- Direct analysis of `furnace.jl` (run_lindbladian, spectral_gap extraction), `structs.jl` (LindbladianResult with Complex{T} spectral_gap), `trajectories.jl` (observable accumulation, TrajectoryWorkspace), `convergence.jl` (build_convergence_observables basis transforms), `qi_tools.jl` (gibbs_state, trace_distance)
- Quick task 20 summary (.planning/quick/20-debug-gns-trotterdomain-0-83-gap-suspect/20-SUMMARY.md) -- Documented basis mismatch causing spurious 0.83 gap (should be 0.0807)
- [LsqFit.jl documentation](https://julianlsolvers.github.io/LsqFit.jl/latest/tutorial/) -- Levenberg-Marquardt, parameter bounds, convergence control
- [LsqFit.jl GitHub](https://github.com/JuliaNLSolvers/LsqFit.jl) -- API reference, weighted fitting
- [Lindbladian Wikipedia](https://en.wikipedia.org/wiki/Lindbladian) -- Spectral gap definition, eigenvalue structure, convergence rate
- [Sandvik (2011) "Excitation Gap from Optimized Correlation Functions in QMC Simulations"](https://ar5iv.labs.arxiv.org/html/1112.2269) -- Signal-to-noise in extracting gaps from Monte Carlo
- [Nachtergaele, Sims (2006) "Spectral Gap and Exponential Decay of Correlations"](https://link.springer.com/article/10.1007/s00220-006-0030-4) -- Theory of gap-controlled correlation decay
- [Mori (2022) "Liouvillian analysis of relaxation time in open quantum systems"](https://www2.yukawa.kyoto-u.ac.jp/~nqs2022/slide/4th/Mori.pdf) -- Complex eigenvalue structure, decay rate vs oscillation frequency
- [Chen, Kastoryano, Gilyen (2025)](https://arxiv.org/abs/2311.09207) -- KMS detailed balance, real eigenvalues under GNS inner product

### Domain knowledge (HIGH confidence, established physics/numerics)
- Exponential fitting of sums of exponentials is an ill-conditioned problem (classic numerical analysis result)
- Late-time fitting extracts the slowest decay mode (standard technique in spectroscopy and QMC)
- Bootstrap resampling is the gold standard for error estimation in Monte Carlo (Efron 1979)
- Signal-to-noise for Monte Carlo averages scales as 1/sqrt(N_samples)
- KMS-detailed-balanced Lindbladians are self-adjoint w.r.t. GNS inner product, giving real eigenvalues
- Total magnetization commutes with isotropic Heisenberg Hamiltonian (SU(2) symmetry)
- Disorder breaks SU(2) symmetry, making all observables generically overlap with all eigenmodes

### Partially verified (MEDIUM confidence)
- [ResearchGate discussion on fitting sum of exponentials](https://www.researchgate.net/post/What_are_good_methods_for_fitting_a_sum_of_exponentials_to_data_without_an_initial_guess) -- Practical advice on multi-exponential fitting
- [Exponential curve fitting numerical conditioning](https://davdata.nl/math/expfitting.html) -- Ill-conditioning analysis
- [Mixing time from Liouvillian spectral gap](https://arxiv.org/html/2411.04454) -- Recent theoretical bounds on mixing time

---
*Pitfalls research for: v1.3 Mixing Time Estimation -- Spectral gap estimation from trajectory observables*
*Researched: 2026-02-16*
