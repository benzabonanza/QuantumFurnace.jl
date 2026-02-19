# Domain Pitfalls: Spectral Gap Refinement Diagnostics

**Domain:** Adding two-exponential fitting, Richardson extrapolation, effective rate plots, bootstrap error bars, anti-Hermitian defect diagnosis, symmetry sector analysis, and automatic fitting window selection to an existing quantum trajectory simulator (QuantumFurnace.jl v1.3)
**Researched:** 2026-02-19
**Confidence:** HIGH (codebase analysis of fitting.jl/gap_estimation.jl/trajectories.jl + established numerical analysis + lessons from Quick-22 through Quick-32 investigations)

**Relationship to prior research:** This document supersedes the v1.3 research PITFALLS.md (2026-02-16) for the v1.4 milestone. Many v1.3 pitfalls remain relevant (basis mismatch, noise floor, symmetry selection) but are now understood empirically. This document focuses on **new** pitfalls specific to the seven features being added in v1.4, plus integration pitfalls from connecting them to the existing system.

---

## Critical Pitfalls

Mistakes that produce silently wrong spectral gap estimates, introduce false conclusions about which error source dominates, or require rewrites.

---

### Pitfall 1: Two-Exponential Fit Parameter Non-Identifiability (Sloppy Direction)

**What goes wrong:**
The two-exponential model `f(t) = A1*exp(-g1*t) + A2*exp(-g2*t) + C` has 5 parameters. When the two decay rates g1 and g2 are within a factor of ~2 of each other, the model becomes structurally non-identifiable: many different (A1,g1,A2,g2) combinations produce nearly identical residuals. The Levenberg-Marquardt optimizer (LsqFit.jl `curve_fit`) converges to whichever local minimum is closest to the initial guess, and the resulting g1 can vary by 50-300% depending on initialization.

**Why it happens:**
This is a fundamental property of sums of exponentials, not a software bug. The condition number of the Jacobian matrix at the solution grows as `~1/(g2-g1)^2` when the rates are close. For the n=4 Heisenberg chain, the exact spectral gap is g1=0.173 and the second decay rate is g2~0.35 (ratio ~2x), which is right at the boundary of identifiability. The Hessian has a "sloppy direction" where A1 and A2 can trade off against g1 and g2 with almost no change in the cost function.

**Specific manifestation in QuantumFurnace:** The existing `fit_exponential_decay` in fitting.jl uses `_log_linear_initial_guess` which estimates parameters from the tail (lines 73-101). For a two-exponential model, this approach gives the SLOW rate reliably but gives an unreliable fast rate, because the tail is dominated by the slow component. The fast component's amplitude and rate are determined almost entirely by the initial transient, which has low signal-to-noise in trajectory data.

**How to avoid:**
1. **Constrain g2 > g1 explicitly.** Without this constraint, the optimizer can swap the two components, doubling the effective number of local minima. Use LsqFit.jl parameter bounds: `lower = [-Inf, 0.0, -Inf, ?, -Inf]` where `?` is `g1 + epsilon` (requires sequential fitting).
2. **Two-stage initialization.** First fit single-exponential to the tail (t > 0.5*T) to get g1_init. Then fit single-exponential to `data - A1_init*exp(-g1_init*t)` in the head (t < 0.3*T) to get g2_init. Use (A1_init, g1_init, A2_init, g2_init, C_init) as the starting point.
3. **Profile likelihood for g1.** Fix g1 at a grid of values, optimize over (A1, A2, g2, C) for each, and plot the profile likelihood. The minimum of this profile gives g1 with reliable uncertainty even when the full 5D landscape is sloppy.
4. **Accept that g2 is unreliable.** For spectral gap estimation, g1 is what matters. Report g2 with an explicit "LOW confidence" flag and do NOT use g2 for Richardson extrapolation or other downstream calculations.
5. **Implement a separation test.** If the fitted g2/g1 < 1.5, flag the fit as "rates too close for reliable two-exponential decomposition." Fall back to single-exponential tail fit.

**Warning signs:**
- LsqFit.jl `stderror(fit)` reports standard errors on g1 and g2 that are comparable to or larger than the gap between them.
- Running the fit from 5 different initial guesses gives 5 different (g1, g2) pairs with similar residuals.
- The covariance matrix `vcov(fit)` has a condition number > 10^6.
- A1 and A2 have opposite signs (the optimizer is using cancellation to achieve a better fit).

**Phase to address:** Two-exponential fitting implementation (early phase). Must be the FIRST thing validated before using two-exp fits for any downstream diagnostic.

---

### Pitfall 2: Anti-Hermitian Defect Blow-Up from rho^{-1/4} at Low Temperature

**What goes wrong:**
The KMS similarity transform uses `rho^{1/4}` and `rho^{-1/4}` to map the Lindbladian to a self-adjoint operator: `L_tilde = rho^{-1/4} * L * rho^{1/4}`. For the Gibbs state `rho = exp(-beta*H) / Z`, the eigenvalues of `rho` are `exp(-beta*E_k) / Z`. At low temperature (beta=10, n=6), the ratio of the smallest to largest eigenvalue of `rho` is `exp(-beta*(E_max - E_min))`. For the Heisenberg chain with n=6, the bandwidth is ~12J, so this ratio is `exp(-120) ~ 10^{-52}`. The `rho^{-1/4}` transform amplifies this by the -1/4 power: `10^{-52*(-1/4)} = 10^{13}`. This means a matrix element of L that couples a high-energy state to a low-energy state gets amplified by a factor of ~10^{13}.

**Why it happens:**
The similarity transform is mathematically well-defined (rho is positive definite, so rho^{-1/4} exists) but numerically catastrophic. In Float64, machine epsilon is ~10^{-16}. A 10^{13} amplification means that numerical noise at the 10^{-16} level becomes noise at the 10^{-3} level in L_tilde. The anti-Hermitian defect `||L_tilde - L_tilde^dagger|| / ||L_tilde||` will be dominated by this amplified noise, NOT by physical violation of KMS detailed balance.

**Specific manifestation in QuantumFurnace:** The existing `gibbs_state_in_eigen` (qi_tools.jl line 198) computes the Gibbs state correctly. But `hamiltonian.gibbs` is stored as a `Hermitian` matrix. Computing `rho^{-1/4}` requires eigendecomposition: `rho = V * D * V'`, `rho^{-1/4} = V * D^{-1/4} * V'` where `D^{-1/4} = diag(d_k^{-1/4})`. For d_k ~ 10^{-52}, d_k^{-1/4} ~ 10^{13}. This is computed correctly in principle but the subsequent matrix multiplications lose all precision for components involving small eigenvalues of rho.

**Consequences:**
- The normality ratio `||[L_tilde, L_tilde^dagger]|| / ||L_tilde||^2` reports a large anti-Hermitian defect even for the exact KMS Lindbladian (BohrDomain), which is known to be exactly self-adjoint in the GNS inner product. This leads to the false conclusion that the Lindbladian construction has a bug.
- Any diagnostic that relies on L_tilde (e.g., verifying that eigenvalues of L_tilde are real) will be corrupted.

**How to avoid:**
1. **Truncate the Gibbs state spectrum.** When computing `rho^{-1/4}`, set a floor on the eigenvalues: `d_k_floor = max(d_k, epsilon_floor)` where `epsilon_floor ~ 1e-12`. This limits the amplification to `(1e-12)^{-1/4} ~ 5600`, which is manageable. Document that this truncation limits the diagnostic to the "thermally accessible" subspace.
2. **Work in the projected subspace.** Identify the k eigenvalues of rho with `d_k > threshold` (e.g., threshold = 1e-10 * max(d_k)). Project L onto this k-dimensional subspace, then compute the similarity transform in the reduced space. This is mathematically cleaner than truncation.
3. **Use the normality ratio relative to the projected norm.** Report `||anti-Hermitian part of L_tilde_projected|| / ||L_tilde_projected||` where both norms are in the projected subspace.
4. **Validate against BohrDomain.** For BohrDomain (exact KMS), the anti-Hermitian defect after projection should be at machine precision (~1e-14). If it is not, the projection threshold is too aggressive. For TrotterDomain, the defect measures the actual KMS violation from Trotterization.
5. **Report the condition number** `max(d_k)/min(d_k)` and warn when it exceeds 10^{10}.

**Warning signs:**
- Anti-Hermitian defect > 1e-6 for BohrDomain (should be ~0 for exact KMS).
- The defect INCREASES when switching from TrotterDomain to BohrDomain (impossible if code is correct).
- `rho^{-1/4}` matrix has entries > 10^{10}.

**Phase to address:** Anti-Hermitian defect diagnosis phase. The truncation/projection strategy must be designed BEFORE computing any similarity transforms.

---

### Pitfall 3: Effective Rate lambda_eff(t) Divergence at Sign Changes

**What goes wrong:**
The effective rate is defined as `lambda_eff(t) = -d/dt log|Delta(t)|` where `Delta(t) = <O>(t) - <O>_ss` is the observable deviation from steady state. In practice this is computed as `lambda_eff(t_i) = -log|Delta(t_{i+1})/Delta(t_i)| / dt`. When `Delta(t)` passes through zero (sign change due to noise or oscillation), the ratio `Delta(t_{i+1})/Delta(t_i)` becomes negative, and `log(negative)` is undefined. Even near a zero crossing, `|Delta(t_{i+1})/Delta(t_i)|` can be astronomically large or small, producing NaN, +Inf, or -Inf in lambda_eff.

**Why it happens:**
For trajectory-averaged observables, `Delta(t)` is the mean over N trajectories. At late times when `|Delta(t)|` is comparable to the statistical noise `sigma/sqrt(N)`, the sign of Delta(t) fluctuates randomly between consecutive time points. This happens at `t ~ (1/gap) * log(A * sqrt(N) / sigma)`. For n=4 with gap~0.17, A~1, N=20000, sigma~1: `t ~ 6 * log(141) ~ 30`. Beyond t~30, lambda_eff is dominated by noise.

Even at earlier times, the observable may have oscillatory components (from non-zero imaginary parts of Lindbladian eigenvalues, see v1.3 Pitfall 1) that cause Delta(t) to change sign. Multi-exponential decay where the fast component has a different sign from the slow component also causes a zero crossing.

**Specific manifestation in QuantumFurnace:** The existing trajectory infrastructure returns `measurements_mean` as a matrix (convergence.jl). The effective rate plot requires computing `Delta(t) = measurements_mean[i, :] .- O_ss` and then the discrete derivative of `log|Delta|`. A naive implementation will produce NaN at every zero crossing.

**Consequences:**
- The effective rate plot has gaps or spikes that look like bugs but are physically meaningful (they indicate the boundary of the usable fitting window).
- Automatic fitting window selection algorithms that scan lambda_eff for a plateau will crash or produce nonsensical windows if they encounter NaN/Inf.
- Plotting lambda_eff with NaN values causes gaps in plots that confuse interpretation.

**How to avoid:**
1. **Guard the log.** Use `lambda_eff(t_i) = -log(max(|Delta(t_{i+1})/Delta(t_i)|, epsilon)) / dt` where `epsilon = 1e-15`. This prevents NaN but produces a large negative lambda_eff at zero crossings.
2. **Mark invalid points explicitly.** Return lambda_eff as a vector of `Union{Float64, Nothing}` or use NaN with a companion boolean mask `valid[i]`. The plotting and window selection code must respect this mask.
3. **Compute lambda_eff from a smoothed Delta(t).** Apply a Savitzky-Golay filter or running median to Delta(t) before computing the log-derivative. This suppresses the single-point noise that causes sign flips. BUT: smoothing introduces bias in the effective rate at early times where Delta(t) changes rapidly. Use a smoothing window that is small compared to 1/gap.
4. **Use absolute value consistently.** Work with `|Delta(t)|` throughout, not `Delta(t)`. The effective rate is well-defined for `|Delta(t)|` as long as it is monotonically decreasing. If `|Delta(t)|` is not monotonically decreasing, that itself is diagnostic information (indicates multi-exponential structure or oscillations).
5. **Define a noise floor cutoff.** Compute `noise_floor(t) = std_over_trajectories(O(t)) / sqrt(N)`. Exclude time points where `|Delta(t)| < 3 * noise_floor(t)` from the lambda_eff computation entirely.

**Warning signs:**
- lambda_eff contains NaN or Inf values.
- lambda_eff has sharp spikes (>10x the median value) at isolated time points.
- lambda_eff oscillates with amplitude comparable to the plateau value.

**Phase to address:** Effective rate plot implementation. The noise floor cutoff and sign-change handling must be built into the lambda_eff computation function, not handled ad-hoc by the caller.

---

### Pitfall 4: Richardson Extrapolation Making Things Worse (Wrong Error Structure)

**What goes wrong:**
Richardson extrapolation assumes the error in the gap estimate has the form `gap_est(delta) = gap_exact + C*delta^p + O(delta^{p+1})` for some known order p (typically p=1 or p=2). The formula then cancels the leading error term. Quick-30 conclusively proved that the gap estimation error does NOT follow this structure: the error/delta ratio varied by 96x across delta=0.1, 0.01, 0.001, and Richardson extrapolation provided only 1.0x improvement. Applying Richardson extrapolation when the assumption is violated can produce results that are WORSE than the un-extrapolated estimate.

**Why it happens:**
The gap estimation error has at least three sources with different delta-dependence:
1. **Trotter discretization error:** O(delta) per trajectory step, accumulates over T/delta steps, but this affects the trajectory-averaged density matrix, NOT the fitted gap rate directly. Quick-32 confirmed this error is small (~1e-3) and does decrease monotonically with delta.
2. **Fitting model error:** The single-exponential fit captures a weighted average of multiple decay modes. Different delta values produce different effective Kraus channels that excite spectral modes with different weights, causing the "effective gap" from fitting to vary non-monotonically with delta. This is delta-independent in magnitude (~37-49% for Mz_stagg) but delta-dependent in direction.
3. **Statistical noise:** Adds O(1/sqrt(N)) noise to the gap estimate, uncorrelated with delta.

Richardson extrapolation cancels only source (1), which is already the smallest contributor. Source (2), which dominates, has NO clean power-law dependence on delta. Applying Richardson amplifies source (2) because the formula `gap_rich = (delta2 * gap(delta1) - delta1 * gap(delta2)) / (delta2 - delta1)` subtracts two similar numbers, amplifying their fitting-error differences.

**Specific manifestation in QuantumFurnace:** Quick-30 data (30-SUMMARY.md lines 66-72): for the (0.01, 0.1) pair, Richardson gave 49.83% error vs 48.72% for the finer estimate alone. For the (0.001, 0.01) pair, Richardson gave 35.88% vs 37.17%. The "improvement" is within noise, and for some observables Richardson makes things noticeably worse.

**How to avoid:**
1. **Do NOT apply Richardson extrapolation to single-exponential gap estimates.** This is a settled conclusion from Quick-30. The code should not even offer this as an automatic option.
2. **For two-exponential fits:** Re-evaluate whether Richardson is viable by checking if the TWO-exponential gap estimate has clean O(delta^p) error structure. Run the same delta-sweep experiment as Quick-30 but with two-exponential fitting. Only if error/delta is approximately constant across 3+ delta values should Richardson be applied.
3. **If Richardson is applied:** Always report the un-extrapolated estimates alongside the extrapolated one, and flag when the extrapolated estimate is OUTSIDE the range of the un-extrapolated estimates (which indicates the extrapolation diverged).
4. **Implement a delta-convergence diagnostic BEFORE Richardson.** Plot the fitted gap vs delta for 3+ delta values. If the curve is monotonic with consistent curvature, Richardson is appropriate. If the curve is non-monotonic or noisy, Richardson is contraindicated. This diagnostic should be a required precondition, not an optional check.

**Warning signs:**
- Richardson-extrapolated gap is outside the range [min(gap_estimates), max(gap_estimates)].
- The extrapolated gap is negative.
- The "improvement factor" (error_rich / error_finest) is > 0.9 (no improvement) or > 1.0 (degradation).

**Phase to address:** Delta-convergence diagnosis phase. The diagnostic must PRECEDE any Richardson extrapolation. The extrapolation itself should be gated on a monotonicity check.

---

### Pitfall 5: Bootstrap Storing Individual Trajectory Data Blows Memory for n=6

**What goes wrong:**
Bootstrap error bars on the spectral gap require access to per-trajectory (or per-batch) observable time series, not just the grand mean. The naive implementation stores per-trajectory data: a matrix of size `(n_obs, n_saves, n_traj)`. For n=6 (dim=64), 8 observables, save_every=10, delta=0.01, mixing_time=50 (5000 steps, 500 saves), and 20,000 trajectories: 8 * 500 * 20000 * 8 bytes = 640 MB. For 50,000 trajectories: 1.6 GB. This exceeds reasonable memory budgets for a diagnostic tool and prevents running multiple bootstrap analyses in a session.

**Why it happens:**
The existing `run_observable_trajectories` (trajectories.jl) only returns `measurements_mean`, the average over all trajectories. It does not store per-trajectory data. Adding bootstrap requires either: (a) storing all per-trajectory data, (b) re-running trajectories with different subsets, or (c) accumulating per-batch statistics during the run. Option (a) is the obvious implementation but is memory-prohibitive. Option (b) is wasteful (re-runs trajectories). Option (c) is efficient but requires modifying the trajectory runner.

**Specific manifestation in QuantumFurnace:** The existing batch infrastructure in `run_trajectories_convergence` (convergence.jl) already computes per-batch density matrices but does NOT store per-batch observable time series. The `_run_batch_no_obs!` function accumulates the density matrix sum. The `_run_chunk_with_obs!` function accumulates measurements into `mean_data_local` but averages across all trajectories in the chunk. There is no existing mechanism to return per-batch measurement arrays.

**How to avoid:**
1. **Block bootstrap over trajectory batches.** Modify the trajectory runner to run in fixed batches (e.g., batch_size=200) and store per-BATCH mean time series: `(n_obs, n_saves, n_batches)`. For 100 batches of 200: 8 * 500 * 100 * 8 bytes = 3.2 MB. Bootstrap then resamples over the 100 batch indices, recomputing the mean for each bootstrap sample as a weighted sum of batch means.
2. **Online sufficient statistics.** During the trajectory run, maintain two running accumulators: `sum_obs[i, t]` and `sum_obs_sq[i, t]` (the latter for variance estimation). This requires zero additional memory beyond 2 * n_obs * n_saves floats, but only gives the normal-approximation confidence interval, not the full bootstrap distribution.
3. **Two-pass approach.** First pass: run all trajectories, store only the grand mean. Second pass: run bootstrap samples (smaller batches) with different seeds. This avoids storing any per-trajectory data but costs ~2x the compute.
4. **Recommendation: Option 1 (block bootstrap) with batch_size=200.** This matches the existing adaptive sampling batch structure and provides genuine bootstrap distributions while keeping memory at ~3 MB. The block bootstrap is valid because trajectories within different batches are independent (different RNG seeds via `Xoshiro(master_seed + traj_id)`).

**Warning signs:**
- Julia process memory exceeds 4 GB during bootstrap (for n=6 simulations that should use ~500 MB).
- Out-of-memory crashes during bootstrap with large trajectory counts.
- Bootstrap resampling takes longer than the original trajectory run (indicates data is being re-loaded from disk).

**Phase to address:** Bootstrap error bar implementation. The per-batch storage mechanism must be designed before any bootstrap code is written. Modify `run_observable_trajectories` or create a new variant.

---

### Pitfall 6: Automatic Window Selection Fails Silently on Multi-Exponential Data

**What goes wrong:**
The "golden window" for single-exponential fitting is the time range where the fast modes have decayed but the signal is still above the noise floor. Automatic window selection looks for this window by scanning the effective rate lambda_eff(t) for a plateau. On clean single-exponential data, this works well: lambda_eff is approximately constant for t > 3/g2 and t < T_noise.

But with multi-exponential data (which is the reality for trajectory observables), lambda_eff has a monotonically DECREASING profile from g2 at early times to g1 at late times, with no clean plateau. The automatic selector either: (a) picks the entire time range (too wide, biased by fast modes), (b) picks the late-time tail (too narrow, dominated by noise), or (c) declares failure and returns no window.

**Why it happens:**
The concept of a "golden window" assumes a single dominant decay mode in the data. For the n=4 Heisenberg chain, the spectral gap mode has overlap with most observables, but the second decay mode also has significant overlap (Quick-30 showed this causes 37-49% single-exponential fitting error). The transition from the "fast mode dominated" regime to the "slow mode dominated" regime is gradual, spanning several decay times. There is no sharp boundary.

**Consequences:**
- Silent failure: the automatic selector returns a window, but it is not the "right" window. The resulting gap estimate is biased, and the user does not know the window was problematic.
- Window selection depends on noise realization: running with different seeds produces different windows, which produce different gap estimates. This inflates the apparent bootstrap uncertainty but does not capture the systematic bias from window placement.
- The two-exponential fit (which is supposed to eliminate the need for careful window selection) is itself sensitive to the fitting window (see Pitfall 1), creating a circularity.

**How to avoid:**
1. **Do NOT rely on automatic window selection as the sole method.** Always offer manual override and report which window was selected so the user can verify.
2. **Implement multiple window strategies and compare.** Run fits with: (a) skip_initial=0.1, (b) skip_initial=0.3, (c) skip_initial=0.5, and (d) the noise-floor-determined window. Report all four gap estimates. If they agree within bootstrap error bars, the result is robust. If they disagree significantly, flag the result as window-dependent.
3. **Use the effective rate plot as a DIAGNOSTIC, not as an input to automatic selection.** Plot lambda_eff(t) and let the user identify the plateau visually. The automatic selector is a convenience, not a substitute for judgment.
4. **For two-exponential fits: use the full signal window.** The whole point of two-exponential fitting is to capture both modes simultaneously, eliminating the need for window selection. Use skip_initial=0.0 (or a very small value to avoid initial-state artifacts) and fit the full time range. If the two-exponential fit converges, it implicitly handles the mode separation.
5. **Plateau detection with explicit quality metric.** If implementing automatic selection: compute lambda_eff in rolling windows of width W. For each window position, compute the coefficient of variation (std/mean) of lambda_eff within the window. The "best" window is where this CV is minimized. Report the CV value so the user knows how flat the plateau actually is. A CV > 0.2 means there is no reliable plateau.

**Warning signs:**
- The automatic selector returns a window spanning less than 2 decay times (not enough data for reliable fitting).
- Different trajectory seeds produce windows that differ by more than 30%.
- The gap estimate changes by more than 20% when the window is shifted by +/- 10%.

**Phase to address:** Automatic window selection implementation. Must be designed alongside, not after, the effective rate plot.

---

### Pitfall 7: Symmetry Sector Labels Ambiguous Near Degeneracies

**What goes wrong:**
Symmetry sector labeling assigns quantum numbers (like total Sz, crystal momentum k, or parity) to Lindbladian eigenvalues. This requires identifying which symmetry sector each eigenvector belongs to by computing overlap with symmetry projectors. When two eigenvalues are nearly degenerate (|lambda_i - lambda_j| / |lambda_i| < 1e-6), the corresponding eigenvectors are arbitrary linear combinations within the degenerate subspace, and the symmetry labels become ambiguous or wrong.

**Why it happens:**
The `eigen()` decomposition of the Liouvillian returns orthogonal eigenvectors, but within a degenerate subspace, any orthonormal basis is valid. Julia's `eigen()` (LAPACK) picks a basis that is numerically convenient, not one that respects the symmetry. If eigenvalues lambda_1 and lambda_2 differ by less than the numerical precision of the eigendecomposition (~1e-12 for a well-conditioned matrix), the returned eigenvectors v1 and v2 will be random linear combinations of the true symmetry eigenstates. Applying a symmetry projector P_k to v1 will give a non-zero overlap for MULTIPLE symmetry sectors, making the assignment ambiguous.

**Specific manifestation in QuantumFurnace:** The existing `eigenbasis_overlap_analysis` (gap_estimation.jl line 279) uses `eigen(L)` and sorts by `abs.(real.(F.values))`. For the n=6 periodic Heisenberg chain, Quick-25 showed the gap mode lives in the k=pi momentum sector. The second gap mode (lambda_3) may be in a different sector. If lambda_2 and lambda_3 are close, their sector labels will be mixed up.

**Consequences:**
- The symmetry sector analysis reports contradictory labels for nearly degenerate eigenvalues.
- The "gap mode sector" diagnostic incorrectly identifies which sector the gap mode belongs to, leading to wrong observable selection recommendations.
- Users may conclude that a symmetry-breaking observable is needed when the symmetry assignment itself is wrong.

**How to avoid:**
1. **Detect near-degeneracies explicitly.** Before labeling, compute pairwise separations `|lambda_i - lambda_j|` for eigenvalues near the gap. Flag any pair with separation < threshold (e.g., 1e-8 * |lambda_gap|) as "ambiguous sector assignment."
2. **Use simultaneous diagonalization.** Instead of diagonalizing L alone, diagonalize L together with the symmetry operator (e.g., the translation operator T). The simultaneous eigendecomposition picks a basis that respects the symmetry even within degenerate subspaces. This requires that [L, T] = 0 (which Quick-25 verified for the translation operator).
3. **Apply symmetry projectors to the degenerate SUBSPACE, not individual eigenvectors.** If lambda_i and lambda_j are nearly degenerate, project the 2D subspace span(v_i, v_j) onto each symmetry sector. Report "this degenerate subspace contains modes from sectors k_a and k_b" rather than "v_i is in sector k_a."
4. **For practical purposes: label the GAP mode first.** The gap mode is typically non-degenerate (there is a unique slowest decay rate). Label it unambiguously. Only flag the degeneracy issue for higher modes.

**Warning signs:**
- Two adjacent eigenvalues differ by less than 1e-8 * |gap|.
- The symmetry projector gives overlap > 0.1 with multiple sectors for a single eigenvector.
- The sector label changes when a small perturbation (random field h ~ 1e-8) is added to the Hamiltonian.

**Phase to address:** Symmetry sector labeling implementation. Degeneracy detection must be the first step, before any sector assignment.

---

## Moderate Pitfalls

Mistakes that cause significant debugging time or wrong intermediate results but are recoverable.

---

### Pitfall 8: Two-Exponential Fit Converges to Negative Amplitude (Cancellation Artifact)

**What goes wrong:**
The two-exponential model `A1*exp(-g1*t) + A2*exp(-g2*t) + C` has no physical constraint requiring A1, A2 > 0. The optimizer may find that the best fit has A1 > 0 and A2 < 0 (or vice versa), where the two exponentials partially cancel. This produces a curve that looks superficially like an exponential decay but is actually a difference of two exponentials. The fitted g1 in this case does NOT correspond to the spectral gap -- it can be arbitrarily different.

**Why it happens:**
For trajectory-averaged observables, the expansion `<O>(t) - <O>_ss = sum_k c_k * exp(lambda_k * t)` can have coefficients c_k of either sign, depending on the observable and initial state. If c_1 > 0 and c_2 < 0, then A1 > 0 and A2 < 0 is the physically correct decomposition. But if both c_1, c_2 > 0 and the data has noise, the optimizer might find a spurious solution with one negative amplitude that exploits cancellation to fit noise features.

**How to avoid:**
1. **Constrain amplitudes to be positive** (`lower = [0.0, 0.0, 0.0, 0.0, -Inf]`) if the physics guarantees positive overlap coefficients. For the observables in `build_preset_trajectory_observables`, this is NOT guaranteed -- overlaps can have either sign.
2. **Instead: flag but do not reject negative amplitudes.** If one amplitude is negative, check whether the other amplitude exceeds the negative one in magnitude at t=0 (which is required for a physically decaying signal). Report "mixed-sign amplitudes" as a diagnostic.
3. **Detect cancellation.** Compute `A1*exp(-g1*t) + A2*exp(-g2*t)` at t=0 and compare to `|A1| + |A2|`. If `|A1 + A2| < 0.3 * (|A1| + |A2|)`, there is significant cancellation. The fitted rates in this case are unreliable.

**Warning signs:**
- A1 and A2 have opposite signs.
- `|A1| + |A2|` is much larger than the observed signal amplitude `|Delta(0)|`.
- The residuals have a characteristic "W" shape (two crossings instead of one).

**Phase to address:** Two-exponential fitting validation tests.

---

### Pitfall 9: Bootstrap Error Bars on lambda_eff(t) Require Per-Batch Time Series

**What goes wrong:**
Computing bootstrap error bars on the effective rate plot `lambda_eff(t)` requires resampling per-batch observable time series and recomputing lambda_eff for each bootstrap sample. This is a different computation from bootstrap on the fitted gap (which resamples the gap fit). If the implementation only supports bootstrap on the gap parameter, the user cannot get uncertainty bands on the lambda_eff plot.

**Why it happens:**
The lambda_eff computation involves a nonlinear transformation (log of the ratio) applied to each time point independently. Bootstrapping the fitted gap gives a scalar confidence interval. Bootstrapping lambda_eff(t) gives a function-valued confidence band. These are different quantities and require different implementations.

**How to avoid:**
1. **Design the bootstrap infrastructure to return per-batch time series**, not just per-batch gap estimates. This is the same data structure needed for Pitfall 5 (memory management). Store `batch_measurements[batch_idx][obs_idx, time_idx]` and use it for both gap bootstrap and lambda_eff bootstrap.
2. **Compute lambda_eff for each bootstrap sample.** For each resampled set of batch indices, compute the resampled mean time series, then compute lambda_eff on that resampled mean. Collect N_bootstrap lambda_eff curves and compute pointwise percentiles for the confidence band.
3. **Handle NaN propagation.** Some bootstrap samples will have more noise than others, producing NaN in lambda_eff at earlier time points. The confidence band computation must handle NaN: use `nanquantile` or equivalent.

**Warning signs:**
- Bootstrap confidence band on lambda_eff is missing (only scalar gap CI is implemented).
- The confidence band has gaps at time points where some bootstrap samples produced NaN.

**Phase to address:** Bootstrap implementation. Design the data flow for both gap bootstrap and lambda_eff bootstrap simultaneously.

---

### Pitfall 10: Testing Two-Exponential Fits Against Stochastic Simulation Output

**What goes wrong:**
Tests for the two-exponential fitting function need test data. Using synthetic data (known amplitudes, rates, noise) validates the fitting algorithm but does NOT test the full pipeline (trajectory simulation -> fitting). Using actual trajectory data makes tests non-deterministic (different seeds give different results) and slow (trajectory simulation dominates test time). The temptation is to test only with synthetic data, which misses integration issues like basis mismatches, wrong steady-state subtraction, or observable normalization errors.

**Why it happens:**
The fitting function (`fit_exponential_decay`, fitting.jl) is a pure numerical function: it takes (times, values) and returns a FitResult. Testing it with synthetic data is straightforward and fast. The integration test (trajectory -> observables -> fitting -> gap) requires building the full system (Hamiltonian, Lindbladian, Kraus operators) and running trajectories, which takes 30-60 seconds for n=4 with 1000 trajectories.

**Specific manifestation in QuantumFurnace:** The existing test suite has both: test_fitting.jl tests the fitting function with clean and noisy synthetic data (fast, deterministic), and test_gap_estimation.jl tests the full pipeline with the SMALL system (n=3, 500 trajectories, fast). But the n=3 system is too small and simple to exercise two-exponential behavior -- its spectral structure is dominated by a single decay mode. The n=4 system is needed but makes tests take ~30 seconds.

**How to avoid:**
1. **Three-tier testing strategy:**
   - **Unit tests (fast, synthetic):** Test the two-exponential fitting function with known parameters. Include edge cases: rates close together, one amplitude zero, negative amplitudes, noisy data. These run in <1 second.
   - **Integration tests (medium, n=3):** Test the full pipeline with the SMALL system. Even though n=3 may not show two-exponential behavior clearly, it validates that the data flows correctly from trajectory simulation to fitting.
   - **Validation tests (slow, n=4, optional/CI-only):** Run the full pipeline with n=4, 1000 trajectories, seed=42. Compare the fitted gap against the exact Liouvillian gap. Mark these as `@testset "slow"` and skip in regular test runs. These are the critical tests but they take 30+ seconds.
2. **Generate reference two-exponential data from the Liouvillian.** For n=4, compute `exp(t*L) * vec(rho0)` exactly and extract the observable time series. This gives clean two-exponential data (no stochastic noise) that can be used in fast deterministic tests. The reference data can be pre-computed and stored.
3. **Use deterministic seeding (seed=42) in integration tests** and test for consistency (same seed -> same result) rather than accuracy (fitted gap close to exact gap). This is already the pattern in test_gap_estimation.jl.

**Warning signs:**
- Two-exponential fitting tests all use synthetic data but fail on actual trajectory data.
- Integration tests pass but validation against exact Liouvillian fails.
- Tests are flaky (pass sometimes, fail other times) because they depend on stochastic simulation output without proper tolerance handling.

**Phase to address:** Every phase that adds a new diagnostic. Each feature needs all three test tiers.

---

### Pitfall 11: Observable-Level vs DM-Level Error Conflation

**What goes wrong:**
The delta-convergence diagnostic compares gap estimates at different delta values. But the "gap estimation error" has two components: (a) the trajectory simulation error (how well the simulated rho matches exp(t*L)*rho0) and (b) the fitting error (how well the single/two-exponential model captures the true multi-exponential decay). Quick-32 proved that component (a) is O(1e-3) and monotonic in delta, while component (b) is O(10-50%) and non-monotonic. If the delta-convergence diagnostic only reports the total error, it will conflate these two components and misattribute fitting errors to simulation errors.

**How to avoid:**
1. **Report both components separately.** For each delta value, compute: (a) `trace_distance(rho_traj, rho_exact)` as the simulation error, and (b) `|gap_fitted - gap_exact|` as the total estimation error. The difference `(b) - (a)` gives the fitting contribution.
2. **Use the Quick-32 approach.** Compare trajectory-averaged OBSERVABLE values (not just fitted gaps) against exact values from exp(t*L). This separates simulation correctness from fitting model adequacy.

**Phase to address:** Delta-convergence diagnosis phase.

---

### Pitfall 12: Steady-State Value Mismatch Between Lindbladian Fixed Point and Gibbs State for Effective Rate

**What goes wrong:**
The effective rate computation requires `Delta(t) = <O>(t) - <O>_ss`. For approximate domains (TrotterDomain), the Lindbladian fixed point differs from the Gibbs state (v1.3 Pitfall 9). Using `<O>_gibbs` instead of `<O>_fixed_point` introduces a constant offset that makes Delta(t) approach a non-zero value at late times, corrupting the effective rate plot. The lambda_eff(t) will NOT plateau at the spectral gap but will instead show a gradual decrease toward zero (as the constant offset dominates over the exponential decay).

**How to avoid:**
1. **Always use the Lindbladian fixed point** from `liouv_result.fixed_point` when available (n=4,6 where exact diagonalization is feasible).
2. **For large systems where the Liouvillian is not available:** use the trajectory-averaged observable at the final time point `<O>(T_final)` as an estimate of `<O>_ss`. This introduces noise but avoids the systematic bias from Gibbs/fixed-point mismatch. Alternatively, use the Gibbs value and document the expected offset.
3. **Diagnostic: plot `Delta(T_final)` for each observable.** If it is not approximately zero, the steady-state value is wrong.

**Phase to address:** Effective rate plot implementation. The steady-state computation must be a required input to the lambda_eff function.

---

## Minor Pitfalls

---

### Pitfall 13: Bootstrap Resampling Indices With Replacement Can Produce Degenerate Samples

**What goes wrong:**
With B batches (e.g., B=100), a bootstrap sample draws B indices with replacement from {1,...,B}. For small B, there is a non-negligible probability that a bootstrap sample consists of repeated copies of very few unique batches (e.g., 80% from one batch). The gap fitted to such a degenerate sample may have extreme values (very high or low), inflating the bootstrap confidence interval.

**How to avoid:**
- Use at least B=50 batches (100 preferred). With B=100 and standard bootstrap theory, the probability of a degenerate sample is negligible.
- Consider the block bootstrap variant: for each bootstrap sample, draw blocks of consecutive batches rather than individual batches. This preserves any temporal correlations between batches (though for this system, batches are independent by construction).

**Phase to address:** Bootstrap implementation (minor concern, just set B large enough).

---

### Pitfall 14: Plotting Effective Rate on Wrong Time Axis

**What goes wrong:**
The effective rate lambda_eff is computed from pairs of consecutive time points: `lambda_eff(t_i) = -log|Delta(t_{i+1})/Delta(t_i)| / (t_{i+1} - t_i)`. This value should be plotted at the midpoint `(t_i + t_{i+1})/2`, not at `t_i` or `t_{i+1}`. Plotting at `t_i` introduces a systematic half-step offset that shifts the apparent plateau position, biasing any automatic window selection that uses the lambda_eff time axis.

**How to avoid:**
- Return the midpoint times alongside lambda_eff values: `times_eff[i] = (times[i] + times[i+1]) / 2`.
- Document the time axis convention in the function docstring.

**Phase to address:** Effective rate implementation (minor, but easy to get wrong).

---

### Pitfall 15: Confusing bootstrap Standard Error with bootstrap Confidence Interval

**What goes wrong:**
The bootstrap gives a distribution of gap estimates. The standard error is the standard deviation of this distribution. The 95% confidence interval is the [2.5th, 97.5th] percentile of this distribution. For non-normal distributions (which is common for gap estimates, since the gap is bounded below by 0), these can differ significantly. Reporting "bootstrap SE" as the uncertainty and then computing "gap +/- 1.96*SE" assumes normality, which may not hold.

**How to avoid:**
- Report the percentile-based confidence interval directly, not a normal approximation.
- Also report the SE for comparison, but label it clearly as "bootstrap SE (normal approximation)" and the CI as "bootstrap CI (percentile method)."

**Phase to address:** Bootstrap error bar reporting.

---

## Technical Debt Patterns

Shortcuts that seem reasonable but create long-term problems.

| Shortcut | Immediate Benefit | Long-term Cost | When Acceptable |
|----------|-------------------|----------------|-----------------|
| Hardcoding `skip_initial=0.1` for all observables | Quick results, no per-observable tuning | Wrong window for observables with different decay profiles; biased gap estimates | Never -- always use per-observable window or the automatic selector |
| Storing bootstrap results only as scalar CIs | Small result structs | Cannot reconstruct per-time-point lambda_eff bands or diagnose non-normal bootstrap distributions | MVP only; full bootstrap distribution needed for publication |
| Using LsqFit.jl `stderror` for two-exponential fit uncertainties | Free with the fit, no extra computation | Dramatically underestimates uncertainty due to ill-conditioning (see Pitfall 1); misleads about parameter identifiability | Never for two-exponential fits; acceptable only for single-exponential fits as a diagnostic |
| Computing anti-Hermitian defect without rho^{-1/4} spectrum truncation | Simpler code, no threshold parameter | Misleading defect values at low temperature (see Pitfall 2); false alarm about KMS violation | Only for high-temperature (beta < 2) where condition number is manageable |
| Selecting best observable by smallest fitted gap without cross-checking | Simple selection criterion | Systematically picks observables that underestimate the gap (Quick-30: Mz_stagg at 37-49% vs YY_avg at 2-6%) | Acceptable for quick estimates; must cross-check for publication-quality results |

## Integration Gotchas

Common mistakes when connecting new v1.4 features to the existing v1.3 infrastructure.

| Integration | Common Mistake | Correct Approach |
|-------------|----------------|------------------|
| Two-exp fitting with existing `FitResult` | Trying to reuse `FitResult` struct (which has 3 params: A, gap, C) for 5 params | Create a new `TwoExpFitResult` struct with fields for both rates, both amplitudes, and offset |
| Bootstrap with `run_observable_trajectories` | Expecting per-batch data from existing function | Modify function to accept a `store_per_batch::Bool` kwarg, or create a new `run_observable_trajectories_batched` variant |
| Lambda_eff using `measurements_mean` from `ObservableTrajectoryResult` | Computing lambda_eff directly from the stored mean (which has already averaged out per-trajectory information) | This is correct for the mean lambda_eff; for bootstrap bands, need per-batch data (see Pitfall 9) |
| Symmetry sector analysis with `eigenbasis_overlap_analysis` | Assuming the existing function's eigendecomposition respects symmetry | The existing `eigen(L)` does NOT respect symmetry. Need a separate simultaneous diagonalization step (see Pitfall 7) |
| Anti-Hermitian defect with `LindbladianResult` | Using `liouv_result.liouvillian` directly for the similarity transform | The Liouvillian is in the vectorized (Liouville) form. The similarity transform `rho^{-1/4} L rho^{1/4}` acts on the operator level, requiring reshaping each column of L into a matrix, transforming, and reshaping back |
| Delta-convergence with `estimate_spectral_gap` | Running `estimate_spectral_gap` at multiple delta values and comparing | This works but conflates simulation and fitting errors (see Pitfall 11). Instead, also compare trajectory-averaged observable values against exact exp(t*L) results, as Quick-32 demonstrated |

## Performance Traps

Patterns that work at small scale but fail at larger system sizes.

| Trap | Symptoms | Prevention | When It Breaks |
|------|----------|------------|----------------|
| Dense eigendecomposition of full Liouvillian for symmetry sector analysis | OOM or multi-hour runtime | Use Arpack shift-invert for just the leading 10-20 eigenvalues; only do full eigen for n<=4 (dim^2 <= 256) | n=6: Liouvillian is 4096x4096 (dense eigen takes ~30s, OK). n=8: 65536x65536 (impossible, 32 GB) |
| Storing the full Liouvillian as a dense matrix for anti-Hermitian defect | Memory: dim^4 * 16 bytes. n=6: 4096^2 * 16 = 268 MB (OK). n=8: 65536^2 * 16 = 69 GB (impossible) | For n>6, use sparse Liouvillian and iterative methods for the defect norm; or restrict defect analysis to n<=6 | n=8 (dim=256, Liouvillian dim=65536) |
| Two-exponential fitting with multi-start (5 initial guesses x 8 observables x 100 bootstrap samples) | 4000 fits per delta value; ~20 seconds for n=4 | Reduce multi-start to 3 guesses; skip fitting for observables that failed single-exponential quality check (R^2 < 0.5) | When bootstrap sample count > 200 or observable count > 10 |
| Per-batch storage for n=8 with 1000 batches | 8 obs * 500 saves * 1000 batches * 8 bytes = 32 MB (manageable) but if save_every is too small (1 instead of 10), it becomes 320 MB | Always use save_every >= 10 for gap estimation; store only the number of batches needed for bootstrap (100 is enough) | When save_every < 5 or n_batches > 500 |

## "Looks Done But Isn't" Checklist

Things that appear complete but are missing critical pieces.

- [ ] **Two-exponential fitting:** Often missing the separation test (g2/g1 > 1.5 check) -- verify that the code flags when rates are too close for reliable decomposition
- [ ] **Bootstrap error bars:** Often missing the per-batch time series storage mechanism -- verify that the trajectory runner actually returns per-batch data, not just the grand mean
- [ ] **Effective rate plot:** Often missing the noise floor cutoff -- verify that lambda_eff is NaN/masked beyond the noise floor, not plotted with spurious values
- [ ] **Richardson extrapolation:** Often missing the precondition check (monotonic delta-convergence) -- verify that the code refuses to extrapolate when the error structure is not O(delta^p)
- [ ] **Anti-Hermitian defect:** Often missing the spectrum truncation for rho^{-1/4} -- verify that the code has a floor on Gibbs eigenvalues and reports the effective condition number
- [ ] **Symmetry sector labels:** Often missing the degeneracy detection -- verify that near-degenerate eigenvalues are flagged as "ambiguous sector"
- [ ] **Automatic window selection:** Often missing the quality metric (CV of lambda_eff in the selected window) -- verify that the code reports HOW FLAT the plateau is, not just WHERE it is
- [ ] **Testing:** Often missing validation against exact Liouvillian -- verify that at least one test compares the diagnostic output against known exact results (not just synthetic data)

## Recovery Strategies

When pitfalls occur despite prevention, how to recover.

| Pitfall | Recovery Cost | Recovery Steps |
|---------|---------------|----------------|
| Two-exp fit non-identifiable (Pitfall 1) | LOW | Fall back to single-exponential tail fit; report g1 from tail with bootstrap CI |
| rho^{-1/4} blow-up (Pitfall 2) | LOW | Recompute with spectrum truncation threshold; compare BohrDomain vs TrotterDomain defect |
| lambda_eff divergence (Pitfall 3) | LOW | Apply noise floor cutoff retroactively; recompute with smoothing |
| Richardson worsens estimate (Pitfall 4) | LOW | Discard Richardson result; use the finest-delta un-extrapolated estimate |
| Memory blow-up from per-trajectory storage (Pitfall 5) | MEDIUM | Kill the process; redesign with block bootstrap; re-run trajectory simulation with per-batch storage |
| Silent window selection failure (Pitfall 6) | MEDIUM | Re-run with multiple manual windows; compare results; identify the actual usable range from lambda_eff plot |
| Ambiguous symmetry labels (Pitfall 7) | LOW | Add small random perturbation to break degeneracy; re-label; verify labels are consistent |
| Negative amplitude in two-exp fit (Pitfall 8) | LOW | Check if the physics allows mixed-sign amplitudes; if spurious, constrain to positive amplitudes and refit |

## Pitfall-to-Phase Mapping

How roadmap phases should address these pitfalls.

| Pitfall | Prevention Phase | Verification |
|---------|------------------|--------------|
| Two-exp non-identifiability (1) | Two-exponential fitting | Multi-start test: 5 initial guesses give same g1 within 10%; separation test passes |
| rho^{-1/4} blow-up (2) | Anti-Hermitian defect diagnosis | BohrDomain defect < 1e-10 after truncation; condition number reported |
| lambda_eff divergence (3) | Effective rate plot | No NaN in returned lambda_eff; noise floor cutoff applied; sign-change points marked |
| Richardson failure (4) | Delta-convergence diagnosis | Monotonicity test runs before Richardson; non-monotonic data blocks extrapolation |
| Memory blow-up (5) | Bootstrap error bars | n=6 with 20k trajectories and 100 bootstrap samples uses < 50 MB for stored batch data |
| Window selection failure (6) | Automatic window selection | Multiple-window comparison shows gap estimates agree within 2x bootstrap SE |
| Symmetry label ambiguity (7) | Symmetry sector analysis | Degeneracy flag triggers for eigenvalue pairs closer than 1e-8 * |gap| |
| Negative amplitude (8) | Two-exponential fitting | Cancellation ratio (|A1+A2|/(|A1|+|A2|)) > 0.3 for all converged fits |
| Bootstrap lambda_eff (9) | Bootstrap error bars | Per-batch time series stored; pointwise CI band computed |
| Testing stochastic output (10) | All diagnostic phases | Each feature has unit tests (synthetic), integration tests (n=3), and validation tests (n=4, optional) |
| Error conflation (11) | Delta-convergence diagnosis | Both simulation error (rho distance) and total estimation error (gap distance) reported separately |
| Steady-state mismatch (12) | Effective rate plot | Uses liouv_result.fixed_point when available; plots Delta(T_final) as diagnostic |

---

## Sources

### Verified from codebase (HIGH confidence)
- QuantumFurnace.jl fitting.jl: single-exponential model `A*exp(-gap*t)+C`, log-linear initial guess, LsqFit.jl `curve_fit` with parameter bounds
- QuantumFurnace.jl gap_estimation.jl: `estimate_spectral_gap` pipeline, `_select_best_observable` smallest-gap criterion, `eigenbasis_overlap_analysis` full dense eigendecomposition
- QuantumFurnace.jl trajectories.jl: `run_observable_trajectories` returns `measurements_mean` (grand average only), no per-batch storage; `_run_chunk_obs_only!` accumulates into shared `mean_data_local`
- QuantumFurnace.jl convergence.jl: `run_trajectories_convergence` has batch structure but stores only per-batch density matrix (via `_run_batch_no_obs!`), not per-batch observable time series
- Quick-30 (30-SUMMARY.md): Gap estimation error NOT O(delta); error/delta varies 96x; Richardson extrapolation 1.0x improvement (ineffective)
- Quick-32 (32-SUMMARY.md): Trajectory simulation correct; observable errors O(1e-3) and monotonic; fitting procedure is sole source of 37-49% gap estimation error
- Quick-25 (25-PLAN.md): n=6 gap mode in k=pi momentum sector; all k=0 observables have zero overlap

### Established numerical analysis (HIGH confidence)
- [Parameter identifiability in two-exponential models](https://bmcsystbiol.biomedcentral.com/articles/10.1186/s12918-015-0219-2): Ill-conditioning of sums-of-exponentials fitting depends on ratio of decay constants and signal-to-noise
- [On the accuracy of Prony's method for recovery of exponential sums with closely spaced exponents](https://www.sciencedirect.com/science/article/abs/pii/S1063520324000642): Numerical conditioning of exponential sum recovery degrades as exponent separation decreases
- [Levenberg-Marquardt algorithm](https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm): Damping regularizes ill-conditioned Jacobian; convergence depends on initialization
- [Richardson extrapolation](https://en.wikipedia.org/wiki/Richardson_extrapolation): Requires error to be a power series in the discretization parameter; fails when assumption is violated
- [Fitting sum of exponentials is ill-conditioned](https://randorithms.com/2020/03/08/exponential-sum-fits.html): On moving from two- to three-exponential models, the condition deteriorates badly
- [Modified Prony Algorithm for Exponential Function Fitting](https://epubs.siam.org/doi/abs/10.1137/0916008): SVD-based approaches improve numerical stability over direct Prony for noisy data
- Bootstrap resampling theory (Efron 1979): Block bootstrap valid for independent blocks; minimum ~50 blocks for reliable percentile CIs

### Domain knowledge (HIGH confidence)
- Similarity transform `rho^{-1/4} L rho^{1/4}` condition number grows as `exp(beta * bandwidth / 4)` -- exponential in inverse temperature and system size
- For KMS-detailed-balanced Lindbladians, L_tilde is self-adjoint (all eigenvalues real) in exact arithmetic; any non-real eigenvalues indicate either approximation error or numerical artifacts
- Effective rate `lambda_eff(t) = -d/dt log|Delta(t)|` is model-free but requires |Delta(t)| > 0 at all computed points
- Near-degenerate eigenvalue subspaces have arbitrary eigenvector orientation under standard eigendecomposition (LAPACK `dsyev`/`zheev`)
- Quick-30 empirically confirmed: single-exponential gap estimation error is dominated by fitting model mismatch (delta-independent), not Trotter discretization (delta-dependent)

---
*Pitfalls research for: v1.4 Spectral Gap Refinement -- Adding diagnostics to QuantumFurnace.jl*
*Researched: 2026-02-19*
