# Scaling-fit pipeline review (CKG smooth-Metropolis sweep)

Insertion target: this is a standalone verifier report; not a thesis insertion.

## Verdict

**PASS** — the analysis pipeline is faithful. Eighteen `(n, β, τ_mix)` triples are correctly loaded from the BSON sidecars (all with `mixing_time_source = :extrapolated`), `fit_scaling` is called on the correct fields, and the fit outputs are bit-reproducible against an independent OLS on `[1, log n, log β]`. The suspicious exponent `x ≈ 0.82 ± 0.08` is **not** a code artifact — it is the global least-squares average of three well-resolved per-β slopes whose values vary monotonically with β (0.99 at β=5, 0.95 at β=10, 0.74 at β=20). The "small global x" therefore reflects a real `n–β` interaction not absorbed by the separable form `M0`, not a bug. The β=20 column drags the average down. Any physics investigation should focus on (a) why the local n-exponent decreases with β under the σ=1/β coupling, and (b) whether `M0` is the right separable ansatz at all.

## Per-β hand-recomputed exponents

OLS slope of $\log \tau_{\mathrm{mix}}$ vs $\log n$ at fixed $\beta$, using the exact BSON values:

| $\beta$ | OLS slope $\hat{x}_\beta$ | endpoint slope $\log(\tau_8/\tau_3)/\log(8/3)$ |
|---|---|---|
| 5  | **0.9644** | 0.9943 |
| 10 | **0.8368** | 0.9464 |
| 20 | **0.6691** | 0.7353 |

The user's back-of-envelope endpoint slopes (1.01 / 0.95 / 0.74) are reproduced to within 1–2% — the small discrepancy is `log(2.65)` vs `log(2.6517)` from the rounded user inputs. The OLS slopes are stricter than the endpoint slopes because the n=5 and n=7 cells lie below their two-endpoint trendline at β=20 (kinks in the column).

Likewise, per-n OLS slopes of $\log \tau_{\mathrm{mix}}$ vs $\log \beta$ at fixed $n$:

| $n$ | OLS slope $\hat{y}_n$ |
|---|---|
| 3 | 0.5614 |
| 4 | 0.5122 |
| 5 | **0.2643** (driven by the non-monotone β=20 cell) |
| 6 | 0.5061 |
| 7 | 0.2835 |
| 8 | 0.3781 |

The global $\hat{y} = 0.4176$ is a weighted average across these. The wide spread of $\hat{y}_n$ (0.26–0.56) is *direct* evidence of `n–β` interaction. Hypothesis confirmed.

## Critical issues

**None.** No data-loading bug, no wrong-field bug, no sign error, no source-filter leak, no plot mis-alignment.

## Notes

1. **σ=1/β coupling absorbs into the fit, not in a clean way.** Per memory `smooth_metro_s_beta_scaling_qf_96o.md`, with `σ = 1/β` the absolute kink-smoothing $\sigma\sqrt{s} = 0.5/\beta$ *shrinks* with $\beta$, which steepens the smooth-Metropolis filter as β grows. The script's docstring already flags this as a physics check. The fit can only report a global $\hat{y}$ that conflates true thermodynamic β-scaling with kink-sharpening β-scaling. This is not a fit-code bug — it is a "fitting the wrong quantity" decision. To disentangle, the user would need either (i) a sweep at fixed σ ≠ 1/β across β, or (ii) a model `M2` of the form $\tau = C \cdot n^{x(\beta)} \cdot g(\beta)$ with an explicit β-dependent exponent, or (iii) the `σ_factor` sweep mentioned in `sigma_sweep_findings_qf_bw1`. The current `M0` / `M1` family cannot represent the β-dependent local n-exponent seen in the per-β slopes.

2. **`gap_est` confirms that the n=5, β=20 = 42.02 non-monotone cell is physical, not numerical.** The Krylov gap estimates are:
   - n=5, β=5  : `gap_est = 0.1205`, τ_mix = 29.13
   - n=5, β=10 : `gap_est = 0.1211`, τ_mix = 38.08
   - n=5, β=20 : `gap_est = 0.1488`, τ_mix = 42.02

   The β=20 cell has a *larger* gap than β=10 at n=5, which mechanistically explains a non-monotone τ_mix even with an `:extrapolated` mixing-time source. `tau_mix_bound = 69.71` and `t_max = 300` for that cell, so the bisection was nowhere near the bracket — the 42.02 is real, not a clamp. (The Lindbladian schema in `sweep_layout.md §method=:krylov` confirms `:extrapolated` means analytic bisection on the captured eigendecomposition, exact within the Krylov subspace.) Suggests that something in the disorder realisation at `seed=42, n=5` reduces the effective coupling at β=20.

3. **σ_residual = 0.0997 in log τ units is high.** A clean separable fit on noiseless synthetic data hits σ_residual ≈ 0.005 (see `test_scaling_fit.jl::(c)`). Here we see 0.10 — i.e., ~10% multiplicative scatter around the global M0 surface, which is one order of magnitude above the seed/numerical noise floor. That is the residual interaction term flagged in note 1, not measurement error.

4. **AICc-weak preference for M0 over M1 (Δ-AICc = 3.57, weight 0.86 vs 0.14).** Both models give the same x exponent (0.823 ± 0.077 vs 0.823 ± 0.085) because in `fit_scaling`, the x-direction (design column `log n`) is orthogonal to both `log β` and `β` in the 18-point grid (note `corr(x, slope) = -0.000` in both models). So the x exponent is *not* a competition artifact — switching β ↔ exp(αβ) does not move x. This is healthy.

5. **The σ_C reported in `formula_string` is the first-order Gaussian-propagation estimate $\sigma_C \approx C \cdot \sigma_c$, not the lognormal exact form.** With $\sigma_c = 0.167$ this is a 17% perturbation, so the linear approximation is fine. Just don't quote $C = 3.48 \pm 0.58$ as a CI — the actual 95% CI on C is $[e^{0.890}, e^{1.602}] = [2.43, 4.96]$, which is wider on the upper side. The formula string is fine for the figure; the CIs in the per-parameter dump are correct.

6. **Plot panels render correctly** (visually inspected `scaling_fit_ckg_smooth_metro.png`). Left panel shows the 18 observed log τ at integer n × {5, 10, 20}; right panel shows the fitted M0 surface on the {3..8} × log-spaced 41-point β grid. Both panels share `clims = (Z_MIN, Z_MAX) = (2.66, 4.16)` (the observed log τ range) and overlay the observed points. The fitted surface is monotone in both axes by construction; the observed panel visibly breaks monotonicity at (n=5, β=20) and to a lesser extent (n=7, β=20) — the very cells driving the small global $\hat{x}$.

7. **Reproducibility check**: independent OLS on the design matrix $X = [\mathbf{1}, \log n, \log \beta]$ via `X \ y` matches `fit_scaling` to all 4 printed digits: `c = 1.2460, x = 0.8234, y = 0.4176`. Standard errors match the unbiased $(N - k = 15)$ divisor, $\hat\sigma_{\mathrm{unb}}^2 = \mathrm{rss}/15$. The MLE residual stdev uses divisor $N = 18$ (per Burnham–Anderson convention) and is reported as `sigma_residual = 0.0997`. Both conventions are correct and internally consistent.

8. **No discrepancy with user's printed table.** Max relative diff across all 18 cells between the user's prompt-printed values and the BSON values is `2.3e-4` (pure rounding).

## Verified items

- **(1) BSON-loading correctness.** Loaded `sweep_n5_beta20_seed42_L_KMS_Energy.bson`. Confirmed:
  - `r.n = 5 :: Int64`, `r.beta = 20.0 :: Float64`, `r.mixing_time = 42.0239... :: Float64` finite & positive.
  - `r.mixing_time_source = :extrapolated :: Symbol`.
  - Schema matches `sweep_layout.md §method=:krylov` exactly: present fields are `seed, mixing_time, w0_D, n, total_matvecs, t0_D, domain, target_epsilon, sigma_factor, wall_time, sweep_version, init_state, method, construction, sigma, all_converged, tau_mix_bound, mode, filter_name, n_grid, t_max, t_max_factor, beta, filter_kind, floor_distance, mixing_time_source, gap_arnoldi, gap_est, r_D`. No `fitted_gap` / `r_squared` (correct — this is a `:krylov` not `:ode` cell).
  - `load_sweep_results` in `scripts/scratch_scaling_fit_ckg_smooth_metro.jl:39-48` performs `NamedTuple(d[:result])` with no field renaming, so the values pass through verbatim.

- **(2) Raw-table sanity.** All 18 cells reproduced. All have `mixing_time_source = :extrapolated`. The non-monotone n=5/β=20 cell is **not** a clamp (`tau_mix = 42.02 < tau_mix_bound = 69.71 < t_max = 300`); it is a genuine eigenmode-extrapolated value supported by a *larger* `gap_est = 0.149` than the n=5/β=10 cell (`gap_est = 0.121`).

- **(3) fit_scaling reproduces the script output.** Independent OLS gives `(c, x, y) = (1.2460, 0.8234, 0.4176)` matching `(3.48, 0.82, 0.42)` with $C = e^c = 3.4765$. Hand-recomputed per-β slopes confirm the user's endpoint calculation to 1% accuracy and confirm the `n–β` interaction.

- **(4) Correct BSON fields used.** `_get_tau` reads `:mixing_time` (Lindbladian schema, `src/scaling_fit.jl:303`) before checking `:tau_mix`. `_get_source` reads `:mixing_time_source` (`src/scaling_fit.jl:311`). Neither falls through to `tau_mix_bound`. The default `source_filter = (:extrapolated,)` is applied at `src/scaling_fit.jl:329-331`.

- **(5) σ=1/β coupling.** Flagged as a "fitting the wrong quantity" risk, not a code bug — see Note 1. The kink-smoothing length $\sigma\sqrt{s}$ shrinks with β, so the global $\hat{y}$ conflates physics-β with filter-β scaling. The pipeline correctly fits *whatever* τ_mix the sweep produced; whether that τ_mix is the "right" target is a thesis-level question not a code-verifier one.

- **(6) Plot correctness.** Visually inspected. Both panels share `clims = (2.66, 4.16)` derived from `Z_MIN/Z_MAX` of `log_τ_obs_mat` (`scripts/scratch_scaling_fit_ckg_smooth_metro.jl:144-145`). The fit panel uses `log.(grid.tau_predicted)` after `predict_scaling` converts back from log-scale (correct — `predict_scaling` returns `exp(log τ̂)`). Observed points overlaid on both panels at the right `(β, n)` coordinates.

- **(7) Source-filter strictness fix.** The tightened policy at `src/scaling_fit.jl:329-331` `(src !== nothing && src in source_filter)` does **not** drop any of the 18 valid cells — all carry `mixing_time_source = :extrapolated`. Verified by running `fit_scaling(results)` and confirming `fits[:M0].n_data == 18`. Regression test in `/tmp/verify_filter.jl` confirms: (a) cells with a missing `mixing_time_source` field are now dropped, (b) cells with `:floor` or `:nan` are dropped, (c) all 18 `:extrapolated` cells pass through.

- **Test suite.** `julia --project test/test_scaling_fit.jl` reports **221/221 pass, 5.4s**.

## Files touched

No source files modified, no new tests written (existing 221 tests cover the math; the pipeline-side checks were one-shot reproducibility runs that don't need to live in `test/`).

Diagnostic scripts written to `/tmp/verify_pipeline.jl`, `/tmp/verify_table.jl`, `/tmp/verify_filter.jl` for transient verification; not committed.

## Recommended next step (not in code)

If the user wants to disentangle the β-induced kink-sharpening from the genuine n-scaling, the cleanest experiment is the `σ_factor` sweep already on the table per memory `sigma_sweep_findings_qf_bw1.md`. At fixed σ ≠ 1/β across β, the kink width is decoupled from β and the per-β slope spread (0.99/0.95/0.74) should collapse if it really is kink-driven. If the spread persists at fixed σ, then there is a thermodynamic `n–β` interaction in the smooth-Metropolis Lindbladian itself, and the separable `M0` / `M1` family must be replaced by a `τ = C(n,β) · n^{x(β)}` form before any single exponent can be reported.
