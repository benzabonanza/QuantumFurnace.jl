# Parameter recommendations for the implemented δ-step channel (qf-b4d.5, REVISED — quadrature recs SUPERSEDED 2026-05-06)

> **\* Update (qf-7xt, 2026-05-06)**: the **quadrature register recommendations** ($r_D$, $r_{b_-}$, $r_{b_+}$ tables in this document) are **superseded** by [`quadrature-convergence-summary.md`](quadrature-convergence-summary.md). The earlier sweeps had two confounders: (i) too-narrow $\omega$-range (truncation cutoff vs grid extent were treated as separate knobs when they should be one); (ii) `gamma_norm_factor` mismatch between test/ref grids inflated raw `‖ΔL‖` by a few %. The new sweeps use a unified $\omega$-range / truncation principle and divide by `gnf` in the difference closure to isolate pure quadrature.
>
> **The Trotter-$M$, generator-splitting, ε-regime, and δ-step recommendations in this document remain authoritative.** Only the per-register $r$-tables should be cross-checked against the new canonical reference.
>
> Numbers measured at one fixture only ($n=4$, $\beta=10$); slopes universal but $K$-prefactors and floor positions can shift by ~1–2 bits with $(n, \beta)$.

This document synthesises the findings from qf-b4d.{1, 2, 3, 4} into a master
recommendation table for the parameters that enter the simulated-algorithm
δ-step channel (TrotterDomain + GQSP, jump-sweep splitting):

```
(t0_b_minus, w0_b_minus, r_b_minus,    # outer coherent register
 t0_b_plus,  w0_b_plus,  r_b_plus,     # inner coherent register
 t0_D,       w0_D,       r_D,          # dissipative register
 M_D,        M_b_minus,  M_b_plus,     # Strang substeps per per-term Trotter cache
 δ_step)                               # outer time step
```

per (n, β, ε_target) ∈ {3, 4, 5, 6} × {5, 10, 20} × {1e-3, 1e-4, 1e-5, 1e-6},
for the **CKG KMS smooth Metropolis** filter (s = 0.25, a = 0).

## REVISIONS from the prior version (review feedback 2026-05-05)

1. **`r_b_minus` is essentially free** — `r_b_minus = 6` saturates the outer
   integral at ε ≈ 6.7e-6 regardless of β/n. The cosh-decay envelope in
   `_compute_b_minus` makes b_-(t) super-algebraic in trapezoidal sums.
   Earlier "uniform `r_b±`" recommendation over-bloated the outer register.
   Verified by `scripts/scratch_coherent_quadrature_split.jl`.

2. **Per-term Trotter caches are required** for the qf-9z0 register split.
   With a single shared `trotter.t0 = t0_D` (current code), the b_+ and b_-
   integration time-grids get rounded to multiples of t0_D, introducing an
   M-INDEPENDENT quantization error that dominates B at moderate r_b±.
   With **three separate Trotter caches** at:
   - `trotter_D.t0 = t0_D`,
   - `trotter_bm.t0 = β · t0_b_minus`,
   - `trotter_bp.t0 = β · t0_b_plus`,
   B-Trotter recovers the predicted **slope -2 in M_b±** (Strang p=2). This
   is a structural code change in `B_trotter`/`construct_lindbladian` that
   has not yet been applied; the analysis emulates it via
   `B_trotter_split` in `scripts/scratch_trotter_M_selection_v2.jl`.

3. **Kinky-Metro slope -2 in 1/N is now verified** (was "partially verified"
   in v1 because qf-b4d.2 coupled w0 with the time-window). The clean check
   in `scripts/scratch_kinky_slope_check.jl` uses **EnergyDomain** (no
   t-window) at fixed total ω-window `2·W_MAX = 10` and varies
   `(w0, r_D)` jointly to keep coverage. Observed local slopes: 2.32, 2.27,
   3.37, 2.60 in the asymptotic regime (matches `eq:per-entry-error-kinky`
   in 2_methods.tex).

4. **Smooth Metropolis saves ~6 bits in `r_D` over kinky** at ε = 1e-7. At
   `r_D = 8`, w0 = π/(5β): smooth ε ≈ 1e-7, kinky ε ≈ 1e-3. To bring kinky
   down to 1e-7 via slope -2, need ω0 ~ 1e-2-fold smaller, i.e., 6 extra
   bits in r_D. This was understated in the prior version.

5. **Trotter prefactor σ^-3 partially verified**: at fixed M=32, the
   `‖ΔL_diss‖` scaling between β=5/10/20 is empirically β^{2.5..2.7}
   instead of the predicted β^3. Discrepancy of ~0.3 in exponent
   attributed to: (a) prefactor in thesis bound being upper bound (not
   tight); (b) operator-norm vs 1→1-norm conversion factors; (c)
   sub-leading O(1/(M³σ⁴)) corrections. Direction of scaling is
   correct, only the exponent is soft.

## Inputs (revised)

### qf-b4d.1: Coherent B quadrature (`coherent-quadrature-table.md`,
`scratch_coherent_quadrature_v2.jl`, `scratch_coherent_quadrature_split.jl`; extended to
n ∈ {3, 4, 5, 6} per qf-dnb)

The b_+ trapezoidal-rule error is **slope -1 in t0' = 2T_+/2^r_b+** with empirical constant
`K = ‖ΔB‖_op · 2^(r_b±−1)` mildly n-favorable: at β=10 we measure
**K ≈ 0.036 (n=3) → 0.026 (n=4) → 0.020 (n=5) → 0.013 (n=6)** — about a 3× reduction across the
sweep. Origin: the `t = 0` L'Hôpital sample contributes `t0' · b_+(0) · K^a(0)` to the discrete
sum — a term the Cauchy P.V. integral defining `B^(s)` excludes. **The η-cutoff branch of
`_compute_b_plus_metro` is dead code under our convention `η < t0'`** (no grid sample falls inside
`(-η, η)` except `t = 0`); kinky and smooth filters are byte-identical in the discretisation. Outer
b_- is Gevrey-class super-algebraic and saturates at r_b_minus = 6 for all n in our sweep.

The recommendation table below uses the conservative `K = 0.036` (n=3) value — at larger n the
recipe is loose by up to 3×, but this only over-provisions r_b± by ~1–2 bits.

| ε_target | r_b_minus (kinky/smooth/Gaussian) | r_b_plus (kinky/smooth) | r_b_plus (Gaussian) |
|----------|-----------------------------------|--------------------------|---------------------|
| 1e-3     | 6                                 | 6                        | 6                   |
| 1e-4     | 6                                 | 7-8                      | 6                   |
| 1e-5     | 6                                 | 9-12                     | 6                   |
| 1e-6     | 6                                 | 12-15+                   | 6                   |

η is set to `ε / (β · ‖H‖ · ‖ΣA†A‖)` to satisfy the continuous-side Background-1 bias bound, but
note the empirical bias in our fixture is ~700× smaller than this bound predicts. The r_b_plus
column above is set by the `K · t0' / T_+ ≤ ε` discretisation requirement, **single-log in 1/ε,
η-independent**. The thesis eq:r_+_metro derivation is for the regime `t0' < η`, which we never
use.

Time-windows: T_- = 10 (outer), T_+ = 5 (inner).
`t0_b_minus = 2·T_- / 2^r_b_minus`, `t0_b_plus = 2·T_+ / 2^r_b_plus`.

### qf-b4d.2: Dissipative quadrature (`dissipative-quadrature-table.md`,
`scratch_dissipative_quadrature.jl`, `scratch_kinky_slope_check.jl`; extended to
n ∈ {3, 4, 5, 6} per qf-dnb, with reduced grids at n=6 for laptop wall-time)

EnergyDomain pure-ω-quadrature errors (no time-window confounding):

| filter | err at w0_D = π/(5β), 2W=10 | slope in 1/N (verified) |
|--------|-----------------------------|-------------------------|
| Gaussian | machine precision (~1e-15) | super-algebraic (slope > 6) |
| Smooth (s=0.25) | ~1e-15 once w0 < 0.04 | super-algebraic (slope ~5-30) |
| Kinky (s=0) | ~1e-3 to 1e-7 region | **algebraic, slope -2** ✓ |

Saturation at fixed `w0_D = π/(5β)` for r_D ≥ 6 (TimeDomain): the time-window
T = π/w0_D is wide enough that t-truncation/discretization is negligible —
only the ω-discretization matters. **Dominant knob: `w0_D`, not `r_D`.**

For smooth Metropolis: r_D = 6 reaches ε ≈ 1e-7 at w0_D = π/(5β). For kinky,
need 6 extra bits to match. Major savings from smoothing in s.

### qf-b4d.3: Trotter M selection — REVISED with split caches (`trotter-M-table-v2.md`,
`scratch_trotter_M_selection_v2.jl`)

With **three Trotter caches** at the natural per-register t0:

| ε_target | M_D (L_diss) | M_b_minus | M_b_plus |
|----------|-------------|-----------|----------|
| 1e-3     | 1           | 1         | 1        |
| 1e-4     | 2           | 1         | 1        |
| 1e-5     | 4           | 1         | 1        |
| 1e-6     | 16-32       | 2         | 2        |

With split caches, B-Trotter is **super-cheap** (M=1 gives B err ≈ 1e-7 at
n=3, β=10). **L_diss is the binding constraint for M.** Predicted slope -2 in
M_D, M_b_minus, M_b_plus all confirmed empirically.

If implemented as a SINGLE shared cache (current code, `trotter.t0 = t0_D`):
B saturates at err ≈ 5e-5, M-INDEPENDENT (rounding error from `round(τ·β /
trotter.t0)` mismatch). To reach ε=1e-6 in B, must increase r_D to 12+
(make t0_D match β·t0_b_plus). This is the work-around in the
unsplit code.

### qf-b4d.4: Generator splitting (`generator-splitting-table.md`)

Slope +2 confirmed perfectly (1.985–1.997). Coh-Diss splitting ~30× smaller
than Jump-wise at β=20 (because `‖L_coh‖ << ‖L_diss‖`). Per-step ε_target →
δ_step:

| ε_per_step | δ_step (jump-wise binding) |
|------------|----------------------------|
| 1e-3       | 5e-2                       |
| 1e-4       | 1.5e-2                     |
| 1e-5       | 5e-3                       |
| 1e-6       | 1.5e-3                     |

## Master table — smooth Metropolis (s = 0.25, REVISED)

Conventions:
- T_- = 10, T_+ = 5 fixed.
- w0_b_minus = π/10 = 0.314, w0_b_plus = π/5 = 0.628.
- t0_b_minus = 2·T_- / 2^r_b_minus = 20/2^r_b_minus.
- t0_b_plus = 2·T_+ / 2^r_b_plus = 10/2^r_b_plus.
- w0_D = π/(5β), t0_D = 2π/(2^r_D · w0_D) = 10β/2^r_D.
- η = ε / (β · ‖H‖ · ‖ΣA†A‖) ≈ ε / (3·β).
- **Three Trotter caches** required (assuming a future code change in
  `B_trotter`):
  - `trotter_D.t0 = t0_D`, M_D substeps.
  - `trotter_bm.t0 = β · t0_b_minus`, M_b_minus substeps.
  - `trotter_bp.t0 = β · t0_b_plus`, M_b_plus substeps.

### β = 5

| ε_target | r_b_minus | r_b_plus | r_D | M_D | M_b- | M_b+ | δ_step | η     |
|----------|-----------|----------|-----|-----|------|------|--------|-------|
| 1e-3     | 6         | 6        | 6   | 1   | 1    | 1    | 5e-2   | 7e-5  |
| 1e-4     | 6         | 7-8      | 6   | 2   | 1    | 1    | 1.5e-2 | 7e-6  |
| 1e-5     | 6         | 10       | 6   | 4   | 1    | 1    | 5e-3   | 7e-7  |
| 1e-6     | 6         | 13-14    | 6   | 16  | 2    | 2    | 1.5e-3 | 7e-8  |

### β = 10

| ε_target | r_b_minus | r_b_plus | r_D | M_D | M_b- | M_b+ | δ_step | η      |
|----------|-----------|----------|-----|-----|------|------|--------|--------|
| 1e-3     | 6         | 6        | 6   | 1   | 1    | 1    | 5e-2   | 3.3e-5 |
| 1e-4     | 6         | 8        | 6   | 2   | 1    | 1    | 1.5e-2 | 3.3e-6 |
| 1e-5     | 6         | 11       | 6   | 4   | 1    | 1    | 5e-3   | 3.3e-7 |
| 1e-6     | 6         | 14+      | 6   | 16  | 2    | 2    | 1.5e-3 | 3.3e-8 |

### β = 20

| ε_target | r_b_minus | r_b_plus | r_D | M_D | M_b- | M_b+ | δ_step | η      |
|----------|-----------|----------|-----|-----|------|------|--------|--------|
| 1e-3     | 6         | 6        | 6   | 1   | 1    | 1    | 5e-2   | 1.7e-5 |
| 1e-4     | 6         | 9        | 6   | 2   | 1    | 1    | 1.5e-2 | 1.7e-6 |
| 1e-5     | 6         | 12       | 6   | 4   | 1    | 1    | 5e-3   | 1.7e-7 |
| 1e-6     | 6         | 14+      | 6   | 16  | 2    | 2    | 1.5e-3 | 1.7e-8 |

**For the cluster (n=11, β=10)**: assumes split-cache implementation.

- r_b_minus = 6, r_b_plus = 11 → ε_quad,B ≈ 1e-5.
- r_D = 6, w0_D = π/50 = 0.063 → ε_quad,L ≈ 1e-7.
- M_D = 4 → ε_Trotter,L ≈ 1e-5.
- M_b- = M_b+ = 1 (super-cheap with split caches).
- δ_step = 5e-3 → per-step splitting ε ≈ 1e-5.
- Total ε_generator ≈ 1e-5 → trace-distance asymp ε_TD ≈ 10-100·ε_gen ≈ 1e-3 to 1e-4. Matches simulated δ/λ floor.

## Discussion

### Slope hierarchy (revised)

Per-source scaling laws and what knob refines each:

| Error source         | Slope               | Refining knob               | Verified in script              |
|----------------------|---------------------|-----------------------------|---------------------------------|
| B b_- quadrature     | super-algebraic     | r_b_minus                   | scratch_coherent_quadrature_split |
| B b_+ quadrature     | -1 in t0' (t=0 L'Hôpital sample) | r_b_plus | scratch_coherent_quadrature_v2 |
| L_diss quadrature (smooth/Gaussian) | super-algebraic | w0_D       | scratch_kinky_slope_check       |
| L_diss quadrature (kinky)           | -2 algebraic    | w0_D       | scratch_kinky_slope_check ✓     |
| L_diss Trotter (Strang)             | -2 in M_D, β^3 prefactor (observed β^{~2.6}) | M_D | scratch_trotter_M_selection_v2 |
| B b_+ Trotter (Strang, split-cache) | -2 in M_b_plus  | M_b_plus, t0_b_plus  | scratch_trotter_M_selection_v2 ✓ |
| B b_- Trotter (Strang, split-cache) | -2 in M_b_minus | M_b_minus, t0_b_minus| scratch_trotter_M_selection_v2 ✓ |
| Generator splitting (jump-wise)     | +2 in δ_step    | δ_step               | scratch_generator_splitting ✓   |
| Generator splitting (coh-diss)      | +2 in δ_step    | δ_step               | scratch_generator_splitting ✓   |

### Surprises / inconsistencies with thesis

1. **r_b_plus slope-(-1) origin is the t=0 L'Hôpital sample, NOT the η-jump**:
   the b_+ trapezoidal-rule error scales as `K · t0' = K · 2T_+/2^r_b+`,
   independent of η, with empirical constant `K ≈ 0.036`. Origin: the
   implementation evaluates `b_+(0)` via the regularised-numerator L'Hôpital
   limit (≈ 0.027), and the trap-rule sum `t0' · b_+(0) · K^a(0)` includes a
   contribution the Cauchy P.V. integral defining `B^(s)` explicitly excludes.
   Under our `η < t0'` convention the η-cutoff branch never fires, so any
   smooth-bump replacement of the indicator (a previously suggested fix, now
   scrapped — see beads `qf-oiq` and `qf-xfa`) computes byte-identical kernel
   values. Thesis eq:r_+_metro derivation is for the opposite regime
   `t0' < η`. The implementation actually achieves single-log
   `r_b+ ~ log_2(2K/ε)` ≈ `log_2(0.07/ε)`, η-independent.

2. **Single-trotter-cache rounding error**: the current
   `B_trotter`/`construct_lindbladian(TrotterDomain)` builds one Trotter
   cache shared across D / b_- / b_+, with `trotter.t0 = t0_D`. The b_-
   and b_+ time-grids round to multiples of t0_D, introducing a hidden
   M-INDEPENDENT quantization error. **Fix**: thread three Trotter caches
   through `_precompute_data` and `B_trotter` / `_precompute_coherent_unitary`
   per qf-9z0 per-term registers. Implementation change in `src/coherent.jl`
   and `src/trotter_domain.jl`.

3. **Trotter β^3 prefactor empirically β^{~2.6}**: the prediction is an
   upper bound (thesis Eq. line 1067), so this is consistent with the bound
   not being tight. No analytic claim to fix; just a quantitative slack.

4. **Kinky in time-domain L_diss**: the slope -2 in 1/w0_D was confirmed at
   FIXED ω-window via EnergyDomain. The earlier qf-b4d.2 TimeDomain w0 sweep
   muddled this with t-window effects; that's why I called it "partially
   verified" — the issue was sweeping methodology, not the prediction.

### Open code changes (recommended)

**Priority 1**: Per-term Trotter caches in TrotterDomain Lindbladian assembly
(`src/coherent.jl::_precompute_coherent_unitary`,
`src/coherent.jl::B_trotter`). Beads `qf-d0w`.

**Thesis-side action (no code change)**: 2_methods.tex Section "Discretized B"
needs a paragraph qualifying eq:r_+_metro. The slope-(-1) we observe is from
the t=0 L'Hôpital sample, not the η-jump; under `η < t0'` the η-cutoff is a
continuous-side accounting device only. Beads `qf-xfa`.

(The previously-suggested "smooth η-regularization in `_compute_b_plus_metro`"
priority is **scrapped** — see beads `qf-oiq` close-out.)

## Cross-check (one combo, from prior version) — STILL VALID

`scripts/scratch_b4d5_crosscheck.jl` runs `predict_channel_trajectory`
with the recommended parameters at (n=3, β=10, ε_gen=1e-4, smooth
Metropolis):
- Final asymptotic trace distance to Gibbs: 1.601e-3 (= 16× ε_gen).
- Channel μ_1 - 1 ≈ 1e-15.

The cross-check used the SHARED-trotter-cache (current code), so the actual
ε_generator is dominated by the B rounding error (~5e-5) rather than the
intended 1e-4 from quadrature. The 1.6e-3 asymptotic distance is consistent
with the loose `20·τ_mix·ε_eff` bound at ε_eff ≈ 5e-5. With the proposed
split-cache implementation, the cross-check would reach a tighter
asymptotic distance.
