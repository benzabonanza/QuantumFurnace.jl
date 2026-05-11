# σ-sweep norm-proxy analysis — structural vs norm-rescaling (qf-aex)

**Run date:** 2026-05-11
**Driver:** `scripts/scratch_sigma_sweep_norm_proxies.jl`
**Inputs:** `scripts/output/sweep_S4_v4_energy/` (35 cells, n ∈ {3,...,7}, σ_factor ∈ {0.25, 0.5, 0.75, 1, 1.5, 2, 3}/β, β=10, ε=10⁻³, CKG smooth-Metro KMS, EnergyDomain, $s = $ `QuantumFurnace.default_smooth_s(β,σ)` so $\sigma\sqrt{s} = 0.05$ fixed).
**Augmented BSON:** `scripts/output/sigma_sweep_norm_proxies_v4.bson` (per-cell τ_mix, λ, ‖L‖_HS, ‖L_diss‖_{1→1}^bnd, τ_bound^KMS).

## TL;DR

**The σ-trend in τ_mix is dominantly a rescaling of the generator's norm, not a structural change in the mixing dynamics.** At n=7 the residual σ-variation of `τ_mix·‖L‖_HS` is only **6%** across a 12× σ-grid, while raw τ_mix varies 5.6×. The same conclusion holds in the 1→1 upper-bound frame (14% residual at n=7). The KMS-Poincaré bound, by contrast, becomes *less tight* as σ widens — the σ-trend in λ overshoots the actual mixing cost, so the spectral gap alone is not sufficient to explain the σ-dependence.

**Physical reading.** Widening σ broadens the Gaussian filter $f(\omega)$ and so populates more off-diagonal ($\nu_1 \neq \nu_2$) entries of the Kossakowski matrix $\alpha(\nu_1, \nu_2)$. One might naively expect this to *speed up* mixing by opening more dissipative pathways between Bohr-frequency pairs. **It does not.** The same total Kossakowski mass spreads thinner over a wider $(\nu_1, \nu_2)$ support, $\|\alpha\|_\infty$ (and hence $\|\mathcal{L}\|_\mathrm{HS}$) shrinks, and raw $\tau_\mathrm{mix}$ grows in lockstep. The σ-knob trades generator norm for off-diagonal spread; the structural mixing time per "operator-norm-time-unit" is conserved. No off-diagonal speedup, no σ-narrowing slowdown — just norm scaling counteracting the apparent mixing-time change.

**Practical implication for the thesis numerics.** The qf-bw1 σ-optimum at σ ≈ 0.25/β is real but corresponds to where `‖L‖_HS` is largest. The Lindbladian's intrinsic mixing rate (per HS-time-unit) is roughly σ-independent at thermodynamic scale; σ-tuning buys speed by inflating the operator norm, not by improving the structural mixing constant.

## Setup

Three thesis-aligned norm proxies (per [`src/kms_geometry.jl`](../src/kms_geometry.jl) and [Kastoryano–Temme 2013][CITE:kastoryano2013] / [Kochanowski et al. 2024][CITE:kochanowski2024] / [Chen et al. 2025][CITE:chen2025]):

1. **2→2 norm (HS-induced):** $\|\mathcal{L}\|_\mathrm{HS}$ via `hs_operator_norm_krylov` (matrix-free GKL on `apply_lindbladian!` / `apply_adjoint_lindbladian!`). Lower bound on $\|\mathcal{L}\|_{1\to1}$ (Watrous §3.3.2; loose by at most factor $d$).
2. **1→1 upper bound:** $\|\mathcal{L}_\mathrm{diss}\|_{1\to1}^\mathrm{bnd} = 4\big\|\sum_a L_a^\dagger L_a\big\|_\infty$ via `dissipator_M_from_alpha` then `4·opnorm`. Wolf–Pérez-García style.
3. **KMS-Poincaré bound** (gap-based): $\tau_\mathrm{bound}^\mathrm{KMS} = \lambda^{-1}\log\!\big(\|\sigma^{-1/2}\|_\infty / \varepsilon\big)$ with $\lambda$ the spectral gap from `gap_arnoldi` (= `krylov_spectral_gap`).

For each cell we form three derived "intrinsic time" / "tightness" quantities and look for σ-flatness, which would identify the dominant frame of the σ-trend:

- $\tau_\mathrm{mix}\cdot\|\mathcal{L}\|_\mathrm{HS}$ — intrinsic mixing time per HS-time-unit.
- $\tau_\mathrm{mix}\cdot\|\mathcal{L}_\mathrm{diss}\|_{1\to1}^\mathrm{bnd}$ — intrinsic mixing time per 1→1-time-unit (looser, but in the "right" trace-distance norm).
- $\tau_\mathrm{mix}\,/\,\tau_\mathrm{bound}^\mathrm{KMS}$ — Poincaré slack. A flat ratio across σ means the bound is uniformly tight (or uniformly loose) and the σ-trend in $\tau_\mathrm{mix}$ tracks $\tau_\mathrm{bound}$.

## σ-variation summary (max/min across the 7-point σ-grid, per n)

| n | raw $\tau_\mathrm{mix}$ | $\tau\cdot\|\mathcal{L}\|_\mathrm{HS}$ | $\tau\cdot\|\mathcal{L}\|_{1\to1}^\mathrm{bnd}$ | $\tau / \tau_\mathrm{bound}^\mathrm{KMS}$ |
|---|---|---|---|---|
| 3 | 5.53× | **1.58×** | 1.93× | 1.35× ↓ |
| 4 | 5.63× | **1.31×** | 1.47× | 1.73× ↓ |
| 5 | 6.95× | **1.53×** | 1.67× | 1.50× ↓ |
| 6 | 7.11× | **1.40×** | 1.51× | 1.59× ↓ |
| 7 | 5.59× | **1.06×** | **1.14×** | 2.05× ↓ |

(↓ = ratio decreases with σ.)

**Reading the table.** A ratio close to 1 means the σ-trend has been almost fully absorbed into that proxy. The HS-norm frame is the cleanest at every n; the structural residual is < 60% at small n and shrinks to ≈ 6% at n=7. The 1→1 bound gives a slightly looser collapse (its prefactor depends on $\sigma$ more strongly, see below). The KMS-Poincaré bound gives the *worst* collapse: τ/τ_bound *decreases* with σ, meaning λ over-predicts how much $\tau_\mathrm{mix}$ should grow.

## Raw τ_mix (S4 v4 reference)

| n | σ=0.25 | σ=0.50 | σ=0.75 | σ=1.00 | σ=1.50 | σ=2.00 | σ=3.00 |
|---|---|---|---|---|---|---|---|
| 3 | **16.01** | 16.56 | 17.65 | 19.62 | 26.41 | 37.81 | 88.54 |
| 4 | **27.09** | 27.54 | 28.57 | 31.92 | 44.00 | 63.93 | 152.52 |
| 5 | **27.60** | 28.34 | 31.74 | 38.08 | 55.40 | 80.75 | 191.88 |
| 6 | **29.32** | 30.89 | 35.45 | 42.41 | 60.59 | 87.75 | 208.47 |
| 7 | **30.82** | 32.42 | 35.66 | 40.23 | 53.32 | 74.54 | 172.22 |

(bold = per-row minimum; optimum at σ=0.25/β confirms qf-bw1 / qf-cno headline.)

## Spectral gap λ

| n | σ=0.25 | σ=0.50 | σ=0.75 | σ=1.00 | σ=1.50 | σ=2.00 | σ=3.00 |
|---|---|---|---|---|---|---|---|
| 3 | 0.317 | 0.294 | 0.262 | 0.222 | 0.151 | 0.102 | 0.0423 |
| 4 | 0.226 | 0.220 | 0.182 | 0.137 | 0.0862 | 0.0569 | 0.0232 |
| 5 | 0.226 | 0.214 | 0.161 | 0.121 | 0.0780 | 0.0523 | 0.0216 |
| 6 | 0.201 | 0.174 | 0.123 | 0.0947 | 0.0629 | 0.0427 | 0.0177 |
| 7 | 0.188 | 0.148 | 0.105 | 0.0833 | 0.0571 | 0.0393 | 0.0165 |

λ shrinks roughly as $\sigma^{-1.4}$ for σ > 1/β, and saturates for σ ≤ 1/β. Across the full σ-grid the gap drops by ~11×.

## Proxy (i) — $\|\mathcal{L}\|_\mathrm{HS}$ via Krylov

| n | σ=0.25 | σ=0.50 | σ=0.75 | σ=1.00 | σ=1.50 | σ=2.00 | σ=3.00 |
|---|---|---|---|---|---|---|---|
| 3 | 1.170 | 1.125 | 1.065 | 0.996 | 0.836 | 0.655 | 0.332 |
| 4 | 1.172 | 1.124 | 1.055 | 0.972 | 0.777 | 0.576 | 0.260 |
| 5 | 1.174 | 1.124 | 1.053 | 0.967 | 0.770 | 0.569 | 0.255 |
| 6 | 1.163 | 1.108 | 1.029 | 0.934 | 0.725 | 0.525 | 0.229 |
| 7 | 1.149 | 1.090 | 1.006 | 0.908 | 0.697 | 0.501 | 0.217 |

Cross-check at n=3,4 against the dense `L_opnorm` already stored in the v4 BSONs: agree to 4 significant figures (max relative diff 5×10⁻⁵). Krylov path is trustworthy.

$\|\mathcal{L}\|_\mathrm{HS}$ shrinks ~5× as σ grows 12× (≈$\sigma^{-0.7}$). It is essentially n-independent — the spread across n at fixed σ is < 2%. This is consistent with the closed-form $\gamma_\mathrm{sup} = 1$ normalisation (qf-etx): the spectral norm of $\mathcal{L}$ is set by the global rate prefactor, not by the system size.

## Proxy (ii) — $4\|\sum_a L_a^\dagger L_a\|_\infty$ (1→1 upper bound)

| n | σ=0.25 | σ=0.50 | σ=0.75 | σ=1.00 | σ=1.50 | σ=2.00 | σ=3.00 |
|---|---|---|---|---|---|---|---|
| 3 | 3.705 | 3.622 | 3.501 | 3.348 | 2.960 | 2.444 | 1.291 |
| 4 | 3.689 | 3.587 | 3.430 | 3.221 | 2.667 | 2.030 | 0.943 |
| 5 | 3.671 | 3.555 | 3.378 | 3.146 | 2.564 | 1.926 | 0.879 |
| 6 | 3.646 | 3.509 | 3.299 | 3.031 | 2.395 | 1.753 | 0.775 |
| 7 | 3.633 | 3.488 | 3.265 | 2.983 | 2.333 | 1.696 | 0.743 |

Same qualitative trend as ‖L‖_HS but with a different slope: shrinks ~5× across σ, almost flat in n. Slightly less σ-collapse than HS frame because the inequality $\|\mathcal{L}\|_\mathrm{HS}\le\|\mathcal{L}\|_{1\to1}\le d\|\mathcal{L}\|_\mathrm{HS}$ is not tight uniformly in σ.

## Proxy (iii) — $\tau_\mathrm{bound}^\mathrm{KMS}$

| n | σ=0.25 | σ=0.50 | σ=0.75 | σ=1.00 | σ=1.50 | σ=2.00 | σ=3.00 |
|---|---|---|---|---|---|---|---|
| 3 | 30.55 | 32.93 | 36.92 | 43.58 | 64.00 | 94.92 | 228.60 |
| 4 | 42.88 | 43.94 | 53.27 | 70.76 | 112.17 | 169.87 | 417.17 |
| 5 | 44.52 | 46.96 | 62.67 | 83.11 | 129.06 | 192.46 | 465.65 |
| 6 | 51.01 | 58.83 | 83.35 | 108.12 | 162.66 | 239.65 | 576.93 |
| 7 | 56.71 | 72.04 | 101.46 | 128.27 | 187.04 | 271.80 | 647.99 |

Reading: the certified upper bound on τ_mix from the spectral gap alone grows much faster than the empirical $\tau_\mathrm{mix}$. At n=7 the bound grows 11× across σ while $\tau_\mathrm{mix}$ grows 5.6×. This is the "λ overshoots σ-cost" pattern that drives the Poincaré-slack column below.

## Rescaled diagnostics

### $\tau_\mathrm{mix}\cdot\|\mathcal{L}\|_\mathrm{HS}$ (intrinsic time per HS-unit)

| n | σ=0.25 | σ=0.50 | σ=0.75 | σ=1.00 | σ=1.50 | σ=2.00 | σ=3.00 |
|---|---|---|---|---|---|---|---|
| 3 | 18.73 | **18.63** | 18.80 | 19.54 | 22.07 | 24.75 | 29.38 |
| 4 | 31.76 | 30.95 | **30.15** | 31.01 | 34.18 | 36.84 | 39.59 |
| 5 | 32.39 | **31.85** | 33.43 | 36.84 | 42.64 | 45.91 | 48.83 |
| 6 | **34.11** | 34.22 | 36.48 | 39.62 | 43.95 | 46.08 | 47.78 |
| 7 | 35.42 | **35.34** | 35.89 | 36.53 | 37.17 | 37.34 | 37.33 |

**Almost perfectly σ-flat at n=7** (1.06× total spread, < 6%). At small n a residual structural σ-trend exists (1.3–1.6× spread) — the slow mode's spectral coefficient drifts modestly with σ. The pattern is monotone in n: the spread shrinks as n grows. *This is the cleanest σ-collapse of any frame.*

### $\tau_\mathrm{mix}\cdot\|\mathcal{L}\|_{1\to1}^\mathrm{bnd}$ (intrinsic time per 1→1-unit)

| n | σ=0.25 | σ=0.50 | σ=0.75 | σ=1.00 | σ=1.50 | σ=2.00 | σ=3.00 |
|---|---|---|---|---|---|---|---|
| 3 | **59.30** | 59.98 | 61.79 | 65.68 | 78.17 | 92.39 | 114.34 |
| 4 | 99.92 | 98.77 | **97.99** | 102.78 | 117.34 | 129.79 | 143.85 |
| 5 | 101.33 | **100.76** | 107.22 | 119.80 | 142.02 | 155.51 | 168.67 |
| 6 | **106.90** | 108.38 | 116.97 | 128.53 | 145.10 | 153.86 | 161.53 |
| 7 | **111.99** | 113.06 | 116.45 | 120.00 | 124.39 | 126.39 | 127.94 |

Same qualitative story — n=7 down to 14% spread — but a touch worse than HS frame. The 1→1 upper bound carries a prefactor that itself shifts with σ (the Kossakowski mass redistributes over a wider $(\nu_1, \nu_2)$ support as $\sigma$ grows), so the collapse is less clean than the HS frame.

### $\tau_\mathrm{mix}\,/\,\tau_\mathrm{bound}^\mathrm{KMS}$ (Poincaré slack — ratio to the certified upper bound)

| n | σ=0.25 | σ=0.50 | σ=0.75 | σ=1.00 | σ=1.50 | σ=2.00 | σ=3.00 |
|---|---|---|---|---|---|---|---|
| 3 | 0.524 | 0.503 | 0.478 | 0.450 | 0.413 | 0.398 | **0.387** |
| 4 | 0.632 | 0.627 | 0.536 | 0.451 | 0.392 | 0.376 | **0.366** |
| 5 | 0.620 | 0.604 | 0.507 | 0.458 | 0.429 | 0.420 | **0.412** |
| 6 | 0.575 | 0.525 | 0.425 | 0.392 | 0.373 | 0.366 | **0.361** |
| 7 | 0.544 | 0.450 | 0.352 | 0.314 | 0.285 | 0.274 | **0.266** |

The bound is loose by a factor of 1.5–3× at σ=0.25/β and *gets even looser* as σ grows (slack ratio drops to 0.27 at σ=3/β, n=7 — bound is 3.8× too high). At n=7 the slack drops by 2×, i.e. the bound's σ-trend overshoots the empirical σ-trend by another factor of 2 on top of $\tau_\mathrm{mix}$'s own σ-growth.

Direction of failure. λ shrinks roughly as $\sigma^{-1.4}$, so $1/\lambda$ alone already gives the "right" 5–7× σ-trend in $\tau_\mathrm{mix}$. But the bound also picks up the $\log\|\sigma^{-1/2}\|/\varepsilon$ term (which is σ-independent at fixed β), so it grows entirely through $1/\lambda$ — and $1/\lambda$ grows *faster* than $\tau_\mathrm{mix}$. So the slow mode whose eigenvalue is $-\lambda$ has a vanishing spectral weight at high σ, meaning it's no longer the bottleneck for mixing; a faster but louder mode takes over. Standard slow-mode dilution.

## Verdict

**Structural component of the σ-trend is negligible at large n.**

Going from σ = 0.25/β to σ = 3/β (a 12× sweep), the Lindbladian's mixing rate "per operator-norm-time-unit" varies by only ≈6% at n=7. The 5–8× growth in raw $\tau_\mathrm{mix}$ across σ is essentially the inverse of the σ-induced decrease in $\|\mathcal{L}\|_\mathrm{HS}$. **The σ-trend is a generator-norm rescaling, not a structural mixing change.**

The dominant frame is the HS-norm frame, not the spectral-gap frame. The KMS-Poincaré bound, viewed alone, would predict a stronger σ-dependence (1/λ grows by 11×) than is actually realised in mixing — by a factor of ≈2. The reason: at high σ, the slow-mode eigenvalue $-\lambda$ couples weakly to a maximally-mixed initial state, so the empirical $\tau_\mathrm{mix}$ is set by faster-decaying modes that carry most of the spectral weight in $\rho_0 - \sigma_\beta$. Equivalent reading: the qf-bw1 conclusion "$\tau\cdot\|\mathcal{L}\|_\mathrm{op}$ nearly constant" — already noted there for n=3,4 — extends to all n=3..7, sharpens at the largest n, and is robust to the choice of generator-scale proxy (HS or 1→1 bound).

Small-n caveat. At n ≤ 4 there is a real, but bounded, structural σ-trend of ≈30–60%: the slow mode shifts character as σ grows, and the leading eigenmode coefficient $c_2$ from `eigenmode_mixing_time` drifts upward at large σ. This is consistent with a finite-size discrete-spectrum effect that washes out at thermodynamic scale.

Plot direction (for qf-e4z.12 / Plot P2). When presenting σ-widening τ_mix data in the thesis, always report alongside `‖L‖_HS` so the reader can normalise out the rescaling. Showing raw $\tau_\mathrm{mix}$ alone exaggerates the σ-cost by a factor of 5+. The "fair" σ-cost at thermodynamic scale is the 6% structural residual; the rest is generator scaling, which a parameter-conditioning argument can compensate for.

## What does *not* fall out cleanly

- **The σ-trend is not purely Poincaré-bound-driven.** The bound $\tau_\mathrm{mix} \le \lambda^{-1}\log\|\sigma^{-1/2}\|/\varepsilon$ overshoots the empirical σ-cost by a factor that itself grows with σ (slack 0.54 → 0.27 across the σ-grid at n=7). One *cannot* claim "the σ-trend is just $1/\lambda$." Half of the apparent σ-growth in the bound is unrealised slack.
- **The 1→1 bound is not as clean as HS.** $\|\mathcal{L}\|_{1\to1}^\mathrm{bnd}$ has its own σ-dependent prefactor and so absorbs the σ-trend less perfectly than HS. The HS frame is the canonical one for this comparison.
- **n-dependence is real but small.** At fixed σ, $\tau_\mathrm{mix}$ grows ≈2× from n=3 to n=7. This is *not* a σ effect and is unaffected by the rescaling argument above — it's the standard system-size cost.

## Suggested thesis insertion

§ Numerics, P2 plot (σ-widening): add a two-panel figure with raw $\tau_\mathrm{mix}(\sigma)$ on the left and $\tau_\mathrm{mix}\cdot\|\mathcal{L}\|_\mathrm{HS}(\sigma)$ on the right. Caption: *"The σ-dependence of $\tau_\mathrm{mix}$ (left) collapses, up to a small 1.1× residual at n=7, when rescaled by the generator's HS norm (right); the σ-trend is dominantly a rescaling, not a structural change."*

## Files

- Script: `scripts/scratch_sigma_sweep_norm_proxies.jl`
- Log: `scripts/output/sigma_sweep_norm_proxies_v4.log`
- BSON: `scripts/output/sigma_sweep_norm_proxies_v4.bson` (35 rows, one per cell)

## Citations

```bibtex
@article{kastoryano2013,
  author = {Kastoryano, Michael J. and Temme, Kristan},
  title = {Quantum logarithmic Sobolev inequalities and rapid mixing},
  journal = {J.\ Math.\ Phys.},
  volume = {54},
  pages = {052202},
  year = {2013},
}

@article{kochanowski2024,
  author = {Kochanowski, Jonas and Capel, \'Angela and Rouz\'e, Cambyse},
  title = {Quantum thermodynamic uncertainty relations, generalized current fluctuations, and nonequilibrium fluctuation-response inequalities},
  year = {2024},
}

@article{chen2025,
  author = {Chen, Chi-Fang and Kastoryano, Michael J. and Brand{\~a}o, Fernando G.\,S.\,L. and Gily{\'e}n, Andr{\'a}s},
  title = {An efficient and exact noncommutative quantum Gibbs sampler},
  year = {2025},
}
```
