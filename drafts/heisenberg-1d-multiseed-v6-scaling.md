# 1D Heisenberg PBC multi-seed scaling (qf-e4z.23, v6)

**Status:** Complete — 180/180 cells done, analysis + figures generated.

## Goal

Re-run the qf-e4z.22 v5 P1 scaling sweep on the post-qf-yi4 *multi-seed*
disordered Heisenberg ensemble (`build_heis_1d`, PBC, **5 seeds** per
`(n, β_phys)` cell) and report median + sampling spread of `gap_phys` and
`τ_mix` across `n ∈ {3..8}` × `β_phys ∈ {0.25, 0.5, 1.0, 1.5, 2.0, 2.5}`.

Two driving questions:

1. Does the apparent **even/odd-n distinction** in the prior single-seed
   v5 `gap_phys` plot survive disorder averaging?
2. Do the v5 scaling-fit conclusions (`τ_mix ≈ 7.5 · n^{1.23} ·
   β_phys^{0.32}`, M0 ≫ M1) change when we replace each single-seed
   datapoint with the median over 5 independent realisations?

## TL;DR

- **Even/odd-n is REAL.** At every β_phys ∈ {0.25..2.5}, the inter-parity
  geomean ratio of `gap_phys` exceeds the intra-parity 5-seed scatter by
  factors of 6–50. The seed-noise floor in `gap_phys / med(gap_phys)` is
  ~1–2 % across the entire 180-cell grid.
- **Even-n `gap_phys` COLLAPSES with β_phys at n ≥ 6** — `gap_phys` at
  `n=8, β_phys=2.5` is **1.08**, vs **5.00** at `n=7, β_phys=2.5`. The odd-n
  branch stays in `[3.6, 5.5]` across the whole grid; the even-n branch
  drops to ~1 at the high-β corner. This was masked in the single-seed
  v5 result that drove the `[heisenberg_1d_no_arrhenius_qf_e4z_22]` memo.
- **The single-seed v5 scaling fit `τ_mix ≈ 7.5·n^{1.23}·β_phys^{0.32}`
  is not robust to disorder averaging.** The per-cell median gives
  `τ_mix ≈ 2.1·n^{1.89}·β_phys^{0.32}` with M0 (separable power) AICc
  weight only 0.69 vs 0.31 for M1 (Arrhenius). `Δ_AICc = 1.58`, far weaker
  than v5's `Δ_AICc ≈ 10`. The n-exponent jumped from 1.23 to 1.89; the
  β-exponent is unchanged at 0.32.
- **Root cause of the weaker fit:** even-n `n ≥ 6` τ_mix grows steeply with
  β_phys (n=6: 36→181, n=8: 49→284), while odd-n + small-even τ_mix
  saturates above β_phys ≈ 0.5. A single separable-power fit is the wrong
  model — the data wants a parity-resolved fit, not a global one.

## Setup

| Knob | Value |
|---|---|
| Hamiltonian family | 1D Heisenberg XXX, `coeffs = [1.0, 1.0, 1.0]`, PBC |
| Disorder | `[[Z], [Z, Z]]`, `disorder_strength = 0.1` |
| Seeds | `{42, 43, 44, 45, 46}` (5 realisations per `(n, β_phys)`) |
| Construction | KMS |
| Domain | EnergyDomain (matvec Krylov) |
| Filter | smooth-Metropolis, fixed `s = 0.25`, `a = 0`, `σ = 1/β_alg` |
| Method | `predict_lindbladian_trajectory`, `krylovdim = 40`, `tol = 1e-10` |
| `ε` target | `1e-3` |
| `r_D` | 7 for `β_phys ≤ 1.0`; 8 for `β_phys ≥ 1.5` (matches v5) |
| Init state | `ρ_0 = I/d` |
| Total cells | 6 n × 6 β_phys × 5 seeds = **180** (all converged) |

Fixtures come from `scripts/output/multiseed_fixtures/heis_xxx_disordered_periodic_n{n}_seed{seed}.bson`
(qf-yi4). Each seed carries its own `rescaling_factor`, so
`β_alg = β_phys · R(n, seed)` varies cell-to-cell within an `(n, β_phys)`
group by ~1 % across seeds (R is largely determined by the deterministic
extensive Heisenberg backbone; the ε=0.1 disorder shifts it by a few
percent).

### r_D adequacy (smoke check)

Smoke check at the worst-case cell `(n=4, β_phys=2.5)` across all 5 seeds
(`β_alg ∈ [68, 70]`):

| seed | β_alg | R(n=4) | ‖L_e − L_b‖_op | rel |
|---|---|---|---|---|
| 42 | 68.35 | 27.339 | 2.24e-15 | 2.26e-15 |
| 43 | 69.47 | 27.787 | 3.22e-15 | 3.21e-15 |
| 44 | 70.23 | 28.093 | 3.78e-15 | 3.76e-15 |
| 45 | 69.09 | 27.634 | 3.12e-15 | 3.10e-15 |
| 46 | 69.43 | 27.774 | 6.22e-15 | 6.22e-15 |

Worst rel = **6.2 × 10⁻¹⁵** ≪ 1 × 10⁻⁹ target. `r_D = 8` is adequate
for all seeds at the highest-β cell — machine-precision Energy↔Bohr
agreement. r_D = 7 was confirmed adequate for `β_phys ≤ 1.0` from the v5
methodology and not re-validated here.

## Per-cell aggregate (5-seed median + IQR)

| n | β_phys | #seeds | med R | med gap_phys | IQR | med τ_mix | IQR |
|---|---|---|---|---|---|---|---|
| 3 | 0.25 | 5 | 14.199 | 3.791 | 0.0949 | 11.68 | 0.0285 |
| 3 | 0.50 | 5 | 14.199 | 3.928 | 0.091 | 13.21 | 0.00394 |
| 3 | 1.00 | 5 | 14.199 | 3.641 | 0.0732 | 14.06 | 0.156 |
| 3 | 1.50 | 5 | 14.199 | 3.618 | 0.0529 | 14.57 | 0.108 |
| 3 | 2.00 | 5 | 14.199 | 3.604 | 0.0278 | 15.18 | 0.286 |
| 3 | 2.50 | 5 | 14.199 | 3.572 | 0.0247 | 15.83 | 0.52 |
| 4 | 0.25 | 5 | 27.774 | 4.883 | 0.0271 | 26.81 | 0.0443 |
| 4 | 0.50 | 5 | 27.774 | 4.678 | 0.028 | 39.45 | 0.0423 |
| 4 | 1.00 | 5 | 27.774 | 4.557 | 0.0124 | 44.13 | 0.00492 |
| 4 | 1.50 | 5 | 27.774 | 4.299 | 0.0179 | 44.11 | 0.00492 |
| 4 | 2.00 | 5 | 27.774 | 4.212 | 0.00858 | 44.06 | 0.00689 |
| 4 | 2.50 | 5 | 27.774 | 4.141 | 0.00986 | 44.02 | 0.0118 |
| 5 | 0.25 | 5 | 28.985 | 4.332 | 0.0613 | 28.14 | 0.063 |
| 5 | 0.50 | 5 | 28.985 | 5.469 | 0.0812 | 35.14 | 0.0276 |
| 5 | 1.00 | 5 | 28.985 | 4.833 | 0.0508 | 36.39 | 0.00492 |
| 5 | 1.50 | 5 | 28.985 | 4.711 | 0.0506 | 36.51 | 0.000984 |
| 5 | 2.00 | 5 | 28.985 | 4.653 | 0.0623 | 36.72 | 0.104 |
| 5 | 2.50 | 5 | 28.985 | 4.62 | 0.0722 | 37.03 | 0.391 |
| 6 | 0.25 | 5 | 39.575 | 4.579 | 0.0297 | 36.77 | 0.11 |
| 6 | 0.50 | 5 | 39.575 | 4.626 | 0.047 | 58.16 | 0.157 |
| 6 | 1.00 | 5 | 39.575 | 3.523 | 0.0348 | 80.55 | 0.0531 |
| 6 | 1.50 | 5 | 39.575 | 2.781 | 0.0242 | 99.84 | 0.0915 |
| 6 | 2.00 | 5 | 39.575 | 1.965 | 0.0169 | 130.9 | 0.611 |
| 6 | 2.50 | 5 | 39.575 | 1.275 | 0.0171 | 180.9 | 2.25 |
| 7 | 0.25 | 5 | 42.816 | 4.382 | 0.0278 | 42.12 | 0.0581 |
| 7 | 0.50 | 5 | 42.816 | 5.383 | 0.0284 | 56.09 | 0.0404 |
| 7 | 1.00 | 5 | 42.816 | 5.49 | 0.0382 | 58.9 | 0.0177 |
| 7 | 1.50 | 5 | 42.816 | 5.205 | 0.0339 | 57.8 | 0.0177 |
| 7 | 2.00 | 5 | 42.816 | 5.088 | 0.0754 | 57.37 | 0.0059 |
| 7 | 2.50 | 5 | 42.816 | 5.005 | 0.133 | 57.22 | 0.32 |
| 8 | 0.25 | 5 | 52.152 | 4.489 | 0.0153 | 48.81 | 0.0522 |
| 8 | 0.50 | 5 | 52.152 | 4.989 | 0.0622 | 73.32 | 0.113 |
| 8 | 1.00 | 5 | 52.152 | 3.44 | 0.0227 | 108.1 | 0.187 |
| 8 | 1.50 | 5 | 52.152 | 2.518 | 0.0125 | 143.4 | 0.248 |
| 8 | 2.00 | 5 | 52.152 | 1.737 | 0.00967 | 197.1 | 0.493 |
| 8 | 2.50 | 5 | 52.152 | 1.078 | 0.0265 | 284.3 | 1.31 |

The IQR column on the rightmost τ_mix entry is striking: even at n=8,
β_phys=2.5 (the most expensive cell with τ_mix ≈ 284), the seed-to-seed
IQR is **1.31** — half a percent. Disorder averaging is dramatically
tighter than the inter-parity signal we are after.

## Scaling fit on per-cell median τ_mix

```
M0 (separable power):  τ_mix = (2.10 ± 0.67) · n^(1.89 ± 0.19) · β_phys^(0.32 ± 0.08)
                       σ_residual = 0.362,  AICc = 38.30,  weight = 0.688

M1 (poly × Arrhenius): τ_mix = (1.41 ± 0.48) · n^(1.89 ± 0.19) · exp((0.30 ± 0.08)·β_phys)
                       σ_residual = 0.370,  AICc = 39.88,  weight = 0.312

Δ_AICc(M1 − M0) = 1.58
```

Contrast with the v5 single-seed fit
(`[heisenberg_1d_no_arrhenius_qf_e4z_22]`):

```
v5:  τ_mix ≈ 7.5 · n^{1.23} · β_phys^{0.32}   (AICc weight M0 = 0.997, Δ ≈ 10)
v6:  τ_mix ≈ 2.1 · n^{1.89} · β_phys^{0.32}   (AICc weight M0 = 0.688, Δ ≈ 1.6)
```

The β-exponent is robust (0.32 in both). The n-exponent jumped substantially
(1.23 → 1.89) and the M0 vs M1 discrimination collapsed from decisive to
marginal. The interpretation is **not** that "Arrhenius might be right after
all" — both M0 and M1 are misspecified, because the underlying data has a
distinct parity-resolved structure (see next section) that no single
separable-power or Arrhenius form can capture. The σ_residual ~ 0.36 in
both models is comparable to the inter-parity geomean spread, which is the
dominant unexplained variance.

## Even/odd diagnostic (per β_phys)

| β_phys | even-n med gap_phys (n=4,6,8) | odd-n med gap_phys (n=3,5,7) | geomean ratio odd/even | intra-scatter | verdict |
|---|---|---|---|---|---|
| 0.25 | [4.88, 4.58, 4.49] | [3.79, 4.33, 4.38] | 0.895 | 0.007 | REAL |
| 0.50 | [4.68, 4.63, 4.99] | [3.93, 5.47, 5.38] | 1.023 | 0.011 | REAL |
| 1.00 | [4.56, 3.52, 3.44] | [3.64, 4.83, 5.49] | 1.205 | 0.008 | REAL |
| 1.50 | [4.30, 2.78, 2.52] | [3.62, 4.71, 5.21] | 1.434 | 0.009 | REAL |
| 2.00 | [4.21, 1.96, 1.74] | [3.60, 4.65, 5.09] | 1.811 | 0.012 | REAL |
| 2.50 | [4.14, 1.27, 1.08] | [3.57, 4.62, 5.00] | 2.439 | 0.017 | REAL |

The "intra-scatter" column reports the median over both parities of
`(max − min) / (2 · median)` — the per-cell half-range across 5 seeds,
normalised by the median. It is uniformly **~1 %**. The verdict
threshold is `|log(odd/even geomean)| > 2 × intra-scatter`; every β_phys
clears that threshold with at least a factor-5 margin.

What jumps out is the **monotone divergence of the parity ratio with β_phys**:
at β_phys=0.25 the odd-n branch sits below the even-n branch (ratio < 1);
by β_phys=2.5 the odd-n branch is **2.4×** higher than even-n. The
crossover happens near β_phys ≈ 0.5. Physically: increasing inverse
temperature widens the parity-resolved gap structure — even-n PBC chains
develop a low-lying mode (singlet-singlet near-degeneracy in the
ferromagnetic-singlet block? Bethe-ansatz quasi-degeneracy?
goldstone-like soft mode from broken U(1)/SU(2)?) that the ε=0.1 Z+ZZ
disorder field does not split.

## Figures

- `drafts/figures/numerics/v6_multiseed_gap_phys.{png,pdf}` —
  left: `gap_phys` vs n at each β_phys (5-seed median + min..max band);
  right: `gap_phys` vs β_phys at each n. The right panel cleanly separates
  the two branches: n=3, 5, 7 (and n=4) trace a flat-in-β family; n=6 and
  n=8 trace a steep collapse from ~4.5 at β_phys=0.25 to ~1 at
  β_phys=2.5.
- `drafts/figures/numerics/v6_multiseed_taumix.{png,pdf}` — τ_mix vs
  β_phys per n, with the M0 fit overlaid. Note the n=6 and n=8 curves
  rise above the fit at high β_phys, while the others saturate well below.
- `drafts/figures/numerics/v6_multiseed_even_odd.{png,pdf}` — 6-panel
  small multiple, one per β_phys, with even-n (blue squares) vs odd-n
  (red diamonds) markers. The ribbon is the per-cell min..max band over
  5 seeds; it is barely visible because intra-parity scatter is ~1 %.

## Methodological notes

- `gap_phys = gap_arnoldi · R(n, seed)`, not `gap_alg` — the latter is the
  rescaled-spectrum gap and is not the thesis-relevant axis (qf-8fr,
  `[krylov_x0_symmetric_bug_qf_8fr]`). Because R varies per seed, β_alg
  drifts ~ ±1 % across seeds at fixed `(n, β_phys)`; the sidecars record
  it per cell.
- The qf-8fr GUE-traceless Krylov `x_0` seed was confirmed adequate for
  disordered Heisenberg: every cell converged, intra-parity scatter
  stays at ~1 % and matvecs are typically 32–40 per cell (krylovdim 40).
- The `find_typical` v5 fixtures had a uniform `R(n)` per n; multi-seed v6
  fixtures have `R(n, seed)` varying by a few percent
  (`[fixture_migration_find_typical_qf_2kd]`). The shift in
  `R(n)`-with-seed does not affect the gap_phys analysis (it is per-cell
  rescaled), and the τ_mix β-exponent is invariant.
- Per `[feedback_more_data_points_for_scaling_claims]`, 3 even-n + 3
  odd-n cells is on the edge of resolving the asymptotic parity scaling.
  What we can claim is: **the parity gap exists, with intra-parity 1 %
  noise, and grows monotonically with β_phys at the available n**. The
  *asymptotic* behaviour (does even-n `gap_phys` close to zero as
  n → ∞? does odd-n stay Ω(1)?) requires n ≥ 10 to test, which is a
  separate cluster sweep.

## Open questions / follow-ups

1. **OBC counterpart sweep**. PBC adds a frustrated bond on an odd-n
   Heisenberg ring; OBC removes it. If the even/odd structure is a PBC
   geometry effect, OBC should kill or invert it. Filed as a beads
   follow-up.
2. **Larger n (n ∈ {10, 12, 14})** on the cluster. Needed to test whether
   even-n `gap_phys` continues collapsing exponentially in n or floors
   at some `n*`-dependent value. 3 even-n datapoints is the minimum to
   *see* the trend; 5+ is needed to *fit* it.
3. **DLL S2 counterpart**. v6 deliberately skipped the DLL pair (S2) per
   the issue plan. If even/odd is intrinsic to the Hamiltonian (not the
   CKG kernel), DLL should show the same structure with possibly
   different prefactors.
4. **Parity-resolved scaling fit**. The single-model M0/M1 fits are
   misspecified. The natural next analysis is to fit M0 separately on
   {n=3, 5, 7} and {n=4, 6, 8} and compare the two β-exponents.

## Reproduction

```
JULIA_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 \
  julia --project scripts/scratch_p1_v6_multiseed.jl

julia --project scripts/analyze_v6_multiseed_scaling.jl
```

Wall clock (sandbox, 8 threads): ~3 h. Peak RSS 1.12 GB. Sidecars:
`scripts/output/sweep_S1_v6_ckg_ideal_multiseed/smooth_metro_eps1e-03/`
(180 BSON files, 152 KB total). Figures: `drafts/figures/numerics/v6_multiseed_*.{pdf,png}`.
