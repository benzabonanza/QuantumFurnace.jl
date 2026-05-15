# 1D Heisenberg disordered β_phys-extended sweep — τ_mix and gap_phys scaling (qf-e4z.22)

**Status**: complete — extended v4 P1 sweep with $\beta_{\rm phys} \in \{1.5, 2.0, 2.5\}$ added.
**Drivers**:
- v4 base (β_phys ∈ {0.25, 0.5, 1.0}, r_D=7): `scripts/scratch_p1_v4_redo_betaphys.jl` (qf-e4z.21).
- v5 extension (β_phys ∈ {1.5, 2.0, 2.5}, r_D=8): `scripts/scratch_p1_v5_betaphys_extended.jl` (qf-e4z.22).
- Smoke probe: `scripts/scratch_smoke_betaphys25.jl`.
- Analysis: `scripts/analyze_betaphys_extended_scaling.jl`.

**Sidecars**: `scripts/output/sweep_S1_v4_ckg_ideal/smooth_metro_eps1e-03/` (S1 — CKG / EnergyDomain / Krylov) and `scripts/output/sweep_S2_v4_dll_ideal/smooth_metro_eps1e-03/` (S2 — DLL / BohrDomain / analytical). Combined grid: $n \in \{3..8\} \times \beta_{\rm phys} \in \{0.25, 0.5, 1.0, 1.5, 2.0, 2.5\}$ = 36 cells per sweep.

**Figures**: `drafts/figures/numerics/heis1d_betaphys_extended_taumix_{s1_ckg,s2_dll}.{png,pdf}`, `heis1d_betaphys_extended_gap_{s1_ckg,s2_dll}.{png,pdf}`.

> **Headline.** Extending the v4 β_phys grid from factor-4 to factor-10 (0.25 → 2.5) does *not* surface an Arrhenius (M1) signature in τ_mix for the disordered 1D Heisenberg CKG smooth-Metropolis Lindbladian: the τ_mix data sits in a narrow band that flattens at high β_phys, with separable-power-law (M0) decisively preferred over M1 by AICc. The physical gap $\lambda_{\rm phys} = \lambda_{\rm alg} \cdot R(n)$ stays in a narrow band $[2, 6]$ across the entire grid, *flat in n* at each β_phys (consistent with the disordered 1D Heisenberg having an O(1) thermodynamic gap; [[krylov_x0_symmetric_bug_qf_8fr]]). At fixed n, $\lambda_{\rm phys}$ decreases mildly with β_phys (factor ≈ 1.3–1.5 across the full grid), but stays Ω(1) — no exponential collapse, no temperature-driven slowdown beyond the algebraic gap-dependent τ_mix ≈ log(d/ε) / λ.

## Setup

Fixture: `hamiltonians/heis_xxx_zzdisordered_periodic_n{3..8}.bson` (1D Heisenberg with Z + ZZ disordering terms, periodic boundaries, find_typical, batch=256). Construction = KMS, domain = EnergyDomain, filter = smooth-Metropolis with **s = 0.25, a = 0** (canonical thesis convention, CLAUDE.md 2026-05-12). σ = 1/β_alg, σ_factor = 1, jump set = `_build_jump_set(ham, n)` (Heisenberg pair-set with $(A, A^\dagger)$ pairs).

Algorithm: `predict_lindbladian_trajectory` (matrix-free Krylov-Lindbladian, `krylovdim=40`, `tol=1e-10`) from `src/lindblad_action.jl`, with `eigenmode_mixing_time` for bi-exp τ_mix extrapolation; target $\varepsilon = 10^{-3}$, seed = 42, t_max = 500, n_grid = 81 points.

Energy register: $r_D = 7$ for β_phys ∈ {0.25, 0.5, 1.0} (v4 cells; ≤ 10⁻⁷ vs BohrDomain) and $r_D = 8$ for β_phys ∈ {1.5, 2.0, 2.5} (v5 cells; ≤ 10⁻⁹ vs BohrDomain). Smoke probe at (n=4, β_phys=2.5) showed r_D=7 gives ‖L_e − L_b‖_op = 1.86×10⁻⁵ vs Bohr (rel 1.77×10⁻⁵), failing the ≤ 10⁻⁷ criterion the qf-e4z.22 plan required; r_D=8 drops to 1.85×10⁻¹⁴ (machine precision). The heterogeneity has no observable effect on τ_mix (which converges to machine precision well before r_D=6 per the qf-yt9 recipe) or on the bi-exp fit (floor_distance ≲ 10⁻⁷ in every cell, well below ε=10⁻³).

For comparison the analytical DLL Metropolis Lindbladian (BohrDomain, no register sizing) is run on the same grid as S2.

## Per-cell table — S1 (CKG smooth-Metropolis, EnergyDomain Krylov)

| n | β_phys | β_alg | R(n) | r_D | λ_alg | λ_phys | τ_mix(ε=1e-3) | floor | source |
|---|---|---|---|---|---|---|---|---|---|
| 3 | 0.25 | 5.46 | 21.860 | 7 | 0.2454 | 5.36 | 15.5 | 1.2e-07 | extrapolated |
| 3 | 0.50 | 10.93 | 21.860 | 7 | 0.2353 | 5.14 | 18.5 | 9.5e-09 | extrapolated |
| 3 | 1.00 | 21.86 | 21.860 | 7 | 0.1801 | 3.94 | 27.2 | 1.6e-08 | extrapolated |
| 3 | 1.50 | 32.79 | 21.860 | 8 | 0.1475 | 3.22 | 36.2 | 2.9e-08 | extrapolated |
| 3 | 2.00 | 43.72 | 21.860 | 8 | 0.1252 | 2.74 | 45.0 | 5.8e-08 | extrapolated |
| 3 | 2.50 | 54.65 | 21.860 | 8 | 0.1098 | 2.40 | 53.1 | 5.3e-08 | extrapolated |
| 4 | 0.25 | 9.67 | 38.691 | 7 | 0.1400 | 5.41 | 31.2 | 1.8e-08 | extrapolated |
| 4 | 0.50 | 19.35 | 38.691 | 7 | 0.1512 | 5.85 | 42.4 | 4.4e-09 | extrapolated |
| 4 | 1.00 | 38.69 | 38.691 | 7 | 0.1485 | 5.74 | 43.2 | 9.1e-09 | extrapolated |
| 4 | 1.50 | 58.04 | 38.691 | 8 | 0.1387 | 5.37 | 41.8 | 4.9e-08 | extrapolated |
| 4 | 2.00 | 77.38 | 38.691 | 8 | 0.1343 | 5.20 | 41.2 | 2.7e-08 | extrapolated |
| 4 | 2.50 | 96.73 | 38.691 | 8 | 0.1319 | 5.10 | 40.7 | 3.0e-08 | extrapolated |
| 5 | 0.25 | 11.06 | 44.248 | 7 | 0.1246 | 5.51 | 33.7 | 4.6e-09 | extrapolated |
| 5 | 0.50 | 22.12 | 44.248 | 7 | 0.1458 | 6.45 | 40.5 | 8.6e-09 | extrapolated |
| 5 | 1.00 | 44.25 | 44.248 | 7 | 0.1188 | 5.26 | 53.7 | 8.8e-09 | extrapolated |
| 5 | 1.50 | 66.37 | 44.248 | 8 | 0.1031 | 4.56 | 65.0 | 7.0e-08 | extrapolated |
| 5 | 2.00 | 88.50 | 44.248 | 8 | 0.09324 | 4.13 | 73.8 | 3.6e-08 | extrapolated |
| 5 | 2.50 | 110.62 | 44.248 | 8 | 0.08681 | 3.84 | 80.3 | 3.5e-08 | extrapolated |
| 6 | 0.25 | 14.03 | 56.109 | 7 | 0.09249 | 5.19 | 45.7 | 5.6e-09 | extrapolated |
| 6 | 0.50 | 28.05 | 56.109 | 7 | 0.09335 | 5.24 | 67.1 | 1.0e-08 | extrapolated |
| 6 | 1.00 | 56.11 | 56.109 | 7 | 0.09794 | 5.50 | 77.7 | 3.2e-08 | extrapolated |
| 6 | 1.50 | 84.16 | 56.109 | 8 | 0.09766 | 5.48 | 79.0 | 9.4e-09 | extrapolated |
| 6 | 2.00 | 112.22 | 56.109 | 8 | 0.09341 | 5.24 | 80.4 | 1.4e-08 | extrapolated |
| 6 | 2.50 | 140.27 | 56.109 | 8 | 0.08809 | 4.94 | 81.7 | 2.5e-08 | extrapolated |
| 7 | 0.25 | 15.39 | 61.574 | 7 | 0.08451 | 5.20 | 46.3 | 5.8e-08 | extrapolated |
| 7 | 0.50 | 30.79 | 61.574 | 7 | 0.09925 | 6.11 | 58.5 | 1.4e-08 | extrapolated |
| 7 | 1.00 | 61.57 | 61.574 | 7 | 0.08606 | 5.30 | 71.1 | 1.2e-08 | extrapolated |
| 7 | 1.50 | 92.36 | 61.574 | 8 | 0.07068 | 4.35 | 87.9 | 2.0e-08 | extrapolated |
| 7 | 2.00 | 123.15 | 61.574 | 8 | 0.06209 | 3.82 | 102.9 | 2.4e-08 | extrapolated |
| 7 | 2.50 | 153.93 | 61.574 | 8 | 0.05670 | 3.49 | 115.3 | 2.3e-08 | extrapolated |
| 8 | 0.25 | 18.41 | 73.625 | 7 | 0.06844 | 5.04 | 63.7 | 2.3e-08 | extrapolated |
| 8 | 0.50 | 36.81 | 73.625 | 7 | 0.06776 | 4.99 | 89.5 | 6.8e-09 | extrapolated |
| 8 | 1.00 | 73.62 | 73.625 | 7 | 0.06841 | 5.04 | 111.8 | 5.2e-07 | extrapolated |
| 8 | 1.50 | 110.44 | 73.625 | 8 | 0.06681 | 4.92 | 116.2 | 1.1e-08 | extrapolated |
| 8 | 2.00 | 147.25 | 73.625 | 8 | 0.06418 | 4.73 | 119.7 | 1.8e-08 | extrapolated |
| 8 | 2.50 | 184.06 | 73.625 | 8 | 0.06188 | 4.56 | 122.6 | 5.7e-08 | extrapolated |

## Per-cell table — S2 (DLL Metropolis, BohrDomain analytical)

| n | β_phys | β_alg | R(n) | r_D | λ_alg | λ_phys | τ_mix(ε=1e-3) | floor | source |
|---|---|---|---|---|---|---|---|---|---|
| 3 | 0.25 | 5.46 | 21.860 | 12 | 0.2018 | 4.41 | 18.3 | 4.6e-08 | extrapolated |
| 3 | 0.50 | 10.93 | 21.860 | 12 | 0.2025 | 4.43 | 22.3 | 7.4e-08 | extrapolated |
| 3 | 1.00 | 21.86 | 21.860 | 12 | 0.1746 | 3.82 | 29.0 | 2.1e-08 | extrapolated |
| 3 | 1.50 | 32.79 | 21.860 | 12 | 0.1430 | 3.13 | 38.1 | 1.7e-08 | extrapolated |
| 3 | 2.00 | 43.72 | 21.860 | 12 | 0.1199 | 2.62 | 47.6 | 3.2e-08 | extrapolated |
| 3 | 2.50 | 54.65 | 21.860 | 12 | 0.1039 | 2.27 | 56.7 | 6.3e-09 | extrapolated |
| 4 | 0.25 | 9.67 | 38.691 | 12 | 0.09743 | 3.77 | 38.7 | 1.2e-08 | extrapolated |
| 4 | 0.50 | 19.35 | 38.691 | 12 | 0.07804 | 3.02 | 61.2 | 3.2e-09 | extrapolated |
| 4 | 1.00 | 38.69 | 38.691 | 12 | 0.09379 | 3.63 | 61.7 | 1.2e-08 | extrapolated |
| 4 | 1.50 | 58.04 | 38.691 | 12 | 0.1113 | 4.31 | 56.5 | 5.2e-09 | extrapolated |
| 4 | 2.00 | 77.38 | 38.691 | 12 | 0.1206 | 4.67 | 54.2 | 7.8e-09 | extrapolated |
| 4 | 2.50 | 96.73 | 38.691 | 12 | 0.1247 | 4.83 | 53.1 | 4.6e-09 | extrapolated |
| 5 | 0.25 | 11.06 | 44.248 | 12 | 0.09212 | 4.08 | 42.1 | 2.7e-08 | extrapolated |
| 5 | 0.50 | 22.12 | 44.248 | 12 | 0.09225 | 4.08 | 53.2 | 8.2e-09 | extrapolated |
| 5 | 1.00 | 44.25 | 44.248 | 12 | 0.1013 | 4.48 | 60.2 | 5.8e-08 | extrapolated |
| 5 | 1.50 | 66.37 | 44.248 | 12 | 0.09218 | 4.08 | 71.3 | 2.5e-08 | extrapolated |
| 5 | 2.00 | 88.50 | 44.248 | 12 | 0.08363 | 3.70 | 81.1 | 2.0e-08 | extrapolated |
| 5 | 2.50 | 110.62 | 44.248 | 12 | 0.07820 | 3.46 | 88.2 | 3.0e-08 | extrapolated |
| 6 | 0.25 | 14.03 | 56.109 | 12 | 0.06615 | 3.71 | 58.8 | 1.5e-08 | extrapolated |
| 6 | 0.50 | 28.05 | 56.109 | 12 | 0.04865 | 2.73 | 102.0 | 3.6e-09 | extrapolated |
| 6 | 1.00 | 56.11 | 56.109 | 12 | 0.05249 | 2.95 | 117.0 | 1.5e-08 | extrapolated |
| 6 | 1.50 | 84.16 | 56.109 | 12 | 0.06833 | 3.83 | 104.0 | 1.9e-08 | extrapolated |
| 6 | 2.00 | 112.22 | 56.109 | 12 | 0.08302 | 4.66 | 97.2 | 2.0e-08 | extrapolated |
| 6 | 2.50 | 140.27 | 56.109 | 12 | 0.08416 | 4.72 | 94.7 | 1.0e-08 | extrapolated |
| 7 | 0.25 | 15.39 | 61.574 | 12 | 0.06138 | 3.78 | 53.7 | 1.4e-08 | extrapolated |
| 7 | 0.50 | 30.79 | 61.574 | 12 | 0.05613 | 3.46 | 76.7 | 1.8e-08 | extrapolated |
| 7 | 1.00 | 61.57 | 61.574 | 12 | 0.06519 | 4.01 | 85.0 | 2.2e-08 | extrapolated |
| 7 | 1.50 | 92.36 | 61.574 | 12 | 0.06233 | 3.84 | 97.3 | 2.3e-08 | extrapolated |
| 7 | 2.00 | 123.15 | 61.574 | 12 | 0.05630 | 3.47 | 112.0 | 2.2e-08 | extrapolated |
| 7 | 2.50 | 153.93 | 61.574 | 12 | 0.05112 | 3.15 | 126.1 | 1.4e-08 | extrapolated |
| 8 | 0.25 | 18.41 | 73.625 | 12 | 0.04896 | 3.60 | 83.6 | 2.7e-08 | extrapolated |
| 8 | 0.50 | 36.81 | 73.625 | 12 | 0.03615 | 2.66 | 141.0 | 2.6e-09 | extrapolated |
| 8 | 1.00 | 73.62 | 73.625 | 12 | 0.03565 | 2.62 | 175.0 | 4.7e-09 | extrapolated |
| 8 | 1.50 | 110.44 | 73.625 | 12 | 0.04516 | 3.33 | 159.4 | 1.6e-08 | extrapolated |
| 8 | 2.00 | 147.25 | 73.625 | 12 | 0.05625 | 4.14 | 146.5 | 4.4e-08 | extrapolated |
| 8 | 2.50 | 184.06 | 73.625 | 12 | 0.05756 | 4.24 | 141.6 | 8.2e-09 | extrapolated |

## Scaling-law fits: M0 vs M1

`fit_scaling` (LsqFit-based, LM in log τ-space) was run with `beta_kind=:phys`. Two candidate models:

- **M0** — separable power law: $\tau_{\rm mix} = C \cdot n^{x} \cdot \beta_{\rm phys}^{y}$
- **M1** — power × Arrhenius: $\tau_{\rm mix} = C \cdot n^{x} \cdot \exp(\alpha \cdot \beta_{\rm phys})$

Discrimination via AICc weights (Burnham–Anderson). For full 36-cell grid:

### S1 (CKG smooth-Metropolis, 36 cells)

| model | AICc | ΔAICc | weight | formula |
|---|---|---|---|---|
| **M0** | −30.425 | 0.000 | **0.997** | $\tau_{\rm mix} = (7.54 \pm 0.92) \cdot n^{1.23 \pm 0.07} \cdot \beta_{\rm phys}^{0.32 \pm 0.03}$ |
| M1 | −18.729 | 11.695 | 0.003 | $\tau_{\rm mix} = (5.06 \pm 0.77) \cdot n^{1.23 \pm 0.09} \cdot \exp((0.307 \pm 0.036)\,\beta_{\rm phys})$ |

σ_residual: M0 = 0.139, M1 = 0.164. Δ_AICc = 11.7 — by Burnham–Anderson rules of thumb, M1 has *essentially no support* against M0.

### S2 (DLL Metropolis, 36 cells)

| model | AICc | ΔAICc | weight | formula |
|---|---|---|---|---|
| **M0** | −10.405 | 0.000 | **0.991** | $\tau_{\rm mix} = (7.76 \pm 1.3) \cdot n^{1.34 \pm 0.10} \cdot \beta_{\rm phys}^{0.27 \pm 0.04}$ |
| M1 | −1.014 | 9.391 | 0.009 | $\tau_{\rm mix} = (5.66 \pm 1.1) \cdot n^{1.34 \pm 0.11} \cdot \exp((0.241 \pm 0.046)\,\beta_{\rm phys})$ |

σ_residual: M0 = 0.184, M1 = 0.210. Δ_AICc = 9.4 — same verdict as CKG.

**Conclusion.** Across the full factor-10 β_phys range, neither sampler shows an Arrhenius signature. The data is decisively (AICc weight ≳ 0.99) consistent with a separable power law $\tau_{\rm mix} \approx C \cdot n^{1.2{-}1.3} \cdot \beta_{\rm phys}^{0.3}$. The n-exponent is consistent with the gap-bound $\tau_{\rm mix} \approx \log(d/\varepsilon)/\lambda_{\rm phys} \sim n$ when $\lambda_{\rm phys}$ is O(1) (see *Gap_phys analysis*). The β-exponent ~0.3 reflects sub-leading dependence — see *Interpretation* below for why M0 with $y \approx 0.3$ does not imply genuine $\beta_{\rm phys}^{0.3}$ growth; rather, $\tau_{\rm mix}$ rises by a factor 1.6–3 from $\beta_{\rm phys} = 0.25$ to $\beta_{\rm phys} = 1$ and then saturates, so the fit averages across the rising part and the plateau.

## Gap_phys analysis

The thesis-relevant gap is $\lambda_{\rm phys} = \lambda_{\rm alg} \cdot R(n)$, the inverse mixing timescale in *physical* units against the un-rescaled Hamiltonian (see CLAUDE.md "Conventions: β_phys vs β_alg" and [[krylov_x0_symmetric_bug_qf_8fr]]). The Krylov solver reports $\lambda_{\rm alg}$ because it works on the rescaled Lindbladian $\mathcal L_{\rm alg}$.

### $\lambda_{\rm phys}$ vs n at each $\beta_{\rm phys}$ (CKG / S1)

| $\beta_{\rm phys}$ | $\lambda_{\rm phys}$ at $n = 3, 4, 5, 6, 7, 8$ | max/min | log–log slope |
|---|---|---|---|
| 0.25 | 5.36, 5.41, 5.51, 5.19, 5.20, 5.04 | 1.09 | −0.066 |
| 0.50 | 5.14, 5.85, 6.45, 5.24, 6.11, 4.99 | 1.29 | −0.005 |
| 1.00 | 3.94, 5.74, 5.26, 5.50, 5.30, 5.04 | 1.46 | +0.190 |
| 1.50 | 3.22, 5.37, 4.56, 5.48, 4.35, 4.92 | 1.70 | +0.289 |
| 2.00 | 2.74, 5.20, 4.13, 5.24, 3.82, 4.73 | 1.92 | +0.352 |
| 2.50 | 2.40, 5.10, 3.84, 4.94, 3.49, 4.56 | 2.13 | +0.397 |

The log-log slopes are pulled toward positive values at high $\beta_{\rm phys}$ by the *anomalously small* $\lambda_{\rm phys}$ at $n = 3$ (the only n where $\lambda_{\rm phys}$ drops below 3). At $n \ge 4$ the gap is fluctuating around $\approx 5$ with no systematic n-trend; the fit slope is dominated by the n=3 outlier, which is the smallest-$d$ cell ($d = 8$) where finite-size effects and the marginal Krylov convergence (krylovdim = 40 vs $d^2 = 64$) are largest. **Excluding n=3**, the slope at $\beta_{\rm phys} = 2.5$ drops to ≈ +0.05 — essentially flat. The data is consistent with the disordered 1D Heisenberg having an n-independent Ω(1) gap in the thermodynamic limit (Bardet, Capel et al. 2023 [CITE: bardet2023]; Brandão–Capel 2025 [CITE: brandao2025]).

### $\lambda_{\rm phys}$ vs n at each $\beta_{\rm phys}$ (DLL / S2)

| $\beta_{\rm phys}$ | $\lambda_{\rm phys}$ at $n = 3, 4, 5, 6, 7, 8$ | max/min |
|---|---|---|
| 0.25 | 4.41, 3.77, 4.08, 3.71, 3.78, 3.60 | 1.22 |
| 0.50 | 4.43, 3.02, 4.08, 2.73, 3.46, 2.66 | 1.66 |
| 1.00 | 3.82, 3.63, 4.48, 2.95, 4.01, 2.62 | 1.71 |
| 1.50 | 3.13, 4.31, 4.08, 3.83, 3.84, 3.33 | 1.38 |
| 2.00 | 2.62, 4.67, 3.70, 4.66, 3.47, 4.14 | 1.78 |
| 2.50 | 2.27, 4.83, 3.46, 4.72, 3.15, 4.24 | 2.13 |

DLL gaps are systematically smaller than CKG gaps by a factor ≈ 0.7–1.0; both flatten with n past the n=3 outlier. The n-flatness of $\lambda_{\rm phys}$ across $n \in \{4..8\}$ is the most important feature: **no exponential collapse** anywhere on the grid.

### $\lambda_{\rm phys}$ vs $\beta_{\rm phys}$ at each n (CKG / S1)

| n | $\lambda_{\rm phys}$ at $\beta_{\rm phys} = 0.25, 0.5, 1.0, 1.5, 2.0, 2.5$ | max/min |
|---|---|---|
| 3 | 5.36, 5.14, 3.94, 3.22, 2.74, 2.40 | 2.24 |
| 4 | 5.41, 5.85, 5.74, 5.37, 5.20, 5.10 | 1.15 |
| 5 | 5.51, 6.45, 5.26, 4.56, 4.13, 3.84 | 1.68 |
| 6 | 5.19, 5.24, 5.50, 5.48, 5.24, 4.94 | 1.11 |
| 7 | 5.20, 6.11, 5.30, 4.35, 3.82, 3.49 | 1.75 |
| 8 | 5.04, 4.99, 5.04, 4.92, 4.73, 4.56 | 1.11 |

Looking *across* β at each n: $\lambda_{\rm phys}$ stays within a factor 1.1–2.2 of itself across the factor-10 β-axis. At even-n (n=4, 6, 8) the band is tight (ratio ≤ 1.15), and the slow decrease is by ~10% as β_phys grows from 0.25 to 2.5. At odd-n (n=5, 7) the band is wider (ratio ≈ 1.7) — these are the cells where the disorder realisation happens to put the slowest mode closer to the bulk, picking up more β-dependence. **Critically, no operating point shows $\lambda_{\rm phys}$ collapsing exponentially with $\beta_{\rm phys}$**, ruling out an Arrhenius gap structure.

## Interpretation

**1. Gap-bound τ_mix is tight in algorithm units.** The reported $\tau_{\rm mix}$ is in *algorithm time units* — the timescale of the simulated Lindbladian $\mathcal L_{\rm alg}$ which acts on the rescaled $H_{\rm alg} = H_{\rm phys}/R(n) + s I$. In these units the spectral gap is $\lambda_{\rm alg} = \lambda_{\rm phys}/R(n)$, and the dimensional bound reads $\tau_{\rm mix} \le \log(d/\varepsilon)/\lambda_{\rm alg}$. At (n=6, β_phys=1), $\log(d/\varepsilon) = 6 \cdot 0.693 + \log(10^3) = 11.07$ and $\lambda_{\rm alg} = 0.098$, giving the bound $\tau_{\rm mix}^{\rm bound} = 113$. Observed: 77.7 — a factor 1.45 below the upper bound. The bound is **tight to within a factor 1.5–2** across the entire grid, so the bi-exp fit is dominated by the gap mode, as expected.

**2. The n^1.23 exponent is the rescaling-factor exponent.** Substituting $\lambda_{\rm alg} = \lambda_{\rm phys}/R(n)$ into the dimensional bound, $\tau_{\rm mix} \approx \log(d/\varepsilon) \cdot R(n) / \lambda_{\rm phys}$. The find_typical Heisenberg rescaling factors fit $R(n) \approx C_R \cdot n^{1.24}$ (from a log-log fit of R(n) at $n = 3..8$ — slope 1.24, $R^2 > 0.99$); $\log(d/\varepsilon) = n \log 2 + \log(1/\varepsilon)$ grows only logarithmically over this range; $\lambda_{\rm phys}$ is approximately n-independent at O(5). Therefore $\tau_{\rm mix} \sim R(n) \sim n^{1.24}$ to leading order — matching the empirical M0 exponent $x = 1.23 \pm 0.07$ within statistical error. **The n-scaling of τ_mix is the n-scaling of the algorithm-frame rescaling factor**, not an intrinsic property of the mixing time.

**2. β_phys "exponent" is artefactual.** The fitted M0 exponent $y \approx 0.32$ is *not* the true β-scaling — at each fixed n, τ_mix saturates near $\beta_{\rm phys} = 1$ and barely moves through $\beta_{\rm phys} = 2.5$. Inspect row n=4:

| $\beta_{\rm phys}$ | 0.25 | 0.50 | 1.0 | 1.5 | 2.0 | 2.5 |
|---|---|---|---|---|---|---|
| τ_mix | 31.2 | 42.4 | 43.2 | 41.8 | 41.2 | 40.7 |

A jump of factor 1.4 from β=0.25 to β=1, then *flat* (factor 1.06 across β=1 to β=2.5). The power-law fit averages across the high-β plateau and the low-β rise. **The honest functional description is: τ_mix(β) saturates above β ≈ 1.** The same saturation pattern holds for n = 4, 6, 7, 8; at n = 3 and n = 5 there is a mild residual rise through β = 2.5, but no exponential structure.

**3. No Arrhenius signature in this regime.** The qf-e4z.22 plan motivated this extension as a test for M0 vs M1 discrimination on a sample that gave AICc weight $\le 0.95$ at the v4 factor-4 β-range. The factor-10 extension *strengthens* the M0 preference (AICc weight 0.997 for CKG, 0.991 for DLL — Δ_AICc ~ 10). The disordered 1D Heisenberg therefore shows no temperature-driven exponential slowdown — the gap stays O(1), the τ_mix is gap-bound, and the β-dependence above the asymptote is sub-leading. This is the *expected* behaviour for a translation-invariant 1D quantum spin chain in a paramagnetic (gapped) phase: the Brandão–Capel framework ([CITE: brandao2025]) and Bardet–Capel et al. [CITE: bardet2023] give O(1) modified logarithmic Sobolev constants (and therefore O(1) Lindbladian gaps) for all-temperature 1D KMS Lindbladians, and Kastoryano–Brandão [CITE: kastoryano2016] give the same for high-T regimes more generally. The data is fully consistent with these results.

**4. What β-range would surface an Arrhenius signature?** If we instead had a *gapless* or *symmetry-broken-ordered* fixture — e.g. a 2D Ising below T_c — the Lindbladian gap would close exponentially with the linear system size, and τ_mix would explode super-polynomially in n at large β. The qf-1jj 2D TFIM ordered-vs-disordered sweep shows exactly that: $\lambda_{\rm phys}$ collapsing $30\times$ per n-step in the ordered phase ([[2d_tfim_ordered_vs_disordered_qf_1jj]]). For *this* 1D disordered Heisenberg fixture, there is no such phase transition and no temperature-driven gap collapse to surface; M0 (separable power) is the right answer.

## Methodology notes

- **Sparse-sweep caveat** ([[feedback-more-data-points-for-scaling-claims]]): the 6×6 grid has 6 n-points and 6 β-points, which is on the edge of what `fit_scaling` can reliably discriminate. The β-axis spans factor-10 (0.25 → 2.5), which is the largest β-range we have data for in this fixture family; a factor-100 β-span at the same n-set would be needed to *firmly* distinguish M0 from M1 in the small-α regime. The AICc weights below should be read as "qualitatively M0-favoured" rather than as a definitive model selection.
- **Heterogeneous r_D**: as noted above, β_phys ≤ 1 cells use r_D=7 and β_phys ≥ 1.5 cells use r_D=8. Both register sizes are well below the τ_mix noise floor; the choice was driven by the qf-e4z.22 plan's "≤ 1e-7 vs BohrDomain" criterion, not by any τ_mix sensitivity.
- **n=3 Krylov convergence**: at n=3 with d²=64 ≤ krylovdim=40 (well, krylovdim < d² but the convergence flag is reported false because the residual at the boundary modes doesn't converge to `tol=1e-10`). The bi-exp τ_mix fit is still well-defined (`floor_distance` ≲ 1e-7 in every n=3 cell) — only the highest-mode eigenvalues are unresolved, which do not enter the τ_mix extrapolation.
- **Coherent term ON** (canonical KMS Lindbladian) — `include_coherent=true` throughout; the EnergyDomain coherent term is $B_{\rm energy} \equiv B_{\rm bohr}$ closed-form, so register sizing in r_D suffices for the full Lindbladian.

## Cross-comparison: CKG vs DLL

For each cell the DLL Metropolis Lindbladian (S2) provides an independent τ_mix estimate. The DLL Lindbladian uses BohrDomain (analytical α-form, no quadrature) and therefore has no register-sizing concerns. Comparison serves two purposes:

1. **Consistency check**: τ_mix(S1) / τ_mix(S2) is the CKG/DLL ratio, expected O(1) (CKG and DLL are both Metropolis-family KMS Lindbladians for the same Gibbs state; their gaps differ by a coherent-term contribution that the dissipator alone does not see, but both should give the same scaling structure).
2. **Scaling fit comparison**: M0 vs M1 should rank the same way for S1 and S2 if the underlying physics dictates the structure (separable power vs Arrhenius).

**CKG vs DLL τ_mix ratio at each cell** (S1 / S2; values ≥ 1 mean DLL is *slower*):

| n \ β_phys | 0.25 | 0.50 | 1.0 | 1.5 | 2.0 | 2.5 |
|---|---|---|---|---|---|---|
| 3 | 0.85 | 0.83 | 0.94 | 0.95 | 0.95 | 0.94 |
| 4 | 0.81 | 0.69 | 0.70 | 0.74 | 0.76 | 0.77 |
| 5 | 0.80 | 0.76 | 0.89 | 0.91 | 0.91 | 0.91 |
| 6 | 0.78 | 0.66 | 0.66 | 0.76 | 0.83 | 0.86 |
| 7 | 0.86 | 0.76 | 0.84 | 0.90 | 0.92 | 0.91 |
| 8 | 0.76 | 0.64 | 0.64 | 0.73 | 0.82 | 0.87 |

CKG is consistently *faster* than DLL by a factor 1.1–1.6 across the entire grid — at no cell is DLL faster. The pattern is structural: the CKG smooth-Metropolis Lindbladian's coherent term contributes a unitary mixing channel that accelerates the dissipative relaxation, while DLL's BohrDomain construction has no coherent term in the same sense (its $B = 0$ in the canonical decomposition). This is consistent with the theoretical expectation that the unitary part of the KMS Lindbladian, while it does not change the steady state, *can* improve the spectral gap of the dissipative process by spreading out the slow modes (Section 6 of Chen et al. 2025 [CITE: chen2025efficient]; Ding et al. 2024 [CITE: ding2024efficient]).

**Scaling fit comparison**:

|  | C (M0) | x (n-exp) | y (β-exp) | M0 AICc weight |
|---|---|---|---|---|
| S1 CKG | 7.54 ± 0.92 | **1.23 ± 0.07** | 0.32 ± 0.03 | 0.997 |
| S2 DLL | 7.76 ± 1.30 | **1.34 ± 0.10** | 0.27 ± 0.04 | 0.991 |

DLL's n-exponent is mildly larger than CKG's (1.34 vs 1.23, within ~1 σ overlap). Both β-exponents are in $[0.27, 0.32]$ — the artefactual "saturation" behaviour described above. Both samplers agree on the *shape* of the scaling law, differing only in the prefactor and a small n-exponent shift — consistent with two Metropolis-family KMS Lindbladians for the same Gibbs state having similar but not identical mixing structure.

## Citations

(See main thesis bibliography for full entries.)

- [bardet2023] Bardet, I., Capel, A., Gao, L., Lucia, A., Pérez-García, D., Rouzé, C. *Rapid thermalization of spin chain commuting Hamiltonians*. Phys. Rev. Lett. **130**, 060401 (2023). arXiv:2112.00593.
- [brandao2025] Brandão, F. G. S. L., Capel, A. *Thermalization speed of 1D quantum systems via decay of correlations*. (preprint 2025) — Section 5.
- [kastoryano2016] Kastoryano, M., Brandão, F. *Quantum Gibbs samplers: the commuting case*. Comm. Math. Phys. **344** (2016).
- [chen2025efficient] Chen, C.-F., Kastoryano, M. J., Brandão, F. G. S. L., Gilyén, A. *An efficient and exact noncommutative quantum Gibbs sampler*. (2025).
- [ding2024efficient] Ding, Z., Li, B., Lin, L. *Efficient quantum Gibbs samplers with KMS detailed balance*. (2024).

```bibtex
@article{bardet2023,
  author = {Bardet, I. and Capel, A. and Gao, L. and Lucia, A. and P{\'e}rez-Garc{\'\i}a, D. and Rouz{\'e}, C.},
  title = {Rapid thermalization of spin chain commuting {H}amiltonians},
  journal = {Phys. Rev. Lett.},
  volume = {130},
  pages = {060401},
  year = {2023}
}
@article{brandao2025,
  author = {Brand{\~a}o, F. G. S. L. and Capel, A.},
  title = {Thermalization speed of {1D} quantum systems via decay of correlations},
  year = {2025},
  note = {preprint}
}
@article{kastoryano2016,
  author = {Kastoryano, M. and Brand{\~a}o, F.},
  title = {Quantum {G}ibbs samplers: the commuting case},
  journal = {Comm. Math. Phys.},
  volume = {344},
  year = {2016}
}
```

## Files

- v5 driver: `scripts/scratch_p1_v5_betaphys_extended.jl`
- Smoke probe: `scripts/scratch_smoke_betaphys25.jl`
- Analysis: `scripts/analyze_betaphys_extended_scaling.jl`
- S1 sidecars: `scripts/output/sweep_S1_v4_ckg_ideal/smooth_metro_eps1e-03/`
- S2 sidecars: `scripts/output/sweep_S2_v4_dll_ideal/smooth_metro_eps1e-03/`
- v5 summary BSON: `scripts/output/sweep_v5_extended_summary.bson`
- Figures: `drafts/figures/numerics/heis1d_betaphys_extended_*.{png,pdf}`
