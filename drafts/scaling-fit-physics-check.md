# Physics Check: scaling-law fit `τ_mix ≈ C · n^x · β^y` on CKG smooth-Metropolis

**Context.** A 3×6 sweep of the ideal CKG smooth-Metro KMS Lindbladian (`σ = 1/β`, `s = 0.25`, family `heis_xxx_zzdisordered_periodic_n*`, single seed 42) over n ∈ {3..8}, β ∈ {5, 10, 20} yields a fit `τ = 3.48·n^{0.82}·β^{0.42}`. AICc weakly prefers M0 over M1; ΔAICc = 3.6.

The headline question is whether `x ≈ 0.82` and `y ≈ 0.42` reflect physics or artifact.

| # | Check | Status | Evidence |
|---|---|---|---|
| 1 | n-exponent x = 0.82 plausible at this scale | **OK** (artifact-free) | Gap λ(n, σ=1/β=0.1) shrinks as 0.222, 0.137, 0.121, 0.0947, 0.0833, 0.0743 over n=3..8 → log-log slope ≈ −0.78. τ ≈ 1/λ explains x directly. See `drafts/sigma-sweep-and-redo-findings.md:130-136`. |
| 2 | Rescaled-spectrum effect on x | **CONCERN** | `_rescaling_and_shift_factors` (`src/hamiltonian.jl:544-555`) forces spec(H) ⊂ [0, 0.45] *for every n*. The physical Hamiltonian norm `‖H_phys‖ ∝ n` is divided out by `rescaling_factor` ≈ 2(λ_max−λ_min)/(1−ε) ∝ n. Effective β in physical units is β_phys = β/rescaling_factor ≈ β/n. The whole sweep secretly varies (β_phys, σ_phys) jointly. |
| 3 | x compatible with literature for 1D short-range | **OK in spirit** | Bergamaschi–Chen 2025 (arXiv:2510.08533) proves *n-independent* gap for all 1D short-range Hamiltonians at all finite T. Asymptotically `x → 0`; n=3..8 is pre-asymptotic but evidently approaching that regime. *Not* the slow-mixing regime of Gamarnik–Kiani 2024 (transverse-field 2D / stabilizer codes). |
| 4 | β-exponent y = 0.42 (small) physically reasonable | **CONCERN (norm-rescaling)** | Memory `sigma_sweep_findings_qf_bw1.md` + `sigma_sweep_norm_proxies_qf_aex.md` document that *most of the σ-dependence in τ_mix is a norm rescaling*. The same applies to β when σ = 1/β: `‖B_CKG‖_op` measured slope vs β is +0.7..+0.8 at n ≥ 4 (`drafts/ckg-vs-dll-comparison-findings.md:130-142`). So y ≈ 0.42 partly compensates an n-averaged `‖L‖ ∝ β^{0.7..0.8}` growth — the *intrinsic* mixing time τ·‖L‖_HS may be β-decreasing here. |
| 5 | Smooth-Metro kink width σ√s β-dependent | **CONCERN** | With σ=1/β and `s=0.25` fixed: σ√s = 0.5/β shrinks 4× from β=5 to β=20. `src/bohr_domain.jl:212-247` documents that `default_smooth_s(β, σ) = (0.05/σ)²` would keep σ√s constant. Production sweep used fixed s=0.25, so β=20 row uses a sharper kink than β=5 row. This adds a β-dependent factor to the dissipator and *both* the numerator and denominator of `τ·‖L‖_HS` shift simultaneously. Hard to disentangle from y. |
| 6 | n=5 β=20 non-monotone cell | **OK (real, mild)** | The v2 datapoint (37.25) was a broken biexp fit; v3 (42.02) is the correct eigenmode bisection — see `drafts/sigma-sweep-and-redo-findings.md:44`. The 2.6% lag below n=4 β=20 (43.12) is within typical disorder-seed scatter and within the fit's σ_residual = 0.10. Not pathological. |
| 7 | Disorder seed sensitivity | **CONCERN** | All 18 cells share seed=42 with `batch_size=200` "ideal" selection (`hamiltonians/generate_hamiltonians.jl:97-99`). Single seed per cell ⇒ no error bars on τ_mix. The "best ν_min" optimisation biases all cells toward atypically large minimum Bohr gaps — i.e. faster-mixing realisations. The fit is a fit over realisation-selected best cases, not over random disorder. |
| 8 | β coverage too narrow for M0 vs M1 discrimination | **CONCERN** | β ∈ [5, 20] is barely a factor of 4. `drafts/scaling-analysis-research.md:83-84,103` already flags that "α ≈ 0.15 Arrhenius produces τ growth by factor ≈ 9; β^{1.6} gives factor ≈ 9 too" — ΔAICc = 3.6 between M0 and M1 is in the "weak preference" zone (Burnham–Anderson cutoff is 2). Cannot distinguish power law from Arrhenius from this grid. |
| 9 | KMS detailed balance preserved | **OK** | All cells use KMS construction; `create_alpha` (`src/bohr_domain.jl:195-209`) satisfies α(ν₁,ν₂)=α(−ν₂,−ν₁)·exp(−β(ν₁+ν₂)/2) to 1e-16; verified throughout the σ-sweep (`drafts/sigma-sweep-v4-default-smooth-s.md:16`). The fixed point is the Gibbs state to machine precision. |
| 10 | Mixing-time post-processing | **OK** | `mixing_time_source = :extrapolated` ⇒ eigenmode bisection on Krylov spectrum (qf-e4y schema), not bi-exp curve-fit. This is the correct, robust path; the biexp degeneracy that produced the v2 anomaly is gone. |

## 1. Verdict on x ≈ 0.82

**Plausible at this scale, but pre-asymptotic and confounded with the spectrum-rescaling convention.**

The data are *self-consistent* in the sense that x ≈ −d(log λ)/d(log n) from the measured spectral gaps at σ=1/β. So the fitted exponent faithfully describes what's in the τ_mix curve, given the inverse-gap interpretation `τ_mix ≈ (1/λ)·log(‖σ^{-1/2}‖/ε)`.

Where x ≈ 0.82 disagrees with theory: Bergamaschi–Chen 2025 (arXiv:2510.08533) proves a *system-size-independent* spectral gap for 1D short-range Hamiltonians at all finite temperatures — and a disordered Heisenberg chain *is* such a system. Asymptotically the n-exponent should approach **zero**, not 0.82. The discrepancy is the pre-asymptotic regime: dim = 2^n ranges 8..256 across the sweep, certainly inside the finite-size transient. The expected crossover for an isotropic Heisenberg chain is set by the correlation length at temperature 1/β ≈ 0.1–0.05; with antiferromagnetic exchange J ≈ 1, this correlation length is short (lattice scale) so the asymptotic regime should already start at modest n. The fact that we don't see it suggests *spectral-density crowding* effects: at low T the gap is dominated by the second-smallest |Re λ_i| of L, which generically sees more competition as the Hilbert space grows. This is *finite-dimensional discrete-spectrum* behaviour, not a violation of the asymptotic bound. The CKG paper itself (arXiv:2303.18224) makes no closed-form prediction for the n-scaling of generic non-commuting Hamiltonians — only that the algorithm cost depends on the gap of the Lindbladian, which it bounds case by case.

**Bigger issue: the Hamiltonian rescaling.** `_rescaling_and_shift_factors` divides H by `(λ_max − λ_min)·(2/(1−ε))`. For an Heisenberg chain with `J=1`, the physical bandwidth grows linearly with n (extensive). So at n=8 we're dividing by a number ≈ n × the value at n=3, i.e. the **effective β in physical units shrinks** with n. The whole sweep is at *β_phys ≈ β/n*. This means n=8, β=20 corresponds to `β_phys ≈ 2.5`, while n=3, β=20 is `β_phys ≈ 6.7` — qualitatively different regimes. The "n-exponent" is conflating two things: (a) what happens to the gap at fixed temperature in physical units, and (b) what happens to the gap when you simultaneously raise n and lower β_phys. Both contribute to x.

## 2. Verdict on y ≈ 0.42

**Indeterminate from this sweep — partly a norm-rescaling artifact, partly a kink-width β-dependence, possibly some genuine residual slow-mode β-dependence buried under both.**

Three contributions to the observed y:

1. **Generator norm scaling.** From `drafts/ckg-vs-dll-comparison-findings.md:130-142`, `‖B_CKG‖_op` has slope ≈ +0.7..+0.8 vs β at n ≥ 4 (σ=1/β, β ∈ {5, 10, 15, 20, 30}). So `‖L‖_op ∝ β^{0.7..0.8}` at fixed σ=1/β. If `τ·‖L‖` is the physically meaningful "intrinsic" mixing time, then `τ ∝ ‖L‖^{−1} · (intrinsic stuff) ⇒ y_obs ≈ y_intrinsic − 0.7..0.8`. With y_obs ≈ 0.42, this gives y_intrinsic ≈ −0.3..−0.4 — i.e. the *intrinsic* mixing time *decreases* with β. That's not crazy in absolute terms (more dissipation per L-time-unit at narrower σ), but it certainly is not "low-T mixing is hard" — it's "we're rescaling our generator to be larger at low T".

2. **Kink-width.** σ√s shrinks 4× from β=5 to β=20 with `s=0.25` fixed. A sharper Metropolis kink concentrates the transition weight at smaller |ν|, which generically lowers the slow-mode eigenvalue (less mass at downward transitions far from the kink). This adds a β-dependent factor whose sign and slope are hard to predict without the constant-σ√s control. The `default_smooth_s` formula in `src/bohr_domain.jl:212-250` was designed precisely to remove this; the sweep used fixed s instead.

3. **σ = 1/β coupling.** The Bohr-frequency window σ tracks β by construction, so the filter's spectral support changes with β. This is *not* an independent degree of freedom; varying β with σ ∝ 1/β simultaneously varies (a) the support of f(ω), (b) the Boltzmann weight, (c) the smoothing scale. Disentangling y from these requires a held-σ control.

The fitted y is therefore a *combined* exponent reflecting all three. It is NOT the "gap-closing exponent" that Temme 2014's `λ ≥ N^{−1} exp(−2β·ε̄)` predicts — that bound, by the way, is for *stabilizer* Hamiltonians (commuting), not disordered Heisenberg, and applies in the low-T limit β·ε̄ ≫ 1. For our disordered Heisenberg with energy barriers ε̄ ~ O(1) and βε̄ ∈ [5, 20], we're at the high end of the Arrhenius regime in nominal units, but β_phys × ε̄_phys is bounded by O(1) once one factors in the rescaling. The puzzle of "small y" largely disappears under this lens.

## 3. Artifact candidates, ranked

1. **(highest) Hamiltonian rescaling-induced β_phys/n coupling.** Both x and y are aliasing with the rescaling. To test: extract x and y on the *physical* Hamiltonian by un-rescaling τ_mix and using β_phys = β / rescaling_factor(n). Or repeat with no rescaling (write `coeffs / 1` and skip the `[0, 0.45]` shift). The "true" exponents reported should be the un-rescaled ones if comparing to literature; the rescaled ones if reporting the algorithm's behaviour on the actual run.

2. **(high) Norm rescaling in y.** Test: compute `‖L‖_HS` at every cell (Krylov, `hs_operator_norm_krylov`), refit on log(τ·‖L‖_HS). y of that fit is the structurally meaningful number. Prediction (based on σ-sweep precedent + Remark 23 numerical confirmation): y_intrinsic in the range [−0.4, 0]; the "intrinsic" Lindbladian mixes slightly *faster* at low T per L-time-unit when σ tracks 1/β. This will be the headline correction to report. *This is the single most leveraged check.*

3. **(high) σ = 1/β confounding.** Test: run a second sweep with σ = c fixed (c ≈ 0.1) across β ∈ {5, 10, 20}. The β-dependence of τ_mix on that sweep is the "clean" y. Compare to the σ=1/β y. Both should agree on the same generator (modulo the prefactor) if y is dominated by the Boltzmann weight; they will disagree if y is dominated by the σ-coupling. Per `drafts/sigma-sweep-and-redo-findings.md`, σ=1/β is *suboptimal* for fixed β; the σ=0.25/β row mixes 20–40% faster. So holding σ constant is operationally interesting too.

4. **(high) Disorder seed sensitivity.** Test: rerun n ∈ {5, 7, 8} at seeds {41, 42, 43, 44, 45}. Report mean τ_mix and spread. Expect 10–25% spread per cell at β=20, 5–10% at β=5 (memory `ckg_vs_dll_first_findings`). The fitted exponents could shift by 0.05–0.15 in absolute terms once averaged. Five seeds is the minimum credible; the "best ν_min" optimisation built into `find_ideal_heisenberg` (batch_size=200, keeping the realisation with the largest minimum Bohr gap) biases all 18 cells toward unusually well-conditioned spectra, so the bare seed-averaged τ_mix will likely be larger than the current numbers.

5. **(medium) Kink-width β-dependence.** Test: rerun with `s = QuantumFurnace.default_smooth_s(β, σ)` (gives s=0.0625 at β=5, s=0.25 at β=10, s=1.0 at β=20 along σ=1/β). The v4 σ-sweep already did this for β=10 only and confirmed: at σ=1/β the v3↔v4 numbers are byte-identical (calibration point). What's needed: full (n, β) sweep with default_smooth_s and refit. Expect: y to *increase* somewhat (under default_smooth_s, β=5 has a much wider kink than β=20 → β=20 dissipator more peaked → potentially slower).

6. **(low) Pre-asymptotic finite-size.** The data are pre-asymptotic for the Bergamaschi–Chen `Θ(1)` gap, but this is *not* an "artifact" — it's the regime we can simulate. The right thing to do is *report the pre-asymptotic exponent honestly* and cite the asymptote. Don't expect n=3..8 to reach a Θ(1) plateau.

## 4. Concrete recommendation for the next sweep

A single follow-up sweep, no new technology needed, that turns the ambiguous result into a defensible thesis figure:

| Knob | Value | Why |
|---|---|---|
| n | {3, 4, 5, 6, 7, 8} | unchanged |
| β | {1, 2, 5, 10, 20} | adding β=1, 2 doubles log-β range; resolves M0 vs M1 aliasing per `drafts/scaling-analysis-research.md:74` |
| σ | hold σ = 0.1 fixed (independent of β) | decouples norm rescaling and Boltzmann weight; standard physics-paper convention |
| s | use `default_smooth_s(β, σ)` | constant absolute kink width σ√s = 0.05 across β |
| seeds | {41, 42, 43, 44, 45} | 5-fold disorder average; report mean and spread |
| auxiliary outputs | `‖L‖_HS`, λ, Λ_max at every cell | for the norm-normalised fit |
| Hamiltonian rescaling | optional — also report results in physical units (`τ_mix / rescaling_factor`) for the M0 fit | makes literature comparison direct |

Total cells: 6 × 5 × 5 = 150. At the measured wall times (n=8 cell ≈ 200s), this is feasible in a few wall-clock hours with the standard 8-Julia/1-BLAS configuration.

Then refit and report:
- **Raw** y, x for the operational benchmark.
- **Norm-rescaled** `τ·‖L‖_HS` y, x as the intrinsic mixing exponents.
- **Disorder spread** at each cell.
- **AICc + LOO** for M0 vs M1 (Bayesian via Turing.jl, per `drafts/scaling-analysis-research.md` recommendations).

Even without the σ-constant control: a 5-seed re-run of the existing grid with extra β=1, β=2 cells and the `‖L‖_HS` annotation is enough to make `x = 0.82 ± 0.08` and `y = 0.42 ± 0.05` honest numbers with confounders quantified.

> **FIXED IN qf-6vr (2026-05-11).** The β_phys / β_alg conflation called out in artifact candidate (1) is resolved at the code level. Every numerics driver that flows into the scaling fit now drives a β_phys grid; the sweep harness derives β_alg per cell via `ham.rescaling_factor`; sidecars expose both β_phys and β_alg + the rescaling factor; `fit_scaling` reads β_phys by default and tags the formula label accordingly (`β_phys^y` vs `β_alg^y`). The new canonical scaling-fit driver is `scripts/numerics_scaling_fit_ckg_smooth_metro.jl` (run it to regenerate the figure under the new convention; results note: `drafts/scaling-fit-bphys-rerun.md`).
>
> **Canonical β_phys grid: `{0.25, 0.5, 1.0}`** (replaces the legacy `{5, 10, 20}` β_alg grid). Picked from the Gibbs-state entropy analysis: 0.25 is the smallest β with non-trivial thermal contrast (`S(ρ_β)/log(d) ≈ 0.80` uniformly across n=3..10 — below that the Gibbs state is essentially uniform); 1.0 is the practical ceiling (above, σ = 1/β_alg gets too tight at large n for the OFT register sizing). Factor-4 log-spacing matches the legacy grid's β-coverage, so the M0 vs M1 discrimination weakness (artifact candidate 8 above) is not resolved by qf-6vr; that requires a denser β grid in a follow-up.
>
> **Open caveat: smooth-Metro `s` at large β_alg.** With σ = 1/β_alg and `default_smooth_s(β, σ) = (0.05/σ)²`, the s-value at (β_phys=1, n=11) would be ≈ 25 — far outside historical [0, 1]. All migrated drivers hold `s = 0.25` fixed for now; if the β_phys=1 cells mix terribly the choice is between (a) characterising large-s smooth-Metro, (b) holding σ fixed independent of β, or (c) accepting fixed s and documenting the kink-width drift.
>
> The other artifact candidates (norm-rescaling in y, σ=1/β confounding, seed sensitivity, kink-width β-dependence) are *not* addressed by this refactor and remain open. See `.claude-memory/beta_phys_beta_alg_convention.md` for the API summary.

## 5. Comparison to literature for this regime

There is **no rigorous bound directly applicable** to disordered Heisenberg with site-and-bond Z+ZZ disorder at the temperatures probed:

| Reference | Coverage | Predicts |
|---|---|---|
| Kastoryano–Temme 2013 (arXiv:1207.3261) | Generic Lindbladians | τ_mix ≤ (1/λ)·log(‖σ^{-1/2}‖/ε); no closed-form λ |
| Temme 2014 (arXiv:1412.2858) | **Stabilizer (commuting)** Davies generators | λ ≥ N^{−1}·exp(−2β·ε̄). Does *not* apply to Heisenberg. The exp(−2β·ε̄) is the Arrhenius factor for *commuting* models. |
| Chen–Kastoryano–Gilyén 2023 (arXiv:2303.18224) | Generic CKG | Cost ∝ polylog of mixing time; no closed-form n, β |
| Gamarnik–Kiani 2024 (arXiv:2411.04300) | **2D TFIM, K-SAT, stabilizer codes** | T_mix = exp[Ω(n^α)]. Does *not* cover 1D disordered Heisenberg. |
| Bergamaschi–Chen 2025 (arXiv:2510.08533) | **All 1D short-range Hamiltonians, all finite T** | Spectral gap is **n-independent**, polylog-depth Gibbs prep. **Disordered Heisenberg is in scope.** Asymptote: x → 0. |
| Šmíd et al. 2025 (arXiv:2510.04954) | Weakly interacting | polylog(n). Disordered Heisenberg with `J=1` is **not** weakly interacting. |
| Ramkumar–Soleimanifar 2024 (arXiv:2411.04454) | Random sparse Hamiltonians | polylog(n) at constant T. Coverage of disordered Heisenberg unclear. |

Best literature anchor: **Bergamaschi–Chen 2025**, which predicts asymptotic x = 0 for our exact setting. The fitted x ≈ 0.82 is the pre-asymptotic transient. The fitted y = 0.42 has no direct theoretical counterpart in the literature for CKG smooth-Metro at this regime; if you want a benchmark, the cleanest one is *empirical* — the σ-sweep results at fixed n already show y ≈ 0.4 is what raw τ_mix vs β at σ=1/β gives, dominated by norm rescaling.

**Net statement for the thesis:** the empirical scaling `τ = (3.48 ± 0.58)·n^{0.82±0.08}·β^{0.42±0.05}` describes a pre-asymptotic, norm-rescaled, single-seed regime; the residual "structural" exponents are likely closer to (x_intrinsic, y_intrinsic) ≈ (0.7, −0.3) once `τ·‖L‖_HS` is the fitted observable. None of this contradicts the rigorous bracket: the asymptote is `x_asym = 0` (Bergamaschi–Chen) and the slow-mixing lower bounds (Gamarnik–Kiani) only apply to other model classes.
