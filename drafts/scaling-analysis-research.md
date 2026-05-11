# Empirical scaling laws from a 2D (n, β) sweep: a research memo

Audience: PhD numerics chapter for KMS detailed-balance Lindbladians (CKG-style "smooth Metropolis"; Ding–Li–Lin DLL variant; Chen 2025-era successors).
Quantity: real-time mixing time τ_mix as a function of system size n ∈ {3,…,11} (9 points) and inverse temperature β ∈ {5, 10, 15, 20} (4 points). Each grid point is a single Lanczos–Krylov estimate with relative accuracy ~1e-5 to 1e-6 but mild heteroscedasticity.

## Bottom line

For 36 noisy points, do **not** start from `τ_mix ≈ C · n^x · β^y` and ask "what are the best exponents?". Treat the analysis as a model-comparison problem in log τ-space, fit a small set of candidate forms with **Bayesian weighted nonlinear regression** (Turing.jl in Julia, or PyMC / emcee in Python), and compare them with **AICc** (small-sample bias correction) plus **PSIS-LOO cross-validation** as an independent check. Use Houdayer–Hartmann data collapse as a *visual* corroboration, not as the primary inference tool — pure FSS collapse routines (pyfssa, BSA) are critical-phenomena tools that assume a critical point and don't quite fit our setting.

The reason for this stance is physical: every existing rigorous bound on τ_mix for KMS Gibbs samplers in the *generic, non-commuting, non-perturbative* regime grows at least exponentially in β times something like an energy barrier — Davies/CKG bounds for stabilizer Hamiltonians are `λ ≥ O(N^{-1} exp(-2β ε̄))` (Temme 2014, [arXiv:1412.2858](https://arxiv.org/abs/1412.2858)); slow-mixing lower bounds for generic low-temperature Gibbs samplers are `T_mix = exp[Ω(n^α)]` (Gamarnik–Kiani 2024, [arXiv:2411.04300](https://arxiv.org/abs/2411.04300)); fast-mixing results are restricted to high-T, 1D, weakly interacting, or commuting cases (Bauerschmidt–Dagallier 2024, [arXiv:2202.02301](https://arxiv.org/abs/2202.02301); Bergamaschi–Chen 2025, [arXiv:2510.08533](https://arxiv.org/abs/2510.08533); Šmíd et al. 2025, [arXiv:2510.04954](https://arxiv.org/abs/2510.04954)). A separable `n^x · β^y` ansatz is *almost certainly* a leading-order fiction for any non-trivial Hamiltonian probed across β ∈ [5, 20]. The right question is "which ansatz family is the data compatible with?", not "what is x?".

## 1. Primary recommendation: Bayesian model comparison in log τ-space

Fit the data in `log τ = f(n, β; θ) + ε`, with `ε ~ N(0, σ²(n,β))` and σ(n,β) inferred jointly. Use Turing.jl (Julia, fits the project's stack — fully compatible with `LsqFit.jl` outputs and `BenchmarkTools` rigs) for HMC/NUTS posterior sampling of the parameters of each candidate form. Candidate families:

- **M0 — separable power law:** `log τ = c + x log n + y log β`
- **M1 — separable power × Arrhenius:** `log τ = c + x log n + α β`
- **M2 — gap-closing in n:** `log τ = c + x log n + α β · g(n)` with `g(n) ∈ {1, n^γ, log n}`
- **M3 — single combined variable:** `log τ = c + a · h(n β)` with `h` either log (gives `(nβ)^a`) or affine (gives Arrhenius in `nβ`)
- **M4 — broken / smooth-crossover:** `log τ = c + x log n + α β · σ(β−β★) + y log β · (1−σ(β−β★))` with `σ` a smooth-step.

Compare with **AICc** (Hurvich–Tsai correction, recommended by Burnham–Anderson when `n/p < 40`; we are at 36/3..6 ≈ 6..12) and **PSIS-LOO** (Vehtari–Gelman–Gabry 2017, [arXiv:1507.02646](https://arxiv.org/abs/1507.02646)) via `Turing` + the `ParetoSmoothedImportanceSampling.jl` interface (or `loo.psis` in R / `arviz.loo` in Python). PSIS-LOO is more robust than naive K-fold at n ≈ 36 and emits a Pareto-k̂ diagnostic that flags individual points whose removal changes the conclusion — exactly the leverage diagnostic you want when each (n, β) cell takes 4^n time.

Why AICc *and* LOO: AICc is closed-form once the model is fit and gives a fast first ranking; LOO is data-driven and catches the case where the AIC-best model is over-fit to a few high-leverage points (here, most likely n=11 at β=20). When AICc and LOO agree, the conclusion is defensible; when they disagree, that itself is the headline finding for the thesis ("our data don't discriminate between M0 and M1; we need more low-β or higher-n points").

This is the right primary tool because it (i) explicitly handles the small-sample regime that breaks classical FSS collapse, (ii) returns honest credible intervals on each parameter from the posterior — no bootstrap needed — (iii) tests the *form* of the ansatz, not just its parameters, and (iv) propagates the heteroscedastic noise model rigorously.

## 2. Backup / cross-check: Houdayer–Hartmann data collapse on `(n, β)` slices

After the Bayesian model comparison selects a form, perform a model-free corroboration using the **Houdayer–Hartmann quality function** (Houdayer–Hartmann 2004; coded in [pyfssa](https://pyfssa.readthedocs.io/) and [autoScale.py](https://github.com/omelchert/autoScale)). The reduced χ² quantity

S = (1/𝒩) Σᵢⱼ (yᵢⱼ − Y(xᵢⱼ))² / (σ²ᵢⱼ + σ²_Y)

attains its minimum ≈ 1 when the rescaled data lie on a single master curve to within the error bars (Bhattacharjee–Seno 2001, [cond-mat/0102515](https://arxiv.org/abs/cond-mat/0102515)). The two diagnostic plots:

- **β-collapse**: plot `τ · n^{−x}` versus β for all n with x fixed at the posterior mean — should overlap if the n-dependence is genuinely a clean power law.
- **n-collapse**: plot `τ · exp(−α β)` (or `· β^{−y}`, whichever M0 vs M1 won) versus n for all β.

A collapse value `S ≈ 1` confirms the ansatz; `S ≫ 1` invalidates it, and the *shape of the residual* (systematic curvature vs random scatter) tells you which direction to extend the model. Houdayer–Hartmann is *independent of* the regression code and uses only the data + errorbars, so it is a genuine cross-check, not a re-fit. Caveat: pyfssa and BSA both assume a critical-point form `(L^{1/ν}(ρ−ρ_c), L^{−ζ/ν})`; you'll need to rewrite the collapse on top of the bare H–H quality function rather than use the canned `autoscale()` routine. ~50 lines of Julia.

## 3. Software

**Recommended (Julia, matches project):**
- `Turing.jl` ([github.com/TuringLang/Turing.jl](https://github.com/TuringLang/Turing.jl)) for the Bayesian regression. Already battle-tested for nonlinear models, supports custom priors and likelihoods. Install: one `Pkg.add`.
- `LsqFit.jl` ([JuliaNLSolvers/LsqFit.jl](https://github.com/JuliaNLSolvers/LsqFit.jl)) for the MLE point estimates as starting values for Turing — Levenberg–Marquardt with `standard_error` and `estimate_covar` for first-cut uncertainties. Already in your dependency tree.
- `MCMCChains.jl` for posterior summaries; `Distributions.jl` for priors. Both Turing-native.
- Roll your own ~50-line H–H collapse — the formula is simple, and there's no Julia FSS package as mature as pyfssa.

**Alternative (Python, if needed):**
- `PyMC` for full Bayesian regression with NUTS, or `emcee` (Foreman-Mackey et al., [arXiv:1202.3665](https://arxiv.org/abs/1202.3665)) for an affine-invariant ensemble sampler that handles modest multi-modality well.
- `pyfssa` ([readthedocs](https://pyfssa.readthedocs.io/)) for the H–H collapse if you want the canned implementation.
- `arviz` for LOO/WAIC and posterior diagnostics.

**Avoid for primary analysis:**
- `PySR` / symbolic regression: with 36 points and a strong physical prior on the *form*, PySR will at best rediscover what you already suspect and at worst over-fit to high-order Taylor terms (a documented failure mode — Cranmer et al. 2024, [link.springer.com/10.1007/s10710-024-09503-4](https://link.springer.com/article/10.1007/s10710-024-09503-4)). It also has reproducibility issues under multi-processing. Useful as a *last-resort exploratory* tool if M0–M4 all fail PSIS-LOO.
- `SINDy` / `pysindy`: governs ODE / PDE *time evolution* — not a scaling-law tool. Wrong instrument.
- BSA / Bayesian scaling analysis: assumes a critical point and a scaling function; mismatched to our smooth, off-critical sweep.

## 4. Diagnostic plots

Six plots will tell the story:

1. **Raw `log τ` vs β at fixed n** (one line per n) and **vs n at fixed β** (one line per β). Linearity in `log β` favours power law; linearity in β itself favours Arrhenius. The shape *across* the n curves tells you whether the slopes are themselves n-dependent (i.e., coupled).
2. **Posterior pair plot** for `(x, y)` from M0, or `(x, α)` from M1 — banana shapes warn of aliasing between exponents (a well-known FSS pitfall, e.g. Sandvik 2010, [arXiv:1101.3281](https://arxiv.org/abs/1101.3281)).
3. **Residual map** `log τ − f(n, β; θ̂)` on the (n, β) grid as a heatmap. Genuine residuals are random; systematic stripes in the n or β direction kill the ansatz.
4. **Houdayer–Hartmann collapse plots** (Section 2) with reported S.
5. **PSIS-LOO Pareto-k̂ diagnostic** ([mc-stan.org/loo](https://mc-stan.org/loo/reference/pareto-k-diagnostic.html)): any k̂ > 0.7 indicates an outlier that the model can't generalize past — usually the corner of the grid. Flagging that corner is itself a result.
6. **Cook's distance / leave-one-out parameter drift** for the LM (LsqFit) fit, to identify high-leverage cells before the Bayesian re-fit.

## 5. Data needs

36 points is *just* enough to fit 2–3 parameters with ≈ 10% error on each, **provided** the data span enough of (n, β) for the exponents not to alias. Concretely:

- **Extend β downwards:** add β ∈ {1, 2, 3} so the full β range is roughly a decade and a half. The CKG complexity is `~ poly(β)` for high-T results (Bauerschmidt 2024; Šmíd 2025) and the Arrhenius regime is `exp(α β)`; without points well below the suspected crossover β★, M0 and M2 are essentially indistinguishable in 36-point regression. **This is the highest-value extension.**
- **Hold off on more n unless β is extended first.** The cost is 4^n; doubling β-range at fixed n is cheap, doubling n at fixed β-range is not.
- **One full diagonal (n, β) = constant** would help: e.g. take β = 5, 10, 20 at n = 7 *and* n = 9 (a 6-point diagonal). Lifted MCMC literature (Diaconis–Holmes–Neal 2000, [semanticscholar.org/...](https://www.semanticscholar.org/paper/ANALYSIS-OF-A-NONREVERSIBLE-MARKOV-CHAIN-SAMPLER-Diaconis-Holmes/e11681da5be491aaf80184f2127ec7ca45d42220)) gives strong reason to also look along the `nβ` combined axis.
- **Replicate corners** of the grid (n=11, β=20 and n=3, β=5) at independent random seeds if there is any stochasticity, to anchor the heteroscedastic noise model. Even if your Krylov estimator is fully deterministic, repeating with different krylovdim values is a cheap proxy for the numerical-floor variance.

Target ≈ 45–50 points if you can afford it. Above that, the AICc penalty becomes irrelevant and standard AIC suffices — but the qubit-cost ratio for 36 → 50 is far better than 11 → 12 in n.

## 6. Pitfalls

1. **Exponent aliasing along a narrow range.** With β over only a factor 4 (5→20) and `log β` over a factor ≈ 2.0, the M0 slope `y` and the M1 slope `α` are *strongly correlated*: both reproduce the same coarse trend. Solution: extend β downwards (Section 5), and explicitly inspect the posterior `corr(y, α)` — values > 0.9 mean the data does not discriminate.
2. **Mistaking exp(αβ) for high-y power law.** Over β ∈ [5, 20], a clean Arrhenius with `α ≈ 0.15` produces τ that grows by a factor ≈ 9; a power law `β^{1.6}` gives a factor ≈ 9 too. AICc on M0 vs M1 within this range will not split them. Either extend β downward or add the combined-variable diagnostic plot 1.
3. **Narrow-n bias.** n = 3 through 11 spans `log n ≈ 1.1 → 2.4` — a factor ~2.2 — which is shorter than the β range in log units. Exponents extracted from such short ranges are biased toward whichever sign of curvature dominates. Sandvik (2010) and the FSS literature recommend at least one decade in the scaling variable; we have less than that. The Bayesian credible intervals will *report* the bias correctly as wide posteriors, but only if you don't fix corrections to scaling at zero.
4. **Heteroscedasticity ignored.** A flat-σ likelihood will give too-tight intervals on the corner points. Always estimate `σ(n, β)` jointly, with a weak prior favoring `σ` growing with n (since Krylov error scales roughly with `δ²/λ`, see memory `project_delta_t0_floor.md`). If you don't have replicated runs to anchor σ, set a hierarchical prior `log σ = σ₀ + σ_n · n + σ_β · β` with weakly informative priors.
5. **Lumping the smooth-Metropolis kink scale into the fit.** Per memory `smooth_metro_s_beta_scaling_qf_96o.md`, the kink smoothing `σ·√s` is *itself* β-dependent through the default `s = (0.05/σ)²`. If σ ∝ 1/β along the sweep diagonal, the "fairness" of the comparison can spuriously add an effective β-dependence to τ_mix. Either freeze (σ, s) explicitly or include them as additional covariates in M0–M4.
6. **Reading the conclusion off the AICc winner.** AICc and BIC frequently disagree on small samples; published guidance is to *report* both, plus the absolute ΔAICc (changes < 4 are weak preference, < 2 are essentially indistinguishable — Burnham–Anderson 2002). For the thesis, list the top 2–3 models with their weights, not a single winner.

## 7. What does the literature predict the answer should be?

This frames the discussion section. Tabulating the rigorous results that bracket our setting:

- **Generic non-commuting, low-T:** No polynomial bound exists. Slow-mixing lower bounds are `T_mix = exp[Ω(n^α)]` for many natural Hamiltonians (Gamarnik–Kiani 2024, [arXiv:2411.04300](https://arxiv.org/abs/2411.04300)). For stabilizer codes the lower bound is `exp[Ω(n)]`; for 2D TFIM `T_mix = exp[n^{1/2−o(1)}]`. Temperature enters through the Arrhenius factor `exp(2β ε̄)` (Temme 2014, [arXiv:1412.2858](https://arxiv.org/abs/1412.2858)).
- **1D short-range:** `gap(L) = Θ(1)` independent of n; depth polylog n (Bergamaschi–Chen 2025, [arXiv:2510.08533](https://arxiv.org/abs/2510.08533)). β-dependence not made explicit in the abstract but implicit in the constants.
- **Weakly interacting:** Rapid mixing `O(polylog n)` at all constant temperature (Šmíd–Meister–Berta–Bondesan 2025, [arXiv:2510.04954](https://arxiv.org/abs/2510.04954)).
- **High-T commuting:** Davies generators rapidly thermalize, `T_mix = O(log n)` (Kastoryano–Temme 2013, [arXiv:1207.3261](https://arxiv.org/abs/1207.3261); Bardet et al. 2024). At low T, the gap can close as `λ ≥ N^{−1} exp(−2β ε̄)`.
- **Classical Glauber (sanity check):** For the d-dim Ising model, polynomial mixing at the critical β (Lubetzky–Sly 2012, [arXiv:1001.1613](https://arxiv.org/abs/1001.1613)); cutoff at high T (Lubetzky–Sly cutoff papers); *exponential* in n at low T due to metastability (mean-field Glauber dynamics literature).
- **Non-reversible / lifted:** Diaconis–Holmes–Neal (2000) gave the first analytic example where lifting cuts mixing time by `√` of the reversible mixing time, hinting that the right combined variable in our setting may be `√(nβ)` rather than `nβ`.

The most honest expectation is: at β = 5 we are likely still close to a high-temperature regime where M0 (separable power law) is OK to within errorbars; at β = 20 we are in the Arrhenius / metastability regime where M2 should win. If the data crosses over within β ∈ [5, 20], M4 (broken-power) will dominate; if it doesn't, β extension (Section 5) is mandatory to even see the crossover.

## 8. Citation table

### Finite-size scaling methodology
- **Privman (1990)**, *Finite Size Scaling and Numerical Simulation of Statistical Systems*, World Scientific, [worldscientific.com/...](https://www.worldscientific.com/worldscibooks/10.1142/1011). Canonical textbook on FSS for Monte Carlo. Background for chapter 5.
- **Sandvik (2010)**, *Computational Studies of Quantum Spin Systems*, [arXiv:1101.3281](https://arxiv.org/abs/1101.3281). Sections on FSS, cutoff and data collapse for QMC; explicit cautions on narrow scaling ranges.
- **Bhattacharjee & Seno (2001)**, "A measure of data-collapse for scaling", [cond-mat/0102515](https://arxiv.org/abs/cond-mat/0102515). Defines the residual-distance metric underlying Houdayer–Hartmann.
- **Harada (2011)**, "Bayesian inference in the scaling analysis of critical phenomena", *Phys. Rev. E* 84, 056704, [arXiv:1102.4149](https://arxiv.org/abs/1102.4149). Gaussian-process regression for scaling functions with credible intervals.
- **Harada (2015)**, "Kernel method for corrections to scaling", *Phys. Rev. E* 92, 012106, [arXiv:1410.3622](https://arxiv.org/abs/1410.3622). Extension handling subleading corrections — relevant if `n=3` produces visible deviations.
- **pyfssa** (Sorge), [readthedocs.io](https://pyfssa.readthedocs.io/) and [github](https://github.com/andsor/pyfssa). Python reference implementation of Houdayer–Hartmann + scaling-function fitting.
- **autoScale.py** (Melchert), [github.com/omelchert/autoScale](https://github.com/omelchert/autoScale). The older, standalone tool that pyfssa builds on.

### Model selection / regression
- **Burnham & Anderson (2002)**, *Model Selection and Multimodel Inference*, Springer; review at [warnercnr.colostate.edu](https://sites.warnercnr.colostate.edu/wp-content/uploads/sites/73/2017/05/Burnham-and-Anderson-2004-SMR.pdf). Definitive reference for AICc weights and multimodel inference when n/p is small.
- **Vehtari, Gelman, Gabry (2017)**, "Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC", *Stat. Comp.*, [arXiv:1507.02646](https://arxiv.org/abs/1507.02646). PSIS-LOO; Pareto-k̂ diagnostic; recommended for small-sample model comparison.

### Gibbs-sampler mixing theory (quantum)
- **Kastoryano & Temme (2013)**, "Quantum logarithmic Sobolev inequalities and rapid mixing", *J. Math. Phys.* 54, 052202, [arXiv:1207.3261](https://arxiv.org/abs/1207.3261). Sets the rapid-mixing framework.
- **Temme (2014)**, "Thermalization time bounds for Pauli stabilizer Hamiltonians", *Comm. Math. Phys.*, [arXiv:1412.2858](https://arxiv.org/abs/1412.2858). Explicit `λ ≥ N^{−1} exp(−2β ε̄)` at low T.
- **Chen, Kastoryano, Brandão, Gilyén (2023)**, "Quantum thermal state preparation", [arXiv:2303.18224](https://arxiv.org/abs/2303.18224). The smooth-Metropolis (CKG) sampler — abstract states cost depends polylogarithmically on mixing time and precision; explicit closed-form n / β scaling is not in the abstract and requires reading the body.
- **Ding, Li, Lin (2024)**, "Efficient quantum Gibbs samplers with KMS detailed balance condition", [arXiv:2404.05998](https://arxiv.org/abs/2404.05998). DLL Lindbladian; complementary to CKG.
- **Gamarnik & Kiani (2024)**, "Slow mixing of quantum Gibbs samplers", [arXiv:2411.04300](https://arxiv.org/abs/2411.04300). Lower bounds `T_mix = exp[Ω(n^α)]` for K-SAT, spin glasses, stabilizer codes, 2D TFIM.
- **Bergamaschi & Chen (2025)**, "Fast mixing of quantum spin chains at all temperatures", [arXiv:2510.08533](https://arxiv.org/abs/2510.08533). 1D, system-size-independent gap, polylog-depth circuit.
- **Šmíd, Meister, Berta, Bondesan (2025)**, "Rapid mixing of quantum Gibbs samplers for weakly-interacting quantum systems", [arXiv:2510.04954](https://arxiv.org/abs/2510.04954). Polylog(n) mixing for weakly-interacting systems.

### Classical mixing for context
- **Lubetzky & Sly (2010)**, "Cutoff for the Ising model on the lattice", *Invent. Math.*, [arXiv:0909.4320](https://arxiv.org/abs/0909.4320). Cutoff at high T.
- **Lubetzky & Sly (2012)**, "Critical Ising on the square lattice mixes in polynomial time", *Comm. Math. Phys.*, [arXiv:1001.1613](https://arxiv.org/abs/1001.1613). Polynomial mixing at criticality.
- **Bauerschmidt & Dagallier (2024)**, "Log-Sobolev inequality for near critical Ising models", *Comm. Pure Appl. Math.* 77, 2568, [arXiv:2202.02301](https://arxiv.org/abs/2202.02301). Polynomial LSI for Ising in d ≥ 5 up to T_c.
- **Diaconis, Holmes, Neal (2000)**, "Analysis of a nonreversible Markov chain sampler", *Ann. Appl. Prob.*, [semantic scholar](https://www.semanticscholar.org/paper/ANALYSIS-OF-A-NONREVERSIBLE-MARKOV-CHAIN-SAMPLER-Diaconis-Holmes/e11681da5be491aaf80184f2127ec7ca45d42220). Lifting; `√` improvement; relevant for combined-variable analyses.
- **Bravyi & Gosset (2017)**, "Polynomial-time classical simulation of quantum ferromagnets", *PRL* 119, 100503, [arXiv:1612.05602](https://arxiv.org/abs/1612.05602). Classical baseline for what's *not* a quantum advantage regime.

### Software (regression / Bayesian / FSS)
- **Turing.jl**, [github.com/TuringLang/Turing.jl](https://github.com/TuringLang/Turing.jl). Primary recommendation; HMC/NUTS in pure Julia.
- **LsqFit.jl**, [github.com/JuliaNLSolvers/LsqFit.jl](https://github.com/JuliaNLSolvers/LsqFit.jl). Levenberg–Marquardt MLE; supplies starting values and quick first-cut SE.
- **PyMC**, [pymc.io](https://github.com/pymc-devs/pymc) — Python backup if Julia tooling is missing a feature.
- **emcee** (Foreman-Mackey et al. 2013), [arXiv:1202.3665](https://arxiv.org/abs/1202.3665), [emcee.readthedocs.io](https://emcee.readthedocs.io/). Affine-invariant ensemble sampler; minimal dependencies.
- **loo / arviz**, [mc-stan.org/loo](https://mc-stan.org/loo/). PSIS-LOO + WAIC + Pareto-k̂ diagnostics.
- **BSA (Bayesian scaling analysis)**, Harada, [kenjiharada.github.io/BSA](https://kenjiharada.github.io/BSA/) and [github.com/KenjiHarada/BSA](https://github.com/KenjiHarada/BSA). Reference if you want the GP-regression alternative — but BSA requires a critical point, mismatched to our use.

### Tools to *avoid* for the primary analysis (kept for transparency)
- **PySR** ([github.com/MilesCranmer/PySR](https://github.com/MilesCranmer/PySR)) — symbolic regression. With 36 points and a strong physical prior, more likely to find a high-order Taylor approximation than the true form (documented in [Springer 2024 review](https://link.springer.com/article/10.1007/s10710-024-09503-4)); also has reproducibility issues under multi-processing.
- **AI Feynman** ([science.org](https://www.science.org/doi/10.1126/sciadv.aay2631)) — also symbolic; designed for *much* larger datasets and clean noise.
- **pysindy / SINDy** ([github.com/dynamicslab/pysindy](https://github.com/dynamicslab/pysindy)) — for governing-equation discovery from *time series*, not scaling laws.

## Action items for the numerics chapter

1. Extend the β sweep down to β ∈ {1, 2, 3} at n ∈ {3, 5, 7, 9} (about 12 new points, cheap). This is the single most leveraged change.
2. Implement the Turing.jl regression for M0–M4 with explicit heteroscedastic noise; report AICc and PSIS-LOO weights for each.
3. Write a ~50-line Houdayer–Hartmann data-collapse routine and produce the two collapse plots after the winning model is identified.
4. Report `(x, α)` posteriors with credible intervals and the corr(x, α) value; if corr > 0.9, state explicitly that the data do not discriminate the form and that the conclusion is conditional on which model family is assumed.
5. Tie the discussion section to the bracketing rigorous results in Section 7: state which regime the data live in, and which conjecture (rapid mixing, Arrhenius, or slow mixing) it is compatible with.
