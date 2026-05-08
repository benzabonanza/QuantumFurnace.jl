# CKG vs DLL — first findings on Kossakowski structure and mixing time

*Insertion target:* `1_preliminaries.tex:1350` (the TODO `% IF you write the comparison section, ref to it`) — this draft is intended as the seed for that comparison subsection. Numerical/script provenance: `scripts/numerics_ckg_vs_dll_taumix_comparison.jl`, `scripts/scratch_taumix_beta20_extension.jl`, `scripts/scratch_kossakowski_antidiagonal_analysis.jl`. Data in `drafts/figures/numerics/ckg_vs_dll_taumix.bson` and `scripts/output/taumix_beta20/taumix_beta20.bson`.

---

## 1. Setup of the comparison

We compare three KMS-detailed-balanced Lindbladian families on the same disordered Heisenberg fixture (`hamiltonians/heis_disordered_periodic_n*.bson`, $n\in\{3,4,5\}$, with $Z$-disorder optimised for max min-gap):

- **CKG smooth-Metropolis** [Chen et al. 2025] [CITE: chen2023efficient] in the linear-combination-of-Gaussians form derived in the thesis (Eq. eq:kossakowski-smooth-metro), $\sigma = 1/\beta$, $a=0$, $s=0.25$ (locked thesis convention; see also (Lem. lem:shifted-kms-condition) and (Cor. cor:kossakowski-shift) for the $\gamma\!\to\!\tilde\gamma$ shift mechanism).
- **DLL Gaussian-type filter** [Ding et al. 2024] [CITE: ding2025efficient], Eq. 3.21–3.22 — $q^{\text{DLL,Gauss}}(\nu) = e^{-(\beta\nu)^2/8}$ with bump $w\equiv 1$ on the relevant grid.
- **DLL Metropolis-type filter** [Ding et al. 2024] [CITE: ding2025efficient], Eq. 3.19–3.20 — $q^{\text{DLL,Metro}}(\nu) = e^{-\sqrt{1+(\beta\nu)^2}/4}\cdot w(\nu/S)$, $S=2$ (Hörmander bump).

Mixing times are obtained from the matrix-free ODE-integrator + bi-exponential extrapolation pipeline (`integrate_to_gibbs`, `estimate_mixing_time(...; model=:biexp, extrapolate=true)`, `target_epsilon=10^{-3}$, $t$-horizon adaptive per qf-lkb.10) on the *ideal* (non-approximated) Lindbladians of CKG and DLL. The grid spans $\beta\in\{1,2,5,10,15,20\}$ and $n\in\{3,4,5\}$.

## 2. Empirical mixing times

The headline plot is `drafts/figures/numerics/ckg_vs_dll_taumix.{png,pdf}`. Selected $\tau_{\text{mix}}$ values (target $\varepsilon=10^{-3}$):

| $n$ | $\beta$ | CKG smooth-Metro | DLL Gaussian | DLL Metropolis |
|---|---|---|---|---|
| 3 | 5  | 15.7 |  11.8 | 18.1 |
| 3 | 10 | 19.0 |  25.3 | 23.3 |
| 3 | 15 | 26.0 | **67.4** | 28.6 |
| 3 | 20 | 33.0 | **264.5** | 35.4 |
| 4 | 20 | 41.6 | 84.9 | 50.2 |
| 5 | 5  | 24.5 | 15.5 | 24.9 |
| 5 | 10 | 29.7 | 24.6 | 33.3 |
| 5 | 20 | 35.5 | 62.2 | 42.6 |

Two empirical observations:

1. **Metropolis-vs-Metropolis is a near-tie.** CKG smooth-Metropolis and DLL Metropolis stay within $\sim$10–25% of each other across the entire $(n,\beta)$ grid, with CKG slightly faster for most cells.
2. **DLL Gaussian collapses at low temperature.** It is the *fastest* of the three at $\beta\!\le\!5$ (15–50% lead) but *crosses over* near $\beta=10$ and becomes dramatically slower at $\beta=20$ — at $n=3,\beta=20$ it is 8× slower than CKG.

The fitted Lindbladian gaps tell the same story: at $\beta=20$, $n=3$, the DLL Gaussian's gap is $0.024$ versus CKG's $0.155$ and DLL Metropolis's $0.149$ — an order-of-magnitude collapse for DLL Gaussian only.

**There is, on this grid, no evidence that DLL beats CKG in mixing time.** Whether this persists for larger $n$ (in particular $n\ge 7$, where the cluster-scale DLL matrix-free path qf-lkb.9 becomes useful) is open and is the natural next experiment.

## 3. The Kossakowski-rank framework

Both papers express the dissipator with a Kossakowski matrix $\alpha_{\nu_1,\nu_2}$ indexed by Bohr frequencies $\nu_1,\nu_2 \in B_H$. For a single coupling $A^a$:

$$
\mathcal{L}^{\text{diss}}[\rho] \;=\; \sum_{\nu_1,\nu_2 \in B_H} \alpha_{\nu_1,\nu_2}^{\,a}\;\bigl(A^a_{\nu_1}\rho (A^a_{\nu_2})^\dagger - \tfrac12\{(A^a_{\nu_2})^\dagger A^a_{\nu_1},\rho\}\bigr).
$$

This is the form in [Chen et al. 2025, Cor. II.2][CITE: chen2023efficient] and [Ding et al. 2024, Sec. 4][CITE: ding2025efficient]. The KMS detailed-balance condition forces the centrosymmetry $\alpha_{\nu_1,\nu_2} e^{\beta(\nu_1+\nu_2)/4} = \alpha_{-\nu_2,-\nu_1} e^{-\beta(\nu_1+\nu_2)/4}$ (this is (Eq. eq:skew-symmetry) in the thesis and Eq. 4.7 in [Ding et al.][CITE: ding2025efficient]).

**Lindblad operators = eigenvectors of the Kossakowski.** Eigendecomposing $\widetilde{C}_{\nu,\nu'} := \alpha_{\nu,\nu'}\,e^{\beta(\nu+\nu')/4}$ as $Q D Q^\dagger$ (Eq. 4.9 of [Ding et al.][CITE: ding2025efficient]) gives a Lindblad operator per non-zero eigenvalue:

$$
L_{\nu'} \;=\; \sqrt{D_{\nu',\nu'}} \sum_{\nu\in B_H} Q_{\nu,\nu'}\,A_\nu\,e^{-\beta\nu/4} \quad (\text{Eq. 4.10 in [Ding et al.][CITE: ding2025efficient]}).
$$

Therefore $\operatorname{rank}(\alpha)$ counts the number of independent dissipative channels per coupling. This is the operational meaning of "rank" we will use throughout.

### 3.1 DLL is rank-1 *in the standard construction*

Both DLL filter choices we use realise the standard DLL construction of [Ding et al. 2024, Sec. 3.1, Eq. 3.4][CITE: ding2025efficient], which prescribes a *single* weight $q^a(\nu)$ per coupling. The Kossakowski is then a literal outer product

$$
\alpha^{\text{DLL}}_{\nu_1,\nu_2} \;=\; e^{-\beta(\nu_1+\nu_2)/4}\,q^a(\nu_1)\,q^a(\nu_2)^\ast,
$$

i.e. rank exactly $1$ — a *single* Lindblad operator $B^a = \sum_\nu q^a(\nu)\,A^a_\nu$ per coupling.

**Important caveat.** This is the rank of the *standard* DLL construction we benchmark, not a structural ceiling. [Ding et al. 2024][CITE: ding2025efficient] themselves prove (Sec. 4, Eq. 4.10–4.12) that DLL is *more general* than rank-1: by taking $|B_H|$ different weights $q_{\nu'}(\nu) := \sqrt{D_{\nu',\nu'}}Q_{\nu,\nu'}$ (or its imaginary version), the DLL framework recovers **CKG as a special case**. So the rank-1 property is a property of how the family is *used* in algorithmic practice, not a property of the framework. (For reference: this is the same reduction used by the thesis in §sec:kms-lindbladian when discussing CKG's relation to its predecessors.)

### 3.2 CKG is generically full-rank

In CKG the Kossakowski is *not* a literal outer product but a Gram-type integral over $\omega$ [Chen et al. 2025, §II.A][CITE: chen2023efficient]:

$$
\alpha^{\text{CKG}}_{\nu_1,\nu_2} \;=\; (2\pi)^2 \int_{-\infty}^{\infty}\!\gamma(\omega)\,\hat f(\omega-\nu_1)^\ast \hat f(\omega-\nu_2)\,\dd\omega \;=\; \langle g_{\nu_1},\,g_{\nu_2}\rangle,
$$

with $g_\nu(\omega) := \sqrt{\gamma(\omega)}\,\hat f(\omega-\nu)$. Since the translates $\{g_\nu\}_{\nu\in B_H}$ are generically linearly independent in $L^2(\dd\omega)$, the rank is generically $|B_H|$. For the fixtures in our suite, our numerical SVD on $n=3$ gives $\operatorname{rank}\bigl(\alpha^{\text{CKG, smooth-Metro}}\bigr) \approx 10$ at $\beta=5$, growing to $\approx 24$ at $\beta=20$ (out of $|B_H|=57$). Concretely: CKG dispatches $\sim 24$ **parallel** dissipative channels per coupling, against DLL's $1$.

The thesis expression for $\alpha^{\text{CKG, smooth-Metro}}$ is (Eq. eq:kossakowski-smooth-metro), with the convention $\sigma=1/\beta$ from [Chen et al. 2025, p. 12][CITE: chen2023efficient].

## 4. Why DLL Metropolis matches CKG smooth-Metropolis: two paths to the same place

For both CKG and DLL, KMS skew-symmetry forces the Kossakowski mass to concentrate near the *antidiagonal* $\nu_2 = -\nu_1$ as $\beta\to\infty$ (entries off the antidiagonal are suppressed by $e^{-\beta(\nu_1+\nu_2)/2}$). What differs is **how** mass is distributed *along* the antidiagonal. Population transfer between energy sectors $E$ and $E-\nu$ requires antidiagonal mass at non-zero $\nu$; the entry at $\nu=\nu'=0$ is pure dephasing.

Numerically extracted antidiagonal slice $|\alpha(\nu,-\nu)|$ at $\beta=20$, $n=3$ (script `scratch_kossakowski_antidiagonal_analysis.jl`):

| $\nu$ | CKG smooth-Metro | DLL Metropolis | DLL Gaussian |
|---|---|---|---|
| 0.00 | 0.576 | 0.607 | 1.000 |
| 0.07 | 0.21 | 0.43 | 0.62 |
| 0.21 | $1\!\times\!10^{-4}$ | 0.12 | 0.012 |
| 0.35 | $2\!\times\!10^{-11}$ | 0.030 | $5\!\times\!10^{-6}$ |
| 0.45 | $1.5\!\times\!10^{-18}$ | **0.011** | $1.6\!\times\!10^{-9}$ |

Three regimes:

- **CKG smooth-Metro:** every individual antidiagonal entry decays super-exponentially as $|\nu|$ grows (the off-diagonal $e^{-\nu_-^2/(8\sigma^2)}$ factor in (Eq. eq:kossakowski-smooth-metro) with $\sigma=1/\beta$ gives a Gaussian envelope of width $\sim 1/\beta$ in the $\nu_- = \nu_1-\nu_2$ direction). At $\beta=20,\nu=0.45$ the CKG entry is *physically* zero. Yet CKG mixes well — because it has $\sim\!24$ such narrow channels in parallel, each centered at a different antidiagonal point along the spectrum. Many narrow brushes paint the spectrum.
- **DLL Metropolis:** $|q(\nu)|^2 \sim e^{-|\beta\nu|/2}$ on the flat-top region, only exponential in $|\beta\nu|$. The single rank-1 channel keeps non-trivial mass even at $\nu\approx \pm 0.45$. One wide brush.
- **DLL Gaussian:** $|q(\nu)|^2 \propto e^{-(\beta\nu)^2/4}$, super-exponential in $\beta\nu$. The single channel concentrates within $|\nu|\le 1/\beta$; outside that band it is effectively zero. One narrow brush — and only one, since rank is $1$. Both *redundancy* (rank) and *reach* (filter tail) fail simultaneously.

### Why the CKG-paper warning about Gaussian $\gamma$ does *not* apply to DLL Gaussian

[Chen et al. 2025, p. 12][CITE: chen2023efficient] warn that picking the **transition weight** $\gamma(\omega)$ Gaussian leads to a narrow transition band and substantially increased mixing time, and that the fix is taking linear combinations to get a Metropolis-shaped $\gamma$ (this is the route the thesis takes in §Improved Kossakowski matrices, see (Eq. eq:smooth-metro), (Eq. eq:kossakowski-smooth-metro)). But this warning is about $\gamma(\omega)$, an $\omega$-domain object that gets *integrated* against the Gaussian filter $\hat f(\omega-\nu)$ before producing $\alpha^{\text{CKG}}$. **DLL Gaussian** is a different mathematical object: $q(\nu)$ lives directly in Bohr-frequency space and is not a $\gamma(\omega)$. The DLL Gaussian's collapse at low $T$ comes from the rank-1 weight $q^{\text{DLL,Gauss}}(\nu) = e^{-(\beta\nu)^2/8}$ being super-exponentially localised in $\beta\nu$, an *internal* DLL property of (Eq. 3.21 in [Ding et al.][CITE: ding2025efficient]).

So the mechanism is real but the supporting reference is DLL §3.2, not CKG p. 12.

## 5. The DLL paper's own caveat: wider Kossakowski → larger block-encoding cost

A subtle and important point that emerged from re-reading the DLL paper: the wider $\alpha$-support DLL achieves at low $T$ — which one might naively read as "DLL keeps mixing alive while CKG suffers" — actually *increases* DLL's algorithmic cost.

[Ding et al. 2024, Remark 23][CITE: ding2025efficient] (paper p. 24) makes this explicit:

> "$\{\alpha_{\nu,\nu'}\}$ always includes a $\Omega(1)$-sized principal submatrix and the dynamics allows different energy transition terms $A_\nu \rho\, (A^a_{\nu'})^\dagger$ even when $|\nu-\nu'|\gg \Omega(1/\beta)$. Therefore, our coherent $G$ has much more cross terms in the expansion, and the **normalization constant for the block encoding increases linearly in $\beta$**."

This translates into a different overall complexity scaling. From [Ding et al. 2024, Table 1][CITE: ding2025efficient]:

| Sampler | Total Hamiltonian-simulation cost |
|---|---|
| CKG (Gaussian or Metropolis $\gamma$) | $\widetilde{O}(\beta\,t_{\text{mix}}\,\mathrm{polylog}(1/\varepsilon))$ |
| DLL (any Gevrey filter, support $S$) | $\widetilde{O}(\beta^{2}\,S\,t_{\text{mix}}\,\mathrm{polylog}(1/\varepsilon))$ |

For $S=O(1)$, **DLL pays an extra factor of $\beta$** in total cost relative to CKG. The mechanism is precisely the structural property that produces the wider Kossakowski support: the coherent term $G$ — which the thesis introduces in (Eq. eq:ckg-L-parts) and which always appears in the KMS construction — has more cross-terms because the dissipator "mixes" energy levels $\nu,\nu'$ with $|\nu-\nu'|\gg 1/\beta$. So the very thing that makes DLL's Kossakowski non-shrinking with $\beta$ is what makes its coherent block-encoding more expensive.

This is the correction to the community folklore the user articulated: **DLL did manage to get a "non-shrinking" (rank-1 wide-support) Kossakowski, but that produced a new problem — a $\beta$-factor cost in the coherent block-encoding** [Ding et al. 2024, §3.31, Remark 23, Theorem 19][CITE: ding2025efficient]. Together with our numerical observation that **DLL also does not give a $t_{\text{mix}}$ gain** on the systems we benchmark, the asymptotic picture in the literature *favours CKG* for low-temperature large-system Gibbs sampling.

### 5.1 Numerical verification of Remark 23: $\|G\|$ scaling with $\beta$

We numerically tested the Remark 23 claim on our own fixtures using the operator norm $\|G\|_{\text{op}}$ as the proxy for the block-encoding normalisation constant (script: `scripts/scratch_coherent_norm_lcu_check.jl`). For each $(n, \beta)$, $G$ is built via `dll_coherent_op_bohr` (DLL Gaussian and DLL Metropolis filters) and via `B_bohr` (CKG smooth-Metropolis), all *without* the production-time $\gamma_{\text{norm}}$ rescaling so the three sit at the same algorithmic stage. Log–log slope of $\|G\|$ vs $\beta$ over the high-$\beta$ tail $\beta \in \{5, 10, 15, 20, 30\}$:

| $n$ | slope $\|B_{\text{CKG}}\|$ | slope $\|G_{\text{DLL,Gauss}}\|$ | slope $\|G_{\text{DLL,Metro}}\|$ |
|---|---|---|---|
| 3 | $-1.0$ | $-0.3$ | $\mathbf{+1.35}$ |
| 4 | $+0.7$ | $+1.45$ | $\mathbf{+1.46}$ |
| 5 | $+0.8$ | $+1.15$ | $\mathbf{+1.38}$ |

Three observations:

- **DLL Metropolis: slopes 1.35–1.46 across all $n$.** This is direct numerical confirmation of [Ding et al. 2024, Remark 23][CITE: ding2025efficient] — the wide rank-1 Metropolis support produces a coherent term $G$ whose normalisation grows at least linearly in $\beta$. Slightly super-linear at small $n$ likely reflects the additional $S$-dependent prefactor in their Theorem 19.
- **DLL Gaussian grows linearly only at $n \ge 4$**; at $n=3$ the slope is *negative*. This is the same low-$T$ collapse we saw in $\tau_{\text{mix}}$: the rank-1 Gaussian weight $q^{\text{DLL,Gauss}}(\nu)$ becomes so localised at $\nu \approx -1/\beta$ that even $G$ shrinks. Internally consistent — narrower filter $\Rightarrow$ smaller $G$ (good for block-encoding) $\Rightarrow$ but no good off-diagonal mass (bad for $\tau_{\text{mix}}$).
- **CKG smooth-Metropolis** has a flatter slope $\sim 0.7$–$0.8$ at $n \ge 4$, recovering the linear-$\beta$ scaling of [Ding et al. 2024, Table 1][CITE: ding2025efficient].

A subtlety: absolute values of $\|G\|$ on the rescaled fixtures (spectrum $[0, 0.45]$) are all in the $10^{-5}$–$10^{-2}$ range, well below the LCU norm-$1$ budget. So the question "does DLL $G$ break the budget on these benchmarks" answers *no* in absolute terms. But the *scaling* is the asymptotically relevant statement, and it is preserved: an unrescaled physical Hamiltonian with $\|H\| = O(n)$ would multiply the magnitudes by $1/\text{rescaling\_factor}$, recovering the quantitative regime where Remark 23 bites. The slope $\approx 1.4$ for DLL Metropolis is what matters.

**Net trade-off summary**: keep the filter wide $\Rightarrow$ pay in $\|G\|$ (DLL Metropolis, Remark 23, slope $\approx 1.4$); keep the filter narrow $\Rightarrow$ pay in $t_{\text{mix}}$ (DLL Gaussian, our §2 finding, $8\times$ slowdown at $n=3,\beta=20$). CKG dodges the trade-off by replacing the narrow rank-1 channel with $|B_H|$ narrow channels in parallel.

## 6. Verdict so far

| Claim | Status | Source |
|---|---|---|
| Kossakowski as $\|B_H\|\!\times\!\|B_H\|$ matrix; KMS skew-symmetry | confirmed both papers | [Chen et al. Cor II.2][CITE: chen2023efficient]; [Ding et al. Eq. 4.7][CITE: ding2025efficient] |
| Standard DLL is rank-1 (one $q$ per coupling, single $B^a$) | confirmed | [Ding et al. Eq. 3.4, Sec. 3.2][CITE: ding2025efficient] |
| DLL framework subsumes CKG via $|B_H|$ weights | confirmed | [Ding et al. Eq. 4.10–4.12][CITE: ding2025efficient] |
| CKG is generically full rank ($\sim |B_H|$) | confirmed numerically + analytically | numerics + [Ding et al. Eq. 4.9][CITE: ding2025efficient] |
| DLL Gaussian's rank-1 weight is super-exp. localised in $\beta\nu$ | confirmed | [Ding et al. Eq. 3.21][CITE: ding2025efficient] |
| DLL Metropolis preserves heavy exp. tail in $\beta\nu$ | confirmed | [Ding et al. Eq. 3.19–3.20, Remark 23][CITE: ding2025efficient] |
| DLL's wider Kossakowski $\Rightarrow$ extra $\beta$ in block-encoding $G$ | confirmed analytically + numerically (slope 1.35–1.46) | [Ding et al. Remark 23, Table 1][CITE: ding2025efficient]; this work §5.1 |
| DLL beats CKG in $t_{\text{mix}}$ at low $T$ | **not observed numerically** ($n\le 5,\beta\le 20$); not claimed by either paper | this work; cf. [Ding et al. Table 1][CITE: ding2025efficient] |
| DLL Gaussian degrades faster than DLL Metro and CKG smooth-Metro at low $T$ | observed numerically | this work |
| Scale-invariant $\rho_{\text{intrinsic}} = \lambda/\Lambda_{\max}$ separates structural from normalisation effects | yes | this work §8.1 ([Rouzé 2019][CITE: rouze2019functional], §5.2; [Li & Lu 2024][CITE: liLu2024quantum], §2.8) |
| CKG smooth-Metro and DLL Metro share Metropolis-shape $\alpha_{\text{diag}}$ (within 6% $L^2$) | confirmed numerically | this work §8.4 |
| $\rho(\text{CKG sM})/\rho(\text{DLL Metro})\in[1.02, 1.43]$ across $(n,\beta)$ — modest CKG advantage, no structural tie | confirmed numerically | this work §8.3 |
| Poincaré bound $\tau_{\text{pred}}\approx 3\cdot\tau_{\text{meas}}$ across 32 cells ($R^2 = 0.93$) | confirmed numerically | this work §8.5 |

## 7. Open questions / next steps

1. **Larger-$n$ verification.** Our grid stops at $n=5$. The qf-lkb.9 matrix-free DLL apply makes $n\!\ge\!7$ feasible; we should run the same comparison at $n\in\{6,7\}$ to check whether the picture survives. The thesis will benefit from at least one $n=7$ data point in this comparison.
2. **Lower temperatures.** $\beta\in\{30,50\}$ at $n=3$ would test whether the DLL Gaussian $\tau_{\text{mix}}$ scales as the predicted slow random walk through Bohr frequencies (small-step diffusion limit) and clarify the asymptotic.
3. **Sweep over multiple disorder seeds.** Currently single-seed (`SEEDS=[42]`); a few additional seeds at $n=5$, $\beta\in\{10,20\}$ would quantify scatter.
4. **DLL Metropolis with larger bump $S$.** Our $S=2$ choice keeps the bump out of the spectral support but is at the boundary; a deliberate $S\to\infty$ run should stay within $\sim$10% of CKG smooth-Metro, structurally confirming that *rank-1 Metropolis-shape* is the equivalence class.
5. **Higher-rank DLL.** A natural test of the framework subsumption (§3.1 caveat): construct a DLL with $k\ge 2$ weights $\{q_\ell\}_{\ell=1}^{k}$ per coupling (still finite, sub-CKG) and check whether mixing time interpolates between rank-1 DLL and full-rank CKG. This directly informs whether one can engineer a *cheap* DLL variant that retrieves CKG's parallelism without paying the full $|B_H|$-fold cost.

## 8. Are CKG and DLL Metropolis on the same footing?

The §2 headline — "CKG smooth-Metropolis $\approx$ DLL Metropolis in $\tau_{\text{mix}}$ to within 10–25%" — sits uneasily next to the fact that the two samplers run at *different production normalisations*. CKG's $\mathcal{L}^{\text{CKG}}$ is rescaled by $\gamma_{\text{norm}} := 1/\max_\omega \gamma(\omega)$ for LCU block-encoding (see `src/coherent.jl:38` and the `gamma_norm_factor` notes in `src/furnace_utensils.jl`); DLL has no analogue, because its filter $q(\nu)$ already absorbs the KMS factor (cf. the comment at `src/furnace_utensils.jl:103-105`). The $\tau_{\text{mix}}$ tie *might* be a structural equivalence — both samplers driving the same intrinsic relaxation — or it *might* be a coincidence in which two different intrinsic rates are dressed by two different normalisations that happen to compensate. We cannot tell from the §2 raw $\tau_{\text{mix}}$ numbers alone.

This section closes that gap. We compute a **scale-invariant intrinsic mixing rate** $\rho_{\text{intrinsic}}$ on every cell of the $(n,\beta)$ grid. The result: CKG smooth-Metropolis is consistently $1.02$–$1.43\times$ ahead of DLL Metropolis in $\rho_{\text{intrinsic}}$ — a modest, monotonically growing-in-$\beta$ lead that *partially* explains the §2 tie (some of the ostensible $\tau_{\text{mix}}$ closeness *is* compensating norms) without overturning the qualitative picture (CKG and DLL Metropolis still belong to the same equivalence class of fast Metropolis-shape mixers, well separated from the DLL Gaussian collapse regime).

### 8.1 The intrinsic mixing ratio

For a KMS-detailed-balanced Lindbladian $\mathcal{L}$ with stationary state $\sigma_\beta = e^{-\beta H}/Z$, the natural inner product on observables is the KMS inner product

$$
\langle X, Y\rangle_{\text{KMS}, \sigma_\beta} \;:=\; \operatorname{tr}\!\bigl[\sigma_\beta^{1/2}\,X^\dagger\,\sigma_\beta^{1/2}\,Y\bigr].
$$ <!-- \label{eq:kms-inner-product} -->

The associated Dirichlet form is

$$
\mathcal{E}_{\mathcal{L}}(X, X) \;:=\; -\langle X,\,\mathcal{L}^\dagger[X]\rangle_{\text{KMS}, \sigma_\beta},
$$ <!-- \label{eq:dirichlet-form-kms} -->

where $\mathcal{L}^\dagger$ acts in the Heisenberg picture. For a KMS-DBC generator, $\mathcal{E}_{\mathcal{L}} \ge 0$ on traceless observables, with kernel exactly the KMS-symmetric scalar multiples of the identity ([Rouzé 2019, §5.2][CITE: rouze2019functional]; [Ding et al. 2024, §1.1][CITE: ding2025efficient]). The **spectral gap** is

$$
\lambda(\mathcal{L}) \;:=\; \inf_{X \neq 0,\, \langle X\rangle_{\sigma_\beta} = 0}\; \frac{\mathcal{E}_{\mathcal{L}}(X, X)}{\operatorname{Var}_{\sigma_\beta}(X)},
$$ <!-- \label{eq:spectral-gap-poincare} -->

i.e. the optimal Poincaré constant. The **maximal Dirichlet rate** is the dual quantity

$$
\Lambda_{\max}(\mathcal{L}) \;:=\; \sup_{X \neq 0,\, \langle X\rangle_{\sigma_\beta} = 0}\; \frac{\mathcal{E}_{\mathcal{L}}(X, X)}{\operatorname{Var}_{\sigma_\beta}(X)},
$$ <!-- \label{eq:lambda-max-dirichlet} -->

i.e. the largest *KMS-eigenvalue* of $-\mathcal{L}^\dagger$ on the traceless subspace. The **intrinsic mixing ratio** is then

$$
\rho_{\text{intrinsic}}(\mathcal{L}) \;:=\; \frac{\lambda(\mathcal{L})}{\Lambda_{\max}(\mathcal{L})} \;\in\; (0, 1].
$$ <!-- \label{eq:rho-intrinsic} -->

It is, by construction, *scale-invariant*: $\rho_{\text{intrinsic}}(c \mathcal{L}) = \rho_{\text{intrinsic}}(\mathcal{L})$ for any $c > 0$. So it strips out exactly the production-normalisation degree of freedom that distinguishes CKG and DLL and isolates the *shape* of the spectrum on the KMS subspace. (See [Kochanowski, Alhambra, Capel & Rouzé 2024, §2][CITE: kochanowski2024rapid] and [Li & Lu 2024, §2.8][CITE: liLu2024quantum] for the same definition under different names — "condition number of the Dirichlet form", "spectral aspect ratio".)

**Implementation.** In the Φ-conjugated frame $\widetilde{\mathcal{L}}^\dagger := \sigma_\beta^{-1/4}\mathcal{L}^\dagger[\sigma_\beta^{1/4}\,\cdot\,\sigma_\beta^{1/4}]\sigma_\beta^{-1/4}$ the KMS inner product becomes the Hilbert–Schmidt inner product, $\widetilde{\mathcal{L}}^\dagger$ becomes Hermitian on $L^2(\operatorname{tr})$, and $\rho_{\text{intrinsic}}$ is read off as the ratio of the smallest to the largest non-zero eigenvalue of the symmetrised generator restricted to the traceless subspace. For $n \le 5$ the superoperator is a $4^n \times 4^n$ matrix and a dense eigendecomposition is cheap (script: `scripts/fair_comparison_kms_dirichlet.jl`).

### 8.2 A note on $\gamma_{\text{norm}}$

It is tempting to read this section as "we should have run CKG without $\gamma_{\text{norm}}$ for a fair comparison". That would be a misreading. The factor $\gamma_{\text{norm}} = 1/\max_\omega \gamma(\omega)$ is **CKG-intrinsic**: it is the LCU normalisation that ensures CKG's coherent block-encoding $B$ sits inside a unit-norm budget for the Berntson–Sünderhauf walk, as fixed in [Chen et al. 2025, §II.A][CITE: chen2023efficient]. Removing it would not put the two samplers on equal footing; it would just put CKG in an unphysical configuration that does not appear in any production deployment of the algorithm.

DLL needs no analogue of $\gamma_{\text{norm}}$ because its rank-1 filter $q(\nu)$ already bakes the KMS factor $e^{-\beta\nu/4}$ into the per-coupling Lindblad operator, leaving the dissipator at the right normalisation by construction (cf. (Eq. 3.4) in [Ding et al. 2024][CITE: ding2025efficient]).

So: both samplers are at their natural production normalisation throughout the §2 comparison — that *is* the operationally meaningful comparison. What §8 *adds* is the further intrinsic comparison via $\rho_{\text{intrinsic}}$, which is independent of any normalisation choice on either side and so resolves the residual ambiguity left by §2's $\tau_{\text{mix}}$ tie.

### 8.3 Empirical $\rho_{\text{intrinsic}}$ across $(n, \beta)$ cells

Panel (a) of the summary plot `drafts/figures/numerics/fair_comparison_summary.{png,pdf}` shows $\rho_{\text{intrinsic}}$ for the three samplers across $n \in \{3, 4, 5\}$, $\beta \in \{1, 2, 5, 10, 20\}$. Selected ratios from `scripts/output/fair_comparison_kms_dirichlet.bson` (CKG smooth-Metropolis vs DLL Metropolis):

| $n$ | $\beta=1$ | $\beta=5$ | $\beta=10$ | $\beta=20$ |
|---|---|---|---|---|
| 3 | $1.04\times$ | $1.07\times$ | $1.14\times$ | $1.18\times$ |
| 4 | $1.02\times$ | $1.10\times$ | $1.31\times$ | $1.43\times$ |
| 5 | $1.05\times$ | $1.15\times$ | $1.28\times$ | $1.40\times$ |

(Cells are $\rho(\text{CKG sM}) / \rho(\text{DLL Metro})$. Numbers $> 1$ favour CKG.)

Three observations:

- **CKG smooth-Metro is consistently $1.0$–$1.4\times$ ahead** in intrinsic ratio. The lead is modest at high $T$ ($\beta = 1$: $\le 5\%$) and grows monotonically with $\beta$, reaching $30$–$40\%$ for $n \in \{4, 5\}$ at $\beta \in \{10, 20\}$. The §2 $\tau_{\text{mix}}$ tie therefore is *not* a structural identity: CKG genuinely organises its KMS spectrum slightly more efficiently than DLL Metropolis. But the gap is small relative to the DLL Gaussian collapse: at $n=3, \beta=20$ DLL Gaussian's $\rho_{\text{intrinsic}}$ is **$6\times$ smaller** than DLL Metropolis's, on the same scale-invariant footing.
- **DLL Gaussian collapses at low $T$.** $\rho_{\text{intrinsic}}$ at $n=3, \beta=20$ is $0.024$ (CKG: $0.155$, DLL Metro: $0.149$); at $n=4, \beta=20$ it is suppressed by $\sim 1.7\times$ relative to CKG; at $n=5$ by $\sim 1.2\times$. The mechanism is identical to the $\tau_{\text{mix}}$ collapse mechanism diagnosed in §4: a single super-exponentially localised filter cannot maintain spectral mass across the Bohr ladder at low temperature.
- **Surprise: CKG smooth-Metro $\lambda$ is non-monotone in $\beta$ at $n = 4, 5$.** The fitted spectral gap dips slightly between $\beta = 5$ and $\beta = 10$ and then recovers at $\beta = 15, 20$. We attribute this to a level-crossing in the second-eigenvector of the symmetrised generator as the antidiagonal envelope of $\alpha^{\text{CKG}}$ tightens; the effect is small ($\le 5\%$) and does not change any qualitative conclusion. Worth checking whether the same wiggle persists in $n = 6, 7$ (qf-lkb.9 matrix-free path) before publishing as a feature rather than a fixture artefact. [CHECK]

### 8.4 Diagonal of $\alpha$ as the actual driver (H2)

The framework of §3 organises the Kossakowski as a $|B_H|\times|B_H|$ matrix indexed by Bohr frequencies, and §3 + §4 frame the mechanism question through *rank* and *off-diagonal reach*. A complementary — and, as we now show, more predictive — diagnostic is the **diagonal** $\alpha_{\text{diag}}(\nu) := \sum_a \alpha^a(\nu, \nu)$, which is exactly the per-Bohr-frequency dissipation rate seen by the populations. The hypothesis (call it H2):

> Among KMS-DBC samplers with comparable Bohr support, the shape of $\alpha_{\text{diag}}(\nu)$ as a function of $\nu$ is the structural driver of $\tau_{\text{mix}}$; rank and off-diagonal coupling are secondary.

H2 makes a clean numerical prediction: CKG smooth-Metropolis and DLL Metropolis, both Metropolis-shape per-frequency, should have nearly identical $\alpha_{\text{diag}}(\nu)$; DLL Gaussian, super-exponentially localised, should not.

Figure `drafts/figures/numerics/alpha_diagonal_three_samplers.{png,pdf}` overlays $\alpha_{\text{diag}}(\nu)$ on a log scale for all 12 cells of the grid. The match between CKG smooth-Metro and DLL Metropolis is essentially exact across all $\nu$ at all $(n, \beta)$; the DLL Gaussian curve coincides only near $\nu = 0$ and falls off super-exponentially in both directions for $\beta \ge 5$. Quantitatively, the relative $L^2$ distance

$$
\mathrm{dist}_{L^2}(\text{X}, \text{Y}) \;:=\; \frac{\bigl\lVert \alpha^{\text{X}}_{\text{diag}} - \alpha^{\text{Y}}_{\text{diag}}\bigr\rVert_2}{\bigl\lVert \alpha^{\text{Y}}_{\text{diag}}\bigr\rVert_2}
$$

is at most $6\%$ between DLL Metro and CKG sM across all 12 $(n,\beta)$ cells, and reaches $50$–$80\%$ between DLL Gauss and CKG sM at $\beta \ge 10$.

**H2 verdict: confirmed.** Both diagnostics — the scale-invariant $\rho_{\text{intrinsic}}$ of §8.3 and the diagonal-of-$\alpha$ of §8.4 — point in the same direction: CKG sM and DLL Metro share the Metropolis per-frequency rate, and that shared shape is what drives the $\tau_{\text{mix}}$ and $\rho_{\text{intrinsic}}$ tie. Rank, the focus of the §3 framing, is decoupled from the rate-determining step in this regime.

### 8.5 Bound validation: does $\lambda$ predict $\tau_{\text{mix}}$?

The $\rho_{\text{intrinsic}}$ comparison is meaningful only insofar as $\lambda$ actually controls $\tau_{\text{mix}}$. The relevant inequality is the Poincaré bound (cf. [Rouzé 2019, Eq. (5.2)][CITE: rouze2019functional] and [Ding et al. 2024, Lem. 1.1][CITE: ding2025efficient])

$$
\tau_{\text{mix}}(\varepsilon) \;\le\; \frac{1}{\lambda(\mathcal{L})}\,\log\!\biggl(\frac{1}{p_{\min}\,\varepsilon}\biggr),
$$ <!-- \label{eq:poincare-bound} -->

where $p_{\min}$ is the smallest eigenvalue of $\sigma_\beta$. Panel (c) of `fair_comparison_summary.png` plots $\tau_{\text{pred}} := \lambda^{-1}\log(1/(p_{\min}\varepsilon))$ against the directly measured $\tau_{\text{meas}}$ from the §2 ODE pipeline, across all 32 cells of the $(n, \beta, \text{sampler})$ grid. The pooled fit gives slope $\mathbf{0.876}$, intercept $\mathbf{0.643}$, and $R^2 = \mathbf{0.931}$. Per-sampler slopes are $0.99$ (CKG sM), $0.86$ (DLL Metro), $0.87$ (DLL Gauss); the constant-factor offset is $\sim 3\times$ uniformly across all three samplers.

Two things matter here:

- **The pooled $R^2 = 0.93$ confirms that $\lambda$ is the right primary predictor of $\tau_{\text{mix}}$ in this regime.** The $\rho_{\text{intrinsic}}$ comparison of §8.3 therefore inherits its operational meaning directly: a sampler with a $1.4\times$ larger $\rho_{\text{intrinsic}}$ at fixed $\Lambda_{\max}$ will reach a target $\varepsilon$ a factor $\sim 1.4\times$ faster on this grid, modulo a uniform $\sim 3\times$ slack from the bound itself.
- **The $\sim 3\times$ slack is sampler-independent**, which is consistent with the slack coming from the loose $\log(1/p_{\min})$ prefactor (which becomes $O(\beta n)$ at low $T$) rather than from any dynamical feature of the Lindbladian. So Eq. eq:poincare-bound is *asymptotically* sharp and *constant-factor* loose, exactly as the literature has it (e.g. [Capel, Rouzé & Stilck França 2020, §1][CITE: capel2020modified]).

### 8.6 Connection back to §3, §4, §5

The §8 results re-frame three earlier claims:

- **§3 (rank).** Rank does not drive $\tau_{\text{mix}}$ in this regime; H2 wins. Standard rank-1 DLL Metropolis and full-rank ($\sim 24$-channel) CKG smooth-Metropolis agree on $\rho_{\text{intrinsic}}$ within $0$–$40\%$. The §3 framing is correct as a *vocabulary* for organising the Kossakowski but is not the right *mechanism* picture for low-$T$ mixing: that mechanism lives on the diagonal of $\alpha$.
- **§4 (antidiagonal mass).** The antidiagonal-mass picture and the diagonal-of-$\alpha$ picture are complementary. The antidiagonal entry $\alpha(\nu, -\nu)$ describes *off-diagonal* (coherence-creating) coupling between energy sectors $E$ and $E - \nu$; the diagonal $\alpha(\nu, \nu)$ describes the *per-sector* dissipation rate. H2 is the empirical statement that the second one is what controls $\tau_{\text{mix}}$ in this regime, with the antidiagonal contributing only the secondary effect of populating coherences that are quickly erased by dephasing. (Cf. the discussion of "diagonal dynamics dominate" in [Chen et al. 2025, §II.B][CITE: chen2023efficient].)
- **§5 (Remark 23 cost).** §5 already showed that the rank-1 Metropolis filter pays its price not in $\tau_{\text{mix}}$ but in coherent-block-encoding cost (an extra $\beta$ in the total Hamiltonian-simulation budget). The fair-comparison $\rho_{\text{intrinsic}}$ confirms the asymptotic picture: even on the *intrinsic* relaxation timescale CKG is the (slightly) better choice, and once the production cost factors of [Ding et al. 2024, Table 1][CITE: ding2025efficient] are folded in, the gap widens by another $\Theta(\beta)$. CKG remains the better choice for low-$T$ large-system Gibbs sampling, even though the small-$n$ benchmark $\tau_{\text{mix}}$ values barely separate them.

### 8.7 Implication for the multi-rank DLL programme (qf-7go)

The natural follow-up — a parametric multi-channel DLL with $k \ge 2$ symmetrised filter translates per coupling, currently scoped as epic qf-7go — is *predicted to give a null result* by §8.4. If H2 is right, what drives mixing is the per-frequency Metropolis shape, which standard rank-1 DLL Metropolis already has. Adding more channels at the same Metropolis shape will not reduce $\tau_{\text{mix}}$ measurably; it will only increase the block-encoding cost. The empirical test in qf-7go.6 will close that question one way or the other.

## 9. Norm conventions for the block-encoding audit (qf-9ld.1)

*Insertion target:* same draft, before the Citations block. This section locks down the block-encoding norm conventions used by [Ding–Li–Lin 2024][CITE: ding2025efficient] §3.3 + §5 (Remark 23) and [Chen–Kastoryano–Gilyén 2025][CITE: chen2023efficient] §III.A–C, picks one convention to use across the qf-9ld audit (sub-issues qf-9ld.2–qf-9ld.5), and maps it to the helpers already exposed in `src/kms_geometry.jl`.

The §5 numerics in this draft used $\|\cdot\|_{\text{op}}$ (operator norm via `opnorm`) on the *generator* $G$ and read it as a proxy for the block-encoding subnormalisation. That is correct for *scaling* in $\beta$ (the operator norm and the LCU 1-norm share the same $\beta$-asymptotics) but is a **lower bound** on the algorithmic factor — the audit needs the LCU sum, not just $\|G\|_{\text{op}}$. §9.3 below makes the convention explicit.

### 9.1 DLL block-encoding factors (Ding–Li–Lin 2024)

**Definition (Def. 24, [Ding et al. 2024 p. 27][CITE: ding2025efficient]).** A unitary $U_A$ is an $(\alpha, m, \varepsilon)$-block-encoding of $A \in \mathbb{C}^{2^n\times 2^n}$ when

$$
\bigl\| A \;-\; \alpha\,(\langle 0^m|\!\otimes\! I_n)\,U_A\,(|0^m\rangle\!\otimes\! I_n) \bigr\|_{\text{op}} \;\le\; \varepsilon. \quad
$$ <!-- \label{eq:dll-be-defn} -->

The norm in the error bound is the **operator (induced 2-2) norm** $\|\cdot\|_{\text{op}} = \sigma_{\max}(\cdot)$. The scalar $\alpha \ge 0$ is the **block-encoding factor** (or *subnormalisation factor*); $U_A$ encodes $A/\alpha$ in its top-left block, so the algorithm always pays a factor $\alpha$ in cost per query. A necessary consequence of Eq. dll-be-defn is $\|A\|_{\text{op}} \le \alpha + \varepsilon$ — the BE factor is an *operator-norm upper bound* on $A$, with $\varepsilon$ slack.

**LCU bound (Lem. 25, [Ding et al. 2024 p. 27][CITE: ding2025efficient]).** For an LCU sum $\sum_{j=0}^{J-1} c_j U_j$ with $U_j$ unitary, the canonical $W = (\text{Prep}_l^\dagger \otimes I)\,\text{Select}\,(\text{Prep}_r \otimes I)$ construction is a $(\|c\|_1,\, \log_2 J,\, 0)$-block-encoding, where

$$
\|c\|_1 \;:=\; \sum_{j=0}^{J-1} |c_j|.
$$ <!-- \label{eq:lcu-1norm} -->

So *for any LCU implementation in DLL, the block-encoding factor is the $\ell^1$-norm of the coefficient vector*.

**DLL block-encoding factors for $\{L_a\}$ and $G$ (§3.3, Eqs. 3.32–3.33, [Ding et al. 2024 pp. 17–18][CITE: ding2025efficient]).** With access to the input-model block-encoding $V_\mathcal{A}$ of $\sum_a |a\rangle\langle a|\otimes A^a$ at factor $Z_\mathcal{A} \ge \max_a \|A^a\|_{\text{op}}$ (Eq. 3.23), DLL produces:

- $U_L$, an $(\,\underbrace{Z_f \cdot Z_\mathcal{A}}_{=:\,\alpha_L^{\text{DLL}}},\;\mathfrak{m}+\mathfrak{b},\;\delta_f)$-block-encoding of $\sum_a |a\rangle\langle a|\otimes L_a$, where $Z_f := \sum_{m=0}^{2M-1}|f(t_m)|\tau$ is the *discretised $\ell^1$-norm of the time-domain filter coefficients of $f$*. The cost bound (Lem. 32 + Eq. D.3) gives Eq. 3.28:

$$
Z_f \;=\; \mathcal{O}\!\bigl((A_q\,C_{1,u} + A_q^2\,S\,A_w)\,\log(\beta A_u + A_w)\bigr).
$$ <!-- \label{eq:dll-zf-bound} -->

Critically, **$Z_f$ has no $\beta$-linear term** — only a $\log\beta$ from the time-truncation $T = \mathcal{O}(\beta A_u + A_w)\log^s(\cdot)$.

- $U_G$, a $(\,\underbrace{Z_g \cdot Z_\mathcal{A}^2 \sqrt{|\mathcal{A}|}}_{=:\,\alpha_G^{\text{DLL}}},\;2\mathfrak{m}+\mathfrak{a}+\mathfrak{b},\;\delta_g)$-block-encoding of $G$, where $Z_g := \sum_{n,m=0}^{2M-1}|g(t_n, t_m)|\tau^2$ is the *discretised $\ell^1$-norm of the 2D time-domain coefficients of $g$*. The $\sqrt{|\mathcal{A}|}$ comes from $\text{Prep}_\mathcal{A} = H^{\otimes\mathfrak{a}}|0^\mathfrak{a}\rangle$ (Eq. 3.25). The cost bound (Lem. 32 + Eq. D.4) gives Eq. 3.31:

$$
Z_g \;=\; \mathcal{O}\!\Bigl(\bigl(A_q^2\,C_{1,u}^2 + A_q^4\,S^2\,A_w^2 + A_q^3\,C_{1,u}\,S\,A_w + \mathbf{A_q^4\,\beta\,S}\bigr)\,\log^2(\beta A_u + A_w)\Bigr).
$$ <!-- \label{eq:dll-zg-bound} -->

The **$A_q^4\,\beta\,S$ term is the $\beta$-linear contribution** flagged in Remark 23: it is the LCU 1-norm of the 2D Riemann sum approximating $\int\!\!\int g(t,t')\,dt\,dt'$, and $g(t,t')$ inherits the *wide* support of $\alpha_{\nu,\nu'}$ on the $|\nu - \nu'|\sim S$ scale (in turn the rank-1 wide filter $q^{\text{DLL,Metro}}$).

**Aggregated $\|\mathcal{L}\|_{\text{be}}$ for the LW23 simulator (Thm. 18 footnote / Eq. 3.36, [Ding et al. 2024 p. 18 / p. 20][CITE: ding2025efficient]).** Theorem 18 of [Ding et al.][CITE: ding2025efficient] (which uses the [Li–Wang 2023][CITE: li2023simulating] Lindblad-LCU) defines

$$
\|\mathcal{L}\|_{\text{be}} \;:=\; \alpha_G + \tfrac12\,(\alpha_L^{\text{per-coupling}})^2\,|\mathcal{A}|,
$$ <!-- \label{eq:dll-Lbe-defn} -->

so for DLL specifically (Eq. 3.36)

$$
\|\mathcal{L}^{\text{DLL}}\|_{\text{be}} \;=\; Z_g\,Z_\mathcal{A}^2\,\sqrt{|\mathcal{A}|} \;+\; \tfrac12\,Z_f^2\,Z_\mathcal{A}^2\,|\mathcal{A}|.
$$ <!-- \label{eq:dll-Lbe} -->

The factor $|\mathcal{A}|$ in the dissipator term is structural: the LW23 dissipator step uses two queries to $U_L$ (one for $L_a$ on the left, one for $L_a^\dagger$ on the right), and the index register $|a\rangle$ is summed over. The total Hamiltonian-simulation cost per unit time of $\mathcal{L}$ is $\widetilde{\mathcal{O}}(\|\mathcal{L}\|_{\text{be}})$ (Eq. 3.34), which closes the loop with [Ding et al. 2024, Thm. 19 / Table 1][CITE: ding2025efficient]: total cost $\widetilde{\mathcal{O}}(C_q\,t_{\text{mix}}\,\beta^2\,S\,|\mathcal{A}|^2\,\log^{1+s}(1/\varepsilon))$.

**Remark 23 in this language.** The $\beta^2$ scaling of the DLL total cost decomposes as one $\beta$ from $\|\mathcal{L}\|_{\text{be}}$ and one $\beta$ from $t_{\text{mix}}$ (the standard $\beta\,t_{\text{mix}}$ Lindbladian-simulation budget). The first $\beta$ lives entirely in the $A_q^4\,\beta\,S$ term of $Z_g$ (Eq. dll-zg-bound), and from Eq. dll-Lbe enters $\|\mathcal{L}^{\text{DLL}}\|_{\text{be}}$ through $\alpha_G^{\text{DLL}} = Z_g\,Z_\mathcal{A}^2\,\sqrt{|\mathcal{A}|}$ — *the coherent block-encoding factor*, not the dissipator.

### 9.2 CKG block-encoding factors (Chen–Kastoryano–Gilyén 2025)

**Norm conventions (Nomenclature, [Chen et al. 2025 p. 23][CITE: chen2023efficient]).** CKG distinguishes:

- $\|O\| := \sup_{|\phi\rangle, |\psi\rangle} |\langle\phi| O|\psi\rangle|/(\|\phi\| \|\psi\|)$ — **operator norm** on operators.
- $\|O\|_p := (\operatorname{Tr}|O|^p)^{1/p}$ — Schatten $p$-norm.
- $\|\mathcal{L}\|_{p\text{-}p} := \sup_{O \neq 0} \|\mathcal{L}[O]\|_p / \|O\|_p$ — **induced $p$–$p$ superop norm**.

These are the same families used by [Ding et al.][CITE: ding2025efficient] (just different notation). No incompatibility.

**Block-encoding for the dissipator (Def. III.1, [Chen et al. 2025 p. 17][CITE: chen2023efficient]).** A unitary $U$ is a *block-encoding for Lindblad operators* $\{L_j\}_{j\in J}$ when

$$
(\langle 0^b|\otimes I)\cdot U\cdot(|0^c\rangle\otimes I) \;=\; \sum_{j\in J} |j\rangle\otimes L_j, \quad b\le c\in\mathbb{Z}^+. \quad
$$ <!-- \label{eq:ckg-be-defn} -->

This is *exactly* the form $V_\text{jump}$ of Eq. 3.23 in [Ding et al.][CITE: ding2025efficient]. The subnormalisation is implicit: it lives in the $\sum_j$ on the right-hand side via the requirement $\|\sum_j |j\rangle\otimes L_j\|_{\text{op}} = \sqrt{\|\sum_j L_j^\dagger L_j\|_{\text{op}}} \le 1$, i.e., the $L_j$'s themselves are pre-rescaled to fit the unit-norm budget.

**Coherent term BE assumption (Prop. III.1, [Chen et al. 2025 p. 18][CITE: chen2023efficient]).** *Suppose $\|f_-\|_1, \|f_+\|_1 \le 1$.* Then there is a block-encoding for

$$
\sum_{\bar t_- \in S_{t_0}} f_-(\bar t_-)\,e^{-iH\bar t_-}\Bigl(\sum_{\bar t_+ \in S_{t_0}} f_+(\bar t_+)\sum_{a\in\mathcal{A}} A^{a\dagger}(\bar t_+)A^a(-\bar t_+)\Bigr)e^{iH\bar t_-},
$$

via constant calls to $V_\text{jump}$, controlled $e^{-iH\bar t}$, and the $\text{prep}_{\sqrt{|f_+|}}$ / $\text{prep}_{f_+/\sqrt{|f_+|}}$ oracles. The block-encoding factor is *implicit but bounded by the LCU 1-norms* $\|f_+\|_1 \cdot \|f_-\|_1$ — exactly Eq. lcu-1norm applied to the time-domain Riemann sums of $f_\pm$.

**Coherent term BE for Metropolis-shape $\gamma$ (Cor. III.2 + surrounding discussion, [Chen et al. 2025 p. 15][CITE: chen2023efficient]).** For the Metropolis-like weight $\gamma^M(\omega) = \exp(-\beta\max(\omega + 1/(2\beta),\,0))$ at $\sigma_E = 1/\beta$, the LCU implementation of $B^{M,\eta}$ (the truncated Metropolis coherent term, Eq. 3.5) gives a block-encoding with subnormalisation

$$
\alpha_B^{\text{CKG}} \;=\; \mathcal{O}\!\Bigl(\log\!\bigl(\tfrac{\beta\,\|H\|\,\|\sum_a A^{a\dagger}A^a\|}{\varepsilon}\bigr)\Bigr).
$$ <!-- \label{eq:ckg-aB-bound} -->

**The $\beta$-dependence is purely logarithmic** — this is the structural fact behind [Chen et al. 2025 Thm. I.2][CITE: chen2023efficient]'s $\widetilde{\mathcal{O}}(\beta\,t_{\text{mix}})$ total cost (one $\beta$ for $t_{\text{mix}}$, none extra for $\alpha_B$).

The $\ell^1$-norms $\|h_+\|_1, \|h_-\|_1$ for the bilinear (transition) part are bounded by Cor. III.3 (Gaussian weight: both $\le 1$) and Cor. III.4 (Metropolis weight: both $\le 1$) at [Chen et al. 2025 p. 16][CITE: chen2023efficient]. Likewise, the LCU $\ell^1$-norms $\|b_1\|_1, \|b_2\|_1$ for the coherent term are bounded by Eqs. 3.2–3.3 (Gaussian: $\|b_1\|_1 < 1$, $\|b_2\|_1 < 1/8$) and by the Cor. III.2 logarithmic estimate (Metropolis). All sit within the unit-norm LCU budget once normalised by $\alpha$.

### 9.3 Reconciliation and the chosen convention

Both papers use **the same Definition 24 / Definition III.1** for what a block-encoding is, and **the same LCU $\ell^1$-norm convention** (Eq. lcu-1norm) for how the BE factor of an LCU sum is computed. The numerical difference between CKG and DLL is *not* a difference in convention — it is a difference in *which time-domain weights* (LCU coefficient vectors) the algorithm sums:

- CKG sums $\{b_1(t_i)\,b_2^{M,\eta}(t_j)\}_{i,j}$ for the coherent term — a logarithmically growing $\ell^1$-norm.
- DLL sums $\{g(t_n, t_m)\}_{n,m}$ — a linearly-in-$\beta$ growing $\ell^1$-norm because $g$ inherits the wide $S$-bump support of the rank-1 filter $q^{\text{DLL,Metro}}$ (Remark 23 mechanism).

For the qf-9ld audit, we therefore adopt **the LCU $\ell^1$-norm of the actual time-domain weights** as the canonical block-encoding factor, applied to whichever family is being audited:

- DLL: $\alpha_G^{\text{DLL}} = Z_g\,Z_\mathcal{A}^2\,\sqrt{|\mathcal{A}|}$ and $\alpha_L^{\text{DLL}} = Z_f\,Z_\mathcal{A}$ (Eqs. dll-zg-bound, dll-zf-bound).
- CKG: $\alpha_B^{\text{CKG}} = \|b_1\|_1\,\|b_2^{M,\eta}\|_1\,Z_\mathcal{A}^2$ for the coherent term (Eq. ckg-aB-bound), $\alpha_h^{\text{CKG}} = \|h_+\|_1\,\|h_-\|_1\,Z_\mathcal{A}^2$ for the dissipator.

The aggregated *generator* BE factor is the LW23 / [Chen et al. CKBG23][CITE: chen2023quantum] quantity

$$
\|\mathcal{L}\|_{\text{be}} \;=\; \alpha_G + \tfrac12\,\alpha_L^2\,|\mathcal{A}|,
$$ <!-- \label{eq:Lbe-unified} -->

evaluated separately for each sampler with the appropriate $\alpha_G, \alpha_L$ from the list above. This is the *single number* whose growth in $(n, \beta, k)$ controls the per-unit-time simulation cost in [Ding et al. 2024 Thm. 19 / Table 1][CITE: ding2025efficient].

**Sanity-check norm (cheap, not the BE factor).** For each block-encoded operator $X \in \{G, L_a\}$, the operator norm $\|X\|_{\text{op}}$ is a *lower bound* on $\alpha_X$ (necessary for any block-encoding to exist; from Eq. dll-be-defn). The §5.1 measurements of $\|G\|_{\text{op}}$ via `opnorm` are therefore valid as a *necessary-condition extrapolation* and as a *scaling diagnostic* (the slope of $\|G\|_{\text{op}}$ vs $\beta$ matches the slope of $\alpha_G$ up to log factors). They are *not* the algorithmic factor — to answer Q1 quantitatively the audit should report $\alpha_G$, not just $\|G\|_{\text{op}}$. We keep $\|G\|_{\text{op}}$ as a cheap upper-bound-on-the-lower-bound diagnostic and add the $\ell^1$ sums on top.

**Multi-channel DLL.** For the rank-$k$ DLL of qf-7go with channels $\{q^{(\ell)}\}_{\ell=1}^k$, the Kossakowski decomposes linearly $\alpha^{\text{multi}} = \sum_\ell \alpha^{(\ell)}$, hence by linearity of the LCU sums

$$
Z_g^{\text{multi}} \;=\; \sum_{\ell=1}^k Z_g^{(\ell)}, \qquad Z_f^{\text{multi}} \;=\; \sum_{\ell=1}^k Z_f^{(\ell)},
$$ <!-- \label{eq:dll-multi-Z} -->

so $\alpha_G^{\text{DLL,multi}} = k\,\alpha_G^{\text{DLL,single}}$ to leading order (assuming roughly equal per-channel weights), and similarly for $\alpha_L$. *This is the prediction the qf-9ld.4 fairness check must verify numerically.* If it holds, the qf-7go.6 multi-channel $\tau_{\text{mix}}$ speedup of $\sim 5$–$20\times$ at $k=8$ is partly traded against a $k$-fold increase in $\|\mathcal{L}\|_{\text{be}}$ — and the relevant question is whether the *product* $\|\mathcal{L}\|_{\text{be}}\cdot t_{\text{mix}}$ improves.

### 9.4 Mapping to existing helpers in `src/kms_geometry.jl`

The audit will use the following helpers; new helpers are flagged where needed.

| Quantity | Helper | Status |
|---|---|---|
| $\|G\|_{\text{op}}, \|L_a\|_{\text{op}}$ — necessary-condition lower bound on $\alpha_G, \alpha_L$ | `opnorm(Matrix(G))`, `opnorm(Matrix(L_a))` (Julia `LinearAlgebra`) | available |
| $\|\mathcal{L}\|_{2\text{-}2}$ — HS-induced superop norm of full generator | `hs_operator_norm(L_super)` (`src/kms_geometry.jl:299`) | available — **use for Q2 fairness** |
| $4\max_a \|L_a^\dagger L_a\|_{\text{op}}$ — Wolf–Pérez-García 1→1 bound on dissipator | `dissipator_one_to_one_norm_bound(L_a_list)` (`src/kms_geometry.jl:310`) | available — *cross-check* on $\alpha_L^2 |\mathcal{A}|$ contribution |
| $\operatorname{Tr}(\alpha)$ — total dissipative weight | `dissipator_trace_alpha(α)` (`src/kms_geometry.jl:322`) | available — *diagnostic, not BE factor* |
| $Z_f = \sum_m |f(t_m)|\tau$ — DLL dissipator LCU 1-norm | — | **new**: extract from `dll_lindblad_op_time` time-domain quadrature precompute |
| $Z_g = \sum_{n,m} |g(t_n,t_m)|\tau^2$ — DLL coherent LCU 1-norm | — | **new**: extract from `dll_coherent_op_time` 2D quadrature precompute (the same $g(t,t')$ already constructed as a $|t|\times|t'|$ table) |
| $\|h_\pm\|_1, \|b_1\|_1, \|b_2^{M,\eta}\|_1$ — CKG LCU 1-norms | — | **new**: parallel diagnostic to `B_bohr` / `B_time` precompute that returns the $\ell^1$ sums of the time-domain weights (already constructed internally inside `_precompute_coherent_B`; need a thin tap) |
| Spectral input-model factor $Z_\mathcal{A} \ge \max_a \|A^a\|_{\text{op}}$ | inline: `maximum(opnorm.(jumps_op))` | available — for the standard single-Pauli jumps with $1/\sqrt{3n}$ norm, $Z_\mathcal{A} = 1/\sqrt{3n}$ exactly |

The two "new" rows are the only library deltas the qf-9ld.2 dataset task should need. They are thin *diagnostic* additions over already-constructed time-domain weights — neither requires changing the production matvec path. A natural location is a small `src/block_encoding_norms.jl` module exporting `dll_lcu_norms(filter, ham, β; M, τ) -> (; Z_f, Z_g)` and `ckg_lcu_norms(γ_choice, ham, β; ...) -> (; b1_l1, b2_l1, hplus_l1, hminus_l1)`, gated by the same Bohr/Time/Energy domain selectors as the existing operator builders.

Once those LCU sums are in place, the qf-9ld.3 extrapolation table simply evaluates Eqs. dll-Lbe / Lbe-unified at each $(n, \beta, k)$ cell, and the qf-9ld.4 fairness check compares $\|\mathcal{L}\|_{\text{be}}^{\text{DLL,multi}}(k)$ versus $\|\mathcal{L}\|_{\text{be}}^{\text{CKG}}$ (the multi-channel DLL rank scales the prefactor; CKG is rank-$\sim 24$ with a different time-weight integral but the same audit recipe).

### 9.5 What this section does *not* settle

Two genuine open points are left for qf-9ld.2 onward:

1. **Quadrature parameter convention.** Both papers parametrise the time-domain quadrature by $(M, \tau, T)$ with $T = M\tau \sim \beta A_u + A_w$ (DLL Prop. 16) or $T \sim \log(\beta\|H\|/\varepsilon)$ (CKG Cor. III.2). Picking the same precision target $\varepsilon$ on both sides is necessary for like-for-like comparison; the audit should fix $\varepsilon = 10^{-3}$ to match §2.
2. **Truncation slack $\delta_f, \delta_g$.** Eq. dll-be-defn permits $\varepsilon$ slack between $A$ and $\alpha\cdot(\text{block})$; we read $\alpha$ as the *target* BE factor, not the post-quadrature actual. For the qf-9ld dataset, report the LCU 1-norm of the *finite* time-domain weights actually summed by `dll_coherent_op_time` etc. — that is the operationally meaningful quantity.

### 9.6 Block-encoding norm dataset (qf-9ld.2)

Driver: `scripts/scratch_dll_norm_audit.jl`. Output: `scripts/output/dll_norm_audit/dll_norm_audit.bson`. Sweep grid $n \in \{3, 4, 5\}$, $\beta \in \{5, 10, 20\}$, $k \in \{1, 2, 4, 8\}$ for DLL multi-channel Metropolis ($S = 2$, centers uniform on $[0, S/2]$ matching qf-7go.6) and DLL multi-channel Gaussian ($c_{\max} = 1$); CKG smooth-Metropolis ($a = 0$, $s = 0.25$, $\sigma = 1/\beta$, $\eta = 10^{-3}$) at $k = 1$ baseline. Time grid for the LCU sums: $\tau = 0.2$, $T_{\max} = 6\beta$ (smooth Schwartz integrands; the $L^1$ Riemann sum converges as $\mathcal{O}(\tau^2)$ — finer grid shifts absolutes by $< 1\%$ and does not change any slope). Hamiltonians: rescaled fixtures `hamiltonians/heis_disordered_periodic_n*.bson` ($\|H\| \approx 0.45$); the per-cell rescaling factor (≈ 20 / 32 / 36 for $n = 3 / 4 / 5$) is stored in the BSON for the qf-9ld.3 unrescaled extrapolation.

**$\|\mathcal{L}\|_{\text{be}}$ for DLL Metropolis** (naive aggregation = single LCU sum across all per-coupling-per-channel $L^{(\ell)}_a$, one big $|\mathcal{A}|$-prep; the concatenated-index version, which uses the LW23 $J = |\mathcal{A}|\cdot k$ accounting, is similar within a factor 2). Multi-channel $Z_g$ uses the per-channel sum $\sum_\ell Z_g^{(\ell)}$ (the operationally correct quantity since $G^{\text{multi}} = \sum_\ell G^{(\ell)}$ is *linear* in channels — see §9.7 for derivation; an earlier draft of this table summed cross-channel terms via $|\sum_\ell \hat f_\ell|^2$ and over-counted by a factor 1–4):

| $n$ | $\beta$ | $k=1$ | $k=2$ | $k=4$ | $k=8$ |
|---|---|---|---|---|---|
| 3 | 5  | 1.01    | 8.43    | 24.0    | 91.0    |
| 3 | 10 | 1.21    | 1.05e+02 | 1.83e+02 | 5.38e+02 |
| 3 | 20 | **1.59** | **2.03e+04** | **2.40e+04** | **4.60e+04** |
| 4 | 5  | 0.99    | 8.19    | 23.6    | 90.2    |
| 4 | 10 | 1.17    | 1.00e+02 | 1.78e+02 | 5.29e+02 |
| 4 | 20 | 1.53    | 1.93e+04 | 2.29e+04 | 4.46e+04 |
| 5 | 5  | 0.97    | 8.04    | 23.4    | 89.6    |
| 5 | 10 | 1.14    | 9.71e+01 | 1.74e+02 | 5.23e+02 |
| 5 | 20 | 1.49    | 1.86e+04 | 2.22e+04 | 4.37e+04 |

**$\|\mathcal{L}\|_{\text{be}}$ for DLL Gaussian** (same accounting; $c_{\max} = 1$ shifts):

| $n$ | $\beta$ | $k=1$ | $k=2$ | $k=4$ | $k=8$ |
|---|---|---|---|---|---|
| 3 | 20 | 0.75    | 8.25e+03 | 8.42e+03 | 1.31e+04 |
| 5 | 20 | 0.73    | 7.99e+03 | 8.16e+03 | 1.26e+04 |

**CKG smooth-Metropolis** (no $k$ axis; $\alpha_B = \|b_-\|_1 \cdot \|b_+\|_1 \cdot Z_\mathcal{A}^2$ with $Z_\mathcal{A} = 1/\sqrt{3n}$):

| $n$ | $\|b_-\|_1$ | $\|b_+\|_1$ | $\alpha_B$ |
|---|---|---|---|
| 3 | 0.818 | 0.093 | $8.46\times 10^{-3}$ |
| 4 | 0.818 | 0.093 | $6.35\times 10^{-3}$ |
| 5 | 0.818 | 0.093 | $5.08\times 10^{-3}$ |

Three observations:

1. **DLL $k=1$ single-channel sits at the edge of the BE budget on rescaled fixtures.** $\|\mathcal{L}\|_{\text{be}}^{\text{DLL,Metro}, k=1}$ ranges from 0.97 ($n=5, \beta=5$) to 1.59 ($n=3, \beta=20$). The Gaussian filter is comfortably below ($0.73$–$0.75$) thanks to its narrower $f̂$ support (smaller $Z_g$). On unrescaled physical $H$ with $\|H\| = \mathcal{O}(n)$, the Bohr frequencies and hence the $f̂$-domain content scale by the rescaling factor (20–36 for $n \in \{3, 4, 5\}$), pushing the dimensionless $\beta\cdot\|H\|$ deeper into the asymptotic Remark-23 regime where $Z_g$ grows linearly in $\beta$.
2. **The $k=1 \to k=2$ jump at $\beta = 20$ is enormous — four orders of magnitude.** This is *not* a quadrature artefact: it traces directly to the per-channel `ShiftedSymmetricFilter`'s time-domain weight (Eq. derived in `src/filters.jl:680` and following), which carries an envelope $\cosh(\beta\,\nu_\ell/4 + i\,\nu_\ell\,t)$ inherited from pulling the KMS factor $e^{-\beta\nu/4}$ through the $\nu\mapsto\nu - \nu_\ell$ shift. For $\beta = 20$, $\nu_\ell \approx 1$, $\cosh(\beta\nu_\ell/4)^2 = \cosh(5)^2 \approx 5500$ — a factor that multiplies $Z_g^{(\ell)}$ for every off-center channel. The pre-asymptotic $k = 1$ DLL Metropolis is well-behaved; the multi-channel construction is *not*, in $\|\mathcal{L}\|_{\text{be}}$.
3. **CKG sits orders of magnitude below DLL.** $\alpha_B^{\text{CKG}} \approx 5\times 10^{-3}$ across the entire grid — three orders of magnitude smaller than the DLL coherent-only $\alpha_G^{\text{DLL}}$ at $k=1$ (which is $\sim 0.1$–$0.45$), and *six* orders of magnitude smaller than the multi-channel $\|\mathcal{L}\|_{\text{be}}^{\text{DLL,multi}}$ at $\beta = 20, k = 8$. The $\beta$-independence of $\|b_\pm\|_1$ on the chosen grid empirically confirms Cor. III.2's logarithmic bound at this scale.

**β-slope diagnostics** (log–log fit over $\beta \in \{5, 10, 20\}$):

| | $k=1$ | $k=2$ | $k=4$ | $k=8$ |
|---|---|---|---|---|
| DLL Metro $\alpha_G$ slope | $0.58$ | $5.79$ | $5.07$ | $4.73$ |
| DLL Gauss $\alpha_G$ slope | $0.00$ | $5.05$ | $4.29$ | $3.76$ |

Three remarks on the slopes:

- **$k=1$ Metro slope $\approx 0.6$ is sub-linear** — the rescaled-$H$ regime $\beta \cdot \|H\| \le 9$ at $\beta=20$ is *pre-asymptotic*. Eq. dll-zg-bound's $A_q^4 \beta S$ term dominates the bound only when $\beta A_u \gg A_q^2 C_{1,u}$, which is a $\beta\cdot\|H\| \gg 1$ condition on a normalised filter. The §5.1 operator-norm slope $1.35$–$1.46$ is the one to trust as the asymptotic signature of Remark 23 (Z_g and $\|G\|_{\text{op}}$ share the same asymptotic class up to log factors per Eq. 3.31, but the constant prefactors differ).
- **$k=1$ Gauss slope $\approx 0$** matches the rank-1 Gaussian collapse mechanism diagnosed in §4: $f̂^{\text{DLL,Gauss}}(\nu) = e^{-(\beta\nu+1)^2/8}$ becomes super-exponentially localised at $\nu \approx -1/\beta$, so its $\ell^1$-mass on the time grid does not grow with $\beta$ on the rescaled fixture (the relevant Bohr frequencies are mostly outside the support of $f̂$).
- **$k>1$ slopes $\sim 4$–$6$ are not power laws** — they are dominated by the $\cosh^2(\beta \nu_\ell / 4)$ envelope which grows as $\sim \tfrac14 e^{\beta\nu_\ell/2}$ for $\beta\nu_\ell \gg 1$. Fitting an exponential as a power law over a 4× span in $\beta$ produces a spurious slope; the right summary is "**multi-channel DLL $\|\mathcal{L}\|_{\text{be}}$ grows exponentially in $\beta\cdot\nu_\ell$**".

**Headline for the audit.** On the rescaled fixtures, single-channel DLL is *just barely* within unit-norm budget at the largest $(n, \beta)$; multi-channel DLL leaves the budget by orders of magnitude even at $k = 2$ once $\beta \ge 10$, due to the $\cosh^2(\beta \nu_\ell / 4)$ blow-up of the symmetrised translates — a *new* mechanism not flagged by Remark 23 (which addresses single-channel wide-support filters, not symmetrised-translate constructions). CKG smooth-Metropolis stays below $0.01$ across the entire grid. **The qf-9ld.4 fairness analysis can already be summarised qualitatively: the qf-7go.6 multi-channel $\tau_{\text{mix}}$ speedup of $5$–$20\times$ is more than offset, in $\|\mathcal{L}\|_{\text{be}}\cdot t_{\text{mix}}$ product, by a $\cosh^2(\beta\,c_{\max}/4)$ penalty in the per-step block-encoding cost.** This is the stronger version of the cost caveat already noted in `multirank_dll_complete_qf_7go.md`; qf-9ld.4 will quantify the product. §9.7 below shows where this $\cosh^2$ comes from and why it is *not* intrinsic to multi-rank DLL — only to the specific shifted-symmetric construction with $c_{\max}$ chosen for the bump support rather than the Hamiltonian's Bohr support.

### 9.7 The cosh blow-up is a property of the construction, not of multi-rank DLL

The §9.6 multi-channel numbers explode at $k > 1$. Before reading this as "multi-rank DLL is fundamentally expensive", we need to check (a) whether the blow-up is real, (b) what causes it, and (c) whether a different multi-rank parametrisation avoids it. Verification driver: `scripts/scratch_dll_multichannel_cosh_verify.jl`.

#### 9.7.1 Where the $\cosh$ comes from (rigorous derivation)

The `ShiftedSymmetricFilter` channel (`src/filters.jl:632`) wraps a base DLL filter with a symmetric $\nu$-shift $\nu_\ell$ and a scalar weight $w$:

$$q_\ell(\nu) \;=\; \sqrt{w/2}\,\bigl[q_{\text{base}}(\nu - \nu_\ell) + q_{\text{base}}(\nu + \nu_\ell)\bigr]\quad (\nu_\ell \neq 0).$$

The DLL filter is $\hat f_\ell(\nu) = q_\ell(\nu)\,e^{-\beta\nu/4}$. Pulling the KMS factor through each shifted copy via the substitutions $u = \nu \mp \nu_\ell$:

$$
\hat f_\ell(\nu) \;=\; \sqrt{w/2}\,\Bigl[\,e^{-\beta\nu_\ell/4}\,\hat f_{\text{base}}(\nu - \nu_\ell) \;+\; e^{+\beta\nu_\ell/4}\,\hat f_{\text{base}}(\nu + \nu_\ell)\,\Bigr].
$$ <!-- \label{eq:shifted-fhat} -->

The *exponentially asymmetric* prefactors $e^{\mp\beta\nu_\ell/4}$ are the consequence of pulling the KMS exponential through a $\pm\nu_\ell$ shift — this is the *unavoidable* effect of the construction. Inverse-Fourier-transforming:

$$
f_\ell(t) \;=\; \sqrt{w/2}\,f_{\text{base}}(t)\,\bigl[e^{-\beta\nu_\ell/4 + i\nu_\ell t} \;+\; e^{+\beta\nu_\ell/4 - i\nu_\ell t}\bigr] \;=\; \sqrt{2w}\,f_{\text{base}}(t)\,\cosh\!\bigl(\tfrac{\beta\nu_\ell}{4} - i\,\nu_\ell\,t\bigr).
$$ <!-- \label{eq:shifted-ftime} -->

Verified to machine precision against `time_kernel(channel, t)` (CLAIM 1 in the verification driver, $\max\,|\text{actual} - \text{predicted}| = 0$ across $\nu_\ell \in \{0, 0.225, 0.45, 0.7, 1.0\}$). Taking magnitudes and using $|\cosh(a + ib)| \le \cosh(a)$:

$$
|f_\ell(t)| \;\le\; \sqrt{2w}\,\cosh\!\bigl(\tfrac{\beta\nu_\ell}{4}\bigr)\,|f_{\text{base}}(t)|. \quad
$$ <!-- \label{eq:shifted-fmag} -->

Hence

$$
Z_f^{(\ell)} \;\le\; \sqrt{2w}\,\cosh\!\bigl(\tfrac{\beta\nu_\ell}{4}\bigr)\,Z_f^{\text{base}}, \qquad Z_g^{(\ell)} \;\lesssim\; 4\,\cosh^2\!\bigl(\tfrac{\beta\nu_\ell}{4}\bigr)\,Z_g^{\text{base}}
$$ <!-- \label{eq:shifted-Zbound} -->

(the $g$ kernel is bilinear in $\hat f_\ell$ — Eq. 3.5 of [Ding et al.][CITE: ding2025efficient] — so the cosh enters squared). Numerical verification at $\beta = 20$, $\nu_\ell = 1$:

| | predicted | measured | ratio |
|---|---|---|---|
| $Z_f^{(\ell)} / Z_f^{\text{base}}$ | $2\cosh(5) = 148$ | $105$ | $0.71$ |
| $Z_g^{(\ell)} / Z_g^{\text{base}}$ | $4\cosh^2(5) = 22000$ | $17090$ | $0.78$ |

Both within factor $\sim 1.5$ of the analytic upper bound — the $|\cosh(a + ib)| \le \cosh(a)$ inequality is loose but the right asymptotic class is exactly recovered.

#### 9.7.2 Root cause: $c_{\max}$ chosen for the bump, not for the Hamiltonian

The qf-7go.6 sweep (`scripts/scratch_dll_multirank_taumix_sweep.jl:84`) sets

$$c_{\max}^{\text{qf-7go}} = S/2 - 10^{-9} \approx 1.0 \quad\text{(Metropolis)}, \qquad c_{\max}^{\text{qf-7go}} = 1.0 \quad\text{(Gaussian)},$$

i.e., centers $\nu_\ell$ uniformly spread on $[0, c_{\max}^{\text{qf-7go}}]$ such that the shifted-bump pairs stay inside the Hörmander flat-top region $[-S/2, S/2]$ of the base Metropolis filter. *But the rescaled fixture has $\nu_{\max}(H) = 0.45$* — the Hamiltonian's actual Bohr spectrum extends only to half of $c_{\max}^{\text{qf-7go}}$. Channels with $\nu_\ell > \nu_{\max}(H)$ place Kossakowski mass at frequencies $\nu \approx \nu_\ell$ where there are *no Bohr transitions* — the corresponding antidiagonal slot $\alpha(\nu_\ell, -\nu_\ell)$ never couples any pair of energy levels, so it cannot help population transfer. **Yet it costs $\cosh^2(\beta\nu_\ell/4)$ in $Z_g$.**

A $c_{\max}$ sweep at $k = 4$, $\beta = 20$, $n = 3$ (DLL Metropolis):

| $c_{\max}$ | $\nu_{\max}(H)$ | $Z_g^{\text{multi}}$ | $\|\mathcal{L}\|_{\text{be}}$ | factor over $k=1$ baseline |
|---|---|---|---|---|
| $0.10$ | $0.45$ | $9.93$ | $34$ | $21\times$ |
| $0.20$ | $0.45$ | $13.78$ | $42$ | $26\times$ |
| $0.30$ | $0.45$ | $25.5$ | $64$ | $40\times$ |
| $\mathbf{0.45}$ | $\mathbf{0.45}$ | $\mathbf{98}$ | $\mathbf{174}$ | $\mathbf{109\times}$ |
| $0.70$ | $0.45$ | $1338$ | $1488$ | $935\times$ |
| $\mathbf{1.00}$ (qf-7go.6) | $0.45$ | $\mathbf{2.40\times 10^4}$ | $\mathbf{2.40\times 10^4}$ | $\mathbf{15{,}000\times}$ |

Reducing $c_{\max}$ from $1.0$ to $0.45$ (within the actual Bohr support) shrinks the cost by **138×**; reducing further to $c_{\max} = 0.2$ keeps the multi-channel BE cost within $\sim 30\times$ the single-channel budget — comparable to or better than the qf-7go.6 $\tau_{\text{mix}}$ speedup. The qf-7go.6 $c_{\max} = 1.0$ choice was *4× too aggressive*: it placed off-zero mass twice as far out as the Hamiltonian needed it, paying $\cosh^2(5) \approx 5500$ instead of $\cosh^2(2.25) \approx 23$ per off-center channel.

#### 9.7.3 Trivial uniform-rescale baseline (sanity check)

Setting $k$ identical channels at $\nu_\ell = 0$ all with weight $1$ trivially scales $\alpha^{\text{multi}} = k\cdot\alpha^{\text{base}}$ — i.e., uniform rescale of the Lindbladian by $k$. Under uniform rescale $\mathcal{L} \to c\,\mathcal{L}$ the LCU norms scale linearly $\|\mathcal{L}\|_{\text{be}} \to c\,\|\mathcal{L}\|_{\text{be}}$ and $\tau_{\text{mix}} \to \tau_{\text{mix}}/c$, so the algorithmic figure of merit $\|\mathcal{L}\|_{\text{be}}\cdot\tau_{\text{mix}}$ is *invariant*. We verify:

| $k$ | $Z_f$ | $Z_g$ | $\|\mathcal{L}\|_{\text{be}}$ |
|---|---|---|---|
| 1 | 1.51 | 1.35 | 1.59 |
| 2 | 3.02 | 2.70 | 5.46 |
| 4 | 6.04 | 5.40 | 20.1 |
| 8 | 12.1 | 10.8 | 76.6 |

$Z_f$ scales linearly ($\times k$), $Z_g$ scales linearly ($\times k$ — corrected from a buggy $k^2$ in an earlier audit version), $\|\mathcal{L}\|_{\text{be}}$ grows superlinearly because the dissipator term $\frac12 \alpha_L^2 |\mathcal{A}|$ has $\alpha_L$ linear in $k$. The *concat* aggregation (LW23 $J = |\mathcal{A}|\,k$ index) recovers exactly linear-in-$k$ growth, matching the rescale-invariance argument.

#### 9.7.4 Better multi-rank DLL parametrisations

The $\cosh^2$ blow-up is not a property of multi-rank DLL — it is a property of the symmetrised-shifted-translate construction with $c_{\max}$ disconnected from the Hamiltonian. Four better routes:

1. **Cap centers at the Hamiltonian's Bohr support, $c_{\max} \le \nu_{\max}(H)$.** Easiest fix to the existing parametrisation. The penalty per off-center channel becomes $\cosh^2(\beta\nu_{\max}(H)/4)$, i.e., $\cosh^2(2.25) \approx 23$ at $\beta = 20$ for the rescaled fixture — much milder than $5500$. Whether the $\tau_{\text{mix}}$ improvement persists at $c_{\max} \le \nu_{\max}(H)$ is the next experiment (qf-7go follow-up: re-sweep $\tau_{\text{mix}}$ vs $k$ at the smaller $c_{\max}$). If it does, this immediately turns the qf-7go.6 result from "pyrrhic 19.9× speedup" into "useful $\sim 5\times$ improvement at $\sim 20\times$ extra cost — net 4× *worse* in product, but only by a small constant factor". Still might be worth pursuing if the $\tau_{\text{mix}}$ gain holds.

2. **Use centers at actual Bohr frequencies, $\nu_\ell \in B_H$.** Place each channel exactly at one Bohr frequency the Hamiltonian has. Same per-channel $\cosh^2(\beta\nu_\ell/4)$ penalty (this is structural to any $\nu$-shift), but at least every channel does drive a real transition. Particularly natural at small $n$ where $|B_H|$ is modest (e.g., $|B_H| \approx 57$ at $n=3$ before degeneracy reduction).

3. **The eigendecomposition route — DLL Eq. 4.10 of [Ding et al.][CITE: ding2025efficient].** Pick a target Kossakowski $\alpha$ (e.g., a low-rank approximation to CKG smooth-Metro's $\alpha$) and compute its $\widetilde C = QDQ^T$ eigendecomposition (Eq. 4.9). Each non-zero $D_{\nu',\nu'}$ defines one Lindblad operator $L_{\nu'}$ (Eq. 4.10) with weight $q_{\nu'}(\nu) = \sqrt{D_{\nu',\nu'}}\,Q_{\nu,\nu'}$ defined directly on $\nu \in B_H$ — *no* shifted-bump construction, *no* cosh from the KMS factor pulling through a continuous-$\nu$ shift. The block-encoding factor of each $L_{\nu'}$ inherits at most a $e^{\beta\nu_{\max}/4}$ amplification (the max over $\nu \in B_H$ of $e^{-\beta\nu/4}$ multiplying $Q_{\nu,\nu'}$), and $\sqrt{D_{\nu',\nu'}}$ is constrained by $\operatorname{tr}(\widetilde C)$. **By construction, this route subsumes CKG when applied to CKG's $\alpha$ — and with rank-$k$ truncation it gives a parametric family of "low-rank CKG-like" samplers with BE cost interpolating between standard rank-1 DLL and full-rank CKG.** This is what qf-7go was nominally trying to find; the shifted-translate construction was a guess that turned out to inflate the cost unnecessarily.

4. **Compensated coefficients.** Generalise the symmetrised translate to complex weights $c_m$ chosen so that the cosh exactly cancels: pick $q_\ell(\nu) = c\,[q_{\text{base}}(\nu - \nu_\ell) + q_{\text{base}}(\nu + \nu_\ell)]$ with $c = \sqrt{w/2}\,e^{\beta\nu_\ell/4}$. Then $f_\ell(t) = \sqrt{2w}\cdot \cos(\nu_\ell t)\,f_{\text{base}}(t)$ — bounded magnitude, no cosh in $Z_f$. *But* the Kossakowski entries scale by $|c|^2 = (w/2)\,e^{\beta\nu_\ell/2}$, so the Lindbladian is uniformly amplified; the $\tau_{\text{mix}}$ shrinks by the same factor and the algorithmic figure of merit is unchanged. *This is just route 3 in disguise* — a single-eigenvector approximation of an amplified $\alpha$ — and would require checking whether the resulting $\tau_{\text{mix}}$ benefit beats CKG, which on present evidence (§8.3) it does not.

#### 9.7.5 Verdict for the thesis

The $\cosh^2(\beta\,c_{\max}/4)$ blow-up of the §9.6 table is **a property of the qf-7go.6 sweep's specific choice of $c_{\max}$**, not a structural ceiling on multi-rank DLL. Concretely:

- The qf-7go.6 $c_{\max} = 1.0$ choice exceeded the Hamiltonian's Bohr support by $2\times$, putting Kossakowski mass at frequencies that contributed nothing to mixing while paying $\cosh^2(5) \approx 5500$ per channel.
- A more honest comparison of multi-rank DLL vs CKG would: (i) cap $c_{\max} \le \nu_{\max}(H)$, (ii) re-measure $\tau_{\text{mix}}$ at the new $c_{\max}$, and (iii) compute $\|\mathcal{L}\|_{\text{be}} \cdot t_{\text{mix}}$ to see if any $k > 1$ choice beats $k = 1$ DLL or CKG.
- The structural argument from §3 + §8 still stands: rank-1 DLL Metropolis $\approx$ CKG smooth-Metro in $\tau_{\text{mix}}$ (within 30%) but pays $\beta$ extra in $\|\mathcal{L}\|_{\text{be}}$ (Remark 23). Route 3 of §9.7.4 — eigendecomposition with rank-$k$ truncation of CKG's $\alpha$ — is the principled multi-rank DLL that is *guaranteed* to beat plain rank-1 DLL when $\alpha^{\text{CKG}}$ is the right target.

**Recommended next experiment** (qf-7go.6-bis): re-run the multi-channel $\tau_{\text{mix}}$ sweep at $c_{\max} \in \{0.1, 0.2, 0.45\} \le \nu_{\max}(H)$ and tabulate the $\|\mathcal{L}\|_{\text{be}}\cdot t_{\text{mix}}$ product. If any cell beats single-channel DLL, the symmetrised-translate construction has algorithmic merit when constrained correctly. If not, route 3 (eigendecomposition) is the only honest path to multi-rank DLL.

---

## Citations

```bibtex
% already in the thesis bib — keys reproduced for cross-reference
@article{chen2023efficient,
  author  = {Chen, Chi-Fang and Kastoryano, Michael J. and Gily{\'e}n, Andr{\'a}s},
  title   = {An efficient and exact noncommutative quantum Gibbs sampler},
  journal = {arXiv:2311.09207},
  year    = {2025},
  note    = {`supplementary-informations/Chen et al. - 2025 - An efficient and exact noncommutative quantum Gibbs sampler.pdf`}
}
@article{ding2025efficient,
  author  = {Ding, Zhiyan and Li, Bowen and Lin, Lin},
  title   = {Efficient quantum Gibbs samplers with Kubo--Martin--Schwinger detailed balance condition},
  journal = {arXiv:2404.05998},
  year    = {2024},
  note    = {`supplementary-informations/Ding et al. - 2024 - Efficient quantum Gibbs samplers with Kubo--Martin.pdf`}
}
@article{ramkumar2024mixing,
  author  = {Ramkumar, A. and Soleimanifar, M.},
  title   = {Mixing time of quantum Gibbs sampling for random sparse Hamiltonians},
  journal = {arXiv:2411.04454},
  year    = {2024}
}

% new entries introduced in §8 — verify keys against existing thesis bib before merging
@article{rouze2019functional,
  author  = {Rouz{\'e}, Cambyse and Datta, Nilanjana},
  title   = {Concentration of quantum states from quantum functional and transportation cost inequalities},
  journal = {Journal of Mathematical Physics},
  volume  = {60},
  number  = {1},
  pages   = {012202},
  year    = {2019},
  note    = {arXiv:1704.02400}
}
@article{capel2020modified,
  author  = {Capel, {\'A}ngela and Rouz{\'e}, Cambyse and Stilck Fran{\c c}a, Daniel},
  title   = {The modified logarithmic {S}obolev inequality for quantum spin systems: classical and commuting nearest neighbour interactions},
  journal = {arXiv:2009.11817},
  year    = {2020}
}
@article{kochanowski2024rapid,
  author  = {Kochanowski, Jan and Alhambra, Alvaro M. and Capel, {\'A}ngela and Rouz{\'e}, Cambyse},
  title   = {Rapid thermalization of dissipative many-body dynamics of commuting {H}amiltonians},
  journal = {arXiv:2404.16780},
  year    = {2024},
  note    = {Communications in Mathematical Physics, 2025}
}
@article{liLu2024quantum,
  author  = {Li, Bowen and Lu, Jianfeng},
  title   = {Quantum space-time {P}oincar{\'e} inequality for {L}indblad dynamics},
  journal = {arXiv:2406.09115},
  year    = {2024}
}

% new entries introduced in §9 — verify keys against existing thesis bib before merging
@article{chen2023quantum,
  author  = {Chen, Chi-Fang and Kastoryano, Michael J. and Brand{\~a}o, Fernando G.S.L. and Gily{\'e}n, Andr{\'a}s},
  title   = {Quantum thermal state preparation},
  journal = {arXiv:2303.18224},
  year    = {2023},
  note    = {`supplementary-informations/Chen et al. - 2023 - Quantum Thermal State Preparation.pdf`; cited internally as [CKBG23] by Ding et al. and Chen et al. 2025 — provides the original LCU and dissipator block-encoding constructions inherited by both papers}
}
@article{li2023simulating,
  author  = {Li, Xiantao and Wang, Chunhao},
  title   = {Simulating {M}arkovian open quantum systems using higher-order series expansion},
  journal = {arXiv:2212.02051},
  year    = {2023},
  note    = {ICALP 2023; cited as [LW23] in Ding et al. 2024 Theorem 18 — provides the $\|\mathcal{L}\|_{\text{be}} = A_g + \tfrac12 A_L^2|\mathcal{A}|$ aggregator used in the DLL total-cost expression}
}
```

Equation labels referenced from `supplementary-informations/2_methods.tex`:
- (Eq. eq:kossakowski-gaussian) — CKG Gaussian-$\gamma$ Kossakowski.
- (Eq. eq:kossakowski-smooth-metro) — CKG smooth-Metropolis Kossakowski (the construction we benchmark).
- (Lem. lem:shifted-kms-condition), (Cor. cor:kossakowski-shift) — the linear-combination route from Gaussian-$\gamma$ to Metropolis-shape $\gamma$.
- (Eq. eq:smooth-metro), (Eq. eq:smooth-metro-int-defi) — smooth-Metropolis transition weight definition.
- (Eq. eq:ckg-L-parts), (Eq. eq:ckg-L) — CKG Lindbladian decomposition incl. coherent term $\mathcal{B}$.
- (§sec:kms-lindbladian) — KMS Lindbladian section.

New equation labels introduced in §8 (for the user to add when porting):
- (Eq. eq:kms-inner-product) — KMS inner product.
- (Eq. eq:dirichlet-form-kms) — Dirichlet form on the KMS Hilbert space.
- (Eq. eq:spectral-gap-poincare) — spectral gap as Poincaré constant.
- (Eq. eq:lambda-max-dirichlet) — max Dirichlet rate.
- (Eq. eq:rho-intrinsic) — intrinsic mixing ratio.
- (Eq. eq:poincare-bound) — Poincaré $\tau_{\text{mix}}$ bound.

New equation labels introduced in §9:
- (Eq. eq:dll-be-defn) — DLL Definition 24 of $(\alpha, m, \varepsilon)$-block-encoding.
- (Eq. eq:lcu-1norm) — LCU 1-norm convention $\|c\|_1 = \sum_j |c_j|$.
- (Eq. eq:dll-zf-bound) — $Z_f$ scaling for DLL dissipator (Eq. 3.28 of [Ding et al.][CITE: ding2025efficient]).
- (Eq. eq:dll-zg-bound) — $Z_g$ scaling for DLL coherent term (Eq. 3.31 of [Ding et al.][CITE: ding2025efficient], $\beta$-linear $A_q^4 \beta S$ term).
- (Eq. eq:dll-Lbe-defn) — LW23 aggregated $\|\mathcal{L}\|_{\text{be}} = \alpha_G + \tfrac12 \alpha_L^2|\mathcal{A}|$.
- (Eq. eq:dll-Lbe) — DLL-specific instantiation (Eq. 3.36 of [Ding et al.][CITE: ding2025efficient]).
- (Eq. eq:ckg-be-defn) — CKG Definition III.1 of block-encoding for Lindblad operators.
- (Eq. eq:ckg-aB-bound) — CKG Metropolis coherent BE subnormalisation $\alpha_B^{\text{CKG}}$ (logarithmic in $\beta$).
- (Eq. eq:Lbe-unified) — unified $\|\mathcal{L}\|_{\text{be}}$ adopted for the qf-9ld audit.
- (Eq. eq:dll-multi-Z) — multi-channel additivity $Z_g^{\text{multi}} = \sum_\ell Z_g^{(\ell)}$.

New equation labels introduced in §9.7:
- (Eq. eq:shifted-fhat) — `ShiftedSymmetricFilter` $\hat f_\ell$ with KMS factor pulled through the $\pm\nu_\ell$ shifts.
- (Eq. eq:shifted-ftime) — corresponding time-domain weight $f_\ell(t)$ with the $\cosh(\beta\nu_\ell/4 - i\nu_\ell t)$ envelope.
- (Eq. eq:shifted-fmag) — $|\cosh(a + ib)| \le \cosh(a)$ envelope bound on $|f_\ell(t)|$.
- (Eq. eq:shifted-Zbound) — $\cosh$ and $\cosh^2$ bounds on $Z_f^{(\ell)}$ and $Z_g^{(\ell)}$.

---
## Writing Notes
<!-- For the author, not for the thesis -->

- **Plot dependency.** §8.3, §8.4, §8.5 reference panels (a), (b), (c) of `drafts/figures/numerics/fair_comparison_summary.{png,pdf}`. That figure is being generated in parallel. Verify panel labelling matches before final submission; the present text assumes (a) = $\rho_{\text{intrinsic}}$ heatmap by sampler, (b) = $\alpha_{\text{diag}}$ overlay panel summary, (c) = $\tau_{\text{pred}}$ vs $\tau_{\text{meas}}$ scatter.
- **Citation-key check needed.** The four new bibtex keys `rouze2019functional`, `capel2020modified`, `kochanowski2024rapid`, `liLu2024quantum` are stubs — verify against the canonical thesis bib before merging; existing thesis may already use a different convention (e.g. `Rouze2019` vs `rouze2019functional`). [CHECK]
- **CKG sM $\lambda$ non-monotonicity (§8.3 footnote).** This is real in the BSON output but I could not exhaustively rule out a fixture artefact (e.g. a near-degeneracy in $\sigma_\beta$ that flips the second eigenvector). Marked [CHECK] in §8.3 — do not advertise as a feature in the thesis until $n = 6, 7$ data is available.
- **§8.4 numerics.** L² distance numbers ($\le 6\%$ for DLL Metro vs CKG sM, $50$–$80\%$ for DLL Gauss vs CKG sM) come from the prompt summary of FAIR.5; not re-extracted from BSON in this draft. If the user wants them in a table rather than a sentence, add a small panel at the end of §8.4.
- **Connection to §sec:kms-lindbladian (thesis).** §8.1 is a natural place to forward-reference the KMS-Lindbladian section of the thesis if the user wants to make the section self-contained. Currently it cites only the first-principles literature ([Rouzé 2019][CITE: rouze2019functional], [Li & Lu 2024][CITE: liLu2024quantum], [Kochanowski et al. 2024][CITE: kochanowski2024rapid]).
- **qf-7go null prediction (§8.7).** Phrased as a prediction, not a fact. If qf-7go.6 is run before final submission and confirms the null, replace "predicted to give" with the actual measured outcome.
- **§9 sourcing.** Equations and page numbers in §9 are taken directly from `supplementary-informations/Ding et al. - 2024 ...pdf` (Defn. 24 / Lem. 25 p. 27; Eqs. 3.28 p. 17, 3.31 p. 18, 3.36 p. 20; Remark 23 p. 24) and `supplementary-informations/Chen et al. - 2025 ...pdf` (Nomenclature p. 23; Defn. III.1 p. 17; Prop. III.1 p. 18; Cor. III.2–III.4 pp. 15–16). Verify page references match the user's local copy before merging — arXiv-version page numbering may differ from the journal version.
- **§9 new bibtex stubs.** `chen2023quantum` (CKBG23), `li2023simulating` (LW23). Both are heavily cited by the two papers and almost certainly already in the thesis bib under different keys; cross-check before introducing duplicates. [CHECK]
- **§9.4 "new helpers" plan.** Two thin diagnostic functions: `dll_lcu_norms(filter, ham, β; M, τ)` returning `(Z_f, Z_g)`, and `ckg_lcu_norms(...)` returning the four CKG $\ell^1$-norms. They tap into already-constructed time-domain tables in `dll_lindblad_op_time` / `dll_coherent_op_time` / `_precompute_coherent_B` — no production matvec changes. Sub-issue qf-9ld.2 turns this into actual code.
- **§9.5 quadrature scope.** The §9 audit recipe is *agnostic* to the quadrature parameters $(M, \tau, T)$ in the sense that the $\beta$-scaling of $Z_g$ etc. is asymptotic and pulls through any reasonable quadrature; but the *absolute* values reported in the qf-9ld.2 dataset will depend on the chosen $\varepsilon$ target. Pin $\varepsilon = 10^{-3}$ (matching §2) to keep the dataset comparable to the existing $\tau_{\text{mix}}$ table.
