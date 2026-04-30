---
paper: "Lubetzky, Sly (2013)"
title: "Cutoff for the Ising model on the lattice"
arxiv: "0909.4320"
year: 2013
venue: "Invent. Math. 191:719-755"
pdf: "supplementary-informations/classical-review/lubetzky-sly-2013-cutoff-ising-lattice.pdf"

temperature: [high-T]
commutativity: [commuting]
locality: [geometric-local]
particle-statistics: [spin]
hamiltonian-models: [Ising]

paradigm: [classical-MCMC]
quantum-or-classical: [classical]

result-type: [mixing-time-upper, mixing-time-lower, structural]
key-scaling: "TV cutoff at t_mix = (d/(2 lambda_infinity)) log n + O(log log n) on the d-dim torus (Z/n Z)^d at any temperature with strong spatial mixing; cutoff window O(log log n)"

related: ["lubetzky-sly-2012-critical-ising-polynomial", "eldan-koehler-zeitouni-2022-spectral-condition-ising", "bauerschmidt-dagallier-2024-near-critical-ising-lsi", "chen-eldan-2022-localization-schemes"]
---

# Cutoff for the Ising model on the lattice — Lubetzky & Sly 2013

**One-sentence takeaway**: For Glauber dynamics on the $d$-dimensional Ising model with periodic boundary conditions at any temperature in the strong-spatial-mixing regime, total-variation mixing exhibits *cutoff* at $t_n = \frac{d}{2\lambda_\infty}\log n$ with window $O(\log\log n)$, where $\lambda_\infty$ is the spectral gap of the dynamics on the infinite-volume lattice — proved via a new $L^1$-to-$L^2$ reduction that lets log-Sobolev / spectral-gap machinery control TV distance and not just $\chi^2$ distance.

## Setting

- **Target**: Ising Gibbs measure $\pi(\sigma) \propto \exp(\beta \sum_{x\sim y}\sigma_x \sigma_y + h\sum_x \sigma_x)$ on the $d$-dimensional discrete torus $(\mathbb{Z}/n\mathbb{Z})^d$ with periodic boundary conditions; $|\Lambda|=n^d$ spins. The natural $n$ here is the *side length* of the torus, not the number of spins.
- **Hamiltonian class**: classical, commuting, ferromagnetic nearest-neighbour Ising on $\mathbb{Z}^d$. $d \ge 1$ arbitrary.
- **Temperature regime**: any $\beta$ at which **strong spatial mixing** (SSM) holds for the infinite-volume measure on $\mathbb{Z}^d$. For $d=2$ this is the entire subcritical regime $\beta < \beta_c$ (by Martinelli–Olivieri–Schonmann); for $d \ge 3$ it is conjecturally subcritical but proved only at high enough temperature.
- **Access**: continuous-time single-site Glauber / heat-bath dynamics on the spin configuration; classical MCMC.
- **Quantity controlled**: total-variation mixing time $t_{\mathrm{mix}}^{(n)}(\varepsilon) = \min\{t : \max_{\sigma_0}\|P_t(\sigma_0,\cdot) - \pi\|_{\mathrm{TV}} \le \varepsilon\}$, and the cutoff phenomenon $t_{\mathrm{mix}}^{(n)}(\varepsilon)/t_{\mathrm{mix}}^{(n)}(1-\varepsilon) \to 1$.

## Main Results

- **Theorem 1 (cutoff at the infinite-volume rate, verbatim spirit)**: continuous-time Glauber dynamics on $(\mathbb{Z}/n\mathbb{Z})^d$ at any $(\beta, h)$ in the SSM regime exhibits cutoff at
  $$ t_n = \frac{d}{2\lambda_\infty}\log n, \qquad \text{window } O(\log\log n), $$
  where $\lambda_\infty := \lim_{L\to\infty}\lambda(\Lambda_L)$ is the spectral gap of Glauber on $\mathbb{Z}^d$ (proved finite and positive in the SSM regime). Concretely $t_{\mathrm{mix}}^{(n)}(\varepsilon) = t_n + O_\varepsilon(\log\log n)$ for every fixed $\varepsilon \in (0,1)$.
- **Corollary (2D up to $\beta_c$)**: on $\mathbb{Z}^2$, SSM holds for all $\beta < \beta_c$ (MOS '94), so cutoff at $\frac{1}{\lambda_\infty}\log n$ with window $O(\log\log n)$ persists *all the way up to* (but not at) the critical temperature.
- **Lower bound (Section 4)**: $t_{\mathrm{mix}}^{(n)}(1-\varepsilon) \ge t_n - C_\varepsilon \log\log n$ via a distinguishing-statistic / second-moment argument on a localized magnetization, using $\lambda_\infty$ as the asymptotic decay rate of correlations.
- **Upper bound (Section 5)**: $t_{\mathrm{mix}}^{(n)}(\varepsilon) \le t_n + C_\varepsilon \log\log n$ via the new $L^1 \to L^2$ reduction (see Method).
- **Structural result ($L^1$-to-$L^2$ reduction, Lemma 3.3)**: for any reversible Markov chain whose stationary measure satisfies a finite-range product structure with SSM, TV distance at time $t = t_0 + s$ is controlled by the $L^2$ distance at time $s$ on a "sparsified" subsystem, where $t_0 = O(\log\log n)$ is a *burn-in* during which most spins are coupled to their stationary values. This is the technical engine.

## Method

Three-step argument tied together by a new $L^1\!\to\!L^2$ reduction. First (burn-in, $\sim \log\log n$ steps), use a censoring / coupling argument à la Peres–Winkler to show that after a short prefix all but an $\varepsilon$-fraction of spins have decoupled from the initial condition and follow the stationary measure conditionally on a sparse "still-coupled" set $S$. Second, on $S$ the conditional measure is approximately a product of small clusters thanks to SSM, so the residual TV distance equals the TV distance of a *low-dimensional* subsystem and can be majorised by its $L^2 = \chi^2$ distance via Cauchy–Schwarz with a tractable prefactor (the cardinality of $S$). Third, $\chi^2$ on the subsystem decays at rate $2\lambda_\infty$ by the standard spectral-gap / log-Sobolev decomposition for Glauber, so the total bound becomes $t_n + O(\log\log n)$. The key novelty is that the reduction *removes* the usual $\sqrt{|\Omega|} = \exp(\Theta(n^d))$ Cauchy–Schwarz overhead between $L^2$ and TV that had previously prevented spectral methods from giving sharp TV cutoff times.

## Quantum vs Classical

- **Baseline**: at the time of writing, the best classical TV mixing bounds in the SSM regime were $O(n^d \log n)$ (Stroock–Zegarlinski / Martinelli–Olivieri via uniform LSI), with the right *order* but a constant factor off and with no cutoff result. Lubetzky–Sly identifies the sharp leading constant $d/(2\lambda_\infty)$ and proves cutoff with an explicit $O(\log\log n)$ window.
- **Gap (vs quantum)**: there is no quantum candidate that beats this on classical commuting Ising in the SSM regime. Quantum Gibbs samplers (Chen et al. KMS DB Lindbladian, Ding et al.) reduce to a near-classical chain on a commuting target and pick up at best polylog block-encoding overhead per step, with no known way to lower the leading constant from $d/(2\lambda_\infty)$. The bar set here for quantum is therefore *constant-tight*: $\frac{d}{2\lambda_\infty}\log n$ Glauber sweeps, equivalently $n^d \cdot \frac{d}{2\lambda_\infty}\log n$ single-site updates in discrete time.
- **Source of the (non-)difference**: SSM is an intrinsically classical correlation-decay condition that has no analogue for non-commuting $H$. The $L^1$-to-$L^2$ trick uses the product structure of the classical configuration space and the *real* spectral decomposition of the reversible Glauber semigroup; the natural quantum lift (KMS-DB Lindbladian) has a non-self-adjoint generator whose $L^1$-to-$L^2$ control is significantly more delicate (Capel et al. 2025 give the noncommutative-OT analogue at the level of spectral-gap-implies-mixing, but no cutoff).
- **Caveat**: the result is conditional on SSM, which is open in $d \ge 3$ in part of the subcritical regime; the constant $\lambda_\infty$ is only computed numerically and depends on $\beta$ in an unknown way (no closed form). Cutoff is *not* claimed at criticality — at $\beta = \beta_c$ the relaxation is polynomial (Lubetzky–Sly 2012) and the leading constant changes.

## Implications for Quantum Advantage

- **Regime cell**: classical Ising × high-T (subcritical) × geometric-local × spin. This is the *canonical anchor* for the corpus's classical-Ising / high-T cell from the lattice / geometric-local angle; complements Eldan–Koehler–Zeitouni (dense / mean-field $J$) and Anari–Liu–Oveis Gharan / Chen–Liu–Vigoda (sparse / bounded-degree antiferromagnetic).
- **What this changes**: tightens the leading constant. Previous bounds gave $O(n^d \log n)$; Lubetzky–Sly gives a *tight* $\frac{d}{2\lambda_\infty} n^d \log n$ (single-site discrete steps) with cutoff. Any quantum-vs-classical comparison in this cell must therefore beat both the leading constant *and* the $\log n$ factor; merely matching the order is no longer competitive evidence for advantage.
- **Promising or not** (for quantum advantage): explicitly *not* promising. The strong-spatial-mixing regime is the sharpest possible "classical wins" cell — cutoff means the classical chain mixes essentially as fast as any reversible chain on this state space ever could, up to lower-order terms. The thesis's punchline cell (non-stoquastic × low-T) sits structurally outside SSM (correlations do not decay at low-$T$ across phase coexistence), so this paper does not constrain the quantum-advantage candidate cell — but it does close the door on advantage in the *high-T classical* anchor.
- **What it does *not* cover**: criticality (handled separately by Lubetzky–Sly 2012 with a polynomial, non-cutoff, non-log-Sobolev mixing bound on 2D); low-$T$ (where Bovier–den Hollander metastability $\sim e^{\beta\Delta}$ takes over); spin glasses and random / antiferromagnetic Ising (no SSM); free / non-periodic boundary conditions (where boundary fluctuations change the leading constant).

## Open Questions / Limitations

- **SSM in $d \ge 3$ all the way to $\beta_c$**: open. Cutoff at the asserted constant is conjectural in the gap between proven SSM and $\beta_c$.
- **Cutoff at criticality**: not addressed. Critical 2D Ising mixes polynomially (Lubetzky–Sly 2012) but cutoff there is open.
- **Low temperature**: the techniques fail at $\beta > \beta_c$ — both because $\lambda_\infty \to 0$ (degenerate) and because SSM fails. The low-$T$ cell is governed by metastability (Bovier–den Hollander), not by cutoff theory.
- **Boundary conditions**: cutoff constant depends on boundary conditions; only periodic torus is treated. Free or plus boundary conditions modify the leading constant via boundary effects.
- **Quantum analogue**: no quantum cutoff theorem is known for KMS-DB Lindbladians in any nontrivial regime. The natural conjecture — cutoff at $\frac{d}{2\lambda_\infty^{\mathrm{KMS}}}\log n$ for high-T quantum Ising — is open even at the level of the right leading constant.

## Connections

- **Lubetzky–Sly 2012** [`lubetzky-sly-2012-critical-ising-polynomial`]: the *critical* sister paper. At $\beta = \beta_c$ on $\mathbb{Z}^2$, Glauber mixes in polynomial (rather than logarithmic) time; cutoff is not established. Together the two papers cover subcritical (cutoff) and critical (polynomial, no cutoff) 2D Ising.
- **Lubetzky–Sly 2016 (information percolation)**: extends cutoff to *all* subcritical $\beta$ in any dimension via a different (information-percolation) cluster expansion, removing the SSM assumption in $d \ge 3$. Where the 2013 paper relies on functional inequalities + $L^1\!\to\!L^2$, the 2016 paper rebuilds the proof on direct space-time coupling.
- **Eldan–Koehler–Zeitouni 2022** [`eldan-koehler-zeitouni-2022-spectral-condition-ising`]: covers the *dense / mean-field* slice of high-T classical Ising via stochastic localization; complementary anchor in the same corpus cell. Both papers ultimately rely on Glauber spectral gaps but apply to disjoint geometries (lattice vs spectral-norm-bounded $J$).
- **Bauerschmidt–Dagallier 2024** [`bauerschmidt-dagallier-2024-near-critical-ising-lsi`]: gives a *uniform-in-volume LSI* up to $\beta_c$ for ferromagnetic Ising on $\mathbb{Z}^d$ ($d \ge 5$). Provides the LSI input needed to push Lubetzky–Sly-style $L^1\!\to\!L^2$ arguments closer to $\beta_c$, but does not by itself give cutoff.
- **Chen–Eldan 2022** [`chen-eldan-2022-localization-schemes`]: localization-schemes framework subsumes much of the spectral-gap-implies-mixing technology; recasts $L^1\!\to\!L^2$ reductions inside a martingale decomposition. Does not currently give a cutoff theorem on its own.
- **Stroock–Zegarlinski 1992 / Martinelli–Olivieri 1994**: the predecessors. Established the equivalence SSM $\Leftrightarrow$ uniform LSI and the $O(\log n)$ relaxation timescale; Lubetzky–Sly upgrades the constant from "some finite $C$" to the sharp $d/(2\lambda_\infty)$ and adds cutoff.
- **Levin–Peres–Wilmer (textbook)**: standard reference for the cutoff phenomenon and $L^2 = \chi^2$ vs TV machinery. Lubetzky–Sly gives the first cutoff for a chain whose stationary measure is itself nontrivial (no exact diagonalisation) — a genuine breakthrough relative to the random-walk / shuffling cases in the textbook.
- **Within the thesis grid**: this is the *Tier-1 anchor for classical Ising × high-T (lattice / geometric-local slice)* in `classical-gibbs-sampling-corpus.md` §2. Sets the constant-tight classical baseline that any quantum Gibbs sampler in this cell must beat — and signals that any honest quantum-advantage claim cannot live in the SSM regime, where classical cutoff makes the comparison unforgiving.
