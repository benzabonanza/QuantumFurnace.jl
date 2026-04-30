---
paper: "Cuff, Ding, Louidor, Lubetzky, Peres, Sly (2012)"
title: "Glauber Dynamics for the Mean-Field Potts Model"
arxiv: "1204.4503"
year: 2012
venue: "J. Stat. Phys. 149:432–477"
pdf: "supplementary-informations/classical-review/cuff-et-al-2012-mean-field-potts-glauber.pdf"

temperature: [high-T, intermediate-T, low-T]
commutativity: [commuting]
locality: [mean-field]
particle-statistics: [spin]
hamiltonian-models: [mean-field-Potts, Curie-Weiss-Potts]

paradigm: [classical-MCMC]
quantum-or-classical: [classical]

result-type: [mixing-time-upper, mixing-time-lower, cutoff, critical-scaling]
key-scaling: "t_mix = alpha_1(beta,q) n log n with cutoff window O(n) for beta < beta_s(q); t_mix = Theta(n^{4/3}) at beta = beta_s(q) with window n^{-2/3}; t_mix >= C_1 exp(C_2 n) for beta > beta_s(q); essential mixing time stays n log n down to beta_c(q) after excluding an exp(-Omega(n))-mass set of starts"

related: ["bovier-den-hollander-2015-metastability", "woodard-schmidler-huber-2009-torpid-tempering", "madras-zheng-2003-swapping-algorithm", "lubetzky-sly-2012-critical-ising-polynomial", "lubetzky-sly-2013-cutoff-ising-lattice", "sly-2010-computational-transition-uniqueness"]
---

# Glauber Dynamics for the Mean-Field Potts Model — Cuff, Ding, Louidor, Lubetzky, Peres, Sly (2012)

**One-sentence takeaway**: Sharp three-regime analysis of single-site Glauber on the Curie–Weiss $q$-state Potts model ($q\ge 3$): cutoff at $\alpha_1(\beta,q)\,n\log n$ for $\beta<\beta_s(q)$, exact $\Theta(n^{4/3})$ critical scaling at the *spinodal* $\beta=\beta_s(q)$ (the dynamical, *not* equilibrium, transition), and $e^{\Theta(n)}$ torpid mixing for $\beta>\beta_s(q)$ — the canonical first-order-transition torpid-mixing benchmark and a strict separation $\beta_s(q)<\beta_c(q)$ that distinguishes Potts from Ising.

## Setting

- **Configuration space.** Complete graph on $V=\{1,\dots,n\}$, single-spin space $Q=\{1,\dots,q\}$, $q\ge 3$ integer; $\Sigma_n = Q^V$.
- **Target measure.** Curie–Weiss Potts distribution
  $$ \mu_n(\sigma) = Z_{\beta,n}^{-1}\,\exp\!\Big\{(\beta/n)\!\!\sum_{u,v\in V}\!\!\mathbf{1}_{\sigma(u)=\sigma(v)}\Big\}, $$
  parametrized by inverse temperature $\beta\ge 0$. Static phase transition at $\beta_c(q)$ (explicitly known); for $q\ge 3$ this is *first-order* — the disordered $1/q$-uniform phase and the $q$ ordered phases coexist.
- **Dynamics.** Discrete-time single-site Glauber: pick $u\in V$ uniformly, resample $\sigma(u)$ from the conditional distribution given the rest. Each step is $O(1)$ work; total-variation mixing time $t_{\mathrm{mix}(\varepsilon)}(n)$.
- **Two thresholds.** $\beta_c(q)$ — *static* / equilibrium critical point (where ordered and disordered phases coexist in $\mu_n$). $\beta_s(q)$ — *dynamical* spinodal, defined (eq. 1.1) as the largest $\beta$ for which the equation $(1+(q-1)e^{2\beta(1-qx)/(q-1)})^{-1} = x$ has only the trivial solution $x=1/q$. Always $\beta_s(q)<\beta_c(q)$ for $q\ge 3$ (in contrast to Ising where $\beta_s(2)=\beta_c(2)=1$).

## Main Results

Verbatim from the abstract and Theorems 1–4:

- **Theorem 1 (subcritical cutoff).** For $q\ge 3$ integer and $\beta<\beta_s(q)$, Glauber exhibits cutoff at
  $$ t_{\mathrm{mix}}(n) = \alpha_1(\beta,q)\,n\log n, \qquad \alpha_1(\beta,q) = [\,2(1-2\beta/q)\,]^{-1}, $$
  with cutoff window $w_\varepsilon(n) = O_\varepsilon(n)$.

- **Theorem 2 (scaling around the spinodal).** Write $\beta(n) = \beta_s(q) - \xi(n)$ with $\xi(n)\to 0$. (i) If $\lim_n n^{2/3}\xi(n) = \infty$ then cutoff with
  $$ t_{\mathrm{mix}}(n) = \alpha_1(\beta(n),q)\,n\log n + \alpha_2(q)\,n/\sqrt{\xi(n)}, \qquad w_\varepsilon(n) = O_\varepsilon\!\big(n + \sqrt{n/\xi(n)^{5/2}}\big). $$
  (ii) If $0\le \liminf n^{2/3}\xi(n) \le \limsup n^{2/3}\xi(n) < \infty$ then **no cutoff** and
  $$ t_{\mathrm{mix}(\varepsilon)}(n) = \Theta_\varepsilon(n^{4/3}). $$
  At criticality $\beta=\beta_s(q)$ exactly, $t_{\mathrm{mix}} = \Theta(n^{4/3})$ with scaling window $n^{-2/3}$.

- **Theorem 3 (supercritical torpidity).** For $q\ge 3$ integer and any fixed $\beta>\beta_s(q)$, every $0<\varepsilon<1$, there exist $C_1,C_2>0$ such that
  $$ t_{\mathrm{mix}(\varepsilon)}(n) \;\ge\; C_1\exp(C_2 n) \quad \text{for all } n. $$

- **Theorem 4 (essential mixing extends past $\beta_s$).** For $q\ge 3$ and any $\beta<\beta_c(q)$ there exist $C_1,C_2>0$ and subsets $\widetilde{\Sigma}_n\subseteq\Sigma_n$ with $\mu_n(\Sigma_n\setminus\widetilde{\Sigma}_n)\le C_1 e^{-C_2 n}$ such that the chain restricted to starts in $\widetilde{\Sigma}_n$ has cutoff at $\alpha_1(\beta,q)\,n\log n$ with window $O_\varepsilon(n)$. The torpidity above $\beta_s$ is therefore caused by an exponentially-rare set of bad initial configurations; from $\mu_n$-typical starts the chain still mixes in $n\log n$ throughout $[\beta_s,\beta_c)$.

These four results give the first complete mixing analysis around the *dynamical* critical temperature for any model with a first-order phase transition, including the critical $n^{4/3}$ power law (compare $n^{3/2}$ window $\sqrt n$ for critical 2D Ising — a different universality class).

## Method

Order-parameter reduction to the proportions vector $S_t \in \mathcal{S} = \{x\in\mathbb{R}_+^q : \|x\|_1=1\}$, on which Glauber acts as a near-deterministic gradient flow on the free-energy surface $F_\beta(s) = -\beta\sum_i s_i^2 - \tfrac{1}{n}\sum_i s_i\log s_i$ plus mean-zero $O(1/\sqrt n)$ noise (Fig. 1, 2 in the paper visualize this). Subcritical cutoff (Thm 1) is proved by an $L^2$-contraction / drift argument on $S_t$; the $\Theta(n^{4/3})$ scaling at $\beta_s$ comes from a saddle-node bifurcation where the leading drift toward $1/q$ vanishes cubically (Fig. 3 — the green middle curve $\beta=\beta_s$ has zero drift to second order), turning the radial coordinate into a critical Ornstein–Uhlenbeck-like process whose hitting time scales as the cube root of the noise variance over the cubic restoring coefficient. Supercritical torpidity (Thm 3) is a conductance / bottleneck lower bound across the metastability barrier separating the disordered minimum from the $q$ ordered minima of $F_\beta$ once they exist, and Thm 4 isolates that bottleneck to an exponentially-rare basin of attraction in $\mu_n$.

## Quantum vs Classical

- **Baseline.** This paper sets the *unconditional classical lower bound* $t_{\mathrm{mix}} \ge e^{\Omega(n)}$ for single-site Glauber on a model with a first-order phase transition, on its native Hamiltonian (mean-field Potts, $q\ge 3$, $\beta\in(\beta_s,\beta_c)$). Unlike Ising's continuous transition where Glauber retains polynomial mixing at criticality (Lubetzky–Sly $n^{3/2}$ in 2D), the first-order transition is a *qualitative* obstacle: an exponential bottleneck appears strictly before the equilibrium transition, so the slowdown is not a mere "critical-slowing-down" power law but a metastability barrier between disordered and ordered phases that exists throughout the entire interval $\beta\in(\beta_s,\beta_c)$.

- **Source of the difference.** Geometric energy barrier in the order-parameter $S_t$. The free-energy surface $F_\beta$ has a single minimum at $1/q$ for $\beta<\beta_s$, develops $q$ new local minima at $\beta=\beta_s$ via a saddle-node bifurcation, and crosses to global minima at $\beta=\beta_c$. The Glauber dynamics is essentially an SDE on $\mathcal{S}$ with $O(1/\sqrt n)$ noise, so crossing a free-energy barrier of height $\Delta F = \Theta(1)$ in $S$ takes $e^{\Theta(n)}$ steps by Eyring–Kramers (cf. `bovier-den-hollander-2015-metastability`).

- **Quantum side: open candidate cell.** The paper makes no quantum claim. The natural question for the thesis: **can a quantum Gibbs sampler (e.g. KMS detailed-balance Lindbladian) avoid the metastability bottleneck on the same Hamiltonian?** Three structural reasons this *might* admit a quantum speedup:
  1. Mean-field Potts has commuting (classical) Hamiltonian, so this is *not* a sign-problem cell — but the hardness here is dynamical (metastability), not sign-problem-related, so it lives in a different axis from the Troyer–Wiese cell.
  2. The bottleneck is a finite-dimensional ($q$-component order-parameter) saddle, in principle accessible to Grover-style $\sqrt{e^{\beta\Delta F}} = e^{\beta\Delta F/2}$ quantum walks if the saddle can be coherently traversed.
  3. Coherent quantum mixing on classical Potts has not been ruled out: no quantum lower bound matching $e^{\Omega(n)}$ is known. [CHECK] whether any quantum Gibbs-sampling lower bound on commuting first-order-transition models exists in the corpus.

- **Caveat — known no-go for adjacent classical accelerations.** Bhatnagar–Randall (2004) and `woodard-schmidler-huber-2009-torpid-tempering` show that *parallel and simulated tempering* are also torpid on this exact model: any classical method that operates by replica-exchange across a temperature ladder fails for the same metastability reason. Galanis–Štefankovič–Vigoda (2015) gives matching $\#\mathrm{BIS}$-hardness for antiferromagnetic Potts in the tree non-uniqueness region, suggesting the difficulty is not just dynamical but partly computational. So a quantum sampler that beats this baseline must do something genuinely non-classical (coherent tunneling), not merely a smarter classical move set.

- **Caveat — essential mixing.** Theorem 4 weakens the worst-case lower bound: from $\mu_n$-typical starts, Glauber *does* mix in $n\log n$ throughout $(\beta_s,\beta_c)$. So the practical question of "sampling from the equilibrium distribution given a typical sample" has classical complexity $\widetilde O(n)$ even in the supposedly-torpid regime. Quantum advantage would have to attack the worst-case cold-start regime (e.g. starting from the all-equal state below $\beta_s$ and warming above $\beta_s$).

## Implications for Quantum Advantage

- **Regime cell (thesis 3×3 grid).** **Classical Ising/Potts × intermediate-T** (the first-order-coexistence interval $(\beta_s,\beta_c)$) and **× low-T** (within the ordered phase, $\beta>\beta_c$, where the $e^{\Theta(n)}$ bound from Thm 3 still applies). Mean-field locality (no geometry); $q$-state spins.
- **What this changes.** Provides a *concrete classical lower bound to beat*: $e^{\Theta(n)}$ for worst-case-start sampling on mean-field Potts at $\beta>\beta_s(q)$ for any $q\ge 3$. Together with Bhatnagar–Randall / Woodard–Schmidler–Huber (tempering torpid) and Galanis–Štefankovič–Vigoda (antiferromagnetic Potts hardness), this triple cements first-order Potts as a *robustly hard* classical sampling cell that no known classical method handles.
- **Promising or not (for quantum advantage).** **Promising — open candidate cell.** First-order-transition metastability is precisely the kind of obstacle a quantum walk could in principle bypass via coherent tunneling, and the underlying Hamiltonian is sign-problem-free, so the difficulty is purely dynamical (no fermionic / non-stoquastic complications muddying the comparison). The cell is small (mean-field Potts is one model) but cleanly delineates the *type* of bottleneck where quantum speedup is conceivable. Conversely, if KMS Lindbladians inherit the same $e^{\Omega(n)}$ obstruction here, that would be strong evidence that the noncommutative log-Sobolev framework cannot resolve first-order transitions either.

## Open Questions / Limitations

- The analysis is *exclusively* mean-field (complete graph). On lattices $\mathbb{Z}^d$ ($d\ge 2$, $q$ large enough) Potts also has a first-order transition but the dynamical picture is expected to be qualitatively different: nucleation has $O(1)$ lifetime so the transition to equilibrium proceeds via local nucleation + droplet growth, and the mixing time is expected to be governed by surface tension $\tau_\beta$ as $e^{(\tau_\beta+o(1))n}$ for low $T$ in a box (related work, p. 5). The mean-field $e^{\Theta(n)}$ bound here is therefore not directly a lattice statement.
- Cluster algorithms (Swendsen–Wang) on the same model: Gore–Jerrum (1999) show SW is *also* torpid on mean-field Potts at the first-order transition, despite its rapid-mixing reputation on ferromagnetic Ising. Closing the algorithmic gap on first-order Potts on *any* graph class is open.
- $\Theta(n^{4/3})$ at $\beta_s$ is sharp in a finite window of size $n^{-2/3}$ around $\beta_s$; outside that window one of the asymptotic regimes of Thm 1 / Thm 3 applies. Universality of the $4/3$ exponent for other first-order-transition models is conjectured but not proved here.
- No quantum analysis. The KMS-detailed-balance Lindbladian on the same Hamiltonian — what does its mixing time look like at the spinodal? Would a coherent quantum walk on the order-parameter simplex realize $\sqrt{e^{\beta\Delta F}}$ scaling, or does the Lindbladian gap still scale as $e^{-\Theta(n)}$? This is exactly the cross-cell question the thesis Review chapter foregrounds.

## Connections

- **Same model, alternative classical samplers (all torpid).** Bhatnagar–Randall 2004 (simulated tempering torpid on mean-field Potts $q\ge 3$); Gore–Jerrum 1999 (Swendsen–Wang torpid on the same model — cluster moves do *not* help here); `woodard-schmidler-huber-2009-torpid-tempering` (general no-go for parallel and simulated tempering, with mean-field Potts as a flagship example) — together these papers form the *classical no-go frontier* for first-order Potts. **No existing review for these in the corpus yet** ([CHECK] and add when reviewed).
- **Same model, positive classical results.** `madras-zheng-2003-swapping-algorithm` proves polynomial PT mixing on mean-field *Ising* (which has no first-order transition); the contrast with mean-field Potts is exactly the structural point — PT works on continuous transitions, fails on first-order ones.
- **Hardness companion.** Galanis–Štefankovič–Vigoda 2015 establishes $\#\mathrm{BIS}$-hardness for *antiferromagnetic* Potts in the tree non-uniqueness region. Together with Cuff et al. this gives a "hard from both sides" picture for $q\ge 3$ Potts: ferromagnetic mean-field is dynamically hard, antiferromagnetic bounded-degree is computationally hard. **No existing review.**
- **Same author school, contrasting model.** `lubetzky-sly-2013-cutoff-ising-lattice` (Glauber-Ising cutoff at $(d/2\lambda_\infty)\log n$ throughout high-T) and `lubetzky-sly-2012-critical-ising-polynomial` (polynomial $n^{3/2}$ mixing of critical 2D Ising) — the contrast crystallizes "Ising = continuous = polynomial at criticality, Potts = first-order = exponential at coexistence".
- **Low-T metastability framework.** `bovier-den-hollander-2015-metastability` provides the Eyring–Kramers $K e^{\beta\Gamma^*}$ machinery that explains the $e^{\Theta(n)}$ in Thm 3 quantitatively. The Bovier–den Hollander book covers Curie–Weiss Ising in Ch. 13–15; Cuff et al. is the natural Potts companion for $q\ge 3$.
- **Hardness anchor on a different axis.** `sly-2010-computational-transition-uniqueness` connects equilibrium phase transitions to *computational* hardness in the antiferromagnetic two-spin setting; Cuff et al. establishes the analogous *dynamical* phase transition in the ferromagnetic mean-field Potts setting — different axis (dynamics vs counting), same physics-meets-complexity flavour.
- **Within the thesis grid.** Cited in `classical-gibbs-sampling-corpus.md` §5 (hardness / no-go) and indirectly in §4 (low-T metastability) and the intermediate-T column of the Tier-1 grid. Quantum-side counterpart to define the comparison: Chen et al. 2023/2025 KMS Lindbladians on mean-field Potts (no published numerics specific to first-order transitions yet — [CHECK] thesis numerics chapter targets this gap).
