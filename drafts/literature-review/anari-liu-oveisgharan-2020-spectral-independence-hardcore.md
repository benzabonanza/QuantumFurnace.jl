---
paper: "Anari, Liu, Oveis Gharan (2020)"
title: "Spectral Independence in High-Dimensional Expanders and Applications to the Hardcore Model"
arxiv: "2001.00303"
year: 2020
venue: "FOCS 2020"
pdf: "supplementary-informations/classical-review/anari-liu-oveisgharan-2020-spectral-independence-hardcore.pdf"

temperature: [high-T, intermediate-T]
commutativity: [commuting]
locality: [sparse, geometric-local]
particle-statistics: [n/a]
hamiltonian-models: [hardcore, antiferromagnetic-2-spin]

paradigm: [classical-MCMC]
quantum-or-classical: [classical]

result-type: [mixing-time-upper, gap-bound, structural]
key-scaling: "t_mix(eps) = O((1+lambda) n)^{1+C(delta)} log(1/(eps mu(tau))), with C(delta) <= exp(O(1/delta)) and lambda < (1-delta) lambda_c(Delta); independent of max-degree Delta"

related: ["chen-liu-vigoda-2021-optimal-mixing-glauber-hde", "weitz-2006-counting-independent-sets-tree-threshold", "sly-2010-computational-transition-uniqueness", "chen-eldan-2022-localization-schemes"]
---

# Spectral Independence in High-Dimensional Expanders and Applications to the Hardcore Model — Anari, Liu, Oveis Gharan 2020

**One-sentence takeaway**: Defines *spectral independence* — a uniform spectral bound on the pairwise-influence matrix $\Psi_\mu$ — and uses high-dimensional-expander local-to-global to prove polynomial-time mixing of single-site Glauber for the hardcore model up to the tree-uniqueness threshold $\lambda_c(\Delta)$ on any bounded-degree graph, matching Weitz/Sly and replacing Weitz's quasi-polynomial FPTAS with an FPRAS independent of $\Delta$.

## Setting

- **Target**: hardcore (independent-set) distribution on a graph $G=(V,E)$ with $|V|=n$, fugacity $\lambda>0$, $\mu(I) \propto \lambda^{|I|}$ for independent sets $I$. Partition function $Z_G(\lambda) = \sum_I \lambda^{|I|}$.
- **Hamiltonian class**: 2-spin antiferromagnetic, *classical* (commuting). The hardcore model is the canonical antiferromagnetic two-state spin system.
- **Temperature regime**: tree-uniqueness regime $\lambda < \lambda_c(\Delta) \;{\stackrel{\rm def}{=}}\; (\Delta-1)^{\Delta-1}/(\Delta-2)^\Delta \approx e/(\Delta-2)$, i.e. above (in $T$) any phase-transition activity.
- **Access**: explicit graph $G$ with maximum degree bounded by $\Delta$; classical Glauber update (single-site heat-bath on the in/out indicator).

## Main Results

- **Theorem 1.3 (spectral independence $\Rightarrow$ spectral gap)**: for a $(\eta_0,\dots,\eta_{n-2})$-spectrally independent $d$-homogeneous distribution $\mu$ on $2^{[n]}$, the natural Glauber chain has spectral gap
$$\mathrm{gap}(P) \ge \frac{1}{n} \prod_{i=0}^{n-2}\!\Big(1-\frac{\eta_i}{n-i-1}\Big).$$
Bounded $\eta_i = O(1)$ gives $\mathrm{gap} \gtrsim 1/n^{1+O(1)}$, hence near-linear mixing.
- **Theorem 1.5 (spectral independence $\Rightarrow$ local spectral expansion)**: $\mu$ being $(\eta_0,\dots,\eta_{n-2})$-spectrally independent implies the auxiliary $n$-partite simplicial complex $X^\mu$ is a $(\eta_0/(n-1),\dots,\eta_{n-2}/1)$-local spectral expander.
- **Theorem 1.8 (hardcore is spectrally independent below uniqueness)**: there is a function $C : [0,1] \to \mathbb{R}_{>0}$ with $C(\delta) \le \exp(O(1/\delta))$ such that for any graph of max degree $\le \Delta$, any $0<\delta<1$, and $\lambda = (1-\delta)\lambda_c(\Delta)$, the hardcore distribution is $(\eta_i)$-spectrally independent with $\eta_i \le \min\{C(\delta),\, \lambda(n-i-1)/(1+\lambda)\}$.
- **Corollary 1.9 (FPRAS)**: for every $\delta>0$ and every graph $G$ of max degree $\le \Delta$, there is an FPRAS for $Z_G(\lambda)$ at $\lambda = (1-\delta)\lambda_c(\Delta)$.
- **Remark 1.10 (mixing time, explicit)**: starting from $\tau$,
$$t_\tau(\varepsilon) \le O\!\Big(\big((1+\lambda) n\big)^{1+C(\delta)} \log\!\big(\tfrac{1}{\varepsilon \mu(\tau)}\big)\Big), \quad C(\delta) \le \exp(O(1/\delta)),$$
giving $t_\emptyset(\varepsilon) \le O\big(n^{2+C(\delta)} \log(1/\varepsilon)\big)$. Crucially **no $\Delta$ dependence**.

## Method

Three-step reduction: (1) introduce the *signed pairwise influence matrix* $\Psi_\mu(i,j) = \Pr[j\mid i] - \Pr[j\mid \bar i]$, define spectral independence as $\lambda_{\max}(\Psi_\mu) \le \eta$; (2) show via the auxiliary $n$-partite simplicial complex $X^\mu$ that spectral independence is equivalent to local spectral expansion of $X^\mu$, then invoke the Alev–Lau local-to-global theorem to bound the gap of the top-level random walk (which is exactly Glauber); (3) bound the row-sums of $\Psi_\mu$ for hardcore via Weitz's self-avoiding-walk tree plus a potential-function correlation-decay estimate, giving $\sum_{u\ne v}|\Psi_\mu(u,v)| \le C(\delta)$.

## Quantum vs Classical

- **Baseline**: Weitz's correlation-decay FPTAS [Wei06] on graphs of max degree $\le\Delta$, which runs in time $(n/\varepsilon)^{C(\delta)\log\Delta}$ — *quasi-polynomial* in $n$ once $\Delta$ is unbounded. Prior best FPRAS (Markov-chain) only reached $\lambda < 2/(\Delta-2)$ [Vig01], i.e. far from the threshold.
- **Gap**: polynomial $\to$ polynomial improvement, but in *unbounded-$\Delta$* graphs the improvement is super-polynomial: $n^{C(\delta)\log\Delta}$ to $n^{O(1)+C(\delta)}$. The only $\Delta$-dependence remaining is through $\lambda_c(\Delta)\le 4$.
- **Source of the difference**: high-dimensional expansion replaces tree-recursion truncation by a global spectral argument; influence sums (not tree recursions) are the only $\Delta$-bounded quantity needed.
- **Quantum side**: paper makes no quantum claim. The relevance for the thesis is that hardcore is a classical commuting model — the "everything commutes, classical sign-free" baseline cell. Quantum Gibbs samplers (Chen et al. KMS, Ding et al.) reduce to a near-optimal classical chain in this regime modulo block-encoding overhead $O(\log(1/\varepsilon)$ in coherent block-encoding cost). No proven quantum speedup is known on hardcore at $\lambda<\lambda_c$, and any quantum advantage would have to beat $\widetilde{O}(n^2)$ classical mixing — the bar set here. [CHECK] whether any quantum lower bound for sampling commuting Gibbs measures appears in the corpus.
- **Caveat**: the constant $C(\delta) = \exp(O(1/\delta))$ blows up as $\lambda \uparrow \lambda_c$; the FPRAS is *not* a uniform $O(n\log n)$ chain (that sharpening is Chen–Liu–Vigoda 2021).

## Implications for Quantum Advantage

- **Regime cell**: high-T (tree-uniqueness regime) $\times$ commuting (classical Hamiltonian) $\times$ sparse / bounded-degree $\times$ classical spin (n/a particle stats).
- **What this changes**: *closes the upper-bound side* of the classical / commuting / high-T cell — together with Sly hardness above $\lambda_c$, this is the canonical "computational threshold = phase-transition threshold" benchmark in the corpus.
- **Promising or not (for quantum advantage)**: explicitly *not promising*. This cell is a high-confidence "classical wins or ties" cell. Any quantum advantage in the thesis chapter must come from non-stoquastic / low-T cells where this entire framework provides no analog.

## Open Questions / Limitations

- $C(\delta) \le \exp(O(1/\delta))$ in influence sum; the paper conjectures (Remark 1.14) that $C(\delta) = O(1/\delta)$ should be tight; tightness on the $\Delta$-regular tree is shown in Appendix B.
- Mixing time exponent $1+C(\delta)$ in $n$ is far from optimal $O(n\log n)$; sharpened to $O(n\log n)$ by Chen–Liu–Vigoda (2021) via modified-log-Sobolev / entropic independence.
- No extension to multi-spin or continuous-spin systems in this paper (later: Chen et al. 2020, Feng et al. 2020 for $q$-colourings; Anari–Jain et al. for entropic independence).
- Pure classical: no statement about quantum spin systems; the spectral-independence framework has no immediate quantum analogue when the Hamiltonian is non-commuting.

## Connections

- **Sharpens** Weitz (2006): replaces a quasi-polynomial deterministic FPTAS by a polynomial-time *randomized* algorithm (with no $\Delta$ dependence in the runtime exponent), using the same self-avoiding-walk tree internally.
- **Matches** Sly (2010) hardness: rapid mixing exactly up to $\lambda_c(\Delta)$, where Sly proves NP-hardness of approximation past $\lambda_c$. Together the two papers cement "computational threshold = uniqueness threshold" for hardcore.
- **Foundation for** Chen–Liu–Vigoda (2021), which promotes spectral independence to a modified log-Sobolev inequality and recovers optimal $O(n\log n)$ mixing.
- **Foundation for** Chen–Eldan (2022) localization-schemes framework: spectral independence is the discrete instance of the same martingale-decomposition idea that drives stochastic localization.
- Within the thesis grid: this is the *classical Ising/Potts $\times$ high-T* anchor cited in `classical-gibbs-sampling-corpus.md` §6 alongside Eldan–Koehler–Zeitouni and Chen–Liu–Vigoda. No quantum-side paper in the corpus directly competes with it (because the model is classical commuting); the relevant quantum-vs-classical question reduces to: can quantum Gibbs samplers ever beat $\widetilde{O}(n^2)$ Glauber on a sparse classical Hamiltonian below uniqueness? Current evidence says no.
