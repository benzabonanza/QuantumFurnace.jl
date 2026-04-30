---
paper: "Chen and Eldan (2022)"
title: "Localization Schemes: A Framework for Proving Mixing Bounds for Markov Chains"
arxiv: "2203.04163"
year: 2022
venue: "FOCS 2022"
pdf: "supplementary-informations/classical-review/chen-eldan-2022-localization-schemes.pdf"

temperature: [high-T, intermediate-T]
commutativity: [commuting]
locality: [n/a]
particle-statistics: [spin]
hamiltonian-models: [Ising, hardcore, log-concave, Sherrington-Kirkpatrick, ferromagnetic-Ising]

paradigm: [classical-MCMC]
quantum-or-classical: [classical]

result-type: [mixing-time-upper, gap-bound, mlsi, structural]
key-scaling: "Hardcore in delta-uniqueness: t_mix <= exp(c/delta)(n log n + 3n log(1/eps)). Ising at ||J||_op < 1/2: t_mix = O(n log n). Graphical Ising in tree-uniqueness: rho_LS >= e^{-8/delta}/n, so t_mix = O(n log n). Strongly log-concave: rho_LS(P^{RGD_eta}) >= mu/(mu+1/eta). Near-critical FM Ising: rho_LS >= n^{-1} exp(-2 ||J||_op int_0^1 chi(lambda) dlambda)"

related: ["anari-liu-oveisgharan-2020-spectral-independence-hardcore", "chen-liu-vigoda-2021-optimal-mixing-glauber-hde", "eldan-koehler-zeitouni-2022-spectral-condition-ising", "bauerschmidt-dagallier-2024-near-critical-ising-lsi", "weitz-2006-counting-independent-sets-tree-threshold", "sly-2010-computational-transition-uniqueness"]
---

# Localization Schemes: A Framework for Proving Mixing Bounds for Markov Chains — Chen and Eldan 2022

**One-sentence takeaway**: Introduces *localization schemes* — a martingale of probability measures $(\nu_t)_t$ with $\nu_0=\nu$ that "zooms in" toward Dirac masses — as a single framework that simultaneously generalizes spectral / entropic independence (discrete) and stochastic localization (continuous), and uses it to prove the first $O(n\log n)$ Glauber mixing for hardcore in tree-uniqueness, optimal $O(n\log n)$ mixing for graphical Ising in uniqueness with no external-field dependence, and KL-decay for log-concave sampling via the Restricted Gaussian Oracle.

## Setting

- **Targets**: probability measures on $\Omega \in \{\{-1,1\}^n, \mathbb{R}^n\}$ — discrete spin systems (hardcore, Ising, ferromagnetic Ising near criticality) and continuous strongly log-concave measures.
- **Hamiltonian class**: classical (commuting) — antiferromagnetic 2-spin (hardcore), graphical Ising $\nu_{G,\beta}(x) \propto \exp(\langle x,v\rangle + \beta \sum_{(i,j)\in E}\mathbf{1}_{x_i=x_j})$, Ising with operator-norm-bounded interaction $\|J\|_{\mathrm{op}}<1/2$ (Sherrington–Kirkpatrick at $\beta\le 1/4$), strongly log-concave $\nu \propto e^{-V}$ with $\nabla^2 V \succeq \mu I_n$.
- **Temperature regime**: tree-uniqueness for hardcore / graphical Ising; spectral-norm regime $\|J\|_{\mathrm{op}}<1/2$ for general Ising; high-T for SK; entire ferromagnetic regime up to criticality (susceptibility $\chi(\lambda)$ controls the bound).
- **Access**: classical Glauber dynamics $P^{\mathrm{GD}}$, or Restricted Gaussian Oracle / restricted Gaussian dynamics for the continuous case.

## Main Results

- **Theorem 49 (general Ising MLSI)**: for $\nu$ on $\{-1,1\}^n$ and PD interaction $J$, if $\|\mathrm{Cov}(\mu_{\lambda,v})\|_{\mathrm{op}} \le \alpha(\lambda)$ for all $\lambda\in[0,1], v\in\mathbb{R}^n$ along the stochastic-localization tilt path, and the product measure at $\lambda=1$ has $\rho_{\mathrm{LS}}\ge\varepsilon$, then $\rho_{\mathrm{LS}}(P^{\mathrm{GD}}(\nu)) \ge \varepsilon\,\exp\!\big(-2\|J\|_{\mathrm{op}}\!\int_0^1\!\alpha(\lambda)\,d\lambda\big)$.
- **Corollary 51 (Ising at $\|J\|_{\mathrm{op}}<1/2$)**: $\rho_{\mathrm{LS}}(P^{\mathrm{GD}}(\nu_{J,v})) \ge n^{-1}(1-2\|J\|_{\mathrm{op}})$, hence Glauber mixes in $O(n\log n)$ for any external field $v$. Recovers the SK result of Eldan–Koehler–Zeitouni / Anari–Jain–Koehler in particular at $\beta\le 1/4$.
- **Corollary 56 (graphical Ising in uniqueness)**: any $\nu_{G,\beta}$ satisfying $e^{|\beta|}<\frac{\Delta(G)-\delta}{\Delta(G)-2+\delta}$ has $\rho_{\mathrm{LS}}(P^{\mathrm{GD}}(\nu)) \ge \exp(-8/\delta)/n$, so $t_{\mathrm{mix}}=O(n\log n)$, **independent of the external field** (improves CFYZ21a).
- **Theorem 62 (hardcore in tree-uniqueness)**: for $\delta$-unique hardcore ($\lambda\le(1-\delta)\lambda_\Delta$ on graph of max degree $\Delta\ge 3$) and any initial $\mu_0$ on $\mathcal{I}_G$, $t_{\mathrm{mix}}(P^{\mathrm{GD}},\varepsilon;\mu_0) \le \exp(c/\delta)\,(n\log n + 3n\log(1/\varepsilon))$. **First optimal $O(n\log n)$ bound for the standard Glauber chain** (prior optima used a tailored balanced Glauber, AJK21a).
- **Corollary 54 (near-critical ferromagnetic Ising)**: $\rho_{\mathrm{LS}}(P^{\mathrm{GD}}(\nu_{J,v})) \ge n^{-1}\exp\!\big(-2\|J\|_{\mathrm{op}}\!\int_0^1\!\chi(\lambda)\,d\lambda\big)$ with $\chi(\lambda):=\|\mathrm{Cov}(\nu_{\lambda J,0})\|_{\mathrm{op}}$, recovering a variant of Bauerschmidt–Dagallier with control by the model's susceptibility up to $T_c$.
- **Theorem 58 (log-concave via Restricted Gaussian Dynamics)**: for $\mu$-strongly log-concave $\nu$ on $\mathbb{R}^n$, the RGD with parameter $\eta$ satisfies $\rho_{\mathrm{LS}}(P^{\mathrm{RGD}_\eta}(\nu)) \ge \mu/(\mu+1/\eta)$. Holds for any $\exp(O(n))$-warm start, resolving an open question of Lee–Shen–Tian.

## Method

To each measure $\nu$ assign a *localization process* $(\nu_t)_t$: a measure-valued martingale with $\nu_0=\nu$ that almost surely converges to a Dirac mass; the natural reverse Markov chain $P_{x\to y}=\mathbb{E}[\nu_\tau(x)\nu_\tau(y)/\nu(x)]$ recovers Glauber (coordinate-by-coordinate scheme), the up-down walk, RGD, etc. The spectral gap / MLSI of $P$ then reduces to *approximate conservation of variance / entropy* along the localization, which for *linear-tilt* localizations follows from boundedness of the covariance / influence matrices of the tilted measures (a Cauchy–Schwarz argument via the log-Laplace transform). Concatenating localizations gives an annealing principle that combines stochastic localization (to reduce $\|J\|_{\mathrm{op}}$) with coordinate-by-coordinate localization (to obtain Glauber gaps), bypassing high-dimensional-expander local-to-global theorems entirely.

## Quantum vs Classical

- **Baseline**: this paper *is* the modern classical ceiling for spin Gibbs sampling in the easy regimes — together with Anari–Liu–Oveis Gharan and Chen–Liu–Vigoda it pins the optimal $O(n\log n)$ Glauber mixing for hardcore and Ising whenever a uniqueness / spectral / log-concavity condition holds. Earlier bounds (spectral independence framework in HDX, Weitz FPTAS, Eldan–Koehler–Zeitouni) are recovered with simpler proofs and (for graphical Ising / hardcore) improved.
- **Gap**: any quantum Gibbs sampler attacking a *classical* Hamiltonian $H_{\mathrm{cl}}$ (diagonal in computational basis) in any of these regimes must beat $\widetilde{O}(n)$ Glauber steps on the relevant operator subroutine just to break even, with the expected $\Omega(\log(1/\varepsilon))$ block-encoding overhead working against any quantum win.
- **Source of the difference**: the localization-schemes argument is fundamentally about *commuting* statistics — covariance and influence of conditioned classical measures. Non-commuting Hamiltonians have no analogue of the coordinate-by-coordinate pinning $\mathcal{R}_u\nu$ that drives the Glauber result, so the framework gives nothing in the quantum regime.
- **Caveat**: the constants $e^{c/\delta}$ (hardcore) and $e^{-8/\delta}/n$ (graphical Ising) blow up at the uniqueness boundary; precise dependence on $\delta$ is not optimized. SK is covered only up to $\beta=1/4$; El Alaoui–Montanari–Sellke push to $\beta=1/2$ via explicit simulation of stochastic localization, which is a *different* quantum-relevant idea (algorithmic localization may be a quantum-algorithm primitive [CHECK]).

## Implications for Quantum Advantage

- **Regime cell**: high-T / intermediate-T (uniqueness regime) $\times$ classical / commuting $\times$ classical spin. Spans hardcore, Ising, log-concave continuous as a unified cell.
- **What this changes**: removes any residual upper-bound slack across the entire commuting-Hamiltonian / above-threshold cell; "classical = $O(n\log n)$ Glauber" is now the quantitative statement for this whole cell, not just for one chain or one regime.
- **Promising or not**: explicitly *not promising* for quantum advantage. This is a high-confidence "classical wins or ties" cell — quantum Gibbs samplers (Chen–Kastoryano–Brandão KMS, Ding et al.) reduce to a near-optimal classical chain modulo coherent block-encoding overhead. The quantum chapter's punchline cell (non-stoquastic $\times$ low-T) sits *outside* the reach of localization schemes by construction, which is the right reason to expect quantum advantage there: there is no martingale-of-classical-conditionings analogue when the Hamiltonian does not diagonalize in any product basis.

## Open Questions / Limitations

- Constants $\exp(c/\delta)$, $\exp(-8/\delta)$ at the uniqueness boundary are not tight; sharper $\mathrm{poly}(1/\delta)$ would close the gap to optimal classical.
- No quantum analogue: the framework rests on linear tilts and pinnings of classical measures; lifting to non-commuting Hamiltonians is open (and likely impossible without genuinely new ideas).
- Ising at $\|J\|_{\mathrm{op}}\in[1/2,1)$ (between the framework reach and SK criticality $\beta=1/2$) is handled by El Alaoui–Montanari–Sellke 2022 via explicit stochastic-localization simulation, not by the static MLSI bounds here.
- Log-concave bound assumes the Restricted Gaussian Oracle, which is itself non-trivial to implement; gives a relative not absolute mixing time.
- Doesn't address discrete distributions outside the linear-tilt class (e.g. matroid bases — handled by ALOG / AJK families separately).

## Connections

- **Unifies** Anari–Liu–Oveis Gharan (spectral independence) and Eldan (stochastic localization) under one martingale-of-measures lens; reproves CLV21 entropic independence without HDX local-to-global.
- **Improves** CFYZ21a graphical-Ising mixing by removing external-field dependence, and AJK21a hardcore optimal-mixing by working with the standard Glauber chain (not balanced Glauber).
- **Recovers** Bauerschmidt–Dagallier near-critical ferromagnetic-Ising LSI in terms of susceptibility (Cor. 54).
- **Recovers** Eldan–Koehler–Zeitouni / AJK21a SK bound at $\beta\le 1/4$ via the operator-norm corollary.
- **Resolves** the $\exp(O(n))$-warm-start question of Lee–Shen–Tian for log-concave RGD sampling.
- **Concurrent** with Chen–Feng–Yin–Zhang (CFYZ22), which independently proves the hardcore optimal bound by a different technique.
- Within the thesis grid: this is the *cross-cutting modern-mixing-framework anchor* in `classical-gibbs-sampling-corpus.md` §6, sitting above the Ising/hardcore high-T cell and feeding all the surrounding entries (ALOG20, CLV21, EKZ22, BD22). The "where can quantum first beat classical?" punchline cell is precisely the cell to which localization schemes do *not* extend — non-stoquastic $\times$ low-T — which is consistent with that being the strongest candidate for quantum advantage.
