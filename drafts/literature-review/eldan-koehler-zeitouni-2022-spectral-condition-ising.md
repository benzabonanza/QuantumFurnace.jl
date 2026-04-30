---
paper: "Eldan, Koehler, Zeitouni (2022)"
title: "A spectral condition for spectral gap: fast mixing in high-temperature Ising models"
arxiv: "2007.08200"
year: 2022
venue: "Probab. Theory Rel. Fields 182:1035-1051"
pdf: "supplementary-informations/classical-review/eldan-koehler-zeitouni-2022-spectral-condition-ising.pdf"

temperature: [high-T]
commutativity: [commuting]
locality: [dense, long-range, geometric-local]
particle-statistics: [spin]
hamiltonian-models: [Ising, Sherrington-Kirkpatrick]

paradigm: [classical-MCMC]
quantum-or-classical: [classical]

result-type: [mixing-time-upper, gap-bound, structural]
key-scaling: "Poincare: (1 - ||J||_op) Var_{nu_0}(phi) <= E_{nu_0}(phi, phi) for any 0 <= J <= Id. Continuous-time Glauber t_mix(eps) <= (1 - ||J||_op)^{-1} ((1 + 2||J||_op) n + 2|h|_1 + log(1/eps)); discrete-time t_mix = O(n^2 + |h|_1 n + n log(1/eps)) / (1 - ||J||_op). Covers Sherrington-Kirkpatrick at any constant beta = O(1) (Dobrushin only reaches beta = O(1/sqrt(n)))."

related: ["chen-eldan-2022-localization-schemes", "anari-liu-oveisgharan-2020-spectral-independence-hardcore", "chen-liu-vigoda-2021-optimal-mixing-glauber-hde", "bauerschmidt-dagallier-2024-near-critical-ising-lsi"]
---

# A spectral condition for spectral gap: fast mixing in high-temperature Ising models — Eldan, Koehler, Zeitouni 2022

**One-sentence takeaway**: Proves a dimension-free Poincaré inequality $(1 - \|J\|_{\mathrm{op}})\,\mathrm{Var}_{\nu_0}(\varphi) \le \mathcal{E}_{\nu_0}(\varphi,\varphi)$ for the Glauber dynamics of any Ising model on $\{\pm 1\}^n$ with positive-semidefinite interaction $J \preceq \mathrm{Id}$, via a stochastic-localization construction that decomposes $\nu_0$ into a mixture of *rank-one* Ising measures — yielding polynomial Glauber mixing whenever $\|J\|_{\mathrm{op}} < 1$, in particular for the Sherrington–Kirkpatrick model at any constant $\beta$ (Dobrushin only reaches $\beta = O(1/\sqrt{n})$).

## Setting

- **Target**: general Ising measure $\dfrac{d\nu_0}{d\mu}(x) = Z^{-1}\exp\!\big(\tfrac{1}{2}\langle x, Jx\rangle + \langle h, x\rangle\big)$ on the hypercube $\{\pm 1\}^n$ with arbitrary symmetric quadratic interaction $J$ and external field $h \in \mathbb{R}^n$; $\mu$ is the uniform measure.
- **Hamiltonian class**: classical (commuting) spin-$1/2$ Ising; *no* sparsity or geometric-locality restriction — $J$ is allowed dense / long-range / mean-field. WLOG $J$ is taken positive semidefinite (a constant shift on the diagonal absorbs into $h$).
- **Temperature regime**: high-T defined *spectrally* by $\|J\|_{\mathrm{op}} < 1$ (operator norm of $J$). Subsumes Dobrushin's $\|J\|_{\infty\to\infty} < 1$ (since $\|J\|_{\mathrm{op}} \le \|J\|_{\infty\to\infty}$) and is strictly weaker on dense matrices.
- **Access**: continuous-time Glauber semigroup $P_t = e^{t\Lambda}$ with generator $(\Lambda\varphi)(x) = \sum_i \big(\mathbb{E}_{\nu_0}[\varphi(X)\,|\,X_{\sim i} = x_{\sim i}] - \varphi(x)\big)$; classical MCMC.

## Main Results

- **Theorem 1 (Poincaré inequality, verbatim)**: for $\nu_0$ as above with $0 \preceq J \preceq \mathrm{Id}$ and any $\varphi : \{\pm 1\}^n \to \mathbb{R}$,
  $$ (1 - \|J\|_{\mathrm{OP}})\,\mathrm{Var}_{\nu_0}(\varphi(X)) \le \mathcal{E}_{\nu_0}(\varphi,\varphi). $$
  Hence the spectral gap $\gamma$ of continuous-time Glauber satisfies $\gamma \ge 1 - \|J\|_{\mathrm{op}}$.
- **Theorem 11 (mixing time, verbatim)**: for $P_t$ continuous-time Glauber on $\nu_0$,
  $$ \max_{x \in \{\pm 1\}^n} \|P_t(x,\cdot) - \nu_0\|_{\mathrm{TV}} \le \varepsilon \quad\text{as long as}\quad t \ge \frac{1}{1-\|J\|_{\mathrm{OP}}}\!\left((1+2\|J\|_{\mathrm{OP}})\,n + 2|h|_1 + \log\!\tfrac{1}{\varepsilon}\right). $$
  In discrete time (one Glauber step per Poisson tick of rate $n$) this gives $t_{\mathrm{mix}} = O\!\big(\tfrac{n^2 + |h|_1\,n + n\log(1/\varepsilon)}{1 - \|J\|_{\mathrm{op}}}\big)$.
- **Sherrington–Kirkpatrick application**: for $J_{ij} \sim \mathcal{N}(0, \beta^2/n)$ symmetric off-diagonal, $\|J\|_{\mathrm{op}} \to 2\beta$ a.s., so Glauber mixes in $\mathrm{poly}(n)$ for any constant $\beta < 1/2$ — Dobrushin only reached $\beta = O(1/\sqrt{n})$.
- **Section 4 — needle decomposition (structural)**: $\nu_0$ admits a measurable decomposition $\nu_0 = \int \nu_T\, d\mathbb{P}(T)$ as a mixture of *rank-one* Ising measures of the form $d\nu_T(x) \propto \exp(\tfrac{1}{2}\langle x, U\rangle^2 + \langle q + h, x\rangle)\,d\mu(x)$ with $|U|^2 \le \|J\|_{\mathrm{op}}$; rank-one analogue of the Kannan–Lovász–Simonovits needle decomposition for log-concave measures.

## Method

Construct an Itô process $\nu_t = F_t \cdot \nu_0$ (a stochastic localization) driven by an $n$-dim Brownian motion, with $J_t = J - \int_0^t C_s^2\, ds$ deterministically *decreasing* until $J_T \succeq 0$ has rank $\le 1$ at a stopping time $T \le \tfrac{1}{2}\mathrm{Tr}(J_0)$; the running matrix-valued coefficient $C_t$ is engineered (Lemma 2) so each step localizes one extra eigendirection of $J$. The Dirichlet form $\mathcal{E}_{\nu_t}(\varphi,\varphi)$ is shown to be a supermartingale (Lemma 9) along the process, so the Poincaré inequality at the rank-one terminal (proved via $\ell_2$-Dobrushin / Perron–Frobenius applied to the influence matrix, Lemmas 7–8) lifts to all of $\nu_0$ with constant $1 - \|J\|_{\mathrm{op}}$. The needle decomposition is precisely this stochastic-localization mixture, viewed structurally rather than dynamically.

## Quantum vs Classical

- **Baseline**: this *is* the modern classical ceiling at high-T for *general* (including dense / spin-glass / mean-field) Ising. The previous best general-purpose criterion was Dobrushin's $\|J\|_{\infty\to\infty} < 1$, which is far from tight on SK and on any dense $J$ where row-$\ell_1$ norms are $\Theta(\sqrt{n})$ times the operator norm. Bauerschmidt–Bodineau 2019 had given a log-Sobolev inequality at high-T for PSD $J$ but with a *different* (gradient-square) Dirichlet form — Eldan–Koehler–Zeitouni show that bound does *not* directly imply polynomial Glauber mixing (the two Dirichlet forms can differ by $e^{\Theta(\beta\sqrt{n})}$ on SK, App. A), and bridges the gap.
- **Gap**: there is no quantum candidate to beat in this cell — Ising on the hypercube is classical and commuting, so any quantum sampler reduces to (at best) a near-optimal classical chain modulo block-encoding overhead. The bar quantum Gibbs samplers must clear here is $\widetilde{O}(n^2/(1-\|J\|_{\mathrm{op}}))$ Glauber sweeps, with no known quantum win.
- **Source of the difference**: stochastic localization is a *commuting* construction — it tilts a classical product-measure base $\mu$ by a martingale of linear forms on the spin coordinates, with no analogue when the Hamiltonian fails to diagonalize in any product basis. Likewise the needle decomposition into rank-one Ising measures relies on the joint diagonalisability of $J$ with itself, not available for non-commuting $H$.
- **Caveat**: the dependence $1 - \|J\|_{\mathrm{op}}$ blows up at the edge $\|J\|_{\mathrm{op}} \to 1$; for SK this corresponds to $\beta \to 1/2$, well below the spin-glass critical $\beta = 1$ (Parisi). El Alaoui–Montanari–Sellke 2022 push to $\beta < 1/2$ via algorithmic stochastic-localization simulation; reaching the conjectured $\beta < 1$ is open. [CHECK] whether SK Glauber mixing is *expected* to be polynomial up to $\beta = 1$ or only up to $\beta = 1/2$.

## Implications for Quantum Advantage

- **Regime cell**: classical Ising $\times$ high-T (dense / mean-field / spin-glass slice). Specifically, the corpus's *classical-Ising / high-T* cell is now closed by this paper for all interactions with $\|J\|_{\mathrm{op}} < 1$, complemented by Anari–Liu–Oveis Gharan / Chen–Liu–Vigoda for sparse antiferromagnetic systems and Bauerschmidt–Dagallier for ferromagnetic Ising up to criticality.
- **What this changes**: tightens the classical baseline at high-T. *Any* claim of quantum advantage on a classical Ising problem at $\|J\|_{\mathrm{op}} < 1$ now competes against $\widetilde{O}(n^2)$ classical Glauber, regardless of how dense or unstructured $J$ is. The result also kills a naïve hope that "spin glasses are hard for classical at any $\beta$" — at constant $\beta$ they are not.
- **Promising or not** (for quantum advantage): explicitly *not promising*. This cell is a high-confidence "classical wins or ties" cell. The thesis chapter's punchline (non-stoquastic $\times$ low-T, where commutativity, sign-freeness, and spatial locality all fail) sits structurally outside the reach of stochastic localization.

## Open Questions / Limitations

- **SK regime**: result reaches $\beta < 1/2$; the conjectured Glauber-rapid-mixing threshold for SK is at the spin-glass transition $\beta = 1$. El Alaoui–Montanari–Sellke and Anari–Jain–Koehler–Pham–Vuong (entropic independence) push variants of the bound; the gap to $\beta = 1$ is open.
- **Functional inequality**: only a Poincaré (spectral-gap) bound is obtained — no log-Sobolev inequality, so the bound on TV mixing is $O(n^2)$ rather than the optimal $O(n \log n)$ that an MLSI would give. Chen–Eldan 2022 promotes this to MLSI / $O(n \log n)$ via the localization-schemes framework at the cost of $\|J\|_{\mathrm{op}} < 1/2$.
- **Asymmetry between Dirichlet forms**: $\mathcal{E}_{\nu_0}$ (Glauber) vs $\mathbb{E}|\nabla\varphi|^2$ (gradient form of Bauerschmidt–Bodineau) can differ by $e^{\Theta(\beta\sqrt{n})}$ on SK (App. A); a single inequality controlling both is open.
- **Antiferromagnetic / sign-indefinite $J$**: PSD assumption is essential to the localization construction. Sign-indefinite (antiferromagnetic) $J$ is reduced via a constant shift, but $\|J\|_{\mathrm{op}}$ can grow and tightness for hardcore-type models is lost relative to spectral-independence.
- **No quantum extension**: the localization construction is intrinsically classical; lifting to non-commuting Hamiltonians is open and likely requires genuinely new ideas.

## Connections

- **Bauerschmidt–Bodineau 2019**: gives a log-Sobolev inequality for the Glauber semigroup of PSD-$J$ Ising at $\|J\|_{\mathrm{op}} < 1$ but with the gradient-square Dirichlet form; Eldan–Koehler–Zeitouni show this LSI does *not* directly imply polynomial mixing on SK and bridge the gap with a true Glauber Poincaré inequality.
- **Chen–Eldan 2022** [`chen-eldan-2022-localization-schemes`]: subsumes and sharpens this paper. Recasts the rank-one decomposition as a *localization scheme* and lifts the bound to a modified-log-Sobolev inequality at $\|J\|_{\mathrm{op}} < 1/2$, recovering optimal $O(n\log n)$ mixing in this sub-regime.
- **Anari, Jain, Koehler, Pham, Vuong (entropic independence)**: parallel route — proves an MLSI for high-T Ising via fractional log-concavity / entropic independence with an explicit gap to the boundary, removing degree-dependence; covers SK at $\beta < 1/4$ before Chen–Eldan extends to $\beta < 1/2$.
- **Anari–Liu–Oveis Gharan 2020** [`anari-liu-oveisgharan-2020-spectral-independence-hardcore`] / **Chen–Liu–Vigoda 2021** [`chen-liu-vigoda-2021-optimal-mixing-glauber-hde`]: the *sparse / bounded-degree* counterpart — spectral / entropic independence in high-dimensional expanders gives $O(n \log n)$ Glauber mixing for hardcore and antiferromagnetic 2-spin up to the tree-uniqueness threshold. Eldan–Koehler–Zeitouni handles dense / mean-field $J$ where bounded-degree HDX arguments do not apply, so the two are *complementary* covers of the classical-Ising / high-T cell.
- **Bauerschmidt–Dagallier 2024** [`bauerschmidt-dagallier-2024-near-critical-ising-lsi`]: pushes the *ferromagnetic* Ising LSI uniformly up to $\beta_c$ (under mean-field-bound on susceptibility), where Eldan–Koehler–Zeitouni stops at $\|J\|_{\mathrm{op}} < 1$; tradeoff is generality of $J$ (Eldan–Koehler–Zeitouni: arbitrary symmetric, including spin-glass; Bauerschmidt–Dagallier: ferromagnetic only).
- **Eldan 2013** (stochastic localization for KLS): the continuous progenitor of the discrete construction here; same martingale-of-tilted-measures idea, applied in $\mathbb{R}^n$ for the KLS conjecture.
- **Kannan–Lovász–Simonovits needle decomposition**: §4 explicitly frames its rank-one mixture as a discrete analogue of the KLS needle decomposition for log-concave measures.
- Within the thesis grid: this is the *Tier-1 anchor for classical Ising × high-T (dense/mean-field slice)* in `classical-gibbs-sampling-corpus.md` §2 / §6, sitting alongside Anari–Liu–Oveis Gharan (sparse slice) and Chen–Liu–Vigoda (optimal $n\log n$). Together with Chen–Eldan it sets the floor every quantum Ising sampler at high-T must beat — and the chapter's quantum-advantage candidate cells (non-stoquastic $\times$ intermediate-/low-T) sit precisely where this framework provides no analogue.
