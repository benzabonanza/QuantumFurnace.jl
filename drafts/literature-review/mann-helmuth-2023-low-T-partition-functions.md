---
paper: "Mann, Helmuth (2023)"
title: "Efficient Algorithms for Approximating Quantum Partition Functions at Low Temperature"
arxiv: "2201.06533"
year: 2023
venue: "Quantum 7:1155"
pdf: "n/a"

temperature: [low-T]
commutativity: [non-commuting]
locality: [geometric-local, any-d, sparse]
particle-statistics: [spin]
hamiltonian-models: [stable-quantum-spin, transverse-ising, perturbation-of-classical]

paradigm: [classical-deterministic, cluster-expansion, polymer-model, contour-expansion]
quantum-or-classical: [classical]

result-type: [FPTAS, polynomial-time-deterministic]
key-scaling: "FPTAS for $\\log Z(\\beta,\\lambda)$ in time $|V(G)|^{O(1)}\\,(1/\\varepsilon)^{O(1)}$ on bounded-degree $\\mathbb{Z}^d$ subgraphs at all $\\beta\\ge\\beta_\\star$, $|\\lambda|\\le\\lambda_\\star$ (model-dependent constants from Peierls $\\alpha_0$, $d$, $|\\Xi|$)"

related: ["bovier-den-hollander-2015-metastability", "mann-helmuth-2021-high-T-partition-functions", "kuwahara-alhambra-anshu-2021-thermal-area-law", "bakshi-liu-moitra-tang-2024-high-T-unentangled"]
---

# Efficient Algorithms for Approximating Quantum Partition Functions at Low Temperature — Mann, Helmuth (2023)

**One-sentence takeaway**: For *stable* quantum perturbations of classical lattice spin systems with finitely many discrete ground states, low-temperature partition functions admit a deterministic FPTAS — closing the low-T cell of the regime grid for the class where a Peierls condition survives the quantum perturbation, while leaving the genuinely hard punchline cells (criticality, frustration, antiferromagnets without symmetry breaking) untouched.

## Setting

- **Object computed**: $\log Z(\beta,\lambda) = \log \mathrm{Tr}\,e^{-\beta H}$ for an $n$-vertex quantum spin Hamiltonian on a bounded-degree subgraph $G \subset \mathbb{Z}^d$, $d\ge 2$.
- **Hamiltonian class**: $H = H_\Phi + \lambda H_\Psi$, where
  - $H_\Phi$ is *diagonal* in a fixed product (classical) spin basis, with finite ground-state set $\Xi$ and a **Peierls condition**: every vertex in an excited configuration pays at least $\alpha_0 > 0$ above the ground-state energy ($\Phi_s(v) \ge e_0 + \alpha_0$ for $s \notin \Xi$).
  - $H_\Psi = \sum_e \Psi(e)$ is a sum of bounded edge operators ($\|\Psi(e)\|\le 1$), generically *non-commuting* with $H_\Phi$ — the genuine quantum perturbation.
  - $\lambda$ is the perturbation strength, controlled by a model-dependent threshold $\lambda_\star$.
- **Temperature regime**: low-T, $\beta \ge \beta_\star = \beta_\star(d,|\Xi|,\alpha_0)$. The thresholds $\beta_\star, \lambda_\star$ are *uniform in system size* — what's required is that the Borgs–Kotecký–Ueltschi quantum Pirogov–Sinai condition holds, i.e. the classical zero-temperature phase diagram survives the quantum perturbation.
- **Access model**: classical input — adjacency list of $G$, classical specifications of $\Phi$ and $\Psi$. Classical output — an $\varepsilon$-approximation to $\log Z$.

## Main Results

- **Theorem 15 (FPTAS for stable quantum spin systems)**: there exist $\beta_\star, \lambda_\star$ depending on $d$, $|\Xi|$, $\alpha_0$, and the maximum degree such that for all $\beta \ge \beta_\star$ and $|\lambda|\le\lambda_\star$, a deterministic algorithm computes $\widehat{F}$ with $|\widehat{F} - \log Z(\beta,\lambda)| \le \varepsilon$ in time
  $$ |V(G)|^{O(1)} \cdot (1/\varepsilon)^{O(1)}. $$
  The exponents depend polynomially on $d$ and on the truncation depth $n = O(\log(|V(G)|/\varepsilon))$ via $\exp(O(n))$ contour-enumeration cost, but $n$ is logarithmic so the overall scaling is polynomial in $|V(G)|$ and $1/\varepsilon$. [CHECK — exact dimension dependence of the $O(1)$ exponent].
- **Sampling**: under the same hypotheses, an FPRAS for the Gibbs measure restricted to a fixed dominant phase (one of the $|\Xi|$ ground-state sectors) follows by self-reducibility from the partition-function FPTAS [CHECK — exact statement].
- The result extends the *high-T* deterministic-FPTAS line of Mann–Helmuth (2021) [arXiv:2004.11568] — partition function in time $|V(G)|^{O(1)}(1/\varepsilon)^{O(1)}$ above a model-dependent $\beta_c$ — to *all* $\beta \ge \beta_\star$, leaving only the *intermediate-T* window (between rapid-mixing high-T and stable low-T) uncovered.

## Method

Quantum Pirogov–Sinai. Apply the **Borgs–Kotecký–Ueltschi 1996** Duhamel/path-integral expansion of $\mathrm{Tr}\,e^{-\beta H}$ as a sum over *contours* in $d{+}1$ dimensions: each contour is a connected space-time region where the quantum perturbation has acted, embedded in a classical-ground-state background. Stability of the classical ground states under the quantum perturbation gives an exponential weight bound on contours ($e^{-\alpha_0 |\gamma|}$ Peierls cost). This converts the quantum partition function into an *abstract polymer model* with controllable activities. The algorithmic step is then **Helmuth–Perkins–Regts 2019** / **Borgs–Chayes–Kolla–Helmuth–Perkins 2020** style: enumerate contours up to size $n = O(\log(|V|/\varepsilon))$ in time $|V|\cdot e^{O(n)}$, evaluate their weights via truncated Dyson series, and sum the cluster expansion. The non-trivial new ingredient over the classical case is showing that *quantum* contour weights — operator-valued objects — can still be computed deterministically to additive accuracy in time exponential only in the contour size, not in $|V|$.

## Quantum vs Classical

- **Baseline (classical)**: Helmuth–Perkins–Regts (2019) and Borgs–Chayes–Kolla–Helmuth–Perkins (2020) give deterministic FPTAS for *classical* low-T Ising / Potts / hard-core via algorithmic Pirogov–Sinai; this paper recovers the *same* polynomial scaling for genuinely quantum (non-commuting) perturbations on top of the same classical scaffolding. Bovier–den Hollander (2015) gives the matching $\Theta(e^{\beta\Gamma^*})$ classical *dynamic* (Glauber) lower bound — orthogonal complement: this paper does not run a Markov chain, it computes $\log Z$ directly.
- **Baseline (quantum)**: Chen–Kastoryano–Gilyén KMS Lindbladians and Ding et al. low-T quantum samplers achieve $\mathrm{poly}(n,\beta,1/\varepsilon)$ for low-T Gibbs states under their own (different) structural hypotheses — but for the *partition-function* problem on stable-perturbation Hamiltonians, the *classical* algorithm here matches them. No quantum advantage is left in this sub-cell.
- **Gap**: **none** — both quantum and classical achieve $\mathrm{poly}(n, 1/\varepsilon)$ at fixed $(\beta, \lambda)$ inside the stability region. The gap is *closed* by the classical side.
- **Source of the difference (or lack thereof)**: the structural reason a stable quantum perturbation of a classical model is classically tractable is exactly the one that makes it *uninteresting* for quantum advantage: the partition function is dominated by an $|\Xi|$-fold sum over classical ground states, with quantum fluctuations confined to small space-time droplets controlled by Peierls. The non-commutativity is *bounded* by the contour expansion. There is no sign problem because $H_\Phi$ is diagonal and $H_\Psi$ is treated perturbatively: contour weights are real and non-negative after the Duhamel expansion under standard ferromagnetic/stoquastic-perturbation assumptions [CHECK — the paper handles a general perturbation but the sign-positivity of weights requires care].
- **Caveat**: the constants $\beta_\star, \lambda_\star$ are explicit but quantitatively far from any physical phase boundary; the FPTAS is asymptotic in $|V|$, not in $\lambda \to \lambda_c$. Critical regimes ($\beta \to \beta_c$, $\lambda \to \lambda_\star$) are not covered.

## Implications for Quantum Advantage

- **Regime cell**: low-T $\times$ non-commuting $\times$ geometric-local $\times$ spin, restricted to **stable quantum perturbations of classical models with discrete symmetry breaking** (e.g. transverse-field ferromagnetic Ising on $\mathbb{Z}^d$ in the ordered phase at small transverse field, ferromagnetic quantum Potts, classical Ising plus weak quantum off-diagonal hopping).
- **What this changes**: **closes a sub-cell**. The low-T row of the regime grid was the strongest candidate for quantum advantage on structural grounds (Bovier–den Hollander $e^{\beta\Gamma^*}$ classical Glauber barrier, Troyer–Wiese sign problem). This paper shows that *partition-function* approximation in the entire stable-perturbation sub-cell is classically poly-time — so the advantage frontier inside the low-T row is **strictly outside** the stability region.
- **Promising or not**: **negative for advantage** in this sub-cell, but **sharpens** the frontier elsewhere. The remaining low-T candidates for quantum speedup are exactly:
  1. **Criticality / phase boundary** ($\lambda \uparrow \lambda_\star$, $\beta \downarrow \beta_\star$) — outside Pirogov–Sinai control;
  2. **Frustrated / antiferromagnetic** Hamiltonians without a stable classical ground-state structure — Peierls condition fails;
  3. **Continuous symmetry breaking** (Heisenberg ferromagnet, $U(1)$) — explicitly excluded by the authors due to slow-decaying correlations;
  4. **Non-stoquastic** systems where the perturbation $H_\Psi$ creates genuinely complex contour weights — the punchline cell of the corpus.
- The companion *high-T* paper (Mann–Helmuth 2021) closes the high-T row for the same model class via the cluster expansion. Together they leave only the *intermediate-T* window (in stable-perturbation models) and the four hard cells above as places where quantum could first beat classical for partition-function approximation.

## Open Questions / Limitations

- **Stability hypothesis is restrictive**: requires a discrete-symmetry-broken classical limit with a quantitative Peierls bound; rules out Heisenberg-type models, frustrated antiferromagnets on triangular/kagome lattices, and any model whose ground-state degeneracy is not finite/discrete.
- **Critical / coexistence regime open**: the algorithm degrades as $\lambda \to \lambda_\star$ or $\beta \to \beta_\star$; the Pirogov–Sinai constants blow up at the phase boundary. Whether quantum methods can extract polynomial-time bounds *near* criticality on these models is open.
- **Dimension dependence of the exponent**: $|V|^{O(1)}\cdot(1/\varepsilon)^{O(1)}$ has $d$-dependent constants from contour enumeration; the $d \to \infty$ behaviour is unaddressed [CHECK].
- **Sampling vs counting**: the partition-function FPTAS yields an FPRAS for the Gibbs state restricted to a single phase, but sampling *across* the $|\Xi|$ symmetry-broken sectors at low T inherits the classical ergodicity barrier — exactly the regime where Bovier–den Hollander's $e^{\beta\Gamma^*}$ time applies to local dynamics.
- **Non-stoquastic perturbations**: the contour-weight analysis does not require $H_\Psi$ to be stoquastic but the absence of catastrophic sign cancellation in the Dyson expansion deserves explicit checking for a non-stoquastic example [CHECK].

## Connections

- **Companion high-T paper**: `mann-helmuth-2021-high-T-partition-functions` — same authors, same algorithmic strategy on the *cluster* (rather than *contour*) expansion; together they cover both temperature extremes of the stable-perturbation class.
- **Classical Pirogov–Sinai algorithmic line**: Helmuth–Perkins–Regts 2019 and Borgs–Chayes–Kolla–Helmuth–Perkins 2020 — the classical scaffolding; this paper lifts their machinery to the quantum setting.
- **Quantum Pirogov–Sinai foundation**: Borgs–Kotecký–Ueltschi (J. Stat. Phys. 1996; CMP 1996) — the contour representation in $d{+}1$ dimensions, this paper's key analytical input.
- **Classical low-T baseline**: `bovier-den-hollander-2015-metastability` — the *dynamic* counterpart at low T (sharp Eyring–Kramers $e^{\beta\Gamma^*}$ for Glauber). Together with Mann–Helmuth this gives "$\log Z$ is easy, but local sampling between phases is hard" — the structural picture for the entire low-T row in stable-perturbation models.
- **Thermal area law**: `kuwahara-alhambra-anshu-2021-thermal-area-law` — orthogonal structural result for any-T 1D Gibbs states; Mann–Helmuth complements this for $d\ge 2$ in a restricted Hamiltonian class.
- **High-T closure**: `bakshi-liu-moitra-tang-2024-high-T-unentangled` — closes the high-T cell *unconditionally* (any local Hamiltonian, not just stable perturbations) via separability; Mann–Helmuth 2023 plays the analogous role at low-T, but only for the stability sub-cell.
- **Punchline cell beyond this paper's reach**: the Review chapter's quantum-advantage frontier lives at non-stoquastic $\times$ low-T or $\times$ intermediate-T (Klassen–Marvian sign-curing NP-hardness, Troyer–Wiese sign problem) — Mann–Helmuth 2023 does *not* enter this cell, which is the central observation for the chapter's framing.
