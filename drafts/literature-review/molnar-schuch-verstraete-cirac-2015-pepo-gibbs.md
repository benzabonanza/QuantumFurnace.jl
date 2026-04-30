---
paper: "Molnar, Schuch, Verstraete & Cirac (2015)"
title: "Approximating Gibbs states of local Hamiltonians efficiently with projected entangled pair states"
arxiv: "1406.2973"
year: 2015
venue: "PRB 91:045138"
pdf: "n/a"

# Regime tags
temperature: [high-T, intermediate-T]
commutativity: [non-commuting]
locality: [geometric-local, 2D, any-d]
particle-statistics: [spin]
hamiltonian-models: [generic-local-spin, heisenberg, ising]

# Algorithmic paradigm
paradigm: [classical-tensor-network, pepo, cluster-expansion]
quantum-or-classical: [classical]

# Result type
result-type: [bond-dimension-bound, polynomial-time, structural]
key-scaling: "$D = (N/\\varepsilon)^{O(\\beta)}$ in any spatial dimension; $D = e^{O(\\log_2^2(N/\\varepsilon))}$ when $\\beta < O(\\log_2 N)$"

# Cross-links
related: ["czarnik-dziarmaga-2015-peps-finite-T", "kuwahara-alhambra-anshu-2021-thermal-area-law", "kuwahara-kato-brandao-2020-cmi-clustering", "bakshi-liu-moitra-tang-2024-high-T-unentangled", "hastings-2006-solving-gapped-locally"]
---

# Approximating Gibbs States of Local Hamiltonians Efficiently with PEPO — Molnar, Schuch, Verstraete & Cirac 2015

**One-sentence takeaway**: There exists a PEPO of bond dimension $D = (N/\varepsilon)^{O(\beta)}$ that approximates the Gibbs state $\rho_\beta = e^{-\beta H}/Z$ of *any* local Hamiltonian on *any* spatial dimension to error $\varepsilon$, providing the canonical structural certificate that classical tensor networks can in principle represent thermal states efficiently outside the low-$T$ regime.

## Setting

- **Hamiltonian**: arbitrary local quantum spin Hamiltonian $H = \sum_i h_i$ with bounded-norm $k$-local terms on a finite lattice $\Lambda \subset \mathbb{Z}^d$ of $N$ sites, in any spatial dimension $d \ge 1$.
- **Object**: the (unnormalized) Gibbs operator $e^{-\beta H}$, viewed as an operator on $(\mathbb{C}^q)^{\otimes N}$.
- **Ansatz**: a projected entangled pair operator (PEPO) — the natural lift of an MPO to general lattice geometry, with one tensor per site and bond indices of dimension $D$ along each lattice edge.
- **Question**: what is the smallest bond dimension $D = D(\beta, N, \varepsilon)$ such that some PEPO $\sigma_D$ satisfies $\|\rho_\beta - \sigma_D\| \le \varepsilon$ in a suitable operator norm?
- **Temperature regime**: any inverse temperature $\beta$, but with bounds that are *only useful* in the high-$T$ to intermediate-$T$ window — at low $T$ ($\beta \gg \log N$) the bond dimension blows up and the certificate becomes vacuous.

## Main Results

- **Polynomial bond-dimension bound (any-$T$, any-$d$)**: there exists a PEPO of bond dimension
  $$ D = \left(\frac{N}{\varepsilon}\right)^{O(\beta)} $$
  approximating $\rho_\beta$ to additive error $\varepsilon$, for any local Hamiltonian on any spatial-dimension lattice. **The exponent is linear in $\beta$**, so this is *quasi-polynomial in $N$* once $\beta = O(\log N)$ and *exponential in $\beta$* otherwise — the crucial structural fact for the regime map.
- **High-$T$ sharpening**: when $\beta < O(\log_2 N)$ the bond dimension can be improved to
  $$ D = \exp\!\big(O(\log_2^2(N/\varepsilon))\big), $$
  i.e. quasi-polynomial in $N/\varepsilon$. In this regime the PEPO is a *genuinely small* tensor network.
- **Method extends to ground states**: when the density of states grows at most polynomially with system size (e.g. gapped local Hamiltonians), the same construction yields a PEPS approximation of the ground state with the same bond-dimension scaling, by taking $\beta = \mathrm{poly}(N)$ low enough that thermal contamination is suppressed.
- **Existence theorem, not algorithm**: the result is a *certificate of existence* — the bond-dimension bound proves a PEPO with the stated $D$ exists, but does not guarantee that contracting it (computing observables) is efficient. PEPS / PEPO contraction is in general $\#$P-hard (Schuch-Wolf-Verstraete-Cirac 2007).

## Method

Two complementary constructions:

1. **Cluster / linked-cluster expansion (after Hastings 2006)**: write $e^{-\beta H} = \sum_C w(C)$ as a sum over clusters $C$ of lattice sites, where the weight of a cluster decays exponentially in its size at high temperature. Truncating at clusters of diameter $\ell$ introduces error controlled by the convergence radius of the expansion; each cluster is encoded locally on the PEPO, giving bond dimension $D = e^{O(\ell^d)}$. Choosing $\ell$ to balance error and cost yields the high-$T$ quasi-polynomial bound.
2. **Suzuki-Trotter / imaginary-time slicing**: split $\beta$ into $M$ Trotter slices $e^{-\beta H/M}$, each represented as a low-bond-dimension PEPO via local exponentiation. Multiplying $M$ such PEPO layers gives bond dimension $D \le D_0^M$. Optimizing $M$ versus the per-slice approximation error and applying compression gives the general $D = (N/\varepsilon)^{O(\beta)}$ bound across all $\beta$ and $d$.

The compression step uses singular-value-decomposition truncation on the bond indices, controlled by an operator-norm error bound that propagates linearly in the lattice size.

## Quantum vs Classical

- **This is the rigorous classical-side existence statement** for tensor-network thermal-state representations in $d \ge 2$. It is the higher-dimensional analogue of Hastings (2006) [1D] and the rigorous companion to the *algorithmic* heuristics of Czarnik-Dziarmaga (2015) [2D PEPS-imaginary-time]. The thesis Review chapter cites this paper as the structural certificate that "above some temperature, classical tensor networks suffice in principle in any dimension."
- **Baseline comparison**: there is no quantum baseline this paper compares against directly — it is purely classical and pre-dates most quantum Gibbs-sampling algorithms. The relevant comparison is downstream: any quantum Gibbs sampler claiming an advantage over classical methods in $d \ge 2$ must contend with the fact that $\rho_\beta$ is *known to admit* a polynomial-bond-dimension classical representation in the high-$T$ / intermediate-$T$ regime.
- **Gap analysis (where this leaves room for quantum)**:
  - **Existence vs computation**: the bound is *non-constructive* in the algorithmic sense — even though a small-$D$ PEPO exists, computing observables from it requires PEPO contraction, which is $\#$P-hard in the worst case (Schuch et al. 2007). Practical contraction relies on heuristic schemes (CTMRG, boundary-MPS) without rigorous guarantees beyond 1D. **This is exactly the gap a quantum Gibbs sampler can exploit**: a quantum algorithm prepares $\rho_\beta$ and lets one *measure* observables directly, sidestepping PEPO contraction entirely.
  - **Low-$T$ blow-up**: $D = \exp(O(\beta))$ becomes super-polynomial once $\beta \gg \log N$. The bound is essentially tight — a thermal area law with $O(\beta)$-scaling boundary entanglement is necessary in the worst case (Wolf-Verstraete-Hastings-Cirac 2008 area law). **This carves out the low-$T$ regime as a structural opening for quantum advantage** in higher dimensions, since classical PEPO/PEPS representations stop being polynomial-size *as a matter of representation*, not just contraction.
- **Source of the difference**: the bond-dimension bound captures the *thermal correlation length*. At high $T$ correlations decay exponentially $\Rightarrow$ small $D$ suffices; at low $T$ critical / ordered correlations require $D = \exp(O(\beta))$. Quantum samplers don't pay this price *in representation* — they store $\rho_\beta$ as a quantum state on $N$ qubits, dimension $q^N$ regardless of correlation length. They pay instead in *mixing time* of the dissipative dynamics, which is the regime where Chen-Eldan / Ding-Chen / KMS bounds do or do not bite.
- **Caveat**: classical tensor-network *contraction* in 2D is heuristic but empirically excellent on stoquastic / unfrustrated 2D models (Czarnik-Dziarmaga); this means the certificate of representation existence is matched by an effective practical algorithm in the cells where physics is also "easy." The genuine quantum-advantage candidate is precisely the cell where contraction heuristics also break down.

## Implications for Quantum Advantage

- **Regime cell**: covers the **high-$T$ + intermediate-$T$ × non-commuting × geometric-local × any-$d$ × spin** cell with a *representation* certificate. The complement — **low-$T$ × any-$d$** — is where this paper's bound goes super-polynomial and quantum advantage becomes structurally plausible.
- **What this changes**: this paper *closes* (in the upper-bound-of-classical sense) the high-$T$ to intermediate-$T$ regime in any dimension as a candidate for an asymptotic quantum advantage *coming from representation alone* — there exists a polynomial-size classical object encoding $\rho_\beta$. Any quantum advantage in this regime must come from (a) the contraction-vs-measurement gap (algorithmic, not structural), or (b) tighter constants / lower polynomial degrees than what classical heuristics achieve in practice.
- **Promising-or-not**: this paper *narrows* the regime where structural quantum advantage is plausible to **low-$T$ × $d \ge 2$**. Combined with Kuwahara-Alhambra-Anshu (2021) sharpening 1D to $\tilde O(\beta^{2/3})$, the structural quantum-advantage frontier sits in 2D / 3D quantum lattice models below the thermal correlation-length regime where tensor networks need $D = \exp(O(\beta))$.
- **Synthesis with cluster-expansion / area-law literature**: paired with Kuwahara-Kato-Brandão (2020) (CMI clustering above threshold $T$) and Bakshi-Liu-Moitra-Tang (2024) (high-$T$ Gibbs states are unentangled), this paper completes a triangle of high-$T$ structural results — *unentangled* (BLMT), *Markovian* (KKB), *PEPO-representable* (this paper) — all reinforcing that high-$T$ is structurally classical-easy.

## Open Questions / Limitations

- **No efficient contraction**: the bond-dimension bound does not imply efficient classical computation of $\langle O \rangle_\beta$. The headline scaling $(N/\varepsilon)^{O(\beta)}$ is the *representation* cost; the *evaluation* cost in 2D and beyond is governed by heuristic boundary-MPS / CTMRG schemes and remains $\#$P-hard in the worst case. This is the gap a quantum algorithm exploits.
- **Tightness in $\beta$**: the linear-in-$\beta$ exponent is essentially forced by the worst-case thermal area law — the result is sharp up to constants for generic local Hamiltonians [CHECK exact tightness statement in paper]. The Kuwahara-Alhambra-Anshu (2021) improvement to $\tilde O(\beta^{2/3})$ in 1D suggests room for sharpening in higher $d$ but no such result yet exists.
- **Implicit constants**: the $O(\beta)$ exponent has constants depending on the locality range $k$, lattice coordination, and local Hilbert-space dimension $q$ — these are not tracked explicitly. For practical $d=2$ models with $\beta \sim 1$ they may be order-unity and the bound informative; for $\beta \sim 10$ they are likely too loose.
- **Norm choice**: the approximation is in operator norm (or 1-norm — [CHECK]). Trace-distance error on the *normalized* Gibbs state $\rho_\beta = e^{-\beta H}/Z$ inherits the bound up to a partition-function-ratio factor that is bounded by the same cluster-expansion machinery in the high-$T$ regime, but this transfer is delicate at low $T$.
- **Ground-state corollary depends on density-of-states assumption**: the extension to ground states assumes polynomial growth of the density of states, which excludes systems with macroscopically degenerate ground manifolds (frustrated systems, topological order with low-energy edge modes).

## Connections

- **Czarnik-Dziarmaga (2015)** — *algorithmic 2D companion*. They demonstrate empirically that a PEPS-with-ancilla evolved in imaginary time gives a workable 2D thermal-state algorithm at $D = 2$, $M \le 32$ for the 2D quantum Ising model. Molnar-SVC supplies the rigorous upper-bound theorem; Czarnik-Dziarmaga supplies the practical algorithm. The two papers are the structural-and-empirical pair for 2D classical Gibbs sampling.
- **Hastings (2006)** — *1D ancestor*. Proves the 1D analogue: thermal MPO with bond dimension polynomial in $1/\varepsilon$ for *gapped* 1D local Hamiltonians. Molnar-SVC lifts the existence statement to *any* lattice and *any* (high-enough) $T$ via a different (cluster-expansion + Suzuki-Trotter) construction.
- **Kuwahara-Alhambra-Anshu (2021)** — *1D sharpening + algorithm*. Improves the thermal area law in 1D from $O(\beta)$ to $\tilde O(\beta^{2/3})$ and turns the existence statement into a *quasi-linear-time* classical algorithm. No analogue exists in 2D. The contrast — rigorous algorithm in 1D, only existence in 2D — is the structural reason 2D quantum Gibbs sampling is the next frontier.
- **Kuwahara-Kato-Brandão (2020)** — *high-$T$ structural twin*. Above a threshold $T$, conditional mutual information clusters exponentially $\Rightarrow$ Markov-network approximation. Combined with Molnar-SVC, the high-$T$ regime in any $d$ has *both* a Markov-network local-recovery structure *and* a PEPO existence proof.
- **Bakshi-Liu-Moitra-Tang (2024) "high-$T$ unentangled"** — *high-$T$ structural climax*. Proves high-$T$ Gibbs states are *unentangled* (separable). This subsumes Molnar-SVC's high-$T$ bound in the sense that an unentangled state needs $D = O(1)$ (product-of-states form), but only above the unentanglement threshold $\beta < \beta^*$, which is generically smaller than $O(\log N)$.
- **Schuch-Wolf-Verstraete-Cirac (2007) "computational complexity of PEPS"** — *the obstruction to making this constructive*. Proves PEPS contraction is $\#$P-hard, which is exactly what prevents Molnar-SVC's existence theorem from immediately becoming a classical algorithm in $d \ge 2$.
- **Thesis Review chapter mapping**: this paper anchors the *higher-dimensional classical structural certificate* row of the Tier-1 grid, complementing the 1D rigorous algorithmic results (Kuwahara-Alhambra-Anshu) and the 2D heuristic algorithm (Czarnik-Dziarmaga).

Sources:
- [arXiv:1406.2973 — Approximating Gibbs states of local Hamiltonians efficiently with projected entangled pair states](https://arxiv.org/abs/1406.2973)
- [PRB 91:045138 (publisher record)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.91.045138)
- [Dissemin record (open-access metadata)](https://dissem.in/p/13161729/approximating-gibbs-states-of-local-hamiltonians-efficiently-with-projected-entangled-pair-states/)
