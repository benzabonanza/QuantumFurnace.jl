---
paper: "Bravyi and Gosset (2017)"
title: "Polynomial-time classical simulation of quantum ferromagnets"
arxiv: "1612.05602"
year: 2017
venue: "PRL"
pdf: "supplementary-informations/classical-review/bravyi-gosset-2017-quantum-ferromagnets.pdf"

temperature: [any-finite-T, ground-state-via-low-T-limit]
commutativity: [stoquastic, non-commuting]
locality: [2-local, geometrically-local-and-non-local, any-graph]
particle-statistics: [spin]
hamiltonian-models: [transverse-field-Ising, XY-ferromagnet, Heisenberg-ferromagnet]

paradigm: [classical-other]
quantum-or-classical: [classical]

result-type: [complexity-upper, FPRAS, partition-function-approximation]
key-scaling: "FPRAS for $\\mathcal{Z}(\\beta,H)$ with runtime $\\mathrm{poly}(n,\\beta,\\varepsilon^{-1})$; explicitly $\\widetilde{O}(n^{115}(1+\\beta^{46})\\varepsilon^{-25})$"

related: ["Jerrum-Sinclair-Vigoda 2004 (permanent FPRAS)", "Jerrum-Sinclair 1993 (ferro-Ising FPRAS)", "Bravyi-DiVincenzo-Oliveira-Terhal 2008 (stoquastic-LH)", "Bravyi-Terhal 2010 (frustration-free stoquastic)", "Bravyi 2015 (Monte Carlo for stoquastic)", "Mann-Helmuth 2021/2023 (polymer-based quantum partition functions)", "Heilmann-Lieb 1972 (matching log-concavity)"]
---

# Polynomial-time classical simulation of quantum ferromagnets — Bravyi & Gosset 2017

**One-sentence takeaway**: For an $n$-qubit *ferromagnetic* XY/Heisenberg model with optional transverse field on *any* graph, the partition function $\mathcal{Z}(\beta,H)$ admits a classical FPRAS with $\mathrm{poly}(n,\beta,\varepsilon^{-1})$ runtime — closing the "stoquastic-ferromagnetic-XY" cell as classically tractable at *all* finite temperatures.

## Setting

- **Hamiltonian family** (Eq. 1):
  $$ H = \sum_{1\le i<j\le n} \bigl(-b_{ij} X_i X_j + c_{ij} Y_i Y_j\bigr) + \sum_{i=1}^n d_i (I+Z_i), $$
  with real coefficients $b_{ij}\ge 0$, $|c_{ij}|\le b_{ij}$, $|d_i|\le 1$.
- **Special cases**: ferromagnetic XY ($c_{ij} = -b_{ij}$), transverse-field ferromagnetic Ising ($c_{ij}=0$), Heisenberg ferromagnet (uniform), plus a continuum interpolating between them.
- **Stoquastic**: each off-diagonal local term has non-positive entries in the computational basis (sign-problem-free).
- **Graph**: arbitrary — no planarity, bipartiteness, or geometric locality required.
- **Task**: approximate $\mathcal{Z}(\beta,H) = \mathrm{Tr}\,e^{-\beta H}$ within relative error $\varepsilon$, hence free energy $\mathcal{F}(\beta) = -\beta^{-1}\log\mathcal{Z}$ to additive error and (by $\beta = O(n\Delta^{-1})$) ground-state energy to additive error $\Delta$.
- **Temperature regime**: any $\beta>0$ (cost polynomial in $\beta$); ground state is the low-$T$ limit, costing $\mathrm{poly}(n,\Delta^{-1})$ for additive accuracy $\Delta$.

## Main Results

- **Theorem 1** (verbatim): $\mathcal{Z}(\beta,H)$ admits a randomized approximation scheme with runtime $\mathrm{poly}(n,\beta,\varepsilon^{-1})$ (relative error $\varepsilon$, success probability $\ge 3/4$).
- **Explicit scaling** (from the proof, p. 3): runtime $\widetilde{O}\bigl(n^{115}(1+\beta^{46})\varepsilon^{-25}\bigr)$ — provably efficient but far from practical.
- **Ground-energy corollary**: setting $\beta=O(n\Delta^{-1})$ yields a $\mathrm{poly}(n,\Delta^{-1})$ classical algorithm for the ground energy to additive error $\Delta$, since $0\le E_0 - \mathcal{F}(\beta)\le nT$.
- **Boundary settled**: this resolves the complexity of a one-parameter family of local Hamiltonian problems studied by Piddock-Montanaro [arXiv:1506.04014] — putting them on the *easy* side.

## Method

A three-step quantum-to-classical reduction. (i) **Suzuki–Trotter discretization** (Lemma 1): approximate $e^{-\beta H}$ by a product of $J = O(n^5(\beta^2+1)\varepsilon^{-1})$ elementary 1- and 2-qubit gates from a fixed set $\mathcal{G}$ whose matrix elements are non-negative reals — the Magnus expansion controls the third-order Trotter error explicitly. (ii) **Gadget mapping**: each gate is realized by a small weighted graph "gadget" so that $\mathcal{Z}_J = \mathrm{Tr}[G_J\cdots G_1] = \mathrm{PerfMatch}(\Gamma)$ for an assembled graph $\Gamma$ with $|V|=O(J)$, $|E|=O(J)$, $w_{\max}\le 2$. (iii) **Permanent FPRAS**: the resulting $\Gamma$ satisfies $\mathrm{NearPerfMatch}(\Gamma)/\mathrm{PerfMatch}(\Gamma) \le O(|V|^2)$, the structural condition required by the Jerrum–Sinclair–Vigoda permanent FPRAS [J. ACM 51, 2004], which then estimates $\mathrm{PerfMatch}(\Gamma)$ in $\widetilde{O}(\varepsilon^{-2}|V|^6|E|^5 w_{\max}^6 q(|V|)^4)$ time.

## Quantum vs Classical

- **Baseline**: prior to this paper, no general classical algorithm beat exact diagonalization for quantum ferromagnets *with* a transverse field; quantum Monte Carlo (SSE, world-line, loop) was empirically successful on bipartite ferromagnets but had no rigorous polynomial-time mixing guarantee. The naive route — sample world-line configurations directly — runs into local autocorrelation issues that no general theorem controlled.
- **Gap**: this work gives the *first* unconditional polynomial-time guarantee for the *finite-temperature* partition function of stoquastic transverse-field XY ferromagnets on arbitrary graphs. Versus quantum Gibbs samplers (Chen et al. KMS Lindbladians, Ding et al.), classical wins outright in this cell — there is no quantum-advantage opening at any temperature.
- **Source of the difference**: ferromagnetic XY/Heisenberg interactions in the standard basis produce *non-negative* off-diagonal matrix elements (stoquastic with the right sign). Combined with the structural fact that the resulting Trotter graph has bounded $\mathrm{NearPerfMatch}/\mathrm{PerfMatch}$ ratio, this exposes a path-integral representation that is not just sign-problem-free but also has a *ratio-bounded* matching structure — exactly the JSV condition for permanent FPRAS. Generic stoquastic Hamiltonians lack the second property.
- **Caveat**: the polynomial is enormous ($n^{115}\beta^{46}\varepsilon^{-25}$), so the result is a *complexity-theoretic* statement, not a practical replacement for SSE/world-line QMC. It is also restricted to ferromagnetic XY-type couplings (with magnitude condition $|c_{ij}|\le b_{ij}$); general stoquastic Hamiltonians remain open.

## Implications for Quantum Advantage

- **Regime cell**: stoquastic × any-finite-T × ferromagnetic-XY/Heisenberg/TFI, on *any* graph — the cell from Tier-1 corpus row "Stoquastic / intermediate-near-critical".
- **What this changes**: closes that cell as classically *provably* polynomial. Even at criticality of the underlying classical phase transition, the quantum algorithm cannot beat classical asymptotically here. Combined with Bravyi-Terhal (frustration-free stoquastic) and Bravyi (2015, Monte Carlo for stoquastic), it tightens the perimeter of "easy classical Hamiltonians" within the stoquastic class.
- **Promising or not**: *not* a quantum-advantage venue. Pushes attention toward (a) *non-stoquastic* quantum spin systems, (b) stoquastic but *non-ferromagnetic* (e.g. antiferromagnetic XY on non-bipartite / frustrated graphs — sign-problem-free becomes representation-dependent), (c) low-T physics where the gap closes faster than $\mathrm{poly}(n,\beta)$ permits — where the $\beta^{46}$ dependence may eventually bite and a quantum sampler with mild $\beta$-dependence could win.

## Open Questions / Limitations

- **Tighten the polynomial**: the $n^{115}\beta^{46}$ scaling is loose; the authors flag this as outstanding. Improvements to JSV (McQuillan 2013; Huang-Lu-Zhang 2016) might extend or sharpen the result.
- **Beyond ferromagnetic XY**: do the techniques extend to (i) non-bipartite ferromagnetic Heisenberg with frustrating local fields, (ii) signed XY where $|c_{ij}|>b_{ij}$ in some bonds, (iii) higher-spin generalizations? The $|c_{ij}|\le b_{ij}$ bound is essential to non-negativity of the gadget weights.
- **Antiferromagnetic case**: the corresponding antiferromagnetic XY/Heisenberg on bipartite graphs *is* stoquastic after a sublattice rotation — does that fall in scope? Not addressed here. Frustrated antiferromagnetic Heisenberg is *not* stoquastic and is a candidate quantum-advantage cell. [CHECK]
- **Critical scaling**: behaviour of the algorithm at the quantum/thermal critical point of the transverse-field Ising model is not analyzed — the polynomial is uniform in $\beta\le \mathrm{poly}(n)$, so criticality is folded in, but the exponents are presumably suboptimal.

## Connections

- **Direct ancestor**: Jerrum-Sinclair (1993) gave the FPRAS for the *classical* ferromagnetic Ising partition function via subgraph-world MCMC. This paper is the quantum/Trotter generalization through a different combinatorial reduction (matchings, not subgraphs).
- **Engine**: Jerrum-Sinclair-Vigoda (2004) permanent FPRAS for $0/1$ matrices is the workhorse — its bipartite-permanent algorithm extends here via the ratio-bounded condition that the gadget construction provably establishes.
- **Stoquastic complexity class**: Bravyi-DiVincenzo-Oliveira-Terhal (2008) defined stoquastic-LH and placed it in AM ∩ MA-hard; Bravyi-Terhal (2010) handled frustration-free stoquastic adiabatic evolution; this paper extends *both* to a finite-temperature partition-function setting for the structured ferromagnetic subclass.
- **Cluster QMC analogue**: structurally cousin to Edwards-Sokal / Swendsen-Wang random-cluster representations — both extract a non-negative auxiliary measure from a sign-problem-free quantum / classical Hamiltonian, but here the dual object is matchings rather than FK clusters, and the combinatorial input is JSV rather than Wolff/Swendsen-Wang.
- **High-T quantum partition functions**: Mann-Helmuth (2021, 2023) handle *any* local Hamiltonian via abstract polymer expansions but only *above a threshold temperature*. Bravyi-Gosset is incomparable — narrower Hamiltonian class but every $\beta$.
- **Quantum-side counterpart**: the Chen et al. 2025 noncommutative Gibbs sampler and Ding et al. 2024 KMS samplers are *quantum* algorithms for general local $H$ with $\mathrm{poly}(n,\beta,\varepsilon^{-1})$ scaling under spectral-gap assumptions; on this paper's Hamiltonian class they offer no asymptotic advantage.
- **Heilmann-Lieb (1972)** matching log-concavity (Theorem 4 of the appendix) is the structural fact that lets the algorithm walk along $Z_k(\Gamma)$ ratios safely.
