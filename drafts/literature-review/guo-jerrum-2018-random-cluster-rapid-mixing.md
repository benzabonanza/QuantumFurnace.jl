---
paper: "Guo and Jerrum (2018)"
title: "Random cluster dynamics for the Ising model is rapidly mixing"
arxiv: "1605.00139"
year: 2018
venue: "Ann. Appl. Probab. 28:1292-1313"
pdf: "supplementary-informations/classical-review/guo-jerrum-2018-random-cluster-rapid-mixing.pdf"

temperature: [all-T]
commutativity: [commuting]
locality: [arbitrary-graph]
particle-statistics: [spin]
hamiltonian-models: [Ising, random-cluster, Potts-q2]

paradigm: [classical-cluster, classical-MCMC]
quantum-or-classical: [classical]

result-type: [mixing-time-upper, FPRAS]
key-scaling: "tau_eps(P_RC) <= 8 n^4 m^2 (m ln(1-p)^{-1} + ln eps^{-1}) for q=2 random-cluster Glauber on any graph with n=|V|, m=|E|, all 0<p<1; transfers polynomially to Swendsen-Wang for ferromagnetic Ising at any beta>0"

related: ["edwards-sokal-1988-fk-representation", "lubetzky-sly-2012-critical-ising-polynomial", "lubetzky-sly-2013-cutoff-ising-lattice", "bauerschmidt-dagallier-2024-near-critical-ising-lsi", "bravyi-gosset-2017-quantum-ferromagnets"]
---

# Random cluster dynamics for the Ising model is rapidly mixing — Guo & Jerrum 2018

**One-sentence takeaway**: Single edge-flip Glauber dynamics on the $q=2$ random-cluster (Fortuin-Kasteleyn) representation mixes in time $O(n^4 m^3)$ on *any* graph at *any* temperature, which by Ullrich's comparison gives the first polynomial mixing-time upper bound for Swendsen-Wang dynamics on ferromagnetic Ising — closing the ferromagnetic-Ising $\times$ all-$T$ cell as classically tractable everywhere with a *dynamic* (Markov-chain) algorithm.

## Setting

- **Sampled object**: random-cluster measure $\pi_{\text{RC};p,q}(S) \propto p^{|S|}(1-p)^{|E\setminus S|}\, q^{\kappa(S)}$ on $S \subseteq E$, where $\kappa(S)$ is the number of connected components of $(V,S)$. At $q=2$ and $p = 1 - 1/\beta$ this is the FK representation of the ferromagnetic Ising model with parameter $\beta > 1$.
- **Hamiltonian class**: ferromagnetic Ising on an arbitrary undirected graph $G=(V,E)$, $n = |V|$, $m = |E|$. No planarity, bipartiteness, geometric locality, or bounded-degree assumption.
- **Temperature regime**: all $\beta > 1$ (equivalently all $0 < p < 1$ with $q=2$). Single result spans high-$T$, near-critical, and low-$T$ uniformly.
- **Access model**: classical MCMC. The dynamics $P_{RC}$ is the *lazy single edge-flip* chain with Metropolis acceptance: pick $e \in E$ uniformly, propose $S \mapsto S \oplus \{e\}$, accept with probability $\min\{1, \pi_{RC}(S \oplus \{e\})/\pi_{RC}(S)\}$ (Eq. (3) of the paper).

## Main Results

- **Theorem 7 (verbatim)**: For the random cluster model with parameters $0 < p < 1$ and $q=2$,
  $$\tau_\varepsilon(P_{RC}) \le 8\, n^4 m^2 \bigl(m \ln(1-p)^{-1} + \ln \varepsilon^{-1}\bigr).$$
  Polynomial in $n, m$ for any fixed $p$; the $\ln(1-p)^{-1}$ factor blows up only logarithmically as $p \to 1$ (low temperature).
- **Lemma 6 (congestion bound)**: $\varrho(\Gamma_{RC}) \le 8 m^2 n^4$ for the canonical-paths flow constructed in Section 4.
- **Corollary (Swendsen-Wang)**: by Ullrich's polynomial comparison [Ullrich 2013], the *Swendsen-Wang* dynamics on ferromagnetic Ising on $G$ inherits a polynomial mixing-time bound at every $\beta > 1$. This is the first such bound for SW; the exponent is far from tight but the *qualitative* "rapid mixing everywhere" is new.

## Method

The proof composes three classical comparison tricks. **(i)** Jerrum-Sinclair's [1993] canonical-paths construction for the (augmented) *even-subgraphs* / worm-process Markov chain bounds congestion of the worm dynamics by $n^4$. **(ii)** The Grimmett-Janson coupling (Theorem 3.5 of [Grimmett-Janson 2009]) shows that adding each edge $e \notin S$ independently with probability $p/(1-p)$ to a sample from the even-subgraph measure $\pi_{\text{even}}(p)$ yields a sample from $\pi_{RC}(2p, 2)$, giving a measure-level bridge between the two state spaces. **(iii)** A new comparison technique (Section 4, "Lifting Canonical Paths") *lifts* each worm canonical path to a flow on random-cluster configurations by mimicking the worm transition while edge-rerandomising the rest of the configuration; this is needed because standard Diaconis-Saloff-Coste / Dyer-Goldberg-Jerrum-Martin [DGJM 2006] comparison requires a single-state-space Markov-chain comparison and cannot handle the worm/RC state-space mismatch on its own. Lifting the flow loses only a polynomial factor (Lemma 5, $O(n^4)$), and Sinclair's canonical-paths bound (Theorem 1) converts this into the mixing-time statement.

## Quantum vs Classical

This is the **strongest "rapid-mixing-everywhere" classical result for any Ising sampler**, and the central reason the *ferromagnetic-Ising* row of the corpus grid is fully closed against quantum advantage:

- **What it closes**: ferromagnetic Ising on *any* graph at *any* $T > 0$ now has a polynomial-time *sampler* (not just an FPRAS for the partition function as in Jerrum-Sinclair 1993). High-$T$, near-critical (Bauerschmidt-Dagallier 2024), critical-2D (Lubetzky-Sly 2012), and low-$T$ are all subsumed under one dynamic algorithm.
- **Implication for quantum advantage**: the quantum-vs-classical regime grid (`drafts/literature-review/classical-gibbs-sampling-corpus.md`) lists the **ferromagnetic-Ising $\times$ all-$T$** cell as classically solved, so any KMS-DB Lindbladian on this Hamiltonian class is competing against a poly-time classical baseline. Quantum advantage — if any — must come from *non-ferromagnetic* (frustrated, antiferromagnetic in non-uniqueness, spin-glass, or non-stoquastic) Hamiltonians. This dovetails with the Sly / Sly-Sun hardness picture for antiferromagnetic two-spin models in the tree non-uniqueness region: classical hardness on the antiferromagnetic side, classical tractability on the ferromagnetic side.
- **Caveat on $q$**: the result is *specific to $q=2$* (Ising). For $q \ge 3$ the random-cluster Glauber chain is *not* polynomial — Gore-Jerrum [1999], Borgs-Chayes-et-al. [1999] and Blanca-Sinclair show first-order phase transitions create exponential bottlenecks for Potts-$q\ge3$ on the complete graph and the lattice. So the quantum-advantage door remains open for ferromagnetic Potts at low $T$ for $q \ge 3$; the present result does *not* close that cell.
- **Caveat on tightness**: the exponent (effectively $n^4 m^3 = O(n^{10})$ on dense graphs) is "well above the true answer" by the authors' own admission and "certainly too high to be of practical interest". This is a complexity-theoretic statement, not a practical sampler — but for the quantum-advantage question only the existence of a polynomial bound matters.
- **Comparison to quantum Gibbs samplers**: Chen-style KMS Lindbladians on the *classical* Ising Hamiltonian (which is diagonal, hence trivially commuting) reduce to classical Glauber-like dynamics with the same mixing barriers and no quantum speedup mechanism; nothing in this cell motivates the quantum machinery. The thesis Review chapter should cite Guo-Jerrum precisely as the closure result against which quantum methods would have to demonstrate value, and conclude that they cannot in this cell.

## Implications

- Closes the **ferromagnetic-Ising $\times$ any-$T$** cell of the corpus grid as polynomial-time classically samplable, on *arbitrary* graphs.
- Provides the first polynomial mixing-time upper bound for the Swendsen-Wang dynamics, partial answer to a question open since $\sim 1990$.
- Demonstrates a useful new comparison technique (path-lifting across distinct state spaces with a measure-level coupling) that may apply elsewhere whenever a coupling between two distributions is known but the dynamics live on different state spaces — potentially relevant to quantum-classical comparisons of Lindbladian and classical dynamics that share a common stationary measure up to a coupling.

## Open Questions

- **Tighter exponent**: the authors explicitly note the $n^4 m^3$ exponent is loose; matching the empirical SW dynamic exponent ($z \approx 0.25$ for 2D Ising) is wide open.
- **Random-cluster at $1 < q < 2$**: the monotone regime where rapid mixing is conjectured but not proven (the authors flag $1 < q < 2$ as "cause to be optimistic"); no quantum analogue exists.
- **Critical exponent**: critical 2D Ising mixing is $O(n^c)$ by Lubetzky-Sly 2012, but the $c$ is unknown; whether random-cluster dynamics is faster than spin Glauber at criticality is open.
- **Non-ferromagnetic extension**: the FK representation has signs for antiferromagnetic / spin-glass interactions, breaking the entire approach. This is precisely the regime where Sly-style hardness kicks in classically, and where quantum samplers might first beat classical — but the present technique gives no leverage there.

## Connections

- **Edwards-Sokal 1988** [edwards-sokal-1988-fk-representation]: provides the joint $(\sigma, \omega)$ representation that defines the random-cluster model and explains exactly when a cluster algorithm exists. The Guo-Jerrum analysis works at $q=2$ because the FK measure is non-negative; outside this regime (antiferromagnetic, complex weights, non-stoquastic quantum extensions) the joint measure has signs and the entire programme breaks. Cross-link this paper as the structural foundation.
- **Jerrum-Sinclair 1993** [classical-gibbs-sampling-corpus.md §5]: the FPRAS for the ferromagnetic Ising partition function via even-subgraph MCMC. Guo-Jerrum upgrade: from approximate counting to approximate *sampling* via the FK chain, keeping the polynomial regime universal.
- **Lubetzky-Sly 2012/2013** [lubetzky-sly-2012-critical-ising-polynomial, lubetzky-sly-2013-cutoff-ising-lattice]: complementary single-site Glauber results — polynomial at $T_c$ for 2D Ising (2012) and cutoff at all $\beta < \beta_c$ on $\mathbb{Z}^d$ (2013). Guo-Jerrum is graph-universal but loose on the exponent; LS gives sharp constants on the lattice but is geometry-specific. Together they cover ferromagnetic Ising completely.
- **Bauerschmidt-Dagallier 2024** [bauerschmidt-dagallier-2024-near-critical-ising-lsi]: log-Sobolev inequality for spin Glauber up to criticality with explicit susceptibility dependence. Strengthens the high-$T$/near-critical part of the cell with sharp LSI constants; Guo-Jerrum extends *down to* low $T$ via the FK detour, which spin Glauber does not (Glauber is exponentially slow on bipartite Ising at low $T$ due to symmetry breaking).
- **Ullrich 2013** [arXiv:1105.3665, ullrich-2013 in corpus]: polynomial comparison between Swendsen-Wang and single-bond / heat-bath spectral gaps on bounded-degree graphs. Guo-Jerrum's RC bound transfers to SW *via* this comparison. Without Ullrich's result the SW corollary would not follow.
- **Gore-Jerrum 1999** (and Borgs-Chayes-et-al. 1999): exponential SW slowdown for Potts $q \ge 3$ at the first-order transition on the complete graph and lattice. Marks the boundary of the present technique: $q=2$ is rapidly mixing, $q \ge 3$ is torpid in general — the FK programme is a $q=2$ phenomenon.
- **Bravyi-Gosset 2017** [bravyi-gosset-2017-quantum-ferromagnets]: classical FPRAS for the partition function of *quantum* ferromagnetic XY/Heisenberg / transverse-field-Ising on arbitrary graphs at any finite $T$, using Jerrum-Sinclair-Vigoda permanents. Together with Guo-Jerrum this closes both the *classical-ferromagnetic-Ising* and *quantum-stoquastic-ferromagnetic-XY* cells against quantum-Gibbs-sampler advantage.
