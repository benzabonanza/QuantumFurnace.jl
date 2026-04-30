---
paper: "Weitz (2006)"
title: "Counting Independent Sets up to the Tree Threshold"
arxiv: "n/a"
year: 2006
venue: "STOC 2006"
pdf: "n/a"

temperature: [high-T, intermediate-T]
commutativity: [commuting]
locality: [sparse, bounded-degree]
particle-statistics: [n/a]
hamiltonian-models: [hardcore, antiferromagnetic-2-spin]

paradigm: [classical-deterministic, correlation-decay]
quantum-or-classical: [classical]

result-type: [FPTAS, deterministic-algorithm]
key-scaling: "FPTAS runtime (n/eps)^{O(log Delta)} on graphs of max degree Delta, valid for any lambda < lambda_c(Delta) = (Delta-1)^{Delta-1}/(Delta-2)^{Delta}; quasi-polynomial in n once Delta is unbounded, polynomial for Delta = O(1)."

related: ["anari-liu-oveisgharan-2020-spectral-independence-hardcore", "chen-liu-vigoda-2021-optimal-mixing-glauber-hde", "sly-2010-computational-transition-uniqueness"]
---

# Counting Independent Sets up to the Tree Threshold — Weitz 2006

**One-sentence takeaway**: Gives the first deterministic FPTAS for the hardcore partition function $Z_G(\lambda)$ on any graph of max degree $\le \Delta$ throughout the entire tree-uniqueness regime $\lambda < \lambda_c(\Delta) = (\Delta-1)^{\Delta-1}/(\Delta-2)^\Delta$, by reducing computation of marginals on $G$ to computation on a *self-avoiding-walk tree* on which strong spatial mixing can be proven directly — establishing $\lambda_c(\Delta)$ as the canonical algorithmic threshold and introducing the SAW-tree reduction reused by essentially all subsequent work.

## Setting

- **Target**: hardcore (independent-set) Gibbs measure on a graph $G=(V,E)$, $|V|=n$, $\mu(I)\propto\lambda^{|I|}$ over independent sets $I\subseteq V$. Partition function $Z_G(\lambda)=\sum_I \lambda^{|I|}$. Equivalently the antiferromagnetic 2-spin / lattice-gas partition function at fugacity $\lambda$.
- **Hamiltonian class**: classical, commuting (diagonal in the configuration basis); the canonical antiferromagnetic 2-spin model.
- **Temperature regime**: tree-uniqueness regime — $\lambda<\lambda_c(\Delta)$ where $\lambda_c(\Delta)=(\Delta-1)^{\Delta-1}/(\Delta-2)^\Delta\approx e/(\Delta-2)$ is the critical activity for uniqueness of the infinite-volume Gibbs measure on the infinite $\Delta$-regular tree $\mathbb{T}_\Delta$. This is "high-T" in the antiferromagnetic-spin-system sense (small $\lambda$ = effective inverse temperature far from any phase transition).
- **Access**: explicit graph $G$ with $\Delta(G)\le\Delta$; the algorithm is deterministic, not MCMC.

## Main Results

- **Theorem (FPTAS, informal)**: for every $\Delta\ge 3$ and every $\lambda<\lambda_c(\Delta)$, there is a deterministic algorithm that, given a graph $G$ of max degree $\le\Delta$ and accuracy $\varepsilon>0$, outputs $\widehat Z$ with $|\widehat Z - Z_G(\lambda)|\le\varepsilon Z_G(\lambda)$ in time $(n/\varepsilon)^{O(\log\Delta)}$.
- **Strong spatial mixing on $\mathbb{T}_\Delta$ up to $\lambda_c(\Delta)$**: for any $\lambda<\lambda_c(\Delta)$, the marginal of the root in the hardcore measure on a depth-$L$ subtree of $\mathbb{T}_\Delta$ depends on the boundary configuration only through an error decaying as $\exp(-\Omega(L))$, *uniformly* in the boundary. This is the technical core; the threshold $\lambda_c(\Delta)$ is sharp (unique fixed point of the tree recursion stops contracting at $\lambda_c$).
- **SAW-tree reduction (Lemma)**: for any vertex $v$ in any graph $G$, the hardcore marginal $\Pr_\mu[v\in I]$ on $G$ equals the corresponding root marginal on the *self-avoiding-walk tree* $T_{\mathrm{SAW}}(G,v)$ — a tree of self-avoiding walks from $v$ in $G$ — with a specific "freezing" rule for vertices that close a cycle. The SAW tree has max degree $\le\Delta$ and depth at most $n$.
- **Algorithm**: compute each marginal $\Pr[v\in I\mid \text{partial config}]$ by truncating $T_{\mathrm{SAW}}(G,v)$ at depth $L=O(\log(n/\varepsilon)/\log(1/\rho))$ where $\rho<1$ is the contraction rate of the tree recursion at $\lambda$; multiply marginals by self-reducibility to recover $Z_G(\lambda)$. Runtime is dominated by enumerating the truncated SAW tree, of size $\Delta^L = (n/\varepsilon)^{O(\log\Delta)}$.

## Method

(1) **SAW-tree reduction**: the hardcore marginal at any vertex of any graph equals the root marginal on a tree of self-avoiding walks from that vertex (with closed-cycle vertices "frozen" alternately to occupied/unoccupied) — this lifts the *graph* problem to a *tree* problem exactly. (2) **Tree recursion + correlation decay**: on a tree, the root's marginal can be computed by a one-dimensional recursion in the children's marginals; one shows the recursion is contractive in a suitable potential (a logarithmic / arctangent reparameterization) iff $\lambda<\lambda_c(\Delta)$, giving exponentially fast strong spatial mixing. (3) **Truncation**: cut the SAW tree at depth $L=O(\log(n/\varepsilon))$ — by SSM the truncation error is $\le\varepsilon$, and the cost is the size of the truncated tree, $\Delta^L=(n/\varepsilon)^{O(\log\Delta)}$.

## Quantum vs Classical

- **Baseline**: prior MCMC-based FPRAS (Luby–Vigoda; Vigoda; Dyer–Greenhill) only reached $\lambda<2/(\Delta-2)$ — far below $\lambda_c(\Delta)\approx e/(\Delta-2)$ — and gave randomized algorithms with no deterministic counterpart. Beyond $\lambda_c$ no FPRAS can exist (Sly 2010 / Sly–Sun 2014, RP$\ne$NP).
- **Gap**: Weitz pushes the algorithmic threshold from $2/(\Delta-2)$ to $\lambda_c(\Delta)$ — sharp, since hardness sets in immediately above. The improvement is qualitative (covers the full uniqueness phase), not in runtime exponent.
- **Source of the difference**: tree spatial mixing (an intrinsic statistical-mechanical property of the infinite tree) is *exactly* the right algorithmic resource — the SAW-tree reduction makes any graph "look like" a tree for marginal computation, so the tree threshold transfers verbatim to general bounded-degree graphs. No MCMC machinery is used; correlation decay replaces canonical-paths / conductance.
- **Caveat — quasi-polynomial in $\Delta$**: the runtime $(n/\varepsilon)^{O(\log\Delta)}$ is polynomial in $n$ only for *bounded* $\Delta$. For graphs with $\Delta=\Theta(\log n)$ or larger, Weitz is super-polynomial in $n$ — the gap that Anari–Liu–Oveis Gharan (2020) close by giving an FPRAS with no $\Delta$ dependence in the runtime exponent and Chen–Liu–Vigoda (2021) sharpen further to optimal $O(n\log n)$ Glauber mixing.
- **Quantum side**: paper makes no quantum claim. For the thesis Review chapter this is a *purely classical, deterministic* benchmark in the **commuting × high-T (uniqueness)** cell. There is no plausible quantum sampler that beats this combination of guarantees in the same regime: the Hamiltonian is diagonal, so any quantum Gibbs sampler degenerates to a quantum implementation of a classical chain plus block-encoding overhead $\Omega(\log(1/\varepsilon))$ — strictly worse, asymptotically, than $\widetilde O(n^2)$ Glauber under spectral-independence (CLV21). Weitz also has the unique virtue of being *deterministic*, which no quantum sampler in the corpus matches.

## Implications for Quantum Advantage

- **Regime cell**: high-T (tree-uniqueness regime) $\times$ commuting (classical Hamiltonian) $\times$ sparse / bounded-degree $\times$ classical spin (n/a particle stats) — top-left of the thesis grid.
- **What this changes**: opens the upper-bound side of the classical / commuting / high-T cell at the *sharp* statistical-physics threshold $\lambda_c(\Delta)$, defining the algorithmic-threshold-equals-uniqueness-threshold paradigm. Together with Sly hardness above $\lambda_c$, this set the agenda that ALOG20 and CLV21 later sharpened.
- **Promising or not (for quantum advantage)**: explicitly *not promising*. The cell admits a deterministic FPTAS (Weitz) and an optimal-time FPRAS (CLV21) matching a hardness lower bound (Sly), with sharp threshold matching the statistical-mechanical phase transition. Any quantum advantage in the thesis chapter must come from cells where this entire paradigm — tree spatial mixing, SAW-tree reduction, classical correlation decay — has no analogue: non-stoquastic / low-T quantum Hamiltonians, where the configuration space does not decompose along a tree and noncommutativity destroys the recursion.

## Open Questions / Limitations

- **Quasi-polynomial in $\Delta$**: $(n/\varepsilon)^{O(\log\Delta)}$ is polynomial only for bounded $\Delta$; closed by ALOG20 (FPRAS with no $\Delta$-dep. in the exponent) and CLV21 (optimal $O(n\log n)$).
- **Deterministic but slow in $\varepsilon$**: $(1/\varepsilon)^{O(\log\Delta)}$ vs the FPRAS $\log(1/\varepsilon)$ — the deterministic-vs-randomized tradeoff for marginal computation; no known way to derandomize the SAW-tree approach to genuine polynomial-in-$\log(1/\varepsilon)$ runtime.
- **Antiferromagnetic 2-spin only in this paper**: extended to general antiferromagnetic 2-spin systems (Li–Lu–Yin) and conjecturally beyond; multi-spin (e.g. $q$-colourings) requires different machinery (Gamarnik–Katz; Liu–Lu).
- **Tree-uniqueness threshold**: the paper covers the entire conjectured tractable region for hardcore — there is no algorithmic slack to exploit on the classical side beyond constants and the runtime exponent.
- **Pure classical**: no statement about quantum spin systems; the SAW-tree reduction crucially uses the classical product structure of configurations and has no obvious operator-valued analogue.

## Connections

- **Sharpened by** Anari–Liu–Oveis Gharan (2020) [`anari-liu-oveisgharan-2020-spectral-independence-hardcore`]: same regime $\lambda<\lambda_c(\Delta)$, but FPRAS via spectral independence on Glauber dynamics, runtime $\widetilde O(n^{2+C(\delta)})$ with no $\Delta$-dep. in the exponent — fixes Weitz's quasi-polynomial-in-$\Delta$ runtime. ALOG20 still uses Weitz's SAW-tree internally to bound the influence sum.
- **Sharpened further by** Chen–Liu–Vigoda (2021) [`chen-liu-vigoda-2021-optimal-mixing-glauber-hde`]: optimal $O(n\log n)$ Glauber mixing in the same regime via MLSI / entropy factorization; replaces Weitz's $(n/\varepsilon)^{O(\log\Delta)}$ FPTAS by an $\widetilde O(n^2)$ FPRAS for $Z_G(\lambda)$.
- **Matched by** Sly (2010) [`sly-2010-computational-transition-uniqueness`]: NP-hardness of approximation for $\lambda>\lambda_c(\Delta)$. Weitz + Sly establish "algorithmic threshold = uniqueness threshold" — the canonical phase-transition-matches-complexity result, mirrored on the multi-spin side by Sly–Sun 2014 and Galanis–Štefankovič–Vigoda 2015.
- **Methodological progenitor**: the SAW-tree + tree-recursion + potential-function correlation-decay recipe is reused throughout the deterministic-counting literature (Bayati–Gamarnik–Katz–Nair–Tetali on monomer-dimer; Li–Lu–Yin for general antiferromagnetic 2-spin; Liu–Lu–Zhang) and inside MCMC analyses (ALOG20's influence-sum bound).
- Within the thesis grid: this is the **classical Ising/hardcore × high-T** Tier-1 anchor in `classical-gibbs-sampling-corpus.md` Section 5, paired with Sly hardness. Together with ALOG20 and CLV21 it pins the entire upper- and lower-bound side of this cell, leaving no room for quantum advantage in the regime — consistent with the chapter's punchline that quantum-vs-classical contests must be staged in non-stoquastic / low-T cells where this whole apparatus has no analogue.
