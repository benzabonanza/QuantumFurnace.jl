---
paper: "Sly (2010)"
title: "Computational Transition at the Uniqueness Threshold"
arxiv: "1005.5584"
year: 2010
venue: "FOCS 2010"
pdf: "n/a"

temperature: [low-T]
commutativity: [commuting]
locality: [sparse, bounded-degree]
particle-statistics: [n/a]
hamiltonian-models: [hardcore, antiferromagnetic-2-spin]

paradigm: [classical-MCMC, hardness]
quantum-or-classical: [classical, lower-bound]

result-type: [hardness, NP-hard, complexity-lower]
key-scaling: "no FPRAS for $Z_G(\\lambda)$ on bounded-degree graphs at $\\lambda \\in (\\lambda_c(\\Delta), \\lambda_c(\\Delta)+\\epsilon(\\Delta))$ unless $\\mathrm{NP}=\\mathrm{RP}$; $\\lambda_c(\\Delta) = (\\Delta-1)^{\\Delta-1}/(\\Delta-2)^\\Delta$"

related: ["anari-liu-oveisgharan-2020-spectral-independence-hardcore", "weitz-2006-counting-independent-sets-tree-threshold", "chen-liu-vigoda-2021-optimal-mixing-glauber-hde", "sly-sun-2014-counting-two-spin-d-regular", "dyer-frieze-jerrum-2002-counting-independent-sets-sparse"]
---

# Computational Transition at the Uniqueness Threshold — Sly (2010)

**One-sentence takeaway**: For the hardcore (independent-set) model on bounded-degree graphs, no FPRAS exists for the partition function at fugacities just above the $\Delta$-regular tree uniqueness threshold $\lambda_c(\Delta) = (\Delta-1)^{\Delta-1}/(\Delta-2)^\Delta$ unless $\mathrm{NP}=\mathrm{RP}$ — pinning the *computational* phase transition exactly to the *statistical-physics* phase transition and closing the lower-bound side of the regime that Weitz / Anari–Liu–Oveis Gharan close on the upper side.

## Setting

- **Target**: hardcore model $\mu_\lambda(I) \propto \lambda^{|I|}$ on independent sets $I$ of a graph $G=(V,E)$ with $|V|=n$ and maximum degree $\le \Delta$. Partition function $Z_G(\lambda) = \sum_I \lambda^{|I|}$.
- **Hamiltonian class**: classical, antiferromagnetic 2-spin / commuting. Equivalent to a hard-constraint Ising model in the limit of infinite repulsion between adjacent occupied vertices.
- **Threshold**: tree non-uniqueness regime, i.e. fugacities $\lambda > \lambda_c(\Delta)$ on the $\Delta$-regular infinite tree, where the Gibbs measure has multiple infinite-volume extensions. *Low-T side* of the transition in the antiferromagnetic statistical-physics convention (high fugacity = strong constraint = low effective temperature).
- **Access model**: explicit graph $G$, decision/approximation oracle for $Z_G(\lambda)$. The hardness target is FPRAS / approximate sampling, not exact counting (which is already $\#$P-hard).

## Main Results

- **Theorem 1 (main)**: For every $\Delta \ge 3$ there exists $\epsilon(\Delta) > 0$ such that, unless $\mathrm{NP}=\mathrm{RP}$, no FPRAS exists for $Z_G(\lambda)$ on graphs of maximum degree $\le \Delta$ when $\lambda \in (\lambda_c(\Delta),\, \lambda_c(\Delta) + \epsilon(\Delta))$. Equivalently, no polynomial-time approximate-sampling algorithm exists for $\mu_\lambda$ in the same window.
- **Concrete instance**: at $\Delta=6$, $\lambda=1$ already lies past $\lambda_c(6)$, so approximately counting independent sets on graphs of max degree $6$ is NP-hard under randomized reductions.
- **Threshold formula**: $\lambda_c(\Delta) = (\Delta-1)^{\Delta-1}/(\Delta-2)^\Delta \approx e/(\Delta-2)$, the classical Kelly threshold for uniqueness of the hardcore Gibbs measure on the $\Delta$-regular tree.
- **Status of conjecture**: closes the matching half of the dichotomy conjectured by Mossel and others — Weitz's FPTAS works iff $\lambda < \lambda_c(\Delta)$, Sly's hardness kicks in at $\lambda > \lambda_c(\Delta)$. (Sharpened later to the entire non-uniqueness region in Sly–Sun (2014) and to general antiferromagnetic 2-spin in Galanis–Štefankovič–Vigoda.)

## Method

A reduction from MAX-CUT through a *random bipartite gadget*. Sly takes a random $\Delta$-regular bipartite graph and shows, via a *second-moment argument applied to a planted-phase model* and a *reconstruction-on-the-tree* analysis, that in the non-uniqueness regime the hardcore Gibbs measure concentrates on two distinguishable "phases" (each side mostly occupied / mostly unoccupied). These two phases act as a binary signal that can be locally decoded, so an FPRAS for $Z_G(\lambda)$ would let one decide MAX-CUT instances embedded as gadgets — turning approximate counting into approximate maximum cut. Formally this rigorizes the physics "replica-symmetric / replica-symmetry-breaking" heuristic for the hardcore antiferromagnet.

## Quantum vs Classical

- **Baseline**: Weitz (2006) gives a deterministic quasi-polynomial FPTAS at $\lambda < \lambda_c(\Delta)$ via the self-avoiding-walk tree; Anari–Liu–Oveis Gharan (2020) and Chen–Liu–Vigoda (2021) sharpen this to a polynomial FPRAS / $O(n\log n)$ Glauber chain. So *below* $\lambda_c(\Delta)$ the classical landscape is fully resolved; Sly's theorem says that crossing $\lambda_c(\Delta)$ from below *immediately* triggers worst-case hardness of approximation.
- **Gap**: this is a *lower bound on classical complexity* — there is no classical "best algorithm" to compare against above $\lambda_c$, since none exists in poly-time. The structural gap is *infinite* (in the sense of a complexity-class barrier) modulo $\mathrm{NP} = \mathrm{RP}$.
- **Source of the difference**: the obstruction is not a sign problem (the model is classical and commuting) — it is a *statistical-physics phase transition* that creates an exponential bottleneck in the configuration space, which Sly converts into a worst-case hardness reduction. The Gibbs measure splits into long-range-correlated phases that any local sampler must distinguish, but doing so encodes a hard combinatorial decision.
- **Caveat for quantum**: the hardness is for the *partition function approximation* problem and the corresponding *approximate sampler*. A quantum Gibbs sampler that prepared $\rho_\beta = e^{-\beta H_{\text{hardcore}}}/Z$ to total-variation precision in poly-time would (via the standard partition-function-from-samples reduction, e.g. Štefankovič–Vempala–Vigoda annealing) yield an FPRAS for $Z_G(\lambda)$ — which Sly forbids unless $\mathrm{NP}=\mathrm{RP}$. So **a quantum Gibbs sampler cannot in general beat classical here either** in the non-uniqueness window: hardness is at the level of the *output*, not the *implementation*. This is the one place in the corpus where a quantum advantage story is *blocked structurally and not by a classical algorithmic existence claim*.

## Implications for Quantum Advantage

- **Regime cell**: low-T (antiferromagnetic non-uniqueness, $\lambda > \lambda_c$) $\times$ commuting (classical Hamiltonian) $\times$ sparse / bounded-degree $\times$ classical spin (n/a particle stats).
- **What this changes**: *closes the lower-bound side* of the classical / commuting / low-T cell for antiferromagnetic 2-spin. Together with Anari–Liu–Oveis Gharan / Weitz on the $\lambda < \lambda_c$ side, this is the canonical "computational threshold = phase-transition threshold" benchmark in the corpus.
- **Promising or not**: **explicitly *not* promising for quantum advantage**, and stronger than the analogous statement for the upper side: any putative quantum Gibbs sampler in this regime — KMS Lindbladian, dissipative-preparation, block-encoded Metropolis — would violate $\mathrm{NP} \ne \mathrm{RP}$ via the partition-function reduction.
- **Where the door is left open for quantum**: the hardness is *worst-case over bounded-degree graphs*, not average-case. A quantum sampler that targets *physically motivated* (e.g. lattice, random-regular, or planted) instances and bypasses the gadget-graphs Sly constructs is not ruled out — the reduction's gadgets are adversarial bipartite expanders, not lattices. Quantum advantage stories in this cell must therefore navigate the partition-function reduction and confine themselves to instance distributions where the Sly gadget cannot be embedded. [CHECK] whether subsequent work (Galanis–Štefankovič–Vigoda, Sly–Sun) closes this loophole on planar / lattice graphs as well.

## Open Questions / Limitations

- Original paper covers only a small window $\lambda_c(\Delta) < \lambda < \lambda_c(\Delta) + \epsilon(\Delta)$; the *entire* non-uniqueness region is closed by Sly–Sun (2014) and Galanis–Štefankovič–Vigoda (2015).
- Hardness is *worst-case* over bounded-degree graphs; it does not rule out polynomial algorithms on structured subclasses (planar graphs, expanders with bounded girth, lattices). Lattice hardness is still partly open in this fugacity window. [CHECK]
- The reduction goes through a randomized argument over random bipartite gadgets, hence the conclusion is conditional on $\mathrm{NP} \ne \mathrm{RP}$ rather than $\mathrm{P} \ne \mathrm{NP}$. Tightening to deterministic hardness would require derandomizing the second-moment / planted-phase analysis.
- No statement about *finite* approximation factors (e.g. $2$-approximation): the FPRAS is ruled out, but constant-factor bounds are not addressed here.
- No quantum lower bound is *proven* — the implication "quantum FPRAS would imply $\mathrm{NP}=\mathrm{RP}$" is via Štefankovič–Vempala–Vigoda annealing, which is classical and assumes the quantum sampler's output can be classically post-processed in poly-time. This is the standard reduction but worth flagging.

## Connections

- **Anari–Liu–Oveis Gharan (2020)** [`anari-liu-oveisgharan-2020-spectral-independence-hardcore`]: matching upper-bound side. Together they cement "computational threshold = uniqueness threshold" for the hardcore model on bounded-degree graphs — Anari et al. give an FPRAS up to $\lambda_c(\Delta)$, Sly forbids one beyond.
- **Weitz (2006)** [`weitz-2006-counting-independent-sets-tree-threshold`]: gives the deterministic FPTAS in the tree-uniqueness regime; combined with Sly this is the original sharp dichotomy paper pair.
- **Chen–Liu–Vigoda (2021)** [`chen-liu-vigoda-2021-optimal-mixing-glauber-hde`]: sharpens Anari et al. to optimal $O(n \log n)$ Glauber mixing in the same regime; Sly's hardness shows this is the best-possible regime.
- **Sly, Sun (2014)** [`sly-sun-2014-counting-two-spin-d-regular`]: extends Sly (2010) to the *entire* tree non-uniqueness region for general antiferromagnetic 2-spin (including Ising with external field), closing the dichotomy.
- **Dyer, Frieze, Jerrum (2002)** [`dyer-frieze-jerrum-2002-counting-independent-sets-sparse`]: classical predecessor — torpid Glauber mixing at $\Delta \ge 6$ and FPRAS hardness at $\Delta \ge 25$. Sly (2010) sharpens the threshold from "$\Delta \ge 25$" to "exactly $\lambda_c(\Delta)$".
- **Bravyi–DiVincenzo–Oliveira–Terhal (2008)** [`bravyi-divincenzo-oliveira-terhal-2008-stoquastic-lh`]: complementary "easy quantum" anchor — together with Sly, they triangulate the corpus: stoquastic = classical-easy (BDOT), antiferromagnetic non-uniqueness = classical-hard *and* quantum-hard (Sly), generic non-stoquastic low-T = the only cell where quantum can plausibly win.
- *Thesis Review chapter:* this is the canonical low-T classical-hardness anchor in the 3×3 grid (`classical-gibbs-sampling-corpus.md` §5). It is the *sharpest available no-go* for any Gibbs sampler — classical or quantum — in the antiferromagnetic 2-spin / commuting / low-T cell, and constrains what a positive quantum-advantage claim in this row can look like (specifically: quantum advantage must target instance distributions where Sly's gadget is not embeddable).
