---
paper: "Woodard, Schmidler, Huber (2009)"
title: "Sufficient conditions for torpid mixing of parallel and simulated tempering"
arxiv: "n/a"
year: 2009
venue: "Electron. J. Probab. 14:780–804, paper 38"
pdf: "n/a"

temperature: [low-T]
commutativity: [commuting]
locality: [mean-field, geometric-local]
particle-statistics: [spin]
hamiltonian-models: [ising, potts, gaussian-mixture]

paradigm: [classical-MCMC, parallel-tempering]
quantum-or-classical: [classical]

result-type: [mixing-time-lower, gap-bound, no-go]
key-scaling: "Spec(PT/ST) <= c_1 e^{-c_2 N} on any target whose hottest distribution has a 'persistent' narrow mode of asymptotically vanishing relative width but non-vanishing mass; consequently t_mix >= e^{Omega(N)} for parallel and simulated tempering with any temperature schedule"

related: ["bovier-den-hollander-2015-metastability", "cuff-et-al-2012-mean-field-potts-glauber", "madras-zheng-2003-swapping-algorithm"]
---

# Sufficient conditions for torpid mixing of parallel and simulated tempering — Woodard, Schmidler, Huber (2009)

**One-sentence takeaway**: Constructs a *generic* persistence criterion — a narrow mode whose relative width shrinks while its mass stays $\Theta(1)$ — under which **every** parallel-tempering (PT) and simulated-tempering (ST) chain, on **any** finite temperature ladder, has spectral gap $\le c_1\, e^{-c_2 N}$, i.e. tempering is *unconditionally* torpid on a wide class of low-T multimodal targets, including ferromagnetic mean-field Potts ($q\ge 3$) and mixtures of Gaussians with unequal covariances.

## Setting

- **Targets.** Sequences of probability measures $\pi_N$ on a state space $\mathcal{X}_N$ (finite or $\mathbb{R}^M$) indexed by a system-size or sharpness parameter $N$. The flagship cases are $\pi_N(\sigma)\propto e^{-\beta H_N(\sigma)}$ for mean-field Potts on $N$ sites with $q\ge 3$, and a fixed mixture of normals on $\mathbb{R}^M$ with the narrowest component's covariance scaling as $\sigma_N \to 0$.
- **Algorithms.** *Parallel tempering* (PT) on a product space $\mathcal{X}^L$ at $L$ inverse temperatures $\beta_0 > \beta_1 > \dots > \beta_{L-1} \ge 0$, with stationary measure $\bigotimes_\ell \pi^{(\beta_\ell)}$ and Metropolis swap moves between adjacent rungs; *simulated tempering* (ST) on $\mathcal{X}\times\{0,\dots,L-1\}$ with adaptive temperature index. In both cases the within-rung kernel is any $\pi^{(\beta_\ell)}$-reversible chain (Glauber, random-walk Metropolis, …); the swap/index move is the only mechanism by which information is transported across the temperature ladder.
- **Question.** Lower bound on the mixing time of PT/ST that holds *uniformly* over the choice of temperature ladder and within-rung kernel — i.e. a structural no-go for the entire tempering family, not just a specific schedule.
- **Comparison reference.** Madras–Zheng (2003) had previously shown a *positive* result: PT mixes in polynomial time on **mean-field Ising** (no field, $q=2$). The present paper exhibits a sharp dividing line: a single structural property of the target — persistence — flips PT/ST from polynomial to exponentially slow, and Potts at $q\ge 3$ falls on the wrong side of that line.

## Main Results

The notation in the abstract and §2 of the paper:

- **Persistence (the structural hypothesis).** A sequence of multimodal distributions $\pi_N$ is *persistent* if it admits a partition $\mathcal{X}_N = A_N \sqcup B_N$ such that, **at the hottest temperature $\beta_{L-1}$ of the ladder**, (i) $\pi^{(\beta_{L-1})}_N(A_N) \ge \delta > 0$ uniformly in $N$ — the "narrow" mode retains $\Theta(1)$ mass even at the highest temperature — and (ii) the bottleneck conductance $\Phi^{(\beta_{L-1})}_N(A_N)\to 0$ exponentially in $N$ for the within-rung kernel at that temperature. Informally: a tall-narrow peak that is too localised for the within-rung sampler to leave even after maximal heating.
- **Theorem (torpid mixing, abstract restated).** Under persistence, for any choice of $L\ge 1$ ladder temperatures $\beta_0>\dots>\beta_{L-1}$ and any within-rung kernels, the second-largest eigenvalue $\lambda_2$ of the PT (resp. ST) chain satisfies
  $$ 1 - \lambda_2 \;\le\; c_1\, e^{-c_2 N}, $$
  hence $t_{\mathrm{mix}}(\varepsilon) \ge \tfrac{1}{2}(1-\lambda_2)^{-1}\log(1/(2\varepsilon)) \ge \tfrac{1}{2}c_1^{-1}e^{c_2 N}\log(1/(2\varepsilon))$ — torpid for both PT and ST.
- **Corollary 1 (mean-field Potts, $q\ge 3$).** For ferromagnetic mean-field Potts with $q\ge 3$ at any $\beta$ in the metastability interval (above the spinodal $\beta_s(q)$), persistence holds with $A_N$ = neighbourhood of an ordered configuration; PT and ST are torpid with rate $c_2 = c_2(\beta,q) > 0$. Recovers and unifies Bhatnagar–Randall (2004).
- **Corollary 2 (mixture of normals with unequal covariances).** Fix a mixture $\pi = \sum_k w_k\,\mathcal{N}(\mu_k,\Sigma_k^{(N)})$ on $\mathbb{R}^M$ with at least one $\Sigma_k^{(N)} = \sigma_N^2 I$ for $\sigma_N\to 0$ (sharp narrow mode) while another component has $\Theta(1)$ covariance and $w_k$ all bounded away from $0$. Then persistence holds for any reasonable Metropolis within-rung kernel and PT/ST is torpid in $N$ (the sharpness parameter, *not* dimension). This is the *new* example of the paper.
- **Models that do *not* satisfy persistence.** Mean-field Ising with no field ($q=2$): symmetric two-mode target, both modes have width $\Theta(1)$ on the spin proportions axis, no persistent narrow component — Madras–Zheng's polynomial PT result is consistent. More generally, multimodal targets whose modes have *comparable widths* under the family of tempered measures escape the hypothesis.

## Method

Conductance / canonical-paths upper bound on the spectral gap, lifted to the product (PT) or extended (ST) state space. The key step is to exhibit a cut $S \subset \mathcal{X}^L$ (resp. $\mathcal{X}\times\{0,\dots,L-1\}$) for which both (a) $\pi^{\otimes}(S)$ is bounded away from $0$ and $1$, and (b) the *flow* across $S$ — bounded by the maximum of within-rung conductances and adjacent-rung swap acceptance — is exponentially small. Persistence makes the within-rung conductance at *every* rung exponentially small (the hottest rung still feels the narrow mode), so no swap schedule can repair the bottleneck. The argument is uniform over $L$ and over the temperature ladder, which is what gives it teeth as a no-go.

## Quantum vs Classical

- **Baseline.** Pre-2009, the standard classical workaround for low-T metastability — when single-site Glauber is torpid by Bovier–den Hollander Eyring–Kramers $e^{\beta\Gamma^*}$ — was to run PT/ST: heat the system above the transition, swap copies across rungs, and exploit polynomial mixing of the hottest rung. Madras–Zheng (2003) had endorsed this for mean-field Ising. This paper closes that escape route on a structurally identified class: whenever the high-temperature target itself contains a tall-narrow peak with $\Theta(1)$ mass, PT/ST is **provably exponential**, regardless of how clever the temperature ladder is.
- **Gap.** *Super-polynomial* (exponential in $N$) lower bound, in the system-size / sharpness parameter $N$. The bound is *unconditional* in the algorithmic axis (any ladder, any within-rung kernel) — only the target's structural property matters.
- **Source of the difference.** **Persistence is a property of the family of tempered measures, not of any one temperature.** Tempering exploits the assumption that heating *flattens* the landscape; if the narrow mode persists at all reachable temperatures, no rung sees a uni-modal target and the swap mechanism cannot help. This is structurally different from the Glauber-only torpidity of Cuff et al. (`cuff-et-al-2012-mean-field-potts-glauber`) and Bovier–den Hollander (`bovier-den-hollander-2015-metastability`): those say the *cold* dynamics is slow, but leave open whether tempering helps. Woodard–Schmidler–Huber says tempering does *not* help on mean-field Potts (recovering Bhatnagar–Randall) and identifies the *general* obstruction.
- **Caveat.** The hypothesis is structural and not always easy to verify: one must exhibit a $\Theta(1)$-mass narrow mode at the *hottest* rung. For models where heating eventually flattens all modes (e.g. mean-field Ising at any field, mixtures of normals with comparable covariances), persistence fails and PT can in principle be polynomial — see the companion paper Woodard–Schmidler–Huber (2009b, *Ann. Appl. Probab.* 19:617–640, arXiv:0906.2341) for matching *rapid*-mixing conditions on the same family.

## Implications for Quantum Advantage

- **Regime cell (thesis 3×3 grid).** **Classical Ising/Potts × low-T** (and mean-field × low-T more generally; also continuous-space mixture targets which the thesis can cite as a side example). The bottleneck is *structural metastability*, not a sign problem — Hamiltonian is classical-commuting, so this is firmly in the no-quantum-axis-issue zone, and the comparison is apples-to-apples.
- **What this changes.** *Closes* a major classical-side escape route for the low-T cell. Together with Bovier–den Hollander (`bovier-den-hollander-2015-metastability`, geometric Glauber lower bound), Cuff et al. (`cuff-et-al-2012-mean-field-potts-glauber`, exponential-in-$n$ Glauber lower bound on mean-field Potts), Bhatnagar–Randall (PT-specific torpidity on mean-field Potts), Gore–Jerrum (Swendsen–Wang torpidity on mean-field Potts), and Galanis–Štefankovič–Vigoda ($\#\mathrm{BIS}$-hardness of antiferromagnetic Potts), Woodard–Schmidler–Huber establishes that the **classical Ising/Potts × low-T** cell of the corpus is a *robust* no-go: no known classical paradigm — local moves, cluster moves, or temperature ladders — handles a persistent narrow mode in polynomial time.
- **Promising or not.** **Strongly promising for quantum.** This is the cleanest argument that the low-T multimodal cell *needs a genuinely non-classical mechanism* to resolve: tempering cannot tunnel through a persistent narrow peak because all of its rungs see the same peak. A coherent quantum sampler — KMS Lindbladian, dissipative engineering, or quantum walk — has an axis (coherent superposition over the configuration space) that classical tempering structurally lacks. If a KMS Lindbladian on mean-field Potts can be shown to escape persistence-style bottlenecks (e.g. via a noncommutative log-Sobolev or detectability-lemma argument), it would constitute a candidate **quantum-first** result on a sign-problem-free Hamiltonian — the cleanest possible benchmark.
- **Counterpoint to keep honest.** The persistence hypothesis is a *width* condition on the target, not a *barrier height* condition. A KMS Gibbs sampler still sees the same target and the same Bohr-frequency structure; nothing in this paper *forces* a quantum slowdown, but nothing rules out a quantum slowdown by an analogous mechanism either. The thesis numerics chapter is the right venue to test whether KMS Lindbladians inherit a persistence-like obstruction on mean-field Potts. [CHECK] whether any quantum analogue of "persistence" has been formulated in the literature.

## Open Questions / Limitations

- **One-sided.** Sufficient conditions for *torpidity*. Matching sufficient conditions for *rapid* PT mixing are in the companion paper (Woodard–Schmidler–Huber 2009b, arXiv:0906.2341) — they require persistence to fail and additional regularity (bounded width-ratio across modes). A *necessary and sufficient* characterisation is open.
- **Mean-field assumption in the Potts example.** The Potts torpidity inherited from this paper is on the complete graph; lattice $\mathbb{Z}^d$ Potts at large $q$ also has a first-order transition but the relevant nucleation picture (surface-tension barrier, droplet growth) is geometrically different and the persistence hypothesis must be re-checked rung-by-rung.
- **Mixture-of-normals example assumes sharpness $\sigma_N\to 0$.** Realistic Bayesian-posterior multimodality (where component widths are fixed but possibly unequal in dimension $M$) is partially covered but the dependence on $M$ is not optimised.
- **Hopper-style algorithms not covered.** Combined PT + cluster moves (e.g. Houdayer cluster + tempering) and continuous-time replica-exchange variants escape the strict generator structure assumed here; whether persistence still implies torpidity for these hybrids is open.
- **No quantum analogue.** No statement about KMS Lindbladians, quantum walks, or any quantum sampler on the same persistent targets. The cleanest open thesis question: *does a KMS detailed-balance Lindbladian on mean-field Potts ($q\ge 3$, $\beta>\beta_s(q)$) inherit a persistence-style bottleneck, or does coherence remove it?*

## Connections

- **Direct positive companion (same authors).** Woodard–Schmidler–Huber 2009b, *Conditions for rapid mixing of parallel and simulated tempering on multimodal distributions*, Ann. Appl. Probab. 19:617–640 (arXiv:0906.2341): matching *upper* bound on the PT/ST gap under persistence-failure + width-ratio control. Together, the two papers form a near-dichotomy for tempering on multimodal targets.
- **`madras-zheng-2003-swapping-algorithm`** — the *positive* result this paper structurally bounds: PT mixes polynomially on mean-field **Ising** (which lacks persistence), but fails on mean-field **Potts** ($q\ge 3$, which has it). The pair is the natural juxtaposition for the Review chapter: same algorithm, different model, opposite outcome — and the explanation is structural (persistence), not algorithmic.
- **`cuff-et-al-2012-mean-field-potts-glauber`** — Glauber on the *same* model (mean-field Potts $q\ge 3$) is exponentially slow above the spinodal $\beta_s(q)$. Woodard–Schmidler–Huber complements this by ruling out the PT/ST workaround. Together they lock in mean-field Potts as a classical low-T no-go cell.
- **`bovier-den-hollander-2015-metastability`** — Eyring–Kramers $K e^{\beta\Gamma^*}$ for *local* low-T classical samplers (Glauber, Kawasaki, PCA). Woodard–Schmidler–Huber extends the no-go to the standard *non-local* classical workaround; together they leave only cluster algorithms and lifted/non-reversible chains as plausible classical escapes — and Gore–Jerrum 1999 closes Swendsen–Wang on the same model.
- **Bhatnagar–Randall 2004** (*Torpid mixing of simulated tempering on the Potts model*): the original PT-torpidity result on mean-field Potts; subsumed and generalised by Woodard–Schmidler–Huber. **No existing review in the corpus** ([CHECK] if it should be added Tier-2).
- **Within the thesis grid.** Cited in `classical-gibbs-sampling-corpus.md` §5 (low-T no-go) and in the **Classical Ising/Potts × low-T** Tier-1 cell. Quantum-side counterpart to define the comparison: Chen et al. 2023/2025 KMS Lindbladians and Ding et al. 2024 low-T quantum Gibbs samplers — none have published numerics on mean-field Potts in the persistence regime; the thesis numerics chapter is the natural place to start.
