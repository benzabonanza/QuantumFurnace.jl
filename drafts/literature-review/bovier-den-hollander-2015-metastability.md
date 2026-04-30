---
paper: "Bovier and den Hollander (2015)"
title: "Metastability: A Potential-Theoretic Approach"
arxiv: "n/a"
year: 2015
venue: "Springer Grundlehren der mathematischen Wissenschaften 351"
pdf: "supplementary-informations/classical-review/Bovier–den Hollander-book.pdf"

temperature: [low-T]
commutativity: [n/a]
locality: [local, geometric-local, mean-field]
particle-statistics: [spin, lattice-gas]
hamiltonian-models: [Ising, Curie-Weiss, lattice-gas, zero-range]

paradigm: [classical-MCMC]
quantum-or-classical: [classical]

result-type: [mixing-time-upper, mixing-time-lower, structural]
key-scaling: "E_m[tau_s] = K e^{beta Gamma^*} (1+o(1)) as beta -> infty (Eyring–Kramers); for 2D Glauber-Ising at field h, Gamma^* = J[4 ell_c] - h[ell_c (ell_c-1)+1] with ell_c = ceil(2J/h)"

related: ["lubetzky-sly-2012-critical-ising-polynomial", "lubetzky-sly-2013-cutoff-ising-lattice", "woodard-schmidler-huber-2009-torpid-tempering", "madras-zheng-2003-swapping-algorithm", "sly-2010-computational-transition-uniqueness", "cuff-et-al-2012-mean-field-potts-glauber"]
---

# Metastability: A Potential-Theoretic Approach — Bovier & den Hollander (2015)

**One-sentence takeaway**: Comprehensive monograph that develops the *potential-theoretic* / capacity framework for metastability and uses it to derive sharp Eyring–Kramers asymptotics $\mathbb{E}_{\mathbf{m}}[\tau_{\mathbf{s}}] = K\,e^{\beta\Gamma^*}(1+o(1))$ — including the explicit Arrhenius prefactor $K$ — for low-temperature reversible Markov dynamics on Ising-like models, fixing the canonical *classical* low-T Gibbs-sampling timescale.

## Setting

- **Configuration space.** Finite $\Lambda \subset \mathbb{Z}^d$ (lattice site set), single-spin space $\Upsilon$ (e.g. $\{-1,+1\}$ for Ising, $\{0,1\}$ for lattice gas); $S = \Upsilon^{\Lambda}$.
- **Target measure.** Gibbs measure $\mu_\beta(\xi) = Z_\beta^{-1} e^{-\beta H(\xi)}$ at inverse temperature $\beta$, with $H$ a local lattice Hamiltonian. The flagship example (Ch. 17) is the 2D ferromagnetic Ising model with external field $h$:
  $$ H(\sigma) = -\tfrac{J}{2}\!\!\sum_{\{x,y\}\in\Lambda^*}\!\!\sigma(x)\sigma(y) - \tfrac{h}{2}\sum_{x\in\Lambda}\sigma(x), \qquad J,h>0. $$
- **Dynamics (access model).** Classical reversible Markov process with single-site moves and Metropolis rates $c_\beta(\xi,\xi') = e^{-\beta[H(\xi')-H(\xi)]_+}$ (Glauber spin-flip in Ch. 17, Kawasaki particle-hop in Ch. 18); heat-bath rates and probabilistic cellular automata are also covered (Sec. 16.4). Each move is $O(1)$ work.
- **Low-T regime.** Asymptotics are taken in $\beta \to \infty$ at *fixed* finite volume (Part VI) or with $|\Lambda_\beta|\to\infty$ slowly enough that $|\Lambda_\beta|\,e^{-\beta\Gamma^*}\to 0$ (Part VII). Curie-Weiss / mean-field at fixed sub-critical $T$ in large volume is treated in Part V.
- **Question.** Mean and law of the *crossover time* $\tau_{\mathbf{s}}$ from a metastable configuration $\mathbf{m}$ (e.g. all-minus) to a stable configuration $\mathbf{s}$ (e.g. all-plus), and the corresponding spectral gap of the generator $-\mathcal{L}_\beta$. This is the Gibbs-sampling mixing time when the chain starts in the wrong phase.

## Main Results

The book proves three "universal" theorems for any low-T Metropolis dynamics satisfying minimal energy-landscape hypotheses (H1 unique metastable/stable pair, H2 protocritical-set uniformity); see Sec. 16.1. Let $\Gamma^* = \Phi(\mathbf{m},\mathbf{s}) - H(\mathbf{m})$ be the *communication height* (energy barrier) between $\mathbf{m}$ and $\mathbf{s}$, and let $\mathscr{C}^*$ be the *critical set* (bottleneck saddles).

- **Theorem 16.4 (gate and uniform entrance).** $\lim_{\beta\to\infty} \mathbb{P}_{\mathbf{m}}(\tau_{\mathscr{C}^*} < \tau_{\mathbf{s}} \mid \tau_{\mathbf{s}} < \tau_{\mathbf{m}}) = 1$, and entrance into $\mathscr{C}^*$ is asymptotically uniform: $\lim_{\beta\to\infty}\mathbb{P}_{\mathbf{m}}(\xi_{\tau_{\mathscr{C}^*}} = \chi) = 1/|\mathscr{C}^*|$ for every $\chi\in\mathscr{C}^*$.

- **Theorem 16.5 (Eyring–Kramers / mean crossover time).** There exists $K \in (0,\infty)$ such that
  $$ \lim_{\beta\to\infty} e^{-\beta\Gamma^*}\,\mathbb{E}_{\mathbf{m}}[\tau_{\mathbf{s}}] = K. $$
  Equivalently, $\mathbb{E}_{\mathbf{m}}[\tau_{\mathbf{s}}] = K e^{\beta \Gamma^*}(1+o(1))$. The barrier $\Gamma^*$ is *model-independent in form* (an energy-landscape quantity) and the prefactor $K$ is given by the variational formula $K = 1/\Theta$ in Lemma 16.17, depending on the geometry of the protocritical/critical sets.

- **Theorem 16.6 (spectrum and exponential law).** $\lim_{\beta\to\infty}\lambda_\beta\,\mathbb{E}_{\mathbf{m}}[\tau_{\mathbf{s}}] = 1$ where $\lambda_\beta$ is the smallest non-zero eigenvalue of $-\mathcal{L}_\beta$, and $\lim_{\beta\to\infty}\mathbb{P}_{\mathbf{m}}(\tau_{\mathbf{s}}/\mathbb{E}_{\mathbf{m}}[\tau_{\mathbf{s}}] > t) = e^{-t}$ for all $t\geq 0$.

For the 2D Ising / Glauber model in a finite torus with $h \in (0, 2J)$ and critical droplet size $\ell_c = \lceil 2J/h\rceil$ (Sec. 17.1), Theorems 17.3–17.4 identify $\Gamma^*$ and $K$ explicitly:
$$ \Gamma^* = J[\,4\ell_c\,] - h\bigl[\ell_c(\ell_c-1)+1\bigr], \qquad K(\Lambda) = \frac{3}{4(2\ell_c-1)}\,\frac{1}{|\Lambda|}. $$
The protocritical set $\mathscr{P}^*$ is the family of $(\ell_c-1)\times \ell_c$ quasi-square droplets and the critical set $\mathscr{C}^*$ adds a single $1\times 1$ protuberance. The 3D extension (Sec. 17.6) gives an explicit recursion for $\Gamma^*$ but only partial information on $\mathscr{C}^*, K$. Kawasaki dynamics in 2D (Ch. 18) gives an analogous picture with droplet size $\ell_c = \lceil U/(2U-\Delta)\rceil$ where $U$ is the binding energy and $\Delta$ the activation energy.

For *large volumes* $|\Lambda_\beta|\to\infty$ with $|\Lambda_\beta|\,e^{-\beta\Gamma^*}\to 0$ (Theorem 19.4), homogeneous nucleation gives
$$ \mathbb{E}_{\nu}[\tau_{\mathscr{S}^c}] = \frac{1}{N_1 |\Lambda_\beta|}\,e^{\beta\Gamma^*}(1+o(1)), \qquad N_1 = 4\ell_c, $$
i.e. the volume factor speeds nucleation by $|\Lambda_\beta|$ (any of $|\Lambda_\beta|$ sites can host the critical droplet) but the dominant exponential timescale $e^{\beta\Gamma^*}$ is unchanged.

## Method

The potential-theoretic approach (Bovier–Eckhoff–Gayrard–Klein 2001) reduces metastable hitting times and the small eigenvalues of the generator to *capacities* of disjoint configuration sets, computed via the Dirichlet variational principle (upper bounds via test potentials) and the Berman–Konsowa flow principle (lower bounds via test flows); see Ch. 7 and Sec. 9.1. The mean exit time is essentially $\mathbb{E}_{\mathbf{m}}[\tau_{\mathbf{s}}] = \mu_\beta(A(\mathbf{m}))/\text{cap}_\beta(\mathbf{m},\mathbf{s})\,(1+o(1))$ (Lemma 16.14, eq. 16.2.28), and the small eigenvalues of $-\mathcal{L}_\beta$ coincide with the inverse mean exit times up to $1+o(1)$ (Theorem 8.34, 8.45). Sharp matching upper and lower bounds on the capacity reduce to a low-dimensional variational problem on the protocritical/critical set, evaluated combinatorially in each model.

## Quantum vs Classical

- **Baseline.** This monograph is the canonical *classical* low-T Gibbs-sampling baseline for Ising-like models. It establishes that any local reversible classical sampler (Glauber, heat-bath, PCA, Kawasaki) on a low-T Hamiltonian with energy barrier $\Gamma^*$ has mixing time $\Theta(e^{\beta\Gamma^*})$, exponential in $\beta\Gamma^*$, with a sharp prefactor.

- **Gap.** No quantum-vs-classical comparison appears: "Throughout the book we only consider classical stochastic dynamics. It would be interesting to consider quantum stochastic dynamics as well, but this is beyond the scope of the book." (Preface, p. vii). The corresponding *quantum* lower bound is unknown for stoquastic Hamiltonians at low T, and the best quantum upper bound for non-stoquastic low-T Hamiltonians (Chen et al. KMS-Lindbladians, Ding et al., Kastoryano–Brandão) is also $e^{O(\beta)}$ in the relevant regimes — so no proven separation either way exists today. The structural reason a separation *might* exist: classical Glauber crosses the barrier sequentially (one spin flip at a time, energy fluctuations of size $\Gamma^*$ are exponentially rare), whereas a coherent quantum walk on the energy landscape could in principle tunnel through the saddle in time $\sim e^{c\sqrt{\beta\Gamma^*}}$ (Grover-style speedup) or $\mathrm{poly}(\beta,\Gamma^*)$ (adiabatic evolution through the gap). Whether KMS detailed-balance Lindbladians actually realize such tunneling at low T is the open question this thesis addresses.

- **Source of the difference.** Energy barriers vs. noncommutativity. The classical bottleneck is *geometric*: the metastable basin is separated from the stable basin by a saddle of height $\Gamma^*$ above $H(\mathbf{m})$ in a configuration-space landscape whose paths are local moves. The quantum bottleneck for a KMS Gibbs sampler is the *Lindbladian spectral gap*, which is governed by Bohr-frequency overlap and noncommutative log-Sobolev / detectability-lemma quantities — different mathematical objects whose low-T behaviour is a priori unrelated to $\Gamma^*$.

- **Caveat.** The Eyring–Kramers timescale $e^{\beta\Gamma^*}$ is *information-theoretic* in the sense that it survives any smart classical post-processing of the local-move chain (Theorem 16.5 covers Glauber, heat-bath and PCA uniformly): a generic classical sampler has no way to "tunnel" across $\Gamma^*$ without spending $e^{\beta\Gamma^*}$ time. Tempering, simulated annealing or cluster moves can break this barrier on specific models (e.g. ferromagnetic Ising via Swendsen–Wang, Madras–Zheng PT) but each method has known torpid-mixing counterexamples (Bhatnagar–Randall, Woodard–Schmidler–Huber, Gore–Jerrum). The Bovier–den Hollander result is thus the *unconditional* low-T classical baseline modulo non-local algorithmic surgery.

## Implications for Quantum Advantage

- **Regime cell.** **Low-T × n/a (commutativity not applicable on the classical side) × geometric-local × spin** (and mean-field via Curie-Weiss in Ch. 13–15; lattice-gas via Kawasaki in Ch. 18). The low-T entry of the corpus's 3×3 Hamiltonian × temperature grid.

- **What this changes.** Defines the *unconditional classical baseline* against which any quantum low-T Gibbs-sampling claim must be benchmarked: any quantum sampler that mixes in time $o(e^{\beta\Gamma^*})$ on a low-T Ising-like model with explicit $\Gamma^*$ given by Theorem 17.3 would constitute a genuine quantum speedup over the optimal classical local sampler. Conversely, a quantum mixing-time *lower* bound matching $e^{\beta\Gamma^*}$ — currently unknown — would close the door on quantum advantage in this cell.

- **Promising or not.** *Mildly promising, not proven.* The exponential-in-$\beta\Gamma^*$ classical timescale is a *real* obstruction with sharp constants, not a worst-case artefact, so the room for quantum speedup in this cell is large. But all known quantum Gibbs samplers (KMS Lindbladians, Davies generators, dissipative engineering) currently inherit some form of low-T slowdown, and no unconditional quantum lower bound exists either. This is exactly the cell where the thesis's KMS-Lindbladian numerics matter most.

## Open Questions / Limitations

- Restricted to *reversible* classical Markov dynamics (Preface); non-reversible variational principles (Sec. 7.4) are weaker. The monograph does not address irreversible / lifted classical samplers (Suwa–Todo) which are sometimes faster.
- Mean-field large-volume models (Curie-Weiss, Ch. 13) are treated at *fixed* sub-critical $T$, not in the joint $T\to 0, N\to\infty$ limit relevant to spin-glass / quantum annealing benchmarks.
- The 3D Ising case (Sec. 17.6) only identifies $\Gamma^*$ explicitly; the prefactor $K_{d=3}$ is given recursively in terms of $K_{d=2}$ and the number $M_{d=3}$ of quasi-cubes inside a critical droplet, but its analytic form is open.
- Kac-type long-range interactions (Sec. 15.6 bibliographical notes) and disordered random-field Curie-Weiss in continuous distributions (Ch. 15) have only partial results.
- Glauber dynamics at small magnetic field $h \to 0$ jointly with $\beta \to \infty$ (Wulff-construction regime) and crystallisation in continuum particle systems are listed as Challenges in Ch. 22–23.
- *No quantum analogue.* The potential-theoretic toolbox does not extend to quantum Lindbladians; capacities of operator-valued Dirichlet forms in noncommutative $L^2$ are not developed here.

## Connections

- **Same low-T cell, alternative samplers.** `madras-zheng-2003-swapping-algorithm` (parallel tempering polynomial mixing on mean-field Ising — concrete *positive* algorithm beating the $e^{\beta\Gamma^*}$ baseline); `woodard-schmidler-huber-2009-torpid-tempering` (general no-go for tempering — confirms PT does not always escape Bovier–den Hollander timescale); `cuff-et-al-2012-mean-field-potts-glauber` (sharp three-regime mixing for Glauber on mean-field Potts, $e^{\Theta(n)}$ in the low-T cell, complementary to Curie-Weiss treatment in Ch. 13).
- **Hardness/no-go.** `sly-2010-computational-transition-uniqueness` (computational hardness at the uniqueness threshold — phase-transition obstruction at the *intermediate*-T side of this row).
- **Cutoff on the "easy" side.** `lubetzky-sly-2013-cutoff-ising-lattice` and `lubetzky-sly-2012-critical-ising-polynomial` cover Glauber-Ising at high and critical $T$ respectively; together with Bovier–den Hollander these three references span the temperature axis for Glauber-Ising (high-T cutoff at $\Theta(\log n)$ — critical polynomial — low-T $e^{\beta\Gamma^*}$).
- **Potential-theoretic underpinnings.** Bovier–Eckhoff–Gayrard–Klein 2001/2002 (foundational paper introducing capacity-based metastability, cited as [33], [34] throughout the book) and Eyring 1935 / Kramers 1940 (origin of the $K e^{\beta E_a}$ Arrhenius formula).
- **Quantum side counterpart this defines the baseline for.** Chen et al. 2023/2025 (KMS Lindbladians), Ding et al. 2024 (low-T quantum Gibbs samplers), Kastoryano–Brandão (cluster-expansion classical/quantum at high T) — comparing the $e^{\beta\Gamma^*}$ classical baseline to their quantum mixing-time bounds is the central numerical question in the thesis Review chapter.
