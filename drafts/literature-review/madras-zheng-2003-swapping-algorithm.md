---
paper: "Madras, Zheng (2003)"
title: "On the swapping algorithm"
arxiv: "n/a"
year: 2003
venue: "Random Struct. Alg. 22:66–97"
pdf: "n/a"

temperature: [low-T]
commutativity: [commuting]
locality: [mean-field]
particle-statistics: [spin]
hamiltonian-models: [mean-field-Ising, Curie-Weiss, exponential-valley]

paradigm: [classical-MCMC, parallel-tempering]
quantum-or-classical: [classical]

result-type: [mixing-time-upper, gap-bound]
key-scaling: "Gap(swap) >= poly(1/n) for symmetric bimodal targets, including the Curie-Weiss Ising model at any fixed beta below or above beta_c, provided the temperature ladder is geometric with ratio 1+c/n and includes O(n) replicas — turning t_mix from exp(Omega(n)) (single-site Glauber at low T) into poly(n)"

related: ["bovier-den-hollander-2015-metastability", "cuff-et-al-2012-mean-field-potts-glauber", "woodard-schmidler-huber-2009-torpid-tempering"]
---

# On the Swapping Algorithm — Madras & Zheng (2003)

**One-sentence takeaway**: The canonical *positive* result for parallel tempering: on the symmetric bimodal mean-field (Curie-Weiss) Ising model — where single-site Glauber is exponentially slow at low temperature — the swapping chain with a geometric temperature ladder of ratio $1+\Theta(1/n)$ and $\Theta(n)$ replicas mixes in time $\mathrm{poly}(n)$, sharp on both the spacing requirement and the symmetry requirement (Bhatnagar-Randall and Woodard-Schmidler-Huber show the same algorithm becomes torpid for $q\ge 3$ Potts).

## Setting

- **Configuration space.** Spins $\sigma \in \{-1,+1\}^n$ on the complete graph (mean-field Ising / Curie-Weiss); the second example is a one-dimensional symmetric "exponential valley" target $\pi(x) \propto e^{-\beta|x|}$ on a bounded interval, used as a clean test bed for the same proof technique.
- **Target measure.** Curie-Weiss Ising $\mu_\beta(\sigma) = Z_\beta^{-1}\exp\!\big(\beta n^{-1}\sum_{i<j}\sigma_i\sigma_j\big)$ at inverse temperature $\beta$. The static phase transition is at $\beta_c = 1$; at $\beta>\beta_c$ the measure is bimodal between the all-$+$ and all-$-$ phases with an $O(\sqrt n)$ free-energy barrier in the magnetisation $m=n^{-1}\sum_i\sigma_i$, so single-site Metropolis / Glauber takes $e^{\Theta(n)}$ steps to cross between modes (this is the prototype Eyring–Kramers low-T classical lower bound; see `bovier-den-hollander-2015-metastability`).
- **Algorithm — swapping (parallel tempering).** Choose a ladder of inverse temperatures $0 = \beta_0 < \beta_1 < \cdots < \beta_M$. The state space is the product $\prod_{k=0}^M \{-1,+1\}^n$, with stationary measure $\bigotimes_k \mu_{\beta_k}$. One step of the chain alternates:
  1. **Local move** — perform a single-site Metropolis step independently within each replica.
  2. **Swap move** — pick a neighbouring pair $(k,k{+}1)$ and propose to swap the configurations of replicas $k$ and $k{+}1$, accepting with the Metropolis probability for the joint product measure.
- **Question.** Spectral gap (and hence mixing time) of the swap chain on the *coldest* replica $\beta_M$, as a function of $n$, $\beta_M$, the ladder spacing, and the number of replicas $M$.

## Main Results

The paper's results are stated as gap (and decomposition) bounds, not explicit mixing-time formulas; the polynomial mixing time for the coldest replica then follows from $t_{\mathrm{mix}} = O(\mathrm{Gap}^{-1}\log(1/\pi_*))$ with $\log(1/\pi_*) = O(n)$.

- **Theorem 1.1 / 2.x — Curie-Weiss Ising rapid mixing.** For the swap chain on the mean-field Ising model with a geometric temperature ladder $\beta_{k+1}/\beta_k = 1 + c/n$ for some absolute constant $c>0$, and $M = \Theta(n\log\beta_M)$ replicas, the spectral gap satisfies
  $$ \mathrm{Gap}(P_{\mathrm{swap}}) \;\ge\; \mathrm{poly}(1/n), $$
  uniformly in $\beta_M$ in any compact range. Hence the total-variation mixing time of the coldest replica is $t_{\mathrm{mix}} = \mathrm{poly}(n)$. **[CHECK]** the precise polynomial degree — the original paper gives an explicit exponent $c'$ (likely $n^{O(1)}$ with $c' \le 6$) but I have not been able to verify it from the paywalled Wiley source; subsequent literature (Woodard–Schmidler–Huber rapid-mixing companion 2009) sharpens it.
- **Theorem — exponential-valley rapid mixing.** Same conclusion ($\mathrm{Gap} \ge \mathrm{poly}(1/n)$) for the symmetric exponential-valley target $\pi(x) \propto e^{-\beta|x|}$ on $[-1,+1]$ discretized on $n$ points, again with a $1+\Theta(1/n)$-geometric ladder.
- **Decomposition lemma (the key technical tool).** Madras–Randall-style Markov-chain decomposition: write the stationary measure as a mixture over disjoint subsets $A_i$, prove a "restricted" gap on each $A_i$ and a "projection" gap between them, then combine. Madras–Zheng prove a clean decomposition statement for *product* state spaces: if the swap chain has a polynomial gap when restricted to each "sign sector" (all replicas in the same metastable basin) *and* the swap moves between sectors at the high-temperature top of the ladder mix with polynomial probability, the global gap is polynomial. They also prove $\mathrm{Gap}(P) \ge n^{-1}\,\mathrm{Gap}(P^{\otimes n})$ — a generic tensorisation comparison used as a building block in the decomposition (this little lemma is the most-cited statement of the paper outside its main theorem).

The key *qualitative* content is two-fold: (i) **sharp ladder spacing** — the $1+c/n$ geometric ratio is *necessary* in the sense that with ratio $1+\omega(1/n)$ adjacent-replica overlap collapses and the swap acceptance probability becomes exponentially small; with ratio $1+o(1/n)$ the ladder needs $\omega(n)$ replicas, increasing the cost. The factor $1/n$ matches the natural fluctuation scale of the magnetisation under $\mu_{\beta_k}$; (ii) **symmetry is essential** — the proof relies on the $\sigma\leftrightarrow -\sigma$ symmetry of the Curie-Weiss measure, which makes the two metastable basins equivalent so that the projection chain on $\{+,-\}$ mixes in $O(1)$ time once both sides have polynomial restricted gaps.

## Method

Markov-chain decomposition (Madras–Randall 2002). Partition each replica's state space into the $+$-basin and $-$-basin of the magnetisation (boundary at $m=0$); restricted to a single basin, single-site Metropolis on Curie-Weiss Ising has $\mathrm{Gap}_{\mathrm{rest}} \ge \mathrm{poly}(1/n)$ at any $\beta$ (the basin is unimodal). The swap moves at high temperature ($\beta_0 = 0$, where $\mu_{\beta_0}$ is uniform) have $\Theta(1)$ probability of straddling the $m=0$ boundary, so propagating "the right basin" up the ladder via swaps gives the projection chain a polynomial gap. Combining the restricted gaps and the projection gap via the Madras–Randall decomposition formula yields a polynomial global gap. The whole argument is a *clean* low-T MCMC proof: it exploits exactly the fact that a hot replica sees a unimodal landscape and feeds the cold replica through a sequence of slow-but-still-bimodal intermediate temperatures that overlap pairwise with $O(1)$ probability.

## Quantum vs Classical

- **Baseline (this paper itself is the classical baseline).** The single-site Glauber chain on the same Hamiltonian has $t_{\mathrm{mix}} = e^{\Theta(n)}$ at any $\beta>\beta_c$ (Eyring–Kramers, `bovier-den-hollander-2015-metastability` Thm 16.5; for Curie-Weiss Ising the barrier is $\Gamma^* = \Theta(n)$). The swap chain with a properly chosen ladder reduces this to $\mathrm{poly}(n)$ — an *exponential classical-on-classical* speedup. So this paper is not a quantum-advantage candidate; it is a *classical algorithmic improvement* that closes one of the most natural quantum-advantage gaps before quantum methods enter the picture.
- **Source of the difference (classical vs classical).** Non-locality of moves. Glauber's locality forces it to climb the $\Theta(n)$ free-energy barrier sequentially. The swap move couples replicas at different $\beta$; effectively the cold replica "borrows" the high-$T$ replica's ability to cross $m=0$ in $O(1)$ time, modulo the cost of propagating the result back through the ladder. The structural ingredient is *symmetry*: the projection chain on basin labels $\{+,-\}$ has spectral gap $\Theta(1)$ only because both basins have equal weight; once that symmetry breaks (e.g. $q\ge 3$ Potts: the disordered phase is *not* symmetric to the $q$ ordered phases at coexistence), the projection chain itself becomes the bottleneck.
- **Sharpness — the failure for $q\ge 3$ Potts.** Bhatnagar–Randall 2004 (mean-field Potts $q\ge 3$ simulated tempering torpid) and `woodard-schmidler-huber-2009-torpid-tempering` (general no-go: "persistence of a narrow mode" criterion applies) prove that the *same algorithm* is torpid on mean-field Potts: $t_{\mathrm{mix}} \ge e^{\Omega(n)}$ for $\beta\in(\beta_s,\beta_c)$ regardless of the temperature ladder. The Madras–Zheng positive result is therefore *sharp to ferromagnetic Ising-style symmetric bimodality*, not a generic feature of replica exchange. This is the most important contextual fact: parallel tempering does *not* generally resolve first-order transitions, and the Madras–Zheng theorem is the canonical example of when it does.
- **Mean-field is essential.** The proof uses that the magnetisation alone is a sufficient statistic — the projection chain has $O(1)$ states (just $\{+,-\}$). On a lattice $\mathbb{Z}^d$, $d\ge 2$, the analogous argument would have to project onto domain-wall configurations and the "projection chain" itself acquires geometric complexity; no analogue of this theorem is known in finite dimension, and the low-T lattice Glauber crossing time is governed by Wulff-construction surface tension (`bovier-den-hollander-2015-metastability` §17.4), which interacts non-trivially with replica exchange.
- **Quantum side: open candidate cell.** The paper makes no quantum statement. The natural cross-axis question for the thesis: a *quantum* Gibbs sampler (KMS Lindbladian, Davies generator) on the same Curie-Weiss Ising Hamiltonian — does it match $\mathrm{poly}(n)$ without needing a ladder? The Hamiltonian is classical (commuting) so KMS-Lindbladian convergence here reduces to classical jump rates with quantum-coherent-superposition structure on the Bohr-frequency side; whether that buys anything over Madras–Zheng's $\mathrm{poly}(n)$ is unclear. **The relevant comparison for a quantum advantage claim on commuting low-T mean-field is therefore not Glauber's $e^{\Theta(n)}$ but Madras–Zheng's $\mathrm{poly}(n)$.** This is the corpus's main reason to flag the paper as Tier-1.

## Implications for Quantum Advantage

- **Regime cell (thesis 3×3 grid).** **Classical Ising/Potts × low-T**, mean-field locality, spin. Within the 3×3 grid this paper sits in the same cell as `cuff-et-al-2012-mean-field-potts-glauber` and `woodard-schmidler-huber-2009-torpid-tempering` and refines it: the cell is *not uniformly hard classically* — the symmetric bimodal subcell (Curie-Weiss Ising) admits a $\mathrm{poly}(n)$ classical algorithm, the asymmetric bimodal subcell ($q\ge 3$ Potts at coexistence) does not.
- **What this changes (for the chapter).** Raises the classical bar from Glauber's $e^{\beta\Gamma^*}$ to parallel tempering's $\mathrm{poly}(n)$ in the *symmetric bimodal* subcell. Any quantum-advantage claim on Curie-Weiss Ising must therefore beat $\mathrm{poly}(n)$ — *not* $e^{\Omega(n)}$ — and must remain efficient *without* using a temperature ladder (or else it competes head-on with classical PT, which already has $\mathrm{poly}(n)$).
- **Promising or not (for quantum advantage).** **Not promising in this subcell.** Once classical PT achieves $\mathrm{poly}(n)$ on Curie-Weiss Ising, the room for quantum advantage shrinks to a polynomial-in-$n$ factor at best, which is unlikely to be detectable past constants. The interesting quantum-advantage question shifts to the *adjacent* cell where PT is torpid: mean-field Potts $q\ge 3$ (Bhatnagar–Randall, Cuff et al., WSH). There the classical baseline is back to $e^{\Omega(n)}$ and a quantum sampler that handles asymmetric metastability — without relying on a tempering ladder it has no analogue of — would constitute a genuine separation. So Madras–Zheng *narrows* the candidate region for quantum advantage but *sharpens its boundary*: it is the asymmetry / first-order character that defeats classical, not bimodality per se.

## Open Questions / Limitations

- **Mean-field only.** No analogue is known on $\mathbb{Z}^d$, $d\ge 2$. The natural conjecture — PT mixes in $\mathrm{poly}(n)$ on lattice Ising at any $\beta\ne\beta_c$ — is open; even rapid mixing of plain Swendsen–Wang on lattice Ising at all $\beta$ is recent (Guo–Jerrum 2018) and uses a different machinery.
- **Symmetric targets only.** Even on the complete graph the proof breaks for the Curie-Weiss Ising model with non-zero external field $h\ne 0$ (asymmetric basins): the projection chain on $\{+,-\}$ acquires a $1-e^{-\Theta(\beta nh)}$ asymmetry, and one of the basins becomes a metastable trap. Woodard–Schmidler–Huber 2009 (rapid mixing) gives a partial extension to "nearly symmetric" multimodal targets but with quantitative degradation in the asymmetry parameter.
- **Polynomial degree not optimised.** The bound is $\mathrm{poly}(n)$ but the explicit exponent is large and almost certainly not tight; subsequent work (Bhatnagar–Randall, Woodard–Schmidler–Huber, Ge et al. 2018, Tawn–Roberts 2019) sharpens both the exponent and the dependence on $\beta$.
- **No quantum analogue.** Whether a noncommutative analogue of replica exchange (e.g. coupling KMS Lindbladians at different $\beta$ via a quantum swap channel) exists and gives any advantage on classical or quantum Hamiltonians is open.

## Connections

- **Paired no-go (mean-field Potts, asymmetric bimodal).** `woodard-schmidler-huber-2009-torpid-tempering` — the direct companion: same algorithm, same model class, opposite conclusion as soon as the $\sigma\leftrightarrow-\sigma$ symmetry is broken (mean-field Potts $q\ge 3$, asymmetric Gaussian mixtures, narrow modes). Together with Bhatnagar–Randall 2004 this delineates the *exact* boundary of where parallel tempering helps classically: symmetric bimodal yes, asymmetric / first-order no.
- **Static (Glauber) baseline.** `cuff-et-al-2012-mean-field-potts-glauber` — on Potts $q\ge 3$ Glauber is $e^{\Theta(n)}$ in the coexistence interval $(\beta_s,\beta_c)$; combined with the WSH no-go, *no known classical method* mixes in $\mathrm{poly}(n)$ on this model. Madras–Zheng's positive Ising result by contrast shows that on the *symmetric* bimodal version of the same physics (Ising = Potts $q=2$, where $\beta_s=\beta_c$ and the transition is continuous so there is no coexistence interval and no first-order metastability anyway — see Cuff et al. for the Ising-vs-Potts distinction), PT achieves $\mathrm{poly}(n)$.
- **Eyring–Kramers anchor.** `bovier-den-hollander-2015-metastability` — defines the $K e^{\beta\Gamma^*}$ classical Glauber baseline that Madras–Zheng beats by switching from local to ladder-coupled moves. The Bovier–den Hollander treatment of Curie-Weiss in Ch. 13–15 makes the barrier $\Gamma^*$ explicit; Madras–Zheng's theorem says that with a $1+c/n$-spaced ladder you don't pay for $\Gamma^*$ at all in the symmetric case.
- **Decomposition technical lineage.** Madras–Randall 2002 (decomposition framework, cited as Tier-2 in the corpus §4) is the immediate predecessor — Madras–Zheng's contribution is the *application* of decomposition to a tempered product state space and the identification of the ladder-spacing/symmetry conditions under which the technique gives polynomial bounds.
- **Modern sharpening.** Woodard–Schmidler–Huber 2009 (rapid mixing companion, arXiv:0906.2341) gives quantitative two-sided gap bounds on tempering for general "nearly symmetric multimodal" targets; Ge–Lee–Risteski 2018 (continuous-state extensions); Tawn–Roberts 2019 (optimal scaling). The corpus tracks these as Tier-2 references; Madras–Zheng remains the foundational Tier-1 entry because it sets the comparison.
- **Within the thesis grid.** Cited in `classical-gibbs-sampling-corpus.md` §4 (low-T metastability / tempering — the *positive* anchor; WSH 2009 is the *negative* anchor in the same row). The mean-field Ising × low-T cell has a unique structural feature: classical PT closes the gap, so quantum advantage in this cell would have to clear a $\mathrm{poly}(n)$ classical bar, not the naive $e^{\Omega(n)}$ Glauber bar. This is the kind of fact the thesis Review chapter must surface explicitly when constructing its regime map.

Sources consulted:
- [On the swapping algorithm — Madras 2003 — Random Structures & Algorithms (Wiley)](https://onlinelibrary.wiley.com/doi/abs/10.1002/rsa.10066)
- [Sufficient Conditions for Torpid Mixing of Parallel and Simulated Tempering — Woodard, Schmidler, Huber 2009 (Project Euclid)](https://projecteuclid.org/journals/electronic-journal-of-probability/volume-14/issue-none/Sufficient-Conditions-for-Torpid-Mixing-of-Parallel-and-Simulated-Tempering/10.1214/EJP.v14-638.full)
- [Conditions for rapid mixing of parallel and simulated tempering — Woodard, Schmidler, Huber 2009 (arXiv:0906.2341)](https://arxiv.org/abs/0906.2341)
- [Torpid Mixing of Simulated Tempering on the Potts Model — Bhatnagar, Randall 2004](https://randall.math.gatech.edu/r-tempsoda.pdf)
- [Simulated Tempering and Swapping on Mean-Field Models — J. Stat. Phys. 2016](https://link.springer.com/article/10.1007/s10955-016-1526-8)
