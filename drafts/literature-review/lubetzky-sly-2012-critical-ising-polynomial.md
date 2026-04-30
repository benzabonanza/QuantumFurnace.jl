---
paper: "Lubetzky, Sly (2012)"
title: "Critical Ising on the square lattice mixes in polynomial time"
arxiv: "1001.1613"
year: 2012
venue: "Comm. Math. Phys. 313:815-836"
pdf: "n/a"

temperature: [intermediate-T]
commutativity: [commuting]
locality: [geometric-local]
particle-statistics: [spin]
hamiltonian-models: [Ising]

paradigm: [classical-MCMC]
quantum-or-classical: [classical]

result-type: [mixing-time-upper, gap-bound]
key-scaling: "1/gap = O(n^c) and t_mix = O(n^c) at beta = beta_c on Z^2 finite box of side n, for an unspecified absolute polynomial c; replaces the previously open exponential / quasi-polynomial gap with a polynomial bound. Restricted to Z^2."

related: ["bauerschmidt-dagallier-2024-near-critical-ising-lsi", "eldan-koehler-zeitouni-2022-spectral-condition-ising", "lubetzky-sly-2013-cutoff-ising-lattice", "cuff-et-al-2012-mean-field-potts-glauber"]
---

# Critical Ising on the square lattice mixes in polynomial time — Lubetzky & Sly 2012

**One-sentence takeaway**: At the critical temperature $\beta = \beta_c$ on a finite box of $\mathbb{Z}^2$, single-site Glauber dynamics for the Ising model has spectral gap and total-variation mixing time that are *polynomial* in the side length $n$, closing the last open temperature for 2D Ising and ruling out the exponential / sub-exponential slowdown that worst-case mean-field heuristics suggest at criticality.

## Setting

- **Sampled object**: ferromagnetic nearest-neighbour Ising measure $\mu_\Lambda^{\beta,\xi}(\sigma) \propto \exp\big(\beta \sum_{x \sim y} \sigma_x \sigma_y\big)$ on $\Lambda = [1,n]^2 \subset \mathbb{Z}^2$, at $\beta = \beta_c = \tfrac12 \log(1+\sqrt{2})$.
- **Boundary conditions**: free, all-plus, all-minus, and periodic, all covered by the theorem; results stable to arbitrary boundary conditions [CHECK — paper covers "any boundary conditions" per §1.2].
- **Hamiltonian class**: classical commuting 2-spin ferromagnet, geometric-local on $\mathbb{Z}^2$. Two dimensions and the *square* lattice are essential to the proof (FK / random-cluster scaling-limit input is dimension-2 specific).
- **Dynamics**: continuous-time single-site heat-bath / Glauber on $\Lambda$. Number of sites $N = n^2$; the paper states bounds in terms of side length $n$.

## Main Results

- **Theorem 1 (main)**: at $\beta = \beta_c$ on $\Lambda = [1,n]^2 \subset \mathbb{Z}^2$, the spectral gap of the Glauber dynamics satisfies $\mathrm{gap} \ge n^{-c}$ for an absolute constant $c < \infty$, hence the total-variation mixing time obeys $t_{\mathrm{mix}} = O(n^c)$ uniformly over the listed boundary conditions. The exponent $c$ is **not given explicitly** in the abstract / theorem statement — it is a sufficiently large absolute constant produced by the proof, *not* the conjectured physics value. [CHECK — explicit numerical exponent]
- Companion bound for the **critical antiferromagnetic** 2D Ising model with the same scaling (by spin-flip symmetry on the bipartite $\mathbb{Z}^2$).
- The paper explicitly *does not* identify the dynamic critical exponent $z$ — the conjectured value from physics is $z \approx 2.1665$ (Nightingale–Blöte and others), and Glauber on $\mathbb{Z}^2$ at $\beta_c$ is expected to satisfy $t_{\mathrm{mix}} \asymp n^z$, but $c$ produced here is far from this conjectured tight value.

## Method

Combines (i) Smirnov / Duminil-Copin / Chelkak–Smirnov-era understanding of the **scaling limit of the critical FK (random-cluster, $q=2$) representation** — in particular conformal-invariance-driven RSW-type crossing estimates and bounded-density spatial mixing of FK clusters at $\beta_c$ — with (ii) classical Markov-chain machinery (block dynamics, canonical paths / decomposition of Glauber on $\Lambda$ into block updates, comparison with the FK heat bath). Concretely: critical-FK spatial mixing gives a polynomial bound on the influence of boundary conditions deep inside $\Lambda$, which is then converted, via a block-dynamics decomposition à la Martinelli–Olivieri (which would give $O(\log n)$ above $T_c$ but breaks at $T_c$), into a polynomial gap for single-site Glauber. The crucial point is that the FK input survives at $\beta_c$, where strong spatial mixing fails, *because* it is a polynomial — not exponential — decay of crossing probabilities.

## Quantum vs Classical

- **Baseline**: this paper *is* the classical baseline for 2D Ising at criticality. Prior to it, polynomial gap was open at $\beta_c$ in 2D — only $\beta < \beta_c$ (Martinelli–Olivieri, log-Sobolev with $O(\log n)$ constant) and $\beta > \beta_c$ (Cesi–Guadagni–Martinelli–Schonmann, $\exp(\Omega(n))$ for $\pm$ phases) were rigorous, and the abstract notes the gap was "rigorous everywhere except at criticality". The corresponding *upper bound* from physics conjectures is $t_{\mathrm{mix}} \asymp n^{2 + z} \approx n^{4.17}$ counting in number of sites $N = n^2$, but this is heuristic.
- **Gap (vs hypothetical quantum)**: none — the regime is already polynomial classically, and a quantum KMS-DB Lindbladian for the 2D Ising Hamiltonian (which is classical and commuting) reduces to a classical sampler modulo block-encoding overhead. There is no sign problem, no noncommutativity, and the model is geometric-local.
- **Source of the (non-)difference**: ferromagnetic + commuting + classical — the FK / Edwards–Sokal joint representation gives a positive-coefficient cluster expansion (Swendsen–Wolff being the algorithmic instance), and conformal invariance of the scaling limit further sharpens estimates. None of the structural obstructions that motivate quantum samplers (sign problem, KMS condition for noncommuting $H$, Lindblad-locality issues) appear here.
- **Caveat**: dimension and lattice are not transferable. The proof uses 2D-specific conformal-invariance / RSW results; in $d \ge 3$ at $\beta_c$ the analogous polynomial bound for Glauber is *open*, and Bauerschmidt–Dagallier 2024 only treats $d \ge 5$ (mean-field regime) and stays *near*, not *at*, $\beta_c$ for $d \le 4$. So while this paper closes the 2D-square cell, it leaves open the higher-dimensional critical cells classically — which still does not open a quantum-advantage door, since stoquastic / commuting reductions apply.

## Implications for Quantum Advantage

- **Regime cell**: classical Ising × intermediate-T (specifically *at* criticality) × geometric-local × $\mathbb{Z}^2$ (square lattice).
- **What this changes**: closes the upper-bound side of the 2D classical Ising × critical cell with an explicit polynomial bound. Together with Bauerschmidt–Dagallier 2024 (functional inequality, $d \ge 5$, near-critical) and Lubetzky–Sly 2013 (cutoff in the high-T phase), the classical-side row "ferromagnetic Ising on $\mathbb{Z}^d$ at intermediate-T" is now polynomially controlled in the cases that matter for the chapter — including the historically worst case ($d=2$, $\beta=\beta_c$) where critical slowing-down was the textbook obstruction.
- **Promising or not**: explicitly *not promising* for quantum advantage. Even at criticality, even with $z \approx 2.17$ slowdown, single-site Glauber on 2D Ising mixes in polynomial time. A quantum Gibbs sampler for the same Hamiltonian (which is diagonal in the computational basis) would inherit a classical reduction at best and pay block-encoding / OFT overhead at worst — there is no asymptotic regime in this cell where quantum could beat classical, and the chapter punchline cell is non-stoquastic × low-T, well away from here. The cell is a "classical wins or ties" anchor.
- **Honest cost**: the polynomial exponent $c$ produced by the proof is not the physical $z$ — it is much larger. So this paper *opens* the cell rather than tightly *closes* it; tightness up to $c = z$ remains an open problem. But the qualitative classical-wins verdict is robust to the value of $c$.

## Open Questions / Limitations

- **Explicit polynomial exponent**: the paper produces a sufficiently large absolute $c$, not the conjectured physical $z \approx 2.1665$. Sharpening $c \to z$ for Glauber on $\mathbb{Z}^2$ at $\beta_c$ is open.
- **Higher dimensions**: the analogous polynomial gap at $\beta_c$ is open for $d \ge 3$. Bauerschmidt–Dagallier handles $d \ge 5$ in the near-critical (not at-critical) sense via mean-field $\chi_\beta$ bounds; $d \in \{3,4\}$ at $\beta_c$ remains unsolved.
- **Other 2D lattices**: the proof leans on FK / conformal-invariance inputs proven for the *square* lattice; triangular and honeycomb at $\beta_c$ are presumably similar but require the corresponding RSW / scaling-limit inputs.
- **No log-Sobolev / cutoff statement here**: the paper bounds the spectral gap and mixing time polynomially but does not give a functional inequality (LSI / MLSI) at $\beta_c$, nor does it identify cutoff. Bauerschmidt–Dagallier 2024 gives the LSI route in higher $d$; Lubetzky–Sly 2013 gives cutoff above $\beta_c$.
- **Antiferromagnetic / spin-glass / non-FKG analogues**: the FK input requires the random-cluster representation and FKG; non-FK extensions (spin glasses, signed couplings) are out of scope.

## Connections

- **Bauerschmidt–Dagallier 2024** [`bauerschmidt-dagallier-2024-near-critical-ising-lsi`]: complementary near-critical result via Polchinski-flow LSI, uniform in volume up to $\beta_c$ in $d \ge 5$ — gives the entropic functional-inequality version where this paper gives the spectral-gap version, in different dimensions. Together they cover (square-lattice 2D $@$ $\beta_c$) and (mean-field-regime $d \ge 5$ near $\beta_c$); $d \in \{3,4\}$ and other lattices remain open.
- **Eldan–Koehler–Zeitouni 2022** [`eldan-koehler-zeitouni-2022-spectral-condition-ising`]: high-T spectral gap from $\|J\|_{\mathrm{op}} < 1$ via stochastic localization. EKZ controls the *strict high-T* regime ($\|J\|_{\mathrm{op}} < 1$, not up to $\beta_c$); Lubetzky–Sly handles exactly the boundary point that EKZ does not reach.
- **Lubetzky–Sly 2013 (cutoff)** [`lubetzky-sly-2013-cutoff-ising-lattice`]: sister paper proving cutoff for Glauber on $\mathbb{Z}^d$ Ising in the *high-T* phase ($\beta < \beta_c$) with mixing time $\frac{1}{2\lambda_\infty} \log n$. Together with this paper, the high-T and critical 2D cells are both classically rigorous; low-T is ruled out by Cesi–Guadagni–Martinelli–Schonmann (exponential mixing in the broken-symmetry phase with mixed boundary conditions).
- **Cuff et al. 2012** [`cuff-et-al-2012-mean-field-potts-glauber`]: mean-field analogue (Curie–Weiss Potts Glauber) where critical slowdown is $n^{4/3}$ (not $n^z$ with $z\approx 2$); shows the 2D-lattice exponent is geometry-driven and substantially larger than mean-field.
- **Within the corpus**: this is the canonical "even at criticality, classical Glauber on 2D Ising mixes in poly time" anchor in the chapter's regime grid. It rules out the strongest classical-side hope a quantum advocate might cite (critical slowing-down) and forces the chapter's quantum-advantage punchline cell into the **non-stoquastic** × **low-T** corner, where no analogous classical polynomial result exists.
