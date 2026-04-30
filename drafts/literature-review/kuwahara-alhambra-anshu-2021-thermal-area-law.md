---
paper: "Kuwahara, Alhambra, Anshu (2021)"
title: "Improved Thermal Area Law and Quasilinear Time Algorithm for Quantum Gibbs States"
arxiv: "2007.11174"
year: 2021
venue: "PRX 11:011047"
pdf: "n/a"

temperature: [all]
commutativity: [non-commuting]
locality: [1D, geometric-local]
particle-statistics: [spin]
hamiltonian-models: [generic-1d-local-spin]

paradigm: [classical-tensor-network, mps, area-law]
quantum-or-classical: [classical]

result-type: [mixing-time-upper, structural, quasi-linear-algorithm, area-law-bound]
key-scaling: "Generic-lattice EoF / Renyi entanglement of purification across a region $\\mathrm{EoF}, R^{(\\alpha)}_{\\mathrm{purif}} = \\widetilde{O}(\\beta^{2/3} |\\partial A|)$; in 1D, MPO bond dimension $\\chi = \\exp\\!\\bigl(\\sqrt{\\widetilde{O}(\\beta\\log(n/\\varepsilon))}\\bigr)$; classical MPS-construction runtime $\\widetilde{O}(n)$ at $\\beta = o(\\log n)$, sub-polynomial in $n/\\varepsilon$."

related: ["czarnik-dziarmaga-2015-peps-finite-T", "molnar-schuch-verstraete-cirac-2015-pepo-gibbs", "kuwahara-kato-brandao-2020-cmi-clustering", "bakshi-liu-moitra-tang-2024-high-T-unentangled"]
---

# Improved Thermal Area Law and Quasilinear Time Algorithm for Quantum Gibbs States — Kuwahara, Alhambra, Anshu (2021)

**One-sentence takeaway**: For *any* finite-range local Hamiltonian on a lattice the thermal state obeys an entanglement-of-formation / Rényi-entanglement-of-purification area law with $\widetilde{O}(\beta^{2/3})$ temperature dependence — a strict improvement over the Wolf–Verstraete–Hastings $O(\beta)$ bound — and in 1D this powers a *quasi-linear* $\widetilde{O}(n)$ classical algorithm for constructing an MPS representation of $\rho_\beta$ at every temperature with $\beta = o(\log n)$.

## Setting

- **Object**: Gibbs state $\rho_\beta = e^{-\beta H}/Z_\beta$ of a finite-range local Hamiltonian $H = \sum_X h_X$ on a $D$-dimensional lattice with bounded local terms $\|h_X\| \le 1$. For the algorithmic part, $D=1$ and $n$ spins.
- **Temperature**: any $\beta>0$ for the area-law bound; the quasi-linear *algorithm* requires $\beta = o(\log n)$, i.e. up to logarithmic-in-$n$ inverse temperature.
- **Quantities controlled**: entanglement of formation $\mathrm{EoF}(\rho_{AB})$, Rényi entanglement of purification $R^{(\alpha)}_{\mathrm{purif}}$, and operator entanglement (Schmidt rank of $\rho_\beta$ across a cut), as a function of the size of the boundary $|\partial A|$ between region $A$ and its complement.
- **Output**: MPS / MPO classical description of $\rho_\beta$ to trace-distance $\varepsilon$, suitable for computing local expectation values $\mathrm{tr}(O \rho_\beta)$ classically in time linear in $n$ up to polylogs.

## Main Results

- **Improved thermal area law (generic lattice)**: for any region $A$ in a $D$-dimensional lattice, the entanglement of formation and Rényi entanglement of purification of $\rho_\beta$ across $\partial A$ obey
  $$ \mathrm{EoF}(\rho_{A:A^c}),\ R^{(\alpha)}_{\mathrm{purif}}(\rho_{A:A^c}) \;\le\; \widetilde{O}\!\bigl(\beta^{2/3}\bigr) \cdot |\partial A|, $$
  improving the Wolf–Verstraete–Hastings–Cirac (2008) mutual-information bound from $O(\beta)$ to $\widetilde{O}(\beta^{2/3})$. Suggests *sub-ballistic / diffusive* propagation of correlations under imaginary-time evolution, in contrast to ballistic real-time Lieb–Robinson spreading.
- **Operator-entanglement / MPO bond dimension in 1D**: for $D=1$ and trace-distance error $\varepsilon$, $\rho_\beta$ admits an MPO approximation with bond dimension
  $$ \chi \;=\; \exp\!\Bigl(\sqrt{\widetilde{O}\bigl(\beta \log(n/\varepsilon)\bigr)}\Bigr), $$
  which is *sub-polynomial* in $n/\varepsilon$ for $\beta = o(\log n)$. [CHECK: paper states $\chi$ smaller than any power of $n/\varepsilon$ in this regime; the explicit $\sqrt{\widetilde{O}(\beta\log(n/\varepsilon))}$ form follows by inverting $\log\chi^2 = \widetilde{O}(\beta\log(n/\varepsilon))$.]
- **Quasi-linear-time MPS construction (1D)**: a classical, deterministic algorithm constructs an MPS / MPO representation of $\rho_\beta$ with the above bond dimension in time $\widetilde{O}(n) \cdot \mathrm{poly}(\chi)$, hence $\widetilde{O}(n)$ — *quasi-linear in the system size* — whenever $\beta = o(\log n)$ and $\varepsilon = 1/\mathrm{poly}(n)$. First rigorous quasi-linear classical Gibbs-sampler for 1D quantum spin chains at all temperatures up to logarithmic-in-$n$.

## Method

The improvement from $O(\beta)$ to $\widetilde{O}(\beta^{2/3})$ rests on a low-degree polynomial approximation of the imaginary-time propagator $e^{-\beta H}$ that maps the $\beta$-dependence of operator spreading onto a *random walk* of length $\widetilde{O}(\beta^{2/3})$ rather than $O(\beta)$ — a Markov-property bootstrap that exposes the diffusive (rather than ballistic) nature of imaginary-time correlation transport. The quasi-linear algorithm uses a new *block decomposition* of $\rho_\beta$ into local pieces, structurally analogous to Haah–Hastings–Kothari–Low's (FOCS '18) decomposition of real-time evolution; each block is computed independently and stitched into an MPO via standard truncation, exploiting the new sub-polynomial bond bound to keep the per-block cost polylogarithmic.

## Quantum vs Classical

- **Baseline classical**:
  - Wolf–Verstraete–Hastings–Cirac (2008): mutual-information thermal area law $I(A:A^c) \le O(\beta)|\partial A|$.
  - Hastings (2006): cluster-expansion MPO approximation in 1D, polynomial-in-$n$ bond dimension at *gapped* temperatures.
  - Molnar–Schuch–Verstraete–Cirac (2015): PEPO with bond dimension $\exp(O(\beta))$ in any dimension — algorithmic but with worse bond growth.
- **Baseline quantum**: dissipative KMS Gibbs samplers (Chen–Kastoryano–Gilyén 2023; Ding–Li–Lin 2024) and dissipative low-T preparation (Ding–Lin 2024). For 1D spin chains at $\beta = o(\log n)$, the Chen et al. high-temperature KMS sampler runs in $\widetilde{O}(n^c)$ with constant $c \ge 1$ on a quantum computer with $\mathrm{poly}(\beta, 1/\varepsilon)$ depth.
- **Gap**: *no quantum advantage* in the 1D, all-T, $\beta = o(\log n)$ regime. Classical Kuwahara-Alhambra-Anshu achieves $\widetilde{O}(n)$ runtime *fully classically* — quasi-linear and likely better than any quantum algorithm in raw cost (the quantum samplers must still simulate the Lindbladian to convergence). The "classical baseline cost" against which any 1D quantum claim is benchmarked is set here.
- **Source of the difference**: the 1D area-law structure (low boundary $|\partial A|$) combined with the diffusive imaginary-time bound makes $\rho_\beta$ classically representable with *sub-polynomial* bond dimension. None of the quantum-specific obstructions (sign problem, KMS-condition implementation cost, noncommuting Lindblad operators) bite in 1D with bounded $\beta$, so the classical method dominates structurally.
- **Caveat**: the result is genuinely 1D — the $D$-dimensional area law gives EoF / Rényi-purification only across $|\partial A|$, but in $D \ge 2$ this still implies bond dimension $\exp(\widetilde{O}(\beta^{2/3} L^{D-1}))$ across an $L \times \cdots \times L$ region, which is *not* sub-polynomial. The quasi-linear algorithm does not lift to 2D.

## Implications for Quantum Advantage

- **Regime cell**: 1D × non-commuting × geometric-local × spin × $\beta = o(\log n)$. This paper *closes* the 1D row of the Tier-1 grid for any inverse temperature growing slower than $\log n$.
- **What this changes**: it sets the strongest rigorous classical baseline for 1D quantum Gibbs sampling. Any 1D quantum-advantage claim — e.g. for Chen et al. 2023, Ding et al. 2024, or Capel et al. 2025 — must either (i) target the regime $\beta = \Theta(\log n)$ or higher (low-T 1D), (ii) target $D \ge 2$ where the $\beta^{2/3}$ bound no longer suppresses bond dimension to sub-polynomial, or (iii) argue a polylog speedup over $\widetilde{O}(n)$ via specific structure (e.g. observable locality, low-depth circuit access). None of these is automatic.
- **Promising or not**: **negative** for 1D quantum advantage at moderate $\beta$. Any honest quantum-vs-classical separation in 1D Gibbs sampling must live at *super-logarithmic* $\beta$ — exactly the low-T 1D regime that this paper explicitly does not cover and that the corpus identifies as one frontier.

## Open Questions / Limitations

- **Sub-logarithmic $\beta$ only**: the algorithm requires $\beta = o(\log n)$. At $\beta = \Theta(\log n)$ or higher in 1D the bond dimension $\chi = \exp(\sqrt{\widetilde{O}(\beta\log n)})$ becomes super-polynomial, and the runtime ceases to be quasi-linear. Whether a different rigorous classical method handles 1D at $\beta = \mathrm{polylog}(n)$ is open. [CHECK: low-T 1D Gibbs sampling at $\beta = n^{O(1)}$ is presumed hard classically — this paper does not pin the threshold.]
- **No 2D / higher-$d$ algorithm**: the area-law improvement holds in $D \ge 2$, but $|\partial A| = \Theta(L^{D-1})$ kills the sub-polynomial-bond-dimension consequence. The natural 2D analogue (PEPO with bond $\exp(\widetilde{O}(\beta^{2/3}))$ rather than $\exp(O(\beta))$) is not constructed, and a 2D analogue of the block decomposition is open.
- **Tightness of $\beta^{2/3}$**: the $2/3$ exponent reflects a particular polynomial approximation; whether the true diffusive scaling is $\beta^{1/2}$ (genuine random-walk $\sqrt{t}$) or stays at $\beta^{2/3}$ is open. Lower bounds on entanglement of formation as a function of $\beta$ are not given.
- **Mutual information vs entanglement of formation**: the paper bounds EoF and Rényi-purification, both upper bounds on quantum mutual information. The original Wolf–Verstraete–Hastings–Cirac mutual-information bound $O(\beta)|\partial A|$ is *not* improved here — the strict gain is on the more entanglement-theoretic quantities. [CHECK: relationship to MI in their Theorem 1 / Theorem 2 statements.]

## Connections

- **Czarnik–Dziarmaga (2015)** [related: `czarnik-dziarmaga-2015-peps-finite-T`]: 2D *heuristic* PEPS imaginary-time analogue. KAA's 1D quasi-linear algorithm is the *rigorous* sibling; whether the diffusive $\beta^{2/3}$ bound can be combined with the PEPS-with-ancilla machinery to give a rigorous 2D algorithm is the natural open follow-up.
- **Molnar–Schuch–Verstraete–Cirac (2015)** [related: `molnar-schuch-verstraete-cirac-2015-pepo-gibbs`]: rigorous PEPO with bond $\exp(O(\beta))$ in any $D$. KAA improves the $\beta$-exponent in 1D from $\beta$ to $\beta^{2/3}$ inside a square root, and crucially constructs a *quasi-linear* algorithm rather than just an existence theorem.
- **Kuwahara–Kato–Brandão (2020)** [related: `kuwahara-kato-brandao-2020-cmi-clustering`]: above-threshold-$T$ CMI clustering and Markov-network structure for $\rho_\beta$ on any lattice. KAA's algorithm covers a strictly larger 1D temperature range (up to $\beta = o(\log n)$, i.e. potentially below the CMI threshold) by exploiting 1D MPS structure rather than CMI.
- **Bakshi–Liu–Moitra–Tang (2024)** [related: `bakshi-liu-moitra-tang-2024-high-T-unentangled`]: high-T separability + classical sampling on *any-dimensional* lattice but only at constant-$\beta$. KAA is the 1D complement — it loses dimensional generality but reaches *all $\beta = o(\log n)$*, including below any high-T separability threshold.
- **Quantum-side counterparts in the corpus**: Chen et al. (2023, 2025), Ding et al. (2024), Capel et al. (2025) — these are exactly the 1D / general-$D$ KMS Gibbs samplers whose 1D performance must be benchmarked against this paper's $\widetilde{O}(n)$ classical bar. In the thesis Review chapter this paper is the canonical "1D × all-T-up-to-logarithmic" classical entry.
- **Haah–Hastings–Kothari–Low (FOCS '18)**: real-time evolution block decomposition. KAA's imaginary-time block decomposition is structurally analogous; the analogy is the load-bearing technical idea, and is worth flagging in the Review chapter as an example of cross-pollination between Hamiltonian-simulation and Gibbs-sampling techniques.
