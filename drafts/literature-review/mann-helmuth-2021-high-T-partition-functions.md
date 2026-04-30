---
paper: "Mann, Helmuth (2021)"
title: "Efficient Algorithms for Approximating Quantum Partition Functions"
arxiv: "2004.11568"
year: 2021
venue: "J. Math. Phys. 62:022201"
pdf: "n/a"

temperature: [high-T]
commutativity: [noncommuting]
locality: [geometric-local, any-d, sparse]
particle-statistics: [spin]
hamiltonian-models: [generic-local-spin, Heisenberg, transverse-field-Ising]

paradigm: [classical-deterministic, cluster-expansion, polymer-model]
quantum-or-classical: [classical]

result-type: [FPTAS, polynomial-time-deterministic, complexity-upper]
key-scaling: "deterministic FPTAS for $Z_G(\\beta)$ in time $\\mathrm{poly}(n, 1/\\varepsilon)$ on any graph $G$ of max degree $\\Delta$ whenever $|\\beta| \\le 1/(e^4 \\Delta)$; threshold $|\\beta| = O(1/\\Delta)$ optimal under RP $\\ne$ NP"

related: ["bakshi-liu-moitra-tang-2024-high-T-unentangled", "kuwahara-kato-brandao-2020-cmi-clustering", "mann-helmuth-2023-low-T-partition-functions"]
---

# Efficient Algorithms for Approximating Quantum Partition Functions — Mann, Helmuth (2021)

**One-sentence takeaway**: For every quantum spin Hamiltonian with pairwise interactions on a graph of maximum degree $\Delta$, the partition function $Z_G(\beta)$ admits a *deterministic* fully polynomial-time approximation scheme whenever $|\beta| \le 1/(e^4 \Delta)$ — sharpening Harrow–Mehraban–Soleimanifar from $n^{\mathrm{polylog}(n)}$ to honest $\mathrm{poly}(n, 1/\varepsilon)$ in the high-temperature regime, and matching the Sly–Sun-style $\Theta(1/\Delta)$ hardness threshold up to constants.

## Setting

- **Object computed**: the quantum partition function $Z_G(\beta) = \mathrm{tr}(e^{-\beta H})$ of an $n$-spin Hamiltonian $H = \sum_{\{u,v\} \in E(G)} h_{uv}$ with one local term per edge of an interaction graph $G = (V,E)$.
- **Hamiltonian class**: pairwise (2-local) interactions on a graph of maximum degree $\Delta$. Each $h_{uv}$ is a bounded-norm Hermitian operator on the two qudits $u,v$ — generic *non-commuting* terms (Heisenberg, $XYZ$, transverse-field Ising, etc.). Strictly more general than geometric locality on $\mathbb{Z}^d$ (any $D$, any periodic / aperiodic / random connectivity with bounded degree).
- **Temperature**: high-T, $|\beta| \le 1/(e^4 \Delta)$, *independent of system size*. Note this is stated for *complex* $\beta$ — the algorithm gives a zero-free disk, hence Taylor truncation and Barvinok-style interpolation.
- **Access model**: classical input (the Hamiltonian as a list of bounded-norm 2-local terms) → classical $\varepsilon$-multiplicative approximation of $\log Z_G(\beta)$.

## Main Results

- **FPTAS for the partition function (Thm 1.1 / Cor.)**: there is a deterministic algorithm that, on input $G$ of max degree $\Delta$, $\beta \in \mathbb{C}$ with $|\beta| \le 1/(e^4 \Delta)$, and $\varepsilon > 0$, outputs $\widehat Z$ with $|\widehat Z / Z_G(\beta) - 1| \le \varepsilon$ in time $\mathrm{poly}(n, 1/\varepsilon)$. Constant exponent in $n$ (not $\mathrm{polylog}(n)$); explicit polynomial dependence on $1/\varepsilon$. [CHECK: exact exponent of $n$ depends on how the polymer truncation depth is set; the leading term is $n \cdot (\log n / \log(1/(\Delta|\beta|)))^{O(1)}$.]
- **Free energy / log-Z approximation**: equivalent additive-error FPTAS for $\log Z_G(\beta)$ at the same threshold, by truncating the cluster expansion of $\log Z$ to $O(\log(n/\varepsilon))$ terms and Barvinok-evaluating each.
- **Optimality of the $1/\Delta$ scaling**: the threshold $|\beta| = \Theta(1/\Delta)$ is tight up to constants; past it, even *classical* antiferromagnetic Ising on $\Delta$-regular graphs is NP-hard to approximate (Sly 2010; Sly–Sun 2014), and quantum hardness inherits.
- **Subclass corollaries**: applies to the antiferromagnetic Heisenberg model, transverse-field Ising, $XXZ$, $XYZ$, etc. on any bounded-degree interaction graph in any spatial dimension.

## Method

The algorithm is the cluster-expansion / abstract-polymer pipeline of Helmuth–Perkins–Regts (STOC 2019), instantiated for the *quantum* cluster expansion of Netočný–Redig: write $\log Z_G(\beta) = \sum_{\Gamma} \phi(\Gamma) \prod_{\gamma \in \Gamma} w(\gamma)$ as a sum over connected sub-multigraphs $\Gamma$ of $G$ with explicit polymer weights $w(\gamma)$ depending on the local Hamiltonian terms. The Kotecký–Preiss convergence condition holds whenever $|\beta| \le 1/(e^4 \Delta)$, giving exponential decay in polymer size; truncating to polymers of size $O(\log(n/\varepsilon))$ and enumerating connected subgraphs of bounded size in $\mathrm{poly}(n)$ time (Patel–Regts) yields the FPTAS. The "simple and slightly sharper" analysis vs Harrow–Mehraban–Soleimanifar comes from working directly in the polymer formalism with the optimal Kotecký–Preiss criterion rather than via complex-zero localisation arguments — this swaps the quasi-polynomial truncation depth $\mathrm{polylog}(n)$ for a *constant-depth-times-$O(\log n)$* enumeration.

## Quantum vs Classical

- **This paper is purely classical** — it computes a quantum quantity ($Z = \mathrm{tr}\, e^{-\beta H}$) on a classical computer. No quantum algorithm is proposed.
- **Predecessor**: Harrow–Mehraban–Soleimanifar 2020 [arXiv:1910.09071] — quasi-polynomial $n^{\mathrm{polylog}(n)}$ algorithm above a similar threshold $|\beta| = O(1/\Delta)$, via complex-zero localisation à la Barvinok.
- **Gap to predecessor**: $n^{\mathrm{polylog}(n)} \to \mathrm{poly}(n)$ — a *super-polynomial-to-polynomial* improvement in $n$ at the *same* temperature regime. Same $\beta$-threshold up to constants; the improvement is purely in the analysis (polymer Kotecký–Preiss vs zero-localisation).
- **Quantum baseline (downstream)**: any quantum Gibbs sampler at $|\beta| < 1/(e^4 \Delta)$ produces samples from $\rho \propto e^{-\beta H}$, from which $Z$ can be estimated. After Mann–Helmuth, *any* such quantum approach must beat $\mathrm{poly}(n)$ in the same regime to claim advantage on partition functions — and beat the Bakshi–Liu–Moitra–Tang depth-1 product-state preparation to claim advantage on the *state*.
- **Source of the difference**: at $|\beta| = O(1/\Delta)$ the imaginary-time propagator $e^{-\beta H}$ is a perturbation of the identity in a precise polymer-norm sense; non-commutativity of the local terms enters only through the operator norms appearing in the Netočný–Redig weights, and is fully absorbed by the Kotecký–Preiss radius. No structural quantum effect (entanglement, sign problem, KMS condition) plays any role at this temperature.
- **Caveat**: the bound is on $Z$, not on $\rho$ or on observables. Estimating $\langle O \rangle = \mathrm{tr}(O e^{-\beta H})/Z$ for a local $O$ is then standard (ratio of two FPTAS evaluations); estimating *non-local* observables, sampling from $\rho$ in a particular basis, or preparing a quantum copy of $\rho$ is not directly given here — those are the Bakshi–Liu–Moitra–Tang and Rouzé–França–Alhambra contributions.

## Implications for Quantum Advantage

- **Regime cell**: noncommuting × high-T × geometric-local (any $D$ via bounded-degree graph) × spin. Tier-1 anchor in the high-T column of the corpus 3×3 grid for *all* local Hamiltonians.
- **What this changes**: this paper *closes* the high-T row for the partition-function task in the strongest sense — there is a deterministic, polynomial-time, $\beta$-threshold-tight classical algorithm. Any quantum Gibbs sampler that needs $\mathrm{poly}(n)$ time to estimate $Z$ in this regime is not an *advantage*, only a parity. Together with Bakshi–Liu–Moitra–Tang (separability of the state itself) and Kuwahara–Kato–Brandão (CMI clustering of the quantum Markov network), the high-T regime is now structurally classical along three independent axes: partition function (here), state structure (BLMT), correlation structure (KKB).
- **Promising or not**: **negative** for advantage in the high-T regime. Reinforces the corpus framing that quantum vs classical separation must move to **intermediate or low temperature**, where polymer convergence fails (correlation length diverges), Barvinok zero-localisation breaks down (zeros approach the real $\beta$ axis at the phase transition), and tensor networks pay $\exp(O(\beta))$ bond dimension.
- This is the canonical *high-T classical baseline* for the thesis Review chapter — every quantum algorithm targeting $Z$ in the high-T regime is benchmarked against this $\mathrm{poly}(n, 1/\varepsilon)$ runtime.

## Open Questions / Limitations

- **Constant in the threshold**: $1/(e^4 \Delta)$ is not the physical critical temperature; the gap to $\beta_c^{\mathrm{phys}}$ remains. Whether polymer methods can be sharpened to *the* phase transition (rather than $\Theta(1/\Delta)$) is an open question for both classical and quantum sides.
- **Low temperature**: not covered. Below the polymer convergence threshold the cluster expansion diverges; the companion paper Mann–Helmuth 2023 [arXiv:2201.06533] handles low-$T$ for *stable* quantum spin systems (dominant ground state, discrete order parameter), but the *intermediate / near-critical* regime remains genuinely hard and is the chapter's punchline cell.
- **Gapless models**: no explicit handling of models with vanishing spectral gap — the polymer expansion only sees $\beta \cdot \|H\|$ via $\Delta$, not the gap, so it works wherever it converges, but gives no information about gapless / critical physics.
- **2D phase transitions in physically interesting models** (e.g. 2D quantum Heisenberg antiferromagnet, transverse-field Ising at $g_c$) sit in the regime *not* covered — convergence fails before $T_c$ in any non-trivial example.
- **Pairwise restriction**: Theorem stated for 2-local interactions; extension to general $k$-local (with $\Delta$ replaced by an analogue counting term-overlaps) is straightforward in principle but is not the form proven in this paper. [CHECK: the Helmuth–Perkins–Regts framework handles general bounded-multiplicity polymer models, and Bakshi–Liu–Moitra–Tang use the $k$-local version explicitly.]
- **Sampling vs counting**: gives a counting FPTAS, not a sampling algorithm. The standard counting-to-sampling reduction (Jerrum–Valiant–Vazirani) is not invoked here, and the BLMT 2024 product-state sampler is the natural companion.

## Connections

- **Direct predecessor — Harrow–Mehraban–Soleimanifar 2020** [arXiv:1910.09071]: same regime, quasi-polynomial time. This paper is the polynomial-time sharpening; the gap is super-poly-to-poly with no temperature improvement.
- **State-side companion — Bakshi–Liu–Moitra–Tang 2024** (`bakshi-liu-moitra-tang-2024-high-T-unentangled`): proves the *state* is separable in roughly the same regime ($\beta < 1/(c\, \mathfrak{d}\,\mathfrak{K})$); together they completely classicalise the high-T row (partition function + state preparation).
- **Correlation-side companion — Kuwahara–Kato–Brandão 2020** (`kuwahara-kato-brandao-2020-cmi-clustering`): proves CMI clustering in the same regime, the structural reason cluster / polymer methods converge.
- **Low-T companion — Mann–Helmuth 2023** (`mann-helmuth-2023-low-T-partition-functions`): polymer methods extended to the deep symmetry-broken phase under stability assumptions; this paper (2021) is the high-T half of the pair.
- **Hardness baseline — Sly 2010, Sly–Sun 2014**: matching $\Theta(1/\Delta)$ classical hardness for antiferromagnetic Ising; sets the constant-factor target for the threshold in this paper.
- **Polymer machinery — Helmuth–Perkins–Regts 2019** (STOC): the abstract-polymer-to-FPTAS reduction (Patel–Regts connected-subgraph enumeration + Barvinok truncation) used here verbatim.
- **Quantum cluster expansion — Netočný–Redig 2004**: defines the polymer weights $w(\gamma)$ for quantum spin systems; this is the quantum input to the otherwise-classical Helmuth–Perkins–Regts pipeline.

Sources:
- [arXiv:2004.11568 — Efficient Algorithms for Approximating Quantum Partition Functions](https://arxiv.org/abs/2004.11568)
- [J. Math. Phys. 62:022201 (2021)](https://pubs.aip.org/aip/jmp/article/62/2/022201/234532/Efficient-algorithms-for-approximating-quantum)
- [Author's hosted PDF — Ryan Mann](https://www.ryanmann.org/media/publications/efficient-algorithms-for-approximating-quantum-partition-functions.pdf)
