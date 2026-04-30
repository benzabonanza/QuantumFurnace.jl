---
paper: "Bakshi, Liu, Moitra, Tang (2024)"
title: "High-Temperature Gibbs States are Unentangled and Efficiently Preparable"
arxiv: "2403.16850"
year: 2024
venue: "Quantum (2024) / arXiv 2403.16850"
pdf: "supplementary-informations/classical-review/bakshi-liu-moitra-tang-2024-high-T-unentangled.pdf"

temperature: [high-T]
commutativity: [noncommuting]
locality: [k-local, geometric-local, sparse]
particle-statistics: [spin]
hamiltonian-models: [generic-local-qubit, generic-local-qudit]

paradigm: [classical-other, structural]
quantum-or-classical: [classical, comparison]

result-type: [structural, complexity-upper, no-go]
key-scaling: "separable for $\\beta < 1/(c\\,\\mathfrak{d}\\,\\mathfrak{K})$; classical preparation in $\\widetilde{O}(n^{6+\\log\\mathfrak{d}/\\log(\\beta_c/\\beta)} \\log^3(1/\\varepsilon) \\cdot \\mathrm{poly}(\\mathfrak{d}, \\mathfrak{K}))$ for $\\beta < \\beta_c = 1/(\\gamma\\mathfrak{d}\\mathfrak{K})$"

related: ["mann-helmuth-2021-high-T-partition-functions", "kuwahara-kato-brandao-2020-cmi-clustering", "rouze-franca-alhambra-2024", "harrow-mehraban-soleimanifar-2020", "anshu-arunachalam-kuwahara-soleimanifar-2024"]
---

# High-Temperature Gibbs States are Unentangled and Efficiently Preparable — Bakshi, Liu, Moitra, Tang (2024)

**One-sentence takeaway**: Above a constant temperature independent of $n$, the Gibbs state of any low-intersection local Hamiltonian is *separable* (a classical mixture of product states) and admits a classical $\mathrm{poly}(n)$-time sampling algorithm whose output is prepared by a depth-1 quantum circuit — closing the high-temperature regime to quantum advantage.

## Setting

- **Object sampled**: Gibbs state $\rho = e^{-\beta H}/\mathrm{tr}(e^{-\beta H})$ of an $n$-qudit Hamiltonian $H = \sum_a H_a$.
- **Hamiltonian class**: $(\mathfrak{d}, \mathfrak{K})$-*low-intersection* — each term $H_a$ acts on at most $\mathfrak{K}$ sites (locality), and in the dual interaction graph each term overlaps at most $\mathfrak{d}$ others (degree). Strictly generalises geometrically local on a $D$-dim lattice, where $(\mathfrak{d}, \mathfrak{K}) = (\mathfrak{K}(2D)^{\mathfrak{K}-1}, \mathfrak{K})$. Term norms $\|H_a\|_\mathrm{op} \le 1$; Pauli-coefficient terms with $|\lambda_a| \le 1$.
- **Temperature**: high-T, $\beta < \beta_c$ for an explicit $\beta_c = \Theta(1/(\mathfrak{d}\mathfrak{K}))$, *independent of system size*.
- **Access model**: classical description of $H$ in; classical description of a stabilizer product state out (then prepared in depth 1 by single-qubit gates if a quantum copy is wanted).

## Main Results

- **Separability (Thm 1.5 / 4.20)**: for $\beta < 1/(100\,\mathfrak{d}\,\mathfrak{K})$, $\rho$ is separable and expressible as a mixture over single-qubit *stabilizer* product states $A_j \in \{\tfrac12(I \pm \sigma_x), \tfrac12(I \pm \sigma_y), \tfrac12(I \pm \sigma_z)\}$. "Sudden death of thermal entanglement" at a system-size-independent temperature.
- **Efficient preparation (Thm 1.7)**: for $\beta < \beta_c = 1/(\gamma \mathfrak{d}\mathfrak{K})$ with universal constant $\gamma$, a classical randomised algorithm outputs a product-state description $\widehat{\rho}$ with $\|\rho - \mathbb{E}[\widehat{\rho}]\|_1 \le \varepsilon$ in time
  $$\widetilde{O}\!\left(n^{6 + \log\mathfrak{d}/\log(\beta_c/\beta)} \cdot \log^3(1/\varepsilon) \cdot \mathrm{poly}(\mathfrak{d}, \mathfrak{K})\right).$$
  Quantum implementation: depth-1 single-qubit-gate circuit on top of the classical sample.
- **Tightness (no-go)**: the constant-temperature threshold cannot be substantially improved — Sly–Sun antiferromagnetic Ising hardness at $\beta = \Theta(1/\mathfrak{d})$ implies *quantum* Gibbs sampling is also NP-hard past the tree-uniqueness threshold.

## Method

A *sampling-to-counting* recursion ("pinning"): peel sites off $e^{-\beta H}$ one at a time by writing the residual operator as $I + cA$ with $A$ a Pauli, $|c| \le (1 - 3/(5\mathfrak{K}))^{|\mathrm{supp}(A)\cap S|}$ — a quasi-local perturbation of the identity, hence a convex combination of stabilizer product states. The coefficient $c$ is sampled term-by-term from the convergent Taylor series $e^{-\beta H} e^{\beta H_{(j)}} = I + \sum_t (\beta^t/t!) f_t$ via a recurrence $f_{t+1} = -[H, f_t] - f_t H_{(1)}$, with normalisation handled by a fast random walk on the pinning tree (Sinclair–Jerrum-style weak-counting reduction).

## Quantum vs Classical

- **Baseline**: dissipative quantum Gibbs samplers (Chen–Kastoryano–Gilyén KMS Lindbladians, and Rouzé–França–Alhambra's concurrent constant-spectral-gap proof at high T).
- **Gap**: *no super-polynomial gap* in the high-T regime. Both classical (this paper) and quantum (RFA24) achieve $\mathrm{poly}(n)$ runtime; classical even achieves a *depth-1* quantum-circuit preparation.
- **Source of the difference**: structural — at high T the Gibbs state is *separable*, so any quantum advantage relying on entanglement of the target state vanishes. The proof exploits a low-degree polynomial approximation to the imaginary-time propagator $e^{-\beta H} e^{\beta H_{(j)}}$, controlled by cluster-expansion / abstract-polymer convergence at $\beta < 1/(\mathrm{const} \cdot \mathfrak{d}\mathfrak{K})$.
- **Caveat**: the threshold has a constant-factor $\sim 1/100$ overhead vs the conjectured physical critical temperature; matches Sly–Sun hardness up to constants. Result is for low-intersection Hamiltonians with bounded one-local terms — Kuwahara–Hatano counterexample shows large external fields scaling as $\log(1/\beta)/\beta$ can still produce entanglement, so the $|\lambda_a| \le 1$ assumption is essential.

## Implications for Quantum Advantage

- **Regime cell**: noncommuting (any) $\times$ high-T (above $\beta_c$). Closes the high-T row of the 3$\times$3 grid for *all* local Hamiltonians, not just stoquastic.
- **What this changes**: removes high-T quantum Gibbs sampling as a candidate for quantum advantage. Together with Harrow–Mehraban–Soleimanifar (partition-function approximation), Mann–Helmuth (poly-time at high-T), Kuwahara–Kato–Brandão (CMI clustering), and Bakshi–Liu–Moitra–Tang's earlier learning paper, the high-T row of the corpus is now structurally classical: *the state itself carries no quantum correlations*.
- **Promising or not**: **negative** for advantage. Reinforces the corpus framing that the search for quantum-vs-classical separation must move to **intermediate or low temperature**, where entanglement is non-trivial and tensor-network / cluster methods break down. Also implies any downstream task derived from high-T Gibbs states (expectation values, sampling in product bases, learning) is classically efficient.

## Open Questions / Limitations

- Sharpness of the constant in $\beta_c = 1/(c\mathfrak{d}\mathfrak{K})$: gap to the Sly–Sun threshold $\beta = \Theta(1/\mathfrak{d})$ leaves the precise *structural* critical temperature open — is there a temperature window where $\rho$ is separable but classical *sampling* is still hard? [CHECK]
- The exponent $6 + \log\mathfrak{d}/\log(\beta_c/\beta)$ blows up as $\beta \uparrow \beta_c$. Whether this is intrinsic or an artefact of the analysis is unclear.
- The depth-1 implementation needs a classical description first — does a fully *quantum* state-preparation circuit (no classical preprocessing) achieve comparable depth?
- Extension beyond bounded one-local terms (Kuwahara–Hatano regime) and to fermionic / non-abelian-symmetry settings (where infinite-T entanglement is possible) is open.
- Generic locality only: long-range / dense Hamiltonians ($\mathfrak{d}$ scaling with $n$) are not covered.

## Connections

- **Concurrent quantum companion**: Rouzé–França–Alhambra 2024 [RFA24] prove the Chen–Kastoryano–Gilyén KMS Lindbladian has constant spectral gap at high-T — equivalent algorithmic conclusion via a dissipative path; this paper bypasses Lindbladians entirely.
- **Cluster-expansion lineage**: Harrow–Mehraban–Soleimanifar (quasi-poly partition function), Mann–Helmuth 2021 (poly partition function high-T), Kuwahara–Kato–Brandão (CMI clustering above threshold), Yin–Lucas (computational-basis sampling) — same regime, this paper strengthens to *full* state preparation in product form.
- **Hardness baseline**: Sly 2010, Sly–Sun 2014 on antiferromagnetic Ising — give the matching lower bound on the threshold temperature.
- **Authors' learning companion**: Bakshi–Liu–Moitra–Tang (STOC 2024) "Learning quantum Hamiltonians at any temperature" uses related cluster-expansion machinery; learning is easy at all $\beta$, sampling/preparation is easy only at high $\beta$.
- **Contrast — entanglement at infinite T**: non-abelian-symmetric systems [LPRS24, MLWHS24] sit outside this paper's scope, hinting at a structural exception worth flagging in the Review chapter.
