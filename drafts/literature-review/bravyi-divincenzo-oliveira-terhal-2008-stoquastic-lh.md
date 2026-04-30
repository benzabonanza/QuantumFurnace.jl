---
paper: "Bravyi, DiVincenzo, Oliveira, Terhal (2008)"
title: "The Complexity of Stoquastic Local Hamiltonian Problems"
arxiv: "quant-ph/0606140"
year: 2008
venue: "Quantum Inf. Comput. 8, 361–385"
pdf: "supplementary-informations/classical-review/bravyi-divincenzo-oliveira-terhal-2008-stoquastic-LH.pdf"

temperature: [n/a]
commutativity: [stoquastic, sign-problem-free]
locality: [k-local, 2-local, 6-local]
particle-statistics: [spin, boson]
hamiltonian-models: [stoquastic-LH, transverse-field-Ising, Heisenberg-ferromagnet, bosonic-Hubbard, Jaynes-Cummings, spin-boson]

paradigm: [structural, complexity-theoretic]
quantum-or-classical: [comparison, classical-side-baseline]

result-type: [complexity-upper, complexity-lower, structural]
key-scaling: "Stoquastic LH-MIN $\\in$ AM; 6-local stoquastic LH-MIN is MA-hard; 2-local stoquastic LH-MIN is MA-hard (via gadgets); with promise of $1/\\mathrm{poly}(n)$ spectral gap, stoquastic LH-MIN $\\in$ PostBPP."

related: ["bravyi-terhal-2010-frustration-free-stoquastic", "bravyi-gosset-2017-quantum-ferromagnets", "klassen-marvian-2020-curing-sign-problem", "hangleiter-roth-nagaj-eisert-2020-easing-sign-problem", "troyer-wiese-2005-sign-problem-hardness"]
---

# The Complexity of Stoquastic Local Hamiltonian Problems — Bravyi, DiVincenzo, Oliveira, Terhal (2008)

**One-sentence takeaway**: Defines the *stoquastic* class — Hamiltonians whose off-diagonal matrix elements in the standard basis are real and non-positive, equivalently those with entry-wise non-negative Gibbs density at every $\beta \ge 0$ — and proves that their local-Hamiltonian problem sits strictly between MA and AM, well below the QMA-completeness of generic local Hamiltonians, formalising the "no sign problem" boundary that determines when classical Monte Carlo can simulate a quantum system.

## Setting

The Local Hamiltonian Problem (LH-MIN): given a $k$-local $n$-qubit Hamiltonian $H = \sum_S H_S$ with $\|H_S\| \le \mathrm{poly}(n)$ and the promise $\lambda(H) \le 0$ or $\lambda(H) \ge \delta = 1/\mathrm{poly}(n)$, decide which. A Hamiltonian is **stoquastic** in a given product basis (typically computational $z$-basis) iff $\langle x | H | y \rangle \le 0$ for all $x \ne y$. Equivalently $\rho_\beta = e^{-\beta H} / Z$ has entry-wise non-negative matrix elements for every $\beta \ge 0$, so the ground state $|\Psi_0\rangle = \sum_i \alpha_i |i\rangle$ has $\alpha_i \ge 0$, defining a probability distribution $\mathbb{P}(i) = \alpha_i / \sum_a \alpha_a$. Examples: ferromagnetic Heisenberg, transverse-field Ising, antiferromagnetic Heisenberg on bipartite graphs (after a sublattice basis change), bosonic Hubbard models, Jaynes–Cummings, spin–boson, and all flux-qubit Josephson-junction Hamiltonians. *No temperature regime* — this is a complexity-theoretic structural result, but it underlies every later "classical QMC works iff stoquastic" claim.

## Main Results

- **Theorem (LH-MIN $\in$ AM)**: stoquastic LH-MIN belongs to AM (Arthur–Merlin, two-round randomized interactive proof). Verbatim: *"stoquastic LH-MIN belongs to the class AM."* Mechanism: the Stoquastic Largest Eigenvalue Problem for a non-negative matrix $G$ reduces to estimating $\mathrm{Tr}(G^L)$, which is a non-negative sum and hence an Approximate Set Size problem solvable by the Goldwasser–Sipser hashing protocol.
- **Theorem (LH-MIN is MA-hard)**: 6-local stoquastic LH-MIN is hard for MA (the probabilistic NP analogue with classical witness). Verbatim: *"6-local stoquastic LH-MIN is MA-hard."* Mechanism: replace Arthur's BPP verifier with a *classical reversible* circuit (Toffoli gates) acting on $|+\rangle$ "coin" qubits and Merlin's witness, then apply Kitaev's clock construction; reversibility makes the resulting clock Hamiltonian stoquastic.
- **Theorem (k-local $\to$ 2-local)**: Any constant $k$-local termwise-stoquastic LH-MIN reduces in polynomial time to 2-local stoquastic LH-MIN. Corollary: 2-local stoquastic LH-MIN is MA-hard. Mechanism: a new three-qubit perturbation-theory gadget that preserves stoquasticity (off-diagonal terms remain non-positive) — note that *k* not bounding completeness suggests stoquastic LH-MIN may be complete for an as-yet-unidentified class strictly between MA and AM.
- **Theorem (gapped $\to$ PostBPP)**: With the additional promise of a $1/\mathrm{poly}(n)$ spectral gap, stoquastic LH-MIN is in PostBPP $=$ BPP$_{\mathrm{path}}$ (BPP with classical post-selection). Corollary: any decision problem solvable by adiabatic quantum computation along a stoquastic path with poly-gap is in PostBPP — i.e., classical post-selected sampling matches stoquastic adiabatic QC. Mechanism: Green's-function QMC + post-selection on a "success" flag.

## Method

The AM upper bound expresses $\mathrm{Tr}(G^L)$ for a non-negative matrix $G = (I - H/C)/2$ as a counting problem (cardinality of a Boolean set) and invokes Goldwasser–Sipser hashing. The MA-hardness uses a Kitaev clock Hamiltonian whose verifier circuit is reversible-classical, making all $H_{\mathrm{prop}}$ off-diagonal elements non-positive. The 2-local reduction is a standard subdivision-plus-mediator perturbation gadget engineered so that the Schrieffer–Wolff effective Hamiltonian remains stoquastic at second order.

## Quantum vs Classical

- **Baseline**: Generic local Hamiltonian problem is QMA-complete (Kitaev), even for 2-local nearest-neighbour qubits on a 2D lattice (Oliveira–Terhal, Kempe–Kitaev–Regev). No classical Monte Carlo can sample $\rho_\beta$ in polynomial time for a generic non-stoquastic $H$ unless QMA collapses.
- **Gap**: For *stoquastic* $H$, the gap to classical hardness narrows dramatically: LH-MIN $\in$ AM $\subseteq \Sigma_2^p$, so unless QMA $\subseteq$ AM (considered very unlikely), generic non-stoquastic LH-MIN is strictly harder than stoquastic LH-MIN. With a poly gap promise, stoquastic LH-MIN lies in PostBPP — almost-classical.
- **Source of the difference**: Non-negativity of $e^{-\beta H}$ in the standard basis means path-integral world-line weights are non-negative, removing the sign problem. Equivalently, the ground-state wavefunction has non-negative amplitudes, so $|\Psi_0|^2$-importance sampling is well-defined classically.
- **Caveat**: "Sign-problem-free in some basis" is a *basis-dependent* property; deciding whether a 2-local Hamiltonian admits a single-qubit basis change to stoquastic form is itself NP-hard ([klassen-marvian-2020-curing-sign-problem]). The class is also not closed under tensor products with general non-stoquastic terms.

## Implications for Quantum Advantage

- **Regime cell**: Defines the *stoquastic row* of the corpus 3×3 grid. Determines the location of the *quantum-vs-classical frontier*: if a Hamiltonian is stoquastic, classical QMC has at least a chance; if non-stoquastic, the Troyer–Wiese sign-problem hardness applies, and quantum Gibbs samplers (KMS Lindbladians, Chen et al., Ding et al.) become the natural candidates.
- **What this changes**: Forces any quantum-advantage claim for Gibbs sampling to either (i) target non-stoquastic Hamiltonians, or (ii) target stoquastic Hamiltonians without the polynomial-spectral-gap promise (since gapped stoquastic adiabatic computation is in PostBPP — a classical-with-post-selection class). Spectrally-gapped stoquastic frustration-free Hamiltonians are *provably* classically polynomial-time simulable ([bravyi-terhal-2010-frustration-free-stoquastic]), removing them from the advantage frontier entirely.
- **Promising or not**: **Not promising for advantage** as a target class — stoquastic Hamiltonians are the canonical "easy" quantum case classically. The paper's value to the thesis is the opposite: it sharply *delineates* the easy-classical regime, leaving non-stoquastic Hamiltonians (especially at low $T$ / past phase transitions) as the residual cell where quantum Gibbs samplers can plausibly outperform classical methods.

## Open Questions / Limitations

- The exact complexity class of stoquastic LH-MIN remains open: known to be MA-hard and in AM, but neither MA-complete (likely not, given the $k$-locality independence) nor known to be AM-hard. A natural candidate class strictly between MA and AM is suggested but unidentified [CHECK: a follow-up by Bravyi (cited as [22] in the paper) reportedly strengthens to SBP].
- The PostBPP result requires *classical post-selection*, which is not BPP and is widely believed strictly stronger; it does not literally imply stoquastic adiabatic QC offers no quantum advantage over BPP, only over PostBPP.
- The result is silent on *finite-temperature Gibbs sampling* per se: it concerns ground-state-energy decisions. The corresponding quantitative classical sampling guarantees come later (Bravyi–Terhal 2010, Bravyi 2015 PI-MCMC, Mann–Helmuth polymer methods) and require additional structure (frustration-freeness, stability, high $T$).
- Stoquasticity is basis-dependent and not preserved under arbitrary local unitaries; the structural significance therefore depends on whether the natural physical basis is stoquastic.

## Connections

- **[bravyi-terhal-2010-frustration-free-stoquastic]**: Strengthens the "easy" side: frustration-free stoquastic Hamiltonians with $1/\mathrm{poly}$ gap admit classical RP-time simulation, sharper than the PostBPP bound here.
- **[bravyi-gosset-2017-quantum-ferromagnets]**: Concrete instance — provides a polynomial-time *partition-function* FPRAS for stoquastic transverse-field XY ferromagnets, demonstrating that for specific stoquastic subclasses classical Gibbs sampling is fully tractable.
- **[klassen-marvian-2020-curing-sign-problem]** and **[hangleiter-roth-nagaj-eisert-2020-easing-sign-problem]**: Address the *boundary* — when can a non-stoquastic Hamiltonian be transformed to stoquastic form by a local basis change? Klassen et al.: NP-hard in general, polynomial-time only in restricted no-one-local-term subclass.
- **[troyer-wiese-2005-sign-problem-hardness]**: The complementary worst-case hardness — a generic polynomial cure of the sign problem would put NP $\subseteq$ BPP. Together with the present paper this fixes the stoquastic-vs-non-stoquastic dichotomy: stoquastic = polynomially-controlled by classical QMC up to spectral-gap issues; non-stoquastic = generically hard for any classical sampler.
- *Thesis Review chapter:* this paper supplies the structural definition that organises the entire stoquastic row of the 3×3 grid; every "classical Gibbs sampler works for $H$" claim in the corpus is implicitly conditional on a stoquastic representation.
