---
paper: "Klassen, Marvian, Piddock, Ioannou, Hen, Terhal (2020)"
title: "Hardness and Ease of Curing the Sign Problem for Two-Local Qubit Hamiltonians"
arxiv: "1906.08800"
year: 2020
venue: "SIAM J. Comput. 49:1332–1362"
pdf: "n/a"

temperature: [low-T, all]
commutativity: [non-commuting]
locality: [2-local, qubit]
particle-statistics: [spin]
hamiltonian-models: [xyz-heisenberg, two-local-qubit, transverse-field-Ising]

paradigm: [classical-QMC, hardness, basis-change]
quantum-or-classical: [classical, comparison, lower-bound]

result-type: [hardness, NP-hard, polynomial-algorithm, structural]
key-scaling: "LocalSignCure is NP-complete with one-local terms (3-SAT reduction); without one-local terms, $O(n^3)$ arithmetic operations over $\\mathbb{R}$ to decide curability."

related: ["hangleiter-roth-nagaj-eisert-2020-easing-sign-problem", "bravyi-divincenzo-oliveira-terhal-2008-stoquastic-lh", "troyer-wiese-2005-sign-problem-hardness", "bravyi-gosset-2017-quantum-ferromagnets"]
---

# Hardness and Ease of Curing the Sign Problem for Two-Local Qubit Hamiltonians — Klassen, Marvian, Piddock, Ioannou, Hen, Terhal (2020)

**One-sentence takeaway**: Establishes a sharp dichotomy for the *single-qubit* basis-change sign-curing problem on two-local qubit Hamiltonians — NP-complete (via 3-SAT) once one-local terms are allowed, but solvable in $O(n^3)$ arithmetic operations when only two-local terms are present (covering pure XYZ Heisenberg in zero field) — pinning the computational boundary of "sign-cure by local rotation" exactly at the addition of a single-qubit field.

## Setting

The sign-curing question asks: given a two-local qubit Hamiltonian $H = \sum_{u<v} H_{uv} + \sum_u h_u$, does there exist an on-site product unitary $U = \bigotimes_i U_i$ such that $UHU^\dagger$ is *stoquastic*, i.e. has real, non-positive off-diagonal entries in the computational basis (a symmetric Z-matrix)?

- **Decision problem (LocalSignCure, Def. 3.1)**: input the local-term coefficients of $H$; output YES if such a $U$ exists, NO otherwise. The transformation set is *single-qubit unitaries* (each $U_i \in U(2)$); equivalently, by an exact gauge argument, real $O(3)$ rotations on the Pauli string of each two-local term.
- **Two regimes by Hamiltonian structure**:
  - **With one-local terms**: arbitrary $\sum_u h_u \cdot \vec\sigma_u$ allowed in addition to two-local couplings.
  - **Without one-local terms ("exactly two-local")**: only $H = \sum_{u<v} H_{uv}$ with $H_{uv} = \sum_{a,b\in\{x,y,z\}} \beta^{(uv)}_{ab}\sigma^a_u\sigma^b_v$.
- **Structural lemma (Prop. 4.3)**: for two-local qubit Hamiltonians, *global* stoquasticity coincides with *termwise* stoquasticity — $H$ is a symmetric Z-matrix iff each $H_{uv}$ is. This locality reduction is what makes the easy side tractable: one can attempt to make every two-local term stoquastic in parallel.

The paper is strictly complexity-theoretic and silent on temperature — the result speaks to *every* world-line / path-integral classical sampler at *any* $\beta$, since stoquasticity is the basis-dependent condition that makes the QMC weights non-negative.

## Main Results

- **Theorem 3.2 (Hardness with one-local terms)**: There exists a family of two-local $n$-qubit Hamiltonians (with one-local terms allowed) for which LocalSignCure is NP-complete. Reduction is from **3-SAT**, not MaxCut. Each clause $(c_1 \vee c_2 \vee c_3)$ is encoded by a clause gadget $H_{\mathrm{OR}} = -(X_d + Z_d + I) \otimes (Z_1 + Z_2 + Z_3 + 2I)$ on an ancilla $d$ and the three clause variables; ancilla penalty terms restrict the allowed single-qubit rotations on each variable to $\{I, W, X, XW\}$ (Hadamard/Pauli group elements), so a sign-curing rotation exists iff the encoded 3-SAT formula is satisfiable (Lemma 5.5).
- **Theorem 3.3 (Polynomial-time algorithm without one-local terms)**: For exactly two-local qubit Hamiltonians (no $h_u$ terms), LocalSignCure is decidable in $O(n^3)$ arithmetic operations over $\mathbb{R}$.
- **XYZ-algorithm (Theorem 2.2, single-edge subroutine)**: For a single two-qubit term $H_{uv}$, decide whether *any* product rotation $U_u \otimes U_v$ makes it stoquastic, by reducing the question to whether the $3\times 3$ Pauli-coefficient matrix $\beta^{(uv)}$ admits an SVD into the form required by a stoquastic XYZ Hamiltonian. Runs in $O(1)$ per edge.
- **Multi-edge consistency (Theorem 3.3 proof)**: edges share qubits, so per-edge XYZ rotations must be globally consistent. The algorithm builds, at each vertex, the intersection of subspaces of admissible $O(3)$ rotations imposed by each incident edge's SVD, and checks non-emptiness via standard linear algebra. The $O(n^3)$ cost is dominated by the SVDs ($O(1)$ each, $O(n^2)$ edges) plus subspace intersections per vertex.
- **Coverage on the easy side**: arbitrary XYZ Heisenberg Hamiltonians $H = \sum_{u<v} (a^{uv}_{xx} X_uX_v + a^{uv}_{yy} Y_uY_v + a^{uv}_{zz} Z_uZ_v)$ with arbitrary couplings (ferromagnetic, antiferromagnetic, frustrated, on any graph). XXZ and XY are special cases. XYZ in zero magnetic field is decidable in $O(n^3)$.

## Method

The hardness side is a 3-SAT-to-LocalSignCure gadget reduction: ancilla qubits are coupled to each variable via penalty terms (Eq. 2 of the paper) that force any sign-curing on-site rotation to lie in a discrete subgroup $\{I, W, X, XW\}$, after which the clause gadget evaluates non-negative iff the corresponding clause is satisfied. The easy side is *not* a Lie-algebraic decomposition — it is a direct linear-algebraic reduction: each two-local Pauli interaction is encoded by a $3 \times 3$ real coefficient matrix $\beta^{(uv)}_{ab}$, the single-edge curability is decided by SVD of $\beta^{(uv)}$ (the XYZ-algorithm), and global consistency across overlapping edges is enforced by intersecting per-vertex orthogonal-rotation subspaces. The crucial structural input is Prop. 4.3 (termwise = global stoquasticity for two-local qubit $H$), which locks the algorithm to per-edge data and stops a global $2^n$-dimensional search.

## Quantum vs Classical

- **Baseline**: Bravyi–DiVincenzo–Oliveira–Terhal [bravyi-divincenzo-oliveira-terhal-2008-stoquastic-lh] defined the stoquastic class; Troyer–Wiese [troyer-wiese-2005-sign-problem-hardness] showed that a *universal* sign-cure cannot exist unless NP $\subseteq$ BPP. The natural intermediate question — can the sign be cured by a single-qubit rotation? — was open at this level of granularity until this paper.
- **Gap**: NP-hard $\Leftrightarrow$ $O(n^3)$, a *super-polynomial* gap (assuming P $\ne$ NP), and the boundary is at *one parameter*: the presence/absence of single-qubit field terms $\sum_u h_u\cdot\vec\sigma_u$.
- **Source of the difference**: with one-local terms, each qubit has independent rotation-space constraints from the field plus additional constraints from every incident two-local edge — these constraints are no longer linearly compatible across the lattice, and one can encode arbitrary 3-SAT into their satisfaction. Without one-local terms, every constraint is between *pairs* of qubits, so the per-edge SVD admits a finite-dimensional consistency check at each vertex — the constraint geometry becomes tractable.
- **Caveat**: the result is for *single-qubit* (on-site) unitaries only. Quasi-local Clifford or matrix-product-unitary basis changes are not addressed and remain open. Likewise, the $O(n^3)$ algorithm only *decides* curability — when the answer is YES it returns the rotations; when it is NO it does not return a "near-stoquastic" approximation (that is the question taken up by Hangleiter et al. [hangleiter-roth-nagaj-eisert-2020-easing-sign-problem]).

## Implications for Quantum Advantage

- **Regime cell**: non-commuting × all-T × 2-local × spin (the canonical non-stoquastic row of the corpus 3×3 grid, with the temperature axis irrelevant since the sign problem is structural).
- **What this changes**:
  - **Refines** the non-stoquastic row's classical-QMC barrier: the boundary between "easy" and "hard" sign-curing-by-local-rotation lies *exactly* at the on/off switch for single-qubit fields. Pure XYZ Heisenberg without external field is classically curable in $O(n^3)$; adding any nonzero $h_u$ flips the problem to NP-hard.
  - **Closes** a tempting classical workaround: it is no longer reasonable to hope that a polynomial-time preprocessor inspects a generic two-local $H$ with field and returns a stoquastic rotation when one exists. (Though Hangleiter et al. show the situation worsens further: even *easing* — relaxing curability to "minimise the non-stoquastic norm" — remains NP-complete on strictly two-local XYZ without one-local terms, under the on-site Clifford or general orthogonal ansatz.)
- **Promising or not**: this paper is a *barrier paper* for classical QMC in the non-stoquastic regime. It does not produce a quantum-advantage construction directly, but it sharpens the structural reason why quantum Gibbs samplers (Chen et al. KMS Lindbladians, Ding et al., Ramkumar–Soleimanifar) become candidate winners on two-local qubit Hamiltonians with magnetic fields — they sample $\rho_\beta = e^{-\beta H}/Z$ in the energy eigenbasis without any basis rotation, sidestepping the entire LocalSignCure obstruction.
- **Headline for the thesis Review chapter**: the natural physical setting — frustrated antiferromagnets in a transverse field, transverse-field Ising, Heisenberg with Zeeman coupling — is *exactly* the field-on side where this paper rules out efficient sign-curing-by-local-rotation. Pure exchange-only XYZ in zero field, where sign-curing is easy, is the structurally exceptional case.

## Open Questions / Limitations

- **Single-qubit ansatz only**: more general bases (quasi-local Clifford, on-site $O(d^k)$ for $k>1$, matrix product unitaries) are not analyzed; their complexity remains open. Hangleiter et al. partially address this by showing easing is hard even under on-site Clifford on strictly 2-local XYZ.
- **No approximation result**: when LocalSignCure returns NO, the algorithm is silent on the approximate ($\nu_1$-easing) version. The Hangleiter et al. follow-up makes this explicit and shows even the approximate variant is NP-complete in a complementary regime.
- **Qubits only**: the dichotomy is for two-level systems; the qudit ($d \ge 3$) case is not treated and structurally different (Prop. 4.3 may fail).
- **Restricted to two-local**: $k$-local for $k \ge 3$ is open. The hardness of LocalSignCure could move further into the easy regime or stay NP-hard; no general statement is available.
- **No quantitative connection to QMC variance**: the paper is binary curable/not-curable. The connection to actual world-line sample complexity (average sign decay rate) goes through Hangleiter et al.'s $\nu_1$ measure, not LocalSignCure directly.
- **[CHECK]** Theorem 3.3's exponent: the paper states "$O(n^3)$ arithmetic operations over $\mathbb{R}$"; verify whether this is tight or merely an upper bound from the proof structure. Also [CHECK] whether bit complexity (rather than arithmetic complexity) introduces extra polynomial factors.

## Connections

- **[bravyi-divincenzo-oliveira-terhal-2008-stoquastic-lh]**: defines the stoquastic class. The current paper makes the basis-change variant of the BDOT question quantitative at the local-rotation granularity, sharpening "is $H$ stoquastic?" to "is $H$ on-site-rotatable to stoquastic?" with a sharp dichotomy.
- **[troyer-wiese-2005-sign-problem-hardness]**: the worst-case "sign-cure NP-hard" result via #P-hard counting. This paper is the natural fine-grained companion, replacing "any cure" with "single-qubit rotation cure" and giving an explicit 3-SAT reduction with explicit polynomial-time solvability on the easy side.
- **[hangleiter-roth-nagaj-eisert-2020-easing-sign-problem]**: the *quantitative* / *easing* analogue. Together, the two papers form the canonical pair for the non-stoquastic row: Klassen et al. = curability dichotomy under single-qubit rotations; Hangleiter et al. = NP-hardness of approximate easing under on-site Clifford / orthogonal — even on the strictly 2-local subclass that is "easy" *exactly* in Klassen et al.'s sense. Note the apparent tension is resolved by the difference in question: LocalSignCure asks $\nu_1 = 0$ feasibility, SignEasing asks $\nu_1 \le A$ optimization, and the latter is harder even where the former is easy.
- **[bravyi-gosset-2017-quantum-ferromagnets]**: a concrete instance on the easy side — stoquastic transverse-field XY ferromagnets admit a classical FPRAS for the partition function. The current paper rules out efficient sign-curing once those Hamiltonians acquire a generic single-qubit field, showing how delicate this easy regime is.
- **Marvian–Lidar–Hen (2019)** and **Ioannou–Terhal (arXiv:2007.11964) on termwise vs global stoquasticity**: extend the picture to extended-Clifford basis changes and to the termwise-versus-globally stoquastic distinction; both build directly on this paper's gadget framework.
- **Thesis Review chapter (non-stoquastic row)**: alongside Hangleiter et al. and Troyer–Wiese, this is a Tier-1 barrier citation for "why classical sign-curing won't save QMC on two-local qubit Hamiltonians with fields" — the canonical motivating obstruction for direct quantum Gibbs sampling on non-stoquastic $H$.
