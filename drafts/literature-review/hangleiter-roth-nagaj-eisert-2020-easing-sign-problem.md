---
paper: "Hangleiter, Roth, Nagaj, Eisert (2020)"
title: "Easing the Monte Carlo sign problem"
arxiv: "1906.02309"
year: 2020
venue: "Sci. Adv. 6, eabb8341"
pdf: "supplementary-informations/classical-review/hangleiter-roth-nagaj-eisert-2020-easing-sign-problem.pdf"

temperature: [n/a]
commutativity: [non-stoquastic, sign-problem-bearing]
locality: [2-local, k-local, translation-invariant]
particle-statistics: [spin, fermionic]
hamiltonian-models: [Heisenberg-ladder, J0-J1-J2-J3, frustrated-Heisenberg, antiferromagnetic-Ising, XYZ]

paradigm: [classical-other, structural, lower-bound]
quantum-or-classical: [classical, comparison]

result-type: [no-go, complexity-lower, structural, heuristic]
key-scaling: "SignEasing for 2-local (XYZ) Hamiltonians is NP-complete under both on-site orthogonal Clifford and on-site general orthogonal transformations (Theorem 2); $\\nu_1$ efficiently computable in $O(1)$ for translation-invariant local $H$."

related: ["bravyi-divincenzo-oliveira-terhal-2008-stoquastic-lh", "klassen-marvian-2020-curing-sign-problem", "troyer-wiese-2005-sign-problem-hardness", "bravyi-terhal-2010-frustration-free-stoquastic", "bravyi-gosset-2017-quantum-ferromagnets"]
---

# Easing the Monte Carlo Sign Problem — Hangleiter, Roth, Nagaj, Eisert (2020)

**One-sentence takeaway**: Replaces the binary "stoquastic vs. not" question with a quantitative, efficiently computable measure of non-stoquasticity $\nu_1(H) = D^{-1}\|H_\neg\|_{\ell_1}$ that bounds the variance of the QMC estimator, gives a heuristic geometric optimizer over local orthogonal bases that practically eases the sign problem on frustrated Heisenberg ladders by factors $2$–$5$, and proves that *optimally* easing the sign problem under on-site orthogonal Clifford or general orthogonal transformations is NP-complete already for 2-local Hamiltonians via reduction from MAXCUT — sharpening the worst-case classical-QMC barrier in the non-stoquastic regime.

## Setting

The world-line QMC sign problem is *basis dependent*: positivity of off-diagonal entries of $-H$ in some product basis suffices to remove it. Bravyi–DiVincenzo–Oliveira–Terhal call such $H$ stoquastic [bravyi-divincenzo-oliveira-terhal-2008-stoquastic-lh]; Klassen et al. show that *deciding* stoquasticity-by-single-qubit-rotations is NP-hard [klassen-marvian-2020-curing-sign-problem]; Troyer–Wiese show that universally curing the sign problem is NP-hard [troyer-wiese-2005-sign-problem-hardness]. This paper asks the more pragmatic question — what if we settle for *easing* rather than *curing*? Concretely:

- **Measure of non-stoquasticity**: for a real $D \times D$ Hamiltonian, define
  $$ \nu_1(H) := D^{-1}\|H_\neg\|_{\ell_1}, $$
  where $H_\neg$ keeps the entries $h_{i,j}$ for $h_{i,j}>0,\ i\ne j$ and zeros otherwise, and $\|\cdot\|_{\ell_1}$ is the entrywise (vector) $\ell_1$-norm. For local Hamiltonians on bounded-degree graphs $\nu_1$ is computable from the local terms; for translation-invariant $H$ it is computable in $O(1)$ effort.
- **Ansatz class**: on-site orthogonal transformations $\mathcal O = O^{\otimes n}$ with $O \in O(d)$, which preserve locality and translation invariance, so $\nu_1(\mathcal O H \mathcal O^\top)$ depends only on the transformed local term.
- **Connection to QMC sample complexity**: if $T_m = \mathbb{1} - \beta H/m$ is the world-line transfer matrix, the relative variance of the partition-function estimator equals $\langle\mathrm{sign}\rangle^{-2}_p - 1$, and the average sign $\langle\mathrm{sign}\rangle = \mathrm{Tr}[T_m^m]/\mathrm{Tr}[|T_m|^m]$ governs sampling cost. The authors argue analytically and numerically that for *generic* 2-local $H$ the average sign decays as $\exp(-c\,D\,\nu_1(H))$, so minimising $\nu_1$ minimises the sample complexity (Sec. II of the paper).

The paper is silent on temperature — its claims are structural / complexity-theoretic and apply at any $\beta$ where world-line QMC is the natural classical sampler.

## Main Results

- **Theorem 2 (NP-completeness of SignEasing)**: SignEasing is NP-complete for 2-local (XYZ) Hamiltonians under (i) on-site orthogonal *Clifford* transformations and (ii) on-site general orthogonal transformations. Verbatim: *"SignEasing is NP-complete for 2-local (XYZ) Hamiltonians under (i) on-site orthogonal Clifford transformations, and (ii) on-site general orthogonal transformations."* Robust under any $\ell_p$-norm choice for the non-stoquastic part. Holds even on subgraphs of the double-layered square lattice (degree six), so hard instances appear at low spatial dimension and constant interaction strength.
- **Decision problem (Definition 1)**: given $H$, $B>A\ge 0$ with $B-A \ge 1/\mathrm{poly}(n)$, and an allowed transformation set $\mathcal U$, decide YES: $\exists U \in \mathcal U: \nu_1(UHU^\dagger) \le A$ vs NO: $\forall U \in \mathcal U: \nu_1(UHU^\dagger) \ge B$.
- **Containment in NP**: SignEasing for arbitrary 2-local $H$ is in NP — given a basis transformation, $\nu_1$ can be approximated to inverse polynomial error from local terms (their Theorem 6).
- **Reduction (sketch)**: encode an antiferromagnetic-Ising MAXCUT instance $H_{\text{AFI}} = \sum_{(i,j)\in E} Z_iZ_j$ into a Hamiltonian $H'$ whose every $ZZ$ edge becomes an $XX$ edge plus penalty terms $C(Z_iZ_j - Z_iZ_{\xi_{i,j}} - Z_jZ_{\xi_{i,j}})$ with ancilla $\xi_{i,j}$ and $C = 2\deg(G)$. Penalty terms force the optimum easing to be of the form $Z^{s_1}\cdots Z^{s_n}$ with $s\in\{0,1\}^n$, so $\nu_1(H')$ falls below the threshold iff the original MAXCUT instance has ground-state energy below the corresponding value. Containment-cum-hardness completes NP-completeness.
- **Indicator independence**: the hardness conclusion is independent of which $\ell_p$-norm of the non-stoquastic part one minimises — the sign-problem barrier is a structural property of the basis-rotation problem, not an artifact of the chosen measure.
- **Heuristic algorithm (proof of principle)**: a conjugate-gradient method on the orthogonal-group manifold $O(d)$, initialised at identity. On random curable Hamiltonians it accurately recovers a stoquastic basis (relative non-stoquasticity dropping by orders of magnitude at $d=2$); on the $J_0$–$J_1$–$J_2$–$J_3$ frustrated Heisenberg model and the frustrated Heisenberg ladder it improves $\nu_1$ by factors $2$–$5$ across the phase diagram, and the inverse average sign drops by factors up to $\sim 10^4$ in the parameter regime studied. Limitations: rugged optimization landscape, frustrated couplings $J_\times \ne J_\parallel$ remain non-trivial, no convergence certificate.
- **Two pathological examples (Sec. II A)**: explicit construction of a sign-problem-free Hamiltonian with $\nu_1 = n$ (Example 3), and of a Hamiltonian with $\nu_1 = bm/(2\beta)$ but average sign $\le C(b-a)/a$ which can be made arbitrarily small (Example 4). Conclusion: a continuous *one-to-one* map $\nu_1 \mapsto \langle\mathrm{sign}\rangle$ does not exist; only a probabilistic / generic correspondence does.

## Method

The methodology has three pillars: (1) a tractable surrogate measure $\nu_1$ that lower-bounds the inverse average sign generically yet is efficiently computable from local terms (avoiding the chicken-and-egg situation that the average sign itself requires QMC sampling, which suffers the very sign problem one is trying to ease); (2) a Riemannian conjugate-gradient optimization over the orthogonal manifold $O(d)$ acting on the local single-site basis, exploiting translation invariance so the optimization scales as $O(1)$ in $n$; (3) a complexity reduction encoding antiferromagnetic-Ising / MAXCUT into the basis-rotation problem via $Z\to X$ swaps with ancilla penalty terms that force the optimum to be on-site Clifford. Proofs of NP-hardness combine the gadget construction in Fig. 3 with a careful enumeration of orthogonal transformations preserved by the penalty.

## Quantum vs Classical

- **Where this lives in the corpus grid**: the *non-stoquastic* row, all temperature columns, but most relevantly intermediate-T and low-T (the stoquastic case is handled by direct QMC; the non-stoquastic case is where the sign problem and hence basis-easing matters).
- **Tightens the worst-case classical-QMC barrier**: Troyer–Wiese (2005) showed *some* universal sign-problem cure cannot exist unless NP $\subseteq$ BPP. Klassen et al. and Marvian–Lidar–Hen showed *single-qubit basis curing* is NP-hard already. This paper extends those barriers to *easing* — even relaxing the goal from "make $\nu_1=0$" to "minimise $\nu_1$" remains NP-hard for 2-local $H$ under on-site Clifford or orthogonal transformations. So no efficient classical preprocessor can be expected to render arbitrary frustrated 2-local QMC tractable in the worst case, even approximately.
- **Where classical QMC retains traction**: the analytical and numerical work shows that *generically* (i.e., on translation-invariant frustrated Heisenberg ladders, away from worst-case constructions) heuristic easing reduces sample complexity by orders of magnitude. The sign problem is structurally hard but practically easable on physically-motivated subclasses — this is the current best classical answer to "what can be salvaged in the non-stoquastic regime?".
- **Quantum-side complement**: the worst-case hardness of basis-easing is exactly the structural reason why quantum Gibbs samplers (Chen et al. KMS Lindbladians, Ding et al., Ramkumar–Soleimanifar) become candidate winners on non-stoquastic $H$ — they sample $\rho_\beta = e^{-\beta H}/Z$ directly without needing a stoquastic basis, side-stepping the entire $\nu_1$-minimisation tower. The thesis Review chapter punchline cell (non-stoquastic $\times$ low-T) is exactly the corner where this paper closes off the most plausible classical workaround.
- **What this paper does *not* show**: it does not exhibit a non-stoquastic Hamiltonian that is provably classically hard *and* quantum-easy. It rules out polynomial-time *worst-case* classical easing under restricted ansätze; it does not rule out polynomial-time easing under more general (e.g. quasi-local Clifford or matrix-product-unitary) ansätze, nor does it imply that every non-stoquastic instance is genuinely hard — many structured instances do admit known cures (auxiliary-field / Jordan–Wigner / dimerisation tricks).

## Implications for Quantum Advantage

- **Refines the boundary of "classical QMC reach"**: the regime grid's non-stoquastic row is hardened — naive classical QMC fails by Troyer–Wiese; sign-curing fails by Klassen et al.; sign-easing fails by Theorem 2 here. Only specialised structural tricks (worm, contour deformation, dual variables, Lefschetz thimbles) or heuristic optimisers without worst-case guarantees survive.
- **Identifies a structural source of quantum advantage**: a quantum Gibbs sampler that operates in the energy / Heisenberg eigenbasis (e.g. the KMS Lindbladians of Chen et al. 2023, 2025; Ding et al.) needs no basis rotation at all, hence is immune to the $\nu_1$-easing barrier. Quantum advantage in non-stoquastic Gibbs sampling is therefore not a "constant-factor speedup" claim but a *worst-case complexity separation* candidate, conditional on standard hierarchy assumptions.
- **Quasi-local circuit caveat (their own outlook)**: if one allows the basis change to be a quasi-local circuit (matrix product unitary, short circuit) that is still efficiently computable, the easing problem may move out of NP-hardness — though the authors point out topological obstructions for some Hamiltonians (Levin–Wen, fixed-point Hamiltonians of non-chiral topologically ordered phases) that may carry an "intrinsic" sign problem. This is the theoretically interesting boundary of what classical methods could in principle still achieve, and it remains open.
- **Promising or not**: this paper is a *barrier* for classical methods; it does not directly produce a quantum-advantage construction. Its value to the thesis is that it rules out the most natural classical workaround (basis-easing-by-polynomial-preprocessing) for the non-stoquastic regime, leaving direct quantum Gibbs sampling as the only remaining candidate for honest poly-time guarantees there.

## Open Questions / Limitations

- The hardness is for *worst-case* 2-local instances; what fraction of physically motivated non-stoquastic Hamiltonians (frustrated antiferromagnets, Hubbard at finite doping, lattice fermions) are easy or hard for the heuristic optimizer is empirical. The Heisenberg-ladder benchmark is suggestive but not a complexity result.
- The link between $\nu_1$ and $\langle\mathrm{sign}\rangle$ is "generically" exponential in $D\nu_1$ but Examples 3 and 4 show it can fail in either direction. A rigorous *typical-case* statement (e.g. for random local Hamiltonians under a natural distribution) is missing.
- All hardness proofs use *real* Hamiltonians and orthogonal transformations; the authors flag the natural extension to complex Hamiltonians (the *complex phase problem*) as an open direction.
- Translation invariance is used to make $\nu_1$ computable in $O(1)$; for disordered or non-translation-invariant systems the optimization landscape grows and even computing $\nu_1$ becomes more delicate.
- Whether quasi-local-Clifford or matrix-product-unitary ansätze admit polynomial-time easing for relevant subclasses is open; the paper conjectures rich structure here but proves nothing.
- Method ≤3 sentence cap: the heuristic optimizer is not certified — convergence to local minima of the rugged $\nu_1$ landscape is observed in the frustrated ladder data and is the empirical bottleneck.
- [CHECK] The factor "$2$–$5$" improvement and "$\sim 10^4$" inverse-average-sign drop are read directly off Figs. 1 and 2 of the paper at $\beta=1$, $m=100$, $2\times 4$ ladder; these are proof-of-principle scales, not asymptotics.

## Connections

- **[bravyi-divincenzo-oliveira-terhal-2008-stoquastic-lh]**: defines the stoquastic class, the binary version of the question this paper makes quantitative. The current paper's $\nu_1$ is the natural metric on the BDOT structural definition, and Theorem 2 here lifts the BDOT "stoquastic-or-not" to "how-far-from-stoquastic and is the optimum efficiently computable?".
- **[klassen-marvian-2020-curing-sign-problem]**: NP-hardness of *curing* (i.e., $\nu_1=0$) under single-qubit unitaries for general 2-local Hamiltonians with one-local terms; polynomial-time *cure* for strictly 2-local. This paper extends to *easing* (i.e., $\nu_1\le A$) and removes the "with one-local terms" caveat — easing is hard already for strictly 2-local (XYZ) Hamiltonians.
- **[troyer-wiese-2005-sign-problem-hardness]**: the original worst-case sign-problem hardness via reduction from #P-hard counting. The present paper's Theorem 2 specialises Troyer–Wiese to a concrete *basis-rotation* class with the matching MAXCUT reduction — explicit, constructive, and tight at 2-local.
- **Marvian–Lidar–Hen (2019)**: extended Clifford-class basis curing of 2-local $H$ is NP-hard. This paper's Theorem 2(i) is the easing analogue under on-site Clifford; the two together establish that no Clifford-class polynomial preprocessor can render frustrated 2-local QMC tractable in the worst case.
- **[bravyi-terhal-2010-frustration-free-stoquastic]** and **[bravyi-gosset-2017-quantum-ferromagnets]**: the easy stoquastic side — frustration-free or ferromagnetic stoquastic $H$ admit polynomial classical Gibbs / partition-function algorithms. Together with the present paper they sharpen the dichotomy: stoquastic + structure $=$ classically tractable; non-stoquastic generic 2-local $=$ classically NP-hard to ease, hence prime quantum-advantage frontier.
- **Thesis Review chapter (corpus row "non-stoquastic")**: this paper is the canonical citation for *why* the non-stoquastic row is the regime where quantum Gibbs sampling could first beat classical — it closes off the obvious classical workaround.
