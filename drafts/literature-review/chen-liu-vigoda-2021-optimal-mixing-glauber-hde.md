---
paper: "Chen, Liu, Vigoda (2021)"
title: "Optimal Mixing of Glauber Dynamics: Entropy Factorization via High-Dimensional Expansion"
arxiv: "2011.02075"
year: 2023
venue: "STOC 2021 / SIAM J. Comput. 2023"
pdf: "supplementary-informations/classical-review/chen-liu-vigoda-2021-optimal-mixing-glauber-hde.pdf"

temperature: [high-T, intermediate-T]
commutativity: [commuting]
locality: [sparse, geometric-local]
particle-statistics: [spin]
hamiltonian-models: [Ising, hardcore, antiferromagnetic-2-spin, ferromagnetic-2-spin, q-colorings, monomer-dimer]

paradigm: [classical-MCMC]
quantum-or-classical: [classical]

result-type: [mixing-time-upper, mlsi, structural]
key-scaling: "t_mix(P_GL, eps) = (Delta/b)^{O(eta/b^2 + 1)} * O(n log(n/eps)); MLSI constant rho_LS >= 1/(C_1 n) with C_1 = (Delta/b)^{O(eta/b^2 + 1)}. For hardcore at lambda <= (1-delta) lambda_c(Delta) and Ising at beta in [(Delta-2+delta)/(Delta-delta), (Delta-delta)/(Delta-2+delta)], mixing time C(Delta, delta) n log n with C ~ Delta^{O(Delta^2/delta)}. For triangle-free q-colorings with q > (alpha* + delta) Delta (alpha* approx 1.763), C n log n. For monomer-dimer, C(Delta, lambda) m log n."

related: ["anari-liu-oveisgharan-2020-spectral-independence-hardcore", "chen-eldan-2022-localization-schemes", "bauerschmidt-dagallier-2024-near-critical-ising-lsi", "weitz-2006-counting-independent-sets-tree-threshold", "sly-2010-computational-transition-uniqueness"]
---

# Optimal Mixing of Glauber Dynamics: Entropy Factorization via High-Dimensional Expansion — Chen, Liu, Vigoda 2021

**One-sentence takeaway**: Promotes the Anari–Liu–Oveis Gharan spectral-independence framework from a spectral-gap bound to a *modified log-Sobolev inequality* via approximate block-factorization of entropy on pure weighted simplicial complexes, yielding the first **optimal $O(n\log n)$** Glauber mixing for hardcore and antiferromagnetic 2-spin systems on any bounded-degree graph throughout the entire tree-uniqueness region — matching Sly's hardness boundary up to constants.

## Setting

- **Targets**: classical $q$-spin systems on a graph $G=(V,E)$, $|V|=n$, defined by a symmetric interaction matrix $A\in\mathbb{R}_{\ge 0}^{q\times q}$ and external field $h\in\mathbb{R}_{>0}^q$, with Gibbs distribution $\mu(\sigma) \propto \prod_{\{u,v\}\in E} A(\sigma_u,\sigma_v)\prod_{v\in V} h(\sigma_v)$. Specializes to: hardcore (independent sets at fugacity $\lambda$), Ising (params $\beta,\gamma,\lambda$ with ferro $\beta\gamma>1$ / antiferro $\beta\gamma<1$), proper $q$-colorings, monomer-dimer matchings.
- **Hamiltonian class**: classical, commuting; hardcore and antiferromagnetic Ising are 2-spin antiferromagnetic, the canonical hard cell on the classical side.
- **Temperature regime**: tree-uniqueness regime — $\lambda<\lambda_c(\Delta)=(\Delta-1)^{\Delta-1}/(\Delta-2)^\Delta$ for hardcore; $\beta\in[(\Delta-2+\delta)/(\Delta-\delta),(\Delta-\delta)/(\Delta-2+\delta)]$ for Ising; $q>(\alpha^*+\delta)\Delta$ with $\alpha^*\approx 1.763$ for triangle-free colourings.
- **Access**: discrete-time single-site Glauber (heat-bath) dynamics $P_{\mathrm{GL}}$; the "Hayes–Sinclair" lower bound $\Omega(n\log n)$ is what optimality means here.

## Main Results

- **Theorem 1.12 (general spin systems, MLSI from spectral independence)**: for $G$ of max degree $\le\Delta$, totally-connected $\mu$ that is $b$-marginally bounded and $\eta$-spectrally independent, the Glauber dynamics satisfies the modified log-Sobolev inequality with constant $1/(C_1 n)$ where $C_1 = (\Delta/b)^{O(\eta/b^2 + 1)}$, hence
$$T_{\mathrm{mix}}(P_{\mathrm{GL}},\varepsilon) = (\Delta/b)^{O(\eta/b^2 + 1)} \cdot O\!\big(n\log(n/\varepsilon)\big).$$
The explicit constant (Remark 1.13) is $C_1 = \frac{18\log(1/b)}{b^4}\big(\frac{24\Delta}{b^2}\big)^{4\eta/b^2 + 1}$.
- **Theorem 1.1 (antiferromagnetic 2-spin)**: for any $\Delta\ge 3$, $\delta\in(0,1)$, and any antiferromagnetic 2-spin system that is *up-to-$\Delta$ unique with gap $\delta$*, mixing time $\le C(\Delta,\delta,\beta,\gamma,\lambda)\,n\log n$.
- **Theorem 1.2 (hardcore)**: for all $\Delta\ge 3$ and $\delta\in(0,1)$, hardcore on any $G$ of max degree $\le\Delta$ at $\lambda\le(1-\delta)\lambda_c(\Delta)$ has $T_{\mathrm{mix}}\le C(\Delta,\delta)\,n\log n$ with $C(\Delta,\delta) \sim \Delta^{O(\Delta^2/\delta)}$ (Remark 1.3).
- **Theorem 1.4 (Ising)**: for all $\Delta\ge 3$, $\delta\in(0,1)$, and $\beta\in[(\Delta-2+\delta)/(\Delta-\delta),(\Delta-\delta)/(\Delta-2+\delta)]$ with any external field $\lambda>0$: $T_{\mathrm{mix}}\le C(\Delta,\delta)\,n\log n$. Remark 1.5: $C=\Delta^{O(1/\delta)}$ suffices, giving polynomial mixing on **unbounded-degree** graphs.
- **Theorem 1.6 ($q$-colorings on triangle-free graphs)**: for $q\ge(\alpha^*+\delta)\Delta$ with $\alpha^*\approx 1.763$ the unique solution to $x=\exp(1/x)$, $T_{\mathrm{mix}}\le Cn\log n$.
- **Theorem 1.7 (monomer-dimer matchings)**: for any $\Delta\ge 3$ and fugacity $\lambda>0$, $T_{\mathrm{mix}}\le C(\Delta,\lambda)\,m\log n$ on graphs with $m$ edges (improves the prior $O(n^2 m\log n)$ of Jerrum–Sinclair).
- **Theorem 1.19 (general simplicial complexes)**: for any pure $n$-dim weighted simplicial complex $(\mathfrak{X},w)$ that is $(b_0,\dots,b_{n-1})$-marginally bounded and $(\zeta_0,\dots,\zeta_{n-2})$-local spectral expander, the order-$(s,r)$ down-up walk satisfies an MLSI with explicit constant $\kappa(r,s)$, giving $T_{\mathrm{mix}}\le\lceil\kappa^{-1}(\log\log(1/\pi_s^*)+\log(1/2\varepsilon^2))\rceil$. Generalizes both Cryan–Guo–Mousa entropy contraction (strongly log-concave) and Alev–Lau Poincaré bounds.
- **Optimality side**: $O(n\log n)$ is tight by Hayes–Sinclair $\Omega(n\log n)$ lower bound; the antiferromagnetic results match Sly–Sun NP-hardness past $\lambda_c(\Delta)$ / $\beta_c(\Delta)$ exactly.

## Method

(1) Reduce Glauber to a $\Theta(n)$-block heat-bath dynamics on the pure weighted simplicial complex $\mathfrak{X}^\mu$ associated with $\mu$, then prove an *approximate tensorization (uniform block factorization) of entropy* for that block dynamics; transfer to single-site Glauber by a comparison argument (Subsection 2.1). (2) Establish the block factorization by adapting the Alev–Lau local-to-global scheme from spectral gap to MLSI: marginal-boundedness $b$ + local spectral expansion $\zeta_k\le\eta/(n-k-1)$ (which follows from $\eta$-spectral independence by Claim 1.18 = AJKPV Theorem 8) imply contraction of relative entropy at every level of the simplicial complex, with rate $\alpha_k = \max\{1-4\zeta_k/(b_k^2(s-k)^2),\,(1-\zeta_k)/(4+2\log(1/(2 b_k b_{k+1})))\}$. (3) Apply known spectral-independence bounds — Weitz SAW-tree + correlation decay (hardcore/Ising), Chen–Galanis–Goldberg–Lapinskas–Lu–Stefankovic–Vigoda for $q$-colourings, a new monomer-dimer bound — to convert spectral independence into the $C_1=(\Delta/b)^{O(\eta/b^2+1)}$ constant.

## Quantum vs Classical

- **Baseline**: Anari–Liu–Oveis Gharan (FOCS 2020) achieved polynomial-time mixing in the same regime but with $\widetilde{O}(n^{2+C(\delta)})$ runtime via a spectral-gap argument and *no* MLSI; Weitz (STOC 2006) gave a deterministic FPTAS at quasi-polynomial $(n/\varepsilon)^{O(\log\Delta)}$. CLV21 closes both: optimal $O(n\log n)$ randomized mixing matching the Hayes–Sinclair lower bound, hence an $\widetilde{O}(n^2)$ FPRAS for the partition function in the same regime.
- **Gap (vs ALOG20)**: $n^{2+C(\delta)}\to n\log n$ — superpolynomial in $n$ to optimal in $n$; the $\Delta$-dependence is now confined to the multiplicative constant $C(\Delta,\delta)\sim\Delta^{O(\Delta^2/\delta)}$, still exponential in $\Delta$ but polynomial for any fixed bounded-degree family.
- **Source of the difference**: replacing the spectral gap (variance contraction) by entropy contraction. The spectral-gap route loses a $\log(1/\pi_*) = \Omega(n)$ factor on huge state spaces; the MLSI route saves it. Mechanistically, this comes from upgrading Alev–Lau's variance local-to-global lemma to an entropy local-to-global lemma on pure weighted simplicial complexes — exactly the step ALOG20 left open.
- **Quantum side**: paper makes no quantum claim. For the thesis Review chapter this is a *purely classical* benchmark in the **commuting × high-T (uniqueness)** cell. Quantum Gibbs samplers for diagonal classical Hamiltonians inherit the same $n\log n$ floor at best and pay a $\Omega(\log(1/\varepsilon))$ block-encoding overhead; there is no known regime in this cell where a quantum sampler beats CLV21. Any plausible quantum advantage in the thesis must come from **non-commuting × low-T** cells where neither spectral nor entropic independence have a known analogue (the matrix MLSI literature for non-commuting Lindbladians is far weaker — Kastoryano–Temme, Bardet et al. — see [CHECK] thesis prelim chapter).
- **Caveat**: the constant $C(\Delta,\delta) \sim \Delta^{O(\Delta^2/\delta)}$ is doubly-bad as $\delta\downarrow 0$ and $\Delta\uparrow\infty$; entropic independence (Anari–Jain–Koehler–Pham–Vuong 2021) and Chen–Eldan localization schemes (2022) later sharpen the $\Delta$-dependence and remove it entirely for Ising. The MLSI is *not* in general a standard log-Sobolev inequality, but the $O(n\log n)$ mixing-time conclusion is the same.

## Implications for Quantum Advantage

- **Regime cell**: classical Ising / hardcore / colourings × **high-T (tree-uniqueness)** × bounded-degree × classical spin. Top-left of the thesis grid (`classical-gibbs-sampling-corpus.md` Section 6).
- **What this changes**: pins the **classical ceiling** in the entire tree-uniqueness cell at the optimal $O(n\log n)$, with constants depending only on $(\Delta,\delta)$ and matching the Hayes–Sinclair lower bound and Sly–Sun hardness boundary on both sides. After this paper, the upper-bound side of this cell is essentially sharp.
- **Promising or not (for quantum advantage)**: explicitly *not promising*. This is the canonical "classical wins or ties" cell — quantum Gibbs samplers for classical commuting Hamiltonians reduce to a quantum implementation of Glauber and inherit the same $n\log n$ scaling, with extra $\log(1/\varepsilon)$ block-encoding overhead working against any speedup. The thesis Review chapter should cite this paper as the *classical benchmark* against which any high-T quantum claim is measured — and use the absence of any classical-side slack here as the structural reason to focus on non-commuting / low-T cells (see `classical-gibbs-sampling-corpus.md` punchline).

## Open Questions / Limitations

- The constant $C_1=(\Delta/b)^{O(\eta/b^2+1)}$ blows up as $\delta\downarrow 0$ at the uniqueness boundary; the conjectured optimum is $\mathrm{poly}(1/\delta)$, achieved later by Anari–Jain–Koehler–Pham–Vuong (entropic independence) for hardcore / Ising and by Chen–Eldan (localization schemes) for graphical Ising with no field dependence.
- Ferromagnetic 2-spin coverage matches prior best (Mossel–Sly, Guo–Liu, Liu–Lu–Yin) at $O(n\log n)$ but does not extend further; the ferromagnetic regime is not the bottleneck.
- $q$-colorings result requires triangle-free graphs; the general $q\ge\alpha^*\Delta$ conjecture remains open.
- Block factorization with constant block size $\Theta(n)$ (not $O(1)$) — the comparison to single-site dynamics costs the $\Delta^{O(\Delta^2/\delta)}$ factor.
- Pure classical: no statement about quantum spin systems. The simplicial-complex / HDX framework has no obvious lift to non-commuting Hamiltonians; trying to define an "operator-valued" pairwise influence runs into ordering / Lindbladian-detailed-balance issues that the matrix-MLSI literature has only partially resolved.
- Remark 1.14: a follow-up by Lau (2022) reformulates the proof without simplicial complexes and brings the constant down to $C_1=(\Delta/b)^{O(\eta/b)+1}$.

## Connections

- **Promotes** Anari–Liu–Oveis Gharan (2020) `anari-liu-oveisgharan-2020-spectral-independence-hardcore`: same spectral-independence hypothesis, gap bound $\to$ MLSI, mixing time $n^{2+C(\delta)}\to n\log n$.
- **Matches** Sly (2010) / Sly–Sun (2014) hardness exactly: rapid mixing up to $\lambda_c(\Delta)$ / $\beta_c(\Delta)$, NP-hardness past it. Together they cement "computational threshold = uniqueness threshold = optimal-mixing threshold" for hardcore and antiferromagnetic Ising.
- **Sharpens** Weitz (2006) `weitz-2006-counting-independent-sets-tree-threshold`: replaces $(n/\varepsilon)^{O(\log\Delta)}$ FPTAS by $\widetilde{O}(n^2)$ FPRAS in the same regime.
- **Extends** Cryan–Guo–Mousa (2019) entropy-factorization for strongly log-concave distributions to general pure weighted simplicial complexes (Theorem 1.19).
- **Generalized by** Chen–Eldan (2022) `chen-eldan-2022-localization-schemes`: localization-schemes framework re-proves CLV21 via martingales of measures without HDX local-to-global, and removes external-field dependence for graphical Ising.
- **Extends to** near-critical ferromagnetic Ising via Bauerschmidt–Dagallier (2024) `bauerschmidt-dagallier-2024-near-critical-ising-lsi`, which uses different (multiscale Bakry–Émery) machinery to reach $T_c$.
- **Companion** to Chen–Liu–Vigoda (2022 SICOMP, arXiv:2004.09083): contraction-based proof of spectral independence for antiferromagnetic 2-spin systems up to uniqueness, supplying the spectral-independence input that Theorem 1.12 here consumes.
- Within the thesis grid (`classical-gibbs-sampling-corpus.md`): this is the **classical Ising/Potts × high-T** Tier-1 anchor in Section 6, alongside Anari–Liu–Oveis Gharan and Eldan–Koehler–Zeitouni. The "where can quantum first beat classical?" punchline cell (non-stoquastic × low-T) sits structurally outside the reach of HDX / spectral independence — consistent with that being the strongest candidate for quantum advantage.
