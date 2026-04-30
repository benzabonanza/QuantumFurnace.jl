---
paper: "Bravyi and Terhal (2010)"
title: "Complexity of Stoquastic Frustration-Free Hamiltonians"
arxiv: "0806.1746"
year: 2010
venue: "SIAM J. Comput. 39:1462-1485"
pdf: "supplementary-informations/classical-review/bravyi-terhal-2010-frustration-free-stoquastic.pdf"

temperature: [ground-state]
commutativity: [stoquastic, frustration-free]
locality: [k-local]
particle-statistics: [spin]
hamiltonian-models: [stoquastic-LH, stoquastic-k-SAT]

paradigm: [classical-MCMC, classical-other]
quantum-or-classical: [classical]

result-type: [complexity-upper, complexity-classification]
key-scaling: "classical sampler runs in $\\mathrm{poly}(n,\\delta^{-1})$ time, given $\\Delta(s) \\ge 1/\\mathrm{poly}(n)$ along the adiabatic path; sample $\\bar\\pi$ with $\\|\\bar\\pi - \\pi\\|_1 \\le \\delta$ where $\\pi(x) = |\\langle x|\\psi\\rangle|^2$"

related: ["bravyi-divincenzo-oliveira-terhal-2008-stoquastic-lh", "bravyi-gosset-2017-quantum-ferromagnets", "mann-helmuth-2021-high-T-partition-functions", "mann-helmuth-2023-low-T-partition-functions", "troyer-wiese-2005-sign-problem-hardness", "hangleiter-roth-nagaj-eisert-2020-easing-sign-problem", "klassen-marvian-2020-curing-sign-problem"]
---

# Complexity of Stoquastic Frustration-Free Hamiltonians — Bravyi, Terhal 2010

**One-sentence takeaway**: Adiabatic evolution along any smooth path of stoquastic frustration-free (SFF) Hamiltonians with $1/\mathrm{poly}(n)$ spectral gap can be classically simulated in randomized polynomial time, so the entire SFF-with-gap class — which includes coherent Gibbs states of any classical local Hamiltonian — is a *proven* classically tractable corner of quantum simulation and serves as the canonical "easy quantum case" benchmark for any quantum-Gibbs-sampling claim.

## Setting

- **Hamiltonian class**: $k$-local *stoquastic frustration-free* (SFF). $H = \sum_a H_a$ with $k=O(1)$, $H_a \succeq 0$, all off-diagonal entries of $H_a$ in the standard basis are real and non-positive (stoquastic), and the ground state $|\psi\rangle$ satisfies $H_a|\psi\rangle = 0$ for every $a$ (frustration-free). Norms $\|H_a\| \le \mathrm{poly}(n)$, $M = \mathrm{poly}(n)$ terms.
- **Target**: sample $x \in \{0,1\}^n$ from $\pi(x) = |\langle x|\psi\rangle|^2$, the *ground-state* probability distribution. By Perron-Frobenius the SFF ground state can be chosen non-negative in the standard basis, so $\pi$ is a genuine probability distribution and not a quasi-probability.
- **Temperature**: ground-state regime ($\beta = \infty$). Note however that *coherent Gibbs states* $|\pi_\beta\rangle = \mathcal{Z}^{-1/2}\sum_x e^{-\beta H_{\rm cl}(x)/2}|x\rangle$ of any classical local Hamiltonian $H_{\rm cl}$ are themselves the unique ground states of an explicit SFF Hamiltonian $H_\beta = \sum_j(\Gamma_j - X_j)$ with $\Gamma_j = X_j e^{-\beta H_{\rm cl}/2} X_j e^{\beta H_{\rm cl}/2}$ (Section 3.1, generalizing Verstraete et al. and Somma et al.) — so the result also covers classical Gibbs sampling at any $\beta$ via SFF ground states.
- **Adiabatic path conditions** (A0)-(A3): (A0) $H(s)$ SFF for all $s\in[0,1]$; (A1) $\|dH/ds\| \le J = \mathrm{poly}(n)$; (A2) gap $\Delta(s) \ge 1/\mathrm{poly}(n)$; (A3) the initial $H_{\rm in}$ admits a $\mathrm{poly}(n)$-time classical algorithm producing a basis vector with overlap $\ge 2^{-\mathrm{poly}(n)}$ on its ground state.

## Main Results

- **Theorem 1.4 (classical simulation of SFF adiabatic evolution)**: under (A0)-(A3), for any $\delta > 0$ there is a classical probabilistic algorithm sampling $x \sim \bar\pi$ with $\|\bar\pi - \pi\|_1 \le \delta$ in time $\mathrm{poly}(n,\delta^{-1})$.
- **Theorem 1.9 (stoquastic $k$-SAT $\in$ MA, for any $k=O(1)$)**: a quantum analogue of Cook-Levin — every problem in MA is Karp-reducible to stoquastic 6-SAT, the first non-trivial MA-complete promise problem. Combined with the prior MA-hardness of stoquastic 6-LH (Bravyi-DiVincenzo-Oliveira-Terhal 2008), stoquastic 6-SAT is **MA-complete**.
- **Proposition 1.5 (universality bound)**: any $\mathrm{poly}(n)$ quantum circuit output can be encoded as the ground state of an SFF Hamiltonian on $\mathrm{poly}(n,\delta^{-1})$ qubits — so SFF adiabatic evolution is BQP-universal *as a model*. Combined with Theorem 1.4 this shows the "frustration-free" hypothesis (not just stoquasticity) is doing the algorithmic work.

## Method

The core tool is a *random walk associated with a SFF Hamiltonian*. Given an SFF $H$ and its non-negative ground state $|\psi\rangle$, define the transition kernel $P_{x \to y} = \frac{\langle y|\psi\rangle}{\langle x|\psi\rangle}\langle y|G|x\rangle$ on the support $\mathcal{S}(\psi)$, where $G = I - \beta H$ is non-negative for $\beta$ small. This is reversible with stationary distribution $\pi(x) = \langle x|\psi\rangle^2$ and spectral gap $\ge \beta\Delta = 1/\mathrm{poly}(n)$. Each step is locally simulable because the spectral projector $\Pi_a$ onto $\ker H_a$ inherits non-negativity (Lemma 4.4, the key technical workhorse) so the ratios $\langle y|\psi\rangle/\langle x|\psi\rangle = \sqrt{\langle y|\Pi_a|y\rangle/\langle x|\Pi_a|x\rangle}$ are computable from $H_a$ alone. To follow the adiabatic path, $T = \mathrm{poly}(n)$ discrete walks $P^{(j)}$ are chained: a $t$-balanced sample for $P^{(j)}$ becomes a *warm start* for $P^{(j+1)}$ provided fidelities $F(\pi^{(j)},\pi^{(j+1)}) \ge 1 - J^2/(T^2\Delta^2)$ are close — i.e. $T \gg J^2/\Delta^2$ steps suffice (Section 2.2).

## Quantum vs Classical

- **Baseline (quantum side)**: SFF adiabatic evolution is BQP-universal as a model (Proposition 1.5). The "naive expectation" would be that it is therefore classically intractable.
- **Gap**: zero. The classical algorithm achieves $\mathrm{poly}(n,\delta^{-1})$ runtime, matching the quantum adiabatic algorithm modulo polynomial factors. There is no exponential separation in this corner.
- **Source of the difference**: stoquasticity gives a non-negative ground state by Perron-Frobenius; frustration-freeness lets one read off ground-state ratios *locally* from each $H_a$ without diagonalizing $H$. Together they reproduce the structural properties of a classical reversible Markov chain — the quantum amplitudes are real, non-negative, and locally accessible, so the simulation reduces to running a polynomial-gap reversible chain.
- **Caveat**: the proof needs *both* stoquasticity and frustration-freeness. Bravyi-DiVincenzo-Oliveira-Terhal showed adiabatic evolution along stoquastic-but-not-frustration-free paths reaches StoqAQC, whose classical complexity remains open (likely AM-hard, possibly StoqMA-hard at the upper end). Dropping stoquasticity collapses to BQP-complete; dropping frustration-freeness loses the local-readout of $\pi$.

## Implications for Quantum Advantage

- **Regime cell**: stoquastic $\times$ frustration-free $\times$ ground-state (effectively any-$\beta$ for coherent Gibbs of *classical* $H$) $\times$ gapped. In the corpus' 3$\times$3 grid this is the *stoquastic / low-T* cell at the ground-state end and equivalently the *classical / any-T* cell via the Section 3.1 mapping.
- **What this changes**: closes the SFF + gap regime as **classically tractable in $\mathrm{poly}(n,1/\Delta)$**. Any quantum-Gibbs-sampling claim of advantage must therefore avoid this corner — i.e. require non-stoquasticity, frustration, or vanishing gap. This is the canonical floor: a quantum sampler that beats classical here is impossible (modulo polynomial factors).
- **Promising or not**: *not promising* for quantum advantage. Together with Bravyi (QIC 2015) — which extends to ground-state expectation values via path-integral MCMC with provable autocorrelation bounds — the entire SFF class is provably classically polynomial. Quantum advantage in Gibbs sampling must come from one of: (i) non-stoquastic Hamiltonians (where Troyer-Wiese, Klassen-Marvian etc. block sign-curing), (ii) frustrated stoquastic with vanishing gap, or (iii) genuinely low-T regimes for non-classical $H$ where Mann-Helmuth-style polymer expansions break down.

## Open Questions / Limitations

- **Restricted to ground states**: the result samples $|\langle x|\psi\rangle|^2$, not finite-$\beta$ Gibbs distributions of *quantum* SFF Hamiltonians. The Section 3.1 trick covers Gibbs states of *classical* $H$ only; finite-$T$ Gibbs sampling of genuinely quantum SFF Hamiltonians is not addressed.
- **General stoquastic adiabatic evolution** (without frustration-freeness): classical complexity remains open. The paper conjectures it is unlikely to be BQP-complete but provides no upper bound.
- **Polynomial degree**: the simulation cost is $\mathrm{poly}(n,1/\Delta,1/\delta)$ but the explicit polynomial degree (in the spectral gap and accuracy) is large and not optimized; later work (Bravyi 2015) sharpens this for ground-state expectation values.
- **Frustration-free verification**: deciding whether a given stoquastic local $H$ is SFF requires testing $H|\psi\rangle = 0$, which is QMA$_1$-related; the paper raises this as an open problem.

## Connections

- **Builds on** Bravyi-DiVincenzo-Oliveira-Terhal (QIC 2008) [qf related: `bravyi-divincenzo-oliveira-terhal-2008-stoquastic-lh`]: defines stoquastic-LH as MA-hard / AM-contained; this paper sharpens to MA-complete in the frustration-free $k$-SAT subproblem.
- **Predecessor of** Bravyi (QIC 2015) [arXiv:1402.2295]: extends from sampling $\pi$ to ground-state expectation values via path-integral MCMC with provable autocorrelation bounds, polynomial in $1/\Delta$ and $n$.
- **Predecessor of** Bravyi-Gosset (PRL 2017) [qf related: `bravyi-gosset-2017-quantum-ferromagnets`]: same "stoquastic + structure $\Rightarrow$ classical poly-time" template, but at finite $T$ for the partition function of stoquastic transverse-field XY ferromagnets.
- **Complementary to** Verstraete et al. and Somma et al. (cited as [18,19] in the paper): they showed coherent classical Gibbs states are SFF ground states; this paper turns that structural fact into a polynomial-time classical sampler.
- **Defines the floor for** Mann-Helmuth (high-T 2021, low-T 2023) [qf related: `mann-helmuth-2021-high-T-partition-functions`, `mann-helmuth-2023-low-T-partition-functions`]: cluster-expansion-based extensions to *finite-T* stoquastic partition functions, which would have to beat this baseline to claim a new tractable corner.
- **Contrast with** Troyer-Wiese (PRL 2005) [qf related: `troyer-wiese-2005-sign-problem-hardness`]: that paper rules out a generic poly-time sign-curing; this paper shows that *when sign-free (stoquastic) and frustration-free hold*, polynomial classical simulation is unconditionally achievable. The two together delineate the stoquastic boundary.
- **Within the thesis grid**: cited in `classical-gibbs-sampling-corpus.md` §7 alongside Bravyi (2015) as the proof of "frustration-free + stoquastic + gap = classically tractable". The companion review for `anari-liu-oveisgharan-2020-spectral-independence-hardcore.md` plays the analogous role on the *commuting / classical* side; together they bracket the regions where quantum-Gibbs-sampling claims must *not* claim advantage.
