---
paper: "Syljuåsen, Sandvik (2002)"
title: "Quantum Monte Carlo with directed loops"
arxiv: "cond-mat/0202316"
year: 2002
venue: "Phys. Rev. E 66:046701"
pdf: "n/a"

temperature: [all-T]
commutativity: [noncommuting, stoquastic, sign-problem-free]
locality: [geometric-local, local]
particle-statistics: [spin]
hamiltonian-models: [Heisenberg, anisotropic-Heisenberg, XXZ, transverse-field-Ising, XY]

paradigm: [classical-MCMC, classical-other]
quantum-or-classical: [classical]

result-type: [structural, runtime-comparison, gap-bound]
key-scaling: "per-MC-step cost $\\sim O(N\\beta)$ (one operator-string sweep + one directed loop of expected length $\\sim N\\beta$); integrated autocorrelation $\\tau_{\\mathrm{int}}$ shown empirically $O(1)$ in $N\\beta$ for the isotropic Heisenberg AFM and reduced by $1$–$2$ orders of magnitude vs prior loop / operator-loop algorithms in the gapless / weak-field regime"

related: ["gull-et-al-2011-ct-qmc-impurity-models", "prokofev-svistunov-tupitsyn-1998-worm-algorithm", "bravyi-gosset-2017-quantum-ferromagnets"]
---

# Quantum Monte Carlo with directed loops — Syljuåsen, Sandvik (2002)

**One-sentence takeaway**: Directed-loop SSE is the workhorse classical Gibbs sampler for sign-problem-free quantum spin lattices — a single linear-program-defined update generalises Sandvik's 1991 loop algorithm to the entire anisotropic Heisenberg / XXZ / TFIM family at arbitrary field, with autocorrelations short enough that 2D Heisenberg systems up to $L=64$ are routinely thermalised at $\beta J \sim 10^3$. Any quantum Gibbs-sampler claim of advantage on stoquastic spin lattices must beat this baseline in practice, not just asymptotically.

## Setting

- **Object sampled**: thermal expectation values $\langle A\rangle = \tfrac{1}{Z}\,\mathrm{tr}\,e^{-\beta H} A$ for sign-problem-free quantum spin Hamiltonians on a finite lattice. The two representations addressed jointly are **stochastic series expansion** (SSE: $Z = \sum_n (-\beta)^n/n!\,\mathrm{tr}\,H^n$ truncated at adaptive cutoff $M$) and **continuous-time world-line / path-integral** (the $M\to\infty$ limit of the same expansion).
- **Model class**: anisotropic spin-$1/2$ Heisenberg antiferromagnet in a longitudinal field,
$$H = \sum_{\langle ij\rangle} \big[J_x(S^x_i S^x_j + S^y_i S^y_j) + J_z S^z_i S^z_j\big] - h\sum_i S^z_i,$$
covering isotropic Heisenberg ($J_x=J_z=J$, $h=0$), XY ($J_z=0$), XXZ, and (by anisotropy limits) transverse-field Ising. The construction extends to any sign-free SSE-representable spin model with two-body bond operators.
- **Sign-free condition**: bipartite lattice with antiferromagnetic $J_x$ allows the sublattice rotation $S^{x,y}_i \to -S^{x,y}_i$ on one sublattice, so all off-diagonal matrix elements of $-H$ become non-negative — *stoquastic in the standard $S^z$ basis*. Frustrated antiferromagnets (triangular, kagome AFM) and fermions away from half-filling fall outside this class and incur a sign problem; the algorithm is silent on them.
- **Temperature regime**: any $T>0$. Truncation cutoff $M \propto N\beta$ is set adaptively; no Trotter error.

## Main Results

- **Directed-loop equations as a vertex-level linear program**: at each vertex of the operator string, an entering loop leg has four exit options (continue-straight, switch, reverse i.e. *back-track*, bounce). The detailed-balance constraints on the four directed-loop transition probabilities $a_{ij}$ form a linear program in the bond weights; the paper states it explicitly (Eqs. 23–27, §III) and gives closed-form solutions for XXZ in arbitrary $h$.
- **"Algorithmic discontinuities" eliminated**: prior loop algorithms (Evertz et al. 1993; Beard–Wiese 1996) and operator-loop SSE (Sandvik 1999) avoided back-tracking only at isolated parameter points (e.g. isotropic Heisenberg, $h=0$). Directed loops give a *continuous family* of valid update rules: the back-tracking probability turns on smoothly as one moves away from the special point and vanishes continuously upon return (Eqs. 31–33).
- **XY model is back-tracking-free for all $h$ up to saturation** (§IV.B) — the linear program admits a no-bounce, no-back-track solution on the entire phase line.
- **Empirical autocorrelations** (Figs. 6–8): on the 2D isotropic Heisenberg AFM at $\beta J = 16$, $L=16$, integrated autocorrelation of the staggered magnetisation drops from $\tau_{\mathrm{int}} \sim 50$ (operator-loop SSE) to $\tau_{\mathrm{int}} \sim 5$ (directed loops) — *one order of magnitude*, growing to two orders in the weak-field anisotropic regime where prior algorithms were nearly stuck.
- **Application: 2D Heisenberg magnetisation curve** (§V): on $L\times L$ lattices up to $L=64$, $\beta J = 2L$, the magnetisation steps between $S^z_{\mathrm{tot}}$ sectors are resolved cleanly enough to extract spin gaps, giving the transverse susceptibility $\chi_\perp = 0.0659 \pm 0.0002$. This is the *quantitative* benchmark a quantum sampler would replace.

## Method

The method represents $\mathrm{tr}\,e^{-\beta H}$ as a stochastic sum over operator strings $\{H_{a_p,b_p}\}_{p=1}^n$ with $H = -\sum_{a,b} H_{a,b}$ split into diagonal ($a=1$) and off-diagonal ($a=2$) bond operators. Updates are: (i) *diagonal update* — Metropolis insertion/removal of diagonal operators along the string; (ii) *directed-loop update* — pick a vertex leg, propagate a "worm head" through the operator string flipping spin states along its path, with the four transition probabilities at each vertex chosen to satisfy detailed balance. The novelty is solving the loop-equation linear program globally rather than at isolated symmetric points; back-tracking (the reverse-direction option) is the slack variable that makes the LP feasible on the whole XXZ + field plane. The combined update is ergodic and its directed-loop step is *non-local* in imaginary time, killing the critical slowing down ($z \approx 2$) that plagues local Metropolis on quantum spin systems.

## Quantum vs Classical

- **Role in the corpus**: defines the *practical classical frontier* for the **stoquastic × all-$T$** cell of the Review chapter's regime grid. SSE / directed-loop QMC is the algorithm against which "quantum advantage on Heisenberg / XXZ / TFIM" must be benchmarked — Hirsch–Fye, ED, and DMRG (outside 1D) are all weaker baselines.
- **Empirical scaling**: per-MC-step cost is dominated by one diagonal sweep (length $M \sim \langle n\rangle \sim N\beta\,\langle H_{\mathrm{bond}}\rangle$, cost $O(N\beta)$) plus one directed loop of expected length $\sim N\beta$. For the 2D Heisenberg AFM $\tau_{\mathrm{int}}$ is empirically $O(1)$ in $N\beta$ in the relevant regime (Fig. 8), giving a *total* sampling cost of $\widetilde{O}(N\beta)$ per independent sample — *linear* in the natural problem size.
- **Comparison to KMS Lindbladian quantum Gibbs samplers**: rigorous quantum samplers [Chen–Kastoryano–Gilyén; Ding et al.] give $\widetilde{O}(\beta\,\mathrm{poly}(N)/\Delta_{\mathrm{LS}})$ on any local Hamiltonian, *including the non-stoquastic ones where directed-loop QMC fails*. On the stoquastic spin lattices that directed loops solve, the quantum Lindbladian asymptotic is competitive in $\beta$ but pays a large constant prefactor (block-encoding overhead, gate counts) and a possibly worse $N$ dependence through $\Delta_{\mathrm{LS}}^{-1}$. **Empirically, the directed-loop algorithm wins on this cell — it is the canonical example of "classical wins inside stoquastic"**.
- **Source of the asymmetry**: the directed-loop update is a *non-local* move in imaginary time tailored to the world-line topology; it bypasses the spectral-gap bottleneck that local Glauber-style dynamics (and many quantum Lindbladian constructions) suffer from at low $T$ and near criticality. The structural reason it beats quantum samplers on stoquastic models is that classical world-line algorithms can directly sample the *positive measure* that stoquasticity exposes, while a quantum Lindbladian must instead simulate dissipative dynamics whose mixing is still controlled by a Hamiltonian gap.
- **Where it fails — sign problem**:
  - Frustrated AFMs (triangular / kagome Heisenberg, $J_1$–$J_2$ on the square lattice with $J_2 > 0$): bipartite trick fails, off-diagonal weights become negative, sign $\langle s\rangle \sim e^{-2\beta\Delta F}$ drops exponentially with $\beta$ and system size.
  - Itinerant fermions away from half filling: same Troyer–Wiese exponential.
  - Real-time / Keldysh: catastrophic dynamical sign problem.
- **Caveat — basis dependence**: the sign-free property is basis-dependent. Klassen–Marvian (2020) showed deciding whether a 2-local Hamiltonian admits a sign-free basis is NP-hard, so the directed-loop algorithm's applicability is in general not algorithmically detectable; practitioners use it on models with manifest stoquasticity.

## Implications for Quantum Advantage

- **Regime cell**: all-$T$ × stoquastic × geometric-local × spin. Directed-loop SSE *closes* the practical front in this cell. Together with the worm algorithm [Prokof'ev–Svistunov–Tupitsyn 1998] for bosonic / world-line lattice models, this is the empirical lower envelope of classical cost — the bar for quantum to clear.
- **What this changes**: makes quantum advantage on **stoquastic spin Hamiltonians at any temperature** unlikely in the foreseeable future. Quantum Gibbs samplers should target the *non-stoquastic* row of the regime grid (frustrated AFMs, sign-problematic fermions, doped Hubbard), where directed loops fail by Troyer–Wiese.
- **Sharpens the meaning of "stoquastic"**: stoquasticity in the standard basis is sufficient for directed-loop QMC to inherit polynomial scaling. The regime where quantum advantage is plausible is therefore *non-stoquastic in every accessible basis* (Klassen–Marvian's worst case), not merely "non-stoquastic in the local $S^z$ basis".
- **Promising or not**: not for the stoquastic cell. The directed-loop algorithm is mature, the constants are small, and the empirical scaling is essentially linear; a quantum sampler beating it would need a constant-factor or polynomial advantage in $\Delta_{\mathrm{LS}}^{-1}$, which is currently out of reach. The Review chapter should explicitly mark the all-$T$ × stoquastic × spin cell as classically dominated.

## Open Questions / Limitations

- **No proven mixing-time bound**: the paper gives empirical autocorrelation data only. There is no rigorous polynomial mixing-time bound on directed-loop SSE for the 2D Heisenberg AFM at $\beta \to \infty$; the absence of critical slowing down is conjectural and confirmed by simulation, not proof. This is the structural gap that quantum Gibbs samplers might exploit if a *worst-case rigorous* lower bound on directed-loop mixing ever appeared.
- **Linear-program degrees of freedom**: when the LP is underdetermined, a choice of solution must be made; the paper picks the no-bounce solution where available but does not optimise over the LP polytope. Extensions (Alet–Wessel–Troyer 2005) explore this systematically.
- **Sign problem unaddressed**: the algorithm is silent on frustrated and fermionic models; the *severity* of the sign problem on, e.g., $J_1$–$J_2$ Heisenberg under directed loops is not quantified here.
- **Higher-spin / multi-body interactions**: the loop equations were originally designed for spin-$1/2$ two-body bonds; spin-$1$ and ring-exchange extensions exist but require model-specific LP solutions. No general theorem guarantees feasibility.
- [CHECK] The claim "back-tracking probability vanishes continuously at the isotropic Heisenberg point" is stated and demonstrated graphically (Fig. 4) but not proven as a theorem in the paper.

## Connections

- **Predecessors / sister classical algorithms**:
  - Sandvik (1991) — first SSE for quantum spin systems; loop updates were limited to the isotropic Heisenberg point. Directed loops generalise this to arbitrary anisotropy and field.
  - Evertz, Lana, Marcu (1993); Beard, Wiese (1996) — continuous-time loop algorithms for the same model class. Directed loops unify the discrete-SSE and continuous-time loop frameworks under one detailed-balance scheme.
  - Sandvik (1999) — operator-loop SSE; directed loops eliminate the bounce-dominated regimes of operator loops.
  - **Prokof'ev–Svistunov–Tupitsyn (1998)** [`prokofev-svistunov-tupitsyn-1998-worm-algorithm`] — worm algorithm for world-line bosonic lattice models; the conceptual cousin (single propagating worm head, non-local imaginary-time update) and direct ancestor of "G-sector" sampling. Directed loops are the spin-lattice / SSE analogue.
- **Sister classical solver in a different cell**:
  - **Gull et al. (2011)** [`gull-et-al-2011-ct-qmc-impurity-models`] — CT-QMC for fermionic impurity models; the directed-loop / worm idea ported to *fermionic* impurities. Same ideology (sample a positive expansion exactly at sign-free corners; sign problem otherwise), different cell of the regime grid (non-commuting × intermediate-T × sparse).
- **Quantum-side counterpart on the same cell**:
  - **Bravyi–Gosset (2017)** [`bravyi-gosset-2017-quantum-ferromagnets`] — provable polynomial-time *classical* algorithm for ferromagnetic stoquastic models (XY ferromagnet). Complements directed loops by giving rigorous worst-case polynomial bounds for a sub-class where directed loops are empirically efficient. The combination of these two papers makes the stoquastic × all-$T$ cell extremely hard for quantum to beat.
- **Sign-problem hardness backbone**:
  - Troyer–Wiese (2005) — NP-hardness of curing the sign problem; explains why directed loops cannot extend to frustrated AFMs and frames the non-stoquastic row of the grid as the quantum-advantage frontier.
  - Klassen–Marvian (2020) — basis-change sign-curing is NP-hard, sharpening the boundary of "stoquastic" in the algorithmic sense.
- **Thesis Review chapter**: directed loops are the canonical "classical wins" entry for the stoquastic × all-$T$ × spin cell. The chapter's narrative — "where can quantum first beat classical?" — should cite this paper as the reason to *exclude* that cell from candidacy and direct attention to the non-stoquastic row.
