---
paper: "Kuwahara, Kato, Brandão (2020)"
title: "Clustering of Conditional Mutual Information for Quantum Gibbs States above a Threshold Temperature"
arxiv: "1910.09425"
year: 2020
venue: "Phys. Rev. Lett. 124:220601"
pdf: "n/a (arXiv only in repo)"

temperature: [high-T]
commutativity: [noncommuting]
locality: [geometric-local, k-local, long-range]
particle-statistics: [spin]
hamiltonian-models: [generic-local-spin, generic-finite-range, power-law-long-range]

paradigm: [structural, cluster-expansion]
quantum-or-classical: [comparison, structural]

result-type: [structural, decay-of-correlations, threshold-temperature, complexity-upper]
key-scaling: "$\\mathcal{I}_\\rho(A:C|B) \\le e\\,\\min(|\\partial A_r|,|\\partial C_r|)\\,(\\beta/\\beta_c)^{d_{A,C}/r}/(1-\\beta/\\beta_c)$ for finite range $r$, with $\\beta_c = 1/(8e^3 k)$; power-law $C_\\beta/d_{A,C}^{\\alpha}$ for $f(R)=R^{-\\alpha}$"

related: ["bakshi-liu-moitra-tang-2024-high-T-unentangled", "kuwahara-alhambra-anshu-2021-thermal-area-law", "mann-helmuth-2021-high-T-partition-functions"]
---

# Clustering of Conditional Mutual Information for Quantum Gibbs States above a Threshold Temperature — Kuwahara, Kato, Brandão (2020)

**One-sentence takeaway**: Above an explicit, system-size-independent threshold temperature $\beta_c = 1/(8e^3 k)$ (for finite-range $k$-local Hamiltonians), the conditional mutual information $\mathcal{I}_\rho(A:C|B)$ of any tripartite split of a Gibbs state decays exponentially in the separation $d_{A,C}$, so high-$T$ quantum Gibbs states are *approximate quantum Markov networks* and admit a quasi-local Petz-style recovery channel — the structural backbone of high-$T$ Gibbs-sampling tractability for arbitrary local quantum Hamiltonians.

## Setting

- **System**: $n$ spins (qudit Hilbert space, dimension $s$ per site) on a graph $G=(V,E)$, with graph-theoretic distance $d_{u,v}$ and surface region $\partial L_l := \{v\in L \mid d_{v,L^c}\le l\}$.
- **Hamiltonian** (Eq. 3): $H=\sum_{|X|\le k} h_X$, generic $k$-body interactions, normalised so that the one-site energy density $g\le 1$ via the locality bound $\sum_{X\ni v,\,\mathrm{diam}(X)\ge R} \|h_X\|\le f(R)$, with $f(1)\le g=1$ (Eq. 4).
- **Two interaction classes treated**: (i) *finite range* $r$, i.e. $f(R)=0$ for $R>r$ (Thm. 1); (ii) *power-law decay* $f(R)=R^{-\alpha}$, $\alpha>0$ (Thm. 3).
- **Object analysed**: Gibbs state $\rho = e^{-\beta H}/Z$ on a tripartite split $V_0 = ABC \subseteq V$ where $A$ and $C$ are *not concatenated* on the graph (i.e. $B$ separates them), with $D=V\setminus V_0$ (see Fig. 1 of the paper).
- **CMI definition** (Eq. 1):
  $$
  \mathcal{I}_\rho(A:C|B) := S(\rho^{AB}) + S(\rho^{BC}) - S(\rho^{ABC}) - S(\rho^B),
  $$
  with $S(\sigma)=-\mathrm{tr}(\sigma\log\sigma)$. This is the quantum analogue of the classical CMI characterising Hammersley–Clifford Markov networks.

## Main Results

- **Theorem 1 (exponential CMI decay, finite range)**: for $f(R)=0$ when $R>r$, the threshold $\beta < \beta_c := 1/(8e^3 k)$ (Eq. 6) implies
  $$
  \mathcal{I}_\rho(A:C|B) \;\le\; e \cdot \min(|\partial A_r|, |\partial C_r|) \cdot \frac{(\beta/\beta_c)^{d_{A,C}/r}}{1-\beta/\beta_c} \qquad (\text{Eq. 7}).
  $$
  Hence exponential decay length $\xi = r/\log(\beta_c/\beta)$, with surface-area prefactor (this strengthens the area law of Wolf et al. by giving its *saturation rate* — Eq. 9).
- **Theorem 2 (quasi-local effective Hamiltonian)**: the effective Hamiltonian $\tilde H_L := -\beta^{-1}\log\mathrm{tr}_{L^c}(e^{-\beta H})$ on a region $L$ admits a localisation $\tilde H_L = H_L + \Phi_L$, and $\Phi_L$ is approximated by a strictly localised $\Phi_{\partial L_l}$ supported within distance $l$ of $\partial L$, with
  $$
  \|\Phi_L - \Phi_{\partial L_l}\|\;\le\; \frac{e}{4\beta}\frac{(\beta/\beta_c)^{l/r}}{1-\beta/\beta_c}\,|\partial L_r|.
  $$
  Computation of $\Phi_L$ to additive error $n\varepsilon$ runs in $n(1/\varepsilon)^{O(k\log(d\,d_G))}$ (Eq. 12), $d_G$ = graph degree.
- **Corollary 1 (poly-time classical Gibbs simulation)**: thermodynamic quantities (local observables, $\log Z$, local entropy $-\mathrm{tr}(\rho^L\log\rho^L)$) are computable classically in $\mathrm{poly}(n)$ at additive error $1/\mathrm{poly}(n)$ — *the first poly-time FPATS for high-$T$ quantum Gibbs states*, refining the FPATS of Harrow–Mehraban–Soleimanifar from quasi-poly to poly.
- **Quantum Gibbs preparation consequence (cited from Brandão–Kastoryano via Thm. 5 of [39])**: on a $D$-dim lattice, there exists a CPTP map $\mathbb{F}=\mathbb{F}_{D+1}\cdots\mathbb{F}_2\mathbb{F}_1$, each $\mathbb{F}_s$ a tensor product of *quasi-local* CPTP channels acting on $O(\log^D n)$ spins, with $\|\mathbb{F}(\psi)-\rho\|_1 = 1/\mathrm{poly}(n)$ and gate count $\exp(O(\log^D n)) = n^{O(\log^{D-1} n)}$ (quasi-polynomial).
- **Theorem 3 (long-range / power law)**: for $f(R)=R^{-\alpha}$, $\alpha>0$, under $\beta < \beta_c/11$ and $d_{A,C}\ge 2\alpha$,
  $$
  \mathcal{I}_\rho(A:C|B)\le \beta\min(|A|,|C|)\,\frac{C_\beta}{d_{A,C}^{\alpha}},\qquad C_\beta = \frac{11e^{1/k}/\beta_c}{1-11\beta/\beta_c}\quad(\text{Eq. 14}).
  $$
  Power-law (rather than exponential) decay; first general clustering theorem at finite $T$ for *quantum* long-range models.
- **Operational meaning (Petz / Fawzi–Renner)**: $\mathcal{I}_\rho(A:C|B)\le \varepsilon$ implies (via Fawzi–Renner) existence of a *recovery channel* $\tau_{B\to BC}$ acting only on $B$ with $\|\tau_{B\to BC}(\rho^{AB})-\rho^{ABC}\|_1\le \sqrt{\varepsilon\log 2}$ (Eq. 5) — i.e. the global Gibbs state is reconstructible from its local marginals by acting on the separator $B$ alone, the quantum analogue of the Hammersley–Clifford Markov property.

## Method

Generalised cluster expansion of $\log\mathrm{tr}_{L^c}(e^{-\beta H_{\vec a}})$ as an analytic function of source parameters $\vec a$ over interaction terms, with the CMI itself parametrised as a multi-derivative of such logarithms. The technical novelty is a systematic bound on multi-derivatives of *logarithms of reduced density matrices* — the standard cluster expansion handles $\log Z$ but not partial-trace logarithms, which is exactly what CMI requires. Convergence of the expansion holds for $\beta < \beta_c = 1/(8e^3 k)$, giving the same regime as Harrow–Mehraban–Soleimanifar zero-freeness but extended from $\log Z$ to local effective Hamiltonians and CMI.

## Quantum vs Classical

- **Classical baseline**: Hammersley–Clifford theorem (1971) — classical Gibbs states on graphs are *exactly* Markov networks, with CMI vanishing identically when $A,C$ are separated by $B$. Classical high-$T$ clustering of correlations: Gross (1979), Park–Yoo (1995), Ueltschi (2004); Kliesch–Gogolin–Kastoryano–Riera–Eisert (2014) for "locality of temperature".
- **Quantum baseline at the time**: Araki (1969) gave 1D exponential clustering at all $T$ for short-range models; Wolf–Verstraete–Hastings–Cirac (2008) gave the thermal area law $\mathcal{I}_\rho(A:B)\le c\beta|\partial A|$ at all $T$ but *no decay rate*; Kato–Brandão (2019) proved CMI decays *subexponentially* in 1D at any $T$. No general higher-dimensional quantum clustering of CMI was known.
- **Gap**: structural — *no algorithmic super-polynomial gap*. The result *enables* polynomial-time classical FPATS at high $T$ (Corollary 1) and quasi-polynomial-time quantum Gibbs preparation (via Brandão–Kastoryano), so quantum and classical complexities are within a quasi-polynomial factor inside this regime cell.
- **Source of the difference vs commuting / classical**: in commuting (or classical) Gibbs states the effective Hamiltonian $\Phi_L$ is *exactly* localised on $\partial L$, giving the *exact* Markov property [36, 37]. For non-commuting $H$, $\Phi_L$ leaks into the bulk of $L$; this paper shows the leakage is *exponentially small* in the bulk distance whenever $\beta<\beta_c$. Below $\beta_c$, generalised cluster expansion fails and topological order (e.g. 4D toric code at finite $T$, Hastings) can produce $\mathcal{I}_\rho(A:C|B)=\Theta(1)$ globally.
- **Caveat**: the threshold $\beta_c=1/(8e^3 k)$ scales with locality $k$ (and via $C_\beta$ with $\alpha,\beta_c$ in the long-range case). It is *not* sharp: physical critical temperatures can be polynomially larger; matching upper and lower thresholds remains open. Long-range result requires $\alpha>0$ and tighter threshold $\beta<\beta_c/11$.

## Implications for Quantum Advantage

- **Regime cell**: high-T $\times$ noncommuting $\times$ (geometric-local *or* power-law long-range) $\times$ spin. Anchors the high-T row of the corpus 3$\times$3 grid as a *structural* tractability cell.
- **What this changes**: *opens* upper bounds. For the first time, a rigorous polynomial-time *classical* algorithm exists for thermodynamic quantities of generic local quantum Hamiltonians at high $T$ (Corollary 1), and a quasi-poly-time *quantum* preparation algorithm follows by combining with Brandão–Kastoryano's finite-correlation-length preparation. Both rest on the same structural fact: high-$T$ Gibbs states are quasi-local Markov networks.
- **Promising or not**: **negative for advantage in this cell.** Together with Bakshi–Liu–Moitra–Tang 2024 (separability + classical poly-time preparation up to $\beta_c \sim 1/(\mathfrak{d}\mathfrak{K})$), Mann–Helmuth (poly partition function), and Harrow–Mehraban–Soleimanifar (zero-free FPATS), this paper is one of four pillars showing the high-T row is structurally classical: *the state itself has only quasi-local quantum correlations, and they can be reconstructed by local recovery channels*. Reinforces the corpus framing that quantum advantage lives below $\beta_c$ — at intermediate / low $T$ where CMI decay is no longer guaranteed.
- **Unique contribution vs Bakshi et al.**: this paper gives *exponential CMI decay* and a *Markov-network / recovery-channel structure*; Bakshi et al. give *unentanglement* (separability into product states). The two are complementary structural statements about the same regime — Bakshi et al. is stronger when it applies (separability $\Rightarrow$ CMI decay), but Kuwahara–Kato–Brandão covers a wider Hamiltonian class (any $k$-local, including long-range $\alpha>0$) at a possibly different threshold.

## Open Questions / Limitations

- **Threshold sharpness**: $\beta_c = 1/(8e^3 k)$ is far below physical phase-transition temperatures. Whether CMI clustering persists up to the *true* one-phase boundary (in 2D and higher) is open. [CHECK]
- **Low-$T$ Markov property**: explicitly flagged in the "Future perspective" — the generalised cluster expansion technique breaks down at low $T$, where topological order in $D\ge 4$ can defeat the Markov property at finite $T$ (Hastings 2011). No replacement technique known.
- **Markov property at non-generic $V_0$**: even for short-range commuting Hamiltonians, the Markov property can fail for special tripartitions (e.g. cluster state on a special $V_0$, Refs [42,43]); the present paper proves *approximate* Markov property holds for *arbitrary* $V_0$ but only at $\beta<\beta_c$.
- **Long-range threshold**: $\beta < \beta_c/11$ is a factor-11 weaker than the short-range threshold; whether this is intrinsic or an artefact is unclear. Power-law decay rate $1/d^\alpha$ in CMI matches the interaction tail — possibly tight.
- **Non-spin / fermionic systems**: setup is bosonic spin; extension to fermionic Gibbs states (Hubbard, lattice gauge) would close one of the gaps in the corpus high-T row.
- **Beyond polynomial classical algorithm**: Corollary 1's classical algorithm gives $\mathrm{poly}(n)$ at fixed $1/\mathrm{poly}(n)$ error but the exponent is $O(k\log(d\,d_G))$, larger than Bakshi et al.'s $\widetilde O(n^{6+\log\mathfrak{d}/\log(\beta_c/\beta)})$. Direct comparison on the same Hamiltonian class would be informative.

## Connections

- **Bakshi–Liu–Moitra–Tang 2024**: stronger structural result (separability) on a slightly different Hamiltonian class (low-intersection rather than $k$-local with norm bound); cluster-expansion lineage shared. The two papers together close the high-T row.
- **Kuwahara–Alhambra–Anshu 2021 (thermal area law)**: improves the area law from $O(\beta)$ to $\tilde O(\beta^{2/3})$ in 1D and gives quasi-linear-time MPS-based 1D Gibbs construction. Complements this paper: KKB gives the *higher-dimensional* CMI clustering, KAA gives the *1D-tight* area-law bond-dimension scaling.
- **Mann–Helmuth 2021 (high-T quantum partition functions)**: same regime, complementary endpoint — MH gives poly-time partition function; KKB Corollary 1 sharpens this to thermodynamic *and* local-observable computation via the effective-Hamiltonian construction.
- **Harrow–Mehraban–Soleimanifar 2020**: same threshold $\beta_c\sim 1/k$ via cluster expansion / zero-freeness, but quasi-poly time and only $\log Z$. KKB upgrades to poly-time via local effective Hamiltonians.
- **Brandão–Kastoryano 2019 ([39])**: provides the dissipative quantum Gibbs sampler whose efficiency on $D$-dim lattices follows from KKB's clustering — directly cited as the route from CMI clustering to a quasi-polynomial-time quantum sampler.
- **Fawzi–Renner 2015 ([33])**: provides the operational meaning of CMI as recoverability error — the bridge from KKB's CMI bound to existence of a local recovery / Petz channel.
- **Kato–Brandão 2019 (1D approximate Markov)**: subexponential 1D CMI clustering at any $T$; KKB strengthens to *exponential* decay in arbitrary dimension at the cost of the high-$T$ assumption.
- **Hastings 2011 (4D topological order at finite $T$)**: explicit obstruction to extending CMI clustering below $\beta_c$ — quantum memories at finite $T$ violate the Markov property globally.
