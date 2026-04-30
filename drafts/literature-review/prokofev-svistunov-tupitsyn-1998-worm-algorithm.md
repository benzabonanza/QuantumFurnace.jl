---
paper: "Prokof'ev, Svistunov, Tupitsyn (1998)"
title: "'Worm' algorithm in quantum Monte Carlo simulations"
arxiv: "cond-mat/9703200"
year: 1998
venue: "Phys. Lett. A 238:253–257 (long version: JETP 87:310, 1998)"
pdf: "n/a"

temperature: [all-T]
commutativity: [noncommuting]
locality: [geometric-local, lattice, sparse]
particle-statistics: [bosonic]
hamiltonian-models: [bose-hubbard, lattice-bosons, quantum-rotor, XXZ-easy-plane]

paradigm: [classical-other]
quantum-or-classical: [classical]

result-type: [structural, runtime-comparison]
key-scaling: "per-update cost $O(1)$ on the extended (broken-worldline) configuration space; integrated autocorrelation time empirically $\\sim O(L^z)$ with $z \\approx 0.5$–$1$ in the superfluid phase and $z \\approx 2$ near the Bose–Hubbard critical point on $L^d$ lattices; wall-clock comparable to or better than directed-loop SSE for $U(1)$-symmetric bosonic / sign-free models"

related: ["gull-et-al-2011-ct-qmc-impurity-models", "syljuasen-sandvik-2002-directed-loops-qmc", "bravyi-gosset-2017-quantum-ferromagnets"]
---

# "Worm" algorithm in quantum Monte Carlo simulations — Prokof'ev, Svistunov, Tupitsyn (1998)

**One-sentence takeaway**: The worm algorithm samples the *extended* configuration space of broken world-lines — propagating two free endpoints ("Ira" and "Masha") through space-time — which makes the off-diagonal Green's function, winding-number sectors, and grand-canonical particle-number changes accessible by purely **local** moves, and is the canonical best-in-class classical Gibbs sampler for sign-free bosonic / $U(1)$-symmetric lattice models (e.g. Bose–Hubbard).

## Setting

- **Object sampled**: thermal expectation values $\langle A\rangle = \tfrac{1}{Z}\,\mathrm{tr}\,e^{-\beta H} A$ for lattice bosons with $U(1)$ particle-number symmetry, primarily the Bose–Hubbard Hamiltonian $H = -t\sum_{\langle i,j\rangle}(b_i^\dagger b_j + \mathrm{h.c.}) + \tfrac{U}{2}\sum_i n_i(n_i-1) - \mu\sum_i n_i$ at finite $T$ (originally demonstrated on the 1D pure and disordered Bose–Hubbard model).
- **Representation**: continuous-imaginary-time **world-line** path integral. Each boson is a closed world-line in $(\text{site} \times [0,\beta))$; the partition function is a sum over closed configurations of these lines, weighted by the exponentiated diagonal energy and the kinetic-hop matrix elements.
- **Configuration space — the trick**: enlarge the partition function to
$$Z_W = Z + \kappa \sum_{i,j}\int_0^\beta d\tau\,\langle b_i^\dagger(\tau)\,b_j(0)\rangle\,e^{-\beta H},$$
i.e. include configurations with **exactly one broken world-line** that ends at two free space-time points (the "Ira" and "Masha" endpoints, called $M$ and $I$ here for "creation" and "annihilation"). The ratio of weights at fixed $i,j,\tau$ vs. closed sector samples directly the imaginary-time single-particle Green's function $G(i-j,\tau)$.
- **Access model**: classical sampler — no oracle, no quantum subroutine. Inputs are the lattice graph, hopping amplitudes, and a stream of pseudo-random numbers; outputs are stochastic estimates of any operator that is diagonal *or* a single creation–annihilation pair (so $\langle n\rangle$, $\langle n^2\rangle$, kinetic energy, $G(\tau)$, superfluid stiffness all read off from a single run).
- **Temperature regime**: any $T > 0$. No Trotter discretisation error (the kink times are continuous in $[0,\beta)$); $\beta\to\infty$ extrapolation is natural and standard.

## Main Results

- **Algorithmic construction**: defines a small move set on $Z_W$ such that (i) when no worm is present, the chain samples the physical partition function $Z$; (ii) when a worm is present, the chain samples the Green-function-extended sector with the correct relative weight; (iii) all moves are **local** in space-time. The move set (essentially universal across worm-style codes) is:
  - **Create / Annihilate worm**: insert a coincident pair $(M,I)$ at a random site and time, or remove them when they meet — couples the two sectors.
  - **Move worm head in time** ($M$ or $I$ shifts by $\Delta\tau$ along its world-line): drives the imaginary-time argument of $G$.
  - **Insert / Remove kink at the worm head**: the worm hops to a neighbouring site by inserting (or absorbing) a hopping vertex $b_i^\dagger b_j$, costing energy $\Delta E$ over the segment and weighted by the hopping amplitude. This is the move that grows or shrinks the spatial separation in $G(i-j,\tau)$ and that **changes the winding number** when the worm wraps the periodic spatial direction.
  - **Reconnect through existing kink**: ergodicity-completing move that lets the worm "swallow" an existing hop on a neighbouring world-line.
- **What this buys** (the headline claims):
  1. **Ergodicity across winding-number / particle-number sectors** — closed-world-line chains (e.g. Beard–Wiese loop) cannot change total winding without global moves; the worm changes winding **one local kink at a time**, removing the topological-sector bottleneck. Critical for measuring the superfluid stiffness $\rho_s = \langle W^2\rangle/(\beta L^{d-2})$.
  2. **Direct grand-canonical sampling** — the create/annihilate worm move changes total particle number by $\pm 1$ when $M$ and $I$ are inserted at different times, so the chemical potential $\mu$ enters the simulation natively (no Legendre transform, no canonical-ensemble bias). Essential for the Bose–Hubbard phase diagram $(\mu, t/U)$ where the Mott lobes are defined in $\mu$.
  3. **Green's function for free** — every move during the open-sector sampling produces a histogram of $(i-j,\tau)$ that *is* $G(i-j,\tau)$ up to normalization. No re-weighting, no separate simulation, no analytic continuation needed for imaginary-time data.
- **Empirical scaling** (1D Bose–Hubbard demonstration in the original Letter; broader literature since): per-update cost $O(1)$, integrated autocorrelation $\tau_{\mathrm{int}} \lesssim L$ in the superfluid phase, growing to $\sim L^2$ near the Mott / superfluid-insulator critical point — matching or beating directed-loop SSE on the same systems and natively continuous-time.
- **Disorder**: the Letter explicitly demonstrates simulations of *disordered* 1D bosons (random $\mu_i$), where local-move ergodicity is critical because non-local cluster moves do not respect random potentials.

## Method

The key idea is to *replace global topological moves by local moves on an enlarged configuration space*. Closed-world-line algorithms must change winding by collective updates (loops, swaps), which become rejected with probability $\to 1$ in interesting parameter regimes; by inserting two free worm endpoints and letting them random-walk through space-time, every elementary local move *can in principle* lead to a winding change once the head has wrapped the lattice. The same enlargement makes off-diagonal correlators $\langle b^\dagger(\tau) b(0)\rangle$ structural quantities of the configuration space rather than estimators that vanish on closed configurations. This is a direct quantum analogue of the **classical** worm algorithm for high-temperature graph expansions of $O(n)$ models (Prokof'ev–Svistunov 2001), and the same structural reason it works there — a Markov chain on broken paths that recombines into closed paths only at the endpoints — drives the quantum case here.

## Quantum vs Classical

- **Role in the corpus**: this paper defines the *classical upper frontier* for the **noncommuting × all-T × bosonic-lattice** cell of the 3$\times$3 grid (rows × temperatures, with bosons sliced inside the noncommuting / sign-free row). Anywhere a quantum lattice problem is sign-free in the world-line representation — Bose–Hubbard, $U(1)$-symmetric bosonic models, $XXZ$ with easy-plane anisotropy in the $S^z$ basis (via the Holstein–Primakoff / hard-core boson map) — worm is the practitioner's gold-standard sampler.
- **Baselines it supersedes**:
  - **Loop algorithm** (Evertz–Lana–Marcu / Beard–Wiese): closed-world-line cluster updates, sign-free in the same regime, but **cannot change winding without a global cluster wrapping the torus** and forces canonical-ensemble simulation. Worm dominates loops on every ergodicity / observable axis except wall-clock at very high $T$.
  - **Directed-loop SSE** [`syljuasen-sandvik-2002-directed-loops-qmc`]: similar performance for spin / hard-core boson systems; SSE works in a high-$T$ series basis rather than continuous imaginary time. Comparable per-update cost; SSE is more natural for spins, worm is more natural for soft-core bosons (no truncation in particle number per site).
  - **Hirsch–Fye / discrete-time world-line**: superseded in continuous time — no $\Delta\tau$ extrapolation, no Trotter error.
- **Gap vs. quantum Gibbs sampling**: in the regime where worm runs (sign-free bosons), there is **no plausible asymptotic quantum advantage** in the worst-case-polynomial sense:
  - Worm runtime is polynomial in $(L^d, \beta, 1/\varepsilon)$ with small constants and a low dynamic exponent ($z \approx 1$ deep in a phase, $z \approx 2$ at criticality).
  - Quantum KMS Lindbladian samplers run in $\mathrm{poly}(\dim\mathcal{H}, \beta, 1/\varepsilon)$ with $\dim\mathcal{H} = (n_{\max}+1)^{L^d}$ for soft-core bosons — exponentially worse in space cost than worm.
  - The relevant comparison is therefore not "asymptotic speedup" but "constant-factor on quantum hardware that natively encodes $\dim\mathcal{H}$ qubits".
- **Source of the difference (why worm exists at all here)**: the Bose–Hubbard kinetic vertex $-t\,b_i^\dagger b_j$ has *positive* off-diagonal matrix elements in the occupation-number basis, so the world-line weights are non-negative — the integrand of the path integral is a probability measure. This is the bosonic / sign-free / *stoquastic-in-the-occupation-basis* structural fact. The same algorithm fails immediately the moment this fails (see below).
- **Caveat**: like CT-QMC [`gull-et-al-2011-ct-qmc-impurity-models`], the polynomial guarantee is empirical / model-class specific. Troyer–Wiese-style worst-case hardness is sidestepped by *staying inside the sign-free corner*, not by curing anything.

## Implications for Quantum Advantage

- **Regime cell**: noncommuting × all-T × geometric-local-lattice × bosonic, sign-free.
- **What this changes**: closes the cell as a candidate for early quantum advantage. Worm runs in low-degree polynomial time on every problem in this cell; quantum samplers cannot beat it on these models without a constant-factor argument that already presupposes scalable fault-tolerant quantum hardware. The chapter should treat the **bosonic sign-free corner as a solved baseline** that quantum advantage is unlikely to penetrate.
- **Where the wall is — failure modes that *do* open quantum advantage**:
  1. **Fermions**: worm's locality is built on positive world-line weights. Generic itinerant fermions have $(-1)^P$ permutation signs from world-line crossings $\Rightarrow$ standard sign problem; worm does not improve over determinantal QMC here. Fermionic worm variants exist but inherit the sign problem.
  2. **Frustrated bosonic / spin systems**: if the off-diagonal Hamiltonian has *negative* matrix elements in any local basis (e.g. antiferromagnetic $XY$ on a triangular lattice, ring-exchange terms), world-line weights become signed and worm degrades exponentially in $\beta$.
  3. **Strong long-range interactions / dipolar gases**: worm still samples, but autocorrelation grows polynomially with the interaction range and the algorithm loses its "local-only" advantage.
  4. **Real-time Keldysh dynamics**: world-line weights become oscillatory, sign $\sim e^{-ct}$.
- **Promising or not (for quantum)**: not in this cell, but the *complement* is the actual frontier — frustrated and fermionic regimes where worm fails. The Pollet (2012) review cites worm + diagrammatic MC as the dominant tool for cold-atom Bose–Hubbard simulations and explicitly identifies the boundary at frustration / finite-density fermions, mirroring the quantum-advantage frontier.
- **Practical role in cold-atom physics**: worm-algorithm Bose–Hubbard simulations are the *gold standard for benchmarking ultracold-atom optical-lattice experiments* (Greiner et al. superfluid–Mott transition, single-site-resolution microscopes, disordered "dirty boson" experiments). Any quantum Gibbs-sampling claim in the bosonic regime must beat worm on the *exact same observables* — $G(i-j,\tau)$, $\rho_s(\beta)$, density profiles in a trap. This sets the experimental yardstick for the row.

## Open Questions / Limitations

- **Critical slowing at the Mott transition**: dynamic exponent $z \approx 2$ near $(t/U)_c$ even with worm; this is the critical-slowing regime where critical exponents demand large $L$ and the simulation becomes expensive. No theoretical reason to expect quantum samplers to do *better* in this regime (KMS gap also closes), but it is the only sub-cell where the constant prefactor matters.
- **Sign-problem-free *only*** — quantitative measure of "how non-stoquastic" before worm collapses is empirical (cf. Hangleiter–Roth–Nagaj–Eisert on sign-easing). No rigorous polynomial mixing-time bound for worm exists in the literature analogous to the spectral-gap / log-Sobolev framework for classical Glauber.
- **Fermionic worm variants**: bold-diagrammatic Monte Carlo (Prokof'ev–Svistunov 2007 PRL, in `supplementary-informations/classical-review/` if present) extends the worm philosophy to skeleton self-energy diagrams for fermions, but convergence relies on (re)summed series and is empirical, not provable.
- **Autocorrelation in disorder**: Letter's 1D disordered demonstration is heuristic; dependence on disorder strength / box size is established empirically, no rigorous bound.
- [CHECK] The exact dynamic exponent at the Bose–Hubbard tip is $z = 1$ in the literature (mean-field-like at the multicritical tip) but $z = 2$ along the generic phase boundary; cite Pollet 2012 / Capogrosso-Sansone–Söyler–Prokof'ev–Svistunov when the chapter quotes a number.

## Connections

- **Sister classical method**: Syljuåsen–Sandvik directed-loop SSE [`syljuasen-sandvik-2002-directed-loops-qmc`] — same regime cell (sign-free quantum lattices), different representation (high-$T$ series vs. continuous imaginary time). Worm dominates for soft-core bosons / off-diagonal observables; SSE dominates for spin systems with no extended occupation basis. Both are sign-locked, and both fail in the same frustration / fermion regimes.
- **Fermionic continuous-time descendant**: Gull et al. (2011) [`gull-et-al-2011-ct-qmc-impurity-models`] — explicitly cites the worm algorithm as conceptual ancestor (RMP §II.A); CT-QMC ports the diagrammatic-sampling-without-Trotter-error idea from bosonic worm to fermionic impurity models, paying the sign-problem price away from particle-hole symmetry.
- **Provable polynomial classical baseline for stoquastic ferromagnets**: Bravyi–Gosset (2017) [`bravyi-gosset-2017-quantum-ferromagnets`] — proves polynomial-time *deterministic* (FPRAS) approximation for partition functions of XY ferromagnets via Suzuki–Trotter to weighted matchings, on a strict subset of the sign-free corner where worm is the empirical practitioner's tool. Worm has no proof of polynomial mixing on the same models but works on a much wider class (general $U(1)$-symmetric bosons, disordered, off-diagonal observables, finite $T$).
- **Continuous-space extension**: Boninsegni–Prokof'ev–Svistunov (2006), continuous-space worm for ${}^4$He / cold-atom path-integral Monte Carlo — same algorithm, no lattice. Establishes that the worm idea works in the limit $a\to 0$ and is the standard for liquid-helium / dilute-Bose-gas first-principles simulations.
- **Bold-diagrammatic / fermion frontier**: Prokof'ev–Svistunov (2007) PRL — same authors, fermionic skeleton-diagram extension; weaker convergence guarantees, but the closest classical attempt at the regime where quantum samplers might first win.
- **Thesis Review chapter placement**: cite as the *bosonic-row classical baseline* paired with Syljuåsen–Sandvik (spin row) and Gull et al. (fermionic impurity row). Three together define the classical sign-free frontier; the row of frustrated / itinerant fermions / non-stoquastic systems beneath them is the quantum-advantage target.

Sources:
- [Worm algorithm in QMC — ScienceDirect (Phys. Lett. A 238:253)](https://www.sciencedirect.com/science/article/abs/pii/S0375960197009572)
- [Worm algorithm in QMC — ADS abstract](https://ui.adsabs.harvard.edu/abs/1998PhLA..238..253P/abstract)
- [Pollet 2012 review — cold-atom QMC](https://pubmed.ncbi.nlm.nih.gov/22885729/)
- [Svistunov lecture notes on worm and diagrammatic MC (QCS 2019)](https://indico.fysik.su.se/event/6664/attachments/4146/4770/QCS2019_Lecture_notes_Svistunov.pdf)
- [Pollet–Prokof'ev–Svistunov disordered Bose–Hubbard (PRL 92:015703)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.92.015703)
