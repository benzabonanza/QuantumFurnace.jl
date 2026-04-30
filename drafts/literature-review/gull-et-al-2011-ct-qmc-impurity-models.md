---
paper: "Gull, Millis, Lichtenstein, Rubtsov, Troyer, Werner (2011)"
title: "Continuous-time Monte Carlo methods for quantum impurity models"
arxiv: "1012.4474"
year: 2011
venue: "Rev. Mod. Phys. 83:349"
pdf: "supplementary-informations/classical-review/gull-et-al-2011-ct-qmc-impurity-models.pdf"

temperature: [intermediate-T, low-T]
commutativity: [noncommuting]
locality: [sparse, k-local]
particle-statistics: [fermionic, bosonic]
hamiltonian-models: [Anderson-impurity, Hubbard, Kondo-lattice, multi-orbital-impurity]

paradigm: [classical-other]
quantum-or-classical: [classical]

result-type: [runtime-comparison, complexity-upper, structural]
key-scaling: "per-sweep cost $\\sim O(\\langle k\\rangle^3)$ where the mean expansion order $\\langle k\\rangle$ is the matrix size; for CT-INT/CT-AUX $\\langle k\\rangle \\sim N\\beta U$ (linear in inverse temperature, system size, and interaction); for CT-HYB segment $\\langle k\\rangle \\sim N\\beta$; sign-problem reweighting cost $\\sim e^{2\\beta\\Delta F}/M$ in the bad regimes"

related: ["syljuasen-sandvik-2002-directed-loops-qmc", "prokofev-svistunov-tupitsyn-1998-worm-algorithm", "troyer-wiese-2005-sign-problem-hardness", "bravyi-divincenzo-oliveira-terhal-2008-stoquastic-lh", "bravyi-terhal-2010-frustration-free-stoquastic", "hangleiter-roth-nagaj-eisert-2020-easing-sign-problem", "klassen-marvian-2020-curing-sign-problem"]
---

# Continuous-time Monte Carlo methods for quantum impurity models — Gull et al. (2011)

**One-sentence takeaway**: CT-QMC (CT-INT, CT-AUX, CT-HYB) is the dominant classical finite-$T$ solver for quantum impurity models — provably polynomial in $\beta$ and system size on sign-problem-free corners (e.g. particle-hole-symmetric Hubbard, segment-representable density-density CT-HYB), but inherits the generic fermionic sign problem, with cost $\sim e^{2\beta\Delta F}$ once particle-hole symmetry is broken or many orbitals couple non-diagonally.

## Setting

- **Object sampled**: thermal expectation values $\langle A\rangle = \tfrac{1}{Z}\,\mathrm{tr}\,e^{-\beta H_\mathrm{QI}} A$ for a *quantum impurity* Hamiltonian $H_\mathrm{QI} = H_\mathrm{loc} + H_\mathrm{bath} + H_\mathrm{hyb}$ — finite (typically $\le 5$ orbitals) interacting impurity coupled to a non-interacting fermionic bath of infinitely many modes (Eq. 1, p.3).
- **Model class**: single- and multi-orbital Anderson, Kondo lattice, Hubbard (via DMFT and quantum-cluster mappings to an impurity), Slater–Kanamori multiplet interactions (Eq. 13). Bosonic (phonon) baths via CT-HYB extensions (§VII).
- **Algorithmic family**: continuous-imaginary-time diagrammatic Monte Carlo, sampling the perturbation expansion $Z = \sum_k (-1)^k \int_0^\beta\!\!\!\cdots\!\int d\tau_1\!\cdots d\tau_k\,\mathrm{tr}[e^{-\beta H_a} H_b(\tau_k)\cdots H_b(\tau_1)]$ (Eq. 19) by Metropolis updates that insert/remove vertices at random imaginary times. Three flavours by choice of $H_a/H_b$ split: **CT-INT** (expand in interaction $H^I_\mathrm{loc}$), **CT-AUX** (CT-INT after a Hubbard–Stratonovich auxiliary-field shift), **CT-HYB** (expand in hybridization $H_\mathrm{hyb}$); plus **CT-J** for Kondo-like exchange (§VI).
- **Key feature**: *no Trotter discretization error* — the "$1/k!$ from the exponential" makes infinite-order summation finite, so the algorithm has only statistical error. Distinguishes CT-QMC from Hirsch–Fye and any path-integral method that discretises $[0,\beta]$ into $M = \beta/\Delta\tau$ slices.
- **Temperature regime**: any $T > 0$ in principle; in practice limited by sign problem and by mean expansion order $\langle k\rangle \propto \beta$.

## Main Results / Capabilities

- **Polynomial-scaling (sign-free) corners** — provable cost cubic in matrix size, polynomial in $(N, \beta)$:
  - **CT-HYB segment** (one-orbital Anderson, density-density multi-orbital with each orbital hybridising to a separate bath): cost $\sim N\beta^3$ per sweep, no sign problem (Table I, p.56). Yoo et al. (2005) *prove* segment Hirsch–Fye / CT-HYB on the Anderson model is sign-free; the proof extends to single-particle-basis-diagonal multi-orbital cases.
  - **CT-INT / CT-AUX** at half filling on bipartite single-orbital Hubbard / particle-hole-symmetric Anderson: bipartite structure makes odd-order diagrams cancel pairwise (CT-AUX trick $K_x = -xU/4$ fixes only even orders, p.29). Cost $\sim (N\beta U)^3$.
- **Where they fail** (sign problem on, exponential cost):
  - CT-INT/CT-AUX with general $U_{ijkl}$ interaction $\Rightarrow$ $(N^2\beta U)^3$ + sign problem (Table I).
  - CT-HYB **matrix** (general non-density-density local interaction such as spin-flip / pair-hopping in Slater–Kanamori): cost $a e^{N}\beta^2 + bN\beta^3$ — *exponential in number of orbitals* $N$, the local-trace evaluation is over the full $4^N$ Hilbert space (5 orbitals at 2024-era HPC limits, p.55). Sign problem present.
  - Cluster DMFT (multi-site impurity): sign problem reappears from inter-site fermionic loops (Hubbard away from half filling, frustrated lattices). "The reduction in severity of the sign problem [vs lattice QMC] has allowed the development of flexible and powerful continuous-time QMC impurity solvers" — but cluster size is the binding constraint, not impurity Hilbert space.
  - Real-time / nonequilibrium (Keldysh): sign grows $\sim e^{ct}$ with simulation time on the real branch (§VIII.D); average sign $10^{-3}$ achievable for $\sim 10$ real-branch operators only.
- **Quantitative comparison vs Hirsch–Fye** (Fig. 18, p.30, single-site DMFT Hubbard at $U/t=4$, $\beta t = 30$): matrix size $\sim 75$ for HF vs $\sim 35$ (CT-INT) and $\sim 10$ (CT-HYB). Cubic scaling means CT-HYB is $\sim 400\times$ faster *and* eliminates the $\Delta\tau \to 0$ extrapolation. CT-HYB matrix size *decreases* with $U/t$ in the strongly correlated / Mott regime, opposite to the weak-coupling expansions — making CT-HYB the unique tractable solver for $U/t \gtrsim 5$.

## Method (1–3 sentences)

Each algorithm samples the same continuous-time expansion (Eq. 19) of $\mathrm{tr}\,e^{-\beta H}$ but differs in the operator chosen as the perturbation: CT-INT/AUX expand in the local interaction (matrix size $\propto U$), CT-HYB expands in the impurity-bath hybridization (matrix size independent of $U$). Updates are Metropolis insertions/removals of vertices at uniformly drawn imaginary times $\tau \in [0,\beta]$; the diagrammatic weight is a determinant of an $\langle k\rangle\times\langle k\rangle$ matrix of bath/hybridization Green's functions (Wick's theorem on $H_\mathrm{bath}$), updated incrementally in $O(k^2)$ time. Continuous time eliminates Suzuki–Trotter error and produces an *adaptive* time grid concentrated where Green's functions vary fastest.

## Quantum vs Classical

- **Role in the corpus**: this paper defines the *classical upper frontier* for the noncommuting × intermediate-T × sparse/impurity cell. Anywhere DMFT or quantum-cluster maps a strongly-correlated lattice problem to an impurity, CT-QMC is the gold-standard solver — Hubbard at $U/t \in [4,7]$, Mott metal–insulator, heavy fermion, multi-orbital correlated materials. Quantum Gibbs samplers must beat CT-QMC's *cubic-in-matrix-size* runtime (where matrix size is $\sim N\beta U$ or $\sim N\beta$) on the same DMFT/quantum-cluster benchmarks, *or* operate in CT-QMC's sign-broken regimes (cluster DMFT off half-filling, multi-orbital with non-diagonal interactions, real-time Keldysh).
- **Where quantum could plausibly beat classical CT-QMC**:
  1. **Multi-orbital matrix CT-HYB**: cost $\sim ae^{N}\beta^2$ scales exponentially in orbitals because the local trace is taken over the full $4^N$ impurity Hilbert space. A KMS-Lindbladian Gibbs sampler running directly on a fault-tolerant quantum simulator of the impurity has $\mathrm{poly}(N)$ Hilbert-space cost — the natural regime for early advantage on transition-metal-oxide multiplet physics.
  2. **Cluster DMFT off half filling / frustrated lattices**: classical sign $\sim e^{-2\beta\Delta F}$ kills runtime exponentially in $\beta$ and cluster size. Quantum Gibbs sampling is sign-blind by construction.
  3. **Nonequilibrium / Keldysh**: real-time CT-QMC has $\langle\mathrm{sign}\rangle \sim e^{-\Gamma t}$, capping $t$ at $\sim 10/\Gamma$. Quantum simulation has no sign problem in real time.
- **Where classical CT-QMC will not be beaten**: single-orbital Anderson, density-density multi-orbital with segment representation, bipartite half-filled Hubbard. These are the "easy" sign-free corners; CT-HYB segment runs in linear-in-$N$, cubic-in-$\beta$ time and the constant prefactor is small enough that any near-term quantum solver will be slower even on tractable instances.
- **Caveat — Troyer–Wiese caveat applies**: the CT-QMC sign problem is the same NP-hard worst-case sign problem [Troyer, Wiese 2005], but its *severity* depends on the basis and partitioning ($H_a$ vs $H_b$). The authors stress (p.7, p.55) that impurity sign problems are typically *milder* than the lattice ones — bath integration "turns off" couplings that make lattice QMC exponentially slow. This is a specific, narrow exemption that classical practitioners exploit aggressively but does not generalise to lattice problems.

## Implications

- **Defines the classical upper-frontier for the noncommuting × intermediate/low-T × sparse cell of the thesis 3$\times$3 grid**: any quantum advantage claim on Hubbard / Anderson / Kondo / DMFT *must* benchmark against CT-QMC, not against Hirsch–Fye (obsolete) or ED (only for $\le 15$ sites).
- **Sharpens "sparse"**: CT-QMC's polynomial regime is fundamentally about *small Hilbert-space impurities coupled to non-interacting baths*. The bath is integrated out analytically (Wick's theorem $\Rightarrow$ determinantal weight), leaving a sampling cost set only by the impurity. Quantum Gibbs samplers should target problems where this analytical bath-integration breaks down — multi-orbital Slater–Kanamori with $N \ge 6$, cluster DMFT with $\ge 8$ sites, finite-density doped Hubbard, and real-time correlators.
- **Sign problem severity is the deciding metric**: Table I gives the precise classification — diagonal hybridisation + density-density gives $N\beta^3$ sign-free; off-diagonal hybridisation gives $(N\beta)^3$ + sign problem; general interaction gives $(N^2\beta U)^3$ + sign problem; matrix CT-HYB scales $e^N\beta^2$ for arbitrary local interactions. Each row of this table is a candidate quantum-advantage cell.

## Open Questions / Limitations

- **Truncation accuracy of large-$N$ Hilbert spaces**: §X.F flags this as an open problem — block diagonalisation by symmetries works in special cases (e.g. Coulomb + spin–orbit) but the accuracy of generic Hilbert-space truncation is "not yet established". Quantum Gibbs sampling on the full $4^N$ Hilbert space provides a clean alternative but the comparison is not quantified.
- **Single-particle-basis assumption**: the entire CT-QMC framework assumes the bath couples through a hybridization $\Delta(\tau)$ that can be written in a single-particle basis. Strongly correlated baths (e.g. impurity-in-Mott-insulator) violate this and require self-consistent (DMFT) iteration; convergence to the self-consistent fixed point has empirical guarantees only.
- **Real-time / Keldysh**: severe sign problem caps simulation times at $t \sim 10\beta$. The authors consider this the most important open quantitative challenge; matches the Review chapter's Keldysh sub-cell.
- **"How severe is the sign problem really?"**: §IX gives runtime numbers but no rigorous bound on $\Delta F$ as a function of $(N, \beta, U)$ — practitioners use heuristics. A rigorous classification of sign-problem-free CT-QMC sub-classes (analogous to Klassen–Marvian's two-local sign-curing dichotomy) is missing.
- [CHECK] CT-J for Kondo-like exchange interactions is reviewed in §VI but is less mature than CT-INT/AUX/HYB; its applicability profile for $f$-electron / heavy-fermion physics is asserted but not benchmarked against KMS Lindbladians.

## Connections

- **Classical predecessors / sister methods**:
  - Hirsch–Fye (1986) — discrete-time predecessor; CT-QMC supersedes it on every metric (Fig. 18).
  - Sandvik–Kurkijärvi (1991) SSE — analogous high-temperature expansion for spin systems; same "sample the perturbation series, no Trotter error" principle. Syljuåsen–Sandvik (2002) [`syljuasen-sandvik-2002-directed-loops-qmc`] is the directed-loop SSE descendant for sign-free spin lattices.
  - Prokof'ev–Svistunov–Tupitsyn (1998) [`prokofev-svistunov-tupitsyn-1998-worm-algorithm`] — worm algorithm; cited as the conceptual foundation for diagrammatic sampling of *bosonic* lattice models. CT-QMC ports the same idea to fermionic impurities.
- **Sign-problem hardness backbone**:
  - Troyer–Wiese (2005) [`troyer-wiese-2005-sign-problem-hardness`] — NP-hardness of generic sign-curing applies to CT-QMC's bad regimes; explicitly cited (p.7) as the structural reason "sign problems are physical and unavoidable, at least in itinerant phases with unpaired fermions".
  - Bravyi–DiVincenzo–Oliveira–Terhal (2008) [`bravyi-divincenzo-oliveira-terhal-2008-stoquastic-lh`] — defines stoquastic Hamiltonians as the sign-free ones in the standard basis; CT-QMC's sign-free corners (segment CT-HYB, particle-hole-symmetric CT-AUX) are *fermionic* analogues that lie outside the qubit-stoquastic dichotomy but share the structural origin (real, non-positive off-diagonal weights in the right basis).
  - Hangleiter–Roth–Nagaj–Eisert (2020) [`hangleiter-roth-nagaj-eisert-2020-easing-sign-problem`] and Klassen–Marvian et al. (2020) [`klassen-marvian-2020-curing-sign-problem`] — modern quantification / hardness of basis-change sign-curing; complementary lens on CT-QMC's particle-hole-symmetric sign-free trick.
- **Stoquastic frustration-free boundary**: Bravyi–Terhal (2010) [`bravyi-terhal-2010-frustration-free-stoquastic`] gives a *provable* polynomial classical simulation for a different "easy" quantum class; CT-QMC is the empirical workhorse on the much broader (non-frustration-free, non-stoquastic-in-qubit-basis) impurity class, with rigour traded for applicability.
- **Quantum-side counterparts** (out-of-corpus but relevant): Chen–Kastoryano–Gilyén KMS Lindbladians, Ding et al. KMS detailed-balance samplers — these run in $\mathrm{poly}(\dim\mathcal{H})$ and so beat matrix CT-HYB asymptotically in $N$; the open question is the *crossover* $N$ where a fault-tolerant quantum sampler becomes practical.
