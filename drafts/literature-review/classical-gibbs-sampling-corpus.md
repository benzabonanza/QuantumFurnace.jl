# Classical Gibbs Sampling — Reference Corpus

This corpus collects the classical-side baseline for the thesis Review chapter on quantum vs classical Gibbs sampling. By "classical Gibbs sampler" we mean both (i) classical algorithms that sample from $p \propto e^{-\beta H}$ for a classical Hamiltonian $H$ on a discrete or continuous configuration space, and (ii) classical algorithms applied to *quantum* Hamiltonians via path-integral / world-line / auxiliary-field quantum Monte Carlo or tensor-network thermal techniques. Quantum-side references (Chen et al. KMS Lindbladians, Gilyén-Su-Low-Wiebe, Rall, Kastoryano, Brandão, etc.) are tracked separately. Within each section the entries are roughly chronological so that classical results, hardness barriers, and modern sharpenings appear in context. References are limited to two sentences each — substantive review of any single paper is left to a later pass.

## Review priority — Tier-1 grid

The chapter is organized as a 3×3 regime grid: **rows = Hamiltonian structure**, **columns = temperature regime**. Modifiers (dimension, system size, particle statistics) live *inside* cells, not as top-level axes. The framing is *current best classical method vs current best quantum method per cell*, not a historic overview — older results are tagged Tier-2 (cite-only) when superseded by sharper modern statements.

|                                 | **High-T** (above transition)                                                                       | **Intermediate / near critical**                                                                    | **Low-T** (below transition)                                                                              |
| ------------------------------- | --------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------- |
| **Classical** Ising / Potts     | Eldan–Koehler–Zeitouni; Chen–Liu–Vigoda; Anari–Liu–Oveis Gharan; Lubetzky–Sly (cutoff, $d$-dep.)    | Bauerschmidt–Dagallier; Lubetzky–Sly (crit. 2D); Guo–Jerrum; Weitz; Cuff et al. (mean-field Potts)  | Sly (hardness); Bovier–den Hollander; Woodard–Schmidler–Huber (PT no-go); Madras–Zheng (PT positive)      |
| **Stoquastic** quantum $H$      | Mann–Helmuth (high-T); Bakshi–Liu–Moitra–Tang; Kuwahara–Kato–Brandão                                | Bravyi–Gosset (XY ferromagnet); Bravyi–DiVincenzo–Oliveira–Terhal (def.)                            | Bravyi–Terhal (frust.-free + gap); Gull et al. (CT-QMC); Mann–Helmuth (low-T polymer)                     |
| **Non-stoquastic** quantum $H$  | (high-T row applies — sign problem mild)                                                            | Klassen–Marvian (sign-cure NP-hard); Hangleiter–Roth–Nagaj–Eisert (sign-easing framework)           | Troyer–Wiese (sign-problem worst-case hard); Klassen–Marvian; Hangleiter et al.                           |

**Cross-cutting Tier-1 anchors** (dimension- or method-defined, span multiple cells):

- **Cluster-algorithm structure:** Edwards–Sokal (random-cluster joint representation — defines exactly when a cluster algorithm exists).
- **1D quantum:** White (METTS); Kuwahara–Alhambra–Anshu (rigorous quasi-linear 1D Gibbs).
- **2D / higher-$d$ quantum:** Czarnik–Dziarmaga (PEPS imaginary-time); Molnar–Schuch–Verstraete–Cirac (PEPO any-$d$, $\exp(O(\beta))$ bond dimension).
- **Sign-free QMC workhorses** (stoquastic / bipartite half-filled): Syljuåsen–Sandvik (directed-loop SSE); Prokof'ev–Svistunov (worm).
- **Modern mixing-time framework:** Chen–Eldan (localization schemes — unifies discrete spectral independence with continuous stochastic localization).

**Tier-1 count: 40 papers** (deep `/lit-review`). **Tier-2: ~75** (cited from the corpus as supporting context, no individual review). Each Tier-1 entry below is tagged `**[T1]**` at the start of its title line; absence of a tag means Tier-2.

**Where quantum could first beat classical (the chapter's punchline cell)**: the **non-stoquastic × low-T** corner is the strongest candidate — classical sign-free QMC fails (Troyer–Wiese), sign-curing is NP-hard (Klassen–Marvian), tensor networks break down past 1D, and Chen / Ding-style KMS Lindbladians give honest poly-time bounds inside this cell. The **non-stoquastic × intermediate-T** cell is the next frontier for the same structural reasons.

## 1. Foundations of classical MCMC mixing theory

### Levin, Peres, Wilmer — *Markov Chains and Mixing Times* (2nd ed., 2017) [AMS]
Canonical textbook covering coupling, strong stationary times, spectral methods, log-Sobolev inequalities, and the cutoff phenomenon, with chapters on lattice spin systems, exclusion processes, and shuffling. Gives the standard $O(n \log n)$ vs $\exp(\Omega(n))$ dichotomy framework that frames every later result and is the single most-cited reference for finite-state mixing.

### Diaconis, Saloff-Coste — *Logarithmic Sobolev inequalities for finite Markov chains* (Ann. Appl. Probab. 6:695–750, 1996)
Introduces log-Sobolev constants $\alpha$ for finite reversible Markov chains and proves $t_{\mathrm{mix}}(\varepsilon) \le \alpha^{-1}(\log\log(1/\pi_*) + \log(1/\varepsilon))$, sharper than the spectral-gap bound when $\pi_*$ is exponentially small. Holds for any reversible chain; the gap-vs-log-Sobolev ratio can be exponential, so this controls rapid mixing on huge state spaces (Ising, hardcore, etc.) where the eigenvalue bound alone is too weak.

### Diaconis, Saloff-Coste — *Nash inequalities for finite Markov chains* (J. Theor. Probab. 9:459–510, 1996)
Develops Nash inequalities as a complement to spectral and log-Sobolev methods, giving polynomial decay rates and sharper short-time bounds for random walks on groups. Applies to any finite Markov chain; Nash bounds dominate when convergence is governed by escape from low-density regions rather than a single bottleneck.

### Aldous, Diaconis — *Strong uniform times and finite random walks* (Adv. Appl. Math. 8:69–97, 1987)
Foundational paper introducing strong stationary times and using them to bound total-variation mixing of card shuffles. The methodology underlies the cutoff phenomenon and is a standard tool whenever the chain has a natural "filling-up" structure (random transpositions, top-to-random, etc.).

### Lyons, Peres — *Probability on Trees and Networks* (Cambridge, 2016)
Comprehensive treatment of random walks on graphs and groups, electrical-network methods, percolation, and isoperimetric inequalities, all of which feed directly into Gibbs-sampler analysis. Covers the regime where geometric structure of the underlying graph determines mixing behaviour, including Galton–Watson trees and transitive graphs.

### Sokal — *Monte Carlo Methods in Statistical Mechanics: Foundations and New Algorithms* (Cargèse Lectures, 1996)
Pedagogical lecture notes laying out the conceptual foundations of MCMC for lattice spin systems, autocorrelation analysis, and dynamic critical exponents $z$ for Metropolis, heat-bath, Swendsen–Wang, and Wolff dynamics. Standard reference for practitioners; pinpoints critical slowing-down ($z \approx 2$ for local algorithms) as the universal obstruction near $T_c$.

### Krauth — *Statistical Mechanics: Algorithms and Computations* (Oxford Master Series, 2006)
Textbook covering Metropolis, heat-bath, cluster, and event-chain algorithms for hard spheres, spins, bosons, and path integrals, with explicit pseudocode. Ranges from finite-T classical statistical mechanics to bosonic path-integral Monte Carlo, providing the practitioner-level baseline for what "classical sampling" actually means in physics simulations.

### **[T1]** Suwa, Todo — *Markov chain Monte Carlo method without detailed balance* (PRL 105:120603, 2010) [arXiv:1007.2262]
Constructs irreversible transition kernels via geometric weight allocation that minimize (often eliminate) rejection while preserving the Gibbs measure. Universal over any finite-state target with the multiplicative structure of an exponentiated Hamiltonian; gives a constant-factor speedup over Metropolis on Potts and similar models, with no impact on the asymptotic dynamic exponent at criticality. Also a cornerstone of §12 (non-reversible / lifted MCMC) — listed here for the lattice spin-system context.

## 2. Glauber / Metropolis on lattice spin systems — log-Sobolev framework

### Holley, Stroock — *Logarithmic Sobolev inequalities and stochastic Ising models* (J. Stat. Phys. 46:1159–1194, 1987)
Proves that a uniform log-Sobolev inequality for the Glauber dynamics implies exponential approach to the Gibbs measure and a complete large-deviation principle for the stochastic Ising model. Covers compact-spin Ising models above the uniqueness threshold; in 1D the convergence is $\exp(-\gamma t / \log t)$, which is the prototype slow-but-rapid behaviour at the uniqueness boundary.

### Stroock, Zegarlinski — *The equivalence of the logarithmic Sobolev inequality and the Dobrushin–Shlosman mixing condition* (Comm. Math. Phys. 144:303–323, 1992)
Proves equivalence between the Dobrushin–Shlosman complete-analyticity / strong spatial-mixing condition and a uniform log-Sobolev inequality for the associated Glauber dynamics on finite-range lattice systems. The equivalence holds for finite spin spaces; this is the canonical "spatial mixing $\Leftrightarrow$ temporal mixing" reduction in the Gibbs-sampling literature.

### Stroock, Zegarlinski — *The logarithmic Sobolev inequality for discrete spin systems on a lattice* (Comm. Math. Phys. 149:175–193, 1992)
Extends the equivalence above to discrete-spin lattice systems and gives the uniform-in-volume log-Sobolev inequality under Dobrushin–Shlosman mixing. Applies in the entire one-phase region for any finite-range interaction; breaks down at coexistence where multiple phases destroy the mixing condition.

### Martinelli, Olivieri — *Approach to equilibrium of Glauber dynamics in the one phase region I, II* (Comm. Math. Phys. 161:447–514, 1994)
Renormalization-group / block-decimation proof that strong spatial mixing implies a $O(\log n)$ log-Sobolev constant for Glauber dynamics on the discrete lattice. Holds throughout the one-phase region (any temperature with unique infinite-volume Gibbs state); fails at and below criticality.

### Martinelli, Olivieri, Schonmann — *For 2-D lattice spin systems weak mixing implies strong mixing* (Comm. Math. Phys. 165:33–47, 1994)
Closes the gap in 2D: weak spatial mixing alone implies the stronger Dobrushin–Shlosman condition for finite-range discrete spins on $\mathbb{Z}^2$. Specific to 2D lattices; combined with Stroock–Zegarlinski this gives uniform log-Sobolev for the 2D Ising model arbitrarily close to (but not at) $T_c$.

### Bodineau, Helffer — *The log-Sobolev inequality for unbounded spins systems* (J. Funct. Anal. 166:168–178, 1999)
Spectral / Witten-Laplacian proof of the log-Sobolev inequality for continuous-spin lattice systems with possibly higher-than-quadratic interactions in the perturbative regime. Holds at high temperature for unbounded continuous spins (e.g. $\phi^4$); does not extend to long-range or critical interactions.

### Bauerschmidt, Bodineau — *A very simple proof of the LSI for high temperature spin systems* (J. Funct. Anal. 276:2582–2588, 2019) [arXiv:1712.10039]
Two-page proof of the log-Sobolev inequality at high temperature using a multiscale Bakry–Émery criterion, dispensing with the renormalization machinery of Stroock–Zegarlinski. High-temperature regime only ($\beta$ less than a model-dependent threshold); the multiscale-Bakry-Émery framework has since been extended to $\phi^4$ and near-critical Ising.

### **[T1]** Eldan, Koehler, Zeitouni — *A spectral condition for spectral gap: fast mixing in high-temperature Ising models* (Probab. Theory Rel. Fields 182:1035–1051, 2022) [arXiv:2007.08200]
Shows that an Ising model with general interaction matrix $J$ on the hypercube has Glauber Poincaré constant bounded by an explicit function of $\|J\|_{\mathrm{op}}$, via a stochastic-localization decomposition into rank-one quadratic measures. Holds whenever $\|J\|_{\mathrm{op}} < 1$ (covers and extends Dobrushin); Glauber mixes in $O(n \log n)$ in this regime.

### **[T1]** Bauerschmidt, Dagallier — *Log-Sobolev inequality for near critical Ising models* (Comm. Pure Appl. Math. 77:2568–2576, 2024) [arXiv:2202.02301]
Proves a log-Sobolev inequality for the ferromagnetic Ising model with constant uniform in the system size up to the critical point, with explicit dependence on the susceptibility. Closes the gap between perturbative high-temperature results and criticality for ferromagnets; does not handle the antiferromagnetic case or low-temperature regime.

### **[T1]** Lubetzky, Sly — *Cutoff for the Ising model on the lattice* (Invent. Math. 191:719–755, 2013) [arXiv:0909.4320]
Establishes total-variation cutoff at $(d/2\lambda_\infty) \log n$ for Glauber dynamics on the $d$-dimensional Ising model with periodic boundary conditions. Holds at any temperature in the high-temperature regime where $\lambda_\infty > 0$; the technique translates $L^1$ to $L^2$ convergence and combines with log-Sobolev.

### Lubetzky, Sly — *Information percolation and cutoff for the stochastic Ising model* (J. AMS 29:729–774, 2016) [arXiv:1401.6065]
Cutoff for stochastic Ising at all subcritical temperatures via a space-time information-percolation cluster expansion. Subcritical regime $\beta < \beta_c$; gives an $O(1)$-window around the natural mixing time and resolves cutoff far beyond the "very high temperature" reach of older log-Sobolev arguments.

### Dobrushin, Shlosman — *Completely analytical interactions: constructive description* (J. Stat. Phys. 46:983–1014, 1987)
Defines the family of equivalent conditions known as complete analyticity (CA) and gives a finite-volume constructive criterion to check it. Foundational for spatial-mixing analysis; CA is equivalent to log-Sobolev for finite-range lattice systems and implies uniqueness of the Gibbs measure plus exponential decay of truncated correlations.

## 3. Critical slowing down and cluster algorithms

### Swendsen, Wang — *Nonuniversal critical dynamics in Monte Carlo simulations* (PRL 58:86–88, 1987)
Introduces the global cluster update for $q$-state Potts via the Fortuin–Kasteleyn random-cluster representation, drastically reducing the dynamic critical exponent $z$ at $T_c$. Empirically $z \approx 0.25$ for 2D Ising vs $z \approx 2$ for Metropolis; works whenever the Fortuin–Kasteleyn representation has a non-negative measure (ferromagnetic Potts/Ising on any graph).

### Wolff — *Collective Monte Carlo updating for spin systems* (PRL 62:361–364, 1989)
Single-cluster variant of Swendsen–Wang in which only one Fortuin–Kasteleyn cluster is grown and flipped per step, generalizable to $O(n)$ models via an embedded-Ising representation. Best-in-class for ferromagnetic Ising/Potts and $O(n)$ models near criticality; like Swendsen–Wang it does not extend cleanly to frustrated or disordered systems where the FK measure has signs.

### Niedermayer — *General cluster updating method for Monte Carlo simulations* (PRL 61:2026–2029, 1988)
Generalizes Swendsen–Wang to a one-parameter family of bond probabilities and arbitrary discrete or continuous global symmetries. Useful when the global symmetry is non-Ising (e.g. $\mathbb{Z}_3$, $O(n)$, sigma models); does not address critical slowing in models without such a symmetry.

### **[T1]** Edwards, Sokal — *Generalization of the Fortuin–Kasteleyn–Swendsen–Wang representation and Monte Carlo algorithm* (Phys. Rev. D 38:2009–2012, 1988)
Joint $(\sigma, n)$ representation that systematizes when a cluster algorithm exists for a given Hamiltonian. Provides the conceptual scaffold (auxiliary-field cluster representations) reused in later quantum cluster algorithms and identifies the obstruction (sign of the joint measure) for non-ferromagnetic systems.

### Gore, Jerrum — *The Swendsen–Wang process does not always mix rapidly* (J. Stat. Phys. 97:67–86, 1999)
Constructs Potts-on-complete-graph instances where Swendsen–Wang has exponentially slow mixing, contradicting the folklore that cluster algorithms always speed up over Metropolis. Shows that cluster speedup is not universal: at first-order transitions in the mean-field Potts model the global cluster move can become a metastable trap.

### Ullrich — *Comparison of Swendsen–Wang and heat-bath dynamics* (Random Struct. Alg. 42:520–535, 2013) [arXiv:1105.3665]
Proves a polynomial comparison between Swendsen–Wang and single-site heat-bath spectral gaps on arbitrary graphs of bounded degree. Holds for arbitrary Potts $q$ and arbitrary graphs; rapid mixing of one implies rapid mixing of the other up to polynomial factors.

### Ullrich — *Rapid mixing of Swendsen–Wang dynamics in two dimensions* (PhD thesis / arXiv:1212.4908, 2012)
Combines the comparison with Martinelli–Olivieri-type heat-bath bounds to prove rapid mixing of Swendsen–Wang on the 2D square lattice at all non-critical temperatures. First proof of rapid mixing of any classical chain for the 2D Ising model at all $T \ne T_c$; does not cover criticality itself.

### **[T1]** Guo, Jerrum — *Random cluster dynamics for the Ising model is rapidly mixing* (Ann. Appl. Probab. 28:1292–1313, 2018) [arXiv:1605.00139]
Shows that single-edge Glauber dynamics for the random-cluster model at $q=2$ has polynomial mixing time on any graph, hence so does Swendsen–Wang for the ferromagnetic Ising model. Universal in graph and temperature for $q=2$ ferromagnetic Ising; this is the strongest "rapid-mixing-everywhere" result for any classical Ising sampler.

### **[T1]** Lubetzky, Sly — *Critical Ising on the square lattice mixes in polynomial time* (Comm. Math. Phys. 313:815–836, 2012) [arXiv:1001.1613]
Polynomial mixing of single-site Glauber for the 2D Ising model at $\beta = \beta_c$. Specifically $O(n^{c})$ for some explicit $c$; shows that critical 2D Ising is not torpid, in contrast to the worst-case mean-field picture.

## 4. Spin glasses, low-T metastability, parallel and simulated tempering

### Geyer — *Markov chain Monte Carlo maximum likelihood* (Computing Science and Statistics 23:156–163, 1991)
Introduces parallel tempering / Metropolis-coupled MCMC with full-configuration swaps between replicas at different temperatures. Generic recipe for any multimodal target; performance depends entirely on the temperature ladder spacing, with no a priori guarantees.

### Marinari, Parisi — *Simulated tempering: a new Monte Carlo scheme* (Europhys. Lett. 19:451–458, 1992)
Introduces simulated tempering, where temperature itself is part of the state space and a single replica diffuses through temperatures. Foundational for spin-glass simulation; rigorous mixing analysis only available for specific models, and no general guarantee of polynomial speedup.

### Madras, Randall — *Markov chain decomposition for convergence rate analysis* (Ann. Appl. Probab. 12:581–606, 2002)
Decomposition technique for spectral gaps of Markov chains in terms of restricted gaps and a projection chain on the partition. Generic technical tool used in essentially every torpid/rapid mixing analysis of tempering; does not by itself give polynomial bounds.

### Bhatnagar, Randall — *Torpid mixing of simulated tempering on the Potts model* (SODA 2004)
Proves exponentially slow mixing of simulated tempering for the ferromagnetic mean-field Potts model with $q \ge 3$ regardless of the temperature schedule. Demonstrates that simulated tempering does not in general resolve first-order transitions; the bottleneck is the entropic gap between the disordered and $q$ ordered phases.

### **[T1]** Woodard, Schmidler, Huber — *Sufficient conditions for torpid mixing of parallel and simulated tempering* (Electron. J. Probab. 14:780–804, 2009)
Gives general sufficient conditions ("persistence" of a narrow mode) under which both parallel and simulated tempering mix exponentially slowly. Applies to mean-field Potts $q\ge 3$, mean-field Ising with bad ladders, and Gaussian mixtures with sharp peaks; together with Bhatnagar–Randall this is the canonical no-go for tempering.

### Woodard, Schmidler, Huber — *Conditions for rapid mixing of parallel and simulated tempering on multimodal distributions* (Ann. Appl. Probab. 19:617–640, 2009)
Companion paper giving sufficient conditions for *rapid* mixing of parallel and simulated tempering. Holds when the modes are quantitatively similar in width and the temperature ladder is geometric; the gap between this regime and the torpid regime is the practical operating envelope of replica exchange.

### **[T1]** Madras, Zheng — *On the swapping algorithm* (Random Struct. Alg. 22:66–97, 2003)
Polynomial mixing-time bound for parallel tempering on the mean-field Ising model with a sufficient number of geometric temperatures. Proof of concept that PT can give polynomial mixing on a model where local Glauber is exponentially slow; tight on $\beta$-spacing requirement.

### Jerrum — *Large cliques elude the Metropolis process* (Random Struct. Alg. 3:347–359, 1992)
Proves that the natural Metropolis chain for finding a planted clique has exponentially slow mixing in the relevant regime, an early canonical "torpid mixing" result. Shows that for planted-clique detection / spin-glass-style problems with overlap structure, Metropolis alone has provable exponential bottleneck.

### **[T1]** Bovier, den Hollander — *Metastability: A Potential-Theoretic Approach* (Springer, 2015)
Comprehensive monograph deriving sharp Eyring–Kramers asymptotics for low-temperature mixing using potential-theoretic / capacity techniques. Standard reference for low-$T$ Glauber on Ising-like models, where $t_{\mathrm{mix}} \sim \exp(\beta \Delta)$ with $\Delta$ the relevant energy barrier.

## 5. Hardness and no-go for classical Gibbs sampling

### **[T1]** Sly — *Computational transition at the uniqueness threshold* (FOCS 2010) [arXiv:1005.5584]
Proves that for the antiferromagnetic 2-spin (Ising / hardcore) model on bounded-degree graphs there is no FPRAS for the partition function at activities just past the tree-uniqueness threshold $\lambda_c(\mathbb{T}_\Delta)$, unless RP=NP. Establishes that the *computational* phase transition coincides with the *statistical-physics* phase transition on the regular tree — the canonical hardness-mirrors-physics result.

### Sly, Sun — *Counting in two-spin models on $d$-regular graphs* (Ann. Probab. 42:2383–2416, 2014) [arXiv:1203.2602]
Sharpens Sly (2010) to all of the tree non-uniqueness region for general antiferromagnetic 2-spin systems, including the Ising model with arbitrary external field. NP-hard to even approximate the partition function or sample from the Gibbs measure throughout the non-uniqueness region; together with Weitz it gives a sharp dichotomy.

### Galanis, Štefankovič, Vigoda — *Inapproximability for antiferromagnetic spin systems in the tree non-uniqueness region* (J. ACM 62:50, 2015) [arXiv:1305.2902]
Extends Sly–Sun-style hardness to general antiferromagnetic multi-spin systems, including $k$-colourings and antiferromagnetic Potts. Proves $\#$BIS-hardness throughout the tree non-uniqueness region for any antiferromagnetic model; matches (up to the threshold) the algorithmic results of Weitz / spectral independence.

### Goldberg, Jerrum — *Inapproximability of the Tutte polynomial* (Inform. and Comput. 206:908–929, 2008)
Establishes a fine-grained inapproximability landscape for the Tutte polynomial / Potts partition function on general graphs. Provides hardness for ferromagnetic Potts at low temperature on certain graphs, complementary to the Jerrum–Sinclair FPRAS for ferromagnetic Ising.

### Dyer, Frieze, Jerrum — *On counting independent sets in sparse graphs* (SIAM J. Comput. 31:1527–1541, 2002)
Shows that Glauber dynamics for hardcore independent sets has exponentially slow mixing on random $\Delta$-regular bipartite graphs for $\Delta \ge 6$, and that no FPRAS exists for $\Delta \ge 25$ unless RP=NP. Classical predecessor to Sly: pinpoints that single-site dynamics on the hardcore model fails well before any computational hardness bound and motivates the spatial-mixing programme.

### Jerrum, Sinclair — *Polynomial-time approximation algorithms for the Ising model* (SIAM J. Comput. 22:1087–1116, 1993)
First FPRAS for the partition function of the *ferromagnetic* Ising model with arbitrary external field via subgraph-world MCMC. Holds for all temperatures in the ferromagnetic regime on any graph; combined with Sly–Sun this completes the dichotomy: ferromagnetic = always tractable, antiferromagnetic = tractable iff in tree-uniqueness region.

### **[T1]** Weitz — *Counting independent sets up to the tree threshold* (STOC 2006)
Deterministic polynomial-time approximation (FPTAS) for the hardcore partition function on bounded-degree graphs in the tree-uniqueness region via the self-avoiding-walk tree reduction. Sharp threshold: $\lambda < \lambda_c(\mathbb{T}_\Delta)$ is in P, $\lambda > \lambda_c(\mathbb{T}_\Delta)$ is NP-hard (Sly); the SAW-tree reduction is the standard trick for converting tree spatial mixing into algorithms on general graphs.

### **[T1]** Cuff et al. — *Glauber dynamics for the mean-field Potts model* (J. Stat. Phys. 149:432–477, 2012) [arXiv:1204.4503]
Sharp three-regime mixing analysis of mean-field Potts Glauber: rapid above $\beta_s$, $\Theta(n^{4/3})$ at $\beta_s$, torpid below $\beta_s$ with $t_{\mathrm{mix}} = e^{\Theta(n)}$. The first-order transition for $q \ge 3$ creates exponentially slow metastability that cluster algorithms (Galanis–Štefankovič–Vigoda) and tempering (Bhatnagar–Randall) also cannot resolve in general.

## 6. Modern probabilistic methods — spectral / entropic independence and stochastic localization

### Anari, Liu, Oveis Gharan, Vinzant — *Log-concave polynomials I-IV* (FOCS 2018, STOC 2019, STOC 2021)
Series proving log-concavity of basis-generating polynomials of matroids and using it to design FPRAS for counting bases of matroids, resolving the Mihail–Vazirani conjecture. Reframes high-dimensional expansion of simplicial complexes as a polynomial inequality, opening the door to all subsequent spectral-independence results.

### **[T1]** Anari, Liu, Oveis Gharan — *Spectral independence in high-dimensional expanders and applications to the hardcore model* (FOCS 2020) [arXiv:2001.00303]
Introduces spectral independence: a uniform spectral bound on a pairwise-influence matrix, which implies optimal $O(n \log n)$ Glauber mixing on bounded-degree graphs. Proves rapid mixing of hardcore Glauber up to the tree-uniqueness threshold $\lambda_c(\mathbb{T}_\Delta)$, matching Weitz / Sly and beating the previous Dobrushin / log-Sobolev bounds.

### **[T1]** Chen, Liu, Vigoda — *Optimal mixing of Glauber dynamics: entropy factorization via high-dimensional expansion* (FOCS 2021) [arXiv:2011.02075]
Promotes spectral independence to a modified-log-Sobolev-inequality bound, giving $O(n \log n)$ Glauber mixing whenever spectral independence holds with bounded degree. Optimal up to constants; covers Ising and hardcore in the entire tree-uniqueness region with no degree-dependence beyond that imposed by Sly's hardness.

### Chen, Liu, Vigoda — *Rapid mixing of Glauber dynamics up to uniqueness via contraction* (SIAM J. Comput. 51:1593–1610, 2022) [arXiv:2004.09083]
Proves $O(n \log n)$ Glauber mixing for antiferromagnetic 2-spin systems up to (and exactly matching) the tree-uniqueness threshold using a contraction-based proof of spectral independence. Closes the gap between Weitz's deterministic FPTAS and Glauber dynamics; the threshold is sharp by Sly–Sun hardness.

### Anari, Jain, Koehler, Pham, Vuong — *Entropic independence I, II* (STOC 2022 / FOCS 2022) [arXiv:2106.04105, arXiv:2111.03247]
Defines entropic independence as the entropic analogue of spectral independence and proves modified log-Sobolev inequalities for fractionally log-concave distributions and high-temperature Ising. Gives nearly-linear-time samplers for hardcore and Ising whose constants depend only on a relative gap to the tree-uniqueness threshold, removing the maximum-degree dependence of Chen–Liu–Vigoda.

### Eldan — *Stochastic localization and the KLS conjecture* (Geom. Funct. Anal. 23:532–569, 2013)
Introduces stochastic localization, which decomposes a high-dimensional log-concave measure into a martingale of progressively localized measures. Foundational machinery later used by Eldan–Koehler–Zeitouni, Anari–Jain–etc., and Chen–Eldan to prove mixing bounds that escape the worst-case spectral-gap regime.

### **[T1]** Chen, Eldan — *Localization schemes: a framework for proving mixing bounds for Markov chains* (FOCS 2022) [arXiv:2203.04163]
Unifies spectral / entropic independence (discrete) and stochastic localization (continuous) under one framework of "localization schemes" attached to martingales of probability measures. Recovers $O(n \log n)$ mixing for hardcore in tree uniqueness, optimal mixing for Ising at any external field with $\|J\|_{\mathrm{op}} < 1$, and KL-divergence decay for log-concave sampling.

### Vempala, Wibisono — *Rapid convergence of the unadjusted Langevin algorithm: isoperimetry suffices* (NeurIPS 2019) [arXiv:1903.08568]
Proves $O(d/\varepsilon)$-step KL convergence of the unadjusted Langevin algorithm assuming only a Poincaré or log-Sobolev inequality on the target. Continuous-state analogue of the discrete spectral-gap framework; covers any non-log-concave continuous Gibbs measure with a functional inequality.

## 7. Sign problem and quantum Monte Carlo for quantum Hamiltonians

### Loh, Gubernatis, Scalettar, White, Sugar, Sugiyama — *Sign problem in the numerical simulation of many-electron systems* (PRB 41:9301–9307, 1990)
Original quantitative statement of the fermion sign problem in determinantal QMC: average sign $\langle s \rangle \sim e^{-\beta N}$ for generic $H$. Universal in $\beta$ and system size $N$ for non-stoquastic fermionic problems away from half filling on bipartite lattices; this is the regime where almost every quantum-vs-classical question lives.

### **[T1]** Troyer, Wiese — *Computational complexity and fundamental limitations to fermionic quantum Monte Carlo simulations* (PRL 94:170201, 2005) [cond-mat/0408370]
Proves that a generic polynomial-time solution to the sign problem would put NP in BPP, ruling out a universal sign-problem cure unless the polynomial hierarchy collapses. Hardness is for the *worst case*; specific models can still be sign-problem-free for structural reasons (stoquasticity, bipartite half filling, etc.).

### **[T1]** Bravyi, DiVincenzo, Oliveira, Terhal — *The complexity of stoquastic local Hamiltonian problems* (Quantum Inf. Comput. 8:361–385, 2008) [quant-ph/0606140]
Defines stoquastic Hamiltonians (real, non-positive off-diagonal entries in the computational basis) and proves their local-Hamiltonian problem is in AM and MA-hard. Stoquastic-LH is strictly easier than the generic local-Hamiltonian problem (which is QMA-complete); these are exactly the Hamiltonians without a sign problem in the standard basis.

### **[T1]** Bravyi, Terhal — *Complexity of stoquastic frustration-free Hamiltonians* (SIAM J. Comput. 39:1462–1485, 2010) [arXiv:0806.1746]
Proves that adiabatic evolution of frustration-free stoquastic Hamiltonians with inverse-polynomial spectral gap can be classically simulated in randomized polynomial time. Frustration-free + stoquastic + gap = classically tractable; this is the canonical "easy quantum case" benchmark for any quantum-Gibbs-sampling claim.

### Bravyi — *Monte Carlo simulation of stoquastic Hamiltonians* (Quantum Inf. Comput. 15:1122–1140, 2015) [arXiv:1402.2295]
Classical algorithm for simulating ground-state expectation values of frustration-free stoquastic Hamiltonians via path-integral MCMC with a provable bound on autocorrelation time. Polynomial in $1/\Delta$ (gap) and system size; demonstrates that the entire frustration-free stoquastic class is *provably* polynomially simulable classically.

### **[T1]** Bravyi, Gosset — *Polynomial-time classical simulation of quantum ferromagnets* (PRL 119:100503, 2017) [arXiv:1612.05602]
Polynomial-time classical algorithm for the partition function of stoquastic transverse-field ferromagnetic XY Hamiltonians via a Suzuki–Trotter mapping to weighted perfect matchings. Specific to ferromagnetic XY / Heisenberg on bipartite graphs; combines the FPRAS of Jerrum–Sinclair–Vigoda for permanents with quantum-to-classical reduction.

### **[T1]** Hangleiter, Roth, Nagaj, Eisert — *Easing the Monte Carlo sign problem* (Sci. Adv. 6:eabb8341, 2020) [arXiv:1906.02309]
Defines a quantitative measure of non-stoquasticity (average sign or related operator-norm distance to stoquastic) and shows how to reduce it via efficiently computable local basis changes. Optimizing over local bases is generically NP-complete (reduction to MaxCut), so curing the sign problem by basis change is provably hard in the worst case but heuristically useful.

### **[T1]** Klassen, Marvian, Piddock, Ioannou, Hen, Terhal — *Hardness and ease of curing the sign problem for two-local qubit Hamiltonians* (SIAM J. Comput. 49:1332–1362, 2020) [arXiv:1906.08800]
Proves NP-hardness of deciding whether a two-local Hamiltonian with one-local terms can be made stoquastic by single-qubit unitaries, and gives a polynomial-time algorithm when no one-local terms are present. Sharp dichotomy: in the "no one-local terms" subclass (e.g. pure XYZ Heisenberg) sign-curing is tractable, otherwise it is NP-hard.

### Marvian, Lidar, Hen — *On the computational complexity of curing non-stoquastic Hamiltonians* (Nat. Commun. 10:1571, 2019) [arXiv:1802.03408]
Strengthens the basis-change hardness picture: shows that even partial sign-easing by extended Clifford-like transformations is intractable for generic Hamiltonians. Confirms that the boundary stoquastic-vs-non-stoquastic is computationally robust — the sign problem is not generically curable by any polynomial-time preprocessing.

### Smith, Hayata — *Easing the sign problem in lattice field theory via complexified path-integral contour deformations* (Phys. Rev. D 95:094501, 2017)
Reformulates lattice field theory on Lefschetz thimbles / complexified contours to suppress the sign problem in finite-density lattice QCD and Hubbard-type models. Reduces the sign problem in specific regimes but does not eliminate it; representative of the contour-deformation literature where worst-case Troyer–Wiese hardness is sidestepped on physically-motivated subclasses.

## 8. Path-integral / world-line / determinantal QMC successes

### Sandvik, Kurkijärvi — *Quantum Monte Carlo simulation method for spin systems* (PRB 43:5950–5961, 1991)
Stochastic series expansion (SSE) representation: high-temperature expansion of $e^{-\beta H}$ truncated at random order, sampled by MCMC. Foundational for sign-problem-free spin systems on bipartite lattices; the SSE basis underlies almost every quantitative numerical result on the antiferromagnetic Heisenberg model.

### **[T1]** Syljuåsen, Sandvik — *Quantum Monte Carlo with directed loops* (PRE 66:046701, 2002) [cond-mat/0202316]
Directed-loop update for SSE and world-line QMC that satisfies detailed balance with fewer rejected moves and full coverage in magnetic field. Standard high-performance QMC method for sign-problem-free quantum spin systems; eliminates ergodicity issues that plagued the earlier loop algorithm and now ships in ALPS / quantum-Monte-Carlo libraries.

### **[T1]** Prokof'ev, Svistunov, Tupitsyn — *"Worm" algorithm in quantum Monte Carlo simulations* (Phys. Lett. A 238:253–257, 1998)
Worm update on the extended configuration space of broken world-lines, allowing direct access to off-diagonal Green's functions, winding numbers, and grand-canonical updates. Best-in-class for bosonic / sign-problem-free models; failure mode is the same as all QMC — fermions, frustration, or strong disorder cause sign or local-trapping issues.

### Boninsegni, Prokof'ev, Svistunov — *Worm algorithm and diagrammatic Monte Carlo for continuous-space path-integral simulations* (PRE 74:036701, 2006)
Continuous-space worm algorithm for bosonic atoms (helium-4, ultracold gases) without lattice discretization. Polynomially scaling for sign-free bosonic simulations; underlies most reliable quantitative results on superfluid-insulator transitions in continuous space.

### **[T1]** Gull, Millis, Lichtenstein, Rubtsov, Troyer, Werner — *Continuous-time Monte Carlo methods for quantum impurity models* (Rev. Mod. Phys. 83:349, 2011) [arXiv:1012.4474]
Comprehensive review of CT-INT, CT-AUX, CT-HYB algorithms for sampling quantum impurity models in continuous imaginary time without Trotter error. The dominant solvers for DMFT and quantum-cluster methods at finite temperature; sign problem still bites away from particle–hole symmetry and at low $T$.

### Prokof'ev, Svistunov — *Bold diagrammatic Monte Carlo* (PRL 99:250201, 2007)
Diagrammatic MCMC sampling of the bold (skeleton) expansion of the self-energy / vertex, with self-consistent renormalization. Important for fermionic systems where standard QMC has a hard sign problem; convergence is empirical (perturbation series may diverge), so guarantees are weaker than for SSE / worm.

### Pollet — *Recent developments in quantum Monte Carlo simulations with applications for cold gases* (Rep. Prog. Phys. 75:094501, 2012) [arXiv:1206.0781]
Modern review of path-integral / worm / determinantal QMC for cold-atom Hubbard and Bose–Hubbard models. Identifies sign-problem-free corners (bosons, half-filled bipartite fermions) where QMC is the gold standard and the sign-problem-bound regimes (frustration, finite-density fermions) where it is unusable.

### Pan, Meng — *The sign problem in quantum Monte Carlo simulations* (Encyclopedia of Condensed Matter Physics, 2024) [arXiv:2204.08777]
Survey of sign-problem manifestations and modern partial cures (basis transformations, reweighting, contour deformation, Majorana / time-reversal symmetry tricks) in fermionic QMC. Up-to-date snapshot of the practical sign-problem landscape, useful for assessing which quantum-vs-classical comparisons are even meaningful.

## 9. Tensor-network finite-temperature methods

### Verstraete, García-Ripoll, Cirac — *Matrix product density operators: simulation of finite-temperature and dissipative systems* (PRL 93:207204, 2004) [cond-mat/0406426]
Introduces matrix product density operators (MPDOs / purification ancilla method) and a variational imaginary-time TDVP algorithm to evolve them. Polynomial in bond dimension for 1D systems satisfying a thermal area law; the foundational tensor-network technique for finite-$T$ in 1D.

### Feiguin, White — *Finite-temperature density matrix renormalization using an enlarged Hilbert space* (PRB 72:220401, 2005) [cond-mat/0510124]
Implements purification / ancilla-DMRG to compute thermal expectation values via $\rho \propto \text{tr}_{\text{anc}}|\beta/2\rangle\langle\beta/2|$. Workhorse for finite-$T$ 1D quantum spin chains, especially when combined with matrix-product compression at each Trotter step.

### **[T1]** White — *Minimally entangled typical thermal states* (PRL 102:190601, 2009) [arXiv:0902.4475]
Algorithm that samples MPS-friendly typical thermal states by alternating imaginary-time evolution and projective measurement in a classical basis. $10^3$–$10^{10}$ speedup over ancilla DMRG in benchmarks; the entanglement of each METTS is much smaller than of the purification, so accuracy is reached at smaller bond dimension.

### Stoudenmire, White — *Minimally entangled typical thermal state algorithms* (New J. Phys. 12:055026, 2010) [arXiv:1002.1305]
Detailed implementation of METTS, including variance estimators, sampling strategies, and benchmark accuracy on quantum spin chains. Establishes METTS as a competitive 1D finite-$T$ method; the regime of utility is the same as MPS more generally — gapped or slightly-correlated thermal states obeying the thermal area law.

### **[T1]** Czarnik, Dziarmaga — *Projected entangled-pair states at finite temperature: imaginary-time evolution with ancillas* (PRB 92:035120, 2015) [arXiv:1209.0454]
2D PEPS with ancillas evolved in imaginary time to represent thermal states of 2D quantum lattices, benchmarked on the 2D quantum Ising model. Polynomially scaling in PEPS bond dimension; cost grows rapidly with $\beta$ but stays tractable above the thermal correlation length scale.

### Czarnik, Dziarmaga — *Variational tensor network renormalization in imaginary time* (PRB 92:035152, 2015) [arXiv:1503.01077]
Variational PEPO + tree-tensor-network ansatz for thermal states of 2D models, optimized to maximize accuracy at fixed bond dimension. Useful for moderately correlated 2D quantum thermal states; no sign problem since this is purely classical tensor-network optimization.

### Hastings — *Solving gapped Hamiltonians locally* (PRB 73:085115, 2006) [cond-mat/0508554]
Cluster-expansion proof that thermal states of gapped local Hamiltonians admit MPO approximations with bond dimension polynomial in $1/\varepsilon$ in 1D. Combined with the area law for 1D ground states, provides the rigorous justification for MPDO / METTS methods on gapped 1D systems.

### **[T1]** Molnar, Schuch, Verstraete, Cirac — *Approximating Gibbs states of local Hamiltonians efficiently with PEPO* (PRB 91:045138, 2015) [arXiv:1406.2973]
Constructs explicit PEPO approximations to thermal states of any local Hamiltonian with bond dimension $\exp(O(\beta))$. Polynomial in $1/\varepsilon$ at fixed $\beta$ and any spatial dimension; the bond-dimension growth in $\beta$ is essentially tight without further (low-correlation) assumptions.

### **[T1]** Kuwahara, Alhambra, Anshu — *Improved thermal area law and quasilinear time algorithm for quantum Gibbs states* (PRX 11:011047, 2021) [arXiv:2007.11174]
Improves the thermal-area-law temperature dependence from $O(\beta)$ to $\tilde{O}(\beta^{2/3})$ and gives a quasi-linear-time classical algorithm for constructing 1D MPS Gibbs states at $\beta = o(\log n)$. Quasi-linear in $n$ for 1D at any temperature below logarithmic-in-$n$; the strongest rigorous classical guarantee for 1D quantum Gibbs sampling.

### **[T1]** Kuwahara, Kato, Brandão — *Clustering of conditional mutual information for quantum Gibbs states above a threshold temperature* (PRL 124:220601, 2020) [arXiv:1910.09425]
Proves exponential decay of conditional mutual information for thermal states of local Hamiltonians above a threshold temperature, hence Markov-network approximation by local quantum recovery maps. Above-threshold regime only; together with Brandão–Kastoryano this provides the structural foundation for efficient classical and quantum Gibbs sampling at high $T$.

### Sugiura, Shimizu — *Canonical thermal pure quantum state* (PRL 111:010401, 2013) [arXiv:1302.3138]
Shows that a single Haar-random "thermal pure quantum" (TPQ) state computes any thermodynamic expectation up to exponentially small error in system size. Provides a classical sampler that beats density-matrix evolution by a Hilbert-space dimension factor; bottleneck is still exponential storage of the quantum state.

## 10. Lattice gauge theory, continuous spins, $O(n)$ models

### Duane, Kennedy, Pendleton, Roweth — *Hybrid Monte Carlo* (Phys. Lett. B 195:216–222, 1987)
Combines molecular dynamics and Metropolis to sample continuous-field theories with non-local proposals, the standard algorithm for lattice QCD. Workhorse for continuous-spin / gauge-field Gibbs sampling; mixing time governed by the integrated autocorrelation of the Hamiltonian flow, which scales polynomially in lattice volume away from criticality.

### Wolff — *Critical slowing down* (Nucl. Phys. B Proc. Suppl. 17:93–102, 1990)
Survey of dynamic critical exponents $z$ for local and cluster algorithms in $O(n)$ sigma models, lattice QCD and Ising-like systems. Establishes the empirical rule $z_{\mathrm{Metropolis}} \approx 2$, $z_{\mathrm{Wolff/SW}} \approx 0$–$1$ near criticality; sets the benchmark targets for any new algorithm.

### Beard, Wiese — *Simulations of discrete quantum systems in continuous Euclidean time* (PRL 77:5130–5133, 1996)
Continuous-Euclidean-time SSE / loop-cluster algorithm for discrete quantum spin systems and lattice gauge theory. Sign-problem-free for many physically-relevant gauge theories and quantum antiferromagnets; representative of the lattice-QCD-style algorithms that dominate finite-$T$ field theory.

### Adams, Chandrasekharan — *Chiral limit of strongly coupled lattice gauge theories* (Nucl. Phys. B 662:220–246, 2003)
Cluster algorithm for strongly coupled lattice gauge theory in a meron / dimer representation that solves the sign problem in this corner. Specific subclass of lattice QCD; shows that the sign problem can be circumvented by switching to dual variables in selected models.

## 11. Modern quantum-side high-temperature classical algorithms

### Harrow, Mehraban, Soleimanifar — *Classical algorithms, correlation decay, and complex zeros of partition functions of quantum many-body systems* (STOC 2020) [arXiv:1910.09071]
Quasi-polynomial-time classical algorithm for the partition function of any local quantum Hamiltonian above a model-dependent threshold temperature, via cluster expansions and zero-free regions. Holds for $\beta < \beta_*$ where $\beta_*$ is below any thermal phase transition; below $\beta_*$ the same problem is in general NP-hard, formalizing the "high-$T$ is easy classically" picture.

### **[T1]** Mann, Helmuth — *Efficient algorithms for approximating quantum partition functions* (J. Math. Phys. 62:022201, 2021) [arXiv:2004.11568]
Polynomial-time approximation algorithm for partition functions of quantum spin Hamiltonians at high temperature using the abstract polymer / cluster expansion framework of Helmuth–Perkins–Regts. Strengthens Harrow–Mehraban–Soleimanifar from quasi-polynomial to polynomial in the high-temperature regime; same regime restriction.

### **[T1]** Mann, Helmuth — *Efficient algorithms for approximating quantum partition functions at low temperature* (Quantum 7:1155, 2023) [arXiv:2201.06533]
Extends polymer-model methods to *low* temperature for stable quantum spin systems with dominant ground-state contribution. Holds in deep symmetry-broken phases with a discrete order; applies to "easy" low-$T$ corners that complement the high-$T$ regime, leaving the critical region as the genuinely hard case.

### Anshu, Arunachalam, Kuwahara, Soleimanifar — *Sample-efficient learning of quantum many-body systems* (Nat. Phys. 20:1027–1031, 2024) [arXiv:2108.04842]
Polynomial sample complexity for learning a local quantum Hamiltonian from copies of its high-temperature Gibbs state. Sample-optimal but with super-polynomial classical time; combined with subsequent work establishes that high-$T$ quantum Gibbs states reveal their generating Hamiltonian classically and efficiently.

### Bakshi, Liu, Moitra, Tang — *Learning quantum Hamiltonians at any temperature in polynomial time* (STOC 2024) [arXiv:2310.02243]
Polynomial-time algorithm to learn a local Hamiltonian at *any* temperature from polynomially many Gibbs-state copies, removing the high-temperature restriction of Anshu et al. Establishes that learning a Hamiltonian from its Gibbs state is classically tractable at all $\beta$, even though *sampling* the Gibbs state classically becomes hard at low $T$.

### **[T1]** Bakshi, Liu, Moitra, Tang — *High-temperature Gibbs states are unentangled and efficiently preparable* (Quantum 2024) [arXiv:2403.16850]
Proves that high-temperature local Hamiltonian Gibbs states are unentangled and admit efficient classical-shadow / product-state representations. Above the threshold temperature only; together with Harrow–Mehraban–Soleimanifar this is the structural reason all high-$T$ quantum Gibbs problems are classically easy.

## 12. Non-reversible / lifted Markov chain Monte Carlo

This section collects the classical (Hilbert-space-free) story of accelerating Gibbs sampling by *breaking detailed balance*. The mechanism is the closest classical analogue to what coherent quantum walks buy: a $\sqrt{\mathrm{gap}}$ mixing-time speedup obtained by enlarging the state space with an auxiliary momentum / direction variable and replacing diffusive exploration with ballistic motion. Without these references the comparison "best classical vs. quantum Gibbs sampler" overstates the quantum advantage, since the relevant classical baseline is the optimal *lift* of the natural reversible chain, not the reversible chain itself. The recent quantum-side companion papers — Apers–Sarlette (cited below at the bridge), Claudon–Piquemal–Monmarché (Nat. Commun. 2025, arXiv:2501.05868), and Li–Lu (arXiv:2505.12187) — are tracked in the quantum corpus; they are listed here only as the natural successors of the classical lift theorem.

### **[T1]** Chen, Lovász, Pak — *Lifting Markov chains to speed up mixing* (STOC 1999, pp. 275–281)
Foundational construction showing that any reversible chain admits a non-reversible lift on an enlarged state space whose mixing time is bounded above by the *square root* of the original mixing time, with a multicommodity-flow characterization of the optimum. Holds for any finite-state reversible chain; the $\sqrt{\mathrm{gap}}$ bound is tight, so the optimal lift gives at most a quadratic mixing-time speedup — exactly the speedup ceiling of the corresponding Szegedy quantum walk.

### **[T1]** Diaconis, Holmes, Neal — *Analysis of a non-reversible Markov chain sampler* (Ann. Appl. Probab. 10:726–752, 2000)
First quantitative analysis of a concrete momentum lift on the cycle / discrete interval, proving a $\Theta(N)$ vs $\Theta(N^2)$ mixing-time gain over the simple reversible random walk. The proof-of-concept that motivated the entire applied non-reversible MCMC programme (lifted Metropolis–Hastings, event-chain, PDMP samplers).

### **[T1]** Suwa, Todo — see §1 above
Geometric-allocation construction of rejection-minimizing non-reversible kernels on finite spin spaces; the canonical physics-side entry to non-reversible MCMC.

### Bernard, Krauth, Wilson — *Event-chain Monte Carlo algorithms for hard-sphere systems* (PRE 80:056704, 2009) [arXiv:0903.2954]
Original event-chain construction for hard spheres: a single update displaces an arbitrarily long chain of particles ballistically along a fixed direction until a collision triggers a swap. Roughly $100\times$ faster than single-particle Metropolis on hard disks and resolves the 2D hard-disk hexatic / liquid–solid transition; subsumed algorithmically by Michel–Kapfer–Krauth 2014 below — cited here as the historical origin and the canonical hard-sphere benchmark.

### **[T1]** Michel, Kapfer, Krauth — *Generalized event-chain Monte Carlo: rejection-free global-balance algorithms from infinitesimal steps* (J. Chem. Phys. 140:054116, 2014) [arXiv:1309.7748]
Generalizes event-chain to arbitrary continuous pairwise interactions via a factorized Metropolis filter and infinitesimal moves, including continuous classical spin systems and chiral magnets in a field. Universal over any continuous-state Hamiltonian with computable factor derivatives; the canonical rejection-free non-reversible algorithm for general particle and continuous-spin systems.

### **[T1]** Bouchard-Côté, Vollmer, Doucet — *The Bouncy Particle Sampler: a non-reversible rejection-free Markov chain Monte Carlo method* (J. Amer. Statist. Assoc. 113:855–867, 2018) [arXiv:1510.02451]
Continuous-time piecewise-deterministic Markov process (PDMP) whose state moves in straight lines with a velocity that reflects at level sets of $-\log\pi$ at state-dependent rates, sampling smooth targets on $\mathbb{R}^d$ without rejection. Polynomial scaling in dimension and rigorous geometric ergodicity for log-concave targets; the de facto general-purpose continuous-state non-reversible sampler in the Bayesian community.

### **[T1]** Bierkens, Fearnhead, Roberts — *The Zig-Zag process and super-efficient sampling for Bayesian analysis of big data* (Ann. Stat. 47:1288–1320, 2019) [arXiv:1607.03188]
PDMP sampler whose velocity components flip sign independently at state-dependent rates, with an exact data-subsampling variant whose per-step cost is $O(1)$ in dataset size. Geometric ergodicity uniform in dimension for log-concave targets and a CLT for ergodic averages; the cleanest theoretical PDMP non-reversible sampler and the natural continuous-time scaling limit of lifted Metropolis–Hastings.

### **[T1]** Andrieu, Livingstone — *Peskun–Tierney ordering for Markovian Monte Carlo: beyond the reversible scenario* (Ann. Stat. 49:1958–1981, 2021) [arXiv:1906.06197]
Extends the classical Peskun–Tierney comparison theorems for asymptotic variance from reversible chains to a broad class of non-reversible chains and PDMPs via $Q$-symmetrisation by an isometric involution. Provides the first general framework for *proving* that a non-reversible sampler beats its reversible parent in asymptotic-variance terms, closing a long-standing gap in non-reversible MCMC theory and giving the right comparison currency for the quantum-vs-classical baseline.

### **[T1]** Apers, Sarlette — *For every quantum walk there is a (classical) lifted Markov chain with faster mixing time* (Quantum Inf. Comput. 18:1109–1126, 2018) [arXiv:1712.02318]
Constructs an explicit classical lift of any reversible Markov chain that exactly reproduces (and slightly outperforms) the quantum-walk mixing distribution at every time, on an enlarged state space of size $n^2 D(G)$. Shows that the quadratic mixing-time advantage of a Szegedy-style quantum walk can be reproduced by a polynomial-cost classical lift — the conceptual hinge between the classical lift theorem and the quantum-walk speedup, and the natural setup paper for the quantum chapter.

### Vucelja — *Lifting — a nonreversible Markov chain Monte Carlo algorithm* (Am. J. Phys. 84:958–968, 2016) [arXiv:1412.8762]
Pedagogical review of lifting with worked examples on the ring, torus, complete-graph Ising, and 1D Ising. Best entry point to the discrete-lift literature for a thesis chapter; covers the mechanics without the proofs.

### Bierkens — *Non-reversible Metropolis-Hastings* (Stat. Comput. 26:1213–1228, 2016) [arXiv:1401.8087]
Modifies the Metropolis–Hastings acceptance ratio with a vorticity-matrix correction to produce non-reversible chains targeting the same stationary distribution. Provides a recipe for converting reversible Metropolis kernels to non-reversible variants; complements Suwa–Todo on the geometric-allocation side.

### Bierkens, Roberts — *A piecewise deterministic scaling limit of Lifted Metropolis–Hastings in the Curie–Weiss model* (Ann. Appl. Probab. 27:846–882, 2017) [arXiv:1509.00302]
Derives the Zig-Zag process as a scaling limit of discrete lifted Metropolis–Hastings on the mean-field Curie–Weiss Ising model. Shows how the discrete-lift framework of Chen–Lovász–Pak and the continuous-time PDMP framework of Bierkens–Fearnhead–Roberts are two sides of the same construction.

### Hwang, Hwang-Ma, Sheu — *Accelerating Gaussian diffusions* (Ann. Appl. Probab. 3:897–913, 1993)
Foundational continuous-state non-reversibility result: adding a divergence-free drift to overdamped Langevin strictly improves the spectral gap of the resulting non-reversible diffusion targeting the same Gaussian. The continuous analogue of the Chen–Lovász–Pak lift; later sharpened by Hwang–Hwang-Ma–Sheu (Ann. Appl. Probab. 2005) and by Lelièvre–Nier–Pavliotis below.

### Lelièvre, Nier, Pavliotis — *Optimal non-reversible linear drift for the convergence to equilibrium of a diffusion* (J. Stat. Phys. 152:237–274, 2013) [arXiv:1212.0876]
Sharpens Hwang's analysis with an explicit construction of the optimal divergence-free drift for Gaussian targets and proves that this drift achieves a square-root mixing-time gain. The continuous-state counterpart to the discrete $\sqrt{\mathrm{gap}}$ bound of Chen–Lovász–Pak.

### Alon, Benjamini, Lubetzky, Sodin — *Non-backtracking random walks mix faster* (Comm. Contemp. Math. 9:585–603, 2007) [arXiv:math/0610550]
Shows that non-backtracking random walks on regular Ramanujan-like expanders mix up to twice as fast as the simple reversible random walk, with the speedup ratio controlled by the spectral-radius gap. Niche in scope but illustrative of the general principle that breaking time-reversal symmetry on a structured graph can give a constant-factor mixing improvement, complementing the lift-based superpolynomial story.

## 13. Surveys, lecture notes, and textbook references

### Sinclair — *Algorithms for Random Generation and Counting: A Markov Chain Approach* (Birkhäuser, 1993)
Foundational monograph linking approximate counting, sampling, and rapidly-mixing Markov chains via the canonical-paths and conductance machinery. Standard reference for the conductance-based approach used as a baseline before spectral independence; covers the FPRAS framework that defines "efficient" classical Gibbs sampling.

### Jerrum — *Counting, Sampling and Integrating: Algorithms and Complexity* (Birkhäuser, 2003)
Pedagogical monograph on counting / sampling algorithms with a focus on Glauber dynamics, canonical paths, comparison, and FPRAS hardness reductions. Companion to Sinclair (1993); accessible introduction to the tools used by Sly, Weitz, and the spectral-independence school.

### Anari, Liu, Oveis Gharan — *Spectral independence and local-to-global techniques for optimal mixing of Markov chains* (lecture notes, 2023) [arXiv:2307.13826]
Modern lecture notes on the spectral / entropic-independence framework, with a self-contained proof that spectral independence implies optimal $O(n \log n)$ Glauber spectral gap. Best entry point to the modern spectral-independence literature; covers Ising, hardcore, matroid bases, and random clusters under one umbrella.

### Levin–Peres–Wilmer — see Section 1 above
(Listed under foundations.)

### Krauth — see Section 1 above
(Listed under foundations.)

### Pollet — *Recent developments in quantum Monte Carlo simulations* — see Section 8 above
(Listed under QMC successes.)

### Gull et al. — *Continuous-time Monte Carlo methods for quantum impurity models* — see Section 8 above
(Listed under QMC successes.)
