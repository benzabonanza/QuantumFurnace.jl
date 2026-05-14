# Literature Synthesis — All

**Scope**: all
**Regime filter**: none
**Reviews ingested**: 99 (39 classical + 60 quantum; out of 39/60 on disk)
**Generated**: 2026-05-14
**Output goal**: structured knowledge dump for the thesis Review chapter (qf-yk4) — *not* polished prose.

The chapter's central question is **"where could quantum Gibbs sampling first beat classical Gibbs sampling?"** This synthesis collates **both sides** of the comparison — the classical-side answers (what classical methods can already do, structural barriers) and the quantum-side answers (what KMS-DB samplers and their alternatives can and cannot do) — and surfaces the cross-scope contact points where the two sides meet on a common axis. The companion files `literature-review-synthesis-classical.md` (May 8) and `literature-review-synthesis-quantum.md` (May 11) are the source per-scope syntheses; this file collates them and adds the cross-scope contact-point and frontier analysis.

The user's framing for this file: **be especially thorough on three categories within the contact-point / frontier analysis** — (1) direct comparisons with a clean axis, (2) "verdict unclear" pairs where the two sides exist on the same model but with incommensurate access / precision / cost metrics, (3) open-frontier cells where neither side has a meaningful result or only one side has been attempted.

---

## Section 1 — Corpus overview

### By regime cell (temperature × commutativity)

| Cell key | Classical count | Quantum count | Cross-scope count |
|---|---|---|---|
| high-T × commuting × sparse / geometric-local × spin | 6 | 4 (Davies / MLSI program) | — |
| high-T × noncommuting × geometric-local / k-local × spin | 4 | 6 + 2 (cross-scope: Yin–Lucas, Zlokapa) | 2 |
| intermediate-T × commuting × geometric-local / mean-field × spin | 3 | 3 (BCL designer, CTH) | — |
| intermediate-T / constant-T × noncommuting × k-local × spin | — | 3 (Bergamaschi, Chen–Lamb, Chen–Markov) | 1 (Zlokapa for SYK) |
| low-T × commuting × geometric-local / mean-field × spin | 5 | 4 (toric code, code-Hamiltonian, BK19) | — |
| low-T × noncommuting × geometric-local / k-local × spin | 1 (stable perturbation) | 6 (Tong–Zhan, Smid, Ding–lowT, Gamarnik, Rakovszky, Becker) | — |
| low-T × non-commuting × fermionic | 2 (Gull, Troyer–Wiese) | 2 (Tong–Zhan, Smid) | — |
| any-T × commuting × arbitrary-graph × spin | 2 (RC, FK / Edwards–Sokal) | — | — |
| any-T × noncommuting × geometric-local / 1D / 2D × spin / boson | 7 (TN, SSE, worm, METTS, PEPO, KAA, Bravyi–Gosset) | 14 (KMS-DB core: CKG / Chen 2025 / Ding / Scandi / Chen–Rouzé / Leng / Fang / Jiang / Lin) | — |
| any-T × bosonic / infinite-dim | — | 1 (Becker) | — |
| ground-state × stoquastic / frustration-free × k-local × spin / boson | 3 (BDOT, Bravyi–Terhal, Bravyi–Gosset) | 4 (Motta, Chowdhury, Poulin–Wocjan, Chiang) | — |
| non-stoquastic / sign-problem-bearing × 2-local × spin / fermion | 2 (Hangleiter, Klassen–Marvian) | — | 2 (Klassen, Hangleiter feed into KMS arguments) |
| high-T / intermediate-T × commuting × bounded-degree × n/a (correlation-decay) | 1 (Weitz) | — | — |
| n/a × n/a (lifting / non-reversible / PDMP — generic targets) | 6 (CLP, DHN, AL, ZZ, BPS, MKK) | 7 (Apers ×2, Claudon, Li–Lu ×2, Lu, Fang–hypocoercivity) | The lifting bridge |

(Cells overlap when papers carry multi-valued temperature/commutativity tags. Counts include cross-listed papers in multiple cells.)

### By paradigm

| Paradigm | Count |
|---|---|
| `classical-MCMC` (Glauber / Metropolis on lattice spin systems) | 12 |
| `classical-MCMC + classical-lifted` / `classical-PDMP` | 6 (CLP, DHN, AL, ZZ, BPS, MKK) |
| `classical-tensor-network` (PEPS / MPS / METTS / PEPO / area law) | 5 (Czarnik–Dziarmaga, KAA, Molnar, White, plus Hastings via KAA refs) |
| `classical-cluster` (FK / Swendsen–Wang / random-cluster) | 2 (Edwards–Sokal, Guo–Jerrum) |
| `classical-deterministic` / `cluster-expansion` / `polymer-model` | 5 (Mann–Helmuth ×2, Weitz, BLMT, KKB) — plus quantum-side Yin–Lucas, Zlokapa |
| `classical-QMC` (CT-QMC / directed-loop SSE / worm / fermionic sign hardness) | 4 (Gull, Syljuåsen–Sandvik, Prokof'ev, Troyer–Wiese) |
| `classical-other` (Bravyi–Gosset, BLMT, Hangleiter) | 3 |
| `structural` / `complexity-theoretic` / `lower-bound` / `hardness` / `basis-change` | 7+ on the classical side; 5+ on quantum side |
| `parallel-tempering` | 2 (Madras–Zheng, Woodard–Schmidler–Huber) |
| `quantum-Lindblad-KMS` (the spine) | 23 |
| `quantum-Lindblad-other` (Davies, hypocoercive, MLSI side) | 12 |
| `quantum-dissipative-other` (system-bath, collision, quasiparticle) | 11 |
| `quantum-block-encoding` | 8 |
| `quantum-walk` | 5 |
| `quantum-phase-estimation-based` | 4 |
| `quantum-lifted` / `quantum-non-reversible` (cross-scope) | 4 (Apers ×2, Claudon, Li–Lu ×2, Lu) |
| `quantum-spin-glass` rigorous | 0 — `[GAP]` Aharonov–Gershoni–Klein 2024 cited but no review |

### By result-type

| Result type | Count |
|---|---|
| `mixing-time-upper` / `gap-bound` / `mlsi` / `cutoff` | 13 (classical) + 30 (quantum) |
| `mixing-time-lower` / `no-go` / `hardness` / `complexity-lower` / `NP-hard` | 8 (classical) + 6 (quantum) |
| `FPTAS` / `FPRAS` / `polynomial-time-deterministic` / `complexity-upper` | 7 (classical) + many quantum upper bounds |
| `structural` / `decay-of-correlations` / `representation-theoretic` / `bond-dimension-bound` / `area-law-bound` | 8 (classical) + 13 (quantum) |
| `runtime-comparison` / `algorithm-heuristic` / `benchmark` / `algorithm` | 6 (classical) + 9 (quantum) |
| `critical-scaling` / `threshold-temperature` | 2 |

### Borderline inclusions

Listed once per cross-scope routing decision.

- **`dervovic-2018-quantum-walk-classical-lift`** (in `quantum-papers-reviewed/`) — frontmatter `paradigm: [classical-lifted, quantum-other, comparison]`, `quantum-or-classical: comparison`. The paper is the structural bridge between classical lifting (CLP) and quantum walks; routed to quantum scope by the rule "paradigm contains quantum-other".
- **`apers-2017-lifting-mix-faster`** (in `quantum-papers-reviewed/`) — frontmatter `paradigm: [classical-MCMC, classical-lifted, comparison]`, `quantum-or-classical: classical`. Per the agent rule, this is *classical* on its tags, but lives in the quantum directory because it is the structural ceiling paper that classical and quantum scopes both need. Listed as borderline in *both* per-scope syntheses; appears in this `all`-scope synthesis in Section 4 (lifting thread).
- **`apers-2018-quantum-walks-classical-lift`** — `quantum-or-classical: comparison`, paradigm contains `quantum-walk`. Quantum-side body content (Szegedy walk simulation). Routed to quantum.
- **`claudon-2025-nonreversible-quantum`** — `comparison`, paradigm contains `quantum-walk, quantum-block-encoding`. Body on the quantum side (GQSVT on classical chain). Routed to quantum.
- **`li-lu-2025-quantum-lifting`**, **`li-lu-2025-spacetime-poincare`**, **`lu-2025-hypocoercivity-lifting`** — all `comparison` paradigms but the quantum lift / KMS-DB hypocoercivity content is quantum-side. Routed to quantum.
- **`bakshi-liu-moitra-tang-2024-high-T-unentangled`** — `[classical, comparison]`, paradigm `[classical-other, structural]`. The algorithmic content (depth-1 product-state preparation) is purely classical — only the *target* is quantum. Routed to classical.
- **`bravyi-divincenzo-oliveira-terhal-2008-stoquastic-lh`** — `[comparison, classical-side-baseline]`, paradigm `[structural, complexity-theoretic]`. Pure complexity-theoretic structural result defining the classical-side baseline (the stoquastic class). Routed to classical.
- **`kuwahara-kato-brandao-2020-cmi-clustering`** — `[comparison, structural]`, paradigm `[structural, cluster-expansion]`. Structural decay-of-correlations + classical poly-time FPATS for thermodynamic quantities. Routed to classical.
- **`hangleiter-roth-nagaj-eisert-2020-easing-sign-problem`**, **`klassen-marvian-2020-curing-sign-problem`** — both `[classical, comparison]`; classical sign-problem hardness results. Routed to classical (both feed quantum-side punchline contact points).
- **`troyer-wiese-2005-sign-problem-hardness`** — `[classical, lower-bound]`. Pure classical-QMC sign-problem hardness. Routed to classical.
- **`sly-2010-computational-transition-uniqueness`** — `[classical, lower-bound]`. Classical hardness result. Routed to classical.
- **`yin-lucas-2023-classical-high-T`** (in `quantum-papers-reviewed/`) — `[classical]`, paradigm `[classical-other, cluster-expansion]`. A *classical* upper bound *for a quantum target*. Strictly classical by the inclusion rule but lives in quantum directory as the canonical "high-T quantum target classical algorithm" entry; we include it as classical here.
- **`zlokapa-2026-syk-classical`** (in `quantum-papers-reviewed/`) — `[classical, comparison]`, paradigm `[classical-other, cluster-expansion, polymer-model, Barvinok-interpolation]`. Same treatment as Yin–Lucas — classical algorithm for a quantum target. Included as classical for this synthesis.
- **`gamarnik-2024-slow-mixing`** — `[lower-bound]`, paradigm contains `lower-bound`. Body bounds *quantum* Gibbs sampler mixing times. Routed to quantum (quantum lower bound).
- **`rakovszky-2024-bottlenecks`** — `[quantum, lower-bound]`, paradigm `[lower-bound, structural, quantum-dissipative-other]`. Routed to quantum.
- **`ilin-lychkovskiy-2021-quantum-speed-limit`** — `[lower-bound]`. Closed-system unitary QSL. Routed to quantum.

### Skipped reviews

None — all 99 review files parsed and ingested.

---

## Section 2 — Regime cells

Cells ordered from most-populated to least-populated. Cell key: `<temperature> × <commutativity> × <locality> × <particle-statistics>`. `n/a` preserved literally. For each cell, papers from both scopes are listed together so the reader can see which side has results and which doesn't.

For brevity, the per-paper entries summarise the per-paper review one-sentence takeaway, key scaling, and the load-bearing "quantum vs classical" line. For full reviews see `drafts/literature-review/classical-papers-reviewed/` and `drafts/literature-review/quantum-papers-reviewed/`.

### all-T × noncommuting × local/k-local/geometric-local × spin

**The KMS-DB Lindbladian core — most populated cell (14 quantum-side reviews + several classical cross-cuts).**

The quantum-side spine: CKG (Chen 2023), Chen 2025 (exact DB), DLL 2024, Ding 2024 KMS, Gilyén 2024 Glauber-Metropolis, Scandi–Alhambra 2025 (physical derivation), Chen 2026 strong-Markov, Chen–Rouzé 2025 (locally Markovian), Fang 2024 (hypocoercivity), Li–Lu 2025 spacetime Poincaré, Fang 2026 detectability, Leng 2026 without-walks, Jiang 2026 single-trajectory, Lin 2025 (perspective).

Classical-side counterparts in this cell: the tensor-network / SSE / worm / METTS / PEPO / KAA family that lives at any-$T$ noncommuting (1D and 2D quantum spin / boson). See "any-T × noncommuting × geometric-local" below.

**Headline quantum scalings** (verbatim from per-paper reviews):
- **chen-2023-thermal-prep**: `Lindbladian sim: \widetilde{O}(beta * t_mix^2 / eps); coherent sampler: \widetilde{O}(beta^2 ||H|| / lambda_gap^{3/2} + beta / (eps * lambda_gap^{3/2}))`. *Promising or not*: extremely promising for **non-stoquastic × low-T** corner. *No $t_{\mathrm{mix}}$ for specific $H$*.
- **chen-2025-exact-noncommutative**: `O~(t_mix * beta)` per Gibbs sample; on $D$-dim lattice gate cost per unit time $\sim \beta(v_{LR}\beta)^D$ system-size-independent. *Promising or not*: very. Necessary precondition for any future MCMC-style quantum advantage — but $t_{\mathrm{mix}}$ not bounded.
- **ding-2024-kms-samplers**: `O~(C_q t_mix beta^2 S |A|^2 log^{1+s}(1/eps))`. Recovers CKG as Gaussian special case.
- **gilyen-2024-glauber-metropolis**: discrete-time channel via QSVT with $\widetilde{O}(\log(1/\varepsilon))$ Ham-sim time when $\|\mathcal{T}'^\dagger[\mathbf{I}]\|\le 1/2$. Closes "exact quantum Metropolis" gap.
- **scandi-alhambra-2025-kms-thermalization**: approximation error $O(\alpha (\sqrt{\Gamma_0\tau}\Gamma\tau + (\alpha\Gamma)^2\beta t)e^{(\alpha\Gamma\beta)^2})$, linear-in-$t$ (vs prior exponential).
- **chen-rouze-2025-locally-markovian**: Gibbs prep depth $e^{O(\log^D(n/\varepsilon))}$ unconditional under uniform clustering.
- **fang-2026-detectability**: stationary-state prep $\widetilde{O}(Mg^2/\mathrm{gap}(\mathcal{L}) \log(1/(\sigma_{\min}\varepsilon)))$; **commuting parent** gives $1/\sqrt{\gamma}$.
- **leng-2026-without-walks**: $\mathcal{O}(\sqrt{J/\Delta}\log(1/\varepsilon))$ queries for any KMS-DB Lindbladian — closes Davies-only restriction.
- **jiang-2026-single-trajectory**: total evolution time $t_{\mathrm{mix}} + N \cdot t_{\mathrm{aut}}$, $N = O(\mathrm{var}_H/\varepsilon^2)$.
- **lin-2025-dissipative-preparation**: framing reference; Lindblad-sim query cost $O(t\|\mathcal{L}\|_{be} \log(t\|\mathcal{L}\|_{be}/\varepsilon))$.

The quantum-side cell is dense; the classical-side counterpart (worm, SSE, METTS, PEPO) is non-empty and gives the bar quantum must clear in sign-free regimes. See cross-scope contact points in Section 5.

---

### any-T × noncommuting × geometric-local / 1D / 2D × spin / boson — the classical TN / SSE / worm cell

**7 classical reviews.** Models: Heisenberg, XXZ, TFIM, XY, Bose–Hubbard, 1D / 2D quantum lattice.

- **bravyi-gosset-2017-quantum-ferromagnets** — 2017 — FPRAS for $\mathcal Z(\beta,H)$ with runtime $\mathrm{poly}(n,\beta,\varepsilon^{-1})$; explicit $\widetilde O(n^{115}(1+\beta^{46})\varepsilon^{-25})$. Closes the "stoquastic-ferromagnetic-XY" cell as classically tractable at *all* finite temperatures.
- **syljuasen-sandvik-2002-directed-loops-qmc** — 2002 — Directed-loop SSE per-MC-step cost $\sim O(N\beta)$; $\tau_{\mathrm{int}}$ empirically $O(1)$ in $N\beta$ for isotropic Heisenberg AFM. The practical classical frontier for stoquastic × all-T.
- **prokofev-svistunov-tupitsyn-1998-worm-algorithm** — 1998 — Per-update cost $O(1)$; integrated autocorrelation empirically $\sim O(L^z)$ with $z\approx 0.5$–$1$ in superfluid, $z\approx 2$ near Bose–Hubbard critical point. Canonical best-in-class for sign-free bosonic / $U(1)$-symmetric.
- **white-2009-metts** — 2009 — Per-sample cost $O(N m_0^3 \beta N_\tau)$ vs ancilla DMRG $O(N m_0^6 \beta)$; speedup $10^3$–$10^{10}$.
- **czarnik-dziarmaga-2015-peps-finite-T** — 2015 — Per-time-step $O(M^3 D^4)$. Heuristic for 2D quantum Ising; classical limit ($h=0$) exact at $D=2$.
- **molnar-schuch-verstraete-cirac-2015-pepo-gibbs** — 2015 — PEPO bond dim $D = (N/\varepsilon)^{O(\beta)}$ in any spatial dimension; $D = e^{O(\log_2^2(N/\varepsilon))}$ when $\beta < O(\log_2 N)$. **Existence theorem, not algorithm** — PEPS contraction is $\#$P-hard in general. *Narrows* candidate region for structural quantum advantage to **low-T × $d\ge 2$**.
- **kuwahara-alhambra-anshu-2021-thermal-area-law** — 2021 — generic-lattice EoF $\widetilde O(\beta^{2/3})\cdot|\partial A|$; 1D MPO bond dim $\chi = \exp(\sqrt{\widetilde O(\beta\log(n/\varepsilon))})$. *Quasi-linear $\widetilde O(n)$ classical algorithm* at $\beta = o(\log n)$. *No quantum advantage* in 1D, $\beta = o(\log n)$.

---

### high-T × commuting × sparse / geometric-local × spin

**6 classical reviews + 4 quantum.**

Classical side (Glauber / log-Sobolev / spectral independence): anari-liu-oveisgharan-2020 (hardcore $\lambda<\lambda_c$, $t_{\mathrm{mix}} = O((1+\lambda)n)^{1+C(\delta)}\log(\cdot)$); chen-liu-vigoda-2021 ($O(n\log n)$ optimal in tree-uniqueness, matching Hayes–Sinclair lower bound); chen-eldan-2022 (Ising at $\|J\|_{\mathrm{op}}<1/2$, $O(n\log n)$ uniformly in field); eldan-koehler-zeitouni-2022 (SK at $\beta < 1/2$, $O(n^2/(1-\|J\|))$); bauerschmidt-dagallier-2024 (LSI from susceptibility, polynomial in $1/(\beta_c-\beta)$ for $d\ge 5$); lubetzky-sly-2013 (cutoff $d/(2\lambda_\infty)\log n$ in SSM regime).

Quantum side (Davies / MLSI on the commuting case): kastoryano-brandao-2016-commuting-case (gap $\Omega(1)$ iff strong $\mathbb L_2$-clustering, $t_{\mathrm{mix}} = O(|\Lambda|/\lambda)$); capel-2021-modified-log-sobolev (system-size-independent MLSI, $t_{\mathrm{mix}} = O(\log(|\Lambda|/\varepsilon))$); kochanowski-2024-rapid-thermalization (1D all-T: $O(\log|\Lambda|\log(1/\varepsilon))$; 2D gap+small $\xi$: $O(\sqrt{|\Lambda|}\log|\Lambda|)$; high-T $\mathbb Z^D$: $O(\log|\Lambda|)$); capel-2025-noncommutative-OT ($o(N^{1+\delta})$ circuit cost under matrix-CMI decay).

*Cell verdict*: classical commuting Gibbs sampling is essentially closed by spectral-independence / localization-schemes machinery (Hayes–Sinclair $\Omega(n\log n)$ matched); quantum-side commuting Davies matches asymptotically. **No quantum advantage** in this cell.

---

### high-T × noncommuting × geometric-local / k-local × spin

**4 classical reviews + 6 quantum + 2 cross-scope.**

Classical side closes the cell:
- **bakshi-liu-moitra-tang-2024-high-T-unentangled** — separable Gibbs for $\beta < 1/(c\mathfrak{d}\mathfrak{K})$; classical preparation in $\widetilde O(n^{6+\log\mathfrak d/\log(\beta_c/\beta)})$ via depth-1 quantum circuit.
- **kuwahara-kato-brandao-2020-cmi-clustering** — exponential CMI decay above $\beta_c = 1/(8e^3 k)$; classical FPATS for thermodynamic quantities.
- **mann-helmuth-2021-high-T-partition-functions** — FPTAS in $\mathrm{poly}(n,1/\varepsilon)$ at $|\beta|\le 1/(e^4\Delta)$; sharpens Harrow–Mehraban–Soleimanifar.
- **kuwahara-alhambra-anshu-2021-thermal-area-law** — 1D, $\widetilde O(n)$ classical algorithm at $\beta = o(\log n)$.

Cross-scope classical-for-quantum-target:
- **yin-lucas-2023-classical-high-T** — `Tr(O rho_beta)` and computational-basis sampling in $O(N^{1+c\delta})$ at $\beta < \beta_\star = [2e^2\mathfrak{d}(\mathfrak{d}+1)]^{-1}$. Closes high-T to MBQC exponential speedups.

Quantum side opens with polynomial bounds; closed asymptotically by the classical side:
- **rouze-2024-thermal-universal** — high-T spectral gap $\ge 1/(2\sqrt{2}e^{1/4})$ constant for $\beta < \beta^*$; $t_{\mathrm{mix}} = O(\log(1/\varepsilon) + |\Lambda|)$. *First unconditional efficient quantum Gibbs sampler for noncommuting*.
- **rouze-2024-optimal-gibbs** — true rapid mixing $O(\log(n/\varepsilon))$ for $\beta < 1/(615^D J)$; end-to-end $\widetilde{O}(n)$.
- **ramkumar-soleimanifar-2024-sparse-mixing** — Random sparse $H$ at $\beta\|H\| = O(1)$, $\lambda_{\mathrm{gap}} = \Omega(1)$, $t_{\mathrm{mix}} = \mathrm{polylog}(n)$.
- **bardet-2023-rapid-thermalization-spin-chains** — 1D translation-invariant commuting Hamiltonian at any $\beta$: Davies MLSI $\alpha_\Lambda = \Omega(1/\log|\Lambda|)$.

*Cell verdict*: **negative for asymptotic quantum advantage**. Classical $\widetilde O(n^7)$ (BLMT) vs quantum $\widetilde O(n)$ (Rouzé 2024) — polynomial parity. The structural row is closed.

---

### intermediate-T × commuting × geometric-local / mean-field × spin

**3 classical reviews + 3 quantum (one BCL designer, two Chen-Markov-style structural).**

Classical: lubetzky-sly-2012 (2D Ising at $\beta_c$, $t_{\mathrm{mix}} = O(n^c)$); cuff-et-al-2012 (mean-field Potts $q\ge 3$, sharp three-regime: $n\log n$ subcritical, $\Theta(n^{4/3})$ at spinodal, $e^{\Theta(n)}$ supercritical); madras-zheng-2003 (PT on Curie–Weiss Ising, $\mathrm{poly}(n)$).

Quantum (designer / structural):
- **bergamaschi-2024-constant-T-advantage** — Designer $O(\log\log n)$-local parent Hamiltonian of shallow IQP circuits: $t_{\mathrm{mix}}(\mathcal{L}_{\mathrm{Davies}}) = O(4^\ell \log n)$; classical sampling intractable under IQP-hardness. *First end-to-end constant-T quantum advantage proof* — but designer Hamiltonian, not natural.

*Cell verdict*: classical first-order-transition torpidity (Cuff et al.) is a candidate for quantum speedup; no quantum analysis exists on the *natural* mean-field model. Bergamaschi opens the constant-T cell for designer Hamiltonians.

---

### low-T × commuting × geometric-local / mean-field × spin

**5 classical reviews + 4 quantum (toric code / stabilizer / BK19).**

Classical (the no-go anchors):
- **bovier-den-hollander-2015-metastability** — $\mathbb E_{\mathbf m}[\tau_{\mathbf s}] = K e^{\beta\Gamma^*}(1+o(1))$; sharp Eyring–Kramers for reversible Markov dynamics on low-T Ising-like.
- **sly-2010-computational-transition-uniqueness** — hardcore past $\lambda_c(\Delta)$, no FPRAS unless NP=RP.
- **woodard-schmidler-huber-2009-torpid-tempering** — PT/ST on persistent-mode targets: $\mathrm{Spec} \le c_1 e^{-c_2 N}$. Mean-field Potts $q\ge 3$, mixtures of Gaussians, unconditional torpidity.
- **cuff-et-al-2012** (cross-listed, $e^{\Theta(n)}$ supercritical mean-field Potts).
- **madras-zheng-2003** (cross-listed, positive PT result on *symmetric* bimodal mean-field Ising).

Quantum: brandao-kastoryano-2019-finite-correlation ($\log(|\Lambda|/\varepsilon)$-depth circuit under uniform clustering + uniform Markov); lucia-2015-rapid-mixing-stability (stability of rapid mixing under local perturbations); kastoryano-brandao-2016 (cross-listed); ding-2024-low-temp (2D toric code, $\mathrm{poly}(N,\beta)$ via four global logical jumps).

*Cell verdict*: the classical no-go anchors (Bovier–den Hollander Eyring–Kramers, WSH persistence) define the bar. Quantum sampling of natural commuting low-T Hamiltonians (e.g. lattice Ising below $T_c$) has *no positive result on natural Hamiltonians* and the Gamarnik / Rakovszky lower bounds say none can exist locally.

---

### low-T × noncommuting × geometric-local / k-local × spin / fermionic

**1 classical review (stable perturbation only) + 6 quantum — the chapter punchline cell.**

Classical side:
- **mann-helmuth-2023-low-T-partition-functions** — Stable quantum perturbations of classical lattice spin systems with finitely many discrete ground states, deterministic FPTAS in $|V(G)|^{O(1)}(1/\varepsilon)^{O(1)}$. *Gap closed by the classical side* in stability region. Leaves criticality, frustration, AFM without symmetry breaking, *non-stoquastic* untouched.

Quantum side — the punchline:
- **tong-zhan-2024-fermionic-mixing** — $D$-dim fermionic $H_0 + V$ with weak $|U| < U_\beta = O(1)$: KMS-DB gap $g = \Omega(1)$ on even-parity sector at *any* $\beta$. $t_{\mathrm{mix}} = O(n + \log(1/\varepsilon))$, gate complexity $O(n^2 \mathrm{polylog})$.
- **smid-2025-fermi-hubbard** — Free-fermion gap $\Delta_0 = 2e^{-4\beta^2\|h\|^2}\cosh(2\beta\|h\|)$; end-to-end $\widetilde{O}(n^3 \mathrm{polylog})$ on Fermi-Hubbard *any $D$*, any constant $\beta$, $|U|/t \le U_{\max}(\beta)$. **First concrete fermionic instance with poly *quantum* + no rigorous polynomial classical**.
- **ding-2024-low-temp** — 2D toric code: $\mathrm{Gap}(-\mathcal{L}_\beta) \ge \max\{e^{-O(\beta)}, \Omega(N^{-3})\}$ via four *global* logical jumps. $t_{\mathrm{mix}} = O(\min\{\beta\mathrm{poly}(N), e^{c\beta}(N+\beta)\})$. Stabilizer-only; "Generalised energy-barrier obstruction" sidestepped, not refuted.
- **gamarnik-2024-slow-mixing** — *Lower bound*: $T_{\mathrm{mix}} = \exp[\Omega(\mathrm{poly}\,n)]$ for KMS/Lindbladian samplers on (i) classically-slow-mixing classical Hamiltonians (random $K$-SAT, $p$-spin glass, lattice Ising, Curie-Weiss), (ii) good stabilizer codes, (iii) 2D TFIM in ferromagnetic phase. *Quantum sampler at least as bad as classical Glauber on same Hamiltonian*.
- **rakovszky-2024-bottlenecks** — *Lower bound*: $t_{\mathrm{mix}} \ge \mathrm{const}\cdot\mathrm{tr}(P_A\rho)\mathrm{tr}(P_C\rho)/\|P_B\rho\|_1$; macroscopic free-energy barrier $\Omega(n)$ $\Rightarrow$ $t_{\mathrm{mix}} = \exp(\Omega(n))$. Closes broad swath of low-T regime for *any local KMS-DB Lindbladian*: SSB, q/LDPC, weakly perturbed stabilizers, classical and quantum LDPC at low $T$.
- **becker-2026-infinite-dim** — Lifts KMS-DB to *unbounded* bosonic Hamiltonians. Single-mode photon-number gap proven; Bose–Hubbard in companion. Interacting many-mode open.

*Cell verdict*: **the chapter's punchline cell**. Sign-problem-hardness barriers (Troyer–Wiese, Klassen–Marvian, Hangleiter) on the classical side make this a worst-case super-polynomial classical-hardness candidate. Quantum upper bounds exist only for (i) weakly interacting fermions (Tong–Zhan, Smid) and (ii) special stabilizer cases (toric code, Ding 2024). Quantum lower bounds (Gamarnik, Rakovszky) close most natural cases for *local* samplers. Where quantum wins on a natural physical model at genuinely low $T$ remains **open**.

---

### low-T × non-commuting × fermionic / spin (sign-problem-bearing)

**2 classical reviews (sign-problem hardness) + 2 quantum (Hubbard).**

Classical:
- **gull-et-al-2011-ct-qmc-impurity-models** — Per-sweep cost $\sim O(\langle k\rangle^3)$ with $\langle k\rangle \sim N\beta U$ (CT-INT) or $\sim N\beta$ (CT-HYB); sign-problem reweighting cost $\sim e^{2\beta\Delta F}/M$. Polynomial on sign-free corners; exponential elsewhere.
- **troyer-wiese-2005-sign-problem-hardness** — $\langle s\rangle = e^{-\beta N \Delta f}$; generic poly-time cure $\Rightarrow$ NP$\subseteq$BPP. **Removes the most obvious classical escape route for non-stoquastic Gibbs sampling**.

Quantum (Tong–Zhan, Smid 2025 in this cell as well — cross-listed).

*Cell verdict*: **highly promising** for quantum advantage. The classical sign barrier is structural; Smid's $\widetilde O(n^3)$ Hubbard upper bound in any $D$ at constant $\beta$ in the *weak-coupling* corner is the cleanest existing fermionic KMS-DB efficiency proof. Strong-coupling / Mott / doped regimes remain open on both sides.

---

### non-stoquastic / sign-problem-bearing × 2-local × spin/fermion

**2 classical reviews (sign-curing hardness).**

- **hangleiter-roth-nagaj-eisert-2020-easing-sign-problem** — SignEasing NP-complete for 2-local under on-site orthogonal Clifford / general orthogonal transformations. $\nu_1(H)$ measure efficiently computable; heuristic optimizer reduces $\nu_1$ by 2–5× on frustrated Heisenberg ladders, $\langle\mathrm{sign}\rangle^{-1}$ drops $\sim 10^4$.
- **klassen-marvian-2020-curing-sign-problem** — LocalSignCure NP-complete with one-local terms (3-SAT reduction); $O(n^3)$ arithmetic operations without one-local terms.

*Cell verdict*: identifies a structural source of quantum advantage. Quantum advantage in non-stoquastic Gibbs sampling is a *worst-case complexity separation* candidate, not a constant-factor speedup claim.

---

### intermediate-T / constant-T × noncommuting × k-local × spin

**3 quantum reviews — the constant-T advantage / hardness trio + Zlokapa SYK.**

- **bergamaschi-2024-constant-T-advantage** (cross-listed — designer commuting CTH).
- **chen-2026-lamb-shift** — $\tau_{\mathrm{mix}} = \widetilde{O}(\beta^5\|H\|^2/(\lambda_{\mathrm{gap}}^3\varepsilon))$; first $\mathcal{O}(1/\varepsilon)$ end-to-end SBI scaling.
- **chen-2026-strong-markov** — clustering + ADB on $AB$ $\Leftrightarrow$ strong local Markov (structural).

Cross-scope:
- **zlokapa-2026-syk-classical** — Classical $n^{O(\log(n/\varepsilon))}$ for SYK local thermal expectations at constant $\beta$ in zero-free disk; proves Lee–Yang zero-free region. **Kills SYK as constant-T quantum-advantage candidate** for local-expectation tasks at sufficiently large constant $\beta$.

*Cell verdict*: split. Bergamaschi exhibits a designer Hamiltonian with proven constant-T quantum advantage. Zlokapa removes SYK from the candidate list (for local expectations at sufficient constant $\beta$). The open question is whether *natural physical noncommuting Hamiltonians* at constant $\beta$ inherit BQP-hardness — chapter's central open question.

---

### ground-state × stoquastic / frustration-free × k-local × spin / boson

**3 classical reviews + 4 quantum.**

Classical (structural baseline): bravyi-divincenzo-oliveira-terhal-2008 (LH-MIN stoquastic $\in$ AM); bravyi-terhal-2010 (SFF + gap classically poly-time); bravyi-gosset-2017 (FPRAS for ferromagnetic XY/Heisenberg).

Quantum (cooling / ITE / QITE): motta-2020-qite ($N_q = k(2C)^d \ln^d(\cdot)$; runtime $T = mn \cdot e^{O(N_q)}$); chowdhury-somma-2016-gibbs-hitting ($\widetilde O(\sqrt{N/Z}\cdot\mathrm{polylog}(1/\varepsilon))$); poulin-wocjan-2009-thermal-sampling ($\sqrt{D/Z(\beta)}\cdot\mathrm{poly}(\beta,1/\varepsilon)$); chiang-wocjan-2010 ($\sqrt{D/Z_\beta}\cdot\mathrm{polylog}(1/\varepsilon)$).

*Cell verdict*: stoquastic + frustration-free + $1/\mathrm{poly}(n)$ gap is *classically tractable* (Bravyi–Terhal); quantum offers Grover/Szegedy $\sqrt{D/Z}$ speedup but no super-polynomial separation in this row.

---

### any-T × commuting × arbitrary-graph × spin (FK / random-cluster)

**2 classical reviews. No quantum counterpart.**

- **edwards-sokal-1988-fk-representation** — structural — non-negative joint $(\sigma,n)$ measure $\Leftrightarrow$ cluster algorithm exists. Same boundary as Bravyi–DiVincenzo–Oliveira–Terhal stoquasticity.
- **guo-jerrum-2018-random-cluster-rapid-mixing** — $\tau_\varepsilon(P_{RC}) \le 8 n^4 m^2 (m \ln(1-p)^{-1} + \ln \varepsilon^{-1})$. First polynomial bound for SW dynamics on ferromagnetic Ising.

*Cell verdict*: classical wins; ferromagnetic Ising × any-T closed. **[OPEN-FRONTIER]** for $q\ge 3$ ferromagnetic Potts at low T — Gore–Jerrum 1999 closes SW; quantum analysis open.

---

### high-T / intermediate-T × commuting × bounded-degree × n/a (deterministic / correlation-decay)

**1 classical review.** weitz-2006-counting-independent-sets-tree-threshold — deterministic FPTAS in $(n/\varepsilon)^{O(\log\Delta)}$ for hardcore at $\lambda < \lambda_c(\Delta)$. Cell admits both deterministic FPTAS (Weitz) and optimal-time FPRAS (CLV21).

---

### any-T × bosonic / infinite-dim   `[OUT-OF-SCOPE: infinite-dim]`

**1 quantum review.** becker-2026-infinite-dim (cross-listed above). **The cell-opening result for bosonic Gibbs sampling.** Truncation $M = \widetilde{O}(\mathrm{polylog}(t\,c\,E_{\mathrm{Gibbs}}\,|A|/\varepsilon))$; single sample $\widetilde{O}((1/\lambda_2)\mathrm{poly}(|A|,\log(c E_{\mathrm{Gibbs}}/\varepsilon)))$.

*Cell verdict*: **structurally promising** but **outside thesis scope** — the thesis works only with finite-dimensional Hilbert spaces (spins, fermions on a lattice with bounded local Fock dim, qubits/qudits). No classical analogue exists for bosonic Gibbs sampling in this generality. Quantitative advantage hinges on gap bounds beyond single-mode and classically-hard target identification.

---

### n/a × n/a (lifting / non-reversible / PDMP — generic targets)

**6 classical reviews + 7 quantum reviews.** Treated in dedicated Section 4.

---

## Section 3 — Cross-cutting paradigms

A second pass over the same papers, grouped by paradigm. For each paradigm with ≥ 2 reviews, the shared idea, papers + headline scaling, regime cells covered and *not* covered, and structural limit.

### `classical-MCMC` (Glauber / Metropolis on lattice spin systems) — 12 reviews

**Shared idea**: single-site reversible Markov chains targeting classical commuting Gibbs measures, analysed via spectral gap, log-Sobolev / modified-log-Sobolev inequalities, conductance, or coupling. Optimal mixing time $O(n\log n)$ (Hayes–Sinclair lower bound) is the target on bounded-degree graphs in spatially-mixing regimes.

**Headline scalings**: see Section 2 cells "high-T × commuting" and "low-T × commuting" — anari-liu-oveisgharan-2020 / chen-liu-vigoda-2021 / chen-eldan-2022 / eldan-koehler-zeitouni-2022 / bauerschmidt-dagallier-2024 / lubetzky-sly-2013 / guo-jerrum-2018 / cuff-et-al-2012 / bovier-den-hollander-2015 / madras-zheng-2003 / woodard-schmidler-huber-2009 / sly-2010.

**Regime cells covered**: high-T × commuting × sparse / dense; intermediate-T × commuting × geometric-local / mean-field; low-T × commuting × geometric-local / mean-field. 

**Cells the paradigm does NOT cover**: any non-commuting Hamiltonian; any cell where SSM / spectral / entropic independence fails (low-T phase coexistence beyond persistence); continuous-state targets (covered by `classical-PDMP`).

**Limit / structural wall**: SSM and spectral / entropic independence are *classical correlation-decay* conditions; constants $\exp(c/\delta)$ at uniqueness boundary, persistence at low T (WSH), Sly hardness past $\lambda_c$. Glauber is exponentially slow at low T on bimodal classical targets (Eyring–Kramers) and torpid for $q\ge 3$ Potts at first-order transitions (Cuff et al.).

### `classical-tensor-network` (DMRG / METTS / PEPS / PEPO) — 5 reviews

**Shared idea**: thermal states $\rho_\beta$ are represented as MPDOs (1D), PEPOs (2D / higher), or sampled as classical product states. Cost depends on bond dimension $D$ needed to capture thermal entanglement.

**Headline scalings**: white-2009-metts ($O(N m_0^3 \beta N_\tau)$ vs $O(N m_0^6 \beta)$); czarnik-dziarmaga-2015 (per-time-step $O(M^3 D^4)$, heuristic); molnar-schuch-verstraete-cirac-2015 ($D = (N/\varepsilon)^{O(\beta)}$ in any $d$; high-T sharpening $D = e^{O(\log_2^2)}$ for $\beta < O(\log_2 N)$); kuwahara-alhambra-anshu-2021 ($\widetilde O(\beta^{2/3})$ thermal area law; 1D quasi-linear $\widetilde O(n)$).

**Regime cells covered**: 1D × non-commuting × any-T (rigorous); 2D × non-commuting × high-T to intermediate-T (heuristic algorithmically; rigorous existence via Molnar et al.).

**Cells not covered**: 2D / 3D × non-commuting × low-T (PEPO bond dim $\exp(O(\beta))$ becomes super-polynomial and contraction is $\#$P-hard); frustrated / non-stoquastic.

**Limit / structural wall**: $D = \exp(O(\beta))$ at low T in any $d\ge 2$ — area law forces super-polynomial bond dimension. PEPS/PEPO contraction is $\#$P-hard worst case.

### `classical-cluster` (FK / Swendsen–Wang / random-cluster) — 2 reviews

**Shared idea**: lift to a non-negative joint $(\sigma, n)$ measure where global (cluster) updates exist.

**Papers**: edwards-sokal-1988; guo-jerrum-2018.

**Cells covered**: ferromagnetic-Ising / Potts $q=2$ × any-T × any-graph.

**Cells not covered**: $q\ge 3$ ferromagnetic Potts at first-order transitions (Gore–Jerrum 1999 torpid SW); antiferromagnetic / frustrated systems (signed FK measure).

**Limit / structural wall**: FK joint measure has signs unless interactions are "ferromagnetic" in the abstract $W_b\in[0,1]$ sense — exact classical mirror of QMC sign problem.

### `classical-deterministic` / `cluster-expansion` / `polymer-model` — 5 reviews

**Shared idea**: write $\log Z$ as a convergent abstract polymer / cluster expansion at high $T$ (or contour expansion at low $T$), truncate to logarithmic depth, enumerate connected subgraphs in $\mathrm{poly}(n)$ time (Helmuth–Perkins–Regts / Patel–Regts pipeline).

**Papers**: mann-helmuth-2021 (high-T, FPTAS at $|\beta|\le 1/(e^4\Delta)$); mann-helmuth-2023 (stable low-T, FPTAS in $|V(G)|^{O(1)}(1/\varepsilon)^{O(1)}$); kuwahara-kato-brandao-2020 (high-T CMI clustering); weitz-2006 (hardcore in tree uniqueness, $(n/\varepsilon)^{O(\log\Delta)}$); bakshi-liu-moitra-tang-2024 ($\widetilde O(n^{6+\log\mathfrak d/\log(\beta_c/\beta)})$). Cross-scope: yin-lucas-2023, zlokapa-2026.

**Cells covered**: high-T × non-commuting (BLMT, KKB, MH21); low-T × non-commuting × stable perturbations (MH23); high-T × commuting × tree-uniqueness (Weitz). Cross-scope: high-T sampling (Yin–Lucas), SYK constant-T (Zlokapa).

**Cells not covered**: intermediate-T (between rapid mixing high-T and stable low-T); frustrated / antiferromagnetic without stable ground state; continuous symmetry breaking; gapless / critical regimes.

**Limit / structural wall**: cluster expansion convergence threshold $\beta < 1/(e^4\Delta)$ or $1/(8e^3 k)$ is far below physical $\beta_c$; Pirogov–Sinai constants blow up at $\lambda \to \lambda_\star$ or $\beta \to \beta_\star$.

### `classical-QMC` (CT-QMC / directed-loop SSE / worm) — 4 reviews

**Shared idea**: sample the perturbation expansion of $e^{-\beta H}$ (continuous time, no Trotter error) by Metropolis updates that insert/remove vertices on imaginary-time strings.

**Papers**: gull-et-al-2011; syljuasen-sandvik-2002; prokofev-svistunov-tupitsyn-1998; troyer-wiese-2005 (hardness).

**Cells covered**: noncommuting × all-T × stoquastic-in-some-basis (bipartite half-filled Hubbard, ferromagnetic Heisenberg, Bose–Hubbard, sign-free $XXZ$).

**Cells not covered**: non-stoquastic in every basis (frustrated AFMs, fermions away from half-filling, doped Hubbard, real-time Keldysh).

**Limit / structural wall**: the sign problem. $\langle s\rangle \sim e^{-\beta N \Delta f}$ kills classical QMC outside sign-free corner; Troyer–Wiese makes worst-case sign-curing NP-hard.

### `classical-PDMP` / `classical-lifted` (Zig-Zag, BPS, ECMC, lifted MH, DHN) — 6 reviews

See dedicated Section 4 on the lifting thread.

### `structural` / `complexity-theoretic` / `lower-bound` / `hardness` (classical) — 7+ reviews

**Shared idea**: classify Hamiltonians or sampling targets by intrinsic structural properties and prove worst-case complexity barriers conditional on standard assumptions.

**Papers**: bravyi-divincenzo-oliveira-terhal-2008; bravyi-terhal-2010; troyer-wiese-2005; klassen-marvian-2020; hangleiter-roth-nagaj-eisert-2020; sly-2010; bakshi-liu-moitra-tang-2024; kuwahara-kato-brandao-2020; edwards-sokal-1988.

**Limit / structural wall**: NP$\ne$RP, NP$\not\subseteq$BPP, P$\ne$NP — barriers conditional on standard assumptions. Hardness is *worst-case*; typical-case statements remain mostly open.

### `quantum-Lindblad-KMS` — 23 reviews (spine of the corpus)

**Shared idea**: continuous-time Lindbladian whose stationary state is the Gibbs state, satisfying KMS detailed balance (or its approximation). Operator Fourier transform of jump operators, Gaussian or Metropolis filter, smooth coherent term to enforce DB despite finite energy resolution. CKG (Chen 2023) is the chassis; Chen 2025 makes DB *exact*; DLL 2024 / Gilyén 2024 give discrete-time and rank-preserving variants; Scandi–Alhambra derives it microscopically; Becker 2026 extends to infinite dim.

**Papers with headline scalings**: see Section 2 cell "all-T × noncommuting × local/k-local". Spine includes: chen-2023-thermal-prep, chen-2025-exact-noncommutative, chen-2026-lamb-shift, chen-rouze-2025-locally-markovian, ding-2024-kms-samplers, ding-2024-low-temp, ding-2025-end-to-end, gilyen-2024-glauber-metropolis, rouze-2024-thermal-universal, rouze-2024-optimal-gibbs, scandi-alhambra-2025-kms-thermalization, ramkumar-soleimanifar-2024-sparse-mixing, tong-zhan-2024-fermionic-mixing, smid-2025-fermi-hubbard, becker-2026-infinite-dim, fang-2026-detectability, leng-2026-without-walks, jiang-2026-single-trajectory, slezak-2026-system-bath, lin-2025-dissipative-preparation.

**Regime cells covered** (settled with rigorous poly bounds):
- high-T noncommuting local (Rouzé 2024 ×2, $\widetilde{O}(\log n)$);
- 1D commuting all-T (Bardet 2023, Kochanowski 2024 $O(\log|\Lambda|)$);
- high-T commuting NN ($\widetilde{O}(\log N)$, Capel 2021);
- weakly-interacting fermionic any-$T$ (Tong–Zhan, Smid; $\widetilde{O}(n^3)$ on Hubbard any $D$);
- random sparse $\beta\|H\|=O(1)$ (Ramkumar–Soleimanifar polylog);
- 2D toric code low-T with logical jumps (Ding 2024);
- constant-T BQP-encoded designer Hamiltonians (Bergamaschi);
- bosonic single-mode photon-number (Becker 2026).

**Cells NOT covered** (no rigorous poly bound): low-T noncommuting generic local (chapter punchline cell); spin glasses (random or structured); non-stoquastic with sign problem at arbitrary $\beta$; gapless critical phases; fermionic strong coupling; topologically ordered models beyond toric code with logical jumps.

**Limit / structural wall**:
- *Low-T bottleneck barrier*: Gamarnik 2024 and Rakovszky 2024 prove $T_{\mathrm{mix}} = \exp[\Omega(\mathrm{poly}\,n)]$ for any local KMS-DB sampler on SSB phases, classical glasses, good LDPC codes, 2D TFIM ferromagnet — locality is precisely what makes the sampler efficient to *implement*, and is what kills it at low T.
- *Square-root ceiling on lifting*: Li–Lu / Lu prove $\sqrt{\lambda_O}$ ceiling on any second-order lift of a KMS-DB Lindbladian — matches CLP classical ceiling; no extra quantum gain from non-reversibility alone.
- *Quasi-locality bottleneck on Szegedy walks*: closed by Leng 2026 (SOS factorization, no Davies-only / no rounding promise) and Fang 2026 (DL for commuting parent). Noncommuting parent Hamiltonians of CKG/DLL still open for direct quadratic speedup.

### `quantum-Lindblad-other` (Davies / hypocoercive / MLSI side) — 12 reviews

**Shared idea**: not strictly KMS-DB OFT, but Davies generator (commuting case), Hyper- or Poincaré-style inequality, MLSI/qLSI, hypocoercivity.

**Papers**: kastoryano-temme-2013-qLSI; temme-2010-chi2-divergence; kastoryano-brandao-2016-commuting-case; capel-2021-modified-log-sobolev; bardet-2023-rapid-thermalization-spin-chains; kochanowski-2024-rapid-thermalization; capel-2025-noncommutative-OT; fang-2024-hypocoercivity; li-lu-2025-spacetime-poincare; lucia-2015-rapid-mixing-stability.

**Cells covered**: commuting Hamiltonian regimes; cells where Davies is implementable.

**Cells not covered**: noncommuting Davies (non-local, exponentially expensive); low-T commuting in $d\ge 2$ generically; topologically ordered.

**Limit / structural wall**: Davies generator requires commutativity for locality; noncommuting Davies needs exponential support of jumps. Hypocoercivity gives $L^2$ decay only — no entropy version yet.

### `quantum-dissipative-other` (system-bath / collision / quasiparticle) — 11 reviews

**Shared idea**: not (necessarily) KMS-DB Lindbladian directly, but realistic physical system-bath dynamics + ancilla reset, with the engineered Lindbladian as the second-order effective generator.

**Papers**: hahn-2025-local-driving; lloyd-abanin-2025-near-term; ong-2026-weak-coupling-bounds; hagan-wiebe-2025-thermo-cost; slezak-2026-system-bath; lloyd-2025-quasiparticle-cooling; chen-2026-lamb-shift; ding-2025-end-to-end; zhang-2023-dissipative-gibbs.

**Cells covered**: same as KMS-DB (with realistic-hardware focus). Slezak bridges *physical* dynamics to *algorithmic* mixing-time guarantees.

**Limit / structural wall**: same as KMS-DB. Hardware-realisability advantage costs polynomial factors in $\tau_{\mathrm{mix}}$ and $1/\varepsilon$.

### `quantum-walk` / `quantum-block-encoding` / `quantum-phase-estimation-based` (coherent-acceleration) — 8 reviews

**Shared idea**: replace incoherent mixing by a coherent Szegedy walk / QSVT eigenstate projection, gaining $\sqrt{\lambda_{\mathrm{gap}}}$ over $\lambda_{\mathrm{gap}}$.

**Papers**: wocjan-temme-2023-szegedy-walks; temme-2011-quantum-metropolis; poulin-wocjan-2009-thermal-sampling; chiang-wocjan-2010-thermal-gibbs; chowdhury-somma-2016-gibbs-hitting; rall-2023-rounding-promises; arunachalam-2022-simpler-faster; leng-2026-without-walks; fang-2026-detectability.

**Cells covered**: classical-Hamiltonian Gibbs (Arunachalam); Davies-generators with rounding promise (WT23); commuting parent Hamiltonians (Fang 2026); arbitrary KMS-DB (Leng 2026).

**Cells not covered**: noncommuting parent Hamiltonians of CKG/DLL (Fang 2026 explicit open problem); generic $H$ without rounding promise.

**Limit / structural wall**: rounding promise unphysical for generic $H$. Boosted shift-invariant in-place QPE *impossible* (Chen 2023 Appendix H).

### `lower-bound` / `no-go` — 5+ reviews

- gamarnik-2024-slow-mixing: $T_{\mathrm{mix}} = \exp(\Omega(\mathrm{poly}\,n))$ for local KMS-DB on classical-slow / stabilizer-code / 2D TFIM-ferromagnet.
- rakovszky-2024-bottlenecks: quantum Cheeger bottleneck — local sampler can't beat free-energy barriers.
- ilin-lychkovskiy-2021-quantum-speed-limit: temperature-aware MT-style QSL for unitary drivings.
- (cross-scope) zlokapa-2026-syk-classical: classical upper bound = quantum no-go for SYK constant-T.
- (cross-scope) yin-lucas-2023-classical-high-T: rules out MBQC speedups at high T.

---

## Section 4 — Non-reversible / lifted Gibbs samplers (dedicated subsection)

This is a thesis subsection on its own. The user wants the lifting / non-reversible thread treated as an independent research direction, with **explicit pairing of classical lifts with their quantum analogues** so the chapter writer can present a clean classical↔quantum table.

### Foundational ceiling

**Chen, Lovász & Pak (1999)** (`chen-lovasz-pak-1999-lifting-markov-chains`): For any reversible Markov chain $M$ with stationary distribution $\pi$, $\pi_0 = \min_i \pi_i$, define $\mathcal C$ = minimum congestion-and-local-cost of a unit multicommodity flow. Then:

- *Set-time / hitting bound*: $\widehat{\mathcal A} \ge \tfrac{\sqrt 5}{40}\sqrt{\mathcal A}$ and $\widehat{\mathcal H} \ge \tfrac{1}{10\sqrt{30}} \sqrt{\mathcal H/\log(1/\pi_0)}$ for any lifting (Theorem 3.1).
- *Optimal-lifting characterisation*: $\tfrac{1}{2}\mathcal C \le \inf_{\text{lifts}} \widehat{\mathcal H} \le 144\,\mathcal C$ (Theorem 3.2). **Square-root speedup ceiling.**
- *Reversible-lifting no-go (Cor. 3.4)*: if $\widehat M$ is reversible, $\widehat{\mathcal H} \ge \mathcal H/(128\log(1/\pi_0))$. *All non-trivial gain requires non-reversibility.*

**Hilbert-space generator-level ceiling** (Lu 2025 + Li–Lu 2025): the **same** $\Theta(\sqrt{\lambda_O})$ ceiling holds for any second-order lift of any KMS-detailed-balanced quantum Markov semigroup.
- **Lu 2025 (`lu-2025-hypocoercivity-lifting`) Theorem 4.5**: $\nu \le (1+\log C)\sqrt{\widetilde{s}_m^{-1}\lambda_O}$ via singular-value-gap-vs-spectral-gap.
- **Lu 2025 Theorem 4.10**: matching flow-Poincaré lower bound at optimal $\gamma_{\max} = \sqrt{C_{1,T}/C_{2,T}}$.
- **Li–Lu 2025 (`li-lu-2025-quantum-lifting`) Theorem 2.16**: $\nu(\mathcal{L}_\gamma) = O(\sqrt{\lambda_O})$ for any second-order lift of a KMS-DB QMS.

### Classical-quantum lift pairing (cleaned for the chapter)

Each classical PDMP / lifted sampler is paired with its quantum analogue, since the user wants this table to read cleanly.

| Classical lift | Quantum analogue | Pairing reference |
|---|---|---|
| **Diaconis–Holmes–Neal (2000)** — cycle $\mathbb{Z}_n$, $\Theta(n)$ vs $\Theta(n^2)$ | **Apers et al. (2018)** — explicit LMC on $\mathbb{Z}_N^d$ giving $\tau(\varepsilon) = O(Md\ln d\ln(1/\varepsilon))$ matching best QW lattice mixing; **Li–Lu 2025 §4.2** — quantum lift of classical random walk on a chain achieving same $\Theta(\sqrt{\lambda_Q})$ rate | DHN saturates CLP ceiling on diameter; Li–Lu 2025 §4.2 constructs the *quantum* lift on a chain reaching the same ceiling |
| **Bierkens–Fearnhead–Roberts (2019)** — Zig-Zag PDMP, $O(1)$ per ESS after $O(n)$ pre-processing | **Claudon, Piquemal & Monmarché (2025)** — quantum walk on flat discriminant of nonreversible Markov kernel, governed by geometric reversibilization gap $\gamma(Q)$; *potentially exponential* speedup over $P$ when $1-\langle\pi\|D^j\|\pi\rangle\ll\gamma(Q)$ | ZZ exploits non-reversibility to beat $\Omega(n)$ floor; Claudon shows the Szegedy quadratic-mixing speedup extends to nonreversible chains by quantizing the flat discriminant |
| **Bouchard-Côté–Vollmer–Doucet (2018)** — BPS, ESS/CPU $\sim d^{-1.47}$ on isotropic Gaussian, refresh rate required | **Fang, Lu, Tong (2024)** — quantum hypocoercivity via twisted Lyapunov; exponential GNS-norm decay even when $\mathcal{L}^D$ has macroscopic kernel | BPS uses refresh to escape non-ergodicity; FLT24 uses twisted Lyapunov as the quantum analogue of velocity-refresh hypocoercivity |
| **Michel, Kapfer & Krauth (2014)** — event-chain Monte Carlo on continuum particles | **(No direct review)** — `[GAP]` no quantum analogue of ECMC has been formulated for discrete-spin Gibbs targets | Discrete-spin / Hamiltonian-Gibbs analogue of PDMP is open |
| **Andrieu–Livingstone (2021)** — Peskun–Tierney variance ordering for $(\mu,Q)$-reversible kernels | **Lu 2025 / Li–Lu 2025 / Scandi–Alhambra 2025** — KMS detailed balance $\mathcal L^* = \Gamma_\beta \mathcal L \Gamma_\beta^{-1}$ where modular involution $\Gamma_\beta$ plays role of $Q$ | The $(\mu,Q)$-reversibility framework and quantum KMS-DB share the same structural device: non-self-adjoint generator with correct fixed point via an involution |
| **Chen–Lovász–Pak (1999) ceiling** — $\inf_{\mathrm{lifts}} H_{\mathrm{mix}} = \Theta(\mathcal C)$ for classical chains | **Apers, Sarlette & Ticozzi (2018)** + **Apers, Ticozzi, Sarlette (2017)** — any $\bar p$-invariant local QW simulable by classical LMC with same mixing time up to $\log(1/\varepsilon)$; conductance bound $\tau \ge 1/(8\Phi_{\bar p})$ binds quantum walks just as CLP binds classical lifts | Quantum walks **are** classical lifts in the CLP sense |
| **(Empirical) lifting via cluster updates** — various / not in corpus | **Apers (2017) `apers-2017-lifting-mix-faster`** — 5-bit taxonomy; lifts gain only when one of (s) or (i) is dropped; under (SI) and absence of (re), $\tau_M \le D_\mathcal{G}+1$ | CLP no-go can be *circumvented* if invariance is relaxed during transient |

### The load-bearing thesis claim

The KMS-DB Lindbladian is itself **a lift** in the formal sense of Lu 2025 Definition 4.4: the dissipative part $\mathcal{L}_s$ (symmetric in KMS inner product) plays the role of the "Ornstein–Uhlenbeck dissipation in momentum"; the coherent Lamb-shift / Hamiltonian drift $\mathcal{L}_a = i[H+G,\cdot]$ plays the role of the "Liouville transport". Both have $\ker\mathcal{L}_s \supsetneq \ker\mathcal{L}$ — the hallmark of hypocoercive lifting.

Three reasons this is a dedicated thesis subsection:

1. **Paradigm-independent ceiling.** CLP applies to any reversible base chain regardless of Hamiltonian regime. Lu 2025 Theorem 4.5 extends the same ceiling to *quantum* Lindbladians. It is a structural ceiling that cuts across every cell of the regime grid.

2. **Quantum-walk bridge.** Apers–Sarlette–Ticozzi (2018) / Apers–Ticozzi–Sarlette (2017) / Claudon: quantum walks **are** classical lifts in the CLP sense — same $\sqrt{\Phi^{-1}}$ structural coincidence, *also* on the quantum-walk side. The KMS-DB system+bath dilation is morally a Hilbert-space analogue of an `(S)`-style classical lift with structured initialization in the bath.

3. **Where quantum can break the ceiling.** Lu 2025 explicitly states (§4.4): the lifting framework gives **no separation** between classical and quantum on the ceiling axis. Quantum advantage from non-reversibility *alone* is bounded. **Quantum advantage in Gibbs sampling must come from outside the lifting framework** — faster overdamped base (Hamiltonian simulation as a primitive that has no classical analogue), genuinely non-local samplers (Ramkumar–Soleimanifar 1-design jumps escaping the locality bottleneck of Gamarnik / Rakovszky), or coherent primitives like Grover/Szegedy with rounding-promise-replacement gadgets (Leng 2026 SOS, Fang 2026 DL for commuting).

### Open questions specific to the lifting thread

- **Constructive gap**: CLP $\mathcal{C}$ is LP-computable but optimal lift is not LP-constructible for physical Gibbs measures; same on the quantum side (Lu 2025 §5).
- **Beyond KMS-DB**: hypocoercivity-as-lifting requires symmetric/antisymmetric splitting under KMS-DB. Non-DB samplers (FLT24 hypocoercivity drops DB at cost of twisted-Lyapunov machinery) need their own ceiling theorems.
- **Discrete-spin / Hamiltonian-Gibbs analogue of PDMP**: classical Zig-Zag, BPS, event-chain rely on $\Psi \in C^1$ on $\mathbb{R}^d$; **no discrete-spin equivalent exists yet** for Ising / spin-glass / quantum-Gibbs lattice problems. **`[GAP]`** This is the structural gap that prevents direct transfer of the classical PDMP machinery to quantum Hamiltonian Gibbs sampling.
- **Time discretisation** (Lu 2025 §5): the continuous-time $\sqrt{\cdot}$ ceiling at Trotter+GQSP discretisation is open.
- **Quantum analogue of "persistence"** (cleanest open thesis question identified in this synthesis): Woodard–Schmidler–Huber identify a structural property (persistence of a narrow mode at the hottest rung) under which classical PT/ST is unconditionally torpid. Whether KMS Lindbladians inherit a persistence-style obstruction on the same models (mean-field Potts $q\ge 3$ in coexistence interval) is **open**.
- **Where physical KMS-DB samplers sit relative to optimal lifts**: Lu's framework certifies *existence* of optimally-tuned lifts; whether CKG, DLL actually saturate $\Theta(\sqrt{\lambda_O})$ on concrete Hamiltonians (and not just on toy bipartite models) is open.
- **`[GAP]`** Suwa–Todo (2010) is named in `classical-gibbs-sampling-corpus.md` §1 as a Tier-1 entry but **no per-paper review exists** yet. Highest-priority next `/lit-review` target on the classical lifted-MCMC side.

### Why this is a thesis subsection

The CLP ceiling is *meta*: it applies to *any* classical Gibbs sampler — Glauber, Metropolis, cluster, tempering, PDMP. The structural mechanism (auxiliary momentum / direction variable that converts diffusive $\sqrt T$ into ballistic $T$) is paradigm-independent and crosses every cell of the regime grid. Three reasons this thread cannot be folded into a single regime cell:

1. *Paradigm-independent ceiling*: CLP applies to any reversible base chain, regardless of Hamiltonian regime.
2. *Quantum-walk bridge*: Apers et al. and Dervovic show quantum walks *are* classical lifts in the CLP sense.
3. *KMS-DB connection*: the modular involution of KMS-DB Lindbladians is the operator-algebraic incarnation of the $Q$ in Andrieu-Livingstone's $(\mu, Q)$-reversibility. The structural device of "non-self-adjoint generator with correct fixed point via an involution" is shared between classical lifting and quantum KMS detailed balance.

---

## Section 5 — Contact points

This is the heart of what the user asked for. Per the explicit framing: **every honest cross-scope contact, categorised** as `[DIRECT]` (clean axis, verdict possible), `[VERDICT-UNCLEAR]` (same model, incommensurate access / precision / cost metrics), or `[OPEN-FRONTIER]` (neither side meaningful, or only one side attempted).

Within-scope contacts (classical↔classical, quantum↔quantum) are listed in the per-scope syntheses (Section 5 of `literature-review-synthesis-classical.md` and 5a of `literature-review-synthesis-quantum.md`); cross-scope contacts are the focus here.

### 5a — Direct cross-scope comparisons `[DIRECT]`

Pairs where same model + same regime + comparable cost metric gives a clean "X beats Y" or "X ties Y" verdict.

#### Contact: rouze-2024-optimal-gibbs ↔ bakshi-liu-moitra-tang-2024-high-T-unentangled
- *Category*: `[DIRECT]`
- *Common axis*: high-T short-range noncommuting local Hamiltonians, polynomial Gibbs prep, constant-$\beta$ threshold.
- *Scaling A* (RFA 2024 quantum): rapid mixing $t_{\mathrm{mix}} = O(\log(n/\varepsilon))$ for $\beta < 1/(615^D J)$; end-to-end $\widetilde{O}(n)$ Hamiltonian-sim + $\widetilde{O}(n)$ two-qubit gates.
- *Scaling B* (BLMT 2024 classical): $\widetilde{O}(n^{6 + \log\mathfrak{d}/\log(\beta_c/\beta)})$ classical depth-1 prep (for $\beta < 1/(c\mathfrak{d}\mathfrak{K})$ separability threshold).
- *Verdict*: **polynomial parity, quantum faster by polynomial factor** — quantum $\widetilde{O}(n)$ vs classical $\widetilde{O}(n^7)$ at $D=O(1)$. Closes the high-T row to *asymptotic* quantum advantage.
- *Caveat*: thresholds $\beta_c$ (classical) and $\beta^* = 1/(615^D J)$ (quantum) may not coincide; both are $\mathcal{O}(1)$, system-size-independent. Quantum prepares a quantum state, classical prepares a measurement-basis sampler.

#### Contact: smid-2025-fermi-hubbard ↔ gull-et-al-2011-ct-qmc-impurity-models
- *Category*: `[DIRECT]`
- *Common axis*: Fermi-Hubbard / fermionic impurity model finite-$T$ sampling, $D\ge 2$, weak-to-intermediate coupling.
- *Scaling A* (Smid quantum): end-to-end $\widetilde{O}(n^3\,\mathrm{polylog}(1/\varepsilon))$ on Hubbard *any $D$*, any constant $\beta$, $|U|/t \le U_{\max}(\beta)$.
- *Scaling B* (Gull CT-QMC classical): per-sweep $O(\langle k\rangle^3)$ with $\langle k\rangle \sim N\beta U$ (CT-INT) or $\sim N\beta$ (CT-HYB). Polynomial *only* on sign-free corners (single-orbital Anderson, segment CT-HYB at density-density, particle-hole symmetric bipartite Hubbard, half-filling). Off half-filling / cluster DMFT / multi-orbital matrix CT-HYB: cost $\sim e^{2\beta\Delta F}/M$, exponential.
- *Verdict*: **super-polynomial quantum vs classical sign-broken regime**. First rigorous polynomial *quantum* upper bound where CT-QMC is exponential. The corner that matters most physically — doped, frustrated, $D\ge 2$ — is precisely where classical fails and quantum succeeds.
- *Caveat*: Smid's $U_{\max}(\beta)$ is implicit (Hastings stability constant), not quantitatively pinned. The Mott / metal-insulator transition at $U/t \sim 4$–8 in $D=2$ is *outside* the proven regime. Empirically (Smid Fig. 1, $n\le 8$) the gap survives to $U/t \sim 6$; asymptotic in $n$ inconclusive.

#### Contact: tong-zhan-2024-fermionic-mixing ↔ troyer-wiese-2005-sign-problem-hardness
- *Category*: `[DIRECT]` (structurally)
- *Common axis*: noncommuting fermionic Gibbs sampling at constant $\beta$.
- *Scaling A* (Tong–Zhan quantum): $g = \Omega(1)$ on even-parity sector for $|U| < U_\beta = O(1)$; $t_{\mathrm{mix}} = O(n + \log(1/\varepsilon))$, gate complexity $O(n^2 \mathrm{polylog})$.
- *Scaling B* (Troyer–Wiese): $\langle s\rangle = e^{-\beta N \Delta f}$; generic poly-time cure $\Rightarrow$ NP$\subseteq$BPP.
- *Verdict*: **super-polynomial worst-case gap conditional on NP$\not\subseteq$BPP**, in the weakly interacting / parity-preserving fermionic corner. Quantum sidesteps the path-integral expansion entirely.
- *Caveat*: TW is worst-case; not every instance is sign-hard. The structural argument (KMS-DB samples $\rho_\beta$ directly, no path integral) is the conceptual basis.

#### Contact: bergamaschi-2024-constant-T-advantage ↔ bremner-jozsa-shepherd / bremner-montanaro-shepherd (IQP hardness, structural classical)
- *Category*: `[DIRECT]`
- *Common axis*: constant-$T$ classically-hard sampling on a designer commuting local Hamiltonian.
- *Scaling A* (BCL 2024 quantum): $t_{\mathrm{mix}}(\mathcal{L}_{\mathrm{Davies}}) = O(4^\ell \log n)$ with $\ell = O(\log\log n)$; end-to-end quantum prep $O(n\,\mathrm{polylog}(n))$.
- *Scaling B* (IQP / BMS classical, structural): classical sampling intractable under standard quantum-supremacy conjectures (anti-concentration + worst-to-average-case + non-collapse of PH).
- *Verdict*: **super-polynomial in $n$, conjectured** — under standard QS hardness, the first end-to-end constant-T quantum advantage proof.
- *Caveat*: designer (CTH-style) Hamiltonian, not natural. Whether *natural* (Heisenberg, Hubbard, transverse-field-Ising, SYK, ...) Hamiltonians at constant $T$ inherit BQP-hardness is the central open question. Hardness is *conjectural* (not unconditional).

#### Contact: rouze-2024-thermal-universal (low-T = $\beta = \Omega(\log n)$) ↔ classical complexity-theoretic baseline
- *Category*: `[DIRECT]` (encoding-based)
- *Common axis*: $\mathsf{BQP} = \mathsf{AdiabQP} = \mathsf{GibbsQP}$ at $\beta = \Omega(\log n)$.
- *Scaling A* (RFA quantum): adiabatic purification $T_{\mathrm{ad}} = O((\beta n)^3/\varepsilon^2)$.
- *Scaling B* (classical): no efficient classical sampler exists unless BPP=BQP (complexity-theoretic).
- *Verdict*: **super-polynomial / complexity-theoretic gap**. No efficient classical sampler under BPP$\ne$BQP.
- *Caveat*: $\beta = \Omega(\log n)$ is not "low temperature" in the physics sense — it is "any super-constant $\beta$". CTH construction, not natural physical Hamiltonians.

#### Contact: arunachalam-2022-simpler-faster ↔ huber-2015 / svv (classical partition function)
- *Category*: `[DIRECT]`
- *Common axis*: classical-Hamiltonian partition function estimation.
- *Scaling A* (Arunachalam quantum): $\widetilde{O}(\ln|\Omega|\ln n\cdot\varepsilon^{-1})$ qsamples; $\widetilde{O}(q\varepsilon^{-1}\Delta^{-1/2})$ quantum-walk steps.
- *Scaling B* (Huber / SVV classical): $O(\ln|\Omega|\ln^2 n\cdot\varepsilon^{-2})$ Gibbs samples.
- *Verdict*: **clean quadratic speedup** in both $\varepsilon$ and $\Delta$ for the partition function task. Same character as Montanaro 2015.
- *Caveat*: Grover/Szegedy speedup, not a regime opener. Both sides polynomial; quantum saves quadratic.

#### Contact: rouze-2024-optimal-gibbs ↔ yin-lucas-2023-classical-high-T
- *Category*: `[DIRECT]`
- *Common axis*: high-T noncommuting local Hamiltonian — *measurement-outcome distribution* sampling task.
- *Scaling A* (RFA 2024 quantum): $\widetilde{O}(n)$ end-to-end quantum Gibbs state preparation.
- *Scaling B* (Yin–Lucas classical): classical $O(N^{1+c\delta})$ sampler of computational-basis distribution $p(\bm{x})$ at $\beta < \beta_\star = [2e^2\mathfrak{d}(\mathfrak{d}+1)]^{-1}$.
- *Verdict*: **parity** — both polynomial. Classical wins on circuit-model independence. Yin–Lucas Thm 3 *closes* high-T to MBQC-style exponential speedups.
- *Caveat*: A gives the quantum state; B gives only the measurement distribution. Different output objects.

---

### 5b — Verdict-unclear cross-scope contacts `[VERDICT-UNCLEAR]`

Same model class, but different access models / different precision targets / one rigorous and the other empirical / different cost metrics.

#### Contact: chen-2023-thermal-prep ↔ troyer-wiese-2005 / klassen-marvian-2020 (chapter punchline contact)
- *Category*: `[VERDICT-UNCLEAR]`
- *Common axis*: sign-problem hardness — what classical methods cannot do, vs what quantum methods *might* do.
- *Scaling A* (CKG 2023 quantum): poly $t_{\mathrm{mix}}$ quantum *if* it exists; algorithm-level $\widetilde{O}(\beta t_{\mathrm{mix}}^2/\varepsilon)$ per Lindbladian step.
- *Scaling B* (Troyer–Wiese, Klassen–Marvian classical): sign problem NP-hard worst case; LocalSignCure NP-complete with one-local terms.
- *Verdict*: chapter punchline — at fixed Hamiltonian where TW/KM apply (non-stoquastic, low-T, frustrated), there is *no* polynomial classical Gibbs sampler; CKG provides a clean *path* to quantum advantage *if* a polynomial $t_{\mathrm{mix}}$ can be exhibited. **Not settled** because $t_{\mathrm{mix}}$ for KMS-DB on sign-problem-hard $H$ is open.
- *Caveat*: Bottleneck lower bounds (Gamarnik 2024 / Rakovszky 2024) suggest a substantial fraction of such cells *also* trap CKG, leaving the open frontier in non-stoquastic *random* / *structurally averaging* cases.

#### Contact: gamarnik-2024-slow-mixing / rakovszky-2024-bottlenecks ↔ cuff-et-al-2012 / woodard-schmidler-huber-2009 / bovier-den-hollander-2015
- *Category*: `[VERDICT-UNCLEAR]` — same model class (low-T metastability), incommensurate (classical Glauber vs quantum local KMS-DB Lindbladian)
- *Common axis*: low-T torpid mixing on classical Hamiltonians (Ising, Potts, persistent narrow-mode targets, mean-field).
- *Scaling A* (Gamarnik / Rakovszky quantum lower bounds): $T_{\mathrm{mix}} = \exp[\Omega(\mathrm{poly}\,n)]$ for any local KMS-DB sampler. Specifically $e^{\Omega(n)}$ for commuting classical models above $\beta_c$, $e^{\Omega(\sqrt{n})}$ for 2D TFIM in ferromagnetic phase.
- *Scaling B* (Cuff et al., WSH, Bovier–den Hollander classical): $e^{\Theta(n)}$ classical Glauber on mean-field Potts / persistent narrow modes / tempering torpidity; Eyring–Kramers $K e^{\beta\Gamma^*}(1+o(1))$.
- *Verdict*: **quantum inherits classical low-T bottlenecks — no asymptotic gap, no separation**. Classical and quantum both fail.
- *Caveat*: cell is *closed* against locally-implementable samplers (classical or quantum). Quantum advantage in this regime, if it exists, *cannot* come from locally constructed KMS-DB samplers and must invoke either non-local jump operators (Ramkumar–Soleimanifar's 1-design jumps), unitary preparation circuits with global structure, or pre-fragmented sampling within a single phase (Brandão–Kastoryano-style local-circuit preparation of one phase).

#### Contact: bauerschmidt-dagallier-2024 ↔ bergamaschi-2024-constant-T-advantage
- *Category*: `[VERDICT-UNCLEAR]` — *different* problem (natural ferromagnetic Ising vs designer CTH), connected via the high-T-to-constant-T axis
- *Common axis*: intermediate-T Glauber-style dynamics, polynomial mixing.
- *Scaling A* (Bauerschmidt–Dagallier classical): LSI from susceptibility; polynomial $1/(\beta_c-\beta)$ for near-critical $\mathbb Z^d$ ($d\ge 5$).
- *Scaling B* (Bergamaschi 2024 quantum): designer commuting CTH at any $\beta = \Theta(1)$ with conjectured classical hardness.
- *Verdict*: **incommensurate** — natural ferromagnetic Ising classically tractable everywhere by Bauerschmidt; the designer Hamiltonian in Bergamaschi is *not* a natural physical model. The comparison is between "classical solves all natural ferromagnetic Ising in poly time" and "quantum solves some designer Hamiltonian in poly time that classical cannot under QS conjectures". These are different sampling problems.
- *Caveat*: open whether *any natural* intermediate-T noncommuting Hamiltonian inherits constant-T BQP-hardness from Bergamaschi's construction.

#### Contact: kuwahara-alhambra-anshu-2021 ↔ rouze-2024-optimal-gibbs / bardet-2023-rapid-thermalization-spin-chains (1D thermal)
- *Category*: `[VERDICT-UNCLEAR]` — 1D thermal: classical wins by polynomial factor with different output type
- *Common axis*: 1D quantum spin chain Gibbs sampling at $\beta = o(\log n)$.
- *Scaling A* (KAA classical): MPO bond dimension $\chi = \exp(\sqrt{\widetilde O(\beta\log(n/\varepsilon))})$; quasi-linear $\widetilde O(n)$ classical algorithm for MPS representation.
- *Scaling B* (RFA / Bardet quantum): $\widetilde O(\log n)$ rapid mixing on 1D commuting (Bardet $\Omega(1/\log|\Lambda|)$ MLSI) / $\widetilde O(n)$ noncommuting high-T (RFA).
- *Verdict*: **classical $\widetilde O(n)$ slightly slower than quantum $\widetilde O(\log n)$ on 1D**, but classical produces a *tensor network* representation; quantum produces a *quantum state*. Different outputs, polynomial gap.
- *Caveat*: KAA is super-logarithmic in $\beta$; RFA / Bardet have explicit $\beta$-dependence. At $\beta \to \infty$ KAA's bond dimension blows up; whether quantum survives depends on Tong–Zhan / Smid weak-$U$ window or Mann–Helmuth stable-perturbation window.

#### Contact: motta-2020-qite ↔ bravyi-terhal-2010-frustration-free-stoquastic
- *Category*: `[VERDICT-UNCLEAR]`
- *Common axis*: SFF (stoquastic frustration-free) ground state preparation along adiabatic path with $1/\mathrm{poly}(n)$ gap.
- *Scaling A* (Motta QITE): $N_q = k(2C)^d\ln^d(\cdot)$; runtime $T = mn\cdot e^{O(N_q)}$ — quasi-polynomial in $n,m$, exponential in $C^d$ (conditional on $C$ bounded).
- *Scaling B* (Bravyi–Terhal): classical sampler in $\mathrm{poly}(n,\delta^{-1})$ given $\Delta(s) \ge 1/\mathrm{poly}(n)$ along path.
- *Verdict*: **classical wins on rigorous guarantees** ($\mathrm{poly}(n,\delta^{-1})$ vs quasi-poly conditional on $C$). QITE is heuristic; SFF + gap is classically tractable rigorously.
- *Caveat*: QITE is NISQ-friendly; Bravyi–Terhal is rigorous-only — different "axes" (heuristic / experimental vs rigorous algorithm).

#### Contact: ding-2024-low-temp ↔ Glauber on toric code stabilizer ground space
- *Category*: `[VERDICT-UNCLEAR]`
- *Common axis*: 2D toric code low-$T$ Gibbs state preparation.
- *Scaling A* (Ding 2024 quantum, with 4 global logical jumps): $\mathrm{Gap}(-\mathcal{L}_\beta) \ge \max\{e^{-O(\beta)}, \Omega(N^{-3})\}$; $t_{\mathrm{mix}} = O(\beta\mathrm{poly}(N))$.
- *Scaling B* (classical, e.g. local Glauber on stabilizer ground space): Toric code is sign-free in computational basis after standard rotation, so classical QMC sampling is in principle possible — but logical-sector ergodicity is still an obstruction. **No formal classical sampler analysis covers this**.
- *Verdict*: **superpolynomial quantum-vs-quantum improvement** from $e^{\Theta(\beta)}$ (AFH09 local Davies) to $\mathrm{poly}(\beta)$. Classical baseline is *ambiguous*: cluster QMC with a tailored cluster move *might* match.
- *Caveat*: Whether classical cluster QMC can match Ding's logical-sector flip with a *non-local* cluster move is a `[CHECK]` open question.

#### Contact: hahn-2025-local-driving ↔ classical lifted samplers (analog NISQ vs digital MCMC)
- *Category*: `[VERDICT-UNCLEAR]`
- *Common axis*: thermalization with no block-encoding, analog Hamiltonian evolution.
- *Scaling A* (Hahn quantum): $T_{\mathrm{tot}} = \widetilde{O}(\beta\tau_{\mathrm{mix}}^2/\varepsilon^2)$.
- *Scaling B* (classical PDMP, e.g. Bouchard-Côté BPS): ESS/CPU $\sim d^{-1.47}$ on isotropic Gaussian; refresh-rate-dependent.
- *Verdict*: **different access models** — Hahn 2025 is analog quantum hardware; PDMP is classical continuous-state. Cannot directly compare without specifying a common target distribution.
- *Caveat*: makes quantum advantage on classical-statistical-mechanics-style continuum systems *less* plausible — obvious "quantum walk speedup over diffusive Metropolis" largely already eaten by lifting + factorized filter.

#### Contact: lubetzky-sly-2012 (2D Ising critical) ↔ quantum at 2D Ising criticality
- *Category*: `[VERDICT-UNCLEAR]`
- *Common axis*: 2D Ising at $\beta = \beta_c$, polynomial classical mixing.
- *Scaling A* (Lubetzky–Sly classical): $t_{\mathrm{mix}} = O(n^c)$ polynomial for some unspecified $c$ (much larger than the conjectured physical $z\approx 2.17$).
- *Scaling B* (quantum): **no quantum analysis exists for 2D Ising at criticality**. Bardet 2023 is 1D; KKB / Kuwahara–Alhambra rule out the cell asymptotically only in the strict high-T (above $\beta_c$). Generic noncommuting Davies fails at critical (long-range correlations).
- *Verdict*: **incommensurate** — classical is polynomial-but-loose; no quantum upper bound on this cell.
- *Caveat*: critical 2D Ising is sign-free / commuting; quantum samplers reduce to Glauber-style at best.

---

### 5c — Open-frontier contacts `[OPEN-FRONTIER]`

Cells where neither side has a meaningful result, OR cells where only one side has results and the other side has not been attempted.

#### Contact: KMS-DB on natural low-T noncommuting × no efficient classical
- *Category*: `[OPEN-FRONTIER]`
- *Common axis*: natural physical noncommuting Hamiltonian at low T (frustrated AFM in transverse field, doped Hubbard at strong coupling, gapped non-stoquastic spin chains).
- *Classical*: Troyer–Wiese / Klassen–Marvian / Hangleiter close basis-curing / SignEasing. No positive classical sampler for natural non-stoquastic models at low T.
- *Quantum*: CKG / Chen 2025 give the Lindbladian; **no $t_{\mathrm{mix}}$ bound** for natural non-stoquastic instances; Gamarnik / Rakovszky close *local* KMS-DB samplers in cases that have classical bottleneck (suggesting the open frontier is non-stoquastic *random* or *structurally averaging* instances where neither classical bottleneck nor classical sign-cure applies).
- *Verdict*: **neither side has a poly-time algorithm**. Chapter's headline open frontier.
- *Suggested closure*: a $t_{\mathrm{mix}} = \mathrm{poly}(n)$ result for CKG / DLL on a *natural physical* non-stoquastic Hamiltonian, OR a quantum lower bound showing it cannot exist.

#### Contact: good qLDPC at low T × classical vs quantum
- *Category*: `[OPEN-FRONTIER]`
- *Common axis*: good qLDPC / expander code Hamiltonians at low T.
- *Quantum*: Rakovszky 2024 Thm 3 closes — extensive energy barrier $\Omega(n)$ $\Rightarrow$ $t_{\mathrm{mix}} = e^{\Omega(n)}$ for any local KMS-DB sampler. Includes classical / quantum LDPC and expander codes.
- *Classical*: **no rigorous classical baseline** for good qLDPC Gibbs sampling. Locally-implementable classical samplers fail by same bottleneck. **Non-local classical sampling not analysed**.
- *Verdict*: **both sides closed for local samplers; non-local sides on both sides completely open**. The combination "is there ANY algorithm — classical or quantum, local or non-local — for low-T good qLDPC Gibbs sampling" is an open question.
- *Caveat*: this is one of the cleanest "no idea if any of it is possible" cells.

#### Contact: gapless critical phases × everything
- *Category*: `[OPEN-FRONTIER]`
- *Common axis*: $H$ at a gapless critical point (continuous symmetry breaking, conformal critical 2D, frustrated AFM at criticality).
- *Classical*: Mann–Helmuth cluster expansion fails (critical regime is outside polymer convergence radius). Lubetzky–Sly only handles classical Ising (commuting). No analysis for gapless quantum critical.
- *Quantum*: KMS-DB samplers have no rigorous $t_{\mathrm{mix}}$ at criticality. Hypocoercivity (Fang 2024) addresses degenerate dissipator but not critical phases.
- *Verdict*: **no result either side**. Open frontier.

#### Contact: quantum spin glasses × neither side
- *Category*: `[OPEN-FRONTIER]`
- *Common axis*: random non-commuting Hamiltonians (Aharonov–Gershoni–Klein 2024 [AGK24] mean-field, or non-mean-field structured).
- *Quantum*: Gamarnik 2024 explicitly flags this as "Major open direction" — bottleneck lemma requires classical-bottleneck lift via Feynman–Kac, which needs $-h\sum_i X_i$ to give positive transition rate; non-stoquastic non-commuting cases not addressed. No quantum upper bound for KMS-DB on random non-commuting.
- *Classical*: no analysis. Spin-glass-like noncommuting models without OGP-style clustering structure are uncharted territory.
- *Verdict*: **no result either side**. Open frontier.

#### Contact: low-T fermionic strong coupling × neither side
- *Category*: `[OPEN-FRONTIER]`
- *Common axis*: Hubbard at $U/t \gg 1$ (Mott / strong coupling), low T, $D\ge 2$.
- *Classical*: CT-QMC has sign problem off half-filling; cluster DMFT off half-filling exponential.
- *Quantum*: Smid 2025 / Tong–Zhan only cover weak coupling $|U| < U_\beta = O(1)$; strong-coupling $U/t \sim 4$–8 in $D=2$ outside proven regime.
- *Verdict*: **no rigorous result either side** in the genuinely strong-coupling, doped, $D\ge 2$ corner.
- *Suggested*: extension of Smid's atomic-limit + perturbation argument (Theorem III.11) to $t > 0$ would close.

#### Contact: continuous symmetry breaking at low T × neither side
- *Category*: `[OPEN-FRONTIER]`
- *Common axis*: Heisenberg ferromagnet at low T, $U(1)$, $O(n)$ models.
- *Classical*: explicitly excluded by Mann–Helmuth 2023 (finitely many ground states required). No review covers this cell from the classical side. Hwang–Hwang-Ma–Sheu (1993) and Lelièvre–Nier–Pavliotis (2013) Tier-2 in the corpus index but no per-paper reviews.
- *Quantum*: no rigorous bound. Hypocoercivity (Fang 2024, Li–Lu 2025) addresses some degenerate-dissipator cases but no concrete continuous-symmetry-breaking application.
- *Verdict*: **no result either side**. Open frontier.

#### Contact: 3D Ising classical at low T × beyond
- *Category*: `[OPEN-FRONTIER]`
- *Common axis*: 3D Ising at $\beta > \beta_c$, low-T Glauber.
- *Classical*: Bovier–den Hollander Sec. 17.6 identifies $\Gamma^*$ but not the prefactor. Partial result.
- *Quantum*: no analysis.
- *Verdict*: classical partial; quantum nothing.

#### Contact: real-time Keldysh (out-of-equilibrium dynamics) × classical / quantum
- *Category*: `[OPEN-FRONTIER]`
- *Common axis*: real-time evolution / Keldysh sampling (oscillatory weights for fermionic / open systems).
- *Classical*: Gull et al. notes real-time Keldysh is classically exponential due to sign problem.
- *Quantum*: not addressed by any Gibbs sampler in the corpus (the KMS-DB framework targets equilibrium $\rho_\beta$, not Keldysh).
- *Verdict*: **out of scope of the corpus**. Open whether KMS-DB extensions can handle this regime.

#### Contact: $q\ge 3$ ferromagnetic Potts at low T × beyond local samplers
- *Category*: `[OPEN-FRONTIER]`
- *Common axis*: low-T ferromagnetic Potts $q\ge 3$ on lattice.
- *Classical*: Glauber torpid ($e^{\Theta(n)}$ on mean field Potts at first-order); SW also torpid (Gore–Jerrum 1999, in corpus index only). $q = 2$ resolved by Guo–Jerrum.
- *Quantum*: Gamarnik 2024 includes $p$-spin glass and lattice Potts $q\ge 3$ above $\beta_c$ as classical-bottleneck-inherited cases — local KMS-DB exponentially slow.
- *Verdict*: **both sides torpid for local samplers**; non-local quantum samplers (e.g. Ramkumar–Soleimanifar 1-design jumps adapted to Potts) not analysed.
- *Caveat*: candidate for "first quantum advantage from non-locality" but no concrete construction.

#### Contact: stoquastic but non-ferromagnetic / frustrated antiferromagnetic on non-bipartite × beyond
- *Category*: `[OPEN-FRONTIER]`
- *Common axis*: stoquastic AFM on non-bipartite (no sublattice rotation to ferromagnetic) at low T.
- *Classical*: Bravyi–Gosset 2017 covers stoquastic ferromagnetic XY/Heisenberg only; frustrated AFM not covered.
- *Quantum*: no rigorous bound. Smid covers weak interactions only.
- *Verdict*: **no result either side** on genuinely frustrated stoquastic AFM at low T.

#### Contact: 4D toric code / topological order with extensive energy barriers × quantum sampler
- *Category*: `[OPEN-FRONTIER]`
- *Common axis*: topologically ordered models with macroscopic energy barriers (4D toric code, Hong–Guo–Lucas quantum memories).
- *Quantum*: Rakovszky 2024 Eq. 9 closes — extensive energy barrier $\Omega(n)$ $\Rightarrow$ $e^{\Omega(n)}$ for any local sampler. Ding 2024 toric code construction only works for 2D toric code (4 global logical jumps); generalisation to 4D not attempted.
- *Classical*: no classical baseline.
- *Verdict*: **closed against local samplers**. Open whether non-local logical-jump samplers (Ding-2024-style) generalize.

---

### 5d — Within-scope contacts (cross-references only)

For complete within-scope contact points (classical↔classical, quantum↔quantum), see Section 5 of `literature-review-synthesis-classical.md` and Section 5a of `literature-review-synthesis-quantum.md`. Key within-scope contacts not duplicated here:

- *Classical*: weitz-2006 ↔ anari-liu-oveisgharan-2020 ↔ chen-liu-vigoda-2021 (FPTAS evolution at tree-uniqueness boundary); molnar-schuch-verstraete-cirac-2015 ↔ czarnik-dziarmaga-2015 (PEPO existence vs PEPO-based heuristic algorithm).
- *Quantum*: rouze-2024-thermal-universal ↔ rouze-2024-optimal-gibbs (gap-only vs rapid mixing); chen-2023-thermal-prep ↔ chen-2025-exact-noncommutative ↔ ding-2024-kms-samplers (CKG evolution); gamarnik-2024-slow-mixing ↔ rakovszky-2024-bottlenecks (complementary lower bounds); leng-2026-without-walks ↔ fang-2026-detectability ↔ wocjan-temme-2023-szegedy-walks (sqrt-gap amplification evolution); tong-zhan-2024 ↔ smid-2025 (fermionic mixing).

---

## Section 6 — Frontier

Per-cell three lines: best classical scaling / best quantum scaling / gap and verdict.

### Frontier table (cell-by-cell)

| Cell | Best classical | Best quantum | Gap | Verdict |
|---|---|---|---|---|
| high-T × commuting × bounded-degree (lattice) | $O(n\log n)$ Glauber (Chen-Liu-Vigoda 2021) | $\widetilde O(\log N)$ Davies (Capel 2021, Kochanowski 2024) | polynomial | quantum advantage **ruled out** asymptotically — classical optimal |
| high-T × noncommuting × geometric-local × spin | $\widetilde O(n^7)$ (BLMT 2024); $O(N^{1+c\delta})$ sampling (Yin-Lucas) | $\widetilde O(\log n)$ rapid (RFA 2024 optimal) | polynomial | quantum advantage **ruled out** — polynomial parity |
| intermediate-T × commuting × geometric-local | poly classical (Lubetzky-Sly 2012, Bauerschmidt 2024) | poly Davies / MLSI | polynomial | quantum advantage **ruled out** for natural models |
| intermediate-T (constant-T) × commuting × designer CTH | classically intractable conjecturally (IQP hardness) | $O(n\,\mathrm{polylog}(n))$ (Bergamaschi 2024) | super-polynomial conjectured | quantum advantage **proven conditionally** — designer Hamiltonian only |
| intermediate-T × noncommuting × generic / SYK | $n^{O(\log)}$ (Zlokapa 2026 for SYK at constant $\beta$ disk) | no $t_{\mathrm{mix}}$ bound | unknown | **VERDICT-UNCLEAR**; SYK constant-$\beta$ killed by Zlokapa; *natural noncommuting* at constant-$T$ open |
| low-T × commuting × geometric-local (lattice Ising etc.) | $e^{\Theta(\beta\Gamma^*)}$ (Bovier–den Hollander); $e^{\Omega(n)}$ classical Glauber on persistence | $e^{\Omega(n)}$ quantum (Gamarnik 2024 / Rakovszky 2024 for local KMS-DB) | no asymptotic gap | quantum advantage **ruled out** for local samplers |
| low-T × commuting × topological order (toric code) | (no formal analysis) | $\mathrm{poly}(N,\beta)$ via 4 global logical jumps (Ding 2024) | super-polynomial in $\beta$ | quantum advantage **proven** for the 4 global logical jump scheme; conditional on hardware implementation |
| low-T × noncommuting × generic / natural | hardness barriers (Troyer-Wiese, Klassen-Marvian, Hangleiter); **no positive classical algorithm** | local KMS-DB closed by Gamarnik / Rakovszky; **no positive quantum algorithm on a natural physical Hamiltonian** | unknown | **OPEN-FRONTIER** — chapter's headline open frontier |
| low-T × noncommuting × fermionic × weak coupling (Hubbard) | exponential CT-QMC in sign-broken regime (Gull et al., Troyer-Wiese) | $\widetilde O(n^3)$ end-to-end (Smid 2025), any $D$, any constant $\beta$, weak $U$ | super-polynomial | quantum advantage **conjectured / plausible from corpus** — first concrete fermionic instance |
| low-T × noncommuting × fermionic × strong coupling (Mott Hubbard, $U/t\gg 1$) | exponential (sign problem) | not covered | unknown | **OPEN-FRONTIER** |
| low-T × noncommuting × bosonic / infinite-dim  `[OUT-OF-SCOPE: infinite-dim]` | (no quantum-target classical baseline) | $\widetilde O((1/\lambda_2)\mathrm{poly}(|A|))$ Becker 2026 | open (gap only for single-mode photon-number) | mentioned only; **outside thesis scope** (thesis is finite-dim only) |
| low-T × noncommuting × stable perturbation | $|V(G)|^{O(1)}(1/\varepsilon)^{O(1)}$ FPTAS (Mann-Helmuth 2023) | (no quantum bound in stability region; quantum could match) | polynomial parity | quantum advantage **ruled out** in stability region |
| low-T × commuting × good qLDPC / expander code | (no analysis) | $e^{\Omega(n)}$ closed by Rakovszky 2024 Thm 3 for local KMS-DB | unknown | **OPEN-FRONTIER** — no algorithm classical or quantum (local) |
| low-T × commuting × $q\ge 3$ ferromagnetic Potts | $e^{\Theta(n)}$ Glauber (Cuff et al. 2012); Gore-Jerrum SW torpid | $e^{\Omega(n)}$ closed by Gamarnik 2024 for local KMS-DB | no asymptotic gap | **OPEN-FRONTIER** for non-local samplers either side |
| ground-state × stoquastic / SFF | $\mathrm{poly}(n,\delta^{-1})$ (Bravyi-Terhal 2010) | $\sqrt{D/Z}\mathrm{polylog}$ (Poulin-Wocjan 2009, Chiang-Wocjan 2010, Chowdhury-Somma 2016) | polynomial (Grover) | quantum advantage **ruled out** in $n$, Grover speedup remains |
| ground-state × frustrated / non-stoquastic | (no efficient classical) | conditional on $\beta = \mathrm{poly}(n)$ RFA / CTH-encoded | super-polynomial conjectured | candidate for quantum advantage |
| any-T × bosonic / Schrödinger / continuous variable  `[OUT-OF-SCOPE: infinite-dim]` | (no general efficient classical) | Becker 2026 (single-mode photon-number proven; many-mode open) | unknown | mentioned only; **outside thesis scope** (thesis is finite-dim only) |
| gapless critical phases (continuous symmetry breaking, conformal) | (no analysis) | (no analysis) | unknown | **OPEN-FRONTIER** — nothing either side |
| quantum spin glasses (random noncommuting) | (no analysis) | (no analysis; Gamarnik 2024 "major open direction") | unknown | **OPEN-FRONTIER** — nothing either side |
| real-time Keldysh / out-of-equilibrium | exponential (sign) | not addressed | unknown | **OPEN-FRONTIER** — out of corpus scope |

### Frontier ranking — top-5 cells where quantum could first beat classical

In priority order, the regime cells the corpus most strongly suggests as the first place quantum could beat classical:

1. **low-T × noncommuting × fermionic × weak coupling Hubbard (any $D$)** — **strongest concrete candidate proven**.
   - Justification: Smid 2025 gives end-to-end $\widetilde{O}(n^3\mathrm{polylog})$ on Hubbard *any $D$*, any constant $\beta$, weak $U$ — in the sign-broken regime that defeats CT-QMC / cluster DMFT. Troyer–Wiese makes worst-case sign-curing NP-hard. No rigorous polynomial classical algorithm known. The corner that matters most physically (doped, frustrated, $D\ge 2$) is in the proven Smid regime up to $U_{\max}(\beta)$.
   - Caveat: $U_{\max}(\beta)$ is implicit; numerics extend to $U/t \sim 6$ at $n\le 8$. Strong coupling open.

2. **constant-T × commuting × designer CTH (Bergamaschi 2024 family)** — **first end-to-end constant-T proof**.
   - Justification: $O(n\,\mathrm{polylog}(n))$ quantum prep + classical hardness under standard QS conjectures. First and only proof that *some* class of physical-thermal-Lindbladians at $\beta = \Theta(1)$ is classically intractable.
   - Caveat: designer Hamiltonian, not natural. The chapter punchline cell — whether natural Hamiltonians inherit any of this BQP-hardness — is **open**.

3. **non-stoquastic × low-T × frustrated / sign-bearing — chapter's headline open frontier**.
   - Justification: Classical sign-cure NP-hard (Klassen–Marvian); SignEasing NP-complete (Hangleiter); worst-case Troyer–Wiese. CKG / Chen 2025 / Ding give KMS-DB Lindbladian for any noncommuting $H$, no path-integral expansion. **But $t_{\mathrm{mix}}$ for KMS-DB on natural non-stoquastic Hamiltonians at low T is open**, and Gamarnik / Rakovszky lower bounds suggest most natural cases are *also* trapped quantum-side. Frontier is the *random* / *structurally averaging* non-stoquastic models where neither classical bottleneck nor classical sign-cure applies.
   - Caveat: not yet exhibited — only structural barriers in place.

4. **all-T × random sparse $\beta\|H\|=O(1)$ noncommuting (Ramkumar–Soleimanifar 2024)** — **average-case opener**.
   - Justification: polylog$(N/\varepsilon)$ mixing for the random-sparse-Hamiltonian model at $\beta\|H\| = O(1)$, structural average-case scaling that no classical algorithm matches in the same regime — the corresponding classical baseline is the cluster-expansion line (Mann–Helmuth 2021 high-T) which loses polylog scaling on random non-stoquastic instances. Genuine finite-dim spin / qubit setting, fits the thesis machinery directly.
   - Caveat: random-instance promise, not worst-case; the regime $\beta\|H\| = O(1)$ is the high-T side of the boundary, and the polylog is in $N$ at fixed $\beta$.

5. **low-T × commuting × topologically ordered with extensive logical jumps (Ding 2024 toric code)** — **proof of principle**.
   - Justification: $\mathrm{poly}(N,\beta)$ via 4 global logical jumps, $\beta$-independent gap on logical sector. Sidesteps the AFH09 $e^{\Theta(\beta)}$ wall for *local* Davies. No genuine classical analogue (no classical "hidden" logical degree of freedom on toric code).
   - Caveat: tailored to toric code's exact logical-vs-syndrome decomposition; whether the template generalises to good qLDPC / 4D toric / general topological order is open.

#### Out-of-scope mention: infinite-dim / bosonic — Becker 2026

The all-T × noncommuting × bosonic / infinite-dim cell (Becker 2026 single-mode photon-number, Bose–Hubbard companion) opens a genuinely new direction with no classical analogue for bosonic Gibbs samplers of noncommuting many-mode systems. **The thesis works only with finite-dimensional Hilbert spaces** — spins, fermions on a lattice with bounded local Fock dim, qubits/qudits — so this cell is **outside thesis scope** [`[OUT-OF-SCOPE: infinite-dim]`]. It is recorded here because it is part of the broader literature, but it is not pursued as a frontier candidate for the chapter.

### Cells where quantum advantage is **ruled out by the corpus**

- **high-T × commuting (lattice)**: spectral-independence / localization-schemes machinery matches Hayes–Sinclair $\Omega(n\log n)$.
- **high-T × noncommuting**: BLMT $\widetilde O(n^7)$ + Yin–Lucas $O(N^{1+c\delta})$ vs quantum $\widetilde O(\log n)$ — polynomial parity.
- **SYK constant-$\beta$ disk**: Zlokapa $n^{O(\log)}$ kills the candidate.
- **low-T × commuting × natural (Ising, Potts, lattice)**: Gamarnik / Rakovszky close any local KMS-DB sampler.
- **ground-state × stoquastic + frustration-free**: Bravyi–Terhal $\mathrm{poly}(n,\delta^{-1})$.
- **low-T × noncommuting × stable perturbation**: Mann–Helmuth 2023 FPTAS in stability region.
- **lifting from non-reversibility alone**: Lu 2025 / Li–Lu 2025 ceiling matches CLP — no quantum gain over best classical lift on the non-reversibility axis.

### Cells where quantum advantage is **fully open** (nothing meaningful either side)

- gapless critical phases (continuous symmetry breaking, conformal critical, frustrated AFM at criticality)
- quantum spin glasses (random non-commuting Hamiltonians)
- strong-coupling fermionic (Mott Hubbard $U/t \gg 1$, $D\ge 2$)
- continuous symmetry breaking at low T (Heisenberg ferromagnet, $U(1)$, $O(n)$)
- good qLDPC / expander codes at low T (closed against local samplers either side; non-local both sides open)
- $q\ge 3$ ferromagnetic Potts at low T beyond local samplers
- 4D toric code / topological order with extensive barriers beyond 2D toric code construction
- real-time Keldysh / out-of-equilibrium

---

## Section 7 — Gaps and open questions (corpus level)

### Empty regime cells

The frontier table in Section 6 lists every cell with results. Cells *fully absent* from the corpus:
- **continuous symmetry breaking at low T** — Mann-Helmuth 2023 excludes; Hwang-Hwang-Ma-Sheu / Lelièvre-Nier-Pavliotis Tier-2 in corpus index but no per-paper reviews. Neither classical nor quantum.
- **gapless critical 2D/3D quantum** — no rigorous analysis either side.
- **quantum spin glasses** — Aharonov-Gershoni-Klein 2024 [AGK24] cited only.
- **3D Ising classical low-T** — only Bovier-den Hollander Sec. 17.6; partial.
- **real-time Keldysh** — out of corpus.
- **fermionic strong coupling at low T** — neither.
- **bosonic interacting many-mode** — Becker single-mode only; Bose-Hubbard companion not reviewed.

### Cells with quantum but no classical, or vice versa

*Quantum upper, no classical baseline*: Becker bosonic infinite-dim (no classical analogue exists in this generality); Ding 2024 toric code logical jumps (no formal classical sampler analysis); Bergamaschi designer CTH (classical conjecturally intractable; no positive classical).

*Classical upper, no quantum upper*: Yin-Lucas sampling at $\beta < \beta_\star$ (quantum gives full state via RFA, but no quantum sampler explicitly compared at the same task); Zlokapa SYK constant-$\beta$ disk (quantum upper bound unknown; KMS-DB $t_{\mathrm{mix}}$ on SYK open); Lubetzky-Sly 2D Ising at criticality (quantum analysis absent).

### `[CHECK]` flags accumulated

Combined from per-scope syntheses; for the full audit list see Section 6 (classical) and Section 7 (quantum) of the respective per-scope syntheses. Summary count: **28 `[CHECK]` flags in classical scope; 80+ `[CHECK]` flags in quantum scope; 5–10 cross-scope ambiguities (paradigm-tag routing decisions for Apers 2017/2018, Claudon, Li-Lu, Yin-Lucas, Zlokapa)**.

### `[GAP]` flags emitted

`[GAP]` papers cited 2+ times in the corpus but lacking per-paper reviews:

1. **Suwa-Todo 2010** (Markov chain without detailed balance, PRL 105:120603) — classical lifted-MH machinery. Highest priority classical-side next `/lit-review`.
2. **Apers-Sarlette-Ticozzi 2018** — has per-paper review (in quantum directory). Listed in `apers-2018-quantum-walks-classical-lift.md`. Resolved.
3. **Davies 1974/1976** (Markovian master equation derivation) — foundational Davies-generator paper. Highest priority quantum-side next `/lit-review`.
4. **Chen-Brandão-Gilyén 2021 (CKBG21)** — early KMS-DB Lindbladian Chen 2023 supersedes. Cited 5+ reviews.
5. **Anshu-Arunachalam-Kuwahara-Soleimanifar 2024** (sample-efficient learning of local Hamiltonians from high-T Gibbs). Cited in classical and quantum syntheses.
6. **Bravyi-Hastings 2017** (efficient algorithm for short-range Gibbs states) — Hastings stability theorem reference. Cited in Becker, Tong-Zhan.
7. **Harrow-Mehraban-Soleimanifar 2020** (quasi-polynomial quantum partition function) — predecessor classical-quantum bridge. Cited in 4+ reviews; would belong in both syntheses.
8. **Yung-Aspuru-Guzik 2012** (quantum-quantum Metropolis, QQMA) — degenerate-spectrum loophole paper. Cited in Wocjan-Temme 2023, Temme 2011.
9. **Aharonov-Chen-Kothari-Klyachko 2025 (ACKK25)** — circuit-complexity lower bound for SYK Gibbs preparation. Cited in Zlokapa.
10. **Anshu-Arad-Vidick 2016** (detectability lemma simple proof) — cited in Fang 2026.
11. **Bilgin-Boixo 2010** (commuting thermal state preparation) — historical commuting-case predecessor.
12. **Bergamaschi-Chen-Vazirani 2025 (BCV25)** — metastable approximate-Markov framework. Cited in Chen 2026 strong-Markov.
13. **Aharonov-Gershoni-Klein 2024 (AGK24)** — quantum spin glasses. Cited in Gamarnik.
14. **Hwang-Hwang-Ma-Sheu 1993** + **Lelièvre-Nier-Pavliotis 2013** — continuous-state analogues of CLP $\sqrt{\mathrm{gap}}$.
15. **Suzuki-Trotter classical analogues** — referenced in cluster expansion context.

### Conjectures presented as proofs / heuristics

- **Lloyd-Abanin 2025**: $O(\theta^2)$ fixed-point bound is *perturbative*; rigorous only at 4th order.
- **Lloyd 2025 (quasiparticle cooling)**: kinetic rate-equation is heuristic.
- **Motta 2020 (QITE)**: rigorous correctness *conditional on uniform $C$ on every intermediate iterate* — fails on GHZ-projector example.
- **Fang 2024 hypocoercivity §5**: conjectured $1/\sqrt{\lambda_M}$ Grover-like speedup, unproven.
- **Zlokapa 2026**: extends only to constant $\beta$ disk; large-$q$ all-$\beta$ extension conjectured from non-rigorous ZK26 replica calculation.
- **Czarnik-Dziarmaga 2015**: heuristic algorithm without rigorous bond-dimension theorem (properly marked `algorithm-heuristic` in frontmatter).

### Suggested next `/lit-review` calls (8 strongest gaps)

1. **`/lit-review davies-1976-markovian-master-equation-derivation`** — closes foundational Davies-generator citation chain. Referenced in 8+ quantum reviews.
2. **`/lit-review chen-brandao-gilyen-2021-early-kms-lindbladian`** — clarifies the `CKBG21 → CKG23 → CKG25` chain that anchors the KMS-DB spine.
3. **`/lit-review suwa-todo-2010-markov-chain-without-detailed-balance`** — closes the classical lifted-MH thread.
4. **`/lit-review harrow-mehraban-soleimanifar-2020-quasi-polynomial-quantum-partition`** — quasi-poly classical $Z$ predecessor; central in 4+ quantum reviews; would belong in both syntheses.
5. **`/lit-review yung-aspuru-guzik-2012-quantum-quantum-metropolis`** — QQMA degenerate-spectrum loophole paper.
6. **`/lit-review aharonov-chen-kothari-klyachko-2025-syk-circuit-complexity`** — quantum lower bound pair to Zlokapa 2026.
7. **`/lit-review anshu-arunachalam-kuwahara-soleimanifar-2024-sample-efficient-learning`** — high-T learning axis.
8. **`/lit-review bergamaschi-chen-vazirani-2025-metastable-approximate-markov`** — pairs with Chen 2026 strong-Markov.

---

## Final report (for the caller)

- **Path**: `/Users/bence/code/QuantumFurnace.jl/drafts/literature-review/literature-review-synthesis-all.md`
- **Reviews ingested**: 99 total / 39 classical + 60 quantum / 99 passed `all` scope / 0 skipped.
- **Populated regime cells**: 18 distinct cells (many cross-listed); 9 cells closed against quantum advantage; 5 cells where quantum advantage is proven or conjectured; 8 cells fully open frontier with nothing either side or only one side attempted.
- **Headline frontier candidate cell**: **low-T × noncommuting × fermionic × weak-coupling Hubbard (any $D$)** — strongest concrete candidate proven (Smid 2025 $\widetilde O(n^3)$ vs CT-QMC exponential off half-filling). With **non-stoquastic × low-T × frustrated** as the chapter's headline *open* frontier.
- **Cross-scope contact points emitted by category**: 7 `[DIRECT]`, 9 `[VERDICT-UNCLEAR]`, 11 `[OPEN-FRONTIER]`. (Plus within-scope contacts cross-referenced to per-scope syntheses.)
- **`[CHECK]` flags**: 28 classical + 80+ quantum, full audit lists in per-scope syntheses.
- **`[GAP]` markers**: 15 explicit `[GAP]` items in Section 7; 8 prioritised next `/lit-review` targets listed.
- **Paradigms / reference papers named in corpus index but lacking per-paper reviews**: Davies 1974/1976 (foundational Davies generator); CKBG21 (early KMS-DB); HMS 2020 (quasi-poly $Z$ classical); Suwa-Todo 2010 (classical lifted-MH); Bilgin-Boixo 2010 (commuting predecessor); Anshu-Arunachalam-Kuwahara-Soleimanifar 2024 (sample-efficient learning); ACKK25 (SYK circuit-complexity lower bound); Yung-Aspuru-Guzik 2012 (QQMA); AGK24 (quantum spin glass); BCV25 (metastable Markov framework); Hwang-Hwang-Ma-Sheu 1993 / Lelièvre-Nier-Pavliotis 2013 (continuous-state $\sqrt{\mathrm{gap}}$ analogues); Bravyi-Hastings 2017 (Hastings stability); Anshu-Arad-Vidick 2016 (detectability lemma); plus the Tier-2 classical foundational references named in `classical-gibbs-sampling-corpus.md`.
