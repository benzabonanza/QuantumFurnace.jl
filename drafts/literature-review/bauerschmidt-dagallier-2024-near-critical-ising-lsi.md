---
paper: "Bauerschmidt, Dagallier (2024)"
title: "Log-Sobolev inequality for near critical Ising models"
arxiv: "2202.02301"
year: 2024
venue: "Comm. Pure Appl. Math. 77:2568-2576"
pdf: "supplementary-informations/classical-review/bauerschmidt-dagallier-2024-near-critical-ising-lsi.pdf"

temperature: [high-T, intermediate-T]
commutativity: [commuting]
locality: [geometric-local, long-range]
particle-statistics: [spin]
hamiltonian-models: [Ising]

paradigm: [classical-MCMC]
quantum-or-classical: [classical]

result-type: [mixing-time-upper, gap-bound, structural]
key-scaling: "1/gamma_{beta,h} <= 1/2 + integral_0^beta exp(2 int_0^t chi_s ds) dt; with mean-field chi_beta <= D/(beta_c - beta), gives 1/gamma ~ beta_c/(2D-1) (1-beta/beta_c)^{1-2D} for D > 1/2 — polynomial in distance to criticality and uniform in volume"

related: ["eldan-koehler-zeitouni-2022-spectral-condition-ising", "lubetzky-sly-2012-critical-ising-polynomial", "lubetzky-sly-2013-cutoff-ising-lattice", "guo-jerrum-2018-random-cluster-rapid-mixing"]
---

# Log-Sobolev inequality for near critical Ising models — Bauerschmidt & Dagallier 2024

**One-sentence takeaway**: For ferromagnetic Ising models with bounded coupling spectral radius, the log-Sobolev constant is bounded purely in terms of the (zero-field) susceptibility $\chi_\beta$, giving uniform-in-volume LSI in the entire high-temperature phase and a *polynomial* dependence on $1/(\beta_c - \beta)$ as criticality is approached on $\mathbb{Z}^d$ for $d \ge 5$.

## Setting

- **Sampled object**: ferromagnetic Ising measure $\mu_{\beta,h}(\sigma) \propto \exp(-\tfrac{\beta}{2}\langle \sigma, A\sigma\rangle + \langle h, \sigma\rangle)$ on $\{\pm 1\}^\Lambda$, finite $\Lambda$.
- **Hamiltonian class**: arbitrary symmetric coupling matrix $A$ with $A_{xy} \le 0$ for $x \ne y$ (ferromagnetic), positive definite, spectral radius $\|A\| \le 1$. No geometric/locality restriction; covers $\mathbb{Z}^d$ nearest-neighbour Ising and more general ferromagnets.
- **Temperature regime**: $\beta \ge 0$ up to (and approaching) the critical $\beta_c$. Low-T phase explicitly out of scope.
- **Access model**: continuous-time Glauber / heat-bath / Metropolis dynamics on the spin configuration; classical MCMC.

## Main Results

- **Theorem 1.1 (verbatim)**: $\dfrac{1}{\gamma_{\beta,h}} \le \dfrac{1}{2} + \displaystyle\int_0^\beta e^{2\int_0^t \chi_s\, ds}\, dt$, where $\gamma_{\beta,h}$ is the LSI constant of Glauber dynamics and $\chi_s = \sup_{x \in \Lambda} \sum_{y \in \Lambda} \mathbb{E}_{\beta=s, 0}(\sigma_x \sigma_y)$ is the zero-field susceptibility.
- **High-T corollary (1.7)**: $\dfrac{1}{\gamma_{\beta,h}} \le \dfrac{1}{2} + \beta\, e^{2\beta \chi_\beta}$ whenever $\chi_\beta$ is bounded (entire high-T phase, uniformly in $\Lambda$, with no spatial-mixing input).
- **Near-critical (Cor. 1.2)**: assuming the mean-field bound $\chi_\beta \le D/(\beta_c - \beta)$ with $D > 1/2$, $\dfrac{1}{\gamma_{\beta,h}} \le \dfrac{1}{2} + \dfrac{\beta_c}{2D-1}\big[(1-\beta/\beta_c)^{1-2D} - 1\big] \sim \dfrac{\beta_c}{2D-1}(1-\beta/\beta_c)^{1-2D}$ as $\beta \to \beta_c$. Finite-volume version with $\chi_\beta \le D/(\beta_c - \beta + L^{-2})$ gives $1/\gamma_{\beta_c,h} \lesssim (L^2 \beta_c)^{2D-1}$.
- **Applies to** $\mathbb{Z}^d$ Ising for $d \ge 5$ (mean-field bound is proven, refs [1,2,12]); to Curie–Weiss ($D=1$, optimal exponent); not yet to $d=2,3,4$ where mean-field bound is open.

## Method

Decompose $\mu_{\beta,h}$ via a Polchinski-type renormalisation flow $C_t = (tA + (\alpha - t)\mathbb{1})^{-1}$ into (i) a product (infinite-T) Bernoulli measure with known LSI and (ii) a continuous Gaussian-type measure $\nu_{0,\beta}$ to which the Polchinski-equation LSI criterion of Bauerschmidt–Bodineau is applied. The key new ingredient is bounding the spectral radius of the covariance $\Sigma_t(f)$ along the flow by $\chi_t$, achieved via FKG, the Ding–Song–Sun correlation inequality for Ising in external field, and Perron–Frobenius positivity — this avoids any Bakry–Émery $\beta < 1$ restriction and any spatial-mixing input.

## Quantum vs Classical

- **Baseline**: this *is* the classical baseline. Glauber on ferromagnetic Ising in the entire high-T phase mixes in $O(n \log n)$ via LSI; near criticality on $\mathbb{Z}^d$ ($d \ge 5$) one gets $t_{\mathrm{mix}} = \mathrm{poly}(n, 1/(\beta_c - \beta))$ uniformly in volume.
- **Gap**: none — there is no quantum advantage candidate in this regime. Ferromagnetic Ising is the *easy* commuting case; classical Glauber, FK / Swendsen–Wang (Guo–Jerrum), and Jerrum–Sinclair FPRAS already handle it polynomially at all $T$. A quantum KMS-DB Lindbladian would at best match this and in practice incurs a block-encoding / OFT / QSP overhead per step that erases any constant-factor win.
- **Source of the difference**: ferromagnetic structure (FKG, Griffiths inequalities, non-negative random-cluster representation) plus commuting Hamiltonian — no sign problem, no noncommutativity, no quantum coherence to exploit.
- **Caveat**: the result requires *ferromagnetic* couplings; antiferromagnetic / spin-glass Ising at intermediate-T is structurally different and is not covered. The mean-field bound on $\chi_\beta$ is open for $\mathbb{Z}^d$ in $d = 2,3,4$, so the polynomial-in-$(\beta_c-\beta)$ statement is conditional there.

## Implications for Quantum Advantage

- **Regime cell**: classical Ising × intermediate / near-critical (and high-T as a corollary).
- **What this changes**: closes the cell more tightly — together with Eldan–Koehler–Zeitouni (high-T spectral condition), Lubetzky–Sly (cutoff / critical 2D polynomial mixing), and Guo–Jerrum (random-cluster rapid mixing), the ferromagnetic Ising near-critical row is now classically tractable with a uniform-in-volume LSI, *not just a spectral gap*. LSI is the entropic functional inequality the quantum-side KMS-DB Lindbladian framework needs to match (modified-LSI / quantum-LSI), so the bar quantum samplers must clear in this cell is high and explicitly polynomial.
- **Promising or not**: not promising for quantum advantage. This is in the "ferromagnetic, classically easy" corner; the chapter punchline cell is non-stoquastic × low-T, far from here.

## Open Questions / Limitations

- Mean-field $\chi_\beta$ bound is open for $\mathbb{Z}^d$ Ising in $d \in \{2,3,4\}$ — exponents 1−2D in (1.8) are not the conjectured optimal ones.
- Method does not extend to the **low-T** phase (boundary-condition-dependent, Lifshitz law regime); explicitly flagged in the paper.
- Antiferromagnetic / non-FKG Ising not covered — FKG and the Ding–Song–Sun correlation inequality are essential.
- For Curie–Weiss the mean-field exponent is recovered with $D=1$ ($1/\gamma \sim (1-\beta/\beta_c)^{-1}$), matching known optimal results.

## Connections

- **Eldan–Koehler–Zeitouni (2022)** [`eldan-koehler-zeitouni-2022-spectral-condition-ising`]: complementary high-T result giving the Poincaré (spectral-gap) constant in terms of $\|J\|_{\mathrm{op}}$ via stochastic localization; Bauerschmidt–Dagallier promotes spectral gap to LSI and pushes uniformly to criticality (under mean-field hypothesis), at the cost of restricting to the ferromagnetic case.
- **Lubetzky–Sly (2012, 2013)** [`lubetzky-sly-2012-critical-ising-polynomial`, `lubetzky-sly-2013-cutoff-ising-lattice`]: 2D polynomial-time mixing at $\beta_c$ via RSW, and cutoff in the high-T phase. Bauerschmidt–Dagallier gives a *functional-inequality* analogue uniformly up to $\beta_c$ in $d \ge 5$ rather than RSW-driven 2D-only results.
- **Bauerschmidt–Bodineau (2019)**: high-T LSI via Polchinski / multiscale Bakry–Émery; this paper removes the $\beta < 1$ restriction by using the full Polchinski-equation criterion plus FKG.
- **Guo–Jerrum (2018)** [`guo-jerrum-2018-random-cluster-rapid-mixing`]: gives rapid mixing of FK / Swendsen–Wang for $q=2$ on any graph at any $T$ by polynomial-time gap; this paper gives the entropic (LSI) refinement on the Glauber side in the FKG regime.
- **Anari–Jain–Koehler–Pham–Vuong (entropic independence)**: parallel route to MLSI for high-T Ising via spectral / entropic independence; Bauerschmidt–Dagallier is geometric/PDE-flavoured (Polchinski) rather than high-dimensional-expander-flavoured.
- **Ding–Song–Sun (2022+)** [ref 15 in paper]: the new correlation inequality for Ising in external field — a key technical input here.
