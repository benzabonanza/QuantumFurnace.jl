# Quantum Signal Processing (v5 — pure Motlagh–Wiebe GQSP; no Low–Chuang framing; ordinary polynomials in $\mathbb{C}[z]$)

> **Insertion target:** `supplementary-informations/1_preliminaries.tex:1068`, body of `\subsection{Quantum Signal Processing.} \label{sec:prelim-QSP}` — between the `\subsection{...}\label{...}` line (line 1068) and `\begin{figure}[ht]` (line 1071). The figure (Fig. `circ:gqsp`, lines 1071–1113) sits immediately below the subsection heading; see Suggestions for the figure-redraw note, which is **required** for v5 to be self-consistent.
>
> **Change from v4:** two substantive corrections.
> 1. The walk operator is no longer $W = R_b\,U_B$ (a Low–Chuang-1019 qubitization reflection about the ancilla $\ket{0^b}$). It is the MW Cor. 8 form $W = -(I - 2\,\text{PREP}\ket{0}\bra{0}\text{PREP}^\dagger)\,\text{SELECT}$ — a reflection about $\text{PREP}\ket{0}$, sandwiched with PREP, composed with a bare SELECT.
> 2. The accessible polynomial class is ordinary $\mathbb{C}[z]$ per MW Cor. 5, **not** Laurent. The Jacobi–Anger Laurent expansion of $\ee^{-\ii\delta\alpha\cos\theta}$ is converted to an ordinary polynomial of degree $2d_{\text{JA}}$ on $\mathbb{T}$ by multiplication by $z^{d_{\text{JA}}}$ (a degree-preserving operation on $\mathbb{T}$). MW Thm. 6 — the $U^\dagger$-bearing "negative powers" route — is dropped entirely, in line with MW's own §III opening: *"our approach eschews the need to use $U^\dagger$ and avoids the use of Laurent polynomials in the analysis"*.
>
> **Existing labels referenced:** `sec:prelim-LCU`, `circ:gqsp`, `alg:coh`, `eq:B-block-encoding`, `eq:b_plus-s-eta`.
> **New labels introduced:** `eq:qsp-walk`, `eq:qsp-poly`, `eq:qsp-degree`. `eq:jacobi-anger` stays as in v4 (defensive label on the raw Jacobi–Anger identity).
> **New citation keys introduced:** `berntsonSunderhauf2025complementary` (Q-finding only). All other keys for this subsection are already in `references.bib`: `motlagh2024generalized`, `chen2023efficient`.
> **Citation keys dropped relative to v4:** `lowChuang2019qubitization`, `haah2019product`, `dong2021efficient`. See Writing Notes for the rationale; the Low–Chuang-framed qubitization walk is not used in v5, and MW's own angle-extraction beats Haah/Dong on every benchmark that matters for us.

---

## Body

The LCU circuit of (§`sec:prelim-LCU`), specialised to the nested kernel construction of Chapter 5, delivers a block encoding $U_B$ on a $b$-qubit block register of the Hermitian coherent term $B$,
$$
\bra{0^b}\,U_B\,\ket{0^b} \;=\; B/\alpha, \qquad \|B/\alpha\| \leq 1, \qquad \alpha \;=\; \|b_-\|_1\,\|b_+\|_1.
$$
Unpacking the LCU, $B$ is realised as an explicit linear combination of unitaries,
$$
B \;=\; \sum_j \alpha_j\,U_j, \qquad \alpha \;=\; \sum_j \alpha_j, \qquad \alpha_j \geq 0,
$$
with a preparation unitary $\text{PREP}\ket{0} = \sum_j \sqrt{\alpha_j/\alpha}\,\ket{j}$ on the block ancilla and a selection unitary $\text{SELECT} = \sum_j \ket{j}\!\bra{j}\otimes U_j$ on the combined block-plus-system register, related to $U_B$ by $U_B = (\text{PREP}^\dagger\otimes I)\,\text{SELECT}\,(\text{PREP}\otimes I)$. What the coherent step of the algorithm actually needs from this block encoding is the unitary $\ee^{-\ii\delta B}$. *Generalized quantum signal processing* (GQSP) [Motlagh & Wiebe 2024] [CITE: motlagh2024generalized] is the tool that turns one into the other: interleaving controlled calls to a walk operator with a short sequence of single-qubit rotations realises any complex polynomial transformation of the walk's block-encoded spectrum, and a judicious polynomial approximates $\ee^{-\ii\delta B}$. The circuit in (Fig. `circ:gqsp`) is the GQSP layer for arbitrary degree $d$.

*Walk operator.* Following [Motlagh & Wiebe 2024, Cor. 8] [CITE: motlagh2024generalized], the unitary that GQSP acts on is
$$
W \;:=\; -\bigl(I \,-\, 2\,\text{PREP}\ket{0}\!\bra{0}\text{PREP}^\dagger\bigr)\,\text{SELECT}
\;=\; \bigl(2\,\text{PREP}\ket{0}\!\bra{0}\text{PREP}^\dagger \,-\, I\bigr)\,\text{SELECT},
$$ <!-- \label{eq:qsp-walk} -->
i.e. a bare SELECT followed by a reflection about the prepared state $\text{PREP}\ket{0}$ (the reflection *sandwiched* with PREP/PREP$^\dagger$; a full $U_B$ is not applied, only one SELECT and one PREP-sandwich). For each eigenvector $\ket{\lambda}$ of $B$ with $B\ket{\lambda} = \lambda\ket{\lambda}$, the two-dimensional subspace spanned by $\text{PREP}\ket{0}\ket{\lambda}$ and an orthogonal state $\ket{\phi_\lambda^-}$ (whose definition is fixed by the action of $W$ on $\text{PREP}\ket{0}\ket{\lambda}$, [Motlagh & Wiebe 2024, Eq. (64)] [CITE: motlagh2024generalized]) is $W$-invariant. In the basis $(\text{PREP}\ket{0}\ket{\lambda},\ket{\phi_\lambda^-})$, $W$ acts as the $\mathrm{SU}(2)$ matrix
$$
W_\lambda \;=\; \begin{pmatrix} \lambda/\alpha & \sqrt{1-(\lambda/\alpha)^2} \\[2pt] -\sqrt{1-(\lambda/\alpha)^2} & \lambda/\alpha \end{pmatrix},
$$
of trace $2\lambda/\alpha$ and determinant $+1$ [Motlagh & Wiebe 2024, Eq. (64)] [CITE: motlagh2024generalized]. Its eigenvalues form a conjugate pair $\ee^{\pm\ii\theta}$ with
$$
\cos\theta \;=\; \lambda/\alpha, \qquad \lambda \in [-\alpha,\alpha], \qquad \theta \in [0,\pi].
$$
Implementing $\ee^{-\ii\delta B}$ on the $\lambda$-eigenspace of $B$ therefore reduces to finding a *polynomial* $P$ that approximates $\ee^{-\ii\delta\alpha\cos\theta}$ on the unit circle $\mathbb{T} = \{z\in\mathbb{C} : |z|=1\}$, evaluated at $z = \ee^{\ii\theta}$.

*GQSP layer and polynomial class.* The circuit in (Fig. `circ:gqsp`) alternates $d$ (open-)controlled calls to $W$ with $d+1$ single-qubit rotations on the QSP ancilla $\ket{0}_\mathrm{QSP}$: a three-parameter opening rotation $R(\theta_0,\phi_0,\lambda_0)$ and $d$ further two-parameter rotations $R(\theta_k,\phi_k,0)$, $k = 1,\dots,d$. (The third parameter $\lambda_0 \in \mathbb{R}$ of the opening rotation, distinct from the spectral variable $\lambda$, absorbs a single global phase per [Motlagh & Wiebe 2024, Thm. 3] [CITE: motlagh2024generalized].) Crucially — and unlike QSVT-style constructions — GQSP requires no controlled $W^\dagger$ gates. Post-selecting both ancillas on $\ket{0}_\mathrm{QSP}\otimes(\text{PREP}\ket{0})$ returns, on each qubitization subspace, the top-left block of the composed circuit, which is an ordinary polynomial $P(\ee^{\ii\theta}) \in \mathbb{C}[z]$ evaluated at the eigenvalue $z = \ee^{\ii\theta}$ of $W$. [Motlagh & Wiebe 2024, Thm. 3 and Cor. 5] [CITE: motlagh2024generalized] characterise the accessible polynomial class: for any $P \in \mathbb{C}[z]$ of degree $\deg P \leq d$ satisfying
$$
\|P\|_{\infty,\mathbb{T}} \;:=\; \max_{z\in\mathbb{T}}\, |P(z)| \;\leq\; 1,
$$ <!-- \label{eq:qsp-poly} -->
there exist rotation angles $(\theta_0,\phi_0,\lambda_0)$ and $(\theta_k,\phi_k)_{k=1}^d$ realising $P$ in that top-left block. Note that the accessible class is **ordinary** polynomials, not Laurent — no negative powers of $\ee^{\ii\theta}$ appear, and the construction uses $W$ only, never $W^\dagger$.

*Target polynomial via Jacobi–Anger.* The target $\ee^{-\ii\delta\alpha\cos\theta}$ is naturally a *Laurent* series: the Jacobi–Anger identity [Motlagh & Wiebe 2024, Eq. (62)] [CITE: motlagh2024generalized] reads
$$
\ee^{-\ii\delta\alpha\cos\theta} \;=\; \sum_{k=-\infty}^{\infty} (-\ii)^k\,J_k(\delta\alpha)\,\ee^{\ii k\theta}, \qquad \theta \in \mathbb{R},
$$ <!-- \label{eq:jacobi-anger} -->
with Bessel coefficients; the sign pattern $(-\ii)^k$ comes from setting $t \to -\delta\alpha$ in MW's identity $\ee^{\ii t\cos\theta} = \sum_k \ii^k J_k(t)\,\ee^{\ii k\theta}$ (and using $J_{-k}(z) = (-1)^k J_k(z)$). Truncating at $|k| \leq d_{\text{JA}}$ yields a Laurent polynomial
$$
L_{d_{\text{JA}}}(\ee^{\ii\theta}) \;:=\; \sum_{k=-d_{\text{JA}}}^{d_{\text{JA}}} (-\ii)^k\,J_k(\delta\alpha)\,\ee^{\ii k\theta},
$$
with uniform tail bound $\|\ee^{-\ii\delta\alpha\cos\theta} - L_{d_{\text{JA}}}\|_{\infty,\mathbb{T}} \leq 2\sum_{|k|>d_{\text{JA}}}|J_k(\delta\alpha)| \leq \mathcal{O}\bigl((\ee\delta\alpha/(2d_{\text{JA}}))^{d_{\text{JA}}}\bigr)$ from the rapid Bessel decay [Motlagh & Wiebe 2024, Eq. (63)] [CITE: motlagh2024generalized]. *This Laurent polynomial is not directly in GQSP's polynomial class.* To land in $\mathbb{C}[z]$, we multiply by the unimodular factor $\ee^{\ii d_{\text{JA}}\theta} = z^{d_{\text{JA}}}$,
$$
P_{2d_{\text{JA}}}(z) \;:=\; z^{d_{\text{JA}}}\, L_{d_{\text{JA}}}(z) \;=\; \sum_{m=0}^{2d_{\text{JA}}} c_m\,z^m, \qquad c_m \;=\; (-\ii)^{m-d_{\text{JA}}}\,J_{m-d_{\text{JA}}}(\delta\alpha).
$$
The multiplication by $z^{d_{\text{JA}}}$ is unit-modulus on $\mathbb{T}$, hence $\|P_{2d_{\text{JA}}}\|_{\infty,\mathbb{T}} = \|L_{d_{\text{JA}}}\|_{\infty,\mathbb{T}} \leq 1 + \mathcal{O}((\ee\delta\alpha/(2d_{\text{JA}}))^{d_{\text{JA}}})$, and a rescaling by $(1-\varepsilon_{\text{QSP}})$ enforces (Eq. `eq:qsp-poly`) strictly. The rescaled $P_{2d_{\text{JA}}}$ is an *ordinary* polynomial of degree $2d_{\text{JA}}$, realisable by GQSP via (Thm. 3, Cor. 5).

*Degree bound and cost.* The Jacobi–Anger truncation index that suffices for $\varepsilon_{\text{QSP}}$-accuracy on $\mathbb{T}$ is, per [Motlagh & Wiebe 2024, Thm. 7] [CITE: motlagh2024generalized],
$$
d_{\text{JA}} \;=\; \mathcal{O}\!\left(\delta\alpha \;+\; \frac{\log(1/\varepsilon_\mathrm{QSP})}{\log\log(1/\varepsilon_\mathrm{QSP})}\right),
$$ <!-- \label{eq:qsp-degree} -->
and the ordinary polynomial $P_{2d_{\text{JA}}}$ that GQSP implements therefore has degree $2d_{\text{JA}} = \mathcal{O}\bigl(\delta\alpha + \log(1/\varepsilon_\mathrm{QSP})/\log\log(1/\varepsilon_\mathrm{QSP})\bigr)$; the factor of $2$ from the Laurent-to-ordinary shift is absorbed into the big-$\mathcal{O}$. The GQSP circuit in (Fig. `circ:gqsp`) makes $2d_{\text{JA}}$ controlled calls to $W$ (one per polynomial degree) interleaved with $2d_{\text{JA}}+1$ single-qubit rotations on the QSP ancilla, matching MW's Cor. 8 bound.

*Small-$\delta\alpha$ default.* In the regime of Chapter 5, the coherent step size $\delta$ is small and $\delta\alpha \lesssim 1$ is the norm. Already $d_{\text{JA}} = 1$ gives, after the $z$-shift,
$$
P_2(z) \;=\; -\ii\,J_1(\delta\alpha)\,z^2 \,+\, J_0(\delta\alpha)\,z \,-\, \ii\,J_1(\delta\alpha),
$$
an ordinary polynomial of degree $2$ implemented with two controlled-$W$ calls and three single-qubit rotations. On each qubitization subspace $P_2(\ee^{\ii\theta}) = \ee^{\ii\theta}(J_0 - 2\ii J_1\cos\theta) \approx \ee^{\ii\theta}\,(1 - \ii\delta\alpha\cos\theta) \approx \ee^{\ii\theta}\,\ee^{-\ii\delta\alpha\cos\theta}$ using the small-argument expansion $J_0(z) = 1 + \mathcal{O}(z^2)$, $2J_1(z) = z + \mathcal{O}(z^3)$, reproducing $\ee^{-\ii\delta\alpha\cos\theta}$ to $\mathcal{O}((\delta\alpha)^2)$ up to the trailing unimodular $\ee^{\ii\theta}$ factor that is absorbed into the qubitization bookkeeping (see Writing Notes on the shift-by-$z^{d_{\text{JA}}}$ identity). Two regimes force $d_{\text{JA}} \geq 2$: (a) when $\alpha$ is large — e.g. the Metropolis kernel $b_+^{(s,\eta)}$ (Eq. `eq:b_plus-s-eta`) with $\ell_1$-norm scaling as $\log(\beta\|H\|\,\|\!\sum_a A_a^\dagger A_a\|/\varepsilon)$ (see `2_methods.tex:1726`), so that $\delta\alpha$ is no longer small; (b) when $\delta$ is coarsened to reduce the total number of coherent steps at fixed overall accuracy. In both cases $d_{\text{JA}}$ grows linearly in $\delta\alpha$ per (Eq. `eq:qsp-degree`) and the polynomial degree $2d_{\text{JA}}$ grows with it.

*Angle finding.* Given the coefficients $\{c_m\}_{m=0}^{2d_{\text{JA}}}$ of $P_{2d_{\text{JA}}}$, the GQSP angle vector $\{(\theta_0,\phi_0,\lambda_0)\}\cup\{(\theta_k,\phi_k)\}_{k=1}^{2d_{\text{JA}}}$ is extracted in two steps. First, a *complementary polynomial* $Q \in \mathbb{C}[z]$ of degree $\deg Q = \deg P_{2d_{\text{JA}}} = 2d_{\text{JA}}$ satisfying
$$
|P_{2d_{\text{JA}}}(z)|^2 + |Q(z)|^2 = 1 \qquad\text{on } \mathbb{T}
$$
is constructed; its existence is [Motlagh & Wiebe 2024, Thm. 4] [CITE: motlagh2024generalized]. Second, the exact recursive extraction of [Motlagh & Wiebe 2024, §IV, Alg. 1] [CITE: motlagh2024generalized] reads off the $(2d_{\text{JA}}+1)$ rotation triples from $(P_{2d_{\text{JA}}}, Q)$ in $\mathcal{O}(d_{\text{JA}}^2)$ arithmetic operations. For the $Q$-finding step, MW propose an FFT-based convolution objective (Eqs. (58)–(60)) solved by nonlinear optimisation. We replace that step with the contour-integral / FFT construction of [Berntson & Sünderhauf 2025] [CITE: berntsonSunderhauf2025complementary], which evaluates $Q$ in the monomial basis in classical time $\tilde{\mathcal{O}}(d_{\text{JA}})$ with explicit a priori error bounds and has been numerically demonstrated up to $\deg P \sim 10^7$. This is a practical gain over MW's own optimisation route — no solver initialisation, no convergence tuning — and delivers the same coefficient vector $\{b_m\}$. The MW angle-extraction recursion (Alg. 1) then proceeds unchanged.

The unitary $\ee^{-\ii\delta B}$ produced by (Fig. `circ:gqsp`) closes the block-encoding pipeline begun in (§`sec:prelim-LCU`): it is the coherent step of (Alg. `alg:coh`). Its cost is $2d_{\text{JA}}$ applications of $W$ (each costing two PREP calls and one SELECT call, plus the PREP-sandwiched reflection) together with $\mathcal{O}(d_{\text{JA}})$ single-qubit gates, the $W$-calls dominating; with $d_{\text{JA}}$ per (Eq. `eq:qsp-degree`), the coherent step contributes a $\tilde{\mathcal{O}}(\delta\alpha)$ factor per $\delta$-step of the Gibbs sampler.

---

## Citations

- **`berntsonSunderhauf2025complementary`** — Berntson, B. K. & Sünderhauf, C., *"Complementary Polynomials in Quantum Signal Processing"*, **Commun. Math. Phys. 406, 161 (2025)**, arXiv:2406.04246. Contour-integral / FFT construction of $Q_d$ in classical time $\tilde{\mathcal{O}}(d)$ with explicit error bounds; numerically demonstrated to $\deg P \sim 10^7$. Replaces MW's §IV FFT-optimisation step for complex $P$.
- **`motlagh2024generalized`** — already in `references.bib`. Thm. 3 and Cor. 5 characterise the achievable GQSP polynomial class ($\mathbb{C}[z]$, $\|P\|_{\infty,\mathbb{T}}\leq 1$). Thm. 7 gives the Jacobi–Anger degree bound. Cor. 8 gives the PREP-sandwiched walk $W = -(I - 2\,\text{PREP}\ket{0}\bra{0}\text{PREP}^\dagger)\text{SELECT}$ and its $2{\times}2$ block action $W_\lambda$ (Eq. 64). Algorithm 1 in §IV gives the $\mathcal{O}(d^2)$ recursive angle extraction.

```bibtex
@article{berntsonSunderhauf2025complementary,
  author        = {Berntson, Bjorn K. and S\"underhauf, Christoph},
  title         = {Complementary Polynomials in Quantum Signal Processing},
  journal       = {Commun. Math. Phys.},
  volume        = {406},
  number        = {7},
  pages         = {161},
  year          = {2025},
  doi           = {10.1007/s00220-025-05302-9},
  eprint        = {2406.04246},
  archivePrefix = {arXiv},
  primaryClass  = {quant-ph},
}
```

---

## Writing Notes

### Why v5 drops the Low–Chuang qubitization framing

v4 cited `lowChuang2019qubitization` in two places: (i) to name the walk $W = R_b\,U_B$ as a "qubitization walk" in the sense of [Low & Chuang 2019, Lemma 10], and (ii) to attribute the additive-optimal $\mathcal{O}(t + \log(1/\varepsilon)/\log\log(1/\varepsilon))$ scaling. Both uses are problematic once we take MW Cor. 8 seriously:

1. MW's walk is $W = -(I - 2\,\text{PREP}\ket{0}\bra{0}\text{PREP}^\dagger)\text{SELECT}$: the reflection is about $\text{PREP}\ket{0}$ (a specific superposition in the block register), *not* about $\ket{0^b}$ (the bare ancilla zero), and the operator in front of it is a bare SELECT, not a full $U_B = \text{PREP}^\dagger\,\text{SELECT}\,\text{PREP}$. Writing "$W = R_b\,U_B$" with $R_b = 2\ket{0^b}\bra{0^b} - I$ is a different walk. Both are qubitization walks in an abstract sense, but MW's Cor. 8 (which we invoke to get $W_\lambda \mapsto F_\lambda$) is stated for the PREP-sandwiched form, not the Low–Chuang form.
2. The scaling can be attributed to MW Thm. 7 directly; citing Low–Chuang for the same scaling is historical padding that does not affect the argument.

So v5 uses MW throughout — one reference covers the walk (Cor. 8), the polynomial class (Thm. 3, Cor. 5), the degree bound (Thm. 7), and the angle extraction (Alg. 1). The figure `circ:gqsp` has to be updated to match; see Suggestions.

### Why Thm. 6 is wrong for v5

MW Thm. 6 ("Polynomials With Negative Powers") uses an auxiliary signal operator $A' = \ket{0}\!\bra{0}\otimes I + \ket{1}\!\bra{1}\otimes U^\dagger$ to extend GQSP to Laurent polynomials in $z$. This is an alternative route to Hamiltonian simulation when one wants to implement the Jacobi–Anger Laurent expansion directly, at the cost of needing controlled $U^\dagger$ (i.e. controlled $W^\dagger$). MW's §III opening sentence explicitly advertises that their construction avoids both $U^\dagger$ and Laurent polynomials; v5 honours this by using ordinary polynomials and the shift trick below.

### The shift-by-$z^{d_{\text{JA}}}$ step and the Cor. 8 identification [CHECK]

The identity invoked in the body — "$P_{2d_{\text{JA}}}(W) = W^{d_{\text{JA}}}\,L_{d_{\text{JA}}}(W)$ acts on the qubitization subspace as $\ee^{-\ii\delta B}$ up to a correctable phase" — deserves a careful statement. For each eigenvector $\ket{\lambda}$ of $B$, in the $W$-invariant basis $\ket{\phi_\lambda^\pm}$ with $W\ket{\phi_\lambda^\pm} = \ee^{\pm\ii\theta_\lambda}\ket{\phi_\lambda^\pm}$,
$$
P_{2d_{\text{JA}}}(W)\ket{\phi_\lambda^\pm} \;=\; \ee^{\pm\ii d_{\text{JA}}\theta_\lambda}\,L_{d_{\text{JA}}}(\ee^{\pm\ii\theta_\lambda})\,\ket{\phi_\lambda^\pm} \;\approx\; \ee^{\pm\ii d_{\text{JA}}\theta_\lambda}\,\ee^{-\ii\lambda\delta}\,\ket{\phi_\lambda^\pm}.
$$
MW Cor. 8 claims the action of GQSP on the 2D subspace is $\operatorname{diag}(\ee^{-\ii\lambda\delta},\ee^{-\ii\lambda\delta})$ (their Eq. 66) — i.e. a $\lambda$-dependent global phase with no $\ee^{\pm\ii d_{\text{JA}}\theta_\lambda}$ mismatch. The reconciliation is that MW implicitly use Cor. 8's statement "apply Cor. 5 for any polynomial $P$ such that $\max\|P(U)\|\leq 1$" to **define** the polynomial $P$ via the desired target $F_\lambda$ — and the polynomial that achieves $P(e^{\ii\theta}) = \ee^{-\ii\delta\alpha\cos\theta}$ and $P(\ee^{-\ii\theta}) = \ee^{-\ii\delta\alpha\cos\theta}$ simultaneously is Laurent, not ordinary. So MW's Cor. 8 actually *uses* Laurent polynomials under the hood, in tension with their §III pitch.

For v5, the operationally honest statement is: the ordinary polynomial $P_{2d_{\text{JA}}}(W)$ implements $\ee^{-\ii\delta B}\cdot W^{d_{\text{JA}}}|_{\text{qub. subspace}}$ rather than $\ee^{-\ii\delta B}$ exactly. The trailing $W^{d_{\text{JA}}}$ is a known, implementable unitary (it acts on the 2D subspace as $W_\lambda^{d_{\text{JA}}}$), so in principle it can be pre-compensated by conjugating the input state — or one re-uses the MW Thm. 6 (Laurent) route. The v5 body states the polynomial target via Cor. 5 cleanly, then the "up to the trailing unimodular $\ee^{\ii\theta}$ factor" is flagged in the $P_2$ paragraph and elaborated here. **[CHECK]** — the most honest version of the Hamiltonian-simulation construction within pure MW Cor. 5 is to apply two GQSP circuits, one for $\cos$ and one for $\sin$, and combine via Berry et al.'s doubling efficiency (LCU with coefficients $\frac{1}{2}(1,\ii)$). The v5 body presents the single-circuit shift route because that is what MW Cor. 8 literally states, but the reader should be aware that this route smuggles in Laurent machinery at the Cor. 8 step. If the thesis needs a fully rigorous single-GQSP-call route without Laurent, the user should reinstate Thm. 6 in v6 and dial back §III pitch accordingly; or switch to the two-call ($\cos + \ii\sin$) LCU route.

### The factor of 2 in the polynomial degree

The "degree" quoted in MW Thm. 7 and in (Eq. `eq:qsp-degree`) is the Jacobi–Anger truncation index $d_{\text{JA}}$. The ordinary polynomial $P_{2d_{\text{JA}}}$ implemented by GQSP has degree $2d_{\text{JA}}$ (because the shift $z^{d_{\text{JA}}}$ converts Laurent range $[-d_{\text{JA}}, d_{\text{JA}}]$ to ordinary range $[0, 2d_{\text{JA}}]$). The number of controlled-$W$ calls in (Fig. `circ:gqsp`) is $2d_{\text{JA}}$, one per polynomial degree, *not* $d_{\text{JA}}$. The factor of 2 matters for the explicit cost accounting in Chapter 5; it is absorbed into the big-$\mathcal{O}$ in the degree-bound statement but should be tracked when the algorithm is actually scheduled. MW Thm. 7's statement "$\mathcal{O}(\alpha t + \log(1/\varepsilon)/\log\log(1/\varepsilon))$ controlled-$U$ operations" is cost-wise correct after the factor-of-2 absorption; the body of v5 calls this out explicitly.

### Derivation of the $P_2$ expansion under $\cos\theta$ and the shift

Starting from (Eq. `eq:jacobi-anger`) truncated at $|k| \leq 1$: $L_1(\ee^{\ii\theta}) = (-\ii)^{-1} J_{-1}\,\ee^{-\ii\theta} + J_0 + (-\ii)^1 J_1\,\ee^{\ii\theta} = \ii\,(-J_1)\,\ee^{-\ii\theta} + J_0 - \ii J_1\,\ee^{\ii\theta} = -\ii J_1(\ee^{-\ii\theta}+\ee^{\ii\theta}) + J_0 = J_0 - 2\ii J_1\cos\theta$ (using $J_{-1}(z) = -J_1(z)$ and $(-\ii)^{-1} = \ii$). Multiplying by $z = \ee^{\ii\theta}$,
$$
P_2(\ee^{\ii\theta}) \;=\; \ee^{\ii\theta}(J_0 - 2\ii J_1\cos\theta) \;=\; J_0\,\ee^{\ii\theta} - \ii J_1\,(\ee^{2\ii\theta} + 1) \;=\; -\ii J_1\,\ee^{2\ii\theta} + J_0\,\ee^{\ii\theta} - \ii J_1,
$$
reading off $c_0 = -\ii J_1$, $c_1 = J_0$, $c_2 = -\ii J_1$ — the degree-2 ordinary polynomial in the body. Leading-order small-$\delta\alpha$ expansion: $J_0(z) = 1 + \mathcal{O}(z^2)$, $2J_1(z) = z + \mathcal{O}(z^3)$, so $P_2(\ee^{\ii\theta}) = \ee^{\ii\theta}(1 - \ii\delta\alpha\cos\theta + \mathcal{O}((\delta\alpha)^2)) = \ee^{\ii\theta}\,\ee^{-\ii\delta\alpha\cos\theta} + \mathcal{O}((\delta\alpha)^2)$. All signs cross-checked against DLMF 10.2.2, 10.12.1, and MW Eq. (62).

### Scope and length

- Six body paragraphs + three numbered displays (`eq:qsp-walk`, `eq:qsp-poly`, `eq:qsp-degree`) + defensively-labelled `eq:jacobi-anger` + one inline $W_\lambda$ matrix + one inline $P_{2d_{\text{JA}}}$ display + one inline $P_2$ display. ~1.3–1.5 pages of typeset body, matching v4.
- Equation labels: `eq:qsp-walk`, `eq:qsp-poly`, `eq:qsp-degree` are all candidates for Chapter 5 cost-analysis cross-references.

### Notation choices

- $\lambda \in [-\alpha,\alpha]$ denotes the eigenvalue of $B$ in physical units (not $B/\alpha$); MW Eq. 64 has $\lambda/\alpha$ everywhere and Eq. 66 has $\ee^{-\ii\lambda t}$, matching.
- $d_{\text{JA}}$ is used for the Jacobi–Anger truncation index to avoid the overload of $d$ (which in Fig. `circ:gqsp` and MW's figures means the GQSP circuit depth = number of controlled-$W$ layers). In v5 the circuit depth is $2d_{\text{JA}}$ and should be documented as such in Chapter 5.
- $\mathbb{T} = \{|z|=1\}$ for the complex unit circle, matching MW and Berntson–Sünderhauf.
- $\text{PREP}\ket{0}$ is the specific superposition $\sum_j \sqrt{\alpha_j/\alpha}\,\ket{j}$ from the LCU section; the "reflection about $\text{PREP}\ket{0}$" $2\,\text{PREP}\ket{0}\bra{0}\text{PREP}^\dagger - I$ is distinct from the raw-ancilla reflection $R_b = 2\ket{0^b}\bra{0^b} - I$ and should not be conflated.
- $\ee$, $\ii$ respected throughout; ancilla labels $\ket{0}_\mathrm{QSP}$, $\ket{0^b}$ match LCU v2 and (Fig. `circ:gqsp`) after its redraw.

### Angle-finding complexity

- The complementary-polynomial step is $\tilde{\mathcal{O}}(d_{\text{JA}})$ via Berntson–Sünderhauf (no optimisation solver required); the MW angle-extraction recursion (§IV, Alg. 1) is $\mathcal{O}(d_{\text{JA}}^2)$ arithmetic operations in exact arithmetic. The body states both explicitly. A single historical-context sentence (MW's §IV comment that traditional QSP angle finding topped out at $\sim 10^4$ degrees, vs BS demonstrating $\sim 10^7$) is worth adding if Chapter 5 later needs a low-temperature benchmark; v5 omits it for length. Dropping `haah2019product` and `dong2021efficient` relative to v4 reflects this.

---

## Suggestions

Items that would polish the subsection but require either a user decision or lie slightly outside v5's scope:

1. **[FIX-NEXT] Redraw `circ:gqsp` (required for v5 self-consistency).** The figure at `1_preliminaries.tex:1071–1113` currently draws $W = R_b\,U_B$ as the sequence "bundled $U_B$ block" followed by a "box labelled $R_b$" — an external reflection about $\ket{0^b}$. Under MW Cor. 8, v5's walk is $W = -(I - 2\,\text{PREP}\ket{0}\bra{0}\text{PREP}^\dagger)\,\text{SELECT}$ — a bare SELECT followed by a PREP-sandwiched reflection about $\text{PREP}\ket{0}$. Two ways to fix:
   - *Option A (honest; what v5's body assumes):* Expand the $U_B$ block into its PREP / SELECT / PREP$^\dagger$ composition, drop the outer PREP$^\dagger$ and PREP that would make a full $U_B$, and instead draw, per GQSP layer: (i) one open-controlled SELECT; (ii) the PREP-sandwiched reflection $2\,\text{PREP}\ket{0}\bra{0}\text{PREP}^\dagger - I$ drawn as a box `PREP$^\dagger$, reflection about $\ket{0^b}$, PREP`; (iii) the $R(\theta_k,\phi_k,0)$ rotation on the QSP ancilla. This exposes the Cor. 8 structure explicitly and removes the mismatch with the body. **This is the recommended redraw.**
   - *Option B (bundled; acceptable but less informative):* Keep the $U_B$ block and rename the external $R_b$ box to reflect the MW Cor. 8 reality, e.g. label it "PREP-sandwiched reflection about $\text{PREP}\ket{0}$", and state in the caption that $U_B$ internally contains the extra PREP$^\dagger$/PREP pair cancellation. The picture still hides the construction but the caption makes the walk unambiguous.
   Either way, the caption text "*degree-$d$ Chebyshev polynomial $P(W)\approx\ee^{-\ii\delta B}$*" (lines 1110–1111) **must be rewritten** in both options — GQSP produces **ordinary polynomials in $\mathbb{C}[z]$ of degree $2d_{\text{JA}}$**, not Chebyshev polynomials. Suggested caption: "Degree-$2d_{\text{JA}}$ GQSP realising $\ee^{-\ii\delta B}$ from the walk operator $W = -(I - 2\,\text{PREP}\ket{0}\bra{0}\text{PREP}^\dagger)\,\text{SELECT}$ (MW Cor. 8), with $\text{PREP}\ket{0}$ the prepared state of the LCU $B = \sum_j \alpha_j U_j$ and SELECT the controlled-$U_j$ oracle; post-selecting the QSP ancilla on $\ket{0}_\mathrm{QSP}$ and the block register on $\text{PREP}\ket{0}$ yields the ordinary polynomial $P_{2d_{\text{JA}}}(W)$ realising $\ee^{-\ii\delta B}$ on the qubitization subspace via the Jacobi–Anger truncation + $z^{d_{\text{JA}}}$ shift".

2. **[FIX-NEXT] Reconcile `2_methods.tex` with v5.** The methods chapter must speak the same convention as v5. Under v5, `2_methods.tex` needs:
   - **Line 1609** (algorithm pseudocode, inside `\textsc{CoherentStep}`): "$P(\ee^{\ii\theta}) \approx \ee^{-\ii\delta\alpha\sin\theta}$" $\to$ "$P_{2d_{\text{JA}}}(\ee^{\ii\theta}) \approx \ee^{\ii\theta d_{\text{JA}}}\ee^{-\ii\delta\alpha\cos\theta}$ (ordinary polynomial in $\mathbb{C}[z]$, degree $2d_{\text{JA}}$)", and update the `\Comment{Degree $d=1$ calls to $U_{B_a}$}` annotation to `\Comment{$d_{\text{JA}}=1$ truncation → $2d_{\text{JA}}=2$ calls to $W$}`.
   - **Line 1731**: "its eigenvalues come in conjugate pairs $\ee^{\pm\ii\theta}$ where $\sin\theta = \lambda$ for each eigenvalue $\lambda$ of $B_a/\alpha$" $\to$ "its eigenvalues come in conjugate pairs $\ee^{\pm\ii\theta}$ where $\cos\theta = \lambda/\alpha$ for each eigenvalue $\lambda$ of $B_a$". Same line: "$P(\ee^{\ii\theta}) \approx \ee^{-\ii\delta\alpha\sin\theta}$" $\to$ "$P_{2d_{\text{JA}}}(\ee^{\ii\theta}) \approx \ee^{\ii\theta d_{\text{JA}}}\ee^{-\ii\delta\alpha\cos\theta}$ on $\mathbb{T}$ (ordinary polynomial)". Also: bring in the walk operator $W$ explicitly, with text pointing to (Eq. `eq:qsp-walk`) in the preliminaries subsection.
   - **Lines 1733–1735** (degree bound): the statement $d = \mathcal{O}(\delta\alpha + \log(1/\varepsilon_{\text{QSP}})/\log\log(1/\varepsilon_{\text{QSP}}))$ is correct if "$d$" means $d_{\text{JA}}$. Rename to $d_{\text{JA}}$ and add a parenthetical "(the GQSP polynomial degree, and hence the number of controlled-$W$ calls per coherent step, is $2d_{\text{JA}}$)".
   - **Line 1737** (default $d = 1$): "degree $d = 1$ is sufficient for the QSP" $\to$ "Jacobi–Anger truncation $d_{\text{JA}} = 1$ is sufficient, giving a degree-2 ordinary polynomial and two controlled-$W$ calls per coherent step".
   - **Cost accounting:** the $R_b$ and "$U_B$" references in the methods cost table should become "$W$-calls", with one $W$ = one SELECT + one PREP-sandwiched reflection + $\mathcal{O}(1)$ single-qubit gates; at $2d_{\text{JA}}$ calls per coherent step.
   - The user memory `.claude-memory/todo_methods_gqsp_reflection.md` tracks this reconciliation; once the methods edit lands, that memo can be closed. The `.claude-memory/thesis_sign_convention_trotter.md` note is unrelated and stays open.

3. **Explicit $d_{\text{JA}} = 1$ GQSP angle formulas.** For the default, the five rotation parameters $(\theta_0, \phi_0, \lambda_0, \theta_1, \phi_1, \theta_2, \phi_2)$ of a degree-2 GQSP circuit reduce to closed-form expressions in $J_0(\delta\alpha)$ and $J_1(\delta\alpha)$ via the MW Alg. 1 recursion applied to $P_2 = -\ii J_1\,z^2 + J_0\,z - \ii J_1$ and its complementary $Q_2$. Listing these explicitly would make the "default in practice" statement self-verifying. Requires separate verification against the numerical implementation; deferred. **Note the $2d_{\text{JA}}+1 = 3$ rotation pairs, not $d_{\text{JA}}+1 = 2$** — another place the factor of 2 is operationally relevant.

4. **Label hooks in Chapter 5.** If the user adds `\label{sec:B-nested-LCU}` at `2_methods.tex:1707` and `\label{sec:B-GQSP}` at `2_methods.tex:1730`, the closing sentence of v5 could point more precisely than `(Alg. alg:coh)`. One-line user edit; not a blocker.

5. **"GQSP" vs "QSP" subsection title.** The existing label `sec:prelim-QSP` and the `\subsection{Quantum Signal Processing.}` heading at `1_preliminaries.tex:1068` use "QSP" as an umbrella. Everything in v5 is strictly GQSP. Consider retitling to `\subsection{Generalized Quantum Signal Processing.}` with label `sec:prelim-GQSP`, and update back-references in `2_methods.tex`. Low priority; current choice is defensible.

6. **[Resolve [CHECK]] The Cor. 8 / shift identity.** See Writing Notes "The shift-by-$z^{d_{\text{JA}}}$ step and the Cor. 8 identification [CHECK]". The cleanest resolution is for v6 to either (a) use Berry et al.'s two-call LCU (`cos + i sin`) approach with explicit factor-of-2 overhead, or (b) reinstate MW Thm. 6's Laurent route and soften §III's "no Laurent" pitch accordingly. v5 presents MW Cor. 8 as written and flags the subtlety, but does not resolve it. Draft-checker should verify that the shift-trick operator identity $P_{2d_{\text{JA}}}(W) = W^{d_{\text{JA}}}\,L_{d_{\text{JA}}}(W)$ on the qubitization subspace matches MW Eq. 66 up to the noted $\ee^{\pm\ii d_{\text{JA}}\theta_\lambda}$ factor — and whether MW's Cor. 8 silently uses Laurent under the hood.

---

## Diff from v4

### Paragraphs rewritten

- **Walk-operator paragraph.** v4: $W = R_b\,U_B$ with $R_b = 2\ket{0^b}\bra{0^b} - I$, cited [Low & Chuang 2019, Lemma 10]. v5: $W = -(I - 2\,\text{PREP}\ket{0}\bra{0}\text{PREP}^\dagger)\,\text{SELECT}$, cited [Motlagh & Wiebe 2024, Cor. 8]. The 2D invariant subspace description is rewritten to use $\text{PREP}\ket{0}\ket{\lambda}$ and $\ket{\phi_\lambda^-}$ as basis vectors (MW's convention), not the abstract Low–Chuang basis. The $W_\lambda$ $2{\times}2$ matrix, trace/det argument, and $\cos\theta = \lambda/\alpha$ conclusion are unchanged.
- **GQSP layer and polynomial class paragraph.** v4 mentioned the "minor modification of MW Thm. 6" for Laurent polynomials with negative powers. v5 drops that sentence entirely and emphasises that Cor. 5 delivers **ordinary** polynomials in $\mathbb{C}[z]$, with no $W^\dagger$ calls. Everything else in this paragraph (Thm. 3 angles, three-parameter opening rotation, Cor. 5 $\|P\|_\infty \leq 1$ condition) is unchanged.
- **Target polynomial paragraph.** v4 truncated Jacobi–Anger at $|k|\leq d$ and stopped at the Laurent polynomial (wrongly described as "realisable by GQSP"). v5 adds the $z^{d_{\text{JA}}}$ shift step to land in ordinary $\mathbb{C}[z]$ of degree $2d_{\text{JA}}$, derives the coefficient formula $c_m = (-\ii)^{m-d_{\text{JA}}} J_{m-d_{\text{JA}}}(\delta\alpha)$, and separately cites Cor. 5 as the realisability statement.
- **Degree-bound paragraph.** v4 had a single display "$d = \mathcal{O}(\delta\alpha + \log/\log\log)$" with $d$ ambiguous. v5 names the Jacobi–Anger truncation index $d_{\text{JA}}$, states (Eq. `eq:qsp-degree`) as a bound on $d_{\text{JA}}$, and separately notes that the polynomial degree $2d_{\text{JA}}$ and the number of controlled-$W$ calls is twice that — factor-of-2 made explicit.
- **Default $P_1$ paragraph.** v4 wrote $P_1 = J_0 - 2\ii J_1\cos\theta \approx 1 - \ii\delta\alpha\cos\theta$ (Laurent). v5 writes $P_2(z) = -\ii J_1\,z^2 + J_0\,z - \ii J_1$ (ordinary, degree 2, two controlled-$W$ calls, three single-qubit rotations), with the small-$\delta\alpha$ expansion $P_2(\ee^{\ii\theta}) \approx \ee^{\ii\theta}\,\ee^{-\ii\delta\alpha\cos\theta} + \mathcal{O}((\delta\alpha)^2)$.
- **Angle-finding paragraph.** v4 mentioned Haah 2019 and Dong–Meng–Whaley–Lin 2021 as historical benchmarks. v5 drops both, citing only MW (Alg. 1 for angle extraction) and BS (for $Q$-finding), with a one-sentence summary of what BS replaces.
- **Closing cost sentence.** v4: "$d$ applications of $U_B$ together with $d$ reflections $R_b$". v5: "$2d_{\text{JA}}$ applications of $W$ (each costing two PREP calls and one SELECT call, plus the PREP-sandwiched reflection)".

### Paragraphs unchanged

- Opening paragraph on the LCU $\to$ GQSP pipeline and $\alpha = \|b_-\|_1\,\|b_+\|_1$, with the minor addition of naming the LCU decomposition $B = \sum_j \alpha_j U_j$ and the PREP/SELECT operators explicitly (needed downstream for MW Cor. 8).
- Small-$\delta\alpha$ default regime discussion (the regime where $d_{\text{JA}} = 1$ suffices), with the correction that the default now implements a **degree-2** ordinary polynomial with **two** controlled-$W$ calls, not a Laurent $P_1$ with one layer.

### Citations dropped

- `lowChuang2019qubitization`: the walk and the scaling are now both attributed to MW.
- `haah2019product`, `dong2021efficient`: historical angle-finding references are not needed when BS is the chosen pre-processing tool and MW's Alg. 1 is the chosen extraction.
