# Quantum Signal Processing (v6 — Motlagh–Wiebe GQSP, Route 2: Thm. 6 Laurent via controlled-$W^\dagger$)

> **Insertion target:** `supplementary-informations/1_preliminaries.tex:1068`, body of `\subsection{Quantum Signal Processing.} \label{sec:prelim-QSP}` — between the `\subsection{...}\label{...}` line (line 1068) and `\begin{figure}[ht]` (line 1071). The figure (Fig. `circ:gqsp`, lines 1071–1113) sits immediately below the subsection heading; see Suggestions for the figure-redraw note, which is **required** for v6 to be self-consistent.
>
> **Change from v5:** commit to Route 2 (MW Thm. 6, Laurent via controlled-$W^\dagger$), dropping the ordinary-polynomial shift trick. This resolves the FATAL issue F1 flagged in `drafts/qsp-subsection-v5-review.md`: the shift trick $P_{2d_\text{JA}}(z) = z^{d_\text{JA}} L_{d_\text{JA}}(z)$ implements $W^{d_\text{JA}}\,\ee^{-\ii\delta B}$, not $\ee^{-\ii\delta B}$, because the trailing $W^{d_\text{JA}}$ factor is $\lambda$-dependent. The resolution is to invoke MW Thm. 6, which replaces $d_{\text{JA}}$ of the $2d_{\text{JA}}$ controlled-$W$ signal-operator slots by controlled-$W^\dagger$, putting $L_{d_{\text{JA}}}(W)$ — the Laurent polynomial directly — in the top-left block. The v5 walk-operator paragraph, the basis/matrix/eigenvalue derivation, the Bessel-tail analysis, and the angle-finding paragraph all carry over verbatim; what changes is the polynomial-class paragraph, the target-polynomial paragraph, and the cost accounting.
>
> **Existing labels referenced:** `sec:prelim-LCU`, `circ:gqsp`, `alg:coh`, `eq:B-block-encoding`, `eq:b_plus-s-eta`.
> **New labels introduced:** `eq:qsp-walk`, `eq:qsp-poly`, `eq:qsp-degree`. `eq:jacobi-anger` stays as in v4/v5.
> **New citation keys introduced:** `berntsonSunderhauf2025complementary` (Q-finding). All other keys are already in `references.bib`: `motlagh2024generalized` (now additionally used for Thm. 6), `chen2023efficient`.
> **Citation keys dropped relative to v4:** `lowChuang2019qubitization`, `haah2019product`, `dong2021efficient` — unchanged from v5.

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
Implementing $\ee^{-\ii\delta B}$ on the $\lambda$-eigenspace of $B$ therefore reduces to finding a polynomial function of $z$ that approximates $\ee^{-\ii\delta\alpha\cos\theta}$ on the unit circle $\mathbb{T} = \{z\in\mathbb{C} : |z|=1\}$, evaluated at $z = \ee^{\ii\theta}$.

*GQSP layer and polynomial class.* The circuit in (Fig. `circ:gqsp`) alternates $d$ (open-)controlled calls to $W$ with $d+1$ single-qubit rotations on the QSP ancilla $\ket{0}_\mathrm{QSP}$: a three-parameter opening rotation $R(\theta_0,\phi_0,\lambda_0)$ and $d$ further two-parameter rotations $R(\theta_k,\phi_k,0)$, $k = 1,\dots,d$ (the third parameter $\lambda_0 \in \mathbb{R}$ of the opening rotation, distinct from the spectral variable $\lambda$, absorbs a single global phase per [Motlagh & Wiebe 2024, Thm. 3] [CITE: motlagh2024generalized]). In the Hamiltonian-simulation application below, this generic circuit depth specialises to $d = 2d_{\text{JA}}$ (see the target-polynomial paragraph). Post-selecting both ancillas on $\ket{0}_\mathrm{QSP}\otimes(\text{PREP}\ket{0})$ returns, on each qubitization subspace, the top-left block of the composed circuit, which for the baseline framework of [Motlagh & Wiebe 2024, Thm. 3 and Cor. 5] [CITE: motlagh2024generalized] is an **ordinary** polynomial $P \in \mathbb{C}[z]$ of degree $\deg P \leq d$ satisfying
$$
\|P\|_{\infty,\mathbb{T}} \;:=\; \max_{z\in\mathbb{T}}\, |P(z)| \;\leq\; 1,
$$ <!-- \label{eq:qsp-poly} -->
evaluated at $z = \ee^{\ii\theta}$. This baseline class, advertised by MW's §III opening as "avoiding $U^\dagger$ and Laurent polynomials", is complete for amplitude-shaping applications (OAA, fixed-point search, error mitigation) where the target is genuinely a polynomial in $\ee^{\ii\theta}$. For Hamiltonian simulation, however, the target $\ee^{-\ii\delta\alpha\cos\theta}$ is a Laurent function of $\ee^{\ii\theta}$ not expressible as an ordinary polynomial in $\mathbb{C}[z]$ of any finite degree — an ordinary $P$ satisfying $P(\ee^{\ii\theta}) = P(\ee^{-\ii\theta})$ on $\mathbb{T}$ has $P(z) = P(1/z)$ on $\mathbb{T}$ and hence (matching Fourier coefficients) must be constant. MW Thm. 6 [Motlagh & Wiebe 2024, Thm. 6] [CITE: motlagh2024generalized] extends the class to Laurent polynomials $P'(z) = z^{-k} P(z)$ with $P \in \mathbb{C}[z]$, $\deg P \leq d$, and $0 \leq k \leq d$, by replacing $k$ of the $d$ controlled-$W$ calls in (Fig. `circ:gqsp`) by controlled-$W^\dagger$ calls — formally, the auxiliary signal operator $A' = \ket{0}\!\bra{0}\otimes I + \ket{1}\!\bra{1}\otimes W^\dagger$ is substituted for $A = \ket{0}\!\bra{0}\otimes I + \ket{1}\!\bra{1}\otimes W$ in $k$ of the slots, and the remaining angles are left unchanged (MW Eqs. 45–53).

*Target polynomial via Jacobi–Anger.* The target $\ee^{-\ii\delta\alpha\cos\theta}$ admits the bilateral Jacobi–Anger expansion [Motlagh & Wiebe 2024, Eq. (62)] [CITE: motlagh2024generalized]
$$
\ee^{-\ii\delta\alpha\cos\theta} \;=\; \sum_{k=-\infty}^{\infty} (-\ii)^k\,J_k(\delta\alpha)\,\ee^{\ii k\theta}, \qquad \theta \in \mathbb{R},
$$ <!-- \label{eq:jacobi-anger} -->
with Bessel coefficients; the sign pattern follows from $t \to -\delta\alpha$ in MW's identity $\ee^{\ii t\cos\theta} = \sum_k \ii^k J_k(t)\,\ee^{\ii k\theta}$ together with $J_n(-z) = (-1)^n J_n(z)$, giving $\ii^n (-1)^n = (-\ii)^n$ directly — no invocation of $J_{-k}$ needed. Truncating at $|k| \leq d_{\text{JA}}$ yields the Laurent polynomial
$$
L_{d_{\text{JA}}}(\ee^{\ii\theta}) \;:=\; \sum_{k=-d_{\text{JA}}}^{d_{\text{JA}}} (-\ii)^k\,J_k(\delta\alpha)\,\ee^{\ii k\theta} \;\in\; \mathbb{C}[z,z^{-1}],
$$
with uniform tail bound $\|\ee^{-\ii\delta\alpha\cos\theta} - L_{d_{\text{JA}}}\|_{\infty,\mathbb{T}} \leq 2\sum_{|k|>d_{\text{JA}}}|J_k(\delta\alpha)| \leq \mathcal{O}\bigl((\ee\delta\alpha/(2(d_{\text{JA}}+1)))^{d_{\text{JA}}+1}\bigr)$ from the rapid Bessel decay [Motlagh & Wiebe 2024, Eq. (63)] [CITE: motlagh2024generalized]. To realise $L_{d_{\text{JA}}}$ via MW Thm. 6, introduce the underlying *ordinary* polynomial
$$
P(z) \;:=\; z^{d_{\text{JA}}}\, L_{d_{\text{JA}}}(z) \;=\; \sum_{m=0}^{2d_{\text{JA}}} c_m\,z^m \;\in\; \mathbb{C}[z], \qquad c_m \;=\; (-\ii)^{m-d_{\text{JA}}}\,J_{m-d_{\text{JA}}}(\delta\alpha),
$$
of degree $2d_{\text{JA}}$. This $P$ is not the operator we implement; it is the bookkeeping polynomial that MW Alg. 1 runs on. Since multiplication by $z^{d_{\text{JA}}}$ is unit-modulus on $\mathbb{T}$, $\|P\|_{\infty,\mathbb{T}} = \|L_{d_{\text{JA}}}\|_{\infty,\mathbb{T}} \leq 1 + \mathcal{O}((\ee\delta\alpha/(2(d_{\text{JA}}+1)))^{d_{\text{JA}}+1})$, and a rescaling $(1-\varepsilon_{\text{QSP}})\,P$ enforces (Eq. `eq:qsp-poly`) strictly without changing the Laurent polynomial that is implemented beyond the same rescaling factor. Applying MW Thm. 6 with $d = 2d_{\text{JA}}$ and $k = d_{\text{JA}}$ converts $P$ to
$$
P'(z) \;=\; z^{-d_{\text{JA}}}\,P(z) \;=\; L_{d_{\text{JA}}}(z),
$$
so the post-selected top-left block directly realises the Jacobi–Anger Laurent truncation. On the qubitization subspace this evaluates to $L_{d_{\text{JA}}}(W_\lambda) = \ee^{-\ii\delta\alpha\cos\theta_\lambda}\,I_{2\times 2} + \mathcal{O}(\varepsilon_{\text{QSP}}) = \ee^{-\ii\delta\lambda}\,I_{2\times 2} + \mathcal{O}(\varepsilon_{\text{QSP}})$, reproducing MW Eq. (66) with no trailing $W^{d_{\text{JA}}}$ factor.

*Degree bound and cost.* The Jacobi–Anger truncation index that suffices for $\varepsilon_{\text{QSP}}$-accuracy on $\mathbb{T}$ is, per [Motlagh & Wiebe 2024, Thm. 7] [CITE: motlagh2024generalized],
$$
d_{\text{JA}} \;=\; \mathcal{O}\!\left(\delta\alpha \;+\; \frac{\log(1/\varepsilon_\mathrm{QSP})}{\log\log(1/\varepsilon_\mathrm{QSP})}\right).
$$ <!-- \label{eq:qsp-degree} -->
The GQSP circuit makes $d_{\text{JA}}$ controlled-$W$ calls and $d_{\text{JA}}$ controlled-$W^\dagger$ calls — $2d_{\text{JA}}$ controlled operations in total — matching MW Cor. 8's $\mathcal{O}(\alpha\delta + \log(1/\varepsilon)/\log\log(1/\varepsilon))$ bound on applications of $W$ and $W^\dagger$. Controlled-$W^\dagger$ costs the same as controlled-$W$: both are implemented as an open-controlled call to the same $2m+n$-qubit block-encoding oracle, with a conjugation applied at the PREP-sandwich level of $W$; no hardware premium arises from the Thm. 6 substitution. The circuit also uses $2d_{\text{JA}}+1$ single-qubit rotations on the QSP ancilla (three parameters on the opening rotation, two on each of the remaining $2d_{\text{JA}}$), contributing $\mathcal{O}(d_{\text{JA}})$ two-qubit gates via standard synthesis.

*Small-$\delta\alpha$ default.* In the regime of Chapter 5, the coherent step size $\delta$ is small and $\delta\alpha \lesssim 1$ is the norm. Already $d_{\text{JA}} = 1$ gives the Laurent polynomial
$$
L_1(\ee^{\ii\theta}) \;=\; J_0(\delta\alpha) \,-\, 2\ii\,J_1(\delta\alpha)\,\cos\theta,
$$
of bilateral degree $1$ (the underlying ordinary polynomial $P(z) = zL_1(z) = -\ii J_1\,z^2 + J_0\,z - \ii J_1$ has degree $2$). The Thm. 6 circuit at $d = 2$, $k = 1$ realises this with **one controlled-$W$ call, one controlled-$W^\dagger$ call, and three single-qubit rotations** $(R_0, R_1, R_2)$ on the QSP ancilla. Under the small-argument Bessel expansions $J_0(z) = 1 + \mathcal{O}(z^2)$ and $2J_1(z) = z + \mathcal{O}(z^3)$,
$$
L_1(\ee^{\ii\theta}) \;=\; 1 \,-\, \ii\,\delta\alpha\,\cos\theta \,+\, \mathcal{O}((\delta\alpha)^2) \;=\; \ee^{-\ii\delta\alpha\cos\theta} \,+\, \mathcal{O}((\delta\alpha)^2),
$$
**exactly** — no residual $\ee^{\ii\theta}$ factor. Two regimes force $d_{\text{JA}} \geq 2$: (a) when $\alpha$ is large — e.g. the Metropolis kernel $b_+^{(s,\eta)}$ (Eq. `eq:b_plus-s-eta`) with $\ell_1$-norm scaling as $\log(\beta\|H\|\,\|\!\sum_a A_a^\dagger A_a\|/\varepsilon)$ (see `2_methods.tex:1726`), so that $\delta\alpha$ is no longer small; (b) when $\delta$ is coarsened to reduce the total number of coherent steps at fixed overall accuracy. In both cases $d_{\text{JA}}$ grows linearly in $\delta\alpha$ per (Eq. `eq:qsp-degree`).

*Angle finding.* Given the coefficients $\{c_m\}_{m=0}^{2d_{\text{JA}}}$ of the underlying ordinary polynomial $P(z) = z^{d_{\text{JA}}}\,L_{d_{\text{JA}}}(z) \in \mathbb{C}[z]$, the GQSP angle vector $\{(\theta_0,\phi_0,\lambda_0)\} \cup \{(\theta_k,\phi_k)\}_{k=1}^{2d_{\text{JA}}}$ is extracted in two steps. (The angles are computed from $P$, the ordinary polynomial of degree $2d_{\text{JA}}$; but the $2d_{\text{JA}}$ signal-operator slots of the circuit are deployed as $d_{\text{JA}}$ controlled-$W$ calls followed by $d_{\text{JA}}$ controlled-$W^\dagger$ calls per MW Thm. 6 — the same angle vector, a different slot assignment, per MW Eqs. 45–53.) First, a *complementary polynomial* $Q \in \mathbb{C}[z]$ of degree $\deg Q = \deg P = 2d_{\text{JA}}$ satisfying
$$
|P(z)|^2 + |Q(z)|^2 = 1 \qquad\text{on } \mathbb{T}
$$
is constructed; its existence is [Motlagh & Wiebe 2024, Thm. 4] [CITE: motlagh2024generalized]. Second, the exact recursive extraction of [Motlagh & Wiebe 2024, §IV, Alg. 1] [CITE: motlagh2024generalized] reads off the $(2d_{\text{JA}}+1)$ rotation triples from $(P, Q)$ in $\mathcal{O}(d_{\text{JA}}^2)$ arithmetic operations. For the $Q$-finding step, MW propose an FFT-based convolution objective (Eqs. (58)–(60)) solved by nonlinear optimisation. We replace that step with the contour-integral / FFT construction of [Berntson & Sünderhauf 2025] [CITE: berntsonSunderhauf2025complementary], which evaluates $Q$ in the monomial basis in classical time $\tilde{\mathcal{O}}(d_{\text{JA}})$ with explicit a priori error bounds and has been numerically demonstrated up to $\deg P \sim 10^7$. This is a practical gain over MW's optimisation route — no solver initialisation, no convergence tuning — and delivers the same coefficient vector $\{b_m\}$. The MW angle-extraction recursion (Alg. 1) then proceeds unchanged.

The unitary $\ee^{-\ii\delta B}$ produced by (Fig. `circ:gqsp`) closes the block-encoding pipeline begun in (§`sec:prelim-LCU`): it is the coherent step of (Alg. `alg:coh`). Its cost is $d_{\text{JA}}$ controlled-$W$ calls and $d_{\text{JA}}$ controlled-$W^\dagger$ calls plus $\mathcal{O}(d_{\text{JA}})$ single-qubit gates; with $d_{\text{JA}}$ per (Eq. `eq:qsp-degree`), the coherent step contributes a $\tilde{\mathcal{O}}(\delta\alpha)$ factor per $\delta$-step of the Gibbs sampler.

---

## Citations

- **`berntsonSunderhauf2025complementary`** — Berntson, B. K. & Sünderhauf, C., *"Complementary Polynomials in Quantum Signal Processing"*, **Commun. Math. Phys. 406, 161 (2025)**, arXiv:2406.04246. Contour-integral / FFT construction of $Q_d$ in classical time $\tilde{\mathcal{O}}(d)$ with explicit error bounds; numerically demonstrated to $\deg P \sim 10^7$. Replaces MW's §IV FFT-optimisation step for complex $P$.
- **`motlagh2024generalized`** — already in `references.bib`. Thm. 3 and Cor. 5 characterise the baseline GQSP polynomial class ($\mathbb{C}[z]$, $\|P\|_{\infty,\mathbb{T}}\leq 1$). **Thm. 6 (Eqs. 45–53)** extends the class to Laurent polynomials $P'(z) = z^{-k}P(z)$ by replacing $k$ of the signal operators in the baseline circuit with controlled-$U^\dagger$; this is the construction used in v6. Thm. 7 gives the Jacobi–Anger degree bound. Cor. 8 gives the PREP-sandwiched walk $W = -(I - 2\,\text{PREP}\ket{0}\bra{0}\text{PREP}^\dagger)\text{SELECT}$ and its $2{\times}2$ block action $W_\lambda$ (Eq. 64); its cost bound is $\mathcal{O}(\alpha\delta + \log/\log\log)$ applications of $W$ (reading "applications of $W$" as encompassing both $W$ and $W^\dagger$ under Thm. 6). Algorithm 1 in §IV gives the $\mathcal{O}(d^2)$ recursive angle extraction, applied to the underlying ordinary polynomial $P$ of degree $2d_{\text{JA}}$.

```bibtex
@article{berntsonSunderhauf2025complementary,
  author        = {Berntson, Bjorn K. and S\"underhauf, Christoph},
  title         = {Complementary Polynomials in Quantum Signal Processing},
  journal       = {Commun. Math. Phys.},
  volume        = {406},
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

### Why Thm. 6 is needed for v6

MW §III opens with *"our approach eschews the need to use $U^\dagger$ and avoids the use of Laurent polynomials in the analysis"*, and v5 took that pitch literally: it used only Cor. 5 and tried to realise the Jacobi–Anger Laurent expansion via the shift trick $P_{2d_{\text{JA}}}(z) = z^{d_{\text{JA}}}L_{d_{\text{JA}}}(z)$. The review `drafts/qsp-subsection-v5-review.md` shows that construction is incorrect: on the qubitization subspace, an ordinary $P \in \mathbb{C}[z]$ acts as $\operatorname{diag}(P(\ee^{\ii\theta_\lambda}), P(\ee^{-\ii\theta_\lambda}))$, so matching MW Eq. (66) (with equal diagonal entries $\ee^{-\ii\lambda\delta}$) would force $P(z) = P(1/z)$ on $\mathbb{T}$. Matching Fourier coefficients of $P(z) = \sum_{m=0}^{2d} c_m z^m$ and $P(1/z) = \sum_{m=0}^{2d} c_m z^{-m}$ on $\mathbb{T}$ gives $c_m = 0$ for all $m \geq 1$ — i.e. $P$ must be constant. Hence **no non-constant ordinary polynomial in $\mathbb{C}[z]$ realises MW Cor. 8's target** within pure Cor. 5.

The §III pitch is therefore accurate about the *baseline* framework (Thm. 3 / Cor. 5 / Thm. 4 — the applications that genuinely need only a polynomial phase such as OAA, fixed-point search, amplitude shaping), but Hamiltonian simulation via Cor. 8 is exactly the canonical place where Thm. 6 (and hence controlled-$U^\dagger$) re-enters. MW's own Cor. 8 proof moves from Eq. (65) to Eq. (66) by choosing "$P(U) = \ee^{-\ii\alpha t H}$" — i.e. they silently invoke the Thm. 6 extension via Thm. 7's Laurent approximation. v6 makes this explicit: Thm. 6 is the construction, $d_{\text{JA}}$ of the $2d_{\text{JA}}$ signal slots become controlled-$W^\dagger$, and the circuit's top-left block is $L_{d_{\text{JA}}}(W) = \ee^{-\ii\delta B} + \mathcal{O}(\varepsilon_{\text{QSP}})$ directly, matching Eq. (66) with no residual factor.

### The shift-by-$z^{d_{\text{JA}}}$ step (resolved and closed)

The v5 `[CHECK]` item "The shift-by-$z^{d_{\text{JA}}}$ step and the Cor. 8 identification" is resolved negatively: the shift trick is **not** the right construction. The review's §"[CHECK] adjudication" verifies the derivation given above (Fourier-coefficient matching forces $P$ constant), and the recommendation is Route 2 — MW Thm. 6 + controlled-$W^\dagger$. v6 implements that recommendation. The alternative Route 1 (two GQSP calls with $\cos + \ii\sin$) does not rescue ordinary polynomials either: each branch $\cos(\delta\alpha\cos\theta)$ and $\sin(\delta\alpha\cos\theta)$ is even in $\theta$, so the same $P(z) = P(1/z)$ obstruction applies per branch, and Laurent polynomials would still be needed per branch (at double the total cost of Route 2). Route 2 is therefore the only option that delivers MW Cor. 8's advertised cost and a clean $\ee^{-\ii\delta B}$ output.

### Thm. 6 circuit-level identity

MW Eqs. 48–53 give the operator identity underlying Thm. 6. Starting from Eq. (49),
$$
\biggl[\prod_{j=1}^{k} R(\theta_{d-k+j}, \phi_{d-k+j}, 0)\,A'\biggr]\,\biggl[\prod_{j=1}^{d-k} R(\theta_j, \phi_j, 0)\,A\biggr]\,R(\theta_0, \phi_0, \lambda),
$$
the substitution $A' = (I \otimes U^\dagger)\,A$ (MW Eq. 48), together with the commutation $(I \otimes U^\dagger) R(\theta, \phi, 0) = R(\theta, \phi, 0)(I \otimes U^\dagger)$ and $(I \otimes U^\dagger) A = A\,(I \otimes U^\dagger)$ (the first because the $\mathrm{SU}(2)$ rotation acts only on the QSP ancilla; the second because $U^\dagger U = UU^\dagger = I$ on the target register and $\ket{k}\!\bra{k}$ are mutually orthogonal projectors, so $(I\otimes U^\dagger)(\ket{0}\!\bra{0}\otimes I + \ket{1}\!\bra{1}\otimes U) = \ket{0}\!\bra{0}\otimes U^\dagger + \ket{1}\!\bra{1}\otimes I = (\ket{0}\!\bra{0}\otimes I + \ket{1}\!\bra{1}\otimes U)(I\otimes U^\dagger)$), moves $(I \otimes U^\dagger)^k$ to the left of the entire rotation–signal product, giving Eqs. (52)–(53): the $2{\times}2$ operator matrix in the $\ket{0}_\mathrm{QSP}, \ket{1}_\mathrm{QSP}$ basis has its original $P(U)$ top-left block replaced by $U^{-k} P(U) = P'(U)$, and the lower row similarly. Thus the same $(\theta_j, \phi_j)$ vector deployed with $k$ of the $d$ slots as $A'$ instead of $A$ realises $P'(U) = U^{-k} P(U)$ in the top-left block. For our application $U \mapsto W$, $d = 2d_{\text{JA}}$, $k = d_{\text{JA}}$, and $P'(z) = z^{-d_{\text{JA}}} z^{d_{\text{JA}}} L_{d_{\text{JA}}}(z) = L_{d_{\text{JA}}}(z)$. The angle vector, computed from $P$ via Alg. 1, is identical between the shift-trick circuit (which implements $P(W) = W^{d_{\text{JA}}} L_{d_{\text{JA}}}(W)$) and the Thm. 6 circuit (which implements $P'(W) = L_{d_{\text{JA}}}(W)$); only the slot assignments differ.

### Factor of 2 in the polynomial degree (updated for v6)

The "degree" quoted in MW Thm. 7 and in (Eq. `eq:qsp-degree`) is the Jacobi–Anger truncation index $d_{\text{JA}}$. The *underlying ordinary polynomial* $P = z^{d_{\text{JA}}}L_{d_{\text{JA}}}$ on which Alg. 1 runs has degree $2d_{\text{JA}}$ (and hence $2d_{\text{JA}}+1$ rotation triples come out of Alg. 1). The *Laurent polynomial actually implemented*, $L_{d_{\text{JA}}}$, has bilateral degree $d_{\text{JA}}$: it runs over $k \in [-d_{\text{JA}}, d_{\text{JA}}]$. The number of controlled operations in (Fig. `circ:gqsp`) is $2d_{\text{JA}}$: $d_{\text{JA}}$ calls to $W$ plus $d_{\text{JA}}$ calls to $W^\dagger$, matching the polynomial degree $d = 2d_{\text{JA}}$ of $P$ slot-for-slot. This is the honest cost accounting; MW Cor. 8's statement "$\mathcal{O}(\alpha t + \log(1/\varepsilon)/\log\log(1/\varepsilon))$ applications of $W$" reads correctly if "applications of $W$" is understood to cover both $W$ and $W^\dagger$ under the Thm. 6 substitution.

### Derivation of $L_1$ (the $d_{\text{JA}} = 1$ default)

Starting from (Eq. `eq:jacobi-anger`) truncated at $|k| \leq 1$: $L_1(\ee^{\ii\theta}) = (-\ii)^{-1} J_{-1}\,\ee^{-\ii\theta} + J_0 + (-\ii)^1 J_1\,\ee^{\ii\theta}$. Using $J_{-1}(z) = -J_1(z)$ and $(-\ii)^{-1} = \ii$: $L_1 = \ii\,(-J_1)\,\ee^{-\ii\theta} + J_0 - \ii J_1\,\ee^{\ii\theta} = -\ii J_1\,(\ee^{-\ii\theta}+\ee^{\ii\theta}) + J_0 = J_0 - 2\ii J_1\cos\theta$. Under the small-argument expansions $J_0(z) = 1 + \mathcal{O}(z^2)$ and $2J_1(z) = z + \mathcal{O}(z^3)$, one gets $L_1 = 1 - \ii\delta\alpha\cos\theta + \mathcal{O}((\delta\alpha)^2) = \ee^{-\ii\delta\alpha\cos\theta} + \mathcal{O}((\delta\alpha)^2)$ directly. For comparison, the underlying ordinary polynomial $P(z) = z L_1(z) = -\ii J_1\,z^2 + J_0\,z - \ii J_1$ has coefficients $(c_0, c_1, c_2) = (-\ii J_1, J_0, -\ii J_1)$ and degree $2$; Alg. 1 runs on $P$ to produce $(\theta_0, \phi_0, \lambda_0), (\theta_1, \phi_1), (\theta_2, \phi_2)$ — three rotation triples — then per MW Thm. 6's canonical "last $k$" pattern (Eq. 46) with $d = 2, k = 1$, the Thm. 6 circuit deploys slot 1 as controlled-$W$ and slot 2 as controlled-$W^\dagger$, giving $P'(z) = z^{-1} P(z) = L_1(z)$ in the top-left block. MW's body text above Thm. 6 ("replacing any $k$ instances of $A$") asserts that any distribution of the $k$ flipped slots realises the same $P'$, though MW's proof covers only the canonical pattern. All signs cross-checked against DLMF 10.2.2, 10.12.1, and MW Eqs. (62), (46).

### Scope and length

- Six body paragraphs + three numbered displays (`eq:qsp-walk`, `eq:qsp-poly`, `eq:qsp-degree`) + defensively-labelled `eq:jacobi-anger` + one inline $W_\lambda$ matrix + one inline $P$ display + one inline $L_1$ display. ~1.3–1.5 pages of typeset body, matching v5.
- Equation labels `eq:qsp-walk`, `eq:qsp-poly`, `eq:qsp-degree` are all candidates for Chapter 5 cost-analysis cross-references.

### Notation choices

- $\lambda \in [-\alpha,\alpha]$ denotes the eigenvalue of $B$ in physical units (not $B/\alpha$); MW Eq. 64 has $\lambda/\alpha$ everywhere and Eq. 66 has $\ee^{-\ii\lambda t}$, matching.
- $d_{\text{JA}}$ is the Jacobi–Anger truncation index; the underlying ordinary polynomial has degree $2d_{\text{JA}}$; the Laurent polynomial actually implemented has bilateral degree $d_{\text{JA}}$; the circuit uses $2d_{\text{JA}}$ controlled operations total ($d_{\text{JA}}$ of $W$, $d_{\text{JA}}$ of $W^\dagger$).
- $\mathbb{T} = \{|z|=1\}$ for the complex unit circle, matching MW and Berntson–Sünderhauf.
- $\text{PREP}\ket{0}$ is the specific superposition $\sum_j \sqrt{\alpha_j/\alpha}\,\ket{j}$ from the LCU section; the "reflection about $\text{PREP}\ket{0}$" $2\,\text{PREP}\ket{0}\bra{0}\text{PREP}^\dagger - I$ is distinct from the raw-ancilla reflection $R_b = 2\ket{0^b}\bra{0^b} - I$ and should not be conflated.
- $\ee$, $\ii$ respected throughout; ancilla labels $\ket{0}_\mathrm{QSP}$, $\ket{0^b}$ match LCU v2 and (Fig. `circ:gqsp`) after its redraw.

### Angle-finding complexity

- The complementary-polynomial step is $\tilde{\mathcal{O}}(d_{\text{JA}})$ via Berntson–Sünderhauf (no optimisation solver required); the MW angle-extraction recursion (§IV, Alg. 1) is $\mathcal{O}(d_{\text{JA}}^2)$ arithmetic operations in exact arithmetic (run on the degree-$2d_{\text{JA}}$ ordinary polynomial $P$, not on the Laurent $L_{d_{\text{JA}}}$). The body states both explicitly. A single historical-context sentence (MW's §IV comment that traditional QSP angle finding topped out at $\sim 10^4$ degrees, vs BS demonstrating $\sim 10^7$) is worth adding if Chapter 5 later needs a low-temperature benchmark; v6 omits it for length.

---

## Suggestions

Items that would polish the subsection but require either a user decision or lie slightly outside v6's scope:

1. **[FIX-NEXT] Redraw `circ:gqsp` (required for v6 self-consistency).** The figure at `1_preliminaries.tex:1071–1113` currently draws $W = R_b\,U_B$ as the sequence "bundled $U_B$ block" followed by a "box labelled $R_b$" — an external reflection about $\ket{0^b}$. Under v6, the walk is MW Cor. 8's PREP-sandwiched form, and *half the signal-operator slots must be flipped in orientation* to realise the Thm. 6 Laurent extension. Two ways to fix:
   - *Option A (honest; what v6's body assumes):* Expand the $U_B$ block into its PREP / SELECT / PREP$^\dagger$ composition, drop the outer PREP$^\dagger$ and PREP that would make a full $U_B$, and instead draw, per GQSP layer: (i) one open-controlled SELECT (or its dagger in the flipped slots); (ii) the PREP-sandwiched reflection $2\,\text{PREP}\ket{0}\bra{0}\text{PREP}^\dagger - I$ drawn as a box `PREP$^\dagger$, reflection about $\ket{0^b}$, PREP`; (iii) the $R(\theta_k,\phi_k,0)$ rotation on the QSP ancilla. For v6, **$d_{\text{JA}}$ of the $2d_{\text{JA}}$ signal-operator slots must be drawn with controlled-$W^\dagger$ (e.g. the last $d_{\text{JA}}$ slots, matching MW Eqs. 45–53's $j = d-k+1, \dots, d$ index range)**. This is the **recommended redraw** — it exposes the Cor. 8 + Thm. 6 structure explicitly.
   - *Option B (bundled; acceptable but less informative):* Keep the $U_B$ block and annotate the signal-operator boxes as "controlled-$W$ / controlled-$W^\dagger$" with the split $d_{\text{JA}}$–$d_{\text{JA}}$ indicated in the caption. The picture still hides the Cor. 8 construction, and the Thm. 6 substitution is only visible in the caption.
   Either way, the caption text "*degree-$d$ Chebyshev polynomial $P(W)\approx\ee^{-\ii\delta B}$*" (lines 1110–1111) **must be rewritten** in both options — GQSP produces the **Laurent polynomial $L_{d_{\text{JA}}}(z) \in \mathbb{C}[z,z^{-1}]$ of bilateral degree $d_{\text{JA}}$**, not Chebyshev polynomials, via MW Thm. 6. Suggested caption: *"Degree-$2d_{\text{JA}}$ GQSP realising $\ee^{-\ii\delta B}$ from the walk operator $W$ of (Eq.~\ref{eq:qsp-walk}) (MW Cor. 8), with $k = d_{\text{JA}}$ of the signal-operator calls promoted to controlled-$W^\dagger$ per MW Thm. 6. Post-selecting both ancillas on $\ket{0}_\mathrm{QSP}\otimes\text{PREP}\ket{0}$ yields the Laurent polynomial $L_{d_{\text{JA}}}(W)$ realising $\ee^{-\ii\delta B}$ on the qubitization subspace."*

2. **[FIX-NEXT] Reconcile `2_methods.tex` with v6.** The methods chapter must speak the same convention. Under v6, `2_methods.tex` needs:
   - **Line 1609** (algorithm pseudocode, inside `\textsc{CoherentStep}`): "$P(\ee^{\ii\theta}) \approx \ee^{-\ii\delta\alpha\sin\theta}$" $\to$ "$L_{d_{\text{JA}}}(\ee^{\ii\theta}) \approx \ee^{-\ii\delta\alpha\cos\theta}$ (Laurent polynomial, bilateral degree $d_{\text{JA}}$, MW Thm. 6)", and update the `\Comment{Degree $d=1$ calls to $U_{B_a}$}` annotation to `\Comment{$d_{\text{JA}}=1$ → 1 call to $W$ + 1 call to $W^\dagger$}`.
   - **Line 1731**: "its eigenvalues come in conjugate pairs $\ee^{\pm\ii\theta}$ where $\sin\theta = \lambda$ for each eigenvalue $\lambda$ of $B_a/\alpha$" $\to$ "its eigenvalues come in conjugate pairs $\ee^{\pm\ii\theta}$ where $\cos\theta = \lambda/\alpha$ for each eigenvalue $\lambda$ of $B_a$". Same line: "$P(\ee^{\ii\theta}) \approx \ee^{-\ii\delta\alpha\sin\theta}$" $\to$ "$L_{d_{\text{JA}}}(\ee^{\ii\theta}) \approx \ee^{-\ii\delta\alpha\cos\theta}$ on $\mathbb{T}$ (Laurent polynomial)". Also: bring in the walk operator $W$ explicitly, with text pointing to (Eq. `eq:qsp-walk`) in the preliminaries subsection.
   - **Lines 1733–1735** (degree bound): the statement $d = \mathcal{O}(\delta\alpha + \log(1/\varepsilon_{\text{QSP}})/\log\log(1/\varepsilon_{\text{QSP}}))$ is correct if "$d$" means $d_{\text{JA}}$. Rename to $d_{\text{JA}}$ and add a parenthetical "(the GQSP polynomial degree, i.e. the count of signal-operator slots, is $2d_{\text{JA}}$ — half of them are controlled-$W^\dagger$ per MW Thm. 6)".
   - **Line 1737** (default $d = 1$): "degree $d = 1$ is sufficient for the QSP" $\to$ "Jacobi–Anger truncation $d_{\text{JA}} = 1$ is sufficient, giving a Laurent polynomial of bilateral degree 1 and one controlled-$W$ call + one controlled-$W^\dagger$ call per coherent step".
   - **Cost accounting:** the $R_b$ and "$U_B$" references in the methods cost table should become "$W$-calls" and "$W^\dagger$-calls", with the split $d_{\text{JA}}$ + $d_{\text{JA}}$ per coherent step; one $W$ (or $W^\dagger$) = one SELECT (or SELECT$^\dagger$) + one PREP-sandwiched reflection + $\mathcal{O}(1)$ single-qubit gates.
   - **Note on $B_a$ vs $\tilde B_a$:** `2_methods.tex:1681`'s quantikz figure `circ:U_B_a` produces $\tilde{B}_a$ (the Trotterized version of $B_a$), while line 1700 (CoherentStep text) refers to "the per-jump Hermitian term $B_a$" (the exact version from `eq:B-bohr-domain`). The v6 preliminaries subsection stays with the generic $B$ (no tilde). The methods chapter will have to pick a convention: either use $\tilde{B}_a$ consistently in the circuit + CoherentStep, or add a sentence clarifying that the coherent step block-encodes $B_a$ up to Trotter error, producing $\ee^{-\ii\delta\tilde{B}_a} = \ee^{-\ii\delta B_a} + \mathcal{O}(\text{Trotter})$.
   - The user memory `.claude-memory/todo_methods_gqsp_reflection.md` tracks this reconciliation; once the methods edit lands, that memo can be updated to reflect v6's punch list (not v4/v5's). The `.claude-memory/thesis_sign_convention_trotter.md` note is unrelated and stays open.

3. **Explicit $d_{\text{JA}} = 1$ GQSP angle formulas.** For the default, the five rotation parameters $(\theta_0, \phi_0, \lambda_0, \theta_1, \phi_1, \theta_2, \phi_2)$ of a degree-2 GQSP circuit reduce to closed-form expressions in $J_0(\delta\alpha)$ and $J_1(\delta\alpha)$ via the MW Alg. 1 recursion applied to $P(z) = -\ii J_1\,z^2 + J_0\,z - \ii J_1$ (the underlying ordinary polynomial) and its complementary $Q$. Listing these explicitly would make the "default in practice" statement self-verifying. Requires separate verification against the numerical implementation; deferred. **Note the $2d_{\text{JA}}+1 = 3$ rotation triples** (seven real parameters with the opening $\lambda_0$), and the slot-assignment choice for which of the two signal-operator slots becomes controlled-$W^\dagger$.

4. **Label hooks in Chapter 5.** If the user adds `\label{sec:B-nested-LCU}` at `2_methods.tex:1707` and `\label{sec:B-GQSP}` at `2_methods.tex:1730`, the closing sentence of v6 could point more precisely than `(Alg. alg:coh)`. One-line user edit; not a blocker.

5. **"GQSP" vs "QSP" subsection title.** The existing label `sec:prelim-QSP` and the `\subsection{Quantum Signal Processing.}` heading at `1_preliminaries.tex:1068` use "QSP" as an umbrella. Everything in v6 is strictly GQSP. Consider retitling to `\subsection{Generalized Quantum Signal Processing.}` with label `sec:prelim-GQSP`, and update back-references in `2_methods.tex`. Low priority; current choice is defensible.

---

## Diff from v5

### Paragraphs rewritten

- **GQSP layer and polynomial class paragraph.** v5 said "Note that the accessible class is ordinary polynomials, not Laurent — no negative powers of $\ee^{\ii\theta}$ appear, and the construction uses $W$ only, never $W^\dagger$." v6 adds: for the *baseline* framework this is true (Thm. 3, Cor. 5); for Hamiltonian simulation, MW Thm. 6 extends the class to Laurent polynomials $P'(z) = z^{-k}P(z)$ by replacing $k$ of the $d$ signal slots with controlled-$W^\dagger$, and states the circuit-level identity $P' = z^{-k}P$. Also adds the "In our application this generic $d$ specialises to $d = 2d_{\text{JA}}$" clarification per review item 10.

- **Target polynomial paragraph.** v5 used the shift trick: truncate Jacobi–Anger at $|k| \leq d_{\text{JA}}$, multiply by $z^{d_{\text{JA}}}$ to land in $\mathbb{C}[z]$ of degree $2d_{\text{JA}}$, realise via Cor. 5. v6 drops the shift trick: truncate at $|k| \leq d_{\text{JA}}$ to get the Laurent $L_{d_{\text{JA}}}$, introduce the *underlying* ordinary polynomial $P = z^{d_{\text{JA}}}L_{d_{\text{JA}}}$ as an Alg. 1 bookkeeping object, then apply MW Thm. 6 with $d = 2d_{\text{JA}}$, $k = d_{\text{JA}}$ to realise $P' = z^{-d_{\text{JA}}}P = L_{d_{\text{JA}}}$ directly in the top-left block. Bessel tail bound tightened from $(\ee\delta\alpha/(2d_{\text{JA}}))^{d_{\text{JA}}}$ to $(\ee\delta\alpha/(2(d_{\text{JA}}+1)))^{d_{\text{JA}}+1}$ per review m1. Jacobi-Anger sign attribution tightened per m2 (direct use of $J_n(-z)$, no separate $J_{-k}$ invocation).

- **Degree bound / cost paragraph.** v5: "$2d_{\text{JA}}$ controlled calls to $W$ ... matching MW's Cor. 8 bound." v6: "$d_{\text{JA}}$ controlled-$W$ calls and $d_{\text{JA}}$ controlled-$W^\dagger$ calls — $2d_{\text{JA}}$ controlled operations in total — matching MW Cor. 8's $\mathcal{O}(\alpha\delta + \log/\log\log)$ bound." Also notes that controlled-$W^\dagger$ costs the same as controlled-$W$.

- **Small-$\delta\alpha$ default paragraph.** v5: $P_2(z) = -\ii J_1\,z^2 + J_0\,z - \ii J_1$ (ordinary, shift-trick, with trailing $\ee^{\ii\theta}$ factor). v6: $L_1(\ee^{\ii\theta}) = J_0 - 2\ii J_1\cos\theta$ (Laurent, implemented via Thm. 6 with 1 controlled-$W$ + 1 controlled-$W^\dagger$), giving $L_1 \approx \ee^{-\ii\delta\alpha\cos\theta} + \mathcal{O}((\delta\alpha)^2)$ **exactly**, no residual $\ee^{\ii\theta}$ factor.

- **Angle-finding paragraph.** v5 and v6 both use the same $(2d_{\text{JA}}+1)$ rotation triples computed by MW Alg. 1 on the underlying ordinary polynomial $P$ of degree $2d_{\text{JA}}$. v6 adds a parenthetical clarifying that the $2d_{\text{JA}}$ signal-operator slots are deployed as $d_{\text{JA}}$ controlled-$W$ calls followed by $d_{\text{JA}}$ controlled-$W^\dagger$ calls per MW Thm. 6 (Eqs. 45–53); the angle vector is identical to v5's, only the slot orientations differ.

- **Closing cost sentence.** v5: "$2d_{\text{JA}}$ applications of $W$ ..." v6: "$d_{\text{JA}}$ controlled-$W$ calls and $d_{\text{JA}}$ controlled-$W^\dagger$ calls plus $\mathcal{O}(d_{\text{JA}})$ single-qubit gates".

### Paragraphs unchanged

- **Walk-operator paragraph** — verbatim from v5, including the $W_\lambda$ matrix, the trace/det argument, and the eigenvalue conclusion $\cos\theta = \lambda/\alpha$.
- **Opening paragraph** on the LCU → GQSP pipeline, the $\alpha = \|b_-\|_1\|b_+\|_1$ norm, and the PREP/SELECT decomposition.

### Citations

- `motlagh2024generalized`: now also used for Thm. 6, which is the central construction in v6. Citation annotation in the Citations section expanded accordingly.
- `berntsonSunderhauf2025complementary`: unchanged. Bibtex `number = {7}` removed per review m3 (CMP 406:161 has no issue number).
- `lowChuang2019qubitization`, `haah2019product`, `dong2021efficient`: remain dropped as in v5.

### Resolved in v6

- **FATAL F1** (shift trick implements $W^{d_{\text{JA}}}\ee^{-\ii\delta B}$, not $\ee^{-\ii\delta B}$): resolved by switching to Thm. 6 Laurent route.
- **MAJOR M1** (cost accounting off by 2 under shift trick): resolved by stating $d_{\text{JA}}$ + $d_{\text{JA}}$ split matching Cor. 8 exactly.
- **MAJOR M2** (trailing $\ee^{\ii\theta}$ factor in $P_2$): resolved by switching to $L_1$ Laurent polynomial, which has no such factor.
- **MAJOR M3** (Route 1 two-call LCU has its own $P(z)=P(1/z)$ obstruction): acknowledged in Writing Notes §"The shift-by-$z^{d_{\text{JA}}}$ step (resolved and closed)" — Route 1 is not a viable alternative.
- **MINOR m1** (Bessel tail $+1$): applied in target-polynomial paragraph.
- **MINOR m2** (Jacobi–Anger attribution): applied — `$J_{-k}$` invocation dropped from the sign derivation; `$J_n(-z) = (-1)^n J_n(z)$` used directly.
- **MINOR m3** (bibtex `number = {7}`): removed.
- **MINOR m4** (memo reference): Suggestion 2 now flags the memo should be updated to reflect v6's methods-side punch list (not v4/v5's).
- **MINOR m5** (caption length): compressed to two sentences per review suggestion.
- **MINOR m6** ($B_a$ vs $\tilde B_a$): added as a bullet in Suggestion 2.
- **Review item 10** (notation consistency — generic $d$ vs $d_{\text{JA}}$): one sentence added to the GQSP layer / polynomial class paragraph tying the generic $d$ to $d = 2d_{\text{JA}}$ explicitly.
