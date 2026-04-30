# Quantum Signal Processing (v7 — generic, Laurent via controlled-$W^\dagger$)

> **Insertion target:** `supplementary-informations/1_preliminaries.tex:1068`, body of `\subsection{Quantum Signal Processing.} \label{sec:prelim-QSP}` — between the `\subsection{...}\label{...}` line (line 1068) and `\begin{figure}[ht]` (line 1071). The figure (Fig. `circ:gqsp`, lines 1071–1113) follows immediately.
>
> **Change from v6 (scope).** Rewritten generically for a Hermitian operator $A$ with block encoding $U_A$, subnormalisation $\alpha$, step size $\delta$, and truncation index $d$. All algorithm-specific objects ($B$, $B_a$, $b_\pm$, $\alpha = \|b_-\|_1\|b_+\|_1$, `alg:coh`, the Metropolis/Gaussian kernels, the Chapter 5 regime) are removed from the body and collected in the **TODOs** block as substitutions for the later §\textsc{CoherentStep} subsection. The v6 symbol $d_{\text{JA}}$ is renamed $d$ throughout; the abstract GQSP polynomial degree is referenced only through $\deg P$ so the two roles no longer collide. All historical commentary about the MW §III framing and earlier draft bugs is removed; the Laurent route (Thm. 6) is presented directly.
>
> **Existing labels referenced:** `sec:prelim-LCU`, `circ:gqsp`.
> **New labels introduced:** `eq:qsp-walk`, `eq:qsp-poly`, `eq:jacobi-anger`, `eq:qsp-degree`.
> **Citation keys:** `motlagh2024generalized` (existing), `berntsonSunderhauf2025complementary` (new).

---

## Body

The LCU circuit of (§`sec:prelim-LCU`) delivers a block encoding $U_A$ on a $b$-qubit block register of any Hermitian operator $A$ on an $n$-qubit system register,
$$
\bra{0^b}\,U_A\,\ket{0^b} \;=\; A/\alpha, \qquad \|A/\alpha\| \leq 1,
$$
equivalently $A = \sum_j \alpha_j U_j$ with $\alpha = \sum_j \alpha_j$, preparation unitary $\text{PREP}\ket{0} = \sum_j \sqrt{\alpha_j/\alpha}\,\ket{j}$, selection unitary $\text{SELECT} = \sum_j \ket{j}\!\bra{j}\otimes U_j$, and $U_A = (\text{PREP}^\dagger\otimes I)\,\text{SELECT}\,(\text{PREP}\otimes I)$. *Generalized quantum signal processing* (GQSP) [Motlagh & Wiebe 2024] [CITE: motlagh2024generalized] converts $U_A$ into any polynomial transformation of $A$'s spectrum by interleaving controlled calls to an associated *walk operator* with single-qubit rotations on a QSP ancilla, post-selecting both ancillas on $\ket{0}$. Choosing the polynomial as a good approximation to $\ee^{-\ii\delta\alpha\cos\theta}$ on the unit circle realises the Hamiltonian-simulation unitary $\ee^{-\ii\delta A}$ up to a QSP truncation error — the application that will be invoked for every coherent step downstream. Figure `circ:gqsp` depicts one GQSP layer.

*Walk operator.* Following [Motlagh & Wiebe 2024, Cor. 8] [CITE: motlagh2024generalized], define
$$
W \;:=\; \bigl(2\,\text{PREP}\ket{0}\!\bra{0}\,\text{PREP}^\dagger \,-\, I\bigr)\,\text{SELECT},
$$ <!-- \label{eq:qsp-walk} -->
a bare SELECT followed by a reflection about the prepared state $\text{PREP}\ket{0}$ (note the PREP-sandwich: a full $U_A$ is *not* applied — only one SELECT and one reflection). For each eigenvector $\ket{\lambda}$ of $A$ with $A\ket{\lambda} = \lambda\ket{\lambda}$, the two-dimensional subspace containing $\text{PREP}\ket{0}\ket{\lambda}$ is $W$-invariant, and in a suitable orthonormal basis $W$ acts as the $\mathrm{SU}(2)$ matrix [Motlagh & Wiebe 2024, Eq. (64)] [CITE: motlagh2024generalized]
$$
W_\lambda \;=\; \begin{pmatrix} \lambda/\alpha & \sqrt{1-(\lambda/\alpha)^2} \\[2pt] -\sqrt{1-(\lambda/\alpha)^2} & \lambda/\alpha \end{pmatrix},
$$
with trace $2\lambda/\alpha$ and determinant $+1$. Its eigenvalues form the conjugate pair $\ee^{\pm\ii\theta}$ with
$$
\cos\theta \;=\; \lambda/\alpha, \qquad \lambda \in [-\alpha,\alpha], \qquad \theta \in [0,\pi].
$$
Any spectral-transformation task on the $\lambda$-eigenspace of $A$ thus reduces to approximating the target function at $z = \ee^{\ii\theta}$ on the unit circle $\mathbb{T} = \{z\in\mathbb{C} : |z|=1\}$.

*GQSP layer and Laurent polynomial class.* A GQSP circuit (Fig. `circ:gqsp`) alternates $(\deg P)$ (open-)controlled calls to $W$ or $W^\dagger$ with $(\deg P) + 1$ single-qubit rotations on the QSP ancilla $\ket{0}_\text{QSP}$: a three-parameter opening rotation $R(\theta_0,\phi_0,\lambda_0)$ (the $\lambda_0\in\mathbb{R}$ is a global-phase parameter, distinct from the spectral variable $\lambda$) and $\deg P$ further two-parameter rotations $R(\theta_j,\phi_j,0)$, $j=1,\dots,\deg P$. [Motlagh & Wiebe 2024, Thm. 6 (Eqs. 45–53)] [CITE: motlagh2024generalized] characterises the top-left block of the post-selected circuit on each qubitization subspace as
$$
P'(W_\lambda) \;=\; W_\lambda^{-k}\,P(W_\lambda), \qquad P \in \mathbb{C}[z], \quad 0 \leq k \leq \deg P, \quad \|P\|_{\infty,\mathbb{T}} := \max_{z\in\mathbb{T}}|P(z)| \leq 1,
$$ <!-- \label{eq:qsp-poly} -->
where $k$ is the number of controlled-$W^\dagger$ slots (and $\deg P - k$ is the number of controlled-$W$ slots) — i.e., any Laurent polynomial $P' \in \mathbb{C}[z,z^{-1}]$ whose ordinary factor $P$ has unit sup-norm on $\mathbb{T}$ is accessible. The rotation angles $\{(\theta_0,\phi_0,\lambda_0)\}\cup\{(\theta_j,\phi_j)\}_{j=1}^{\deg P}$ depend on $P$ only; the slot assignment (which $k$ slots are flipped) encodes the shift $z^{-k}$ independently. Controlled-$W^\dagger$ and controlled-$W$ carry equal hardware cost: both are open-controlled calls to the same block-encoding oracle, so the Thm. 6 substitution incurs no gate premium.

*Target polynomial via Jacobi–Anger.* The Hamiltonian-simulation target $\ee^{-\ii\delta\alpha\cos\theta}$ is a Laurent function of $\ee^{\ii\theta}$ not expressible as a non-constant ordinary polynomial: an ordinary $P\in\mathbb{C}[z]$ satisfying $P(\ee^{\ii\theta}) = P(\ee^{-\ii\theta})$ on $\mathbb{T}$ has $P(z) = P(1/z)$ on $\mathbb{T}$, and matching Fourier coefficients forces $P$ constant. The Thm. 6 Laurent class is exactly what is needed. The bilateral Jacobi–Anger expansion [Motlagh & Wiebe 2024, Eq. (62)] [CITE: motlagh2024generalized]
$$
\ee^{-\ii\delta\alpha\cos\theta} \;=\; \sum_{k=-\infty}^{\infty} (-\ii)^k\,J_k(\delta\alpha)\,\ee^{\ii k\theta}, \qquad \theta \in \mathbb{R},
$$ <!-- \label{eq:jacobi-anger} -->
(sign pattern from $t\to -\delta\alpha$ in $\ee^{\ii t\cos\theta} = \sum_k \ii^k J_k(t)\,\ee^{\ii k\theta}$ together with $J_n(-z) = (-1)^n J_n(z)$) truncated at $|k|\leq d$ yields the Laurent polynomial
$$
L_d(\ee^{\ii\theta}) \;:=\; \sum_{k=-d}^{d} (-\ii)^k\,J_k(\delta\alpha)\,\ee^{\ii k\theta} \;\in\; \mathbb{C}[z,z^{-1}],
$$
of bilateral degree $d$, with Bessel tail bound
$$
\bigl\|\ee^{-\ii\delta\alpha\cos\theta} - L_d\bigr\|_{\infty,\mathbb{T}} \;\leq\; 2\!\!\sum_{|k|>d}\!|J_k(\delta\alpha)| \;\leq\; \mathcal{O}\!\left(\!\bigl(\ee\delta\alpha/(2(d+1))\bigr)^{d+1}\right)
$$
[Motlagh & Wiebe 2024, Eq. (63)] [CITE: motlagh2024generalized]. To realise $L_d$ via Thm. 6, introduce the *underlying* ordinary polynomial
$$
P(z) \;:=\; z^{d}\,L_{d}(z) \;=\; \sum_{m=0}^{2d} c_m\,z^m, \qquad c_m \;=\; (-\ii)^{m-d}\,J_{m-d}(\delta\alpha),
$$
of degree $2d$ — a bookkeeping object on which the angle-extraction (below) runs, *not* an operator that is implemented. Since $|z^d|=1$ on $\mathbb{T}$, $\|P\|_{\infty,\mathbb{T}} = \|L_d\|_{\infty,\mathbb{T}} \leq 1 + \mathcal{O}((\ee\delta\alpha/(2(d+1)))^{d+1})$, and a rescaling $(1-\varepsilon_\text{QSP})P$ enforces (Eq. `eq:qsp-poly`) strictly. Applying Thm. 6 with $\deg P = 2d$ and $k = d$ — i.e. $d$ of the $2d$ signal slots deployed as controlled-$W^\dagger$, the remaining $d$ as controlled-$W$ — produces
$$
P'(z) \;=\; z^{-d}\,P(z) \;=\; L_d(z),
$$
so the post-selected top-left block directly realises the Jacobi–Anger Laurent truncation, giving $L_d(W_\lambda) = \ee^{-\ii\delta\alpha\cos\theta_\lambda}\,I_{2\times 2} + \mathcal{O}(\varepsilon_\text{QSP}) = \ee^{-\ii\delta\lambda}\,I_{2\times 2} + \mathcal{O}(\varepsilon_\text{QSP})$ on each qubitization subspace, and hence $\ee^{-\ii\delta A} + \mathcal{O}(\varepsilon_\text{QSP})$ on the system register.

*Degree, cost, small-$\delta\alpha$ default.* The truncation index $d$ that suffices for $\varepsilon_\text{QSP}$-accuracy on $\mathbb{T}$ is [Motlagh & Wiebe 2024, Thm. 7] [CITE: motlagh2024generalized]
$$
d \;=\; \mathcal{O}\!\left(\delta\alpha \;+\; \frac{\log(1/\varepsilon_\text{QSP})}{\log\log(1/\varepsilon_\text{QSP})}\right).
$$ <!-- \label{eq:qsp-degree} -->
The resulting GQSP circuit uses $2d$ controlled signal slots in total — half of them ($d$ slots) deployed as controlled-$W$ and the other half as controlled-$W^\dagger$ per Thm. 6 with $k = d$ — plus $2d+1$ single-qubit rotations on the QSP ancilla, matching the cost bound of [Motlagh & Wiebe 2024, Cor. 8] [CITE: motlagh2024generalized]. Explicitly: $d=1$ needs 1 controlled-$W$ and 1 controlled-$W^\dagger$; $d=2$ needs 2 and 2; $d=3$ needs 3 and 3; etc. In the small-$\delta\alpha$ regime $\delta\alpha \lesssim 1$, already $d = 1$ suffices:
$$
L_1(\ee^{\ii\theta}) \;=\; J_0(\delta\alpha) \,-\, 2\ii\,J_1(\delta\alpha)\,\cos\theta,
$$
realised by one controlled-$W$ call, one controlled-$W^\dagger$ call, and three rotation triples $(R_0, R_1, R_2)$. Under $J_0(z) = 1 + \mathcal{O}(z^2)$ and $2J_1(z) = z + \mathcal{O}(z^3)$,
$$
L_1(\ee^{\ii\theta}) \;=\; 1 \,-\, \ii\,\delta\alpha\,\cos\theta \,+\, \mathcal{O}((\delta\alpha)^2) \;=\; \ee^{-\ii\delta\alpha\cos\theta} \,+\, \mathcal{O}((\delta\alpha)^2),
$$
with no residual factor. Larger $\alpha$ or a larger $\delta$ at fixed total accuracy push $d \geq 2$ via (Eq. `eq:qsp-degree`); the cost then grows linearly in $\delta\alpha$.

*Angle finding.* Given $\{c_m\}_{m=0}^{2d}$, the angle vector $\{(\theta_0,\phi_0,\lambda_0)\}\cup\{(\theta_j,\phi_j)\}_{j=1}^{2d}$ is extracted in two classical steps. First, a *complementary polynomial* $Q\in\mathbb{C}[z]$ with $\deg Q = 2d$ and
$$
|P(z)|^2 + |Q(z)|^2 \;=\; 1 \qquad \text{on } \mathbb{T}
$$
is constructed; existence is [Motlagh & Wiebe 2024, Thm. 4] [CITE: motlagh2024generalized]. We use the contour-integral / FFT construction of [Berntson & Sünderhauf 2025] [CITE: berntsonSunderhauf2025complementary], which evaluates $Q$ in classical time $\tilde{\mathcal{O}}(d)$ with explicit a priori error bounds. Second, the recursive extraction of [Motlagh & Wiebe 2024, §IV, Alg. 1] [CITE: motlagh2024generalized] reads off the $2d+1$ rotation triples from $(P,Q)$ in $\mathcal{O}(d^2)$ arithmetic operations. Both steps are performed once, offline, for each value of $\delta\alpha$; the same angle vector is reused at every quantum invocation of the layer. The $2d$ signal slots are deployed as $d$ controlled-$W$ followed by $d$ controlled-$W^\dagger$ per [Motlagh & Wiebe 2024, Eq. (46)] [CITE: motlagh2024generalized].

---

## TODOs — specialisations for the later §\textsc{CoherentStep} subsection (`2_methods.tex:1606–1740`)

The coherent-step subsection consumes the preliminaries above via the substitutions listed here. Each item replaces a generic object of the prelim by an algorithm-specific one; no new GQSP facts are introduced.

1. **Hermitian operator.** $A \to B_a$, the per-jump coherent term of (Eq. `eq:B-bohr-domain`). The block encoding $U_A$ is the nested-LCU circuit $U_{B_a}$ of Fig. `circ:U_B_a`, constructed in Step 2 of `alg:coh` on outer/inner time registers $T_\mp$.

2. **Subnormalisation.** $\alpha \to \|b_-\|_1\,\|b_+\|_1$ per (Eq. `eq:B-block-encoding`); its scaling — bounded for the Gaussian weight, logarithmic for the Metropolis weight $b_+^{(s,\eta)}$ (Eq. `eq:b_plus-s-eta`) — is what determines whether $d=1$ is sufficient.

3. **Step size.** $\delta \to \delta$ of `alg:coh` (the coherent step of the splitting; used to match the weak-measurement dissipative error).

4. **Default truncation index.** $d = 1$ in the Chapter 5 regime ($\delta\alpha \lesssim 1$); cite $L_1$ (body) and the three rotation triples. State explicitly that (a) low-temperature Metropolis (large $\alpha$) and (b) high-precision requirements (small $\delta$) force $d \geq 2$, with cost growing linearly in $\delta\alpha$.

5. **Walk $W$.** The walk of (Eq. `eq:qsp-walk`) is built from $U_{B_a}$'s combined PREP$\,=\,$PREP$'_-\otimes\text{PREP}_+$ (and the corresponding SELECT) — the reflection is about the nested prepared state, *not* the raw-ancilla reflection $R_b = 2\ket{0^b}\!\bra{0^b}-I$ currently shown in Fig. `circ:gqsp` (see Suggestion 1 for the required figure redraw).

6. **Cost accounting.** The $2d$ controlled signal slots become $d$ calls to $U_{B_a}$ and $d$ calls to $U_{B_a}^\dagger$, each one an outer+inner PREP + controlled Strang $S_2$ + outer+inner $\text{PREP}^\dagger$ (or its reverse). Plug into `alg:coh`'s coherent-step total; double the block-encoding-oracle count relative to a single $U_{B_a}$ application.

7. **$B_a$ vs $\tilde B_a$.** `alg:coh`'s coherent step produces $\ee^{-\ii\delta\tilde B_a} = \ee^{-\ii\delta B_a} + \mathcal{O}(\text{Trotter})$ with $\tilde B_a$ the Trotterized version. A single sentence in §\textsc{CoherentStep} should fix the convention: either use $\tilde B_a$ uniformly in the block encoding and the QSP output, or flag the Trotter approximation explicitly.

8. **Pseudocode update in `alg:coh` (line 1609).** "$P(\ee^{\ii\theta}) \approx \ee^{-\ii\delta\alpha\sin\theta}$" $\to$ "$L_d(\ee^{\ii\theta}) \approx \ee^{-\ii\delta\alpha\cos\theta}$ (Laurent polynomial, bilateral degree $d$, MW Thm. 6)". Update `\Comment{Degree $d=1$ calls to $U_{B_a}$}` at line 1611 to `\Comment{$d=1$: one call to $U_{B_a}$ + one call to $U_{B_a}^\dagger$}`.

9. **Trig branch (line 1731).** "its eigenvalues come in conjugate pairs $\ee^{\pm\ii\theta}$ where $\sin\theta = \lambda$ for each eigenvalue $\lambda$ of $B_a/\alpha$" $\to$ "$\ldots$ where $\cos\theta = \lambda/\alpha$ for each eigenvalue $\lambda$ of $B_a$". Also: update "$P(\ee^{\ii\theta}) \approx \ee^{-\ii\delta\alpha\sin\theta}$" $\to$ "$L_d(\ee^{\ii\theta}) \approx \ee^{-\ii\delta\alpha\cos\theta}$". Introduce $W$ by reference to (Eq. `eq:qsp-walk`).

10. **Degree bound paragraph (lines 1733–1737).** The bound's "$d$" matches this prelim exactly; add a parenthetical "(the circuit uses $2d$ signal-operator slots: $d$ controlled-$U_{B_a}$ + $d$ controlled-$U_{B_a}^\dagger$ per MW Thm. 6)". Restate the $d=1$ default in Laurent-polynomial language: "Jacobi–Anger truncation $d=1$ is sufficient, giving a Laurent polynomial of bilateral degree 1 and one controlled-$U_{B_a}$ + one controlled-$U_{B_a}^\dagger$ per coherent step."

---

## Citations

- **`motlagh2024generalized`** (already in `references.bib`). Thm. 6 (Eqs. 45–53) gives the Laurent class, realised by flipping $k$ of the $\deg P$ signal slots to controlled-$W^\dagger$. Thm. 7 (Eq. 62) gives the Jacobi–Anger degree bound. Cor. 8 (Eqs. 64–66) gives the PREP-sandwiched walk and its $2{\times}2$ block action. Alg. 1 (§IV, p.9) gives the $\mathcal{O}(d^2)$ recursive angle extraction. Thm. 4 (p.6) guarantees the complementary polynomial $Q$. MW Eq. 46 fixes the canonical slot assignment (last $k$).
- **`berntsonSunderhauf2025complementary`** (new). Berntson, B. K. & Sünderhauf, C., *"Complementary Polynomials in Quantum Signal Processing"*, **Commun. Math. Phys. 406, 161 (2025)**, arXiv:2406.04246. Contour-integral / FFT construction of the complementary $Q$ in classical time $\tilde{\mathcal{O}}(\deg P)$ with explicit a priori error bounds; demonstrated numerically to $\deg P \sim 10^7$.

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

### Notation

- $A$: generic Hermitian operator on the $n$-qubit system register; $\alpha$: LCU subnormalisation ($\|A/\alpha\|\leq 1$); $\lambda\in[-\alpha,\alpha]$: eigenvalue of $A$ in physical units (not $A/\alpha$); $\delta$: step size; $d$: Jacobi–Anger truncation index (= bilateral degree of $L_d$); $\deg P = 2d$: degree of the underlying ordinary polynomial; $k = d$: number of flipped (controlled-$W^\dagger$) slots in the Hamiltonian-simulation specialisation.
- $\mathbb{T} = \{z\in\mathbb{C} : |z|=1\}$.
- $\text{PREP}\ket{0} = \sum_j\sqrt{\alpha_j/\alpha}\ket{j}$; the reflection about $\text{PREP}\ket{0}$ is $2\,\text{PREP}\ket{0}\bra{0}\text{PREP}^\dagger - I$, distinct from $R_b = 2\ket{0^b}\bra{0^b} - I$.
- The $\lambda_0\in\mathbb{R}$ of the opening rotation $R(\theta_0,\phi_0,\lambda_0)$ is a global-phase parameter and is **not** the spectral variable $\lambda$.

### Derivation of $L_1$

From (Eq. `eq:jacobi-anger`) truncated at $|k|\leq 1$: $L_1(\ee^{\ii\theta}) = (-\ii)^{-1}J_{-1}\ee^{-\ii\theta} + J_0 + (-\ii)^1 J_1\ee^{\ii\theta}$. With $J_{-1}(z) = -J_1(z)$ and $(-\ii)^{-1} = \ii$:
$$
L_1 \;=\; \ii(-J_1)\ee^{-\ii\theta} + J_0 - \ii J_1\ee^{\ii\theta} \;=\; -\ii J_1(\ee^{-\ii\theta}+\ee^{\ii\theta}) + J_0 \;=\; J_0 - 2\ii J_1\cos\theta.
$$
Under $J_0(z) = 1 + \mathcal{O}(z^2)$ and $2J_1(z) = z + \mathcal{O}(z^3)$, $L_1 = 1 - \ii\delta\alpha\cos\theta + \mathcal{O}((\delta\alpha)^2)$. The underlying ordinary polynomial is $P(z) = zL_1(z) = -\ii J_1\,z^2 + J_0\,z - \ii J_1$, of degree $2$; Alg. 1 run on $P$ yields three rotation triples, and MW Eq. 46 assigns slot 1 as controlled-$W$, slot 2 as controlled-$W^\dagger$ (last-$k$ canonical pattern).

### Polynomial-degree / slot bookkeeping

Three numbers are easy to confuse; the "factor of 2" runs as follows:

| quantity | value |
|----------|-------|
| bilateral degree of the Laurent $L_d$ actually implemented | $d$ |
| degree of the underlying ordinary $P = z^d L_d$ (on which Alg. 1 runs) | $2d$ |
| number of controlled signal slots in the GQSP layer | $2d$ |
| number of controlled-$W$ calls | $d$ |
| number of controlled-$W^\dagger$ calls | $d$ |
| number of rotation triples | $2d + 1$ |
| $\deg Q$ of the complementary polynomial | $2d$ |

"$\mathcal{O}(\alpha\delta + \log/\log\log)$ applications of $W$" in MW Cor. 8 reads as "applications of $W$ **or** $W^\dagger$" under Thm. 6.

### Sign derivation in (Eq. `eq:jacobi-anger`)

From the standard identity $\ee^{\ii t\cos\theta} = \sum_k \ii^k J_k(t)\ee^{\ii k\theta}$, substitute $t \to -\delta\alpha$ and use $J_n(-z) = (-1)^n J_n(z)$ (DLMF 10.2.2, 10.12.1): the factor $\ii^n(-1)^n = (-\ii)^n$ appears directly. No separate invocation of $J_{-k}$ is needed at the level of the generic expansion (the $J_{-1}$ step in the derivation of $L_1$ above is a coefficient-level simplification at $d=1$, not a general expansion fact).

---

## Suggestions

1. **[FIX-NEXT] Redraw `circ:gqsp` (`1_preliminaries.tex:1071–1113`).** The current figure draws $W = R_b\,U_B$ with $U_B$ as a bundled block and $R_b$ as an external $\ket{0^b}$-reflection; this does not match the PREP-sandwich $W$ of (Eq. `eq:qsp-walk`), nor does it expose the Thm. 6 $W/W^\dagger$ split.

   *Option A (recommended, matches the body).* Expand the $U_A$ box into PREP, SELECT, PREP$^\dagger$; drop the outer PREP$^\dagger$ / PREP that would make a full $U_A$; per GQSP layer, draw: (i) one open-controlled SELECT (or SELECT$^\dagger$ in a flipped slot); (ii) the PREP-sandwiched reflection $2\,\text{PREP}\ket{0}\bra{0}\text{PREP}^\dagger - I$; (iii) the $R(\theta_j,\phi_j,0)$ rotation on the QSP ancilla. In the Hamiltonian-simulation specialisation, draw the last $d$ of the $2d$ slots with a dagger on SELECT (controlled-$W^\dagger$) per MW Eq. 46.

   *Option B (bundled, acceptable).* Keep the $U_A$ box; annotate signal boxes as "controlled-$W$ / controlled-$W^\dagger$"; indicate the $d$–$d$ split in the caption.

   **Caption rewrite is required in both options.** Current caption line 1111, *"degree-$d$ Chebyshev polynomial $P(W)\approx\ee^{-\ii\delta B}$"*, is wrong on two counts (GQSP produces a Laurent polynomial, not a Chebyshev polynomial; no $B$ or $\delta$ enters the generic prelim). Suggested:
   > *"GQSP layer of depth $\deg P$. With $k$ of the $\deg P$ controlled signal slots deployed as controlled-$W^\dagger$ (Thm. 6 of~\cite{motlagh2024generalized}), the post-selected top-left block on each qubitization subspace realises the Laurent polynomial $P'(W) = W^{-k}P(W)$. In the Hamiltonian-simulation application (§\textsc{CoherentStep}), $\deg P = 2d$ and $k = d$ reproduce $\ee^{-\ii\delta A}$ up to Jacobi–Anger truncation error."*

2. **Section title.** `sec:prelim-QSP` and the heading "Quantum Signal Processing" are an umbrella; this subsection is strictly GQSP. Consider retitling *"Generalized Quantum Signal Processing"* with label `sec:prelim-GQSP`, updating cross-references in `2_methods.tex`. Low priority.

3. **Explicit $d=1$ angle formulas.** The seven-parameter rotation vector for the default circuit admits closed-form expressions in $J_0(\delta\alpha), J_1(\delta\alpha)$ via Alg. 1 run on $P(z) = -\ii J_1 z^2 + J_0 z - \ii J_1$ and its $Q$. Listing these would make the "default in practice" claim self-verifying; requires separate numerical verification; deferred.
