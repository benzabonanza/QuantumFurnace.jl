# Draft Review: QSP Subsection v5 (GQSP, ordinary-polynomial route)

**Reviewed**: `drafts/qsp-subsection-v5.md`
**Primary source**: `supplementary-informations/Motlagh and Wiebe - 2024 - Generalized Quantum Signal Processing.pdf` (Thm. 3, Cor. 5, Thm. 6, Thm. 7, Cor. 8, Eqs. 62–68; §IV Alg. 1; §V.A)
**Secondary source**: `supplementary-informations/Berntson and Sünderhauf - 2025 - Complementary Polynomials in Quantum Signal Processing.pdf` (Thm. 2, Cor. 1.2, Alg. 1)
**Date**: 2026-04-23
**Overall assessment**: v5 is a substantial and mostly correct rewrite. All the routine mathematical checks pass (Jacobi–Anger signs, the $P_2$ derivation, the coefficient formula $c_m = (-\ii)^{m-d_{\text{JA}}} J_{m-d_{\text{JA}}}(\delta\alpha)$, the $W_\lambda$ matrix entries, the trace/det/eigenvalue statements, the Berntson–Sünderhauf runtime claim). The draft's own flagged `[CHECK]` is the single real mathematical issue in the body, and the review below confirms it is a genuine problem with the single-circuit construction as presented. Not ready to insert as-is: **the `[CHECK]` tension must be resolved before commitment**, and the figure `circ:gqsp` redraw (Suggestion 1 in the draft) is required for self-consistency. With those two items addressed, v5 can land.

---

## FATAL

### F1: The shift-trick single-circuit construction does not implement $\ee^{-\ii\delta B}$ — it implements $W^{d_{\text{JA}}}\,\ee^{-\ii\delta B}$ on the qubitization subspace.

- **Location**: Body paragraph *"Small-$\delta\alpha$ default"* (line 73 of the draft), and the closing claim at line 81 — *"The unitary $\ee^{-\ii\delta B}$ produced by (Fig. `circ:gqsp`) closes the block-encoding pipeline..."*.
- **Problem**: The draft's body sentence
  > "reproducing $\ee^{-\ii\delta\alpha\cos\theta}$ to $\mathcal{O}((\delta\alpha)^2)$ up to the trailing unimodular $\ee^{\ii\theta}$ factor that is absorbed into the qubitization bookkeeping (see Writing Notes on the shift-by-$z^{d_{\text{JA}}}$ identity)"

  and the closing sentence
  > "The unitary $\ee^{-\ii\delta B}$ produced by (Fig. `circ:gqsp`) closes the block-encoding pipeline..."

  under-state what the Writing Notes `[CHECK]` correctly flags. On the qubitization subspace for eigenvalue $\lambda$, the ordinary polynomial $P_{2d_{\text{JA}}}$ built from the shift trick implements
  $$P_{2d_{\text{JA}}}(W)\ket{\phi^\pm_\lambda} \;=\; \ee^{\pm\ii d_{\text{JA}}\theta_\lambda}\,L_{d_{\text{JA}}}(\ee^{\pm\ii\theta_\lambda})\,\ket{\phi^\pm_\lambda} \;\approx\; \ee^{\pm\ii d_{\text{JA}}\theta_\lambda}\,\ee^{-\ii\lambda\delta}\,\ket{\phi^\pm_\lambda}$$
  — i.e. the map $W_\lambda \mapsto \operatorname{diag}(\ee^{\ii d_{\text{JA}}\theta_\lambda - \ii\lambda\delta},\, \ee^{-\ii d_{\text{JA}}\theta_\lambda - \ii\lambda\delta}) = W_\lambda^{d_{\text{JA}}}\cdot\operatorname{diag}(\ee^{-\ii\lambda\delta}, \ee^{-\ii\lambda\delta}) = W_\lambda^{d_{\text{JA}}}\,F_\lambda$ in MW's notation. This is **not** MW Eq. (66)'s target $F_\lambda = \operatorname{diag}(\ee^{-\ii\lambda\delta}, \ee^{-\ii\lambda\delta})$. The trailing factor $W_\lambda^{d_{\text{JA}}}$ is **$\lambda$-dependent** — in particular it depends on $\theta_\lambda = \arccos(\lambda/\alpha)$ — so it is not a global phase that can be "absorbed into bookkeeping". It couples to the spectrum and, post-projection with $\text{PREP}|0\rangle$, leaves a residual unitary $(W|_\text{subspace})^{d_{\text{JA}}}$ on top of $\ee^{-\ii\delta B}$.
- **Independent verification that no ordinary polynomial can realise MW Eq. (66) directly.** An ordinary $P(z) = \sum_{m=0}^{2d} c_m z^m \in \mathbb{C}[z]$ acts on the 2D qubitization subspace as
  $$P(W_\lambda) \;=\; \operatorname{diag}(P(\ee^{\ii\theta_\lambda}), P(\ee^{-\ii\theta_\lambda}))$$
  in the $\{|\phi^+_\lambda\rangle,|\phi^-_\lambda\rangle\}$ eigenbasis. MW Eq. (66) asks for the two diagonal entries to be $\ee^{-\ii\lambda\delta}$, i.e. **equal** — which forces $P(z) = P(1/z)$ on $\mathbb{T}$. Writing $P(z) = \sum_{m=0}^{2d} c_m z^m$ and $P(1/z) = \sum_{m=0}^{2d} c_m z^{-m}$ and matching Fourier coefficients on $\mathbb{T}$ forces $c_m = 0$ for all $m \geq 1$ — i.e. $P$ constant. Hence MW Cor. 8 cannot be realised within pure Cor. 5 / ordinary $\mathbb{C}[z]$ by a non-constant polynomial. Either (i) one accepts the $\lambda$-dependent trailing $W^{d_{\text{JA}}}$ factor (what the shift trick actually gives); or (ii) one invokes MW Thm. 6 (controlled $U^\dagger$, Laurent), contradicting MW §III's pitch; or (iii) one uses LCU with two GQSP calls (one for an even function of $\theta$, one for an odd, but see the [CHECK] adjudication below — that route has its own subtlety).
- **Why it matters**: as written, the draft claims a Hamiltonian-simulation routine that does not exist — the subalgorithm `CoherentStep` would not produce $\ee^{-\ii\delta B_a}$ as stated in `2_methods.tex:1705`, it would produce $W^{d_{\text{JA}}}\ee^{-\ii\delta B_a}$ (on the post-selected qubitization subspace). Any downstream Chapter-5 claim that the coherent step realises $\mathcal{L}_{a,\text{coh}}(\cdot) = -\ii[B_a,\cdot]$ is then off by a $W^{d_{\text{JA}}}$ conjugation that doesn't commute with $B_a$.
- **Suggested fix (minimal)**: the cleanest fix is to commit to an explicit route and rewrite the *Target polynomial via Jacobi–Anger* paragraph and the closing sentence accordingly. See the **[CHECK] adjudication** section below for the recommendation; in short, **Route 2 (reinstate MW Thm. 6, Laurent, controlled-$U^\dagger$, soften §III pitch)** is the right commitment for the thesis, and the textual changes are local. Draft patch:
  - In the *Target polynomial* paragraph, replace
    > "This Laurent polynomial is not directly in GQSP's polynomial class. To land in $\mathbb{C}[z]$, we multiply by the unimodular factor $\ee^{\ii d_{\text{JA}}\theta} = z^{d_{\text{JA}}}$, $P_{2d_{\text{JA}}}(z) \;:=\; z^{d_{\text{JA}}}\, L_{d_{\text{JA}}}(z) \;=\; \sum_{m=0}^{2d_{\text{JA}}} c_m\,z^m$..."

    with
    > "MW's Thm. 6 [Motlagh & Wiebe 2024, Thm. 6] [CITE: motlagh2024generalized] extends the Cor. 5 class to Laurent polynomials in $z = \ee^{\ii\theta}$ — i.e. polynomials $P'(z) = z^{-k} P(z)$ with $P \in \mathbb{C}[z]$, $\deg P \leq d$, $0\leq k\leq d$ — by replacing $k$ of the controlled-$W$ calls in Fig. `circ:gqsp` by controlled-$W^\dagger$ calls. Choosing $d = 2d_{\text{JA}}$, $k = d_{\text{JA}}$, and $P(z) = z^{d_{\text{JA}}}L_{d_{\text{JA}}}(z) \in \mathbb{C}[z]$ gives $P'(z) = L_{d_{\text{JA}}}(z)$ — the Jacobi–Anger Laurent truncation directly, realised via $d_{\text{JA}}$ calls to $W$ and $d_{\text{JA}}$ calls to $W^\dagger$."
  - In the *Small-$\delta\alpha$ default* paragraph, rewrite the $P_2$ derivation around the direct Laurent $L_1(\ee^{\ii\theta}) = J_0 - 2\ii J_1\cos\theta$ (which genuinely equals $\ee^{-\ii\delta\alpha\cos\theta} + \mathcal{O}((\delta\alpha)^2)$ on $\mathbb{T}$, no trailing $\ee^{\ii\theta}$ factor, but needs one $W$ call and one $W^\dagger$ call per layer under Thm. 6).
  - Add one sentence to the *GQSP layer and polynomial class* paragraph acknowledging that §III's "no Laurent" pitch is about Cor. 5 / Thm. 3 alone, and that Thm. 6 re-introduces $U^\dagger$ specifically for phase-function applications — including Cor. 8's Hamiltonian simulation.

If the user prefers to keep the shift trick, the fix is larger: the body claim "$\ee^{-\ii\delta B}$" must be weakened to "$W^{d_{\text{JA}}}\,\ee^{-\ii\delta B}$" everywhere, and a paragraph must be added explaining how the $W^{d_{\text{JA}}}$ factor is implemented / compensated (it is not, without further machinery, a benign correction).

---

## MAJOR

### M1: Cost accounting is inconsistent with MW Cor. 8 as written.

- **Location**: Paragraph *"Degree bound and cost"*, line 67: *"The GQSP circuit in (Fig. `circ:gqsp`) makes $2d_{\text{JA}}$ controlled calls to $W$ (one per polynomial degree) ... matching MW's Cor. 8 bound."*
- **Problem**: MW Cor. 8's stated cost is "$\mathcal{O}(\alpha t + \log(1/\epsilon)/\log\log(1/\epsilon))$ applications of $W$" (MW pg. 11, end of Cor. 8's proof), with the polynomial degree $d$ from MW Thm. 7 equal to the Jacobi–Anger index $d_{\text{JA}}$ in the draft's notation. MW's cost is therefore exactly $d_{\text{JA}}$ applications of $W$, not $2d_{\text{JA}}$. The factor-of-2 the draft introduces is the **honest** cost of the shift-trick ordinary-polynomial construction: under Cor. 5 with an ordinary degree-$2d_{\text{JA}}$ polynomial, Fig. 2 of MW indeed has $2d_{\text{JA}}$ controlled-$W$ boxes. MW's Cor. 8 cost of $d_{\text{JA}}$ calls is achieved only under Thm. 6 (Laurent), where $d_{\text{JA}}$ of the $2d_{\text{JA}}$ signal-operator calls are promoted to controlled-$U^\dagger$, saving nothing in controlled-operation count but preserving the Laurent degree. In neither accounting is the cost literally "$2d_{\text{JA}}$ calls of $W$ and this matches MW's Cor. 8 bound."
- **Why it matters**: the thesis's Chapter 5 cost tables depend on this factor; if the user carries the `$2d_{\text{JA}}$` count forward for per-step cost analysis, they are over-counting relative to MW Cor. 8's literal statement by a factor of 2 (and hence over-estimating total coherent-step gate count by that factor).
- **Suggested fix**: tie cost accounting to the route commitment from F1.
  - Under **Route 2 (Thm. 6, Laurent)**: $d_{\text{JA}}$ controlled-$W$ calls + $d_{\text{JA}}$ controlled-$W^\dagger$ calls = $2d_{\text{JA}}$ controlled operations in total, split evenly. This literally matches MW's Thm. 7 and Cor. 8, and the "$2d_{\text{JA}}$ matches MW's bound" sentence becomes **correct** if we read "applications of $W$" as "applications of $W$ or its inverse". Rewrite the line as: "The GQSP circuit makes $d_{\text{JA}}$ controlled-$W$ calls and $d_{\text{JA}}$ controlled-$W^\dagger$ calls, $2d_{\text{JA}}$ controlled operations total, matching MW Cor. 8's $\mathcal{O}(\alpha\delta + \log/\log\log)$ bound."
  - Under **Route 1 (two GQSP circuits)**: each circuit contributes $2d_{\text{JA}}$ calls (since each implements a shift-trick ordinary polynomial for one of $\cos$ / $\sin$), so total is $4d_{\text{JA}}$ controlled-$W$ calls. This is double the MW count, which the text should acknowledge explicitly as the cost of staying Laurent-free.

### M2: The "$\approx \ee^{\ii\theta}\,\ee^{-\ii\delta\alpha\cos\theta}$" approximation in the $P_2$ line is being read by the draft as "$\approx \ee^{-\ii\delta\alpha\cos\theta}$ up to bookkeeping" — but the trailing $\ee^{\ii\theta}$ factor is spectrally dependent.

- **Location**: Line 73: *"$P_2(\ee^{\ii\theta}) = \ee^{\ii\theta}(J_0 - 2\ii J_1\cos\theta) \approx \ee^{\ii\theta}\,(1 - \ii\delta\alpha\cos\theta) \approx \ee^{\ii\theta}\,\ee^{-\ii\delta\alpha\cos\theta}$ ... reproducing $\ee^{-\ii\delta\alpha\cos\theta}$ to $\mathcal{O}((\delta\alpha)^2)$ up to the trailing unimodular $\ee^{\ii\theta}$ factor that is absorbed into the qubitization bookkeeping"*.
- **Problem**: Parallel to F1 at the $d_{\text{JA}}=1$ default. The factor $\ee^{\ii\theta}$ is $\ee^{\ii\arccos(\lambda/\alpha)}$ in the $\lambda$-eigenspace — this is $W_\lambda$ itself (acting on $|\phi^+_\lambda\rangle$). It is not a global phase and it is not "absorbed into bookkeeping"; it is a concrete unitary error that rotates the $|\phi^\pm_\lambda\rangle$ eigenvectors by an eigenvalue-dependent phase, which after projection with $\text{PREP}|0\rangle$ gives a nontrivial residual $W$ factor multiplying $\ee^{-\ii\delta B}$.
- **Suggested fix**: under Route 2, replace the $P_2$ derivation by $L_1(\ee^{\ii\theta}) = J_0 - 2\ii J_1\cos\theta$ directly (no shift), implemented by one controlled-$W$ and one controlled-$W^\dagger$ per layer. Under Route 1, state explicitly that the $\ee^{\ii\theta}$ factor survives as a residual unitary that must be cancelled in the LCU combiner.

### M3: The draft's Writing Notes flags Berry-style two-call LCU as an escape hatch, but that escape hatch has its own problem that needs acknowledging.

- **Location**: Writing Notes *"The shift-by-$z^{d_{\text{JA}}}$ step and the Cor. 8 identification [CHECK]"*, closing sentence: *"the most honest version of the Hamiltonian-simulation construction within pure MW Cor. 5 is to apply two GQSP circuits, one for $\cos$ and one for $\sin$, and combine via Berry et al.'s doubling efficiency (LCU with coefficients $\frac{1}{2}(1,\ii)$)."*
- **Problem**: the target decomposition $\ee^{-\ii\delta\alpha\cos\theta} = \cos(\delta\alpha\cos\theta) - \ii\sin(\delta\alpha\cos\theta)$ gives two real targets, $\cos(\delta\alpha\cos\theta)$ and $\sin(\delta\alpha\cos\theta)$, each an **even** function of $\theta$ (i.e. invariant under $\theta \to -\theta$). The symmetry obstruction from F1 applies to each: even ordinary polynomials in $z$ satisfying $P(z) = P(1/z)$ on $\mathbb{T}$ are still forced to be constant. So even this "two-call" LCU, if one demands ordinary polynomials per branch, runs into the same problem. What the draft has in mind — each branch being an "ordinary function of $\cos\theta$" — is true (the Chebyshev expansion of $\cos(\delta\alpha\cos\theta)$ in $\cos\theta$ uses $T_n(\cos\theta) = \cos(n\theta)$) but the corresponding *polynomial* $P \in \mathbb{C}[z]$ with $P(\ee^{\ii\theta}) = \cos(\delta\alpha\cos\theta)$ must have $P(z) = P(1/z)$, which as above forces $P$ constant.
- **Resolution**: the two-call LCU does work, but only after the same Laurent machinery is invoked per branch. It saves a factor of 2 in total cost relative to implementing a single $\ee^{-\ii\delta\alpha\cos\theta}$ as a Laurent polynomial of twice the length, but it is not a genuine escape from Laurent / $U^\dagger$.
- **Suggested fix**: in the [CHECK] adjudication section's v6 recommendation (Suggestion 6), drop Route 1 as a non-solution to the Cor. 5 issue and commit to Route 2 outright. See the dedicated [CHECK] section below.

---

## MINOR

### m1: Bessel tail bound, implicit $+1$ shift.

- **Location**: Line 57: *"$\|\ee^{-\ii\delta\alpha\cos\theta} - L_{d_{\text{JA}}}\|_{\infty,\mathbb{T}} \leq 2\sum_{|k|>d_{\text{JA}}}|J_k(\delta\alpha)| \leq \mathcal{O}((\ee\delta\alpha/(2d_{\text{JA}}))^{d_{\text{JA}}})$"*.
- **Issue**: MW Eq. (63) gives $J_n(t) \in \mathcal{O}((\ee t/(2n))^n)$, so the leading tail term at $|k| = d_{\text{JA}}+1$ is $\mathcal{O}((\ee\delta\alpha/(2(d_{\text{JA}}+1)))^{d_{\text{JA}}+1})$. The draft's bound has $d_{\text{JA}}$ where $d_{\text{JA}}+1$ would be tighter. At big-$\mathcal{O}$ this is fine, but v4 wrote the precise form $(\ee\delta\alpha/(2(d+1)))^{d+1}$ — v5 lost that precision.
- **Fix**: write $\mathcal{O}((\ee\delta\alpha/(2(d_{\text{JA}}+1)))^{d_{\text{JA}}+1})$.

### m2: Sign pattern attribution in Jacobi–Anger derivation (nitpick).

- **Location**: Line 53: *"the sign pattern $(-\ii)^k$ comes from setting $t \to -\delta\alpha$ in MW's identity $\ee^{\ii t\cos\theta} = \sum_k \ii^k J_k(t)\,\ee^{\ii k\theta}$ (and using $J_{-k}(z) = (-1)^k J_k(z)$)."*
- **Issue**: $J_{-k}(z) = (-1)^k J_k(z)$ is needed for the $k < 0$ terms of the *original* identity as it is often stated with $k \geq 0$; MW Eq. (62) already runs $n \in \mathbb{Z}$, so $J_n(-\delta\alpha) = (-1)^n J_n(\delta\alpha)$ directly gives $\sum_n \ii^n (-1)^n J_n = \sum_n (-\ii)^n J_n$ without invoking $J_{-k}$. The $J_{-k}$ identity is used later for the $P_2$ coefficient at $k=-1$, which is a separate step.
- **Fix**: drop "(and using $J_{-k}(z) = (-1)^k J_k(z)$)" here; it is distracting. The substitution $t \to -\delta\alpha$ with MW's bilateral sum and $J_n(-z) = (-1)^n J_n(z)$ is all that is used.

### m3: `berntsonSunderhauf2025complementary` bibtex `number = {7}`.

- **Location**: Citations section, bibtex block.
- **Issue**: CMP 406:161 (2025) uses the article-identifier format: the "161" in the journal string *is* the article number. There is no issue/number 7 — the paper's DOI `10.1007/s00220-025-05302-9` and the CMP front page show only volume 406 and article number 161. Similarly, arXiv 2406.04246 front matter shows "Commun. Math. Phys. (2025) 406:161" without an issue number.
- **Fix**: remove the `number = {7}` line from the bibtex entry (or replace with `articleno = {161}` if the user's bibtex style uses it). Keep `pages = {161}` since that is the article ID in CMP.

### m4: `.claude-memory/todo_methods_gqsp_reflection.md` is already superseded.

- **Location**: Writing Notes line 179: *"The user memory `.claude-memory/todo_methods_gqsp_reflection.md` tracks this reconciliation; once the methods edit lands, that memo can be closed."*
- **Issue**: that memo was opened under v3/v4 which still referenced $R_b U_B$. Under v5 the reconciliation items are different (walk operator rewrite, $\sin \to \cos$ in `2_methods.tex:1731`, degree $d \to d_{\text{JA}}$ rename, cost rewrite). The memo should be updated, not just "closed once the methods edit lands".
- **Fix**: minor note. The user should update the memo to reflect the v5 methods-side punch list in Suggestion 2 of the draft, which is correct and more comprehensive than the memo.

### m5: "Chebyshev polynomial" caption word on line 1111 of `1_preliminaries.tex`.

- **Location**: Suggestion 1 in the draft ("Suggested caption" at line 171).
- **Issue**: the suggested caption rewrite is fine but unusually long for a figure caption — it sets the pipeline, names the walk, names PREP and SELECT, and states what post-selection yields, all in one sentence. This is heavy for a figure caption.
- **Fix**: compress to two sentences. Suggested: *"Degree-$2d_{\text{JA}}$ GQSP realising $\ee^{-\ii\delta B}$ from the walk operator $W$ of (Eq.~\ref{eq:qsp-walk}) (MW Cor. 8). Post-selection on $\ket{0}_\mathrm{QSP}\otimes\text{PREP}\ket{0}$ yields the ordinary polynomial $P_{2d_{\text{JA}}}(W)$ which, via the Jacobi–Anger shift trick, implements $\ee^{-\ii\delta B}$ on the qubitization subspace (see `sec:prelim-QSP`)."*

### m6: "$\tilde{B}_a$ vs $B_a$" in Methods Ch. 5.

- **Location**: `2_methods.tex:1681` shows the quantikz figure producing $\tilde{B}_a\ket{\psi}/\alpha$, while line 1700 (CoherentStep text) says *"the per-jump Hermitian term $B_a$"*. The "tilde" is the Trotterized version of $B_a$, distinct from the exact one from \eqref{eq:B-bohr-domain}. v5's body uniformly uses $B$ (no tilde). Methods-side reconciliation (draft Suggestion 2) should also clarify whether `CoherentStep` operates on $B_a$ or $\tilde{B}_a$ — the block encoding \eqref{eq:B-block-encoding} is of $B_a$, but the quantikz of Fig. `circ:U_B_a` produces $\tilde{B}_a$. Under v5 the QSP preliminaries subsection should remain conventional ($B$, no tilde); the methods chapter will have to pay the tilde price.
- **Fix**: no change needed in v5. Add a one-line note to Suggestion 2 so the user doesn't forget the tilde when reconciling `2_methods.tex`.

### m7: Small typography — subscripts on $d_\text{JA}$.

- **Location**: throughout.
- **Issue**: `d_{\text{JA}}` is the correct form and the draft uses it consistently — no issue, flagging as "nothing to fix" to confirm the `d` / `d_{JA}` rename came through cleanly. The methods-side renaming (Suggestion 2) is where $d$ must become $d_{\text{JA}}$.

---

## [CHECK] adjudication: the shift-by-$z^{d_{\text{JA}}}$ step and MW Cor. 8

**Summary of the tension.** MW §III opens with "*our approach eschews the need to use $U^\dagger$ and avoids the use of Laurent polynomials in the analysis*." MW Cor. 5 delivers ordinary polynomials $P \in \mathbb{C}[z]$ only. MW Cor. 8 then claims to implement $\ee^{-\ii\delta H}$ using only Cor. 5 and the walk $W$. The body's `[CHECK]` flag asks whether Cor. 8 is literally consistent with Cor. 5 (ordinary polynomials, no $U^\dagger$), or whether it smuggles in Laurent / $U^\dagger$ under the hood.

**The draft's reading is correct.** An ordinary polynomial $P \in \mathbb{C}[z]$ of degree $2d$ acts on the 2D qubitization subspace as
$$P(W_\lambda) = \operatorname{diag}\!\bigl(P(\ee^{\ii\theta_\lambda}),\,P(\ee^{-\ii\theta_\lambda})\bigr)$$
in the $\{|\phi^+_\lambda\rangle, |\phi^-_\lambda\rangle\}$ eigenbasis of $W_\lambda$. MW Eq. (66)'s target $F_\lambda = \operatorname{diag}(\ee^{-\ii\lambda\delta},\ee^{-\ii\lambda\delta})$ has the two diagonal entries equal, so the polynomial must satisfy $P(\ee^{\ii\theta}) = P(\ee^{-\ii\theta})$ on $\mathbb{T}$, i.e. $P(z) = P(1/z)$ on $\mathbb{T}$. For $P(z) = \sum_{m=0}^{2d} c_m z^m$ this forces $c_m z^m = c_m z^{-m}$ coefficient-by-coefficient in the Fourier expansion, hence $c_m = 0$ for $m \geq 1$ — $P$ must be constant.

Hence MW Cor. 8 **cannot** be realised within pure Cor. 5 (ordinary polynomials) by a non-constant polynomial. The shift trick $P_{2d_{\text{JA}}}(z) = z^{d_{\text{JA}}} L_{d_{\text{JA}}}(z)$ is a legitimate ordinary polynomial of degree $2d_{\text{JA}}$, but it implements $W^{d_{\text{JA}}} \cdot \ee^{-\ii\delta B}$ on the qubitization subspace, not $\ee^{-\ii\delta B}$. The $W^{d_{\text{JA}}}$ factor is $\lambda$-dependent and does not go away.

MW's Cor. 8 proof, read carefully: MW writes Eq. (65) as "apply Cor. 5 ... for any degree $d$ polynomial $P$" and Eq. (66) as "choose $P(U) = \ee^{-\ii\alpha t H}$". The step between (65) and (66) silently lets $P$ be Laurent — MW does not invoke Cor. 5 on an explicit ordinary $P$, they invoke it on $P(U) = \ee^{-\ii\alpha t H}$, a phase function that is not a polynomial at all and must be approximated by a Laurent polynomial via Thm. 7. So MW Cor. 8 is operationally a Thm. 6 + Thm. 7 construction, not a pure Cor. 5 + Thm. 7 construction. The §III opening sentence applies to the GQSP *framework* (Thm. 3, Cor. 5) and to applications that do not need $U^\dagger$ (Thm. 4 for general $|P|\leq 1$ polynomial synthesis) — but Hamiltonian simulation via Cor. 8 is exactly the canonical place where $U^\dagger$ re-enters via Thm. 6.

**No resolution was missed.** There is no operator identity involving PREP, no basis choice, and no parity constraint on $P$ that collapses the asymmetry $P(\ee^{\ii\theta}) = P(\ee^{-\ii\theta})$ within $\mathbb{C}[z]$. The draft's framing is honest.

---

### Route adjudication

**Route 1 (two-GQSP-call LCU with $\cos + \ii\sin$):** does not rescue ordinary polynomials. Each branch $\cos(\delta\alpha\cos\theta)$ and $\sin(\delta\alpha\cos\theta)$ is even in $\theta$, and again requires $P(z) = P(1/z)$ — so any ordinary-polynomial implementation of either branch runs into the same obstruction and needs a shift trick (or Thm. 6) per branch. Net cost: **$4d_{\text{JA}}$** controlled-$W$ calls under shift tricks with doubled residual $W^{d_{\text{JA}}}$ factors (one per branch, different phases), or $2d_{\text{JA}}$ + $2d_{\text{JA}}$ = $4d_{\text{JA}}$ controlled operations under Thm. 6, plus LCU overhead. **Not recommended**: the draft's intuition ("each branch is a real function, so ordinary polynomial should suffice") is incorrect.

**Route 2 (reinstate MW Thm. 6; Laurent directly; controlled-$W$ plus controlled-$W^\dagger$):** mathematically clean. The Jacobi–Anger truncation $L_{d_{\text{JA}}}(z) \in \mathbb{C}[z,z^{-1}]$ is implemented directly via Fig. 2 of MW modified with $d_{\text{JA}}$ of the $2d_{\text{JA}}$ signal-operator boxes being controlled-$W^\dagger$ instead of controlled-$W$ (per MW Thm. 6's proof, Eqs. 45–53). Total controlled operations: $2d_{\text{JA}}$ ($d_{\text{JA}}$ of $W$, $d_{\text{JA}}$ of $W^\dagger$). This matches MW Cor. 8's cost bound literally. The only cost is textual: the §III pitch ("no $U^\dagger$") must be softened in the thesis subsection with a one-line acknowledgment that Thm. 6 is used specifically for Hamiltonian simulation.

**Recommendation**: **commit to Route 2 in v6.** The thesis does not depend on the "no-$U^\dagger$" property of GQSP anywhere in Chapter 5, and the hardware cost of controlled-$W^\dagger$ is identical to controlled-$W$ (both are calls to the same $2m+n$-qubit oracle with a conjugation at the PREP-sandwich level). Route 2 is the only route that delivers MW Cor. 8's advertised cost and a clean "$\ee^{-\ii\delta B_a}$" output — both prerequisites for Chapter 5's cost analysis and for the `CoherentStep` correctness claim.

**What this changes in v5 → v6:**

1. **Body, *Walk operator* paragraph** (unchanged): keep as is — the walk $W$ is MW's PREP-sandwiched walk regardless of route.
2. **Body, *GQSP layer and polynomial class* paragraph**: add one sentence: *"For Hamiltonian simulation, which needs a function of $W$ of mixed parity in $\theta$, we additionally invoke MW Thm. 6, which extends the class to Laurent polynomials $P'(z) = z^{-k} P(z)$ at the cost of replacing $k$ of the $d$ controlled-$W$ calls in Fig. `circ:gqsp` with controlled-$W^\dagger$ calls."*
3. **Body, *Target polynomial* paragraph**: drop the shift trick and the $z^{d_{\text{JA}}}$ argument entirely. Write: *"The GQSP circuit implements the Laurent polynomial $L_{d_{\text{JA}}}(z) \in \mathbb{C}[z,z^{-1}]$ directly via MW Thm. 6 with $d = 2d_{\text{JA}}$ and $k = d_{\text{JA}}$; $P(z) := z^{d_{\text{JA}}} L_{d_{\text{JA}}}(z)$ is the underlying ordinary polynomial, and the Thm. 6 conversion $P'(z) = z^{-d_{\text{JA}}}P(z) = L_{d_{\text{JA}}}(z)$ is the Laurent target."*
4. **Body, *Default $d_{\text{JA}}=1$* paragraph**: $P_1(\ee^{\ii\theta}) = L_1(\ee^{\ii\theta}) = J_0(\delta\alpha) - 2\ii J_1(\delta\alpha)\cos\theta$; this is a Laurent polynomial of degree 1, implemented via 1 controlled-$W$ call and 1 controlled-$W^\dagger$ call per layer, plus three single-qubit rotations. Small-$\delta\alpha$ expansion: $P_1 \approx 1 - \ii\delta\alpha\cos\theta \approx \ee^{-\ii\delta\alpha\cos\theta} + \mathcal{O}((\delta\alpha)^2)$ **exactly** — no trailing $\ee^{\ii\theta}$ factor.
5. **Body, *Closing cost sentence***: "$2d_{\text{JA}}$ controlled operations ($d_{\text{JA}}$ controlled-$W$ + $d_{\text{JA}}$ controlled-$W^\dagger$)", matching MW Cor. 8.
6. **Writing Notes, *Why Thm. 6 is wrong for v5* section**: delete or rewrite as *"Why Thm. 6 is needed for v6"*, acknowledging that §III's "no-$U^\dagger$" pitch applies to the pure GQSP framework (Thm. 3, Cor. 5) but Hamiltonian simulation (Cor. 8) needs Thm. 6.
7. **Figure `circ:gqsp`**: the redraw (draft Suggestion 1) gains a new element: $k = d_{\text{JA}}$ of the per-layer signal operators become controlled-$W^\dagger$ instead of controlled-$W$. Option A of the draft's suggested redraw ("expose SELECT + PREP-sandwich") is still the right choice; in the Laurent variant, $k$ of the layers flip the sign of the innermost SELECT.
8. **`2_methods.tex:1731`**: unchanged from Suggestion 2's current recommendation, modulo the degree rename and the $\cos$ swap. Cost accounting: "$d_{\text{JA}}$ controlled-$U_{B_a}$ + $d_{\text{JA}}$ controlled-$U_{B_a}^\dagger$ per coherent step", or equivalently "$d_{\text{JA}}$ calls to $W$ + $d_{\text{JA}}$ calls to $W^\dagger$".

---

## Cross-check of items 1–12 in the task

1. **Walk operator sign, MW "1" vs "I" notation, $U_B$ identity.** MW Cor. 8 reads *"$W = -(1 - 2\,\text{PREPARE}|0\rangle\langle 0|\text{PREPARE}^\dagger)\text{SELECT}$"* (pg. 11). Draft: $W = -(I - 2\,\text{PREP}|0\rangle\langle 0|\text{PREP}^\dagger)\text{SELECT}$. MW's "1" is identity; draft's "$I$" is the standard cleanup. **Sign ✓, notation ✓.** The identity $U_B = (\text{PREP}^\dagger\otimes I)\,\text{SELECT}\,(\text{PREP}\otimes I)$ is standard for block encodings of $B = \sum_j \alpha_j U_j$ with $\alpha_j \geq 0$, and gives $\bra{0^b}U_B\ket{0^b} = \sum_j (\alpha_j/\alpha) U_j = B/\alpha$. **Identity ✓.** SELECT lives on ancilla $\otimes$ system; PREP lives on ancilla only — tensor structure in the draft is correct.

2. **$W_\lambda$ matrix and basis.** MW Eq. (64): $W_{\lambda_j} = \begin{pmatrix}\lambda_j/\alpha & \sqrt{1-\lambda_j^2/\alpha^2}\\ -\sqrt{1-\lambda_j^2/\alpha^2} & \lambda_j/\alpha\end{pmatrix}$, with MW's basis convention "taking the former [PREP$|0\rangle|\lambda_j\rangle$] to $|0\rangle$ and the latter [$|\phi^-_j\rangle$] to $|1\rangle$" (pg. 11, above Eq. 64). Draft: exact same matrix in basis $(\text{PREP}|0\rangle|\lambda\rangle,|\phi^-_\lambda\rangle)$. **Basis ordering ✓, matrix ✓.** Eigenvalue derivation: trace $= 2\lambda/\alpha$, det $= 1$, so $\chi(x) = x^2 - (2\lambda/\alpha)x + 1$, roots $x = \lambda/\alpha \pm \ii\sqrt{1-\lambda^2/\alpha^2} = \ee^{\pm\ii\theta}$ with $\cos\theta = \lambda/\alpha$. **Eigenvalues ✓.**

3. **Jacobi–Anger identity.** MW Eq. (62): $\ee^{\ii t\cos\theta} = \sum_n \ii^n J_n(t) \ee^{\ii n\theta}$. Substituting $t = -\delta\alpha$: $\ee^{-\ii\delta\alpha\cos\theta} = \sum_n \ii^n J_n(-\delta\alpha)\ee^{\ii n\theta}$. Using $J_n(-z) = (-1)^n J_n(z)$: $= \sum_n \ii^n (-1)^n J_n(\delta\alpha)\ee^{\ii n\theta} = \sum_n (-\ii)^n J_n(\delta\alpha)\ee^{\ii n\theta}$. **Draft ✓.** (Minor: see m2 for the attribution tidy-up.)

4. **$P_2$ derivation.** Full check:
   - $(-\ii)^{-1} = 1/(-\ii) = \ii$ **✓**.
   - $J_{-1}(z) = -J_1(z)$ **✓** (standard Bessel parity).
   - $L_1(\ee^{\ii\theta}) = \ii\cdot(-J_1)\ee^{-\ii\theta} + J_0 + (-\ii)J_1\ee^{\ii\theta} = -\ii J_1(\ee^{-\ii\theta}+\ee^{\ii\theta}) + J_0 = J_0 - 2\ii J_1\cos\theta$ **✓**.
   - Multiply by $z = \ee^{\ii\theta}$: $J_0\ee^{\ii\theta} - 2\ii J_1\cos\theta \cdot\ee^{\ii\theta} = J_0\ee^{\ii\theta} - \ii J_1(\ee^{2\ii\theta}+1) = -\ii J_1\ee^{2\ii\theta} + J_0\ee^{\ii\theta} - \ii J_1$. Coefficients $(c_0, c_1, c_2) = (-\ii J_1, J_0, -\ii J_1)$ **✓**.
   - Small-$\delta\alpha$: $J_0(z) = 1 + \mathcal{O}(z^2)$, $2J_1(z) = z + \mathcal{O}(z^3)$; $P_2(\ee^{\ii\theta}) = \ee^{\ii\theta}(1 + \mathcal{O}((\delta\alpha)^2) - \ii\delta\alpha\cos\theta + \mathcal{O}((\delta\alpha)^3)) = \ee^{\ii\theta}\,(1 - \ii\delta\alpha\cos\theta) + \mathcal{O}((\delta\alpha)^2) = \ee^{\ii\theta}\,\ee^{-\ii\delta\alpha\cos\theta} + \mathcal{O}((\delta\alpha)^2)$ **✓** (the approximation is correct; what the draft misreads is whether the $\ee^{\ii\theta}$ factor is benign — see F1/M2).

5. **Coefficient formula.** At general $d_{\text{JA}}$: $c_m = [z^m]\,z^{d_{\text{JA}}}L_{d_{\text{JA}}}(z) = [z^{m - d_{\text{JA}}}]\,L_{d_{\text{JA}}}(z) = (-\ii)^{m-d_{\text{JA}}} J_{m-d_{\text{JA}}}(\delta\alpha)$ for $-d_{\text{JA}} \leq m - d_{\text{JA}} \leq d_{\text{JA}}$, i.e. $0 \leq m \leq 2d_{\text{JA}}$. **✓**. At $d_{\text{JA}} = 1$: $c_0 = (-\ii)^{-1}J_{-1} = \ii\cdot(-J_1) = -\ii J_1$ ✓; $c_1 = J_0$ ✓; $c_2 = -\ii J_1$ ✓.

6. **Bessel tail.** See m1 above: draft's bound is asymptotically correct at big-$\mathcal{O}$, slightly loose relative to MW Eq. (63).

7. **Degree bound.** MW Thm. 7: $\mathcal{O}(t + \log(1/\epsilon)/\log\log(1/\epsilon))$ controlled-$U$ for $\ee^{\ii t\cos H}$ or $\ee^{\ii t\sin H}$ (pg. 11). With $t \to \delta\alpha$: $d_{\text{JA}} = \mathcal{O}(\delta\alpha + \log(1/\varepsilon_{\text{QSP}})/\log\log(1/\varepsilon_{\text{QSP}}))$. **Draft ✓.**

8. **Cost accounting.** See M1 above. The draft's "$2d_{\text{JA}}$ matches MW's Cor. 8 bound" is imprecise; MW's cost is $d_{\text{JA}}$ calls to $W$ under Thm. 6 (Laurent), while the ordinary-polynomial shift trick costs $2d_{\text{JA}}$ calls to $W$. This is exactly the trade-off F1 flagged.

9. **Citations.** 
   - `motlagh2024generalized` and `chen2023efficient`: confirmed in use in `2_methods.tex:1608`, `2_methods.tex:1733`, etc. — pre-existing in the user's `references.bib` (not visible in the local repo, but grep confirms existing `\cite{...}` uses). **✓.**
   - `berntsonSunderhauf2025complementary`: new. ArXiv 2406.04246 confirmed on the PDF; DOI 10.1007/s00220-025-05302-9 confirmed on the CMP front page; volume 406, article 161 confirmed. **`number = {7}` is incorrect** — see m3.
   - `lowChuang2019qubitization`, `haah2019product`, `dong2021efficient` dropped: no orphan references in v5 body; draft's accounting is correct. **✓ (dropped cleanly).**

10. **Notation consistency — $d_{\text{JA}}$ vs $d$.** `d_{\text{JA}}` is used throughout v5 body; $d$ (bare) appears only in two places: the *GQSP layer and polynomial class* paragraph at line 43 ("degree $d$ … $d+1$ single-qubit rotations … $R(\theta_k,\phi_k,0)$, $k=1,\dots,d$") where $d$ is the GQSP circuit depth (still generic, not yet specialised); and in the degree-bound paragraph where it transitions cleanly to $d_{\text{JA}}$. **Rename is clean** except that one line's $d$ is ambiguous — it's being used as "the circuit depth = polynomial degree" before the Jacobi–Anger specialisation, but in v5 the final depth is $2d_{\text{JA}}$. The draft's Suggestion 3 ("Note the $2d_{\text{JA}}+1 = 3$ rotation pairs, not $d_{\text{JA}}+1 = 2$") covers this, but the body could benefit from one sentence tying the generic $d$ of the *GQSP layer* paragraph to $2d_{\text{JA}}$ explicitly: "In our application below, this generic $d$ specialises to $d = 2d_{\text{JA}}$".

11. **Thesis-writing rule compliance.** Markdown headers ✓, inline `$math$` ✓, display `$$...$$` with HTML-comment labels ✓ (`<!-- \label{eq:qsp-walk} -->` etc.), `[Author Year] [CITE: bibkey]` format ✓, Citations section at bottom ✓, Suggestions section ✓, Writing Notes section ✓, Diff from v4 section ✓. **Fully compliant.**

12. **Cross-reference labels.** 
    - `sec:prelim-LCU`: `1_preliminaries.tex:1037` **✓**.
    - `circ:gqsp`: `1_preliminaries.tex:1112` **✓**.
    - `alg:coh`: `2_methods.tex:1540` **✓**.
    - `eq:B-block-encoding`: `2_methods.tex:1702` **✓**.
    - `eq:b_plus-s-eta`: `2_methods.tex:523` **✓**.

---

## What's good

- The draft has the right amount of self-criticism: it flags its own `[CHECK]` item instead of hiding it. The honesty is preserved.
- The $P_2$ derivation in Writing Notes is a model of clarity: $(-\ii)^{-1}$, $J_{-1} = -J_1$, Euler sum, multiply by $z$, read off coefficients. Every sign is justified inline. Keep this style for v6.
- Dropping `lowChuang2019qubitization`, `haah2019product`, `dong2021efficient` is the right call — v5's citation graph is now lean and single-source.
- The Notation choices bullet list is exemplary: every symbol that could clash is flagged ($\lambda$ vs $\lambda_0$; $d$ vs $d_{\text{JA}}$; $\text{PREP}|0\rangle$ vs $\ket{0^b}$; $\mathbb{T}$).
- Suggestion 1 (figure redraw) correctly identifies the most impactful change required for self-consistency and gives two implementation options.
- Suggestion 2 (methods-side reconciliation) is the right level of detail, with line numbers and exact string replacements.

---

## Summary Scorecard

| Category | Rating | Notes |
|----------|--------|-------|
| Logical correctness | MAJOR | The shift-trick single-circuit route implements $W^{d_{\text{JA}}}\ee^{-\ii\delta B}$, not $\ee^{-\ii\delta B}$. Fixable by committing to Route 2 (Thm. 6). |
| Notation consistency | OK | $d_{\text{JA}}$ rename is clean; one line (43) could tie generic $d$ to $2d_{\text{JA}}$ explicitly. |
| Citations | MINOR | Berntson–Sünderhauf `number = {7}` is spurious; drop. All other citations verified. |
| Thesis consistency | OK | Methods-side reconciliation is correctly scoped in Suggestion 2. Tilde on $B_a$ vs $\tilde{B}_a$ is a minor Methods concern. |
| Exposition quality | OK | Body reads well; Writing Notes are comprehensive and helpful. |

---

## Verdict

**v5 is NOT ready to insert.** Two blockers:

1. **Must-fix (FATAL F1 / MAJOR M1–M3)**: commit to **Route 2 (MW Thm. 6, Laurent + controlled-$W^\dagger$)** in v6 and rewrite the body paragraphs accordingly (see [CHECK] adjudication, v5→v6 items 1–7). This is the only route that delivers the literal MW Cor. 8 cost bound and a clean `$\ee^{-\ii\delta B}$` output. The textual cost is two sentences acknowledging that Thm. 6 is used specifically for Hamiltonian simulation.

2. **Must-fix (self-consistency)**: redraw `circ:gqsp` per Suggestion 1 in the draft, updated to reflect $k = d_{\text{JA}}$ of the controlled-$W$ calls being controlled-$W^\dagger$ under Route 2. Without the redraw the figure and body disagree on the walk operator.

3. **Should-fix (MINOR m1–m6)**: Bessel tail $+1$ shift (m1), bibtex `number` removal (m3), caption compression (m5), memo update (m4). All cosmetic; can be batched in a polish pass.

With items 1 and 2 addressed, v6 should be ready to insert. Suggestion 2 (methods reconciliation) is a follow-on and does not block v6 insertion into `1_preliminaries.tex`.
