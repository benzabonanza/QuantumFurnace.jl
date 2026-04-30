# Final Review: QSP Subsection v6 (Route 2, MW Thm. 6 Laurent)

**Reviewed**: `drafts/qsp-subsection-v6.md`
**Against**: `drafts/qsp-subsection-v5-review.md` §"What this changes in v5 → v6"
**Primary source verified**: `supplementary-informations/Motlagh and Wiebe - 2024 - Generalized Quantum Signal Processing.pdf` — Thm. 3 (p.4, Eq. 7), Cor. 5 (p.6, Eqs. 31–32), Thm. 6 (p.8, Eqs. 45–53), Thm. 7 (p.11, Eq. 62), Cor. 8 (p.11, Eqs. 64–66, p.12 Eqs. 67–68), Alg. 1 (p.9), §IV Eqs. 58–60.
**Secondary source verified**: `supplementary-informations/Berntson and Sünderhauf - 2025 - Complementary Polynomials in Quantum Signal Processing.pdf` (Commun. Math. Phys. 406:161).
**Insertion target**: `supplementary-informations/1_preliminaries.tex:1068–1113`.
**Methods-side target**: `supplementary-informations/2_methods.tex:1606–1740`.
**Date**: 2026-04-23

**Overall assessment**: v6 resolves the FATAL F1 and all MAJOR (M1–M3) issues from the v5 review cleanly; the MW Thm. 6 construction is invoked honestly with correct indexing, the Laurent polynomial is genuinely what's implemented, the cost accounting matches MW Cor. 8 literally, and the $L_1$ default expansion no longer carries the trailing $\ee^{\ii\theta}$ factor. All m1–m6 minor fixes are applied, all cross-references resolve, and thesis-writing rule compliance is clean. Two issues remain, both MINOR: a loose symmetric-to-which-of-the-two-slot-patterns claim at Writing Notes line 131, and a mild imprecision in the Thm. 6 "commutation" sentence at line 131 that is worth tightening. Neither blocks insertion. **READY TO INSERT**, with two tiny polish items (listed at the end) optional.

---

## FATAL

None.

---

## MAJOR

None.

---

## MINOR

### m_v6_1: The "or vice versa" claim about flipping first-vs-last $d_{\text{JA}}$ slots is slightly looser than MW's construction.

- **Location**: Writing Notes, "Derivation of $L_1$" section, line 139:
  > "the Thm. 6 circuit deploys slot 1 as controlled-$W^\dagger$ and slot 2 as controlled-$W$ (or vice versa; the two arrangements differ by a global phase under the PREP-sandwich, both implement $L_1(W)$)"

  Also in the "Thm. 6 circuit-level identity" section, line 131:
  > "The angle vector, computed from $P$ via Alg. 1, is identical between the shift-trick circuit (which implements $P(W) = W^{d_{\text{JA}}} L_{d_{\text{JA}}}(W)$) and the Thm. 6 circuit"

- **Problem**: MW Thm. 6's canonical construction (Eq. 46) has the flipped signal operators at the **last** $k$ of the $d$ slots — specifically, slots $j = d-k+1, \dots, d$. MW's body text above Thm. 6 (p.8) does say "replacing any $k$ instances of $A$", suggesting a distributed choice also works, but MW's proof (Eqs. 49–53) only formally derives the result for the specific "last $k$" pattern because it needs $(I\otimes U^\dagger)^k$ to pull through the entire rotation–signal product to the left. Whether distributed patterns give the same $P'(z) = z^{-k}P(z)$ in the top-left block (exactly, without phase or coefficient modifications) requires a short check that v6 does not carry out.

  The "differ by a global phase" claim is therefore plausible (commutations would absorb into the diagonal structure) but not rigorously established within MW's proof. For the $d_{\text{JA}} = 1$, $d = 2$ default, both orderings (slot 1 = $W^\dagger$, slot 2 = $W$) versus (slot 1 = $W$, slot 2 = $W^\dagger$) likely give the same top-left block up to a diagonal conjugation, but this deserves either (a) a brief direct computation or (b) a citation of MW's exact wording ("any $k$ instances").

- **Why it matters**: low — the body text of v6 consistently uses the canonical MW pattern ($d_{\text{JA}}$ controlled-$W$ calls followed by $d_{\text{JA}}$ controlled-$W^\dagger$ calls), which is exactly MW Eq. 46. Only the Writing Notes ambiguate; since Writing Notes are for the user's own reference and not part of the thesis body, this is purely a cleanliness issue.

- **Fix**: replace "or vice versa; the two arrangements differ by a global phase under the PREP-sandwich, both implement $L_1(W)$" with "following MW Eq. 46 the canonical choice is slot 1 $= W$, slot 2 $= W^\dagger$ (matching MW's 'last $k$' pattern with $d=2, k=1$); MW's preceding body text asserts that any distribution of $k$ controlled-$W^\dagger$ gates across the $d$ slots realises the same $P'$, though the proof is given only for the canonical pattern." This is one sentence and removes the ambiguity.

### m_v6_2: The Thm. 6 commutation sentence at Writing Notes line 131 is phrased confusingly.

- **Location**: Writing Notes, "Thm. 6 circuit-level identity" section, line 129:
  > "the substitution $A' = (I \otimes U^\dagger)\,A$ (MW Eq. 48), together with the commutation $(I \otimes U^\dagger) R(\theta, \phi, 0) = R(\theta, \phi, 0)(I \otimes U^\dagger)$ and $(I \otimes U^\dagger) A = A\,(I \otimes U^\dagger)$ (since the SU(2) rotation acts only on the QSP ancilla and $A$ acts trivially on $I \otimes U^\dagger$ except through the $\ket{1}\!\bra{1}\otimes U$ term, where the $U\cdot U^\dagger = I$ commutation holds up to $U^\dagger$ passing through trivially)"

- **Problem**: the "up to $U^\dagger$ passing through trivially" parenthetical is not correct as written. The actual mechanism is: $A = \ket{0}\!\bra{0}\otimes I + \ket{1}\!\bra{1}\otimes U$. Then $(I\otimes U^\dagger)A = \ket{0}\!\bra{0}\otimes U^\dagger + \ket{1}\!\bra{1}\otimes U^\dagger U = \ket{0}\!\bra{0}\otimes U^\dagger + \ket{1}\!\bra{1}\otimes I$, while $A(I\otimes U^\dagger) = \ket{0}\!\bra{0}\otimes U^\dagger + \ket{1}\!\bra{1}\otimes UU^\dagger = \ket{0}\!\bra{0}\otimes U^\dagger + \ket{1}\!\bra{1}\otimes I$. They are equal because $U$ and $U^\dagger$ commute with each other on the target register (a unitary commutes with its inverse), and tensor with $\ket{k}\bra{k}$ factors separately. The sentence as written reads like there's some subtle "trivial passing" involved, but there isn't: standard $UU^\dagger = U^\dagger U = I$ on a single register.

- **Why it matters**: low — the conclusion is correct and matches MW's Eqs. 49–51. Only the verbal justification of why commutation holds is hand-wavy.

- **Fix**: replace the parenthetical with "the commutation holds because $U^\dagger U = UU^\dagger = I$ on the target register and $\ket{k}\bra{k}$ are mutually orthogonal projectors".

### m_v6_3: Notation polish — generic $d$ in the polynomial-class paragraph is now tied to $d = 2d_{\text{JA}}$ correctly, but the tie-in sentence is a parenthetical and reads slightly disjointed.

- **Location**: Body, "GQSP layer and polynomial class" paragraph, lines 41–45:
  > "the opening rotation $R(\theta_0,\phi_0,\lambda_0)$ and $d$ further two-parameter rotations $R(\theta_k,\phi_k,0)$, $k = 1,\dots,d$ (the third parameter $\lambda_0 \in \mathbb{R}$ of the opening rotation, distinct from the spectral variable $\lambda$, absorbs a single global phase per [Motlagh & Wiebe 2024, Thm. 3] [CITE: motlagh2024generalized]); in our application below, this generic $d$ specialises to $d = 2d_{\text{JA}}$."

- **Problem**: the specialisation sentence is tacked on at the end of a long parenthetical-chain sentence, making it easy to miss. This is a readability issue, not a correctness one.

- **Why it matters**: very low — this is the kind of thing a copy-editor would tweak; the current text is understandable on careful reading.

- **Fix** (optional): promote the specialisation to its own sentence immediately after the rotation description, e.g. after "two-parameter rotations $R(\theta_k, \phi_k, 0)$, $k = 1, \dots, d$.", insert "In the Hamiltonian-simulation application below, this generic circuit depth specialises to $d = 2d_{\text{JA}}$ (see the target-polynomial paragraph)." Defensible to keep as-is.

---

## Item-by-Item Verification Table

| # | Item (from task) | Verdict | Justification |
|---|-----|--------|---------------|
| 1 | FATAL F1 resolution — Laurent $L_{d_{\text{JA}}}$ via MW Thm. 6, no residual $W^{d_{\text{JA}}}$ | **PASS** | Body §"Target polynomial via Jacobi–Anger" line 59–63 explicitly applies Thm. 6 with $d=2d_{\text{JA}}, k=d_{\text{JA}}$ to produce $P'(z) = z^{-d_{\text{JA}}}P(z) = L_{d_{\text{JA}}}(z)$ directly; line 63 confirms $L_{d_{\text{JA}}}(W_\lambda) = \ee^{-\ii\delta\alpha\cos\theta_\lambda} I_{2\times 2} + \mathcal{O}(\varepsilon_{\text{QSP}})$ with no trailing $W^{d_{\text{JA}}}$ factor. $P$ is clearly flagged as "the bookkeeping polynomial that MW Alg. 1 runs on" (line 59). |
| 2 | MAJOR M1 cost accounting — $d_{\text{JA}}$ controlled-$W$ + $d_{\text{JA}}$ controlled-$W^\dagger$ = $2d_{\text{JA}}$ ops | **PASS** | Body §"Degree bound and cost" line 69 says exactly "$d_{\text{JA}}$ controlled-$W$ calls and $d_{\text{JA}}$ controlled-$W^\dagger$ calls — $2d_{\text{JA}}$ controlled operations in total". Same split stated in closing sentence (line 87). Controlled-$W^\dagger$ cost equivalence explicitly justified: "both are implemented as an open-controlled call to the same $2m+n$-qubit block-encoding oracle, with a conjugation applied at the PREP-sandwich level of $W$; no hardware premium arises from the Thm. 6 substitution" (line 69). MW Cor. 8's $\mathcal{O}(\alpha\delta + \log/\log\log)$ citation is paired with the reading "applications of $W$ or $W^\dagger$" (line 69 and Citations §, line 94). |
| 3 | MAJOR M2 $L_1$ Laurent, no trailing $\ee^{\ii\theta}$ factor | **PASS** | Body §"Small-$\delta\alpha$ default" line 73 writes $L_1(\ee^{\ii\theta}) = J_0 - 2\ii J_1\cos\theta$ (Laurent of bilateral degree 1), implemented by "one controlled-$W$ call, one controlled-$W^\dagger$ call, and three single-qubit rotations" (line 75). Small-$\delta\alpha$ expansion at line 77 gives $L_1 = \ee^{-\ii\delta\alpha\cos\theta} + \mathcal{O}((\delta\alpha)^2)$ with explicit "**exactly** — no residual $\ee^{\ii\theta}$ factor" emphasis. The underlying ordinary $P$ is parenthetically mentioned (line 75) but not mistaken for what's implemented. |
| 4 | MAJOR M3 Route 1 acknowledgment | **PASS** | Writing Notes §"The shift-by-$z^{d_{\text{JA}}}$ step (resolved and closed)" line 123 explicitly addresses Route 1: "The alternative Route 1 (two GQSP calls with $\cos + \ii\sin$) does not rescue ordinary polynomials either: each branch $\cos(\delta\alpha\cos\theta)$ and $\sin(\delta\alpha\cos\theta)$ is even in $\theta$, so the same $P(z) = P(1/z)$ obstruction applies per branch, and Laurent polynomials would still be needed per branch (at double the total cost of Route 2)". |
| 5a | MINOR m1 — tightened Bessel tail bound $(\ee\delta\alpha/(2(d_{\text{JA}}+1)))^{d_{\text{JA}}+1}$ | **PASS** | Body line 55: "$\|\ee^{-\ii\delta\alpha\cos\theta} - L_{d_{\text{JA}}}\|_{\infty,\mathbb{T}} \leq 2\sum_{|k|>d_{\text{JA}}}|J_k(\delta\alpha)| \leq \mathcal{O}\bigl((\ee\delta\alpha/(2(d_{\text{JA}}+1)))^{d_{\text{JA}}+1}\bigr)$" — exact form requested. Same bound appears in the $\|P\|_{\infty,\mathbb{T}}$ statement at line 59. |
| 5b | MINOR m2 — Jacobi–Anger uses $J_n(-z)=(-1)^n J_n(z)$ directly, not $J_{-k}$ | **PASS** | Body line 51: "the sign pattern follows from $t \to -\delta\alpha$ in MW's identity $\ee^{\ii t\cos\theta} = \sum_k \ii^k J_k(t)\,\ee^{\ii k\theta}$ together with $J_n(-z) = (-1)^n J_n(z)$, giving $\ii^n (-1)^n = (-\ii)^n$ directly — **no invocation of $J_{-k}$ needed**." Note the $J_{-1}$ reference still appears in the $L_1$ Writing-Notes derivation (line 139), which is a separate step about constructing $L_1$ from the truncated sum; this is consistent with the v5 review's guidance (m2 is about the general Jacobi–Anger sign derivation, not the coefficient-level computation at $d_{\text{JA}}=1$). |
| 5c | MINOR m3 — `number = {7}` removed from Berntson–Sünderhauf bibtex | **PASS** | Bibtex block (lines 97–109): volume 406, pages 161, no number field. DOI 10.1007/s00220-025-05302-9. Clean. |
| 5d | MINOR m4 — memo reference updated for v6 | **PASS** | Suggestion 2 last bullet (line 176): "once the methods edit lands, that memo can be updated to reflect v6's punch list (not v4/v5's)." This is the exact wording requested. |
| 5e | MINOR m5 — figure caption compressed to two sentences | **PASS** | Suggestion 1 (line 167) suggested caption: *"Degree-$2d_{\text{JA}}$ GQSP realising $\ee^{-\ii\delta B}$ from the walk operator $W$ of (Eq.~\ref{eq:qsp-walk}) (MW Cor. 8), with $k = d_{\text{JA}}$ of the signal-operator calls promoted to controlled-$W^\dagger$ per MW Thm. 6. Post-selecting both ancillas on $\ket{0}_\mathrm{QSP}\otimes\text{PREP}\ket{0}$ yields the Laurent polynomial $L_{d_{\text{JA}}}(W)$ realising $\ee^{-\ii\delta B}$ on the qubitization subspace."* — two sentences, matching the v5 review's m5 suggestion in spirit and length. |
| 5f | MINOR m6 — $B_a$ vs $\tilde B_a$ note added to Suggestion 2 | **PASS** | Suggestion 2 bullet "Note on $B_a$ vs $\tilde B_a$" (line 175) spells out the Trotterized-vs-exact distinction, references `2_methods.tex:1681` and `:1700` by line number, and offers two remediation paths (use $\tilde{B}_a$ consistently, or add a clarifying sentence that the coherent step block-encodes $B_a$ up to Trotter error). |
| 5g | Notation — generic $d$ in polynomial-class paragraph tied to $d = 2d_{\text{JA}}$ | **PASS** | Line 41: "in our application below, this generic $d$ specialises to $d = 2d_{\text{JA}}$." (though the phrasing is slightly buried in a long sentence — see m_v6_3 above for an optional polish). |
| 6 | No new errors introduced; v5 correct math preserved | **PASS** | Walk operator (Eq. `eq:qsp-walk`, lines 28–32) verbatim from v5. $W_\lambda$ matrix (lines 33–34), trace/det, $\cos\theta = \lambda/\alpha$ (line 37) verbatim. Jacobi–Anger (Eq. `eq:jacobi-anger`, lines 49–50) unchanged modulo m2's attribution edit. $c_m = (-\ii)^{m-d_{\text{JA}}} J_{m-d_{\text{JA}}}(\delta\alpha)$ appears at line 57 as part of the definition of $P$. MW Cor. 8 cost statement (line 69) unchanged. Berntson–Sünderhauf use (Q-finding only) unchanged (line 85). |
| 7 | Thm. 6 construction check ($k = d_{\text{JA}}$ of $d = 2d_{\text{JA}}$, resulting Laurent range) | **PASS** | Against MW Eqs. 45–46 (p.8): $P'(z) = z^{-k} P(z)$ with $\deg P \leq d$, $0 \leq k \leq d$. MW's body text (p.8, above Thm. 6) specifies that if $P(x) = \sum_{n=0}^d a_n e^{inx}$, then $P'(x) = \sum_{n=-k}^{d-k} a_{n+k} e^{inx}$. For v6: $d = 2d_{\text{JA}}$, $k = d_{\text{JA}}$, so $P'$ range is $n \in [-d_{\text{JA}}, d_{\text{JA}}]$ — exactly $L_{d_{\text{JA}}}$'s bilateral range. Verified. v6 body line 59: "Applying MW Thm. 6 with $d = 2d_{\text{JA}}$ and $k = d_{\text{JA}}$ converts $P$ to $P'(z) = z^{-d_{\text{JA}}}\,P(z) = L_{d_{\text{JA}}}(z)$" — correct. |
| 8 | $L_{d_{\text{JA}}}(W)$ on qubitization subspace genuinely equals $\ee^{-\ii\delta B} I + \mathcal{O}(\varepsilon)$ | **PASS** | The Laurent-palindromic check: $L_{d_{\text{JA}}}(1/z) = \sum_k (-\ii)^k J_k(\delta\alpha) z^{-k}$. Relabel $k \to -k$: $\sum_k (-\ii)^{-k} J_{-k}(\delta\alpha) z^k = \sum_k \ii^k \cdot (-1)^k J_k z^k = \sum_k (-\ii)^k J_k z^k = L_{d_{\text{JA}}}(z)$. ✓. Therefore on the 2D qubitization subspace $L_{d_{\text{JA}}}(W_\lambda) = \operatorname{diag}(L_{d_{\text{JA}}}(\ee^{\ii\theta_\lambda}), L_{d_{\text{JA}}}(\ee^{-\ii\theta_\lambda})) = L_{d_{\text{JA}}}(\ee^{\ii\theta_\lambda}) I_{2\times 2}$, matching MW Eq. 66 without approximation error beyond the Bessel truncation. v6 body line 63 states this cleanly and correctly. Additional small-$\delta\alpha$ check: $L_1 = 1 - \ii(\delta\alpha)\cos\theta + \mathcal{O}((\delta\alpha)^2)$ (using $J_0 = 1 + O(z^2)$, $2J_1 = z + O(z^3)$), matching line 77. |
| 9 | Citations and cross-references | **PASS** | `motlagh2024generalized` appears with appropriate theorem/cor numbers at lines 24 (GQSP intro), 26 (Cor. 8), 31 (Eq. 64), 35 (Eq. 64), 45 (Thm. 6, Eqs. 45–53), 47 (Eq. 62), 55 (Eq. 63), 65 (Thm. 7), 81 (Thm. 4), 85 (§IV Alg. 1). `berntsonSunderhauf2025complementary` appears only at line 85 for Q-finding. `chen2023efficient` not mentioned in v6 body (matches v5, where it was also absent from the QSP body — it's used only in Chapter 5). No orphan references to `lowChuang2019qubitization`, `haah2019product`, or `dong2021efficient`. All cross-refs resolve: `sec:prelim-LCU` → `1_preliminaries.tex:1037`, `circ:gqsp` → `1_preliminaries.tex:1112`, `alg:coh` → `2_methods.tex:1540`, `eq:B-block-encoding` → `2_methods.tex:1702`, `eq:b_plus-s-eta` → `2_methods.tex:523`. |
| 10 | Writing Notes and Suggestions consistency with v6 | **PASS** | Writing Notes: `[CHECK]` flag is replaced by "The shift-by-$z^{d_{\text{JA}}}$ step (resolved and closed)" (line 121), which explicitly references the v5 review. "Why Thm. 6 is needed for v6" section exists (line 115, replacing v5's "wrong for v5"). Factor-of-2 section updated for Route 2 (line 133). $L_1$ derivation is clean Laurent, with underlying ordinary $P$ listed for bookkeeping only (line 139). Suggestions: Suggestion 1 updated for Route 2 with $k = d_{\text{JA}}$ slots flipped (line 165). Suggestion 2 updated for Route 2 cost and expanded to include the $B_a$/$\tilde{B}_a$ bullet (line 175). v5's Suggestion 6 (resolve `[CHECK]`) is correctly dropped. Suggestions 3, 4, 5 preserved. |
| 11 | Thesis-writing rule compliance | **PASS** | Markdown headers (`#`, `##`, `###`) used for structure; prose in plain Markdown; inline math `$...$` and display math `$$...$$`; display labels via HTML comments (`<!-- \label{eq:qsp-walk} -->`, etc.); citations as `[Author Year] [CITE: bibkey]`; Citations section at bottom with bibtex fence; cross-references as `(Fig. `circ:gqsp`)`, `(Eq. `eq:qsp-walk`)`, `(§`sec:prelim-LCU`)`; opening "Insertion target" line; Writing Notes and Suggestions sections present. No LaTeX-only constructs misused. |

---

## New `[CHECK]` flags

None. The v5 `[CHECK]` is explicitly closed (negatively — the shift trick was not the right construction, v6 uses Thm. 6 instead). No new tensions have been introduced.

---

## What's good (worth preserving through any further edits)

- The "Thm. 6 circuit-level identity" Writing-Notes section (line 127) walks through MW Eqs. 48–53 explicitly, showing exactly how the substitution moves $(I\otimes U^\dagger)^k$ to the left and yields $P'(U) = U^{-k}P(U)$ in the top-left block. This is the clearest justification the draft could offer for the Thm. 6 application, and is much stronger than v5's equivalent (which didn't explain the construction at all because v5 wasn't using Thm. 6).
- The deliberate separation of "what's implemented" (the Laurent polynomial $L_{d_{\text{JA}}}$, a bilateral function of $\ee^{\ii\theta}$) from "what Alg. 1 runs on" (the ordinary polynomial $P = z^{d_{\text{JA}}} L_{d_{\text{JA}}}$ of degree $2d_{\text{JA}}$, which is the bookkeeping device for Alg. 1's rotation-angle recursion) is handled cleanly in the body (lines 55–63) and reinforced in the angle-finding paragraph (line 81) with the explicit parenthetical clarifying the distinction.
- The small-$\delta\alpha$ default paragraph's use of *"**exactly** — no residual $\ee^{\ii\theta}$ factor"* (line 77) is a model of precise, self-checking writing — it preemptively contrasts v6 with v5's $P_2$ expansion that had the trailing factor.
- The Bessel tail-bound tightening ($(\ee\delta\alpha/(2(d_{\text{JA}}+1)))^{d_{\text{JA}}+1}$) is applied consistently wherever the tail bound appears (lines 55 and 59).
- The Route 1 rebuttal at Writing Notes line 123 is tight: one sentence explaining why each $\cos/\sin$ branch has the same $P(z) = P(1/z)$ obstruction.
- The Citations-section annotation for `motlagh2024generalized` (line 94) lists every theorem/corollary invoked, with the note that "applications of $W$" in MW Cor. 8 should be read as "$W$ and $W^\dagger$". This is exactly the reading v6's cost accounting requires, and flagging it in the Citations section is the right place.
- The `[FIX-NEXT]` tags in Suggestions 1 and 2 correctly signal the non-negotiable follow-ups (figure redraw, methods reconciliation); both are outside v6's scope but required for downstream self-consistency.

---

## Summary Scorecard

| Category | Rating | Notes |
|----------|--------|-------|
| Logical correctness | **OK** | Route 2 (Thm. 6 Laurent) correctly invoked; $L_{d_{\text{JA}}}(W)$ genuinely equals $\ee^{-\ii\delta B} I_{2\times 2} + \mathcal{O}(\varepsilon)$ on the qubitization subspace (palindromic Laurent verified). No residual $W^{d_{\text{JA}}}$. |
| Notation consistency | **OK** | $d = 2d_{\text{JA}}$, $k = d_{\text{JA}}$ indexing consistent throughout; $d_{\text{JA}}$ used uniformly in body; generic $d$ flagged and tied to $2d_{\text{JA}}$. |
| Citations | **OK** | `motlagh2024generalized` with theorem/cor precision; `berntsonSunderhauf2025complementary` bibtex clean; no orphans. |
| Thesis consistency | **OK** | Methods-side reconciliation (Suggestion 2) is the correct scope; no insertion-blocking issues in `1_preliminaries.tex`. |
| Exposition quality | **OK** | Body reads clearly; Writing Notes thorough and self-critical; two small polish items (m_v6_1, m_v6_2) are purely cosmetic. |

---

## Verdict: **READY TO INSERT**

v6 is ready for insertion into `supplementary-informations/1_preliminaries.tex` lines 1068–1113. All FATAL / MAJOR items from v5 are resolved, all minor items from the v5 review (m1–m6 and notation item 10) are applied, and no new critical issues have arisen.

### Nice-to-haves (non-blocking, can be applied in a polish pass):

1. **[Optional, m_v6_1]** Tighten the "or vice versa" claim in Writing Notes line 139 (and similar in line 131) to explicitly follow MW Eq. 46's canonical "last $k$" ordering, or cite MW's "any $k$ instances" body-text claim explicitly.
2. **[Optional, m_v6_2]** Replace the parenthetical at Writing Notes line 129 ("up to $U^\dagger$ passing through trivially") with the cleaner justification: "because $U^\dagger U = UU^\dagger = I$ and $\ket{k}\bra{k}$ are mutually orthogonal projectors".
3. **[Optional, m_v6_3]** Move the "this generic $d$ specialises to $d = 2d_{\text{JA}}$" tie-in (line 41) out of the parenthetical chain and into its own sentence for readability.

### Required-but-out-of-scope-for-v6 follow-ups (covered in Suggestions 1–2):

- **Figure `circ:gqsp` redraw** at `1_preliminaries.tex:1071–1113`: Option A (expose SELECT + PREP-sandwich per Cor. 8, with $d_{\text{JA}}$ slots flipped to controlled-$W^\dagger$) is the recommended choice; caption rewrite is required either way. Without this, the figure and the v6 body disagree on the walk operator's internal structure.
- **Methods reconciliation** at `2_methods.tex:1606–1740`: lines 1609, 1731, 1733–1735, 1737, and the cost accounting all need the $\cos$/Laurent/$d_{\text{JA}}$ updates per Suggestion 2. Also the $B_a$ vs $\tilde{B}_a$ convention needs a decision (Suggestion 2 last bullet).

Neither follow-up blocks inserting v6 into the preliminaries file; they become necessary as soon as the reader cross-references Chapter 5.
