# Draft Review v3: Trotterization (Preliminaries subsection)

**Reviewed**: `/Users/bence/code/QuantumFurnace.jl/drafts/trotterization.md` (v3)
**Target insertion**: body of `\subsection{Trotterization.}` in `/Users/bence/code/QuantumFurnace.jl/supplementary-informations/1_preliminaries.tex` (between lines 1034 and 1037).
**Primary source verified against**: `/Users/bence/code/QuantumFurnace.jl/supplementary-informations/Childs et al. - 2021 - Theory of Trotter Error with Commutator Scaling.pdf` (PRX 11, 011020 (2021)).
**Downstream consumer**: `/Users/bence/code/QuantumFurnace.jl/supplementary-informations/2_methods.tex` line 881 onwards (`prop:trotter-diss` at 908, Strang bound at 929).
**Date**: 2026-04-21

**Overall assessment**: v2's core is preserved clean. The v3 trims are all clean — no orphan references — and the new unitarity sentence is correct (modulo one minor qualification about Hermiticity). However, the new $M$-count display (`eq:trotter-M-count-heisenberg`) carries an **algebraic error in the explicit constant** (factor-of-2), and the new *Principle* remark's quoted combinatorial constant $(2d)^p$ mixes two different graphs ("graph on summands" vs "graph on bonds") in a way that will confuse a careful reader. Both are fixable with a one-line edit.

---

## Severity: one [BLOCKER], one [MAJOR], a handful of [MINOR]/[NIT]s.

---

## [BLOCKER] B1. Algebraic error in the explicit $M$-count constant (`eq:trotter-M-count-heisenberg`)

- **Location**: draft lines 131–135, display `eq:trotter-M-count-heisenberg`.
- **Claim under review**: $M \ge \sqrt{15\cdot 54}\,|J|^{3/2}\,n^{1/2}/(\sigma^{3/2}\,\varepsilon^{1/2})$.
- **Actual algebra**: starting from the Strang bound of `prop:trotter-diss` at `2_methods.tex:929`,
  $$
  \varepsilon \;\ge\; \frac{\sqrt{15}\,\tilde\alpha^{(2)}_{\mathrm{comm}}}{M^2\,\sigma^3}
  \quad\xrightarrow{\ \tilde\alpha^{(2)}_{\mathrm{comm}}\le 54|J|^3 n\ }\quad
  \varepsilon \;\ge\; \frac{\sqrt{15}\,\cdot\,54\,|J|^3\,n}{M^2\,\sigma^3}.
  $$
  Solving for $M^2$ and then $M$:
  $$
  M^2 \;\ge\; \frac{\sqrt{15}\cdot 54\,|J|^3\,n}{\sigma^3\,\varepsilon}
  \quad\Longrightarrow\quad
  M \;\ge\; \underbrace{\sqrt{\sqrt{15}\cdot 54}}_{\approx\,14.46}\,\frac{|J|^{3/2}\,n^{1/2}}{\sigma^{3/2}\,\varepsilon^{1/2}}.
  $$
  The draft's stated constant is $\sqrt{15\cdot 54}=\sqrt{810}\approx 28.46$. The correct constant is $\sqrt{\sqrt{15}\cdot 54}=(15)^{1/4}\sqrt{54}\approx 14.46$. **Almost exactly a factor of 2 too large.**

  The slip is in the square root: the draft pulled the $\sqrt{15}$ *inside* a single square root, but the $\sqrt{15}$ is already a square root. Equivalently, $M^2 \ge \sqrt{15}\cdot K$ implies $M \ge 15^{1/4}\sqrt{K}$, not $M \ge \sqrt{15}\cdot\sqrt{K}$.

- **Why it matters**: this is the one explicit numerical constant the subsection provides. Getting it wrong is the worst kind of error in a preliminaries chapter: readers trust the explicit number and will propagate it into any resource-counting estimate downstream. The $\mathcal{O}$-form on the right-hand side of the display *is* correct — the big-O hides the constant — but the **equality with an explicit constant on the left is wrong**.
- **Suggested fix**: replace
  > $M \;\ge\; \frac{\sqrt{15 \cdot 54}\;|J|^{3/2}\,n^{1/2}}{\sigma^{3/2}\,\varepsilon^{1/2}} \;=\; \mathcal{O}\!\left(\frac{|J|^{3/2}\,n^{1/2}}{\sigma^{3/2}\,\varepsilon^{1/2}}\right)$

  with either
  > $M \;\ge\; \frac{15^{1/4}\sqrt{54}\;|J|^{3/2}\,n^{1/2}}{\sigma^{3/2}\,\varepsilon^{1/2}} \;=\; \mathcal{O}\!\left(\frac{|J|^{3/2}\,n^{1/2}}{\sigma^{3/2}\,\varepsilon^{1/2}}\right)$

  or, cleaner and slightly loose (if the explicit constant is not needed),
  > $M \;\ge\; \sqrt{\frac{\sqrt{15}\cdot 54\,|J|^3\,n}{\sigma^3\,\varepsilon}} \;=\; \mathcal{O}\!\left(\frac{|J|^{3/2}\,n^{1/2}}{\sigma^{3/2}\,\varepsilon^{1/2}}\right).$

  The second form shows the square root structure plainly and is what I'd recommend. Either way, flag the mistake item 7 in the Revision-notes block before the user copies it to `.tex`.

- **Sanity-check update to the Revision-notes block (line 168)**: change "$M \ge \sqrt{15\cdot 54}\,|J|^{3/2}\,n^{1/2}/(\sigma^{3/2}\varepsilon^{1/2})$" to "$M \ge 15^{1/4}\sqrt{54}\,|J|^{3/2}\,n^{1/2}/(\sigma^{3/2}\varepsilon^{1/2})$".

---

## [MAJOR] M1. *Principle: geometrically-local commutator scaling* — the $(2d)^p$ bound conflates the summand graph and the "inside-$H_k$" graph

- **Location**: draft line 78 (the new *Principle* remark).
- **Claim under review**: "on a graph of max degree $d$ the number of connected $(p+1)$-tuples anchored at a pivot vertex is at most $(2d)^p$, giving $n \cdot (2d)^p = \mathcal{O}(n)$ non-vanishing nested commutators in total."
- **Problem**: the Principle starts by fixing a decomposition $H = \sum_k H_k$ (indexed by $k$), with interaction graph of max degree $d$ on that index set, and then speaks of the $(p+1)$-nested commutator sum. The number of non-vanishing terms in $\sum_{k_1,\dots,k_{p+1}}\|\dots\|$ is a count of connected $(p+1)$-tuples on the graph of summands. So the natural claim is
  $$ \#\{\text{non-vanishing }(p+1)\text{-tuples}\} \;\le\; K\cdot C_{d,p}, \qquad C_{d,p}\le (2d)^p\ \text{(say)}, $$
  where $K$ is the number of summands. If one sticks with the coarse two-group split $H = H_E + H_O$ (as the Heisenberg derivation actually does), $K = 2$, **not $n$**, and the Principle's "$\cdot n$" does not come from the pivot vertex choice on the summand graph.
  The $\cdot n$ appears only when one further expands each $H_k = \sum_{b\in E} h_b$ into $\Theta(n)$ commuting bonds *inside* the summand — which is what the Heisenberg derivation actually does (line 98). The "graph" for that count is the graph of *bonds*, which has max degree $2$ on a 1D chain.
  The Principle as written elides this two-step structure (outer $K=\mathcal{O}(1)$ summands; inner $\Theta(n)$ commuting bonds within each), and a careful reader will be left unsure whether "$n$" in $n\cdot(2d)^p$ counts pivot summands or pivot bonds.
- **Concrete failure mode**: a reader tries to apply the Principle directly to the all-to-all TFIM (Childs Eq. 82), which has $K=2$ summands but $O(n^2)$ bonds inside each. The Principle predicts $\mathcal{O}(n)$ nested commutators; Childs proves $\mathcal{O}(n^3)$. The discrepancy is because "all-to-all" is *not* geometrically local — but the Principle as stated never names the true locality assumption (bounded *bond* degree, not bounded summand degree).
- **Side verification of the $(2d)^p$ bound on the concrete Heisenberg case**: $p=2$, $d=2$ (each bond has 2 site-sharing neighbours on each side), so $(2d)^p = 16$. The draft's own derivation (line 98) counts $2$ triples per pivot (per parity), times $2$ parities, giving $4$ connected triples per pivot bond, $\le 16$. So $(2d)^p$ is a valid but loose upper bound in this case. Good — but only because $d$ here is bond-degree, not summand-degree.
- **Why it matters**: the Principle is pitched as the umbrella that justifies all three specialisations (Heisenberg, TFIM, 2D Heisenberg). For that to read cleanly, it has to name *bond-level* locality, not the top-level summand decomposition. As written, it is slightly hand-wavy and will attract a reviewer's pen during defence.
- **Suggested fix**: refactor the Principle into a two-layer statement. Something like:
  > **Principle: geometrically-local commutator scaling.** A Hamiltonian is *geometrically local* if it can be written as $H = \sum_{b} h_b$ where each *local term* $h_b$ is supported on a bounded-size subset of sites and the *interaction hypergraph* on these terms has bounded maximum degree $d$ (each local term overlaps $\mathcal{O}(1)$ others). Grouping the local terms into $K = \mathcal{O}(1)$ commuting classes $H_k := \sum_{b\in\text{class}_k} h_b$ — each class being a sum of pairwise-commuting local terms — the nested commutator $[H_{k_{p+1}},\dots[H_{k_2},H_{k_1}]\dots]$ expands into a sum over ordered $(p+1)$-tuples of local terms $(b_1,\dots,b_{p+1})$, only *connected* tuples contributing. The number of such connected tuples anchored at a pivot $b_1$ is at most $(2d)^p$; summing over the $\Theta(n)$ pivots gives
  > $$ \tilde\alpha^{(p)}_{\mathrm{comm}} \;\le\; c_{d,p}\,n\cdot\bigl(\max_b\|h_b\|\bigr)^{p+1}, \qquad c_{d,p} \le (2d)^p\cdot\text{(ordering factor)}. $$
  > All specialisations below are instances of this two-layer counting: outer commuting-class grouping, inner bounded-degree bond counting.

  This makes plain that the $\cdot n$ is the pivot-bond count (not the pivot-summand count) and that the Principle's "bounded degree" is at the local-term level (so all-to-all TFIM is excluded, as it should be).
- **Lighter alternative (one-sentence patch)**: keep the current Principle but insert one parenthetical: after "… each site participates in $\mathcal{O}(1)$ summands" add "*— equivalently, each local bond overlaps $\mathcal{O}(d)$ others, which is what actually drives the $(2d)^p$ count below*". This prevents the misreading without a full refactor.

---

## [MINOR] m1. Unitarity sentence does not explicitly require Hermitian $H_k$

- **Location**: draft line 35, new sentence: "Each factor $\ee^{-\ii t H_k}$ is unitary, so $S_p(t)$ and its iterate $S_p^{(M)}(t) := S_p(t/M)^M$ are unitary for every choice of $p$ and $M$".
- **Issue**: the statement $\ee^{-\ii t H_k}$ is unitary **requires $H_k$ Hermitian**. The draft writes $H$ as "a sum of $K$ non-commuting summands" (line 16) without stating Hermiticity of each summand. It is implicit (because $H$ is Hermitian and each $H_k$ is "a single Pauli string or a commuting family of such strings," which are Hermitian by construction), but a careful reader would want this spelled out.
- **Fix**: prepend one short clause. Either
  > "Since each $H_k$ is Hermitian by construction, each factor $\ee^{-\ii t H_k}$ is unitary, so $S_p(t)$ and its iterate …"

  or, slightly earlier, rework line 21 to say "*non-commuting Hermitian summands*".
- **Placement**: the sentence is currently well-placed (after the Suzuki paragraph, before the "Commutator-scaling bound" section header), and the logical flow is fine.

## [MINOR] m2. Suzuki retention — the "why we keep it" is implicit, not stated

- **Location**: draft line 35: "Higher-order Suzuki formulas $S_{2k}$ can be built recursively from $S_2$ by the fractal construction … giving error $\mathcal{O}(t^{2k+1})$. The numerics of Part IV use $p=2$ throughout, so the displayed error analysis below is specialised to Strang."
- **Issue**: the draft tells the reader (i) that higher-order Suzuki exists and (ii) that Part IV doesn't use it. It does *not* explicitly say *why we bother to state the recursion at all*. The user's self-reassessment note implies it's for the Chapter 5 props, which are stated for general $p$. That's a legitimate reason, but only implicit in the current phrasing.
- **Fix (one-liner)**: insert between the two sentences:
  > "… giving error $\mathcal{O}(t^{2k+1})$. The Chapter 5 error analysis (`prop:trotter-diss`, `prop:trotter-B`) is stated for general $p$ and therefore carries the Suzuki recursion implicitly; the numerics of Part IV use $p=2$ throughout, so the displayed error analysis below is specialised to Strang."
- **Priority**: cosmetic. Do not block on this; merge-with-minor if the algebraic fix is applied.

## [MINOR] m3. "Other geometrically-local models" — the $\mathcal{O}(J^3 n)$ for 2D Heisenberg is stated without explicit bond-count justification

- **Location**: draft line 121, the collapsed closing paragraph: "…for 2D square-lattice Heisenberg (degree $d=4$) a four-colour edge split gives $\tilde\alpha^{(2)}_{\mathrm{comm}} = \mathcal{O}(J^3 n)$ with a larger constant prefactor."
- **Issue**: in 2D the bond-interaction graph has max degree $d = 6$ per bond (each bond touches $2\cdot 3 = 6$ other bonds at its two endpoints), not $d = 4$. The $d = 4$ is the *site* degree on the 2D square lattice, which is not the same as the bond-level $d$ that drives the Principle's $(2d)^p$ count.
- **Verdict**: the scaling $\mathcal{O}(J^3 n)$ is correct (still linear in $n$), so the conclusion stands. But quoting "degree $d = 4$" while using the bond-level Principle is inconsistent — same misalignment as M1.
- **Fix**: replace "(degree $d = 4$)" with either "(site degree $4$; bond degree $6$)" or just drop the parenthetical entirely and say "the same counting applies: bonds split into four commuting colour classes, pivot-bond count $\Theta(n)$, connected-triple count $\mathcal{O}(1)$ per pivot".

## [MINOR] m4. Sign-convention cut left one residual `\pm`

- **Location**: draft line 12, first paragraph: "the unitary Hamiltonian evolution $\ee^{\pm \ii H t}$ itself, used inside the controlled cascade of (Fig. circ:qpe) — each controlled-$U^{2^k}$ is implemented via Trotter steps of $S_p(t_0)$."
- **Issue**: the rest of the draft uses $\ee^{-\ii H t}$ exclusively (lines 24, 30, 41, 51, 59, 108). The $\pm$ in the motivational first sentence is fine *once* (it advertises that forward and backward evolutions both appear in the dissipator), but the draft's own Revision note 1 says "the draft silently uses $\ee^{-\ii Ht}$". That's not strictly true — there's still one $\pm$ on line 12.
- **Fix (pick one)**:
  - (a) change line 12's $\ee^{\pm \ii H t}$ to $\ee^{-\ii H t}$ and say "we show only the forward sign; the backward case is identical by $t \to -t$";
  - (b) leave line 12 as is and amend Revision note 1 to "the body uses $\ee^{-\ii Ht}$ throughout; the motivating first sentence retains $\ee^{\pm \ii Ht}$ to signal that the Lindbladian generator uses both signs."
- **Priority**: truly cosmetic — this one does not block.

---

## Orphan-reference sweep (Item 4 of the request)

Explicit verification:

| Target | Status |
|---|---|
| Removed sign-convention paragraph | No residual reference in body. ✓ |
| Removed "gain over blind triangle" meta-paragraph | No residual reference in body. The equivalent $\mathcal{O}(n^3)\to\mathcal{O}(n)$ comparison is now *inlined* in the Consequence section (line 131: "bounding each nested commutator by $4\|H_k\|^3$ with $\|H_E\|,\|H_O\|=\mathcal{O}(|J|n)$ and summing $K^2=4$ ordering terms"), so the baseline is self-contained. ✓ |
| Removed `eq:ising-alpha2` | `grep` confirms no reference to `eq:ising-alpha2` anywhere in the body (only in the Revision-notes block, which is to be deleted). ✓ |
| Removed `eq:lie-trotter-prefactor` | No reference anywhere in the body. ✓ |
| Removed 2D Heisenberg stand-alone paragraph | Collapsed into the final "Other geometrically-local models" paragraph; no dangling reference. ✓ |
| New `eq:trotter-M-count-heisenberg` label | Currently a **standalone display** (the display is self-referenced only by its own defining equation). There is no `(Eq. eq:trotter-M-count-heisenberg)` back-reference in the body. If the numerics chapter later cites the $M$-count by this label, good; if not, the label is harmless but unused. Acceptable. ✓ |
| `eq:heisenberg-evenodd-split` label | Also effectively standalone — not back-referenced. Same verdict: harmless. ✓ |

All orphan checks pass.

---

## Citation-table pruning sweep (Item 7 of the request)

Checked the body against the citation table (draft lines 149–156). The body cites:

- Thm 6 (line 39, `eq:childs-comm-bound` and surrounding text)
- Cor. 7 (line 59, Trotter-number bound)
- Prop. 9 Eq. 120 (line 49, "triangle-inequality bootstrap of Prop. 9")
- Prop. 10 Eq. 121 (lines 49, 112, "Strang analysis of Prop. 10 … a direct specialisation of Prop. 10")
- Eq. (122) (line 76, OBC Heisenberg)
- Fig. 3 (line 112, "$p=4$ version…")

Those are exactly the 6 rows in the table. I searched for `Sec. IV C`, `App. K`, `Eq. (82)`, `Eq. (90)`, `Eq. (K2)`, `Prop. M.1`, `Prop. M.2`, `Eq. M13`, `Eq. M14`, `Eq. (123)`, `Eq. (125)`, `Table II`, and `Fig. 3 inset` — none appear in the body. The pruning is clean. ✓

---

## Items confirmed accurate (preserve during revision)

- **$\sqrt{n}$ vs $n^{3/2}$ framing (line 137)**: correct. Plugging the generic $\tilde\alpha^{(2)}_{\mathrm{comm}} = \mathcal{O}(n^3)$ into Cor. 7 gives $M = \mathcal{O}(\sqrt{n^3}) = \mathcal{O}(n^{3/2})$; the structural $\tilde\alpha^{(2)}_{\mathrm{comm}} = \mathcal{O}(n)$ gives $M = \mathcal{O}(\sqrt{n})$. The quadratic sharpening claim is honest.
- **$\mathcal{O}(n^3)$ baseline derivation (line 131)**: now in-line and explicit — "bounding each nested commutator by $4\|H_k\|^3$ with $\|H_E\|,\|H_O\| = \mathcal{O}(|J|n)$ and summing $K^2 = 4$ ordering terms". This fully absorbs the old "gain over blind triangle" meta-paragraph. ✓
- **Unitarity sentence's placement** (line 35, after Suzuki paragraph, before "Commutator-scaling bound"): flows well. It sets up the Parseval remark referenced in the Consequence section without interrupting the main narrative.
- **Heisenberg constant $54|J|^3 n$ and the breakdown $54 = 4\cdot 108\cdot(1/12+1/24)$**: arithmetic re-verified ($108/12 + 108/24 = 9 + 4.5 = 13.5$; then $4 \cdot 13.5 = 54$). ✓
- **Prop. 9 inline acknowledgement** (line 49): the "triangle-inequality bootstrap of Prop. 9 (Eq. 120) gives the analogous tight first-order prefactor" reads cleanly and cedes the $p = 1$ analogue to the paper without a separate display. Good choice.
- **No new citations introduced** — verified: the body cites only `childs2021theory`, which is already in `references.bib`. The author's claim in the header holds. ✓

---

## Summary Scorecard

| Category | Rating | Notes |
|----------|--------|-------|
| Logical correctness | **MAJOR** | `eq:trotter-M-count-heisenberg` explicit constant off by factor of 2 (B1); Principle $(2d)^p$ count conflates summand vs bond graph (M1). Scaling exponents ($\sqrt n$, $n^3\to n$) all correct. |
| Notation consistency | OK | $K$ vs $\Gamma$ harmonised; $\ee^{-\ii Ht}$ uniform except one motivational $\pm$ in line 12 (m4, cosmetic). |
| Citations | OK | Table pruned cleanly; body cites only Thm 6 / Cor 7 / Prop 9 / Prop 10 / Eq. 122 / Fig. 3, all verified against the PDF. No new bibtex entries. |
| Thesis consistency | OK | Strang bound `prop:trotter-diss:929` picks up the new `eq:heisenberg-alpha2` substitution cleanly; the only issue is the arithmetic slip that will propagate downstream if not caught (B1). |
| Exposition quality | MINOR | Unitarity sentence needs Hermiticity qualifier (m1); Suzuki retention rationale is implicit (m2); 2D Heisenberg degree typo (m3). All small. |

---

## Priority-ordered fix list

**Must fix before merging (blocker + major):**

1. **B1**: replace $\sqrt{15\cdot 54}$ with $15^{1/4}\sqrt{54}$ (or rewrite as $\sqrt{\sqrt{15}\cdot 54\,|J|^3 n /(\sigma^3\varepsilon)}$) in `eq:trotter-M-count-heisenberg`.
2. **M1**: refactor the *Principle* to name *bond-level* locality (not summand-level). Either the full two-layer rewrite above or the one-sentence parenthetical patch.

**Should fix (minor):**

3. **m1**: add "Since each $H_k$ is Hermitian" to the unitarity sentence (one clause).
4. **m3**: change the "(degree $d = 4$)" in the 2D Heisenberg aside to "(site degree 4; bond degree 6)" or drop the degree parenthetical entirely.

**Optional polish (nit):**

5. **m2**: add one half-sentence explaining why Suzuki is retained.
6. **m4**: either drop the stray `\pm` on line 12 or amend the Revision note 1 to acknowledge it.

---

## Verdict

**One more revision needed.**

The algebraic slip in `eq:trotter-M-count-heisenberg` is a [BLOCKER] by any reasonable standard — an off-by-two in the one explicit numerical constant the subsection provides. That single character-level edit, plus the one-sentence Principle refactor (M1), gets the draft to merge-ready.

Turnaround is low-effort: changing $\sqrt{15\cdot 54}$ to $15^{1/4}\sqrt{54}$ is a one-character fix in the display plus a mirror fix in Revision note 7. The Principle refactor is three–four sentences. With those two edits in v4, this subsection is clean to copy into `1_preliminaries.tex`.
