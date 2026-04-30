# Draft Review: Weak-measurement subsection (v4) — focus on the new "Choice of Lindblad simulator" section

**Reviewed**: `drafts/weak-measurement-subsection-v4.md` (v3 → v4 delta is the new section spanning lines 43–51, plus updated header line 7 and bibtex stubs lines 71–95).
**Date**: 2026-04-27
**Overall assessment**: The new section makes the right architectural point and the source-paper attributions are essentially correct, but the body sentence misstates Li-Wang Theorem 1's polylog argument (`polylog(1/ε)` vs. the correct `polylog(t/ε)`), uses bare `t` in place of `τ`, and introduces a notation inconsistency between the new section's `$\hat{A}_a(\omega)$` and the rest of the v4 draft / methods chapter's `$\tilde{A}_a(\bar\omega)`. None of the issues are conceptual; all are local fixes. Verdict at the bottom.

---

## Fatal

### F1. Body sentence misstates Li-Wang Theorem 1's polylog argument and drops $\|\mathcal{L}\|_{\mathrm{be}}$
- **Location**: line 45, body of the new section. The sentence reads:

  > *"$\widetilde{\mathcal{O}}(t \,\mathrm{polylog}(1/\varepsilon))$ block-encoding queries (Theorem 1, with $\tau := t \|\mathcal{L}\|_{\mathrm{be}}$)"*

- **Issue**: Li-Wang Theorem 1 (informal, p. 4 of arXiv:2212.02051v2) says, verbatim, the algorithm uses $\widetilde{O}(\tau\,\mathrm{polylog}(t/\varepsilon))$ queries to $U_H$ and $U_{L_j}$, with $\tau \mathrel{:=} t\|\mathcal{L}\|_{\mathrm{be}}$ and $\|\mathcal{L}\|_{\mathrm{be}} \mathrel{:=} \alpha_0 + \tfrac{1}{2}\sum_{j=1}^m \alpha_j^2$ (their Eq. 4). The v4 body has *two* misstatements packed in one line: (a) the polylog argument is $t/\varepsilon$, not $1/\varepsilon$; (b) the leading factor should be $\tau$, not $t$. The parenthetical immediately after defines $\tau$ but then never uses it — the reader is told $\tau := t\|\mathcal{L}\|_{\mathrm{be}}$ and *also* sees the cost written with bare $t$, which is internally inconsistent.
- **Why it matters**: this is the only quantitative claim about Li-Wang in the section, and it is what the architectural argument turns on (a polylog-in-precision algorithm). Stating it sloppily undermines the very point the section is making, and any reader who looks up the Li-Wang paper will find a precise mismatch.
- **Suggested fix**: replace the offending clause with

  > "$\widetilde{\mathcal{O}}(\tau\,\mathrm{polylog}(t/\varepsilon))$ block-encoding queries to $U_H$ and to the $U_{L_j}$, where $\tau := t\,\|\mathcal{L}\|_{\mathrm{be}}$ (Theorem 1)"

  Note the BibTex stub at line 81 already has the correct `polylog(t/ε)`, so this is just a body-sentence sync.

### F2. The body sentence omits the linear-$m$ qualifier, which is the load-bearing part of the contrast
- **Location**: line 45–49, body of the new section.
- **Issue**: the section's argument is "Li-Wang gives polylog-in-$\varepsilon$ but cost is *linear in $m$*, the number of jump operators; weak-measurement gives polynomial in $1/\varepsilon$ but absorbs the OFT integral into a *single* $U_L$." The linear-$m$ part is the hinge of the argument and it is what makes the comparison non-trivial. But the body sentence quotes only $\widetilde{\mathcal{O}}(\tau\,\mathrm{polylog}(t/\varepsilon))$ block-encoding *queries* and never mentions Li-Wang's $\widetilde{\mathcal{O}}(m\tau\,\mathrm{polylog}(t/\varepsilon))$ additional 1- and 2-qubit gates. The next sentence ("their cost is linear in $m$") then asserts linear-in-$m$ scaling that the previous sentence's complexity expression does not display.
- **Why it matters**: the reader sees a polylog cost that has no $m$ in it and is then told the cost is linear in $m$ — these look contradictory unless the reader already knows that the queries-vs-gates split is where the $m$ lives. Without the $m$ visible somewhere in a cost expression, the architectural argument has no quantitative anchor.
- **Suggested fix**: after the cost statement add a parenthetical:

  > "(plus $\widetilde{\mathcal{O}}(m\tau\,\mathrm{polylog}(t/\varepsilon))$ additional 1- and 2-qubit gates, which is what makes the total cost linear in $m$)"

  This matches the precise breakdown given in Li-Wang Theorem 1 and motivates the next sentence.

---

## Important

### I1. Notation collision: `$\hat A_a(\omega)$` (new section) vs. `$\tilde A_a(\bar\omega)$` (everywhere else in the same v4 file and in the methods chapter)
- **Location**: new section line 47 uses `$\hat A_a(\omega)$`, but line 49 also uses `$\hat A_a(\omega)$`; the §"Generality of the control structure" at line 55 uses `$\tilde A_a(\bar\omega)$`. The methods chapter uses `$\hat A^a(\omega)$` (note: superscript, not subscript) in its continuous OFT-Lindbladian display at `2_methods.tex:64, 109, 210` and switches to `$\tilde A_a(\bar\omega)$` (subscript, tilde, discretized) in the algorithm at `2_methods.tex:1778, 1887, 1893, 1895`.
- **Issue**: there are now *three* OFT notations floating in the immediate vicinity: `$\hat A_a(\omega)$` (v4 new section, subscript + hat + continuous $\omega$), `$\tilde A_a(\bar\omega)$` (v4 generality section + methods algorithm, subscript + tilde + discretized $\bar\omega$), `$\hat A^a(\omega)$` (methods CKBG/CKG Lindbladians, *superscript* + hat + continuous $\omega$). The new section's choice is the closest analogue to the methods CKBG form except that it has a subscript $a$ instead of superscript $a$. The note at draft line 110 ("$L$ versus $L_a$") commits the prelims subsection to subscript-$a$ throughout, so the new section is internally consistent with v3's `$L_a$` choice but breaks symmetry with how `$\hat A^a$` is written in the methods continuous form.
- **Why it matters**: the section sits one chapter before the methods' CKG Lindbladian and is intended to introduce the OFT Lindbladian shape to a reader who has not yet seen `eq:L-ckbg`. If they later see `$\hat A^a(\omega)$` in the methods, they have to do mental gymnastics to recognize it as the same object. Worse, the v4 generality section talks about "the Euler-step identity persists with $L$ replaced by the OFT-filtered jump $\sqrt{\gamma(\bar\omega)}\,\tilde A_a(\bar\omega)$" — so within v4 *itself* the OFT object is named two different ways (hat + $\omega$ in the new section, tilde + $\bar\omega$ in the generality section), with no remark explaining that they are the continuous and discretized versions of the same map.
- **Suggested fix**: pick one of the following:
  1. (Cheapest) Add one inline phrase to the new section's display equation: "with the continuously $\omega$-parametrized OFT $\hat A_a(\omega)$, of which the discretized $\tilde A_a(\bar\omega)$ used in §sec:algorithm is the quadrature on the $\Omega$-register grid". This costs one clause and resolves the apparent collision.
  2. (More invasive) Drop $\hat A_a$ entirely from the new section and write the display as

      $$\mathcal{L}(\rho) \;=\; \sum_{a \in \mathcal{A}} \int \dd\omega\;\gamma(\omega)\,\mathcal{D}_{\hat A^a(\omega)}(\rho),$$

      matching the methods superscript-$a$ convention exactly. Then the line 49 textual reference also becomes "$\hat A^a(\omega)$" and the only remaining v4-internal mismatch is hat-vs-tilde (continuous vs discretized), which fix (1) handles.

### I2. The label `\label{sec:prelim-QPE}` does not exist in `1_preliminaries.tex` — references will dangle
- **Location**: new section line 49 ("addressed coherently through the $\Omega$ register built up in §sec:prelim-QPE"); line 51 ("the OFT machinery of §sec:prelim-QPE is only useful if it is exercised"). Also references in v3's body that are now consumed by v4 unchanged: insertion-target header line 5, generality section line 55 (no `sec:prelim-QPE` ref there but the chain of cross-refs depends on it).
- **Issue**: I checked `1_preliminaries.tex` directly. The Phase Estimation subsection at line 911 has *no* `\label{sec:prelim-QPE}` — only `\label{circ:qpe}` for the figure (line 907). The `\label{sec:prelim-LCU}` is at line 1079 (so that ref *does* resolve). The Trotter subsection at line 968 has `\label{sec:prelim-trotter}`. So `\ref{sec:prelim-QPE}` will produce `??` at compile time.
- **Why it matters**: the new section uses the `sec:prelim-QPE` label twice and the body of v3 also uses it; once the .md is translated to .tex, both refs will be broken. This is a v3 carry-over but the new section makes the dangling reference more conspicuous (it now appears in the rhetorically loaded closing sentence at line 51).
- **Suggested fix**: ask the user to add `\label{sec:prelim-QPE}` to the `\subsection{Phase Estimation.}` line at `1_preliminaries.tex:911`, e.g. change

  ```latex
  \subsection{Phase Estimation.} One of the key subroutines ...
  ```

  to

  ```latex
  \subsection{Phase Estimation.} \label{sec:prelim-QPE} One of the key subroutines ...
  ```

  This is a one-line patch and is also implicitly required by every prior v3 reference. Flag it in the insertion-target header at line 5 of the draft (where labels-already-in-the-chapter are listed) — currently line 5 lists `sec:prelim-QPE` as if it exists.

### I3. Constraint compliance: the `$\gamma(\omega)$` "Kossakowski function" framing is borderline pre-supposing KMS DB
- **Location**: line 47 (display) and line 51 ("replacing $\gamma(\omega)$ with a discrete sum of $\delta$-functions"). The section deliberately stays name-free on Chen 2025 KMS / Ding-Li-Lin in the body, which honours the user's hard constraint. The brief context note in the draft (line 112) explicitly says "the new section deliberately stays in fully general Lindblad-simulation terms". So far so good.
- **Issue**: but the framing of $\gamma(\omega)$ as a continuously parametrized weight on an OFT-decomposed jump set is *almost a verbatim* setup for the KMS-DB Lindbladian of the methods chapter (compare to `eq:L-ckbg` and `eq:ckg-L`); a reader who has not yet seen the KMS construction will read this as "the Lindbladian *I will be working with*" rather than "a generic Lindbladian shape". The closing sentence at line 51 — "Since the construction the rest of this thesis follows keeps $\gamma(\omega)$ continuous (Chapter 5; the OFT machinery of §sec:prelim-QPE is only useful if it is exercised)" — names "Chapter 5" and explicitly forward-points to the rest of the thesis, which softens the "general" framing intentionally. Whether this is OK depends on how strict the user wants the constraint to be.
- **Why it matters**: the user's hard constraint as restated in the brief is "no naming of Chen 2025 KMS or Ding et al in the body". Naming "Chapter 5" is not a violation but it does break the "fully general" pretence and the reader will realise the section is implicitly motivating the choice of weak-measurement for the rest of the thesis. If the user wants a self-contained technical observation (the brief says "self-contained"), the closing should not name the chapter.
- **Suggested fix**: two options:
  1. (Lighter) Change the parenthetical "(Chapter 5; the OFT machinery of §sec:prelim-QPE is only useful if it is exercised)" to "(see the chapter on the dissipative Gibbs sampler; the OFT machinery of §sec:prelim-QPE is only useful if it is exercised)" — keeps the forward pointer but doesn't lock to a specific chapter number.
  2. (Stricter — preferred for self-contained reading) Drop the closing parenthetical entirely and end the section with: "Since the construction the rest of this thesis follows keeps $\gamma(\omega)$ continuous, the weak-measurement primitive is the only one of the two that applies." This removes both "Chapter 5" and the slightly cute "only useful if it is exercised" phrasing (see I4 below).

### I4. The closing sentence at line 51 is too cute and reads as a non-sequitur after a clean technical comparison
- **Location**: line 51 last sentence: "*Since the construction the rest of this thesis follows keeps $\gamma(\omega)$ continuous (Chapter 5; the OFT machinery of §sec:prelim-QPE is only useful if it is exercised), the weak-measurement primitive is the only one of the two that applies.*"
- **Issue**: the parenthetical "the OFT machinery of §sec:prelim-QPE is only useful if it is exercised" is a rhetorical flourish — it does not justify the conclusion (the conclusion is justified by the preceding "keeps $\gamma(\omega)$ continuous"), and reads as a slightly winking justification of a thesis-organisational choice rather than a technical observation. The brief flagged this exact concern, and I agree.
- **Why it matters**: the rest of the section is sober declarative prose; this sentence breaks register and the reader notices.
- **Suggested fix**: as above (I3 fix 2): drop the parenthetical and the final clause becomes a clean technical conclusion: "Since the construction the rest of this thesis follows keeps $\gamma(\omega)$ continuous, the weak-measurement primitive is the only one of the two that applies." Five words shorter and the architectural punchline is preserved.

### I5. `[CITE: ding2024efficient]` is on the wrong sentence — it currently points at the *workaround* but the cited paper *is* the workaround
- **Location**: line 51, mid-sentence: "*typical workarounds for the continuous-$\omega$ case — e.g. replacing $\gamma(\omega)$ with a discrete sum of $\delta$-functions so that the integral becomes a finite sum [CITE: ding2024efficient] — give back exactly the finite-$m$ regime in which the higher-order route applies, at the cost of working with a different jump ensemble.*"
- **Issue**: the citation is correctly attached to the technique (discrete-$\delta$ Kossakowski). I verified Ding-Li-Lin p. 3, first paragraph: *"In particular, $\gamma(\omega)$ can be chosen to be a discrete sum of $\delta$ functions, leading to a finite number of jump operators ... can be efficiently simulated using any high-order simulation algorithms, including those in [LW23, DLL24]."* So the citation is correctly placed. **Not a fatal issue, just a phrasing point**: calling Ding-Li-Lin's construction a "workaround" is mildly pejorative — their paper *is* the family of KMS-DB Lindbladians built around the discrete-$\delta$ Kossakowski, and the discrete-jump option is one of its central design choices, not a workaround. The brief's draft of the contribution paragraph (line 155) used the more neutral "the natural way to enter that finite-$m$ regime from a KMS-DB Lindbladian", which is better.
- **Suggested fix**: change "*typical workarounds for the continuous-$\omega$ case — e.g. replacing $\gamma(\omega)$ with a discrete sum of $\delta$-functions*" to "*reformulating the KMS-DB Lindbladian so that $\gamma(\omega)$ is a discrete sum of $\delta$-functions and the integral becomes a finite sum [CITE: ding2024efficient]*". This is more accurate and doesn't editorialise.

### I6. The statement "weak-measurement Euler step ... cost ... in $1/\varepsilon$ — polynomial, not poly-logarithmic" is correct in spirit but the use of $\widetilde{\mathcal{O}}$ is unmotivated
- **Location**: line 45, opening of the body: "*$T/\delta = \widetilde{\mathcal{O}}(T^2/\varepsilon)$ calls to $U_L$ — polynomial, not poly-logarithmic, in $1/\varepsilon$.*"
- **Issue**: the standard derivation is: per-step error is $\mathcal{O}(\delta^2)$ from `eq:weak-meas-euler`; over $T/\delta$ steps the cumulative diamond-norm error is $\mathcal{O}(T\delta)$; setting this $\leq \varepsilon$ gives $\delta = \mathcal{O}(\varepsilon/T)$, so $T/\delta = \mathcal{O}(T^2/\varepsilon)$. There is no polylog factor anywhere — the bare $\mathcal{O}$ is the right symbol, not $\widetilde{\mathcal{O}}$. The tilde absorbs subdominant polylog factors, and there are none here.
- **Why it matters**: the reader will pause to wonder what hidden log factors are being absorbed — there aren't any. The line is also drawing a contrast with Li-Wang's $\widetilde{\mathcal{O}}(\tau\,\mathrm{polylog}(t/\varepsilon))$ where the tilde is doing real work. Using the same symbol for both elides the qualitative difference.
- **Suggested fix**: change "$T/\delta = \widetilde{\mathcal{O}}(T^2/\varepsilon)$" to "$T/\delta = \mathcal{O}(T^2/\varepsilon)$". Two characters; sharpens the polynomial-vs-polylog contrast that is the section's whole point.

---

## Suggestions

### S1. Tighten the bibtex stubs to match `references.bib` style
- **Location**: lines 71–95 (bibtex block).
- **Verified**: the arXiv IDs (2212.02051 and 2404.05998) and authorship are correct. Li & Wang's full title in the actual paper / arXiv is "Simulating Markovian open quantum systems using higher-order series expansion" (note the hyphen in "higher-order"); the v4 stub omits the hyphen. The Ding-Li-Lin paper title in the PDF first page is "EFFICIENT QUANTUM GIBBS SAMPLERS WITH KUBO–MARTIN–SCHWINGER DETAILED BALANCE CONDITION" with regular en-dashes between Kubo, Martin, Schwinger; v4 already uses double-hyphens which is biblatex-canonical. Both arXiv IDs check out: 2212.02051 (Li-Wang, Dec 2022) and 2404.05998 (Ding-Li-Lin, Apr 2024).
- **Suggestion**: change `Simulating {M}arkovian open quantum systems using higher order series expansion` to `Simulating {M}arkovian open quantum systems using higher-order series expansion` (insert hyphen between "higher" and "order" to match the canonical title). For Ding-Li-Lin the v4 stub is already correct.

### S2. Two paragraphs in the new section: paragraph 2 does add a distinct point, but the "trade-off is architectural rather than numerical" framing could be sharpened
- **Location**: lines 51 (paragraph 2 of the new section).
- **Observation**: the brief asks whether paragraph 2 repeats paragraph 1 or adds a distinct point. My read is that it *does* add a distinct point: paragraph 1 establishes the technical fact (continuous-$\omega$ requires single-block-encoded $L$ with sub-OFT-grid resolution); paragraph 2 reframes this as a pre-compile-time vs run-time choice (Li-Wang requires the jump list to be collapsed *at compile time*; weak-measurement defers it). But the connecting phrase "trade-off is therefore architectural rather than numerical" doesn't unpack what "architectural" means — the rest of the paragraph then does unpack it implicitly. A tighter version would say "the choice is between collapsing the jump list at compile time (Li-Wang) and absorbing it coherently into a single $U_L$ at run time (weak-measurement)".
- **Suggestion**: reword "*The trade-off between the two routes is therefore architectural rather than numerical.*" to "*The trade-off between the two routes is therefore one of compile-time vs run-time jump-set resolution.*" — same length, more concrete.

### S3. The cost-comparison table mentioned in Suggestion 6 (line 145) really would help a future reader, and is not "premature" given the new section's load
- **Location**: line 145 (Suggestions in writing notes).
- **Observation**: with the new section, three quantities are now relevant in the body: (a) Li-Wang query count $\widetilde{\mathcal{O}}(\tau\,\mathrm{polylog}(t/\varepsilon))$; (b) Li-Wang gate count $\widetilde{\mathcal{O}}(m\tau\,\mathrm{polylog}(t/\varepsilon))$; (c) weak-measurement step count $\mathcal{O}(T^2/\varepsilon)$. Asking the reader to keep three quantities in their head while parsing two paragraphs of prose is a real ask. A 2 × 3 table (Euler vs Li-Wang × queries / gates / scaling-with-jumps) would be ~6 lines of LaTeX and would crystallize the architectural argument that the prose is making.
- **Suggestion**: re-evaluate Suggestion 6 (line 145) — the table is now a real win, not just a future-growth contingency. If the user prefers prose, leave as-is; if they will accept a small table, adding one would compress the two-paragraph section into one short paragraph + one table.

### S4. The display equation at line 47 should pick up an HTML-comment label so the surrounding text can `\eqref{}` it later
- **Location**: line 47, the unlabelled display $\mathcal{L}(\rho) = \sum_a \int \dd\omega \gamma(\omega) \mathcal{D}_{\hat A_a(\omega)}(\rho)$.
- **Observation**: this is the *only* unlabelled display equation in the new section; all of v3's displays carry inline label hints (or are explicitly noted as unlabelled). The note at line 116 says one display is labelled (`eq:weak-meas-euler`) and "one is the unlabelled $\sum_a \int \dd\omega$ Lindbladian display in the new section". Currently nothing in the body refers back to it, so the choice not to label is defensible. But if the v4 → v5 process refines the section and the body needs to point back to "the OFT-Lindbladian shape we considered", the absence of a label will be a small friction.
- **Suggestion**: add an inline label hint: change line 47's display block to

  $$\mathcal{L}(\rho) \;=\; \sum_{a \in \mathcal{A}} \int \dd\omega\;\gamma(\omega)\,\mathcal{D}_{\hat A_a(\omega)}(\rho),$$ <!-- \label{eq:oft-lindbladian-shape} -->

  with name `eq:oft-lindbladian-shape` (or similar). This costs nothing now and gives later revisions a hook.

### S5. Particularly well-done: keep these on revision
- The section's *strategic* placement between §"One Euler step" and §"Generality of the control structure" is correct: the reader has just absorbed the Euler-step cost (line 41), and now naturally asks "but couldn't I do better with a higher-order simulator?" — the new section answers exactly this.
- The architectural framing (collapse-at-compile-time vs absorb-into-single-$U_L$-at-run-time) is *the* right framing and is not in the literature in this compact form. This is a real expository contribution.
- The hard constraint to stay name-free in the body (no Chen 2025 KMS, no Ding-Li-Lin) is honoured (modulo I3 above), and the citation pattern follows v3 cleanly.
- The bibtex stubs are *almost* perfect (just S1's hyphen) and the BibTex stub for Li-Wang correctly carries the precise polylog($t/\varepsilon$) form even though the body misstates it (F1).

---

## Summary Scorecard

| Category | Rating | Notes |
|---|---|---|
| Logical correctness | MAJOR | F1 (polylog argument is wrong) and F2 (linear-$m$ qualifier missing) are quantitative bugs that the architectural argument turns on; the rest of the section's logic is correct. |
| Notation consistency | MAJOR | I1 — the new section's `$\hat A_a(\omega)$` collides with the v4 generality section's `$\tilde A_a(\bar\omega)$` and the methods chapter's `$\hat A^a(\omega)$`. Needs at minimum a one-line bridge. |
| Citations | MINOR | All three citations (chen2023quantum, chen2023efficient, liwang2023, ding2024efficient) verified against the source PDFs / arXiv pages; the BibTex stubs are correct modulo S1's missing hyphen in Li-Wang title. The `[CITE: ding2024efficient]` placement is correct (I5 is phrasing only). |
| Thesis consistency | MAJOR | I2 — `\label{sec:prelim-QPE}` does not exist in `1_preliminaries.tex`. The new section's two cross-refs to it will dangle, as will v3's reference. Needs a one-line patch to `1_preliminaries.tex:911`. |
| Exposition quality | MINOR | I3 (constraint compliance softness via "Chapter 5"), I4 (cute closing parenthetical), I6 ($\widetilde{\mathcal{O}}$ overuse) are all small. Paragraph structure is sound (S2 confirms paragraph 2 adds a distinct point). |

---

## Verdict

**Not ready for insertion** as currently written, but the gap to ready-for-insertion is small.

Required fixes before insertion (all are local, no new technical content needed):
1. Fix F1 (correct the polylog argument and use $\tau$ not $t$) — line 45.
2. Fix F2 (add the $\widetilde{\mathcal{O}}(m\tau\,\mathrm{polylog}(t/\varepsilon))$ gate-count parenthetical to make the linear-$m$ claim visible) — line 45.
3. Fix I1 (one-line bridge between `$\hat A_a(\omega)$` and `$\tilde A_a(\bar\omega)$`) — line 47 or 49.
4. Fix I2 (ask the user to add `\label{sec:prelim-QPE}` to `1_preliminaries.tex:911`; flag in the insertion-target header at draft line 5).

Recommended fixes:
5. I3 + I4 (clean up the closing sentence — drop "(Chapter 5; ... only useful if it is exercised)").
6. I5 (re-word "workaround" — change to "reformulating the KMS-DB Lindbladian so that $\gamma(\omega)$ is a discrete sum of $\delta$-functions").
7. I6 ($\mathcal{O}$ not $\widetilde{\mathcal{O}}$ for the Euler-step count).
8. S1 (insert hyphen in Li-Wang bibtex title).

After these eight one-line patches, the new section is ready. The rest of v4 (the v3-inherited body) is unaffected by this review per the brief's instruction.
