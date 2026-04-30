# Weak-measurement (v2)

> **Insertion target:** `supplementary-informations/1_preliminaries.tex`, body of `\subsection{Weak-measurement.} \label{sec:prelim-weak-meas}` — text goes between the `\subsection{...}\label{...}` line at line 1116 and `\begin{figure}[ht]` at line 1117 (the figure `circ:weak-measurement` sits at lines 1117–1167).
> **Depends on labels already in the chapter:** `sec:prelim-LCU`, `sec:prelim-QPE`, `eq:kraus` (Kraus form of CPTP maps at `1_preliminaries.tex:47`-ish), `circ:weak-measurement`, and the $Y_\theta := R_Y(2\arcsin\sqrt\theta)$ shorthand at `1_preliminaries.tex:825`.
> **Forward reference:** `alg:diss` (the `DissipativeStep` algorithm) in Chapter 5 at `2_methods.tex:1744`.
> **Bib keys used, already in `references.bib`:** `chen2023quantum`, `chen2023efficient`. No new keys introduced.

---

A Lindblad dissipator $\mathcal{D}_L(\rho) := L\rho L^\dagger - \tfrac{1}{2}\{L^\dagger L,\rho\}$ combines a *transition* term $L\rho L^\dagger$ with an anti-Hermitian *drain* $-\tfrac12\{L^\dagger L,\rho\}$ whose interplay is what makes $\rho \mapsto \ee^{t\mathcal{D}_L}\rho$ trace-preserving. The LCU of §sec:prelim-LCU delivers the transition term: its post-selection on $\ket{0^b}$ of the block-encoded $U_L$ imprints $L\ket\psi/\alpha$ on $S$ (with $\|L/\alpha\|\leq 1$ by §sec:prelim-LCU). What LCU alone does *not* deliver is the drain. The coherent-undo construction of [Chen et al. 2023] [CITE: chen2023quantum] fills this gap by weaving a single weak-measurement ancilla $\ket{0}_\delta$ between $U_L$ and a conditional $U_L^\dagger$, in such a way that post-selecting $b$ on $\ket{0^b}$ and measuring $\ket{0}_\delta$ implements exactly one Euler step of $\ee^{\delta\mathcal{D}_L/\alpha^2}$ on $\rho$, up to $\mathcal{O}(\delta^2)$ (Fig. circ:weak-measurement). It is this primitive that reappears as the dissipative step of Chapter 5's Lindblad simulation (Alg. alg:diss). (The present $\mathcal{D}_L$ is distinct from the KMS discriminant $\mathcal{D}(\rho,\mathcal{L})$ of (Eq. eq:discriminant) in §sec:davies; the jump argument keeps the two apart.)

## Circuit walkthrough

Initialise the three registers in $\ket{0}_\delta \otimes \ket{0^b} \otimes \rho$ and apply the four operations of Fig. circ:weak-measurement in turn.

*(i) Block encoding.* $U_L$ acts on $(b,S)$; on the $\ket{0^b}$-block of the combined state it imprints the amplitude $L/\alpha$ on $S$, leaving an orthogonal component $\ket{\phi^\perp}_{bS} := (\idm - \ket{0^b}\!\bra{0^b}\!\otimes\idm)\,U_L\,\ket{0^b}\!\ket\psi$ outside the block.

*(ii) Ancilla rotation.* The $Y_\delta$ of the thesis's $Y_\theta := R_Y(2\arcsin\sqrt\theta)$ shorthand at line 825 is *applied conditional on $\ket{0^b}$*, i.e. it fires only on the post-selectable branch of (i):

$$Y_\delta\ket{0}_\delta\ket{0^b} \;=\; \bigl(\sqrt{1-\delta}\,\ket{0}_\delta + \sqrt{\delta}\,\ket{1}_\delta\bigr)\otimes\ket{0^b},$$

so the weak-measurement ancilla acquires amplitude $\sqrt{\delta}$ on $\ket{1}_\delta$ precisely where $U_L$ succeeded. The $\ket{\phi^\perp}_{bS}$ branch carries the original $\ket{0}_\delta$.

*(iii) Open-controlled undo.* The open-controlled $U_L^\dagger$ fires only on the $\ket{0}_\delta$ branch. On the $\ket{0^b}$-block of $U_L\ket{0^b}\ket\psi$ this uncomputes to $\ket{0^b}\ket\psi$; the leaked $\ket{\phi^\perp}_{bS}$ is *not* uncomputed and survives as an out-of-block residual. On the $\ket{1}_\delta$ branch $U_L^\dagger$ is absent, and $\ket{0^b}\otimes(L/\alpha)\ket\psi$ stands.

*(iv) Post-selection and measurement.* Post-select $b$ on $\ket{0^b}$ and measure $\ket{0}_\delta$.

## Kraus operators

The two outcomes of step (iv) give two Kraus operators on $S$ (in the sense of (Eq. eq:kraus)). In the *no-jump* branch ($\ket{0}_\delta$ outcome, $\ket{0^b}$ post-selected), the uncomputation $U_L^\dagger U_L = \idm$ acts on the block, while the $\ket{\phi^\perp}_{bS}$ residual has $\ket{0^b}$ post-selection failure probability $\|L\ket\psi\|^2/\alpha^2$; bookkeeping the trace loss reconstitutes the anti-Hermitian drain. In the *jump* branch ($\ket{1}_\delta$ outcome, $\ket{0^b}$ post-selected), the full amplitude $L\ket\psi/\alpha$ survives, weighted by the $\sqrt\delta$ of the rotation. The Chen et al. 2023 amplitude bookkeeping [CITE: chen2023quantum][CHEN23-EQ3.2] gives

$$M_0 \;=\; \sqrt{1-\delta}\,\idm \;-\; \tfrac{1}{2}\,\tfrac{\delta}{\alpha^2}\,L^\dagger L \;+\; \mathcal{O}(\delta^2), \qquad M_1 \;=\; \sqrt{\delta}\;\frac{L}{\alpha}.$$

## One Euler step

Summing the two branches and expanding in $\delta$,

$$\rho \;\mapsto\; M_0\,\rho\,M_0^\dagger + M_1\,\rho\,M_1^\dagger \;=\; \rho \;+\; \frac{\delta}{\alpha^{2}}\Bigl(L\rho L^\dagger - \tfrac{1}{2}\{L^\dagger L,\rho\}\Bigr) \;+\; \mathcal{O}(\delta^{2}),$$ <!-- \label{eq:weak-meas-euler} -->

which is one Euler step of $\dot\rho = (1/\alpha^2)\,\mathcal{D}_L(\rho)$ — the rstick of Fig. circ:weak-measurement. Identifying $\delta = \gamma_{\mathrm{eff}}\,\dd t$ with $\gamma_{\mathrm{eff}} := 1/\alpha^2$ the effective Lindblad rate of the block-encoded generator (per unit time), the post-selected channel implements one Euler step of $\dot\rho = \gamma_{\mathrm{eff}}\,\mathcal{D}_L(\rho)$. The post-selection on $b$ succeeds with probability $1 - \delta\,\|L\ket\psi\|^2/\alpha^2 + \mathcal{O}(\delta^2) = 1 - \mathcal{O}(\delta)$, so — unlike the $1/\alpha^2$ post-selection loss of §sec:prelim-LCU — no oblivious amplitude amplification is needed: $\delta$-smallness replaces amplification as the deterministic-acceptance mechanism. In effect, the coherent-undo construction turns the LCU's $1-1/\alpha^2$ block-failure amplitude into signal: the leaked projection becomes the Lindblad drain. The extension to the multi-jump sum $\mathcal{L} = \sum_a \mathcal{D}_{L_a}$ and to the KMS-weighted filtered jumps of [Chen et al. 2025] [CITE: chen2023efficient] that replace $L$ in Chapter 5 is the dissipative step of (Alg. alg:diss).

---

## Citations

Both keys already exist in `references.bib` — no new bibtex entries are introduced.

- **`chen2023quantum`** — Chen, C.-F., Kastoryano, M. J., Brandão, F. G. S. L., Gilyén, A., *"Quantum Thermal State Preparation"*, arXiv:2303.18224 (2023). Theorem III.1 and Eq. (3.2) in their proof establish the weak-measurement implementation of a Lindblad dissipator from a block-encoded jump; the four-step circuit of Fig. circ:weak-measurement is the single-jump ($|J|=1$) restriction of their construction.
- **`chen2023efficient`** — Chen, C.-F., Kastoryano, M. J., Gilyén, A., *"An Efficient and Exact Noncommutative Quantum Gibbs Sampler"*, arXiv:2311.09207 (2023). The same weak-measurement primitive survives in the KMS-exact variant; only the jumps fed into $U_L$ change.

---

## Writing notes

<!-- For the author, not for the thesis -->

### Style / notation decisions

- **`\ket{0}_\delta`** single-qubit label and **`\ket{0^b}`** bundled block-register label: verbatim from the circuit figure's `\lstick` labels at `1_preliminaries.tex:1125, 1139`.
- **$Y_\delta$**: used throughout as the shorthand $R_Y(2\arcsin\sqrt\delta)$ per `1_preliminaries.tex:825`. The amplitude imprinted on $\ket{1}_\delta$ is $\sqrt\delta$ (not $\sin(\delta/2)\approx\delta/2$), which is why the dissipator prefactor is $\delta$ — not $\delta^2$ — in the rstick. The convention is made explicit inline when $Y_\delta$ is first used, with a pointer back to line 825.
- **Dissipator symbol**: renamed from $\mathcal{D}[L]$ to $\mathcal{D}_L$ to (a) differ visibly from the KMS discriminant $\mathcal{D}(\rho,\mathcal{L})$ already used at `1_preliminaries.tex:647, 655, 656` and `2_methods.tex:1940`, and (b) keep the jump argument as a subscript, parallel to $L \to L_a$ in Chapter 5 where the specialisation $\mathcal{D}_{L_a}$ writes itself naturally. The figure's rstick currently uses $\mathcal{D}(\rho)$ without argument — recommended to change to $\mathcal{D}_L(\rho)$ in the `.tex` for consistency; noted in Suggestions below.
- **Macros used**: `\ee`, `\ii`, `\dd`, `\idm` (all defined at `main.tex:83–90`).
- **$L$ versus $L_a$**: the figure and preliminaries discuss a single abstract jump $L$; Chapter 5's $L_a = \sqrt{\gamma(\bar\omega)}\,\tilde{A}_a(\bar\omega)$ is the concrete instantiation. The subsection sticks to $L$ throughout to match the figure caption.
- **"First Chen citation" pedagogical pointer**: the $M_0/M_1$ display carries a footnote-style tag `[CHEN23-EQ3.2]` in the draft, which expands in the `.tex` to `\cite[Eq.~(3.2)]{chen2023quantum}` (or the biblatex equivalent). This gives the reader a specific equation handle rather than a fifty-page-paper citation.

### Length

5 short paragraphs + 3 display equations (one labelled `eq:weak-meas-euler`). Rendered length is just under 1 page at the thesis font size — the shortest of the three preliminaries circuit subsections, as requested.

### Forward link

The subsection closes with an explicit pointer to (Alg. alg:diss), which at `2_methods.tex:1744` is labelled `alg:diss` and whose commentary at `2_methods.tex:1884–1931` already uses the coherent-undo pattern established here.

### Suggestions (figure-level edits the user should apply)

1. **\[CHECK — important\] Thesis figure vs Chen 2023 Fig. 3: control on $Y_\delta$.** The body text above follows Chen 2023 Theorem III.1 / their Figure 3 *as Chen drew it*: the $Y_\delta$ rotation is **block-controlled on $\ket{0^b}$**, so that it fires only on the post-selectable branch of $U_L$. The current thesis figure at `1_preliminaries.tex:1125–1131` draws $Y_\delta$ as **unconditional** on the $\ket{0}_\delta$ wire (no filled control dot from the block register). With $Y_\delta$ unconditional, the Kraus derivation above does *not* go through: one instead gets $M_0 = \sqrt{1-\delta}\,\idm$ (without the $L^\dagger L$ correction) and $M_1 = \sqrt\delta\,L/\alpha$, whose post-selected channel is $\rho\mapsto(1-\delta)\rho + (\delta/\alpha^2)\,L\rho L^\dagger$ — not a trace-preserving Lindblad-dissipator Euler step unless $L^\dagger L = \alpha^2\,\idm$.
   - **Recommended one-line TikZ fix at `1_preliminaries.tex:1125–1141`**: add a filled control dot from the block register (row 2) onto the $Y_\delta$ column (currently column 5), so that the $Y_\delta$ on row 1 becomes the target of a control carried by the block wire. Concretely, on row 1 the cell `\gate{Y_\delta}` in column 5 remains a target, and on row 2 one inserts a `\ctrl{-1}` (or the appropriate quantikz directive for a control on a *bundled* wire — likely the half-filled SELECT-style dot the thesis uses at `1_preliminaries.tex:827`) in that same column. The exact TikZ incantation depends on the bundle style the user is comfortable with; the mathematical requirement is a block-wire control on $Y_\delta$.
   - **Alternative (not recommended):** if the user prefers to keep $Y_\delta$ unconditional as currently drawn, the body text must be rewritten: the rstick would need to read $\rho + (\delta/\alpha^2)L\rho L^\dagger - \delta\,\rho + \mathcal{O}(\delta^2)$ (not a Lindblad-dissipator Euler step), and the subsection's own purpose — the weak-measurement implementation of $\mathcal{D}_L$ — would have to be abandoned for this figure. Not advisable, and inconsistent with Chapter 5's `alg:diss` step 8 at `2_methods.tex:1789–1791`, which is the block-controlled convention.

2. **Caption typo: $\delta^2 \to \delta$.** The caption at `1_preliminaries.tex:1165` currently writes $\gamma = \delta^2/(\alpha^2\,\dd t)$, while the rstick at line 1158 writes $\rho + (\delta/\alpha^2)\mathcal{D}(\rho) + \mathcal{O}(\delta^2)$ (linear in $\delta$), and Chapter 5 at `2_methods.tex:1910, 1929` reads $\rho + \delta\,\mathcal{L}_{a,\mathrm{diss}}(\rho) + \mathcal{O}(\delta^2)$. Four independent cross-checks (Chen 2023 Thm III.1 Eq. (3.2)–(3.3), `alg:diss` step 8's $R_Y(2\arcsin\sqrt\delta)$, `2_methods.tex:1910`, `2_methods.tex:1929`) confirm the Euler reading $\delta \sim \dd t$, which forces $\gamma = \delta/(\alpha^2\,\dd t)$. **Recommended one-character edit at `1_preliminaries.tex:1165`**: replace $\delta^2$ by $\delta$.

3. **Figure rstick: $\mathcal{D}(\rho)\to\mathcal{D}_L(\rho)$.** To align with the renamed dissipator $\mathcal{D}_L$ in the body (distinct from the KMS discriminant $\mathcal{D}(\rho,\mathcal{L})$), edit the rstick at `1_preliminaries.tex:1158` from $\mathcal{D}(\rho)$ to $\mathcal{D}_L(\rho)$, and the caption's $\mathcal{D}(\rho) = L\rho L^\dagger - \tfrac12\{L^\dagger L,\rho\}$ at line 1163 to $\mathcal{D}_L(\rho) = L\rho L^\dagger - \tfrac12\{L^\dagger L,\rho\}$. One-character edit in each.

4. **Optionally swap the labelled display to `eq:weak-meas`** or similar, to match the naming of sibling labels (`eq:qpe-unified` in QPE, `eq:B-block-encoding` in Chapter 5). The current choice `eq:weak-meas-euler` emphasises that this is specifically the Euler-step identity; either works, but the shorter form is more consistent with sibling labels.

---

## Review-driven choices applied

The following fixes from `drafts/weak-measurement-subsection-review.md` are applied to this v2 draft (not reproduced as issues):

- **Critical Issue 1 (figure-vs-text mismatch, Chen's block-controlled $Y_\delta$):** the body text has been rewritten so that the block-control is *explicit* in prose — step (ii) says "applied conditional on $\ket{0^b}$" and makes the two-branch split visible, and the $\ket{\phi^\perp}_{bS}$ residual is named in step (i) and re-entered in step (iii) so the Kraus derivation is traceable through the circuit as intended. The figure-level resolution is raised prominently as `[CHECK — important]` in Suggestion 1 above, with a concrete one-line TikZ fix and a spelled-out statement of the failure mode if the figure is kept as drawn.
- **Critical Issue 2 (walkthrough language for the $\ket{0}_\delta$ branch ignores leakage):** the previous "cleanly uncomputes $U_L$ and restores $\rho$" has been replaced by step (iii)'s "On the $\ket{0^b}$-block … this uncomputes to $\ket{0^b}\ket\psi$; the leaked $\ket{\phi^\perp}_{bS}$ is *not* uncomputed and survives as an out-of-block residual." This primes the Kraus derivation correctly (the reader knows where the $L^\dagger L$ term will come from).
- **N1 (single-index $L$):** kept $L$ without index; the $|J|=1$ specialisation of Chen's multi-jump formula is called out in the Citations section.
- **N2 (Chen 2023 vs 2025 bibkey):** kept `chen2023efficient`; internally consistent with thesis usage at `2_methods.tex:18`.
- **N3 ($\mathcal{D}$ notation clash with discriminant):** renamed $\mathcal{D}[L]\to\mathcal{D}_L$ and added a one-line disambiguating parenthetical pointing at `eq:discriminant` in §sec:davies. Also flagged in Suggestion 3 to update the figure rstick.
- **N4 (Theorem III.1 pointer sharpness):** the body's $M_0/M_1$ display carries a `[CHEN23-EQ3.2]` tag to pin the citation to the exact equation.
- **N5 (alg:diss label):** confirmed and kept.
- **N6 (macros):** `\ee`, `\ii`, `\dd`, `\idm` all preserved.
- **[CHECK A] caption typo $\delta^2 \to \delta$:** sharpened from the drafter's apologetic "corresponds to the diffusive-scaling reading" phrasing into a direct typo statement in Suggestion 2, with the four independent cross-checks enumerated.
- **[CHECK B] multi-jump defer:** kept deferred; single forward-pointing sentence at the end of the Euler-step paragraph.
- **[CHECK C] $M_0$ justification:** kept as one-liner + pin-pointed Chen citation, per the reviewer's preferred option under the "figure fixed" scenario (which is what the body's prose assumes).
- **Exposition Suggestion 1 (why-before-what opening):** adopted. Paragraph 1 now opens with what a Lindblad dissipator *is* (transition + drain) before identifying the drain as the missing piece.
- **Exposition Suggestion 2 ($\gamma$ vs $\gamma_{\mathrm{eff}}$):** adopted; the rate is written as $\gamma_{\mathrm{eff}} = 1/\alpha^2$ and the name $\gamma$ is left free for Chapter 5's transition weight $\gamma(\bar\omega)$ of the Boltzmann filter.
- **Exposition Suggestion 3 (OAA closure):** adopted as the phrase "no oblivious amplitude amplification is needed: $\delta$-smallness replaces amplification as the deterministic-acceptance mechanism".
- **Exposition Suggestion 4 ("failed post-selection" phrasing):** adopted; the sentence "bookkeeping the trace loss reconstitutes the anti-Hermitian drain" replaces the original "the failed post-selection on $b$ when $U_L^\dagger$ does not act" phrasing, with the $\ket{\phi^\perp}_{bS}$ residual named earlier.
- **Exposition Suggestion 5 (tighten the punchline):** adopted; "turns the LCU's $1-1/\alpha^2$ block-failure amplitude into signal: the leaked projection becomes the Lindblad drain."
- **Exposition Suggestion 6 (collapse $M_0/M_1$ to one display):** adopted; the two Kraus operators now sit on one display line.
- **Drafter's Sug. 2 ($\|L/\alpha\|\leq 1$ subnormalisation):** adopted as the six-word parenthetical "(with $\|L/\alpha\|\leq 1$ by §sec:prelim-LCU)".
- **E3 (link to Kraus form of §prelim.kraus):** adopted as "(in the sense of (Eq. eq:kraus))" in the Kraus operators paragraph.
- **Drafter's Sug. 1, 3, 4:** applied per the reviewer's verdicts (caption fix → Suggestion 2; `eq:weak-meas-euler`→`eq:weak-meas` rename → Suggestion 4 optional; composition-forward link → not added, per reviewer's advice to keep chapter division of labour clean).
- **E1, E2 (Markov-chain rejection, Lin/Fang connection):** not added; the subsection is the shortest of three and these would dilute focus without paying off in a preliminaries setting.
