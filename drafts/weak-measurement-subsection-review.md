# Draft Review: Weak-measurement subsection

**Reviewed**: `drafts/weak-measurement-subsection.md`
**Date**: 2026-04-23
**Overall assessment**: The prose is crisp, the structure is tight, and the choice to defer the multi-jump/KMS-weighted generalisation to Chapter 5 is right. However, there is one genuine *critical* issue that the drafter did not flag: **the thesis's Figure `circ:weak-measurement` deviates from Chen 2023 Figure 3 in a load-bearing way** — the $Y_\delta$ rotation is drawn *unconditional* in the thesis figure but is *controlled on $\ket{0^b}$* in Chen's original. Under the thesis's circuit as drawn, the derivation of $M_0 = \sqrt{1-\delta}\,I - (\delta/2\alpha^2)L^\dagger L + O(\delta^2)$ given in the draft (borrowed verbatim from Chen) does **not** go through, and the post-selected channel is not a Lindblad-dissipator Euler step. This has to be resolved before the subsection lands. Otherwise: the draft's verdict on the caption $\delta^2 \to \delta$ typo is correct and well-supported; notation/macros are clean; length is appropriately short.

---

## Critical Issues

### Issue 1: The thesis figure differs from Chen 2023 Fig. 3 — $Y_\delta$ is unconditional, breaking the Kraus derivation

- **Location**: Draft §"Circuit walkthrough" and §"Kraus operators"; underlying problem at `1_preliminaries.tex:1129`.
- **Problem**: The draft describes the thesis circuit correctly as "a $Y_\delta$-rotation on $\ket{0}_\delta$ … and an open-controlled $U_L^\dagger$", i.e. $Y_\delta$ is unconditional. But the Kraus formula $M_0 = \sqrt{1-\delta}\,I - (\delta/2\alpha^2)L^\dagger L + O(\delta^2)$ is lifted from Chen 2023 Theorem III.1, whose circuit (their Figure 3) has $Y_\delta$ **controlled on $\ket{0^b}$**. Tracing through the thesis circuit as literally drawn:

  1. $U_L|0^b\rangle|\psi\rangle = |0^b\rangle\otimes(L/\alpha)|\psi\rangle + |\phi^\perp\rangle$.
  2. Applying $Y_\delta$ unconditionally: $(\sqrt{1-\delta}|0\rangle_\delta + \sqrt\delta|1\rangle_\delta)\otimes(|0^b\rangle(L/\alpha)|\psi\rangle + |\phi^\perp\rangle)$.
  3. Applying $U_L^\dagger$ open-controlled on $q_\delta$: on the $|0\rangle_\delta$ branch, $U_L^\dagger$ fires and cleanly uncomputes to $|0^b\rangle|\psi\rangle$; on $|1\rangle_\delta$, nothing happens.
  4. Post-selecting $|0^b\rangle$ and measuring $q_\delta$:

     $$M_0 = \sqrt{1-\delta}\,I, \qquad M_1 = \sqrt{\delta}\,L/\alpha.$$

     There is **no** $-(\delta/2\alpha^2)L^\dagger L$ term in $M_0$. The post-selected channel is then
     $$\rho \mapsto (1-\delta)\rho + \frac{\delta}{\alpha^2}L\rho L^\dagger,$$
     which is not a trace-preserving Lindblad Euler step of $\mathcal{D}[L]$ unless $L^\dagger L = \alpha^2 I$ (i.e. $L/\alpha$ unitary).
- **Why it matters**: This is a logical error at the core of the subsection — the body text's claim that the circuit implements $\rho + (\delta/\alpha^2)\mathcal{D}[L]\rho + O(\delta^2)$ is false for the circuit as drawn. It also means the figure does not match Chen 2023 and does not match the usage in Chapter 5: `alg:diss` step 8 at `2_methods.tex:1789-1791` has $Y_\delta$ *controlled* on $q_\gamma = |0\rangle_\gamma$ (and implicitly, via the Boltzmann filter, on the block-encoding branch of the OFT), consistent with Chen's block-controlled $Y_\delta$.

  Chen's derivation (Chen 2023, proof of Theorem III.1, eq. (3.2)) gets the $-(\delta/2\alpha^2)L^\dagger L$ term precisely *because* the $(I - |0^b\rangle\langle 0^b|\otimes I)U|0^c\rangle|\psi\rangle$ component sits in the $|0\rangle$-branch only, and partial-tracing that branch over the block-failure outcomes yields the missing $\{L^\dagger L,\rho\}$ drain. Remove the $|0^b\rangle$-control on $Y_\delta$ and that bookkeeping vanishes.

- **Suggested fix**: The cleanest resolution is to **edit the thesis figure to match Chen Fig. 3**: add a block-register control dot on $Y_\delta$ so $Y_\delta$ fires only on the $|0^b\rangle$ branch. At `1_preliminaries.tex:1125-1131` change

  ```
  \lstick{$\ket{0}_\delta$}
    & \qw & \qw & \qw
    & \gate{Y_\delta}          <-- currently unconditional
    & \octrl{1} & \meter{}
  ```

  to a block-controlled form (using quantikz's `\ctrl` from the block register wire to the $Y_\delta$ gate), and add a corresponding receiver control on the block line so that $Y_\delta$'s cell on row 1 is a target and the block wire carries a $\bullet$ in that column. Then the draft's body text and Kraus derivation are correct as written.

  If the user prefers to keep the figure as drawn (unconditional $Y_\delta$) — which is a legitimate *simpler* circuit, one just not sufficient to realise $\mathcal{D}[L]$ — then the body text and rstick are wrong and must be rewritten. The draft's current derivation is not self-consistent with the figure as it stands.

- **Drafter's note**: this issue is not flagged in the writing notes, but it eclipses the $\delta$-vs-$\delta^2$ caption issue in importance. Recommend fixing both in the same figure-edit pass.

### Issue 2: Internal language inconsistency — "uncomputes $U_L$ and restores $\rho$" on the $|0\rangle_\delta$ branch ignores the failure leakage

- **Location**: Draft §"Circuit walkthrough" step (iii): "there it cleanly uncomputes $U_L$ and restores $\rho$ on $S$".
- **Problem**: Even under Chen's block-controlled $Y_\delta$, the $|0\rangle_\delta$ branch does *not* cleanly return $\rho$ on $S$. It returns $|0^b\rangle\otimes(L/\alpha)|\psi\rangle + |\phi^\perp\rangle$-style leakage weighted by $(1 - \sqrt{1-\delta})$, which after post-selecting $|0^b\rangle$ is what produces the $L^\dagger L$ correction in $M_0$. The sentence "cleanly uncomputes $U_L$ and restores $\rho$" reads as if no correction is needed — which then jars against the $M_0 = \sqrt{1-\delta}\,I - (\delta/2\alpha^2)L^\dagger L$ formula that appears two paragraphs later.
- **Why it matters**: The Kraus derivation is the technical heart of the subsection. If the walkthrough says "clean uncomputation" and the Kraus formula says "plus an $L^\dagger L$ correction from failed post-selection on $b$", the reader is left to reconcile the two.
- **Suggested fix**: Replace "there it cleanly uncomputes $U_L$ and restores $\rho$ on $S$" with something like: "there $U_L^\dagger$ fires: on the $\ket{0^b}$-block of $U_L|0^b\rangle|\psi\rangle$ this uncomputes to $\rho$, but the orthogonal component $(I - \ket{0^b}\!\bra{0^b})U_L|0^b\rangle|\psi\rangle$ is not uncomputed — it survives as a sub-normalised term whose failure to post-select on $|0^b\rangle$ contributes the $L^\dagger L$ drain in $M_0$." This also sets up the Kraus derivation naturally.

---

## Notation / Citation Fixes

### N1: $L$ vs $L_j$ — single-index notation

- **Location**: Draft §"Kraus operators", in the $M_0$ formula and $M_1 = \sqrt\delta\,L/\alpha$.
- **Thesis convention**: thesis figure at `1_preliminaries.tex:1161-1164` uses a single abstract $L$ (no index) and $\mathcal{D}(\rho) = L\rho L^\dagger - \tfrac12\{L^\dagger L,\rho\}$ (single-operator dissipator). Draft matches this.
- **Chen 2023 usage**: Chen's eq. (3.2) writes $\sum_{j\in J}L_j\rho L_j^\dagger - \tfrac12\{\sum_j L_j^\dagger L_j,\rho\}$, i.e. a multi-jump form that comes from $U$ being the block-encoding of a **set** of Lindblad operators via $(\langle 0^b|\otimes I)U(|0^c\rangle\otimes I) = \sum_j |j\rangle\otimes L_j$.
- **Fix**: None needed in the draft — the thesis figure correctly specialises to single $L$ and the draft follows. Just flagging so the reader is aware the Chen formula $M_0 = \sqrt{1-\delta}\,I - (\delta/2\alpha^2)L^\dagger L$ is the $|J|=1$ specialisation of Chen's $M_0 = (I - \frac{1-\sqrt{1-\delta}}{\alpha^2}\sum_j L_j^\dagger L_j)(\sqrt{\cdot})$ — the draft doesn't misattribute.

### N2: Chen 2023 vs Chen 2025 bibkey

- **Location**: Draft Citations section and inline.
- **Draft claims**: `chen2023efficient` is the "2023 KMS-exact" paper, citing arXiv:2311.09207 as 2023.
- **Checked**: `arXiv:2311.09207` was originally submitted in Nov 2023 (hence `chen2023efficient`) but the published PRX article carries a 2025 date, and the PDF file in `supplementary-informations/` is actually titled "Chen et al. - 2025 - An efficient and exact noncommutative quantum Gibbs sampler.pdf". The thesis itself at `2_methods.tex:18` uses `chen2023efficient` and the corresponding bibtex key is well-established in the thesis — no fix needed in the draft itself; the *.bib* entry already exists. But the reference-list prose in the draft ("Chen, Kastoryano, Gilyén, 2023") is fine because the bibkey year (2023) is the arXiv year.
- **Fix**: None needed; noting in case a future reviewer wants to rename the key to `chen2025efficient`. Current draft is internally consistent.

### N3: $\mathcal{D}(\rho)$ vs $\mathcal{D}[L](\rho)$

- **Location**: Draft's opening display ("$\mathcal{D}[L](\rho) := L\rho L^\dagger - \tfrac12\{L^\dagger L,\rho\}$"); thesis figure rstick at line 1158 uses $\mathcal{D}(\rho)$ without the $[L]$; thesis caption at line 1163 uses $\mathcal{D}(\rho)$ without the $[L]$.
- **Thesis convention**: Note that `$\mathcal{D}$` is already heavily used in the thesis as **discriminant** (`1_preliminaries.tex:647, 655, 656`; `2_methods.tex:1940`) — `$\mathcal{D}(\rho, \mathcal{L})$` and `$\mathcal{D}(\rho_\beta, \mathcal{L})$`. Using `$\mathcal{D}[L](\rho)$` for the Lindblad dissipator, plain `$\mathcal{D}(\rho)$` for the same object inside the figure, **and** `$\mathcal{D}(\rho, \mathcal{L})$` for the discriminant is three different meanings of a single symbol in close proximity.
- **Fix**: Either (a) rename the dissipator with different notation, e.g. `$\mathscr{L}_\mathrm{diss}[L](\rho)$` or simply `$\mathcal{D}_L(\rho)$` throughout; or (b) keep `$\mathcal{D}[L]$` but explicitly remark that this *dissipator* `$\mathcal{D}$` is distinct from the *discriminant* `$\mathcal{D}(\rho, \mathcal{L})$` of `eq. 655`. The draft currently does (b) implicitly by always writing `$\mathcal{D}[L]$` with the bracketed argument, and matches the figure's `$\mathcal{D}(\rho)$` which has a parenthetical. A one-sentence footnote or parenthetical — "(distinguished from the KMS discriminant of Eq. 655 by the explicit jump argument)" — would go a long way. Low priority, but worth doing before the section goes to print.

### N4: Citation correctness — Theorem III.1 attribution

- **Location**: Draft: "a standard amplitude-amplification bookkeeping [CITE: chen2023quantum]".
- **Verified**: Chen 2023 Theorem III.1 proves the weak-measurement construction for **general** purely-irreversible Lindbladians $\mathcal{L}[\cdot] = \sum_{j\in J} L_j \cdot L_j^\dagger - \tfrac12\{\sum_j L_j^\dagger L_j,\cdot\}$, via the Kraus derivation at eq. (3.2). The draft's "single-jump restriction" phrasing is exactly right. Citation is accurate.
- **Fix**: None. But consider adding a more specific pointer: "see the proof of Theorem III.1 at Chen et al. 2023, eq. (3.2)" rather than the bare citation, because the Kraus-by-Kraus bookkeeping is exactly what's relevant and it's easy to miss in a fifty-page paper.

### N5: Forward reference — `alg:diss` is correct

- **Verified**: `2_methods.tex:1744` has `\label{alg:diss}` inside `\begin{algorithm}`. Draft uses `(Alg. alg:diss)` which is the right target.
- **Fix**: None.

### N6: Macros `\ee`, `\ii`, `\dd`, `\idm`

- **Checked**: All four appear in the draft; `\idm` is used in place of `$I$` for the identity in $M_0$ and $M_1$; `\ee^{\delta\mathcal{D}[L]/\alpha^2}$` uses `\ee` correctly; `\dd t` in the $\gamma$ formula. All match thesis macros.
- **Fix**: None.

---

## [CHECK] resolutions

### [CHECK A]: $\delta$ vs $\delta^2$ in figure caption — drafter's verdict is correct

**Verdict**: The drafter's conclusion is correct. The figure **caption** at `1_preliminaries.tex:1165` has a typo: `$\gamma = \delta^2/(\alpha^2\,\dd t)$` should read `$\gamma = \delta/(\alpha^2\,\dd t)$`. The rstick at line 1158 and the body of Chapter 5 are both in the Euler-reading ($\delta \sim \dd t$).

**Evidence**:
1. **Chapter 5 `alg:diss` step 8** (`2_methods.tex:1789-1791`): "rotate $q_\delta$: $R_Y(2\arcsin\sqrt\delta)$" — this is exactly the $Y_\delta$ of the preliminaries with the thesis's $Y_\theta := R_Y(2\arcsin\sqrt\theta)$ convention at line 825. Amplitude on $|1\rangle$ is $\sqrt\delta$ (not $\sqrt{\delta^2} = \delta$, and not $\sin(\delta/2) \approx \delta/2$).
2. **`2_methods.tex:1910`**: "The net probability of the 'jump' outcome ($q_\gamma = |1\rangle_\gamma$) is $\gamma(\bar\omega)\cdot\delta$, i.e. **first-order in $\delta$** which is required by the weak-measurement approximation $e^{\delta\mathcal{L}_{a,\mathrm{diss}}} = \mathrm{id} + \delta\mathcal{L}_{a,\mathrm{diss}} + \mathcal{O}(\delta^2)$." — explicitly first-order in $\delta$.
3. **`2_methods.tex:1929`**: "$\rho \mapsto \rho + \delta\,\mathcal{L}_{a,\mathrm{diss}}(\rho) + \mathcal{O}(\delta^2)$" — exactly the Euler reading, matching the thesis figure's rstick `$\rho + (\delta/\alpha^2)\mathcal{D}(\rho) + \mathcal{O}(\delta^2)$`.
4. **Chen 2023 Thm III.1** (see Chen pages 19-20 eq. (3.2)-(3.3)): the per-step Kraus sum gives $(\mathcal{I} + \delta\mathcal{L})[|\psi\rangle\langle\psi|] + \mathcal{O}(\delta^2)$. First-order in $\delta$. Chen then sets $\delta = \Theta(\epsilon/t)$ and repeats $\Theta(t^2/\epsilon)$ times — the step-size $\delta$ is the Euler step, not its square.

All four independent cross-checks agree: **the caption has the typo, not the body**. The caption should read `$\gamma = \delta/(\alpha^2\,\dd t)$`.

Drafter's phrasing in the body ("see [CHECK] in the notes: the figure caption writes $\gamma = \delta^2/(\alpha^2\,\dd t)$, which corresponds to the diffusive-scaling reading") is accurate but somewhat apologetic for the caption-bug. Recommend sharpening: the in-body parenthetical should say "(the figure caption at `1165` contains a typo: $\delta^2\to\delta$)".

### [CHECK B]: Multi-jump extension deferred — good call

**Verdict**: Defer is correct. Inserting the multi-jump generalisation in the preliminaries would duplicate Chapter 5's `alg:diss` commentary verbatim. The single forward-pointing sentence at the end of the Euler-step paragraph is the right touch.

### [CHECK C]: $M_0$ formula justification

**Verdict**: The one-line justification + Chen citation is borderline. Given the derivation in Chen eq. (3.2) is exactly the bookkeeping that gives the $L^\dagger L$ term from the $(I - |0^b\rangle\langle 0^b|\otimes I)U|0^c\rangle|\psi\rangle$ residual, a half-paragraph unpacking would be pedagogically worth it — but is **only if the figure is fixed to match Chen's block-controlled $Y_\delta$** (Critical Issue 1). If the thesis figure stays as drawn (unconditional $Y_\delta$), the one-liner + Chen citation is actually *misleading*, because Chen's derivation does not apply.

Two sub-options under "figure fixed" scenario:
- **Keep as is** (one line, Chen citation): fine for a preliminaries subsection that is the shortest of three. Reader is expected to chase Chen Thm III.1 for the bookkeeping.
- **Unpack to half-paragraph**: four sentences of the form "Writing $U_L|0^b\rangle|\psi\rangle = |0^b\rangle\otimes(L/\alpha)|\psi\rangle + (I-|0^b\rangle\langle 0^b|\otimes I)U_L|0^b\rangle|\psi\rangle =: |0^b\rangle\otimes(L/\alpha)|\psi\rangle + |\phi^\perp\rangle$, ..., after post-selecting $b$ on $|0^b\rangle$ the $\sqrt{1-\delta}|0\rangle_\delta$-branch contributes $\sqrt{1-\delta}\,\langle 0^b|U_L^\dagger U_L|0^b\rangle = \sqrt{1-\delta}\,I$, while the $\sqrt\delta|1\rangle_\delta$-branch has $|0^b\rangle$ post-selection failure probability $\|(I-|0^b\rangle\langle 0^b|)U_L|0^b\rangle|\psi\rangle\|^2 = 1 - \|L|\psi\rangle\|^2/\alpha^2$, and bookkeeping the latter into the $M_0$ Kraus map via the $Y_\delta$-rotation's $|1\rangle$-branch gives the $-(\delta/2\alpha^2)L^\dagger L$ correction." This is precise but adds ~4 sentences; probably tips the subsection over the "shortest of the three" goal. Depends on the user's taste.

Recommend: **keep the one-liner** for now (shortest subsection goal), but pin a `\todo{consider half-paragraph derivation if reviewers push back}` in the draft.

---

## Consistency Issues

### C1: Figure caption and rstick disagree with each other (the [CHECK A] above)

Already covered. Affects `1_preliminaries.tex:1158` and `:1165` — one typo fix, one-line edit.

### C2: Thesis figure diverges from Chen 2023 Fig. 3 (Critical Issue 1)

Already covered. Affects `1_preliminaries.tex:1125-1131` — either fix the figure or fix the body text + caption.

### C3: `qpe-subsection.md` references `\subsection{Quantum Phase Estimation.}` inside "Chapter Quantum Circuits"

This is a QPE draft issue, not a weak-meas issue, but flagging: the three preliminaries subsections (LCU, QSP, weak-measurement) all live inside the Quantum Circuits chapter of `1_preliminaries.tex`. All three drafts are internally consistent about this, and all use `sec:prelim-*` label prefixes uniformly. No fix.

### C4: `$\alpha$` as subnormalisation is consistent with `lcu-subsection.md`

Draft's "$\|L/\alpha\|\leq 1$" implicit from LCU — matches `drafts/lcu-subsection.md` Eq. (lcu-block-encoding) and `2_methods.tex:1724`. No fix.

### C5: `$Y_\theta := R_Y(2\arcsin\sqrt\theta)$` convention

Draft uses the thesis convention correctly. Amplitude on $|1\rangle$ is $\sqrt\delta$. Verified at `1_preliminaries.tex:825`. No fix.

---

## Exposition Suggestions (ranked by impact)

### Suggestion 1: Open with the why-before-what — first state what a Lindblad semigroup needs, then the block encoding delivers one half

- **Location**: Draft paragraph 1.
- **Current**: "The Linear Combination of Unitaries of §sec:prelim-LCU block-encodes a non-unitary operator $L$ … and so delivers the amplitude $L|\psi\rangle/\alpha$ … What it does *not* deliver is the accompanying anti-Hermitian drain …".
- **Suggested**: Flip the order: open with "A Lindblad dissipator $\mathcal{D}[L](\rho)$ consists of a *transition* term $L\rho L^\dagger$ and an anti-Hermitian *drain* $-\tfrac12\{L^\dagger L,\rho\}$ that together form a trace-preserving channel. The LCU of §sec:prelim-LCU gives us the transition term via the block encoding $U_L$; the drain is what the weak-measurement construction of Chen et al. 2023 supplies." Then state the Kraus-pair result as the single-sentence summary. Rationale: the reader is told what the dissipator is (structural definition) before the specific missing piece (the drain) is named.
- Impact: medium. Sharpens motivation; adds ~20 words.

### Suggestion 2: Sharpen the rate-rewrite at the end of the Euler paragraph

- **Location**: Draft §"One Euler step", last sentence: "Identifying $\delta$ as the Euler step-size of the semigroup then gives the effective Lindblad rate $\gamma = \delta/(\alpha^2\,\dd t)$".
- **Current**: The sentence reads as a definition of $\gamma$, but $\gamma$ is already overloaded in Chapter 5 as the *transition weight* $\gamma(\bar\omega)$ of the Boltzmann filter — a quite different object.
- **Suggested**: "Identifying $\delta = \gamma_{\mathrm{eff}}\,\dd t$ with $\gamma_{\mathrm{eff}} := 1/\alpha^2$ the effective Lindblad *rate* of the semigroup (per unit time), the channel implements one Euler step of $\dot\rho = \gamma_{\mathrm{eff}}\,\mathcal{D}[L](\rho)$." Rationale: avoids overloading $\gamma$, makes the role of $\alpha$ clearer (it's what sets the natural timescale of the block-encoded generator).
- Impact: low-medium. Saves a reader hiccup. Combined with a note that Chapter 5's $\gamma(\bar\omega)$ multiplies $\delta$ to give the *jump* probability $\gamma(\bar\omega)\cdot\delta$, it also foreshadows the Boltzmann filter.

### Suggestion 3: One-sentence connection to oblivious amplitude amplification

- **Location**: New sentence near the end of §"One Euler step" or in the final forward-pointing sentence.
- **Current**: No reference to amplification.
- **Suggested**: "(Unlike the LCU post-selection of §sec:prelim-LCU, this one is near-deterministic per step because $\delta\ll 1$; no oblivious amplitude amplification is needed.)" Rationale: §sec:prelim-LCU motivates amplification as the cure for small post-selection probabilities; the weak-measurement subsection should close the loop by pointing out that here, $\delta$-smallness replaces amplification as the deterministic-acceptance mechanism.
- Impact: low. Nice-to-have for continuity.

### Suggestion 4: The "failed post-selection on $b$ when $U_L^\dagger$ does not act" phrasing is slightly awkward

- **Location**: Draft §"Kraus operators": "the second term comes from the failed post-selection on $b$ when $U_L^\dagger$ does not act".
- **Problem**: "failed post-selection" is the outcome, not the cause; the cause is that $(I - |0^b\rangle\langle 0^b|\otimes I)U_L|0^b\rangle|\psi\rangle$ component of $U_L|0^b\rangle|\psi\rangle$ stays entangled in the $|1\rangle_\delta$-branch (under Chen's circuit) and partial-tracing over its block-failure part gives the $L^\dagger L$ contribution back to the $|0\rangle_\delta$-branch's $M_0$.
- **Suggested**: "the second term arises from the $(I - |0^b\rangle\!\langle 0^b|\otimes I)\,U_L|0^b\rangle|\psi\rangle$ residual on the $|1\rangle_\delta$-branch: its post-selection on $|0^b\rangle$ fails with probability $\|L|\psi\rangle\|^2/\alpha^2$, and bookkeeping the trace loss back into the $|0\rangle_\delta$-branch reconstitutes the anti-Hermitian drain."
- Impact: medium (and this is really only applicable if Critical Issue 1 is resolved by fixing the figure to match Chen's block-controlled $Y_\delta$).

### Suggestion 5: Make the "trades one failure mode for the dissipator" punchline more punchy

- **Location**: Last sentence of §"One Euler step": "the entire construction trades one failure mode of the LCU block-encoding (the $1-1/\alpha^2$ rejection probability) for the dissipator that the simulation needs".
- **Current**: Correct but slightly tangled.
- **Suggested**: "In effect, the coherent-undo construction turns the LCU's $1-1/\alpha^2$ post-selection loss into signal: the failed-projection amplitude becomes the Lindblad drain." Rationale: keeps the same message in fewer words, puts "signal" as the payoff.
- Impact: low. Cosmetic polish.

### Suggestion 6: Length — the subsection could trim 10-15%

- **Current**: 5 paragraphs plus 3 displays; ~1 page (estimate).
- **Observations**: Paragraph 3 ("The two branches of the ancilla measurement…") and paragraph 4 ("Summing the two branches…") overlap significantly. The $M_0$ and $M_1$ displays are each preceded by a full sentence; the drafter's own note (Length section, last line) offers a one-line collapse `$M_0 = \sqrt{1-\delta}\,\idm + \mathcal{O}(\delta/\alpha^2)$, $M_1 = \sqrt\delta\,L/\alpha$`.
- **Suggested**: Adopt the collapsed form. Then the paragraph reads: "The two branches give Kraus operators $M_0 = \sqrt{1-\delta}\,\idm + \mathcal{O}(\delta/\alpha^2)$ (the $O(\delta/\alpha^2)$ is the $-(\delta/2\alpha^2)L^\dagger L$ failed-post-selection correction, [Chen et al. 2023, Thm III.1]) and $M_1 = \sqrt\delta\,L/\alpha$." Saves half a paragraph, doesn't lose content.
- Impact: medium.

---

## Connections and Enrichments (optional)

### E1: Connect to Markov-chain rejection sampling

The draft has a latent parallel to classical Metropolis-Hastings (acceptance-rejection): "no-jump" ($|0\rangle_\delta$) = accept-current-state; "jump" ($|1\rangle_\delta$) = accept-proposed-transition. This parallel is half-made in Chapter 5's `alg:diss` ("net jump probability $\gamma(\bar\omega)\cdot\delta$"), and Chen 2023 Fig. 1 explicitly draws the Markov-chain picture. A one-sentence remark in the weak-meas subsection could reinforce this — "Operationally, the two ancilla outcomes implement a quantum analogue of accept ($|0\rangle_\delta$)/jump ($|1\rangle_\delta$), matching the classical MH kernel of (Eq. X.X)." Optional; strengthens the pedagogical arc from §classical-MCMC to §weak-meas.

### E2: Lin 2025 / Fang et al. 2026 connection

[Lin 2025] and [Fang et al. 2026] use the weak-measurement primitive for *dissipative* state preparation with detectability-lemma arguments. A half-sentence at the end of the subsection — "This primitive also appears in recent work on dissipative ground-state preparation [CITE: lin2025dissipative, fang2026detectability]" — would place the construction in the broader landscape. Optional; the subsection is already citing the Chen lineage and can't carry too many forward references to stay short.

### E3: Relate to Kraus form of quantum channels from §prelim.kraus

The draft implicitly uses the Kraus representation $\rho \mapsto \sum_k M_k\rho M_k^\dagger$ but does not reference the earlier thesis introduction of Kraus operators (§post-measurement state / eq:kraus at `1_preliminaries.tex:47`-ish range). A one-word pointer — "gives two Kraus operators (see Eq. eq:kraus)" — is near-free and strengthens the thesis's internal cohesion.

---

## Drafter's suggestions — independent assessment

### Drafter's Sug. 1: Fix the caption $\delta^2 \to \delta$ in `1_preliminaries.tex:1165`

**Agree strongly.** Confirmed as a typo (see [CHECK A]). One-character edit. Remove a genuine internal inconsistency.

### Drafter's Sug. 2: Mention $\|L/\alpha\|\leq 1$ (block-encoding subnormalisation) in the opening paragraph

**Lean in favour**, for a different reason than the drafter: the subnormalisation is what makes the $M_0 = \sqrt{1-\delta}\,I + O(\delta/\alpha^2)$ bound *hold uniformly in the state* (because $\|L|\psi\rangle\|^2/\alpha^2 \leq 1$ automatically). Six words — "(with $\|L/\alpha\|\leq 1$ by §prelim-LCU)" — after the introduction of $U_L$ is enough.

### Drafter's Sug. 3: Rename label `eq:weak-meas-euler` to `eq:weak-meas`

**Indifferent.** Both names are fine; `eq:weak-meas-euler` is more informative because it emphasises that the Euler reading is what the equation says. If the user wants a shorter label matching `eq:qpe-unified` and `eq:B-block-encoding`, renaming is trivial. Low priority.

### Drafter's Sug. 4: Add a DissipativeStep-composition forward link

**Disagree, for now.** The current division of labour (preliminaries = three circuit primitives; Chapter 5 = their composition) is cleaner. Adding a composition-forward sentence here anticipates structure the reader has not yet met and duplicates what Chapter 5 does in the opening of `alg:diss` commentary.

---

## Summary Scorecard

| Category | Rating | Notes |
|----------|--------|-------|
| Logical correctness | **MAJOR** | Critical Issue 1: thesis figure's unconditional $Y_\delta$ contradicts Chen-style Kraus derivation. Either fix figure or rewrite body. |
| Notation consistency | MINOR | `$\mathcal{D}[L]$` collides with the thesis's discriminant `$\mathcal{D}(\rho,\mathcal{L})$`; minor. |
| Citations | OK | `chen2023quantum` and `chen2023efficient` both present and used correctly; `alg:diss` label confirmed. |
| Thesis consistency | MINOR | Caption $\delta^2/(\alpha^2 \dd t)$ at line 1165 is a typo (confirmed against Chapter 5 and Chen); fix in the figure fix pass. |
| Exposition quality | MINOR | Tightly written overall; suggest minor reorder in the opening paragraph (Suggestion 1) and the Kraus-formula phrasing (Suggestion 4). |

---

**Verdict and per-bucket counts**: Critical 2, Notation/citation 4 (2 fixes needed, 2 informational), Exposition 6 (all non-blocking), [CHECK] 3 (all resolved: A = caption typo confirmed, B = defer correct, C = one-liner OK conditional on Critical 1 fix). The drafter's write-up is technically sharp and structurally correct *for Chen's circuit*, but the thesis Figure at `1_preliminaries.tex:1125-1131` draws $Y_\delta$ unconditional rather than block-controlled; this has to be reconciled before the subsection is merged. The $\delta$-vs-$\delta^2$ caption issue is a genuine typo (body text is right), and the fix is a one-character edit at `1_preliminaries.tex:1165`.
