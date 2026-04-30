# Weak-measurement (v3)

> **Insertion target:** `supplementary-informations/1_preliminaries.tex`, body of `\subsection{Weak-measurement.} \label{sec:prelim-weak-meas}` — text goes between the `\subsection{...}\label{...}` line at line 1116 and `\begin{figure}[ht]` at line 1117 (the figure `circ:weak-measurement` sits at lines 1117–1167).
> **Figure replacement:** the thesis figure at `1_preliminaries.tex:1117–1167` must be replaced with `drafts/circuits/weak-measurement-v2.tex` — the current figure draws $Y_\delta$ as unconditional (inconsistent with Chen 2023 Thm III.1, the construction the body follows), and the new file adds the missing block-register open control and fixes the caption typo.
> **Depends on labels already in the chapter:** `sec:prelim-LCU`, `sec:prelim-QPE`, `eq:kraus` (Kraus form of CPTP maps at `1_preliminaries.tex:47`-ish), `circ:weak-measurement`, and the $Y_\theta := R_Y(2\arcsin\sqrt\theta)$ shorthand at `1_preliminaries.tex:825`.
> **Forward references:** `alg:diss` (the `DissipativeStep` algorithm) in Chapter 5 at `2_methods.tex:1744`, and `sec:algorithm` at `2_methods.tex:1620`-ish.
> **Bib keys used, already in `references.bib`:** `chen2023quantum`, `chen2023efficient`. No new keys introduced.

---

A Lindblad dissipator $\mathcal{D}_L(\rho) := L\rho L^\dagger - \tfrac{1}{2}\{L^\dagger L,\rho\}$ combines a *transition* term $L\rho L^\dagger$ with an anti-Hermitian *drain* $-\tfrac12\{L^\dagger L,\rho\}$ whose interplay is what makes $\rho \mapsto \ee^{t\mathcal{D}_L}\rho$ trace-preserving. The LCU of §sec:prelim-LCU delivers the transition term: its post-selection on $\ket{0^b}$ of the block-encoded $U_L$ imprints $L\ket\psi/\alpha$ on $S$ (with $\|L/\alpha\|\leq 1$ by §sec:prelim-LCU). What LCU alone does *not* deliver is the drain. The coherent-undo construction of [Chen et al. 2023] [CITE: chen2023quantum] fills this gap by weaving a single weak-measurement ancilla $\ket{0}_\delta$ between $U_L$ and a conditional $U_L^\dagger$, in such a way that post-selecting $b$ on $\ket{0^b}$ and measuring $\ket{0}_\delta$ implements exactly one Euler step of $\ee^{\delta\mathcal{D}_L/\alpha^2}$ on $\rho$, up to $\mathcal{O}(\delta^2)$ (Fig. circ:weak-measurement). It is this primitive that reappears as the dissipative step of Chapter 5's Lindblad simulation (Alg. alg:diss). (The present $\mathcal{D}_L$ is distinct from the KMS discriminant $\mathcal{D}(\rho,\mathcal{L})$ of (Eq. eq:discriminant) in §sec:davies; the jump argument keeps the two apart.)

## Circuit walkthrough

Initialise the three registers in $\ket{0}_\delta \otimes \ket{0^b} \otimes \rho$ and apply the four operations of Fig. circ:weak-measurement in turn.

*(i) Block encoding.* $U_L$ acts on $(b,S)$; on the $\ket{0^b}$-block of the combined state it imprints the amplitude $L/\alpha$ on $S$, leaving an orthogonal component $\ket{\phi^\perp}_{bS} := (\idm - \ket{0^b}\!\bra{0^b}\!\otimes\idm)\,U_L\,\ket{0^b}\!\ket\psi$ outside the block.

*(ii) Ancilla rotation, controlled on the block.* The $Y_\delta$ of the thesis's $Y_\theta := R_Y(2\arcsin\sqrt\theta)$ shorthand at line 825 is *open-controlled on the block register*, i.e. it fires only on the $\ket{0^b}$-branch of (i):

$$Y_\delta\,\bigl(\ket{0}_\delta\otimes\ket{0^b}\bigr) \;=\; \bigl(\sqrt{1-\delta}\,\ket{0}_\delta + \sqrt{\delta}\,\ket{1}_\delta\bigr)\otimes\ket{0^b},\qquad Y_\delta\,\bigl(\ket{0}_\delta\otimes\ket{\phi^\perp}_{bS}\bigr) \;=\; \ket{0}_\delta\otimes\ket{\phi^\perp}_{bS},$$

so the weak-measurement ancilla acquires amplitude $\sqrt{\delta}$ on $\ket{1}_\delta$ precisely where $U_L$ succeeded, and carries $\ket{0}_\delta$ untouched on the leaked $\ket{\phi^\perp}_{bS}$ branch.

*(iii) Open-controlled undo.* The open-controlled $U_L^\dagger$ fires only on the $\ket{0}_\delta$ branch. On the $\ket{0^b}$-block of $U_L\ket{0^b}\ket\psi$ this uncomputes to $\ket{0^b}\ket\psi$; the leaked $\ket{\phi^\perp}_{bS}$ is *not* uncomputed and survives as an out-of-block residual. On the $\ket{1}_\delta$ branch $U_L^\dagger$ is absent, and $\ket{0^b}\otimes(L/\alpha)\ket\psi$ stands.

*(iv) Post-selection and measurement.* Post-select $b$ on $\ket{0^b}$ and measure $\ket{0}_\delta$.

## Kraus operators

The two outcomes of step (iv) give two Kraus operators on $S$ (in the sense of (Eq. eq:kraus)). In the *no-jump* branch ($\ket{0}_\delta$ outcome, $\ket{0^b}$ post-selected), the uncomputation $U_L^\dagger U_L = \idm$ acts on the block, while the $\ket{\phi^\perp}_{bS}$ residual has $\ket{0^b}$ post-selection failure probability $\|L\ket\psi\|^2/\alpha^2$; bookkeeping the trace loss reconstitutes the anti-Hermitian drain. In the *jump* branch ($\ket{1}_\delta$ outcome, $\ket{0^b}$ post-selected), the full amplitude $L\ket\psi/\alpha$ survives, weighted by the $\sqrt\delta$ of the rotation. The Chen et al. 2023 amplitude bookkeeping [CITE: chen2023quantum][CHEN23-EQ3.2] gives

$$M_0 \;=\; \sqrt{1-\delta}\,\idm \;-\; \tfrac{1}{2}\,\tfrac{\delta}{\alpha^2}\,L^\dagger L \;+\; \mathcal{O}(\delta^2), \qquad M_1 \;=\; \sqrt{\delta}\;\frac{L}{\alpha}.$$

## One Euler step

Summing the two branches and expanding in $\delta$,

$$\rho \;\mapsto\; M_0\,\rho\,M_0^\dagger + M_1\,\rho\,M_1^\dagger \;=\; \rho \;+\; \frac{\delta}{\alpha^{2}}\Bigl(L\rho L^\dagger - \tfrac{1}{2}\{L^\dagger L,\rho\}\Bigr) \;+\; \mathcal{O}(\delta^{2}),$$ <!-- \label{eq:weak-meas-euler} -->

which is one Euler step of $\dot\rho = (1/\alpha^2)\,\mathcal{D}_L(\rho)$ — the rstick of Fig. circ:weak-measurement. Identifying $\delta = \gamma_{\mathrm{eff}}\,\dd t$ with $\gamma_{\mathrm{eff}} := 1/\alpha^2$ the effective Lindblad rate of the block-encoded generator (per unit time), the post-selected channel implements one Euler step of $\dot\rho = \gamma_{\mathrm{eff}}\,\mathcal{D}_L(\rho)$. The post-selection on $b$ succeeds with probability $1 - \delta\,\|L\ket\psi\|^2/\alpha^2 + \mathcal{O}(\delta^2) = 1 - \mathcal{O}(\delta)$, so — unlike the $1/\alpha^2$ post-selection loss of §sec:prelim-LCU — no oblivious amplitude amplification is needed: $\delta$-smallness replaces amplification as the deterministic-acceptance mechanism. In effect, the coherent-undo construction turns the LCU's $1-1/\alpha^2$ block-failure amplitude into signal: the leaked projection becomes the Lindblad drain. The extension to the multi-jump sum $\mathcal{L} = \sum_a \mathcal{D}_{L_a}$ and to the KMS-weighted filtered jumps of [Chen et al. 2025] [CITE: chen2023efficient] that replace $L$ in Chapter 5 is the dissipative step of (Alg. alg:diss).

## Generality of the control structure

The figure above depicts the *minimal* LCU-block-encoded jump: a single block register $\ket{0^b}$ plays the role of "$U_L$ succeeded", and $Y_\delta$ is open-controlled on $\ket{0^b}$ to restrict the weak rotation to the post-selection branch. Chapter 5's dissipative step (Alg. alg:diss, §sec:algorithm) replaces this minimal block register with a richer ancilla structure: the frequency register $\Omega$ and the Boltzmann qubit $q_\gamma$ together encode the KMS-weighted jump $\sqrt{\gamma(\bar\omega)}\,\tilde{A}_a(\bar\omega)$ on the system, and $Y_\delta$ is then open-controlled on $q_\gamma = \ket{0}_\gamma$ in place of the generic $\ket{0^b}$ (Alg. alg:diss step 8, `2_methods.tex:1789`). The mechanism is the same — $Y_\delta$ always fires on the ancilla branch flagging successful application of the effective jump — and the Euler-step identity (Eq. eq:weak-meas-euler) persists with $L$ replaced by the OFT-filtered jump $\sqrt{\gamma(\bar\omega)}\,\tilde{A}_a(\bar\omega)$. In Chapter 5's picture, $\Omega$ is not post-selected but coherently uncomputed together with $q_\gamma$ on the no-jump branch (Alg. alg:diss steps 9–14), so the role of the block register's post-selection is taken over by the controlled uncomputation $U_{\mathrm{diss}}^\dagger$; the $Y_\delta$-control structure is inherited verbatim.

---

## Citations

Both keys already exist in `references.bib` — no new bibtex entries are introduced.

- **`chen2023quantum`** — Chen, C.-F., Kastoryano, M. J., Brandão, F. G. S. L., Gilyén, A., *"Quantum Thermal State Preparation"*, arXiv:2303.18224 (2023). Theorem III.1 and Eq. (3.2) in their proof establish the weak-measurement implementation of a Lindblad dissipator from a block-encoded jump; the circuit of Fig. circ:weak-measurement is the single-jump ($|J|=1$) restriction of their construction, and their Fig. 3 is the canonical drawing of the block-controlled $Y_\delta$ convention.
- **`chen2023efficient`** — Chen, C.-F., Kastoryano, M. J., Gilyén, A., *"An Efficient and Exact Noncommutative Quantum Gibbs Sampler"*, arXiv:2311.09207 (2023). The same weak-measurement primitive survives in the KMS-exact variant; only the jumps fed into $U_L$ change.

---

## Writing notes

<!-- For the author, not for the thesis -->

### Style / notation decisions

- **`\ket{0}_\delta`** single-qubit label and **`\ket{0^b}`** bundled block-register label: verbatim from the circuit figure's `\lstick` labels at `1_preliminaries.tex:1125, 1139`.
- **$Y_\delta$**: used throughout as the shorthand $R_Y(2\arcsin\sqrt\delta)$ per `1_preliminaries.tex:825`. The amplitude imprinted on $\ket{1}_\delta$ is $\sqrt\delta$ (not $\sin(\delta/2)\approx\delta/2$), which is why the dissipator prefactor is $\delta$ — not $\delta^2$ — in the rstick. The convention is made explicit inline when $Y_\delta$ is first used, with a pointer back to line 825.
- **Block-control as open control on a bundle.** The block register's control on $Y_\delta$ is the *open* control (`\octrl`), firing on $\ket{0^b}$ — matching the thesis's open-dot convention at line 827. On the bundled wire this reads operationally as "all qubits of $b$ are zero", i.e. the post-selectable branch. Chen 2023 Fig. 3 draws the same control as a filled disk because their register labels are ket-carrying (they write "controlled on the $\ket{0^b}$ state" in the proof); both conventions are standard, and the thesis's open-dot notation is the internally consistent one.
- **Dissipator symbol**: $\mathcal{D}_L(\rho)$ throughout, to (a) differ visibly from the KMS discriminant $\mathcal{D}(\rho,\mathcal{L})$ already used at `1_preliminaries.tex:647, 655, 656` and `2_methods.tex:1940`, and (b) keep the jump argument as a subscript, parallel to $L \to L_a$ in Chapter 5 where the specialisation $\mathcal{D}_{L_a}$ writes itself naturally. The corrected figure's rstick now writes $\mathcal{D}_L(\rho)$ (matching the body); the caption does too.
- **Macros used**: `\ee`, `\ii`, `\dd`, `\idm` (all defined at `main.tex:83–90`).
- **$L$ versus $L_a$**: the figure and preliminaries discuss a single abstract jump $L$; Chapter 5's $L_a = \sqrt{\gamma(\bar\omega)}\,\tilde{A}_a(\bar\omega)$ is the concrete instantiation. The subsection sticks to $L$ throughout to match the figure caption.
- **"First Chen citation" pedagogical pointer**: the $M_0/M_1$ display carries a footnote-style tag `[CHEN23-EQ3.2]` in the draft, which expands in the `.tex` to `\cite[Eq.~(3.2)]{chen2023quantum}` (or the biblatex equivalent). This gives the reader a specific equation handle rather than a fifty-page-paper citation.

### Length

5 short body paragraphs + 1 generality paragraph + 3 display equations (one labelled `eq:weak-meas-euler`). Rendered length is approximately 1 page at the thesis font size — still the shortest of the three preliminaries circuit subsections, as requested. The new generality paragraph is short (one paragraph); it replaces the prominent `[CHECK — important]` callout in the v2 Suggestions section, so the overall length is roughly flat.

### Forward link

The subsection closes its body with an explicit pointer to (Alg. alg:diss), which at `2_methods.tex:1744` is labelled `alg:diss` and whose commentary at `2_methods.tex:1884–1931` already uses the coherent-undo pattern established here. The new generality paragraph adds a forward pointer to `sec:algorithm` and to Alg. 1b steps 8 and 9–14.

### Figure-replacement note

The thesis figure `circ:weak-measurement` at `1_preliminaries.tex:1117–1167` should be *replaced* by the contents of `drafts/circuits/weak-measurement-v2.tex` (verbatim, minus the standalone-document preamble). The fixes included in the new file are:

1. **Block-register open control on $Y_\delta$** — `\octrl{-1}` on row 2 at the $Y_\delta$ column (new).
2. **Caption typo fix** — $\gamma = \delta/(\alpha^{2}\,\dd t)$ (was $\delta^{2}/(\alpha^{2}\,\dd t)$).
3. **Rstick notation** — $\mathcal{D}_L(\rho)$ (was $\mathcal{D}(\rho)$), matching the body's dissipator-vs-discriminant disambiguation.
4. **Caption sentence on the control** — one sentence added to the caption documenting that $Y_\delta$ is block-controlled on $\ket{0^b}$.

No other changes to the drawing; `row sep`, `column sep`, bundling conventions, and label placement are identical to the thesis's original.

### Suggestions

1. **Rstick alignment between figure and body**: the figure's rstick uses $\mathcal{D}_L(\rho)$ per the corrected file; the body uses $\mathcal{D}_L(\rho)$ throughout; and the body's opening line defines $\mathcal{D}_L$ explicitly. These are all aligned in v3. No action required beyond adopting the new figure.

2. **Multi-jump extension deferred**: the generalisation $\mathcal{L} = \sum_a \mathcal{D}_{L_a}$ is stated in one forward-pointing sentence at the end of §"One Euler step", not developed here. Developing it would duplicate Chapter 5's Alg. 1b commentary verbatim; defer stands.

3. **Optional label rename `eq:weak-meas-euler` → `eq:weak-meas`**: to match the naming style of `eq:qpe-unified` and `eq:B-block-encoding`. Both labels are fine; the longer one flags the Euler reading explicitly. Low priority.

4. **Optional link to classical Metropolis-Hastings** (not adopted in v3): the two ancilla outcomes $\ket{0}_\delta, \ket{1}_\delta$ map to MH's accept/propose; a one-sentence remark would reinforce the pedagogical arc from classical MCMC to quantum dissipation. Chapter 5's Alg. 1b commentary already makes this parallel, so duplicating it here would add words without new information.

5. **Optional reference to Lin 2025 / Fang et al. 2026** (not adopted in v3): these works apply the weak-measurement primitive to dissipative ground-state preparation with detectability-lemma arguments. Adding them here would place the construction in the broader landscape but dilute the preliminaries' focus; their natural home is Chapter 5 or a literature-review section.

---

## Review-driven choices applied (v2 → v3)

The following changes were made from `drafts/weak-measurement-subsection-v2.md` to v3, prompted by the user's directive to remove the outstanding figure-fix callout and to add a generality paragraph:

- **Removed** the `[CHECK — important]` Suggestion 1 from the v2 Suggestions section. The underlying issue — the thesis figure draws $Y_\delta$ as unconditional, contrary to Chen 2023 Thm III.1 — is now addressed by providing the corrected figure as `drafts/circuits/weak-measurement-v2.tex`, and a short note in the **Figure-replacement note** of the Writing notes lists the three one-line fixes.
- **Removed** the v2 Suggestions 2 and 3 (the $\delta^2 \to \delta$ caption fix and the $\mathcal{D}(\rho) \to \mathcal{D}_L(\rho)$ rstick fix), because both are now included in the new figure file and thus no longer pending author action.
- **Retained** the body's v2 prose for steps (i)–(iv), the Kraus display, and the Euler-step paragraph — no changes, because v2's body already describes the correct block-controlled $Y_\delta$ construction.
- **Tightened** step (ii) of the Circuit walkthrough: added the explicit second half of the $Y_\delta$ action — *"and carries $\ket{0}_\delta$ untouched on the leaked $\ket{\phi^\perp}_{bS}$ branch"* — to make the two-branch split visible, and wrote the rotation as two displayed identities (one for the block branch, one for the leaked branch) rather than one display. This closes a small gap between step (i)'s introduction of $\ket{\phi^\perp}_{bS}$ and step (iii)'s re-entry.
- **Renamed** step (ii)'s section heading from "Ancilla rotation" to "Ancilla rotation, controlled on the block", to foreground the control structure that is the subject of the new figure fix.
- **Added** §"Generality of the control structure" as a new short paragraph at the end of the body, before the Citations section. The paragraph situates the single-block-register $\ket{0^b}$ control within the richer $(\Omega, q_\gamma)$ ancilla structure of Chapter 5's Alg. 1b, explaining that $q_\gamma = \ket{0}_\gamma$ plays the role of $\ket{0^b}$ in the generic construction. This answers the user's question about how the same control mechanism generalises beyond the minimal LCU case.
- **Clarified** in the Writing notes' style section that the block-control is drawn as `\octrl` (open control), explaining the notational difference from Chen 2023 Fig. 3's filled-disk drawing (they label the state $\ket{0^b}$; we use the thesis's open-dot convention).
- **Updated** the insertion-target header to flag the figure replacement explicitly.
- **Kept** all v2 notation decisions (single-index $L$, `chen2023efficient` bibkey, $\mathcal{D}_L$ symbol, macros, pedagogical `[CHEN23-EQ3.2]` pointer, `eq:weak-meas-euler` label).
- **Kept** all v2 exposition choices (why-before-what opening, $\gamma_{\mathrm{eff}}$ naming, no-OAA-needed closure, "leaked projection becomes the Lindblad drain" punchline, collapsed $M_0/M_1$ display, $\|L/\alpha\|\leq 1$ subnormalisation parenthetical, (Eq. eq:kraus) pointer).
