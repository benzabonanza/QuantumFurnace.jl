# Weak-measurement (v4)

> **Insertion target:** `supplementary-informations/1_preliminaries.tex`, body of `\subsection{Weak-measurement.} \label{sec:prelim-weak-meas}` — text goes between the `\subsection{...}\label{...}` line at line 1116 and `\begin{figure}[ht]` at line 1117 (the figure `circ:weak-measurement` sits at lines 1117–1167).
> **Figure replacement:** the thesis figure at `1_preliminaries.tex:1117–1167` must be replaced with `drafts/circuits/weak-measurement-v2.tex` — the current figure draws $Y_\delta$ as unconditional (inconsistent with Chen 2023 Thm III.1, the construction the body follows), and the new file adds the missing block-register open control and fixes the caption typo.
> **Depends on labels already in the chapter:** `sec:prelim-LCU`, `sec:prelim-QPE`, `eq:kraus` (Kraus form of CPTP maps at `1_preliminaries.tex:47`-ish), `circ:weak-measurement`, and the $Y_\theta := R_Y(2\arcsin\sqrt\theta)$ shorthand at `1_preliminaries.tex:825`.
> **Forward references:** `alg:diss` (the `DissipativeStep` algorithm) in Chapter 5 at `2_methods.tex:1744`, and `sec:algorithm` at `2_methods.tex:1620`-ish.
> **Bib keys used, already in `references.bib`:** `chen2023quantum`, `chen2023efficient`. **New keys introduced in v4:** `liwang2023`, `ding2024efficient` — bibtex stubs at the bottom.

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

## Choice of Lindblad simulator: continuous vs discrete jumps

The Euler-step accuracy of (Eq. eq:weak-meas-euler) is only first order in $\delta$, so reaching diamond-norm error $\varepsilon$ over total time $T$ costs $T/\delta = \widetilde{\mathcal{O}}(T^2/\varepsilon)$ calls to $U_L$ — polynomial, not poly-logarithmic, in $1/\varepsilon$. One might therefore ask why we do not instead pick a higher-order Lindblad simulator, such as the Duhamel-principle / scaled-Gaussian-quadrature scheme of [Li and Wang 2023] [CITE: liwang2023], which simulates a generator $\mathcal{L}(\rho) = -\ii[H,\rho] + \sum_{j=1}^m \mathcal{D}_{L_j}(\rho)$ with $\widetilde{\mathcal{O}}(t \,\mathrm{polylog}(1/\varepsilon))$ block-encoding queries (Theorem 1, with $\tau := t \|\mathcal{L}\|_{\mathrm{be}}$). The constraint hidden in the "$j = 1, \dots, m$" of that theorem is the answer: higher-order series simulators ingest a *finite* list $\{L_j\}_{j=1}^m$ of jump operators, each supplied as its own block-encoding $U_{L_j}$, and their cost is linear in $m$. They do not absorb a generator of the form

$$\mathcal{L}(\rho) \;=\; \sum_{a \in \mathcal{A}} \int \dd\omega\;\gamma(\omega)\,\mathcal{D}_{\hat A_a(\omega)}(\rho),$$

whose jump operators $\hat A_a(\omega)$ are *continuously parametrized* by a Bohr-frequency variable $\omega \in \mathbb{R}$ — even after a fine quadrature discretisation, $m$ is at least the size of the frequency grid, and each grid point would need its own block encoding rather than being addressed coherently through the $\Omega$ register built up in §sec:prelim-QPE. The weak-measurement primitive is precisely what one uses when the integral cannot be peeled off the block encoding: a single $U_L$ call already block-encodes the full OFT-coherent jump $\sum_a \int \dd\omega\,\sqrt{\gamma(\omega)}\,\hat A_a(\omega)$, and one Euler step of (Eq. eq:weak-meas-euler) realises one step of the *whole* Lindbladian per call, with cost that scales with the resolution of the OFT rather than with the number of effective jumps.

The trade-off between the two routes is therefore architectural rather than numerical. The higher-order series scheme buys $\mathrm{polylog}(1/\varepsilon)$-accurate channel-as-LCU implementation, but only for jump sets that have been collapsed to a small finite list at compile time; typical workarounds for the continuous-$\omega$ case — e.g. replacing $\gamma(\omega)$ with a discrete sum of $\delta$-functions so that the integral becomes a finite sum [CITE: ding2024efficient] — give back exactly the finite-$m$ regime in which the higher-order route applies, at the cost of working with a different jump ensemble. The weak-measurement Euler step, conversely, accepts an arbitrary block-encoded $L$ — including the OFT-coherent superposition over $\omega$ — in exchange for first-order accuracy in $\delta$. Since the construction the rest of this thesis follows keeps $\gamma(\omega)$ continuous (Chapter 5; the OFT machinery of §sec:prelim-QPE is only useful if it is exercised), the weak-measurement primitive is the only one of the two that applies.

## Generality of the control structure

The figure above depicts the *minimal* LCU-block-encoded jump: a single block register $\ket{0^b}$ plays the role of "$U_L$ succeeded", and $Y_\delta$ is open-controlled on $\ket{0^b}$ to restrict the weak rotation to the post-selection branch. Chapter 5's dissipative step (Alg. alg:diss, §sec:algorithm) replaces this minimal block register with a richer ancilla structure: the frequency register $\Omega$ and the Boltzmann qubit $q_\gamma$ together encode the KMS-weighted jump $\sqrt{\gamma(\bar\omega)}\,\tilde{A}_a(\bar\omega)$ on the system, and $Y_\delta$ is then open-controlled on $q_\gamma = \ket{0}_\gamma$ in place of the generic $\ket{0^b}$ (Alg. alg:diss step 8, `2_methods.tex:1789`). The mechanism is the same — $Y_\delta$ always fires on the ancilla branch flagging successful application of the effective jump — and the Euler-step identity (Eq. eq:weak-meas-euler) persists with $L$ replaced by the OFT-filtered jump $\sqrt{\gamma(\bar\omega)}\,\tilde{A}_a(\bar\omega)$. In Chapter 5's picture, $\Omega$ is not post-selected but coherently uncomputed together with $q_\gamma$ on the no-jump branch (Alg. alg:diss steps 9–14), so the role of the block register's post-selection is taken over by the controlled uncomputation $U_{\mathrm{diss}}^\dagger$; the $Y_\delta$-control structure is inherited verbatim.

---

## Citations

`chen2023quantum` and `chen2023efficient` already exist in `references.bib` — no new bibtex entries for them.

- **`chen2023quantum`** — Chen, C.-F., Kastoryano, M. J., Brandão, F. G. S. L., Gilyén, A., *"Quantum Thermal State Preparation"*, arXiv:2303.18224 (2023). Theorem III.1 and Eq. (3.2) in their proof establish the weak-measurement implementation of a Lindblad dissipator from a block-encoded jump; the circuit of Fig. circ:weak-measurement is the single-jump ($|J|=1$) restriction of their construction, and their Fig. 3 is the canonical drawing of the block-controlled $Y_\delta$ convention.
- **`chen2023efficient`** — Chen, C.-F., Kastoryano, M. J., Gilyén, A., *"An Efficient and Exact Noncommutative Quantum Gibbs Sampler"*, arXiv:2311.09207 (2023). The same weak-measurement primitive survives in the KMS-exact variant; only the jumps fed into $U_L$ change.

**New for v4:**

- **`liwang2023`** — Li, X., Wang, C., *"Simulating Markovian open quantum systems using higher order series expansion"*, arXiv:2212.02051 (2023). Theorem 1: given block-encodings $U_H$ of $H$ and $U_{L_j}$ of $L_j$ ($j=1,\dots,m$), simulates $\ee^{\mathcal{L} t}$ in diamond norm to error $\varepsilon$ using $\widetilde{\mathcal{O}}(\tau\, \mathrm{polylog}(t/\varepsilon))$ block-encoding queries with $\tau = t \|\mathcal{L}\|_{\mathrm{be}}$, via Duhamel's principle and scaled Gaussian quadrature.
- **`ding2024efficient`** — Ding, Z., Li, B., Lin, L., *"Efficient quantum Gibbs samplers with Kubo–Martin–Schwinger detailed balance condition"*, arXiv:2404.05998 (2024). Constructs a family of KMS-DB Lindbladians whose Kossakowski function $\gamma(\omega)$ may be a *discrete* sum of $\delta$ functions, yielding a finite jump set that can then be simulated by any high-order Lindblad simulator (cf. their Sec. 1.1, p. 2 last paragraph and Sec. 1.2, p. 3 first paragraph).

```bibtex
@article{liwang2023,
  author        = {Li, Xiantao and Wang, Chunhao},
  title         = {Simulating {M}arkovian open quantum systems using higher order series expansion},
  journal       = {arXiv preprint},
  year          = {2023},
  eprint        = {2212.02051},
  archivePrefix = {arXiv},
  primaryClass  = {quant-ph},
  url           = {https://arxiv.org/abs/2212.02051},
  note          = {Theorem~1: $\widetilde{\mathcal{O}}(\tau\,\mathrm{polylog}(t/\varepsilon))$ queries to block-encodings of $H$ and finitely many $L_j$.}
}

@article{ding2024efficient,
  author        = {Ding, Zhiyan and Li, Bowen and Lin, Lin},
  title         = {Efficient quantum {G}ibbs samplers with {K}ubo--{M}artin--{S}chwinger detailed balance condition},
  journal       = {arXiv preprint},
  year          = {2024},
  eprint        = {2404.05998},
  archivePrefix = {arXiv},
  primaryClass  = {quant-ph},
  url           = {https://arxiv.org/abs/2404.05998},
  note          = {KMS-DB Lindbladian family; $\gamma(\omega)$ may be a discrete sum of $\delta$-functions, yielding finitely many jump operators amenable to high-order Lindblad simulation.}
}
```

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
- **Continuous-vs-discrete language in the new §"Choice of Lindblad simulator"**: the new section deliberately stays in fully general Lindblad-simulation terms ("continuously parametrized jumps", "finite list of $L_j$", "discrete sum of $\delta$-functions"). It does not name "Chen 2025 KMS" or "Ding–Li–Lin" inside the body — these names appear only in the Citations entries — because at the position of the prelims subsection the reader has not yet been introduced to either Lindbladian variant. The two BibKeys it does use (`liwang2023` for the high-order simulator, `ding2024efficient` for the discrete-jump option) are matched to the referenced *technique*, not to any thesis-bookkeeping forward pointer.

### Length

5 short body paragraphs + 1 generality paragraph + 1 new "Choice of simulator" two-paragraph section + 4 display equations (one labelled `eq:weak-meas-euler`, one is the unlabelled $\sum_a \int \dd\omega$ Lindbladian display in the new section). Rendered length is approximately 1.3 pages at the thesis font size — still the shortest of the three preliminaries circuit subsections, as the new section is intentionally compact (two paragraphs and one display).

### Forward link

The subsection closes its body with an explicit pointer to (Alg. alg:diss), which at `2_methods.tex:1744` is labelled `alg:diss` and whose commentary at `2_methods.tex:1884–1931` already uses the coherent-undo pattern established here. The new generality paragraph adds a forward pointer to `sec:algorithm` and to Alg. 1b steps 8 and 9–14. The new §"Choice of Lindblad simulator" deliberately does *not* add a forward pointer — it stands as a self-contained technical observation about block-encoding granularity.

### Figure-replacement note

The thesis figure `circ:weak-measurement` at `1_preliminaries.tex:1117–1167` should be *replaced* by the contents of `drafts/circuits/weak-measurement-v2.tex` (verbatim, minus the standalone-document preamble). The fixes included in the new file are:

1. **Block-register open control on $Y_\delta$** — `\octrl{-1}` on row 2 at the $Y_\delta$ column (new).
2. **Caption typo fix** — $\gamma = \delta/(\alpha^{2}\,\dd t)$ (was $\delta^{2}/(\alpha^{2}\,\dd t)$).
3. **Rstick notation** — $\mathcal{D}_L(\rho)$ (was $\mathcal{D}(\rho)$), matching the body's dissipator-vs-discriminant disambiguation.
4. **Caption sentence on the control** — one sentence added to the caption documenting that $Y_\delta$ is block-controlled on $\ket{0^b}$.

No other changes to the drawing; `row sep`, `column sep`, bundling conventions, and label placement are identical to the thesis's original.

### Suggestions

1. **Rstick alignment between figure and body**: the figure's rstick uses $\mathcal{D}_L(\rho)$ per the corrected file; the body uses $\mathcal{D}_L(\rho)$ throughout; and the body's opening line defines $\mathcal{D}_L$ explicitly. These are all aligned in v3/v4. No action required beyond adopting the new figure.

2. **Multi-jump extension deferred**: the generalisation $\mathcal{L} = \sum_a \mathcal{D}_{L_a}$ is stated in one forward-pointing sentence at the end of §"One Euler step", not developed here. Developing it would duplicate Chapter 5's Alg. 1b commentary verbatim; defer stands.

3. **Optional label rename `eq:weak-meas-euler` → `eq:weak-meas`**: to match the naming style of `eq:qpe-unified` and `eq:B-block-encoding`. Both labels are fine; the longer one flags the Euler reading explicitly. Low priority.

4. **Optional link to classical Metropolis-Hastings** (not adopted in v3/v4): the two ancilla outcomes $\ket{0}_\delta, \ket{1}_\delta$ map to MH's accept/propose; a one-sentence remark would reinforce the pedagogical arc from classical MCMC to quantum dissipation. Chapter 5's Alg. 1b commentary already makes this parallel, so duplicating it here would add words without new information.

5. **Optional reference to Lin 2025 / Fang et al. 2026** (not adopted in v3/v4): these works apply the weak-measurement primitive to dissipative ground-state preparation with detectability-lemma arguments. Adding them here would place the construction in the broader landscape but dilute the preliminaries' focus; their natural home is Chapter 5 or a literature-review section.

6. **Cost-comparison table** (not adopted in v4): the new §"Choice of Lindblad simulator" could be tightened into a two-row table (Euler vs Li–Wang, columns: jump-set type / $\varepsilon$-scaling / block-encoding count). The two-paragraph prose form was kept to match the surrounding declarative voice; a table would jut out visually as the only such object in the subsection. Re-evaluate if the new section grows.

---

## Review-driven choices applied (v3 → v4)

The following changes were made from `drafts/weak-measurement-subsection-v3.md` to v4, prompted by the user's directive to add a self-contained comparison between the weak-measurement primitive and a higher-order Lindblad simulator:

- **Added** a new short section §"Choice of Lindblad simulator: continuous vs discrete jumps", inserted between §"One Euler step" and §"Generality of the control structure". Two paragraphs + one display equation. The section makes two technical points: (a) the weak-measurement Euler step is the only way to feed an OFT-coherent, continuously $\omega$-parametrized jump set into a Lindblad simulation via a *single* block encoding $U_L$, with cost scaling in the OFT resolution rather than in the (effectively infinite) jump count; (b) when the jump set is finite, a higher-order Duhamel-principle / scaled-Gaussian-quadrature simulator [CITE: liwang2023] gives $\widetilde{\mathcal{O}}(t\,\mathrm{polylog}(1/\varepsilon))$ scaling instead of the Euler step's polynomial $1/\varepsilon$, but it ingests the jumps as a finite list $\{L_j\}_{j=1}^m$ of separate block encodings and is linear in $m$. The discrete-$\delta$-function Kossakowski variant of [CITE: ding2024efficient] is the natural way to enter that finite-$m$ regime from a KMS-DB Lindbladian.
- **Wrote the new section in fully general terms** — "continuously parametrized jumps", "finite list of $L_j$", "discrete sum of $\delta$-functions" — and **never named Chen 2025 KMS or Ding–Li–Lin in the body**. This is a hard constraint from the user: the prelims subsection sits at line 1268 of `1_preliminaries.tex`, before either Lindbladian variant has been introduced. The two new bibkeys `liwang2023` and `ding2024efficient` appear *only* via `[CITE: ...]` annotations next to the relevant techniques (the high-order simulator and the discrete-$\delta$ Kossakowski option), matching the citation pattern already used in v3 (`chen2023quantum` for the weak-measurement primitive itself, `chen2023efficient` for the KMS-weighted jumps).
- **Verified the source claims against the papers before writing.** Li & Wang 2023 (`supplementary-informations/Li and Wang - Simulating Markovian open quantum systems using higher order series expansion.pdf`) Theorem 1 (informal, p. 4): $\widetilde{\mathcal{O}}(\tau\,\mathrm{polylog}(t/\varepsilon))$ queries to $U_H$ and to the block-encodings $U_{L_j}$, with $\tau := t\|\mathcal{L}\|_{\mathrm{be}}$ and an *additional* $\widetilde{\mathcal{O}}(m\tau\,\mathrm{polylog}(t/\varepsilon))$ 1- and 2-qubit gates — the linear-in-$m$ dependence is what motivates the "finite list" framing. Ding–Li–Lin 2024 KMS (`supplementary-informations/Ding et al. - 2024 - Efficient quantum Gibbs samplers with Kubo--Martin.pdf`) p. 2 last paragraph confirms the obstruction for the continuous-$\omega$ Chen 2025 case (*"To the best of our knowledge, high-order Lindblad simulators designed for a finite number of jump operators ... are not suitable for this task"*); p. 3 first paragraph confirms the discrete-$\delta$-function construction (*"$\gamma(\omega)$ can be chosen to be a discrete sum of $\delta$ functions, leading to a finite number of jump operators ... can be efficiently simulated using any high-order simulation algorithms, including those in [LW23, DLL24]"*). The new section's two technical claims and the citation pairing track these passages directly.
- **Updated the Insertion-target header** to flag the two new bibtex keys (`liwang2023`, `ding2024efficient`).
- **Added a Citations subsection** for the two new keys with full author/title/journal/arXiv data, plus a fenced `bibtex` block holding stub entries (the user can replace with canonical entries from `references.bib` upstream if needed).
- **Added a §Style note** to the Writing notes documenting the deliberate choice to keep the new section name-free (no Chen 2025, no Ding–Li–Lin in the body).
- **Added a §Length entry update** to reflect the new section's contribution (~0.3 page).
- **Added a §Forward link entry update** noting that the new section, unlike §"Generality of the control structure", deliberately omits any forward pointer, so as to read as a self-contained technical observation.
- **Added Suggestion 6** flagging the table form as a future option if the new section grows.
- **Kept everything else verbatim from v3**: the entire body of §"Circuit walkthrough", §"Kraus operators", §"One Euler step" up to its closing sentence, and §"Generality of the control structure" are character-for-character identical to v3. The v2 → v3 review-driven choices remain valid and are not repeated here.
