# Weak-measurement

> **Insertion target:** `supplementary-informations/1_preliminaries.tex`, body of `\subsection{Weak-measurement.} \label{sec:prelim-weak-meas}` — text goes between the `\subsection{...}\label{...}` line at line 1116 and `\begin{figure}[ht]` at line 1117 (the figure `circ:weak-measurement` is already in place at lines 1117–1167).
> **Depends on labels already in the chapter:** `sec:prelim-LCU`, `sec:prelim-QPE`, `circ:weak-measurement`, and the $Y_\theta := R_Y(2\arcsin\sqrt\theta)$ shorthand introduced at `1_preliminaries.tex:825`.
> **Forward reference:** `alg:diss` (Algorithm 1b, `\textsc{DissipativeStep}`) in Chapter 5 at `2_methods.tex:1744`.
> **Bib keys used, already in `references.bib`:** `chen2023quantum`, `chen2023efficient`. No new keys introduced.

---

The Linear Combination of Unitaries of §sec:prelim-LCU block-encodes a non-unitary operator $L$ into the $\ket{0^b}$-block of a unitary $U_L$, and so delivers the amplitude $L\ket\psi/\alpha$ on the system register $S$ after post-selecting the block register $b$ on $\ket{0^b}$. What it does *not* deliver is the accompanying anti-Hermitian drain $-\tfrac{1}{2}\{L^\dagger L,\rho\}$ that promotes the bare action $\rho \mapsto L\rho L^\dagger$ into a trace-preserving Lindblad dissipator 
$$\mathcal{D}[L](\rho) \;:=\; L\rho L^\dagger - \tfrac{1}{2}\bigl\{L^\dagger L,\,\rho\bigr\}.$$
A single additional ancilla --- the *weak-measurement qubit* $\ket{0}_\delta$ --- closes this gap. The coherent-undo construction of [Chen et al. 2023] [CITE: chen2023quantum] threads $U_L$, a $Y_\delta$-rotation on $\ket{0}_\delta$, and an open-controlled $U_L^\dagger$ in such a way that post-selecting $b$ on $\ket{0^b}$ and measuring $\ket{0}_\delta$ implements exactly one Euler step of the semigroup $\ee^{\delta\mathcal{D}[L]/\alpha^2}$ on $\rho$, up to $\mathcal{O}(\delta^2)$ (Fig. circ:weak-measurement). It is this primitive that appears as the dissipative step of Chapter 5's Lindblad simulation (Alg. alg:diss).

*Circuit walkthrough.* Initialise the three registers in $\ket{0}_\delta \otimes \ket{0^b} \otimes \rho$ and apply the four operations of Fig. circ:weak-measurement in turn. (i) The block encoding $U_L$ acts on $(b,S)$; projected onto $\ket{0^b}$ the amplitude on $S$ is $L/\alpha$, so on the mixed state this imprints the unnormalised Kraus component $(L/\alpha)\,\rho\,(L^\dagger/\alpha)$ on the $\ket{0^b}$-branch. (ii) The rotation $Y_\delta$ on $\ket{0}_\delta$, using the thesis shorthand $Y_\theta := R_Y(2\arcsin\sqrt\theta)$ of `1_preliminaries.tex:825`, acts as
$$Y_\delta\ket{0}_\delta \;=\; \sqrt{1-\delta}\,\ket{0}_\delta \;+\; \sqrt{\delta}\,\ket{1}_\delta,$$
so that the joint state splits into two branches of amplitudes $\sqrt{1-\delta}$ and $\sqrt\delta$. (iii) The open-controlled $U_L^\dagger$ on $(b,S)$ fires only on the $\ket{0}_\delta$ branch: there it cleanly uncomputes $U_L$ and restores $\rho$ on $S$ with the ancillas returned to $\ket{0^b}\ket{0}_\delta$; on the $\ket{1}_\delta$ branch $U_L$ stands, and the system still carries the amplitude $L\rho L^\dagger/\alpha^2$ inside the $\ket{0^b}$-component of $(b,S)$. (iv) Post-select $b$ on $\ket{0^b}$ and measure $\ket{0}_\delta$.

*Kraus operators.* The two branches of the ancilla measurement give two Kraus operators on $S$. In the *no-jump* branch ($\ket{0}_\delta$ outcome, $\ket{0^b}$ post-selected), the uncomputation $U_L^\dagger U_L = \idm$ acts inside the $\ket{0^b}$-block, and a standard amplitude-amplification bookkeeping [CITE: chen2023quantum] shows
$$M_0 \;=\; \sqrt{1-\delta}\;\idm \;-\; \tfrac{1}{2}\,\tfrac{\delta}{\alpha^2}\,L^\dagger L \;+\; \mathcal{O}(\delta^2),$$
where the second term comes from the failed post-selection on $b$ when $U_L^\dagger$ does *not* act — it is $-(\delta/(2\alpha^2))L^\dagger L$, precisely the anti-Hermitian drain of $\mathcal{D}[L]$. In the *jump* branch ($\ket{1}_\delta$ outcome, $\ket{0^b}$ post-selected), $U_L^\dagger$ is absent on $S$ and the full amplitude $(L/\alpha)\rho(L^\dagger/\alpha)$ survives, weighted by the $\sqrt\delta$ of the rotation:
$$M_1 \;=\; \sqrt{\delta}\;\frac{L}{\alpha}.$$

*One Euler step.* Summing the two branches and expanding in $\delta$,
$$\rho \;\mapsto\; M_0\,\rho\,M_0^\dagger + M_1\,\rho\,M_1^\dagger \;=\; \rho \;+\; \frac{\delta}{\alpha^{2}}\Bigl(L\rho L^\dagger - \tfrac{1}{2}\{L^\dagger L,\rho\}\Bigr) \;+\; \mathcal{O}(\delta^{2}),$$ <!-- \label{eq:weak-meas-euler} -->
which is one Euler step of $\dot\rho = (1/\alpha^2)\,\mathcal{D}[L](\rho)$ — the rstick of Fig. circ:weak-measurement. Identifying $\delta$ as the Euler step-size of the semigroup then gives the effective Lindblad rate $\gamma = \delta/(\alpha^2\,\dd t)$ (see [CHECK] in the notes: the figure caption writes $\gamma = \delta^2/(\alpha^2\,\dd t)$, which corresponds to the diffusive-scaling reading $\delta \sim \sqrt{\dd t}$ instead). The post-selection on $b$ succeeds with probability $1 - \delta\,\|L\ket\psi\|^2/\alpha^2 + \mathcal{O}(\delta^2) = 1 - \mathcal{O}(\delta)$, so the channel is near-deterministic per step; the entire construction trades one failure mode of the LCU block-encoding (the $1-1/\alpha^2$ rejection probability) for the dissipator that the simulation needs. Full derivation of (Eq. eq:weak-meas-euler), including the extension to the multi-jump sum $\mathcal{L} = \sum_a \mathcal{D}[L_a]$ and to the KMS-weighted filtered jumps of [Chen et al. 2025] [CITE: chen2023efficient] that replace $L$ in Chapter 5, is deferred to the commentary around (Alg. alg:diss).

---

## Citations

Both keys already exist in `references.bib` — no new bibtex entries are introduced.

- **`chen2023quantum`** — Chen, C.-F., Kastoryano, M. J., Brandão, F. G. S. L., Gilyén, A., *"Quantum Thermal State Preparation"*, arXiv:2303.18224 (2023). Theorem III.1 establishes the weak-measurement implementation of a Lindblad dissipator from a block-encoded jump; the four-step circuit of Fig. circ:weak-measurement is the single-jump restriction of their construction.
- **`chen2023efficient`** — Chen, C.-F., Kastoryano, M. J., Gilyén, A., *"An Efficient and Exact Noncommutative Quantum Gibbs Sampler"*, arXiv:2311.09207 (2023). The same weak-measurement primitive survives in the KMS-exact variant; only the jumps fed into $U_L$ change.

---

## Writing notes

<!-- For the author, not for the thesis -->

### [CHECK] items

- **Caption vs rstick: $\delta^2$ vs $\delta$.** Fig. circ:weak-measurement's rstick reads $\rho + \tfrac{\delta}{\alpha^{2}}\,\mathcal{D}(\rho) + \mathcal{O}(\delta^{2})$ (linear in $\delta$, consistent with the $Y_\theta := R_Y(2\arcsin\sqrt\theta)$ convention used throughout the thesis and with Chapter 5's Euler expansion $\ee^{\delta\mathcal{L}_{a,\mathrm{diss}}} = \mathrm{id} + \delta\mathcal{L}_{a,\mathrm{diss}} + \mathcal{O}(\delta^2)$ at `2_methods.tex:1910,1929`). The caption, however, writes the rate as $\gamma = \delta^2/(\alpha^2\,\dd t)$. The two are reconcilable under *different* scalings of $\delta$ against $\dd t$:
  - *Euler reading* ($\delta \sim \dd t$, i.e. $\delta$ is the step-size directly): the per-step Kraus sum gives $\rho + (\delta/\alpha^2)\mathcal{D}(\rho) + \mathcal{O}(\delta^2)$, so matching $\rho + \gamma\,\dd t\,\mathcal{D}(\rho)$ forces $\gamma = \delta/(\alpha^2\,\dd t)$.
  - *Diffusive reading* ($\delta \sim \sqrt{\dd t}$, i.e. the ancilla rotation angle is $\sqrt{\delta}$ which is $\propto (\dd t)^{1/4}$): then $\delta^2 \sim \dd t$ and $\gamma = \delta^2/(\alpha^2\,\dd t)$ is the physical rate, matching the caption, but the rstick should then read $(\delta^2/\alpha^2)\mathcal{D}(\rho)$.

  The thesis's Chapter 5 usage is unambiguously the *Euler reading* (step 8 of `alg:diss` rotates by $R_Y(2\arcsin\sqrt\delta)$, and the channel has prefactor $\delta$ not $\delta^2$). I therefore suspect the caption is the typo and should read $\gamma = \delta/(\alpha^2\,\dd t)$. I have written the body text in the Euler convention (matching the rstick and Chapter 5), and flagged the caption discrepancy inline. Please confirm the intended reading and fix the caption accordingly.

- **Multi-jump extension.** I deliberately deferred the multi-jump generalization $\mathcal{L}(\rho) = \sum_a \mathcal{D}[L_a](\rho)$ and the KMS-weighted filtered jumps $\tilde{A}_a(\bar\omega)$ to Chapter 5, where they are the natural next step (the Boltzmann filter and OFT-register dressing of `alg:diss` add nothing new at the level of Fig. circ:weak-measurement's mechanism). A single forward-pointing sentence at the end of (Eq. eq:weak-meas-euler)'s paragraph makes the connection explicit. If you want an additional half-sentence in the preliminaries remarking that one block-encodes $\sum_a L_a/\alpha$ instead of $L/\alpha$ by an extra LCU layer, I can add it, but it would duplicate material that Chapter 5 already has — recommend leaving it out.

- **$M_0$ and the failed-post-select correction.** The formula $M_0 = \sqrt{1-\delta}\,\idm - (\delta/(2\alpha^2))L^\dagger L + \mathcal{O}(\delta^2)$ is the standard amplitude-bookkeeping result: the $\sqrt{1-\delta}\,\idm$ comes from the uncomputed $\ket{0}_\delta$-branch, and the $-(\delta/(2\alpha^2))L^\dagger L$ contribution is the leakage from the $\ket{1}_\delta$-branch that also has amplitude in $\ket{0^b}$ (i.e. the $L^\dagger L/\alpha^2$-weighted rejection that, after Taylor-expanding the post-selection probability, reconstitutes the anti-Hermitian drain). I have given this as a one-line statement with a Chen CKBG citation rather than expanding the bookkeeping, because (a) the full derivation is in [Chen et al. 2023, Thm III.1] and (b) the subsection is the last of three short preliminaries subsections, not a chapter. If a one-paragraph derivation is preferred, it can be inserted without changing the rest of the structure.

### Style / notation decisions

- **`\ket{0}_\delta`** single-qubit label, **`\ket{0^b}`** bundled block-register label: verbatim from the circuit figure's `\lstick` labels at lines 1125 and 1139. The external circuit file `drafts/circuits/weak-measurement.tex` uses `q_\delta` for the ancilla label, but the in-thesis figure uses `\ket{0}_\delta` — I followed the thesis figure.
- **$Y_\delta$**: used throughout as the shorthand $R_Y(2\arcsin\sqrt\delta)$ per `1_preliminaries.tex:825`, matching the thesis convention. Crucially the amplitude imprinted on $\ket{1}_\delta$ is $\sqrt\delta$, not $\sin(\delta/2)\approx\delta/2$ — this is why the dissipator prefactor is $\delta$ (not $\delta^2$) in the rstick. I made the convention explicit inline when $Y_\delta$ is first used, with a pointer back to line 825.
- **Macros used**: `\ee`, `\ii`, `\dd`, `\idm` (all defined at `main.tex:83–90`).
- **$L$ versus $L_a$**: the figure and preliminaries discuss a single abstract jump $L$; Chapter 5's $L_a = \sqrt{\gamma(\bar\omega)}\,\tilde{A}_a(\bar\omega)$ is the concrete instantiation. The subsection sticks to $L$ throughout to match the figure caption.

### Length

6 short paragraphs + 3 display equations (two unlabelled, one labelled `eq:weak-meas-euler`). Final rendered length is comfortably under the 1.5-page cap — closer to 1 page once set in the thesis font size, consistent with the user's ask for the *shortest* of the three circuit subsections. If the extra Kraus-operator paragraph feels heavier than a preliminary needs, the $M_0$ and $M_1$ display equations can be merged into one line: `$M_0 = \sqrt{1-\delta}\,\idm + \mathcal{O}(\delta/\alpha^2),\;M_1 = \sqrt\delta\,L/\alpha$` without changing the derivation that follows.

### Forward link

The subsection closes with an explicit pointer to (Alg. alg:diss), which at `2_methods.tex:1744` is labelled `alg:diss` and whose in-text commentary at `2_methods.tex:1884–1931` already uses exactly the coherent-undo pattern introduced here. The wording "deferred to the commentary around (Alg. alg:diss)" matches the cross-reference style used in the QPE subsection (`drafts/qpe-subsection.md`, last paragraph).

### Suggestions (rounding out the subsection)

1. **Fix the caption $\delta^2 \to \delta$** in Fig. circ:weak-measurement to match the rstick and Chapter 5 (see [CHECK] above). This is a one-character edit at `1_preliminaries.tex:1165` and removes a genuine internal inconsistency.
2. **Optionally mention $\|L/\alpha\| \leq 1$** (block-encoding subnormalisation) in the opening paragraph, since it is the reason the expansion in $\delta$ converges uniformly and the $\mathcal{O}(\delta^2)$-error is meaningful. I chose to leave this implicit because §sec:prelim-LCU already defines the subnormalisation $\alpha$, but a six-word parenthetical would not hurt.
3. **Optionally swap the labelled display to `eq:weak-meas`** or similar, to match the naming of sibling labels (`eq:qpe-unified` in QPE, `eq:B-block-encoding` in Chapter 5). The current choice `eq:weak-meas-euler` emphasises that this is specifically the Euler-step identity, which I think is informative, but either works.
4. **If a future revision adds a figure for the full $\textsc{DissipativeStep}$** at the end of the preliminaries (so that the preliminaries chapter contains all three primitives *and* their composition), the weak-measurement subsection could close with a one-sentence "composed with the OFT of §sec:prelim-QPE this gives the full KMS-weighted dissipator $\ee^{\delta\mathcal{L}_{a,\mathrm{diss}}}$" forward link. I chose not to include that because the current preliminaries chapter stops at the three individual primitives and defers the composition to Chapter 5; adding it would shift the division of labour between chapters. Flagging it for consideration.
