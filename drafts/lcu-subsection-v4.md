# Linear Combination of Unitaries (v4 — $\cos\theta = \lambda/\alpha$ convention, aligned with QSP v4)

> **Insertion target:** `supplementary-informations/1_preliminaries.tex`, body of `\subsection{Linear Combination of Unitaries.}\label{sec:prelim-LCU}`, between the `\subsection{...}\label{...}` line at 1037 and the `\begin{figure}` at 1039. The circuit figure `circ:lcu` (lines 1039–1064) stays immediately below the body.
> **Change from v3:** three `sin` → `cos` surgical fixes in the GQSP forward pointer (paragraph 5, the Citations entry for `motlagh2024generalized`, and the v2 → v3 change summary), to match QSP v4's Motlagh–Wiebe convention $\cos\theta = \lambda/\alpha$ for the qubitization walk $W = R_b\,U_B$. All other v3 content is preserved verbatim.
> **Existing labels referenced:** `circ:lcu`, `sec:prelim-QSP`, `eq:B-block-encoding`, `eq:b_minus`, `eq:b_plus-s-eta`, `circ:U_B_a`.
> **New labels introduced:** `eq:lcu-block-encoding`.
> **New citation keys introduced:** `lowChuang2019qubitization` (bibtex stub at the bottom). All other keys (`childs2012hamiltonian`, `berry2015simulating`, `chen2023efficient`, `motlagh2024generalized`) already appear in the thesis's `references.bib`.

---

Trotterization gives us Hamiltonian evolution; phase estimation gives us energy labels. The remaining circuit primitive needed for Chapter 5's coherent correction $B$ is a way to realise a *weighted sum of unitaries* as a single block-encoded unitary circuit, at a controlled cost in ancilla qubits and success probability. This is the *linear combination of unitaries* (LCU) scheme of [Childs and Wiebe 2012] [CITE: childs2012hamiltonian] and [Berry et al. 2015] [CITE: berry2015simulating], which we recall here in its modern block-encoding form [CITE: lowChuang2019qubitization] because that is the interface on which the GQSP subsection (see §sec:prelim-QSP) builds.

Fix a unitary family $\{U(t)\}_{t\in S}$ indexed by a discrete grid $S$ and a scalar weight $f : S \to \mathbb{R}_{\geq 0}$ with finite $\ell_1$ norm $\|f\|_1 = \sum_t f(t)$. The target operator
$$
A \;=\; \sum_{t \in S} f(t)\,U(t)
$$
is in general **Hermitian but not unitary** when the unitaries $U(t)$ come in conjugate pairs $U(-t) = U(t)^\dagger$ and $f$ is real and symmetric; its operator norm is at most $\|f\|_1$ by the triangle inequality and typically well below it. The LCU circuit in (Fig. circ:lcu) realises $A$ up to a known rescaling using one *block ancilla* register of $b = \lceil \log_2|S|\rceil$ qubits and three ingredients:

1. $\textsc{Prep}$ on the ancilla, loading the filter state $\ket{f} := \frac{1}{\sqrt{\|f\|_1}}\sum_{t \in S} \sqrt{f(t)}\,\ket{t}$;
2. $\textsc{Select}$ on system$\otimes$ancilla, the controlled family $\sum_{t}\ket{t}\!\bra{t}\otimes U(t)$;
3. $\textsc{Prep}^\dagger$, which uncomputes the ancilla.

Writing the composite circuit as $U_A := (\textsc{Prep}^\dagger \otimes \idm)\, \textsc{Select}\, (\textsc{Prep} \otimes \idm)$, a direct calculation shows that its $\ket{0^b}$–$\ket{0^b}$ matrix element on the ancilla is exactly
$$
\bra{0^b}\,U_A\,\ket{0^b} \;=\; \frac{A}{\|f\|_1},
$$ <!-- \label{eq:lcu-block-encoding} -->
i.e. $U_A$ is a *block encoding* of $A$ with subnormalisation $\alpha := \|f\|_1$ [CITE: lowChuang2019qubitization]. The subnormalisation is the price for packing a non-unitary $A$ inside a unitary: $U_A$ acts as $A/\alpha$ on the block projector $\ket{0^b}\!\bra{0^b}$ and as an unspecified completion on its orthogonal complement, and all downstream complexity statements (GQSP polynomial degree, simulation time, success probability) scale with $\alpha$ rather than with $\|A\|$.

If $f$ takes signed or complex values, two equivalent implementations are in common use. The first absorbs $\mathrm{sgn}(f(t)) := f(t)/|f(t)|$ into $\textsc{Select}$ by redefining $\tilde U(t) := \mathrm{sgn}(f(t))\,U(t)$ and letting $\textsc{Prep}$ load the magnitude state $\ket{|f|} := \frac{1}{\sqrt{\|f\|_1}}\sum_t \sqrt{|f(t)|}\,\ket{t}$. The second splits $\textsc{Prep}$ into an **asymmetric left/right pair** $\textsc{Prep}'$ / $\textsc{Prep}$: $\textsc{Prep}'$ writes the signed / phased amplitude $\sum_t \mathrm{sgn}(f(t))\sqrt{|f(t)|/\|f\|_1}\,\ket{t}$ on the left, and $\textsc{Prep}$ writes the magnitude-only amplitude $\sum_t \sqrt{|f(t)|/\|f\|_1}\,\ket{t}$ on the right, so that the two amplitudes multiply to $f(t)/\|f\|_1$. The subnormalisation $\alpha = \|f\|_1 = \sum_t |f(t)|$ is the same in both conventions. Chapter 5 adopts the asymmetric variant for the signed kernels $b_-$ (see \eqref{eq:b_minus}) and $b_+$ (see \eqref{eq:b_plus-s-eta}), where $\mathrm{sgn}(\cdot)$ is to be read as a unit-modulus phase in the complex case (see Fig. circ:U_B_a, following [Chen et al. 2023, Fig. 4] [CITE: chen2023efficient]).

Running $U_A$ on $\ket{0^b}\otimes\ket{\psi}$ and *post-selecting* the ancilla on $\ket{0^b}$ leaves, up to normalisation, the state $A\ket{\psi}/\alpha$ on the system, with success probability $\|A\ket{\psi}\|^2/\alpha^2$. When $\alpha \leq 1$, the block encoding can double as a deterministic implementation of $A$. In the typical case $\alpha > 1$, however, post-selection succeeds only with probability $\mathcal{O}(1/\alpha^2)$ and LCU is *not* deterministic. Two canonical textbook remedies exist. *(i) Repeat-until-success.* Running the circuit $\mathcal{O}(\alpha^2)$ times in expectation yields one useful sample; acceptable when $A$ feeds a sampling subroutine whose outputs are anyway averaged, but quadratic in $\alpha$ and incoherent. *(ii) Oblivious amplitude amplification.* When $A/\alpha$ acts as a near-isometry on the relevant input subspace — a property that holds for the truncated Taylor series approximation of $\ee^{-\ii Ht}$ [CITE: berry2015simulating] — a Grover-style reflection on the block ancilla boosts the success probability to $\Theta(1)$ using only $\mathcal{O}(\alpha)$ invocations of $U_A$ and $U_A^\dagger$. The quadratic overhead of (i) is replaced by a linear overhead, and the boosted circuit is itself a unitary acting coherently on system$\otimes$ancilla.

Chapter 5 takes a different route. The coherent correction $B$ is Hermitian and is built as a *nested* LCU (see the block-encoding relation \eqref{eq:B-block-encoding}) with two time registers $T_-, T_+$ and subnormalisations $\alpha_\pm = \|b_\pm\|_1$, so that the composite subnormalisation is
$$
\alpha \;=\; \|b_-\|_1\,\|b_+\|_1 \;=\; \mathcal{O}\!\left(\log\frac{\beta\|H\|\,\|\textstyle\sum_a A^{a\dagger}A^a\|}{\varepsilon}\right)
$$
in the regularised Metropolis case [CITE: chen2023efficient], i.e. polylogarithmic in the problem parameters. For a Hermitian $B$ with mixed-sign spectrum in $[-\alpha,\alpha]$, $B/\alpha$ is **not** close to an isometry — eigenvectors with small eigenvalue are sent to vectors of small norm — so the prerequisite for oblivious amplitude amplification fails on $U_B$ directly. Instead, the $\alpha$ factor is absorbed *at the next circuit layer*: feeding $U_B$ together with the block reflection $R_b := 2\ket{0^b}\bra{0^b} - \idm$ into the qubitization walk operator $W := R_b\,U_B$ [CITE: lowChuang2019qubitization], the *generalized QSP* (GQSP) construction of [Motlagh and Wiebe 2024] [CITE: motlagh2024generalized] realises $\ee^{-\ii\delta B}$ via a Laurent polynomial $P(\ee^{\ii\theta}) \approx \ee^{-\ii\delta\alpha\cos\theta}$ of degree $d = \mathcal{O}(\delta\alpha + \log(1/\varepsilon_\mathrm{QSP})/\log\log(1/\varepsilon_\mathrm{QSP}))$ [Motlagh–Wiebe 2024, Thm. 7] [CITE: motlagh2024generalized]. Because the target is unit-modulus on the qubitization subspace, the composite circuit is already deterministic up to GQSP truncation error [CITE: motlagh2024generalized]. The polylog $\alpha$ above therefore enters *linearly* in the GQSP query count without requiring a separate amplification pass.

The output of the LCU circuit of (Fig. circ:lcu), specialised to the nested two-register kernel construction of Chapter 5, is therefore a unitary $U_B$ acting on system $\otimes\, T_- \otimes T_+$ with
$$
\bra{0^b}\,U_B\,\ket{0^b} \;=\; \frac{B}{\alpha}, \qquad B \text{ Hermitian},
$$
where $\ket{0^b}$ denotes the fiducial all-zero state on the combined block register $T_-\otimes T_+$ and $\alpha = \|b_-\|_1\,\|b_+\|_1$. It is $U_B$, together with the reflection $R_b$ on the block register, that enters the qubitization walk $W := R_b\,U_B$ of §sec:prelim-QSP. The qubitization picture exploits Hermiticity of $B$ to parametrise the spectrum of $W$ in the eigenvalues of $B/\alpha$ [CITE: lowChuang2019qubitization], and Motlagh–Wiebe GQSP [CITE: motlagh2024generalized] then delivers the $\delta$-step of coherent evolution $\ee^{-\ii\delta B}$.

---

## Citations

- **`childs2012hamiltonian`** — Childs, A. M. and Wiebe, N., *"Hamiltonian simulation using linear combinations of unitary operations"*, **Quantum Inf. Comput. 12, 901–924 (2012)**, arXiv:1202.5822. *Already in `references.bib`.*
- **`berry2015simulating`** — Berry, D. W., Childs, A. M., Cleve, R., Kothari, R. and Somma, R. D., *"Simulating Hamiltonian dynamics with a truncated Taylor series"*, **Phys. Rev. Lett. 114, 090502 (2015)**, arXiv:1412.4687. Oblivious amplitude amplification for LCU. *Already in `references.bib`.*
- **`chen2023efficient`** — Chen, C.-F., Kastoryano, M. J. and Gilyén, A., *"An efficient and exact noncommutative quantum Gibbs sampler"*, arXiv:2311.09207 (2023). Nested LCU block encoding of $B$, Proposition III.1; Fig. 4 gives the asymmetric $\textsc{Prep}' / \textsc{Prep}$ circuit. *Already in `references.bib`.*
- **`motlagh2024generalized`** — Motlagh, D. and Wiebe, N., *"Generalized Quantum Signal Processing"*, arXiv:2308.01501 (2024). Theorem 7 (Laurent-polynomial degree bound $d = \mathcal{O}(\delta\alpha + \log(1/\varepsilon)/\log\log(1/\varepsilon))$ for $\ee^{-\ii\delta\alpha\cos\theta}$, equivalently $\ee^{-\ii\delta\alpha\sin\theta}$; Thm. 7 covers both) and Corollary 8 (deterministic success for unit-modulus targets; walk $W = R_b\,U_B$ with $\cos\theta = \lambda/\alpha$). *Already in `references.bib`.*
- **`lowChuang2019qubitization`** — Low, G. H. and Chuang, I. L., *"Hamiltonian Simulation by Qubitization"*, **Quantum 3, 163 (2019)**, arXiv:1610.06546. Modern block-encoding and subnormalisation framing; qubitization of a Hermitian $B/\alpha$ via $W = R_b\,U_B$.

```bibtex
@article{lowChuang2019qubitization,
  author        = {Low, Guang Hao and Chuang, Isaac L.},
  title         = {Hamiltonian Simulation by Qubitization},
  journal       = {Quantum},
  volume        = {3},
  pages         = {163},
  year          = {2019},
  doi           = {10.22331/q-2019-07-12-163},
  eprint        = {1610.06546},
  archivePrefix = {arXiv},
  primaryClass  = {quant-ph},
}
```

---

## Writing Notes

- **Length.** Six body paragraphs + three labelled / unlabelled displays. Expected to sit just under the 2-page cap at thesis line-spacing; the figure `circ:lcu` fills roughly the next half page. The six-paragraph arc is: (1) motivation + scope, (2) primitive + block-encoding identity, (3) signed-$f$ conventions (both variants), (4) success probability + two canonical $\alpha > 1$ remedies, (5) Chapter-5 route (Motlagh–Wiebe GQSP on the qubitization walk absorbs $\alpha$; OAA prerequisite fails on Hermitian $B$), (6) hand-off to §sec:prelim-QSP naming $U_B$, $R_b$, $W$. This structure is what the user's original three-point spec asked for, corrected so that Chapter 5 is described honestly (Motlagh–Wiebe GQSP, not OAA) while the two canonical textbook remedies are still taught.
- **Notation.** `\ket{0^b}` / `\bra{0^b}` throughout, matching (Fig. circ:lcu) and (Fig. circ:gqsp) captions and the GQSP draft. `\textsc{Prep}` / `\textsc{Select}` match the thesis figure caption. `\idm` for the identity in body text (thesis preamble `main.tex:83`), `\ee`, `\ii` for $e$ and $i$. Subnormalisation (no hyphen) per `2_methods.tex:1705`. The asymmetric `PREP' / PREP` pair matches `2_methods.tex:1618–1695` and the figure caption at `circ:U_B_a`.
- **Two signed-$f$ conventions.** Review C2 asked for one of {describe both; describe the symmetric one and forward-point}. I chose to describe both in a single paragraph and name the asymmetric one as "the Chapter-5 choice". This keeps the preliminaries subsection self-contained (a reader never flipping to Chapter 5 still sees the canonical textbook form first), while being honest about the specific variant Chapter 5 consumes. A parenthetical in the paragraph points at `circ:U_B_a` directly.
- **Scaling form.** Took review CS1 / [CHECK] 1 verdict: use the full thesis expression $\alpha = \mathcal{O}(\log(\beta\|H\|\,\|\sum_a A^{a\dagger}A^a\|/\varepsilon))$ copied verbatim from `2_methods.tex:1726`, including the jump-sum factor. No $\tilde O$ shorthand.
- **Hermiticity declaration.** Review cross-draft recommendation (GQSP review §"Things the GQSP draft uses from LCU that LCU actually sets up"): $B$ was previously Hermitian implicitly. v2+ declares this explicitly in two places: (i) the signed-$A$ paragraph now names $A$ as "Hermitian but not unitary" when $f$ is real-symmetric and $U(-t) = U(t)^\dagger$, and (ii) the closing block-encoding identity is tagged "$B$ Hermitian" inline. The GQSP subsection's opening "For Hermitian $B$..." is now properly set up.
- **GQSP hand-off.** The review's C3 fix is applied: replaced "singular values of $U_B$" with the qubitization-eigenvalue language ("qubitization picture exploits Hermiticity of $B$ to parametrise the spectrum of $W$ in the eigenvalues of $B/\alpha$"). The closing sentence names $W$ explicitly as the walk operator and hands it to §sec:prelim-QSP, matching drafter S3 of the original Suggestions list and closing GQSP review item §Cross-reference verdict → Forward refs from LCU to GQSP. Paragraph 5 now names Motlagh–Wiebe GQSP explicitly and cites Thm. 7 on the first introduction of the degree bound, not only as an afterthought.
- **OAA attribution.** Review C1 fix: paragraph 4 now presents OAA as the *generic* remedy (valid for truncated-Taylor Hamiltonian simulation, the use case in Berry et al. 2015), without the false "applies to our Hermitian $B$-block" claim. Paragraph 5 explicitly states why OAA fails on $U_B$ (Hermitian mixed-sign spectrum is not near-isometric) and credits Motlagh–Wiebe GQSP with the $\alpha$-absorption.
- **$\cos\theta$ target (new in v4).** The target polynomial named in paragraph 5 is $P(\ee^{\ii\theta}) \approx \ee^{-\ii\delta\alpha\cos\theta}$, matching QSP v4 (`drafts/qsp-subsection-v4.md`) and MW's own convention for the walk $W = R_b\,U_B$ (their Cor. 8, Eq. 64: $W_\lambda$'s trace is $2\lambda/\alpha$, eigenvalues $\ee^{\pm\ii\arccos(\lambda/\alpha)}$). MW Thm. 7 gives the same additive degree bound for both $\ee^{\pm\ii t\cos\theta}$ and $\ee^{\pm\ii t\sin\theta}$, so the scaling $d = \mathcal{O}(\delta\alpha + \log(1/\varepsilon)/\log\log(1/\varepsilon))$ stands unchanged. The methods chapter currently uses $\sin\theta$ at `2_methods.tex:1731, 1733, 1609` and needs a separate `sin` → `cos` pass to match (user memory `todo_methods_gqsp_reflection.md` tracks this).
- **No Gilyén–Su–Low–Wiebe 2019.** Same reasoning as v1/v2: that key is properly invoked in (§sec:prelim-QSP) for QSVT context; introducing it here would duplicate the framing citation `lowChuang2019qubitization` already covers.
- **Citation key for OAA primary reference.** Took review [CHECK] 2 verdict: keep `berry2015simulating` as the single citation for OAA. After the C1 rewrite OAA is no longer load-bearing for $B$, so a second key (`berry2014exponential`) would clutter the bibliography without pedagogical benefit.

---

## Suggestions

*Items not requested by the user but worth considering for polish. Items 2 and 5 of the v1 list were dropped per the review's assessment; items 1, 3, 4 are retained (with the C2 fix folded into item 1).*

1. **Numbered display for signed-$f$.** If the user wants Chapter 5's `2_methods.tex:1713–1718` derivation to cite back instead of repeating the convention, promote the signed-$f$ paragraph (or the asymmetric-pair formula for $\textsc{Prep}'$ / $\textsc{Prep}$) to a labelled display `eq:lcu-signed-prep`. Because v2+ now describes *both* conventions, the labelled equation should be the asymmetric one that Chapter 5 actually uses: $\textsc{Prep}'_\pm \ket{0} = \sum_t \mathrm{sgn}(f(t))\sqrt{|f(t)|/\|f\|_1}\,\ket{t}$, $\textsc{Prep}_\pm \ket{0} = \sum_t \sqrt{|f(t)|/\|f\|_1}\,\ket{t}$. Low effort, saves Chapter 5 an equation.

2. **Forward link to a labelled $B$-block paragraph in Chapter 5.** Currently paragraph 5 refers to \eqref{eq:B-block-encoding}. If the user adds `\label{sec:B-nested-LCU}` at `2_methods.tex:1707` (the paragraph "The block encoding of the coherent term..."), the subsection can cite it directly as a fuller pointer. One-line user edit.

3. **Connection to the weak-measurement subsection (`sec:prelim-weak-meas`).** The weak-measurement figure at `1_preliminaries.tex:1139–1167` uses $U_L$ to block-encode a single jump $L$ with subnormalisation $\alpha_L$. A half-sentence — "The same primitive also block-encodes the jump operators $L$ of §sec:prelim-weak-meas, with subnormalisation $\alpha_L$ depending on $L$'s construction" — would make (§sec:prelim-LCU) the single source of truth for every block encoding in the thesis. Low priority; post-polish cosmetic.

---

## Review-driven choices applied

The following fixes from `drafts/lcu-subsection-review.md` are applied to this draft (not reproduced as open issues):

- **C1 (OAA vs GQSP retargeting, critical).** Paragraph 4 now presents OAA as the canonical textbook remedy for the *generic* $\alpha > 1$ case (applied via [Berry et al. 2015] to truncated-Taylor Hamiltonian simulation), without the false "applies to our Hermitian $B$" claim. Paragraph 5 explicitly states that $B/\alpha$ is *not* a near-isometry (mixed-sign spectrum sends small-eigenvalue vectors to small norm), so OAA's prerequisite fails on $U_B$; the $\alpha$-absorption is then attributed correctly to Motlagh–Wiebe GQSP on the qubitization walk $W = R_b U_B$ via `motlagh2024generalized` (Thm 7 degree bound + Cor. 8 deterministic-success for unit-modulus targets). The phrase "makes the $B$-block encoding efficient rather than exponentially costly" is removed; the efficient-route story is now told through GQSP degree scaling.
- **C2 (signed-$f$ convention, critical).** Paragraph 3 now describes **both** canonical implementations — the symmetric phase-into-$\textsc{Select}$ variant and the asymmetric $\textsc{Prep}'$ / $\textsc{Prep}$ pair — and names the asymmetric one as the one Chapter 5's `circ:U_B_a` (following [Chen et al. 2023, Fig. 4]) actually uses. The previous false claim "this is exactly the convention used for ... Chapter 5" is removed. This choice (option (a) of the review's two alternatives) keeps the preliminaries self-contained while being honest about Chapter 5's variant.
- **C3 (singular values → qubitization eigenvalues, critical).** Closing paragraph rewritten: "a Laurent-polynomial transformation of $U_B$'s singular values" is replaced by "The qubitization picture exploits Hermiticity of $B$ to parametrise the spectrum of $W$ in the eigenvalues of $B/\alpha$ [CITE: lowChuang2019qubitization], and Motlagh–Wiebe GQSP [CITE: motlagh2024generalized] then delivers the $\delta$-step $\ee^{-\ii\delta B}$." $W$ is named explicitly, matching the GQSP subsection's walk-operator framing.
- **N1 (joint ground state → fiducial all-zero state).** Final display now reads "where $\ket{0^b}$ denotes the fiducial all-zero state on the combined block register $T_-\otimes T_+$".
- **N3 (hyphenation).** All instances of "sub-normalisation" → "subnormalisation" (zero hyphens), matching `2_methods.tex:1685, 1705`.
- **CS1 / [CHECK] 1 ($\alpha$-scaling form).** Display now uses the exact thesis form $\mathcal{O}(\log(\beta\|H\|\,\|\sum_a A^{a\dagger}A^a\|/\varepsilon))$ verbatim from `2_methods.tex:1726`, replacing the $\tilde O(\log(\beta\|H\|/\varepsilon))$ shorthand. No $\tilde O$-hiding.
- **CS4 (spelling).** "generalised quantum signal processing" → "generalized quantum signal processing" throughout, matching the GQSP draft and `2_methods.tex:1705`.
- **E1 (opening sentence "deterministic operation" → "block-encoded unitary circuit").** Fixed in paragraph 1.
- **E2 (drop `eq:lcu-success-prob`).** Removed: paragraph 4's inline form $p_\mathrm{succ} = \|A\ket\psi\|^2/\alpha^2$ is no longer a labelled display. Only `eq:lcu-block-encoding` remains numbered.
- **E3 (price-is-paid polish).** "is the price paid for" → "is the price for" in paragraph 2.
- **E5 (composite $\alpha$ display label).** Deliberately *not* labelled in v2/v3/v4; the display sits inside paragraph 5's argument flow and does not need to be cross-referenced by later chapters (Chapter 5 already has the same formula at `2_methods.tex:1726`, and the Numerics chapter can cite that). The review flagged this as optional.
- **Cross-draft from GQSP review ($B$ Hermitian explicit).** Declared explicitly twice: (i) in paragraph 2's description of $A$ when $f$ is real-symmetric and $U(-t) = U(t)^\dagger$, giving the generic Hermitian-but-non-unitary setting; (ii) inline in the closing block-encoding identity "$B$ Hermitian". The GQSP subsection's opening "For Hermitian $B$..." is now properly prepared.

**Items not applied, with reason:**

- **N2 ($\idm$ vs $I$ inconsistency across drafts).** The LCU body uses `\idm` (two places — both body-text macros, per thesis preamble convention). The thesis figure captions at `1_preliminaries.tex:1109` use plain `I` inside `\ket{0^b}\bra{0^b} - I`. I followed the review's tolerance and kept `\idm` in body text, accepting that the figure caption uses `I`. No change here.
- **MC1 / [CHECK] 2 (BCCKS 2014 OAA primary reference).** Kept `berry2015simulating` as the single OAA citation, per the review's own verdict. After C1, OAA is no longer the load-bearing tool and adding a second key would clutter the bibliography.
- **MC4 (`eq:B-block-encoding-alpha` label).** The review flagged that no such label exists. I reference `\eqref{eq:B-block-encoding}` (which does exist at `2_methods.tex:1702`) and do not introduce the nonexistent label. If the user adds `\label{eq:alpha-metropolis-scaling}` at `2_methods.tex:1726` later, I can cite it in a future revision.
- **CS3 ($U_A$ vs $U_B$ naming across LCU and GQSP).** Review flagged this as "actually good" — no change.
- **CE1 (weak-measurement connection).** Review rated low priority / post-polish; captured in Suggestion 3 above.
- **CE3 (single sentence on why subnormalisation exists).** Review rated "very optional polish"; the current paragraph-2 wording ("the price for packing a non-unitary $A$ inside a unitary") already carries the intuition.

---

## v3 → v4 change summary

**Three surgical `sin` → `cos` fixes to match QSP v4's Motlagh–Wiebe convention $\cos\theta = \lambda/\alpha$.** No other content changes.

1. **Paragraph 5, GQSP forward pointer.** "$P(\ee^{\ii\theta}) \approx \ee^{-\ii\delta\alpha\sin\theta}$" → "$P(\ee^{\ii\theta}) \approx \ee^{-\ii\delta\alpha\cos\theta}$". The degree bound $d = \mathcal{O}(\delta\alpha + \log(1/\varepsilon_\mathrm{QSP})/\log\log(1/\varepsilon_\mathrm{QSP}))$ is unchanged — MW Thm. 7 on p. 11 states that bound for both $\ee^{\pm\ii t\sin H}$ and $\ee^{\pm\ii t\cos H}$ simultaneously.
2. **Citations entry `motlagh2024generalized`.** "Laurent-polynomial degree bound ... for $\ee^{-\ii\delta\alpha\sin\theta}$" → "... for $\ee^{-\ii\delta\alpha\cos\theta}$, equivalently $\ee^{-\ii\delta\alpha\sin\theta}$; Thm. 7 covers both". Corollary 8 gloss expanded to name the walk convention: "walk $W = R_b\,U_B$ with $\cos\theta = \lambda/\alpha$".
3. **New Writing Note.** Added a "$\cos\theta$ target (new in v4)" bullet explaining the convention match with QSP v4 and flagging that `2_methods.tex:1731, 1733, 1609` still use $\sin\theta$ and need a separate pass (tracked in user memory `todo_methods_gqsp_reflection.md`).

No v3 content is otherwise reworded. The three-point user spec, the six-paragraph arc, the OAA-vs-GQSP reframing, the asymmetric $\textsc{Prep}'/\textsc{Prep}$ discussion, the $\alpha$-scaling form, the Hermiticity declaration, and the GQSP hand-off are all preserved verbatim from v3.

---

## v2 → v3 change summary (retained for history)

**Single surgical change in v3: the forward pointer to §sec:prelim-QSP was retargeted to Motlagh–Wiebe GQSP applied to the qubitization walk $W = R_b U_B$, reflecting that the QSP subsection actually uses Motlagh–Wiebe GQSP (not Low–Chuang QSP).**

The v2 draft already named GQSP in paragraphs 5 and 6 and already cited `motlagh2024generalized` at the degree-bound display. v3 tightened the attribution in three places:

- **Paragraph 1, final sentence.** "...because that is the interface on which the *QSP* subsection builds" → "...because that is the interface on which the *GQSP* subsection builds". The label `sec:prelim-QSP` is unchanged (it is the user's real label in `1_preliminaries.tex:1068`).
- **Paragraph 2, final sentence.** "*QSP* polynomial degree, simulation time, success probability" → "*GQSP* polynomial degree, simulation time, success probability". Consistent with what Chapter 5 actually uses.
- **Paragraph 5, the "absorbed at the next circuit layer" sentence.** Rewritten so that the walk operator $W = R_b U_B$ appears **inside** the sentence that names Motlagh–Wiebe GQSP (not only in paragraph 6), and [Motlagh and Wiebe 2024] is cited by name (not only as `motlagh2024generalized` in the degree formula). The Laurent-polynomial target $P(\ee^{\ii\theta}) \approx \ee^{-\ii\delta\alpha\sin\theta}$ from `2_methods.tex:1731` was explicitly stated here in v3 (v2 only had it implicit via the degree formula); v4 now replaces $\sin$ with $\cos$ in this sentence (see v3 → v4 summary above).
- **Paragraph 6, final clause.** "GQSP [CITE: motlagh2024generalized] then delivers…" → "Motlagh–Wiebe GQSP [CITE: motlagh2024generalized] then delivers…". Named attribution on the final forward hand-off.

The three-point user spec is untouched: (i) block encoding and $\|0^b\rangle$–$\|0^b\rangle$ projector identity in paragraphs 2–3, (ii) success probability $\|A\|\psi\rangle\|^2/\alpha^2$ plus the $\alpha \leq 1$ special case in paragraph 4, (iii) the $\alpha > 1$ remedies — repeat-until-success and oblivious amplitude amplification as canonical textbook answers in paragraph 4, followed by Motlagh–Wiebe GQSP on the qubitization walk as the thesis's actual choice in paragraph 5 — is preserved verbatim from v2. The OAA-vs-GQSP reframing (GQSP absorbs $\alpha$ into polynomial degree; OAA fails on Hermitian $B/\alpha$ because it is not a near-isometry) is preserved. The asymmetric $\textsc{Prep}'$/$\textsc{Prep}$ convention discussion in paragraph 3 is preserved.

No new citation keys introduced in v3 beyond those already in v2. `motlagh2024generalized` (already in `references.bib` per `2_methods.tex:1705, 1733`) is cited in paragraph 5 as well as paragraph 6; `lowChuang2019qubitization` — the Low–Chuang 2019 qubitization reference for the walk $W = R_b U_B$ — is kept as v2 already introduced it, because Low–Chuang 2019 is the correct attribution for the qubitization construction even though the signal-processing layer is Motlagh–Wiebe GQSP. The Berntson–Sünderhauf angle-finding reference belongs to the GQSP subsection, not here, and is not cited in v3/v4.
