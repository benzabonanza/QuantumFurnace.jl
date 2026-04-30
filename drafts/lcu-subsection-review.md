# Draft Review: LCU Subsection

**Reviewed**: `drafts/lcu-subsection.md`
**Target insertion**: `supplementary-informations/1_preliminaries.tex:1037`, body of `\subsection{Linear Combination of Unitaries.}\label{sec:prelim-LCU}`.
**Date**: 2026-04-23
**Overall assessment**: Ship with fixes. Exposition is clean and the three stated goals are hit, but **the closing two paragraphs contain a substantive logical error** about why the $\alpha$-scaling for $B$ is not fatal: the thesis does *not* use oblivious amplitude amplification on $U_B$ — GQSP already absorbs $\alpha$ into the polynomial degree, and `drafts/coherent-step.md:120, 144` explicitly states "No amplitude amplification is needed." The signed-$f$ paragraph also describes a different convention than the one Chapter 5 actually uses, and the signal-processing description elides "singular values" where the thesis specifically works with Hermitian $B$ (eigenvalues) via qubitization. One minor notation issue ("joint ground state" for the all-zero computational-basis ancilla).

---

## Summary Scorecard

| Category | Rating | Notes |
|----------|--------|-------|
| Logical correctness | **MAJOR** | Paragraph 5 misattributes the thesis's $\alpha$-absorption to OAA; it is actually GQSP. |
| Notation consistency | **MINOR** | "Joint ground state" is off for $\ket{0}^{\otimes b}$; "singular values of $U_B$" is wrong (should be eigenvalues of $B/\alpha$, via qubitization on Hermitian $B$). |
| Citations | **MINOR** | `berry2015simulating` is acceptable for OAA but BCCKS 2014 (arXiv:1312.1414) is the primary reference. `lowChuang2019qubitization` stub is structurally valid. |
| Thesis consistency | **MAJOR** | The signed-$f$ paragraph states "the circuit structure is unchanged: absorb the phase into `\textsc{Select}`"; Chapter 5 actually uses an asymmetric PREP'/PREP convention (distinct left/right preparations) — the draft's claim "This is exactly the convention used for the signed kernels $b_-$ and $b_+$ of Chapter 5" is false in letter (same subnormalisation $\|f\|_1$, but different circuit). |
| Exposition quality | OK | Six-paragraph structure is tight; primitive / signed-$f$ / success / remedies / Chapter-5 hand-off / $U_B$ closing paragraph all land cleanly modulo the corrections below. |

---

## Critical Issues

### C1: Paragraph 5 misattributes $\alpha$-absorption to OAA; the thesis actually uses GQSP

- **Location**: draft lines 35 (sentence starting "*(ii) Oblivious amplitude amplification.*") and 41 ("the $\Theta(1)$-success oblivious amplitude amplification of [CITE: berry2015simulating] keeps that factor coherent and linear rather than quadratic, and is what makes the $B$-block encoding efficient rather than exponentially costly once ported onto the GQSP walk").
- **Problem**: The thesis does *not* apply OAA to $U_B$. The QSP/GQSP subsection (see `drafts/qsp-subsection.md` lines 36–41, and the already-committed `drafts/coherent-step.md:120, 144`) takes the view that the target function $\ee^{-\ii\delta\alpha\sin\theta}$ has unit modulus, so the QSP-side "success probability is $1 - \mathcal{O}(\varepsilon_\mathrm{QSP})$" *before* any amplification, by Motlagh–Wiebe Corollary 8. The $\alpha$ factor enters the GQSP polynomial degree as $d = \mathcal{O}(\delta\alpha + \log(1/\varepsilon)/\log\log(1/\varepsilon))$, i.e. linearly in $\alpha$, through the Jacobi–Anger construction — *not* through OAA on $U_B$. The draft's current phrasing suggests that OAA is the reason the $B$-block encoding is efficient, which contradicts the thesis's own Coherent-Step draft and its eventual algorithm.
- **Why it matters**: This is the load-bearing pedagogical hand-off from the LCU subsection to the QSP subsection. Claiming OAA does the work sets up an inconsistency that the QSP subsection will then have to contradict, and a reviewer would catch it immediately.
- **Secondary concern**: OAA requires the block-encoded operator to be close to a unitary (more precisely, close to an isometry on the relevant subspace). For a Hermitian $B$ with eigenvalues in $[-\alpha, \alpha]$ (the thesis uses $\alpha = \|b_-\|_1\|b_+\|_1$ and $B/\alpha$ has eigenvalues in $[-1,1]$), $B/\alpha$ is *not* close to unitary in general: an eigenvector with eigenvalue near $0$ gets sent to a vector of small norm. So even if one wanted to apply OAA directly to $U_B$, the prerequisite would fail for general $\ket{\psi}$. The draft's parenthetical "which happens whenever $A$ is close to a unitary, as is the case for the Hamiltonian-simulation Taylor truncation and for our Hermitian $B$-block" is wrong on the "Hermitian $B$-block" half: Hamiltonian-simulation Taylor truncation is close to a unitary because it approximates $\ee^{-\ii Ht}$; a Hermitian $B$ with mixed-sign spectrum is not.
- **Suggested fix**: rewrite the last sentence of paragraph 5 (currently "the $\Theta(1)$-success oblivious amplitude amplification of [CITE: berry2015simulating] keeps that factor coherent and linear rather than quadratic, and is what makes the $B$-block encoding efficient rather than exponentially costly once ported onto the GQSP walk of §sec:prelim-QSP") to correctly credit GQSP:

  > In Chapter 5, the $\alpha$ factor is absorbed not via oblivious amplitude amplification of $U_B$ itself, but at the next circuit layer: the GQSP polynomial implementing $\ee^{-\ii\delta B}$ from the qubitized walk operator has degree $d = \mathcal{O}(\delta\alpha + \log(1/\varepsilon_\mathrm{QSP})/\log\log(1/\varepsilon_\mathrm{QSP}))$ (see §sec:prelim-QSP), and the target function $\ee^{-\ii\delta\alpha\sin\theta}$ has unit modulus on the qubitization subspace, so the composite circuit is already deterministic up to QSP truncation error. The polylog $\alpha$ of (Eq. eq:B-block-encoding-alpha) therefore enters linearly in the GQSP query count and does not require a separate amplification pass.

  And loosen paragraph 5's parenthetical so that OAA is presented as the *generic* remedy (true for truncated-Taylor Hamiltonian simulation) without claiming it applies to the $B$-block:

  > *(ii) Oblivious amplitude amplification.* When $A/\alpha$ acts as a near-isometry on the relevant input subspace — a property that holds for the truncated Taylor series approximation of $\ee^{-\ii Ht}$, as used in [CITE: berry2015simulating] — a Grover-style reflection on the block ancilla boosts the success probability to $\Theta(1)$ using only $O(\alpha)$ invocations of $U_A$ and $U_A^\dagger$. The overhead of (i) is replaced by a linear overhead, and the boosted circuit is itself a unitary acting coherently on system $\otimes$ ancilla.

### C2: Signed-$f$ paragraph describes the wrong Chapter-5 convention

- **Location**: draft line 28.
- **Problem**: The draft writes
  > absorb the phase into `\textsc{Select}` by redefining $\tilde U(t) := \mathrm{sgn}(f(t))\,U(t)$ with $\mathrm{sgn}(f(t)) := f(t)/|f(t)|$ … and let `\textsc{Prep}` load the amplitude-magnitude state $\ket{f} = \frac{1}{\sqrt{\|f\|_1}}\sum_t \sqrt{|f(t)|}\,\ket{t}$; … This is exactly the convention used for the signed kernels $b_-$ and $b_+$ of Chapter 5 (see \eqref{eq:b_minus}, \eqref{eq:b_plus-s-eta}).
- This is **not** the Chapter-5 convention. The nested-LCU circuit at `2_methods.tex:1618–1695` (also `drafts/circuits/block-encoding-ba-halfctrl.tex`) uses an **asymmetric** pair of preparations: `PREP'_±` (signed / phased amplitude $b_\pm(\bar t)/\sqrt{|b_\pm(\bar t)|}$) on the left-hand side, and `PREP_±` (magnitude-only amplitude $\sqrt{|b_\pm(\bar t)|}$) on the right-hand side. The sign ends up in one PREP, not in SELECT. The two conventions give the same final operator (both block-encode $A/\|f\|_1$) but they are genuinely different circuits.
- **Why it matters**: The preliminaries subsection claims to set up the convention Chapter 5 consumes. If the two convention stories disagree, the reader who flips back from `2_methods.tex:1707` to the preliminaries will see conflicting circuits.
- **Suggested fix**: either (a) describe both conventions, or (b) describe the symmetric PREP / PREP$^\dagger$ variant as the textbook form and then point the reader at Chapter 5's specific signed convention. Concretely, replace the sentence in draft line 28 with:

  > If $f$ takes signed or complex values, two equivalent implementations are in common use: one absorbs $\mathrm{sgn}(f(t)) := f(t)/|f(t)|$ into `\textsc{Select}` (redefining $\tilde U(t) := \mathrm{sgn}(f(t))\,U(t)$ and letting `\textsc{Prep}` load the magnitude state $\ket{|f|} = \frac{1}{\sqrt{\|f\|_1}}\sum_t \sqrt{|f(t)|}\,\ket{t}$), the other splits `\textsc{Prep}` into an asymmetric left/right pair `\textsc{Prep}' / \textsc{Prep}` carrying phase and magnitude respectively. The subnormalisation $\alpha = \|f\|_1 = \sum_t |f(t)|$ is the same in both. Chapter 5 adopts the asymmetric variant (Fig. circ:U_B_a, following [CITE: chen2023efficient, Fig. 4]) so that each of the signed kernels $b_-$ (Eq. eq:b_minus) and $b_+$ (Eq. eq:b_plus-s-eta) is carried by its own PREP' / PREP pair on two separate time registers.

### C3: "Singular values" is the wrong word for QSP on Hermitian $B$

- **Location**: draft line 47, final sentence of the subsection: "a Laurent-polynomial transformation of $U_B$'s singular values turns the block encoding of $B/\alpha$ into a $\delta$-step of coherent evolution $\ee^{-\ii\delta B}$."
- **Problem**: The QSP / qubitization picture used in `drafts/qsp-subsection.md` (lines 18–22) is specifically the Hermitian-case picture: the walk operator $W = R_b\,U_B$ has eigenvalues $\ee^{\pm\ii\arccos\lambda}$ where $\lambda$ is an *eigenvalue* of $B/\alpha$. That is a transformation of the *eigenvalues* of $B/\alpha$, not the *singular values* of $U_B$ (those are all 1, since $U_B$ is unitary). The singular-value picture is QSVT (Gilyén–Su–Low–Wiebe 2019) for general non-Hermitian block encodings; the thesis deliberately stays in the simpler qubitization picture because $B$ is Hermitian.
- **Why it matters**: Minor mathematical sloppiness that the next subsection immediately contradicts — the QSP draft at line 18 opens "For Hermitian $B$, a single extra reflection promotes $U_B$ to a walk operator" and is explicit that eigenvalues are the target. The LCU draft's language here pre-commits to the wrong framework.
- **Suggested fix**: replace "a Laurent-polynomial transformation of $U_B$'s singular values" with "a Laurent-polynomial transformation of the walk operator $W = R_b\,U_B$'s spectrum (which for Hermitian $B$ is parametrised by the eigenvalues of $B/\alpha$)", or more briefly: "a polynomial transformation of the eigenvalues of $B/\alpha$ (inherited via qubitization of $W = R_b\,U_B$)".

---

## Notation Issues

### N1: "Joint ground state" for $\ket{0^b}$ is non-standard

- **Location**: draft line 47.
- **Problem**: Calling $\ket{0^b} := \ket{0}^{\otimes b}$ the "joint ground state" of the $T_-, T_+$ ancillas suggests $T_\pm$ have a Hamiltonian of which $\ket{0}$ is the lowest-energy eigenstate — they are just fiducial ancilla registers. The LCU figure at `1_preliminaries.tex:1042–1047` and the QSP figure at `1_preliminaries.tex:1090–1095` just say "ancilla in $\ket{0^b}$".
- **Fix**: replace "joint ground state" with "fiducial all-zero state" (or simply "initial state of the combined ancilla registers"). E.g.:

  > where $\ket{0^b}$ denotes the fiducial all-zero state on the combined block register $T_- \otimes T_+$, and $\alpha = \|b_-\|_1\,\|b_+\|_1$.

### N2: `\idm` (thesis-defined) used correctly; `\ee`, `\ii` used correctly

Confirmed via `main.tex:83,88–90`. No action needed. Minor observation: the draft uses `I` (not `\idm`) in `R_b = 2\ket{0^b}\bra{0^b} - I` in `drafts/qsp-subsection.md:20` (for consistency with the thesis, which also writes `R_b = 2\ket{0^b}\!\bra{0^b} - I` in `1_preliminaries.tex:1109`), but in the LCU draft `\idm` is used in two places (draft lines 22, 47). This is internally inconsistent within the two drafts — not actually a problem because both appear in figure captions / body text respectively, but it is worth flagging so the user picks one convention. Recommended: use `\idm` in body text (per thesis preamble convention) and tolerate `I` inside `\gate{}` / `\bra{}\ket{}` constructions where `\idm` would render less cleanly.

### N3: `sub-normalisation` hyphenation

- **Location**: throughout draft (lines 26, 28, 37).
- The thesis body uses `subnormalisation` without a hyphen (e.g. `2_methods.tex:1685, 1705`). The draft writes `sub-normalisation` with a hyphen. Very minor; pick one and be consistent. Recommended: drop the hyphen to match the thesis.

### N4: `\textsc{Prep}` / `\textsc{Select}` convention

- **Location**: draft lines 18–22.
- **Thesis**: The LCU figure caption at `1_preliminaries.tex:1057–1059` uses `\textsc{Prep}` / `\textsc{Select}` in body text and `\gate{\text{PREP}}` inside the quantikz gate label. The draft mirrors this exactly. ✓ No action needed.

### N5: `\ket{0^b}`-`\bra{0^b}` (draft) vs `\ket{0}_b` / `\bra{0}_b` alternatives

Verified: the draft uses `\ket{0^b}` / `\bra{0^b}` exclusively (lines 22, 24, 25 (in the HTML label), 30, 32, 45, 47). This exactly matches the thesis's established convention from `1_preliminaries.tex:1042, 1047, 1090, 1095, 1139, 1145`, `drafts/qsp-subsection.md:14`, and `drafts/weak-measurement-subsection.md`. ✓ The drafter's "Notation harmony" note is correct.

---

## Missing Citations

### MC1: OAA primary reference is BCCKS 2014, not Berry et al. 2015

- **Location**: draft lines 35, 41 (both citing `berry2015simulating` for OAA).
- **Observation**: Oblivious amplitude amplification was *introduced* in Berry, Childs, Cleve, Kothari, Somma 2014 "Exponential improvement in precision for simulating sparse Hamiltonians" (arXiv:1312.1414), where the abstract explicitly introduces "a new form of oblivious amplitude amplification". The 2015 PRL (arXiv:1412.4687) uses it as a tool but is not the introducing paper. The thesis already cites `berry2015simulating` elsewhere (`1_preliminaries.tex:1062`, `2_methods.tex:1709`) alongside `childs2012hamiltonian`, so using it is not wrong, but the cleaner primary reference is BCCKS 2014.
- **Recommendation**: either keep `berry2015simulating` (acceptable — the thesis already uses it as the LCU-era citation and the typical community practice conflates the two Berry papers), or add a second key `berry2014exponential` (stub below) for the specific OAA citation. The drafter's [CHECK] note at the bottom of the draft flagged this correctly; my recommendation is **keep `berry2015simulating`** for consistency with the existing thesis citations *but* drop the "observed" language in paragraph 5. After C1's rewrite the attribution becomes lower-stakes since OAA is no longer doing the hard work.

  Optional bibtex stub if the user wants to add the 2014 paper:
  ```bibtex
  @article{berry2014exponential,
    author  = {Berry, Dominic W. and Childs, Andrew M. and Cleve, Richard and Kothari, Robin and Somma, Rolando D.},
    title   = {Exponential Improvement in Precision for Simulating Sparse Hamiltonians},
    journal = {Forum of Mathematics, Sigma},
    volume  = {5},
    pages   = {e8},
    year    = {2017},
    doi     = {10.1017/fms.2017.2},
    eprint  = {1312.1414},
    archivePrefix = {arXiv},
    primaryClass  = {quant-ph},
  }
  ```
  (conference version in STOC 2014; journal version in *Forum of Mathematics, Sigma* 2017).

### MC2: `lowChuang2019qubitization` bibtex stub is structurally valid

Verified the stub at draft lines 58–70 has all mandatory fields (`author`, `title`, `journal`, `volume`, `pages`, `year`, `doi`, `eprint`, `archivePrefix`, `primaryClass`). DOI `10.22331/q-2019-07-12-163` is correct for the Quantum journal version. Volume `3` and pages `163` are correct. ✓ Note: the QSP draft introduces the same key with the *same* stub at `drafts/qsp-subsection.md:73–84` — the user only needs to commit it once. This is *not* a duplicate bug in the LCU draft; it is consistent cross-draft.

### MC3: Missing cross-reference to Figure `circ:lcu` in the body text

- **Location**: paragraph 2 (draft line 16) introduces `(Fig. circ:lcu)` once. Paragraph 5 and 6 do not re-reference the figure.
- Observation: this is a stylistic choice. Given that the figure is *right below* the body text, one reference is sufficient. ✓ No action needed.

### MC4: Equation label `eq:B-block-encoding-alpha` does not exist

- **Location**: my C1 fix suggestion above refers to "(Eq. eq:B-block-encoding-alpha)". 
- That label is *not* defined anywhere in the thesis. `2_methods.tex:1726` carries an unlabelled display. If the user wants to cite the $\alpha$-scaling formula later, they would need to add `\label{eq:alpha-metropolis-scaling}` (or similar) at that display. This is a housekeeping note on my C1 patch, not a defect of the draft.

---

## Consistency Issues

### CS1: `\tilde O(\log(\beta\|H\|/\varepsilon))` omits the $\|\sum_a A^{a\dagger}A^a\|$ factor of the thesis expression

- **Location**: draft line 39.
- **Thesis expression** (`2_methods.tex:1726`): $\alpha = \mathcal{O}\!\left(\log\frac{\beta\|H\|\,\|\sum_a A^{a\dagger}A^a\|}{\varepsilon}\right)$ — includes the squared-jump sum inside the log.
- **Draft expression** (line 39): $\alpha = \tilde{O}\!\left(\log\frac{\beta\|H\|}{\varepsilon}\right)$ — drops the jump-sum factor and uses `\tilde O` (polylog-hiding) rather than plain `\mathcal{O}`.
- **Assessment**: both forms are defensible — the jump-sum is $O(1)$ per jump for Pauli-string jumps and the `\tilde O` absorbs it. The drafter's [CHECK] note correctly flagged this. For a preliminaries subsection that prepares the reader for Chapter 5, I recommend matching Chapter 5's exact form to avoid later surprise when the reader sees the full expression on the first pass:

  > Suggested replacement for draft lines 37–41:
  > $$\alpha = \|b_-\|_1\,\|b_+\|_1 = \mathcal{O}\!\left(\log\frac{\beta\|H\|\,\|\sum_a A^{a\dagger}A^a\|}{\varepsilon}\right)$$
  > in the regularised Metropolis case [CITE: chen2023efficient] (see \eqref{eq:alpha-metropolis-scaling} in Chapter 5), i.e. polylogarithmic in the problem parameters.

  This matches `2_methods.tex:1726` to the character and avoids introducing a $\tilde O$ shorthand that the reader would then need to unpack on reaching Chapter 5. (It also resolves the drafter's [CHECK] 6.)

### CS2: Draft paragraph 5 says OAA works for "Hermitian $B$-block"; conflicts with `drafts/coherent-step.md:120, 144`

Already covered in C1. Flagging again here because it is a cross-draft contradiction: the coherent-step draft's [RESOLVED] note ("No amplitude amplification is needed") is authoritative for Chapter 5.

### CS3: Minor inconsistency between $U_A$ / $U_B$ and $U$ across circuits

- **Location**: draft lines 22, 30, 43; the thesis's `circ:lcu` caption at `1_preliminaries.tex:1057` uses the generic "$A = \sum_t f(t)\,U(t)$" — no `$U_A$` symbol. The draft introduces `$U_A$` in body text (line 22, 30) and then `$U_B$` in the closing paragraph (line 43). ✓ This is actually *good* — it gives the composite circuit a name before specialising to $B$. The only concern: the QSP subsection figure uses `$U_B$` as a single two-qubit gate label (`\gate[2]{U_B}` at `1_preliminaries.tex:1092`) without further decomposition, so the reader transitions from the LCU's three-piece circuit $U_A = \textsc{Prep}^\dagger \cdot \textsc{Select} \cdot \textsc{Prep}$ to the QSP's black-box $U_B$ cleanly. No action needed.

### CS4: "generalised" vs "generalized" spelling

- **Location**: draft line 47 uses "generalised".
- The thesis uses both spellings; `2_methods.tex:1705` writes "generalized" (American). `drafts/qsp-subsection.md:16` uses "generalized". Recommend switching the draft's "generalised" (line 47) to "generalized" for cross-draft consistency with the QSP draft.

---

## Exposition Improvements

### E1: Opening sentence flows but "deterministic operation" is slightly misleading

- **Location**: draft line 10, second sentence: "a way to realise a weighted sum of unitaries as a single deterministic operation".
- LCU per se is not deterministic — it succeeds only with probability $1/\alpha^2$, which the subsection itself then explains. "Deterministic" here belongs with the OAA-boosted version in paragraph 5, not with the block-encoding primitive. Reword to:

  > … as a single block-encoded unitary circuit, at a controlled cost in ancilla qubits and success probability.

- Low-priority stylistic nit; the reader who finishes the subsection understands the caveat.

### E2: Label `eq:lcu-success-prob` is introduced but not cited

- **Location**: draft line 33 introduces `\label{eq:lcu-success-prob}`.
- Neither the rest of this subsection nor the QSP subsection draft refers back to `eq:lcu-success-prob`. If no Chapter 5 argument needs this as a labelled equation, drop the label (Eq. eq:lcu-block-encoding at line 25 is the one Chapter 5 needs).
- If the user wants to keep it for forward citations from the algorithms chapter (where post-selection cost shows up), keep it.
- My recommendation: drop `eq:lcu-success-prob`; keep `eq:lcu-block-encoding`.

### E3: Paragraph 4 phrase "which is a price paid for packing a non-unitary $A$ inside a unitary"

- **Location**: draft line 26.
- Reads well. Minor polish: "the sub-normalisation is the price for packing a non-unitary $A$ inside a unitary" — drops the existential "is paid". Not critical.

### E4: Length is well within the 2-page cap

- 6 body paragraphs + 4 display equations (2 labelled). At the thesis's formatting this comes to ~1.1 pages of body text, plus the figure. ✓ Well under cap; the signed-$f$ demotion suggestion in the drafter's Writing notes (bullet 3) is not needed.

### E5: The Chapter-5 transition paragraph could carry an equation label for the composite $\alpha$

- **Location**: draft lines 37–39.
- Suggestion: label the $\alpha$-scaling display `\label{eq:alpha-metropolis-scaling}` so that later chapters (including the Numerics chapter) can cite it by equation number rather than by the Chapter-5 paragraph. Small one-liner, high future-proofing value.

---

## [CHECK] Resolutions

### [CHECK] 1: $\alpha$-scaling phrasing — $\tilde O(\beta\log(1/\varepsilon))$ vs $\mathcal{O}(\log(\beta\|H\|/\varepsilon))$

**Verdict**: Use the thesis's full $\mathcal{O}(\log(\beta\|H\|\,\|\sum_a A^{a\dagger}A^a\|/\varepsilon))$ form exactly as written at `2_methods.tex:1726`. Rationale:

1. The drafter's intuition is right that the two forms are inequivalent shorthands for different objects: $\tilde O(\beta\log(1/\varepsilon))$ would hide the single-log structure and the full $\|H\|$ dependence; $\mathcal{O}(\log(\beta\|H\|/\varepsilon))$ is the correct Chen 2023 expression but drops the squared-jump factor.
2. For the preliminaries subsection, the priority is consistency with the chapter the reader will hit next. The full form is what `2_methods.tex:1726` commits to; copying it verbatim into the preliminaries avoids the reader's "wait, which form is it?" moment on first contact with Chapter 5.
3. User's spec ($\tilde O(\beta\log(1/\varepsilon))$) appears to be from an earlier bibliographic note — the thesis body has moved past it.

Patch: see CS1 above.

### [CHECK] 2: OAA primary reference — `berry2015simulating` vs BCCKS 2014

**Verdict**: Keep `berry2015simulating` as the primary citation for OAA. Rationale:

1. The thesis already uses `berry2015simulating` at `1_preliminaries.tex:1062` and `2_methods.tex:1709` — introducing a second OAA citation for a one-paragraph context is weaker than staying consistent.
2. After C1's rewrite, OAA is no longer the load-bearing remedy for the $B$-block, so the strength of the citation matters less.
3. The community practice is to cite Berry et al. 2015 or Berry et al. 2014 interchangeably for OAA; both are acceptable.

Optional: if the user decides the OAA-introducing paper warrants its own key, the stub in MC1 above is ready to drop in.

---

## Suggestions Assessment (the drafter's own "Suggestions" section)

The drafter listed six suggestions at draft lines 86–95. My independent opinions:

1. **Numbered display for signed-$f$** — *Worth doing.* Chapter 5 at `2_methods.tex:1713–1718` re-derives the signed convention from scratch. A labelled preliminaries equation `\eqref{eq:lcu-signed-prep}` that Chapter 5 can cite would remove that duplication. However, make sure the displayed identity matches the *asymmetric* PREP'/PREP convention Chapter 5 actually uses (see C2), not the symmetric phase-into-SELECT variant the current signed-$f$ paragraph states. Low effort, but hinges on fixing C2 first.

2. **Gilyén–Su–Low–Wiebe 2019 citation** — *Optional.* I agree with the drafter's judgment that `gilyen2019quantum` belongs properly to §sec:prelim-QSP (where QSVT would be introduced if the thesis used it; see the QSP draft line 168 for the same thinking). Adding it alongside `lowChuang2019qubitization` in the LCU subsection would introduce a third reference for a framing point that already has two. Skip.

3. **Preview $W = R_b U_B$ as the walk operator** — *Worth doing.* The LCU subsection already mentions $R_b$ in its closing paragraph (draft line 47) but does not name $W$. A one-sentence preview — "We will call the composite $W = R_b\,U_B$ the *walk operator*; it is the primitive the QSP subsection consumes." — would streamline the transition to §sec:prelim-QSP without duplicating QSP's derivation. Low effort, clean hand-off.

4. **Forward-link to Chapter 5's labelled $B$-block paragraph** — *Worth doing if* the user is willing to add a `\label{sec:B-nested-LCU}` in `2_methods.tex:1707` (a one-line edit). Otherwise skip: \eqref{eq:B-block-encoding} is already the best available handle. The drafter's existing cross-references to \eqref{eq:b_minus}, \eqref{eq:b_plus-s-eta} are sufficient for the preliminaries reader.

5. **"Chapter 5's GQSP walk consumes $U_B$ coherently, so the relevant regime is (ii)"** — *Drop after C1.* This suggestion assumes paragraph 5 correctly presents OAA as Chapter 5's tool. After C1 rewrites paragraph 5 to credit GQSP instead, the relevant statement is "Chapter 5 consumes $U_B$ coherently through GQSP, which absorbs the $\alpha$ factor into the polynomial degree without a separate amplification pass" — and that belongs in the *GQSP* context sentence, not as a coda to paragraph 5. If the rewrite from C1 is adopted, this suggestion becomes redundant.

6. **Scaling clarification** — *Already resolved* by CS1 / [CHECK] 1 above. Use the thesis's exact $\mathcal{O}(\log(\beta\|H\|\,\|\sum_a A^{a\dagger}A^a\|/\varepsilon))$ form; no footnote needed.

---

## Connections and Enrichments

### CE1: Hook to the Weak-Measurement subsection

- The weak-measurement figure at `1_preliminaries.tex:1139–1167` uses $U_L$ to block-encode a single jump $L$ with subnormalisation $\alpha$. This is another use of LCU within the same chapter. A half-sentence at the end of the subsection acknowledging this reuse — "The same primitive also block-encodes the jump operators $L$ of §sec:prelim-weak-meas with subnormalisation $\alpha_L$ depending on $L$'s construction" — would make the LCU subsection the single source of truth for every block encoding in the thesis.
- Low priority; the user can add this during polish.

### CE2: Connection to Chapter 2's statistical-mechanics framing

- Not needed; Chapter 4's preliminaries are deliberately decoupled from the Lindbladian-side chapters 2–3. ✓

### CE3: Single sentence on "why sub-normalisation exists"

- The drafter handles this implicitly in paragraph 4 ("the price paid for packing a non-unitary $A$ inside a unitary"). Could be sharpened: "The block-encoding identity $U_A^\dagger U_A = I$ forces $\|A/\alpha\| \leq 1$, so no choice of block-encoding circuit can reach $\alpha < \|A\|$; the $\ell_1$ bound $\alpha = \|f\|_1 \geq \|A\|$ of the LCU construction is usually not far from optimal when $f$ is a well-chosen kernel."
- Very optional polish.

---

## Headline Verdict

**Ship with fixes.** Counts: **3 Critical** (C1 OAA misattribution, C2 signed-$f$ convention, C3 singular-values-vs-eigenvalues), **5 Notation / Consistency** (N1 "joint ground state", N3 hyphenation, CS1 $\alpha$ form, CS3 inconsistency with GQSP, CS4 spelling), **5 Exposition** (E1–E5), **2 [CHECK] resolutions** (scaling form = thesis form; OAA primary reference = keep 2015). Both [CHECK]s flagged by the drafter were good instincts, and one led to C1's discovery that the paragraph downstream is actually saying the wrong thing about how $\alpha$ gets absorbed. Once C1–C3 are patched, the subsection is a clean, tight fit for `1_preliminaries.tex:1037`.
