# Draft Review v4: Trotterization (Preliminaries subsection)

**Reviewed**: `/Users/bence/code/QuantumFurnace.jl/drafts/trotterization-v4.md`
**Compared against**: `/Users/bence/code/QuantumFurnace.jl/drafts/trotterization.md` (v3) and `/Users/bence/code/QuantumFurnace.jl/drafts/trotterization-review.md` (v3 review)
**Target insertion**: body of `\subsection{Trotterization.}` in `/Users/bence/code/QuantumFurnace.jl/supplementary-informations/1_preliminaries.tex` (between line 1034 and 1079).
**Primary source verified against**: `/Users/bence/code/QuantumFurnace.jl/supplementary-informations/Childs et al. - 2021 - Theory of Trotter Error with Commutator Scaling.pdf`.
**Downstream consumer**: `/Users/bence/code/QuantumFurnace.jl/supplementary-informations/2_methods.tex` (lines 881–1245).
**Date**: 2026-04-27

---

## Overall assessment

V4 is a clean, well-executed trim of v3. The two known v3 issues are fully resolved: B1 (factor-of-2 in the explicit M-count) is eliminated by deletion as planned, and M1 (Principle conflating summand vs bond graph) is properly refactored as a two-layer count with the all-to-all TFIM correctly flagged as the converse-failure case. The label-name harmonisations all verify against `1_preliminaries.tex` at the cited line numbers. The citation table is pruned and every Childs citation in the body is anchored to a verified PDF page/equation.

There is **one [BLOCKER]**: the `eq:H_heis` label that v4 uses (twice) is **not defined anywhere in `1_preliminaries.tex` or `2_methods.tex`** — and was not defined in v3 either. The v4 metadata claims it lives at `1_preliminaries.tex:1048`, but that line **only references** the label (`\eqref{eq:H_heis}`); the defining equation is missing. This is inherited from the surrounding tex, not introduced by v4, but v4 propagates the broken reference without flagging it.

The remaining items are minor or nits.

---

## [BLOCKER] B1. `eq:H_heis` is referenced but never defined

- **Location**: v4 line 6 (metadata: "Existing labels referenced: ... `eq:H_heis` (the Heisenberg-Hamiltonian display already typed into `1_preliminaries.tex`)"), v4 line 86 ("the 1D periodic Heisenberg chain (Eq. eq:H_heis)"), v4 Writing Notes line 129 ("v3's `eq:heisenberg-hamiltonian` → `eq:H_heis`").
- **Claim under review**: that `eq:H_heis` already exists in `1_preliminaries.tex` and v4 is harmonising to a label that is "already typed".
- **Actual state of the .tex** (verified by direct reading of `1_preliminaries.tex` lines 1–1079 and `2_methods.tex` lines 1–1976):
  - `1_preliminaries.tex:1048` reads: "the periodic Heisenberg chain $H_\mathrm{Heis}$ \eqref{eq:H_heis} is given in \hyperref[circ:trotter-strang]{...}". **This is a `\eqref{}`, not a `\label{}`**.
  - I searched lines 1–1079 of `1_preliminaries.tex` and lines 1–1976 of `2_methods.tex` for `\label{eq:H_heis}` and there is no such label. The Heisenberg Hamiltonian is never written as a labelled display in either file. v4 assumes it exists; it doesn't.
- **Why it matters**: v4 introduces this `(Eq. eq:H_heis)` cross-reference at line 86 (`"The 1D periodic Heisenberg chain (Eq. eq:H_heis) is the canonical instance..."`), expecting the user to type a `\eqref{eq:H_heis}` against an already-typed display. If the user copies v4 verbatim, both v4's cross-reference *and* the existing `1_preliminaries.tex:1048` `\eqref{eq:H_heis}` will compile to `??` until the user separately decides where to define the label. The v4 metadata's "Existing labels referenced" section (line 6) mis-states the situation: this label does **not** exist yet.
- **Suggested fix**: pick one of these three options and execute it as part of the v4 → .tex transcription:
  - (a) v4 itself defines the Heisenberg display under `\label{eq:H_heis}` at line 70–74 of v4. The body already has the unlabelled display
    ```
    H_{\mathrm{Heis}} = J \sum_{i=1}^{n} (X_i X_{i+1} + Y_i Y_{i+1} + Z_i Z_{i+1}) + \sum_{i=1}^{n} h_i Z_i,
    ```
    — change the metadata line 6 from "`eq:H_heis` (the Heisenberg-Hamiltonian display already typed into `1_preliminaries.tex`)" to "`eq:H_heis` (introduced by this draft on its $H_{\mathrm{Heis}}$ display, line 70)" and add the HTML-comment label `<!-- \label{eq:H_heis} -->` after the display. The metadata "New labels introduced: none" claim (line 6) then has to change to "introduces `eq:H_heis`".
  - (b) keep v4's metadata claim, but update the Writing Notes line 129 to "v4 silently switches v3's `eq:heisenberg-hamiltonian` → `eq:H_heis`, which **the user must define manually** before transcribing v4 — `1_preliminaries.tex:1048` already contains a forward-reference `\eqref{eq:H_heis}` that is currently dangling, and v4's body adds another such reference."
  - (c) leave the v4 body as is, but state explicitly in the Writing Notes that the Heisenberg display has to be added to `1_preliminaries.tex` before the trotterization subsection (e.g. as part of a future "Standard models" subsection). This would make v4's *use* of `eq:H_heis` parasitic on a future addition.
- **Recommendation**: option (a). The v4 body is the natural home for the Heisenberg display, since it is the canonical instance of the Principle. Adding the HTML-comment label is one line of work and removes the only remaining dangling reference in the chain.

**Sanity-check note on the same label**: `1_preliminaries.tex:1048` already has `\eqref{eq:H_heis}` predating v4. Whether v4 or some other future addition supplies the label, this is a real existing bug in the .tex source that v4 inherits but does not introduce. Flagging it here so the user is aware before pasting v4 in.

---

## [MINOR] m1. v4 retains `\pm` in the motivational first sentence (v3-review m4 was flagged but unfixed)

- **Location**: v4 line 12, first sentence: "the unitary Hamiltonian evolution $\ee^{\pm \ii H t}$ itself".
- **Issue**: v3-review m4 flagged this `\pm` as a stray (cosmetic) that contradicts the rest of the draft's exclusive use of $\ee^{-\ii H t}$ (lines 24, 30, 41, 51, 53, 60, 70). v4 has not changed line 12, but **also** has not amended the metadata or Writing Notes to acknowledge the stray. The v3-review's two suggestions were either (a) drop the $\pm$ or (b) add a one-sentence note. Neither was applied.
- **Comparison with `1_preliminaries.tex:1034`**: the existing .tex prose also uses $\ee^{\pm \ii H t}$ in the same motivational position, so v4's choice is actually **consistent with the existing source**. This is a softer status than the v3 review described — it is not "stray", it matches the intro sentence already in `1_preliminaries.tex:1034`. Re-classifying from v3-review m4 (cosmetic) to "intentional inheritance, OK".
- **Verdict**: I am downgrading this from MINOR-issue to MINOR-acknowledgement. No fix needed if the user is content with the same convention as `1_preliminaries.tex:1034`. If the user wants to harmonise, the matching prose in `1_preliminaries.tex:1034` would also need editing — and that's a sign-convention pass that the memory file `thesis_sign_convention_trotter.md` already flags as a future TODO. Not a blocker for v4.

---

## [MINOR] m2. Forward-reference `(§ref-to-numerics)` is a placeholder (acknowledged, but worth a louder flag)

- **Location**: v4 line 86: "The explicit prefactor and its empirical comparison with the Suzuki-$\mathcal{S}_4$ data of [CITE: childs2021theory, Fig. 3] are deferred to the numerics chapter (§ref-to-numerics)."
- **Issue**: the placeholder `(§ref-to-numerics)` is not a real label. The v4 Writing Notes (line 128) acknowledge this as something the user must replace later. Flagging it here for visibility because it is the **only forward-reference in the entire draft body** that does not resolve to an existing label, and it is easy to miss while transcribing.
- **Suggested fix**: after the user has chosen / created the numerics-chapter section label, replace `(§ref-to-numerics)` with the real `\ref{sec:numerics-trotter-prefactor}` (or whatever it ends up being). Alternatively, until that section exists, flag the placeholder more prominently — e.g. wrap it as `[CHECK: placeholder]` or use a dedicated `\todo{}` macro so it doesn't quietly slip through.
- **Severity**: cosmetic, but it would print as `??` in a compiled draft until fixed.

---

## [NIT] n1. The Heisenberg paragraph uses both $H_{\mathrm{Heis}}$ and $H^{\mathrm{iso}}_{\mathrm{Heis}}$ implicitly

- **Location**: v4 line 86: "the $K=2$ commuting-class split $H_E + H_O$ along the parity of bond indices yields $\tilde\alpha^{(2)}_{\mathrm{comm}}(H_{\mathrm{Heis}}) = \mathcal{O}(|J|^3 n)$".
- **Issue**: in v3 the Heisenberg derivation distinguished the *isotropic* part $H^{\mathrm{iso}}_{\mathrm{Heis}} := H_{\mathrm{Heis}} - \sum_i h_i Z_i$ from the full $H_{\mathrm{Heis}}$, and the bound $54|J|^3 n$ was for $H^{\mathrm{iso}}_{\mathrm{Heis}}$ specifically. The disorder contribution was added separately (v3 line 117–123). v4 collapses this distinction into one sentence: "$H_E + H_O$ ... yields $\tilde\alpha^{(2)}_{\mathrm{comm}}(H_{\mathrm{Heis}}) = \mathcal{O}(|J|^3 n)$; the $Z$-disorder field $\sum_i h_i Z_i$ contributes at the same $\mathcal{O}(n)$ order, with a constant prefactor polynomial of total degree three in $(|J|, h)$".
- **Risk**: a reader could read the first clause as saying the bound is $\mathcal{O}(|J|^3 n)$ for the full disordered $H_{\mathrm{Heis}}$, then be confused by the second clause that adds the disorder contribution separately. The clause "with a constant prefactor polynomial of total degree three in $(|J|, h)$" is doing a lot of work — it is a one-line summary of v3's lines 117–123.
- **Fix (one-clause patch)**: Change "yields $\tilde\alpha^{(2)}_{\mathrm{comm}}(H_{\mathrm{Heis}}) = \mathcal{O}(|J|^3 n)$" to "yields $\tilde\alpha^{(2)}_{\mathrm{comm}}(H^{\mathrm{iso}}_{\mathrm{Heis}}) = \mathcal{O}(|J|^3 n)$ for the isotropic part" — then the second clause that brings in disorder reads more naturally. Equivalently, write: "$\tilde\alpha^{(2)}_{\mathrm{comm}}(H_{\mathrm{Heis}}) = \mathcal{O}((|J| + h)^3 n)$ with the constant absorbing the disorder field" if the user prefers a single bound.
- **Severity**: nit. The current phrasing is technically correct but slightly compressed. Either patch makes the dependence cleaner.

---

## [NIT] n2. The Principle's combinatorial constant: $K^{p+1}(2d)^p$ vs the v3-review's suggestion of $(2d)^p \cdot \text{ordering factor}$

- **Location**: v4 line 84: "with $c_{d,p} \le K^{p+1}\,(2d)^p$ a combinatorial constant depending only on the bond-graph degree $d$ and the order $p$".
- **Issue**: the v3-review M1 suggested fix used $c_{d,p} \le (2d)^p \cdot \text{(ordering factor)}$ — the ordering factor is the count of $(p+1)$-tuples with specified position labels, which is a $K^{p+1}$ multiplicative factor coming from the *outer summand-class sum*. v4 collapses both into $K^{p+1}(2d)^p$. This is fine on the surface, but the v4 body sentence directly preceding the bound is "It is the inner bond layer ... that carries the geometric locality assumption". Calling the *bond-level* count $(2d)^p$ but stamping a $K^{p+1}$ multiplier on top can read as if the outer-class loop is still doing the geometric work, which it isn't. The two $\mathcal{O}(1)$ factors are conceptually different (outer class-index loop = "which classes does each $H_{k_j}$ come from?" vs inner pivot-and-connected-tuple loop = "given the classes, which bonds are connected?") and conflating them in one product is a stylistic choice that may slightly muddy the very distinction v4 is trying to make.
- **Sanity-check**: for 1D Heisenberg even–odd ($p=2$, $d=2$, $K=2$): $K^{p+1}(2d)^p = 2^3 \cdot 4^2 = 128$. v3 derivation showed actual count $= 4$ per pivot per parity $\times 2$ parities $\times \Theta(n) = \Theta(n) \times 8 \le 128 \cdot \Theta(n) / n$, so the bound is loose but valid, as expected.
- **Fix**: optional. If you want to preserve the conceptual separation, write "$c_{d,p} \le K^{p+1}\,(2d)^p$, where the $K^{p+1}$ is the count of class-index orderings ($\mathcal{O}(1)$ since $K = \mathcal{O}(1)$) and the $(2d)^p$ is the bond-pivot connected-tuple count". One extra clause; still clean.
- **Severity**: nit. Both are valid upper bounds; the issue is purely expository.

---

## [NIT] n3. Suzuki recursion paragraph: "specialised to Strang" → tense ambiguity

- **Location**: v4 line 35: "the numerics of Part IV use $p=2$ throughout, so the displayed error analysis below is specialised to Strang."
- **Issue**: this sentence reads OK but "specialised to Strang" is a slight oddity — the analysis isn't specialised to Strang per se, it's *displayed* in Strang form to match the Part IV numerics. The Chapter 5 props (Prop. prop:trotter-diss, Prop. prop:trotter-B) are general-$p$ statements, so "the displayed error analysis below" is a localised reference to the explicit Eq. eq:strang-explicit-error display, not the broader thesis chain.
- **Suggested rephrase**: "the displayed error analysis in (Eq. eq:strang-explicit-error) is the $p=2$ case of (Eq. eq:childs-comm-bound) needed for the Part IV numerics; the Chapter 5 propositions remain stated for general $p$." Slightly more verbose but unambiguous.
- **Severity**: pure cosmetic. Skip if you want to preserve the current pacing.

---

## Orphan-reference sweep

Verified: no surviving back-references in v4 to any of the v3 labels that v4 dropped:

| Dropped label (from v3) | Searched in v4 body | Result |
|---|---|---|
| `eq:trotter-M-count-heisenberg` | no occurrences in v4 lines 12–101 | Clean ✓ |
| `eq:heisenberg-alpha2` | no occurrences | Clean ✓ |
| `eq:heisenberg-evenodd-split` | no occurrences | Clean ✓ |
| `eq:tfim-hamiltonian` | no occurrences | Clean ✓ |
| `eq:heisenberg-2d-hamiltonian` | no occurrences | Clean ✓ |
| `eq:trotter-general` | no occurrences | Clean ✓ (was self-referential in v3) |
| `eq:strang` | no occurrences | Clean ✓ (kept as comment-label `<!-- \label{eq:strang} -->` on the display itself, but never back-referenced) |
| `eq:heisenberg-hamiltonian` | no occurrences | Clean ✓ (renamed to `eq:H_heis`, see B1) |
| `eq:strang-prefactor` | no occurrences | Clean ✓ (renamed to `eq:strang-explicit-error`) |
| `eq:trotter-count` | no occurrences | Clean ✓ (renamed to `eq:trott-num-bound`) |

All harmonised labels (`eq:H_heis`, `eq:strang-explicit-error`, `eq:trott-num-bound`) verified by direct reading of `1_preliminaries.tex`:
- `eq:strang-explicit-error` appears at `1_preliminaries.tex:1062` ✓ (matches v4 metadata)
- `eq:trott-num-bound` appears at `1_preliminaries.tex:1071` ✓ (matches v4 metadata)
- `eq:H_heis` **does not appear anywhere** as a `\label{}` in `1_preliminaries.tex` or `2_methods.tex` — see B1 above.

---

## Citation-table pruning sweep

Verified each row in v4 lines 113–121 against the body:

| Row | Body location | PDF location verified |
|---|---|---|
| Thm 6 (Eq. 44) | line 39 ("Specialising Theorem 6...") | p.12, Eq. (44) — anti-Hermitian form ✓ |
| Cor. 7 (Eq. 45) | line 59 ("Trotter-number bound of [CITE..., Cor. 7]") | p.13, Eq. (45) ✓ |
| Prop. 9 (Eq. 120) | line 49 ("triangle-inequality bootstrap of Prop. 9") | p.22, Eq. (120) ✓ |
| Prop. 10 (Eq. 121) | line 49 ("Strang analysis of Prop. 10") | p.22, Eq. (121) ✓ |
| Eq. 122 (OBC Heis.) | line 76 ("Childs et al.\ analyse an OBC variant ... [CITE, Eq.~(122)]") | p.23, Eq. (122) ✓ — flagged as OBC vs v4's PBC, with $h_j \in [-1,1]$ in Childs vs v4's $h_i \in [-h, h]$ |
| Fig. 3 | line 86 ("empirical comparison with the Suzuki-$\mathcal{S}_4$ data of [CITE..., Fig. 3]") | p.24, paragraph after Fig. 3 ✓ |
| Sec. IV C | line 84 ("the all-to-all TFIM analysed in [CITE..., Sec. IV C]") | p.19, Eq. (82) defines all-to-all TFIM; p.20, Eq. (93) shows $\|\mathrm{ad}_A^2 B\| = \mathcal{O}(n^3 j^3)$ ✓ |

All 7 rows match the body and the PDF. Only `childs2021theory` is cited; no new bibtex entries are introduced.

**Minor table observation**: row 6 says "(referenced as the empirical-comparison anchor for the numerics chapter)". The body actually mentions Fig. 3 once at line 86 ("the Suzuki-$\mathcal{S}_4$ data of ... [CITE..., Fig. 3] are deferred to the numerics chapter"). The table description is slightly stronger ("anchor") than the body warrants ("deferred to"); cosmetic.

---

## Items confirmed accurate (preserve during revision)

- **B1 (v3) fully fixed by deletion**: the explicit numerical constant $\sqrt{15\cdot 54}$ / $15^{1/4}\sqrt{54}$ is gone from the body. The v4 Consequence section (lines 88–101) now has only the $\mathcal{O}$-form, with the symbolic substitution $\tilde\alpha^{(2)}_{\mathrm{comm}} = \mathcal{O}(n)$ producing the $\mathcal{O}(|J|^{3/2}\,n^{1/2}/(\sigma^{3/2}\,\varepsilon^{1/2}))$ M-count cleanly. ✓

- **M1 (v3) fully addressed**: the Principle paragraph (v4 line 80–84) now correctly distinguishes the two layers — outer commuting-class grouping ($K = \mathcal{O}(1)$) and inner bounded-degree bond graph (max degree $d$). The all-to-all TFIM is named explicitly as the converse-failure case ("an all-to-all interaction (e.g. the all-to-all TFIM analysed in [CITE..., Sec. IV C]) has $K=2$ summands but unbounded bond degree $d$, yielding the sharper $\mathcal{O}(n^3)$ commutator constant"). The phrase "two-layer count" appears explicitly. ✓

- **m1 (v3) addressed**: the Hermiticity qualifier "Since each $H_k$ is Hermitian by construction, each factor $\ee^{-\ii t H_k}$ is unitary" appears at v4 line 35. ✓ (Note: the same Hermiticity sentence is now also typed into `1_preliminaries.tex:1048`, so this is a clean match.)

- **m3 (v3) addressed**: 2D Heisenberg paragraph in v3 is now collapsed to a single sentence (v4 line 86 last sentence) with no degree-typo. ✓

- **m4 (v3) inherited from .tex source**: see m1 above.

- **m2 (v3) addressed**: v4 line 35 reads "The Chapter 5 error analysis (Prop. prop:trotter-diss, Prop. prop:trotter-B) is stated for general $p$ and therefore carries the Suzuki recursion implicitly; the numerics of Part IV use $p=2$ throughout..." — explicit motivation for retaining the Suzuki recursion. ✓

- **Logical flow end-to-end**: no jumps. Product formulas → commutator-scaling bound → geometric-local Principle → CKG consequence. The Principle does pull its weight in the Consequence section: the deletion of the explicit Heisenberg arithmetic does not leave a hole because the Principle's $\mathcal{O}(n)$ scaling feeds straight into Prop. prop:trotter-diss substitution, and the in-line $\mathcal{O}(n^3)$ baseline derivation at line 96 ("bounding each nested commutator by $4\,(\max_k\|H_k\|)^3$ with $\|H_k\| = \mathcal{O}(|J|n)$ for each commuting class and summing $K^2 = 4$ ordering terms") supplies the contrast that justifies the $\sqrt{n}$ payoff. ✓

- **Density** (v3 was too dense, v4 trim risk over-trim): v4 reads at the right level — denser than typical preliminary prose but not encyclopedic. The Principle paragraph is a single dense paragraph but does include a worked sanity-check (the parenthetical "$K=2$ summands but unbounded bond degree $d$ ... $\mathcal{O}(n^3)$" example), so even a careful first-time reader can confirm the count works. The Consequence section's $\sqrt{n}$ vs $n^{3/2}$ contrast lands cleanly in the last paragraph. ✓

- **Markdown discipline**: v4 follows the `.claude/rules/thesis-writing.md` style — `#`/`##` headers, `$math$`/`$$math$$`, HTML-comment labels, `[CITE: bibkey]` annotations. No `\subsection{}` or LaTeX dumps. ✓

- **Algebra spot-checks**:
  - $54 = 4 \cdot 108 \cdot (1/12 + 1/24) = 4 \cdot 108 \cdot 1/8 = 54$ ✓ (the v3 prefactor that has been deleted from v4 body but lives in Writing-Notes line 130; arithmetic re-verified for completeness)
  - $K^{p+1}(2d)^p$ for $K=2, d=2, p=2$: $8 \cdot 16 = 128$, valid loose upper bound for the actual count of $4$ per pivot $\times 2$ parities $= 8$ per pivot ✓
  - $\sqrt{n} \leftrightarrow n^{3/2}$ exponents under the Cor. 7 substitution: $M = \mathcal{O}(\sqrt{\tilde\alpha^{(2)}_{\mathrm{comm}} \cdot t^3 / \varepsilon})$, so $\tilde\alpha^{(2)}_{\mathrm{comm}} = \mathcal{O}(n^3) \Rightarrow M = \mathcal{O}(n^{3/2})$ vs $\tilde\alpha^{(2)}_{\mathrm{comm}} = \mathcal{O}(n) \Rightarrow M = \mathcal{O}(n^{1/2})$ ✓

- **Citation count audit**: the body cites only `childs2021theory`, exactly as the metadata claims. The Writing-Notes claim "No new citations introduced" verified ✓.

---

## Summary Scorecard

| Category | Rating | Notes |
|----------|--------|-------|
| Logical correctness | OK | All explicit constants now in $\mathcal{O}$-form; B1 (v3) and M1 (v3) both correctly resolved. Two-layer counting argument is precise. |
| Notation consistency | OK | $\tilde\alpha^{(2)}_{\mathrm{comm}}$, $\tilde\alpha_{\mathrm{comm}}$, $H_k$, $h_b$, $K$, $d$ all used consistently across v4. The motivational $\ee^{\pm \ii Ht}$ in line 12 matches `1_preliminaries.tex:1034` (intentional). |
| Citations | OK | Table fully verified against PDF; 7 rows, 7 verified. Only `childs2021theory` cited; no new bibtex entries. |
| Thesis consistency | **MINOR** | `eq:H_heis` is referenced but not defined in `1_preliminaries.tex` (B1 above). v4 inherits this from the existing .tex which already has a dangling `\eqref{eq:H_heis}`, but v4's metadata claim "already typed into `1_preliminaries.tex`" is wrong. Easy fix (one HTML-comment label on the v4 body display). |
| Exposition quality | OK | Density is appropriate. Suzuki retention is now motivated. Principle paragraph does the load-bearing work that the deleted Heisenberg arithmetic used to do. |

---

## Priority-ordered fix list

**Must fix before merging (blocker):**

1. **B1**: define `eq:H_heis`. Recommended: put `<!-- \label{eq:H_heis} -->` on the v4 body's $H_{\mathrm{Heis}}$ display (after the comma at line 73, on line 74 after the `\qquad X_{n+1} \equiv X_1, \ldots,` close), and update v4 metadata line 6 + Writing Notes line 129 to say "v4 introduces `eq:H_heis`". This *also* fixes the pre-existing dangling `\eqref{eq:H_heis}` in `1_preliminaries.tex:1048` once v4 is transcribed in.

**Should fix (minor):**

2. **m2**: replace the placeholder `(§ref-to-numerics)` at v4 line 86 with either a real `\ref{...}` or a more visible `[CHECK]` flag so it doesn't compile to `??`.

**Optional polish (nit):**

3. **n1**: add "for the isotropic part" or similar in v4 line 86 to disambiguate $H_{\mathrm{Heis}}$ vs $H^{\mathrm{iso}}_{\mathrm{Heis}}$ in the Heisenberg-as-canonical-instance paragraph.
4. **n2**: optionally split $K^{p+1}(2d)^p$ into its outer/inner contributions in the Principle paragraph, with a short clause naming each factor.
5. **n3**: optionally rephrase "the displayed error analysis below is specialised to Strang" in v4 line 35.

---

## Verdict

**Merge-ready after one [BLOCKER] fix.**

The B1 fix is a single-line edit (add `<!-- \label{eq:H_heis} -->` on v4 line 74, plus a small metadata update). After that, v4 is the cleanest pass yet: every Childs citation is verified to the PDF page, every label-name harmonisation matches the existing `.tex`, the M1 refactor correctly distinguishes summand-graph from bond-graph locality, and the deleted explicit-constant arithmetic (B1 from v3) is replaced by the symbolic Principle without losing any of the conceptual content.

The trim is well-judged. v3 was too dense and over-detailed for a preliminaries section; v4 keeps exactly the load-bearing pieces (commutator-scaling bound, Strang tight prefactor, Trotter-number count, Principle, Consequence) and defers all numerical specifics to Part IV via the symbolic constant $\tilde\alpha^{(2)}_{\mathrm{comm}}$. The downstream consumer in `2_methods.tex` (Prop. prop:trotter-diss, Prop. prop:trotter-B) consumes the symbolic constant directly, so v4's choice keeps the chain consistent.

If the user opts for fix-suggestion (a) for B1, the v4 → .tex transcription is now: copy the body, add the HTML-comment label, edit the metadata to claim the new label, fix the placeholder `(§ref-to-numerics)`. Total effort: under 5 minutes.
