# Linear Combination of Unitaries

> **Insertion target:** `supplementary-informations/1_preliminaries.tex`, body of `\subsection{Linear Combination of Unitaries.}\label{sec:prelim-LCU}`, between the `\subsection{...}\label{...}` line at 1037 and the `\begin{figure}` at 1039. The circuit figure `circ:lcu` (lines 1039–1064) stays immediately below the body.
> **Existing labels referenced:** `circ:lcu`, `sec:prelim-QSP`, `eq:B-block-encoding`, `eq:b_minus`, `eq:b_plus-s-eta`.
> **New label introduced:** `eq:lcu-block-encoding`, `eq:lcu-success-prob`.
> **New citation keys introduced:** `lowChuang2019qubitization` (bibtex stub at the bottom). All other keys (`childs2012hamiltonian`, `berry2015simulating`, `chen2023efficient`) already appear in the thesis's `references.bib`.

---

Trotterization gives us Hamiltonian evolution; phase estimation gives us energy labels. The remaining circuit primitive needed for Chapter 5's coherent correction $B$ is a way to realise a *weighted sum of unitaries* as a single deterministic operation, at a controlled cost in ancilla qubits and success probability. This is the \textit{linear combination of unitaries} (LCU) scheme of [Childs and Wiebe 2012] [CITE: childs2012hamiltonian] and [Berry et al. 2015] [CITE: berry2015simulating], which we recall here in its modern block-encoding form [CITE: lowChuang2019qubitization] because that is the interface on which the QSP subsection (see §sec:prelim-QSP) builds.

Fix a unitary family $\{U(t)\}_{t\in S}$ indexed by a discrete grid $S$ and a scalar weight $f : S \to \mathbb{R}_{\geq 0}$ with finite $\ell_1$ norm $\|f\|_1 = \sum_t f(t)$. The target operator
$$
A \;=\; \sum_{t \in S} f(t)\,U(t)
$$
is in general not unitary; its operator norm is at most $\|f\|_1$ by the triangle inequality and typically well below it. The LCU circuit in (Fig. circ:lcu) realises $A$ up to a known rescaling using one *block ancilla* register of $b = \lceil \log_2|S|\rceil$ qubits and three ingredients:

1. $\textsc{Prep}$ on the ancilla, loading the filter state $\ket{f} := \frac{1}{\sqrt{\|f\|_1}}\sum_{t \in S} \sqrt{f(t)}\,\ket{t}$;
2. $\textsc{Select}$ on system$\otimes$ancilla, the controlled family $\sum_{t}\ket{t}\!\bra{t}\otimes U(t)$;
3. $\textsc{Prep}^\dagger$, which uncomputes the ancilla.

Writing the composite circuit as $U_A := (\textsc{Prep}^\dagger \otimes \idm)\, \textsc{Select}\, (\textsc{Prep} \otimes \idm)$, a direct calculation shows that its $\ket{0^b}$–$\ket{0^b}$ matrix element on the ancilla is exactly
$$
\bra{0^b}\,U_A\,\ket{0^b} \;=\; \frac{A}{\|f\|_1},
$$ <!-- \label{eq:lcu-block-encoding} -->
i.e.\ $U_A$ is a *block encoding* of $A$ with sub-normalisation $\alpha := \|f\|_1$ [CITE: lowChuang2019qubitization]. The sub-normalisation is the price paid for packing a non-unitary $A$ inside a unitary: $U_A$ acts as $A/\alpha$ on the block projector $\ket{0^b}\!\bra{0^b}$ and as an unspecified completion on its orthogonal complement, and all downstream complexity statements (QSP polynomial degree, simulation time, success probability) scale with $\alpha$ rather than with $\|A\|$.

If $f$ takes signed or complex values, the circuit structure is unchanged: absorb the phase into \textsc{Select} by redefining $\tilde U(t) := \mathrm{sgn}\bigl(f(t)\bigr)\,U(t)$ with $\mathrm{sgn}\bigl(f(t)\bigr) := f(t)/|f(t)|$ a unit-modulus phase, and let $\textsc{Prep}$ load the amplitude-magnitude state $\ket{f} = \frac{1}{\sqrt{\|f\|_1}}\sum_t \sqrt{|f(t)|}\,\ket{t}$; the sub-normalisation $\alpha = \|f\|_1 = \sum_t |f(t)|$ is now the usual $\ell_1$ norm. This is exactly the convention used for the signed kernels $b_-$ and $b_+$ of Chapter 5 (see \eqref{eq:b_minus}, \eqref{eq:b_plus-s-eta}).

Running $U_A$ on $\ket{0^b}\otimes\ket{\psi}$ and *post-selecting* the ancilla on $\ket{0^b}$ leaves, up to normalisation, the state $A\ket{\psi}/\alpha$ on the system, with success probability
$$
p_{\mathrm{succ}} \;=\; \frac{\|A\ket{\psi}\|^2}{\alpha^2}.
$$ <!-- \label{eq:lcu-success-prob} -->

When $\alpha \leq 1$, the block encoding can double as a deterministic implementation of $A$ (since then $p_{\mathrm{succ}}\geq\|A\ket{\psi}\|^2$ and no sub-normalisation amplification is needed). In the typical case $\alpha > 1$, however, post-selection succeeds only with probability $O(1/\alpha^2)$ and LCU is *not* deterministic. Two canonical remedies exist. *(i) Repeat-until-success.* Running the circuit $O(\alpha^2)$ times in expectation yields one useful sample; this is acceptable when $A$ is applied to a sampling subroutine whose outputs are anyway averaged, but quadratic in $\alpha$ and incoherent. *(ii) Oblivious amplitude amplification.* [Berry et al. 2015] [CITE: berry2015simulating] observed that when $A/\alpha$ has nearly unit norm on the relevant subspace --- which happens whenever $A$ is close to a unitary, as is the case for the Hamiltonian-simulation Taylor truncation and for our Hermitian $B$-block --- a Grover-style reflection on the block ancilla boosts the success probability to $\Theta(1)$ using only $O(\alpha)$ invocations of $U_A$ and $U_A^\dagger$. The quadratic overhead of (i) is replaced by a linear overhead, and the boosted circuit is itself a unitary acting coherently on system$\otimes$ancilla.

The $\alpha$-scaling is decisive for the coherent correction $B$ of Chapter 5. The operator $B$ is built as a *nested* LCU (see the block-encoding relation \eqref{eq:B-block-encoding}) with two time registers $T_-, T_+$ and sub-normalisations $\alpha_\pm = \|b_\pm\|_1$, so the composite sub-normalisation is
$$
\alpha \;=\; \|b_-\|_1\,\|b_+\|_1 \;=\; \tilde{O}\!\left(\log\frac{\beta\|H\|}{\varepsilon}\right)
$$
in the regularised Metropolis case [CITE: chen2023efficient], i.e.\ polylogarithmic in the problem parameters. The na\"ive success probability $1/\alpha^2$ would therefore be a $\mathrm{polylog}$ factor only, but the $\Theta(1)$-success oblivious amplitude amplification of [CITE: berry2015simulating] keeps that factor coherent and linear rather than quadratic, and is what makes the $B$-block encoding efficient rather than exponentially costly once ported onto the GQSP walk of §sec:prelim-QSP.

The output of the LCU circuit of (Fig. circ:lcu), specialised to the nested two-register kernel construction of Chapter 5, is therefore a unitary $U_B$ acting on system $\otimes T_- \otimes T_+$ with
$$
\bra{0^b}\,U_B\,\ket{0^b} \;=\; \frac{B}{\alpha},
$$
where $\ket{0^b}$ denotes the joint ground state of the $T_-, T_+$ block ancillas and $\alpha = \|b_-\|_1\,\|b_+\|_1$. It is $U_B$, together with a reflection $R_b = 2\ket{0^b}\!\bra{0^b} - \idm$ on the block register, that enters the qubitization walk and the generalised quantum signal processing of §sec:prelim-QSP, where a Laurent-polynomial transformation of $U_B$'s singular values turns the block encoding of $B/\alpha$ into a $\delta$-step of coherent evolution $\ee^{-\ii\delta B}$.

---

## Citations

- **`childs2012hamiltonian`** — Childs, A. M. and Wiebe, N., *"Hamiltonian simulation using linear combinations of unitary operations"*, **Quantum Inf. Comput. 12, 901–924 (2012)**, arXiv:1202.5822. *Already in `references.bib`.*
- **`berry2015simulating`** — Berry, D. W., Childs, A. M., Cleve, R., Kothari, R. and Somma, R. D., *"Simulating Hamiltonian dynamics with a truncated Taylor series"*, **Phys. Rev. Lett. 114, 090502 (2015)**, arXiv:1412.4687. Oblivious amplitude amplification for LCU. *Already in `references.bib`.*
- **`chen2023efficient`** — Chen, C.-F., Kastoryano, M. J. and Gilyén, A., *"An efficient and exact noncommutative quantum Gibbs sampler"*, arXiv:2311.09207 (2023). Nested LCU block encoding of $B$, Proposition III.1. *Already in `references.bib`.*
- **`lowChuang2019qubitization`** — Low, G. H. and Chuang, I. L., *"Hamiltonian simulation by qubitization"*, **Quantum 3, 163 (2019)**, arXiv:1610.06546. Modern block-encoding and sub-normalisation framing (Def. 1, Lemma 2).

```bibtex
@article{lowChuang2019qubitization,
  author  = {Low, Guang Hao and Chuang, Isaac L.},
  title   = {Hamiltonian Simulation by Qubitization},
  journal = {Quantum},
  volume  = {3},
  pages   = {163},
  year    = {2019},
  doi     = {10.22331/q-2019-07-12-163},
  eprint  = {1610.06546},
  archivePrefix = {arXiv},
  primaryClass  = {quant-ph},
}
```

---

## Writing notes

- **Length check.** Body text is 6 paragraphs plus 4 short display equations (two labelled, two unlabelled). At the thesis's line-spacing and font size this sits comfortably under 2 pages; the figure `circ:lcu` then takes roughly another half-page. The three-way structure (primitive + success probability + $\alpha > 1$ remedies) is the one explicitly requested, and the forward link to $U_B$ in the final paragraph closes the loop with §sec:prelim-QSP.
- **Notation harmony.** `\ket{0^b}` and `\bra{0^b}` throughout, matching both the figure `circ:lcu` caption at `1_preliminaries.tex:1042–1047` and the QSP figure caption at `1_preliminaries.tex:1090–1095`. `\textsc{Prep}` / `\textsc{Select}` small-caps match the CKG-style figure caption already in the thesis. The block-encoding relation $\bra{0^b} U_B \ket{0^b} = B/\alpha$ is stated identically to the expression in (Fig. circ:gqsp) caption and in \eqref{eq:B-block-encoding}, so the reader can trace $U_B$ across the three subsections without notational translation.
- **Signed-$f$ display.** The user's spec said the signed/complex $f$ case "can live in a single sentence"; I chose to give it its own full sentence (not a display equation), because it is a notational correction and not a new theorem, and because Chapter 5's $b_\pm$ kernels rely on exactly this convention. If the user prefers it demoted to a parenthetical, the sentence beginning "If $f$ takes signed or complex values..." can be compressed by half.
- **Scaling caveat.** [CHECK] The user wrote "$\tilde O(\beta \log(1/\varepsilon))$" for the nested $\alpha$. In the regularised Metropolis case, `2_methods.tex:1726` gives $\alpha = \mathcal{O}(\log(\beta\|H\|\|{\sum_a A^{a\dagger}A^a}\|/\varepsilon))$, i.e.\ the $\beta$ and $\log(1/\varepsilon)$ factors enter *inside* a single logarithm, not as a product. I used $\tilde O(\log(\beta\|H\|/\varepsilon))$ to match the existing thesis text; if the user really wants the "$\tilde O(\beta\log(1/\varepsilon))$" shorthand (which is correct only if $\beta$ is treated as a separate polylog factor coming from the Gaussian $\|b_-\|_1 \leq (\pi/\sqrt 2)\,\ee^{\beta^2\sigma^2/8}$, not from the Metropolis $\|b_+\|_1$), the adjustment is a one-line edit. Either form is self-consistent; the thesis's own `2_methods.tex` uses the $\log(\beta\|H\|/\varepsilon)$ form.
- **What I did not cite.** I deliberately did *not* introduce `gilyen2019quantum` here, although it is the canonical modern block-encoding reference. The thesis already cites it once, at `2_methods.tex:1916`, for QSVT. Introducing it here would duplicate a framing citation that `lowChuang2019qubitization` already covers, and the thesis body in (§sec:prelim-QSP) would be the natural place to pull it in. Flagged in Suggestions below.
- **Oblivious amplitude amplification provenance.** [CHECK] The original clean statement of oblivious amplitude amplification for LCU is in [Berry et al. 2015] (Lemma 3.6 / Sec.~3). It was already implicit in [Berry, Childs, Cleve, Kothari, Somma 2014] (the "exponential improvement" precursor, arXiv:1312.1414), but the 2015 PRL is the accepted primary reference and is the bibkey already in the thesis, so I cited that one. If the user's `references.bib` also has the 2014 companion, it can be added in the "[CITE: berry2015simulating]" slot as a second key without changing the prose.

---

## Suggestions

*Items not requested by the user but worth considering for a polished thesis.*

1. **Numbered display for signed-$f$.** Consider promoting the signed-$f$ remark to a short numbered remark (`\begin{remark}[Signed and complex $f$]...\end{remark}`) or to its own display equation showing the replacement $\tilde U(t) = \mathrm{sgn}(f(t))\,U(t)$. This would make it cross-referenceable from Chapter 5, where `2_methods.tex:1714–1718` uses the same identity without a pointer back to the preliminaries. Low effort, and it gives (§sec:prelim-LCU) a second `\label` for Chapter 5 to cite.
2. **Gilyén–Su–Low–Wiebe 2019 citation.** The modern block-encoding formalism (Lemma 48, the `(\alpha, a, \varepsilon)`-block-encoding definition used throughout the literature) is cleanest in [Gilyén, Su, Low, Wiebe 2019] rather than [Low, Chuang 2019] --- the former is cited already at `2_methods.tex:1916`. A one-line enrichment would change "recall here in its modern block-encoding form [CITE: lowChuang2019qubitization]" to "...form [CITE: lowChuang2019qubitization, gilyen2019quantum]". Optional and stylistic.
3. **Explicit statement that QSP consumes only $U_B$ and $R_b$.** The final paragraph already names $U_B$ and mentions the reflection $R_b$. If the next subsection (§sec:prelim-QSP) is already written, a one-sentence preview naming $W = R_b\,U_B$ as the *walk operator* would make the handoff seamless; this is exactly the language of the figure `circ:gqsp` caption. Trivial to add; whether it reads as redundant depends on how the QSP prose begins.
4. **Forward links to Chapter 5 $B$-block encoding.** I link to \eqref{eq:B-block-encoding}, \eqref{eq:b_minus}, \eqref{eq:b_plus-s-eta}. A richer forward reference would also point at the nested-LCU derivation at `2_methods.tex:1707–1737` --- but that derivation does not yet carry a labelled block (the `alg:coh` algorithm environment does). If the user adds a label to the paragraph beginning "The block encoding of the coherent term..." at `2_methods.tex:1707` (e.g.\ `\label{sec:B-nested-LCU}`), the final paragraph of this subsection can cite it directly rather than leaving a bare reference to \eqref{eq:B-block-encoding}.
5. **Post-selection vs coherent use.** The prose distinguishes repeat-until-success from amplitude amplification, but does not explicitly say that the latter is the one used in Chapter 5. A half-sentence --- "Chapter 5's GQSP walk consumes $U_B$ coherently, so the relevant regime is (ii)" --- would make the pedagogical flow tighter. Low-risk addition; see paragraph 5.
6. **Scaling clarification** (already flagged as [CHECK]). If the user prefers the shorthand "$\tilde O(\beta \log(1/\varepsilon))$", rewrite the last paragraph's display as $\alpha = \|b_-\|_1\|b_+\|_1 = \tilde{O}(\beta \log(1/\varepsilon))$ with a footnote pointing at the precise $\log(\beta\|H\|/\varepsilon)$ form of [CITE: chen2023efficient]. Choice is purely stylistic.
