# Trotterization

> **Insertion target:** `supplementary-informations/1_preliminaries.tex`, body of `\subsection{Trotterization.}` inside `\chapter{Quantum Circuits}` — between the figure `circ:trotter-strang` (ends at line 1034) and `\subsection{Linear Combination of Unitaries.}` at line 1079. The subsection label `sec:prelim-trotter` at line 968 is already referenced by `2_methods.tex:887`.
> **Figure referenced:** `circ:trotter-strang` (line 1033). The circuit is **not** modified; the draft only cites it.
> **Existing labels referenced:** `sec:prelim-trotter`, `circ:trotter-strang`, `eq:paulis`, `rem:palindromic`, `eq:childs-comm-bound` (already in `1_preliminaries.tex:1053`), `eq:strang-explicit-error` (already at `1_preliminaries.tex:1062`), `eq:trott-num-bound` (already at `1_preliminaries.tex:1071`), `eq:childs-trotter-bound` at `2_methods.tex:888`.
> **New labels introduced:** `eq:H_heis` — the Heisenberg-Hamiltonian display below is its definition. Note: `1_preliminaries.tex:1048` already contains a forward `\eqref{eq:H_heis}` (currently dangling); typing this draft in supplies the missing definition. (v3 had introduced `eq:trotter-general`, `eq:strang`, `eq:childs-comm-bound`, `eq:strang-prefactor`, `eq:trotter-count`, `eq:heisenberg-hamiltonian`, `eq:heisenberg-evenodd-split`, `eq:heisenberg-alpha2`, `eq:trotter-M-count-heisenberg`, `eq:tfim-hamiltonian`, `eq:heisenberg-2d-hamiltonian`. v4 deletes the explicit Heisenberg-arithmetic displays and the TFIM / 2D-Heisenberg displays, so most of the new-label list collapses; only `eq:H_heis` survives.)
> **User action required in `2_methods.tex`:** add `\label{sec:hamiltonian-trotterization}` on the line immediately after `\subsection{Hamiltonian Trotterization.}` at `2_methods.tex:881`. Without this, the `\ref{sec:hamiltonian-trotterization}` invocation below prints `??`.
> **New citation keys introduced:** none. `childs2021theory` is already in `references.bib`.
>
> **TODO — propagate notation to `2_methods.tex`.** This draft frees $k$ for $k$-locality by adopting Childs's notation verbatim: summand index $\ell$, summand count $\Gamma$, Suzuki recursion depth $r$. The Hamiltonian-Trotterization section in `2_methods.tex` (`\subsection{Hamiltonian Trotterization.}` at line 881) still uses the old convention. Concrete patches required there:
> - line 894–896: the master commutator-bound display `\sum_{k_1, k_2, \ldots, k_{p+1} = 1}^{K} \|[H_{k_{p+1}}, \cdots [H_{k_2}, H_{k_1}]\cdots]\|` → rename indices and count: $k_i \to \ell_i$, $K \to \Gamma$.
> - line 918–924 (`eq:childs-trotter-bound`, the Strang prefactor): $\sum_{k_1 < k_2}$ and the $H_{k_1}, H_{k_2}$ subscripts inside → $\ell_1, \ell_2$.
> - line 1031–1037: the same Strang prefactor invoked again inside the dissipator-error proof (`\sum_{k_1 < k_2}\ldots H_{k_2}, [H_{k_2}, H_{k_1}]\ldots`) → same rename.
> - any "$H = \sum_{k=1}^K H_k$" decomposition statement around the Trotter-section preamble → $\sum_{\ell=1}^\Gamma H_\ell$.
> - Suzuki-recursion mentions (if the depth letter $k$ appears in $S_{2k}$, $u_k$, $4^{1/(2k-1)}$, $\mathcal{O}(t^{2k+1})$) → switch to $r$.
> - Algorithm pseudocode at line 1385 (`\For{$\ell = 1, \ldots, L$}`) is *unaffected*: that $\ell$ is already the outer Lindbladian-step counter, scope-local to the algorithm box, and does not collide with the Trotter summand index introduced here.
> - Leave the QPE-register bit index $k$ in `controlled-$U^{2^k}$` displays alone — that is standard QPE notation and does not collide with $k$-locality (it indexes register bits, not Hamiltonian terms).

---

Up to this point the quantum circuits we discussed performed either read-out of an energy register (phase estimation) or state preparation on an ancilla. The dissipative generators of Chapter 5 require a third primitive: the unitary Hamiltonian evolution $\ee^{\pm \ii H t}$ itself, used inside the controlled cascade of (Fig. circ:qpe) — each controlled-$U^{2^k}$ is implemented via Trotter steps of $S_p(t_0)$. For the local Hamiltonians we care about — each term acting on a bounded number of qubits — there is no known shallow exact circuit implementing $\ee^{\pm \ii H t}$; we approximate it by a *product formula*, the subject of this subsection. The error analysis follows [Childs et al. 2021] [CITE: childs2021theory], whose commutator-scaling bound is the main analytical handle. The material here is a preliminary for the sharpened Trotter error bound of (§sec:hamiltonian-trotterization), where the generic commutator constant $\tilde\alpha^{(2)}_{\mathrm{comm}}$ enters Prop. prop:trotter-diss and Prop. prop:trotter-B as a black box; an explicit Heisenberg-specific evaluation is deferred to the numerics of Part IV.

## Product formulas

Write the Hamiltonian as a sum of $\Gamma$ non-commuting Hermitian summands
$$
H \;=\; \sum_{\ell=1}^{\Gamma} H_\ell,
$$

where each $H_\ell$ is either a single Pauli string (Eq. eq:paulis) or a commuting family of such strings whose exponential $\ee^{-\ii t H_\ell}$ is easy to implement (typically because all terms in $H_\ell$ share support disjointly or pairwise commute). A *$p$-th order product formula* $S_p(t)$ is an ordered product of the single-summand evolutions $\ee^{-\ii t H_\ell}$ that approximates $\ee^{-\ii H t}$ up to error $\mathcal{O}(t^{p+1})$ for small $t$ [CITE: childs2021theory]. The two formulas that appear throughout the thesis are

$$
S_1(t) \;:=\; \ee^{-\ii t H_\Gamma}\cdots \ee^{-\ii t H_1}
\qquad\text{(Lie--Trotter, $p=1$)},
$$

and the symmetric Strang splitting
$$
S_2(t) \;:=\; \Bigl(\prod_{\ell=1}^{\Gamma} \ee^{-\ii (t/2) H_\ell}\Bigr)\Bigl(\prod_{\ell=\Gamma}^{1} \ee^{-\ii (t/2) H_\ell}\Bigr),
$$ <!-- \label{eq:strang} -->

with errors $\|S_1(t) - \ee^{-\ii H t}\| = \mathcal{O}(t^2)$ and $\|S_2(t) - \ee^{-\ii H t}\| = \mathcal{O}(t^3)$ respectively. By construction $S_2(t)^\dagger = S_2(-t)$ — Strang is *palindromic* — a structural property that the Chapter 5 adjoint-symmetry argument (Rem. rem:palindromic) relies on to preserve KMS detailed balance under Trotterization. The circuit realisation of $S_2$ for the $n=3$ periodic Heisenberg chain is given in (Fig. circ:trotter-strang): three Pauli-bond groups ($XX+YY+ZZ$ per colour) arranged as A-B-C-B-A, with the central C-block carrying the odd-$n$ wrap-around bond at doubled angle.

Higher-order Suzuki formulas $S_{2r}$ can be built recursively from $S_2$ by the fractal construction $S_{2r}(t) = S_{2r-2}(u_r t)^2\,S_{2r-2}((1-4u_r)t)\,S_{2r-2}(u_r t)^2$ with $u_r = 1/(4 - 4^{1/(2r-1)})$, giving error $\mathcal{O}(t^{2r+1})$ [CITE: childs2021theory]. The Chapter 5 error analysis (Prop. prop:trotter-diss, Prop. prop:trotter-B) is stated for general $p$ and therefore carries the Suzuki recursion implicitly; the numerics of Part IV use $p=2$ throughout, so the displayed error analysis below is specialised to Strang. Since each $H_\ell$ is Hermitian by construction, each factor $\ee^{-\ii t H_\ell}$ is unitary, so $S_p(t)$ and its iterate $S_p^{(M)}(t) := S_p(t/M)^M$ are unitary for every choice of $p$ and $M$ — a property that the CKG error propagation (§sec:hamiltonian-trotterization) relies on when applying Parseval's identity to the Trotterised OFT.

## Commutator-scaling bound

A naive Taylor truncation only controls the error through $\|H_\ell\|$, ignoring commutativity of the summands. The key result of [Childs et al. 2021] [CITE: childs2021theory] is that the error of any $p$-th order formula is in fact controlled by *nested commutators* of the summands. Specialising Theorem 6 of [CITE: childs2021theory] to anti-Hermitian $\ii H_\ell$ and real $t$, one has
$$
\bigl\|S_p(t) - \ee^{-\ii H t}\bigr\|
\;=\; \mathcal{O}\!\bigl(\tilde\alpha_{\mathrm{comm}}\,|t|^{p+1}\bigr),
\qquad
\tilde\alpha_{\mathrm{comm}}
\;:=\; \sum_{\ell_1,\dots,\ell_{p+1}=1}^{\Gamma}
\bigl\|\,[H_{\ell_{p+1}},\,\cdots\,[H_{\ell_2},H_{\ell_1}]\cdots]\,\bigr\|,
$$ <!-- \label{eq:childs-comm-bound} -->

where the asymptotic constant is universal (independent of $\Gamma$, $\|H_\ell\|$, and the ordering fixed by the particular formula). The triangle-inequality bootstrap of [CITE: childs2021theory, Prop. 9 (Eq. 120)] gives the analogous tight first-order prefactor, which the Strang analysis of [CITE: childs2021theory, Prop. 10 (Eq. 121)] bootstraps to the tight second-order form we use below. Following the cleaner notation already adopted in `2_methods.tex:1032`, the latter reads
$$
\bigl\|S_2(t) - \ee^{-\ii H t}\bigr\|
\;\leq\; \sum_{\ell_1 < \ell_2}\!\left(
\frac{t^3}{12}\bigl\|[H_{\ell_2},[H_{\ell_2},H_{\ell_1}]]\bigr\|
\;+\; \frac{t^3}{24}\bigl\|[H_{\ell_1},[H_{\ell_1},H_{\ell_2}]]\bigr\|
\right)
\;+\; \mathcal{O}(t^4).
$$ <!-- \label{eq:strang-explicit-error} -->

The coefficients $1/12$ and $1/24$ match the lowest-order BCH expansion and are tight [CITE: childs2021theory]. Splitting the target evolution $\ee^{-\ii H t}$ into $M$ Trotter steps of size $\delta t = t/M$ and applying the triangle inequality over steps yields the *Trotter-number bound* of [CITE: childs2021theory, Cor. 7]: to achieve total error $\varepsilon$, it suffices to take
$$
M \;\geq\; \Omega\!\left(\frac{\tilde\alpha_{\mathrm{comm}}^{1/p}\;t^{\,1+1/p}}{\varepsilon^{\,1/p}}\right).
$$ <!-- \label{eq:trott-num-bound} -->

For $p=1$ this is the $\mathcal{O}(\tilde\alpha_{\mathrm{comm}}\, t^2/\varepsilon)$ worst-case count, and for Strang ($p=2$) the $\mathcal{O}(\sqrt{\tilde\alpha_{\mathrm{comm}}}\, t^{3/2}/\sqrt{\varepsilon})$ bound that governs the Chapter 5 analysis.

## Geometrically-local Hamiltonians

The Hamiltonian used in the numerical chapter and sketched in (Fig. circ:trotter-strang) is the isotropic periodic Heisenberg chain on $n$ qubits,
$$
H_{\mathrm{Heis}}
\;=\; J \sum_{i=1}^{n} \bigl( X_i X_{i+1} + Y_i Y_{i+1} + Z_i Z_{i+1} \bigr)
\;+\; \sum_{i=1}^{n} h_i Z_i,
\qquad X_{n+1} \equiv X_1,\ \ldots,
$$ <!-- \label{eq:H_heis} -->

with coupling $J>0$ and disorder fields $h_i \in [-h,h]$. Childs et al.\ analyse an OBC variant of this model with $j=1,\ldots,n-1$ in [CITE: childs2021theory, Eq.~(122)]; the PBC extension adds only a constant number of boundary-bond pairs and so preserves the $n$-scaling. Rather than evaluate $\tilde\alpha^{(2)}_{\mathrm{comm}}(H_{\mathrm{Heis}})$ in closed form here, we abstract the structural feature that drives the Chapter 5 bounds — geometric locality — into the following principle, of which the Heisenberg chain is one instance among many.

**Principle: geometrically-local commutator scaling.** A Hamiltonian is *geometrically local* if it can be written as
$$ H \;=\; \sum_b h_b, \qquad \|h_b\| \le h_*, $$
where each *local term* $h_b$ is supported on at most $k$ sites — equivalently, $H$ is *$k$-local* with $k = \mathcal{O}(1)$ — and the *interaction graph* on the local terms — vertex set $\{h_b\}$, edges between terms with overlapping supports — has bounded maximum degree $d$ (each $h_b$ shares support with $\mathcal{O}(d)$ others). Such an $H$ admits a partition into $\Gamma = \mathcal{O}(1)$ commuting classes $H_\ell := \sum_{b\in\mathrm{class}_\ell} h_b$, each class a sum of pairwise-commuting local terms. Expanding the $(p+1)$-nested commutator $[H_{\ell_{p+1}}, \cdots [H_{\ell_2}, H_{\ell_1}]]$ term-by-term over local terms, one obtains a sum over ordered $(p+1)$-tuples of bonds $(b_1,\dots,b_{p+1})$, of which only *connected* tuples — each $b_{i+1}$ overlapping with at least one of $b_1,\dots,b_i$ — survive, since disjoint-support terms commute. The count of connected $(p+1)$-tuples anchored at a pivot bond $b_1$ is at most $(2d)^p$, and there are $\Theta(n)$ choices of pivot bond, giving
$$
\tilde\alpha^{(p)}_{\mathrm{comm}} \;\leq\; c_{d,p}\,n\,h_*^{p+1} \;=\; \mathcal{O}(n),
$$
with $c_{d,p} \le \Gamma^{p+1}\,(2d)^p$ a combinatorial constant depending only on the bond-graph degree $d$ and the order $p$. The argument is a *two-layer* count: outer commuting-class grouping fixes $\Gamma = \mathcal{O}(1)$ and selects which classes contribute, while inner bounded-degree bond counting yields the $\Theta(n)\cdot(2d)^p$ pivot-and-connected-tuple bound. It is the inner bond layer — not the outer summand layer — that carries the geometric locality assumption: an all-to-all interaction (e.g. the all-to-all TFIM analysed in [CITE: childs2021theory, Sec. IV C]) has $\Gamma=2$ summands but unbounded bond degree $d$, yielding the sharper $\mathcal{O}(n^3)$ commutator constant rather than $\mathcal{O}(n)$.

The 1D periodic Heisenberg chain (Eq. eq:H_heis) is the canonical instance of this Principle: its bond set is $2$-colourable for even $n$ (a third group closes the boundary for odd $n$, as carried by the C-block of Fig. circ:trotter-strang), the bond-interaction degree on the chain is $d=2$, and the $\Gamma=2$ commuting-class split $H_E + H_O$ of the isotropic part along the parity of bond indices yields $\tilde\alpha^{(2)}_{\mathrm{comm}}(H^{\mathrm{iso}}_{\mathrm{Heis}}) = \mathcal{O}(|J|^3 n)$; the $Z$-disorder field $\sum_i h_i Z_i$ contributes at the same $\mathcal{O}(n)$ order, with a constant prefactor polynomial of total degree three in $(|J|, h)$, so $\tilde\alpha^{(2)}_{\mathrm{comm}}(H_{\mathrm{Heis}}) = \mathcal{O}(n)$ overall. The explicit prefactor and its empirical comparison with the Suzuki-$\mathcal{S}_4$ data of [CITE: childs2021theory, Fig. 3] are deferred to the numerics chapter (§ref-to-numerics) [CHECK: replace placeholder with the numerics-chapter section label once it exists]. The same Principle applies to the periodic 1D transverse-field Ising model and the 2D Heisenberg model on a square lattice, again with $\tilde\alpha^{(2)}_{\mathrm{comm}} = \mathcal{O}(n)$ scaling and constant prefactors that change only with the bond-interaction degree $d$.

## Consequence for the CKG Lindbladian

The Hamiltonian-Trotterization analysis of (§sec:hamiltonian-trotterization) states its main propositions (Prop. prop:trotter-diss, Prop. prop:trotter-B) in terms of the *generic* constant $\tilde\alpha^{(2)}_{\mathrm{comm}}$ of (Eq. eq:strang-explicit-error), making no structural use of $H$. Substituting the geometrically-local scaling $\tilde\alpha^{(2)}_{\mathrm{comm}} = \mathcal{O}(n)$ into the Strang-splitting dissipator bound
$$
\bigl\|\bar{\mathcal{L}}_{\mathrm{diss}} - \tilde{\mathcal{L}}_{\mathrm{diss}}\bigr\|_{1\to 1}
\;\leq\; \frac{\sqrt{15}\,\tilde\alpha^{(2)}_{\mathrm{comm}}}{M^2\,\sigma^3}
$$

of Prop. prop:trotter-diss sharpens the $n$-dependence from the generic $\mathcal{O}(n^3)$ of a blind triangle inequality on (Eq. eq:childs-comm-bound) — bounding each nested commutator by $4\,(\max_\ell\|H_\ell\|)^3$ with $\|H_\ell\| = \mathcal{O}(|J|n)$ for each commuting class and summing $\Gamma^2 = 4$ ordering terms — to the $\mathcal{O}(n)$ scaling of the Principle. Solving for the number of Strang steps needed to reach dissipator error $\varepsilon$ then gives
$$
M \;=\; \mathcal{O}\!\left(\frac{|J|^{3/2}\,n^{1/2}}{\sigma^{3/2}\,\varepsilon^{1/2}}\right).
$$

The payoff of the specialisation is visible in the $n$-exponent: the required Trotter-step count scales as $\sqrt{n}$ for any geometrically-local model, not the $n^{3/2}$ that a generic $\mathcal{O}(n^3)$ commutator constant would yield. The same substitution propagates through Prop. prop:trotter-B for the coherent term, with no structural change to either proof. Numerical evaluation of the explicit Heisenberg constant — together with its empirical comparison against the simulated dissipator error — is carried out in Part IV; the rest of the thesis keeps $\tilde\alpha^{(2)}_{\mathrm{comm}}$ symbolic so the bounds remain applicable to any decomposition $H = \sum_\ell H_\ell$.

---

## Citations

- **`childs2021theory`** — Childs, A.~M., Su, Y., Tran, M.~C., Wiebe, N., Zhu, S., *"Theory of Trotter Error with Commutator Scaling"*, **Phys. Rev. X 11, 011020 (2021)**, arXiv:1912.08854. Already in `references.bib` and cited throughout `2_methods.tex`.

No new bibtex entries introduced.

**Page-level verification of every Childs citation used above** (against `supplementary-informations/Childs et al. - 2021 - Theory of Trotter Error with Commutator Scaling.pdf`):

| Citation in draft | PDF location | Verified |
|---|---|---|
| Thm 6 (master commutator bound, anti-Hermitian form) | p. 12, Eq. (44) | ✓ |
| Cor. 7 (Trotter-number count) | p. 13, Eq. (45) | ✓ |
| Prop. 9 (first-order tight bound) | p. 22, Eq. (120) | ✓ |
| Prop. 10 (Strang tight bound with 1/12, 1/24) | p. 22, Eq. (121) | ✓ |
| Eq. (122) (OBC Heisenberg + random field) | p. 23 | ✓ (flagged as OBC vs our PBC) |
| Fig. 3 (Suzuki-$\mathcal{S}_4$ even–odd, $n=10$ tightness) | p. 24, paragraph after figure | ✓ (referenced as the empirical-comparison anchor for the numerics chapter) |
| Sec. IV C (all-to-all TFIM, $\mathcal{O}(n^3)$ commutator) | p. 19, Eq. (82) and Eq. (93) | ✓ (cited only as converse to the geometric-locality principle) |

---

## Writing Notes

- **Notation update (v4 → v4-relabel).** v4 used $H = \sum_{k=1}^{K} H_k$ (with a footnote acknowledging Childs's $\Gamma$). To free $k$ for $k$-locality across the thesis, this revision adopts Childs's notation verbatim: summand index $\ell$, summand count $\Gamma$, Suzuki recursion-depth letter $r$ (so $S_{2r}$, $u_r$, $4^{1/(2r-1)}$). The $K \to \Gamma$ rename also applies to every count appearance in the body ($\Gamma=2$ commuting-class split, $\Gamma^2=4$ ordering terms in the Consequence section, $c_{d,p}\le \Gamma^{p+1}(2d)^p$). The footnote about $K$-vs-$\Gamma$ is dropped. The Principle paragraph now explicitly introduces "$k$-local with $k=\mathcal{O}(1)$" so the reserved letter has a clear definitional anchor in the same chapter; this is the only intentional appearance of $k$ in the body. The QPE-register bit $k$ in `controlled-$U^{2^k}$` is *not* renamed (standard QPE notation, indexes register bits not Hamiltonian terms; flagged in the TODO so the reader can confirm it is intentional).
- **Self-orphan-check (post-trim).** Searched the v4 body for back-references to every label I dropped from v3: `eq:heisenberg-evenodd-split` (gone — never back-referenced even in v3), `eq:heisenberg-alpha2` (gone — back-referenced once in the v3 Consequence section, which I rewrote in $\mathcal{O}$-form), `eq:trotter-M-count-heisenberg` (gone — was only self-referential), `eq:tfim-hamiltonian` (gone — never referenced anywhere in `1_preliminaries.tex` or `2_methods.tex`, confirmed by `grep`), `eq:heisenberg-2d-hamiltonian` (gone — same `grep` check). No orphan references survive in the v4 body. The citation table is also pruned: Eq. 122 is mentioned in body once (the OBC remark), Sec. IV C once (converse-failure parenthetical), Fig. 3 once (deferral to numerics chapter); all three are kept in the table for verification purposes.
- **Connection to numerics chapter.** The deferred explicit Heisenberg evaluation goes into the numerics chapter (item G.21 of `thesis_numerics_plan.md`). The reference `(§ref-to-numerics)` in the body is a placeholder — the user should replace it with the correct `\ref{sec:...}` once the numerics-chapter section is in place. If the numerics section will instead reference *back* to (Eq. eq:H_heis) and the symbolic statement here, that's the cleaner direction; either way, the thesis has a single home for the explicit constant arithmetic.
- **Label-name harmonisation.** v4 silently switches v3's `eq:strang-prefactor` → `eq:strang-explicit-error` and `eq:trotter-count` → `eq:trott-num-bound` to match labels already typed in `1_preliminaries.tex:1062` and `1_preliminaries.tex:1071`. v4 also adopts the name `eq:H_heis` (previously `eq:heisenberg-hamiltonian` in v3) because `1_preliminaries.tex:1048` already references this label via `\eqref{eq:H_heis}`; that reference is currently dangling because the label is *not* yet defined anywhere in the existing `.tex`. Typing this draft in (specifically, the `<!-- \label{eq:H_heis} -->` comment on the $H_{\mathrm{Heis}}$ display) supplies the missing definition and resolves the existing forward-reference. So unlike the other two harmonisations, `eq:H_heis` is *introduced* by this draft, not inherited from the existing `.tex`.
- **v3 → v4 Revision notes.** The trim removed (i) the explicit even–odd bond-count derivation (`54 = 4 \cdot 108 \cdot (1/12+1/24)` arithmetic, ~30 lines), (ii) the disorder-contribution sub-paragraph (compressed to a single clause), (iii) the standalone TFIM and 2D-Heisenberg Hamiltonian displays plus their derivation paragraphs (compressed to a single sentence), and (iv) the explicit-constant $M$-count display $M \ge 15^{1/4}\sqrt{54}\,|J|^{3/2}n^{1/2}/(\sigma^{3/2}\varepsilon^{1/2})$ (now stated only in $\mathcal{O}$-form). The v3-review B1 algebraic slip is therefore *eliminated by deletion*: no explicit numerical constant survives in the preliminaries, so there is nothing for a reader to mis-propagate. The v3-review M1 *Principle* refactor (bond-level locality, two-layer count) is applied verbatim. The v3-review m1 unitarity-Hermiticity qualifier is also already in v3 and is preserved. v4 weighs in at roughly 60% of v3's length, with the Principle as the conceptual centerpiece.
- **What v4 deliberately keeps (load-bearing, do not re-edit):** the product-formula introduction, Suzuki recursion paragraph, palindromic / unitarity sentences, Childs commutator-scaling bound (Eq. eq:childs-comm-bound), Strang tight prefactor (Eq. eq:strang-explicit-error), Trotter-number bound (Eq. eq:trott-num-bound), and the qualitative $\sqrt{n}$-vs-$n^{3/2}$ contrast in the Consequence section. These all flow into Prop. prop:trotter-diss and Prop. prop:trotter-B in `2_methods.tex` and must remain stated in the symbolic form those propositions consume.
