# Trotterization

> **Insertion target:** `supplementary-informations/1_preliminaries.tex`, body of `\subsection{Trotterization.}` inside `\chapter{Quantum Circuits}` — between the figure `circ:trotter-strang` (ends at line 1034) and `\subsection{Linear Combination of Unitaries.}` at line 1037. The subsection label `sec:prelim-trotter` at line 968 is already referenced by `2_methods.tex:887`.
> **Figure referenced:** `circ:trotter-strang` (line 1033). The circuit is **not** modified; the draft only cites it.
> **Existing labels referenced:** `sec:prelim-trotter`, `circ:trotter-strang`, `eq:paulis`, `rem:palindromic`, `eq:childs-trotter-bound` (at `2_methods.tex:888`).
> **New labels introduced:** `eq:trotter-general`, `eq:strang`, `eq:childs-comm-bound`, `eq:strang-prefactor`, `eq:trotter-count`, `eq:heisenberg-hamiltonian`, `eq:heisenberg-evenodd-split`, `eq:heisenberg-alpha2`, `eq:trotter-M-count-heisenberg`, `eq:tfim-hamiltonian`, `eq:heisenberg-2d-hamiltonian`.
> **User action required in `2_methods.tex`:** add `\label{sec:hamiltonian-trotterization}` on the line immediately after `\subsection{Hamiltonian Trotterization.}` at `2_methods.tex:881`. Without this, the `\ref{sec:hamiltonian-trotterization}` invocations below will print `??`.
> **New citation keys introduced:** none. `childs2021theory` is already in `references.bib` and cited in `2_methods.tex`.

---

Up to this point the quantum circuits we discussed performed either read-out of an energy register (phase estimation) or state preparation on an ancilla. The dissipative generators of Chapter 5 require a third primitive: the unitary Hamiltonian evolution $\ee^{\pm \ii H t}$ itself, used inside the controlled cascade of (Fig. circ:qpe) — each controlled-$U^{2^k}$ is implemented via Trotter steps of $S_p(t_0)$. For the local Hamiltonians we care about — each term acting on a bounded number of qubits — there is no known shallow exact circuit implementing $\ee^{\pm \ii H t}$; we approximate it by a *product formula*, the subject of this subsection. The error analysis follows [Childs et al. 2021] [CITE: childs2021theory], whose commutator-scaling bound is the main analytical handle. The material here is a preliminary for the sharpened Trotter error bound of (§sec:hamiltonian-trotterization), where the generic commutator constant $\tilde\alpha^{(2)}_{\mathrm{comm}}$ is substituted by its explicit 1D-Heisenberg value.

## Product formulas

Write the Hamiltonian as a sum of $K$ non-commuting summands
$$
H \;=\; \sum_{k=1}^{K} H_k,
$$ <!-- \label{eq:trotter-general} -->

where each $H_k$ is either a single Pauli string (Eq. eq:paulis) or a commuting family of such strings whose exponential $\ee^{-\ii t H_k}$ is easy to implement (typically because all terms in $H_k$ share support disjointly or pairwise commute).\footnote{We follow the thesis convention $K$ for the number of summands; Childs et al.\ denote the same quantity $\Gamma$.} A *$p$-th order product formula* $S_p(t)$ is an ordered product of the single-summand evolutions $\ee^{-\ii t H_k}$ that approximates $\ee^{-\ii H t}$ up to error $\mathcal{O}(t^{p+1})$ for small $t$ [CITE: childs2021theory]. The two formulas that appear throughout the thesis are

$$
S_1(t) \;:=\; \ee^{-\ii t H_K}\cdots \ee^{-\ii t H_1}
\qquad\text{(Lie--Trotter, $p=1$)},
$$

and the symmetric Strang splitting
$$
S_2(t) \;:=\; \Bigl(\prod_{k=1}^{K} \ee^{-\ii (t/2) H_k}\Bigr)\Bigl(\prod_{k=K}^{1} \ee^{-\ii (t/2) H_k}\Bigr),
$$ <!-- \label{eq:strang} -->

with error $\|S_1(t) - \ee^{-\ii H t}\| = \mathcal{O}(t^2)$ and $\|S_2(t) - \ee^{-\ii H t}\| = \mathcal{O}(t^3)$ respectively. By construction $S_2(t)^\dagger = S_2(-t)$ — Strang is *palindromic* — a structural property that the Chapter 5 adjoint-symmetry argument (Rem. rem:palindromic) relies on to preserve KMS detailed balance under Trotterization. The circuit realisation of $S_2$ for the $n=3$ periodic Heisenberg chain is given in (Fig. circ:trotter-strang): three Pauli-bond groups ($XX+YY+ZZ$ per colour) arranged as A-B-C-B-A, with the central C-block carrying the odd-$n$ wrap-around bond at doubled angle.

Higher-order Suzuki formulas $S_{2k}$ can be built recursively from $S_2$ by the fractal construction $S_{2k}(t) = S_{2k-2}(u_k t)^2\,S_{2k-2}((1-4u_k)t)\,S_{2k-2}(u_k t)^2$ with $u_k = 1/(4 - 4^{1/(2k-1)})$, giving error $\mathcal{O}(t^{2k+1})$ [CITE: childs2021theory]. The Chapter 5 error analysis (Prop. prop:trotter-diss, Prop. prop:trotter-B) is stated for general $p$ and therefore carries the Suzuki recursion implicitly; the numerics of Part IV use $p=2$ throughout, so the displayed error analysis below is specialised to Strang. Since each $H_k$ is Hermitian by construction, each factor $\ee^{-\ii t H_k}$ is unitary, so $S_p(t)$ and its iterate $S_p^{(M)}(t) := S_p(t/M)^M$ are unitary for every choice of $p$ and $M$ — a property that the CKG error propagation (§sec:hamiltonian-trotterization) relies on when applying Parseval's identity to the Trotterised OFT.

## Commutator-scaling bound

A naive Taylor truncation only controls the error through $\|H_k\|$, ignoring commutativity of the summands. The key result of [Childs et al. 2021] [CITE: childs2021theory] is that the error of any $p$-th order formula is in fact controlled by *nested commutators* of the summands. Specialising Theorem 6 of [CITE: childs2021theory] to anti-Hermitian $\ii H_k$ and real $t$, one has
$$
\bigl\|S_p(t) - \ee^{-\ii H t}\bigr\|
\;=\; \mathcal{O}\!\bigl(\tilde\alpha_{\mathrm{comm}}\,|t|^{p+1}\bigr),
\qquad
\tilde\alpha_{\mathrm{comm}}
\;:=\; \sum_{k_1,\dots,k_{p+1}=1}^{K}
\bigl\|\,[H_{k_{p+1}},\,\cdots\,[H_{k_2},H_{k_1}]\cdots]\,\bigr\|,
$$ <!-- \label{eq:childs-comm-bound} -->

where the asymptotic constant is universal (independent of $K$, $\|H_k\|$, and the ordering fixed by the particular formula). The triangle-inequality bootstrap of [CITE: childs2021theory, Prop. 9 (Eq. 120)] gives the analogous tight first-order prefactor, which the Strang analysis of [CITE: childs2021theory, Prop. 10 (Eq. 121)] bootstraps to the tight second-order form we use below. Following the cleaner notation already adopted in `2_methods.tex:1032`, the latter reads
$$
\bigl\|S_2(t) - \ee^{-\ii H t}\bigr\|
\;\leq\; \sum_{k_1 < k_2}\!\left(
\frac{t^3}{12}\bigl\|[H_{k_2},[H_{k_2},H_{k_1}]]\bigr\|
\;+\; \frac{t^3}{24}\bigl\|[H_{k_1},[H_{k_1},H_{k_2}]]\bigr\|
\right)
\;+\; \mathcal{O}(t^4).
$$ <!-- \label{eq:strang-prefactor} -->

The coefficients $1/12$ and $1/24$ match the lowest-order BCH expansion and are tight [CITE: childs2021theory]. Splitting the target evolution $\ee^{-\ii H t}$ into $r$ Trotter steps of size $\delta t = t/r$ and applying the triangle inequality over steps yields the *Trotter-number bound* of [CITE: childs2021theory, Cor. 7]: to achieve total error $\varepsilon$, it suffices to take
$$
r \;\geq\; \Omega\!\left(\frac{\tilde\alpha_{\mathrm{comm}}^{1/p}\;t^{\,1+1/p}}{\varepsilon^{\,1/p}}\right).
$$ <!-- \label{eq:trotter-count} -->

For $p=1$ this is the $\mathcal{O}(\tilde\alpha_{\mathrm{comm}}\, t^2/\varepsilon)$ worst-case count, and for Strang ($p=2$) the $\mathcal{O}(\sqrt{\tilde\alpha_{\mathrm{comm}}}\, t^{3/2}/\sqrt{\varepsilon})$ bound that governs the Chapter 5 analysis.

## Specialisations

The Hamiltonian used in the numerical chapter and sketched in (Fig. circ:trotter-strang) is the isotropic periodic Heisenberg chain on $n$ qubits,
$$
H_{\mathrm{Heis}}
\;=\; J \sum_{i=1}^{n} \bigl( X_i X_{i+1} + Y_i Y_{i+1} + Z_i Z_{i+1} \bigr)
\;+\; \sum_{i=1}^{n} h_i Z_i,
\qquad X_{n+1} \equiv X_1,\ \ldots,
$$ <!-- \label{eq:heisenberg-hamiltonian} -->

with coupling $J>0$ and disorder fields $h_i \in [-h,h]$. Childs et al.\ analyse an OBC variant of this model with $j=1,\ldots,n-1$ [CITE: childs2021theory, Eq.~(122)]; the PBC extension below adds only a constant number of boundary-bond pairs and so preserves the $n$-scaling.

**Principle: geometrically-local commutator scaling.** All Hamiltonians of interest in this thesis are *geometrically local*: $H = \sum_b h_b$ where each *local term* $h_b$ is supported on a bounded-size subset of sites, is uniformly norm-bounded $\|h_b\| \le h_*$, and the interaction graph on the local terms has bounded maximum degree $d$ (each $h_b$ shares support with $\mathcal{O}(d)$ others). Such an $H$ admits a partition into $K = \mathcal{O}(1)$ commuting classes $H_k := \sum_{b\in\mathrm{class}_k} h_b$, each class a sum of pairwise-commuting local terms. Expanding the $(p+1)$-nested commutator $[H_{k_{p+1}}, \cdots [H_{k_2}, H_{k_1}]]$ term-by-term over local terms, only *connected* $(p+1)$-tuples $(b_1,\dots,b_{p+1})$ survive — disjoint-support terms commute — and on a graph of max degree $d$ the number of such tuples anchored at a pivot local term $b_1$ is at most $(2d)^{p}$. The $\Theta(n)$ pivot choices then give
$$
\tilde\alpha^{(p)}_{\mathrm{comm}} \;\leq\; c_{d,p}\,n\,h_*^{p+1} \;=\; \mathcal{O}(n),
$$
with $c_{d,p}$ a combinatorial constant depending only on $d$ and $p$. All the specialisations below are instances of this two-layer counting: outer commuting-class grouping (fixes $K$), inner bounded-degree local-term counting (gives the $\Theta(n)$ pivots and $\mathcal{O}(1)$ connected tuples per pivot). The converse failure is also instructive: an all-to-all interaction (e.g. the all-to-all TFIM analysed in [CITE: childs2021theory, Sec. IV C]) has unbounded $d$ and consequently the sharper $\mathcal{O}(n^3)$ commutator constant rather than $\mathcal{O}(n)$.

**1D periodic Heisenberg — two-colourable bond split.** For even $n$ the bond set $\{(i,i+1)\}_{i=1}^{n}$ is 2-colourable (as an edge colouring, not a site bipartition): write $E := \{(1,2),(3,4),\dots,(n{-}1,n)\}$ (*odd* bonds in the chain indexing) and $O := \{(2,3),(4,5),\dots,(n,1)\}$ (*even* bonds, closing the ring), and group the isotropic part as
$$
H^{\mathrm{iso}}_{\mathrm{Heis}} \;=\; H_E + H_O,
\qquad
H_E := J\!\!\sum_{(i,i+1)\in E}\!\! \vec\sigma_i\!\cdot\!\vec\sigma_{i+1},
\qquad
H_O := J\!\!\sum_{(i,i+1)\in O}\!\! \vec\sigma_i\!\cdot\!\vec\sigma_{i+1},
$$ <!-- \label{eq:heisenberg-evenodd-split} -->

with $\vec\sigma_i\!\cdot\!\vec\sigma_j := X_iX_j + Y_iY_j + Z_iZ_j$. Bonds within $E$ share no sites and therefore commute, as do bonds within $O$; only bonds *across* the two colours share one site. For odd $n$ the ring is not bipartite at the edge level and one needs a third group to close the boundary — the $R^2_{\sigma_3\sigma_1}$ gate spanning three wires in (Fig. circ:trotter-strang) is exactly this wrap-around bond, promoted to its own layer. The three-group Strang bound (Eq. eq:strang-prefactor) adds $\mathcal{O}(1)$ new inter-group triples beyond the two-group case; each still reduces to $\mathcal{O}(n)$ adjacent-bond pairs, so the $n$-scaling is unchanged and only the constant prefactor grows.

Applying the counting principle with $K=2$, $d=2$ to $H = H_E + H_O$: only the cross-pair $k_1=1, k_2=2$ contributes to the sum in (Eq. eq:strang-prefactor), leaving
$$
\tilde\alpha^{(2)}_{\mathrm{comm}}(H^{\mathrm{iso}}_{\mathrm{Heis}})
\;=\; \tfrac{1}{12}\bigl\|[H_O,[H_O,H_E]]\bigr\|
\;+\; \tfrac{1}{24}\bigl\|[H_E,[H_E,H_O]]\bigr\|.
$$

Expand each group as a sum over its bonds, $H_E = \sum_{b\in E} h_b$ with $h_b := J\,\vec\sigma_b$, and analogously for $H_O$; $\|h_b\| \le 3|J|$ by the triangle inequality on the three Pauli-Pauli terms. A doubly nested commutator $[h_{b_3},[h_{b_2},h_{b_1}]]$ vanishes unless $b_1, b_2, b_3$ form a *connected* chain of adjacent bonds (each sharing a site with the next), because disjoint bonds commute. The number of such connected triples on a periodic chain is $\Theta(n)$: fix the pivot bond $b_1 \in E$ in $n/2$ ways, pick $b_2 \in O$ adjacent to $b_1$ (two choices), then $b_3 \in E$ adjacent to $b_2$ (two choices). Each nested commutator satisfies
$$
\bigl\|[h_{b_3},[h_{b_2},h_{b_1}]]\bigr\| \;\leq\; 4\,\|h_{b_3}\|\,\|h_{b_2}\|\,\|h_{b_1}\| \;\leq\; 4\cdot (3|J|)^3 \;=\; 108\,|J|^3,
$$

and the same bound holds with $E,O$ swapped. Collecting:
$$
\tilde\alpha^{(2)}_{\mathrm{comm}}(H^{\mathrm{iso}}_{\mathrm{Heis}})
\;\leq\; 54\,|J|^3\,n,
\qquad
\bigl\|S_2(\delta t) - \ee^{-\ii H_{\mathrm{Heis}}^{\mathrm{iso}}\delta t}\bigr\|
\;=\; \mathcal{O}\bigl(|J|^3\,n\,|\delta t|^3\bigr),
$$ <!-- \label{eq:heisenberg-alpha2} -->

with the explicit constant $54 = 4 \cdot 108 \cdot (\tfrac{1}{12}+\tfrac{1}{24})$ (the factor $4$ absorbs the adjacent-bond count at each pivot and is loose). This is a direct specialisation of [CITE: childs2021theory, Prop. 10 (Eq. 121)] to the even–odd partition; Childs et al.\ do not tabulate the $p=2$ Heisenberg prefactor themselves, but they plot the $p=4$ version in [CITE: childs2021theory, Fig. 3], where the corresponding Suzuki-$\mathcal{S}_4$ bound is loose by only a factor $\approx 5$ at $n=10$ for the even–odd ordering.

**Disorder contribution.** The disorder field $\sum_i h_i Z_i$ commutes with each $Z_j Z_{j+1}$ but not with $X_j X_{j+1}$ or $Y_j Y_{j+1}$. Re-running the same connected-triple count with one or two factors replaced by $h_i Z_i$ (norm $\leq h$) and the rest by Heisenberg bonds (norm $\leq 3|J|$), the non-vanishing triples are again $\Theta(n)$ and each is bounded by $4\cdot(3|J|)^2\cdot h$ or $4\cdot 3|J|\cdot h^2$. The disorder therefore contributes
$$
\tilde\alpha^{(2)}_{\mathrm{comm}}\bigl[\text{disorder}\bigr] \;=\; \mathcal{O}\bigl(n(h^2 J + h J^2)\bigr),
$$

which is of the same $n$-scaling as the isotropic contribution and can be absorbed into the prefactor of (Eq. eq:heisenberg-alpha2) at the price of replacing $|J|^3$ by a polynomial of total degree three in $(|J|,h)$. For the numerics of Part IV we keep $h \lesssim J$, so this is a bounded correction.

**Other geometrically-local models.** The same principle specialises cleanly to two other canonical models that will reappear later in the thesis. The *1D periodic transverse-field Ising model* (TFIM) on $n$ qubits is
$$
H_{\mathrm{TFIM}} \;=\; J\sum_{i=1}^{n} Z_i Z_{i+1} \;+\; h\sum_{i=1}^{n} X_i,
\qquad Z_{n+1} \equiv Z_1,
$$ <!-- \label{eq:tfim-hamiltonian} -->

with coupling $J$ and transverse field $h$; with the two-group split $H_{ZZ} + H_X$, the connected-triple count gives $\tilde\alpha^{(2)}_{\mathrm{comm}}(H_{\mathrm{TFIM}}) = \mathcal{O}(n(J^2 h + J h^2))$. The *2D periodic Heisenberg model* on the $L\times L$ square lattice ($n = L^2$ qubits) is
$$
H^{2\mathrm{D}}_{\mathrm{Heis}}
\;=\; J \sum_{\langle i,j\rangle}\bigl( X_i X_j + Y_i Y_j + Z_i Z_j \bigr)
\;+\; \sum_{i} h_i Z_i,
$$ <!-- \label{eq:heisenberg-2d-hamiltonian} -->

where $\langle i,j\rangle$ runs over nearest-neighbour bonds on the square lattice with periodic boundaries in both directions. The bonds split into four commuting colour classes (horizontal-even, horizontal-odd, vertical-even, vertical-odd) with bond-level interaction degree $d = 6$, yielding $\tilde\alpha^{(2)}_{\mathrm{comm}}(H^{2\mathrm{D}}_{\mathrm{Heis}}) = \mathcal{O}(J^3 n)$ with a larger constant prefactor than in 1D. Both are routine applications of the counting argument above.

## Consequence for the CKG Lindbladian

The Hamiltonian-Trotterization analysis of (§sec:hamiltonian-trotterization) states its main propositions (Prop. prop:trotter-diss, Prop. prop:trotter-B) in terms of the *generic* constant $\tilde\alpha^{(2)}_{\mathrm{comm}}$ of (Eq. eq:strang-prefactor), making no structural use of $H$. For the 1D periodic Heisenberg chain of (Eq. eq:heisenberg-hamiltonian) used in the numerics of Part IV, substituting the specialised constant (Eq. eq:heisenberg-alpha2) into the Strang-splitting dissipator bound
$$
\bigl\|\bar{\mathcal{L}}_{\mathrm{diss}} - \tilde{\mathcal{L}}_{\mathrm{diss}}\bigr\|_{1\to 1}
\;\leq\; \frac{\sqrt{15}\,\tilde\alpha^{(2)}_{\mathrm{comm}}}{M^2\,\sigma^3}
$$

of Prop. prop:trotter-diss sharpens the $n$-dependence from the generic $\mathcal{O}(n^3)$ of a blind triangle inequality on (Eq. eq:childs-comm-bound) — bounding each nested commutator by $4\|H_k\|^3$ with $\|H_E\|,\|H_O\| = \mathcal{O}(|J|n)$ and summing $K^2 = 4$ ordering terms — to the $\mathcal{O}(n)$ scaling of (Eq. eq:heisenberg-alpha2). Solving for the number of Strang steps needed to reach dissipator error $\varepsilon$ and substituting $\tilde\alpha^{(2)}_{\mathrm{comm}} \le 54\,|J|^3\,n$ from (Eq. eq:heisenberg-alpha2) gives
$$
M \;\ge\; \sqrt{\frac{\sqrt{15}\cdot 54\,|J|^3\,n}{\sigma^3\,\varepsilon}}
\;=\; 15^{1/4}\sqrt{54}\;\frac{|J|^{3/2}\,n^{1/2}}{\sigma^{3/2}\,\varepsilon^{1/2}}
\;=\; \mathcal{O}\!\left(\frac{|J|^{3/2}\,n^{1/2}}{\sigma^{3/2}\,\varepsilon^{1/2}}\right).
$$ <!-- \label{eq:trotter-M-count-heisenberg} -->

The payoff of the specialisation is visible in the $n$-exponent: the required Trotter-step count scales as $\sqrt{n}$ for the 1D Heisenberg chain, not the $n^{3/2}$ that the generic $\mathcal{O}(n^3)$ commutator constant would yield. The same substitution propagates through Prop. prop:trotter-B for the coherent term, with no structural change to either proof. We carry out the substitution and its numerical cross-check in Part IV; the rest of the thesis keeps $\tilde\alpha^{(2)}_{\mathrm{comm}}$ generic so the bounds remain applicable to any decomposition $H = \sum_k H_k$.

---

## Citations

- **`childs2021theory`** — Childs, A.~M., Su, Y., Tran, M.~C., Wiebe, N., Zhu, S., *"Theory of Trotter Error with Commutator Scaling"*, **Phys. Rev. X 11, 011020 (2021)**, arXiv:1912.08854. Already in `references.bib` and cited throughout `2_methods.tex`.

No new bibtex entries introduced.

**Page-level verification of every Childs citation used above** (against `supplementary-informations/Childs et al. - 2021 - Theory of Trotter Error with Commutator Scaling.pdf`):

| Citation in draft | PDF location | Verified |
|---|---|---|
| Thm 6 (master commutator bound) | p. 12, Eq. (44) — anti-Hermitian form | ✓ |
| Cor. 7 (Trotter-number count) | p. 13, Eq. (45) | ✓ |
| Prop. 9 (first-order tight bound) | p. 22, Eq. (120) | ✓ |
| Prop. 10 (Strang tight bound with 1/12, 1/24) | p. 22, Eq. (121) | ✓ |
| Eq. (122) (OBC Heisenberg + random field) | p. 23 | ✓ (flagged as OBC vs our PBC) |
| Fig. 3 (factor-$\approx 5$ at $n=10$, even–odd, $p=4$) | p. 24, text paragraph after figure | ✓ (explicitly $\mathcal{S}_4$, not Strang) |
| Sec. IV C (all-to-all TFIM, $\mathcal{O}(n^3)$ commutator) | p. 19, Eq. (82) and Eq. (93) | ✓ (cited only as converse to the geometric-locality principle) |

---

## Revision notes (v3 vs v2 — delete before merging)

- **Item 1 — Sign-convention paragraph cut.** The standalone paragraph beginning "The sign convention $S_p(t) \approx \ee^{-\ii Ht}$ used here differs from…" is removed entirely; the draft silently uses $\ee^{-\ii Ht}$.
- **Item 2 — "Gain over blind triangle" meta-paragraph cut.** Removed; the same quantitative point is already made in the Consequence section.
- **Item 3 — 2D Heisenberg paragraph collapsed.** Standalone "**2D Heisenberg on a square lattice.**" sub-paragraph removed; its content is folded into the one-paragraph "Other geometrically-local models" at the end of the Specialisations section.
- **Item 4 — Lie–Trotter tight-prefactor display cut.** Eq. `eq:lie-trotter-prefactor` removed; replaced by a single inline sentence acknowledging Prop. 9 as the $p=1$ analogue that bootstraps Prop. 10.
- **Item 5 — Suzuki recursion retained.** Unchanged from v2.
- **Item 6 — Geometrically-local principle added.** New "**Principle: geometrically-local commutator scaling.**" remark states the bond-level two-layer counting argument (outer commuting-class grouping $K = \mathcal{O}(1)$, inner bounded-degree local-term graph with max degree $d$; $\Theta(n)$ pivots times $(2d)^p$ connected tuples per pivot), giving $\tilde\alpha^{(p)}_{\mathrm{comm}} \le c_{d,p}\,n\,h_*^{p+1}$. All-to-all TFIM is called out as the converse-failure case. 1D Heisenberg derivation is presented as the canonical instance, with the explicit $54|J|^3 n$ bound preserved; 1D TFIM and 2D Heisenberg are collapsed into a single closing paragraph.
- **Item 7 — Explicit $M$-count for CKG Heisenberg.** New display (Eq. eq:trotter-M-count-heisenberg) in the Consequence section: $M \ge 15^{1/4}\sqrt{54}\,|J|^{3/2}\,n^{1/2}/(\sigma^{3/2}\varepsilon^{1/2})$, showing $\sqrt{n}$ rather than the generic $n^{3/2}$. The display is written in the $M \ge \sqrt{\sqrt{15}\cdot 54\,|J|^3 n/(\sigma^3 \varepsilon)}$ form to make the square-root structure plainly visible.
- **Item 8 — Unitarity of $S_p^{(M)}$.** One sentence added after the Suzuki paragraph: since each $H_k$ is Hermitian, each factor $\ee^{-\ii t H_k}$ is unitary, hence $S_p(t)$ and $S_p^{(M)}(t)$ are unitary, which is what the CKG Parseval argument relies on.
