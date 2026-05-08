## Insertion target

`drafts/ckg-vs-dll-comparison-findings.md`, append as a new top-level
section §9 *after* §8.7 (which already gestures at the multi-rank
question). The new §9 also retroactively contradicts the *predicted
null-result* sentence in §8.7 — that paragraph should get a one-line
addendum noting that §9 supersedes it.

# 9. Empirical $\tau_{\text{mix}}$ vs channel count $k$ (qf-7go)

§8.7 predicted that adding rank to DLL — a parametric multi-channel
filter with $k > 1$ symmetrised translates of the same base shape — would
*not* shorten $\tau_{\text{mix}}$, because the per-frequency Metropolis
shape on the diagonal of $\alpha$ already determines the measured rate
(H2). This section reports the qf-7go.6 numerical sweep that tests that
prediction directly.

## 9.1 Construction

For each base filter $f̂_{\text{base}}$ (DLL Metropolis-type Eq. 3.20 or
DLL Gaussian-type Eq. 3.22) and integer $k \geq 1$, build $k$
**symmetrised-translate channels** indexed by $\nu_\ell$:

$$
q_\ell(\nu) \;=\; \sqrt{\tfrac{w_\ell}{2}} \,
\bigl[\, q_{\text{base}}(\nu - \nu_\ell) + q_{\text{base}}(\nu + \nu_\ell) \,\bigr]
\quad (\nu_\ell \neq 0),
\qquad
q_\ell(\nu) \;=\; \sqrt{w_\ell} \, q_{\text{base}}(\nu)
\quad (\nu_\ell = 0).
$$ <!-- \label{eq:dll-symmetrised-translate} -->

The associated Lindblad operator and Kossakowski matrix per coupling
$A^a$ are

$$
L_a^{(\ell)} \;=\; \sum_{\nu \in B_H} q_\ell(\nu) e^{-\beta\nu/4} A^a_{\nu},
\qquad
\alpha^{(\ell)}(\nu, \nu') \;=\; q_\ell(\nu) q_\ell(\nu')
\, e^{-\beta(\nu+\nu')/4},
$$ <!-- \label{eq:dll-multirank-perchannel} -->

and the multi-channel dissipator and coherent term are the linear sum
over $\ell$:

$$
\mathcal{L}^{\text{multi}}(X) \;=\; i [G^{\text{multi}}, X]
\;+\; \sum_a \sum_{\ell=1}^k \Bigl( (L_a^{(\ell)})^\dagger X L_a^{(\ell)}
- \tfrac{1}{2} \{ (L_a^{(\ell)})^\dagger L_a^{(\ell)}, X \} \Bigr),
\qquad
G^{\text{multi}} \;=\; \sum_\ell G^{(\ell)},
$$ <!-- \label{eq:dll-multirank-lindbladian} -->

where each per-channel $G^{(\ell)}$ is the standard DLL coherent operator
(Eq. eq:dll-coherent / Ding–Li–Lin Eq. 3.7) built from the channel's
$f̂_\ell$. The decomposition $G^{\text{multi}} = \sum_\ell G^{(\ell)}$
follows because the canonical-form map $V \mapsto -i\,\tanh \circ
\log(\Delta^{1/4}_{\sigma_\beta}) \cdot V$ from the joint $V = \tfrac{1}{2}
\sum_{a,\ell} (L_a^{(\ell)})^\dagger L_a^{(\ell)}$ is *linear*
([Ding et al. 2024][CITE: ding2025efficient], Theorem 10 / Eq. 2.33).

By construction each $q_\ell$ is real-even, so $q_\ell(-\nu) = q_\ell(\nu)
= \overline{q_\ell(\nu)}$ — the Eq. 3.2 symmetry condition holds per
channel, and the per-channel Kossakowski $\alpha^{(\ell)}$ satisfies the
KMS-DBC skew-symmetry $\alpha(\nu, \nu') = \alpha(-\nu', -\nu)
e^{-\beta(\nu+\nu')/2}$ (Eq. 4.7 of [Ding et al. 2024]). The sum
$\alpha^{\text{multi}} = \sum_\ell \alpha^{(\ell)}$ inherits both
properties, and the full $\mathcal{L}^{\text{multi}}$ is verified to be
KMS-DB by the discriminant test
$\|\mathcal{D} - \mathcal{D}^\dagger\|_{2 \to 2} / \|\mathcal{D}\|_{2 \to 2}
\leq 10^{-10}$ on the $n = 3$ disordered Heisenberg fixture across
$\beta \in \{1, 5, 10\}$, $k \in \{1, 2, 4\}$.

The construction adds rank to the Kossakowski (each $\alpha^{(\ell)}$
is rank-$1$, so $\text{rank}(\alpha^{\text{multi}}) \leq k$), recovering
the spirit of the §3 framing without invoking the SVD-of-CKG route the
DLL paper sketches in Eq. 4.10. Centers $\{\nu_\ell\}$ for the sweep
below are uniform on $[0, S/2]$ with $S = 2$ (Metropolis flat-top
constraint); the same center grid is used for the Gaussian base.

## 9.2 Sweep grid

| Parameter | Values |
|-----------|--------|
| $n$ | $3$ |
| $\beta$ | $\{1, 5, 10, 20\}$ |
| $k$ | $\{1, 2, 4, 8\}$ |
| Base filter | DLL Metropolis ($S = 2$); DLL Gaussian |
| Centers | uniform on $[0, S/2]$ |
| $\varepsilon_{\text{target}}$ | $10^{-3}$ |
| Domain | BohrDomain (matrix-free) |

Each cell uses the qf-lkb adaptive horizon
$t_{\max} = T \cdot \text{gap}^{-1}$ with $T = \max(5, 1.5 \log_{10}
\varepsilon^{-1}) = 5$, the qf-lkb.10 bi-exponential extrapolation, and
the qf-lkb.10 observed-mixing fallback when extrapolation fails. The
spectral-gap pre-pass uses matrix-free Krylov on
`apply_lindbladian!` (qf-lkb.9 path) so the cost stays linear in
$|B_H| \cdot d^2$; total wall time across the 32 cells was 17 s on a
single Julia thread.

## 9.3 Results

The sweep contradicts the §8.7 prediction: $\tau_{\text{mix}}$ shortens
substantially with $k$ on every $(\beta, \text{shape})$ cell.

**DLL Metropolis ($n = 3$):**

| $\beta$ | $k=1$ | $k=2$ | $k=4$ | $k=8$ | $k{=}1 \to k{=}8$ ratio |
|---:|---:|---:|---:|---:|---:|
| $1$  | $11.2$  | $4.56$  | $1.80$  | $0.82$  | $13.7\times$ |
| $5$  | $18.1$  | $13.6$  | $5.06$  | $2.22$  | $8.2\times$  |
| $10$ | $23.3$  | $22.2$  | $10.3$  | $4.52$  | $5.2\times$  |
| $20$ | $35.4$  | $35.3$  | $26.4$  | $12.9$  | $2.7\times$  |

**DLL Gaussian ($n = 3$):**

| $\beta$ | $k=1$ | $k=2$ | $k=4$ | $k=8$ | $k{=}1 \to k{=}8$ ratio |
|---:|---:|---:|---:|---:|---:|
| $1$  | $6.79$  | $2.65$  | $1.08$  | $0.50$  | $13.7\times$ |
| $5$  | $11.8$  | $11.4$  | $4.70$  | $2.02$  | $5.9\times$  |
| $10$ | $25.3$  | $25.3$  | $9.69$  | $4.36$  | $5.8\times$  |
| $20$ | $\mathbf{264}$ | $\mathbf{264}$ | $\mathbf{20.8}$ | $\mathbf{13.3}$ | $\mathbf{19.9\times}$ |

The full picture in `drafts/figures/numerics/dll_multirank_taumix.{png,pdf}`
plots $\tau_{\text{mix}}(k)$ on a log-$y$ scale, one panel per base shape,
with one line per $\beta$.

## 9.4 Two headline observations

**(i) The $\beta = 20$ DLL Gaussian collapse is fully recovered by adding
rank.** The single-channel DLL Gaussian $\tau_{\text{mix}}$ at
$(n = 3, \beta = 20)$ is $\approx 264$ — an order of magnitude worse than
the corresponding CKG smooth-Metropolis baseline of $\approx 28$
(qf-mto fair-comparison data) — because the Eq. 3.22 weight
$\hat{f}_{\text{base}}(\nu) \propto e^{-(\beta\nu+1)^2/8}$ has width
$O(\beta^{-1})$, leaving a single near-dephasing channel at low
temperature. With $k = 4$ channels at centers $\{0, \tfrac{1}{3},
\tfrac{2}{3}, 1\}$, $\tau_{\text{mix}}$ drops to $20.8$; at $k = 8$, to
$13.3$. The Kossakowski rank, not the per-channel filter shape, is the
bottleneck *in this regime*. The collapse from §3 of the comparison
draft is therefore an artefact of insisting on rank $1$, not a structural
limitation of the Gaussian-type weighting.

**(ii) DLL Metropolis at $k = 8$ is comparable to CKG even before any
shape-tuning.** Across $\beta \in \{5, 10, 20\}$ the multi-rank
Metropolis $\tau_{\text{mix}}$ at $k = 8$ stays in the same ballpark as
or slightly below the corresponding CKG smooth-Metropolis baselines
(qf-mto data at $n = 3$: $\beta = 5$: CKG $\sim 12$ vs DLL $2.22$;
$\beta = 10$: CKG $\sim 18$ vs DLL $4.52$; $\beta = 20$: CKG $\sim 28$ vs
DLL $12.9$). This is *not* a like-for-like comparison — CKG sM uses
$\sim 24$ effective channels per coupling (full-rank Kossakowski), while
the multi-rank DLL here uses $k = 8$ — but it shows that the rank-1
Metropolis "one wide brush" of §4 is leaving substantial relaxation
budget on the table and that the absolute floor of attainable
$\tau_{\text{mix}}$ within the DLL framework is well below the
single-channel value.

## 9.5 What revises §8.7's prediction

§8.4 / H2 argued that the *diagonal* $\alpha(\nu, \nu)$ — the
per-frequency dissipation rate — drives $\tau_{\text{mix}}$. That claim
is consistent with §9: the multi-channel construction *does* enrich the
diagonal, because $\alpha^{\text{multi}}(\nu, \nu) = \sum_\ell |q_\ell(\nu)|^2
e^{-\beta\nu/2}$ aggregates over all channels' weight at $\nu$. With
shifted $q_\ell$, more $\nu$ values inside the relevant Bohr-frequency
range receive non-negligible diagonal mass, while the antidiagonal
also gains contributions $|q_\ell(\nu)|^2 e^{-\beta\nu/2}$ from
channels whose centers $\nu_\ell$ straddle $\nu$. The §8.7 prediction
implicitly assumed that *adding channels does not change the diagonal
shape* — true only if the centers are all $0$ (in which case $k$
identical copies just rescale by $\sqrt{k}$ and produce the same
Lindbladian up to a $k$-fold trivial rescaling). The qf-7go.5
construction breaks that degeneracy by spreading channels across the
spectrum.

A cleaner re-statement of H2 in light of §9 is: *the diagonal-of-$\alpha$
shape that minimises $\tau_{\text{mix}}$ is not unique to CKG — it is
recoverable inside the DLL framework as a finite combination of
ν-shifted base filters, and matching it tightens both DLL Metropolis
and DLL Gaussian to within a small factor of the CKG sM benchmark.* The
follow-up question, deferred to a future epic, is whether the qf-mto
Dirichlet-form analysis can predict the optimal centers $\{\nu_\ell\}$
analytically rather than by sweep.

## 9.6 Implication for the cost analysis of §5

§5 / Remark 23 argued that the rank-1 DLL Metropolis pays its price not
in $\tau_{\text{mix}}$ but in coherent-block-encoding cost — a
$\Theta(\beta)$ blow-up of $\|G\|$ relative to CKG sM. Multi-rank DLL
inherits the same cost picture *per channel*: each $G^{(\ell)}$ has the
same Remark 23 scaling, so $\|G^{\text{multi}}\| \leq \sum_\ell
\|G^{(\ell)}\| = O(k \cdot \beta^{1.4})$ in the worst case. The
$\tau_{\text{mix}}$ improvements above are therefore offset by a factor
$k$ in the per-step block-encoding cost: at $\beta = 20$, the
$k = 1 \to k = 8$ DLL Gaussian gain of $\sim 20\times$ in
$\tau_{\text{mix}}$ comes at $\sim 8\times$ in coherent cost, for a net
$\sim 2.5\times$ total-cost reduction. This is still a real win at the
*total* level, but more modest than the $\tau_{\text{mix}}$ table alone
suggests, and it does not flip the §6 verdict that CKG remains the
better choice for low-$T$ large-system Gibbs sampling under the cost
metric of [Ding et al. 2024][CITE: ding2025efficient], Table 1.
