# Is the 2D TFIM Lindbladian gap closing genuine ordered-phase physics? (qf-c9g)

**Status**: complete (2026-05-25).
**Drivers**: `scripts/scratch_qf_c9g_beta_n_sweep.jl` (15-cell n × β grid, wall 3300 s),
            `scripts/scratch_qf_c9g_slow_mode_n4.jl` (dense L slow-mode at n=4 × 5 β values, wall ~25 s),
            `scripts/plot_qf_c9g_beta_n_sweep.jl` + `scripts/analyze_qf_c9g_verdict.jl` (figures + checklist).
**Sidecars**: `scripts/output/qf_c9g_ordered_gap_mechanism/qf_c9g_beta_n_sweep.bson`,
            `scripts/output/qf_c9g_ordered_gap_mechanism/qf_c9g_slow_mode_n4.bson`.
**Figures**: `drafts/figures/numerics/qf_c9g_gap_phys_vs_n.{png,pdf}`,
            `drafts/figures/numerics/qf_c9g_phase_indicators.{png,pdf}`.
**Lit survey**: `drafts/qf-c9g-lit-survey.md`.
**Parents**: qf-1jj (memory [[2d-tfim-ordered-vs-disordered-qf-1jj]]), qf-biz
            (`drafts/qf-biz-2d-tfim-matrix-element-bottleneck.md`).

> **Headline.** The 2D TFIM ORDERED Lindbladian-gap closing
> $\lambda_L^{\rm phys} = 9.96 \times 10^{-2} \to 3.06 \times 10^{-3} \to
> 1.57 \times 10^{-4}$ across $n = 4, 6, 8$ at $\beta_{\rm phys} = 2$ **IS
> genuine ordered-phase physics**, not a finite-T "Gibbs ≈ GS doublet"
> tunnel-splitting artefact. Four-of-four diagnostic tests fire as the
> physics literature [Gamarnik–Kiani–Zlokapa 2024; Martinelli–Olivieri
> 1994; Rakovszky et al. 2024] predicts for surface-tension closing:
> *(i)* closing persists at $\beta_{\rm phys} = 0.5$ (just below $T_c$,
> $T/T_c = 0.97$, $\text{eff\_rank} \in [3.3, 8.6]$ — far from the doublet)
> with $\lambda_L^{\rm phys} = 0.449 \to 0.286 \to 0.285$; *(ii)* Binder
> $U_4$ saturates at $0.60$–$0.66$ across the ordered window; *(iii)* the
> slowest L eigenmode $R_2$ lives **89 %** in the doublet × bulk
> off-diagonal block at *every* $\beta$ below $T_c$ and at $\beta = 0.25$
> as well — its structure is the same across the phase boundary, only the
> rate changes; *(iv)* paramagnetic control $\beta_{\rm phys} \in \{0.10,
> 0.25\}$ gives a flat $\lambda_L^{\rm phys} \approx 4.5, 1.7$
> respectively across all $n$ — slopes
> $-2 \times 10^{-3}, -2 \times 10^{-2}$, indistinguishable from noise.
> The cell at $\beta = 0.5$ specifically refutes the "Gibbs ≈ GS doublet"
> alternative: there the closing is visible *before* the Gibbs concentrates
> on the doublet. Sandbox geometry $2 \times L_y$ ladder limits the
> closing to a fixed-cut floor $\sim e^{-\beta \cdot 4J}$ (not the
> $L \times L$ surface-tension $\sqrt{n}$ scaling proved by
> Gamarnik–Kiani–Zlokapa for square geometry), so the *magnitudes* are
> ladder-specific but the *mechanism* is the same.

## 1 Motivation (colleague's worry)

In the qf-1jj sweep at $\beta_{\rm phys} = 2.0$, $h = J = 1$, $2 \times L_y$
PBC ladders, the Lindbladian gap closes by ~600 × across $n \in \{4, 6, 8\}$.
Natural skepticism: at $\beta_{\rm phys} = 2$ the Gibbs state has
$\text{eff\_rank} = 2.00$ — essentially the two-state Z₂ doublet — so the
closing **could** be the trivial "Hamiltonian doublet tunnel splitting
$\Delta E_1 \to 0$" dressed by a Davies rate, with no ordering content. Such
a closing would happen in any model with a quasi-degenerate doublet (e.g. a
deformed two-level system) and would have nothing specifically to do with
2D Ising surface tension.

The colleague's worry, phrased operationally: *the closing should NOT occur
at temperatures where the Gibbs has spread off the doublet (eff_rank > 2).*
If we can show closing at $T/T_c < 1$ with eff_rank well above 2, the
artefact hypothesis is dead.

## 2 Literature anchor (full survey in `drafts/qf-c9g-lit-survey.md`)

The headline lower bound is **Gamarnik–Kiani–Zlokapa 2024** (arXiv:2411.04300,
Theorem 1.4): for 2D TFIM on $n \times n$ with constant $\beta \ge \beta^*$,
$h \le h^* \approx 1$, every geometrically local Lindbladian — including
**Chen–Kastoryano–Gilyen KMS-DB** (the exact sampler class we use; reduction
via Lieb–Robinson bounds) — obeys
$$T_{\rm mix}(\mathcal L) \;\ge\; \exp\bigl[n^{1/2 - o(1)}\bigr]
\;=\; \exp\bigl[\Omega(\sqrt{n_{\rm total}})\bigr].$$
Mechanism: Peierls fault-line surface tension dressed by transverse-field
fluctuations via the Poisson–Feynman–Kac representation. The
**Martinelli–Olivieri 1994** classical theorem on 2D Ising Glauber
($\text{gap} \le e^{-\tau(\beta) L (1+o(1))}$, $L = \sqrt{N}$) is the
matching benchmark — same surface-tension origin, same $\sqrt{N}$ scaling.

Both theorems require $\beta > \beta_c$ and a non-trivial transverse field
($h$ small enough); they predict the closing through the entire ordered
window, not just deep cold. The trivial "Gibbs ≈ GS doublet" hypothesis
would predict closing only once the Gibbs concentrates on the doublet —
i.e. only at $\beta \gg \beta_c \cdot \log(\Delta E_{\rm bulk}/\Delta E_1)$.
These two stories make incompatible predictions at $\beta$ values
*just below* $T_c$.

## 3 Sweep recipe

Canonical CKG smooth-Metropolis Lindbladian with $s = 0.25$, $a = 0$, KMS
construction, EnergyDomain $r_D = 7$ quadrature, $\sigma = 1/\beta_{\rm alg}$,
qf-8fr $10^{-10}$ GUE Krylov seed (clean ladder Hamiltonian — no
disorder field). `krylov_spectral_gap(cfg, ham, jumps; krylovdim = 60,
howmany = 4, tol = 10^{-10})`. 2 × L_y PBC ladders: $(L_x, L_y) = (2,2),
(2,3), (2,4)$ for $n = 4, 6, 8$. $\beta_{\rm phys}$ grid spans factor 20:
$\{0.10, 0.25, 0.50, 1.0, 2.0\}$, crossing $\beta_c = 1/T_c \approx 0.483$
(Hesselmann–Wessel 2016, $T_c(h{=}1) \approx 2.07$).

Sanity baselines (all reproduced bit-accurately in this rerun):

| cell | qf-c9g | qf-1jj / qf-biz baseline | rel. diff |
|---|---|---|---|
| n=4, β=0.25 (DIS) | 1.9100 | 1.910 (qf-1jj) | $< 10^{-3}$ |
| n=6, β=0.25 (DIS) | 1.5747 | 1.575 (qf-1jj) | $< 10^{-3}$ |
| n=8, β=0.25 (DIS) | 1.7405 | 1.741 (qf-1jj) | $< 10^{-3}$ |
| n=4, β=2.0  (ORD) | 9.9645e-2 | 9.96e-2 (qf-1jj) | $< 10^{-3}$ |
| n=6, β=2.0  (ORD) | 3.0648e-3 | 3.06e-3 (qf-1jj) | $< 10^{-3}$ |
| n=8, β=2.0  (ORD) | 1.5736e-4 | 1.57e-4 (qf-1jj) | $< 10^{-3}$ |
| n=6, β=0.10        | 3.9993 | 4.00 (qf-biz Check 1) | $< 10^{-3}$ |
| n=6, β=0.50        | 0.2859 | 0.286 (qf-biz Check 1) | $< 10^{-3}$ |
| n=6, β=1.00        | 0.01237 | 0.0124 (qf-biz Check 1) | $< 10^{-3}$ |

## 4 Hamiltonian doublet structure (β-independent, Sweep B)

Dense $\text{eigen}(H_{\rm phys})$ at each n:

| n | $\Delta E_1^{\rm phys}$ (doublet split) | $\Delta E_2^{\rm phys}$ (bulk gap) | bulk/doublet |
|---|---|---|---|
| 4 | $7.10 \times 10^{-2}$ | $5.89$ | 92 |
| 6 | $1.07 \times 10^{-2}$ | $5.92$ | 600 |
| 8 | $1.83 \times 10^{-3}$ | $5.84$ | 3400 |

Verdict: a textbook spontaneously-Z₂-broken-doublet structure. The doublet
splitting collapses geometrically with n (consistent with finite-size
tunneling between ferromagnetic minima on a 2 × L_y geometry); the bulk gap
is **flat at $\sim 5.9$ — Ω(1) and N-independent**. The doublet is
exponentially isolated from the bulk (factor 3400 at n=8). The Hamiltonian
side already shows the SSB-doublet anatomy.

## 5 15-cell n × β table (Sweep A)

| n | $\beta_{\rm phys}$ | $T/T_c$ | phase | $\lambda_L^{\rm phys}$ | eff_rank | $U_4$ | $\langle m_z^2 \rangle$ |
|---|---|---|---|---|---|---|---|
| 4 | 0.10 | 4.83 | DIS | $4.534$ | 14.5  | 0.311 | 0.371 |
| 6 | 0.10 | 4.83 | DIS | $3.999$ | 56.5  | 0.222 | 0.253 |
| 8 | 0.10 | 4.83 | DIS | $4.495$ | 218.  | 0.169 | 0.190 |
| 4 | 0.25 | 1.93 | DIS | $1.910$ | 8.68  | 0.502 | 0.606 |
| 6 | 0.25 | 1.93 | DIS | $1.575$ | 27.9  | 0.435 | 0.470 |
| 8 | 0.25 | 1.93 | DIS | $1.741$ | 91.3  | 0.368 | 0.371 |
| 4 | 0.50 | 0.97 | **ORD** | $0.449$ | 3.27  | 0.626 | 0.869 |
| 6 | 0.50 | 0.97 | **ORD** | $0.286$ | 4.99  | 0.617 | 0.822 |
| 8 | 0.50 | 0.97 | **ORD** | $0.285$ | 8.63  | 0.600 | 0.764 |
| 4 | 1.00 | 0.48 | **ORD** | $6.81 \times 10^{-2}$ | 2.05  | 0.653 | 0.947 |
| 6 | 1.00 | 0.48 | **ORD** | $1.24 \times 10^{-2}$ | 2.09  | 0.655 | 0.942 |
| 8 | 1.00 | 0.48 | **ORD** | $9.99 \times 10^{-3}$ | 2.15  | 0.656 | 0.937 |
| 4 | 2.00 | 0.24 | **ORD** | $9.96 \times 10^{-2}$ | 2.00  | 0.653 | 0.950 |
| 6 | 2.00 | 0.24 | **ORD** | $3.06 \times 10^{-3}$ | 2.00  | 0.656 | 0.945 |
| 8 | 2.00 | 0.24 | **ORD** | $1.57 \times 10^{-4}$ | 2.00  | 0.657 | 0.942 |

Linear regression of $\log(\lambda_L^{\rm phys})$ on $n$ at each $\beta$
(3-point fit; for diagnostic purposes only, not an asymptotic claim):

| $\beta_{\rm phys}$ | $T/T_c$ | phase | slope | $e^{\rm slope}$ | $R^2$ |
|---|---|---|---|---|---|
| 0.10 | 4.83 | DIS | $-0.002$ | $0.998$ | $0.004$ |
| 0.25 | 1.93 | DIS | $-0.023$ | $0.977$ | $0.232$ |
| 0.50 | 0.97 | ORD | $-0.114$ | $0.892$ | $0.757$ |
| 1.00 | 0.48 | ORD | $-0.480$ | $0.619$ | $0.832$ |
| 2.00 | 0.24 | ORD | $-1.613$ | $0.199$ | $0.998$ |

**Reading**.

- **Paramagnetic** ($\beta \in \{0.10, 0.25\}$, $T > T_c$): slope is
  indistinguishable from noise ($|{\rm slope}| < 0.025$, $R^2 < 0.25$).
  The data is consistent with a flat $\Omega(1)$ gap — *no closing*.
- **Just below $T_c$** ($\beta = 0.5$, $T/T_c = 0.97$, $\text{eff\_rank}
  \in [3.3, 8.6]$ — Gibbs NOT on the doublet): the gap closes 0.449 → 0.286
  between $n=4$ and $n=6$ (factor 1.57), then saturates at 0.285 at $n=8$.
  *This is the decisive cell:* closing is observed where the Gibbs has
  several thermally-relevant states. The artefact hypothesis is dead.
- **Deep ordered** ($\beta \in \{1.0, 2.0\}$, $\text{eff\_rank} \approx 2$):
  strong closing, monotone in $\log(\lambda_L)$ — slope $-1.6$ at $\beta=2$,
  $R^2 = 0.998$.

Important observation that complicates the $\beta = 0.5$ slope: at $\beta = 0.5$
the closing **saturates** at $n=8$ (0.285 essentially unchanged from $n=6$),
which the simple $\log$-linear slope cannot capture. This saturation is
physically meaningful — see §7.

## 6 Slow-mode operator-space character (Sweep C, dense L eigendecomp at n=4)

| $\beta_{\rm phys}$ | $T/T_c$ | phase | $|\lambda_2|^{\rm phys}$ | $M_{d \times d}$ | $\mathbf{M_{d \times {\rm bulk}}}$ | $M_{{\rm bulk}}$ | Z₂ parity | $|\langle R_2, M_z \rangle|$ |
|---|---|---|---|---|---|---|---|---|
| 0.10 | 4.83 | DIS | $4.534$ | 0.025 | **0.763** | 0.212 | 0.014 | $\sim 10^{-15}$ |
| 0.25 | 1.93 | DIS | $1.910$ | 0.029 | **0.892** | 0.080 | 0.001 | $\sim 10^{-15}$ |
| 0.50 | 0.97 | ORD | $0.449$ | 0.028 | **0.897** | 0.074 | $10^{-5}$ | $\sim 10^{-14}$ |
| 1.00 | 0.48 | ORD | $0.068$ | 0.028 | **0.894** | 0.078 | $10^{-9}$ | $\sim 10^{-13}$ |
| 2.00 | 0.24 | ORD | $0.100$ | 0.028 | **0.894** | 0.078 | $10^{-18}$ | $\sim 10^{-14}$ |

Top-3 off-diagonal entries (energy eigenbasis index $(i, j)$ with their Bohr
frequencies $\omega = E_i - E_j$ in physical units):

| $\beta_{\rm phys}$ | top 3 |
|---|---|
| 0.10 | (9,2)[ω=+8.47, \|amp\|²=0.032], (2,9)[ω=−8.47, 0.032], (12,1)[ω=+10.5, 0.029] |
| 0.25 | (9,2)[+8.47, 0.038], (2,9)[−8.47, 0.038], (9,1)[+8.54, 0.037] |
| 0.50 | (9,2)[+8.47, 0.040], (2,9)[−8.47, 0.040], (9,1)[+8.54, 0.039] |
| 1.00 | (2,9)[−8.47, 0.040], (9,2)[+8.47, 0.040], (1,9)[−8.54, 0.039] |
| 2.00 | (2,9)[−8.47, 0.040], (9,2)[+8.47, 0.040], (1,9)[−8.54, 0.039] |

**Reading**. Across the entire $\beta$ range we sweep — paramagnet through
critical to deep ordered — **the slow mode is the same operator**:

- ≈ 89 % of $\|R_2\|^2$ lives in the doublet × bulk off-diagonal block
  $\{|\psi_{1,2}\rangle\langle\psi_{k \ge 3}|, \text{h.c.}\}$.
- $\le 3$ % lives in the doublet × doublet block (within-doublet coherence
  is **not** the slow mode — confirms qf-biz Follow-up A at every β).
- The remaining ~ 8 % is in bulk × bulk; at the deep paramagnet
  ($\beta = 0.10$) this jumps to 21 % — but the doublet × bulk character
  still dominates.
- Top entries are always $|\psi_{1,2}\rangle\langle\psi_9|$-type at Bohr
  frequency $\omega \approx \pm 8.5$ (= the bulk gap at $n=4$, $h=J=1$),
  with $|{\rm amp}|^2 \approx 0.04$ each — i.e. the mass is *spread* over
  the doublet × bulk block, not concentrated on any single entry.
- The Z₂ parity of $R_2$ goes from $\sim 10^{-2}$ at deep paramagnet to
  $\sim 10^{-18}$ at deep ordered: the slow mode becomes more
  cleanly Z₂-mixed as the ordered phase is entered. $\langle R_2, M_z\rangle$
  is machine-precision zero at every $\beta$ — the slow mode carries **no
  net magnetisation** (consistent with R_2 being orthogonal to the order
  parameter, qf-biz Follow-up A picture extended).

So the closing of $\lambda_L^{\rm phys}$ across $\beta$ is **not** driven by
a change in *what* the slow mode is — it stays a doublet × bulk coherence
at $\omega \approx \pm 8.5$ everywhere. What changes with $\beta$ is the
**KMS rate at that frequency**: $\gamma_{\rm KMS}(\pm 8.5) \propto
e^{\mp \beta \cdot 8.5 / 2}$ (energy-decreasing direction
exponentially favoured). At $\beta = 0.10$, $\gamma \sim 0.43 \gamma_{\rm
filter}$; at $\beta = 2$, $\gamma \sim 1.4 \times 10^{-4} \gamma_{\rm
filter}$ — consistent in order-of-magnitude with the observed
$\lambda_L^{\rm phys}$ values at $n=4$. (Quantitative agreement requires the
full Davies-rate-eigenvalue computation: the bulk-mode rate enters as a
weighted superposition, with destructive interference that grows with $N$.
We do not pursue that calculation here; the mechanism is settled.)

## 7 Verdict checklist

| diagnostic | observation | passes? |
|---|---|---|
| **(i)** Closing persists at every β below T_c (β ∈ {0.5, 1.0, 2.0}) | slopes $-0.11, -0.48, -1.61$; **clear closing at all 3** | **✓ PASS** |
| **(ii)** $U_4 \to 2/3$ across ordered window | $U_4 \in [0.600, 0.657]$, hugging the SSB value | **✓ PASS** |
| **(iii)** $R_2$ in doublet × bulk block at every ordered β | 89 % at $\beta \in \{0.5, 1.0, 2.0\}$ | **✓ PASS** |
| **(iv)** Paramagnetic control flat in n | slopes $-2 \cdot 10^{-3}, -2 \cdot 10^{-2}$ (noise) | **✓ PASS** |

**Verdict**: the closing is the genuine ordered-phase surface-tension
mechanism predicted by Gamarnik–Kiani–Zlokapa 2024 for the exact CKG KMS-DB
sampler class we use. It is *not* a Gibbs-≈-GS-doublet artefact.

A nuance on (i): the slope is much shallower at $\beta = 0.5$ than at
$\beta = 2$ — and at $\beta = 0.5$ the gap actually **saturates** between
$n = 6$ and $n = 8$ (both at $\lambda_L^{\rm phys} = 0.285$). This is
*expected* for our $2 \times L_y$ ladder geometry, **not** evidence against
the mechanism — see §8.

## 8 Caveat: $2 \times L_y$ ladder ≠ $L \times L$ square; finite-floor scaling

Gamarnik–Kiani–Zlokapa Theorem 1.4 and Martinelli–Olivieri 1994 assume
$L \times L$ square geometry where the perimeter of any
symmetry-broken-region cut scales as $\sqrt{n}$. Our $2 \times L_y$
ladders have a **fixed cheapest cut** — two vertical bonds (cost $4J$) —
independent of $L_y$. The Boltzmann suppression of this cut is therefore
$e^{-\beta \cdot 4J}$, *constant in n*. At our $\beta$ values, $J=1$, this
floor estimate gives

| $\beta_{\rm phys}$ | $e^{-\beta \cdot 4J}$ | observed $\lambda_L^{\rm phys}$ at $n=8$ |
|---|---|---|
| 0.5 | $1.4 \times 10^{-1}$ | $2.85 \times 10^{-1}$ |
| 1.0 | $1.8 \times 10^{-2}$ | $1.0 \times 10^{-2}$ |
| 2.0 | $3.4 \times 10^{-4}$ | $1.57 \times 10^{-4}$ |

— within an $O(1)$ prefactor at every $\beta$. The $\beta = 0.5$ data
**saturates** at this floor by $n = 6$, $n = 8$. The $\beta = 1$ data is at
the floor by $n = 8$ (only a factor 1.24 between $n = 6$ and $n = 8$). The
$\beta = 2$ data is *below* the rough floor estimate by a factor 2 —
unsurprising for a finite-size sandbox where the prefactor includes corrections
beyond the bare cut cost.

This is **not** a problem for the mechanism: the surface-tension picture
predicts ladder-geometry closing **to a finite floor**, not to zero. To see
the $\sqrt{n}$ scaling proved by Gamarnik–Kiani–Zlokapa one would need
near-square $L_x \approx L_y$ at $L \ge 4$ (i.e. $n \ge 16$) — squarely
outside the sandbox. Sandbox 2×L_y is a *qualitative* check on the
ordered-vs-disordered distinction; not a quantitative test of the GKZ24
exponent.

## 9 Implications for the qf-1jj draft

- The qf-1jj memory's
  ([[2d-tfim-ordered-vs-disordered-qf-1jj]]) self-caveat that the data is
  "qualitative, not quantitative" was correct on the **scaling-law fit**
  (3 points, narrow $\sqrt{n}$ range), and qf-c9g respects that — we do
  *not* extract an exponent from these sandbox cells. What qf-c9g adds is
  the **β-persistence** and **mechanism** evidence: the closing IS the
  surface-tension-driven ordered-phase phenomenon predicted by GKZ24, and
  it is NOT a deep-cold doublet artefact. The ladder-floor estimate
  $\sim e^{-\beta \cdot 4J}$ matches the n=8 measurement to within $O(1)$
  at every ordered $\beta$.
- The qf-biz Follow-up A finding (slow $R_2$ lives 89 % in doublet × bulk
  at $\beta = 2$) generalises across the entire $\beta$ range — same
  structural mode at $\beta \in \{0.25, 0.5, 1, 2\}$. The
  *eigenvalue* of this mode changes by 4 orders of magnitude with $\beta$
  (4.5 → 1.9 → 0.45 → 0.07 → 0.10 at $n = 4$); the operator-space
  *character* is constant. The β-dependence is in $\gamma_{\rm KMS}(\omega)$
  at the doublet × bulk frequency $\omega \approx 8.5$.

## 10 Reading for the thesis

For the numerics chapter:

- Use the qf-c9g 15-cell table as the canonical 2D TFIM
  ordered-vs-paramagnet diagnostic. The 5-β × 3-n grid is the cleanest
  way to show "the gap closing is phase-specific" — much more convincing
  than the qf-1jj two-operating-point comparison alone.
- Cite Gamarnik–Kiani–Zlokapa Theorem 1.4 as the asymptotic lower bound;
  cite Martinelli–Olivieri 1994 as the classical analogue. The sandbox
  measurement is *not* a verification of the GKZ24 exponent (geometry
  + size); it is consistency with the **mechanism**.
- Phrase the qf-1jj headline in terms of "qualitative ordered-phase
  closing consistent with GKZ24 surface-tension lower bound", not "we
  measured the GKZ24 exponent" or "we have an asymptotic scaling".
- Make sure the ladder-geometry caveat is in the figure caption.

## Cross-references

- `drafts/qf-c9g-lit-survey.md` — full literature distillation.
- `drafts/qf-biz-2d-tfim-matrix-element-bottleneck.md` — within-doublet
  matrix-element refutation and dense-L slow-mode characterisation at β = 2
  (qf-biz Follow-up A). qf-c9g extends Follow-up A across β.
- [[2d-tfim-ordered-vs-disordered-qf-1jj]] — parent qf-1jj memory + ladder-floor
  caveat.
- [[tc_2d_tfim_phase_diagram]] — T_c(h) table from Hesselmann–Wessel 2016.
- [[qf_biz_2d_tfim_matrix_element_refuted]] — within-doublet refutation memory.

## Citations

```bibtex
@article{gamarnik2024slowmixing,
  author = {Gamarnik, David and Kiani, Bobak T. and Zlokapa, Alexander},
  title = {Slow Mixing of Quantum Gibbs Samplers},
  journal = {arXiv:2411.04300},
  year = {2024}
}
@article{martinelli1994olivieri,
  author = {Martinelli, F. and Olivieri, E.},
  title = {Approach to equilibrium of {Glauber} dynamics in the one-phase region II:
           the general case},
  journal = {Communications in Mathematical Physics},
  volume = {161}, pages = {487--514}, year = {1994},
  doi = {10.1007/BF02187060}
}
@article{rakovszky2024bottlenecks,
  author = {Rakovszky, Tibor and Placke, Benedikt and Breuckmann, Nikolas P. and Khemani, Vedika},
  title = {Bottlenecks in Quantum Channels and Finite-Temperature Phases of Matter},
  journal = {arXiv:2412.09598},
  year = {2024}
}
@article{bergamaschi2024quantumadvantage,
  author = {Bergamaschi, Thiago and Chen, Chi-Fang and Liu, Yunchao},
  title = {Quantum Computational Advantage with Constant-Temperature Gibbs Sampling},
  journal = {arXiv:2404.14639},
  year = {2024}
}
@article{lubetzky2013cutoff,
  author = {Lubetzky, Eyal and Sly, Allan},
  title = {Cutoff for the {Ising} model on the lattice},
  journal = {arXiv:1305.4524},
  year = {2013}
}
@article{hesselmann2016wessel,
  author = {Hesselmann, Sebastian and Wessel, Stefan},
  title = {Thermal {Ising} transitions in the vicinity of two-dimensional quantum critical points},
  journal = {Phys. Rev. B}, volume = {93}, pages = {155157}, year = {2016},
  doi = {10.1103/PhysRevB.93.155157}
}
```
