# qf-e4z.35 σ-sweep — structural mixing pattern (n=3..7)

**Question:** at what σ value does the smooth-Metropolis CKG sampler
have the best structural mixing for the project-canonical initial state
$\rho_0 = |+\rangle\!\langle +|^{\otimes N}$, after normalising the
trajectory mixing time by a $\sigma$-dependent norm of $\mathcal L$?

**Setup:** the qf-e4z.35 σ-sweep produces per-cell `mixing_time_phys`,
`gap_phys`, `‖L‖_HS_phys`, and `d_{1→1}_phys`. The last is the operator
norm of the Kossakowski-weighted jump-overlap matrix $M_\alpha[k,j] =
\sum_a \sum_i \overline{A^a_{ik}} A^a_{ij}\, \alpha(\nu_{ik}, \nu_{ij})$
(see `scripts/scratch_qf_e4z_35_sigma_sweep_plot.jl:294`–298); it
upper-bounds the dissipator's induced 1→1 superoperator norm via
$\|\mathcal L_{\mathrm{diss}}\|_{1\to 1} \le 4\,\|M_\alpha\|_{\mathrm{op}}$
but is NOT itself that induced norm, NOR a $1/\sigma$ rate scale —
empirically it varies sub-$1/\sigma$ in $c$. We use two normalisations
($\|\mathcal L\|_{\mathrm{HS}}$ and $d_{1\to 1}$) precisely *because*
neither is canonical: their agreement on $c^\star$ tells us the trend
is robust under the choice of norm, not that we have identified "the"
rate. **A separate follow-up, qf-j4a, is filed to characterise the
right σ-independent norm for rate-removed comparisons.** Until that
lands, $\tau_{\mathrm{mix}} \cdot \mathcal{N}$ should be read as a
geometric *quantity per unit of $\mathcal{N}$*, not as the rate-removed
mixing time.

Analysis script: `scripts/scratch_qf_e4z_35_structural_pattern.jl`.
Full tabular output: `drafts/qf-e4z-35-structural-pattern-n3-7.txt`.

## c* (argmin of τ_mix · ‖L‖) tables

Both L-norms agree to within ±1 grid step everywhere — the structural
optimum is robust under norm choice.

**τ_mix · ‖L‖_HS:**

| β_phys \\ n | n=3 | n=4 | n=5 | n=6 | n=7 |
|------------|-----|-----|-----|-----|-----|
| 0.25 (hot) | 0.25 | 0.75 | 0.75 | 0.75 | 0.75 |
| 0.5        | 0.25 | 1.50 | 1.00 | 1.50 | 1.50 |
| 1.0 (cold) | 0.25 | 2.00 | 1.50 | 2.00 | 2.00 |

**τ_mix · d_{1→1}:**

| β_phys \\ n | n=3 | n=4 | n=5 | n=6 | n=7 |
|------------|-----|-----|-----|-----|-----|
| 0.25 (hot) | 0.25 | 0.75 | 0.75 | 0.75 | 0.75 |
| 0.5        | 0.25 | 1.50 | 1.00 | 1.50 | 1.50 |
| 1.0 (cold) | 0.25 | 2.00 | 1.50 | 1.50 | 1.50 |

## Findings

1. **n=3 is anomalous** — every β_phys cell prefers c*=0.25. Likely because
   the d=8 spectrum is too small / dense for a wider smoothing kernel
   to be helpful: every nontrivial Bohr frequency is already
   well-separated, so the narrowest kink wins.

2. **For n ≥ 4: optimal σ grows with β_phys** (colder ⇒ wider σ is
   structurally better):
   - β=0.25:  c* ≈ 0.75
   - β=0.5:   c* ≈ 1.0–1.5
   - β=1.0:   c* ≈ 1.5–2.0 (often at the upper edge of our c-grid — the
     true optimum may lie at c > 2 for β=1)
   The optimum is essentially **n-independent** within each β bin —
   this is a temperature effect, not a finite-size effect.

3. **Structural spread is modest:** the max/min of τ_mix·‖L‖ across the
   c-grid is **1.25–1.7×** at n ≥ 4. Picking the right σ buys roughly
   25–70% structural improvement, not orders of magnitude.

4. **The canonical σ = 1/β_alg (c=1) is well-placed at β=0.5 but
   sub-optimal at both temperature extremes** — slightly too wide at
   β=0.25 (optimum c=0.75) and too narrow at β=1.0 (optimum c=1.5–2.0).

## Why the two normalisations agree

`‖L‖_HS` is the Hilbert–Schmidt norm of the full generator (dissipator
+ coherent commutator); `d_{1→1}` is the Kossakowski opnorm
$\|M_\alpha\|_{\mathrm{op}}$ defined above, an upper bound on the
dissipator's induced 1→1 norm. They agree on $c^\star$ because (a) the
coherent part is approximately σ-independent (the Hamiltonian
commutator $i[H + B, \cdot]$ with $B$ from a σ-dependent integrand that
integrates out to a near-σ-independent operator), and (b) the
dissipator-side σ-dependence dominates the σ-scaling at the (n, β)
cells we sweep. Both norms therefore track the same underlying
dissipator-geometry quantity within an O(1) constant; their $c^\star$
agreement is corroborating evidence, not independent normalisations.

## Different metrics, different answers

- **τ_mix (absolute mixing time):** prefers c=0.25 everywhere because
  smaller σ ⇒ larger $\|\mathcal L\|$ (in either norm) ⇒ faster
  wall-clock mixing. This is the σ-dependence of the generator's
  overall scale, which the σ-normalised quantities below subtract off.
- **τ_mix · ‖L‖_HS / τ_mix · d_{1→1} (σ-normalised):** see above —
  optimum migrates from c=0.75 (hot) to c=1.5–2.0 (cold) for n ≥ 4.
  Neither norm is canonical — see qf-j4a — but their agreement
  shows the geometric trend is robust under the choice.
- **τ_mix · gap (rate-resolved):** picks c*=2.0 almost uniformly. This
  measures the $|c_2|$ coefficient of the $|+\rangle\!\langle+|^{\otimes N}$
  initial state on the slowest mode — at wider σ, that overlap shrinks
  (the swap), so τ_mix per gap drops. Different question from "is L
  structurally better"; see also
  [`drafts/qf-e4z-35-spectral-mechanism.md`](qf-e4z-35-spectral-mechanism.md).
- **gap / ‖L‖_HS (relative spectral gap):** peaks at c ≈ 0.25–0.75 at
  n ≥ 4 — smaller σ ⇒ slow mode more separated from the bulk in
  relative terms — but this isn't what determines τ_mix (because $|c_2|$
  for the chosen initial state matters, not just the gap).

So **for a claim of "best σ-normalised mixing for $\rho_0 = |+\rangle^{\otimes N}$"**, the
two L-norm normalisations (HS or $d_{1\to 1}$) point at the same
$c^\star$. The answer:

> For the canonical 1D Heisenberg PBC fixtures at canonical β_phys ≤ 1
> and n ≥ 4: the σ-normalised optimum (under either $\|L\|_{\mathrm{HS}}$
> or $d_{1\to 1}$) is a wider-than-canonical σ that grows with β_phys
> (from c ≈ 0.75 at the hottest cell to c ≈ 2 at the coldest), nearly
> system-size independent.

**The V-shape is a property of $(\mathcal L, \rho_0)$, not of $\mathcal L$
alone.** The mode-swap that drives the V-shape (see the companion
spectral-mechanism draft) is intrinsic to the $\mathcal L(c)$ family,
but the $|c_2|$ collapse that makes the slowest mode irrelevant for
trace distance requires $\rho_0$ to be near-orthogonal to the
post-swap coherence mode — which is exactly what $|+\rangle^{\otimes N}$
is. A different $\rho_0$ with weight on both branches would smooth out
the V-shape. The $c^\star$ values above should therefore be read as
"best σ-normalised $c$ *for $\rho_0 = |+\rangle\!\langle+|^{\otimes N}$*",
matching the project's canonical-$\rho_0$ policy, not as a property of
the Lindbladian family on its own.

## What's next (when n=8 resumes)

The n=8 cells will pin down whether the n-independence of c* persists
into the larger system. The β=1 optimum at c=2.0 is at the grid edge —
a follow-up sweep extending the c-grid (e.g. c ∈ {2, 3, 4}) at β=1
would tell us whether the cold-cell optimum sits inside the grid or
keeps growing.

## Why the V-shape exists — mechanism

The structural V-shape in $\tau_{\mathrm{mix}}\cdot d_{1\to 1}$ (for
$n \ge 4$) is an *avoided eigenvalue crossing* in the L spectrum at
fixed $\beta$, threaded by $c$. Below $c^\star$ the slowest
L-eigenmode is a population mode (diagonal in the energy eigenbasis)
that $\rho_0 = |+\rangle^{\otimes N}$ maximally excites; above
$c^\star$ a coherence mode becomes slowest and $\rho_0$ has near-zero
overlap with it, so the operationally effective mixing rate jumps to
the next-slowest mode. The argmin of $\tau_{\mathrm{mix}}\cdot
d_{1\to 1}$ in $c$ is co-located with this swap. See
[`drafts/qf-e4z-35-spectral-mechanism.md`](qf-e4z-35-spectral-mechanism.md)
for the direct evidence (qf-e4z.35.1: $|c_2|^2$ drops 3–6 orders of
magnitude across one grid step at every $n \ge 4$ cell;
$\mathrm{off\text{-}diag}(R_2)$ is bimodal $\{0, 1\}$). The $n=3$ cell
is excluded from that mechanism analysis — at $d=8$ no separated
population sub-band can exist, the slow mode is already a coherence
mode at every $c$ in the grid, and the V-shape is governed by
different sub-dominant-mode physics that the mechanism diagnostic
isn't probing here.
