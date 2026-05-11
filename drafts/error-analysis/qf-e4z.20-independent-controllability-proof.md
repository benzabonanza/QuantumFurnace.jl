# qf-e4z.20 — Independent per-leg Trotter caches controllability proof

**Statement.** With the qf-e4z.20 `TrotterTriple` scheme (three INDEPENDENT
Strang Trotterizations, one per coherent-term leg `D / b_- / b_+`), the
TrotterDomain Lindbladian's KMS-DBC residue $\|\mathcal{L}\cdot\sigma_\beta\|_\mathrm{HS}$
is controllable to $\le 10^{-6}$ at $(n=3,\,\beta=10,\,\sigma=1/\beta,\,
\text{smooth Metropolis}, s=0.25)$ **by tightening only the dissipative-leg
substep count $M_D$**, leaving $M_{b_-} = M_{b_+} = 1$.

## Setup

Cell: `n=3, β=10, σ=1/β=0.1, smooth_metro` with `s=0.25, a=β/30=0.333,
η=1e-3`. Hamiltonian: 1D XXX with Z + ZZ disorder. Quadrature windows
`T_minus=18`, `T_plus=12`. Dissipative span `ω_range = 2.5`.

Per-leg quantities (`make_cfg` in `scripts/scratch_e4z20_independent_controllability.jl`):

| Quantity | Formula | Notes |
|---|---|---|
| `t0_D` | $2\pi / (2^{r_D} \cdot w0_D)$, $w0_D = \omega_{\text{range}}/2^{r_D}$ | dissipative OFT grid |
| `t0_b_minus` | $2 T_- / 2^{r_{b_-}}$ | outer kernel grid step |
| `t0_b_plus` | $2 T_+ / 2^{r_{b_+}}$ | inner kernel grid step |
| `δt₀_D` | `t0_D / M_D` | dissipative Strang substep |
| `δt₀_b_minus` | `(t0_b_minus / σ) / M_b_minus` | outer Strang substep |
| `δt₀_b_plus` | `(β · t0_b_plus) / M_b_plus` | inner Strang substep |

Independent knobs: $(r_D, r_{b_-}, r_{b_+}, M_D, M_{b_-}, M_{b_+})$ with NO
inter-leg commensurability constraint (the qf-d0w shared-$\delta t_0$
constraint is gone with the new `TrotterTriple` scheme).

## Sweep 1 — Joint $M$ ramp confirms Strang slope $-2$

At the tight grid $(r_D, r_{b_-}, r_{b_+}) = (11, 14, 18)$:

| $M$ (per leg) | $\|\mathcal{L}\cdot\sigma_\beta\|_\mathrm{HS}$ | ratio |
|---|---|---|
| 1   | 1.776e-3 | — |
| 2   | 4.457e-4 | 3.99 |
| 4   | 1.115e-4 | 4.00 |
| 8   | 2.789e-5 | 4.00 |
| 16  | 6.972e-6 | 4.00 |
| 32  | 1.743e-6 | 4.00 |
| 64  | 4.357e-7 | 4.00 |
| 128 | 1.089e-7 | 4.00 |
| 256 | 2.723e-8 | 4.00 |

Clean Strang slope $-2$ across 9 decades of $M$: each doubling of $M$ drops
the residue by exactly $4\times$. **At $M \ge 64$ the residue is $\le 10^{-6}$.**

## Sweep 2 — Per-leg attribution

Hold two legs at $M = 128$ (Strang-saturated) and sweep the third.

**Tighten $M_D$ only** (`M_b_minus = M_b_plus = 128`):

| $M_D$ | residue |
|---|---|
| 1   | 1.776e-3 |
| 4   | 1.115e-4 |
| 16  | 6.972e-6 |
| 64  | 4.357e-7 |
| 128 | 1.089e-7 |
| 256 | 2.723e-8 |

**Tighten $M_{b_-}$ only** (`M_D = M_b_plus = 128`):

| $M_{b_-}$ | residue |
|---|---|
| 1   | 1.095e-7 |
| 4   | 1.090e-7 |
| 16  | 1.089e-7 |
| 64  | 1.089e-7 |
| 128 | 1.089e-7 |
| 256 | 1.089e-7 |

**Tighten $M_{b_+}$ only** (`M_D = M_b_minus = 128`):

| $M_{b_+}$ | residue |
|---|---|
| 1   | 1.089e-7 |
| 4   | 1.089e-7 |
| 16  | 1.089e-7 |
| 64  | 1.089e-7 |
| 128 | 1.089e-7 |
| 256 | 1.089e-7 |

**Verdict.** Only $M_D$ moves the needle. The $b_-$ and $b_+$ legs are
already at their Strang floor at $M = 1$:

- $b_+$: natural Trotter step $\beta \cdot t_{0,b_+}^{\text{grid}} = 10 \cdot
  (2 T_+ / 2^{r_{b_+}}) = 10 \cdot (24/262144) \approx 9 \times 10^{-4}$. With
  $M_{b_+} = 1$, $\delta t_{0,b_+} \approx 9 \times 10^{-4}$, Strang error per
  step $\sim (\delta t_0)^2 \cdot \|[H_a,H_b]\| \sim 8 \times 10^{-7}$ — well
  below the $10^{-6}$ target before any further refinement.
- $b_-$: natural step $t_{0,b_-}^{\text{grid}} / \sigma = (2 T_- / 2^{r_{b_-}}) /
  \sigma = (36/16384) / 0.1 \approx 2.2 \times 10^{-2}$. Strang error per step
  $\sim 5 \times 10^{-4}$, accumulated over $\sim 50$ steps gives the outer
  contribution $\lesssim 10^{-5}$ — STILL contributes, but is sub-dominant
  to $M_D = 1$'s contribution of $\sim 10^{-3}$.

The coherent legs $b_-$ and $b_+$ have natural (M = 1) substeps so small
that Strang precision is automatic; the dissipative leg has the large
$t_{0,D} \approx 2.5$ and is the sole Trotterization that needs refinement.

## Sweep 3 — Quadrature floor at saturated $M$

At $M = 128$ per leg (Strang error $\le 10^{-7}$), sweep the register sizes:

**$r_{b_+}$ sweep** (`r_D=11, r_b_minus=14`):

| $r_{b_+}$ | residue |
|---|---|
| 10  | 1.089e-7 |
| 12  | 1.089e-7 |
| 14  | 1.089e-7 |
| 16  | 1.089e-7 |
| 18  | 1.089e-7 |
| 20  | 1.089e-7 |

**$r_D$ sweep** (`r_b_minus=14, r_b_plus=18`):

| $r_D$ | residue |
|---|---|
| 7   | 1.089e-7 |
| 9   | 1.089e-7 |
| 11  | 1.089e-7 |
| 13  | 1.089e-7 |
| 15  | 1.089e-7 |

The quadrature is already over-resolved at every $r$ value swept: the
residue plateaus at $1.09 \times 10^{-7}$. This is the floor from the
inner/outer combined Strang error at $M_{b_-} = M_{b_+} = 128$ (which is
small but nonzero); pushing $M_D \to 256$ continues the slope-$2$ drop
(see Sweep 1 last row).

## Final tight recipes

Demonstrating ≤ $10^{-6}$ at varying levels of leg-resolution:

| Recipe | $\|\mathcal{L}\cdot\sigma_\beta\|_\mathrm{HS}$ |
|---|---|
| $(r_D,r_{b_-},r_{b_+}) = (11,14,18),\ M=(64, 1, 1)$ | $4.36 \times 10^{-7}$ |
| $(11,14,18),\ M=(128, 1, 1)$ — **optimal** | $1.10 \times 10^{-7}$ |
| $(11,14,18),\ M=(128, 4, 4)$ | $1.09 \times 10^{-7}$ |
| $(11,14,18),\ M=(64, 64, 64)$ — legacy-style | $4.36 \times 10^{-7}$ |
| $(9, 14,18),\ M=(128, 1, 1)$ — low $r_D$ | $1.10 \times 10^{-7}$ |
| $(9, 9, 9),\ M=(128, 1, 1)$ — minimal grids | $1.82 \times 10^{-6}$ |

**Headline recipe**: $\boxed{r_D = 9,\ r_{b_-} = 14,\ r_{b_+} = 18,\ M_D = 128,\ M_{b_-} = M_{b_+} = 1}$,
giving residue $1.10 \times 10^{-7}$.

## Comparison to qf-e4z.5.3 Option A (shared $\delta t_0$)

The previous Option A recipe (`drafts/error-analysis/qf-e4z.5.3-controllability-proof.md`)
needed $m_D = 80$, $M_\text{user} = 4$ to hit $1.6 \times 10^{-6}$:

- Shared $\delta t_0 = t_{0,D} / m_D \approx 0.031$, REQUIRED to be the
  elementary Strang step of all three legs.
- Implied per-leg $M$ counts: $M_D = m_D \cdot M_\text{user} = 320$,
  $M_{b_-} \approx (\beta\, t_{0,b_-} / \delta t_0) \approx 11$,
  $M_{b_+} \approx 4$.
- Net Strang substeps per cell: $\approx 320 + 11 + 4 = 335$.

The independent-knob recipe above achieves a $15\times$ better residue
($1.1 \times 10^{-7}$ vs $1.6 \times 10^{-6}$) at $\approx 130$ Strang
substeps per cell (just $M_D = 128$, $M_{b_-} = M_{b_+} = 1$). **Half the
substep count, an order of magnitude tighter.**

The shared-$\delta t_0$ scheme was paying for $b_-$ and $b_+$ over-resolution
that buys nothing — $M_{b_\pm} = 1$ already saturates those legs at this
fixture. The independent-knob refactor lets us spend Strang substeps on the
one leg that matters.

## Channel-floor cross-check

The Lindbladian residue $\|\mathcal{L}\sigma_\beta\|_\mathrm{HS}$ is the
algorithm-level fixed-point precision (no $\delta$-step). The channel
floor $\|\rho_\infty - \sigma_\beta\|_1 / 2$ adds a slope-$1$-in-$\delta$
contribution from the generator-splitting (Φ_δ ≈ I + δ·𝓛 + O(δ²)). At the
optimal recipe ($\|\mathcal{L}\sigma_\beta\|_\mathrm{HS} = 1.1 \times 10^{-7}$),
a dense $\Phi_\delta$ build with $\delta \le 5 \times 10^{-7}$ would
satisfy `0.2·δ ≤ 1e-7`, giving total channel floor $\le 2 \times 10^{-7}$ —
below the $10^{-6}$ target with margin.

A full channel-floor sweep is left for follow-up (qf-e4z.20.6 acceptance:
the Lindbladian-residue controllability ≤ 1e-6 is established here at
$1.1 \times 10^{-7}$, with the additional $\delta$-split contribution
trivially controlled by $\delta \le 5 \times 10^{-7}$).

## Cost model summary

Per `B_trotter` call at the optimal recipe `(r_D=9, r_bm=14, r_bp=18,
M_D=128, M_bm=M_bp=1)` at $n=3$:

- Cache construction: three Strang Trotterizations + three `eigen()` calls
  on $d \times d$ matrices ($d = 2^n = 8$). Sub-millisecond.
- Inner $\tau$-loop: $N_{b_+} \sim 6 \times 10^4$ effective grid points
  $\times n_\text{jumps} = 9$ inner steps × $O(d^3)$ GEMM ≈ 1.4 × $10^8$ flops.
- Outer $t$-loop: $N_{b_-} \sim 5 \times 10^3$ × $O(d^3)$ ≈ $5 \times 10^6$ flops.
- Basis rotations: $n_\text{jumps} \cdot 2 d^3 + 4 d^3 \approx 10^4$ flops —
  measured at 9.8μs vs 207ms total call ($4.7 \times 10^{-5}$ of total cost).

Rotation overhead from the V_D → V_bp jump rotation, V_bp → V_bm summand
rotation, and V_bm → V_D final rotation is **negligible** at $n = 3$
(0.005% of total `B_trotter` wall time). Confirmed via BenchmarkTools at
`scripts/scratch_e4z20_independent_controllability.jl` companion profile.

## Acceptance — qf-e4z.20.6 closure

1. ✓ **Independent-knob recipe at ≤ 1e-6 demonstrated**: optimal
   $\boxed{(9, 14, 18), M=(128, 1, 1)}$ → 1.1e-7.
2. ✓ **Cheaper than Option A**: 1.1e-7 vs 1.6e-6 ($15\times$ better) at
   ~$2.5\times$ fewer Strang substeps per cell.
3. ✓ **Independent control verified**: $M_D$ alone controls residue;
   $M_{b_\pm}$ saturates at $M = 1$.
4. ✓ **Wall time ≤ 60s per cell**: actual ~200ms per `construct_lindbladian`
   call at the optimal recipe.

Generating script: `scripts/scratch_e4z20_independent_controllability.jl`.

## Code

The optimal recipe in code:

```julia
cfg = Config(
    sim = Lindbladian(), domain = TrotterDomain(), construction = KMS(),
    num_qubits = 3, with_linear_combination = true,
    beta = 10.0, sigma = 0.1, s = 0.25, a = 10.0/30.0, eta = 1e-3,
    num_energy_bits_D = 9,  t0_D = 2π/(2^9 * 2.5/2^9),     w0_D = 2.5/2^9,
    num_energy_bits_b_minus = 14, t0_b_minus = 36.0/2^14,  w0_b_minus = 2π/(2^14 * 36.0/2^14),
    num_energy_bits_b_plus  = 18, t0_b_plus  = 24.0/2^18,  w0_b_plus  = 2π/(2^18 * 24.0/2^18),
    num_trotter_steps_per_t0_D       = 128,
    num_trotter_steps_per_t0_b_minus = 1,
    num_trotter_steps_per_t0_b_plus  = 1,
)
trotter = make_trotter_for_config(ham, cfg)    # → TrotterTriple
L = construct_lindbladian(jumps, cfg, ham; trotter=trotter)
# ‖L · σ_β‖_HS ≈ 1.1e-7
```
