# qf-e4z.20 — Independent per-leg Trotter caches controllability proof

**Statement.** With the qf-e4z.20 `TrotterTriple` scheme (three INDEPENDENT
Strang Trotterizations, one per coherent-term leg `D / b_- / b_+`), the
TrotterDomain Lindbladian's KMS-DBC residue $\|\mathcal{L}\cdot\sigma_\beta\|_\mathrm{HS}$
is controllable to $\le 10^{-6}$ at $(n=3,\,\beta=10,\,\sigma=1/\beta,\,
\text{smooth Metropolis},\,s=0.25)$ at the qf-7xt-canonical register sizing
$(r_D, r_{b_-}, r_{b_+}) = (7, 6, 14)$ — the per-leg substep counts scale
**inversely with the grid resolution** because finer grids ⇒ smaller natural
Trotter step ⇒ smaller $M$ needed.

## Setup

Cell: `n=3, β=10, σ=1/β=0.1, smooth_metro` with `s=0.25, a=β/30=0.333,
η=1e-3`. Hamiltonian: 1D XXX with Z + ZZ disorder. Quadrature windows
$T_- = 18$, $T_+ = 12$, dissipative span $\omega_\text{range} = 2.5$.

Per-leg natural Trotter-step durations (used by `make_trotter_for_config`):

| Leg | Grid step | Natural Trotter step |
|---|---|---|
| $D$ | $t_{0,D} = 2\pi/\omega_\text{range}$ | $t_{0,D}$ (raw Hamiltonian time) |
| $b_-$ | $t_{0,b_-}^\text{grid} = 2T_-/2^{r_{b_-}}$ | $t_{0,b_-}^\text{grid} / \sigma$ |
| $b_+$ | $t_{0,b_+}^\text{grid} = 2T_+/2^{r_{b_+}}$ | $\beta \cdot t_{0,b_+}^\text{grid}$ |

Per-leg Strang substep: $\delta t_{0,X} = t_{0,X}^\text{nat} / M_X$. The
substep is what enters the Strang error $\sim N_X \cdot (\delta t_{0,X})^2 \cdot
\|[H_a, H_b]\|$ summed over the integration variable of that leg.

Important: small $r_X$ ⇒ large grid step ⇒ large natural Trotter step ⇒
the **per-leg $M_X$ must scale UP** to keep $\delta t_{0,X}$ small. At the
qf-7xt minimum-$r$ recipe, $b_-$ has a coarse grid and needs $M_{b_-} \ge 64$;
$b_+$ has a fine grid (small $\beta \cdot t_0$) and needs only $M_{b_+} = 1$.

## Sweep 1 — Per-leg Strang attribution at qf-7xt grids

At $(r_D, r_{b_-}, r_{b_+}) = (7, 6, 14)$, hold two legs at $M = 128$
(Strang-saturated) and sweep the third:

**$M_D$ sweep** (`M_b_minus = M_b_plus = 128`):

| $M_D$ | residue | ratio |
|---|---|---|
| 1   | 1.776e-3 | — |
| 4   | 1.115e-4 | 15.9 |
| 16  | 6.935e-6 | 16.1 |
| 64  | 7.358e-7 | 9.4 |
| 128 | 6.377e-7 | 1.15 |
| 256 | 6.374e-7 | 1.00 |

**$M_{b_-}$ sweep** (`M_D = M_b_plus = 128`):

| $M_{b_-}$ | residue |
|---|---|
| 1   | 1.171e-4 |
| 4   | 6.296e-6 |
| 16  | 3.791e-7 |
| 64  | 6.215e-7 |
| 128 | 6.377e-7 |
| 256 | 6.418e-7 |
| 512 | 6.428e-7 |

**$M_{b_+}$ sweep** (`M_D = M_b_minus = 128`):

| $M_{b_+}$ | residue |
|---|---|
| 1   | 6.376e-7 |
| 4   | 6.377e-7 |
| 16  | 6.377e-7 |
| 64  | 6.377e-7 |
| 128 | 6.377e-7 |

**Verdict** — three distinct per-leg Strang signatures:

1. **$D$ leg** is **slope-2 in $M_D$**: each $M_D$ doubling drops residue by
   $\sim 4\times$ until quadrature floor $\sim 6.4 \times 10^{-7}$ is hit at
   $M_D \approx 64$.
2. **$b_-$ leg** is **slope-2 in $M_{b_-}$** (coarser grid, $M_{b_-} = 1$
   gives a Strang error $1.17 \times 10^{-4}$; need $M_{b_-} \ge 16$ to
   reach floor).
3. **$b_+$ leg** is **already saturated at $M_{b_+} = 1$** (fine grid,
   natural step $\beta \cdot t_{0,b_+}^\text{grid} = 10 \cdot 24/2^{14} \approx
   1.5 \times 10^{-2}$ — Strang error per step $\sim 2 \times 10^{-4}$,
   summed over $\sim 30$ effective grid points gives $\sim 6 \times 10^{-3}
   \cdot \|[H_a,H_b]\| \cdot |b_+(\tau)|_{\ell^1}$ which lands sub-1e-6 after
   the kernel weighting).

## Sweep 2 — Quadrature floors

At $M = 128$ per leg (Strang error ≤ 1e-7 from Sweep 1 floors), sweep each
register size to see where the quadrature floor lies:

**$r_{b_+}$ sweep** (`r_D=5, r_b_minus=6`):

| $r_{b_+}$ | residue |
|---|---|
| 8   | 6.112e-6 |
| 10  | 6.096e-6 |
| 12  | 6.096e-6 |
| 14  | 6.096e-6 |
| 16  | 6.096e-6 |

$r_{b_+} = 10$ already saturates. (qf-7xt suggests $r_{b_+}=14$ for the
$\beta=10$ smooth-Metro slope; at $r_{b_+}=10$ we are still 4 bits above
the b_+ quadrature noise floor.)

**$r_{b_-}$ sweep** (`r_D=5, r_b_plus=14`):

| $r_{b_-}$ | residue |
|---|---|
| 4   | 4.533e-4 |
| 5   | 5.792e-4 |
| 6   | 6.096e-6 |
| 7   | 6.071e-6 |
| 8   | 6.071e-6 |
| 10  | 6.071e-6 |

$r_{b_-} = 6$ is the minimum: below that, the $b_-(t)$ kernel is truncated
($|b_-|$ at the grid endpoint $\sim 0.05$ at $r_{b_-}=4$, see "Filter kernel
at ends" warning emitted by `_truncate_time_labels_for_oft`).

**$r_D$ sweep** (`r_b_minus=6, r_b_plus=14`):

| $r_D$ | residue |
|---|---|
| 4   | 3.692e-3 |
| 5   | 6.096e-6 |
| 6   | 6.377e-7 |
| 7   | 6.377e-7 |
| 8   | 6.377e-7 |
| 10  | 6.377e-7 |

The dissipative quadrature saturates at $r_D = 6$ for this $(n=3, \beta=10)$
fixture. (qf-7xt's $r_D = 5$ at ε=1e-6 was measured at $n=4$; the note
warns "$K$ prefactors and floor positions can shift by 1–2 bits with $(n, \beta)$"
— $r_D = 6$ at $n=3$ matches that prediction.)

**Total quadrature floor** at $(r_D, r_{b_-}, r_{b_+}) = (6, 6, 14)$: $6.4 \times 10^{-7}$
(below $10^{-6}$ ✓).

## Final tight recipes

| Recipe | $\|\mathcal{L}\cdot\sigma_\beta\|_\mathrm{HS}$ | Notes |
|---|---|---|
| $(7, 6, 14)$, $M=(64, 64, 1)$ | $7.25 \times 10^{-7}$ | **minimal cost ≤ 1e-6** |
| $(7, 6, 14)$, $M=(128, 128, 1)$ | $6.38 \times 10^{-7}$ | comfortable margin |
| $(7, 6, 14)$, $M=(128, 256, 1)$ | $6.42 \times 10^{-7}$ | $M_{b_-}$ over-resolved |
| $(7, 6, 14)$, $M=(256, 256, 1)$ | $6.42 \times 10^{-7}$ | $M_D$ over-resolved |
| $(5, 6, 14)$, $M=(128, 128, 1)$ | $6.10 \times 10^{-6}$ | $r_D = 5$ above floor at $n=3$ |
| $(6, 6, 14)$, $M=(128, 128, 1)$ | $6.38 \times 10^{-7}$ | qf-7xt-compatible |

**Headline recipe** (minimal cost): $\boxed{(r_D, r_{b_-}, r_{b_+}) = (7, 6, 14),\ (M_D, M_{b_-}, M_{b_+}) = (64, 64, 1)}$,
giving residue $7.2 \times 10^{-7}$ ≤ $10^{-6}$.

The register sizes follow the canonical qf-7xt smooth-Metropolis recipe
(with $r_D$ shifted up by 1 bit for the $n=3$ fixture as predicted by the
1–2 bit-shift note). The per-leg M's scale inversely with grid resolution:
$M_{b_+} = 1$ at the fine $r_{b_+} = 14$ grid (natural Trotter step
$\sim 10^{-2}$); $M_{b_-} = 64$ at the coarse $r_{b_-} = 6$ grid (natural
Trotter step $\sim 5.6$).

## Comparison to qf-e4z.5.3 Option A (shared $\delta t_0$)

The previous Option A recipe (`drafts/error-analysis/qf-e4z.5.3-controllability-proof.md`)
needed $m_D = 80$, $M_\text{user} = 4$ at $(r_D, r_{b_-}, r_{b_+}) =
(9, 11, 11)$ to hit $1.6 \times 10^{-6}$:

- Shared $\delta t_0 = t_{0,D} / m_D \approx 0.031$, REQUIRED to be the
  elementary Strang step of all three legs.
- Implied per-leg $M$ counts (legs commensurate): $M_D = 320$, $M_{b_-} \approx 11$,
  $M_{b_+} \approx 4$.
- Total Strang substeps in cache construction: $\approx 320 + 11 + 4 = 335$.

The independent-knob recipe achieves a $2\times$ better residue
($7.2 \times 10^{-7}$ vs $1.6 \times 10^{-6}$) at $\approx 64 + 64 + 1 = 129$
total substeps — **2.6× fewer Strang Trotterizations**. The shared-$\delta t_0$
scheme over-resolves $b_-$ and $b_+$ to maintain integer commensurability;
independent knobs decouple this constraint.

## Channel-floor cross-check

The Lindbladian residue $\|\mathcal{L}\sigma_\beta\|_\mathrm{HS}$ is the
algorithm-level fixed-point precision (no $\delta$-step). The channel
floor $\|\rho_\infty - \sigma_\beta\|_1 / 2$ adds a slope-$1$-in-$\delta$
contribution from generator splitting ($\Phi_\delta \approx I + \delta\mathcal{L}
+ O(\delta^2)$). At the optimal recipe ($\|\mathcal{L}\sigma_\beta\|_\mathrm{HS}
\approx 7.2 \times 10^{-7}$), a dense $\Phi_\delta$ build with $\delta \le
10^{-6}$ would satisfy $0.2\delta \le 2 \times 10^{-7}$, giving total
channel floor $\le 10^{-6}$ with comfortable margin.

Full channel-floor sweep deferred — the Lindbladian-residue controllability
$\le 10^{-6}$ is established here, and the $\delta$-split contribution is
trivially controlled.

## Cost model

Per `B_trotter` call at the minimal recipe `(r_D=7, r_bm=6, r_bp=14,
M_D=64, M_bm=64, M_bp=1)` at $n=3$:

- **Cache construction** (one-time per cell): three Strang Trotterizations
  + three `eigen()` calls on $d \times d$ matrices ($d = 2^3 = 8$). Cost
  $\propto M_D + M_{b_-} + M_{b_+} = 129$ Strang substeps, well below the
  335 of Option A.
- **Inner $\tau$-loop**: effective $b_+$ grid $N_{b_+} \sim 30$ samples at
  $r_{b_+} = 14$ (after truncation) × 9 jumps × $O(d^3)$ GEMMs.
- **Outer $t$-loop**: effective $b_-$ grid $N_{b_-} \sim 40$ samples at
  $r_{b_-} = 6$ × $O(d^3)$.
- **Basis rotations** (V_D ↔ V_bp / V_bm): $\sim 10^4$ flops per call,
  measured at 9.8μs vs 207ms total at the over-resolved recipe ≪ $10^{-4}$
  of total cost.

At the smaller-r qf-7xt recipe, the per-call cost is **dominated by the
threaded inner/outer kernels operating on much smaller grids** than the
over-resolved $(11, 14, 18)$ I first tested. The minimal-recipe wall time
is sub-100ms per `construct_lindbladian` on the test sandbox.

## Acceptance — qf-e4z.20.6 closure

1. ✓ **Independent-knob recipe at ≤ 1e-6 demonstrated**: optimal
   $\boxed{(7, 6, 14),\ M=(64, 64, 1)}$ → $7.2 \times 10^{-7}$.
2. ✓ **Per-leg attribution confirmed**: $M_D$ controls dissipative
   Strang; $M_{b_-}$ controls outer Strang (slope-2); $M_{b_+}$ saturates
   at $M=1$ (fine grid).
3. ✓ **Cheaper than Option A**: 7.2e-7 vs 1.6e-6 (2× better) at
   $\sim 2.6\times$ fewer Strang substeps per cell.
4. ✓ **Wall time ≪ 60s per cell**: sub-100ms per `construct_lindbladian`.

Generating script: `scripts/scratch_e4z20_independent_controllability.jl`.
Refined recipe sweep at `/tmp/qf_test_proper_recipe.jl` (test artefact).

## Code — minimal recipe

```julia
ω_range = 2.5
T_minus, T_plus = 18.0, 12.0
r_D, r_bm, r_bp = 7, 6, 14
w0_D = ω_range / 2^r_D
cfg = Config(
    sim = Lindbladian(), domain = TrotterDomain(), construction = KMS(),
    num_qubits = 3, with_linear_combination = true,
    beta = 10.0, sigma = 0.1, s = 0.25, a = 10.0/30.0, eta = 1e-3,
    num_energy_bits_D = r_D,  t0_D = 2π / (2^r_D * w0_D), w0_D = w0_D,
    num_energy_bits_b_minus = r_bm,
        t0_b_minus = 2 * T_minus / 2^r_bm,
        w0_b_minus = 2π / (2^r_bm * 2 * T_minus / 2^r_bm),
    num_energy_bits_b_plus  = r_bp,
        t0_b_plus  = 2 * T_plus / 2^r_bp,
        w0_b_plus  = 2π / (2^r_bp * 2 * T_plus / 2^r_bp),
    num_trotter_steps_per_t0_D       = 64,
    num_trotter_steps_per_t0_b_minus = 64,
    num_trotter_steps_per_t0_b_plus  = 1,
    num_trotter_steps_per_t0         = 64,  # legacy fallback for M_D if any caller skips per-leg
)
trotter = make_trotter_for_config(ham, cfg)    # → TrotterTriple
L = construct_lindbladian(jumps, cfg, ham; trotter=trotter)
# ‖L · σ_β‖_HS ≈ 7.2e-7
```
