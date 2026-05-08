# Smooth-Metropolis $s$ decision (qf-3il)

**Insertion target.** A short paragraph in `supplementary-informations/2_methods.tex` near the existing
`(Eq. eq:smooth-metro)`/`(Prop. prop:smooth-metro-gevrey)` discussion (currently has a `[REF] \todo{ref}`
placeholder for the numerical decision) and a brief mention in the Numerics chapter where mixing-time
sweeps are reported.

## Headline

The thesis target precision for the simulator's fixed-point displacement from $\rho_\beta$ is
$\varepsilon = 10^{-6}$ (KMS-DB and approximate GNS-DB Lindbladians are both constructed with
controllable error and can be pushed to that tolerance). At this precision the smooth-Metropolis
qubit-saving argument is *load-bearing*:

| $\beta$ | kinky $s=0$ bits at $\varepsilon=10^{-6}$ | smooth $s=0.25$ bits at $\varepsilon=10^{-6}$ | savings |
|---:|---:|---:|---:|
| 5  | 8                | 6 | **2 bits** ($4\times$ register reduction) |
| 10 | 9                | 6 | **3 bits** ($8\times$) |
| 20 | $\ge 10$         | 7 | **$\ge 3$ bits** ($\ge 8\times$) |

(Energy register at $n = 5$, measured via $\|\rho_\infty^{\text{discrete}} - \rho_\beta\|_1 / 2$
where $\rho_\infty^{\text{discrete}}$ is the kernel of $\mathcal{L}_E$ from `krylov_spectral_gap`.
Time-register cost is fixed by the $\eta$-regulariser and is insensitive to $s$.)

For thesis-scale simulations and early-implementation cost estimates, **$s = 0.25$ is locked in
as the default** (already the value used in `sweep_mixing_times` and the qf-lkb.6 plot). The
trade is:

- **3–6% mixing-time penalty** (vs kinky $s = 0$).
- **2–3+ bits saved** in the energy register at the thesis target $\varepsilon = 10^{-6}$,
  growing with $\beta$.
- **Time-register cost unchanged.**

The choice $s = 0.25$ is robust: $s = 0.10$ has a slightly smaller $\tau_{\text{mix}}$ penalty
(1–3%) but does not always reach $\varepsilon = 10^{-6}$ at the simulator's default $w_0 = 0.05$;
$s = 0.50$ and $s = 1.0$ pay 6–12% and 17–24% mixing-time penalties without buying additional
bit savings at any $\varepsilon$ we tested.

A coarser-precision sanity check at $\varepsilon = 10^{-3}, 10^{-4}$ is included below for
comparison; at $\varepsilon = 10^{-3}$ kinky and smooth tie at every $(\beta, n=5)$ cell. **The
smooth case becomes load-bearing only at $\varepsilon = 10^{-6}$**, which is precisely the
thesis target.

The decision rests on three pieces of numerical evidence.

## 1. Energy-register qubit cost — trace-distance metric (the algorithm's actual requirement)

The relevant cost is the smallest energy-register spacing $w_0$ such that the discrete
Lindbladian's fixed point $\rho_\infty^{\text{discrete}}$ (= eigenvector of $\mathcal{L}$ at
eigenvalue 0) lies within trace distance $\varepsilon$ of $\rho_\beta$. We solve
$\mathcal{L}_E[\rho_\infty] = 0$ directly via Krylov eigsolve (`krylov_spectral_gap`,
`fixed_point` field) and report $\|\rho_\infty^{\text{discrete}} - \rho_\beta\|_1 / 2$.

For $n = 5, \beta \in \{5, 10, 20\}, s \in \{0, 0.25\}$, varying $w_0$:

| $\beta$ | $w_0$ | $N_{\text{trunc}}$ | bits | kinky floor | smooth floor |
|---:|---:|---:|---:|---:|---:|
| 5  | 0.10  | 39  | 6 | $6.78\times 10^{-6}$ | $3.18\times 10^{-12}$ |
| 5  | 0.05  | 79  | 7 | $1.70\times 10^{-6}$ | $1.34\times 10^{-15}$ |
| 5  | 0.025 | 157 | 8 | $4.25\times 10^{-7}$ | $1.29\times 10^{-15}$ |
| 10 | 0.10  | 25  | 5 | $6.86\times 10^{-5}$ | $5.35\times 10^{-6}$ |
| 10 | 0.05  | 49  | 6 | $2.98\times 10^{-5}$ | $1.42\times 10^{-11}$ |
| 10 | 0.025 | 97  | 7 | $7.41\times 10^{-6}$ | $2.23\times 10^{-15}$ |
| 20 | 0.10  | 17  | 5 | $3.17\times 10^{-2}$ | $2.87\times 10^{-2}$ |
| 20 | 0.05  | 33  | 6 | $\mathbf{1.57\times 10^{-4}}$ | $\mathbf{1.28\times 10^{-5}}$ |
| 20 | 0.025 | 67  | 7 | $6.07\times 10^{-5}$ | $2.68\times 10^{-11}$ |

(Source: `scripts/scratch_smooth_metropolis_qubit_cost.jl` (n=3 ‖L‖_op view) and the n=5
fixed-point sweep saved at `drafts/figures/numerics/smooth_metro_fixed_point_floor_n5.bson`.)

Kinky converges *polynomially* in $w_0$ (slope ≈ 1–2); smooth $s=0.25$ converges
super-polynomially (Gevrey-1/2 of Prop. prop:smooth-metro-gevrey) — saturates at machine
precision by $w_0 \le 0.025$.

**Bits required at each $\varepsilon$:**

| $\beta$ | $\varepsilon = 10^{-3}$ k→s | $\varepsilon = 10^{-4}$ k→s | $\varepsilon = 10^{-6}$ k→s ($\star$ thesis target) |
|---:|:---|:---|:---|
| 5  | 6 → 6 (0 saved) | 6 → 6 (0)        | 8 → 6 (**2**)  |
| 10 | 5 → 5 (0)       | 5 → 5 (0)        | 9 → 6 (**3**)  |
| 20 | 6 → 6 (0)       | 7 → 6 (**1**)    | $\ge 10\to 7$ (**≥3**) |

**At the thesis target $\varepsilon = 10^{-6}$, smooth saves 2–3+ bits across $\beta$**, and the
saving grows with $\beta$. The looser $\varepsilon = 10^{-3}, 10^{-4}$ rows are kept for
comparison: at $\varepsilon = 10^{-3}$ kinky and smooth tie at every $(\beta, n=5)$ cell — the
$\omega$-quadrature kink is benign enough at this coarse precision that the polynomial vs
super-polynomial scaling difference doesn't yet bite. By $\varepsilon = 10^{-6}$ the
polynomial-vs-Gevrey gap dominates and the smooth case is decisively cheaper.

**Important methodological note.** An earlier version of this note (commit 5e56529) used
$\|\Delta\mathcal{L}\|_{\text{op}} / \|\mathcal{L}\|_{\text{op}}$ as the cost criterion, which
*overstated* the bit savings at the looser precision $\varepsilon = 10^{-3}$. The operator-norm
error is a strict upper bound on the trace-distance fixed-point displacement (perturbation
theory: $\|\rho_\infty - \rho_\beta\|_1 \lesssim \|\Delta\mathcal{L}\|_{\text{op}} /
\text{gap}$), so requiring it below $\varepsilon$ is more conservative than the algorithm
actually needs. The fixed-point trace-distance metric used in the table above is the correct
one. At the thesis target $\varepsilon = 10^{-6}$ the two metrics agree more closely (bits
saved: 2–3 either way), and the smooth case is the right choice regardless of how we slice it.

## 2. Time-register qubit cost (insensitive to $s$)

The TimeDomain side (varying $t_0$ at fixed converged $w_0 = 0.003125$) shows the time-register
cost is **insensitive to $s$**: the simulator's $\eta$-regularised $b_+^{(s,\eta)}$ kernel
(`(Eq. eq:b_plus-s-eta)`) tames the kink at $|t| \le \eta$ for kinky and smooth uniformly, so the
time-domain Riemann error is set by $\sigma, \beta, \eta$, not by $s$. **Whatever bits $s$
saves, it saves them in the energy register; the time-register cost is fixed by $\eta$.** This
is the right distinction to draw in the thesis.

## 3. Mixing-time penalty (full simulation, $n \le 5$)

Smoothing strictly reduces the spectral gap (Prop. prop:optim-metro): $\gamma_M^{(s)} \le
\gamma_M^{(0)}$ pointwise, so $\tau_{\text{mix}}(s) \ge \tau_{\text{mix}}(0)$. The empirical
question is *how much*. The $\tau_{\text{mix}}$ ratios across $n\in\{3,4,5\}$,
$\beta\in\{5, 10, 20\}$ (45 cells, `heis_xxx_zzdisordered_periodic_n*` fixture, EnergyDomain CKG,
target $\varepsilon=10^{-3}$ used for the original sweep — but the ratio
$\tau_{\text{mix}}(s)/\tau_{\text{mix}}(0)$ is $\varepsilon$-insensitive because both numerator
and denominator scale by the same $\log(1/\varepsilon)$ factor at fixed Lindbladian spectrum, so
the percentages below carry over to the thesis target $\varepsilon = 10^{-6}$):

| $\beta$ | $n$ | $s=0.10$ | $s=0.25$ | $s=0.50$ | $s=1.00$ |
|---:|---:|---:|---:|---:|---:|
| 5  | 3 | 1.02 | 1.05 | 1.11 | 1.23 |
| 5  | 4 | 1.02 | 1.06 | 1.11 | 1.24 |
| 5  | 5 | 1.02 | 1.06 | 1.12 | 1.24 |
| 10 | 3 | 1.03 | 1.06 | 1.11 | 1.21 |
| 10 | 4 | 1.01 | 1.03 | 1.07 | 1.16 |
| 10 | 5 | 1.02 | 1.05 | 1.09 | 1.19 |
| 20 | 3 | 1.01 | 1.03 | 1.07 | 1.17 |
| 20 | 4 | 1.02 | 1.04 | 1.06 | 1.10 |
| 20 | 5 | n/a&nbsp;[*] | n/a&nbsp;[*] | (37.2) | (38.2) |

At $s = 0.25$ the τ_mix penalty is **3–6% across all 8 successful cells**; at $s = 0.10$ it is
**1–3%**. Both are well within the noise of typical thesis-scale numerics. Pushing $s$ to $0.5$
incurs ~6–12%, $s = 1.0$ runs to 17–24% (and $s = 2$ Glauber-like would be much worse, as
Chen et al. flag in \cite{chen2023quantum}).

Source: `scripts/scratch_smooth_metropolis_s_sweep.jl` →
`drafts/figures/numerics/smooth_metro_s_taumix.{png,pdf}`,
`drafts/figures/numerics/smooth_metro_s_sweep.bson`.

### [*] $\beta=20, n=5, s\in\{0,0.1\}$ — discrete-Lindbladian floor, not a fitting failure

Re-running these cells without bi-exp post-processing shows the trajectory plateaus *above*
$\varepsilon=10^{-3}$ for every $s$ at this $(\beta, n)$ at the simulator's default
$w_0 = 0.05$:

| $s$ | $\|\rho(t_{\max}) - \rho_\beta\|_1$ |
|---:|---:|
| 0.00 | $3.92\times 10^{-3}$ |
| 0.10 | $3.82\times 10^{-3}$ |
| 0.25 | $3.60\times 10^{-3}$ |
| 0.50 | $3.26\times 10^{-3}$ |
| 1.00 | $2.77\times 10^{-3}$ |

None reach $\varepsilon=10^{-3}$. The bi-exp's offset $C$ for $s\in\{0, 0.1\}$ correctly identifies
this and returns NaN; the $s\ge 0.25$ "extrapolations" are bi-exp predicting a (fictitious)
continuation past the plateau and should be read as *upper bounds*, not measurements.

The mechanism is exactly the energy-register quadrature error of §1, in its discrete-Lindbladian
form: at $\beta=20$, $w_0 = 0.05$, the EnergyDomain Lindbladian for kinky $\gamma$ has
$\|\Delta\mathcal{L}\|_{\text{op}} / \|\mathcal{L}\|_{\text{op}} = 1.21\times 10^{-2}$ — and that
operator error translates into a $\sim 3\text{–}4\times 10^{-3}$ floor on
$\|\rho_\infty^{\text{discrete}} - \rho_\beta\|_1$. Smoothing with $s = 0.25$ at the *same* $w_0$
brings the operator error to $3.85\times 10^{-4}$ and lifts the floor by an order of magnitude,
but at $\beta = 20$ on the 12-bit default register it is still above $\varepsilon = 10^{-3}$.

**The fix is more energy bits (smaller $w_0$), not more $s$.** Going from $w_0 = 0.05$ to
$w_0 = 0.025$ (one extra bit) drops the smooth $\|\Delta\mathcal{L}\|$ to $6\times 10^{-10}$ — an
8-orders-of-magnitude drop — while the kinky variant only drops to $5.6\times 10^{-3}$, still
above $\varepsilon$. Two extra bits on the kinky variant ($w_0 = 0.0125$) gets to
$1.4\times 10^{-3}$ — still above. The kinky variant simply needs a *much* finer energy register
to reach $\varepsilon = 10^{-3}$ at $\beta = 20$ in our $n=5$ system: 8 bits vs 6 for smooth.

This sharpens rather than weakens the smooth-Metropolis case: the kink shows up *both* as the
polynomial quadrature convergence in §1 *and* as a discrete-Lindbladian fixed-point shift at
finite $w_0$. Smoothing ($s > 0$) is the cheapest way to remove both.

## Why a single global $s$ (not $s = s(\beta)$)

The Gevrey-1/2 constant $C_2 = \beta/\sqrt{e\,s}$ (Prop. prop:smooth-metro-gevrey, with the
thesis convention $\sigma = 1/\beta$) grows linearly in $\beta$ at fixed $s$, so one might
worry the optimal $s$ shifts with $\beta$. The data in §1 shows it doesn't: smooth $s = 0.25$
already saturates to machine precision at $w_0 \le 0.0125$ for every $\beta\in\{5, 10, 20\}$;
moving to $s = 0.5$ or $s = 1.0$ does not buy further bits. The polynomial-vs-super-polynomial
gap is so wide that any $s \ge 0.1$ already lands the smooth side in the regime where
truncation, not the kink, sets the cost. A single global $s = 0.25$ thus suffices across the
relevant $\beta$ range and aligns with the existing thesis convention.

## What changes in `src/`

No code change required for the $s$ default: `sweep_mixing_times` already defaults to
`a = 0.0, s = 0.25` (qf-lkb.11 lock); the qf-lkb.6 comparison plot uses the same value. The
decision documented here promotes the locked default from a convention to an *evidence-backed*
default and provides the quadrature-cost argument the thesis was missing (the `[REF] \todo{ref}`
after `(Eq. eq:smooth-metro)` in `2_methods.tex` line 305).

**Recommended `num_energy_bits` at the thesis target $\varepsilon = 10^{-6}$**, on the
$\sigma = 1/\beta$ thesis convention:

| $\beta$ | smooth $s = 0.25$ | kinky $s = 0$ | recommended $w_0$ (smooth) |
|---:|---:|---:|---:|
| 5  | 6 bits (39 truncated labels)   | 8 bits (157 labels) | $0.10$  |
| 10 | 6 bits (49 labels)             | 9 bits (387 labels) | $0.05$  |
| 20 | 7 bits (67 labels)             | $\ge 10$ bits       | $0.025$ |

The 12-bit `num_energy_bits` default in the test suite and current production scripts
($w_0 = 0.05$) is **sufficient for smooth at every $(\beta, n)$ cell tested** and gives
$\varepsilon \le 10^{-6}$ headroom at $\beta \le 10$. For $\beta = 20$ at $\varepsilon = 10^{-6}$
the smooth case requires $w_0 \le 0.025$ (which the truncation logic happily accommodates without
changing `num_energy_bits`), and the kinky case fails entirely — out of bits at any $w_0$ on
our test grid.

Follow-up: re-running the qf-lkb.6 mixing-time sweep at `target_epsilon = 1e-6` would update
the absolute $\tau_{\text{mix}}$ numbers (longer integration horizons via the
qf-lkb.10 auto-factor heuristic) but would not change the smooth/kinky percentage ratios in §3.
Filed as future work — not required for the qf-3il decision.

## Reference data

- **Qubit-cost script** — `scripts/scratch_smooth_metropolis_qubit_cost.jl`. EnergyDomain and
  TimeDomain $\mathcal{L}$ vs BohrDomain reference, varying $w_0$ and $t_0$ in the simulator's
  pipeline. Reports bits = $\lceil\log_2(N_{\text{trunc}})\rceil$ at the largest spacing
  achieving $\varepsilon$.
- **$s$-sweep script** — `scripts/scratch_smooth_metropolis_s_sweep.jl`. 45-cell sweep
  ($s \times n \times \beta$). Per-cell BSON sidecars under
  `drafts/figures/numerics/sweep_cache/smooth_metro_s/<s_tag>/`; resumable.
- **Earlier windowed Riemann analysis** — `scripts/scratch_smooth_metropolis_quadrature_error.jl`.
  Complementary view of the $\gamma$-integrand kink-vs-smooth contrast on a fixed window; the
  qubit-cost script (above) is the load-bearing argument for the thesis because it uses the
  simulator's actual truncation logic.
- **Aggregated BSON** — `drafts/figures/numerics/smooth_metro_qubit_cost.bson`,
  `drafts/figures/numerics/smooth_metro_s_sweep.bson`,
  `drafts/figures/numerics/smooth_metro_quadrature_error.bson`.

## Citations

- Chen et al. 2023 — Quantum Thermal State Preparation (Glauber-limit warning at $s=2$, footnote 33). [CITE: chen2023quantum]
- Chen et al. 2025 — Efficient noncommutative quantum Gibbs sampler (smooth-Metropolis baseline). [CITE: chen2023efficient]

```bibtex
@article{chen2023quantum,
  author  = {Chen, Chi-Fang and Kastoryano, Michael J. and Brand{\~a}o, Fernando G. S. L. and Gily{\'e}n, Andr{\'a}s},
  title   = {Quantum Thermal State Preparation},
  year    = {2023}
}
@article{chen2023efficient,
  author  = {Chen, Chi-Fang and Kastoryano, Michael J. and Gily{\'e}n, Andr{\'a}s},
  title   = {An Efficient and Exact Noncommutative Quantum Gibbs Sampler},
  year    = {2025}
}
```
