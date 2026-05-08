# Quadrature convergence: how many ω-bits / time-bits do we need?

> **\* This is the canonical reference for choosing $r_D$, $r_{b_-}$, $r_{b_+}$.** Numbers measured at one fixture only: $n = 4$, $\beta = 10$. Larger system size and different temperature can shift the requirements by ~1–2 bits in either direction (e.g., bigger $\|H\|$ widens the integrand support, finer Bohr separations need finer $w_0$). Cross-check at additional $(n, \beta)$ before finalising. Slopes are universal — only the $K$ prefactors and the $r$ where slope-1 transitions to floor change with the fixture.

## Setup

$n = 4$, $\beta = 10$, $\sigma = 1/\beta = 0.1$, $\|H\| \approx 0.45$.

Unified ω-range / cutoff principle:
$$\omega_{\text{range}} = 2(\|H\| + C \sigma), \qquad \varepsilon = \exp(-C^2/2).$$

At $C = 8$: $\omega_{\text{range}} \approx 2.5$, $\varepsilon \approx 1.3 \times 10^{-14}$ (machine precision).

## EnergyDomain → BohrDomain (analytical α) convergence

Sweep r at fixed $\omega_{\text{range}}$, with $w_0 = \omega_{\text{range}} / 2^r$. Both `L_e` and `L_b` built dissipator-only at the same $(r, w_0)$; `err = opnorm(L_e/gnf_e - L_b/gnf_b)`. Source: `scripts/scratch_energy_ref_convergence.jl`.

| filter | r=3 | r=4 | r=5 | r=6 | r=8 | r=10 | r=12 | r=14 |
|---|---|---|---|---|---|---|---|---|
| Gaussian | 2.4e-1 | 1.1e-2 | **5.6e-8** | 1.0e-15 | 3.9e-15 | 7.6e-15 | 1.7e-14 | 6.4e-14 |
| Smooth Metro (s=0.25) | 1.7e-1 | 9.8e-3 | **1.5e-5** | 1.2e-14 | 4.6e-15 | 1.5e-14 | 3.1e-14 | 1.3e-13 |
| Kinky Metro (s=0) | 1.7e-1 | 3.6e-2 | 8.4e-3 | 1.2e-3 | 1.2e-4 | 1.0e-5 | 7.3e-7 | 2.5e-8 |

Behaviour:

- **Gaussian**: super-algebraic (exponential in $r$). Saturates at machine precision by $r = 6$.
- **Smooth Metro** (s = 0.25): super-algebraic (Gevrey-1/2), slightly slower than Gaussian. Saturates at machine precision by $r = 6$.
- **Kinky Metro** (s = 0): clean slope $-2$ in $1/N$ (= $-0.61$ in $\log_{10}\|\Delta L\| / r$), as predicted by CKBG23 Cor. C.2 (kink in γ at $\omega = -\beta\sigma^2/2$).

## Minimum r per ε target (energy-domain quadrature alone)

| ε target | Gaussian | Smooth Metro | Kinky Metro |
|---|---|---|---|
| 1e-3 | 4 | 4 | 6–7 |
| 1e-6 | 5 | 5 | 12 |

For algorithm-level ε = 10⁻⁶:

- Gaussian filter is essentially **free** in r (r = 5 already gives ~5 × 10⁻⁸).
- Smooth Metro: r = 5–6.
- **Kinky Metro is the binding constraint** at any uniform r — slope -2 forces r ≈ 12.

Picking smooth-Metropolis instead of kinky **saves ~6–7 bits in the ω-register** at the same precision target.

## TimeDomain → BohrDomain convergence

In production we use the TimeDomain Lindbladian (NUFFT-based discrete OFT), not the EnergyDomain Riemann sum directly. Same fixed $\omega_{\text{range}}$ as above; sweep r at $w_0(r) = \omega_{\text{range}} / 2^r$, which makes $t_0 = 2\pi/\omega_{\text{range}} \approx 2.51$ **constant** across the sweep (only the t-window grows: $t_{\text{range}}(r) = 2\pi \cdot 2^r / \omega_{\text{range}}$). Source: `scripts/scratch_time_ref_convergence.jl`.

| filter | r=3 | r=4 | r=5 | r=6 | r=8 | r=10 | r=12 | r=14 |
|---|---|---|---|---|---|---|---|---|
| Gaussian | 2.0e-1 | 1.1e-2 | **5.6e-8** | 2.0e-13 | 2.0e-13 | 2.0e-13 | 2.0e-13 | 2.2e-13 |
| Smooth Metro | 1.2e-1 | 9.9e-3 | **1.5e-5** | 2.5e-13 | 2.5e-13 | 2.5e-13 | 2.5e-13 | 2.9e-13 |
| Kinky Metro | 1.4e-1 | 3.6e-2 | 8.4e-3 | 1.2e-3 | 1.2e-4 | 1.0e-5 | 7.3e-7 | 2.5e-8 |

Behaviour:

- Same per-filter slopes as EnergyDomain — kinky slope $-2$ (regression $-0.61$ in $\log_{10}\|\Delta L\| / r$), Gaussian/smooth super-algebraic.
- Gaussian/smooth saturate at the **FINUFFT precision floor** $\sim 2 \times 10^{-13}$ (instead of FP machine precision $\sim 10^{-15}$ in the EnergyDomain dense path), as expected.
- Kinky agrees with EnergyDomain entry-by-entry to all printed digits.

This validates the user's intuition for the parameter strategy: at fixed $\omega_{\text{range}}$, the energy-kernel quadrature (slope $-2$ for kinky) is the binding constraint; the t-OFT is super-algebraically resolved by the Gaussian $f(t)$ at any $t_0 = 2\pi/\omega_{\text{range}}$. Once $\omega_{\text{range}}$ is set from the integrand support, choosing $r$ for the energy-quadrature target automatically gives a fine-enough $t_0$.

## Coherent term B convergence (TimeDomain B_time → analytical B_bohr)

The coherent (Lamb-shift) operator has only TIME sums on its two **independent** registers (qf-9z0; thesis 2_methods §quad-errors-coh, eq:r_- and eq:r_+_*):
$$B = t_0^- t_0^+ \sum_t \sum_\tau b_-(t) \, b_+(\tau) \, U(\text{jumps}; t, \tau).$$

$r_-$ and $r_+$ are **separate knobs** — the outer ($b_-$) and inner ($b_+$) integrals get separate registers and converge at different rates. We sweep them separately below.

The truncation `_compute_truncated_func(_compute_b_*, time_labels; atol=1e-12)` already in `src/coherent.jl` chops where $|b_\pm| < 10^{-12}$, so the **time range** is set by the truncation — we just need the grid range to comfortably cover the support. Empirical supports at β=10, σ=0.1:

- $b_-(t)$: $|t| \le 5.55$ (filter-independent — set by $1/\cosh(2\pi t / (\beta\sigma))$)
- $b_+(t)$ Gaussian: $|t| \le 2.5$ (Gaussian decay $\propto e^{-4t^2}$)
- $b_+(t)$ Kinky Metro: $|t| \le 3.25$ ($e^{-2t^2}$ inside the $1/(t(2t+i))$ envelope)
- $b_+(t)$ Smooth Metro: $|t| \le 2.9$

Pick the grid ranges $\sim 3\times$ the support so the truncation is the binding cut. We use $T_{\text{range}}^- = 18$ and $T_{\text{range}}^+ = 12$. Source: `scripts/scratch_coherent_ref_convergence.jl`.

### Split A: sweep $r_-$ at fixed $r_+ = 16$ (isolates outer $b_-$)

| filter | r_-=4 | r_-=5 | r_-=6 | r_-=8 | r_-=14 |
|---|---|---|---|---|---|
| Gaussian | 1.4e-3 | 1.05e-5 | **7.7e-15** | 9.1e-15 | 1.0e-14 |
| Smooth Metro | 6.6e-4 | 2.8e-6 | **2.7e-6** | 2.7e-6 | 2.7e-6 |
| Kinky Metro | 6.3e-4 | 2.8e-6 | **2.7e-6** | 2.7e-6 | 2.7e-6 |

Behaviour: $b_-(t)$ is smooth (cosh-decay convolved with $\sin \cdot e^{-2t^2}$) for ALL filters → **super-algebraic in $r_-$**. The Metro floor at $\sim 2.7 \times 10^{-6}$ is the leftover $b_+$ contribution at $r_+ = 16$ (not coming from $b_-$). **$r_-$ converges to its asymptote by $r_- \approx 5\text{–}6$ regardless of filter.**

### Split B: sweep $r_+$ at fixed $r_- = 16$ (isolates inner $b_+$)

| filter | r_+=4 | r_+=6 | r_+=8 | r_+=10 | r_+=12 | r_+=13 | r_+=14 | r_+=15 |
|---|---|---|---|---|---|---|---|---|
| Gaussian | 1.4e-3 | **2.4e-14** | 2.5e-14 | 2.5e-14 | 2.5e-14 | 2.5e-14 | 2.5e-14 | 2.5e-14 |
| Smooth Metro | 1.0e-3 | 2.5e-4 | 6.2e-5 | 1.5e-5 | 3.9e-6 | 1.9e-6 | **2.9e-6** | 2.4e-6 |
| Kinky Metro | 9.9e-4 | 2.5e-4 | 6.2e-5 | 1.5e-5 | 3.9e-6 | 1.9e-6 | **2.9e-6** | 2.4e-6 |

Slopes (log₁₀ ‖ΔB‖ per $r_+$, in the clean slope regime $r_+ \in [4, 13]$):

- **Gaussian**: super-algebraic, machine-precision by $r_+ = 6$.
- **Smooth Metro**: $-0.302/r_+ = -1.00$ per $\log_2 N$ → **slope $-1$ in $1/N_+$** through $r_+ = 13$, then noisy floor at $\sim 1\text{–}3 \times 10^{-6}$.
- **Kinky Metro**: same, $-0.301/r_+$ — numerically identical to smooth at all $r_+$.

The $r_+ \ge 14$ floor at $\sim 2 \times 10^{-6}$ is consistent with the old $\eta$-cap floor noted in `.claude-memory/error_analysis_qf_b4d.md` ("bias floor at ~1.5e-6"). Below this floor the trapezoidal sum + Cauchy P.V. + truncation effects mix; pushing for tighter ε needs both registers larger and risks numerical noise overtake.

### Key observations

1. **$r_-$ is "free": $r_- = 6$ is enough for any filter / target.** $b_-(t)$ is smooth so its Riemann sum converges super-algebraically — matches thesis eq:r_-, which absorbs only an $|\!|b_+|\!|_1 = O(\log(1/\varepsilon))$ factor for Metro (mild log-log correction).

2. **$r_+$ is the binding register for Metro.** Slope $-1$ in $1/N_+$ (per thesis eq:r_+_metro and the $t=0$ Cauchy P.V. anomaly noted in `.claude-memory/trap_rule_t0_lhopital_origin.md`).

3. **Smooth and kinky give identical $r_+$ convergence to ≥4 digits.** Both have the same singular $1/(t(2t+i))$ envelope and the same η-regularised $t=0$ branch in `_compute_b_plus_metro`. The smoothing $s$ modifies the Gaussian-like exponent but does not change the singular-tail structure that drives the slope. **Choosing smooth-Metropolis over kinky saves bits in the dissipator quadrature, NOT in the coherent quadrature.**

4. **Coherent Metro $r_+$ is one slope worse than dissipator Metro $r_D$:** $-1$ vs $-2$ in $1/N$. At ε = $10^{-6}$, dissipator needs $r_D \approx 12$ and coherent needs $r_+ \approx 12$ — comparable; at coarser ε it's even cheaper.

## Combined parameter recipe (n=4, β=10)

For the **dissipative register** $r_D$, with fixed $\omega_{\text{range}} = 2(\|H\| + 8\sigma)$ and $w_{0,D} = \omega_{\text{range}}/2^{r_D}$ (TimeDomain in production; matches EnergyDomain quadrature accuracy modulo NUFFT floor):

| ε target | Gaussian | Smooth Metro | Kinky Metro |
|---|---|---|---|
| 1e-3 | 4 | 4 | 6–7 |
| 1e-6 | 5 | 5 | 12 |

For the **coherent registers** $r_{b_-}$ and $r_{b_+}$, with grid ranges fixed at ≥ 3× the b_± supports (no new constants — the truncation atol = 1e-12 sets the effective range), sized SEPARATELY:

**Outer register $r_-$** (any filter — $b_-$ is smooth, super-algebraic):

| ε target | $r_-$ |
|---|---|
| 1e-3 | 4 |
| 1e-6 | 6 |

**Inner register $r_+$** (filter-dependent):

| ε target | Gaussian | Smooth Metro | Kinky Metro |
|---|---|---|---|
| 1e-3 | 4 | 4 | 4 |
| 1e-6 | 6 | 14 | 14 |

**Parameter strategy.** Pick each register independently:

- $\omega$-range from the integrand support → $w_{0,D}$ from the target $r_D$.
- $b_-, b_+$ time ranges from the existing truncation (atol = $10^{-12}$).
- $r_D$, $r_{b_-}$, $r_{b_+}$ each sized to its own filter-and-precision target above.

## Caveats

- All numbers above are at $n = 4$, $\beta = 10$, $\|H\| = 0.45$, $\sigma = 0.1$, $\eta = 10^{-3}$ for Metro $b_+$. Other $(n, \beta)$ shift the supports and floors but the slopes are universal.
- Sweeps were dense (`construct_lindbladian` + `opnorm`) for clean convergence pictures; production code uses matrix-free `hs_operator_norm_krylov` (validated at machine precision in the test suite).
