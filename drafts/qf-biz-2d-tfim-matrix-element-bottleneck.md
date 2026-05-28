# 2D TFIM ORDERED gap closing: phase check + matrix-element diagnostic (qf-biz)

**Status**: complete (n=6 β-sweep + n=4,6 matrix elements + dense-L eigenmode at n=4 + kinky-Metro non-monotonicity check). Verified by `code-verifier` (bit-exact agreement with qf-1jj) and `physics-checker` (claims 1–4 OK; claim 5 amended below).
**Driver**: `scripts/scratch_qf_biz_phase_and_matrix_elements.jl` (Checks 1–2, wall 60.9 s) + `scripts/scratch_qf_biz_followup.jl` (Follow-ups A–B, wall 34.3 s).
**Sidecars**: `scripts/output/qf_biz_phase_and_matrix_elements/check{1,2}_*.bson` + `followup_a_b.bson`.
**Reviews**: `drafts/qf-biz-code-verifier-report.md`.
**Parent**: qf-1jj (`drafts/2d-tfim-ordered-vs-disordered.md`).

> **Headline.** The ORDERED β_phys=2 cell is genuinely deep in the ordered phase (⟨M_z²⟩≈0.95, U_4≈0.66, eff_rank=2.00 across the entire β_phys ≥ 1 range — **distinct** from the 1D-Heisenberg-even-n "p_0 ≈ 99 % pure-GS" artefact). **But the within-doublet matrix-element bottleneck I proposed earlier is wrong**: $M^2 = \sum_a |\langle \psi_2 | A_a | \psi_1\rangle|^2$ is **constant** at 0.311 from n=4 to n=6 (sum rule: 1/(3n)·∑_i ⟨ψ_2|σ_z^(i)|ψ_1⟩² ≈ n/(3n) = 1/3). Follow-up A (dense L eigendecomp at n=4) pinpoints why: **the slow mode is not within-doublet at all — 89 % of $\|R_2\|^2$ lives in the doublet × bulk off-diagonal block** (largest entries are $|\psi_{1,2}\rangle\langle\psi_8|$-type coherences), only 2.8 % in the doublet × doublet block. The mode is a Z₂-mixed-parity, M_z-orthogonal coherence between the doublet and a specific bulk state, decohering at a rate involving energy-increasing transitions out of |ψ_1,2⟩ to bulk. The 2D-TFIM ORDERED gap closing IS thermal-order-driven, but the kinetic bottleneck is in **doublet-to-bulk** matrix elements / Boltzmann factors, not in the bare within-doublet ones. Follow-up B (kinky-Metro spot check) confirms the β-non-monotonicity in λ_L(β) is **physical**, not a smooth-Metro filter artefact.

## Check 1 — β sweep at n=6 (2×3 ladder, h=1)

Spans deep-disordered (T/T_c ≈ 4.8) through T_c (β_c ≈ 0.48) to deep-ordered (T/T_c ≈ 0.16). All numerics use the canonical CKG smooth-Metro (s=0.25, a=0, r_D=7, kdim=40), no Bohr cross-check (validated for n=6 in qf-1jj), no trajectory τ_mix.

| β_phys | T/T_c | β_alg | ⟨M_z²⟩ | U_4 | S/log d | eff_rank | w_+ | w_- | Δ_1^phys | λ_L^phys | λ_L/ΔE_1 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.10 | 4.83 |   5.04 | 0.253 | 0.222 | 0.970 | 56.5 | 0.500 | 0.500 | 1.07e-2 | 4.00 | 373 |
| 0.25 | 1.93 |  12.6  | 0.470 | 0.435 | 0.800 | 27.9 | 0.500 | 0.500 | 1.07e-2 | 1.57 | 147 |
| 0.50 | 0.97 |  25.2  | 0.822 | 0.617 | 0.387 |  4.99 | 0.501 | 0.499 | 1.07e-2 | 0.286 | 26.7 |
| 1.00 | 0.48 |  50.4  | 0.942 | 0.655 | 0.177 |  2.09 | 0.503 | 0.497 | 1.07e-2 | 1.24e-2 | 1.15 |
| 1.50 | 0.32 |  75.7  | 0.945 | 0.656 | 0.167 |  2.00 | 0.504 | 0.496 | 1.07e-2 | 2.65e-3 | 0.247 |
| 2.00 | 0.24 | 101    | 0.945 | 0.656 | 0.167 |  2.00 | 0.505 | 0.495 | 1.07e-2 | 3.06e-3 | 0.286 |
| 2.50 | 0.19 | 126    | 0.945 | 0.656 | 0.167 |  2.00 | 0.507 | 0.493 | 1.07e-2 | 4.18e-3 | 0.390 |
| 3.00 | 0.16 | 151    | 0.945 | 0.656 | 0.167 |  2.00 | 0.508 | 0.492 | 1.07e-2 | 5.57e-3 | 0.520 |

**Phase verdict**: β=2 is deep in the ordered phase — already at the saturated values of every order parameter from β=1.5 onward. eff_rank = 2.00 exactly is the cleanest signature: the Gibbs state is supported on the Z₂ doublet, not on a single GS (as in the 1D Heisenberg even-n p_0 ≈ 99% case). Crossover happens at β ≈ 0.5 (T/T_c ≈ 1), in agreement with Hesselmann–Wessel 2016. **No need to go colder for the ordered-phase signature** — saturated already.

**Unexpected observation — L gap is non-monotone in β** (resolved by Follow-up B below). The Lindbladian gap λ_L^phys reaches a minimum near β_phys ≈ 1.5 and then *increases* as we go colder: 1.24e-2 (β=1) → 2.65e-3 (β=1.5) → 3.06e-3 (β=2) → 4.18e-3 (β=2.5) → 5.57e-3 (β=3). Naively one expects a monotone decrease at low T. Follow-up B (kinky Metropolis spot check, see below) shows this is **not** a smooth-Metro filter artefact — the kinky-Metro gap follows the same shape to within 4 % at every β. The non-monotonicity is therefore physical (most likely a slow-mode crossover, e.g. from "thermal Glauber barrier crossing" to "kinetic tunnelling" as β grows, well-documented for exponentially-small-gap quantum-annealing problems [Bando et al. 2018, PMC6060131 / arXiv:1801.02006]).

## Check 2 — matrix elements between Z₂ doublet, n=4 and n=6 (ORDERED β_phys=2)

Dense `eigen(H_phys)`; |ψ_1⟩, |ψ_2⟩ are the lowest two eigenvectors. Z₂ parities verified at machine precision (⟨ψ_1|P|ψ_1⟩ = +1, ⟨ψ_2|P|ψ_2⟩ = -1; P = X^⊗n).

Matrix element of each canonical 1/√(3n)-normalised jump operator A_a:

$$M^2 = \sum_a |\langle \psi_2 | A_a | \psi_1\rangle|^2,
\qquad
M^2_\text{bare} = \sum_{i,\alpha} |\langle \psi_2 | \sigma_\alpha^{(i)} | \psi_1\rangle|^2$$

with α ∈ {X, Y, Z}, i ∈ {1..n}, 3n jump operators total.

| n | ΔE_1^phys | ΔE_1^alg | M²_norm | M²_bare | \|X\|²_sum | \|Y\|²_sum | \|Z\|²_sum | γ(+ΔE_1) | γ(-ΔE_1) | λ_L^phys | λ_L^alg |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 4 | 7.10e-2 | 1.87e-3 | 0.311 | 3.73 | ~10⁻³⁰ | 4.69e-3 | 3.73 | 0.510 | 0.642 | 9.96e-2 | 2.62e-3 |
| 6 | 1.07e-2 | 2.12e-4 | 0.311 | 5.60 | ~10⁻²⁸ | 1.61e-4 | 5.60 | 0.565 | 0.585 | 3.06e-3 | 6.08e-5 |
| **ratio n6/n4** | **0.151** | — | **1.001** | **1.501** | — | **0.034** | **1.503** | — | — | **0.0308** | **0.0232** |

### Interpretation

1. **σ_x matrix elements are EXACTLY zero** (~10⁻³⁰, machine precision). Mechanism: P = X^⊗n commutes with σ_x^(i) (Z₂-even), so σ_x can't connect P=+1 |ψ_1⟩ to P=-1 |ψ_2⟩.
2. **σ_z matrix elements are O(1)**. Each |⟨ψ_2|σ_z^(i)|ψ_1⟩|² ≈ 0.93 at both n. Mechanism: σ_z is Z₂-odd (anti-commutes with P at the matching site, commutes elsewhere), and σ_z just changes the sign of the |↓...↓⟩ component in the dressed-FM-minimum decomposition, so the leading-order matrix element saturates at 1 in the limit |ψ_{1,2}⟩ → (|↑...↑⟩ ± |↓...↓⟩)/√2. The per-site value of 0.93 (vs 1.00) is the small dressing correction at h/J = 1.
3. **σ_y matrix elements are small and DO close with n** — 4.7e-3 → 1.6e-4 from n=4 to n=6. This is the σ_y ~ i σ_x σ_z perturbative tail, which involves both spin-flip and sign-flip and gets exponentially suppressed in the dressed FM minima.
4. **M²_bare scales linearly with n** (3.73 → 5.60 = 1.50× = 6/4 exactly). **Sum-rule explanation**: M²_bare is dominated by σ_z, each contributing ≈ 1 per site (claim 2), summed over n sites gives M²_bare ≈ n. The 1/(3n) normalisation absorbs the linear-in-n growth, so M²_norm is constant at ≈ 1/3 independent of system size. *This is exactly the right invariance to test for a within-doublet bottleneck — if it broke, it'd be evidence of dressing corrections to the bare FM-minimum picture, not of a kinetic bottleneck.*

### The matrix-element bottleneck story fails

If λ_L ≈ M² · γ_KMS(±ΔE_1) · (n-independent prefactor) held, the L gap should track M² × γ — predicted ratio 0.911 from n=4 to n=6. The actual ratio is 0.023, **off by a factor of 40**. The closing factor of 43× in λ_L^alg cannot be attributed to either γ (ratio 0.91) or M² (ratio 1.00).

So the qualitative statement "L gap closes faster than ΔE_1 → matrix-element bottleneck" was **wrong as a mechanism**. The phenomenon (λ closes faster than ΔE_1) is real; the explanation isn't the bare-doublet matrix elements.

## Follow-up A — what IS the slow mode? (dense L eigendecomp at n=4)

`scripts/scratch_qf_biz_followup.jl` — `construct_lindbladian(jumps, cfg, ham; include_coherent=true)` at n=4, ORDERED β=2, d=16, d²=256. `eigen(L_dense)` in 0.13 s. Sort by |Re λ|.

**Spectrum**: lowest six |Re λ_alg|: `6.4e-17  2.62e-3  0.149  0.203  0.203  0.219`. The λ_1 ≈ 0 fixed point checks out at floating-point noise; λ_2 = -2.6243e-3 + i·O(10⁻²⁸) reproduces the Krylov gap to all displayed digits.

**Operator-basis decomposition of the slowest mode R_2** (Frobenius-normalised, viewed as a d×d matrix):

| Quantity | Value | Interpretation |
|---|---|---|
| Z₂ parity `tr(R_2† P R_2 P) / ‖R_2‖²` | **+0.000** | NOT a pure Z₂-parity eigenmode — mixed sector. |
| \|⟨R_2, (\|ψ_1⟩⟨ψ_2\| + h.c.)/√2⟩\| | 0.119 | Small — within-doublet symmetric coherence is *not* the slow mode. |
| \|⟨R_2, i(\|ψ_1⟩⟨ψ_2\| − h.c.)/√2⟩\| | 0.000 | Zero — neither is the anti-symmetric coherence. |
| \|⟨R_2, (\|ψ_1⟩⟨ψ_1\| − \|ψ_2⟩⟨ψ_2\|)/√2⟩\| | 0.007 | Negligible — and not the doublet imbalance. |
| \|⟨R_2, M_z⟩\| (Frob-normalised) | 0.000 | Zero — slow mode has zero overlap with the magnetisation. |

**Block decomposition of \|R_2\|² in the energy eigenbasis** (i, j ∈ {1, 2} ≡ doublet; > 2 ≡ bulk):

| Block | Mass fraction |
|---|---|
| doublet × doublet | **2.8 %** |
| **doublet × bulk** | **89.4 %** |
| bulk × bulk | 7.8 % |

**Top eight \|R_2\|-entries in the energy eigenbasis**: all of the form `(2, 8), (8, 2), (1, 8), (8, 1), (2, 13), (13, 2), (2, 9), (9, 2)` — i.e. `|ψ_{1,2}⟩⟨ψ_{8,9,13}|`-type matrix elements.

### So the slow mode is

**a coherence between the Z₂-doublet ground states and a specific bulk excited state |ψ_8⟩** (with smaller contributions from |ψ_9⟩, |ψ_13⟩). Magnetisation- and Z₂-orthogonal. The associated Bohr frequencies are ω ≈ E_8 − E_{1,2} > 0, which are energy-increasing transitions; KMS-DB Boltzmann-suppresses γ(+ω) ~ exp(−β·ω). At β_phys=2, n=4, E_3 − E_1 ≈ 5.9 (phys), so β·ω ≈ 12 and γ ~ 10⁻⁶ for the bulk-excitation rate — *but* the L gap is ~ 10⁻³, so the L-gap rate is **not** just γ × ME² — there's a residual scale (likely from the destructive sum over multiple bulk states, or from a Z₂-mixed slow direction that decouples from M_z).

### Mechanism take-away

Both candidate stories I sketched in the discussion earlier are partially wrong:

- *Within-doublet matrix-element bottleneck* — refuted (M² constant, slow mode lives 89 % outside doublet×doublet).
- *Boltzmann barrier crossing on a single classical path* — too simple; doesn't quantitatively predict the gap.

The right framing is: in the ORDERED phase the slowest L mode is an **off-diagonal doublet-to-bulk coherence**, magnetisation-orthogonal and Z₂-mixed. Its decay rate is set by the full Davies-rate-equation eigenvalue, not by any single bare matrix element. The closing with n is real thermal-order physics, but quantitative scaling needs the full eigenvalue, not the bare doublet matrix element.

## Follow-up B — kinky-Metropolis spot check at n=6

`scripts/scratch_qf_biz_followup.jl` — repeat Check 1 with **kinky Metropolis (s=0, a=0)** at the same five β_phys ∈ {1.0, 1.5, 2.0, 2.5, 3.0}. Goal: verify the non-monotonicity in λ_L(β) is not a smooth-Metro filter artefact.

| β_phys | λ_phys kinky | λ_phys smooth | ratio |
|---|---|---|---|
| 1.0 | 1.265e-2 | 1.237e-2 | 1.02 |
| 1.5 | 2.752e-3 | 2.645e-3 | 1.04 |
| 2.0 | 3.175e-3 | 3.065e-3 | 1.04 |
| 2.5 | 4.297e-3 | 4.182e-3 | 1.03 |
| 3.0 | 5.683e-3 | 5.574e-3 | 1.02 |

Kinky and smooth filters agree to 4 % at every β, with **identical non-monotonicity shape** (both bottom out at β=1.5 and rise to β=3). The β-non-monotonicity is **physical**, not a filter artefact. Most likely the slow-mode crossover (thermal-Glauber barrier crossing → kinetic tunnelling) that has been documented for finite-T quantum-annealing problems on exponentially-small-gap Hamiltonians [Bando et al. 2018, PMC6060131]. Detailed mechanism beyond the scope of qf-biz.

## Implications for the qf-1jj draft

The qf-1jj draft already claimed correctly that λ_L^phys closes faster than ΔE_1^phys (ratio 1.4 → 0.29 → 0.086 across n=4,6,8). That part stands. The qf-1jj draft does **not** make any inline claim that the closing is due to ⟨ψ_2|A|ψ_1⟩ (physics-checker confirmed this in the sweep) — so no edits to qf-1jj headlines are required. The "Methodology notes" section can usefully gain a brief footnote:

> *qf-biz follow-up (2026-05-24): the within-doublet matrix element $\sum_a |\langle \psi_2|A_a|\psi_1\rangle|^2$ is constant in n (≈ 0.31) and does **not** account for the gap closing. The slow mode of L lives 89 % in the doublet × bulk off-diagonal block at n=4 — a Z₂-mixed coherence between |ψ_1,2⟩ and bulk excited states like |ψ_8⟩. The thermal-order signature is genuine but the bottleneck is at doublet ↔ bulk transitions, not bare within-doublet matrix elements.*

Add the qf-biz draft to the cross-references list at the bottom of qf-1jj.

## Follow-up C — canonical $\rho_0 = |+\rangle\langle+|^{\otimes N}$ rerun (qf-b4i / qf-65e, 2026-05-24)

The qf-1jj sweep ran the trajectory predictor from $\rho_0 = I/d$, which was the default at the time. The project-wide canonical choice (qf-e4z.30, decreed 2026-05-24) is $\rho_0 = (|+\rangle\langle+|)^{\otimes N}$. qf-b4i re-runs the six cells with the canonical $\rho_0$ keeping every other knob fixed; gap_phys is $\rho_0$-independent and unchanged (verified to rel-diff $\le 1.1\,\text{e}{-10}$ vs qf-1jj at all six cells), but the trajectory τ_mix moves substantially in the ordered phase.

| $n$ | $\tau_{\rm mix}^{\rm ORD}\ (|+\rangle)$ | $\tau_{\rm mix}^{\rm ORD}\ (I/d)$ | $|c_2|^{\rm ORD}\ (|+\rangle)$ | mechanism |
|---|---|---|---|---|
| 4 | $1349$ | $47.0$ | $4.88\,\text{e}{-2}$ | $|+\rangle$ couples to tunnelling mode — $\tau \sim 1/\lambda_2$ |
| 6 | $6.23\,\text{e}{4}$ | $104.7$ | $6.23\,\text{e}{-2}$ | same — $\tau$ tracks the Z₂-tunnelling collapse |
| 8 | $251.8$ | $190.1$ | $1.34\,\text{e}{-3}$ | $|+\rangle$ symmetry-decouples from L's slow mode — $\tau \sim 1/\lambda_3$ |

**Mechanism.** Under $P = X^{\otimes N}$, both $\rho_0$'s are parity-EVEN ($P|+\rangle\langle+|^{\otimes N} P^\dagger = |+\rangle\langle+|^{\otimes N}$ since $X|+\rangle = |+\rangle$; $P \cdot I/d \cdot P^\dagger = I/d$). The slow Z₂-odd eigenmode lives 89 % in the $|\psi_{1,2}\rangle\langle\psi_8|$-type off-diagonal block (Follow-up A above) — that block has substantial overlap with the X-only Pauli content of $|+\rangle\langle+|^{\otimes N}$ at small $n$ where the doublet states have spatial breadth, but at $n = 8$ the doublet states are essentially Hamming-weight-pure on the $|0...0\rangle$ / $|1...1\rangle$ classical ground states (up to the exponentially small tunnelling). The off-diagonal block becomes nearly Y/Z-string-orthogonal to the broad X-string content, $|c_2|$ collapses $46\times$ between $n = 6$ and $n = 8$, and $\tau_{\rm mix}(|+\rangle)$ joins $\tau_{\rm mix}(I/d)$ at the bulk-mode timescale.

**Reading**: τ_mix($|+\rangle$) is the canonical mixing time for the canonical initial state. The qf-1jj numbers (τ_mix($I/d$) modest at all $n$) reflected $I/d$'s exact decoupling from the parity-odd slow mode by symmetry, not the "true" Lindbladian mixing rate. Neither symmetric $\rho_0$ tracks $1/\lambda_2$ at large $n$; in a hardware run one needs a parity-asymmetric initial state (e.g. one of the doublet states $|\psi_1\rangle\langle\psi_1|$, or a thermal-perturbed seed) for $\tau$ to track the full gap collapse.

**qf-65e Krylov-size correction.** The initial qf-b4i rerun used krylovdim $= 40$ (matching qf-1jj) and exposed the ORD $n = 8$ trajectory to a Pass-1 truncation artefact: at kdim = 40 the trajectory predictor reported a "slowest captured mode" at $\lambda^{\rm phys} = 4.77\,\text{e}{-3}$, $30\times$ larger than L's true slowest mode ($\lambda_2^{\rm phys} = 1.57\,\text{e}{-4}$), with floor $= 9.3\,\text{e}{-4}$ (uncomfortably close to $\varepsilon = 10^{-3}$). The qf-65e rerun bumped Pass-1 to the Heisenberg-recipe krylovdim$_{p1} = 80$ (Pass-2 fixed at $30$ per qf-e4z.33); ORD $n = 8$ floor dropped $170\times$ to $5.5\,\text{e}{-6}$, the Pass-1↔Pass-2 gap ratio is now $1.0$ exactly, and τ_mix corrected from $9228$ (spurious) to $251.8$ (the bulk-relaxation regime above). A saturation check at krylovdim$_{p1} = 100$ reproduces the krylovdim$_{p1} = 80$ result bit-identically (Pass-2 gap, traj-eff gap, floor, τ_mix, $|c_2|$ all match to 7 figures). **Verdict: undersizing, not parity-trap** — Pass-1 at kdim = 80 does capture L's slowest mode in its Krylov subspace; the $|c_2| = 1.3\,\text{e}{-3}$ coupling is real and tiny, but non-zero. The kdim sensitivity sidecar is at `scripts/output/sweep_2d_tfim_ordered_vs_disordered/qf_b4i_n8_kdim_diagnostic.bson`.

**Direct $d(t)$ verification of the n=6 ↔ n=8 crossover.** Script `scripts/scratch_qf_65e_verify_ord_n8.jl` directly evaluates $d(t) = \tfrac{1}{2}\|\rho_{\rm eigen}(t) - \sigma_\beta\|_1$ at many $t$ on ORD n=8 (krylovdim$_{p1} = 80$). $d(t)$ is monotone-decreasing: $d(0) = 0.997 \to d(30) = 0.405 \to d(200) = 1.29\,\text{e}{-3} \to \mathbf{d(252) = 1.00\,\text{e}{-3} = \varepsilon} \to$ "plateau" $d(300 \le t \le 10^4) \approx 9.5\,\text{e}{-4} \to d(5 \times 10^5) = 3.4\,\text{e}{-4}$ (tunnelling decay onset, $1/\lambda_2 \approx 4.8 \times 10^5$) $\to d(5 \times 10^6) = 5.5\,\text{e}{-6}$ (floor). The plateau value is $|c_2| \cdot \|R_2\|_1 / 2 = 1.34\,\text{e}{-3} \cdot \sqrt{2} / 2 = 9.51\,\text{e}{-4}$, **just below $\varepsilon$** — narrow-escape regime. At n=6 the plateau is $\sim 4\,\text{e}{-2}$ (well above $\varepsilon$) → tunnelling-limited, τ_mix $\approx 6 \times 10^4$. At n=8 it just clears $\varepsilon$ → bulk-limited, τ_mix = 252. The sharp non-monotonicity is a **genuine physical crossover** at the $|c_2| \cdot \|R_2\|_1 / 2 \approx \varepsilon$ threshold, not a numerical artefact.

## Follow-up D — $\rho_0 = |0\rangle\langle 0|^{\otimes N}$ rerun (qf-b4i.1, 2026-05-25)

Hypothesis (issue-description): on 2D TFIM the Z₂ generator is $P = X^{\otimes N}$, so $X|+\rangle = |+\rangle$ ⇒ $|+\rangle\langle+|^{\otimes N}$ is X-parity-EVEN under conjugation by $P$ (same parity character as $I/d$). $X|0\rangle = |1\rangle \ne |0\rangle$ ⇒ $|0\rangle\langle 0|^{\otimes N}$ is X-parity-MIXED *as a density operator* ($P \rho_0 P^\dagger = |1\rangle\langle 1|^{\otimes N} \ne \rho_0$). The doublet states approach $|\psi_{1,2}\rangle \approx (|0...0\rangle \pm |1...1\rangle)/\sqrt{2}$ in the deep ordered phase, giving $\langle\psi_{1,2}|0...0\rangle \to 1/\sqrt{2}$ — a macroscopic overlap with both doublet states at every $N$. So *if* the slow Lindbladian mode were the within-doublet coherence $|\psi_1\rangle\langle\psi_2|$, $|0\rangle\langle 0|^{\otimes N}$ would couple to it with $|c_2| \sim 1/2$ at every $N$ — no orthogonality catastrophe — and $\tau_{\rm mix}$ would track $1/\lambda_2$ at the bottom of the gap collapse (e.g. $\sim 4 \times 10^5$ phys at ORD $n = 8$).

Method: `scripts/scratch_qf_b4i_1_2d_tfim_taumix_zero.jl` reruns the six cells via `_krylov_spectral_decomposition` (Pass-1 only; krylovdim$_{p1} = 80$). Pass-2 is skipped — we cross-check the Pass-1 traj-effective gap_phys directly against the qf-65e Pass-2 gap_phys ($\rho_0$-independent). Audit threshold: $\le 1\%$. Result:

| phase | $n$ | gap_phys(qf-b4i.1, traj-eff) | gap_phys(qf-65e, Pass-2) | rel-diff | $\tau_{\rm mix}(|0\rangle)$ alg | $\tau_{\rm mix}(|+\rangle)$ alg | $\tau_{\rm mix}(I/d)$ alg | $|c_2|(|0\rangle)$ | floor |
|---|---|---|---|---|---|---|---|---|---|
| DISORDERED | 4 | 1.910031 | 1.910031 | $7.0 \times 10^{-15}$ | 33.3 | 86.6 | 33.9 | $9.3 \times 10^{-15}$ | $2.5 \times 10^{-16}$ |
| DISORDERED | 6 | 1.574719 | 1.574719 | $2.0 \times 10^{-14}$ | 63.8 | 143.2 | 59.7 | $5.5 \times 10^{-13}$ | $3.6 \times 10^{-16}$ |
| DISORDERED | 8 | 1.740497 | 1.740497 | $2.4 \times 10^{-14}$ | 102.3 | 89.7 | 88.9 | $1.5 \times 10^{-11}$ | $1.2 \times 10^{-15}$ |
| ORDERED | 4 | $9.9645 \times 10^{-2}$ | $9.9645 \times 10^{-2}$ | $1.9 \times 10^{-13}$ | 17.1 | 1349 | 47.0 | $2.9 \times 10^{-15}$ | $3.1 \times 10^{-7}$ |
| ORDERED | 6 | $3.0648 \times 10^{-3}$ | $3.0648 \times 10^{-3}$ | $1.5 \times 10^{-11}$ | 17.3 | $6.23 \times 10^4$ | 104.7 | $8.9 \times 10^{-13}$ | $1.2 \times 10^{-6}$ |
| ORDERED | 8 | $1.5736 \times 10^{-4}$ | $1.5736 \times 10^{-4}$ | $3.6 \times 10^{-10}$ | 17.4 | 251.8 | 190.1 | $2.9 \times 10^{-11}$ | $5.5 \times 10^{-6}$ |

All six audits pass at $\sim 10^{-14}$ to $10^{-10}$ — Pass-1 from $|0\rangle$ captures L's slowest eigenvalue. **But the hypothesis is refuted: $|c_2|(|0\rangle) \le 10^{-11}$ at every cell**, vastly smaller than $|c_2|(|+\rangle)$ ($\sim 10^{-1}$ at small-$n$ ORD), and $\tau_{\rm mix}(|0\rangle)^{\rm ORD}$ is *constant at $\sim 17$ alg time units* across $n \in \{4, 6, 8\}$ — orders of magnitude shorter than any other column.

### Why the hypothesis fails (mechanism)

The hypothesis assumed L's slow mode is the within-doublet coherence $R_2 \propto |\psi_1\rangle\langle\psi_2|$. Follow-up A above already refuted this: at ORD $n=4$, $R_2$ lives **89% in the doublet × bulk off-diagonal block** (top entries $|\psi_{1,2}\rangle\langle\psi_8|$-type), only 2.8% in the doublet × doublet block.

The deep-ordered-phase decomposition of $|0\rangle\langle 0|^{\otimes N}$ in the energy eigenbasis is

$$|0...0\rangle\langle 0...0| \approx \tfrac{1}{2}\bigl(|\psi_1\rangle\langle\psi_1| + |\psi_1\rangle\langle\psi_2| + |\psi_2\rangle\langle\psi_1| + |\psi_2\rangle\langle\psi_2|\bigr) + \mathcal{O}(\text{tunnelling tails}),$$

i.e. *entirely* on the doublet × doublet block plus exponentially small bulk tails. R_2 lives in the doublet × bulk block — orthogonal to where $|0\rangle\langle 0|^{\otimes N}$ has weight. Hence $|c_2| \to 0$ even faster than for $|+\rangle\langle+|^{\otimes N}$ (which has *some* doublet × bulk content via its broad X-string Pauli mass).

Mechanically the picture for $|0\rangle\langle 0|^{\otimes N}$ at ORDERED is:

1. $\rho_0$ starts as a rank-1 projector concentrated on the doublet × doublet block.
2. The "fast" relaxation modes ($R_k$ for $k \ge 3$ with $|c_k|$ noticeable) **flatten the doublet diagonal** from $\{p_{|\psi_1\rangle} = 1, p_{|\psi_2\rangle} = 0\}$ to $\{1/2, 1/2\}$ and **zero out the within-doublet coherence** $\langle\psi_1|\rho|\psi_2\rangle$. Both happen at the bulk-mode rate $\lambda_3 \sim 10^{-2}$ phys (or faster — for $|0\rangle$ specifically, $|c_3|$ is itself tiny, and the actual relaxation is driven by $R_k$ for $k$ in the bulk).
3. Once the doublet content has thermalised, $\rho \approx \sigma_\beta$ (since $\sigma_\beta \approx \text{diag}(1/2, 1/2)$ on the doublet at $\beta = 2$, eff_rank = 2.00 from Check 1), and $d(t) < \varepsilon$ is reached. R_2 (the doublet × bulk coherence) was never excited; it neither needs to decay nor matters for τ_mix.

Quantitatively at ORD $n = 8$: $\tau_{\rm mix}(|0\rangle)^{\rm alg} = 17.4$, and the dominant decay mode in the eigenvalue-amplitude product is some $R_k$ with $|c_k| \sim |c_3|^{|0\rangle}$ at amplitude $\mathcal{O}(10^{-5})$ and eigenvalue $\lambda_k \sim 10^{-2}$ alg, giving $\tau \sim \log(|c_k| / \varepsilon) / \lambda_k = \log(10^{-5} / 10^{-3}) / 10^{-2} \approx 5 / 10^{-2}$ — order-of-magnitude $\sim 500$, *but in the limit where the doublet × doublet block content is most of $\rho_0$*, the trace distance is dominated by the doublet projection itself which thermalises faster*. The constant-with-$n$ value $\sim 17$ alg suggests a single $n$-independent local bulk-relaxation rate dominates — consistent with the doublet-block content being independent of $n$ (always a rank-1 projector on the doublet × doublet $2 \times 2$ block).

### Implication for the canonical-$\rho_0$ choice on 2D TFIM

The issue's hypothesis — "choose a single-qubit factor that is not a $P$-eigenstate" — captures the *necessary* condition for X-parity breaking at the operator level but NOT the *sufficient* condition for coupling to L's slow mode. The relevant criterion is whether $\rho_0$ has weight on the same operator-space block as $R_2$, which on this 2D TFIM model is the **doublet × bulk off-diagonal block**.

- *None of the three audit columns satisfies that.* $|+\rangle\langle+|^{\otimes N}$ has only-X-Pauli-string content (orthogonal to Y/Z-string mass of R_2 at large $n$); $|0\rangle\langle 0|^{\otimes N}$ has only-diagonal-Z content (lives entirely in the doublet × doublet block); $I/d$ has weight only on the identity Pauli (decoupled from every $R_k$). All three τ_mix's revert to bulk-mode timescales at large $n$.
- The natural initial states that *would* track $1/\lambda_2$ are (a) one of the doublet states $|\psi_1\rangle\langle\psi_1|$ (broken X-parity *and* doublet × bulk content via the small bulk tails), (b) a thermal-perturbed Gibbs seed $\sigma_\beta + \alpha \cdot R_2$, or (c) any random Gaussian-perturbed state with O(1) projection on R_2.

For the project-canonical-$\rho_0$ memory entry ($\rho_0 = |+\rangle\langle+|^{\otimes N}$, decreed 2026-05-24): the **2D TFIM result does not motivate a model-specific switch to $|0\rangle\langle 0|^{\otimes N}$**. The 2D TFIM mixing-time anomaly is fundamentally a property of how concentrated $R_2$'s support is in operator space (doublet × bulk, low-rank, hard to reach from any natural product state), not a parity-protection artefact removable by choosing a different $\rho_0$ in the same product-state family. The canonical $|+\rangle\langle+|^{\otimes N}$ rule stays — with the documented understanding that on 2D TFIM (and any other model where $R_2$ lives in a low-rank off-diagonal block far from the product-state manifold), $\tau_{\rm mix}$ from the canonical $\rho_0$ is the bulk-mode timescale, not $1/\lambda_2$.

## Acceptance (qf-biz)

1. ✓ β-sweep at n=6 confirms ORDERED at β ≥ 1 with eff_rank = 2.00; crossover at β ≈ 0.5 (T_c per HW 2016).
2. ✓ At β_phys=2 deep ordered (⟨M_z²⟩ = 0.945, U_4 = 0.656, eff_rank = 2.00). **No need to go colder** — saturated.
3. ✓ Matrix-element ratio M²(n=6)/M²(n=4) = 1.001 — within-doublet matrix elements are NOT the bottleneck (refuted).
4. ✓ Follow-up A: slow mode identified — Z₂-mixed coherence in the doublet × bulk block (89 % of |R_2|²), magnetisation-orthogonal.
5. ✓ Follow-up B: β-non-monotonicity in λ_L(β) is **physical** (kinky-Metro reproduces shape to 4 %), not a filter artefact.
6. ✓ Driver + BSON sidecars written; `drafts/qf-biz-code-verifier-report.md` PASS bit-exact against qf-1jj; physics-checker validated claims 1–4 algebraically and amended claim 5 (single-classical-path Glauber → "Davies-rate-eigenvalue with doublet-to-bulk dominance" — captured in Follow-up A).
7. ✓ Follow-up D (qf-b4i.1): $\rho_0 = |0\rangle\langle 0|^{\otimes N}$ rerun confirms the qf-biz Follow-up A picture independently. The X-parity-breaking-density-operator hypothesis (issue qf-b4i.1) predicted $|c_2| \sim 1/2$ at ORD; observed $|c_2| \le 10^{-11}$ at every cell. The mechanism is that $R_2$ lives in the doublet × bulk block, while $|0\rangle\langle 0|^{\otimes N}$ at deep-ORDERED is supported on the doublet × doublet block — orthogonal. **Confirms Follow-up A independently** (the slow mode is NOT within-doublet) and **closes the symmetry-protection investigation** for this 2D TFIM model: no symmetric product-state ρ_0 in $\{I/d, |+\rangle, |0\rangle\}$ couples to $R_2$, regardless of operator-level parity character.

## Cross-references

- qf-1jj — parent draft `drafts/2d-tfim-ordered-vs-disordered.md`.
- qf-8fr — Krylov x_0 GUE seed (`.claude-memory/krylov_x0_symmetric_bug_qf_8fr.md`).
- qf-b4i — canonical $\rho_0 = |+\rangle\langle+|^{\otimes N}$ trajectory rerun. `scripts/scratch_qf_b4i_rerun_2d_tfim_taumix_plus.jl`.
- qf-65e — Krylov-size correction at $n = 8$ ORDERED: krylovdim$_{p1} = 80$, krylovdim$_{p2} = 30$. `scripts/scratch_qf_b4i_n8_kdim_diagnostic.jl`.
- qf-b4i.1 — $\rho_0 = |0\rangle\langle 0|^{\otimes N}$ X-parity-breaking trial (Follow-up D). `scripts/scratch_qf_b4i_1_2d_tfim_taumix_zero.jl`.
- `.claude-memory/canonical_taumix_setup_qf_e4z_30.md` — canonical $\rho_0$ policy (project-wide).
- `.claude-memory/qf_e4z_33_kdim_convergence.md` — Heisenberg-recipe krylovdim$_{p1} \ge 80$ for $n = 7$ low-$\beta$ cells.
- `.claude-memory/heisenberg_1d_even_odd_mechanism.md` — the 1D Heisenberg even-n p_0 ≈ 99% artefact (distinct mechanism, not at play here).
- `.claude-memory/feedback_check_gibbs_when_simulation_off.md` — why Subtask B' Gibbs check is the right first move before invoking algorithmic explanations.
