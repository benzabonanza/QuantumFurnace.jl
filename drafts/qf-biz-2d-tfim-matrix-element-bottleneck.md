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

## Acceptance (qf-biz)

1. ✓ β-sweep at n=6 confirms ORDERED at β ≥ 1 with eff_rank = 2.00; crossover at β ≈ 0.5 (T_c per HW 2016).
2. ✓ At β_phys=2 deep ordered (⟨M_z²⟩ = 0.945, U_4 = 0.656, eff_rank = 2.00). **No need to go colder** — saturated.
3. ✓ Matrix-element ratio M²(n=6)/M²(n=4) = 1.001 — within-doublet matrix elements are NOT the bottleneck (refuted).
4. ✓ Follow-up A: slow mode identified — Z₂-mixed coherence in the doublet × bulk block (89 % of |R_2|²), magnetisation-orthogonal.
5. ✓ Follow-up B: β-non-monotonicity in λ_L(β) is **physical** (kinky-Metro reproduces shape to 4 %), not a filter artefact.
6. ✓ Driver + BSON sidecars written; `drafts/qf-biz-code-verifier-report.md` PASS bit-exact against qf-1jj; physics-checker validated claims 1–4 algebraically and amended claim 5 (single-classical-path Glauber → "Davies-rate-eigenvalue with doublet-to-bulk dominance" — captured in Follow-up A).

## Cross-references

- qf-1jj — parent draft `drafts/2d-tfim-ordered-vs-disordered.md`.
- qf-8fr — Krylov x_0 GUE seed (`.claude-memory/krylov_x0_symmetric_bug_qf_8fr.md`).
- `.claude-memory/heisenberg_1d_even_odd_mechanism.md` — the 1D Heisenberg even-n p_0 ≈ 99% artefact (distinct mechanism, not at play here).
- `.claude-memory/feedback_check_gibbs_when_simulation_off.md` — why Subtask B' Gibbs check is the right first move before invoking algorithmic explanations.
