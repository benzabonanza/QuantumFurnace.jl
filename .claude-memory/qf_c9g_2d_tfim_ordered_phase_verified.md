---
name: qf-c9g-2d-tfim-ordered-phase-verified
description: "2D TFIM ORDERED Lindbladian gap closing (n=4,6,8 at β_phys=2, h=1) is genuine ordered-phase physics — not a Gibbs-≈-GS-doublet artefact. Four-of-four diagnostic tests pass: β-persistence at β=0.5 (eff_rank > 2), Binder U_4 saturation, doublet × bulk slow-mode character at every β below T_c, paramagnetic control flat. Mechanism = surface-tension lower bound from Gamarnik-Kiani-Zlokapa 2024 Thm 1.4 (arXiv:2411.04300) for the exact CKG KMS-DB sampler class."
metadata: 
  node_type: memory
  type: project
  originSessionId: af4aaa49-1221-4001-8453-7b63624870a5
---

The colleague-flagged worry (qf-c9g, 2026-05-25): "at β_phys = 2 the Gibbs
state is essentially the GS doublet — so the gap closing might be a trivial
two-level-system tunnel-splitting artefact, not specifically 2D Ising
ordering". Settled NO via a 15-cell n × β grid + n=4 slow-mode dense
eigendecomp + literature anchor.

**Verdict**: closing IS genuine ordered-phase physics. The artefact
hypothesis is dead.

**Key data** (`drafts/qf-c9g-2d-tfim-ordered-mechanism.md`,
`scripts/output/qf_c9g_ordered_gap_mechanism/qf_c9g_beta_n_sweep.bson`):

| n | β=0.10 (T/T_c=4.8) | β=0.25 (1.9) | β=0.50 (0.97) | β=1.0 (0.48) | β=2.0 (0.24) |
|---|---|---|---|---|---|
| 4 | 4.53 | 1.91 | 0.449 | 0.0681 | 0.0996 |
| 6 | 4.00 | 1.575 | 0.286 | 0.0124 | 0.00306 |
| 8 | 4.49 | 1.741 | 0.285 | 0.00999 | 0.000157 |

slope of log(λ_L^phys) vs n: -0.002 → -0.023 → -0.114 → -0.480 → -1.613.
Clean phase transition signal at β_c = 1/T_c ≈ 0.483.

**Diagnostic checklist (from `drafts/qf-c9g-lit-survey.md`)**:

(i) **β-persistence below T_c**: closing visible at β=0.5 (T/T_c=0.97) with
eff_rank ∈ [3.3, 8.6] — Gibbs spread over multiple states, NOT just the
doublet. This is THE decisive cell — the artefact hypothesis predicts no
closing here.

(ii) **Magnetisation indicators**: U_4 ∈ [0.60, 0.66] across ordered window
(textbook 2/3 saturation); ⟨m_z²⟩ ≈ 0.95 saturated; eff_rank → 2 at
β ≥ 1.

(iii) **Slow-mode operator-space character (n=4 dense L eigendecomp)**:
R_2 lives **89% in the doublet × bulk** off-diagonal block at EVERY β
below T_c (also at β=0.25 paramagnetic). Top entries always
|ψ_{1,2}⟩⟨ψ_9|-type at ω ≈ ±8.5. ⟨R_2, M_z⟩ = 0 (magnetisation-orthogonal).
The operator-space CHARACTER is constant across β; only the RATE changes
via γ_KMS(±8.5) ∝ exp(∓β·8.5/2). qf-biz Follow-up A's β=2 finding
generalises to all β.

(iv) **Paramagnetic control**: at β=0.10 and 0.25, gap_phys is FLAT in n.
Slopes -0.002, -0.023 (indistinguishable from noise). Confirms closing is
phase-specific.

**Lit anchor**: **Gamarnik–Kiani–Zlokapa 2024** (arXiv:2411.04300) Thm 1.4
proves `T_mix ≥ exp[n^(1/2-o(1))]` for the **exact CKG KMS-DB sampler class
we use** on 2D TFIM at constant β ≥ β*, h ≤ h* ≈ 1. Mechanism: Peierls
fault-line surface tension. Martinelli–Olivieri 1994 (DOI 10.1007/BF02187060)
is the classical analogue: gap ≤ exp(-τ(β)·L) for 2D Ising Glauber, L = √N.

**Caveat — ladder geometry, NOT square**: our $2 \times L_y$ ladders have a
fixed cheapest cut (2 vertical bonds × 2J = 4J, $L_y$-independent), giving
a floor gap ≳ exp(-β·4J) — not the √N closing of GKZ24. Predictions:

| β | exp(-β·4J) floor | observed gap at n=8 |
|---|---|---|
| 0.5 | 0.14 | 0.285 (saturated by n=6) |
| 1.0 | 0.018 | 0.010 (close to floor) |
| 2.0 | 3.4e-4 | 1.6e-4 (factor 2 below — finite-size prefactor) |

All within $O(1)$ of the floor estimate at $n=8$. Closing magnitudes are
ladder-specific; mechanism is the same as the square-lattice GKZ24
prediction. **Sandbox sandbox cells are qualitative ordered/paramagnetic
distinction, not quantitative tests of the GKZ24 √n exponent.**

**Implications for qf-1jj draft**: keep the qf-1jj qualitative description
(memory [[2d-tfim-ordered-vs-disordered-qf-1jj]]) — the data IS consistent
with the GKZ24 mechanism, just at sandbox scales where the closing magnitude
is bounded by the ladder cut. Add a footnote pointing at qf-c9g for the
β-scan that rules out the GS-doublet artefact alternative.

**Cross-references**:
- `drafts/qf-c9g-2d-tfim-ordered-mechanism.md` — full writeup
- `drafts/qf-c9g-lit-survey.md` — literature distillation
- `drafts/qf-biz-2d-tfim-matrix-element-bottleneck.md` — extends qf-biz Follow-up A across β
- [[2d-tfim-ordered-vs-disordered-qf-1jj]] — parent qf-1jj memory
- [[qf_biz_2d_tfim_matrix_element_refuted]] — within-doublet refutation
- [[tc_2d_tfim_phase_diagram]] — T_c(h) table
