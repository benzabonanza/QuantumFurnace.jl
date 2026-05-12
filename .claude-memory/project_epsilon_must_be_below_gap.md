---
name: project-epsilon-must-be-below-gap
description: "ε (total operator-norm error of L_discrete) must be smaller than the smallest spectral gap |λ_1| for precise Krylov gap / τ_mix estimates. Pick registers relative to the gap you want to resolve, not the other way around."
metadata: 
  node_type: memory
  type: project
  originSessionId: 14338e97-2944-44d2-96b9-d88a1f6145c9
---

# ε must be below the spectral gap for precise Krylov gap/τ_mix estimates

For any Krylov-based gap or τ_mix predictor
(`predict_lindbladian_trajectory`, `predict_channel_trajectory`,
`krylov_spectral_gap`) the relative error on the gap (and hence on
τ_mix = 1/|λ_1|) is bounded by

  δ|λ_1| / |λ_1| ≲ κ · ‖ΔL‖_op / |λ_1|

where ‖ΔL‖_op is the **total** operator-norm error between the discrete
Lindbladian we built and the continuum L_∞. The "ε" label on any
register table is shorthand for the **sum of all quadrature
contributions** (dissipator quadrature, coherent quadrature, FINUFFT
precision, Trotter M, δ-step on the channel side, etc.). The Krylov
factorisation itself reads the spectrum of whatever operator we built
to machine precision — all of the error is in the operator.

**Why:** Bauer-Fike (with κ caveat for non-normal Lindbladians):
eigenvalues shift by at most ‖ΔL‖_op in absolute terms. When the gap is
small (low T, large n) the same absolute shift becomes a large relative
error on τ_mix. Verified data point at n=6, β_phys=2, smooth_def_s:
TRUE gap (EnergyDomain, B exact, r_D=8) = 3.97e-3; TimeDomain
"ε=10⁻⁶" register gave 4.06e-3 (+2.2% off because ‖ΔB‖ floored at
2.4e-5); TimeDomain "ε=10⁻³" register gave 3.40e-3 (-14.5% off
because ‖ΔB‖ ≈ 1.15e-3 from undersized r_b_plus=7).

**How to apply:**

- Do not pick ε as a free parameter. Pick it relative to the smallest
  gap you intend to resolve: ε ≤ 10⁻² · min|λ_1| for 1% τ_mix
  accuracy, ε ≤ 10⁻¹ · min|λ_1| for 10%.
- For the *ideal Lindbladian τ_mix* plot (P1 / thesis S1): use
  `Config{Lindbladian, EnergyDomain, KMS}` with r_D ≥ 7. Coherent is
  automatically exact (B_energy ≡ B_bohr). Dissipator quadrature at
  r_D=7 reaches ~10⁻⁹ to 10⁻¹¹ at every (n=3..6, β_phys=0.25..2)
  cell — i.e. ‖ΔL‖ much smaller than any reasonable gap. So
  EnergyDomain r_D=7 is essentially gap-exact at our parameter range.
- For the *implemented-channel τ_mix* plot (compare CPTP δ-step
  channel vs ideal Lindbladian): the channel has its own (δ, GQSP,
  Trotter M) error budget which contributes to ‖ΔΦ_δ‖. Same rule —
  effective ε is the channel's intended error budget, must satisfy
  ε < gap for trustworthy τ_mix.
- At low T (β_phys ≥ 2, n ≥ 5) the gap drops to 10⁻³–10⁻⁴ range.
  TimeDomain has a β-amplified coherent-term floor here (2.4e-5 at
  n=6, β_phys=2 — scales like β_alg³ from the t=0 L'Hopital anomaly
  combined with smooth_def_s `s ∝ β_alg²` — see
  [[trap_rule_t0_lhopital_origin]]). **TimeDomain cannot resolve a
  gap of order 10⁻³ at low T to better than ~1% regardless of
  register sizing.** Use EnergyDomain.
- First pass for a new sweep: do a cheap probe at coarse registers,
  read off the minimum gap across the grid, then pick the production
  register sizing to satisfy ε ≤ 10⁻² × min(gap). Adaptive.
- The τ_mix you measure is the τ_mix of *whatever operator you
  built*, not of the continuum L_∞. They coincide only when ε ≪ gap.

## Pointers

- Detailed write-up with full numerical data: `drafts/error-analysis/epsilon-vs-gap-rule.md`
- v2 register recipe (the "ε" labels): [[quadrature_register_recipe_v2]]
- TimeDomain coherent t=0 floor origin: [[trap_rule_t0_lhopital_origin]]
- Krylov routes (Lindbladian / Channel, Energy / Time): [[krylov_two_routes]]
- include_coherent default: [[feedback_coherent_term_on_by_default]]
- Beads: qf-1au, qf-f45 (closed); epic qf-yt9 thread
