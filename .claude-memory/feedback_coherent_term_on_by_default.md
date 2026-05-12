---
name: feedback-coherent-term-on-by-default
description: "Simulations must keep include_coherent=true by default. Dissipator-only runs (include_coherent=false) must be called out explicitly in headers, outputs, and summaries — they do NOT size coherent registers and the resulting operator is non-physical."
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 14338e97-2944-44d2-96b9-d88a1f6145c9
---

# Coherent term ON by default — flag any dissipator-only run

When calling `apply_lindbladian!`, `apply_adjoint_lindbladian!`,
`construct_lindbladian`, `predict_lindbladian_trajectory`,
`lindblad_action_integrate`, `krylov_spectral_gap`, or any other
Lindbladian-facing API: **keep the default `include_coherent = true`**.
The default IS the physical KMS Lindbladian. Setting it to `false`
yields a **dissipator-only** operator that does not satisfy KMS detailed
balance — its spectrum / fixed point are not the physical ones.

## Why

User asked me to extend the qf-yt9 quadrature sweep to β_phys=2 and then
to TimeDomain at n=5. I built both campaigns with
`include_coherent = false` (copying the qf-yt9 dissipator-only S1 setup)
but never flagged in the markdown summary or the side-by-side tables
that those numbers were dissipator-only. The user spotted it and asked
"so you chose r_bp r_bm as well and had B created in the TimeDomain and
everything?" — confirming that the framing I produced was misleading.
Reporting "the Lindbladian matches at r_D = 7" when the matvec was
dissipator-only conflates two register-sizing questions (the dissipator
quadrature vs the coherent-term quadrature) that decouple completely:
the TimeDomain coherent term caps at ε ≈ 10⁻⁶ from the t = 0
L'Hôpital sample, so the full-Lindbladian register triple is *not*
`r_D = r_b_minus = r_b_plus = 7`.

## How to apply

- Default in any new script: do not pass `include_coherent` at all —
  use the canonical default.
- When `include_coherent = false` is genuinely needed (isolating
  dissipator quadrature, debugging R_total assembly, etc.):
  - Say so in the script header docstring.
  - Name the script and output directories accordingly (e.g.
    `scratch_quad_S1_dissipator_*` rather than `scratch_quad_*`).
  - Open the markdown summary with a one-line caveat block: "These
    numbers measure the dissipator term only — the coherent commutator
    $i[B, \rho]$ is masked. Full-Lindbladian register sizing requires a
    separate coherent-term sweep."
  - Never present a dissipator-only `r_D` cutoff as the full-Lindbladian
    register recipe.
- B_time / B_bohr / B_trotter is still BUILT in the workspace even with
  `include_coherent = false` (it lives in `G_left`/`G_right` and
  cancels in the matvec). So the per-matvec wall is not a clean
  dissipator-only measurement — it includes B construction. Adjust wall-
  time expectations accordingly when extrapolating to n + 1.
- For Krylov spectral-gap / fixed-point work the coherent term is
  spectrum-shifting and Hermitian → contributes anti-Hermitian
  $-iB^\top$. Never drop it to "simplify the eigenproblem"; the gap and
  steady state change.

Rule lives in `.claude/rules/julia-code.md` ("Coherent Term: ON by
Default") and `.claude/rules/scripts.md`.

## Pointers

- API surfaces: `apply_lindbladian!` etc. all have
  `include_coherent::Bool = true` keyword.
- Workspace internal: `G_left = +iB^\top - R/2`,
  `G_right = -iB^\top - R/2` (src/krylov_workspace.jl:49-50); with
  `include_coherent=false` matvec uses
  `(G_left + G_right)/2 = -R/2`, dropping the unitary part.
- Related: [[quadrature_register_recipe_v2]] — TimeDomain coherent ε
  cap at slope-(−1) from t=0 L'Hôpital sample, see
  [[trap_rule_t0_lhopital_origin]].
