---
name: qf-e4z-31-v6-plus-parity-trap-quantified
description: qf-e4z.31 PARTIAL re-run with rho_0=|+⟩⟨+|^⊗N + krylovdim=60 settled the v6 (I/d, kdim=40) parity trap; even/odd splitting in 1D Heisenberg PBC is REAL physics but v6 had inflated odd-n gap by 14-21%
metadata:
  type: project
---

## qf-e4z.31 v6_plus sweep: rho_0=|+⟩⟨+|^⊗N + krylovdim=60 (n=3..7, 150 cells, 2026-05-16)

Re-ran the qf-e4z.23 multiseed Heisenberg PBC sweep on the SAME 30 fixtures
(scripts/output/multiseed_fixtures/heis_xxx_disordered_periodic_n{3..7}_seed{42..46}.bson)
with the qf-e4z.30 canonical recipe: `rho_0 = |+⟩⟨+|^⊗N` + krylovdim=60 +
single-pass via `predict_lindbladian_trajectory`. Saved Pass-1 gap
(`eigenvalues[2]` from rho_0-seeded Arnoldi) AND Pass-2 gap (post-qf-e4z.27
patch via `krylov_spectral_gap` with `_krylov_default_x0` + KrylovKit
thick-restart at kdim=30) per cell as an internal cross-check.

**Why:** qf-e4z.23 reported a striking even/odd-n splitting (even-n gap_phys
collapsing with β, odd-n flat). v6 had used `I/d` Arnoldi at kdim=40, which is
parity-trapped on Z+ZZ-disordered fixtures (P = Z^⊗N commutes with both H and
rho_0). Needed to know whether the splitting was real physics or an I/d artefact.

**How to apply:** Cite the v6_plus sidecars at
`scripts/output/sweep_S1_v6_plus_ckg_ideal_multiseed/smooth_metro_eps1e-03/`
(150 BSONs) and figure `drafts/figures/numerics/v6_plus_vs_v6_gap_phys.{png,pdf}`.
Driver `scripts/scratch_p1_v6_multiseed_plus.jl`; analysis `scripts/analyze_v6_plus_vs_v6.jl`.
The summary BSON is `scripts/output/sweep_v6_plus_multiseed_summary.bson`.

### Findings

(1) **Parity-trap correction is real and quantified** (v6_plus / v6 gap_phys ratios):

|     n     | β=0.25 | β=0.5 | β=1.0 | β=1.5 | β=2.0 | β=2.5 |
|-----------|--------|-------|-------|-------|-------|-------|
| 3 (odd)   | 1.000  | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| 4 (even)  | 1.000  | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| 5 (odd)   | 1.000  | 0.931 | 0.854 | 0.858 | 0.867 | 0.874 |
| 6 (even)  | 1.000  | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| 7 (odd)   | 1.000  | 1.000 | 0.791 | 0.798 | 0.811 | 0.823 |

→ v6 overestimated odd-n gap (n=5: 14-15%, n=7: 18-21%) at β ≥ 1 via I/d parity
trap. Even-n and n=3 unaffected. n=5,7 at β=0.25 also unaffected (low β: parity
sectors not yet split). Matches qf-e4z.28 spot measurements (n=5: 16.6%, n=7: 24.3%).

(2) **Even/odd splitting SURVIVES the fix** (v6_plus gap_phys, TRUE L gap):

|     n     | β=0.25 | β=2.5 | min/max | regime |
|-----------|--------|-------|---------|--------|
| 3 (odd)   | 3.79   | 3.57  | 0.942   | flat |
| 4 (even)  | 4.88   | 4.14  | 0.848   | modest decrease |
| 5 (odd)   | 4.33   | 4.04  | 0.932   | flat |
| 6 (even)  | 4.58   | 1.27  | **0.278** | **monotone collapse** |
| 7 (odd)   | 4.38   | 4.12  | 0.939   | flat |

→ Even-n monotonely collapses with β (especially n=6: 3.6× drop across the grid);
odd-n stays roughly flat (~7% variation). Inter-parity contrast at β=2.5: n=6=1.27
vs n=7=4.12 → ratio **3.2×** (down from v6's 3.6× over-estimate, but still strongly
inter-parity-resolved against ~1% intra-parity seed noise from 5 fixtures).

(3) **Pass1↔Pass2 internal cross-check** at v6_plus kdim=60 + |+⟩⟨+|^⊗N:
- n=3-4: 1e-14 (machine precision, all β)
- n=5: 5e-8 at β=0.25, 1e-11 at β≥1
- n=6: 5e-7 at β=0.25, 1e-12 at β≥0.5
- **n=7: 4.5e-4 at β=0.5, 2e-9 at β≥1.0** — low-β floor visible
→ At n=7 low β the spectrum has a wide continuum of slow modes that kdim=60
can't resolve to 1e-9. Still 4-decimal accurate (well below thesis ε=1e-3),
but qf-e4z.33 will harden to <1e-6 via kdim=80, 100 sweep + dense ref at n≤6.

### Policy decisions

POLICY (set 2026-05-16): **ALL future Heisenberg multiseed runs use
`rho_0 = |+⟩⟨+|^⊗N`** (qf-e4z.30 canonical). I/d is forbidden — wastes
compute reporting parity-trapped gaps. See [[canonical-taumix-setup-qf-e4z-30]].

Outstanding work:
- [[qf-e4z-32]] n=8 follow-up using same |+⟩⟨+|^⊗N + kdim=60 setup
- [[qf-e4z-33]] dense-ref validation at n≤6 + Pass-1 self-saturation at n=7 low-β
  (Pass 1 at kdim ∈ {60, 80, 100}, Pass 2 at kdim ∈ {30, 60, 80}).
  Cap at kdim=100; if not converged investigate root cause (algorithm bias,
  Krylov-resistant spectrum, or matvec floor).

### Caveats (per [[feedback-more-data-points-for-scaling-claims]])

3 even-n (n=4, 6) + 3 odd-n (n=3, 5, 7) is on the edge for parity-resolved
claims. The qualitative even-n β-collapse / odd-n flat structure is robust
(inter-parity log-spread 0.49→0.91 across β >> intra-parity scatter ~1%),
but extracting an asymptotic scaling exponent (e.g. `gap^even ~ exp(-γ·N)`)
from 3 even-n points is NOT defensible. n=8 from [[qf-e4z-32]] will give 3
even-n at one more system size, still in the qualitative regime.

The n=6 β=2.5 collapse to gap_phys=1.27 has 3-way Krylov agreement (v6,
v6_plus Pass 1, v6_plus Pass 2 — three different algorithms / kdims) to 4
decimals. Pending dense `eigvals` cross-check from [[qf-e4z-33]] part (A)
for the 4th-method confirmation.
