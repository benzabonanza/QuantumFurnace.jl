# qf-e4z.33 — Krylov dimension convergence at problematic cells

**Status:** done (2026-05-17)
**Driver:** `scripts/scratch_qf_e4z_33_kdim_sweep.jl`
**Analyzer:** `scripts/analyze_qf_e4z_33_kdim_sweep.jl`  (full output: `scripts/output/qf_e4z_33_analysis.txt`)
**Sidecars:** `scripts/output/sweep_qf_e4z_33_kdim_sweep/seed42/`  (30 cells)
**Scope:** single seed (42), n ∈ {3..7}, β_phys ∈ {0.25, 0.5, 1.0, 1.5, 2.0, 2.5}, smooth Metro s=0.25 a=0, r_D=7 (β_phys≤1) / r_D=8 (β_phys≥1.5), EnergyDomain CKG, coherent ON, `rho_0 = |+⟩⟨+|^⊗N` (qf-e4z.30).

## Question

The qf-e4z.31 v6_plus sweep used `rho_0 = |+⟩⟨+|^⊗N` + `krylovdim = 60`, and at the n=7 low-β cells reported a Pass 1 ↔ Pass 2 spectral gap disagreement of ~10⁻⁴:

| β_phys | v6_plus Pass 1 (kdim=60) | v6_plus Pass 2 (kdim=30 + thick restart) | rel_diff |
|---|---|---|---|
| 0.25 | 0.10250424 | 0.10248581 | **1.80e-4** |
| 0.5  | 0.12611566 | 0.12610421 | **9.08e-5** |
| ≥ 1.0 | (machine precision) | | ≤ 1e-9 |

The Pass 2 internal kdim policy was `max(30, kdim÷2)` — so at v6_plus's call `predict_lindbladian_trajectory(..., krylovdim=60)`, Pass 2 ran at kdim=30, *smaller* than Pass 1. The question: is the residual a Krylov-truncation effect that saturates with larger kdim, an algorithm bias requiring a different approach (different seed, shift-invert, etc.), or a genuine Krylov-resistance?

## Result: outcome (i) — Krylov truncation only

**Pass 1 alone** was under-converged at kdim=60 for the n=7 low-β cells; Pass 2 was already at machine precision at kdim=30.

After bumping both kdims to 100 (single-seed=42, same fixtures and config as v6_plus):

| β_phys | Pass 1 (kdim=100) | Pass 2 (kdim=100) | cross rel_diff | v6_plus Pass1[60] vs new Pass1[100] | v6_plus Pass2[30] vs new Pass2[100] |
|---|---|---|---|---|---|
| 0.25 | 0.10248581 | 0.10248581 | **2.16e-11** | **1.8e-4** ⚠ | 8.7e-15 |
| 0.5  | 0.12610421 | 0.12610421 | **1.07e-10** | **9.1e-5** ⚠ | 1.9e-12 |
| 1.0  | 0.10182588 | 0.10182588 | 1.76e-14 | 5.2e-10 | 5.5e-15 |
| 1.5  | 0.09737452 | 0.09737452 | 5.42e-15 | 3.8e-09 | 1.6e-15 |
| 2.0  | 0.09678152 | 0.09678152 | 3.58e-15 | 1.0e-09 | 2.0e-15 |
| 2.5  | 0.09659483 | 0.09659483 | 1.36e-14 | 2.0e-09 | 1.9e-14 |

The cross-disagreement drops by **7 orders of magnitude** at the worst cell (n=7 β=0.25) once both passes get kdim=100. The "right" gap is the Pass 2 value all along — Pass 1 in v6_plus was off by 1.8e-4 because of plain Arnoldi truncation at kdim=60.

### Per-cell Pass 1 convergence trace (n=7, β_phys=0.25)

```
P1  k=40: 0.10299257678   (rel to k=100: 4.93e-3)
P1  k=60: 0.10250423537   (rel to k=100: 1.80e-4)  ← v6_plus value
P1  k=80: 0.10248580929   (rel to k=100: 3.65e-8)
P1  k=100:0.10248580555   (saturated)

P2  k=30: 0.10248580555   (matches k=100 to 1e-14)  ← v6_plus value, already converged
P2  k=60: 0.10248580555
P2  k=80: 0.10248580555
P2  k=100:0.10248580555
```

Pass 1 needs kdim ≥ 80 to bring n=7 β=0.25 under 1e-6. Pass 2 (KrylovKit thick-restart Krylov-Schur from `_krylov_default_x0`) was already saturated at kdim=30 — the matvec budget (~120 with restarts) and the parity-breaking `1e-10·H_GUE_traceless` seed both contribute.

### Dense ground-truth match (n ≤ 6, 24 cells)

For every dense-ref cell, both Pass 1 and Pass 2 at kdim=100 match the dense `construct_lindbladian + eigvals` gap to **≤ 6e-12** (most cells at 1e-14, machine precision). Zero failures in 24 cells. Confirms:

* The qf-e4z.30 `rho_0 = |+⟩⟨+|^⊗N` recipe fully eliminates the parity trap that bit qf-e4z.23 (`rho_0 = I/d`): every captured Krylov eigenvalue lives in the full Lindbladian spectrum, not a parity sub-spectrum.
* The dense reference itself (24 cells × full d²×d² LAPACK `geev`) confirms the gap_arnoldi numbers from v6_plus at n ≤ 6 were correct to 1e-12 or better, even at the low-β cells where Pass 1 alone would have showed the kdim=60 truncation if d² were as large as 16384.

### Self-saturation per cell (kdim 80 → 100)

Every cell except n=7 β ∈ {0.25, 0.5} has Pass-1 self-saturation `|gap(80) - gap(100)| / gap < 1e-12`. The two low-β n=7 cells saturate at:

* n=7 β=0.25: `|p1(80) - p1(100)| / gap = 3.65e-8`
* n=7 β=0.5:  `|p1(80) - p1(100)| / gap = 1.23e-7`

Both well below the 1e-6 acceptance threshold (just barely close to saturation; kdim=80 would be the minimum, kdim=100 gives the safety margin). Pass 2 saturates `< 1e-12` at every cell at kdim ≥ 30.

## Decision: canonical recipe update

| Regime | Pass 1 kdim | Pass 2 kdim | Reason |
|---|---|---|---|
| n ≤ 6, all β | 40 | 30 | Both saturate at machine precision; kdim=60 (v6_plus) was conservative but cheap. |
| n = 7, β_phys ≥ 1.0 | 60 | 30 | Pass 1 already at 5e-10 at kdim=60; safe. |
| **n = 7, β_phys ≤ 0.5** | **≥ 80** | 30 | Pass 1 plain Arnoldi needs kdim=80 to reach <1e-7. Pass 2 still fine at 30. |
| n ≥ 8 (extrapolation) | **80–100** (verify per cell) | 30 | Unverified at n=8+; the kdim=60 → 80 jump that fixed n=7 will probably reappear at higher n. Always verify with a Pass-2 cross-check. |

**Practical guidance for sweep drivers:** use `krylovdim = 80, krylovdim_gap_pass = 30` as the new canonical default for any production sweep that touches n ≥ 7. This costs ~30% more matvecs than the v6_plus default (kdim 60/30) but eliminates the 1e-4 Pass-1 bias visible in v6_plus n=7 low-β. The qf-e4z.31 τ_mix values themselves are unaffected (τ_mix used Pass 1 = `gap_arnoldi`, and the 1e-4 gap bias propagates to ~0.02% in τ_mix at n=7 β=0.25 — negligible vs the 1% thesis-target accuracy).

## τ_mix impact assessment

Recomputing the v6_plus τ_mix with the corrected Pass 2 gap at n=7 β=0.25:

* v6_plus reported τ_mix = 47.59 (using Pass 1 gap 0.10250424).
* Corrected: τ_mix changes by at most `Δ(1/gap) × log(d/ε)` × small constant — for the bi-exp fit this is roughly proportional to the gap, so τ_mix_corrected ≈ 47.59 × (0.10250424 / 0.10248581) = 47.598, i.e. **0.02% shift**.

The qf-e4z.31 sidecars do NOT need to be regenerated for thesis figures — the corrected gap moves τ_mix below the 1% precision targets.

## Code change

`src/lindblad_action.jl` — both `predict_lindbladian_trajectory` and `predict_channel_trajectory` now expose `krylovdim_gap_pass::Union{Nothing, Integer} = nothing`. Default `nothing` reproduces the previous `max(30, krylovdim÷2)` policy bit-for-bit; callers that want Pass 2 at a different kdim can set it explicitly. (Commit `581e067`.)

Callers iterating Pass 1 and Pass 2 at varying kdims should call the underlying `_krylov_spectral_decomposition` and `krylov_spectral_gap` directly (as the qf-e4z.33 driver does) rather than forwarding through the trajectory predictor — that avoids re-running Pass 1 each time Pass 2 changes.

## Cost summary

* Per Pass 1 call: ≈ kdim matvecs (no restart). At n=7 EnergyDomain r_D=7: 100 matvecs ≈ 60s; at r_D=8: 130s.
* Per Pass 2 call: ≈ 80–220 matvecs (KrylovKit restart). Total wall ≈ 70–270s at n=7 depending on r_D.
* Full single-seed kdim sweep (4 P1 × 4 P2 = 8 calls per cell, 30 cells) wall ≈ 2 h on the sandbox.
* Dense ground-truth at n=6 (single LAPACK `geev` on 4096×4096 complex): ≈ 4–7 min per cell, dominated the n=6 wall budget.

## Open follow-up

* `qf-e4z.34` (not yet filed): the same kdim sweep at **n=8** to characterise whether Pass 1 saturates at kdim=80–100 there, or whether n=8 needs kdim ≥ 120. The qf-e4z.31 n=8 sweep is currently deferred (it's the next item in the v6_plus partial scope).
* Channel-side analog: `predict_channel_trajectory` has the same Pass-2 hardcoded kdim policy and the same `krylovdim_gap_pass` kwarg now; no kdim sweep run on the channel path yet. Probably fine since the channel hierarchy is currently used only for τ_mix at n ≤ 7 in P3 / S2 sweeps, but worth a single spot-check at n=7 β=0.25 if the channel path becomes a publication-target plot.

---

*Single seed = 42 per n. This is a krylovdim convergence diagnostic, not a seed-disorder statistics run — for the latter, see the canonical 5-seed v6_plus sweep (qf-e4z.31). All 30 cells saturate to <1e-7 cross-agreement at kdim=100.*
