---
name: qf-e4z-33-kdim-convergence
description: kdim=60 was the bias source for v6_plus n=7 low-β 1e-4 Pass1/Pass2 disagreement; bump Pass 1 to kdim=80 for n≥7. krylovdim_gap_pass kwarg added.
metadata:
  type: project
---

# qf-e4z.33 kdim convergence — Pass 1 plain Arnoldi was the bias source (2026-05-17)

**Result (single-seed=42, n ∈ {3..7}, β_phys ∈ {0.25,…,2.5}, EnergyDomain CKG, `rho_0 = |+⟩⟨+|^⊗N`, dense ground-truth at n ≤ 6):**

The v6_plus `gap_arnoldi_pass1` ↔ `gap_arnoldi_pass2` disagreement of 1.80e-4 (n=7 β=0.25) and 9.08e-5 (n=7 β=0.5) is **fully** explained by Pass 1 plain-Arnoldi under-convergence at kdim=60. Pass 2 (KrylovKit thick-restart Krylov-Schur from `_krylov_default_x0`) was already at machine precision at its hardcoded kdim=30 — its number is the correct gap. Bumping Pass 1 to kdim=80 brings n=7 β=0.25 self-saturation under 1e-7; kdim=100 to 1e-11. Cross-agreement at kdim=100 reaches 2.16e-11 (β=0.25) / 1.07e-10 (β=0.5), down from 1.8e-4 / 9.1e-5 in v6_plus.

**Dense `construct_lindbladian + eigvals` cross-check (n ≤ 6, 24 cells)**: zero failures. Both passes match dense gap to ≤ 6e-12 (most cells 1e-14) — the qf-e4z.30 `|+⟩⟨+|^⊗N` recipe really does eliminate parity-trap artefacts at every n we can dense-check. The 1e-4 effect at n=7 was Pass 1 under-convergence on the plain-Arnoldi side, NOT a residual parity-sub-spectrum bias.

**Why:** plain Arnoldi from `vec(|+⟩⟨+|^⊗N)` captures only the top `kdim` modes by Ritz-value magnitude, with no restart. At low β (β_alg ~ 10–20), the Lindbladian spectrum at n=7 is dense near the slowest mode, so kdim=60 misses by 1e-4. KrylovKit's thick-restart amplifies the slowest modes regardless of where they sit, so kdim=30 + restart converges with ~120 effective matvecs.

**Canonical recipe update (replaces [[canonical-taumix-setup-qf-e4z-30]] kdim=60):**

| n | Pass 1 kdim | Pass 2 kdim |
|---|---|---|
| ≤ 6 | 40 | 30 |
| 7, β_phys ≥ 1.0 | 60 | 30 |
| 7, β_phys ≤ 0.5 | **≥ 80** | 30 |
| ≥ 8 | 80–100 (verify per cell) | 30 |

For production sweep drivers: use `krylovdim = 80, krylovdim_gap_pass = 30` as the new canonical default for any sweep touching n ≥ 7. Costs ~30% more matvecs than v6_plus default.

**τ_mix impact on qf-e4z.31 sidecars: negligible.** The 1e-4 gap shift at n=7 β=0.25 propagates to ≈ 0.02 % in τ_mix (47.59 → 47.598). Below the 1 % thesis-precision target — sidecars do NOT need regeneration.

**Src change:** `src/lindblad_action.jl` commit `581e067` adds `krylovdim_gap_pass::Union{Nothing, Integer} = nothing` kwarg on `predict_lindbladian_trajectory` / `predict_channel_trajectory`. Default `nothing` preserves `max(30, krylovdim÷2)` policy. Drivers iterating Pass 1 and Pass 2 at varying kdims should call `_krylov_spectral_decomposition` / `krylov_spectral_gap` directly (as qf-e4z.33 driver does) — the trajectory predictor re-runs Pass 1 each time Pass 2 changes.

**Open:** [[qf-e4z-32]] n=8 sweep deferred — applying the new kdim=80 / 30 recipe is recommended; verify Pass 1 self-saturation at n=8 once the sweep runs. Channel-side `predict_channel_trajectory` has identical Pass 2 wiring + kwarg; not kdim-swept yet.

**Driver:** `scripts/scratch_qf_e4z_33_kdim_sweep.jl`. Analyzer: `scripts/analyze_qf_e4z_33_kdim_sweep.jl`. Sidecars: `scripts/output/sweep_qf_e4z_33_kdim_sweep/seed42/`. Draft: `drafts/qf-e4z-33-kdim-convergence.md`. Analysis text: `scripts/output/qf_e4z_33_analysis.txt`.
