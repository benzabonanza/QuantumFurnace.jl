---
name: qf-e4z-34-ckg-vs-dll-plot-sweep
description: "qf-e4z.34 plot-grade CKG-vs-DLL sweep 1D Heisenberg PBC n=3..8 β_phys∈{0.25,0.5,1.0} 5 seeds completed 2026-05-17 — CKG faster than DLL uniformly, even-n gap does not collapse at n=8 on β≤1"
metadata: 
  node_type: memory
  type: project
  originSessionId: a66686ad-a81e-4448-bb55-c211572e55de
---

180/180 cells (90 CKG = KMS+Energy+smooth-Metro s=0.25, 90 DLL = DLL+Bohr+DLLMetropolisFilter), each with 5 disorder seeds, single-pass `predict_lindbladian_trajectory` from `rho_0 = |+⟩⟨+|^⊗N` ([[canonical-taumix-setup-qf-e4z-30]]). Flat kdim=80 at n=7,8; kdim=60 at n≤6. ε=1e-3, t_max=500.

Result tables in `scripts/output/sweep_qf_e4z_34_ckg_vs_dll_plot_summary{,_stats}.bson` (the latter has p25/median/p75 across seeds for the 36 (sampler, n, β_phys) cells). all_converged=true everywhere, flagged_cells=0, mv/cell at n=8 is 168-216 (CKG) or 199-216 (DLL).

**Why:** First plot-grade pipeline run combining all qf-e4z resolutions ([[predict-trajectory-two-pass-qf-e4z-27]] → [[canonical-taumix-setup-qf-e4z-30]] for parity-broken seed, [[qf-e4z-33-kdim-convergence]] for kdim recipe, [[feedback-thesis-numerics-grid-canonical]] for the β grid). Single canonical dataset for the CKG-vs-DLL comparison figure (Figure A: gap_phys; Figure B: τ_mix_phys; both 2 cols × 3 rows = even/odd × β).

**How to apply:**
- Plotting is deferred — `scripts/plot_qf_e4z_34_ckg_vs_dll.jl` is the next step; loads the `_summary_stats.bson` and produces the IQR-band figures.
- For any extension (n ≥ 9 cluster run, alternative β cells, new sampler), re-use the canonical pipeline parameters from `scripts/scratch_qf_e4z_34_ckg_vs_dll_plot_sweep.jl` lines 70-96 (kdim recipe, β grid, register layout, |+⟩⟨+|^⊗N seed).
- **CKG/DLL ratio** (canonical regime): CKG τ_mix is **smaller** at every cell (ratio CKG/DLL ∈ [0.79, 1.00], median ≈ 0.89). CKG advantage is strongest at low β (~0.80 at β=0.25, ~0.94 at β=1.0). The smooth-Metro kink kernel is genuinely faster than DLL Metropolis on this grid.
- **gap_phys**: CKG gap is **larger** at most cells (ratio 1.11-1.72), modulo 3 cells where odd-n high-β has DLL slightly ahead (n=3,5,7 at β=1.0, ratios 0.95-1.00). Gap and τ_mix do not move together — the τ_mix advantage is mostly the gap advantage, but for the odd-n β=1 cells τ_mix still favours CKG even though the gap is essentially tied.
- **Even-n n=8 gap_phys NOT collapsing** on this grid: CKG even-n at β=1 has gap_phys 3.44 (vs n=6 → 3.52, n=4 → 4.56) — modest decline, no exp(-cn) collapse. The dramatic even-n collapse from [[heisenberg-1d-multiseed-even-odd-qf-e4z-23]] / [[heisenberg-1d-even-odd-mechanism]] is a β_phys ≥ 1.5 effect (Gibbs concentrates on GS, slowest mode becomes GS↔first-excited coherence). At β ≤ 1 the Gibbs is still spread enough that the mechanism does not kick in. The thesis plot at the canonical β grid will therefore show "qualitative even/odd split with CKG and DLL both staying Ω(1)", not the dramatic collapse story — consistent with [[feedback-numerics-grid-canonical]].
- All n=8 cells report `mixing_time_source = :extrapolated`. **The name is misleading**: `eigenmode_mixing_time` (src/mixing.jl:536) does an HONEST BISECTION on `d(t) − ε = 0` using the closed-form eigenmode formula; `:extrapolated` means "bisection found a crossing", NOT "fit-extrapolated past data". `τ_mix` is genuinely `τ : d(τ) = 1e-3` to `atol=1e-3` in t. `floor_distance ≈ 1e-15 ≪ ε` at every cell — the asymptotic floor is irrelevant.

**Post-hoc analysis (2026-05-17 same session):**

1. **HS-norm / d_{1→1} diagnostic** (`scripts/scratch_qf_e4z_34_norm_diagnostic.jl`, output `scripts/output/qf_e4z_34_norm_diagnostic.bson`): at every (n, β_phys, seed=42) cell n ∈ {3..7}:
   - `‖L_CKG‖_HS / ‖L_DLL‖_HS = 0.99–1.01` (within ±0.4%)
   - `‖L_diss^CKG‖_{1→1} / ‖L_diss^DLL‖_{1→1} = 0.99–1.03` (within ±1%)
   - So the CKG/DLL speedup is **~100% structural** (a wider gap at the same total operator scale), NOT rate-scale. Same total rates, better spectral structure. β_phys-aligned comparison (vs the older β_alg-aligned comparison) is what makes this clean — at the same physical temperature both samplers have the same total transition-rate magnitude.

2. **Gap accuracy verified** (`scripts/scratch_qf_e4z_34_verify_b05.jl`): independent `krylov_spectral_gap` (KrylovKit thick restart from `_krylov_default_x0`, kdim=80, tol=1e-12) reproduces sidecar `gap_phys` to **rel_diff ≤ 1.2e-7** at every β_phys=0.5 cell n ∈ {4, 6, 7, 8} (typical 1e-11). Sidecar Pass-1 gap (from `traj.eigenvalues[2]` aligned with `|+⟩⟨+|^⊗N`) is reliable.

3. **β_phys=0.5 trend revisited** (was suspicious — gap ratio widens with n while τ ratio stays flat):
   - At β_phys=0.5 the gap ratio λ_C/λ_D widens with n (even: 1.43 → 1.56 → 1.72 from n=4 to n=8; odd: 1.12 → 1.41 → 1.62 from n=3 to n=7) — confirmed by independent cross-check, NOT a Krylov artefact.
   - Observed τ_mix-speedup stays flat at ~1.10–1.20 across the same n range.
   - **Resolution from τ·λ diagnostic** (alg frame at β_phys=0.5): for CKG `τ·λ ≈ 7` and `exp(τ·λ)·ε ≈ |c₂| ≈ 1` (slow mode dominates the decay); for DLL `τ·λ ≈ 5` and `|c₂| ≈ 0.1` (slow mode barely populated by `|+⟩⟨+|^⊗N`). So DLL's small gap doesn't hurt its mixing as much as it "should" — faster modes carry the bulk of the trace distance until τ_mix, and the slow mode hasn't even decayed by 1 e-folding at the crossing.
   - **Interpretation**: the gap ratio is the **asymptotic / structural** speedup potential; the τ_mix ratio is what's REALIZED for the specific initial state. They are not the same when `|c₂|` differs strongly between samplers, as it does here. For thesis presentation, report **τ_mix from `|+⟩⟨+|^⊗N`** as the observed mixing performance; report **gap_phys** as a spectral diagnostic — do NOT call the gap ratio "the speedup". The gap-widens-with-n story at β_phys=0.5 is real but applies to the *asymptotic* mixing timescale, not the ε=1e-3 crossing time.

**How to apply (added 2026-05-17):**
- When presenting the comparison: **τ_mix ratio τ_D/τ_C ∈ [0.97, 1.28] (median ~1.15)** is the observed plot story; the gap ratio is a structural diagnostic but overstates the effective speedup at high n / mid β because `|c₂|_DLL ≪ |c₂|_CKG` from `|+⟩⟨+|^⊗N`.
- For "where does CKG help more" claims: lean on **β_phys-axis** (CKG advantage decays from ~25% at β=0.25 to ~7% at β=1.0) — robust, parity-stable, monotone. The system-size axis at this grid is essentially flat in τ.
- For `bd note qf-e4z.34`-style follow-ups: never quote "gap ratio widens with n at β=0.5 → CKG advantage grows with N" without the `|c₂|_DLL` caveat. Initially I made that claim and had to retract it within an hour.

**Cluster follow-ups parked under [[qf-e4z-32]] / [[qf-e4z-23]] cluster notes:** n ∈ {10,12,14} on the same canonical grid, both PBC and OBC, same pipeline parameters. Re-verify kdim=80 saturation at n=10+ first ([[canonical-taumix-setup-qf-e4z-30]] flagged this as still-to-confirm).
