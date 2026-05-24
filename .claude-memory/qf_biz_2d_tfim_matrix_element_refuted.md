---
name: qf-biz-2d-tfim-matrix-element-refuted
description: Within-doublet matrix elements are NOT the bottleneck for 2D TFIM ORDERED L gap closing; slow mode lives 89% in the doublet-to-bulk coherence block (2026-05-24).
metadata: 
  node_type: memory
  type: project
  originSessionId: 2e0ec557-e0b4-44c1-92de-7a08e13be5d7
---

The "matrix-element bottleneck" story for the 2D TFIM ORDERED-phase L gap closing (qf-1jj, gap_phys = 9.96e-2 → 3.06e-3 → 1.57e-4 across n=4,6,8 at h=J=1, β_phys=2) — which I floated 2026-05-24 in response to user skepticism about a "Gibbs-on-GS" artefact — **is wrong**. Refuted in qf-biz (`drafts/qf-biz-2d-tfim-matrix-element-bottleneck.md`).

**Why:** M²_norm = Σ_a |⟨ψ_2|A_a|ψ_1⟩|² (1/√(3n)-normalised canonical jump set) is **constant** at 0.311 from n=4 to n=6, while λ_L^alg ratio is 0.023 (43× closing). Sum-rule: σ_z dominates (~1 per site by FM-minimum dressed-state algebra, σ_y small, σ_x = 0 by Z₂), summed over n sites gives M²_bare ∝ n, but the 1/(3n) normalisation kills the n-growth.

**Actual mechanism (Follow-up A, dense L at n=4):** the slow mode R_2 is **not within-doublet** — 89% of ‖R_2‖² lives in the doublet × bulk off-diagonal block (largest entries are |ψ_{1,2}⟩⟨ψ_8|-type coherences in the energy eigenbasis), only 2.8% in doublet × doublet. R_2 is Z₂-mixed-parity (tr(R_2† P R_2 P)/‖R_2‖² = 0.000), magnetisation-orthogonal (|⟨R_2, M_z⟩| = 0.000), and has small (0.119) overlap with the symmetric within-doublet coherence. The closing mechanism involves the full Davies-rate-equation eigenvalue over doublet-to-bulk coherences (energy-increasing transitions ω ≈ 5.9 phys → β·ω ≈ 12 → Boltzmann-suppressed), not bare ⟨ψ_2|A_a|ψ_1⟩.

**Independent finding (Follow-up B):** λ_L^phys at fixed n=6 is **non-monotone** in β_phys — minimum at β=1.5, rises as we go colder (5.57e-3 at β=3). Kinky-Metropolis spot check (s=0, a=0) matches the smooth-Metro shape to 4% at every β — the non-monotonicity is **physical**, not a filter artefact. Consistent with a slow-mode crossover (thermal-Glauber-barrier-crossing → kinetic-tunnelling) of the kind documented for finite-T quantum-annealing problems on exponentially-small-gap Hamiltonians [Bando et al. 2018, PMC6060131].

**Phase confirmation (Check 1):** β_phys=2 IS deep in the ordered phase (eff_rank=2.00, ⟨M_z²⟩=0.945, U_4=0.656). Crossover at β≈0.5 = β_c per HW2016 T_c(h=1)≈2.07. Distinct from [[heisenberg-1d-even-odd-mechanism]] (1D Heis even-n: eff_rank ≈ 1, p_0 ≈ 99%, gap drifts mildly with n — the artefact regime). 2D TFIM ORDERED is the genuine thermal-order regime.

**qf-1jj draft (`drafts/2d-tfim-ordered-vs-disordered.md`) implications:** no headline edits needed — physics-checker confirmed no inline matrix-element claim exists. Add a 2-sentence footnote to "Methodology notes" pointing to qf-biz for the refuted bottleneck and the doublet-to-bulk slow-mode finding.

**Driver scripts:** `scripts/scratch_qf_biz_phase_and_matrix_elements.jl` (Checks 1+2, 60s) and `scripts/scratch_qf_biz_followup.jl` (Follow-ups A+B, 34s). Sidecars in `scripts/output/qf_biz_phase_and_matrix_elements/`. Verifier reports: `drafts/qf-biz-{code-verifier,physics-checker}-report.md`.

Open question for a future session if needed: precise prefactor and n-scaling of the Davies-rate-equation eigenvalue for the dominant doublet × bulk coherence mode. Would require dense L eigenmode at n=6 (d²=4096, ~1 min) and possibly n=8 (d²=65536, heavy) to see how the doublet × bulk block mass evolves with n.
