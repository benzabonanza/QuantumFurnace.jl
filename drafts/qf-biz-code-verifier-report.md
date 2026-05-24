# qf-biz code-verifier report

**Verdict: PASS.** The qf-biz diagnostic script is correct on every axis I checked. No bugs found. The published numbers (M²_norm = 0.311 at n=4 *and* n=6, L gap at n=6 β_phys=2 = 3.06e-3) reproduce exactly on a fresh run and match the qf-1jj sidecar **bit-for-bit**.

## Numerical agreement with qf-1jj

| field | qf-biz (re-run) | qf-1jj sidecar | rel-err |
|---|---|---|---|
| `gap_phys` (n=6, β_phys=2) | 3.0648095915e-03 | 3.0648095915e-03 | **0** (identical bytes) |
| `gap_alg`  (n=6, β_phys=2) | 6.0769489796e-05 | 6.0769489796e-05 | **0** |
| `gap_phys` (n=4, β_phys=2) | 9.9645e-02 | 9.9644545416e-02 | < 1e-4 (stdout rounding only) |
| `gap_alg`  (n=4, β_phys=2) | 2.6243e-03 | 2.6243376030e-03 | < 1e-4 (stdout rounding) |
| ⟨M_z²⟩, U_4, S/log d, eff_rank (n=6, β=2) | 0.9454, 0.6557, 0.1667, 2.0001 | identical (draft claims 0.945, 0.656, 0.167, 2.00) | match |

The exact agreement of `gap_phys` is expected: same `krylov_spectral_gap` call, same `_build_jump_set`, same CKG cfg (smooth-Metro s=0.25, a=0, r_D=7, kdim=40, howmany=4). The qf-biz Config builder reproduces qf-1jj's `t0_D=5.93523`, `w0_D=0.00827051`, `rescaling_factor=50.4334` exactly.

Total wall: 60.6 s (matches the draft's "60.9 s").

## Matrix-element computation — independently verified

I rebuilt `H_phys` for n=4 (2×2, h=1, J=1), did a fresh `eigen` → |ψ_1⟩, |ψ_2⟩, and computed `|⟨ψ_2|σ_α^(i)|ψ_1⟩|²` via *explicit* `kron`-based padding (independent of `pad_term`). Result matches the script's `pad_term` path to machine precision:

- `|X|²_sum_bare = 9.60e-31` (machine zero; symmetry-protected)
- `|Y|²_sum_bare = 4.69e-03`
- `|Z|²_sum_bare = 3.7267`
- `M²_norm = M²_bare / (3n) = 0.3109` (matches script's 0.3109; matches n=6 value 0.3112 to <0.1%)

Z₂ commutator structure verified explicitly: `‖[P, σ_x^(i)]‖ = 0`, `‖{P, σ_{y,z}^(i)}‖ = 0` on every site — so σ_x ME = 0 by symmetry (the draft's claim). Parity expectation values `⟨ψ_1|P|ψ_1⟩ = +1.0`, `⟨ψ_2|P|ψ_2⟩ = -1.0` to 1e-16.

Re. frame independence: `eigen(H_phys)` vs `eigen(H_alg)` give *identical* eigenvectors (`H_alg = H_phys/R + s·I` differs by scaling and uniform shift). The script does the matrix elements in the PHYS frame, which is the cleanest convention. `Δ1_phys = raw.eigvals[2:1] · R` reproduces `eigen(H_phys)` directly to 1e-13.

## Issues found

| # | Severity | Location | Issue |
|---|---|---|---|
| 1 | NOTE | `matrix_element_diagnostic`:355 | `eigvals_phys_low5 = eigvals_phys_sorted[1:min(5, end)]` uses `end` inside an indexing expression on a slice — works as written but reads as if `end` refers to the slice; clearer as `min(5, length(eigvals_phys_sorted))`. Cosmetic only. |
| 2 | NOTE | draft table line 14 (`Δ_1^phys` column for β_sweep) | Δ_1^phys is constant (≈ 1.07e-2) across all β rows because it depends only on H, not β. Worth a 1-line footnote in the draft saying "ΔE_1^phys is β-independent; column kept for cross-reference to qf-1jj." |
| 3 | NOTE | line 27 of draft | The "L gap is non-monotone in β" observation is real (script confirms minimum at β≈1.5: 1.24e-2 → 2.65e-3 → 3.06e-3 → 4.18e-3 → 5.57e-3). One likely cause not mentioned: in the algorithm frame the slowest mode rate scales as `γ(ΔE_1^alg)` with `σ_alg = 1/β_alg`; at large β_alg, σ shrinks and the smooth-Metro kernel narrows. Worth flagging since this connects to the canonical `s=0.25` decision (.claude-memory: "smooth-Metropolis s convention"). Orthogonal to the bottleneck question; fine to defer. |

No `BLOCKER` or `WARNING`. All physics annotations in the script are accurate, and the comment headers correctly flag deliberate skips (BohrDomain cross-check, trajectory τ_mix, multiseed).

## Suggestions for clarifying the draft

1. The headline "matrix-element bottleneck story is wrong" is correct as stated, but the draft's "What's actually closing?" section (lines 60-67) speculates on Davies/Glauber paths without citing the closest evidence in the qf-biz output itself — namely that `γ` is ~constant in n (0.51→0.56) and `M²` is ~constant in n, so the closing must come from the *eigenvector of L's slowest mode* (not the doublet coherence). Mention this projection mismatch as the cleanest next diagnostic (the draft does suggest it in line 67, just emphasize that this is the only remaining knob).
2. Line 8 phrasing: "the matrix-element bottleneck I proposed earlier... is WRONG." Suggest softening to "...is NOT the mechanism" — the matrix element of |ψ_1⟩↔|ψ_2⟩ really is O(1), but the *L gap* still has a "matrix-element-like" decomposition through its actual slowest left eigenvector. The bottleneck framing isn't wrong as a concept, just wrong about which doublet matters.
3. The "by Z₂ symmetry" σ_x = 0 argument (item 1 in draft Interpretation) is correct but technically the relevant statement is "P = X^⊗n commutes with σ_x^(i), so σ_x is Z₂-EVEN; an operator can connect parities ±1 ↔ ∓1 only if it ANTI-commutes with P." The current text says "P commutes with σ_x^(i), so σ_x is Z₂-even, can't connect" — true but the logical chain skips one step (Z₂-even ⇒ same-sector). 30-second clarification.

## Final notes

- The script obeys all scripts.md rules: standalone runnable, header docstring, `# PHYSICS CHECK` flags, no new deps, BLAS=1 set, JULIA_NUM_THREADS respected.
- `include_coherent` is left at default `true` everywhere (correct per julia-code.md "Coherent Term: ON by Default").
- The `cfg.beta_phys` path is used consistently; β_alg derived via `beta_alg(ham, β_phys)`.
- Sidecar BSON contains all rows + diags + gap eigenvalues; re-runnable for audit.

Verification: **PASS**.

Driver: `/Users/bence/code/QuantumFurnace.jl/scripts/scratch_qf_biz_phase_and_matrix_elements.jl`
Sidecars: `/Users/bence/code/QuantumFurnace.jl/scripts/output/qf_biz_phase_and_matrix_elements/check{1,2}_*.bson`
Run log: `/tmp/qf_biz_run.log`
