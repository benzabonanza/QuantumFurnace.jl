---
name: predict-trajectory-parity-bug-qf-e4z-26
description: "predict_lindbladian_trajectory had the qf-8fr symmetric-x0 bug too — qf-e4z.26 patched with a 1e-6 GUE perturbation (fails at n=7). Superseded by qf-e4z.27 two-pass architecture — see [[predict-trajectory-two-pass-qf-e4z-27]]."
metadata:
  type: feedback
  node_type: memory
---

# predict_lindbladian_trajectory parity-trapped Arnoldi — qf-e4z.26

> **SUPERSEDED 2026-05-16.** The qf-e4z.26 1e-6-GUE-perturbed seed in
> `_krylov_spectral_decomposition` caught the parity-projection bug at
> n ≤ 5 but FAILED at n = 7 (gap 0.1203 reported vs true 0.0966 from
> `krylov_spectral_gap`). The proper architectural fix is
> [[predict-trajectory-two-pass-qf-e4z-27]] (two-pass Krylov: single-
> seed Arnoldi for trajectory + separate `krylov_spectral_gap` call
> for the gap). Below is the historical write-up — kept for context.

**Bug.** `predict_lindbladian_trajectory` and `predict_channel_trajectory`
(`src/lindblad_action.jl`) both went through `_krylov_spectral_decomposition`,
which seeded Arnoldi with `x_0 = vec(rho_0)` (default `rho_0 = I/d`).
When the Lindbladian commutes with a parity superoperator `P̂[ρ] = PρP`
(weak symmetry, in the Buča–Prosen sense), Arnoldi from a parity-even
seed stays trapped in the parity-even sector and silently misses every
parity-odd eigenmode. The reported `spectral_gap` is then the slowest
parity-even rate, not the true Lindbladian gap.

**Why:** The qf-8fr fix (2026-05-13, see
[[krylov-x0-symmetric-bug-qf-8fr]]) addressed this trap in
`krylov_spectral_gap` but only patched `src/krylov_eigsolve.jl`. The
sister code path inside `_krylov_spectral_decomposition` was left with
the unpatched seed, with a docstring that explicitly rationalised "pure
I/d would also work for the forward problem" — wrong whenever
`rho_0` is symmetry-protected. The qf-8fr write-up further claimed
disordered Heisenberg fixtures "break the symmetry"; that claim is
false for the production fixture family
(`build_heis_1d`, base XX+YY+ZZ, disorder [[Z], [Z,Z]]) because every
listed term commutes with `P = Z^⊗N`.

**Impact on qf-e4z.23:** The 180-cell multi-seed sweep ran
`predict_lindbladian_trajectory(rho_0 = I/d, ...)` and read
`gap_arnoldi = traj.spectral_gap`. Dense LAPACK at the suspect cells
shows:
- `n=5, β_phys=2.5, seed=42`: sidecar `0.1607`, dense `0.1378` (17%
  gap overestimate; missed parity-odd doublet at -0.1378).
- `n=7, β_phys=2.5, seed=42`: sidecar `0.1203`, krylov_spectral_gap
  `0.0966` (25% overestimate; missed parity-odd doublet at -0.0966).
- n=3 seed=46 cells broke Arnoldi at exactly 32 = d²/2 matvecs —
  numerical confirmation that Krylov stayed inside the parity-even
  sector of dimension 32.

**τ_mix is unchanged.** Because `rho_0 = I/d` and `σ_β` are both
parity-even, the residual `rho_0 - σ_β` lives entirely in the
parity-even sector. The biorthogonal projection
`c_i = ⟨L_i, rho_0 - σ_β⟩_HS` is identically zero on parity-odd `L_i`,
so parity-odd eigenmodes contribute nothing to the trajectory. The
`mixing_time` field in the sidecars is physically correct — it is the
mixing time *from `I/d`*, which is governed by the slowest *parity-even*
rate, which is what the buggy Arnoldi reports.

**How to apply (POST qf-e4z.27):**
- `_krylov_spectral_decomposition` is single-seed Arnoldi from
  `vec(rho_0)` again. `eigenvalues[2]` is the P̂-EVEN-sub-spectrum
  gap (CORRECT decay rate for the trajectory of `rho_0`).
- `traj.spectral_gap` is sourced from a separate `krylov_spectral_gap`
  pass and is the TRUE Lindbladian gap at every n.
- See [[predict-trajectory-two-pass-qf-e4z-27]] for the full
  architectural rationale and verification.
- Sidecars from qf-e4z.{20-25} are STILL pre-fix on `gap_arnoldi`;
  τ_mix data stands.
- When designing a Lindbladian fixture / cross-check / sweep, FIRST
  determine whether the Lindbladian inherits a weak symmetry from the
  Hamiltonian + jump set, by checking commutation of each term in
  `_build_jump_set` with the candidate symmetry. Z^⊗N is conserved by
  any pure-Z disorder; X^⊗N is conserved by any pure-X disorder.
  Better still: use a parity-breaking disorder family like
  [[X],[Z,Z]] — see [[xz-disorder-fixtures-qf-e4z-27]].

**Cells written:**
- `drafts/qf-e4z-23-audit-findings.md` — full audit write-up
- `scripts/scratch_audit_e4z23_spectrum.jl` — dense vs Krylov gap
- `scripts/scratch_audit_e4z23_parity_gap.jl` — krylov_spectral_gap vs
  sidecar comparison + symmetry verification
- `scripts/scratch_audit_e4z23_parity_decomp.jl` — parity-signature
  decomposition of dense eigenmodes
- `scripts/scratch_audit_e4z23_verify_fix.jl` — post-fix verification

**Open**: regression test for the patched
`predict_lindbladian_trajectory` path should be added (mirror the
existing qf-8fr regression in `test/test_krylov_eigsolve.jl` on
classical Ising n=3). Re-run qf-e4z.23 to refresh `gap_arnoldi` values
in sidecars.
