---
name: predict-trajectory-two-pass-qf-e4z-27
description: "qf-e4z.27 architectural fix for the predict_lindbladian_trajectory / predict_channel_trajectory parity-projection bug — two-pass Krylov (single-seed Arnoldi for trajectory + krylov_spectral_gap for the true gap). Supersedes the qf-e4z.26 transitional eps=1e-6 band-aid."
metadata:
  type: feedback
  node_type: memory
---

# qf-e4z.27 fix: two-pass Krylov for trajectory predictor

**Supersedes** the qf-e4z.26 band-aid (`vec(rho_0 + 1e-6 H_GUE)` seed in
`_krylov_spectral_decomposition`), which caught the parity-projection
bug at n ≤ 5 but failed at n = 7 (gap 0.1203 reported vs true 0.0966
from `krylov_spectral_gap`). See [[predict-trajectory-parity-bug-qf-e4z-26]]
for the bug history.

**Architecture.** `predict_lindbladian_trajectory` and
`predict_channel_trajectory` (`src/lindblad_action.jl`) now use TWO
complementary Krylov passes:

1. **Trajectory pass** — `_krylov_spectral_decomposition` runs single-seed
   `_arnoldi_factorize` from `x_0 = vec(rho_0)`. The captured Krylov
   subspace is tightly aligned with rho_0's orbit ⇒ rho_0 ∈ span(Q)
   exactly, trajectory reconstructs to ~1e-7 at every t. The eigenvalues
   returned reflect the symmetry sector rho_0 lives in (P̂-EVEN-only for
   rho_0 = I/d on a parity-symmetric L); this is the CORRECT decay rate
   for the trajectory of that rho_0, because P̂-odd modes have c_i ≡ 0
   by symmetry.
2. **Gap pass** — `krylov_spectral_gap` on the same `(config, ham, jumps)`
   with `_krylov_default_x0(dim) = vec(I/d + 1e-10 H_GUE)` + KrylovKit
   thick-restart Krylov–Schur. Captures parity-odd modes reliably at
   every n (restart re-amplifies the 1e-10 perturbation, vs unrestarted
   MGS Arnoldi which loses content < O(eps^k) over k steps). Reports
   the TRUE Lindbladian gap.

The returned NamedTuple's `spectral_gap` field comes from pass 2.
Trajectory, `eigenvalues`, `R_modes`, `c`, `rho_inf` come from pass 1.
`total_matvecs = decomp.matvec_count + gap_res.matvec_count` reports the
cumulative cost (typically +50–100 matvecs vs the pre-qf-e4z.27 single
pass).

**Why not a single dual-seed Arnoldi?** Tested in qf-e4z.27 dev: a
dual-seed Arnoldi with `(vec(rho_0), vec(H_GUE))` gives the true gap
AND keeps rho_0 in span(Q), BUT the captured subspace's *budget* is
split between rho_0's orbit and GUE's orbit. For rho_0 = I/d on
parity-symmetric L, half the Krylov dims are wasted on parity-odd
modes that rho_0 cannot populate, degrading parity-even trajectory
accuracy from ~1e-7 to ~1e-3. The two-pass approach pays a single
extra `krylov_spectral_gap` call (~100 matvecs) to recover both.

**Why not KrylovKit eigsolve for both?** Tested in qf-e4z.27 dev:
KrylovKit's eigsolve with Gram-inverse biorthogonality breaks for the
near-degenerate slow modes characteristic of KMS-DB Lindbladians.
Captured right eigvecs of degenerate eigenvalues are nearly parallel
⇒ `G = R'R` is nearly singular ⇒ `L = R G⁻¹` is unreliable ⇒
biorthogonality violated ⇒ trajectory at t=0 has O(1) error. The
plain MGS Arnoldi + Q-lift in `_krylov_spectral_decomposition`
preserves biorthogonality exactly via the orthonormal Q, regardless
of degeneracy.

**Why:** The qf-e4z.26 single-pass band-aid fundamentally cannot work
at large n because plain MGS Arnoldi has no restart — a small 1e-6
seed component cannot survive 40 GS sweeps when its overlap with the
target eigenmode is exponentially small in the spectral separation
times the number of steps.

**How to apply:**
- New trajectory work: trust `traj.spectral_gap` post-qf-e4z.27. It is
  the true Lindbladian gap (matches `krylov_spectral_gap` to KrylovKit
  tol = 1e-10).
- Sidecars from qf-e4z.{20-25} (pre-qf-e4z.27) have wrong
  `gap_arnoldi` values; τ_mix data stands. Re-run the sweep if you
  need accurate gap_arnoldi values.
- For new fixture families: prefer X-type (or any non-Z) disorder over
  Z-type for AFM Heisenberg. [[X],[Z,Z]] structurally breaks Z^⊗N parity
  at the Hamiltonian level. NOTE per
  [[single-pass-xz-disorder-decorative-qf-e4z-28]]: the 0.1 X-disorder
  in our X+ZZ fixtures is too weak to remove the algorithmic parity
  trap from plain Arnoldi at krylovdim=40 — Pass 2 still required on
  those fixtures. The X+ZZ family is the right physics fixture (no
  even/odd-n parity artefact in qf-e4z.23-style sweeps), not an
  algorithmic substitute. See [[xz-disorder-fixtures-qf-e4z-27]].
- Cost: the qf-e4z.27 two-pass approach adds ~50–100 matvecs vs the
  pre-qf-e4z.27 single pass. For n ≥ 7 sweeps this is meaningful;
  for thesis-grade work this is the right default.

**Regression tests:**
- `test/test_predict_lindbladian.jl::(f)` — classical 1D Ising n=3
  rtol tightened from 1e-6 → 1e-8 (krylov_spectral_gap path converges
  to its `tol=1e-10` floor by construction, no MGS noise to fight).
- `test/test_predict_channel.jl::(d)` (NEW) — sister test for
  `predict_channel_trajectory` on the same classical Ising fixture.
  Tolerance is 1e-3 because channel μ → λ = (μ-1)/δ conversion has
  O(δ) error from the channel-vs-generator first-order approximation.

**Files touched:**
- `src/lindblad_action.jl`
  - `_krylov_spectral_decomposition` docstring updated to flag the
    parity caveat and point at `predict_*_trajectory` for the true gap.
  - `predict_lindbladian_trajectory` + `predict_channel_trajectory`
    add a `krylov_spectral_gap` call after the decomp; replace
    `spectral_gap` field with the true Lindbladian gap.
- `test/test_predict_lindbladian.jl::(f)` — rtol 1e-6 → 1e-8.
- `test/test_predict_channel.jl::(d)` — new regression.
- `scripts/scratch_qf_e4z_27_verification.jl` — Part D spot-cell
  verification driver.
- `scripts/scratch_multiseed_disordered_fixtures.jl` — adds
  `run_1d_heisenberg_XZ_disorder` for the new [[X],[Z,Z]] family.

**Spot-cell verification.** qf-e4z.27 spot cell n=7, β_phys=2.5,
seed=42 on the legacy Z+ZZ fixture: `predict_lindbladian_trajectory`
gap = `krylov_spectral_gap` reference to relative error 0.0 (= machine
precision). Pre-qf-e4z.27 reported 0.1203 (24 % off from true 0.0966).
Driver: `scripts/scratch_qf_e4z_27_verification.jl`.
