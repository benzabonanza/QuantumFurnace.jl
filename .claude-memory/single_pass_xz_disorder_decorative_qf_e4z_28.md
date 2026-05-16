---
name: single-pass-xz-disorder-decorative-qf-e4z-28
description: "qf-e4z.28 finding: X+ZZ structural parity-break does NOT eliminate the qf-e4z.27 parity-trap bug. The 0.1 X-disorder is too weak to survive 40-step MGS Arnoldi from vec(I/d) — single-pass gives the SAME wrong gap as Z+ZZ. Even a 1e-10 GUE-perturbed seed fails in plain Arnoldi; KrylovKit thick-restart in Pass 2 is the actual fix. Keep the qf-e4z.27 two-pass."
metadata: 
  node_type: memory
  type: feedback
  originSessionId: eb034d2e-3830-4584-91f6-0d52f53ce2fe
---

# qf-e4z.28: X+ZZ structural fix is decorative — keep two-pass

**Rule:** Do NOT drop the qf-e4z.27 two-pass in `predict_lindbladian_trajectory`
/ `predict_channel_trajectory` on the strength of "X+ZZ breaks Z^⊗N parity
structurally". The two-pass with KrylovKit thick-restart in Pass 2 is
the algorithmic fix; X+ZZ is a physics-level convenience, not an
algorithmic substitute.

**Why:** Empirical test on the historically problematic spot cells
(n=5 and n=7, β_phys=2.5, seed=42; qf-e4z.27 verification drivers).
At krylovdim=40, `_krylov_spectral_decomposition` seeded with vec(I/d):

| Cell | gap_1pass (X+ZZ) | gap_1pass (Z+ZZ) | true gap |
|---|---|---|---|
| n=5  | 0.1608 | 0.1607 | 0.1378 |
| n=7  | 0.1201 | 0.1203 | 0.0966 |

Identical 16-24% errors on both fixtures despite ‖[P, H]‖ = 0.30 (n=5)
/ 0.45 (n=7) on X+ZZ — the X-disorder *does* break Z^⊗N parity at the
Hamiltonian level (verified by direct commutator). But the disorder
strength (0.1) is small relative to the dominant ZZ exchange, so the
parity-odd content the X-term injects into the Krylov orbit of vec(I/d)
is below the MGS orthogonalisation tol after a handful of steps. The
captured 40-d Krylov subspace stays inside the *approximately*
parity-even sector and reports the leading parity-even eigenvalue,
exactly as on Z+ZZ.

Even substituting Pass 2's seed `vec(I/d + 1e-10·H_GUE_traceless)` into
plain Arnoldi (Test E) fails identically on both fixtures — the 1e-10
GUE-perturbed seed only works with KrylovKit's thick-restart Krylov-Schur,
which re-amplifies the rare parity-odd content via implicit restart at
each cycle. Plain MGS Arnoldi has no restart; once a small component
falls below ~`tol` it is gone for good.

A `|+⟩⟨+|^⊗N` initial state DOES break the trap at n=5 (rel_err 2.5e-7)
but only partially at n=7 (rel_err 2.6e-5 X+ZZ / 9.2e-5 Z+ZZ — above the
1e-6 threshold for thesis-grade work). Not a robust replacement.

**How to apply:**
- New `predict_lindbladian_trajectory` / `predict_channel_trajectory`
  callers: trust `traj.spectral_gap` (= Pass 2). The matvec overhead
  (~50–100 vs single-pass) is required.
- X+ZZ-disordered fixtures (`heis_xxx_XZdisordered_*`) remain the
  *physics-correct* family per [[xz-disorder-fixtures-qf-e4z-27]]: they
  remove the parity artefact at the Hamiltonian level and avoid the
  qf-e4z.23 even/odd-n splitting interpretation. But they do not buy
  any matvec saving in the trajectory predictor — Pass 2 is still
  required for the true gap.
- Do not propose "drop Pass 2 when ‖[P, H]‖ > 0" optimisations. The
  threshold for "parity-broken enough to survive plain Arnoldi" would
  scale with the disorder strength times the number of Krylov steps —
  not a robust API contract.
- If a future caller really needs to skip Pass 2 (e.g. perf-critical
  inner-loop usage where the rho_0 trajectory is the only output and
  the user doesn't care about `spectral_gap`), expose a kwarg
  `compute_true_gap::Bool = true` rather than auto-detecting. The
  trajectory itself (eigenvalues, c, R_modes from Pass 1) is correct
  for any P̂-even rho_0 even on parity-symmetric L — only the reported
  `spectral_gap` is affected. Currently no such caller exists.

**Driver:** `scripts/scratch_qf_e4z_28_single_pass.jl` (5 tests:
A–E spanning rho_0 ∈ {I/d, |+⟩^⊗N, I/d + 1e-10·H_GUE} × fixture ∈
{X+ZZ, Z+ZZ} at n ∈ {5, 7}, β_phys = 2.5, seed = 42). Full output
preserved at `scripts/output/qf_e4z_28_single_pass_log.txt`.

**Files touched:** none in `src/` — investigation only.
