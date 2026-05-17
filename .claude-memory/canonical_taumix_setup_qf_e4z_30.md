---
name: canonical-taumix-setup-qf-e4z-30
description: "Canonical 1D τ_mix / spectral-gap setup for the thesis numerics. Z+ZZ-disordered Heisenberg + rho_0 = |+⟩⟨+|^⊗N + single-pass predict_lindbladian_trajectory at krylovdim=60. Eliminates the parity-sector ambiguity — reported gap equals true L gap equals trajectory decay rate, all one number. Replaces the rho_0 = I/d convention. Validated dense at n=5; krylovdim=60 confirmed sufficient at n=7."
metadata:
  type: project
  node_type: memory
---

# Canonical 1D τ_mix / gap setup (qf-e4z.29 + qf-e4z.30)

**Rule:** For all thesis τ_mix and spectral-gap plots on 1D
Heisenberg, use this single configuration:

| Choice | Value |
|---|---|
| Hamiltonian | 1D AFM Heisenberg, **Z+ZZ disorder, strength 0.1** |
| Fixture file | `heis_xxx_disordered_periodic_n{n}_seed{seed}.bson` |
| Initial state | **`rho_0 = \|+⟩⟨+\|^⊗N`** (one Hadamard layer on `\|0⟩^⊗N`) — replace every `Matrix(I(d)/d)` call site |
| Algorithm | `predict_lindbladian_trajectory(cfg, ham, jumps, rho_0, t_grid; krylovdim, tol=1e-10)` — **single-pass effectively** because Pass 2's reported gap matches Pass 1's `eigenvalues[2]` to ~1e-9 once rho_0 is parity-broken (the Pass 2 call inside the predictor becomes a redundant cross-check; remove only if perf matters) |
| Reported gap | `traj.eigenvalues[2]` (= `traj.spectral_gap` here, by construction) |
| krylovdim | **60** for n ∈ {3..9}. Bump to 80 if pushing for machine precision, or at n ≥ 10 if convergence rates degrade (see "n-dependence" below) |
| β_phys grid | Canonical `{0.25, 0.5, 1.0}` per [[β-phys-vs-β-alg-convention-qf-6vr]] |
| n grid | **Odd-n preferred** (`{3, 5, 7, 9}`) to avoid the qf-e4z.23 even/odd-n splitting ambiguity. Keep open for revision once n ≥ 10 data lands and we know whether the splitting is real physics or a parity-sub-spectrum artefact (see [[heisenberg-1d-multiseed-even-odd-qf-e4z-23]]) |

**Why this works.** `\|+⟩⟨+\|^⊗N` is a pure product state with full
off-diagonal coherences in the computational basis (`H^⊗N` applied to
`\|0⟩^⊗N`). As a superoperator-input vector `vec(\|+⟩⟨+\|^⊗N)`, it is
**not** a P̂-eigenstate for P̂[ρ] = Z^⊗N ρ Z^⊗N — in fact
P̂(\|+⟩⟨+\|^⊗N) = \|-⟩⟨-\|^⊗N. So Arnoldi's orbit
`{rho_0, L(rho_0), L²(rho_0), …}` naturally spans both parity sectors of
the Liouvillian, and `_krylov_spectral_decomposition`'s captured
eigenmodes include the true gap mode regardless of `[P, L]`. The
trajectory `e^{tL}(\|+⟩⟨+\|^⊗N)` and the reported gap
`\|eigenvalues[2]\|` then both reflect the true L spectrum, with no
sector-projection caveat.

`I/d` (maximally mixed) is the polar opposite — it is unitary-invariant,
so `P̂(I/d) = I/d` exactly. Arnoldi from `vec(I/d)` is permanently
symmetry-trapped in the parity-even sector on any L with `[P, L] = 0`,
which is why qf-e4z.26/27 had to add the [[predict-trajectory-two-pass-qf-e4z-27]]
patch. `\|+⟩⟨+\|^⊗N` solves the problem at the input, removing the
need for the patch entirely (though we leave the two-pass code as a
safety net for any future caller that does pass `rho_0 = I/d`).

`I/d` is also NOT preparable on hardware with Hadamards: H^⊗N applied
to `\|0⟩^⊗N` gives `\|+⟩⟨+\|^⊗N`, not `I/d`. `I/d` requires noise / a
random coin / partial trace.

**Empirical validation.** Two scripts at the spot cells used in the
qf-e4z.26/27 bug history:

* `scripts/scratch_qf_e4z_29_trajectory_vs_dense.jl` at n=5, β_phys=2.5,
  seed=42. Reference = dense `exp(L_dense · t) · vec(rho_0)` at every
  t in the trajectory grid. With `rho_0 = \|+⟩⟨+\|^⊗N`:
  err_max = **3.9e-14** (Z+ZZ) / 3.7e-14 (X+ZZ) at krylovdim=40 —
  machine precision. `gap_traj = gap_dense` exactly.
  With `rho_0 = I/d`: err_max = 7.1e-9 (the trajectory is fine —
  qf-e4z.27 symmetry argument verified) BUT `gap_traj` is 17% off
  (= leading parity-even mode, not L's true gap).

* `scripts/scratch_qf_e4z_30_kdim_sweep_n7.jl` at n=7, β_phys=2.5,
  seed=42, krylovdim ∈ {40, 60, 80, 100, 120}. With
  `rho_0 = \|+⟩⟨+\|^⊗N`, single-pass `_krylov_spectral_decomposition`
  rel_err vs `krylov_spectral_gap` reference:

  | krylovdim | X+ZZ rel_err | Z+ZZ rel_err |
  |---|---|---|
  | 40  | 2.59e-5  | 9.17e-5  |
  | **60**  | **2.45e-9**  | **2.02e-9**  |
  | 80  | 1.31e-13 | 3.43e-14 |
  | 100 | 4.47e-14 | 1.90e-14 |
  | 120 | 4.63e-14 | 2.41e-14 |

  krylovdim=60 already 4 orders of magnitude under the 1e-7 thesis
  threshold; krylovdim=80 reaches machine precision. **Both fixtures
  behave identically** — the X+ZZ vs Z+ZZ distinction is irrelevant
  once `rho_0` is parity-broken.

**Cost.** Single-pass at krylovdim=60 = 60 matvecs. The legacy
two-pass at krylovdim=40 was ~166 matvecs (Pass 1 = 40, Pass 2 ≈ 126).
~2.75× faster *and* 4 orders of magnitude more accurate. The Pass 2
call inside `predict_lindbladian_trajectory` is now redundant for any
caller using `rho_0 = \|+⟩⟨+\|^⊗N`; the ~50-100 extra Pass 2 matvecs are
a sunk cost not yet optimised out at the API level.

**n-dependence caveat.** krylovdim = 60 is sufficient for n ≤ 9 at our
fixtures (extrapolating from the n=5,7 trend with a small safety margin).
At n ≥ 10 the captured fraction of d² = 4^n drops further and the
Krylov truncation may start mattering: krylovdim=60 captures
60/4^10 ≈ 6e-5 of operator space at n=10. **Empirically verify with a
krylovdim sweep before signing off on thesis plots at n ≥ 10.** The
qf-e4z.30 driver is the template — re-run at the largest n and pick the
smallest krylovdim that reaches rel_err < 1e-7.

**Seed policy (added 2026-05-17 from qf-e4z.34 canonical-seed analysis):**

- n ≤ 8 (cheap sandbox runs): **5 seeds** `{42, 43, 44, 45, 46}` + median + IQR bands — the qf-e4z.34 cost-default.
- **n ≥ 9 (cluster runs): SINGLE SEED = 46.** Validated by qf-e4z.34: across all 72 (sampler × metric × cell) comparisons over n=3..8 × β_phys ∈ {0.25, 0.5, 1.0}, seed=46 is the closest seed to the 5-seed median **in 33/72 cases (46%, far above 1/5 = 20% by chance)** with mean deviation 0.35% from the cell median; the worst-case deviation across ALL cells is 1.99% (τ_CKG) / 1.37% (gap_DLL) / 1.32% (gap_CKG) / 1.26% (τ_DLL). The next-best seed (43) has max τ deviations of 5.4% (CKG) / 9.85% (DLL) — significantly worse.
- This ±2% single-seed proxy is **WELL BELOW** any signal we care about now that the numerics chapter is pivoting to QUALITATIVE understanding rather than asymptotic scaling fits (per [[feedback-more-data-points-for-scaling-claims]] + [[feedback-numerics-grid-canonical]]). Exact numerical scaling exponents would need 5+ seeds; for the "where does CKG help more, what's the parity story, does the gap collapse" qualitative picture, 1 seed is plenty.
- If extra robustness is wanted at a single n=10..14 spot, run 2 seeds `{46, 43}`. Going beyond 2 seeds at cluster scale is wasted compute.

**How to apply:**

- For all NEW τ_mix / gap drivers, write `rho_0 = rho_plus_tensor(n) =
  prod kron of [0.5 0.5; 0.5 0.5]`. Replace every `Matrix(I(d)/d)` in
  existing sweep scripts during their next edit. Use `krylovdim=60`.
- Both `predict_lindbladian_trajectory` and `predict_channel_trajectory`
  benefit equally — `\|+⟩⟨+\|^⊗N` is parity-broken regardless of which
  superoperator (Lindbladian L or channel Φ_δ) it's fed to.
- Keep the existing Z+ZZ canonical fixture family
  [[fixture-migration-find-typical-qf-2kd]]. The X+ZZ family
  [[xz-disorder-fixtures-qf-e4z-27]] becomes redundant — preserved on
  disk for historical audit but not needed for new work.
- Open follow-up: qf-e4z.32 (filed) — re-run qf-e4z.23's 180-cell sweep
  with `\|+⟩⟨+\|^⊗N` to determine whether the even/odd-n splitting is
  real physics or a parity-sub-spectrum artefact from the I/d
  convention. If real physics, the thesis n grid stays odd-only; if
  artefact, even-n becomes available for the larger system-size data
  points and the n grid widens to `{3, 5, 7, 9, ...}` or `{3, 4, 5, ..., 9}`.

**Drivers (validation):**
- `scripts/scratch_qf_e4z_29_trajectory_vs_dense.jl` — n=5 dense reference.
- `scripts/scratch_qf_e4z_30_kdim_sweep_n7.jl` — n=7 krylovdim convergence.

Output BSON: `scripts/output/qf_e4z_29/`.
