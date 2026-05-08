# Dissipative L_diss quadrature error — TimeDomain vs EnergyDomain (qf-7xt)

## Method

Compare TimeDomain (NUFFT-based discrete ω-sum on `[-N/2, N/2-1]·w0_D`,
`N = 2^r_D`) against EnergyDomain reference at r_ref=16 with
w0_ref = (2·‖H‖+5σ)/2^r_ref ≪ all test w0_D. EnergyDomain at this
resolution agrees with BohrDomain to ~2e-14 (Gaussian) and ~4e-6
(smooth/kinky Metro) — see `scripts/scratch_energy_ref_convergence.jl`.

Operator-norm error is computed matrix-free via
`hs_operator_norm_krylov` (KrylovKit.svdsolve, GKL bidiagonalization)
on the difference closure, avoiding the OOM that the BohrDomain dense
build hits at n=6 and would hit catastrophically at n=7-8. Cross-checked
vs the dense path at n=4 (testset 8 in `test_kms_geometry.jl`) to
machine precision.  η = 0.001 for Metropolis-like b_+.

## Section 1: r_D sweep at w0_D = π/(5β) (β-adaptive)

Predicted slopes (thesis §quad-errors-diss):
- Gaussian: super-algebraic / exponential in r_D.
- Smooth Metro (s=0.25): super-algebraic, slower than Gaussian.
- Kinky Metro (s=0): algebraic slope -2 in `2^r_D` (∝ ω0_D²) due to ω-kink.

| n | β | filter | w0_D | r_D | ‖ΔL‖_op | ‖L_ref‖ |
|---|---|---|---|---|---|---|
| 3 | 5.0 | gaussian | 0.1257 | 6 | 5.140e-12 | 7.759e-01 |
| 3 | 5.0 | gaussian | 0.1257 | 7 | 4.588e-12 | 7.759e-01 |
| 3 | 5.0 | gaussian | 0.1257 | 8 | 5.236e-12 | 7.759e-01 |
| 3 | 5.0 | gaussian | 0.1257 | 9 | 4.858e-12 | 7.759e-01 |
| 3 | 5.0 | gaussian | 0.1257 | 10 | 4.817e-12 | 7.759e-01 |
| 3 | 5.0 | gaussian | 0.1257 | 11 | 4.261e-12 | 7.759e-01 |
| 3 | 5.0 | gaussian | 0.1257 | 12 | 4.656e-12 | 7.759e-01 |
| 3 | 5.0 | gaussian | 0.1257 | 13 | 4.223e-12 | 7.759e-01 |
| 3 | 5.0 | kinky | 0.1257 | 6 | 1.080e-03 | 9.636e-01 |
| 3 | 5.0 | kinky | 0.1257 | 7 | 1.080e-03 | 9.636e-01 |
| 3 | 5.0 | kinky | 0.1257 | 8 | 1.080e-03 | 9.636e-01 |
| 3 | 5.0 | kinky | 0.1257 | 9 | 1.080e-03 | 9.636e-01 |
| 3 | 5.0 | kinky | 0.1257 | 10 | 1.080e-03 | 9.636e-01 |
| 3 | 5.0 | kinky | 0.1257 | 11 | 1.080e-03 | 9.636e-01 |
| 3 | 5.0 | kinky | 0.1257 | 12 | 1.080e-03 | 9.636e-01 |
| 3 | 5.0 | kinky | 0.1257 | 13 | 1.080e-03 | 9.636e-01 |
| 3 | 5.0 | smooth | 0.1257 | 6 | 2.719e-07 | 9.145e-01 |
| 3 | 5.0 | smooth | 0.1257 | 7 | 2.719e-07 | 9.145e-01 |
| 3 | 5.0 | smooth | 0.1257 | 8 | 2.719e-07 | 9.145e-01 |
| 3 | 5.0 | smooth | 0.1257 | 9 | 2.719e-07 | 9.145e-01 |
| 3 | 5.0 | smooth | 0.1257 | 10 | 2.719e-07 | 9.145e-01 |
| 3 | 5.0 | smooth | 0.1257 | 11 | 2.719e-07 | 9.145e-01 |
| 3 | 5.0 | smooth | 0.1257 | 12 | 2.719e-07 | 9.145e-01 |
| 3 | 5.0 | smooth | 0.1257 | 13 | 2.719e-07 | 9.145e-01 |
| 3 | 10.0 | gaussian | 0.0628 | 6 | 2.948e-12 | 6.145e-01 |
| 3 | 10.0 | gaussian | 0.0628 | 7 | 2.228e-12 | 6.145e-01 |
| 3 | 10.0 | gaussian | 0.0628 | 8 | 2.746e-12 | 6.145e-01 |
| 3 | 10.0 | gaussian | 0.0628 | 9 | 2.767e-12 | 6.145e-01 |
| 3 | 10.0 | gaussian | 0.0628 | 10 | 2.282e-12 | 6.145e-01 |
| 3 | 10.0 | gaussian | 0.0628 | 11 | 2.400e-12 | 6.145e-01 |
| 3 | 10.0 | gaussian | 0.0628 | 12 | 2.396e-12 | 6.145e-01 |
| 3 | 10.0 | gaussian | 0.0628 | 13 | 2.801e-12 | 6.145e-01 |
| 3 | 10.0 | kinky | 0.0628 | 6 | 8.412e-04 | 9.869e-01 |
| 3 | 10.0 | kinky | 0.0628 | 7 | 8.412e-04 | 9.869e-01 |
| 3 | 10.0 | kinky | 0.0628 | 8 | 8.412e-04 | 9.869e-01 |
| 3 | 10.0 | kinky | 0.0628 | 9 | 8.412e-04 | 9.869e-01 |
| 3 | 10.0 | kinky | 0.0628 | 10 | 8.412e-04 | 9.869e-01 |
| 3 | 10.0 | kinky | 0.0628 | 11 | 8.412e-04 | 9.869e-01 |
| 3 | 10.0 | kinky | 0.0628 | 12 | 8.412e-04 | 9.869e-01 |
| 3 | 10.0 | kinky | 0.0628 | 13 | 8.412e-04 | 9.869e-01 |
| 3 | 10.0 | smooth | 0.0628 | 6 | 2.487e-07 | 9.504e-01 |
| 3 | 10.0 | smooth | 0.0628 | 7 | 2.487e-07 | 9.504e-01 |
| 3 | 10.0 | smooth | 0.0628 | 8 | 2.487e-07 | 9.504e-01 |
| 3 | 10.0 | smooth | 0.0628 | 9 | 2.487e-07 | 9.504e-01 |
| 3 | 10.0 | smooth | 0.0628 | 10 | 2.487e-07 | 9.504e-01 |
| 3 | 10.0 | smooth | 0.0628 | 11 | 2.487e-07 | 9.504e-01 |
| 3 | 10.0 | smooth | 0.0628 | 12 | 2.487e-07 | 9.504e-01 |
| 3 | 10.0 | smooth | 0.0628 | 13 | 2.487e-07 | 9.504e-01 |
| 3 | 20.0 | gaussian | 0.0314 | 6 | 2.664e-12 | 4.698e-01 |
| 3 | 20.0 | gaussian | 0.0314 | 7 | 2.891e-12 | 4.698e-01 |
| 3 | 20.0 | gaussian | 0.0314 | 8 | 2.899e-12 | 4.698e-01 |
| 3 | 20.0 | gaussian | 0.0314 | 9 | 2.576e-12 | 4.698e-01 |
| 3 | 20.0 | gaussian | 0.0314 | 10 | 2.814e-12 | 4.698e-01 |
| 3 | 20.0 | gaussian | 0.0314 | 11 | 2.736e-12 | 4.698e-01 |
| 3 | 20.0 | gaussian | 0.0314 | 12 | 3.195e-12 | 4.698e-01 |
| 3 | 20.0 | gaussian | 0.0314 | 13 | 2.872e-12 | 4.698e-01 |
| 3 | 20.0 | kinky | 0.0314 | 6 | 7.730e-04 | 9.838e-01 |
| 3 | 20.0 | kinky | 0.0314 | 7 | 7.730e-04 | 9.838e-01 |
| 3 | 20.0 | kinky | 0.0314 | 8 | 7.730e-04 | 9.838e-01 |
| 3 | 20.0 | kinky | 0.0314 | 9 | 7.730e-04 | 9.838e-01 |
| 3 | 20.0 | kinky | 0.0314 | 10 | 7.730e-04 | 9.838e-01 |
| 3 | 20.0 | kinky | 0.0314 | 11 | 7.730e-04 | 9.838e-01 |
| 3 | 20.0 | kinky | 0.0314 | 12 | 7.730e-04 | 9.838e-01 |
| 3 | 20.0 | kinky | 0.0314 | 13 | 7.730e-04 | 9.838e-01 |
| 3 | 20.0 | smooth | 0.0314 | 6 | 1.518e-07 | 9.516e-01 |
| 3 | 20.0 | smooth | 0.0314 | 7 | 1.518e-07 | 9.516e-01 |
| 3 | 20.0 | smooth | 0.0314 | 8 | 1.518e-07 | 9.516e-01 |
| 3 | 20.0 | smooth | 0.0314 | 9 | 1.518e-07 | 9.516e-01 |
| 3 | 20.0 | smooth | 0.0314 | 10 | 1.518e-07 | 9.516e-01 |
| 3 | 20.0 | smooth | 0.0314 | 11 | 1.518e-07 | 9.516e-01 |
| 3 | 20.0 | smooth | 0.0314 | 12 | 1.518e-07 | 9.516e-01 |
| 3 | 20.0 | smooth | 0.0314 | 13 | 1.518e-07 | 9.516e-01 |

## Section 2: w0_D sweep (r_D = 10, w0 = mult × π/(5β))

| n | β | filter | mult | w0_D | ‖ΔL‖_op | ‖L_ref‖ |
|---|---|---|---|---|---|---|
| 3 | 5.0 | gaussian | 0.50 | 0.0628 | 1.420e-13 | 7.759e-01 |
| 3 | 5.0 | gaussian | 0.75 | 0.0942 | 1.234e-13 | 7.759e-01 |
| 3 | 5.0 | gaussian | 1.00 | 0.1257 | 5.295e-12 | 7.759e-01 |
| 3 | 5.0 | gaussian | 1.50 | 0.1885 | 1.285e-05 | 7.759e-01 |
| 3 | 5.0 | gaussian | 2.00 | 0.2513 | 1.355e-03 | 7.759e-01 |
| 3 | 5.0 | gaussian | 4.00 | 0.5027 | 1.574e-01 | 7.759e-01 |
| 3 | 5.0 | kinky | 0.50 | 0.0628 | 1.519e-03 | 9.636e-01 |
| 3 | 5.0 | kinky | 0.75 | 0.0942 | 4.980e-03 | 9.636e-01 |
| 3 | 5.0 | kinky | 1.00 | 0.1257 | 1.080e-03 | 9.636e-01 |
| 3 | 5.0 | kinky | 1.50 | 0.1885 | 1.549e-02 | 9.636e-01 |
| 3 | 5.0 | kinky | 2.00 | 0.2513 | 2.818e-02 | 9.636e-01 |
| 3 | 5.0 | kinky | 4.00 | 0.5027 | 8.411e-02 | 9.636e-01 |
| 3 | 5.0 | smooth | 0.50 | 0.0628 | 2.199e-13 | 9.145e-01 |
| 3 | 5.0 | smooth | 0.75 | 0.0942 | 3.550e-11 | 9.145e-01 |
| 3 | 5.0 | smooth | 1.00 | 0.1257 | 2.719e-07 | 9.145e-01 |
| 3 | 5.0 | smooth | 1.50 | 0.1885 | 1.938e-04 | 9.145e-01 |
| 3 | 5.0 | smooth | 2.00 | 0.2513 | 2.815e-03 | 9.145e-01 |
| 3 | 5.0 | smooth | 4.00 | 0.5027 | 7.130e-02 | 9.145e-01 |
| 3 | 10.0 | gaussian | 0.50 | 0.0314 | 1.202e-13 | 6.145e-01 |
| 3 | 10.0 | gaussian | 0.75 | 0.0471 | 1.353e-13 | 6.145e-01 |
| 3 | 10.0 | gaussian | 1.00 | 0.0628 | 2.214e-12 | 6.145e-01 |
| 3 | 10.0 | gaussian | 1.50 | 0.0942 | 1.134e-05 | 6.145e-01 |
| 3 | 10.0 | gaussian | 2.00 | 0.1257 | 1.197e-03 | 6.145e-01 |
| 3 | 10.0 | gaussian | 4.00 | 0.2513 | 1.161e-01 | 6.145e-01 |
| 3 | 10.0 | kinky | 0.50 | 0.0314 | 1.174e-03 | 9.869e-01 |
| 3 | 10.0 | kinky | 0.75 | 0.0471 | 3.862e-03 | 9.869e-01 |
| 3 | 10.0 | kinky | 1.00 | 0.0628 | 8.412e-04 | 9.869e-01 |
| 3 | 10.0 | kinky | 1.50 | 0.0942 | 1.212e-02 | 9.869e-01 |
| 3 | 10.0 | kinky | 2.00 | 0.1257 | 2.090e-02 | 9.869e-01 |
| 3 | 10.0 | kinky | 4.00 | 0.2513 | 6.777e-02 | 9.869e-01 |
| 3 | 10.0 | smooth | 0.50 | 0.0314 | 2.090e-13 | 9.504e-01 |
| 3 | 10.0 | smooth | 0.75 | 0.0471 | 2.359e-11 | 9.504e-01 |
| 3 | 10.0 | smooth | 1.00 | 0.0628 | 2.487e-07 | 9.504e-01 |
| 3 | 10.0 | smooth | 1.50 | 0.0942 | 1.838e-04 | 9.504e-01 |
| 3 | 10.0 | smooth | 2.00 | 0.1257 | 2.323e-03 | 9.504e-01 |
| 3 | 10.0 | smooth | 4.00 | 0.2513 | 5.682e-02 | 9.504e-01 |
| 3 | 20.0 | gaussian | 0.50 | 0.0157 | 5.381e-14 | 4.698e-01 |
| 3 | 20.0 | gaussian | 0.75 | 0.0236 | 5.107e-14 | 4.698e-01 |
| 3 | 20.0 | gaussian | 1.00 | 0.0314 | 3.197e-12 | 4.698e-01 |
| 3 | 20.0 | gaussian | 1.50 | 0.0471 | 1.042e-05 | 4.698e-01 |
| 3 | 20.0 | gaussian | 2.00 | 0.0628 | 1.297e-03 | 4.698e-01 |
| 3 | 20.0 | gaussian | 4.00 | 0.1257 | 1.159e-01 | 4.698e-01 |
| 3 | 20.0 | kinky | 0.50 | 0.0157 | 1.008e-03 | 9.838e-01 |
| 3 | 20.0 | kinky | 0.75 | 0.0236 | 3.322e-03 | 9.838e-01 |
| 3 | 20.0 | kinky | 1.00 | 0.0314 | 7.730e-04 | 9.838e-01 |
| 3 | 20.0 | kinky | 1.50 | 0.0471 | 1.049e-02 | 9.838e-01 |
| 3 | 20.0 | kinky | 2.00 | 0.0628 | 1.795e-02 | 9.838e-01 |
| 3 | 20.0 | kinky | 4.00 | 0.1257 | 6.053e-02 | 9.838e-01 |
| 3 | 20.0 | smooth | 0.50 | 0.0157 | 1.263e-13 | 9.516e-01 |
| 3 | 20.0 | smooth | 0.75 | 0.0236 | 1.934e-11 | 9.516e-01 |
| 3 | 20.0 | smooth | 1.00 | 0.0314 | 1.518e-07 | 9.516e-01 |
| 3 | 20.0 | smooth | 1.50 | 0.0471 | 1.555e-04 | 9.516e-01 |
| 3 | 20.0 | smooth | 2.00 | 0.0628 | 1.996e-03 | 9.516e-01 |
| 3 | 20.0 | smooth | 4.00 | 0.1257 | 5.069e-02 | 9.516e-01 |

## Recommended r_D per ε_target at w0_D = π/(5β)

| n | β | filter | ε=1e-3 | ε=1e-4 | ε=1e-5 | ε=1e-6 |
|---|---|---|---|---|---|---|
| 3 | 5.0 | gaussian | 6 | 6 | 6 | 6 |
| 3 | 5.0 | kinky | — | — | — | — |
| 3 | 5.0 | smooth | 6 | 6 | 6 | 6 |
| 3 | 10.0 | gaussian | 6 | 6 | 6 | 6 |
| 3 | 10.0 | kinky | 6 | — | — | — |
| 3 | 10.0 | smooth | 6 | 6 | 6 | 6 |
| 3 | 20.0 | gaussian | 6 | 6 | 6 | 6 |
| 3 | 20.0 | kinky | 6 | — | — | — |
| 3 | 20.0 | smooth | 6 | 6 | 6 | 6 |

---

## Section 3: Dissipative-only Krylov-SVD sweep — unified ω-range (canonical)

> **Supersedes sections 1 and 2**, which used a confounded fixed-w0 methodology.

Method: `include_coherent=false` on both `apply_lindbladian!` matvecs, so only the dissipator
$\mathcal{L}_{\text{diss}}(\rho) = \sum_a L_a \rho L_a^\dagger - \tfrac{1}{2}\{R, \rho\}$
is compared. Norm computed matrix-free via `hs_operator_norm_krylov` (GKL bidiagonalization).

Reference: EnergyDomain at $r_{\text{ref}} = 8$ (smooth Metro saturates by $r = 6$).
Test: TimeDomain at varying $r_D$ (NUFFT-based discrete $\omega$-sum).
Both at the same unified $\omega_{\text{range}} = 2(\|H\| + 8\sigma)$, $\beta = 10$, $\sigma = 0.1$, TAIL_C = 8.

Source: `scripts/scratch_dissipative_krylov_quadrature.jl`.

### Smooth Metropolis (s=0.25), β=10, σ=0.1

Cross-referenced against canonical $n=4$ dense data from `quadrature-convergence-summary.md`
(which used TimeDomain vs BohrDomain, dense both sides).

| $r_D$ | n=4 (dense, Time→Bohr) | n=5 (Krylov, Time→Energy) | n=6 (Krylov) | n=7 (Krylov) |
|---|---|---|---|---|
| 3 | 1.2e-1 | 1.20e-1 | 1.25e-1 | 1.18e-1 |
| 4 | 9.9e-3 | 9.73e-3 | 1.15e-2 | 1.23e-2 |
| 5 | 1.5e-5 | 1.27e-5 | 1.48e-5 | 1.32e-5 |
| 6 | 2.5e-13 | 6.80e-14 | 5.09e-14 | 4.88e-14 |

**Finding**: super-algebraic convergence confirmed n-independent across n=4,5,6,7.
At $r_D = 5$: $\|\Delta L_{\text{diss}}\| \approx 1.3{-}1.5 \times 10^{-5}$ uniformly (spread < 1.2×).
At $r_D = 6$: saturates at the NUFFT precision floor ($\sim 5{-}7 \times 10^{-14}$ via Krylov;
$\sim 2.5 \times 10^{-13}$ via dense at n=4). Floor difference is reference-method dependent
(EnergyDomain $r=8$ vs BohrDomain), not physically meaningful — both are $\ll 10^{-6}$.

Wall times (8 Julia threads, 4 BLAS): n=5: 1–15s per point; n=6: 14–37s;
n=7: 233–502s (reference matvec at $r_{\text{ref}}=8$ with $d=128$ dominates).
n=7 total: ~26 min for 4 points.
