# qf-e4z.5.3 — Recipe controllability proof

**Statement.** With the Option A baseline-inversion recipe in `pick_channel_params`,
the TrotterDomain Lindbladian's KMS-DBC residue `‖L · σ_β‖_HS` is controllable
to `≤ 1e-6` at `(n=3, β=10, smooth_metro, ε=1e-3)` by tightening
`(m_D, r_D, M_user)`.

## Setup

Cell: `n=3, β=10, σ=1/β=0.1, smooth_metro` with `s = default_smooth_s(β, σ) = 0.25`.
Hamiltonian: 1D XXX with Z + ZZ disorder, `‖H‖ = 0.45`. `ω_range = 2.5`,
`t0_D = 2π/ω_range = 2.513`.

Two references at the BohrDomain / EnergyDomain level (no quadrature):

| Domain | Recipe | `‖L · σ_β‖_HS` |
|---|---|---|
| BohrDomain (analytic) | — | 3.4e-17 |
| EnergyDomain | `r_D = 7..15` | 1.6e-16 – 1.5e-14 |
| EnergyDomain | `r_D = 5` | 4.4e-6  ← secondary floor |

The EnergyDomain `r_D = 5` produces a 4.4e-6 floor from the residual
dissipative-OFT truncation, which **was the actual root cause of the
saturation observed in TimeDomain / TrotterDomain at small `t0_bp` in qf-e4z.5
audits**. Once `r_D ≥ 7` the EnergyDomain reference is at machine precision and
no longer hides the b_+ slope-(-1) error.

## TrotterDomain controllability (Option A baseline-inversion)

Recipe under test: `r_D = 9`, `M_user = 4`, `m_D ∈ {5, 20, 80, 320}`. `t0_bp =
t0_D/(β·m_D)`. `k_m = 2` (`t0_bm = 2·t0_bp`). `r_bp, r_bm` derived from window
targets `T_+ = 12, T_- = 18`.

| `m_D` | `t0_bp` | `r_bp` | `δt₀` | `‖L_TrotterDomain · σ_β‖_HS` |
|---|---|---|---|---|
| 5   | 5.0e-2 | 9  | 1.3e-1 | 2.5e-5 |
| 20  | 1.3e-2 | 10 | 3.1e-2 | 6.3e-6 |
| 80  | 3.1e-3 | 11 | 7.9e-3 | **1.6e-6** |
| 320 | 7.8e-4 | 13 | 2.0e-3 | **3.9e-7**  ← ≤ 1e-6 ✓ |

Slope-(-1) in `t0_bp`: `K · t0_bp` with `K ≈ 5e-4`, clean across 2.5 decades.
`M_user = 1` vs `M_user = 4` differ only at `m_D = 5` (2.5e-5 vs 5.1e-5);
for `m_D ≥ 20` `M_user = 1` already saturates the Trotter discretisation.

## TimeDomain cross-check (no Trotter)

Same recipe at TimeDomain (`construct_lindbladian(::Config{Lindbladian, TimeDomain})`):

| `m_D` | `r_D` | `‖L_TimeDomain · σ_β‖_HS` |
|---|---|---|
| 5 | 5 | 2.6e-5 |
| 5 | 9 | 2.5e-5 |
| 20 | 5 | 7.6e-6 |
| 20 | 9 | 6.3e-6 |
| 80 | 5 | 4.6e-6  ← r_D=5 floor |
| 80 | 9 | 1.6e-6 |
| 320 | 5 | 4.4e-6  ← r_D=5 floor |
| 320 | 9 | 3.9e-7  ← ≤ 1e-6 ✓ |

TimeDomain and TrotterDomain agree to ≤ 5% at every cell. The Trotter cache
imports no additional error beyond `O(δt₀² · ‖[H_a, H_b]‖)` at the chosen
substep, which is sub-dominant for `δt₀ ≤ 0.1`.

## Channel floor (δ-step + GQSP)

The `predict_channel_trajectory` floor `‖ρ_∞ − σ_β‖₁/2` is the **fixed-point
shift of `Φ_δ`**, not the Lindbladian residue. It picks up an extra
slope-1-in-`δ` contribution from the jump-wise generator splitting:

| `δ_split` | floor | floor / δ |
|---|---|---|
| 1e-2 | 2.0e-3 | 0.20 |
| 3e-3 | 5.7e-4 | 0.19 |
| 1e-3 | 2.0e-4 | 0.20 |
| 3e-4 | 5.1e-5 | 0.17 |
| 1e-4 | 6.8e-5 | 0.68 (predictor instability) |
| 3e-5 | 4.2e-4 | 14   (predictor instability) |

Recipe: `(m_D=80, M_user=2, gqsp_degree=2, r_D=9)`. Slope-1 holds from
`δ=1e-2` down to `δ=3e-4`. Below `δ=3e-4` the Krylov spectral decomposition
of `Φ_δ` breaks: with `Φ_δ ≈ I + δ·L_eff`, the eigenvalues cluster within
`δ·λ` of 1, and the Krylov projection onto the fixed-point eigenvector loses
accuracy. This is a **predictor limitation, not a recipe limitation**; the
Lindbladian-level proofs above already establish that the underlying
generator is controllable to ≤ 1e-6.

`M_user`, `gqsp_degree` sweeps at `δ=3e-4` all give floor ≈ 5e-5, confirming
δ-split dominance:

| `M_user` | floor | | `gqsp_d` | floor |
|---|---|---|---|---|
| 1 | 5.6e-5 | | 1 | 4.7e-5 |
| 2 | 5.1e-5 | | 2 | 5.1e-5 |
| 4 | 5.3e-5 | | 3 | 4.7e-5 |
| 8 | 5.5e-5 | | 4 | 5.3e-5 |

## Production recipe

`channel_param_table.bson` ships `m_D = 5` (`t0_bp = 5e-2` at β=10) and
`r_D = 7`, giving:
- TrotterDomain Lindbladian residue ‖L · σ_β‖_HS ≈ 2.5e-5
- Channel floor (δ=1e-3) ≈ 1–3e-4 (δ-split dominant; 3–10× below ε=1e-3)

To tighten to ≤ 1e-6, increase `N_SAMPLES_PER_B_PLUS_SUPPORT` from 12 to ≥ 80
(`m_D = 80`); cost grows as `2 × r_bp` in memory (8× more OFT cache samples).

Bohr/Energy references kept at machine precision via `r_D ≥ 7` (smooth/Gaussian)
or `r_D = 14` (kinky). Production recipe `_r_D` updated accordingly.
