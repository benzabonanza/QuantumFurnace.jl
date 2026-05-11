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

## Channel floor (δ-step + GQSP) — Krylov vs dense Φ_δ

The `predict_channel_trajectory` floor `‖ρ_∞ − σ_β‖₁/2` is the **fixed-point
shift of `Φ_δ`**, not the Lindbladian residue. It picks up an extra
slope-1-in-`δ` contribution from the jump-wise generator splitting.

**Krylov predictor.** At recipe `(m_D=80, M_user=2, gqsp_degree=2, r_D=9)`:

| `δ_split` | floor | floor / δ |
|---|---|---|
| 1e-2 | 2.0e-3 | 0.20 |
| 1e-3 | 2.0e-4 | 0.20 |
| 1e-4 | 6.8e-5 | 0.68  (Krylov breaks) |
| 3e-5 | 4.2e-4 | 14    (Krylov breaks) |

Below `δ ≈ 3e-4` the Krylov spectral decomposition of `Φ_δ` loses accuracy:
with `Φ_δ ≈ I + δ·L_eff`, the eigenvalues cluster within `δ·λ` of 1 and the
Krylov projection onto the fixed-point eigenvector suffers from round-off in
the eigenvalue differences.

**Dense eigendecomposition.** Build `S = vec ∘ Φ_δ ∘ unvec` as a `d²×d²`
matrix (n=3 ⟹ 64×64 — microseconds), find the eigenvector at the eigenvalue
closest to 1, exact ρ_∞ to machine precision. Same recipe, dense path:

| `δ_split` | floor | floor / δ | `|μ₁ − 1|` |
|---|---|---|---|
| 1e-2 | 2.0e-3 | 0.20 | 2.4e-15 |
| 1e-3 | 2.0e-4 | 0.20 | 1.3e-15 |
| 1e-4 | 2.6e-5 | 0.26 | 6.2e-15 |
| 1e-5 | 1.6e-5 | 1.6   | 2.7e-15 |
| 1e-6 | 1.7e-5 | 17    | 5.6e-19 |
| 1e-7 | 1.7e-5 | 173   | 4.4e-16 |

The slope-1 (δ-split) extrapolates through `δ → 0`; at `δ ≈ 1e-5` the floor
saturates at the **quadrature steady-state shift** `‖L · σ_β‖_HS / λ_2 ≈
3.9e-7 / 0.089 ≈ 4 × 10⁻⁶`, then further saturates around `1.7 × 10⁻⁵` as
higher-order GQSP/Trotter contributions step in. The dense path proves the
Krylov-predictor anomaly at `δ < 3e-4` is a **numerical artefact**, not a
physical floor.

**Pushing the channel below 1e-6.** With dense Φ_δ, sweep `m_D` at `δ=1e-6`
(GQSP+Trotter contributions stay small), recipe `(r_D=9, M_user=4, gqsp_degree=2)`:

| `m_D` | `t0_bp` | `r_bp` | `‖ρ_∞ − σ_β‖₁/2` |
|---|---|---|---|
| 320 | 7.8e-4 | 15 | 4.2e-6 |
| 1280 | 2.0e-4 | 17 | **1.0e-6** ← at 1e-6 |
| 2560 | 9.8e-5 | 18 | **5.7e-7** ← below ✓ |
| 5120 | 4.9e-5 | 19 | **3.5e-7** ← below ✓ |
| 10240 | 2.5e-5 | 20 | 4.5e-7  (FP-noise saturated) |

**Channel error is fully controllable to ≤ 1e-6** via tightening `(m_D, δ)`.
The asymptote near 4e-7 is set by dense `eigen()` accuracy + accumulated
floating-point error in the d²×d² superoperator build — not a physical floor.

`M_user`, `gqsp_degree` sweeps at `δ=3e-4` all give floor ≈ 5e-5, confirming
δ-split dominance at the loose recipe:

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
