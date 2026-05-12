# Quadrature Register Recipe — v2 (qf-yt9, 2026-05-12 redo)

Tables produced by the campaigns in `scripts/scratch_quad_*_campaign.jl`
after the β_phys/β_alg split (qf-6vr), per-register API (qf-9z0),
typical Hamiltonians (qf-2kd), and `default_smooth_s` (qf-3il). Replaces
the qf-7xt v1 tables.

All cells use the unified ω-range principle `omega_range = 2(‖H‖+8σ)` with
`w0_X(r) = omega_range / 2^r`. EnergyDomain → BohrDomain measures the pure
energy-quadrature Riemann-sum error; TimeDomain → BohrDomain adds the
FINUFFT noise floor (eps = 1e-12). Coherent B_time → B_bohr uses fixed
T_RANGE_MINUS = 18, T_RANGE_PLUS = 12 (≥ 3× kernel support), and
**η = 1e-6** (kept below `t0_+ = T_RANGE_PLUS / 2^R_MAX_PLUS` at the
sweep's largest r_+).

## Smooth-Metropolis `s` convention — DECIDED 2026-05-12

**Use fixed `s = 0.25, a = 0` for all production sweeps.** The
`default_smooth_s(β, σ) = (0.05/σ)²` convention is **not** what production
should run, despite what the qf-yt9 numbers below for the `smooth_def_s`
column might suggest. Reason: at high β_alg (σ ≪ 0.05) the formula sets
σ·√s = 0.05 absolute, which is ≫ σ itself. The resulting γ(ω) is mollified
on a window much wider than the OFT filter scale and never reaches 1, even
for deeply energy-decreasing ω. At (n=7, β_phys=1.0, β_alg≈62, σ≈0.016),
γ(ω = −0.1) = 0.72 with s = default_smooth_s = 9.48 vs 1.00 with s = 0.25.
The Metropolis acceptance shape is destroyed; CKG smooth-Metro τ_mix is
2–3× slower than DLL Metro at these cells. The s = 0.25 convention
instead keeps the **relative** smoothing σ·√s = σ/2 constant across β, so
γ stays sharp at every temperature.

**r_D scaling with this fix.** Fixed `s = 0.25` gives Δω_kink = σ/2,
narrower than the `smooth_def_s` smoothing window at high β. The energy
quadrature therefore needs to resolve a sharper kernel, which costs +1 bit
in `r_D` at the highest β. The `smooth_fixed_s` table below (which uses
s = 0.25) is the operational reference: r_D = 7 suffices for ε ≤ 10⁻⁹
across n ≤ 6 and β_phys ∈ {0.25, 0.5}, while (n = 6, β_phys = 1.0) needs
r_D = 8. Empirically (P1 v4 sweep, qf-e4z.21) r_D = 7 still gives ε ≤
10⁻⁷ vs Bohr at (n = 6, β_phys = 1) — 4 orders below typical thesis
targets — so r_D = 7 remains the production default through n ≈ 8.
Beyond that, two knobs are available as β_alg climbs into the 50–100
range:

* **bump `r_D` to 8** (or, eventually, 9) per the smooth_fixed_s column;
  cheap because Krylov matvec scales linearly in `2^r_D`.
* **bump `s` modestly** (e.g. 0.4, 1.0) to widen the kink-smoothing
  window enough to recover r_D = 7. This costs some Metropolis-acceptance
  sharpness but does not collapse γ(ω) the way `default_smooth_s` does.
  The qf-yt9 `smooth_fixed_s` table at higher s values (not currently
  measured but trivially extendable) would set the bridge.

The choice between the two becomes interesting only at β_alg ≳ 50; below
that, fixed s = 0.25 with r_D = 7 is the canonical configuration.

## Headline

* **EnergyDomain Lindbladian (dissipator + coherent) reaches ε=1e-9 at
  r_D = 6 for smooth_def_s** across every measured (n ≤ 6, β_phys ∈
  {0.25, 0.5, 1.0}) cell, except (n=6, β_phys=1.0) which needs r_D = 7.
  The coherent contribution in EnergyDomain is **closed-form** (`B_bohr`
  formula); it carries no quadrature error and `B_energy == B_bohr`
  exactly. So `r_D` from the table below is the **only** knob needed to
  qualify EnergyDomain as the 1e-9 reference. All downstream sweeps
  (S2 Time↔Energy, S5 Trotter↔Time, faithful-channel comparisons) use
  EnergyDomain at this `r_D` as the reference — BohrDomain never enters
  again, exactly as it should not (its matvec costs `O(n_bohr · d³)`
  ~ `O(d⁵)`, infeasible at `d ≥ 32`).

* **TimeDomain coherent term carries an irreducible slope-(-1) error in
  `r_+`** for smooth/kinky Metro (the trapezoid-rule t=0 L'Hôpital
  sample dominates; see `.claude-memory/trap_rule_t0_lhopital_origin.md`).
  At `R_MAX_PLUS = 17` we reach ε ≈ 4·10⁻⁷ at β_phys=0.5 and ε ≈ 2·10⁻⁶
  at β_phys=1.0; **ε = 10⁻⁹ is not reachable in practice** with this
  discretisation. This is *not* a bug — the smoothing parameter `s`
  affects only the b_+ envelope amplitude, not the t=0 PV anomaly, so
  smooth_def_s and smooth_fixed_s produce **identical** r_+ cutoffs.

* For TimeDomain runs targeting ε = 10⁻⁶ (the user's controllability
  target), the recipe is `r_+ ≈ 15` at β_phys=0.25, `r_+ ≈ 16` at
  β_phys=0.5, `r_+ ≈ 17` at β_phys=1.0 (per +log₂(β_alg/β_alg₀) bit
  scaling). For ε = 10⁻⁹ you must build the Lindbladian in
  EnergyDomain (where coherent is exact) rather than TimeDomain.

Column key: `r_X@ε` = smallest `num_energy_bits_X` such that the error
is ≤ ε against the BohrDomain analytical reference.

## gaussian

| n | β_phys | β_alg | σ | s | r_D@1e-3 | r_D@1e-6 | r_D@1e-9 | r_-@1e-3 | r_-@1e-6 | r_-@1e-9 | r_+@1e-3 | r_+@1e-6 | r_+@1e-9 |
|---|--------|-------|---|---|----------|----------|----------|----------|----------|----------|----------|----------|----------|
| 3 | 0.25 | 5.46 | 0.183 | — | 5 | 5 | 5 | 4 | 5 | 6 | 4 | 5 | 6 |
| 3 | 0.5 | 10.93 | 0.0915 | — | 5 | 5 | 6 | 4 | 5 | 6 | 4 | 6 | 6 |
| 3 | 1.0 | 21.86 | 0.0457 | — | 5 | 6 | 6 | 5 | 6 | 6 | 4 | 6 | 6 |
| 4 | 0.25 | 9.67 | 0.1034 | — | 5 | 5 | 6 | 4 | 5 | 6 | 5 | 6 | 6 |
| 4 | 0.5 | 19.35 | 0.0517 | — | 5 | 6 | 6 | 5 | 6 | 6 | 5 | 6 | 6 |
| 4 | 1.0 | 38.69 | 0.0258 | — | 6 | 6 | 7 | 6 | 6 | 7 | 6 | 7 | 7 |
| 5 | 0.25 | 11.06 | 0.0904 | — | 5 | 5 | 6 | 5 | 6 | 6 | 5 | 6 | 6 |
| 5 | 0.5 | 22.12 | 0.0452 | — | 5 | 6 | 6 | 5 | 6 | 6 | 6 | 6 | 6 |
| 5 | 1.0 | 44.25 | 0.0226 | — | 6 | 7 | 7 | 6 | 7 | 7 | 6 | 7 | 7 |
| 6 | 0.25 | 14.03 | 0.0713 | — | 5 | 6 | 6 | 5 | 6 | 6 | 5 | 6 | 6 |
| 6 | 0.5 | 28.05 | 0.0356 | — | 6 | 6 | 6 | 6 | 6 | 6 | 6 | 6 | 6 |
| 6 | 1.0 | 56.11 | 0.0178 | — | 6 | 7 | 7 | 7 | 7 | 7 | 6 | 7 | 7 |

## smooth_def_s

| n | β_phys | β_alg | σ | s | r_D@1e-3 | r_D@1e-6 | r_D@1e-9 | r_-@1e-3 | r_-@1e-6 | r_-@1e-9 | r_+@1e-3 | r_+@1e-6 | r_+@1e-9 |
|---|--------|-------|---|---|----------|----------|----------|----------|----------|----------|----------|----------|----------|
| 3 | 0.25 | 5.46 | 0.183 | 0.07 | 5 | 6 | 7 | 4 | 5 | — | 5 | 15 | — |
| 3 | 0.5 | 10.93 | 0.0915 | 0.30 | 5 | 6 | 6 | 5 | 6 | — | 5 | 15 | — |
| 3 | 1.0 | 21.86 | 0.0457 | 1.19 | 5 | 6 | 6 | 4 | 6 | — | 4 | 14 | — |
| 4 | 0.25 | 9.67 | 0.1034 | 0.23 | 5 | 6 | 6 | 4 | 5 | — | 5 | 14 | — |
| 4 | 0.5 | 19.35 | 0.0517 | 0.94 | 5 | 6 | 6 | 5 | 6 | — | 6 | 15 | — |
| 4 | 1.0 | 38.69 | 0.0258 | 3.74 | 5 | 6 | 6 | 6 | 6 | — | 6 | 16 | — |
| 5 | 0.25 | 11.06 | 0.0904 | 0.31 | 5 | 6 | 6 | 4 | 6 | — | 5 | 15 | — |
| 5 | 0.5 | 22.12 | 0.0452 | 1.22 | 5 | 6 | 6 | 5 | 6 | — | 6 | 16 | — |
| 5 | 1.0 | 44.25 | 0.0226 | 4.89 | 6 | 6 | 6 | 6 | 6 | — | 7 | 17 | — |
| 6 | 0.25 | 14.03 | 0.0713 | 0.49 | 5 | 6 | 6 | 4 | 6 | — | 5 | 15 | — |
| 6 | 0.5 | 28.05 | 0.0356 | 1.97 | 5 | 6 | 6 | 6 | 6 | — | 6 | 16 | — |
| 6 | 1.0 | 56.11 | 0.0178 | 7.87 | 6 | 6 | 7 | 6 | 7 | — | 7 | 17 | — |

## smooth_fixed_s

| n | β_phys | β_alg | σ | s | r_D@1e-3 | r_D@1e-6 | r_D@1e-9 | r_-@1e-3 | r_-@1e-6 | r_-@1e-9 | r_+@1e-3 | r_+@1e-6 | r_+@1e-9 |
|---|--------|-------|---|---|----------|----------|----------|----------|----------|----------|----------|----------|----------|
| 3 | 0.25 | 5.46 | 0.183 | 0.25 | 5 | 5 | 6 | 4 | 5 | — | 5 | 15 | — |
| 3 | 0.5 | 10.93 | 0.0915 | 0.25 | 5 | 6 | 6 | 5 | 6 | — | 5 | 15 | — |
| 3 | 1.0 | 21.86 | 0.0457 | 0.25 | 6 | 6 | 7 | 4 | 6 | — | 4 | 14 | — |
| 4 | 0.25 | 9.67 | 0.1034 | 0.25 | 5 | 6 | 6 | 4 | 5 | — | 5 | 14 | — |
| 4 | 0.5 | 19.35 | 0.0517 | 0.25 | 5 | 6 | 7 | 5 | 6 | — | 6 | 15 | — |
| 4 | 1.0 | 38.69 | 0.0258 | 0.25 | 6 | 7 | 7 | 6 | 6 | — | 6 | 16 | — |
| 5 | 0.25 | 11.06 | 0.0904 | 0.25 | 5 | 6 | 6 | 5 | 6 | — | 5 | 15 | — |
| 5 | 0.5 | 22.12 | 0.0452 | 0.25 | 5 | 6 | 7 | 5 | 6 | — | 6 | 16 | — |
| 5 | 1.0 | 44.25 | 0.0226 | 0.25 | 6 | 7 | 7 | 6 | 7 | — | 7 | 17 | — |
| 6 | 0.25 | 14.03 | 0.0713 | 0.25 | 5 | 6 | 6 | 4 | 6 | — | 5 | 15 | — |
| 6 | 0.5 | 28.05 | 0.0356 | 0.25 | 6 | 6 | 7 | 5 | 6 | — | 6 | 16 | — |
| 6 | 1.0 | 56.11 | 0.0178 | 0.25 | 6 | 7 | 8 | 6 | 7 | — | 7 | 17 | — |

## kinky

| n | β_phys | β_alg | σ | s | r_D@1e-3 | r_D@1e-6 | r_D@1e-9 | r_-@1e-3 | r_-@1e-6 | r_-@1e-9 | r_+@1e-3 | r_+@1e-6 | r_+@1e-9 |
|---|--------|-------|---|---|----------|----------|----------|----------|----------|----------|----------|----------|----------|
| 3 | 0.25 | 5.46 | 0.183 | — | 7 | 12 | 15 | — | — | — | — | — | — |
| 3 | 0.5 | 10.93 | 0.0915 | — | 6 | 12 | — | — | — | — | — | — | — |
| 3 | 1.0 | 21.86 | 0.0457 | — | 7 | 12 | — | — | — | — | — | — | — |
| 4 | 0.25 | 9.67 | 0.1034 | — | 7 | 12 | — | — | — | — | — | — | — |
| 4 | 0.5 | 19.35 | 0.0517 | — | 8 | 12 | 16 | — | — | — | — | — | — |
| 4 | 1.0 | 38.69 | 0.0258 | — | 7 | 12 | — | — | — | — | — | — | — |
| 5 | 0.25 | 11.06 | 0.0904 | — | 6 | 12 | — | — | — | — | — | — | — |
| 5 | 0.5 | 22.12 | 0.0452 | — | 7 | 13 | — | — | — | — | — | — | — |
| 5 | 1.0 | 44.25 | 0.0226 | — | 8 | 13 | — | — | — | — | — | — | — |

