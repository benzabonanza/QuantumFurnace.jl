---
name: β_phys vs β_alg convention (qf-6vr)
description: Hamiltonian fixtures store a rescaled spectrum, so any inverse temperature `β` has two meanings — physical (against `H_phys`) and algorithm-side (against the rescaled `H` stored in the HamHam). qf-6vr / Phase qf-bphys splits them explicitly.
type: project
---

## Why this convention exists

`_rescaling_and_shift_factors` (`src/hamiltonian.jl`) divides every
Hamiltonian fixture by `rescaling_factor ≈ 2(λ_max − λ_min)/(1 − ε)` so the
stored spectrum sits in `[0, 0.45]` for *every* `n`. For extensive
Hamiltonians (Heisenberg, TFIM) this factor scales roughly linearly with the
qubit count `n`. Consequence: a sweep "at fixed β_alg across n" silently
varies the *physical* β by the same factor (and vice versa). The qf-now
scaling-fit `τ ≈ 3.48·n^0.82·β^0.42` was diagnosed as partly an artifact of
this β_phys / n coupling (`drafts/scaling-fit-physics-check.md`).

## Two β scales

- `β_phys` — against the un-rescaled `H_phys`. The thing a physicist types.
- `β_alg = β_phys · ham.rescaling_factor` — against the rescaled spectrum
  that the simulator actually sees. Equals the legacy `cfg.beta`.

## API surface (post-qf-6vr)

| Object | API | Convention |
|---|---|---|
| `HamHam` | `HamHam(raw, β)` (positional) | β = β_alg, legacy semantics |
| `HamHam` | `HamHam(raw; beta_phys=…)` (keyword) | derives β_alg internally |
| `HamHam` | `beta_alg(ham, β_phys)` / `beta_phys(ham, β_alg)` | scalar conversions |
| `Config` | `cfg.beta` (= `beta_alg(cfg)`) | β_alg, required field |
| `Config` | `cfg.beta_phys` (= `beta_phys(cfg)`) | optional, `Union{T, Nothing}` |
| `validate_config!(cfg, ham)` | 2-arg method | enforces `cfg.beta ≈ cfg.beta_phys · ham.rescaling_factor` when both set |
| `sweep_mixing_times` / `sweep_channel_mixing` | positional `beta_values` | β_alg list (legacy) |
| `sweep_mixing_times` / `sweep_channel_mixing` | kwarg `beta_phys_values` | β_phys list (qf-6vr); mutually exclusive with positional |
| Sidecar BSON | `:beta_phys`, `:beta_alg`, `:rescaling_factor` | always emitted (both modes) |
| Sidecar filename tag | `beta<β_alg>` vs `betaphys<β_phys>` | switched on input mode |
| `fit_scaling(::Vector{<:NamedTuple})` | `beta_kind = :auto` | prefers `:beta_phys`, falls back to `:beta_alg`/`:beta` |
| `ScalingFit.beta_kind` | `:phys` or `:alg` | sets the label in `formula_string` (`β_phys^y` vs `β_alg^y`) |
| Test constants | `BETA` (= `BETA_ALG`) | β_alg, legacy preserved |
| Test constants | `BETA_PHYS` (= `BETA / TEST_HAM.rescaling_factor`), `N3_BETA_PHYS` | β_phys for n=4 / n=3 fixtures |

## Default β grids

- **Drivers calibrated at β_alg** (legacy): `[5.0, 10.0, 20.0]`
- **Canonical qf-6vr β_phys grid** (decided 2026-05-11, after the
  Gibbs-state entropy analysis): `[0.25, 0.5, 1.0]`. Used by every
  migrated production driver:
  `numerics_ckg_vs_dll_taumix_comparison.jl`,
  `numerics_sweep_S3.jl`,
  `numerics_scaling_fit_ckg_smooth_metro.jl` (Task 8).
- **Auxiliary drivers** kept on the legacy β_alg path (calibrated at the
  β_alg ranges that fed the existing parameter table):
  `numerics_p1_*.jl`, `numerics_S3_audit_floors.jl`, `numerics_param_table.jl`.

### Why this grid

Choice constraints from the Gibbs-state entropy analysis (xxx_zzdis fixture,
seed 42, n ∈ {3..10}):

- **Lower 0.25**: smallest β_phys with non-trivial thermal contrast.
  S/log(d) ≈ 0.80 nearly *uniformly across n=3..10* — the same fractional
  entropy at every system size. Below 0.25 (e.g. β_phys=0.1) one gets
  S/log(d) > 0.96, essentially infinite-temperature; the Gibbs state is
  maximally mixed and τ_mix(β) carries no β-signal worth fitting.
- **Upper 1.0**: practical ceiling. β_alg at the corners:
  n=3 → 26.8, n=6 → 52.9, n=10 → 90.3, n=11 → ~99. The OFT filter width
  σ = 1/β_alg drops to ~0.01 at n=10. The legacy quadrature recipe
  (`quadrature_register_recipe_qf_7xt`) is calibrated up to β_alg ≈ 20;
  the register sizing audit for β_alg ≈ 90 is feasible but ugly. Above
  β_phys=1 the audit becomes pathological without further σ-tightening.
- **Factor-4 log-spacing**: same span as the legacy grid `{5, 10, 20}`,
  which the original physics-check (`drafts/scaling-fit-physics-check.md`)
  already noted as "barely enough" for M0 vs M1 discrimination (ΔAICc≈3.6
  in the "weak preference" zone). A wider span isn't available under the
  current constraints; if M0/M1 discrimination matters, add a 4th/5th
  point with denser log-spacing later.

### `s` blowup caveat at β_phys=1, large n

The qf-96o rule `default_smooth_s(β, σ) = (0.05/σ)²` is designed to keep
the *absolute* kink width `σ·√s = 0.05` constant across β-sweeps. With
σ = 1/β_alg, the rule reads `s = (0.05·β_alg)²`. At the high corner
(β_phys=1, n=11): β_alg ≈ 99, so `s ≈ 25` — far outside the historical
[0, 1] range of s values we have tested.

**Status (2026-05-11):** all migrated drivers hold `s = 0.25` fixed (the
legacy thesis convention), NOT the `default_smooth_s` rule. We do not yet
know:

1. Whether a smooth-Metro γ-rate with s = O(10) is physically meaningful
   (the smoothing kernel may suppress the transition rate at positive ν
   so heavily that the dissipator's spectral coverage collapses).
2. How much s = O(10) suppresses γ vs the optimal kinky-Metropolis
   value γ_kinky ≈ exp(-βν/2) at low-T.
3. Whether `τ_mix` at β_phys=1 cells is dominated by the kink-resolution
   issue (in which case `default_smooth_s` would help) or by the
   σ-narrowness (filter spectral coverage issue, which neither s rule
   addresses).

**Decision deferred** until the first full β_phys ∈ {0.25, 0.5, 1.0}
sweep lands. Mitigation choices then are:

- (a) Switch to `default_smooth_s` and accept large s, characterising the
  γ-rate suppression at one cell first.
- (b) Hold σ = 0.1 (or σ = c/β with c < 1, per `sigma_sweep_findings_qf_bw1`
  optimum at σ ≈ 0.25/β) fixed independent of β, decoupling the kink
  width from β.
- (c) Reinstate the fixed s = 0.25 value (current driver default) and
  document the kink-width drift as a known artifact of the β-sweep.

## Why the minimal-churn strategy

Strict rename `Config.beta → Config.beta_alg` was rejected in favour of
keeping `cfg.beta` = β_alg with a new `cfg.beta_phys` companion. The rename
would have forced ~80 src callsites and ~30 test files to change for zero
new physics. The minimal-churn version preserves every existing test
verbatim while making β_phys explicit at the user-script boundary, which is
where the operational confusion was.

## Risk register

- **β_alg at n=11, β_phys=3** — `ham.rescaling_factor` for the
  `xxx_zzdisordered` family is ≈ 30–35 at n ≥ 8, so β_phys=3 gives
  β_alg ≈ 90–105. Smoke test at (n=3, β_phys=1) already produced
  `mixing_time_source = :floor` with `target_ε = 1e-3` against the
  EnergyDomain CKG path — the asymptotic floor of the captured Krylov
  subspace exceeds `target_ε` at high β_alg. Mitigations the user can
  apply before launching the full Task 8 sweep:
  1. Bump `target_epsilon` to `1e-2`.
  2. Increase `spectral_krylovdim` from 60.
  3. Switch low-β_phys cells to BohrDomain (exact, no register sizing).
- **Param table cell coverage** — `numerics_param_table.jl` is calibrated
  on β_alg ∈ {5, 10, 20}. β_phys-first channel sweeps at β_phys=3 at n=8
  need β_alg ≈ 90 cells that don't exist in the table; rebuild the table
  before running `numerics_sweep_S3.jl` Phase B (n ≥ 7).

## Migration utility

`scripts/migrate_bson_beta_phys.jl` walks legacy sidecar trees and writes
`*.betaphys.bson` companions carrying `:beta_phys`, `:beta_alg`,
`:rescaling_factor` derived from the matching fixture. Idempotent.
Originals are preserved unchanged. The default scan covers
`scripts/output/` and `drafts/figures/numerics/sweep_cache/`; pass `--dir`
to target a specific tree.
