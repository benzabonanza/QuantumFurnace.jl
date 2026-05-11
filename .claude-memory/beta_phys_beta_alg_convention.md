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
- **Drivers calibrated at β_phys** (qf-6vr): `[1.0, 2.0, 3.0]` — the
  operationally visible payoff of the refactor. Migrated drivers:
  `numerics_ckg_vs_dll_taumix_comparison.jl`,
  `numerics_sweep_S3.jl`,
  `numerics_scaling_fit_ckg_smooth_metro.jl` (Task 8).
- **Auxiliary drivers** kept on the legacy β_alg path (calibrated at the
  β_alg ranges that fed the existing parameter table):
  `numerics_p1_*.jl`, `numerics_S3_audit_floors.jl`, `numerics_param_table.jl`.

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
