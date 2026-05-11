# Sweep BSON layout — channel & Lindbladian (qf-e4z)

Per-cell sidecar conventions for the numerics-chapter plot suite (`qf-e4z`).
One BSON per `(n, β, seed, [ε, filter,] construction, domain)` cell, written by
`sweep_mixing_times` (Lindbladian / S1 + S2) and `sweep_channel_mixing`
(channel / S3 + S5 + S6) under `output_dir`.  Plot scripts compose snapshots
without re-running anything.

## β convention (qf-6vr)

Each sidecar carries three β fields:

| field | meaning |
|---|---|
| `beta_phys` | physical inverse temperature against the un-rescaled Hamiltonian (`H_phys`) |
| `beta_alg`  | algorithm-side inverse temperature against the rescaled spectrum (`H_rescaled` stored in `ham.eigvals`) |
| `beta`      | back-compat alias = `beta_alg` (legacy readers continue to work) |
| `rescaling_factor` | `ham.rescaling_factor` for the cell's HamHam (≈ `2(λ_max − λ_min)/(1 − ε)`) — relates the two via `beta_alg = beta_phys · rescaling_factor` |

The driver writes a sidecar through either of two entry points:

- **β_phys-first (qf-6vr):** call `sweep_*(n_values; beta_phys_values = [β_phys, …])`. The sweep harness loads `ham` per cell, derives `beta_alg = beta_phys · ham.rescaling_factor`, and tags the sidecar filename with `betaphys<β_phys>`.
- **β_alg-first (legacy):** call `sweep_*(n_values, [β_alg, …])`. The sidecar filename keeps the historic `beta<β_alg>` tag; `beta_phys` is derived from `ham.rescaling_factor` and recorded for forward compatibility.

`fit_scaling(::Vector{<:NamedTuple}; beta_kind = :auto)` reads `:beta_phys` preferentially, falling back to `:beta_alg`/`:beta`. The returned `ScalingFit` carries `beta_kind ∈ {:phys, :alg}` so `formula_string` prints the correct β label (`β_phys^y` vs `β_alg^y`).

## Channel sweeps (`sweep_channel_mixing` — qf-e4z.2)

### Sidecar filename

```
channel_n<n>_beta<β_alg>_seed<seed>_eps<ε>_<filter>_<construction>_<domain>.bson      (legacy β_alg-first)
channel_n<n>_betaphys<β_phys>_seed<seed>_eps<ε>_<filter>_<construction>_<domain>.bson  (qf-6vr β_phys-first)
```

- `<β_alg>` / `<β_phys>` formatted as up to 6 decimals, trailing zeros stripped (`5`, `10`, `20.5`).
- The `betaphys` prefix prevents collisions between legacy β_alg-keyed caches and β_phys-first re-runs in the same `output_dir`.
- `<ε>` in scientific notation `%.0e` (e.g. `1e-03`, `1e-06`).
- `<filter>` ∈ {`gaussian`, `smooth_metro`, `kinky_metro`}.
- `<construction>` ∈ {`KMS`, `GNS`}.
- `<domain>` ∈ {`Trotter`, `Time`}.

### Cell schema

Each sidecar contains `Dict(:result => Dict(...))` whose entries reconstruct
the per-cell `NamedTuple` returned by `sweep_channel_mixing`. Keys:

| field | description |
|---|---|
| `n, beta, seed, eps, filter` | identifying tuple (`beta` = β_alg; back-compat) |
| `beta_phys, beta_alg, rescaling_factor` | qf-6vr β-pair: `beta_alg = beta_phys · rescaling_factor` |
| `family, construction, domain` | string tags (e.g. `xxx_zzdisordered`, `KMS`, `Trotter`) |
| `r_D, w0_D, t0_D` | dissipative-register triple |
| `r_bm, w0_bm, t0_bm` | outer-coherent register (b_-) |
| `r_bp, w0_bp, t0_bp` | inner-coherent register (b_+) |
| `M_D, M_bm, M_bp` | Trotter step counts (per-term, qf-d0w shared-δt₀ scheme) |
| `delta, eta` | splitting δ + Metropolis regularisation |
| `with_gqsp, gqsp_degree` | GQSP cost-model flags (qf-e4z.18) |
| `tau_mix, tau_mix_source` | mixing time + source ∈ {`:extrapolated`, `:floor`, `:nan`} (eigenmode schema, qf-e4y) |
| `lambda_gap_channel` | Lindbladian-equivalent gap from the channel's leading non-trivial eigenvalue (`(1−|μ_2|)/δ`) |
| `floor_distance` | asymptotic ``\| \rho_\infty - \sigma_\beta \|_1 / 2`` from the captured Krylov subspace — exposes the O(δ) channel shift |
| `n_steps_to_target, k_max, t_max` | trajectory bookkeeping |
| `achieved_dist_at_kmax` | trace distance at end of `k_grid` (= channel-shift floor as observed) |
| `total_matvecs, all_converged_predict` | predictor convergence flags |
| `oft_time_per_step, b_per_be_per_step, b_time_per_step, per_step_time` | per-step Hamiltonian-simulation time decomposition |
| `n_steps_total, total_ham_sim_time` | aggregated Ham-sim time over `T = τ_mix` |
| `wall_time_seconds, init_state, family_tag` | bookkeeping |

`tau_mix_source` resolution order (post-qf-e4y eigenmode schema):
1. **`:extrapolated`** — `eigenmode_mixing_time` bisected the closed-form trace distance `d(t) = ‖(ρ_∞ − σ_β) + Σ c_i e^{λ_eff_i t} R_i‖_1 / 2 = ε` (continuous-time, with `λ_eff_i = log(μ_i) / δ`). `tau_mix` is the bisection result.
2. **`:floor`** — `ε` is below the asymptotic floor `floor_distance = ‖ρ_∞ − σ_β‖_1 / 2`, so no `t` solves `d(t) = ε`. `tau_mix` is populated with the conservative `log(d / ε) / λ_gap_channel` bound (matches the prior `:gap` source's value for plot continuity).
3. **`:nan`** — predictor produced fewer than 2 eigenvalues, or the bisection bracket failed even after one 3× expansion.

The previous `:observed` (trajectory crossing) and `:gap` (gap-based bound) sources were collapsed into the `:floor` branch — the eigenmode formula is exact on the captured Krylov subspace, so a "trajectory observation" carries no extra information beyond the eigendecomposition.

### Parameter source

Every per-cell parameter (`r_D, w0_D, t0_D, r_b±, …, M_D, M_b±, δ, η, with_gqsp, gqsp_degree`) is sourced from `scripts/output/channel_param_table.bson` (qf-e4z.1 / P0a). Editing the recipe means re-running `scripts/numerics_param_table.jl`; downstream sweeps then re-load the new table on next invocation.

### Hamiltonian-simulation time cost model

`oft_time_per_step` and `b_per_be_per_step` come straight out of `compute_simulation_time`. `b_time_per_step` applies the GQSP multiplier (qf-e4z.18 cost model = MW2024 Eq. 46 / Form B):

```
b_time_per_step = (with_gqsp && coherent) ? 2 · gqsp_degree · b_per_be_per_step
                                          : b_per_be_per_step
```

The `2d` coefficient matches the live GQSP circuit (Form B, `d` controlled-`W` + `d` closed-controlled-`W†` slots, no uncontrolled tail) — see `drafts/sim-time-gqsp-audit.md` for the form derivation and the qf-e4z.19 history.

`per_step_time = 2 · oft_time_per_step + b_time_per_step` (forward + backward OFT plus one CoherentStep). `total_ham_sim_time = n_steps_total · per_step_time` over `T = τ_mix`.

## Lindbladian sweeps (`sweep_mixing_times`)

### Sidecar filename

```
sweep_n<n>_beta<β_alg>_seed<seed>_<mode>_<construction>_<domain>.bson         (legacy β_alg-first)
sweep_n<n>_betaphys<β_phys>_seed<seed>_<mode>_<construction>_<domain>.bson    (qf-6vr β_phys-first)
```

`mode` ∈ {`:L`, `:K`}; `domain` ∈ {`Bohr`, `Energy`}; `construction` ∈ {`KMS`, `DLL`}. The β tag is `beta` when called with the positional `beta_values` (β_alg list), `betaphys` when called with `beta_phys_values` (β_phys list, qf-6vr).

### Schema — split by `method` (qf-e4y)

The runner emits two distinct NamedTuple shapes depending on the `method` keyword. Plot scripts must branch on the `method` field (which is always emitted) before reading method-specific columns.

#### `method = :krylov` (production / thesis numerics)

Single-pass spectral expansion via `predict_lindbladian_trajectory`. τ_mix(ε) is bisected analytically on the captured eigendecomposition.

| field | description |
|---|---|
| `n, beta, seed, init_state, mode, method, construction, domain` | identifying tuple (`method = :krylov` here; `beta` = β_alg, back-compat) |
| `beta_phys, beta_alg, rescaling_factor` | qf-6vr β-pair (always emitted; β_phys derived from `ham.rescaling_factor` in legacy mode) |
| `filter_name, filter_kind, target_epsilon` | filter identifying tuple |
| `r_D, w0_D, t0_D` | dissipative-register triple |
| `gap_est` | smallest `|Re(λ_i)|` over non-steady eigenvalues from the predictor's Arnoldi |
| `t_max, t_max_factor, tau_mix_bound` | seed-bracket bookkeeping |
| `n_grid, total_matvecs, all_converged` | predictor diagnostics |
| `mixing_time, mixing_time_source` | τ_mix + source ∈ {`:extrapolated`, `:floor`, `:nan`} (same enum + resolution as the channel sweep above) |
| `floor_distance` | asymptotic `‖ρ_∞ − σ_β‖_1 / 2` |
| `wall_time` | per-cell wall time (s) |

NO `fitted_gap`, `r_squared`, or `converged_fit` on this path.

#### `method = :ode` (legacy / debug)

Matrix-free ODE integrator over `t_grid` followed by a bi-exponential curve fit on the trace-distance trajectory.

| field | description |
|---|---|
| `n, beta, seed, init_state, mode, method, construction, domain` | identifying tuple (`method = :ode` here; `beta` = β_alg, back-compat) |
| `beta_phys, beta_alg, rescaling_factor` | qf-6vr β-pair (always emitted) |
| `filter_name, filter_kind, target_epsilon` | filter identifying tuple |
| `r_D, w0_D, t0_D` | dissipative-register triple |
| `gap_est` | spectral gap from a separate `krylov_spectral_gap` Arnoldi pre-pass |
| `t_max, t_max_factor, tau_mix_bound, n_grid, total_matvecs, all_converged` | trajectory bookkeeping |
| `fitted_gap, r_squared, converged_fit` | bi-exponential curve-fit diagnostics |
| `mixing_time, mixing_time_source` | τ_mix + source ∈ {`:extrapolated`, `:observed`, `:nan`} |
| `wall_time` | per-cell wall time (s) |

The two paths are intentionally distinct: `:ode` is preserved for benchmark scripts that compare the eigenmode and biexp extrapolators (`scripts/benchmark_lmode_vs_kmode.jl`); `:krylov` is the default for thesis numerics.

## Loading sidecars

```julia
using BSON, QuantumFurnace
d = BSON.load("output/channel_n3_beta10_seed42_eps1e-03_smooth_metro_KMS_Time.bson",
              QuantumFurnace)
result = NamedTuple(d[:result])     # NamedTuple matching the schema above
```

## Composing plots

The plot scripts (`drafts/figures/numerics/*.jl`, to be added) glob the relevant sidecars, build a `DataFrame` keyed by the identifying tuple, and join with whatever else they need (e.g. ideal-Lindbladian `τ_mix` from the Lindbladian sweep, ham-sim time scaling). No simulation work in the plot scripts — they read sidecars only.
