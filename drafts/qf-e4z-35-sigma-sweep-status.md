# qf-e4z.35 σ-sweep — status (paused at n=8)

**Date:** 2026-05-18 (sweep paused after ~5 hours; user opted to finish n=8 later)

## What ran

- Driver: `scripts/scratch_qf_e4z_35_sigma_sweep_plot.jl`
- Analyzer: `scripts/analyze_qf_e4z_35_sigma_sweep_plot.jl`
- Floor-cell diagnostic: `scripts/scratch_qf_e4z_35_diagnostic.jl`
- Pipeline: CKG = KMS + EnergyDomain + smooth-Metro (s=0.25, a=0),
  `rho_0 = |+⟩⟨+|^⊗N`, single Arnoldi pass via `predict_lindbladian_trajectory`,
  τ_mix via `eigenmode_mixing_time` at ε=1e-3, t_max=500, n_grid=81.
- Grid: n ∈ {3..8}, β_phys ∈ {0.25, 0.5, 1.0}, σ_factor c ∈ {0.25, 0.5, 0.75, 1.0, 1.5, 2.0}, seed=46
- **r_D = 8** (raised from issue spec's r_D=7; see decision below)
- HS norm via matrix-free GKL (krylovdim=20) and d_{1→1} via Kossakowski M
  baked into every sidecar.

## Completion state

- **94 / 108 cells done** (saved as BSON sidecars under
  `scripts/output/sweep_qf_e4z_35_sigma_sweep_plot/ckg/`).
- **n = 3..7 fully complete** (90 cells).
- **n = 8 partial**: only β_phys=0.25 c ∈ {0.25, 0.5, 0.75, 1.0} done (4 cells).
- 14 n=8 cells remain. Resume = re-run the same driver; `skip_existing=true`
  skips the 94 cached cells automatically.
- **No `:floor` cells anywhere** in the 94 completed; all `extrapolated`,
  all Arnoldi converged. r_D=8 was the right call.

## r_D decision (7 → 8)

Initial smoke test at r_D=7 (the issue spec's baseline) hit a `:floor` cell
at (n=5, β_phys=1, c=0.25) — τ_mix=Inf because the Krylov-captured fixed
point sat 1.98e-3 away from the analytic Gibbs state, above ε=1e-3.

Focused diagnostic (`scratch_qf_e4z_35_diagnostic.jl`) at the same cell:

| variant | gap_alg | τ_mix_alg | floor_dist | source |
|---|---|---|---|---|
| r_D=7, kdim ∈ {60..200} | 0.21058 | Inf | 1.98e-3 | floor |
| r_D=8, kdim=100 | 0.21063 | 36.39 | 1.15e-5 | extrapolated |
| r_D=9, kdim=100 | 0.21063 | 36.39 | 2.56e-12 | extrapolated |
| r_D=10, kdim=100 | 0.21063 | 36.39 | 1.78e-15 | extrapolated |
| **BohrDomain dense** (exact ref) | **0.21063** | — | 0 | exact |

Reading: r_D=7 with shared-ω_range from σ_max=2/β_alg under-resolves the
σ=c·1/β_alg kink at the smallest c. Bumping kdim from 60 to 200 does
**nothing** — Krylov is not the bottleneck. r_D=8 doubles grid pts per kink
width (from 0.275 to 0.55) and recovers the Bohr-exact gap to <1e-5 and a
clean τ_mix extrapolation.

## n=8 wall-time issue

At n=8 (d=256), per-cell wall jumped from ~10 s (n=7) to **25–55 minutes**:

| n | wall typ. | dominant cost |
|---|---|---|
| 7 | 200–470 s | apply_lindbladian at d=128, r_D=8 (≈256 energy bits) |
| 8 | 1500–3300 s | apply_lindbladian at d=256, r_D=8 — **memory-bandwidth wall**: working set d²·2^r_D ≈ 16 MComplexF64 per matvec = 256 MB read/write, well past L2 |

Each n=8 cell also exhibits **GC drift**: within β_phys=0.25, wall grew
1557 → 2400 → 2810 → 3240 s across c ∈ {0.25, 0.5, 0.75, 1.0}, even though
the matvec count stayed flat (167–231). The trajectory + HS-norm Krylov
+ Kossakowski-M build inside one Julia session accumulate enough live
memory to slow the GC tail.

Projected total n=8 time at current pace: **~12 hours for the remaining
14 cells**. Paused.

## Patterns in the n=3..7 data (single-seed=46)

Reading the analyzer table (`drafts/qf-e4z-35-sigma-sweep-partial-n3-7.txt`):

- **τ_mix monotonically increases in c at every (n, β_phys)** — wider σ means
  slower mixing, even though it widens the Metropolis kink window.
  At β_phys=0.25 the trend is steep (τ grows 3–5× from c=0.25 to c=2);
  at β_phys=1.0 the small-c plateau is essentially flat (n=3 example: 1.04,
  1.14, 1.29, 1.50, 2.15, 3.27) — most of the σ-growth is concentrated
  at c ≥ 1.
- **gap_phys monotonically decreases in c** — small σ gives the largest
  spectral gap. At β=0.25 the gap drops ~5× from c=0.25 to c=2; at β=1.0
  the gap is largely flat for c ≤ 1 then drops at c ≥ 1.5.
- **‖L‖_HS_phys also decreases monotonically in c** — wider σ smears α(ν, ν′)
  out, lowering the operator norm.
- **τ_mix · ‖L‖_HS generally increases in c** — the σ-trend in τ_mix is NOT
  pure rate-scale change. There is a structural slow-down at wider σ
  on top of the ‖L‖ shrinkage.
- **gap / ‖L‖_HS decreases in c** — relative spectral gap shrinks as σ
  widens; the slowest mode sits closer to the bulk at wider σ.
- **τ_mix · gap (rate-resolved mixing)** is approximately constant in c
  across most cells, drifting downward at large c (esp. at colder β).
  This says the slowest-mode coefficient |c_2| in `eigenmode_mixing_time`
  is roughly σ-independent — τ_mix tracks 1/gap fairly well.

## Resume plan (for the deferred n=8 finish)

```bash
JULIA_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 \
  julia --project scripts/scratch_qf_e4z_35_sigma_sweep_plot.jl
```

`skip_existing=true` in the driver makes this idempotent — the 94 cached
sidecars are skipped, only the 14 missing n=8 cells run. Expected wall
time: **~10–12 hours** at the current per-cell rate. Better to run on
the cluster or split into chunks (one β_phys at a time).

If the n=8 GC-drift problem looks bad even on the cluster, two mitigations:
- Restart Julia between (n, β_phys) blocks (each fresh process pays AOT
  compile cost ~10 s but starts at GC zero). The cleanest version is a
  short bash launcher that loops over (n, β_phys) blocks; or
- Drop kdim_p2 use (we already only do a single Arnoldi pass per qf-e4z.30) —
  no obvious savings without changing the algorithm.

## Files

- `scripts/scratch_qf_e4z_35_sigma_sweep_plot.jl` — driver
- `scripts/analyze_qf_e4z_35_sigma_sweep_plot.jl` — analyzer (writes
  `scripts/output/sweep_qf_e4z_35_sigma_sweep_plot_summary_stats.bson`)
- `scripts/scratch_qf_e4z_35_diagnostic.jl` — r_D / kdim floor-cell diagnostic
- `scripts/output/sweep_qf_e4z_35_sigma_sweep_plot/ckg/` — 94 sidecars
- `scripts/output/sweep_qf_e4z_35_sigma_sweep_plot_summary_stats.bson` —
  analyzer output (single-row-per-cell, all metadata)
- `drafts/qf-e4z-35-sigma-sweep-partial-n3-7.txt` — printed analyzer table
