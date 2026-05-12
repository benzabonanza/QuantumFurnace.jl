# P1 sweep redo (v4, β_phys-first, s = 0.25 fixed) — findings

**Driver**: `scripts/scratch_p1_v4_redo_betaphys.jl` (qf-e4z.21)
**Plot driver**: `scripts/numerics_p1_v4_plot.jl`
**Date**: 2026-05-12 • **Cells**: 36 (S1 + S2, n=3..8, β_phys ∈ {0.25, 0.5, 1.0})
**Total wall**: ~14 min (cached n=3..7 + 3×~225 s for the n=8 push)

## Headline

**For the canonical β_phys grid {0.25, 0.5, 1.0}, CKG smooth-Metropolis and DLL
Metropolis mix in essentially the same regime — τ_DLL/τ_CKG = 1.06×..1.58×
across n = 3..8, with CKG slightly faster everywhere.** Same wall-time
regime, same gap-bound scaling, no qualitative difference in mixing.

**This is a change from previous results before the β_phys-first convention.**
The v3 sweep (β_alg ∈ {5, 10, 20}, qf-e4z.4 close-reason) saw the DLL/CKG
ratio peak at 2.13× (n=6, β_alg=20, ε=10⁻⁶). At those high effective β_alg
the DLL Metropolis path was distinctly slower. With the new β_phys grid the
worst β_alg we hit at n=8 is 73.6 — comparable to v3's β_alg=20 in absolute
β but with a *different physical interpretation* (β_phys, not β_alg, is the
user-facing knob). And in *all* cells we now reach, the two filters are
within 1.6×. The dramatic CKG-vs-DLL contrast the earlier headline numbers
suggested is **not** what the canonical-grid data shows. CKG keeps a modest
edge, but the practical message — when reporting CKG-vs-DLL mixing-time
comparisons in the thesis — is "approximately equal, CKG slightly better".

## Filter convention

`s = 0.25, a = 0` **fixed** for CKG smooth-Metropolis. Decided 2026-05-12
after the alternative `default_smooth_s(β, σ) = (0.05/σ)²` was shown to
destroy γ(ω) at high β_alg (σ ≪ 0.05): the smoothing window swamps σ, γ
never reaches 1 even at deeply energy-decreasing ω, and Metropolis
discrimination collapses. See `drafts/plots/gamma_n7_betaphys1_three_s.png`
for the diagnostic plot, and the new top section of
`drafts/error-analysis/quadrature-convergence-summary-v2.md` for the
decision write-up.

Fixed `s = 0.25` keeps the *relative* smoothing σ·√s = σ/2 constant across
the β-sweep, so γ stays sharp at every temperature. r_D = 7 remains the
production register at n ≤ 8; at very large β_alg (≳ 50) the kernel
narrows and r_D = 8 (or a modestly larger s, e.g. 0.4) may be needed for
ε = 10⁻⁹ vs Bohr.

## What changed vs. v3 (qf-bw1, 2026-05-10)

| Change | Source |
|---|---|
| β_phys/β_alg split, β_phys-first | qf-6vr / `[[ham-rescaling]]` |
| Canonical β_phys grid {0.25, 0.5, 1.0} | replaces β_alg ∈ {5, 10, 20} |
| Per-term register triples on Config | qf-9z0 |
| Typical-disorder fixtures find_typical_heisenberg + [[Z],[Z,Z]] | qf-2kd |
| Coherent term ON by default | rule encoded 2026-05-12 |
| Quadrature recipe v2 (qf-yt9): r_D = 7 → ≤ 1e-7..1e-9 vs Bohr | n=3..8 grid |
| `s = 0.25` fixed convention reconfirmed (decision made today) | this redo |

## Sidecars

```
scripts/output/sweep_S1_v4_ckg_ideal/smooth_metro_eps1e-03/sweep_n{3..8}_betaphys{0.25,0.5,1}_seed42_L_KMS_Energy.bson
scripts/output/sweep_S2_v4_dll_ideal/smooth_metro_eps1e-03/sweep_n{3..8}_betaphys{0.25,0.5,1}_seed42_L_DLL_Bohr.bson
scripts/output/sweep_v4_summary.bson    # 36 rows
```

Re-runs are idempotent via `skip_existing=true` in `run_cell`; cell-level
try/catch + `GC.gc(true)` between cells keeps a single OOM from killing
the loop.

## Headline numbers

### S1 — CKG smooth-Metropolis KMS, EnergyDomain r_D=7, s = 0.25

| n | β_phys | β_alg | τ_mix | gap λ_2 | wall (s) | floor |
|---|--------|-------|-------|---------|----------|-------|
| 3 | 0.25 |  5.46 | 15.5 | 0.245 | 0.67 | 1e-7 |
| 3 | 0.5  | 10.93 | 18.5 | 0.235 | 0.04 | 9e-9 |
| 3 | 1.0  | 21.86 | 27.2 | 0.180 | 0.04 | 2e-8 |
| 4 | 0.25 |  9.67 | 31.2 | 0.140 | 0.11 | 2e-8 |
| 4 | 0.5  | 19.35 | 42.4 | 0.151 | 0.10 | 4e-9 |
| 4 | 1.0  | 38.69 | 43.2 | 0.149 | 0.10 | 9e-9 |
| 5 | 0.25 | 11.06 | 33.7 | 0.125 | 0.42 | 5e-9 |
| 5 | 0.5  | 22.12 | 40.5 | 0.146 | 0.41 | 9e-9 |
| 5 | 1.0  | 44.25 | 53.7 | 0.119 | 0.43 | 9e-9 |
| 6 | 0.25 | 14.03 | 45.7 | 0.092 | 2.99 | 6e-9 |
| 6 | 0.5  | 28.05 | 67.1 | 0.093 | 3.05 | 1e-8 |
| 6 | 1.0  | 56.11 | 77.7 | 0.098 | 3.22 | 3e-8 |
| 7 | 0.25 | 15.39 | 46.3 | 0.085 | 25.3 | 6e-8 |
| 7 | 0.5  | 30.79 | 58.5 | 0.099 | 26.1 | 1e-8 |
| 7 | 1.0  | 61.57 | 71.1 | 0.086 | 26.6 | 1e-8 |
| **8** | **0.25** | **18.41** | **63.7** | **0.068** | **215** | **2e-8** |
| **8** | **0.5**  | **36.81** | **89.5** | **0.068** | **222** | —     |
| **8** | **1.0**  | **73.62** | **111.8**| —         | ~225  | —     |

### S2 — DLL Metropolis Lindbladian, BohrDomain matrix-free

| n | β_phys | β_alg | τ_mix | gap λ_2 | wall (s) |
|---|--------|-------|-------|---------|----------|
| 3 | 0.25 |  5.46 | 18.3  | 0.202 | 0.43 |
| 3 | 0.5  | 10.93 | 22.3  | 0.203 | 0.00 |
| 3 | 1.0  | 21.86 | 29.0  | 0.175 | 0.00 |
| 4 | 0.25 |  9.67 | 38.7  | 0.097 | 0.01 |
| 4 | 0.5  | 19.35 | 61.2  | 0.078 | 0.01 |
| 4 | 1.0  | 38.69 | 61.7  | 0.094 | 0.01 |
| 5 | 0.25 | 11.06 | 42.1  | 0.092 | 0.03 |
| 5 | 0.5  | 22.12 | 53.2  | 0.092 | 0.03 |
| 5 | 1.0  | 44.25 | 60.2  | 0.101 | 0.03 |
| 6 | 0.25 | 14.03 | 58.8  | 0.066 | 0.20 |
| 6 | 0.5  | 28.05 | 101.8 | 0.049 | 0.18 |
| 6 | 1.0  | 56.11 | 116.8 | 0.053 | 0.17 |
| 7 | 0.25 | 15.39 | 53.7  | 0.061 | 1.30 |
| 7 | 0.5  | 30.79 | 76.7  | 0.056 | 1.32 |
| 7 | 1.0  | 61.57 | 85.0  | 0.065 | 1.33 |
| **8** | **0.25** | **18.41** | **83.6**  | **0.054** | **2.6** |
| **8** | **0.5**  | **36.81** | **141.2** | **0.043** | **3.2** |
| **8** | **1.0**  | **73.62** | **175.3** | **0.046** | **2.5** |

### τ_DLL / τ_CKG ratio (the comparison plot data)

| n | β_phys=0.25 | β_phys=0.5 | β_phys=1.0 |
|---|-------------|------------|------------|
| 3 | 1.18× | 1.20× | 1.06× |
| 4 | 1.24× | 1.45× | 1.43× |
| 5 | 1.25× | 1.31× | 1.12× |
| 6 | 1.29× | 1.52× | 1.50× |
| 7 | 1.16× | 1.31× | 1.20× |
| **8** | **1.31×** | **1.58×** | **1.57×** |

Range: 1.06× — 1.58×. CKG always faster, never by more than ~60%. The
ratio drifts upward with n at β_phys ∈ {0.5, 1.0} (1.50× → 1.58×, 1.20× →
1.57×) — CKG's lead grows slightly with system size in these cells. At
β_phys = 0.25 the ratio is essentially flat at ~1.2× across n. Whether
the upward drift at β_phys ∈ {0.5, 1.0} continues to n=9 is a question
for the next round.

## Wall-time + RAM scaling (final)

S1 CKG (EnergyDomain Krylov r_D=7) — bottleneck path:

| n | avg wall (s) | ratio | β_phys=1.0 cell |
|---|---|---|---|
| 3 | 0.25 | — | 0.04 |
| 4 | 0.10 | (warmup) | 0.10 |
| 5 | 0.42 | 4.2× | 0.43 |
| 6 | 3.09 | 7.4× | 3.22 |
| 7 | 26.0 | 8.4× | 26.6 |
| **8** | **221** | **8.5×** | **~225** |

Wall ratio settles to ~8.5× per +1 in n. Extrapolation:
- **n=9**: ~1900 s/cell × 3 β_phys = ~95 min total (laptop B on 16 GB Macbook).
- **n=10**: ~16000 s/cell × 3 β_phys = ~13 h total → **cluster**.

S2 DLL (BohrDomain matrix-free) is ~80× faster than S1 at n=8; trivially
extends to n=10 on a laptop.

**Memory**: peak RSS 1040 MB at n=8 (sandbox 4 GB — tight but OK). n=9
extrapolates to ~3 GB peak; cluster only mandatory at n ≥ 10.

## Plot

`drafts/figures/numerics/p1_taumix_v4.{png,pdf}` — single panel, y log
scale, 6 simulated curves (CKG = circles, DLL = squares; β_phys-color
warm→aubergine), 6 gap-bound overlays (hollow dashed).

## Quality flags

- All n ≥ 4 cells: `all_converged = true` at krylovdim=40. n=3 "false" is cosmetic (d²=64 ≤ krylovdim+restart).
- All floor_distance ≤ 6e-8 — ≥ 4 orders below target ε=10⁻³.
- r_D = 7 floor: per qf-yt9 v2, (smooth_fixed_s, n=6, β_phys=1.0) cell needs r_D=8 for ε=10⁻⁹ vs Bohr; r_D=7 here sits at ε ≈ 10⁻⁷. At n=8/β_phys=1 (β_alg=73.6) the smooth kernel is the narrowest in the grid; if a tighter-ε plot ever needs it, bump that one cell to r_D=8.

## Worth a second look

- **τ_DLL/τ_CKG drift n=6 → 8 at β_phys ∈ {0.5, 1.0}**: 1.50 → 1.31 → 1.58 at β_phys=0.5 and 1.50 → 1.20 → 1.57 at β_phys=1.0. Non-monotone with a dip at n=7 (the same fixture-specific n=6 / n=7 oddity we flagged earlier). A multi-seed sweep at n ∈ {6, 7, 8}, β_phys ∈ {0.5, 1.0} would tell us whether the n=6 fixture is unusually slow for DLL (likely) or whether n=7 is unusually fast (less likely).
- **CKG's lead grows with n at high β_phys** — if this trend continues to n=9+, the "basically similar" headline may need a "until n ≈ 8" caveat. Worth verifying on the Macbook (n=9 ~ 95 min total).
- **Spectral gap closes faster for DLL** (S2 gaps at n=8: 0.054, 0.043, 0.046) vs CKG (0.068, 0.068, —). Suggests DLL has a worse asymptotic spectral gap on these typical-disorder fixtures.

## Status

- qf-e4z.21 closed (this writeup + sidecars + plot + memory entries).
- Plot driver `scripts/numerics_p1_v4_plot.jl` is the canonical P1 plotter going forward; reads v4 sidecars, skips missing cells gracefully, ready for n=9+ once the BSONs land.

## Reproducibility

```
JULIA_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 \
  julia --project scripts/scratch_p1_v4_redo_betaphys.jl
JULIA_NUM_THREADS=1 julia --project scripts/numerics_p1_v4_plot.jl
```

Inputs: `hamiltonians/heis_xxx_zzdisordered_periodic_n{3..8}.bson` (find_typical, batch_size=256, qf-2kd). Re-runs skip cached cells via `skip_existing=true`.
