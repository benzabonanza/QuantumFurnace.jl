# `predict_lindbladian_trajectory` TimeDomain bench — ε = 10⁻³ vs ε = 10⁻⁶

**Task**: qf-1au follow-up. Same grid (n = 3..6, β_phys ∈ {1, 2},
smooth_def_s, coherent ON, krylovdim = 30, multithreaded), now run at
two register triples corresponding to two target precision levels.

| Target ε | r_D | r_b_minus | r_b_plus | v2 recipe source |
|----------|-----|-----------|----------|------------------|
| 10⁻⁶ | 7 | 7 | 17 | smooth_def_s table @ 10⁻⁶ |
| 10⁻³ | 6 | 6 | 7  | smooth_def_s table @ 10⁻³ |

The killer knob is `r_b_plus`: 17 → 7 cuts the inner FINUFFT of
`B_time` by `2¹⁰ = 1024×`. `r_D` 7 → 6 halves the per-matvec ω-loop.

## Wall-time comparison

| n | β_phys | d | ws_1e-6 | ws_1e-3 | ws speed-up | krylov_1e-6 | krylov_1e-3 | krylov speed-up | total_1e-6 | total_1e-3 |
|---|--------|---|---------|---------|-------------|-------------|-------------|------------------|------------|------------|
| 3 | 1.0 | 8  | 3.691 *(JIT)* | 3.457 *(JIT)* | — | 1.500 | 1.505 | — | 5.191 | 4.963 |
| 3 | 2.0 | 8  | 0.200 | 0.025 | **8×** | 0.035 | 0.018 | 1.9× | 0.234 | 0.043 |
| 4 | 1.0 | 16 | 0.374 | 0.004 | **104×** | 0.087 | 0.042 | 2.1× | 0.461 | 0.046 |
| 4 | 2.0 | 16 | 0.229 | 0.009 | **25×** | 0.074 | 0.040 | 1.9× | 0.303 | 0.049 |
| 5 | 1.0 | 32 | 1.211 | 0.007 | **176×** | 0.252 | 0.149 | 1.7× | 1.463 | 0.156 |
| 5 | 2.0 | 32 | 0.654 | 0.007 | **95×** | 0.292 | 0.158 | 1.9× | 0.946 | 0.165 |
| 6 | 1.0 | 64 | 6.879 | 0.033 | **207×** | 2.060 | 1.047 | 2.0× | **8.939** | **1.080** |
| 6 | 2.0 | 64 | 3.748 | 0.034 | **111×** | 2.092 | 1.068 | 2.0× | **5.840** | **1.102** |

All cells timed at 30 matvecs (krylovdim). Wall numbers in seconds.

**Headline**: at n = 6 the cold-start
`predict_lindbladian_trajectory` cost drops from ~9 s to ~1 s — an
8× speedup at β_phys = 1 (5× at β_phys = 2 where the ε = 10⁻⁶ workspace
was already cheaper due to faster-decaying kernels). Workspace build,
which dominated ε = 10⁻⁶ at small n, collapses to milliseconds at
ε = 10⁻³.

## Does the loosened precision matter physically?

Spectral gap (the quantity that sets τ_mix) at the two precision levels:

| n | β_phys | gap @ ε=10⁻⁶ | gap @ ε=10⁻³ | Δ |
|---|--------|---------------|---------------|---|
| 3 | 1.0 | 1.478e−1 | 1.478e−1 | 0.00% |
| 3 | 2.0 | 5.830e−2 | 5.836e−2 | 0.10% |
| 4 | 1.0 | 1.135e−1 | 1.135e−1 | 0.03% |
| 4 | 2.0 | 3.816e−2 | 3.833e−2 | 0.45% |
| 5 | 1.0 | 5.205e−2 | 5.207e−2 | 0.03% |
| 5 | 2.0 | 5.626e−3 | 5.921e−3 | **5.23%** |
| 6 | 1.0 | 3.902e−2 | 3.901e−2 | 0.02% |
| 6 | 2.0 | 4.060e−3 | 3.398e−3 | **16.3%** |

At **β_phys = 1, all n**: gap shifts by < 0.05%. The loosened registers
are effectively invisible to τ_mix.

At **β_phys = 2, n ≥ 5**: gap shifts by 5–16%. This is the actual
ε = 10⁻³ quadrature noise hitting the Lindbladian spectrum — and
because the gap itself is small at low T (4×10⁻³ at n=6, β=2), a fixed
10⁻³ perturbation becomes a *relative* O(10%) on the gap. So
**τ_mix = 1/gap is reliable to ~10% at ε = 10⁻³, exact at ε = 10⁻⁶**.

So: if you care about precise τ_mix at low T (β_phys ≥ 2), use the
ε = 10⁻⁶ triple. If you just want order-of-magnitude or eigenmode
*shapes* (R_modes), ε = 10⁻³ is fine and ~8× cheaper. Most of our
plot work targets ε = 10⁻³ to 10⁻⁵ ([thesis target precision
split](../../.claude-memory/thesis_target_precision_1e6.md)), so the
ε = 10⁻³ recipe is the right default unless explicitly proving
algorithm-level controllability at 10⁻⁶ to 10⁻⁹.

## What's left of the wall at ε = 10⁻³

For n = 6 (1.08 s total at β = 1):

- **97% Krylov factorisation** (1.05 s). 30 `apply_lindbladian!` matvecs,
  each 35 ms. d² scaling dominated by the per-jump sandwich GEMM.
- **3% workspace build** (0.03 s). `B_time` 2D FINUFFT now has 2⁶ × 2⁷ =
  8 192 grid points instead of 2⁷ × 2¹⁷ = 16.8 M — trivially fast.

At ε = 10⁻³ the `workspace = ws` reuse-knob path is essentially
irrelevant: building a fresh workspace per call costs 30 ms, while the
Krylov is 1 s. Just build it inline.

## Files

- Script (unchanged): `scripts/scratch_predict_lindbladian_time_bench.jl`
- BSON sidecars (ε=10⁻³): `scripts/output/predict_lindbladian_time_bench_eps1e3/n{3..6}_bphys{1.0,2.0}_smooth_def_s_time.bson`
- BSON sidecars (ε=10⁻⁶): `scripts/output/predict_lindbladian_time_bench/n{3..6}_bphys{1.0,2.0}_smooth_def_s_time.bson`
- Beads: qf-1au (closed; this is a follow-up comparison)
